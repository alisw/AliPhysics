/* Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSEHFResonanceBuilder
// \brief Analysis task to produce trees of D-meson candidates for ML analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include "yaml-cpp/yaml.h"

#include <TRandom3.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliNeutralTrackParam.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliHFMLResponseD0toKpi.h"
#include "AliHFMLResponseDplustoKpipi.h"
#include "AliHFMLResponseDstartoD0pi.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisUtils.h"
#include "AliAODHandler.h"
#include "AliAODExtension.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSECharmHadronMLSelector.h"

#include "AliAnalysisTaskSEHFResonanceBuilder.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEHFResonanceBuilder);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEHFResonanceBuilder::AliAnalysisTaskSEHFResonanceBuilder() : AliAnalysisTaskSE()
{
    /// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEHFResonanceBuilder::AliAnalysisTaskSEHFResonanceBuilder(const char *name, int decayChannel, AliRDHFCuts *analysisCuts) :
    AliAnalysisTaskSE(name),
    fDecChannel(decayChannel)
{
    /// Standard constructor
    SetAnalysisCuts(analysisCuts);
    SetDecayChannel(decayChannel);

    DefineOutput(1, TList::Class());
    switch(fDecChannel){
        case kD0toKpi:
            DefineOutput(2,AliRDHFCutsD0toKpi::Class());       //Cut object for D0
        break;
        case kDplustoKpipi:
            DefineOutput(2,AliRDHFCutsDplustoKpipi::Class());  //Cut object for Dplus
        break;
        case kDstartoD0pi:
            DefineOutput(2,AliRDHFCutsDStartoKpipi::Class());  //Cut object for D*
        break;
    }
    DefineOutput(3, AliNormalizationCounter::Class());
    DefineOutput(4, TNtuple::Class());
}

//________________________________________________________________________
AliAnalysisTaskSEHFResonanceBuilder::~AliAnalysisTaskSEHFResonanceBuilder()
{
    // Destructor
    delete fOutput;
    delete fCounter;
    delete fListCuts;
    delete fRDCuts;
    delete fNtupleCharmReso;

    if (fApplyML && fMLResponse)
        delete fMLResponse;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFResonanceBuilder::LocalInit()
{
    // Initialization

    switch (fDecChannel)
    {
        case kD0toKpi:
        {
            AliRDHFCutsD0toKpi *copycut = new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi *>(fRDCuts)));
            PostData(2, copycut);
        }
        break;
        case kDplustoKpipi:
        {
            AliRDHFCutsDplustoKpipi *copycut = new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi *>(fRDCuts)));
            PostData(2, copycut);
        }
        break;
        case kDstartoD0pi:
        {
            AliRDHFCutsDStartoKpipi *copycut = new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi *>(fRDCuts)));
            PostData(2, copycut);
        }
        break;
    }

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFResonanceBuilder::UserCreateOutputObjects()
{

    if (fInvMassResoPiMin.size() != fInvMassResoPiMax.size())
        AliFatal("Different size of fInvMassResoPiMin and fInvMassResoPiMax");
    if (fInvMassResoKaMin.size() != fInvMassResoKaMax.size())
        AliFatal("Different size of fInvMassResoKaMin and fInvMassResoKaMax");
    if (fInvMassResoPrMin.size() != fInvMassResoPrMax.size())
        AliFatal("Different size of fInvMassResoPrMin and fInvMassResoPrMax");

    if (std::accumulate(fEnableBachelor.begin(), fEnableBachelor.end(), 0) && std::accumulate(fEnableV0.begin(), fEnableV0.end(), 0))
        AliFatal("Combination with bachelors and V0s simultenously is not supported");

    /// Create the output container
    //

    // Several histograms are more conveniently managed in a TList
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos");

    fHistNEvents = new TH1F("fHistNEvents", "number of events ", 15, -0.5, 14.5);
    fHistNEvents->GetXaxis()->SetBinLabel(1, "no. read events");
    fHistNEvents->GetXaxis()->SetBinLabel(2, "no. matched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(3, "no. mismatched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(4, "no. analysed events");
    fHistNEvents->GetXaxis()->SetBinLabel(5, "no. passing IsEvSelected");
    fHistNEvents->GetXaxis()->SetBinLabel(6, "no. rejected due to trigger");
    fHistNEvents->GetXaxis()->SetBinLabel(7, "no. rejected due to not reco vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(8, "no. rejected for contr vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(9, "no. rejected for vertex out of accept");
    fHistNEvents->GetXaxis()->SetBinLabel(10, "no. rejected for pileup events");
    fHistNEvents->GetXaxis()->SetBinLabel(11, "no. of out centrality events");
    fHistNEvents->GetXaxis()->SetBinLabel(12, "no. rejected for phys sel");
    fHistNEvents->GetXaxis()->SetBinLabel(13, "no. of D candidates");
    fHistNEvents->GetXaxis()->SetBinLabel(14, "no. of D after filtering cuts");
    fHistNEvents->GetXaxis()->SetBinLabel(15, "no. of D after selection cuts");
    fHistNEvents->GetXaxis()->SetNdivisions(1, false);
    fHistNEvents->SetMinimum(0);
    fOutput->Add(fHistNEvents);

    // QA histograms
    std::array<std::string, kNumBachIDs> partLabel = {"Pi", "Ka", "Pr", "De"};
    for (int iHypo{0}; iHypo<kNumBachIDs; ++iHypo) {
        fHistNsigmaTPCSelBach[iHypo] = new TH2F(Form("fHistNsigmaTPC%sSelBach", partLabel[iHypo].data()), Form(";#it{p} (GeV/#it{c});#it{N}_{#sigma}^{ TPC} (%s)", partLabel[iHypo].data()), 100, 0., 10., 200, -10., 10.);
        fHistNsigmaTOFSelBach[iHypo] = new TH2F(Form("fHistNsigmaTOF%sSelBach", partLabel[iHypo].data()), Form(";#it{p} (GeV/#it{c});#it{N}_{#sigma}^{ TOF} (%s)", partLabel[iHypo].data()), 100, 0., 10., 200, -10., 10.);
        fOutput->Add(fHistNsigmaTOFSelBach[iHypo]);
        fOutput->Add(fHistNsigmaTPCSelBach[iHypo]);
    }
    fHistBDTOutputScore[0] = new TH1F("fHistBDTOutputScoreBkg", ";ML output score for bkg;counts", 1000, 0., 1.);
    fHistBDTOutputScore[1] = new TH1F("fHistBDTOutputScorePrompt", ";ML output score for prompt D;counts", 1000, 0., 1.);
    fHistBDTOutputScore[2] = new TH1F("fHistBDTOutputScoreNonPrompt", ";ML output score for nonprompt D;counts", 1000, 0., 1.);
    fOutput->Add(fHistBDTOutputScore[0]);
    fOutput->Add(fHistBDTOutputScore[1]);
    fOutput->Add(fHistBDTOutputScore[2]);
    float minMass = -1.;
    float maxMass = -1.;
    switch(fDecChannel)
    {
        case kD0toKpi: {
            minMass = 1.7;
            maxMass = 2.1;
            break;
        }
        case kDplustoKpipi: {
            minMass = 1.7;
            maxMass = 2.1;
            break;
        }
        case kDstartoD0pi: {
            minMass = 0.14;
            maxMass = 0.20;
            break;
        }
    }

    std::array<std::string, kNumV0IDs> v0Label = {"Kz", "La"};
    for (int iHypo{0}; iHypo<kNumV0IDs; ++iHypo) {
        auto massV0 = TDatabasePDG::Instance()->GetParticle(kPdgV0IDs[iHypo])->Mass();
        fHistMassSelV0[iHypo] = new TH2F(Form("fHistMass%sSel", v0Label[iHypo].data()), ";#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2})", 100, 0., 50., 100, massV0-0.1, massV0+0.1);
        fOutput->Add(fHistMassSelV0[iHypo]);
    }

    fInvMassVsPt = new TH2F("fInvMassVsPt", ";#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2})", 100, 0., 50., 200., minMass, maxMass);
    fOutput->Add(fInvMassVsPt);

    //Counter for Normalization
    fCounter = new AliNormalizationCounter("NormalizationCounter");
    fCounter->Init();
    PostData(3, fCounter);

    //Loading of ML models
    if(fApplyML) {
        if(!fDependOnMLSelector)
        {
            switch (fDecChannel)
            {
                case kD0toKpi:
                    fMLResponse = new AliHFMLResponseD0toKpi("D0toKpiMLResponse", "D0toKpiMLResponse", fConfigPath.data());
                    break;
                case kDplustoKpipi:
                    fMLResponse = new AliHFMLResponseDplustoKpipi("DplustoKpipiMLResponse", "DplustoKpipiMLResponse", fConfigPath.data());
                    break;
                case kDstartoD0pi:
                    fMLResponse = new AliHFMLResponseDstartoD0pi("DstartoD0piMLResponse", "DstartoD0piMLResponse", fConfigPath.data());
                    break;
            }
            fMLResponse->MLResponseInit();
        }
        else {
            std::string configLocalPath = AliMLModelHandler::ImportFile(fConfigPath.data());
            YAML::Node nodeList;
            try
            {
                nodeList = YAML::LoadFile(configLocalPath);
            }
            catch (std::exception &e)
            {
                AliFatal(Form("Yaml-ccp error: %s! Exit", e.what()));
            }
            fPtLimsML = nodeList["BINS"].as<vector<float> >();

            for (const auto &model : nodeList["MODELS"])
            {
                fMLScoreCuts.push_back(model["cut"].as<std::vector<double> >());
                fMLOptScoreCuts.push_back(model["cut_opt"].as<std::vector<std::string> >());
            }
        }
    }

    // Creat MC gen acc histos
    if (fReadMC) {
        switch (fDecChannel) {
            case kDplustoKpipi:
            {
                std::array<std::vector<int>, 2> pdgReso = {std::vector<int>{435, 10433}, std::vector<int>{}}; 
                std::set<int> pdgResoAllDecays{};
                for (auto &array: pdgReso) {
                    pdgResoAllDecays.insert(array.begin(), array.end());
                }
                for (auto &pdg: pdgResoAllDecays) {
                    fHistMCGenAccAllDecaysPrompt.push_back(new TH2F(Form("hPromptMCGenPtVsY_%d", pdg), Form("%d all decays ;#it{p}_{T} (GeV/#it{c});#it{y}", pdg), 100, 0., 50., 100., -1., 1.));
                    fHistMCGenAccAllDecaysNonPrompt.push_back(new TH2F(Form("hNonPromptMCGenPtVsY_%d", pdg), Form("%d all decays ;#it{p}_{T} (GeV/#it{c});#it{y}", pdg), 100, 0., 50., 100., -1., 1.));
                }
                for (auto iV0{0u}; iV0<fEnableV0.size(); ++iV0) {
                    if (fEnableV0[iV0]) {
                        for (auto iReso{0u}; iReso<pdgReso[iV0].size(); ++iReso) {
                            fHistMCGenAccPrompt.push_back(new TH2F(Form("hPromptMCGenPtVsY_%d_to_411_%d", pdgReso[iV0][iReso], kPdgV0IDs[iV0]), Form("%d #rightarrow 411 %d;#it{p}_{T} (GeV/#it{c});#it{y}", pdgReso[iV0][iReso], kPdgV0IDs[iV0]), 100, 0., 50., 100., -1., 1.));
                            fHistMCGenAccNonPrompt.push_back(new TH2F(Form("hNonPromptMCGenPtVsY_%d_to_411_%d", pdgReso[iV0][iReso], kPdgV0IDs[iV0]), Form("%d #rightarrow 411 %d;#it{p}_{T} (GeV/#it{c});#it{y}", pdgReso[iV0][iReso], kPdgV0IDs[iV0]), 100, 0., 50., 100., -1., 1.));
                        }
                    }
                }
                break;
            }
            case kDstartoD0pi:
            {
                std::array<std::vector<int>, 2> pdgReso = {std::vector<int>{435, 10433}, std::vector<int>{}}; 
                std::set<int> pdgResoAllDecays{};
                for (auto &array: pdgReso) {
                    pdgResoAllDecays.insert(array.begin(), array.end());
                }
                for (auto &pdg: pdgResoAllDecays) {
                    fHistMCGenAccAllDecaysPrompt.push_back(new TH2F(Form("hPromptMCGenPtVsY_%d", pdg), Form("%d all decays ;#it{p}_{T} (GeV/#it{c});#it{y}", pdg), 100, 0., 50., 100., -1., 1.));
                    fHistMCGenAccAllDecaysNonPrompt.push_back(new TH2F(Form("hNonPromptMCGenPtVsY_%d", pdg), Form("%d all decays ;#it{p}_{T} (GeV/#it{c});#it{y}", pdg), 100, 0., 50., 100., -1., 1.));
                }
                for (auto iV0{0u}; iV0<fEnableV0.size(); ++iV0) {
                    if (fEnableV0[iV0]) {
                        for (auto iReso{0u}; iReso<pdgReso[iV0].size(); ++iReso) {
                            fHistMCGenAccPrompt.push_back(new TH2F(Form("hPromptMCGenPtVsY_%d_to_413_%d", pdgReso[iV0][iReso], kPdgV0IDs[iV0]), Form("%d #rightarrow 413 %d;#it{p}_{T} (GeV/#it{c});#it{y}", pdgReso[iV0][iReso], kPdgV0IDs[iV0]), 100, 0., 50., 100., -1., 1.));
                            fHistMCGenAccNonPrompt.push_back(new TH2F(Form("hNonPromptMCGenPtVsY_%d_to_413_%d", pdgReso[iV0][iReso], kPdgV0IDs[iV0]), Form("%d #rightarrow 413 %d;#it{p}_{T} (GeV/#it{c});#it{y}", pdgReso[iV0][iReso], kPdgV0IDs[iV0]), 100, 0., 50., 100., -1., 1.));
                        }
                    }
                }
                break;
            }
        }
        for (auto iHist{0u}; iHist<fHistMCGenAccPrompt.size(); ++iHist) {
            fOutput->Add(fHistMCGenAccPrompt[iHist]);
            fOutput->Add(fHistMCGenAccNonPrompt[iHist]);
        }
        for (auto iHist{0u}; iHist<fHistMCGenAccAllDecaysPrompt.size(); ++iHist) {
            fOutput->Add(fHistMCGenAccAllDecaysPrompt[iHist]);
            fOutput->Add(fHistMCGenAccAllDecaysNonPrompt[iHist]);
        }
    }

    PostData(1, fOutput);

    fHistMCGenAccPrompt.clear();
    fHistMCGenAccNonPrompt.clear();
    fHistMCGenAccAllDecaysPrompt.clear();
    fHistMCGenAccAllDecaysNonPrompt.clear();

    if (std::accumulate(fEnableBachelor.begin(), fEnableBachelor.end(), 0))
        fNtupleCharmReso = new TNtuple("fNtupleCharmReso", "fNtupleCharmReso", "delta_inv_mass_reso:pt_reso:signal_reso:inv_mass_D:pt_D:charge_D:origin_D:pt_track:charge_track:id_track:signal_track:nsigma_tpc_track:nsigma_tof_track:outputscore_bkg_D:outputscore_prompt_D:outputscore_fd_D");
    else 
        fNtupleCharmReso = new TNtuple("fNtupleCharmReso", "fNtupleCharmReso", "delta_inv_mass_reso:pt_reso:signal_reso:inv_mass_D:pt_D:charge_D:origin_D:inv_mass_v0:pt_v0:charge_v0:id_v0:signal_v0:outputscore_bkg_D:outputscore_prompt_D:outputscore_fd_D");

    PostData(4, fNtupleCharmReso);

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFResonanceBuilder::UserExec(Option_t * /*option*/)
{
    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());

    fHistNEvents->Fill(0); // all events
    if (fAODProtection >= 0)
    {
        //   Protection against different number of events in the AOD and deltaAOD
        //   In case of discrepancy the event is rejected.
        int matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
        if (matchingAODdeltaAODlevel < 0 || (matchingAODdeltaAODlevel == 0 && fAODProtection == 1))
        {
            // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
            fHistNEvents->Fill(2);
            PostData(1, fOutput);
            return;
        }
        fHistNEvents->Fill(1);
    }

    TClonesArray *arrayCand = nullptr;
    TClonesArray *arrayCandDDau = nullptr;
    if (!fAOD && AODEvent() && IsStandardAOD())
    {
        // In case there is an AOD handler writing a standard AOD, use the AOD
        // event in memory rather than the input (ESD) event.
        fAOD = dynamic_cast<AliAODEvent *>(AODEvent());
        // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
        // have to taken from the AOD event hold by the AliAODExtension
        AliAODHandler *aodHandler = dynamic_cast<AliAODHandler *>((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
        if (aodHandler->GetExtensions())
        {
            AliAODExtension *ext = dynamic_cast<AliAODExtension *>(aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root"));
            AliAODEvent *aodFromExt = ext->GetAOD();
            switch (fDecChannel)
            {
                case kD0toKpi:
                    arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("D0toKpi"));
                    break;
                case kDplustoKpipi:
                    arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Charm3Prong"));
                    break;
                case kDstartoD0pi:
                    arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Dstar"));
                    arrayCandDDau = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("D0toKpi"));
                    break;
            }
        }
    }
    else if (fAOD)
    {
        switch (fDecChannel)
        {
            case kD0toKpi:
                arrayCand = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("D0toKpi"));
                break;
            case kDplustoKpipi:
                arrayCand = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Charm3Prong"));
                break;
            case kDstartoD0pi:
                arrayCand = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Dstar"));
                arrayCandDDau = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("D0toKpi"));
                break;
        }
    }

    if (!fAOD || !arrayCand || (fDecChannel == kDstartoD0pi && !arrayCandDDau))
    {
        AliError("Candidate branch not found!\n");
        PostData(1, fOutput);
        return;
    }

    // fix for temporary bug in ESDfilter
    // the AODs with null vertex pointer didn't pass the PhysSel
    if (!fAOD->GetPrimaryVertex() || TMath::Abs(fAOD->GetMagneticField()) < 0.001)
    {
        PostData(1, fOutput);
        return;
    }

    fHistNEvents->Fill(3); // count event

    fCounter->StoreEvent(fAOD, fRDCuts, fReadMC);

    bool isEvSel = fRDCuts->IsEventSelected(fAOD);

    if (fRDCuts->IsEventRejectedDueToTrigger())
        fHistNEvents->Fill(5);
    if (fRDCuts->IsEventRejectedDueToNotRecoVertex())
        fHistNEvents->Fill(6);
    if (fRDCuts->IsEventRejectedDueToVertexContributors())
        fHistNEvents->Fill(7);
    if (fRDCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion())
        fHistNEvents->Fill(8);
    if (fRDCuts->IsEventRejectedDueToPileup())
        fHistNEvents->Fill(9);
    if (fRDCuts->IsEventRejectedDueToCentrality())
        fHistNEvents->Fill(10);
    if (TESTBIT(fRDCuts->GetEventRejectionBitMap(), AliRDHFCuts::kPhysicsSelection))
        fHistNEvents->Fill(11);

    TClonesArray *arrayMC = nullptr;
    AliAODMCHeader *mcHeader = nullptr;

    // load MC particles
    if (fReadMC)
    {
        arrayMC = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
        if (!arrayMC)
        {
            AliWarning("MC particles branch not found!");
            PostData(1, fOutput);
            return;
        }

        // load MC header
        mcHeader = dynamic_cast<AliAODMCHeader *>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
        if (!mcHeader)
        {
            AliWarning("MC header branch not found!");
            PostData(1, fOutput);
            return;
        }
        FillMCGenHistos(arrayMC, mcHeader);
    }

    if (!isEvSel)
    {
        PostData(1, fOutput);
        return;
    }

    fHistNEvents->Fill(4); // accepted event

    // check if the train includes the common ML selector for the given charm-hadron species
    AliAnalysisTaskSECharmHadronMLSelector *taskMLSelect = nullptr;
    std::vector<int> chHadIdx{};
    if(fDependOnMLSelector)
    {
        taskMLSelect = dynamic_cast<AliAnalysisTaskSECharmHadronMLSelector*>(AliAnalysisManager::GetAnalysisManager()->GetTask(fMLSelectorName.data()));
        if(!taskMLSelect)
        {
            AliFatal("ML Selector not present in train and ML models not compiled!");
            return;
        }
        chHadIdx = taskMLSelect->GetSelectedCandidates();
        fScoresFromMLSelector = taskMLSelect->GetMLSCores();
        fScoresFromMLSelectorSecond = taskMLSelect->GetMLSCoresSecond();
    }
    else
    {
        for (int iCand = 0; iCand < arrayCand->GetEntriesFast(); iCand++)
            chHadIdx.push_back(iCand);
    }

    if (chHadIdx.size() == 0) // we don't have charm hadrons
        return;

    AliAODPidHF *pidHF = fRDCuts->GetPidHF();

    // prepare vector of selected tracks
    std::vector<int> selectedTrackIndices{};
    std::vector<int> selectedTrackIds{};
    std::vector<std::array<double, kNumBachIDs>> nSigmaTPC{}, nSigmaTOF{};
    std::vector<std::array<bool, kNumBachIDs>> selectedTrackSignal{};
  
    if (std::accumulate(fEnableBachelor.begin(), fEnableBachelor.end(), 0)) {
        for (int iTrack{0}; iTrack < fAOD->GetNumberOfTracks(); ++iTrack) {
            AliAODTrack *track = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(iTrack));
            std::array<double, kNumBachIDs> nSigmaTrkTPC = {-999., -999., -999., -999.}; // nsigma TPC and TOF
            std::array<double, kNumBachIDs> nSigmaTrkTOF = {-999., -999., -999., -999.}; // nsigma TPC and TOF
            int id = IsBachelorSelected(track, pidHF, nSigmaTrkTPC, nSigmaTrkTOF);
            if (id > 0)
            {
                selectedTrackIndices.push_back(iTrack);
                selectedTrackIds.push_back(id);
                nSigmaTPC.push_back(nSigmaTrkTPC);
                nSigmaTOF.push_back(nSigmaTrkTOF);
                std::array<bool, kNumBachIDs> isSignal = {false, false, false, false};
                if (fReadMC) { // match to MC signals
                    if (track->GetLabel() >= 0) {
                        AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(arrayMC->At(track->GetLabel()));
                        int pdgBach = part->GetPdgCode();
                        for (int iHypo{0}; iHypo<kNumBachIDs; ++iHypo) {
                            if (TMath::Abs(pdgBach) == kPdgBachIDs[iHypo]) {
                                isSignal[iHypo] = true;
                                break;
                            }
                        }
                    }
                }
                selectedTrackSignal.push_back(isSignal);
            }
        }
    }

    // prepare vector of selected V0s
    std::vector<int> selectedV0Indices{};
    std::vector<int> selectedV0Ids{};
    std::vector<std::array<bool, kNumV0IDs>> selectedV0Signal{};
    std::vector<std::array<int, kNumV0IDs>> selectedV0Labels{};
    if (std::accumulate(fEnableV0.begin(), fEnableV0.end(), 0)) {
        for (int iV0{0}; iV0 < fAOD->GetNumberOfV0s(); ++iV0) {
            AliAODv0 *v0 = fAOD->GetV0(iV0);
            int id = IsV0Selected(v0);
            if (id > 0)
            {
                selectedV0Indices.push_back(iV0);
                selectedV0Ids.push_back(id);
                std::array<bool, kNumV0IDs> isSignal = {false, false};
                std::array<int, kNumV0IDs> mcLab = {-1, -1};
                if (fReadMC) { // match to MC signals
                    if (TESTBIT(id, kK0S)) {
                        int pdgDaus[2] = {211, 211};
                        mcLab[kK0S] = v0->MatchToMC(kPdgV0IDs[kK0S], arrayMC, 2, pdgDaus);
                        if (mcLab[kK0S] >= 0)
                            isSignal[kK0S] = true;
                    } 
                    if (!isSignal[kK0S] && TESTBIT(id, kLambda)) {
                        int pdgDaus[2] = {2212, 211};
                        mcLab[kLambda] = v0->MatchToMC(kPdgV0IDs[kLambda], arrayMC, 2, pdgDaus);
                        if (mcLab[kLambda] >= 0)
                            isSignal[kLambda] = true;
                    }
                }
                selectedV0Signal.push_back(isSignal);
                selectedV0Labels.push_back(mcLab);
            }
        }
    }

    // vHF object is needed to call the method that refills the missing info of the candidates
    // if they have been deleted in dAOD reconstruction phase
    // in order to reduce the size of the file
    AliAnalysisVertexingHF vHF = AliAnalysisVertexingHF();
    for (size_t iCand = 0; iCand < chHadIdx.size(); ++iCand)
    {
        AliAODRecoDecayHF *dMeson = dynamic_cast<AliAODRecoDecayHF *>(arrayCand->UncheckedAt(chHadIdx[iCand]));
        AliAODRecoDecayHF *dMesonWithVtx;
        if(fDecChannel == kDstartoD0pi)
        {
            if(dMeson->GetIsFilled()<1)
                dMesonWithVtx = dynamic_cast<AliAODRecoDecayHF *>(arrayCandDDau->UncheckedAt(dMeson->GetProngID(1)));
            else
                dMesonWithVtx = dynamic_cast<AliAODRecoCascadeHF *>(dMeson)->Get2Prong();
        }
        else
            dMesonWithVtx = dMeson;

        bool unsetVtx = false;
        bool recVtx = false;
        AliAODVertex *origOwnVtx = nullptr;
        std::vector<double> scores{}, scoresSecond{};

        int isSelected = IsCandidateSelected(dMeson, dMesonWithVtx, &vHF, unsetVtx, recVtx, origOwnVtx, pidHF, iCand, scores, scoresSecond);
        if (!isSelected)
        {
            if (unsetVtx)
                dMesonWithVtx->UnsetOwnPrimaryVtx();
            if (recVtx)
                fRDCuts->CleanOwnPrimaryVtx(dMesonWithVtx, fAOD, origOwnVtx);
            continue;
        }

        fHistNEvents->Fill(14); // candidate selected

        // get MC truth
        AliAODMCParticle *partD = nullptr;
        int labD = -1;
        int pdgCode0 = -999;
        int orig = 0;
        int pdgD0Dau[2] = {321, 211};
        int pdgDplusDau[3] = {321, 211, 211};
        int pdgDstarDau[2] = {421, 211};

        if (fReadMC)
        {
            switch (fDecChannel)
            {
                case kD0toKpi:
                    labD = (dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson))->MatchToMC(fPdgD, arrayMC, 2, pdgD0Dau);
                break;
                case kDplustoKpipi:
                    labD = (dynamic_cast<AliAODRecoDecayHF3Prong *>(dMeson))->MatchToMC(fPdgD, arrayMC, 3, pdgDplusDau);
                break;
                case kDstartoD0pi:
                    labD = (dynamic_cast<AliAODRecoCascadeHF *>(dMeson))->MatchToMC(fPdgD, 421, pdgDstarDau, pdgD0Dau, arrayMC, false);
                break;
            }

            if (labD >= 0)
            {
                partD = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labD));
                if (fDecChannel == kD0toKpi) // check if signal is reflected
                {
                    int labDau0 = dynamic_cast<AliAODTrack *>(dMeson->GetDaughter(0))->GetLabel();
                    AliAODMCParticle *dau0 = dynamic_cast<AliAODMCParticle *>(arrayMC->UncheckedAt(TMath::Abs(labDau0)));
                    pdgCode0 = std::abs(dau0->GetPdgCode());
                }
            }
            if (partD)
            {
                orig = AliVertexingHFUtils::CheckOrigin(arrayMC, partD, true);
            }
        }

        int chargeD[2] = {0, 0};
        std::array<ROOT::Math::PxPyPzMVector, 2> fourVecD{};
        double massD[2] = {-1, -1};
        double massD4Delta[2] = {-1, -1};
        switch(fDecChannel) {
            case kD0toKpi:
            {
                double massDau[2] = {TDatabasePDG::Instance()->GetParticle(211)->Mass(), TDatabasePDG::Instance()->GetParticle(321)->Mass()};
                for (int iProng=0; iProng<2; ++iProng) {
                    fourVecD[0] += ROOT::Math::PxPyPzMVector(dMeson->PxProng(iProng), dMeson->PyProng(iProng), dMeson->PzProng(iProng), massDau[iProng]);
                    fourVecD[1] += ROOT::Math::PxPyPzMVector(dMeson->PxProng(iProng), dMeson->PyProng(iProng), dMeson->PzProng(iProng), massDau[1-iProng]);
                }
                if (isSelected == 1 || isSelected == 3) {
                    massD4Delta[0] = fourVecD[0].M();
                    massD[0] = dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson)->InvMassD0();
                    chargeD[0] = 1;
                }
                if (isSelected >=2) {
                    massD4Delta[1] = fourVecD[1].M();
                    massD[1] = dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson)->InvMassD0bar();
                    chargeD[1] = -1;
                }
                break;
            }
            case kDplustoKpipi:
            {
                double massDau[3] = {TDatabasePDG::Instance()->GetParticle(211)->Mass(), TDatabasePDG::Instance()->GetParticle(321)->Mass(), TDatabasePDG::Instance()->GetParticle(211)->Mass()};
                for (int iProng=0; iProng<3; ++iProng) {
                    fourVecD[0] += ROOT::Math::PxPyPzMVector(dMeson->PxProng(iProng), dMeson->PyProng(iProng), dMeson->PzProng(iProng), massDau[iProng]);
                }
                massD4Delta[0] = fourVecD[0].M();
                massD[0] = dynamic_cast<AliAODRecoDecayHF3Prong *>(dMeson)->InvMassDplus();
                chargeD[0] = dMeson->Charge();
                break;
            }
            case kDstartoD0pi:
            {
                fourVecD[0] += ROOT::Math::PxPyPzMVector(dMeson->PxProng(0), dMeson->PyProng(0), dMeson->PzProng(0), TDatabasePDG::Instance()->GetParticle(211)->Mass());
                AliAODRecoDecayHF2Prong* dZero = nullptr;
                if(dMeson->GetIsFilled()<1)
                    dZero = dynamic_cast<AliAODRecoDecayHF2Prong *>(arrayCandDDau->UncheckedAt(dMeson->GetProngID(1)));
                else
                    dZero = dynamic_cast<AliAODRecoCascadeHF *>(dMeson)->Get2Prong();
                double massDau[2];
                if (dMeson->Charge() > 0) {
                    massDau[1] = TDatabasePDG::Instance()->GetParticle(321)->Mass();
                    massDau[0] = TDatabasePDG::Instance()->GetParticle(211)->Mass();
                }
                else {
                    massDau[0] = TDatabasePDG::Instance()->GetParticle(321)->Mass();
                    massDau[1] = TDatabasePDG::Instance()->GetParticle(211)->Mass();
                }
                for (int iProng=0; iProng<2; ++iProng) {
                    fourVecD[0] += ROOT::Math::PxPyPzMVector(dZero->PxProng(iProng), dZero->PyProng(iProng), dZero->PzProng(iProng), massDau[iProng]);
                }
                massD4Delta[0] = fourVecD[0].M();
                massD[0] = dynamic_cast<AliAODRecoCascadeHF *>(dMeson)->DeltaInvMass();
                chargeD[0] = dMeson->Charge();
                break;
            }
        }

        // loop over tracks
        for (std::size_t iTrack{0}; iTrack < selectedTrackIndices.size(); ++iTrack) {
            AliAODTrack* track = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(selectedTrackIndices[iTrack]));

            for (int iHypo{0}; iHypo<kNumBachIDs; ++iHypo) {
                if(!fEnableBachelor[iHypo])
                    continue;
                if (IsDaughterTrack(track, dMeson, arrayCandDDau, &vHF))
                    continue;
                if (!TESTBIT(selectedTrackIds[iTrack], iHypo))
                    continue;

                int signalReso = -1;
                if (fReadMC) {
                    if (orig >= 4) { // D is signal
                        if (selectedTrackSignal[iTrack][iHypo]) { // bachelor is signal
                            AliAODMCParticle *partD = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labD));
                            AliAODMCParticle *partBach = dynamic_cast<AliAODMCParticle *>(arrayMC->At(track->GetLabel()));
                            signalReso = MatchResoToMC(partD, partBach, arrayMC);
                        }
                    }
                }

                // propagate bachelor track to PV to compute invariant mass
                std::unique_ptr<AliESDtrack> trackESD(new AliESDtrack(track));
                trackESD.get()->PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), kVeryBig);

                double massBachelor = (iHypo != kDeuteron) ? TDatabasePDG::Instance()->GetParticle(kPdgBachIDs[iHypo])->Mass() : 1.87561294257;
                auto fourVecReso = ROOT::Math::PxPyPzMVector(trackESD.get()->Px(), trackESD.get()->Py(), trackESD.get()->Pz(), massBachelor);
                if (fDecChannel != kD0toKpi || ((isSelected == 1 || isSelected == 3))) {
                    fourVecReso += fourVecD[0];
                    double deltaInvMassReso = fourVecReso.M() - massD4Delta[0];
                    if (std::abs(pdgCode0) == 321)
                        orig *= -1.; //refelcted signal
                    if (IsInvMassResoSelected(deltaInvMassReso, iHypo, -1)) {
                        std::vector<float> arr4Tuple{};
                        if(fDependOnMLSelector)
                            arr4Tuple = std::vector<float>{float(deltaInvMassReso), float(fourVecReso.Pt()), float(signalReso), float(massD[0]), float(dMeson->Pt()), float(chargeD[0]), float(orig), float(track->Pt()), float(track->Charge()), float(kPdgBachIDs[iHypo]), float(selectedTrackSignal[iTrack][iHypo]), float(nSigmaTPC[iTrack][iHypo]), float(nSigmaTOF[iTrack][iHypo]), float(fScoresFromMLSelector[iCand][0]), float(fScoresFromMLSelector[iCand][1]), float(fScoresFromMLSelector[iCand][2])};
                        else
                            arr4Tuple = std::vector<float>{float(deltaInvMassReso), float(fourVecReso.Pt()), float(signalReso), float(massD[0]), float(dMeson->Pt()), float(chargeD[0]), float(orig), float(track->Pt()), float(track->Charge()), float(kPdgBachIDs[iHypo]), float(selectedTrackSignal[iTrack][iHypo]), float(nSigmaTPC[iTrack][iHypo]), float(nSigmaTOF[iTrack][iHypo]), float(scores[0]), float(scores[1]), float(scores[2])};
                        fNtupleCharmReso->Fill(arr4Tuple.data());
                    }
                }
                if (fDecChannel == kD0toKpi && isSelected >= 2) {
                    fourVecReso += fourVecD[1];
                    double deltaInvMassReso = fourVecReso.M() - massD4Delta[1];
                    if (std::abs(pdgCode0) == 211)
                        orig *= -1.; //refelcted signal
                    if (IsInvMassResoSelected(deltaInvMassReso, iHypo, -1)) {
                        std::vector<float> arr4Tuple{};
                        if(fDependOnMLSelector)
                            arr4Tuple = std::vector<float>{float(deltaInvMassReso), float(fourVecReso.Pt()), float(signalReso), float(massD[1]), float(dMeson->Pt()), float(chargeD[1]), float(orig), float(track->Pt()), float(track->Charge()), float(kPdgBachIDs[iHypo]), float(selectedTrackSignal[iTrack][iHypo]), float(nSigmaTPC[iTrack][iHypo]), float(nSigmaTOF[iTrack][iHypo]), float(fScoresFromMLSelectorSecond[iCand][0]), float(fScoresFromMLSelectorSecond[iCand][1]), float(fScoresFromMLSelectorSecond[iCand][2])};
                        else
                            arr4Tuple = std::vector<float>{float(deltaInvMassReso), float(fourVecReso.Pt()), float(signalReso), float(massD[1]), float(dMeson->Pt()), float(chargeD[1]), float(orig), float(track->Pt()), float(track->Charge()), float(kPdgBachIDs[iHypo]), float(selectedTrackSignal[iTrack][iHypo]), float(nSigmaTPC[iTrack][iHypo]), float(nSigmaTOF[iTrack][iHypo]), float(scoresSecond[0]), float(scoresSecond[1]), float(scoresSecond[2])};
                        fNtupleCharmReso->Fill(arr4Tuple.data());
                    }
                }
            }
        }

        // loop over V0s
        for (std::size_t iV0{0}; iV0 < selectedV0Indices.size(); ++iV0) {
            AliAODv0* v0 = dynamic_cast<AliAODv0 *>(fAOD->GetV0(selectedV0Indices[iV0]));
            for (int iHypo{0}; iHypo<kNumV0IDs; ++iHypo) {
                if(!fEnableV0[iHypo])
                    continue;
                if (!TESTBIT(selectedV0Ids[iV0], iHypo))
                    continue;

                // propagate V0 track to PV to compute invariant mass
                const AliVTrack *trackVV0 = dynamic_cast<const AliVTrack*>(v0);
                if (!trackVV0)
                    continue;
                std::unique_ptr<AliNeutralTrackParam> trackV0(new AliNeutralTrackParam(trackVV0));
                trackV0.get()->PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), kVeryBig);

                double massV0 = TDatabasePDG::Instance()->GetParticle(kPdgV0IDs[iHypo])->Mass();
                auto fourVecReso = ROOT::Math::PxPyPzMVector(trackV0.get()->Px(), trackV0.get()->Py(), trackV0.get()->Pz(), massV0);
                double deltaInvMassReso[2] = {0., 0.};
                bool isRefl[2] = {false, false};
                if (fDecChannel != kD0toKpi || (isSelected == 1 || isSelected ==3)) {
                    fourVecReso += fourVecD[0];
                    deltaInvMassReso[0] = fourVecReso.M() - massD4Delta[0];
                    if (std::abs(pdgCode0) == 321)
                        isRefl[0] = true;
                }
                if (fDecChannel == kD0toKpi && isSelected >= 2) {
                    fourVecReso += fourVecD[1];
                    deltaInvMassReso[1] = fourVecReso.M() - massD4Delta[1];
                    if (std::abs(pdgCode0) == 211)
                        isRefl[1] = true;
                }
                for (int iMass=0; iMass<2; ++iMass) {
                    if (IsInvMassResoSelected(deltaInvMassReso[iMass], -1, iHypo)) {
                        if (isRefl[iMass])
                            orig *= -1.;
                        double invMassV0 = -999.;
                        double chargeV0 = 0.;
                        if (iHypo == kK0S) {
                            invMassV0 = v0->MassK0Short();
                        }
                        else if (iHypo == kLambda) {
                            if (TESTBIT(selectedV0Ids[iV0], kNumV0IDs)) {
                                invMassV0 = v0->MassAntiLambda();
                                chargeV0 = -1.;
                            }
                            else {
                                invMassV0 = v0->MassLambda();
                                chargeV0 = 1.;
                            }
                        }

                        int signalReso = -1;
                        if (fReadMC) {
                            if (orig >= 4) { // D is signal
                                if (selectedV0Signal[iV0][iHypo]) { // V0 is signal
                                    AliAODMCParticle *partD = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labD));
                                    AliAODMCParticle *partV0 = dynamic_cast<AliAODMCParticle *>(arrayMC->At(selectedV0Labels[iV0][iHypo]));
                                    signalReso = MatchResoToMC(partD, partV0, arrayMC);
                                }
                            }
                        }

                        if(fDependOnMLSelector)
                            fNtupleCharmReso->Fill(deltaInvMassReso[iMass], fourVecReso.Pt(), signalReso, massD[iMass], dMeson->Pt(), chargeD[iMass], orig, invMassV0, v0->Pt(), chargeV0, kPdgV0IDs[iHypo], selectedV0Signal[iV0][iHypo], fScoresFromMLSelector[iCand][0], fScoresFromMLSelector[iCand][1], fScoresFromMLSelector[iCand][2]);
                        else
                            fNtupleCharmReso->Fill(deltaInvMassReso[iMass], fourVecReso.Pt(), signalReso, massD[iMass], dMeson->Pt(), chargeD[iMass], orig, invMassV0, v0->Pt(), chargeV0, kPdgV0IDs[iHypo], selectedV0Signal[iV0][iHypo], scores[0], scores[1], scores[2]);
                    }
                }
            }
        }

        if (unsetVtx)
            dMesonWithVtx->UnsetOwnPrimaryVtx();
        if (recVtx)
            fRDCuts->CleanOwnPrimaryVtx(dMesonWithVtx, fAOD, origOwnVtx);
    }

    PostData(1, fOutput);
    PostData(3, fCounter);
    PostData(4, fNtupleCharmReso);
}

//________________________________________________________________________
int AliAnalysisTaskSEHFResonanceBuilder::IsCandidateSelected(AliAODRecoDecayHF *&dMeson, AliAODRecoDecayHF *&dMesonWithVtx, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx, AliAODPidHF *&pidHF, std::size_t &iCand, std::vector<double> &modelPred, std::vector<double> &modelPredSecond)
{
    if (!dMeson || !dMesonWithVtx || !vHF)
        return 0;
    fHistNEvents->Fill(12);

    // Preselection to speed up task
    TObjArray arrDauTracks(3);
    int nDau = 3;
    if (fDecChannel == kD0toKpi)
        nDau = 2;

    for (int iDau = 0; iDau < nDau; iDau++)
    {
        AliAODTrack *track;
        if ((fDecChannel != kDstartoD0pi) || (iDau == 0))
            track = vHF->GetProng(fAOD, dMeson, iDau);
        else
            track = vHF->GetProng(fAOD, dMesonWithVtx, iDau-1); //D0<-D* daughters
        arrDauTracks.AddAt(track, iDau);
    }

    if (!fRDCuts->PreSelect(arrDauTracks))
    {
        return 0;
    }

    bool isSelBit = true;
    switch (fDecChannel)
    {
        case kD0toKpi:
        {
            isSelBit = dMeson->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts);
            if (!isSelBit || !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson)))
            {
                return 0;
            }
            break;
        }
        case kDplustoKpipi:
        {
            isSelBit = dMeson->HasSelectionBit(AliRDHFCuts::kDplusCuts);
            if (!isSelBit || !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF3Prong *>(dMeson)))
            {
                return 0;
            }
            break;
        }
        case kDstartoD0pi:
        {
            if (!vHF->FillRecoCasc(fAOD, dynamic_cast<AliAODRecoCascadeHF *>(dMeson), true))
            {
                return 0;
            }
            break;
        }
    }

    fHistNEvents->Fill(13);

    unsetVtx = false;
    if (!dMesonWithVtx->GetOwnPrimaryVtx())
    {
        dMesonWithVtx->SetOwnPrimaryVtx(dynamic_cast<AliAODVertex *>(fAOD->GetPrimaryVertex()));
        unsetVtx = true;
        // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
        // Pay attention if you use continue inside this loop!!!
    }

    double ptD = dMeson->Pt();
    double yD = dMeson->Y(fPdgD);

    int ptBin = fRDCuts->PtBin(ptD);
    if (ptBin < 0)
    {
        if (unsetVtx)
            dMesonWithVtx->UnsetOwnPrimaryVtx();
        return 0;
    }

    bool isFidAcc = fRDCuts->IsInFiducialAcceptance(ptD, yD);
    if (!isFidAcc)
    {
        if (unsetVtx)
            dMesonWithVtx->UnsetOwnPrimaryVtx();
        return 0;
    }

    int isSelected = fRDCuts->IsSelected(dMeson, AliRDHFCuts::kAll, fAOD);
    if (!isSelected)
    {
        if (unsetVtx)
            dMesonWithVtx->UnsetOwnPrimaryVtx();
        return 0;
    }

    recVtx = false;
    origOwnVtx = nullptr;

    if (fRDCuts->GetIsPrimaryWithoutDaughters())
    {
        if (dMesonWithVtx->GetOwnPrimaryVtx())
            origOwnVtx = new AliAODVertex(*dMesonWithVtx->GetOwnPrimaryVtx());
        if (fRDCuts->RecalcOwnPrimaryVtx(dMesonWithVtx, fAOD))
            recVtx = true;
        else
            fRDCuts->CleanOwnPrimaryVtx(dMesonWithVtx, fAOD, origOwnVtx);
    }

    float massD[2] = {-1., -1.};
    switch(fDecChannel)
    {
        case kD0toKpi:
            if (isSelected == 1 || isSelected == 3) {
                massD[0] = dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson)->InvMassD0();
            }
            if (isSelected >= 2) {
                massD[1] = dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson)->InvMassD0bar();
            }
            break;
        case kDplustoKpipi:
            massD[0] = dynamic_cast<AliAODRecoDecayHF3Prong *>(dMeson)->InvMassDplus();
            break;
        case kDstartoD0pi:
            massD[0] = dynamic_cast<AliAODRecoCascadeHF *>(dMeson)->DeltaInvMass();
            break;
    }

    if(fApplyML)
    {
        //variables for ML application
        int isMLsel = 0;
        modelPred = {};
        modelPredSecond = {};
        double ptCand = dMeson->Pt();

        if((fDecChannel == kD0toKpi && (isSelected == 1 || isSelected == 3)) || fDecChannel == kDplustoKpipi || fDecChannel == kDstartoD0pi)
        {
            if(fDependOnMLSelector)
            {
                std::vector<float>::iterator low = std::lower_bound(fPtLimsML.begin(), fPtLimsML.end(), ptCand);
                unsigned int bin = low - fPtLimsML.begin() - 1;
                if(bin < 0)
                    bin = 0;
                else if(bin > fPtLimsML.size()-2)
                    bin = fPtLimsML.size()-2;

                isMLsel += 1;
                for(auto iScore{0u}; iScore < fScoresFromMLSelector[iCand].size(); ++iScore) {
                    if((fMLOptScoreCuts[bin][iScore] == "upper" && fScoresFromMLSelector[iCand][iScore] > fMLScoreCuts[bin][iScore]) ||
                       (fMLOptScoreCuts[bin][iScore] == "lower" && fScoresFromMLSelector[iCand][iScore] < fMLScoreCuts[bin][iScore]))
                    {
                        isMLsel -= 1;
                        break;
                    }
                }
            }
            else if (fMLResponse->IsSelectedMultiClass(modelPred, dMeson, fAOD->GetMagneticField(), pidHF, 0))
            {
                isMLsel += 1;
            }

            if (isMLsel >= 1) {
                fInvMassVsPt->Fill(dMeson->Pt(), massD[0]);
                std::size_t nClasses = fDependOnMLSelector ? fScoresFromMLSelector[iCand].size() : modelPred.size();
                for(size_t iScore = 0; iScore < nClasses; iScore++) {
                    if(fDependOnMLSelector)
                        fHistBDTOutputScore[iScore]->Fill(fScoresFromMLSelector[iCand][iScore]);
                    else
                        fHistBDTOutputScore[iScore]->Fill(modelPred[iScore]);
                }
            }
        }
        if(fDecChannel == kD0toKpi && isSelected >= 2)
        {
            if(fDependOnMLSelector)
            {
                std::vector<float>::iterator low = std::lower_bound(fPtLimsML.begin(), fPtLimsML.end(), ptCand);
                int bin = low - fPtLimsML.begin() - 1;
                if(bin < 0)
                    bin = 0;
                else if(bin > fPtLimsML.size()-2)
                    bin = fPtLimsML.size()-2;

                isMLsel += 2;
                for(auto iScore{0u}; iScore < fScoresFromMLSelectorSecond[iCand].size(); ++iScore) {
                    if((fMLOptScoreCuts[bin][iScore] == "upper" && fScoresFromMLSelectorSecond[iCand][iScore] > fMLScoreCuts[bin][iScore]) ||
                       (fMLOptScoreCuts[bin][iScore] == "lower" && fScoresFromMLSelectorSecond[iCand][iScore] < fMLScoreCuts[bin][iScore]))
                    {
                        isMLsel -= 2;
                        break;
                    }
                }
            }
            else if(fMLResponse->IsSelectedMultiClass(modelPredSecond, dMeson, fAOD->GetMagneticField(), pidHF, 1)){
                isMLsel += 2;
            }

            if (isMLsel >= 2) {
                fInvMassVsPt->Fill(dMeson->Pt(), massD[1]);
                std::size_t nClasses = fDependOnMLSelector ? fScoresFromMLSelector[iCand].size() : modelPredSecond.size();
                for(size_t iScore = 0; iScore < nClasses; iScore++) {
                    if(fDependOnMLSelector)
                        fHistBDTOutputScore[iScore]->Fill(fScoresFromMLSelector[iCand][iScore]);
                    else
                        fHistBDTOutputScore[iScore]->Fill(modelPredSecond[iScore]);
                }
            }
        }
      
        if(modelPred.size() > modelPredSecond.size())
            for(int iScore=0; iScore<modelPred.size(); iScore++)
                modelPredSecond.push_back(-9999.);
        else if(modelPred.size() < modelPredSecond.size())
            for(int iScore=0; iScore<modelPredSecond.size(); iScore++)
                modelPred.push_back(-9999.);

        return isMLsel;
    }
    if (isSelected == 1 || isSelected == 3) {
        fInvMassVsPt->Fill(dMeson->Pt(), massD[0]);
    }
    else if (isSelected >= 2) {
        fInvMassVsPt->Fill(dMeson->Pt(), massD[1]);
    }


    return isSelected;
}

//________________________________________________________________________
int AliAnalysisTaskSEHFResonanceBuilder::IsBachelorSelected(AliAODTrack *&track, AliAODPidHF *&pidHF, std::array<double, kNumBachIDs>& nSigmaTPC, std::array<double, kNumBachIDs>& nSigmaTOF)
{
    int retVal = 0;

    if (!track)
        return retVal;

    if (!track->TestFilterBit(fFilterBitBachelor))
        return retVal;
    if (track->Pt() < fPtTrackMin)
        return retVal;
    if (std::abs(track->Eta()) < 0.8)
        return retVal;

    AliPID::EParticleType parthypo[kNumBachIDs] = {AliPID::kPion, AliPID::kKaon, AliPID::kProton, AliPID::kDeuteron};
    for (int iHypo{0}; iHypo<kNumBachIDs; ++iHypo)
    {
        if (!fEnableBachelor[iHypo])
            continue;

        nSigmaTPC[iHypo] = -999.;
        int okTPC = pidHF->GetnSigmaTPC(track, parthypo[iHypo], nSigmaTPC[iHypo]);
        nSigmaTOF[iHypo] = -999.;
        int okTOF = pidHF->GetnSigmaTOF(track, parthypo[iHypo], nSigmaTOF[iHypo]);
        if (((std::abs(nSigmaTPC[iHypo]) < fNsigmaBachelorTPC[iHypo]) || okTPC<0) && ((std::abs(nSigmaTOF[iHypo]) < fNsigmaBachelorTOF[iHypo]) || okTOF<0)) {
            retVal |= BIT(iHypo);
            fHistNsigmaTPCSelBach[iHypo]->Fill(track->P(), nSigmaTPC[iHypo]);
            fHistNsigmaTOFSelBach[iHypo]->Fill(track->P(), nSigmaTOF[iHypo]);
        }
    }

    return retVal;
}

//________________________________________________________________________
int AliAnalysisTaskSEHFResonanceBuilder::IsV0Selected(AliAODv0 *&v0)
{
    int retVal = 0;

    if (!v0)
        return retVal;

    if (v0->GetOnFlyStatus())
        return retVal;

    // FIXME: hard-coded cut values

    AliAODTrack *pTrack=dynamic_cast<AliAODTrack *>(v0->GetDaughter(0)); //0->Positive Daughter
    AliAODTrack *nTrack=dynamic_cast<AliAODTrack *>(v0->GetDaughter(1)); //1->Negative Daughter

    if (!pTrack || !nTrack) {
      return retVal;
    }

    if (pTrack->GetID()<0 || nTrack->GetID()<0)
        return retVal;
    if (pTrack->Charge() == nTrack->Charge())
        return retVal;


    if (( pTrack->GetTPCClusterInfo(2,1) < 70 ) || ( nTrack->GetTPCClusterInfo(2,1) ) < 70 )
        return retVal;

    //GetKinkIndex condition
    if( pTrack->GetKinkIndex(0)>0 || nTrack->GetKinkIndex(0)>0 )
        return retVal;

    //Findable clusters > 0 condition
    if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 )
        return retVal;

    //Compute ratio Crossed Rows / Findable clusters
    double lPosTrackCrossedRowsOverFindable = pTrack->GetTPCClusterInfo(2,1) / ((double)(pTrack->GetTPCNclsF()));
    double lNegTrackCrossedRowsOverFindable = nTrack->GetTPCClusterInfo(2,1) / ((double)(nTrack->GetTPCNclsF()));
    if (lPosTrackCrossedRowsOverFindable < 0.8 || lNegTrackCrossedRowsOverFindable < 0.8)
        return retVal;

    const AliVVertex* vtx = fAOD->GetPrimaryVertex();
    double posPrimVtx[3];
    vtx->GetXYZ(posPrimVtx);

    if (v0->DcaPosToPrimVertex() < 0.06 || v0->DcaNegToPrimVertex() < 0.06)
        return retVal;

    // we want only primary V0s
    if(v0->DcaV0ToPrimVertex() > 0.05)
        return retVal;

    double dca = v0->DcaV0Daughters();
    if (dca > 1.)
        return retVal;

    double rad = std::sqrt(v0->Xv()*v0->Xv() + v0->Yv()*v0->Yv());
    if (rad < 0.5 || rad > 100.)
        return retVal;

    double invMasses[kNumV0IDs] = {v0->MassK0Short(), v0->MassLambda()};
    double invMassesAntiPart[kNumV0IDs] = {-999., v0->MassAntiLambda()};
    double y[kNumV0IDs] = {v0->RapK0Short(), v0->RapLambda()};
    double cpaMin[kNumV0IDs] = {0.97, 0.99};
    for (int iHypo{0}; iHypo<kNumV0IDs; ++iHypo)
    {
        if (!fEnableV0[iHypo])
            continue;

        if (v0->CosPointingAngle(posPrimVtx) < cpaMin[iHypo])
            continue;

        if (std::abs(y[iHypo]) > 0.8)
            continue;

        auto expMass = TDatabasePDG::Instance()->GetParticle(kPdgV0IDs[iHypo])->Mass();
        bool isPart = false, isAntiPart = false;
        if (std::abs(invMasses[iHypo] - expMass) < 0.1)
            isPart = true;
        if (std::abs(invMassesAntiPart[iHypo] - expMass) < 0.1)
            isAntiPart = true;

        if (!isPart && !isAntiPart)
            continue;

        fHistMassSelV0[iHypo]->Fill(v0->Pt(), invMasses[iHypo]);
        retVal |= BIT(iHypo);
        if (isAntiPart)
            retVal |= BIT(kNumV0IDs); // N+1 V0 tells us that is selected as antiparticle
    }

    return retVal;
}

//________________________________________________________________________
bool AliAnalysisTaskSEHFResonanceBuilder::IsInvMassResoSelected(double &mass, int bachHypo, int V0hypo)
{
    if (bachHypo >= 0) {
        switch (bachHypo)
        {
            case kPion:
            {
                for (std::size_t iMass{0}; iMass<fInvMassResoPiMin.size(); ++iMass)
                {
                    if (mass > fInvMassResoPiMin[iMass] && mass < fInvMassResoPiMax[iMass])
                        return true;
                }
                break;
            }
            case kKaon:
            {
                for (std::size_t iMass{0}; iMass<fInvMassResoKaMin.size(); ++iMass)
                {
                    if (mass > fInvMassResoKaMin[iMass] && mass < fInvMassResoKaMax[iMass])
                        return true;
                }
                break;
            }
            case kProton:
            {
                for (std::size_t iMass{0}; iMass<fInvMassResoPrMin.size(); ++iMass)
                {
                    if (mass > fInvMassResoPrMin[iMass] && mass < fInvMassResoPrMax[iMass])
                        return true;
                }
                break;
            }
            case kDeuteron:
            {
                for (std::size_t iMass{0}; iMass<fInvMassResoDeMin.size(); ++iMass)
                {
                    if (mass > fInvMassResoDeMin[iMass] && mass < fInvMassResoDeMax[iMass])
                        return true;
                }
                break;
            }
        }
    }
    else if (V0hypo >= 0) {
        switch(V0hypo)
        {
            case kK0S:
            {
                for (std::size_t iMass{0}; iMass<fInvMassResoKzMin.size(); ++iMass)
                {
                    if (mass > fInvMassResoKzMin[iMass] && mass < fInvMassResoKzMax[iMass])
                        return true;
                }
                break;
            }
            case kLambda:
            {
                for (std::size_t iMass{0}; iMass<fInvMassResoLaMin.size(); ++iMass)
                {
                    if (mass > fInvMassResoLaMin[iMass] && mass < fInvMassResoLaMax[iMass])
                        return true;
                }
                break;
            }
        }
    }

    return false;
}

//________________________________________________________________________
bool AliAnalysisTaskSEHFResonanceBuilder::IsDaughterTrack(AliAODTrack *&track, AliAODRecoDecayHF *&dMeson, TClonesArray *&arrayCandDDau, AliAnalysisVertexingHF *vHF) {

    int trkIdx = track->GetID();
    switch(fDecChannel) {
        case kD0toKpi:
        {
            for (int iProng=0; iProng<2; ++iProng) {
                AliAODTrack* dauD = dynamic_cast<AliAODTrack *>(vHF->GetProng(fAOD, dMeson, dMeson->GetProngID(iProng)));
                if (trkIdx == dauD->GetID())
                    return true;
            }
            break;
        }
        case kDplustoKpipi:
        {
            for (int iProng=0; iProng<3; ++iProng) {
                AliAODTrack* dauD = dynamic_cast<AliAODTrack *>(vHF->GetProng(fAOD, dMeson, dMeson->GetProngID(iProng)));
                if (trkIdx == dauD->GetID())
                    return true;
            }
            break;
        }
        case kDstartoD0pi:
        {
            AliAODTrack* dauD = dynamic_cast<AliAODTrack *>(vHF->GetProng(fAOD, dMeson, dMeson->GetProngID(0)));
            if (trkIdx == dauD->GetID())
                return true;

            AliAODRecoDecayHF2Prong* dZero = nullptr;
            if(dMeson->GetIsFilled()<1)
                dZero = dynamic_cast<AliAODRecoDecayHF2Prong *>(arrayCandDDau->UncheckedAt(dMeson->GetProngID(1)));
            else
                dZero = dynamic_cast<AliAODRecoCascadeHF *>(dMeson)->Get2Prong();
            for (int iProng=0; iProng<2; ++iProng) {
                dauD = dynamic_cast<AliAODTrack *>(vHF->GetProng(fAOD, dZero, dZero->GetProngID(iProng)));
                if (trkIdx == dauD->GetID())
                    return true;
            }
            break;
        }
    }

    return false;
}

//________________________________________________________________________
int AliAnalysisTaskSEHFResonanceBuilder::MatchResoToMC(AliAODMCParticle *partD, AliAODMCParticle *partLight, TClonesArray* arrayMC) {
    std::vector<int> modthersD{};
    std::vector<int> modthersLight{};
    int motherD = 0;
    while(motherD >= 0) {
        motherD = partD->GetMother();
        modthersD.push_back(motherD);
    }
    int motherLight = 0;
    while(motherLight >= 0) {
        motherLight = partD->GetMother();
        modthersLight.push_back(motherLight);
    }
    std::sort(modthersD.begin(), modthersD.end());
    std::sort(modthersLight.begin(), modthersLight.end());

    double momSumDaughters[3] = {partD->Px()+partLight->Px(), partD->Py()+partLight->Py(), partD->Pz()+partLight->Pz()};
    std::vector<int> commonMothers{};
    std::set_intersection(modthersD.begin(), modthersD.end(), modthersLight.begin(), modthersLight.end(), std::back_inserter(commonMothers));
    for (auto &mother: commonMothers) {
        AliAODMCParticle *partMother = dynamic_cast<AliAODMCParticle *>(arrayMC->At(mother));
        int pdgMother = partMother->GetPdgCode();
        if ((std::abs(pdgMother)/100 == 4) || (std::abs(pdgMother)/1000 == 4) || ((std::abs(pdgMother)-10000)/100 == 4) || ((std::abs(pdgMother)-20000)/100 == 4)) { // we are interested in charm resonances
            // let's check also momentum conservation
            double momMother[3] = {partMother->Px(), partMother->Py(), partMother->Pz()};
            bool isMomConserved = true;
            for (int iEl{0}; iEl<3; ++iEl) {
                if (std::abs(momMother[iEl]-momSumDaughters[iEl]) / (std::abs(momMother[iEl]) + 1.e-13) > 0.00001) {
                    isMomConserved = false;
                    break;
                }
            }
            if (isMomConserved) {
                return pdgMother;
            }
        }
    }

    return -1;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFResonanceBuilder::FillMCGenHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader) {

    std::array<std::vector<int>, 2> pdgReso{};
    switch(fDecChannel) {
        case kDplustoKpipi:
        {
            pdgReso[kK0S] = std::vector<int>{435, 10433};
            pdgReso[kLambda] =  std::vector<int>{};
            break;
        }
        case kDstartoD0pi:
        {
            pdgReso[kK0S] = std::vector<int>{435, 10433};
            pdgReso[kLambda] =  std::vector<int>{};
            break;
        }
        case kD0toKpi:
        {
            return;
            break;
        }
    }

    std::set<int> pdgResoAllDecays{};
    for (auto &array: pdgReso) {
        pdgResoAllDecays.insert(array.begin(), array.end());
    }

    double zMCVertex = mcHeader->GetVtxZ(); // vertex MC
    if (TMath::Abs(zMCVertex) <= fRDCuts->GetMaxVtxZ()) {
        for (int iPart{0}; iPart < arrayMC->GetEntriesFast(); ++iPart) {
            AliAODMCParticle *mcPart = dynamic_cast<AliAODMCParticle *>(arrayMC->At(iPart));
            std::array<int, 2> decay{}; 
            bool isGoodReso = false;
            int pdgPart = mcPart->GetPdgCode();
            int labDau[5] = {-1, -1, -1, -1, -1};

            for (auto &pdg: pdgResoAllDecays) {
                if (TMath::Abs(pdgPart) == pdg) {
                    isGoodReso = true;
                    break;
                }
            }

            if (isGoodReso) {
                if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iPart, mcHeader, arrayMC))
                    continue;

                int orig = AliVertexingHFUtils::CheckOrigin(arrayMC, mcPart, true);
                double pt = mcPart->Pt();
                double rapid = mcPart->Y();

                if (orig < 4)
                    continue;

                if (orig == 4) {
                    dynamic_cast<TH2F*>(fOutput->FindObject(Form("hPromptMCGenPtVsY_%d", TMath::Abs(pdgPart))))->Fill(pt, rapid);
                }
                else if (orig == 5) {
                    dynamic_cast<TH2F*>(fOutput->FindObject(Form("hNonPromptMCGenPtVsY_%d", TMath::Abs(pdgPart))))->Fill(pt, rapid);
                }

                switch(fDecChannel) {
                    case kDplustoKpipi:
                    {
                        if (fEnableV0[kK0S]) {
                            decay[kK0S] = AliVertexingHFUtils::CheckResoToDplusK0SDecay(arrayMC, mcPart, labDau);
                        }
                        decay[kLambda] = -1; //TODO: implement check for Xic* -> D+Lambda
                        break;
                    }
                    case kDstartoD0pi:
                    {
                        if (fEnableV0[kK0S]) {
                            decay[kK0S] = AliVertexingHFUtils::CheckResoToDstarK0SDecay(arrayMC, mcPart, labDau);
                        }
                        decay[kLambda] = -1; //TODO: implement check for Xic* -> D*Lambda
                        break;
                    }
                }

                int whichV0 = -1;
                for (auto iV0{0u}; iV0<fEnableV0.size(); ++iV0) {
                    if (fEnableV0[iV0] && decay[iV0]) {
                        whichV0 = iV0;
                    }
                }

                if (decay[whichV0] > 0) {
                    if (orig == 4) {
                        dynamic_cast<TH2F*>(fOutput->FindObject(Form("hPromptMCGenPtVsY_%d_to_%d_%d", TMath::Abs(decay[whichV0]), fPdgD, kPdgV0IDs[whichV0])))->Fill(pt, rapid);
                    }
                    else if (orig == 5) {
                        dynamic_cast<TH2F*>(fOutput->FindObject(Form("hNonPromptMCGenPtVsY_%d_to_%d_%d", TMath::Abs(decay[whichV0]), fPdgD, kPdgV0IDs[whichV0])))->Fill(pt, rapid);
                    }
                }
            }
        }
    }
}
