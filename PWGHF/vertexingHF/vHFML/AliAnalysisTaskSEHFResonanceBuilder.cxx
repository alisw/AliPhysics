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
#include "AliAODMCParticle.h"
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

    /// Create the output container
    //

    // Several histograms are more conveniently managed in a TList
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos");

    fHistNEvents = new TH1F("fHistNEvents", "number of events ", 16, -0.5, 15.5);
    fHistNEvents->GetXaxis()->SetBinLabel(1, "nEventsRead");
    fHistNEvents->GetXaxis()->SetBinLabel(2, "nEvents Matched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(3, "nEvents Mismatched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(4, "nEventsAnal");
    fHistNEvents->GetXaxis()->SetBinLabel(5, "n. passing IsEvSelected");
    fHistNEvents->GetXaxis()->SetBinLabel(6, "n. rejected due to trigger");
    fHistNEvents->GetXaxis()->SetBinLabel(7, "n. rejected due to not reco vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(8, "n. rejected for contr vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(9, "n. rejected for vertex out of accept");
    fHistNEvents->GetXaxis()->SetBinLabel(10, "n. rejected for pileup events");
    fHistNEvents->GetXaxis()->SetBinLabel(11, "no. of out centrality events");
    fHistNEvents->GetXaxis()->SetBinLabel(12, "no. of D candidates");
    fHistNEvents->GetXaxis()->SetBinLabel(13, "no. of D after filtering cuts");
    fHistNEvents->GetXaxis()->SetBinLabel(14, "no. of D after selection cuts");
    fHistNEvents->GetXaxis()->SetBinLabel(15, "no. of not on-the-fly rec D");
    fHistNEvents->GetXaxis()->SetBinLabel(16, "no. of D rejected by preselect");
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
        }
        case kDplustoKpipi: {
            minMass = 1.7;
            maxMass = 2.1;
        }
        case kDstartoD0pi: {
            minMass = 0.;
            maxMass = 0.2;
        }
    }

    fInvMassVsPt = new TH2F("fInvMassVsPt", ";#it{M} (GeV/#it{c});counts", 100, 0., 50., 200., minMass, maxMass);

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

    PostData(1, fOutput);

    OpenFile(4);
    fNtupleCharmReso = new TNtuple("fNtupleCharmReso", "fNtupleCharmReso", "inv_mass_reso:pt_reso:inv_mass_D_1:inv_mass_D_2:pt_D:charge_D:origin_D:pt_track:charge_track:id_track");
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

    AliAODPidHF *pidHF = fRDCuts->GetPidHF();

    // prepare vector of selected tracks
    std::vector<int> selectedTrackIndices{};
    std::vector<int> selectedTrackIds{};
    for (int iTrack{0}; iTrack < fAOD->GetNumberOfTracks(); ++iTrack) {
        AliAODTrack *track = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(iTrack));
        int id = IsBachelorSelected(track, pidHF);
        if (id > 0)
        {
            selectedTrackIndices.push_back(iTrack);
            selectedTrackIds.push_back(id);
        }
    }

    // vHF object is needed to call the method that refills the missing info of the candidates
    // if they have been deleted in dAOD reconstruction phase
    // in order to reduce the size of the file
    AliAnalysisVertexingHF vHF = AliAnalysisVertexingHF();
    for (size_t iCand = 0; iCand < chHadIdx.size(); iCand++)
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

        int isSelected = IsCandidateSelected(dMeson, dMesonWithVtx, &vHF, unsetVtx, recVtx, origOwnVtx, pidHF, iCand);
        if (!isSelected)
        {
            if (unsetVtx)
                dMesonWithVtx->UnsetOwnPrimaryVtx();
            if (recVtx)
                fRDCuts->CleanOwnPrimaryVtx(dMesonWithVtx, fAOD, origOwnVtx);
            continue;
        }

        fHistNEvents->Fill(13); // candidate selected

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
                    pdgCode0 = TMath::Abs(dau0->GetPdgCode());
                }
            }
            if (partD)
            {
                orig = AliVertexingHFUtils::CheckOrigin(arrayMC, partD, true);
            }
        }

        float massD[2] = {-1., -1.};
        int chargeD = 0.;
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
                chargeD = vHF.GetProng(fAOD, dMeson, 0)->Charge();
                break;
            case kDstartoD0pi:
                massD[0] = dynamic_cast<AliAODRecoCascadeHF *>(dMeson)->DeltaInvMass();
                chargeD = dynamic_cast<AliAODRecoCascadeHF *>(dMeson)->GetBachelor()->Charge();
                break;
        }

        // loop over tracks
        for (std::size_t iTrack{0}; iTrack < selectedTrackIndices.size(); ++iTrack) {
            AliAODTrack* track = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(selectedTrackIndices[iTrack]));
            auto fourVecD = ROOT::Math::PxPyPzMVector(dMeson->Px(), dMeson->Py(), dMeson->Pz(), TDatabasePDG::Instance()->GetParticle(fPdgD)->Mass());
            for (int iHypo{0}; iHypo<kNumBachIDs; ++iHypo) {
                if (!TESTBIT(selectedTrackIds[iTrack], iHypo))
                    continue;
                double massBachelor = (iHypo != kDeuteron) ? TDatabasePDG::Instance()->GetParticle(kPdgBachIDs[iHypo])->Mass() : 1.87561294257;
                auto fourVecBach = ROOT::Math::PxPyPzMVector(track->Px(), track->Py(), track->Pz(), massBachelor);
                auto fourVecReso = fourVecD + fourVecBach;
                auto invMassReso = fourVecReso.M();
                if (IsInvMassResoSelected(invMassReso, iHypo)) {
                    fNtupleCharmReso->Fill(invMassReso, fourVecReso.Pt(), massD[0], massD[1], dMeson->Pt(), chargeD, orig, track->Pt(), track->Charge(), kPdgBachIDs[iHypo]);
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
int AliAnalysisTaskSEHFResonanceBuilder::IsCandidateSelected(AliAODRecoDecayHF *&dMeson, AliAODRecoDecayHF *&dMesonWithVtx, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx, AliAODPidHF *&pidHF, std::size_t &iCand)
{
    if (!dMeson || !dMesonWithVtx || !vHF)
        return 0;
    fHistNEvents->Fill(11);

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
        fHistNEvents->Fill(15);
        return 0;
    }

    bool isSelBit = true;
    switch (fDecChannel)
    {
        case kD0toKpi:
            isSelBit = dMeson->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts);
            if (!isSelBit || !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson)))
            {
                fHistNEvents->Fill(14);
                return 0;
            }
            break;
        case kDplustoKpipi:
            isSelBit = dMeson->HasSelectionBit(AliRDHFCuts::kDplusCuts);
            if (!isSelBit || !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF3Prong *>(dMeson)))
            {
                fHistNEvents->Fill(14);
                return 0;
            }
            break;
        case kDstartoD0pi:
            if (!vHF->FillRecoCasc(fAOD, dynamic_cast<AliAODRecoCascadeHF *>(dMeson), true))
            {
                fHistNEvents->Fill(14);
                return 0;
            }
            break;
    }

    fHistNEvents->Fill(12);

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
        std::vector<double> modelPred = {};
        int isMLsel = 0;
        double ptCand = dMeson->Pt();

        if((fDecChannel == kD0toKpi && (isSelected == 1 || isSelected == 3)) || fDecChannel == kDplustoKpipi || fDecChannel == kDstartoD0pi)
        {
            if(fDependOnMLSelector)
            {
                std::vector<float>::iterator low = std::lower_bound(fPtLimsML.begin(), fPtLimsML.end(), ptCand);
                int bin = low - fPtLimsML.begin() - 1;
                if(bin < 0)
                    bin = 0;
                else if(bin > fPtLimsML.size()-2)
                    bin = fPtLimsML.size()-2;

                isMLsel += 1;
                for(size_t iScore = 0; iScore < fScoresFromMLSelector[iCand].size(); iScore++) {
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
                for(size_t iScore = 0; iScore < fScoresFromMLSelectorSecond[iCand].size(); iScore++) {
                    if((fMLOptScoreCuts[bin][iScore] == "upper" && fScoresFromMLSelectorSecond[iCand][iScore] > fMLScoreCuts[bin][iScore]) ||
                       (fMLOptScoreCuts[bin][iScore] == "lower" && fScoresFromMLSelectorSecond[iCand][iScore] < fMLScoreCuts[bin][iScore]))
                    {
                        isMLsel -= 2;
                        break;
                    }
                }
            }
            else if(fMLResponse->IsSelectedMultiClass(modelPred, dMeson, fAOD->GetMagneticField(), pidHF, 1)){
                isMLsel += 2;                
            }

            if (isMLsel >= 2) {
                fInvMassVsPt->Fill(dMeson->Pt(), massD[1]);
                std::size_t nClasses = fDependOnMLSelector ? fScoresFromMLSelector[iCand].size() : modelPred.size();
                for(size_t iScore = 0; iScore < nClasses; iScore++) {
                    if(fDependOnMLSelector)
                        fHistBDTOutputScore[iScore]->Fill(fScoresFromMLSelector[iCand][iScore]);
                    else
                        fHistBDTOutputScore[iScore]->Fill(modelPred[iScore]);
                }
            }
        }

        return isMLsel;
    }

    if (isSelected == 1 || isSelected == 3)
        fInvMassVsPt->Fill(dMeson->Pt(), massD[0]);
    if (isSelected >= 2)
        fInvMassVsPt->Fill(dMeson->Pt(), massD[1]);

    return isSelected;
}

//________________________________________________________________________
int AliAnalysisTaskSEHFResonanceBuilder::IsBachelorSelected(AliAODTrack *&track, AliAODPidHF *&pidHF)
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
    for (int iHypo{0}; iHypo<kNumBachIDs; iHypo++)
    {
        double nSigmaTPC = -999.;
        int okTPC = pidHF->GetnSigmaTPC(track, parthypo[iHypo], nSigmaTPC);
        double nSigmaTOF = -999.;
        int okTOF = pidHF->GetnSigmaTOF(track, parthypo[iHypo], nSigmaTOF);
        if (((std::abs(nSigmaTPC) < fNsigmaBachelorTPC[iHypo]) || okTPC<0) && ((std::abs(nSigmaTOF) < fNsigmaBachelorTOF[iHypo]) || okTOF<0)) {
            retVal |= BIT(iHypo);
            fHistNsigmaTPCSelBach[iHypo]->Fill(track->P(), nSigmaTPC);
            fHistNsigmaTOFSelBach[iHypo]->Fill(track->P(), nSigmaTOF);
        }
    }

    return retVal;
}

//________________________________________________________________________
bool AliAnalysisTaskSEHFResonanceBuilder::IsInvMassResoSelected(double &mass, int &bachHypo)
{
    switch (bachHypo)
    {
        case kPion:
        {
            for (std::size_t iMass{0}; iMass<fInvMassResoPiMin.size(); ++iMass)
            {
                if (mass > fInvMassResoPiMin[iMass] && mass < fInvMassResoPiMax[iMass])
                    return true;   
            }
        }
        case kKaon:
        {
            for (std::size_t iMass{0}; iMass<fInvMassResoKaMin.size(); ++iMass)
            {
                if (mass > fInvMassResoKaMin[iMass] && mass < fInvMassResoKaMax[iMass])
                    return true;   
            }
        }
        case kProton:
        {
            for (std::size_t iMass{0}; iMass<fInvMassResoPrMin.size(); ++iMass)
            {
                if (mass > fInvMassResoPrMin[iMass] && mass < fInvMassResoPrMax[iMass])
                    return true;   
            }
        }
        case kDeuteron:
        {
            for (std::size_t iMass{0}; iMass<fInvMassResoDeMin.size(); ++iMass)
            {
                if (mass > fInvMassResoDeMin[iMass] && mass < fInvMassResoDeMax[iMass])
                    return true;   
            }
        }
    }

    return false;
}