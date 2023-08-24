/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSEDstarPolarization
// \brief Analysis task to perform D*+ polarization analysis
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
// S. Kundu, sourav.kundu@cern.ch
/////////////////////////////////////////////////////////////

#include "yaml-cpp/yaml.h"

#include <TRandom3.h>

#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisUtils.h"
#include "AliAODHandler.h"
#include "AliAODExtension.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliNeutralTrackParam.h"
#include "AliAnalysisTaskSECharmHadronMLSelector.h"
#include "AliHFQnVectorHandler.h"
#include "AliAnalysisTaskSEHFTenderQnVectors.h"

#include "AliAnalysisTaskSEDstarPolarization.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEDstarPolarization);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEDstarPolarization::AliAnalysisTaskSEDstarPolarization() : AliAnalysisTaskSE()
{
    /// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEDstarPolarization::AliAnalysisTaskSEDstarPolarization(const char *name, AliRDHFCuts *analysisCuts) :
    AliAnalysisTaskSE(name)
{
    /// Standard constructor
    SetAnalysisCuts(analysisCuts);

    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskSEDstarPolarization::~AliAnalysisTaskSEDstarPolarization()
{
    // Destructor
    delete fOutput;
    delete fListCuts;
    delete fRDCuts;
    if (fApplyML && fMLResponse)
        delete fMLResponse;
    if(fEsdTrackCutsSoftPi)
        delete fEsdTrackCutsSoftPi;
    if(fTrkFilterSoftPi)
        delete fTrkFilterSoftPi;
    if (fApplyTrackCutVariations) {
        for (int iTrkCut{0}; iTrkCut<4; ++iTrkCut) {
            if (fRDCutsTrackVariations[iTrkCut]) {
                delete fRDCutsTrackVariations[iTrkCut];
            }
        }
    }
}

//________________________________________________________________________
void AliAnalysisTaskSEDstarPolarization::LocalInit()
{
    // Initialization

    if (fDecChannel == kDstartoD0pi) {
        AliRDHFCutsDStartoKpipi *copycut = new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi *>(fRDCuts)));
        PostData(2, copycut);
    }
    else {
        AliRDHFCutsD0toKpi *copycut = new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi *>(fRDCuts)));
        PostData(2, copycut);
    }

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDstarPolarization::UserCreateOutputObjects()
{
    /// Create the output container
    //

    // Several histograms are more conveniently managed in a TList
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos");

    fHistNEvents = new TH1F("hNEvents", "number of events ", 16, -0.5, 15.5);
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

    fHistEvPlane[0] = new TH1F("fHistEvPlaneV0M", "event plane angle V0M", 180, 0., TMath::Pi());
    fHistEvPlane[1] = new TH1F("fHistEvPlaneV0A", "event plane angle V0A", 180, 0., TMath::Pi());
    fHistEvPlane[2] = new TH1F("fHistEvPlaneV0C", "event plane angle V0C", 180, 0., TMath::Pi());
    fOutput->Add(fHistEvPlane[0]);
    fOutput->Add(fHistEvPlane[1]);
    fOutput->Add(fHistEvPlane[2]);

    fHistEvPlaneResol[0] = new TH2F("fHistEvPlaneResolV0MTPCpos", ";centrality;cos2(#psi_{V0M}-#psi_{TPCpos})", 100, 0., 100., 220, -1.1, 1.1);
    fHistEvPlaneResol[1] = new TH2F("fHistEvPlaneResolV0MTPCneg", ";centrality;cos2(#psi_{V0M}-#psi_{TPCneg})", 100, 0., 100., 220, -1.1, 1.1);
    fHistEvPlaneResol[2] = new TH2F("fHistEvPlaneResolTPCposTPCneg", ";centrality;cos2(#psi_{TPCpos}-#psi_{TPCneg})", 100, 0., 100., 220, -1.1, 1.1);
    fOutput->Add(fHistEvPlaneResol[0]);
    fOutput->Add(fHistEvPlaneResol[1]);
    fOutput->Add(fHistEvPlaneResol[2]);

    // Sparses for efficiencies (only gen)
    if (fReadMC)
        CreateEffSparses();

    //Loading of ML models
    if (fRecomputeDstarCombinatorial)
        fDependOnMLSelector = false; // if we want to use the ML selector, we cannot recompute the combinatorial D*

    if (fApplyML) {
        if (!fDependOnMLSelector)
        {
            switch (fDecChannel)
            {
                case kDstartoD0pi:
                    fMLResponse = new AliHFMLResponseDstartoD0pi("DstartoD0piMLResponse", "DstartoD0piMLResponse", fConfigPath.data());
                    break;
                case kD0toKpi:
                    fMLResponse = new AliHFMLResponseD0toKpi("D0toKpiMLResponse", "D0toKpiMLResponse", fConfigPath.data());
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

    CreateRecoSparses();

    if (fRecomputeDstarCombinatorial) {
        fEsdTrackCutsSoftPi = new AliESDtrackCuts("AliESDtrackCuts", "default");
        fEsdTrackCutsSoftPi->SetRequireITSRefit(true);
        fEsdTrackCutsSoftPi->SetMinNClustersITS(2);
        fEsdTrackCutsSoftPi->SetMaxDCAToVertexXY(1.);  
        fEsdTrackCutsSoftPi->SetMaxDCAToVertexZ(1.);
        fEsdTrackCutsSoftPi->SetPtRange(0.1, 1.e10);
        fEsdTrackCutsSoftPi->SetEtaRange(-0.8, +0.8);
        fTrkFilterSoftPi = new AliAnalysisFilter("fTrkFilterSoftPi");
        fTrkFilterSoftPi->AddCuts(fEsdTrackCutsSoftPi);
    }

    PostData(1, fOutput);

    if (fApplyTrackCutVariations && fDecChannel == kDstartoD0pi) {
        for (int iTrkCut{0}; iTrkCut<4; ++iTrkCut) {
            fRDCutsTrackVariations[iTrkCut] = new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi *>(fRDCuts)));
        }
        fRDCutsTrackVariations[0]->GetTrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC();
        fRDCutsTrackVariations[1]->SetMinCrossedRowsTPCPtDep("120-(5/pt)");
        fRDCutsTrackVariations[2]->SetMinRatioClsOverCrossRowsTPC(0.65);
        fRDCutsTrackVariations[3]->GetTrackCutsSoftPi()->SetMinNClustersITS(4);
    }

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDstarPolarization::UserExec(Option_t * /*option*/)
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
    std::vector<AliAODRecoCascadeHF*> arrayCandRecomputed{};
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
            arrayCandDDau = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("D0toKpi"));
            if (!fRecomputeDstarCombinatorial)
                arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Dstar"));
        }
    }
    else if (fAOD)
    {
        arrayCandDDau = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("D0toKpi"));
        if (!fRecomputeDstarCombinatorial) 
            arrayCand = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Dstar"));
    }

    if (!fAOD || (!fRecomputeDstarCombinatorial && !arrayCand) || !arrayCandDDau)
    {
        AliWarning("Candidate branch not found!\n");
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

    double centrality = -999.;
    AliMultSelection *multSelection = dynamic_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
    if (multSelection)
        centrality = multSelection->GetMultiplicityPercentile("V0M");

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

        // fill MC acceptance histos
        FillMCGenAccHistos(arrayMC, mcHeader, centrality);
    }

    if (!isEvSel)
    {
        PostData(1, fOutput);
        return;
    }

    fHistNEvents->Fill(4); // accepted event

    if (fRecomputeDstarCombinatorial) {
        RecomputeDstarCombinatorial(arrayCandDDau, fAOD->GetTracks(), arrayCandRecomputed);
    }

    //Get Qn-vectors from tender task
    AliHFQnVectorHandler *HFQnVectorHandler = nullptr;
    double QnFullV0[2], QnV0A[2], QnV0C[2];
    double PsinFullV0 = -1., PsinV0A = -1., PsinV0C = -1.;
    double PsinFullTPC = -1., PsinPosTPC = -1., PsinNegTPC = -1.;

    if (fComputeQnVectors && !fReadMC) {
        bool isFromTender = false;
        AliAnalysisTaskSEHFTenderQnVectors *HFQnVectorTask = dynamic_cast<AliAnalysisTaskSEHFTenderQnVectors*>(AliAnalysisManager::GetAnalysisManager()->GetTask(fTenderTaskName.data()));
        if(HFQnVectorTask) {
            HFQnVectorHandler = HFQnVectorTask->GetQnVectorHandler();
        }

        if(HFQnVectorHandler) {
            isFromTender = true;
            if(HFQnVectorHandler->GetCalibrationsOADBFileName() != fQnCalibFileName) {
                AliWarning("OADB file name for calibrations of task and Qn-vector handler not consistent!");
                return;
            }
        }
        else { //create a new handler if not found in tender task
            AliWarning("Qn-vector tender task not found! Create a new one");
            HFQnVectorHandler = new AliHFQnVectorHandler(AliHFQnVectorHandler::kQnCalib, AliHFQnVectorHandler::kQoverQlength, 2, fQnCalibFileName);
            HFQnVectorHandler->SetAODEvent(fAOD);
            HFQnVectorHandler->ComputeCalibratedQnVectorTPC();
            HFQnVectorHandler->ComputeCalibratedQnVectorV0();
        }

        //get the unnormalised Qn-vectors --> normalisation can be done in the task
        HFQnVectorHandler->GetQnVecV0(QnFullV0, QnV0A, QnV0C);
        HFQnVectorHandler->GetEventPlaneAngleV0(PsinFullV0, PsinV0A, PsinV0C);
        HFQnVectorHandler->GetEventPlaneAngleTPC(PsinFullTPC, PsinPosTPC, PsinNegTPC);

        fHistEvPlane[0]->Fill(PsinFullV0);
        fHistEvPlane[1]->Fill(PsinV0A);
        fHistEvPlane[2]->Fill(PsinV0C);

        fHistEvPlaneResol[0]->Fill(centrality, TMath::Cos(2*GetDeltaPsiSubInRange(PsinFullV0, PsinPosTPC)));
        fHistEvPlaneResol[1]->Fill(centrality, TMath::Cos(2*GetDeltaPsiSubInRange(PsinFullV0, PsinNegTPC)));
        fHistEvPlaneResol[2]->Fill(centrality, TMath::Cos(2*GetDeltaPsiSubInRange(PsinPosTPC, PsinNegTPC)));

        if (!isFromTender)
            delete HFQnVectorHandler;
    }

    // check if the train includes the common ML selector for the given charm-hadron species
    AliAnalysisTaskSECharmHadronMLSelector *taskMLSelect = nullptr;
    std::vector<int> chHadIdx{};
    std::vector<std::vector<double> > scoresFromMLSelector{}, scoresFromMLSelectorSecond{};
    if (fDependOnMLSelector) 
    {
        taskMLSelect = dynamic_cast<AliAnalysisTaskSECharmHadronMLSelector*>(AliAnalysisManager::GetAnalysisManager()->GetTask(fMLSelectorName.data()));
        if (!taskMLSelect)
        {
            AliFatal("ML Selector not present in train and ML models not compiled!");
            return;
        }
        chHadIdx = taskMLSelect->GetSelectedCandidates();
        scoresFromMLSelector = taskMLSelect->GetMLSCores();
        scoresFromMLSelectorSecond = taskMLSelect->GetMLSCoresSecond();
    }
    else if(!fRecomputeDstarCombinatorial)
    {
        for (int iCand = 0; iCand < ((fDecChannel == kDstartoD0pi) ? arrayCand->GetEntriesFast() : arrayCandDDau->GetEntriesFast()); iCand++)
        {
            chHadIdx.push_back(iCand);
            scoresFromMLSelector.push_back({});
            scoresFromMLSelectorSecond.push_back({});
        }
    }

    // vHF object is needed to call the method that refills the missing info of the candidates
    // if they have been deleted in dAOD reconstruction phase
    // in order to reduce the size of the file
    AliAnalysisVertexingHF vHF = AliAnalysisVertexingHF();

    size_t nEntries = (fRecomputeDstarCombinatorial ? arrayCandRecomputed.size() : chHadIdx.size());
    for (size_t iCand = 0; iCand < nEntries; iCand++)
    {
        AliAODRecoDecayHF *dMeson = nullptr;
        AliAODRecoCascadeHF *dStar = nullptr;
        AliAODRecoDecayHF2Prong *dZeroDau = nullptr;
        if (!fRecomputeDstarCombinatorial) {
            if (fDecChannel == kDstartoD0pi) {
                dMeson = dynamic_cast<AliAODRecoDecayHF *>(arrayCand->UncheckedAt(chHadIdx[iCand]));
                dStar = dynamic_cast<AliAODRecoCascadeHF *>(dMeson);
                if (dMeson->GetIsFilled()<1)
                    dZeroDau = dynamic_cast<AliAODRecoDecayHF2Prong *>(arrayCandDDau->UncheckedAt(dStar->GetProngID(1)));
                else
                    dZeroDau = dynamic_cast<AliAODRecoDecayHF2Prong *>(dStar->Get2Prong());
            }
            else
                dMeson = dynamic_cast<AliAODRecoDecayHF *>(arrayCandDDau->UncheckedAt(chHadIdx[iCand]));
        } 
        else {
            dMeson = arrayCandRecomputed[iCand];
            dStar = arrayCandRecomputed[iCand];
            if (dMeson->GetIsFilled()<1)
                dZeroDau = dynamic_cast<AliAODRecoDecayHF2Prong *>(arrayCandDDau->UncheckedAt(dStar->GetProngID(1)));
            else {
                dZeroDau = dynamic_cast<AliAODRecoDecayHF2Prong *>(dStar->Get2Prong());
            }
            scoresFromMLSelector.push_back({});
            scoresFromMLSelectorSecond.push_back({});
        }

        bool unsetVtx = false;
        bool recVtx = false;
        AliAODVertex *origOwnVtx = nullptr;

        int trackCutFlags[4] = {0, 0, 0, 0};
        std::vector<double> scores{}, scoresSecond{};
        int isSelected = IsCandidateSelected(dMeson, dZeroDau, &vHF, unsetVtx, recVtx, origOwnVtx, scoresFromMLSelector[iCand], scoresFromMLSelectorSecond[iCand], scores, scoresSecond, trackCutFlags);
        if (!isSelected)
        {
            if (fDecChannel == kDstartoD0pi) {
                if (unsetVtx)
                    dZeroDau->UnsetOwnPrimaryVtx();
                if (recVtx)
                    fRDCuts->CleanOwnPrimaryVtx(dZeroDau, fAOD, origOwnVtx);
            }
            else {
                if (unsetVtx)
                    dMeson->UnsetOwnPrimaryVtx();
                if (recVtx)
                    fRDCuts->CleanOwnPrimaryVtx(dMeson, fAOD, origOwnVtx);
            }
            continue;
        }

        fHistNEvents->Fill(13); // candidate selected

        // get MC truth
        AliAODMCParticle *partD = nullptr;
        int labD = -1;
        int orig = 0;
        int pdgD0Dau[2] = {321, 211};
        int pdgDstarDau[2] = {421, 211};

        if (fReadMC)
        {
            if (fDecChannel == kDstartoD0pi)
                labD = dStar->MatchToMC(413, 421, pdgDstarDau, pdgD0Dau, arrayMC, false);
            else
                labD = dMeson->MatchToMC(421, arrayMC, 2, pdgD0Dau);

            if (labD >= 0)
                partD = dynamic_cast<AliAODMCParticle *>(arrayMC->UncheckedAt(labD));

            if (partD)
                orig = AliVertexingHFUtils::CheckOrigin(arrayMC, partD, true);
        }

        // actual analysis
        double ptCand = dMeson->Pt();
        double pCand = dMeson->P();
        double yCand = (fDecChannel == kDstartoD0pi) ? dMeson->Y(413) : dMeson->Y(421);

        // random axis to test null hypothesis
        double phiRandom = gRandom->Uniform(0., 2*TMath::Pi());
        double thetaRandom = gRandom->Uniform(0., TMath::Pi());
        ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(TMath::Sin(thetaRandom) * TMath::Cos(phiRandom), TMath::Sin(thetaRandom) * TMath::Sin(phiRandom), TMath::Cos(thetaRandom));
        ROOT::Math::XYZVector randomVecXY = ROOT::Math::XYZVector(-TMath::Sin(phiRandom), TMath::Cos(phiRandom), 0.);
        ROOT::Math::XYZVector randomVecX = ROOT::Math::XYZVector(TMath::Cos(phiRandom), TMath::Sin(phiRandom), 0.);

        if (fDecChannel == kDstartoD0pi) {
            AliAODTrack* dauPi = dynamic_cast<AliAODTrack *>(dStar->GetBachelor());
            AliAODRecoDecayHF2Prong* dauD0 = dynamic_cast<AliAODRecoDecayHF2Prong *>(dStar->Get2Prong());
            fourVecPi = ROOT::Math::PxPyPzMVector(dauPi->Px(), dauPi->Py(), dauPi->Pz(), TDatabasePDG::Instance()->GetParticle(211)->Mass());
            fourVecD0 = ROOT::Math::PxPyPzMVector(dauD0->Px(), dauD0->Py(), dauD0->Pz(), TDatabasePDG::Instance()->GetParticle(421)->Mass());
            fourVecDstar = fourVecPi + fourVecD0;

            ROOT::Math::Boost boostv12{fourVecDstar.BoostToCM()};
            fourVecPiCM = boostv12(fourVecPi);

            ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(dStar->Py() / ptCand, -dStar->Px() / ptCand, 0.);
            ROOT::Math::XYZVector helicityVec = ROOT::Math::XYZVector(dStar->Px() / pCand, dStar->Py() / pCand, dStar->Pz() / pCand);
            ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0., 0., 1.);
            ROOT::Math::XYZVector Q2VecNorm = ROOT::Math::XYZVector(QnFullV0[1], -QnFullV0[0], 0.);
            ROOT::Math::XYZVector Q2Vec = ROOT::Math::XYZVector(QnFullV0[0], QnFullV0[1], 0.);

            ROOT::Math::XYZVector threeVecPiCM = fourVecPiCM.Vect();
            ROOT::Math::XYZVector threeVecPiCMXY = ROOT::Math::XYZVector(threeVecPiCM.X(), threeVecPiCM.Y(), 0.);

            double cosThetaStarProd = TMath::Abs(normalVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
            double cosThetaStarHelicity = TMath::Abs(helicityVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
            double cosThetaStarBeam = TMath::Abs(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
            double cosThetaStarRandom = TMath::Abs(randomVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
            double cosThetaStarEvPlane = fReadMC ? TMath::Abs(randomVecXY.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2())) : TMath::Abs(Q2VecNorm.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
            double cosThetaQvector = fReadMC ? TMath::Abs(randomVecX.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2())) : TMath::Abs(Q2Vec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
            double cosPhiStar = fReadMC ? TMath::Abs(randomVecX.Dot(threeVecPiCMXY) / TMath::Sqrt(threeVecPiCMXY.Mag2())) : TMath::Abs(Q2Vec.Dot(threeVecPiCMXY) / TMath::Sqrt(threeVecPiCMXY.Mag2()));
            double thetaStarBeam = TMath::ACos(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
            double phiStarBeam = TMath::ATan2(threeVecPiCM.Y(), threeVecPiCM.X());

            double mass = dStar->DeltaInvMass();

            double deltaPhi = fReadMC ? GetPhiInRange(dMeson->Phi() - phiRandom) : GetPhiInRange(dMeson->Phi() - PsinFullV0);

            std::vector<double> var4nSparse = {mass, ptCand, yCand, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, cosThetaStarEvPlane, cosThetaStarRandom, deltaPhi, centrality, scores[0], scores[1], scores[2], cosThetaQvector, cosPhiStar};
            if (fApplyTrackCutVariations) {
                for (int iTrkCut{0}; iTrkCut<4; ++iTrkCut) {
                    var4nSparse.push_back((double)trackCutFlags[iTrkCut]);
                }
            }
            std::vector<double> var4nSparseThetaPhiStar = {mass, ptCand, thetaStarBeam, phiStarBeam};

            if (!fReadMC) {
                fnSparseReco[0]->Fill(var4nSparse.data());
                fnSparseRecoThetaPhiStar[0]->Fill(var4nSparseThetaPhiStar.data());
            }
            else
            {
                if (labD > 0) {
                    if (orig == 4) {
                        fnSparseReco[1]->Fill(var4nSparse.data());
                        fnSparseRecoThetaPhiStar[1]->Fill(var4nSparseThetaPhiStar.data());
                    }
                    else if (orig == 5) {
                        fnSparseReco[2]->Fill(var4nSparse.data());
                        fnSparseRecoThetaPhiStar[2]->Fill(var4nSparseThetaPhiStar.data());
                    }
                }
                else {
                    if(fFillBkgSparse) {
                        fnSparseReco[3]->Fill(var4nSparse.data());
                        fnSparseRecoThetaPhiStar[3]->Fill(var4nSparseThetaPhiStar.data());
                    }
                }
            }

            if (unsetVtx)
                dZeroDau->UnsetOwnPrimaryVtx();
            if (recVtx)
                fRDCuts->CleanOwnPrimaryVtx(dZeroDau, fAOD, origOwnVtx);
        }
        else {
            if (isSelected == 1 || isSelected == 3) {
                AliAODTrack* dauPi = dynamic_cast<AliAODTrack *>(dMeson->GetDaughter(0));
                AliAODTrack* dauK = dynamic_cast<AliAODTrack *>(dMeson->GetDaughter(1));
                fourVecPi = ROOT::Math::PxPyPzMVector(dauPi->Px(), dauPi->Py(), dauPi->Pz(), TDatabasePDG::Instance()->GetParticle(211)->Mass());
                fourVecD0 = ROOT::Math::PxPyPzMVector(dauK->Px(), dauK->Py(), dauK->Pz(), TDatabasePDG::Instance()->GetParticle(321)->Mass()); // it's a kaon
                fourVecDstar = fourVecPi + fourVecD0;

                ROOT::Math::Boost boostv12{fourVecDstar.BoostToCM()};
                fourVecPiCM = boostv12(fourVecPi);

                ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(dMeson->Py() / ptCand, -dMeson->Px() / ptCand, 0.);
                ROOT::Math::XYZVector helicityVec = ROOT::Math::XYZVector(dMeson->Px() / pCand, dMeson->Py() / pCand, dMeson->Pz() / pCand);
                ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0., 0., 1.);
                ROOT::Math::XYZVector Q2VecNorm = ROOT::Math::XYZVector(QnFullV0[1], -QnFullV0[0], 0.);
                ROOT::Math::XYZVector Q2Vec = ROOT::Math::XYZVector(QnFullV0[0], QnFullV0[1], 0.);

                ROOT::Math::XYZVector threeVecPiCM = fourVecPiCM.Vect();
                ROOT::Math::XYZVector threeVecPiCMXY = ROOT::Math::XYZVector(threeVecPiCM.X(), threeVecPiCM.Y(), 0.);

                double cosThetaStarProd = TMath::Abs(normalVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarHelicity = TMath::Abs(helicityVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarBeam = TMath::Abs(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarEvPlane = fReadMC ? TMath::Abs(randomVecXY.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2())) : TMath::Abs(Q2VecNorm.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaQvector = fReadMC ? TMath::Abs(randomVecX.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2())) : TMath::Abs(Q2Vec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosPhiStar = fReadMC ? TMath::Abs(randomVecX.Dot(threeVecPiCMXY) / TMath::Sqrt(threeVecPiCMXY.Mag2())) : TMath::Abs(Q2Vec.Dot(threeVecPiCMXY) / TMath::Sqrt(threeVecPiCMXY.Mag2()));
                double cosThetaStarRandom = TMath::Abs(randomVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double thetaStarBeam = TMath::ACos(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double phiStarBeam = TMath::ATan2(threeVecPiCM.Y(), threeVecPiCM.X());

                double mass = dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson)->InvMassD0();

                double deltaPhi = fReadMC ? GetPhiInRange(dMeson->Phi() - phiRandom) : GetPhiInRange(dMeson->Phi() - PsinFullV0);

                std::vector<double> var4nSparse = {mass, ptCand, yCand, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, cosThetaStarEvPlane, cosThetaStarRandom, deltaPhi, centrality, scores[0], scores[1], scores[2], cosThetaQvector, cosPhiStar};
                if (fApplyTrackCutVariations) {
                    for (int iTrkCut{0}; iTrkCut<4; ++iTrkCut) {
                        var4nSparse.push_back((double)trackCutFlags[iTrkCut]);
                    }
                }
                std::vector<double> var4nSparseThetaPhiStar = {mass, ptCand, thetaStarBeam, phiStarBeam};

                if (!fReadMC) {
                    fnSparseReco[0]->Fill(var4nSparse.data());
                    fnSparseRecoThetaPhiStar[0]->Fill(var4nSparseThetaPhiStar.data());
                }
                else
                {
                    if (labD >= 0) {
                        //check if reflected signal
                        int labDauFirst = dauPi->GetLabel();
                        AliAODMCParticle* dauFirst = nullptr;
                        if (labDauFirst >= 0)
                            dauFirst = dynamic_cast<AliAODMCParticle *>(arrayMC->UncheckedAt(labDauFirst));
                        if (dauFirst && TMath::Abs(dauFirst->GetPdgCode()) == 211) {
                            if (orig == 4) {
                                fnSparseReco[1]->Fill(var4nSparse.data());
                                fnSparseRecoThetaPhiStar[1]->Fill(var4nSparseThetaPhiStar.data());
                            }
                            else if (orig == 5) {
                                fnSparseReco[2]->Fill(var4nSparse.data());
                                fnSparseRecoThetaPhiStar[2]->Fill(var4nSparseThetaPhiStar.data());
                            }
                        }
                    }
                    else {
                        if(fFillBkgSparse) {
                            fnSparseReco[3]->Fill(var4nSparse.data());
                            fnSparseRecoThetaPhiStar[3]->Fill(var4nSparseThetaPhiStar.data());
                        }
                    }
                }
            }
            if (isSelected >= 2) {
                AliAODTrack* dauPi = dynamic_cast<AliAODTrack *>(dMeson->GetDaughter(1));
                AliAODTrack* dauK = dynamic_cast<AliAODTrack *>(dMeson->GetDaughter(0));
                fourVecPi = ROOT::Math::PxPyPzMVector(dauPi->Px(), dauPi->Py(), dauPi->Pz(), TDatabasePDG::Instance()->GetParticle(211)->Mass());
                fourVecD0 = ROOT::Math::PxPyPzMVector(dauK->Px(), dauK->Py(), dauK->Pz(), TDatabasePDG::Instance()->GetParticle(321)->Mass()); // it's a kaon
                fourVecDstar = fourVecPi + fourVecD0;

                ROOT::Math::Boost boostv12{fourVecDstar.BoostToCM()};
                fourVecPiCM = boostv12(fourVecPi);

                ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(dMeson->Py() / ptCand, -dMeson->Px() / ptCand, 0.);
                ROOT::Math::XYZVector helicityVec = ROOT::Math::XYZVector(dMeson->Px() / pCand, dMeson->Py() / pCand, dMeson->Pz() / pCand);
                ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0., 0., 1.);
                ROOT::Math::XYZVector Q2VecNorm = ROOT::Math::XYZVector(QnFullV0[1], -QnFullV0[0], 0.);
                ROOT::Math::XYZVector Q2Vec = ROOT::Math::XYZVector(QnFullV0[0], QnFullV0[1], 0.);

                ROOT::Math::XYZVector threeVecPiCM = fourVecPiCM.Vect();
                ROOT::Math::XYZVector threeVecPiCMXY = ROOT::Math::XYZVector(threeVecPiCM.X(), threeVecPiCM.Y(), 0.);

                double cosThetaStarProd = TMath::Abs(normalVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarHelicity = TMath::Abs(helicityVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarBeam = TMath::Abs(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarRandom = TMath::Abs(randomVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarEvPlane = fReadMC ? TMath::Abs(randomVecXY.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2())) : TMath::Abs(Q2VecNorm.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaQvector = fReadMC ? TMath::Abs(randomVecX.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2())) : TMath::Abs(Q2Vec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosPhiStar = fReadMC ? TMath::Abs(randomVecX.Dot(threeVecPiCMXY) / TMath::Sqrt(threeVecPiCMXY.Mag2())) : TMath::Abs(Q2Vec.Dot(threeVecPiCMXY) / TMath::Sqrt(threeVecPiCMXY.Mag2()));
                double thetaStarBeam = TMath::ACos(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double phiStarBeam = TMath::ATan2(threeVecPiCM.Y(), threeVecPiCM.X());

                double mass = dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson)->InvMassD0bar();

                double deltaPhi = fReadMC ? GetPhiInRange(dMeson->Phi() - phiRandom) : GetPhiInRange(dMeson->Phi() - PsinFullV0);

                std::vector<double> var4nSparse = {mass, ptCand, yCand, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, cosThetaStarEvPlane, cosThetaStarRandom, deltaPhi, centrality, scoresSecond[0], scoresSecond[1], scoresSecond[2], cosThetaQvector, cosPhiStar};
                if (fApplyTrackCutVariations) {
                    for (int iTrkCut{0}; iTrkCut<4; ++iTrkCut) {
                        var4nSparse.push_back((double)trackCutFlags[iTrkCut]);
                    }
                }
                std::vector<double> var4nSparseThetaPhiStar = {mass, ptCand, thetaStarBeam, phiStarBeam};

                if (!fReadMC) {
                    fnSparseReco[0]->Fill(var4nSparse.data());
                    fnSparseRecoThetaPhiStar[0]->Fill(var4nSparseThetaPhiStar.data());
                }
                else
                {
                    if (labD >= 0) {
                        //check if reflected signal
                        int labDauFirst = dauPi->GetLabel();
                        AliAODMCParticle* dauFirst = nullptr;
                        if (labDauFirst >= 0)
                            dauFirst = dynamic_cast<AliAODMCParticle *>(arrayMC->UncheckedAt(labDauFirst));
                        if (dauFirst && TMath::Abs(dauFirst->GetPdgCode()) == 211) {
                            if (orig == 4) {
                                fnSparseReco[1]->Fill(var4nSparse.data());
                                fnSparseRecoThetaPhiStar[1]->Fill(var4nSparseThetaPhiStar.data());
                            }
                            else if (orig == 5) {
                                fnSparseReco[2]->Fill(var4nSparse.data());
                                fnSparseRecoThetaPhiStar[2]->Fill(var4nSparseThetaPhiStar.data());
                            }
                        }
                    }
                    else {
                        if(fFillBkgSparse) {
                            fnSparseReco[3]->Fill(var4nSparse.data());
                            fnSparseRecoThetaPhiStar[3]->Fill(var4nSparseThetaPhiStar.data());
                        }
                    }
                }
            }
            if (unsetVtx)
                dMeson->UnsetOwnPrimaryVtx();
            if (recVtx)
                fRDCuts->CleanOwnPrimaryVtx(dMeson, fAOD, origOwnVtx);
        }
    }

    if (fRecomputeDstarCombinatorial) {
        for (auto &dStar: arrayCandRecomputed) {
            if (dStar) {
                AliAODVertex *vtxDS = (AliAODVertex*)dStar->GetSecondaryVtx();
                if(vtxDS) {
                    delete vtxDS;
                    vtxDS = nullptr;
                }
                delete dStar;
                dStar = nullptr;
            }
        }
        arrayCandRecomputed.clear();
    }

    PostData(1, fOutput);
}

//________________________________________________________________________
int AliAnalysisTaskSEDstarPolarization::IsCandidateSelected(AliAODRecoDecayHF *&d, AliAODRecoDecayHF2Prong *&dZeroDau, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx, std::vector<double> scoresFromMLSelector, std::vector<double> scoresFromMLSelectorSecond, std::vector<double> &scores, std::vector<double> &scoresSecond, int trackCutFlags[4])
{
    if (!d || (!dZeroDau && fDecChannel == kDstartoD0pi) || !vHF)
        return 0;
    fHistNEvents->Fill(11);

    AliAODRecoCascadeHF* dStar = nullptr;
    int nDau = 0;
    if (fDecChannel == kDstartoD0pi) {
        dStar = dynamic_cast<AliAODRecoCascadeHF*>(d);
        nDau = 3;
    }
    else {
        nDau = 2;
    }

    if(!fRecomputeDstarCombinatorial) {
        // Preselection to speed up task
        TObjArray arrDauTracks(nDau);

        for (int iDau = 0; iDau < nDau; iDau++) {
            AliAODTrack *track = nullptr;
            if (fDecChannel == kDstartoD0pi) {
                if (iDau == 0)
                    track = vHF->GetProng(fAOD, dStar, iDau);
                else
                    track = vHF->GetProng(fAOD, dZeroDau, iDau-1); //D0<-D* daughters
            }
            else {
                track = vHF->GetProng(fAOD, d, iDau);
            }
            arrDauTracks.AddAt(track, iDau);
        }

        if (!fRDCuts->PreSelect(arrDauTracks)) {
            fHistNEvents->Fill(15);
            return 0;
        }
        if (fDecChannel == kDstartoD0pi && !vHF->FillRecoCasc(fAOD, dStar, false)) {
            fHistNEvents->Fill(14);
            return 0;
        }
        else if (fDecChannel == kD0toKpi && !d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts) && !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF2Prong *>(d))) {
            fHistNEvents->Fill(14);
            return 0;
        }
    }

    fHistNEvents->Fill(12);

    unsetVtx = false;
    if (fDecChannel == kDstartoD0pi && !dZeroDau->GetOwnPrimaryVtx()) {
        dZeroDau->SetOwnPrimaryVtx(dynamic_cast<AliAODVertex *>(fAOD->GetPrimaryVertex()));
        unsetVtx = true;
        // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
        // Pay attention if you use continue inside this loop!!!
    }
    else if (fDecChannel == kD0toKpi && !d->GetOwnPrimaryVtx()) {
        d->SetOwnPrimaryVtx(dynamic_cast<AliAODVertex *>(fAOD->GetPrimaryVertex()));
        unsetVtx = true;
        // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
        // Pay attention if you use continue inside this loop!!!
    }

    double ptD = d->Pt();

    int ptBin = fRDCuts->PtBin(ptD);
    if (ptBin < 0) {
        if (unsetVtx) {
            if (fDecChannel == kDstartoD0pi)
                dZeroDau->UnsetOwnPrimaryVtx();
            else
                d->UnsetOwnPrimaryVtx();
        }
        return 0;
    }

    int isSelected = 1;
    if (!fRecomputeDstarCombinatorial) {
        isSelected = (fDecChannel == kDstartoD0pi) ? fRDCuts->IsSelected(dStar, AliRDHFCuts::kAll, fAOD) : fRDCuts->IsSelected(d, AliRDHFCuts::kAll, fAOD);
    }
    if (!isSelected) {
        if (unsetVtx) {
            if (fDecChannel == kDstartoD0pi)
                dZeroDau->UnsetOwnPrimaryVtx();
            else
                d->UnsetOwnPrimaryVtx();
        }
        return isSelected;
    }

    recVtx = false;
    origOwnVtx = nullptr;

    if (fRDCuts->GetIsPrimaryWithoutDaughters()) {
        if (fDecChannel == kDstartoD0pi) {
            if (dZeroDau->GetOwnPrimaryVtx())
                origOwnVtx = new AliAODVertex(*dZeroDau->GetOwnPrimaryVtx());
            if (fRDCuts->RecalcOwnPrimaryVtx(dZeroDau, fAOD))
                recVtx = true;
            else
                fRDCuts->CleanOwnPrimaryVtx(dZeroDau, fAOD, origOwnVtx);
        }
        else {
            if (d->GetOwnPrimaryVtx())
                origOwnVtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
            if (fRDCuts->RecalcOwnPrimaryVtx(d, fAOD))
                recVtx = true;
            else
                fRDCuts->CleanOwnPrimaryVtx(d, fAOD, origOwnVtx);
        }
    }

    if (fApplyTrackCutVariations) {
        for (int iTrkCut{0}; iTrkCut<4; ++iTrkCut) {
            if (fRDCutsTrackVariations[iTrkCut]->IsSelected(dStar, AliRDHFCuts::kAll, fAOD)) {
                trackCutFlags[iTrkCut] = 1;
            } else {
                trackCutFlags[iTrkCut] = 0;
            }
        }
    }

    if (!fApplyML) {
        scores.push_back(-999.);
        scores.push_back(-999.);
        scores.push_back(-999.);
        scoresSecond.push_back(-999.);
        scoresSecond.push_back(-999.);
        scoresSecond.push_back(-999.);
        return isSelected;
    }
    else {
        //variables for ML application
        std::vector<double> modelPred = {};
        int isMLsel = isSelected;
        double ptCand = d->Pt();

        AliAODPidHF *pidHF = fRDCuts->GetPidHF();

        if (fDependOnMLSelector) {
            std::vector<float>::iterator low = std::lower_bound(fPtLimsML.begin(), fPtLimsML.end(), ptCand);
            int bin = low - fPtLimsML.begin() - 1;
            if (bin < 0)
                bin = 0;
            else if (bin > fPtLimsML.size()-2)
                bin = fPtLimsML.size()-2;

            scores = scoresFromMLSelector;
            scoresSecond = scoresFromMLSelectorSecond;

            if ((fDecChannel == kD0toKpi && (isSelected == 1 || isSelected == 3)) || isSelected) {
                for(size_t iScore = 0; iScore < scoresFromMLSelector.size(); iScore++) {
                    if ((fMLOptScoreCuts[bin][iScore] == "upper" && scoresFromMLSelector[iScore] > fMLScoreCuts[bin][iScore]) ||
                       (fMLOptScoreCuts[bin][iScore] == "lower" && scoresFromMLSelector[iScore] < fMLScoreCuts[bin][iScore]))
                    {
                        isMLsel -= 1;
                        break;
                    }
                }
            }
            if (fDecChannel == kD0toKpi && isSelected >= 2) {
                for(size_t iScore = 0; iScore < scoresFromMLSelectorSecond.size(); iScore++) {
                    if ((fMLOptScoreCuts[bin][iScore] == "upper" && scoresFromMLSelectorSecond[iScore] > fMLScoreCuts[bin][iScore]) ||
                    (fMLOptScoreCuts[bin][iScore] == "lower" && scoresFromMLSelectorSecond[iScore] < fMLScoreCuts[bin][iScore]))
                    {
                        isMLsel -= 2;
                        break;
                    }
                }
            }
        }
        else {
            if ((fDecChannel == kD0toKpi && (isSelected == 1 || isSelected == 3)) || isSelected) {
                if (!fMLResponse->IsSelectedMultiClass(modelPred, d, fAOD->GetMagneticField(), pidHF, 0))
                    isMLsel -= 1;
                scores = modelPred;
            }
            if (fDecChannel == kD0toKpi && isSelected >= 2) {
                if (!fMLResponse->IsSelectedMultiClass(modelPred, d, fAOD->GetMagneticField(), pidHF, 1))
                    isMLsel -= 2;
                scoresSecond = modelPred;
            }
        }
        return isMLsel;
    }
}

//________________________________________________________________________
void AliAnalysisTaskSEDstarPolarization::FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader, double centrality)
{
    /// Fill MC histos for cuts study
    ///    - at GenLimAccStep and AccStep (if fFillAcceptanceLevel=false)
    ///    - at AccStep (if fFillAcceptanceLevel=true)

    double zMCVertex = mcHeader->GetVtxZ(); //vertex MC
    if (TMath::Abs(zMCVertex) <= fRDCuts->GetMaxVtxZ())
    {
        for (int iPart = 0; iPart < arrayMC->GetEntriesFast(); iPart++)
        {
            AliAODMCParticle *mcPart = dynamic_cast<AliAODMCParticle *>(arrayMC->UncheckedAt(iPart));
            auto pdgCode = TMath::Abs(mcPart->GetPdgCode());
            if ((fDecChannel == kDstartoD0pi && pdgCode == 413) || (fDecChannel == kD0toKpi && pdgCode == 421))
            {
                int orig = AliVertexingHFUtils::CheckOrigin(arrayMC, mcPart, true); //Prompt = 4, FeedDown = 5
                bool isParticleFromOutOfBunchPileUpEvent = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iPart, mcHeader, arrayMC);

                int deca = 0;
                bool isGoodDecay = false;
                int labDau[3] = {-1, -1, -1};
                bool isFidAcc = false;
                bool isDaugInAcc = false;
                int nDau = 0;
                if (fDecChannel == kDstartoD0pi) {
                    nDau = 3;
                    deca = AliVertexingHFUtils::CheckDstarDecay(arrayMC, mcPart, labDau);
                }
                else {
                    nDau = 2;
                    deca = AliVertexingHFUtils::CheckD0Decay(arrayMC, mcPart, labDau);
                }

                if (labDau[0] == -1)
                    continue; //protection against unfilled array of labels
                if (deca == 1)
                    isGoodDecay = true;

                if (isGoodDecay)
                {
                    double pt = mcPart->Pt();
                    double p = mcPart->P();
                    double rapid = mcPart->Y();
                    isFidAcc = fRDCuts->IsInFiducialAcceptance(pt, rapid);
                    isDaugInAcc = CheckDaugAcc(arrayMC, nDau, labDau);

                    if ((fFillAcceptanceLevel && isFidAcc && isDaugInAcc) || (!fFillAcceptanceLevel && TMath::Abs(rapid) < 1))
                    {
                        int labDauFirst = mcPart->GetDaughterFirst();
                        if (labDauFirst < 0)
                            continue;
                        AliAODMCParticle* dauFirst = dynamic_cast<AliAODMCParticle *>(arrayMC->UncheckedAt(labDauFirst));
                        fourVecDstar = ROOT::Math::PxPyPzMVector(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcPart->M());
                        fourVecPi = ROOT::Math::PxPyPzMVector(dauFirst->Px(), dauFirst->Py(), dauFirst->Pz(), dauFirst->M());

                        ROOT::Math::Boost boostv12{fourVecDstar.BoostToCM()};
                        fourVecPiCM = boostv12(fourVecPi);

                        ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(mcPart->Py() / pt, -mcPart->Px() / pt, 0.);
                        ROOT::Math::XYZVector helicityVec = ROOT::Math::XYZVector(mcPart->Px() / p, mcPart->Py() / p, mcPart->Pz() / p);
                        ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0., 0., 1.);

                        ROOT::Math::XYZVector threeVecPiCM = fourVecPiCM.Vect();
                        ROOT::Math::XYZVector threeVecPiCMXY = ROOT::Math::XYZVector(threeVecPiCM.X(), threeVecPiCM.Y(), 0.);

                        double phiRandom = gRandom->Uniform(0., 2*TMath::Pi());
                        double thetaRandom = gRandom->Uniform(0., TMath::Pi());
                        ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(TMath::Sin(thetaRandom) * TMath::Cos(phiRandom), TMath::Sin(thetaRandom) * TMath::Sin(phiRandom), TMath::Cos(thetaRandom));
                        ROOT::Math::XYZVector randomVecXY = ROOT::Math::XYZVector(-TMath::Sin(phiRandom), TMath::Cos(phiRandom), 0.);
                        ROOT::Math::XYZVector randomVecX = ROOT::Math::XYZVector(TMath::Cos(phiRandom), TMath::Sin(phiRandom), 0.);

                        double cosThetaStarProd = TMath::Abs(normalVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double cosThetaStarHelicity = TMath::Abs(helicityVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double cosThetaStarBeam = TMath::Abs(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double cosThetaStarRandom = TMath::Abs(randomVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double cosThetaStarRandomXY = TMath::Abs(randomVecXY.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double cosThetaQvector = TMath::Abs(randomVecX.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double cosPhiStar = TMath::Abs(randomVecX.Dot(threeVecPiCMXY) / TMath::Sqrt(threeVecPiCMXY.Mag2()));
                        double thetaStarBeam = TMath::ACos(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double phiStarBeam = TMath::ATan2(threeVecPiCM.Y(), threeVecPiCM.X());

                        double deltaPhi = GetPhiInRange(mcPart->Phi() - phiRandom);

                        if (orig == 4 && !isParticleFromOutOfBunchPileUpEvent)
                        {
                            double var4nSparseAcc[knVarForSparseAcc] = {pt, rapid, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, cosThetaStarRandomXY, cosThetaStarRandom, deltaPhi, centrality, cosThetaQvector, cosPhiStar};
                            double var4nSparseAccThetaPhiStar[3] = {pt, thetaStarBeam, phiStarBeam};
                            fnSparseMC[0]->Fill(var4nSparseAcc);
                            fnSparseMCThetaPhiStar[0]->Fill(var4nSparseAccThetaPhiStar);
                        }
                        else if (orig == 5 && !isParticleFromOutOfBunchPileUpEvent)
                        {
                            double var4nSparseAcc[knVarForSparseAcc] = {pt, rapid, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, cosThetaStarRandomXY, cosThetaStarRandom, deltaPhi, centrality, cosThetaQvector, cosPhiStar};
                            double var4nSparseAccThetaPhiStar[3] = {pt, thetaStarBeam, phiStarBeam};
                            fnSparseMC[1]->Fill(var4nSparseAcc);
                            fnSparseMCThetaPhiStar[1]->Fill(var4nSparseAccThetaPhiStar);
                        }
                    }
                }
            }
        }
    }
}

//________________________________________________________________________
bool AliAnalysisTaskSEDstarPolarization::CheckDaugAcc(TClonesArray *arrayMC, int nProng, int *labDau)
{
    /// check if the decay products are in the good eta and pt range

    for (int iProng = 0; iProng < nProng; iProng++)
    {
        bool isSoftPion = false;
        AliAODMCParticle *mcPartDaughter = dynamic_cast<AliAODMCParticle *>(arrayMC->UncheckedAt(labDau[iProng]));
        if (!mcPartDaughter)
            return false;

        AliAODMCParticle *mother = dynamic_cast<AliAODMCParticle *>(arrayMC->UncheckedAt(mcPartDaughter->GetMother()));
        if (mother && TMath::Abs(mother->GetPdgCode()) == 413)
            isSoftPion = true;

        double eta = mcPartDaughter->Eta();
        double pt = mcPartDaughter->Pt();
        double minPt = (!isSoftPion) ? 0.1 : 0.06;

        if (TMath::Abs(eta) > 0.9 || pt < minPt)
            return false;
    }
    return true;
}

//_________________________________________________________________________
void AliAnalysisTaskSEDstarPolarization::CreateEffSparses()
{
    /// use sparses to be able to add variables if needed (multiplicity, Zvtx, etc)

    int nPtBinsCutObj = fRDCuts->GetNPtBins();
    float *ptLims = fRDCuts->GetPtBinLimits();
    int nPtBins = (int)ptLims[nPtBinsCutObj];
    if (fUseFinPtBinsForSparse)
        nPtBins = nPtBins * 10;

    int nBinsAcc[knVarForSparseAcc] = {nPtBins, 20, 5, 5, 5, 5, 5, 90, 100, 5, 5};
    double xminAcc[knVarForSparseAcc] = {0., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    double xmaxAcc[knVarForSparseAcc] = {ptLims[nPtBinsCutObj], 1., 1., 1., 1., 1., 1., TMath::Pi(), 100., 1., 1.};

    int nBinsThetaPhiAcc[3] = {nPtBins, 100, 100};
    double xminThetaPhiAcc[3] = {0., 0., 0.};
    double xmaxThetaPhiAcc[3] = {ptLims[nPtBinsCutObj], TMath::Pi(), TMath::Pi()};

    TString label[2] = {"fromC", "fromB"};
    for (int iHist = 0; iHist < 2; iHist++)
    {
        TString titleSparse = Form("MC nSparse (%s)- %s", fFillAcceptanceLevel ? "Acc.Step" : "Gen.Acc.Step", label[iHist].Data());
        fnSparseMC[iHist] = new THnSparseF(Form("fnSparseAcc_%s", label[iHist].Data()), titleSparse.Data(), knVarForSparseAcc, nBinsAcc, xminAcc, xmaxAcc);
        fnSparseMC[iHist]->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
        fnSparseMC[iHist]->GetAxis(1)->SetTitle("#it{y}");
        fnSparseMC[iHist]->GetAxis(2)->SetTitle("|cos(#theta*)| (beam)");
        fnSparseMC[iHist]->GetAxis(3)->SetTitle("|cos(#theta*)| (production)");
        fnSparseMC[iHist]->GetAxis(4)->SetTitle("|cos(#theta*)| (helicity)");
        fnSparseMC[iHist]->GetAxis(5)->SetTitle("|cos(#theta*)| (random XY)");
        fnSparseMC[iHist]->GetAxis(6)->SetTitle("|cos(#theta*)| (random)");
        fnSparseMC[iHist]->GetAxis(7)->SetTitle("#varphi - #psi_{2}");
        fnSparseMC[iHist]->GetAxis(8)->SetTitle("centrality");
        fnSparseMC[iHist]->GetAxis(9)->SetTitle("|cos(#theta_{Q-vector})|");
        fnSparseMC[iHist]->GetAxis(10)->SetTitle("|cos(#varphi*)|");
        fOutput->Add(fnSparseMC[iHist]);

        fnSparseMCThetaPhiStar[iHist] = new THnSparseF(Form("fnSparseMCThetaPhiStar_%s", label[iHist].Data()), titleSparse.Data(), 3, nBinsThetaPhiAcc, xminThetaPhiAcc, xmaxThetaPhiAcc);
        fnSparseMCThetaPhiStar[iHist]->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
        fnSparseMCThetaPhiStar[iHist]->GetAxis(1)->SetTitle("#theta* (beam)");
        fnSparseMCThetaPhiStar[iHist]->GetAxis(2)->SetTitle("#varphi* (beam)");
        fOutput->Add(fnSparseMCThetaPhiStar[iHist]);
    }
}

//_________________________________________________________________________
void AliAnalysisTaskSEDstarPolarization::CreateRecoSparses()
{
    int nPtBinsCutObj = fRDCuts->GetNPtBins();
    float *ptLims = fRDCuts->GetPtBinLimits();
    int nPtBins = (int)ptLims[nPtBinsCutObj];
    if (fUseFinPtBinsForSparse)
        nPtBins = nPtBins * 10;

    int nMassBins = 500;
    double massMin = 0.138, massMax = 0.160;
    TString massTitle = "#it{M}(K#pi#pi) #minus #it{M}(K#pi)";
    if (fDecChannel == kD0toKpi) {
        massMin = 1.65;
        massMax = 2.15;
        massTitle = "#it{M}(K#pi)";
    }

    int nCosThetaBins = 5;

    int nBinsReco[knVarForSparseReco] = {nMassBins, nPtBins, 20, nCosThetaBins, nCosThetaBins, nCosThetaBins, nCosThetaBins, nCosThetaBins, 90, 100, fNBinsML[0], fNBinsML[1], fNBinsML[2], 5, 5, 2, 2, 2, 2};
    double xminReco[knVarForSparseReco] = {massMin, 0., -1., 0., 0., 0., 0., 0., 0., 0., fMLOutputMin[0], fMLOutputMin[1], fMLOutputMin[2], 0., 0., -0.5, -0.5, -0.5, -0.5};
    double xmaxReco[knVarForSparseReco] = {massMax, ptLims[nPtBinsCutObj], 1., 1., 1., 1., 1., 1., TMath::Pi(), 100., fMLOutputMax[0], fMLOutputMax[1], fMLOutputMax[2], 1., 1., 1.5, 1.5, 1.5, 1.5};

    int nBinsThetaPhiReco[4] = {nMassBins, nPtBins, 100, 100};
    double xminThetaPhiReco[4] = {massMin, 0., 0., 0.};
    double xmaxThetaPhiReco[4] = {massMax, ptLims[nPtBinsCutObj], TMath::Pi(), TMath::Pi()};

    TString label[4] = {"all", "fromC", "fromB", "bkg"};
    for (int iHist = 0; iHist < 4; iHist++)
    {
        int nVars = (fApplyTrackCutVariations) ? knVarForSparseReco : knVarForSparseReco - 4;
        TString titleSparse = Form("Reco nSparse - %s", label[iHist].Data());
        fnSparseReco[iHist] = new THnSparseF(Form("fnSparseReco_%s", label[iHist].Data()), titleSparse.Data(), nVars, nBinsReco, xminReco, xmaxReco);
        fnSparseReco[iHist]->GetAxis(0)->SetTitle(Form("%s (GeV/#it{c}^{2})", massTitle.Data()));
        fnSparseReco[iHist]->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
        fnSparseReco[iHist]->GetAxis(2)->SetTitle("#it{y}");
        fnSparseReco[iHist]->GetAxis(3)->SetTitle("|cos(#theta*)| (beam)");
        fnSparseReco[iHist]->GetAxis(4)->SetTitle("|cos(#theta*)| (production)");
        fnSparseReco[iHist]->GetAxis(5)->SetTitle("|cos(#theta*)| (helicity)");
        if (fReadMC)
            fnSparseReco[iHist]->GetAxis(6)->SetTitle("|cos(#theta*)| (random XY)");
        else
            fnSparseReco[iHist]->GetAxis(6)->SetTitle("|cos(#theta*)| (event plane V0M)");
        fnSparseReco[iHist]->GetAxis(7)->SetTitle("|cos(#theta*)| (random)");
        fnSparseReco[iHist]->GetAxis(8)->SetTitle("#varphi - #psi_{2}");
        fnSparseReco[iHist]->GetAxis(9)->SetTitle("centrality %");
        fnSparseReco[iHist]->GetAxis(10)->SetTitle("ML bkg output score");
        fnSparseReco[iHist]->GetAxis(11)->SetTitle("ML prompt output score");
        fnSparseReco[iHist]->GetAxis(12)->SetTitle("ML non-prompt output score");
        fnSparseReco[iHist]->GetAxis(13)->SetTitle("|cos(#theta_{Q-vector})|");
        fnSparseReco[iHist]->GetAxis(14)->SetTitle("|cos(#varphi*)|");
        if (fApplyTrackCutVariations) {
            for (int iTrkCut{0}; iTrkCut<4; ++iTrkCut) {
                fnSparseReco[iHist]->GetAxis(15 + iTrkCut)->SetTitle(Form("track cut %d", iTrkCut));
            }
        }
        fOutput->Add(fnSparseReco[iHist]);

        fnSparseRecoThetaPhiStar[iHist] = new THnSparseF(Form("fnSparseRecoThetaPhiStar_%s", label[iHist].Data()), titleSparse.Data(), 4, nBinsThetaPhiReco, xminThetaPhiReco, xmaxThetaPhiReco);
        fnSparseRecoThetaPhiStar[iHist]->GetAxis(0)->SetTitle(Form("%s (GeV/#it{c}^{2})", massTitle.Data()));
        fnSparseRecoThetaPhiStar[iHist]->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
        fnSparseRecoThetaPhiStar[iHist]->GetAxis(2)->SetTitle("#theta* (beam)");
        fnSparseRecoThetaPhiStar[iHist]->GetAxis(3)->SetTitle("#varphi* (beam)");
        fOutput->Add(fnSparseRecoThetaPhiStar[iHist]);
    }
}

//_________________________________________________________________________
bool AliAnalysisTaskSEDstarPolarization::RecomputeDstarCombinatorial(TClonesArray *array2Prongs, TClonesArray *arrayTracks, std::vector<AliAODRecoCascadeHF*> &arrayDstar) {

    AliAnalysisVertexingHF vHF = AliAnalysisVertexingHF();

    fBzkG = (double)fAOD->GetMagneticField();
    const AliVVertex *vprimary = fAOD->GetPrimaryVertex();
    double pos[3];
    double cov[6];
    vprimary->GetXYZ(pos);
    vprimary->GetCovarianceMatrix(cov);
    AliESDVertex *fV1 = new AliESDVertex(pos,cov,100.,100,vprimary->GetName());
    fV1->GetCovMatrix(cov);

    for (int iTrack=0; iTrack<arrayTracks->GetEntriesFast(); iTrack++) {
        auto track = static_cast<AliAODTrack*>(arrayTracks->At(iTrack));
        AliESDtrack *trackESD = new AliESDtrack(track);
        if (!trackESD->PropagateToDCA(fV1, fBzkG, kVeryBig)) {
            delete trackESD;
            trackESD = nullptr;
            continue;
        }
        trackESD->RelateToVertex(fV1, fBzkG, kVeryBig);
        if (!fTrkFilterSoftPi->IsSelected(trackESD)) {
            delete trackESD;
            trackESD = nullptr;
            continue;
        }

        for (int iCand=0; iCand<array2Prongs->GetEntriesFast(); iCand++) {
            auto dZero = static_cast<AliAODRecoDecayHF2Prong*>(array2Prongs->At(iCand));
            if(!vHF.FillRecoCand(fAOD, dZero) || !dZero->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts))
                continue;
            AliAODRecoCascadeHF *dStar = MakeCascade(dZero, track, trackESD, fV1);
            if (dStar) {
                if (fRDCuts->IsSelected(dStar, AliRDHFCuts::kAll, fAOD))
                    arrayDstar.push_back(dStar);
                else {
                    AliAODVertex *vtxDS = (AliAODVertex*)dStar->GetSecondaryVtx();
                        if(vtxDS) {
                            delete vtxDS;
                            vtxDS = nullptr;
                        }
                    delete dStar;
                    dStar = nullptr;
                }
            }
        }
        delete trackESD;
        trackESD = nullptr;
    }

    delete fV1;
    fV1 = nullptr;

    return true;
}

//----------------------------------------------------------------------------
AliAODRecoCascadeHF *AliAnalysisTaskSEDstarPolarization::MakeCascade(AliAODRecoDecayHF2Prong* trackD0, AliAODTrack *track, AliESDtrack *esdTrackPi, AliESDVertex *fV1) {

    double pxNoVtx[2] = {esdTrackPi->Px(), trackD0->Px()};
    double pyNoVtx[2] = {esdTrackPi->Py(), trackD0->Py()};
    double pzNoVtx[2] = {esdTrackPi->Pz(), trackD0->Pz()};
    if (!SelectInvMassAndPtDstarD0pi(pxNoVtx, pyNoVtx, pzNoVtx))
       return nullptr;

    TObjArray* twoTrackArrayCasc = new TObjArray(2);
    AliNeutralTrackParam *trackV0 = new AliNeutralTrackParam(trackD0);

    twoTrackArrayCasc->AddAt(esdTrackPi, 0);
    twoTrackArrayCasc->AddAt(trackV0, 1);

    const AliVVertex *vprimary = fAOD->GetPrimaryVertex();
    double pos[3];
    double cov[6];
    fV1->GetXYZ(pos);
    fV1->GetCovMatrix(cov);

    double dca = 0.;
    AliAODVertex *vtxCasc = nullptr;
    double chi2perNDF = fV1->GetChi2toNDF();
    vtxCasc = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,2);
    if(!vtxCasc) {
        twoTrackArrayCasc->Clear();
        twoTrackArrayCasc->Delete();
        delete twoTrackArrayCasc;
        delete vtxCasc;
        vtxCasc=nullptr;
        delete trackV0;
        trackV0=nullptr;
        return nullptr;
    }

    double px[2],py[2],pz[2],d0[2],d0err[2];
    // propagate tracks to secondary vertex, to compute inv. mass
    esdTrackPi->PropagateToDCA(vtxCasc,fBzkG,kVeryBig);
    trackV0->PropagateToDCA(vtxCasc,fBzkG,kVeryBig);
    double momentum[3];
    esdTrackPi->GetPxPyPz(momentum);
    px[0] = momentum[0]; py[0] = momentum[1]; pz[0] = momentum[2];
    trackV0->GetPxPyPz(momentum);
    px[1] = momentum[0]; py[1] = momentum[1]; pz[1] = momentum[2];

    double d0z0[2],covd0z0[3];
    esdTrackPi->PropagateToDCA(vprimary,fBzkG,kVeryBig,d0z0,covd0z0);
    d0[0] = d0z0[0];
    d0err[0] = TMath::Sqrt(covd0z0[0]);
    trackV0->PropagateToDCA(vprimary,fBzkG,kVeryBig,d0z0,covd0z0);
    d0[1] = d0z0[0];
    d0err[1] = TMath::Sqrt(covd0z0[0]);
    AliAODRecoCascadeHF *rCasc = new AliAODRecoCascadeHF(0x0, 1, px, py, pz, d0, d0err, 0.);
    rCasc->SetPxPyPzProngs(2,px,py,pz);
    rCasc->SetDCA(dca);
    rCasc->Setd0Prongs(2,d0);
    rCasc->Setd0errProngs(2,d0err);
    rCasc->SetOwnPrimaryVtx((AliAODVertex*)fAOD->GetPrimaryVertex());
    vtxCasc->SetParent(rCasc);
    rCasc->SetSecondaryVtx(vtxCasc);
    vtxCasc->AddDaughter(track);
    vtxCasc->AddDaughter(trackD0);
    rCasc->SetPrimaryVtxRef((AliAODVertex*)fAOD->GetPrimaryVertex());
    rCasc->SetCharge(esdTrackPi->Charge());
    // get PID info from ESD
    double esdpid0[5]={0.,0.,0.,0.,0.};
    if(esdTrackPi->GetStatus()&AliESDtrack::kESDpid) esdTrackPi->GetESDpid(esdpid0);
    double esdpid1[5]={0.,0.,0.,0.,0.};
    double esdpid[10];
    for(int i=0;i<5;i++) {
        esdpid[i]   = esdpid0[i];
        esdpid[5+i] = esdpid1[i];
    }
    rCasc->SetPID(2,esdpid);
    rCasc->SetIsFilled(2);

    twoTrackArrayCasc->Clear();
    twoTrackArrayCasc->Delete();
    delete twoTrackArrayCasc;
    delete trackV0;
    trackV0=nullptr;

    return rCasc;
}

//-----------------------------------------------------------------------------
bool AliAnalysisTaskSEDstarPolarization::SelectInvMassAndPtDstarD0pi(double *px, double *py, double *pz){
  unsigned int pdg2[2];
  int nprongs=2;
  double minv2,mrange;
  double lolim,hilim;
  double minPt=0;
  bool retval=false;

  double d02[2]={0.,0.};
  AliAODRecoDecay fMassCalc2(0x0,2,0,d02);
  fMassCalc2.SetPxPyPzProngs(nprongs,px,py,pz);
  // pt cut
  double ptcand = TMath::Sqrt(fMassCalc2.Pt2());
  minPt = fRDCuts->GetMinPtCandidate();
  if(minPt > 0.1 && ptcand < minPt)
    return retval;

  // mass cut
  int jPtBinStar=fRDCuts->PtBin(ptcand);
  if(jPtBinStar<0) jPtBinStar=0;
  mrange=((AliRDHFCutsDStartoKpipi*)fRDCuts)->GetMassCut(jPtBinStar);
  double massDstar=TDatabasePDG::Instance()->GetParticle(413)->Mass();
  lolim=massDstar-mrange;
  hilim=massDstar+mrange;
  pdg2[0]=211; pdg2[1]=421; // in twoTrackArrayCasc we put the pion first
  minv2 = fMassCalc2.InvMass2(nprongs, pdg2);
  if(minv2>lolim*lolim && minv2<hilim*hilim ){
    retval=true;
  }

  return retval;
}

//________________________________________________________________________
double AliAnalysisTaskSEDstarPolarization::GetPhiInRange(double phi)
{
    // Sets the phi angle in the range [0,2*pi/harmonic]

    double result = phi;
    while(result < 0) {
        result = result + 2. * TMath::Pi() / 2;
    }
    while(result > 2.*TMath::Pi() / 2){
        result = result - 2. * TMath::Pi() / 2;
    }
    return result;
}

//________________________________________________________________________
double AliAnalysisTaskSEDstarPolarization::GetDeltaPsiSubInRange(double psi1, double psi2)
{
    // difference of subevents reaction plane angle cannot be bigger than pi / n

    double delta = psi1 - psi2;
    if(TMath::Abs(delta) > TMath::Pi() / 2) {
        if(delta>0.) delta -= 2.*TMath::Pi() / 2;
        else delta += 2.*TMath::Pi() / 2;
    }

    return delta;
}
