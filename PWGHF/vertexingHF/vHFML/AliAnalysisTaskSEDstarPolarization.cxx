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
#include "AliAnalysisTaskSECharmHadronMLSelector.h"

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
    if(fApplyML && fMLResponse)
        delete fMLResponse;
}

//________________________________________________________________________
void AliAnalysisTaskSEDstarPolarization::LocalInit()
{
    // Initialization

    if(fDecChannel == kDstartoD0pi) {
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

    // Sparses for efficiencies (only gen)
    if(fReadMC)
        CreateEffSparses();

    //Loading of ML models
    if(fApplyML) {
        if(!fDependOnMLSelector)
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

    PostData(1, fOutput);

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
            arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Dstar"));
            arrayCandDDau = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("D0toKpi"));
        }
    }
    else if (fAOD)
    {
        arrayCand = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Dstar"));
        arrayCandDDau = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("D0toKpi"));
    }

    if (!fAOD || !arrayCand || !arrayCandDDau)
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
    if(multSelection)
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

    // check if the train includes the common ML selector for the given charm-hadron species
    AliAnalysisTaskSECharmHadronMLSelector *taskMLSelect = nullptr;
    std::vector<int> chHadIdx{};
    std::vector<std::vector<double> > scoresFromMLSelector{}, scoresFromMLSelectorSecond{};
    if(fDependOnMLSelector) 
    {
        taskMLSelect = dynamic_cast<AliAnalysisTaskSECharmHadronMLSelector*>(AliAnalysisManager::GetAnalysisManager()->GetTask(fMLSelectorName.data()));
        if(!taskMLSelect)
        {
            AliFatal("ML Selector not present in train and ML models not compiled!");
            return;
        }
        chHadIdx = taskMLSelect->GetSelectedCandidates();
        scoresFromMLSelector = taskMLSelect->GetMLSCores();
        scoresFromMLSelectorSecond = taskMLSelect->GetMLSCoresSecond();
    }
    else
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

    for (size_t iCand = 0; iCand < chHadIdx.size(); iCand++)
    {
        AliAODRecoDecayHF *dMeson = nullptr;
        AliAODRecoCascadeHF *dStar = nullptr;
        AliAODRecoDecayHF2Prong *dZeroDau = nullptr;
        if(fDecChannel == kDstartoD0pi) {
            dMeson = dynamic_cast<AliAODRecoDecayHF *>(arrayCand->UncheckedAt(chHadIdx[iCand]));
            dStar = dynamic_cast<AliAODRecoCascadeHF *>(dMeson);
            if(dMeson->GetIsFilled()<1)
                dZeroDau = dynamic_cast<AliAODRecoDecayHF2Prong *>(arrayCandDDau->UncheckedAt(dStar->GetProngID(1)));
            else
                dZeroDau = dynamic_cast<AliAODRecoDecayHF2Prong *>(dStar->Get2Prong());
        }
        else
            dMeson = dynamic_cast<AliAODRecoDecayHF *>(arrayCandDDau->UncheckedAt(chHadIdx[iCand]));

        bool unsetVtx = false;
        bool recVtx = false;
        AliAODVertex *origOwnVtx = nullptr;

        int isSelected = IsCandidateSelected(dMeson, dZeroDau, &vHF, unsetVtx, recVtx, origOwnVtx, scoresFromMLSelector[iCand], scoresFromMLSelectorSecond[iCand]);
        if (!isSelected)
        {
            if(fDecChannel == kDstartoD0pi) {
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
            if(fDecChannel == kDstartoD0pi)
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

        if(fDecChannel == kDstartoD0pi) {
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

            ROOT::Math::XYZVector threeVecPiCM = fourVecPiCM.Vect();

            double cosThetaStarProd = TMath::Abs(normalVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
            double cosThetaStarHelicity = TMath::Abs(helicityVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
            double cosThetaStarBeam = TMath::Abs(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
            double cosThetaStarRandom = TMath::Abs(randomVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
            double thetaStarBeam = TMath::ACos(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
            double phiStarBeam = TMath::ATan2(threeVecPiCM.Y(), threeVecPiCM.X());

            double mass = dStar->DeltaInvMass();

            std::vector<double> var4nSparse = {mass, ptCand, yCand, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, cosThetaStarRandom, centrality};
            std::vector<double> var4nSparseThetaPhiStar = {mass, ptCand, thetaStarBeam, phiStarBeam};

            if(!fReadMC) {
                fnSparseReco[0]->Fill(var4nSparse.data());
                fnSparseRecoThetaPhiStar[0]->Fill(var4nSparseThetaPhiStar.data());
            }
            else
            {
                if(labD > 0) {
                    if(orig == 4) {
                        fnSparseReco[1]->Fill(var4nSparse.data());
                        fnSparseRecoThetaPhiStar[1]->Fill(var4nSparseThetaPhiStar.data());
                    }
                    else if(orig == 5) {
                        fnSparseReco[2]->Fill(var4nSparse.data());
                        fnSparseRecoThetaPhiStar[2]->Fill(var4nSparseThetaPhiStar.data());
                    }
                }
                else {
                    fnSparseReco[3]->Fill(var4nSparse.data());
                    fnSparseRecoThetaPhiStar[3]->Fill(var4nSparseThetaPhiStar.data());
                }
            }

            if (unsetVtx)
                dZeroDau->UnsetOwnPrimaryVtx();
            if (recVtx)
                fRDCuts->CleanOwnPrimaryVtx(dZeroDau, fAOD, origOwnVtx);
        }
        else {
            if(isSelected == 1 || isSelected == 3) {
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

                ROOT::Math::XYZVector threeVecPiCM = fourVecPiCM.Vect();

                double cosThetaStarProd = TMath::Abs(normalVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarHelicity = TMath::Abs(helicityVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarBeam = TMath::Abs(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarRandom = TMath::Abs(randomVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double thetaStarBeam = TMath::ACos(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double phiStarBeam = TMath::ATan2(threeVecPiCM.Y(), threeVecPiCM.X());

                double mass = dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson)->InvMassD0();

                std::vector<double> var4nSparse = {mass, ptCand, yCand, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, cosThetaStarRandom, centrality};
                std::vector<double> var4nSparseThetaPhiStar = {mass, ptCand, thetaStarBeam, phiStarBeam};

                if(!fReadMC) {
                    fnSparseReco[0]->Fill(var4nSparse.data());
                    fnSparseRecoThetaPhiStar[0]->Fill(var4nSparseThetaPhiStar.data());
                }
                else
                {
                    if(labD >= 0) {
                        //check if reflected signal
                        int labDauFirst = dauPi->GetLabel();
                        AliAODMCParticle* dauFirst = nullptr;
                        if(labDauFirst >= 0)
                            dauFirst = dynamic_cast<AliAODMCParticle *>(arrayMC->UncheckedAt(labDauFirst));
                        if(dauFirst && TMath::Abs(dauFirst->GetPdgCode()) == 211) {
                            if(orig == 4) {
                                fnSparseReco[1]->Fill(var4nSparse.data());
                                fnSparseRecoThetaPhiStar[1]->Fill(var4nSparseThetaPhiStar.data());
                            }
                            else if(orig == 5) {
                                fnSparseReco[2]->Fill(var4nSparse.data());
                                fnSparseRecoThetaPhiStar[2]->Fill(var4nSparseThetaPhiStar.data());
                            }
                        }
                    }
                    else {
                        fnSparseReco[3]->Fill(var4nSparse.data());
                        fnSparseRecoThetaPhiStar[3]->Fill(var4nSparseThetaPhiStar.data());
                    }
                }
            }
            if(isSelected >= 2) {
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

                ROOT::Math::XYZVector threeVecPiCM = fourVecPiCM.Vect();

                double cosThetaStarProd = TMath::Abs(normalVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarHelicity = TMath::Abs(helicityVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarBeam = TMath::Abs(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double cosThetaStarRandom = TMath::Abs(randomVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double thetaStarBeam = TMath::ACos(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                double phiStarBeam = TMath::ATan2(threeVecPiCM.Y(), threeVecPiCM.X());

                double mass = dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson)->InvMassD0bar();

                std::vector<double> var4nSparse = {mass, ptCand, yCand, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, cosThetaStarRandom, centrality};
                std::vector<double> var4nSparseThetaPhiStar = {mass, ptCand, thetaStarBeam, phiStarBeam};

                if(!fReadMC) {
                    fnSparseReco[0]->Fill(var4nSparse.data());
                    fnSparseRecoThetaPhiStar[0]->Fill(var4nSparseThetaPhiStar.data());
                }
                else
                {
                    if(labD >= 0) {
                        //check if reflected signal
                        int labDauFirst = dauPi->GetLabel();
                        AliAODMCParticle* dauFirst = nullptr;
                        if(labDauFirst >= 0)
                            dauFirst = dynamic_cast<AliAODMCParticle *>(arrayMC->UncheckedAt(labDauFirst));
                        if(dauFirst && TMath::Abs(dauFirst->GetPdgCode()) == 211) {
                            if(orig == 4) {
                                fnSparseReco[1]->Fill(var4nSparse.data());
                                fnSparseRecoThetaPhiStar[1]->Fill(var4nSparseThetaPhiStar.data());
                            }
                            else if(orig == 5) {
                                fnSparseReco[2]->Fill(var4nSparse.data());
                                fnSparseRecoThetaPhiStar[2]->Fill(var4nSparseThetaPhiStar.data());
                            }
                        }
                    }
                    else {
                        fnSparseReco[3]->Fill(var4nSparse.data());
                        fnSparseRecoThetaPhiStar[3]->Fill(var4nSparseThetaPhiStar.data());
                    }
                }
            }
            if (unsetVtx)
                dMeson->UnsetOwnPrimaryVtx();
            if (recVtx)
                fRDCuts->CleanOwnPrimaryVtx(dMeson, fAOD, origOwnVtx);
        }
    }

    PostData(1, fOutput);
}

//________________________________________________________________________
int AliAnalysisTaskSEDstarPolarization::IsCandidateSelected(AliAODRecoDecayHF *&d, AliAODRecoDecayHF2Prong *&dZeroDau, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx, std::vector<double> scoresFromMLSelector, std::vector<double> scoresFromMLSelectorSecond)
{

    if (!d || (!dZeroDau && fDecChannel == kDstartoD0pi) || !vHF)
        return 0;
    fHistNEvents->Fill(11);

    AliAODRecoCascadeHF* dStar = nullptr;
    int nDau = 0;
    if(fDecChannel == kDstartoD0pi) {
        dStar = dynamic_cast<AliAODRecoCascadeHF*>(d);
        nDau = 3;
    }
    else {
        nDau = 2;
    }

    // Preselection to speed up task
    TObjArray arrDauTracks(nDau);

    for (int iDau = 0; iDau < nDau; iDau++) {
        AliAODTrack *track = nullptr;
        if(fDecChannel == kDstartoD0pi) {
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

    int isSelected = (fDecChannel == kDstartoD0pi) ? fRDCuts->IsSelected(dStar, AliRDHFCuts::kAll, fAOD) : fRDCuts->IsSelected(d, AliRDHFCuts::kAll, fAOD);
    if (!isSelected) {
        if (unsetVtx) {
            if (fDecChannel == kDstartoD0pi)
                dZeroDau->UnsetOwnPrimaryVtx();
            else
                d->UnsetOwnPrimaryVtx();
        }
        return 0;
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

    if(!fApplyML) {
        return isSelected;
    }
    else {
        //variables for ML application
        std::vector<double> modelPred = {};
        int isMLsel = isSelected;
        double ptCand = d->Pt();

        AliAODPidHF *pidHF = fRDCuts->GetPidHF();

        if(fDependOnMLSelector) {
            std::vector<float>::iterator low = std::lower_bound(fPtLimsML.begin(), fPtLimsML.end(), ptCand);
            int bin = low - fPtLimsML.begin() - 1;
            if(bin < 0)
                bin = 0;
            else if(bin > fPtLimsML.size()-2)
                bin = fPtLimsML.size()-2;

            if(isSelected == 1 || isSelected == 3) {
                for(size_t iScore = 0; iScore < scoresFromMLSelector.size(); iScore++) {
                    if((fMLOptScoreCuts[bin][iScore] == "upper" && scoresFromMLSelector[iScore] > fMLScoreCuts[bin][iScore]) ||
                       (fMLOptScoreCuts[bin][iScore] == "lower" && scoresFromMLSelector[iScore] < fMLScoreCuts[bin][iScore]))
                    {
                        isMLsel -= 1;
                        break;
                    }
                }
            }
            if(isSelected >= 2) {
                for(size_t iScore = 0; iScore < scoresFromMLSelectorSecond.size(); iScore++) {
                    if((fMLOptScoreCuts[bin][iScore] == "upper" && scoresFromMLSelectorSecond[iScore] > fMLScoreCuts[bin][iScore]) ||
                    (fMLOptScoreCuts[bin][iScore] == "lower" && scoresFromMLSelectorSecond[iScore] < fMLScoreCuts[bin][iScore]))
                    {
                        isMLsel -= 2;
                        break;
                    }
                }
            }
        }
        else {
            if(isSelected == 1 || isSelected == 3) {
                if(!fMLResponse->IsSelectedMultiClass(modelPred, d, fAOD->GetMagneticField(), pidHF, 0))
                    isMLsel -= 1;
            }
            if(isSelected >= 2) {
                if(!fMLResponse->IsSelectedMultiClass(modelPred, d, fAOD->GetMagneticField(), pidHF, 1))
                    isMLsel -= 2;
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
                if(fDecChannel == kDstartoD0pi) {
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
                        if(labDauFirst < 0)
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

                        double phiRandom = gRandom->Uniform(0., 2*TMath::Pi());
                        double thetaRandom = gRandom->Uniform(0., TMath::Pi());
                        ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(TMath::Sin(thetaRandom) * TMath::Cos(phiRandom), TMath::Sin(thetaRandom) * TMath::Sin(phiRandom), TMath::Cos(thetaRandom));

                        double cosThetaStarProd = TMath::Abs(normalVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double cosThetaStarHelicity = TMath::Abs(helicityVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double cosThetaStarBeam = TMath::Abs(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double cosThetaStarRandom = TMath::Abs(randomVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double thetaStarBeam = TMath::ACos(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double phiStarBeam = TMath::ATan2(threeVecPiCM.Y(), threeVecPiCM.X());

                        if (orig == 4 && !isParticleFromOutOfBunchPileUpEvent)
                        {
                            double var4nSparseAcc[knVarForSparseAcc] = {pt, rapid, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, cosThetaStarRandom, centrality};
                            double var4nSparseAccThetaPhiStar[3] = {pt, thetaStarBeam, phiStarBeam};
                            fnSparseMC[0]->Fill(var4nSparseAcc);
                            fnSparseMCThetaPhiStar[0]->Fill(var4nSparseAccThetaPhiStar);
                        }
                        else if (orig == 5 && !isParticleFromOutOfBunchPileUpEvent)
                        {
                            double var4nSparseAcc[knVarForSparseAcc] = {pt, rapid, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, cosThetaStarRandom, centrality};
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
        if(mother && TMath::Abs(mother->GetPdgCode()) == 413)
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

    int nBinsAcc[knVarForSparseAcc] = {nPtBins, 100, 5, 5, 5, 5, 100};
    double xminAcc[knVarForSparseAcc] = {0., -1., 0., 0., 0., 0., 0.};
    double xmaxAcc[knVarForSparseAcc] = {ptLims[nPtBinsCutObj], 1., 1., 1., 1., 1., 100.};

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
        fnSparseMC[iHist]->GetAxis(5)->SetTitle("|cos(#theta*)| (random)");
        fnSparseMC[iHist]->GetAxis(6)->SetTitle("centrality");
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
    if(fDecChannel == kD0toKpi) {
        massMin = 1.65;
        massMax = 2.15;
        massTitle = "#it{M}(K#pi)";
    }

    int nCosThetaBins = 5;

    int nBinsReco[knVarForSparseReco] = {nMassBins, nPtBins, 100, nCosThetaBins, nCosThetaBins, nCosThetaBins, nCosThetaBins, 100};
    double xminReco[knVarForSparseReco] = {massMin, 0., -1., 0., 0., 0., 0., 0.};
    double xmaxReco[knVarForSparseReco] = {massMax, ptLims[nPtBinsCutObj], 1., 1., 1., 1., 1., 100.};

    int nBinsThetaPhiReco[4] = {nMassBins, nPtBins, 100, 100};
    double xminThetaPhiReco[4] = {massMin, 0., 0., 0.};
    double xmaxThetaPhiReco[4] = {massMax, ptLims[nPtBinsCutObj], TMath::Pi(), TMath::Pi()};

    TString label[4] = {"all", "fromC", "fromB", "bkg"};
    for (int iHist = 0; iHist < 4; iHist++)
    {
        TString titleSparse = Form("Reco nSparse - %s", label[iHist].Data());
        fnSparseReco[iHist] = new THnSparseF(Form("fnSparseReco_%s", label[iHist].Data()), titleSparse.Data(), knVarForSparseReco, nBinsReco, xminReco, xmaxReco);
        fnSparseReco[iHist]->GetAxis(0)->SetTitle(Form("%s (GeV/#it{c}^{2})", massTitle.Data()));
        fnSparseReco[iHist]->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
        fnSparseReco[iHist]->GetAxis(2)->SetTitle("#it{y}");
        fnSparseReco[iHist]->GetAxis(3)->SetTitle("|cos(#theta*)| (beam)");
        fnSparseReco[iHist]->GetAxis(4)->SetTitle("|cos(#theta*)| (production)");
        fnSparseReco[iHist]->GetAxis(5)->SetTitle("|cos(#theta*)| (helicity)");
        fnSparseReco[iHist]->GetAxis(6)->SetTitle("|cos(#theta*)| (random)");
        fnSparseReco[iHist]->GetAxis(7)->SetTitle("centrality %");
        fOutput->Add(fnSparseReco[iHist]);

        fnSparseRecoThetaPhiStar[iHist] = new THnSparseF(Form("fnSparseRecoThetaPhiStar_%s", label[iHist].Data()), titleSparse.Data(), 4, nBinsThetaPhiReco, xminThetaPhiReco, xmaxThetaPhiReco);
        fnSparseRecoThetaPhiStar[iHist]->GetAxis(0)->SetTitle(Form("%s (GeV/#it{c}^{2})", massTitle.Data()));
        fnSparseRecoThetaPhiStar[iHist]->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
        fnSparseRecoThetaPhiStar[iHist]->GetAxis(2)->SetTitle("#theta* (beam)");
        fnSparseRecoThetaPhiStar[iHist]->GetAxis(3)->SetTitle("#varphi* (beam)");
        fOutput->Add(fnSparseRecoThetaPhiStar[iHist]);
    }
}
