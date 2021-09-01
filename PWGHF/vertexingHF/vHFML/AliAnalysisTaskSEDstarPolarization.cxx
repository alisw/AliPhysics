/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSEDstarPolarization
// \brief Analysis task to perform D*+ polarization analysis
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
// S. Kundu, sourav.kundu@cern.ch
/////////////////////////////////////////////////////////////

#include <TRandom3.h>

#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisUtils.h"
#include "AliAODHandler.h"
#include "AliAODExtension.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"

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

    AliRDHFCutsDStartoKpipi *copycut = new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi *>(fRDCuts)));
    PostData(2, copycut);

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
    if(fApplyML)
    {
        fMLResponse = new AliHFMLResponseDstartoD0pi("DstartoD0piMLResponse", "DstartoD0piMLResponse", fConfigPath.data());
        fMLResponse->MLResponseInit();
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
        }
    }
    else if (fAOD)
    {
        arrayCand = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Dstar"));
    }

    if (!fAOD || !arrayCand)
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

    // vHF object is needed to call the method that refills the missing info of the candidates
    // if they have been deleted in dAOD reconstruction phase
    // in order to reduce the size of the file
    AliAnalysisVertexingHF vHF = AliAnalysisVertexingHF();

    for (int iCand = 0; iCand < arrayCand->GetEntriesFast(); iCand++)
    {
        AliAODRecoCascadeHF *dStar = dynamic_cast<AliAODRecoCascadeHF *>(arrayCand->UncheckedAt(iCand));

        bool unsetVtx = false;
        bool recVtx = false;
        AliAODVertex *origOwnVtx = nullptr;

        int isSelected = IsCandidateSelected(dStar, &vHF, unsetVtx, recVtx, origOwnVtx);
        if (!isSelected)
        {
            if (unsetVtx)
                dStar->UnsetOwnPrimaryVtx();
            if (recVtx)
                fRDCuts->CleanOwnPrimaryVtx(dStar, fAOD, origOwnVtx);
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
            labD = dStar->MatchToMC(413, 421, pdgDstarDau, pdgD0Dau, arrayMC, false);
            partD = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labD));

            if (partD)
                orig = AliVertexingHFUtils::CheckOrigin(arrayMC, partD, true);
        }

        // actual analysis
        double mass = dStar->DeltaInvMass();
        double ptCand = dStar->Pt();
        double pCand = dStar->P();

        AliAODTrack* dauPi = dynamic_cast<AliAODTrack *>(dStar->GetBachelor());
        AliAODRecoDecayHF2Prong* dauD0 = dynamic_cast<AliAODRecoDecayHF2Prong *>(dStar->Get2Prong());
        fourVecPi = ROOT::Math::PxPyPzMVector(dauPi->Px(), dauPi->Py(), dauPi->Pz(), TDatabasePDG::Instance()->GetParticle(211)->Mass());
        fourVecD0 = ROOT::Math::PxPyPzMVector(dauD0->Px(), dauD0->Py(), dauD0->Pz(), TDatabasePDG::Instance()->GetParticle(421)->Mass());
        fourVecDstar = fourVecPi + fourVecD0;

        ROOT::Math::Boost boostv12{fourVecDstar.BoostToCM()};
        fourVecPiCM = boostv12(fourVecPi);
        fourVecD0CM = boostv12(fourVecD0);

        ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(dStar->Py() / ptCand, -dStar->Px() / ptCand, 0.);
        ROOT::Math::XYZVector helicityVec = ROOT::Math::XYZVector(dStar->Px() / pCand, dStar->Py() / pCand, dStar->Pz() / pCand);
        ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0., 0., 1.);

        ROOT::Math::XYZVector threeVecPiCM = fourVecPiCM.Vect();

        double cosThetaStarProd = TMath::Abs(normalVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
        double cosThetaStarHelicity = TMath::Abs(helicityVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
        double cosThetaStarBeam = TMath::Abs(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));

        std::vector<double> var4nSparse = {mass, ptCand, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, centrality};

        if(!fReadMC)
            fnSparseReco[0]->Fill(var4nSparse.data());
        else
        {
            if(labD > 0) {
                if(orig == 4)
                    fnSparseReco[1]->Fill(var4nSparse.data());
                else if(orig == 5)
                    fnSparseReco[2]->Fill(var4nSparse.data());
            }
            else
                fnSparseReco[3]->Fill(var4nSparse.data());
        }

        if (unsetVtx)
            dStar->UnsetOwnPrimaryVtx();
        if (recVtx)
            fRDCuts->CleanOwnPrimaryVtx(dStar, fAOD, origOwnVtx);
    }

    PostData(1, fOutput);
}

//________________________________________________________________________
int AliAnalysisTaskSEDstarPolarization::IsCandidateSelected(AliAODRecoCascadeHF *&dStar, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx)
{

    if (!dStar || !vHF)
        return 0;
    fHistNEvents->Fill(11);

    // Preselection to speed up task
    TObjArray arrDauTracks(3);
    int nDau = 3;

    for (int iDau = 0; iDau < nDau; iDau++)
    {
        AliAODTrack *track = vHF->GetProng(fAOD, dStar, iDau);
        arrDauTracks.AddAt(track, iDau);
    }

    if (!fRDCuts->PreSelect(arrDauTracks))
    {
        fHistNEvents->Fill(15);
        return 0;
    }

    if (!vHF->FillRecoCasc(fAOD, dStar, false))
    {
        fHistNEvents->Fill(14);
        return 0;
    }

    fHistNEvents->Fill(12);

    unsetVtx = false;
    if (!dStar->GetOwnPrimaryVtx())
    {
        dStar->SetOwnPrimaryVtx(dynamic_cast<AliAODVertex *>(fAOD->GetPrimaryVertex()));
        unsetVtx = true;
        // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
        // Pay attention if you use continue inside this loop!!!
    }

    double ptD = dStar->Pt();
    double yD = dStar->Y(413);

    int ptBin = fRDCuts->PtBin(ptD);
    if (ptBin < 0)
    {
        if (unsetVtx)
            dStar->UnsetOwnPrimaryVtx();
        return 0;
    }

    bool isFidAcc = fRDCuts->IsInFiducialAcceptance(ptD, yD);
    if (!isFidAcc)
    {
        if (unsetVtx)
            dStar->UnsetOwnPrimaryVtx();
        return 0;
    }

    int isSelected = fRDCuts->IsSelected(dStar, AliRDHFCuts::kAll, fAOD);
    if (!isSelected)
    {
        if (unsetVtx)
            dStar->UnsetOwnPrimaryVtx();
        return 0;
    }

    recVtx = false;
    origOwnVtx = nullptr;

    if (fRDCuts->GetIsPrimaryWithoutDaughters())
    {
        if (dStar->GetOwnPrimaryVtx())
            origOwnVtx = new AliAODVertex(*dStar->GetOwnPrimaryVtx());
        if (fRDCuts->RecalcOwnPrimaryVtx(dStar, fAOD))
            recVtx = true;
        else
            fRDCuts->CleanOwnPrimaryVtx(dStar, fAOD, origOwnVtx);
    }


    if(fApplyML)
    {
        //variables for ML application
        std::vector<double> modelPred = {};
        bool isMLsel = false;

        AliAODPidHF *pidHF = fRDCuts->GetPidHF();
        isMLsel = fMLResponse->IsSelectedMultiClass(modelPred, dStar, fAOD->GetMagneticField(), pidHF, 0);
        if(isMLsel)
            isSelected = 1;
    }

    return isSelected;
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
            AliAODMCParticle *mcPart = dynamic_cast<AliAODMCParticle *>(arrayMC->At(iPart));

            if (TMath::Abs(mcPart->GetPdgCode()) == 413)
            {
                int orig = AliVertexingHFUtils::CheckOrigin(arrayMC, mcPart, true); //Prompt = 4, FeedDown = 5
                bool isParticleFromOutOfBunchPileUpEvent = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iPart, mcHeader, arrayMC);

                int deca = 0;
                bool isGoodDecay = false;
                int labDau[3] = {-1, -1, -1};
                bool isFidAcc = false;
                bool isDaugInAcc = false;
                int nDau = 3;
                deca = AliVertexingHFUtils::CheckDstarDecay(arrayMC, mcPart, labDau);

                if (labDau[0] == -1)
                    continue; //protection against unfilled array of labels
                if (deca > 0)
                    isGoodDecay = true;

                if (isGoodDecay)
                {
                    double pt = mcPart->Pt();
                    double p = mcPart->P();
                    double rapid = mcPart->Y();
                    isFidAcc = fRDCuts->IsInFiducialAcceptance(pt, rapid);
                    isDaugInAcc = CheckDaugAcc(arrayMC, nDau, labDau);

                    if ((fFillAcceptanceLevel && isFidAcc && isDaugInAcc) || (!fFillAcceptanceLevel && TMath::Abs(rapid) < 0.5))
                    {
                        int labDauFirst = mcPart->GetDaughterFirst();
                        AliAODMCParticle* dauFirst = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labDauFirst));
                        fourVecDstar = ROOT::Math::PxPyPzMVector(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcPart->M());
                        fourVecPi = ROOT::Math::PxPyPzMVector(dauFirst->Px(), dauFirst->Py(), dauFirst->Pz(), dauFirst->M());

                        ROOT::Math::Boost boostv12{fourVecDstar.BoostToCM()};
                        fourVecPiCM = boostv12(fourVecPi);

                        ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(mcPart->Py() / pt, -mcPart->Px() / pt, 0.);
                        ROOT::Math::XYZVector helicityVec = ROOT::Math::XYZVector(mcPart->Px() / p, mcPart->Py() / p, mcPart->Pz() / p);
                        ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0., 0., 1.);

                        ROOT::Math::XYZVector threeVecPiCM = fourVecPiCM.Vect();

                        double cosThetaStarProd = TMath::Abs(normalVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double cosThetaStarHelicity = TMath::Abs(helicityVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));
                        double cosThetaStarBeam = TMath::Abs(beamVec.Dot(threeVecPiCM) / TMath::Sqrt(threeVecPiCM.Mag2()));

                        if (orig == 4 && !isParticleFromOutOfBunchPileUpEvent)
                        {
                            double var4nSparseAcc[knVarForSparseAcc] = {pt, rapid, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, centrality};
                            fnSparseMC[0]->Fill(var4nSparseAcc);
                        }
                        else if (orig == 5 && !isParticleFromOutOfBunchPileUpEvent)
                        {
                            double var4nSparseAcc[knVarForSparseAcc] = {pt, rapid, cosThetaStarBeam, cosThetaStarProd, cosThetaStarHelicity, centrality};
                            fnSparseMC[1]->Fill(var4nSparseAcc);
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
        AliAODMCParticle *mcPartDaughter = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labDau[iProng]));
        if (!mcPartDaughter)
            return false;

        AliAODMCParticle *mother = dynamic_cast<AliAODMCParticle *>(arrayMC->At(mcPartDaughter->GetMother()));
        if(TMath::Abs(mother->GetPdgCode()) == 413)
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

    int nBinsAcc[knVarForSparseAcc] = {nPtBins, 20, 5, 5, 5, 100};
    double xminAcc[knVarForSparseAcc] = {0., -1., 0., 0., 0., 0.};
    double xmaxAcc[knVarForSparseAcc] = {ptLims[nPtBinsCutObj], 1., 1., 1., 1., 100.};

    TString label[2] = {"fromC", "fromB"};
    for (int iHist = 0; iHist < 2; iHist++)
    {
        TString titleSparse = Form("MC nSparse (%s)- %s", fFillAcceptanceLevel ? "Acc.Step" : "Gen.Acc.Step", label[iHist].Data());
        fnSparseMC[iHist] = new THnSparseF(Form("fnSparseAcc_%s", label[iHist].Data()), titleSparse.Data(), knVarForSparseAcc, nBinsAcc, xminAcc, xmaxAcc);
        fnSparseMC[iHist]->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/c)");
        fnSparseMC[iHist]->GetAxis(1)->SetTitle("#it{y}");
        fnSparseMC[iHist]->GetAxis(2)->SetTitle("|cos(#theta*)| (beam)");
        fnSparseMC[iHist]->GetAxis(3)->SetTitle("|cos(#theta*)| (production)");
        fnSparseMC[iHist]->GetAxis(4)->SetTitle("|cos(#theta*)| (helicity)");
        fnSparseMC[iHist]->GetAxis(5)->SetTitle("centrality");
        fOutput->Add(fnSparseMC[iHist]);
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

    int nCosThetaBins = 5;

    int nBinsReco[knVarForSparseReco] = {nMassBins, nPtBins, nCosThetaBins, nCosThetaBins, nCosThetaBins, 100};
    double xminReco[knVarForSparseReco] = {massMin, 0., 0., 0., 0., 0.};
    double xmaxReco[knVarForSparseReco] = {massMax, ptLims[nPtBinsCutObj], 1., 1., 1., 100.};

    TString label[4] = {"all", "fromC", "fromB", "bkg"};
    for (int iHist = 0; iHist < 4; iHist++)
    {
        TString titleSparse = Form("Reco nSparse - %s", label[iHist].Data());
        fnSparseReco[iHist] = new THnSparseF(Form("fnSparseReco_%s", label[iHist].Data()), titleSparse.Data(), knVarForSparseReco, nBinsReco, xminReco, xmaxReco);
        fnSparseReco[iHist]->GetAxis(0)->SetTitle("#it{M}(K#pi#pi) #minus #it{M}(K#pi) (GeV/#it{c})");
        fnSparseReco[iHist]->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/c)");
        fnSparseReco[iHist]->GetAxis(2)->SetTitle("|cos(#theta*)| (beam)");
        fnSparseReco[iHist]->GetAxis(3)->SetTitle("|cos(#theta*)| (production)");
        fnSparseReco[iHist]->GetAxis(4)->SetTitle("|cos(#theta*)| (helicity)");
        fnSparseReco[iHist]->GetAxis(5)->SetTitle("centrality %");

        fOutput->Add(fnSparseReco[iHist]);
    }
}
