#include <Riostream.h>

#include <stdio.h>
#include <iostream>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"

#include "AliAnalysisManager.h"

#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliMultiplicity.h"

#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODMCHeader.h"
#include "AliAODInputHandler.h"
#include "AliPIDResponse.h"

#include "AliAnalysisTaskFlatenicityLambdaK0s.h"

class AliAnalysisTaskFlatenicityLambdaK0s;

using namespace std;

Int_t gTrack;
Float_t invmassK0s;

ClassImp(AliAnalysisTaskFlatenicityLambdaK0s)

    AliAnalysisTaskFlatenicityLambdaK0s::AliAnalysisTaskFlatenicityLambdaK0s() : AliAnalysisTaskSE(),
                                                     fAOD(0), fOutputList(0), hinvmassK0s(0), fESD(0), fPIDResponse(0), hinvmassLambda(0), hinvmassAntiLambda(0), hflat(0), fESDpid(0x0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskFlatenicityLambdaK0s::AliAnalysisTaskFlatenicityLambdaK0s(const char *name) : AliAnalysisTaskSE(name),
                                                                 fAOD(0), fOutputList(0), hinvmassK0s(0), fESD(0), fPIDResponse(0), hinvmassLambda(0), hinvmassAntiLambda(0), hflat(0), fESDpid(0x0)
{
    // constructor
    DefineInput(0, TChain::Class()); // define the input of the analysis: in this case we take a 'chain' of events
                                     // this chain is created by the analysis manager, so no need to worry about it,
                                     // it does its work automatically
    DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this case it's a list of histograms
                                     // you can add more output objects by calling DefineOutput(2, classname::Class())
                                     // if you add more output objects, make sure to call PostData for all of them, and to
                                     // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskFlatenicityLambdaK0s::~AliAnalysisTaskFlatenicityLambdaK0s()
{
    // destructor
    if (fOutputList)
    {
        delete fOutputList; // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityLambdaK0s::UserCreateOutputObjects()
{

    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    // example of a histogram
    hinvmassK0s = new TH1F("hinvmassK0s", "hinvmassK0s", 100, 0.41, 0.58);
    hinvmassLambda = new TH1F("hinvmassLambda", "hinvmassLambda", 100, 1.08, 1.15);
    hinvmassAntiLambda = new TH1F("hinvmassAntiLambda", "hinvmassAntiLambda", 100, 1.08, 1.15);
    hflat = new TH1F("hflat", "hflat", 101, -0.005, 1.005);


    fOutputList->Add(hinvmassK0s);
    fOutputList->Add(hinvmassLambda);
    fOutputList->Add(hinvmassAntiLambda);
    fOutputList->Add(hflat);

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityLambdaK0s::UserExec(Option_t *)
{
    fESD = (AliESDEvent *)InputEvent();

    if (!fESD)
    {
        Printf("ERROR: fESD not available");
        return;
    }

    Double_t lPLambda = 0;
    Double_t lPAntiLambda = 0;
    Double_t lPK0s = 0;
    Double_t lMagneticField = 999;

    // Multiplcity:
    Int_t nv0sTot = 0, nv0s = 0;
    //  Int_t nv0sMI =0;
    // Variables:
    Double_t lV0Position[3];

    Double_t lDcaPosToPrimVertex = 0;
    Double_t lDcaNegToPrimVertex = 0;
    Double_t lDcaV0Daughters = 0;
    Double_t lV0cosPointAngle = 0;
    Double_t lChi2V0 = 0;
    Double_t lV0DecayLength = 0;
    Double_t lV0Radius = 0;
    Double_t lDcaV0ToPrimVertex = 0;
    Double_t lcTauLambda = 0;
    Double_t lcTauAntiLambda = 0;
    Double_t lcTauK0s = 0;
    Int_t lOnFlyStatus = 0;
    // Float_t   tdcaPosToPrimVertexXYZ[2], tdcaNegToPrimVertexXYZ[2]; // ..[0] = Impact parameter in XY plane and ..[1] = Impact parameter in Z
    // Double_t  tdcaDaughterToPrimVertex[2];                          // ..[0] = Pos and ..[1] = Neg

    Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
    Double_t lPtK0s = 0, lPtLambda = 0, lPtAntiLambda = 0;
    Double_t lRapK0s = 0, lRapLambda = 0, lRapAntiLambda = 0;
    //  Double_t lEtaK0s     = 0, lEtaLambda     = 0, lEtaAntiLambda     = 0;
    Double_t lAlphaV0 = 0, lPtArmV0 = 0;

    Double_t lPzK0s = 0, lPzLambda = 0, lPzAntiLambda = 0;

    Double_t lV0Eta = 999;

    // to study Associated V0s:
    Int_t lIndexTrackPos = 0, lIndexTrackNeg = 0;
    UInt_t lLabelTrackPos = 0, lLabelTrackNeg = 0;
    Int_t lCheckPIdK0Short = 0, lCheckMcK0Short = 0;
    Int_t lCheckPIdLambda = 0, lCheckMcLambda = 0;
    Int_t lCheckPIdAntiLambda = 0, lCheckMcAntiLambda = 0;
    Int_t lCheckSecondaryK0s = 0, lCheckSecondaryLambda = 0, lCheckSecondaryAntiLambda = 0;
    Int_t lCheckGamma = 0;
    Double_t mcPosMotherX = 0, mcPosMotherY = 0, mcPosMotherZ = 0;
    Double_t mcPosMotherR = 0;
    Double_t mcMotherPt = 0, mcMotherRap = 0;

    Int_t lIndexPosMother = 0;
    Int_t lIndexNegMother = 0;
    Int_t lIndexMotherOfMother = 0;
    Int_t lPDGCodePosDaughter = 0;
    Int_t lPDGCodeNegDaughter = 0;
    Int_t lPdgcodeMother = 0;
    Int_t lPdgcodeMotherOfMother = 0;

    // Reconstructed position
    // Double_t rcPosXK0s        = 0,  rcPosYK0s        = 0, rcPosZK0s        = 0;
    Double_t rcPosRK0s = 0;
    // Double_t rcPosXLambda     = 0,  rcPosYLambda     = 0, rcPosZLambda     = 0;
    Double_t rcPosRLambda = 0;
    //  Double_t rcPosXAntiLambda = 0,  rcPosYAntiLambda = 0, rcPosZAntiLambda = 0;
    Double_t rcPosRAntiLambda = 0;

    // Pt resolution
    // Double_t deltaPtK0s  = 0, deltaPtLambda  = 0, deltaPtAntiLambda  = 0;
    AliESDtrack *myTrackPos = NULL;
    AliESDtrack *myTrackNeg = NULL;

    // Daughters' momentum:
    Double_t lMomPos[3] = {999, 999, 999};
    Double_t lMomNeg[3] = {999, 999, 999};
    Double_t lPtPos = 999, lPtNeg = 999;
    Double_t lPPos = 999, lPNeg = 999;

    // Inner Wall parameters:
    Double_t lMomInnerWallPos = 999, lMomInnerWallNeg = 999;

    // AliKF Chi2 and Armenteros variables
    //  Double_t lChi2KFK0s  = 0, lChi2KFLambda = 0,  lChi2KFAntiLambda = 0;
    //  Double_t lAlphaV0K0s = 0, lAlphaV0Lambda = 0,  lAlphaV0AntiLambda = 0;
    // Double_t lPtArmV0K0s = 0, lPtArmV0Lambda = 0,  lPtArmV0AntiLambda = 0;
    //  Double_t lQlPos   = 0, lQlNeg   = 0;

    // PID
    Float_t nSigmaPosPion = 0;
    Float_t nSigmaNegPion = 0;

    Float_t nSigmaPosProton = 0;
    Float_t nSigmaNegProton = 0;

    Int_t lCheckPIDK0sPosDaughter = 0, lCheckPIDK0sNegDaughter = 0;
    Int_t lCheckPIDLambdaPosDaughter = 0, lCheckPIDLambdaNegDaughter = 0;
    Int_t lCheckPIDAntiLambdaPosDaughter = 0, lCheckPIDAntiLambdaNegDaughter = 0;

    //****************************************************
    // Primary Vertex cuts &
    // Magnetic field and Quality tracks cuts
    //****************************************************

    Double_t lPrimaryVtxPosition[3];
    Double_t lPrimaryVtxCov[6];
    Double_t lPrimaryVtxChi2 = 999;
    Double_t lResPrimaryVtxX = 999;
    Double_t lResPrimaryVtxY = 999;
    Double_t lResPrimaryVtxZ = 999;

    AliAODVertex *myPrimaryVertex = NULL;

    const AliESDVertex *myBestPrimaryVertex = ((AliESDEvent *)fESD)->GetPrimaryVertex();
    myBestPrimaryVertex = ((AliESDEvent *)fESD)->GetPrimaryVertex();
    if (!myBestPrimaryVertex)
        return;
    if (!myBestPrimaryVertex->GetStatus())
        return;
    myBestPrimaryVertex->GetXYZ(lPrimaryVtxPosition);
    myBestPrimaryVertex->GetCovMatrix(lPrimaryVtxCov);

    if ((TMath::Abs(lPrimaryVtxPosition[2])) > 10.0)
        return; //// cut on z of prim. vertex!!!!!

    lPrimaryVtxChi2 = myBestPrimaryVertex->GetChi2toNDF();
    lResPrimaryVtxX = myBestPrimaryVertex->GetXRes();
    lResPrimaryVtxY = myBestPrimaryVertex->GetYRes();
    lResPrimaryVtxZ = myBestPrimaryVertex->GetZRes();

    // const AliESDVertex *mySPDPrimaryVertex = ((AliESDEvent *)fESD)->GetPrimaryVertexSPD();
    // if (!mySPDPrimaryVertex)
    //     return;
    // const AliESDVertex *myPrimaryVertexTracking = ((AliESDEvent *)fESD)->GetPrimaryVertexTracks();
    // if (!myPrimaryVertexTracking)
    //     return;

    // if (!mySPDPrimaryVertex->GetStatus() && !myPrimaryVertexTracking->GetStatus())
    //     return;

    myPrimaryVertex = new AliAODVertex(lPrimaryVtxPosition, lPrimaryVtxCov, lPrimaryVtxChi2, NULL, -1, AliAODVertex::kPrimary);
    if (!myPrimaryVertex)
        return;

    lMagneticField = ((AliESDEvent *)fESD)->GetMagneticField();
    // AliKFVertex primaryVtxKF(*myPrimaryVertex);
    // AliKFParticle::SetField(lMagneticField);

    //=========================================================================================================
    fESDpid = new AliESDpid;
    fESDpid->GetTPCResponse().SetBetheBlochParameters(1.41543 / 50.0, 2.63394E1, 5.0411E-11, 2.12543, 4.88663);
    //=========================================================================================================

    if (fESD->IsPileupFromSPD())
        return;
    Double_t flat = GetFlatenicityV0();
    ((TH1D *)(fOutputList->FindObject("hflat")))->Fill(1.0-flat);

    nv0sTot = fESD->GetNumberOfV0s();
    Int_t lComeFromSigma = 0;
    Int_t cnt = 0;
    for (Int_t iV0 = 0; iV0 < nv0sTot; iV0++)
    {

        lIndexPosMother = 0;
        lIndexNegMother = 0;
        lIndexMotherOfMother = 0;
        lCheckPIdK0Short = 0;
        lCheckMcK0Short = 0;
        lCheckSecondaryK0s = 0;
        lCheckPIdLambda = 0;
        lCheckMcLambda = 0;
        lCheckSecondaryLambda = 0;
        lCheckPIdAntiLambda = 0;
        lCheckMcAntiLambda = 0;
        lCheckSecondaryAntiLambda = 0;
        lComeFromSigma = 0;
        lCheckGamma = 0;

        AliESDv0 *v0 = ((AliESDEvent *)fESD)->GetV0(iV0);
        if (!v0)
            continue;
        // V0's Daughters
        lIndexTrackPos = TMath::Abs(v0->GetPindex());
        lIndexTrackNeg = TMath::Abs(v0->GetNindex());
        AliESDtrack *myTrackPosTest = ((AliESDEvent *)fESD)->GetTrack(lIndexTrackPos);
        AliESDtrack *myTrackNegTest = ((AliESDEvent *)fESD)->GetTrack(lIndexTrackNeg);
        if (!myTrackPosTest || !myTrackNegTest)
        {
            Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
            continue;
        }
        // Remove like-sign
        if ((Int_t)myTrackPosTest->GetSign() == (Int_t)myTrackNegTest->GetSign())
        {
            continue;
        }
        // VO's main characteristics to check the reconstruction cuts
        lOnFlyStatus = v0->GetOnFlyStatus();
        lChi2V0 = v0->GetChi2V0();
        lDcaV0Daughters = v0->GetDcaV0Daughters();
        lDcaV0ToPrimVertex = v0->GetD(lPrimaryVtxPosition[0], lPrimaryVtxPosition[1], lPrimaryVtxPosition[2]);
        lV0cosPointAngle = v0->GetV0CosineOfPointingAngle(lPrimaryVtxPosition[0], lPrimaryVtxPosition[1], lPrimaryVtxPosition[2]);
        v0->GetXYZ(lV0Position[0], lV0Position[1], lV0Position[2]);

        lV0Radius = TMath::Sqrt(lV0Position[0] * lV0Position[0] + lV0Position[1] * lV0Position[1]);
        lV0DecayLength = TMath::Sqrt(TMath::Power(lV0Position[0] - lPrimaryVtxPosition[0], 2) +
                                     TMath::Power(lV0Position[1] - lPrimaryVtxPosition[1], 2) +
                                     TMath::Power(lV0Position[2] - lPrimaryVtxPosition[2], 2));
        if (lV0DecayLength < 0.5)
            continue;

        if (myTrackPosTest->GetSign() == 1)
        {

            myTrackPos = ((AliESDEvent *)fESD)->GetTrack(lIndexTrackPos);
            myTrackNeg = ((AliESDEvent *)fESD)->GetTrack(lIndexTrackNeg);

            // Daughters' momentum;
            v0->GetPPxPyPz(lMomPos[0], lMomPos[1], lMomPos[2]);
            v0->GetNPxPyPz(lMomNeg[0], lMomNeg[1], lMomNeg[2]);
        }

        if (myTrackPosTest->GetSign() == -1)
        {

            myTrackPos = ((AliESDEvent *)fESD)->GetTrack(lIndexTrackNeg);
            myTrackNeg = ((AliESDEvent *)fESD)->GetTrack(lIndexTrackPos);

            // Daughters' momentum;
            v0->GetPPxPyPz(lMomNeg[0], lMomNeg[1], lMomNeg[2]);
            v0->GetNPxPyPz(lMomPos[0], lMomPos[1], lMomPos[2]);
        }

        Float_t lPosTrackCrossedRowsOverFindable = myTrackPos->GetTPCClusterInfo(2, 1) / ((double)(myTrackPos->GetTPCNclsF()));
        Float_t lNegTrackCrossedRowsOverFindable = myTrackNeg->GetTPCClusterInfo(2, 1) / ((double)(myTrackNeg->GetTPCNclsF()));

        Float_t fTreeVariableLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
        if (lNegTrackCrossedRowsOverFindable < fTreeVariableLeastRatioCrossedRowsOverFindable)
            fTreeVariableLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;

        // Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
        if ((fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8))
            continue;

        lLabelTrackPos = (UInt_t)TMath::Abs(myTrackPos->GetLabel());
        lLabelTrackNeg = (UInt_t)TMath::Abs(myTrackNeg->GetLabel());

        // Daughters Pt and P:
        lPtPos = TMath::Sqrt(lMomPos[0] * lMomPos[0] + lMomPos[1] * lMomPos[1]);
        lPtNeg = TMath::Sqrt(lMomNeg[0] * lMomNeg[0] + lMomNeg[1] * lMomNeg[1]);

        lPPos = TMath::Sqrt(lMomPos[0] * lMomPos[0] + lMomPos[1] * lMomPos[1] + lMomPos[2] * lMomPos[2]);
        lPNeg = TMath::Sqrt(lMomNeg[0] * lMomNeg[0] + lMomNeg[1] * lMomNeg[1] + lMomNeg[2] * lMomNeg[2]);

        // DCA between daughter and Primary Vertex:
        if (myTrackPos)
            lDcaPosToPrimVertex = TMath::Abs(myTrackPos->GetD(lPrimaryVtxPosition[0], lPrimaryVtxPosition[1], lMagneticField));

        if (myTrackNeg)
            lDcaNegToPrimVertex = TMath::Abs(myTrackNeg->GetD(lPrimaryVtxPosition[0], lPrimaryVtxPosition[1], lMagneticField));

        if (TMath::Abs(lDcaPosToPrimVertex) < 0.06 || TMath::Abs(lDcaNegToPrimVertex) < 0.06)
            continue;

        // Armenteros variables:
        // lAlphaV0 = v0->AlphaV0();
        // lPtArmV0 = v0->PtArmV0();

        // Pseudorapidity:
        lV0Eta = v0->Eta();

        if ((((myTrackPos->GetTPCClusterInfo(2, 1)) < 70) || ((myTrackNeg->GetTPCClusterInfo(2, 1)) < 70)))
            continue;

        // GetKinkIndex condition
        if (myTrackPos->GetKinkIndex(0) < 0 || myTrackNeg->GetKinkIndex(0) < 0)
            continue;

        // Findable clusters > 0 condition
        if (myTrackPos->GetTPCNclsF() <= 0 || myTrackNeg->GetTPCNclsF() <= 0)
            continue;

        Float_t fTreeVariableNegEta = myTrackNeg->Eta();
        Float_t fTreeVariablePosEta = myTrackPos->Eta();

        if (TMath::Abs(fTreeVariableNegEta) > 0.8 || TMath::Abs(fTreeVariablePosEta) > 0.8)
            continue;

        nSigmaPosPion = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackPos, AliPID::kPion));
        nSigmaNegPion = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackNeg, AliPID::kPion));
        nSigmaPosProton = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackPos, AliPID::kProton));
        nSigmaNegProton = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackNeg, AliPID::kProton));
        //////////////////////////////////////////////////////////////////////////
        // Invariant mass
        v0->ChangeMassHypothesis(310);
        lInvMassK0s = v0->GetEffMass();
        if (TMath::Abs(lInvMassK0s - 0.497611) < 0.2)
        {
            if (lV0cosPointAngle < 0.97)
                continue;
            if (TMath::Abs(nSigmaPosPion) > 6 && TMath::Abs(nSigmaNegPion) > 6)
                continue;
            if (v0->GetDcaV0Daughters() > 1)
                continue;
            lPtK0s = v0->Pt();
            lPzK0s = v0->Pz();
            if (lPtK0s == 0)
                continue;
            lRapK0s = v0->Y(310);
            if (TMath::Abs(lRapK0s) > 0.5)
                continue;

            Float_t ctauK0s = lV0DecayLength * lInvMassK0s / v0->P();
            if (ctauK0s > 20)
                continue;
            ((TH1D *)(fOutputList->FindObject("hinvmassK0s")))->Fill(lInvMassK0s);
        }
        Float_t massLambda = 1.11568;

        v0->ChangeMassHypothesis(3122);
        lInvMassLambda = v0->GetEffMass();
        if (TMath::Abs(lInvMassLambda - massLambda) < 0.2)
        {
            if (lV0cosPointAngle < 0.995)
                continue;
            if (TMath::Abs(nSigmaPosProton) > 6 && TMath::Abs(nSigmaNegPion) > 6)
                continue;
            if (v0->GetDcaV0Daughters() > 1)
                continue;

            lPtLambda = v0->Pt();
            lPzLambda = v0->Pz();
            if (lPtLambda == 0)
                continue;
            lRapLambda = v0->Y(3122);
            if (TMath::Abs(lRapLambda) > 0.5)
                continue;
            Float_t ctauLambda = lV0DecayLength * lInvMassLambda / v0->P();
            if (ctauLambda > 30)
                continue;
            ((TH1D *)(fOutputList->FindObject("hinvmassLambda")))->Fill(lInvMassLambda);
        }

        v0->ChangeMassHypothesis(-3122);

        lInvMassAntiLambda = v0->GetEffMass();

        if (TMath::Abs(lInvMassAntiLambda - massLambda) < 0.2)
        {
            if (lV0cosPointAngle < 0.995)
                continue;
            if (TMath::Abs(nSigmaPosPion) > 6 && TMath::Abs(nSigmaNegProton) > 6)
                continue;
            if (v0->GetDcaV0Daughters() > 1)
                continue;

            lPtAntiLambda = v0->Pt();
            lPzAntiLambda = v0->Pz();

            if (lPtAntiLambda == 0)
                continue;
            lRapAntiLambda = v0->Y(-3122);
            if (TMath::Abs(lRapAntiLambda) > 0.5)
                continue;
            Float_t ctauAntiLambda = lV0DecayLength * lInvMassAntiLambda / v0->P();
            if (ctauAntiLambda > 30)
                continue;
            ((TH1D *)(fOutputList->FindObject("hinvmassAntiLambda")))->Fill(lInvMassAntiLambda);
        }
    }

    // if (TestTrackCuts)
    //     delete TestTrackCuts;
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityLambdaK0s::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskFlatenicityLambdaK0s::GetFlatenicityV0()
{
    AliVVZERO *lVV0 = 0x0;
    AliVEvent *lVevent = 0x0;
    lVevent = dynamic_cast<AliVEvent *>(InputEvent());
    if (!lVevent)
    {
        AliWarning("ERROR: ESD / AOD event not available \n");
        return -1;
    }
    // Get VZERO Information for multiplicity later
    lVV0 = lVevent->GetVZEROData();
    if (!lVV0)
    {
        AliError("AliVVZERO not available");
        return 9999;
    }
    // Flatenicity calculation
    const Int_t nRings = 4;
    const Int_t nSectors = 8;
    Float_t minEtaV0C[nRings] = {-3.7, -3.2, -2.7, -2.2};
    Float_t maxEtaV0C[nRings] = {-3.2, -2.7, -2.2, -1.7};
    Float_t maxEtaV0A[nRings] = {5.1, 4.5, 3.9, 3.4};
    Float_t minEtaV0A[nRings] = {4.5, 3.9, 3.4, 2.8};
    // Grid
    const Int_t nCells = nRings * 2 * nSectors;
    Float_t RhoLattice[nCells];
    for (Int_t iCh = 0; iCh < nCells; iCh++)
    {
        RhoLattice[iCh] = 0.0;
    }

    Int_t nringA = 0;
    Int_t nringC = 0;
    for (Int_t iCh = 0; iCh < nCells; iCh++)
    {
        Float_t detaV0 = -1;
        Float_t mult = lVV0->GetMultiplicity(iCh);
        if (iCh < 32)
        { // V0C
            if (iCh < 8)
            {
                nringC = 0;
            }
            else if (iCh >= 8 && iCh < 16)
            {
                nringC = 1;
            }
            else if (iCh >= 16 && iCh < 24)
            {
                nringC = 2;
            }
            else
            {
                nringC = 3;
            }
            detaV0 = maxEtaV0C[nringC] - minEtaV0C[nringC];
        }
        else
        { // V0A
            if (iCh < 40)
            {
                nringA = 0;
            }
            else if (iCh >= 40 && iCh < 48)
            {
                nringA = 1;
            }
            else if (iCh >= 48 && iCh < 56)
            {
                nringA = 2;
            }
            else
            {
                nringA = 3;
            }
            detaV0 = maxEtaV0A[nringA] - minEtaV0A[nringA];
        }
        RhoLattice[iCh] =
            mult / detaV0; // needed to consider the different eta coverage
    }
    // Filling histos with mult info
    // for (Int_t iCh = 0; iCh < nCells; iCh++) {
    // hActivityV0DataSect->Fill(iCh, RhoLattice[iCh]);
    //}
    Float_t mRho = 0;
    Float_t flatenicity = -1;
    for (Int_t iCh = 0; iCh < nCells; iCh++)
    {
        mRho += RhoLattice[iCh];
    }
    Float_t multiplicityV0M = mRho;
    // average activity per cell
    mRho /= (1.0 * nCells);
    // get sigma
    Double_t sRho_tmp = 0;
    for (Int_t iCh = 0; iCh < nCells; iCh++)
    {
        sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
    }
    sRho_tmp /= (1.0 * nCells * nCells);
    Float_t sRho = TMath::Sqrt(sRho_tmp);
    Bool_t fRemoveTrivialScaling = kFALSE;
    if (mRho > 0)
    {
        if (fRemoveTrivialScaling)
        {
            flatenicity = TMath::Sqrt(multiplicityV0M) * sRho / mRho;
        }
        else
        {
            flatenicity = sRho / mRho;
        }
    }
    else
    {
        flatenicity = 9999;
    }
    return flatenicity;
}
