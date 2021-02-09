/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Modified version of AliAnalysisTaskCheckCascade.cxx.
// This is a 'hybrid' output version, in that it uses a classic TTree
// ROOT object to store the candidates, plus a couple of histograms filled on
// a per-event basis for storing variables too numerous to put in a tree. 
//
// --- Added bits of code for checking V0s 
//      (from AliAnalysisTaskCheckStrange.cxx)
//
//  --- Algorithm Description 
//   1. Perform Physics Selection
//   2. Perform Primary Vertex |z|<10cm selection
//   3. Perform Primary Vertex NoTPCOnly vertexing selection (>0 contrib.)
//   4. Perform Pileup Rejection
//   5. Analysis Loops: 
//    5a. Fill TTree object with V0 information, candidates
//
//  Please Report Any Bugs! 
//
//   --- David Dobrigkeit Chinellato
//        (david.chinellato@gmail.com)
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TTree;
class TParticle;
class TVector3;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "AliLog.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"

#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"

#include "AliCFContainer.h"
#include "AliMultiplicity.h"

#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliESDHeader.h"
#include "AliAODTrack.h"
#include "AliAnalysisTaskhStrCorr.h"

//debugging purposes
#include "TObjectTable.h"

ClassImp(AliAnalysisTaskhStrCorr)

AliAnalysisTaskhStrCorr::AliAnalysisTaskhStrCorr() 
: AliAnalysisTaskSE(),
//Output lists
fOutput(0),

fHistTriggerPt(0),

fHistK0ShortMassVsPt(0),
fHistLambdaMassVsPt(0),
fHistAntiLambdaMassVsPt(0),
fHistXiMinusMassVsPt(0),
fHistXiPlusMassVsPt(0),
fHistOmegaMinusMassVsPt(0),
fHistOmegaPlusMassVsPt(0),

//Particle correlation analysis: same-event correlation functions
fHistCorrFuncPeakK0ShortVsPt(0),
fHistCorrFuncPeakLambdaVsPt(0),
fHistCorrFuncPeakAntiLambdaVsPt(0),
fHistCorrFuncPeakXiMinusVsPt(0),
fHistCorrFuncPeakXiPlusVsPt(0),
fHistCorrFuncPeakOmegaMinusVsPt(0),
fHistCorrFuncPeakOmegaPlusVsPt(0),

fHistCorrFuncSideBandK0ShortVsPt(0),
fHistCorrFuncSideBandLambdaVsPt(0),
fHistCorrFuncSideBandAntiLambdaVsPt(0),
fHistCorrFuncSideBandXiMinusVsPt(0),
fHistCorrFuncSideBandXiPlusVsPt(0),
fHistCorrFuncSideBandOmegaMinusVsPt(0),
fHistCorrFuncSideBandOmegaPlusVsPt(0),

//Functions to calculate relevant peak info for stat treatment of 2pc
fParametricK0ShortMean(0),
fParametricK0ShortSigma(0),
fParametricLambdaMean(0),
fParametricLambdaSigma(0),
fParametricXiMean(0),
fParametricXiSigma(0),

//Task Control / Utils
fPIDResponse(0)
{
    // Dummy Constructor
    for(Int_t iV0selIdx   = 0; iV0selIdx   < 7; iV0selIdx++   ) { fV0Sels          [iV0selIdx   ] = -1.; }
}

AliAnalysisTaskhStrCorr::AliAnalysisTaskhStrCorr(const char *name) 
: AliAnalysisTaskSE(name),
//Output lists
fOutput(0),

fHistTriggerPt(0),

fHistK0ShortMassVsPt(0),
fHistLambdaMassVsPt(0),
fHistAntiLambdaMassVsPt(0),
fHistXiMinusMassVsPt(0),
fHistXiPlusMassVsPt(0),
fHistOmegaMinusMassVsPt(0),
fHistOmegaPlusMassVsPt(0),

//Particle correlation analysis: same-event correlation functions
fHistCorrFuncPeakK0ShortVsPt(0),
fHistCorrFuncPeakLambdaVsPt(0),
fHistCorrFuncPeakAntiLambdaVsPt(0),
fHistCorrFuncPeakXiMinusVsPt(0),
fHistCorrFuncPeakXiPlusVsPt(0),
fHistCorrFuncPeakOmegaMinusVsPt(0),
fHistCorrFuncPeakOmegaPlusVsPt(0),

fHistCorrFuncSideBandK0ShortVsPt(0),
fHistCorrFuncSideBandLambdaVsPt(0),
fHistCorrFuncSideBandAntiLambdaVsPt(0),
fHistCorrFuncSideBandXiMinusVsPt(0),
fHistCorrFuncSideBandXiPlusVsPt(0),
fHistCorrFuncSideBandOmegaMinusVsPt(0),
fHistCorrFuncSideBandOmegaPlusVsPt(0),

//Functions to calculate relevant peak info for stat treatment of 2pc
fParametricK0ShortMean(0),
fParametricK0ShortSigma(0),
fParametricLambdaMean(0),
fParametricLambdaSigma(0),
fParametricXiMean(0),
fParametricXiSigma(0),

//Task Control / Utils
fPIDResponse(0)
{
    // Constructor
    // V0 Selection cuts as per default (pp-like but based on K0Short rather than Lambda)
    fV0Sels[0] =  33.  ;  // never used in ALICE, actually. Boh
    fV0Sels[1] =   0.04;  // DCA neg/pos
    fV0Sels[2] =   0.04;  // DCA neg/pos
    fV0Sels[3] =   1.5 ;  // DCA Daughters (sigh, effective centimeters)
    fV0Sels[4] =   0.97;  // Cosine of pointing angle and all that jazz
    fV0Sels[5] =    0.5;  // R2D: fine for now
    fV0Sels[6] = 200.  ;  // R2D-max, but beware finder
    //Cascades
    fCascSels[0] =  33.;     // never used
    fCascSels[1] =   0.07;   // min allowed V0 impact parameter
    fCascSels[2] =   0.008;  // "window" around the Lambda mass
    fCascSels[3] =   0.05;   // min allowed bachelor's impact parameter
    fCascSels[4] =   1.5;    // max allowed DCA between the V0 and the bachelor
    fCascSels[5] =   0.97;   // min allowed cosine of the cascade pointing angle
    fCascSels[6] =   0.8;    // min radius of the fiducial volume
    fCascSels[7] =   100.;     // max radius of the fiducial volume
    
    
    // Output slot #0 writes into a TList container (Lambda Histos and fTree)
    DefineOutput(1, TList::Class());
}


AliAnalysisTaskhStrCorr::~AliAnalysisTaskhStrCorr()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------
    
    if (fOutput){
        delete fOutput;
        fOutput = 0x0;
    }
}



//________________________________________________________________________
void AliAnalysisTaskhStrCorr::UserCreateOutputObjects()
{
    
    //Define Output Lists
    fOutput = new TList();
    fOutput->SetOwner();
    
    
    //Histogram Output: Event-by-Event
    fHistEvent = new TH1D( "fHistEvent", ";Evt. Sel. Step;Count",4,0,4);
    fHistEvent->GetXaxis()->SetBinLabel(1, "Processed");
    fHistEvent->GetXaxis()->SetBinLabel(2, "Phys-Sel");
    fHistEvent->GetXaxis()->SetBinLabel(3, "Has Vtx");
    fHistEvent->GetXaxis()->SetBinLabel(4, "Vtx |z|<10cm");
    fOutput->Add(fHistEvent);
    
    Double_t lK0Mass = 0.497, lLamMass = 1.116, lXiMass = 1.322, lOmMass = 1.672;
    Double_t lK0Win = 0.150, lLamWin = 0.100, lXiWin = 0.100, lOmWin = 0.100;
    Double_t lPtMax = 20; Int_t lNPtBins = 200;
    
    //Create invariant mass histograms
    if(!fHistK0ShortMassVsPt){
        fHistK0ShortMassVsPt = new TH2D("fHistK0ShortMassVsPt", "", lNPtBins, 0, lPtMax, 600, lK0Mass-lK0Win, lK0Mass+lK0Win);
        fOutput->Add(fHistK0ShortMassVsPt);
    }
    if(!fHistLambdaMassVsPt){
        fHistLambdaMassVsPt = new TH2D("fHistLambdaMassVsPt", "", lNPtBins, 0, lPtMax, 400, lLamMass-lLamWin, lLamMass+lLamWin);
        fOutput->Add(fHistLambdaMassVsPt);
    }
    if(!fHistAntiLambdaMassVsPt){
        fHistAntiLambdaMassVsPt = new TH2D("fHistAntiLambdaMassVsPt", "", lNPtBins, 0, lPtMax, 400, lLamMass-lLamWin, lLamMass+lLamWin);
        fOutput->Add(fHistAntiLambdaMassVsPt);
    }
    if(!fHistXiMinusMassVsPt){
        fHistXiMinusMassVsPt = new TH2D("fHistXiMinusMassVsPt", "", lNPtBins, 0, lPtMax, 400, lXiMass-lXiWin, lXiMass+lXiWin);
        fOutput->Add(fHistXiMinusMassVsPt);
    }
    if(!fHistXiPlusMassVsPt){
        fHistXiPlusMassVsPt = new TH2D("fHistXiPlusMassVsPt", "", lNPtBins, 0, lPtMax, 400, lXiMass-lXiWin, lXiMass+lXiWin);
        fOutput->Add(fHistXiPlusMassVsPt);
    }
    if(!fHistOmegaMinusMassVsPt){
        fHistOmegaMinusMassVsPt = new TH2D("fHistOmegaMinusMassVsPt", "", lNPtBins, 0, lPtMax, 400, lOmMass-lOmWin, lOmMass+lOmWin);
        fOutput->Add(fHistOmegaMinusMassVsPt);
    }
    if(!fHistOmegaPlusMassVsPt){
        fHistOmegaPlusMassVsPt = new TH2D("fHistOmegaPlusMassVsPt", "", lNPtBins, 0, lPtMax, 400, lOmMass-lOmWin, lOmMass+lOmWin);
        fOutput->Add(fHistOmegaPlusMassVsPt);
    }
    
    //Create 2PC objects
    if(!fHistCorrFuncPeakK0ShortVsPt){
        fHistCorrFuncPeakK0ShortVsPt = new TH3D("fHistCorrFuncPeakK0ShortVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncPeakK0ShortVsPt);
    }
    if(!fHistCorrFuncPeakLambdaVsPt){
        fHistCorrFuncPeakLambdaVsPt = new TH3D("fHistCorrFuncPeakLambdaVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncPeakLambdaVsPt);
    }
    if(!fHistCorrFuncPeakAntiLambdaVsPt){
        fHistCorrFuncPeakAntiLambdaVsPt = new TH3D("fHistCorrFuncPeakAntiLambdaVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncPeakAntiLambdaVsPt);
    }
    if(!fHistCorrFuncPeakXiMinusVsPt){
        fHistCorrFuncPeakXiMinusVsPt = new TH3D("fHistCorrFuncPeakXiMinusVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncPeakXiMinusVsPt);
    }
    if(!fHistCorrFuncPeakXiPlusVsPt){
        fHistCorrFuncPeakXiPlusVsPt = new TH3D("fHistCorrFuncPeakXiPlusVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncPeakXiPlusVsPt);
    }
    if(!fHistCorrFuncPeakOmegaMinusVsPt){
        fHistCorrFuncPeakOmegaMinusVsPt = new TH3D("fHistCorrFuncPeakOmegaMinusVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncPeakOmegaMinusVsPt);
    }
    if(!fHistCorrFuncPeakOmegaPlusVsPt){
        fHistCorrFuncPeakOmegaPlusVsPt = new TH3D("fHistCorrFuncPeakOmegaPlusVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncPeakOmegaPlusVsPt);
    }
    
    if(!fHistCorrFuncSideBandK0ShortVsPt){
        fHistCorrFuncSideBandK0ShortVsPt = new TH3D("fHistCorrFuncSideBandK0ShortVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncSideBandK0ShortVsPt);
    }
    if(!fHistCorrFuncSideBandLambdaVsPt){
        fHistCorrFuncSideBandLambdaVsPt = new TH3D("fHistCorrFuncSideBandLambdaVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncSideBandLambdaVsPt);
    }
    if(!fHistCorrFuncSideBandAntiLambdaVsPt){
        fHistCorrFuncSideBandAntiLambdaVsPt = new TH3D("fHistCorrFuncSideBandAntiLambdaVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncSideBandAntiLambdaVsPt);
    }
    if(!fHistCorrFuncSideBandXiMinusVsPt){
        fHistCorrFuncSideBandXiMinusVsPt = new TH3D("fHistCorrFuncSideBandXiMinusVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncSideBandXiMinusVsPt);
    }
    if(!fHistCorrFuncSideBandXiPlusVsPt){
        fHistCorrFuncSideBandXiPlusVsPt = new TH3D("fHistCorrFuncSideBandXiPlusVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncSideBandXiPlusVsPt);
    }
    if(!fHistCorrFuncSideBandOmegaMinusVsPt){
        fHistCorrFuncSideBandOmegaMinusVsPt = new TH3D("fHistCorrFuncSideBandOmegaMinusVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncSideBandOmegaMinusVsPt);
    }
    if(!fHistCorrFuncSideBandOmegaPlusVsPt){
        fHistCorrFuncSideBandOmegaPlusVsPt = new TH3D("fHistCorrFuncSideBandOmegaPlusVsPt", "", lNPtBins, 0, lPtMax, 32, -1.6, 1.6, 32, -0.5*TMath::Pi(), +1.5*-0.5*TMath::Pi());
        fOutput->Add(fHistCorrFuncSideBandOmegaPlusVsPt);
    }
    
    //TF1's: don't initialize here.
    //User has to provide, will stream together with task
    
    //------------------------------------------------
    // Particle Identification Setup
    //------------------------------------------------
    
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    //Regular output: Histograms
    PostData(1, fOutput);
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskhStrCorr::UserExec(Option_t *) 
{
    
    // Main loop
    // Called for each event
    //gObjectTable->Print();
    AliAODEvent *lAODevent = 0x0;
    
    //AliAODEvent *lAODevent = 0x0;
    Int_t    nV0s                        = -1;
    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    Double_t lMagneticField                 = -10.;
    
    // Connect to the InputEvent
    // After these lines, we should have an ESD/AOD event + the number of cascades in it.
    
    // Appropriate for ESD analysis!
    
    lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
    if (!lAODevent) {
        AliWarning("ERROR: lAODevent not available from InputEvent() trying with AODEvent()");
        
        //  assume that the AOD is in the general output...
        lAODevent  = AODEvent();
        if(!lAODevent){
            AliWarning("ERROR: lAODevent not available from AODEvent() Aborting event!");
            return;
        }
    }
    fHistEvent->Fill(0.5);
    
    //------------------------------------------------
    // Physics Selection
    //------------------------------------------------
    
    // new method
    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
    
    //Standard Min-Bias Selection
    if ( ! isSelected ) {
        PostData(1, fOutput);
    }
    fHistEvent->Fill(1.5);
    
    //------------------------------------------------
    // After Trigger Selection
    //------------------------------------------------
    
    nV0s = lAODevent->GetNumberOfV0s();
    
    //------------------------------------------------
    // Getting: Primary Vertex + MagField Info
    //------------------------------------------------
    
    const AliAODVertex *lPrimaryBestAODVtx = lAODevent->GetPrimaryVertex();
    // get the best primary vertex available for the event
    // As done in AliCascadeVertexer, we keep the one which is the best one available.
    // between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
    // This one will be used for next calculations (DCA essentially)
    lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );
    
    Double_t tPrimaryVtxPosition[3];
    const AliVVertex *primaryVtx = lAODevent->GetPrimaryVertex();
    tPrimaryVtxPosition[0] = primaryVtx->GetX();
    tPrimaryVtxPosition[1] = primaryVtx->GetY();
    tPrimaryVtxPosition[2] = primaryVtx->GetZ();
    
    //------------------------------------------------
    // Only look at events with well-established PV
    //------------------------------------------------
    
    const AliAODVertex *lPrimaryTrackingAODVtxCheck = lAODevent->GetPrimaryVertex();
    const AliAODVertex *lPrimarySPDVtx = lAODevent->GetPrimaryVertexSPD();
    if (!lPrimarySPDVtx && !lPrimaryTrackingAODVtxCheck ){
        AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
        PostData(1, fOutput);
        return;
    }
    fHistEvent->Fill(2.5);
    
    //------------------------------------------------
    // Primary Vertex Z position: SKIP
    //------------------------------------------------
    
    if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0 ) {
        AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !");
        PostData(1, fOutput);
        return;
    }
    
    lMagneticField = lAODevent->GetMagneticField( );
    fHistEvent->Fill(3.5);
    
    //------------------------------------------------
    // MAIN LAMBDA LOOP STARTS HERE
    //------------------------------------------------
    
    //Variable definition
    Int_t    lOnFlyStatus = 0;// nv0sOn = 0, nv0sOff = 0;
    Double_t lChi2V0 = 0;
    Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
    Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
    Double_t lV0CosineOfPointingAngle = 0;
    Double_t lV0Radius = 0, lPt = 0;
    Double_t lRapK0Short = 0, lRapLambda = 0;
    Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
    Double_t lAlphaV0 = 0, lPtArmV0 = 0;
    
    Double_t fMinV0Pt = 0;
    Double_t fMaxV0Pt = 100;
    
    Int_t nv0s = 0;
    nv0s = lAODevent->GetNumberOfV0s();
    
    //Define candidate indices of interest
    TArrayI lIndK0Short(nv0s);
    TArrayI lIndLambda(nv0s);
    TArrayI lIndAntiLambda(nv0s);
    TArrayI lIndXiMinus(nv0s);
    TArrayI lIndXiPlus(nv0s);
    TArrayI lIndOmegaMinus(nv0s);
    TArrayI lIndOmegaPlus(nv0s);
    Long_t lNK0Short=0;
    Long_t lNLambda=0;
    Long_t lNAntiLambda=0;
    Long_t lNXiMinus=0;
    Long_t lNXiPlus=0;
    Long_t lNOmegaMinus=0;
    Long_t lNOmegaPlus=0;
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    // This loop establishes good V0 indices to be used for the correlation analysis
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
    for (Int_t iV0 = 0; iV0 < nv0s; iV0++)
    {// This is the begining of the V0 loop
        AliAODv0 *v0 = lAODevent->GetV0(iV0);
        if (!v0) continue;
        
        //Only use Offline Candidates for QA
        lOnFlyStatus = v0->GetOnFlyStatus();
        if( lOnFlyStatus == kTRUE ) continue;
        
        Double_t tDecayVertexV0[3]; v0->GetXYZ(tDecayVertexV0);
        Double_t tV0mom[3];
        v0->GetPxPyPz( tV0mom );
        Double_t lV0TotalMomentum = TMath::Sqrt(
                                                tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );
        
        lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);
        lPt = v0->Pt();
        lRapK0Short = v0->RapK0Short();
        lRapLambda  = v0->RapLambda();
        
        Double_t lMomPos[3]; //v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
        Double_t lMomNeg[3]; //v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);
        lMomPos[0] = v0->MomPosX();
        lMomPos[1] = v0->MomPosY();
        lMomPos[2] = v0->MomPosZ();
        lMomNeg[0] = v0->MomNegX();
        lMomNeg[1] = v0->MomNegY();
        lMomNeg[2] = v0->MomNegZ();
        
        AliAODTrack *pTrack=(AliAODTrack *)v0->GetDaughter(0); //0->Positive Daughter
        AliAODTrack *nTrack=(AliAODTrack *)v0->GetDaughter(1); //1->Negative Daughter
        if (!pTrack || !nTrack) {
            Printf("ERROR: Could not retrieve one of the daughter track");
            continue;
        }
        
        //Daughter Eta for Eta selection, afterwards
        Float_t lNegEta = nTrack->Eta();
        Float_t lPosEta = pTrack->Eta();
        
        // Filter like-sign V0 (next: add counter and distribution)
        if ( pTrack->Charge() == nTrack->Charge()){
            continue;
        }

        //________________________________________________________________________
        // Track quality cuts
        Float_t lPosTrackCrossedRows = pTrack->GetTPCClusterInfo(2,1);
        Float_t lNegTrackCrossedRows = nTrack->GetTPCClusterInfo(2,1);
        Int_t lLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
        if( lNegTrackCrossedRows < lLeastNbrCrossedRows )
            lLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;
        
        // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
        if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        
        if ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;
        
        //Findable clusters > 0 condition
        if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;
        
        //Compute ratio Crossed Rows / Findable clusters
        //Note: above test avoids division by zero!
        Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(pTrack->GetTPCNclsF()));
        Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(nTrack->GetTPCNclsF()));
        
        Float_t lLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
        if( lNegTrackCrossedRowsOverFindable < lLeastRatioCrossedRowsOverFindable )
            lLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;
        
        //Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
        if ( lLeastRatioCrossedRowsOverFindable < 0.8 ) continue;
        
        //End track Quality Cuts
        //________________________________________________________________________
        
        lDcaPosToPrimVertex = v0->DcaPosToPrimVertex();
        lDcaNegToPrimVertex = v0->DcaNegToPrimVertex();
        
        lOnFlyStatus = v0->GetOnFlyStatus();
        lChi2V0 = v0->Chi2V0();
        lDcaV0Daughters = v0->DcaV0Daughters();
        lDcaV0ToPrimVertex = v0->DcaV0ToPrimVertex();
        lV0CosineOfPointingAngle = v0->CosPointingAngle(tPrimaryVtxPosition);
        
        // Getting invariant mass infos directly from ESD
        lInvMassK0s        = v0->MassK0Short();
        lInvMassLambda     = v0->MassLambda();
        lInvMassAntiLambda = v0->MassAntiLambda();
        lAlphaV0 = v0->AlphaV0();
        lPtArmV0 = v0->PtArmV0();
        
        //Official means of acquiring N-sigmas
        Float_t lNSigmasPosProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
        Float_t lNSigmasPosPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
        Float_t lNSigmasNegProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
        Float_t lNSigmasNegPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );
        
        //This requires an Invariant Mass Hypothesis afterwards
        Float_t lDistOverTotMom = TMath::Sqrt(
                                              TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
                                              TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
                                              TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
                                              );
        lDistOverTotMom /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure
        
        //------------------------------------------------
        // Determine valid candidate index for later use
        // and fill invariant mass histograms
        //------------------------------------------------
        
        if( lDcaNegToPrimVertex > fV0Sels[1] && lDcaPosToPrimVertex > fV0Sels[2]      &&
           lDcaV0Daughters     < fV0Sels[3] && lV0CosineOfPointingAngle > fV0Sels[4] &&
           lV0Radius           > fV0Sels[5] && lV0Radius < fV0Sels [6] ){
            
            //Specific fV0Sel selection level, dE/dx applied
            if ( TMath::Abs(lNSigmasPosPion)   < 5 && TMath::Abs(lNSigmasNegPion)   < 5 ){
                lIndK0Short[lNK0Short++] = iV0;
                fHistK0ShortMassVsPt->Fill( lPt, lInvMassK0s );
            }
            if ( TMath::Abs(lNSigmasPosProton) < 5 && TMath::Abs(lNSigmasNegPion)   < 5 ){
                lIndK0Short[lNLambda++] = iV0;
                fHistLambdaMassVsPt->Fill( lPt, lInvMassLambda );
            }
            if ( TMath::Abs(lNSigmasPosPion)   < 5 && TMath::Abs(lNSigmasNegProton) < 5 ){
                lIndK0Short[lNAntiLambda++] = iV0;
                fHistAntiLambdaMassVsPt->Fill( lPt, lInvMassAntiLambda );
            }
        }
    }// This is the end of the V0 loop
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    // This loop establishes good cascade indices to be used for the corr analysis
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
    // Post output data.
    PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskhStrCorr::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query
    // This will draw the V0 candidate multiplicity, whose
    // number of entries corresponds to the number of triggered events.
    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList){
        Printf("ERROR - AliAnalysisTaskhStrCorr : ouput data container list not available\n");
        return;
    }
    fHistEvent = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEvent")  );
    if (!fHistEvent) {
        Printf("ERROR - AliAnalysisTaskhStrCorr : fHistEvent not available");
        return;
    }
    TCanvas *canCheck = new TCanvas("AliAnalysisTaskhStrCorr","V0 Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();
    fHistEvent->SetMarkerStyle(22);
    fHistEvent->DrawCopy("E");
}
