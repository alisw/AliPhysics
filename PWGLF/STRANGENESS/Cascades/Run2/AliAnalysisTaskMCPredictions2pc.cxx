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
// This task is meant to acquire MC-level predictions for 2pc analyses
//
// Please report any bugs, complaints, suggestions to:
// --- david.dobrigkeit.chinellato@cern.ch
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
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
//#include "AliLog.h"

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
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliPPVsMultUtils.h"
#include "AliPWG0Helper.h"
#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskMCPredictions2pc.h"
////
// #include "AliMultEstimator.h"
// #include "AliMultVariable.h"
// #include "AliMultInput.h"
// #include "AliMultSelection.h"

#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenPythiaEventHeader.h"
//#include "AliPythia8.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMCPredictions2pc)

AliAnalysisTaskMCPredictions2pc::AliAnalysisTaskMCPredictions2pc()
: AliAnalysisTaskSE(),
fListHist(0),

// relevant for task configuration
fSmallMultRange(1000),
fLargeMultRange(2000),
fRebinFactor(1), 
fkNIntervals(1),
fkSelectINELgtZERO(kTRUE),
fkMinimumMultiplicity(-1),
fkDo2pc(kTRUE),
fkMinPtTrigger(2.0),
fkMaxPtTrigger(4.0),
fkMinEta(-0.8),
fkMaxEta(+0.8),
fkNEtaBins(80),

// Histograms
fHistForwardMult(0),
fHistNchVsForwardMult(0),
fHistEventCounter(0),
fHistChargedEta(0),

// 2pc histograms
fHistEtaVsPtTrigger(0x0), 

fHist3d2pcK0Short(0x0),
fHist3d2pcLambda(0x0),
fHist3d2pcAntiLambda(0x0),
fHist3d2pcXiMinus(0x0),
fHist3d2pcXiPlus(0x0),
fHist3d2pcOmegaMinus(0x0),
fHist3d2pcOmegaPlus(0x0), 

fHist3d2pcMixedK0Short(0x0),
fHist3d2pcMixedLambda(0x0),
fHist3d2pcMixedAntiLambda(0x0),
fHist3d2pcMixedXiMinus(0x0),
fHist3d2pcMixedXiPlus(0x0),
fHist3d2pcMixedOmegaMinus(0x0),
fHist3d2pcMixedOmegaPlus(0x0)

{
  for(Int_t ii=0; ii<50; ii++){
    fEMBufferEta[ii]=0;
    fEMBufferPhi[ii]=0;
  }
  for(Int_t ii=0; ii<10; ii++){
    fkIntervalMinEta[ii]=0;
    fkIntervalMaxEta[ii]=0;
  }
  fkIntervalMinEta[0] = -1.4;
  fkIntervalMaxEta[0] = +1.4;
}

AliAnalysisTaskMCPredictions2pc::AliAnalysisTaskMCPredictions2pc(const char *name, Int_t lNSmallBinning, Int_t lNLargeBinning, Int_t lRebinFactor, Int_t lNEtaBins)
: AliAnalysisTaskSE(name),
fListHist(0),

// relevant for task configuration
fSmallMultRange(1000),
fLargeMultRange(2000),
fRebinFactor(lRebinFactor), 
fkNIntervals(1),
fkSelectINELgtZERO(kTRUE),
fkMinimumMultiplicity(-1),
fkDo2pc(kTRUE),
fkMinPtTrigger(2.0),
fkMaxPtTrigger(4.0),
fkMinEta(-0.8),
fkMaxEta(+0.8),
fkNEtaBins(lNEtaBins),

// Histograms
fHistForwardMult(0),
fHistNchVsForwardMult(0),
fHistEventCounter(0),
fHistChargedEta(0),

// 2pc histograms
fHistEtaVsPtTrigger(0x0), 

fHist3d2pcK0Short(0x0),
fHist3d2pcLambda(0x0),
fHist3d2pcAntiLambda(0x0),
fHist3d2pcXiMinus(0x0),
fHist3d2pcXiPlus(0x0),
fHist3d2pcOmegaMinus(0x0),
fHist3d2pcOmegaPlus(0x0), 

fHist3d2pcMixedK0Short(0x0),
fHist3d2pcMixedLambda(0x0),
fHist3d2pcMixedAntiLambda(0x0),
fHist3d2pcMixedXiMinus(0x0),
fHist3d2pcMixedXiPlus(0x0),
fHist3d2pcMixedOmegaMinus(0x0),
fHist3d2pcMixedOmegaPlus(0x0)

{
  for(Int_t ii=0; ii<50; ii++){
    fEMBufferEta[ii]=0;
    fEMBufferPhi[ii]=0;
  }
  for(Int_t ii=0; ii<10; ii++){
    fkIntervalMinEta[ii]=0;
    fkIntervalMaxEta[ii]=0;
  }
  fkIntervalMinEta[0] = -1.4;
  fkIntervalMaxEta[0] = +1.4;
}


AliAnalysisTaskMCPredictions2pc::~AliAnalysisTaskMCPredictions2pc()
{
  //------------------------------------------------
  // DESTRUCTOR
  //------------------------------------------------
  
  if (fListHist) {
    delete fListHist;
    fListHist = 0x0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskMCPredictions2pc::UserCreateOutputObjects()
{
  //------------------------------------------------
  // Histograms: Basic Analysis Output
  //------------------------------------------------
  // Create histograms
  fListHist = new TList();
  fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
  
  //Settings for transverse momentum
  Int_t lNPtBins = 250;
  Double_t lMaxPt = 25.0;
  
  Int_t lNEtaBins = fkNEtaBins;
  Double_t lMaxAbsEta = 4;
  
  //Settings for charged particle counters (integers!)
  Int_t lNNchBins = fSmallMultRange/fRebinFactor;
  Double_t lLowNchBound  = -0.5;
  Double_t lHighNchBound = -0.5 + ((double)(fSmallMultRange));
  
  Int_t lNNchBinsForward = fLargeMultRange/fRebinFactor;
  Double_t lLowNchBoundForward  = -0.5;
  Double_t lHighNchBoundForward = -0.5 + ((double)(fLargeMultRange));
  
  if(! fHistEventCounter ) {
    //Histogram Output: Event-by-Event
    fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",1,0,1);
    //Keeps track of some basics
    fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
    fListHist->Add(fHistEventCounter);
  }
  //___________________________________________________
  if(! fHistForwardMult ) {
    //Histogram Output: Event-by-Event
    fHistForwardMult = new TH1D( "fHistForwardMult", ";Forward Mult;Count",lNNchBinsForward,lLowNchBoundForward,lHighNchBoundForward);
    //Keeps track of some basics
    fListHist->Add(fHistForwardMult);
  }
  if(! fHistNchVsForwardMult ) {
    //Histogram Output: Event-by-Event
    fHistNchVsForwardMult = new TH2D( "fHistNchVsForwardMult", ";Forward Mult;Count",
                                 lNNchBinsForward,lLowNchBoundForward,lHighNchBoundForward,
                                 lNNchBins,lLowNchBound,lHighNchBound);
    //Keeps track of some basics
    fListHist->Add(fHistNchVsForwardMult);
  }
  //___________________________________________________
  if(! fHistEtaVsPtTrigger ) {
    //Histogram Output: Event-by-Event
    fHistEtaVsPtTrigger = new TH2D( "fHistPtTrigger", ";p_{T};Count",128, -0.8, 0.8, 200,0,20);
    fListHist->Add(fHistEtaVsPtTrigger);
  }
  //___________________________________________________
  if(! fHist3d2pcK0Short ) {
    fHist3d2pcK0Short = new TH3D("fHist3d2pcK0Short","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcK0Short);
  }
  if(! fHist3d2pcLambda ) {
    fHist3d2pcLambda = new TH3D("fHist3d2pcLambda","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcLambda);
  }
  if(! fHist3d2pcAntiLambda ) {
    fHist3d2pcAntiLambda = new TH3D("fHist3d2pcAntiLambda","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcAntiLambda);
  }
  if(! fHist3d2pcXiMinus ) {
    fHist3d2pcXiMinus = new TH3D("fHist3d2pcXiMinus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcXiMinus);
  }
  if(! fHist3d2pcXiPlus ) {
    fHist3d2pcXiPlus = new TH3D("fHist3d2pcXiPlus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcXiPlus);
  }
  if(! fHist3d2pcOmegaMinus ) {
    fHist3d2pcOmegaMinus = new TH3D("fHist3d2pcOmegaMinus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcOmegaMinus);
  }
  if(! fHist3d2pcOmegaPlus ) {
    fHist3d2pcOmegaPlus = new TH3D("fHist3d2pcOmegaPlus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcOmegaPlus);
  }
  //___________________________________________________
  if(! fHist3d2pcMixedK0Short ) {
    fHist3d2pcMixedK0Short = new TH3D("fHist3d2pcMixedK0Short","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcMixedK0Short);
  }
  if(! fHist3d2pcMixedLambda ) {
    fHist3d2pcMixedLambda = new TH3D("fHist3d2pcMixedLambda","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcMixedLambda);
  }
  if(! fHist3d2pcMixedAntiLambda ) {
    fHist3d2pcMixedAntiLambda = new TH3D("fHist3d2pcMixedAntiLambda","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcMixedAntiLambda);
  }
  if(! fHist3d2pcMixedXiMinus ) {
    fHist3d2pcMixedXiMinus = new TH3D("fHist3d2pcMixedXiMinus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcMixedXiMinus);
  }
  if(! fHist3d2pcMixedXiPlus ) {
    fHist3d2pcMixedXiPlus = new TH3D("fHist3d2pcMixedXiPlus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcMixedXiPlus);
  }
  if(! fHist3d2pcMixedOmegaMinus ) {
    fHist3d2pcMixedOmegaMinus = new TH3D("fHist3d2pcMixedOmegaMinus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcMixedOmegaMinus);
  }
  if(! fHist3d2pcMixedOmegaPlus ) {
    fHist3d2pcMixedOmegaPlus = new TH3D("fHist3d2pcMixedOmegaPlus","",2*lNEtaBins,-2*lMaxAbsEta,+2*lMaxAbsEta,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),20,0,10);
    fListHist->Add(fHist3d2pcMixedOmegaPlus);
  }

  //List of Histograms: Normal
  PostData(1, fListHist);  
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskMCPredictions2pc::UserExec(Option_t *)
{
  // Main loop
  // Called for each event
  
  AliMCEvent  *lMCevent  = 0x0;
  AliStack    *lMCstack  = 0x0;
  
  // Appropriate for ESD analysis!
  
  lMCevent = MCEvent();
  if (!lMCevent) {
    Printf("ERROR: Could not retrieve MC event \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }
  
  lMCstack = lMCevent->Stack();
  if (!lMCstack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }
  
  //------------------------------------------------
  // Multiplicity Information Acquistion
  //------------------------------------------------
  
  //Monte Carlo Level information !
  //--------- GENERATED NUMBER OF CHARGED PARTICLES
  // ---> Variable Definition
  
  Long_t lNchEta5   = 0;
  Long_t lNchEta8   = 0;
  
  //this keeps multiplicity over a wide range
  //(multiple intervals as configured by the fkInterval... vars) 
  Double_t lNchForward  = 0;
  
  Bool_t lEvSel_INELgtZEROStackPrimaries=kFALSE;
  
  
  //----- Loop on Stack ----------------------------------------------------------------
  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
  {   // This is the begining of the loop on tracks
    TParticle* particleOne = lMCstack->Particle(iCurrentLabelStack);
    if(!particleOne) continue;
    if(!particleOne->GetPDG()) continue;
    Double_t lThisCharge = particleOne->GetPDG()->Charge()/3.;
    if(TMath::Abs(lThisCharge)<0.001) continue;
    if(! (lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) ) continue;
    
    Double_t gpt = particleOne -> Pt();
    Double_t geta = particleOne -> Eta();
    
    if( TMath::Abs(geta) < 0.5 ) lNchEta5++;
    if( TMath::Abs(geta) < 0.8 ) lNchEta8++;
    if( TMath::Abs(geta) < 1.0 ) lEvSel_INELgtZEROStackPrimaries = kTRUE;
    
    //Special treatment: multiple intervals 
    for(Int_t ii = 0; ii<fkNIntervals; ii++ )
      if( fkIntervalMinEta[ii] < geta && geta < fkIntervalMaxEta[ii] ) lNchForward+=fkIntervalWeight[ii];
  }//End of loop on tracks
  //----- End Loop on Stack ------------------------------------------------------------
  
  //Reject non-INEL>0 if requested
  if( !lEvSel_INELgtZEROStackPrimaries && fkSelectINELgtZERO ) return;
  
  if( lNchForward < fkMinimumMultiplicity+1e-10 ) return;
  
  //------------------------------------------------
  // Fill Event Counters
  //------------------------------------------------
  
  //Basics: All Processed
  if( !fHistEventCounter ) {
    Printf("ERROR: Could not retrieve fHistEventCounter! This will crash!\n");
  }
  fHistEventCounter->Fill(0.5);
  
  if(fHistForwardMult)      fHistForwardMult        -> Fill ( lNchForward );
  if(fHistNchVsForwardMult) fHistNchVsForwardMult   -> Fill ( lNchForward, lNchEta5  );
  
  //------------------------------------------------
  // Fill Spectra as Needed
  //------------------------------------------------
  
  //===== Start 2pc =================
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // step 1: find particle indices to be correlated

  std::vector<uint32_t> triggerIndices;
  std::vector<std::vector<uint32_t>> associatedIndices;

  std::vector<uint32_t> k0ShortIndices;
  std::vector<uint32_t> lambdaIndices;
  std::vector<uint32_t> antiLambdaIndices;
  std::vector<uint32_t> xiMinusIndices;
  std::vector<uint32_t> xiPlusIndices;
  std::vector<uint32_t> omegaMinusIndices;
  std::vector<uint32_t> omegaPlusIndices;

  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
  {
    // Determine if within acceptance, otherwise fully reject from list
    // done such that this check is done O(N) and not O(N^2)
    TParticle* lThisParticle = lMCstack->Particle(iCurrentLabelStack);
    AliMCParticle* lMCPart = (AliMCParticle*)lMCevent->GetTrack(iCurrentLabelStack);
    if(!lThisParticle) continue;
    Bool_t lIsPhysicalPrimary = lMCstack->IsPhysicalPrimary(iCurrentLabelStack);
    Double_t geta = lThisParticle -> Eta();
    Double_t gpt = lThisParticle -> Pt();

    // kick out stuff not at midrapidity
    if( ( geta < fkMinEta || geta > fkMaxEta) ) continue; 
    
    if(!lThisParticle->GetPDG()) continue;
    Double_t lThisCharge = lThisParticle->GetPDG()->Charge()/3.;

    // populate triggers with charged particles within desired pT window
    if (fkMinPtTrigger<gpt && gpt<fkMaxPtTrigger && TMath::Abs(lThisCharge)>0.001){
      //valid trigger 
      fHistEtaVsPtTrigger->Fill(gpt, geta);
      triggerIndices.emplace_back(iCurrentLabelStack);
    }

    if ( lThisParticle->GetPdgCode() ==  310 ) k0ShortIndices.emplace_back(iCurrentLabelStack);
    if ( lThisParticle->GetPdgCode() ==  3122 ) lambdaIndices.emplace_back(iCurrentLabelStack);
    if ( lThisParticle->GetPdgCode() == -3122 ) antiLambdaIndices.emplace_back(iCurrentLabelStack);
    if ( lThisParticle->GetPdgCode() ==  3312 ) xiMinusIndices.emplace_back(iCurrentLabelStack);
    if ( lThisParticle->GetPdgCode() == -3312 ) xiPlusIndices.emplace_back(iCurrentLabelStack);
    if ( lThisParticle->GetPdgCode() ==  3334 ) omegaMinusIndices.emplace_back(iCurrentLabelStack);
    if ( lThisParticle->GetPdgCode() == -3334 ) omegaPlusIndices.emplace_back(iCurrentLabelStack);
  }

  associatedIndices.emplace_back(k0ShortIndices);
  associatedIndices.emplace_back(lambdaIndices);
  associatedIndices.emplace_back(antiLambdaIndices);
  associatedIndices.emplace_back(xiMinusIndices);
  associatedIndices.emplace_back(xiPlusIndices);
  associatedIndices.emplace_back(omegaMinusIndices);
  associatedIndices.emplace_back(omegaPlusIndices);

  TH3D *h3dAssociated[7] = {fHist3d2pcK0Short, fHist3d2pcLambda, fHist3d2pcAntiLambda, fHist3d2pcXiMinus, fHist3d2pcXiPlus, fHist3d2pcOmegaMinus, fHist3d2pcOmegaPlus};
  TH3D *h3dMixed[7] = {fHist3d2pcMixedK0Short, fHist3d2pcMixedLambda, fHist3d2pcMixedAntiLambda, fHist3d2pcMixedXiMinus, fHist3d2pcMixedXiPlus, fHist3d2pcMixedOmegaMinus, fHist3d2pcMixedOmegaPlus};

  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  //  Actually correlate stuff with stuff
  for (Int_t iTrigger = 0;  iTrigger < triggerIndices.size(); iTrigger++){   // trigger loop
    TParticle* lTriggerParticle = lMCstack->Particle(triggerIndices[iTrigger]);
    
    Double_t geta = lTriggerParticle -> Eta();
    Double_t gphi = lTriggerParticle -> Phi();

    for (Int_t iassocSpecies = 0;  iassocSpecies < associatedIndices.size(); iassocSpecies++){   // associated loop
      for (Int_t iassoc = 0;  iassoc < associatedIndices[iassocSpecies].size(); iassoc++){   // associated loop
        TParticle* lAssociatedParticle = 0x0;
        lAssociatedParticle = lMCstack->Particle( associatedIndices[iassocSpecies][iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        Double_t lThisPt    = lAssociatedParticle->Pt();
        h3dAssociated[iassocSpecies]->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
      }
    }
  }
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
  //Event mixing loop
  if( fEMBufferFull ){ //require also that a trigger exists
    for (Int_t iTrigger = 0;  iTrigger < fEMBufferSize; iTrigger++){   // trigger loop
      Double_t geta = fEMBufferEta[iTrigger]; //from previous events
      Double_t gphi = fEMBufferPhi[iTrigger]; //from previous events

      for (Int_t iassocSpecies = 0;  iassocSpecies < associatedIndices.size(); iassocSpecies++){   // associated loop
        for (Int_t iassoc = 0;  iassoc < associatedIndices[iassocSpecies].size(); iassoc++){   // associated loop
          TParticle* lAssociatedParticle = 0x0;
          lAssociatedParticle = lMCstack->Particle( associatedIndices[iassocSpecies][iassoc] );
          if(!lAssociatedParticle) {
            Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
            continue;
          }
          Double_t geta2 = lAssociatedParticle -> Eta();
          Double_t gphi2 = lAssociatedParticle -> Phi();
          Double_t lThisPt    = lAssociatedParticle->Pt();
          h3dMixed[iassocSpecies]->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt);
        }
      }
    }
  }
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // fill EM buffer
  if( fEMBufferSize > 0 ){
    for (Int_t iTrigger = 0;  iTrigger < triggerIndices.size(); iTrigger++){   // trigger loop
      TParticle* lThisParticle = lMCstack->Particle(triggerIndices[iTrigger]);
      AliMCParticle* lMCPart = (AliMCParticle*)lMCevent->GetTrack(triggerIndices[iTrigger]);
      if(!lThisParticle) continue;

      //Add to buffer
      fEMBufferEta[ fEMBufferCycle ] = lThisParticle->Eta();
      fEMBufferPhi[ fEMBufferCycle ] = lThisParticle->Phi();
      fEMBufferCycle++;
      if(fEMBufferCycle>=fEMBufferSize) fEMBufferFull = kTRUE;
      fEMBufferCycle = fEMBufferCycle%fEMBufferSize;
    }
  }
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  
  // Post output data.
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskMCPredictions2pc::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
  TList *cRetrievedList = 0x0;
  cRetrievedList = (TList*)GetOutputData(1);
  if(!cRetrievedList) {
    Printf("ERROR - AliAnalysisTaskMCPredictions2pc : ouput data container list not available\n");
    return;
  }
  
  fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
  if (!fHistEventCounter) {
    Printf("ERROR - AliAnalysisTaskMCPredictions2pc : fHistEventCounter not available");
    return;
  }
  
  TCanvas *canCheck = new TCanvas("AliAnalysisTaskMCPredictions2pc","Event Multiplicity",10,10,510,510);
  canCheck->cd(1)->SetLogy();
  
  fHistEventCounter->SetMarkerStyle(22);
  fHistEventCounter->DrawCopy("E");
}

//______________________________________________________________________
Double_t AliAnalysisTaskMCPredictions2pc::MyRapidity(Double_t rE, Double_t rPz) const
{
  // Local calculation for rapidity
  Double_t ReturnValue = -100;
  if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
    ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
  }
  return ReturnValue;
}

//______________________________________________________________________
Double_t AliAnalysisTaskMCPredictions2pc::ComputeDeltaPhi( Double_t phi1, Double_t phi2) const {
  //To be completely sure, use inner products
  Double_t x1, y1, x2, y2;
  x1 = TMath::Cos( phi1 );
  y1 = TMath::Sin( phi1 );
  x2 = TMath::Cos( phi2 );
  y2 = TMath::Sin( phi2 );
  Double_t lInnerProd  = x1*x2 + y1*y2;
  Double_t lVectorProd = x1*y2 - x2*y1;
  
  Double_t lReturnVal = 0;
  if( lVectorProd > 1e-8 ){
    lReturnVal = TMath::ACos(lInnerProd);
  }
  if( lVectorProd < -1e-8 ){
    lReturnVal = -TMath::ACos(lInnerProd);
  }
  if( lReturnVal < -TMath::Pi()/2 ) lReturnVal += 2*TMath::Pi();
  return lReturnVal;
}

//______________________________________________________________________
void AliAnalysisTaskMCPredictions2pc::PrintEtaIntervals(){
  for(Int_t ii=0; ii<fkNIntervals; ii++)
    cout<<"Interval #"<<ii<<"\tmin = "<<fkIntervalMinEta[ii]<<"\tmax = "<<fkIntervalMaxEta[ii]<<"\tweight = "<<fkIntervalWeight[ii]<<endl;
}