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
#include "THn.h"
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
fkNPhiBins(80),
fkVerboseMode(kFALSE),
fkDoEventMixing(kFALSE),

// Histograms
fHistForwardMult(0),
fHistNchVsForwardMult(0),
fHistEventCounter(0),
fHistChargedEta(0),

// 2pc histograms
fHist3dTrigger(0x0),
fHist3dAssoPions(0x0),
fHist3dAssoK0Short(0x0),
fHist3dAssoLambda(0x0),
fHist3dAssoAntiLambda(0x0),
fHist3dAssoXiMinus(0x0),
fHist3dAssoXiPlus(0x0),
fHist3dAssoOmegaMinus(0x0),
fHist3dAssoOmegaPlus(0x0),

fHist4d2pcPions(0x0),
fHist4d2pcK0Short(0x0),
fHist4d2pcLambda(0x0),
fHist4d2pcAntiLambda(0x0),
fHist4d2pcXiMinus(0x0),
fHist4d2pcXiPlus(0x0),
fHist4d2pcOmegaMinus(0x0),
fHist4d2pcOmegaPlus(0x0),

fHist4d2pcMixedPions(0x0),
fHist4d2pcMixedK0Short(0x0),
fHist4d2pcMixedLambda(0x0),
fHist4d2pcMixedAntiLambda(0x0),
fHist4d2pcMixedXiMinus(0x0),
fHist4d2pcMixedXiPlus(0x0),
fHist4d2pcMixedOmegaMinus(0x0),
fHist4d2pcMixedOmegaPlus(0x0),
fNMultBins(0),
fNPtBins(0)
{
  for(Int_t ii=0; ii<100; ii++){
    fMultBinBounds[ii] = 0.0;
    fPtBinBounds[ii] = 0.0;
  }
  for(Int_t ii=0; ii<20; ii++){
    fEMBufferFull[ii]=0;
  }
  for(Int_t ii=0; ii<20; ii++){
    fEMBufferCycle[ii]=0;
  }
  for(Int_t ii=0; ii<20; ii++){
    fEMBufferSize[ii]=10;
  }
  for(Int_t ii=0; ii<50; ii++){
    for(Int_t jj=0; jj<20; jj++){
      fEMBufferEta[ii][jj]=0;
      fEMBufferPhi[ii][jj]=0;
    }
  }
  for(Int_t ii=0; ii<10; ii++){
    fkIntervalMinEta[ii]=0;
    fkIntervalMaxEta[ii]=0;
  }
  fkIntervalMinEta[0] = -1.4;
  fkIntervalMaxEta[0] = +1.4;
}

AliAnalysisTaskMCPredictions2pc::AliAnalysisTaskMCPredictions2pc(const char *name, Int_t lNSmallBinning, Int_t lNLargeBinning, Int_t lRebinFactor, Int_t lNEtaBins, Int_t lNPhiBins)
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
fkNPhiBins(lNPhiBins),
fkVerboseMode(kFALSE),
fkDoEventMixing(kFALSE),

// Histograms
fHistForwardMult(0),
fHistNchVsForwardMult(0),
fHistEventCounter(0),
fHistChargedEta(0),

// 2pc histograms
fHist3dTrigger(0x0),
fHist3dAssoPions(0x0),
fHist3dAssoK0Short(0x0),
fHist3dAssoLambda(0x0),
fHist3dAssoAntiLambda(0x0),
fHist3dAssoXiMinus(0x0),
fHist3dAssoXiPlus(0x0),
fHist3dAssoOmegaMinus(0x0),
fHist3dAssoOmegaPlus(0x0),

fHist4d2pcPions(0x0),
fHist4d2pcK0Short(0x0),
fHist4d2pcLambda(0x0),
fHist4d2pcAntiLambda(0x0),
fHist4d2pcXiMinus(0x0),
fHist4d2pcXiPlus(0x0),
fHist4d2pcOmegaMinus(0x0),
fHist4d2pcOmegaPlus(0x0),

fHist4d2pcMixedPions(0x0),
fHist4d2pcMixedK0Short(0x0),
fHist4d2pcMixedLambda(0x0),
fHist4d2pcMixedAntiLambda(0x0),
fHist4d2pcMixedXiMinus(0x0),
fHist4d2pcMixedXiPlus(0x0),
fHist4d2pcMixedOmegaMinus(0x0),
fHist4d2pcMixedOmegaPlus(0x0),
fNMultBins(0),
fNPtBins(0)
{
  for(Int_t ii=0; ii<100; ii++){
    fMultBinBounds[ii] = 0.0;
    fPtBinBounds[ii] = 0.0;
  }
  for(Int_t ii=0; ii<20; ii++){
    fEMBufferFull[ii]=0;
  }
  for(Int_t ii=0; ii<20; ii++){
    fEMBufferCycle[ii]=0;
  }
  for(Int_t ii=0; ii<20; ii++){
    fEMBufferSize[ii]=10;
  }
  for(Int_t ii=0; ii<50; ii++){
    for(Int_t jj=0; jj<20; jj++){
      fEMBufferEta[ii][jj]=0;
      fEMBufferPhi[ii][jj]=0;
    }
  }
  for(Int_t ii=0; ii<10; ii++){
    fkIntervalMinEta[ii]=0;
    fkIntervalMaxEta[ii]=0;
  }
  fkIntervalMinEta[0] = -1.4;
  fkIntervalMaxEta[0] = +1.4;

  DefineOutput(1, TList::Class()); // Event Counter Histo
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
  //Int_t lNPtBins = 250;
  //Double_t lMaxPt = 25.0;
  
  Int_t lNEtaBins = fkNEtaBins;
  Double_t lMaxAbsEta = fkMaxEta;
  
  //Settings for charged particle counters (integers!)
  Int_t lNNchBins = fSmallMultRange/fRebinFactor;
  Double_t lLowNchBound  = -0.5;
  Double_t lHighNchBound = -0.5 + ((double)(fSmallMultRange));
  
  Int_t lNNchBinsForward = fLargeMultRange/fRebinFactor;
  Double_t lLowNchBoundForward  = -0.5;
  Double_t lHighNchBoundForward = -0.5 + ((double)(fLargeMultRange));

  std::vector<std::vector<double>> expandedAxes;
  std::vector<double> edgesEta;
  std::vector<double> edgesDeltaEta;
  std::vector<double> edgesDeltaPhi;
  std::vector<double> edgesPtTrigger;
  std::vector<double> edgesMult;
  Double_t etaBoundary,phiBoundary = 0.0;
  for (int i=0; i<lNEtaBins+1; i++)
  {
    Double_t etaBinWidth = 2*lMaxAbsEta/lNEtaBins;
    etaBoundary = -1*lMaxAbsEta+static_cast<float>(i)*etaBinWidth;
    edgesEta.emplace_back(etaBoundary);
  }
  for (int i=0; i<lNEtaBins+1; i++)
  {
    Double_t etaBinWidth = 4*lMaxAbsEta/lNEtaBins;
    etaBoundary = -2*lMaxAbsEta+static_cast<float>(i)*etaBinWidth;
    edgesDeltaEta.emplace_back(etaBoundary);
  }
  for (int i=0; i<fkNPhiBins+1; i++)
  {
    Double_t phiBinWidth = 2.0*TMath::Pi()/fkNPhiBins;
    phiBoundary = -0.5*TMath::Pi()+static_cast<float>(i)*phiBinWidth;
    edgesDeltaPhi.emplace_back(phiBoundary);
  }
  for (int i=0; i<fNPtBins+1; i++)
  {
    edgesPtTrigger.emplace_back(fPtBinBounds[i]);
  }
  for (int i=0; i<fNMultBins+1; i++)
  {
    edgesMult.emplace_back(fMultBinBounds[i]);
  }
  const Int_t nBins[4] = {lNEtaBins, fkNPhiBins, fNPtBins, fNMultBins};
  expandedAxes.emplace_back(edgesDeltaEta);
  expandedAxes.emplace_back(edgesDeltaPhi);
  expandedAxes.emplace_back(edgesPtTrigger);
  expandedAxes.emplace_back(edgesMult);
  
  if(! fHistEventCounter ) {
    //Histogram Output: Event-by-Event
    fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",fNMultBins, fMultBinBounds);
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
    if(! fHist3dTrigger ) {
    //Histogram Output: Event-by-Event
    fHist3dTrigger = new TH3D( "fHist3dTrigger", " ",lNEtaBins,edgesEta.data(), fNPtBins, edgesPtTrigger.data(),fNMultBins,edgesMult.data());
    fListHist->Add(fHist3dTrigger);
  }
  //___________________________________________________
    if(! fHist3dAssoPions ) {
    //Histogram Output: Event-by-Event
    fHist3dAssoPions = new TH3D( "fHist3dAssoPions", " ",lNEtaBins,edgesEta.data(), fNPtBins, edgesPtTrigger.data(),fNMultBins,edgesMult.data());
    fListHist->Add(fHist3dAssoPions);
  }
  //___________________________________________________
      if(! fHist3dAssoK0Short ) {
    //Histogram Output: Event-by-Event
    fHist3dAssoK0Short = new TH3D( "fHist3dAssoK0Short", " ",lNEtaBins,edgesEta.data(), fNPtBins, edgesPtTrigger.data(),fNMultBins,edgesMult.data());
    fListHist->Add(fHist3dAssoK0Short);
  }
  //___________________________________________________
      if(! fHist3dAssoLambda ) {
    //Histogram Output: Event-by-Event
    fHist3dAssoLambda = new TH3D( "fHist3dAssoLambda", " ",lNEtaBins,edgesEta.data(), fNPtBins, edgesPtTrigger.data(),fNMultBins,edgesMult.data());
    fListHist->Add(fHist3dAssoLambda);
  }
  //___________________________________________________
      if(! fHist3dAssoAntiLambda ) {
    //Histogram Output: Event-by-Event
    fHist3dAssoAntiLambda = new TH3D( "fHist3dAssoAntiLambda", " ",lNEtaBins,edgesEta.data(), fNPtBins, edgesPtTrigger.data(),fNMultBins,edgesMult.data());
    fListHist->Add(fHist3dAssoAntiLambda);
  }
  //___________________________________________________
      if(! fHist3dAssoXiMinus ) {
    //Histogram Output: Event-by-Event
    fHist3dAssoXiMinus = new TH3D( "fHist3dAssoXiMinus", " ",lNEtaBins,edgesEta.data(), fNPtBins, edgesPtTrigger.data(),fNMultBins,edgesMult.data());
    fListHist->Add(fHist3dAssoXiMinus);
  }
  //___________________________________________________
      if(! fHist3dAssoXiPlus ) {
    //Histogram Output: Event-by-Event
    fHist3dAssoXiPlus = new TH3D( "fHist3dAssoXiPlus", " ",lNEtaBins,edgesEta.data(), fNPtBins, edgesPtTrigger.data(),fNMultBins,edgesMult.data());
    fListHist->Add(fHist3dAssoXiPlus);
  }
  //___________________________________________________
      if(! fHist3dAssoOmegaMinus ) {
    //Histogram Output: Event-by-Event
    fHist3dAssoOmegaMinus = new TH3D( "fHist3dAssoOmegaMinus", " ",lNEtaBins,edgesEta.data(), fNPtBins, edgesPtTrigger.data(),fNMultBins,edgesMult.data());
    fListHist->Add(fHist3dAssoOmegaMinus);
  }
  //___________________________________________________
      if(! fHist3dAssoOmegaPlus ) {
    //Histogram Output: Event-by-Event
    fHist3dAssoOmegaPlus = new TH3D( "fHist3dAssoOmegaPlus", " ",lNEtaBins,edgesEta.data(), fNPtBins, edgesPtTrigger.data(),fNMultBins,edgesMult.data());
    fListHist->Add(fHist3dAssoOmegaPlus);
  }
  //___________________________________________________

  if(! fHist4d2pcPions ) {
    fHist4d2pcPions = new THnF("fHist4d2pcPions","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcPions);
  }
  if(! fHist4d2pcK0Short ) {
    fHist4d2pcK0Short = new THnF("fHist4d2pcK0Short","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcK0Short);
  }
  if(! fHist4d2pcLambda ) {
    fHist4d2pcLambda = new THnF("fHist4d2pcLambda","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcLambda);
  }
  if(! fHist4d2pcAntiLambda ) {
    fHist4d2pcAntiLambda = new THnF("fHist4d2pcAntiLambda","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcAntiLambda);
  }
  if(! fHist4d2pcXiMinus ) {
    fHist4d2pcXiMinus = new THnF("fHist4d2pcXiMinus","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcXiMinus);
  }
  if(! fHist4d2pcXiPlus ) {
    fHist4d2pcXiPlus = new THnF("fHist4d2pcXiPlus","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcXiPlus);
  }
  if(! fHist4d2pcOmegaMinus ) {
    fHist4d2pcOmegaMinus = new THnF("fHist4d2pcOmegaMinus","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcOmegaMinus);
  }
  if(! fHist4d2pcOmegaPlus ) {
    fHist4d2pcOmegaPlus = new THnF("fHist4d2pcOmegaPlus","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcOmegaPlus);
  }
  //___________________________________________________
  if(! fHist4d2pcMixedPions && fkDoEventMixing ) {
    fHist4d2pcMixedPions = new THnF("fHist4d2pcMixedPions","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcMixedPions);
  }
  if(! fHist4d2pcMixedK0Short && fkDoEventMixing) {
    fHist4d2pcMixedK0Short = new THnF("fHist4d2pcMixedK0Short","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcMixedK0Short);
  }
  if(! fHist4d2pcMixedLambda && fkDoEventMixing ) {
    fHist4d2pcMixedLambda = new THnF("fHist4d2pcMixedLambda","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcMixedLambda);
  }
  if(! fHist4d2pcMixedAntiLambda && fkDoEventMixing ) {
    fHist4d2pcMixedAntiLambda = new THnF("fHist4d2pcMixedAntiLambda","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcMixedAntiLambda);
  }
  if(! fHist4d2pcMixedXiMinus && fkDoEventMixing ) {
    fHist4d2pcMixedXiMinus = new THnF("fHist4d2pcMixedXiMinus","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcMixedXiMinus);
  }
  if(! fHist4d2pcMixedXiPlus && fkDoEventMixing ) {
    fHist4d2pcMixedXiPlus = new THnF("fHist4d2pcMixedXiPlus","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcMixedXiPlus);
  }
  if(! fHist4d2pcMixedOmegaMinus && fkDoEventMixing ) {
    fHist4d2pcMixedOmegaMinus = new THnF("fHist4d2pcMixedOmegaMinus","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcMixedOmegaMinus);
  }
  if(! fHist4d2pcMixedOmegaPlus && fkDoEventMixing ) {
    fHist4d2pcMixedOmegaPlus = new THnF("fHist4d2pcMixedOmegaPlus","",4, nBins, expandedAxes);
    fListHist->Add(fHist4d2pcMixedOmegaPlus);
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
  Long_t Nall = 0; 
  Long_t Ncharged = 0; 
  Long_t Nprimary = 0;

  Bool_t lEvSel_INELgtZEROStackPrimaries=kFALSE;
  
  
  //----- Loop on Stack ----------------------------------------------------------------
  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCevent->GetNumberOfTracks()); iCurrentLabelStack++)
  {   // This is the begining of the loop on tracks
    AliMCParticle* particleOne = (AliMCParticle*)lMCevent->GetTrack(iCurrentLabelStack);
    if(!particleOne) continue;
    Nall++; 
    Double_t lThisCharge = particleOne->Charge();
    if(TMath::Abs(lThisCharge)<0.001) continue;
    Ncharged++; 
    if(! (particleOne->IsPhysicalPrimary()) ) continue;
    Nprimary++;
    
    Double_t gpt = particleOne -> Pt();
    Double_t geta = particleOne -> Eta();
    
    if( TMath::Abs(geta) < 0.5 ) lNchEta5++;
    if( TMath::Abs(geta) < 0.8 ) lNchEta8++;
    if( TMath::Abs(geta) < 1.0 ) lEvSel_INELgtZEROStackPrimaries = kTRUE;
    
    //Special treatment: multiple intervals
    for(Int_t ii = 0; ii<fkNIntervals; ii++ )
      if( fkIntervalMinEta[ii] < geta && geta < fkIntervalMaxEta[ii] ) lNchForward++;
  }//End of loop on tracks
  //----- End Loop on Stack ------------------------------------------------------------

  if( fkVerboseMode ){ 
    Printf("Particle counters - all = %li, charged = %li, primary = %li", Nall, Ncharged, Nprimary); 
    Printf("Particle counters - |eta|<0.5: %li; |eta|<0.8: %li; custom counter: %f", lNchEta5, lNchEta8, lNchForward); 
  }

  //Reject non-INEL>0 if requested
  if( !lEvSel_INELgtZEROStackPrimaries && fkSelectINELgtZERO ) return;
  
  if( lNchForward < fkMinimumMultiplicity+1e-10 ) return;
  Int_t multiplicityIndex = fHistEventCounter->FindBin(lNchForward)-1; //use variable binning
  
  //------------------------------------------------
  // Fill Event Counters
  //------------------------------------------------
  
  //Basics: All Processed
  if( !fHistEventCounter ) {
    Printf("ERROR: Could not retrieve fHistEventCounter! This will crash!\n");
  }
  fHistEventCounter->Fill(lNchForward);
  
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

  std::vector<uint32_t> piIndices;
  std::vector<uint32_t> k0ShortIndices;
  std::vector<uint32_t> lambdaIndices;
  std::vector<uint32_t> antiLambdaIndices;
  std::vector<uint32_t> xiMinusIndices;
  std::vector<uint32_t> xiPlusIndices;
  std::vector<uint32_t> omegaMinusIndices;
  std::vector<uint32_t> omegaPlusIndices;

  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCevent->GetNumberOfTracks()); iCurrentLabelStack++)
  {
    // Determine if within acceptance, otherwise fully reject from list
    // done such that this check is done O(N) and not O(N^2)
    // TParticle* lThisParticle = lMCstack->Particle(iCurrentLabelStack);
    AliMCParticle* lThisParticle = (AliMCParticle*)lMCevent->GetTrack(iCurrentLabelStack);
    if(!lThisParticle) continue;
    Bool_t lIsPhysicalPrimary = lThisParticle->IsPhysicalPrimary();
    Double_t geta = lThisParticle -> Eta();
    Double_t gpt = lThisParticle -> Pt();

    // kick out stuff not at midrapidity
    if( ( geta < fkMinEta || geta > fkMaxEta) ) continue;
    if( !( lIsPhysicalPrimary ) ) continue;
    
    Double_t lThisCharge = lThisParticle->Charge();

    // populate triggers with charged particles within desired pT window
    if (fkMinPtTrigger<gpt && gpt<fkMaxPtTrigger && TMath::Abs(lThisCharge)>0.001){
      //valid trigger -- select also on PDG species 

      if( TMath::Abs(lThisParticle->PdgCode()) ==  211 ||
          TMath::Abs(lThisParticle->PdgCode()) ==  321 ||
          TMath::Abs(lThisParticle->PdgCode()) == 2212 || 
          TMath::Abs(lThisParticle->PdgCode()) ==   11 ||
          TMath::Abs(lThisParticle->PdgCode()) ==   13 )
      {
        fHist3dTrigger->Fill(geta, gpt,lNchForward);
        triggerIndices.emplace_back(iCurrentLabelStack);
      }
    }

    if ( TMath::Abs(lThisParticle->PdgCode()) ==  211 ){
      piIndices.emplace_back(iCurrentLabelStack);
      fHist3dAssoPions->Fill(geta, gpt,lNchForward);
    }
    if ( lThisParticle->PdgCode() ==  310 ){
      k0ShortIndices.emplace_back(iCurrentLabelStack);
      fHist3dAssoK0Short->Fill(geta, gpt,lNchForward);
    }
    if ( lThisParticle->PdgCode() ==  3122 ){
      lambdaIndices.emplace_back(iCurrentLabelStack);
      fHist3dAssoLambda->Fill(geta, gpt,lNchForward);
    }
    if ( lThisParticle->PdgCode() == -3122 ){
      antiLambdaIndices.emplace_back(iCurrentLabelStack);
      fHist3dAssoAntiLambda->Fill(geta, gpt,lNchForward);
    }
    if ( lThisParticle->PdgCode() ==  3312 ){
      xiMinusIndices.emplace_back(iCurrentLabelStack);
      fHist3dAssoXiMinus->Fill(geta, gpt,lNchForward);
    }
    if ( lThisParticle->PdgCode() == -3312 ){
      xiPlusIndices.emplace_back(iCurrentLabelStack);
      fHist3dAssoXiPlus->Fill(geta, gpt,lNchForward);
    }
    if ( lThisParticle->PdgCode() ==  3334 ){
      omegaMinusIndices.emplace_back(iCurrentLabelStack);
      fHist3dAssoOmegaMinus->Fill(geta, gpt,lNchForward);
    }
    if ( lThisParticle->PdgCode() == -3334 ){
      omegaPlusIndices.emplace_back(iCurrentLabelStack);
      fHist3dAssoOmegaPlus->Fill(geta, gpt,lNchForward);
    }
  }

  associatedIndices.emplace_back(piIndices);
  associatedIndices.emplace_back(k0ShortIndices);
  associatedIndices.emplace_back(lambdaIndices);
  associatedIndices.emplace_back(antiLambdaIndices);
  associatedIndices.emplace_back(xiMinusIndices);
  associatedIndices.emplace_back(xiPlusIndices);
  associatedIndices.emplace_back(omegaMinusIndices);
  associatedIndices.emplace_back(omegaPlusIndices);

  THnF *h4dSame[8] = {fHist4d2pcPions, fHist4d2pcK0Short, fHist4d2pcLambda, fHist4d2pcAntiLambda, fHist4d2pcXiMinus, fHist4d2pcXiPlus, fHist4d2pcOmegaMinus, fHist4d2pcOmegaPlus};
  THnF *h4dMixed[8] = {fHist4d2pcMixedPions, fHist4d2pcMixedK0Short, fHist4d2pcMixedLambda, fHist4d2pcMixedAntiLambda, fHist4d2pcMixedXiMinus, fHist4d2pcMixedXiPlus, fHist4d2pcMixedOmegaMinus, fHist4d2pcMixedOmegaPlus};

  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  //  Actually correlate stuff with stuff
  for (Int_t iTrigger = 0;  iTrigger < triggerIndices.size(); iTrigger++){   // trigger loop
    //TParticle* lTriggerParticle = lMCstack->Particle(triggerIndices[iTrigger]);
    AliMCParticle* lTriggerParticle = (AliMCParticle*)lMCevent->GetTrack(triggerIndices[iTrigger]);
    
    Double_t geta = lTriggerParticle -> Eta();
    Double_t gphi = lTriggerParticle -> Phi();

    for (Int_t iassocSpecies = 0;  iassocSpecies < associatedIndices.size(); iassocSpecies++){   // associated loop
      for (Int_t iassoc = 0;  iassoc < associatedIndices[iassocSpecies].size(); iassoc++){   // associated loop
        if( triggerIndices[iTrigger] == associatedIndices[iassocSpecies][iassoc] ) continue; // avoid self
        AliMCParticle* lAssociatedParticle = (AliMCParticle*)lMCevent->GetTrack( associatedIndices[iassocSpecies][iassoc] );
        if(!lAssociatedParticle) {
          Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
          continue;
        }
        Double_t geta2 = lAssociatedParticle -> Eta();
        Double_t gphi2 = lAssociatedParticle -> Phi();
        Double_t lThisPt    = lAssociatedParticle->Pt();
        h4dSame[iassocSpecies]->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt,lNchForward);
      }
    }
  }
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  //Event mixing loop
  if( fkDoEventMixing ){ //require also that a trigger exists
    for (Int_t iTrigger = 0;  iTrigger < TMath::Min(fEMBufferSize[multiplicityIndex], fEMBufferFull[multiplicityIndex]); iTrigger++){   // trigger loop
      Double_t geta = fEMBufferEta[iTrigger][multiplicityIndex]; //from previous events
      Double_t gphi = fEMBufferPhi[iTrigger][multiplicityIndex]; //from previous events

      for (Int_t iassocSpecies = 0;  iassocSpecies < associatedIndices.size(); iassocSpecies++){   // associated loop
        for (Int_t iassoc = 0;  iassoc < associatedIndices[iassocSpecies].size(); iassoc++){   // associated loop
          //lAssociatedParticle = lMCstack->Particle( associatedIndices[iassocSpecies][iassoc] );
          AliMCParticle* lAssociatedParticle = (AliMCParticle*)lMCevent->GetTrack(associatedIndices[iassocSpecies][iassoc] );
          if(!lAssociatedParticle) {
            Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iassoc );
            continue;
          }
          Double_t geta2 = lAssociatedParticle -> Eta();
          Double_t gphi2 = lAssociatedParticle -> Phi();
          Double_t lThisPt    = lAssociatedParticle->Pt();
          h4dMixed[iassocSpecies]->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt,lNchForward);
        }
      }
    }
  }
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // fill EM buffer
  if( fEMBufferSize[multiplicityIndex] > 0 && fkDoEventMixing){
    for (Int_t iTrigger = 0;  iTrigger < triggerIndices.size(); iTrigger++){   // trigger loop
      //TParticle* lThisParticle = lMCstack->Particle(triggerIndices[iTrigger]);
      AliMCParticle* lThisParticle = (AliMCParticle*)lMCevent->GetTrack(triggerIndices[iTrigger]);
      if(!lThisParticle) continue;

      //Add to buffer
      fEMBufferEta[fEMBufferCycle[multiplicityIndex]][multiplicityIndex] = lThisParticle->Eta();
      fEMBufferPhi[fEMBufferCycle[multiplicityIndex]][multiplicityIndex] = lThisParticle->Phi();
      fEMBufferCycle[multiplicityIndex]++;
      fEMBufferFull[multiplicityIndex]++;
      if(fEMBufferFull[multiplicityIndex]>fEMBufferSize[multiplicityIndex]) fEMBufferFull[multiplicityIndex] = fEMBufferSize[multiplicityIndex];
      fEMBufferCycle[multiplicityIndex] = fEMBufferCycle[multiplicityIndex]%fEMBufferSize[multiplicityIndex];
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
