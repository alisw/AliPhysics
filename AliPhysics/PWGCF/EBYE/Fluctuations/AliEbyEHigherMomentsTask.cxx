/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Satyajit Jena.                                                 *
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


//=========================================================================//
//                                                                         //
//           Analysis Task for Net-Charge Higher Moment Analysis           //
//              Author: Satyajit Jena || Nirbhay K. Behera                 //
//                      sjena@cern.ch || nbehera@cern.ch                   //
//           Last Modified: 30/01/2014: only for net-charge part           //
//=========================================================================//
#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TCanvas.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliStack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliPID.h"
#include "AliAODPid.h"                
#include "AliPIDResponse.h"
#include "AliAODpidUtil.h"        
#include "AliPIDCombined.h"
#include "AliHelperPID.h"

#include "AliEbyEHigherMomentsTask.h"

using std::cout;
using std::endl;

ClassImp(AliEbyEHigherMomentsTask)


//-----------------------------------------------------------------------
AliEbyEHigherMomentsTask::AliEbyEHigherMomentsTask( const char *name )
: AliAnalysisTaskSE( name ),
  fListOfHistos(0),
  fArrayMC(0),
  fAnalysisType(0),
  fCentralityEstimator("V0M"),
  fCentrality(0),
  fVxMax(3.),
  fVyMax(3.),
  fVzMax(10.),
  fPtLowerLimit(0.3),
  fPtHigherLimit(1.5),
  fNptBins(0),
  fBin(0),
  fEtaLowerLimit(-0.8),
  fEtaHigherLimit(0.8),
  fAODtrackCutBit(128),
  fEventCounter(0),
  fHistVxVy(0),
  fTHnCentNplusNminusCh(0),
  fTHnCentNplusNminusChTruth(0),
  fPtBinNplusNminusCh(0),
  fPtBinNplusNminusChTruth(0)
{ 
  
  for ( Int_t i = 0; i < 4; i++) { 
    fHistQA[i] = NULL;
  }
 
  DefineOutput(1, TList::Class()); // Outpput....
  //DefineOutput(2, TList::Class()); 
}

AliEbyEHigherMomentsTask::~AliEbyEHigherMomentsTask() {

  if(fListOfHistos) delete fListOfHistos;
}

//---------------------------------------------------------------------------------
void AliEbyEHigherMomentsTask::UserCreateOutputObjects() {
  
  fListOfHistos = new TList();
  fListOfHistos->SetOwner(kTRUE);
  
  fEventCounter = new TH1D("fEventCounter","EventCounter", 10, 0.5,10.5);
  fEventCounter->GetXaxis()->SetBinLabel(1,"Event Accesed");
  fEventCounter->GetXaxis()->SetBinLabel(2,"Within 0-80% centrality");
  fEventCounter->GetXaxis()->SetBinLabel(5,"Have a vertex");
  fEventCounter->GetXaxis()->SetBinLabel(6,"After vertex Cut");
  fEventCounter->GetXaxis()->SetBinLabel(7,"Event Analyzed");
  fListOfHistos->Add(fEventCounter);
  
  //For QA-Histograms
  fHistQA[0] = new TH1D("fHistQAvz", "Histo Vz After Cut", 500, -25., 25.);  
  fHistQA[1] = new TH1D("fHistQAPt","p_{T} distribution",600,0,6);
  fHistQA[2] = new TH1D("fHistQAEta","#eta distribution",240,-1.2,1.2);
  fHistQA[3] = new TH1D("fHistQAPhi","#phi distribution",340,0,6.8);
 
  for(Int_t i = 0; i < 4; i++)
    {
      fListOfHistos->Add(fHistQA[i]);
    }
  
  fHistVxVy = new TH2D("fHistVxVy","Vertex-x Vs Vertex-y", 200, -1., 1., 200, -1., 1.);
  fListOfHistos->Add(fHistVxVy);
 
 
  const Int_t nDim = 3;
  
  Int_t fBinsCh[nDim] = {100, 1500, 1500};
  Double_t fMinCh[nDim] = { -0.5, -0.5, -0.5 };
  Double_t fMaxCh[nDim] = { 99.5, 1499.5, 1499.5};
  
  const Int_t dim = fNptBins*2 + 1;
  //const Int_t dim = ;
  Int_t bin[dim];
  Double_t min[dim];
  Double_t max[dim];
  bin[0] = 100; min[0] = -0.5; max[0] = 99.5;

  for(Int_t i = 1; i < dim; i++){
    
    bin[i] = 900;
    min[i] = -0.5;
    max[i] = 899.5;
 
  } 
  
  fTHnCentNplusNminusCh = new THnSparseD("fTHnCentNplusNminusCh","Cent-NplusChrg-NminusChrg", nDim, fBinsCh, fMinCh, fMaxCh); 
  fTHnCentNplusNminusCh->GetAxis(0)->SetTitle("Centrality");
  fTHnCentNplusNminusCh->GetAxis(1)->SetTitle("Nplus");
  fTHnCentNplusNminusCh->GetAxis(2)->SetTitle("Nminus");
  fListOfHistos->Add(fTHnCentNplusNminusCh);

  fPtBinNplusNminusCh = new THnSparseI("fPtBinNplusNminusCh","cent-nplus-nminus", dim, bin, min, max);
  
  fListOfHistos->Add(fPtBinNplusNminusCh);
  
  if( fAnalysisType == "MCAOD" ){   
    
    fTHnCentNplusNminusChTruth = new THnSparseD("fTHnCentNplusNminusChTruth","Cent-NplusChrg-NminusChrg", nDim, fBinsCh, fMinCh, fMaxCh); 
    fTHnCentNplusNminusChTruth->GetAxis(0)->SetTitle("Centrality");
    fTHnCentNplusNminusChTruth->GetAxis(1)->SetTitle("Nplus");
    fTHnCentNplusNminusChTruth->GetAxis(2)->SetTitle("Nminus");
    fListOfHistos->Add(fTHnCentNplusNminusChTruth);
    
    fPtBinNplusNminusChTruth = new THnSparseI("fPtBinNplusNminusChTruth","cent-nplus-nminus", dim, bin, min, max);
    fListOfHistos->Add(fPtBinNplusNminusChTruth);
    
  }//MCAOD---
  
 
  
  //PostData(1, fListOfHistosQA);
  PostData(1, fListOfHistos);   
  
  
}


//----------------------------------------------------------------------------------
void AliEbyEHigherMomentsTask::UserExec( Option_t * ){
  
  fEventCounter->Fill(1);
  
  if(fAnalysisType == "AOD") {
    
    doAODEvent();
    
  }//AOD--analysis-----
  
  else if(fAnalysisType == "MCAOD") {
    
    doMCAODEvent();
    
  }
  
  else return;
 
  
  
}

//--------------------------------------------------------------------------------------
void AliEbyEHigherMomentsTask::doAODEvent(){
  
  Double_t positiveSum = 0.;
  Double_t negativeSum = 0.;
  const Int_t dim = fNptBins*2;
  Double_t PtCh[dim];

  for(Int_t idx = 0; idx < dim; idx++){
    PtCh[idx] = 0.;
  }


  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    cout<<"ERROR: Analysis manager not found."<<endl;
    return;
  }
  //coneect to the inputHandler------------
  AliAODInputHandler* inputHandler = dynamic_cast<AliAODInputHandler*> (manager->GetInputEventHandler());
  if (!inputHandler) {
    cout<<"ERROR: Input handler not found."<<endl;
    return;
  }
  
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());

  if (!fAOD)
    {
      cout<< "ERROR: AOD not found " <<endl;
      return;
    }
  
  if(!ProperVertex(fAOD)) return;   
  
  AliAODHeader *aodHeader = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
  if(!aodHeader) AliFatal("Not a standard AOD");
  
  fCentrality = (Int_t)aodHeader->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
  
  if(fCentrality < 0 || fCentrality >= 81) return;
  
  fEventCounter->Fill(2);
  
  Int_t nTracks = fAOD->GetNumberOfTracks();
  
  for(Int_t i = 0; i < nTracks; i++) {
    
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(i)); 
    
    if(!aodTrack) {
      AliError(Form("ERROR: Could not retrieve AODtrack %d",i));
      continue;
    }
    
    if(!AcceptTrack(aodTrack) ) continue;

    
    fBin = GetPtBin(aodTrack->Pt());
    
    Short_t gCharge = aodTrack->Charge();
    if(gCharge > 0){
      positiveSum += 1.;
      PtCh[fBin] += 1;
    }
    if(gCharge < 0){
      negativeSum += 1.;
      PtCh[fNptBins+fBin] += 1;
    }
   
   
  }//--------- Track Loop to select with filterbit
  
  Double_t fContainerCh[3] = { static_cast<Double_t>(fCentrality), positiveSum, negativeSum};
  
  fTHnCentNplusNminusCh->Fill(fContainerCh);

  //cout << fCentrality << " " << positiveSum <<" " << negativeSum << endl;
  
  Double_t ptContainer[dim+1];

  ptContainer[0] = fCentrality;

  for(Int_t i = 1; i <= dim; i++){
    ptContainer[i] = PtCh[i-1];
    //cout << PtCh[i-1] <<" ,";
  }
  //cout << endl;

  fPtBinNplusNminusCh->Fill(ptContainer);
  
  
  fEventCounter->Fill(7);
  PostData(1, fListOfHistos);
  return;
  
}
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
void AliEbyEHigherMomentsTask::doMCAODEvent(){
  
  
  //---------
  Double_t positiveSumMCRec = 0.;
  Double_t negativeSumMCRec = 0.;
 
  Double_t positiveSumMCTruth = 0.;
  Double_t negativeSumMCTruth = 0.;
 
  const Int_t dim = fNptBins*2;
  Double_t PtCh[dim];
  Double_t ptChMC[dim];
  for(Int_t idx = 0; idx < dim; idx++){
    PtCh[idx] = 0;
    ptChMC[idx] = 0;
  }
  
  //Connect to Anlaysis manager------
  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    cout<<"ERROR: Analysis manager not found."<<endl;
    return;
  }
  
  AliAODInputHandler* inputHandler = dynamic_cast<AliAODInputHandler*> (manager->GetInputEventHandler());
  if (!inputHandler) {
    cout<<"ERROR: Input handler not found."<<endl;
    return;
  }
  
  //AOD event------
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  
  if (!fAOD)
    {
      cout<< "ERROR: AOD not found " <<endl;
      return;
    }

  // -- Get MC header
  // ------------------------------------------------------------------
  
  fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!fArrayMC) {
    AliFatal("Error: MC particles branch not found!\n");
    return;
  }
  
  AliAODMCHeader *mcHdr=NULL;
  mcHdr=(AliAODMCHeader*)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());  
  if(!mcHdr) {
    printf("MC header branch not found!\n");
    return;
  }
  
  if(!ProperVertex(fAOD)) return;

  AliAODHeader *aodHeader = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
  if(!aodHeader) AliFatal("Not a standard AOD");
  
  fCentrality = (Int_t)aodHeader->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());

  if( fCentrality < 0 || fCentrality >= 81) return;
  
  fEventCounter->Fill(2);
 
  Int_t nTracks = fAOD->GetNumberOfTracks();
  
  for(Int_t i = 0; i < nTracks; i++) {
    
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(i)); 
    
    if(!aodTrack) {
      AliError(Form("ERROR: Could not retrieve AODtrack %d",i));
      continue;
    }
    
    if(!AcceptTrack(aodTrack) ) continue;

    
    fBin = GetPtBin(aodTrack->Pt());
    
    //cout << "Pt Bin " << fBin << endl;

    Short_t gCharge = aodTrack->Charge();
    if(gCharge > 0){
      positiveSumMCRec += 1.;
     PtCh[fBin] += 1;
    }
    if(gCharge < 0){
      negativeSumMCRec += 1.;    
      PtCh[fNptBins+fBin] += 1;
    }
   
  }//--------- Track Loop to select with filterbit
  
  //===============================================================================
  //---------------------MC Truth--------------------------------------------------
  //===============================================================================
  
  Int_t nMCTrack = fArrayMC->GetEntriesFast();
  
  for (Int_t iMC = 0; iMC < nMCTrack; iMC++) {
    
    AliAODMCParticle *partMC = (AliAODMCParticle*) fArrayMC->At(iMC);
    
    if(!partMC){
      AliError(Form("ERROR: Could not retrieve AODMCtrack %d",iMC));
      continue;
    }
    
    if(partMC->Charge() == 0) continue;
    if(!partMC->IsPhysicalPrimary()) continue;
    
    if (partMC->Eta() < fEtaLowerLimit || partMC->Eta() > fEtaHigherLimit) continue;
    if (partMC->Pt() < fPtLowerLimit ||  partMC->Pt() > fPtHigherLimit) continue;
    
    Short_t gCharge = partMC->Charge();
    Int_t mcbin = GetPtBin(partMC->Pt());
    
    if(gCharge > 0){
      positiveSumMCTruth += 1.;
      ptChMC[mcbin] += 1.;
    }
    if(gCharge < 0){
      negativeSumMCTruth += 1.;
      ptChMC[fNptBins+mcbin] += 1.;
    }
   
   
  }//MC-Truth Track loop--

  Double_t fContainerCh[3] = { static_cast<Double_t>(fCentrality), positiveSumMCRec, negativeSumMCRec};//Reco. values ch. hadrons
  Double_t fContainerChTruth[3] = { static_cast<Double_t>(fCentrality), positiveSumMCTruth, negativeSumMCTruth};  
 
  fTHnCentNplusNminusCh->Fill(fContainerCh);//Fill the rec. ch. particles---
  fTHnCentNplusNminusChTruth->Fill(fContainerChTruth);//MC -Truth ch. particles
  
  //cout << fCentrality << " " << positiveSumMCRec << " " << negativeSumMCRec <<endl;
  //cout <<fCentrality<<"   " << positiveSumMCTruth << " " << negativeSumMCTruth << endl;
  //cout <<"   " << posPidSumMCRec << " " << negPidSumMCRec << endl;
  //cout <<"   " << posPidSumMCTruth << " " << negPidSumMCTruth << endl;
  //cout <<"---------------------------------" << endl;
  
  Double_t ptContainer[dim+1];
  Double_t ptContainerMC[dim+1];
  ptContainer[0] = fCentrality;
  ptContainerMC[0] = fCentrality;

  for(Int_t i = 1; i <= dim; i++){
    ptContainer[i] = PtCh[i-1];
    ptContainerMC[i] = ptChMC[i-1];
    //cout <<" "<< PtCh[i-1]<<endl;
    //  cout<< " MC=" << ptChMC[i-1];
  }

  //cout << endl;

  fPtBinNplusNminusCh->Fill(ptContainer);
  fPtBinNplusNminusChTruth->Fill(ptContainerMC);

  fEventCounter->Fill(7);
  PostData(1, fListOfHistos);
  return;
  
}

//---------------------------------------------------------------------------------------
Bool_t AliEbyEHigherMomentsTask::ProperVertex(AliAODEvent *fAOD) const{
  
  Bool_t IsVtx = kFALSE;
  
  const AliAODVertex *vertex = fAOD->GetPrimaryVertex();
  
  if(vertex) {
    Double32_t fCov[6];
    vertex->GetCovarianceMatrix(fCov);
    if(vertex->GetNContributors() > 0) {
      if(fCov[5] != 0) {
	
	Double_t lvx = vertex->GetX();
	Double_t lvy = vertex->GetY();
	Double_t lvz = vertex->GetZ();
	
	fEventCounter->Fill(5);
      
	if(TMath::Abs(lvx) < fVxMax) {
	  if(TMath::Abs(lvy) < fVyMax) {
	    if(TMath::Abs(lvz) < fVzMax) {
	      
	      fEventCounter->Fill(6);
	      fHistQA[0]->Fill(lvz);
	      fHistVxVy->Fill(lvx,lvy);
	      IsVtx = kTRUE; 
	      
	    }//Z-Vertex cut---
 	  }//Y-vertex cut--
	}//X-vertex cut---
      }//Covariance------
    }//Contributors check---
  }//If vertex-----------

  AliCentrality *centrality = fAOD->GetCentrality();
  if (centrality->GetQuality() != 0) IsVtx = kFALSE;
  
  return IsVtx;
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
Bool_t AliEbyEHigherMomentsTask::AcceptTrack(AliAODTrack* track) const{
  
  if(!track) return kFALSE;
  if( track->Charge() == 0 ) return kFALSE;

  Double_t pt = track->Pt();
  Double_t eta = track->Eta();
  if(!track->TestFilterBit(fAODtrackCutBit)) return kFALSE;
  if( pt < fPtLowerLimit || pt > fPtHigherLimit ) return kFALSE;
  if( eta < fEtaLowerLimit || eta > fEtaHigherLimit ) return kFALSE;
  
  fHistQA[1]->Fill(pt);
  fHistQA[2]->Fill(eta);
  fHistQA[3]->Fill(track->Phi());
 
  return kTRUE;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
Int_t AliEbyEHigherMomentsTask::GetPtBin(Double_t pt){
  
  Int_t bin = 0;

  Double_t BinSize = (fPtHigherLimit - fPtLowerLimit)/fNptBins;
  
  for(Int_t iBin = 0; iBin < fNptBins; iBin++){
    
    Double_t xLow = fPtLowerLimit + iBin*BinSize;
    Double_t xHigh = fPtLowerLimit + (iBin+1)*BinSize;
    
    if( pt >= xLow && pt < xHigh){
      bin = iBin;
      //cout << "Pt Bin #" << bin <<"pt=" << pt << endl;
      break;
    }
    
  }//for 
  
  
  return bin;
}
//------------------------------------------------------------------------
//------------------------------------------------------------------------
void AliEbyEHigherMomentsTask::Terminate( Option_t * ){
  
  Info("AliEbyEHigherMomentTask"," Task Successfully finished");
  
}

//------------------------------------------------------------------------
