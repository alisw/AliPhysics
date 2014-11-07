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
//                                                                         //
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

#include "AliEbyEHigherMomentsTaskPID.h"

using std::cout;
using std::endl;

ClassImp(AliEbyEHigherMomentsTaskPID)


//-----------------------------------------------------------------------
AliEbyEHigherMomentsTaskPID::AliEbyEHigherMomentsTaskPID( const char *name )
: AliAnalysisTaskSE( name ),
  fListOfHistos(0),
  fArrayMC(0),
  fPIDResponse(0),
  fParticleSpecies(AliPID::kProton),
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
  fRapidityCut(0.5),
  fNSigmaCut(3.),
  fAODtrackCutBit(128),
  fHelperPID(0),
  fEventCounter(0),
  fHistVxVy(0),
  fTHnCentNplusNminusPid(0),
  fTHnCentNplusNminusPidTruth(0),
  fPtBinNplusNminusPid(0),
  fPtBinNplusNminusPidTruth(0)
{ 
  
  for ( Int_t i = 0; i < 4; i++) { 
    fHistQA[i] = NULL;
  }
  
  DefineOutput(1, TList::Class()); // Outpput...
}

AliEbyEHigherMomentsTaskPID::~AliEbyEHigherMomentsTaskPID(){
  
  if(fListOfHistos) delete fListOfHistos;
  if(fHelperPID) delete fHelperPID;
}

//---------------------------------------------------------------------------------
void AliEbyEHigherMomentsTaskPID::UserCreateOutputObjects() {

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
  const Int_t nPid = 5;
  TString Species[nPid] = {"Electron","Muon","Pion","Kaon","Proton"};
  
  Int_t fBins[nPid][nDim];
  Double_t fMin[nPid][nDim];
  Double_t fMax[nPid][nDim];
  
  for( Int_t i = 0; i < nPid; i++ ){
    for( Int_t j = 0; j < nDim; j++ ){
      fBins[i][j] = 0;
      fMin[i][j] = 0.;
      fMax[i][j] = 0.;
    }
  }  
  

  const Int_t dim = fNptBins*2 + 1;
  Int_t bin[nPid][dim];
  Double_t min[nPid][dim];
  Double_t max[nPid][dim];

  bin[2][0] = 100; min[2][0] = -0.5; max[2][0] = 99.5;
  bin[3][0] = 100; min[3][0] = -0.5; max[3][0] = 99.5;
  bin[4][0] = 100; min[4][0] = -0.5; max[4][0] = 99.5;

  for(Int_t i = 1; i < dim; i++){
    
    bin[2][i] = 900;
    min[2][i] = -0.5;
    max[2][i] = 899.5;
    
    bin[3][i] = 700;
    min[3][i] = -0.5;
    max[3][i] = 699.5;

    bin[4][i] = 400;
    min[4][i] = -0.5;
    max[4][i] = 399.5;

 
  } 


  
  
 
  TString hname1, hname11;
  TString htitle1, axisTitle1, axisTitle2; 
  
  Int_t gPid = (Int_t)fParticleSpecies;
  
  //Pion----
  fBins[2][0] = 100; fBins[2][1] = 1000; fBins[2][2] = 1000;
  fMin[2][0] = -0.5; fMin[2][1] = -0.5; fMin[2][2] = -0.5;
  fMax[2][0] = 99.5; fMax[2][1] = 999.5; fMax[2][2] = 999.5;
  //Kaon------
  fBins[3][0] = 100; fBins[3][1] = 700; fBins[3][2] = 700;
  fMin[3][0] = -0.5; fMin[3][1] = -0.5; fMin[3][2] = -0.5;
  fMax[3][0] = 99.5; fMax[3][1] = 699.5; fMax[3][2] = 699.5;
  //Proton-----
  fBins[4][0] = 100; fBins[4][1] = 400; fBins[4][2] = 400;
  fMin[4][0] = -0.5; fMin[4][1] = -0.5; fMin[4][2] = -0.5;
  fMax[4][0] = 99.5; fMax[4][1] = 399.5; fMax[4][2] = 399.5;
  
  hname1 = "fCentNplusNminusPid"; hname1 +=gPid;
  htitle1 = Species[gPid]; htitle1 +=" And Neg-"; htitle1 +=Species[gPid];
  axisTitle1 ="Pos-"; axisTitle1 += Species[gPid];
  axisTitle2 ="Neg-"; axisTitle2 += Species[gPid];
  
  fTHnCentNplusNminusPid = new THnSparseD(hname1.Data(),htitle1.Data(), nDim, fBins[gPid], fMin[gPid],fMax[gPid]);
  
  fTHnCentNplusNminusPid->GetAxis(0)->SetTitle("Centrality");
  fTHnCentNplusNminusPid->GetAxis(1)->SetTitle(axisTitle1.Data());
  fTHnCentNplusNminusPid->GetAxis(2)->SetTitle(axisTitle2.Data());
  fListOfHistos->Add(fTHnCentNplusNminusPid);

  TString hname2 = "fPtBinNplusNminusPid"; 
  TString htitle2 = "cent-nplus-nminus-ptbinwise";
  hname2 += gPid;
  htitle2 += gPid;
  fPtBinNplusNminusPid = new THnSparseI(hname2.Data(),htitle2.Data(), dim, bin[gPid], min[gPid], max[gPid]);  
  fListOfHistos->Add(fPtBinNplusNminusPid);
  
  
  
  
  if( fAnalysisType == "MCAOD" ){
   
    hname11 = "fCentNplusNminusPidTruth"; hname11 +=gPid;
    fTHnCentNplusNminusPidTruth = new THnSparseD(hname11.Data(),htitle1.Data(), nDim, fBins[gPid], fMin[gPid],fMax[gPid]);
    
    fTHnCentNplusNminusPidTruth->GetAxis(0)->SetTitle("Centrality");
    fTHnCentNplusNminusPidTruth->GetAxis(1)->SetTitle(axisTitle1.Data());
    fTHnCentNplusNminusPidTruth->GetAxis(2)->SetTitle(axisTitle2.Data());
    fListOfHistos->Add(fTHnCentNplusNminusPidTruth);
    
    TString hname22 = "fPtBinNplusNminusPidTruth"; 
    TString htitle22 = "cent-nplus-nminus-ptbinwise";
    hname22 += gPid;
    htitle22 += gPid;    
    
    fPtBinNplusNminusPidTruth = new THnSparseI(hname22.Data(),htitle22.Data(), dim, bin[gPid], min[gPid], max[gPid]);
    fListOfHistos->Add(fPtBinNplusNminusPidTruth);
    
  }//fUsePid-------

  PostData(1, fListOfHistos);   
  
  
}


//----------------------------------------------------------------------------------
void AliEbyEHigherMomentsTaskPID::UserExec( Option_t * ){
  
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
void AliEbyEHigherMomentsTaskPID::doAODEvent(){
 
  Double_t posPidSum = 0.;
  Double_t negPidSum = 0.;

  const Int_t dim = fNptBins*2;
  Double_t PtCh[dim];
  
  for(Int_t idx = 0; idx < dim; idx++){
    PtCh[idx] = 0.;
  }
  
 
  Int_t gPid = (Int_t)fParticleSpecies;  

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
  
  fPIDResponse = inputHandler->GetPIDResponse();
  
  
  if (!fPIDResponse){
    AliFatal("This Task needs the PID response attached to the inputHandler");
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
    
    
    Double_t rap =  aodTrack->Y(AliPID::ParticleMass(fParticleSpecies));
    if( fabs(rap ) > fRapidityCut ) continue;//Rapidity cut
    
    Int_t PID = fHelperPID->GetParticleSpecies((AliVTrack*)aodTrack, kFALSE);     
    
    if( (PID+2) == gPid ){
      
      if(gCharge > 0){
	posPidSum += 1.;
	PtCh[fBin] += 1;
      }
      if(gCharge < 0){
	negPidSum += 1.;
	PtCh[fNptBins+fBin] += 1;
      }
      
    }//chek the PID
    
    
  }//--------- Track Loop to select with filterbit
  
  
  Double_t fContainerPid[3] = { static_cast<Double_t>(fCentrality), posPidSum, negPidSum};
  fTHnCentNplusNminusPid->Fill(fContainerPid);
  
  //cout << fCentrality <<" "<< posPidSum <<" " << negPidSum << endl;

  Double_t ptContainer[dim+1];
  
  ptContainer[0] = fCentrality;
  
  for(Int_t i = 1; i <= dim; i++){
    ptContainer[i] = PtCh[i-1];
    //cout << PtCh[i-1] <<" ,";
  }
  //cout << endl;

  fPtBinNplusNminusPid->Fill(ptContainer);
  
  
  fEventCounter->Fill(7);
  PostData(1, fListOfHistos);
  return;
  
}
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
void AliEbyEHigherMomentsTaskPID::doMCAODEvent(){
  
  
  //---------
  Double_t posPidSumMCRec = 0.;
  Double_t negPidSumMCRec = 0.;

  Double_t posPidSumMCTruth = 0.;
  Double_t negPidSumMCTruth = 0.;
  
  const Int_t dim = fNptBins*2;
  Double_t PtCh[dim];
  Double_t ptChMC[dim];

  for(Int_t idx = 0; idx < dim; idx++){
    PtCh[idx] = 0;
    ptChMC[idx] = 0;
  }
  
  
  Int_t gPid = (Int_t)fParticleSpecies; 
  
  Int_t gPdgCode = AliPID::ParticleCode(fParticleSpecies);
  
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
  
  fPIDResponse =inputHandler->GetPIDResponse();
  
  
  if(!fPIDResponse){
    AliFatal("This Task needs the PID response attached to the inputHandler");
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
    
    Double_t rap =  aodTrack->Y(AliPID::ParticleMass(fParticleSpecies));    
    if( fabs(rap ) > fRapidityCut ) continue;//Rapidity cut
    
    Int_t PID = fHelperPID->GetParticleSpecies((AliVTrack*)aodTrack, kFALSE);  
    
    if( (PID+2) == gPid  ){
      if(gCharge > 0){
	posPidSumMCRec += 1;
	PtCh[fBin] += 1;
      }
      if( gCharge < 0 ){
	negPidSumMCRec += 1.;
	PtCh[fNptBins+fBin] += 1;
      }
    }//nSigmaCut-----
  
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
  
    Double_t rap =  partMC->Y();
    if( fabs(rap ) > fRapidityCut ) continue;//Rapidity cut
    if(TMath::Abs(partMC->GetPdgCode()) != gPdgCode) continue;
    
    if(gCharge > 0){
      posPidSumMCTruth += 1.;
      ptChMC[mcbin] += 1;
    }
    if(gCharge < 0){
      negPidSumMCTruth += 1.;
      ptChMC[fNptBins+mcbin] += 1;
    }
    
    
  }//MC-Truth Track loop--
  
  
  //cout << fCentrality << " " << positiveSumMCRec << " " << negativeSumMCRec <<endl;
  //cout <<"   " << positiveSumMCTruth << " " << negativeSumMCTruth << endl;
  //cout <<fCentrality<<"   " << posPidSumMCRec << " " << negPidSumMCRec << endl;
  //cout <<fCentrality<<"   " << posPidSumMCTruth << " " << negPidSumMCTruth << endl;
  //cout <<"---------------------------------" << endl;
    
  Double_t fContainerPid[3] = { static_cast<Double_t>(fCentrality), posPidSumMCRec, negPidSumMCRec};//Reco.
  Double_t fContainerPidTruth[3] = { static_cast<Double_t>(fCentrality), posPidSumMCTruth, negPidSumMCTruth};
  
  fTHnCentNplusNminusPid->Fill(fContainerPid);//Fill the rec. pid tracks
  fTHnCentNplusNminusPidTruth->Fill(fContainerPidTruth);//MC-Truth pid
  
  Double_t ptContainer[dim+1];
  Double_t ptContainerMC[dim+1];
  ptContainer[0] = fCentrality;
  ptContainerMC[0] = fCentrality;
  
  for(Int_t i = 1; i <= dim; i++){
    ptContainer[i] = PtCh[i-1];
    ptContainerMC[i] = ptChMC[i-1];
    //cout <<" "<< PtCh[i-1]<<endl;
    //cout<< " Rec=" << PtCh[i-1];
  }

  //cout << endl;
  
  fPtBinNplusNminusPid->Fill(ptContainer);
  fPtBinNplusNminusPidTruth->Fill(ptContainerMC);

  fEventCounter->Fill(7);
  PostData(1, fListOfHistos);
  return;
  
}

//---------------------------------------------------------------------------------------
Bool_t AliEbyEHigherMomentsTaskPID::ProperVertex(AliAODEvent *fAOD) const{
  
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
Bool_t AliEbyEHigherMomentsTaskPID::AcceptTrack(AliAODTrack* track) const{
  
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
Int_t AliEbyEHigherMomentsTaskPID::GetPtBin(Double_t pt){
  
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
void AliEbyEHigherMomentsTaskPID::Terminate( Option_t * ){
  
  Info("AliEbyEHigherMomentTask"," Task Successfully finished");
  
}

//------------------------------------------------------------------------
