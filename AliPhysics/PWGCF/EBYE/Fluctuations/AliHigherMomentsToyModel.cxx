/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Nirbhay K. Behera                                              *
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
//           Toy Model for Net-Charge Higher Moment Analysis               //
//              Author: Satyajit Jena || Nirbhay K. Behera                 //
//                      sjena@cern.ch || nbehera@cern.ch                   //
//                                                                         //
//=========================================================================//

#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TParticle.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliPID.h"
#include "AliAODPid.h"                
#include "AliPIDResponse.h"
#include "AliAODpidUtil.h"        
#include "AliPIDCombined.h"

#include "AliHigherMomentsToyModel.h"

using std::cout;
using std::endl;


ClassImp(AliHigherMomentsToyModel)


//-----------------------------------------------------------------------
AliHigherMomentsToyModel::AliHigherMomentsToyModel( const char *name )
: AliAnalysisTaskSE( name ),
  fListOfHistosQA(0),
  fListOfHistos(0),
  fAOD(0),
  fArrayMC(0),
  fPIDResponse(0),
  fParticleSpecies(AliPID::kProton),
  fAnalysisType("AOD"),
  fCentralityEstimator("V0M"),
  fCentrality(0),
  fVxMax(3.),
  fVyMax(3.),
  fVzMax(10.),
  fDCAxy(2.4),
  fDCAz(3.),
  fPtLowerLimit(0.3),
  fPtHigherLimit(1.5),
  fEtaLowerLimit(-0.8),
  fEtaHigherLimit(0.8),
  fRapidityCut(0.5),
  fNSigmaCut(3.),
  fTPCNClus(80),
  fChi2perNDF(4.),
  fAODtrackCutBit(1024),
  fUsePid(kFALSE),
  fEventCounter(0),
  fHistDCA(0),
  fTPCSig(0),
  fTPCSigA(0),
  fTHnCentNplusNminusCh(0),
  fTHnCentNplusNminusChTruth(0),
  fTHnCentNplusNminus(0)
{ 
  
  for ( Int_t i = 0; i < 13; i++) { 
    fHistQA[i] = NULL;
  }
  
  for ( Int_t i = 0; i < 5; i++ ){
    fTHnCentNplusNminusPid[i] = NULL;
    fTHnCentNplusNminusPidTruth[i] = NULL;
  }

  DefineOutput(1, TList::Class()); // Outpput....
  DefineOutput(2, TList::Class()); 

}

AliHigherMomentsToyModel::~AliHigherMomentsToyModel() {
  
  if(fListOfHistosQA) delete fListOfHistosQA;
  if(fListOfHistos) delete fListOfHistos;

}

//---------------------------------------------------------------------------------
void AliHigherMomentsToyModel::UserCreateOutputObjects() {
  //For QA-Histograms
  fListOfHistosQA = new TList();
  fListOfHistosQA->SetOwner(kTRUE);
  
  fListOfHistos = new TList();
  fListOfHistos->SetOwner(kTRUE);
  
  fEventCounter = new TH1D("fEventCounter","EventCounter", 10, 0.5,10.5);
  fEventCounter->GetXaxis()->SetBinLabel(1,"Event Accesed");
  fEventCounter->GetXaxis()->SetBinLabel(2,"Within 0-90% centrality");
  fEventCounter->GetXaxis()->SetBinLabel(5,"Have a vertex");
  fEventCounter->GetXaxis()->SetBinLabel(6,"After vertex Cut");
  fEventCounter->GetXaxis()->SetBinLabel(7,"Event Analyzed");
  fEventCounter->GetXaxis()->SetBinLabel(8,"Event Analysis finished");
  fListOfHistosQA->Add(fEventCounter);
  
  fHistQA[0] = new TH1D("fHistQAvx", "Histo Vx After Cut", 400, -4., 4.);
  fHistQA[1] = new TH1D("fHistQAvy", "Histo Vy After Cut", 400, -4., 4.);
  fHistQA[2] = new TH1D("fHistQAvz", "Histo Vz After Cut", 500, -25., 25.);  
  fHistQA[3] = new TH1D("fHistQAvxA", "Histo Vx All", 500, -5., 5.);
  fHistQA[4] = new TH1D("fHistQAvyA", "Histo Vy All", 500, -5., 5.);
  fHistQA[5] = new TH1D("fHistQAvzA", "Histo Vz All", 500, -25., 25.);
  fHistQA[6] = new TH1D("fHistQADcaXyC", "Histo DCAxy after cut", 500, -5., 5.);
  fHistQA[7] = new TH1D("fHistQADcaZC", "Histo DCAz after cut", 500, -5., 5.);   
  fHistQA[8] = new TH1D("fHistQAPt","p_{T} distribution",600,0,6);
  fHistQA[9] = new TH1D("fHistQAEta","#eta distribution",240,-1.2,1.2);
  fHistQA[10] = new TH1D("fHistQAPhi","#phi distribution",340,0,6.8);
  fHistQA[11] = new TH1D("fHistQANCls","Number of TPC cluster",200,0,200);
  fHistQA[12] = new TH1D("fHistQAChi2","Chi2 per NDF",100,0,10);
  
  for(Int_t i = 0; i < 13; i++)
    {
      fListOfHistosQA->Add(fHistQA[i]);
    }
  
  fHistDCA = new TH2D("fHistDCA","DCAxy Vs DCAz", 500, -5., 5., 500, -5., 5.);
  fTPCSig = new TH2D("fTPCSig","TPC signal",200, 0.0, 10. ,1000,0.,1000);
  fTPCSig->SetMarkerColor(kRed);
  fTPCSigA = new TH2D("fTPCSigA","TPC signal all ",200, 0.0, 10. ,1000,0.,1000);
  fListOfHistosQA->Add(fHistDCA);
  fListOfHistosQA->Add(fTPCSig);
  fListOfHistosQA->Add(fTPCSigA);
  
  
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
  
  Int_t fBinsCh[nDim] = {100, 1500, 1500};
  Double_t fMinCh[nDim] = { -0.5, -0.5, -0.5 };
  Double_t fMaxCh[nDim] = { 99.5, 1499.5, 1499.5};
  
  fTHnCentNplusNminusCh = new THnSparseD("fTHnCentNplusNminusCh","Cent-NplusChrg-NminusChrg", nDim, fBinsCh, fMinCh, fMaxCh); 
  fTHnCentNplusNminusCh->GetAxis(0)->SetTitle("Centrality");
  fTHnCentNplusNminusCh->GetAxis(1)->SetTitle("Nplus");
  fTHnCentNplusNminusCh->GetAxis(2)->SetTitle("Nminus");
  fListOfHistos->Add(fTHnCentNplusNminusCh);

  if( fAnalysisType == "MCAOD" ){   
    
    fTHnCentNplusNminusChTruth = new THnSparseD("fTHnCentNplusNminusChTruth","Cent-NplusChrg-NminusChrg", nDim, fBinsCh, fMinCh, fMaxCh); 
    fTHnCentNplusNminusChTruth->GetAxis(0)->SetTitle("Centrality");
    fTHnCentNplusNminusChTruth->GetAxis(1)->SetTitle("Nplus");
    fTHnCentNplusNminusChTruth->GetAxis(2)->SetTitle("Nminus");
    fListOfHistos->Add(fTHnCentNplusNminusChTruth);
  }//MCAOD---
 
  
  TString hname1, hname11;
  TString htitle1, axisTitle1, axisTitle2;
  
  
  
  if( fUsePid ){
    
    Int_t gPid = (Int_t)fParticleSpecies;
    
    if( gPid > 1 && gPid < 5 ){ 
      //Pion----
      fBins[2][0] = 100; fBins[2][1] = 1000; fBins[2][2] = 1000;
      fMin[2][0] = -0.5; fMin[2][1] = -0.5; fMin[2][2] = -0.5;
      fMax[2][0] = 99.5; fMax[2][1] = 999.5; fMax[2][2] = 999.5;
      //Kaon------
      fBins[3][0] = 100; fBins[3][1] = 600; fBins[3][2] = 600;
      fMin[3][0] = -0.5; fMin[3][1] = -0.5; fMin[3][2] = -0.5;
      fMax[3][0] = 99.5; fMax[3][1] = 599.5; fMax[3][2] = 599.5;
      //Proton-----
      fBins[4][0] = 100; fBins[4][1] = 400; fBins[4][2] = 400;
      fMin[4][0] = -0.5; fMin[4][1] = -0.5; fMin[4][2] = -0.5;
      fMax[4][0] = 99.5; fMax[4][1] = 399.5; fMax[4][2] = 399.5;
      
      hname1 = "fCentNplusNminusPid"; hname1 +=gPid;
      htitle1 = Species[gPid]; htitle1 +=" And Neg-"; htitle1 +=Species[gPid];
      axisTitle1 = Species[gPid];
      axisTitle2 ="Neg-"; axisTitle2 += Species[gPid];
      
      fTHnCentNplusNminusPid[gPid] = new THnSparseD(hname1.Data(),htitle1.Data(), nDim, fBins[gPid], fMin[gPid],fMax[gPid]);
      
      fTHnCentNplusNminusPid[gPid]->GetAxis(0)->SetTitle("Centrality");
      fTHnCentNplusNminusPid[gPid]->GetAxis(1)->SetTitle(axisTitle1.Data());
      fTHnCentNplusNminusPid[gPid]->GetAxis(1)->SetTitle(axisTitle2.Data());
      
      fListOfHistos->Add(fTHnCentNplusNminusPid[gPid]);
      
      if( fAnalysisType == "MCAOD" ){
	
	hname11 = "fCentNplusNminusPidTruth"; hname11 +=gPid;
	fTHnCentNplusNminusPidTruth[gPid] = new THnSparseD(hname11.Data(),htitle1.Data(), nDim, fBins[gPid], fMin[gPid],fMax[gPid]);
	
	fTHnCentNplusNminusPidTruth[gPid]->GetAxis(0)->SetTitle("Centrality");
	fTHnCentNplusNminusPidTruth[gPid]->GetAxis(1)->SetTitle(axisTitle1.Data());
	fTHnCentNplusNminusPidTruth[gPid]->GetAxis(1)->SetTitle(axisTitle2.Data());
	
	fListOfHistos->Add(fTHnCentNplusNminusPidTruth[gPid]);
	
      }//MCAOD-----
    }//Pion, Koan and Proton-------
    
    else{
      
      Int_t fBinsX[nDim] = {100, 1500, 1500};
      Double_t fMinX[nDim] = { -0.5, -0.5, -0.5 };
      Double_t fMaxX[nDim] = { 99.5, 1499.5, 1499.5};
      
      fTHnCentNplusNminus = new THnSparseD("fTHnCentNplusNminus","Cent-NplusChrg-NminusChrg", nDim, fBinsX, fMinX, fMaxX); 
      fTHnCentNplusNminus->GetAxis(0)->SetTitle("Centrality");
      fTHnCentNplusNminus->GetAxis(1)->SetTitle("UnwantedPlus");
      fTHnCentNplusNminus->GetAxis(2)->SetTitle("UnwantedMinus");
      fListOfHistos->Add(fTHnCentNplusNminus);
      
    }//Unwanted particle -------
    
  }//fUsePid-------
  
  
  
  PostData(1, fListOfHistosQA);
  PostData(2, fListOfHistos);
  
}


//----------------------------------------------------------------------------------
void AliHigherMomentsToyModel::UserExec( Option_t * ){
  
  fEventCounter->Fill(1);

  
  if(fAnalysisType == "AOD") {
    
    doAODEvent();
    
  }//AOD--analysis-----

  else if(fAnalysisType == "MCAOD") {
  
    doMCAODEvent();
    
  }
  
  else return;
  
  
  
  
  fEventCounter->Fill(8);
  
  PostData(1, fListOfHistosQA);
  PostData(2, fListOfHistos);
  
}

//--------------------------------------------------------------------------------------
void AliHigherMomentsToyModel::doAODEvent(){
 
  //-------------------
  Double_t nPlusCharge = 0.;
  Double_t nMinusCharge = 0.;
  Double_t nPartile = 0.;
  Double_t nAntiParticle = 0.;
  Int_t gPid = 0.;


  //connect to the analysis mannager-----
 
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
  
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD)
    {
      cout<< "ERROR 01: AOD not found " <<endl;
      return;
    }
  
  fPIDResponse =inputHandler->GetPIDResponse();
  
  
  if (!fPIDResponse){
    AliFatal("This Task needs the PID response attached to the inputHandler");
    return;
  }  
  
  
  AliAODHeader *aodHeader = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
  if(!aodHeader) AliFatal("Not a standard AOD");
  
  fCentrality = (Int_t)aodHeader->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
  /* Int_t cent = -1;
  cent =  aodHeader->GetCentralityP()->GetCentralityClass10(fCentralityEstimator.Data());
  
  if (cent == 0)
    fCentrality = aodHeader->GetCentralityP()->GetCentralityClass5(fCentralityEstimator.Data());
  else if (cent == 10 || cent == -1.)
    fCentrality = -1;
  else if (cent > 0 && cent < 9)
    fCentrality = cent + 1;
  */
  if(fCentrality < 0 || fCentrality >= 91) return;
  
  fEventCounter->Fill(2);
  
  if(!ProperVertex()) return;   
  
  Int_t nTracks = fAOD->GetNumberOfTracks();
  
  TExMap *trackMap = new TExMap();//Mapping matrix----
  
  for(Int_t i = 0; i < nTracks; i++) {
    
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(i)); 
    
    if(!aodTrack) {
      AliError(Form("ERROR: Could not retrieve AODtrack %d",i));
      continue;
    }
    
    Double_t tpcSignalAll = aodTrack->GetTPCsignal();
    fTPCSigA->Fill(aodTrack->GetTPCmomentum(),tpcSignalAll);
    
    Int_t gID = aodTrack->GetID();
    
    if (aodTrack->TestFilterBit(1)) trackMap->Add(gID, i);//Global tracks-----
  }//1st track loop----
 
 AliAODTrack* newAodTrack; 
  
  for( Int_t j = 0; j < nTracks; j++ )
    {
      
      AliAODTrack* aodTrack1 = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(j));
      
      if(!aodTrack1) {
	AliError(Form("ERROR: Could not retrieve AODtrack %d",j));
	continue;
      }
      
      
      if(!aodTrack1->TestFilterBit(fAODtrackCutBit)) continue;
      
      Int_t gID = aodTrack1->GetID();
      
      //if( aodTrack->GetID() != (-aodTrack1->GetID() -1) ) continue;
      newAodTrack = gID >= 0 ? aodTrack1 : dynamic_cast<AliAODTrack*>(fAOD->GetTrack(trackMap->GetValue(-1-gID))); //Take those global track who corresponds to TPC only track
      if(!newAodTrack) AliFatal("Not a standard AOD");
      
      Float_t dxy = 0., dz = 0.;
      
      dxy = aodTrack1->DCA();
      dz  = aodTrack1->ZAtDCA();
      
      Double_t pt = aodTrack1->Pt();
      Double_t eta = aodTrack1->Eta();
      Double_t nclus = aodTrack1->GetTPCClusterInfo(2,1);
      Double_t chi2ndf = aodTrack1->Chi2perNDF();
      
      /* 
      if( fabs(dxy) > fDCAxy ) continue; 
      if( fabs(dz) > fDCAz ) continue;
      //Extra cut on DCA---( Similar to BF Task.. )	      
      if( fDCAxy !=0 && fDCAz !=0 ){
	if( TMath::Sqrt((dxy*dxy)/(fDCAxy*fDCAxy)+(dz*dz)/(fDCAz*fDCAz)) > 1. ) continue;
      }
      */
      if( pt < fPtLowerLimit || pt > fPtHigherLimit ) continue;
      if( eta < fEtaLowerLimit || eta > fEtaHigherLimit ) continue;
      if( nclus < fTPCNClus ) continue;
      if( chi2ndf > fChi2perNDF ) continue;
      
      
      fHistQA[6]->Fill(dxy);
      fHistQA[7]->Fill(dz);
      fHistQA[8]->Fill(pt);
      fHistQA[9]->Fill(eta);
      fHistQA[10]->Fill(aodTrack1->Phi());
      fHistQA[11]->Fill(nclus);
      fHistQA[12]->Fill(chi2ndf);
      fHistDCA->Fill(dxy,dz);
      
      Short_t gCharge = aodTrack1->Charge();
      
      if(gCharge > 0) nPlusCharge += 1.;
      if(gCharge < 0) nMinusCharge += 1.;
      
      if( fUsePid ) {
	
	gPid = (Int_t)fParticleSpecies;

	Double_t rap =  newAodTrack->Y(AliPID::ParticleMass(fParticleSpecies));
	Double_t tpcSignal = newAodTrack->GetTPCsignal();
	//Double_t rap =  aodTrack1->Y(AliPID::ParticleMass(fParticleSpecies));
	//Double_t tpcSignal = aodTrack1->GetTPCsignal();
	
	if( fabs(rap ) > fRapidityCut ) continue;//Rapidity cut
	
	fTPCSig->Fill(newAodTrack->GetTPCmomentum(),tpcSignal);
	
	Float_t nsigmaTPCPID = -999.;
	//Float_t nsigmaTOFPID = -999.;
	//Float_t nsigmaTPCTOFPID = -999.;
	
	nsigmaTPCPID = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(newAodTrack,fParticleSpecies));
	//nsigmaTOFPID = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(newAodTrack,fParticleSpecies));
	
	if ( nsigmaTPCPID < fNSigmaCut  ){
	  
	  if (gCharge > 0) nPartile +=1.;
	  if( gCharge < 0 ) nAntiParticle +=1.;
	  
	}
      }//fUsepid-----
      
    }//--------- Track Loop to select with filterbit
  
  
  
  Double_t fContainerCh[3] = { static_cast<Double_t>(fCentrality), nPlusCharge, nMinusCharge};
  Double_t fContainerPid[3] = { static_cast<Double_t>(fCentrality), nPartile, nAntiParticle};
  
  
  fTHnCentNplusNminusCh->Fill(fContainerCh);
  
  if( fUsePid ){
    gPid = (Int_t)fParticleSpecies;
    fTHnCentNplusNminusPid[gPid]->Fill(fContainerPid);
    
    // cout << "nCentrality "<< fCentrality <<", nParticle="<< nPartile << ", nMinusParticle=" << nAntiParticle << endl;
  }
  
  
  //cout << "nCentrality "<< fCentrality <<", nPlusCharge="<< nPlusCharge << ", nMinusCharge=" << nMinusCharge << endl;
  
  fEventCounter->Fill(7);
  
  return;
  
}
//--------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------
void AliHigherMomentsToyModel::doMCAODEvent(){
  
  
  //---------
  Double_t nPlusCharge = 0.;
  Double_t nMinusCharge = 0.;

  Double_t nPlusChargeTruth = 0.;
  Double_t nMinusChargeTruth = 0.;

  Double_t nPartile = 0.;
  Double_t nAntiParticle = 0.;
  Double_t nPartileTruth = 0.;
  Double_t nAntiParticleTruth = 0.;

  Int_t gPid = 0;
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
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD)
    {
      cout<< "ERROR 01: AOD not found " <<endl;
      return;
    }
  
  fPIDResponse =inputHandler->GetPIDResponse();
  
  
  if (!fPIDResponse){
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
  
  AliAODHeader *aodHeader = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
  if(!aodHeader) AliFatal("Not a standard AOD");
  
  fCentrality = (Int_t)aodHeader->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
  
  /*
    Int_t cent = -1;
    cent =  aodHeader->GetCentralityP()->GetCentralityClass10(fCentralityEstimator.Data());
    
    if (cent == 0)
    fCentrality = aodHeader->GetCentralityP()->GetCentralityClass5(fCentralityEstimator.Data());
    else if (cent == 10 || cent == -1.)
    fCentrality = -1;
    else if (cent > 0 && cent < 9)
    fCentrality = cent + 1;
  */
  if( fCentrality < 0 || fCentrality >= 91) return;
  
  fEventCounter->Fill(2);
  
  
  
  if(!ProperVertex()) return;   
  
  Int_t nTracks = fAOD->GetNumberOfTracks();
 
  
  TExMap *trackMap = new TExMap();//Mapping matrix----
  
  for(Int_t i = 0; i < nTracks; i++) {
    
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(i)); 
    
    if(!aodTrack) {
      AliError(Form("ERROR: Could not retrieve AODtrack %d",i));
      continue;
    }
    
    Double_t tpcSignalAll = aodTrack->GetTPCsignal();
    fTPCSigA->Fill(aodTrack->GetTPCmomentum(),tpcSignalAll);
    
    Int_t gID = aodTrack->GetID();
    
    if (aodTrack->TestFilterBit(1)) trackMap->Add(gID, i);//Global tracks-----
    
  }//1st track loop----
  
  AliAODTrack* newAodTrack; 
  
  for( Int_t j = 0; j < nTracks; j++ ){
    
    AliAODTrack* aodTrack1 = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(j));
    
    if(!aodTrack1) {
      AliError(Form("ERROR: Could not retrieve AODtrack %d",j));
      continue;
    }
    
   
    if(!aodTrack1->TestFilterBit(fAODtrackCutBit)) continue;
    
    
    Int_t gID = aodTrack1->GetID();
    
    //if( aodTrack->GetID() != (-aodTrack1->GetID() -1) ) continue;
    
    newAodTrack = gID >= 0 ? aodTrack1 : dynamic_cast<AliAODTrack*>(fAOD->GetTrack(trackMap->GetValue(-1-gID))); //Take those global track who corresponds to TPC only track
    if(!newAodTrack) AliFatal("Not a standard AOD");
    
  
    //cout << dxy << endl;
    Double_t pt = aodTrack1->Pt();
    Double_t eta = aodTrack1->Eta();
    Double_t nclus = aodTrack1->GetTPCClusterInfo(2,1);
    Double_t chi2ndf = aodTrack1->Chi2perNDF();
  
    Float_t dxy = 0., dz = 0.;
    
    if( fAODtrackCutBit == 128 ){
     
      dxy = aodTrack1->DCA();
      dz  = aodTrack1->ZAtDCA();      
      
      if( fabs(dxy) > fDCAxy ) continue; 
      if( fabs(dz) > fDCAz ) continue;
      //Extra cut on DCA---( Similar to BF Task.. )	      
      if( fDCAxy !=0 && fDCAz !=0 ){
	if( TMath::Sqrt((dxy*dxy)/(fDCAxy*fDCAxy)+(dz*dz)/(fDCAz*fDCAz)) > 1. ) continue;
      }

      fHistQA[6]->Fill(dxy);
      fHistQA[7]->Fill(dz);
      fHistDCA->Fill(dxy,dz);

    }
    

    if( pt < fPtLowerLimit || pt > fPtHigherLimit ) continue;
    if( eta < fEtaLowerLimit || eta > fEtaHigherLimit ) continue;
    if( nclus < fTPCNClus ) continue;
    if( chi2ndf > fChi2perNDF ) continue;
    
    
    
    fHistQA[8]->Fill(pt);
    fHistQA[9]->Fill(eta);
    fHistQA[10]->Fill(aodTrack1->Phi());
    fHistQA[11]->Fill(nclus);
    fHistQA[12]->Fill(chi2ndf);
    
    
    Short_t gCharge = aodTrack1->Charge();
    
    if( gCharge == 0 ) continue;
    
    
    if(gCharge > 0) nPlusCharge += 1.;
    if(gCharge < 0) nMinusCharge += 1.;
    
    if( fUsePid ) {
      
      gPid = (Int_t)fParticleSpecies;;
      Double_t rap =  newAodTrack->Y(AliPID::ParticleMass(fParticleSpecies));
      Double_t tpcSignal = newAodTrack->GetTPCsignal();
   
      if( fabs(rap ) > fRapidityCut ) continue;//Rapidity cut
      
      fTPCSig->Fill(newAodTrack->GetTPCmomentum(),tpcSignal);
      
      Float_t nsigmaTPCPID = -999.;
      //Float_t nsigmaTOFPID = -999.;
      //Float_t nsigmaTPCTOFPID = -999.;
      
      nsigmaTPCPID = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(newAodTrack,fParticleSpecies));
      //nsigmaTOFPID = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(newAodTrack,fParticleSpecies));
      
      if ( nsigmaTPCPID < fNSigmaCut  ){

	
	if (gCharge > 0) nPartile +=1.;
	if( gCharge < 0 ) nAntiParticle +=1.;
	
      }//nSigmaCut-----
    }//fUsepid-----
   
  }//--------- Track Loop to select with filterbit
  
  //cout << "Cent=" << fCentrality << " MC-PlusChrg=" << nPlusCharge << " MC-MinusChrg=" << nMinusCharge << endl;
  
  //cout << "nCentrality "<< fCentrality <<", nParticle="<< nPartile << ", nMinusParticle=" << nAntiParticle << endl; 
  
  Double_t fContainerCh[3] = { static_cast<Double_t>(fCentrality), nPlusCharge, nMinusCharge};
  Double_t fContainerPid[3] = { static_cast<Double_t>(fCentrality), nPartile, nAntiParticle};
  
  fTHnCentNplusNminusCh->Fill(fContainerCh);  
  
  if( fUsePid ){

    gPid = (Int_t)fParticleSpecies;
    fTHnCentNplusNminusPid[gPid]->Fill(fContainerPid);
    
  }
  
  //===========================================
  //--------MC Truth---------------------------
  //===========================================

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

    if(gCharge > 0) nPlusChargeTruth += 1.;
    if(gCharge < 0) nMinusChargeTruth += 1.;
   
    if(fUsePid){
      
      Double_t rap =  partMC->Y();
      
      if( fabs(rap ) > fRapidityCut ) continue;//Rapidity cut
      
      if(TMath::Abs(partMC->GetPdgCode()) != gPdgCode) continue;
      
      if(gCharge > 0) nPartileTruth += 1.;
      if(gCharge < 0) nAntiParticleTruth += 1.;
      
    }//if(fUsePid) ----
  
  }//MC-Truth Track loop--
  
  Double_t fContainerChTruth[3] = { static_cast<Double_t>(fCentrality), nPlusChargeTruth, nMinusChargeTruth };
  Double_t fContainerPidTruth[3] = { static_cast<Double_t>(fCentrality), nPartileTruth, nAntiParticleTruth };
  
  //cout << "Cent=" << fCentrality << " MC-PlusChrgT=" << nPlusChargeTruth << " MC-MinusChrgT=" << nMinusChargeTruth << endl;
  
  //cout << "nCentrality "<< fCentrality <<", nParticleT="<< nPartileTruth << ", nMinusParticleT=" << nAntiParticleTruth << endl; 
  
  //cout << "----------------------------" << endl;
  fTHnCentNplusNminusChTruth->Fill(fContainerChTruth);  
  
  if( fUsePid ){
    gPid = (Int_t)fParticleSpecies;
    fTHnCentNplusNminusPidTruth[gPid]->Fill(fContainerPidTruth);
    
  }
  
  
  fEventCounter->Fill(7);
  
  return;
  
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
Bool_t AliHigherMomentsToyModel::ProperVertex(){
  
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
	
	fHistQA[3]->Fill(lvx);fHistQA[4]->Fill(lvy);fHistQA[5]->Fill(lvz);
	
	if(TMath::Abs(lvx) < fVxMax) {
	  if(TMath::Abs(lvy) < fVyMax) {
	    if(TMath::Abs(lvz) < fVzMax) {
	      
	      fEventCounter->Fill(6);
	      fHistQA[0]->Fill(lvx);fHistQA[1]->Fill(lvy);fHistQA[2]->Fill(lvz);
	      
	      IsVtx = kTRUE; 
	      
	    }//Z-Vertex cut---
	  }//Y-vertex cut--
	}//X-vertex cut---
      }//Covariance------
    }//Contributors check---
  }//If vertex-----------
  
  return IsVtx;
}
//---------------------------------------------------------------------------------------



//------------------------------------------------------------------------
void AliHigherMomentsToyModel::Terminate( Option_t * ){
  
  Info("AliHigherMomentsToyModel"," Task Successfully finished");
  
}

//------------------------------------------------------------------------
