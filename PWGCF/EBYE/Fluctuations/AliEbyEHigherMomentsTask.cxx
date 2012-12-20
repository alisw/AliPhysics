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
#include "TH2F.h"
#include "TCanvas.h"

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
#include "AliStack.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliPID.h"
#include "AliAODPid.h"                
#include "AliPIDResponse.h"
#include "AliAODpidUtil.h"        
#include "AliPIDCombined.h"

#include "AliEbyEHigherMomentsTask.h"

using std::cout;
using std::endl;

ClassImp(AliEbyEHigherMomentsTask)


//-----------------------------------------------------------------------
AliEbyEHigherMomentsTask::AliEbyEHigherMomentsTask( const char *name )
: AliAnalysisTaskSE( name ),
  fListOfHistos(0),
  fAnalysisType(0),
  fCentralityEstimator("V0M"),
  fVxMax(3.),
  fVyMax(3.),
  fVzMax(15.),
  fDCAxy(2.4),
  fDCAz(3.),
  fPtLowerLimit(0.3),
  fPtHigherLimit(1.5),
  fEtaLowerLimit(-0.8),
  fEtaHigherLimit(0.8),
  fTPCNClus(80),
  fChi2perNDF(4.),
  nAODtrackCutBit(128),
  fEventCounter(0)
{ 
  
  for(Int_t i = 0; i < 91; i++)
    {
      fhNplusNminus[i] = NULL;
    }
  
  for ( Int_t i = 0; i < 13; i++) { 
    fHistQA[i] = NULL;
  }
  
  DefineOutput(1, TList::Class()); // Outpput....
}

AliEbyEHigherMomentsTask::~AliEbyEHigherMomentsTask() {
  if(fListOfHistos) delete fListOfHistos;
}

//---------------------------------------------------------------------------------
void AliEbyEHigherMomentsTask::UserCreateOutputObjects() {
  //For QA-Histograms
  fListOfHistos = new TList();
  fListOfHistos->SetOwner(kTRUE);
  fEventCounter = new TH1D("fEventCounter","EventCounter", 150, 0.5,150.5);
  fListOfHistos->Add(fEventCounter);
  
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
      fListOfHistos->Add(fHistQA[i]);
    }
  
  
  TString hname;
  TString htitle;
  
  for(Int_t i = 0; i < 25; i++)  {
    hname  = "hPlusMinusCentBin"; hname+=i;
    htitle = "N+ and N- in Cent Bin "; htitle+=i;
    fhNplusNminus[i] = new TH2F(hname.Data(),htitle.Data(),1400,0.5,1400.5,1400,0.5,1400.5);
    fListOfHistos->Add(fhNplusNminus[i]);
  }
  
  for(Int_t i = 25; i < 50; i++) {
    hname  = "hPlusMinusCentBin"; hname+=i;
    htitle = "N+ and N- in Cent Bin "; htitle+=i;
    fhNplusNminus[i] = new TH2F(hname.Data(),htitle.Data(),800,0.5,800.5,800,0.5,800.5);
    fListOfHistos->Add(fhNplusNminus[i]);
  }
  
  for(Int_t i = 50; i < 91; i++) {
    hname  = "hPlusMinusCentBin"; hname+=i;
    htitle = "N+ and N- in Cent Bin "; htitle+=i;
    fhNplusNminus[i] = new TH2F(hname.Data(),htitle.Data(),500,0.5,500.5,500,0.5,500.5);
    fListOfHistos->Add(fhNplusNminus[i]);
  }
  
  PostData(1, fListOfHistos);
}


//----------------------------------------------------------------------------------
void AliEbyEHigherMomentsTask::UserExec( Option_t * ){
  
  fEventCounter->Fill(1);
  Double_t nPlusCharge = 0.;
  Double_t nMinusCharge = 0.;
  
  AliAODEvent* gAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!gAOD) {
    cout<< "ERROR 01: AOD not found " <<endl;
    return;
  }
  
  
  if(fAnalysisType == "AOD") {
    AliAODHeader *aodHeader = gAOD->GetHeader();
    Double_t cent =  aodHeader->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
    Double_t count = (Int_t)cent + 20;
    fEventCounter->Fill(count);
    
    if(cent < 0 || cent >= 91) return;
    
    fEventCounter->Fill(2);
    
    Int_t nCentrality = (Int_t)cent;
    
    const AliAODVertex *vertex = gAOD->GetPrimaryVertex();
    
    if(vertex) {
      Double32_t fCov[6];
      vertex->GetCovarianceMatrix(fCov);
      if(vertex->GetNContributors() > 0){
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
		
		Int_t ntracks = gAOD->GetNumberOfTracks();
		
		for(Int_t i = 0; i < ntracks; i++) {
		  AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(gAOD->GetTrack(i)); 
		  if(!aodTrack) {
		    AliError(Form("ERROR: Could not retrieve AODtrack %d",i));
		    continue;
		  }
		  
		  if(!aodTrack->TestFilterBit(nAODtrackCutBit)) continue;
		  
		  Float_t dxy, dz;
		  
		  dxy = aodTrack->DCA();
		  dz  = aodTrack->ZAtDCA();
		  
		  Double_t pt = aodTrack->Pt();
		  Double_t eta = aodTrack->Eta();
		  //Double_t nclus = aodTrack->GetTPCNcls();
		  Double_t nclus = aodTrack->GetTPCClusterInfo(2,1);//important for 2011 data
		  Double_t chi2ndf = aodTrack->Chi2perNDF(); //
		  
		  //Set the track cuts--------------
		  if( dxy > fDCAxy ) continue; 
		  if( dz > fDCAz ) continue;
		  if( pt < fPtLowerLimit || pt > fPtHigherLimit ) continue;
		  if( eta < fEtaLowerLimit || eta > fEtaHigherLimit ) continue;
		  if( nclus < fTPCNClus ) continue;
		  if( chi2ndf > fChi2perNDF ) continue;
		  
		  //Fill the Trackwise QA-histograms-----
		  fHistQA[6]->Fill(dxy);
		  fHistQA[7]->Fill(dz);
		  fHistQA[8]->Fill(pt);
		  fHistQA[9]->Fill(eta);
		  fHistQA[10]->Fill(aodTrack->Phi());
		  fHistQA[11]->Fill(nclus);
		  fHistQA[12]->Fill(chi2ndf);
		  
		  //Count the charge particles-----
		  Short_t gCharge = aodTrack->Charge();
		  
		  if(gCharge > 0) nPlusCharge += 1.;
		  if(gCharge < 0) nMinusCharge += 1.;
		  
		}//--------- Track Loo
		//cout << "nCentrality "<<nCentrality<<" nPlusCharge "<< nPlusCharge << " nMinusCharge " << nMinusCharge << endl;
		fhNplusNminus[nCentrality]->Fill(nPlusCharge,nMinusCharge);
		fEventCounter->Fill(7);
	      }//Z-vertex cut-------
	    }//Vetecx-Y cut---
	  }//Vertex-X cut
	}//Covariant matrix cuts---
      }//N-Contributors-----------
    }//Check for vertex--------
  }//AOD--analysis-----
  
  else if(fAnalysisType == "MCAOD") {
    TClonesArray *arrayMC = dynamic_cast<TClonesArray*>(gAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!arrayMC) {
      AliError("Error: MC particles branch not found!\n");
      return;
    }
    AliAODMCHeader *mcHdr=0;
    mcHdr=(AliAODMCHeader*)gAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());  
    if(!mcHdr) {
      printf("MC header branch not found!\n");
      return;
    }
    
    AliAODHeader *aodHeader = gAOD->GetHeader();
    AliCentrality* fcentrality = aodHeader->GetCentralityP();
    Double_t cent =fcentrality->GetCentralityPercentile(fCentralityEstimator.Data());
    
    Int_t nCentrality = (Int_t)cent;
    if(cent < 0 || cent >= 91) return;
    
    Bool_t ver = kFALSE;
    const AliAODVertex *vertex = gAOD->GetPrimaryVertex();
    if(vertex) {
      Double32_t fCov[6];
      vertex->GetCovarianceMatrix(fCov);
      if(vertex->GetNContributors() > 0) {
	if(fCov[5] != 0) {
	  
	  if(TMath::Abs(vertex->GetX()) < fVxMax) {
	    if(TMath::Abs(vertex->GetY()) < fVyMax) {
	      if(TMath::Abs(vertex->GetZ()) < fVzMax) {
		ver = kTRUE;
	      }
	    }
	  }
	}
      }
    }
    
    if(!ver) return;
    Int_t nMC = arrayMC->GetEntries();
    for (Int_t iMC = 0; iMC < nMC; iMC++) {
      AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(iMC);
      if(!partMC){
	AliError(Form("ERROR: Could not retrieve AODMCtrack %d",iMC));
	continue;
      }
      if(!partMC->IsPhysicalPrimary())continue;
      if (partMC->Eta() < fEtaLowerLimit || partMC->Eta() > fEtaHigherLimit) continue;
      if (partMC->Pt() < fPtLowerLimit ||  partMC->Pt() > fPtHigherLimit) continue;
      
      if(partMC->Charge() > 0) nPlusCharge += 1.;
      if(partMC->Charge() < 0) nMinusCharge += 1.;
      
    }//MC Track loop-----
    fhNplusNminus[nCentrality]->Fill(nPlusCharge,nMinusCharge);
    //cout << "Cent=" << nCentrality << " MC-PlusChrg=" << nPlusCharge << " MC-MinusChrg=" << nMinusCharge << endl;
    
  }//MCAOD--analysis---- 
  
  else return;
    
  fEventCounter->Fill(8);
  PostData(1, fListOfHistos);
  
}


//------------------------------------------------------------------------
void AliEbyEHigherMomentsTask::Terminate( Option_t * ){
  
  Info("AliEbyEHigherMomentTask"," Task Successfully finished");
  
}

//------------------------------------------------------------------------
