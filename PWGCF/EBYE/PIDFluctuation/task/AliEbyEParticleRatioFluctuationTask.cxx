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
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//
//             V0.1 2013/03/25 Using THnSparse
//             V0.2 2013/04/03 Cleanup
//             V1.0 2013/04/10 Cleanup Bins for Less Memory
//             V1.1 2013/04/15 Bins Added 
//   Todo: pp and pA, Mix Events
//=========================================================================//

#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliPIDResponse.h"
#include "AliAODHeader.h"
#include "AliAODpidUtil.h"
#include "AliHelperPID.h"

#include "AliEbyEParticleRatioFluctuationTask.h"

ClassImp(AliEbyEParticleRatioFluctuationTask)

//-----------------------------------------------------------------------
AliEbyEParticleRatioFluctuationTask::AliEbyEParticleRatioFluctuationTask( const char *name ) : AliAnalysisTaskSE( name ), 
  fThnList(NULL), 
  fAnalysisType("AOD"), 
  fAnalysisData("PbPb"), 
  fCentralityEstimator("V0M"), 
  fVxMax(3.), 
  fVyMax(3.), 
  fVzMax(10.), 
  fDCAxy(2.4),  
  fDCAz(3.2), 
  fPtLowerLimit(0.2), 
  fPtHigherLimit(5.), 
  fEtaLowerLimit(-1.), 
  fEtaHigherLimit(1.), 
  fTPCNClus(80),
  fAODtrackCutBit(128), 
  isQA(kFALSE), 
  fDebug(kFALSE), 
  fHelperPID(0x0),
  fEventCounter(NULL), 
  fHistoCorrelation(NULL) { 
  for(Int_t i = 0; i < 14; i++ ) fHistQA[i] = NULL;
  DefineOutput(1, TList::Class()); //! Connect Outpput....
}

AliEbyEParticleRatioFluctuationTask::~AliEbyEParticleRatioFluctuationTask() {
  //!   Cleaning up
  if (fThnList)   delete fThnList;
  if (fHelperPID) delete fHelperPID;
}

//---------------------------------------------------------------------------------
void AliEbyEParticleRatioFluctuationTask::UserCreateOutputObjects() {
  fThnList = new TList();
  fThnList->SetOwner(kTRUE);

  fEventCounter = new TH1D("fEventCounter","EventCounter", 300, 0.5,300.5);
  if (isQA) fThnList->Add(fEventCounter);
  
  fHistQA[0] = new TH2F("fHistQAvx", "Histo Vx Selected;Centrality;Vx", 100,0,100, 5000, -5., 5.);
  fHistQA[1] = new TH2F("fHistQAvy", "Histo Vy Selected;Centrality;Vy", 100,0,100, 5000, -5., 5.);
  fHistQA[2] = new TH2F("fHistQAvz", "Histo Vz Selected;Centrality;Vz", 100,0,100, 5000, -25., 25.);  
  fHistQA[3] = new TH2F("fHistQAvxA", "Histo Vx;Centrality;Vx", 100,0,100, 5000, -5., 5.);
  fHistQA[4] = new TH2F("fHistQAvyA", "Histo Vy;Centrality;Vy", 100,0,100, 5000, -5., 5.);
  fHistQA[5] = new TH2F("fHistQAvzA", "Histo Vz;Centrality;Vz", 100,0,100, 5000, -25., 25.);

  fHistQA[6] = new TH2F("fHistQADcaXyA", "Histo DCAxy;Centrality;DCAxy",100,0,100, 600, -15., 15.);
  fHistQA[7] = new TH2F("fHistQADcaZA", "Histo DCAz;Centrality;DCAz ",100,0,100, 600, -15., 15.);   
  fHistQA[8] = new TH2F("fHistQAPtA","p_{T} distribution;Centrality;p_{T}",100,0,100,1000,0,10);
  fHistQA[9] = new TH2F("fHistQAEtaA","#eta distribution;Centrality;#eta",100,0,100,240,-1.2,1.2);

  fHistQA[10] = new TH2F("fHistQADcaXy", "Histo DCAxy after Selected;Centrality;DCAxy", 100,0,100,600, -15., 15.);
  fHistQA[11] = new TH2F("fHistQADcaZ", "Histo DCAz Selected;Centrality;DCAz", 100,0,100,600, -15., 15.);   
  fHistQA[12] = new TH2F("fHistQAPt","p_{T} distribution Selected;Centrality;p_{T}",100,0,100,1000,0,10);
  fHistQA[13] = new TH2F("fHistQAEta","#eta distribution Selected;Centrality;#eta",100,0,100, 240,-1.2,1.2);
  
  if (isQA) for(Int_t i = 0; i < 14; i++) fThnList->Add(fHistQA[i]);
  
  Int_t fgSparseDataBins[kNSparseData]   = {100, 5000, 5000, 2500, 2500, 3000, 1500, 1500, 1000, 500, 500, 500, 250, 250};
  Double_t fgSparseDataMin[kNSparseData] = {0.,  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  0.,  0.,  0.,  0.};
  Double_t fgSparseDataMax[kNSparseData] = {100.,5000.,5000.,2500.,2500.,3000.,1500.,1500.,1000.,500.,500.,500.,250.,250.};
  
  const Char_t *fgkSparseDataTitle[] = {"centrality","RefMult","N_{ch}", "N_{+}","N_{-}","N_{#pi}", "N_{#pi^{+}}","N_{#pi^{-}}","N_{K}","N_{K^{+}}", "N_{K^{-}}","N_{p}","N_{p}","N_{#bar{p}}"};
    
  fHistoCorrelation = new THnSparseI("fThnCorr", "", kNSparseData, fgSparseDataBins, fgSparseDataMin, fgSparseDataMax);
  for (Int_t iaxis = 0; iaxis < kNSparseData; iaxis++)
    fHistoCorrelation->GetAxis(iaxis)->SetTitle(fgkSparseDataTitle[iaxis]);
  
  if(!isQA) fThnList->Add(fHistoCorrelation);
  if(isQA)  
    if (fHelperPID)
      fThnList->Add(fHelperPID);
  
  PostData(1, fThnList);
}

//----------------------------------------------------------------------------------
void AliEbyEParticleRatioFluctuationTask::UserExec( Option_t * ){

  if (isQA) fEventCounter->Fill(1);

  AliAODEvent* event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event) {
    cout<< "ERROR 01: AOD not found " <<endl;
    return;
  }

  Int_t gCent   = -1;
  Float_t gRefMul = -1;
  
  AliAODHeader *aodHeader = event->GetHeader();
  gCent = (Int_t)aodHeader->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
  gRefMul = aodHeader->GetRefMultiplicity();
  if (gCent < 0 || gCent > 100) return;
  if (isQA) fEventCounter->Fill(2);  

  if (!AcceptEvent(event,gCent)) return;
  
  Int_t icharge = -1;
  Int_t gCharge[2];
  Int_t gPid[3][2];
  
  for (Int_t i = 0; i < 2; i++) {
    gCharge[i] = 0;
    for (Int_t ii = 0; ii < 3; ii++) {
      gPid[ii][i] = 0;
    }
  }


  if(fAnalysisType == "AOD") {
    if (isQA) {
	fEventCounter->Fill(5);
	fEventCounter->Fill(50+gCent);
    }
    
    for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {
      AliAODTrack* track = dynamic_cast<AliAODTrack *>(event->GetTrack(itrk));
      if (!track) continue;
      
      if (!AcceptTrack(track, gCent)) continue;
            
      Int_t a = fHelperPID->GetParticleSpecies((AliVTrack*) track,kTRUE);

      if(a < 0 || a > 2) continue;
      icharge = track->Charge() > 0 ? 0 : 1;
      gCharge[icharge]++;
      gPid[a][icharge]++;
    }
  }
  else if(fAnalysisType == "MCAOD") {
    TClonesArray *arrayMC= 0; 
    arrayMC = dynamic_cast<TClonesArray*> (event->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
    if (!arrayMC) {
      Printf("Error: MC particles branch not found!\n");
      return;
    }
    AliAODMCHeader *mcHdr=0;
    mcHdr=(AliAODMCHeader*)event->GetList()->FindObject(AliAODMCHeader::StdBranchName());  
    if(!mcHdr) {
      Printf("MC header branch not found!\n");
      return;
    }
    
    if (isQA) {
      fEventCounter->Fill(5);
      fEventCounter->Fill(50+gCent);
    }
    
    Int_t nMC = arrayMC->GetEntries();
    for (Int_t iMC = 0; iMC < nMC; iMC++) {
      AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(iMC);
      if(!AcceptMCTrack(partMC, gCent)) continue;
      Int_t a = fHelperPID->GetMCParticleSpecie((AliVEvent*) event,(AliVTrack*)partMC,1);
      if(a < 0 || a > 2) continue;
      icharge = partMC->Charge() > 0 ? 0 : 1;
      gCharge[icharge]++;
      gPid[a][icharge]++;
    }
  }
  else {
    printf(" No Event Type is Defined ");
    return;
  }
  
  if( (gCharge[0] + gCharge[1]) != 0 ) {
    if (isQA) {
      fEventCounter->Fill(6); 
      fEventCounter->Fill(160 + gCent);
    }
    else { 
      Double_t vsparse[kNSparseData];
      vsparse[0]   = gCent;
      vsparse[1]   = gRefMul;
      vsparse[2]   = gCharge[0] + gCharge[1];
      vsparse[3]   = gCharge[0];
      vsparse[4]   = gCharge[1];
      vsparse[5]   = gPid[0][0] + gPid[0][1];
      vsparse[6]   = gPid[0][0];
      vsparse[7]   = gPid[0][1];
      vsparse[8]   = gPid[1][0] + gPid[1][0];
      vsparse[9]   = gPid[1][0];
      vsparse[10]  = gPid[1][1];
      vsparse[11]  = gPid[2][0] + gPid[2][1];
      vsparse[12]  = gPid[2][0];
      vsparse[13]  = gPid[2][1];
      fHistoCorrelation->Fill(vsparse);
    }
  }
  
  if(fDebug && isQA)  Printf(" %6d  %6d %6d  %6d %6d  %6d %6d %6d %6d %6d %6d %6d %6d", 
			     (Int_t)fEventCounter->GetBinContent(1), 
			     (Int_t)fEventCounter->GetBinContent(2), 
			     (Int_t)fEventCounter->GetBinContent(3), 
			     (Int_t)gCent, (Int_t)gRefMul, 
			     gCharge[0], gCharge[1], 
			     gPid[0][0], gPid[0][1],  gPid[1][0], 
			     gPid[1][1],  gPid[2][0], gPid[2][1]);
  
  PostData(1, fThnList);
}

void AliEbyEParticleRatioFluctuationTask::Terminate( Option_t * ){
  Info("AliEbyEParticleRatioFluctuationTask"," Task Successfully finished");
}

//___________________________________________________________
Bool_t AliEbyEParticleRatioFluctuationTask::AcceptEvent(AliAODEvent *event, Int_t cent) const {
  Bool_t ver = kFALSE;
  const AliAODVertex *vertex = event->GetPrimaryVertex();
  if(vertex) {
    Double32_t fCov[6];
    vertex->GetCovarianceMatrix(fCov);
    if(vertex->GetNContributors() > 0) {
      if(fCov[5] != 0) {
	
	if(isQA) {	
	  fEventCounter->Fill(3);
	  fHistQA[3]->Fill(cent,vertex->GetX());
	  fHistQA[4]->Fill(cent,vertex->GetY());
	  fHistQA[5]->Fill(cent,vertex->GetZ());
	}
	
	if(TMath::Abs(vertex->GetX()) < fVxMax) {
	  if(TMath::Abs(vertex->GetY()) < fVyMax) {
	    if(TMath::Abs(vertex->GetZ()) < fVzMax) {
	      ver = kTRUE;
	      if(isQA) {	
		fEventCounter->Fill(4);
		fHistQA[0]->Fill(cent,vertex->GetX());
		fHistQA[1]->Fill(cent,vertex->GetY());
		fHistQA[2]->Fill(cent,vertex->GetZ());
	      }
	    }
	  }
	}
      }
    }
  }
  
  AliCentrality *centrality = event->GetCentrality();
  if (centrality->GetQuality() != 0) ver = kFALSE;
  return ver;
}


//___________________________________________________________
Bool_t AliEbyEParticleRatioFluctuationTask::AcceptTrack(AliAODTrack *track, Int_t cent) const {
  if(!track)                                  return kFALSE;
  if (track->Charge() == 0 )                  return kFALSE;
  
  if(isQA) {
    fHistQA[6]->Fill(cent,track->DCA());
    fHistQA[7]->Fill(cent,track->ZAtDCA());
    fHistQA[8]->Fill(cent,track->Pt());
    fHistQA[9]->Fill(cent,track->Eta());
  }
  
  if (!track->TestFilterBit(fAODtrackCutBit)) return kFALSE;

  if (track->Eta() < fEtaLowerLimit ||
      track->Eta() > fEtaHigherLimit)         return kFALSE;
  if (track->Pt() < fPtLowerLimit ||
      track->Pt() > fPtHigherLimit)           return kFALSE;  
  if ( track->DCA() > fDCAxy )                return kFALSE; 
  if ( track->ZAtDCA() > fDCAz )              return kFALSE;
   
  if(isQA) {
    fHistQA[10]->Fill(cent,track->DCA());
    fHistQA[11]->Fill(cent,track->ZAtDCA());
    fHistQA[12]->Fill(cent,track->Pt());
    fHistQA[13]->Fill(cent,track->Eta());
  }
 
  return kTRUE;
}


//___________________________________________________________
Bool_t AliEbyEParticleRatioFluctuationTask::AcceptMCTrack(AliAODMCParticle *track, Int_t cent) const {
  if(!track) return kFALSE;
  if(!track->IsPhysicalPrimary()) return kFALSE;
  if (track->Charge() == 0 ) return kFALSE;
  if(isQA) {
    fHistQA[8]->Fill(cent,track->Pt());
    fHistQA[9]->Fill(cent,track->Eta());
  }

  if (track->Eta() < fEtaLowerLimit ||
      track->Eta() > fEtaHigherLimit) return kFALSE;
  if (track->Pt() < fPtLowerLimit ||
      track->Pt() > fPtHigherLimit) return kFALSE;  
    
  if(isQA) {
    fHistQA[12]->Fill(cent,track->Pt());
    fHistQA[13]->Fill(cent,track->Eta());
  }
 
  return kTRUE;
}
