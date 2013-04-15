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
AliEbyEParticleRatioFluctuationTask::AliEbyEParticleRatioFluctuationTask( const char *name ) : AliAnalysisTaskSE( name ), fThnList(NULL), fAnalysisType("AOD"), fAnalysisData("PbPb"), fCentralityEstimator("V0M"), fVxMax(3.), fVyMax(3.), fVzMax(10.), fDCAxy(2.4),  fDCAz(3.2), fPtLowerLimit(0.2), fPtHigherLimit(5.), fEtaLowerLimit(-1.), fEtaHigherLimit(1.), fAODtrackCutBit(128), isQA(kFALSE), fDebug(kFALSE), fHelperPID(0x0),fEventCounter(NULL), fHistoCorrelation(NULL) { 
  for(Int_t i = 0; i < 14; i++ ) fHistQA[i] = NULL;
  DefineOutput(1, TList::Class()); // Outpput....
}

AliEbyEParticleRatioFluctuationTask::~AliEbyEParticleRatioFluctuationTask() {
  //clean up
  if (fThnList) delete fThnList;
  if (fHelperPID) delete fThnList;
}

//---------------------------------------------------------------------------------
void AliEbyEParticleRatioFluctuationTask::UserCreateOutputObjects() {
  fThnList = new TList();
  fThnList->SetOwner(kTRUE);

  fEventCounter = new TH1D("fEventCounter","EventCounter", 250, 0.5,250.5);
  if (isQA) fThnList->Add(fEventCounter);
  
  fHistQA[0] = new TH1D("fHistQAvx", "Histo Vx Selected", 5000, -5., 5.);
  fHistQA[1] = new TH1D("fHistQAvy", "Histo Vy Selected", 5000, -5., 5.);
  fHistQA[2] = new TH1D("fHistQAvz", "Histo Vz Selected", 5000, -25., 25.);  
  fHistQA[3] = new TH1D("fHistQAvxA", "Histo Vx", 5000, -5., 5.);
  fHistQA[4] = new TH1D("fHistQAvyA", "Histo Vy", 5000, -5., 5.);
  fHistQA[5] = new TH1D("fHistQAvzA", "Histo Vz", 5000, -25., 25.);

  fHistQA[6] = new TH1D("fHistQADcaXyA", "Histo DCAxy", 600, -15., 15.);
  fHistQA[7] = new TH1D("fHistQADcaZA", "Histo DCAz ", 600, -15., 15.);   
  fHistQA[8] = new TH1D("fHistQAPtA","p_{T} distribution",1000,0,10);
  fHistQA[9] = new TH1D("fHistQAEtaA","#eta distribution",240,-1.2,1.2);

  fHistQA[10] = new TH1D("fHistQADcaXy", "Histo DCAxy after Selected", 600, -15., 15.);
  fHistQA[11] = new TH1D("fHistQADcaZ", "Histo DCAz Selected", 600, -15., 15.);   
  fHistQA[12] = new TH1D("fHistQAPt","p_{T} distribution Selected",1000,0,10);
  fHistQA[13] = new TH1D("fHistQAEta","#eta distribution Selected",240,-1.2,1.2);
  
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

  if (!AcceptEvent(event)) return;
  if (isQA) fEventCounter->Fill(2);  
  
  Int_t icharge = -1;
  Int_t gCharge[2];
  Int_t gPid[3][2];
  
  for (Int_t i = 0; i < 2; i++) {
    gCharge[i] = 0;
    for (Int_t ii = 0; ii < 3; ii++) {
      gPid[ii][i] = 0;
    }
  }

  Float_t gCent   = -1;
  Float_t gRefMul = -1;

  if(fAnalysisType == "AOD") {
    AliAODHeader *aodHeader = event->GetHeader();
    gCent = aodHeader->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
    gRefMul = aodHeader->GetRefMultiplicity();
    
    if (gCent < 0 || gCent > 100) return;
    if (isQA) fEventCounter->Fill(50+gCent);
    
    for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {
      AliAODTrack* track = dynamic_cast<AliAODTrack *>(event->GetTrack(itrk));
      if (!track) continue;
      
      if (!AcceptTrack(track)) continue;
            
      Int_t a = fHelperPID->GetParticleSpecies((AliVTrack*) track,kTRUE);
      if(a < 0 || a > 2) continue;
      icharge = track->Charge() > 0 ? 0 : 1;
      gCharge[icharge]++;
      gPid[a][icharge]++;
    }
  }
  else if(fAnalysisType == "MCAOD") {
    TClonesArray *arrayMC = (TClonesArray*) event->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if (!arrayMC) {
      AliFatal("Error: MC particles branch not found!\n");
    }
    AliAODMCHeader *mcHdr=0;
    mcHdr=(AliAODMCHeader*)event->GetList()->FindObject(AliAODMCHeader::StdBranchName());  
    if(!mcHdr) {
      printf("MC header branch not found!\n");
      return;
    }
    
    AliAODHeader *aodHeader = event->GetHeader();
    AliCentrality* fcentrality = aodHeader->GetCentralityP();
    gCent = fcentrality ->GetCentralityPercentile(fCentralityEstimator.Data());
    if (gCent < 0 || gCent > 100) return;
    if (isQA) fEventCounter->Fill(50+gCent);      

    Int_t nMC = arrayMC->GetEntries();
    for (Int_t iMC = 0; iMC < nMC; iMC++) {
      AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(iMC);
      if(!AcceptMCTrack(partMC)) continue;
      Int_t a = fHelperPID->GetMCParticleSpecie((AliVEvent*) event,(AliVTrack*)partMC,1);
      if(a < 0 || a > 2) continue;
      icharge = partMC->Charge() > 0 ? 0 : 1;
      gCharge[icharge]++;
      gPid[a][icharge]++;
    }
  }
  else return;
  
    
  if( (gCharge[0] + gCharge[1]) != 0 ) {
    if (isQA) {
      fEventCounter->Fill(3); 
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
Bool_t AliEbyEParticleRatioFluctuationTask::AcceptEvent(AliAODEvent *event) const {
  Bool_t ver = kFALSE;
  const AliAODVertex *vertex = event->GetPrimaryVertex();
  if(vertex) {
    Double32_t fCov[6];
    vertex->GetCovarianceMatrix(fCov);
    if(vertex->GetNContributors() > 0) {
      if(fCov[5] != 0) {
	
	if(isQA) {	
	  fEventCounter->Fill(5);
	  fHistQA[3]->Fill(vertex->GetX());
	  fHistQA[4]->Fill(vertex->GetY());
	  fHistQA[5]->Fill(vertex->GetZ());
	}
	
	if(TMath::Abs(vertex->GetX()) < fVxMax) {
	  if(TMath::Abs(vertex->GetY()) < fVyMax) {
	    if(TMath::Abs(vertex->GetZ()) < fVzMax) {
	      ver = kTRUE;
	      if(isQA) {	
		fEventCounter->Fill(6);
		fHistQA[0]->Fill(vertex->GetX());
		fHistQA[1]->Fill(vertex->GetY());
		fHistQA[2]->Fill(vertex->GetZ());
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
Bool_t AliEbyEParticleRatioFluctuationTask::AcceptTrack(AliAODTrack *track) const {
  if(!track) return kFALSE;
  if (track->Charge() == 0 ) return kFALSE;
  
  if(isQA) {
    fHistQA[6]->Fill(track->DCA());
    fHistQA[7]->Fill(track->ZAtDCA());
    fHistQA[8]->Fill(track->Pt());
    fHistQA[9]->Fill(track->Eta());
  }
  
  if (!track->TestFilterBit(fAODtrackCutBit)) return kFALSE;

  if (track->Eta() < fEtaLowerLimit ||
      track->Eta() > fEtaHigherLimit) return kFALSE;
  if (track->Pt() < fPtLowerLimit ||
      track->Pt() > fPtHigherLimit) return kFALSE;  
  if ( track->DCA() > fDCAxy ) return kFALSE; 
  if ( track->ZAtDCA() > fDCAz ) return kFALSE;
   
  if(isQA) {
    fHistQA[10]->Fill(track->DCA());
    fHistQA[11]->Fill(track->ZAtDCA());
    fHistQA[12]->Fill(track->Pt());
    fHistQA[13]->Fill(track->Eta());
  }
 
  return kTRUE;
}


//___________________________________________________________
Bool_t AliEbyEParticleRatioFluctuationTask::AcceptMCTrack(AliAODMCParticle *track) const {
  if(!track) return kFALSE;
  if(!track->IsPhysicalPrimary()) return kFALSE;
  if (track->Charge() == 0 ) return kFALSE;
  if(isQA) {
    fHistQA[8]->Fill(track->Pt());
    fHistQA[9]->Fill(track->Eta());
  }

  if (track->Eta() < fEtaLowerLimit ||
      track->Eta() > fEtaHigherLimit) return kFALSE;
  if (track->Pt() < fPtLowerLimit ||
      track->Pt() > fPtHigherLimit) return kFALSE;  
    
  if(isQA) {
    fHistQA[12]->Fill(track->Pt());
    fHistQA[13]->Fill(track->Eta());
  }
 
  return kTRUE;
}
