/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: ALICE Offline.                                                 *
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
//          Dealing with Wide pT Addition (Test Only - not for use)
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
#include "AliPIDCombined.h"
#include "AliAODHeader.h"
#include "AliAODpidUtil.h"
#include "AliHelperPID.h"
using std::endl;
using std::cout;
#include "AliEbyEPidTTask.h"

ClassImp(AliEbyEPidTTask)

//-----------------------------------------------------------------------
AliEbyEPidTTask::AliEbyEPidTTask( const char *name ) : AliAnalysisTaskSE( name ), 
  fThnList(NULL), 
  isMC(kFALSE), 
  fCentralityEstimator("V0M"), 
  fAODtrackCutBit(128),
  fHelperPID(0x0),
  fEventCounter(NULL), 
  fEventTree(NULL),
  isQA(kFALSE), 
  fDebug(kFALSE), 
  fRunNumber(0),
  fNumberOfTracks(8000),
  fNumberOfTracksM(8000),
  fFilterBit(99),
  fCentPercentile(-1),
  fVertexX(-99),
  fVertexY(-99),
  fVertexZ(-99) {
  /*
    fTrackPt[kTrack];
    fTrackEta[kTrack];
    fTrackPhi[kTrack];
    fTrackPtM[kTrack];
    fTrackEtaM[kTrack];
    fTrackPhiM[kTrack];
    fTrackCharge[kTrack];
    fTrackPid[kTrack];
    fTrackChargeM[kTrack];
    fTrackPidM[kTrack];
  */
 
  DefineOutput(1, TList::Class()); //! Connect Outpput....
  DefineOutput(2, TTree::Class()); //! Connect Outpput....
}

AliEbyEPidTTask::~AliEbyEPidTTask() {
  //!   Cleaning up
  if (fThnList)   delete fThnList;
  if (fHelperPID) delete fHelperPID;
  if (fEventTree) delete fEventTree;
}

//---------------------------------------------------------------------------------
void AliEbyEPidTTask::UserCreateOutputObjects() {
  fThnList = new TList();
  fThnList->SetOwner(kTRUE);

  fEventCounter = new TH1D("fEventCounter","EventCounter", 70, 0.5,70.5);
  fThnList->Add(fEventCounter);

  TList *ll = (TList*)fHelperPID->GetOutputList();
  for (Int_t ikey = 0; ikey < ll->GetEntries(); ikey++) {
    fThnList->Add(ll->At(ikey));
  }
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fEventTree = new TTree("fEventTree","fEventTree");
  owd->cd();

  fEventTree->Branch("fRunNumber",      &fRunNumber,     "fRunNumber/I");
  fEventTree->Branch("fFilterBit",      &fFilterBit,     "fFilterBit/I");
  fEventTree->Branch("fNumberOfTracks", &fNumberOfTracks,"fNumberOfTracks/I");
  fEventTree->Branch("fCentPercentile", &fCentPercentile,"fCentPercentile/F");
  
  fEventTree->Branch("fVertexX",        &fVertexX,       "fVertexX/F");
  fEventTree->Branch("fVertexY",        &fVertexY,       "fVertexY/F");
  fEventTree->Branch("fVertexZ",        &fVertexZ,       "fVertexZ/F");
  
  fEventTree->Branch("fTrackPt",        fTrackPt,        "fTrackPt[fNumberOfTracks]/F");
  fEventTree->Branch("fTrackPhi",       fTrackPhi,       "fTrackPhi[fNumberOfTracks]/F");
  fEventTree->Branch("fTrackEta",       fTrackEta,       "fTrackEta[fNumberOfTracks]/F");
  fEventTree->Branch("fTrackCharge",    fTrackCharge,    "fTrackCharge[fNumberOfTracks]/I");
  fEventTree->Branch("fTrackPid",       fTrackPid,       "fTrackPid[fNumberOfTracks]/I");
  
  if(isMC) {
    fEventTree->Branch("fNumberOfTracksM", &fNumberOfTracksM, "fNumberOfTracksM/I");
    fEventTree->Branch("fTrackPtM",        fTrackPtM,        "fTrackPtM[fNumberOfTracksM]/F");
    fEventTree->Branch("fTrackPhiM",       fTrackPhiM,       "fTrackPhiM[fNumberOfTracksM]/F");
    fEventTree->Branch("fTrackEtaM",       fTrackEtaM,       "fTrackEtaM[fNumberOfTracksM]/F");
    fEventTree->Branch("fTrackChargeM",    fTrackChargeM,    "fTrackChargeM[fNumberOfTracksM]/I");
    fEventTree->Branch("fTrackPidM",       fTrackPidM,       "fTrackPidM[fNumberOfTracksM]/I");
  }

  PostData(1, fThnList);
  PostData(2, fEventTree);  
}

//----------------------------------------------------------------------------------
void AliEbyEPidTTask::UserExec( Option_t * ){

   fEventCounter->Fill(1);

  AliAODEvent* event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event) {
    Printf("ERROR 01: AOD not found ");
    return;
  }

  fEventCounter->Fill(2);

  Int_t gCent   = -1;
  Float_t gRefMul = -1;
  
  AliAODHeader *aodHeader = event->GetHeader();
  gCent = (Int_t)aodHeader->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
  gRefMul = aodHeader->GetRefMultiplicity();
  if (gCent < 0 || gCent > 100) return;
  fEventCounter->Fill(3);  

  const AliAODVertex *vertex = event->GetPrimaryVertex();
  if(!vertex) return;
  fEventCounter->Fill(4);
  Bool_t vtest = kFALSE;
  Double32_t fCov[6];
  vertex->GetCovarianceMatrix(fCov);
  if(vertex->GetNContributors() > 0) {
      if(fCov[5] != 0) {
	vtest = kTRUE;
      }
  }
  if(!vtest)return;
  
  fEventCounter->Fill(5);

  AliCentrality *centrality = event->GetCentrality();
  if (centrality->GetQuality() != 0) return;

  fEventCounter->Fill(6);

  fRunNumber = event->GetRunNumber();
  fFilterBit = fAODtrackCutBit;
  fCentPercentile = gCent;
  fVertexX = vertex->GetX();
  fVertexY = vertex->GetY();
  fVertexZ = vertex->GetZ();
  
//Default TPC priors
  if(fHelperPID->GetPIDType()==kBayes)fHelperPID->GetPIDCombined()->SetDefaultTPCPriors();//FIXME maybe this can go in the UserCreateOutputObject?




  Int_t iTracks = 0; 
  for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {
    AliAODTrack* track = dynamic_cast<AliAODTrack *>(event->GetTrack(itrk));
    fEventCounter->Fill(10);
    if (!track) continue;
    fEventCounter->Fill(11);
    if (!AcceptTrack(track)) continue;
    fEventCounter->Fill(12);
    Int_t a = fHelperPID->GetParticleSpecies((AliVTrack*) track,kTRUE);
    //  if(a < 0 || a > 2) continue;
    fEventCounter->Fill(13);
    
    Int_t b = -999;
    if (a == 0 ) b = 1;
    else if (a == 1 ) b = 2;
    else if (a == 2 ) b = 3;
    else b = 4;
    
    if (track->Charge()  < 0 ) b = -1*b;
    
    //    Int_t icharge = track->Charge() > 0 ? 0 : 1;

    // cout << icharge << "  " << track->Charge() << endl;

    fTrackPt[iTracks]     = (Float_t)track->Pt();
    fTrackPhi[iTracks]    = (Float_t)track->Phi();
    fTrackEta[iTracks]    = (Float_t)track->Eta();
    fTrackCharge[iTracks] = track->Charge();
    fTrackPid[iTracks] = b;
    iTracks++;
  }
  fNumberOfTracks = iTracks;

  if(isMC) {
    fEventCounter->Fill(21);
    TClonesArray *arrayMC= 0; 
    arrayMC = dynamic_cast<TClonesArray*> (event->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
    if (!arrayMC) {
      Printf("Error: MC particles branch not found!\n");
      return;
    }
    fEventCounter->Fill(22);
    AliAODMCHeader *mcHdr=0;
    mcHdr=(AliAODMCHeader*)event->GetList()->FindObject(AliAODMCHeader::StdBranchName());  
    if(!mcHdr) {
      Printf("MC header branch not found!\n");
      return;
    }
    
    fEventCounter->Fill(23);
    
    Int_t nMC = arrayMC->GetEntries();
    Int_t mTracks = 0; 
    for (Int_t iMC = 0; iMC < nMC; iMC++) {
      fEventCounter->Fill(24);
      AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(iMC);
      if(!AcceptMCTrack(partMC)) continue;
      
      fEventCounter->Fill(25);
      Int_t a  = partMC->GetPdgCode();
      
      // Int_t a = fHelperPID->GetMCParticleSpecie((AliVEvent*) event,(AliVTrack*)partMC,1);
      
      //if(a < 0 || a > 2) continue;
      Int_t icharge = a > 0 ? 0 : 1;
      
      //  cout << a << "  " << icharge << endl;
      
      fTrackPtM[mTracks]     = (Float_t)partMC->Pt();
      fTrackPhiM[mTracks]    = (Float_t)partMC->Phi();
      fTrackEtaM[mTracks]    = (Float_t)partMC->Eta();
      fTrackChargeM[mTracks] = icharge;
      fTrackPidM[mTracks] = a;
      
      mTracks++;
    }
    fEventCounter->Fill(26);
    fNumberOfTracksM = mTracks;
  }
  
  fEventCounter->Fill(30);
  fEventTree->Fill();
  
  PostData(1, fThnList);  
  PostData(2, fEventTree);
}

//___________________________________________________________
void AliEbyEPidTTask::Terminate( Option_t * ){
  Info("AliEbyEPidTTask"," Task Successfully finished");
}

//___________________________________________________________
Bool_t AliEbyEPidTTask::AcceptTrack(AliAODTrack *track) const {
  if(!track)                                  return kFALSE;
  if (track->Charge() == 0 )                  return kFALSE;
  if (!track->TestFilterBit(fAODtrackCutBit)) return kFALSE;
  return kTRUE;
}


//___________________________________________________________
Bool_t AliEbyEPidTTask::AcceptMCTrack(AliAODMCParticle *track) const {
  if(!track) return kFALSE;
  if(!track->IsPhysicalPrimary()) return kFALSE;
  if (track->Charge() == 0 ) return kFALSE;
  return kTRUE;
}
