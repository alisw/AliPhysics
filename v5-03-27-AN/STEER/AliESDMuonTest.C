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

/* $Id$ */

//
// Macro for checking aliroot output and associated files contents
// Gines Martinez, Subatech June 2003
//

// ROOT includes
#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1.h"
#include "TParticle.h"
#include "TTree.h"

// STEER includes
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliLoader.h"
#include "AliStack.h"
#include "AliESDEvent.h"

// MUON includes
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONConstants.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrack.h"
#include "AliMUONTriggerTrack.h"
#include "AliESDMuonTrack.h"


void AliESDMuonTest(char * filename="galice.root", Int_t run=0){

  TClonesArray* recTracksArray;
  TClonesArray* recTrigTracksArray;
  TTree* treeE;
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  
  // creating ESD file
  TFile *ef=TFile::Open("AliESD.Muon.root","RECREATE");
  if (!ef->IsOpen()) cerr<<"Can't open AliESD.root file!\n";

  treeE = new TTree("ESD", "ESD");

  AliLoader* MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadTracks("READ");

  // Creating MUON data container
  AliMUONData muondata(MUONLoader,"MUON","MUON");

  // declaration  
  Int_t ievent, nevents;
  Int_t ntrackhits;
  Double_t fitfmin;
  Int_t nrectracks;

  Double_t bendingSlope, nonBendingSlope, fInverseBendingMomentum;
  Double_t fXRec, fYRec, fZRec, chi2MatchTrigger;
  Bool_t matchTrigger;

  nevents = RunLoader->GetNumberOfEvents();

  // setting pointer for tracks, triggertracks& trackparam at vertex
  AliMUONTrack* rectrack;
  AliMUONTriggerTrack* rectriggertrack;
  AliMUONTrackParam* trackParam;

  for (ievent = 0; ievent < nevents; ievent++) {

    RunLoader->GetEvent(ievent);

    // setting ESD class pointer
    AliESDEvent* event = new AliESDEvent();
    event->CreateStdContent();
    char name[255];
    event->WriteToTree(treeE);

    event->SetRunNumber(run);
    event->SetEventNumber(ievent);

    // setting ESD MUON class
    AliESDMuonTrack* ESDTrack = new  AliESDMuonTrack() ;

 // ---------------- tracks ----------------
    muondata.SetTreeAddress("RT");
    muondata.GetRecTracks();
    recTracksArray = muondata.RecTracks();
    nrectracks = (Int_t) recTracksArray->GetEntriesFast();

// --------------- trigger tracks ----------
    Long_t trigPat = 0;

    muondata.SetTreeAddress("RL");
    muondata.GetRecTriggerTracks();
    recTrigTracksArray = muondata.RecTriggerTracks();
    rectriggertrack = (AliMUONTriggerTrack*) recTrigTracksArray->First();
    trigPat = rectriggertrack->GetGTPattern();

    printf(">>> Event %d Number of Recconstructed tracks %d \n",ievent, nrectracks);
 
    for (Int_t irectracks = 0; irectracks <  nrectracks;  irectracks++) {

      // ---------------- tracks ----------------
	rectrack = (AliMUONTrack*) recTracksArray->At(irectracks);    
	trackParam = rectrack->GetTrackParamAtVertex();

	bendingSlope            = trackParam->GetBendingSlope();
	nonBendingSlope         = trackParam->GetNonBendingSlope();
	//	printf(" SlopeX %f SlopeY %f\n",bendingSlope ,nonBendingSlope);

	fInverseBendingMomentum = trackParam->GetInverseBendingMomentum();
	fXRec  = trackParam->GetNonBendingCoor();
	fYRec  = trackParam->GetBendingCoor();
	//      printf(" X %f Y %f\n", fXRec, fYRec);

	fZRec  = trackParam->GetZ();
	//      printf(" Z %f\n", fZRec);

	ntrackhits = rectrack->GetNTrackHits();
	fitfmin = rectrack->GetFitFMin();
	matchTrigger     = rectrack->GetMatchTrigger();
	chi2MatchTrigger = rectrack->GetChi2MatchTrigger();

	// setting data member of ESD MUON
	ESDTrack->SetInverseBendingMomentum(fInverseBendingMomentum);
	ESDTrack->SetThetaX(TMath::ATan(nonBendingSlope));
	ESDTrack->SetThetaY(TMath::ATan(bendingSlope));
	ESDTrack->SetZ(fZRec);
	ESDTrack->SetBendingCoor(fYRec);
	ESDTrack->SetNonBendingCoor(fXRec);
	ESDTrack->SetChi2(fitfmin);
	ESDTrack->SetNHit(ntrackhits);
	ESDTrack->SetMatchTrigger(matchTrigger);
	ESDTrack->SetChi2MatchTrigger(chi2MatchTrigger);

      // storing ESD MUON Track into ESD Event & reset muondata
      event->AddMuonTrack(ESDTrack);
    }// track loop
 
    event->SetTrigger(trigPat);

    for (Int_t iTrack = 0; iTrack < event->GetNumberOfMuonTracks(); iTrack++) {
      AliESDMuonTrack* muonTrack = event->GetMuonTrack(iTrack);
      Double_t ptInv = TMath::Abs(muonTrack->GetInverseBendingMomentum());
      cout << "  ptInv: "<<ptInv <<"  nb track: "<< event->GetNumberOfMuonTracks() << endl;
    }
    treeE->Fill();
    event->Reset();
        
    muondata.ResetRecTracks();
    muondata.ResetRecTriggerTracks();

  } // end loop on event
  treeE->GetUserInfo()->Add(event);
  ef->Write();
  ef->Close();
  MUONLoader->UnloadTracks();
  delete event;

}










