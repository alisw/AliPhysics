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
#include "AliESD.h"

// MUON includes
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONConstants.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrack.h"
#include "AliMUONTriggerTrack.h"
#include "AliESDMuonTrack.h"


void AliESDMuonTest(char * filename="galice.root", Int_t run=0){

  TClonesArray * recTracksArray;
  
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  
  // creating ESD file
  TFile *ef=TFile::Open("AliESD.Muon.root","RECREATE");
  if (!ef->IsOpen()) cerr<<"Can't open AliESD.root file!\n";

  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadTracks("READ");

  // Creating MUON data container
  AliMUONData muondata(MUONLoader,"MUON","MUON");

  // declaration  
  Int_t ievent, nevents;
  Int_t ntrackhits;
  Double_t fitfmin;
  Int_t nrectracks;

  Double_t bendingSlope, nonBendingSlope, fInverseBendingMomentum;
  Double_t fXRec, fYRec, fZRec;

  Float_t  x11, y11, thetaX,thetaY ;

  nevents = RunLoader->GetNumberOfEvents();
  
  // setting pointer for tracks, triggertracks& trackparam at vertex
  AliMUONTrack * rectrack;
  AliMUONTriggerTrack * rectriggertrack;
  AliMUONTrackParam *trackParam;

  for (ievent = 0; ievent < nevents; ievent++) {
    RunLoader->GetEvent(ievent);

    //   cerr<<"\n\nProcessing event number : "<<ievent<<endl;

    // setting ESD class pointer
    AliESD *event = new AliESD(); 
    event->SetRunNumber(run);
    event->SetEventNumber(ievent);

    // -------------------- tracks-------------

    // setting ESD MUON class
    AliESDMuonTrack* ESDTrack = new  AliESDMuonTrack() ;
    muondata.SetTreeAddress("RT");
    muondata.GetRecTracks();
    recTracksArray = muondata.RecTracks();
        
    nrectracks = (Int_t) recTracksArray->GetEntriesFast(); //
 
    printf(">>> Event %d Number of Recconstructed tracks %d \n",ievent, nrectracks);
   
    // read track infos
    for (Int_t irectracks = 0; irectracks <  nrectracks;  irectracks++) {

      rectrack = (AliMUONTrack*) recTracksArray->At(irectracks);

      trackParam = rectrack->GetTrackParamAtVertex();

      bendingSlope            = trackParam->GetBendingSlope();
      nonBendingSlope         = trackParam->GetNonBendingSlope();
      fInverseBendingMomentum = trackParam->GetInverseBendingMomentum();
      fXRec  = trackParam->GetNonBendingCoor();
      fYRec  = trackParam->GetBendingCoor();
      fZRec  = trackParam->GetZ();

      ntrackhits = rectrack->GetNTrackHits();
      fitfmin = rectrack->GetFitFMin();

      // setting data member of ESD MUON
      ESDTrack->SetInverseBendingMomentum(fInverseBendingMomentum);
      ESDTrack->SetThetaX(TMath::ATan(nonBendingSlope));
      ESDTrack->SetThetaY(TMath::ATan(bendingSlope));
      ESDTrack->SetZ(fZRec);
      ESDTrack->SetBendingCoor(fYRec);
      ESDTrack->SetNonBendingCoor(fXRec);
      ESDTrack->SetChi2(fitfmin);
      ESDTrack->SetNHit(ntrackhits);
    }

    // -------------------- trigger tracks-------------
    muondata.SetTreeAddress("RL");
    muondata.GetRecTriggerTracks();
    recTracksArray = muondata.RecTriggerTracks();
        
    nrectracks = (Int_t) recTracksArray->GetEntriesFast(); //
 
    printf(">>> Event %d Number of Recconstructed tracks %d \n",ievent, nrectracks);
   
    // read trigger track infos
    for (Int_t irectracks = 0; irectracks <  nrectracks;  irectracks++) {

      rectriggertrack = (AliMUONTriggerTrack*) recTracksArray->At(irectracks);
    
      x11 = rectriggertrack->GetY11();
      y11 = rectriggertrack->GetY11();
      thetaX = rectriggertrack->GetThetax();
      thetaY = rectriggertrack->GetThetay();

      // setting data member of ESD MUON trigger
      ESDTrack->SetThetaX11(thetaX);
      ESDTrack->SetThetaY11(thetaY);
      ESDTrack->SetX11(x11);
      ESDTrack->SetY11(y11);
    }

    // storing ESD MUON Track into ESD Event & reset muondata
    event->AddMuonTrack(ESDTrack);
    muondata.ResetRecTracks();
    muondata.ResetRecTriggerTracks();

    // writting ESD event
    Char_t ename[100]; 
    sprintf(ename,"%d",ievent);
    ef->cd();
    if (!event->Write(ename)) cerr<<"Something bad happened...\n";
    delete event;

  } // end loop on event  
  ef->Close();
  MUONLoader->UnloadTracks();
  MUONLoader->UnloadTriggerTracks();

}










