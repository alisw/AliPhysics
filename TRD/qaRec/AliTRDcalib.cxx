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

/* $Id: AliTRDcalib.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Calibration                                                           //
//                                                                        //
//  Authors:                                                              //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TProfile2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <TParticle.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"

#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliTRDcluster.h"
#include "AliTRDgeometry.h"
#include "AliTRDCalibraFillHisto.h"
#include "AliCDBManager.h"
#include "AliTRDcalibDB.h"
#include "TTreeStream.h"

#include <cstdio>
#include <cstring>

#include "AliTRDcalib.h"
#include "AliTRDtrackInfo/AliTRDtrackInfo.h"

ClassImp(AliTRDcalib)


//____________________________________________________________________
AliTRDcalib::AliTRDcalib(const Char_t *name):
  AliAnalysisTask(name, "")
  ,fESD(0x0)
  ,fESDfriend(0x0)
  ,fListHist(0x0)
  ,fTRDCalibraFillHisto(0x0)
  ,flow(0)
  ,fhigh(30)
  ,ffillZero(kFALSE)
  ,fdebugLevel(0)
  ,fspecificstorage("local://$ALICE_ROOT")
{
  //
  // Default constructor
  //

  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());
}
//____________________________________________________________________
void AliTRDcalib::ConnectInputData(Option_t *)
{
  //
  // Link the Input Data
  //
  TTree *tree = dynamic_cast<TChain*>(GetInputData(0));
  if(!tree){
    printf("ERROR - ESD event not found");
  } else {
    tree->SetBranchStatus("Tracks", 1);
    tree->SetBranchStatus("ESDfriend*",1);
  }
  
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!esdH){
		printf("ERROR - ESD input handler not found");
  } else {
    fESD = esdH->GetEvent();
    if(!fESD){
      printf("ERROR - ESD event not found");
    } else {
      esdH->SetActiveBranches("ESDfriend*");
      fESDfriend = (AliESDfriend *)fESD->FindListObject("AliESDfriend");
      printf("fESDfriend = %p\n", (void*)fESDfriend);
    }
  }
}

//____________________________________________________________________
void AliTRDcalib::CreateOutputObjects()
{	
  //
  // Create Output Containers (TObjectArray containing 1D histograms)
  //

  fListHist = new TList();

  // Number of time bins
  AliCDBManager *cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("local://$ALICE_ROOT");
  cdbManager->SetSpecificStorage("TRD/Calib/FEE",fspecificstorage);
  cdbManager->SetRun(0);
  
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  Int_t numberoftimebins = cal->GetNumberOfTimeBins();

  // instance calibration
  fTRDCalibraFillHisto = AliTRDCalibraFillHisto::Instance();
  fTRDCalibraFillHisto->SetHisto2d(); // choose to use histograms
  fTRDCalibraFillHisto->SetCH2dOn();  // choose to calibrate the gain
  fTRDCalibraFillHisto->SetPH2dOn();  // choose to calibrate the drift velocity
  fTRDCalibraFillHisto->Init2Dhistos(); // initialise the histos
  fTRDCalibraFillHisto->SetNumberClusters(flow); // At least flow clusters
  fTRDCalibraFillHisto->SetNumberClustersf(fhigh); // Maximum fhigh clusters
  fTRDCalibraFillHisto->SetFillWithZero(ffillZero); // Fill with Zero or not
  fTRDCalibraFillHisto->SetDebugLevel(fdebugLevel); //debug stuff

  TH1F *nbClusters = new TH1F("NbClusters","",35,0,35);
  nbClusters->Sumw2();
  //
  TProfile2D *pHSum = new TProfile2D("PH2dSum","Nz0Nrphi0"
				     ,numberoftimebins,-0.05,((Double_t)(numberoftimebins/10.0-0.05))
				     ,540,0,540);
  pHSum->SetYTitle("Det/pad groups");
  pHSum->SetXTitle("time [#mus]");
  pHSum->SetZTitle("<PH> [a.u.]");
  pHSum->SetStats(0);
  //
  TH2I *cHSum = new TH2I("CH2dSum","Nz0Nrphi0",100,0,300,540,0,540);
  cHSum->SetYTitle("Det/pad groups");
  cHSum->SetXTitle("charge deposit [a.u]");
  cHSum->SetZTitle("counts");
  cHSum->SetStats(0);
  cHSum->Sumw2();

  
  fListHist->Add(fTRDCalibraFillHisto->GetCH2d()); //TH2I
  fListHist->Add(fTRDCalibraFillHisto->GetPH2d()); //TProfile2D
  fListHist->Add(nbClusters);
  fListHist->Add(pHSum);
  fListHist->Add(cHSum);
      
}
//____________________________________________________________________
void AliTRDcalib::Exec(Option_t *){
  //
  // Run the Analysis
  //
  if(!fESD){
    puts("Error: ESD not found");
    return;
  }
  if(!fESDfriend){
    puts("Error: ESD friend not found");
    return;
  }
  fESD->SetESDfriend(fESDfriend);


  TH1F * nbClusters = (TH1F *) fListHist->At(2);
  if(!nbClusters) {
    puts("Error: nbClusters not found");
    return;
  }
   TProfile2D * pHSum = (TProfile2D *) fListHist->At(3);
   if(!pHSum) {
     puts("Error: pHSum not found");
     return;
   }
   TH2I * cHSum = (TH2I *) fListHist->At(4);
   if(!cHSum) {
     puts("Error: cHSum not found");
     return;
   }
   
   AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
   Int_t numberoftimebins = cal->GetNumberOfTimeBins();
   
   Int_t nTracks = fESD->GetNumberOfTracks();
   
   AliESDtrack *esdTrack = 0x0;
  AliESDfriendTrack *esdFriendTrack = 0x0;
  TObject *calObject = 0x0;
  AliTRDtrackV1 *track = 0x0;
  AliTRDseedV1  *tracklet = 0x0;
  AliTRDcluster *cl = 0x0;
  
  for(Int_t itrk = 0; itrk < nTracks; itrk++){
    
    esdTrack = fESD->GetTrack(itrk);
    
    // read REC info
    esdFriendTrack = fESDfriend->GetTrack(itrk);
    if(esdFriendTrack){
      Int_t icalib = 0;
      while((calObject = esdFriendTrack->GetCalibObject(icalib++))){
	if(strcmp(calObject->IsA()->GetName(),"AliTRDtrackV1") != 0) continue; // Look for the TRDtrack
	track = dynamic_cast<AliTRDtrackV1*>(calObject);
	if(!track) continue;
	fTRDCalibraFillHisto->UpdateHistogramsV1(track);
	for(Int_t itr = 0; itr < 6; itr++){
	  if(!(tracklet = track->GetTracklet(itr))) continue;
	   if(!tracklet->IsOK()) continue;
	   Int_t nbclusters = 0;
	   // For PH
	   Double_t *phtb = new Double_t[numberoftimebins];
	   for(Int_t k=0; k < numberoftimebins; k++){
	     phtb[k] = 0.0;
	   }
	   // For CH
	   Double_t sum = 0.0;
	   // normalisation
	   Float_t normalisation = 6.67;
	   Int_t detector = 0;
	   for(int ic=0; ic<AliTRDseed::knTimebins; ic++){
	     if(!(cl = tracklet->GetClusters(ic))) continue;
	     nbclusters++;
	     Int_t time = cl->GetPadTime();
	     Float_t ch =  TMath::Abs(tracklet->GetdQdl(ic));
	     detector = cl->GetDetector();	  
	     if((time>-1) && (time<numberoftimebins)) phtb[time]=ch/normalisation;
	     sum += ch/normalisation;
	   }
	   nbClusters->Fill(nbclusters);
	   if((nbclusters > flow) && (nbclusters < fhigh)){
	     cHSum->Fill(sum/20.0,0.0);
	     for(int ic=0; ic<numberoftimebins; ic++){
	       if(ffillZero) pHSum->Fill((Double_t)ic/10.0,0.0,(Double_t)phtb[ic]);
	       else {
		 if(phtb[ic] > 0.0) pHSum->Fill((Double_t)ic/10.0,0.0,(Double_t)phtb[ic]);
	       }
	     }
	   }
	}
      }
      
    }
  }
  
  PostData(0, fListHist);
}
//____________________________________________________________________
void AliTRDcalib::Terminate(Option_t *)
{
  //
  // Stays empty because we are only interested in the tree
  //
  printf("terminate\n");
  if(fTRDCalibraFillHisto) fTRDCalibraFillHisto->DestroyDebugStreamer();
}
