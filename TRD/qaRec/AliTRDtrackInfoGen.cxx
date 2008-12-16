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

/* $Id: AliTRDtrackInfoGen.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
//  Authors:                                                              //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TH1F.h>
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
#include "AliESDHeader.h"
#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliTRDtrackV1.h"
#include "AliTrackReference.h"
#include "AliTRDgeometry.h"
#include "TTreeStream.h"

#include <cstdio>
#include <cstring>
#include <iostream>

#include "AliTRDtrackInfoGen.h"
#include "AliTRDtrackInfo/AliTRDtrackInfo.h"
#include "AliTRDtrackInfo/AliTRDeventInfo.h"

ClassImp(AliTRDtrackInfoGen)

const Float_t AliTRDtrackInfoGen::xTPC = 290.;
const Float_t AliTRDtrackInfoGen::xTOF = 365.;

//____________________________________________________________________
AliTRDtrackInfoGen::AliTRDtrackInfoGen():
  AliTRDrecoTask("InfoGen", "Track List Generator")
  ,fESD(0x0)
  ,fMC(0x0)
  ,fESDfriend(0x0)
  ,fTrackInfo(0x0)
  ,fEventInfo(0x0)
{
  //
  // Default constructor
  //

  DefineInput(0, TChain::Class());
  DefineOutput(0, TObjArray::Class());
  DefineOutput(1, AliTRDeventInfo::Class());
  //DefineOutput(1, TTree::Class());
}

//____________________________________________________________________
AliTRDtrackInfoGen::~AliTRDtrackInfoGen()
{
  if(fTrackInfo) delete fTrackInfo; fTrackInfo = 0x0;
  if(fEventInfo) delete fEventInfo; fEventInfo = 0x0;
  if(fContainer){ 
    fContainer->Delete(); delete fContainer;
    fContainer = 0x0;
  }
}

//____________________________________________________________________
void AliTRDtrackInfoGen::ConnectInputData(Option_t *)
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
    }
  }
  if(HasMCdata()){
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!mcH){ 
      AliError("MC input handler not found");
    } else {
      fMC = mcH->MCEvent();
    }
  }
}

//____________________________________________________________________
void AliTRDtrackInfoGen::CreateOutputObjects()
{	
  //
  // Create Output Containers (TObjectArray containing 1D histograms)
  //
  fTrackInfo = new AliTRDtrackInfo();
  fEventInfo = new AliTRDeventInfo();
  fContainer = new TObjArray(1000);
  fContainer->SetOwner(kTRUE);

/*  OpenFile(1, "RECREATE");
  fTree = new TTree("trd", "extract of the TRD detector");
  fTree->Branch("info",  &fTrackInfo);
  printf("output tree build in %s\n", fTree->GetDirectory()->GetName());*/
}

//____________________________________________________________________
void AliTRDtrackInfoGen::Exec(Option_t *){
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
  if(HasMCdata() && !fMC){
    puts("Error: Monte Carlo Event not available");
    return;
  }
  fContainer->Delete();
  fEventInfo->Delete("");
  fESD->SetESDfriend(fESDfriend);
  new(fEventInfo)AliTRDeventInfo(fESD->GetHeader(), const_cast<AliESDRun *>(fESD->GetESDRun()));
  
  Bool_t *trackMap = 0x0;
  AliStack * mStack = 0x0;
  if(HasMCdata()){
    mStack = fMC->Stack();
    if(!mStack){
      puts("Error: Cannot get the Monte Carlo Stack");
      return;
    }
    trackMap = new Bool_t[fMC->GetNumberOfTracks()];
    memset(trackMap, 0, sizeof(Bool_t) * fMC->GetNumberOfTracks());
  }
  
  Int_t nTRD = 0, nTPC = 0, nclsTrklt;
  Int_t nTracks = fESD->GetNumberOfTracks();
  if(fDebugLevel>=1){ 
    printf("Entry[%3d] Tracks: ESD[%d] MC[%d]\n", (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), nTracks, HasMCdata() ? mStack->GetNtrack() : 0);
  }
  AliESDtrack *esdTrack = 0x0;
  AliESDfriendTrack *esdFriendTrack = 0x0;
  TObject *calObject = 0x0;
  AliTRDtrackV1 *track = 0x0;
  for(Int_t itrk = 0; itrk < nTracks; itrk++){
    esdTrack = fESD->GetTrack(itrk);
    if(fDebugLevel>=2) printf("\n%3d ITS[%d] TPC[%d] TRD[%d]\n", itrk, esdTrack->GetNcls(0), esdTrack->GetNcls(1), esdTrack->GetNcls(2));
    if(esdTrack->GetNcls(1)) nTPC++;
    if(esdTrack->GetNcls(2)) nTRD++;

    // look at external track param
    const AliExternalTrackParam *op = esdTrack->GetOuterParam();
    Double_t xyz[3];
    if(op){
        op->GetXYZ(xyz);
        op->Global2LocalPosition(xyz, op->GetAlpha());
        if(fDebugLevel>=2) printf("op @ X[%7.3f]\n", xyz[0]);
    }

    // read MC info
    Int_t fPdg = -1;
    Int_t label = -1;
    if(HasMCdata()){
      label = esdTrack->GetLabel();
      if(label < fMC->GetNumberOfTracks()) trackMap[TMath::Abs(label)] = kTRUE; // register the track
      //if (TMath::Abs(label) > mStack->GetNtrack()) continue; 
      AliMCParticle *mcParticle = 0x0; 
      if(!(mcParticle = fMC->GetTrack(TMath::Abs(label)))){
        printf("E - AliTRDtrackInfoGen::Exec() : MC particle missing for ESD label %d\n", label);
        continue;
      }
      fPdg = mcParticle->Particle()->GetPdgCode();
      Int_t nRefs = mcParticle->GetNumberOfTrackReferences();
      Int_t iref = 0; AliTrackReference *ref = 0x0; 
      while(iref<nRefs){
        ref = mcParticle->GetTrackReference(iref);
        if(ref->LocalX() > xTPC) break;
        //printf("\ttrackRef[%2d] @ %7.3f\n", iref, ref->LocalX());
        iref++;
      }

      new(fTrackInfo) AliTRDtrackInfo();
      fTrackInfo->SetPDG(fPdg);
      fTrackInfo->SetPrimary(mcParticle->Particle()->IsPrimary());
      Int_t jref = iref;//, kref = 0;
      while(jref<nRefs){
        ref = mcParticle->GetTrackReference(jref);
        if(ref->LocalX() > xTOF) break;
        if(fDebugLevel>=3) printf("\ttrackRef[%2d (%2d)] @ %7.3f OK\n", jref-iref, jref, ref->LocalX());
        fTrackInfo->AddTrackRef(ref);
        jref++;
      }
      if(fDebugLevel>=2) printf("NtrackRefs[%d(%d)]\n", fTrackInfo->GetNTrackRefs(), nRefs);
    } else {
      new (fTrackInfo) AliTRDtrackInfo();
      fTrackInfo->SetPDG(fPdg);
    }

    // copy some relevant info to TRD track info
    fTrackInfo->SetStatus(esdTrack->GetStatus());
    fTrackInfo->SetTrackId(esdTrack->GetID());
    Double_t p[AliPID::kSPECIES]; esdTrack->GetTRDpid(p);
    fTrackInfo->SetESDpid(p);
    fTrackInfo->SetESDpidQuality(esdTrack->GetTRDpidQuality());
    fTrackInfo->SetLabel(label);
    fTrackInfo->SetNumberOfClustersRefit(esdTrack->GetNcls(2));
    // some other Informations which we may wish to store in order to find problematic cases
    fTrackInfo->SetKinkIndex(esdTrack->GetKinkIndex(0));
    fTrackInfo->SetTPCncls(static_cast<UShort_t>(esdTrack->GetNcls(1)));
    nclsTrklt = 0;
  

    // read REC info
    esdFriendTrack = fESDfriend->GetTrack(itrk);
    if(esdFriendTrack){
      Int_t icalib = 0;
      while((calObject = esdFriendTrack->GetCalibObject(icalib++))){
        if(strcmp(calObject->IsA()->GetName(),"AliTRDtrackV1") != 0) continue; // Look for the TRDtrack
        if(!(track = dynamic_cast<AliTRDtrackV1*>(calObject))) break;
        nTRD++;
        if(fDebugLevel>=3) printf("TRD track OK\n");
        fTrackInfo->SetTrack(track);
        break;
      }
      if(fDebugLevel>=2) printf("Ntracklets[%d]\n", fTrackInfo->GetNTracklets());
    } else if(fDebugLevel>=2) printf("No ESD friends\n");
    if(op) fTrackInfo->SetOuterParam(op);

    if(fDebugLevel >= 1){
      AliTRDtrackInfo info(*fTrackInfo);
      (*fDebugStream) << "trackInfo"
      << "TrackInfo.=" << &info
      << "\n";
      info.Delete("");
    }
  
    fContainer->Add(new AliTRDtrackInfo(*fTrackInfo));
    fTrackInfo->Delete("");
  }
  if(fDebugLevel>=2) printf("%3d Tracks: TPC[%d] TRD[%d]\n", (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), nTPC, nTRD);

  // Insert also MC tracks which are passing TRD where the track is not reconstructed
  if(HasMCdata()){
    if(fDebugLevel > 10){
      printf("Output of the MC track map:\n");
      for(Int_t itk = 0; itk < fMC->GetNumberOfTracks();  itk++)
        printf("trackMap[%d] = %s\n", itk, trackMap[itk] == kTRUE ? "TRUE" : "kFALSE");
    }
  
    for(Int_t itk = 0; itk < fMC->GetNumberOfTracks(); itk++){
      if(fDebugLevel >=2 ) printf("Number of MC tracks: %d\n", fMC->GetNumberOfTracks());
      if(trackMap[itk]) continue;
      AliMCParticle *mcParticle = fMC->GetTrack(TMath::Abs(itk));
      Int_t fPdg = mcParticle->Particle()->GetPdgCode();
      Int_t nRefs = mcParticle->GetNumberOfTrackReferences();
      Int_t iref = 0; AliTrackReference *ref = 0x0; 
      Int_t nRefsTRD = 0;
      new(fTrackInfo) AliTRDtrackInfo();
      fTrackInfo->SetPDG(fPdg);
      while(iref<nRefs){
        ref = mcParticle->GetTrackReference(iref);
        if(fDebugLevel > 3) printf("\ttrackRef[%2d] @ %7.3f", iref, ref->LocalX());
        if(ref->LocalX() > 250. && ref->LocalX() < 370.){
          if(fDebugLevel > 3) printf(" OK\n");
          fTrackInfo->AddTrackRef(ref);
          nRefsTRD++;
        }
        else
          if(fDebugLevel > 3) printf("\n");
        iref++;
      }
      if(!nRefsTRD){
        // In this stage we at least require 1 hit inside TRD. What will be done with this tracks is a task for the 
        // analysis job
        fTrackInfo->Delete("");
        continue;
      }
      fTrackInfo->SetPrimary(mcParticle->Particle()->IsPrimary());
      fTrackInfo->SetLabel(itk);
      if(fDebugLevel >= 1){
        AliTRDtrackInfo info(*fTrackInfo);
        (*fDebugStream) << "trackInfo"
        << "TrackInfo.=" << &info
        << "\n";
        info.Delete("");
      }
      if(fDebugLevel > 2)printf("Registering rejected MC track with label %d\n", itk);
      fContainer->Add(new AliTRDtrackInfo(*fTrackInfo));
      fTrackInfo->Delete("");
    }
    delete[] trackMap;
  }
  PostData(0, fContainer);
  PostData(1, fEventInfo);
}


//____________________________________________________________________
void AliTRDtrackInfoGen::Terminate(Option_t *)
{
  //
  // Stays empty because we are only interested in the tree
  //
  if(fDebugLevel>=1) printf("Terminate:\n");
  //TFile *f =((TFile*)gROOT->FindObject("TRD.TrackInfo.root"));
  //f->cd(); f->Write(); f->Close();

  if(fDebugStream){ 
    delete fDebugStream;
    fDebugStream = 0x0;
  }
}
