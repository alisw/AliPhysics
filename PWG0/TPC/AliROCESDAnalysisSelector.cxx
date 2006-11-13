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
//  This class analyses TPC cosmics data from the ESD and the ESDfriend
//
//  Authors: Jan.Fiete.Grosse-Oetringhaus@cern.ch, Claus.Jorgensen@cern.ch
//

#include "AliROCESDAnalysisSelector.h"

#include <AliLog.h>
#include <AliESD.h>
#include <AliESDfriend.h>
#include <../TPC/AliTPCclusterMI.h>
#include <../TPC/AliTPCseed.h>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>

#include "TPC/AliTPCClusterHistograms.h"

ClassImp(AliROCESDAnalysisSelector)

AliROCESDAnalysisSelector::AliROCESDAnalysisSelector() :
  AliSelector(),
  fESDfriend(0)
{
  //
  // Constructor. Initialization of pointers
  //
  
  for (Int_t i=0; i<kTPCHists; i++)
    fClusterHistograms[i] = 0;
}

AliROCESDAnalysisSelector::~AliROCESDAnalysisSelector()
{
  //
  // Destructor
  //
}

void AliROCESDAnalysisSelector::SlaveBegin(TTree* tree)
{
  //
  
  AliSelector::SlaveBegin(tree);
} 

void AliROCESDAnalysisSelector::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.

  AliSelector::Init(tree);
  
  printf("Init called %p\n", (void*) fESDfriend);

  // Set branch address
  if (tree) 
  {
    tree->SetBranchAddress("ESDfriend", &fESDfriend);
  
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("fTracks.*", 1);
    tree->SetBranchStatus("fTimeStamp", 1);
    //tree->SetBranchStatus("fTracks.fCalibContainer", 0);
  }

  if (fESDfriend != 0)
    AliDebug(AliLog::kInfo, "INFO: Found ESDfriend branch in chain.");
}

Bool_t AliROCESDAnalysisSelector::Process(Long64_t entry)
{
  //
  // Implement your analysis here. Do not forget to call the parent class Process by
  // if (AliSelector::Process(entry) == kFALSE)
  //   return kFALSE;
  //

  if (AliSelector::Process(entry) == kFALSE)
    return kFALSE;

  // Check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  // Check prerequisites
  if (!fESDfriend)
  {
    AliDebug(AliLog::kError, "ESDfriend branch not available");
    return kFALSE;
  }
  
  fESD->SetESDfriend(fESDfriend);

  Int_t nTracks = fESD->GetNumberOfTracks();
  
  Int_t nSkippedSeeds = 0;
  
  // loop over esd tracks
  for (Int_t t=0; t<nTracks; t++)
  {

    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (fESD->GetTrack(t));
    if (!esdTrack)
    {
      AliDebug(AliLog::kError, Form("ERROR: Could not retrieve track %d.", t));
      continue;
    }
    
    AliESDfriendTrack* friendtrack = const_cast<AliESDfriendTrack*> (dynamic_cast<const AliESDfriendTrack*> (esdTrack->GetFriendTrack()));
    if (!friendtrack)
    {
      AliDebug(AliLog::kError, Form("ERROR: Could not retrieve friend of track %d.", t));
      continue;
    }
    
    const AliTPCseed* seed = dynamic_cast<const AliTPCseed*> (friendtrack->GetCalibObject(0));
    if (!seed)
    {
      AliDebug(AliLog::kDebug, Form("ERROR: Could not retrieve seed of track %d.", t));
      nSkippedSeeds++;
      continue;
    }
    
    if (!AcceptTrack(seed)) 
      continue;

    for (Int_t clusterID = 0; clusterID < 160; clusterID++)
    {
      AliTPCclusterMI* cluster = seed->GetClusterPointer(clusterID);
      if (!cluster)
      {
        //AliDebug(AliLog::kError, Form("ERROR: Could not retrieve cluster %d of track %d.", clusterID, t));
        continue;
      }
      
      //AliDebug(AliLog::kDebug, Form("We found a cluster from sector %d", cluster->GetDetector()));

      Int_t detector = cluster->GetDetector();
      
      if (detector < 0 || detector >= kTPCSectors) {
	AliDebug(AliLog::kDebug, Form("We found a cluster from invalid sector %d", detector));
	continue;
      }
      
      // TODO: find a clever way to handle the time      
      Int_t time = 0;

      if (fESD->GetTimeStamp()>1160000000)
	time = fESD->GetTimeStamp();      

      if (!fClusterHistograms[detector])
        fClusterHistograms[detector] = new AliTPCClusterHistograms(detector,"",time,time+7*60*60);
      
      if (!fClusterHistograms[detector+kTPCSectors])
        fClusterHistograms[detector+kTPCSectors] = new AliTPCClusterHistograms(detector,"",time,time+7*60*60, kTRUE);

      fClusterHistograms[detector]->FillCluster(cluster, time);
      fClusterHistograms[detector+kTPCSectors]->FillCluster(cluster, time);
    }
  }
  
  if (nSkippedSeeds > 0)
    printf("WARNING: The seed was not found for %d out of %d tracks.\n", nSkippedSeeds, nTracks);

  // TODO This should not be needed, the TTree::GetEntry() should take care of this, maybe because it has a reference member, to be analyzed
  // if the ESDfriend is not deleted we get a major memory leak
  // here the esdfriend seems to be also deleted, very weird behaviour....
  delete fESD;
  fESD = 0;    
  
  //delete fESDfriend;
  //fESDfriend = 0;
   
  return kTRUE;
}


Bool_t AliROCESDAnalysisSelector::AcceptTrack(const AliTPCseed* track) {
  //
  //
  //

  // TODO : implement min number of rows to accept track

  const Int_t   kMinClusters = 20;
  const Float_t kMinRatio    = 0.75;
  const Float_t kMax1pt      = 0.5;


  if (track->GetNumberOfClusters()<kMinClusters) return kFALSE;
  Float_t ratio = track->GetNumberOfClusters()/(track->GetNFoundable()+1.);
  if (ratio<kMinRatio) return kFALSE;
  Float_t mpt = track->Get1Pt();
  if (TMath::Abs(mpt)>kMax1pt) return kFALSE;

  //if (TMath::Abs(track->GetZ())>240.) return kFALSE;
  //if (TMath::Abs(track->GetZ())<10.) return kFALSE;
  //if (TMath::Abs(track->GetTgl())>0.03) return kFALSE;
  
  return kTRUE;
}

void AliROCESDAnalysisSelector::SlaveTerminate()
{
  //
  
  if (fOutput)
  {
    for (Int_t i=0; i<kTPCHists; i++)
      if (fClusterHistograms[i])
        fOutput->Add(fClusterHistograms[i]);
  }
} 

void AliROCESDAnalysisSelector::Terminate()
{
  // 
  // read the objects from the output list and write them to a file
  // the filename is modified by the object comment passed in the tree info or input list
  //

  if (fOutput)
  {  
    fOutput->Print();
        
    for (Int_t i=0; i<kTPCSectors; i++)
      fClusterHistograms[i] = dynamic_cast<AliTPCClusterHistograms*> (fOutput->FindObject(AliTPCClusterHistograms::FormDetectorName(i, kFALSE)));
    for (Int_t i=0; i<kTPCSectors; i++)
      fClusterHistograms[kTPCSectors+i] = dynamic_cast<AliTPCClusterHistograms*> (fOutput->FindObject(AliTPCClusterHistograms::FormDetectorName(i, kTRUE)));
  }
  
  TNamed* comment = 0;
  if (fTree && fTree->GetUserInfo())
    comment = dynamic_cast<TNamed*>(fTree->GetUserInfo()->FindObject("comment"));
  if (!comment && fInput)
    comment = dynamic_cast<TNamed*>(fInput->FindObject("comment"));

  if (comment)
  {
    AliDebug(AliLog::kInfo, Form("INFO: Found comment in input list: %s \n", comment->GetTitle()));
  }
  else
    return;

  TFile* file = TFile::Open(Form("rocESD_%s.root",comment->GetTitle()), "RECREATE");
  
  for (Int_t i=0; i<kTPCHists; i++)
    if (fClusterHistograms[i]) {
      fClusterHistograms[i]->SaveHistograms();
      TCanvas* c = fClusterHistograms[i]->DrawHistograms();
      c->SaveAs(Form("plots_%s_%s.eps",comment->GetTitle(),c->GetName()));
      c->SaveAs(Form("plots_%s_%s.gif",comment->GetTitle(),c->GetName()));

      c->Close();
      delete c;
    }
  file->Close();
} 
