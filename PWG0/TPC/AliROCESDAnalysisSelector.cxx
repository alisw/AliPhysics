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

// The ESD is available as member fESD
//
// The Process function is nearly empty. Implement your analysis there and look at the other listed below functions you
// might need.
//
// The following methods can be overrriden. Please do not forgot to call the base class function.
//
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Init():         called for each new tree. Enable/Disable branches here.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
//  Author: Jan.Fiete.Grosse-Oetringhaus@cern.ch

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
  fESDfriend(0),
  fClusterHistograms(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AliROCESDAnalysisSelector::~AliROCESDAnalysisSelector()
{
  //
  // Destructor
  //

  if (fClusterHistograms) {
    delete fClusterHistograms;
    fClusterHistograms = 0;    
  }

}

void AliROCESDAnalysisSelector::SlaveBegin(TTree* tree)
{
  //
  
  AliSelector::SlaveBegin(tree);

  fClusterHistograms = new AliTPCClusterHistograms("test","test");

} 

void AliROCESDAnalysisSelector::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.

  AliSelector::Init(tree);

  // Set branch address
  if (tree)
    tree->SetBranchAddress("ESDfriend", &fESDfriend);
    
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
      AliDebug(AliLog::kError, Form("ERROR: Could not retrieve seed of track %d.", t));
      continue;
    }
    
    for (Int_t clusterID = 0; clusterID < 160; clusterID++)
    {
      AliTPCclusterMI* cluster = seed->GetClusterPointer(clusterID);
      if (!cluster)
      {
        //AliDebug(AliLog::kError, Form("ERROR: Could not retrieve cluster %d of track %d.", clusterID, t));
        continue;
      }
      
      //AliDebug(AliLog::kDebug, Form("We found a cluster from sector %d", cluster->GetDetector()));

      if (cluster->GetDetector() != 5)
        continue;
      
      fClusterHistograms->FillCluster(cluster);
    }
  }
   
  return kTRUE;
}

void AliROCESDAnalysisSelector::SlaveTerminate()
{
  //
  
  fOutput->Add(fClusterHistograms);
} 

void AliROCESDAnalysisSelector::Terminate()
{
  //
    
  TFile* file = TFile::Open("rocESD.root", "RECREATE");
  
  fClusterHistograms->SaveHistograms();
  file->Close();
} 
