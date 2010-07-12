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

/* $Id: AliTRDTriggerL1.cxx 31904 2009-04-08 16:42:03Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TRD trigger L1 (GTU) simulation steering                                  //
// currently the Trigger() method calls the GTU tracking simulation and      //
// runs two example triggers, namely on a single high pt particle and        //
// on a jet.                                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TObjArray.h"

#include "AliLog.h"
#include "AliTriggerInput.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliTRDTriggerL1.h"
#include "AliTRDgtuSim.h"
#include "AliTRDtrackGTU.h"

AliTRDTriggerL1::AliTRDTriggerL1()
{
  SetName("TRD");
}

AliTRDTriggerL1::~AliTRDTriggerL1()
{

}

void AliTRDTriggerL1::CreateInputs()
{
  // create the trigger inputs for TRD

  if (fInputs.GetEntriesFast() > 0)
    return;

  fInputs.AddLast(new AliTriggerInput("1HSH", "TRD", 1));
  fInputs.AddLast(new AliTriggerInput("1HJT", "TRD", 1));
}

void AliTRDTriggerL1::Trigger()
{
  // run the trigger algorithms

  AliRunLoader *runLoader = AliRunLoader::Instance();
  if (!runLoader)
    return;
  AliLoader *trdLoader = runLoader->GetLoader("TRDLoader");
  if (!trdLoader)
    return;
  
  // now running the GTU tracking;
  AliTRDgtuSim *gtusim = new AliTRDgtuSim();
  gtusim->RunGTU(trdLoader, 0x0);
  gtusim->WriteTracksToLoader();
  
  TTree *trackTree = trdLoader->GetDataLoader("gtutracks")->Tree();
  if (!trackTree) {
    AliDebug(1,"Did not find track tree");
    return;
  }
  TBranch *branch = trackTree->GetBranch("TRDtrackGTU");
  AliDebug(1,Form("TRD trigger: found %d tracks", trackTree->GetEntriesFast()));
  
  // trigger thresholds should go elsewhere
  Float_t ptThreshold1 = 2;
  Float_t ptThreshold2 = 9.9;
  Int_t trackThreshold1 = 6;
  Int_t trackThreshold2 = 2;
  
  // trigger algorithms to come, e.g.
  Bool_t triggeredHighPt = kFALSE;
  Bool_t triggeredJet = kFALSE;
  
  if (branch) {
    AliTRDtrackGTU *trk = 0x0;
    branch->SetAddress(&trk);

    // high pt trigger
    for (Int_t iTrack = 0; iTrack < trackTree->GetEntriesFast(); iTrack++) {
      trackTree->GetEntry(iTrack);
      if (TMath::Abs(trk->GetPt()) > 3.0) {
        AliDebug(1, Form("Found track in sector %2i, stack %i with pt = %3.1f, triggered", 
                         trk->GetSector(), trk->GetStack(), trk->GetPt()));
        triggeredHighPt = kTRUE;
      }
    }

    // jet trigger
    Int_t nTracks1[90]; // tracks above lower pt threshold
    Int_t nTracks2[90]; // tracks above higher pt threshold
    for (Int_t iTrack = 0; iTrack < trackTree->GetEntriesFast(); iTrack++) {
      trackTree->GetEntry(iTrack);
      if (TMath::Abs(trk->GetPt()) > ptThreshold1)
        nTracks1[5*trk->GetSector() + trk->GetStack()]++;
      if (TMath::Abs(trk->GetPt()) > ptThreshold2)
        nTracks2[5*trk->GetSector() + trk->GetStack()]++;
    }
    for (Int_t iStack = 0; iStack < 90; iStack++) {
      if ((nTracks1[iStack] >= trackThreshold1) || (nTracks2[iStack] >= trackThreshold2))
        triggeredJet = kTRUE;
    }
  }
  else {
    AliWarning("GTU Branch not found");
  }

  if (triggeredHighPt) { 
    AliInfo("Fired high-pt trigger");
    SetInput("1HSH");
  }

  if (triggeredJet) {
    AliInfo("Fired jet trigger");
    SetInput("1HJT");
  }

  // cleaning up
  delete gtusim;
}
