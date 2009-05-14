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

/* $Id: AliTRDTrigger.cxx 31904 2009-04-08 16:42:03Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TRD trigger interface class to CTP                                        //
// currently the Trigger() method calls the GTU tracking simulation and      //
// runs two example triggers, namely on a single high pt particle and        //
// on a jet.                                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TClonesArray.h"

#include "AliLog.h"
#include "AliTriggerInput.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliTRDgtuSim.h"
#include "AliTRDtrackGTU.h"
#include "AliTRDTrigger.h"

AliTRDTrigger::AliTRDTrigger()
{
  // defautl constructor

  SetName("TRD");
}

AliTRDTrigger::~AliTRDTrigger()
{
  // destructor
}

void AliTRDTrigger::CreateInputs()
{
  // Create the inputs to CTP for the TRD

  // if inputs already created return
  if (fInputs.GetEntriesFast() > 0)
    return;

  AliInfo("Creating TRD trigger inputs");
  fInputs.AddLast(new AliTriggerInput("TRD_HIGHPT_L1", "TRD", 1));
  fInputs.AddLast(new AliTriggerInput("TRD_JET_L1", "TRD", 1));
}

void AliTRDTrigger::Trigger()
{
  // TRD trigger steering
  // currently the L1 trigger is directly put here
  // lateron can be separated such that from here the
  // pretrigger (generating an L0) and L1 can be called

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
  AliDebug(1,Form("TRD trigger: found %i tracks", trackTree->GetEntriesFast()));

  // trigger thresholds should go elsewhere
  Float_t ptThreshold1 = 2;
  Float_t ptThreshold2 = 9.9;
  Int_t trackThreshold1 = 6;
  Int_t trackThreshold2 = 2;

  // trigger algorithms to come, e.g.
  Bool_t triggered_highpt = kFALSE;
  Bool_t triggered_jet = kFALSE;

  if (branch) {
    AliTRDtrackGTU *trk = 0x0;
    branch->SetAddress(&trk);

    // high pt trigger
    for (Int_t iTrack = 0; iTrack < trackTree->GetEntriesFast(); iTrack++) {
      trackTree->GetEntry(iTrack);
      if (TMath::Abs(trk->GetPt()) > 3.0) {
        AliInfo(Form("Found track in sector %2i, stack %i with pt = %3.1f, triggered", 
                     trk->GetSector(), trk->GetStack(), trk->GetPt()));
        triggered_highpt = kTRUE;
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
        triggered_jet = kTRUE;
    }
  }
  else {
    AliWarning("GTU Branch not found");
  }

  if (triggered_highpt) { 
    AliInfo("Fired high-pt trigger");
    SetInput("TRD_HIGHPT_L1");
  }

  if (triggered_jet) {
    AliInfo("Fired jet trigger");
    SetInput("TRD_JET_L1");
  }

  // cleaning up
  delete gtusim;
}
