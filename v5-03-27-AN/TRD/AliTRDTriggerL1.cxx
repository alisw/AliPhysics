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
#include <TTree.h>

#include "AliLog.h"
#include "AliTriggerInput.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliTRDTriggerL1.h"
#include "AliTRDgtuSim.h"
#include "AliTRDtrackGTU.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCalDCSGTU.h"

AliTRDTriggerL1::AliTRDTriggerL1() :
  AliTriggerDetector(),
  fPtThresholdA(3.),
  fPtThresholdB(5.),
  fPidThresholdA(0),
  fPidThresholdB(185),
  fNoThreshold(1),
  fNoThresholdA(3),
  fNoThresholdB(200),
  fNoThresholdJetA(3),
  fNoThresholdJetB(200),
  fNoThresholdElA(1),
  fNoThresholdElB(1)
{
  // ctor

  SetName("TRD");
}

AliTRDTriggerL1::~AliTRDTriggerL1()
{
  // dtor
}

void AliTRDTriggerL1::CreateInputs()
{
  // create the trigger inputs for TRD

  if (fInputs.GetEntriesFast() > 0)
    return;

  fInputs.AddLast(new AliTriggerInput("1HCO", "TRD", 1));
  fInputs.AddLast(new AliTriggerInput("1HJT", "TRD", 1));
  fInputs.AddLast(new AliTriggerInput("1HSE", "TRD", 1));
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

  TTree *trackTree = trdLoader->GetDataLoader("gtutracks")->Tree();
  if (!trackTree) {
    AliDebug(1,"Did not find track tree");
    return;
  }
  TBranch *branch = trackTree->GetBranch("TRDtrackGTU");
  AliDebug(1,Form("TRD trigger: found %lld tracks", trackTree->GetEntriesFast()));

  // trigger algorithms to come, e.g.
  Bool_t triggered1HCO    = kFALSE;
  Bool_t triggered1HJT    = kFALSE;
  Bool_t triggered1HSE    = kFALSE;

  if (branch) {
    AliTRDtrackGTU *trk = 0x0;
    branch->SetAddress(&trk);

    Int_t nTracks[90]      = { 0 }; // number of tracks
    Int_t nTracksA[90]     = { 0 }; // number of tracks above pt threshold A
    Int_t nTracksB[90]     = { 0 }; // number of tracks above pt threshold B
    Int_t nTracksElA[90]   = { 0 }; // number of tracks above pt threshold A and PID threshold A
    Int_t nTracksElB[90]   = { 0 }; // number of tracks above pt threshold B and PID threshold B

    for (Int_t iTrack = 0; iTrack < trackTree->GetEntriesFast(); iTrack++) {
      trackTree->GetEntry(iTrack);

      nTracks[5*trk->GetSector() + trk->GetStack()]++;

      if (TMath::Abs(trk->GetPt()) > fPtThresholdA) {
        nTracksA[5*trk->GetSector() + trk->GetStack()]++;
	if (trk->GetPID() > fPidThresholdA)
	  nTracksElA[5*trk->GetSector() + trk->GetStack()]++;
      }

      if (TMath::Abs(trk->GetPt()) > fPtThresholdB) {
        nTracksB[5*trk->GetSector() + trk->GetStack()]++;
	if (trk->GetPID() > fPidThresholdB)
	  nTracksElB[5*trk->GetSector() + trk->GetStack()]++;
      }
    }

    for (Int_t iStack = 0; iStack < 90; iStack++) {
      if ((nTracksA[iStack] >= fNoThresholdJetA) || (nTracksB[iStack] >= fNoThresholdJetB))
        triggered1HJT = kTRUE;

      if ((nTracksElA[iStack] >= fNoThresholdElA))
        triggered1HCO = kTRUE;

      if ((nTracksElB[iStack] >= fNoThresholdElB))
        triggered1HSE = kTRUE;
    }
  }
  else {
    AliWarning("GTU Branch not found");
  }

  if (triggered1HCO) {
    AliDebug(1, "Fired cosmic trigger");
    SetInput("1HCO");
  }

  if (triggered1HJT) {
    AliDebug(1, "Fired jet trigger");
    SetInput("1HJT");
  }

  if (triggered1HSE) {
    AliDebug(1, "Fired single electron trigger");
    SetInput("1HSE");
  }

  // cleaning up
  delete gtusim;
}
