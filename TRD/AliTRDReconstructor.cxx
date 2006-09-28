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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TRD reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TFile.h>

#include "AliRunLoader.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliESDTrdTrack.h"
#include "AliESD.h"

#include "AliTRDReconstructor.h"
#include "AliTRDclusterizerV1.h"
#include "AliTRDtracker.h"
#include "AliTRDpidESD.h"
#include "AliTRDtrigger.h"
#include "AliTRDtrigParam.h"
#include "AliTRDgtuTrack.h"

ClassImp(AliTRDReconstructor)

Bool_t AliTRDReconstructor::fgkSeedingOn  = kFALSE;
Int_t  AliTRDReconstructor::fgStreamLevel = 0;      // Stream (debug) level

//_____________________________________________________________________________
void AliTRDReconstructor::Reconstruct(AliRunLoader *runLoader) const
{
  //
  // Reconstruct clusters
  //

  AliLoader *loader = runLoader->GetLoader("TRDLoader");
  loader->LoadRecPoints("recreate");

  runLoader->CdGAFile();
  Int_t nEvents = runLoader->GetNumberOfEvents();

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    AliTRDclusterizerV1 clusterer("clusterer","TRD clusterizer");
    clusterer.Open(runLoader->GetFileName(),iEvent);
    clusterer.ReadDigits();
    clusterer.MakeClusters();
    clusterer.WriteClusters(-1);
  }

  loader->UnloadRecPoints();

  //
  // Trigger (tracklets, LTU)
  //
  loader->LoadTracks("RECREATE");
  AliInfo("Trigger tracklets will be produced");

  AliTRDtrigger trdTrigger("Trigger","Trigger class"); 

  AliTRDtrigParam *trigp = new AliTRDtrigParam("TRDtrigParam"
                                              ,"TRD Trigger parameters");

  Float_t field = AliTracker::GetBz() * 0.1; // Tesla
  AliInfo(Form("Trigger set for magnetic field = %f Tesla \n",field));
  trigp->SetField(field);
  trigp->Init();
  trdTrigger.SetParameter(trigp);

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    trdTrigger.Open(runLoader->GetFileName(),iEvent);
    trdTrigger.ReadDigits();
    trdTrigger.MakeTracklets();
    trdTrigger.WriteTracklets(-1);
  }

  loader->UnloadTracks();

}

//_____________________________________________________________________________
void AliTRDReconstructor::Reconstruct(AliRunLoader *runLoader
                                    , AliRawReader *rawReader) const
{
  //
  // Reconstruct clusters
  //

  AliInfo("Reconstruct TRD clusters from RAW data");

  AliLoader *loader = runLoader->GetLoader("TRDLoader");
  loader->LoadRecPoints("recreate");

  runLoader->CdGAFile();
  Int_t nEvents = runLoader->GetNumberOfEvents();

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    if (!rawReader->NextEvent()) break;
    AliTRDclusterizerV1 clusterer("clusterer","TRD clusterizer");
    clusterer.Open(runLoader->GetFileName(),iEvent);
    clusterer.ReadDigits(rawReader);
    clusterer.MakeClusters();
    clusterer.WriteClusters(-1);
  }

  loader->UnloadRecPoints();

  //
  // Trigger (tracklets, LTU)
  //
  loader->LoadTracks();
  if (loader->TreeT()) {
    AliError("Tracklets already exist");
    return;
  }
  AliInfo("Trigger tracklets will be produced");

  AliTRDtrigger trdTrigger("Trigger","Trigger class"); 

  AliTRDtrigParam *trigp = new AliTRDtrigParam("TRDtrigParam"
                                              ,"TRD Trigger parameters");

  Float_t field = AliTracker::GetBz() * 0.1; // Tesla
  AliInfo(Form("Trigger set for magnetic field = %f Tesla \n",field));
  trigp->SetField(field);
  trigp->Init();
  trdTrigger.SetParameter(trigp);

  rawReader->RewindEvents();

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    if (!rawReader->NextEvent()) break;
    trdTrigger.Open(runLoader->GetFileName(),iEvent);
    trdTrigger.ReadDigits(rawReader);
    trdTrigger.MakeTracklets();
    trdTrigger.WriteTracklets(-1);
  }

  loader->UnloadTracks();

}

//_____________________________________________________________________________
AliTracker *AliTRDReconstructor::CreateTracker(AliRunLoader *runLoader) const
{
  //
  // Create a TRD tracker
  //

  runLoader->CdGAFile();

  return new AliTRDtracker(gFile);

}

//_____________________________________________________________________________
void AliTRDReconstructor::FillESD(AliRunLoader *runLoader
				, AliESD *esd) const
{
  //
  // Make PID
  //

  AliTRDpidESD trdPID;
  trdPID.MakePID(esd);

  //
  // Trigger (tracks, GTU)
  //
  AliTRDtrigger trdTrigger("Trigger","Trigger class"); 

  AliTRDtrigParam *trigp = new AliTRDtrigParam("TRDtrigParam"
                                              ,"TRD Trigger parameters");

  Float_t field = AliTracker::GetBz() * 0.1; // Tesla
  AliInfo(Form("Trigger set for magnetic field = %f Tesla \n",field));
  trigp->SetField(field);
  trigp->Init();

  trdTrigger.SetParameter(trigp);
  trdTrigger.SetRunLoader(runLoader);
  trdTrigger.Init();

  Int_t iEvent = runLoader->GetEventNumber(); 
  runLoader->GetEvent(iEvent);
  trdTrigger.ReadTracklets(runLoader);

  AliESDTrdTrack *TrdTrack = new AliESDTrdTrack();
  AliTRDgtuTrack *GtuTrack;

  Int_t nTracks = trdTrigger.GetNumberOfTracks();
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

    GtuTrack = trdTrigger.GetTrack(iTrack);

    TrdTrack->SetYproj(GtuTrack->GetYproj());
    TrdTrack->SetZproj(GtuTrack->GetZproj());
    TrdTrack->SetSlope(GtuTrack->GetSlope());
    TrdTrack->SetDetector(GtuTrack->GetDetector());
    TrdTrack->SetTracklets(GtuTrack->GetTracklets());
    TrdTrack->SetPlanes(GtuTrack->GetPlanes());
    TrdTrack->SetClusters(GtuTrack->GetClusters());
    TrdTrack->SetPt(GtuTrack->GetPt());
    TrdTrack->SetPhi(GtuTrack->GetPhi());
    TrdTrack->SetEta(GtuTrack->GetEta());
    TrdTrack->SetLabel(GtuTrack->GetLabel());
    TrdTrack->SetPID(GtuTrack->GetPID());
    TrdTrack->SetIsElectron(GtuTrack->IsElectron());

    esd->AddTrdTrack(TrdTrack);

  }

}
