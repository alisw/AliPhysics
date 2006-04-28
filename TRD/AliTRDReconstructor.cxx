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

#include "AliTRDReconstructor.h"
#include "AliRunLoader.h"
#include "AliTRDclusterizerV1.h"
#include "AliTRDtracker.h"
#include "AliTRDpidESD.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliTRDtrigger.h"
#include "AliTRDtrigParam.h"
#include "AliTRDgtuTrack.h"
#include "AliRun.h"
#include "AliESDTrdTrack.h"
#include "AliESD.h"

ClassImp(AliTRDReconstructor)

Bool_t  AliTRDReconstructor::fgkSeedingOn = kFALSE;
Int_t   AliTRDReconstructor::fgStreamLevel     = 0;      // stream (debug) level

//_____________________________________________________________________________
void AliTRDReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
// reconstruct clusters

  AliLoader *loader=runLoader->GetLoader("TRDLoader");
  loader->LoadRecPoints("recreate");

  AliTRDclusterizerV1 clusterer("clusterer", "TRD clusterizer");
  runLoader->CdGAFile();
  Int_t nEvents = runLoader->GetNumberOfEvents();

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    clusterer.Open(runLoader->GetFileName(), iEvent);
    clusterer.ReadDigits();
    clusterer.MakeClusters();
    clusterer.WriteClusters(-1);
  }

  loader->UnloadRecPoints();

  // Trigger (tracklets, LTU)

  loader->LoadTracks("UPDATE");
  if (loader->TreeT()) {
    Info("Reconstruct","Tracklets already exist");
    return;
  }
  Info("Reconstruct","Trigger tracklets will be produced");

  AliTRDtrigger trdTrigger("Trigger","Trigger class"); 

  AliTRDtrigParam *trigp = new AliTRDtrigParam("TRDtrigParam","TRD Trigger parameters");

  if (runLoader->GetAliRun() == 0x0) runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  Double_t x[3] = { 0.0, 0.0, 0.0 };
  Double_t b[3];
  gAlice->Field(x,b);  // b[] is in kilo Gauss
  Float_t field = b[2] * 0.1; // Tesla
  Info("Reconstruct","Trigger set for magnetic field = %f Tesla \n",field);

  trigp->SetField(field);
  trigp->Init();
  trdTrigger.SetParameter(trigp);

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    trdTrigger.Open(runLoader->GetFileName(), iEvent);
    trdTrigger.ReadDigits();
    trdTrigger.MakeTracklets();
    trdTrigger.WriteTracklets(-1);
  }

  loader->UnloadTracks();

}

//_____________________________________________________________________________
void AliTRDReconstructor::Reconstruct(AliRunLoader* runLoader,
                                      AliRawReader* rawReader) const
{
// reconstruct clusters

  AliInfo("Reconstruct TRD clusters from RAW data");

  AliLoader *loader=runLoader->GetLoader("TRDLoader");
  loader->LoadRecPoints("recreate");

  AliTRDclusterizerV1 clusterer("clusterer", "TRD clusterizer");
  runLoader->CdGAFile();
  Int_t nEvents = runLoader->GetNumberOfEvents();

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    if (!rawReader->NextEvent()) break;
    clusterer.Open(runLoader->GetFileName(), iEvent);
    clusterer.ReadDigits(rawReader);
    clusterer.MakeClusters();
    clusterer.WriteClusters(-1);
  }

  loader->UnloadRecPoints();

  // Trigger (tracklets, LTU)

  loader->LoadTracks();
  if (loader->TreeT()) {
    Info("Reconstruct","Tracklets already exist");
    return;
  }
  Info("Reconstruct","Trigger tracklets will be produced");

  AliTRDtrigger trdTrigger("Trigger","Trigger class"); 

  AliTRDtrigParam *trigp = new AliTRDtrigParam("TRDtrigParam","TRD Trigger parameters");

  if (runLoader->GetAliRun() == 0x0) runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  Double_t x[3] = { 0.0, 0.0, 0.0 };
  Double_t b[3];
  gAlice->Field(x,b);  // b[] is in kilo Gauss
  Float_t field = b[2] * 0.1; // Tesla
  Info("Reconstruct","Trigger set for magnetic field = %f Tesla \n",field);

  trigp->SetField(field);
  trigp->Init();
  trdTrigger.SetParameter(trigp);

  rawReader->RewindEvents();

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    if (!rawReader->NextEvent()) break;
    trdTrigger.Open(runLoader->GetFileName(), iEvent);
    trdTrigger.ReadDigits(rawReader);
    trdTrigger.MakeTracklets();
    trdTrigger.WriteTracklets(-1);
  }

  loader->UnloadTracks();

}

//_____________________________________________________________________________
AliTracker* AliTRDReconstructor::CreateTracker(AliRunLoader* runLoader) const
{
// create a TRD tracker

  runLoader->CdGAFile();
  return new AliTRDtracker(gFile);
}

//_____________________________________________________________________________
void AliTRDReconstructor::FillESD(AliRunLoader* runLoader, 
				  AliESD* esd) const
{
// make PID

  Double_t parTRD[] = {
    280., // Min. Ionizing Particle signal.  Check it !!!
    0.23, // relative resolution             Check it !!!
    10.   // PID range (in sigmas)
  };
  AliTRDpidESD trdPID(parTRD);
  trdPID.MakePID(esd);

  // Trigger (tracks, GTU)

  AliTRDtrigger trdTrigger("Trigger","Trigger class"); 

  AliTRDtrigParam *trigp = new AliTRDtrigParam("TRDtrigParam","TRD Trigger parameters");

  if (runLoader->GetAliRun() == 0x0) runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  Double_t x[3] = { 0.0, 0.0, 0.0 };
  Double_t b[3];
  gAlice->Field(x,b);  // b[] is in kilo Gauss
  Float_t field = b[2] * 0.1; // Tesla
  Info("FillESD","Trigger set for magnetic field = %f Tesla \n",field);

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



