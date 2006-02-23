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


#include "AliTRDReconstructor.h"
#include "AliRunLoader.h"
#include "AliTRDparameter.h"
#include "AliTRDclusterizerV1.h"
#include "AliTRDtracker.h"
#include "AliTRDpidESD.h"
#include <TFile.h>
#include "AliRawReader.h"
#include "AliLog.h"

ClassImp(AliTRDReconstructor)


//_____________________________________________________________________________
void AliTRDReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
// reconstruct clusters

  AliLoader *loader=runLoader->GetLoader("TRDLoader");
  loader->LoadRecPoints("recreate");

  AliTRDclusterizerV1 clusterer("clusterer", "TRD clusterizer");
  runLoader->CdGAFile();
  AliTRDparameter* trdParam = GetTRDparameter(runLoader); 
  if (!trdParam) {
    Error("Reconstruct", "no TRD parameters found");
    return;
  }
  trdParam->ReInit();
  clusterer.SetParameter(trdParam);
  Int_t nEvents = runLoader->GetNumberOfEvents();

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    clusterer.Open(runLoader->GetFileName(), iEvent);
    clusterer.ReadDigits();
    clusterer.MakeClusters();
    clusterer.WriteClusters(-1);
  }

  loader->UnloadRecPoints();
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
  AliTRDparameter* trdParam = GetTRDparameter(runLoader); 
  if (!trdParam) {
    Error("Reconstruct", "no TRD parameters found");
    return;
  }
  trdParam->ReInit();
  clusterer.SetParameter(trdParam);
  Int_t nEvents = runLoader->GetNumberOfEvents();

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    if (!rawReader->NextEvent()) break;
    clusterer.Open(runLoader->GetFileName(), iEvent);
    clusterer.ReadDigits(rawReader);
    clusterer.MakeClusters();
    clusterer.WriteClusters(-1);
  }

  loader->UnloadRecPoints();
}

//_____________________________________________________________________________
AliTracker* AliTRDReconstructor::CreateTracker(AliRunLoader* runLoader) const
{
// create a TRD tracker

  runLoader->CdGAFile();
  return new AliTRDtracker(gFile);
}

//_____________________________________________________________________________
void AliTRDReconstructor::FillESD(AliRunLoader* /*runLoader*/, 
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
}


//_____________________________________________________________________________
AliTRDparameter* AliTRDReconstructor::GetTRDparameter(AliRunLoader* runLoader) const
{
// get the TRD parameters

  runLoader->CdGAFile();
  AliTRDparameter* trdParam = (AliTRDparameter*) gFile->Get("TRDparameter"); 
  if (!trdParam) {
    Error("GetTRDparameter", "no TRD parameters available");
    return NULL;
  }
  return trdParam;
}


