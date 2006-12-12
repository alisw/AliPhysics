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
// class for TOF reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"

#include "AliLog.h"
#include "AliRawReader.h"
#include "AliRunLoader.h"

#include "AliTOFClusterFinder.h"
#include "AliTOFGeometry.h"
#include "AliTOFGeometryV5.h"
#include "AliTOFtrackerMI.h"
#include "AliTOFtracker.h"
#include "AliTOFReconstructor.h"

class TTree;

class AliESD;

extern TDirectory *gDirectory;
extern TFile *gFile;

ClassImp(AliTOFReconstructor)

//_____________________________________________________________________________
  void AliTOFReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
// reconstruct clusters from digits

  AliTOFClusterFinder tofClus(runLoader);
  tofClus.Load();
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++)
    {
      tofClus.Digits2RecPoints(iEvent);
    }
  tofClus.UnLoad();

}

//_____________________________________________________________________________
void AliTOFReconstructor::Reconstruct(AliRunLoader* runLoader,
                                      AliRawReader *rawReader) const
{
// reconstruct clusters from Raw Data

  AliTOFClusterFinder tofClus(runLoader);
  tofClus.LoadClusters();
  Int_t iEvent = 0;
  while (rawReader->NextEvent()) {
    tofClus.Digits2RecPoints(iEvent,rawReader);
    //tofClus.Raw2Digits(iEvent,rawReader); // temporary solution
    iEvent++;
  }
  tofClus.UnLoadClusters();

}

//_____________________________________________________________________________
void AliTOFReconstructor::Reconstruct(AliRawReader *rawReader,
                                      TTree *clustersTree) const
{
// reconstruct clusters from Raw Data

  AliTOFClusterFinder tofClus;
  tofClus.Digits2RecPoints(rawReader, clustersTree);

}

//_____________________________________________________________________________
AliTracker* AliTOFReconstructor::CreateTracker(AliRunLoader* runLoader) const
{
// create a TOF tracker

//  AliTOFGeometry* geom = GetTOFGeometry(runLoader);
  AliTOFGeometry* geom = new AliTOFGeometryV5();
  if (!geom) return NULL;
  //  Double_t parPID[] = {130., 5.};
  Double_t parPID[] = {80., 5.};
  TString selectedTracker = GetOption();
  // use MI tracker if selected
  if (selectedTracker.Contains("MI")) return new AliTOFtrackerMI(geom,parPID);

  return new AliTOFtracker(geom, parPID);
}

//_____________________________________________________________________________
void AliTOFReconstructor::FillESD(AliRunLoader* /*runLoader*/, 
				  AliESD* /*esd*/) const
{
// nothing to be done

}

//_____________________________________________________________________________
AliTOFGeometry* AliTOFReconstructor::GetTOFGeometry(AliRunLoader* runLoader) const
{
// get the TOF parameters

  AliTOFGeometry *tofGeom;

  runLoader->CdGAFile();
  TDirectory *savedir=gDirectory; 
  TFile *in=(TFile*)gFile;  
  if (!in->IsOpen()) {
    AliWarning("Geometry file is not open default  TOF geometry will be used");
    tofGeom = new AliTOFGeometryV5();
  }
  else {
    in->cd();  
    tofGeom = (AliTOFGeometry*) in->Get("TOFgeometry");
  }

  savedir->cd();  

  if (!tofGeom) {
    AliError("no TOF geometry available");
    return NULL;
  }
  return tofGeom;
}
