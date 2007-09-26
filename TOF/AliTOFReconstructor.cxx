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

#include "AliTOFClusterFinder.h"
#include "AliTOFGeometry.h"
#include "AliTOFGeometryV5.h"
#include "AliTOFcalib.h"
#include "AliTOFtrackerMI.h"
#include "AliTOFtracker.h"
#include "AliTOFReconstructor.h"

class TTree;

class AliESDEvent;

extern TDirectory *gDirectory;
extern TFile *gFile;

ClassImp(AliTOFReconstructor)

 //____________________________________________________________________
AliTOFReconstructor::AliTOFReconstructor() 
  : AliReconstructor(),
    fTOFGeometry(0),
    fTOFcalib(0)
{
//
// ctor
//
  //Retrieving the TOF calibration info  
  fTOFGeometry = new AliTOFGeometryV5();
  fTOFcalib    = new AliTOFcalib(fTOFGeometry);
  if(!fTOFcalib->ReadParFromCDB("TOF/Calib",-1)) {AliFatal("Exiting, no CDB object found!!!");exit(0);}  
}

//------------------------------------------------------------------------
AliTOFReconstructor::AliTOFReconstructor(const AliTOFReconstructor &source)
  : AliReconstructor(),
    fTOFGeometry(0),
    fTOFcalib(0)
{
//
// copy ctor
//
  this->fTOFGeometry=source.fTOFGeometry;
  this->fTOFcalib=source.fTOFcalib;
}

//------------------------------------------------------------------------
AliTOFReconstructor & AliTOFReconstructor::operator=(const AliTOFReconstructor &source)
{
//
// assignment op.
//
  this->fTOFGeometry=source.fTOFGeometry;
  this->fTOFcalib=source.fTOFcalib;
  return *this;
}
//_____________________________________________________________________________
AliTOFReconstructor::~AliTOFReconstructor() 
{
//
// dtor
//
  delete fTOFGeometry;
  delete fTOFcalib;
}

//_____________________________________________________________________________
void AliTOFReconstructor::Reconstruct(AliRawReader *rawReader,
                                      TTree *clustersTree) const
{
// reconstruct clusters from Raw Data

  AliTOFClusterFinder tofClus(fTOFcalib);
  tofClus.Digits2RecPoints(rawReader, clustersTree);

}

//_____________________________________________________________________________
void AliTOFReconstructor::Reconstruct(TTree *digitsTree,
                                      TTree *clustersTree) const
{
// reconstruct clusters from Raw Data

  AliDebug(2,Form("Global Event loop mode: Creating Recpoints from Digits Tree")); 
  AliTOFClusterFinder tofClus(fTOFcalib);
  tofClus.Digits2RecPoints(digitsTree, clustersTree);

}
//_____________________________________________________________________________
  void AliTOFReconstructor::ConvertDigits(AliRawReader* reader, TTree* digitsTree) const
{
// reconstruct clusters from digits

  AliDebug(2,Form("Global Event loop mode: Converting Raw Data to a Digits Tree")); 
  AliTOFClusterFinder tofClus(fTOFcalib);
  tofClus.Raw2Digits(reader, digitsTree);

}

//_____________________________________________________________________________
AliTracker* AliTOFReconstructor::CreateTracker() const
{
// create a TOF tracker

  TString selectedTracker = GetOption();
  // use MI tracker if selected
  if (selectedTracker.Contains("MI")) return new AliTOFtrackerMI();
  return new AliTOFtracker();
}

//_____________________________________________________________________________
AliTOFGeometry* AliTOFReconstructor::GetTOFGeometry() const
{
  // get the TOF geometry
  return fTOFGeometry;
}
