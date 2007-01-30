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
//  TRD cluster finder base class                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TObjArray.h>

#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliLog.h"

#include "AliTRDclusterizer.h"
#include "AliTRDcluster.h"
#include "AliTRDrecPoint.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"

ClassImp(AliTRDclusterizer)

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer()
  :TNamed()
  ,fRunLoader(NULL)
  ,fClusterTree(NULL)
  ,fRecPoints(NULL)
{
  //
  // AliTRDclusterizer default constructor
  //

}

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer(const Text_t* name, const Text_t* title)
  :TNamed(name,title)
  ,fRunLoader(NULL)
  ,fClusterTree(NULL)
  ,fRecPoints(NULL)
{
  //
  // AliTRDclusterizer constructor
  //

}

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer(const AliTRDclusterizer &c)
  :TNamed(c)
  ,fRunLoader(NULL)
  ,fClusterTree(NULL)
  ,fRecPoints(NULL)
{
  //
  // AliTRDclusterizer copy constructor
  //

}

//_____________________________________________________________________________
AliTRDclusterizer::~AliTRDclusterizer()
{
  //
  // AliTRDclusterizer destructor
  //

  if (fRecPoints) {
    fRecPoints->Delete();
    delete fRecPoints;
  }

}

//_____________________________________________________________________________
AliTRDclusterizer &AliTRDclusterizer::operator=(const AliTRDclusterizer &c)
{
  //
  // Assignment operator
  //

  if (this != &c) {
    ((AliTRDclusterizer &) c).Copy(*this);
  }
  return *this;

}

//_____________________________________________________________________________
void AliTRDclusterizer::Copy(TObject &c) const
{
  //
  // Copy function
  //

  ((AliTRDclusterizer &) c).fClusterTree = NULL;
  ((AliTRDclusterizer &) c).fRecPoints   = NULL;  

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::Open(const Char_t *name, Int_t nEvent)
{
  //
  // Opens the AliROOT file. Output and input are in the same file
  //

  TString evfoldname = AliConfig::GetDefaultEventFolderName();
  fRunLoader         = AliRunLoader::GetRunLoader(evfoldname);

  if (!fRunLoader) {
    fRunLoader = AliRunLoader::Open(name);
  }

  if (!fRunLoader) {
    AliError(Form("Can not open session for file %s.",name));
    return kFALSE;
  }

  OpenInput(nEvent);
  OpenOutput();

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::OpenOutput()
{
  //
  // Open the output file
  //

  TObjArray *ioArray = 0;

  AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
  loader->MakeTree("R");

  fClusterTree = loader->TreeR();
  fClusterTree->Branch("TRDcluster","TObjArray",&ioArray,32000,0);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::OpenInput(Int_t nEvent)
{
  //
  // Opens a ROOT-file with TRD-hits and reads in the digits-tree
  //

  // Import the Trees for the event nEvent in the file
  fRunLoader->GetEvent(nEvent);
  
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::WriteClusters(Int_t det)
{
  //
  // Fills TRDcluster branch in the tree with the clusters 
  // found in detector = det. For det=-1 writes the tree. 
  //

  if ((det <                      -1) || 
      (det >= AliTRDgeometry::Ndet())) {
    AliError(Form("Unexpected detector index %d.\n",det));
    return kFALSE;
  }
 
  TBranch *branch = fClusterTree->GetBranch("TRDcluster");
  if (!branch) {
    TObjArray *ioArray = 0;
    branch = fClusterTree->Branch("TRDcluster","TObjArray",&ioArray,32000,0);
  }

  if ((det >=                      0) && 
      (det <  AliTRDgeometry::Ndet())) {

    Int_t nRecPoints = RecPoints()->GetEntriesFast();
    TObjArray *detRecPoints = new TObjArray(400);

    for (Int_t i = 0; i < nRecPoints; i++) {
      AliTRDcluster *c = (AliTRDcluster *) RecPoints()->UncheckedAt(i);
      if (det == c->GetDetector()) {
        detRecPoints->AddLast(c);
      }
      else {
        AliError("Attempt to write a cluster with unexpected detector index\n");
      }
    }

    branch->SetAddress(&detRecPoints);
    fClusterTree->Fill();

    delete detRecPoints;

    return kTRUE;

  }

  if (det == -1) {

    AliInfo(Form("Writing the cluster tree %s for event %d."
	        ,fClusterTree->GetName(),fRunLoader->GetEventNumber()));

    if (fRecPoints) {

      branch->SetAddress(&fRecPoints);

      AliLoader *loader = fRunLoader->GetLoader("TRDLoader");
      loader->WriteRecPoints("OVERWRITE");
  
    }
    else {

      AliError("Cluster tree does not exist. Cannot write clusters.\n");
      return kFALSE;

    }

    return kTRUE;  

  }

  AliError(Form("Unexpected detector index %d.\n",det));
 
  return kFALSE;  
  
}


//_____________________________________________________________________________
AliTRDcluster* AliTRDclusterizer::AddCluster(Double_t *pos, Int_t timebin
                                           , Int_t det, Double_t amp
				           , Int_t *tracks, Double_t *sig
                                           , Int_t iType, Float_t center)
{
  //
  // Add a cluster for the TRD
  //

  AliTRDcluster *c = new AliTRDcluster();

  c->SetDetector(det);
  c->SetQ(amp);
  c->SetX(pos[2]);
  c->SetY(pos[0]);
  c->SetZ(pos[1]);
  c->SetSigmaY2(sig[0]);   
  c->SetSigmaZ2(sig[1]);
  c->SetLocalTimeBin(timebin);
  c->SetCenter(center);

  if (tracks) {
    c->AddTrackIndex(tracks);
  }

  switch (iType) {
  case 0:
    c->Set2pad();
    break;
  case 1:
    c->Set3pad();
    break;
  case 2:
    c->Set4pad();
    break;
  case 3:
    c->Set5pad();
    break;
  case 4:
    c->SetLarge();
    break;
  };

  RecPoints()->Add(c);
  return c;

}

//_____________________________________________________________________________
Double_t AliTRDclusterizer::CalcXposFromTimebin(Float_t timebin, Int_t idet
                                              , Int_t col, Int_t row)
{
  //
  // Calculates the local x position in the detector from the timebin, 
  // depends on the drift velocity and t0
  //
  
  AliTRDcalibDB     *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("Cannot find calibration object");
    return -1;
  }
  AliTRDCommonParam *parCom      = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliError("Could not get common parameters\n");
    return kFALSE;
  }

  Float_t vdrift            = calibration->GetVdrift(idet,col,row);  
  Float_t t0                = calibration->GetT0(idet,col,row);
  Float_t samplingFrequency = parCom->GetSamplingFrequency();

  timebin -= t0;

  return timebin / samplingFrequency * vdrift;

}

//_____________________________________________________________________________
void AliTRDclusterizer::ResetRecPoints() 
{
  //
  // Resets the list of rec points
  //

  if (fRecPoints) {
    fRecPoints->Delete();
  }

}

//_____________________________________________________________________________
TObjArray* AliTRDclusterizer::RecPoints() 
{
  //
  // Returns the list of rec points
  //

  if (!fRecPoints) {
    fRecPoints = new TObjArray(400);
  }
 
  return fRecPoints;

}
