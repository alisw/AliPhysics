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

/*
$Log$
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD cluster finder base class                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliRun.h"

#include "AliTRD.h"
#include "AliTRDclusterizer.h"

ClassImp(AliTRDclusterizer)

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer():TNamed()
{
  //
  // AliTRDclusterizer default constructor
  //

  fInputFile = NULL;
  fEvent     = 0;

}

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer(const Text_t* name, const Text_t* title)
                  :TNamed(name, title)
{
  //
  // AliTRDclusterizer default constructor
  //

  fInputFile = NULL;
  fEvent     = 0;

  Init();

}

//_____________________________________________________________________________
AliTRDclusterizer::~AliTRDclusterizer()
{

  if (fInputFile) {
    fInputFile->Close();
    delete fInputFile;
  }

}

//_____________________________________________________________________________
void AliTRDclusterizer::Init()
{
  //
  // Initializes the cluster finder
  //

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::Open(const Char_t *name, Int_t nEvent)
{
  //
  // Opens a ROOT-file with TRD-hits and reads in the digits-tree
  //

  // Connect the AliRoot file containing Geometry, Kine, and Hits
  fInputFile = (TFile*) gROOT->GetListOfFiles()->FindObject(name);
  if (!fInputFile) {
    printf("AliTRDclusterizer::Open -- ");
    printf("Open the ALIROOT-file %s.\n",name);
    fInputFile = new TFile(name,"UPDATE");
  }
  else {
    printf("AliTRDclusterizer::Open -- ");
    printf("%s is already open.\n",name);
  }

  // Get AliRun object from file or create it if not on file
  //if (!gAlice) {
    gAlice = (AliRun*) fInputFile->Get("gAlice");
    if (gAlice) {
      printf("AliTRDclusterizer::Open -- ");
      printf("AliRun object found on file.\n");
    }
    else {
      printf("AliTRDclusterizer::Open -- ");
      printf("Could not find AliRun object.\n");
      return kFALSE;
    }
  //}

  fEvent = nEvent;

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(fEvent);
  if (nparticles <= 0) {
    printf("AliTRDclusterizer::Open -- ");
    printf("No entries in the trees for event %d.\n",fEvent);
    return kFALSE;
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::WriteCluster()
{
  //
  // Writes out the TRD-cluster
  //

  // Write the new tree into the input file (use overwrite option)
  Char_t treeName[7];
  sprintf(treeName,"TreeR%d",fEvent);
  printf("AliTRDclusterizer::WriteCluster -- ");
  printf("Write the cluster tree %s for event %d.\n"
        ,treeName,fEvent);
  gAlice->TreeR()->Write(treeName,2);

  return kTRUE;

}

ClassImp(AliTRDcluster)

//_____________________________________________________________________________
AliTRDcluster::AliTRDcluster(Int_t *tracks, Int_t *cluster, Float_t energy, Float_t* position)
              :TObject()
{
  //
  // Create a TRD cluster
  //

  fDetector  = cluster[0];

  fTimeSlice = cluster[1];
  fEnergy    = energy;

  fX         = position[0];
  fY         = position[1];
  fZ         = position[2];

  fTracks[0] = tracks[0];
  fTracks[1] = tracks[1];
  fTracks[2] = tracks[2];

}
