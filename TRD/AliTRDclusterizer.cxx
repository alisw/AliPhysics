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

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliTRD.h"
#include "AliTRDclusterizer.h"
#include "AliTRDcluster.h"
#include "AliTRDrecPoint.h"
#include "AliTRDgeometry.h"
#include "AliTRDparameter.h"

ClassImp(AliTRDclusterizer)

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer():TNamed()
{
  //
  // AliTRDclusterizer default constructor
  //

  fClusterTree = NULL;
  fTRD         = 0;
  fEvent       = 0;
  fVerbose     = 0;
  fPar         = 0;

}

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer(const Text_t* name, const Text_t* title)
                  :TNamed(name, title)
{
  //
  // AliTRDclusterizer default constructor
  //

  fClusterTree = NULL;
  fEvent       = 0;
  fVerbose     = 0;
  fPar         = 0;

}

//_____________________________________________________________________________
AliTRDclusterizer::AliTRDclusterizer(const AliTRDclusterizer &c):TNamed(c)
{
  //
  // AliTRDclusterizer copy constructor
  //

  ((AliTRDclusterizer &) c).Copy(*this);

}

//_____________________________________________________________________________
AliTRDclusterizer::~AliTRDclusterizer()
{
  //
  // AliTRDclusterizer destructor
  //

}

//_____________________________________________________________________________
AliTRDclusterizer &AliTRDclusterizer::operator=(const AliTRDclusterizer &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDclusterizer &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDclusterizer::Copy(TObject &c)
{
  //
  // Copy function
  //

  ((AliTRDclusterizer &) c).fClusterTree = NULL;
  ((AliTRDclusterizer &) c).fEvent       = 0;  
  ((AliTRDclusterizer &) c).fVerbose     = fVerbose;  
  ((AliTRDclusterizer &) c).fPar         = 0;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::Open(const Char_t *name, Int_t nEvent)
{
  //
  // Opens the AliROOT file. Output and input are in the same file
  //
  fRunLoader = AliRunLoader::Open(name);
  if (!fRunLoader)
   {
     Error("Open","Can not open session for file %s.",name);
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

  // Connect the AliRoot file containing Geometry, Kine, and Hits
  fRunLoader->LoadgAlice();
  gAlice = fRunLoader->GetAliRun();

  if (!(gAlice)) {
    fRunLoader->LoadgAlice();
    gAlice = fRunLoader->GetAliRun();
      if (!(gAlice)) {
        printf("AliTRDclusterizer::OpenInput -- ");
        printf("Could not find AliRun object.\n");
        return kFALSE;
      }
  }

  fEvent = nEvent;

  // Import the Trees for the event nEvent in the file
  fRunLoader->GetEvent(fEvent);
  
  // Get the TRD object
  fTRD = (AliTRD*) gAlice->GetDetector("TRD"); 
  if (!fTRD) {
    printf("AliTRDclusterizer::OpenInput -- ");
    printf("No TRD detector object found\n");
    return kFALSE;
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::WriteClusters(Int_t det)
{
  //
  // Fills TRDcluster branch in the tree with the clusters 
  // found in detector = det. For det=-1 writes the tree. 
  //

  if ((det < -1) || (det >= AliTRDgeometry::Ndet())) {
    printf("AliTRDclusterizer::WriteClusters -- ");
    printf("Unexpected detector index %d.\n",det);
    return kFALSE;
  }
 

  TBranch *branch = fClusterTree->GetBranch("TRDcluster");
  if (!branch) {
    TObjArray *ioArray = 0;
    branch = fClusterTree->Branch("TRDcluster","TObjArray",&ioArray,32000,0);
  }

  if ((det >= 0) && (det < AliTRDgeometry::Ndet())) {

    Int_t nRecPoints = fTRD->RecPoints()->GetEntriesFast();
    TObjArray *detRecPoints = new TObjArray(400);

    for (Int_t i = 0; i < nRecPoints; i++) {
      AliTRDcluster *c = (AliTRDcluster *) fTRD->RecPoints()->UncheckedAt(i);
      if (det == c->GetDetector()) {
        detRecPoints->AddLast(c);
      }
      else {
        printf("AliTRDclusterizer::WriteClusters --");
        printf("Attempt to write a cluster with unexpected detector index\n");
      }
    }

    branch->SetAddress(&detRecPoints);
    fClusterTree->Fill();

    delete detRecPoints;

    return kTRUE;

  }

  if (det == -1) {

    printf("AliTRDclusterizer::WriteClusters -- ");
    printf("Writing the cluster tree %-18s for event %d.\n"
	  ,fClusterTree->GetName(),fEvent);
    /*
    fClusterTree->Write();
    AliTRDgeometry *geo = fTRD->GetGeometry();
    geo->SetName("TRDgeometry");
    geo->Write();
    fPar->Write();
    */
    AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
    loader->WriteRecPoints("OVERWRITE");
  
    return kTRUE;  

  }
  /*
  AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
  loader->WriteDigits("OVERWRITE");
  */
  printf("AliTRDclusterizer::WriteClusters -- ");
  printf("Unexpected detector index %d.\n",det);
 
  return kFALSE;  
  
}



