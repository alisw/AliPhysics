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

Revision 1.1.4.5  2000/10/15 23:40:01  cblume
Remove AliTRDconst

Revision 1.1.4.4  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.4.3  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.1.4.2  2000/09/22 14:49:49  cblume
Adapted to tracking code

Revision 1.5  2000/10/02 21:28:19  fca
Removal of useless dependecies via forward declarations

Revision 1.4  2000/06/09 11:10:07  cblume
Compiler warnings and coding conventions, next round

Revision 1.3  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.1  2000/05/08 15:08:03  cblume
Remove the class AliTRDcluster

Revision 1.4  2000/06/09 11:10:07  cblume
Compiler warnings and coding conventions, next round

Revision 1.3  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.1  2000/05/08 15:08:03  cblume
Remove the class AliTRDcluster

Revision 1.1  2000/02/28 18:57:58  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD cluster finder base class                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>

#include "AliRun.h"
#include "AliTRD.h"
#include "AliTRDclusterizer.h"
#include "AliTRDrecPoint.h"
#include "AliTRDgeometry.h"

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
AliTRDclusterizer::AliTRDclusterizer(const AliTRDclusterizer &c)
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

  if (fInputFile) {
    fInputFile->Close();
    delete fInputFile;
  }

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

  ((AliTRDclusterizer &) c).fInputFile = NULL;
  ((AliTRDclusterizer &) c).fEvent     = 0;  

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
Bool_t AliTRDclusterizer::WriteClusters(Int_t det)
{
  //
  // Fills TRDrecPoints branch in TRDrecPoints## tree with rec. points 
  // found in detector = det. For det=-1 writes the tree. 
  // For det=-2 recreates the tree.

  Char_t treeName[14];
  sprintf(treeName,"TRDrecPoints%d", fEvent);

  if (det == -2) {
    fInputFile->Delete(treeName);
    TTree *tree = new TTree(treeName,"Tree with TRD rec. points");
    tree->Write();
    return kTRUE;
  }

  TTree *tree=(TTree*)fInputFile->Get(treeName);
  TBranch *branch=tree->GetBranch("TRDrecPoints");

  if(!branch) {
    TObjArray *ioArray = 0;
    branch = tree->Branch("TRDrecPoints","TObjArray",&ioArray,32000,0);
  }

  if ((det >= 0) && (det < AliTRDgeometry::Ndet())) {

    AliTRD *TRD = (AliTRD*) gAlice->GetDetector("TRD");
    Int_t nRecPoints = TRD->RecPoints()->GetEntriesFast();
    TObjArray *fDetRecPoints = new TObjArray(400);

    for (Int_t i=0; i<nRecPoints; i++) {
      AliTRDrecPoint *p=(AliTRDrecPoint*)TRD->RecPoints()->UncheckedAt(i);
      if(det == p->GetDetector()) fDetRecPoints->AddLast(p);
      else printf("attempt to write a RecPoint with unexpected detector index");
    }

    branch->SetAddress(&fDetRecPoints);
    tree->Fill();
    return kTRUE;
  }

  if (det == -1) {

    printf("\rAliTRDclusterizer::WriteClusters -- ");
    printf("Writing the cluster tree %-18s for event %d.\n"
	   ,tree->GetName(),fEvent);

    tree->Write();     
    return kTRUE;  
  }
  
  printf("\rAliTRDclusterizer::WriteClusters -- ");
  printf("Unexpected detector index %d.\n", det); 
  return kFALSE;  

}



