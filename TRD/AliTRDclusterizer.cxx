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
Revision 1.17  2003/02/10 11:09:26  cblume
Further modifications in OpenOutput and WriteCluster

Revision 1.16  2003/02/07 09:45:24  cblume
Modification in OpenOutput() for many events in single file

Revision 1.15  2002/11/07 15:52:09  cblume
Update of tracking code for tilted pads

Revision 1.14  2002/10/14 14:57:43  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.11.6.2  2002/10/11 07:26:37  hristov
Updating VirtualMC to v3-09-02

Revision 1.13  2002/09/18 09:20:53  cblume
Write the parameter class into the cluster file

Revision 1.12  2002/06/12 09:54:35  cblume
Update of tracking code provided by Sergei

Revision 1.11  2001/11/27 08:50:33  hristov
BranchOld replaced by Branch

Revision 1.10  2001/11/14 10:50:45  cblume
Changes in digits IO. Add merging of summable digits

Revision 1.9  2001/10/21 18:30:02  hristov
Several pointers were set to zero in the default constructors to avoid memory management problems

Revision 1.8  2001/05/07 08:06:44  cblume
Speedup of the code. Create only AliTRDcluster

Revision 1.7  2001/03/30 14:40:14  cblume
Update of the digitization parameter

Revision 1.6  2000/11/01 14:53:20  cblume
Merge with TRD-develop


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

  fInputFile   = NULL;
  fInputFileCreated   = kFALSE;
  fOutputFile  = NULL;
  fOutputFileCreated  = kFALSE;
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

  fInputFile   = NULL;
  fOutputFile  = NULL;
  fClusterTree = NULL;
  fEvent       = 0;
  fVerbose     = 0;
  fPar         = 0;

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

  if (fInputFile && fInputFileCreated) {
    fInputFile->Close();
    delete fInputFile;
  }

  if (fOutputFile && fOutputFileCreated) {
    fOutputFile->Close();
    delete fOutputFile;
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

  ((AliTRDclusterizer &) c).fInputFile   = NULL;
  ((AliTRDclusterizer &) c).fOutputFile  = NULL;
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

  OpenInput(name,nEvent);
  OpenOutput(name);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::Open(const Char_t *inname, const Char_t *outname
                              , Int_t nEvent)
{
  //
  // Opens the AliROOT file. Output and input are in different files
  //

  OpenInput(inname,nEvent);
  OpenOutput(outname);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::OpenOutput(const Char_t *name)
{
  //
  // Open the output file
  //

  TDirectory *savedir = NULL;

  if (!fInputFile) return kFALSE;

  if (strcmp(name,fInputFile->GetName()) != 0) {
    savedir = gDirectory;
    TFile *file = (TFile *) gROOT->FindObject(name);
    if (file) {
      fOutputFile = file;
      fOutputFile->cd();
    }
    else {
      printf("AliTRDclusterizer::OpenOutput -- ");
      printf("Open the output file %s.\n",name);
      fOutputFile = new TFile(name,"RECREATE");
      fOutputFileCreated = kTRUE;
    }
  }

  // Create a tree for the cluster
  Char_t treeName[12];
  sprintf(treeName,"TreeR%d_TRD",fEvent);
  fClusterTree = new TTree(treeName,"TRD cluster");
  TObjArray *ioArray = 0;
  fClusterTree->Branch("TRDcluster","TObjArray",&ioArray,32000,0);

  if (savedir) {
    savedir->cd();
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizer::OpenInput(const Char_t *name, Int_t nEvent)
{
  //
  // Opens a ROOT-file with TRD-hits and reads in the digits-tree
  //

  // Connect the AliRoot file containing Geometry, Kine, and Hits
  fInputFile = (TFile*) gROOT->GetListOfFiles()->FindObject(name);
  if (!fInputFile) {
    printf("AliTRDclusterizer::OpenInput -- ");
    printf("Open the ALIROOT-file %s.\n",name);
    fInputFile = new TFile(name,"UPDATE");
    fInputFileCreated = kTRUE;
  }
  else {
    printf("AliTRDclusterizer::OpenInput -- ");
    printf("%s is already open.\n",name);
  }

  // Get AliRun object from file
  if (!(gAlice)) gAlice = (AliRun *) fInputFile->Get("gAlice");
  if (!(gAlice)) {
    printf("AliTRDclusterizer::OpenInput -- ");
    printf("Could not find AliRun object.\n");
    return kFALSE;
  }

  fEvent = nEvent;

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(fEvent);
  if (nparticles <= 0) {
    printf("AliTRDclusterizer::OpenInput -- ");
    printf("No entries in the trees for event %d.\n",fEvent);
    return kFALSE;
  }

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
 
  TDirectory *savedir = gDirectory;

  if (fOutputFile) {
    fOutputFile->cd();
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

    return kTRUE;

  }

  if (det == -1) {

    printf("AliTRDclusterizer::WriteClusters -- ");
    printf("Writing the cluster tree %-18s for event %d.\n"
	  ,fClusterTree->GetName(),fEvent);

    fClusterTree->Write();

    AliTRDgeometry *geo = fTRD->GetGeometry();
    if (!fOutputFile->Get("TRDgeometry")) {
      geo->SetName("TRDgeometry");
      geo->Write();
    }
    if (!fOutputFile->Get("TRDparameter")) {
      fPar->Write();     
    }    

    return kTRUE;  

  }
  
  savedir->cd();

  printf("AliTRDclusterizer::WriteClusters -- ");
  printf("Unexpected detector index %d.\n",det);
 
  return kFALSE;  

}



