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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Alice segment manager class                                           //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TDirectory.h>
#include <TError.h>
#include <TClass.h>

#include "AliLog.h"

#include "AliTRDgeometry.h"
#include "AliTRDsegmentArray.h"
#include "AliTRDsegmentID.h"
#include "AliTRDdataArray.h"
#include "AliTRDarrayI.h"

ClassImp(AliTRDsegmentArray)

//_____________________________________________________________________________
AliTRDsegmentArray::AliTRDsegmentArray()
  :TNamed()
  ,fSegment(0) 
  ,fTreeIndex(0)
  ,fNSegment(0)
  ,fTree(0)
  ,fBranch(0)
  ,fClass(0)
{
  //
   // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDsegmentArray::AliTRDsegmentArray(const char *classname, Int_t n)
  :TNamed()
  ,fSegment(0) 
  ,fTreeIndex(0)
  ,fNSegment(0)
  ,fTree(0)
  ,fBranch(0)
  ,fClass(0)
{
  //
  //  Create an array of objects of <classname>. The class must inherit from
  //  AliTRDsegmentID. The second argument sets the number of entries in 
  //  the array.
  //

  AliTRDdataArray *dataArray;  

  SetClass(classname);

  if (MakeArray(n) == kFALSE) {
    AliError(Form("Cannot allocate %d segments in memory",n));
    return;
  }

  for (Int_t i = 0; i < n; i++) {
    dataArray = (AliTRDdataArray *) AddSegment(i);
  }

}

//_____________________________________________________________________________
AliTRDsegmentArray::AliTRDsegmentArray(AliTRDsegmentArray &a)
  :TNamed(a)
  ,fSegment(a.fSegment) 
  ,fTreeIndex(a.fTreeIndex)
  ,fNSegment(a.fNSegment)
  ,fTree(a.fTree)
  ,fBranch(a.fBranch)
  ,fClass(a.fClass)
{
  //
  // AliTRDsegmentArray copy constructor
  //

  a.Copy(*this);

}

//_____________________________________________________________________________
AliTRDsegmentArray::~AliTRDsegmentArray()
{
  //
  // AliTRDsegmentArray destructor
  //

  Delete();

  if (fNSegment) {
    fSegment->Delete();
    delete fSegment;
  }

  if (fTreeIndex) {
    delete fTreeIndex;
  }

}

//_____________________________________________________________________________
AliTRDsegmentArray &AliTRDsegmentArray::operator=(const AliTRDsegmentArray &a)
{
  //
  // Assignment operator
  //

  if (this != &a) ((AliTRDsegmentArray &) a).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDsegmentArray::Copy(TObject &a) const
{
  //
  // Copy function
  //

  TNamed::Copy(a);

  fSegment->Copy(*((AliTRDsegmentArray &) a).fSegment);
  fTreeIndex->Copy(*((AliTRDsegmentArray &) a).fTreeIndex);
  fClass->Copy(*((AliTRDsegmentArray &) a).fClass);

  ((AliTRDsegmentArray &) a).fNSegment = fNSegment;

}

//_____________________________________________________________________________
Bool_t AliTRDsegmentArray::SetClass(const Char_t *classname)
{
  //
  // Sets the classname of the stored object
  //

  if (fTree    != 0) {
    delete fTree;
    fTree      = 0;
    fBranch    = 0;
    delete fTreeIndex;
    fTreeIndex = 0;
  } 
  if (fSegment != 0) {
    fSegment->Delete();
    delete fSegment;
    fSegment   = 0;
  }

  if (!gROOT) {
    AliFatal("ROOT system not initialized");
    exit(1);
  }   

  fClass = gROOT->GetClass(classname);
  if (!fClass) {
    AliError(Form("%s is not a valid class name",classname));
    return kFALSE;
  }
  if (!fClass->InheritsFrom(AliTRDsegmentID::Class())) {
    AliError(Form("%s does not inherit from AliTRDsegmentID",classname));
    return kFALSE;
  }
  
  return kTRUE;

}

//_____________________________________________________________________________
AliTRDsegmentID *AliTRDsegmentArray::NewSegment()
{
  //
  // Create a new object according to the class information
  //

  if (fClass  == 0) {
    return 0;
  }

  AliTRDsegmentID *segment = (AliTRDsegmentID *) fClass->New();

  if (segment == 0) {
    return 0;
  }
  else {
    return segment;
  }

}

//_____________________________________________________________________________
Bool_t AliTRDsegmentArray::AddSegment(AliTRDsegmentID *segment)
{
  //
  // Add a segment to the array
  //

  if (segment  == 0) {
    return kFALSE;
  }
  if (fSegment == 0) {
    return kFALSE;
  }
  if (fClass   == 0) {
    return kFALSE;
  }

  if (!(segment->IsA()->InheritsFrom(fClass))) {
    AliError(Form("added class %s is not of proper type"
                 ,segment->IsA()->GetName()));
    return kFALSE;
  }

  fSegment->AddAt(segment,segment->GetID());
  fNSegment = fSegment->GetLast() + 1;

  return kTRUE;

}

//_____________________________________________________________________________
AliTRDsegmentID *AliTRDsegmentArray::AddSegment(Int_t index)
{
  //
  // Add a segment to the array
  //

  if (fSegment == 0) {
    return 0;
  }
  if (fClass   == 0) {
    return 0;
  }

  AliTRDsegmentID *segment = NewSegment();
  if (segment  == 0) {
    return 0;
  }

  fSegment->AddAt(segment,index);
  segment->SetID(index);
  fNSegment = fSegment->GetLast() + 1;

  return segment;

}

//_____________________________________________________________________________
Bool_t AliTRDsegmentArray::MakeArray(Int_t n)
{
  //
  // Create an array of pointers to the segments
  //

  if (fSegment) {
    fSegment->Delete();
    delete fSegment;
  }
  if (fTreeIndex) delete fTreeIndex;  

  fSegment   = new TObjArray(n);
  fTreeIndex = new AliTRDarrayI();
  fTreeIndex->Set(n);
  fNSegment  = n;
  if ((fSegment) && (fTreeIndex)) {
    return kTRUE;
  }
  else { 
    return kFALSE;
  }
		  
}

//_____________________________________________________________________________
void AliTRDsegmentArray::ClearSegment(Int_t index)
{
  //
  // Remove a segment from the active memory    
  //

  if (fSegment->At(index)) {
    delete fSegment->RemoveAt(index);
  }

}

//_____________________________________________________________________________
void AliTRDsegmentArray::MakeTree(char *file)
{
  //
  // Create a tree for the segment
  //

  AliTRDsegmentID *psegment = NewSegment();  

  if (fTree) {
    delete fTree;
  }

  fTree   = new TTree("Segment Tree","Tree with segments");
  fBranch = fTree->Branch("Segment",psegment->IsA()->GetName(),&psegment,64000);

  if (file) {
    fBranch->SetFile(file);      
  }

  delete psegment;

}              

//_____________________________________________________________________________
Bool_t AliTRDsegmentArray::ConnectTree(const char *treeName)
{
  //
  // Connect a tree from current directory  
  //

  if (fTree) {
    delete fTree;
    fTree   = 0;
    fBranch = 0;
  }

  fTree = (TTree *) gDirectory->Get(treeName);
  if (fTree   == 0) {
    return kFALSE;
  }
  fBranch = fTree->GetBranch("Segment");
  if (fBranch == 0) {
    return kFALSE;
  }

  MakeDictionary(TMath::Max(fNSegment,Int_t(fTree->GetEntries())));

  return kTRUE;

}

//_____________________________________________________________________________
AliTRDsegmentID *AliTRDsegmentArray::LoadSegment(Int_t index)
{
  //
  // Load a segment with index <index> into the memory
  //

  if (fTreeIndex == 0) {
    MakeDictionary(3000);
  }

  // First try to load dictionary 
  if (fTreeIndex == 0) {
    return 0;
  }
  if (fBranch    == 0) {
    return 0;
  }
  if (index > fTreeIndex->fN) {
    return 0;
  }

  AliTRDsegmentID *s = (AliTRDsegmentID *) fSegment->At(index);
  if (s == 0) {
    s = NewSegment();
  }
  s->SetID(index);
  
  if (s != 0) {
    Int_t treeIndex = (*fTreeIndex)[index];
    if (treeIndex < 1) {
      return 0;
    }
    else { 
      treeIndex--;
    }   
    fBranch->SetAddress(&s);
    fTree->GetEvent(treeIndex);
    fSegment->AddAt((TObject*) s, index);
  }
  else { 
    return 0;
  }

  return s;

}

//_____________________________________________________________________________
AliTRDsegmentID *AliTRDsegmentArray::LoadEntry(Int_t index)
{
  //
  // Load a segment at position <index> in the tree into the memory
  //

  if (fBranch == 0) {
    return 0;
  }
  if (index > fTree->GetEntries()) {
    return 0;
  }

  AliTRDsegmentID *s = NewSegment();  
  if (s) {
    fBranch->SetAddress(&s);
    fTree->GetEvent(index);
  }
  else {
    return 0;
  }

  Int_t nindex = s->GetID();
  ClearSegment(nindex);
  fSegment->AddAt((TObject *) s, nindex);

  return s;

}

//_____________________________________________________________________________
void AliTRDsegmentArray::StoreSegment(Int_t index)
{
  //
  // Make a segment persistent 
  //

  const AliTRDsegmentID *kSegment = (*this)[index];
  if (kSegment == 0) {
    return;
  }
  if (fTree    == 0) {
    MakeTree();
  }
  fBranch->SetAddress(&kSegment);
  fTree->Fill();

}

//_____________________________________________________________________________
Bool_t  AliTRDsegmentArray::MakeDictionary(Int_t size)
{
  //
  // Create an index table for the tree
  //  

  if (size < 1) {
    return kFALSE;
  }
  if (fTreeIndex) {
    delete fTreeIndex;
  }

  fTreeIndex = new AliTRDarrayI(); 
  fTreeIndex->Set(size);
  
  AliTRDsegmentID   segment;
  AliTRDsegmentID *psegment = &segment;

  fBranch->SetAddress(&psegment);
  TBranch *brindix = fTree->GetBranch("fSegmentID");

  Int_t nevent = (Int_t) fTree->GetEntries();  
  for (Int_t i = 0; i < nevent; i++) {
    brindix->GetEvent(i);
    Int_t treeIndex = segment.GetID();
    if (fTreeIndex->fN < treeIndex) {
      fTreeIndex->Expand(Int_t (Float_t(treeIndex) * 1.5) + 1);
    }
    (*fTreeIndex)[treeIndex] = i + 1; 
  }

  return kTRUE;

}

//_____________________________________________________________________________
const AliTRDsegmentID * AliTRDsegmentArray::operator[](Int_t i) const
{
  //
  // Returns a segment with the given index <i>
  //

  if ((i <          0) || 
      (i >= fNSegment)) {
    return 0; 
  }

  return (AliTRDsegmentID *) fSegment->At(i);

}

//_____________________________________________________________________________
const AliTRDsegmentID *AliTRDsegmentArray::At(Int_t i) const
{
  //
  // Returns a segment with the given index <i>
  //

  if ((i <          0) || 
      (i >= fNSegment)) {
    return 0; 
  }

  return (AliTRDsegmentID *) fSegment->At(i);

}

//_____________________________________________________________________________
void AliTRDsegmentArray::Delete()
{
  //
  // Deletes all detector segments from the array
  //

  for (Int_t iDet = 0; iDet < fNSegment; iDet++) {
    ClearSegment(iDet);
  }

}

//_____________________________________________________________________________
Bool_t AliTRDsegmentArray::LoadArray(const Char_t *branchname, TTree *tree)
{
  //
  // Loads all segments of the array from the branch <branchname> of
  // the digits tree <tree>
  //

  fTree = tree;

  if (!fTree) {
    AliError("Digits tree is not defined\n");
    return kFALSE;
  }

  // Get the branch
  fBranch = fTree->GetBranch(branchname);
  if (!fBranch) {
    AliError(Form("Branch %s is not defined\n",branchname));
    return kFALSE;
  }

  // Loop through all segments and read them from the tree
  Bool_t status = kTRUE;
  for (Int_t iSegment = 0; iSegment < fNSegment; iSegment++) {
    AliTRDdataArray *dataArray = (AliTRDdataArray *) fSegment->At(iSegment);
    if (!dataArray) {
      status = kFALSE;
      break;    
    }
    fBranch->SetAddress(&dataArray);
    fBranch->GetEntry(iSegment);
  }

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDsegmentArray::StoreArray(const Char_t *branchname, TTree *tree)
{
  //
  // Stores all segments of the array in the branch <branchname> of 
  // the digits tree <tree>
  //

  fTree = tree;

  if (!fTree) {
    AliError("Digits tree is not defined\n");
    return kFALSE;
  }

  // Get the branch
  fBranch = fTree->GetBranch(branchname);
  if (!fBranch) {
    AliError(Form("Branch %s is not defined\n",branchname));
    return kFALSE;
  }

  // Loop through all segments and fill them into the tree
  Bool_t status = kTRUE;
  for (Int_t iSegment = 0; iSegment < fNSegment; iSegment++) {
    const AliTRDdataArray *kDataArray = 
         (AliTRDdataArray *) AliTRDsegmentArray::At(iSegment);
    if (!kDataArray) {
      status = kFALSE;
      break;
    }
    fBranch->SetAddress(&kDataArray);
    fBranch->Fill();
  }

  return status;

}

//_____________________________________________________________________________
AliTRDdataArray *AliTRDsegmentArray::GetDataArray(Int_t det) const
{
  //
  // Returns the data array for a given detector
  //

  return ((AliTRDdataArray *) AliTRDsegmentArray::At(det));

}

//_____________________________________________________________________________
AliTRDdataArray *AliTRDsegmentArray::GetDataArray(Int_t pla
                                                , Int_t cha
                                                , Int_t sec) const
{
  //
  // Returns the data array for a given detector
  //

  Int_t det = AliTRDgeometry::GetDetector(pla,cha,sec);
  return GetDataArray(det);

}
