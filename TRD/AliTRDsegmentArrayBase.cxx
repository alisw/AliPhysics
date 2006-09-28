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
//  Alice segment manager base class                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TDirectory.h>
#include <TError.h>
#include <TClass.h>

#include "AliLog.h"

#include "AliTRDsegmentArrayBase.h"
#include "AliTRDarrayI.h"
#include "AliTRDsegmentID.h"

ClassImp(AliTRDsegmentArrayBase)
  
//_____________________________________________________________________________
AliTRDsegmentArrayBase::AliTRDsegmentArrayBase()
  :TNamed()
  ,fSegment(0) 
  ,fTreeIndex(0)
  ,fNSegment(0)
  ,fTree(0)
  ,fBranch(0)
  ,fClass(0)
{
  //
  // AliTRDsegmentArrayBase default constructor
  //

}

//_____________________________________________________________________________
AliTRDsegmentArrayBase::AliTRDsegmentArrayBase(const char *classname, Int_t n)
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

  SetClass(classname);

  if (MakeArray(n) == kFALSE) {
    Error("AliTRDsegmentArrayBase","Cannot allocate %d segments in memory",n);
    return;
  }

}

//_____________________________________________________________________________
AliTRDsegmentArrayBase::AliTRDsegmentArrayBase(const AliTRDsegmentArrayBase &a)
  :TNamed(a)
  ,fSegment(a.fSegment) 
  ,fTreeIndex(a.fTreeIndex)
  ,fNSegment(a.fNSegment)
  ,fTree(a.fTree)
  ,fBranch(a.fBranch)
  ,fClass(a.fClass)
{
  //
  // AliTRDsegmentArrayBase copy constructor
  //

}

//_____________________________________________________________________________
AliTRDsegmentArrayBase::~AliTRDsegmentArrayBase()
{
  //
  // AliTRDsegmentArrayBase destructor
  //

  if (fNSegment) {
    fSegment->Delete();
    delete fSegment;
  }

  if (fTreeIndex) {
    delete fTreeIndex;
  }

}

//_____________________________________________________________________________
AliTRDsegmentArrayBase &AliTRDsegmentArrayBase
                        ::operator=(const AliTRDsegmentArrayBase &a)
{
  //
  // Assignment operator
  //

  if (this != &a) ((AliTRDsegmentArrayBase &) a).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDsegmentArrayBase::Copy(TObject &a) const
{
  //
  // Copy function
  //

  TNamed::Copy(a);

  fSegment->Copy(*((AliTRDsegmentArrayBase &) a).fSegment);
  fTreeIndex->Copy(*((AliTRDsegmentArrayBase &) a).fTreeIndex);
  fClass->Copy(*((AliTRDsegmentArrayBase &) a).fClass);

  ((AliTRDsegmentArrayBase &) a).fNSegment = fNSegment;

}

//_____________________________________________________________________________
Bool_t AliTRDsegmentArrayBase::SetClass(const char *classname)
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
AliTRDsegmentID *AliTRDsegmentArrayBase::NewSegment()
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
Bool_t AliTRDsegmentArrayBase::AddSegment(AliTRDsegmentID *segment)
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
AliTRDsegmentID *AliTRDsegmentArrayBase::AddSegment(Int_t index)
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
Bool_t AliTRDsegmentArrayBase::MakeArray(Int_t n)
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
void AliTRDsegmentArrayBase::ClearSegment(Int_t index)
{
  //
  // Remove a segment from the active memory    
  //

  if (fSegment->At(index)) {
    delete fSegment->RemoveAt(index);
  }

}

//_____________________________________________________________________________
void AliTRDsegmentArrayBase::MakeTree(char *file)
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
Bool_t AliTRDsegmentArrayBase::ConnectTree(const char *treeName)
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
AliTRDsegmentID *AliTRDsegmentArrayBase::LoadSegment(Int_t index)
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
AliTRDsegmentID *AliTRDsegmentArrayBase::LoadEntry(Int_t index)
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
void AliTRDsegmentArrayBase::StoreSegment(Int_t index)
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
Bool_t  AliTRDsegmentArrayBase::MakeDictionary(Int_t size)
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
const AliTRDsegmentID * AliTRDsegmentArrayBase::operator[](Int_t i) const
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
const AliTRDsegmentID *AliTRDsegmentArrayBase::At(Int_t i) const
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
