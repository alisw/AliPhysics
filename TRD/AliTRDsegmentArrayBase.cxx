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
Revision 1.7  2000/11/20 08:56:07  cblume
Cleanup of data arrays

Revision 1.6  2000/11/01 14:53:21  cblume
Merge with TRD-develop

Revision 1.1.4.3  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.4.2  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.5  2000/06/09 11:10:07  cblume
Compiler warnings and coding conventions, next round

Revision 1.4  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.3  2000/06/07 16:27:01  cblume
Try to remove compiler warnings on Sun and HP

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.1  2000/05/08 14:55:03  cblume
Bug fixes

Revision 1.1  2000/02/28 19:02:56  cblume
Add new TRD classes

*/

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

#include "AliTRDarrayI.h"
#include "AliTRDsegmentID.h"
#include "AliTRDsegmentArrayBase.h"

ClassImp(AliTRDsegmentArrayBase)
  
//_____________________________________________________________________________
AliTRDsegmentArrayBase::AliTRDsegmentArrayBase():TNamed()
{
  //
  // AliTRDsegmentArrayBase default constructor
  //

  fNSegment  = 0;
  fSegment   = 0; 
  fTreeIndex = 0;
  fTree      = 0;
  fClass     = 0;
  fBranch    = 0;

}

//_____________________________________________________________________________
AliTRDsegmentArrayBase::AliTRDsegmentArrayBase(Text_t *classname, Int_t n)
{
  //
  //  Create an array of objects of <classname>. The class must inherit from
  //  AliTRDsegmentID. The second argument sets the number of entries in 
  //  the array.
  //

  fNSegment  = 0;
  fSegment   = 0; 
  fTreeIndex = 0;
  fTree      = 0;
  fClass     = 0;
  fBranch    = 0;

  SetClass(classname);

  if (MakeArray(n) == kFALSE) {
     Error("AliTRDsegmentArrayBase","Cannot allocate %d segments in memory",n);
     return;
  }

}

//_____________________________________________________________________________
AliTRDsegmentArrayBase::AliTRDsegmentArrayBase(const AliTRDsegmentArrayBase &a)
{
  //
  // AliTRDsegmentArrayBase copy constructor
  //
  
  ((AliTRDsegmentArrayBase &) a).Copy(*this);

}

//_____________________________________________________________________________
AliTRDsegmentArrayBase::~AliTRDsegmentArrayBase()
{
  //
  // AliTRDsegmentArrayBase destructor
  //

  if (fNSegment){
    fSegment->Delete();
    delete fSegment;
  }

  if (fTree)      delete fTree;
  if (fTreeIndex) delete fTreeIndex;
  if (fClass)     delete fClass;

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
void AliTRDsegmentArrayBase::Copy(TObject &a)
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
Bool_t AliTRDsegmentArrayBase::SetClass(Text_t *classname)
{
  //
  // Sets the classname of the stored object
  //

  if (fClass   != 0) {
    delete fClass;
    fClass = 0;
  }
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

  if (!gROOT) ::Fatal("AliTRDsegmentArrayBase::AliTRDsegmentArrayBase"
                     ,"ROOT system not initialized");
   
   fClass = gROOT->GetClass(classname);
   if (!fClass) {
     Error("AliTRDsegmentArrayBase","%s is not a valid class name",classname);
     return kFALSE;
   }
   if (!fClass->InheritsFrom(AliTRDsegmentID::Class())) {
     Error("AliTRDsegmentArrayBase"
          ,"%s does not inherit from AliTRDsegmentID",classname);
     return kFALSE;
   }
  
   return kTRUE;

}

//_____________________________________________________________________________
AliTRDsegmentID * AliTRDsegmentArrayBase::NewSegment()
{
  //
  // Create a new object according to the class information
  //

  if (fClass  == 0) return 0;

  AliTRDsegmentID *segment = (AliTRDsegmentID *) fClass->New();
  if (segment == 0) return 0;

  return segment;

}

//_____________________________________________________________________________
Bool_t AliTRDsegmentArrayBase::AddSegment(AliTRDsegmentID *segment)
{
  //
  // Add a segment to the array
  //

  if (segment  == 0) return kFALSE;
  if (fSegment == 0) return kFALSE;
  if (fClass   == 0) return kFALSE;

  if (!(segment->IsA()->InheritsFrom(fClass))) {
    Error("AliTRDsegmentArrayBase","added class %s is not of proper type",
	  segment->IsA()->GetName());
    return kFALSE;
  }

  fSegment->AddAt(segment,segment->GetID());
  fNSegment = fSegment->GetLast() + 1;

  return kTRUE;

}

//_____________________________________________________________________________
AliTRDsegmentID * AliTRDsegmentArrayBase::AddSegment(Int_t index)
{
  //
  // Add a segment to the array
  //

  if (fSegment == 0) return 0;
  if (fClass   == 0) return 0;

  AliTRDsegmentID *segment = NewSegment();
  if (segment  == 0) return 0;

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
  if ((fSegment) && (fTreeIndex)) 
    return kTRUE;
  else 
    return kFALSE;
		  
}

//_____________________________________________________________________________
void AliTRDsegmentArrayBase::ClearSegment(Int_t index)
{
  //
  // Remove a segment from the active memory    
  //

  if ((*fSegment)[index]){
    delete (*fSegment)[index]; // because problem with deleting TClonesArray
    fSegment->RemoveAt(index);
  }

}

//_____________________________________________________________________________
void AliTRDsegmentArrayBase::MakeTree(char *file)
{
  //
  // Create a tree for the segment
  //

  AliTRDsegmentID *psegment = NewSegment();  

  if (fTree) delete fTree;
  fTree   = new TTree("Segment Tree","Tree with segments");

  fBranch = fTree->Branch("Segment",psegment->IsA()->GetName(),&psegment,64000,1);
  if (file) 
      fBranch->SetFile(file);      

  delete psegment;

}              

//_____________________________________________________________________________
Bool_t AliTRDsegmentArrayBase::ConnectTree(const char * treeName)
{
  //
  // Connect a tree from current directory  
  //

  if (fTree){
    delete fTree;
    fTree   = 0;
    fBranch = 0;
  }

  fTree   = (TTree*) gDirectory->Get(treeName);
  if (fTree   == 0) return kFALSE;
  fBranch = fTree->GetBranch("Segment");
  if (fBranch == 0) return kFALSE;

  MakeDictionary(TMath::Max(fNSegment,Int_t(fTree->GetEntries())));

  return kTRUE;

}

//_____________________________________________________________________________
AliTRDsegmentID *AliTRDsegmentArrayBase::LoadSegment(Int_t index)
{
  //
  // Load a segment with index <index> into the memory
  //

  if (fTreeIndex == 0) MakeDictionary(3000);

  // First try to load dictionary 
  if (fTreeIndex == 0)        return 0;
  if (fBranch    == 0)        return 0;
  if (index > fTreeIndex->fN) return 0;
  AliTRDsegmentID *s = (AliTRDsegmentID*) (*fSegment)[index];
  if (s == 0) s = NewSegment();
  s->SetID(index);
  
  if (s != 0) {
    Int_t treeIndex = (*fTreeIndex)[index];
    if (treeIndex < 1) 
      return 0;
    else 
      treeIndex--;   
    fBranch->SetAddress(&s);
    fTree->GetEvent(treeIndex);
    (*fSegment)[index] = (TObject*) s;
  }
  else 
    return 0;

  return s;

}

//_____________________________________________________________________________
AliTRDsegmentID *AliTRDsegmentArrayBase::LoadEntry(Int_t index)
{
  //
  // Load a segment at position <index> in the tree into the memory
  //

  if (fBranch == 0)                return 0;
  if (index > fTree->GetEntries()) return 0;

  AliTRDsegmentID *s = NewSegment();  
  if (s) {
    fBranch->SetAddress(&s);
    fTree->GetEvent(index);
  }
  else 
    return 0;

  Int_t nindex = s->GetID();
  ClearSegment(nindex);
  (*fSegment)[nindex] = (TObject *) s;

  return s;

}

//_____________________________________________________________________________
void AliTRDsegmentArrayBase::StoreSegment(Int_t index)
{
  //
  // Make a segment persistent 
  //

  const AliTRDsegmentID *kSegment = (*this)[index];
  if (kSegment == 0) return;
  if (fTree    == 0) MakeTree();
  fBranch->SetAddress(&kSegment);
  fTree->Fill();

}

//_____________________________________________________________________________
Bool_t  AliTRDsegmentArrayBase::MakeDictionary(Int_t size)
{
  //
  // Create an index table for the tree
  //  

  if (size < 1)   return kFALSE;
  if (fTreeIndex) delete fTreeIndex;

  fTreeIndex = new AliTRDarrayI(); 
  fTreeIndex->Set(size);
  
  AliTRDsegmentID   segment;
  AliTRDsegmentID *psegment = &segment;

  fBranch->SetAddress(&psegment);
  TBranch *brindix = fTree->GetBranch("fSegmentID");

  Int_t nevent = (Int_t) fTree->GetEntries();  
  for (Int_t i = 0; i < nevent; i++){
    brindix->GetEvent(i);
    Int_t treeIndex = segment.GetID();
    if (fTreeIndex->fN < treeIndex) 
      fTreeIndex->Expand(Int_t (Float_t(treeIndex) * 1.5) + 1);
    (*fTreeIndex)[treeIndex] = i + 1; 
  }

  return kTRUE;

}

//_____________________________________________________________________________
const AliTRDsegmentID * AliTRDsegmentArrayBase::operator[](Int_t i)
{
  //
  // Returns a segment with the given index <i>
  //

  if ((i < 0) || (i >= fNSegment)) return 0; 
  return (AliTRDsegmentID *) fSegment->At(i);

}

//_____________________________________________________________________________
const AliTRDsegmentID *AliTRDsegmentArrayBase::At(Int_t i) const
{
  //
  // Returns a segment with the given index <i>
  //

  if ((i < 0) || (i >= fNSegment)) return 0; 
  return (AliTRDsegmentID *)((*fSegment)[i]);

}
