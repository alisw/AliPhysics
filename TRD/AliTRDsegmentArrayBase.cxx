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
Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.1  2000/05/08 14:55:03  cblume
Bug fixes

Revision 1.1  2000/02/28 19:02:56  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Alice segment manager object                                             //
//                                                                           //
//  AliTRDsegmentIDArray object  is array of pointers to object derived from //
//  AliTRDsegmentID object                                                   //
//  AliTRDsegmentID - object in comparison with TObject enhalt               //
//  additional information fSegmentID                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include  <TROOT.h>
#include <TTree.h>
#include "TClonesArray.h"
#include "TDirectory.h"
#include "AliTRDarrayI.h"
#include "TError.h"
#include "TClass.h"

#include "AliTRDsegmentID.h"
#include "AliTRDsegmentArrayBase.h"

ClassImp(AliTRDsegmentArrayBase)
  
//_____________________________________________________________________________
AliTRDsegmentArrayBase::AliTRDsegmentArrayBase()
{
  //
  //
  //

  fNSegment  = 0;
  fSegment   = 0; 
  fTreeIndex = 0;
  fTree      = 0;
  fClass     = 0;

}

//_____________________________________________________________________________
AliTRDsegmentArrayBase::AliTRDsegmentArrayBase(Text_t *classname, Int_t n)
{
  //
  //constructor which 
  // 
  //  Create an array of objects of classname. The class must inherit from
  //  AliTRDsegmentID .  The second argument adjust number of entries in 
  //  the array.
  fNSegment=0;
  fSegment =0; 
  fTreeIndex = 0;
  fTree  = 0;
  fClass = 0;
  SetClass(classname);
  if (MakeArray(n)==kFALSE){
     Error("AliTRDsegmentArrayBase", "can't allocate %d segments in memory",n);
     return;
   }
}

//_____________________________________________________________________________
Bool_t AliTRDsegmentArrayBase:: SetClass(Text_t *classname)
{
  //
  //set class of stored object
  if ( fClass !=0 ) {
    delete fClass;
    fClass = 0;
  }
  if (fTree !=0) {
    delete fTree;
    fTree = 0;
    fBranch = 0;
    delete fTreeIndex;
    fTreeIndex = 0;
  } 
  if (fSegment != 0) {
    fSegment->Delete();
    delete fSegment;
    fSegment = 0;
  }
  if (!gROOT)
      ::Fatal("AliTRDsegmentArrayBase::AliTRDsegmentArrayBase", "ROOT system not initialized");
   
   fClass = gROOT->GetClass(classname);
   if (!fClass) {
      Error("AliTRDsegmentArrayBase", "%s is not a valid class name", classname);
      return kFALSE;
   }
   if (!fClass->InheritsFrom(AliTRDsegmentID::Class())) {
      Error("AliTRDsegmentArrayBase", "%s does not inherit from AliTRDsegmentID", classname);
      return kFALSE;
   }  
   return kTRUE;
}

//_____________________________________________________________________________
//Bool_t AliTRDsegmentArrayBase::ClassError( )
//{
  //signalize class error 
  //  if (!fClass) {
  //    Error("AliTRDsegmentArrayBase", "%s is not a valid class name", classname);
  //    return kFALSE;
  // }
////  return kFALSE;
//}

//_____________________________________________________________________________
AliTRDsegmentArrayBase::~AliTRDsegmentArrayBase()
{
  if (fNSegment>0){
    fSegment->Delete();
    delete fSegment;
  }
  if (fTree) delete fTree;
  if (fTreeIndex) delete fTreeIndex;
  if (fClass!=0) delete fClass;
}

//_____________________________________________________________________________
AliTRDsegmentID * AliTRDsegmentArrayBase::NewSegment()
{
  //
  //create object according class information
  if (fClass==0) return 0;
  AliTRDsegmentID * segment = (AliTRDsegmentID * )fClass->New();
  if (segment == 0) return 0;
  return segment;
}

//_____________________________________________________________________________
Bool_t AliTRDsegmentArrayBase::AddSegment(AliTRDsegmentID *segment)
{
  //
  // add segment to array
  //
  if (segment==0) return kFALSE;
  if (fSegment==0) return kFALSE;
  if (fClass==0) return kFALSE;
  if (!(segment->IsA()->InheritsFrom(fClass))){
    Error("AliTRDsegmentArrayBase", "added class %s  is not of proper type ",
	  segment->IsA()->GetName());
      return kFALSE;
  }
  fSegment->AddAt(segment,segment->GetID());
  fNSegment = fSegment->GetLast()+1;
  return kTRUE;
}

//_____________________________________________________________________________
AliTRDsegmentID * AliTRDsegmentArrayBase::AddSegment(Int_t index)
{
  //
  // add segment to array
  //
  if (fSegment==0) return 0;
  if (fClass==0) return 0;
  //  AliTRDsegmentID * segment = (AliTRDsegmentID * )fClass->New();
  AliTRDsegmentID * segment = NewSegment();
  if (segment == 0) return 0;
  fSegment->AddAt(segment,index);
  segment->SetID(index);
  fNSegment = fSegment->GetLast()+1;
  return segment;
}

//_____________________________________________________________________________
Bool_t AliTRDsegmentArrayBase::MakeArray(Int_t n)
{
  //
  //make array of pointers to Segments
  //
  if (fSegment) {
    fSegment->Delete();
    delete fSegment;
  }
  if (fTreeIndex) delete   fTreeIndex;  
  fSegment = new TObjArray(n);
  fTreeIndex = new AliTRDarrayI;
  fTreeIndex->Set(n);
  fNSegment=n;
  if ( (fSegment) && (fTreeIndex)) return kTRUE;
  else return kFALSE;		  
}

//_____________________________________________________________________________
void AliTRDsegmentArrayBase::ClearSegment(Int_t index)
{
  //
  //remove segment from active memory    
  //
  if ((*fSegment)[index]){
    //    (*fSegment)[index]->Delete(); //not working for TClonesArray
    delete (*fSegment)[index]; //because problem with deleting TClonesArray
    fSegment->RemoveAt(index);
  }
}

//_____________________________________________________________________________
void AliTRDsegmentArrayBase::MakeTree()
{
  //  AliTRDsegmentID  segment;
  AliTRDsegmentID * psegment = NewSegment();  
  if (fTree) delete fTree;
  fTree = new TTree("Segment Tree","Tree with segments");
  fBranch = fTree->Branch("Segment",psegment->IsA()->GetName(),&psegment,64000,1);
  delete psegment;
}              

//_____________________________________________________________________________
Bool_t AliTRDsegmentArrayBase::ConnectTree(const char * treeName)
{
  //connect tree from current directory  
  if (fTree){
    delete fTree;
    fTree = 0;
    fBranch = 0;
  }
  fTree =(TTree*)gDirectory->Get(treeName);
  if (fTree == 0)    return kFALSE;
  fBranch = fTree->GetBranch("Segment");
  if (fBranch==0) return kFALSE;
  MakeDictionary(TMath::Max(fNSegment,Int_t(fTree->GetEntries())));
  return kTRUE;
}

//_____________________________________________________________________________
AliTRDsegmentID *AliTRDsegmentArrayBase::LoadSegment(Int_t index)
{
  //
  //load segment with index to the memory
  //
  //
  if (fTreeIndex ==0 ) MakeDictionary(3000);
  //firstly try to load dictionary 
  if (fTreeIndex ==0 ) return 0;
  if (fBranch==0) return 0;
  if (index>fTreeIndex->fN) return 0;
  AliTRDsegmentID *s = (AliTRDsegmentID*)(*fSegment)[index];
  if (s==0)  s=  NewSegment();
  s->SetID(index);
  //  new AliTRDsegmentID(index);
  
  if (s!=0) {
    Int_t treeIndex =(*fTreeIndex)[index];
    if (treeIndex<1) return 0;
    else treeIndex--;   //I don't like it Int table I have index shifted by 1		       
    fBranch->SetAddress(&s);
    fTree->GetEvent(treeIndex);
    (*fSegment)[index] = (TObject*) s;
  }
  else 
    return 0;
  return s;
  //  AbstractMethod("LoadSegment");
}

//_____________________________________________________________________________
AliTRDsegmentID *AliTRDsegmentArrayBase::LoadEntry(Int_t index)
{
  //
  //load segment at position inex in tree  to the memory
  //
  //
  if (fBranch==0) return 0;
  if (index>fTree->GetEntries()) return 0;
  AliTRDsegmentID * s =  NewSegment();
  
  if (s) {
    fBranch->SetAddress(&s);
    fTree->GetEvent(index);
  }
  else 
    return 0;
  Int_t nindex = s->GetID();
  ClearSegment(nindex);
  (*fSegment)[nindex] = (TObject*) s;
  return s;
  //  AbstractMethod("LoadSegment");
}

void AliTRDsegmentArrayBase::StoreSegment(Int_t index)
{
  //
  //make segment persistent 
  //
  const AliTRDsegmentID *  segment = (*this)[index];
  if (segment == 0 ) return;
  if (fTree==0) MakeTree();
  fBranch->SetAddress(&segment);
  fTree->Fill();
}

//_____________________________________________________________________________
Bool_t  AliTRDsegmentArrayBase::MakeDictionary(Int_t size)
{
  //
  //create index table for tree
  //  
  if (size<1) return kFALSE;
  if (fTreeIndex) delete fTreeIndex;
  fTreeIndex = new AliTRDarrayI(); 
  fTreeIndex->Set(size);
  
  AliTRDsegmentID  segment;
  AliTRDsegmentID * psegment = &segment;
  fBranch->SetAddress(&psegment);
  TBranch * brindix = fTree->GetBranch("fSegmentID");
  Int_t nevent = (Int_t)fTree->GetEntries();  
  for (Int_t i = 0; i<nevent; i++){
    brindix->GetEvent(i);
    Int_t treeIndex=segment.GetID();
    if (fTreeIndex->fN<treeIndex) fTreeIndex->Expand(Int_t(Float_t(treeIndex)*1.5)+1);
    //    Int_t index = segment.GetID(); 
    (*fTreeIndex)[treeIndex]=i+1; // MI 19.5. I'm sorry  -index 0 couldn't be use in AliTRDarrayI   
  }
  return kTRUE;
}

//_____________________________________________________________________________
const AliTRDsegmentID*  AliTRDsegmentArrayBase::operator[](Int_t i)
{
  //
  //return segment with given index
  //
  if ( (i<0) || (i>=fNSegment)) return 0; 
  return (AliTRDsegmentID *)fSegment->At(i);
}

//_____________________________________________________________________________
const AliTRDsegmentID*  AliTRDsegmentArrayBase::At(Int_t i)
{
  //
  //return segment with given index
  //
  if ( (i<0) || (i>=fNSegment)) return 0; 
  return (AliTRDsegmentID *)((*fSegment)[i]);
}
