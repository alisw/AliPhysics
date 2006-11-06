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
//  Alice segment manager object                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include <Riostream.h>

#include <TTree.h>
#include <TROOT.h>
#include "AliRun.h"

#include "TClonesArray.h"
#include "TDirectory.h"
#include <TArrayI.h>
#include "TError.h"
#include "TClass.h"
#include "TFile.h"

#include "AliSegmentID.h"
#include "AliSegmentArray.h"
#include "TObjString.h"


//_____________________________________________________________________________
ClassImp(AliSegmentArray)
  
  AliSegmentArray::AliSegmentArray()
                  :TNamed(),
                   fSegment(0),
                   fTreeIndex(0),
                   fNSegment(0),
                   fTree(0),
                   fTreeOwner(kFALSE),
                   fBranch(0),
                   fClass(0)                                      
{
  //
  // constructor
  //

}

AliSegmentArray::AliSegmentArray(const char *classname, Int_t n)
                :TNamed("SegmentArray","SegmentArray"),
                 fSegment(0),
                 fTreeIndex(0),
                 fNSegment(0),
                 fTree(0),
                 fTreeOwner(kFALSE),
                 fBranch(0),
                 fClass(0) 
{
  //
  //constructor which 
  // 
  //  Create an array of objects of classname. The class must inherit from
  //  AliSegmentID .  The second argument adjust number of entries in 
  //  the array.
 

  SetClass(classname);
  if (MakeArray(n)==kFALSE){
     Error("AliSegmentArray", "can't allocate %d segments in memory",n);
     return;
   }
}

AliSegmentArray::AliSegmentArray(const AliSegmentArray &segment)
                :TNamed(segment),
                 fSegment(0),
                 fTreeIndex(0),
                 fNSegment(0),
                 fTree(0),
                 fTreeOwner(kFALSE),
                 fBranch(0),
                 fClass(0)                                      
               
{
  //
  //copy constructor
  // to be later implemented
}

AliSegmentArray &AliSegmentArray::operator = (const AliSegmentArray & /*segment*/)
{
  //assignment operator
  //to be later implemented
  return (*this);
}

AliSegmentArray::~AliSegmentArray()
{
  //
  // default destructor
  if (fNSegment>0){
    fSegment->Delete();
    delete fSegment;
  }
  if (fTree) 
   if (fTreeOwner) delete fTree;
   else fTree->Reset();

  if (fTreeIndex) delete fTreeIndex;
  //  if (fClass!=0) delete fClass;
}


Bool_t AliSegmentArray::SetClass(const char *classname)
{
  //
  //set class of stored object
  if ( fClass !=0 ) {
    //delete fClass; not ower of fClass
    fClass = 0;
  }
  if (fTree !=0) {
    if (fTreeOwner) delete fTree;
    else fTree->Reset();
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
      ::Fatal("AliSegmentArray::AliSegmentArray", "ROOT system not initialized");
   
   fClass = gROOT->GetClass(classname);
   if (!fClass) {
      Error("AliSegmentArray", "%s is not a valid class name", classname);
      return kFALSE;
   }
   if (!fClass->InheritsFrom(AliSegmentID::Class())) {
      Error("AliSegmentArray", "%s does not inherit from AliSegmentID", classname);
      return kFALSE;
   }  
   return kTRUE;
}


AliSegmentID * AliSegmentArray::NewSegment()
{
  //
  //create object according class information
  if (fClass==0) return 0;
  AliSegmentID * segment = (AliSegmentID * )fClass->New();
  if (segment == 0) return 0;
  return segment;
}


Bool_t AliSegmentArray::AddSegment(AliSegmentID *segment)
{
  //
  // add segment to array
  //
  if (segment==0) return kFALSE;
  if (fSegment==0) return kFALSE;
  if (fClass==0) return kFALSE;
  if (!(segment->IsA()->InheritsFrom(fClass))){
    Error("AliSegmentArray", "added class %s  is not of proper type ",
	  segment->IsA()->GetName());
      return kFALSE;
  }
  fSegment->AddAt(segment,segment->GetID());
  fNSegment = fSegment->GetLast()+1;
  return kTRUE;
}

AliSegmentID * AliSegmentArray::AddSegment(Int_t index)
{
  //
  // add segment to array
  //
  if (fSegment==0) return 0;
  if (fClass==0) return 0;
  AliSegmentID * segment = NewSegment();
  if (segment == 0) return 0;
  fSegment->AddAt(segment,index);
  segment->SetID(index);
  fNSegment = fSegment->GetLast()+1;
  return segment;
}


void AliSegmentArray::ClearSegment(Int_t index)
{
  //
  //remove segment from active memory    
  //
  //PH  if ((*fSegment)[index]){
  if (fSegment->At(index)){
    //    (*fSegment)[index]->Delete(); //not working for TClonesArray
    //PH    delete (*fSegment)[index]; //because problem with deleting TClonesArray
    //PH    fSegment->RemoveAt(index);
    delete fSegment->RemoveAt(index);
  }
}


Bool_t AliSegmentArray::MakeArray(Int_t n)
{
  //
  //make array of pointers to Segments
  //
  if (fSegment) {
    fSegment->Delete();
    delete fSegment;
  }  
  fSegment = new TObjArray(n);  
  fNSegment=n;
  if (fSegment) return kTRUE;  
  else return kFALSE;		  
}
void AliSegmentArray::MakeTree(TTree* tree)
{
             //Make tree with the name
  AliSegmentID * psegment = NewSegment();  
  fTree = tree;
  //PH  fBranch = fTree->Branch("Segment",psegment->IsA()->GetName(),&psegment,64000);
  fBranch = fTree->Branch("Segment",psegment->IsA()->GetName(),&psegment,64000,99);

}

void AliSegmentArray::MakeTree(char *file)
{
  //  AliSegmentID  segment;
  AliSegmentID * psegment = NewSegment();  
  if (fTree) 
    if (fTreeOwner) 
     {
       delete fTree;
       fTree = new TTree("Segment Tree","Tree with segments");     
     }
    else fTree->Reset();

  
  //PH  fBranch = fTree->Branch("Segment",psegment->IsA()->GetName(),&psegment,64000);
   fBranch = fTree->Branch("Segment",psegment->IsA()->GetName(),&psegment,64000,99);
 
  if (file) {
        TString outFile = gAlice->GetBaseFile();
        outFile = outFile + "/" + file;
        fBranch->SetFile(outFile.Data());
        TDirectory *wd = gDirectory;
        TBranch *b = fBranch;
        TIter next( b->GetListOfBranches());
        while ((b=(TBranch*)next())) {
           b->SetFile(outFile.Data());
        }
        cout << "Diverting branch " << "Segment" << " to file " << outFile << endl;  
        wd->cd(); 
    }
  delete psegment;
}              


Bool_t  AliSegmentArray::MakeDictionary(Int_t size)
{
  //
  //create index table for tree
  //  
  if (size<1) return kFALSE;
  if (fTreeIndex) delete fTreeIndex;
  fTreeIndex = new TArrayI(); 
  fTreeIndex->Set(size);
  
  AliSegmentID * psegment = NewSegment(); //MI change
  fBranch->SetAddress(&psegment);
  TBranch * brindix = fTree->GetBranch("fSegmentID");
  Int_t nevent = (Int_t)fTree->GetEntries();  
  for (Int_t i = 0; i<nevent; i++){
    brindix->GetEvent(i);
    Int_t treeIndex=psegment->GetID();
    if (fTreeIndex->fN<treeIndex) fTreeIndex->Set(Int_t(Float_t(treeIndex)*1.5)+1);
    //    Int_t index = segment.GetID(); 
    (*fTreeIndex)[treeIndex]=i+1; //  
  }
  if (psegment) delete psegment;
  return kTRUE;
}

Bool_t AliSegmentArray::ConnectTree(TTree* tree)
{
  fTree =tree;
  if (fTree == 0)    return kFALSE;
  fBranch = fTree->GetBranch("Segment");
  if (fBranch==0) return kFALSE;
  MakeDictionary(TMath::Max(fNSegment,Int_t(fTree->GetEntries())));
  MakeArray(fTreeIndex->fN);
  return kTRUE;
}


Bool_t AliSegmentArray::ConnectTree(const char * treeName)
{
  //connect tree from current directory  
  if (fTree){
   if (fTreeOwner) 
    {
     delete fTree;
     fTree = 0;
    }
   else fTree->Reset();
   fBranch = 0;
  }
  fTree =(TTree*)gDirectory->Get(treeName);
  
  if (fTree == 0)    return kFALSE;
  fBranch = fTree->GetBranch("Segment");
  if (fBranch==0) return kFALSE;
  MakeDictionary(TMath::Max(fNSegment,Int_t(fTree->GetEntries())));
  MakeArray(fTreeIndex->fN);
  return kTRUE;
}

AliSegmentID *AliSegmentArray::LoadSegment(Int_t index)
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
  //PH  AliSegmentID *s = (AliSegmentID*)(*fSegment)[index];
  AliSegmentID *s = (AliSegmentID*)fSegment->At(index);
  if (s==0)  s=  NewSegment();
  s->SetID(index);
  //  new AliSegmentID(index);
  
  if (s!=0) {
    Int_t treeIndex =(*fTreeIndex)[index];
    if (treeIndex<1) return 0;
    else treeIndex--;   //I don't like it Int table I have index shifted by 1		       
    fBranch->SetAddress(&s);
    fTree->GetEvent(treeIndex);
    //PH    (*fSegment)[index] = (TObject*) s;
    fSegment->AddAt((TObject*) s, index);
  }
  else 
    return 0;
  return s;

}
AliSegmentID *AliSegmentArray::LoadEntry(Int_t index)
{
  //
  //load segment at position inex in tree  to the memory
  //
  //
  if (fBranch==0) return 0;
  if (index>fTree->GetEntries()) return 0;
  AliSegmentID * s =  NewSegment();
  
  if (s) {
    fBranch->SetAddress(&s);
    fTree->GetEvent(index);
  }
  else 
    return 0;
  Int_t nindex = s->GetID();
  ClearSegment(nindex);
  //PH  (*fSegment)[nindex] = (TObject*) s;
  fSegment->AddAt((TObject*) s, nindex);
  return s;
}

void AliSegmentArray::StoreSegment(Int_t index)
{
  //
  //make segment persistent 
  //
  const AliSegmentID *  ksegment = (*this)[index];
  if (ksegment == 0 ) return;
  if (fTree==0) MakeTree();
  fBranch->SetAddress(&ksegment);
  fTree->Fill();
}


void AliSegmentArray::Streamer(TBuffer &R__b)
{
  TObjString treeName, * ptreeName=&treeName;
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TNamed::Streamer(R__b);
    R__b>>ptreeName;
    if (fTree && fTreeOwner) delete fTree;
    ConnectTree(ptreeName->String());   
  } else {
    R__b.WriteVersion(AliSegmentArray::IsA());
    TNamed::Streamer(R__b);      
    //  char  ch[200];
    //  sprintf(ch,"%s",fTrre->GetTitle());
    treeName.String() = fTree->GetTitle();
    R__b<<ptreeName;
    fTree->Write();
  }
}

