
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFast virtual base class for Makers                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TTree.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TBrowser.h>

#include "AliFMaker.h"
#include "AliFast.h"

ClassImp(AliFMaker)

//_____________________________________________________________________________
AliFMaker::AliFMaker()
{
   fBranchName = "";
   fSave       = 0;
   fHistograms = 0;
   fFruits     = 0;
   fClones     = 0;
   fIsClonable = kTRUE;
}

//_____________________________________________________________________________
AliFMaker::AliFMaker(const char *name, const char *title)
       :TNamed(name,title)
{
   fBranchName = "";
   fSave       = 0;
   fHistograms = new TList();
   fClones     = 0;
   fIsClonable = kTRUE;
   
   gAliFast->Makers()->Add(this);
}

//_____________________________________________________________________________
AliFMaker::~AliFMaker()
{
  delete fFruits;
  delete fClones;
}

//______________________________________________________________________________
void AliFMaker::Browse(TBrowser *b)
{
//  Insert Maker objects in the list of objects to browsed.

  char name[64];
  if( b == 0  || fFruits == 0) return;
  TObject *obj;

// If fFruits is a ClonesArray, insert all the objects in the list
// of browsable objects
  if (fFruits->InheritsFrom("TClonesArray")) {
     TClonesArray *clones = (TClonesArray*)fFruits;
     Int_t nobjects = clones->GetEntries();
     for (Int_t i=0;i<nobjects;i++) {
        obj = clones->At(i);
        sprintf(name,"%s_%d",obj->GetName(),i);
        if (strstr(name,"AliF")) b->Add(obj, &name[4]);
        else                     b->Add(obj, &name[0]);
     }
// fFruits points to an object in general. Insert this object in the browser
  } else {
      b->Add( fFruits, fFruits->GetName());
  }
}

//_____________________________________________________________________________
void AliFMaker::Clear(Option_t *option)
{
  if (fFruits) fFruits->Clear(option);
  delete fClones;
  fClones = 0;
}

//_____________________________________________________________________________
void AliFMaker::Draw(Option_t *)
{
//    Insert products of this maker in graphics pad list

  TObject *obj;

// If fFruits is a ClonesArray, insert all the objects in the list
// of objects to be painted
  if (fFruits->InheritsFrom("TClonesArray")) {
     TClonesArray *clones = (TClonesArray*)fFruits;
     Int_t nobjects = clones->GetEntries();
     for (Int_t i=0;i<nobjects;i++) {
        obj = clones->At(i);
        if (obj) obj->AppendPad();
     }
// fFruits points to an object in general. Insert this object in the pad
  } else {
     fFruits->AppendPad();
  }
}

//_____________________________________________________________________________
void AliFMaker::FillClone()
{
//   Copy original fruits in a separate list (clones)

   if (!fIsClonable || fFruits == 0) return;
   fClones = fFruits->Clone();
}

//_____________________________________________________________________________
void AliFMaker::Init()
{
   //dummy
}

//_____________________________________________________________________________
void AliFMaker::Finish()
{

   //dummy
}

//_____________________________________________________________________________
void AliFMaker::Make()
{

   Warning("Make","Dummy function called");
}

//_____________________________________________________________________________
void AliFMaker::PrintInfo()
{
   printf("*************************************************************\n");
   printf("*                                                           *\n");
   printf("*            %25s                      *\n",GetName());
   printf("*                                                           *\n");
   printf("*************************************************************\n");

   Dump();
}

//_____________________________________________________________________________
void AliFMaker::MakeBranch()
{
//   Adds the list of physics objects to the AliFast tree as a new branch

   if (fSave == 0) return;

   TTree *tree = gAliFast->Tree();
   if (tree == 0  || fFruits == 0  || fBranchName.Length() == 0) return;

//  Make a branch tree if a branch name has been set
   Int_t buffersize = 4000;
   if (fFruits->InheritsFrom("TClonesArray")) {
      tree->Branch(fBranchName.Data(), &fFruits, buffersize);
   } else {
      tree->Branch(fBranchName.Data(),fFruits->ClassName(), &fFruits, buffersize);
   }
}

//_____________________________________________________________________________
void AliFMaker::SetChainAddress(TChain *chain)
{
//   Set branch address in a chain of files

   if (chain == 0) return;

   chain->SetBranchAddress(fBranchName.Data(), &fFruits);
}

//______________________________________________________________________________
void AliFMaker::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliFMaker.

   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      
      AliFMaker::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
          //this is an addition to the standard rootcint version of Streamer
          //branch address for this maker is set automatically
      TTree *tree = gAliFast->Tree();
      if (tree == 0  || fFruits == 0  || fBranchName.Length() == 0) return;
      TBranch *branch = tree->GetBranch(fBranchName.Data());
      if (branch)  branch->SetAddress(&fFruits);
   } else {
      AliFMaker::Class()->WriteBuffer(R__b,this);
   }
}














