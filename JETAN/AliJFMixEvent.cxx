// $Id$

#include <Riostream.h>

#include <TParticle.h>
#include <TTree.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TChain.h>

#include "AliJFMixEvent.h"

AliJFMixEvent::AliJFMixEvent(Char_t *file1, Char_t *file2) : 
  fStatus(-1),fEvent1(-1),fEvent2(-1),fMaxEvent1(-1),fMaxEvent2(-1),
  fMarkPythia(0),fTree1(0),fTree2(0),fFile1(0),
  fFile2(0),fBranch1(0),fBranch2(0),
  fParticles1(0),fParticles2(0),fMixedParticles(0)
{
  Init(file1,file2);
}

AliJFMixEvent::AliJFMixEvent(Char_t *files1, Char_t *tname1, Char_t *files2=0, Char_t *tname2=0) :
  fStatus(-1),fEvent1(-1),fEvent2(-1),fMaxEvent1(-1),fMaxEvent2(-1),
  fMarkPythia(0),fTree1(0),fTree2(0),fFile1(0),
  fFile2(0),fBranch1(0),fBranch2(0),
  fParticles1(0),fParticles2(0),fMixedParticles(0)
{
  InitChains(files1,tname1,files2,tname2);
}

AliJFMixEvent::~AliJFMixEvent() 
{
  Clean();
}

void AliJFMixEvent::Init(Char_t *file1, Char_t *file2)
{
  Clean();

  if (!file1) {
    cerr << "Error AliJFMixEvent:: Must give at least one filename!" << endl;
    return;
  }
  
  fFile1 = NULL;
  fFile1 = new TFile(file1,"READ");
  if (!fFile1 || !(fFile1->IsOpen())) {
    cerr << "Error AliJFMixEvent:: Couldn't open input file " << file1 << endl;
    return;
  }

  fFile2 = NULL;
  if(file2!=0){
    fFile2 = new TFile(file2,"READ");
    if (!fFile2 || !(fFile2->IsOpen())) {
      cerr << "Error AliJFMixEvent:: Couldn't open input file " << file2 << endl;
      return;
    }
  }

  Int_t prealloc1=100000;
  fTree1 = NULL;
  fTree1 = (TTree *)fFile1->Get("hijing");
  if (!fTree1){
    fTree1 = (TTree *)fFile1->Get("pythia");
    prealloc1=10000;
  }
  if (!fTree1) {
    cerr << "Error AliJFMixEvent:: Didn't find a TParticle tree for " << file1 << endl;
    return;
  }
  fMaxEvent1=(Int_t)fTree1->GetEntries();

  Int_t prealloc2=0;
  fTree2 = NULL;
  if(file2!=0){
    prealloc2=100000;
    fTree2 = (TTree *)fFile2->Get("hijing");
    if (!fTree2){
      fTree2 = (TTree *)fFile2->Get("pythia");
      prealloc2=10000;
    }
    if (!fTree2) {
      cerr << "Error AliJFMixEvent:: Didn't find a TParticle tree for " << file2 << endl;
      return;
    }
    fMaxEvent2=(Int_t)fTree2->GetEntries();
  }

  fParticles1=new TClonesArray("TParticle",prealloc1);
  if(file2!=0) fParticles2=new TClonesArray("TParticle",prealloc2);
  fMixedParticles=new TClonesArray("TParticle",prealloc1+prealloc2);

  fBranch1=fTree1->GetBranch("particles");
  fBranch1->SetAddress(&fParticles1);

  if(file2!=0){
    fBranch2=fTree2->GetBranch("particles");
    fBranch2->SetAddress(&fParticles2);
  }

  fStatus=0;
}

void AliJFMixEvent::InitChains(Char_t *files1, Char_t *tname1, Char_t *files2, Char_t *tname2)
{
  Clean();

  if (!files1 || !tname1) {
    cerr << "Error AliJFMixEvent:: Must give at least one file/tree name!" << endl;
    return;
  }
  
  Int_t prealloc1=100000;
  Int_t prealloc2=0;

  fTree1=new TChain(tname1);
  ((TChain*)fTree1)->Add(files1);
  fParticles1=new TClonesArray("TParticle",prealloc1);
  fTree1->SetBranchAddress("particles",&fParticles1);
  fMaxEvent1=(Int_t)fTree1->GetEntries();

  if (files2 && tname2) {
    prealloc2=100000;
    fTree2=new TChain(tname2);
    ((TChain*)fTree2)->Add(files2);
    fParticles2=new TClonesArray("TParticle",prealloc2);
    fTree2->SetBranchAddress("particles",&fParticles2);
    fMaxEvent2=(Int_t)fTree2->GetEntries();
  }

  fMixedParticles=new TClonesArray("TParticle",prealloc1+prealloc2);
  fStatus=0;
}

void AliJFMixEvent::Clean()
{
  fStatus=-1; //not initialized!

  fEvent1=-1;
  fEvent2=-1;
  fMaxEvent1=-1;
  fMaxEvent2=-1;

  if(fParticles1) delete fParticles1;
  fParticles1=NULL;

  if(fParticles2) delete fParticles2;
  fParticles2=NULL;

  if(fMixedParticles) delete fMixedParticles;
  fMixedParticles=NULL;

  if(fTree1) delete fTree1;
  fTree1=NULL;

  if(fTree2) delete fTree2;
  fTree2=NULL;

  if(fFile1) delete fFile1;
  fFile1=NULL;

  if(fFile2) delete fFile2;
  fFile2=NULL;
}

Int_t AliJFMixEvent::CreateMixedEvent(Int_t i,Int_t j)
{
  if(fStatus<0) return -1;
  if(i<0 || j<0) return -1;
  if((i>=fMaxEvent1)||((fTree2)&&(j>=fMaxEvent2))) return -1;

  if(fEvent1!=i){
    fTree1->GetEvent(i);
    fEvent1=i;
  }

  if(fTree2){
    if(fEvent2!=j){
      fTree2->GetEvent(j);
      fEvent2=j;
    }
  }
  return MixEvent();
}

Int_t AliJFMixEvent::CreateNextMixedEvent()
{
  if(fStatus<0) return -1;

  if(fStatus==0) { //load first data samples
    fEvent1=0;
    fTree1->GetEvent(0);
    fEvent2=0;
    if(fTree2) fTree2->GetEvent(0);
  } else { //continue loading
    fEvent1++;
    if(fEvent1==fMaxEvent1){
      fEvent1=0;

      if(fTree2){
	fEvent2++;
	if(fEvent2==fMaxEvent2) 
	  fEvent2=0;
	fTree2->GetEvent(fEvent2);
      }
    }

    fTree1->GetEvent(fEvent1);
  }

  return MixEvent();
}

Int_t AliJFMixEvent::MixEvent()
{
  Int_t n1=fParticles1->GetEntries();
  Int_t n2=0;
  if(fTree2) n2=fParticles2->GetEntries();

  fMixedParticles->ExpandCreateFast(n1+n2);
  Int_t n=0;
  TParticle *particle;
  TIterator *i=fParticles1->MakeIterator();
  while ((particle = (TParticle *) i->Next()) != NULL) {
    if(fMarkPythia){
      particle->SetWeight(-123); //mark pythia particles
    }
    new ((*fMixedParticles)[n]) TParticle(*particle);
    n++;
  }
  if(fTree2){
    i=fParticles2->MakeIterator();
    while ((particle = (TParticle *) i->Next()) != NULL) {
      new((*fMixedParticles)[n]) TParticle(*particle);
      n++;
    }
  }
  fStatus++; 

  return n;
}

void AliJFMixEvent::Debug()
{
  if(!fMixedParticles) return;

  TParticle *p;
  TIterator *iter=fMixedParticles->MakeIterator();
  while ((p = (TParticle *) iter->Next()) != NULL) {
    p->Print();
  }
}
