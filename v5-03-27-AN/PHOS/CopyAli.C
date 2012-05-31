////////////////////////////////////////////////////////////////////////
//
// name: CopyAli.C
// date: 17.3.2002
// last update: 21.3.2002
// author: Jiri Chudoba
// version: 1.0
//
// description: 
//    copy some alice objects from 1 file to another
//
// ToDo: 
//    add support for more events in 1 file
//
//
// Note:
//    copied objects are not deleted, I assume that the root
//    session is exited after this copy
//
// Example:
//     aliroot -b -q CopyAli.C\("galice.root","TreeK.root",0,1,1,1\)
//
// History:
//
// 21.3.02 - first version
//
////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "iostream.h"
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "AliRun.h"
#include "TParticle.h"
#include "AliHeader.h"
#include "TObjArray.h"
#include "TString.h"
#endif

void CopyAli(TString fnOrig="rfio:galice.root", TString fnNew="galice_new.root",Int_t iEvent = 0, Bool_t copygAlice=kTRUE,Bool_t copyTreeK = kFALSE) 
{

  TFile *fileOrig = TFile::Open(fnOrig);
  if (!fileOrig->IsOpen()) {
    cerr<<"Cannot open input file "<<fnOrig.Data()<<endl;
    return;
  }
//  AliRun *gAlice;
  if (gAlice) {delete gAlice; gAlice = 0;}  
  gAlice = (AliRun*)(fileOrig->Get("gAlice"));
  if (!gAlice) {
    cerr<<"Cannot read gAlice from the input file"<<endl;
    return;
  }

  Int_t nAllTracks = gAlice->GetEvent(iEvent);
  cout<<"nAllTracks = "<<nAllTracks<<endl;

// Open the new file

  TFile *fileNew = TFile::Open(fnNew,"update");
  if (!fileNew->IsOpen()) {
    cerr<<"Cannot open output file "<<fnNew.Data()<<endl;
    return;
  }
  if (copygAlice) {
    cout<<"Copy gAlice: ";
    gAlice->Write();
    cout<<"done"<<endl;

    TTree *treeE  = gAlice->TreeE();
    if (!treeE) {
      cerr<<"No TreeE found for event "<<iEvent<<endl;
      return;
    }      
    cout<<"Copy TreeE: ";
    AliHeader *header = new AliHeader();
    treeE->SetBranchAddress("Header", &header);
    treeE->SetBranchStatus("*",1);
    TTree *treeENew =  treeE->CloneTree();
    treeENew->Write();
    cout<<"done"<<endl;

  if (copyTreeK) {
    cout<<"Copy TreeK: ";
    TTree *treeK  = gAlice->TreeK();
    if (!treeK) {
      cerr<<"No TreeK found for event "<<iEvent<<endl;
      return;
    }
    TParticle *particle = new TParticle();
    treeK->SetBranchAddress("Particles",&particle);
    treeK->SetBranchStatus("*",1);
    TTree *treeKNew =  treeK->CloneTree();
    treeKNew->Write();
    cout<<"done"<<endl;
  }

 
  delete gAlice;
  gAlice = 0;
  fileNew->Close();
  fileOrig->Close();
  delete fileNew;
  delete fileOrig;
}
