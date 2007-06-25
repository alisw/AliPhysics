/*

 /

  STANDALONE - no PROOF
  .x ~/rootlogon.C
  .L $ALICE_ROOT/TPC/TPCcalib/AliTPCcalibTracks.cxx+
  .L $ALICE_ROOT/TPC/TPCcalib/AliTPCSelectorTracks.cxx+
  .L $ALICE_ROOT/TPC/TPCcalib/tpcSelectorTracks.C+
  gMaxFiles = 50;
  //  TChain * chain = makeChain("list.list", kTRUE,kFALSE, kTRUE);
  TChain * chain = makeChain("listTr.list", kTRUE,kTRUE, kFALSE);
  chain->SetBranchStatus("*",1);
  chain->Process("$ALICE_ROOT/TPC/TPCcalib/AliTPCSelectorTracks.cxx+"); 

 
  
*/


#include <fstream.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TDSet.h>
#include <TVectorF.h>
#include <TLinearFitter.h>
#include "AliESD.h"




char  *prefix = "root://lxfs35.gsi.de:1094//alice/testtpc/rec0606";


Int_t gMaxFiles  =10000;
TChain * gChain = 0;
void MakeSet( char * ifile, TChain *cESDTree=0, TChain *cESDFriend=0, Bool_t check=kTRUE);
void MakeSetZip( char * ifile, TChain *cESDTree=0, TChain *cESDFriend=0, Bool_t check=kTRUE);



TChain * makeChain(char * input, Bool_t check, Bool_t bFriend=kFALSE, Bool_t bZip = kTRUE ){
  //
  //
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");
  
  TChain *cESDTree = new TChain("esdTree");
  TChain *cESDFriend = 0;
  if (bFriend) cESDFriend = new TChain("esdFriendTree");
  if (bFriend) {
    MakeSet(input,cESDTree,cESDFriend, check);
  }else{
    if (bZip){
      MakeSetZip(input,cESDTree,cESDFriend, check);
    }
    else{
      MakeSet(input,cESDTree,cESDFriend, check);
    }
  }
  if (bFriend) cESDTree->AddFriend(cESDFriend,"kokot");  
  cESDTree->Lookup();
  gChain=cESDTree;
  return cESDTree;
}


void MakeSet( char * ifile, TChain *cESDTree, TChain *cESDFriend, Bool_t check){
  //
  //
  // 
  char dir[100];
  char esdFile[100];
  char friendFile[100];
  ifstream in(ifile);
  Int_t count =0;
  while (!in.eof()) {
   if (in.eof()) break;
   in >> dir;
   if (in.eof()) break;
   sprintf(esdFile,"%s/%s/AliESDs.root",prefix,gSystem->DirName(Form("%s/",dir)));
   sprintf(friendFile,"%s/%s/AliESDfriends.root",prefix,gSystem->DirName(Form("%s/",dir)));
   //
   // check file
   if (check){
     Int_t entries0=0;
     Int_t entries1=0;
     TFile * fESD = TFile::Open(esdFile);
     if (!fESD) continue;
     TTree * treeESD = (TTree*)fESD->Get("esdTree");
     if (!treeESD) {delete fESD; continue;}
     else{
       entries0 = treeESD->GetEntries();
     }
     delete fESD;     
     TFile * fFriend = TFile::Open(friendFile);
     if (!fFriend) continue;
     //
     if (cESDFriend){
       TTree * tFriend = (TTree*)fFriend->Get("esdFriendTree");
       if (!tFriend) { delete tFriend; continue;}
       
       else{
	 entries1 = tFriend->GetEntries();
       }
       delete fFriend;
       if (entries0!=entries1) continue;
       if (entries0==0) continue;
     }
   }
   //
   //  
   if ( cESDTree)   cESDTree->Add(esdFile);
   if (cESDFriend)   cESDFriend->Add(friendFile);
   printf ("%s\n",dir);
   count++;
   if (count>gMaxFiles) break;
 }
}


void MakeSetZip( char * ifile, TChain *cESDTree, TChain */*cESDFriend*/, Bool_t check){
  //
  //
  // 
  char dir[500];
  char esdFile[500];
  char friendFile[500];
  ifstream in(ifile);
  Int_t count=0;
  while (!in.eof()) {
   if (in.eof()) break;
   in >> dir;
   if (in.eof()) break;
   sprintf(esdFile,"%s/root_archive.zip#AliESDs.root",dir);
   sprintf(friendFile,"%s/root_archive.zip#AliESDfriends.root",dir);
   printf("%s\n", esdFile);
   //
   // check file
   if (check){
     Int_t entries0=0;
     //     Int_t entries1=0;
     TFile * fESD = TFile::Open(esdFile);
     if (!fESD) continue;
     TTree * treeESD = (TTree*)fESD->Get("esdTree");
     if (!treeESD) {delete fESD; continue;}
     else{
       entries0 = treeESD->GetEntries();
     }
     delete fESD;     
     //TFile * fFriend = TFile::Open(friendFile);
     //if (!fFriend) continue;
   }
   //
   //  
   if ( cESDTree)   cESDTree->Add(esdFile);
   count++;
   if (count>gMaxFiles) break;
 }
}



