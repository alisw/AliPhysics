/////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
/////////////////////////////////////////////////////////////
//
// Description: AliRoot macro to make AliFlowEvents from KineTree (new way) 
//
/////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TFile.h"
#include "TObjArray"
#include "TStopwatch.h"

using namespace std; //required for resolving the 'cout' symbol

TTree* kOpen(int evtN=0) ;
//TChain* CreateKineChain(const char* aDataDir = "MCfiles.txt", Int_t aRuns = 10, Int_t offset = 0) ;
//void LookupWrite(TChain* chain, const char* target) ;

void testKiner(TString output = "flowKevts.root")
{
 cout << " . Here the new flow kinemaker (2007b) ... " << endl ;
 cout << endl ;

 bool kOne = kFALSE ;

 TStopwatch timer;
 timer.Start();

 gSystem->Load("AliFlow_All.so");

 // output file //

 TFile * fFlowfile = new TFile(output.Data(),"RECREATE") ;
 //fFlowfile->cd() ; 

 // esd chain //

// // TString fESDfileName = "AliESDs.root" ; TString fESDtree = "esdTree" ; 
// TChain* pMCchain = CreateESDChain(".",10,0);
// Int_t fNumberOfEvents = (Int_t)pMCchain->GetEntries() ;
// cout << " tot. " << fNumberOfEvents << " events in the TChain ... " << endl ; cout << endl ;

// TString fKineBranch = "TreeK" ;
// TTree treeK = 0 ;
// pMCchain->SetBranchAddress(fKineBranch.Data(),&treeK) ;

 Int_t fNumberOfEvents = 1 ;
 TTree* treeK = 0 ;

 // flow maker //

 AliFlowKineMaker * flowKiner = new  AliFlowKineMaker() ;
 // cuts, etc.
 flowKiner->SetAbsEtaCut(2.1) ;
 flowKiner->SetECut(0.01,100.) ;
 //flowKiner->SetLabelCut(..,..) ;
 flowKiner->SetPrimaryCut(Bool_t prim = kTRUE) ;
 flowKiner->PrintCutList() ;

 // loop //

 Int_t evtN = 0 ;
 AliFlowEvent * flowEvt = 0 ;
 for(evtN=0;evtN<fNumberOfEvents;evtN++)
 {
  treeK = kOpen(0) ;
  
  Int_t evtNN = -1 ;
  Int_t nTrk = treeK->GetEntries() ;

  cout << endl ; cout << " Event " << evtN << "  ( " << evtNN << " )  : " << nTrk << " tracks  ." << endl ;

  flowEvt = flowKiner->FillFlowEvent(treeK) ;

  cout << " Event filled " << flowEvt << " ... " << endl ;
  // cout << endl ; cout << " trks : " << flowEvt->TrackCollection()->GetEntries() << endl ;
  // flowEvt->Dump() ; cout << endl ;

  TString evtID = "" ; evtID += evtN ; 
  fFlowfile->cd() ; 
  flowEvt->Write(evtID.Data()) ;
  cout <<  " Event " << evtN << "  ( " << evtID.Data() << " )  -  written on disk (" << output << ") ." << endl;
  delete flowEvt ;
 }
 
 fFlowfile->Close() ; 
 
 cout <<  endl ;
 cout << " Finished ... " << endl ;
 cout << "  nTracks:  " << flowKiner->GetNgoodTracks() << endl ;   
 cout << "  nV0s:  " << flowKiner->GetNgoodV0s()  << endl ;  	     
 cout << "  nTracks (|eta|<0.5):  " << flowKiner->GetNgoodTracksEta() << endl ; 
 cout << "  nTracks+:  " << flowKiner->GetNposiTracks() << endl ; 	     
 cout << "  nTracks-:  " << flowKiner->GetNnegaTracks() << endl ; 	     
 cout << "  nTracks unconstrained:  " << flowKiner->GetNunconstrained() << endl ; 	 
 //cout << "  Bayesian :  " ; 
 //for(int ii=0;ii<5;ii++) { cout << flowKiner->GetBayesianNorm(ii) << "   " ; } 
 cout << " . " << endl ; 

 timer.Stop() ;
 cout << endl ;
 timer.Print() ;
 cout << " . here it was (kiner) ... " << endl ;  //juice!
 cout << endl ;

 // break ;

}

TTree* kOpen(int evtN)
{
 TString fileName = "./0/galice.root" ;
 AliRunLoader* rl = AliRunLoader::Open(fileName.Data(),"MyEvent","read");
 rl->LoadgAlice();
 AliRun* gAlice = rl->GetAliRun();
 rl->LoadHeader();
 rl->LoadKinematics();
 Int_t fNumberOfEvents = rl->GetNumberOfEvents() ;
 cout << " Found :  " << fNumberOfEvents << "  event(s) ... " << endl ; 

 Int_t exitStatus = rl->GetEvent(evtN) ; if(exitStatus!=0) { return 0 ; }

 TTree* kTree = (TTree*)rl->TreeK();	      // Particles TTree (KineTree)
 AliStack* kStack = gAlice->Stack();	      // Particles Stack (use "Label()" to get the number in the stack)

 Int_t fNumberOfParticles = kTree->GetEntries() ;
 Int_t nPart = kStack->GetNtrack() ;
 cout << " Event n. " << evtN << "  contains :  " << fNumberOfParticles << "  particles in the TTree  ( =  " << nPart << "  in the stack ) . " << endl ; 

 return kTree ; 
}


// // Helper macros for creating chains (from: CreateESDChain.C,v 1.10 jgrosseo Exp)
// TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset)
// {
//   // creates chain of files in a given directory or file containing a list.
//   // In case of directory the structure is expected as:
//   // <aDataDir>/<dir0>/AliESDs.root
//   // <aDataDir>/<dir1>/AliESDs.root
//   // ...
// 
//   if (!aDataDir)
//     return 0;
// 
//   Long_t id, size, flags, modtime;
//   if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
//   {
//     printf("%s not found.\n", aDataDir);
//     return 0;
//   }
// 
//   TChain* chain = new TChain("esdTree");
//   TChain* chaingAlice = 0;
// 
//   if (flags & 2)
//   {
//     TString execDir(gSystem->pwd());
//     TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
//     TList* dirList            = baseDir->GetListOfFiles();
//     Int_t nDirs               = dirList->GetEntries();
//     gSystem->cd(execDir);
// 
//     Int_t count = 0;
// 
//     for (Int_t iDir=0; iDir<nDirs; ++iDir)
//     {
//       TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
//       if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
//         continue;
// 
//       if (offset > 0)
//       {
//         --offset;
//         continue;
//       }
// 
//       if (count++ == aRuns)
//         break;
// 
//       TString presentDirName(aDataDir);
//       presentDirName += "/";
//       presentDirName += presentDir->GetName();
// 
//       chain->Add(presentDirName + "/AliESDs.root/esdTree");
//     }
//   }
//   else
//   {
//     // Open the input stream
//     ifstream in;
//     in.open(aDataDir);
// 
//     Int_t count = 0;
// 
//     // Read the input list of files and add them to the chain
//     TString esdfile;
//     while(in.good()) {
//       in >> esdfile;
//       if (!esdfile.Contains("root")) continue; // protection
// 
//       if (offset > 0)
//       {
//         --offset;
//         continue;
//       }
// 
//       if (count++ == aRuns)
//         break;
// 
//         // add esd file
//       chain->Add(esdfile);
//     }
// 
//     in.close();
//   }
// 
//   return chain;
// }
// 
// void LookupWrite(TChain* chain, const char* target)
// {
//   // looks up the chain and writes the remaining files to the text file target
// 
//   chain->Lookup();
// 
//   TObjArray* list = chain->GetListOfFiles();
//   TIterator* iter = list->MakeIterator();
//   TObject* obj = 0;
// 
//   ofstream outfile;
//   outfile.open(target);
// 
//   while ((obj = iter->Next()))
//     outfile << obj->GetTitle() << "#AliESDs.root" << endl;
// 
//   outfile.close();
// 
//   delete iter;
// }
