/////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
/////////////////////////////////////////////////////////////
//
// Description: AliRoot macro to make AliFlowEvents from AliESDs (new way) 
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

TChain* CreateESDChain(const char* aDataDir = "ESDfiles.txt", Int_t aRuns = -1, Int_t offset = 0) ;
void LookupWrite(TChain* chain, const char* target) ;

//////////////////////////////////////////////////////////////////////////////////////////////////////

int testMaker(int cen = -1)
{
 cout << " . Here the new flow maker (2007 rs) ... " << endl ;
 cout << endl ;

 bool kOne = kFALSE ;
 int limit = -1 ;

 if(limit>0)
 {
  cout << " . limited to " << limit << "events . " << endl ;
  cout << endl ;
 }
 
 // flowEvents file (output) //

 TString output = "flowEvts" ;
 if(cen >= 0) { output += cen ; }
 output += ".root" ;
 
 // start, load libs //

 TStopwatch timer;
 timer.Start();

 //gSystem->Load("libPhysics.so");
 gSystem->Load("libPWG2flow.so");

 // open output file //

 TFile * fFlowfile = new TFile(output.Data(),"RECREATE") ;
 //fFlowfile->cd() ; 

 // esd chain (imput) //

 TString baseDir = "." ; 
 TChain* pESDchain = CreateESDChain(baseDir.Data(),limit,0);
 Int_t fNumberOfEvents = (Int_t)pESDchain->GetEntries() ;
 cout << " tot. " << fNumberOfEvents << " events in the TChain ... " << endl ; cout << endl ;

 TString fESDbranch = "ESD" ;
 AliESD * pEsd = 0 ;
 pESDchain->SetBranchAddress(fESDbranch.Data(),&pEsd) ;

 // flow maker //

 AliFlowMaker * flowMaker = new  AliFlowMaker() ;

 // flags 
 flowMaker->DoTracks(kTRUE) ;
 flowMaker->DoV0s(kFALSE) ;     // kTRUE

 // cuts, etc.
 flowMaker->SetNHitsCut(1) ;
 flowMaker->SetECut(0.01,100.) ;
 //flowMaker->SetLabelCut(..,..) ;
 flowMaker->PrintCutList() ;

 // loop //

 Int_t evtN = 0 ;
 Int_t count = 0 ;
 AliFlowEvent * flowEvt = 0 ;
 for(evtN=0;evtN<fNumberOfEvents;evtN++)
 {
  pESDchain->GetEntry(evtN,1) ;

  Int_t evtNN = -1 ;
  // Int_t evtNN = pEsd->GetEventNumber() ;
  Int_t nTrk = pEsd->GetNumberOfTracks() ;
  Int_t nV0s = pEsd->GetNumberOfV0s() ;
  cout << endl ; cout << " Event " << evtN << "  ( " << evtNN << " )  : " << nTrk << " tracks  &  " << nV0s << " v0s ." << endl ;

  flowEvt = flowMaker->FillFlowEvent(pEsd) ;
  // cout << " Event filled " << flowEvt << " ... " << endl ; // flowEvt->Dump() ; cout << endl ;

  TString evtID = "" ; evtID += evtN ; 
  fFlowfile->cd() ; flowEvt->Write(evtID.Data()) ;
  cout <<  " Event " << evtN << "  ( " << evtID.Data() << " )  -  written on disk (" << output << ") ." << endl;
  delete flowEvt ;

  if(count == limit) { break ; }  
  count++ ;  
 }
 
 fFlowfile->Close() ; 
 
 cout <<  endl ;
 cout << " Finished ... " << endl ;
 cout << "  nTracks:  " << flowMaker->GetNgoodTracks() << " (" << flowMaker->GetNposiTracks() << "+ , " << flowMaker->GetNnegaTracks() << "-) " << endl ;   
 cout << "  nV0s:  " << flowMaker->GetNgoodV0s() << endl ;  	     
 cout << "  nTracks (|eta|<0.5):  " << flowMaker->GetNgoodTracksEta() << endl ; 
 // cout << "  nTracks unconstrained:  " << flowMaker->GetNunconstrained() << endl ; 	 
 cout << "  Bayesian :  " ; 
 for(int ii=0;ii<5;ii++) { cout << flowMaker->GetBayesianNorm(ii) << "   " ; } 
 cout << " . " << endl ; 

 timer.Stop() ;
 cout << endl ;
 timer.Print() ;
 cout << " . here it was (maker rs) ... " << endl ;  //juice!
 cout << endl ;
 
 // cout << endl ; cout << " Memory Check (from Paul)" << endl ; 
 // gObjectTable->Print();
 // cout << endl ; cout << endl ;

 return cen ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

// Helper macros for creating chains (from: CreateESDChain.C,v 1.10 jgrosseo Exp)
TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliESDs.root
  // <aDataDir>/<dir1>/AliESDs.root
  // ...

  if (!aDataDir)
    return 0;

  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
  {
    printf("%s not found.\n", aDataDir);
    return 0;
  }

  TChain* chain = new TChain("esdTree");
  TChain* chaingAlice = 0;

  if (flags & 2)
  {
    TString execDir(gSystem->pwd());
    TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
    TList* dirList            = baseDir->GetListOfFiles();
    Int_t nDirs               = dirList->GetEntries();
    gSystem->cd(execDir);

    cout << "* Reading :  " << aDataDir << "   (+" << nDirs << ")" << endl ;
      
    Int_t count = 0;

    for (Int_t iDir=0; iDir<nDirs; ++iDir)
    {
      TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
      if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
      {
        cout << " Skip : " << presentDir->GetName() << endl ;
      
        continue;
      }

      if (offset > 0)
      {
        cout << " OffSet : " << offset << endl ;
      
        --offset;
        continue;
      }

      if (count++ == aRuns)
      {
        cout << " aRuns : " << aRuns << endl ;
	
        break;
      }

      TString presentDirName(aDataDir);
      presentDirName += "/";
      presentDirName += presentDir->GetName();

      chain->Add(presentDirName + "/AliESDs.root/esdTree");
    }
  }
  else
  {
    // Open the input stream
    ifstream in;
    in.open(aDataDir);

    Int_t count = 0;

    // Read the input list of files and add them to the chain
    TString esdfile;
    while(in.good()) {
      in >> esdfile;
      if (!esdfile.Contains("root")) continue; // protection

      if (offset > 0)
      {
        --offset;
        continue;
      }

      if (count++ == aRuns)
        break;

        // add esd file
      chain->Add(esdfile);
    }

    in.close();
  }

  return chain;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void LookupWrite(TChain* chain, const char* target)
{
  // looks up the chain and writes the remaining files to the text file target

  chain->Lookup();

  TObjArray* list = chain->GetListOfFiles();
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  ofstream outfile;
  outfile.open(target);

  while ((obj = iter->Next()))
    outfile << obj->GetTitle() << "#AliESDs.root" << endl;

  outfile.close();

  delete iter;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
