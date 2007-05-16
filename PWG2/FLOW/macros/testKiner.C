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

const char* aDataDir = "./" ; Int_t aRuns = -1 ; Int_t offset = 0 ;

//////////////////////////////////////////////////////////////////////////////////////////////////////

int testKiner(int cen = -1)
{
 cout << " . Here the new flow kinemaker (2007 noRL) ... " << endl ;
 cout << endl ;

 int limit = -1 ;
 if(limit>0)
 {
  cout << " . limited to " << limit << "events . " << endl ;
  cout << endl ;
 }
 
 // flowEvents file (output) //

 TString output = "flowKevts" ;
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

 // flow maker //

 AliFlowKineMaker * flowKiner = new  AliFlowKineMaker() ;
 // cuts, etc.
 flowKiner->SetAbsEtaCut(10.) ;
 flowKiner->SetECut(0.001,100.) ;
 //flowKiner->SetLabelCut(..,..) ;
 flowKiner->SetPrimaryCut(kTRUE) ;
 flowKiner->PrintCutList() ;

 AliFlowEvent * flowEvt = 0 ;
 Int_t count = 0 ;

 // loop (folders) //

 TString execDir(gSystem->pwd());
 TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
 TList* dirList 	   = baseDir->GetListOfFiles();
 Int_t nDirs		   = dirList->GetEntries();
 gSystem->cd(execDir);
 
 for(Int_t iDir=0; iDir<nDirs; ++iDir)
 {
  TSystemFile* presentDir = (TSystemFile*)dirList->At(iDir) ;
  if(!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0) 
  {
   cout << endl ; 
   cout << "Directory (" << iDir << "):  " << presentDir->GetName() << " - Skipping ... " << endl ;
   continue ;   
  }
  if(offset > 0)  { --offset ; continue ; }
  if((aRuns > 0) && (count >= aRuns)) { break ; }
 
  TString presentDirName(aDataDir);
  presentDirName += presentDir->GetName();
  presentDirName += "/";

  TString fileName = presentDirName ; 
  fileName += "galice.root" ;
  Long_t *id, *size, *flags, *modtime ;
  if(gSystem->GetPathInfo(fileName.Data(),id,size,flags,modtime)) 
  { 
   cout << " File : " << fileName << " does NOT exist ! - Skipping ... " << endl ; 
   continue ; 
  }
  cout << endl ; cout << "Directory (" << iDir << "):  " << presentDirName << "  ... " << endl ;

 // loop (simulations in the present dir) //

  TSystemDirectory* evtsDir = new TSystemDirectory(".", presentDirName.Data());
  TList* fileList 	    = evtsDir->GetListOfFiles();
  Int_t nFiles		    = fileList->GetEntries();
  gSystem->cd(execDir);

  for(Int_t iFiles=0; iFiles<nFiles; ++iFiles)
  {
   TSystemFile* presentFile = (TSystemFile*) fileList->At(iFiles);

   TString presentFileName(presentDirName);
   presentFileName += presentFile->GetName();

   if(!(presentFileName.Contains("Kinematics") && presentFileName.Contains("root"))) { continue ; }

   cout << " found: " << presentFileName.Data() << endl ; 
  
   TFile* kineFile = new TFile(presentFileName.Data(), "READ") ; 
   // kineFile->ls() ;
   Int_t nEvts = kineFile->GetNkeys() ; 
   cout << "  . found: " << nEvts << " KineTree(s) in " << presentFileName.Data() << endl ;
   TList* kineEventsList = (TList*)kineFile->GetListOfKeys() ; 
   TTree* kTree ;
   TIter next(kineEventsList); 
   TKey* key ;

   // Loop over the events
   while( key=(TKey *)next() ) 
   {
    TDirectory* tDir = (TDirectory*)key->ReadObj() ;
    if(!tDir) break;
 
    TString evtDir(tDir->GetName()) ; 
    cout << "  . . found: " << tDir->GetName() << endl ;

    kTree = (TTree *)tDir->Get("TreeK");
    if(!kTree) break;

    Int_t nPart = kTree->GetEntries() ;
    cout << "  . . . kTree " << count << " has " << nPart << " particles " << endl ;
    
   // fill and save the flow event
    flowEvt = flowKiner->FillFlowEvent(kTree) ;
    cout << "  . . . flowEvent " << flowEvt << " filled from ttree " << kTree << " ... " << endl ; 
    // flowEvt->Dump() ; cout << endl ;
    TString evtID = "" ; evtID += iDir ; evtID += "-" ; evtID += evtDir ; 
    fFlowfile->cd() ; flowEvt->Write(evtID.Data()) ;
    cout <<  "  . . . flowEvent " << flowEvt << "  ( " << evtID.Data() << " )  -  written on disk (" << output << ") ... " << endl;
    delete flowEvt ; cout << endl ;
   // -

    if(count == limit) { break ; }
    count ++ ;
    
    delete kTree ;
   }
   delete kineFile ;
  }
  delete evtsDir ;
 }
 
 fFlowfile->Close() ; 
 
 cout <<  endl ;
 cout << " Finished ... " << endl ;
 cout << "  nParticles:  " << (flowKiner->GetNgoodTracks() + flowKiner->GetNgoodV0s()) << " (" << flowKiner->GetNposiTracks() << "+ , " << flowKiner->GetNnegaTracks() << "- , " << flowKiner->GetNgoodV0s() << " neutral) " << endl ;   
 cout << "  <nCharged> (|eta|<0.5):  " << (Int_t)(flowKiner->GetNgoodTracksEta()/count) << endl ; 
 cout << "  Bayesian :  " ; 
 for(int ii=0;ii<5;ii++) { cout << flowKiner->GetBayesianNorm(ii) << "   " ; } 
 cout << " . " << endl ; 

 timer.Stop() ;
 cout << endl ;
 timer.Print() ;
 cout << " . here it was (kiner noRL) ... " << endl ;  //juice!
 cout << endl ;

 // cout << endl ; cout << " Memory Check (from Paul)" << endl ; 
 // gObjectTable->Print();
 // cout << endl ; cout << endl ;

 return cen ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

