/* AliSelectorLYZ.cxx, v1.0 30/07/2007 kolk Exp */
/* derived from AliSelectorFoF.cxx, v1.1 01/02/2007 esimili Exp */
/* derived from AliSelector.cxx,v 1.17 2006/08/31 jgrosseo Exp */

// The class definition in esdV0.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("AliSelector.C")
// Root > T->Process("AliSelector.C","some options")
// Root > T->Process("AliSelector.C+")
//

#include "AliSelectorLYZ.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TRegexp.h>
#include <TTime.h>
#include <TFriendElement.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TTimeStamp.h>
#include <TMath.h>
#include <TVector2.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TSelector.h>

#include "AliLog.h"		  
#include "AliESD.h"		  
#include "AliESDtrack.h"  
#include "AliESDv0.h"		   

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

#include "AliFlowEvent.h"
#include "AliFlowTrack.h"
#include "AliFlowV0.h"
#include "AliFlowConstants.h"
#include "AliFlowSelection.h"
#include "AliFlowMaker.h"
#include "AliFlowLYZConstants.h"
#include "AliFlowLeeYangZerosMaker.h"
#include "AliFlowTrack.h"


ClassImp(AliSelectorLYZ)

//-----------------------------------------------------------------------

AliSelectorLYZ::AliSelectorLYZ() :
  TSelector(),
  fTree(0),
  fESD(0),
  fCountFiles(0),
  fFirstRunLYZ(kTRUE), //set boolean for firstrun to initial value
  fUseSumLYZ(kTRUE),   //set boolean for use sum to initial value
  fKineFile(0)
 {
  //
  // Constructor. Initialization of pointers
  //
  
  fFlowEventFileName       = "fof_flowEvts.root" ;    
  
}

//-----------------------------------------------------------------------

AliSelectorLYZ::~AliSelectorLYZ()
{
  //
  // Destructor
  //

 if (fTree) { fTree->ResetBranchAddresses() ; }

 if (fESD)
 {
   delete fESD;
   fESD = 0;
 }
}

//-----------------------------------------------------------------------

void AliSelectorLYZ::CheckOptions()
{
  // checks the option string for the debug flag
  
  AliLog::SetClassDebugLevel(ClassName(), AliLog::kInfo);

  TString option = GetOption();

  if (option.Contains("moredebug"))
  {
    printf("Enabling verbose debug mode for %s\n", ClassName());
    AliLog::SetClassDebugLevel(ClassName(), AliLog::kDebug+1);
    AliInfo(Form("Called with option %s.", option.Data()));
  }
  else if (option.Contains("debug"))
  {
    printf("Enabling debug mode for %s\n", ClassName());
    AliLog::SetClassDebugLevel(ClassName(), AliLog::kDebug);
    AliInfo(Form("Called with option %s.", option.Data()));
  }
}

//-----------------------------------------------------------------------

void AliSelectorLYZ::Begin(TTree*)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  cerr << " HERE I begin !!! " << endl ; cout << endl ;
  

  CheckOptions();

  AliDebug(AliLog::kDebug, "============BEGIN===========");

 // flags
 fDoNothing     = kFALSE ;
 fOnFlyAnalysis = kTRUE ; 
 fSaveFlowEvents= kFALSE ;

// Maker part :
 fFlowMaker = new AliFlowMaker() ;
 // ESD Cuts
 fFlowMaker->SetNHitsCut(1) ;
 fFlowMaker->SetECut(0.01,100.) ; 
 fFlowMaker->PrintCutList() ;

// Opens flow event (output) file
 if(fSaveFlowEvents)  
 { 
 // flow events file (output)
  fFlowfile = new TFile(fFlowEventFileName.Data(),"RECREATE") ;
  fFlowfile->cd() ; 
  //cout << " . Writing AliFlowEvents on  : " << fFlowEventFileName.Data()   << "  . " << endl ;
 }
 
 // Analysis part :
 if(fOnFlyAnalysis)  
 { 
   //cout << " . Here the flow selection ... " << endl ;
  cout << endl ;

 // AliFlowSelection...
  fFlowSelect = new AliFlowSelection() ;
  // Event Cuts
  fFlowSelect->SetCentralityCut(-1) ;
  fFlowSelect->SetRunIdCut(-1) ;
  // R.P. calculation cuts
  for(int j=0;j<AliFlowConstants::kHars;j++)
    {
      fFlowSelect->SetEtaCut(0., 2., j, 1) ;
      fFlowSelect->SetPtCut(0.1, 10. , j, 1);  
    }
  fFlowSelect->SetConstrainCut(kTRUE) ;
  fFlowSelect->SetDcaGlobalCut(0.,0.1);
  // Correlation analysis cuts (not all of them)
  fFlowSelect->SetEtaPart(-1.1,1.1);
  fFlowSelect->SetPtPart(0.1,10.);   
  fFlowSelect->SetConstrainablePart(kTRUE);
  fFlowSelect->SetDcaGlobalPart(0.,0.1);
  // V0 analysis cuts (not all of them ... they are useless anyway)
  fFlowSelect->SetV0Mass(0.4875,0.5078) ;	 // Mk0 = 0.49765
  fFlowSelect->SetV0SideBands(0.1) ;
  fFlowSelect->SetV0Pt(0.1,10.) ;
  fFlowSelect->SetV0Eta(-2.1,2.1) ;
  // print list :
  //cout << " . Selection for R.P. calculation: " << endl ;
  fFlowSelect->PrintSelectionList() ;
  //cout << " . Selection for correlation analysis: " << endl ;
  fFlowSelect->PrintList() ;
  //cout << " . Selection for V0 analysis: " << endl ;
  fFlowSelect->PrintV0List() ;

  //cout << " . Here the flow analysis ... " << endl ;
  //cout << endl ;
  
  //AliFlowLYZAnalyser...
  fFlowLYZ = new AliFlowLeeYangZerosMaker(fFlowSelect) ;
  //fFlowLYZ -> SetFirstRun(fFirstRunLYZ);   //kFALSE //kTRUE 
  fFlowLYZ -> SetFirstRun(GetFirstRunLYZ()); //set first run true or false
  fFlowLYZ -> SetUseSum(GetUseSumLYZ());     //set use sum true or false
  fFlowLYZ -> SetDebug(kTRUE) ;              //for debugging porposes

  // first run file (read)
  if (!fFirstRunLYZ){
    TString firstRunFileName = "fof_flowLYZAnal_firstrun.root" ;
    fFirstRunFile = new TFile(firstRunFileName.Data(),"READ");
    if(!fFirstRunFile || fFirstRunFile->IsZombie()) { cerr << " ERROR: NO first Run file... " << endl ; }
    else { fFlowLYZ -> SetFirstRunFile(fFirstRunFile); }
  }

  // analysis file (output)
  if (fFirstRunLYZ) {fFlowLYZAnalysisFileName = "fof_flowLYZAnal_firstrun.root" ;}
  else {fFlowLYZAnalysisFileName = "fof_flowLYZAnal_secondrun.root" ;}
  fFlowLYZ->SetHistFileName(fFlowLYZAnalysisFileName.Data()) ;
  cout << " . Writing Analysis Histograms on  : " << fFlowLYZ->GetHistFileName()   << "  . " << endl ;
 

  fFlowLYZ->Init() ;
  cout << "fFlowLYZ->Init() ;"<<endl;
 }
 cout << "                ... init is done " << endl ; 
}

//-----------------------------------------------------------------------

void AliSelectorLYZ::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  CheckOptions();

  AliDebug(AliLog::kDebug, "=======SLAVEBEGIN========");
  AliDebug(AliLog::kDebug, Form("Hostname: %s", gSystem->HostName()));
  AliDebug(AliLog::kDebug, Form("Time: %s", gSystem->Now().AsString()));

  if (tree != 0) { Init(tree) ; }
}

//-----------------------------------------------------------------------

void AliSelectorLYZ::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.

  AliDebug(AliLog::kDebug, "=========Init==========");

  fTree = tree;

  if (fTree == 0)
  {
   AliDebug(AliLog::kError, "ERROR: tree argument is 0.");
   return;
  }

  // Set branch address
  fTree->SetBranchAddress("ESD", &fESD);
  if (fESD != 0) { AliDebug(AliLog::kInfo, "INFO: Found ESD branch in chain.") ; }
}

//-----------------------------------------------------------------------

Bool_t AliSelectorLYZ::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. Typically here the branch pointers
  // will be retrieved. It is normaly not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed.

  AliDebug(AliLog::kDebug, "=========NOTIFY==========");
  AliDebug(AliLog::kDebug, Form("Hostname: %s", gSystem->HostName()));
  AliDebug(AliLog::kDebug, Form("Time: %s", TTimeStamp(time(0)).AsString()));

  ++fCountFiles;
  if (fTree)
  {
    TFile *f = fTree->GetCurrentFile();
    AliDebug(AliLog::kInfo, Form("Processing %d. file %s", fCountFiles, f->GetName()));
  }
  else
  {
    AliDebug(AliLog::kError, "fTree not available");
  }

  DeleteKinematicsFile();

  return kTRUE;
}

//-----------------------------------------------------------------------

Bool_t AliSelectorLYZ::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.

  // WARNING when a selector is used with a TChain, you must use
  //  the pointer to the current TTree to call GetEntry(entry).
  //  The entry is always the local entry number in the current tree.
  //  Assuming that fTree is the pointer to the TChain being processed,
  //  use fTree->GetTree()->GetEntry(entry).

  AliDebug(AliLog::kDebug, Form("=========PROCESS========== Entry %lld", entry));

  if(!fTree)
  {
   AliDebug(AliLog::kError, "ERROR: fTree is 0.") ;
   return kFALSE ;
  }

  fEventNumber = entry ;
  fTree->GetTree()->GetEntry(fEventNumber) ;

  if(fESD) { AliDebug(AliLog::kDebug, Form("ESD: We have %d tracks.", fESD->GetNumberOfTracks())); }
  cout << " event !!! " << entry << endl ;

  fRunID = fESD->GetRunNumber() ;
  //  fEventNumber = fESD->GetEventNumber() ;
  fEventNumber = -1 ;
  fNumberOfTracks = fESD->GetNumberOfTracks() ;
  fNumberOfV0s = fESD->GetNumberOfV0s() ;

  cout << " *evt n. " << fEventNumber << " (run " << fRunID << ") " << endl ;
  cout << "  tracks: " << fNumberOfTracks << " ,   v0s " << fNumberOfV0s << endl ;

 // Dummy Loop
  if(fDoNothing)  { cout << " ... doing nothing ... " << endl ; return kTRUE; } 

 // Instantiate a new AliFlowEvent
  cout << " filling the flow event :| " << endl ;
  fFlowEvent = fFlowMaker->FillFlowEvent(fESD) ;
  if(!fFlowEvent) { cout << "! something bad occurred !" << endl ; return kFALSE ; }
  else 		  { cout << "# event done :) " << entry << "     # ok ! #" << endl ; }  
     

 // Saves the AliFlowEvent
  if(fSaveFlowEvents)
  { 
   cout << " saving flow event :| " << endl ;
   TString strID = "" ; strID += entry ; 
   fFlowfile->cd() ; fFlowEvent->Write(strID.Data()) ;
   cout << "# event saved :) " << strID.Data() << "     # ok ! #" << endl ; cout << endl ;  
  }

   
  // On-Fly Analysis
  if(fOnFlyAnalysis)  
  {
    Bool_t donelyz = fFlowLYZ->Make(fFlowEvent) ;
    if(donelyz) { cout << "# LYZ analysis done :) " << entry << "         # ok ! #       " << endl ; 
    delete fFlowEvent; }

  }
   
  return kTRUE;
}

//-----------------------------------------------------------------------

void AliSelectorLYZ::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliDebug(AliLog::kDebug, "=======SLAVETERMINATE=======");

  DeleteKinematicsFile();
}

//-----------------------------------------------------------------------

void AliSelectorLYZ::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliDebug(AliLog::kDebug, "=========TERMINATE==========");

  cout << " Finished ... " << endl ;

  if(fDoNothing) { cout << "             ... & nothing was done ! " << endl ; cout << endl ; return ; }

  cout << endl ;
  cout << "  nTracks:  " << fFlowMaker->GetNgoodTracks() << endl ;   
  cout << "  nV0s:  " << fFlowMaker->GetNgoodV0s()  << endl ;  	     
  cout << "  nTracks (|eta|<0.5):  " << fFlowMaker->GetNgoodTracksEta() << endl ; 
  cout << "  nTracks+:  " << fFlowMaker->GetNposiTracks() << endl ; 	     
  cout << "  nTracks-:  " << fFlowMaker->GetNnegaTracks() << endl ; 	     
  cout << "  nTracks unconstrained:  " << fFlowMaker->GetNunconstrained() << endl ; 	 
  cout << "  Bayesian (e,mu,pi,k,p) :  " ; 
  for(int ii=0;ii<5;ii++) { cout << fFlowMaker->GetBayesianNorm(ii) << "   " ; } 
  cout << " . " << endl ; cout << endl ;
  
  if(fOnFlyAnalysis) 
    { 
      fFlowLYZ->Finish() ; cout << endl;
      delete fFlowLYZ;
    }
     
  if(fSaveFlowEvents) { fFlowfile->Close() ; cout << " file closed . " << endl ; }
  delete fFlowMaker ;

  delete fFlowSelect ;
  
  cout << endl ; 
  return ;
} 

//-----------------------------------------------------------------------

TTree* AliSelectorLYZ::GetKinematics()
{
  // Returns kinematics tree corresponding to current ESD active in fTree
  // Loads the kinematics from the kinematics file, the file is identified by replacing "AliESDs" to
  // "Kinematics" in the file path of the ESD file. This is a hack, to be changed!

  if (!fKineFile)
  {
    if(!fTree->GetCurrentFile()) { return 0 ; }

    TString fileName(fTree->GetCurrentFile()->GetName());
    fileName.ReplaceAll("AliESDs", "Kinematics");

    // temporary workaround for PROOF bug #18505
    fileName.ReplaceAll("#Kinematics.root#Kinematics.root", "#Kinematics.root");

    AliDebug(AliLog::kInfo, Form("Opening %s", fileName.Data()));

    fKineFile = TFile::Open(fileName);
    if(!fKineFile) { return 0 ; }
  }

  return dynamic_cast<TTree*> (fKineFile->Get(Form("Event%d/TreeK", fTree->GetTree()->GetReadEntry())));
}

//-----------------------------------------------------------------------

void AliSelectorLYZ::DeleteKinematicsFile()
{
  //
  // Closes the kinematics file and deletes the pointer.
  //

  if (fKineFile)
  {
    fKineFile->Close();
    delete fKineFile;
    fKineFile = 0;
  }
}

//-----------------------------------------------------------------------


