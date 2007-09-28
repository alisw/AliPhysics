#define AliAnalysisTaskRLLYZ_cxx
 
#include <iostream>


#include "TChain.h"
#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TParticle.h"

#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliStack.h"
#include <AliHeader.h>
#include <AliGenEventHeader.h>

#include "AliAnalysisTaskRL.h"
#include "AliAnalysisTaskRLLYZ.h"
#include "AliFlowConstants.h"
#include "AliFlowLYZConstants.h"
#include "AliFlowSelection.h"
#include "AliFlowEvent.h"
#include "AliFlowMaker.h"
#include "AliFlowTrack.h"
#include "AliFlowLeeYangZerosMaker.h"

//#include "TObjectTable.h"



ClassImp(AliAnalysisTaskRLLYZ)

 //-----------------------------------------------------------------------
 
 AliAnalysisTaskRLLYZ::AliAnalysisTaskRLLYZ(const char *name, Bool_t firstrun) :
   AliAnalysisTaskRL(name,""),
   fESD(0),
   fFirstRunLYZ(firstrun), //set boolean for firstrun to initial value
   fUseSumLYZ(kTRUE)    //set boolean for use sum to initial value
{

  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  if (!firstrun) DefineInput(1, TList::Class()); //for second loop 
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  
}

 
 //-----------------------------------------------------------------------


 AliAnalysisTaskRLLYZ::~AliAnalysisTaskRLLYZ() 
 {
   //destructor
   
 }
 
//-----------------------------------------------------------------------


void AliAnalysisTaskRLLYZ::ConnectInputData(Option_t *) {
  // Initialize branches.
  printf("   ConnectInputData of task %s\n", GetName());
  //cerr<<"fESD ("<<fESD<<")"<<endl;
  if (!fESD) {
    //cerr<<"no fESD"<<endl;
    char ** address = (char **)GetBranchAddress(0, "ESD");
    if (address) fESD = (AliESD*)(*address);
    if (!fESD) {
      //cerr<<"still no fESD"<<endl;
      fESD = new AliESD();
      SetBranchAddress(0, "ESD", &fESD);
      cerr<<"new fESD"<<endl;
    }
  }
}

//-----------------------------------------------------------------------
void AliAnalysisTaskRLLYZ::CreateOutputObjects() {

  
  fFlowMaker = new AliFlowMaker();
  cerr<<"create fFlowMaker ("<<fFlowMaker<<")"<<endl;
  fFlowMaker->SetNHitsCut(1);
  fFlowMaker->SetECut(0.01,100.);
  fFlowMaker->PrintCutList();

  fFlowSelect = new AliFlowSelection();
  cerr<<"create fFlowSelect ("<<fFlowSelect<<")"<<endl;
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

  //AliFlowLYZAnalyser...
  fFlowLYZ = new AliFlowLeeYangZerosMaker(fFlowSelect);
  cerr<<"fFlowLYZ ("<<fFlowLYZ<<")"<<endl;
  fFlowLYZ -> SetFirstRun(GetFirstRunLYZ());   //set first run true or false
  fFlowLYZ -> SetUseSum(GetUseSumLYZ());     //set use sum true or false
  //fFlowLYZ -> SetDebug(kTRUE) ;     //for debugging porposes
  
     
  //output file
  if (fFirstRunLYZ) fFlowLYZ->SetHistFileName("testTaskLYZ_firstrun.root");
  else fFlowLYZ->SetHistFileName("testTaskLYZ_secondrun.root");


  // Get data from input slot 1
  if (GetNinputs() == 2) {                   //if there are two input slots
    fFirstRunFile = (TFile*)GetInputData(1);
    cerr<<"fFirstRunFile ("<<fFirstRunFile<<")"<<endl;
    cerr<<"fFirstRunFile -> IsOpen() = "<<fFirstRunFile -> IsOpen()<<endl;

    fFlowLYZ -> SetFirstRunFile(fFirstRunFile);
  }
  fFlowLYZ -> Init() ; cerr<<"fFlowLYZ->Init()"<<endl;
      
} 
 
//-----------------------------------------------------------------------
 
void AliAnalysisTaskRLLYZ::Exec(Option_t *) {

  
  // Get data from input slot 0
  TTree *tinput = (TTree*)GetInputData(0);
  Long64_t ientry = tinput->GetReadEntry();
  if (AliAnalysisTaskRL::GetEntry(ientry) == kFALSE) {
    printf("Couldn't get event from the runLoader\n");
    return;
  }
  
  if (!fESD) {
    cout << "No ESD branch available" << endl;
    return;
  }
  
  cerr<<"fESD ("<<fESD <<") in begin Exec"<<endl;
  cerr<<"number of tracks: "<<fESD->GetNumberOfTracks()<<endl;

  AliStack* stack = GetStack();
  if (!stack) {
    AliDebug(AliLog::kError, "Stack not available");
    // return kFALSE;
  }

  cerr<<"fFlowMaker ("<<fFlowMaker<<")"<<endl;
  cerr<<"fFlowSelect ("<<fFlowSelect<<")"<<endl;
  
  fFlowEvent = new AliFlowEvent() ;
  cerr<<"create fFlowEvent ("<<fFlowEvent<<")"<<endl;
  
  if (!fFlowMaker){cerr<<"no fFlowMaker: nullpointer"<<endl;}
 
  else { 
    if (!fESD) { cerr<<"no fESD: NULL pointer"<<endl;}
    
    else {
      cerr<<"fFlowMaker ("<<fFlowMaker<<"), fFlowEvent ("<<fFlowEvent<<") and fESD ("<<fESD<<") available"<<endl;
      fFlowEvent = fFlowMaker->FillFlowEvent(fESD);
      
      if (!fFlowEvent){ cerr<<"no fFlowEvent: NULL pointer"<<endl; }
      else {
	if(fFlowSelect->Select(fFlowEvent)){	 // event selected 
	 
	  cerr<<"event selected"<<endl;
	  fFlowEvent->SetSelections(fFlowSelect) ;
	  cerr<<"fFlowLYZ ("<<fFlowLYZ<<")"<<endl;
	  fFlowLYZ->Make(fFlowEvent);
	}
      }
    }
  }
  
  delete fFlowEvent;
  cerr<<"delete fFlowEvent ("<<fFlowEvent<<")"<<endl;
  //post data
  //PostData(0,fFlowLYZ->GetHistFile());

}

  //--------------------------------------------------------------------    
void AliAnalysisTaskRLLYZ::Terminate(Option_t *) {
   
  if (GetNinputs() == 2) cerr<<"fFirstRunFile -> IsOpen() = "<<fFirstRunFile -> IsOpen()<<endl;
  cerr<<"fFlowLYZ->GetHistFile() -> IsOpen() = "<<fFlowLYZ->GetHistFile() -> IsOpen()<<endl;
  fFlowLYZ->Finish();
  PostData(0,fFlowLYZ->GetHistFile());

  delete fFlowLYZ;
  delete fFlowMaker;
  delete fFlowSelect;

  cout<<".....finished"<<endl;
 }

 

