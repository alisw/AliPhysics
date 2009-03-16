//------------------------------------------------------
// AliAnalysisTaskFemto - A task for the analysis framework
// from the FEMTOSCOPY analysis of PWG2. Creates the necessary
// connection between the ESD or AOD input and the femtoscopic
// code.
// Author: Adam Kisiel, OSU; Adam.Kisiel@cern.ch
//------------------------------------------------------
#include "TROOT.h"
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"

#include "AliAnalysisTask.h"

#include "AliESDEvent.h"

#include "AliFemtoAnalysis.h"
#include "AliAnalysisTaskFemto.h"
#include "AliVHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"

ClassImp(AliAnalysisTaskFemto)

// Default name for the setup macro of femto analysis  
// This function MUST be defined in the separate file !!!
extern AliFemtoManager *ConfigFemtoAnalysis();

//________________________________________________________________________
  AliAnalysisTaskFemto::AliAnalysisTaskFemto(const char *name): 
    AliAnalysisTask(name,""), 
    fESD(0), 
    fAOD(0),
    fStack(0),
    fOutputList(0), 
    fReader(0x0),
    fManager(0x0),
    fAnalysisType(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TList::Class());
}

AliAnalysisTaskFemto::AliAnalysisTaskFemto(const AliAnalysisTaskFemto& aFemtoTask){
  // copy constructor
  fESD = aFemtoTask.fESD; 
  fAOD = aFemtoTask.fAOD; 
  fStack = aFemtoTask.fStack;
  fOutputList = aFemtoTask.fOutputList;   
  fReader = aFemtoTask.fReader;       
  fManager = aFemtoTask.fManager;      
  fAnalysisType = aFemtoTask.fAnalysisType; 
}


AliAnalysisTaskFemto& AliAnalysisTaskFemto::operator=(const AliAnalysisTaskFemto& aFemtoTask){
  // assignment operator
  if (this == &aFemtoTask)
    return *this;

  fESD = aFemtoTask.fESD; 
  fAOD = aFemtoTask.fAOD; 
  fStack = aFemtoTask.fStack;
  fOutputList = aFemtoTask.fOutputList;   
  fReader = aFemtoTask.fReader;       
  fManager = aFemtoTask.fManager;      
  fAnalysisType = aFemtoTask.fAnalysisType; 

  return *this;
}

//________________________________________________________________________
void AliAnalysisTaskFemto::ConnectInputData(Option_t *) {
  printf("   ConnectInputData %s\n", GetName());

  fESD = 0;
  fAOD = 0;
  fAnalysisType = 0;

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if(esdH) {
      cout << "Selected ESD analysis" << endl;
      fAnalysisType = 1;
      
      if (!esdH) {
	Printf("ERROR: Could not get ESDInputHandler");
      } 
      else {
	fESD = esdH->GetEvent();
      }
    }
    else {
      AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      
      if (!aodH) {
	Printf("ERROR: Could not get AODInputHandler");
      } 
      else {
	cout << "Selected AOD analysis" << endl;
	fAnalysisType = 2;

	fAOD = aodH->GetEvent();
      }
    }
    if ((!fAOD) && (!fESD)) {
      Printf("Wrong analysis type: Only ESD and AOD types are allowed!");
    }
  }
  
  
}

//________________________________________________________________________
void AliAnalysisTaskFemto::CreateOutputObjects() {
  printf("Creating Femto Analysis objects\n");

  SetFemtoManager(ConfigFemtoAnalysis());

  TList *tOL;
  fOutputList = fManager->Analysis(0)->GetOutputList();

  for (unsigned int ian = 1; ian<fManager->AnalysisCollection()->size(); ian++) {
    tOL = fManager->Analysis(ian)->GetOutputList();

    TIter nextListCf(tOL);
    while (TObject *obj = nextListCf()) {
      fOutputList->Add(obj);
    }

    delete tOL;
  }
}

//________________________________________________________________________
void AliAnalysisTaskFemto::Exec(Option_t *) {
  // Task making a femtoscopic analysis.

  if (fAnalysisType==1) {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }

    //Get MC data
    AliMCEventHandler*    mctruth = (AliMCEventHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    
    AliGenHijingEventHeader *hdh = 0;
    if(mctruth) {
      fStack = mctruth->MCEvent()->Stack();

      AliGenCocktailEventHeader *hd = dynamic_cast<AliGenCocktailEventHeader *> (mctruth->MCEvent()->GenEventHeader());
      
      if (hd) {
	
	printf ("Got MC cocktail event header %p\n", (void *) hd);
	TList *lhd = hd->GetHeaders();
	printf ("Got list of headers %d\n", lhd->GetEntries());
	
	for (int iterh=0; iterh<lhd->GetEntries(); iterh++) 
	  {
	    hdh = dynamic_cast<AliGenHijingEventHeader *> (lhd->At(iterh));
	    printf ("HIJING header at %i is %p\n", iterh, (void *) hdh);
	  }
      }    
    }

    printf("Tracks in ESD: %d \n",fESD->GetNumberOfTracks());

    if (fESD->GetNumberOfTracks() >= 0) {
    
      if (!fReader) {
 	printf("ERROR: No ESD reader for ESD analysis !\n");
      }
      
      AliFemtoEventReaderESDChain* fesdc = dynamic_cast<AliFemtoEventReaderESDChain *> (fReader);
      if (fesdc)
	{
	  // Process the event with no Kine information
	  fesdc->SetESDSource(fESD);
	  fManager->ProcessEvent();
	}
      AliFemtoEventReaderESDChainKine* fesdck = dynamic_cast<AliFemtoEventReaderESDChainKine *> (fReader);
      if (fesdck) 
	{
	  // Process the event with Kine information
	  fesdck->SetESDSource(fESD);
	  fesdck->SetStackSource(fStack);
	  
	  fesdck->SetGenEventHeader(hdh);
	  fManager->ProcessEvent();
	}
    } 

    // Post the output histogram list
    PostData(0, fOutputList);
  }
  
  if (fAnalysisType==2) {    
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    }
    printf("Tracks in AOD: %d \n",fAOD->GetNumberOfTracks());
    
    if (fAOD->GetNumberOfTracks() > 0) {
      if (!fReader) {
	printf("ERROR: No AOD reader for AOD analysis! \n");
      }
      else {
	AliFemtoEventReaderAODChain* faodc = dynamic_cast<AliFemtoEventReaderAODChain *> (fReader);

	if (faodc) {
	  // Process the event
	  faodc->SetAODSource(fAOD);
	  fManager->ProcessEvent();
	}
      }
    } 

    // Post the output histogram list
    PostData(0, fOutputList);
  }
}      

//________________________________________________________________________
void AliAnalysisTaskFemto::Terminate(Option_t *) {
  // Do the final processing
  if (fManager) {
    fManager->Finish();
  }
}
//________________________________________________________________________
void AliAnalysisTaskFemto:: FinishTaskOutput() {
  // Do the final processing
  if (fManager) {
    fManager->Finish();
  }
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoReaderESD(AliFemtoEventReaderESDChain *aReader)
{
  printf("Selectring Femto reader for ESD\n");
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoReaderESDKine(AliFemtoEventReaderESDChainKine *aReader)
{
  printf("Selectring Femto reader for ESD with Kinematics information\n");
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoReaderAOD(AliFemtoEventReaderAODChain *aReader)
{
  printf("Selecting Femto reader for AOD\n");
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoManager(AliFemtoManager *aManager)
{
  fManager = aManager;
  printf("Got reader %p\n", (void *) aManager->EventReader());
  AliFemtoEventReaderESDChain     *tReaderESDChain     = dynamic_cast<AliFemtoEventReaderESDChain *> (aManager->EventReader());
  AliFemtoEventReaderESDChainKine *tReaderESDChainKine = dynamic_cast<AliFemtoEventReaderESDChainKine *> (aManager->EventReader());
  AliFemtoEventReaderAODChain     *tReaderAODChain     = dynamic_cast<AliFemtoEventReaderAODChain *> (aManager->EventReader());

  if ((!tReaderESDChain) && (!tReaderESDChainKine) && (!tReaderAODChain)) {
    printf("No AliFemto event reader created. Will not run femto analysis.\n");
    return;
  }
  if (tReaderESDChain) SetFemtoReaderESD(tReaderESDChain);
  if (tReaderESDChainKine) SetFemtoReaderESDKine(tReaderESDChainKine);
  if (tReaderAODChain) SetFemtoReaderAOD(tReaderAODChain);
}

