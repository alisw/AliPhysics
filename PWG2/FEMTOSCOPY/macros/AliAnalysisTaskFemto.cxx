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

ClassImp(AliAnalysisTaskFemto)

//________________________________________________________________________
  AliAnalysisTaskFemto::AliAnalysisTaskFemto(const char *name): 
    AliAnalysisTask(name,""), 
    fESD(0), 
    fAOD(0),
    fOutputList(0), 
    fReaderESD(0x0),
    fReaderAOD(0x0),
    fManager(0x0),
    fAnalysisType(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TList::Class());
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
}

//________________________________________________________________________
void AliAnalysisTaskFemto::Exec(Option_t *) {
  // Task making a femtoscopic analysis.
  if (fAnalysisType==1) {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }

    printf("Tracks in ESD: %d \n",fESD->GetNumberOfTracks());
  
    if (fESD->GetNumberOfTracks() >= 0) {
    
      if (!fReaderESD) {
	printf("ERROR: No ESD reader for ESD analysis !\n");
      }
      else {
	fReaderESD->SetESDSource(fESD);
	fManager->ProcessEvent();
      }
    } 
    fOutputList = fManager->Analysis(0)->GetOutputList();
    PostData(0, fOutputList);
  }
  
  if (fAnalysisType==2) {    
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    }
    printf("Tracks in AOD: %d \n",fAOD->GetNumberOfTracks());
    
    if (fAOD->GetNumberOfTracks() > 0) {
      Double_t pxyz[3];
      
      if (!fReaderAOD) {
	printf("ERROR: No AOD reader for AOD analysis! \n");
      }
      else {
	fReaderAOD->SetAODSource(fAOD);
	fManager->ProcessEvent();
      }
    } 
    fOutputList = fManager->Analysis(0)->GetOutputList();
    PostData(0, fOutputList);
  }
}      

//________________________________________________________________________
void AliAnalysisTaskFemto::Terminate(Option_t *) {
  // Do the final processing
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoReaderESD(AliFemtoEventReaderESDChain *aReader)
{
  fReaderESD = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoReaderAOD(AliFemtoEventReaderAODChain *aReader)
{
  fReaderAOD = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoManager(AliFemtoManager *aManager)
{
  fManager = aManager;
}

