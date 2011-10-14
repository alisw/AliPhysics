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
#include "TInterpreter.h"

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
// extern AliFemtoManager *ConfigFemtoAnalysis();

//________________________________________________________________________
AliAnalysisTaskFemto::AliAnalysisTaskFemto(const char *name, const char *aConfigMacro, const char *aConfigParams):
    AliAnalysisTask(name,""), 
    fESD(0), 
    fESDpid(0),
    fAOD(0),
    fAODpidUtil(0),
    fStack(0),
    fOutputList(0), 
    fReader(0x0),
    fManager(0x0),
    fAnalysisType(0),
    fConfigMacro(0),
    fConfigParams(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TList::Class());
  fConfigMacro = (char *) malloc(sizeof(char) * strlen(aConfigMacro));
  strcpy(fConfigMacro, aConfigMacro);
  fConfigParams = (char *) malloc(sizeof(char) * strlen(aConfigParams));
  strcpy(fConfigParams, aConfigParams);
}
//________________________________________________________________________
AliAnalysisTaskFemto::AliAnalysisTaskFemto(const char *name, const char *aConfigMacro="ConfigFemtoAnalysis.C"): 
    AliAnalysisTask(name,""), 
    fESD(0), 
    fESDpid(0),
    fAOD(0),
    fAODpidUtil(0),
    fStack(0),
    fOutputList(0), 
    fReader(0x0),
    fManager(0x0),
    fAnalysisType(0),
    fConfigMacro(0),
    fConfigParams(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TList::Class());
  fConfigMacro = (char *) malloc(sizeof(char) * strlen(aConfigMacro));
  strcpy(fConfigMacro, aConfigMacro);
  fConfigParams = (char *) malloc(sizeof(char) * 2);
  strcpy(fConfigParams, "");
}

AliAnalysisTaskFemto::AliAnalysisTaskFemto(const AliAnalysisTaskFemto& aFemtoTask):
    AliAnalysisTask(aFemtoTask), 
    fESD(0), 
    fESDpid(0),
    fAOD(0),
    fAODpidUtil(0),
    fStack(0),
    fOutputList(0), 
    fReader(0x0),
    fManager(0x0),
    fAnalysisType(0),
    fConfigMacro(0),
    fConfigParams(0)
{
  // copy constructor
  fESD = aFemtoTask.fESD; 
  fESDpid = aFemtoTask.fESDpid; 
  fAOD = aFemtoTask.fAOD; 
  fAODpidUtil = aFemtoTask.fAODpidUtil;
  fStack = aFemtoTask.fStack;
  fOutputList = aFemtoTask.fOutputList;   
  fReader = aFemtoTask.fReader;       
  fManager = aFemtoTask.fManager;      
  fAnalysisType = aFemtoTask.fAnalysisType; 
  fConfigMacro = (char *) malloc(sizeof(char) * strlen(aFemtoTask.fConfigMacro));
  strcpy(fConfigMacro, aFemtoTask.fConfigMacro);
  fConfigParams = (char *) malloc(sizeof(char) * strlen(aFemtoTask.fConfigParams));
  strcpy(fConfigParams, aFemtoTask.fConfigParams);
}


AliAnalysisTaskFemto& AliAnalysisTaskFemto::operator=(const AliAnalysisTaskFemto& aFemtoTask){
  // assignment operator
  if (this == &aFemtoTask)
    return *this;

  fESD = aFemtoTask.fESD; 
  fESDpid = aFemtoTask.fESDpid;
  fAOD = aFemtoTask.fAOD; 
  fAODpidUtil = aFemtoTask.fAODpidUtil;
  fStack = aFemtoTask.fStack;
  fOutputList = aFemtoTask.fOutputList;   
  fReader = aFemtoTask.fReader;       
  fManager = aFemtoTask.fManager;      
  fAnalysisType = aFemtoTask.fAnalysisType; 
  if (fConfigMacro) free(fConfigMacro);
  fConfigMacro = (char *) malloc(sizeof(char) * strlen(aFemtoTask.fConfigMacro));
  strcpy(fConfigMacro, aFemtoTask.fConfigMacro);
  if (fConfigParams) free(fConfigParams);
  fConfigParams = (char *) malloc(sizeof(char) * strlen(aFemtoTask.fConfigParams));
  strcpy(fConfigParams, aFemtoTask.fConfigParams);

  return *this;
}

AliAnalysisTaskFemto::~AliAnalysisTaskFemto() 
{
  if (fConfigMacro) free(fConfigMacro);
  if (fConfigParams) free(fConfigParams);
}


//________________________________________________________________________
void AliAnalysisTaskFemto::ConnectInputData(Option_t *) {
  AliInfo(Form("   ConnectInputData %s\n", GetName()));

  fESD = 0;
  fESDpid = 0;
  fAOD = 0;
  fAODpidUtil = 0;
  fAnalysisType = 0;

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    AliWarning("Could not read chain from input slot 0");
    return;
  } 

  AliFemtoEventReaderESDChain *femtoReader = dynamic_cast<AliFemtoEventReaderESDChain *> (fReader);
  if ((dynamic_cast<AliFemtoEventReaderESDChain *> (fReader))) {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if(esdH) {
      AliInfo("Selected ESD analysis");
      fAnalysisType = 1;
      
//       if (!esdH) {
// 	AliWarning("Could not get ESDInputHandler");
//       } 
//       else {
	fESD = esdH->GetEvent();
        fESDpid = esdH->GetESDpid();
        femtoReader->SetESDPid(fESDpid);
//       }
    }
}
  else if ((dynamic_cast<AliFemtoEventReaderKinematicsChain *> (fReader))) {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if(esdH) {
      AliInfo("Selected ESD analysis");
      fAnalysisType = 1;
      
//       if (!esdH) {
// 	AliWarning("Could not get ESDInputHandler");
//       } 
//       else {
	fESD = esdH->GetEvent();
        //fESDpid = esdH->GetESDpid();
        //femtoReader->SetESDPid(fESDpid);
//       }
    }
  }


  AliFemtoEventReaderESDChainKine *femtoReaderESDKine = dynamic_cast<AliFemtoEventReaderESDChainKine *> (fReader);
  if ((dynamic_cast<AliFemtoEventReaderESDChainKine *> (fReader))) {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if(esdH) {
      AliInfo("Selected ESD analysis");
      fAnalysisType = 1;
      
//       if (!esdH) {
// 	AliWarning("Could not get ESDInputHandler");
//       } 
//       else {
	fESD = esdH->GetEvent();
        fESDpid = esdH->GetESDpid();
        femtoReaderESDKine->SetESDPid(fESDpid);
//       }
    }
}

  AliFemtoEventReaderKinematicsChain *femtoReaderKine = dynamic_cast<AliFemtoEventReaderKinematicsChain *> (fReader);
if ((dynamic_cast<AliFemtoEventReaderKinematicsChain *> (fReader))) {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if(esdH) {
      AliInfo("Selected ESD analysis");
      fAnalysisType = 1;
      
//       if (!esdH) {
// 	AliWarning("Could not get ESDInputHandler");
//       } 
//       else {
	fESD = esdH->GetEvent();
        //fESDpid = esdH->GetESDpid();
        //femtoReader->SetESDPid(fESDpid);
//       }
    }
 }

  
    AliFemtoEventReaderAODChain *femtoReaderAOD = dynamic_cast<AliFemtoEventReaderAODChain *> (fReader);
  if (dynamic_cast<AliFemtoEventReaderAODChain *> (fReader)) {
    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	
    if (!aodH) {
      TObject *handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
      AliInfo("Has output handler ");
      if( handler && handler->InheritsFrom("AliAODHandler") ) {
	AliInfo("Selected AOD analysis");

	fAOD = ((AliAODHandler*)handler)->GetAOD();
	fAnalysisType = 2;
      }
      else {
	AliWarning("Selected AOD reader but no AOD handler found");
      }
    } 
    else {
      AliInfo("Selected AOD analysis");
      fAnalysisType = 2;
      
      fAOD = aodH->GetEvent();

      fAODpidUtil = aodH->GetAODpidUtil();
      //      printf("aodH->GetAODpidUtil(): %x",aodH->GetAODpidUtil());
      femtoReaderAOD->SetAODpidUtil(fAODpidUtil);
    }
  }

  if ((!fAOD) && (!fESD)) {
    AliWarning("Wrong analysis type: Only ESD and AOD types are allowed!");
  }
}

//________________________________________________________________________
void AliAnalysisTaskFemto::CreateOutputObjects() {
  AliInfo("Creating Femto Analysis objects\n");

  gSystem->SetIncludePath("-I$ROOTSYS/include -I./STEERBase/ -I./ESD/ -I./AOD/ -I./ANALYSIS/ -I./ANALYSISalice/ -I./PWG2AOD/AOD -I./PWG2femtoscopy/FEMTOSCOPY/AliFemto -I./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser");
  //  char fcm[2000];
//   sprintf(fcm, "%s++", fConfigMacro);
//   gROOT->LoadMacro(fcm);
  gROOT->LoadMacro(fConfigMacro);
  //  fJetFinder = (AliJetFinder*) gInterpreter->ProcessLine("ConfigJetAnalysis()");
  if (!fConfigParams)
    SetFemtoManager((AliFemtoManager *) gInterpreter->ProcessLine("ConfigFemtoAnalysis()"));
  else
    SetFemtoManager((AliFemtoManager *) gInterpreter->ProcessLine(Form("ConfigFemtoAnalysis(%s)", fConfigParams)));

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

  PostData(0, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskFemto::Exec(Option_t *) {
  // Task making a femtoscopic analysis.

  if (fAnalysisType==1) {
    if (!fESD) {
      AliWarning("fESD not available");
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
	
	//	AliInfo ("Got MC cocktail event header %p\n", (void *) hd);
	TList *lhd = hd->GetHeaders();
	//	AliInfo ("Got list of headers %d\n", lhd->GetEntries());
	
	for (int iterh=0; iterh<lhd->GetEntries(); iterh++) 
	  {
	    hdh = dynamic_cast<AliGenHijingEventHeader *> (lhd->At(iterh));
	    //	    AliInfo ("HIJING header at %i is %p\n", iterh, (void *) hdh);
	  }
      }    
    }

    // Get ESD
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if (!esdH) {
      AliWarning("Could not get ESDInputHandler");
      return;
    } 
    else {
      fESD = esdH->GetEvent();
      fESDpid = esdH->GetESDpid();   
    }

    AliInfo(Form("Tracks in ESD: %d \n",fESD->GetNumberOfTracks()));

    if (fESD->GetNumberOfTracks() >= 0) {
    
      if (!fReader) {
 	AliWarning("No ESD reader for ESD analysis !\n");
      }
      
      AliFemtoEventReaderESDChain* fesdc = dynamic_cast<AliFemtoEventReaderESDChain *> (fReader);
      if (fesdc)
	{
	  // Process the event with no Kine information
	  fesdc->SetESDSource(fESD);
	  fManager->ProcessEvent();
	}

	AliFemtoEventReaderKinematicsChain* fkinec = dynamic_cast<AliFemtoEventReaderKinematicsChain *> (fReader);
      if (fkinec)
	{
	  // Process the event with Kine information only
	  fkinec->SetStackSource(fStack);
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
      AliFemtoEventReaderStandard* fstd = dynamic_cast<AliFemtoEventReaderStandard *> (fReader);
      if (fstd) 
	{
	  // Process the event with Kine information
	  fstd->SetESDSource(fESD);
	  if (mctruth) {
	    fstd->SetStackSource(fStack);
	    fstd->SetGenEventHeader(hdh);
	    fstd->SetInputType(AliFemtoEventReaderStandard::kESDKine);
	  }
	  else
	    fstd->SetInputType(AliFemtoEventReaderStandard::kESD);
	  fManager->ProcessEvent();
	}
    } 

    // Post the output histogram list
    PostData(0, fOutputList);
  }

  if (fAnalysisType==2) {    
    if (!fAOD) {
      AliWarning("fAOD not available");
      return;
    }

    // Get AOD
//     AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      
//     if (!aodH) {
//       AliWarning("Could not get AODInputHandler");
//       return;
//     } 
//     else {

//       fAOD = aodH->GetEvent();
//     }

    AliInfo(Form("Tracks in AOD: %d \n",fAOD->GetNumberOfTracks()));
    
    if (fAOD->GetNumberOfTracks() > 0) {
      if (!fReader) {
	AliWarning("No AOD reader for AOD analysis! \n");
      }
      else {
	AliFemtoEventReaderAODChain* faodc = dynamic_cast<AliFemtoEventReaderAODChain *> (fReader);

	if (faodc) {
	  // Process the event
	  faodc->SetAODSource(fAOD);
	  fManager->ProcessEvent();
	}
	AliFemtoEventReaderStandard* fstd = dynamic_cast<AliFemtoEventReaderStandard *> (fReader);

	if (fstd) {
	  // Process the event
	  fstd->SetAODSource(fAOD);
	  fstd->SetInputType(AliFemtoEventReaderStandard::kAOD);
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
  AliInfo("Selecting Femto reader for ESD\n");
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoReaderESDKine(AliFemtoEventReaderESDChainKine *aReader)
{
  AliInfo("Selecting Femto reader for ESD with Kinematics information\n");
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoReaderAOD(AliFemtoEventReaderAODChain *aReader)
{
  AliInfo("Selecting Femto reader for AOD\n");
  fReader = aReader;
}
void AliAnalysisTaskFemto::SetFemtoReaderStandard(AliFemtoEventReaderStandard *aReader)
{
  AliInfo("Selecting Standard all-purpose Femto reader\n");
  fReader = aReader;
}
void AliAnalysisTaskFemto::SetFemtoReaderKinematics(AliFemtoEventReaderKinematicsChain *aReader)
{
  printf("Selecting Femto reader for Kinematics (Monte Carlo) information\n");
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoManager(AliFemtoManager *aManager)
{
  fManager = aManager;
  AliInfo(Form("Got reader %p\n", (void *) aManager->EventReader()));
  AliFemtoEventReaderESDChain     *tReaderESDChain     = dynamic_cast<AliFemtoEventReaderESDChain *> (aManager->EventReader());
  AliFemtoEventReaderESDChainKine *tReaderESDChainKine = dynamic_cast<AliFemtoEventReaderESDChainKine *> (aManager->EventReader());
  AliFemtoEventReaderAODChain     *tReaderAODChain     = dynamic_cast<AliFemtoEventReaderAODChain *> (aManager->EventReader());
  AliFemtoEventReaderStandard     *tReaderStandard     = dynamic_cast<AliFemtoEventReaderStandard *> (aManager->EventReader());
  AliFemtoEventReaderKinematicsChain *tReaderKineChain = dynamic_cast<AliFemtoEventReaderKinematicsChain *> (aManager->EventReader());

 if ((!tReaderESDChain) && (!tReaderESDChainKine) && (!tReaderAODChain) && (!tReaderStandard) && (!tReaderKineChain)) {
    AliWarning("No AliFemto event reader created. Will not run femto analysis.\n");
    return;
  }
  if (tReaderESDChain) SetFemtoReaderESD(tReaderESDChain);
  if (tReaderESDChainKine) SetFemtoReaderESDKine(tReaderESDChainKine);
  if (tReaderAODChain) SetFemtoReaderAOD(tReaderAODChain);
  if (tReaderStandard) SetFemtoReaderStandard(tReaderStandard);
  if (tReaderKineChain) SetFemtoReaderKinematics(tReaderKineChain);
}

