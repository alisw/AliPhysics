///
/// \file AliAnalysisTaskFemto.cxx
///

#include "TROOT.h"
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TInterpreter.h"

//#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"

#include "AliESDEvent.h"

#include "AliFemtoAnalysis.h"
#include "AliAnalysisTaskFemtoMJ.h"
#include "AliVHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliAnalysisTaskFemtoMJ);
  /// \endcond
#endif


// Default name for the setup macro of femto analysis
// This function MUST be defined in the separate file !!!
// extern AliFemtoManager *ConfigFemtoAnalysis();

//________________________________________________________________________
AliAnalysisTaskFemtoMJ::AliAnalysisTaskFemtoMJ(TString name, TString aConfigMacro, TString aConfigParams, Bool_t aVerbose):
  AliAnalysisTaskSE(name), //AliAnalysisTask(name,""),
  fESD(NULL),
  fESDpid(NULL),
  fAOD(NULL),
  fAODpidUtil(NULL),
  fAODheader(NULL),
  fStack(NULL),
  fOutputList(NULL),
  fReader(NULL),
  fManager(NULL),
  fAnalysisType(0),
  fConfigMacro(aConfigMacro),
  fConfigParams(aConfigParams),
  fVerbose(aVerbose)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  //DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TList::Class());

}
//________________________________________________________________________
AliAnalysisTaskFemtoMJ::AliAnalysisTaskFemtoMJ(TString name, TString aConfigMacro = "ConfigFemtoAnalysis.C", Bool_t aVerbose):
  AliAnalysisTaskSE(name), //AliAnalysisTask(name,""),
  fESD(NULL),
  fESDpid(NULL),
  fAOD(NULL),
  fAODpidUtil(NULL),
  fAODheader(NULL),
  fStack(NULL),
  fOutputList(NULL),
  fReader(NULL),
  fManager(NULL),
  fAnalysisType(0),
  fConfigMacro(aConfigMacro),
  fConfigParams(""),
  fVerbose(aVerbose)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  //DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TList::Class());

}

AliAnalysisTaskFemtoMJ::AliAnalysisTaskFemtoMJ(const AliAnalysisTaskFemtoMJ &aFemtoTask):
  AliAnalysisTaskSE(aFemtoTask), //AliAnalysisTask(aFemtoTask),
  fESD(aFemtoTask.fESD),
  fESDpid(aFemtoTask.fESDpid),
  fAOD(aFemtoTask.fAOD),
  fAODpidUtil(aFemtoTask.fAODpidUtil),
  fAODheader(aFemtoTask.fAODheader),
  fStack(aFemtoTask.fStack),
  fOutputList(aFemtoTask.fOutputList),
  fReader(aFemtoTask.fReader),
  fManager(aFemtoTask.fManager),
  fAnalysisType(aFemtoTask.fAnalysisType),
  fConfigMacro(aFemtoTask.fConfigMacro),
  fConfigParams(aFemtoTask.fConfigParams),
  fVerbose(aFemtoTask.fVerbose)
{
  // copy constructor
}


AliAnalysisTaskFemtoMJ &AliAnalysisTaskFemtoMJ::operator=(const AliAnalysisTaskFemtoMJ &aFemtoTask)
{
  // assignment operator
  if (this == &aFemtoTask)
    return *this;

  fESD = aFemtoTask.fESD;
  fESDpid = aFemtoTask.fESDpid;
  fAOD = aFemtoTask.fAOD;
  fAODpidUtil = aFemtoTask.fAODpidUtil;
  fAODheader = aFemtoTask.fAODheader;
  fStack = aFemtoTask.fStack;
  fOutputList = aFemtoTask.fOutputList;
  fReader = aFemtoTask.fReader;
  fManager = aFemtoTask.fManager;
  fAnalysisType = aFemtoTask.fAnalysisType;

  fConfigMacro = aFemtoTask.fConfigMacro;
  fConfigParams = aFemtoTask.fConfigParams;
  fVerbose = aFemtoTask.fVerbose;

  return *this;
}

AliAnalysisTaskFemtoMJ::~AliAnalysisTaskFemtoMJ()
{
}


//________________________________________________________________________
void AliAnalysisTaskFemtoMJ::ConnectInputData(Option_t *)
{
  AliInfo(Form("   ConnectInputData %s\n", GetName()));

  fESD = 0;
  fESDpid = 0;
  fAOD = 0;
  fAODpidUtil = 0;
  fAODheader = 0;
  fAnalysisType = 0;

  TTree *tree = dynamic_cast<TTree *>(GetInputData(0));
  if (!tree) {
    AliWarning("Could not read chain from input slot 0");
    return;
  }

  AliFemtoEventReaderESDChain *femtoReader = dynamic_cast<AliFemtoEventReaderESDChain *>(fReader);
  if ((dynamic_cast<AliFemtoEventReaderESDChain *>(fReader))) {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (esdH) {
      if (fVerbose)
        AliInfo("Selected ESD analysis");
      fAnalysisType = 1;

      fESD = (AliESDEvent*)esdH->GetEvent();
      fESDpid = esdH->GetESDpid();
      femtoReader->SetESDPid(fESDpid);
    }
  } else if ((dynamic_cast<AliFemtoEventReaderKinematicsChain *>(fReader))) {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (esdH) {
      if (fVerbose)
        AliInfo("Selected ESD analysis");
      fAnalysisType = 1;
      fESD = (AliESDEvent*)esdH->GetEvent();
    }
  } else if ((dynamic_cast<AliFemtoEventReaderKinematicsChainESD *>(fReader))) {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (esdH) {
      if (fVerbose)
        AliInfo("Selected ESD analysis");
      fAnalysisType = 1;
      fESD = (AliESDEvent*)esdH->GetEvent();
    }
  }

  AliFemtoEventReaderESDChainKine *femtoReaderESDKine = dynamic_cast<AliFemtoEventReaderESDChainKine *>(fReader);
  if ((dynamic_cast<AliFemtoEventReaderESDChainKine *>(fReader))) {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (esdH) {
      if (fVerbose)
        AliInfo("Selected ESD analysis");
      fAnalysisType = 1;
      fESD = (AliESDEvent*)esdH->GetEvent();
      fESDpid = esdH->GetESDpid();
      femtoReaderESDKine->SetESDPid(fESDpid);

    }
  }

  //  AliFemtoEventReaderKinematicsChain *femtoReaderKine = dynamic_cast<AliFemtoEventReaderKinematicsChain *> (fReader);
  if ((dynamic_cast<AliFemtoEventReaderKinematicsChain *>(fReader))) {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (esdH) {
      if (fVerbose)
        AliInfo("Selected ESD analysis");
      fAnalysisType = 1;

      //       if (!esdH) {
      //  AliWarning("Could not get ESDInputHandler");
      //       }
      //       else {
      fESD = (AliESDEvent*)esdH->GetEvent();
      //fESDpid = esdH->GetESDpid();
      //femtoReader->SetESDPid(fESDpid);
      //       }
    }
  }


  AliFemtoEventReaderAODChain *femtoReaderAOD = dynamic_cast<AliFemtoEventReaderAODChain *>(fReader);
  if (dynamic_cast<AliFemtoEventReaderAODChain *>(fReader)) {
    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!aodH) {
      TObject *handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
      if (fVerbose)
        AliInfo("Has output handler ");
      if (handler && handler->InheritsFrom("AliAODHandler")) {
        if (fVerbose)
          AliInfo("Selected AOD analysis");

        fAOD = ((AliAODHandler *)handler)->GetAOD();
        fAnalysisType = 2;
      } else {
        if (fVerbose)
          AliWarning("Selected AOD reader but no AOD handler found");
      }
    } else {
      if (fVerbose)
        AliInfo("Selected AOD analysis");
      fAnalysisType = 2;

      fAOD = aodH->GetEvent();

      fAODpidUtil = aodH->GetAODpidUtil(); //correct way
      //fAODpidUtil = new AliAODpidUtil(); //not correct way
      //      printf("aodH->GetAODpidUtil(): %x",aodH->GetAODpidUtil());
      if (fVerbose)
        cout << "AliAnalysisTaskFemtoMJ::AodpidUtil:" << fAODpidUtil << endl;
      femtoReaderAOD->SetAODpidUtil(fAODpidUtil);

      fAODheader = dynamic_cast<AliAODHeader *>(fAOD->GetHeader());
      if (!fAODheader) AliFatal("Not a standard AOD");
      femtoReaderAOD->SetAODheader(fAODheader);

    }
  }

  if ((!fAOD) && (!fESD)) {
    if (fVerbose)
      AliWarning("Wrong analysis type: Only ESD and AOD types are allowed!");
  }
}

//________________________________________________________________________
void AliAnalysisTaskFemtoMJ::CreateOutputObjects()
{
  if (fVerbose)
    AliInfo("Creating Femto Analysis objects\n");

  gSystem->SetIncludePath("-I$ROOTSYS/include -I./STEERBase/ -I./ESD/ -I./AOD/ -I./ANALYSIS/ -I./ANALYSISalice/ -I./PWG2AOD/AOD -I./PWG2femtoscopy/FEMTOSCOPY/AliFemto -I./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser");
  //  char fcm[2000];
//   sprintf(fcm, "%s++", fConfigMacro);
//   gROOT->LoadMacro(fcm);
  cout<< "AliAnalysisTaskFemtoMJ: Loading config of the analysis: "<<fConfigMacro<<endl; //!!!
  gROOT->LoadMacro(fConfigMacro);
  //  fJetFinder = (AliJetFinder*) gInterpreter->ProcessLine("ConfigJetAnalysis()");
  if (!fConfigParams)
    SetFemtoManager((AliFemtoManager *) gInterpreter->ProcessLine("ConfigFemtoAnalysis()"));
  else
    SetFemtoManager((AliFemtoManager *) gInterpreter->ProcessLine(Form("ConfigFemtoAnalysis(%s)", fConfigParams.Data())));

  fOutputList = fManager->Analysis(0)->GetOutputList();
  fOutputList->SetOwner(kTRUE);

  for (UInt_t ian = 1; ian < fManager->AnalysisCollection()->size(); ian++) {
    TList* tOL = fManager->Analysis(ian)->GetOutputList();

    TIter nextListCf(tOL);
    while (TObject *obj = nextListCf()) {
      fOutputList->Add(obj);
    }

    delete tOL;
  }

  PostData(0, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskFemtoMJ::Exec(Option_t *)
{
  //cout<<"AliAnalysisTaskFemtoMJ::Exec()"<<endl;
  // Task making a femtoscopic analysis.
  if (fOfflineTriggerMask) {
    Bool_t isSelected = (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fOfflineTriggerMask);
    if (!isSelected) {
      if (fVerbose)
	cout << "AliAnalysisTaskFemtoMJ: is not selected" << endl;
      //return;//!!!
    }
  }

  //cout<<"AliAnalysisTaskFemtoMJ:: retreiving MCEvent()"<<endl;//!!!



  //jedna metoda

  /*AliMCEvent* mcEvent = MCEvent(); //!!!
    if (!mcEvent) {//!!!
      Printf("ERROR: Could not retrieve MC event");//!!!
      return;//!!!
    }//!!!

    AliStack *stack = mcEvent->Stack();//!!!

    AliMCEventHandler *mctruth = (AliMCEventHandler *)
    ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());

    AliGenHijingEventHeader *hdh = 0;*/

    //druga metoda
    AliMCEventHandler *mctruth = (AliMCEventHandler *)
                                    ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());


    if(!mctruth)
      {
	Printf("ERROR: Could not retrieve MC handler");//!!!
      }

    AliMCEvent* mcEvent = mctruth->MCEvent();
    if (!mcEvent) {//!!!
      Printf("ERROR: Could not retrieve MC event");//!!!
      return;//!!!
    }//!!!

      fStack = mcEvent->Stack();//!!!

      AliGenHijingEventHeader *hdh = 0;
	  AliGenEventHeader *header = 0;








    if (mctruth) {
    //fStack = mctruth->MCEvent()->Stack(); //tak bylo//!!!


      AliGenCocktailEventHeader *hd = dynamic_cast<AliGenCocktailEventHeader *>(mctruth->MCEvent()->GenEventHeader());
	  header = dynamic_cast<AliGenEventHeader *>(mctruth->MCEvent()->Header()->GenEventHeader());


      if (hd) {

        //  AliInfo ("Got MC cocktail event header %p\n", (void *) hd);
        TList *lhd = hd->GetHeaders();
        //  AliInfo ("Got list of headers %d\n", lhd->GetEntries());

        for (int iterh = 0; iterh < lhd->GetEntries(); iterh++) {
          hdh = dynamic_cast<AliGenHijingEventHeader *>(lhd->At(iterh));
          //      AliInfo ("HIJING header at %i is %p\n", iterh, (void *) hdh);
        }
      }
    }


    AliFemtoEventReaderKinematicsChain *fkinec = dynamic_cast<AliFemtoEventReaderKinematicsChain *>(fReader);
    if (fkinec) {
      // Process the event with Kine information only
      fkinec->SetStackSource(fStack);
      fkinec->SetGenEventHeader(header);
      fManager->ProcessEvent();
    }


    PostData(0, fOutputList);

  if (fAnalysisType == 1) {
    if (!fESD) {
      if (fVerbose)
	AliWarning("fESD not available");
      //return;//!!!-
    }
    //Get MC data


    // Get ESD
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      if (fVerbose)
        AliWarning("Could not get ESDInputHandler");
      return;
    } else {
      fESD = (AliESDEvent*)esdH->GetEvent();
      fESDpid = esdH->GetESDpid();
    }

    if (fVerbose)
      AliInfo(Form("Tracks in ESD: %d \n", fESD->GetNumberOfTracks()));

    if (fESD->GetNumberOfTracks() >= 0) {

      if (!fReader && fVerbose) {
          AliWarning("No ESD reader for ESD analysis !\n");
      }

      AliFemtoEventReaderESDChain *fesdc = dynamic_cast<AliFemtoEventReaderESDChain *>(fReader);
      if (fesdc) {
        // Process the event with no Kine information
        fesdc->SetESDSource(fESD);
        fManager->ProcessEvent();
      }
    }



    AliFemtoEventReaderESDChainKine *fesdck = dynamic_cast<AliFemtoEventReaderESDChainKine *>(fReader);
    if (fesdck) {
      // Process the event with Kine information
      fesdck->SetESDSource(fESD);
      fesdck->SetStackSource(fStack);
      cout << "Set Stack:" << fStack << endl;
      fesdck->SetGenEventHeader(hdh);
      fManager->ProcessEvent();
    }


    AliFemtoEventReaderKinematicsChainESD *fkcesd = dynamic_cast<AliFemtoEventReaderKinematicsChainESD *>(fReader);
    if (fkcesd) {
      // Process the event with Kine information
      fkcesd->SetESDSource(fESD);
      fkcesd->SetStackSource(fStack);
      fkcesd->SetGenEventHeader(hdh);
      fManager->ProcessEvent();
    }

    AliFemtoEventReaderStandard *fstd = dynamic_cast<AliFemtoEventReaderStandard *>(fReader);
    if (fstd) {
      // Process the event with Kine information
      fstd->SetESDSource(fESD);
      if (mctruth) {
        fstd->SetStackSource(fStack);
        fstd->SetGenEventHeader(hdh);
        fstd->SetInputType(AliFemtoEventReaderStandard::kESDKine);
      } else
        fstd->SetInputType(AliFemtoEventReaderStandard::kESD);
      fManager->ProcessEvent();
    }


    // Post the output histogram list
    PostData(0, fOutputList);
  }

  if (fAnalysisType == 2) {
    if (!fAOD) {
      if (fVerbose)
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




    if (fVerbose)
      AliInfo(Form("Tracks in AOD: %d \n", fAOD->GetNumberOfTracks()));

    if (fAOD->GetNumberOfTracks() > 0) {
      if (!fReader) {
        if (fVerbose)
          AliWarning("No AOD reader for AOD analysis! \n");
      } else {
        AliFemtoEventReaderAODChain *faodc = dynamic_cast<AliFemtoEventReaderAODChain *>(fReader);

        if (faodc) {
          // Process the event
          faodc->SetAODSource(fAOD);
          fManager->ProcessEvent();
        }
        AliFemtoEventReaderStandard *fstd = dynamic_cast<AliFemtoEventReaderStandard *>(fReader);

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
void AliAnalysisTaskFemtoMJ::Terminate(Option_t *)
{
  // Do the final processing
  if (fManager) {
    fManager->Finish();
  }
}
//________________________________________________________________________
void AliAnalysisTaskFemtoMJ:: FinishTaskOutput()
{
  // Do the final processing
  if (fManager) {
    fManager->Finish();
  }
}
//________________________________________________________________________
void AliAnalysisTaskFemtoMJ::SetFemtoReaderESD(AliFemtoEventReaderESDChain *aReader)
{
  if (fVerbose)
    AliInfo("Selecting Femto reader for ESD\n");
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemtoMJ::SetFemtoReaderESDKine(AliFemtoEventReaderESDChainKine *aReader)
{
  if (fVerbose)
    AliInfo("Selecting Femto reader for ESD with Kinematics information\n");
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemtoMJ::SetFemtoReaderAOD(AliFemtoEventReaderAODChain *aReader)
{
  if (fVerbose)
    AliInfo("Selecting Femto reader for AOD\n");
  fReader = aReader;
}
void AliAnalysisTaskFemtoMJ::SetFemtoReaderStandard(AliFemtoEventReaderStandard *aReader)
{
  if (fVerbose)
    AliInfo("Selecting Standard all-purpose Femto reader\n");
  fReader = aReader;
}
void AliAnalysisTaskFemtoMJ::SetFemtoReaderKinematics(AliFemtoEventReaderKinematicsChain *aReader)
{
  if (fVerbose)
    printf("Selecting Femto reader for Kinematics (Monte Carlo) information\n");
  fReader = aReader;
}
void AliAnalysisTaskFemtoMJ::SetFemtoReaderKinematicsESD(AliFemtoEventReaderKinematicsChainESD *aReader)
{
  if (fVerbose)
    printf("Selecting Femto reader for Kinematics (Monte Carlo) information + ESD\n");
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemtoMJ::SetFemtoManager(AliFemtoManager *aManager)
{
  fManager = aManager;

  AliFemtoEventReader *eventReader = aManager->EventReader();

  if (fVerbose)
    AliInfo(Form("Got reader %p\n", (void *)eventReader));

  // Determine the subclass of the AliFemtoManager's event reader. If none of
  // the dynamic casts resove the event reader, we do NOT add it to the task
  // and optionally print a warning. The analysis will presumably end.
  //
  // I do not know if this is the expected behavior, as the someone could write
  // their own subclass of AliFemtoEventReader, and not be able to use it in an
  // analysis as there is no other apparent way to set fReader.
  if (eventReader == NULL) {
    if (fVerbose)
      AliWarning("No AliFemto event reader created. Will not run femto analysis.\n");
  }
  else if (dynamic_cast<AliFemtoEventReaderESDChain*>(eventReader) != NULL) {
    SetFemtoReaderESD((AliFemtoEventReaderESDChain*) eventReader);
  }
  else if (dynamic_cast<AliFemtoEventReaderESDChainKine*>(eventReader) != NULL) {
    SetFemtoReaderESDKine((AliFemtoEventReaderESDChainKine *) eventReader);
  }
  else if (dynamic_cast<AliFemtoEventReaderAODChain*>(eventReader) != NULL) {
    SetFemtoReaderAOD((AliFemtoEventReaderAODChain *) eventReader);
  }
  else if (dynamic_cast<AliFemtoEventReaderStandard*>(eventReader) != NULL) {
    SetFemtoReaderStandard((AliFemtoEventReaderStandard *) eventReader);
  }
  else if (dynamic_cast<AliFemtoEventReaderKinematicsChain*>(eventReader) != NULL) {
    SetFemtoReaderKinematics((AliFemtoEventReaderKinematicsChain *) eventReader);
  }
  else if (dynamic_cast<AliFemtoEventReaderKinematicsChainESD*>(eventReader) != NULL) {
    SetFemtoReaderKinematicsESD((AliFemtoEventReaderKinematicsChainESD *) eventReader);
  }
  else {
    if (fVerbose)
      AliWarning("Specified AliFemto event reader does *not* inherit from an "
                 "approved AliFemtoEventReader subclass. Will not run femto analysis.\n");\
  }
}
