///
/// \file AliAnalysisTaskFemto.cxx
///
#include <sstream>

#include "TROOT.h"
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TInterpreter.h"
#include "TMacro.h"
#include "TFile.h"
#include "TGrid.h"

//#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"

#include "AliESDEvent.h"

#include "AliFemtoAnalysis.h"
#include "AliAnalysisTaskFemto.h"
#include "AliVHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliFemtoEventReaderNanoAODChain.h"

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliAnalysisTaskFemto);
  /// \endcond
#endif

// Default name for the setup macro of femto analysis
// This function MUST be defined in the separate file !!!
// extern AliFemtoManager *ConfigFemtoAnalysis();

//________________________________________________________________________
AliAnalysisTaskFemto::AliAnalysisTaskFemto(TString name,
                                           TString aConfigMacro,
                                           TString aConfigParams,
                                           Bool_t aVerbose,
					   Bool_t aGridConfig,
					   TString aUserName,
					   TString aConfigFunName):
  AliAnalysisTaskSE(name), //AliAnalysisTask(name,""),
  fESD(NULL),
  fESDpid(NULL),
  fVEvent(NULL),
  fAOD(NULL),
  fAODpidUtil(NULL),
  fAODheader(NULL),
  fNanoAODheader(NULL),
  fStack(NULL),
  fOutputList(NULL),
  fReader(NULL),
  fManager(NULL),
  fAnalysisType(0),
  fConfigMacro(aConfigMacro),
  fConfigParams(aConfigParams),
  fVerbose(aVerbose),
  f1DcorrectionsPions(NULL),
  f1DcorrectionsKaons(NULL),
  f1DcorrectionsProtons(NULL),
  f1DcorrectionsPionsMinus(NULL),
  f1DcorrectionsKaonsMinus(NULL),
  f1DcorrectionsProtonsMinus(NULL),
  f1DcorrectionsAll(NULL),
  f1DcorrectionsLambdas(NULL),
  f1DcorrectionsLambdasMinus(NULL),
  f1DcorrectionsXiMinus(NULL),
  f1DcorrectionsXiPlus(NULL),
  f4DcorrectionsPions(NULL),
  f4DcorrectionsKaons(NULL),
  f4DcorrectionsProtons(NULL),
  f4DcorrectionsPionsMinus(NULL),
  f4DcorrectionsKaonsMinus(NULL),
  f4DcorrectionsProtonsMinus(NULL),
  f4DcorrectionsAll(NULL),
  f4DcorrectionsLambdas(NULL),
  f4DcorrectionsLambdasMinus(NULL),
  fGridConfig(aGridConfig),
  fConfigTMacro(NULL),
  fSaveConfigTMacro(NULL),
  fUserName(aUserName),
  fconfigFunName(aConfigFunName)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  //DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TList::Class());

}
//________________________________________________________________________
AliAnalysisTaskFemto::AliAnalysisTaskFemto(TString name,
                                           TString aConfigMacro,
                                           Bool_t aVerbose,
					   Bool_t aGridConfig,
					   TString aUserName,
					   TString aConfigFunName):
  AliAnalysisTaskSE(name), //AliAnalysisTask(name,""),
  fESD(NULL),
  fESDpid(NULL),
  fVEvent(NULL),
  fAOD(NULL),
  fAODpidUtil(NULL),
  fAODheader(NULL),
  fNanoAODheader(NULL),
  fStack(NULL),
  fOutputList(NULL),
  fReader(NULL),
  fManager(NULL),
  fAnalysisType(0),
  fConfigMacro(aConfigMacro),
  fConfigParams(""),
  fVerbose(aVerbose),
  f1DcorrectionsPions(NULL),
  f1DcorrectionsKaons(NULL),
  f1DcorrectionsProtons(NULL),
  f1DcorrectionsPionsMinus(NULL),
  f1DcorrectionsKaonsMinus(NULL),
  f1DcorrectionsProtonsMinus(NULL),
  f1DcorrectionsAll(NULL),
  f1DcorrectionsLambdas(NULL),
  f1DcorrectionsLambdasMinus(NULL),
  f1DcorrectionsXiMinus(NULL),
  f1DcorrectionsXiPlus(NULL),
  f4DcorrectionsPions(NULL),
  f4DcorrectionsKaons(NULL),
  f4DcorrectionsProtons(NULL),
  f4DcorrectionsPionsMinus(NULL),
  f4DcorrectionsKaonsMinus(NULL),
  f4DcorrectionsProtonsMinus(NULL),
  f4DcorrectionsAll(NULL),
  f4DcorrectionsLambdas(NULL),
  f4DcorrectionsLambdasMinus(NULL),
  fGridConfig(aGridConfig),
  fConfigTMacro(NULL),
  fSaveConfigTMacro(false),
  fUserName(aUserName),
  fconfigFunName(aConfigFunName)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  //DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TList::Class());

}

AliAnalysisTaskFemto::AliAnalysisTaskFemto(const AliAnalysisTaskFemto &aFemtoTask):
  AliAnalysisTaskSE(aFemtoTask), //AliAnalysisTask(aFemtoTask),
  fESD(aFemtoTask.fESD),
  fESDpid(aFemtoTask.fESDpid),
  fVEvent(aFemtoTask.fVEvent),
  fAOD(aFemtoTask.fAOD),
  fAODpidUtil(aFemtoTask.fAODpidUtil),
  fAODheader(aFemtoTask.fAODheader),
  fNanoAODheader(aFemtoTask.fNanoAODheader),
  fStack(aFemtoTask.fStack),
  fOutputList(aFemtoTask.fOutputList),
  fReader(aFemtoTask.fReader),
  fManager(aFemtoTask.fManager),
  fAnalysisType(aFemtoTask.fAnalysisType),
  fConfigMacro(aFemtoTask.fConfigMacro),
  fConfigParams(aFemtoTask.fConfigParams),
  fVerbose(aFemtoTask.fVerbose),
  f1DcorrectionsPions(aFemtoTask.f1DcorrectionsPions),
  f1DcorrectionsKaons(aFemtoTask.f1DcorrectionsKaons),
  f1DcorrectionsProtons(aFemtoTask.f1DcorrectionsProtons),
  f1DcorrectionsPionsMinus(aFemtoTask.f1DcorrectionsPionsMinus),
  f1DcorrectionsKaonsMinus(aFemtoTask.f1DcorrectionsKaonsMinus),
  f1DcorrectionsProtonsMinus(aFemtoTask.f1DcorrectionsProtonsMinus),
  f1DcorrectionsAll(aFemtoTask.f1DcorrectionsAll),
  f1DcorrectionsLambdas(aFemtoTask.f1DcorrectionsLambdas),
  f1DcorrectionsLambdasMinus(aFemtoTask.f1DcorrectionsLambdasMinus),
  f1DcorrectionsXiMinus(aFemtoTask.f1DcorrectionsXiMinus),
  f1DcorrectionsXiPlus(aFemtoTask.f1DcorrectionsXiPlus),
  f4DcorrectionsPions(aFemtoTask.f4DcorrectionsPions),
  f4DcorrectionsKaons(aFemtoTask.f4DcorrectionsKaons),
  f4DcorrectionsProtons(aFemtoTask.f4DcorrectionsProtons),
  f4DcorrectionsPionsMinus(aFemtoTask.f4DcorrectionsPionsMinus),
  f4DcorrectionsKaonsMinus(aFemtoTask.f4DcorrectionsKaonsMinus),
  f4DcorrectionsProtonsMinus(aFemtoTask.f4DcorrectionsProtonsMinus),
  f4DcorrectionsAll(aFemtoTask.f4DcorrectionsAll),
  f4DcorrectionsLambdas(aFemtoTask.f4DcorrectionsLambdas),
  f4DcorrectionsLambdasMinus(aFemtoTask.f4DcorrectionsLambdasMinus),
  fGridConfig(aFemtoTask.fGridConfig),
  fConfigTMacro(aFemtoTask.fConfigTMacro),
  fSaveConfigTMacro(aFemtoTask.fSaveConfigTMacro),
  fUserName(aFemtoTask.fUserName),
  fconfigFunName(aFemtoTask.fconfigFunName)
{
  // copy constructor
}


AliAnalysisTaskFemto &AliAnalysisTaskFemto::operator=(const AliAnalysisTaskFemto &aFemtoTask)
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

  f1DcorrectionsPions = aFemtoTask.f1DcorrectionsPions;
  f1DcorrectionsKaons = aFemtoTask.f1DcorrectionsKaons;
  f1DcorrectionsProtons = aFemtoTask.f1DcorrectionsProtons;
  f1DcorrectionsPionsMinus = aFemtoTask.f1DcorrectionsPionsMinus;
  f1DcorrectionsKaonsMinus = aFemtoTask.f1DcorrectionsKaonsMinus;
  f1DcorrectionsProtonsMinus = aFemtoTask.f1DcorrectionsProtonsMinus;
  f1DcorrectionsAll = aFemtoTask.f1DcorrectionsAll;
  f1DcorrectionsLambdas = aFemtoTask.f1DcorrectionsLambdas;
  f1DcorrectionsLambdasMinus = aFemtoTask.f1DcorrectionsLambdasMinus;
  f1DcorrectionsXiMinus = aFemtoTask.f1DcorrectionsXiMinus;
  f1DcorrectionsXiPlus = aFemtoTask.f1DcorrectionsXiPlus;


  f4DcorrectionsPions = aFemtoTask.f4DcorrectionsPions;
  f4DcorrectionsKaons = aFemtoTask.f4DcorrectionsKaons;
  f4DcorrectionsProtons = aFemtoTask.f4DcorrectionsProtons;
  f4DcorrectionsPionsMinus = aFemtoTask.f4DcorrectionsPionsMinus;
  f4DcorrectionsKaonsMinus = aFemtoTask.f4DcorrectionsKaonsMinus;
  f4DcorrectionsProtonsMinus = aFemtoTask.f4DcorrectionsProtonsMinus;
  f4DcorrectionsAll = aFemtoTask.f4DcorrectionsAll;
  f4DcorrectionsLambdas = aFemtoTask.f4DcorrectionsLambdas;
  f4DcorrectionsLambdasMinus = aFemtoTask.f4DcorrectionsLambdasMinus;

  fGridConfig = aFemtoTask.fGridConfig;
  fConfigTMacro = aFemtoTask.fConfigTMacro;
  fSaveConfigTMacro = aFemtoTask.fSaveConfigTMacro;
  fUserName = aFemtoTask.fUserName;
  fconfigFunName = aFemtoTask.fconfigFunName;

  return *this;
}

AliAnalysisTaskFemto::~AliAnalysisTaskFemto()
{
}


//________________________________________________________________________
void AliAnalysisTaskFemto::ConnectInputData(Option_t *)
{
  AliInfo(Form("   ConnectInputData %s\n", GetName()));
  fESD = nullptr;
  fESDpid = nullptr;
  fAOD = nullptr;
  fAODpidUtil = nullptr;
  fAODheader = nullptr;
  fNanoAODheader = nullptr;
  fAnalysisType = 0;

  TTree *tree = dynamic_cast<TTree *>(GetInputData(0));
  if (!tree) {
    AliWarning("Could not read chain from input slot 0");
    return;
  }

  if (auto *femtoReader = dynamic_cast<AliFemtoEventReaderESDChain *>(fReader)) {
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
  else if (auto *femtoReaderESDKine = dynamic_cast<AliFemtoEventReaderESDChainKine *>(fReader)) {
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

  if (auto *femtoReaderAOD = dynamic_cast<AliFemtoEventReaderAODChain *>(fReader)) {
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
        cout << "AliAnalysisTaskFemto::AodpidUtil:" << fAODpidUtil << endl;
      femtoReaderAOD->SetAODpidUtil(fAODpidUtil);

      //Applying 1D corrections

      if(f1DcorrectionsPions) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections pions"<<f1DcorrectionsPions;
        femtoReaderAOD->Set1DCorrectionsPions(f1DcorrectionsPions);
      }
      if(f1DcorrectionsKaons) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections kaons"<<f1DcorrectionsKaons;
        femtoReaderAOD->Set1DCorrectionsKaons(f1DcorrectionsKaons);
      }
      if(f1DcorrectionsProtons) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections protons"<<f1DcorrectionsProtons;
        femtoReaderAOD->Set1DCorrectionsProtons(f1DcorrectionsProtons);
      }
      if(f1DcorrectionsPionsMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections pions Minus"<<f1DcorrectionsPionsMinus;
        femtoReaderAOD->Set1DCorrectionsPionsMinus(f1DcorrectionsPionsMinus);
      }
      if(f1DcorrectionsKaonsMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections kaons Minus"<<f1DcorrectionsKaonsMinus;
        femtoReaderAOD->Set1DCorrectionsKaonsMinus(f1DcorrectionsKaonsMinus);
      }
      if(f1DcorrectionsProtonsMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections protons Minus"<<f1DcorrectionsProtonsMinus;
        femtoReaderAOD->Set1DCorrectionsProtonsMinus(f1DcorrectionsProtonsMinus);
      }
      if(f1DcorrectionsAll) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections all"<<f1DcorrectionsAll;
        femtoReaderAOD->Set1DCorrectionsAll(f1DcorrectionsAll);
      }
      if(f1DcorrectionsLambdas) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections lambas"<<f1DcorrectionsLambdas;
        femtoReaderAOD->Set1DCorrectionsLambdas(f1DcorrectionsLambdas);
      }
      if(f1DcorrectionsLambdasMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections lambas Minus"<<f1DcorrectionsLambdasMinus;
        femtoReaderAOD->Set1DCorrectionsLambdasMinus(f1DcorrectionsLambdasMinus);
      }

      if(f1DcorrectionsXiMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections xi Minus"<< f1DcorrectionsXiMinus;
        femtoReaderAOD->Set1DCorrectionsXiMinus(f1DcorrectionsXiMinus);
      }
      if(f1DcorrectionsXiPlus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections xi plus"<<f1DcorrectionsXiPlus;
        femtoReaderAOD->Set1DCorrectionsXiPlus(f1DcorrectionsXiPlus);
      }


      //Applying 4D corrections
      if(f4DcorrectionsPions) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections pions"<<f4DcorrectionsPions;
        femtoReaderAOD->Set4DCorrectionsPions(f4DcorrectionsPions);
      }
      if(f4DcorrectionsKaons) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections kaons"<<f4DcorrectionsKaons;
        femtoReaderAOD->Set4DCorrectionsKaons(f4DcorrectionsKaons);
      }
      if(f4DcorrectionsProtons) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections protons"<<f4DcorrectionsProtons;
        femtoReaderAOD->Set4DCorrectionsProtons(f4DcorrectionsProtons);
      }
      if(f4DcorrectionsPionsMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections pions Minus"<<f4DcorrectionsPionsMinus;
        femtoReaderAOD->Set4DCorrectionsPionsMinus(f4DcorrectionsPionsMinus);
      }
      if(f4DcorrectionsKaonsMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections kaons Minus"<<f4DcorrectionsKaonsMinus;
        femtoReaderAOD->Set4DCorrectionsKaonsMinus(f4DcorrectionsKaonsMinus);
      }
      if(f4DcorrectionsProtonsMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections protons Minus"<<f4DcorrectionsProtonsMinus;
        femtoReaderAOD->Set4DCorrectionsProtonsMinus(f4DcorrectionsProtonsMinus);
      }
      if(f4DcorrectionsAll) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections all"<<f4DcorrectionsAll;
        femtoReaderAOD->Set4DCorrectionsAll(f4DcorrectionsAll);
      }
      if(f4DcorrectionsLambdas) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections lambas"<<f4DcorrectionsLambdas;
        femtoReaderAOD->Set4DCorrectionsLambdas(f4DcorrectionsLambdas);
      }
      if(f4DcorrectionsLambdasMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections lambas Minus"<<f4DcorrectionsLambdasMinus;
        femtoReaderAOD->Set4DCorrectionsLambdasMinus(f4DcorrectionsLambdasMinus);
      }

      fAODheader = dynamic_cast<AliAODHeader *>(fAOD->GetHeader());
      if (!fAODheader) AliFatal("Not a standard AOD");
      femtoReaderAOD->SetAODheader(fAODheader);

    }
  }

  if (auto *femtoReaderNanoAOD = dynamic_cast<AliFemtoEventReaderNanoAODChain *>(fReader)) {
    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    AliAnalysisTaskSE::ConnectInputData();
    if (!aodH) {
      TObject *handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
      if (fVerbose) {
        AliInfo("Output handler ");
      }
      if (handler && handler->InheritsFrom("AliAODHandler")) {
        if (fVerbose)
          AliInfo("Selected NanoAOD analysis");
      }
      else {
        if (fVerbose)
          AliWarning("Selected NanoAOD reader but no AOD handler found");
      }
    }
    else{
      //fAOD = ((AliAODHandler *)handler)->GetAOD();

      fAnalysisType = 3;

      //Applying 1D corrections
      if (fVerbose)
	cout << "AliAnalysisTaskFemto::Applying 1D corrections:" << f1DcorrectionsPions << endl;

      if(f1DcorrectionsPions) {
	if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections pions"<<f1DcorrectionsPions;
	femtoReaderNanoAOD->Set1DCorrectionsPions(f1DcorrectionsPions);
      }
      if(f1DcorrectionsKaons) {
	if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections kaons"<<f1DcorrectionsKaons;
	femtoReaderNanoAOD->Set1DCorrectionsKaons(f1DcorrectionsKaons);
      }
      if(f1DcorrectionsProtons) {
	if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections protons"<<f1DcorrectionsProtons;
	femtoReaderNanoAOD->Set1DCorrectionsProtons(f1DcorrectionsProtons);
      }
      if(f1DcorrectionsPionsMinus) {
	if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections pions Minus"<<f1DcorrectionsPionsMinus;
	femtoReaderNanoAOD->Set1DCorrectionsPionsMinus(f1DcorrectionsPionsMinus);
      }
      if(f1DcorrectionsKaonsMinus) {
	if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections kaons Minus"<<f1DcorrectionsKaonsMinus;
	femtoReaderNanoAOD->Set1DCorrectionsKaonsMinus(f1DcorrectionsKaonsMinus);
      }
      if(f1DcorrectionsProtonsMinus) {
	if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections protons Minus"<<f1DcorrectionsProtonsMinus;
	femtoReaderNanoAOD->Set1DCorrectionsProtonsMinus(f1DcorrectionsProtonsMinus);
      }
      if(f1DcorrectionsAll) {
	if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections all"<<f1DcorrectionsAll;
	femtoReaderNanoAOD->Set1DCorrectionsAll(f1DcorrectionsAll);
      }
      if(f1DcorrectionsLambdas) {
	if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections lambas"<<f1DcorrectionsLambdas;
	femtoReaderNanoAOD->Set1DCorrectionsLambdas(f1DcorrectionsLambdas);
      }
      if(f1DcorrectionsLambdasMinus) {
	if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections lambas Minus"<<f1DcorrectionsLambdasMinus;
	femtoReaderNanoAOD->Set1DCorrectionsLambdasMinus(f1DcorrectionsLambdasMinus);
      }

  if(f1DcorrectionsXiMinus) {
     if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections xi Minus"<< f1DcorrectionsXiMinus;
  femtoReaderNanoAOD->Set1DCorrectionsXiMinus(f1DcorrectionsXiMinus);
      }
  if(f1DcorrectionsXiPlus) {
     if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 1d corrections xi plus"<<f1DcorrectionsXiPlus;
  femtoReaderNanoAOD->Set1DCorrectionsXiPlus(f1DcorrectionsXiPlus);
      }
//Applying 4D corrections
      if(f4DcorrectionsPions) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections pions"<<f4DcorrectionsPions;
        femtoReaderNanoAOD->Set4DCorrectionsPions(f4DcorrectionsPions);
      }
      if(f4DcorrectionsKaons) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections kaons"<<f4DcorrectionsKaons;
        femtoReaderNanoAOD->Set4DCorrectionsKaons(f4DcorrectionsKaons);
      }
      if(f4DcorrectionsProtons) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections protons"<<f4DcorrectionsProtons;
        femtoReaderNanoAOD->Set4DCorrectionsProtons(f4DcorrectionsProtons);
      }
      if(f4DcorrectionsPionsMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections pions Minus"<<f4DcorrectionsPionsMinus;
        femtoReaderNanoAOD->Set4DCorrectionsPionsMinus(f4DcorrectionsPionsMinus);
      }
      if(f4DcorrectionsKaonsMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections kaons Minus"<<f4DcorrectionsKaonsMinus;
        femtoReaderNanoAOD->Set4DCorrectionsKaonsMinus(f4DcorrectionsKaonsMinus);
      }
      if(f4DcorrectionsProtonsMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections protons Minus"<<f4DcorrectionsProtonsMinus;
        femtoReaderNanoAOD->Set4DCorrectionsProtonsMinus(f4DcorrectionsProtonsMinus);
      }
      if(f4DcorrectionsAll) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections all"<<f4DcorrectionsAll;
        femtoReaderNanoAOD->Set4DCorrectionsAll(f4DcorrectionsAll);
      }
      if(f4DcorrectionsLambdas) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections lambas"<<f4DcorrectionsLambdas;
        femtoReaderNanoAOD->Set4DCorrectionsLambdas(f4DcorrectionsLambdas);
      }
      if(f4DcorrectionsLambdasMinus) {
        if (fVerbose)	cout<<"AliAnalysisTaskFemto::Setting 4d corrections lambas Minus"<<f4DcorrectionsLambdasMinus;
        femtoReaderNanoAOD->Set4DCorrectionsLambdasMinus(f4DcorrectionsLambdasMinus);
      }

      // cout<<"nano AOD header: "<<fVEvent<<" "<<fVEvent->GetHeader()<<endl;
      // fNanoAODheader = dynamic_cast<AliNanoAODHeader *>(fVEvent->GetHeader());
      // if (!fNanoAODheader) AliFatal("Not a standard NanoAOD");
      // femtoReaderNanoAOD->SetAODheader(fNanoAODheader);

    }

  }

  if (auto *femtoReaderAODKine = dynamic_cast<AliFemtoEventReaderAODKinematicsChain *>(fReader)) {
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

      fAODheader = dynamic_cast<AliAODHeader *>(fAOD->GetHeader());
      if (!fAODheader) AliFatal("Not a standard AOD");
      femtoReaderAODKine->SetAODheader(fAODheader);

    }
  }

  if ((!fAOD) && (!fESD)) {
    if (fVerbose)
      AliWarning("Wrong analysis type: Only ESD and AOD types are allowed!");
  }
}

//________________________________________________________________________
void AliAnalysisTaskFemto::CreateOutputObjects()
{
  if (fVerbose)
    AliInfo("Creating Femto Analysis objects\n");

  gSystem->SetIncludePath("-I$ROOTSYS/include -I./STEERBase/ -I./ESD/ -I./AOD/ -I./ANALYSIS/ -I./ANALYSISalice/ -I./PWG2AOD/AOD -I./PWG2femtoscopy/FEMTOSCOPY/AliFemto -I./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser");
  if(!fGridConfig)
    {
      gROOT->LoadMacro(fConfigMacro);
      if(fSaveConfigTMacro)
	fConfigTMacro = new TMacro(fConfigMacro);
    }
  else
    {
      printf("*** Connect to AliEn ***\n");
      TGrid::Connect("alien://");
      TFile *fileConfig = TFile::Open(fConfigMacro.Data());
      fConfigTMacro = dynamic_cast<TMacro*>(fileConfig->Get(fconfigFunName.Data())->Clone());
      LoadMacro(fConfigTMacro);
      fileConfig->Close();
    }


  TString cmd = Form("%s(%s)", fconfigFunName.Data(),fConfigParams.Data());
  auto *femto_manager = reinterpret_cast<AliFemtoManager*>(gInterpreter->ProcessLine(cmd));
  if (femto_manager == nullptr) {
    AliError(Form("ConfigFemtoAnalysis function returned NULL (i.e. no manager)\n--- invoked function ---\n%s\n---", cmd.Data()));
  }

  SetFemtoManager(femto_manager);

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


  if(fSaveConfigTMacro && fConfigTMacro)
    fOutputList->Add(fConfigTMacro);

  PostData(0, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskFemto::Exec(Option_t *)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  auto *event_handler = static_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());

  // Task making a femtoscopic analysis.
  if (fOfflineTriggerMask) {
    Bool_t isSelected = event_handler->IsEventSelected() & fOfflineTriggerMask;
    if (!isSelected) {
      if (fVerbose) {
        // std::cout << "AliAnalysisTaskFemto: is not selected" << endl;
        AliInfo("Event is not selected");
      }
      return;
    }
  }

  if (fAnalysisType == 1) {
    if (!fESD) {
      if (fVerbose) {
        AliWarning("fESD not available");
      }
      return;
    }
    //Get MC data
    AliMCEventHandler *mctruth = static_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());

    AliGenHijingEventHeader *hdh = 0;
    AliGenCocktailEventHeader *hd = 0;
    AliGenEventHeader *header = 0;
    if (mctruth) {
      fStack = mctruth->MCEvent()->Stack();

      hd = dynamic_cast<AliGenCocktailEventHeader *>(mctruth->MCEvent()->GenEventHeader());
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

      if (auto *fesdc = dynamic_cast<AliFemtoEventReaderESDChain *>(fReader)) {
        // Process the event with no Kine information
        fesdc->SetESDSource(fESD);
        fManager->ProcessEvent();
      }
    }

    if (auto *fkinec = dynamic_cast<AliFemtoEventReaderKinematicsChain *>(fReader)) {
      // Process the event with Kine information only
      fkinec->SetStackSource(fStack);
      fkinec->SetGenEventHeader(header);
      fManager->ProcessEvent();
    }
    else if (auto *fesdck = dynamic_cast<AliFemtoEventReaderESDChainKine *>(fReader)) {
      // Process the event with Kine information
      fesdck->SetESDSource(fESD);
      fesdck->SetStackSource(fStack);
      cout << "Set Stack:" << fStack << endl;
      fesdck->SetGenEventHeader(hdh);
      fManager->ProcessEvent();
    }
    else if (auto *fkcesd = dynamic_cast<AliFemtoEventReaderKinematicsChainESD *>(fReader)) {
      // Process the event with Kine information
      fkcesd->SetESDSource(fESD);
      fkcesd->SetStackSource(fStack);
      fkcesd->SetGenEventHeader(hdh);
      fManager->ProcessEvent();
    }
    else if (auto *fstd = dynamic_cast<AliFemtoEventReaderStandard *>(fReader)) {
      // Process the event with Kine information
      fstd->SetESDSource(fESD);
      if (mctruth) {
        fstd->SetStackSource(fStack);
        fstd->SetGenEventHeader(hdh);
        fstd->SetInputType(AliFemtoEventReaderStandard::kESDKine);
      } else {
        fstd->SetInputType(AliFemtoEventReaderStandard::kESD);
      }
      fManager->ProcessEvent();
    }

    // Post the output histogram list
    if(fSaveConfigTMacro && fConfigTMacro)
      fOutputList->Add(fConfigTMacro);
    PostData(0, fOutputList);
  }

  if (fAnalysisType == 2) {

    if (!fAOD) {
      if (fVerbose)
	AliWarning("fAOD not available");
      return;
    }

    if (fVerbose) {
      AliInfo(Form("Tracks in AOD: %d \n", fAOD->GetNumberOfTracks()));
    }

    if (fAOD->GetNumberOfTracks() > 0) {
      if (!fReader) {
        if (fVerbose) {
          AliWarning("No AOD reader for AOD analysis! \n");
        }
      } else {
        if (auto *faodc = dynamic_cast<AliFemtoEventReaderAODChain *>(fReader)) {
          // Process the event
          faodc->SetAODSource(fAOD);
          fManager->ProcessEvent();
        }

        else if (auto *fstd = dynamic_cast<AliFemtoEventReaderStandard *>(fReader)) {
          // Process the event
          fstd->SetAODSource(fAOD);
          fstd->SetInputType(AliFemtoEventReaderStandard::kAOD);
          fManager->ProcessEvent();
        }

        else if (auto *faodkine = dynamic_cast<AliFemtoEventReaderAODKinematicsChain *>(fReader)) {
          // Process the event
          faodkine->SetAODSource(fAOD);
          fManager->ProcessEvent();
        }
      }
    }

    // Post the output histogram list
    if(fSaveConfigTMacro && fConfigTMacro)
      fOutputList->Add(fConfigTMacro);
    PostData(0, fOutputList);
  }

    if (fAnalysisType == 3) {

      if (auto *faodc = dynamic_cast<AliFemtoEventReaderNanoAODChain *>(fReader)) {
	// Process the event
	if (!fInputEvent)
	  {
	    return;
	  }

	faodc->SetInputEvent(fInputEvent);
	AliNanoAODHeader* nanoHeader = dynamic_cast<AliNanoAODHeader*>(fInputEvent->GetHeader());
	faodc->SetAODheader(nanoHeader);
        fManager->ProcessEvent();
      }
      // Post the output histogram list
      if(fSaveConfigTMacro && fConfigTMacro)
	fOutputList->Add(fConfigTMacro);
      PostData(0, fOutputList);
    }

}

//________________________________________________________________________
void AliAnalysisTaskFemto::Terminate(Option_t *)
{
  // Do the final processing
  if (fManager) {
    fManager->Finish();
  }
}
//________________________________________________________________________
void AliAnalysisTaskFemto:: FinishTaskOutput()
{
  // Do the final processing
  if (fManager) {
    fManager->Finish();
  }
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoReaderESD(AliFemtoEventReaderESDChain *aReader)
{
  if (fVerbose)
    AliInfo("Selecting Femto reader for ESD\n");
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoReaderESDKine(AliFemtoEventReaderESDChainKine *aReader)
{
  if (fVerbose)
    AliInfo("Selecting Femto reader for ESD with Kinematics information\n");
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoReaderAOD(AliFemtoEventReaderAODChain *aReader)
{
  if (fVerbose)
    AliInfo("Selecting Femto reader for AOD\n");
  fReader = aReader;
}

//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoReaderNanoAOD(AliFemtoEventReaderNanoAODChain *aReader)
{
  if (fVerbose)
    AliInfo("Selecting Femto reader for NanoAOD\n");
  fReader = aReader;
}

void AliAnalysisTaskFemto::SetFemtoReaderStandard(AliFemtoEventReaderStandard *aReader)
{
  if (fVerbose)
    AliInfo("Selecting Standard all-purpose Femto reader\n");
  fReader = aReader;
}
void AliAnalysisTaskFemto::SetFemtoReaderKinematics(AliFemtoEventReaderKinematicsChain *aReader)
{
  if (fVerbose)
    printf("Selecting Femto reader for Kinematics (Monte Carlo) information\n");
  fReader = aReader;
}
void AliAnalysisTaskFemto::SetFemtoReaderKinematicsESD(AliFemtoEventReaderKinematicsChainESD *aReader)
{
  if (fVerbose)
    printf("Selecting Femto reader for Kinematics (Monte Carlo) information + ESD\n");
  fReader = aReader;
}
void AliAnalysisTaskFemto::SetFemtoReaderAODKinematics(AliFemtoEventReaderAODKinematicsChain *aReader)
{
  if (fVerbose)
    printf("Selecting Femto reader for AOD kinematics (Monte Carlo) information (from AOD)\n");
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoManager(AliFemtoManager *aManager)
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
  else if (dynamic_cast<AliFemtoEventReaderNanoAODChain*>(eventReader) != NULL) {
    SetFemtoReaderNanoAOD((AliFemtoEventReaderNanoAODChain *) eventReader);
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
  else if (dynamic_cast<AliFemtoEventReaderAODKinematicsChain*>(eventReader) != NULL) {
    SetFemtoReaderAODKinematics((AliFemtoEventReaderAODKinematicsChain *) eventReader);
  }

  else {
    if (fVerbose)
      AliWarning("Specified AliFemto event reader does *not* inherit from an "
                 "approved AliFemtoEventReader subclass. Will not run femto analysis.\n");\
  }
}

void AliAnalysisTaskFemto::Set1DCorrectionsPions(TH1D *h1)
{
  if (fVerbose)
    printf("Reading corrections plus\n");
  f1DcorrectionsPions = h1;
}

void AliAnalysisTaskFemto::Set1DCorrectionsKaons(TH1D *h1)
{
  f1DcorrectionsKaons = h1;
}


void AliAnalysisTaskFemto::Set1DCorrectionsProtons(TH1D *h1)
{
  f1DcorrectionsProtons = h1;
}


void AliAnalysisTaskFemto::Set1DCorrectionsPionsMinus(TH1D *h1)
{
  if (fVerbose)
    printf("Reading corrections minus\n");
  f1DcorrectionsPionsMinus = h1;
}

void AliAnalysisTaskFemto::Set1DCorrectionsKaonsMinus(TH1D *h1)
{
  f1DcorrectionsKaonsMinus = h1;
}


void AliAnalysisTaskFemto::Set1DCorrectionsProtonsMinus(TH1D *h1)
{
  f1DcorrectionsProtonsMinus = h1;
}



void AliAnalysisTaskFemto::Set1DCorrectionsAll(TH1D *h1)
{
  f1DcorrectionsAll = h1;
}

void AliAnalysisTaskFemto::Set1DCorrectionsLambdas(TH1D *h1)
{
  f1DcorrectionsLambdas = h1;
}


void AliAnalysisTaskFemto::Set1DCorrectionsLambdasMinus(TH1D *h1)
{
  f1DcorrectionsLambdasMinus = h1;
}

void AliAnalysisTaskFemto::Set1DCorrectionsXiPlus(TH1D *h1)
{
  f1DcorrectionsXiPlus = h1;
}


void AliAnalysisTaskFemto::Set1DCorrectionsXiMinus(TH1D *h1)
{
  f1DcorrectionsXiMinus = h1;
}

void AliAnalysisTaskFemto::Set4DCorrectionsPions(THnSparse *h1)
{
  if (fVerbose)
    printf("Reading corrections\n");
  f4DcorrectionsPions = h1;
}

void AliAnalysisTaskFemto::Set4DCorrectionsKaons(THnSparse *h1)
{
  f4DcorrectionsKaons = h1;
}

void AliAnalysisTaskFemto::Set4DCorrectionsProtons(THnSparse *h1)
{
  f4DcorrectionsProtons = h1;
}

void AliAnalysisTaskFemto::Set4DCorrectionsPionsMinus(THnSparse *h1)
{
  f4DcorrectionsPionsMinus = h1;
}

void AliAnalysisTaskFemto::Set4DCorrectionsKaonsMinus(THnSparse *h1)
{
  f4DcorrectionsKaonsMinus = h1;
}

void AliAnalysisTaskFemto::Set4DCorrectionsProtonsMinus(THnSparse *h1)
{
  f4DcorrectionsProtonsMinus = h1;
}

void AliAnalysisTaskFemto::Set4DCorrectionsAll(THnSparse *h1)
{
  f4DcorrectionsAll = h1;
}

void AliAnalysisTaskFemto::Set4DCorrectionsLambdas(THnSparse *h1)
{
  f4DcorrectionsLambdas = h1;
}


void AliAnalysisTaskFemto::Set4DCorrectionsLambdasMinus(THnSparse *h1)
{
  f4DcorrectionsLambdasMinus = h1;
}


 ////////////////////////////////////////////////////////////////////////////////
 /// Load the macro into the interpreter.
 /// Function copied from TMacro class of ROOT 6, not present in ROOT 5.34
void AliAnalysisTaskFemto::LoadMacro(TMacro *macro)
 {
    if (macro == nullptr) {
      return;
    }

    std::stringstream ss;

    TList *fLines = macro->GetListOfLines();
    TIter next(fLines);
    TObjString *obj;
    while ((obj = (TObjString*) next())) {
       ss << obj->GetName() << std::endl;
    }
    gInterpreter->LoadText(ss.str().c_str());

    if(fSaveConfigTMacro)
      fConfigTMacro = dynamic_cast<TMacro*>(macro->Clone());
 }

void AliAnalysisTaskFemto::SaveConfigTMacro(Bool_t save)
{
  fSaveConfigTMacro = save;
}

void AliAnalysisTaskFemto::SetGRIDUserName(TString aUserName)
{
  fUserName = aUserName;
}
