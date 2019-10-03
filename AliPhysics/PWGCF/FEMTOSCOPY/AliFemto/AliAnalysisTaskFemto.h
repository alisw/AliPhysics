///
/// \file AliAnalysisTaskFemto.h
///
/// \class AliAnalysisTaskFemto
/// \brief An AliAnalysisTaskSE class for handling and running AliFemtoAnalysis objects
///
/// Creates the necessary connection between the ESD or AOD input and the
/// femtoscopic code.
///
/// As with all [AliAnalysisTasks](@ref AliAnalysisTask), a AliAnalysisTaskFemto
/// should be created from your macro and added to the AliAnalysisManager. The
/// proper way to construct an AliAnalysisTaskFemto is through the use of a
/// configuration macro file which defines a function with signature
///    `AliFemtoManager* ConfigFemtoAnalysis()`
/// whose name is specified to the constructor. The default name of this file
/// is 'ConfigFemtoAnalysis.C'. The macro is responsible for creating an
/// AliFemtoAnalysisManager and loading your custom AliFemtoAnalyses. Upon the
/// call to ::CreateOutputObjects, the macro will be loaded, and the analysis
/// manager created and added to the AliAnalysisTaskFemto. Parameters can be
/// specified in the constructor via the TString `aConfigParams`.
///
/// \author Adam Kisiel <Adam.Kisiel@cern.ch>, OSU
///

#pragma once

#ifndef ALIANALYSISTASKFEMTO_H
#define ALIANALYSISTASKFEMTO_H

#include "TH1.h"
#include "TString.h"
#include "TMacro.h"

#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"

//#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODChain.h"
#include "AliFemtoEventReaderNanoAODChain.h"
#include "AliFemtoEventReaderStandard.h"
#include "AliFemtoEventReaderKinematicsChain.h"
#include "AliFemtoEventReaderKinematicsChainESD.h"
#include "AliFemtoEventReaderAODKinematicsChain.h"
#include "AliFemtoManager.h"

#include "AliESDpid.h"
#include "AliAODpidUtil.h"
#include "AliAODHeader.h"


class AliAnalysisTaskFemto : public AliAnalysisTaskSE { //AliAnalysisTask
public:

  /// Default parameters initialize all values to zero (NULL) and empty strings
  ///
  /// Do not use this constructor - objects created with it have no way to set
  /// the configuration macro or parameters. Use an alternative constructor.
  AliAnalysisTaskFemto();

  /// Full Constructor - Set the name of the task, configuration macro filename
  /// and paramters, and optional verbosity flag.
  AliAnalysisTaskFemto(TString name, TString aConfigMacro, TString aConfigParams, Bool_t aVerbose=kFALSE, Bool_t aGridConfig=kFALSE, TString aUserName = "", TString aConfigUserName ="ConfigFemtoAnalysis");

  /// Construct with task name, configuration filename, and verbosity flag.
  ///
  /// The paramters are set to the empty string.
  AliAnalysisTaskFemto(TString name, TString aConfigMacro="ConfigFemtoAnalysis.C", Bool_t aVerbose=kFALSE, Bool_t aGridConfig=kFALSE, TString aUserName = "", TString aConfigUserName ="ConfigFemtoAnalysis");

  /// Copy Constructor - should not be used
  AliAnalysisTaskFemto(const AliAnalysisTaskFemto& aFemtoTask);
  virtual ~AliAnalysisTaskFemto();

  AliAnalysisTaskFemto& operator=(const AliAnalysisTaskFemto& aFemtoTask);

  /// ConnectInputData
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *);
  virtual void FinishTaskOutput();

  /// Set the femtomanager containing this task's analyses.
  void SetFemtoManager(AliFemtoManager *aManager);

  void SetFemtoReaderESD(AliFemtoEventReaderESDChain *aReader);
  void SetFemtoReaderESDKine(AliFemtoEventReaderESDChainKine *aReader);
  void SetFemtoReaderAOD(AliFemtoEventReaderAODChain *aReader);
  void SetFemtoReaderNanoAOD(AliFemtoEventReaderNanoAODChain *aReader);
  void SetFemtoReaderStandard(AliFemtoEventReaderStandard *aReader);
  void SetFemtoReaderKinematics(AliFemtoEventReaderKinematicsChain *aReader);
  void SetFemtoReaderKinematicsESD(AliFemtoEventReaderKinematicsChainESD *aReader);
  void SetFemtoReaderAODKinematics(AliFemtoEventReaderAODKinematicsChain *aReader);

  void Set1DCorrectionsPions(TH1D *h1);
  void Set1DCorrectionsKaons(TH1D *h1);
  void Set1DCorrectionsProtons(TH1D *h1);
  void Set1DCorrectionsPionsMinus(TH1D *h1);
  void Set1DCorrectionsKaonsMinus(TH1D *h1);
  void Set1DCorrectionsProtonsMinus(TH1D *h1);
  void Set1DCorrectionsAll(TH1D *h1);
  void Set1DCorrectionsLambdas(TH1D *h1);
  void Set1DCorrectionsLambdasMinus(TH1D *h1);
  void Set1DCorrectionsXiMinus(TH1D *h1);
  void Set1DCorrectionsXiPlus(TH1D *h1);

  void Set4DCorrectionsPions(THnSparse *h1);
  void Set4DCorrectionsKaons(THnSparse *h1);
  void Set4DCorrectionsProtons(THnSparse *h1);
  void Set4DCorrectionsPionsMinus(THnSparse *h1);
  void Set4DCorrectionsKaonsMinus(THnSparse *h1);
  void Set4DCorrectionsProtonsMinus(THnSparse *h1);
  void Set4DCorrectionsAll(THnSparse *h1);
  void Set4DCorrectionsLambdas(THnSparse *h1);
  void Set4DCorrectionsLambdasMinus(THnSparse *h1);
  void LoadMacro(TMacro *macro);
  void SaveConfigTMacro(Bool_t save);
  void SetGRIDUserName(TString aUserName);

protected:
  AliESDEvent          *fESD;          //!<! ESD object
  AliESDpid            *fESDpid;       //!<! ESDpid object
  AliVEvent            *fVEvent;       //!<! AliVEvent object
  AliAODEvent          *fAOD;          //!<! AOD object
  AliAODpidUtil        *fAODpidUtil;   ///<  AliAODpidUtil object
  AliAODHeader         *fAODheader;     ///< AliAODHeader object (to get reference multiplicity in pp)
  AliNanoAODHeader     *fNanoAODheader; //!<! AliNanoAODHeader object

  AliStack             *fStack;        //!<! Stack from Kinematics
  TList                *fOutputList;   ///<  AliFemto results list
  AliFemtoEventReader  *fReader;       //!<! Reference to the reader
  AliFemtoManager      *fManager;      //!<! AliFemto top-level manager
  Int_t                fAnalysisType;  ///<  Mark ESD of AOD analysis
  TString              fConfigMacro;   ///<  Config macro location
  TString              fConfigParams;  ///<  Config macro parameters
  Bool_t               fVerbose;

  TH1D                 *f1DcorrectionsPions; //file with corrections, pT dependant
  TH1D                 *f1DcorrectionsKaons; //file with corrections, pT dependant
  TH1D                 *f1DcorrectionsProtons; //file with corrections, pT dependant
  TH1D                 *f1DcorrectionsPionsMinus; //file with corrections, pT dependant
  TH1D                 *f1DcorrectionsKaonsMinus; //file with corrections, pT dependant
  TH1D                 *f1DcorrectionsProtonsMinus; //file with corrections, pT dependant
  TH1D                 *f1DcorrectionsAll; //file with corrections, pT dependant
  TH1D                 *f1DcorrectionsLambdas; //file with corrections, pT dependant
  TH1D                 *f1DcorrectionsLambdasMinus; //file with corrections, pT dependant
  TH1D                 *f1DcorrectionsXiMinus; //file with corrections, pT dependant
  TH1D                 *f1DcorrectionsXiPlus; //file with corrections, pT dependant

  THnSparse            *f4DcorrectionsPions; //file with corrections, pT dependant
  THnSparse            *f4DcorrectionsKaons; //file with corrections, pT dependant
  THnSparse            *f4DcorrectionsProtons; //file with corrections, pT dependant
  THnSparse            *f4DcorrectionsPionsMinus; //file with corrections, pT dependant
  THnSparse            *f4DcorrectionsKaonsMinus; //file with corrections, pT dependant
  THnSparse            *f4DcorrectionsProtonsMinus; //file with corrections, pT dependant
  THnSparse            *f4DcorrectionsAll; //file with corrections, pT dependant
  THnSparse            *f4DcorrectionsLambdas; //file with corrections, pT dependant
  THnSparse            *f4DcorrectionsLambdasMinus; //file with corrections, pT dependant

  Bool_t fGridConfig; //use config stored on the grid
  TMacro *fConfigTMacro; //macro of the config file
  Bool_t fSaveConfigTMacro; //flag to save config TMacro in output list
  TString fUserName; //GRID user name
  TString fconfigFunName; //name of the config fucntion (like "ConfigFemtoAnalysis")


  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskFemto, 3);
  /// \endcond
};

inline
AliAnalysisTaskFemto::AliAnalysisTaskFemto():
  AliAnalysisTaskSE(),
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
  fConfigMacro(),
  fConfigParams(),
  fVerbose(kTRUE),
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
  fGridConfig(false),
  fConfigTMacro(NULL),
  fSaveConfigTMacro(false),
  fUserName(),
  fconfigFunName()
{
  /* no-op */
}

#endif
