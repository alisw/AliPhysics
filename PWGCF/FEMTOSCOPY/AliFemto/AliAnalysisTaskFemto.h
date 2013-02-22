//------------------------------------------------------
// AliAnalysisTaskFemto - A task for the analysis framework
// from the FEMTOSCOPY analysis of PWG2. Creates the necessary
// connection between the ESD or AOD input and the femtoscopic
// code.
// Author: Adam Kisiel, OSU; Adam.Kisiel@cern.ch
//------------------------------------------------------
#ifndef ALIANALYSISTASKFEMTO_H
#define ALIANALYSISTASKFEMTO_H

#include "TH1.h"
#include "TString.h"

#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"

//#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODChain.h"
#include "AliFemtoEventReaderStandard.h"
#include "AliFemtoEventReaderKinematicsChain.h"
#include "AliFemtoManager.h"

#include "AliESDpid.h"
#include "AliAODpidUtil.h"
#include "AliAODHeader.h"


class AliAnalysisTaskFemto : public AliAnalysisTaskSE { //AliAnalysisTask
 public:
AliAnalysisTaskFemto() : AliAnalysisTaskSE(), fESD(0), fESDpid(0), fAOD(0), fAODpidUtil(0), fAODheader(0), fStack(0), fOutputList(0), fReader(0x0), fManager(0x0), fAnalysisType(0), fConfigMacro(0), fConfigParams(0), fVerbose(kTRUE) {}
    AliAnalysisTaskFemto(TString name, TString aConfigMacro, TString aConfigParams, Bool_t aVerbose=kFALSE);
    AliAnalysisTaskFemto(TString name, TString aConfigMacro, Bool_t aVerbose=kFALSE);
  AliAnalysisTaskFemto(const AliAnalysisTaskFemto& aFemtoTask);
  virtual ~AliAnalysisTaskFemto();
  
  AliAnalysisTaskFemto& operator=(const AliAnalysisTaskFemto& aFemtoTask);

  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   FinishTaskOutput();

  void SetFemtoReaderESD(AliFemtoEventReaderESDChain *aReader);
  void SetFemtoReaderESDKine(AliFemtoEventReaderESDChainKine *aReader);
  void SetFemtoReaderAOD(AliFemtoEventReaderAODChain *aReader);
  void SetFemtoReaderStandard(AliFemtoEventReaderStandard *aReader);
  void SetFemtoManager(AliFemtoManager *aManager);
  void SetFemtoReaderKinematics(AliFemtoEventReaderKinematicsChain *aReader);

 private:
  AliESDEvent                 *fESD;          //! ESD object
  AliESDpid                   *fESDpid;       //! ESDpid object
  AliAODEvent                 *fAOD;          //! AOD object
  AliAODpidUtil               *fAODpidUtil;   // AliAODpidUtil object
  AliAODHeader                *fAODheader;    //AliAODHeader object (to get reference multiplicity in pp)

  AliStack                    *fStack;        //! Stack from Kinematics
  TList                       *fOutputList;   //  AliFemto results list
  AliFemtoEventReader         *fReader;       //! Reference to the reader
  AliFemtoManager             *fManager;      //! AliFemto top-level manager 
  int                          fAnalysisType; //  Mark ESD of AOD analysis
    TString                      fConfigMacro;  //  Config macro location
    TString                      fConfigParams; //  Config macro parameters
  Bool_t fVerbose;

  ClassDef(AliAnalysisTaskFemto, 3); // example of analysis
};

#endif

