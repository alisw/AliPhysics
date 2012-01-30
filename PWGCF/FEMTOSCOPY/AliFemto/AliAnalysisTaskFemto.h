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

#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"

#include "AliAnalysisTask.h"
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


class AliAnalysisTaskFemto : public AliAnalysisTask {
 public:
  AliAnalysisTaskFemto() : AliAnalysisTask(), fESD(0), fESDpid(0), fAOD(0), fAODpidUtil(0), fStack(0), fOutputList(0), fReader(0x0), fManager(0x0), fAnalysisType(0), fConfigMacro(0), fConfigParams(0) {}
  AliAnalysisTaskFemto(const char *name, const char *aConfigMacro, const char *aConfigParams);
  AliAnalysisTaskFemto(const char *name, const char *aConfigMacro);
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

  AliStack                    *fStack;        //! Stack from Kinematics
  TList                       *fOutputList;   //  AliFemto results list
  AliFemtoEventReader         *fReader;       //! Reference to the reader
  AliFemtoManager             *fManager;      //! AliFemto top-level manager 
  int                          fAnalysisType; //  Mark ESD of AOD analysis
  char                        *fConfigMacro;  //  Config macro location
  char                        *fConfigParams; //  Config macro parameters


  ClassDef(AliAnalysisTaskFemto, 3); // example of analysis
};

#endif

