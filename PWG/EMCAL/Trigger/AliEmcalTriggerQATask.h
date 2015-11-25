/**
 * \file AliEmcalTriggerQATask.h
 * \brief Class to do some fast QA of the EMCal trigger. Useful also to tune trigger thresholds.
 *
 * Class to do some fast QA of the EMCal trigger. Useful also to tune trigger thresholds.
 * The input for the process are the trigger patches AliEmcalTriggerPatchInfoAPV1 produced by the AliEmcalTriggerMaker class.
 *
 * The output is a bunch of histograms
 *
 * \author Salvatore Aiola <>, Yale University
 * \date Oct 22, 2015
 */
#ifndef ALIEMCALTRIGGERQATASK_H
#define ALIEMCALTRIGGERQATASK_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TClonesArray;
class THistManager;
class TString;
class AliEmcalTriggerQAAP;

#include "AliLog.h"
#include "AliAnalysisTaskEmcal.h"

/**
 * \class AliEmcalTriggerQATask
 * \brief EMCAL trigger QA task
 *
 * This Class does some fast QA of the EMCal trigger. Useful also to tune trigger thresholds.
 */
class AliEmcalTriggerQATask : public AliAnalysisTaskEmcal {
 public:

  AliEmcalTriggerQATask();
  AliEmcalTriggerQATask(const char *name);
  virtual ~AliEmcalTriggerQATask();

  void SetCaloTriggersInName(const char *name)      { fCaloTriggersInName      = name; }

  AliEmcalTriggerQAAP* GetTriggerQA()               { return fEMCALTriggerQA         ; }

 protected:
  void                                      UserCreateOutputObjects();
  void                                      ExecOnce();
  Bool_t                                    Run();
  Bool_t                                    FillHistograms();
  
  TString                                   fCaloTriggersInName;         ///< name of input trigger array
  AliEmcalTriggerQAAP                      *fEMCALTriggerQA;             ///< produces the QA histograms

  TClonesArray                             *fCaloTriggersIn;             //!<! trigger array out

 private:
  AliEmcalTriggerQATask(const AliEmcalTriggerQATask&);            // not implemented
  AliEmcalTriggerQATask &operator=(const AliEmcalTriggerQATask&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerQATask, 1) // Task to make QA of EMCAL trigger
  /// \endcond
};

#endif
