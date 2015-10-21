#ifndef ALIEMCALTRIGGERMAKERTASK_H
#define ALIEMCALTRIGGERMAKERTASK_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcal.h"
#include <TString.h>

class TClonesArray;
class THistManager;
class AliEmcalTriggerBitConfig;
class AliEmcalTriggerMakerKernel;
class AliVVZERO;

class AliEmcalTriggerMakerTask : public AliAnalysisTaskEmcal {
public:
  /***
   * \enum TriggerMakerTriggerBitConfig_t
   * \brief Definition of trigger bit configurations
   *
   * This enumeration handles different trigger bit configurations for the
   * EMCAL Level1 triggers (with and without different thresholds) applied
   * in the reconstruction of different samples.
   */
  enum TriggerMakerBitConfig_t {
    kOldConfig = 0,///< Old configuration, no distinction between high and low threshold
    kNewConfig = 1 ///< New configuration, distiction between high and low threshold
  };

  AliEmcalTriggerMakerTask();
  AliEmcalTriggerMakerTask(const char *name, Bool_t doQA = kFALSE);
  virtual ~AliEmcalTriggerMakerTask();

  virtual void      UserCreateOutputObjects();
  virtual void      ExecOnce();
  virtual Bool_t    Run();

  void SetRunQA(Bool_t doQA = kTRUE) { fDoQA = doQA; }
  void SetCaloTriggersOutName(const char *name)     { fCaloTriggersOutName      = name; }
  void SetCaloTriggerSetupOutName(const char *name) { fCaloTriggerSetupOutName  = name; }
  void SetV0InName(const char *name) { fV0InName      = name; }
  void SetUseTriggerBitConfig(TriggerMakerBitConfig_t bitConfig) { fUseTriggerBitConfig = bitConfig; }


protected:

  AliEmcalTriggerMakerKernel              *fTriggerMaker;             ///< The actual trigger maker kernel

  AliEmcalTriggerBitConfig                *fTriggerBitConfig;         //!<! Trigger bit configuration
  AliVVZERO                               *fV0;                       //!<! VZERO data

  TString                                 fCaloTriggersOutName;       ///< name of output track array
  TString                                 fCaloTriggerSetupOutName;   ///< name of output track array
  TString                                 fV0InName;                  ///< name of output track array
  TriggerMakerBitConfig_t                 fUseTriggerBitConfig;       ///< type of trigger config
  TClonesArray                            *fCaloTriggersOut;          //!<! trigger array out

  Bool_t                                  fDoQA;                      ///< Fill QA histograms
  THistManager                            *fQAHistos;                 //!<! Histograms for QA

private:
  AliEmcalTriggerMakerTask(const AliEmcalTriggerMakerTask &);
  AliEmcalTriggerMakerTask &operator=(const AliEmcalTriggerMakerTask &);

  ClassDef(AliEmcalTriggerMakerTask, 1);
};

#endif
