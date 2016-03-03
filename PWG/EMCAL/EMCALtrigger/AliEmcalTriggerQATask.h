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
class TObjArray;
class THistManager;
class TString;
class THnSparse;

#include "AliEMCALTriggerQA.h"
#include "AliLog.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliEMCALTriggerChannelContainer.h"

/**
 * \class AliEmcalTriggerQATask
 * \brief EMCAL trigger QA task
 *
 * This Class does some fast QA of the EMCal trigger. Useful also to tune trigger thresholds.
 */
class AliEmcalTriggerQATask : public AliAnalysisTaskEmcal {
 public:

  enum CaloTriggers_t {
    kEMCalL0,
    kEMCalL1G1,
    kEMCalL1G2,
    kEMCalL1J1,
    kEMCalL1J2,
    kDCalL0,
    kDCalL1G1,
    kDCalL1G2,
    kDCalL1J1,
    kDCalL1J2,
    kPHOSL0,
    kPHOSL1H,
    kPHOSL1M,
    kPHOSL1L,
    kMinBias,
    kLastCaloTrigger
  };

  enum CaloTriggerBits_t {
    kEMCalL0bit = BIT(kEMCalL0),
    kEMCalL1G1bit = BIT(kEMCalL1G1),
    kEMCalL1G2bit = BIT(kEMCalL1G2),
    kEMCalL1J1bit = BIT(kEMCalL1J1),
    kEMCalL1J2bit = BIT(kEMCalL1J2),
    kDCalL0bit = BIT(kDCalL0),
    kDCalL1G1bit = BIT(kDCalL1G1),
    kDCalL1G2bit = BIT(kDCalL1G2),
    kDCalL1J1bit = BIT(kDCalL1J1),
    kDCalL1J2bit = BIT(kDCalL1J2),
    kPHOSL0bit = BIT(kPHOSL0),
    kPHOSL1Hbit = BIT(kPHOSL1H),
    kPHOSL1Mbit = BIT(kPHOSL1M),
    kPHOSL1Lbit = BIT(kPHOSL1L),
    kMinBiasbit = BIT(kMinBias),

    kEMCalL1Anybit = kEMCalL1G1bit | kEMCalL1G2bit | kEMCalL1J1bit | kEMCalL1J2bit,
    kEMCalAnybit = kEMCalL0bit | kEMCalL1Anybit,

    kDCalL1Anybit = kDCalL1G1bit | kDCalL1G2bit | kDCalL1J1bit | kDCalL1J2bit,
    kDCalAnybit = kDCalL0bit | kDCalL1Anybit,

    kEMCalDCalL0bit = kEMCalL0bit | kDCalL0bit,
    kEMCalDCalL1G1bit = kEMCalL1G1bit | kDCalL1G1bit,
    kEMCalDCalL1G2bit = kEMCalL1G2bit | kDCalL1G2bit,
    kEMCalDCalL1J1bit = kEMCalL1J1bit | kDCalL1J1bit,
    kEMCalDCalL1J2bit = kEMCalL1J2bit | kDCalL1J2bit,

    kEMCalDCalL1Anybit = kEMCalDCalL1G1bit | kEMCalDCalL1G2bit | kEMCalDCalL1J1bit | kEMCalDCalL1J2bit,
    kEMCalDCalAnybit = kEMCalDCalL0bit | kEMCalDCalL1Anybit,

    kPHOSL1Anybit = kPHOSL1Hbit | kPHOSL1Mbit | kPHOSL1Lbit,
    kPHOSAnybit = kPHOSL0bit | kPHOSL1Anybit,

    kCALOL0bit = kEMCalDCalL0bit | kPHOSL0bit,
    kCALOL1bit = kEMCalDCalL1Anybit | kPHOSL1Anybit,
    kCALOAnybit = kCALOL0bit | kCALOL1bit,

    kCALOMinBias = kCALOAnybit | kMinBiasbit
  };

  AliEmcalTriggerQATask();
  AliEmcalTriggerQATask(const char *name);
  virtual ~AliEmcalTriggerQATask();

  void SetTriggerPatchesName(const char *name)      { fTriggerPatchesName      = name; }
  void SetBkgPatchType(Int_t t);
  void SetADCperBin(Int_t n);
  void SetNCentBins(Int_t n);

  AliEMCALTriggerQA* GetTriggerQA(Int_t i = 0)    { return i >= 0 && i < fNcentBins ? static_cast<AliEMCALTriggerQA*>(fEMCALTriggerQA->At(i)) : 0; }

  void Set2015CaloTriggerNames();

 protected:
  void                                      UserCreateOutputObjects();
  void                                      ExecOnce();
  Bool_t                                    Run();
  Bool_t                                    FillHistograms();
  void                                      FillEventQA();
  
  UInt_t                                    SteerFiredTriggers(const TString& firedTriggersStr) const;

  TString                                   fCaloTriggerNames[kLastCaloTrigger]; ///< names of the calo trigger classes

  TString                                   fTriggerPatchesName;         ///< name of input trigger array
  TObjArray                                *fEMCALTriggerQA;             ///< produces the QA histograms
  Int_t                                     fADCperBin;                  ///< ADC counts per bin
  Int_t                                     fBkgPatchType;               ///< Background patch type
  AliEMCALTriggerChannelContainer           fBadChannels;                ///< Container of bad channels

  TClonesArray                             *fTriggerPatches;             //!<! trigger array in
  TH1                                      *fHistEMCalTriggers;          //!<! EMCal triggers
  THnSparse                                *fHistEventQA;                //!<! Event QA

 private:
  AliEmcalTriggerQATask(const AliEmcalTriggerQATask&);            // not implemented
  AliEmcalTriggerQATask &operator=(const AliEmcalTriggerQATask&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerQATask, 1) // Task to make QA of EMCAL trigger
  /// \endcond
};

#endif
