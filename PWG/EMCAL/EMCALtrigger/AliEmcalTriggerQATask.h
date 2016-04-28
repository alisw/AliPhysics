/// \class AliEmcalTriggerQATask
/// \brief Class to do some fast QA of the EMCal trigger. Useful also to tune trigger thresholds.
///
/// Class to do some fast QA of the EMCal trigger. Useful also to tune trigger thresholds.
/// The input for the process are the trigger patches AliEMCALTriggerPatchInfo produced by the AliEmcalTriggerMaker class.
///
/// The output is a bunch of histograms
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Apr 4, 2016

#ifndef ALIEMCALTRIGGERQATASK_H
#define ALIEMCALTRIGGERQATASK_H

/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

class TClonesArray;
class TObjArray;
class THistManager;
class TString;
class AliEMCALTriggerQA;
class THnSparse;
class AliESDEvent;

#include <AliEMCALTriggerQA.h>
#include <AliLog.h>
#include <AliEMCALTriggerChannelContainer.h>

#include "AliAnalysisTaskEmcalLight.h"

/**
 * \class AliEmcalTriggerQATask
 * \brief EMCAL trigger QA task
 *
 * This Class does some fast QA of the EMCal trigger. Useful also to tune trigger thresholds.
 */
class AliEmcalTriggerQATask : public AliAnalysisTaskEmcalLight {
 public:

  AliEmcalTriggerQATask();
  AliEmcalTriggerQATask(const char *name, UInt_t nCentBins=0, Bool_t online=kFALSE);
  virtual ~AliEmcalTriggerQATask();

  void SetTriggerPatchesName(const char *name)      { fTriggerPatchesName      = name; }
  void SetADCperBin(Int_t n);
  void SetMinAmplitude(Int_t m)                     { fMinAmplitude            = m   ; }
  void EnableDCal(Bool_t e = kTRUE)                 { fDCalPlots               = e   ; }
  void SetTimeStampRange(UInt_t min, UInt_t max)    { fMinTimeStamp            = min ; fMaxTimeStamp = max; }
  void EnableHistogramsByTimeStamp(UInt_t binWidth = 600){ fTimeStampBinWidth  = binWidth   ; }

  AliEMCALTriggerQA* GetTriggerQA(Int_t i = 0)    { return fEMCALTriggerQA && i >= 0 && i < fEMCALTriggerQA->GetEntriesFast() ? static_cast<AliEMCALTriggerQA*>(fEMCALTriggerQA->At(i)) : 0; }

 protected:
  void                                      UserCreateOutputObjects();
  void                                      ExecOnce();
  Bool_t                                    Run();
  Bool_t                                    FillHistograms();
  void                                      FillEventQA();

  TString                                   fTriggerPatchesName;         ///< name of input trigger array
  TObjArray                                *fEMCALTriggerQA;             ///< produces the QA histograms
  Int_t                                     fADCperBin;                  ///< ADC counts per bin
  Int_t                                     fMinAmplitude;               ///< Minimum trigger patch amplitude
  Bool_t                                    fDCalPlots;                  ///< Whether to add DCal QA plots
  UInt_t                                    fMinTimeStamp;               ///< Minimum event time stamp (only ESD)
  UInt_t                                    fMaxTimeStamp;               ///< Maximum event time stamp (only ESD)
  UInt_t                                    fTimeStampBinWidth;          ///< Time stamp bin width

  AliESDEvent                              *fESDEvent;                   //!<! current ESD event
  TClonesArray                             *fTriggerPatches;             //!<! trigger array in

 private:
  AliEmcalTriggerQATask(const AliEmcalTriggerQATask&);            // not implemented
  AliEmcalTriggerQATask &operator=(const AliEmcalTriggerQATask&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerQATask, 4);
  /// \endcond
};

#endif
