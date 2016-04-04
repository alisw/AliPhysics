/// \class AliEmcalTriggerQATaskPP
/// \brief Class to do some fast QA of the EMCal trigger (pp collisions). Useful also to tune trigger thresholds.
///
/// Class to do some fast QA of the EMCal trigger. Useful also to tune trigger thresholds.
/// The input for the process are the trigger patches AliEmcalTriggerPatchInfo produced by the AliEmcalTriggerMaker class.
/// This task is optimized for pp collisions.
///
/// The output is a bunch of histograms
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Feb 20, 2016

#ifndef ALIEMCALTRIGGERQATASKPP_H
#define ALIEMCALTRIGGERQATASKPP_H

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
class AliEmcalTriggerQAPP;
class THnSparse;
class AliESDEvent;

#include "AliLog.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliEMCALTriggerChannelContainer.h"

/**
 * \class AliEmcalTriggerQATask
 * \brief EMCAL trigger QA task
 *
 * This Class does some fast QA of the EMCal trigger. Useful also to tune trigger thresholds.
 */
class AliEmcalTriggerQATaskPP : public AliAnalysisTaskEmcal {
 public:

  AliEmcalTriggerQATaskPP();
  AliEmcalTriggerQATaskPP(const char *name);
  virtual ~AliEmcalTriggerQATaskPP();

  void SetTriggerPatchesName(const char *name)      { fTriggerPatchesName      = name; }
  void SetADCperBin(Int_t n);
  void SetMinAmplitude(Int_t m)                     { fMinAmplitude            = m   ; }
  void EnableDCal(Bool_t e = kTRUE)                 { fDCalPlots               = e   ; }
  void SetTimeStampRange(UInt_t min, UInt_t max)    { fMinTimeStamp            = min ; fMaxTimeStamp = max; }
  void EnableHistogramsByTimeStamp(UInt_t binWidth = 600){ fTimeStampBinWidth  = binWidth   ; }

  AliEmcalTriggerQAPP* GetTriggerQA(Int_t i = 0)    { return i >= 0 && i < fNcentBins ? static_cast<AliEmcalTriggerQAPP*>(fEMCALTriggerQA->At(i)) : 0; }

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
  AliEmcalTriggerQATaskPP(const AliEmcalTriggerQATaskPP&);            // not implemented
  AliEmcalTriggerQATaskPP &operator=(const AliEmcalTriggerQATaskPP&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerQATaskPP, 3)
  /// \endcond
};

#endif
