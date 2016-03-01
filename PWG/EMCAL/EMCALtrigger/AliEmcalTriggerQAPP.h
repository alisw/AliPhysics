/// \class AliEmcalTriggerQAPP
/// \brief Class to generate EMCal trigger QA plots in pp collisions
///
/// This class generates QA plots for EMCal trigger in pp collisions
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Feb. 20, 2016

#ifndef ALIEMCALTRIGGERQAPP_H
#define ALIEMCALTRIGGERQAPP_H

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

#include <set>
#include <iostream>
#include <TNamed.h>
#include "THistManager.h"

class AliEMCALTriggerPatchInfo;
class TObjArray;
class THashList;
class AliEMCALTriggerFastOR;

class AliEmcalTriggerQAPP : public TNamed {
public:

  typedef EMCALTrigger::EMCalTriggerType_t EMCalTriggerType_t;

  enum PatchTypes_t {
    kOnlinePatch,
    kRecalcPatch,
    kOfflinePatch
  };

  struct AliEmcalCellInfo {
    AliEmcalCellInfo() : fAbsId(-1), fEnergy(0.) {}
    void Set(Int_t absId, Double_t e) { fAbsId = absId; fEnergy = e; }

    Short_t  fAbsId;
    Double_t fEnergy;
  };

  AliEmcalTriggerQAPP();
  AliEmcalTriggerQAPP(const char* name);
  AliEmcalTriggerQAPP(const AliEmcalTriggerQAPP& triggerQA);
  virtual ~AliEmcalTriggerQAPP();

  void   SetDebugLevel(Int_t l)               { fDebugLevel = l; }
  void   SetADCperBin(Int_t i)                { fADCperBin  = i; }
  void   SetL0TimeRange(Int_t min, Int_t max) { fL0MinTime = min; fL0MaxTime = max; }
  void   AddOfflineBadChannel(Short_t absId)  { fOfflineBadChannels.insert(absId)   ; }
  void   ReadOfflineBadChannelFromFile(const char* fname);
  void   ReadOfflineBadChannelFromStream(std::istream& stream);

  Int_t  GetDebugLevel()        const { return fDebugLevel    ; }

  void   EnablePatchType(PatchTypes_t type, Bool_t e = kTRUE);
  void   EnableTriggerType(EMCalTriggerType_t type, Bool_t e = kTRUE);

  void   Init();
  void   ProcessPatch(AliEMCALTriggerPatchInfo* patch);
  void   ProcessFastor(AliEMCALTriggerFastOR* fastor);
  void   ProcessCell(const AliEmcalCellInfo& cell);
  void   EventCompleted();
  void   ComputeBackground();



  THashList* GetListOfHistograms()  { return fHistManager.GetListOfHistograms(); }

  static const Int_t      fgkMaxPatchAmp[6];            ///< Maximum patch amplitude for the histograms
  static const TString    fgkPatchTypes[3];             ///< Patch type names

protected:

  std::set<Short_t>       fOfflineBadChannels;          ///< Abd ID of offline bad channels
  Bool_t                  fEnabledPatchTypes[3];        ///< Patch types to be plotted
  Bool_t                  fEnabledTriggerTypes[6];      ///< Trigger types to be plotted
  Int_t                   fFastorL0Th;                  ///< FastOR L0 threshold
  Int_t                   fFastorL1Th;                  ///< FastOR L1 threshold
  Int_t                   fADCperBin;                   ///< ADC counts per bin
  Int_t                   fDebugLevel;                  ///< Debug level
  Int_t                   fL0MinTime;                   ///< Minimum L0 time
  Int_t                   fL0MaxTime;                   ///< Maximum L0 time

  THistManager            fHistManager;                 ///<  Histogram manager
  Int_t                   fMaxPatchEMCal[6][3];         //!<! EMCal max ADC amplitude (will be reset each event)
  Int_t                   fMaxPatchDCal[6][3];          //!<! DCal max ADC amplitude (will be reset each event)

private:
  AliEmcalTriggerQAPP &operator=(const AliEmcalTriggerQAPP &);

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerQAPP, 1);
  /// \endcond
};

#endif
