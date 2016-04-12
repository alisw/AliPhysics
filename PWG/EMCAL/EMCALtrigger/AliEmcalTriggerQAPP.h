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
#include <TArrayF.h>
#include "THistManager.h"

class AliEMCALTriggerPatchInfo;
class TObjArray;
class THashList;
class AliEMCALTriggerFastOR;
class AliEMCALGeometry;
class AliVCaloCells;

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
  void   AddFastORBadChannel(Short_t absId) { fBadChannels.insert(absId); }
  void   ReadFastORBadChannelFromStream(std::istream& stream);
  void   ReadFastORBadChannelFromFile(const char* fname);
  void   AddOfflineBadChannel(Short_t absId)  { fOfflineBadChannels.insert(absId)   ; }
  void   ReadOfflineBadChannelFromFile(const char* fname);
  void   ReadOfflineBadChannelFromStream(std::istream& stream);
  void   ReadFastORPedestalFromStream(std::istream& stream);
  void   ReadFastORPedestalFromFile(const char* fname);
  void   SetFastORPedestal(Short_t absId, Float_t ped);
  void   ResetFastORPedestal() { fFastORPedestal.Reset(); }

  Int_t  GetDebugLevel()        const { return fDebugLevel    ; }
  void   SetFastORandCellThresholds(Int_t l0, Int_t l1, Double_t cell) { fMinL0FastORAmp = l0; fMinL1FastORAmp = l1; fMinCellAmp = cell; }

  void   EnablePatchType(PatchTypes_t type, Bool_t e = kTRUE);
  void   EnableTriggerType(EMCalTriggerType_t type, Bool_t e = kTRUE);
  void   EnableDCal(Bool_t e = kTRUE) { fDCalPlots = e; }
  void   EnableHistogramsByTimeStamp(UInt_t binWidth = 600){ fTimeStampBinWidth  = binWidth   ; }

  void   Init();
  void   ExecOnce();
  void   ProcessPatch(AliEMCALTriggerPatchInfo* patch);
  void   ProcessFastor(AliEMCALTriggerFastOR* fastor, AliVCaloCells* cells);
  void   ProcessCell(const AliEmcalCellInfo& cell);
  void   EventCompleted();
  void   EventTimeStamp(UInt_t timeStamp);

  static Int_t  GetAmplitude(AliEMCALTriggerPatchInfo* patch, Int_t itype);

  THashList* GetListOfHistograms()  { return fHistManager.GetListOfHistograms(); }

  static const Int_t      fgkMaxPatchAmp[6];            ///< Maximum patch amplitude for the histograms
  static const TString    fgkPatchTypes[3];             ///< Patch type names

protected:

  std::set<Short_t>          fOfflineBadChannels;          ///< Abs ID of offline bad channels
  std::set<Short_t>          fBadChannels;                 ///< Container of bad channels
  TArrayF                    fFastORPedestal;              ///< FastOR pedestal
  Bool_t                     fEnabledPatchTypes[3];        ///< Patch types to be plotted
  Bool_t                     fEnabledTriggerTypes[6];      ///< Trigger types to be plotted
  Int_t                      fFastorL0Th;                  ///< FastOR L0 threshold
  Int_t                      fFastorL1Th;                  ///< FastOR L1 threshold
  Int_t                      fADCperBin;                   ///< ADC counts per bin
  Int_t                      fDebugLevel;                  ///< Debug level
  Int_t                      fL0MinTime;                   ///< Minimum L0 time
  Int_t                      fL0MaxTime;                   ///< Maximum L0 time
  Double_t                   fMinCellAmp;                  ///< Minimum offline amplitude of the cells
  Int_t                      fMinL0FastORAmp;              ///< Minimum L0 amplitude of the FastORs
  Int_t                      fMinL1FastORAmp;              ///< Minimum L1 amplitude of the FastORs
  THistManager               fHistManager;                 ///< Histogram manager
  Bool_t                     fDCalPlots;                   ///< Whether to add DCal QA plots
  UInt_t                     fTimeStampBinWidth;           ///< Time stamp bin width

  AliEMCALGeometry          *fGeom;                        //!<! EMCal geometry
  AliEMCALTriggerPatchInfo  *fMaxPatchEMCal[6][3];         //!<! EMCal max patch (will be reset each event)
  AliEMCALTriggerPatchInfo  *fMaxPatchDCal[6][3];          //!<! DCal max patch (will be reset each event)
  Double_t                   fSumOfflineEMCal;             //!<! EMCal sum of all offline energy deposition (will be reset each event)
  Int_t                      fSumL0EMCal;                  //!<! EMCal sum of all online energy deposition (will be reset each event)
  Int_t                      fSumL1EMCal;                  //!<! EMCal sum of all online energy deposition (will be reset each event)
  Double_t                   fSumOfflineDCal;              //!<! DCal sum of all offline energy deposition (will be reset each event)
  Int_t                      fSumL0DCal;                   //!<! DCal sum of all online energy deposition (will be reset each event)
  Int_t                      fSumL1DCal;                   //!<! DCal sum of all online energy deposition (will be reset each event)
  Int_t                      fNCellEMCal;                  //!<! EMCal number of offline cells (will be reset each event)
  Int_t                      fNL0EMCal;                    //!<! EMCal number of L0 FastORs (will be reset each event)
  Int_t                      fNL1EMCal;                    //!<! EMCal number of L1 FastORs (will be reset each event)
  Int_t                      fNCellDCal;                   //!<! DCal number of offline cells (will be reset each event)
  Int_t                      fNL0DCal;                     //!<! DCal number of L0 FastORs (will be reset each event)
  Int_t                      fNL1DCal;                     //!<! DCal number of L1 FastORs (will be reset each event)
  UInt_t                     fEventTimeStamp;              //!<! Time stamp of the current event
  UInt_t                     fEventTimeStampBin;           //!<! Time stamp bin
  Int_t                      fNTotTRU;                     //!<! Total number of TRUs
  Int_t                      fMaxFORabsId;                 //!<! Maximum FastOR abs id

private:
  AliEmcalTriggerQAPP &operator=(const AliEmcalTriggerQAPP &);

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerQAPP, 4);
  /// \endcond
};

#endif
