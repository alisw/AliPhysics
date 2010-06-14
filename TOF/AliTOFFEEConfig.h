#ifndef ALITOFFEECONFIG_H
#define ALITOFFEECONFIG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFDecoder.h,v 1.2 2007/05/08 11:55:24 arcelli Exp $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF FEE config object.         //
//                                                           //
///////////////////////////////////////////////////////////////

#include <time.h>
#include "TROOT.h"

class AliTOFDRMConfig
{

 private:

  Int_t fDRMId; // DRM id
  Int_t fOptLinkId; // opt link id
  Int_t fFineDelayTTC; // fine delay TTC
  Int_t fCoarseDelayTTC; // coarse delay TTC
  Int_t fSelectTTC; // select TTC
  Int_t fGhostDDL; // ghost DDL
  Int_t fActiveTRM; // active TRM
  Int_t fSlotMask; // slot mask
  Int_t fPulserPolarity; // pulser polarity
  Int_t fPulserFlag; // pulser flag
  Int_t fPulserSection0; // pulser section 0
  Int_t fPulserSection1; // pulser section 1
  Int_t fPulserSection2; // pulser section 2
  Int_t fPrePulseEnable; // pre pulse enable
  Int_t fSelectMode; // select mode
  Int_t fBLTMask; // BLT mask

 public:

  Int_t GetDRMId() const {return fDRMId;}; // get DRM id
  Int_t GetOptLink() const {return fOptLinkId;}; // get opt link id
  Int_t GetFineDelayTTC() const {return fFineDelayTTC;}; // get fine delay TTC
  Int_t GetCoarseDelayTTC() const {return fCoarseDelayTTC;}; // get coarse delay TTC
  Int_t GetSelectTTC() const {return fSelectTTC;}; // get select TTC
  Int_t GetGhostDDL() const {return fGhostDDL;}; // get ghost DDL
  Int_t GetActiveTRM() const {return fActiveTRM;}; // get active TRM
  Int_t GetSlotMask() const {return fSlotMask;}; // get slot mask
  Int_t GetPulserPolarity() const {return fPulserPolarity;}; // get pulser polarity
  Int_t GetPulserFlag() const {return fPulserFlag;}; // get pulser flag
  Int_t GetPulserSection0() const {return fPulserSection0;}; // get pulser section 0
  Int_t GetPulserSection1() const {return fPulserSection1;}; // get pulser section 1
  Int_t GetPulserSection2() const {return fPulserSection2;}; // get pulser section 2
  Int_t GetPrePulserEnable() const {return fPrePulseEnable;}; // get pre pulse enable
  Int_t GetSelectMode() const {return fSelectMode;}; // get select mode
  Int_t GetBLTMask() const {return fBLTMask;}; // get BLT mask

};

class AliTOFLTMConfig
{

 private:

  Int_t fSlotId; // slot id
  Int_t fSuperModuleId; // super module id
  Int_t fOptLinkId; // opt link id
  Int_t fVMEAddress; // VME address
  Int_t fThreshold; // threshold
  Int_t fTriggerEnable; // trigger enable
  Int_t fTriggerMask1; // trigger mask1
  Int_t fTriggerMask2; // trigger mask2
#if 0
  Int_t fDelayLine0; // delay line 0
  Int_t fDelayLine1; // delay line 1
  Int_t fDelayLine2; // delay line 2
  Int_t fDelayLine3; // delay line 3
  Int_t fDelayLine4; // delay line 4
  Int_t fDelayLine5; // delay line 5
  Int_t fDelayLine6; // delay line 6
  Int_t fDelayLine7; // delay line 7
#endif

 public:

  Int_t GetSlotId() const {return fSlotId;}; // get slot id
  Int_t GetSuperModuleId() const {return fSuperModuleId;}; // get super module id
  Int_t GetOptLinkId() const {return fOptLinkId;}; // get opt link id
  Int_t GetVMEAddress() const {return fVMEAddress;}; // get VME address
  Int_t GetThreshold() const {return fThreshold;}; // get threshold
  Int_t GetTriggerEnable() const {return fTriggerEnable;}; // get trigger enable
  Int_t GetTriggerMask1() const {return fTriggerMask1;}; // get trigger mask1
  Int_t GetTriggerMask2() const {return fTriggerMask2;}; // get trigger mask2
#if 0
  Int_t GetDelayLine0() const {return fDelayLine0;}; // get delay line 0
  Int_t GetDelayLine1() const {return fDelayLine1;}; // get delay line 1
  Int_t GetDelayLine2() const {return fDelayLine2;}; // get delay line 2
  Int_t GetDelayLine3() const {return fDelayLine3;}; // get delay line 3
  Int_t GetDelayLine4() const {return fDelayLine4;}; // get delay line 4
  Int_t GetDelayLine5() const {return fDelayLine5;}; // get delay line 5
  Int_t GetDelayLine6() const {return fDelayLine6;}; // get delay line 6
  Int_t GetDelayLine7() const {return fDelayLine7;}; // get delay line 7
#endif

};

class AliTOFTRMConfig
{

 private:

  Int_t fSlotId; // slot id
  Int_t fSuperModuleId; // super module id
  Int_t fOptLinkId; // opt link id
  Int_t fVMEAddress; // VME address
  Int_t fMatchingWindow; // matching window
  Int_t fLatencyWindow; // latency window
  Int_t fBunchCrossingAdjust; // bunch crossing adjust
  Int_t fTriggerLevelConfig; // trigger level config
  Int_t fTriggerSubtraction; // trigger subtracion
  Int_t fEdgeDetection; // edge detection
  Int_t fPackingFlag; // packing flag
  Int_t fLVStatus; // LV status
  Int_t fChainAFlag; // chain A flag
  Int_t fChainBFlag; // chain B flag
  Int_t fActiveChipA; // active chip A
  Int_t fActiveChipB; // active chip B
  Int_t fMaskPB0; // mask PB 0
  Int_t fMaskPB1; // mask PB 1
  Int_t fMaskPB2; // mask PB 2
  Int_t fMaskPB3; // mask PB 3
  Int_t fMaskPB4; // mask PB 4
  Int_t fMaskPB5; // mask PB 5
  Int_t fMaskPB6; // mask PB 6
  Int_t fMaskPB7; // mask PB 7
  Int_t fMaskPB8; // mask PB 8
  Int_t fMaskPB9; // mask PB 9

 public:

  Int_t GetSlotId() const {return fSlotId;}; // get slot id
  Int_t GetSuperModuleId() const {return fSuperModuleId;}; // get super module id
  Int_t GetOptLinkId() const {return fOptLinkId;}; // get opt link id
  Int_t GetVMEAddress() const {return fVMEAddress;}; // get VME address
  Int_t GetMatchingWindow() const {return fMatchingWindow;}; // get matching window
  Int_t GetLatencyWindow() const {return fLatencyWindow;}; // get latency window
  Int_t GetBunchCrossingAdjust() const {return fBunchCrossingAdjust;}; // get bunch crossing adjust
  Int_t GetTriggerLevelConfig() const {return fTriggerLevelConfig;}; // get trigger level config
  Int_t GetTriggerSubtraction() const {return fTriggerSubtraction;}; // get trigger subtracion
  Int_t GetEdgeDetection() const {return fEdgeDetection;}; // get edge detection
  Int_t GetPackingFlag() const {return fPackingFlag;}; // get packing flag
  Int_t GetLVStatus() const {return fLVStatus;}; // get LV status
  Int_t GetChainAFlag() const {return fChainAFlag;}; // get chain A flag
  Int_t GetChainBFlag() const {return fChainBFlag;}; // get chain B flag
  Int_t GetActiveChipA() const {return fActiveChipA;}; // get active chip A
  Int_t GetActiveChipB() const {return fActiveChipB;}; // get active chip B
  Int_t GetMaskPB0() const {return fMaskPB0;}; // get mask PB 0
  Int_t GetMaskPB1() const {return fMaskPB1;}; // get mask PB 1
  Int_t GetMaskPB2() const {return fMaskPB2;}; // get mask PB 2
  Int_t GetMaskPB3() const {return fMaskPB3;}; // get mask PB 3
  Int_t GetMaskPB4() const {return fMaskPB4;}; // get mask PB 4
  Int_t GetMaskPB5() const {return fMaskPB5;}; // get mask PB 5
  Int_t GetMaskPB6() const {return fMaskPB6;}; // get mask PB 6
  Int_t GetMaskPB7() const {return fMaskPB7;}; // get mask PB 7
  Int_t GetMaskPB8() const {return fMaskPB8;}; // get mask PB 8
  Int_t GetMaskPB9() const {return fMaskPB9;}; // get mask PB 9

};


class AliTOFCrateConfig
{

 private:
  
  static const Int_t fgkNumberOfTRMs = 10; // number of TRMs

  enum EStatus_t {
    kStatusDisabled = 0,
    kStatusEnabled = 1
  };

  Int_t fDRMEnabled; // DRM enabled
  AliTOFDRMConfig fDRMConfig; // DRM config
  Int_t fLTMEnabled; // LTM enabled
  AliTOFLTMConfig fLTMConfig; // LTM config
  Int_t fCPDMEnabled; // CPDM enabled
  Int_t fTRMEnabled[fgkNumberOfTRMs]; // TRM enabled array
  AliTOFTRMConfig fTRMConfig[fgkNumberOfTRMs]; // TRM config array

 public:
  
  static Int_t GetNumberOfTRMs() {return fgkNumberOfTRMs;}; // get number of TRMs
  Int_t GetDRMEnabled() const {return fDRMEnabled;}; // get DRM enabled
  Bool_t IsDRMEnabled() const {return fDRMEnabled == kStatusEnabled;}; // is DRM enabled
  AliTOFDRMConfig *GetDRMConfig() {return &fDRMConfig;}; // get DRM config
  Int_t GetLTMEnabled() const {return fLTMEnabled;}; // get LTM enabled
  Bool_t IsLTMEnabled() const {return fLTMEnabled == kStatusEnabled;}; // is LTM enabled
  AliTOFLTMConfig *GetLTMConfig() {return &fLTMConfig;}; // get LTM config
  Int_t GetCPDMEnabled() const {return fCPDMEnabled;}; // get CPDM enabled
  Bool_t IsCPDMEnabled() const {return fCPDMEnabled == kStatusEnabled;}; // is CPDM enabled
  Int_t GetTRMEnabled(UShort_t iTRM) const {return iTRM < GetNumberOfTRMs() ? fTRMEnabled[iTRM] : 0;}; // get TRM enabled
  Bool_t IsTRMEnabled(UShort_t iTRM) const {return iTRM < GetNumberOfTRMs() ? fTRMEnabled[iTRM] ==  kStatusEnabled : 0;}; // is TRM enabled
  AliTOFTRMConfig *GetTRMConfig(UShort_t iTRM) {return iTRM < GetNumberOfTRMs() ? &fTRMConfig[iTRM] : NULL;}; // get TRM config

};

class AliTOFFEEConfig
{

 private:

  static const Int_t fgkNumberOfCrates = 72; // number of crates

  Int_t fVersion; // version
  time_t fDumpTime; // dump time
  Int_t fRunNumber; // run number
  Int_t fRunType; // run type
  Int_t fBytes; // bytes
  Int_t fCTTMTriggerMask[fgkNumberOfCrates]; // CTTM trigger mask
  AliTOFCrateConfig fCrateConfig[fgkNumberOfCrates]; // crate config array

 public:

  static Int_t GetNumberOfCrates() {return fgkNumberOfCrates;}; // get number of crates
  Int_t GetVersion() const {return fVersion;}; // get version
  time_t GetDumpTime() const {return fDumpTime;}; // get dump time
  Int_t GetRunNumber() const {return fRunNumber;}; // get run number
  Int_t GetRunType() const {return fRunType;}; // get run type
  Int_t GetBytes() const {return fBytes;}; // get bytes
  Int_t GetCTTMTriggerMask(UShort_t iCrate) const {return iCrate < GetNumberOfCrates() ? fCTTMTriggerMask[iCrate] : NULL;}; // get CTTM trigger mask
  AliTOFCrateConfig *GetCrateConfig(UShort_t iCrate) {return iCrate < GetNumberOfCrates() ? &fCrateConfig[iCrate] : NULL;}; // get crate config

};

#endif /* ALITOFFEE_H */
