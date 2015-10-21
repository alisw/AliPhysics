#ifndef ALIEMCALTRIGGERMAKERKERNEL_H
#define ALIEMCALTRIGGERMAKERKERNEL_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

#include "AliEmcalTriggerChannelContainer.h"

class TObjArray;
class AliEmcalTriggerPatchInfo;
class AliEmcalTriggerRawPatch;
class AliEMCALGeometry;
class AliVCaloCells;
class AliVCaloTrigger;
class AliVEvent;
class AliVVZERO;
template<class T> class AliEmcalTriggerDataGrid;

class AliEmcalTriggerMakerKernel : public TObject {
public:
  /**
   * \enum TriggerMakerTriggerType_t
   * \brief Definition of different trigger patch types
   *
   * This enumeration defines the different trigger patch types
   * processed by the trigger maker. Each trigger patch type has
   * a certain patch size and therefore a certain length and
   * geometric center
   */
  enum TriggerMakerTriggerType_t {
    kTMEMCalJet = 0,            ///< EMCAL Jet trigger
    kTMEMCalGamma = 1,          ///< EMCAL Gamma trigger
    kTMEMCalLevel0 = 2,         ///< EMCAL Level0 patches
    kTMEMCalRecalcJet = 3,      ///< EMCAL Jet patches, recalculated
    kTMEMCalRecalcGamma = 4,    ///< EMCAL Gamma patches, recalculated
    kTMUndefined = -1
  };

  AliEmcalTriggerMakerKernel();
  virtual ~AliEmcalTriggerMakerKernel();

  void Init();
  TObjArray *CreateTriggerPatches(const AliVEvent *inputevent);

  inline Bool_t IsEJE(Int_t tBits) const;
  inline Bool_t IsEGA(Int_t tBits) const;
  inline Bool_t IsLevel0(Int_t tBits) const;

  void SetTriggerThresholdJetLow   ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[2][0] = a; fThresholdConstants[2][1] = b; fThresholdConstants[2][2] = c; }
  void SetTriggerThresholdJetHigh  ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[0][0] = a; fThresholdConstants[0][1] = b; fThresholdConstants[0][2] = c; }
  void SetTriggerThresholdGammaLow ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[3][0] = a; fThresholdConstants[3][1] = b; fThresholdConstants[3][2] = c; }
  void SetTriggerThresholdGammaHigh( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[3][0] = a; fThresholdConstants[1][1] = b; fThresholdConstants[1][2] = c; }

  void SetCaloTriggers(AliVCaloTrigger *triggers) { fCaloTriggers = triggers; }
  void SetCaloCells(AliVCaloCells *cells) { fCaloCells = cells; }
  void SetVZERO(AliVVZERO *vzero) { fV0 = vzero; }

  void SetMC(Bool_t isMC) { fIsMC = isMC; }
  void SetRunNumber(Int_t runnumber) { fRunNumber = runnumber; }

  void SetGeometry(const AliEMCALGeometry *const geo) { fGeometry = geo; }
  void SetTriggerBitConfig(const AliEmcalTriggerBitConfig *const config) { fTriggerBitConfig = config; }
  void SetRunTriggerType(TriggerMakerTriggerType_t type, Bool_t doTrigger = kTRUE) { fRunTriggerType[type] = doTrigger; }

  /**
   * Switch on rejection of patches which leave the EMCAL acceptance in \f$ \eta \f$ and \f$ \phi \f$
   * \param doReject If true we reject patches outside the EMCAL acceptance
   */
  void SetRejectOffAcceptancePatches(Bool_t doReject = kTRUE) { fRejectOffAcceptancePatches = doReject; }

  static const TString &GetTriggerTypeName(int index) { return fgkTriggerTypeNames[index]; }

  AliEmcalTriggerPatchInfo *ConvertRawPatch(const AliEmcalTriggerRawPatch *patch) const;

protected:
  static const TString                      fgkTriggerTypeNames[5];       ///< Histogram name tags
  static const int                          kColsEta;                     ///< Number of columns in eta direction

  void                                      RunSimpleOfflineTrigger();
  Bool_t                                    NextTrigger( Bool_t &isOfflineSimple );
  AliEmcalTriggerPatchInfo*                 ProcessPatch(TriggerMakerTriggerType_t type, Bool_t isOfflineSimple);
  Bool_t                                    CheckForL0(const AliVCaloTrigger &trg) const;

  AliEmcalTriggerChannelContainer           fBadChannels;                 ///< Container of bad channels
  const AliEmcalTriggerBitConfig            *fTriggerBitConfig;           ///< Trigger bit configuration, aliroot-dependent
  const AliEMCALGeometry                    *fGeometry;                   //!<! Underlying EMCAL geometry
  AliEmcalTriggerSetupInfo                  *fCaloTriggerSetupOut;        //!<! trigger setup

  AliEmcalTriggerDataGrid<float>            *fPatchAmplitudes;            //!<! TRU Amplitudes (for L0)
  AliEmcalTriggerDataGrid<double>           *fPatchADCSimple;             //!<! patch map for simple offline trigger
  AliEmcalTriggerDataGrid<int>              *fPatchADC;                   //!<! ADC values map
  AliEmcalTriggerDataGrid<char>             *fLevel0TimeMap;              //!<! Map needed to store the level0 times

  // Temporary objects from AliAnalysisTaskEmcal - will be refactored
  AliVCaloCells                             *fCaloCells;
  AliVCaloTrigger                           *fCaloTriggers;
  AliVVZERO                                 *fV0;
  AliAODCaloTrigger                         *fSimpleOfflineTriggers;      //!<! simple offline trigger
  Double_t                                  fVertex[3];

  Int_t                                     fThresholdConstants[4][3];    ///< simple offline trigger thresholds constants
  Bool_t                                    fIsMC;                        ///< Switch between data and MC mode
  Int_t                                     fRunNumber;                   ///< Run number
  Bool_t                                    fRunTriggerType[5];           ///< Run patch maker for a given trigger type
  Bool_t                                    fRejectOffAcceptancePatches;  ///< Switch for rejection of patches outside the acceptance
  Int_t                                     fDebugLevel;                  ///< Debug lebel;


  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerMakerKernel, 1);
  /// \endcond
};

/**
 * Check whehter trigger is jet patch trigger according to trigger bits
 * @param tBits Trigger bits of the fastor
 * @return True if fastor belongs to a jet patch trigger, false otherwise
 */
Bool_t AliEmcalTriggerMakerKernel::IsEJE(Int_t tBits) const {
  if( tBits & ( 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetJetLowBit()) | 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetJetHighBit()) | 1 << (fTriggerBitConfig->GetJetLowBit()) | 1 << (fTriggerBitConfig->GetJetHighBit()) ))
    return kTRUE;
  else
    return kFALSE;
}

/**
 * Check whether trigger is gamma patch trigger according to trigger bits
 * @param tBits Trigger bits of the fastor
 * @return True if fastor belongs to a gamma trigger, false otherwise
 */
Bool_t AliEmcalTriggerMakerKernel::IsEGA(Int_t tBits) const {
  if( tBits & ( 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetGammaLowBit()) | 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetGammaHighBit()) | 1 << (fTriggerBitConfig->GetGammaLowBit()) | 1 << (fTriggerBitConfig->GetGammaHighBit()) ))
    return kTRUE;
  else
    return kFALSE;
}

/**
 * Check whether trigger is level0 trigger according to trigger bits
 * @param tBits Trigger bits of the fastor
 * @return True if fastor belongs to a level0 trigger, false otherwise
 */
Bool_t AliEmcalTriggerMakerKernel::IsLevel0(Int_t tBits) const {
  if( tBits & (1 << (fTriggerBitConfig->GetLevel0Bit() + fTriggerBitConfig->GetLevel0Bit()) | (1 << fTriggerBitConfig->GetLevel0Bit())))
    return kTRUE;
  return kFALSE;
}

#endif
