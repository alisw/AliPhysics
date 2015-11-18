#ifndef ALIEMCALTRIGGERMAKERKERNEL_H
#define ALIEMCALTRIGGERMAKERKERNEL_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

#include "AliEmcalTriggerChannelContainerAP.h"

class TObjArray;
class AliEmcalTriggerPatchInfo;
class AliEMCALGeometry;
class AliVCaloCells;
class AliVCaloTrigger;
class AliVEvent;
class AliVVZERO;
template<class T> class AliEmcalTriggerDataGridAP;
template<class T> class AliEmcalTriggerAlgorithmAP;
template<class T> class AliEmcalTriggerPatchFinderAP;

class AliEmcalTriggerMakerKernel : public TObject {
public:
  AliEmcalTriggerMakerKernel();
  virtual ~AliEmcalTriggerMakerKernel();

  void Init();
  TObjArray *CreateTriggerPatches(const AliVEvent *inputevent);

  void SetTriggerThresholdJetLow   ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[2][0] = a; fThresholdConstants[2][1] = b; fThresholdConstants[2][2] = c; }
  void SetTriggerThresholdJetHigh  ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[0][0] = a; fThresholdConstants[0][1] = b; fThresholdConstants[0][2] = c; }
  void SetTriggerThresholdGammaLow ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[3][0] = a; fThresholdConstants[3][1] = b; fThresholdConstants[3][2] = c; }
  void SetTriggerThresholdGammaHigh( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[3][0] = a; fThresholdConstants[1][1] = b; fThresholdConstants[1][2] = c; }

  void SetMC(Bool_t isMC) { fIsMC = isMC; }
  void SetRunNumber(Int_t runnumber) { fRunNumber = runnumber; }

  void SetGeometry(const AliEMCALGeometry *const geo) { fGeometry = geo; }
  void SetTriggerBitConfig(const AliEmcalTriggerBitConfig *const config) { fTriggerBitConfig = config; }

  /**
   * Switch on rejection of patches which leave the EMCAL acceptance in \f$ \eta \f$ and \f$ \phi \f$
   * \param doReject If true we reject patches outside the EMCAL acceptance
   */
  void SetRejectOffAcceptancePatches(Bool_t doReject = kTRUE) { fRejectOffAcceptancePatches = doReject; }

  /**
   * Reset data grids
   */
  void Reset();

  /**
   * Read the calo trigger data
   * @param trigger Input calo trigger data
   */
  void ReadTriggerData(AliVCaloTrigger *trigger);

  /**
   * Read the EMCAL cell data
   * @param cells EMCAL cell data
   */
  void ReadCellData(AliVCaloCells *cells);

  /**
   * Build VZERO-dependent thresholds for the offline trigger
   * @param vzdata VERO charges
   */
  void BuildL1ThresholdsOffline(const AliVVZERO *vzdata);

protected:
  enum{
    kColsEta = 48
  };

  Bool_t                                    CheckForL0(Int_t col, Int_t row) const;

  AliEmcalTriggerChannelContainerAP         fBadChannels;                 ///< Container of bad channels
  const AliEmcalTriggerBitConfig            *fTriggerBitConfig;           ///< Trigger bit configuration, aliroot-dependent
  const AliEMCALGeometry                    *fGeometry;                   //!<! Underlying EMCAL geometry

  AliEmcalTriggerDataGridAP<double>         *fPatchAmplitudes;            //!<! TRU Amplitudes (for L0)
  AliEmcalTriggerDataGridAP<double>         *fPatchADCSimple;             //!<! patch map for simple offline trigger
  AliEmcalTriggerDataGridAP<double>         *fPatchADC;                   //!<! ADC values map
  AliEmcalTriggerDataGridAP<char>           *fLevel0TimeMap;              //!<! Map needed to store the level0 times
  AliEmcalTriggerDataGridAP<int>            *fTriggerBitMap;              //!<! Map of trigger bits

  AliEmcalTriggerPatchFinderAP<double>      *fPatchFinder;                //!<! The actual patch finder
  AliEmcalTriggerAlgorithmAP<double>        *fLevel0PatchFinder;          //!<! Patch finder for Level0 patches

  Int_t                                     fThresholdConstants[4][3];    ///< simple offline trigger thresholds constants
  ULong64_t                                 fL1ThresholdsOffline[4];      ///< container for V0-dependent offline thresholds
  Bool_t                                    fIsMC;                        ///< Switch between data and MC mode
  Int_t                                     fRunNumber;                   ///< Run number
  Bool_t                                    fRejectOffAcceptancePatches;  ///< Switch for rejection of patches outside the acceptance
  Int_t                                     fDebugLevel;                  ///< Debug lebel;


  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerMakerKernel, 1);
  /// \endcond
};

#endif
