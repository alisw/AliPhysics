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

/**
 * @class AliEmcalTriggerMakerKernel
 * @brief Kernel of the EMCAL trigger patch maker
 *
 * The trigger maker kernel contains the core functionality of
 * the trigger maker:
 * - Filling of the data grids
 * - Steering and running the patch finders
 * - Conversion of the raw patches obtained by the patch finders
 *   to full EMCAL trigger patch info objects.
 * I/O, which means interaction with the ALICE analysis system,
 * is done outside in a dedicated steering task.
 */
class AliEmcalTriggerMakerKernel : public TObject {
public:

  /**
   * Constructor
   */
  AliEmcalTriggerMakerKernel();

  /**
   * Destructor
   */
  virtual ~AliEmcalTriggerMakerKernel();

  /**
   * Initialize the trigger maker Kernel:
   * - Allocates the data grids according to the EMCAL setup for
   *   the given run (obtained from the geometry)
   * - Initialize patch finders: Different patch finders are created
   *   for EMCAL and DCAL both for gamma and jet patches, using the
   *   data grids only in the active area of the EMCAL or DCAL
   */
  void Init();

  /**
   * Run patch finders on input data. Patches are converted from raw patches into AliEMCALTriggerPatchInfo data.
   * Trigger patches contain all information of the given category;
   * - Bit selection map from the STU (bits for non-matching patch types are removed)
   * - Offline trigger selection status
   * - Recalculation trigger selection status
   * There are no more several patches for the three categories. The only difference left over is the distinction of
   * patch types according to
   * - Jet patches (16x16 FAST-ors)
   * - Gamma patches (2x2 FAST-ors)
   * - Level0 patches (2x2 FAST-ors, using L0 amplitude and L0 times for the selection)
   * @param inputevent Input ESD/AOD event, used for kinematics calculation
   * @return Array of reconstructed trigger patches
   */
  TObjArray *CreateTriggerPatches(const AliVEvent *inputevent);

  void SetTriggerThresholdJetLow   ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[2][0] = a; fThresholdConstants[2][1] = b; fThresholdConstants[2][2] = c; }
  void SetTriggerThresholdJetHigh  ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[0][0] = a; fThresholdConstants[0][1] = b; fThresholdConstants[0][2] = c; }
  void SetTriggerThresholdGammaLow ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[3][0] = a; fThresholdConstants[3][1] = b; fThresholdConstants[3][2] = c; }
  void SetTriggerThresholdGammaHigh( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[1][0] = a; fThresholdConstants[1][1] = b; fThresholdConstants[1][2] = c; }

  /**
   * Define whether running on MC or not (for offset)
   * @param isMC Flag for MC
   */
  void SetIsMC(Bool_t isMC) { fIsMC = isMC; }

  /**
   * Provide the EMCAL geometry to the trigger maker Kernel
   * @param geo
   */
  void SetGeometry(const AliEMCALGeometry *const geo) { fGeometry = geo; }

  /**
   * Set the trigger bit configuration applied for the given data set
   * @param config Trigger bit config applied
   */
  void SetTriggerBitConfig(const AliEmcalTriggerBitConfig *const config) { fTriggerBitConfig = config; }

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

  /**
   * Accept trigger patch as Level0 patch. Level0 patches are identified as 2x2 FASTOR patches
   * in the same TRU
   * @param trg Triggers object with the pointer set to the patch to inspect
   * @return True if the patch is accepted, false otherwise.
   */
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
  Bool_t                                    fIsMC;                        ///< Set MC offset
  Int_t                                     fDebugLevel;                  ///< Debug lebel;

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerMakerKernel, 1);
  /// \endcond
};

#endif
