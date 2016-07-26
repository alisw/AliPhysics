#ifndef ALIEMCALTRIGGERMAKERKERNEL_H
#define ALIEMCALTRIGGERMAKERKERNEL_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <set>
#include <iostream>

#include <TObject.h>

#include "AliEMCALTriggerChannelContainer.h"

class TObjArray;
class AliEMCALTriggerPatchInfo;
class AliEMCALTriggerRawPatch;
class AliEMCALGeometry;
class AliVCaloCells;
class AliVCaloTrigger;
class AliVEvent;
class AliVVZERO;
template<class T> class AliEMCALTriggerDataGrid;
template<class T> class AliEMCALTriggerAlgorithm;
template<class T> class AliEMCALTriggerPatchFinder;

// To be moved to AliRoot in AliEMCALTriggerConstants.h at the first occasion
namespace EMCALTrigger {
const Double_t kEMCL0ADCtoGeV_AP = 0.018970588*4;  // 0.075882352;             ///< Conversion from EMCAL Level0 ADC to energy
}

/**
 * @class AliEmcalTriggerMakerKernel
 * @brief Kernel of the EMCAL trigger patch maker
 * @ingroup EMCALTRGFW
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

  enum ELevel0TriggerStatus_t { kNotLevel0, kLevel0Candidate, kLevel0Fired };

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
  TObjArray *CreateTriggerPatches(const AliVEvent *inputevent, Bool_t useL0amp=kFALSE);

  void SetTriggerThresholdJetLow   ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[2][0] = a; fThresholdConstants[2][1] = b; fThresholdConstants[2][2] = c; }
  void SetTriggerThresholdJetHigh  ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[0][0] = a; fThresholdConstants[0][1] = b; fThresholdConstants[0][2] = c; }
  void SetTriggerThresholdGammaLow ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[3][0] = a; fThresholdConstants[3][1] = b; fThresholdConstants[3][2] = c; }
  void SetTriggerThresholdGammaHigh( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[1][0] = a; fThresholdConstants[1][1] = b; fThresholdConstants[1][2] = c; }
  void SetBackgroundThreshold(Int_t t)                           { fBkgThreshold             = t; }
  void SetL0Threshold(Int_t t)                                   { fL0Threshold              = t; }
  /**
   * Set the size of jet patches
   * @param patchsize Size of jet patches
   */
  void SetJetPatchsize(Int_t patchsize)                          { fJetPatchsize = patchsize; }
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
  void SetTriggerBitConfig(const AliEMCALTriggerBitConfig *const config) { fTriggerBitConfig = config; }

  /**
   * Set range for L0 time
   * @param min Minimum L0 time (default is 7)
   * @param max Maximum L0 time (default is 10)
   */
  void SetL0TimeRange(Int_t min, Int_t max) { fL0MinTime = min; fL0MaxTime = max; }

  /**
   * Add an offline bad channel to the set
   * @param absId Absolute ID of the bad channel
   */
  void AddOfflineBadChannel(Short_t absId) { fOfflineBadChannels.insert(absId); }

  /**
   * Read the offline bad channel map from a standard stream
   * @param stream A reference to a standard stream to read from (can be a file stream)
   */
  void ReadOfflineBadChannelFromStream(std::istream& stream);

  /**
   * Read the offline bad channel map from a text file
   * @param fname Path and name of the file
   */
  void ReadOfflineBadChannelFromFile(const char* fname);

  /**
   * Reset data grids
   */
  void Reset();

  /**
   * Read the calo trigger data
   * @param trigger Input calo trigger data
   * @param timeMin Minimum L0 time
   * @param timeMax Maximum L0 time
   * @param applyTimeCut If true, only L0 amoplitudes in the range timeMin, timeMax are considered
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
   * @return the status of the patch (not L0, candidate, fired)
   */
  ELevel0TriggerStatus_t CheckForL0(Int_t col, Int_t row) const;

  /**
   * Create trigger algorithm for gamma triggers
   * @param rowmin Minimum row the trigger algorithm operates on
   * @param rowmax Maximum row the trigger algorithm operates on
   * @return The gamma trigger algorithm
   */
  AliEMCALTriggerAlgorithm<double> *CreateGammaTriggerAlgorithm(Int_t rowmin, Int_t rowmax) const;
  /**
   * Create trigger algorithm for jet triggers
   * @param rowmin Minimum row the trigger algorithm operates on
   * @param rowmax Maximum row the trigger algorithm operates on
   * @return The jet trigger algorithm
   */
  AliEMCALTriggerAlgorithm<double> *CreateJetTriggerAlgorithm(Int_t rowmin, Int_t rowmax) const;
  /**
   * Create trigger algorithm for jet triggers
   * @param rowmin Minimum row the trigger algorithm operates on
   * @param rowmax Maximum row the trigger algorithm operates on
   * @return The jet trigger algorithm
   */
  AliEMCALTriggerAlgorithm<double> *CreateBkgTriggerAlgorithm(Int_t rowmin, Int_t rowmax) const;

  /**
   * Check from the bitmask whether the patch is a gamma patch
   * @param patch Patch to check
   * @return True if patch is a gamma patch
   */
  Bool_t IsGammaPatch(const AliEMCALTriggerRawPatch &patch) const;
  /**
   * Check from the bitmask whether the patch is a jet patch
   * @param patch Patch to check
   * @return True if patch is a jet patch
   */
  Bool_t IsJetPatch(const AliEMCALTriggerRawPatch &patch) const;
  /**
   * Check from the bitmask whether the patch is a background patch
   * @param patch Patch to check
   * @return True if patch is a background patch
   */
  Bool_t IsBkgPatch(const AliEMCALTriggerRawPatch &patch) const;

  AliEMCALTriggerChannelContainer           fBadChannels;                 ///< Container of bad channels
  std::set<Short_t>                         fOfflineBadChannels;          ///< Abd ID of offline bad channels
  const AliEMCALTriggerBitConfig            *fTriggerBitConfig;           ///< Trigger bit configuration, aliroot-dependent
  const AliEMCALGeometry                    *fGeometry;                   //!<! Underlying EMCAL geometry

  AliEMCALTriggerDataGrid<double>           *fPatchAmplitudes;            //!<! TRU Amplitudes (for L0)
  AliEMCALTriggerDataGrid<double>           *fPatchADCSimple;             //!<! patch map for simple offline trigger
  AliEMCALTriggerDataGrid<double>           *fPatchADC;                   //!<! ADC values map
  AliEMCALTriggerDataGrid<char>             *fLevel0TimeMap;              //!<! Map needed to store the level0 times
  AliEMCALTriggerDataGrid<int>              *fTriggerBitMap;              //!<! Map of trigger bits

  AliEMCALTriggerPatchFinder<double>        *fPatchFinder;                //!<! The actual patch finder
  AliEMCALTriggerAlgorithm<double>          *fLevel0PatchFinder;          //!<! Patch finder for Level0 patches
  Int_t                                     fL0MinTime;                   ///< Minimum L0 time
  Int_t                                     fL0MaxTime;                   ///< Maximum L0 time
  Double_t                                  fADCtoGeV;                    //!<! Conversion factor from ADC to GeV

  Int_t                                     fJetPatchsize;                ///< Size of a jet patch
  Int_t                                     fThresholdConstants[4][3];    ///< simple offline trigger thresholds constants
  ULong64_t                                 fL1ThresholdsOffline[4];      ///< container for V0-dependent offline thresholds
  Int_t                                     fBkgThreshold;                ///< threshold for the background patches (8x8)
  Int_t                                     fL0Threshold;                 ///< threshold for the L0 patches (2x2)
  Bool_t                                    fIsMC;                        ///< Set MC offset

  Int_t                                     fDebugLevel;                  ///< Debug lebel;

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerMakerKernel, 1);
  /// \endcond
};

#endif
