#ifndef ALIEMCALTRIGGERMAKERKERNEL_H
#define ALIEMCALTRIGGERMAKERKERNEL_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <set>
#include <iostream>

#include <TObject.h>
#include <TArrayF.h>

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
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @author Salvatore Aiola, Yale University
 * @since Oct. 22nd, 2015
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
   * @param[in] inputevent Input ESD/AOD event, used for kinematics calculation
   * @param[in] useL0amp if true the Level0 amplitude is used
   * @return Array of reconstructed trigger patches
   */
  TObjArray *CreateTriggerPatches(const AliVEvent *inputevent, Bool_t useL0amp=kFALSE);

  /**
   * Get the list of online masked FastOR's used in the trigger maker
   * @return Absolute FastOR IDs (full EMCAL + DCAL in run2+) of masked channels
   */
  const std::set<Short_t> &GetListOfBadFastORAbsIDs() const { return fBadChannels; }

  /**
   * Get the list of bad cells used when calculating offline patch energy. These
   * cells were not masked at online level, instead they were determined by a calibration
   * algorithm.
   * @return Absolute Cell IDs (full EMCAL + DCAL in run2+) of bad cells
   */
  const std::set<Short_t> &GetListOfOfflineBadCellAbsIDs() const { return fOfflineBadChannels; }

  void SetTriggerThresholdJetLow   ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[2][0] = a; fThresholdConstants[2][1] = b; fThresholdConstants[2][2] = c; }
  void SetTriggerThresholdJetHigh  ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[0][0] = a; fThresholdConstants[0][1] = b; fThresholdConstants[0][2] = c; }
  void SetTriggerThresholdGammaLow ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[3][0] = a; fThresholdConstants[3][1] = b; fThresholdConstants[3][2] = c; }
  void SetTriggerThresholdGammaHigh( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[1][0] = a; fThresholdConstants[1][1] = b; fThresholdConstants[1][2] = c; }
  void SetBackgroundThreshold(Int_t t)                           { fBkgThreshold             = t; }
  void SetL0Threshold(Int_t t)                                   { fL0Threshold              = t; }

  /**
   * Define whether running on MC or not (for offset)
   * @param isMC Flag for MC
   */
  void SetIsMC(Bool_t isMC) { fIsMC = isMC; }

  /**
   * Provide the EMCAL geometry to the trigger maker Kernel
   * @param[in] geo EMCAL geometry
   */
  void SetGeometry(const AliEMCALGeometry *const geo) { fGeometry = geo; }

  /**
   * Set the trigger bit configuration applied for the given data set
   * @param[in] config Trigger bit config applied
   */
  void SetTriggerBitConfig(const AliEMCALTriggerBitConfig *const config);

  /**
   * Set range for L0 time
   * @param min Minimum L0 time (default is 7)
   * @param max Maximum L0 time (default is 10)
   */
  void SetL0TimeRange(Int_t min, Int_t max) { fL0MinTime = min; fL0MaxTime = max; }

  /**
   * Set thresholds applied to FastORs and offline cells before patch reconstruction
   * @param[in] l0 Threshold for L0 FastOR amplitudes
   * @param[in] l1 Threshold for L1 FastOR amplitudes
   * @param[in] cell Threshold for cell amplitudes
   */
  void SetFastORandCellThresholds(Int_t l0, Int_t l1, Int_t cell) { fMinL0FastORAmp = l0; fMinL1FastORAmp = l1; fMinCellAmp = cell; }

  /**
   * Add a FastOR bad channel to the list
   * @param[in] absId Absolute ID of the bad channel
   */
  void AddFastORBadChannel(Short_t absId) { fBadChannels.insert(absId); }

  /**
   * Read the FastOR bad channel map from a standard stream
   * @param[in] stream A reference to a standard stream to read from (can be a file stream)
   */
  void ReadFastORBadChannelFromStream(std::istream& stream);

  /**
   * Clear FastOR bad channel list
   */
  void ClearFastORBadChannels();

  /**
   * Read the FastOR bad channel map from a text file
   * @param[in] fname Path and name of the file
   */
  void ReadFastORBadChannelFromFile(const char* fname);

  /**
   * Add an offline bad channel to the set
   * @param[in] absId Absolute ID of the bad channel
   */
  void AddOfflineBadChannel(Short_t absId) { fOfflineBadChannels.insert(absId); }

  /**
   * Read the offline bad channel map from a standard stream
   * @param[in] stream A reference to a standard stream to read from (can be a file stream)
   */
  void ReadOfflineBadChannelFromStream(std::istream& stream);

  /**
   * Read the offline bad channel map from a text file
   * @param[in] fname Path and name of the file
   */
  void ReadOfflineBadChannelFromFile(const char* fname);

  /**
   * Clear offline bad channel list.
   */
  void ClearOfflineBadChannels();

  /**
   * Read the FastOR pedestals from a standard stream
   * @param[in] stream A reference to a standard stream to read from (can be a file stream)
   */
  void ReadFastORPedestalFromStream(std::istream& stream);

  /**
   * Read the FastOR pedestals from a text file
   * @param[in] fname Path and name of the file
   */
  void ReadFastORPedestalFromFile(const char* fname);

  /**
   * Set the pedestal value for a FastOR
   * @param[in] absId Absolute ID of a FastOR
   * @param[in] ped   Pedestal value
   */
  void SetFastORPedestal(Short_t absId, Float_t ped);

  /**
   * Reset the FastOR pedestal array
   */
  void ResetFastORPedestal() { fFastORPedestal.Reset(); }

  /**
   * Set symmmetric limit of the cell time allowed to accept the cell contributing to cell offline energy / offline patches.
   * @param[in] maxCellTime Maximum accepted value of the absolute cell time
   */
  void SetMaxAbsCellTime(Double_t maxCellTime) { fCellTimeLimits[0] = -maxCellTime; fCellTimeLimits[1] = maxCellTime; }

  /**
   * Set symmmetric limit of the cell time allowed to accept the cell contributing to cell offline energy / offline patches.
   * @param[in] minCellTime Minimum allowed cell time
   * @param[in] maxCellTime Maximum allowed cell time
   */
  void SetAcceptCellTimeRange(Double_t minCellTime, Double_t maxCellTime) { fCellTimeLimits[0] = minCellTime; fCellTimeLimits[1] = maxCellTime; }

  /**
   * Set lower cut on the cell time allowed to accept the cell contributing to cell offline energy / offline patches.
   * @param[in] minCellTime Minimum allowed cell time
   */
  void SetLowerCellTimeCut(Double_t minCellTime) { fCellTimeLimits[0] = minCellTime; }

  /**
   * Set upper cut on the cell time allowed to accept the cell contributing to cell offline energy / offline patches.
   * @param[in] minCellTime Minimum allowed cell time
   */
  void SetUpperCellTimeCut(Double_t maxCellTime) { fCellTimeLimits[1] = maxCellTime; }

  /**
   * Set cut on the minimum cell amplitude
   * @param[in] minamp Min cell amplitude accepted
   */
  void SetMinFEEAmplitude(Double_t minamp) { fMinCellAmplitude = minamp; }

  /**
   * Apply online bad channel masking to offline channel energies. This means that
   * cell energies in cells within FastOrs that are masked online are ignored. By
   * applying this the online trigger acceptance can be applied to offline patches
   * as well.
   * @param[in] doApply If true the online masking is applied to offline patch energies
   */
  void SetApplyOnlineBadChannelMaskingToOffline(Bool_t doApply = kTRUE) { fApplyOnlineBadChannelsToOffline = doApply; }

  /**
   * Reset data grids
   */
  void Reset();

  /**
   * Read the calo trigger data
   * @param[in] trigger Input calo trigger data
   */
  void ReadTriggerData(AliVCaloTrigger *trigger);

  /**
   * Read the EMCAL cell data
   * @param[in] cells EMCAL cell data
   */
  void ReadCellData(AliVCaloCells *cells);

  /**
   * Build VZERO-dependent thresholds for the offline trigger
   * @param vzdata VERO charges
   */
  void BuildL1ThresholdsOffline(const AliVVZERO *vzdata);

  /**
   * Add a L1 trigger algorithm
   * @param[in] rowmin Minimum row value
   * @param[in] rowmax Maximum row value
   * @param[in] bitmask Offline bit mask to be applied to the patches
   * @param[in] patchSize Size of the patches
   * @param[in] subregionSize Size of the sliding sub region
   */
  void AddL1TriggerAlgorithm(Int_t rowmin, Int_t rowmax, UInt_t bitmask, Int_t patchSize, Int_t subregionSize);

  /**
   * Set the L0 algorithm
   * @param[in] rowmin Minimum row value
   * @param[in] rowmax Maximum row value
   * @param[in] bitmask Offline bit mask to be applied to the patches
   * @param[in] patchSize Size of the patches
   * @param[in] subregionSize Size of the sliding sub region
   */
  void SetL0TriggerAlgorithm(Int_t rowmin, Int_t rowmax, UInt_t bitmask, Int_t patchSize, Int_t subregionSize);

  /**
   * Check whether the trigger maker has been specially configured. Status has to
   * be set in the functions ConfigureForXX.
   * @return True if the trigger maker kernel is configured, false otherwise
   */
  Bool_t IsConfigured() const { return fConfigured; }

  /**
   * Configure the class for 2015 PbPb
   */
  void ConfigureForPbPb2015();

  /**
   * Configure the class for 2015 pp
   */
  void ConfigureForPP2015();

  /**
   * Configure the class for 2013 pPb
   */
  void ConfigureForPPb2013();

  /**
   * Configure the class for 2012 pp
   */
  void ConfigureForPP2012();

  /**
   * Configure the class for 2011 PbPb
   */
  void ConfigureForPbPb2011();

  /**
   * Configure the class for 2011 pp
   */
  void ConfigureForPP2011();

protected:
  enum{
    kColsEta = 48
  };

  /**
   * Accept trigger patch as Level0 patch. Level0 patches are identified as 2x2 FASTOR patches
   * in the same TRU
   * @param[in] col Starting coloumn of the patch
   * @param[in] row Starting row of the patch
   * @return the status of the patch (not L0, candidate, fired)
   */
  ELevel0TriggerStatus_t CheckForL0(Int_t col, Int_t row) const;

  /**
   * Check from the bitmask whether the patch is a gamma patch
   * @param[in] patch Patch to check
   * @return True if patch is a gamma patch
   */
  Bool_t IsGammaPatch(const AliEMCALTriggerRawPatch &patch) const;
  /**
   * Check from the bitmask whether the patch is a jet patch
   * @param[in] patch Patch to check
   * @return True if patch is a jet patch
   */
  Bool_t IsJetPatch(const AliEMCALTriggerRawPatch &patch) const;
  /**
   * Check from the bitmask whether the patch is a background patch
   * @param[in] patch Patch to check
   * @return True if patch is a background patch
   */
  Bool_t IsBkgPatch(const AliEMCALTriggerRawPatch &patch) const;

  /**
   * Check according to geometrical cuts whether the patch has
   * overlap with the PHOS region. Partial overlap is sufficient
   * for this check.
   * @param[in] patch Patch to check
   * @return True if the patch has at least partial overlap with PHOS
   */
  bool HasPHOSOverlap(const AliEMCALTriggerRawPatch &patch) const;

  std::set<Short_t>                         fBadChannels;                 ///< Container of bad channels
  std::set<Short_t>                         fOfflineBadChannels;          ///< Abd ID of offline bad channels
  TArrayF                                   fFastORPedestal;              ///< FastOR pedestal
  const AliEMCALTriggerBitConfig           *fTriggerBitConfig;            ///< Trigger bit configuration, aliroot-dependent

  AliEMCALTriggerPatchFinder<double>       *fPatchFinder;                 ///< The actual patch finder
  AliEMCALTriggerAlgorithm<double>         *fLevel0PatchFinder;           ///< Patch finder for Level0 patches
  Int_t                                     fL0MinTime;                   ///< Minimum L0 time
  Int_t                                     fL0MaxTime;                   ///< Maximum L0 time
  Int_t                                     fMinCellAmp;                  ///< Minimum offline amplitude of the cells used to generate the patches
  Int_t                                     fMinL0FastORAmp;              ///< Minimum L0 amplitude of the FastORs used to generate the patches
  Int_t                                     fMinL1FastORAmp;              ///< Minimum L1 amplitude of the FastORs used to generate the patches
  Int_t                                     fThresholdConstants[4][3];    ///< simple offline trigger thresholds constants
  ULong64_t                                 fL1ThresholdsOffline[4];      ///< container for V0-dependent offline thresholds
  Int_t                                     fBkgThreshold;                ///< threshold for the background patches (8x8)
  Int_t                                     fL0Threshold;                 ///< threshold for the L0 patches (2x2)
  Bool_t                                    fIsMC;                        ///< Set MC offset
  Int_t                                     fDebugLevel;                  ///< Debug lebel;
  Double_t                                  fCellTimeLimits[2];           ///< Maximum allowed abs cell time (default [0] = - 10000, [1] = 10000)
  Double_t                                  fMinCellAmplitude;            ///< Minimum amplitude in cell required to be considered for filling the data grid
  Bool_t                                    fApplyOnlineBadChannelsToOffline;   ///< Apply online bad channels to offline ADC values
  Bool_t                                    fConfigured;                  ///< Switch specifying whether the trigger maker kernel has been configured for a given data set

  const AliEMCALGeometry                    *fGeometry;                   //!<! Underlying EMCAL geometry
  AliEMCALTriggerDataGrid<double>           *fPatchAmplitudes;            //!<! TRU Amplitudes (for L0)
  AliEMCALTriggerDataGrid<double>           *fPatchADCSimple;             //!<! patch map for simple offline trigger
  AliEMCALTriggerDataGrid<double>           *fPatchADC;                   //!<! ADC values map
  AliEMCALTriggerDataGrid<char>             *fLevel0TimeMap;              //!<! Map needed to store the level0 times
  AliEMCALTriggerDataGrid<int>              *fTriggerBitMap;              //!<! Map of trigger bits

  Double_t                                  fADCtoGeV;                    //!<! Conversion factor from ADC to GeV

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerMakerKernel, 4);
  /// \endcond
};

#endif
