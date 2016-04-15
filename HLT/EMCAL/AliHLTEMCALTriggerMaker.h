/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Markus Fasel                                          *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTEMCALTRIGGERMAKER_H
#define ALIHLTEMCALTRIGGERMAKER_H

/**
 * @file   AliHLTEMCALTriggerMaker.h
 * @author Markus Fasel
 * @date   Oct. 30, 2015
 * @brief  Trigger maker kernel used by the EMCAL HLT trigger maker component
 */

#include <TObject.h>
#include "AliHLTLogging.h"

struct AliHLTCaloTriggerPatchDataStruct;
class AliHLTEMCALGeometry;

class AliEMCALTriggerBitConfig;
class AliEMCALTriggerRawPatch;
template<typename T> class AliEMCALTriggerDataGrid;
template<typename T> class AliEMCALTriggerPatchFinder;

/**
 * @class AliHLTEMCALTriggerMaker
 * @brief EMCAL Trigger Maker implementation runnin on the HLT
 *
 * This class implements the trigger maker for the EMCAL in the
 * HLT.
 */
class AliHLTEMCALTriggerMaker : public TObject, public AliHLTLogging {
public:
  enum ELevel0TriggerStatus_t { kNotLevel0, kLevel0Candidate, kLevel0Fired };

  enum ThresholdType_t{
	  kHighThreshold = 0,
	  kLowThreshold = 1,
	  kNthresholds = 2
  };
  /**
   * Constructor
   */
  AliHLTEMCALTriggerMaker();

  /**
   * Destructor
   */
  virtual ~AliHLTEMCALTriggerMaker();

  /**
   * Initializes the trigger maker kernel. Sets the geometry ptr,
   * allocates space for the trigger channels and configures the
   * trigger patch finder.
   *
   * Based on the geometry versions the trigger channel data grid is
   * allocated and the trigger patch finders in EMCAL (and DCAL if available)
   * are configured.
   *
   * @param geo Pointer to the geometry used in the trigger maker.
   */
  void Initialise(const AliHLTEMCALGeometry *geo);

  /**
   * Rest ADC values
   */
  void ResetADC();

  /**
   * Add digit structure to the data grid
   * @param digit Input digit data
   */
  void AddDigit(const AliHLTCaloDigitDataStruct *digit);

  /**
   * Set the adc value for a given col / row combination
   * @param col Column
   * @param row Row
   * @param adc ADC value
   */
  void SetADC(Int_t col, Int_t row, Float_t adc);

  /**
   * Set the L0 amplitude for a given col / row combination
   * @param col Column
   * @param row Row
   * @param amp L0 Amplitude
   */
  void SetL0Amplitude(Int_t col, Int_t row, Float_t amp);

  /**
   * Set the L0 trigger time for a given col / row combination
   * @param col Column
   * @param row Row
   * @param amp L0 trigger time
   */
  void SetL0Time(Int_t col, Int_t row, UChar_t time);

    /**
   * Set the bit mask from the STU for a given col / row combination
   * @param col Column
   * @param row Row
   * @param bitMask Bit mask
   */
  void SetBitMask(Int_t col, Int_t row, Int_t bitMask);

  /**
   * Set the pointer to the writeout buffer
   * @param outputcont Pointer to the writeout buffer
   */
  void SetTriggerPatchDataPtr(AliHLTCaloTriggerPatchDataStruct *outputcont, AliHLTUInt32_t buffersize) { fTriggerPatchDataPtr = outputcont; fBufferSize = buffersize; }

  /**
   * Find EMCAL trigger patches.
   * The trigger patches are internally converted into HLT EMCAL trigger patches and
   * pushed to the writeout buffer.
   * @return Number of found patches (-1 in case of buffer overflow)
   */
  Int_t FindPatches();

  /**
   * Set the thresholds for jet trigger patches
   * @param onlineTh Online threshold to be applied
   * @param offlineTh Offline threshold to be applied
   */
  void SetJetThresholds(ThresholdType_t threshtype, Float_t onlineTh, Float_t offlineTh) { fJetThresholdOnline[threshtype] = onlineTh; fJetThresholdOffline[threshtype] = offlineTh; }
  
  /**
   * Set the thresholds for gamma trigger patches
   * @param onlineTh Online threshold to be applied
   * @param offlineTh Offline threshold to be applied
   */
  void SetGammaThresholds(ThresholdType_t threshtype, Float_t onlineTh, Float_t offlineTh) { fGammaThresholdOnline[threshtype] = onlineTh; fGammaThresholdOffline[threshtype] = offlineTh; }

  /**
   * Set the thresholds for bkg trigger patches
   * @param onlineTh Online threshold to be applied
   * @param offlineTh Offline threshold to be applied
   */
  void SetBkgThresholds(Float_t onlineTh, Float_t offlineTh) { fBkgThresholdOnline = onlineTh; fBkgThresholdOffline = offlineTh; }

  /**
   * Define whether we run the algorithm for background patches
   * @param doRun Switch for the background algorithm
   */
  void SetRunBkgAlgorithm(Bool_t doRun) { fRunBkgAlgorithm = doRun; }

  /**
   * Set the thresholds for level0 trigger patches
   * @param onlineTh Online Threshold to be applied
   * @param offlineTh Offline Threshold to be applied
   */
  void SetLevel0Thresholds(Float_t onlineTh, Float_t offlineTh) { fLevel0ThresholdOnline = onlineTh; fLevel0ThresholdOffline = offlineTh; }

  /**
   * Set the patch size/subregion for jet trigger patches
   * @param size Size of the patch in number of FastORs
   * @param subregion Subregion of the sliding window
   */
  void SetJetPatch(Int_t size, Int_t subregion)    { fJetPatchSize = size  ; fJetSubregionSize = subregion  ; }

  /**
   * Set the patch size/subregion for gamma trigger patches
   * @param size Size of the patch in number of FastORs
   * @param subregion Subregion of the sliding window
   */
  void SetGammaPatch(Int_t size, Int_t subregion)  { fGammaPatchSize = size; fGammaSubregionSize = subregion; }

  /**
   * Set the patch size/subregion for background trigger patches
   * @param size Size of the patch in number of FastORs
   * @param subregion Subregion of the sliding window
   */
  void SetBkgPatch(Int_t size, Int_t subregion)    { fBkgPatchSize = size  ; fBkgSubregionSize = subregion  ; }

protected:
  /**
   * Convert raw patches found by the trigger patch finders into HLT EMCAL trigger patches.
   * @param inputpatch EMCAL raw patch to be converted into an EMCAL HLT patch
   * @param output HLT trigger patch obtaied using the information in the EMCAL raw patch
   */
   void MakeHLTPatch(const AliEMCALTriggerRawPatch &inputpatch, AliHLTCaloTriggerPatchDataStruct &output, UInt_t offlinebits, UInt_t onlinebitmask, UInt_t level0bits) const;

  /**
   * Initialise trigger patch finders in the EMCAL or DCAL at Level1
   * @param isDCAL Switch distinguishing between DCAL (true) and EMCAL (false)
   */
  void InitializeLevel1PatchFinders(Bool_t isDCAL);

  /** Initialize L0 Patch finders in EMCAL and DCAL */
  void InitializeLevel0PatchFinders(Bool_t isDCAL);

  /**
   * Initialize the lookup tables used by the trigger maker
   */
  void InitializeLookupTables();

  /**
   * Check whether all fastors are in the same TRU. This
   * is a condition to accept the patch as valid Level0
   * patch.
   * @param patch Patch to check
   * @return True if all fastors are in the same TRU, false
   * otherwise
   */
  ELevel0TriggerStatus_t CheckForL0(Int_t col, Int_t row) const;

  Bool_t IsSameTRU(const AliEMCALTriggerRawPatch &patch) const;

private:
  /** Pointer to the output buffer */
  AliHLTCaloTriggerPatchDataStruct              *fTriggerPatchDataPtr;
  /** Underlying EMCAL geometry*/
  const AliHLTEMCALGeometry                     *fkGeometryPtr;

  /** EMCAL trigger algorithms used to find trigger patches */
  AliEMCALTriggerPatchFinder<float>             *fPatchFinder;
  /** EMCAL Level0 patch finder (operating on different data */
  AliEMCALTriggerPatchFinder<float>				*fL0PatchFinder;
  /** Grid with ADC values used for the trigger patch finding */
  AliEMCALTriggerDataGrid<float>                *fADCValues;
  /** Grid with ADC values used for the trigger patch finding */
  AliEMCALTriggerDataGrid<float>                *fADCOfflineValues;
  /** Grid with L0 Amplitude values used for L0 trigger patch finding */
  AliEMCALTriggerDataGrid<float>				*fL0Amplitudes;
  /** Grid with trigger bit mask from STU */
  AliEMCALTriggerDataGrid<int>                  *fTriggerBitMasks;
  /** Grid with L0 trigger time values used to retrieve L0 decision */
  AliEMCALTriggerDataGrid<unsigned char>        *fLevel0TimeMap;
  /** Lookup table with TRU indices */
  AliEMCALTriggerDataGrid<int>                  *fTRUIndexMap;
  /** Trigger bit configurtion */
  AliEMCALTriggerBitConfig                      *fTriggerBitConfig;
  /** Jet patch size **/
  Int_t                                         fJetPatchSize;
  /** Jet subregion size **/
  Int_t                                         fJetSubregionSize;
  /** Gamma patch size **/
  Int_t                                         fGammaPatchSize;
  /** Gamme subregion size **/
  Int_t                                         fGammaSubregionSize;
  /** Background patch size **/
  Int_t                                         fBkgPatchSize;
  /** Background subregion size **/
  Int_t                                         fBkgSubregionSize;
  /** Minimum time bin for a valid L0 trigger decision **/
  Char_t                                        fL0MinTime;
  /** Maximum time bin for a valid L0 trigger decision **/
  Char_t                                        fL0MaxTime;

  /** Available space in buffer */
  AliHLTUInt32_t                                fBufferSize;
  /** online threshold for gamma patches */
  Float_t                                       fGammaThresholdOnline[2];
  /** offline threshold for gamma patches */
  Float_t                                       fGammaThresholdOffline[2];
  /** online threshold for jet patches */
  Float_t                                       fJetThresholdOnline[2];
  /** offline threshold for jet patches */
  Float_t                                       fJetThresholdOffline[2];
  /** online threshold for bkg patches */
  Float_t                                       fBkgThresholdOnline;
  /** offline threshold for bkg patches */
  Float_t                                       fBkgThresholdOffline;
  /** Switch for running the background trigger algorithm */
  Bool_t                                        fRunBkgAlgorithm;
  /** Threshold for accepting level0 online patches */
  Float_t                                       fLevel0ThresholdOnline;
  /** Threshold for accepting level0 offline patches */
  Float_t                                       fLevel0ThresholdOffline;

  ClassDef(AliHLTEMCALTriggerMaker, 3);
};

#endif
