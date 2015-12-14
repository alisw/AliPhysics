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
   void MakeHLTPatch(const AliEMCALTriggerRawPatch &inputpatch, AliHLTCaloTriggerPatchDataStruct &output, UInt_t offlinebits) const;

  /**
   * Initialise trigger patch finders in the EMCAL
   */
  void InitializeEMCALPatchFinders();
  /**
   * Initializse trigger patch finders in the DCAL
   */
  void InitializeDCALPatchFinders();

private:
  /** Pointer to the output buffer */
  AliHLTCaloTriggerPatchDataStruct              *fTriggerPatchDataPtr;
  /** Underlying EMCAL geometry*/
  const AliHLTEMCALGeometry                     *fkGeometryPtr;

  /** EMCAL trigger algorithms used to find trigger patches */
  AliEMCALTriggerPatchFinder<float>             *fPatchFinder;
  /** Grid with ADC values used for the trigger patch finding */
  AliEMCALTriggerDataGrid<float>                *fADCValues;
  /** Grid with ADC values used for the trigger patch finding */
  AliEMCALTriggerDataGrid<float>                *fADCOfflineValues;
  /** Grid with trigger bit mask from STU */
  AliEMCALTriggerDataGrid<int>                  *fTriggerBitMasks;
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

  ClassDef(AliHLTEMCALTriggerMaker, 2);
};

#endif
