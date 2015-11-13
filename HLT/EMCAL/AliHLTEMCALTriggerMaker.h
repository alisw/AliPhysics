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
  enum{
    kOriginSTU = 0,
    kOriginRECAL = 1,
    kOriginCELLS = 1
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
   * Set the origin of the trigger patches (STU, RECAL, CELLS)
   * @param origin Origin of the trigger patches
   */
  void SetOrigin(UChar_t origin) { fOrigin = origin; }

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
   * Set the pointer to the writeout buffer
   * @param outputcont Pointer to the writeout buffer
   */
  void SetTriggerPatchDataPtr(AliHLTCaloTriggerPatchDataStruct *outputcont) { fTriggerPatchDataPtr = outputcont; }

  /**
   * Find EMCAL trigger patches.
   * The trigger patches are internally converted into HLT EMCAL trigger patches and
   * pushed to the writeout buffer.
   * @return Number of found patches (-1 in case of buffer overflow)
   */
  Int_t FindPatches();

  /**
   * Set the threshold for jet trigger patches
   * @param threshold Threshold to be applied
   */
  void SetJetThreshold(Float_t threshold) { fJetThreshold = threshold; }

  /**
   *
   */
  void SetGammaThreshold(Float_t threshold) { fGammaThreshold = threshold; }

protected:
  /**
   * Convert raw patches found by the trigger patch finders into HLT EMCAL trigger patches.
   * @param inputpatch EMCAL raw patch to be converted into an EMCAL HLT patch
   * @param output HLT trigger patch obtaied using the information in the EMCAL raw patch
   */
   void MakeHLTPatch(const AliEMCALTriggerRawPatch &inputpatch, AliHLTCaloTriggerPatchDataStruct &output, UChar_t origin) const;

  /**
   * Initialise trigger patch finders in the EMCAL
   */
  void InitializeEMCALPatchFinders();
  /**
   * Initializse trigger patch finders in the DCAL
   */
  void InitializeDCALPatchFinders();

private:
  UChar_t                                       fOrigin;
  /** Pointer to the output buffer */
  AliHLTCaloTriggerPatchDataStruct              *fTriggerPatchDataPtr;
  /** Underlying EMCAL geometry*/
  const AliHLTEMCALGeometry                     *fkGeometryPtr;

  /** EMCAL trigger algorithms used to find trigger patches */
  AliEMCALTriggerPatchFinder<float>             *fPatchFinder;
  /** Grid with ADC values used for the trigger patch finding */
  AliEMCALTriggerDataGrid<float>                *fADCValues;

  /** offline threshold for gamma patches */
  Float_t                                       fGammaThreshold;
  /** offline threshold for jet patches */
  Float_t                                       fJetThreshold;

  ClassDef(AliHLTEMCALTriggerMaker, 1);
};

#endif
