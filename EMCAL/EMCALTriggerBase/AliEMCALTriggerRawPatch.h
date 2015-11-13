#ifndef ALIEMCALTRIGGERRAWPATCH_H
#define ALIEMCALTRIGGERRAWPATCH_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <Riosfwd.h>
#include <TObject.h>
#include "AliEMCALTriggerConstants.h"

/**
 * @class AliEmcalTriggerRawPatch
 * @brief Raw patch information used inside the trigger maker kernel
 * for the offline trigger and trigger recalculator
 *
 * Within the EMCAL trigger maker patches can be found by the offline trigger
 * or by the trigger patch recalculator. The trigger raw patch is supposed to
 * keep a minimum information needed with trigger patch information in order
 * to calculated the more detailed AliEmcalTriggerPatchInfo from it.
 */
class AliEMCALTriggerRawPatch : public TObject {
public:
  /**
   * Dummy constructor
   */
  AliEMCALTriggerRawPatch();

  /**
   * Main constructor
   * @param col0 Starting column
   * @param row0 Starting row
   * @param size Patch size
   * @param adc ADC value
   */
  AliEMCALTriggerRawPatch(Int_t col0, Int_t row0, Int_t size, Double_t adc, Double_t offlineAdc);

  /**
   * Destructor
   */
  virtual ~AliEMCALTriggerRawPatch() {}

  /**
   * Comparison operator for equalness: Patches are equal if they have the same position and
   * the same trigger bit mask.
   * @param other Patch to compare to
   * @return True if the patches share the same position and trigger bit mask, false otherwise
   */
  bool operator==(const AliEMCALTriggerRawPatch &other) const;

  /**
   * Comparison operator for smaller. As this is used in sorting algorithms, the comparison
   * is made based on the patch ADC.
   * @param other Patch to compate to
   * @return True if the patch ADC of this patch is smaller, false otherwise
   */
  bool operator<(const AliEMCALTriggerRawPatch &other) const;

  /**
   * Set the starting column of the patch
   * @param col0 Starting column of the patch
   */
  void SetColStart(Int_t col0) { fCol0 = col0; }

  /**
   * Set the starting row of the patch
   * @param row0 Starting row of the patch
   */
  void SetRowStart(Int_t row0) { fRow0 = row0; }

  /**
   * Set the patch size
   * @param patchsize Patch size
   */
  void SetPatchSize(Int_t patchsize) { fSize = patchsize; }

  /**
   * Set the patch ADC
   * @param adc Patch ADC
   */
  void SetADC(Double_t adc) { fADC = adc; }

  /**
   * Set the patch offline ADC
   * @param adc Patch offline ADC
   */
  void SetOfflineADC(Double_t adc) { fOfflineADC = adc; }

  /**
   * Set the patch trigger bit mask
   * @param bitmask Patch trigger bit mask
   */
  void SetBitmask(ULong_t bitmask) { fBitMask = bitmask; }

  /**
   * Get the starting column of the patch
   * @return Starting column of the patch
   */
  Int_t GetColStart() const { return fCol0; }

  /**
   * Get the starting row of the patch
   * @return Starting row of the patch
   */
  Int_t GetRowStart() const { return fRow0; }

  /**
   * Get the size of then patch (in number of FAST-ors per direction)
   * @return Patch size
   */
  Int_t GetPatchSize() const { return fSize; }

  /**
   * Get the patch ADC
   * @return patch ADC
   */
  Double_t GetADC() const { return fADC; }

  /**
   * Get the patch offline ADC
   * @return patch offline ADC
   */
  Double_t GetOfflineADC() const { return fOfflineADC; }

  /**
   * Get the patch trigger bit mask
   * @return Patch trigger bit mask
   */
  ULong_t GetBitmask() const { return fBitMask; }

  /**
   * Print trigger patch information to a stream
   * @param stream Output stream
   */
  void PrintStream(std::ostream &stream) const;

protected:
  ULong_t                       fBitMask;         ///< Trigger bit mask
  Int_t                         fCol0;            ///< Start column of the patch
  Int_t                         fRow0;            ///< Start row of the patch
  Int_t                         fSize;            ///< Patch size in number of FAST-ors
  Double_t                      fADC;             ///< Patch ADC
  Double_t                      fOfflineADC;      ///< Patch ADC

  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerRawPatch, 2);
  /// \endcond
};

/**
 * output stream operator
 */
std::ostream &operator<<(std::ostream &in, const AliEMCALTriggerRawPatch &patch);

#endif
