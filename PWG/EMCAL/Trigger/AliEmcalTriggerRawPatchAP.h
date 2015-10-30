/**
 * @file AliEmcalTriggerRawPatchAP.h
 * @since Oct 23, 2015
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#ifndef AliEmcalTriggerRawPatchAPAP_H
#define AliEmcalTriggerRawPatchAPAP_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <Riosfwd.h>
#include <TObject.h>

/**
 * @class AliEmcalTriggerRawPatchAP
 * @brief Raw patch information used inside the trigger maker kernel
 * for the offline trigger and trigger recalculator
 *
 * Within the EMCAL trigger maker patches can be found by the offline trigger
 * or by the trigger patch recalculator. The trigger raw patch is supposed to
 * keep a minimum information needed with trigger patch information in order
 * to calculated the more detailed AliEmcalTriggerPatchInfo from it.
 */
class AliEmcalTriggerRawPatchAP : public TObject {
public:
  /**
   * Dummy constructor
   */
  AliEmcalTriggerRawPatchAP();

  /**
   * Main constructor
   * @param col0 Starting column
   * @param row0 Starting row
   * @param size Patch size
   * @param adc ADC value
   */
  AliEmcalTriggerRawPatchAP(Int_t col0, Int_t row0, Int_t size, Double_t adc);

  /**
   * Destructor
   */
  virtual ~AliEmcalTriggerRawPatchAP() {}

  /**
   * Comparison operator for equalness: Patches are equal if they have the same position and
   * the same trigger bit mask.
   * @param other Patch to compare to
   * @return True if the patches share the same position and trigger bit mask, false otherwise
   */
  bool operator==(const AliEmcalTriggerRawPatchAP &other) const;

  /**
   * Comparison operator for smaller. As this is used in sorting algorithms, the comparison
   * is made based on the patch ADC.
   * @param other Patch to compate to
   * @return True if the patch ADC of this patch is smaller, false otherwise
   */
  bool operator<(const AliEmcalTriggerRawPatchAP &other) const;

  void SetColStart(Int_t col0) { fCol0 = col0; }
  void SetRowStart(Int_t row0) { fRow0 = row0; }
  void SetPatchSize(Int_t patchsize) { fSize = patchsize; }
  void SetADC(Double_t adc) { fADC = adc; }
  void SetBitmask(ULong_t bitmask) { fBitMask = bitmask; }

  Int_t GetColStart() const { return fCol0; }
  Int_t GetRowStart() const { return fRow0; }
  Int_t GetPatchSize() const { return fSize; }
  Int_t GetADC() const { return fADC; }
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

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerRawPatchAP, 1);
  /// \endcond
};

/**
 * output stream operator
 */
std::ostream &operator<<(std::ostream &in, const AliEmcalTriggerRawPatchAP &patch);

#endif
