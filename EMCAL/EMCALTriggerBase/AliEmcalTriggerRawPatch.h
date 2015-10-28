#ifndef ALIEMCALTRIGGERRAWPATCH_H
#define ALIEMCALTRIGGERRAWPATCH_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <Riosfwd.h>
#include <TObject.h>

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
class AliEmcalTriggerRawPatch : public TObject {
public:
  AliEmcalTriggerRawPatch();
  AliEmcalTriggerRawPatch(Int_t col0, Int_t row0, Int_t size, Double_t adc);
  virtual ~AliEmcalTriggerRawPatch() {}

  bool operator==(const AliEmcalTriggerRawPatch &other) const;
  bool operator<(const AliEmcalTriggerRawPatch &other) const;

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

  void PrintStream(std::ostream &stream) const;

protected:
  ULong_t                       fBitMask;         ///< Trigger bit mask
  Int_t                         fCol0;            ///< Start column of the patch
  Int_t                         fRow0;            ///< Start row of the patch
  Int_t                         fSize;            ///< Patch size in number of FAST-ors
  Double_t                      fADC;             ///< Patch ADC

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerRawPatch, 1);
  /// \endcond
};

std::ostream &operator<<(std::ostream &in, const AliEmcalTriggerRawPatch &patch);

#endif
