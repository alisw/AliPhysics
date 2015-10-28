/**
 * @file AliEmcalTriggerAlgorithm.h
 * @date Oct. 23, 2015
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#ifndef ALIEMCALTRIGGERALGORITHM_H
#define ALIEMCALTRIGGERALGORITHM_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <vector>

#include "AliEmcalTriggerRawPatch.h"


template<typename T> class AliEmcalTriggerDataGrid;

/**
 * @class AliEmcalTriggerAlgorithm
 * @brief Base class for EMCAL Level1 trigger algorithms
 */
template<typename T>
class AliEmcalTriggerAlgorithm : public TObject {
public:
  AliEmcalTriggerAlgorithm();
  AliEmcalTriggerAlgorithm(Int_t rowmin, Int_t rowmax, ULong_t bitmask);
  virtual ~AliEmcalTriggerAlgorithm();

  void SetRowMin(Int_t rowmin);
  void SetRowMax(Int_t rowmax);
  void SetThreshold(Double_t threshold);
  void SetBitMask(ULong_t bitmask);
  void SetPatchSize(Int_t patchsize);

  virtual std::vector<AliEmcalTriggerRawPatch> FindPatches(const AliEmcalTriggerDataGrid<T> &adc) const;

protected:
  int                               fRowMin;
  int                               fRowMax;
  int                               fPatchSize;
  ULong_t                           fBitMask;
  double                            fThreshold;

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerAlgorithm, 1);
  /// \endcond
};

/**
 * @class AliEmcalJetTriggerAlgorithm
 * @brief Implementation of the EMCAL jet trigger algorithm
 *
 * A jet
 */
template<typename T>
class AliEmcalJetTriggerAlgorithm : public AliEmcalTriggerAlgorithm<T> {
public:
  /**
   * Constructor
   */
  AliEmcalJetTriggerAlgorithm();
  /**
   * Constructor, setting also range limits and bit mask
   * @param rowmin Min. row used for patch finding
   * @param rowmax Max. row used for patch finding
   * @param bitmask Bitmask stored in the raw patches
   */
  AliEmcalJetTriggerAlgorithm(Int_t rowmin, Int_t rowmax, ULong_t bitmask);
  /**
   * Destructor
   */
  virtual ~AliEmcalJetTriggerAlgorithm();

  /// \cond CLASSIMP
  ClassDef(AliEmcalJetTriggerAlgorithm, 1);
  /// \endcond
};

template<typename T>
class AliEmcalGammaTriggerAlgorithm : public AliEmcalTriggerAlgorithm<T> {
public:
  AliEmcalGammaTriggerAlgorithm();
  AliEmcalGammaTriggerAlgorithm(Int_t rowmin, Int_t rowmax, ULong_t Bitmask);
  virtual ~AliEmcalGammaTriggerAlgorithm();

  /// \cond CLASSIMP
  ClassDef(AliEmcalGammaTriggerAlgorithm, 1);
  /// \endcond
};

#endif
