/**
 * @file AliEMCALTriggerAlgorithm.h
 * @date Oct. 23, 2015
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#ifndef AliEMCALTRIGGERALGORITHM_H
#define AliEMCALTRIGGERALGORITHM_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliEMCALTriggerRawPatch.h"
#include <TObject.h>
#include <vector>

template<typename T> class AliEMCALTriggerDataGrid;

/**
 * @class AliEMCALTriggerAlgorithm
 * @brief Base class for EMCAL Level1 trigger algorithms
 */
template<typename T>
class AliEMCALTriggerAlgorithm : public TObject {
public:
  AliEMCALTriggerAlgorithm();
  AliEMCALTriggerAlgorithm(Int_t rowmin, Int_t rowmax, UInt_t bitmask);
  virtual ~AliEMCALTriggerAlgorithm();

  void SetRowMin(Int_t rowmin) { fRowMin = rowmin; }
  void SetRowMax(Int_t rowmax) { fRowMax = rowmax; }
  void SetThresholds(Float_t th, Float_t offTh) { fThreshold = th; fOfflineThreshold = offTh; }
  void SetBitMask(UInt_t bitmask) { fBitMask = bitmask; }
  void SetPatchSize(Int_t patchsize) { fPatchSize = patchsize; }
  void SetSubregionSize(Int_t subregionsize) { fSubregionSize = subregionsize; }

  virtual std::vector<AliEMCALTriggerRawPatch> FindPatches(const AliEMCALTriggerDataGrid<T> &adc, const AliEMCALTriggerDataGrid<T> &offlineAdc) const;

protected:
  Int_t                             fRowMin;
  Int_t                             fRowMax;
  Int_t                             fPatchSize;
  Int_t                             fSubregionSize;
  UInt_t                            fBitMask;
  Float_t                           fThreshold;
  Float_t                           fOfflineThreshold;

  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerAlgorithm, 2);
  /// \endcond
};

/**
 * @class AliEMCALJetTriggerAlgorithm
 * @brief Implementation of the EMCAL jet trigger algorithm
 *
 * A jet
 */
template<typename T>
class AliEMCALJetTriggerAlgorithm : public AliEMCALTriggerAlgorithm<T> {
public:
  /**
   * Constructor
   */
  AliEMCALJetTriggerAlgorithm();
  /**
   * Constructor, setting also range limits and bit mask
   * @param rowmin Min. row used for patch finding
   * @param rowmax Max. row used for patch finding
   * @param bitmask Bitmask stored in the raw patches
   */
  AliEMCALJetTriggerAlgorithm(Int_t rowmin, Int_t rowmax, UInt_t bitmask);
  /**
   * Destructor
   */
  virtual ~AliEMCALJetTriggerAlgorithm();

  /// \cond CLASSIMP
  ClassDef(AliEMCALJetTriggerAlgorithm, 1);
  /// \endcond
};

template<typename T>
class AliEMCALGammaTriggerAlgorithm : public AliEMCALTriggerAlgorithm<T> {
public:
  AliEMCALGammaTriggerAlgorithm();
  AliEMCALGammaTriggerAlgorithm(Int_t rowmin, Int_t rowmax, UInt_t Bitmask);
  virtual ~AliEMCALGammaTriggerAlgorithm();

  /// \cond CLASSIMP
  ClassDef(AliEMCALGammaTriggerAlgorithm, 1);
  /// \endcond
};

template<typename T>
class AliEMCALBkgTriggerAlgorithm : public AliEMCALTriggerAlgorithm<T> {
public:
  AliEMCALBkgTriggerAlgorithm();
  AliEMCALBkgTriggerAlgorithm(Int_t rowmin, Int_t rowmax, UInt_t Bitmask);
  virtual ~AliEMCALBkgTriggerAlgorithm();

  /// \cond CLASSIMP
  ClassDef(AliEMCALBkgTriggerAlgorithm, 1);
  /// \endcond
};


#endif
