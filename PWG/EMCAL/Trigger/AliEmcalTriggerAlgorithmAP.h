/**
 * @file AliEmcalTriggerAlgorithmAP.h
 * @date Oct. 23, 2015
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#ifndef AliEmcalTriggerAlgorithmAPAP_H
#define AliEmcalTriggerAlgorithmAPAP_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliEmcalTriggerRawPatchAP.h"
#include <TObject.h>
#include <vector>



template<typename T> class AliEmcalTriggerDataGridAP;

/**
 * @class AliEmcalTriggerAlgorithmAP
 * @brief Base class for EMCAL Level1 trigger algorithms
 */
template<typename T>
class AliEmcalTriggerAlgorithmAP : public TObject {
public:
  AliEmcalTriggerAlgorithmAP();
  AliEmcalTriggerAlgorithmAP(Int_t rowmin, Int_t rowmax, ULong_t bitmask);
  virtual ~AliEmcalTriggerAlgorithmAP();

  void SetRowMin(Int_t rowmin) { fRowMin = rowmin; }
  void SetRowMax(Int_t rowmax) { fRowMax = rowmax; }
  void SetThresholds(Float_t th, Float_t offTh) { fThreshold = th; fOfflineThreshold = offTh; }
  void SetBitMask(UInt_t bitmask) { fBitMask = bitmask; }
  void SetPatchSize(Int_t patchsize) { fPatchSize = patchsize; }

  virtual std::vector<AliEmcalTriggerRawPatchAP> FindPatches(const AliEmcalTriggerDataGridAP<T> &adc, const AliEmcalTriggerDataGridAP<T> &offlineAdc) const;

protected:
  int                               fRowMin;
  int                               fRowMax;
  int                               fPatchSize;
  ULong_t                           fBitMask;
  double                            fThreshold;
  double                            fOfflineThreshold;

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerAlgorithmAP, 1);
  /// \endcond
};

/**
 * @class AliEmcalJetTriggerAlgorithmAP
 * @brief Implementation of the EMCAL jet trigger algorithm
 *
 * A jet
 */
template<typename T>
class AliEmcalJetTriggerAlgorithmAP : public AliEmcalTriggerAlgorithmAP<T> {
public:
  /**
   * Constructor
   */
  AliEmcalJetTriggerAlgorithmAP();
  /**
   * Constructor, setting also range limits and bit mask
   * @param rowmin Min. row used for patch finding
   * @param rowmax Max. row used for patch finding
   * @param bitmask Bitmask stored in the raw patches
   */
  AliEmcalJetTriggerAlgorithmAP(Int_t rowmin, Int_t rowmax, ULong_t bitmask);
  /**
   * Destructor
   */
  virtual ~AliEmcalJetTriggerAlgorithmAP();

  /// \cond CLASSIMP
  ClassDef(AliEmcalJetTriggerAlgorithmAP, 1);
  /// \endcond
};

template<typename T>
class AliEmcalGammaTriggerAlgorithmAP : public AliEmcalTriggerAlgorithmAP<T> {
public:
  AliEmcalGammaTriggerAlgorithmAP();
  AliEmcalGammaTriggerAlgorithmAP(Int_t rowmin, Int_t rowmax, ULong_t Bitmask);
  virtual ~AliEmcalGammaTriggerAlgorithmAP();

  /// \cond CLASSIMP
  ClassDef(AliEmcalGammaTriggerAlgorithmAP, 1);
  /// \endcond
};

#endif
