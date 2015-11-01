/**
 * @file AliEmcalTriggerPatchFinderAP.h
 * @date Oct. 23, 2015
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#ifndef AliEmcalTriggerPatchFinderAPAP_H
#define AliEmcalTriggerPatchFinderAPAP_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliEmcalTriggerRawPatchAP.h"
#include <TObject.h>
#include <vector>


template<typename T> class AliEmcalTriggerDataGridAP;
template<typename T> class AliEmcalTriggerAlgorithmAP;

/**
 * @class AliEmcalTriggerPatchFinderAP
 * @brief Steering class for patch finder
 *
 * This class steers the trigger patch finding by calling
 * the trigger algorithms assigned to the patch finder. Trigger
 * algorithms are added via the function AddTriggerAlgorithm.
 * Patches are found via the function FindPatches.
 */
template<typename T>
class AliEmcalTriggerPatchFinderAP : public TObject {
public:
  /**
   * Constructor
   */
  AliEmcalTriggerPatchFinderAP();
  /**
   * Destructor
   */
  virtual ~AliEmcalTriggerPatchFinderAP();

  /**
   * Add trigger algorithm to the trigger patch finder
   * @param trigger Trigger algorithm assigned to the patch finder
   */
  void AddTriggerAlgorithm(AliEmcalTriggerAlgorithmAP<T> *trigger) { fTriggerAlgorithms.push_back(trigger); }

  /**
   * Find trigger patches usin the grid of adc values. All trigger patch finders are called one after each other.
   * The result contains the vector of patches from all trigger algorithms. The trigger patches are sorted in energy.
   * @param adc Data grid with ADC values
   * @return List of trigger patches found by all trigger algorithms assigned to this trigger patch finder.
   */
  std::vector<AliEmcalTriggerRawPatchAP> FindPatches(const AliEmcalTriggerDataGridAP<T> &adc) const;

protected:
  std::vector<AliEmcalTriggerAlgorithmAP<T> *>                 fTriggerAlgorithms;      ///< Trigger algoritms to be used

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerPatchFinderAP, 1);
  /// \endcond
};

#endif
