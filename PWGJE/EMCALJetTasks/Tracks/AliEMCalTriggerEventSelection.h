/**
 * \file AliEMCalTriggerEventSelection.h
 * \brief Class for event selection (apart from trigger selection)
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#ifndef ALIEMCALTRIGGEREVENTSELECTION_H
#define ALIEMCALTRIGGEREVENTSELECTION_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include "AliCutValueRange.h"

class AliVEvent;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliEMCalTriggerEventData;

/**
 * \class AliEMCalTriggerEventSelection
 * \brief Class for event selection in the analysis of triggered events
 *
 * Basic event selection component: Selects events according to the pA cut and a vertex-z cut
 * For more sophisticated event selection the method IsEventSelected has to be overwritten
 */
class AliEMCalTriggerEventSelection: public TObject {
public:
  AliEMCalTriggerEventSelection();
  /**
   * Destructor
   */
  virtual ~AliEMCalTriggerEventSelection() {}

  void SetVertexCut(double zmin, double zmax) { fVertexCut.SetLimits(zmin, zmax); }
  void SetOldPileupSelection(bool doOld = true) { fOldPileupSelection = doOld; }
  void SetOldVertexSelection(bool doOld = true) { fOldVertexSelection = doOld; }

  virtual bool IsEventSelected(const AliEMCalTriggerEventData * const ev) const;

protected:
  Bool_t FalseVertexSelectionPA2013(const AliVEvent *const ev) const;
  AliCutValueRange<double>    fVertexCut;                         ///< cut range for the vertex selection
  Bool_t                      fOldPileupSelection;                ///< apply old pileup selection
  Bool_t                      fOldVertexSelection;                ///< apply old vertex selection

  ClassDef(AliEMCalTriggerEventSelection, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGEREVENTSELECTION_H */
