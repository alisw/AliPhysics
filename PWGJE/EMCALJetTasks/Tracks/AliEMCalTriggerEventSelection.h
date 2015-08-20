#ifndef ALIEMCALTRIGGEREVENTSELECTION_H
#define ALIEMCALTRIGGEREVENTSELECTION_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#include <TObject.h>
#include "AliCutValueRange.h"

namespace EMCalTriggerPtAnalysis {

/**
 * \namespace EMCalTriggerPtAnalysis
 * \brief Analysis of high-\f$ p_{t} \f$ tracks in triggered events
 *
 * This namespace contains classes for the analysis of high-\f$ p_{t} \f$ tracks in
 * triggered events.
 */
class AliEMCalTriggerEventData;

class AliEMCalTriggerEventSelection: public TObject {
public:
  AliEMCalTriggerEventSelection();
  AliEMCalTriggerEventSelection(const AliEMCalTriggerEventSelection &ref);
  AliEMCalTriggerEventSelection &operator=(const AliEMCalTriggerEventSelection &ref);
  virtual ~AliEMCalTriggerEventSelection() {}

  void SetVertexCut(double zmin, double zmax) { fVertexCut.SetLimits(zmin, zmax); }

  virtual bool IsEventSelected(const AliEMCalTriggerEventData * const ev) const;

protected:
  AliCutValueRange<double>   fVertexCut;

  ClassDef(AliEMCalTriggerEventSelection, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGEREVENTSELECTION_H */
