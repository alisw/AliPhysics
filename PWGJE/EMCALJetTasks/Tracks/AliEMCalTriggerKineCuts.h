/**
 * \file AliEMCalTriggerKineCuts.h
 * \brief
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \since Dec 12, 2014
 */
#ifndef ALIEMCALTRIGGERKINECUTS_H
#define ALIEMCALTRIGGERKINECUTS_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include "AliCutValueRange.h"

class AliVParticle;

namespace PWGJE {
  
namespace EMCALJetTasks {

/**
 * \class AliEMCalTriggerKineCuts
 * \brief Basic kinematic cuts for single track selection
 *
 * Class implementing basic kinematic track selection. Following selection
 * cuts are implemented:
 * - \f$ p_{t} \$f
 * - \f$ \eta \f$
 * - \f$ \phi \f$
 * For each cut variable a range (min and max) is expected.
 */
class AliEMCalTriggerKineCuts: public TObject {
public:
  AliEMCalTriggerKineCuts();
  /**
   * Destructor
   */
  virtual ~AliEMCalTriggerKineCuts() {}

  /**
   *
   * @param ptmin
   * @param ptmax
   */
  void SetPtRange(double ptmin, double ptmax) { fPtCut.SetLimits(ptmin, ptmax); }
  /**
   *
   * @param etamin
   * @param etamax
   */
  void SetEtaRange(double etamin, double etamax) { fEtaCut.SetLimits(etamin, etamax); }
  /**
   *
   * @param phimin
   * @param phimax
   */
  void SetPhiRange(double phimin, double phimax) { fPhiCut.SetLimits(phimin, phimax); }

  bool IsSelected(const AliVParticle *const track) const;

protected:
  AliCutValueRange<double>     fPtCut;            ///< Cut range in \f$ p_{t} \f$
  AliCutValueRange<double>     fEtaCut;           ///< Cut range in \f$ \eta \f$
  AliCutValueRange<double>     fPhiCut;           ///< Cut range in \f$ \phi \f$

  ClassDef(AliEMCalTriggerKineCuts, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERKINECUTS_H */
