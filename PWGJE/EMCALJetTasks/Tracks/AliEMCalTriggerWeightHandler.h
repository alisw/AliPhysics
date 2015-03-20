/**
 * \file AliEMCalTriggerWeightHandler.h
 * \brief Weight handler for the analysis of high-\f$ p_{t} \f$ tracks in EMCAL-triggered events
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Mar 15, 2015
 */
#ifndef ALIEMCALTRIGGERWEIGHTHANDLER_H
#define ALIEMCALTRIGGERWEIGHTHANDLER_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class TF1;

class AliMCEvent;

/**
 * \namespace EMCalTriggerPtAnalysis
 * \brief Analysis of high-\f$ p_{t} \f$ tracks in triggered events
 *
 * This namespace contains classes for the analysis of high-\f$ p_{t} \f$ tracks in
 * triggered events.
 */
namespace EMCalTriggerPtAnalysis {

/**
 * \class AliEMCalTriggerWeightHandler
 * \brief Weight handler
 *
 * Weight handler, assigning an event-dependent weight. The weight is coming from a weight model,
 * which is based on an analytic description. For the moment it is assumed that the event depends
 * on the \f$ p_{t} \f$ of the hard interaction.
 */
class AliEMCalTriggerWeightHandler : public TObject {
public:
  AliEMCalTriggerWeightHandler();
  virtual ~AliEMCalTriggerWeightHandler() {}

  /**
   * Defines whether we use the \f$ p_{t} \f$ hard value of an event.
   * \param usePtHard
   */
  void SetUsePtHard(bool usePtHard) { fUsePtHard = usePtHard; }

  /**
   * Set the weight model
   * \param model The weight model
   */
  void SetWeightModel(TF1 *model) { fWeightModel = model; }

  double GetEventWeight(const AliMCEvent *const event) const;

private:
  TF1               *fWeightModel;    ///< Weight model
  bool               fUsePtHard;      ///< Calculate weight using pt-hard

  /// \cond CLASSIMP
  ClassDef(AliEMCalTriggerWeightHandler, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERWEIGHTHANDLER_H */
