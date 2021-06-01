/**
 * \file AliEMCalTriggerMCJetAnalysisComponent.h
 * \brief Definition of class AliEMCalTriggerMCJetAnalysisComponent
 *
 * This class contains the definition of the component analysing jets and
 * particles in jets at generator level.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec. 12, 2014
 */
#ifndef ALIEMCALTRIGGERMCJETANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERMCJETANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliEMCalTriggerTracksAnalysisComponent.h"

class TString;
class AliVParticle;
class AliEmcalJet;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliEMCalTriggerEventData;

/**
 * \class AliEMCalTriggerMCJetAnalysisComponent
 * \brief Analysis component for particles in jets at generator level.
 *
 * Analysis component for tracks in jets of MC particles where the jet has a given
 * minimum pt
 */
class AliEMCalTriggerMCJetAnalysisComponent: public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerMCJetAnalysisComponent();
  AliEMCalTriggerMCJetAnalysisComponent(const char * name);
  virtual ~AliEMCalTriggerMCJetAnalysisComponent() {}

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

  /**
   * Specify minimum \f$ p_{t} \f$ for selected jets.
   * \param minpt The minimum \f$ p_{t} \f$
   */
  void SetMinimumJetPt(Double_t minpt) { fMinimumJetPt = minpt; }

protected:
  void FillHistogram(const TString &histname, const AliVParticle *track, const AliEmcalJet *jet, double vz, double weight);
  void FillJetHistogram(const TString &histname, const AliEmcalJet *recjet, double vz, double weight);

  Double_t                fMinimumJetPt;                      ///< Min. \f$ p_{t} \f$ request for the jet

  ClassDef(AliEMCalTriggerMCJetAnalysisComponent, 1);         // Analysis component for MC Jets
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERMCJETANALYSISCOMPONENT_H */
