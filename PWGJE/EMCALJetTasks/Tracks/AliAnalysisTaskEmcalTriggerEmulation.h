#ifndef ALIANALYSISTASKEMCALTRIGGEREMULATION_H
#define ALIANALYSISTASKEMCALTRIGGEREMULATION_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcal.h"
#include "AliCutValueRange.h"
#include <TString.h>

class THistManager;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliEmcalTriggerEmulation;
class AliEMCalTriggerWeightHandler;

/**
 * @class AliAnalysisTaskEmcalTriggerEmulation
 * @brief Emulate EMCAL trigger offline in simulation, and monitor reference distribution
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since May 4, 2016
 */
class AliAnalysisTaskEmcalTriggerEmulation : public AliAnalysisTaskEmcal {
public:
  /**
   * Dummy constructor (for ROOT I/O)
   */
  AliAnalysisTaskEmcalTriggerEmulation();
  /**
   * Constructor, initializing the task with name
   * @param[in] name Name of the task
   */
  AliAnalysisTaskEmcalTriggerEmulation(const char *name);
  /**
   * Destructor
   */
  virtual ~AliAnalysisTaskEmcalTriggerEmulation();

  /**
   * Specify handler with bin weights for a certain \f$ p_{t} \f$-hard production
   * @param[in] handler Handler with bin weights
   */
  void SetWeightHandler(const AliEMCalTriggerWeightHandler * handler) { fkWeightHandler = handler; }
  /**
   * Set offline trigger emulation (used only in case of running on MC simulation)
   * @param emul EMCAL offline trigger emulation
   */
  void SetTriggerEmulation(const AliEmcalTriggerEmulation *emul) { fkTriggerEmulation = emul; }

  /**
   * Set the name of the MC particle container
   * @param[in] name Name of the MC particle container
   */
  void SetNameMCParticles(const TString &name)  { fNameMCParticles = name; }
  /**
   * Set the name of the cluster container
   * @param[in] name Name of the cluster container
   */
  void SetNameClusters(const TString &name)     { fNameClusters = name; }
  /**
   * Set the name of the track container
   * @param name Name of the track container
   */
  void SetNameTracks(const TString &name)       { fNameTracks = name; }

  /**
   * Specify whether task runs on MC events
   * @param[in] isMC Flag specifying whether input events come from simulations
   */
  void SetIsMC(bool isMC) {fIsMC = isMC; }
  /**
   * Set the kinemactic acceptance for MC particles and tracks in \f$ \eta \f$-direction
   * @param[in] min Minimum of the kinematic acceptance
   * @param[in] max Maximum of the kinematic acceptance
   */
  void SetEtaRange(double min, double max) { fEtaRange.SetLimits(min, max); }
  /**
   * Set the kinemactic acceptance for MC particles and tracks in \f$ \phi \f$-direction
   * @param[in] min Minimum of the kinematic acceptance
   * @param[in] max Maximum of the kinematic acceptance
   */
  void SetPhiRange(double min, double max) { fPhiRange.SetLimits(min, max); }

protected:
  virtual void UserCreateOutputObjects();
  virtual bool Run();

private:
  const AliEMCalTriggerWeightHandler          *fkWeightHandler;         ///< Handler containing bin weights for \f$ p_{t} \f$-hard productions
  const AliEmcalTriggerEmulation              *fkTriggerEmulation;      ///< Offline emulation of the EMCAL trigger
  THistManager                                *fHistManager;            ///< Handler for output histograms

  TString                                     fNameMCParticles;         ///< Name of the MC particle container
  TString                                     fNameClusters;            ///< Name of the cluster container
  TString                                     fNameTracks;              ///< Name of the track container

  Bool_t                                      fIsMC;                    ///< Specify whether running on MC events
  AliCutValueRange<double>                    fEtaRange;                ///< Acceptance range for particles and tracks in \f$ \eta \f$-direction
  AliCutValueRange<double>                    fPhiRange;                ///< Acceptance range for particles and tracks in \f$ \phi \f$-direction

  ClassDef(AliAnalysisTaskEmcalTriggerEmulation, 1)
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKEMCALTRIGGEREMULATION_H */
