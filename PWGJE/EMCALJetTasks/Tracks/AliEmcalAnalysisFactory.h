#ifndef ALIEMCALANALYSISFACTORY_H
#define ALIEMCALANALYSISFACTORY_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TString.h>
#include "AliEmcalTriggerOfflineSelection.h"

class AliEmcalTrackSelection;

namespace EMCalTriggerPtAnalysis {

/**
 * @class AliEmcalAnalysisFactory
 * @brief Collection of helper functions used to configure the analysis
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @ingroup PWGJETASKS
 * @since Feb 23, 2016
 *
 * Helpers implemented in this class are:
 * - Configurator for track selection
 * - Configurator for trigger selection
 * - Name of the default track container (differs between ESDs and AODs)
 */
class AliEmcalAnalysisFactory : public TObject {
public:
  AliEmcalAnalysisFactory() {}
  virtual ~AliEmcalAnalysisFactory(){}

  /**
   * @brief Get name of the default cluster container
   *
   * In case of usedefault the default cluster container is
   * used in the analysis. The name differs between ESDs and
   * AODs. The function returns the proper name for the
   * input data type.
   *
   * @param[in] isAOD Switch between ESDs and AODs
   * @return Name of the default cluster container
   */
  static TString ClusterContainerNameFactory(Bool_t isAOD);

  /**
   * @brief Fully-configure EMCAL track selection independent of the data type
   *
   * Default configurations can be used to configure the virtual track selection
   * independent of the input data type (ESDs or AODs) - data type dependent handling
   * is hidden from the user. Selections for AOD are of course only handled if AODs
   * support this. Available configurations:
   * - standard (RAA) track cuts
   * - hybrid track cuts
   *
   * @param name Name of the track cuts, used as cutstring to configure the track selection
   * @param isAOD Switch whether to create cuts for ESDs or AODs
   * @return Fully configured EMCAL track selection
   */
  static AliEmcalTrackSelection *TrackCutsFactory(TString name, Bool_t isAOD);

  /**
   * @brief Configures EMCAL trigger offline selection used to restrict EMCAL triggered sample
   *
   * Defines the threshold at Level0 and Level1 for the various
   * EMCAL triggers, and specifies the type of energy measurement
   * used to select EMCAL trigger patches
   *
   * @param[in] el0 Energy threshold for EMCAL Level0 trigger
   * @param[in] eg1 Energy threshold for EMCAL Level1 G1 trigger
   * @param[in] eg2 Energy threshold for EMCAL Level1 G2 trigger
   * @param[in] ej1 Energy threshold for EMCAL Level1 J1 trigger
   * @param[in] ej2 Energy threshold for EMCAL Level1 J2 trigger
   * @param[in] endef Energy type (ADC, offline energy, transverse energy / ADC) used for the patch selection
   * @return Fully configured EMCAL trigger offline selection
   */
  static AliEmcalTriggerOfflineSelection *TriggerSelectionFactory(Double_t el0, Double_t eg1, Double_t eg2, Double_t ej1, Double_t ej2, AliEmcalTriggerOfflineSelection::EmcalEnergyDefinition_t endef = AliEmcalTriggerOfflineSelection::kFEEEnergy);

  /// \cond CLASSIMP
  ClassDef(AliEmcalAnalysisFactory, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALANALYSISFACTORY_H */
