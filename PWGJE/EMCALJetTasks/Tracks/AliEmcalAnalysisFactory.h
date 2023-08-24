/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIEMCALANALYSISFACTORY_H
#define ALIEMCALANALYSISFACTORY_H

#include <TString.h>
#include "AliEmcalTriggerOfflineSelection.h"

class AliESDtrackCuts;
class AliEmcalTrackSelection;


namespace PWGJE {
  
namespace EMCALJetTasks {

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
class AliEmcalAnalysisFactory {
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
   * @brief Get name of the default track container
   *
   * In case of usedefault the default track container is
   * used in the analysis. The name differs between ESDs and
   * AODs. The function returns the proper name for the
   * input data type.
   *
   * @param[in] isAOD Switch between ESDs and AODs
   * @return Name of the default track container
   */
  static TString TrackContainerNameFactory(Bool_t isAOD);

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

  /**
   * Generate default cut settings for the analysis
   * @return Default cut settings for the analysis
   */
  static AliESDtrackCuts *GenerateDefaultCutsESD();

  /**
   * Generate cut settings with loose DCA (used for the hybrid cuts)
   * @return Cut settings with loose DCA cuts
   */
  static AliESDtrackCuts *GenerateLooseDCACutsESD();

  /**
   * Definition: ITSchi2XXXX
   * - 3 Digits before . (to be filled with 0)
   * - 1 Digit after .
   */
  static double ValueDecoder(const char *string, const char *tag);

};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALANALYSISFACTORY_H */
