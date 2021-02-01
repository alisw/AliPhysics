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
#ifndef ALIANALYSISTASKEMCALTRIGGERSELECTION_H
#define ALIANALYSISTASKEMCALTRIGGERSELECTION_H

#include <iosfwd>
#include <exception>
#include <string>
#include <vector>
#include <TList.h>
#include <TNamed.h>
#include <TString.h>
#include "AliAnalysisTaskEmcal.h"
#include "AliEmcalTriggerSelectionCuts.h"

// operator<< has to be forward declared carefully to stay in the global namespace so that it works with CINT.
// For generally how to keep the operator in the global namespace, See: https://stackoverflow.com/a/38801633
namespace PWG { namespace EMCAL { class AliAnalysisTaskEmcalTriggerSelection; } }
std::ostream & operator<< (std::ostream &in, const PWG::EMCAL::AliAnalysisTaskEmcalTriggerSelection &task);

namespace PWG{
namespace EMCAL {

class AliEmcalTriggerDecision;
class AliEmcalTriggerDecisionContainer;
class AliEmcalTriggerSelection;

/**
 * @class AliAnalysisTaskEmcalTriggerSelection
 * @brief Task providing an event selection for EMCAL-triggered events based on the
 * reconstructed EMCAL trigger patches
 * @ingroup EMCALFWTASKS
 * @author Markus Fasel <markus.fasel@cern.ch> Oak Ridge National Laboratory
 * @since Dec 17, 2014
 * 
 * # Perparation wagon for transparent trigger selection by users
 * 
 * Trigger event selection is centralized in this task. The trigger selecion
 * is implemented in trigger selecion object, selecting all trigger patches
 * firing the trigger using a trigger patch cut object fully describing the trigger.
 * Trigger selections have to be configured by the user and should fully describe 
 * all Level1 triggers supported in the corresponding dataset (for simulation the 
 * dataset it is anchored to). Typically different modes are used on data and
 * simulation:
 * 
 * - Data: Selection based on FastOR ADC, used to cleanup noisy triggers or 
 *   equalize triggers for different periods with different settings
 * - Simulation: Selection based on FEE energy with an energy threshold corresponding
 *   to the ADC threshold applied in data, used to simulate the trigger response
 * 
 * Users have to configure the task acording to the different modes and the settings
 * matching the dataset. For several datasets default configurations are already 
 * implemented. Missing configurations are implemented by the author on user request.
 * 
 * The trigger selection tasks appends an object to the event where users can easily
 * query the trigger selection so the selection based on trigger primitives (patches)
 * is centralized and not needed in the user analysis. The object can be retrieved
 * by the users from the event.
 * 
 * The task has also a builtin QA for each trigger class, monitoring the energy spectra
 * of all trigger patches and of the main trigger patches. It is recommended to always
 * run the trigger selection QA - resource consumption should not be an issue.
 * 
 * ## Adding the wagon to the train
 * 
 * An add macro is provided under PWG/EMCAL/AddEmcalTriggerSelectionTask.C which can be
 * used to configure the train wagon. The macro customization should include the following
 * lines:
 * 
 * ~~~{.cxx}
 * __R_ADDTASK__->SetGlobalDecisionContainerName("EmcalTriggerDecision");
 * __R_ADDTASK__->AutoConfigure(kDataset);
 * ~~~
 *
 * kDataset should be the name of a supported dataset. Please refer to the description of class 
 * PWG::EMCAL::AliAnalysisTaskEmcalTriggerSelection for supported datasets.
 * 
 * ## Configuring the trigger selection task for the corresponding period
 * 
 * Default configurations exists for run 1 pp (2012, 8 TeV) and for run2 pp (2016-2018, 13 TeV)
 * and the corresponding datasets. For existing configurations and the corresponing cut
 * settings please refer to the matching Configure function.
 * 
 * An AutoConfigure function provides a frontend for the user to configure the task based
 * on the dataset name, and internally mapping datasets to configurations. On the data side
 * pp datasets for 2012 and 2016-2018 are fully supported. Unfortunately the situation is 
 * more complicated on the simulation side: Due to the absence of a logical naming scheme
 * where an algorithm can extract whether it is a simulation dataset or what period it is 
 * anchored to, supported datasets must be handled by a lookup table, and this can by far
 * not be complete. Users should contact the author in case a new dataset should be added
 * to the configuration.
 * 
 * Configuration by hand is also possible. This is of relevance for systematic studies of
 * the trigger response via cut variations. Users have to provide a trigger selection via
 * the function AddTriggerSelection. Each trigger selection object represents one Level1
 * trigger. For such configurations it is recommended to use a config macro.
 * 
 * ## Accessing the trigger selection in the user task
 * 
 * After processing the trigger selection the task appends a container with the selection
 * results for each trigger supported for the dataset to the input event. The trigger decision
 * container provides information about
 * 
 * - Trigger selection status
 * - Main patch according to the energy and patch specification
 * - All trigger patches firing the trigger
 * - The trigger selection cuts
 * 
 * Users can use the trigger selection container in order to select events as triggered events
 * for the given trigger in case the selection status is true. The name of the trigger must
 * match the name in the configuration. Names of trigger classes follow the EMCAL naming convention.
 * The following example reads the trigger decision container from the event and selects events
 * as EG1 event:
 * 
 * ~~~{.cxx}
 * auto trgsel = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(input->FindListObject("EmcalTriggerDecision"))
 * if(trgsel->IsEventSelected("EG1")) std::cout << "Event selected as EG1 event" << std::endl; 
 * ~~~
 */
class AliAnalysisTaskEmcalTriggerSelection: public AliAnalysisTaskEmcal {
public:

  /**
   * @class ConfigValueException
   * @brief Handling of incorrect values in YAML configuration files
   * 
   * Many information (Acceptance, patch type, ...) are represented in the 
   * YAML configuration file as strings. Thus they correspond to a finite
   * set of values, usually handled as enumeration type. This class handles
   * the error raised for improper configuration values;
   */
  class ConfigValueException : public std::exception {
  public:
    /**
     * @brief Construct a new ConfigValueException object
     * 
     * Exception is thrown when decoding a configuration string with an unknown value
     * 
     * @param key   Key for which an improper value was set
     * @param value Improper value
     */
    ConfigValueException(const char *key, const char *value);

    /**
     * @brief Destructor
     */
    virtual ~ConfigValueException() throw() {}

    /**
     * @brief Display error message
     * @return Error message string
     */
    const char *what() const throw() { return fMessage.data(); }

    const std::string &getKey() const { return fKey; }
    const std::string &getValue() const { return fValue; }

  private:
    std::string           fKey;       ///< Key for which an unknown value was assigned
    std::string           fValue;     ///< Improper value raising the exception
    std::string           fMessage;   ///< Error message shown in what()
  };
  /**
   * @brief Dummy constructor
   *
   * Only for I/O, not to be used by the user
   */
  AliAnalysisTaskEmcalTriggerSelection();

  /**
   * @brief Main constructor
   *
   * Initializing an empty task with default settings 
   * (i.e. name of the trigger decision container). 
   * The task is not yet configured and has no trigger
   * selection object attached - users must use the 
   * function AutoConfigure or an appropriate Configure
   * function. The name of the task is used to retrieve
   * the task from the AliAnlysisManager and must be
   * unique.
   * 
   * @param[in] name Name of the task
   */
  AliAnalysisTaskEmcalTriggerSelection(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskEmcalTriggerSelection() {}

  /**
   * @brief Add trigger selection to the trigger selection task
   * 
   * This function is used to configure the trigger maker manually
   * by providing trigger selection objects configured by hand.
   * Each trigger selection object corresponds to one trigger
   * class. The method is forseen for systematic settings. For 
   * default settings please refer to the AutoConfigure function
   * or the appropriate Configure function.
   *
   * @param[in] selection the trigger selection to be added
   */
  void AddTriggerSelection(AliEmcalTriggerSelection * const selection);

  /**
   * @brief Access to the trigger selection objects attached to the task
   * @return List of trigger selections
   */
  const TList &GetListOfTriggerSelections() const { return fTriggerSelections; }

  /**
   * @brief Set the name of the global trigger decision container
   *
   * The AliAnalysisTaskEmcalTriggerSelection appends an object
   * of type PWG::EMCAL::AliEmcalTriggerDecisionContainer with
   * the status of the trigger selecion for all trigger classes 
   * supported by the configuration to the event. Users can query
   * the object in their task via
   * 
   * ~~~{.cxx}
   * auto trgcont = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(nametriggercontainer));
   * ~~~
   * 
   * The name must match the name specified here. If no name
   * is specified the default one ("EmcalTriggerDecision") is 
   * used.
   *
   * @param[in] name Name of the trigger decision container
   */
  void SetGlobalDecisionContainerName(const char *name) { fGlobalDecisionContainerName = name; }

  /**
   * @brief Provide access to the name of the global trigger decision container name
   * @return Name of the trigger decision container
   */
  const TString  &GetGlobalDecisionContainerName() const { return fGlobalDecisionContainerName; }

  /**
   * @brief Configure task using YAML configuration file
   * 
   * This interface allows to setup a complicated trigger scheme
   * using a singe YAML file. The YAML file consists of a global
   * part with settings in common for all triggers and a trigger
   * specific part. Global configurations are (keyname):
   * 
   * - containername: Name of the output container
   * - energydef: Energy definition
   * - energysource: Offline (FEE) or Recalc (ADC, ignoring online STU decision)
   * - triggerclasses: Array with the names of the trigger classes
   * 
   * Supported values for key energydef:
   * 
   * | Value name    | Energy definition                |
   * |---------------|----------------------------------|
   * | ADC           | ADC from FastORs                 |
   * | Energy        | FEE Energy                       |
   * | EnergyRough   | Energy estimated from FastORs    |
   * | EnergySmeared | FEE Energy with offline smearing |
   * 
   * The key for the trigger class is the name of the trigger class.
   * Configurations for the specific triggers are (keyname):
   * 
   * - acceptance: Acceptance type
   * - patchtype: Patch type
   * - Threshold: ADC / Energy threshold
   * 
   * Supported values for acceptance:
   * 
   * | Value name | Acceptance type      | 
   * |------------|----------------------|
   * | EMCAL      | EMCAL acceptance     |
   * | DCAL       | DCAL+PHOS acceptance |
   * 
   * Supported values for patchtype:
   * 
   * | Value type  | Type of the patch                  |
   * |-------------|------------------------------------|
   * | L1Gamma     | Level1 gamma patch, any threshold  |
   * | L1GammaHigh | Level1 gamma patch, high threshold |
   * | L1GammaLow  | Level1 gamma patch, low threshold  |
   * | L1Jet       | Level1 jet patch, any threshold    |
   * | L1JetHigh   | Level1 jet patch, high threshold   |
   * | L1JetLow    | Level1 jet patch, low threshold    |
   * 
   * Here is an example for a valid config file
   * 
   * ~~~{.yaml}
   * containername: "EmcalTriggerDecisionV1"
   * energysource: "Recalc"
   * energydef: "ADC"
   * triggerclasses:
   * - EGA
   * - EJE
   * EGA:
   *     acceptance: "EMCAL"
   *     patchtype" "L1Gamma"
   *     threshold: 130
   * EJE:
   *     acceptance: "EMCAL"
   *     patchtype: "L1Jet"
   *     threshold: 200
   * ~~~
   * 
   * @param yamlconfig Name of the YAML configuration file
   */
  void ConfigureFromYAML(const char *yamlconfig);

  /**
   * @brief Automatically configure trigger decision handler for different periods
   */
  void AutoConfigure(const char *period);

  /**
   * @brief Trigger configuration for data anchored to run1 pp 7TeV (2011)
   *
   * Configuration is representing the Level0 trigger available in the
   * run1 data taking for pp in 2011. This configuration is for the 
   * simulation mode, correspondingly the trigger patch selection
   * is applied on trigger patches calulated from FEE energies with energy 
   * smearing applied. The following table lists the trigger classes supported 
   * together with the corresponding settings:
   * 
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (ADC) |
   * |---------------|-----------------|----------------|-----------------|
   * |       L0      |     EMCAL       |   EMCL0 (2x2)  |       292       |
   * 
   * The size of the jet patch and the subregion size are defined in the configuration
   * of the trigger maker kernel.
   * 
   * The trigger thresholds are tuned to describe the data.
   */
  void ConfigurePP7TeV2011();

  /**
   * @brief Trigger configuration for MC anchored to run1 pp 7TeV (2011)
   *
   * Configuration is representing the Level9 trigger available in the
   * run1 data taking for pp in 2011. This configuration is for the 
   * simulation mode, correspondingly the trigger patch selection
   * is applied on trigger patches calulated from FEE energies with energy 
   * smearing applied. The following table lists the trigger classes supported 
   * together with the corresponding settings:
   * 
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (GeV) |
   * |---------------|-----------------|----------------|-----------------|
   * |       L0      |     EMCAL       |   EMCL0 (2x2)  |       5.5       |
   * 
   * The size of the jet patch and the subregion size are defined in the configuration
   * of the trigger maker kernel.
   * 
   * The trigger thresholds are tuned to describe the data.
   */
  void ConfigureMCPP7TeV2011();
  /**
   * @brief Trigger configuration for run1 pp (2012) - data mode
   * 
   * Configuration is representing all Level1 triggers available in the
   * run1 data taking for pp in 2012. This configuration is for the 
   * online mode (data events), correspondingly the trigger patch selection
   * is applied on trigger patches fulfilling the recalc type. The following
   * table lists the trigger classes supported together with the corresponding
   * settings:
   * 
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (ADC) |
   * |---------------|-----------------|----------------|-----------------|
   * |     EGA       |     EMCAL       |   EGA (2x2)    |       130       |
   * |     EJE       |     EMCAL       |   EJE (16x16)  |       200       |
   * 
   * The size of the jet patch and the subregion size are defined in the configuration
   * of the trigger maker kernel.
   * 
   * The trigger thresholds match the ones in https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalTriggerOffline
   */
  void ConfigurePP2012();

  /**
   * @brief Trigger configuration for MC anchored to run1 pp (2012)
   *
   * Configuration is representing all Level1 triggers available in the
   * run2 data taking for pp in 2012. This configuration is for the 
   * simulation mode, correspondingly the trigger patch selection
   * is applied on trigger patches calulated from FEE energies with energy 
   * smearing applied. The following table lists the trigger classes supported 
   * together with the corresponding settings:
   * 
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (GeV) |
   * |---------------|-----------------|----------------|-----------------|
   * |     EG1       |     EMCAL       |   EGA (2x2)    |       10        |
   * |     EJ1       |     EMCAL       |   EJE (16x16)  |       15.5      |
   * 
   * The size of the jet patch and the subregion size are defined in the configuration
   * of the trigger maker kernel.
   * 
   * The trigger thresholds are tuned to describe the data.
   */
  void ConfigureMCPP2012();
  /**
   * @brief Trigger configuration for run1 pPb (2013) - data mode
   * 
   * Configuration is representing all Level1 triggers available in the
   * run1 data taking for pPb in 2013. This configuration is for the 
   * online mode (data events), correspondingly the trigger patch selection
   * is applied on trigger patches fulfilling the recalc type. The following
   * table lists the trigger classes supported together with the corresponding
   * settings:
   * 
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (ADC) |
   * |---------------|-----------------|----------------|-----------------|
   * |      L0       |     EMCAL       |   L0 (2x2)     |       158       |
   * |     EG1       |     EMCAL       |   EGA (2x2)    |       140       |
   * |     EG2       |     EMCAL       |   EGA (2x2)    |        89       |
   * |     EJ1       |     EMCAL       |   EJE (16x16)  |       260       |
   * |     EJ2       |     EMCAL       |   EJE (16x16)  |       127       |
   * 
   * The size of the jet patch and the subregion size are defined in the configuration
   * of the trigger maker kernel.
   * 
   * The trigger thresholds match the ones in https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalTriggerOffline
   */
  void ConfigurePPB5TeV2013();

  /**
   * @brief Trigger configuration for MC anchored to run1 pPb (2013)
   *
   * Configuration is representing all Level1 triggers available in the
   * run2 data taking for pPb in 2013. This configuration is for the 
   * simulation mode, correspondingly the trigger patch selection
   * is applied on trigger patches calulated from FEE energies with energy 
   * smearing applied. The following table lists the trigger classes supported 
   * together with the corresponding settings:
   * 
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (GeV) |
   * |---------------|-----------------|----------------|-----------------|
   * |      L0       |     EMCAL       |   L0 (2x2)     |       3.0       |
   * |     EG1       |     EMCAL       |   EGA (2x2)    |       11.       |
   * |     EG2       |     EMCAL       |   EGA (2x2)    |       7.0       |
   * |     EJ1       |     EMCAL       |   EJE (16x16)  |       20.       |
   * |     EJ2       |     EMCAL       |   EJE (16x16)  |       10.       |
   * 
   * The size of the jet patch and the subregion size are defined in the configuration
   * of the trigger maker kernel.
   * 
   * The trigger thresholds are tuned to describe the data.
   */
  void ConfigureMCPPB5TeV2013();

  /**
   * @brief Trigger configuration for run2 pp 5TeV (2015) - data mode
   * 
   * Configuration is representing all Level0 triggers available in the
   * run2 data taking for pp 2015. This configuration is for the 
   * online mode (data events), correspondingly the trigger patch selection
   * is applied on trigger patches fulfilling the recalc type. The following
   * table lists the trigger classes supported together with the corresponding
   * settings:
   * 
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (ADC) |
   * |---------------|-----------------|----------------|-----------------|
   * |      L0       |     EMCAL       |   L0 (2x2)     |       263       |
   * |      L0       |     DCAL        |   L0 (2x2)     |       263       |
   * 
   * The trigger thresholds match the ones in https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalTriggerOffline
   */
  void ConfigurePP5TeV2015();

  /**
   * @brief Trigger configuration for MC anchored to run2 pp 5TeV (2015)
   *
   * Configuration is representing all Level0 triggers available in the
   * run2 data taking for pp 5TeV 2015. This configuration is for the 
   * simulation mode, correspondingly the trigger patch selection
   * is applied on trigger patches calulated from FEE energies with energy 
   * smearing applied. The following table lists the trigger classes supported 
   * together with the corresponding settings:
   * 
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (GeV) |
   * |---------------|-----------------|----------------|-----------------|
   * |      L0       |     EMCAL       |   L0 (2x2)     |       5.0       |
   * |      L0       |     DCAL        |   L0 (2x2)     |       5.0       |
   * 
   * The trigger thresholds are tuned to describe the data.
   */
  void ConfigureMCPP5TeV2015();

  /**
   * @brief Trigger configuration for run2 pp (2016 - 2018) - data mode
   * 
   * Configuration is representing all Level1 triggers available in the
   * run2 data taking for pp 2016 to 2018. This configuration is for the 
   * online mode (data events), correspondingly the trigger patch selection
   * is applied on trigger patches fulfilling the recalc type. The following
   * table lists the trigger classes supported together with the corresponding
   * settings:
   * 
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (ADC) |
   * |---------------|-----------------|----------------|-----------------|
   * |     EG1       |     EMCAL       |   EGA (2x2)    |       115       |
   * |     EG2       |     EMCAL       |   EGA (2x2)    |        51       |
   * |     DG1       |     DCAL        |   EGA (2x2)    |       115       |
   * |     DG2       |     DCAL        |   EGA (2x2)    |        51       |
   * |     EJ1       |     EMCAL       |   EJE (16x16)  |       255       |
   * |     EJ2       |     EMCAL       |   EJE (16x16)  |       204       |
   * |     DJ1       |     DCAL        |   EJE (8x8)    |       255       |
   * |     DJ2       |     DCAL        |   EJE (8x8)    |       204       |
   * 
   * The size of the jet patch and the subregion size are defined in the configuration
   * of the trigger maker kernel.
   * 
   * The trigger thresholds match the ones in https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalTriggerOffline
   * 
   * Attention: Selection of jet triggered events will not work for LHC16o+p, incorrect
   * patchsize was applied at hardware level, therefore too strong thresholds.
   */
  void ConfigurePP2016();

  /**
   * @brief Trigger configuration for MC anchored to run2 pp (2016-2018)
   *
   * Configuration is representing all Level1 triggers available in the
   * run2 data taking for pp 2016 to 2018. This configuration is for the 
   * simulation mode, correspondingly the trigger patch selection
   * is applied on trigger patches calulated from FEE energies with energy 
   * smearing applied. The following table lists the trigger classes supported 
   * together with the corresponding settings:
   * 
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (GeV) |
   * |---------------|-----------------|----------------|-----------------|
   * |     EG1       |     EMCAL       |   EGA (2x2)    |         9       |
   * |     EG2       |     EMCAL       |   EGA (2x2)    |         4       |
   * |     DG1       |     DCAL        |   EGA (2x2)    |         9       |
   * |     DG2       |     DCAL        |   EGA (2x2)    |         4       |
   * |     EJ1       |     EMCAL       |   EJE (16x16)  |        14       |
   * |     EJ2       |     EMCAL       |   EJE (16x16)  |        14       |
   * |     DJ1       |     DCAL        |   EJE (8x8)    |        19       |
   * |     DJ2       |     DCAL        |   EJE (8x8)    |        14       |
   * 
   * The size of the jet patch and the subregion size are defined in the configuration
   * of the trigger maker kernel.
   * 
   * The trigger thresholds are tuned to describe the data.
   */
  void ConfigureMCPP2016();

  /**
   * @brief Trigger configuration for run2 pPb 5TeV (2016) - data mode
   * 
   * Configuration is representing all Level1 triggers available in the
   * run2 data taking for pPb 5TeV 2016. This configuration is for the 
   * online mode (data events), correspondingly the trigger patch selection
   * is applied on trigger patches fulfilling the recalc type. The following
   * table lists the trigger classes supported together with the corresponding
   * settings:
   * 
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (ADC) |
   * |---------------|-----------------|----------------|-----------------|
   * |     EG1       |     EMCAL       |   EGA (2x2)    |       140       |
   * |     EG2       |     EMCAL       |   EGA (2x2)    |        83       |
   * |     DG1       |     DCAL        |   EGA (2x2)    |       140       |
   * |     DG2       |     DCAL        |   EGA (2x2)    |        83       |
   * |     EJ1       |     EMCAL       |   EJE (16x16)  |       318       |
   * |     EJ2       |     EMCAL       |   EJE (16x16)  |       255       |
   * |     DJ1       |     DCAL        |   EJE (8x8)    |       318       |
   * |     DJ2       |     DCAL        |   EJE (8x8)    |       255       |
   * 
   * The trigger thresholds match the ones in https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalTriggerOffline
   */
  void ConfigurePPB5TeV2016();

  /**
   * @brief Trigger configuration for MC anchored to run2 pPb 5 TeV (2016)
   *
   * Configuration is representing all Level1 triggers available in the
   * run2 data taking for pPb 5TeV 2016. This configuration is for the
   * simulation mode, correspondingly the trigger patch selection
   * is applied on trigger patches calulated from FEE energies with energy
   * smearing applied. The following table lists the trigger classes supported
   * together with the corresponding settings:
   *
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (GeV) |
   * |---------------|-----------------|----------------|-----------------|
   * |     EG1       |     EMCAL       |   EGA (2x2)    |        11.      |
   * |     EG2       |     EMCAL       |   EGA (2x2)    |        6.5      |
   * |     DG1       |     DCAL        |   EGA (2x2)    |        11.      |
   * |     DG2       |     DCAL        |   EGA (2x2)    |        6.5      |
   * |     EJ1       |     EMCAL       |   EJE (16x16)  |        25.      |
   * |     EJ2       |     EMCAL       |   EJE (16x16)  |        20.      |
   * |     DJ1       |     DCAL        |   EJE (8x8)    |        25.      |
   * |     DJ2       |     DCAL        |   EJE (8x8)    |        20.      |
   *
   * The trigger thresholds are tuned to describe the data.
   */
  void ConfigureMCPPB5TeV2016();

  /**
   * @brief Trigger configuration for run2 pPb 8TeV (2016) - data mode
   *
   * Configuration is representing all Level1 triggers available in the
   * run2 data taking for pPb 8TeV 2016. This configuration is for the
   * online mode (data events), correspondingly the trigger patch selection
   * is applied on trigger patches fulfilling the recalc type. The following
   * table lists the trigger classes supported together with the corresponding
   * settings:
   *
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (ADC) |
   * |---------------|-----------------|----------------|-----------------|
   * |     EG1       |     EMCAL       |   EGA (2x2)    |       102       |
   * |     EG2       |     EMCAL       |   EGA (2x2)    |        70       |
   * |     DG1       |     DCAL        |   EGA (2x2)    |       102       |
   * |     DG2       |     DCAL        |   EGA (2x2)    |        70       |
   * |     EJ1       |     EMCAL       |   EJE (16x16)  |       293       |
   * |     EJ2       |     EMCAL       |   EJE (16x16)  |       229       |
   * |     DJ1       |     DCAL        |   EJE (8x8)    |       293       |
   * |     DJ2       |     DCAL        |   EJE (8x8)    |       229       |
   *
   * The trigger thresholds match the ones in https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalTriggerOffline
   */
  void ConfigurePPB8TeV2016();

  /**
   * @brief Trigger configuration for MC anchored to run2 pPb 8 TeV (2016)
   *
   * Configuration is representing all Level1 triggers available in the
   * run2 data taking for pPb 8TeV 2016. This configuration is for the
   * simulation mode, correspondingly the trigger patch selection
   * is applied on trigger patches calulated from FEE energies with energy
   * smearing applied. The following table lists the trigger classes supported
   * together with the corresponding settings:
   *
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (GeV) |
   * |---------------|-----------------|----------------|-----------------|
   * |     EG1       |     EMCAL       |   EGA (2x2)    |        8.0      |
   * |     EG2       |     EMCAL       |   EGA (2x2)    |        5.5      |
   * |     DG1       |     DCAL        |   EGA (2x2)    |        8.0      |
   * |     DG2       |     DCAL        |   EGA (2x2)    |        5.5      |
   * |     EJ1       |     EMCAL       |   EJE (16x16)  |        23.      |
   * |     EJ2       |     EMCAL       |   EJE (16x16)  |        18.      |
   * |     DJ1       |     DCAL        |   EJE (8x8)    |        23.      |
   * |     DJ2       |     DCAL        |   EJE (8x8)    |        18.      |
   *
   * The trigger thresholds are tuned to describe the data.
   */
  void ConfigureMCPPB8TeV2016();

  /**
   * @brief Trigger configuration for run2 pp 5TeV (2017) - data mode
   *
   * Configuration is representing all Level1 triggers available in the
   * run2 data taking for pp 5TeV 2017. This configuration is for the
   * online mode (data events), correspondingly the trigger patch selection
   * is applied on trigger patches fulfilling the recalc type. The following
   * table lists the trigger classes supported together with the corresponding
   * settings:
   *
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (ADC) |
   * |---------------|-----------------|----------------|-----------------|
   * |     EG1       |     EMCAL       |   EGA (2x2)    |       115       |
   * |     EG2       |     EMCAL       |   EGA (2x2)    |        51       |
   * |     DG1       |     DCAL        |   EGA (2x2)    |       115       |
   * |     DG2       |     DCAL        |   EGA (2x2)    |        51       |
   * |     EJ1       |     EMCAL       |   EJE (16x16)  |       255       |
   * |     EJ2       |     EMCAL       |   EJE (16x16)  |       204       |
   * |     DJ1       |     DCAL        |   EJE (8x8)    |       255       |
   * |     DJ2       |     DCAL        |   EJE (8x8)    |       204       |
   *
   * The trigger thresholds match the ones in https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalTriggerOffline
   */
  void ConfigurePP5TeV2017();

  /**
   * @brief Trigger configuration for MC anchored to run2 pp 5 TeV (2017)
   *
   * Configuration is representing all Level1 triggers available in the
   * run2 data taking for pp 5TeV 2017. This configuration is for the
   * simulation mode, correspondingly the trigger patch selection
   * is applied on trigger patches calulated from FEE energies with energy
   * smearing applied. The following table lists the trigger classes supported
   * together with the corresponding settings:
   *
   * | Trigger class | Acceptance type |   Patch Type   | Threshold (GeV) |
   * |---------------|-----------------|----------------|-----------------|
   * |     EG1       |     EMCAL       |   EGA (2x2)    |        9.0      |
   * |     EG2       |     EMCAL       |   EGA (2x2)    |        4.0      |
   * |     DG1       |     DCAL        |   EGA (2x2)    |        9.0      |
   * |     DG2       |     DCAL        |   EGA (2x2)    |        4.0      |
   * |     EJ1       |     EMCAL       |   EJE (16x16)  |        20.      |
   * |     EJ2       |     EMCAL       |   EJE (16x16)  |        16.      |
   * |     DJ1       |     DCAL        |   EJE (8x8)    |        20.      |
   * |     DJ2       |     DCAL        |   EJE (8x8)    |        16.      |
   *
   * The trigger thresholds are tuned to describe the data.
   */
  void ConfigureMCPP5TeV2017();

  /**
   * @brief Output stream operator
   * 
   * Logging all settings (output container, trigger classes, trigger selection cuts) to the
   * output stream
   * 
   * @param stream Output stream used for logging
   * @param task Trigger selection task to be put on the stream
   * @return stream after logging
   */
  friend std::ostream& ::operator<<(std::ostream &stream, const AliAnalysisTaskEmcalTriggerSelection &task);
  
protected:

  /**
   * @class AliEmcalTriggerSelectionQA
   * @brief Helper class for the trigger selection
   * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
   * @since Oct 3, 2017
   */
  class AliEmcalTriggerSelectionQA : public TNamed {
  public:

    /**
     * @brief Dummy constructor
     */
    AliEmcalTriggerSelectionQA();

    /**
     * @brief Constructor
     *
     * Initializing component matching to the trigger selection
     *
     * @param[in] sel Corresponding trigger selection to be monitored by this QA component
     */
    AliEmcalTriggerSelectionQA(const AliEmcalTriggerSelection * const sel);

    /**
     * @brief Copy constructor
     *
     * Histograms will be shared between the two QA components
     * @param[in] ref Reference for copy
     */
    AliEmcalTriggerSelectionQA(const AliEmcalTriggerSelectionQA &ref);

    /**
     * @brief Assignment operator
     *
     * Histograms will be shared between the two QA components
     * @param[in] ref Reference for assignment
     * @return This object after assignment
     */
    AliEmcalTriggerSelectionQA &operator=(const AliEmcalTriggerSelectionQA &ref);

    /**
     * @brief Destructor
     */
    virtual ~AliEmcalTriggerSelectionQA() {}

    /**
     * Fill histograms for the main patch for the given trigger decision
     * @param[in] decision
     */
    void Fill(const AliEmcalTriggerDecision * const decision);

    /**
     * @brief Fill histograms of this QA component into the targetlist
     *
     * The outputlist will take ownership over the histograms
     *
     * @param[out] targetlist list toe
     */
    void GetHistos(TList *targetlist) const;

  private:
     TH1 *fMaxPatchADC;                 ///< Histogram with patch ADC of the max patch
     TH1 *fMaxPatchEnergy;              ///< Histogram with patch energy of the max patch
     TH1 *fMaxPatchEnergySmeared;       ///< Histogram with smeared patch energy of the max patch
  };

  /**
   * @brief Initialization of output container
   * 
   * Initialization will take care of the QA histograms specified by the 
   * various trigger selections in the corresponding QA objects.
   * 
   * Function overrides UserCreateOutputObjects from AliAnalysisTaskSE
   */
  virtual void UserCreateOutputObjects();

  /**
   * @brief Initializations performed when the first event is created
   * 
   * Initializing common output container for trigger decision
   * 
   * Function overrides UserExecOnce from AliAnalysisTaskEmcal
   */
  virtual void UserExecOnce();

  /**
   * @brief User event loop.
   * 
   * The user event loop performs the trigger selection. The selection
   * process is implemented in AliEmcalTriggerSelection objects. The 
   * event loop processes all trigger selection objects configured in
   * this task and collects the result in the trigger selection container
   * which is attached to the event.
   * 
   * Function overrides Run from AliAnalysisTaskEmcal
   * 
   * @return Always true
   */
  virtual Bool_t Run();

  /**
   * @brief Filling basic QA Histograms of the trigger selection task
   *
   * The QA histograms are connected to the main patch and monitor
   * - ADC Amplitude
   * - Energy
   * - Smeared energy
   * @return Always true
   */
  virtual Bool_t FillHistograms();

  /**
   * @brief Find the main trigger container in the input event.
   *
   * If not available, create it and add it to the input event
   */
  AliEmcalTriggerDecisionContainer *GetGlobalTriggerDecisionContainer() const { return fTriggerDecisionContainer; }

  /**
   * @brief Fill QA histograms for the event
   * @param cont
   */
  void MakeQA(const AliEmcalTriggerDecisionContainer *cont);

  /**
   * @brief Initialize QA histograms for trigger selection
   * @param sel
   */
  void InitQA(const AliEmcalTriggerSelection *const sel);

  Bool_t Is2011PP7TeV(const char *dataset) const;
  Bool_t Is2011MCPP7TeV(const char *dataset) const;
  Bool_t Is2012PP(const char *dataset) const;
  Bool_t Is2012MCPP(const char *dataset) const;
  Bool_t Is2013PPB(const char *dataset) const;
  Bool_t Is2013MCPPB(const char *dataset) const;
  Bool_t Is2015PP5TeV(const char *dataset) const;
  Bool_t Is2015MCPP5TeV(const char *dataset) const;
  Bool_t Is2016PPB5TeV(const char *dataset) const;
  Bool_t Is2016MCPPB5TeV(const char *dataset) const;
  Bool_t Is2016PPB8TeV(const char *dataset) const;
  Bool_t Is2016MCPPB8TeV(const char *dataset) const;
  Bool_t Is2016PP(const char *dataset) const;
  Bool_t Is2016MCPP(const char *dataset) const;
  Bool_t Is2017PP5TeV(const char *dataset) const;
  Bool_t Is2017MCPP5TeV(const char *dataset) const;
  Bool_t IsSupportedMCSample(const char *period, std::vector<TString> &supportedProductions) const;
 
  AliEmcalTriggerSelectionCuts::AcceptanceType_t  DecodeAcceptanceString(const std::string &acceptancestring);
  AliEmcalTriggerSelectionCuts::PatchType_t       DecodePatchTypeString(const std::string &patchtypestring);
  AliEmcalTriggerSelectionCuts::SelectionMethod_t DecodeEnergyDefinition(const std::string &energydefstring);

  /**
   * @brief Calculate rho value from reconstructed trigger patches
   * @param patches Input patch container
   * @param recalc If true the rho is calculated from recalc patches (as ADC), otherwise from offline patches
   * @param isEmcal  If true the rho is calculated from EMCAL patches, otherwise from DCAL patches
   * @return rho value for the given combination
   */
  double CalculateRho(const TClonesArray *const patches, Bool_t recalc, Bool_t isEmcal);

  AliEmcalTriggerDecisionContainer          *fTriggerDecisionContainer;        ///< Trigger decision container objects
  TString                                    fGlobalDecisionContainerName;     ///< Name of the global trigger selection
  TList                                      fTriggerSelections;               ///< List of trigger selections
  TList                                      fSelectionQA;                     ///< Trigger selection QA
  Bool_t                                     fUseRho;                          //!<! Check if any of teh trigger selections is requiring rho calculation

private:

  /**
   * @brief Print information about the trigger decision container to the output stream
   * @param stream Output stream used from printing
   */
  void PrintStream(std::ostream &stream) const;

  ClassDef(AliAnalysisTaskEmcalTriggerSelection, 1);    // Task running different EMCAL trigger selections
};


}
}

#endif /* ALIANALYSISTASKEMCALTRIGGERSELECTION_H */
