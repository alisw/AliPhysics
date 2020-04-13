/************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef __ALIANALYSISTASKEMCALTRIGGERNORMALIZATION_H__
#define __ALIANALYSISTASKEMCALTRIGGERNORMALIZATION_H__

#include "AliAnalysisTaskEmcal.h"
#include <exception>
#include <string>
#include <vector>

class THistManager;

namespace PWG {

namespace EMCAL {

/**
 * @class AliAnalysisTaskEmcalTriggerNormalization
 * @brief Task for trigger normalization 
 * @ingroup EMCALFWTASKS
 * @author Markus Fasel <markus.fasel@cern.ch> Oak Ridge National Laboratory
 * @date June 12, 2019
 * 
 * # Trigger normalization based on min. bias L0 counters
 * 
 * EMCAL triggers require a coincidence with the min. bias L0 trigger. Knowing the
 * amount of min. bias L0 triggers in a trigger cluster after event selection the 
 * luminosity of the cluster and the trigger rejection factors can be calculated.
 * 
 * Trigger clusters are configured in the way that the cluster is read out if any
 * of the trigger classes in the cluster fires. Firing the trigger class depends on
 * - Trigger selection defined by the trigger class
 * - Condition (busy status, downscaling)
 * As the busy status is a feature of the cluster, it applies to all triggers. Thus
 * the luminosity of the cluster can be estimated by correcting back the number of
 * min. bias triggers recorded by the downscale factor. The trigger rejection factor
 * can then be calculated as 
 * 
 * \f$ R = N(MB) / N(EMCAL)\f$
 * 
 * The task fills a histogram hTriggerNorm with the downscale-corrected event counts
 * for the min. bias reference trigger and the various EMCAL L1 triggers. The histogram
 * is filled as function of the centrality percentile. Users have to rebin into centrality
 * classes according to their definition of the centrality classes.
 * 
 * # Basic task configuration
 * 
 * A static Add function performs the basic task configuration. However for the 
 * luminosity determination the trigger cluster (and if different the EMCAL trigger
 * cluster), as well as the L0 trigger and the min. bias ref trigger have to be
 * defined by the user. In case one of them not set the Run method will throw 
 * an exception signaling the missing setting.
 * 
 * The following example configures the task for pp, \f$\sqrt{s} = 13 TeV\f$:
 * 
 * ~~~.{cxx}
 * auto normtask = PWG::EMCAL::AliAnalysisTaskEmcalTriggerNormalization::AddTaskEmcalTriggerNormalization("normatask");
 * normtask->SetTriggerCluster("CENT");       // Matching for INT7 and EMCAL triggers
 * normtask->SetL0TriggerEMCAL("EMC7");
 * normtask->AddMBTriggerClass("INT7");
 * ~~~
 * 
 * Attention: The task must run only on runs where the EMCAL was included and was 
 * a trigger detector.
 */
class AliAnalysisTaskEmcalTriggerNormalization : public AliAnalysisTaskEmcal {
public:

  /**
   * @brief Constructor (for ROOT I/O), not to be used by the user
   */
  AliAnalysisTaskEmcalTriggerNormalization();

  /**
   * @brief Named constructor, fully initializing the object
   * 
   * Attention: Constructor does not set default values for the trigger
   * cluster, the min. bias reference trigger, and the L0 trigger used
   * for EMCAL L0. Users must defined them. If not set, the Run method
   * will throw exceptions for each of them not set.
   */
  AliAnalysisTaskEmcalTriggerNormalization(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskEmcalTriggerNormalization() {}


  /**
   * @brief Set the trigger cluster for which the luminosity is estimated
   * 
   * Function must be called by the user. If not called, Run will throw
   * a TriggerClusterNotSetException.
   * 
   * @param triggercluster Name of the trigger cluster
   */
  void SetTriggerCluster(const char *triggercluster) { fTriggerCluster = triggercluster; } 

  /**
   * @brief Set the trigger cluster used for EMCAL trigger classes
   * 
   * This function is optional and needed only in rare cases (i.e. PbPb 2015 where the EMCAL was 
   * in the CENTNOPMD cluster, the luminosity however is estimated for the CENT cluster, and both
   * clusters have the same livetime). If not specified, the cluster specified in SetTriggerCluster
   * is used for both min. bias and EMCAL triggers.
   * 
   * @param triggercluster Name of the trigger cluster
   */
  void SetTriggerClusterEMCAL(const char *triggercluster) { fTriggerClusterEMCAL = triggercluster; }

  /**
   * @brief Set the L0 trigger used for the EMCAL L1
   * 
   * Function must be called by the user. If not called, Run will throw a 
   * L0TriggerNotSetException
   * 
   * @param trigger Name of the L0 trigger used for EMCAL L1
   */
  void SetL0TriggerEMCAL(const char *trigger) { fEMCALL0trigger = trigger; }

  /**
   * @brief Specify min. bias reference trigger
   * 
   * Function must be called by the user. If not called, Run will throw a 
   * MBTriggerNotSetException. If different min. bias triggers were used
   * during a period, multiple triggers can be defined. Any of the triggers
   * will be used to calculate the downscale corrected number of min. bias
   * events.
   * 
   * @param triggerclass Min. bias reference trigger class
   */
  void AddMBTriggerClass(const char *triggerclass) { fMBTriggerClasses.emplace_back(triggerclass); }

  /**
   * @brief Request normalization as function of mulitplicity for p-Pb
   * 
   * Attention: Relies on AliMultiplicitySelection. If not found an exception 
   * is raised.
   * 
   * @param doRequest If true the multiplicity percentile is requested for pPb
   */
  void SetRequestCentralityForpPb(Bool_t doRequest) { fUseCentralityForpPb = doRequest; }

  /**
   * @brief Configure normalization task and add it to the analysis manager
   * 
   * Configuring input and output slot, and adding the task to the analysis manager.
   * Attention: The function does not set default values for the trigger
   * cluster, the min. bias reference trigger, and the L0 trigger used
   * for EMCAL L0. Users must defined them. If not set, the Run method
   * will throw exceptions for each of them not set. 
   * 
   * @param name Name of the task
   */
  static AliAnalysisTaskEmcalTriggerNormalization *AddTaskEmcalTriggerNormalization(const char *name);

  /**
   * @class TriggerClusterNotSetException
   * @brief Exception handling trigger cluster not set error
   */
  class TriggerClusterNotSetException : public std::exception {
  public:

    /**
     * @brief Constructor
     */
    TriggerClusterNotSetException() {}

    /**
     * @brief Destructor
     */
    virtual ~TriggerClusterNotSetException() throw() {}

    /**
     * @brief Create error message
     * @return error message
     */
    virtual const char *what() const throw() { return "Trigger cluster not set."; }
  };

  /**
   * @class L0TriggerNotSetException
   * @brief Excepiton handling L0 trigger for EMCAL not set error
   */
  class L0TriggerNotSetException : public std::exception {
  public:

    /**
     * @brief Constructor
     */
    L0TriggerNotSetException() {}

    /**
     * @brief Destructor
     */
    virtual ~L0TriggerNotSetException() {}

    /**
     * @brief Create error message
     * @return error message
     */
    virtual const char *what() const throw() { return "L0 trigger not defined."; }
  };

  /**
   * @class MBTriggerNotSetException
   * @brief Exception handling MB reference trigger not set error
   */
  class MBTriggerNotSetException : public std::exception {
  public:

    /**
     * @brief Constructor
     */
    MBTriggerNotSetException() {}

    /**
     * @brief Destructor
     */
    virtual ~MBTriggerNotSetException() throw() {}

    /**
     * @brief Create error message
     * @return error message
     */
    virtual const char *what() const throw() { return "No min. bias trigger defined."; }
  };

  /**
   * @class CentralityNotSetException
   * @brief Expection handling errors reading the multiplicity object from the event
   */
  class CentralityNotSetException : public std::exception {
  public:

    /**
     * @brief Constructor
     */
    CentralityNotSetException() {}
    
    /**
     * @brief Destructor
     */
    virtual ~CentralityNotSetException() throw() {}
    
    /**
     * @brief Create error message
     * @return error message
     */
    virtual const char *what() const throw() { return "Centrality estimator not available"; }
  };

protected:

  /**
   * @brief Creating normalization histogram
   * 
   * Creating also general histograms (from AliAnalysisTaskEmcal::UserCreateOutputObjects)
   */
  virtual void UserCreateOutputObjects();

  /**
   * @brief Perform action for selected events
   * 
   * Calculating the trigger weights for the min. bias reference trigger
   * and the EMCAL triggers and filling the normalization histogram.
   * 
   * @throw TriggerClusterNotSetException if the trigger cluster is not specified
   * @throw L0TriggerNotSetException if the L0 trigger for EMCAL is not defined
   * @throw MBTriggerNotSetException if the min. bias reference trigger is not specified
   * @throw CentralityNotSetException if PbPb and the tasks fails reading the multiplicity object
   */
  virtual bool Run();

  /**
   * @brief Performing run-dependent initialization
   * 
   * Loading downscale factors from the OCDB for the new run
   * 
   * @param newrun New run
   */
  virtual void RunChanged(int newrun);

  std::string MatchTrigger(EMCAL_STRINGVIEW triggerstring, const std::vector<std::string> &triggers, EMCAL_STRINGVIEW triggercluster) const;

private:
  THistManager             *fHistos;                ///< List of histograms
  std::string               fTriggerCluster;        ///< Trigger cluster       
  std::string               fTriggerClusterEMCAL;   ///< Cluster of the EMCAL triggers (if different)
  std::string               fEMCALL0trigger;        ///< L0 trigger required for L1
  std::vector<std::string>  fMBTriggerClasses;      ///< List of valid min. bias trigger classes
  Bool_t                    fUseCentralityForpPb;   ///< Request centrality also for p-Pb (from AliMultiplicitySelection)

  AliAnalysisTaskEmcalTriggerNormalization(const AliAnalysisTaskEmcalTriggerNormalization &);
  AliAnalysisTaskEmcalTriggerNormalization &operator=(const AliAnalysisTaskEmcalTriggerNormalization &);

  ClassDef(AliAnalysisTaskEmcalTriggerNormalization, 1);
};

}

}


#endif