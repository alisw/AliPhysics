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
   * @enum TriggerCluster_t
   * @brief Numerical index for the trigger cluster selection
   */
  enum TriggerCluster_t {
    kTrgClusterANY,             ///< Any trigger cluster
    kTrgClusterCENT,            ///< CENT cluster (non-exclusive)
    kTrgClusterCENTNOTRD,       ///< CENTNOTRD cluster (non-exclusive)
    kTrgClusterCALO,            ///< CALO cluster (non-exclusive)
    kTrgClusterCALOFAST,        ///< CALOFAST cluster (non-exclusive)
    kTrgClusterCENTBOTH,        ///< CENT and CENTNOTRD
    kTrgClusterOnlyCENT,        ///< CENT cluster (exclusive)
    kTrgClusterOnlyCENTNOTRD,   ///< CENTNOTRD cluster (exclusive)
    kTrgClusterCALOBOTH,        ///< CALO and CALOFAST cluster
    kTrgClusterOnlyCALO,        ///< CALO cluster (exclusive)
    kTrgClusterOnlyCALOFAST,    ///< CALOFAST cluster (exclusive)
    kTrgClusterCENTNOPMD,       ///< CENTNOPMD cluster (non-exclusive)
    kTrgClusterALL,             ///< ALL cluster (CENT+MUON, non-exclusive)
    kTrgClusterALLNOTRD,        ///< ALLNOTRD cluster
    kTrgClusterALLBOTH,         ///< ALL and ALLNOTRD cluster
    kTrgClusterOnlyALL,         ///< ALL cluster (exclusive)
    kTrgClusterOnlyALLNOTRD,    ///< ALLNOTRD cluster (non-exclusive)
    kTrgClusterN                ///< Number of trigger clusters
  };

  /**
   * @brief Get the trigger cluster labels for a given trigger cluster index
   * @param index Index representation of the trigger cluster
   * @return std::string label of the trigger cluster
   */
  static std::string GetTriggerClusterLabels(TriggerCluster_t &index) { return fgkTriggerClusterLabels[int(index)]; }

  /**
   * @brief Get the trigger cluster labels for a given trigger cluster index
   * @param index Index representation of the trigger cluster
   * @return std::string label of the trigger cluster
   */
  static std::string GetTriggerClusterLabels(int &index) { return fgkTriggerClusterLabels[index]; }

  /**
   * @brief get the trigger cluster index from a given trigger cluster label
   * @param triggerclustername Label of the trigger cluster
   * @return Index of the trigger cluster (-1 if label is unknown) 
   */
  static int GetIndexFromTriggerClusterLabel(EMCAL_STRINGVIEW triggerclustername);

  /**
   * @brief Constructor (for ROOT I/O), not to be used by the user
   */
  AliAnalysisTaskEmcalTriggerNormalization();

  /**
   * @brief Named constructor, fully initializing the object
   * @param name Name of the task
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
   * @brief Specify min. bias reference trigger
   * @param triggerclass Min. bias reference trigger class
   * 
   * Function must be called by the user. If not called, Run will throw a 
   * MBTriggerNotSetException. If different min. bias triggers were used
   * during a period, multiple triggers can be defined. Any of the triggers
   * will be used to calculate the downscale corrected number of min. bias
   * events.
   */
  void AddMBTriggerClass(const char *triggerclass) { fMBTriggerClasses.emplace_back(triggerclass); }

  /**
   * @brief Request normalization as function of mulitplicity for p-Pb
   * @param doRequest If true the multiplicity percentile is requested for pPb
   * 
   * Attention: Relies on AliMultiplicitySelection. If not found an exception 
   * is raised.
   */
  void SetRequestCentralityForpPb(Bool_t doRequest) { fUseCentralityForpPb = doRequest; }

  /**
   * @brief Configure normalization task and add it to the analysis manager
   * @param name Name of the task
   * 
   * Configuring input and output slot, and adding the task to the analysis manager.
   * Attention: The function does not set default values for the trigger
   * cluster, the min. bias reference trigger, and the L0 trigger used
   * for EMCAL L0. Users must defined them. If not set, the Run method
   * will throw exceptions for each of them not set. 
   */
  static AliAnalysisTaskEmcalTriggerNormalization *AddTaskEmcalTriggerNormalization(const char *name);

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
   * @throw TriggerClusterNotSetException if the trigger cluster is not specified
   * @throw MBTriggerNotSetException if the min. bias reference trigger is not specified
   * @throw CentralityNotSetException if PbPb and the tasks fails reading the multiplicity object
   * 
   * Calculating the trigger weights for the min. bias reference trigger
   * and the EMCAL triggers and filling the normalization histogram.
   */
  virtual bool Run();

  /**
   * @brief Performing run-dependent initialization
   * @param newrun New run
   * 
   * Loading downscale factors from the OCDB for the new run
   */
  virtual void RunChanged(int newrun);

#if (defined(__CINT_) && !defined(__CLING__)) || (defined(__MAKECINT__) && !defined(__ROOTCLING__))
  // ROOT5 function headers
  std::vector<PWG::EMCAL::AliAnalysisTaskEmcalTriggerNormalization::TriggerCluster_t> GetTriggerClusterIndices(EMCAL_STRINGVIEW triggerstring) const;
  std::vector<PWG::EMCAL::AliAnalysisTaskEmcalTriggerNormalization::TriggerCluster_t> GetTriggerClustersANY() const { return {kTrgClusterANY}; }
#else
  // ROOT6 function headers
  std::vector<TriggerCluster_t> GetTriggerClusterIndices(EMCAL_STRINGVIEW triggerstring) const;
  std::vector<TriggerCluster_t> GetTriggerClustersANY() const { return {kTrgClusterANY}; }
#endif

  /**
   * @brief Find fired trigger(es) class within an event trigger string
   * @param triggerstring Event trigger string to be parsed
   * @param triggers Trigger classes to be matched within 
   * @return Full trigger class string matched 
   */
  std::string MatchTrigger(EMCAL_STRINGVIEW triggerstring, const std::vector<std::string> &triggers) const;

  /**
   * @brief Reset downscale factor cache for all EMCAL triggers
   */
  void ResetDownscaleFactors();

private:
  static const std::vector<std::string> fgkTriggerClusterLabels;   ///< Labels of the trigger cluster
  THistManager                         *fHistos;                    ///< List of histograms
  std::vector<std::string>              fMBTriggerClasses;          ///< List of valid min. bias trigger classes
  std::map<std::string, std::string>    fCacheTriggerClasses;       ///< List of trigger classes available for run
  std::map<std::string, double>         fCacheDownscaleFactors;     ///< Cache for downscale factor for a given run
  Bool_t                                fUseCentralityForpPb;       ///< Request centrality also for p-Pb (from AliMultiplicitySelection)

  AliAnalysisTaskEmcalTriggerNormalization(const AliAnalysisTaskEmcalTriggerNormalization &);
  AliAnalysisTaskEmcalTriggerNormalization &operator=(const AliAnalysisTaskEmcalTriggerNormalization &);

  ClassDef(AliAnalysisTaskEmcalTriggerNormalization, 1);
};

}

}


#endif