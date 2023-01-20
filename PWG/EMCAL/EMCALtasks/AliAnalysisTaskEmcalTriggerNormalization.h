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
 * # Trigger normalization based on trigger luminosity
 * 
 * In addtion the task creates a histogram hTriggerLuminosity which measures the luminosity
 * for the different triggers based on the CENT cluster inspected luminosity and the trigger
 * livetime from downscaling. In order to obtain the full luminosity of the trigger the
 * trigger luminosity needs to be corrected to the various trigger clusters. For each trigger
 * class a corresponding counter counts the events in the various trigger classes. Triggers
 * which are either not downscaled (i.e. high-threshold triggers) or synchronized with min. bias
 * (i.e. EMC7) can be used to calculate the trigger cluster luminosities in other trigger
 * clusters based on the luminosity of the CENT cluster. Which trigger class to use depends 
 * on the various datasets. The following examples can be used for the various datasets:
 * 
 * | Dataset       | Conversion trigger cluster | Trigger classes    |
 * |---------------|----------------------------|--------------------|
 * | pp 13 TeV     | CENT->CENTNOTRD            | EJ1, EG1, DJ1, DG1 |
 * | p-Pb 8.16 TeV | CENT->CENTNOPMD            | EMC7, DMC7         |
 * | p-Pb 8.16 TeV | CENTNOPMD->CENTNOTRD       | EJ1, EG1, DJ1, DG1 |
 * 
 * # Trigger correlation 
 * 
 * Further histograms are filled for each trigger cluster monitoring the correlation between
 * different trigger classes. The diagonal elements are self-correlations and can therefore
 * be used to as absolute event counter. The off-diagonal elements of the matrix monitor
 * which triggers fire together with the trigger class and show the correlation of a trigger
 * with other triggers, where the other trigger is always a subset.
 * 
 * Downscale factors can be checked comparing the event count of the off-diagonal element for
 * the trigger to the diagonal eleemnt for a reference trigger. In this case the reference 
 * process must not be downscaled.
 * 
 * # Basic task configuration
 * 
 * A static Add function performs the basic task configuration. The tasks configures itself
 * automatically based on the run number. Only datasets which were taken in the rare
 * trigger mode are supported. These inlude all periods staring from p-Pb 2013, with the 
 * exception of pp 2015 as in these datasets the trigger was in development but still 
 * enabled. Running the task on unsupported datasets will lead to an AliFatal.
 * 
 * The following example configures the task for any known dataset:
 * 
 * ~~~.{cxx}
 * auto normtask = PWG::EMCAL::AliAnalysisTaskEmcalTriggerNormalization::AddTaskEmcalTriggerNormalization("normatask");
 * ~~~
 * 
 * Attention: The task will skip calculating the luminosity for runs for which the EMCAL
 * was not a trigger detector. Also the normalization histograms and the trigger cluster
 * counters are not filled for these runs. Only the min. bias luminosity (non-downscaled)
 * is counted.
 * 
 * Attention: Task needs access to the OCDB, therefore a the CDB connect task is required
 * to be included in train runs or local / grid analyses.
 * 
 * Attention: The task must run on data only - MC is not supported.
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

  /**
   * @brief Initialize task settings at the beginning of the event loop
   * 
   * Initialization is done based on the run number of the first event. The
   * following datasets are supported:
   * 
   * - 2013 p-Pb at 5.02 TeV (LHC13b-f)
   * - 2015 Pb-Pb at 5.02 TeV (LHC15o)
   * - 2016-2018 pp at 13 TeV (LHC16d-p, LHC17c-o, LHC17r, LHC18b-p)
   * - 2016 p-Pb at 5.02 TeV (LHC16q+t, note: EMCAL triggers were just a subset of min. bias)
   * - 2016 p-Pb at 8.16 TeV (LHC16r-s)
   * - 2017 pp at 5.02 TeV (LHC17q-r, note: EMCAL triggers were running only in CALOFAST cluster)
   * 
   * In case an unsupported dataset is requested the processing stops with AliFatal
   */
  virtual void UserExecOnce(); 

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
   * @brief Check whether run number is from a p-Pb run from 2013
   * @param runnumber Run number to check
   * @return true in case runs from p-Pb 2013 are processed, false otherwise 
   */
  inline bool IsRun1pPb5TeV(int runnumber) const;

  /**
   * @brief Check whether the run number is from a pp run at 13 TeV from 2016-2018
   * @param runnumber Run number to check
   * @return true in case runs from pp 2016-2018 at 13 TeV are processed, false otherwise
   */
  inline bool IsRun2pp13TeV(int runnumber) const;

  /**
   * @brief Check whether the run number is from a pp run at 5.02 TeV from 2017
   * @param runnunber Run number to check
   * @return true in case runs from pp 2017 at 5.02 TeV are processed, false otherwise
   */
  inline bool IsRun2pp5TeV(int runnunber) const;

  /**
   * @brief Check whether the run number is from a p-Pb run at 5.02 TeV from 2016
   * @param runnunber Run number to check
   * @return true in case runs from p-Pb 2016 at 5.02 TeV are processed, false otherwise
   */
  inline bool IsRun2pPb5TeV(int runnumber) const;

  /**
   * @brief Check whether the run number is from a p-Pb run at 8.16 TeV from 2016
   * @param runnunber Run number to check
   * @return true in case runs from p-Pb 2016 at 8.16 TeV are processed, false otherwise
   */
  inline bool IsRun2pPb8TeV(int runnumber) const;

  /**
   * @brief Check whether the run number is from Pb-Pb run at 5.02 TeV from 2015 or 2018
   * @param runnunber Run number to check
   * @return true in case runs from Pb-Pb 2015 or 2018 at 5.02 TeV are processed, false otherwise
   */
  inline bool IsRun2PbPb5TeV(int runnumber) const;

  /**
   * @brief Get the year in which a run was taken
   * @param runnumber Run number to check
   * @return Year in which the run was taken
   */
  inline int getYear(int runnumber) const;


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

bool AliAnalysisTaskEmcalTriggerNormalization::IsRun1pPb5TeV(int runnumber) const {
  return runnumber >= 195344 && runnumber <= 197388;              // LHC13b-f
}

bool AliAnalysisTaskEmcalTriggerNormalization::IsRun2pp13TeV(int runnumber) const {
  return (runnumber >= 252235 && runnumber <= 264347) ||                                                  // LHC16d - LHC16r
         (runnumber >= 270581 && runnumber <= 281961) || (runnumber >= 282528 && runnumber <= 282704) ||  // LHC17c - LHC17o, LHC17r 
         (runnumber >= 285009 && runnumber <= 294925);                                                    // LHC18b - LHC18p
}

bool AliAnalysisTaskEmcalTriggerNormalization::IsRun2pp5TeV(int runnunber) const {
  return (runnunber >= 282008 && runnunber <= 282441);
}

bool AliAnalysisTaskEmcalTriggerNormalization::IsRun2pPb5TeV(int runnumber) const {
  return (runnumber >= 265309 && runnumber <= 265525) ||           // LHC16q
         (runnumber >= 267163 && runnumber <= 267166);             // LHC16t
}

bool AliAnalysisTaskEmcalTriggerNormalization::IsRun2pPb8TeV(int runnumber) const {
  return runnumber >= 265589 && runnumber <= 267131;               // LHC16r-s
}

bool AliAnalysisTaskEmcalTriggerNormalization::IsRun2PbPb5TeV(int runnumber) const {
  return (runnumber >= 244824 && runnumber <= 246994) ||            // LHC15o
         (runnumber >= 295585 && runnumber <= 297624);              // LHC18q+r
}



int AliAnalysisTaskEmcalTriggerNormalization::getYear(int runnumber) const {
  if(runnumber >= 66719  && runnumber <= 105523) return 2009;
  if(runnumber >= 105524 && runnumber <= 139667) return 2010;
  if(runnumber >= 140390 && runnumber <= 170718) return 2011;
  if(runnumber >= 170730 && runnumber <= 194306) return 2012;
  if(runnumber >= 194481 && runnumber <= 199162) return 2013;
  if(runnumber >= 208402 && runnumber <= 247167) return 2015;
  if(runnumber >= 247656 && runnumber <= 267252) return 2016;
  if(runnumber >= 267402 && runnumber <= 282843) return 2017;
  if(runnumber >= 282908 && runnumber <= 297635) return 2018;
  return -1;
}

}

}


#endif
