#ifndef ALIANALYSISTASKEMCALJETHUTILS_H
#define ALIANALYSISTASKEMCALJETHUTILS_H

/**
 * @class AliAnalysisTaskEmcalJetHUtils
 * @brief Jet-hadron correlations utilities class
 *
 * Contains functionality that is shared between the various classes. Could have been
 * a namespace except it wouldn't play nice with ROOT.
 *
 * @author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
 * @date 23 Feb 2018
 */

#include <string>
#include <map>

#include "AliEventCuts.h"
#include "AliYAMLConfiguration.h"
// NOTE: AliAnalysisTaskSE is just needed for AliQnCorrectionsVarManagerTask...
#include "AliAnalysisTaskSE.h"
#include "AliQnCorrectionsVarManagerTask.h"
#include "AliEmcalJet.h"
class AliEmcalContainer;
class AliParticleContainer;
class AliTrackContainer;
class AliClusterContainer;
class AliAnalysisTaskFlowVectorCorrections;
class AliQnCorrectionsManager;
class AliQnCorrectionsHistos;

namespace PWGJE {
namespace EMCALJetTasks {

class AliAnalysisTaskEmcalJetHUtils {
 public:
  /**
   * @enum ELeadingHadronBiasType_t
   * @brief Determine the jet leading hadron bias type
   */
  enum ELeadingHadronBiasType_t {
    kCharged = 0,   //!<! Charged hadron bias
    kNeutral = 1,   //!<! Neutral hadron bias
    kBoth = 2       //!<! Leading is from the leader of both
  };
  static const std::map<std::string, ELeadingHadronBiasType_t> fgkLeadingHadronBiasMap; //!<! Map from name to leading hadron bias used with the YAML config
  static double GetLeadingHadronPt(AliEmcalJet * jet, ELeadingHadronBiasType_t leadingHadronType);

  static double RelativeEPAngle(double jetAngle, double epAngle);

  static void ConfigureEventCuts(AliEventCuts & eventCuts, PWG::Tools::AliYAMLConfiguration & yamlConfig, const UInt_t offlineTriggerMask, const std::string & baseName, const std::string & taskName);
  static AliParticleContainer * CreateParticleOrTrackContainer(const std::string & collectionName);
  static void ConfigureEMCalContainersFromYAMLConfig(std::vector<std::string> baseName,
                            std::string containerName,
                            AliEmcalContainer* cont,
                            PWG::Tools::AliYAMLConfiguration& yamlConfig,
                            std::string taskName);
  static void ConfigureTrackContainersFromYAMLConfig(std::vector<std::string> baseNameWithContainer,
                            AliTrackContainer* trackCont,
                            PWG::Tools::AliYAMLConfiguration& yamlConfig,
                            std::string taskName);
  static void ConfigureClusterContainersFromYAMLConfig(std::vector<std::string> baseNameWithContainer,
                             AliClusterContainer* clusterCont,
                             PWG::Tools::AliYAMLConfiguration& yamlConfig,
                             std::string taskName);
  // Determine jet acceptance from YAML
  static const std::map<std::string, AliEmcalJet::JetAcceptanceType> fgkJetAcceptanceMap;
  static UInt_t DetermineJetAcceptanceFromYAML(const std::vector<std::string> & selections);

  // AddTask for Qn flow vector corrections
  static AliAnalysisTaskFlowVectorCorrections * AddTaskFlowQnVectorCorrections(const std::string & configFilename);

  static double GetJetPt(const AliEmcalJet * jet, const double rho);

  /**
   * @enum EEfficiencyPeriodIdentifier_t
   * @brief Identify the beam type and period that is being analyzed.
   */
  enum EEfficiencyPeriodIdentifier_t {
    kDisableEff = -1,   //!<! No identifier
    kLHC11h = 0,        //!<! Run1 PbPb - LHC11h
    kLHC15o = 1,        //!<! Run2 PbPb - LHC15o
    kLHC18qr = 2,       //!<! Run2 PbPb - LHC18{q,r}
    kLHC11a = 3,        //!<! Run1 pp - LHC11a (2.76 TeV)
    kpA = 4,            //!<! Generic pA
    kpp = 5             //!<! Generic pp
  };
  static const std::map<std::string, EEfficiencyPeriodIdentifier_t> fgkEfficiencyPeriodIdentifier; //!<! Map from name to efficiency period identifier for use with the YAML config
  static double DetermineTrackingEfficiency(const double trackPt, const double trackEta, const int centralityBin,
                       const EEfficiencyPeriodIdentifier_t efficiencyPeriodIdentifier,
                       const std::string& taskName);

 protected:
  // Private interface for tracking efficiency functions
  // LHC11a
  static double LHC11aTrackingEfficiency(const double trackPt, const double trackEta, const int centralityBin,
                      const std::string& taskName);
  // LHC11h
  static double LHC11hTrackingEfficiency(const double trackPt, const double trackEta, const int centralityBin,
                      const std::string& taskName);
  // LHC15o
  static double LHC15oTrackingEfficiency(const double trackPt, const double trackEta, const int centralityBin,
                      const std::string& taskName);

  // LHC11a efficiency parameters are stored directly in the function
  // LHC11h Efficiency parameters
  static const double LHC11hParam_0_10[17];                     //!<! 0-10% good runs
  static const double LHC11hParam_10_30[17];                    //!<! 10-30% good runs
  static const double LHC11hParam_30_50[17];                    //!<! 30-50% good runs
  static const double LHC11hParam_50_90[17];                    //!<! 50-90% good runs
  // LHC15o efficiency parameters
  // pt parameters
  static const double LHC15oParam_0_10_pt[10];                    //!<! 0-10% pt parameters
  static const double LHC15oParam_10_30_pt[10];                   //!<! 10-30% pt parameters
  static const double LHC15oParam_30_50_pt[10];                   //!<! 30-50% pt parameters
  static const double LHC15oParam_50_90_pt[10];                   //!<! 50-90% pt parameters
  // Eta parameters
  static const double LHC15oParam_0_10_eta[13];                   //!<! 0-10% eta parameters
  static const double LHC15oParam_10_30_eta[13];                  //!<! 10-30% eta parameters
  static const double LHC15oParam_30_50_eta[13];                  //!<! 30-50% eta parameters
  static const double LHC15oParam_50_90_eta[13];                  //!<! 50-90% eta parameters
  // Helper functions for determining the LHC15o tracking efficiency
  static double LHC15oPtEfficiency(const double trackPt, const double params[10]);
  static double LHC15oLowPtEfficiencyImpl(const double trackPt, const double params[10], const int index);
  static double LHC15oHighPtEfficiencyImpl(const double trackPt, const double params[10], const int index);
  static double LHC15oEtaEfficiency(const double trackEta, const double params[13]);
  static double LHC15oEtaEfficiencyImpl(const double trackEta, const double params[13], const int index);

};

/**
 * All functions in the FlowVectorCorrections static class are copied diretly from the AddTask. I wish I didn't have
 * to do this, but I need it in compiled code, so I have no choice.
 */
class FlowVectorCorrections {
 public:
  static void AddVZERO(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager, AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity);
  static void AddTPC(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager, AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity);
  static void AddSPD(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager, AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity);
  static void AddTZERO(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager, AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity);
  static void AddZDC(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager, AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity);
  static void AddFMD(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager, AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity);
  static void AddRawFMD(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager, AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity);
  static void DefineHistograms(AliQnCorrectionsManager* QnManager, AliQnCorrectionsHistos* histos, TString histClass);
};

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

#endif /* AliAnalysisTaskEmcalJetHUtils.h */
