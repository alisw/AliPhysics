#ifndef ALIANALYSISTASKEMCALJETHPERFORMANCE_H
#define ALIANALYSISTASKEMCALJETHPERFORMANCE_H

/**
 * @class AliAnalysisTaskEmcalJetHPerformance
 * @brief Jet-hadron correlations task dedicated to performance of the analysis
 *
 * @author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
 * @date 23 Feb 2018
 */

#include <string>

class AliJetContainer;
class AliEmcalJet;
class AliVCaloCells;
class AliQnCorrectionsManager;
#include "THistManager.h"
#include "AliYAMLConfiguration.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliEmcalEmbeddingQA.h"
#include "AliAnalysisTaskEmcalJetHUtils.h"

// operator<< has to be forward declared carefully to stay in the global namespace so that it works with CINT.
// For generally how to keep the operator in the global namespace, See: https://stackoverflow.com/a/38801633
namespace PWGJE { namespace EMCALJetTasks { class AliAnalysisTaskEmcalJetHPerformance; } }
std::ostream & operator<< (std::ostream &in, const PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance &myTask);
void swap(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance & first, PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance & second);

namespace PWGJE {
namespace EMCALJetTasks {

class AliAnalysisTaskEmcalJetHPerformance : public AliAnalysisTaskEmcalJet {
 public:
  /**
   * Wrapper to contain the properties for a jet to simplify filling the response matrix histogram
   */
  struct ResponseMatrixFillWrapper {
    double fPt;                 //!<! Jet pt
    double fArea;               //!<! Jet area
    double fPhi;                //!<! jet phi
    double fRelativeEPAngle;    //!<! Angle relative to the event plane
    double fLeadingHadronPt;    //!<! Leading hadron pt
    double fDistance;           //!<! Distance to the matched jet
    double fCentrality;         //!<! Centrality of the given jet (an event level property, but useful to here available here)
  };

  AliAnalysisTaskEmcalJetHPerformance();
  AliAnalysisTaskEmcalJetHPerformance(const char * name);
  // Additional constructors
  AliAnalysisTaskEmcalJetHPerformance(const AliAnalysisTaskEmcalJetHPerformance & other);
  AliAnalysisTaskEmcalJetHPerformance& operator=(AliAnalysisTaskEmcalJetHPerformance other);
  friend void ::swap(AliAnalysisTaskEmcalJetHPerformance & first, AliAnalysisTaskEmcalJetHPerformance & second);
  // Avoid implementing move since c++11 is not allowed in the header

  void UserCreateOutputObjects();

  // Initialize the task
  // Configuration is handled via the YAML configuration file
  bool Initialize();
  void AddConfigurationFile(const std::string & configurationPath, const std::string & configName = "") { fYAMLConfig.AddConfiguration(configurationPath, configName); }

  // Utility functions
  // AddTask
  static AliAnalysisTaskEmcalJetHPerformance * AddTaskEmcalJetHPerformance(const char * suffix = "");
  // Printing
  std::string toString() const;
  friend std::ostream & ::operator<<(std::ostream &in, const AliAnalysisTaskEmcalJetHPerformance &myTask);
  void Print(Option_t* opt = "") const;
  std::ostream & Print(std::ostream &in) const;

 private:

  Bool_t Run();

  // Helper functions
  double DetermineTrackingEfficiency(double trackPt, double trackEta);

  // Configuration
  void RetrieveAndSetTaskPropertiesFromYAMLConfig();
  void SetupJetContainersFromYAMLConfig();
  void SetupParticleContainersFromYAMLConfig();
  void SetupClusterContainersFromYAMLConfig();

  // QA histograms
  void SetupQAHists();
  void QAHists();
  void FillQAHists();
  // Cell QA
  void SetupCellQAHistsWithPrefix(const std::string & prefix);
  void SetupCellQAHists();
  void FillCellQAHists(const std::string & prefix, AliVCaloCells * cells);
  void FillCellQAHists();

  // Jet matching
  void SetupJetMatchingQA();
  void ResetMatching(const AliJetContainer & jetCont) const;
  bool PerformGeometricalJetMatching(AliJetContainer& contBase, AliJetContainer& contTag, double maxDist) const;
  void FillJetMatchingQA(AliJetContainer& contBase, AliJetContainer& contTag, const std::string& prefix);

  // Response matrix functions
  void SetupResponseMatrixHists();
  void ResponseMatrix();
  void FillResponseMatrix(AliEmcalJet* jet1, AliEmcalJet* jet2, const double jet1Rho);
  ResponseMatrixFillWrapper CreateResponseMatrixFillWrapper(AliEmcalJet * jet, const double rho) const;

  // Basic configuration
  PWG::Tools::AliYAMLConfiguration fYAMLConfig; ///< YAML configuration file.
  bool fConfigurationInitialized;     ///<  True if the task configuration has been successfully initialized.

  // Histograms
  THistManager fHistManager;          ///<  Histogram manager.
  AliEmcalEmbeddingQA fEmbeddingQA;   //!<! Embedding QA hists (will only be added if embedding).

  // Configuration options
  bool fCreateQAHists;                ///<  If true, create QA histograms.
  bool fCreateCellQAHists;            ///<  If true, create the Cell QA histograms. It doesn't gracefully turn off when not configured like the containers, so we have a switch.
  bool fPerformJetMatching;           ///<  If true, enables jet matching.
  bool fCreateResponseMatrix;         ///<  If true, create a response matrix with the available jet collections.

  // QA variables
  std::string fEmbeddedCellsName;                 ///<  Set the embedded cells collection name
  UInt_t fPreviousEventTrigger;                   ///<  Physics selection (offline trigger) of the previous event for determine why a small number of embedded event are double counted.
  bool fPreviousEmbeddedEventSelected;            ///<  True if the previous embedded event was selected. Used to determine why a small number of embedded event are double counted.
  AliAnalysisTaskEmcalJetHUtils::EEfficiencyPeriodIdentifier_t fEfficiencyPeriodIdentifier;  ///<  Identifies the period for determining the efficiency correction to apply
  AliQnCorrectionsManager *fFlowQnVectorManager;  //!<! Qn corrections framework manager.

  // Jet matching
  double fMaxJetMatchingDistance;                 ///<  Matx jet matching distance.

  // Response matrix variables
  // Response matrix fill map
  // For whatever reason, it appears that CINT can't handle this definition
  #if !(defined(__CINT__) || defined(__MAKECINT__))
  // For some reason, root can't handle this typedef
  //typedef double ResponseMatrixFillWrapper::*MP;
  //std::map <std::string, std::pair<int, MP>> fResponseMatrixFillMap;   //!<!
  std::map <std::string, std::pair<int, double ResponseMatrixFillWrapper::*>> fResponseMatrixFillMap;   //!<! Map from axis title to pair of (jet number, function to retrieve fill value)
  #endif
  bool fResponseFromThreeJetCollections; ///<  If true, the det level jets in collection 2 are only an intermediate step. They are used to get part level jets to match to hybrid jets
  double fMinFractionShared;             ///<  Minimum fraction of shared jet pt required for matching a hybrid jet to detector level
  AliAnalysisTaskEmcalJetHUtils::ELeadingHadronBiasType_t fLeadingHadronBiasType; ///<  Leading hadron in jet bias type (either charged, neutral, or both)

  ClassDef(AliAnalysisTaskEmcalJetHPerformance, 7);
};

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

#endif /* AliAnalysisTaskEmcalJetHPerformance.h */
