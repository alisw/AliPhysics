#ifndef ALIANALYSISTASKTRACKSKIM_H
#define ALIANALYSISTASKTRACKSKIM_H

/**
 * @class AliAnalysisTaskTrackSkim
 * @brief Minimal track skim.
 *
 * Track skimming for jet analysis. Borrows from HFTreeCreator.
 *
 * @author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
 * @date 1 August 2021
 */

#include <map>
#include <string>

class TH1;
class TH2;
class TH3;
class TH3F;
class TTree;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAnalysisManager;
class AliJetContainer;
class AliEmcalJetFinder;

#include "AliAnalysisTaskEmcalJet.h"
#include "AliYAMLConfiguration.h"
#include "AliMCAnalysisUtils.h"

// operator<< has to be forward declared carefully to stay in the global namespace so that it works with CINT.
// For generally how to keep the operator in the global namespace, See: https://stackoverflow.com/a/38801633
namespace PWGJE {
namespace EMCALJetTasks {
  class AliAnalysisTaskTrackSkim;
} // namespace EMCALJetTasks
} // namespace PWGJE
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::AliAnalysisTaskTrackSkim& myTask);
void swap(PWGJE::EMCALJetTasks::AliAnalysisTaskTrackSkim& first,
     PWGJE::EMCALJetTasks::AliAnalysisTaskTrackSkim& second);

namespace PWGJE {
namespace EMCALJetTasks {

class AliAnalysisTaskTrackSkim : public AliAnalysisTaskEmcalJet
{
 public:
  enum EncodedInformation_t {
    kNone = 0,                  //!<! No additional encoded information
    kPhotonPrompt = 1,          //!<! Prompt photons (only meaningful at MC level)
    kPhotonFragmentation = 2    //!<! Fragmentation photons (only meaningful at MC level)
  };
  static const std::map<std::string, AliMCAnalysisUtils::generator> fgkMCUtilsGeneratorMap; //!<! Map from name to MC utils generator used with the YAML config

  AliAnalysisTaskTrackSkim();
  AliAnalysisTaskTrackSkim(const char* name);
  // Additional constructors
  AliAnalysisTaskTrackSkim(const AliAnalysisTaskTrackSkim & other);
  AliAnalysisTaskTrackSkim& operator=(AliAnalysisTaskTrackSkim other);
  friend void ::swap(AliAnalysisTaskTrackSkim & first, AliAnalysisTaskTrackSkim & second);
  // Avoid implementing move since c++11 is not allowed in the header
  virtual ~AliAnalysisTaskTrackSkim() {}

  void UserCreateOutputObjects();

  // Initialize the task
  // Configuration is handled via the YAML configuration file
  bool Initialize();
  void AddConfigurationFile(const std::string & configurationPath, const std::string & configName = "") { fYAMLConfig.AddConfiguration(configurationPath, configName); }

  static AliAnalysisTaskTrackSkim* AddTaskTrackSkim(const std::string& suffix = "");

  // Printing
  std::string toString() const;
  friend std::ostream & ::operator<<(std::ostream &in, const AliAnalysisTaskTrackSkim &myTask);
  void Print(Option_t* opt = "") const;
  std::ostream & Print(std::ostream &in) const;

 protected:
  Bool_t IsEventSelected();
  Bool_t Run();
  Bool_t FillHistograms();

  void RetrieveAndSetTaskPropertiesFromYAMLConfig();
  void AddParticleKinematicsToMap(std::map<std::string, std::vector<float>> & map);
  void AddParticleLabelsToMap(std::map<std::string, std::vector<int>> & map);
  void SetupTree();
  void ClearTree();

  // Helpers
  std::string GetMCHeaderName();
  // Implemented as a singleton for convenience (approach borrowed from Florian)
  AliMCAnalysisUtils * GetMCAnalysisUtils();
  void EncodeAdditionalParticleLevelInformation(AliVParticle * particle);

  // Basic configuration
  PWG::Tools::AliYAMLConfiguration fYAMLConfig; ///<  YAML configuration file.
  bool fConfigurationInitialized;               ///<  True if the task configuration has been successfully initialized.
  bool fChargeScaling;                          ///<  Flag explicitly whether charge scaling is enabled.
  bool fEncodeAdditionalParticleInformation;    ///<  Flag whether to encode additional particle information.
  // Analysis helpers
  AliMCAnalysisUtils::generator fMCAnalysisUtilsGenerator;  ///<  MC analysis tools generator
  AliMCAnalysisUtils* fMCAnalysisUtils;                     //!<!  MC analysis tools

  TH1F* fNAcceptedFirstTrackCollection;         //!<! Keep track of number of tracks in the first collection (which also tracks events with zero accepted tracks)

  // Tree variables
  // Event level variables
  bool fTriggerBitINT7;                       //!<! Flag explicitly whether trigger bitmap contains INT7
  float fCentralityForTree;                   //!<! Centrality (as a float to reduce storage requirements).
  float fEventPlaneV0MForTree;                //!<! V0M event plane (as a float to reduce storage requirements).
  bool fTriggerBitCentral;                    //!<! Flag explicitly whether trigger bitmap contains kCentral
  bool fTriggerBitSemiCentral;                //!<! Flag explicitly whether trigger bitmap contains kSemiCentral
  // Track level properties
  // We need at most two collections:
  std::map<std::string, std::vector<float>> fFirstTrackCollectionKinematics;   //!<! Kinematics for first track collection (stored in floats).
  std::map<std::string, std::vector<int>> fFirstTrackCollectionLabels;         //!<! Labels for first track collection (optional - only used for MC).
  std::map<std::string, std::vector<float>> fSecondTrackCollectionKinematics;  //!<! Kinematics for second track collection (stored in floats).
  std::map<std::string, std::vector<int>> fSecondTrackCollectionLabels;        //!<! Labels for second track collection (optional - only used for MC).
  // The tree itself.
  TTree* fTreeSkim;                           //!<! Tree containing the track skim.

  // And version the output, so we can store it in the tree name
  const static unsigned int fOutputVersion = 3;     ///<  Version of the track skim output. Be sure to increment with any meaningful changes!
 private:
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskTrackSkim, 3)  // Track skim
  /// \endcond
};

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

#endif  /* ALIANALYSISTASKTRACKSKIM_H */
