//
// Utilities class for Jet-Hadron correlation analysis
//

#include "AliAnalysisTaskEmcalJetHUtils.h"

// Require to use AliLog streams with some types
#include <iostream>

#include <TMath.h>

#include <AliLog.h>

#include "AliEmcalJet.h"
#include "AliEmcalContainerUtils.h"
#include "AliEmcalContainer.h"
#include "AliParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalParticleJetConstituent.h"
#include "AliEmcalClusterJetConstituent.h"

namespace PWGJE {
namespace EMCALJetTasks {

const std::map<std::string, AliAnalysisTaskEmcalJetHUtils::ELeadingHadronBiasType_t> AliAnalysisTaskEmcalJetHUtils::fgkLeadingHadronBiasMap = {
  { "kCharged", AliAnalysisTaskEmcalJetHUtils::kCharged},
  { "kNeutral", AliAnalysisTaskEmcalJetHUtils::kNeutral},
  { "kBoth", AliAnalysisTaskEmcalJetHUtils::kBoth}
};

/**
 * Determine leading hadron pt in a jet. This is inspired by AliJetContainer::GetLeadingHadronMomentum(), but
 * that particular function is avoided because the cluster energy retrieved is always the raw E while the
 * cluster energy used in creating the jet would be preferred. One could create a cluster container and go
 * through all of those steps, but there is a simpler approach: the leading charged and neutral momenta
 * are stored in AliEmcalJet while performing jet finding.
 *
 * @param[in] jet Jet from which the leading hadron pt should be extracted
 * @param[in] leadingHadronType Type of leading hadron pt to retrieve
 *
 * @return Value of the leading hadron pt
 */
double AliAnalysisTaskEmcalJetHUtils::GetLeadingHadronPt(AliEmcalJet * jet, AliAnalysisTaskEmcalJetHUtils::ELeadingHadronBiasType_t leadingHadronType)
{
  double maxTrackPt = 0;
  double maxClusterPt = 0;

  if (leadingHadronType == kCharged || leadingHadronType == kBoth) {
    auto particle = jet->GetLeadingParticleConstituent();
    if (particle) {
      maxTrackPt = particle->Pt();
    }
  }
  if (leadingHadronType == kNeutral || leadingHadronType == kBoth) {
    // NOTE: We don't want to use jet->MaxNeutralPt() because this uses energy
    //       from the neutral particles at the particle level. While this is not
    //       strictly wrong, it can be rather misleading to have a leading neutral
    //       particle value when we are really interested in the cluster pt that is
    //       only meaningful at detector level.
    auto cluster = jet->GetLeadingClusterConstituent();
    if (cluster) {
      // Uses the energy definition that was used when the constituent was created
      // to calculate the Pt(). Usually, this would be the hadronic corrected energy
      maxClusterPt = cluster->Pt();
    }
  }

  // The max value will be 0 unless it was filled. Thus, it will only be greater if
  // it was requested.
  return (maxTrackPt > maxClusterPt) ? maxTrackPt : maxClusterPt;
}

/**
 * Function to calculate angle between jet and EP in the 1st quadrant (0,Pi/2).
 * Adapted from AliAnalysisTaskEmcalJetHadEPpid.
 *
 * @param jetAngle Phi angle of the jet (could be any particle)
 * @param epAngle Event plane angle
 *
 * @return Angle between jet and EP in the 1st quadrant (0,Pi/2)
 */
double AliAnalysisTaskEmcalJetHUtils::RelativeEPAngle(double jetAngle, double epAngle)
{
  double dphi = (epAngle - jetAngle);

  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi<-1*TMath::Pi() ) {
    dphi = dphi + 1*TMath::Pi();
  } // this assumes we are doing full jets currently

  if( (dphi>0) && (dphi<1*TMath::Pi()/2) ) {
    // Do nothing! we are in quadrant 1
  } else if ( (dphi>1*TMath::Pi()/2) && (dphi<1*TMath::Pi()) ) {
    dphi = 1*TMath::Pi() - dphi;
  } else if ( (dphi<0) && (dphi>-1*TMath::Pi()/2) ) {
    dphi = std::abs(dphi);
  } else if ( (dphi<-1*TMath::Pi()/2) && (dphi>-1*TMath::Pi()) ) {
    dphi = dphi + 1*TMath::Pi();
  }

  // Warn if we are not in the proper range
  if ( dphi < 0 || dphi > TMath::Pi()/2 ) {
    AliWarningGeneralStream("AliAnalysisTaskEmcalJetHUtils") << ": dPHI not in range [0, 0.5*Pi]!\n";
  }

  return dphi;   // dphi in [0, Pi/2]
}

/**
 * Configure an AliEventCuts object with the options in the given AliYAMLConfiguration object and the task trigger mask.
 *
 * @param[in] eventCuts AliEventCuts object to configure.
 * @param[in] yamlConfig %YAML configuration object to be used in configuring the event cuts object.
 * @param[in] offlineTriggerMask Trigger mask (set via SelectCollisionCandidates()) from the task. The value can be updated by values in the YAML config.
 * @param[in] baseName Name under which the settings should be looked for in the %YAML config.
 * @param[in] taskName Name of the analysis task for which this function was called. This is to make it clear which task is being configured.
 */
void AliAnalysisTaskEmcalJetHUtils::ConfigureEventCuts(AliEventCuts & eventCuts, PWG::Tools::AliYAMLConfiguration & yamlConfig, const UInt_t offlineTriggerMask, const std::string & baseName, const std::string & taskName)
{
  // The trigger can be set regardless of event cuts settings.
  // Event cuts trigger selection.
  bool useEventCutsAutomaticTriggerSelection = false;
  bool res = yamlConfig.GetProperty(std::vector<std::string>({baseName, "useAutomaticTriggerSelection"}), useEventCutsAutomaticTriggerSelection, false);
  if (res && useEventCutsAutomaticTriggerSelection) {
    // Use the automatic selection. Nothing to be done.
    AliInfoGeneralStream(taskName.c_str()) << "Using the automatic trigger selection from AliEventCuts.\n";
  }
  else {
    // Use the cuts selected by SelectCollisionCandidates() (or via YAML)
    AliInfoGeneralStream(taskName.c_str()) << "Using the trigger selection specified with SelectCollisionCandidates() or via YAML. Value: " << offlineTriggerMask << "\n";
    eventCuts.OverrideAutomaticTriggerSelection(offlineTriggerMask, true);
  }

  // Manual mode
  bool manualMode = false;
  yamlConfig.GetProperty({baseName, "manualMode"}, manualMode, false);
  if (manualMode) {
    AliInfoGeneralStream(taskName.c_str()) << "Configuring manual event cuts.\n";
    eventCuts.SetManualMode();
    // Confgure manual mode via YAML
    // Select the period
    typedef void (AliEventCuts::*MFP)();
    std::map<std::string, MFP> eventCutsPeriods = { std::make_pair("LHC11h", &AliEventCuts::SetupRun1PbPb),
                            std::make_pair("LHC15o", &AliEventCuts::SetupLHC15o) };
    std::string manualCutsPeriod = "";
    yamlConfig.GetProperty({ baseName, "cutsPeriod" }, manualCutsPeriod, true);
    auto eventCutsPeriod = eventCutsPeriods.find(manualCutsPeriod);
    if (eventCutsPeriod != eventCutsPeriods.end()) {
      // Call event cuts period setup.
      (eventCuts.*eventCutsPeriod->second)();
      AliDebugGeneralStream(taskName.c_str(), 3) << "Configuring event cuts with period \"" << manualCutsPeriod << "\"\n";
    } else {
      AliFatalGeneralF(taskName.c_str(), "Period %s was not found in the event cuts period map.", manualCutsPeriod.c_str());
    }

    // Additional settings must be after setting the period to ensure that the settings aren't overwritten.
    // Centrality
    std::pair<double, double> centRange;
    res = yamlConfig.GetProperty("centralityRange", centRange, false);
    if (res) {
      AliDebugGeneralStream(taskName.c_str(), 3) << "Setting centrality range of (" << centRange.first << ", " << centRange.second << ").\n";
      eventCuts.SetCentralityRange(centRange.first, centRange.second);
    }

    // MC
    bool mc = false;
    yamlConfig.GetProperty({ baseName, "MC" }, mc, false);
    eventCuts.fMC = mc;

    // Set the 15o pileup cuts. Defaults to on.
    if (manualCutsPeriod == "LHC15o") {
      bool enablePileupCuts = true;
      yamlConfig.GetProperty({ baseName, "enablePileupCuts" }, enablePileupCuts, false);
      AliDebugGeneralStream(taskName.c_str(), 3) << "Setting 15o pileup cuts to " << std::boolalpha << enablePileupCuts << ".\n";
      eventCuts.fUseVariablesCorrelationCuts = enablePileupCuts;
    }
  }
}

/**
 * Utility function to create a particle or track container given the collection name of the desired container.
 *
 * @param[in] collectionName Name of the particle or track collection name.
 *
 * @return A newly created particle or track container.
 */
AliParticleContainer * AliAnalysisTaskEmcalJetHUtils::CreateParticleOrTrackContainer(const std::string & collectionName)
{
  AliParticleContainer * partCont = nullptr;
  if (collectionName == AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kTrack)) {
    AliTrackContainer * trackCont = new AliTrackContainer(collectionName.c_str());
    partCont = trackCont;
  }
  else if (collectionName != "") {
    partCont = new AliParticleContainer(collectionName.c_str());
  }

  return partCont;
}

/**
 * Configure an EMCal container according to the specified YAML configuration.
 *
 * @param[in] baseName Name under which the config is stored, with the container name included.
 * @param[in] containerName Name of the container.
 * @param[in] cont Existing particle container.
 * @param[in] yamlConfig YAML configuration to be used.
 * @param[in] taskName Name of the task which is calling this function (for debugging purposes).
 */
void AliAnalysisTaskEmcalJetHUtils::ConfigureEMCalContainersFromYAMLConfig(std::vector<std::string> baseName,
                                      std::string containerName,
                                      AliEmcalContainer* cont,
                                      PWG::Tools::AliYAMLConfiguration& yamlConfig,
                                      std::string taskName)
{
  // Initial setup
  cont->SetName(containerName.c_str());

  // Set the properties
  double tempDouble = -1.0;
  bool tempBool = false;
  std::vector<double> tempRange;
  // Min Pt
  bool result = yamlConfig.GetProperty(baseName, "minPt", tempDouble, false);
  if (result) {
    AliDebugGeneralStream(taskName.c_str(), 2) << cont->GetName() << ": Setting minPt of " << tempDouble << "\n";
    cont->SetMinPt(tempDouble);
  }
  // Min E
  result = yamlConfig.GetProperty(baseName, "minE", tempDouble, false);
  if (result) {
    AliDebugGeneralStream(taskName.c_str(), 2) << cont->GetName() << ": Setting minE of " << tempDouble << "\n";
    cont->SetMinE(tempDouble);
  }
  // Eta min, max
  result = yamlConfig.GetProperty(baseName, "etaLimits", tempRange, false);
  if (result) {
    if (tempRange.size() != 2) {
      AliErrorGeneralStream(taskName.c_str()) << "Passed eta range with " << tempRange.size()
                          << " entries, but 2 values are required. Ignoring values.\n";
    } else {
      AliDebugGeneralStream(taskName.c_str(), 2)
       << "Setting eta range to [" << tempRange.at(0) << ", " << tempRange.at(1) << "]\n";
      cont->SetEtaLimits(tempRange.at(0), tempRange.at(1));
    }
  }
  // Phi min, max
  result = yamlConfig.GetProperty(baseName, "phiLimits", tempRange, false);
  if (result) {
    if (tempRange.size() != 2) {
      AliErrorGeneralStream(taskName.c_str()) << "Passed phi range with " << tempRange.size()
                          << " entries, but 2 values are required. Ignoring values.\n";
    } else {
      AliDebugGeneralStream(taskName.c_str(), 2)
       << "Setting phi range to [" << tempRange.at(0) << ", " << tempRange.at(1) << "]\n";
      cont->SetPhiLimits(tempRange.at(0), tempRange.at(1));
    }
  }
  // Embedded
  result = yamlConfig.GetProperty(baseName, "embedding", tempBool, false);
  if (result) {
    AliDebugGeneralStream(taskName.c_str(), 2)
     << cont->GetName() << ": Setting embedding to " << (tempBool ? "enabled" : "disabled") << "\n";
    cont->SetIsEmbedding(tempBool);
  }
}

/**
 * Configure a track container according to the specified YAML configuration.
 *
 * @param[in] baseNameWithhContainer Name under which the config is stored, with the container name included.
 * @param[in] trackCont Existing particle container.
 * @param[in] yamlConfig YAML configuration to be used.
 * @param[in] taskName Name of the task which is calling this function (for debugging purposes).
 */
void AliAnalysisTaskEmcalJetHUtils::ConfigureTrackContainersFromYAMLConfig(
 std::vector<std::string> baseNameWithContainer, AliTrackContainer* trackCont,
 PWG::Tools::AliYAMLConfiguration& yamlConfig, std::string taskName)
{
  // Initial setup
  std::string tempString = "";

  // Track selection
  // AOD Filter bits as a sequence
  std::vector<UInt_t> filterBitsVector;
  bool result = yamlConfig.GetProperty(baseNameWithContainer, "aodFilterBits", filterBitsVector, false);
  if (result) {
    UInt_t filterBits = 0;
    for (int filterBit : filterBitsVector) {
      filterBits += filterBit;
    }
    AliDebugGeneralStream(taskName.c_str(), 2)
     << trackCont->GetName() << ": Setting filterBits of " << filterBits << std::endl;
    trackCont->SetAODFilterBits(filterBits);
  }

  // SetTrackFilterType enum
  result = yamlConfig.GetProperty(baseNameWithContainer, "trackFilterType", tempString, false);
  if (result) {
    // Need to get the enumeration
    AliEmcalTrackSelection::ETrackFilterType_t trackFilterType =
     AliTrackContainer::fgkTrackFilterTypeMap.at(tempString);
    AliDebugGeneralStream(taskName.c_str(), 2)
     << trackCont->GetName() << ": Setting trackFilterType of " << trackFilterType << " (" << tempString << ")\n";
    trackCont->SetTrackFilterType(trackFilterType);
  }

  // Track cuts period
  result = yamlConfig.GetProperty(baseNameWithContainer, "trackCutsPeriod", tempString, false);
  if (result) {
    // Need to get the enumeration
    AliDebugGeneralStream(taskName.c_str(), 2)
     << trackCont->GetName() << ": Setting track cuts period to " << tempString << std::endl;
    trackCont->SetTrackCutsPeriod(tempString.c_str());
  }
}

/**
 * Configure a cluster container according to the specified YAML configuration.
 *
 * @param[in] baseNameWithhContainer Name under which the config is stored, with the container name included.
 * @param[in] clusterCont Existing cluster container.
 * @param[in] yamlConfig YAML configuration to be used.
 * @param[in] taskName Name of the task which is calling this function (for debugging purposes).
 */
void AliAnalysisTaskEmcalJetHUtils::ConfigureClusterContainersFromYAMLConfig(
 std::vector<std::string> baseNameWithContainer, AliClusterContainer* clusterCont,
 PWG::Tools::AliYAMLConfiguration& yamlConfig, std::string taskName)
{
  // Initial setup
  double tempDouble = 0;
  std::string tempString = "";
  bool tempBool = false;

  // Default energy
  bool result = yamlConfig.GetProperty(baseNameWithContainer, "defaultClusterEnergy", tempString, false);
  if (result) {
    // Need to get the enumeration
    AliVCluster::VCluUserDefEnergy_t clusterEnergyType =
     AliClusterContainer::fgkClusterEnergyTypeMap.at(tempString);
    AliDebugGeneralStream(taskName.c_str(), 2)
     << clusterCont->GetName() << ": Setting cluster energy type to " << clusterEnergyType << std::endl;
    clusterCont->SetDefaultClusterEnergy(clusterEnergyType);
  }

  // NonLinCorrEnergyCut
  result = yamlConfig.GetProperty(baseNameWithContainer, "clusNonLinCorrEnergyCut", tempDouble, false);
  if (result) {
    AliDebugGeneralStream(taskName.c_str(), 2)
     << clusterCont->GetName() << ": Setting clusNonLinCorrEnergyCut of " << tempDouble << std::endl;
    clusterCont->SetClusNonLinCorrEnergyCut(tempDouble);
  }

  // HadCorrEnergyCut
  result = yamlConfig.GetProperty(baseNameWithContainer, "clusHadCorrEnergyCut", tempDouble, false);
  if (result) {
    AliDebugGeneralStream(taskName.c_str(), 2)
     << clusterCont->GetName() << ": Setting clusHadCorrEnergyCut of " << tempDouble << std::endl;
    clusterCont->SetClusHadCorrEnergyCut(tempDouble);
  }

  // SetIncludePHOS
  result = yamlConfig.GetProperty(baseNameWithContainer, "includePHOS", tempBool, false);
  if (result) {
    AliDebugGeneralStream(taskName.c_str(), 2) << clusterCont->GetName() << ": Setting Include PHOS to "
                          << (tempBool ? "enabled" : "disabled") << std::endl;
    clusterCont->SetIncludePHOS(tempBool);
  }
}

/**
 * Get the background subtracted jet pt.
 *
 * @param[in] jet Jet to be subtracted.
 * @param[in] rho Rho value for the jet collection.
 */
double AliAnalysisTaskEmcalJetHUtils::GetJetPt(const AliEmcalJet* jet, const double rho)
{
  double pT = jet->Pt() - rho * jet->Area();
  return pT;
}

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */
