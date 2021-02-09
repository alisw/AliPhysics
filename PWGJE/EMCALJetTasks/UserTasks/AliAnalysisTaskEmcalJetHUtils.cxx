//
// Utilities class for Jet-Hadron correlation analysis
//

#include "AliAnalysisTaskEmcalJetHUtils.h"

// Require to use AliLog streams with some types
#include <iostream>
#include <cmath>

#include <TObjArray.h>
#include <TObjString.h>
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

// Flow vector corrections
#include "AliQnCorrectionsProfileCorrelationComponents.h"
#include "AliQnCorrectionsProfile3DCorrelations.h"
#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsCutWithin.h"
#include "AliQnCorrectionsDataVector.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliQnCorrectionsDetector.h"
#include "AliQnCorrectionsDetectorConfigurationTracks.h"
#include "AliQnCorrectionsDetectorConfigurationChannels.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsInputGainEqualization.h"
#include "AliQnCorrectionsQnVectorRecentering.h"
#include "AliQnCorrectionsQnVectorAlignment.h"
#include "AliQnCorrectionsQnVectorTwistAndRescale.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"

namespace PWGJE {
namespace EMCALJetTasks {

const std::map<std::string, AliAnalysisTaskEmcalJetHUtils::ELeadingHadronBiasType_t> AliAnalysisTaskEmcalJetHUtils::fgkLeadingHadronBiasMap = {
  { "kCharged", AliAnalysisTaskEmcalJetHUtils::kCharged},
  { "kNeutral", AliAnalysisTaskEmcalJetHUtils::kNeutral},
  { "kBoth", AliAnalysisTaskEmcalJetHUtils::kBoth}
};

const std::map<std::string, AliEmcalJet::JetAcceptanceType> AliAnalysisTaskEmcalJetHUtils::fgkJetAcceptanceMap = {
  {"kTPC", AliEmcalJet::kTPC},
  {"kTPCfid", AliEmcalJet::kTPCfid},
  {"kEMCAL", AliEmcalJet::kEMCAL},
  {"kEMCALfid", AliEmcalJet::kEMCALfid},
  {"kDCAL", AliEmcalJet::kDCAL},
  {"kDCALfid", AliEmcalJet::kDCALfid},
  {"kDCALonly", AliEmcalJet::kDCALonly},
  {"kDCALonlyfid", AliEmcalJet::kDCALonlyfid},
  {"kPHOS", AliEmcalJet::kPHOS},
  {"kPHOSfid", AliEmcalJet::kPHOSfid},
  {"kUser", AliEmcalJet::kUser}
};

const std::map<std::string, AliAnalysisTaskEmcalJetHUtils::EEfficiencyPeriodIdentifier_t> AliAnalysisTaskEmcalJetHUtils::fgkEfficiencyPeriodIdentifier = {
  { "DisableEff", AliAnalysisTaskEmcalJetHUtils::kDisableEff},
  { "LHC11h", AliAnalysisTaskEmcalJetHUtils::kLHC11h},
  { "LHC15o", AliAnalysisTaskEmcalJetHUtils::kLHC15o},
  { "LHC18q", AliAnalysisTaskEmcalJetHUtils::kLHC18qr},
  { "LHC18r", AliAnalysisTaskEmcalJetHUtils::kLHC18qr},
  { "LHC11a", AliAnalysisTaskEmcalJetHUtils::kLHC11a},
  { "pA", AliAnalysisTaskEmcalJetHUtils::kpA },
  { "pp", AliAnalysisTaskEmcalJetHUtils::kpp }
};

// LHC11h efficiency parameters for good runs
// 0-10% centrality
const double AliAnalysisTaskEmcalJetHUtils::LHC11hParam_0_10[17] = {
  0.971679, 0.0767571, 1.13355,  -0.0274484, 0.856652,  0.00536795, 3.90795e-05, 1.06889, 0.011007,
  0.447046, -0.146626, 0.919777, 0.192601,   -0.268515, 1.00243,    0.00620849,  0.709477
};
// 10-30% centrality
const double AliAnalysisTaskEmcalJetHUtils::LHC11hParam_10_30[17] = {
  0.97929,  0.0776039, 1.12213,  -0.0300645, 0.844722,  0.0134788, -0.0012333, 1.07955, 0.0116835,
  0.456608, -0.132743, 0.930964, 0.174175,   -0.267154, 0.993118,  0.00574892, 0.765256
};
// 30-50% centrality
const double AliAnalysisTaskEmcalJetHUtils::LHC11hParam_30_50[17] = {
  0.997696, 0.0816769,  1.14341,  -0.0353734, 0.752151,  0.0744259, -0.0102926, 1.01561, 0.00713274,
  0.57203,  -0.0640248, 0.947747, 0.102007,   -0.194698, 0.999164,  0.00568476, 0.7237
};
// 50-90% centrality
const double AliAnalysisTaskEmcalJetHUtils::LHC11hParam_50_90[17] = {
  0.97041, 0.0813559,  1.12151,  -0.0368797, 0.709327, 0.0701501, -0.00784043, 1.06276, 0.00676173,
  0.53607, -0.0703117, 0.982534, 0.0947881,  -0.18073, 1.03229,   0.00580109,  0.737801
};

// For pt parameters, first 5 are low pt, next 5 are high pt
// For eta parameters, first 6 are eta =< -0.04 (eta left in Eliane's def), next 6 are => -0.04 (eta right
// in Eliane's def). The last parameter normalizes the eta values such that their maximum is 1. This was apparently
// part of their definition, but was implementing by normalizing a TF1 afterwards. My implementation approach here
// is more useful when not using a TF1.
// 0-10% centrality
const double AliAnalysisTaskEmcalJetHUtils::LHC15oParam_0_10_pt[10] = { 0.8350,         0.0621,         0.0986, 0.2000,
                                    1.0124,         0.7568,         0.0277, -0.0034,
                                    0.1506 * 0.001, -0.0023 * 0.001 };
const double AliAnalysisTaskEmcalJetHUtils::LHC15oParam_0_10_eta[13] = { 1.0086,  0.0074, 0.2404, -0.1230, -0.0107,
                                     0.0427,  0.8579, 0.0088, 0.4697,  0.0772,
                                     -0.0352, 0.0645, 0.7716 };
// 10-30% centrality
const double AliAnalysisTaskEmcalJetHUtils::LHC15oParam_10_30_pt[10] = {
  0.8213, 0.0527, 0.0867, 0.1970, 1.1518, 0.7469, 0.0300, -0.0038, 0.1704 * 0.001, -0.0026 * 0.001
};
const double AliAnalysisTaskEmcalJetHUtils::LHC15oParam_10_30_eta[13] = { 0.9726,  0.0066, 0.2543, -0.1167, -0.0113,
                                     0.0400,  0.8729, 0.0122, 0.4537,  0.0965,
                                     -0.0328, 0.0623, 0.7658 };
// 30-50% centrality
const double AliAnalysisTaskEmcalJetHUtils::LHC15oParam_30_50_pt[10] = {
  0.8381, 0.0648, 0.1052, 0.1478, 1.0320, 0.7628, 0.0263, -0.0032, 0.1443 * 0.001, -0.0023 * 0.001
};
const double AliAnalysisTaskEmcalJetHUtils::LHC15oParam_30_50_eta[13] = { 0.9076,  0.0065, 0.3216, -0.1130, -0.0107,
                                     0.0456,  0.8521, 0.0073, 0.4764,  0.0668,
                                     -0.0363, 0.0668, 0.7748 };
// 50-90% centrality
const double AliAnalysisTaskEmcalJetHUtils::LHC15oParam_50_90_pt[10] = {
  0.8437, 0.0668, 0.1083, 0.2000, 0.9741, 0.7677, 0.0255, -0.0030, 0.1260 * 0.001, -0.0019 * 0.001
};
const double AliAnalysisTaskEmcalJetHUtils::LHC15oParam_50_90_eta[13] = { 1.1259,  0.0105, 0.1961, -0.1330, -0.0103,
                                     0.0440,  0.8421, 0.0066, 0.5061,  0.0580,
                                     -0.0379, 0.0651, 0.7786 };

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
    // Configure manual mode via YAML
    // Select the period
    typedef void (AliEventCuts::*MFP)();
    std::map<std::string, MFP> eventCutsPeriods = { std::make_pair("LHC11h", &AliEventCuts::SetupRun1PbPb),
                            std::make_pair("LHC15o", &AliEventCuts::SetupLHC15o),
                            std::make_pair("LHC18qr", &AliEventCuts::SetupPbPb2018) };
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
    res = yamlConfig.GetProperty({ baseName, "centralityRange" }, centRange, false);
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
 * @param[in] taskName Name of the task which is calling this function (for logging purposes).
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
 * @param[in] taskName Name of the task which is calling this function (for logging purposes).
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
     << trackCont->GetName() << ": Setting filterBits of " << filterBits << "\n";
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
     << trackCont->GetName() << ": Setting track cuts period to " << tempString << "\n";
    trackCont->SetTrackCutsPeriod(tempString.c_str());
  }
}

/**
 * Configure a cluster container according to the specified YAML configuration.
 *
 * @param[in] baseNameWithhContainer Name under which the config is stored, with the container name included.
 * @param[in] clusterCont Existing cluster container.
 * @param[in] yamlConfig YAML configuration to be used.
 * @param[in] taskName Name of the task which is calling this function (for logging purposes).
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
     << clusterCont->GetName() << ": Setting cluster energy type to " << clusterEnergyType << "\n";
    clusterCont->SetDefaultClusterEnergy(clusterEnergyType);
  }

  // NonLinCorrEnergyCut
  result = yamlConfig.GetProperty(baseNameWithContainer, "clusNonLinCorrEnergyCut", tempDouble, false);
  if (result) {
    AliDebugGeneralStream(taskName.c_str(), 2)
     << clusterCont->GetName() << ": Setting clusNonLinCorrEnergyCut of " << tempDouble << "\n";
    clusterCont->SetClusNonLinCorrEnergyCut(tempDouble);
  }

  // HadCorrEnergyCut
  result = yamlConfig.GetProperty(baseNameWithContainer, "clusHadCorrEnergyCut", tempDouble, false);
  if (result) {
    AliDebugGeneralStream(taskName.c_str(), 2)
     << clusterCont->GetName() << ": Setting clusHadCorrEnergyCut of " << tempDouble << "\n";
    clusterCont->SetClusHadCorrEnergyCut(tempDouble);
  }

  // SetIncludePHOS
  result = yamlConfig.GetProperty(baseNameWithContainer, "includePHOS", tempBool, false);
  if (result) {
    AliDebugGeneralStream(taskName.c_str(), 2) << clusterCont->GetName() << ": Setting Include PHOS to "
                          << (tempBool ? "enabled" : "disabled") << "\n";
    clusterCont->SetIncludePHOS(tempBool);
  }
}

/**
 * Determines the jet acceptance that is retrieved from a YAML configuration. Note that the result is an OR of
 * all of the individual acceptances selected in the input.
 *
 * @return The desired jet acceptance. Note that a `UInt_t` is explicitly passed for the acceptance, so it's fine to return it here.
 */
UInt_t AliAnalysisTaskEmcalJetHUtils::DetermineJetAcceptanceFromYAML(const std::vector<std::string> & selections)
{
  UInt_t jetAcceptance = 0;
  for (auto selection : selections) {
    auto sel = fgkJetAcceptanceMap.find(selection);
    AliDebugGeneralStream("AliAnalysisTaskEmcalJetHUtils", 3) << "Adding jet acceptance: " << selection << "\n";
    if (sel != fgkJetAcceptanceMap.end()) {
      jetAcceptance |= sel->second;
    } else {
      AliFatalGeneralF("AliAnalysisTaskEmcalJetHUtils", "Could not find jet acceptance with key \"%s\"",
               selection.c_str());
    }
  }
  return jetAcceptance;
}

/**
 * AddTask for Qn flow vector corrections. The AddTask `AddTaskFlowQnVectorCorrectionsNewDetConfig.C`
 * provides the right options, but it enables the calibration and event histograms, apparently without
 * the ability to turn them off. This is problematic, because the output size is quite large. This AddTask
 * will allow for the additional histograms to be turned off. The only purpose of this AddTask is to make
 * those options configurable.
 *
 * In order to make this AddTask compileable, a number of other minor changes were required.
 *
 * Note that this function uses a YAML configuration file instead of the frameworks configuration method.
 * This is done solely for convenience and to cut down on the options that don't matter when just running
 * the task on the LEGO train.
 *
 * @param[in] configFilename Filename and path of the YAML configuration file.
 * @returns AliAnalysisTaskFlowVectorCorrections object configured to provide corrected Qn vectors.
 */
AliAnalysisTaskFlowVectorCorrections* AliAnalysisTaskEmcalJetHUtils::AddTaskFlowQnVectorCorrections(
 const std::string& configFilename)
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFlowQnVectorCorrections", "No analysis manager found.");
    return 0;
  }

  // Create the correction task and manager
  AliQnCorrectionsManager* QnManager = new AliQnCorrectionsManager();
  AliAnalysisTaskFlowVectorCorrections* taskQnCorrections =
   new AliAnalysisTaskFlowVectorCorrections("FlowQnVectorCorrections");

  // Determine the task configuration
  PWG::Tools::AliYAMLConfiguration yamlConfig;
  yamlConfig.AddConfiguration(configFilename, "config");
  std::string baseName = "";

  // General configuration
  // Use VZERO centrality or multiplicity percentile for centrality determination
  // Use centrality for 2010 (and 2011?), use multiplicity of 2015 (ie. run 2)
  bool useMultiplicityPercentileForCentralityDetermination = true;
  yamlConfig.GetProperty("useMultiplicityPercentileForCentralityDetermination",
              useMultiplicityPercentileForCentralityDetermination, true);
  AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity = AliQnCorrectionsVarManagerTask::kCentVZERO;
  if (useMultiplicityPercentileForCentralityDetermination) {
    varForEventMultiplicity = AliQnCorrectionsVarManagerTask::kVZEROMultPercentile;
  }
  // Select the Z vertex, centrality for when to calibrate (and correct?).
  std::pair<double, double> zVertexRange;
  yamlConfig.GetProperty("zVertex", zVertexRange, true);
  std::pair<double, double> centRange;
  yamlConfig.GetProperty("centrality", centRange, true);
  // Select only events validated for centrality calibration
  // Check information about your runs of interest in
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliMultSelectionCalibStatus.
  // Learn more about its usage in https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentralityCodeSnippets
  bool useOnlyCentralityCalibratedEvents = false;
  yamlConfig.GetProperty("useOnlyCentralityCalibratedEvents", useOnlyCentralityCalibratedEvents, true);
  // Select runs to use when gather calibration data.
  std::vector<std::string> listOfRuns;
  yamlConfig.GetProperty("runsToUseDuringCalibration", listOfRuns, false);
  // Physics selection. Defaults to kAnyINT
  std::vector<std::string> physicsSelection;
  bool res = yamlConfig.GetProperty("physicsSelection", physicsSelection, false);
  if (res) {
    taskQnCorrections->SelectCollisionCandidates(
     AliEmcalContainerUtils::DeterminePhysicsSelectionFromYAML(physicsSelection));
  } else {
    // Defaults to using kAnyINT
    taskQnCorrections->SelectCollisionCandidates(AliVEvent::kAnyINT);
  }

  // Location of correction histograms
  baseName = "correctionHistograms";
  std::string correctionsSource = "";
  yamlConfig.GetProperty({ baseName, "source" }, correctionsSource, true);
  std::string correctionsFilePath = "";
  yamlConfig.GetProperty({ baseName, "path" }, correctionsFilePath, true);
  std::string correctionsFileName = "";
  yamlConfig.GetProperty({ baseName, "filename" }, correctionsFileName, true);

  // Detector configuration (optional - all are off by default)
  baseName = "detectors";
  bool useTPC = false;
  yamlConfig.GetProperty({ baseName, "TPC" }, useTPC, false);
  bool useSPD = false;
  yamlConfig.GetProperty({ baseName, "SPD" }, useSPD, false);
  bool useVZERO = false;
  yamlConfig.GetProperty({ baseName, "VZERO" }, useVZERO, false);
  bool useTZERO = false;
  yamlConfig.GetProperty({ baseName, "TZERO" }, useTZERO, false);
  bool useFMD = false;
  yamlConfig.GetProperty({ baseName, "FMD" }, useFMD, false);
  bool useRawFMD = false;
  yamlConfig.GetProperty({ baseName, "RawFMD" }, useRawFMD, false);
  bool useZDC = false;
  yamlConfig.GetProperty({ baseName, "ZDC" }, useZDC, false);

  // Outputs configuration
  baseName = "outputs";
  bool fillQVectorTree = false;
  yamlConfig.GetProperty({ baseName, "QVectorTree" }, fillQVectorTree, false);
  bool fillQAHistograms = true;
  yamlConfig.GetProperty({ baseName, "QAHistograms" }, fillQAHistograms, false);
  bool fillNveQAHistograms = true;
  yamlConfig.GetProperty({ baseName, "NveQAHistograms" }, fillNveQAHistograms, false);
  bool fillOutputHistograms = true;
  yamlConfig.GetProperty({ baseName, "OutputHistograms" }, fillOutputHistograms, false);
  bool fillExchangeContainerWithQvectors = true;
  yamlConfig.GetProperty({ baseName, "ExchangeContainerWithQvectors" }, fillExchangeContainerWithQvectors, false);
  bool fillEventQA = true;
  yamlConfig.GetProperty({ baseName, "EventQA" }, fillEventQA, false);

  // Create a string to describe the configuration.
  // It allows the user to verify that the configuration was accessed successfully and that the
  // values were extracted into the proper variables.
  std::stringstream tempSS;
  tempSS << std::boolalpha;
  tempSS << "Flow Qn vector corrections configuration:\n";
  tempSS << "Use multiplicity percentile for centrality determination: "
      << useMultiplicityPercentileForCentralityDetermination << "\n";
  tempSS << "Z vertex: [" << zVertexRange.first << ", " << zVertexRange.second << "]\n";
  tempSS << "Centrality: [" << centRange.first << ", " << centRange.second << "]\n";
  tempSS << "Use only centrality calibrated events: " << useOnlyCentralityCalibratedEvents << "\n";
  tempSS << "Runs to use during calibration:\n";
  bool atLeastOneRun = false;
  for (const auto& run : listOfRuns) {
    atLeastOneRun = true;
    tempSS << "\t" << run << "\n";
  }
  if (!atLeastOneRun) {
    tempSS << "\tNone\n";
  }
  tempSS << "Physics selection: " << taskQnCorrections->GetCollisionCandidates() << "\n";

  tempSS << "Correction histograms:\n";
  tempSS << "\tSource: " << correctionsSource << "\n";
  tempSS << "\tPath: " << correctionsFilePath << "\n";
  tempSS << "\tName: " << correctionsFileName << "\n";

  tempSS << "Detectors:\n";
  tempSS << "\tTPC: " << useTPC << "\n";
  tempSS << "\tSPD: " << useSPD << "\n";
  tempSS << "\tVZERO: " << useVZERO << "\n";
  tempSS << "\tTZERO: " << useTZERO << "\n";
  tempSS << "\tFMD: " << useFMD << "\n";
  tempSS << "\tRaw FMD: " << useRawFMD << "\n";
  tempSS << "\tZDC: " << useZDC << "\n";

  tempSS << "Outputs:\n";
  tempSS << "\tQ vector tree: " << fillQVectorTree << "\n";
  tempSS << "\tQA histograms: " << fillQAHistograms << "\n";
  tempSS << "\tNve QA histograms: " << fillNveQAHistograms << "\n";
  tempSS << "\tOutput histograms: " << fillOutputHistograms << "\n";
  tempSS << "\tExchange containers with Q vectors: " << fillExchangeContainerWithQvectors << "\n";
  tempSS << "\tEvent QA: " << fillEventQA << "\n";

  // Print the task configuration
  // Print outside of the ALICE Log system to ensure that it is always available!
  std::cout << tempSS.str();

  // Convert list of runs from vector to TObjArray so they can be passed to the task.
  // Allocate using new so we can pass it into the class. Otherwise, it will segfault because the task attempts
  // to access the object after it's gone out of scope (instead of copying the TObjArray).
  TObjArray* listOfRunsTOBJ = new TObjArray();
  listOfRunsTOBJ->SetOwner(kTRUE);
  for (const auto& run : listOfRuns) {
    listOfRunsTOBJ->Add(new TObjString(run.c_str()));
  }

  /* let's establish the event cuts for event selection */
  AliQnCorrectionsCutsSet* eventCuts = new AliQnCorrectionsCutsSet();
  eventCuts->Add(
   new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kVtxZ, zVertexRange.first, zVertexRange.second));
  eventCuts->Add(new AliQnCorrectionsCutWithin(varForEventMultiplicity, centRange.first, centRange.second));
  taskQnCorrections->SetEventCuts(eventCuts);
  taskQnCorrections->SetUseOnlyCentCalibEvents(useOnlyCentralityCalibratedEvents);

  TString histClass = "";
  histClass += "Event_NoCuts;";
  histClass += "Event_Analysis;";
  histClass += "TrackQA_NoCuts;";

  /* add the selected detectors */
  if (useTPC) {
    FlowVectorCorrections::AddTPC(taskQnCorrections, QnManager, varForEventMultiplicity);
    histClass += "TrackQA_TPC;";
  }
  if (useSPD) {
    FlowVectorCorrections::AddSPD(taskQnCorrections, QnManager, varForEventMultiplicity);
    histClass += "TrackletQA_SPD;";
  }
  if (useVZERO) {
    FlowVectorCorrections::AddVZERO(taskQnCorrections, QnManager, varForEventMultiplicity);
  }
  if (useTZERO) {
    FlowVectorCorrections::AddTZERO(taskQnCorrections, QnManager, varForEventMultiplicity);
  }
  if (useFMD) {
    FlowVectorCorrections::AddFMD(taskQnCorrections, QnManager, varForEventMultiplicity);
  }
  if (useRawFMD) {
    FlowVectorCorrections::AddRawFMD(taskQnCorrections, QnManager, varForEventMultiplicity);
  }
  if (useZDC) {
    FlowVectorCorrections::AddZDC(taskQnCorrections, QnManager, varForEventMultiplicity);
  }

  QnManager->SetShouldFillQnVectorTree(fillQVectorTree);
  QnManager->SetShouldFillQAHistograms(fillQAHistograms);
  QnManager->SetShouldFillNveQAHistograms(fillNveQAHistograms);
  QnManager->SetShouldFillOutputHistograms(fillOutputHistograms);

  taskQnCorrections->SetFillExchangeContainerWithQvectors(fillExchangeContainerWithQvectors);
  taskQnCorrections->SetFillEventQA(fillEventQA);

  taskQnCorrections->SetAliQnCorrectionsManager(QnManager);
  taskQnCorrections->DefineInOutput();
  taskQnCorrections->SetRunsLabels(listOfRunsTOBJ);

  /* let's handle the calibration file */
  AliInfoGeneralStream("PWGJE::EMCALJetTasks::AddTaskFlowQnVectorCorrections")
   << "=================== CALIBRATION FILE =============================================\n";
  std::string inputCalibrationFilename = correctionsFilePath + "/" + correctionsFileName;
  if (correctionsSource == "local") {
    AliInfoGeneralStream("PWGJE::EMCALJetTasks::AddTaskFlowQnVectorCorrections")
     << "\t File " << inputCalibrationFilename << "\n\t being taken locally when building the task object.\n";
    taskQnCorrections->SetCalibrationHistogramsFile(AliAnalysisTaskFlowVectorCorrections::CALIBSRC_local,
                            inputCalibrationFilename.c_str());
  } else if (correctionsSource == "aliensingle") {
    AliInfoGeneralStream("PWGJE::EMCALJetTasks::AddTaskFlowQnVectorCorrections")
     << "\t File " << inputCalibrationFilename << " being taken from alien in the execution nodes\n";
    taskQnCorrections->SetCalibrationHistogramsFile(AliAnalysisTaskFlowVectorCorrections::CALIBSRC_aliensingle,
                            inputCalibrationFilename.c_str());
  } else if (correctionsSource == "alienmultiple") {
    AliInfoGeneralStream("PWGJE::EMCALJetTasks::AddTaskFlowQnVectorCorrections")
     << "\t File " << inputCalibrationFilename
     << " being taken from alien in the execution nodes on a per run basis\n";
    taskQnCorrections->SetCalibrationHistogramsFile(AliAnalysisTaskFlowVectorCorrections::CALIBSRC_alienmultiple,
                            inputCalibrationFilename.c_str());
  } else if (correctionsSource == "OADBsingle") {
    AliInfoGeneralStream("PWGJE::EMCALJetTasks::AddTaskFlowQnVectorCorrections")
     << "\t File " << inputCalibrationFilename << " being taken from OADB in the execution nodes\n";
    taskQnCorrections->SetCalibrationHistogramsFile(AliAnalysisTaskFlowVectorCorrections::CALIBSRC_OADBsingle,
                            inputCalibrationFilename.c_str());
  } else if (correctionsSource == "OADBmultiple") {
    AliInfoGeneralStream("PWGJE::EMCALJetTasks::AddTaskFlowQnVectorCorrections")
     << "\t File " << inputCalibrationFilename
     << " being taken from OADB in the execution nodes on a per run basis\n";
    taskQnCorrections->SetCalibrationHistogramsFile(AliAnalysisTaskFlowVectorCorrections::CALIBSRC_OADBmultiple,
                            inputCalibrationFilename.c_str());
  } else {
    AliErrorGeneralStream("AddTaskFlowQnVectorCorrections")
     << "\t CALIBRATION FILE SOURCE NOT SUPPORTED. ABORTING!!!";
    return NULL;
  }
  AliInfoGeneralStream("PWGJE::EMCALJetTasks::AddTaskFlowQnVectorCorrections")
   << "==================================================================================\n";

  AliQnCorrectionsHistos* hists = taskQnCorrections->GetEventHistograms();
  FlowVectorCorrections::DefineHistograms(QnManager, hists, histClass);

  mgr->AddTask(taskQnCorrections);
  mgr->ConnectInput(taskQnCorrections, 0, mgr->GetCommonInputContainer());

  // create output containers
  if (QnManager->GetShouldFillOutputHistograms()) {
    AliAnalysisDataContainer* cOutputHist =
     mgr->CreateContainer(QnManager->GetCalibrationHistogramsContainerName(), TList::Class(),
                AliAnalysisManager::kOutputContainer, "CalibrationHistograms.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotHistQn(), cOutputHist);
  }

  if (QnManager->GetShouldFillQnVectorTree()) {
    AliAnalysisDataContainer* cOutputQvec = mgr->CreateContainer(
     "CalibratedQvector", TTree::Class(), AliAnalysisManager::kOutputContainer, "QvectorsTree.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotTree(), cOutputQvec);
  }

  if (QnManager->GetShouldFillQAHistograms()) {
    AliAnalysisDataContainer* cOutputHistQA =
     mgr->CreateContainer(QnManager->GetCalibrationQAHistogramsContainerName(), TList::Class(),
                AliAnalysisManager::kOutputContainer, "CalibrationQA.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotHistQA(), cOutputHistQA);
  }

  if (QnManager->GetShouldFillNveQAHistograms()) {
    AliAnalysisDataContainer* cOutputHistNveQA =
     mgr->CreateContainer(QnManager->GetCalibrationNveQAHistogramsContainerName(), TList::Class(),
                AliAnalysisManager::kOutputContainer, "CalibrationQA.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotHistNveQA(), cOutputHistNveQA);
  }

  if (taskQnCorrections->GetFillEventQA()) {
    AliAnalysisDataContainer* cOutputQnEventQA =
     mgr->CreateContainer("QnEventQA", TList::Class(), AliAnalysisManager::kOutputContainer, "QnEventQA.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotEventQA(), cOutputQnEventQA);
  }

  AliAnalysisDataContainer* cOutputQvecList = mgr->CreateContainer(
   "CalibratedQvectorList", TList::Class(), AliAnalysisManager::kExchangeContainer, "QvectorsList.root");

  if (taskQnCorrections->GetFillExchangeContainerWithQvectors())
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotGetListQnVectors(), cOutputQvecList);

  return taskQnCorrections;
}

/**
 * Get the background subtracted jet pt.
 *
 * @param[in] jet Jet to be subtracted.
 * @param[in] rho Rho value for the jet collection.
 * @returns The rho corrected jet pt (or just the raw jet pt if no rho value is provided).
 */
double AliAnalysisTaskEmcalJetHUtils::GetJetPt(const AliEmcalJet* jet, const double rho)
{
  double pT = jet->Pt() - rho * jet->Area();
  return pT;
}

/**
 * Calculate the tracking efficiency for LHC11h - run 1 PbPb at 2.76 TeV. See Joel's jet-hadron analysis note.
 *
 * @param[in] trackPt Track pt
 * @param[in] trackEta Track eta
 * @param[in] centralityBin Centrality bin of the current event.
 * @param[in] taskName Name of the task which is calling this function (for logging purposes).
 * @returns The efficiency of measuring the given single track.
 */
double AliAnalysisTaskEmcalJetHUtils::LHC11hTrackingEfficiency(const double trackPt, const double trackEta, const int centralityBin, const std::string & taskName)
{
  // Setup
  double etaAxis = 0;
  double ptAxis = 0;
  double efficiency = 1;

  // Assumes that the centrality bins follow (as defined in AliAnalysisTaskEmcal)
  // 0 = 0-10%
  // 1 = 10-30%
  // 2 = 30-50%
  // 3 = 50-90%
  switch (centralityBin) {
    case 0 :
      // Parameter values for GOOD TPC (LHC11h) runs (0-10%):
      ptAxis =
       (trackPt < 2.9) * (LHC11hParam_0_10[0] * exp(-pow(LHC11hParam_0_10[1] / trackPt, LHC11hParam_0_10[2])) +
                 LHC11hParam_0_10[3] * trackPt) +
       (trackPt >= 2.9) *
        (LHC11hParam_0_10[4] + LHC11hParam_0_10[5] * trackPt + LHC11hParam_0_10[6] * trackPt * trackPt);
      etaAxis =
       (trackEta < 0.0) *
        (LHC11hParam_0_10[7] * exp(-pow(LHC11hParam_0_10[8] / std::abs(trackEta + 0.91), LHC11hParam_0_10[9])) +
         LHC11hParam_0_10[10] * trackEta) +
       (trackEta >= 0.0 && trackEta <= 0.4) *
        (LHC11hParam_0_10[11] + LHC11hParam_0_10[12] * trackEta + LHC11hParam_0_10[13] * trackEta * trackEta) +
       (trackEta > 0.4) * (LHC11hParam_0_10[14] *
                 exp(-pow(LHC11hParam_0_10[15] / std::abs(-trackEta + 0.91), LHC11hParam_0_10[16])));
      efficiency = ptAxis * etaAxis;
      break;

    case 1:
      // Parameter values for GOOD TPC (LHC11h) runs (10-30%):
      ptAxis = (trackPt < 2.9) *
            (LHC11hParam_10_30[0] * exp(-pow(LHC11hParam_10_30[1] / trackPt, LHC11hParam_10_30[2])) +
            LHC11hParam_10_30[3] * trackPt) +
           (trackPt >= 2.9) * (LHC11hParam_10_30[4] + LHC11hParam_10_30[5] * trackPt +
                     LHC11hParam_10_30[6] * trackPt * trackPt);
      etaAxis =
       (trackEta < 0.0) * (LHC11hParam_10_30[7] *
                  exp(-pow(LHC11hParam_10_30[8] / std::abs(trackEta + 0.91), LHC11hParam_10_30[9])) +
                 LHC11hParam_10_30[10] * trackEta) +
       (trackEta >= 0.0 && trackEta <= 0.4) * (LHC11hParam_10_30[11] + LHC11hParam_10_30[12] * trackEta +
                           LHC11hParam_10_30[13] * trackEta * trackEta) +
       (trackEta > 0.4) * (LHC11hParam_10_30[14] *
                 exp(-pow(LHC11hParam_10_30[15] / std::abs(-trackEta + 0.91), LHC11hParam_10_30[16])));
      efficiency = ptAxis * etaAxis;
      break;

    case 2:
      // Parameter values for GOOD TPC (LHC11h) runs (30-50%):
      ptAxis = (trackPt < 2.9) *
            (LHC11hParam_30_50[0] * exp(-pow(LHC11hParam_30_50[1] / trackPt, LHC11hParam_30_50[2])) +
            LHC11hParam_30_50[3] * trackPt) +
           (trackPt >= 2.9) * (LHC11hParam_30_50[4] + LHC11hParam_30_50[5] * trackPt +
                     LHC11hParam_30_50[6] * trackPt * trackPt);
      etaAxis =
       (trackEta < 0.0) * (LHC11hParam_30_50[7] *
                  exp(-pow(LHC11hParam_30_50[8] / std::abs(trackEta + 0.91), LHC11hParam_30_50[9])) +
                 LHC11hParam_30_50[10] * trackEta) +
       (trackEta >= 0.0 && trackEta <= 0.4) * (LHC11hParam_30_50[11] + LHC11hParam_30_50[12] * trackEta +
                           LHC11hParam_30_50[13] * trackEta * trackEta) +
       (trackEta > 0.4) * (LHC11hParam_30_50[14] *
                 exp(-pow(LHC11hParam_30_50[15] / std::abs(-trackEta + 0.91), LHC11hParam_30_50[16])));
      efficiency = ptAxis * etaAxis;
      break;

    case 3:
      // Parameter values for GOOD TPC (LHC11h) runs (50-90%):
      ptAxis = (trackPt < 2.9) *
            (LHC11hParam_50_90[0] * exp(-pow(LHC11hParam_50_90[1] / trackPt, LHC11hParam_50_90[2])) +
            LHC11hParam_50_90[3] * trackPt) +
           (trackPt >= 2.9) * (LHC11hParam_50_90[4] + LHC11hParam_50_90[5] * trackPt +
                     LHC11hParam_50_90[6] * trackPt * trackPt);
      etaAxis =
       (trackEta < 0.0) * (LHC11hParam_50_90[7] *
                  exp(-pow(LHC11hParam_50_90[8] / std::abs(trackEta + 0.91), LHC11hParam_50_90[9])) +
                 LHC11hParam_50_90[10] * trackEta) +
       (trackEta >= 0.0 && trackEta <= 0.4) * (LHC11hParam_50_90[11] + LHC11hParam_50_90[12] * trackEta +
                           LHC11hParam_50_90[13] * trackEta * trackEta) +
       (trackEta > 0.4) * (LHC11hParam_50_90[14] *
                 exp(-pow(LHC11hParam_50_90[15] / std::abs(-trackEta + 0.91), LHC11hParam_50_90[16])));
      efficiency = ptAxis * etaAxis;
      break;

    default:
      AliErrorGeneralStream(taskName.c_str()) << "Invalid centrality for determine tracking efficiency.\n";
      efficiency = 0;
  }

  return efficiency;
}

/**
 * Determine the pt efficiency axis for LHC15o. This is the main interface
 * for getting the efficiency.
 *
 * @param[in] trackEta Track eta.
 * @param[in] params Parameters for use with the function.
 * @returns The efficiency associated with the eta parameterization.
 */
double AliAnalysisTaskEmcalJetHUtils::LHC15oPtEfficiency(const double trackPt, const double params[10])
{
  return ((trackPt <= 3.5) * LHC15oLowPtEfficiencyImpl(trackPt, params, 0) +
      (trackPt > 3.5) * LHC15oHighPtEfficiencyImpl(trackPt, params, 5));
}

/**
 * Determine the pt efficiency axis for low pt tracks in LHC15o. Implementation function.
 *
 * @param[in] trackEta Track eta.
 * @param[in] params Parameters for use with the function.
 * @param[in] index Index where it should begin accessing the parameters.
 * @returns The efficiency associated with the eta parameterization.
 */
double AliAnalysisTaskEmcalJetHUtils::LHC15oLowPtEfficiencyImpl(const double trackPt, const double params[10], const int index)
{
  return (params[index + 0] + -1.0 * params[index + 1] / trackPt) +
      params[index + 2] * TMath::Gaus(trackPt, params[index + 3], params[index + 4]);
}

/**
 * Determine the pt efficiency axis for high pt tracks in LHC15o. Implementation function.
 *
 * @param[in] trackEta Track eta.
 * @param[in] params Parameters for use with the function.
 * @param[in] index Index where it should begin accessing the parameters.
 * @returns The efficiency associated with the eta parameterization.
 */
double AliAnalysisTaskEmcalJetHUtils::LHC15oHighPtEfficiencyImpl(const double trackPt, const double params[10], const int index)
{
  return params[index + 0] + params[index + 1] * trackPt + params[index + 2] * std::pow(trackPt, 2) +
      params[index + 3] * std::pow(trackPt, 3) + params[index + 4] * std::pow(trackPt, 4);
}

/**
 * Determine the eta efficiency axis for LHC15o.
 *
 * @param[in] trackEta Track eta.
 * @param[in] params Parameters for use with the function.
 * @returns The efficiency associated with the eta parameterization.
 */
double AliAnalysisTaskEmcalJetHUtils::LHC15oEtaEfficiency(const double trackEta, const double params[13])
{
  // Just modify the arguments - the function is the same.
  return ((trackEta <= -0.04) * LHC15oEtaEfficiencyImpl(trackEta, params, 0) +
      (trackEta > -0.04) * LHC15oEtaEfficiencyImpl(trackEta, params, 6));
}

/**
 * Determine the eta efficiency axis for LHC15o. Implementation function.
 *
 * @param[in] trackEta Track eta.
 * @param[in] params Parameters for use with the function.
 * @param[in] index Index where it should begin accessing the parameters.
 * @returns The efficiency associated with the eta parameterization.
 */
double AliAnalysisTaskEmcalJetHUtils::LHC15oEtaEfficiencyImpl(const double trackEta, const double params[13],
                               const int index)
{
  // We need to multiply the track eta by -1 if we are looking at eta > 0 (which corresponds to
  // the second set of parameters, such that the index is greater than 0).
  int sign = index > 0 ? -1 : 1;
  return (params[index + 0] *
       std::exp(-1.0 * std::pow(params[index + 1] / std::abs(sign * trackEta + 0.91), params[index + 2])) +
      params[index + 3] * trackEta + params[index + 4] * TMath::Gaus(trackEta, -0.04, params[index + 5])) /
      params[12];
}

/**
 * Calculate the track efficiency for LHC15o - PbPb at 5.02 TeV. See the gamma-hadron analysis (from Eliane via Michael).
 *
 * @param[in] trackPt Track pt
 * @param[in] trackEta Track eta
 * @param[in] centralityBin Centrality bin of the current event.
 * @param[in] taskName Name of the task which is calling this function (for logging purposes).
 * @returns The efficiency of measuring the given single track.
 */
double AliAnalysisTaskEmcalJetHUtils::LHC15oTrackingEfficiency(const double trackPt, const double trackEta, const int centralityBin, const std::string & taskName)
{
  // We use the switch to determine the parameters needed to call the functions.
  // Assumes that the centrality bins follow (as defined in AliAnalysisTaskEmcal)
  // 0 = 0-10%
  // 1 = 10-30%
  // 2 = 30-50%
  // 3 = 50-90%
  const double* ptParams = nullptr;
  const double* etaParams = nullptr;
  switch (centralityBin) {
    case 0:
      ptParams = LHC15oParam_0_10_pt;
      etaParams = LHC15oParam_0_10_eta;
      break;
    case 1:
      ptParams = LHC15oParam_10_30_pt;
      etaParams = LHC15oParam_10_30_eta;
      break;
    case 2:
      ptParams = LHC15oParam_30_50_pt;
      etaParams = LHC15oParam_30_50_eta;
      break;
    case 3:
      ptParams = LHC15oParam_50_90_pt;
      etaParams = LHC15oParam_50_90_eta;
      break;
    default:
      AliFatalGeneral(taskName.c_str(), "Invalid centrality for determine tracking efficiency.\n");
  }

  // Calculate the efficiency using the parameters.
  double ptAxis = LHC15oPtEfficiency(trackPt, ptParams);
  double etaAxis = LHC15oEtaEfficiency(trackEta, etaParams);
  double efficiency = ptAxis * etaAxis;

  return efficiency;
}

/**
 * Calculate the track efficiency for LHC11a - pp at 2.76 TeV. Calculated using LHC12f1a. See the jet-hadron
 * analysis note for more details.
 *
 * @param[in] trackPt Track pt
 * @param[in] trackEta Track eta
 * @param[in] centralityBin Centrality bin of the current event.
 * @param[in] taskName Name of the task which is calling this function (for logging purposes).
 * @returns The efficiency of measuring the given single track.
 */
double AliAnalysisTaskEmcalJetHUtils::LHC11aTrackingEfficiency(const double trackPt, const double trackEta, const int centralityBin, const std::string & taskName)
{
  // Pt axis
  // If the trackPt > 6 GeV, then all we need is this coefficient
  Double_t coefficient = 0.898052; // p6
  if (trackPt < 6) {
    coefficient = (1 + -0.442232 * trackPt              // p0
            + 0.501831 * std::pow(trackPt, 2)   // p1
            + -0.252024 * std::pow(trackPt, 3)  // p2
            + 0.062964 * std::pow(trackPt, 4)   // p3
            + -0.007681 * std::pow(trackPt, 5)  // p4
            + 0.000365 * std::pow(trackPt, 6)); // p5
  }

  // Eta axis
  double efficiency = coefficient * (1 + 0.402825 * std::abs(trackEta)            // p7
                    + -2.213152 * std::pow(trackEta, 2)          // p8
                    + 4.311098 * std::abs(std::pow(trackEta, 3)) // p9
                    + -2.778200 * std::pow(trackEta, 4));        // p10

  return efficiency;
}

/**
 * Note that this function relies on the user taking advantage of the centrality binning in ``AliAnalysisTaskEmcal``. If it's
 * not used, the proper centrality dependent efficiency may not be applied.
 *
 * @param[in] trackPt Track pt
 * @param[in] trackEta Track eta
 * @param[in] centralityBin Centrality bin of the current event.
 * @param[in] efficiencyPeriodIdentifier Identifies the efficiency which should be applied.
 * @param[in] taskName Name of the task which is calling this function (for logging purposes).
 * @returns The efficiency of measuring the given single track.
 */
double AliAnalysisTaskEmcalJetHUtils::DetermineTrackingEfficiency(
 const double trackPt, const double trackEta, const int centralityBin,
 const EEfficiencyPeriodIdentifier_t efficiencyPeriodIdentifier, const std::string& taskName)
{
  // Efficiency is determined entirely based on the given efficiency period.
  double efficiency = 1;
  switch (efficiencyPeriodIdentifier) {
    case AliAnalysisTaskEmcalJetHUtils::kDisableEff:
      efficiency = 1;
      break;
    case AliAnalysisTaskEmcalJetHUtils::kLHC11h:
      efficiency = LHC11hTrackingEfficiency(trackPt, trackEta, centralityBin, taskName);
      break;
    case AliAnalysisTaskEmcalJetHUtils::kLHC15o:
      efficiency = LHC15oTrackingEfficiency(trackPt, trackEta, centralityBin, taskName);
      break;
    case AliAnalysisTaskEmcalJetHUtils::kLHC11a:
      efficiency = LHC11aTrackingEfficiency(trackPt, trackEta, centralityBin, taskName);
      break;
    case AliAnalysisTaskEmcalJetHUtils::kLHC18qr:
    case AliAnalysisTaskEmcalJetHUtils::kpA:
    case AliAnalysisTaskEmcalJetHUtils::kpp:
      // Intetionally fall through for LHC18, kpA
      AliFatalGeneral(taskName.c_str(),
              TString::Format("Tracking efficiency for period identifier %d is not yet implemented.",
                      efficiencyPeriodIdentifier));
      break;
    default:
      // No efficiency period option selected. Notify the user.
      // ie. The efficiency correction is disabled.
      AliErrorGeneralStream(taskName.c_str())
       << "No single track efficiency setting selected! Please select one.\n";
      efficiency = 0.;
  }

  return efficiency;
}

/**
 * Add the VZERO configuration. Copied directly from ``AddTaskFlowQnVectorCorrectionsNewDetConfig.h``.
 */
void FlowVectorCorrections::AddVZERO(AliAnalysisTaskFlowVectorCorrections* task, AliQnCorrectionsManager* QnManager,
                   AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity)
{
  Bool_t VZEROchannels[4][64];
  for (Int_t iv0 = 0; iv0 < 4; iv0++)
    for (Int_t ich = 0; ich < 64; ich++)
      VZEROchannels[iv0][ich] = kFALSE;

  for (Int_t ich = 32; ich < 64; ich++)
    VZEROchannels[0][ich] = kTRUE; // channel list: kTRUE if channel should be used
  for (Int_t ich = 0; ich < 32; ich++)
    VZEROchannels[1][ich] = kTRUE;
  for (Int_t ich = 0; ich < 64; ich++)
    VZEROchannels[2][ich] = kTRUE;

  Int_t channelGroups[64];
  for (Int_t ich = 0; ich < 64; ich++)
    channelGroups[ich] = Int_t(ich / 8);

  //-----------------------------------------------------------
  // Our event classes for V0
  //
  const Int_t nVZEROdim = 2;
  AliQnCorrectionsEventClassVariablesSet* CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nVZEROdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4 }, { -7.0, 1 }, { 7.0, 8 }, { 10.0, 1 } };
  Double_t Ctbinning[][2] = { { 0.0, 2 }, { 100.0, 100 } };
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   AliQnCorrectionsVarManagerTask::kVtxZ, task->VarName(AliQnCorrectionsVarManagerTask::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   varForEventMultiplicity, Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the VZERO detector */
  AliQnCorrectionsDetector* VZERO = new AliQnCorrectionsDetector("VZERO", AliQnCorrectionsVarManagerTask::kVZERO);

  /* the VZEROA detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* VZEROAconf =
   new AliQnCorrectionsDetectorConfigurationChannels("VZEROAQoverM", CorrEventClasses, 64, /* number of channels */
                            4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROAconf->SetChannelsScheme(VZEROchannels[0], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROAconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization* eqA = new AliQnCorrectionsInputGainEqualization();
  eqA->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqA->SetShift(1.0);
  eqA->SetScale(0.1);
  eqA->SetUseChannelGroupsWeights(kTRUE);
  VZEROAconf->AddCorrectionOnInputData(eqA);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROAconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* alignA = new AliQnCorrectionsQnVectorAlignment();
  alignA->SetHarmonicNumberForAlignment(2);
  alignA->SetReferenceConfigurationForAlignment("TPCQoverM");
  VZEROAconf->AddCorrectionOnQnVector(alignA);
  /* lets configrure the QA histograms */
  VZEROAconf->SetQACentralityVar(varForEventMultiplicity);
  VZEROAconf->SetQAMultiplicityAxis(100, 0.0, 500.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleA = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleA->SetApplyTwist(kTRUE);
  twScaleA->SetApplyRescale(kTRUE);
  twScaleA->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleA->SetReferenceConfigurationsForTwistAndRescale("TPCQoverM", "VZEROCQoverM");
  /* now we add it to the detector configuration */
  VZEROAconf->AddCorrectionOnQnVector(twScaleA);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROAconf);

  /* the VZEROC detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* VZEROCconf =
   new AliQnCorrectionsDetectorConfigurationChannels("VZEROCQoverM", CorrEventClasses, 64, /* number of channels */
                            4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROCconf->SetChannelsScheme(VZEROchannels[1], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROCconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization* eqC = new AliQnCorrectionsInputGainEqualization();
  eqC->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqC->SetShift(1.0);
  eqC->SetScale(0.1);
  eqC->SetUseChannelGroupsWeights(kTRUE);
  VZEROCconf->AddCorrectionOnInputData(eqC);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROCconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* alignC = new AliQnCorrectionsQnVectorAlignment();
  alignC->SetHarmonicNumberForAlignment(2);
  alignC->SetReferenceConfigurationForAlignment("TPCQoverM");
  VZEROCconf->AddCorrectionOnQnVector(alignC);
  /* lets configrure the QA histograms */
  VZEROCconf->SetQACentralityVar(varForEventMultiplicity);
  VZEROCconf->SetQAMultiplicityAxis(100, 0.0, 500.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleC = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleC->SetApplyTwist(kTRUE);
  twScaleC->SetApplyRescale(kTRUE);
  twScaleC->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleC->SetReferenceConfigurationsForTwistAndRescale("TPCQoverM", "VZEROAQoverM");
  /* now we add it to the detector configuration */
  VZEROCconf->AddCorrectionOnQnVector(twScaleC);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROCconf);

  /* the full VZERO detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* VZEROconf =
   new AliQnCorrectionsDetectorConfigurationChannels("VZEROQoverM", CorrEventClasses, 64, /* number of channels */
                            4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROconf->SetChannelsScheme(VZEROchannels[2], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization* eq = new AliQnCorrectionsInputGainEqualization();
  eq->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eq->SetShift(1.0);
  eq->SetScale(0.1);
  eq->SetUseChannelGroupsWeights(kTRUE);
  VZEROconf->AddCorrectionOnInputData(eq);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* align = new AliQnCorrectionsQnVectorAlignment();
  align->SetHarmonicNumberForAlignment(2);
  align->SetReferenceConfigurationForAlignment("TPCQoverM");
  VZEROconf->AddCorrectionOnQnVector(align);
  /* lets configrure the QA histograms */
  VZEROconf->SetQACentralityVar(varForEventMultiplicity);
  VZEROconf->SetQAMultiplicityAxis(100, 0.0, 500.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScale = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScale->SetApplyTwist(kTRUE);
  twScale->SetApplyRescale(kTRUE);
  twScale->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScale->SetReferenceConfigurationsForTwistAndRescale("TPCQoverM", "VZEROCQoverM");
  /* now we add it to the detector configuration */
  VZEROconf->AddCorrectionOnQnVector(twScale);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROconf);

  /* the VZEROA detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* VZEROAconfQoverQlength =
   new AliQnCorrectionsDetectorConfigurationChannels("VZEROAQoverQlength", CorrEventClasses,
                            64, /* number of channels */
                            4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROAconfQoverQlength->SetChannelsScheme(VZEROchannels[0], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROAconfQoverQlength->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverQlength);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization* eqAQoverQlength = new AliQnCorrectionsInputGainEqualization();
  eqAQoverQlength->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqAQoverQlength->SetShift(1.0);
  eqAQoverQlength->SetScale(0.1);
  eqAQoverQlength->SetUseChannelGroupsWeights(kTRUE);
  VZEROAconfQoverQlength->AddCorrectionOnInputData(eqAQoverQlength);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROAconfQoverQlength->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* alignAQoverQlength = new AliQnCorrectionsQnVectorAlignment();
  alignAQoverQlength->SetHarmonicNumberForAlignment(2);
  alignAQoverQlength->SetReferenceConfigurationForAlignment("TPCQoverQlength");
  VZEROAconfQoverQlength->AddCorrectionOnQnVector(alignAQoverQlength);
  /* lets configrure the QA histograms */
  VZEROAconfQoverQlength->SetQACentralityVar(varForEventMultiplicity);
  VZEROAconfQoverQlength->SetQAMultiplicityAxis(100, 0.0, 500.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleAQoverQlength = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleAQoverQlength->SetApplyTwist(kTRUE);
  twScaleAQoverQlength->SetApplyRescale(kTRUE);
  twScaleAQoverQlength->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleAQoverQlength->SetReferenceConfigurationsForTwistAndRescale("TPCQoverQlength", "VZEROCQoverQlength");
  /* now we add it to the detector configuration */
  VZEROAconfQoverQlength->AddCorrectionOnQnVector(twScaleAQoverQlength);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROAconfQoverQlength);

  /* the VZEROC detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* VZEROCconfQoverQlength =
   new AliQnCorrectionsDetectorConfigurationChannels("VZEROCQoverQlength", CorrEventClasses,
                            64, /* number of channels */
                            4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROCconfQoverQlength->SetChannelsScheme(VZEROchannels[1], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROCconfQoverQlength->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverQlength);
  /* lets configure the equalization of input data */
  eqAQoverQlength = new AliQnCorrectionsInputGainEqualization();
  eqAQoverQlength->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqAQoverQlength->SetShift(1.0);
  eqAQoverQlength->SetScale(0.1);
  eqAQoverQlength->SetUseChannelGroupsWeights(kTRUE);
  VZEROCconfQoverQlength->AddCorrectionOnInputData(eqAQoverQlength);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROCconfQoverQlength->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  alignAQoverQlength = new AliQnCorrectionsQnVectorAlignment();
  alignAQoverQlength->SetHarmonicNumberForAlignment(2);
  alignAQoverQlength->SetReferenceConfigurationForAlignment("TPCQoverQlength");
  VZEROCconfQoverQlength->AddCorrectionOnQnVector(alignAQoverQlength);
  /* lets configrure the QA histograms */
  VZEROCconfQoverQlength->SetQACentralityVar(varForEventMultiplicity);
  VZEROCconfQoverQlength->SetQAMultiplicityAxis(100, 0.0, 500.0);
  /* let's configure the twist and rescale correction step */
  twScaleAQoverQlength = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleAQoverQlength->SetApplyTwist(kTRUE);
  twScaleAQoverQlength->SetApplyRescale(kTRUE);
  twScaleAQoverQlength->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleAQoverQlength->SetReferenceConfigurationsForTwistAndRescale("TPCQoverQlength", "VZEROAQoverQlength");
  /* now we add it to the detector configuration */
  VZEROCconfQoverQlength->AddCorrectionOnQnVector(twScaleAQoverQlength);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROCconfQoverQlength);

  /* the full VZERO detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* VZEROconfQoverQlength =
   new AliQnCorrectionsDetectorConfigurationChannels("VZEROQoverQlength", CorrEventClasses,
                            64, /* number of channels */
                            4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROconfQoverQlength->SetChannelsScheme(VZEROchannels[2], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROconfQoverQlength->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverQlength);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization* eqQoverQlength = new AliQnCorrectionsInputGainEqualization();
  eqQoverQlength->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqQoverQlength->SetShift(1.0);
  eqQoverQlength->SetScale(0.1);
  eqQoverQlength->SetUseChannelGroupsWeights(kTRUE);
  VZEROconfQoverQlength->AddCorrectionOnInputData(eqQoverQlength);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROconfQoverQlength->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* alignQoverQlength = new AliQnCorrectionsQnVectorAlignment();
  alignQoverQlength->SetHarmonicNumberForAlignment(2);
  alignQoverQlength->SetReferenceConfigurationForAlignment("TPCQoverQlength");
  VZEROconfQoverQlength->AddCorrectionOnQnVector(alignQoverQlength);
  /* lets configrure the QA histograms */
  VZEROconfQoverQlength->SetQACentralityVar(varForEventMultiplicity);
  VZEROconfQoverQlength->SetQAMultiplicityAxis(100, 0.0, 500.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleQoverQlength = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleQoverQlength->SetApplyTwist(kTRUE);
  twScaleQoverQlength->SetApplyRescale(kTRUE);
  twScaleQoverQlength->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleQoverQlength->SetReferenceConfigurationsForTwistAndRescale("TPCQoverQlength", "VZEROCQoverQlength");
  /* now we add it to the detector configuration */
  VZEROconfQoverQlength->AddCorrectionOnQnVector(twScaleQoverQlength);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROconfQoverQlength);

  /* the VZEROA detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* VZEROAconfQoverSqrtM =
   new AliQnCorrectionsDetectorConfigurationChannels("VZEROAQoverSqrtM", CorrEventClasses,
                            64, /* number of channels */
                            4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROAconfQoverSqrtM->SetChannelsScheme(VZEROchannels[0], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROAconfQoverSqrtM->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverSqrtM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization* eqAQoverSqrtM = new AliQnCorrectionsInputGainEqualization();
  eqAQoverSqrtM->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqAQoverSqrtM->SetShift(1.0);
  eqAQoverSqrtM->SetScale(0.1);
  eqAQoverSqrtM->SetUseChannelGroupsWeights(kTRUE);
  VZEROAconfQoverSqrtM->AddCorrectionOnInputData(eqAQoverSqrtM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROAconfQoverSqrtM->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* alignAQoverSqrtM = new AliQnCorrectionsQnVectorAlignment();
  alignAQoverSqrtM->SetHarmonicNumberForAlignment(2);
  alignAQoverSqrtM->SetReferenceConfigurationForAlignment("TPCQoverSqrtM");
  VZEROAconfQoverSqrtM->AddCorrectionOnQnVector(alignAQoverSqrtM);
  /* lets configrure the QA histograms */
  VZEROAconfQoverSqrtM->SetQACentralityVar(varForEventMultiplicity);
  VZEROAconfQoverSqrtM->SetQAMultiplicityAxis(100, 0.0, 500.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleAQoverSqrtM = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleAQoverSqrtM->SetApplyTwist(kTRUE);
  twScaleAQoverSqrtM->SetApplyRescale(kTRUE);
  twScaleAQoverSqrtM->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleAQoverSqrtM->SetReferenceConfigurationsForTwistAndRescale("TPCQoverSqrtM", "VZEROCQoverSqrtM");
  /* now we add it to the detector configuration */
  VZEROAconfQoverSqrtM->AddCorrectionOnQnVector(twScaleAQoverSqrtM);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROAconfQoverSqrtM);

  /* the VZEROC detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* VZEROCconfQoverSqrtM =
   new AliQnCorrectionsDetectorConfigurationChannels("VZEROCQoverSqrtM", CorrEventClasses,
                            64, /* number of channels */
                            4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROCconfQoverSqrtM->SetChannelsScheme(VZEROchannels[1], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROCconfQoverSqrtM->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverSqrtM);
  /* lets configure the equalization of input data */
  eqAQoverSqrtM = new AliQnCorrectionsInputGainEqualization();
  eqAQoverSqrtM->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqAQoverSqrtM->SetShift(1.0);
  eqAQoverSqrtM->SetScale(0.1);
  eqAQoverSqrtM->SetUseChannelGroupsWeights(kTRUE);
  VZEROCconfQoverSqrtM->AddCorrectionOnInputData(eqAQoverSqrtM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROCconfQoverSqrtM->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  alignAQoverSqrtM = new AliQnCorrectionsQnVectorAlignment();
  alignAQoverSqrtM->SetHarmonicNumberForAlignment(2);
  alignAQoverSqrtM->SetReferenceConfigurationForAlignment("TPCQoverSqrtM");
  VZEROCconfQoverSqrtM->AddCorrectionOnQnVector(alignAQoverSqrtM);
  /* lets configrure the QA histograms */
  VZEROCconfQoverSqrtM->SetQACentralityVar(varForEventMultiplicity);
  VZEROCconfQoverSqrtM->SetQAMultiplicityAxis(100, 0.0, 500.0);
  /* let's configure the twist and rescale correction step */
  twScaleAQoverSqrtM = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleAQoverSqrtM->SetApplyTwist(kTRUE);
  twScaleAQoverSqrtM->SetApplyRescale(kTRUE);
  twScaleAQoverSqrtM->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleAQoverSqrtM->SetReferenceConfigurationsForTwistAndRescale("TPCQoverSqrtM", "VZEROAQoverSqrtM");
  /* now we add it to the detector configuration */
  VZEROCconfQoverSqrtM->AddCorrectionOnQnVector(twScaleAQoverSqrtM);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROCconfQoverSqrtM);

  /* the full VZERO detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* VZEROconfQoverSqrtM =
   new AliQnCorrectionsDetectorConfigurationChannels("VZEROQoverSqrtM", CorrEventClasses,
                            64, /* number of channels */
                            4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROconfQoverSqrtM->SetChannelsScheme(VZEROchannels[2], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROconfQoverSqrtM->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverSqrtM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization* eqQoverSqrtM = new AliQnCorrectionsInputGainEqualization();
  eqQoverSqrtM->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqQoverSqrtM->SetShift(1.0);
  eqQoverSqrtM->SetScale(0.1);
  eqQoverSqrtM->SetUseChannelGroupsWeights(kTRUE);
  VZEROconfQoverSqrtM->AddCorrectionOnInputData(eqQoverSqrtM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROconfQoverSqrtM->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* alignQoverQoverSqrtM = new AliQnCorrectionsQnVectorAlignment();
  alignQoverQoverSqrtM->SetHarmonicNumberForAlignment(2);
  alignQoverQoverSqrtM->SetReferenceConfigurationForAlignment("TPCQoverSqrtM");
  VZEROconfQoverSqrtM->AddCorrectionOnQnVector(alignQoverQoverSqrtM);
  /* lets configrure the QA histograms */
  VZEROconfQoverSqrtM->SetQACentralityVar(varForEventMultiplicity);
  VZEROconfQoverSqrtM->SetQAMultiplicityAxis(100, 0.0, 500.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleQoverSqrtM = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleQoverSqrtM->SetApplyTwist(kTRUE);
  twScaleQoverSqrtM->SetApplyRescale(kTRUE);
  twScaleQoverSqrtM->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleQoverSqrtM->SetReferenceConfigurationsForTwistAndRescale("TPCQoverSqrtM", "VZEROCQoverM");
  /* now we add it to the detector configuration */
  VZEROconfQoverSqrtM->AddCorrectionOnQnVector(twScaleQoverSqrtM);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROconfQoverSqrtM);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(VZERO);
}

/**
 * Add the TPC configuration. Copied directly from ``AddTaskFlowQnVectorCorrectionsNewDetConfig.h``.
 */
void FlowVectorCorrections::AddTPC(AliAnalysisTaskFlowVectorCorrections* task, AliQnCorrectionsManager* QnManager,
                  AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity)
{
  /////////////// Add TPC subdetectors ///////////////////

  //-----------------------------------------------------------
  // Our event classes for TPC
  //
  const Int_t nTPCdim = 2;
  AliQnCorrectionsEventClassVariablesSet* CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nTPCdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4 }, { -7.0, 1 }, { 7.0, 8 }, { 10.0, 1 } };
  Double_t Ctbinning[][2] = { { 0.0, 2 }, { 100.0, 100 } };
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   AliQnCorrectionsVarManagerTask::kVtxZ, task->VarName(AliQnCorrectionsVarManagerTask::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   varForEventMultiplicity, Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the TPC  detector */
  AliQnCorrectionsDetector* TPC = new AliQnCorrectionsDetector("TPC", AliQnCorrectionsVarManagerTask::kTPC);

  /* the TPC detector configuration */
  AliQnCorrectionsDetectorConfigurationTracks* TPCconf = new AliQnCorrectionsDetectorConfigurationTracks(
   "TPCQoverM", CorrEventClasses, 4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  TPCconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TPCconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector twist correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScale = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScale->SetApplyTwist(kTRUE);
  twScale->SetApplyRescale(kFALSE);
  twScale->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
  TPCconf->AddCorrectionOnQnVector(twScale);

  /* define the cuts to apply */
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD = mgr->GetInputEventHandler()->IsA() == AliESDInputHandler::Class();
  AliQnCorrectionsCutsSet* cutsTPC = new AliQnCorrectionsCutsSet();
  if (!isESD) {
    cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFilterBitMask768, 0.5, 1.5));
    cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.8));
    cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
  } else {
    Bool_t UseTPConlyTracks = kFALSE; // Use of TPC standalone tracks or Global tracks (only for ESD analysis)
    task->SetUseTPCStandaloneTracks(UseTPConlyTracks);
    if (UseTPConlyTracks) {
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -3.0, 3.0));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -3.0, 3.0));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.8));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCnclsIter1, 70.0, 161.0));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2Iter1, 0.2, 4.0));
    } else {
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -0.3, 0.3));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -0.3, 0.3));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.8));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCncls, 70.0, 161.0));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2, 0.2, 4.0));
    }
  }
  TPCconf->SetCuts(cutsTPC);

  /* add the configuration to the detector */
  TPC->AddDetectorConfiguration(TPCconf);

  /* the TPC -- negative eta -- detector configuration */
  AliQnCorrectionsDetectorConfigurationTracks* TPCNegEtaconf = new AliQnCorrectionsDetectorConfigurationTracks(
   "TPCNegEtaQoverM", CorrEventClasses, 4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  TPCNegEtaconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TPCNegEtaconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector twist correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleNegEta = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleNegEta->SetApplyTwist(kTRUE);
  twScaleNegEta->SetApplyRescale(kFALSE);
  twScaleNegEta->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
  TPCNegEtaconf->AddCorrectionOnQnVector(twScaleNegEta);

  /* define the cuts to apply */
  isESD = mgr->GetInputEventHandler()->IsA() == AliESDInputHandler::Class();
  AliQnCorrectionsCutsSet* cutsTPCNegEta = new AliQnCorrectionsCutsSet();
  if (!isESD) {
    cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFilterBitMask768, 0.5, 1.5));
    cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.));
    cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
  } else {
    Bool_t UseTPConlyTracks = kFALSE; // Use of TPC standalone tracks or Global tracks (only for ESD analysis)
    task->SetUseTPCStandaloneTracks(UseTPConlyTracks);
    if (UseTPConlyTracks) {
      cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -3.0, 3.0));
      cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -3.0, 3.0));
      cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.));
      cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCNegEta->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCnclsIter1, 70.0, 161.0));
      cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2Iter1, 0.2, 4.0));
    } else {
      cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -0.3, 0.3));
      cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -0.3, 0.3));
      cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.));
      cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCncls, 70.0, 161.0));
      cutsTPCNegEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2, 0.2, 4.0));
    }
  }
  TPCNegEtaconf->SetCuts(cutsTPCNegEta);

  /* add the configuration to the detector */
  TPC->AddDetectorConfiguration(TPCNegEtaconf);

  /* the TPC -- negative eta -- detector configuration */
  AliQnCorrectionsDetectorConfigurationTracks* TPCPosEtaconf = new AliQnCorrectionsDetectorConfigurationTracks(
   "TPCPosEtaQoverM", CorrEventClasses, 4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  TPCPosEtaconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TPCPosEtaconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector twist correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScalePosEta = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScalePosEta->SetApplyTwist(kTRUE);
  twScalePosEta->SetApplyRescale(kFALSE);
  twScalePosEta->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
  TPCPosEtaconf->AddCorrectionOnQnVector(twScalePosEta);

  /* define the cuts to apply */
  isESD = mgr->GetInputEventHandler()->IsA() == AliESDInputHandler::Class();
  AliQnCorrectionsCutsSet* cutsTPCPosEta = new AliQnCorrectionsCutsSet();
  if (!isESD) {
    cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFilterBitMask768, 0.5, 1.5));
    cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, 0., 0.8));
    cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
  } else {
    Bool_t UseTPConlyTracks = kFALSE; // Use of TPC standalone tracks or Global tracks (only for ESD analysis)
    task->SetUseTPCStandaloneTracks(UseTPConlyTracks);
    if (UseTPConlyTracks) {
      cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -3.0, 3.0));
      cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -3.0, 3.0));
      cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, 0., 0.8));
      cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCPosEta->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCnclsIter1, 70.0, 161.0));
      cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2Iter1, 0.2, 4.0));
    } else {
      cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -0.3, 0.3));
      cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -0.3, 0.3));
      cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, 0., 0.8));
      cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCncls, 70.0, 161.0));
      cutsTPCPosEta->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2, 0.2, 4.0));
    }
  }
  TPCPosEtaconf->SetCuts(cutsTPCPosEta);

  /* add the configuration to the detector */
  TPC->AddDetectorConfiguration(TPCPosEtaconf);

  /* the TPC detector configuration */
  AliQnCorrectionsDetectorConfigurationTracks* TPCconfQoverQlength = new AliQnCorrectionsDetectorConfigurationTracks(
   "TPCQoverQlength", CorrEventClasses, 4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  TPCconfQoverQlength->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverQlength);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TPCconfQoverQlength->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector twist correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleQoverQlength = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleQoverQlength->SetApplyTwist(kTRUE);
  twScaleQoverQlength->SetApplyRescale(kFALSE);
  twScaleQoverQlength->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
  TPCconfQoverQlength->AddCorrectionOnQnVector(twScaleQoverQlength);

  /* define the cuts to apply */
  isESD = mgr->GetInputEventHandler()->IsA() == AliESDInputHandler::Class();
  AliQnCorrectionsCutsSet* cutsTPCQoverQlength = new AliQnCorrectionsCutsSet();
  if (!isESD) {
    cutsTPCQoverQlength->Add(
     new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFilterBitMask768, 0.5, 1.5));
    cutsTPCQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.8));
    cutsTPCQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
  } else {
    Bool_t UseTPConlyTracks = kFALSE; // Use of TPC standalone tracks or Global tracks (only for ESD analysis)
    task->SetUseTPCStandaloneTracks(UseTPConlyTracks);
    if (UseTPConlyTracks) {
      cutsTPCQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -3.0, 3.0));
      cutsTPCQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -3.0, 3.0));
      cutsTPCQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.8));
      cutsTPCQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCnclsIter1, 70.0, 161.0));
      cutsTPCQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2Iter1, 0.2, 4.0));
    } else {
      cutsTPCQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -0.3, 0.3));
      cutsTPCQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -0.3, 0.3));
      cutsTPCQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.8));
      cutsTPCQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCncls, 70.0, 161.0));
      cutsTPCQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2, 0.2, 4.0));
    }
  }
  TPCconfQoverQlength->SetCuts(cutsTPCQoverQlength);

  /* add the configuration to the detector */
  TPC->AddDetectorConfiguration(TPCconfQoverQlength);

  /* the TPC eta < 0 detector configuration */
  AliQnCorrectionsDetectorConfigurationTracks* TPCNegEtaconfQoverQlength =
   new AliQnCorrectionsDetectorConfigurationTracks("TPCNegEtaQoverQlength", CorrEventClasses,
                           4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  TPCNegEtaconfQoverQlength->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverQlength);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TPCNegEtaconfQoverQlength->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector twist correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleNegEtaQoverQlength = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleNegEtaQoverQlength->SetApplyTwist(kTRUE);
  twScaleNegEtaQoverQlength->SetApplyRescale(kFALSE);
  twScaleNegEtaQoverQlength->SetTwistAndRescaleMethod(
   AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
  TPCNegEtaconfQoverQlength->AddCorrectionOnQnVector(twScaleNegEtaQoverQlength);

  /* define the cuts to apply */
  isESD = mgr->GetInputEventHandler()->IsA() == AliESDInputHandler::Class();
  AliQnCorrectionsCutsSet* cutsTPCNegEtaQoverQlength = new AliQnCorrectionsCutsSet();
  if (!isESD) {
    cutsTPCNegEtaQoverQlength->Add(
     new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFilterBitMask768, 0.5, 1.5));
    cutsTPCNegEtaQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.));
    cutsTPCNegEtaQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
  } else {
    Bool_t UseTPConlyTracks = kFALSE; // Use of TPC standalone tracks or Global tracks (only for ESD analysis)
    task->SetUseTPCStandaloneTracks(UseTPConlyTracks);
    if (UseTPConlyTracks) {
      cutsTPCNegEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -3.0, 3.0));
      cutsTPCNegEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -3.0, 3.0));
      cutsTPCNegEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.));
      cutsTPCNegEtaQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCNegEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCnclsIter1, 70.0, 161.0));
      cutsTPCNegEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2Iter1, 0.2, 4.0));
    } else {
      cutsTPCNegEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -0.3, 0.3));
      cutsTPCNegEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -0.3, 0.3));
      cutsTPCNegEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.));
      cutsTPCNegEtaQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCNegEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCncls, 70.0, 161.0));
      cutsTPCNegEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2, 0.2, 4.0));
    }
  }
  TPCNegEtaconfQoverQlength->SetCuts(cutsTPCNegEtaQoverQlength);
  TPC->AddDetectorConfiguration(TPCNegEtaconfQoverQlength);

  /* the TPC eta > 0 detector configuration */
  AliQnCorrectionsDetectorConfigurationTracks* TPCPosEtaconfQoverQlength =
   new AliQnCorrectionsDetectorConfigurationTracks("TPCPosEtaQoverQlength", CorrEventClasses,
                           4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  TPCPosEtaconfQoverQlength->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverQlength);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TPCPosEtaconfQoverQlength->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector twist correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScalePosEtaQoverQlength = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScalePosEtaQoverQlength->SetApplyTwist(kTRUE);
  twScalePosEtaQoverQlength->SetApplyRescale(kFALSE);
  twScalePosEtaQoverQlength->SetTwistAndRescaleMethod(
   AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
  TPCPosEtaconfQoverQlength->AddCorrectionOnQnVector(twScalePosEtaQoverQlength);

  /* define the cuts to apply */
  isESD = mgr->GetInputEventHandler()->IsA() == AliESDInputHandler::Class();
  AliQnCorrectionsCutsSet* cutsTPCPosEtaQoverQlength = new AliQnCorrectionsCutsSet();
  if (!isESD) {
    cutsTPCPosEtaQoverQlength->Add(
     new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFilterBitMask768, 0.5, 1.5));
    cutsTPCPosEtaQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, 0., 0.8));
    cutsTPCPosEtaQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
  } else {
    Bool_t UseTPConlyTracks = kFALSE; // Use of TPC standalone tracks or Global tracks (only for ESD analysis)
    task->SetUseTPCStandaloneTracks(UseTPConlyTracks);
    if (UseTPConlyTracks) {
      cutsTPCPosEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -3.0, 3.0));
      cutsTPCPosEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -3.0, 3.0));
      cutsTPCPosEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, 0., 0.8));
      cutsTPCPosEtaQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCPosEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCnclsIter1, 70.0, 161.0));
      cutsTPCPosEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2Iter1, 0.2, 4.0));
    } else {
      cutsTPCPosEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -0.3, 0.3));
      cutsTPCPosEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -0.3, 0.3));
      cutsTPCPosEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, 0., 0.8));
      cutsTPCPosEtaQoverQlength->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCPosEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCncls, 70.0, 161.0));
      cutsTPCPosEtaQoverQlength->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2, 0.2, 4.0));
    }
  }
  TPCPosEtaconfQoverQlength->SetCuts(cutsTPCPosEtaQoverQlength);

  /* add the configuration to the detector */
  TPC->AddDetectorConfiguration(TPCPosEtaconfQoverQlength);

  /* the TPC detector configuration */
  AliQnCorrectionsDetectorConfigurationTracks* TPCconfQoverSqrtM = new AliQnCorrectionsDetectorConfigurationTracks(
   "TPCQoverSqrtM", CorrEventClasses, 4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  TPCconfQoverSqrtM->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverSqrtM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TPCconfQoverSqrtM->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector twist correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleQoverSqrtM = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleQoverSqrtM->SetApplyTwist(kTRUE);
  twScaleQoverSqrtM->SetApplyRescale(kFALSE);
  twScaleQoverSqrtM->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
  TPCconfQoverSqrtM->AddCorrectionOnQnVector(twScaleQoverSqrtM);

  /* define the cuts to apply */
  isESD = mgr->GetInputEventHandler()->IsA() == AliESDInputHandler::Class();
  AliQnCorrectionsCutsSet* cutsTPCQoverSqrtM = new AliQnCorrectionsCutsSet();
  if (!isESD) {
    cutsTPCQoverSqrtM->Add(
     new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFilterBitMask768, 0.5, 1.5));
    cutsTPCQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.8));
    cutsTPCQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
  } else {
    Bool_t UseTPConlyTracks = kFALSE; // Use of TPC standalone tracks or Global tracks (only for ESD analysis)
    task->SetUseTPCStandaloneTracks(UseTPConlyTracks);
    if (UseTPConlyTracks) {
      cutsTPCQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -3.0, 3.0));
      cutsTPCQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -3.0, 3.0));
      cutsTPCQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.8));
      cutsTPCQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCnclsIter1, 70.0, 161.0));
      cutsTPCQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2Iter1, 0.2, 4.0));
    } else {
      cutsTPCQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -0.3, 0.3));
      cutsTPCQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -0.3, 0.3));
      cutsTPCQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.8));
      cutsTPCQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCncls, 70.0, 161.0));
      cutsTPCQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2, 0.2, 4.0));
    }
  }
  TPCconfQoverSqrtM->SetCuts(cutsTPCQoverSqrtM);

  /* add the configuration to the detector */
  TPC->AddDetectorConfiguration(TPCconfQoverSqrtM);

  /* the TPC detector configuration */
  AliQnCorrectionsDetectorConfigurationTracks* TPCNegEtaconfQoverSqrtM =
   new AliQnCorrectionsDetectorConfigurationTracks("TPCNegEtaQoverSqrtM", CorrEventClasses,
                           4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  TPCNegEtaconfQoverSqrtM->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverSqrtM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TPCNegEtaconfQoverSqrtM->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector twist correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleNegEtaQoverSqrtM = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleNegEtaQoverSqrtM->SetApplyTwist(kTRUE);
  twScaleNegEtaQoverSqrtM->SetApplyRescale(kFALSE);
  twScaleNegEtaQoverSqrtM->SetTwistAndRescaleMethod(
   AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
  TPCNegEtaconfQoverSqrtM->AddCorrectionOnQnVector(twScaleNegEtaQoverSqrtM);

  /* define the cuts to apply */
  isESD = mgr->GetInputEventHandler()->IsA() == AliESDInputHandler::Class();
  AliQnCorrectionsCutsSet* cutsTPCNegEtaQoverSqrtM = new AliQnCorrectionsCutsSet();
  if (!isESD) {
    cutsTPCNegEtaQoverSqrtM->Add(
     new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFilterBitMask768, 0.5, 1.5));
    cutsTPCNegEtaQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.));
    cutsTPCNegEtaQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
  } else {
    Bool_t UseTPConlyTracks = kFALSE; // Use of TPC standalone tracks or Global tracks (only for ESD analysis)
    task->SetUseTPCStandaloneTracks(UseTPConlyTracks);
    if (UseTPConlyTracks) {
      cutsTPCNegEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -3.0, 3.0));
      cutsTPCNegEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -3.0, 3.0));
      cutsTPCNegEtaQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.));
      cutsTPCNegEtaQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCNegEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCnclsIter1, 70.0, 161.0));
      cutsTPCNegEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2Iter1, 0.2, 4.0));
    } else {
      cutsTPCNegEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -0.3, 0.3));
      cutsTPCNegEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -0.3, 0.3));
      cutsTPCNegEtaQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, -0.8, 0.));
      cutsTPCNegEtaQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCNegEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCncls, 70.0, 161.0));
      cutsTPCNegEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2, 0.2, 4.0));
    }
  }
  TPCNegEtaconfQoverSqrtM->SetCuts(cutsTPCNegEtaQoverSqrtM);

  /* add the configuration to the detector */
  TPC->AddDetectorConfiguration(TPCNegEtaconfQoverSqrtM);

  /* the TPC detector configuration */
  AliQnCorrectionsDetectorConfigurationTracks* TPCPosEtaconfQoverSqrtM =
   new AliQnCorrectionsDetectorConfigurationTracks("TPCPosEtaQoverSqrtM", CorrEventClasses,
                           4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  TPCPosEtaconfQoverSqrtM->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverSqrtM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TPCPosEtaconfQoverSqrtM->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector twist correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScalePosEtaQoverSqrtM = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScalePosEtaQoverSqrtM->SetApplyTwist(kTRUE);
  twScalePosEtaQoverSqrtM->SetApplyRescale(kFALSE);
  twScalePosEtaQoverSqrtM->SetTwistAndRescaleMethod(
   AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
  TPCPosEtaconfQoverSqrtM->AddCorrectionOnQnVector(twScalePosEtaQoverSqrtM);

  /* define the cuts to apply */
  isESD = mgr->GetInputEventHandler()->IsA() == AliESDInputHandler::Class();
  AliQnCorrectionsCutsSet* cutsTPCPosEtaQoverSqrtM = new AliQnCorrectionsCutsSet();
  if (!isESD) {
    cutsTPCPosEtaQoverSqrtM->Add(
     new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFilterBitMask768, 0.5, 1.5));
    cutsTPCPosEtaQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, 0., 0.8));
    cutsTPCPosEtaQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
  } else {
    Bool_t UseTPConlyTracks = kFALSE; // Use of TPC standalone tracks or Global tracks (only for ESD analysis)
    task->SetUseTPCStandaloneTracks(UseTPConlyTracks);
    if (UseTPConlyTracks) {
      cutsTPCPosEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -3.0, 3.0));
      cutsTPCPosEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -3.0, 3.0));
      cutsTPCPosEtaQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, 0., 0.8));
      cutsTPCPosEtaQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCPosEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCnclsIter1, 70.0, 161.0));
      cutsTPCPosEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2Iter1, 0.2, 4.0));
    } else {
      cutsTPCPosEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY, -0.3, 0.3));
      cutsTPCPosEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ, -0.3, 0.3));
      cutsTPCPosEtaQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta, 0., 0.8));
      cutsTPCPosEtaQoverSqrtM->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt, 0.2, 5.));
      cutsTPCPosEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCncls, 70.0, 161.0));
      cutsTPCPosEtaQoverSqrtM->Add(
       new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2, 0.2, 4.0));
    }
  }
  TPCPosEtaconfQoverSqrtM->SetCuts(cutsTPCPosEtaQoverSqrtM);

  /* add the configuration to the detector */
  TPC->AddDetectorConfiguration(TPCPosEtaconfQoverSqrtM);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(TPC);
}

/**
 * Add the SPD configuration. Copied directly from ``AddTaskFlowQnVectorCorrectionsNewDetConfig.h``.
 */
void FlowVectorCorrections::AddSPD(AliAnalysisTaskFlowVectorCorrections* task, AliQnCorrectionsManager* QnManager,
                  AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity)
{
  /////////////// Add SPD subdetectors ///////////////////

  //-----------------------------------------------------------
  // Our event classes for SPD
  //
  const Int_t nSPDdim = 2;
  AliQnCorrectionsEventClassVariablesSet* CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nSPDdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4 }, { -7.0, 1 }, { 7.0, 8 }, { 10.0, 1 } };
  Double_t Ctbinning[][2] = { { 0.0, 2 }, { 100.0, 100 } };
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   AliQnCorrectionsVarManagerTask::kVtxZ, task->VarName(AliQnCorrectionsVarManagerTask::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   varForEventMultiplicity, Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the SPD detector */
  AliQnCorrectionsDetector* SPD = new AliQnCorrectionsDetector("SPD", AliQnCorrectionsVarManagerTask::kSPD);

  /* the SPD detector configuration */
  AliQnCorrectionsDetectorConfigurationTracks* SPDconf = new AliQnCorrectionsDetectorConfigurationTracks(
   "SPD", CorrEventClasses, 4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  SPDconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  SPDconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector twist correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScale = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScale->SetApplyTwist(kTRUE);
  twScale->SetApplyRescale(kFALSE);
  twScale->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
  SPDconf->AddCorrectionOnQnVector(twScale);

  /* add the configuration to the detector */
  SPD->AddDetectorConfiguration(SPDconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(SPD);
}

/**
 * Add the TZERO configuration. Copied directly from ``AddTaskFlowQnVectorCorrectionsNewDetConfig.h``.
 */
void FlowVectorCorrections::AddTZERO(AliAnalysisTaskFlowVectorCorrections* task, AliQnCorrectionsManager* QnManager,
                   AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity)
{
  /////////////// Add TZERO subdetectors ///////////////////

  Bool_t TZEROchannels[2][24];
  for (Int_t iv0 = 0; iv0 < 2; iv0++)
    for (Int_t ich = 0; ich < 24; ich++)
      TZEROchannels[iv0][ich] = kFALSE;

  for (Int_t ich = 12; ich < 24; ich++)
    TZEROchannels[0][ich] = kTRUE; // channel list: value 1 if channel should be used
  for (Int_t ich = 0; ich < 12; ich++)
    TZEROchannels[1][ich] = kTRUE;

  Int_t channelGroups[24];
  for (Int_t ich = 0; ich < 24; ich++)
    channelGroups[ich] = Int_t(ich / 12);

  //-----------------------------------------------------------
  // Our event classes for TZERO
  //
  const Int_t nTZEROdim = 2;
  AliQnCorrectionsEventClassVariablesSet* CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nTZEROdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4 }, { -7.0, 1 }, { 7.0, 8 }, { 10.0, 1 } };
  Double_t Ctbinning[][2] = { { 0.0, 2 }, { 100.0, 100 } };
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   AliQnCorrectionsVarManagerTask::kVtxZ, task->VarName(AliQnCorrectionsVarManagerTask::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   varForEventMultiplicity, Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the TZERO detector */
  AliQnCorrectionsDetector* TZERO = new AliQnCorrectionsDetector("TZERO", AliQnCorrectionsVarManagerTask::kTZERO);

  /* the TZEROA detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* TZEROAconf =
   new AliQnCorrectionsDetectorConfigurationChannels("TZEROA", CorrEventClasses, 24, /* number of channels */
                            4); /* number of harmonics: 1, 2, 3 and 4 */
  TZEROAconf->SetChannelsScheme(TZEROchannels[0], channelGroups);
  /* let's configure the Q vector calibration */
  TZEROAconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization* eqA = new AliQnCorrectionsInputGainEqualization();
  eqA->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqA->SetShift(1.0);
  eqA->SetScale(0.1);
  eqA->SetUseChannelGroupsWeights(kFALSE);
  TZEROAconf->AddCorrectionOnInputData(eqA);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TZEROAconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* alignA = new AliQnCorrectionsQnVectorAlignment();
  alignA->SetHarmonicNumberForAlignment(2);
  alignA->SetReferenceConfigurationForAlignment("TPC");
  TZEROAconf->AddCorrectionOnQnVector(alignA);
  /* let's configure the QA histograms */
  TZEROAconf->SetQACentralityVar(varForEventMultiplicity);
  TZEROAconf->SetQAMultiplicityAxis(100, 0.0, 150.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleA = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleA->SetApplyTwist(kTRUE);
  twScaleA->SetApplyRescale(kTRUE);
  twScaleA->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleA->SetReferenceConfigurationsForTwistAndRescale("TPC", "TZEROC");
  /* now we add it to the detector configuration */
  TZEROAconf->AddCorrectionOnQnVector(twScaleA);

  /* add the configuration to the detector */
  TZERO->AddDetectorConfiguration(TZEROAconf);

  /* the TZEROC detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* TZEROCconf =
   new AliQnCorrectionsDetectorConfigurationChannels("TZEROC", CorrEventClasses, 24, /* number of channels */
                            4); /* number of harmonics: 1, 2, 3 and 4 */
  TZEROCconf->SetChannelsScheme(TZEROchannels[1], channelGroups);
  /* let's configure the Q vector calibration */
  TZEROCconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization* eqC = new AliQnCorrectionsInputGainEqualization();
  eqC->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqC->SetShift(1.0);
  eqC->SetScale(0.1);
  eqC->SetUseChannelGroupsWeights(kFALSE);
  TZEROCconf->AddCorrectionOnInputData(eqC);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TZEROCconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* alignC = new AliQnCorrectionsQnVectorAlignment();
  alignC->SetHarmonicNumberForAlignment(2);
  alignC->SetReferenceConfigurationForAlignment("TPC");
  TZEROCconf->AddCorrectionOnQnVector(alignC);
  /* let's configure the QA histograms */
  TZEROCconf->SetQACentralityVar(varForEventMultiplicity);
  TZEROCconf->SetQAMultiplicityAxis(100, 0.0, 150.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleC = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleC->SetApplyTwist(kTRUE);
  twScaleC->SetApplyRescale(kTRUE);
  twScaleC->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleC->SetReferenceConfigurationsForTwistAndRescale("TPC", "TZEROA");
  /* now we add it to the detector configuration */
  TZEROCconf->AddCorrectionOnQnVector(twScaleC);

  /* add the configuration to the detector */
  TZERO->AddDetectorConfiguration(TZEROCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(TZERO);
}

/**
 * Add the ZDC configuration. Copied directly from ``AddTaskFlowQnVectorCorrectionsNewDetConfig.h``.
 */
void FlowVectorCorrections::AddZDC(AliAnalysisTaskFlowVectorCorrections* task, AliQnCorrectionsManager* QnManager,
                  AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity)
{
  /////////////// Add ZDC subdetectors ///////////////////

  Bool_t ZDCchannels[2][10];
  for (Int_t iv0 = 0; iv0 < 2; iv0++)
    for (Int_t ich = 0; ich < 10; ich++)
      ZDCchannels[iv0][ich] = kFALSE;

  for (Int_t ich = 6; ich < 10; ich++)
    ZDCchannels[0][ich] = kTRUE;
  for (Int_t ich = 1; ich < 5; ich++)
    ZDCchannels[1][ich] = kTRUE;

  //-----------------------------------------------------------
  // Our event classes for ZDC
  //
  const Int_t nZDCdim = 3;
  AliQnCorrectionsEventClassVariablesSet* CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nZDCdim);
  Double_t VtxXbinning[][2] = { { -0.3, 2 }, { 0.3, 10 } };
  Double_t VtxYbinning[][2] = { { -0.3, 2 }, { 0.3, 10 } };
  Double_t Ctbinning[][2] = { { 0.0, 2 }, { 100.0, 100 } };
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   AliQnCorrectionsVarManagerTask::kVtxX, task->VarName(AliQnCorrectionsVarManagerTask::kVtxX), VtxXbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   AliQnCorrectionsVarManagerTask::kVtxY, task->VarName(AliQnCorrectionsVarManagerTask::kVtxY), VtxYbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   varForEventMultiplicity, Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the ZDC detector */
  AliQnCorrectionsDetector* ZDC = new AliQnCorrectionsDetector("ZDC", AliQnCorrectionsVarManagerTask::kZDC);

  /* the ZDCA detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* ZDCAconf =
   new AliQnCorrectionsDetectorConfigurationChannels("ZDCA", CorrEventClasses, 10, /* number of channels */
                            3); /* number of harmonics: 1, 2 and 3 */
  ZDCAconf->SetChannelsScheme(ZDCchannels[0], NULL /* no groups */);
  /* let's configure the Q vector calibration */
  ZDCAconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  ZDCAconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());

  /* add the configuration to the detector */
  ZDC->AddDetectorConfiguration(ZDCAconf);

  /* the ZDCC detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* ZDCCconf =
   new AliQnCorrectionsDetectorConfigurationChannels("ZDCC", CorrEventClasses, 10, /* number of channels */
                            3); /* number of harmonics: 1, 2 and 3 */
  ZDCCconf->SetChannelsScheme(ZDCchannels[1], NULL /* no groups */);
  /* let's configure the Q vector calibration */
  ZDCCconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  ZDCCconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());

  /* add the configuration to the detector */
  ZDC->AddDetectorConfiguration(ZDCCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(ZDC);
}

/**
 * Add the FMD configuration. Copied directly from ``AddTaskFlowQnVectorCorrectionsNewDetConfig.h``.
 */
void FlowVectorCorrections::AddFMD(AliAnalysisTaskFlowVectorCorrections* task, AliQnCorrectionsManager* QnManager,
                  AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity)
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD = mgr->GetInputEventHandler()->IsA() == AliESDInputHandler::Class();
  if (isESD) {
    AliFatalGeneral("PWGJE::EMCALJetTasks::FlowVectorCorrections::AddFMD",
            "ESD analysis is not supported for the FMD.");
  }

  Bool_t FMDchannels[2][4000];
  for (Int_t iv0 = 0; iv0 < 2; iv0++)
    for (Int_t ich = 0; ich < 4000; ich++)
      FMDchannels[iv0][ich] = kFALSE;

  for (Int_t ich = 2000; ich < 4000; ich++)
    FMDchannels[0][ich] = kTRUE; // channel list: value 1 if channel should be used
  for (Int_t ich = 0; ich < 2000; ich++)
    FMDchannels[1][ich] = kTRUE;

  //-----------------------------------------------------------
  // Our event classes for FMD
  //
  const Int_t nFMDdim = 2;
  AliQnCorrectionsEventClassVariablesSet* CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nFMDdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4 }, { -7.0, 1 }, { 7.0, 8 }, { 10.0, 1 } };
  Double_t Ctbinning[][2] = { { 0.0, 2 }, { 100.0, 100 } };
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   AliQnCorrectionsVarManagerTask::kVtxZ, task->VarName(AliQnCorrectionsVarManagerTask::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   varForEventMultiplicity, Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the FMD detector */
  AliQnCorrectionsDetector* FMD = new AliQnCorrectionsDetector("FMD", AliQnCorrectionsVarManagerTask::kFMD);

  /* the FMDA detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* FMDAconf =
   new AliQnCorrectionsDetectorConfigurationChannels("FMDA", CorrEventClasses, 4000, /* number of channels */
                            4); /* number of harmonics: 1, 2 and 3 */
  FMDAconf->SetChannelsScheme(FMDchannels[0], NULL /* no groups */);
  /* let's configure the Q vector calibration */
  FMDAconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  FMDAconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* alignA = new AliQnCorrectionsQnVectorAlignment();
  alignA->SetHarmonicNumberForAlignment(2);
  alignA->SetReferenceConfigurationForAlignment("TPC");
  FMDAconf->AddCorrectionOnQnVector(alignA);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleA = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleA->SetApplyTwist(kTRUE);
  twScaleA->SetApplyRescale(kTRUE);
  twScaleA->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleA->SetReferenceConfigurationsForTwistAndRescale("TPC", "FMDC");
  /* now we add it to the detector configuration */
  FMDAconf->AddCorrectionOnQnVector(twScaleA);

  /* add the configuration to the detector */
  FMD->AddDetectorConfiguration(FMDAconf);

  /* the FMDC detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* FMDCconf =
   new AliQnCorrectionsDetectorConfigurationChannels("FMDC", CorrEventClasses, 4000, /* number of channels */
                            4); /* number of harmonics: 1, 2, 3 and 4 */
  FMDCconf->SetChannelsScheme(FMDchannels[1], NULL /* no groups */);
  /* let's configure the Q vector calibration */
  FMDCconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  FMDCconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* alignC = new AliQnCorrectionsQnVectorAlignment();
  alignC->SetHarmonicNumberForAlignment(2);
  alignC->SetReferenceConfigurationForAlignment("TPC");
  FMDCconf->AddCorrectionOnQnVector(alignC);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale* twScaleC = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleC->SetApplyTwist(kTRUE);
  twScaleC->SetApplyRescale(kTRUE);
  twScaleC->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleC->SetReferenceConfigurationsForTwistAndRescale("TPC", "FMDA");
  /* now we add it to the detector configuration */
  FMDCconf->AddCorrectionOnQnVector(twScaleC);

  /* add the configuration to the detector */
  FMD->AddDetectorConfiguration(FMDCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(FMD);
}

/**
 * Add the raw FMD configuration. Copied directly from ``AddTaskFlowQnVectorCorrectionsNewDetConfig.h``.
 */
void FlowVectorCorrections::AddRawFMD(AliAnalysisTaskFlowVectorCorrections* task, AliQnCorrectionsManager* QnManager,
                   AliQnCorrectionsVarManagerTask::Variables varForEventMultiplicity)
{
  /////////////// Add FMD subdetectors ///////////////////
  /* FMD1 and FMD2 make FMDA and FMD3 make FMDC */

  const Int_t nNoOfDetectors = 3;             ///< the number of FMD detectors
  const Int_t detectorNumber[] = { 1, 2, 3 }; ///< the number of the FMD detector
  const Int_t nNoOfRings[] = { 1, 2, 2 };     ///< the number of rings for each detector
  const Int_t ringNoOfSectors[] = { 20, 40 }; ///< ring number of sectors
  const Int_t nTotalNoOfChannels =
   ringNoOfSectors[0] * 3 + ringNoOfSectors[1] * 2; ///< three inner sectors plus two outer ones
  // const Int_t nTotalNoOfGroups     = nNoOfRings[0] + nNoOfRings[1] + nNoOfRings[2]; ///< each ring one channel
  // group
  const Int_t FMDCdetectorNumber = 3; ///< the number of the detector associated to FMDC
  Int_t nSectorId = 0;                ///< the sector id used as channel number
  Int_t nRingId = 0;                  ///< the ring id (0..4) used as group number

  Bool_t FMDchannels[2][nTotalNoOfChannels];  ///< the assignment of channels to each subdetector
  Int_t FMDchannelGroups[nTotalNoOfChannels]; ///< the group associated to each channel
  for (Int_t i = 0; i < 2; i++)
    for (Int_t c = 0; c < nTotalNoOfChannels; c++)
      FMDchannels[i][c] = kFALSE;

  for (Int_t detector = 0; detector < nNoOfDetectors; detector++) {
    for (Int_t ring = 0; ring < nNoOfRings[detector]; ring++) {
      for (Int_t sector = 0; sector < ringNoOfSectors[ring]; sector++) {
        FMDchannels[0][nSectorId] = ((detectorNumber[detector] != FMDCdetectorNumber) ? kTRUE : kFALSE);
        FMDchannels[1][nSectorId] = ((detectorNumber[detector] != FMDCdetectorNumber) ? kFALSE : kTRUE);
        FMDchannelGroups[nSectorId] = nRingId;
        nSectorId++;
      }
      nRingId++;
    }
  }

  //-----------------------------------------------------------
  // Our event classes for FMD
  //
  const Int_t nFMDdim = 2;
  AliQnCorrectionsEventClassVariablesSet* CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nFMDdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4 }, { -7.0, 1 }, { 7.0, 8 }, { 10.0, 1 } };
  Double_t Ctbinning[][2] = { { 0.0, 2 }, { 100.0, 100 } };
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   AliQnCorrectionsVarManagerTask::kVtxZ, task->VarName(AliQnCorrectionsVarManagerTask::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(
   varForEventMultiplicity, Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  AliQnCorrectionsCutsSet* cutFMDA = new AliQnCorrectionsCutsSet();
  cutFMDA->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFMDEta, 0.0, 6.0));

  AliQnCorrectionsCutsSet* cutFMDC = new AliQnCorrectionsCutsSet();
  cutFMDC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFMDEta, -6.0, 0.0));

  /* the FMD detector */
  AliQnCorrectionsDetector* FMDraw = new AliQnCorrectionsDetector("FMDraw", AliQnCorrectionsVarManagerTask::kFMDraw);

  /* the FMDAraw detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* FMDArawconf = new AliQnCorrectionsDetectorConfigurationChannels(
   "FMDAraw", CorrEventClasses, nTotalNoOfChannels, 4); /* number of harmonics: 1, 2, 3 and 4 */
  FMDArawconf->SetChannelsScheme(FMDchannels[0], FMDchannelGroups);
  /* let's configure the Q vector calibration */
  FMDArawconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization* eqA = new AliQnCorrectionsInputGainEqualization();
  eqA->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqA->SetShift(1.0);
  eqA->SetScale(0.1);
  eqA->SetUseChannelGroupsWeights(kTRUE);
  FMDArawconf->AddCorrectionOnInputData(eqA);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  FMDArawconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* alignA = new AliQnCorrectionsQnVectorAlignment();
  alignA->SetHarmonicNumberForAlignment(2);
  alignA->SetReferenceConfigurationForAlignment("TPC");
  FMDArawconf->AddCorrectionOnQnVector(alignA);
  /* and add the cuts */
  FMDArawconf->SetCuts(cutFMDA);

  /* add the configuration to the detector */
  FMDraw->AddDetectorConfiguration(FMDArawconf);

  /* the FMDCraw detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels* FMDCrawconf = new AliQnCorrectionsDetectorConfigurationChannels(
   "FMDCraw", CorrEventClasses, nTotalNoOfChannels, 4); /* number of harmonics: 1, 2, 3 and 4 */
  FMDCrawconf->SetChannelsScheme(FMDchannels[1], FMDchannelGroups);
  /* let's configure the Q vector calibration */
  FMDCrawconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization* eqC = new AliQnCorrectionsInputGainEqualization();
  eqC->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqC->SetShift(1.0);
  eqC->SetScale(0.1);
  eqC->SetUseChannelGroupsWeights(kTRUE);
  FMDCrawconf->AddCorrectionOnInputData(eqC);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  FMDCrawconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment* alignC = new AliQnCorrectionsQnVectorAlignment();
  alignC->SetHarmonicNumberForAlignment(2);
  alignC->SetReferenceConfigurationForAlignment("TPC");
  FMDCrawconf->AddCorrectionOnQnVector(alignC);
  /* and add the cuts */
  FMDCrawconf->SetCuts(cutFMDC);

  /* add the configuration to the detector */
  FMDraw->AddDetectorConfiguration(FMDCrawconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(FMDraw);
}

/**
 * Define Qn vector histograms. Copied directly from ``AddTaskFlowQnVectorCorrectionsNewDetConfig.h``.
 */
void FlowVectorCorrections::DefineHistograms(AliQnCorrectionsManager* QnManager, AliQnCorrectionsHistos* histos,
                       TString histClass)
{
  //
  // define the histograms
  //
  const Char_t* histClasses = histClass.Data();

  AliInfoGeneralStream("PWGJE::EMCALJetTasks::FlowVectorCorrections::DefineHistograms")
   << "Defining histograms ...\n";
  AliInfoGeneralStream("PWGJE::EMCALJetTasks::FlowVectorCorrections::DefineHistograms")
   << "histogram classes: " << histClass << "\n";

  // fHistosFile=new TFile(output,"RECREATE");

  TString classesStr(histClasses);
  TObjArray* arr = classesStr.Tokenize(";");

  const Int_t kNRunBins = 3000;
  Double_t runHistRange[2] = { 137000., 140000. };

  for (Int_t iclass = 0; iclass < arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    AliInfoGeneralStream("PWGJE::EMCALJetTasks::FlowVectorCorrections::DefineHistograms")
     << "hist class: " << classStr.Data() << "\n";

    // Event wise histograms
    if (classStr.Contains("Event")) {
      histos->AddHistClass(classStr.Data());
      histos->AddHistogram(classStr.Data(), "RunNo", "Run numbers;Run", kFALSE, kNRunBins, runHistRange[0],
                 runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo);
      histos->AddHistogram(classStr.Data(), "BC", "Bunch crossing;BC", kFALSE, 3000, 0., 3000.,
                 AliQnCorrectionsVarManagerTask::kBC);
      histos->AddHistogram(classStr.Data(), "IsPhysicsSelection", "Physics selection flag;;", kFALSE, 2, -0.5,
                 1.5, AliQnCorrectionsVarManagerTask::kIsPhysicsSelection, 0, 0.0, 0.0,
                 AliQnCorrectionsVarManagerTask::kNothing, 0, 0.0, 0.0,
                 AliQnCorrectionsVarManagerTask::kNothing, "off;on");

      histos->AddHistogram(classStr.Data(), "VtxZ", "Vtx Z;vtx Z (cm)", kFALSE, 300, -30.0, 30.0,
                 AliQnCorrectionsVarManagerTask::kVtxZ);
      // histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)",
      // kFALSE,300,-15.,15.,AliQnCorrectionsVarManagerTask::kVtxZ);
      histos->AddHistogram(classStr.Data(), "VtxX", "Vtx X;vtx X (cm)", kFALSE, 300, -1., 1.,
                 AliQnCorrectionsVarManagerTask::kVtxX);
      histos->AddHistogram(classStr.Data(), "VtxY", "Vtx Y;vtx Y (cm)", kFALSE, 300, -1., 1.,
                 AliQnCorrectionsVarManagerTask::kVtxY);

      histos->AddHistogram(
       classStr.Data(), "CentVZEROvsMultPVZERO",
       "Multiplicity percentile (VZERO);multiplicity VZERO (percents);centrality VZERO (percents)", kFALSE, 100,
       0.0, 100.0, AliQnCorrectionsVarManagerTask::kVZEROMultPercentile, 100, 0.0, 100.0,
       AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(), "CentVZEROvsCentSPD",
                 "Centrality(VZERO);centrality VZERO (percents);centrality SPD (percents)", kFALSE, 100,
                 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 100, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentSPD);
      histos->AddHistogram(classStr.Data(), "CentTPCvsCentSPD",
                 "Centrality(TPC);centrality TPC (percents);centrality SPD (percents)", kFALSE, 100,
                 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentTPC, 100, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentSPD);
      histos->AddHistogram(classStr.Data(), "CentTPCvsCentVZERO",
                 "Centrality(TPC);centrality TPC (percents);centrality VZERO (percents)", kFALSE, 100,
                 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentTPC, 100, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(), "CentTPCvsCentZDC",
                 "Centrality(TPC);centrality TPC (percents);centrality ZDC (percents)", kFALSE, 100,
                 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentTPC, 100, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentZDC);
      histos->AddHistogram(classStr.Data(), "CentZDCvsCentVZERO",
                 "Centrality(ZDC);centrality ZDC (percents);centrality VZERO (percents)", kFALSE, 100,
                 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentZDC, 100, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(), "CentZDCvsCentSPD",
                 "Centrality(ZDC);centrality ZDC (percents);centrality SPD (percents)", kFALSE, 100,
                 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentZDC, 100, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentSPD);

      histos->AddHistogram(classStr.Data(), "MultVZEROvsCentVZERO",
                 "Multiplicity;multiplicity VZERO;VZERO centrality", kFALSE, 100, 0.0, 32000.0,
                 AliQnCorrectionsVarManagerTask::kVZEROTotalMult, 100, 0., 100.,
                 AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(), "MultSPDvsCentPSD", "Multiplicity;SPD tracklets;SPD centrality",
                 kFALSE, 100, 0.0, 3000.0, AliQnCorrectionsVarManagerTask::kSPDntracklets, 100, 0.,
                 100., AliQnCorrectionsVarManagerTask::kCentSPD);
      histos->AddHistogram(classStr.Data(), "MultTPCvsCentSPD", "Multiplicity;TPC selected tracks;TPC centrality",
                 kFALSE, 100, 0.0, 3500.0, AliQnCorrectionsVarManagerTask::kNtracksSelected, 100, 0.,
                 100., AliQnCorrectionsVarManagerTask::kCentTPC);
      histos->AddHistogram(classStr.Data(), "MultZDCvsCentZDC", "Multiplicity;multiplicity ZDC;ZDC centrality",
                 kFALSE, 100, 0.0, 300000.0, AliQnCorrectionsVarManagerTask::kZDCTotalEnergy, 100, 0.,
                 100., AliQnCorrectionsVarManagerTask::kCentZDC);

      histos->AddHistogram(classStr.Data(), "MultTPCvsMultVZERO", "Multiplicity;tracks TPC;multiplicity VZERO",
                 kFALSE, 100, 0.0, 3500.0, AliQnCorrectionsVarManagerTask::kNtracksSelected, 100, 0.0,
                 32000.0, AliQnCorrectionsVarManagerTask::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(), "MultTPCvsMultSPD", "Multiplicity;tracklets SPD;tracks TPC", kFALSE,
                 100, 0.0, 3500.0, AliQnCorrectionsVarManagerTask::kNtracksSelected, 100, 0.0, 3000.0,
                 AliQnCorrectionsVarManagerTask::kSPDntracklets);
      histos->AddHistogram(classStr.Data(), "MultSPDvsMultVZERO", "Multiplicity;tracklets SPD;multiplicity VZERO",
                 kFALSE, 100, 0.0, 32000.0, AliQnCorrectionsVarManagerTask::kVZEROTotalMult, 100, 0.0,
                 3000.0, AliQnCorrectionsVarManagerTask::kSPDntracklets);
      histos->AddHistogram(classStr.Data(), "MultTPCvsMultZDC", "Multiplicity;tracks TPC;energy ZDC", kFALSE, 100,
                 0.0, 3500.0, AliQnCorrectionsVarManagerTask::kNtracksSelected, 100, 0.0, 300000.0,
                 AliQnCorrectionsVarManagerTask::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(), "MultVZEROvsMultZDC", "Multiplicity;multiplicity VZERO;energy ZDC",
                 kFALSE, 100, 0.0, 32000.0, AliQnCorrectionsVarManagerTask::kVZEROTotalMult, 100, 0.0,
                 300000.0, AliQnCorrectionsVarManagerTask::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(), "MultSPDvsMultZDC", "Multiplicity;tracklets SPD;energy ZDC", kFALSE,
                 100, 0.0, 3000.0, AliQnCorrectionsVarManagerTask::kSPDntracklets, 100, 0.0, 300000.0,
                 AliQnCorrectionsVarManagerTask::kZDCTotalEnergy);

      histos->AddHistogram(classStr.Data(), "MultVZERO", "Multiplicity;multiplicity VZERO", kFALSE, 320, 0.0,
                 25000.0, AliQnCorrectionsVarManagerTask::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(), "MultVZEROA", "Multiplicity;multiplicity VZEROA", kFALSE, 250, 0.0,
                 9500.0, AliQnCorrectionsVarManagerTask::kVZEROATotalMult); // 10000.0
      histos->AddHistogram(classStr.Data(), "MultVZEROC", "Multiplicity;multiplicity VZEROC", kFALSE, 250, 0.0,
                 16000.0, AliQnCorrectionsVarManagerTask::kVZEROCTotalMult); // 15000.0
      histos->AddHistogram(classStr.Data(), "MultZDC", "Multiplicity;multiplicity ZDC", kFALSE, 200, 0.0,
                 300000.0, AliQnCorrectionsVarManagerTask::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(), "MultZDCA", "Multiplicity;multiplicity ZDCA", kFALSE, 200, 0.0,
                 150000.0, AliQnCorrectionsVarManagerTask::kZDCATotalEnergy);
      histos->AddHistogram(classStr.Data(), "MultZDCC", "Multiplicity;multiplicity ZDCC", kFALSE, 200, 0.0,
                 150000.0, AliQnCorrectionsVarManagerTask::kZDCCTotalEnergy);
      histos->AddHistogram(classStr.Data(), "MultFMD1", "Multiplicity;multiplicity FMD1", kFALSE, 300, 0.0,
                 10000.0, AliQnCorrectionsVarManagerTask::kFMD1TotalMult);
      histos->AddHistogram(classStr.Data(), "MultFMD2I", "Multiplicity;multiplicity FMD2I", kFALSE, 300, 0.0,
                 10000.0, AliQnCorrectionsVarManagerTask::kFMD2ITotalMult);
      histos->AddHistogram(classStr.Data(), "MultFMD2O", "Multiplicity;multiplicity FMD2O", kFALSE, 300, 0.0,
                 10000.0, AliQnCorrectionsVarManagerTask::kFMD2OTotalMult);
      histos->AddHistogram(classStr.Data(), "MultFMD3I", "Multiplicity;multiplicity FMD3I", kFALSE, 300, 0.0,
                 10000.0, AliQnCorrectionsVarManagerTask::kFMD3ITotalMult);
      histos->AddHistogram(classStr.Data(), "MultFMD3O", "Multiplicity;multiplicity FMD3O", kFALSE, 300, 0.0,
                 10000.0, AliQnCorrectionsVarManagerTask::kFMD3OTotalMult);
      histos->AddHistogram(classStr.Data(), "MultTZEROA", "Multiplicity;multiplicity TZEROA", kFALSE, 300, 0.0,
                 3000.0, AliQnCorrectionsVarManagerTask::kTZEROATotalMult);
      histos->AddHistogram(classStr.Data(), "MultTZEROC", "Multiplicity;multiplicity TZEROC", kFALSE, 300, 0.0,
                 3000.0, AliQnCorrectionsVarManagerTask::kTZEROCTotalMult);

      histos->AddHistogram(classStr.Data(), "MultPercentVZERO",
                 "Multiplicity percentile (VZERO);multiplicity VZERO (percents)", kFALSE, 100, 0.0,
                 100.0, AliQnCorrectionsVarManagerTask::kVZEROMultPercentile);
      histos->AddHistogram(classStr.Data(), "CentVZERO", "Centrality(VZERO);centrality VZERO (percents)", kFALSE,
                 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(), "CentSPD", "Centrality(SPD);centrality SPD (percents)", kFALSE, 100,
                 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentSPD);
      histos->AddHistogram(classStr.Data(), "CentTPC", "Centrality(TPC);centrality TPC (percents)", kFALSE, 100,
                 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentTPC);
      histos->AddHistogram(classStr.Data(), "CentZDC", "Centrality(ZDC);centrality ZDC (percents)", kFALSE, 100,
                 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentZDC);

      histos->AddHistogram(classStr.Data(), "CentQuality", "Centrality quality;centrality quality", kFALSE, 100,
                 -50.5, 49.5, AliQnCorrectionsVarManagerTask::kCentQuality);
      histos->AddHistogram(classStr.Data(), "CentVZERO_Run_prof",
                 "<Centrality(VZERO)> vs run;Run; centrality VZERO (%)", kTRUE, kNRunBins,
                 runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0.0,
                 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(), "CentSPD_Run_prof",
                 "<Centrality(SPD)> vs run;Run; centrality SPD (%)", kTRUE, kNRunBins, runHistRange[0],
                 runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentSPD);
      histos->AddHistogram(classStr.Data(), "CentTPC_Run_prof",
                 "<Centrality(TPC)> vs run;Run; centrality TPC (%)", kTRUE, kNRunBins, runHistRange[0],
                 runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentTPC);
      histos->AddHistogram(classStr.Data(), "CentZDC_Run_prof",
                 "<Centrality(ZDC)> vs run;Run; centrality ZDC (%)", kTRUE, kNRunBins, runHistRange[0],
                 runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentZDC);

      histos->AddHistogram(classStr.Data(), "NV0sTotal", "Number of V0 candidates per event;# pairs", kFALSE,
                 1000, 0., 30000., AliQnCorrectionsVarManagerTask::kNV0total);
      histos->AddHistogram(classStr.Data(), "NV0sSelected", "Number of selected V0 candidates per event;# pairs",
                 kFALSE, 1000, 0., 10000., AliQnCorrectionsVarManagerTask::kNV0selected);
      histos->AddHistogram(classStr.Data(), "NPairs", "Number of candidates per event;# pairs", kFALSE, 5000, 0.,
                 5000., AliQnCorrectionsVarManagerTask::kNdielectrons);
      histos->AddHistogram(classStr.Data(), "NPairsSelected", "Number of selected pairs per event; #pairs",
                 kFALSE, 5000, 0., 5000., AliQnCorrectionsVarManagerTask::kNpairsSelected);
      histos->AddHistogram(classStr.Data(), "NTracksTotal", "Number of total tracks per event;# tracks", kFALSE,
                 1000, 0., 30000., AliQnCorrectionsVarManagerTask::kNtracksTotal);
      histos->AddHistogram(classStr.Data(), "NTracksSelected", "Number of selected tracks per event;# tracks",
                 kFALSE, 1000, 0., 30000., AliQnCorrectionsVarManagerTask::kNtracksSelected);
      histos->AddHistogram(classStr.Data(), "SPDntracklets", "SPD #tracklets; tracklets", kFALSE, 3000, -0.5,
                 2999.5, AliQnCorrectionsVarManagerTask::kSPDntracklets);
      histos->AddHistogram(classStr.Data(), "SPDnSingleClusters", "SPD #single clusters; tracklets", kFALSE, 3000,
                 -0.5, 2999.5, AliQnCorrectionsVarManagerTask::kSPDnSingleClusters);

      histos->AddHistogram(classStr.Data(), "NV0total_Run_prof", "<Number of total V0s> per run; Run; #tracks",
                 kTRUE, kNRunBins, runHistRange[0], runHistRange[1],
                 AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000.,
                 AliQnCorrectionsVarManagerTask::kNV0total);
      histos->AddHistogram(classStr.Data(), "NV0selected_Run_prof",
                 "<Number of selected V0s> per run; Run; #tracks", kTRUE, kNRunBins, runHistRange[0],
                 runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000.,
                 AliQnCorrectionsVarManagerTask::kNV0selected);
      histos->AddHistogram(classStr.Data(), "Ndielectrons_Run_prof",
                 "<Number of dielectrons> per run; Run; #tracks", kTRUE, kNRunBins, runHistRange[0],
                 runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000.,
                 AliQnCorrectionsVarManagerTask::kNdielectrons);
      histos->AddHistogram(classStr.Data(), "NpairsSelected_Run_prof",
                 "<Number of selected pairs> per run; Run; #tracks", kTRUE, kNRunBins, runHistRange[0],
                 runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000.,
                 AliQnCorrectionsVarManagerTask::kNpairsSelected);
      histos->AddHistogram(classStr.Data(), "NTracksTotal_Run_prof", "<Number of tracks> per run; Run; #tracks",
                 kTRUE, kNRunBins, runHistRange[0], runHistRange[1],
                 AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000.,
                 AliQnCorrectionsVarManagerTask::kNtracksTotal);
      histos->AddHistogram(classStr.Data(), "NTracksSelected_Run_prof",
                 "<Number of selected tracks> per run; Run; #tracks", kTRUE, kNRunBins, runHistRange[0],
                 runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000.,
                 AliQnCorrectionsVarManagerTask::kNtracksSelected);
      histos->AddHistogram(classStr.Data(), "SPDntracklets_Run_prof", "<SPD ntracklets> per run; Run; #tracks",
                 kTRUE, kNRunBins, runHistRange[0], runHistRange[1],
                 AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000.,
                 AliQnCorrectionsVarManagerTask::kSPDntracklets);

      histos->AddHistogram(classStr.Data(), "VtxZ_CentVZERO",
                 "Centrality(VZERO) vs vtx. Z;vtx Z (cm); centrality VZERO (%)", kFALSE, 300, -15., 15.,
                 AliQnCorrectionsVarManagerTask::kVtxZ, 100, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(), "VtxZ_CentSPD",
                 "Centrality(SPD) vs vtx. Z;vtx Z (cm); centrality SPD (%)", kFALSE, 300, -15., 15.,
                 AliQnCorrectionsVarManagerTask::kVtxZ, 100, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentSPD);
      histos->AddHistogram(classStr.Data(), "VtxZ_CentTPC",
                 "Centrality(TPC) vs vtx. Z;vtx Z (cm); centrality TPC (%)", kFALSE, 300, -15., 15.,
                 AliQnCorrectionsVarManagerTask::kVtxZ, 100, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentTPC);
      continue;
    } // end if className contains "Event"

    // Offline trigger histograms
    if (classStr.Contains("OfflineTriggers")) {
      histos->AddHistClass(classStr.Data());

      TString triggerNames = "";
      for (Int_t i = 0; i < 64; ++i) {
        triggerNames += Form("%s", AliQnCorrectionsVarManagerTask::fOfflineTriggerNames[i]);
        triggerNames += ";";
      }

      histos->AddHistogram(classStr.Data(), "Triggers", "Offline triggers fired; ; ;", kFALSE, 64, -0.5, 63.5,
                 AliQnCorrectionsVarManagerTask::kOfflineTrigger, 2, -0.5, 1.5,
                 AliQnCorrectionsVarManagerTask::kOfflineTriggerFired, 0, 0.0, 0.0,
                 AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data(), "off;on");
      histos->AddHistogram(classStr.Data(), "Triggers2", "Offline triggers fired; ; ;", kFALSE, 64, -0.5, 63.5,
                 AliQnCorrectionsVarManagerTask::kOfflineTriggerFired2, 0, 0.0, 0.0,
                 AliQnCorrectionsVarManagerTask::kNothing, 0, 0.0, 0.0,
                 AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentVZERO_Triggers2",
                 "Offline triggers fired vs centrality VZERO; ; centrality VZERO;", kFALSE, 64, -0.5,
                 63.5, AliQnCorrectionsVarManagerTask::kOfflineTriggerFired2, 20, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentVZERO, 0, 0.0, 0.0,
                 AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentTPC_Triggers2",
                 "Offline triggers fired vs centrality TPC; ; centrality TPC;", kFALSE, 64, -0.5, 63.5,
                 AliQnCorrectionsVarManagerTask::kOfflineTriggerFired2, 20, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentTPC, 0, 0.0, 0.0,
                 AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentSPD_Triggers2",
                 "Offline triggers fired vs centrality SPD; ; centrality SPD;", kFALSE, 64, -0.5, 63.5,
                 AliQnCorrectionsVarManagerTask::kOfflineTriggerFired2, 20, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentSPD, 0, 0.0, 0.0,
                 AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentZDC_Triggers2",
                 "Offline triggers fired vs centrality ZDC; ; centrality ZDC;", kFALSE, 64, -0.5, 63.5,
                 AliQnCorrectionsVarManagerTask::kOfflineTriggerFired2, 20, 0.0, 100.0,
                 AliQnCorrectionsVarManagerTask::kCentZDC, 0, 0.0, 0.0,
                 AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "VtxZ_Triggers2", "Offline triggers fired vs vtxZ; ; vtx Z (cm.);",
                 kFALSE, 64, -0.5, 63.5, AliQnCorrectionsVarManagerTask::kOfflineTriggerFired2, 200,
                 -20.0, 20.0, AliQnCorrectionsVarManagerTask::kVtxZ, 0, 0.0, 0.0,
                 AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data());
      continue;
    }

    // Track histograms
    if (classStr.Contains("Tracks")) {
      histos->AddHistClass(classStr.Data());
      for (Int_t ih = 0; ih < 6; ++ih) {
        histos->AddHistogram(
         classStr.Data(), Form("Cos%dPhi_CentVZERO", ih + 1),
         Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih + 1, ih + 1), kTRUE, 20,
         0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 500, -1., 1.,
         AliQnCorrectionsVarManagerTask::kCosNPhi + ih);
        histos->AddHistogram(
         classStr.Data(), Form("Sin%dPhi_CentVZERO", ih + 1),
         Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih + 1, ih + 1), kTRUE, 20,
         0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 500, -1.0, 1.0,
         AliQnCorrectionsVarManagerTask::kSinNPhi + ih);
      }
    }

    // Track histograms
    if (classStr.Contains("TrackQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution; p_{T} (GeV/c^{2});", kFALSE, 1000, 0.0,
                 50.0, AliQnCorrectionsVarManagerTask::kPt);
      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE, 1000, -1.5, 1.5,
                 AliQnCorrectionsVarManagerTask::kEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE, 1000, 0.0, 6.3,
                 AliQnCorrectionsVarManagerTask::kPhi);
      histos->AddHistogram(classStr.Data(), "DCAxy", "DCAxy; DCAxy (cm.)", kFALSE, 1000, -10.0, 10.0,
                 AliQnCorrectionsVarManagerTask::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz", "DCAz; DCAz (cm.)", kFALSE, 1000, -10.0, 10.0,
                 AliQnCorrectionsVarManagerTask::kDcaZ);
      histos->AddHistogram(classStr.Data(), "TPCncls", "TPCncls; TPCncls", kFALSE, 160, 0.0, 160.0,
                 AliQnCorrectionsVarManagerTask::kTPCncls);
      histos->AddHistogram(classStr.Data(), "TPCsa_TPCncls", "TPC standalone TPCncls; TPCncls", kFALSE, 160, 0.0,
                 160.0, AliQnCorrectionsVarManagerTask::kTPCnclsIter1);

      // run dependence
      histos->AddHistogram(classStr.Data(), "Pt_Run", "<p_{T}> vs run; run;", kTRUE, kNRunBins, runHistRange[0],
                 runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 1000, 0.0, 50.0,
                 AliQnCorrectionsVarManagerTask::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Run", "<#eta> vs run; run;", kTRUE, kNRunBins, runHistRange[0],
                 runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 1000, -1.5, 1.5,
                 AliQnCorrectionsVarManagerTask::kEta);
      histos->AddHistogram(classStr.Data(), "Phi_Run", "<#varphi> vs run; run;", kTRUE, kNRunBins,
                 runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 1000, 0.0,
                 6.3, AliQnCorrectionsVarManagerTask::kPhi);
      histos->AddHistogram(classStr.Data(), "DCAxy_Run", "<DCAxy> vs run; run;", kTRUE, kNRunBins,
                 runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 1000, -10.0,
                 10.0, AliQnCorrectionsVarManagerTask::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Run", "<DCAz> vs run; run;", kTRUE, kNRunBins, runHistRange[0],
                 runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 1000, -10.0, 10.0,
                 AliQnCorrectionsVarManagerTask::kDcaZ);

      // correlations between parameters
      histos->AddHistogram(classStr.Data(), "Eta_Pt_prof", "<p_{T}> vs #eta; #eta; p_{T} (GeV/c);", kTRUE, 300,
                 -1.5, +1.5, AliQnCorrectionsVarManagerTask::kEta, 100, 0.0, 10.0,
                 AliQnCorrectionsVarManagerTask::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt", "p_{T} vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kFALSE,
                 300, -0.01, 6.3, AliQnCorrectionsVarManagerTask::kPhi, 100, 0.0, 2.2,
                 AliQnCorrectionsVarManagerTask::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt_prof", "<p_{T}> vs #varphi; #varphi (rad.); p_{T} (GeV/c)",
                 kTRUE, 300, 0.0, 6.3, AliQnCorrectionsVarManagerTask::kPhi, 100, 0.0, 10.0,
                 AliQnCorrectionsVarManagerTask::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE, 200,
                 -1.0, +1.0, AliQnCorrectionsVarManagerTask::kEta, 100, 0.0, 6.3,
                 AliQnCorrectionsVarManagerTask::kPhi);
      histos->AddHistogram(
       classStr.Data(), "TPCncls_Eta_Phi_prof", "<TPC ncls> vs #varphi vs #eta; #eta; #varphi (rad.);TPC ncls",
       kTRUE, 200, -1.0, +1.0, AliQnCorrectionsVarManagerTask::kEta, 100, 0.0, 6.3,
       AliQnCorrectionsVarManagerTask::kPhi, 10, 0.0, 200., AliQnCorrectionsVarManagerTask::kTPCncls);
      histos->AddHistogram(
       classStr.Data(), "DCAxy_Eta_Phi_prof", "<DCAxy> vs #varphi vs #eta; #eta; #varphi (rad.);DCAxy (cm)",
       kTRUE, 200, -1.0, +1.0, AliQnCorrectionsVarManagerTask::kEta, 100, 0.0, 6.3,
       AliQnCorrectionsVarManagerTask::kPhi, 10, 0.0, 200., AliQnCorrectionsVarManagerTask::kDcaXY);
      histos->AddHistogram(
       classStr.Data(), "DCAz_Eta_Phi_prof", "<DCAz> vs #varphi vs #eta; #eta; #varphi (rad.);DCAz (cm)", kTRUE,
       200, -1.0, +1.0, AliQnCorrectionsVarManagerTask::kEta, 100, 0.0, 6.3,
       AliQnCorrectionsVarManagerTask::kPhi, 10, 0.0, 200., AliQnCorrectionsVarManagerTask::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Pt_DCAxy", "DCAxy vs p_{T}; p_{T} (GeV/c); DCA_{xy} (cm)", kFALSE,
                 100, 0.0, 10.0, AliQnCorrectionsVarManagerTask::kPt, 500, -2.0, 2.0,
                 AliQnCorrectionsVarManagerTask::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Pt_DCAz", "DCAz vs p_{T}; p_{T} (GeV/c); DCA_{z} (cm)", kFALSE, 100,
                 0.0, 10.0, AliQnCorrectionsVarManagerTask::kPt, 500, -2.0, 2.0,
                 AliQnCorrectionsVarManagerTask::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Eta_DCAxy", "DCAxy vs #eta; #eta; DCA_{xy} (cm)", kFALSE, 100, -1.0,
                 1.0, AliQnCorrectionsVarManagerTask::kEta, 500, -2.0, 2.0,
                 AliQnCorrectionsVarManagerTask::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Eta_DCAz", "DCAz vs #eta; #eta; DCA_{z} (cm)", kFALSE, 100, -1.0,
                 1.0, AliQnCorrectionsVarManagerTask::kEta, 500, -2.0, 2.0,
                 AliQnCorrectionsVarManagerTask::kDcaZ);

      for (Int_t ih = 0; ih < 6; ++ih) {
        // histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs
        // (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 500, -1., 1.,
        //                   AliQnCorrectionsVarManagerTask::kCosNPhi+ih);
        // histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs
        // (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 500, -1.0, 1.0,
        //                   AliQnCorrectionsVarManagerTask::kSinNPhi+ih);
        histos->AddHistogram(
         classStr.Data(), Form("Cos%dPhi_Pt_Eta", ih + 1),
         Form("<cos%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <cos%d #varphi>", ih + 1, ih + 1), kTRUE,
         20, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kEta, 30, 0.0, 3.0,
         AliQnCorrectionsVarManagerTask::kPt, 500, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kCosNPhi + ih);
        histos->AddHistogram(
         classStr.Data(), Form("Sin%dPhi_Pt_Eta", ih + 1),
         Form("<sin%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <sin%d #varphi>", ih + 1, ih + 1), kTRUE,
         20, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kEta, 30, 0.0, 3.0,
         AliQnCorrectionsVarManagerTask::kPt, 500, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kSinNPhi + ih);
        histos->AddHistogram(
         classStr.Data(), Form("Cos%dPhi_CentVZERO_VtxZ", ih + 1),
         Form("<cos%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <cos%d #varphi>", ih + 1,
            ih + 1),
         kTRUE, 30, -15.0, 15.0, AliQnCorrectionsVarManagerTask::kVtxZ, 20, 0.0, 100.0,
         AliQnCorrectionsVarManagerTask::kCentVZERO, 500, -1., 1.,
         AliQnCorrectionsVarManagerTask::kCosNPhi + ih);
        histos->AddHistogram(
         classStr.Data(), Form("Sin%dPhi_CentVZERO_VtxZ", ih + 1),
         Form("<sin%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <sin%d #varphi>", ih + 1,
            ih + 1),
         kTRUE, 30, -15.0, 15.0, AliQnCorrectionsVarManagerTask::kVtxZ, 20, 0.0, 100.0,
         AliQnCorrectionsVarManagerTask::kCentVZERO, 500, -1.0, 1.0,
         AliQnCorrectionsVarManagerTask::kSinNPhi + ih);
      }
    }

    // Tracklet histograms
    if (classStr.Contains("TrackletQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE, 1000, -3.0, 3.0,
                 AliQnCorrectionsVarManagerTask::kSPDtrackletEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE, 300, -0.01, 6.3,
                 AliQnCorrectionsVarManagerTask::kSPDtrackletPhi);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE, 200,
                 -3.0, +3.0, AliQnCorrectionsVarManagerTask::kSPDtrackletEta, 100, 0.0, 6.3,
                 AliQnCorrectionsVarManagerTask::kSPDtrackletPhi);
    }
  }

  AliInfoGeneralStream("PWGJE::EMCALJetTasks::FlowVectorCorrections::DefineHistograms") << " done\n";
}

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */
