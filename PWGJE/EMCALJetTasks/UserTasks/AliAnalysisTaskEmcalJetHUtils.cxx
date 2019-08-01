//
// Utilities class for Jet-Hadron correlation analysis
//

#include "AliAnalysisTaskEmcalJetHUtils.h"

// Require to use AliLog streams with some types
#include <iostream>
#include <cmath>

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

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */
