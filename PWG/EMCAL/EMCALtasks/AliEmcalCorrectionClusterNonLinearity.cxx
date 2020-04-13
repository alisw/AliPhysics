// AliEmcalCorrectionClusterNonLinearity
//

#include "AliEmcalCorrectionClusterNonLinearity.h"

#include <TList.h>

#include "AliClusterContainer.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionClusterNonLinearity);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionClusterNonLinearity> AliEmcalCorrectionClusterNonLinearity::reg("AliEmcalCorrectionClusterNonLinearity");

const std::map <std::string, AliEMCALRecoUtils::NonlinearityFunctions> AliEmcalCorrectionClusterNonLinearity::fgkNonlinearityFunctionMap = {
    { "kPi0MC", AliEMCALRecoUtils::kPi0MC },
    { "kPi0GammaGamma", AliEMCALRecoUtils::kPi0GammaGamma },
    { "kPi0GammaConversion", AliEMCALRecoUtils::kPi0GammaConversion },
    { "kNoCorrection", AliEMCALRecoUtils::kNoCorrection },
    { "kBeamTest", AliEMCALRecoUtils::kBeamTest },
    { "kBeamTestCorrected", AliEMCALRecoUtils::kBeamTestCorrected },
    { "kPi0MCv2", AliEMCALRecoUtils::kPi0MCv2 },
    { "kPi0MCv3", AliEMCALRecoUtils::kPi0MCv3 },
    { "kBeamTestCorrectedv2", AliEMCALRecoUtils::kBeamTestCorrectedv2 },
    { "kSDMv5", AliEMCALRecoUtils::kSDMv5 },
    { "kPi0MCv5", AliEMCALRecoUtils::kPi0MCv5 },
    { "kSDMv6", AliEMCALRecoUtils::kSDMv6 },
    { "kPi0MCv6", AliEMCALRecoUtils::kPi0MCv6 },
    { "kBeamTestCorrectedv3", AliEMCALRecoUtils::kBeamTestCorrectedv3 },
    { "kPCMv1", AliEMCALRecoUtils::kPCMv1 },
    { "kPCMplusBTCv1", AliEMCALRecoUtils::kPCMplusBTCv1 },
    { "kPCMsysv1", AliEMCALRecoUtils::kPCMsysv1 },
    { "kBeamTestCorrectedv4", AliEMCALRecoUtils::kBeamTestCorrectedv4 },
    { "kBeamTestNS", AliEMCALRecoUtils::kBeamTestNS },
    { "kPi0MCNS", AliEMCALRecoUtils::kPi0MCNS }
};

/**
 * Default constructor
 */
AliEmcalCorrectionClusterNonLinearity::AliEmcalCorrectionClusterNonLinearity() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionClusterNonLinearity"),
  fEnergyDistBefore(0),
  fEnergyTimeHistBefore(0),
  fEnergyDistAfter(0),
  fEnergyTimeHistAfter(0),
  fSetForceClusterE(kFALSE)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionClusterNonLinearity::~AliEmcalCorrectionClusterNonLinearity()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionClusterNonLinearity::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  
  std::string nonLinFunctStr = "";
  GetProperty("nonLinFunct", nonLinFunctStr);
  UInt_t nonLinFunct = fgkNonlinearityFunctionMap.at(nonLinFunctStr);

  GetProperty("setForceClusterE", fSetForceClusterE);
  
  // init reco utils
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;
  fRecoUtils->SetNonLinearityFunction(nonLinFunct);
  
  if (fRecoUtils) {
    fRecoUtils->InitNonLinearityParam();
    fRecoUtils->Print("");
  }

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionClusterNonLinearity::UserCreateOutputObjects()
{   
  AliEmcalCorrectionComponent::UserCreateOutputObjects();
  
  // Create my user objects.
  if (fCreateHisto){
    fEnergyDistBefore = new TH1F("hEnergyDistBefore","hEnergyDistBefore;E_{clus} (GeV)",1500,0,150);
    fOutput->Add(fEnergyDistBefore);
    fEnergyTimeHistBefore = new TH2F("hEnergyTimeDistBefore","hEnergyTimeDistBefore;E_{clus} (GeV);time (s)",1500,0,150,500,-1e-6,1e-6);
    fOutput->Add(fEnergyTimeHistBefore);
    fEnergyDistAfter = new TH1F("hEnergyDistAfter","hEnergyDistAfter;E_{clus} (GeV)",1500,0,150);
    fOutput->Add(fEnergyDistAfter);
    fEnergyTimeHistAfter = new TH2F("hEnergyTimeDistAfter","hEnergyTimeDistAfter;E_{clus} (GeV);time (s)",1500,0,150,500,-1e-6,1e-6);
    fOutput->Add(fEnergyTimeHistAfter);
    
    // Take ownership of output list
    fOutput->SetOwner(kTRUE);
  }
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionClusterNonLinearity::Run()
{
  AliEmcalCorrectionComponent::Run();
  
  // loop over clusters
  AliVCluster *clus = 0;
  AliClusterContainer * clusCont = 0;
  TIter nextClusCont(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(nextClusCont()))) {

    if (!clusCont) return kFALSE;
    auto clusItCont = clusCont->all_momentum();

    for (AliClusterIterableMomentumContainer::iterator clusIterator = clusItCont.begin(); clusIterator != clusItCont.end(); ++clusIterator) {
      clus = static_cast<AliVCluster *>(clusIterator->second);

      if (!clus->IsEMCAL()) continue;

      if (fCreateHisto) {
        fEnergyDistBefore->Fill(clus->E());
        fEnergyTimeHistBefore->Fill(clus->E(), clus->GetTOF());
      }

      if (fRecoUtils) {
        if (fRecoUtils->GetNonLinearityFunction() != AliEMCALRecoUtils::kNoCorrection) {
          Double_t energy = fRecoUtils->CorrectClusterEnergyLinearity(clus);
          clus->SetNonLinCorrEnergy(energy);
          if ( fSetForceClusterE ) clus->SetE(energy);
        }
      }

      // Fill histograms only if cluster is not exotic, as in ClusterMaker (the clusters are flagged, not removed)
      if (fCreateHisto && !clus->GetIsExotic()) {
        Float_t energy = clus->GetNonLinCorrEnergy();
        if(fSetForceClusterE) energy = clus->E();
        
        fEnergyDistAfter->Fill(energy);
        fEnergyTimeHistAfter->Fill(energy, clus->GetTOF());
      }
    }
  }
  
  return kTRUE;
}
