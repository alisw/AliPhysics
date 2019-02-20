// AliEmcalCorrectionClusterEnergyVariation
//

#include "AliEmcalCorrectionClusterEnergyVariation.h"

#include "AliClusterContainer.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionClusterEnergyVariation);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionClusterEnergyVariation> AliEmcalCorrectionClusterEnergyVariation::reg("AliEmcalCorrectionClusterEnergyVariation");

/**
 * Default constructor
 */
AliEmcalCorrectionClusterEnergyVariation::AliEmcalCorrectionClusterEnergyVariation() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionClusterEnergyVariation"),
  fEnergyScaleShift(0.),
  fSmearingWidth(0.),
  fRandom(0)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionClusterEnergyVariation::~AliEmcalCorrectionClusterEnergyVariation()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionClusterEnergyVariation::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  
  GetProperty("energyScaleShift", fEnergyScaleShift);
  GetProperty("smearingWidth", fSmearingWidth);

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionClusterEnergyVariation::UserCreateOutputObjects()
{   
  AliEmcalCorrectionComponent::UserCreateOutputObjects();
}

/**
 * Called before the first event to initialize the correction.
 */
void AliEmcalCorrectionClusterEnergyVariation::ExecOnce()
{
  fRandom.SetSeed(0);
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionClusterEnergyVariation::Run()
{
  AliEmcalCorrectionComponent::Run();
  
  // Loop over all clusters
  AliTLorentzVector clusFourVec;
  AliVCluster *cluster = 0;
  AliClusterContainer * clusCont = 0;
  TIter nextClusCont(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(nextClusCont()))) {
    auto clusItCont = clusCont->all_momentum();
    for (AliClusterIterableMomentumContainer::iterator clusIterator = clusItCont.begin(); clusIterator != clusItCont.end(); ++clusIterator) {
      
      // Get the cluster energy (using the default energy type set in the cluster container)
      clusFourVec.Clear();
      clusFourVec = clusIterator->first;
      Double_t energyclus = clusFourVec.E();

      // Scale energy (if applicable)
      if (TMath::Abs(fEnergyScaleShift) > 1e-6) {
        energyclus *= (1 + fEnergyScaleShift);
      }
      
      // Smear energy (if applicable)
      if (fSmearingWidth > 1e-6) {
        energyclus *= (1 + fRandom.Gaus(0, fSmearingWidth) );
      }
      
      // Set the shifted/smeared energy to the user-defined energy field
      cluster = static_cast<AliVCluster *>(clusIterator->second);
      if (cluster->IsEMCAL() && energyclus > 0) {
        cluster->SetUserDefEnergy(AliVCluster::kUserDefEnergy1, energyclus);
      }
      
    }
  }

  return kTRUE;
}
