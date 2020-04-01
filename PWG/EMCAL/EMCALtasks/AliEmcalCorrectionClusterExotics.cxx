// AliEmcalCorrectionClusterExotics
//

#include "AliEmcalCorrectionClusterExotics.h"

#include <TH2F.h>
#include <TList.h>

#include "AliClusterContainer.h"
#include "AliEMCALRecoUtils.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionClusterExotics);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionClusterExotics> AliEmcalCorrectionClusterExotics::reg("AliEmcalCorrectionClusterExotics");

/**
 * Default constructor
 */
AliEmcalCorrectionClusterExotics::AliEmcalCorrectionClusterExotics() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionClusterExotics"),
  fEtaPhiDistBefore(0),
  fEtaPhiDistAfter(0),
  fEtaPhiDistAfterNDiffCut(0),
  fEnergyExoticClusters(0),
  fEnergyExoticClustersNDiffCut(0),
  fExoticMinCellAmplitude(0),
  fMaxFcross(0),
  fCellCrossMaxTimeDiff(0),
  fHighEnergyNdiffCut(0),
  fMinCellEnNdiffCut(0)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionClusterExotics::~AliEmcalCorrectionClusterExotics()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionClusterExotics::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  
  fExoticMinCellAmplitude = 4.;
  GetProperty("fExoticMinCellAmplitude", fExoticMinCellAmplitude);
  fMaxFcross = 0.97;
  GetProperty("fMaxFcross", fMaxFcross);
  fCellCrossMaxTimeDiff = 1e6;
  GetProperty("fCellCrossMaxTimeDiff", fCellCrossMaxTimeDiff);
  fHighEnergyNdiffCut = 300; // open, use ~50
  GetProperty("fHighEnergyNdiffCut", fHighEnergyNdiffCut);
  fMinCellEnNdiffCut = 0;
  GetProperty("fMinCellEnNdiffCut",fMinCellEnNdiffCut);

  // init reco utils
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;
  fRecoUtils->SwitchOnRejectExoticCluster();
  fRecoUtils->SetExoticCellMinAmplitudeCut(fExoticMinCellAmplitude);
  fRecoUtils->SetExoticCellFractionCut(fMaxFcross);
  fRecoUtils->SetExoticCellDiffTimeCut(fCellCrossMaxTimeDiff);
  if (fRecoUtils)
    fRecoUtils->Print("");

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionClusterExotics::UserCreateOutputObjects()
{   
  AliEmcalCorrectionComponent::UserCreateOutputObjects();
  
  // Create my user objects.
  if (fCreateHisto){
    fEtaPhiDistBefore = new TH2F("hEtaPhiDistBefore","hEtaPhiDistBefore;#eta;#phi",280,-0.7,0.7,200*3.14,0,2*3.14);
    fOutput->Add(fEtaPhiDistBefore);
    fEtaPhiDistAfter = new TH2F("hEtaPhiDistAfter","hEtaPhiDistAfter;#eta;#phi",280,-0.7,0.7,200*3.14,0,2*3.14);
    fOutput->Add(fEtaPhiDistAfter);
    fEtaPhiDistAfterNDiffCut = new TH2F("hEtaPhiDistAfterNDiffCut","hEtaPhiDistAfterNDiffCut;#eta;#phi",280,-0.7,0.7,200*3.14,0,2*3.14);
    fOutput->Add(fEtaPhiDistAfterNDiffCut);
    fEnergyExoticClusters = new TH1F("fEnergyExoticClusters","fEnergyExoticClusters;E_{ex clus} (GeV)",2000,0,200);
    fOutput->Add(fEnergyExoticClusters);
    fEnergyExoticClustersNDiffCut = new TH1F("fEnergyExoticClustersNDiffCut","fEnergyExoticClustersNDiffCut;E_{ex clus} (GeV)",2000,0,200);
    fOutput->Add(fEnergyExoticClustersNDiffCut);
    
    // Take ownership of output list
    fOutput->SetOwner(kTRUE);
  }
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionClusterExotics::Run()
{
  AliEmcalCorrectionComponent::Run();
  
  // loop over clusters
  AliVCluster *clus = 0;
  AliClusterContainer * clusCont = 0;
  TIter nextClusCont(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(nextClusCont()))) {

    if (!clusCont) continue;
    auto clusItCont = clusCont->all_momentum();

    for (AliClusterIterableMomentumContainer::iterator clusIterator = clusItCont.begin(); clusIterator != clusItCont.end(); ++clusIterator) {
      clus = static_cast<AliVCluster *>(clusIterator->second);

      if (!clus->IsEMCAL()) continue;

      if (fCreateHisto) {
        Float_t pos[3] = {0.};
        clus->GetPosition(pos);
        TVector3 vec(pos);
        // Phi needs to be in 0 to 2 Pi
        fEtaPhiDistBefore->Fill(vec.Eta(), TVector2::Phi_0_2pi(vec.Phi()));
      }

      Bool_t exResult = kFALSE;
      Bool_t exResult2= kFALSE;

      if (fRecoUtils) {
        if (fRecoUtils->IsRejectExoticCluster()) {
          Bool_t exRemoval = fRecoUtils->IsRejectExoticCell();
          fRecoUtils->SwitchOnRejectExoticCell();                  //switch on temporarily
          
          //exResult = fRecoUtils->IsExoticCluster(clus, fCaloCells);
         
          // Copy fRecoUtils->IsExoticCluster(), to avoid double calling of absIdMax finding
          //
          AliEMCALGeometry* geom = AliEMCALGeometry::GetInstance();
          
          if(!geom)
          {
            AliError("No instance of the geometry is available");
            return kFALSE;
          }
          
          Int_t iSupMod = -1, absIdMax = -1, ieta = -1, iphi = -1;
          Bool_t shared = kFALSE;
          fRecoUtils->GetMaxEnergyCell(geom, fCaloCells, clus, 
                                       absIdMax, iSupMod, ieta, iphi, shared);
          
          exResult = fRecoUtils->IsExoticCell(absIdMax, fCaloCells);
          
          if (!exRemoval) fRecoUtils->SwitchOffRejectExoticCell(); //switch back off

          clus->SetIsExotic(exResult);
          if ( !exResult && clus->E() >  fHighEnergyNdiffCut )
          {
            Int_t   nDiff = 0, nSame = 0; 
            Float_t eDiff = 0, eSame = 0;
            fRecoUtils->GetEnergyAndNumberOfCellsInTCard(clus, absIdMax, fCaloCells, 
                                                         nDiff, nSame, eDiff, eSame, 
                                                         fMinCellEnNdiffCut) ;
            
            if ( nDiff == 0 ) exResult2 = kTRUE;
            
            clus->SetIsExotic(exResult2);
          }
        }
      }

      if (fCreateHisto) {
        if (exResult) {
          fEnergyExoticClusters->Fill(clus->E());
        }
        else {
          Float_t pos[3] = {0.};
          clus->GetPosition(pos);
          TVector3 vec(pos);
          // Phi needs to be in 0 to 2 Pi
          fEtaPhiDistAfter->Fill(vec.Eta(), TVector2::Phi_0_2pi(vec.Phi()));
          if(!exResult2) fEtaPhiDistAfterNDiffCut->Fill(vec.Eta(), TVector2::Phi_0_2pi(vec.Phi()));
        }
        if(exResult2) fEnergyExoticClustersNDiffCut->Fill(clus->E());
      }
    }
  }
  
  return kTRUE;
}
