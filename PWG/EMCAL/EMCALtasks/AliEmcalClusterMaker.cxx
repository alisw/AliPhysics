//
// Cluster maker task.
//
// Author: C.Loizides, S.Aiola

#include <TChain.h>
#include <TClonesArray.h>
#include "AliAODCaloCluster.h"
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliEmcalClusterMaker.h"
#include "AliClusterContainer.h"

ClassImp(AliEmcalClusterMaker)

//________________________________________________________________________
AliEmcalClusterMaker::AliEmcalClusterMaker() : 
  AliAnalysisTaskEmcal("AliEmcalClusterMaker", kFALSE),
  fOutCaloName(),
  fRecoUtils(0),
  fEsdMode(kTRUE),
  fOutClusters(0),
  fEnergyDistBefore(0),
  fEtaPhiDistBefore(0),
  fEnergyTimeHistBefore(0),
  fEnergyDistAfter(0),
  fEtaPhiDistAfter(0),
  fEnergyTimeHistAfter(0),
  fEnergyExoticClusters(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliEmcalClusterMaker::AliEmcalClusterMaker(const char *name, Bool_t histo) : 
   AliAnalysisTaskEmcal(name, histo),
   fOutCaloName("EmcClusters"),
   fRecoUtils(0),
   fEsdMode(kTRUE),
   fOutClusters(0),
   fEnergyDistBefore(0),
   fEtaPhiDistBefore(0),
   fEnergyTimeHistBefore(0),
   fEnergyDistAfter(0),
   fEtaPhiDistAfter(0),
   fEnergyTimeHistAfter(0),
   fEnergyExoticClusters(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(histo);

  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";
}

//________________________________________________________________________
AliEmcalClusterMaker::~AliEmcalClusterMaker()
{
  // Destructor
}

//________________________________________________________________________
void AliEmcalClusterMaker::UserCreateOutputObjects()
{
  // Create my user objects.

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if (fRecoUtils) {
    fRecoUtils->InitNonLinearityParam();
    fRecoUtils->Print("");
  }

  if (!fCreateHisto) return;

  fEnergyDistBefore = new TH1F("hEnergyDistBefore","hEnergyDistBefore;E_{clus} (GeV)",1500,0,150);
  fOutput->Add(fEnergyDistBefore);
  fEtaPhiDistBefore = new TH2F("hEtaPhiDistBefore","hEtaPhiDistBefore;#eta;#phi",280,-0.7,0.7,800,1.3,3.3);
  fOutput->Add(fEtaPhiDistBefore);
  fEnergyTimeHistBefore = new TH2F("hEnergyTimeDistBefore","hEnergyTimeDistBefore;E_{clus} (GeV);time",1500,0,150,500,0,1e-6);
  fOutput->Add(fEnergyTimeHistBefore);
  fEnergyDistAfter = new TH1F("hEnergyDistAfter","hEnergyDistAfter;E_{clus} (GeV)",1500,0,150);
  fOutput->Add(fEnergyDistAfter);
  fEtaPhiDistAfter = new TH2F("hEtaPhiDistAfter","hEtaPhiDistAfter;#eta;#phi",280,-0.7,0.7,800,1.3,3.3);
  fOutput->Add(fEtaPhiDistAfter);
  fEnergyTimeHistAfter = new TH2F("hEnergyTimeDistAfter","hEnergyTimeDistAfter;E_{clus} (GeV);time",1500,0,150,500,0,1e-6);
  fOutput->Add(fEnergyTimeHistAfter);
  fEnergyExoticClusters = new TH1F("fEnergyExoticClusters","fEnergyExoticClusters;E_{ex clus} (GeV)",1500,0,150);
  fOutput->Add(fEnergyExoticClusters);
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliEmcalClusterMaker::ExecOnce() 
{
  // Initialize the analysis.

  // Do base class initializations and if it fails -> bail out
  AliAnalysisTaskEmcal::ExecOnce();
  if (!fInitialized) {
    return;
  }

  if (dynamic_cast<AliAODEvent*>(InputEvent())) {
    fEsdMode = kFALSE;
  }

  if (!fOutCaloName.IsNull()) {  // if empty updates old clusters instead of creating a new collection
    if (fEsdMode) {
      fOutClusters = new TClonesArray("AliESDCaloCluster");
    }
    else {
      fOutClusters = new TClonesArray("AliAODCaloCluster");
    }
    fOutClusters->SetName(fOutCaloName);

    // post output in event if not yet present
    if (!(InputEvent()->FindListObject(fOutCaloName))) {
      InputEvent()->AddObject(fOutClusters);
    }
    else {
      fInitialized = kFALSE;
      AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fOutCaloName.Data()));
      return;
    }
  }
}

//________________________________________________________________________
Bool_t AliEmcalClusterMaker::Run() 
{
  // Run the cluster maker

  // delete output
  if (fOutClusters) fOutClusters->Delete();

  AliClusterContainer* clusters = GetClusterContainer(0);
  if (!clusters) return kFALSE;

  // loop over clusters
  Int_t clusCount = 0;
  //Int_t entries   = fCaloClusters->GetEntries();
  AliVCluster *clus = 0;
  clusters->ResetCurrentID();
  while ((clus = clusters->GetNextCluster())) {
    if (!clus->IsEMCAL()) continue;

    if (fCreateHisto) {
      fEnergyDistBefore->Fill(clus->E());
      Float_t pos[3] = {0.};
      clus->GetPosition(pos);
      TVector3 vec(pos);
      fEtaPhiDistBefore->Fill(vec.Eta(),vec.Phi());
      fEnergyTimeHistBefore->Fill(clus->E(), clus->GetTOF());
    }

    Bool_t exResult = kFALSE;

    if (fRecoUtils) {
      if (fRecoUtils->IsRejectExoticCluster()) {
        Bool_t exRemoval = fRecoUtils->IsRejectExoticCell();
        fRecoUtils->SwitchOnRejectExoticCell();                  //switch on temporarily
        exResult = fRecoUtils->IsExoticCluster(clus, fCaloCells);
        if (!exRemoval) fRecoUtils->SwitchOffRejectExoticCell(); //switch back off

        clus->SetIsExotic(exResult);
      }
      if (fRecoUtils->GetNonLinearityFunction() != AliEMCALRecoUtils::kNoCorrection) {
        Double_t energy = fRecoUtils->CorrectClusterEnergyLinearity(clus);
        clus->SetNonLinCorrEnergy(energy);
      }
    }

    if (fCreateHisto) {
      if (exResult) {
        fEnergyExoticClusters->Fill(clus->E());
      }
      else {
        fEnergyDistAfter->Fill(clus->GetNonLinCorrEnergy());
        Float_t pos[3] = {0.};
        clus->GetPosition(pos);
        TVector3 vec(pos);
        fEtaPhiDistAfter->Fill(vec.Eta(), vec.Phi());
        fEnergyTimeHistAfter->Fill(clus->GetNonLinCorrEnergy(), clus->GetTOF());
      }
    }

    UInt_t rejectionReason = 0;
    if (fOutClusters && clusters->AcceptCluster(clus, rejectionReason) && !exResult) {

      AliVCluster *oc = 0;
      if (fEsdMode) {
        AliESDCaloCluster *ec = dynamic_cast<AliESDCaloCluster*>(clus);
        if (!ec) continue;
        oc = new ((*fOutClusters)[clusCount]) AliESDCaloCluster(*ec);
      }
      else {
        AliAODCaloCluster *ac = dynamic_cast<AliAODCaloCluster*>(clus);
        if (!ac) continue;
        oc = new ((*fOutClusters)[clusCount]) AliAODCaloCluster(*ac);
      }

      oc->SetE(clus->GetNonLinCorrEnergy());
      oc->SetNonLinCorrEnergy(clus->GetNonLinCorrEnergy());
      oc->SetHadCorrEnergy(0);//just to make sure that in this container version there is no hadronic correction applied - can be set at a later stage
      clusCount++;
    }
  }

  return kTRUE;
}
