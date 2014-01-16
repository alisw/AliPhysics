// $Id$
//
// Cluster maker task.
//
// Author: C.Loizides

#include <TChain.h>
#include <TClonesArray.h>
#include "AliAODCaloCluster.h"
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliEmcalClusterMaker.h"

ClassImp(AliEmcalClusterMaker)

//________________________________________________________________________
AliEmcalClusterMaker::AliEmcalClusterMaker() : 
  AliAnalysisTaskEmcal("AliEmcalClusterMaker", kFALSE),
  fOutCaloName(),
  fRecoUtils(0),
  fEsdMode(kTRUE),
  fOutClusters(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliEmcalClusterMaker::AliEmcalClusterMaker(const char *name, Bool_t histo) : 
  AliAnalysisTaskEmcal(name, histo),
  fOutCaloName("EmcClusters"),
  fRecoUtils(0),
  fEsdMode(kTRUE),
  fOutClusters(0)
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

  if (fRecoUtils)
    fRecoUtils->Print("");
  //  PostData(1, fOutput);
}

//________________________________________________________________________
void AliEmcalClusterMaker::ExecOnce() 
{
  // Initialize the analysis.

  // Do base class initializations and if it fails -> bail out
  AliAnalysisTaskEmcal::ExecOnce();
  if (!fInitialized) 
    return;

  if (dynamic_cast<AliAODEvent*>(InputEvent())) 
    fEsdMode = kFALSE;

  if (fEsdMode) 
    fOutClusters = new TClonesArray("AliESDCaloCluster");
  else 
    fOutClusters = new TClonesArray("AliAODCaloCluster");

  fOutClusters->SetName(fOutCaloName);

  // post output in event if not yet present
  if (!(InputEvent()->FindListObject(fOutCaloName))) {
    InputEvent()->AddObject(fOutClusters);
  } else {
    fInitialized = kFALSE;
    AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fOutCaloName.Data()));
    return;
  }
}

//________________________________________________________________________
Bool_t AliEmcalClusterMaker::Run() 
{
  // Run the hadronic correction

  // delete output
  fOutClusters->Delete();

  // loop over clusters
  Int_t clusCount = 0;
  Int_t entries   = fCaloClusters->GetEntries();
  for (Int_t i=0; i<entries; ++i) {
    AliVCluster *clus = static_cast<AliVCluster*>(fCaloClusters->At(i));
    if (!clus || !clus->IsEMCAL())
      continue;
    AliVCluster *oc = 0;
    if (fEsdMode) {
      AliESDCaloCluster *ec = dynamic_cast<AliESDCaloCluster*>(clus);
      if (!ec) continue;
      oc = new ((*fOutClusters)[clusCount]) AliESDCaloCluster(*ec);
    } else { 
      AliAODCaloCluster *ac = dynamic_cast<AliAODCaloCluster*>(clus);
      if (!ac) continue;
      oc = new ((*fOutClusters)[clusCount]) AliAODCaloCluster(*ac);
    }
    if (fRecoUtils) {
      if (fRecoUtils->IsRejectExoticCluster()) {
	if (fRecoUtils->IsExoticCluster(oc,fCaloCells))
	continue;
      }
      if (fRecoUtils->GetNonLinearityFunction()!=AliEMCALRecoUtils::kNoCorrection) {
	Double_t energy = fRecoUtils->CorrectClusterEnergyLinearity(oc);
	oc->SetE(energy);
      }
    }
    if (!AcceptCluster(oc))
      continue;
    clusCount++;
  }
  if ((clusCount>0) && (clusCount==fOutClusters->GetEntries()))
    fOutClusters->RemoveAt(clusCount);
  return kTRUE;
}
