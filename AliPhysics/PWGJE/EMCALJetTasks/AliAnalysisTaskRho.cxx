#include "AliAnalysisTaskRho.h"

#include <TClonesArray.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliVEventHandler.h"
#include "AliAnalysisDataContainer.h"

ClassImp(AliAnalysisTaskRho)

//________________________________________________________________________
AliAnalysisTaskRho::AliAnalysisTaskRho() : 
  AliAnalysisTaskRhoBase("AliAnalysisTaskRho"),
  fNExclLeadJets(0)
{
  // Constructor.
}

//________________________________________________________________________
AliAnalysisTaskRho::AliAnalysisTaskRho(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoBase(name, histo),
  fNExclLeadJets(0)
{
  // Constructor.
}


//________________________________________________________________________
Bool_t AliAnalysisTaskRho::Run() 
{
  // Run the analysis.

  fOutRho->SetVal(0);
  if (fOutRhoScaled)
    fOutRhoScaled->SetVal(0);

  if (!fJets)
    return kFALSE;

  const Int_t Njets   = fJets->GetEntries();

  Int_t maxJetIds[]   = {-1, -1};
  Float_t maxJetPts[] = { 0,  0};

  if (fNExclLeadJets > 0) {
    for (Int_t ij = 0; ij < Njets; ++ij) {
      AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(ij));
      if (!jet) {
	AliError(Form("%s: Could not receive jet %d", GetName(), ij));
	continue;
      } 

      if (!AcceptJet(jet))
        continue;

      if (jet->Pt() > maxJetPts[0]) {
	maxJetPts[1] = maxJetPts[0];
	maxJetIds[1] = maxJetIds[0];
	maxJetPts[0] = jet->Pt();
	maxJetIds[0] = ij;
      } else if (jet->Pt() > maxJetPts[1]) {
	maxJetPts[1] = jet->Pt();
	maxJetIds[1] = ij;
      }
    }
    if (fNExclLeadJets < 2) {
      maxJetIds[1] = -1;
      maxJetPts[1] = 0;
    }
  }

  static Double_t rhovec[999];
  Int_t NjetAcc = 0;

  // push all jets within selected acceptance into stack
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {

    // exlcuding lead jets
    if (iJets == maxJetIds[0] || iJets == maxJetIds[1])
      continue;

    AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(iJets));
    if (!jet) {
      AliError(Form("%s: Could not receive jet %d", GetName(), iJets));
      continue;
    } 

    if (!AcceptJet(jet))
      continue;

    rhovec[NjetAcc] = jet->Pt() / jet->Area();
    ++NjetAcc;
  }


  if (NjetAcc > 0) {
    //find median value
    Double_t rho = TMath::Median(NjetAcc, rhovec);
    fOutRho->SetVal(rho);

    if (fOutRhoScaled) {
      Double_t rhoScaled = rho * GetScaleFactor(fCent);
      fOutRhoScaled->SetVal(rhoScaled);
    }
  }

  return kTRUE;
} 

//________________________________________________________________________
AliAnalysisTaskRho* AliAnalysisTaskRho::AddTaskRhoNew (
    const char* nTracks, const char* nClusters, const char* nRho,
    Double_t jetradius, UInt_t acceptance,  AliJetContainer::EJetType_t jetType, const Bool_t histo,
    AliJetContainer::ERecoScheme_t rscheme, const char* suffix
)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { ::Error("AddTaskRho", "No analysis manager to connect to."); return NULL; }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) { ::Error("AddTaskEmcalJetSpectraQA", "This task requires an input event handler"); return NULL; }

  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(nTracks);
  TString clusName(nClusters);

  if (trackName == "usedefault") {
    if (dataType == kESD) {
      trackName = "Tracks";
    }
    else if (dataType == kAOD) {
      trackName = "tracks";
    }
    else {
      trackName = "";
    }
  }

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    }
    else if (dataType == kAOD) {
      clusName = "caloClusters";
    }
    else {
      clusName = "";
    }
  }

  TString name("AliAnalysisTaskRho");
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskRho* mgrTask = (AliAnalysisTaskRho*)(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;

  AliAnalysisTaskRho *rhotask = new AliAnalysisTaskRho(name, histo);
  rhotask->SetOutRhoName(nRho);

  if (trackName == "mcparticles") {
    AliMCParticleContainer* mcpartCont = rhotask->AddMCParticleContainer(trackName);
    mcpartCont->SelectPhysicalPrimaries(kTRUE);
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    AliTrackContainer* trackCont = rhotask->AddTrackContainer(trackName);
  }
  else if (!trackName.IsNull()) {
    rhotask->AddParticleContainer(trackName);
  }
  AliParticleContainer* partCont = rhotask->GetParticleContainer(0);

  AliClusterContainer *clusterCont = rhotask->AddClusterContainer(clusName);
  if (clusterCont) {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(0.3);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  AliJetContainer *jetCont = rhotask->AddJetContainer(jetType, AliJetContainer::kt_algorithm, rscheme, jetradius, acceptance, partCont, clusterCont);
  if (jetCont) jetCont->SetJetPtCut(0);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(rhotask);

  // Create containers for input/output
  mgr->ConnectInput(rhotask, 0, mgr->GetCommonInputContainer());
  if (histo) {
    TString contname(name);
    contname += "_histos";
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
        TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(rhotask, 1, coutput1);
  }

  return rhotask;
}

