// $Id$
//
// Calculation of rho from a collection of jets.
// If scale function is given the scaled rho will be exported
// with the name as "fOutRhoName".apppend("_Scaled").
//
// Authors: R.Reed, S.Aiola, M.Connors

#include "AliAnalysisTaskRhoSparse.h"

#include <TClonesArray.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include <AliVEventHandler.h>
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"

ClassImp(AliAnalysisTaskRhoSparse)

AliAnalysisTaskRhoSparse::AliAnalysisTaskRhoSparse() : 
  AliAnalysisTaskRhoBase("AliAnalysisTaskRhoSparse"),
  fNExclLeadJets(0),
  fExcludeOverlaps(0),
  fRhoCMS(0),
  fUseTPCArea(0),
  fExcludeAreaExcludedJets(0),
  fHistOccCorrvsCent(0)
{
}

AliAnalysisTaskRhoSparse::AliAnalysisTaskRhoSparse(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoBase(name, histo),
  fNExclLeadJets(0),
  fExcludeOverlaps(0),
  fRhoCMS(0),
  fUseTPCArea(0),
  fExcludeAreaExcludedJets(0),
  fHistOccCorrvsCent(0)
{
}

void AliAnalysisTaskRhoSparse::UserCreateOutputObjects()
{
  if (!fCreateHisto) return;

  AliAnalysisTaskRhoBase::UserCreateOutputObjects();

  fHistOccCorrvsCent = new TH2F("OccCorrvsCent", "OccCorrvsCent", 101, -1, 100, 2000, 0 , 2);
  fOutput->Add(fHistOccCorrvsCent);
}

Bool_t AliAnalysisTaskRhoSparse::IsJetOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2)
{
  for (Int_t i = 0; i < jet1->GetNumberOfTracks(); ++i)
  {
    Int_t jet1Track = jet1->TrackAt(i);
    for (Int_t j = 0; j < jet2->GetNumberOfTracks(); ++j)
    {
      Int_t jet2Track = jet2->TrackAt(j);
      if (jet1Track == jet2Track)
        return kTRUE;
    }
  }
  return kFALSE;
}

Bool_t AliAnalysisTaskRhoSparse::IsJetSignal(AliEmcalJet* jet)
{
  if(jet->Pt()>5){
    return kTRUE;
  }else{
    return kFALSE;
  }
}

Bool_t AliAnalysisTaskRhoSparse::Run() 
{

  fOutRho->SetVal(0);
  if (fOutRhoScaled)
    fOutRhoScaled->SetVal(0);

  if (!fJets)
    return kFALSE;
  const Int_t Njets = fJets->GetEntries();

  AliJetContainer *sigjets = static_cast<AliJetContainer*>(fJetCollArray.At(1));
  Int_t NjetsSig = 0;
  if (sigjets) NjetsSig = sigjets->GetNJets();

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

      //Search the two leading KT jets
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
  Double_t TotaljetAreaPhys=0;
  Double_t TotalAreaCovered=0;
  Double_t TotalTPCArea=2*TMath::Pi()*0.9;

  // push all jets within selected acceptance into stack
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {

    AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(iJets));
    if (!jet) {
      AliError(Form("%s: Could not receive jet %d", GetName(), iJets));
      continue;
    } 

    // Take into account the area of real jets (no pure ghost jets)
    // and also the area of later excluded jets.
    // Some of these jets do not contribute to the rho calculation
    // but to the factor with which this rho is scaled down
    if (fExcludeAreaExcludedJets==0)
    {
      // Total area of physical jets that are used for the rho calculation
      if(jet->GetNumberOfTracks()>0)
      {
       TotaljetAreaPhys+=jet->Area();
      }
      // Total area of all found jets ghost+phyical jets. This is a proxy
      // for the available detector area in which the rho could have been calculated
      TotalAreaCovered+=jet->Area();
    }
    // Exclude leading background jets (could be signal)
    if (iJets == maxJetIds[0] || iJets == maxJetIds[1])
      continue;
    // Exclude background jets that do not fullfill basic cuts defined in AliJetContainer
    if (!AcceptJet(jet))
      continue;

    // Search for overlap with signal jets
    Bool_t isOverlapping = kFALSE;
    if (sigjets) {
      for(Int_t j=0;j<NjetsSig;j++)
      {
        AliEmcalJet* signalJet = sigjets->GetAcceptJet(j);
        if(!signalJet)
          continue;
        if(!IsJetSignal(signalJet))
          continue;

        if(IsJetOverlapping(signalJet, jet))
        {
          isOverlapping = kTRUE;
          break;
        }
      }
    }
    // Exclude background jets that overlap with anti-kT signal jets
    if(fExcludeOverlaps && isOverlapping)
      continue;

    // Take into account only real jets (no pure ghost jets) that also
    // contribute to the rho calculation
    if (fExcludeAreaExcludedJets==1)
    {
      // Total area of physical jets that are used for the rho calculation
      if(jet->GetNumberOfTracks()>0)
      {
       TotaljetAreaPhys+=jet->Area();
      }
      // Total area of all found jets ghost+phyical jets. This is a proxy
      // for the available detector area in which the rho could have been calculated
      TotalAreaCovered+=jet->Area();
    }

    // Exclude pure ghost jets from the rho calculation.
    // Use only jets that fulfill your background jet selection
    // for the rho calculation.
    // Eg. real signal jets should not bias the background rho
    if(jet->GetNumberOfTracks()>0)
    {
      rhovec[NjetAcc] = jet->Pt() / jet->Area();
      ++NjetAcc;
    }
  }

  Double_t OccCorr=1;
  //Use the total TPC area in which rho is calculated as denominater
  if(fUseTPCArea)
  {
    OccCorr = TotaljetAreaPhys/TotalTPCArea;
  }
  //Use the total area covered by all ghost+physical jets that pass the selection in the detector as denominator
  else if (TotalAreaCovered>0)
  {
    OccCorr = TotaljetAreaPhys/TotalAreaCovered;
  }

  if (fCreateHisto)
    fHistOccCorrvsCent->Fill(fCent, OccCorr);

  if (NjetAcc > 0) {
    //find median value
    Double_t rho = TMath::Median(NjetAcc, rhovec);

    if(fRhoCMS){
      rho = rho * OccCorr;
    }

    fOutRho->SetVal(rho);

    if (fOutRhoScaled) {
      Double_t rhoScaled = rho * GetScaleFactor(fCent);
      fOutRhoScaled->SetVal(rhoScaled);
    }
  }

  return kTRUE;
} 

AliAnalysisTaskRhoSparse* AliAnalysisTaskRhoSparse::AddTaskRhoSparse(
    const char    *nTracks,
    const char    *nClusters,
    const char    *nRho,
    Double_t       jetradius,
    UInt_t         acceptance,
    AliJetContainer::EJetType_t jetType,
    AliJetContainer::ERecoScheme_t rscheme,
    const Bool_t   histo,
    const char    *nJetsSig,
    const char    *cutType,
    Double_t       jetptcut,
    Double_t       jetareacut,
    Double_t       emcareacut,
    const char    *suffix
)
{

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRhoSparse", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalJetPerformance", "This task requires an input event handler");
    return 0;
  }

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

  TString name("AliAnalysisTaskRhoSparse");
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskRhoSparse* mgrTask = (AliAnalysisTaskRhoSparse*)(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;

  AliAnalysisTaskRhoSparse *rhotask = new AliAnalysisTaskRhoSparse(name, histo);
  rhotask->SetHistoBins(1000,-0.1,9.9);
  rhotask->SetRhoCMS(1);
  rhotask->SetSmallSystem(1);
  rhotask->SetOutRhoName(nRho);

  AliTrackContainer* trackCont;
  AliParticleContainer* partCont;
  if (trackName == "mcparticles") {
    AliMCParticleContainer* mcpartCont = rhotask->AddMCParticleContainer(trackName);
    mcpartCont->SelectPhysicalPrimaries(kTRUE);
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    trackCont = rhotask->AddTrackContainer(trackName);
  }
  else if (!trackName.IsNull()) {
    rhotask->AddParticleContainer(trackName);
  }
  partCont = rhotask->GetParticleContainer(0);

  AliClusterContainer *clusterCont = rhotask->AddClusterContainer(clusName);
  if (clusterCont) {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(0.3);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  //
  AliJetContainer *bkgJetCont = rhotask->AddJetContainer(jetType, AliJetContainer::kt_algorithm, rscheme, jetradius, acceptance, partCont, clusterCont);
  if (bkgJetCont) {
    //why?? bkgJetCont->SetJetAreaCut(jetareacut);
    //why?? bkgJetCont->SetAreaEmcCut(emcareacut);
    bkgJetCont->SetJetPtCut(0.);
  }

  //This is only for excluding overlaps with signal jets in case signal jets are provided
  AliJetContainer *sigJetCont;
  if(nJetsSig!=0)
  {
    sigJetCont = rhotask->AddJetContainer(nJetsSig,cutType,jetradius);
    if (sigJetCont)
    {
      sigJetCont->SetJetAreaCut(jetareacut);
      sigJetCont->SetAreaEmcCut(emcareacut);
      sigJetCont->SetJetPtCut(jetptcut);
      sigJetCont->ConnectParticleContainer(trackCont);
      sigJetCont->ConnectClusterContainer(clusterCont);
    }
  }
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
        TList::Class(),AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(rhotask, 1, coutput1);
  }

  return rhotask;
}

