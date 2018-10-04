/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliAnalysisTaskConvJet.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskConvJet);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskConvJet::AliAnalysisTaskConvJet() : 
  AliAnalysisTaskEmcalJet(),
//   fHistManager(),
  fNJets(0),
  fVectorJetPt(0),
  fVectorJetEta(0),
  fVectorJetPhi(0),
  fVectorJetR(0)
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskConvJet::AliAnalysisTaskConvJet(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
//   fHistManager(name),
  fNJets(0),
  fVectorJetPt(0),
  fVectorJetEta(0),
  fVectorJetPhi(0),
  fVectorJetR(0)
{
  SetMakeGeneralHistograms(kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskConvJet::~AliAnalysisTaskConvJet()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskConvJet::UserCreateOutputObjects()
{
   AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskConvJet::FillHistograms()
{
  DoJetLoop();

  return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskConvJet::DoJetLoop()
{
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    UInt_t count = 0;
    fNJets = 0 ;
    fVectorJetPt.clear();
    fVectorJetEta.clear();
    fVectorJetPhi.clear();
    fVectorJetR.clear();
    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      count++;
      fVectorJetPt.push_back(jet->Pt());
      fVectorJetEta.push_back(jet->Eta());
      fVectorJetPhi.push_back(jet->Phi());
      fVectorJetR.push_back(jet->Area());
    }
    fNJets = count ;
  }
}

/**
 * This function performs a loop over the reconstructed tracks
 * in the current event and fills the relevant histograms.
 */
// void AliAnalysisTaskConvJet::DoTrackLoop()
// {
//   AliClusterContainer* clusCont = GetClusterContainer(0);
// 
//   TString histname;
//   TString groupname;
//   UInt_t sumAcceptedTracks = 0;
//   AliParticleContainer* partCont = 0;
//   TIter next(&fParticleCollArray);
//   while ((partCont = static_cast<AliParticleContainer*>(next()))) {
//     groupname = partCont->GetName();
//     UInt_t count = 0;
//     for(auto part : partCont->accepted()) {
//       if (!part) continue;
//       count++;
// 
//       histname = TString::Format("%s/histTrackPt", groupname.Data());
//       fHistManager.FillTH1(histname, part->Pt());
// 
//       histname = TString::Format("%s/histTrackPhi", groupname.Data());
//       fHistManager.FillTH1(histname, part->Phi());
// 
//       histname = TString::Format("%s/histTrackEta", groupname.Data());
//       fHistManager.FillTH1(histname, part->Eta());
// 
//       if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
//         const AliVTrack* track = static_cast<const AliVTrack*>(part);
// 
//         histname = TString::Format("%s/fHistDeltaEtaPt", groupname.Data());
//         fHistManager.FillTH1(histname, track->Pt(), track->Eta() - track->GetTrackEtaOnEMCal());
// 
//         histname = TString::Format("%s/fHistDeltaPhiPt", groupname.Data());
//         fHistManager.FillTH1(histname, track->Pt(), track->Phi() - track->GetTrackPhiOnEMCal());
// 
//         histname = TString::Format("%s/fHistDeltaPtvsPt", groupname.Data());
//         fHistManager.FillTH1(histname, track->Pt(), track->Pt() - track->GetTrackPtOnEMCal());
// 
//         if (clusCont) {
//           Int_t iCluster = track->GetEMCALcluster();
//           if (iCluster >= 0) {
//             AliVCluster* cluster = clusCont->GetAcceptCluster(iCluster);
//             if (cluster) {
//               histname = TString::Format("%s/fHistEoverPvsP", groupname.Data());
//               fHistManager.FillTH2(histname, track->P(), cluster->GetNonLinCorrEnergy() / track->P());
//             }
//           }
//         }
//       }
//     }
//     sumAcceptedTracks += count;
// 
//     histname = TString::Format("%s/histNTracks", groupname.Data());
//     fHistManager.FillTH1(histname, count);
//   }
// 
//   histname = "fHistSumNTracks";
//   fHistManager.FillTH1(histname, sumAcceptedTracks);
// }

/**
 * This function performs a loop over the reconstructed EMCal clusters
 * in the current event and fills the relevant histograms.
 */
// void AliAnalysisTaskConvJet::DoClusterLoop()
// {
//   TString histname;
//   TString groupname;
//   UInt_t sumAcceptedClusters = 0;
//   AliClusterContainer* clusCont = 0;
//   TIter next(&fClusterCollArray);
//   while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
//     groupname = clusCont->GetName();
// 
//     for(auto cluster : clusCont->all()) {
//       if (!cluster) continue;
// 
//       if (cluster->GetIsExotic()) {
//         histname = TString::Format("%s/histClusterEnergyExotic", groupname.Data());
//         fHistManager.FillTH1(histname, cluster->E());
//       }
//     }
// 
//     UInt_t count = 0;
//     for(auto cluster : clusCont->accepted()) {
//       if (!cluster) continue;
//       count++;
// 
//       AliTLorentzVector nPart;
//       cluster->GetMomentum(nPart, fVertex);
// 
//       histname = TString::Format("%s/histClusterEnergy", groupname.Data());
//       fHistManager.FillTH1(histname, cluster->E());
// 
//       histname = TString::Format("%s/histClusterNonLinCorrEnergy", groupname.Data());
//       fHistManager.FillTH1(histname, cluster->GetNonLinCorrEnergy());
// 
//       histname = TString::Format("%s/histClusterHadCorrEnergy", groupname.Data());
//       fHistManager.FillTH1(histname, cluster->GetHadCorrEnergy());
// 
//       histname = TString::Format("%s/histClusterPhi", groupname.Data());
//       fHistManager.FillTH1(histname, nPart.Phi_0_2pi());
// 
//       histname = TString::Format("%s/histClusterEta", groupname.Data());
//       fHistManager.FillTH1(histname, nPart.Eta());
//     }
//     sumAcceptedClusters += count;
// 
//     histname = TString::Format("%s/histNClusters", groupname.Data());
//     fHistManager.FillTH1(histname, count);
//   }
// 
//   histname = "fHistSumNClusters";
//   fHistManager.FillTH1(histname, sumAcceptedClusters);
// }
// 
// /**
//  * This function performs a loop over the reconstructed EMCal cells
//  * in the current event and fills the relevant histograms.
//  */
// void AliAnalysisTaskConvJet::DoCellLoop()
// {
//   if (!fCaloCells) return;
// 
//   TString histname;
// 
//   const Short_t ncells = fCaloCells->GetNumberOfCells();
// 
//   histname = TString::Format("%s/histNCells", fCaloCellsName.Data());
//   fHistManager.FillTH1(histname, ncells);
// 
//   histname = TString::Format("%s/histCellEnergy", fCaloCellsName.Data());
//   for (Short_t pos = 0; pos < ncells; pos++) {
//     Double_t amp   = fCaloCells->GetAmplitude(pos);
//     
//     fHistManager.FillTH1(histname, amp);
//   }
// }

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskConvJet::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskConvJet::Run()
{
  return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskConvJet::Terminate(Option_t *) 
{
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskConvJet * AliAnalysisTaskConvJet::AddTask_GammaConvJet(
  const char *ntracks,
  const char *nclusters,
  const char* ncells,
  const char *suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTask_GammaConvJet", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTask_GammaConvJet", "This task requires an input event handler");
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

  TString trackName(ntracks);
  TString clusName(nclusters);
  TString cellName(ncells);

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

  if (cellName == "usedefault") {
    if (dataType == kESD) {
      cellName = "EMCALCells";
    }
    else if (dataType == kAOD) {
      cellName = "emcalCells";
    }
    else {
      cellName = "";
    }
  }

  TString name("AliAnalysisTaskConvJet");
//   if (!trackName.IsNull()) {
//     name += "_";
//     name += trackName;
//   }
//   if (!clusName.IsNull()) {
//     name += "_";
//     name += clusName;
//   }
//   if (!cellName.IsNull()) {
//     name += "_";
//     name += cellName;
//   }
//   if (strcmp(suffix,"") != 0) {
//     name += "_";
//     name += suffix;
//   }

  AliAnalysisTaskConvJet* sampleTask = new AliAnalysisTaskConvJet(name);
  sampleTask->SetCaloCellsName(cellName);
  sampleTask->SetVzRange(-10,10);

  if (trackName == "mcparticles") {
    sampleTask->AddMCParticleContainer(trackName);
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    sampleTask->AddTrackContainer(trackName);
  }
  else if (!trackName.IsNull()) {
    sampleTask->AddParticleContainer(trackName);
  }
  sampleTask->AddClusterContainer(clusName);
  
  sampleTask->GetClusterContainer(0)->SetClusECut(0.);
  sampleTask->GetClusterContainer(0)->SetClusPtCut(0.);
  sampleTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
  sampleTask->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
  sampleTask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  sampleTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
  sampleTask->GetParticleContainer(0)->SetParticleEtaLimits(-0.8,0.8);
  sampleTask->GetTrackContainer(0)->SetParticlePtCut(0.15);
  sampleTask->GetTrackContainer(0)->SetParticleEtaLimits(-0.8,0.8);
  sampleTask->GetTrackContainer(0)->SetFilterHybridTracks(kTRUE);
//   sampleTask->SetHistoBins(600, 0, 300);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(sampleTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
      TList::Class(),AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (sampleTask, 0,  cinput1 );
  mgr->ConnectOutput (sampleTask, 1, coutput1 );

  return sampleTask;
}
