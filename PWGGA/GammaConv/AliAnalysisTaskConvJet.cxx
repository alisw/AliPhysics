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


AliAnalysisTaskConvJet::AliAnalysisTaskConvJet() :
  AliAnalysisTaskEmcalJet(),
  fNJets(0),
  fVectorJetPt(0),
  fVectorJetPx(0),
  fVectorJetPy(0),
  fVectorJetPz(0),
  fVectorJetEta(0),
  fVectorJetPhi(0),
  fVectorJetR(0),
  fTrueNJets(0),
  fTrueVectorJetPt(0),
  fTrueVectorJetPx(0),
  fTrueVectorJetPy(0),
  fTrueVectorJetPz(0),
  fTrueVectorJetEta(0),
  fTrueVectorJetPhi(0),
  fTrueVectorJetR(0)
{
}

AliAnalysisTaskConvJet::AliAnalysisTaskConvJet(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fNJets(0),
  fVectorJetPt(0),
  fVectorJetPx(0),
  fVectorJetPy(0),
  fVectorJetPz(0),
  fVectorJetEta(0),
  fVectorJetPhi(0),
  fVectorJetR(0),
  fTrueVectorJetPt(0),
  fTrueVectorJetPx(0),
  fTrueVectorJetPy(0),
  fTrueVectorJetPz(0),
  fTrueVectorJetEta(0),
  fTrueVectorJetPhi(0),
  fTrueVectorJetR(0)
{
  SetMakeGeneralHistograms(kTRUE);
}


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
 * The body of this function should contain instructions to fill the vectors.
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
 * in the current event and fills the vectors.
 */
void AliAnalysisTaskConvJet::DoJetLoop()
{
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    TString JetName = jetCont->GetTitle();
    TObjArray *arr = JetName.Tokenize("__");
    TObjString* testObjString = (TObjString*)arr->At(2);
    if(testObjString->GetString() != "mcparticles"){
      UInt_t count = 0;
      fNJets = 0 ;
      fVectorJetPt.clear();
      fVectorJetPx.clear();
      fVectorJetPy.clear();
      fVectorJetPz.clear();
      fVectorJetEta.clear();
      fVectorJetPhi.clear();
      fVectorJetR.clear();
      for(auto const& jet : jetCont->accepted()) {
        if (!jet) continue;
        count++;
        fVectorJetPt.push_back(jet->Pt());
        fVectorJetPx.push_back(jet->Px());
        fVectorJetPy.push_back(jet->Py());
        fVectorJetPz.push_back(jet->Pz());
        fVectorJetEta.push_back(jet->Eta());
        fVectorJetPhi.push_back(jet->Phi());
        fVectorJetR.push_back(jet->Area());
      }
      fNJets = count ;
    }else{
      UInt_t count = 0;
      fTrueNJets = 0 ;
      fTrueVectorJetPt.clear();
      fTrueVectorJetPx.clear();
      fTrueVectorJetPy.clear();
      fTrueVectorJetPz.clear();
      fTrueVectorJetEta.clear();
      fTrueVectorJetPhi.clear();
      fTrueVectorJetR.clear();
      for(auto const& jet : jetCont->accepted()) {
        if (!jet) continue;
        count++;
        fTrueVectorJetPt.push_back(jet->Pt());
        fTrueVectorJetPx.push_back(jet->Px());
        fTrueVectorJetPy.push_back(jet->Py());
        fTrueVectorJetPz.push_back(jet->Pz());
        fTrueVectorJetEta.push_back(jet->Eta());
        fTrueVectorJetPhi.push_back(jet->Phi());
        fTrueVectorJetR.push_back(jet->Area());
      }
      fTrueNJets = count ;
    }
  }
}

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

  if(trackName != "mcparticles"){
      sampleTask->GetTrackContainer(0)->SetFilterHybridTracks(kTRUE);
      sampleTask->GetTrackContainer(0)->SetParticlePtCut(0.15);
      sampleTask->GetTrackContainer(0)->SetParticleEtaLimits(-0.8,0.8);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(sampleTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(trackName);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
      TList::Class(),AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (sampleTask, 0,  cinput1 );
  mgr->ConnectOutput (sampleTask, 1, coutput1 );

  return sampleTask;
}
