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

#include "AliAnalysisTaskJetOutlierRemoval.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetOutlierRemoval);
/// \endcond


AliAnalysisTaskJetOutlierRemoval::AliAnalysisTaskJetOutlierRemoval() :
  AliAnalysisTaskEmcalJet(),
  fMaxJetPt(0)
{
}

AliAnalysisTaskJetOutlierRemoval::AliAnalysisTaskJetOutlierRemoval(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fMaxJetPt(0)
{
  SetMakeGeneralHistograms(kTRUE);
}


AliAnalysisTaskJetOutlierRemoval::~AliAnalysisTaskJetOutlierRemoval()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskJetOutlierRemoval::UserCreateOutputObjects()
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
Bool_t AliAnalysisTaskJetOutlierRemoval::FillHistograms()
{
  fMaxJetPt = -1;
  DoJetLoop();

  return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the vectors.
 */
void AliAnalysisTaskJetOutlierRemoval::DoJetLoop()
{
  AliDebugStream(1) << "Using custom MC outlier rejection" << std::endl;
  auto partjets = GetJetContainer("partjets");
  if(!partjets){
    AliFatal("MC jet container not found for jet pt determination used in outlier removal!");
  }

  // Loop over particle level jets and find max pT jet in TPC acceptance
  auto jetiter = partjets->accepted();
  auto max = std::max_element(jetiter.begin(), jetiter.end(), [](const AliEmcalJet *lhs, const AliEmcalJet *rhs ) { return lhs->Pt() < rhs->Pt(); });
  if(max != jetiter.end())  {
    fMaxJetPt = (*max)->Pt();
    // std::cout << "Found max jet with pt " << (*max)->Pt() << " GeV/c" << std::endl;
    // std::cout << "jetfinder max jet eta " << (*max)->Eta() << std::endl;
  }
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskJetOutlierRemoval::ExecOnce()
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
Bool_t AliAnalysisTaskJetOutlierRemoval::Run()
{
  return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskJetOutlierRemoval::Terminate(Option_t *) 
{
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskJetOutlierRemoval * AliAnalysisTaskJetOutlierRemoval::AddTask_GammaOutlierRemoval()
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTask_GammaOutlierRemoval", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTask_GammaOutlierRemoval", "This task requires an input event handler");
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
  TString name("AliAnalysisTaskJetOutlierRemoval");

  AliAnalysisTaskJetOutlierRemoval* task = new AliAnalysisTaskJetOutlierRemoval(name);

  AliTrackContainer *tracks(nullptr);
  AliClusterContainer *clusters(nullptr);
    tracks = task->AddTrackContainer(dataType == kAOD ? "tracks" : "Tracks");
    tracks->SetMinPt(0.15);

    clusters = task->AddClusterContainer(dataType == kAOD ? "caloClusters" : "CaloClusters");
    clusters->SetDefaultClusterEnergy(0.2);
    clusters->SetClusUserDefEnergyCut(0.2, 0.3);




  auto partcont = task->AddMCParticleContainer("mcparticlesSelected");
  partcont->SetMinPt(0.);
  
  auto pjcont = task->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::E_scheme, 0.2, AliJetContainer::kTPCfid, partcont, nullptr);
  pjcont->SetName("partjets");
  pjcont->SetMinPt(0);
  pjcont->SetMaxTrackPt(1000.);
  std::cout << "Adding jet container with underlying array:" << pjcont->GetArrayName() << std::endl;


  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname("mcparticlejet");
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
      TList::Class(),AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (task, 0,  cinput1 );
  mgr->ConnectOutput (task, 1, coutput1 );

  return task;
}
