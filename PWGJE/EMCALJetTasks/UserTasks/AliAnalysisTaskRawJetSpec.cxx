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

#include "AliAnalysisTaskRawJetSpec.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskRawJetSpec);
/// \endcond

// Default constructor. Needed by ROOT I/O
AliAnalysisTaskRawJetSpec::AliAnalysisTaskRawJetSpec() : 
  AliAnalysisTaskEmcalJet(),
  fHistManager()
{
}

/**
 * Standard constructor. Should be used by the user.
 * @param[in] name Name of the task
 */
AliAnalysisTaskRawJetSpec::AliAnalysisTaskRawJetSpec(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistManager(name)
{
  SetMakeGeneralHistograms(kTRUE);
}

// Destructor
AliAnalysisTaskRawJetSpec::~AliAnalysisTaskRawJetSpec(){}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskRawJetSpec::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AllocateJetHistograms();

  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskRawJetSpec::AllocateJetHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    //for (Int_t cent = 0; cent < fNcentBins; cent++) {
    for (Int_t cent = 0; cent < 10; cent++) {
      histname = TString::Format("%s/histJetPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt);

      histname = TString::Format("%s/histJetArea_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{A}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, 3);

      histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histJetEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      histname = TString::Format("%s/histNJets_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of jets;events", histname.Data());
      if (fForceBeamType != kpp) fHistManager.CreateTH1(histname, histtitle, 500, 0, 500);
      else                                   fHistManager.CreateTH1(histname, histtitle, 100, 0, 100);

      if (!jetCont->GetRhoName().IsNull()) {
        histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH1(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
      }
    }
  }
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskRawJetSpec::FillHistograms()
{
  DoJetLoop();

  return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskRawJetSpec::DoJetLoop()
{
  TString histname;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);

  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    UInt_t count = 0;
    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      count++;

      //histname = TString::Format("%s/histJetPt_%d", groupname.Data(), fCentBin);
      histname = TString::Format("%s/histJetPt_%d", groupname.Data(), int(fCent/10));
      fHistManager.FillTH1(histname, jet->Pt());

      histname = TString::Format("%s/histJetArea_%d", groupname.Data(), int(fCent/10));
      fHistManager.FillTH1(histname, jet->Area());

      histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), int(fCent/10));
      fHistManager.FillTH1(histname, jet->Phi());

      histname = TString::Format("%s/histJetEta_%d", groupname.Data(), int(fCent/10));
      fHistManager.FillTH1(histname, jet->Eta());

      if (jetCont->GetRhoParameter()) {
        histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), int(fCent/10));
        fHistManager.FillTH1(histname, jet->Pt() - jetCont->GetRhoVal() * jet->Area());
      }
    }
    histname = TString::Format("%s/histNJets_%d", groupname.Data(), int(fCent/10));
    fHistManager.FillTH1(histname, count);
  }
}

void AliAnalysisTaskRawJetSpec::ExecOnce()
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
Bool_t AliAnalysisTaskRawJetSpec::Run()
{
  return kTRUE;
}


 // This function is called once at the end of the analysis.
void AliAnalysisTaskRawJetSpec::Terminate(Option_t *){}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskRawJetSpec * AliAnalysisTaskRawJetSpec::AddTaskRawJetSpec(
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
    ::Error("AddTaskRawJetSpec", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskRawJetSpec", "This task requires an input event handler");
    return 0;
  }

  enum EDataType_t { kAOD};
  EDataType_t dataType = kAOD;

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString trackName(ntracks);
  TString clusName(nclusters);
  TString cellName(ncells);

  if (trackName == "usedefault" && dataType == kAOD) trackName = "tracks";
  if (clusName == "usedefault" && dataType == kAOD) clusName = "caloClusters";
  if (cellName == "usedefault" && dataType == kAOD) cellName = "emcalCells";

  TString name("AliAnalysisTaskRawJetSpec");

  AliAnalysisTaskRawJetSpec* sampleTask = new AliAnalysisTaskRawJetSpec(name);
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

  // Final settings, pass to manager and set the containers  ---------------------------------------------------
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
