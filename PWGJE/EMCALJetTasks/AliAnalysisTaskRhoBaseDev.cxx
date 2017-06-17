/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
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

#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TClonesArray.h>
#include <TGrid.h>

#include <AliLog.h>
#include <AliVEventHandler.h>
#include <AliAnalysisManager.h>

#include "AliRhoParameter.h"
#include "AliEmcalJet.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliJetContainer.h"

#include "AliAnalysisTaskRhoBaseDev.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskRhoBaseDev);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskRhoBaseDev::AliAnalysisTaskRhoBaseDev() :
  AliAnalysisTaskEmcalJetLight("AliAnalysisTaskRhoBaseDev", kFALSE),
  fOutRhoName(),
  fOutRhoScaledName(),
  fRhoFunction(0),
  fScaleFunction(0),
  fAttachToEvent(kTRUE),
  fNbins(250),
  fMinBinPt(0),
  fMaxBinPt(500),
  fOutRho(0),
  fOutRhoScaled(0),
  fHistJetPtvsCent(0),
  fHistJetAreavsCent(0),
  fHistJetRhovsCent(0),
  fHistNjetvsCent(0),
  fHistJetPtvsNtrack(0),
  fHistJetAreavsNtrack(0),
  fHistNjetvsNtrack(0),
  fHistRhovsCent(0),
  fHistRhoScaledvsCent(0),
  fHistRhovsNtrack(0),
  fHistRhoScaledvsNtrack(0),
  fHistRhovsNcluster(0),
  fHistRhoScaledvsNcluster(0)
{
  for (Int_t i = 0; i < 4; i++) {
    fHistJetNconstVsPt[i] = 0;
    fHistJetRhovsEta[i] = 0;
  }
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name  Name of the task
 * @param[in] histo If kTRUE, the task will also produce QA histograms
 */
AliAnalysisTaskRhoBaseDev::AliAnalysisTaskRhoBaseDev(const char *name, Bool_t histo) :
    AliAnalysisTaskEmcalJetLight(name, histo),
  fOutRhoName(),
  fOutRhoScaledName(),
  fRhoFunction(0),
  fScaleFunction(0),
  fAttachToEvent(kTRUE),
  fNbins(250),
  fMinBinPt(0),
  fMaxBinPt(500),
  fOutRho(0),
  fOutRhoScaled(0),
  fHistJetPtvsCent(0),
  fHistJetAreavsCent(0),
  fHistJetRhovsCent(0),
  fHistNjetvsCent(0),
  fHistJetPtvsNtrack(0),
  fHistJetAreavsNtrack(0),
  fHistNjetvsNtrack(0),
  fHistRhovsCent(0),
  fHistRhoScaledvsCent(0),
  fHistRhovsNtrack(0),
  fHistRhoScaledvsNtrack(0),
  fHistRhovsNcluster(0),
  fHistRhoScaledvsNcluster(0)
{
  for (Int_t i = 0; i < 4; i++) {
    fHistJetNconstVsPt[i] = 0;
    fHistJetRhovsEta[i] = 0;
  }
  SetMakeGeneralHistograms(histo);
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskRhoBaseDev::UserCreateOutputObjects()
{
  if (!fCreateHisto) return;

  AliAnalysisTaskEmcalJetLight::UserCreateOutputObjects();
  //ranges for PbPb
  Float_t Ntrackrange[2] = {0, 6000};
  //set multiplicity related axes to a smaller max value
  if (fBeamType != kAA) Ntrackrange[1] = 200.;
  
  fHistRhovsCent = new TH2F("fHistRhovsCent", "fHistRhovsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt);
  fHistRhovsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhovsCent->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
  fOutput->Add(fHistRhovsCent);

  if (fParticleCollArray.size() > 0) {
    fHistRhovsNtrack = new TH2F("fHistRhovsNtrack", "fHistRhovsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt);
    fHistRhovsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistRhovsNtrack->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhovsNtrack);
  }

  if (fClusterCollArray.size()>0) {
    fHistRhovsNcluster = new TH2F("fHistRhovsNcluster", "fHistRhovsNcluster", 50, 0, 1500, fNbins, fMinBinPt, fMaxBinPt);
    fHistRhovsNcluster->GetXaxis()->SetTitle("No. of clusters");
    fHistRhovsNcluster->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhovsNcluster);
  }

  if (fJetCollArray.size() > 0) {
    fHistJetPtvsCent = new TH2F("fHistJetPtvsCent", "fHistJetPtvsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt);
    fHistJetPtvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetPtvsCent->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
    fOutput->Add(fHistJetPtvsCent);

    fHistJetAreavsCent = new TH2F("fHistJetAreavsCent", "fHistJetAreavsCent", 101, -1, 100, 100, 0, 1);
    fHistJetAreavsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetAreavsCent->GetYaxis()->SetTitle("Jet area");
    fOutput->Add(fHistJetAreavsCent);

    fHistJetRhovsCent = new TH2F("fHistJetRhovsCent", "fHistJetRhovsCent", 101, -1, 100, fNbins, fMinBinPt, fMaxBinPt);
    fHistJetRhovsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetRhovsCent->GetYaxis()->SetTitle("Jet #rho (GeV/c)");
    fOutput->Add(fHistJetRhovsCent);

    fHistNjetvsCent = new TH2F("fHistNjetvsCent",  "fHistNjetvsCent", 101, -1, 100, 150, -0.5, 149.5);
    fHistNjetvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistNjetvsCent->GetYaxis()->SetTitle("No. of jets");
    fOutput->Add(fHistNjetvsCent);


    if (fParticleCollArray.size() > 0) {
      fHistJetPtvsNtrack = new TH2F("fHistJetPtvsNtrack", "fHistJetPtvsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt);
      fHistJetPtvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistJetPtvsNtrack->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
      fOutput->Add(fHistJetPtvsNtrack);

      fHistJetAreavsNtrack = new TH2F("fHistJetAreavsNtrack", "fHistJetAreavsNtrack", 150, Ntrackrange[0], Ntrackrange[1], 100, 0, 1);
      fHistJetAreavsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistJetAreavsNtrack->GetYaxis()->SetTitle("Jet area");
      fOutput->Add(fHistJetAreavsNtrack);

      fHistNjetvsNtrack = new TH2F("fHistNjetvsNtrack", "fHistNjetvsNtrack", 150, Ntrackrange[0], Ntrackrange[1], 150, -0.5, 149.5);
      fHistNjetvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistNjetvsNtrack->GetYaxis()->SetTitle("No. of jets");
      fOutput->Add(fHistNjetvsNtrack);
    }


    TString name;
    for (Int_t i = 0; i < 4; i++) {
      name = Form("fHistJetNconstVsPt_%d",i);
      fHistJetNconstVsPt[i] = new TH2F(name, name, 150, -0.5, 149.5, fNbins, fMinBinPt, fMaxBinPt);
      fHistJetNconstVsPt[i]->GetXaxis()->SetTitle("No. of constituents");
      fHistJetNconstVsPt[i]->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");
      fOutput->Add(fHistJetNconstVsPt[i]);

      name = Form("fHistJetRhovsEta_%d",i);
      fHistJetRhovsEta[i] = new TH2F(name, name, fNbins, fMinBinPt, fMaxBinPt, 16, -0.8, 0.8);
      fHistJetRhovsEta[i]->GetXaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
      fHistJetRhovsEta[i]->GetYaxis()->SetTitle("#eta");
      fOutput->Add(fHistJetRhovsEta[i]);
    }
  }

  if (fScaleFunction) {
    fHistRhoScaledvsCent = new TH2F("fHistRhoScaledvsCent", "fHistRhoScaledvsCent", 101, -1, 100, fNbins, fMinBinPt , fMaxBinPt);
    fHistRhoScaledvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistRhoScaledvsCent->GetYaxis()->SetTitle("#rho_{scaled} (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhoScaledvsCent);

    if (fParticleCollArray.size() > 0) {
      fHistRhoScaledvsNtrack = new TH2F("fHistRhoScaledvsNtrack", "fHistRhoScaledvsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt);
      fHistRhoScaledvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistRhoScaledvsNtrack->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
      fOutput->Add(fHistRhoScaledvsNtrack);
    }

    if (fClusterCollArray.size() > 0) {
      fHistRhoScaledvsNcluster = new TH2F("fHistRhoScaledvsNcluster", "fHistRhoScaledvsNcluster", 50, 0, 1500, fNbins, fMinBinPt, fMaxBinPt);
      fHistRhoScaledvsNcluster->GetXaxis()->SetTitle("No. of clusters");
      fHistRhoScaledvsNcluster->GetYaxis()->SetTitle("#rho_{scaled} (GeV/c * rad^{-1})");
      fOutput->Add(fHistRhoScaledvsNcluster);
    }
  }
}

/**
 * Run the analysis.
 */
Bool_t AliAnalysisTaskRhoBaseDev::Run()
{
  Double_t rho = GetRhoFactor(fCent);
  fOutRho->SetVal(rho);

  if (fScaleFunction) {
    Double_t rhoScaled = rho * GetScaleFactor(fCent);
    fOutRhoScaled->SetVal(rhoScaled);
  }

  return kTRUE;
}

/**
 * Fill histograms.
 */
Bool_t AliAnalysisTaskRhoBaseDev::FillHistograms()
{
  Int_t Ntracks   = 0;
  Int_t Nclusters = 0;
  Int_t NjetAcc   = 0;

  // Loop over all possible containers
  for (auto partCont : fParticleCollArray) {
    Ntracks += partCont.second->GetNAcceptedParticles();
  }
  for (auto clusCont : fClusterCollArray) {
    Nclusters += clusCont.second->GetNAcceptedClusters();
  }

  AliJetContainer* jetCont = fJetCollArray["Background"];
  for (auto jet : jetCont->accepted()) {
    fHistJetPtvsCent->Fill(fCent, jet->Pt());
    fHistJetAreavsCent->Fill(fCent, jet->Area());
    fHistJetRhovsCent->Fill(fCent, jet->Pt() / jet->Area());
    fHistJetRhovsEta[fCentBin]->Fill(jet->Pt() / jet->Area(), jet->Eta());

    if (fHistJetPtvsNtrack) fHistJetPtvsNtrack->Fill(Ntracks, jet->Pt());
    if (fHistJetAreavsNtrack) fHistJetAreavsNtrack->Fill(Ntracks, jet->Area());

    fHistJetNconstVsPt[fCentBin]->Fill(jet->GetNumberOfConstituents(), jet->Pt());
    NjetAcc++;
  }

  if (fHistNjetvsCent) fHistNjetvsCent->Fill(fCent, NjetAcc);
  if (fHistNjetvsNtrack) fHistNjetvsNtrack->Fill(Ntracks, NjetAcc);

  fHistRhovsCent->Fill(fCent, fOutRho->GetVal());

  if (fHistRhovsNtrack) fHistRhovsNtrack->Fill(Ntracks, fOutRho->GetVal());
  if (fHistRhovsNcluster) fHistRhovsNcluster->Fill(Nclusters, fOutRho->GetVal());

  if (fOutRhoScaled) {
    fHistRhoScaledvsCent->Fill(fCent, fOutRhoScaled->GetVal());
    if (fHistRhoScaledvsNtrack) fHistRhoScaledvsNtrack->Fill(Ntracks, fOutRhoScaled->GetVal());
    if (fHistRhoScaledvsNcluster) fHistRhoScaledvsNcluster->Fill(Nclusters,  fOutRhoScaled->GetVal());
  }

  return kTRUE;
}      

/**
 * Init the analysis.
 */
void AliAnalysisTaskRhoBaseDev::ExecOnce()
{
  if (!fOutRho) {
    fOutRho = new AliRhoParameter(fOutRhoName, 0);

    if (fAttachToEvent) {
      if (!(InputEvent()->FindListObject(fOutRhoName))) {
        InputEvent()->AddObject(fOutRho);
      } else {
        AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fOutRhoName.Data()));
        return;
      }
    }
  }

  if (fScaleFunction && !fOutRhoScaled) {
    fOutRhoScaled = new AliRhoParameter(fOutRhoScaledName, 0);

    if (fAttachToEvent) {
      if (!(InputEvent()->FindListObject(fOutRhoScaledName))) {
        InputEvent()->AddObject(fOutRhoScaled);
      } else {
        AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fOutRhoScaledName.Data()));
        return;
      }
    }
  }

  AliAnalysisTaskEmcalJetLight::ExecOnce();
}

/**
 * Return rho per centrality.
 * @param cent Centrality percentile.
 * @return The average background from the parametrization loaded in memory.
 */
Double_t AliAnalysisTaskRhoBaseDev::GetRhoFactor(Double_t cent)
{
  Double_t rho = 0;
  if (fRhoFunction) rho = fRhoFunction->Eval(cent);
  return rho;
}

/**
 * Get scale factor.
 * @param cent Centrality percentile.
 * @return The scale factor (charged -> full) from the parametrization loaded in memory.
 */
Double_t AliAnalysisTaskRhoBaseDev::GetScaleFactor(Double_t cent)
{
  Double_t scale = 1;
  if (fScaleFunction) scale = fScaleFunction->Eval(cent);
  return scale;
}

/**
 * Load the scale function from a file.
 * @param path Path to the file
 * @param name Name of the object inside the file
 * @return On success, the pointer to the TF1 object
 */
TF1* AliAnalysisTaskRhoBaseDev::LoadRhoFunction(const char* path, const char* name)
{
  TString fname(path);
  if (fname.BeginsWith("alien://")) {
    TGrid::Connect("alien://");
  }

  TFile* file = TFile::Open(path);

  if (!file || file->IsZombie()) {
    ::Error("AddTaskRho", "Could not open scale function file");
    return 0;
  }

  TF1* sfunc = dynamic_cast<TF1*>(file->Get(name));

  if (sfunc) {
    ::Info("AliAnalysisTaskRhoBaseDev::LoadRhoFunction", "Scale function %s loaded from file %s.", name, path);
  }
  else {
    ::Error("AliAnalysisTaskRhoBaseDev::LoadRhoFunction", "Scale function %s not found in file %s.", name, path);
    return 0;
  }

  fScaleFunction = static_cast<TF1*>(sfunc->Clone());

  file->Close();
  delete file;

  return fScaleFunction;
}

/**
 * Create an instance of this class and add it to the analysis manager
 * @param trackName name of the track collection
 * @param clusName name of the calorimeter cluster collection
 * @param nRho name of the output rho object
 * @param jetradius Radius of the kt jets used to calculate the background
 * @param acceptance Fiducial acceptance of the kt jets
 * @param jetType Jet type (full/charged)
 * @param rscheme Recombination scheme
 * @param histo If kTRUE the task will also produce QA histograms
 * @param suffix additional suffix that can be added at the end of the task name
 * @return pointer to the new AliAnalysisTaskRhoDev task
 */
AliAnalysisTaskRhoBaseDev* AliAnalysisTaskRhoBaseDev::AddTaskRhoBaseDev(TString trackName, TString clusName, TString nRho, Double_t jetradius, UInt_t acceptance, AliJetContainer::EJetType_t jetType , AliJetContainer::ERecoScheme_t rscheme, Bool_t histo, TString suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliAnalysisTaskRhoBaseDev::AddTaskRhoBaseDev", "No analysis manager to connect to.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AliAnalysisTaskRhoBaseDev::AddTaskRhoBaseDev", "This task requires an input event handler");
    return nullptr;
  }

  EDataType_t dataType = kUnknownDataType;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  // Init the task and do settings
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

  TString name("AliAnalysisTaskRhoBaseDev");
  if (!suffix.IsNull()) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskRhoBaseDev* mgrTask = dynamic_cast<AliAnalysisTaskRhoBaseDev*>(mgr->GetTask(name.Data()));
  if (mgrTask) {
    ::Warning("AliAnalysisTaskRhoBaseDev::AddTaskRhoBaseDev", "Not adding the task again, since a task with the same name '%s' already exists", name.Data());
    return mgrTask;
  }

  AliAnalysisTaskRhoBaseDev* rhotask = new AliAnalysisTaskRhoBaseDev(name, histo);
  rhotask->SetOutRhoName(nRho);

  AliParticleContainer* partCont = rhotask->AddParticleContainer(trackName.Data());
  AliClusterContainer *clusterCont = rhotask->AddClusterContainer(clusName.Data());
  if (clusterCont) {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(0.3);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  AliJetContainer *jetCont = new AliJetContainer(jetType, AliJetContainer::kt_algorithm, rscheme, jetradius, partCont, clusterCont);
  if (jetCont) {
    jetCont->SetJetPtCut(0);
    jetCont->SetJetAcceptanceType(acceptance);
    jetCont->SetName("Background");
    rhotask->AdoptJetContainer(jetCont);
  }

  // Final settings, pass to manager and set the containers
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
