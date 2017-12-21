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
  AliAnalysisTaskJetUE(),
  fOutRhoName(),
  fOutRhoScaledName(),
  fRhoFunction(0),
  fScaleFunction(0),
  fAttachToEvent(kTRUE),
  fTaskConfigured(kFALSE),
  fOutRho(0),
  fOutRhoScaled(0),
  fHistRhoVsCent(0),
  fHistRhoVsLeadJetPt(),
  fHistLeadJetPtVsCent(),
  fHistLeadJetPtDensityVsCent(),
  fHistTotJetAreaVsCent(),
  fHistLeadJetNconstVsCent(),
  fHistLeadJetNconstVsPt(),
  fHistNjetVsCent(),
  fHistNjetVsNtrack(),
  fHistRhoVsLeadTrackPt(0),
  fHistRhoVsNtrack(0),
  fHistLeadTrackPtVsCent(0),
  fHistNtrackVsCent(0),
  fHistRhoVsLeadClusterE(0),
  fHistRhoVsNcluster(0),
  fHistLeadClusterEVsCent(0),
  fHistNclusterVsCent(0),
  fHistRhoScaledVsCent(0),
  fHistRhoScaledVsNtrack(0),
  fHistRhoScaledVsNcluster(0)
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name  Name of the task
 * @param[in] histo If kTRUE, the task will also produce QA histograms
 */
AliAnalysisTaskRhoBaseDev::AliAnalysisTaskRhoBaseDev(const char *name, Bool_t histo) :
  AliAnalysisTaskJetUE(name, histo),
  fOutRhoName(),
  fOutRhoScaledName(),
  fRhoFunction(0),
  fScaleFunction(0),
  fAttachToEvent(kTRUE),
  fTaskConfigured(kFALSE),
  fOutRho(0),
  fOutRhoScaled(0),
  fHistRhoVsCent(0),
  fHistRhoVsLeadJetPt(),
  fHistLeadJetPtVsCent(),
  fHistLeadJetPtDensityVsCent(),
  fHistTotJetAreaVsCent(),
  fHistLeadJetNconstVsCent(),
  fHistLeadJetNconstVsPt(),
  fHistNjetVsCent(),
  fHistNjetVsNtrack(),
  fHistRhoVsLeadTrackPt(0),
  fHistRhoVsNtrack(0),
  fHistLeadTrackPtVsCent(0),
  fHistNtrackVsCent(0),
  fHistRhoVsLeadClusterE(0),
  fHistRhoVsNcluster(0),
  fHistLeadClusterEVsCent(0),
  fHistNclusterVsCent(0),
  fHistRhoScaledVsCent(0),
  fHistRhoScaledVsNtrack(0),
  fHistRhoScaledVsNcluster(0)
{
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

  ::Info("UserCreateOutputObjects", "CreateOutputObjects of task %s", GetName());

  TString name;

  Int_t maxTracks = 6000;
  Double_t maxRho = 500;
  Int_t nRhoBins = 500;

  if (fForceBeamType == kpp) {
    maxRho = 50;
    maxTracks = 200;
  }
  else if (fForceBeamType == kpA) {
    maxRho = 200;
    maxTracks = 500;
  }

  Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);
  
  fHistRhoVsCent = new TH2F("fHistRhoVsCent", "fHistRhoVsCent", 100, 0,  100, nRhoBins, 0, maxRho);
  fHistRhoVsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhoVsCent->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
  fOutput->Add(fHistRhoVsCent);

  if (fParticleCollArray.size() > 0) {
    fHistRhoVsNtrack = new TH2F("fHistRhoVsNtrack", "fHistRhoVsNtrack", 200, 0, maxTracks, nRhoBins, 0, maxRho);
    fHistRhoVsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistRhoVsNtrack->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoVsNtrack);

    fHistNtrackVsCent = new TH2F("fHistNtrackVsCent", "fHistNtrackVsCent", 100, 0,  100, 200, 0, maxTracks);
    fHistNtrackVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistNtrackVsCent->GetYaxis()->SetTitle("No. of tracks");
    fOutput->Add(fHistNtrackVsCent);

    fHistRhoVsLeadTrackPt = new TH2F("fHistRhoVsLeadTrackPt", "fHistRhoVsLeadTrackPt", nPtBins, 0, fMaxPt, nRhoBins, 0, maxRho);
    fHistRhoVsLeadTrackPt->GetXaxis()->SetTitle("#it{p}_{T,track} (GeV/c)");
    fHistRhoVsLeadTrackPt->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoVsLeadTrackPt);

    fHistLeadTrackPtVsCent = new TH2F("fHistLeadTrackPtVsCent", "fHistLeadTrackPtVsCent", 100, 0,  100, nPtBins, 0, fMaxPt);
    fHistLeadTrackPtVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistLeadTrackPtVsCent->GetYaxis()->SetTitle("#it{p}_{T,track} (GeV/c)");
    fOutput->Add(fHistLeadTrackPtVsCent);
  }

  if (fClusterCollArray.size()>0) {
    fHistRhoVsNcluster = new TH2F("fHistRhoVsNcluster", "fHistRhoVsNcluster", 50, 0, maxTracks / 4, nRhoBins, 0, maxRho);
    fHistRhoVsNcluster->GetXaxis()->SetTitle("No. of clusters");
    fHistRhoVsNcluster->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoVsNcluster);

    fHistNclusterVsCent = new TH2F("fHistNclusterVsCent", "fHistNclusterVsCent", 100, 0,  100, 50, 0, maxTracks / 4);
    fHistNclusterVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistNclusterVsCent->GetYaxis()->SetTitle("No. of clusters");
    fOutput->Add(fHistNclusterVsCent);

    fHistRhoVsLeadClusterE = new TH2F("fHistRhoVsLeadClusterE", "fHistRhoVsLeadClusterE", nPtBins, 0, fMaxPt, nRhoBins, 0, maxRho);
    fHistRhoVsLeadClusterE->GetXaxis()->SetTitle("#it{p}_{T,track} (GeV/c)");
    fHistRhoVsLeadClusterE->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoVsLeadClusterE);

    fHistLeadClusterEVsCent = new TH2F("fHistLeadClusterEVsCent", "fHistLeadClusterEVsCent", 100, 0,  100, nPtBins, 0, fMaxPt);
    fHistLeadClusterEVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistLeadClusterEVsCent->GetYaxis()->SetTitle("#it{p}_{T,track} (GeV/c)");
    fOutput->Add(fHistLeadClusterEVsCent);
  }

  for (auto jetCont : fJetCollArray) {
    name = TString::Format("%s_fHistRhoVsLeadJetPt", jetCont.first.c_str());
    fHistRhoVsLeadJetPt[jetCont.first] = new TH2F(name, name, nPtBins, 0, fMaxPt, nRhoBins, 0, maxRho);
    fHistRhoVsLeadJetPt[jetCont.first]->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
    fHistRhoVsLeadJetPt[jetCont.first]->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoVsLeadJetPt[jetCont.first]);

    name = TString::Format("%s_fHistLeadJetPtVsCent", jetCont.first.c_str());
    fHistLeadJetPtVsCent[jetCont.first] = new TH2F(name, name, 100, 0,  100, nPtBins, 0, fMaxPt);
    fHistLeadJetPtVsCent[jetCont.first]->GetXaxis()->SetTitle("Centrality (%)");
    fHistLeadJetPtVsCent[jetCont.first]->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
    fOutput->Add(fHistLeadJetPtVsCent[jetCont.first]);

    name = TString::Format("%s_fHistLeadJetPtDensityVsCent", jetCont.first.c_str());
    fHistLeadJetPtDensityVsCent[jetCont.first] = new TH2F(name, name, 100, 0,  100, nPtBins, 0, fMaxPt*2);
    fHistLeadJetPtDensityVsCent[jetCont.first]->GetXaxis()->SetTitle("Centrality (%)");
    fHistLeadJetPtDensityVsCent[jetCont.first]->GetYaxis()->SetTitle("#it{p}_{T,jet} / #it{A}_{jet} (GeV/#it{c})");
    fOutput->Add(fHistLeadJetPtDensityVsCent[jetCont.first]);

    name = TString::Format("%s_fHistLeadJetNconstVsCent", jetCont.first.c_str());
    fHistLeadJetNconstVsCent[jetCont.first] = new TH2F(name, name, 100, 0,  100, 150, -0.5, 149.5);
    fHistLeadJetNconstVsCent[jetCont.first]->GetXaxis()->SetTitle("Centrality (%)");
    fHistLeadJetNconstVsCent[jetCont.first]->GetYaxis()->SetTitle("No. of constituents");
    fOutput->Add(fHistLeadJetNconstVsCent[jetCont.first]);

    if (fCentBins.size() > 1) {
      fHistLeadJetNconstVsPt[jetCont.first] = new TH2*[fCentBins.size()-1];
      for (Int_t i = 0; i < fCentBins.size()-1; i++) {
        name = TString::Format("%s_fHistJetNconstVsPt_Cent%d_%d", jetCont.first.c_str(), TMath::FloorNint(fCentBins[i]), TMath::FloorNint(fCentBins[i+1]));
        fHistLeadJetNconstVsPt[jetCont.first][i] = new TH2F(name, name, nPtBins, 0, fMaxPt, 150, -0.5, 149.5);
        fHistLeadJetNconstVsPt[jetCont.first][i]->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/#it{c})");
        fHistLeadJetNconstVsPt[jetCont.first][i]->GetYaxis()->SetTitle("No. of constituents");
        fOutput->Add(fHistLeadJetNconstVsPt[jetCont.first][i]);
      }
    }

    name = TString::Format("%s_fHistTotJetAreaVsCent", jetCont.first.c_str());
    fHistTotJetAreaVsCent[jetCont.first] = new TH2F(name, name, 100, 0, 100, 500, 0, 15);
    fHistTotJetAreaVsCent[jetCont.first]->GetXaxis()->SetTitle("Centrality (%)");
    fHistTotJetAreaVsCent[jetCont.first]->GetYaxis()->SetTitle("Jet area");
    fOutput->Add(fHistTotJetAreaVsCent[jetCont.first]);

    name = TString::Format("%s_fHistNjetVsCent", jetCont.first.c_str());
    fHistNjetVsCent[jetCont.first] = new TH2F(name, name, 100, 0, 100, 150, -0.5, 149.5);
    fHistNjetVsCent[jetCont.first]->GetXaxis()->SetTitle("Centrality (%)");
    fHistNjetVsCent[jetCont.first]->GetYaxis()->SetTitle("No. of jets");
    fOutput->Add(fHistNjetVsCent[jetCont.first]);

    if (fParticleCollArray.size() > 0) {
      name = TString::Format("%s_fHistNjetVsNtrack", jetCont.first.c_str());
      fHistNjetVsNtrack[jetCont.first] = new TH2F(name, name, 200, 0, maxTracks, 150, -0.5, 149.5);
      fHistNjetVsNtrack[jetCont.first]->GetXaxis()->SetTitle("No. of tracks");
      fHistNjetVsNtrack[jetCont.first]->GetYaxis()->SetTitle("No. of jets");
      fOutput->Add(fHistNjetVsNtrack[jetCont.first]);
    }
  }

  if (fScaleFunction) {
    fHistRhoScaledVsCent = new TH2F("fHistRhoScaledVsCent", "fHistRhoScaledVsCent", 100, 0, 100, nRhoBins, 0, maxRho);
    fHistRhoScaledVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistRhoScaledVsCent->GetYaxis()->SetTitle("#rho_{scaled} (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoScaledVsCent);

    if (fParticleCollArray.size() > 0) {
      fHistRhoScaledVsNtrack = new TH2F("fHistRhoScaledVsNtrack", "fHistRhoScaledVsNtrack", 200, 0, maxTracks, nRhoBins, 0, maxRho);
      fHistRhoScaledVsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistRhoScaledVsNtrack->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
      fOutput->Add(fHistRhoScaledVsNtrack);
    }

    if (fClusterCollArray.size() > 0) {
      fHistRhoScaledVsNcluster = new TH2F("fHistRhoScaledVsNcluster", "fHistRhoScaledVsNcluster", 50, 0, maxTracks / 4, nRhoBins, 0, maxRho);
      fHistRhoScaledVsNcluster->GetXaxis()->SetTitle("No. of clusters");
      fHistRhoScaledVsNcluster->GetYaxis()->SetTitle("#rho_{scaled} (GeV/#it{c} #times rad^{-1})");
      fOutput->Add(fHistRhoScaledVsNcluster);
    }
  }
}


/**
 * Calculates the average background using a given parametrization
 * as a function of centrality.
 */
void AliAnalysisTaskRhoBaseDev::CalculateRho()
{
  Double_t rho = GetRhoFactor(fCent);
  fOutRho->SetVal(rho);
}

/**
 * Run the analysis.
 */
Bool_t AliAnalysisTaskRhoBaseDev::Run()
{
  fOutRho->SetVal(0);
  if (fOutRhoScaled) fOutRhoScaled->SetVal(0);

  CalculateEventProperties();

  CalculateRho();

  if (fScaleFunction) {
    Double_t rhoScaled = fOutRho->GetVal() * GetScaleFactor(fCent);
    fOutRhoScaled->SetVal(rhoScaled);
  }

  return kTRUE;
}

/**
 * Fill histograms.
 */
Bool_t AliAnalysisTaskRhoBaseDev::FillHistograms()
{
  fHistRhoVsCent->Fill(fCent, fOutRho->GetVal());

  if (fLeadingParticle) {
    fHistLeadTrackPtVsCent->Fill(fCent, fLeadingParticle->Pt());
    fHistRhoVsLeadTrackPt->Fill(fLeadingParticle->Pt(), fOutRho->GetVal());
  }

  if (fLeadingCluster) {
    fHistLeadClusterEVsCent->Fill(fCent, fLeadingCluster->E());
    fHistRhoVsLeadClusterE->Fill(fLeadingCluster->E(), fOutRho->GetVal());
  }

  if (fHistNtrackVsCent) fHistNtrackVsCent->Fill(fCent, fNtracks);
  if (fHistNclusterVsCent) fHistNclusterVsCent->Fill(fCent, fNclusters);


  if (fHistRhoVsNtrack) fHistRhoVsNtrack->Fill(fNtracks, fOutRho->GetVal());
  if (fHistRhoVsNcluster) fHistRhoVsNcluster->Fill(fNclusters, fOutRho->GetVal());

  if (fOutRhoScaled) {
    fHistRhoScaledVsCent->Fill(fCent, fOutRhoScaled->GetVal());
    if (fHistRhoScaledVsNtrack) fHistRhoScaledVsNtrack->Fill(fNtracks, fOutRhoScaled->GetVal());
    if (fHistRhoScaledVsNcluster) fHistRhoScaledVsNcluster->Fill(fNclusters,  fOutRhoScaled->GetVal());
  }

  for (auto jetCont : fJetCollArray) {
    fHistTotJetAreaVsCent[jetCont.first]->Fill(fCent, fTotJetArea[jetCont.first]);
    if (fLeadingJet[jetCont.first]) {
      fHistLeadJetPtVsCent[jetCont.first]->Fill(fCent, fLeadingJet[jetCont.first]->Pt());
      fHistLeadJetPtDensityVsCent[jetCont.first]->Fill(fCent, fLeadingJet[jetCont.first]->Pt() / fLeadingJet[jetCont.first]->Area());
      fHistLeadJetNconstVsCent[jetCont.first]->Fill(fCent, fLeadingJet[jetCont.first]->GetNumberOfConstituents());
      fHistRhoVsLeadJetPt[jetCont.first]->Fill(fLeadingJet[jetCont.first]->Pt(), fOutRho->GetVal());
      if (fCentBin >=0 && fCentBin < fCentBins.size()-1 && fHistLeadJetNconstVsPt[jetCont.first]) fHistLeadJetNconstVsPt[jetCont.first][fCentBin]->Fill(fLeadingJet[jetCont.first]->Pt(), fLeadingJet[jetCont.first]->GetNumberOfConstituents());
    }
    if (fHistNjetVsCent[jetCont.first]) fHistNjetVsCent[jetCont.first]->Fill(fCent, fNjets[jetCont.first]);

    if (fHistNjetVsNtrack[jetCont.first]) fHistNjetVsNtrack[jetCont.first]->Fill(fNtracks, fNjets[jetCont.first]);
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

  fTaskConfigured = VerifyContainers();

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
