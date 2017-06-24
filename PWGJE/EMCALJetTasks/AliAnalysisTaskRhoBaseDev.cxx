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
  AliAnalysisTaskEmcalJetLight(),
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
  fHistLeadJetPtVsCent(0),
  fHistLeadTrackPtVsCent(0),
  fHistLeadClusterEVsCent(0),
  fHistLeadJetPtDensityVsCent(0),
  fHistSubLeadJetPtVsCent(0),
  fHistTotJetAreaVsCent(0),
  fHistLeadJetNconstVsCent(0),
  fHistLeadJetNconstVsPt(0),
  fHistNtrackVsCent(0),
  fHistNclusterVsCent(0),
  fHistNjetVsCent(0),
  fHistNjetVsNtrack(0),
  fHistRhoVsCent(0),
  fHistRhoVsLeadJetPt(0),
  fHistRhoVsLeadTrackPt(0),
  fHistRhoVsLeadClusterE(0),
  fHistRhoVsNtrack(0),
  fHistRhoVsNcluster(0),
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
  fHistLeadJetPtVsCent(0),
  fHistLeadTrackPtVsCent(0),
  fHistLeadClusterEVsCent(0),
  fHistLeadJetPtDensityVsCent(0),
  fHistSubLeadJetPtVsCent(0),
  fHistTotJetAreaVsCent(0),
  fHistLeadJetNconstVsCent(0),
  fHistLeadJetNconstVsPt(0),
  fHistNtrackVsCent(0),
  fHistNclusterVsCent(0),
  fHistNjetVsCent(0),
  fHistNjetVsNtrack(0),
  fHistRhoVsCent(0),
  fHistRhoVsLeadJetPt(0),
  fHistRhoVsLeadTrackPt(0),
  fHistRhoVsLeadClusterE(0),
  fHistRhoVsNtrack(0),
  fHistRhoVsNcluster(0),
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

  ::Info("UserCreateOutputObjects", "CreateOutputObjects of task %s", GetName());

  AliAnalysisTaskEmcalJetLight::UserCreateOutputObjects();
  //ranges for PbPb
  Float_t Ntrackrange[2] = {0, 6000};
  //set multiplicity related axes to a smaller max value
  if (fBeamType != kAA) Ntrackrange[1] = 200.;
  
  fHistRhoVsCent = new TH2F("fHistRhoVsCent", "fHistRhoVsCent", 100, 0,  100, fNbins, fMinBinPt, fMaxBinPt);
  fHistRhoVsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhoVsCent->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
  fOutput->Add(fHistRhoVsCent);

  if (fParticleCollArray.size() > 0) {
    fHistRhoVsNtrack = new TH2F("fHistRhoVsNtrack", "fHistRhoVsNtrack", 200, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt);
    fHistRhoVsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistRhoVsNtrack->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoVsNtrack);

    fHistNtrackVsCent = new TH2F("fHistNtrackVsCent", "fHistNtrackVsCent", 100, 0,  100, 200, Ntrackrange[0], Ntrackrange[1]);
    fHistNtrackVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistNtrackVsCent->GetYaxis()->SetTitle("No. of tracks");
    fOutput->Add(fHistNtrackVsCent);

    fHistRhoVsLeadTrackPt = new TH2F("fHistRhoVsLeadTrackPt", "fHistRhoVsLeadTrackPt", fNbins, fMinBinPt, fMaxBinPt*2, fNbins, fMinBinPt, fMaxBinPt);
    fHistRhoVsLeadTrackPt->GetXaxis()->SetTitle("#it{p}_{T,track} (GeV/c)");
    fHistRhoVsLeadTrackPt->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoVsLeadTrackPt);

    fHistLeadTrackPtVsCent = new TH2F("fHistLeadTrackPtVsCent", "fHistLeadTrackPtVsCent", 100, 0,  100, fNbins, fMinBinPt, fMaxBinPt*2);
    fHistLeadTrackPtVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistLeadTrackPtVsCent->GetYaxis()->SetTitle("#it{p}_{T,track} (GeV/c)");
    fOutput->Add(fHistLeadTrackPtVsCent);
  }

  if (fClusterCollArray.size()>0) {
    fHistRhoVsNcluster = new TH2F("fHistRhoVsNcluster", "fHistRhoVsNcluster", 50, Ntrackrange[0] / 4, Ntrackrange[1] / 4, fNbins, fMinBinPt, fMaxBinPt);
    fHistRhoVsNcluster->GetXaxis()->SetTitle("No. of clusters");
    fHistRhoVsNcluster->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoVsNcluster);

    fHistNclusterVsCent = new TH2F("fHistNclusterVsCent", "fHistNclusterVsCent", 100, 0,  100, 50, Ntrackrange[0] / 4, Ntrackrange[1] / 4);
    fHistNclusterVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistNclusterVsCent->GetYaxis()->SetTitle("No. of clusters");
    fOutput->Add(fHistNclusterVsCent);

    fHistRhoVsLeadClusterE = new TH2F("fHistRhoVsLeadClusterE", "fHistRhoVsLeadClusterE", fNbins, fMinBinPt, fMaxBinPt*2, fNbins, fMinBinPt, fMaxBinPt);
    fHistRhoVsLeadClusterE->GetXaxis()->SetTitle("#it{p}_{T,track} (GeV/c)");
    fHistRhoVsLeadClusterE->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoVsLeadClusterE);

    fHistLeadClusterEVsCent = new TH2F("fHistLeadClusterEVsCent", "fHistLeadClusterEVsCent", 100, 0,  100, fNbins, fMinBinPt, fMaxBinPt*2);
    fHistLeadClusterEVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistLeadClusterEVsCent->GetYaxis()->SetTitle("#it{p}_{T,track} (GeV/c)");
    fOutput->Add(fHistLeadClusterEVsCent);
  }

  if (fJetCollArray.size() > 0) {
    fHistRhoVsLeadJetPt = new TH2F("fHistRhoVsLeadJetPt", "fHistRhoVsLeadJetPt", fNbins, fMinBinPt, fMaxBinPt*2, fNbins, fMinBinPt, fMaxBinPt);
    fHistRhoVsLeadJetPt->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
    fHistRhoVsLeadJetPt->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoVsLeadJetPt);

    fHistLeadJetPtVsCent = new TH2F("fHistLeadJetPtVsCent", "fHistLeadJetPtVsCent", 100, 0,  100, fNbins, fMinBinPt, fMaxBinPt*2);
    fHistLeadJetPtVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistLeadJetPtVsCent->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
    fOutput->Add(fHistLeadJetPtVsCent);

    fHistSubLeadJetPtVsCent = new TH2F("fHistSubLeadJetPtVsCent", "fHistSubLeadJetPtVsCent", 100, 0,  100, fNbins, fMinBinPt, fMaxBinPt*2);
    fHistSubLeadJetPtVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistSubLeadJetPtVsCent->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/#it{c})");
    fOutput->Add(fHistSubLeadJetPtVsCent);

    fHistLeadJetPtDensityVsCent = new TH2F("fHistLeadJetPtDensityVsCent", "fHistLeadJetPtDensityVsCent", 100, 0,  100, fNbins, fMinBinPt, fMaxBinPt*4);
    fHistLeadJetPtDensityVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistLeadJetPtDensityVsCent->GetYaxis()->SetTitle("#it{p}_{T,jet} / #it{A}_{jet} (GeV/#it{c})");
    fOutput->Add(fHistLeadJetPtDensityVsCent);

    fHistLeadJetNconstVsCent = new TH2F("fHistLeadJetNconstVsCent", "fHistLeadJetNconstVsCent", 100, 0,  100, 150, -0.5, 149.5);
    fHistLeadJetNconstVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistLeadJetNconstVsCent->GetYaxis()->SetTitle("No. of constituents");
    fOutput->Add(fHistLeadJetNconstVsCent);

    if (fCentBins.size() > 1) {
      TString name;
      fHistLeadJetNconstVsPt = new TH2*[fCentBins.size()-1];
      for (Int_t i = 0; i < fCentBins.size()-1; i++) {
        name = Form("fHistJetNconstVsPt_Cent%d_%d", TMath::FloorNint(fCentBins[i]), TMath::FloorNint(fCentBins[i+1]));
        fHistLeadJetNconstVsPt[i] = new TH2F(name, name, fNbins, fMinBinPt, fMaxBinPt*2, 150, -0.5, 149.5);
        fHistLeadJetNconstVsPt[i]->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/#it{c})");
        fHistLeadJetNconstVsPt[i]->GetYaxis()->SetTitle("No. of constituents");
        fOutput->Add(fHistLeadJetNconstVsPt[i]);
      }
    }

    fHistTotJetAreaVsCent = new TH2F("fHistTotJetAreaVsCent", "fHistTotJetAreaVsCent", 100, 0, 100, 500, 0, 15);
    fHistTotJetAreaVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistTotJetAreaVsCent->GetYaxis()->SetTitle("Jet area");
    fOutput->Add(fHistTotJetAreaVsCent);

    fHistNjetVsCent = new TH2F("fHistNjetVsCent",  "fHistNjetVsCent", 100, 0, 100, 150, -0.5, 149.5);
    fHistNjetVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistNjetVsCent->GetYaxis()->SetTitle("No. of jets");
    fOutput->Add(fHistNjetVsCent);

    if (fParticleCollArray.size() > 0) {
      fHistNjetVsNtrack = new TH2F("fHistNjetVsNtrack",  "fHistNjetVsNtrack", 200, Ntrackrange[0], Ntrackrange[1], 150, -0.5, 149.5);
      fHistNjetVsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistNjetVsNtrack->GetYaxis()->SetTitle("No. of jets");
      fOutput->Add(fHistNjetVsNtrack);
    }
  }

  if (fScaleFunction) {
    fHistRhoScaledVsCent = new TH2F("fHistRhoScaledVsCent", "fHistRhoScaledVsCent", 100, 0, 100, fNbins, fMinBinPt , fMaxBinPt);
    fHistRhoScaledVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistRhoScaledVsCent->GetYaxis()->SetTitle("#rho_{scaled} (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoScaledVsCent);

    if (fParticleCollArray.size() > 0) {
      fHistRhoScaledVsNtrack = new TH2F("fHistRhoScaledVsNtrack", "fHistRhoScaledVsNtrack", 200, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt);
      fHistRhoScaledVsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistRhoScaledVsNtrack->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
      fOutput->Add(fHistRhoScaledVsNtrack);
    }

    if (fClusterCollArray.size() > 0) {
      fHistRhoScaledVsNcluster = new TH2F("fHistRhoScaledVsNcluster", "fHistRhoScaledVsNcluster", 50, Ntrackrange[0] / 4, Ntrackrange[1] / 4, fNbins, fMinBinPt, fMaxBinPt);
      fHistRhoScaledVsNcluster->GetXaxis()->SetTitle("No. of clusters");
      fHistRhoScaledVsNcluster->GetYaxis()->SetTitle("#rho_{scaled} (GeV/#it{c} #times rad^{-1})");
      fOutput->Add(fHistRhoScaledVsNcluster);
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
  AliVParticle* leadPart = nullptr;
  for (auto partCont : fParticleCollArray) {
    for (auto track : partCont.second->accepted()) {
      if (!leadPart || track->Pt() > leadPart->Pt()) leadPart = track;
      Ntracks++;
    }
  }

  // Loop over all possible containers
  AliVCluster* leadClus = nullptr;
  for (auto clusCont : fClusterCollArray) {
    for (auto clus : clusCont.second->accepted()) {
      if (!leadClus || clus->E() > leadClus->E()) leadClus = clus;
      Nclusters++;
    }
  }

  AliEmcalJet* leadJet = nullptr;
  AliEmcalJet* subLeadJet = nullptr;
  Double_t totJetArea = 0;

  AliJetContainer* jetCont = fJetCollArray["Background"];
  for (auto jet : jetCont->accepted()) {
    if (!leadJet || jet->Pt() > leadJet->Pt()) {
      subLeadJet = leadJet;
      leadJet = jet;
    }
    else if (!subLeadJet || jet->Pt() > subLeadJet->Pt()) {
      subLeadJet = jet;
    }

    if (!jet->IsGhost()) {
      NjetAcc++;
      totJetArea += jet->Area();
    }
  }

  fHistTotJetAreaVsCent->Fill(fCent, totJetArea);
  fHistRhoVsCent->Fill(fCent, fOutRho->GetVal());

  if (leadJet) {
    fHistLeadJetPtVsCent->Fill(fCent, leadJet->Pt());
    fHistLeadJetPtDensityVsCent->Fill(fCent, leadJet->Pt() / leadJet->Area());
    fHistLeadJetNconstVsCent->Fill(fCent, leadJet->GetNumberOfConstituents());
    fHistRhoVsLeadJetPt->Fill(leadJet->Pt(), fOutRho->GetVal());
    if (fCentBin >=0 && fCentBin < fCentBins.size()-1 && fHistLeadJetNconstVsPt) fHistLeadJetNconstVsPt[fCentBin]->Fill(leadJet->Pt(), leadJet->GetNumberOfConstituents());
  }

  if (subLeadJet) {
    fHistSubLeadJetPtVsCent->Fill(fCent, subLeadJet->Pt());
  }

  if (leadPart) {
    fHistLeadTrackPtVsCent->Fill(fCent, leadPart->Pt());
    fHistRhoVsLeadTrackPt->Fill(leadPart->Pt(), fOutRho->GetVal());
  }

  if (leadClus) {
    fHistLeadClusterEVsCent->Fill(fCent, leadClus->E());
    fHistRhoVsLeadClusterE->Fill(leadClus->E(), fOutRho->GetVal());
  }

  if (fHistNtrackVsCent) fHistNtrackVsCent->Fill(fCent, Ntracks);
  if (fHistNclusterVsCent) fHistNclusterVsCent->Fill(fCent, Nclusters);
  if (fHistNjetVsCent) fHistNjetVsCent->Fill(fCent, NjetAcc);

  if (fHistNjetVsNtrack) fHistNjetVsNtrack->Fill(Ntracks, NjetAcc);

  if (fHistRhoVsNtrack) fHistRhoVsNtrack->Fill(Ntracks, fOutRho->GetVal());
  if (fHistRhoVsNcluster) fHistRhoVsNcluster->Fill(Nclusters, fOutRho->GetVal());

  if (fOutRhoScaled) {
    fHistRhoScaledVsCent->Fill(fCent, fOutRhoScaled->GetVal());
    if (fHistRhoScaledVsNtrack) fHistRhoScaledVsNtrack->Fill(Ntracks, fOutRhoScaled->GetVal());
    if (fHistRhoScaledVsNcluster) fHistRhoScaledVsNcluster->Fill(Nclusters,  fOutRhoScaled->GetVal());
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
