// Track QA task (efficiency and pt resolution)
//
// Author: S.Aiola

#include <TH3F.h>
#include <THnSparse.h>
#include <TMath.h>
#include <TString.h>
#include <Riostream.h>

#include "AliPicoTrack.h"
#include "AliAODMCParticle.h"
#include "AliParticleContainer.h"
#include "AliLog.h"

#include "AliEmcalTrackingQATask.h"

ClassImp(AliEmcalTrackingQATask)

//________________________________________________________________________
AliEmcalTrackingQATask::AliEmcalTrackingQATask() : 
  AliAnalysisTaskEmcal("AliEmcalTrackingQA", kTRUE),
  fGeneratorLevel(0),
  fDetectorLevel(0),
  fTracksAll(0),
  fTracksSelected(0),
  fParticlesAllPhysPrim(0),
  fParticlesSelected(0),
  fParticlesFindable(0),
  fParticlesMatched(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliEmcalTrackingQATask::AliEmcalTrackingQATask(const char *name) : 
  AliAnalysisTaskEmcal(name, kTRUE),
  fGeneratorLevel(0),
  fDetectorLevel(0),
  fTracksAll(0),
  fTracksSelected(0),
  fParticlesAllPhysPrim(0),
  fParticlesSelected(0),
  fParticlesFindable(0),
  fParticlesMatched(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliEmcalTrackingQATask::~AliEmcalTrackingQATask()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalTrackingQATask::UserCreateOutputObjects()
{
  // Create my user objects.

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if (!fCreateHisto) return;

  fTracksAll = new TH3**[fNcentBins];
  fTracksSelected = new TH3**[fNcentBins];
  fParticlesAllPhysPrim = new TH3*[fNcentBins];
  fParticlesSelected = new TH3*[fNcentBins];  

  TString histname;

  for (Int_t i = 0; i < fNcentBins; i++) {

    fTracksAll[i] = new TH3*[3];
    fTracksSelected[i] = new TH3*[3];
    for (Int_t j = 0; j < 3; j++) {
      histname = Form("fTracksAll_%d_%d",i,j);
      fTracksAll[i][j] = new TH3F(histname,histname, 100, -1, 1, 101, 0, TMath::Pi() * 2.02, fNbins, fMinBinPt, fMaxBinPt);
      fTracksAll[i][j]->GetXaxis()->SetTitle("#eta");
      fTracksAll[i][j]->GetYaxis()->SetTitle("#phi");
      fTracksAll[i][j]->GetZaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fOutput->Add(fTracksAll[i][j]);

      histname = Form("fTracksSelected_%d_%d",i,j);
      fTracksSelected[i][j] = new TH3F(histname,histname, 100, -1, 1, 101, 0, TMath::Pi() * 2.02, fNbins, fMinBinPt, fMaxBinPt);
      fTracksSelected[i][j]->GetXaxis()->SetTitle("#eta");
      fTracksSelected[i][j]->GetYaxis()->SetTitle("#phi");
      fTracksSelected[i][j]->GetZaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fOutput->Add(fTracksSelected[i][j]);
    }

    histname = Form("fParticlesAllPhysPrim_%d",i);
    fParticlesAllPhysPrim[i] = new TH3F(histname,histname, 100, -1, 1, 101, 0, TMath::Pi() * 2.02, fNbins, fMinBinPt, fMaxBinPt);
    fParticlesAllPhysPrim[i]->GetXaxis()->SetTitle("#eta");
    fParticlesAllPhysPrim[i]->GetYaxis()->SetTitle("#phi");
    fParticlesAllPhysPrim[i]->GetZaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutput->Add(fParticlesAllPhysPrim[i]);

    histname = Form("fParticlesSelected_%d",i);
    fParticlesSelected[i] = new TH3F(histname,histname, 100, -1, 1, 101, 0, TMath::Pi() * 2.02, fNbins, fMinBinPt, fMaxBinPt);
    fParticlesSelected[i]->GetXaxis()->SetTitle("#eta");
    fParticlesSelected[i]->GetYaxis()->SetTitle("#phi");
    fParticlesSelected[i]->GetZaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutput->Add(fParticlesSelected[i]);
  }

  AllocateFindableParticlesTHnSparse();
  AllocateMatchedParticlesTHnSparse();

  PostData(1, fOutput);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::AllocateFindableParticlesTHnSparse()
{
  Int_t dim = 0;
  TString title[20];
  Int_t nbins[20] = {0};
  Double_t min[20] = {0};
  Double_t max[20] = {0};
  
  if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
    title[dim] = "Centrality %";
    nbins[dim] = 101;
    min[dim] = 0;
    max[dim] = 101;
    dim++;
  }

  title[dim] = "#it{p}_{T} (GeV/#it{c})";
  nbins[dim] = fNbins;
  min[dim] = fMinBinPt;
  max[dim] = fMaxBinPt;
  dim++;

  title[dim] = "#eta";
  nbins[dim] = 100;
  min[dim] = -1;
  max[dim] = 1;
  dim++;

  title[dim] = "#phi";
  nbins[dim] = 101;
  min[dim] = 0;
  max[dim] = TMath::Pi() * 2.02;
  dim++;
 
  fParticlesFindable = new THnSparseF("fParticlesFindable","fParticlesFindable",dim,nbins,min,max);
  for (Int_t i = 0; i < dim; i++)
    fParticlesFindable->GetAxis(i)->SetTitle(title[i]);
  fOutput->Add(fParticlesFindable);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::AllocateMatchedParticlesTHnSparse()
{
  Int_t dim = 0;
  TString title[20];
  Int_t nbins[20] = {0};
  Double_t min[20] = {0};
  Double_t max[20] = {0};
  
  if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
    title[dim] = "Centrality %";
    nbins[dim] = 101;
    min[dim] = 0;
    max[dim] = 101;
    dim++;
  }

  title[dim] = "#it{p}_{T}^{gen} (GeV/#it{c})";
  nbins[dim] = fNbins;
  min[dim] = fMinBinPt;
  max[dim] = fMaxBinPt;
  dim++;

  title[dim] = "#eta^{gen}";
  nbins[dim] = 100;
  min[dim] = -1;
  max[dim] = 1;
  dim++;

  title[dim] = "#phi^{gen}";
  nbins[dim] = 101;
  min[dim] = 0;
  max[dim] = TMath::Pi() * 2.02;
  dim++;

  title[dim] = "#it{p}_{T}^{det} (GeV/#it{c})";
  nbins[dim] = fNbins;
  min[dim] = fMinBinPt;
  max[dim] = fMaxBinPt;
  dim++;

  title[dim] = "#eta^{det}";
  nbins[dim] = 100;
  min[dim] = -1;
  max[dim] = 1;
  dim++;

  title[dim] = "#phi^{det}";
  nbins[dim] = 101;
  min[dim] = 0;
  max[dim] = TMath::Pi() * 2.02;
  dim++;

  title[dim] = "(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{gen}";
  nbins[dim] = fNbins;
  min[dim] = -1;
  max[dim] = 1;
  dim++;

  title[dim] = "track type";
  nbins[dim] = 3;
  min[dim] = -0.5;
  max[dim] = 2.5;
  dim++;

  fParticlesMatched = new THnSparseF("fParticlesMatched","fParticlesMatched",dim,nbins,min,max);
  for (Int_t i = 0; i < dim; i++)
    fParticlesMatched->GetAxis(i)->SetTitle(title[i]);
  fOutput->Add(fParticlesMatched);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::SetGeneratorLevelName(const char* name)
{
  if (!fGeneratorLevel) {  // first check if the generator level array is set
    fGeneratorLevel = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
    if (fGeneratorLevel) {  // now check if the first collection array has been added already
      fGeneratorLevel->SetArrayName(name);
    }
    else {
      fGeneratorLevel = AddParticleContainer(name);
    }
    fGeneratorLevel->SetClassName("AliAODMCParticle");
    fGeneratorLevel->SelectPhysicalPrimaries(kTRUE);
  }
  fGeneratorLevel->SetArrayName(name);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::SetDetectorLevelName(const char* name)
{
  if (!fGeneratorLevel) {
    AliError("Please, first set the generatol level array!");
    return;
  }
  if (!fDetectorLevel) {  // first check if the detector level array is set
    fDetectorLevel = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
    if (fDetectorLevel) {  // now check if the second collection array has been added already
      fDetectorLevel->SetArrayName(name);
    }
    else {
      fDetectorLevel = AddParticleContainer(name);
    }
    fDetectorLevel->SetClassName("AliPicoTrack");
    fDetectorLevel->SelectPhysicalPrimaries(kTRUE);
  }
  fDetectorLevel->SetArrayName(name);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::ExecOnce()
{
  // Init the analysis.

  if (fParticleCollArray.GetEntriesFast() < 2) {
    AliFatal("This task needs at least two particle containers!");
  }

  if (!fGeneratorLevel) {
    fGeneratorLevel = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
    fGeneratorLevel->SetClassName("AliAODMCParticle");
  }

  if (!fDetectorLevel) {
    fDetectorLevel = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
    fDetectorLevel->SetClassName("AliPicoTrack");
  }

  AliAnalysisTaskEmcal::ExecOnce();
}

//________________________________________________________________________
void AliEmcalTrackingQATask::FillFindableParticlesTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt)
{
  Double_t contents[20]={0};

  for (Int_t i = 0; i < fParticlesFindable->GetNdimensions(); i++) {
    TString title(fParticlesFindable->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = cent;
    else if (title=="#it{p}_{T} (GeV/#it{c})")
      contents[i] = partPt;
    else if (title=="#eta")
      contents[i] = partEta;
    else if (title=="#phi")
      contents[i] = partPhi;
    else 
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), fParticlesFindable->GetName()));
  }

  fParticlesFindable->Fill(contents);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::FillMatchedParticlesTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt,
							   Double_t trackEta, Double_t trackPhi, Double_t trackPt, Byte_t trackType)
{
  Double_t contents[20]={0};

  for (Int_t i = 0; i < fParticlesMatched->GetNdimensions(); i++) {
    TString title(fParticlesMatched->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = cent;
    else if (title=="#it{p}_{T}^{gen} (GeV/#it{c})")
      contents[i] = partPt;
    else if (title=="#eta^{gen}")
      contents[i] = partEta;
    else if (title=="#phi^{gen}")
      contents[i] = partPhi;
    else if (title=="#it{p}_{T}^{det} (GeV/#it{c})")
      contents[i] = trackPt;
    else if (title=="#eta^{det}")
      contents[i] = trackEta;
    else if (title=="#phi^{det}")
      contents[i] = trackPhi;
    else if (title=="(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{gen}")
      contents[i] = (partPt - trackPt) / partPt;
    else if (title=="track type")
      contents[i] = trackType;
    else 
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), fParticlesMatched->GetName()));
  }

  fParticlesMatched->Fill(contents);
}

//________________________________________________________________________
Bool_t AliEmcalTrackingQATask::FillHistograms()
{
  // Fill the histograms.

  AliPicoTrack *track = static_cast<AliPicoTrack*>(fDetectorLevel->GetNextAcceptParticle(0));
  while (track != 0) {
    Byte_t type = track->GetTrackType();
    if (type<= 2 && type >= 0) {
      fTracksAll[fCentBin][type]->Fill(track->Eta(), track->Phi(), track->Pt());

      Int_t label = TMath::Abs(track->GetLabel());

      if (label==0 || track->GetGeneratorIndex() == 0) {  // reject particles generated from other generator in the cocktail but keep fake tracks (label == 0)
	fTracksSelected[fCentBin][type]->Fill(track->Eta(), track->Phi(), track->Pt());
      }
      
      if (label > 0) {
	AliAODMCParticle *part =  static_cast<AliAODMCParticle*>(fGeneratorLevel->GetAcceptParticleWithLabel(label));
	if (part) {
	  if (part->GetGeneratorIndex() == 0) {
	    Int_t pdg = TMath::Abs(part->PdgCode());
	    // select charged pions, protons, kaons , electrons, muons
	    if (pdg == 211 || pdg == 2212 || pdg == 321 || pdg == 11 || pdg == 13) {
	      FillMatchedParticlesTHnSparse(fCent, part->Eta(), part->Phi(), part->Pt(), track->Eta(), track->Phi(), track->Pt(), type);
	    }
	  }
	}
      }
    }
    else {
      AliError(Form("Track %d has type %d not recognized!", fDetectorLevel->GetCurrentID(), type));
    }

    track = static_cast<AliPicoTrack*>(fDetectorLevel->GetNextAcceptParticle());
  }

  AliAODMCParticle *part = static_cast<AliAODMCParticle*>(fGeneratorLevel->GetNextAcceptParticle(0));
  while (part != 0) {
    fParticlesAllPhysPrim[fCentBin]->Fill(part->Eta(), part->Phi(), part->Pt());

    if (part->GetGeneratorIndex() == 0) {
      fParticlesSelected[fCentBin]->Fill(part->Eta(), part->Phi(), part->Pt());

      Int_t pdg = TMath::Abs(part->PdgCode());
      // select charged pions, protons, kaons , electrons, muons
      if (pdg == 211 || pdg == 2212 || pdg == 321 || pdg == 11 || pdg == 13) {
	FillFindableParticlesTHnSparse(fCent, part->Eta(), part->Phi(), part->Pt());
      }
    }
    
    part = static_cast<AliAODMCParticle*>(fGeneratorLevel->GetNextAcceptParticle());
  }

  return kTRUE;
}
