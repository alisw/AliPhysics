// Track QA task (efficiency and pt resolution)
//
// Author: S.Aiola

#include <TH3F.h>
#include <THnSparse.h>
#include <TMath.h>
#include <TString.h>
#include <Riostream.h>

#include "AliPicoTrack.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliMCParticleContainer.h"
#include "AliLog.h"

#include "AliEmcalTrackingQATask.h"

ClassImp(AliEmcalTrackingQATask)

//________________________________________________________________________
AliEmcalTrackingQATask::AliEmcalTrackingQATask() : 
  AliAnalysisTaskEmcal("AliEmcalTrackingQA", kTRUE),
  fDoSigma1OverPt(kFALSE),
  fDoSigmaPtOverPtGen(kFALSE),
  fGeneratorLevel(0),
  fDetectorLevel(0),
  fNPtHistBins(0),
  fPtHistBins(0),
  fNEtaHistBins(0),
  fEtaHistBins(0),
  fNPhiHistBins(0),
  fPhiHistBins(0),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtRelDiffHistBins(0),
  fPtRelDiffHistBins(0),
  fNPtResHistBins(0),
  fPtResHistBins(0),
  f1OverPtResHistBins(0),
  fN1OverPtResHistBins(0),
  fNIntegerHistBins(0),
  fIntegerHistBins(0),
  fTracks(0),
  fParticlesPhysPrim(0),
  fParticlesMatched(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);

  GenerateHistoBins();
}

//________________________________________________________________________
AliEmcalTrackingQATask::AliEmcalTrackingQATask(const char *name) : 
  AliAnalysisTaskEmcal(name, kTRUE),
  fDoSigma1OverPt(kFALSE),
  fDoSigmaPtOverPtGen(kFALSE),
  fGeneratorLevel(0),
  fDetectorLevel(0),
  fNPtHistBins(0),
  fPtHistBins(0),
  fNEtaHistBins(0),
  fEtaHistBins(0),
  fNPhiHistBins(0),
  fPhiHistBins(0),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtRelDiffHistBins(0),
  fPtRelDiffHistBins(0),
  fNPtResHistBins(0),
  fPtResHistBins(0),
  f1OverPtResHistBins(0),
  fN1OverPtResHistBins(0),
  fNIntegerHistBins(0),
  fIntegerHistBins(0),
  fTracks(0),
  fParticlesPhysPrim(0),
  fParticlesMatched(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(kTRUE);

  GenerateHistoBins();
}

//________________________________________________________________________
AliEmcalTrackingQATask::~AliEmcalTrackingQATask()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalTrackingQATask::GenerateHistoBins()
{
  fNPtHistBins = 82;
  fPtHistBins = new Double_t[fNPtHistBins+1];
  GenerateFixedBinArray(6, 0, 0.3, fPtHistBins);
  GenerateFixedBinArray(7, 0.3, 1, fPtHistBins+6);
  GenerateFixedBinArray(10, 1, 3, fPtHistBins+13);
  GenerateFixedBinArray(14, 3, 10, fPtHistBins+23);
  GenerateFixedBinArray(10, 10, 20, fPtHistBins+37);
  GenerateFixedBinArray(15, 20, 50, fPtHistBins+47);
  GenerateFixedBinArray(20, 50, 150, fPtHistBins+62);

  fNEtaHistBins = 100;
  fEtaHistBins = new Double_t[fNEtaHistBins+1];
  GenerateFixedBinArray(fNEtaHistBins, -1, 1, fEtaHistBins);

  fNPhiHistBins = 101;
  fPhiHistBins = new Double_t[fNPhiHistBins+1];
  GenerateFixedBinArray(fNPhiHistBins, 0, TMath::Pi() * 2.02, fPhiHistBins);

  fNCentHistBins = 4;
  fCentHistBins = new Double_t[fNCentHistBins+1];
  fCentHistBins[0] = 0;
  fCentHistBins[1] = 10;
  fCentHistBins[2] = 30;
  fCentHistBins[3] = 50;
  fCentHistBins[4] = 90;

  fNPtResHistBins = 175;
  fPtResHistBins = new Double_t[fNPtResHistBins+1];
  GenerateFixedBinArray(50, 0, 0.05, fPtResHistBins);
  GenerateFixedBinArray(25, 0.05, 0.10, fPtResHistBins+50);
  GenerateFixedBinArray(25, 0.10, 0.20, fPtResHistBins+75);
  GenerateFixedBinArray(30, 0.20, 0.50, fPtResHistBins+100);
  GenerateFixedBinArray(25, 0.50, 1.00, fPtResHistBins+130);
  GenerateFixedBinArray(20, 1.00, 2.00, fPtResHistBins+155);

  fNPtRelDiffHistBins = 200;
  fPtRelDiffHistBins = new Double_t[fNPtRelDiffHistBins+1];
  GenerateFixedBinArray(fNPtRelDiffHistBins, -2, 2, fPtRelDiffHistBins);

  fN1OverPtResHistBins = 385;
  f1OverPtResHistBins = new Double_t[fN1OverPtResHistBins+1];
  GenerateFixedBinArray(100, 0, 0.02, f1OverPtResHistBins);
  GenerateFixedBinArray(60, 0.02, 0.05, f1OverPtResHistBins+100);
  GenerateFixedBinArray(50, 0.05, 0.1, f1OverPtResHistBins+160);
  GenerateFixedBinArray(50, 0.1, 0.2, f1OverPtResHistBins+210);
  GenerateFixedBinArray(75, 0.2, 0.5, f1OverPtResHistBins+260);
  GenerateFixedBinArray(50, 0.5, 1.5, f1OverPtResHistBins+335);

  fNIntegerHistBins = 10;
  fIntegerHistBins = new Double_t[fNIntegerHistBins+1];
  GenerateFixedBinArray(fNIntegerHistBins, -0.5, 9.5, fIntegerHistBins);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::UserCreateOutputObjects()
{
  // Create my user objects.

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if (fParticleCollArray.GetEntriesFast() < 1) {
    AliFatal("This task needs at least one particle container!");
  }

  if (!fDetectorLevel) {
    fDetectorLevel = static_cast<AliTrackContainer*>(fParticleCollArray.At(0));
  }

  if (!fGeneratorLevel && fParticleCollArray.GetEntriesFast() > 1) {
    fGeneratorLevel = static_cast<AliMCParticleContainer*>(fParticleCollArray.At(1));
  }

  AllocateDetectorLevelTHnSparse();

  if (fGeneratorLevel) {
    AllocateGeneratorLevelTHnSparse();
    AllocateMatchedParticlesTHnSparse();
  }
}

//________________________________________________________________________
void AliEmcalTrackingQATask::AllocateDetectorLevelTHnSparse()
{
  Int_t dim = 0;
  TString title[20];
  Int_t nbins[20] = {0};
  Double_t *binEdges[20] = {0};

  if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
    title[dim] = "Centrality %";
    nbins[dim] = fNCentHistBins;
    binEdges[dim] = fCentHistBins;
    dim++;
  }

  title[dim] = "#it{p}_{T} (GeV/#it{c})";
  nbins[dim] = fNPtHistBins;
  binEdges[dim] = fPtHistBins;
  dim++;

  title[dim] = "#eta";
  nbins[dim] = fNEtaHistBins;
  binEdges[dim] = fEtaHistBins;
  dim++;

  title[dim] = "#phi";
  nbins[dim] = fNPhiHistBins;
  binEdges[dim] = fPhiHistBins;
  dim++;

  title[dim] = "MC Generator";
  nbins[dim] = 2;
  binEdges[dim] = fIntegerHistBins;
  dim++;

  title[dim] = "track type";
  nbins[dim] = 3;
  binEdges[dim] = fIntegerHistBins;
  dim++;

  if (fIsEsd) {
    if (fDoSigma1OverPt) {
      title[dim] = "#sigma(1/#it{p}_{T}) (GeV/#it{c})^{-1}";
      nbins[dim] = fN1OverPtResHistBins;
      binEdges[dim] = f1OverPtResHistBins;
      dim++;
    }
    else {
      title[dim] = "#sigma(#it{p}_{T}) / #it{p}_{T}";
      nbins[dim] = fNPtResHistBins;
      binEdges[dim] = fPtResHistBins;
      dim++;
    }    
  }

  fTracks = new THnSparseF("fTracks","fTracks",dim,nbins);
  for (Int_t i = 0; i < dim; i++) {
    fTracks->GetAxis(i)->SetTitle(title[i]);
    fTracks->SetBinEdges(i, binEdges[i]);
  }

  fOutput->Add(fTracks);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::AllocateGeneratorLevelTHnSparse()
{
  Int_t dim = 0;
  TString title[20];
  Int_t nbins[20] = {0};
  Double_t *binEdges[20] = {0};

  if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
    title[dim] = "Centrality %";
    nbins[dim] = fNCentHistBins;
    binEdges[dim] = fCentHistBins;
    dim++;
  }

  title[dim] = "#it{p}_{T} (GeV/#it{c})";
  nbins[dim] = fNPtHistBins;
  binEdges[dim] = fPtHistBins;
  dim++;

  title[dim] = "#eta";
  nbins[dim] = fNEtaHistBins;
  binEdges[dim] = fEtaHistBins;
  dim++;

  title[dim] = "#phi";
  nbins[dim] = fNPhiHistBins;
  binEdges[dim] = fPhiHistBins;
  dim++;

  title[dim] = "MC Generator";
  nbins[dim] = 2;
  binEdges[dim] = fIntegerHistBins;
  dim++;

  title[dim] = "Findable";
  nbins[dim] = 2;
  binEdges[dim] = fIntegerHistBins;
  dim++;

  fParticlesPhysPrim = new THnSparseF("fParticlesPhysPrim","fParticlesPhysPrim",dim,nbins);
  for (Int_t i = 0; i < dim; i++) {
    fParticlesPhysPrim->GetAxis(i)->SetTitle(title[i]);
    fParticlesPhysPrim->SetBinEdges(i, binEdges[i]);
  }

  fOutput->Add(fParticlesPhysPrim);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::AllocateMatchedParticlesTHnSparse()
{
  Int_t dim = 0;
  TString title[20];
  Int_t nbins[20] = {0};
  Double_t *binEdges[20] = {0};

  if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
    title[dim] = "Centrality %";
    nbins[dim] = fNCentHistBins;
    binEdges[dim] = fCentHistBins;
    dim++;
  }

  title[dim] = "#it{p}_{T}^{gen} (GeV/#it{c})";
  nbins[dim] = fNPtHistBins;
  binEdges[dim] = fPtHistBins;
  dim++;

  title[dim] = "#eta^{gen}";
  nbins[dim] = fNEtaHistBins;
  binEdges[dim] = fEtaHistBins;
  dim++;

  title[dim] = "#phi^{gen}";
  nbins[dim] = fNPhiHistBins;
  binEdges[dim] = fPhiHistBins;
  dim++;

  title[dim] = "#it{p}_{T}^{det} (GeV/#it{c})";
  nbins[dim] = fNPtHistBins;
  binEdges[dim] = fPtHistBins;
  dim++;

  title[dim] = "#eta^{det}";
  nbins[dim] = fNEtaHistBins;
  binEdges[dim] = fEtaHistBins;
  dim++;

  title[dim] = "#phi^{det}";
  nbins[dim] = fNPhiHistBins;
  binEdges[dim] = fPhiHistBins;
  dim++;

  if (fDoSigmaPtOverPtGen) {
    title[dim] = "(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{gen}";
    nbins[dim] = fNPtRelDiffHistBins;
    binEdges[dim] = fPtRelDiffHistBins;
    dim++;
  }
  else {
    title[dim] = "(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{det}";
    nbins[dim] = fNPtRelDiffHistBins;
    binEdges[dim] = fPtRelDiffHistBins;
    dim++;
  }

  title[dim] = "track type";
  nbins[dim] = 3;
  binEdges[dim] = fIntegerHistBins;
  dim++;

  fParticlesMatched = new THnSparseF("fParticlesMatched","fParticlesMatched",dim,nbins);
  for (Int_t i = 0; i < dim; i++) {
    fParticlesMatched->GetAxis(i)->SetTitle(title[i]);
    fParticlesMatched->SetBinEdges(i, binEdges[i]);
  }

  fOutput->Add(fParticlesMatched);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::SetGeneratorLevelName(const char* name)
{
  if (!fDetectorLevel) {
    AliError("Please, first set the detector level array!");
    return;
  }
  if (!fGeneratorLevel) {  // first check if the generator level array is set
    fGeneratorLevel = dynamic_cast<AliMCParticleContainer*>(fParticleCollArray.At(1));
    if (fGeneratorLevel) {  // now check if the first collection array has been added already
      fGeneratorLevel->SetArrayName(name);
    }
    else {
      fGeneratorLevel = AddMCParticleContainer(name);
    }
    fGeneratorLevel->SelectPhysicalPrimaries(kTRUE);
    fGeneratorLevel->SetParticlePtCut(0);
  }
  fGeneratorLevel->SetArrayName(name);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::SetDetectorLevelName(const char* name)
{
  if (!fDetectorLevel) {  // first check if the detector level array is set
    fDetectorLevel = static_cast<AliTrackContainer*>(fParticleCollArray.At(0));
    if (fDetectorLevel) {  // now check if the second collection array has been added already
      fDetectorLevel->SetArrayName(name);
    }
    else {
      fDetectorLevel = AddTrackContainer(name);
    }
  }
  fDetectorLevel->SetArrayName(name);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::ExecOnce()
{
  // Init the analysis.

  AliAnalysisTaskEmcal::ExecOnce();
}

//________________________________________________________________________
void AliEmcalTrackingQATask::FillDetectorLevelTHnSparse(Double_t cent, Double_t trackEta, Double_t trackPhi, Double_t trackPt, 
    Double_t sigma1OverPt, Int_t mcGen, Byte_t trackType)
{
  Double_t contents[20]={0};

  for (Int_t i = 0; i < fTracks->GetNdimensions(); i++) {
    TString title(fTracks->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = cent;
    else if (title=="#it{p}_{T} (GeV/#it{c})")
      contents[i] = trackPt;
    else if (title=="#eta")
      contents[i] = trackEta;
    else if (title=="#phi")
      contents[i] = trackPhi;
    else if (title=="#sigma(1/#it{p}_{T}) (GeV/#it{c})^{-1}")
      contents[i] = sigma1OverPt;
    else if (title=="#sigma(#it{p}_{T}) / #it{p}_{T}")
      contents[i] = sigma1OverPt*trackPt;
    else if (title=="MC Generator")
      contents[i] = mcGen;
    else if (title=="track type")
      contents[i] = trackType;
    else 
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), fTracks->GetName()));
  }

  fTracks->Fill(contents);
}

//________________________________________________________________________
void AliEmcalTrackingQATask::FillGeneratorLevelTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt, Int_t mcGen, Byte_t findable)
{
  Double_t contents[20]={0};

  for (Int_t i = 0; i < fParticlesPhysPrim->GetNdimensions(); i++) {
    TString title(fParticlesPhysPrim->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = cent;
    else if (title=="#it{p}_{T} (GeV/#it{c})")
      contents[i] = partPt;
    else if (title=="#eta")
      contents[i] = partEta;
    else if (title=="#phi")
      contents[i] = partPhi;
    else if (title=="MC Generator")
      contents[i] = mcGen;
    else if (title=="Findable")
      contents[i] = findable;
    else 
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), fParticlesPhysPrim->GetName()));
  }

  fParticlesPhysPrim->Fill(contents);
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
    else if (title=="(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{det}")
      contents[i] = (partPt - trackPt) / trackPt;
    else if (title=="track type")
      contents[i] = (Double_t)trackType;
    else 
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), fParticlesMatched->GetName()));
  }

  fParticlesMatched->Fill(contents);
}

//________________________________________________________________________
Bool_t AliEmcalTrackingQATask::FillHistograms()
{
  // Fill the histograms.

  fDetectorLevel->ResetCurrentID();
  AliVTrack *track = fDetectorLevel->GetNextAcceptTrack();
  while (track != 0) {
    Byte_t type = fDetectorLevel->GetTrackType(fDetectorLevel->GetCurrentID());
    if (type <= 2) {
      Double_t sigma = 0;
      if (fIsEsd) {
        AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track);
        if (esdTrack) sigma = TMath::Sqrt(esdTrack->GetSigma1Pt2());
      }

      Int_t label = TMath::Abs(track->GetLabel());
      Int_t mcGen = 1;
      // reject particles generated from other generators in the cocktail but keep fake tracks (label == 0)
      if (label==0 || track->GetGeneratorIndex() == 0) mcGen = 0;

      FillDetectorLevelTHnSparse(fCent, track->Eta(), track->Phi(), track->Pt(), sigma, mcGen, type);

      if (fGeneratorLevel && label > 0) {
        AliAODMCParticle *part =  fGeneratorLevel->GetAcceptMCParticleWithLabel(label);
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

    track = fDetectorLevel->GetNextAcceptTrack();
  }

  if (fGeneratorLevel) {
    fGeneratorLevel->ResetCurrentID();
    AliAODMCParticle *part = fGeneratorLevel->GetNextAcceptMCParticle();
    while (part != 0) {
      Int_t mcGen = 1;
      Byte_t findable = 0;

      if (part->GetGeneratorIndex() == 0) mcGen = 0;

      Int_t pdg = TMath::Abs(part->PdgCode());
      // select charged pions, protons, kaons , electrons, muons
      if (pdg == 211 || pdg == 2212 || pdg == 321 || pdg == 11 || pdg == 13) findable = 1;

      FillGeneratorLevelTHnSparse(fCent, part->Eta(), part->Phi(), part->Pt(), mcGen, findable);    
      part = fGeneratorLevel->GetNextAcceptMCParticle();
    }
  }

  return kTRUE;
}
