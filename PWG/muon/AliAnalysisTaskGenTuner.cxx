/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

#include <Riostream.h>
#include <cassert>

// ROOT includes
#include <TSystem.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1D.h>
#include <TF1.h>
#include <TArrayD.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TObjArray.h>
#include <TLegend.h>
#include <THashList.h>
#include <TParameter.h>

// STEER includes
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODDimuon.h"
#include "AliCounterCollection.h"

// ANALYSIS includes
#include "AliMultSelection.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskGenTuner.h"

using std::ifstream;


ClassImp(AliAnalysisTaskGenTuner)

//________________________________________________________________________
AliAnalysisTaskGenTuner::AliAnalysisTaskGenTuner() :
AliAnalysisTaskSE(),
fList(0x0),
fEventCounters(0x0),
fCentMin(-FLT_MAX),
fCentMax(FLT_MAX),
fMuonTrackCuts(0x0),
fPtCut(-1.),
fGenPtCut(-1.),
fCharge(0),
fWeight(kFALSE),
fDataFile(""),
fPtFuncOld(0x0),
fPtFuncNew(0x0),
fPtFunc(0x0),
fPtFuncMC(0x0),
fYFuncOld(0x0),
fYFuncNew(0x0),
fYFunc(0x0),
fYFuncMC(0x0),
fMuPlusFracOld(-1.),
fMuPlusFracNew(-1.),
fMuPlusFrac(-1.),
fRunWeights(0x0),
fRunWeight(1.),
fcRes(0x0),
fcRat(0x0)
{
  /// Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskGenTuner::AliAnalysisTaskGenTuner(const char *name) :
AliAnalysisTaskSE(name),
fList(0x0),
fEventCounters(0x0),
fCentMin(-FLT_MAX),
fCentMax(FLT_MAX),
fMuonTrackCuts(0x0),
fPtCut(-1.),
fGenPtCut(-1.),
fCharge(0),
fWeight(kFALSE),
fDataFile(""),
fPtFuncOld(0x0),
fPtFuncNew(0x0),
fPtFunc(0x0),
fPtFuncMC(0x0),
fYFuncOld(0x0),
fYFuncNew(0x0),
fYFunc(0x0),
fYFuncMC(0x0),
fMuPlusFracOld(-1.),
fMuPlusFracNew(-1.),
fMuPlusFrac(-1.),
fRunWeights(0x0),
fRunWeight(1.),
fcRes(0x0),
fcRat(0x0)
{
  /// Constructor
  
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,TObjArray::Class());
  DefineOutput(2,AliCounterCollection::Class());
  
}

//________________________________________________________________________
AliAnalysisTaskGenTuner::~AliAnalysisTaskGenTuner()
{
  /// Destructor
  
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fList;
    delete fEventCounters;
  }
  delete fMuonTrackCuts;
  delete fPtFuncOld;
  delete fPtFuncNew;
  delete fPtFunc;
  delete fPtFuncMC;
  delete fYFuncOld;
  delete fYFuncNew;
  delete fYFunc;
  delete fYFuncMC;
  delete fRunWeights;
  delete fcRes;
  delete fcRat;
  
}

//___________________________________________________________________________
void AliAnalysisTaskGenTuner::UserCreateOutputObjects()
{
  /// Create histograms
  
  // initialize histos
  fList = new TObjArray(2000);
  fList->SetOwner();
  
  TH1* hPtGen = new TH1D("hPtGen","generated p_{T} distribution;p_{T} (GeV/c);dN/dp_{T}", 300, 0., 30.);
  hPtGen->Sumw2();
  fList->AddAtAndExpand(hPtGen, kPtGen);
  TH1* hPtRec = new TH1D("hPtRec","reconstructed p_{T} distribution;p_{T} (GeV/c);dN/dp_{T}", 300, 0., 30.);
  hPtRec->Sumw2();
  fList->AddAtAndExpand(hPtRec, kPtRec);
  
  TH1* hYGen = new TH1D("hYGen","generated y distribution;y;dN/dy", 250, -4.5, -2.);
  hYGen->Sumw2();
  fList->AddAtAndExpand(hYGen, kYGen);
  TH1* hYRec = new TH1D("hYRec","reconstructed y distribution;y;dN/dy", 250, -4.5, -2.);
  hYRec->Sumw2();
  fList->AddAtAndExpand(hYRec, kYRec);
  
  TH1* hPhiGen = new TH1D("hPhiGen","generated #phi distribution;#phi (deg);dN/d#phi", 360, 0., 360.);
  hPhiGen->Sumw2();
  fList->AddAtAndExpand(hPhiGen, kPhiGen);
  TH1* hPhiRec = new TH1D("hPhiRec","reconstructed #phi distribution;#phi (deg);dN/d#phi", 360, 0., 360.);
  hPhiRec->Sumw2();
  fList->AddAtAndExpand(hPhiRec, kPhiRec);
  
  TH1* hSignGen = new TH1D("hSignGen","generated sign distribution;sign;dN/dsign", 2, -2., 2.);
  hSignGen->Sumw2();
  fList->AddAtAndExpand(hSignGen, kSignGen);
  TH1* hSignRec = new TH1D("hSignRec","reconstructed sign distribution;sign;dN/dsign", 2, -2., 2.);
  hSignRec->Sumw2();
  fList->AddAtAndExpand(hSignRec, kSignRec);
  
  // initialize event counters
  fEventCounters = new AliCounterCollection(GetOutputSlot(2)->GetContainer()->GetName());
  fEventCounters->AddRubric("run", 100000);
  fEventCounters->Init();
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fList);
  PostData(2, fEventCounters);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::NotifyRun()
{
  /// Prepare processing of new run: load corresponding OADB objects...
  
  // get the trackCuts for this run
  if (!fMuonTrackCuts) AliFatal("You must specify the requested selections (AliMuonTrackCut obj is missing)");
  fMuonTrackCuts->SetRun(fInputHandler);
  
  // set the current run weight
  if (fRunWeights) {
    TParameter<Double_t> *weight = static_cast<TParameter<Double_t>*>(fRunWeights->FindObject(Form("%d", fCurrentRunNumber)));
    fRunWeight = weight ? weight->GetVal() : 0.;
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::UserExec(Option_t *)
{
  /// process event
  
  // get AOD event
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());
  if ( !aod ) return;
  
  // trigger selection
  //if (!aod->GetFiredTriggerClasses().Contains("CINT7-B-NOPF-MUFAST")) return;
  
  // select the centrality range
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(aod->FindListObject("MultSelection"));
  Float_t centrality = multSelection ? multSelection->GetMultiplicityPercentile("V0M") : -1.;
  if (centrality < fCentMin || centrality > fCentMax) return;
  
  // analyzed events
  fEventCounters->Count(Form("run:%d", fCurrentRunNumber));
  
  // fill the MC part if running on MC
  TArrayD weight;
  if (MCEvent()) {
    
    weight.Set(MCEvent()->GetNumberOfTracks());
    for (Int_t i = 0; i < MCEvent()->GetNumberOfTracks(); i++) {
      
      AliAODMCParticle *mctrack = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(i));
      weight[i] = 0.;
      
      if (!mctrack->IsPrimary()) continue;
      
      if (fCharge*mctrack->Charge() < 0) continue;
      
      // compute the weights for all primary particles (other weights are 0)
      Double_t y = mctrack->Y();
      Double_t pT = mctrack->Pt();
      if (fGenPtCut > 0. && pT < fGenPtCut) weight[i] = 0.;
      else if (fWeight && fPtFuncOld && fPtFuncNew && fYFuncOld && fYFuncNew) {
	weight[i] = fPtFuncNew->Eval(pT) / fPtFuncOld->Eval(pT) * fYFuncNew->Eval(y) / fYFuncOld->Eval(y);
        if (fMuPlusFracOld > 0. && fMuPlusFracOld < 1. && fMuPlusFracNew >= 0. && fMuPlusFracNew <= 1.) {
          weight[i] *= (mctrack->Charge() > 0.) ? fMuPlusFracNew / fMuPlusFracOld : (1.-fMuPlusFracNew) / (1.-fMuPlusFracOld);
        }
	if (weight[i] < 0.) {
	  AliError(Form("negative weight: y = %g, pT = %g: w = %g", y, pT, weight[i]));
	  weight[i] = 0.;
	}
      } else weight[i] = 1.;
      
      Double_t w = fRunWeight * weight[i];
      ((TH1*)fList->UncheckedAt(kPtGen))->Fill(pT, w);
      ((TH1*)fList->UncheckedAt(kYGen))->Fill(y, w);
      ((TH1*)fList->UncheckedAt(kPhiGen))->Fill(mctrack->Phi()*TMath::RadToDeg(), w);
      ((TH1*)fList->UncheckedAt(kSignGen))->Fill((mctrack->Charge() < 0.) ? -1. : 1., w);
      
    }
    
  }
  
  // fill the reconstructed part
  for (Int_t i = 0; i < aod->GetNumberOfTracks(); i++){
    
    AliAODTrack *track = static_cast<AliAODTrack*>(aod->GetTrack(i));
    
    Double_t pT = track->Pt();
    Int_t mcLabel = track->GetLabel();
    if (!fMuonTrackCuts->IsSelected(track) || (fPtCut > 0. && pT < fPtCut) ||
	fCharge*track->Charge() < 0 || (MCEvent() && mcLabel < 0)) continue;
    
    Double_t w = (weight.GetSize() > 0) ? fRunWeight * weight[mcLabel] : fRunWeight;
    ((TH1*)fList->UncheckedAt(kPtRec))->Fill(pT, w);
    ((TH1*)fList->UncheckedAt(kYRec))->Fill(track->Y(), w);
    ((TH1*)fList->UncheckedAt(kPhiRec))->Fill(track->Phi()*TMath::RadToDeg(), w);
    ((TH1*)fList->UncheckedAt(kSignRec))->Fill(track->Charge(), w);
    
  }
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fList);
  PostData(2, fEventCounters);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::Terminate(Option_t *)
{
  /// post-processing
  
  // get current results
  Int_t hIndex[8] = {kPtGen, kYGen, kPhiGen, kSignGen, kPtRec, kYRec, kPhiRec, kSignRec};
  fList = static_cast<TObjArray*>(GetOutputData(1));
  TH1 *h[8];
  for (Int_t i = 0; i < 8; i++) {
    h[i] = static_cast<TH1*>(fList->UncheckedAt(hIndex[i])->Clone());
    h[i]->SetDirectory(0);
  }
  
  // get the fit ranges
  Double_t fitRangeMC[4][2];
  fitRangeMC[0][0] = GetFitLowEdge(*(h[0]));
  fitRangeMC[0][1] = 999.;
  fitRangeMC[1][0] = GetFitLowEdge(*(h[1]));
  fitRangeMC[1][1] = GetFitUpEdge(*(h[1]));
  fitRangeMC[2][0] = h[2]->GetXaxis()->GetXmin();
  fitRangeMC[2][1] = h[2]->GetXaxis()->GetXmax();
  fitRangeMC[3][0] = h[3]->GetXaxis()->GetXmin();
  fitRangeMC[3][1] = h[3]->GetXaxis()->GetXmax();
  Double_t fitRange[4][2];
  fitRange[0][0] = (fPtCut > 0.) ? TMath::Max(fitRangeMC[0][0], fPtCut) : fitRangeMC[0][0];
  fitRange[0][1] = h[4]->GetXaxis()->GetXmax();
  fitRange[1][0] = -3.98; // not -4. because of the influence of the eta cut
  fitRange[1][1] = -2.51;
  fitRange[2][0] = h[6]->GetXaxis()->GetXmin();
  fitRange[2][1] = h[6]->GetXaxis()->GetXmax();
  fitRange[3][0] = h[7]->GetXaxis()->GetXmin();
  fitRange[3][1] = h[7]->GetXaxis()->GetXmax();
  
  // compute acc*eff corrections if it is simulated data
  TH1 *hAccEff[4] = {0x0, 0x0, 0x0, 0x0};
  for (Int_t i = 0; i < 4 && h[i]->GetEntries() > 0; i++) {
    hAccEff[i] = ComputeAccEff(*(h[i]), *(h[i+4]), Form("%sOverGen",h[i+4]->GetName()), "Acc#{times}Eff");
  }
  
  // get reference data if provided
  TH1 *hRef[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  if (!fDataFile.IsNull()) {
    TFile* dataFile = TFile::Open(fDataFile.Data(),"READ");
    if (!dataFile || !dataFile->IsOpen()) return;
    TObjArray* data = static_cast<TObjArray*>(dataFile->FindObjectAny("Histograms"));
    if (!data) return;
    for (Int_t i = 4; i < 8; i++) {
      if (hIndex[i] < data->GetEntries()) {
        hRef[i] = static_cast<TH1*>(data->UncheckedAt(hIndex[i])->Clone());
        hRef[i]->SetDirectory(0);
      }
    }
    dataFile->Close();
  }
  
  // compute corrected data
  for (Int_t i = 0; i < 4 && hRef[i+4] && hAccEff[i]; i++) {
    hRef[i] = static_cast<TH1*>(hRef[i+4]->Clone(Form("%sCorr",hRef[i+4]->GetName())));
    hRef[i]->SetTitle("corrected data");
    hRef[i]->Divide(hAccEff[i]);
  }
  
  // get the (weighted) number of generated(corrected)/reconstructed muons
  Double_t nMu[2] = {0., 0.};
  nMu[0] = h[3]->Integral();
  nMu[1] = h[7]->Integral();
  Double_t nMuRef[2] = {0., 0.};
  if (hRef[3]) nMuRef[0] = hRef[3]->Integral();
  if (hRef[7]) nMuRef[1] = hRef[7]->Integral();
  
  // normalize histograms
  Bool_t normalized = kFALSE;
  for (Int_t i = 0; i < 4 && hRef[i]; i++) {
    Double_t integral = h[i]->Integral(h[i]->FindBin(fitRange[i][0]), h[i]->FindBin(fitRange[i][1]), "width");
    Double_t norm = (integral != 0.) ? 1./integral : 1.;
    h[i]->Scale(norm);
    h[i+4]->Scale(norm);
    integral = hRef[i]->Integral(hRef[i]->FindBin(fitRange[i][0]), hRef[i]->FindBin(fitRange[i][1]), "width");
    norm = (integral != 0.) ? 1./integral : 1.;
    hRef[i]->Scale(norm);
    hRef[i+4]->Scale(norm);
    normalized = kTRUE;
  }
  
  // compute data/MC ratios
  TH1 *hRat[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  for (Int_t i = 0; i < 8; i++) {
    if (hRef[i]) {
      hRat[i] = static_cast<TH1*>(hRef[i]->Clone(Form("%sDataOverMC",hRef[i]->GetName())));
      hRat[i]->SetTitle("data / MC");
      hRat[i]->Divide(h[i]);
    }
  }
  
  // prepare fitting functions depending whether the original distributions are weighted or not
  if (hAccEff[0]) {
    if (fWeight && fPtFuncNew) {
      fPtFunc = new TF1(*fPtFuncNew);
      fPtFuncMC = new TF1(*fPtFuncNew);
    } else if (!fWeight && fPtFuncOld) {
      fPtFunc = new TF1(*fPtFuncOld);
      fPtFuncMC = new TF1(*fPtFuncOld);
    }
    if (fPtFunc) {
      fPtFunc->SetName("fPtFunc");
      fPtFunc->SetRange(fitRange[0][0], fitRange[0][1]);
      NormFunc(fPtFunc, fitRange[0][0], fitRange[0][1]);
    }
    if (fPtFuncMC) {
      fPtFuncMC->SetName("fPtFuncMC");
      fPtFuncMC->SetRange(fitRangeMC[0][0], fitRangeMC[0][1]);
      NormFunc(fPtFuncMC, fitRange[0][0], fitRange[0][1]);
    }
    if (hRef[0] && fPtFuncNew) {
      fPtFuncNew->SetRange(fitRange[0][0], fitRange[0][1]);
      NormFunc(fPtFuncNew, fitRange[0][0], fitRange[0][1]);
    }
    if (!normalized) {
      Double_t integral = h[0]->Integral(h[0]->FindBin(fitRange[0][0]), h[0]->FindBin(fitRange[0][1]), "width");
      if (fPtFunc) fPtFunc->SetParameter(0, fPtFunc->GetParameter(0)*integral);
      if (fPtFuncMC) fPtFuncMC->SetParameter(0, fPtFuncMC->GetParameter(0)*integral);
      if (fPtFuncNew) {
	integral = hRef[0]->Integral(hRef[0]->FindBin(fitRange[0][0]), hRef[0]->FindBin(fitRange[0][1]), "width");
	fPtFuncNew->SetParameter(0, fPtFuncNew->GetParameter(0)*integral);
      }
    }
  }
  if (hAccEff[1]) {
    if (fWeight && fYFuncNew) {
      fYFunc = new TF1(*fYFuncNew);
      fYFuncMC = new TF1(*fYFuncNew);
    } else if (!fWeight && fYFuncOld) {
      fYFunc = new TF1(*fYFuncOld);
      fYFuncMC = new TF1(*fYFuncOld);
    }
    if (fYFunc) {
      fYFunc->SetName("fYFunc");
      fYFunc->SetRange(fitRange[1][0], fitRange[1][1]);
      NormFunc(fYFunc, fitRange[1][0], fitRange[1][1]);
    }
    if (fYFuncMC) {
      fYFuncMC->SetName("fYFuncMC");
      fYFuncMC->SetRange(fitRangeMC[1][0], fitRangeMC[1][1]);
      NormFunc(fYFuncMC, fitRange[1][0], fitRange[1][1]);
    }
    if (hRef[1] && fYFuncNew) {
      fYFuncNew->SetRange(fitRange[1][0], fitRange[1][1]);
      NormFunc(fYFuncNew, fitRange[1][0], fitRange[1][1]);
    }
    if (!normalized) {
      Double_t integral = h[1]->Integral(h[1]->FindBin(fitRange[1][0]), h[1]->FindBin(fitRange[1][1]), "width");
      if (fYFunc) fYFunc->SetParameter(0, fYFunc->GetParameter(0)*integral);
      if (fYFuncMC) fYFuncMC->SetParameter(0, fYFuncMC->GetParameter(0)*integral);
      if (fYFuncNew) {
	integral = hRef[1]->Integral(hRef[1]->FindBin(fitRange[1][0]), hRef[1]->FindBin(fitRange[1][1]), "width");
	fYFuncNew->SetParameter(0, fYFuncNew->GetParameter(0)*integral);
      }
    }
  }
  
  // plot and fit results
  TLegend *lRes = new TLegend(0.5,0.55,0.85,0.75);
  fcRes = new TCanvas("cRes", "results", 1200, 600);
  fcRes->Divide(4,2);
  fcRes->cd(1);
  gPad->SetLogy();
  if (hAccEff[0]) {
    if (fPtFuncMC) {
      fPtFuncMC->SetLineColor(3);
      fPtFuncMC->SetLineWidth(3);
      h[0]->Fit(fPtFuncMC, "GWLMR", "e0sames");
      lRes->AddEntry(fPtFuncMC,"MC range MC","l");
    } else h[0]->Draw("e0");
    if (fPtFunc) {
      fPtFunc->SetLineColor(4);
      h[0]->Fit(fPtFunc, "GWLMR+");
      lRes->AddEntry(fPtFunc,"MC range Data","l");
    }
  }
  if (hRef[0]) {
    hRef[0]->SetLineColor(2);
    if (fPtFuncNew) {
      fPtFuncNew->SetLineColor(2);
      hRef[0]->Fit(fPtFuncNew, "GWLMR", "e0sames");
      lRes->AddEntry(fPtFuncNew,"Data","l");
    } else hRef[0]->Draw("e0sames");
  }
  lRes->Draw("same");
  fcRes->cd(2);
  if (hAccEff[1]) {
    if (fYFuncMC) {
      fYFuncMC->SetLineColor(3);
      fYFuncMC->SetLineWidth(3);
      h[1]->Fit(fYFuncMC, "GWLMR", "e0sames");
    } else h[1]->Draw("e0");
    if (fYFunc) {
      fYFunc->SetLineColor(4);
      h[1]->Fit(fYFunc, "GWLMR+");
    }
  }
  if (hRef[1]) {
    hRef[1]->SetLineColor(2);
    if (fYFuncNew) {
      fYFuncNew->SetLineColor(2);
      hRef[1]->Fit(fYFuncNew, "GWLMR", "e0sames");
    } else hRef[1]->Draw("e0sames");
  }
  for (Int_t i = 2; i < 8; i++) {
    fcRes->cd(i+1);
    if (i == 4) gPad->SetLogy();
    h[i]->Draw("e0");
    if (hRef[i]) {
      hRef[i]->SetLineColor(2);
      hRef[i]->Draw("e0sames");
    }
  }
  
  // normalize functions to their integral in the range used in MC
  if (hAccEff[0] && fPtFunc) {
    fPtFunc->SetRange(fitRangeMC[0][0], fitRangeMC[0][1]);
    NormFunc(fPtFunc, fitRangeMC[0][0], fitRangeMC[0][1]);
  }
  if (hAccEff[0] && fPtFuncMC) {
    NormFunc(fPtFuncMC, fitRangeMC[0][0], fitRangeMC[0][1]);
  }
  if (hRef[0] && fPtFuncNew) {
    fPtFuncNew->SetRange(fitRangeMC[0][0], fitRangeMC[0][1]);
    NormFunc(fPtFuncNew, fitRangeMC[0][0], fitRangeMC[0][1]);
  }
  if (hAccEff[1] && fYFunc) {
    fYFunc->SetRange(fitRangeMC[1][0], fitRangeMC[1][1]);
    NormFunc(fYFunc, fitRangeMC[1][0], fitRangeMC[1][1]);
  }
  if (hAccEff[1] && fYFuncMC) {
    NormFunc(fYFuncMC, fitRangeMC[1][0], fitRangeMC[1][1]);
  }
  if (hRef[1] && fYFuncNew) {
    fYFuncNew->SetRange(fitRangeMC[1][0], fitRangeMC[1][1]);
    NormFunc(fYFuncNew, fitRangeMC[1][0], fitRangeMC[1][1]);
  }
  
  // Get fractions of mu+
  if (hAccEff[3]) {
    Double_t sum  = h[3]->GetBinContent(1) + h[3]->GetBinContent(2);
    if (sum > 0.) fMuPlusFrac = h[3]->GetBinContent(2) / sum;
  }
  if (hRef[3]) {
    Double_t sum  = hRef[3]->GetBinContent(1) + hRef[3]->GetBinContent(2);
    if (sum > 0.) fMuPlusFracNew = hRef[3]->GetBinContent(2) / sum;
  }
  
  // prepare data/MC function ratios
  TF1 *ptRat = (hRat[0] && fPtFunc && fPtFuncNew) ? new TF1("ptRat", this, &AliAnalysisTaskGenTuner::PtRat, fitRangeMC[0][0], hRat[0]->GetXaxis()->GetXmax(), 0, "AliAnalysisTaskGenTuner", "PtRat") : 0x0;
  TF1 *yRat = (hRat[1] && fYFunc && fYFuncNew) ? new TF1("yRat", this, &AliAnalysisTaskGenTuner::YRat, fitRangeMC[1][0], fitRangeMC[1][1], 0, "AliAnalysisTaskGenTuner", "YRat") : 0x0;
  
  // plot ratios
  fcRat = new TCanvas("cRat", "ratios", 1200, 600);
  fcRat->Divide(4,2);
  for (Int_t i = 0; i < 8; i++) {
    if (hRat[i]) {
      fcRat->cd(i+1);
      hRat[i]->Draw("e0");
      if (i == 0 && ptRat) ptRat->Draw("sames");
      else if (i == 1 && yRat) yRat->Draw("sames");
    }
  }
  
  // print fitting ranges
  if (fPtFuncMC && fYFuncMC) {
    printf("\npT fitting range MC = [%g, %g]\n", fitRangeMC[0][0], fitRangeMC[0][1]);
    printf("y fitting range MC = [%g, %g]\n\n", fitRangeMC[1][0], fitRangeMC[1][1]);
  }
  if (fPtFunc && fYFunc) {
    printf("pT fitting range = [%g, %g]\n", fitRange[0][0], fitRange[0][1]);
    printf("y fitting range = [%g, %g]\n\n", fitRange[1][0], fitRange[1][1]);
  }
  
  // print parameters
  if (fPtFuncMC) {
    printf("Double_t oldPtParamMC[%d] = {", fPtFuncMC->GetNpar());
    for (Int_t i = 0; i < fPtFuncMC->GetNpar()-1; i++) printf("%g, ", fPtFuncMC->GetParameter(i));
    printf("%g};\n", fPtFuncMC->GetParameter(fPtFuncMC->GetNpar()-1));
  }
  if (fYFuncMC) {
    printf("Double_t oldYParamMC[%d] = {", fYFuncMC->GetNpar());
    for (Int_t i = 0; i < fYFuncMC->GetNpar()-1; i++) printf("%g, ", fYFuncMC->GetParameter(i));
    printf("%g};\n\n", fYFuncMC->GetParameter(fYFuncMC->GetNpar()-1));
  }
  if (fPtFunc) {
    printf("Double_t oldPtParam[%d] = {", fPtFunc->GetNpar());
    for (Int_t i = 0; i < fPtFunc->GetNpar()-1; i++) printf("%g, ", fPtFunc->GetParameter(i));
    printf("%g};\n", fPtFunc->GetParameter(fPtFunc->GetNpar()-1));
  }
  if (fYFunc) {
    printf("Double_t oldYParam[%d] = {", fYFunc->GetNpar());
    for (Int_t i = 0; i < fYFunc->GetNpar()-1; i++) printf("%g, ", fYFunc->GetParameter(i));
    printf("%g};\n\n", fYFunc->GetParameter(fYFunc->GetNpar()-1));
  }
  if (fPtFuncNew) {
    printf("Double_t newPtParam[%d] = {", fPtFuncNew->GetNpar());
    for (Int_t i = 0; i < fPtFuncNew->GetNpar()-1; i++) printf("%g, ", fPtFuncNew->GetParameter(i));
    printf("%g};\n", fPtFuncNew->GetParameter(fPtFuncNew->GetNpar()-1));
  }
  if (fYFuncNew) {
    printf("Double_t newYParam[%d] = {", fYFuncNew->GetNpar());
    for (Int_t i = 0; i < fYFuncNew->GetNpar()-1; i++) printf("%g, ", fYFuncNew->GetParameter(i));
    printf("%g};\n\n", fYFuncNew->GetParameter(fYFuncNew->GetNpar()-1));
  }
  if (hAccEff[3]) printf("Double_t oldMuPlusFrac = %f\n", fMuPlusFrac);
  if (hRef[3]) printf("Double_t newMuPlusFrac = %f\n\n", fMuPlusFracNew);
  if (hAccEff[3]) printf("MC  muons: generated (MC range): %g ; reconstructed (reco range): %g --> ratio = %g\n", nMu[0], nMu[1], (nMu[0]>0.) ? nMu[1]/nMu[0] : 0.);
  if (hRef[7]) printf("Ref muons: corrected (MC range): %g ; reconstructed (reco range): %g\n\n", nMuRef[0], nMuRef[1]);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::SetOriginPtFunc(TString formula, const Double_t *param, const Bool_t *fixParam,
                                              Double_t xMin, Double_t xMax)
{
  /// Create the original function with the parameters used in simulation to generate the pT distribution.
  /// It must be in the form [0]*(...) to allow for global normalization.
  /// The [xMin,xMax] range is used to normalized the function.
  /// Some parameters can be fixed when fitting the generated distribution for cross-check.
  
  assert(param);
  
  delete fPtFuncOld;
  fPtFuncOld = new TF1("fPtFuncOld", formula.Data(), xMin, xMax);
  
  fPtFuncOld->SetParameters(param);
  if (fixParam) for (Int_t i = 0; i < fPtFuncOld->GetNpar(); ++i)
    if (fixParam[i]) fPtFuncOld->FixParameter(i, fPtFuncOld->GetParameter(i));
  
  NormFunc(fPtFuncOld, xMin, xMax);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::SetNewPtFunc(TString formula, const Double_t *param, const Bool_t *fixParam,
                                           Double_t xMin, Double_t xMax)
{
  /// Create the new function with its initial parameters to fit the generated/weighted pT distribution.
  /// It must be in the form [0]*(...) to allow for global normalization.
  /// The [xMin,xMax] range is used to normalized the function.
  /// Some parameters can be fixed when fitting the generated distribution.
  
  assert(param);
  
  delete fPtFuncNew;
  fPtFuncNew = new TF1("fPtFuncNew", formula.Data(), xMin, xMax);
  
  fPtFuncNew->SetParameters(param);
  if (fixParam) for (Int_t i = 0; i < fPtFuncNew->GetNpar(); ++i)
    if (fixParam[i]) fPtFuncNew->FixParameter(i, fPtFuncNew->GetParameter(i));
  
  NormFunc(fPtFuncNew, xMin, xMax);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::SetOriginYFunc(TString formula, const Double_t *param, const Bool_t *fixParam,
                                             Double_t xMin, Double_t xMax)
{
  /// Create the original function with the parameters used in simulation to generate the y distribution.
  /// It must be in the form [0]*(...) to allow for global normalization.
  /// The [xMin,xMax] range is used to normalized the function.
  /// Some parameters can be fixed when fitting the generated distribution for cross-check.
  
  assert(param);
  
  delete fYFuncOld;
  fYFuncOld = new TF1("fYFuncOld", formula.Data(), xMin, xMax);
  
  fYFuncOld->SetParameters(param);
  if (fixParam) for (Int_t i = 0; i < fYFuncOld->GetNpar(); ++i)
    if (fixParam[i]) fYFuncOld->FixParameter(i, fYFuncOld->GetParameter(i));
  
  NormFunc(fYFuncOld, xMin, xMax);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::SetNewYFunc(TString formula, const Double_t *param, const Bool_t *fixParam,
                                          Double_t xMin, Double_t xMax)
{
  /// Create the new function with its initial parameters to fit the generated/weighted y distribution.
  /// It must be in the form [0]*(...) to allow for global normalization.
  /// The [xMin,xMax] range is used to normalized the function.
  /// Some parameters can be fixed when fitting the generated distribution.
  
  assert(param);
  
  delete fYFuncNew;
  fYFuncNew = new TF1("fYFuncNew", formula.Data(), xMin, xMax);
  
  fYFuncNew->SetParameters(param);
  if (fixParam) for (Int_t i = 0; i < fYFuncNew->GetNpar(); ++i)
    if (fixParam[i]) fYFuncNew->FixParameter(i, fYFuncNew->GetParameter(i));
  
  NormFunc(fYFuncNew, xMin, xMax);
  
}

//________________________________________________________________________
TH1* AliAnalysisTaskGenTuner::ComputeAccEff(TH1 &hGen, TH1 &hRec, const Char_t *name, const Char_t *title)
{
  /// Compute acc*eff and binomial errors by hand, i.e. not using TGraphAsymmErrors
  /// Result is identical to divide histograms with option "B", except here error is forced > 1/gen
  
  Int_t nbins = hGen.GetNbinsX();
  TH1* hAcc = new TH1D(name,title, nbins, hGen.GetXaxis()->GetXmin(), hGen.GetXaxis()->GetXmax());
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t accEff = 0.;
    Double_t accEffErr = 0.;
    Double_t gen = hGen.GetBinContent(i);
    if (gen > 0.) {
      Double_t rec = hRec.GetBinContent(i);
      Double_t genErr = hGen.GetBinError(i);
      Double_t recErr = hRec.GetBinError(i);
      accEff = rec/gen;
      //Double_t accEffErr2 = ((1.-2.*accEff)*recErr*recErr + accEff*accEff*genErr*genErr)/(gen*gen);
      //accEffErr = TMath::Max(accEff*genErr/gen, TMath::Sqrt(TMath::Abs(accEffErr2)));
      accEffErr = TMath::Max(recErr/gen, accEff*genErr/gen);
    }
    hAcc->SetBinContent(i, accEff);
    hAcc->SetBinError(i, accEffErr);
  }
  
  return hAcc;
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTuner::GetFitLowEdge(TH1 &h)
{
  /// adjust the lower edge of the fit range according to the content of the histogram
  Int_t binAbove0 = h.FindFirstBinAbove(0.);
  if (h.GetBinContent(binAbove0) < 0.1*h.GetBinContent(binAbove0+1)) binAbove0++;
  return h.GetBinLowEdge(binAbove0);
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTuner::GetFitUpEdge(TH1 &h)
{
  /// adjust the upper edge of the fit range according to the content of the histogram
  Int_t binAbove0 = h.FindLastBinAbove(0.);
  if (h.GetBinContent(binAbove0) < 0.1*h.GetBinContent(binAbove0-1)) binAbove0--;
  return h.GetBinLowEdge(binAbove0+1);
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::NormFunc(TF1 *f, Double_t min, Double_t max)
{
  /// normalize the function to its integral in the given range
  Double_t integral = f->Integral(min, max);
  if (integral != 0.) f->SetParameter(0, f->GetParameter(0)/integral);
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTuner::PtRat(const Double_t *x, const Double_t */*p*/)
{
  /// generated pT fit function ratio
  return (fPtFunc && fPtFuncNew) ? fPtFuncNew->Eval(*x) / fPtFunc->Eval(*x) : 0.;
}

//________________________________________________________________________
Double_t AliAnalysisTaskGenTuner::YRat(const Double_t *x, const Double_t */*p*/)
{
  /// generated y fit function ratio
  return (fYFunc && fYFuncNew) ? fYFuncNew->Eval(*x) / fYFunc->Eval(*x) : 0.;
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::RunWeight(TString &fileName)
{
  /// weight each run using the values in the given text file
  
  if (!gSystem->GetFromPipe(Form("file %s", fileName.Data())).Contains("text"))
    AliFatal("weights must be in a text file");
  
  delete fRunWeights;
  fRunWeights = LoadRunWeights(fileName);
  
}

//________________________________________________________________________
void AliAnalysisTaskGenTuner::RunWeight(TString &fileNameOrigin, TString &fileNameNew)
{
  /// weight each run using the original/new number of events in the given text of root files
  
  delete fRunWeights;
  fRunWeights = new THashList(1000);
  fRunWeights->SetOwner();
  
  THashList *nEvsOld = LoadRunWeights(fileNameOrigin);
  if (!nEvsOld) return;
  nEvsOld->Sort();
  
  THashList *nEvsNew = LoadRunWeights(fileNameNew);
  if (!nEvsNew) return;
  nEvsNew->Sort();
  
  if (nEvsOld->GetSize() != nEvsNew->GetSize()) AliFatal("the 2 files do not contain the same number of runs");
  
  TIter nextRunOld(nEvsOld);
  TIter nextRunNew(nEvsNew);
  TParameter<Double_t> *nEvOld = 0x0, *nEvNew = 0x0;
  while ((nEvOld = static_cast<TParameter<Double_t>*>(nextRunOld())) &&
         (nEvNew = static_cast<TParameter<Double_t>*>(nextRunNew()))) {
    
    const char *run = nEvNew->GetName();
    if (!strstr(nEvOld->GetName(), run)) AliFatal("the 2 files do not contain the same runs");
    
    Double_t weight = (nEvOld->GetVal() > 0) ? nEvNew->GetVal() / nEvOld->GetVal() : 0.;
    fRunWeights->Add(new TParameter<Double_t>(run, weight));
    
  }
  
  delete nEvsOld;
  delete nEvsNew;
  
}

//________________________________________________________________________
THashList* AliAnalysisTaskGenTuner::LoadRunWeights(const TString &fileName)
{
  /// Load the weights (or number of events) per run
  
  TString fileType = gSystem->GetFromPipe(Form("file %s", fileName.Data()));
  
  if (fileType.Contains("ROOT")) return LoadRunWeightsFromRootFile(fileName);
  else if (fileType.Contains("text")) return LoadRunWeightsFromTextFile(fileName);
  else if (fileType.Contains("No such file")) AliFatal(Form("file %s not found", fileName.Data()));
  else AliFatal("number of events must be either in a text file or in an AliCounterCollection in a root file");
  
  return 0x0;
  
}

//________________________________________________________________________
THashList* AliAnalysisTaskGenTuner::LoadRunWeightsFromTextFile(const TString &fileName)
{
  /// Load the weights (or number of events) per run from the given text file
  
  THashList *weights = new THashList(1000);
  weights->SetOwner();
  
  ifstream inFile(fileName);
  if (!inFile.is_open()) AliFatal(Form("cannot open file %s", fileName.Data()));
  
  TString line;
  while (! inFile.eof() ) {
    
    line.ReadLine(inFile,kTRUE);
    if(line.IsNull()) continue;
    
    TObjArray *param = line.Tokenize(" ");
    if (param->GetEntries() != 2) AliFatal(Form("bad input line %s", line.Data()));
    
    Int_t run = ((TObjString*)param->UncheckedAt(0))->String().Atoi();
    if (run < 0) AliFatal(Form("invalid run number: %d", run));
    
    Float_t weight = ((TObjString*)param->UncheckedAt(1))->String().Atof();
    if (weight < 0.) AliFatal(Form("invalid weight: %g", weight));
    
    if (weights->FindObject(Form("%d",run))) AliFatal(Form("weight for run %d already exist", run));
    
    weights->Add(new TParameter<Double_t>(Form("%d",run), weight));
    
    delete param;
  }
  
  inFile.close();
  
  return weights;
  
}

//________________________________________________________________________
THashList* AliAnalysisTaskGenTuner::LoadRunWeightsFromRootFile(const TString &fileName)
{
  /// Load the number of events per run from the AliCounterCollection in the given root file
  
  THashList *weights = new THashList(1000);
  weights->SetOwner();
  
  TFile* inFile = TFile::Open(fileName.Data(),"READ");
  if (!inFile || !inFile->IsOpen()) {
    AliFatal(Form("cannot open file %s", fileName.Data()));
    return 0x0;
  }
  
  AliCounterCollection* eventCounters = static_cast<AliCounterCollection*>(inFile->FindObjectAny("eventCounters"));
  if (!eventCounters) {
    AliFatal(Form("cannot find eventCounters in file %s", fileName.Data()));
    return 0x0;
  }
  
  TH1D *hnEvs = eventCounters->Get("run", "");
  if (!hnEvs) {
    AliFatal(Form("invalid eventCounters in file %s", fileName.Data()));
    return 0x0;
  }
  hnEvs->SetDirectory(0);
  
  for (Int_t irun = 1; irun <= hnEvs->GetNbinsX(); ++irun) {
    
    Int_t run = atoi(hnEvs->GetXaxis()->GetBinLabel(irun));
    if (run < 0) AliFatal(Form("invalid run number: %d", run));
    
    Double_t weight = hnEvs->GetBinContent(irun);
    if (weight < 0.) AliFatal(Form("invalid weight: %g", weight));
    
    if (weights->FindObject(Form("%d",run))) AliFatal(Form("weight for run %d already exist", run));
    
    weights->Add(new TParameter<Double_t>(Form("%d",run), weight));
    
  }
  
  delete hnEvs;
  inFile->Close();
  
  return weights;
  
}

