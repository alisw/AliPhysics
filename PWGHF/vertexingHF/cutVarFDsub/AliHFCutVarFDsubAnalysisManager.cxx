#include "AliHFCutVarFDsubAnalysisManager.h"

#include <iostream>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"

#include "AliHFMassFitter.h"

#include "AliHFCutVarFDsubAxis.h"
#include "AliHFCutVarFDsubCut.h"
#include "AliHFCutVarFDsubCutSet.h"
#include "AliHFCutVarFDsubEfficiency.h"
#include "AliHFCutVarFDsubMassFitter.h"
#include "AliHFCutVarFDsubMinimiser.h"

/// \cond CLASSIMP
ClassImp(AliHFCutVarFDsubAnalysisManager);
/// \endcond

using std::cerr;
using std::endl;
using std::vector;

AliHFCutVarFDsubAnalysisManager::AliHFCutVarFDsubAnalysisManager()
  : TObject()
  , fData(0x0)
  , fAxes(0x0)
  , fCuts(0x0)
  , fRawYields(0x0)
  , fCorrYieldPrompt(0x0)
  , fCorrYieldFD(0x0)
  , fResiduals(0x0)
  , fPulls(0x0)
  , fFprompt(0x0)
  , fFpromptRaw(0x0)
  , fIncDistError(0x0)
  , fxAxisTitle("#it{p}_{T} (GeV/#it{c})")
  , fBinsX(0x0)
{
  /// Default constructor

  fMCgenLevel[0]  = 0x0;
  fMCgenLevel[1]  = 0x0;
  fMCafterCuts[0] = 0x0;
  fMCafterCuts[1] = 0x0;
  fEffListVsCutSets[0] = 0x0;
  fEffListVsCutSets[1] = 0x0;
  fEffListVsBins[0] = 0x0;
  fEffListVsBins[1] = 0x0;
}


//_________________________________________________________________________________________________
AliHFCutVarFDsubAnalysisManager::~AliHFCutVarFDsubAnalysisManager() {
  /// Destructor

  for (UInt_t iOrigin=kPrompt; iOrigin<=kFD; ++iOrigin) {
    if (fMCgenLevel[iOrigin])       { delete fMCgenLevel[iOrigin];       fMCgenLevel[iOrigin]       = 0x0; }
    if (fMCgenLevel[iOrigin])       { delete fMCgenLevel[iOrigin];       fMCgenLevel[iOrigin]       = 0x0; }
    if (fEffListVsBins[iOrigin])    { delete fEffListVsBins[iOrigin];    fEffListVsBins[iOrigin]    = 0x0; }
    if (fEffListVsCutSets[iOrigin]) { delete fEffListVsCutSets[iOrigin]; fEffListVsCutSets[iOrigin] = 0x0; }
  }
  if (fRawYields)       { delete   fRawYields;       fRawYields       = 0x0; }
  if (fCorrYieldPrompt) { delete   fCorrYieldPrompt; fCorrYieldPrompt = 0x0; }
  if (fCorrYieldFD)     { delete   fCorrYieldFD;     fCorrYieldFD     = 0x0; }
  if (fResiduals)       { delete   fResiduals;       fResiduals       = 0x0; }
  if (fPulls)           { delete   fPulls;           fPulls           = 0x0; }
  if (fFprompt)         { delete   fFprompt;         fFprompt         = 0x0; }
  if (fFpromptRaw)      { delete   fFpromptRaw;      fFpromptRaw      = 0x0; }
  if (fIncDistError)    { delete   fIncDistError;    fIncDistError    = 0x0; }
  if (fBinsX)           { delete[] fBinsX;           fBinsX           = 0x0; }
}


//_________________________________________________________________________________________________
void AliHFCutVarFDsubAnalysisManager::DrawDistributions(TString strOutputFolder) {

  // release pT axis
  for(Int_t iOrigin=0; iOrigin<kFD+1; ++iOrigin) {
    fMCafterCuts[iOrigin]->GetAxis(0)->SetRange(-1, -1); //pT axes should be 0;
  }

  TCanvas *cDist = new TCanvas("cDist", "", 1920, 1080);
  cDist->Clear();
  Int_t nAxes = 0;
  for(Int_t iAxis=0; iAxis<fAxes->GetEntries(); ++iAxis) {
    AliHFCutVarFDsubAxis *axis = (AliHFCutVarFDsubAxis*)fAxes->At(iAxis);
    if(axis->GetAxisNo(AliHFCutVarFDsubAxis::kMCafterCuts)<(UInt_t)-1)
      nAxes++;
    if(axis->GetAxisName()=="PID")
      nAxes = nAxes-1;
  }
  
  Int_t nPads = 1;

  if(nAxes-1<=3) {
    nPads = nAxes-1;
    cDist->Divide(nPads,1);
  }
  else if(nAxes-1>3 && nAxes-1<=6){
    if((nAxes-1)%2==0) {
      nPads = (nAxes-1)/2;
      cDist->Divide(nPads,2);
    }
    else {
      nPads = (nAxes-1)/2+1;
      cDist->Divide(nPads,2);
    }
  }
  else {
    if((nAxes-1)%3==0) {
      nPads = (nAxes-1)/3;
      cDist->Divide(nPads,3);
    }
    else {
      nPads = (nAxes-1)/3+1;
      cDist->Divide(nPads,3);
    }
  }

  //select the pT bins
  for(Int_t iBin=0; iBin<fCuts->GetEntries(); ++iBin){
    TList* cutSets = (TList*)fCuts->At(iBin);
    //takes the first cut set (every cut set should have the same pt range
    TLegend* l = new TLegend(0.65,0.8,0.9,0.9);
    l->SetTextSize(0.04);
    for(UInt_t iOrigin=kPrompt; iOrigin<=kFD; ++iOrigin) {
      AliHFCutVarFDsubCutSet* cutSet = (AliHFCutVarFDsubCutSet*)cutSets->At(0);
      AliHFCutVarFDsubCut* Ptcut = (AliHFCutVarFDsubCut*)cutSet->GetCut(0);//takes the cut on pT
      AliHFCutVarFDsubAxis* Ptaxis = (AliHFCutVarFDsubAxis*)fAxes->At(Ptcut->fAxisId);
      UInt_t PtaxisNo = Ptaxis->GetAxisNo(AliHFCutVarFDsubAxis::kMCafterCuts);
      if (PtaxisNo<(UInt_t)-1) {
        TAxis* ax = fMCafterCuts[iOrigin]->GetAxis(PtaxisNo);
        Int_t PtBinMin = ax->FindBin(Ptcut->fLow*1.0001);
        Int_t PtBinMax = ax->FindBin(Ptcut->fHigh*0.9999);
        ax->SetRange(PtBinMin,PtBinMax);
      }
      Int_t axisCounter = 0;
      for (Int_t iCut=1; iCut<cutSet->GetEntries(); ++iCut) {
        AliHFCutVarFDsubCut* cut = (AliHFCutVarFDsubCut*)cutSet->GetCut(iCut);//takes the cut on pT
        AliHFCutVarFDsubAxis* axis = (AliHFCutVarFDsubAxis*)fAxes->At(cut->fAxisId);
        UInt_t axisNo = axis->GetAxisNo(AliHFCutVarFDsubAxis::kMCafterCuts);
        if (axisNo<(UInt_t)-1 && axis->GetAxisName()!="PID") {
          axisCounter++;
          if(iOrigin==kPrompt) {
            TH1F* hProjPrompt = (TH1F*)fMCafterCuts[iOrigin]->Projection(axisNo);
            hProjPrompt->Scale(1./hProjPrompt->Integral());
            if(axisCounter==1)
              l->AddEntry(hProjPrompt,"Prompt","l");
            hProjPrompt->SetStats(kFALSE);
            hProjPrompt->GetYaxis()->SetTitle("Normalised Entries");
            Double_t max = hProjPrompt->GetMaximum()*5;
            Double_t min = 0.0001;
            hProjPrompt->GetYaxis()->SetRangeUser(min,max);
            cDist->cd(axisCounter)->SetLogy();
            hProjPrompt->Draw();
          }
          else {
            TH1F* hProjFD = (TH1F*)fMCafterCuts[iOrigin]->Projection(axisNo);
            hProjFD->Scale(1./hProjFD->Integral());
            hProjFD->SetStats(kFALSE);
            hProjFD->SetLineColor(kRed);
            if(axisCounter==1)
              l->AddEntry(hProjFD,"FD","l");
            cDist->cd(axisCounter);
            hProjFD->Draw("same");
            l->Draw("same");
          }
        }
      }
    }

    cDist->Update();

    cDist->SaveAs(Form("%s/CutVarDist_%d.pdf",  strOutputFolder.Data(), iBin));
    cDist->SaveAs(Form("%s/CutVarDist_%d.root", strOutputFolder.Data(), iBin));

    delete l;
    l=0x0;
  }

  delete cDist;
  cDist=0x0;

}


//_________________________________________________________________________________________________
void AliHFCutVarFDsubAnalysisManager::GetEfficiencies(TString strOutputFolder/*="."*/, Bool_t ptWeight/*=kFALSE*/, TF1* funcWeightsD/*=0x0*/, TF1* funcWeightsB/*=0x0*/) {
  /// Obtain the efficiencies from the THnSparses

  Int_t* nCutSets = new Int_t[fCuts->GetEntries()];
  Int_t nMaxCutSets = 0;
  for (Int_t iBin=0; iBin<fCuts->GetEntries(); ++iBin) nCutSets[iBin] = 0;
  
  if (!fBinsX) GetXaxisInformation();
  TList* set = (TList*)fCuts->At(0);
  
  TH1F*** hEff = new TH1F**[2];
  hEff[0] = new TH1F*[set->GetEntries()];
  hEff[1] = new TH1F*[set->GetEntries()];
  
  for(Int_t iCutSet=0; iCutSet<set->GetEntries(); ++iCutSet) {
    hEff[0][iCutSet] = new TH1F(Form("hEffPrompt_Set%d",iCutSet+1),"",fCuts->GetEntries(),fBinsX);
    hEff[1][iCutSet] = new TH1F(Form("hEffFD_Set%d",iCutSet+1),"",fCuts->GetEntries(),fBinsX);
  }

  std::vector< std::vector< std::vector<Double_t> > > eff;
  std::vector< std::vector< std::vector<Double_t> > > effErr;
  eff.resize(kFD+1);
  effErr.resize(kFD+1);
  for (UInt_t iOrigin=kPrompt; iOrigin<=kFD; ++iOrigin) {
    eff[iOrigin].resize(fCuts->GetEntries());
    effErr[iOrigin].resize(fCuts->GetEntries());
  }

  for (Int_t iBin=0; iBin<fCuts->GetEntries(); ++iBin) {
    TList* cutSets = (TList*)fCuts->At(iBin);
    if (nCutSets[iBin]<cutSets->GetEntries()) nCutSets[iBin] = cutSets->GetEntries();
    if (nCutSets[iBin]>nMaxCutSets) nMaxCutSets=nCutSets[iBin];
    for (Int_t iCutSet=0; iCutSet<cutSets->GetEntries(); ++iCutSet) {
      AliHFCutVarFDsubCutSet* cutSet = (AliHFCutVarFDsubCutSet*)cutSets->At(iCutSet);

      for (UInt_t iOrigin=kPrompt; iOrigin<=kFD; ++iOrigin) {
        AliHFCutVarFDsubEfficiency* efficiency = 0x0;
        if(iOrigin==kPrompt) {
          efficiency = new AliHFCutVarFDsubEfficiency(fMCgenLevel[iOrigin], fMCafterCuts[iOrigin], cutSet, fAxes, ptWeight, funcWeightsD);
        }
        else {
          efficiency = new AliHFCutVarFDsubEfficiency(fMCgenLevel[iOrigin], fMCafterCuts[iOrigin], cutSet, fAxes, ptWeight, funcWeightsB);
        }
          
        eff[iOrigin][iBin].push_back(efficiency->GetEfficiency());
        effErr[iOrigin][iBin].push_back(efficiency->GetEfficiencyError());

        hEff[iOrigin][iCutSet]->SetBinContent(iBin+1,efficiency->GetEfficiency());
        hEff[iOrigin][iCutSet]->SetBinError(iBin+1,efficiency->GetEfficiencyError());
        
        delete efficiency;
        efficiency = 0x0;
      }
    }
  }

  for (UInt_t iOrigin=kPrompt; iOrigin<=kFD; ++iOrigin) {

    TString suffix = (iOrigin==kPrompt) ? "Prompt" : "FD";

    // preparte the output lists
    if (fEffListVsBins[iOrigin]) {
      delete fEffListVsBins[iOrigin];
      fEffListVsBins[iOrigin] = 0x0;
    }
    if (fEffListVsCutSets[iOrigin]) {
      delete fEffListVsCutSets[iOrigin];
      fEffListVsCutSets[iOrigin] = 0x0;
    }
    fEffListVsBins[iOrigin] = new TList();
    fEffListVsCutSets[iOrigin] = new TList();
    fEffListVsBins[iOrigin]->SetOwner();
    fEffListVsCutSets[iOrigin]->SetOwner();

    for (Int_t iBin=0; iBin<fCuts->GetEntries(); ++iBin) {
      TH1F* h = new TH1F(Form("hEff%s_%0.1f-%0.1f", suffix.Data(), fBinsX[iBin], fBinsX[iBin+1]),
                         ";Cut Set;#varepsilon", nMaxCutSets, 0., (Double_t)nMaxCutSets);
      for (UInt_t iCutSet=0; iCutSet<eff[iOrigin][iBin].size(); ++iCutSet) {
        h->SetBinContent(iCutSet+1, eff[iOrigin][iBin][iCutSet]   );
        h->SetBinError(  iCutSet+1, effErr[iOrigin][iBin][iCutSet]);
      }
      fEffListVsCutSets[iOrigin]->Add((TObject*)h);
    }

    for (Int_t iCutSet=0; iCutSet<nMaxCutSets; ++iCutSet) {
      TH1F* h = new TH1F(Form("hEff%s_%d", suffix.Data(), iCutSet),
                         Form(";%s;#varepsilon_{%d}", fxAxisTitle.Data(), iCutSet),
                         fCuts->GetEntries(), fBinsX);
      for (Int_t iBin=0; iBin<fCuts->GetEntries(); ++iBin) {
        h->SetBinContent(iBin+1, (iCutSet<nCutSets[iBin]) ? eff[iOrigin][iBin][iCutSet]    : 0.);
        h->SetBinError(  iBin+1, (iCutSet<nCutSets[iBin]) ? effErr[iOrigin][iBin][iCutSet] : 0.);
      }
      fEffListVsBins[iOrigin]->Add((TObject*)h);
    }
  }

  TFile outfile(Form("%s/Efficiencies.root",strOutputFolder.Data()),"RECREATE");
  for(Int_t iCutSet=0; iCutSet<set->GetEntries(); ++iCutSet) {
    for(UInt_t iOrigin=kPrompt; iOrigin<=kFD; ++iOrigin)
      hEff[iOrigin][iCutSet]->Write();
  }
  outfile.Close();
  for (Int_t iCutSet=0; iCutSet<set->GetEntries(); ++iCutSet) {
    for(UInt_t iOrigin=kPrompt; iOrigin<=kFD; ++iOrigin)
      delete hEff[iOrigin][iCutSet];
  }
  delete[] hEff;
}


//_________________________________________________________________________________________________
void AliHFCutVarFDsubAnalysisManager::DrawEfficiencies(TString strOutputFolder,
                                                       TString prefix/*="eff"*/,
                                                       UInt_t xAxis/*=0*/) {
  /// Drawing efficiency plots (prompt and feed-down combined)

  // xAxis: 0 = vs Cut set, 1 = vs pt

  TList* effListPrompt = (xAxis==0) ? fEffListVsCutSets[kPrompt] : fEffListVsBins[kPrompt];
  TList* effListFD     = (xAxis==0) ? fEffListVsCutSets[kFD]     : fEffListVsBins[kFD];

  TCanvas* c = new TCanvas("cEff", "", 1920, 1080);
  if (!effListPrompt || !effListFD) {
    std::cerr << "Determine efficiencies before trying to plot them" << std::endl;
  }

  if (effListPrompt->GetEntries() != effListFD->GetEntries()) {
    std::cerr << "Size of the efficiency lists don't match! Not drawing anything." << std::endl;
  }
  for (Int_t iEff=0; iEff<effListPrompt->GetEntries(); ++iEff) {
    c->Clear();
    TLegend* l = new TLegend(0.75,0.75,0.9,0.9);
    TH1* hPrompt = (TH1*)effListPrompt->At(iEff);
    TH1* hFD     = (TH1*)effListFD->At(iEff);
    l->AddEntry(hPrompt, "Prompt", "l");
    l->AddEntry(hFD,     "FD",     "l");
    hPrompt->SetStats(kFALSE);
    hFD->SetStats(kFALSE);
    hFD->SetLineColor(kRed);
    hPrompt->Draw();
    hFD->Draw("same");
    Double_t max = (hFD->GetMaximum()>hPrompt->GetMaximum()) ?
      hFD->GetMaximum() : hPrompt->GetMaximum();
    max *= 1.1;
    hPrompt->GetYaxis()->SetRangeUser(0, max);
    l->Draw();
    c->Update();

    c->SaveAs(Form("%s/%s_%d.pdf",  strOutputFolder.Data(), prefix.Data(), iEff));
    c->SaveAs(Form("%s/%s_%d.root", strOutputFolder.Data(), prefix.Data(), iEff));
  }
  delete c;
  c = 0x0;
}


//_________________________________________________________________________________________________
void AliHFCutVarFDsubAnalysisManager::GetRawYields(Bool_t drawFit/*=kFALSE*/,TString strOutputFolder/*="."*/) {
  /// Obtain the raw yields from the THnSparse

  Int_t* nCutSets = new Int_t[fCuts->GetEntries()];
  Int_t nMaxCutSets = 0;

  for (Int_t iBin=0; iBin<fCuts->GetEntries(); ++iBin) nCutSets[iBin] = 0;

  std::vector< std::vector<Double_t> > rawYield;
  std::vector< std::vector<Double_t> > rawYieldErr;
  rawYield.resize(fCuts->GetEntries());
  rawYieldErr.resize(fCuts->GetEntries());

  if (!fBinsX) GetXaxisInformation();
  TList* set = (TList*)fCuts->At(0);
  
  TH1F** hRawYields = new TH1F*[set->GetEntries()];
  TH1F** hRawYieldSigmas = new TH1F*[set->GetEntries()];
  TH1F** hRawYieldMeans = new TH1F*[set->GetEntries()];
  TH1F** hRawYieldChisq = new TH1F*[set->GetEntries()];
  for(Int_t iCutSet=0; iCutSet<set->GetEntries(); ++iCutSet) {
    hRawYields[iCutSet] = new TH1F(Form("hRawYields_Set%d",iCutSet+1),"",fCuts->GetEntries(),fBinsX);
    hRawYieldSigmas[iCutSet] = new TH1F(Form("hRawYieldSigmas_Set%d",iCutSet+1),"",fCuts->GetEntries(),fBinsX);
    hRawYieldMeans[iCutSet] = new TH1F(Form("hRawYieldMeans_Set%d",iCutSet+1),"",fCuts->GetEntries(),fBinsX);
    hRawYieldChisq[iCutSet] = new TH1F(Form("hRawYieldChisq_Set%d",iCutSet+1),"",fCuts->GetEntries(),fBinsX);
  }
  
  TCanvas *cFitter = new TCanvas("cFitter","",1920,1080);

  for (Int_t iBin=0; iBin<fCuts->GetEntries(); ++iBin) {
    TList* cutSets = (TList*)fCuts->At(iBin);
    if (nCutSets[iBin]<cutSets->GetEntries()) nCutSets[iBin] = cutSets->GetEntries();
    if (nCutSets[iBin]>nMaxCutSets) nMaxCutSets=nCutSets[iBin];
    AliHFMassFitter** fitter = new AliHFMassFitter*[cutSets->GetEntries()];
    if(drawFit) {
      cFitter->Clear();
      cFitter->Divide(nCutSets[iBin],1);
    }
    for (Int_t iCutSet=0; iCutSet<cutSets->GetEntries(); ++iCutSet) {
      AliHFCutVarFDsubCutSet* cutSet = (AliHFCutVarFDsubCutSet*)cutSets->At(iCutSet);

      AliHFCutVarFDsubMassFitter* massFitter =
        new AliHFCutVarFDsubMassFitter(fData, fAxes, cutSet);

      rawYield[iBin].push_back(massFitter->GetSig());
      rawYieldErr[iBin].push_back(massFitter->GetSigErr());

      hRawYields[iCutSet]->SetBinContent(iBin+1,massFitter->GetSig());
      hRawYields[iCutSet]->SetBinError(iBin+1,massFitter->GetSigErr());
      hRawYieldSigmas[iCutSet]->SetBinContent(iBin+1,massFitter->GetSigma());
      hRawYieldSigmas[iCutSet]->SetBinError(iBin+1,massFitter->GetSigmaErr());
      hRawYieldMeans[iCutSet]->SetBinContent(iBin+1,massFitter->GetMean());
      hRawYieldMeans[iCutSet]->SetBinError(iBin+1,massFitter->GetMeanErr());
      hRawYieldChisq[iCutSet]->SetBinContent(iBin+1,massFitter->GetRedChiSquare());

        if(drawFit) {
        fitter[iCutSet]=(AliHFMassFitter*)massFitter->GetFitter();
        cFitter->cd(iCutSet+1);
        fitter[iCutSet]->DrawHere(gPad);
      }

      delete massFitter;
      massFitter = 0x0;
    }

    if(drawFit) {
      cFitter->SaveAs(Form("%s/Fitter_%d.pdf",strOutputFolder.Data(),iBin));
      cFitter->SaveAs(Form("%s/Fitter_%d.root",strOutputFolder.Data(),iBin));
      cFitter->Update();
    }
    delete fitter;
    fitter = 0x0;
  }
  if (fRawYields){
    delete fRawYields;
    fRawYields = 0x0;
  }
  fRawYields = new TList();
  fRawYields->SetOwner();
  for (Int_t iBin=0; iBin<fCuts->GetEntries(); ++iBin) {
    TH1F* h = new TH1F(Form("hRawYield_%d", iBin), ";Cut Set;Raw yield (a.u.)",
                       rawYield[iBin].size(), 0, (Double_t)rawYield[iBin].size());
    for (UInt_t iCutSet=0; iCutSet<rawYield[iBin].size(); ++iCutSet) {
      h->SetBinContent(iCutSet+1, rawYield[iBin][iCutSet]   );
      h->SetBinError(  iCutSet+1, rawYieldErr[iBin][iCutSet]);
    }
    fRawYields->Add((TObject*)h);
  }

  TFile outfile(Form("%s/RawYields.root",strOutputFolder.Data()),"RECREATE");
  for (Int_t iCutSet=0; iCutSet<set->GetEntries(); ++iCutSet) {
    hRawYields[iCutSet]->Write();
    hRawYieldSigmas[iCutSet]->Write();
    hRawYieldMeans[iCutSet]->Write();
    hRawYieldChisq[iCutSet]->Write();
  }
  outfile.Close();
  for (Int_t iCutSet=0; iCutSet<set->GetEntries(); ++iCutSet) {
    delete hRawYields[iCutSet];
    delete hRawYieldSigmas[iCutSet];
    delete hRawYieldMeans[iCutSet];
    delete hRawYieldChisq[iCutSet];
  }
  delete[] hRawYields;
  delete[] hRawYieldSigmas;
  delete[] hRawYieldMeans;
  delete[] hRawYieldChisq;

  delete cFitter;
  cFitter=0x0;
}


//_________________________________________________________________________________________________
void AliHFCutVarFDsubAnalysisManager::GetXaxisInformation() {
  /// Obtain x axis information

  if (fBinsX) {
    delete[] fBinsX;
    fBinsX = 0x0;
  }
  fBinsX = new Double_t[fCuts->GetEntries()+1];
  Bool_t xContinuous = kTRUE;
  for (Int_t iBin=0; iBin<fCuts->GetEntries(); ++iBin) {
    TList* cutSets = (TList*)fCuts->At(iBin);
    AliHFCutVarFDsubCutSet* cutSet = (AliHFCutVarFDsubCutSet*)cutSets->At(0);
    AliHFCutVarFDsubCut* firstCut = (AliHFCutVarFDsubCut*)cutSet->GetCut(0);
    if (iBin>0) {
      if (TMath::Abs((fBinsX[iBin]-firstCut->fLow)/(fBinsX[iBin]+1.e-3))>1.e-3) xContinuous = kFALSE;
    }
    fBinsX[iBin]   = firstCut->fLow;
    fBinsX[iBin+1] = firstCut->fHigh;
  }
  if (!xContinuous) {
    Printf("x-axis bins not continuous, using bin numbers instead of values!");
    for (Int_t iBin=0; iBin<fCuts->GetEntries()+1; ++iBin) {
      fBinsX[iBin] = (Double_t)iBin;
    }
  }
}


//_________________________________________________________________________________________________
Bool_t AliHFCutVarFDsubAnalysisManager::Minimise(UInt_t method/*=0*/,
                                                 UInt_t nIterations/*=10*/,
                                                 Bool_t useWeights/*=kTRUE*/,
                                                 Double_t relSysteEffErr/*=0.*/,
                                                 Int_t nSim/*=1000*/) {
  /// Obtain the corrected yields

  if (fCorrYieldPrompt) { delete fCorrYieldPrompt; fCorrYieldPrompt = 0x0; }
  if (fCorrYieldFD)     { delete fCorrYieldFD;     fCorrYieldFD     = 0x0; }
  if (fResiduals)       { delete fResiduals;       fResiduals       = 0x0; }
  if (fPulls)           { delete fPulls;           fPulls           = 0x0; }
  if (fIncDistError)    { delete fIncDistError;    fIncDistError    = 0x0; }
  if (fFprompt)         { delete fFprompt;         fFprompt         = 0x0; }
  if (fFpromptRaw)      { delete fFpromptRaw;      fFpromptRaw      = 0x0; }

  fCorrYieldPrompt = new TH1F("hCorrYieldPrompt", Form(";%s;Prompt D^{0} corrected yield #frac{d#it{N}_{Prompt}}{d#it{p}_{T}} (GeV^{-1}#it{c})", fxAxisTitle.Data()), fCuts->GetEntries(), fBinsX);
  fCorrYieldFD =     new TH1F("hCorrYieldFD",     Form(";%s;FD D^{0} corrected yield #frac{d#it{N}_{FD}}{d#it{p}_{T}} (GeV^{-1}#it{c})", fxAxisTitle.Data()), fCuts->GetEntries(), fBinsX);
  fResiduals  =   new TList();
  fPulls      =   new TList();
  fFprompt    =   new TList();
  fFpromptRaw =   new TList();
  fIncDistError = new TList();

  Bool_t result = kTRUE;

  for (Int_t iBin=0; iBin<fCuts->GetEntries(); ++iBin) {
    AliHFCutVarFDsubMinimiser* min
      = new AliHFCutVarFDsubMinimiser((TH1F*)fRawYields->At(iBin),
                                      (TH1F*)fEffListVsCutSets[kPrompt]->At(iBin),
                                      (TH1F*)fEffListVsCutSets[kFD]->At(iBin),
                                      method, nIterations, useWeights, relSysteEffErr, nSim);

    Bool_t minimised = min->GetStatus();
    if (!minimised) result = kFALSE;

    fCorrYieldPrompt->SetBinContent(iBin+1, (minimised) ? min->GetPromptYield()    : 0.);
    fCorrYieldPrompt->SetBinError(  iBin+1, (minimised) ? min->GetPromptYieldErr() : 0.);
    fCorrYieldFD->SetBinContent(iBin+1, (minimised) ? min->GetFDYield()    : 0.);
    fCorrYieldFD->SetBinError(  iBin+1, (minimised) ? min->GetFDYieldErr() : 0.);
    
    fResiduals->Add( (minimised) ? min->GetResiduals()->Clone()  : (TObject*)0x0);
    fPulls->Add(     (minimised) ? min->GetPulls()->Clone()      : (TObject*)0x0);
    
    if(method==1) {
      fIncDistError->Add( (minimised) ? min->GetIncDistError()->Clone() : (TObject*)0x0);
    }
    
    fFprompt->Add(   (minimised) ? min->GetFprompt()->Clone()    : (TObject*)0x0);
    fFpromptRaw->Add((minimised) ? min->GetFpromptRaw()->Clone() : (TObject*)0x0);

    delete min;
    min = 0x0;
  }
  return result;
}


//_________________________________________________________________________________________________
void AliHFCutVarFDsubAnalysisManager::DrawLines(TString strOutputFolder, Bool_t IsIncentre) {

  TList* cutSets = (TList*)fCuts->At(0);
  Int_t nCutSets = cutSets->GetEntries();

  TF1 **line = new TF1*[nCutSets];
  TF1 **lineErrLow = new TF1*[nCutSets];
  TF1 **lineErrHigh = new TF1*[nCutSets];
  const char *formula = "([2]-[1]*x)/[0]";//[0]->EffPrompt, [1]->EffFD, [2]->RawYield

  Double_t xmin = -1e+5;
  Double_t xmax = 1e+5;
  Double_t ymin = -1e+5;
  Double_t ymax = 1e+5;

  if (!fEffListVsCutSets[kPrompt] || !fEffListVsCutSets[kFD]) {
    std::cerr << "Determine efficiencies before trying to plot them" << std::endl;
  }

  if(!fRawYields) {
    std::cerr << "Determine raw yields before trying to plot them" << std::endl;
  }

  if (fEffListVsCutSets[kPrompt]->GetEntries() != fEffListVsCutSets[kFD]->GetEntries()) {
    std::cerr << "Size of the efficiency lists don't match! Not drawing anything." << std::endl;
  }

  if (fEffListVsCutSets[kPrompt]->GetEntries() != fRawYields->GetEntries()) {
    std::cerr << "Size of the efficiency lists and the raw yield list don't match! Not drawing anything." << std::endl;
  }

  Color_t color[10] = {kRed,kBlue,kGreen+3,kOrange+7,kBlack,kMagenta,kCyan,kOrange+3,kGray+2,kViolet+3};//10 is the max number of cut sets

  TCanvas *cLines = new TCanvas("cLines","", 1920, 1080);

  for (Int_t iBin=0; iBin<fCuts->GetEntries(); ++iBin) {
    cLines->Clear();

    TLegend* l = new TLegend(0.75,0.75,0.9,0.9);

    TH1F *hEffPrompt = (TH1F*)fEffListVsCutSets[kPrompt]->At(iBin);
    TH1F *hEffFD = (TH1F*)fEffListVsCutSets[kFD]->At(iBin);
    TH1F *hRawYields = (TH1F*)fRawYields->At(iBin);

    Double_t xmin = -1e+5;
    Double_t xmax = 1e+5;
    Double_t ymin = -1e+5;
    Double_t ymax = 1e+5;

    if(fCorrYieldFD!=0x0 && fCorrYieldFD->GetBinContent(iBin+1)!=0) {
      xmin = fCorrYieldFD->GetBinContent(iBin+1)-4*TMath::Abs(fCorrYieldFD->GetBinContent(iBin+1));
      xmax = fCorrYieldFD->GetBinContent(iBin+1)+4*TMath::Abs(fCorrYieldFD->GetBinContent(iBin+1));
    }

    if(fCorrYieldPrompt!=0x0 && fCorrYieldPrompt->GetBinContent(iBin+1)!=0) {
      ymin = fCorrYieldPrompt->GetBinContent(iBin+1)-4*TMath::Abs(fCorrYieldPrompt->GetBinContent(iBin+1));
      ymax = fCorrYieldPrompt->GetBinContent(iBin+1)+4*TMath::Abs(fCorrYieldPrompt->GetBinContent(iBin+1));
    }

    for(Int_t iCutSet=0; iCutSet<nCutSets; ++iCutSet) {

      line[iCutSet] = new TF1(Form("lineSet%d_PtBin%d",iCutSet+1,iBin),formula,xmin,xmax);
      line[iCutSet]->SetParameters(hEffPrompt->GetBinContent(iCutSet+1),
                                   hEffFD->GetBinContent(iCutSet+1),
                                   hRawYields->GetBinContent(iCutSet+1));
      line[iCutSet]->SetLineColor(color[iCutSet]);
      line[iCutSet]->SetLineWidth(1);
      l->AddEntry(line[iCutSet],Form("Set %d",iCutSet+1),"l");

      lineErrLow[iCutSet] = new TF1(Form("lineErrLowSet%d_PtBin%d",iCutSet+1,iBin),formula,xmin,xmax);
      lineErrLow[iCutSet]->SetParameters(hEffPrompt->GetBinContent(iCutSet+1)+hEffPrompt->GetBinError(iCutSet+1),
                                         hEffFD->GetBinContent(iCutSet+1)+hEffFD->GetBinError(iCutSet+1),
                                         hRawYields->GetBinContent(iCutSet+1)-hRawYields->GetBinError(iCutSet+1));
      lineErrLow[iCutSet]->SetLineColor(color[iCutSet]);
      lineErrLow[iCutSet]->SetLineWidth(1);
      lineErrLow[iCutSet]->SetLineStyle(7);

      lineErrHigh[iCutSet] = new TF1(Form("lineErrHighSet%d_PtBin%d",iCutSet+1,iBin),formula,xmin,xmax);
      lineErrHigh[iCutSet]->SetParameters(hEffPrompt->GetBinContent(iCutSet+1)-hEffPrompt->GetBinError(iCutSet+1),
                                          hEffFD->GetBinContent(iCutSet+1)+-hEffFD->GetBinError(iCutSet+1),
                                          hRawYields->GetBinContent(iCutSet+1)+hRawYields->GetBinError(iCutSet+1));
      lineErrHigh[iCutSet]->SetLineColor(color[iCutSet]);
      lineErrHigh[iCutSet]->SetLineWidth(1);
      lineErrHigh[iCutSet]->SetLineStyle(7);

      line[iCutSet]->GetXaxis()->SetTitle("N_{FD}");
      line[iCutSet]->GetYaxis()->SetTitle("N_{Prompt}");

      if(iCutSet==0) {
        line[iCutSet]->SetTitle("");
        line[iCutSet]->GetYaxis()->SetRangeUser(ymin,ymax);
        line[iCutSet]->Draw();
      }
      else
        line[iCutSet]->Draw("same");

      lineErrLow[iCutSet]->Draw("same");
      lineErrHigh[iCutSet]->Draw("same");
    }

    l->Draw("same");

    if(IsIncentre) {
      TH2F *hIncDist = (TH2F*)fIncDistError->At(iBin);
      hIncDist->SetStats(kFALSE);
      if(hIncDist)
        hIncDist->Draw("same");
    }

    cLines->SaveAs(Form("%s/Lines_%d.pdf",  strOutputFolder.Data(), iBin));
    cLines->SaveAs(Form("%s/Lines_%d.root", strOutputFolder.Data(), iBin));
    cLines->Update();

    delete l;
    l=0x0;
  }

  for(Int_t iCutSet=0; iCutSet<nCutSets; ++iCutSet) {
    delete line[iCutSet];
    line[iCutSet] = 0x0;
    delete line[iCutSet];
    lineErrLow[iCutSet] = 0x0;
    delete lineErrHigh[iCutSet];
    lineErrHigh[iCutSet] = 0x0;
  }

  delete cLines;
  cLines=0x0;
  delete[] line;
  line=0x0;
  delete[] lineErrLow;
  lineErrLow=0x0;
  delete[] lineErrHigh;
  lineErrHigh=0x0;
}


//_________________________________________________________________________________________________
AliHFCutVarFDsubAnalysisManager::AliHFCutVarFDsubAnalysisManager(const AliHFCutVarFDsubAnalysisManager& am)
  : TObject()
  , fData(am.fData)
  , fAxes(am.fAxes)
  , fCuts(am.fCuts)
  , fRawYields(am.fRawYields)
  , fCorrYieldPrompt(am.fCorrYieldPrompt)
  , fCorrYieldFD(am.fCorrYieldFD)
  , fResiduals(am.fResiduals)
  , fPulls(am.fPulls)
  , fFprompt(am.fFprompt)
  , fFpromptRaw(am.fFpromptRaw)
  , fxAxisTitle(am.fxAxisTitle)
  , fBinsX(am.fBinsX)
{
  /// Copy constructor

  // TODO: do proper deep copies
  Printf("Do not use this copy constructor!");
}

//_________________________________________________________________________________________________
AliHFCutVarFDsubAnalysisManager AliHFCutVarFDsubAnalysisManager::operator=(const AliHFCutVarFDsubAnalysisManager& am)
{
  /// Assignment operator

  // TODO: do proper deep copies
  Printf("Do not use this assignment operator!");

  if (this != &am) {
    fData = am.fData;
    fAxes = am.fAxes;
    fCuts = am.fCuts;
    fRawYields = am.fRawYields;
    fCorrYieldPrompt = am.fCorrYieldPrompt;
    fCorrYieldFD = am.fCorrYieldFD;
    fResiduals = am.fResiduals;
    fPulls = am.fPulls;
    fFprompt = am.fFprompt;
    fFpromptRaw = am.fFpromptRaw;
    fxAxisTitle = am.fxAxisTitle;
    fBinsX = am.fBinsX;
  }
  return *this;
}
