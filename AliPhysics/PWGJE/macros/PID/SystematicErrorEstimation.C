#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"

#include "AliPID.h"

#include <iostream>
#include <iomanip>

#include "SystematicErrorUtils.h"

const Int_t numSpecies = AliPID::kSPECIES;

enum resultType { kFraction = 0, kYield = 1, kToPiRatio = 2 };

//________________________________________________________
TCanvas* calculateSystematics(TString canvName, TString canvTitle, TH1F** histos, Int_t numHistos, Int_t speciesID, Double_t /*nSigma*/,
                              const TString* systematicsHistosName, Int_t reference, TH1F** hSystematics, TGraphAsymmErrors** gr,
                              Bool_t ignoreSigmaErrors, Bool_t setMean, Double_t setMeanLowerThreshold, Double_t setMeanUpperThreshold, 
                              Bool_t toPiRatio = kFALSE)
{
  // For every bin:
  // Since the method with the root finding already takes into account the statistical error,
  // there is no need to use nSigma > 0.
  // If the statistical error is ignored, nevertheless don't use nSigma > 0 because this might
  // give zero systematic error for high pT, which is usually not accepted by people, although
  // the natural point of view "no systematic visible for given statistical error" is reasonable to me.
  
  if (!histos || !hSystematics || !gr)
    return 0x0;
  
  Double_t ymax = 0;
  Double_t ymin = 0;
  
  
  // Just for drawing
  for (Int_t j = 0; j < numHistos; j++) {
    hSystematics[j] = new TH1F(*histos[j]);
    hSystematics[j]->SetName(Form("%s%s_%s", systematicsHistosName[j].Data(), toPiRatio ? "toPiRatio" : "",
                                  AliPID::ParticleName(speciesID)));
    hSystematics[j]->Reset(); 
    hSystematics[j]->GetXaxis()->SetRange(0, -1);
    
    for (Int_t bin = 1; bin <= histos[j]->GetNbinsX(); bin++) {
      hSystematics[j]->SetBinContent(bin, histos[reference]->GetBinContent(bin) - histos[j]->GetBinContent(bin));
      hSystematics[j]->SetBinError(bin, TMath::Sqrt(TMath::Abs(TMath::Power(histos[reference]->GetBinError(bin), 2) - 
                                                               TMath::Power(histos[j]->GetBinError(bin), 2))));
  
      if (hSystematics[j]->GetBinError(bin) == 0)
        hSystematics[j]->SetBinError(bin, 1e-10);
      Double_t temp = hSystematics[j]->GetBinContent(bin) + hSystematics[j]->GetBinError(bin);
      if (temp > ymax)
        ymax = temp;
        
      temp = hSystematics[j]->GetBinContent(bin) - hSystematics[j]->GetBinError(bin);
      if (temp < ymin)
        ymin = temp;
    }
  }
  
  TCanvas* canv = new TCanvas(canvName.Data(), canvTitle.Data(),100,10,1200,800);
  canv->SetGridy(1);
  
  hSystematics[reference]->Draw("e p");
  hSystematics[reference]->GetYaxis()->SetRangeUser(ymin, ymax);
  for (Int_t j = 0; j < numHistos; j++) {
    if (j == reference)
      continue;
  
    hSystematics[j]->SetMarkerStyle(20 + j);
    hSystematics[j]->Draw("e p same");
  }
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  
  for (Int_t j = 0; j < numHistos; j++) {
    legend->AddEntry(hSystematics[j], Form("%s", systematicsHistosName[j].Data()), "p");
  }
  legend->Draw();
             
             
  const Int_t nBins = histos[reference]->GetNbinsX();
  Double_t x[nBins];
  Double_t y[nBins];
  Double_t xerr[nBins];
  Double_t yerrl[nBins];
  Double_t yerrh[nBins];
  
  Double_t meansForFit[numHistos];
  Double_t sigmasForFit[numHistos];
  
  for (Int_t bin = 0; bin < nBins; bin++) {
    x[bin] = histos[reference]->GetBinCenter(bin + 1);
    xerr[bin] = histos[reference]->GetBinWidth(bin + 1) / 2.;
    y[bin] = histos[reference]->GetBinContent(bin + 1);
    
    Double_t ymaxOriginal = -1;
    Double_t yminOriginal = -1;
    Bool_t problem = kFALSE;
    
    for (Int_t j = 0; j < numHistos; j++) {
      meansForFit[j] = histos[j]->GetBinContent(bin + 1);
      sigmasForFit[j] = histos[j]->GetBinError(bin + 1);
      
      if (toPiRatio && (meansForFit[j] < 0 || meansForFit[j] > 10))
        problem = kTRUE;
      
      Double_t tempOriginal = histos[j]->GetBinContent(bin + 1);
      if (tempOriginal > ymaxOriginal || ymaxOriginal < 0)
        ymaxOriginal = tempOriginal;
      if (tempOriginal < yminOriginal || yminOriginal < 0)
        yminOriginal = tempOriginal;
    }
    
    Double_t weightedMean = 0;
    
    if (problem) {
      weightedMean = y[bin];
      yerrl[bin] = yerrh[bin] = 1.;
      
      printf("Found problem in bin %d (x = %e) (toPiRatio %d), %s\n", bin, histos[reference]->GetXaxis()->GetBinCenter(bin), toPiRatio,
        canvName.Data());
      for (Int_t j = 0; j < numHistos; j++) 
        printf("meansForFit[%d] = %e\n", j, meansForFit[j]);
      printf("weightedMean set to ref value, error set to 1.\n\n");
    }
    // The root finder needs some given range. For yields, the error should lie within [0, maxMean - minMean]
    yerrl[bin] = yerrh[bin] = findSystematicError(numHistos, meansForFit, sigmasForFit, ignoreSigmaErrors, 0, ymaxOriginal - yminOriginal,
                                                  yminOriginal, ymaxOriginal, &weightedMean);
    
    // If the mean should be set to the weighted mean and the bin is in the desired range, do it for the reference histo and the graph!
    // NOTE: The reference histo is only changed for the summary plot, but it is not saved separately!
    if (setMean && x[bin] >= setMeanLowerThreshold && x[bin] <= setMeanUpperThreshold) {
      histos[reference]->SetBinContent(bin + 1, weightedMean);
      y[bin] = weightedMean;
    }
  }
  
  TGraphAsymmErrors* gTemp = new TGraphAsymmErrors(nBins, x, y, xerr, xerr, yerrl, yerrh);
  *gr = gTemp;
  (*gr)->SetName(Form("systematicError%s_%s", toPiRatio ? "ToPiRatio" : "", AliPID::ParticleName(speciesID)));
  (*gr)->SetLineColor(hSystematics[0]->GetMarkerColor());
  //(*gr)->SetFillColor(kGray);
  (*gr)->SetFillStyle(0);//3004 + reference);
  
  return canv;
}


/*OLD
//________________________________________________________
TCanvas* calculateSystematics(TString canvName, TString canvTitle, TH1F** histos, Int_t numHistos, Int_t speciesID, Double_t nSigma,
                              const TString* systematicsHistosName, Int_t reference, TH1F** hSystematics, TGraphAsymmErrors** gr)
{
  Double_t ymax = 0;
  Double_t ymin = 0;
  
  for (Int_t j = 0; j < numHistos; j++) {
    hSystematics[j] = new TH1F(*histos[j]);
    hSystematics[j]->SetName(Form("%s_%s", systematicsHistosName[j].Data(), AliPID::ParticleName(speciesID)));
    hSystematics[j]->Reset(); 
    hSystematics[j]->GetXaxis()->SetRange(0, -1);
    
    for (Int_t bin = 1; bin <= histos[j]->GetNbinsX(); bin++) {
      hSystematics[j]->SetBinContent(bin, histos[reference]->GetBinContent(bin) - histos[j]->GetBinContent(bin));
      hSystematics[j]->SetBinError(bin, TMath::Sqrt(TMath::Abs(TMath::Power(histos[reference]->GetBinError(bin), 2) - 
                                                               TMath::Power(histos[j]->GetBinError(bin), 2))));
  
      if (hSystematics[j]->GetBinError(bin) == 0)
        hSystematics[j]->SetBinError(bin, 1e-10);
      Double_t temp = hSystematics[j]->GetBinContent(bin) + hSystematics[j]->GetBinError(bin);
      if (temp > ymax)
        ymax = temp;
        
      temp = hSystematics[j]->GetBinContent(bin) - hSystematics[j]->GetBinError(bin);
      if (temp < ymin)
        ymin = temp;
    }
  }
  
  TCanvas* canv = new TCanvas(canvName.Data(), canvTitle.Data(),100,10,1200,800);
  canv->SetGridy(1);
  
  hSystematics[reference]->Draw("e p");
  hSystematics[reference]->GetYaxis()->SetRangeUser(ymin, ymax);
  for (Int_t j = 0; j < numHistos; j++) {
    if (j == reference)
      continue;
  
    hSystematics[j]->SetMarkerStyle(20 + j);
    hSystematics[j]->Draw("e p same");
  }
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  
  for (Int_t j = 0; j < numHistos; j++) {
    legend->AddEntry(hSystematics[j], Form("%s", systematicsHistosName[j].Data()), "p");
  }
  legend->Draw();
			       
			       
  const Int_t nBins = histos[reference]->GetNbinsX();
  Double_t x[nBins];
  Double_t y[nBins];
  Double_t xerr[nBins];
  Double_t yerrl[nBins];
  Double_t yerrh[nBins];
  
  for (Int_t bin = 0; bin < nBins; bin++) {
    x[bin] = histos[reference]->GetBinCenter(bin + 1);
    xerr[bin] = histos[reference]->GetBinWidth(bin + 1) / 2.;
    y[bin] = histos[reference]->GetBinContent(bin + 1);
    
    // Take all points that are more than nSigma sigma away from 0.
    // If there are at least 2 such points, take the difference between
    // the extreme values (i.e. maximum and minimum) as a measure of
    // the systematics
    Int_t count = 0;
    Double_t deltaMin = 0;
    Double_t deltaMax = 0;
    
    for (Int_t j = 0; j < numHistos; j++) {
      if (hSystematics[j]->GetBinError(bin + 1) == 0) // Per definition always true for reference histo
        continue;
        
      Double_t delta = hSystematics[j]->GetBinContent(bin + 1);
      if (TMath::Abs(delta / hSystematics[j]->GetBinError(bin + 1)) > nSigma)  {
        //if (count == 0) {
        //  deltaMin = delta;
        //  deltaMax = delta;
        //}
        //else {
          if (delta < deltaMin)
            deltaMin = delta;
          if (delta > deltaMax)
            deltaMax = delta;
        //}
        count++;
      }
    }
    
    //if (deltaMax > 0.) 
    //  yerrh[bin] = deltaMax;
    //else
    //  yerrh[bin] = 0.;
    //  
    //if (deltaMin < 0.)
    //  yerrl[bin] = -deltaMin;
    //else
    //  yerrl[bin] = 0.;
    
    if (count < 1) // Reference histo is not counted. One can only do systematics if there is at least one other histogram
      yerrl[bin] = yerrh[bin] = 0.;
    else
      yerrl[bin] = yerrh[bin] = (deltaMax - deltaMin) / TMath::Sqrt(2);
    
  }
  
  TGraphAsymmErrors* gTemp = new TGraphAsymmErrors(nBins, x, y, xerr, xerr, yerrl, yerrh);
  *gr = gTemp;
  (*gr)->SetName(Form("systematicError_%s", AliPID::ParticleName(speciesID)));
  (*gr)->SetLineColor(hSystematics[0]->GetMarkerColor());
  //(*gr)->SetFillColor(kGray);
  (*gr)->SetFillStyle(0);//3004 + reference);
  
  return canv;
}*/


//________________________________________________________
TCanvas* DrawResultHistos(TString canvName, TString canvTitle, Double_t pLow, Double_t pHigh, TH1F*** hist, Int_t reference,
                          TGraphAsymmErrors** gr, Int_t resultType, Int_t mode = -1)
{
  if (!hist || !gr)
    return 0x0;
  
  TCanvas* canv = new TCanvas(canvName.Data(), canvTitle.Data(),100,10,1200,800);
  canv->SetGridx(1);
  canv->SetGridy(1);
  if (mode <= 0)
    canv->SetLogx(1);
  for (Int_t i = 0; i < numSpecies; i++) {
    if (!hist[i][reference])
      continue;
    
    if (i == AliPID::kPion && resultType == kToPiRatio)
      continue;
    
    if (resultType == kYield) {
      if (mode < 0)
        hist[i][reference]->GetYaxis()->SetRangeUser(2e-12, 18);
      else  
        hist[i][reference]->GetYaxis()->SetRangeUser(2e-7, 6);
    }
    else
      hist[i][reference]->GetYaxis()->SetRangeUser(0.0, 1.0);
    //hist[i][reference]->GetYaxis()->SetTitle(canvTitle.Data());
    if (mode <= 0)
      hist[i][reference]->GetXaxis()->SetRangeUser(pLow, pHigh);
    //hist[i][reference]->SetFillStyle(3004 + i);
    //hist[i][reference]->SetFillColor(kGray);
    hist[i][reference]->SetFillStyle(0);
    hist[i][reference]->SetFillColor(hist[i][reference]->GetMarkerColor());
    hist[i][reference]->SetLineColor(hist[i][reference]->GetMarkerColor());
  }
  if (resultType != kToPiRatio) {
    hist[AliPID::kPion][reference]->SetMarkerStyle(20);
    hist[AliPID::kPion][reference]->Draw("e p");
  }
  hist[AliPID::kElectron][reference]->SetMarkerStyle(21);
  hist[AliPID::kElectron][reference]->Draw(Form("e p%s", (resultType == kToPiRatio) ? "" : " same"));
  hist[AliPID::kKaon][reference]->SetMarkerStyle(22);
  hist[AliPID::kKaon][reference]->Draw("e p same");
  hist[AliPID::kProton][reference]->SetMarkerStyle(29);
  hist[AliPID::kProton][reference]->Draw("e p same");
  hist[AliPID::kMuon][reference]->SetMarkerStyle(30);
  hist[AliPID::kMuon][reference]->Draw("e p same");
  
  gr[AliPID::kElectron]->Draw("2 same");
  gr[AliPID::kKaon]->Draw("2 same");
  if (resultType != kToPiRatio)
    gr[AliPID::kPion]->Draw("2 same");
  gr[AliPID::kProton]->Draw("2 same");
  gr[AliPID::kMuon]->Draw("2 same");
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  
  TString suffix = (resultType == kToPiRatio) ? "/#pi" : "";
  if (resultType != kToPiRatio)
    legend->AddEntry(hist[AliPID::kPion][reference], Form("#pi%s", suffix.Data()), "flp");
  legend->AddEntry(hist[AliPID::kElectron][reference], Form("e%s", suffix.Data()), "flp");
  legend->AddEntry(hist[AliPID::kKaon][reference], Form("K%s", suffix.Data()), "flp");
  legend->AddEntry(hist[AliPID::kProton][reference], Form("p%s", suffix.Data()), "flp");
  legend->AddEntry(hist[AliPID::kMuon][reference], Form("#mu%s", suffix.Data()), "flp");
  legend->Draw();
  
  return canv;
}


//________________________________________________________
TH1F* loadHisto(const TString histName, TFile* f, Bool_t warningIfHistoNotExists = kTRUE)
{
  if (!f) {
    std::cout << "No file. Cannot load hist \"" << histName.Data() << "\" from file \"" << f->GetName() << "\"!" << std::endl;
    return 0x0;
  }
  
  TH1F* hTemp = (TH1F*)(f->Get(histName.Data()));
  if (!hTemp) {
    if (warningIfHistoNotExists)
      std::cout << "Failed to load histo \"" << histName.Data() << "\" from file \"" << f->GetName() << "\"!" << std::endl;
    return 0x0;
  } 
  
  return hTemp;
}


//________________________________________________________
Int_t SystematicErrorEstimation(const TString path, const TString outFileTitle, const TString* fileNames, const TString* histTitles,
                                const Int_t numFiles, const Double_t nSigma, const Bool_t ignoreSigmaErrors, Bool_t setMean, 
                                const Double_t setMeanLowerThreshold, const Double_t setMeanUpperThreshold, const Int_t mode = -1) 
{ 
  // NOTE: Normally, weighted average should be done either for the yields or should be WEIGHTED with the total yield
  // (must be done anyway for the to-pion-ratios!). BUT: The whole approach relies on the assumption that the overall statistics
  // is the same for all compared data points. This must be checked anyway and is true within percent typically!
  if (!fileNames || numFiles < 1)
    return -1;
  
  TFile* f[numFiles];
  TH1F** hFractions[numSpecies];
  for (Int_t i = 0; i < numSpecies; i++)  {
    hFractions[i] = new TH1F*[numFiles];
    for (Int_t j = 0; j < numFiles; j++) 
      hFractions[i][j] = 0x0;
  }
  
  TH1F** hRatioToPi[numSpecies];
  for (Int_t i = 0; i < numSpecies; i++) {
    hRatioToPi[i] = new TH1F*[numFiles];
    for (Int_t j = 0; j < numFiles; j++) 
      hRatioToPi[i][j] = 0x0;
  }
  
  const Int_t reference = 0;
  TH1F* hYields[numSpecies] = { 0x0, }; // Only the reference yields
  
  // Histos for inclusive case (mode = -1)
  TString histNames[numSpecies] = {"hFractionElectrons", "hFractionMuons", "hFractionPions", "hFractionKaons", "hFractionProtons" };
  TString histNamesYields[numSpecies] = {"hYieldElectrons", "hYieldMuons", "hYieldPions", "hYieldKaons", "hYieldProtons" };
  TString histNamesToPiRatios[numSpecies] = {"hRatioToPiElectrons", "hRatioToPiMuons", "", "hRatioToPiKaons", "hRatioToPiProtons" };
  
  // Histos for jet case (mode = 0,1,2 = pT,z,xi)
  if (mode == 0) {
    for (Int_t i = 0; i < numSpecies; i++) {
      histNames[i] = Form("hFractionIDFFtrackPt_%s", AliPID::ParticleShortName(i));
      histNamesYields[i] = Form("hIDFFtrackPt_%s", AliPID::ParticleShortName(i));
      
      if (i != AliPID::kPion)
        histNamesToPiRatios[i] = Form("hRatioToPiIDFFtrackPt_%s", AliPID::ParticleShortName(i));
    }
  }
  else if (mode == 1) {
    for (Int_t i = 0; i < numSpecies; i++) {
      histNames[i] = Form("hFractionIDFFz_%s", AliPID::ParticleShortName(i));
      histNamesYields[i] = Form("hIDFFz_%s", AliPID::ParticleShortName(i));
      
      if (i != AliPID::kPion)
        histNamesToPiRatios[i] = Form("hRatioToPiIDFFz_%s", AliPID::ParticleShortName(i));
    }
  }
  else if (mode == 2) {
    for (Int_t i = 0; i < numSpecies; i++) {
      histNames[i] = Form("hFractionIDFFxi_%s", AliPID::ParticleShortName(i));
      histNamesYields[i] = Form("hIDFFxi_%s", AliPID::ParticleShortName(i));
      
      if (i != AliPID::kPion)
        histNamesToPiRatios[i] = Form("hRatioToPiIDFFxi_%s", AliPID::ParticleShortName(i));
    }
  }
  
  
  TH1* hNumEvents = 0x0;
  TH1* hNumEventsTriggerSel = 0x0;
  TH1* hNumEventsTriggerSelVtxCut = 0x0;
  TH1* hNumEventsTriggerSelVtxCutNoPileUpRejection = 0x0;
  
  TH1* hNjetsGen = 0x0;
  TH1* hNjetsRec = 0x0;
  
  // For backward compatibility
  Bool_t hasToPiRatios = kFALSE;
  Bool_t hasYield = kFALSE;
  
  for (Int_t iFile = 0; iFile < numFiles; iFile++) {
    f[iFile] = TFile::Open(fileNames[iFile].Data());
    if (!f[iFile])  {
      std::cout << "Failed to open file \"" << fileNames[iFile].Data() << "\"!" << std::endl;
      return -1;
    }
        
    // Extract the data histograms
    for (Int_t i = 0; i < numSpecies; i++) {
      hFractions[i][iFile] = loadHisto(histNames[i], f[iFile]);
      if (!hFractions[i][iFile])
        return -1;
      
      if (histNamesToPiRatios[i] != "") {
        hRatioToPi[i][iFile] = loadHisto(histNamesToPiRatios[i], f[iFile]);
        if (!hRatioToPi[i][iFile]) {
          if (hasToPiRatios)
            return -1;
        }
        else
          hasToPiRatios = kTRUE;
      }
      
      if (iFile == reference) {
        hYields[i] = loadHisto(histNamesYields[i], f[iFile]);
        if (!hYields[i]) {
          if (hasYield)
            return -1;
        }
        else
          hasYield = kTRUE;
        
        // Load numEvent histos from reference file and store them (should ideally be the same for all,
        // but since mainly a ratio of events is used (normalisation of yields already done) for a correction later on,
        // it doesn't matter if it is not the same. What matters is that the numbers are correct for the reference!
        // NOTE: Anyway, for the sum of the errors the correct histos will be stored only in the corresponding file.
        // These histos are just, if processing of this intermediate results is done...
        hNumEvents = loadHisto("fhEventsProcessed", f[iFile]);
        hNumEventsTriggerSel = loadHisto("fhEventsTriggerSel", f[iFile]);
        hNumEventsTriggerSelVtxCut = loadHisto("fhEventsTriggerSelVtxCut", f[iFile]);
        hNumEventsTriggerSelVtxCutNoPileUpRejection = loadHisto("fhEventsProcessedNoPileUpRejection", f[iFile]);
        
        // Same for numJet histos
        hNjetsGen = loadHisto("fh2FFJetPtGen", f[iFile], kFALSE);
        hNjetsRec = loadHisto("fh2FFJetPtRec", f[iFile], kFALSE);
      }
    }
  }
  
  if (!hasToPiRatios) {
    printf("No to-pi-ratios found...\n");
  }
  
  TGraphAsymmErrors* grSysErrors[numSpecies] = {0x0,};
  TGraphAsymmErrors* grSysErrorsToPiRatios[numSpecies] = {0x0,};
  TGraphAsymmErrors* grSysErrorsYields[numSpecies] = {0x0,};
  
  // Fractions
  TH1F* hSystematicsPions[numFiles];
  TH1F* hSystematicsElectrons[numFiles];
  TH1F* hSystematicsKaons[numFiles];
  TH1F* hSystematicsProtons[numFiles]; 
  TH1F* hSystematicsMuons[numFiles];
  
  for (Int_t i = 0; i < numFiles; i++) {
    hSystematicsPions[i] = 0x0;
    hSystematicsElectrons[i] = 0x0;
    hSystematicsKaons[i] = 0x0;
    hSystematicsProtons[i] = 0x0;
    hSystematicsMuons[i] = 0x0;
  }

  TCanvas* cSystematicsPions = calculateSystematics("cSystematicsPions", "Systematics Pions", hFractions[AliPID::kPion], numFiles,
                                                    AliPID::kPion, nSigma, histTitles, reference, hSystematicsPions,
                                                    &grSysErrors[AliPID::kPion], ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                                    setMeanUpperThreshold);

  TCanvas* cSystematicsElectrons = calculateSystematics("cSystematicsElectrons", "Systematics Electrons", hFractions[AliPID::kElectron], 
                                                        numFiles, AliPID::kElectron, nSigma, histTitles, reference, hSystematicsElectrons,
                                                        &grSysErrors[AliPID::kElectron], ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                                        setMeanUpperThreshold);
  
  TCanvas* cSystematicsKaons = calculateSystematics("cSystematicsKaons", "Systematics Kaons", hFractions[AliPID::kKaon], numFiles, 
                                                    AliPID::kKaon, nSigma, histTitles, reference, hSystematicsKaons, 
                                                    &grSysErrors[AliPID::kKaon], ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                                    setMeanUpperThreshold);
  
  TCanvas* cSystematicsProtons = calculateSystematics("cSystematicsProtons", "Systematics Protons", hFractions[AliPID::kProton], numFiles,
                                                      AliPID::kProton, nSigma, histTitles, reference, hSystematicsProtons,
                                                      &grSysErrors[AliPID::kProton], ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                                      setMeanUpperThreshold);
  
  TCanvas* cSystematicsMuons = calculateSystematics("cSystematicsMuons", "Systematics Muons", hFractions[AliPID::kMuon], numFiles,
                                                    AliPID::kMuon, nSigma, histTitles, reference, hSystematicsMuons, 
                                                    &grSysErrors[AliPID::kMuon], ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                                    setMeanUpperThreshold);
  
  // To-pi-ratios
  TH1F* hSystematicsToPiRatioElectrons[numFiles];
  TH1F* hSystematicsToPiRatioKaons[numFiles];
  TH1F* hSystematicsToPiRatioProtons[numFiles]; 
  TH1F* hSystematicsToPiRatioMuons[numFiles];
  
  for (Int_t i = 0; i < numFiles; i++) {
    hSystematicsToPiRatioElectrons[i] = 0x0;
    hSystematicsToPiRatioKaons[i] = 0x0;
    hSystematicsToPiRatioProtons[i] = 0x0;
    hSystematicsToPiRatioMuons[i] = 0x0;
  }
  
  TCanvas* cSystematicsToPiRatioElectrons = calculateSystematics("cSystematicsToPiRatioElectrons", "Systematics Electrons",
                                                                 hRatioToPi[AliPID::kElectron], numFiles, AliPID::kElectron, nSigma, 
                                                                 histTitles, reference, hSystematicsToPiRatioElectrons, 
                                                                 &grSysErrorsToPiRatios[AliPID::kElectron], ignoreSigmaErrors, setMean, setMeanLowerThreshold, setMeanUpperThreshold, kTRUE);
  
  TCanvas* cSystematicsToPiRatioKaons = calculateSystematics("cSystematicsToPiRatioKaons", "Systematics Kaons", hRatioToPi[AliPID::kKaon], 
                                                             numFiles, AliPID::kKaon, nSigma, histTitles, reference, 
                                                             hSystematicsToPiRatioKaons, &grSysErrorsToPiRatios[AliPID::kKaon], 
                                                             ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                                             setMeanUpperThreshold, kTRUE);
  
  TCanvas* cSystematicsToPiRatioProtons = calculateSystematics("cSystematicsToPiRatioProtons", "Systematics Protons", 
                                                               hRatioToPi[AliPID::kProton], numFiles, AliPID::kProton, nSigma, histTitles, 
                                                               reference, hSystematicsToPiRatioProtons, 
                                                               &grSysErrorsToPiRatios[AliPID::kProton], ignoreSigmaErrors, setMean, setMeanLowerThreshold, setMeanUpperThreshold, kTRUE);
  
  TCanvas* cSystematicsToPiRatioMuons = calculateSystematics("cSystematicsToPiRatioMuons", "Systematics Muons", hRatioToPi[AliPID::kMuon], 
                                                             numFiles, AliPID::kMuon, nSigma, histTitles, reference, hSystematicsToPiRatioMuons, &grSysErrorsToPiRatios[AliPID::kMuon], 
                                                             ignoreSigmaErrors, setMean, setMeanLowerThreshold,
                                                             setMeanUpperThreshold, kTRUE);
  
  // Yields
  
  // The error of the fractions and the yield is just a constant factor (number of tracks in that bin).
  // Thus, the relative errors are the same for fractions and yields and I can just use this fact to
  // transform the errors from the fractions to those of the yields.
  // However, this causes trouble in case of fraction = 0. Therefore, sum up the yields to the total yield and use this for scaling
  // NOTE: See comment at beginning of this function concerning problems with different weights!
  
  if (hasYield) {
    for (Int_t i = 0; i < numSpecies; i++) {
      grSysErrorsYields[i] = new TGraphAsymmErrors(*grSysErrors[i]);
      TString name = grSysErrors[i]->GetName();
      name.ReplaceAll("systematicError_", "systematicErrorYields_");
      grSysErrorsYields[i]->SetName(name.Data());
      
      for (Int_t ind = 0; ind < grSysErrorsYields[i]->GetN(); ind++) {
        Double_t totalYield = 0;
        for (Int_t j = 0; j < numSpecies; j++)
          totalYield += hYields[j]->GetBinContent(ind + 1);
        
        const Double_t yield = hYields[i]->GetBinContent(ind + 1);
        const Double_t sysErrorLow = grSysErrors[i]->GetErrorYlow(ind);
        const Double_t sysErrorHigh = grSysErrors[i]->GetErrorYhigh(ind);
        
        grSysErrorsYields[i]->SetPoint(ind, grSysErrorsYields[i]->GetX()[ind], yield);
        grSysErrorsYields[i]->SetPointEYhigh(ind, totalYield * sysErrorHigh);
        grSysErrorsYields[i]->SetPointEYlow(ind, totalYield * sysErrorLow);
        
        /*
        Double_t totalYield = 0;
        for (Int_t j = 0; j < numSpecies; j++)
          totalYield += hYields[j]->GetBinContent(ind + 1);
        
        const Double_t yield = hYields[i]->GetBinContent(ind + 1);
        const Double_t fraction = hFractions[i][reference]->GetBinContent(ind + 1);
        const Double_t sysErrorLow = grSysErrors[i]->GetErrorYlow(ind);
        const Double_t sysErrorHigh = grSysErrors[i]->GetErrorYhigh(ind);
        
        if (fraction <= 0.) {
          printf("Error: Fraction = 0 for species %d. Cannot transform error....\n", i);
          return -1;
        }
        const Double_t sysErrorLowRel = sysErrorLow / fraction;
        const Double_t sysErrorHighRel = sysErrorHigh / fraction;
        
        grSysErrorsYields[i]->SetPoint(ind, grSysErrorsYields[i]->GetX()[ind], yield);
        grSysErrorsYields[i]->SetPointEYhigh(ind, sysErrorHighRel * yield);
        grSysErrorsYields[i]->SetPointEYlow(ind, sysErrorLowRel * yield);
        */
      }
    }
  }
  
  // Summary plots
  Double_t pLow = 0.15;
  Double_t pHigh = 50.;
  TCanvas* cFractionsWithSystematicError = DrawResultHistos("cFractionsWithSystematicError", "Particle fraction", pLow, pHigh, 
                                                              hFractions, reference, grSysErrors, kFraction, mode);
  
  TCanvas* cToPiRatiosWithSystematicError = DrawResultHistos("cToPiRatiosWithSystematicError", "Ratio", pLow, pHigh, 
                                                               hRatioToPi, reference, grSysErrorsToPiRatios, kToPiRatio, mode);
  
  TCanvas* cYieldWithSystematicError = 0x0;
  
  if (hasYield) {
    TH1F** hYieldsDummyArray[numSpecies];
    for (Int_t i = 0; i < numSpecies; i++) {
      hYieldsDummyArray[i] = new TH1F*[numFiles];
      hYieldsDummyArray[i][reference] = hYields[i];
    }
    cYieldWithSystematicError = DrawResultHistos("cYieldWithSystematicError", "Ratio", pLow, pHigh, 
                                                 hYieldsDummyArray, reference, grSysErrorsYields, kYield, mode);
  }
  
  // Output file
  TFile* fSave = 0x0;
  TDatime daTime;
  TString saveFileName;
  
  TString modeString[4] = { "", "_pT", "_z", "_xi" };
    
  saveFileName = Form("outputSystematics_%s%s_nSigma%.1f__%04d_%02d_%02d.root", outFileTitle.Data(), modeString[mode + 1].Data(),
                      nSigma, daTime.GetYear(),
                      daTime.GetMonth(), daTime.GetDay());
    
  fSave = TFile::Open(Form("%s/%s", path.Data(), saveFileName.Data()), "recreate");
  if (!fSave) {
    std::cout << "Failed to open save file \"" << Form("%s/%s", path.Data(), saveFileName.Data()) << "\"!" << std::endl;
    return -1;
  }
  
  // Save final results
  fSave->cd();
  
  for (Int_t i = 0; i < numFiles; i++) {
    if (hSystematicsElectrons[i])
      hSystematicsElectrons[i]->Write();
    
    if (hSystematicsPions[i])
      hSystematicsPions[i]->Write();
    
    if (hSystematicsKaons[i])
      hSystematicsKaons[i]->Write();
    
    if (hSystematicsProtons[i])
      hSystematicsProtons[i]->Write();
    
    if (hSystematicsMuons[i])
      hSystematicsMuons[i]->Write();
  }
  
  if (cSystematicsElectrons)
      cSystematicsElectrons->Write();
  
  if (cSystematicsPions)
      cSystematicsPions->Write();
  
  if (cSystematicsKaons)
      cSystematicsKaons->Write();
  
  if (cSystematicsProtons)
      cSystematicsProtons->Write();
  
  if (cSystematicsMuons)
      cSystematicsMuons->Write();
  
  if (cFractionsWithSystematicError)
      cFractionsWithSystematicError->Write();
  
  for (Int_t i = 0; i < numFiles; i++) {
    if (hSystematicsToPiRatioElectrons[i])
      hSystematicsToPiRatioElectrons[i]->Write();
    
    if (hSystematicsToPiRatioKaons[i])
      hSystematicsToPiRatioKaons[i]->Write();
    
    if (hSystematicsToPiRatioProtons[i])
      hSystematicsToPiRatioProtons[i]->Write();
    
    if (hSystematicsToPiRatioMuons[i])
      hSystematicsToPiRatioMuons[i]->Write();
  }
  
  if (cSystematicsToPiRatioElectrons)
      cSystematicsToPiRatioElectrons->Write();
  
  if (cSystematicsToPiRatioKaons)
      cSystematicsToPiRatioKaons->Write();
  
  if (cSystematicsToPiRatioProtons)
      cSystematicsToPiRatioProtons->Write();
  
  if (cSystematicsToPiRatioMuons)
      cSystematicsToPiRatioMuons->Write();
  
  if (cToPiRatiosWithSystematicError)
      cToPiRatiosWithSystematicError->Write();
  
  for (Int_t i = 0; i < numSpecies; i++) {
    if (hFractions[i][reference])
      hFractions[i][reference]->Write();
    
    if (grSysErrors[i])
      grSysErrors[i]->Write();
    
    if (hYields[i])
      hYields[i]->Write();
    
    if (grSysErrorsYields[i])
      grSysErrorsYields[i]->Write();
    
    if (hRatioToPi[i][reference])
      hRatioToPi[i][reference]->Write();
    
    if (grSysErrorsToPiRatios[i])
      grSysErrorsToPiRatios[i]->Write();
  }
  
  if (cYieldWithSystematicError)
      cYieldWithSystematicError->Write();
  
  
  if (hNumEvents)
    hNumEvents->Write();
  
  if (hNumEventsTriggerSel)
    hNumEventsTriggerSel->Write();
  
  if (hNumEventsTriggerSelVtxCut)
    hNumEventsTriggerSelVtxCut->Write();
  
  if (hNumEventsTriggerSelVtxCutNoPileUpRejection)
    hNumEventsTriggerSelVtxCutNoPileUpRejection->Write();
  
  if (hNjetsGen)
    hNjetsGen->Write();
  
  if (hNjetsRec)
    hNjetsRec->Write();
  
  // Save list of file names in output file
  TString listOfFileNames = "";
  for (Int_t i = 0; i < numFiles; i++) {
    listOfFileNames.Append(Form("%s%d: %s", i == 0 ? "" : ", ", i, fileNames[i].Data()));
  }
  
  TNamed* settings = new TNamed(Form("Used files for systematics: %s\n", listOfFileNames.Data()),
                                Form("Used files for systematics: %s\n", listOfFileNames.Data()));
  settings->Write();
    
  fSave->Close();
  
  delete cSystematicsElectrons;
  delete cSystematicsPions;
  delete cSystematicsKaons;
  delete cSystematicsMuons;
  delete cSystematicsProtons;
  delete cFractionsWithSystematicError;
  
  delete cSystematicsToPiRatioElectrons;
  delete cSystematicsToPiRatioKaons;
  delete cSystematicsToPiRatioMuons;
  delete cSystematicsToPiRatioProtons;
  delete cToPiRatiosWithSystematicError;
  
  delete cYieldWithSystematicError;
  
  for (Int_t iFile = 0; iFile < numFiles; iFile++)
    f[iFile]->Close();
  
  return 0;
}
