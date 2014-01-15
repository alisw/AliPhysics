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

const Int_t numSpecies = 5;

//________________________________________________________
TCanvas* calculateSystematics(TString canvName, TString canvTitle, TH1F** histos, Int_t numHistos, Int_t speciesID, Double_t /*nSigma*/,
                              const TString* systematicsHistosName, Int_t reference, TH1F** hSystematics, TGraphAsymmErrors** gr,
                              Bool_t ignoreSigmaErrors)
{
  // For every bin:
  // Since the method with the root finding already takes into account the statistical error,
  // there is no need to use nSigma > 0.
  // If the statistical error is ignored, nevertheless don't use nSigma > 0 because this might
  // give zero systematic error for high pT, which is usually not accepted by people, although
  // the natural point of view "no systematic visible for given statistical error" is reasonable to me.
  
  Double_t ymax = 0;
  Double_t ymin = 0;
  
  
  // Just for drawing
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
  
  Double_t meansForFit[numHistos];
  Double_t sigmasForFit[numHistos];
  
  for (Int_t bin = 0; bin < nBins; bin++) {
    x[bin] = histos[reference]->GetBinCenter(bin + 1);
    xerr[bin] = histos[reference]->GetBinWidth(bin + 1) / 2.;
    y[bin] = histos[reference]->GetBinContent(bin + 1);
    
    for (Int_t j = 0; j < numHistos; j++) {
      meansForFit[j] = histos[j]->GetBinContent(bin + 1);
      sigmasForFit[j] = histos[j]->GetBinError(bin + 1);
    }
  
    yerrl[bin] = yerrh[bin] = findSystematicError(numHistos, meansForFit, sigmasForFit, ignoreSigmaErrors);
  }
  
  TGraphAsymmErrors* gTemp = new TGraphAsymmErrors(nBins, x, y, xerr, xerr, yerrl, yerrh);
  *gr = gTemp;
  (*gr)->SetName(Form("systematicError_%s", AliPID::ParticleName(speciesID)));
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
TCanvas* DrawFractionHistos(TString canvName, TString canvTitle, Double_t pLow, Double_t pHigh, TH1F*** hist, Int_t reference,
                            TGraphAsymmErrors** gr)
{
  TCanvas* canv = new TCanvas(canvName.Data(), canvTitle.Data(),100,10,1200,800);
  canv->SetGridx(1);
  canv->SetGridy(1);
  canv->SetLogx(1);
  for (Int_t i = 0; i < numSpecies; i++) {
    hist[i][reference]->GetYaxis()->SetRangeUser(0.0, 1.0);
    hist[i][reference]->GetYaxis()->SetTitle(canvTitle.Data());
    hist[i][reference]->GetXaxis()->SetRangeUser(pLow, pHigh);
    //hist[i][reference]->SetFillStyle(3004 + i);
    //hist[i][reference]->SetFillColor(kGray);
    hist[i][reference]->SetFillStyle(0);
    hist[i][reference]->SetFillColor(hist[i][reference]->GetMarkerColor());
    hist[i][reference]->SetLineColor(hist[i][reference]->GetMarkerColor());
  }
  hist[2][reference]->SetMarkerStyle(20);
  hist[2][reference]->Draw("e p");
  hist[0][reference]->SetMarkerStyle(21);
  hist[0][reference]->Draw("e p same");
  hist[1][reference]->SetMarkerStyle(22);
  hist[1][reference]->Draw("e p same");
  hist[3][reference]->SetMarkerStyle(29);
  hist[3][reference]->Draw("e p same");
  hist[4][reference]->SetMarkerStyle(30);
  hist[4][reference]->Draw("e p same");
  
  gr[0]->Draw("2 same");
  gr[1]->Draw("2 same");
  gr[2]->Draw("2 same");
  gr[3]->Draw("2 same");
  gr[4]->Draw("2 same");
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(hist[2][reference], "#pi", "flp");
  legend->AddEntry(hist[0][reference], "e", "flp");
  legend->AddEntry(hist[1][reference], "K", "flp");
  legend->AddEntry(hist[3][reference], "p", "flp");
  legend->AddEntry(hist[4][reference], "#mu", "flp");
  legend->Draw();
  
  return canv;
}


//________________________________________________________
TH1F* loadHisto(const TString histName, TFile* f)
{
  if (!f) {
    std::cout << "No file. Cannot load hist \"" << histName.Data() << "\n!" << std::endl;
    return 0x0;
  }
  
  TH1F* hTemp = dynamic_cast<TH1F*>(f->Get(histName.Data()));
  if (!hTemp) {
    std::cout << "Failed to load histo \"" << histName.Data() << "\"!" << std::endl;
    return 0x0;
  } 
  
  return hTemp;
}


//________________________________________________________
Int_t SystematicErrorEstimation(const TString path, const TString outFileTitle, const TString* fileNames, const TString* histTitles,
                                const Int_t numFiles, const Double_t nSigma, const Bool_t ignoreSigmaErrors) 
{ 
  if (!fileNames || numFiles < 1)
    return -1;
  
  TFile* f[numFiles];
  TH1F** hFractions[numSpecies];
  for (Int_t i = 0; i < numSpecies; i++) 
    hFractions[i] = new TH1F*[numFiles];
  
  const Int_t reference = 0;
  TH1F* hYields[numSpecies]; // Only the reference yields
  
  const TString histNames[numSpecies] = {"hFractionElectrons", "hFractionKaons", "hFractionPions", "hFractionProtons", "hFractionMuons" };
  
  const TString histNamesYields[numSpecies] = {"hYieldElectrons", "hYieldKaons", "hYieldPions", "hYieldProtons", "hYieldMuons" };
    
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
      
      if (iFile == reference) {
        hYields[i] = loadHisto(histNamesYields[i], f[iFile]);
        if (!hYields[i])
          return -1;
      }
    }
  }
    
  
  TGraphAsymmErrors* grSysErrors[numSpecies] = {0x0,};
  TGraphAsymmErrors* grSysErrorsYields[numSpecies] = {0x0,};
  
  TH1F* hSystematicsPions[numFiles];
  TCanvas* cSystematicsPions = calculateSystematics("cSystematicsPions", "Systematics Pions", hFractions[2], numFiles,
                                                    AliPID::kPion, nSigma,
                                                    histTitles, reference, hSystematicsPions, &grSysErrors[2], ignoreSigmaErrors);

  TH1F* hSystematicsElectrons[numFiles];
  TCanvas* cSystematicsElectrons = calculateSystematics("cSystematicsElectrons", "Systematics Electrons", hFractions[0], numFiles,
                                                        AliPID::kElectron,
                                                        nSigma, histTitles, reference, hSystematicsElectrons,
                                                        &grSysErrors[0], ignoreSigmaErrors);
  
  TH1F* hSystematicsKaons[numFiles];
  TCanvas* cSystematicsKaons = calculateSystematics("cSystematicsKaons", "Systematics Kaons", hFractions[1], numFiles, AliPID::kKaon, nSigma,
                                                    histTitles, reference, hSystematicsKaons, &grSysErrors[1], ignoreSigmaErrors);
  
  TH1F* hSystematicsProtons[numFiles];  
  TCanvas* cSystematicsProtons = calculateSystematics("cSystematicsProtons", "Systematics Protons", hFractions[3], numFiles,
                                                      AliPID::kProton, nSigma,
                                                      histTitles, reference, hSystematicsProtons, &grSysErrors[3], ignoreSigmaErrors);
  
  TH1F* hSystematicsMuons[numFiles];
  TCanvas* cSystematicsMuons = calculateSystematics("cSystematicsMuons", "Systematics Muons", hFractions[4], numFiles,
                                                    AliPID::kMuon, nSigma,
                                                    histTitles, reference, hSystematicsMuons, &grSysErrors[4], ignoreSigmaErrors);
  
  Double_t pLow = 0.15;
  Double_t pHigh = 50.;
  TCanvas* cFractionsWithSystematicError = DrawFractionHistos("cFractionsWithSystematicError", "Particle fractions", pLow, pHigh, hFractions, reference, 
                                                              grSysErrors);
  
  
  //TODO At the moment, the error of the fractions and the yield is just a constant factor (number of tracks in that bin)
  // (-> But this can change in future (I have to think about it a little bit more carefully)).
  // Thus, the relative errors are the same for fractions and yields and I can just use this fact to
  // transform the errors from the fractions to those of the yields.
  // However, this causes trouble in case of fraction = 0. Therefore, sum up the yields to the total yield and use this for scaling
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
  
  
  // Output file
  TFile* fSave = 0x0;
  TDatime daTime;
  TString saveFileName;
    
  saveFileName = Form("outputSystematics_%s_nSigma%.1f__%04d_%02d_%02d.root", outFileTitle.Data(), nSigma, daTime.GetYear(),
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
  
  for (Int_t i = 0; i < numSpecies; i++) {
    if (hFractions[i][reference])
      hFractions[i][reference]->Write();
    
    if (grSysErrors[i])
      grSysErrors[i]->Write();
    
    if (hYields[i])
      hYields[i]->Write();
    
    if (grSysErrorsYields[i])
      grSysErrorsYields[i]->Write();
  }
  
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
  
  return 0;
}
