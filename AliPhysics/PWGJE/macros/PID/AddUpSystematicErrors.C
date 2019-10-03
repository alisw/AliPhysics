#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TMath.h"

#include "AliPID.h"

#include "THnSparseDefinitions.h"

#include <iostream>
#include <iomanip>

const Int_t numSpecies = AliPID::kSPECIES;

enum resultType { kFraction = 0, kYield = 1, kToPiRatio = 2 };

//________________________________________________________
TCanvas* DrawFractionHistos(TString canvName, TString canvTitle, Double_t pLow, Double_t pHigh, TGraphAsymmErrors** gr,  TH1F** histRef, 
                            Int_t resultType, Int_t mode)
{
  if (!histRef || !gr)
    return 0x0;
  
  TCanvas* canv = new TCanvas(canvName.Data(), canvTitle.Data(),100,10,1200,800);
  canv->SetGridx(1);
  canv->SetGridy(1);
  if (mode <= 0)
    canv->SetLogx(1);
  if (resultType == kYield)
    canv->SetLogy(1);
  
  for (Int_t i = 0; i < numSpecies; i++) {
    if (!histRef[i])
      continue;
    if (resultType == kYield) {
      if (mode < 0)
        histRef[i]->GetYaxis()->SetRangeUser(2e-12, 19);
      else  
        histRef[i]->GetYaxis()->SetRangeUser(2e-7, 9);
    }
    else
      histRef[i]->GetYaxis()->SetRangeUser(0.0, 1.0);
    //histRef[i]->GetYaxis()->SetTitle(canvTitle.Data());
    if (mode <= 0)
      histRef[i]->GetXaxis()->SetRangeUser(pLow, pHigh);
    //histRef[i]->SetFillStyle(3004 + i);
    //histRef[i]->SetFillColor(kGray);
    histRef[i]->SetFillStyle(0);
    histRef[i]->SetFillColor(histRef[i]->GetMarkerColor());
    histRef[i]->SetLineColor(histRef[i]->GetMarkerColor());
  }
  if (resultType != kToPiRatio) {
    histRef[AliPID::kPion]->SetMarkerStyle(20);
    histRef[AliPID::kPion]->Draw("e p");
  }
  histRef[AliPID::kElectron]->SetMarkerStyle(21);
  histRef[AliPID::kElectron]->Draw(Form("e p%s", (resultType == kToPiRatio) ? "" : " same"));
  histRef[AliPID::kKaon]->SetMarkerStyle(22);
  histRef[AliPID::kKaon]->Draw("e p same");
  histRef[AliPID::kProton]->SetMarkerStyle(29);
  histRef[AliPID::kProton]->Draw("e p same");
  histRef[AliPID::kMuon]->SetMarkerStyle(30);
  histRef[AliPID::kMuon]->Draw("e p same");
  
  if (mode <= 0)
    gr[AliPID::kElectron]->GetHistogram()->GetXaxis()->SetRangeUser(pLow, pHigh);
  gr[AliPID::kElectron]->GetHistogram()->GetYaxis()->SetRangeUser(0., 1.0);
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
    legend->AddEntry(histRef[AliPID::kPion], Form("#pi%s", suffix.Data()), "flp");
  legend->AddEntry(histRef[AliPID::kElectron], Form("e%s", suffix.Data()), "flp");
  legend->AddEntry(histRef[AliPID::kKaon], Form("K%s", suffix.Data()), "flp");
  legend->AddEntry(histRef[AliPID::kProton], Form("p%s", suffix.Data()), "flp");
  legend->AddEntry(histRef[AliPID::kMuon], Form("#mu%s", suffix.Data()), "flp");
  legend->Draw();
  
  ClearTitleFromHistoInCanvas(canv);
  
  return canv;
}


//________________________________________________________
TGraphAsymmErrors* loadGraph(const TString graphName, TFile* f)
{
  if (!f) {
    std::cout << "No file. Cannot load graph \"" << graphName.Data() << "\n!" << std::endl;
    return 0x0;
  }
  
  TGraphAsymmErrors* grTemp = dynamic_cast<TGraphAsymmErrors*>(f->Get(graphName.Data()));
  if (!grTemp) {
    std::cout << "Failed to load histo \"" << graphName.Data() << "\"!" << std::endl;
    return 0x0;
  } 
  
  return grTemp;
}


//________________________________________________________
TH1F* loadHisto(const TString histName, TFile* f)
{
  if (!f) {
    std::cout << "No file. Cannot load hist \"" << histName.Data() << "\n!" << std::endl;
    return 0x0;
  }
  
  TH1F* hTemp = (TH1F*)(f->Get(histName.Data()));
  if (!hTemp) {
    std::cout << "Failed to load histo \"" << histName.Data() << "\"!" << std::endl;
    return 0x0;
  } 
  
  return hTemp;
}


//________________________________________________________
Int_t AddUpSystematicErrors(const TString path, const TString outFileTitle, const TString* fileNames, const Int_t numFiles,
                            const TString fileNameReference, const Int_t mode = -1, const Bool_t useDirectFit = kTRUE) 
{ 
  if (!fileNames || numFiles < 1)
    return -1;
  
  // NOTE: the names of the systematic errors are always the same, only the names of the reference histos change!
  // Histo names for inclusive case (mode = -1)
  TFile* f[numFiles];
  TGraphAsymmErrors** grSysError[numSpecies] = { 0x0, };
  TString grNames[numSpecies] = {"systematicError_electron", "systematicError_muon", "systematicError_pion", "systematicError_kaon", 
                                 "systematicError_proton" };
  
  TString histNames[numSpecies] = {"hFractionElectrons", "hFractionMuons", "hFractionPions", "hFractionKaons", "hFractionProtons" };
  
  TGraphAsymmErrors** grSysErrorToPiRatios[numSpecies] = { 0x0, };
  TString grNamesToPiRatios[numSpecies] = {"systematicErrorToPiRatio_electron", "systematicErrorToPiRatio_muon", "",          
                                           "systematicErrorToPiRatio_kaon", "systematicErrorToPiRatio_proton",  };
  
  TString histNamesToPiRatios[numSpecies] = {"hRatioToPiElectrons", "hRatioToPiMuons", "", "hRatioToPiKaons", "hRatioToPiProtons" };
  
  TGraphAsymmErrors** grSysErrorYields[numSpecies] = { 0x0, };
  TString grNamesYields[numSpecies] = {"systematicErrorYields_electron", "systematicErrorYields_muon", "systematicErrorYields_pion",
                                       "systematicErrorYields_kaon", "systematicErrorYields_proton"  };
                                 
  TString histNamesYields[numSpecies] = {"hYieldElectrons", "hYieldMuons", "hYieldPions", "hYieldKaons", "hYieldProtons" };
  
  
  // Hist names for jet case (mode = 0,1,2 = pT,z,xi) && !useDirectFit
  if (!useDirectFit) {
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
  }
  
  for (Int_t i = 0; i < numSpecies; i++) {
    grSysError[i] = new TGraphAsymmErrors*[numFiles];
    grSysErrorToPiRatios[i] = new TGraphAsymmErrors*[numFiles];
    grSysErrorYields[i] = new TGraphAsymmErrors*[numFiles];
    
    for (Int_t j = 0; j < numFiles; j++) {
      grSysError[i][j] = 0x0;
      grSysErrorToPiRatios[i][j] = 0x0;
      grSysErrorYields[i][j] = 0x0;
    }
  }
  
  Bool_t hasToPiRatios = kFALSE;
  Bool_t hasYield = kFALSE;
    
  for (Int_t iFile = 0; iFile < numFiles; iFile++) {
    f[iFile] = TFile::Open(fileNames[iFile].Data());
    if (!f[iFile])  {
      std::cout << "Failed to open file \"" << fileNames[iFile].Data() << "\"!" << std::endl;
      return -1;
    }
        
    // Extract the data graphs
    for (Int_t i = 0; i < numSpecies; i++) {
      grSysError[i][iFile] = loadGraph(grNames[i], f[iFile]);
      if (!grSysError[i][iFile])
        return -1;
      
      grSysError[i][iFile]->SetName(Form("%s_fileID%d", grSysError[i][iFile]->GetName(), iFile));
      
      if (grNamesToPiRatios[i] != "") {
        grSysErrorToPiRatios[i][iFile] = loadGraph(grNamesToPiRatios[i], f[iFile]);
        if (!grSysErrorToPiRatios[i][iFile]) {
          if (hasToPiRatios)
            return -1;
        }
        else
          hasToPiRatios = kTRUE;
        
        grSysErrorToPiRatios[i][iFile]->SetName(Form("%s_fileID%d", grSysErrorToPiRatios[i][iFile]->GetName(), iFile));
      }
      
      grSysErrorYields[i][iFile] = loadGraph(grNamesYields[i], f[iFile]);
      if (!grSysErrorYields[i][iFile])
        continue;
      else
        hasYield = kTRUE;
      
      grSysErrorYields[i][iFile]->SetName(Form("%s_fileID%d", grSysErrorYields[i][iFile]->GetName(), iFile));
    }
  }
  
  // Load data points with statistical errors
  TFile* fReferenceData = TFile::Open(fileNameReference.Data());
  if (!fReferenceData) {
    std::cout << "Failed to open file \"" << fileNameReference.Data() << "\"!" << std::endl;
    return -1;
  }
    
  TH1F* hReferenceFractions[numSpecies] = { 0x0, };
  TH1F* hReferenceToPiRatios[numSpecies] = { 0x0, };
  TH1F* hReferenceYields[numSpecies] = { 0x0, };
  
  for (Int_t i = 0; i < numSpecies; i++) {
    hReferenceFractions[i] = loadHisto(histNames[i], fReferenceData);
    if (!hReferenceFractions[i])
      return -1;
    
    hReferenceYields[i] = loadHisto(histNamesYields[i], fReferenceData);
    if (!hReferenceYields[i] && hasYield)
      return -1;
    
    if (histNamesToPiRatios[i] != "") {
      hReferenceToPiRatios[i] = loadHisto(histNamesToPiRatios[i], fReferenceData);
      if (!hReferenceToPiRatios[i] && hasToPiRatios)
        return -1;
    }
  }
  
  // Forward number of events for further processing
  TH1* hNumEvents = loadHisto("fhEventsProcessed", fReferenceData);
  TH1* hNumEventsTriggerSel = loadHisto("fhEventsTriggerSel", fReferenceData);
  TH1* hNumEventsTriggerSelVtxCut = loadHisto("fhEventsTriggerSelVtxCut", fReferenceData);
  TH1* hNumEventsTriggerSelVtxCutNoPileUpRejection = loadHisto("fhEventsProcessedNoPileUpRejection", fReferenceData);
  
  // Same for number of jet histos
  TH1* hNjetsGen = loadHisto("fh2FFJetPtGen", fReferenceData);
  TH1* hNjetsRec = loadHisto("fh2FFJetPtRec", fReferenceData);
  
  const Int_t reference = 0;
  
  // The x,y coordinates should be those of the reference graph. Then, all corresponding sys. errors are added in quadrature.
  TGraphAsymmErrors* grTotSysError[numSpecies] = { 0x0, };
  TGraphAsymmErrors* grTotSysErrorToPiRatios[numSpecies] = { 0x0, };
  TGraphAsymmErrors* grTotSysErrorYields[numSpecies] = { 0x0, };
  Double_t sysErrorTotalSquared = 0;
  Double_t temp = 0;
  
  for (Int_t i = 0; i < numSpecies; i++) {
    grTotSysError[i] = new TGraphAsymmErrors(*grSysError[i][reference]);
    grTotSysError[i]->SetName(grNames[i]);
    
    for (Int_t iPoint = 0; iPoint < grTotSysError[i]->GetN(); iPoint++) {
      sysErrorTotalSquared = 0;
      for (Int_t iFile = 0; iFile < numFiles; iFile++) {
        // Already averages high and low value -> Since they are by now the same, this is ok.
        temp = grSysError[i][iFile]->GetErrorY(iPoint); 
        if (temp > 0) {
          sysErrorTotalSquared += temp * temp;
        }
      }
      
      grTotSysError[i]->SetPointEYlow(iPoint, TMath::Sqrt(sysErrorTotalSquared));
      grTotSysError[i]->SetPointEYhigh(iPoint, TMath::Sqrt(sysErrorTotalSquared));
    }
    
    // Same for the yields
    if (hasYield) {
      grTotSysErrorYields[i] = new TGraphAsymmErrors(*grSysErrorYields[i][reference]);
      grTotSysErrorYields[i]->SetName(grNamesYields[i]);
      
      for (Int_t iPoint = 0; iPoint < grTotSysErrorYields[i]->GetN(); iPoint++) {
        sysErrorTotalSquared = 0;
        for (Int_t iFile = 0; iFile < numFiles; iFile++) {
          // Already averages high and low value -> Since they are by now the same, this is ok.
          temp = grSysErrorYields[i][iFile]->GetErrorY(iPoint); 
          if (temp > 0) {
            sysErrorTotalSquared += temp * temp;
          }
        }
        
        grTotSysErrorYields[i]->SetPointEYlow(iPoint, TMath::Sqrt(sysErrorTotalSquared));
        grTotSysErrorYields[i]->SetPointEYhigh(iPoint, TMath::Sqrt(sysErrorTotalSquared));
      }
    }
    
    // Same for the to-pi-ratios
    if (hasToPiRatios && grSysErrorToPiRatios[i][reference]) {
      grTotSysErrorToPiRatios[i] = new TGraphAsymmErrors(*grSysErrorToPiRatios[i][reference]);
      grTotSysErrorToPiRatios[i]->SetName(grNamesToPiRatios[i]);
      
      for (Int_t iPoint = 0; iPoint < grTotSysErrorToPiRatios[i]->GetN(); iPoint++) {
        sysErrorTotalSquared = 0;
        for (Int_t iFile = 0; iFile < numFiles; iFile++) {
          // Already averages high and low value -> Since they are by now the same, this is ok.
          temp = grSysErrorToPiRatios[i][iFile]->GetErrorY(iPoint); 
          if (temp > 0) {
            sysErrorTotalSquared += temp * temp;
          }
        }
        
        grTotSysErrorToPiRatios[i]->SetPointEYlow(iPoint, TMath::Sqrt(sysErrorTotalSquared));
        grTotSysErrorToPiRatios[i]->SetPointEYhigh(iPoint, TMath::Sqrt(sysErrorTotalSquared));
      }
    }
  }
  
  const Double_t pLow = 0.15;
  const Double_t pHigh = 50.;
  TCanvas* cFractionsWithTotalSystematicError = DrawFractionHistos("cFractionsWithTotalSystematicError", "Particle fraction", pLow,
                                                                   pHigh, grTotSysError, hReferenceFractions, kFraction, mode);
  
  TCanvas* cToPiRatiosWithTotalSystematicError = DrawFractionHistos("cToPiRatiosWithTotalSystematicError", "Ratio", pLow,
                                                                    pHigh, grTotSysErrorToPiRatios, hReferenceToPiRatios, kToPiRatio, 
                                                                    mode);
  
  TCanvas* cYieldsWithTotalSystematicError = 0x0;
  if (hasYield) 
    cYieldsWithTotalSystematicError = DrawFractionHistos("cYieldsWithTotalSystematicError", "Yield", pLow,
                                                         pHigh, grTotSysErrorYields, hReferenceYields, kYield, mode);
  
  // Output file
  TFile* fSave = 0x0;
  TDatime daTime;
  TString saveFileName;
  
  saveFileName = Form("outputSystematicsTotal_%s__%04d_%02d_%02d.root", outFileTitle.Data(),
                      daTime.GetYear(), daTime.GetMonth(), daTime.GetDay());
    
  fSave = TFile::Open(Form("%s/%s", path.Data(), saveFileName.Data()), "recreate");
  if (!fSave) {
    std::cout << "Failed to open save file \"" << Form("%s/%s", path.Data(), saveFileName.Data()) << "\"!" << std::endl;
    return -1;
  }
  
  // Save final results
  fSave->cd();
  
  if (cFractionsWithTotalSystematicError)
    cFractionsWithTotalSystematicError->Write();
  
  if (cToPiRatiosWithTotalSystematicError)
    cToPiRatiosWithTotalSystematicError->Write();
  
  if (cYieldsWithTotalSystematicError)
    cYieldsWithTotalSystematicError->Write();
  
  for (Int_t i = 0; i < numSpecies; i++) {
    if (grTotSysError[i])
      grTotSysError[i]->Write();
    if (hReferenceFractions[i])
      hReferenceFractions[i]->Write();
    
    if (grTotSysErrorYields[i])
      grTotSysErrorYields[i]->Write();
    if (hReferenceYields[i])
      hReferenceYields[i]->Write();
    
    if (grTotSysErrorToPiRatios[i])
      grTotSysErrorToPiRatios[i]->Write();
    if (hReferenceToPiRatios[i])
      hReferenceToPiRatios[i]->Write();
  }
  
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
  
  delete cFractionsWithTotalSystematicError;
  
  delete cToPiRatiosWithTotalSystematicError;
  
  delete cYieldsWithTotalSystematicError;

  fSave->Close();
     
  return 0;
}
