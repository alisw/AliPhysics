#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"

#include "AliPID.h"

#include <iostream>

#include "THnSparseDefinitions.h"

void normaliseHist(TH2* h)
{
  if (!h)
    return;
  
  for (Int_t binX = 1; binX <= h->GetNbinsX(); binX++) {
    for (Int_t binY = 1; binY <= h->GetNbinsY(); binY++) {
      Double_t widthFactor = h->GetXaxis()->GetBinWidth(binX) * h->GetYaxis()->GetBinWidth(binY);
      if (widthFactor > 0) {
        h->SetBinContent(binX, binY, h->GetBinContent(binX, binY) / widthFactor);
        h->SetBinError(binX, binY, h->GetBinError(binX, binY) / widthFactor);
      }
    }
  }
}


//___________________________________________________________________
Int_t extractPtResolution(TString pathNameData, TString listName,
                          Int_t chargeMode /*kNegCharge = -1, kAllCharged = 0, kPosCharge = 1*/,
                          Double_t lowerCentrality /*= -2*/, Double_t upperCentrality /*= -2*/,
                          Double_t lowerJetPt /*= -1*/ , Double_t upperJetPt/* = -1*/)
{
  if (listName == "") {
    listName = pathNameData;
    listName.Replace(0, listName.Last('/') + 1, "");
    listName.ReplaceAll(".root", "");
  }
  
  TString pathData = pathNameData;
  pathData.Replace(pathData.Last('/'), pathData.Length(), "");
  
  TH1D* hPtResolutionFit[AliPID::kSPECIES] = {0x0, };
  TH2D* hPtResolution[AliPID::kSPECIES] = {0x0, };
  THnSparse* hPtResolutionRaw[AliPID::kSPECIES] = {0x0, };

  TFile* fileData = TFile::Open(pathNameData.Data());
  if (!fileData) {
    printf("Failed to open data file \"%s\"\n", pathNameData.Data());
    return -1;
  }
  
  TObjArray* histList = (TObjArray*)(fileData->Get(listName.Data()));
  
  if (!histList) {
    printf("Failed to load list!\n");
    return -1;
  }
  
  Double_t actualLowerCentrality = -2;
  Double_t actualUpperCentrality = -2;
  
  Double_t actualLowerJetPt = -1.;
  Double_t actualUpperJetPt = -1.;
  
  Bool_t restrictJetPtAxis = (lowerJetPt >= 0 && upperJetPt >= 0);
  const Bool_t restrictCentrality = ((lowerCentrality >= -1) && (upperCentrality >= -1));
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    const TString sparseName = Form("fPtResolution_%s", AliPID::ParticleShortName(species));
    hPtResolutionRaw[species] = (THnSparse*)histList->FindObject(sparseName.Data());
    
    if (!hPtResolutionRaw[species]) {
      printf("Failed to load THnSparse for %s: %s!\n", AliPID::ParticleShortName(species), sparseName.Data());
      return -1;
    }
    
    // Set proper errors, if not yet calculated
    if (!hPtResolutionRaw[species]->GetCalculateErrors()) {
      std::cout << "Re-calculating errors of " << hPtResolutionRaw[species]->GetName() << "..." << std::endl;
      
      hPtResolutionRaw[species]->Sumw2();
      
      Long64_t nBinsPtResolutionRaw = hPtResolutionRaw[species]->GetNbins();
      Double_t binContent = 0;
      for (Long64_t bin = 0; bin < nBinsPtResolutionRaw; bin++) {
        binContent = hPtResolutionRaw[species]->GetBinContent(bin);
        hPtResolutionRaw[species]->SetBinError(bin, TMath::Sqrt(binContent));
      }
    }
    
    // Integral(lowerCentBinLimit, uppCentBinLimit) will not be restricted if these values are kept
    const Int_t lowerCentralityBinLimit = restrictCentrality ? hPtResolutionRaw[species]->GetAxis(kPtResCentrality)->FindBin(lowerCentrality + 0.001) 
                                                                    : -1;
    const Int_t upperCentralityBinLimit = restrictCentrality ? hPtResolutionRaw[species]->GetAxis(kPtResCentrality)->FindBin(upperCentrality - 0.001) 
                                                                    : -2;
    
    
    if (restrictCentrality) {
      actualLowerCentrality = hPtResolutionRaw[species]->GetAxis(kPtResCentrality)->GetBinLowEdge(lowerCentralityBinLimit);
      actualUpperCentrality = hPtResolutionRaw[species]->GetAxis(kPtResCentrality)->GetBinLowEdge(upperCentralityBinLimit);
      
      hPtResolutionRaw[species]->GetAxis(kPtResCentrality)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
    }
    
    const Bool_t restrictCharge = (chargeMode != kAllCharged);
    
    Int_t lowerChargeBinLimit = -1;
    Int_t upperChargeBinLimit = -2;
      
    if (restrictCharge) {
      // Add subtract a very small number to avoid problems with values right on the border between to bins
      if (chargeMode == kNegCharge) {
        lowerChargeBinLimit = hPtResolutionRaw[species]->GetAxis(kPtResCharge)->FindBin(-1. + 0.001);
        upperChargeBinLimit = hPtResolutionRaw[species]->GetAxis(kPtResCharge)->FindBin(0. - 0.001);
      }
      else if (chargeMode == kPosCharge) {
        lowerChargeBinLimit = hPtResolutionRaw[species]->GetAxis(kPtResCharge)->FindBin(0. + 0.001);
        upperChargeBinLimit = hPtResolutionRaw[species]->GetAxis(kPtResCharge)->FindBin(1. - 0.001);
      }
      
      // Check if the values look reasonable
      if (lowerChargeBinLimit <= upperChargeBinLimit && lowerChargeBinLimit >= 1
          && upperChargeBinLimit <= hPtResolutionRaw[species]->GetAxis(kPtResCharge)->GetNbins()) {
        // OK
      }
      else {
        std::cout << std::endl;
        std::cout << "Requested charge range out of limits or upper and lower limit are switched!" << std::endl;
        return -1;
      }
      
      hPtResolutionRaw[species]->GetAxis(kPtResCharge)->SetRange(lowerChargeBinLimit, upperChargeBinLimit);
    }
    
    // If desired, restrict jetPt axis
    Int_t lowerJetPtBinLimit = -1;
    Int_t upperJetPtBinLimit = -1;
    
    if (restrictJetPtAxis) {
      // Add subtract a very small number to avoid problems with values right on the border between to bins
      lowerJetPtBinLimit = hPtResolutionRaw[species]->GetAxis(kPtResJetPt)->FindBin(lowerJetPt + 0.001);
      upperJetPtBinLimit = hPtResolutionRaw[species]->GetAxis(kPtResJetPt)->FindBin(upperJetPt - 0.001);
      
      // Check if the values look reasonable
      if (lowerJetPtBinLimit <= upperJetPtBinLimit && lowerJetPtBinLimit >= 1 &&
          upperJetPtBinLimit <= hPtResolutionRaw[species]->GetAxis(kPtResJetPt)->GetNbins()) {
        actualLowerJetPt = hPtResolutionRaw[species]->GetAxis(kPtResJetPt)->GetBinLowEdge(lowerJetPtBinLimit);
        actualUpperJetPt = hPtResolutionRaw[species]->GetAxis(kPtResJetPt)->GetBinUpEdge(upperJetPtBinLimit);

        restrictJetPtAxis = kTRUE;
      }
      else {
        std::cout << std::endl;
        std::cout << "Requested jet pT range out of limits or upper and lower limit are switched!" << std::endl;
        return -1;
      }
    }
    
    std::cout << "jet pT: ";
    if (restrictJetPtAxis) {
      std::cout << actualLowerJetPt << " - " << actualUpperJetPt << std::endl;
      hPtResolutionRaw[species]->GetAxis(kPtResJetPt)->SetRange(lowerJetPtBinLimit, upperJetPtBinLimit);
    }
    else {
      std::cout << "All" << std::endl;
    }
    
    hPtResolution[species] = hPtResolutionRaw[species]->Projection(kPtResGenPt, kPtResRecPt, "e");
    hPtResolution[species]->SetName(Form("hPtResolution_%s", AliPID::ParticleShortName(species)));
    hPtResolution[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
    hPtResolution[species]->SetStats(kFALSE);
    hPtResolution[species]->SetLineColor(getLineColorAliPID(species));
    hPtResolution[species]->SetMarkerColor(getLineColorAliPID(species));
    
    normaliseHist(hPtResolution[species]);
    
    TObjArray aSlices;
    hPtResolution[species]->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices);
    TH1D* hMean = (TH1D*)(aSlices.At(1));
    TH1D* hSigma = (TH1D*)(aSlices.At(2));
    
    hPtResolutionFit[species] = new TH1D(*hSigma);
    hPtResolutionFit[species]->SetName(Form("hPtResolutionFit_%s", AliPID::ParticleShortName(species)));
    hPtResolutionFit[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
    hPtResolutionFit[species]->SetLineColor(getLineColorAliPID(species));
    hPtResolutionFit[species]->SetMarkerColor(getLineColorAliPID(species));
    
    hPtResolutionFit[species]->Divide(hSigma, hMean);
  }
  
  
  
  // Save results to file
  TString chargeString = "";
  if (chargeMode == kPosCharge)
    chargeString = "_posCharge";
  else if (chargeMode == kNegCharge)
    chargeString = "_negCharge";
  
  TString saveFileName = pathNameData;
  saveFileName.Replace(0, pathNameData.Last('/') + 1, "");
  
  TString savePath = pathNameData;
  savePath.ReplaceAll(Form("/%s", saveFileName.Data()), "");
  
  saveFileName.Prepend("output_extractedPTResolution_");
  TString centralityString = restrictCentrality ? Form("_centrality_%.0f_%.0f.root", actualLowerCentrality,
                                                       actualUpperCentrality)
                                                    : "_centrality_all";
  TString jetPtString = restrictJetPtAxis ? Form("_jetPt_%.0f_%.0f.root", actualLowerJetPt, actualUpperJetPt)
                                          : "";  
  saveFileName.ReplaceAll(".root", Form("%s%s%s.root", centralityString.Data(), jetPtString.Data(), chargeString.Data()));
  
  TString saveFilePathName = Form("%s/%s", savePath.Data(), saveFileName.Data());
  TFile* saveFile = TFile::Open(saveFilePathName.Data(), "RECREATE");
  
  if (!saveFile) {
    printf("Failed to save results to file \"%s\"!\n", saveFilePathName.Data());
    return -1;
  }
  
  saveFile->cd();
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (hPtResolution[species])
      hPtResolution[species]->Write();
    
    if (hPtResolutionFit[species])
      hPtResolutionFit[species]->Write();
  }
  
  TNamed* settings = new TNamed(
      Form("Settings: Data file \"%s\", lowerCentrality %.3f, upperCentrality %.3f, lowerJetPt %.1f, upperJetPt %.1f\n",
           pathNameData.Data(), lowerCentrality, upperCentrality, lowerJetPt, upperJetPt), "");
  settings->Write();
  
  saveFile->Close();
  
  
  return 0;
}