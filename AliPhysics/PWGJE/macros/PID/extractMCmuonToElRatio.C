#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "THnSparse.h"
#include "TObjArray.h"

#include <iostream>

#include "THnSparseDefinitions.h"

Int_t extractMCmuonToElRatio(const TString pathNameData, const TString listName, const Double_t lowerJetPt, const Double_t upperJetPt,
                             const Double_t lowerCentrality, const Double_t upperCentrality, const Bool_t jetHandling)
{
  TFile* fileData = TFile::Open(pathNameData.Data());
  if (!fileData) {
    printf("Failed to open data file \"%s\"\n", pathNameData.Data());
    return -1;
  }
  
  TObjArray* histList = (TObjArray*)(fileData->Get(listName.Data()));
  if (!histList) {
    std::cout << "Failed to load list \"" << listName.Data() << "\"!" << std::endl;
    return -1;
  }
  
  // Extract the data histogram
  THnSparse* hPIDdata = dynamic_cast<THnSparse*>(histList->FindObject("hPIDdataAll"));
  if (!hPIDdata) {
    std::cout << "Failed to load data histo!" << std::endl;
    return -1;
  }
  
  // Set proper errors, if not yet calculated
  if (!hPIDdata->GetCalculateErrors()) {
    std::cout << "Re-calculating errors of " << hPIDdata->GetName() << "..." << std::endl;
    hPIDdata->Sumw2();
    Long64_t nBinsTHnSparse = hPIDdata->GetNbins();
    Double_t binContent = 0;
    
    for (Long64_t bin = 0; bin < nBinsTHnSparse; bin++) {
      binContent = hPIDdata->GetBinContent(bin);
      hPIDdata->SetBinError(bin, TMath::Sqrt(binContent));
    }
  }
  
  
  Int_t lowerJetPtBinLimit = -1;
  Int_t upperJetPtBinLimit = -1;
  Bool_t restrictJetPtAxis = kFALSE;
  Double_t actualLowerJetPt = -1.;
  Double_t actualUpperJetPt = -1.;
  
  if (lowerJetPt >= 0 && upperJetPt >= 0) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerJetPtBinLimit = hPIDdata->GetAxis(kPidJetPt)->FindBin(lowerJetPt + 0.001);
    upperJetPtBinLimit = hPIDdata->GetAxis(kPidJetPt)->FindBin(upperJetPt - 0.001);
    
    // Check if the values look reasonable
    if (lowerJetPtBinLimit <= upperJetPtBinLimit && lowerJetPtBinLimit >= 1 &&
        upperJetPtBinLimit <= hPIDdata->GetAxis(kPidJetPt)->GetNbins()) {
      actualLowerJetPt = hPIDdata->GetAxis(kPidJetPt)->GetBinLowEdge(lowerJetPtBinLimit);
      actualUpperJetPt = hPIDdata->GetAxis(kPidJetPt)->GetBinUpEdge(upperJetPtBinLimit);

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
    hPIDdata->GetAxis(kPidJetPt)->SetRangeUser(actualLowerJetPt, actualUpperJetPt); 
  }
  else {
    std::cout << "All" << std::endl;
  }
  
  
  // If desired, restrict centrality axis
  Int_t lowerCentralityBinLimit = -1;
  Int_t upperCentralityBinLimit = -1;
  Bool_t restrictCentralityAxis = kFALSE;
  Double_t actualLowerCentrality = -1.;
  Double_t actualUpperCentrality = -1.;
  
  if (lowerCentrality >= -1 && upperCentrality >= -1) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerCentralityBinLimit = hPIDdata->GetAxis(kPidCentrality)->FindBin(lowerCentrality + 0.001);
    upperCentralityBinLimit = hPIDdata->GetAxis(kPidCentrality)->FindBin(upperCentrality - 0.001);
    
    // Check if the values look reasonable
    if (lowerCentralityBinLimit <= upperCentralityBinLimit && lowerCentralityBinLimit >= 1
        && upperCentralityBinLimit <= hPIDdata->GetAxis(kPidCentrality)->GetNbins()) {
      actualLowerCentrality = hPIDdata->GetAxis(kPidCentrality)->GetBinLowEdge(lowerCentralityBinLimit);
      actualUpperCentrality = hPIDdata->GetAxis(kPidCentrality)->GetBinUpEdge(upperCentralityBinLimit);

      restrictCentralityAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested centrality range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  std::cout << "centrality: ";
  if (restrictCentralityAxis) {
    std::cout << actualLowerCentrality << " - " << actualUpperCentrality << std::endl;
     hPIDdata->GetAxis(kPidCentrality)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
  }
  else {
    std::cout << "All" << std::endl;
  }
  
  // Obtain MC information about particle yields
  hPIDdata->GetAxis(kSelectSpecies)->SetRange(1, 1); // Do not count each particle more than once
  
  hPIDdata->GetAxis(kPidMCpid)->SetRange(3, 3);
  TH1D* hMCmuonYield = (TH1D*)hPIDdata->Projection(kPidPt, "e");
  hMCmuonYield->SetName("hMCmuonYield");
  
  TCanvas* cYields = new TCanvas;
  cYields->SetLogx(kTRUE);
  cYields->SetLogy(kTRUE);
  hMCmuonYield->SetLineColor(kOrange);
  hMCmuonYield->GetXaxis()->SetRangeUser(0.15, 50.);
  hMCmuonYield->Draw();
  
  hPIDdata->GetAxis(kPidMCpid)->SetRange(1, 1);
  TH1D* hMCelectronYield = (TH1D*)hPIDdata->Projection(kPidPt, "e");
  hMCelectronYield->SetName("hMCelectronYield");
  
  hMCelectronYield->SetLineColor(kMagenta);
  hMCelectronYield->Draw("same");
  
  TH1D* hMCmuonToElRatio = new TH1D(*hMCmuonYield);
  hMCmuonToElRatio->Divide(hMCelectronYield);
  
  TCanvas* cFit = new TCanvas;
  cFit->SetLogx(kTRUE);
  hMCmuonToElRatio->GetYaxis()->SetRangeUser(0., 1.);
  hMCmuonToElRatio->Draw();
  
  TF1 f("f", "[0]+[1]/TMath::Min(x, [4])+[2]*TMath::Min(x, [4])+[3]*TMath::Min(x, [4])*TMath::Min(x, [4])+[5]*TMath::Min(x, [4])*TMath::Min(x, [4])*TMath::Min(x, [4])+[6]*(x>[7])*TMath::Min(x-[7], [8]-[7])",
        0.01, 50.);
  
  f.SetParameters(-0.688, 0.042, 4.52, -6.2, 0.53, 2.89, 0.035, 2.3, 6.);
  restrictJetPtAxis = kTRUE;
  if (jetHandling)
    f.FixParameter(8, 6.);
  
  f.SetParLimits(4, 0.2, 2.0);
  f.SetParLimits(7, 1.0, 6.0);
  
  hMCmuonToElRatio->Fit(&f, jetHandling ? "+" : "+W", "", 0.15, restrictJetPtAxis ? actualUpperJetPt : 15.0);
  
  return 0;
}