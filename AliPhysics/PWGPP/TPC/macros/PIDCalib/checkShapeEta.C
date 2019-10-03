#include "THnSparse.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TVirtualFitter.h"
#include "TList.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMath.h"
#include "TDatime.h"

#include <iostream>
#include <iomanip>

#include "THnSparseDefinitions.h"
#include "ProgressBar.h"

Bool_t correctedData = kFALSE;

//--------------------------------------------------
Double_t getMedianOfNonZeros(Double_t* input, const Int_t dim)
{
  Double_t values[dim];
  for (Int_t i = 0; i < dim; i++)
    values[i] = 0.0;
    
  Int_t numNonZero = 0;
  
  for (Int_t i = 0; i < dim; i++) {
    if (input[i] > 0) {
      values[numNonZero] = input[i];
      numNonZero++;
    }
  }
  
  return ((numNonZero > 0) ? TMath::Median(numNonZero, values) : 0);
}


//--------------------------------------------------
void prepareHisto(TH2* h)
{
  for (Int_t binX = 1; binX <= h->GetXaxis()->GetNbins(); binX++) {
    
    // Find maximum for fixed x
    Double_t maxBinContent = -1;
    Int_t maxBin = -1;
    for (Int_t binY = 1; binY <= h->GetYaxis()->GetNbins(); binY++) {
      Double_t cont = h->GetBinContent(binX, binY);
      if (cont > maxBinContent) {
        maxBinContent = cont;
        maxBin = binY;
      }
    }
    
    Double_t prevBinContent = -1;
    Double_t currBinContent = -1;
    
    // Start at the maximum, go to the left and remove all bins that fall too steeply (due to TPC cut) by setting large errors
    for (Int_t binY = maxBin; binY >= 1; binY--) {
      currBinContent = h->GetBinContent(binX, binY);
      
      if (currBinContent > 0)   {
        Double_t ratio = prevBinContent / currBinContent;
        
        if (ratio > 5)    {
          for (Int_t binY2 = binY; binY2 >= 1; binY2--) {
            if (h->GetBinContent(binX, binY2) > 0) {
              h->SetBinContent(binX, binY2, 0);
              h->SetBinError(binX, binY2, maxBinContent);
            }
          }
          break;
        }
      }
      
      prevBinContent = currBinContent;
    }
    
    prevBinContent = -1;
    currBinContent = -1;
    
    // Start at the maximum, go to the right and remove all bins that fall too steeply (due to TPC cut) by setting large errors
    for (Int_t binY = maxBin; binY <= h->GetYaxis()->GetNbins(); binY++) {
      currBinContent = h->GetBinContent(binX, binY);
      
      if (currBinContent > 0)   {
        Double_t ratio = prevBinContent / currBinContent;
        
        if (ratio > 5)    {
          for (Int_t binY2 = binY; binY2 <= h->GetYaxis()->GetNbins(); binY2++) {
            if (h->GetBinContent(binX, binY2) > 0) {
              h->SetBinContent(binX, binY2, 0);
              h->SetBinError(binX, binY2, maxBinContent);
            }
          }
          break;
        }
      }
      
      prevBinContent = currBinContent;
    }    
  }
}


//--------------------------------------------------
Bool_t FitHist(TH1 *h, Double_t heightFractionForRange, TString fitOption, Double_t *results, Double_t *resultErrors)
{
  if (!h) return kFALSE;
  if (!results || !resultErrors) return kFALSE;
  
  // Average around maximum bin -> Might catch some outliers
  Int_t maxBin = h->GetMaximumBin();
  Double_t maxVal = h->GetBinContent(maxBin);
  
  if (maxVal < 2) { // It could happen that all entries are in overflow/underflow; don't fit in this case
    return kFALSE;
  }
  
  UChar_t usedBins = 1;
  if (maxBin > 1) {
    maxVal += h->GetBinContent(maxBin - 1);
    usedBins++;
  }
  if (maxBin < h->GetNbinsX()) {
    maxVal += h->GetBinContent(maxBin + 1);
    usedBins++;
  }
  maxVal /= usedBins;
  
  Double_t thresholdFraction = heightFractionForRange * maxVal; 
  Int_t lowThrBin = TMath::Max(1, h->FindFirstBinAbove(thresholdFraction));
  Int_t highThrBin = TMath::Min(h->GetNbinsX(), h->FindLastBinAbove(thresholdFraction));
  
  Double_t lowThreshold = h->GetBinCenter(lowThrBin);
  Double_t highThreshold = h->GetBinCenter(highThrBin);
  
  TFitResultPtr res = h->Fit("gaus", Form("%sS", fitOption.Data()), "", lowThreshold, highThreshold);
  
  if ((Int_t)res == 0) {
    for (Int_t i = 0; i < 3; i++) {
      results[i] = res->GetParams()[i];
      resultErrors[i] = res->GetErrors()[i];
    }
    results[3] = res->Ndf()>0 ? res->Chi2()/res->Ndf() : 0;
    resultErrors[3] = 0;
    
    return kTRUE;
  }
  
  return kFALSE;
}


//--------------------------------------------------
void FitSlicesY(TH2 *hist, Double_t heightFractionForRange, Int_t cutThreshold, TString fitOption, TObjArray *arr)
{
  if (!hist) return;
  if (!arr) return;

  // If old style is to be used
  /*
  hist->FitSlicesY(0, 0, -1, cutThreshold, fitOption.Data(), &array);
  return;
  */
  
  
  arr->Clear();
  arr->SetOwner();
  arr->Expand(4);
  
  TAxis *axis=hist->GetXaxis();
  const TArrayD *bins = axis->GetXbins();
  TH1D** hList = new TH1D*[4];
  
  for (Int_t i = 0; i < 4; i++) {
    delete gDirectory->FindObject(Form("%s_%d", hist->GetName(), i));
    
    if (bins->fN == 0) {
      hList[i] = new TH1D(Form("%s_%d", hist->GetName(), i), i < 3 ? Form("Fitted value of par[%d]", i) : "Chi2/NDF",
                          hist->GetNbinsX(), axis->GetXmin(), axis->GetXmax());
    } else {
      hList[i] = new TH1D(Form("%s_%d", hist->GetName(), i), i < 3 ? Form("Fitted value of par[%d]", i) : "Chi2/NDF", hist->GetNbinsX(), bins->fArray);
    }
    
    (*arr)[i] = hList[i];
  }
  
  Double_t results[4] = {0.0, 0.0, 0.0, 0.0 };
  Double_t resultErrors[4] = {0.0, 0.0, 0.0, 0.0 };
  
  for (Int_t ibin=0, ibin2=axis->GetFirst(); ibin2<=axis->GetLast(); ++ibin2, ++ibin) {
    TH1 *h=hist->ProjectionY("_temp",ibin2,ibin2);
    if (!h)
      continue;
    
    if (h->GetEntries() < cutThreshold) {
      delete h;
      continue;
    }
    
    Bool_t fitSuccessful = FitHist(h, heightFractionForRange, fitOption, results, resultErrors);
    
    if (fitSuccessful) {
      Int_t resBin = ibin2;
      hList[0]->SetBinContent(resBin,results[0]);
      hList[0]->SetBinError(resBin,resultErrors[0]);
      hList[1]->SetBinContent(resBin,results[1]);
      hList[1]->SetBinError(resBin,resultErrors[1]);
      hList[2]->SetBinContent(resBin,results[2]);
      hList[2]->SetBinError(resBin,resultErrors[2]);
      hList[3]->SetBinContent(resBin,results[3]);
    }
    
    delete h;
  }
}


//--------------------------------------------------
/*Obsolete: dEdx is NOT corrected!
TH1D* getSeparation(THnSparse* hist, Bool_t fromV0s, TFile* fSave)
{
  TObjArray arr;
  
  Int_t binIDel = 1;
  if (fromV0s)
    binIDel = 7;
  
  Int_t binIDpi = 4;
  if (fromV0s)
    binIDpi = 8;
  
  //hist->GetAxis(kEta)->SetRangeUser(-0.199, 0.199);
  
  // Select and fit electrons
  hist->GetAxis(kMCpid)->SetRange(binIDel,binIDel); // (V0) electrons
  hist->GetAxis(kSelectSpecies)->SetRange(1,1); // Delta_electron
  
  Int_t referenceAxis = kDeDx;
  
  TH2D* histEl = (TH2D*)hist->Projection(referenceAxis, kPtpcInner);
  histEl->SetName(fromV0s ? "histV0El" : "histEl");
  histEl->GetXaxis()->SetTitle(hist->GetAxis(kPtpcInner)->GetTitle());
  histEl->GetYaxis()->SetTitle(hist->GetAxis(referenceAxis)->GetTitle());
  
  if (!fromV0s)
    prepareHisto(histEl);
  histEl->FitSlicesY(0,0,-1,0,"QNR",&arr); 
  TH1D* hMeanEl = (TH1D*)arr.At(1)->Clone(fromV0s ? "hMeanV0El" : "hMeanEl");
  TH1D* hSigmaEl = (TH1D*)arr.At(2)->Clone(fromV0s ? "hSigmaV0El" : "hSigmaEl");
  TH1D* hChi2El = (TH1D*)arr.At(3)->Clone(fromV0s ? "hChi2V0El" : "hChi2El");

  hMeanEl->GetXaxis()->SetTitle(hist->GetAxis(kPtpcInner)->GetTitle());
  hMeanEl->GetYaxis()->SetTitle("Mean (Gauss)");
  hSigmaEl->GetXaxis()->SetTitle(hist->GetAxis(kPtpcInner)->GetTitle());
  hSigmaEl->GetYaxis()->SetTitle("#sigma (Gauss)");
  hChi2El->GetXaxis()->SetTitle(hist->GetAxis(kPtpcInner)->GetTitle());
  hChi2El->GetYaxis()->SetTitle("#chi^{2}/NDF (Gauss)");
  
  hMeanEl->SetMarkerStyle(20);
  hSigmaEl->SetMarkerStyle(20);
  
  if (fSave)    {
    fSave->cd();
    
    histEl->Write();
  }
  
  delete histEl;

  // Select and fit pions
  hist->GetAxis(kMCpid)->SetRange(binIDpi,binIDpi); // (V0) pions
  hist->GetAxis(kSelectSpecies)->SetRange(3,3); // Delta_pion
  
  TH2D* histPi = (TH2D*)hist->Projection(referenceAxis, kPtpcInner);
  histPi->SetName(fromV0s ? "histV0Pi" : "histPi");
  histPi->GetXaxis()->SetTitle(hist->GetAxis(kPtpcInner)->GetTitle());
  histPi->GetYaxis()->SetTitle(hist->GetAxis(referenceAxis)->GetTitle());
  
  if (!fromV0s)
    prepareHisto(histPi);
  histPi->FitSlicesY(0,0,-1,0,"QNR",&arr);
  TH1D* hMeanPi = (TH1D*)arr.At(1)->Clone(fromV0s ? "hMeanV0Pi" : "hMeanPi");
  TH1D* hSigmaPi = (TH1D*)arr.At(2)->Clone(fromV0s ? "hSigmaV0Pi" : "hSigmaPi");
  TH1D* hChi2Pi = (TH1D*)arr.At(3)->Clone(fromV0s ? "hChi2V0Pi" : "hChi2Pi");
  
  hMeanPi->GetXaxis()->SetTitle(hist->GetAxis(kPtpcInner)->GetTitle());
  hMeanPi->GetYaxis()->SetTitle("Mean (Gauss)");
  hSigmaPi->GetXaxis()->SetTitle(hist->GetAxis(kPtpcInner)->GetTitle());
  hSigmaPi->GetYaxis()->SetTitle("#sigma (Gauss)");
  hChi2Pi->GetXaxis()->SetTitle(hist->GetAxis(kPtpcInner)->GetTitle());
  hChi2Pi->GetYaxis()->SetTitle("#chi^{2}/NDF (Gauss)");
  
  hMeanPi->SetMarkerStyle(20);
  hMeanPi->SetLineColor(kRed);
  hMeanPi->SetMarkerColor(kRed);
  hSigmaPi->SetMarkerStyle(20);
  hSigmaPi->SetLineColor(kRed);
  hSigmaPi->SetMarkerColor(kRed);
  hChi2Pi->SetMarkerStyle(20);
  hChi2Pi->SetLineColor(kRed);
  hChi2Pi->SetMarkerColor(kRed);
  
  // Separation
  
  TH1D* hSeparation= (TH1D*)hMeanEl->Clone(fromV0s ? "hSeparationV0" : "hSeparation"); //to get same binning
  hSeparation->SetMarkerStyle(20);
  hSeparation->GetYaxis()->SetTitle("Separation");

  const Int_t nBins = hMeanEl->GetNbinsX();
  
  Double_t deltaMean[nBins] ;
  Double_t deltaSigma[nBins];
  
  for(Int_t i = 0 ;i < nBins; i++) {
    deltaMean[i] = TMath::Abs(hMeanEl->GetBinContent(i) - hMeanPi->GetBinContent(i));
    deltaSigma[i] = TMath::Abs((hSigmaEl->GetBinContent(i) + hSigmaPi->GetBinContent(i))) / 2.;
    
    if(TMath::Abs(deltaSigma[i]) < 0.000001)
      continue;

    hSeparation->SetBinContent(i, deltaMean[i] / deltaSigma[i]);
  }


  // Reset ranges
  hist->GetAxis(kMCpid)->SetRange(0,-1); 
  hist->GetAxis(kSelectSpecies)->SetRange(0,-1);
  hist->GetAxis(kEta)->SetRange(0, -1);

  if (fSave)    {
    fSave->cd();

    histPi->Write();
    
    hMeanEl->Write();
    hMeanPi->Write();

    hSigmaEl->Write();
    hSigmaPi->Write();
    
    hChi2El->Write();
    hChi2Pi->Write();
    
    hSeparation->Write();
  }
  
  delete histPi;

  return hSeparation;
}
*/

//--------------------------------------------------
Int_t checkShapeEta(TString path = ".", Double_t multiplicityStepSize = -1/*-1 for all multiplicities*/,
                    Double_t maxMultiplicity = 20000,
                    Int_t onlyEtaShapeComparison = 0, TString fileName = "bhess_PIDetaAdv.root", 
                    TString listName = "bhess_PIDetaAdv", Int_t momentumAxis = kPtpcInner/* = kPt*/, Int_t binTypePt = kPtBinTypePPMult) { 
			       
  Int_t nPtBins = -1;
  const Double_t* binsPt = GetPtBins(binTypePt, nPtBins);
  
  Double_t binWidthsPt[nPtBins];
  for (Int_t i = 0; i < nPtBins; i++) 
    binWidthsPt[i] = (binsPt[i + 1] - binsPt[i]) / 2.;

  Int_t cutForFitting = 10;
  Double_t heightFractionForFittingRange = 0.1;
  
  TList* histList = 0x0;
  
	TFile* f = TFile::Open(Form("%s/%s", path.Data(), fileName.Data()));
  if (!f)  {
    std::cout << "Failed to open file \"" << Form("%s/%s", path.Data(), fileName.Data()) << "\"!" << std::endl;
    return -1;
  }
  
  THnSparse* hPIDdata = 0x0;
  
  if (correctedData) {
    hPIDdata = dynamic_cast<THnSparse*>(f->Get("hPIDdataCorrected"));
    if (!hPIDdata) {
      std::cout << "Failed to load (extracted) data histo!" << std::endl;
      return -1;
    }
  }
  else  {
    histList = (TList*)(f->Get(listName.Data()));
    if (!histList) {
      std::cout << "Failed to load list \"" << listName.Data() << "\"!" << std::endl;
      return -1;
    }
    
    // Extract the data histogram
    hPIDdata = dynamic_cast<THnSparse*>(histList->FindObject("hPIDdataAll"));
    if (!hPIDdata) {
      hPIDdata = dynamic_cast<THnSparse*>(f->Get(Form("%s/hPIDdataAll", listName.Data())));
      if (!hPIDdata)  {
        std::cout << "Failed to load data histo!" << std::endl;
        return -1;
      }
    }
  }
  
  // Set proper errors
  hPIDdata->Sumw2();
  Long64_t nBinsTHnSparse = hPIDdata->GetNbins();
  Double_t binContent = 0;
  
  for (Long64_t bin = 0; bin < nBinsTHnSparse; bin++) {
    binContent = hPIDdata->GetBinContent(bin);
    hPIDdata->SetBinError(bin, TMath::Sqrt(binContent));
  }
  
  // Output file
  TDatime daTime;
  TString saveFileName = Form("outputCheckShapeEtaAdv_%04d_%02d_%02d__%02d_%02d.root", daTime.GetYear(), daTime.GetMonth(), 
                              daTime.GetDay(), daTime.GetHour(), daTime.GetMinute());
            
  TFile* fSave = TFile::Open(Form("%s/%s", path.Data(), saveFileName.Data()), "recreate");
  if (!fSave) {
    std::cout << "Failed to open save file \"" << Form("%s/%s", path.Data(), saveFileName.Data()) << "\"!" << std::endl;
    return -1;
  }
  
  // p dependece of eta dependence
  Bool_t first = kTRUE;
  
  Int_t momLow = 1;// = 1 or 7 or 21
  Int_t momHigh = hPIDdata->GetAxis(momentumAxis)->GetNbins();//;24 or 40 or hPIDdata->GetAxis(momentumAxis)->GetNbins()
  std::cout << "Momentum range: ";
  std::cout << hPIDdata->GetAxis(momentumAxis)->GetBinLowEdge(momLow) << " - " << hPIDdata->GetAxis(momentumAxis)->GetBinUpEdge(momHigh);
  std::cout << std::endl;
  
  //TF1* constantFunc = new TF1("constantFunc", "1", -2, +2);
  
  const Int_t nModes = 3;
  TString mode[nModes] = {"dE/dx", "Delta_{Species}", "Delta'_{Species}"};
  Int_t dataAxis[nModes] = {0/*kDeDx*/, 0/*kDelta*/, kDeltaPrime};
  
  
  
  if (onlyEtaShapeComparison > 0) {
    // Check, if there is a dependence on the species and/or p directly by choosing momentum bins with equal dE/dx for
    // different species and plotting Delta vs eta etc.  
    
    
    for (onlyEtaShapeComparison = 1; onlyEtaShapeComparison <= 3; onlyEtaShapeComparison++) { 
      TCanvas* canvSpecific = new TCanvas(Form("EtaDependenceOfDifferentSpeciesAtSamedEdx_%d", onlyEtaShapeComparison), 
                                          "#eta dependence of different species at same dE/dx",
                                          100, 10, 1380 / 2, 700);
      
      canvSpecific->SetRightMargin(0.001);
      canvSpecific->SetTopMargin(0.01);
      canvSpecific->SetLeftMargin(0.18);
      canvSpecific->SetBottomMargin(0.11);
      Int_t speciesColour[6] = { 4, 1, 3, 2, 6, 7 };
      
      Int_t momPions = -1;
      //Int_t momPions2 = hPIDdata->GetAxis(momentumAxis)->FindBin(4.0);
      //Int_t momPions3 = hPIDdata->GetAxis(momentumAxis)->FindBin(7.0);
      Int_t momElectrons = -1; 
      
      Int_t momProtons = -1;
      Int_t momKaons = -1;
      
      Int_t dEdx = -1;
      
      if (onlyEtaShapeComparison == 1)  {// For dE/dx about 50
        momProtons = hPIDdata->GetAxis(momentumAxis)->FindBin(2.7);
        momKaons = hPIDdata->GetAxis(momentumAxis)->FindBin(1.3);
        momPions = hPIDdata->GetAxis(momentumAxis)->FindBin(0.45);
        dEdx = 50;
      }
      else if (onlyEtaShapeComparison == 2)  {// For dE/dx about 75
        momProtons = hPIDdata->GetAxis(momentumAxis)->FindBin(1.0);
        momKaons = hPIDdata->GetAxis(momentumAxis)->FindBin(0.5);
        momElectrons = hPIDdata->GetAxis(momentumAxis)->FindBin(0.6);
        dEdx = 75;
      }
      else if (onlyEtaShapeComparison == 3)  {// For dE/dx about 60
        momProtons = hPIDdata->GetAxis(momentumAxis)->FindBin(1.4);
        momKaons = hPIDdata->GetAxis(momentumAxis)->FindBin(0.735);
        // Sample not clean and V0 pi don't have sufficient statistics momPions = hPIDdata->GetAxis(momentumAxis)->FindBin(3.9);
        dEdx = 60;
      }
      else if (onlyEtaShapeComparison == 4)  {// For dE/dx about 100
        momProtons = hPIDdata->GetAxis(momentumAxis)->FindBin(0.8);
        momKaons = hPIDdata->GetAxis(momentumAxis)->FindBin(0.4);
        dEdx = 100;
      }
      else  {
        std::cout << "Eta shape comparison not defined for this value!" << std::endl;
        return -1;
      }
      
      
      /*
      canvSpecific->Divide(3,1, 0.01, 0.01);
      
      Double_t x1, x2, y1, y2;
      
      for (Int_t  i = 1; i <= 3; i++) {
        x1 = 0. + (i-1) * 1./3.;
        y1 = 0.85;
        x2 = x1 + 1./3.;
        y2 = 0.; 
            
        canvSpecific->GetPad(i)->SetPad(x1, y1, x2, y2);
      }*/
      
      
      /*Definitions from ADV task (18.04.14)
      hist->GetAxis(0)->SetBinLabel(1, "e");
      hist->GetAxis(0)->SetBinLabel(2, "K");
      hist->GetAxis(0)->SetBinLabel(3, "#mu");
      hist->GetAxis(0)->SetBinLabel(4, "#pi");
      hist->GetAxis(0)->SetBinLabel(5, "p");
      hist->GetAxis(0)->SetBinLabel(6, "V0+TOF e");
      hist->GetAxis(0)->SetBinLabel(7, "V0 e");
      hist->GetAxis(0)->SetBinLabel(8, "V0 #pi");
      hist->GetAxis(0)->SetBinLabel(9, "V0 p");

      hist->GetAxis(1)->SetTitle("Select Species");
      hist->GetAxis(1)->SetBinLabel(1, "e");
      hist->GetAxis(1)->SetBinLabel(2, "K");
      hist->GetAxis(1)->SetBinLabel(3, "#pi");
      hist->GetAxis(1)->SetBinLabel(4, "p");
      */
      
      TLegend* legendSpecific = new TLegend(0.62, 0.76, 0.95, 0.95);    
      //legendSpecific->SetNColumns(2);
      legendSpecific->SetBorderSize(0);
      legendSpecific->SetFillColor(0);
        
      Bool_t isFirst = kTRUE;
      
      Int_t i = 2; //Delta'
      
      for (Int_t species = 1; species <= 4/*6*/; species++) {
        TString histTitle = "";
        // Select particle species and corresponding deltaSpecies
        if (species == 1) {// Electrons
          if (momElectrons < 0)
            continue;
          hPIDdata->GetAxis(kMCpid)->SetRange(1,1); // Only particles of species X
          hPIDdata->GetAxis(momentumAxis)->SetRange(momElectrons, momElectrons);
          hPIDdata->GetAxis(kSelectSpecies)->SetRange(1,1); 
          histTitle = Form("e, %.2f #leq p (GeV/c) #leq %.2f", 
                            hPIDdata->GetAxis(momentumAxis)->GetBinLowEdge(momElectrons),
                            hPIDdata->GetAxis(momentumAxis)->GetBinUpEdge(momElectrons));
        }
        else if (species == 2) {// Kaons
          hPIDdata->GetAxis(kMCpid)->SetRange(2,2); // Only particles of species X
          hPIDdata->GetAxis(momentumAxis)->SetRange(momKaons, momKaons);
          hPIDdata->GetAxis(kSelectSpecies)->SetRange(2,2); 
          histTitle = Form("K, %.2f #leq p (GeV/c) #leq %.2f", 
                            hPIDdata->GetAxis(momentumAxis)->GetBinLowEdge(momKaons),
                            hPIDdata->GetAxis(momentumAxis)->GetBinUpEdge(momKaons));
        }
        else if (species == 3) {// Pions
          if (momPions < 0)
            continue;
          hPIDdata->GetAxis(kMCpid)->SetRange(4,4); // Only particles of species X
          hPIDdata->GetAxis(momentumAxis)->SetRange(momPions, momPions);
          hPIDdata->GetAxis(kSelectSpecies)->SetRange(3,3);
          histTitle = Form("#pi, %.2f #leq p (GeV/c) #leq %.2f", 
                            hPIDdata->GetAxis(momentumAxis)->GetBinLowEdge(momPions),
                            hPIDdata->GetAxis(momentumAxis)->GetBinUpEdge(momPions));
        }
        else if (species == 4) {// Protons
          hPIDdata->GetAxis(kMCpid)->SetRange(5,5); // Only particles of species X
          hPIDdata->GetAxis(momentumAxis)->SetRange(momProtons, momProtons);
          hPIDdata->GetAxis(kSelectSpecies)->SetRange(4,4);
          histTitle = Form("p, %.2f #leq p (GeV/c) #leq %.2f", 
                            hPIDdata->GetAxis(momentumAxis)->GetBinLowEdge(momProtons),
                            hPIDdata->GetAxis(momentumAxis)->GetBinUpEdge(momProtons));
        }/*
        else if (species == 5) {// Pions at higher momentum
          hPIDdata->GetAxis(kMCpid)->SetRange(4,4); // Only particles of species X
          hPIDdata->GetAxis(momentumAxis)->SetRange(momPions2, momPions2);
          hPIDdata->GetAxis(kSelectSpecies)->SetRange(3,3);
        }
        else if (species == 6) {// Pions at even higher momentum
          hPIDdata->GetAxis(kMCpid)->SetRange(4,4); // Only particles of species X
          hPIDdata->GetAxis(momentumAxis)->SetRange(momPions3, momPions3);
          hPIDdata->GetAxis(kSelectSpecies)->SetRange(3,3);
        }*/
        
        
        TH2D* hEta = hPIDdata->Projection(dataAxis[i], kEta);
        hEta->Rebin2D(2, 1);
        hEta->SetName(Form("hEta_%s_%d_%d", mode[i].Data(), onlyEtaShapeComparison, species));
        hEta->SetTitle(Form("Eta dependence of %s for %s", mode[i].Data(),  hPIDdata->GetAxis(kMCpid)->GetBinLabel(species)));
      
        TObjArray aSlices;
        FitSlicesY(hEta, heightFractionForFittingRange, cutForFitting, "QNR", &aSlices);
        TH1D* hEtaProj = (TH1D*)(aSlices.At(1)->Clone(Form("hEtaProj_%s_%d_%d", mode[i].Data(), onlyEtaShapeComparison, species)));
        hEtaProj->SetMarkerStyle(20 + species - 1);
        hEtaProj->SetTitle(histTitle.Data());
        hEtaProj->SetLineColor(speciesColour[species - 1]);
        hEtaProj->SetMarkerColor(speciesColour[species - 1]);
        
        hEtaProj->GetXaxis()->SetTitle(hPIDdata->GetAxis(kEta)->GetTitle());
        
        // Scale to unity
        Double_t integral = 0;
        Double_t unityIntegral = 0;
        for (Int_t j = 1; j <= hEtaProj->GetNbinsX(); j++) {
          if (hEtaProj->GetBinContent(j) > 0) {
            // Remove bins with too large error
            if (hEtaProj->GetBinError(j) > 0.005) {
              hEtaProj->SetBinContent(j, -1);
              continue;
            }
            
            unityIntegral += 1;//hEtaProj->GetXaxis()->GetBinWidth(j);
            integral += hEtaProj->GetBinContent(j);
          }
        }
        hEtaProj->Scale(unityIntegral/integral);
        
        hEtaProj->SetStats(kFALSE);
        hEtaProj->GetYaxis()->SetRangeUser(0.97,1.059);
        hEtaProj->GetYaxis()->SetTitle("#Delta'_{#lower[-0.5]{scaled}}");
        hEtaProj->GetYaxis()->SetTitleOffset(1.4);
        
        hEtaProj->DrawCopy(Form("%s", (isFirst ? "" : "same")));
        isFirst = kFALSE;
        
        legendSpecific->AddEntry(hEtaProj, hEtaProj->GetTitle(), "flp");
        
        //hEtaProj->Fit("pol3", "+", "same", -1, 0);
        //hEtaProj->Fit("pol3", "+", "same", 0, 1);
        //((TF1*)hEtaProj->GetListOfFunctions()->At(0))->SetLineColor(species);
        //((TF1*)hEtaProj->GetListOfFunctions()->At(1))->SetLineColor(species);
      }
      
      legendSpecific->Draw();
      
      ClearTitleFromHistoInCanvas(canvSpecific);
      
      //canvSpecific->cd(0);
      /*
      TLegend* legendSpecific = new TLegend(0.1, 0.85, 0.9, 0.99);    
      legendSpecific->SetNColumns(2);
      legendSpecific->SetBorderSize(0);
      legendSpecific->SetFillColor(0);
      
      
      TLegendEntry* ent = legendSpecific->AddEntry((TObject*)0, path.Data(), "");
      ent->SetTextColor(1);
      
      ent = legendSpecific->AddEntry((TObject*)0, "", "");
      
      ent = legendSpecific->AddEntry((TObject*)0, Form("K, %.2f #leq p (GeV/c) #leq %.2f", 
                                                      hPIDdata->GetAxis(momentumAxis)->GetBinLowEdge(momKaons),
                                                      hPIDdata->GetAxis(momentumAxis)->GetBinUpEdge(momKaons)),
                                    "p");
      ent->SetLineColor(speciesColour[1]);
      ent->SetTextColor(speciesColour[1]);
      ent->SetMarkerColor(speciesColour[1]);
      ent->SetMarkerStyle(2);
      ent->SetMarkerSize(2);
      
      ent = legendSpecific->AddEntry((TObject*)0, Form("e, %.2f #leq p (GeV/c) #leq %.2f", 
                                                      hPIDdata->GetAxis(momentumAxis)->GetBinLowEdge(momElectrons),
                                                      hPIDdata->GetAxis(momentumAxis)->GetBinUpEdge(momElectrons)),
                                    "p");
      ent->SetLineColor(speciesColour[0]);
      ent->SetTextColor(speciesColour[0]);
      ent->SetMarkerColor(speciesColour[0]);
      ent->SetMarkerStyle(2);
      ent->SetMarkerSize(2);
        
      ent = legendSpecific->AddEntry((TObject*)0, Form("p, %.2f #leq p (GeV/c) #leq %.2f", 
                                                      hPIDdata->GetAxis(momentumAxis)->GetBinLowEdge(momProtons),
                                                      hPIDdata->GetAxis(momentumAxis)->GetBinUpEdge(momProtons)),
                                    "p");
      ent->SetLineColor(speciesColour[3]);
      ent->SetTextColor(speciesColour[3]);
      ent->SetMarkerColor(speciesColour[3]);
      ent->SetMarkerStyle(2);
      ent->SetMarkerSize(2);
      
      ent = legendSpecific->AddEntry((TObject*)0, Form("#pi, %.2f #leq p (GeV/c) #leq %.2f", 
                                                      hPIDdata->GetAxis(momentumAxis)->GetBinLowEdge(momPions),
                                                      hPIDdata->GetAxis(momentumAxis)->GetBinUpEdge(momPions)),
                                    "p");
      ent->SetLineColor(speciesColour[2]);
      ent->SetTextColor(speciesColour[2]);
      ent->SetMarkerColor(speciesColour[2]);
      ent->SetMarkerStyle(2);
      ent->SetMarkerSize(2);
      
      ent = legendSpecific->AddEntry((TObject*)0, Form("V0 #pi, %.2f #leq p (GeV/c) #leq %.2f", 
                                                      hPIDdata->GetAxis(momentumAxis)->GetBinLowEdge(momPions2),
                                                      hPIDdata->GetAxis(momentumAxis)->GetBinUpEdge(momPions2)),
                                    "p");
      ent->SetLineColor(speciesColour[4]);
      ent->SetTextColor(speciesColour[4]);
      ent->SetMarkerColor(speciesColour[4]);
      ent->SetMarkerStyle(2);
      ent->SetMarkerSize(2);
      
      ent = legendSpecific->AddEntry((TObject*)0, Form("V0 #pi, %.2f #leq p (GeV/c) #leq %.2f", 
                                                      hPIDdata->GetAxis(momentumAxis)->GetBinLowEdge(momPions3),
                                                      hPIDdata->GetAxis(momentumAxis)->GetBinUpEdge(momPions3)),
                                    "p");
      ent->SetLineColor(speciesColour[5]);
      ent->SetTextColor(speciesColour[5]);
      ent->SetMarkerColor(speciesColour[5]);
      ent->SetMarkerStyle(2);
      ent->SetMarkerSize(2);
        
      legendSpecific->Draw("");
      */
      fSave->cd();
      canvSpecific->Write();
      TString pdfName = saveFileName;
      pdfName = pdfName.ReplaceAll(".root", Form("__dEdx_%d.pdf", dEdx));
      canvSpecific->SaveAs(Form("%s/%s", path.Data(), pdfName.Data()));
    }
    
    return 0;
  }
  
  // Reset ranges
  hPIDdata->GetAxis(kMCpid)->SetRange(0, -1);
  hPIDdata->GetAxis(momentumAxis)->SetRange(0, -1);
  hPIDdata->GetAxis(kSelectSpecies)->SetRange(0, -1);
  
  
  
  
  
  
  /*
  getSeparation(hPIDdata, kFALSE, fSave);
  getSeparation(hPIDdata, kTRUE, fSave);
  return 0;*/
  
  // If V0s are available in histo, add them! Also: Backward compatibility to data type with "all id. primaries" bin instead of "V0plusTOF el"
  Int_t additionalV0bins = 0;
  
  Bool_t oldStyleWithAllIDprimaries = strcmp(hPIDdata->GetAxis(kMCpid)->GetBinLabel(6), "all id. primaries") == 0;
  printf("Backward compatibility: Checking bin label of bin 6 of MCpid axis: %s -> %s\n",
         hPIDdata->GetAxis(kMCpid)->GetBinLabel(6), oldStyleWithAllIDprimaries ? "Old style" : "New style");
  
  if (hPIDdata->GetAxis(kMCpid)->GetNbins() >= 8)
    additionalV0bins = oldStyleWithAllIDprimaries ? 3 : 4;
  
  // If valid multiplicityStepSize, fill each multiplicity bin + 1 bin with all multiplicities. Otherwise: Only bin with all multiplicities
  const Int_t nMultBins = (multiplicityStepSize <= 0) ? 1 :
                          (((Int_t)((Double_t)maxMultiplicity / multiplicityStepSize)) + (((Int_t)maxMultiplicity % TMath::Max(1, (Int_t)multiplicityStepSize)) != 0) + 1);
  
  const Int_t nSpeciesBins = (oldStyleWithAllIDprimaries ? 6 : 5) + additionalV0bins;
  TH2D* hPEtaDependence[nMultBins][nSpeciesBins];
  
  TH1D* hPMeandEdxDependence[nMultBins][nSpeciesBins];
  
  const Int_t numTotalBins = (momHigh - momLow + 1) * (4 + additionalV0bins) * nMultBins; // nMomBins * nSpecies * nMultBins
  
  for (Int_t i = 0; i < nMultBins; i++) {
    TString currDir = (i == nMultBins - 1) ? "AllMultiplicites" : 
                      Form("Multiplicity%g_%g", i * multiplicityStepSize, TMath::Min(maxMultiplicity, (i + 1) * multiplicityStepSize));
    
    TString multTitle = (i == nMultBins - 1) ? "all multiplicites" :
                        Form("multiplicity %g - %g", i * multiplicityStepSize, TMath::Min(maxMultiplicity, (i + 1) * multiplicityStepSize));
                        
    fSave->cd();
    fSave->mkdir(currDir.Data());
    // nSpeciesBins + 1 to take into account bin "all primaries"
    for (Int_t species = 1, deltaSpecies = 1; species < nSpeciesBins + 1; species++, deltaSpecies++) {
      if (species == 3) {
        deltaSpecies--;
        continue; // skip muons
      }
      else if (species == 6)    {
        if (oldStyleWithAllIDprimaries) {
          continue; // skip bin "all primaries"
        }
        else {
          deltaSpecies = 1; // V0+TOF el
        }
      }
      else if (species == 7)    {
        deltaSpecies = 1; // V0 el
      }
      else if (species == 8)    {
        deltaSpecies = 3; // V0 pi
      }
      else if (species == 9)    {
        deltaSpecies = 4; // V0 pr
      }
      
      first = kTRUE;
      
      hPMeandEdxDependence[i][species - 1] = new TH1D(Form("hPMeandEdxDependenceDeltaPrime_%s_%s", partShortName[species - 1].Data(), currDir.Data()), 
                                                      Form("Momentum dependence of mean dEdx (Delta'_{Species}) of %s, %s",
                                                           partShortName[species - 1].Data(), multTitle.Data()),
                                                           nPtBins, binsPt); 
      
      hPMeandEdxDependence[i][species - 1]->GetXaxis()->SetTitle(hPIDdata->GetAxis(momentumAxis)->GetTitle());
      hPMeandEdxDependence[i][species - 1]->GetYaxis()->SetTitle("#Delta'");
      
      
      
      
      hPEtaDependence[i][species - 1] = new TH2D(Form("hPEtaDependence_species_DeltaPrime_%s_%s", partShortName[species - 1].Data(), currDir.Data()), 
                                                 Form("Momentum dependence of eta dependence (#Delta'_{Species}) of %s, %s",
                                                      partShortName[species - 1].Data(), multTitle.Data()),
                                                 nPtBins, binsPt, 
                                                 hPIDdata->GetAxis(kEta)->GetNbins(), hPIDdata->GetAxis(kEta)->GetBinLowEdge(1),
                                                 hPIDdata->GetAxis(kEta)->GetBinUpEdge(hPIDdata->GetAxis(kEta)->GetNbins())); 
      
      hPEtaDependence[i][species - 1]->GetXaxis()->SetTitle(hPIDdata->GetAxis(momentumAxis)->GetTitle());
      hPEtaDependence[i][species - 1]->GetYaxis()->SetTitle(hPIDdata->GetAxis(kEta)->GetTitle());
      hPEtaDependence[i][species - 1]->GetYaxis()->SetLabelSize(0.03);
      hPEtaDependence[i][species - 1]->GetYaxis()->SetTitleOffset(0.9);
      hPEtaDependence[i][species - 1]->GetXaxis()->SetTitleOffset(1.2);
      hPEtaDependence[i][species - 1]->GetZaxis()->SetTitle("#Delta'");
      hPEtaDependence[i][species - 1]->GetZaxis()->SetTitleOffset(0.8);
      
      TH2D* hPEtaDependenceScaled = (TH2D*)hPEtaDependence[i][species - 1]->Clone(Form("%s_scaled", hPEtaDependence[i][species - 1]->GetName())); 
      hPEtaDependenceScaled->SetTitle(Form("%s (scaled)", hPEtaDependenceScaled->GetTitle()));
      
      hPIDdata->GetAxis(kMCpid)->SetRange(species,species); // Only particles of type X
      hPIDdata->GetAxis(kSelectSpecies)->SetRange(deltaSpecies,deltaSpecies); // Seclect corresponding deltaSpecies
      // Seclect corresponding multiplicity range (last bin (or the only bin) = all multiplicities)
      if (i == nMultBins - 1)
        hPIDdata->GetAxis(kMultiplicity)->SetRange(0, -1); // All multiplicities
      else
        hPIDdata->GetAxis(kMultiplicity)->SetRangeUser(i * multiplicityStepSize, TMath::Min(maxMultiplicity, (i + 1) * multiplicityStepSize)); 
      
      hPIDdata->GetAxis(momentumAxis)->SetRange(0, -1); // Full momentum range
      TH3D* h3Dproj = hPIDdata->Projection(momentumAxis, kDeltaPrime, kEta);
      
      for (Int_t momentum_Bin = momLow; momentum_Bin <= momHigh; momentum_Bin++)  {
        h3Dproj->GetXaxis()->SetRange(momentum_Bin, momentum_Bin); // Select p/pT bin
        
        TH1D* hMeandEdx = (TH1D*)h3Dproj->Project3D("y");
        hMeandEdx->SetName(Form("hMeandEdx_DeltaPrime_%s_%d_%d", partShortName[species - 1].Data(), momentum_Bin, i));
        hMeandEdx->SetTitle(Form("p dependence of mean dEdx of #Delta'_{Species} for %s (p %.3f), %s",
                            hPIDdata->GetAxis(kMCpid)->GetBinLabel(species), hPIDdata->GetAxis(momentumAxis)->GetBinCenter(momentum_Bin),
                            multTitle.Data()));
        
        Double_t results[4] = {0.0, 0.0, 0.0, 0.0 };
        Double_t resultErrors[4] = {0.0, 0.0, 0.0, 0.0 };
  
        Bool_t meandEdxFitSuccessful = FitHist(hMeandEdx, heightFractionForFittingRange, "QN", results, resultErrors);
        if (meandEdxFitSuccessful) {
          hPMeandEdxDependence[i][species - 1]->SetBinContent(momentum_Bin, results[1]);
          hPMeandEdxDependence[i][species - 1]->SetBinError(momentum_Bin, resultErrors[1]);
        }
        
        TH2D* hEta = (TH2D*)h3Dproj->Project3D("yz");
        hEta->SetName(Form("hEta_DeltaPrime_%d_%d_%d", species, momentum_Bin, i));
        
        hEta->SetTitle(Form("p dependence of eta dependence of #Delta' for %s (p %.3f), %s",
                            hPIDdata->GetAxis(kMCpid)->GetBinLabel(species), hPIDdata->GetAxis(momentumAxis)->GetBinCenter(momentum_Bin),
                            multTitle.Data()));
                    
        TObjArray aSlices;
        FitSlicesY(hEta, heightFractionForFittingRange, cutForFitting, "QN", &aSlices);
        
        TH1D* hEtaProj = (TH1D*)(aSlices.At(1)->Clone(Form("hEtaProj_DeltaPrime_%s_p_%.3f", partShortName[species - 1].Data(), 
                                                           hPIDdata->GetAxis(momentumAxis)->GetBinCenter(momentum_Bin))));
        
        if (!hEtaProj)  {
          std::cout << "Failed to load eta projection from FitSlicesY!" << std::endl;
          return -1;
        }
        
        first = kFALSE;
        
        for (Int_t etaBin = 1; etaBin <= hEtaProj->GetXaxis()->GetNbins(); etaBin++)  {
          hPEtaDependence[i][species - 1]->SetBinContent(momentum_Bin, etaBin, hEtaProj->GetBinContent(etaBin));
          hPEtaDependence[i][species - 1]->SetBinError(momentum_Bin, etaBin, hEtaProj->GetBinError(etaBin));
        }
                
        Int_t nBinsFinishedMomentum = (momentum_Bin - momLow + 1);
        Int_t numSpecies = species;
        if (species > 3)
          numSpecies--; // Muons not used
          
        if (oldStyleWithAllIDprimaries && species > 6)
          numSpecies--; // All primaries not used
        Int_t nBinsFinishedSpecies = (numSpecies - 1) * (momHigh - momLow + 1);
        Int_t nBinsFinishedMult = i * (4 + additionalV0bins) * (momHigh - momLow + 1);
        
        delete hMeandEdx;
        delete hEta;
        delete hEtaProj;
        
        progressbar(100. * ( ((Double_t) (nBinsFinishedMomentum + nBinsFinishedSpecies + nBinsFinishedMult)) / numTotalBins));
        
      }
      
      delete h3Dproj;
      
      hPEtaDependence[i][species - 1]->GetXaxis()->SetRangeUser(0.15, 5.0);
      hPEtaDependence[i][species - 1]->GetZaxis()->SetRangeUser(0.9, 1.1);
      
      fSave->cd(currDir.Data());
      if (hPEtaDependence[i][species - 1]) {
        hPEtaDependence[i][species - 1]->SetObjectStat(0);
        hPEtaDependence[i][species - 1]->Write();
      }
      
      if (hPMeandEdxDependence[i][species - 1]) {
        hPMeandEdxDependence[i][species - 1]->GetXaxis()->SetRangeUser(0.15, 5.0);
        hPMeandEdxDependence[i][species - 1]->GetYaxis()->SetRangeUser(0.9, 1.1);
        hPMeandEdxDependence[i][species - 1]->Write();
      }
      
      // Scale histograms 
      const Int_t nEtaBins = hPEtaDependence[i][species - 1]->GetYaxis()->GetNbins();
                              
      for (Int_t momentum_Bin = momLow; momentum_Bin <= momHigh; momentum_Bin++)  {
                
        // Scale such that radially emitted particles have deltaPrime equal 1!!
        // To be robust against fluctuations: Take median of a few eta bins around eta=0 for scaling
        Double_t values[nEtaBins];
        for (Int_t etaBin = 1; etaBin <= nEtaBins; etaBin++)
          values[etaBin - 1] = 0;
          
        for (Int_t etaBin = hPEtaDependence[i][species - 1]->GetYaxis()->FindBin(-0.175); 
                  etaBin <= hPEtaDependence[i][species - 1]->GetYaxis()->FindBin(+0.175); etaBin++)  {
          Double_t temp = hPEtaDependence[i][species - 1]->GetBinContent(momentum_Bin, etaBin);
          values[etaBin - 1] = (temp > 0) ? temp : 0;
        }
        
        Double_t temp = getMedianOfNonZeros(values, nEtaBins);
        
        if (temp <= 0) 
          continue;
        
        Double_t scale = 1. / temp;
        
        for (Int_t etaBin = 1; etaBin <= nEtaBins; etaBin++)  {
          hPEtaDependenceScaled->SetBinContent(momentum_Bin, etaBin, 
                                                scale * hPEtaDependence[i][species - 1]->GetBinContent(momentum_Bin, etaBin));
          hPEtaDependenceScaled->SetBinError(momentum_Bin, etaBin, 
                                              scale * hPEtaDependence[i][species - 1]->GetBinError(momentum_Bin, etaBin));
        }
      }
      
      hPEtaDependenceScaled->GetXaxis()->SetRangeUser(0.15, 5.0);
      fSave->cd(currDir.Data());
      hPEtaDependenceScaled->SetDrawOption("surf1");
      hPEtaDependenceScaled->Write();
      
      delete hPEtaDependenceScaled;
    }
  }
  
  // Reset ranges
  hPIDdata->GetAxis(kMCpid)->SetRange(0, -1);
  hPIDdata->GetAxis(momentumAxis)->SetRange(0, -1);
  hPIDdata->GetAxis(kSelectSpecies)->SetRange(0, -1);
  hPIDdata->GetAxis(kMultiplicity)->SetRange(0, -1);
  
  progressbar(100.);
  printf("\n");

  /*
  printf("Determining separation...\n");

  getSeparation(hPIDdata, kFALSE, fSave);
  getSeparation(hPIDdata, kTRUE, fSave);
  */
  return 0; 
}
