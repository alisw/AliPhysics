#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TDatime.h"

#include "TSpline.h"
#include "AliPID.h"

#include <iostream>
#include "ProgressBar.h"
#include "THnSparseDefinitions.h"

//________________________________________________________________________
void FitSlicesY(TH2 *hist, Double_t heightFractionForRange, Int_t cutThreshold, TString fitOption, TObjArray *arr)
{
  //heightPercentageForRange
  // custom slices fit
  //

  if (!hist) return;
  if (!arr) return;

  // If old style is to be used
  /*
  hist->FitSlicesY(0, 0, -1, cutThreshold, fitOption.Data(), arr);
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
  
  for (Int_t ibin=axis->GetFirst(); ibin<=axis->GetLast(); ++ibin){
    TH1 *h=hist->ProjectionY("_temp",ibin,ibin);
    if (!h)
      continue;
    
    if (h->GetEntries() < cutThreshold) {
      delete h;
      continue;
    }
    
    // Average around maximum bin -> Might catch some outliers
    Int_t maxBin = h->GetMaximumBin();
    Double_t maxVal = h->GetBinContent(maxBin);
    
    if (maxVal < 2) { // It could happen that all entries are in overflow/underflow; don't fit in this case
      delete h;
      continue;
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
      Int_t resBin = ibin;
      hList[0]->SetBinContent(resBin,res->GetParams()[0]);
      hList[0]->SetBinError(resBin,res->GetErrors()[0]);
      hList[1]->SetBinContent(resBin,res->GetParams()[1]);
      hList[1]->SetBinError(resBin,res->GetErrors()[1]);
      hList[2]->SetBinContent(resBin,res->GetParams()[2]);
      hList[2]->SetBinError(resBin,res->GetErrors()[2]);
      hList[3]->SetBinContent(resBin,res->Ndf()>0?res->Chi2()/res->Ndf():0);
    }
    
    delete h;
  }
  
  delete [] hList;
}

//__________________________________________________________________________________________
Int_t extractMultiplicityDependence(TString pathTree, Double_t widthFactor /*=1*/, Int_t multStepSize /*= 400*/, Bool_t isMC,
                                    Bool_t correct = kFALSE, TString fileNameTree = "bhess_PIDetaTree.root", TString treeName = "fTree",
                                    TString pathNameSplines = "splines_10h.pass2.root", TString splineName = "TSPLINE3_DATA_PROTON_LHC10H_PASS2_PBPB_MEAN") 
{ 
  const Double_t massProton = AliPID::ParticleMass(AliPID::kProton);
  const Double_t massPion = AliPID::ParticleMass(AliPID::kPion);
  
  Double_t mass = massProton;
  
  
  // Extract the data tree               
  TFile* fTree = TFile::Open(Form("%s/%s", pathTree.Data(), fileNameTree.Data()));
  if (!fTree)  {
    std::cout << "Failed to open file \"" << Form("%s/%s", pathTree.Data(), fileNameTree.Data()) << "\"!" << std::endl;
    return -1;
  }
  
  TTree* tree = dynamic_cast<TTree*>(fTree->Get(treeName.Data()));
  if (!tree) {
    std::cout << "Failed to load data tree!" << std::endl;
    return -1;
  }
  
  // Output file
  TDatime daTime;
  TString savefileName = Form("%s_multiplicityDependence_%04d_%02d_%02d__%02d_%02d.root", fileNameTree.ReplaceAll(".root", "").Data(),
                              daTime.GetYear(), daTime.GetMonth(), daTime.GetDay(), daTime.GetHour(), daTime.GetMinute());
  
  TFile* fSave = TFile::Open(Form("%s/%s", pathTree.Data(), savefileName.Data()), "recreate");
  if (!fSave) {
    std::cout << "Failed to open save file \"" << Form("%s/%s", pathTree.Data(), savefileName.Data()) << "\"!" << std::endl;
    return -1;
  }
  
  Long64_t nTreeEntries = tree->GetEntriesFast();
  
  Double_t pTPC = 0.; 
  Double_t dEdx = 0.; 
  //Double_t dEdxExpected = 0.;
  Double_t tanTheta = 0.; 
  UShort_t tpcSignalN = 0.;
  Int_t fMultiplicity = 0.; 
  UChar_t pidType = 0;
  
  tree->SetBranchStatus("*", 0); // Disable all branches
  tree->SetBranchStatus("pTPC", 1);
  tree->SetBranchStatus("dEdx", 1);
  //tree->SetBranchStatus("dEdxExpected", 1);
  tree->SetBranchStatus("tanTheta", 1);
  tree->SetBranchStatus("fMultiplicity", 1);
  tree->SetBranchStatus("tpcSignalN", 1);
  tree->SetBranchStatus("pidType", 1);
    
  
  tree->SetBranchAddress("pTPC", &pTPC);
  tree->SetBranchAddress("dEdx", &dEdx);
  //tree->SetBranchAddress("dEdxExpecttanThetaCentered", &dEdxExpected);
  tree->SetBranchAddress("tanTheta", &tanTheta);
  tree->SetBranchAddress("fMultiplicity", &fMultiplicity);
  tree->SetBranchAddress("tpcSignalN", &tpcSignalN);
  tree->SetBranchAddress("pidType", &pidType);
  
  const Int_t numMultBins = TMath::Ceil((20000. - 0.) / multStepSize);
  
  const Int_t numPbins = 24;
  Double_t pBinEdges[numPbins * 2] = {
    0.3, 0.3 + 0.005*widthFactor,
    0.35, 0.35 + 0.005*widthFactor,
    0.4, 0.4 + 0.005*widthFactor,
    0.45, 0.45 + 0.005*widthFactor,
    0.5, 0.5 + 0.005*widthFactor,
    0.55, 0.55 + 0.005*widthFactor,
    0.6, 0.6 + 0.005*widthFactor,
    0.65, 0.65 + 0.005*widthFactor,
    0.7, 0.7 + 0.005*widthFactor,
    0.8, 0.8 + 0.01*widthFactor,
    0.9, 0.9 + 0.01*widthFactor,
    1.0, 1.0 + 0.01*widthFactor,
    1.1, 1.1 + 0.01*widthFactor,
    1.2, 1.2 + 0.01*widthFactor,
    1.3, 1.3 + 0.01*widthFactor,
    1.4, 1.4 + 0.01*widthFactor,
    1.5, 1.5 + 0.01*widthFactor,
    1.6, 1.6 + 0.02*widthFactor,
    1.7, 1.7 + 0.02*widthFactor,
    1.8, 1.8 + 0.02*widthFactor,
    1.9, 1.9 + 0.02*widthFactor,
    2.0, 2.0 + 0.1*widthFactor,
    2.5, 2.5 + 0.1*widthFactor,
    3.0, 3.0 + 0.15*widthFactor };
    
  const Int_t numTanThetaBins = 10;
  
  Double_t tanThetaBinEdges[numTanThetaBins * 2] = {
    -1.0, -0.8,
    -0.8, -0.6,
    -0.6, -0.4,
    -0.4, -0.2,
    -0.2, 0.0,
    0.0, 0.2,
    0.2, 0.4,
    0.4, 0.6,
    0.6, 0.8,
    0.8, 1.0 };
    /*TODO tanThetaAbs
  Double_t tanThetaBinEdges[numTanThetaBins * 2] = {
    0.0, 0.1,
    0.1, 0.2,
    0.2, 0.3,
    0.3, 0.4,
    0.4, 0.5,
    0.5, 0.6,
    0.6, 0.7,
    0.7, 0.8,
    0.8, 0.9,
    0.9, 1.0 };*/
    
  TH2D* h2[numPbins][numTanThetaBins];
  TH2D* h2AllEta[numPbins];
  
  for (Int_t i = 0; i < numPbins; i++) {
    h2AllEta[i] = new TH2D(Form("h2AllEta_%d", i), Form("p %.2f - %.2f GeV/c, all #eta; Multiplicity; dEdx", pBinEdges[2*i], pBinEdges[2*i+1]),
                           numMultBins, 0, 20000, 590, 10.0, 600.0);
    
    for (Int_t j = 0; j < numTanThetaBins; j++) {
      h2[i][j] = new TH2D(Form("h2_%d_%d", i, j), Form("p %.2f - %.2f GeV/c & %.2f < tanTheta #leq %.2f; Multiplicity; dEdx",
      //TODO tanThetaAbs h2[i][j] = new TH2D(Form("h2_%d_%d", i, j), Form("p %.2f - %.2f GeV/c & %.2f < |tanTheta| #leq %.2f; Multiplicity; dEdx",
                                                       pBinEdges[2*i], pBinEdges[2*i+1], tanThetaBinEdges[2*j], tanThetaBinEdges[2*j+1]),
                          numMultBins, 0, 20000, 590, 10.0, 600.0);
    }
  }
  
  TSpline3* splProton = 0x0;
  
  if (correct) {
    printf("CORRECTING data for multiplicity dependence...TODO OBSOLETE - DO NOT USE!!!!\n");
    
    TFile* f = 0x0;
    f = TFile::Open(pathNameSplines.Data());
    if (!f)  {
      std::cout << "Failed to open spline file \"" << pathNameSplines.Data() << "\"!" << std::endl;
      return -1;
    }
    
    splProton = (TSpline3*)f->Get(splineName.Data());
    if (!splProton) {
      std::cout << "Failed to load splines: " << splineName.Data() << "!" << std::endl;
      return -1;
    }    
    
    if (splineName.Contains("PION")){
      printf("Pion mass assumed: %f...\n", massPion);
      mass = massPion;
    }
    else {
      printf("Proton mass assumed: %f...\n", massProton);
      mass = massProton;
    }
  }
  else
    printf("NOT CORRECTING data for multiplicity dependence...\n");

  progressbar(0.);
  
  TF1 cutFunc("cutFunc","50./TMath::Power(x,1.3)", 0.05, 6);
  
  TF1 corrFunc("corrFunc", "pol2", 0, 0.01);
  corrFunc.SetParameter(0, 5.5e-6);
  corrFunc.SetParameter(1, -0.00436);
  corrFunc.SetParameter(2, 0.103);
  
  for (Long64_t i = 0; i < nTreeEntries; i++) {
    tree->GetEntry(i);
    
    if (i % 100000 == 0)
      progressbar(100. * (((Double_t)i) / nTreeEntries));

    // Filter tracks according to shape recognition
    //if (tpcSignalN < 60 || (!isMC && dEdx < cutFunc.Eval(pTPC))) continue;
    if (dEdx <= 0 || tpcSignalN < 60) {
      continue;
    }
    
    if (!isMC) {
      if ((pTPC < 0.6 && pidType != kTPCid) ||               // No V0s/TOF below 0.6 due to bad statistics
      //if ((pTPC < 0.6 && pidType == kTPCandTOFid) ||         // No TOF below 0.6 GeV/c due to bad statistics
          //TODO NOW Maybe don't use this line since the tails get fewer V0's than the main peak, so the tails are overall reduced (pTPC >= 0.6 && pTPC < pThreholdV0cut && pidType != kTPCandTOFid) || // Do not mix V0's and TOF in this range
          (pTPC >= 0.6 && pTPC < 0.9 && pidType != kTPCandTOFid) ||
          (pTPC >= 0.9 && pidType != kV0idPlusTOFaccepted && pidType != kV0idNoTOF))      // Only V0's WITH TOF above pThreholdV0cut due to contamination
        continue;
    }
  
    // CORRECT multiplicity dependence -> WARNING: Better take into account the eta dependence for dEdxSplines for the determination of corrFactor!!!
    if (correct) {
      Double_t dEdxSplines = 50. * splProton->Eval(pTPC / mass);
      Double_t dEdxSplinesInv = 1. / dEdxSplines;
      Double_t relSlope = corrFunc.Eval(dEdxSplinesInv);
      Double_t corrFactor = relSlope * fMultiplicity;
      dEdx /= 1. + corrFactor;
    }
    
    
    
    //TODO tanThetaAbs Double_t absTanTheta = TMath::Abs(tanTheta);
    
    for (Int_t j = 0; j < numPbins; j++) {
      if (pTPC >= pBinEdges[2*j] && pTPC < pBinEdges[2*j+1]) {
        h2AllEta[j]->Fill(fMultiplicity, dEdx);
        
        for (Int_t k = 0; k < numTanThetaBins; k++) {
          if (tanTheta >= tanThetaBinEdges[2*k] && tanTheta < tanThetaBinEdges[2*k+1]) {
          //TODO tanThetaAbs if (absTanTheta >= tanThetaBinEdges[2*k] && absTanTheta < tanThetaBinEdges[2*k+1]) {
            h2[j][k]->Fill(fMultiplicity, dEdx);
            break;
          }
        }
        break;
      }
    }
  }

  progressbar(100.);
  
  printf("\n\n");

  
  const Double_t MCfitThreshold = 0.006;
  const Double_t heightFractionForRange = 0.1;
  TObjArray arr;
  
  { // Fitting for all eta
    TGraphErrors* gSlvsOffset = new TGraphErrors(numPbins);
    TGraphErrors* gSlvsInvOffset = new TGraphErrors(numPbins);
    
    TGraphErrors* gSigmaSlopevsdEdx = new TGraphErrors(numPbins);
    TGraphErrors* gSigmaSlopevsInvdEdx = new TGraphErrors(numPbins);
    
    fSave->mkdir("allTanTheta");
    
    for (Int_t i = 0; i < numPbins; i++) {
      gSlvsOffset->SetPoint(i, -10, 0);
      gSlvsOffset->SetPointError(i, 0, 0);
          
      gSlvsInvOffset->SetPoint(i, -10, 0);
      gSlvsInvOffset->SetPointError(i, 0, 0);
      
      gSigmaSlopevsdEdx->SetPoint(i, -10, 0);
      gSigmaSlopevsdEdx->SetPointError(i, 0, 0);
      
      gSigmaSlopevsInvdEdx->SetPoint(i, -10, 0);
      gSigmaSlopevsInvdEdx->SetPointError(i, 0, 0);
      
      if (h2AllEta[i]->GetEntries() < 10)
        continue;
      
      printf("Momentum %.2f - %.2f GeV/c, all eta:\n", pBinEdges[2*i], pBinEdges[2*i+1]);
      TCanvas* c = new TCanvas(Form("canv_pBin%d_allEta", i), Form("canv_pBin%d_allEta", i), 100,10,1380,800);
      FitSlicesY(h2AllEta[i], heightFractionForRange, 10, "QN", &arr);
      h2AllEta[i]->GetYaxis()->SetRange(h2AllEta[i]->FindFirstBinAbove(1, 2), h2AllEta[i]->FindLastBinAbove(1,2));
      h2AllEta[i]->Draw("colz");
      
      TH1D* hMean = (TH1D*)arr.At(1);//new TH1D((TH1D&)*arr.At(1));
      TH1D* hSigma = (TH1D*)arr.At(2);//new TH1D((TH1D&)*arr.At(2));
      hMean->SetName(Form("hMean_allEta_pBin%d", i));
      hSigma->SetName(Form("hSigma_allEta_pBin%d", i));
      hSigma->SetLineColor(kMagenta);
      
      TH1D* hSigmaRel = new TH1D(*hSigma);
      hSigmaRel->SetName(Form("hSigmaRel_allEta_pBin%d", i));
      hSigmaRel->Divide(hMean);
      hSigmaRel->SetLineColor(kMagenta + 2);
      
      TFitResultPtr fitRes = hMean->Fit("pol1", "s", "same", 0, 12000);
      TFitResultPtr fitResSigma = hSigmaRel->Fit("pol1", "s", "same", 0, 12000);
      
      hMean->DrawClone("same");
      
      if (((Int_t)fitRes) != 0) {
        printf("Fit failed!\n");
      }
      else {
        Double_t p0 = fitRes->GetParams()[0];
        Double_t p1 = fitRes->GetParams()[1];
        Double_t err0 = fitRes->GetErrors()[0];
        Double_t err1 = fitRes->GetErrors()[1];
        
        Double_t totError = 99999;
        if (p0 != 0 && p1 != 0) {
          totError = TMath::Abs(p1 / p0 * TMath::Sqrt(TMath::Power(err0 / p0, 2) + TMath::Power(err1 / p1, 2)));
        }
        printf("Result: p1 / p0 = %f / %f = %e\n\n", p1, p0, (p0 != 0 ? p1 / p0 : -99999));
        if (p0 > 0) {
          gSlvsOffset->SetPoint(i, p0, p1/p0);
          gSlvsOffset->SetPointError(i, err0, totError);
          
          gSlvsInvOffset->SetPoint(i, 1./p0, p1/p0);
          gSlvsInvOffset->SetPointError(i, TMath::Abs(err0/p0/p0), totError);
          
          // Only makes sense, if valid p0 available
          if (((Int_t)fitResSigma) != 0) {
            printf("Sigma fit failed!\n");
          }
          else {
            Double_t pSigma0 = fitResSigma->GetParams()[0];
            Double_t errSigma0 = fitResSigma->GetErrors()[0];
            
            Double_t pSigma1 = fitResSigma->GetParams()[1];
            Double_t errSigma1 = fitResSigma->GetErrors()[1];
            
            printf("Result sigma: pSigma1 / pSigma0 = %f / %f = %e\n\n", pSigma1, pSigma0, (pSigma0 != 0 ? pSigma1 / pSigma0 : -99999));
            
            Double_t totErrorSigma = 99999;
            if (pSigma0 != 0 && pSigma1 != 0) {
              totErrorSigma = TMath::Abs(pSigma1 / pSigma0 * TMath::Sqrt(TMath::Power(errSigma0 / pSigma0, 2) + TMath::Power(errSigma1 / pSigma1, 2)));
            }
            
            if (pSigma0 > 0) {          
              gSigmaSlopevsdEdx->SetPoint(i, p0, pSigma1/pSigma0);
              gSigmaSlopevsdEdx->SetPointError(i, err0, totErrorSigma);
                
              gSigmaSlopevsInvdEdx->SetPoint(i, 1./p0, pSigma1/pSigma0);
              gSigmaSlopevsInvdEdx->SetPointError(i, TMath::Abs(err0/p0/p0), totErrorSigma);
              
               // Set artificially large errors for points with positive slope vs. multiplicity
              // (negative slope makes physically no sense and is most likely due to uncertainties)
              if (pSigma1 < 0) {
                gSigmaSlopevsdEdx->SetPointError(i, err0, 1000);
                gSigmaSlopevsInvdEdx->SetPointError(i, TMath::Abs(err0/p0/p0), 1000);
              }
            }
          }
        }
      }
      
      fSave->cd("allTanTheta");
      c->Write();
      hSigma->Write();
      hSigmaRel->Write();
      arr.Clear();
    }
  
    gSlvsOffset->SetName("slopes_AllEta");
    gSlvsOffset->SetTitle("Rel. slope vs. offset (all #eta)");
    gSlvsOffset->SetMarkerColor(kBlack);
    gSlvsOffset->SetLineColor(kBlack);
    gSlvsOffset->SetMarkerStyle(20);
    gSlvsOffset->Draw("Ap");
    gSlvsOffset->GetHistogram()->GetXaxis()->SetTitle("dEdx (offset)");
    gSlvsOffset->GetHistogram()->GetYaxis()->SetTitle("Multiplicity slope/offset");
    gSlvsOffset->Fit("pol2", "W EX0", "same", 50, 300);
    TF1* fitFuncResult = dynamic_cast<TF1*>(gSlvsOffset->GetListOfFunctions()->At(0));
    if (fitFuncResult)
      fitFuncResult->SetLineColor(kBlack);
    
    gSlvsInvOffset->SetName("slopes2_AllEta");
    gSlvsInvOffset->SetTitle("Rel. slope vs. 1./offset (all #eta)");
    gSlvsInvOffset->SetMarkerColor(kBlack);
    gSlvsInvOffset->SetLineColor(kBlack);
    gSlvsInvOffset->SetMarkerStyle(20);
    gSlvsInvOffset->Draw("Ap");
    gSlvsInvOffset->GetHistogram()->GetXaxis()->SetTitle("1./dEdx (1./offset)");
    gSlvsInvOffset->GetHistogram()->GetYaxis()->SetTitle("Multiplicity slope/offset");
    
    TFitResultPtr fitRes = gSlvsInvOffset->Fit("pol2", "S 0 W EX0", "", 0, 0.02);//isMC ? MCfitThreshold : 0., 0.018); // Estimate pol2 params for subsequent fit
    //TODO NEW Introduced parameter [4]
    TF1 fitFunc("fitFunc", "[0] + [1]*TMath::Max([4], TMath::Min(x, [3])) + [2] * TMath::Power(TMath::Max([4], TMath::Min(x, [3])), 2)", 0., 0.2);
    //TF1 fitFunc("fitFunc", "[0] + [1]*TMath::Min(x, [3]) + [2] * TMath::Power(TMath::Min(x, [3]), 2)", 0., 0.1);

    
    //TF1* fitFuncResult2 = dynamic_cast<TF1*>(gSlvsInvOffset->GetListOfFunctions()->At(0));
    if ((Int_t)fitRes == 0)
      fitFunc.SetParameters(fitRes->GetParams());
    fitFunc.SetParameter(3, isMC ? 0.1 : 0.018);
    fitFunc.SetParLimits(3, isMC ? 0.014 : 0.01, 0.2);//TODO NEW
    // No lower threshold for data
    if (isMC) {
      fitFunc.SetParameter(4, MCfitThreshold);
      fitFunc.SetParLimits(4, 0., 0.01);
    }
    else
      fitFunc.FixParameter(4, 0.);
    
    fitFunc.SetLineColor(kBlack);
    
    gSlvsInvOffset->Fit(&fitFunc, "EX0", "same", /*TODO NEW isMC ? MCfitThreshold : */0., 0.2);
    
    
    
    gSigmaSlopevsdEdx->SetName("sigmaSlopes_AllEta");
    gSigmaSlopevsdEdx->SetTitle("Sigma rel slope vs. dEdx offset (all #eta)");
    gSigmaSlopevsdEdx->SetMarkerColor(kBlack);
    gSigmaSlopevsdEdx->SetLineColor(kBlack);
    gSigmaSlopevsdEdx->SetMarkerStyle(20);
    gSigmaSlopevsdEdx->Draw("Ap");
    gSigmaSlopevsdEdx->GetHistogram()->GetXaxis()->SetTitle("dEdx (offset)");
    gSigmaSlopevsdEdx->GetHistogram()->GetYaxis()->SetTitle("Multiplicity sigma rel. slope");
    
    gSigmaSlopevsInvdEdx->SetName("sigmaSlopes2_AllEta");
    gSigmaSlopevsInvdEdx->SetTitle("Sigma rel slope vs. 1./dEdx offset (all #eta)");
    gSigmaSlopevsInvdEdx->SetMarkerColor(kBlack);
    gSigmaSlopevsInvdEdx->SetLineColor(kBlack);
    gSigmaSlopevsInvdEdx->SetMarkerStyle(20);
    gSigmaSlopevsInvdEdx->Draw("Ap");
    gSigmaSlopevsInvdEdx->GetHistogram()->GetXaxis()->SetTitle("1./dEdx (1./offset)");
    gSigmaSlopevsInvdEdx->GetHistogram()->GetYaxis()->SetTitle("Multiplicity sigma rel. slope");
    
    TFitResultPtr fitRes2 = gSigmaSlopevsInvdEdx->Fit("pol2", "S 0 W EX0", "same", 0., 0.012); // Estimate pol2 params for subsequent fit
    TF1 fitFuncSigma("fitFuncSigma_AllEta", "TMath::Max(0, [0] + [1]*TMath::Min(x, [3]) + [2] * TMath::Power(TMath::Min(x, [3]), 2))", 0., 0.2);
    if ((Int_t)fitRes2 == 0)
      fitFuncSigma.SetParameters(fitRes2->GetParams());
    fitFuncSigma.SetParameter(3, 0.012);
    fitFuncSigma.SetParLimits(3, 0.01, 0.2);
    
    fitFuncSigma.SetLineColor(kBlack);
    
    gSigmaSlopevsInvdEdx->Fit(&fitFuncSigma, "EX0", "same", 0., 0.2);
    
    
    fSave->cd("allTanTheta");
    gSlvsOffset->Write();
    gSlvsInvOffset->Write();
    
    gSigmaSlopevsdEdx->Write();
    gSigmaSlopevsInvdEdx->Write();
  }
  
  
  TGraphErrors* gSlvsTanTheta[numPbins];
  for (Int_t j = 0; j < numPbins; j++) {
    gSlvsTanTheta[j] = new TGraphErrors(numTanThetaBins);
    gSlvsTanTheta[j]->SetName(Form("tanThetaSlopes_pBin%d", j));
    gSlvsTanTheta[j]->SetTitle(Form("Rel. slope vs. tanTheta (%.4f #leq pTPC/(GeV/c) < %.4f)", pBinEdges[2*j], pBinEdges[2*j+1]));
    gSlvsTanTheta[j]->SetMarkerColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));
    gSlvsTanTheta[j]->SetLineColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));
    gSlvsTanTheta[j]->SetMarkerStyle(20);
  }
  
  TGraphErrors* gSigmaSlopevsTanTheta[numPbins];
  for (Int_t j = 0; j < numPbins; j++) {
    gSigmaSlopevsTanTheta[j] = new TGraphErrors(numTanThetaBins);
    gSigmaSlopevsTanTheta[j]->SetName(Form("tanThetaSigmaSlopes_pBin%d", j));
    gSigmaSlopevsTanTheta[j]->SetTitle(Form("Rel. sigma slope vs. tanTheta (%.4f #leq pTPC/(GeV/c) < %.4f)", pBinEdges[2*j], pBinEdges[2*j+1]));
    gSigmaSlopevsTanTheta[j]->SetMarkerColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));
    gSigmaSlopevsTanTheta[j]->SetLineColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));
    gSigmaSlopevsTanTheta[j]->SetMarkerStyle(20);
  }
  
  printf("\n\n\n");
  
  // Fitting for eta bins
  for (Int_t j = 0; j < numTanThetaBins; j++) {
    printf("%.2f < tanTheta <= %.2f:\n", tanThetaBinEdges[2*j], tanThetaBinEdges[2*j+1]); //TODO absTanTheta
    
    TGraphErrors* gSlvsOffset = new TGraphErrors(numPbins);
    TGraphErrors* gSlvsInvOffset = new TGraphErrors(numPbins);
    
    TGraphErrors* gSigmaSlopevsdEdx = new TGraphErrors(numPbins);
    TGraphErrors* gSigmaSlopevsInvdEdx = new TGraphErrors(numPbins);
    
    fSave->mkdir(Form("tanThetaBin%d", j));
    
    for (Int_t i = 0; i < numPbins; i++) {
      gSlvsOffset->SetPoint(i, -10, 0);
      gSlvsOffset->SetPointError(i, 0, 0);
          
      gSlvsInvOffset->SetPoint(i, -10, 0);
      gSlvsInvOffset->SetPointError(i, 0, 0);
      
      gSlvsTanTheta[i]->SetPoint(j, -10, 0);
      gSlvsTanTheta[i]->SetPointError(j, 0, 0);
      
      gSigmaSlopevsTanTheta[i]->SetPoint(j, -10, 0);
      gSigmaSlopevsTanTheta[i]->SetPointError(j, 0, 0);
      
      gSigmaSlopevsdEdx->SetPoint(i, -10, 0);
      gSigmaSlopevsdEdx->SetPointError(i, 0, 0);
      
      gSigmaSlopevsInvdEdx->SetPoint(i, -10, 0);
      gSigmaSlopevsInvdEdx->SetPointError(i, 0, 0);
      
      if (h2[i][j]->GetEntries() < 10)
        continue;
      
      printf("Momentum %.2f - %.2f GeV/c:\n", pBinEdges[2*i], pBinEdges[2*i+1]);
      TCanvas* c = new TCanvas(Form("canv_pBin%d_tanThetaBin%d", i, j), Form("canv_pBin%d_tanThetaBin_%d", i, j), 100,10,1380,800);
      FitSlicesY(h2[i][j], heightFractionForRange, 10, "QN", &arr);
      h2[i][j]->GetYaxis()->SetRange(h2[i][j]->FindFirstBinAbove(1, 2), h2[i][j]->FindLastBinAbove(1,2));
      h2[i][j]->Draw("colz");
      
      TH1D* hMean = (TH1D*)arr.At(1);//new TH1D((TH1D&)*arr.At(1));
      TH1D* hSigma = (TH1D*)arr.At(2);//new TH1D((TH1D&)*arr.At(2));
      hMean->SetName(Form("hMean_pBin%d_tanThetaBin_%d", i, j));
      hSigma->SetName(Form("hSigma_pBin%d_tanThetaBin_%d", i, j));
      hSigma->SetLineColor(kMagenta);
      
      TH1D* hSigmaRel = new TH1D(*hSigma);
      hSigmaRel->SetName(Form("hSigmaRel_pBin%d_tanThetaBin_%d", i, j));
      hSigmaRel->Divide(hMean);
      hSigmaRel->SetLineColor(kMagenta + 2);
      TFitResultPtr fitRes = hMean->Fit("pol1", "s", "same", 0, 12000);
      TFitResultPtr fitResSigma = hSigmaRel->Fit("pol1", "s", "same", 0, 12000);
      
      hMean->DrawClone("same");
      
      if (((Int_t)fitRes) != 0) {
        printf("Fit failed!\n");
      }
      else {
        Double_t p0 = fitRes->GetParams()[0];
        Double_t p1 = fitRes->GetParams()[1];
        Double_t err0 = fitRes->GetErrors()[0];
        Double_t err1 = fitRes->GetErrors()[1];
        
        Double_t totError = 99999;
        if (p0 != 0 && p1 != 0) {
          totError = TMath::Abs(p1 / p0 * TMath::Sqrt(TMath::Power(err0 / p0, 2) + TMath::Power(err1 / p1, 2)));
        }
        printf("Result: p1 / p0 = %f / %f = %e\n\n", p1, p0, (p0 != 0 ? p1 / p0 : -99999));
        if (p0 > 0) {
          gSlvsOffset->SetPoint(i, p0, p1/p0);
          gSlvsOffset->SetPointError(i, err0, totError);
          
          gSlvsInvOffset->SetPoint(i, 1./p0, p1/p0);
          gSlvsInvOffset->SetPointError(i, TMath::Abs(err0/p0/p0), totError);
          
          
          gSlvsTanTheta[i]->SetPoint(j, (tanThetaBinEdges[2*j] + tanThetaBinEdges[2*j+1]) / 2., p1/p0);
          gSlvsTanTheta[i]->SetPointError(j, (tanThetaBinEdges[2*j+1] - tanThetaBinEdges[2*j]) / 2., totError); 
          
          // Only makes sense, if valid p0 available
          if (((Int_t)fitResSigma) != 0) {
            printf("Sigma fit failed!\n");
          }
          else {
            Double_t pSigma0 = fitResSigma->GetParams()[0];
            Double_t errSigma0 = fitResSigma->GetErrors()[0];
            
            Double_t pSigma1 = fitResSigma->GetParams()[1];
            Double_t errSigma1 = fitResSigma->GetErrors()[1];
            
            printf("Result sigma: pSigma1 / pSigma0 = %f / %f = %e\n\n", pSigma1, pSigma0, (pSigma0 != 0 ? pSigma1 / pSigma0 : -99999));
            
            Double_t totErrorSigma = 99999;
            if (pSigma0 != 0 && pSigma1 != 0) {
              totErrorSigma = TMath::Abs(pSigma1 / pSigma0 * TMath::Sqrt(TMath::Power(errSigma0 / pSigma0, 2) + TMath::Power(errSigma1 / pSigma1, 2)));
            }
            
            if (pSigma0 > 0) {          
              gSigmaSlopevsdEdx->SetPoint(i, p0, pSigma1/pSigma0);
              gSigmaSlopevsdEdx->SetPointError(i, err0, totErrorSigma);
                
              gSigmaSlopevsInvdEdx->SetPoint(i, 1./p0, pSigma1/pSigma0);
              gSigmaSlopevsInvdEdx->SetPointError(i, TMath::Abs(err0/p0/p0), totErrorSigma);
              
              gSigmaSlopevsTanTheta[i]->SetPoint(j, (tanThetaBinEdges[2*j] + tanThetaBinEdges[2*j+1]) / 2., pSigma1/pSigma0);
              gSigmaSlopevsTanTheta[i]->SetPointError(j, (tanThetaBinEdges[2*j+1] - tanThetaBinEdges[2*j]) / 2., totErrorSigma); 
              
              // Set artificially large errors for points with positive slope vs. multiplicity
              // (negative slope makes physically no sense and is most likely due to uncertainties)
              if (pSigma1 < 0) {
                gSigmaSlopevsdEdx->SetPointError(i, err0, 1000);
                gSigmaSlopevsInvdEdx->SetPointError(i, TMath::Abs(err0/p0/p0), 1000);
                gSigmaSlopevsTanTheta[i]->SetPointError(j, (tanThetaBinEdges[2*j+1] - tanThetaBinEdges[2*j]) / 2., 1000); 
              }
            }
          }
        }
      }
      
      fSave->cd(Form("tanThetaBin%d", j));
      c->Write();
      hSigma->Write();
      hSigmaRel->Write();
      arr.Clear();
    }
    
    
    gSlvsOffset->SetName(Form("slopes_tanTheta%d", j));
    gSlvsOffset->SetTitle(Form("Rel. slope vs. offset (%.2f < tanTheta #leq %.2f)", tanThetaBinEdges[2*j], tanThetaBinEdges[2*j+1]));//TODO tanThetaAbs
    gSlvsOffset->SetMarkerColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));//   (2 + ((j > 2) ? j + 1 : j));
    gSlvsOffset->SetLineColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));//(2 + ((j > 2) ? j + 1 : j));
    gSlvsOffset->SetMarkerStyle(20);
    gSlvsOffset->Draw("Ap");
    gSlvsOffset->GetHistogram()->GetXaxis()->SetTitle("dEdx (offset)");
    gSlvsOffset->GetHistogram()->GetYaxis()->SetTitle("Multiplicity slope/offset");
    gSlvsOffset->Fit("pol2", "W Ex0", "same", 10, 2000);
    TF1* fitFuncResult = dynamic_cast<TF1*>(gSlvsOffset->GetListOfFunctions()->At(0));
    if (fitFuncResult)
      fitFuncResult->SetLineColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));//(2 + ((j > 2) ? j + 1 : j));
    
    gSlvsInvOffset->SetName(Form("slopes2_tanTheta%d", j));
    gSlvsInvOffset->SetTitle(Form("Rel. slope vs. 1./offset (%.2f < tanTheta #leq %.2f)", tanThetaBinEdges[2*j], tanThetaBinEdges[2*j+1]));//TODO tanThetaAbs
    gSlvsInvOffset->SetMarkerColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));//(2 + ((j > 2) ? j + 1 : j));
    gSlvsInvOffset->SetLineColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));//(2 + ((j > 2) ? j + 1 : j));
    gSlvsInvOffset->SetMarkerStyle(20);
    gSlvsInvOffset->Draw("Ap");
    gSlvsInvOffset->GetHistogram()->GetXaxis()->SetTitle("1./dEdx (1./offset)");
    gSlvsInvOffset->GetHistogram()->GetYaxis()->SetTitle("Multiplicity slope/offset");
    
    TFitResultPtr fitRes = gSlvsInvOffset->Fit("pol2", "S 0 W EX0", "same", isMC ? MCfitThreshold : 0., 0.018); // Estimate pol2 params for subsequent fit
    //TODO NEW Introduced parameter [4]
    TF1 fitFunc(Form("fitFunc_tanTheta%d", j), "[0] + [1]*TMath::Max([4], TMath::Min(x, [3])) + [2] * TMath::Power(TMath::Max([4], TMath::Min(x, [3])), 2)", 0., 0.2);
    //TF1 fitFunc(Form("fitFunc_tanTheta%d", j), "[0] + [1]*TMath::Min(x, [3]) + [2] * TMath::Power(TMath::Min(x, [3]), 2)", 0., 0.1);
    if ((Int_t)fitRes == 0)
      fitFunc.SetParameters(fitRes->GetParams());
    fitFunc.SetParameter(3, isMC ? 0.1 : 0.018);
    fitFunc.SetParLimits(3, isMC ? 0.014 : 0.01, 0.2);//TODO NEW
    // No lower threshold for data
    if (isMC) {
      fitFunc.SetParameter(4, MCfitThreshold);
      fitFunc.SetParLimits(4, 0., 0.01);
    }
    else
      fitFunc.FixParameter(4, 0.);
    
    fitFunc.SetLineColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));//(2 + ((j > 2) ? j + 1 : j));
    
    gSlvsInvOffset->Fit(&fitFunc, "EX0", "same", /*TODO NEW isMC ? MCfitThreshold : */0., 0.2);
    
    
    
    gSigmaSlopevsdEdx->SetName(Form("sigmaSlopes_tanTheta%d", j));
    gSigmaSlopevsdEdx->SetTitle(Form("Sigma rel slope vs. dEdx offset (%.2f < tanTheta #leq %.2f)", tanThetaBinEdges[2*j], 
                                     tanThetaBinEdges[2*j+1]));//TODO tanThetaAbs
    gSigmaSlopevsdEdx->SetMarkerColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));
    gSigmaSlopevsdEdx->SetLineColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));
    gSigmaSlopevsdEdx->SetMarkerStyle(20);
    gSigmaSlopevsdEdx->Draw("Ap");
    gSigmaSlopevsdEdx->GetHistogram()->GetXaxis()->SetTitle("dEdx (offset)");
    gSigmaSlopevsdEdx->GetHistogram()->GetYaxis()->SetTitle("Multiplicity sigma rel. slope");
    
    gSigmaSlopevsInvdEdx->SetName(Form("sigmaSlopes2_tanTheta%d", j));
    gSigmaSlopevsInvdEdx->SetTitle(Form("Sigma rel slope vs. 1./dEdx offset (%.2f < tanTheta #leq %.2f)", tanThetaBinEdges[2*j],
                                        tanThetaBinEdges[2*j+1]));//TODO tanThetaAbs
    gSigmaSlopevsInvdEdx->SetMarkerColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));
    gSigmaSlopevsInvdEdx->SetLineColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));
    gSigmaSlopevsInvdEdx->SetMarkerStyle(20);
    gSigmaSlopevsInvdEdx->Draw("Ap");
    gSigmaSlopevsInvdEdx->GetHistogram()->GetXaxis()->SetTitle("1./dEdx (1./offset)");
    gSigmaSlopevsInvdEdx->GetHistogram()->GetYaxis()->SetTitle("Multiplicity sigma rel. slope");
    
    
    TFitResultPtr fitRes2 = gSigmaSlopevsInvdEdx->Fit("pol2", "S 0 W EX0", "same", 0., 0.012); // Estimate pol2 params for subsequent fit
    TF1 fitFuncSigma(Form("fitFuncSigma_tanTheta%d", j), "TMath::Max(0, [0] + [1]*TMath::Min(x, [3]) + [2] * TMath::Power(TMath::Min(x, [3]), 2))",
                     0., 0.2);
    if ((Int_t)fitRes2 == 0)
      fitFuncSigma.SetParameters(fitRes2->GetParams());
    fitFuncSigma.SetParameter(3, 0.012);
    fitFuncSigma.SetParLimits(3, 0.006, 0.2); //TODO was 3, 0.01, 0.2
    
    fitFuncSigma.SetLineColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));//(2 + ((j > 2) ? j + 1 : j));
    
    gSigmaSlopevsInvdEdx->Fit(&fitFuncSigma, "EX0", "same", 0., 0.2);
    
    
    fSave->cd(Form("tanThetaBin%d", j));
    gSlvsOffset->Write();
    gSlvsInvOffset->Write();
    
    gSigmaSlopevsdEdx->Write();
    gSigmaSlopevsInvdEdx->Write();
    
    printf("\n\n\n");
  }
  
  
  
  
  fSave->mkdir("tanThetaFits");
  fSave->cd("tanThetaFits");

  for (Int_t j = 0; j < numPbins; j++) {
    /*
    // Scale graphs, just that first bin sits at +-1
    Double_t scale = 1.0;
    for (Int_t k = 0; k < gSlvsTanTheta[j]->GetN(); k++) {
      if (k == 0) {
        scale = 1. / TMath::Max(1e-10, TMath::Abs(gSlvsTanTheta[j]->GetY()[k]));
      }
      
      gSlvsTanTheta[j]->SetPoint(k, gSlvsTanTheta[j]->GetX()[k], gSlvsTanTheta[j]->GetY()[k] * scale);
      gSlvsTanTheta[j]->SetPointError(k, gSlvsTanTheta[j]->GetEX()[k], gSlvsTanTheta[j]->GetEY()[k] * scale);
    }
    */
    
    // Shift graphs, just that reference bin sits at 0
    Double_t shift = 0.0;
    const Int_t refBin = 7;
    shift = gSlvsTanTheta[j]->GetY()[refBin];
    
    for (Int_t k = 0; k < gSlvsTanTheta[j]->GetN(); k++) {      
      gSlvsTanTheta[j]->SetPoint(k, gSlvsTanTheta[j]->GetX()[k], gSlvsTanTheta[j]->GetY()[k] - shift);
    }
    
    gSlvsTanTheta[j]->Draw("Ap");
    gSlvsTanTheta[j]->GetHistogram()->GetXaxis()->SetTitle("tan(#Theta)");//TODO tanThetaAbs
    gSlvsTanTheta[j]->GetHistogram()->GetYaxis()->SetTitle("Multiplicity slope/offset");
    
    gSlvsTanTheta[j]->Fit("pol2", "EX0", "same", tanThetaBinEdges[0], tanThetaBinEdges[numTanThetaBins * 2 - 1]);
    
    TF1* fitFuncResult = dynamic_cast<TF1*>(gSlvsTanTheta[j]->GetListOfFunctions()->At(0));
    if (fitFuncResult)
      fitFuncResult->SetLineColor((j == 7) ? kBlack : (2 + ((j > 2) ? j + 1 : j)));
    gSlvsTanTheta[j]->Write();
    
    gSigmaSlopevsTanTheta[j]->Draw("Ap");
    gSigmaSlopevsTanTheta[j]->GetHistogram()->GetXaxis()->SetTitle("tanTheta");//TODO tanThetaAbs
    gSigmaSlopevsTanTheta[j]->GetHistogram()->GetYaxis()->SetTitle("Multiplicity rel. sigma slope");
    gSigmaSlopevsTanTheta[j]->Write();
  }
  
  fSave->cd();
  if (correct)
    corrFunc.Write();
  
  
  TNamed* settings = new TNamed(Form("Settings: widthFactor %f, multStepSize %d, isMC %d, correct %d, pathNameSplines %s, splineName %s\n",
                                     widthFactor, multStepSize, isMC, correct, pathNameSplines.Data(), splineName.Data()), "");
  settings->Write();
  fSave->Close();
  
  
  return 0;
}
