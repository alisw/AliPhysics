#include "TH1.h"
#include "TH2F.h"
//#include "THnSparseDefinitions.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "./AddUpSystematicErrors.C"
#include <iostream>
using namespace std;

enum sysErrorAdditionType {kQuadratic = 0, kLinear = 1, kNosysErrorCalculation = 2};

TString getStringFromTObjStrArray(TObjArray *arr, Int_t position)
{
  if (position > arr->GetEntriesFast()-1)
    return "";
  
  TObjString* objStr = (TObjString*)(arr->At(position));
  if (!objStr)
    return "";
  
  return objStr->GetString();
}

Double_t CalculateJointSystematicError(Double_t sysErrorOriginal, Double_t sysErrorUE, Int_t sysErrorAddition);
Double_t CalculateToPiRatioSysError(Double_t yieldk,Double_t sysErrorSpecies,Double_t yieldPion,Double_t sysErrorPion);
  
void SubtractUnderlyingEvent(TString jetFilePattern, TString ueFilePattern, TString jetPtStepsString, TString centStepsString, TString modesInputString, TString outputFilePattern, Int_t sysErrorAddition = kNosysErrorCalculation, TString fileToRebin = "", Bool_t compare = kFALSE) { 
  
  TObjArray* jetPtBins = jetPtStepsString.Tokenize(";");
	Int_t numJetPtBins = jetPtBins->GetEntriesFast();
	Int_t* jetPt = new Int_t[numJetPtBins];
	for (Int_t i=0;i<numJetPtBins;i++) {
		jetPt[i] = getStringFromTObjStrArray(jetPtBins, i).Atoi();
	}
	
  TObjArray* centBins = centStepsString.Tokenize(";");
	Int_t numCentralities = centBins->GetEntriesFast();
	Int_t* centralities = new Int_t[numCentralities];
	for (Int_t i=0;i<numCentralities;i++) {
		centralities[i] = getStringFromTObjStrArray(centBins, i).Atoi();
	}
  
  TObjArray* modeArray = modesInputString.Tokenize(";");
	Int_t maxModes = modeArray->GetEntriesFast();
	TString* modeString = new TString[maxModes];
	for (Int_t i=0;i<maxModes;i++) {
    modeString[i] = getStringFromTObjStrArray(modeArray, i);
	}  
  
  Bool_t kDrawElectrons = kFALSE;
  Bool_t kdoNotUseUEErrorForFractions = kTRUE;
  
  const Int_t nOfYieldHistograms = 5;
  const Int_t nOfTotalYieldHistograms = 3;  
  
  const Int_t nOfEventHistograms = 4;
  const Int_t nOfJetsHistograms = 2;
  
  Bool_t efficiencyCorrected = jetFilePattern.Contains("EfficiencyCorrection");
  
  //Names
  TString yieldNames[nOfYieldHistograms] = { "hYieldElectrons", "hYieldMuons", "hYieldPions", "hYieldKaons", "hYieldProtons" };
  if (efficiencyCorrected) {
      for (Int_t i=0;i<nOfYieldHistograms;++i) {
          yieldNames[i] = yieldNames[i] + "_corrected";
      }
  }
  TString yieldMCNames[nOfYieldHistograms] = { "hMCYield_e", "hMCYield_mu", "hMCYield_pi", "hMCYield_K", "hMCYield_p" };
  TString toPiMCNames[nOfYieldHistograms] = { "", "", "", "hMCToPi_K", "hMCToPi_p" };
  TString yieldNamessysError[nOfYieldHistograms] = { "hYieldElectrons_sysError_corrected", "hYieldMuons_sysError_corrected", "hYieldPions_sysError_corrected", "hYieldKaons_sysError_corrected", "hYieldProtons_sysError_corrected" };
  
  const TString histNames[nOfYieldHistograms] = {"hFractionElectrons", "hFractionMuons", "hFractionPions", "hFractionKaons", "hFractionProtons" };
  const TString histNamesToPiRatios[nOfYieldHistograms] = {"hRatioToPiElectrons", "hRatioToPiMuons", "", "hRatioToPiKaons", "hRatioToPiProtons" }; 
  TString eventNames[nOfEventHistograms] = { "fhEventsProcessed", "fhEventsTriggerSel", "fhEventsTriggerSelVtxCut", "fhEventsProcessedNoPileUpRejection" };
  TString jetNames[nOfJetsHistograms] = { "fh2FFJetPtGen", "fh2FFJetPtRec" };
  const Int_t nOfPionHistogram = 2;
  
  const TString sysErrorYieldsGraphNames[nOfYieldHistograms] = {"systematicErrorYields_electron", "systematicErrorYields_muon", "systematicErrorYields_pion", "systematicErrorYields_kaon", "systematicErrorYields_proton"};
  const TString sysErrorGraphNames[nOfYieldHistograms] = {"systematicError_electron", "systematicError_muon", "systematicError_pion", "systematicError_kaon", "systematicError_proton"};
  const TString sysErrorToPiRatioGraphNames[nOfYieldHistograms] = {"systematicErrorToPiRatio_electron", "systematicErrorToPiRatio_muon", "systematicErrorToPiRatio_pion", "systematicErrorToPiRatio_kaon", "systematicErrorToPiRatio_proton"};
  
  //Yield Histograms
  TH1 *yieldHistogramsOriginal[nOfYieldHistograms] = {0x0};
  TH1 *yieldHistogramsOriginal_sysErrors[nOfYieldHistograms] = {0x0};
  
  TH1 *yieldHistogramsOriginalMC[nOfYieldHistograms] = {0x0};
  
  TH1 *yieldHistogramsUE[nOfYieldHistograms] = {0x0};
  TH1 *yieldHistogramsUE_sysErrors[nOfYieldHistograms] = {0x0};
  
  TH1 *yieldHistogramsUEMC[nOfYieldHistograms] = {0x0};
  
  TH1 *yieldHistogramsUEsubtracted[nOfYieldHistograms] = {0x0};
  
  TH1 *yieldHistogramsUEsubtractedMC[nOfYieldHistograms] = {0x0};
  TH1 *ratioToPiHistogramsUEsubtractedMC[nOfYieldHistograms] = {0x0};
  
  TH1 *totalYieldHistograms[nOfTotalYieldHistograms] = {0x0};
	
  TH1* yieldHistogramFractions[nOfYieldHistograms] = {0x0};
  TH1* ratioToPiHistograms[nOfYieldHistograms] = {0x0};
  
  //Systematic Errors 
  TGraphAsymmErrors* systematicErrorYieldsOriginal[nOfYieldHistograms] = {0x0};
  TGraphAsymmErrors* systematicErrorYieldsUE[nOfYieldHistograms] = {0x0};
  TGraphAsymmErrors* systematicErrorYieldsUEsubtracted[nOfYieldHistograms] = {0x0};
  TGraphAsymmErrors* systematicErrorUEsubtracted[nOfYieldHistograms] = {0x0};
  TGraphAsymmErrors* systematicErrorToPiRatioUEsubtracted[nOfYieldHistograms] = {0x0};

  for (Int_t i=0;i<nOfYieldHistograms;++i) {
    systematicErrorYieldsUEsubtracted[i] = 0x0;
    systematicErrorUEsubtracted[i] = 0x0;
    systematicErrorToPiRatioUEsubtracted[i] = 0x0;
  }
  
  Bool_t recalculateSysErrors[nOfYieldHistograms] = {0x0};
	
  //Event Histograms
  
  TH1D *eventHistogramsOriginal[nOfEventHistograms] = {0x0};
  
  TH2D *jetHistogramsOriginal[nOfJetsHistograms] = {0x0};
	
  TH1D *eventHistogramsUEsubtracted[nOfEventHistograms] = {0x0};
  
  TH2D *jetHistogramsUEsubtracted[nOfJetsHistograms] = {0x0};
  
  TH1::SetDefaultSumw2();
  
  for (Int_t iCent = 0; iCent < numCentralities/2; iCent++) {
    TString centralityString = Form("_centrality%d_%d", centralities[2*iCent], centralities[2*iCent + 1]);
    for (Int_t iMode = 0; iMode < maxModes; iMode++) {
      for (Int_t iJetPt = 0; iJetPt < numJetPtBins/2; iJetPt++) {
        TString jetPtString = Form("_jetPt%d.0_%d.0", jetPt[2*iJetPt], jetPt[2*iJetPt + 1]);
        for (Int_t k=0;k<nOfYieldHistograms;++k) {
          if (sysErrorAddition == kNosysErrorCalculation)
            recalculateSysErrors[k] = kFALSE;
          else
            recalculateSysErrors[k] = kTRUE;
        }
    
        //Set File names
        TString inputName = Form(jetFilePattern, modeString[iMode].Data(), centralityString.Data(), jetPtString.Data());
        TString ueName = Form(ueFilePattern, modeString[iMode].Data(), centralityString.Data(), jetPtString.Data());
        TString outputName = Form(outputFilePattern, modeString[iMode].Data(), centralityString.Data(), jetPtString.Data());
    
        //Open file with complete yield
        TFile *all = new TFile(inputName,"READ");
        if (!all || all->IsZombie()) {
          std::cout << "Input file " << inputName << " could not be opened" << std::endl;
          continue;
        }
        
        if (!efficiencyCorrected) {
          //Loading Event Yield Histograms
          for (Int_t k=0;k<nOfEventHistograms;++k) {
              all->GetObject(eventNames[k], eventHistogramsOriginal[k]);
              if (!eventHistogramsOriginal[k])
              std::cout << "Could not retrieve " << eventNames[k] << " of file " << inputName << std::endl;
          }
          
          //Loading Number of Jets Histograms
          for (Int_t k=0;k<nOfJetsHistograms;++k) {
              all->GetObject(jetNames[k], jetHistogramsOriginal[k]);
              if (!jetHistogramsOriginal[k])
              std::cout << "Could not retrieve " << jetNames[k] << " of file " << inputName << std::endl;
          }
        }
        //Loading Original Yield Histograms and Errors
        for (Int_t k=0;k<nOfYieldHistograms;++k) {
          all->GetObject(yieldNames[k], yieldHistogramsOriginal[k]);
          if (!yieldHistogramsOriginal[k])
            std::cout << "Could not retrieve original " << yieldNames[k] << " of file " << inputName << std::endl;
            if (efficiencyCorrected) {
              all->GetObject(yieldNamessysError[k], yieldHistogramsOriginal_sysErrors[k]);
            }
            else {
              all->GetObject(sysErrorYieldsGraphNames[k], systematicErrorYieldsOriginal[k]);
              if (!systematicErrorYieldsOriginal[k]) {
                cout << "Could not load original " << sysErrorYieldsGraphNames[k] << " of file " << inputName << endl;
                recalculateSysErrors[k] = kFALSE;
              } 
              all->GetObject(sysErrorGraphNames[k],systematicErrorUEsubtracted[k]);
              all->GetObject(sysErrorToPiRatioGraphNames[k], systematicErrorToPiRatioUEsubtracted[k]);
            }
          yieldHistogramsOriginalMC[k] = (TH1D*)all->FindObjectAny(yieldMCNames[k].Data());
        }    
        
        //Opening underlying event file
        TFile *ueFile = new TFile(ueName,"READ");
        if (!ueFile || ueFile->IsZombie()) {
          std::cout << "Underlying event File " << ueName << " could not be opened" << std::endl;     
          continue;
        }
        
        cout << "Loading UE yield" << endl;
        //Loading UE Yield Histograms and Errors
        for (Int_t k=0;k<nOfYieldHistograms;++k) {
          ueFile->GetObject(yieldNames[k], yieldHistogramsUE[k]);
          if (!yieldHistogramsUE[k])
            std::cout << "Could not retrieve underlying event " << yieldNames[k] << " of file " << ueName << std::endl;
          
          if (efficiencyCorrected) {
            ueFile->GetObject(yieldNamessysError[k], yieldHistogramsUE_sysErrors[k]);
          }
          else {
            ueFile->GetObject(sysErrorYieldsGraphNames[k], systematicErrorYieldsUE[k]);
            if (!systematicErrorYieldsUE[k]) {
              cout << "Could not load original " << sysErrorYieldsGraphNames[k] << " of file " << ueName << endl;
              recalculateSysErrors[k] = kFALSE;
            }      
          }
          yieldHistogramsUEMC[k] = (TH1D*)ueFile->FindObjectAny(yieldMCNames[k].Data());
        }
        
        if (!efficiencyCorrected) {
          //Clone Event Histograms
          for (Int_t k=0;k<nOfEventHistograms;++k) {
              eventHistogramsUEsubtracted[k] = dynamic_cast<TH1D*>(eventHistogramsOriginal[k]->Clone(eventNames[k]));
              if (!eventHistogramsUEsubtracted[k])
              std::cout << "Could not clone " << eventNames[k] << std::endl;
          }
          
          //Clone Number of Jets Histograms
          for (Int_t k=0;k<nOfJetsHistograms;++k) {
              jetHistogramsUEsubtracted[k] = dynamic_cast<TH2D*>(jetHistogramsOriginal[k]->Clone(jetNames[k]));
              if (!jetHistogramsUEsubtracted[k])
              std::cout << "Could not clone " << jetNames[k] << std::endl;
          }
        }
      
        //Clone Yield Histograms
        for (Int_t k=0;k<nOfYieldHistograms;++k) {
          yieldHistogramsUEsubtracted[k] = dynamic_cast<TH1*>(yieldHistogramsOriginal[k]->Clone(yieldNames[k]));
          if (!yieldHistogramsUEsubtracted[k])
            std::cout << "Could not clone " << yieldNames[k] << std::endl;
        }
        
        cout << "subtract underlying event" << endl;
        //Subtract Underlying event
        for (Int_t k=0;k<nOfYieldHistograms;++k) {
          for (Int_t l=1;l<=yieldHistogramsUEsubtracted[k]->GetNbinsX();++l) {
            Double_t originalYield = yieldHistogramsOriginal[k]->GetBinContent(l);
            Double_t ueYield = yieldHistogramsUE[k]->GetBinContent(l);
            if (originalYield >= 0.0 && ueYield >= 0.0)
              yieldHistogramsUEsubtracted[k]->SetBinContent(l, originalYield - ueYield);
            else if (originalYield >= 0.0)
              yieldHistogramsUEsubtracted[k]->SetBinContent(l, originalYield);
            else if (ueYield >= 0.0)
              yieldHistogramsUEsubtracted[k]->SetBinContent(l, -ueYield);
            else
              yieldHistogramsUEsubtracted[k]->SetBinContent(l, 0.0);
          }
          if (yieldHistogramsOriginalMC[k] && yieldHistogramsUEMC[k]) {
            yieldHistogramsUEsubtractedMC[k] = (TH1*)(yieldHistogramsOriginalMC[k]->Clone(yieldMCNames[k]));
            yieldHistogramsUEsubtractedMC[k]->Add(yieldHistogramsUEMC[k],-1);
          }
        }
  
        cout << "construct total yield" << endl;
        //Construct total Yield Histograms
        totalYieldHistograms[0] = dynamic_cast<TH1*>(yieldHistogramsOriginal[0]->Clone(yieldNames[0]));
        for (Int_t k=1;k<nOfYieldHistograms;++k) {
          totalYieldHistograms[0]->Add(yieldHistogramsOriginal[k],1.0);
        }
        totalYieldHistograms[0]->SetNameTitle("hTotalYield","Total Yield");
        
        //Construct total UEYields Histograms
        totalYieldHistograms[1] = dynamic_cast<TH1*>(yieldHistogramsUE[0]->Clone(yieldNames[0]));
        for (Int_t k=1;k<nOfYieldHistograms;++k) {
          totalYieldHistograms[1]->Add(yieldHistogramsUE[k],1.0);
        } 
        
        totalYieldHistograms[1]->SetNameTitle("hTotalYieldUE","Total Yield Underlying Event");
        
        //Construct total subtracted Yield Histograms
        totalYieldHistograms[2] = dynamic_cast<TH1*>(yieldHistogramsUEsubtracted[0]->Clone(yieldNames[0]));
        for (Int_t k=1;k<nOfYieldHistograms;++k) {
          totalYieldHistograms[2]->Add(yieldHistogramsUEsubtracted[k],1.0);
        } 
        
        totalYieldHistograms[2]->SetNameTitle("hTotalYieldUEsubtracted","Total Yield with subtracted UE");
        
        for (Int_t k=0;k<3;k++) {
          totalYieldHistograms[k]->SetLineColor(kBlack);
          totalYieldHistograms[k]->SetMarkerColor(kBlack);
          totalYieldHistograms[k]->SetMarkerStyle(kFullSquare);
        }
        
        cout << "construct fractions" << endl;
        for (Int_t iSpecies=0;iSpecies<nOfYieldHistograms;++iSpecies) {
          yieldHistogramFractions[iSpecies] = dynamic_cast<TH1*>(yieldHistogramsUEsubtracted[iSpecies]->Clone(histNames[iSpecies]));
          yieldHistogramFractions[iSpecies]->Divide(totalYieldHistograms[2]);
          if (iSpecies != nOfPionHistogram) {
            ratioToPiHistograms[iSpecies] = dynamic_cast<TH1*>(yieldHistogramsUEsubtracted[iSpecies]->Clone(histNamesToPiRatios[iSpecies]));
            ratioToPiHistograms[iSpecies]->Divide(yieldHistogramsUEsubtracted[nOfPionHistogram]);
          }
        }
        
        TCanvas* cMCYields = 0x0;
        TCanvas* cMCToPiRatios = 0x0;

        if (yieldHistogramsUEsubtractedMC[0]) {
          std::cout << "Construct MC Truth ratios and canvases" << std::endl;
        
          TH1D* yieldHistogramUEsubtractedMCTotal = (TH1D*)yieldHistogramsUEsubtractedMC[0]->Clone("hUEsubtractedMCTotalYield");
          for (Int_t k=1;k<nOfYieldHistograms;++k) {
            yieldHistogramUEsubtractedMCTotal->Add(yieldHistogramsUEsubtractedMC[k]);
          }
          yieldHistogramUEsubtractedMCTotal->SetLineColor(kBlack);
          yieldHistogramUEsubtractedMCTotal->SetMarkerColor(kBlack);
          yieldHistogramUEsubtractedMCTotal->SetMarkerStyle(kFullSquare);
          
          cMCYields = new TCanvas("cMCYields","MC truth yields");
          cMCYields->SetLogy(kTRUE);
          if (iMode == 0)
            cMCYields->SetLogx(kTRUE);  
          
          for (Int_t k=0;k<nOfYieldHistograms;k++) {
            yieldHistogramsUEsubtractedMC[k]->GetYaxis()->SetRangeUser(yieldHistogramsUEsubtractedMC[AliPID::kElectron]->GetBinContent(yieldHistogramsUEsubtractedMC[AliPID::kElectron]->FindLastBinAbove(0.)) * 0.1, yieldHistogramUEsubtractedMCTotal->GetBinContent(yieldHistogramUEsubtractedMCTotal->GetMaximumBin()) * 10.);
            yieldHistogramsUEsubtractedMC[k]->Draw(k==0 ? "E1" : "E1 same");
            if (yieldHistogramsUEsubtracted[k]) {
              yieldHistogramsUEsubtracted[k]->Draw("E1 same");
  //             hYieldCorrectedSysError[k]->Draw("E2 same");
            }
          }
          
          cMCToPiRatios = new TCanvas("cMCToPiRatios", "MC truth ratios");
          if (iMode == 0)
            cMCToPiRatios->SetLogx(kTRUE);
          
          for (Int_t k=3;k<5;k++) {
            ratioToPiHistogramsUEsubtractedMC[k] = (TH1*)yieldHistogramsUEsubtractedMC[k]->Clone(toPiMCNames[k].Data());
            ratioToPiHistogramsUEsubtractedMC[k]->SetNameTitle(toPiMCNames[k].Data(), toPiMCNames[k].Data());
            ratioToPiHistogramsUEsubtractedMC[k]->Divide(yieldHistogramsUEsubtractedMC[2]);
            ratioToPiHistogramsUEsubtractedMC[k]->GetYaxis()->SetRangeUser(0,1);
            ratioToPiHistogramsUEsubtractedMC[k]->Draw(k==0 ? "E1" : "E1 same");
            if (ratioToPiHistograms[k]) {
              ratioToPiHistograms[k]->Draw("E1 same");
  //             hRatioToPiCorrectedSysError[k]->Draw("E2 same");
            }
          }
        }
    
        cout << "Construct systematic errors for yields" << endl;
        const Int_t nBins = yieldHistogramsOriginal_sysErrors[0]->GetNbinsX();
        Double_t* x = 0x0;
        Double_t* y = 0x0;
        Double_t* xerr = 0x0;
        Double_t* yerrl = 0x0;
        Double_t* yerrh = 0x0;
        x = new Double_t[nBins+1];
        y = new Double_t[nBins+1];
        xerr = new Double_t[nBins+1];
        yerrl = new Double_t[nBins+1];
        yerrh = new Double_t[nBins+1];
        for (Int_t iSpecies=0;iSpecies<nOfYieldHistograms;++iSpecies) {
          for (Int_t bin=1;bin<=nBins;bin++) {
            x[bin] = yieldHistogramsOriginal_sysErrors[iSpecies]->GetBinCenter(bin);
            xerr[bin] = yieldHistogramsOriginal_sysErrors[iSpecies]->GetBinWidth(bin) / 2.;
            y[bin] = yieldHistogramsOriginal_sysErrors[iSpecies]->GetBinContent(bin);
            yerrl[bin] = yieldHistogramsOriginal_sysErrors[iSpecies]->GetBinErrorLow(bin);
            yerrh[bin] = yieldHistogramsOriginal_sysErrors[iSpecies]->GetBinErrorUp(bin);
          }
          systematicErrorYieldsOriginal[iSpecies] = new TGraphAsymmErrors(nBins, x, y, xerr, xerr, yerrl, yerrh);
          systematicErrorYieldsOriginal[iSpecies]->SetLineColor(yieldHistogramsOriginal_sysErrors[iSpecies]->GetMarkerColor());
          systematicErrorYieldsOriginal[iSpecies]->SetFillStyle(0);
          for (Int_t bin=1;bin<=nBins;bin++) {
            x[bin] = yieldHistogramsUE_sysErrors[iSpecies]->GetBinCenter(bin);
            xerr[bin] = yieldHistogramsUE_sysErrors[iSpecies]->GetBinWidth(bin) / 2.;
            y[bin] = yieldHistogramsUE_sysErrors[iSpecies]->GetBinContent(bin);
            yerrl[bin] = yieldHistogramsUE_sysErrors[iSpecies]->GetBinErrorLow(bin);
            yerrh[bin] = yieldHistogramsUE_sysErrors[iSpecies]->GetBinErrorUp(bin);
          }
          systematicErrorYieldsUE[iSpecies] = new TGraphAsymmErrors(nBins, x, y, xerr, xerr, yerrl, yerrh);
          systematicErrorYieldsUE[iSpecies]->SetLineColor(yieldHistogramsUE_sysErrors[iSpecies]->GetMarkerColor());
          systematicErrorYieldsUE[iSpecies]->SetFillStyle(0);          
        }
        delete[] x;
        delete[] y;
        delete[] xerr;
        delete[] yerrl;
        delete[] yerrh;
        x = y = xerr = yerrl = yerrh = 0x0;
        
        //Cloning 
        for (Int_t iSpecies=0;iSpecies<nOfYieldHistograms;++iSpecies) {
          systematicErrorYieldsUEsubtracted[iSpecies] = dynamic_cast<TGraphAsymmErrors*>(systematicErrorYieldsOriginal[iSpecies]->Clone());
          systematicErrorUEsubtracted[iSpecies] = dynamic_cast<TGraphAsymmErrors*>(systematicErrorYieldsOriginal[iSpecies]->Clone());
          systematicErrorToPiRatioUEsubtracted[iSpecies] = dynamic_cast<TGraphAsymmErrors*>(systematicErrorYieldsOriginal[iSpecies]->Clone());
        }
        
        //Calculate systematic Errors. First for pions, because their systematic errors are necessary for the calculation of the sys errors of the topi histograms
        cout << "Construct systematic errors for the pion " << endl;
        if (systematicErrorYieldsOriginal[nOfPionHistogram]->GetN() == systematicErrorYieldsUE[nOfPionHistogram]->GetN()) {
          for (Int_t n=0;n<systematicErrorYieldsUEsubtracted[nOfPionHistogram]->GetN();++n) {
            Double_t x = yieldHistogramsUEsubtracted[nOfPionHistogram]->GetBinCenter(n);
            Double_t y = yieldHistogramsUEsubtracted[nOfPionHistogram]->GetBinContent(n);
            Double_t ytotal = totalYieldHistograms[2]->GetBinContent(n);
            Double_t eyhigh_originalYield = systematicErrorYieldsOriginal[nOfPionHistogram]->GetErrorYhigh(n);
            Double_t eylow_originalYield = systematicErrorYieldsOriginal[nOfPionHistogram]->GetErrorYlow(n);
            Double_t eyhigh = CalculateJointSystematicError(eyhigh_originalYield,systematicErrorYieldsUE[nOfPionHistogram]->GetErrorYhigh(n),sysErrorAddition);
            Double_t eylow = CalculateJointSystematicError(eylow_originalYield,systematicErrorYieldsUE[nOfPionHistogram]->GetErrorYlow(n),sysErrorAddition);
            systematicErrorYieldsUEsubtracted[nOfPionHistogram]->SetPoint(n,x,y);
            systematicErrorYieldsUEsubtracted[nOfPionHistogram]->SetPointEYhigh(n,eyhigh);
            systematicErrorYieldsUEsubtracted[nOfPionHistogram]->SetPointEYlow(n,eylow);
            if (kdoNotUseUEErrorForFractions) {
              eyhigh = eyhigh_originalYield;
              eylow = eylow_originalYield;
            }
            if (ytotal > 0.0) {
              systematicErrorUEsubtracted[nOfPionHistogram]->SetPoint(n,x,y/ytotal);
              systematicErrorUEsubtracted[nOfPionHistogram]->SetPointEYhigh(n,eyhigh/ytotal);
              systematicErrorUEsubtracted[nOfPionHistogram]->SetPointEYlow(n,eylow/ytotal);
            }
            else {
              systematicErrorUEsubtracted[nOfPionHistogram]->SetPoint(n,x,0.0);
              systematicErrorUEsubtracted[nOfPionHistogram]->SetPointEYhigh(n,0.0);
              systematicErrorUEsubtracted[nOfPionHistogram]->SetPointEYlow(n,0.0);            
            }
          }
        } else {
          cout << "Don't recalculate systematic Errors, number of bins does not match" << endl;
        }
        
        cout << "Construct missing systematic errors" << endl;
        for (Int_t k=0;k<nOfYieldHistograms;++k) {
          if (k == nOfPionHistogram || (systematicErrorYieldsOriginal[k]->GetN() != systematicErrorYieldsUE[k]->GetN()))
            continue;
                    
  //         cout << "Construct error for " << k << endl;
          systematicErrorToPiRatioUEsubtracted[k]->SetTitle("(without respect to correlations between yield k and yield pion)");
          for (Int_t n=0;n<systematicErrorYieldsUEsubtracted[k]->GetN();++n) {
            Double_t eyhigh_pion, eylow_pion;
            if (kdoNotUseUEErrorForFractions) {
              eyhigh_pion = systematicErrorYieldsOriginal[nOfPionHistogram]->GetErrorYhigh(n);
              eylow_pion = systematicErrorYieldsOriginal[nOfPionHistogram]->GetErrorYlow(n);  
            }
            else {
              eyhigh_pion = systematicErrorYieldsUEsubtracted[nOfPionHistogram]->GetErrorYhigh(n);
              eylow_pion = systematicErrorYieldsUEsubtracted[nOfPionHistogram]->GetErrorYlow(n);       
            }
            Double_t x = yieldHistogramsUEsubtracted[k]->GetBinCenter(n);
            Double_t y = yieldHistogramsUEsubtracted[k]->GetBinContent(n);
            Double_t ytotal = totalYieldHistograms[2]->GetBinContent(n);
            Double_t eyhigh_originalYield = systematicErrorYieldsOriginal[k]->GetErrorYhigh(n);
            Double_t eylow_originalYield = systematicErrorYieldsOriginal[k]->GetErrorYlow(n);
            Double_t eyhigh = CalculateJointSystematicError(eyhigh_originalYield,systematicErrorYieldsUE[k]->GetErrorYhigh(n),sysErrorAddition);
            Double_t eylow = CalculateJointSystematicError(eylow_originalYield,systematicErrorYieldsUE[k]->GetErrorYlow(n),sysErrorAddition);
            systematicErrorYieldsUEsubtracted[k]->SetPoint(n,x,y);
            systematicErrorYieldsUEsubtracted[k]->SetPointEYhigh(n,eyhigh);
            systematicErrorYieldsUEsubtracted[k]->SetPointEYlow(n,eylow);
            if (kdoNotUseUEErrorForFractions) {
              eyhigh = eyhigh_originalYield;
              eylow = eylow_originalYield;
            }         
            if (ytotal > 0.0) {
              systematicErrorUEsubtracted[k]->SetPoint(n,x,y/ytotal);
              systematicErrorUEsubtracted[k]->SetPointEYhigh(n,eyhigh/ytotal);
              systematicErrorUEsubtracted[k]->SetPointEYlow(n,eylow/ytotal);
            }
            else {
              systematicErrorUEsubtracted[k]->SetPoint(n,x,0.0);
              systematicErrorUEsubtracted[k]->SetPointEYhigh(n,0.0);
              systematicErrorUEsubtracted[k]->SetPointEYlow(n,0.0);            
            }
            Double_t yPi = yieldHistogramsUEsubtracted[nOfPionHistogram]->GetBinContent(n);
            if (yPi > 0.0) {
              systematicErrorToPiRatioUEsubtracted[k]->SetPoint(n,x,y/yPi);
              systematicErrorToPiRatioUEsubtracted[k]->SetPointEYhigh(n,CalculateToPiRatioSysError(y,eyhigh,yPi,eyhigh_pion));
              systematicErrorToPiRatioUEsubtracted[k]->SetPointEYlow(n,CalculateToPiRatioSysError(y,eylow,yPi,eylow_pion));
            }
            else {
              systematicErrorToPiRatioUEsubtracted[k]->SetPoint(n,x,0.0);
              systematicErrorToPiRatioUEsubtracted[k]->SetPointEYhigh(n,CalculateToPiRatioSysError(y,eyhigh,0.0,eyhigh_pion));
              systematicErrorToPiRatioUEsubtracted[k]->SetPointEYlow(n,CalculateToPiRatioSysError(y,eylow,0.0,eylow_pion));            
            }
          }
        }

        const Double_t pLow = 0.15;
        const Double_t pHigh = 50.;
        
        TCanvas* cFractionsWithTotalSystematicError_Total = 0x0;
        TCanvas* cToPiRatiosWithTotalSystematicError_Total = 0x0; 
        TCanvas* cYieldsWithTotalSystematicError_Total = 0x0;
        TCanvas* cFractionsWithTotalSystematicError_UE = 0x0;
        TCanvas* cToPiRatiosWithTotalSystematicError_UE = 0x0; 
        TCanvas* cYieldsWithTotalSystematicError_UE = 0x0;
        TCanvas* cFractionsWithTotalSystematicError_UEsubtracted = 0x0;
        TCanvas* cToPiRatiosWithTotalSystematicError_UEsubtracted = 0x0; 
        TCanvas* cYieldsWithTotalSystematicError_UEsubtracted = 0x0;     
         
          cout << "Construct canvases" << endl;
          cout << systematicErrorUEsubtracted << endl << yieldHistogramFractions << endl;
//           Total
//           cFractionsWithTotalSystematicError_Total = DrawFractionHistos("cFractionsWithTotalSystematicError_Total", "Particle fraction", pLow,
//                                                                             pHigh, systematicErrorUEsubtracted, yieldHistogramFractions, kFraction, iMode, kDrawElectrons);
// 
//           cToPiRatiosWithTotalSystematicError_Total = DrawFractionHistos("cToPiRatiosWithTotalSystematicError_Total", "Ratio", pLow,
//                                                                             pHigh, systematicErrorToPiRatioUEsubtracted, ratioToPiHistograms , kToPiRatio, 
//                                                                             iMode, kDrawElectrons);
//           
//           cYieldsWithTotalSystematicError_Total = DrawFractionHistos("cYieldsWithTotalSystematicError_Total", "Yield", pLow, pHigh, systematicErrorYieldsOriginal, yieldHistogramsOriginal, kYield, iMode, kDrawElectrons);
//           
// //           UE
//           cFractionsWithTotalSystematicError_UE = DrawFractionHistos("cFractionsWithTotalSystematicError_UE", "Particle fraction", pLow,
//                                                                             pHigh, systematicErrorUEsubtracted, yieldHistogramFractions, kFraction, iMode, kDrawElectrons);
// 
//           cToPiRatiosWithTotalSystematicError_UE = DrawFractionHistos("cToPiRatiosWithTotalSystematicError_UE", "Ratio", pLow,
//                                                                             pHigh, systematicErrorToPiRatioUEsubtracted, ratioToPiHistograms , kToPiRatio, 
//                                                                             iMode, kDrawElectrons);
//           
//           cYieldsWithTotalSystematicError_UE = DrawFractionHistos("cYieldsWithTotalSystematicError_UE", "Yield", pLow, pHigh, systematicErrorYieldsUE, yieldHistogramsUE, kYield, iMode, kDrawElectrons); 
//           
        //UE subtracted
        cFractionsWithTotalSystematicError_UEsubtracted = DrawFractionHistos("cFractionsWithTotalSystematicError_UEsubtracted", "Particle fraction", pLow, pHigh, systematicErrorUEsubtracted, yieldHistogramFractions, kFraction, iMode, kDrawElectrons);

        cToPiRatiosWithTotalSystematicError_UEsubtracted = DrawFractionHistos("cToPiRatiosWithTotalSystematicError_UEsubtracted", "Ratio", pLow, pHigh, systematicErrorToPiRatioUEsubtracted, ratioToPiHistograms, kToPiRatio, iMode, kDrawElectrons);
        
        cYieldsWithTotalSystematicError_UEsubtracted = DrawFractionHistos("cYieldsWithTotalSystematicError_UEsubtracted", "Yield", pLow, pHigh, systematicErrorYieldsUEsubtracted, yieldHistogramsUEsubtracted, kYield, iMode, kDrawElectrons);
        totalYieldHistograms[2]->Draw("same E2");
        
        
        cout << "Canvases constructed" << endl;
        TH1* hTotalYieldRebinned = 0x0;
        TH1* hTotalYieldUERebinned = 0x0;
        TH1* hTotalYieldUEsubtractedRebinned = 0x0;
        
        TFile *rebinFile = 0x0;
        
        TH1F* ratioSignal = 0x0;
        TH1F* ratioUE = 0x0;
        
        if (!(fileToRebin == "")) {
          std::cout << "Rebinning." << std::endl;
          rebinFile = new TFile(fileToRebin,"READ");
          
          TString modus = iMode==0 ? "TrackPt" : "Z";
          TString zero = iJetPt==0 ? "0" : "";
          
          TString histogramToRebin = TString::Format("fh1FF%sRecCuts_%s%i_%i", modus.Data(), zero.Data(), jetPt[iJetPt], jetPt[iJetPt+1]);
          TString histogramToCompareSignal = TString::Format("fh1FF%sRecCuts_%s%i_%i", modus.Data(), zero.Data(), jetPt[iJetPt], jetPt[iJetPt+1]);
          TString histogramToCompareUE = TString::Format("fh1Bgr%sBckgPerpRecCuts_%s%i_%i", modus.Data(), zero.Data(), jetPt[iJetPt], jetPt[iJetPt+1]);
          
          if (rebinFile) {
            
            TDirectoryFile* df = 0x0;
            rebinFile->GetObject("clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00_noGenJets_trackTypeUndef_jetTypeUndef_cl0_trackRefs_BckgPerp",df);
            
            if (!df)
              break;
            
            TH1* rebinHistogram = 0x0;
            df->GetObject(histogramToRebin, rebinHistogram);
            
            if (!rebinHistogram)
              df->GetObject(TString::Format("fh1FF%sRecCuts_%s%i_%i", modus.Data(), "", 20, 30), rebinHistogram);

            if (rebinHistogram) {
            
              TAxis* axeForRebinning = rebinHistogram->GetXaxis();
              
              Int_t numberOfNewBins = axeForRebinning->GetXbins()->GetSize() - 1;
              const Double_t* xbinsRebinned = axeForRebinning->GetXbins()->GetArray();
              
              TH1* intermediate = 0x0;
              
              intermediate = dynamic_cast<TH1*>(totalYieldHistograms[0]->Clone("hintermediate"));
              for (Int_t a=1;a<=intermediate->GetNbinsX();++a) {
                intermediate->SetBinContent(a,intermediate->GetBinContent(a) * intermediate->GetBinWidth(a));
              }
              hTotalYieldRebinned = intermediate->Rebin(numberOfNewBins, "hTotalYieldRebinned", xbinsRebinned); 
              for (Int_t a=1;a<=intermediate->GetNbinsX();++a) {
                hTotalYieldRebinned->SetBinContent(a,hTotalYieldRebinned->GetBinContent(a) / hTotalYieldRebinned->GetBinWidth(a));
              }

              intermediate = dynamic_cast<TH1*>(totalYieldHistograms[1]->Clone("hintermediate"));
              for (Int_t a=1;a<=intermediate->GetNbinsX();++a) {
                intermediate->SetBinContent(a,intermediate->GetBinContent(a) * intermediate->GetBinWidth(a));
              }
              hTotalYieldUERebinned = intermediate->Rebin(numberOfNewBins, "hTotalYieldUERebinned", xbinsRebinned); 
              for (Int_t a=1;a<=intermediate->GetNbinsX();++a) {
                hTotalYieldUERebinned->SetBinContent(a,hTotalYieldUERebinned->GetBinContent(a) / hTotalYieldUERebinned->GetBinWidth(a));
              }  
              
              intermediate = dynamic_cast<TH1*>(totalYieldHistograms[2]->Clone("hintermediate"));
              for (Int_t a=1;a<=intermediate->GetNbinsX();++a) {
                intermediate->SetBinContent(a,intermediate->GetBinContent(a) * intermediate->GetBinWidth(a));
              }
              hTotalYieldUEsubtractedRebinned = intermediate->Rebin(numberOfNewBins, "hTotalYieldUEsubtractedRebinned", xbinsRebinned); 
              for (Int_t a=1;a<=intermediate->GetNbinsX();++a) {
                hTotalYieldUEsubtractedRebinned->SetBinContent(a,hTotalYieldUEsubtractedRebinned->GetBinContent(a) / hTotalYieldUEsubtractedRebinned->GetBinWidth(a));
              }
              
              if (compare) {
                TH1* compareHistogramSignal = 0x0;
                df->GetObject(histogramToCompareSignal, compareHistogramSignal);
                if (compareHistogramSignal) {
                  ratioSignal = new TH1F("h1ratioSignal", "Ratio of Signal results",numberOfNewBins,xbinsRebinned);
                  ratioSignal->Divide(hTotalYieldRebinned, compareHistogramSignal);
                }
              
                TH1* compareHistogramUE = 0x0;
                df->GetObject(histogramToCompareUE, compareHistogramUE);
                if (compareHistogramUE) {
                  ratioUE = new TH1F("h1ratioUE", "Ratio of UE results",numberOfNewBins,xbinsRebinned);
                  ratioUE->Divide(hTotalYieldUERebinned, compareHistogramUE);
                }
              }
            }
            else {
              std::cout << "Histogram " << histogramToRebin << " used for Rebinning not found." << std::endl;
            }
          }
          else {
            std::cout << "File " << fileToRebin << " containing the histogram for rebinning not found." << std::endl;
          }
          
          std::cout << "Finished Rebinning." << std::endl;
        }
        
        //Opening output file
        TFile *output = new TFile(outputName,"RECREATE");     //NOTE Change to update?
        if (output->IsZombie())
          std::cout << "Could not create output file " << outputName << std::endl;    

        cout << "Write Canvases" << endl;
        if (cFractionsWithTotalSystematicError_Total)
          cFractionsWithTotalSystematicError_Total->Write();
        
        if (cFractionsWithTotalSystematicError_UE)
        cFractionsWithTotalSystematicError_UE->Write();
        
        if (cFractionsWithTotalSystematicError_UEsubtracted)
        cFractionsWithTotalSystematicError_UEsubtracted->Write();
        
        if (cToPiRatiosWithTotalSystematicError_Total)
          cToPiRatiosWithTotalSystematicError_Total->Write();
        
        if (cToPiRatiosWithTotalSystematicError_UE)
        cToPiRatiosWithTotalSystematicError_UE->Write();
        
        if (cToPiRatiosWithTotalSystematicError_UEsubtracted)
        cToPiRatiosWithTotalSystematicError_UEsubtracted->Write();

        if (cYieldsWithTotalSystematicError_Total)
          cYieldsWithTotalSystematicError_Total->Write();
        
        if (cYieldsWithTotalSystematicError_UE)
        cYieldsWithTotalSystematicError_UE->Write();
        
        if (cYieldsWithTotalSystematicError_UEsubtracted)
        cYieldsWithTotalSystematicError_UEsubtracted->Write();      
        
        if (cMCYields) {
          output->mkdir("MC truth");
          output->cd("MC truth");
          
          cMCYields->Write();
          for (Int_t k=0;k<nOfYieldHistograms;++k)
            yieldHistogramsUEsubtractedMC[k]->Write();
          
          cMCToPiRatios->Write();
          for (Int_t k=3;k<nOfYieldHistograms;++k)
            yieldHistogramsUEsubtractedMC[k]->Write();
          
          output->cd();
        }
        
        cout << "Write histograms and systematic errors" << endl;
        //Write histograms
        for (Int_t k = 0;k<nOfYieldHistograms;++k) {
          if (systematicErrorUEsubtracted[k])
            systematicErrorUEsubtracted[k]->Write();
          
          yieldHistogramFractions[k]->Write();
          if (systematicErrorYieldsUEsubtracted[k])
            systematicErrorYieldsUEsubtracted[k]->Write();
          yieldHistogramsUEsubtracted[k]->Write();
          if (k != nOfPionHistogram) {
            if (systematicErrorToPiRatioUEsubtracted[k])
              systematicErrorToPiRatioUEsubtracted[k]->Write();
            
            ratioToPiHistograms[k]->Write();
          }
          
//           if (yieldHistogramsUEsubtractedMC[k])
//             yieldHistogramsUEsubtractedMC[k]->Write();
//           
//           if (ratioToPiHistogramsUEsubtractedMC[k])
//             ratioToPiHistogramsUEsubtractedMC[k]->Write();
        }
        
        cout << "Write total yield histograms" << endl;
        //Write total Yield Hisotgrams
        for (Int_t k=0;k<nOfTotalYieldHistograms;++k) {
          totalYieldHistograms[k]->Write();
        }
        
        if (!efficiencyCorrected) {
          cout << "Write event and jet histograms" << endl;
          //Write Event Histograms
          for (Int_t k=0;k<nOfEventHistograms;++k) {
            if (eventHistogramsOriginal[k])
              eventHistogramsUEsubtracted[k]->Write();
          }
          
          //Write Number of Jets Histograms
          for (Int_t k=0;k<nOfJetsHistograms;++k) {
            if (jetHistogramsUEsubtracted[k])
              jetHistogramsUEsubtracted[k]->Write();
          }   
        }
        
        //Write rebinned Histograms
        if (hTotalYieldRebinned) 
          hTotalYieldRebinned->Write();
        
        if (hTotalYieldUERebinned) 
          hTotalYieldUERebinned->Write();
        
        if (hTotalYieldUEsubtractedRebinned) 
          hTotalYieldUEsubtractedRebinned->Write();
        //Write ratio
        if (ratioSignal)
          ratioSignal->Write();
        
        if (ratioUE)
          ratioUE->Write();
        
        //Close files
        all->Close();
        ueFile->Close();
        output->Close();
        if (rebinFile)
          rebinFile->Close();
      }
    }
  }
  
  return;
}

      
Double_t CalculateJointSystematicError(Double_t sysErrorOriginal, Double_t sysErrorUE, Int_t sysErrorAddition) {
  if (sysErrorAddition == kQuadratic) {
    return TMath::Sqrt(sysErrorOriginal * sysErrorOriginal + sysErrorUE * sysErrorUE);
  }
  else if (sysErrorAddition == kLinear) {
    return sysErrorOriginal + sysErrorUE;
  }
  return sysErrorOriginal;
}

Double_t CalculateToPiRatioSysError(Double_t yieldk,Double_t sysErrorSpecies,Double_t yieldPion,Double_t sysErrorPion) {
  if (yieldPion > 0.0)
    return TMath::Sqrt(sysErrorSpecies * sysErrorSpecies/(yieldPion * yieldPion) + TMath::Power(yieldk * sysErrorPion/(yieldPion * yieldPion),2));
  else
    return 0.0;
}
