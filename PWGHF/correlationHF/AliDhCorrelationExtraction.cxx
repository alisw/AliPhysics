
/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

//
//  Base class to extract D-h correlation distribution from analysis task
//
//-----------------------------------------------------------------------
//  Author F.Colamaria
//  INFN Bari
//  fabio.colamaria@cern.ch
//-----------------------------------------------------------------------

#include "AliDhCorrelationExtraction.h"

//___________________________________________________________________________________________
AliDhCorrelationExtraction::AliDhCorrelationExtraction():
// default constructor
fFile(0x0),
fDirSE(0x0),
fDirME(0x0),
fMassList(0x0),
fTracksList(0x0),
fSECorrelationList(0x0),
fMECorrelationList(0x0),
fDmesonSpecies(kD0toKpi),
fSandBextraction(kBfromBinCount),
fSBscaling(kBinCountScaling),
fDmesonLabel(""),
fFileName(""), 
fDirNameSE(""), 
fListNameSE(""), 
fDirNameME(""), 
fListNameME(""),
fListMassName(""),
fMassHistoName(""),
fSECorrelHistoName(""),
fSECorrelHistoName_DstarBkg(""),
fMEsuffix(""),
fReadMassesOnly(kFALSE),
fRebinMassPlots(1),    
fNpTbins(3),
fFirstpTbin(1),
fLastpTbin(3),
fNumberOfSigmasFitter(3),
fCorrectPoolsSeparately(kTRUE),
fNpools(9),
fDmesonFitterSignal(0x0),
fDmesonFitterSignalError(0x0),
fDmesonFitterBackground(0x0),
fDmesonFitterBackgroundError(0x0),
fDMesonFitterSBCand(0x0),
fDMesonFitterSBCandErr(0x0),
fDmesonFitterMean(0x0),
fDmesonFitterMeanError(0x0),
fDmesonFitterSigma(0x0),
fDmesonFitterSigmaError(0x0),
fDmesonFitterSignificance(0x0),
fDmesonFitterSignificanceError(0x0),
fDmesonFitterSOverB(0x0),
fSignalSigmas(2),
fAutoSBRange(kFALSE),
fSBOuterSigmas(8),
fSBInnerSigmas(4),
fSBSingleBin(kTRUE),
fDeltaEtaMin(-1.0),
fDeltaEtaMax(1.0),
fSignalCorrel(0x0),
fBackgrCorrel(0x0),
fRangesSignL(0x0),
fRangesSignR(0x0),
fRangesSB1L(0x0),
fRangesSB1R(0x0),
fRangesSB2L(0x0),
fRangesSB2R(0x0),
fScaleFactor(0x0),
fIntegratePtBins(kFALSE),
fDebug(0),
fMassFit(0x0),
fBkgFit(0x0),
fMassHisto(0x0)
{

}

//___________________________________________________________________________________________
AliDhCorrelationExtraction::AliDhCorrelationExtraction(const AliDhCorrelationExtraction &source):
// copy constructor
fFile(source.fFile),
fDirSE(source.fDirSE),
fDirME(source.fDirME),
fMassList(source.fMassList),
fTracksList(source.fTracksList),
fSECorrelationList(source.fSECorrelationList),
fMECorrelationList(source.fMECorrelationList),
fDmesonSpecies(source.fDmesonSpecies),
fSandBextraction(source.fSandBextraction),
fSBscaling(source.fSBscaling),
fDmesonLabel(source.fDmesonLabel),
fFileName(source.fFileName), 
fDirNameSE(source.fDirNameSE), 
fListNameSE(source.fListNameSE), 
fDirNameME(source.fDirNameME), 
fListNameME(source.fListNameME),
fListMassName(source.fListMassName),
fMassHistoName(source.fMassHistoName),
fSECorrelHistoName(source.fSECorrelHistoName),
fSECorrelHistoName_DstarBkg(source.fSECorrelHistoName_DstarBkg),
fMEsuffix(source.fMEsuffix),
fReadMassesOnly(source.fReadMassesOnly),
fRebinMassPlots(source.fRebinMassPlots),    
fNpTbins(source.fNpTbins),
fFirstpTbin(source.fFirstpTbin),
fLastpTbin(source.fLastpTbin),
fNumberOfSigmasFitter(source.fNumberOfSigmasFitter),
fCorrectPoolsSeparately(source.fCorrectPoolsSeparately),
fNpools(source.fNpools),
fDmesonFitterSignal(source.fDmesonFitterSignal),
fDmesonFitterSignalError(source.fDmesonFitterSignalError),
fDmesonFitterBackground(source.fDmesonFitterBackground),
fDmesonFitterBackgroundError(source.fDmesonFitterBackgroundError),
fDMesonFitterSBCand(source.fDMesonFitterSBCand),
fDMesonFitterSBCandErr(source.fDMesonFitterSBCandErr),
fDmesonFitterMean(source.fDmesonFitterMean),
fDmesonFitterMeanError(source.fDmesonFitterMeanError),
fDmesonFitterSigma(source.fDmesonFitterSigma),
fDmesonFitterSigmaError(source.fDmesonFitterSigmaError),
fDmesonFitterSignificance(source.fDmesonFitterSignificance),
fDmesonFitterSignificanceError(source.fDmesonFitterSignificanceError),
fDmesonFitterSOverB(source.fDmesonFitterSOverB),
fSignalSigmas(source.fSignalSigmas),
fAutoSBRange(source.fAutoSBRange),
fSBOuterSigmas(source.fSBOuterSigmas),
fSBInnerSigmas(source.fSBInnerSigmas),
fSBSingleBin(source.fSBSingleBin),
fDeltaEtaMin(source.fDeltaEtaMin),
fDeltaEtaMax(source.fDeltaEtaMax),
fSignalCorrel(source.fSignalCorrel),
fBackgrCorrel(source.fBackgrCorrel),
fRangesSignL(source.fRangesSignL),
fRangesSignR(source.fRangesSignR),
fRangesSB1L(source.fRangesSB1L),
fRangesSB1R(source.fRangesSB1R),
fRangesSB2L(source.fRangesSB2L),
fRangesSB2R(source.fRangesSB2R),
fScaleFactor(source.fScaleFactor),
fIntegratePtBins(source.fIntegratePtBins),
fDebug(source.fDebug),
fMassFit(source.fMassFit),
fBkgFit(source.fBkgFit),
fMassHisto(source.fMassHisto)
{

}

//___________________________________________________________________________________________
AliDhCorrelationExtraction::~AliDhCorrelationExtraction() {
//destructor

}

//___________________________________________________________________________________________
Bool_t AliDhCorrelationExtraction::SetDmesonSpecie(DMesonSpecies k){

  if(k<0 || k>2) {
    printf("Error! D meson specie not correctly set!\n");
    return kFALSE;
  } else if(k==0) fDmesonLabel="Dzero";
  else if(k==1) fDmesonLabel="Dplus";  
  else fDmesonLabel="Dstar";

  fDmesonSpecies=k;
  return kTRUE;

}

//___________________________________________________________________________________________
Bool_t AliDhCorrelationExtraction::ReadInputs(){

  fFile = TFile::Open(fFileName.Data());
  if(!fFile->IsOpen()){
    std::cout << "File " << fFileName << " cannot be opened! check your file path!"; return kFALSE;
  }

  fDirSE = (TDirectoryFile*)fFile->Get(fDirNameSE.Data());
  if(!fReadMassesOnly ) fDirME = (TDirectoryFile*)fFile->Get(fDirNameME.Data());
    
  std::cout << "=================== " <<std::endl;
  std::cout << "Read inputs " <<std::endl;
  std::cout << "TDir SE    = " << fDirNameSE <<std::endl;
  std::cout << "TDir ME    = " << fDirNameME <<std::endl;
  std::cout << "TList SE   = " << fListNameSE <<std::endl;
  std::cout << "TList ME   = " << fListNameME <<std::endl;
  std::cout << "TList mass = " << fListMassName <<std::endl;
  std::cout << "=================== " <<std::endl;
  std::cout << " " <<std::endl;
    
  if(!fDirSE) {
    std::cout << "Cannot open the TDirectory for SE" << "\n" << "Your input: " << fDirNameSE << "\n" << "File content " << std::endl;
    fFile->ls();
    return kFALSE;
  }
    
  if(!fDirME && !fReadMassesOnly) {
    std::cout << "Cannot open the TDirectory for ME" << "\n" << "Your input: " << fDirNameME << "\n" << "File content " << std::endl;

    fFile->ls();
    return kFALSE;
  }
    
  fMassList = (TList*)fDirSE->Get(fListMassName);
  if(!fReadMassesOnly)   fSECorrelationList = (TList*)fDirSE->Get(fListNameSE);
  if(!fReadMassesOnly)    fMECorrelationList = (TList*)fDirME->Get(fListNameME);
    
  if(!fMassList) {
    std::cout << "Cannot open the TList for Mass" << "\n" << "Your input: " << fListMassName << "\n" << "TDirectory content " << std::endl;
    fDirSE->ls();
    return kFALSE;
  }
    
  if(!fSECorrelationList && !fReadMassesOnly) {
    std::cout << "Cannot open the TList for SE correlation" << "\n" << "Your input: " << fListNameSE << "\n" << "TDirectory content " << std::endl;
    fDirSE->ls();
    return kFALSE;
  }
    
  if(!fMECorrelationList && !fReadMassesOnly) {
    std::cout << "Cannot open the TList for ME correlation" << "\n" << "Your input: " << fListNameME << "\n" << "TDirectory content " << std::endl;
    fDirME->ls();
    return kFALSE;
  }

  return kTRUE;
       
}

//___________________________________________________________________________________________
Bool_t AliDhCorrelationExtraction::FitInvariantMass() {

  TCanvas *cDMass;

  Int_t nbinsDraw = fNpTbins;
  if(fIntegratePtBins) nbinsDraw = 1;

  if(nbinsDraw == 1) cDMass = new TCanvas("cDMass","Dmeson invariant mass in pt bins",0,0,1400,800);
  if(nbinsDraw == 2) {cDMass = new TCanvas("cDMass","Dmeson invariant mass in pt bins",0,0,1400,800);
    cDMass->Divide(2,1);
  }
  if(nbinsDraw == 3) {cDMass = new TCanvas("cDMass","Dmeson invariant mass in pt bins",0,0,1750,800);
    cDMass->Divide(3,1);
  }
  if(nbinsDraw == 4) {cDMass = new TCanvas("cDMass","Dmeson invariant mass in pt bins",0,0,1200,900);
    cDMass->Divide(2,2);
  }
  if(nbinsDraw == 5 || nbinsDraw == 6) {cDMass = new TCanvas("cDMass","Dmeson invariant mass in pt bins",0,0,1000,800);
    cDMass->Divide(3,2);
  }
  if(nbinsDraw == 7 || nbinsDraw == 8) {cDMass = new TCanvas("cDMass","Dmeson invariant mass in pt bins",0,0,1200,900);
    cDMass->Divide(4,2);
  }
  if(nbinsDraw == 9) {cDMass = new TCanvas("cDMass","Dmeson invariant mass in pt bins",0,0,1200,1000);
    cDMass->Divide(3,3);
  }
  if(nbinsDraw == 10 || nbinsDraw == 11 || nbinsDraw == 12) {cDMass = new TCanvas("cDMass","Dmeson invariant mass in pt bins",0,0,1200,1000);
    cDMass->Divide(3,4);
  }
  if(!cDMass){std::cout << "Cannot create canvas " <<std::endl; return kFALSE;}

  // define arrays containing the outputs of the d meson fits - the last member contains the output of the merged ptbin
  fDmesonFitterSignal = new Double_t[fNpTbins];
  fDmesonFitterSignalError = new Double_t[fNpTbins];
  fDmesonFitterBackground = new Double_t[fNpTbins];
  fDmesonFitterBackgroundError = new Double_t[fNpTbins];
  fDMesonFitterSBCand = new Double_t[fNpTbins];
  fDMesonFitterSBCandErr = new Double_t[fNpTbins];
  fDmesonFitterMean = new Double_t[fNpTbins];
  fDmesonFitterMeanError = new Double_t[fNpTbins];
  fDmesonFitterSigma = new Double_t[fNpTbins];
  fDmesonFitterSigmaError = new Double_t[fNpTbins];
  fDmesonFitterSignificance = new Double_t[fNpTbins];
  fDmesonFitterSignificanceError = new Double_t[fNpTbins];
  fDmesonFitterSOverB = new Double_t[fNpTbins];

  fSignalCorrel = new Double_t[fNpTbins]; 
  fBackgrCorrel = new Double_t[fNpTbins]; 
  fRangesSignL = new Double_t[fNpTbins];
  fRangesSignR = new Double_t[fNpTbins];
  fScaleFactor = new Double_t[fNpTbins];

  fMassFit = new TF1*[fNpTbins];
  fBkgFit = new TF1*[fNpTbins];

  fMassHisto = new TH1F*[fNpTbins];

  AliHFMassFitter *fitter = NULL;

  if(fLastpTbin-fFirstpTbin+1 != fNpTbins) { // if fits fails - display warning and keep going
      std::cout << "Error in the definition of the pT bins! Exiting..." << std::endl;
      return kFALSE;
  }

  Double_t LSBLowLim[fNpTbins], LSBUppLim[fNpTbins], RSBLowLim[fNpTbins], RSBUppLim[fNpTbins];
 
  //Loop on pT bins
  for(int i=0; i<fNpTbins; i++) {
    if(!fIntegratePtBins) { //Regular input extraction
      if(fDmesonSpecies==kDplusKpipi) {
        THnSparse *h = (THnSparse*)fMassList->FindObject(Form("%s%d",fMassHistoName.Data(),i+fFirstpTbin));
        fMassHisto[i] = (TH1F*)h->Projection(0);
      }
      else fMassHisto[i] = (TH1F*)fMassList->FindObject(Form("%s%d",fMassHistoName.Data(),i+fFirstpTbin));
    } else { //Integrate mass pT bins in wider correlation pT bin
      if(i>0) continue;
      MergeMassPlotsVsPt();
    }

    //Rebinning of the mass histogram - PAY ATTENTION! Now bins are different w.r.t. THnSparse 
    if(fRebinMassPlots>1) fMassHisto[i]->Rebin(fRebinMassPlots);

    //Settings of mass fitter and fit masses
    if(fDmesonSpecies==kD0toKpi) {
      fitter = new AliHFMassFitter(fMassHisto[i],fLeftFitRange,fRightFitRange,1,fBkgFitFunction,0);
      fitter->SetInitialGaussianMean(1.864);
      fitter->SetInitialGaussianSigma(0.010);
    } else if(fDmesonSpecies==kDplusKpipi) {
      fitter = new AliHFMassFitter(fMassHisto[i],fLeftFitRange,fRightFitRange,1,fBkgFitFunction,0);
      fitter->SetInitialGaussianMean(1.869);
      fitter->SetInitialGaussianSigma(0.010);
    } else if(fDmesonSpecies==kDStarD0pi) {
      fitter = new AliHFMassFitter(fMassHisto[i],fLeftFitRange,fRightFitRange,1,fBkgFitFunction,0);
      fitter->SetInitialGaussianMean(0.1454);
      fitter->SetInitialGaussianSigma(0.0005);
    } else {
      printf("Error! Decay channel not supported!\n"); 
      return kFALSE;
    }

    //Do the fit!
    std::cout << "Fitting mass plots... " <<std::endl;
    Bool_t isFitted = fitter->MassFitter(kFALSE);
        
    if(!isFitted) { // if fits fails - display warning and keep going
      std::cout << ">>>>>>>>>>> Fit at bin " << i << " not successful - skipping " << std::endl;
      return kFALSE;
    }

    //Get the outputs of the Fit
    fMassFit[i] = fitter->GetMassFunc();
    fMassFit[i]->SetRange(fLeftFitRange,fRightFitRange);
    fMassFit[i]->SetLineColor(4);

    fBkgFit[i] = fitter->GetBackgroundRecalcFunc();
    fBkgFit[i]->SetRange(fLeftFitRange,fRightFitRange);
    fBkgFit[i]->SetLineColor(2);

    fitter->Signal(fNumberOfSigmasFitter,fDmesonFitterSignal[i],fDmesonFitterSignalError[i]);
    fitter->Background(fNumberOfSigmasFitter,fDmesonFitterBackground[i],fDmesonFitterBackgroundError[i]);
    fitter->Significance(fNumberOfSigmasFitter,fDmesonFitterSignificance[i],fDmesonFitterSignificanceError[i]);
        
    if(fDmesonFitterBackground[i]) fDmesonFitterSOverB[i] = fDmesonFitterSignal[i]/fDmesonFitterBackground[i];
    else fDmesonFitterSOverB[i] = -1;
   
    fDmesonFitterMean[i] = fitter->GetMean();
    fDmesonFitterMeanError[i] = fitter->GetMeanUncertainty();
    fDmesonFitterSigma[i] = fitter->GetSigma();
    fDmesonFitterSigmaError[i]= fitter->GetSigmaUncertainty();

    //Draw mass plots
    cDMass->cd(i+1);
    SetHistoStyle(fMassHisto[i],1);
    fMassHisto[i]->GetXaxis()->SetRangeUser(fLeftFitRange,fRightFitRange);
    fMassHisto[i]->SetMinimum(0);
    fMassHisto[i]->SetMaximum(fMassHisto[i]->GetMaximum()*1.3);
    fMassHisto[i]->Sumw2();
    fMassHisto[i]->Draw();
    fBkgFit[i]->Draw("same");
    fMassFit[i]->Draw("same");
 
    TPaveText *paveTextDMass;
    DefinePaveText(paveTextDMass,0.57,0.5,0.88,0.88,Form("masspavetext_%d",i+fFirstpTbin));
       
    paveTextDMass->AddText(Form("Bin #%d",i+fFirstpTbin));
    paveTextDMass->AddText(Form("Signal (%.1f #sigma) =  %.1f #pm %.1f",fNumberOfSigmasFitter,fDmesonFitterSignal[i],fDmesonFitterSignalError[i]));
    paveTextDMass->AddText(Form("Bkg (%.1f #sigma) =  %.1f #pm %.1f",fNumberOfSigmasFitter,fDmesonFitterBackground[i],fDmesonFitterBackgroundError[i]));
    paveTextDMass->AddText(Form("Significance(%.1f #sigma) =  %.1f #pm %.1f",fNumberOfSigmasFitter,fDmesonFitterSignificance[i],fDmesonFitterSignificanceError[i]));
    paveTextDMass->AddText(Form("S/B(%.1f #sigma) =  %.1f ",fNumberOfSigmasFitter,fDmesonFitterSOverB[i]));
    paveTextDMass->AddText(Form("Mean = %.2f #pm %.2f (MeV/c^{2}) ",fDmesonFitterMean[i]*1000,fDmesonFitterMeanError[i]*1000));
    paveTextDMass->AddText(Form("Sigma = %.2f #pm %.2f (MeV/c^{2})  ",fDmesonFitterSigma[i]*1000,fDmesonFitterSigmaError[i]*1000));

    paveTextDMass->Draw("same");

    //Set Signal region and sideband region (if no external range is passed)
    // *** WARNING! Ranges are rounded so that they match bin edges of mass histograms!
    fRangesSignL[i] = fMassHisto[i]->GetXaxis()->GetBinLowEdge(fMassHisto[i]->FindBin(fDmesonFitterMean[i] - fSignalSigmas*fDmesonFitterSigma[i]));
    fRangesSignR[i] = fMassHisto[i]->GetXaxis()->GetBinUpEdge(fMassHisto[i]->FindBin(fDmesonFitterMean[i] + fSignalSigmas*fDmesonFitterSigma[i]));
    if(fAutoSBRange) {
      if (fDmesonSpecies!=kDStarD0pi) {
	LSBLowLim[i] = fMassHisto[i]->GetXaxis()->GetBinLowEdge(fMassHisto[i]->FindBin(fDmesonFitterMean[i] - fSBOuterSigmas*fDmesonFitterSigma[i]));
        LSBUppLim[i] = fMassHisto[i]->GetXaxis()->GetBinUpEdge(fMassHisto[i]->FindBin(fDmesonFitterMean[i] - fSBInnerSigmas*fDmesonFitterSigma[i]));
      }
      RSBLowLim[i] = fMassHisto[i]->GetXaxis()->GetBinLowEdge(fMassHisto[i]->FindBin(fDmesonFitterMean[i] + fSBInnerSigmas*fDmesonFitterSigma[i]));
      RSBUppLim[i] = fMassHisto[i]->GetXaxis()->GetBinUpEdge(fMassHisto[i]->FindBin(fDmesonFitterMean[i] + fSBOuterSigmas*fDmesonFitterSigma[i]));
      if(LSBLowLim[i]<fLeftFitRange || RSBUppLim[i]>fRightFitRange) {
        printf("Error! Sideband ranges go outside the fit range of the background fitting function!\nDepending on the approach, this could still be ok, but I'll exit anyway!\n");
	return kFALSE;
      }
    }

  }  //end of for cycle

  if(fAutoSBRange) {
    if (fDmesonSpecies==kDStarD0pi) SetSBRanges(RSBLowLim,RSBUppLim);
    else SetSBRanges(LSBLowLim,LSBUppLim,RSBLowLim,RSBUppLim);
  }

  for(int i=0; i<fNpTbins; i++) {
    //Extract signal and background in 'fSignalSigmas', for trigger normalization and SB correlation rescaling; then evaluate SB scaling factor
    if(fIntegratePtBins && i>0) continue;
    GetSignalAndBackgroundForNorm(i,fMassHisto[i]);
    GetSBScalingFactor(i,fMassHisto[i]);
  }

  std::cout << "Mass distributions for the bins considered were fitted" <<std::endl;

  cDMass->SaveAs(Form("Output_png/InvMassDistributions_%s_Bins%dto%d.png",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin));
  cDMass->SaveAs(Form("Output_Root/InvMassDistributions_%s_Bins%dto%d.root",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin));

   //HERE WE SHALL DELETE ALL DYNAMIC OBJECTS (ARRAYS, FUNCTIONS...)


  return kTRUE;
}

//___________________________________________________________________________________________
Bool_t AliDhCorrelationExtraction::ExtractCorrelations(Double_t thrMin, Double_t thrMax) {

  if(!fCorrectPoolsSeparately) fNpools = 1; //single histogram with integrated pools

  //Histograms definition
  TH2D* hSE_Sign[fNpools][fNpTbins];
  TH2D* hME_Sign[fNpools][fNpTbins];  
  TH2D* hSE_Sideb[fNpools][fNpTbins];
  TH2D* hME_Sideb[fNpools][fNpTbins];

  TH2D* hSE_Sign_PtInt[fNpools];
  TH2D* hME_Sign_PtInt[fNpools];
  TH2D* hSE_Sideb_PtInt[fNpools];
  TH2D* hME_Sideb_PtInt[fNpools];

  TH2D* hCorr_Sign_PtInt[fNpools];
  TH2D* hCorr_Sideb_PtInt[fNpools];

  TH2D* h2D_Sign;
  TH2D* h2D_Sideb;
  TH2D* h2D_Subtr;

  TH1D* h1D_Sign;
  TH1D* h1D_Sideb;
  TH1D* h1D_Subtr;
  TH1D* h1D_SubtrNorm;

  for(int iPool=0; iPool<fNpools; iPool++) {

    for(int iBin=0; iBin<fNpTbins; iBin++) {

      if(fIntegratePtBins && iBin>0) continue;

      //Retrieve 2D plots for SE and ME, signal and bkg regions, for each pTbin and pool
      hSE_Sign[iPool][iBin] = GetCorrelHisto(kSE,kSign,iPool,iBin,thrMin,thrMax);
      hME_Sign[iPool][iBin] = GetCorrelHisto(kME,kSign,iPool,iBin,thrMin,thrMax);
      hSE_Sideb[iPool][iBin] = GetCorrelHisto(kSE,kSideb,iPool,iBin,thrMin,thrMax);
      hME_Sideb[iPool][iBin] = GetCorrelHisto(kME,kSideb,iPool,iBin,thrMin,thrMax);

      //Scale bkg plots by ratio of signal region/sidebands
      hSE_Sideb[iPool][iBin]->Scale(fScaleFactor[iBin]);
      hME_Sideb[iPool][iBin]->Scale(fScaleFactor[iBin]);
      hSE_Sideb[iPool][iBin]->SetEntries(hSE_Sideb[iPool][iBin]->GetEntries()*fScaleFactor[iBin]);
      hME_Sideb[iPool][iBin]->SetEntries(hME_Sideb[iPool][iBin]->GetEntries()*fScaleFactor[iBin]);

      if(fDebug>=1) { 
	TCanvas *c = new TCanvas(Form("cInput_%d_%d_%1.1fto%1.1f",iPool,iBin,thrMin,thrMax),Form("InputCorr_%s_pool%d_bin%d_pTassoc%1.1fto%1.1f",fDmesonLabel.Data(),iPool,iBin,thrMin,thrMax),100,100,1200,900);
	c->Divide(2,2);
	c->cd(1);
	hSE_Sign[iPool][iBin]->SetMinimum(0);
	hSE_Sign[iPool][iBin]->Draw("lego2");
	c->cd(2);
	hME_Sign[iPool][iBin]->SetMinimum(0);
	hME_Sign[iPool][iBin]->Draw("lego2");
	c->cd(3);
	hSE_Sideb[iPool][iBin]->SetMinimum(0);
	hSE_Sideb[iPool][iBin]->Draw("lego2");
	c->cd(4);
	hME_Sideb[iPool][iBin]->SetMinimum(0);
	hME_Sideb[iPool][iBin]->Draw("lego2");
	c->SaveAs(Form("Output_png/InputCorr_%s_Canvas_pool%d_pTbin%d_thr%1.1fto%1.1f.png",fDmesonLabel.Data(),iPool,fFirstpTbin+iBin,thrMin,thrMax));
	c->SaveAs(Form("Output_Root/InputCorr_%s_Canvas_pool%d_pTbin%d_thr%1.1fto%1.1f.root",fDmesonLabel.Data(),iPool,fFirstpTbin+iBin,thrMin,thrMax));
      }

      if(iBin==0) {
        hSE_Sign_PtInt[iPool] = (TH2D*)hSE_Sign[iPool][iBin]->Clone(Form("hSE_Sign_PtInt_p%d",iPool));
        hME_Sign_PtInt[iPool] = (TH2D*)hME_Sign[iPool][iBin]->Clone(Form("hME_Sign_PtInt_p%d",iPool));
        hSE_Sideb_PtInt[iPool] = (TH2D*)hSE_Sideb[iPool][iBin]->Clone(Form("hSE_Sideb_PtInt_p%d",iPool));
        hME_Sideb_PtInt[iPool] = (TH2D*)hME_Sideb[iPool][iBin]->Clone(Form("hME_Sideb_PtInt_p%d",iPool));
      }
      else {
        hSE_Sign_PtInt[iPool]->Add(hSE_Sign[iPool][iBin]);
        hME_Sign_PtInt[iPool]->Add(hME_Sign[iPool][iBin]);
        hSE_Sideb_PtInt[iPool]->Add(hSE_Sideb[iPool][iBin]);
        hME_Sideb_PtInt[iPool]->Add(hME_Sideb[iPool][iBin]);
      }

    } // end of pT for

    NormalizeMEplot(hME_Sign_PtInt[iPool]);
    NormalizeMEplot(hME_Sideb_PtInt[iPool]);

    //Apply Event Mixing Correction
    hCorr_Sign_PtInt[iPool] = (TH2D*)hSE_Sign_PtInt[iPool]->Clone(Form("hCorr_Sign_PtInt_p%d",iPool));
    hCorr_Sign_PtInt[iPool]->Divide(hME_Sign_PtInt[iPool]);   
    hCorr_Sideb_PtInt[iPool] = (TH2D*)hSE_Sideb_PtInt[iPool]->Clone(Form("hCorr_Sideb_PtInt_p%d",iPool));
    hCorr_Sideb_PtInt[iPool]->Divide(hME_Sideb_PtInt[iPool]); 
    Double_t N_SEsign = 0, N_SEsideb = 0, N_sign = 0, N_sideb = 0;  
    for(int i=1;i<=hCorr_Sign_PtInt[iPool]->GetXaxis()->GetNbins();i++) {
      for(int j=1;j<=hCorr_Sign_PtInt[iPool]->GetYaxis()->GetNbins();j++) {
	N_SEsign += hSE_Sign_PtInt[iPool]->GetBinContent(i,j);
	N_SEsideb += hSE_Sideb_PtInt[iPool]->GetBinContent(i,j);
	N_sign += hCorr_Sign_PtInt[iPool]->GetBinContent(i,j);
	N_sideb += hCorr_Sideb_PtInt[iPool]->GetBinContent(i,j);
      }
    }
    hSE_Sign_PtInt[iPool]->SetEntries(N_SEsign); 
    hSE_Sideb_PtInt[iPool]->SetEntries(N_SEsideb); 
    hCorr_Sign_PtInt[iPool]->SetEntries(N_sign); 
    hCorr_Sideb_PtInt[iPool]->SetEntries(N_sideb); 

    if(fDebug>=1) {
      TCanvas *c = new TCanvas(Form("cSEME_%d_%1.1fto%1.1f",iPool,thrMin,thrMax),Form("cSEME_%s_pool%d_pTassoc%1.1fto%1.1f",fDmesonLabel.Data(),iPool,thrMin,thrMax),100,100,1600,900);
      c->Divide(3,2);
      c->cd(1);
      hSE_Sign_PtInt[iPool]->SetMinimum(0);
      hSE_Sign_PtInt[iPool]->Draw("lego2");
      c->cd(2);
      hME_Sign_PtInt[iPool]->SetMinimum(0);
      hME_Sign_PtInt[iPool]->Draw("lego2");
      c->cd(3);
      hCorr_Sign_PtInt[iPool]->SetMinimum(0);
      hCorr_Sign_PtInt[iPool]->Draw("lego2");
      c->cd(4);
      hSE_Sideb_PtInt[iPool]->SetMinimum(0);
      hSE_Sideb_PtInt[iPool]->Draw("lego2");
      c->cd(5);
      hME_Sideb_PtInt[iPool]->SetMinimum(0);
      hME_Sideb_PtInt[iPool]->Draw("lego2");
      c->cd(6);
      hCorr_Sideb_PtInt[iPool]->SetMinimum(0);
      hCorr_Sideb_PtInt[iPool]->Draw("lego2");
      c->SaveAs(Form("Output_png/CorrSEandME_%s_Canvas_PtIntBins%dto%d_pool%d_thr%1.1fto%1.1f.png",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin,iPool,thrMin,thrMax));
      c->SaveAs(Form("Output_Root/CorrSEandME_%s_Canvas_PtIntBins%dto%d_pool%d_thr%1.1fto%1.1f.root",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin,iPool,thrMin,thrMax));
    }

    //Pools integration
    if(iPool==0) {
      h2D_Sign = (TH2D*)hCorr_Sign_PtInt[0]->Clone("h2D_Sign");
      h2D_Sideb = (TH2D*)hCorr_Sideb_PtInt[0]->Clone("h2D_Sideb");
    } else {
      h2D_Sign->Add(hCorr_Sign_PtInt[iPool]);
      h2D_Sideb->Add(hCorr_Sideb_PtInt[iPool]);
    }

  } //end of pool for

  //Draw 2D plots (Signal region and Sidebands)
  TCanvas *c2D = new TCanvas(Form("c2D_IntPools_%1.1fto%1.1f",thrMin,thrMax),Form("c2D_%s_IntPools_pTassoc%1.1fto%1.1f",fDmesonLabel.Data(),thrMin,thrMax),100,100,1500,800);
  c2D->Divide(2,1);
  c2D->cd(1);
  h2D_Sign->SetMinimum(0);
  h2D_Sign->Draw("lego2");
  c2D->cd(2);
  h2D_Sideb->SetMinimum(0);
  h2D_Sideb->Draw("lego2");
  c2D->SaveAs(Form("Output_png/h2D_%s_Canvas_PtIntBins%dto%d_PoolInt_thr%1.1fto%1.1f.png",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin,thrMin,thrMax));
  c2D->SaveAs(Form("Output_Root/h2D_%s_Canvas_PtIntBins%dto%d_PoolInt_thr%1.1fto%1.1f.root",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin,thrMin,thrMax));
	
  //Bkg subtraction (2D plot)
  TCanvas *c2D_Sub = new TCanvas(Form("c2D_Subtr_IntPools_%1.1fto%1.1f",thrMin,thrMax),Form("c2D_%s_Subtr_IntPools_pTassoc%1.1fto%1.1f",fDmesonLabel.Data(),thrMin,thrMax),100,100,1500,800);
  h2D_Subtr = (TH2D*)h2D_Sign->Clone("h2D_Subtr");
  h2D_Subtr->Add(h2D_Sideb,-1);
  h2D_Subtr->SetEntries(h2D_Sign->GetEntries()-h2D_Sideb->GetEntries());
  h2D_Subtr->SetTitle("Signal region after sideb. subt. corr. - 2D");
  h2D_Subtr->Draw("lego2");
  c2D_Sub->SaveAs(Form("Output_png/h2D_%s_Subtr_Canvas_PtIntBins%dto%d_PoolInt_thr%1.1fto%1.1f.png",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin,thrMin,thrMax));
  c2D_Sub->SaveAs(Form("Output_Root/h2D_%s_Subtr_Canvas_PtIntBins%dto%d_PoolInt_thr%1.1fto%1.1f.root",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin,thrMin,thrMax));

  //Evaluate total number of triggers
  Double_t Snorm = 0, Bnorm = 0;
  for(int iBin=0; iBin<fNpTbins; iBin++) {
    if(fIntegratePtBins && iBin>0) continue;
    Snorm += fSignalCorrel[iBin];
    Bnorm += fBackgrCorrel[iBin];
  }

  //1D projection
  h1D_Sign = (TH1D*)h2D_Sign->ProjectionX("h1D_Sign");
  h1D_Sideb = (TH1D*)h2D_Sideb->ProjectionX("h1D_Sideb");
  h1D_Sign->SetTitle("Signal region correlations");
  h1D_Sideb->SetTitle("Sidebands correlations");
  h1D_Sign->Scale(1./h1D_Sign->GetXaxis()->GetBinWidth(1));
  h1D_Sideb->Scale(1./h1D_Sideb->GetXaxis()->GetBinWidth(1));

  //Bkg subtraction (1D plot)
  h1D_Subtr = (TH1D*)h1D_Sign->Clone("h1D_Subtr");
  h1D_Subtr->Add(h1D_Sideb,-1);
  h1D_Subtr->SetEntries(h1D_Sign->GetEntries()-h1D_Sideb->GetEntries());
  h1D_Subtr->SetTitle("Signal region after sideb. subt. corr.");

  //Draw 1D plots (Signal region, Sidebands, S-SB)
  TCanvas *c1D = new TCanvas(Form("c1D_IntPools_%1.1fto%1.1f",thrMin,thrMax),Form("c1D_%s_IntPools_pTassoc%1.1fto%1.1f",fDmesonLabel.Data(),thrMin,thrMax),100,100,1600,700);
  c1D->Divide(3,1);
  SetHistoCorrStyle(h1D_Sign);
  SetHistoCorrStyle(h1D_Sideb);
  SetHistoCorrStyle(h1D_Subtr);
  c1D->cd(1);
  h1D_Sign->Draw();
  c1D->cd(2);
  h1D_Sideb->Draw();
  c1D->cd(3);
  h1D_Subtr->Draw();
  c1D->SaveAs(Form("Output_png/h1D_%s_Canvas_PtIntBins%dto%d_PoolInt_thr%1.1fto%1.1f.png",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin,thrMin,thrMax));
  c1D->SaveAs(Form("Output_Root/h1D_%s_Canvas_PtIntBins%dto%d_PoolInt_thr%1.1fto%1.1f.root",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin,thrMin,thrMax));

  if(fDebug>=1) {
    //Draw 1D plots (Signal region, Sidebands, S-SB)
    TCanvas *c1D_Norm = new TCanvas(Form("c1D_IntPools_Normalized_%1.1fto%1.1f",thrMin,thrMax),Form("c1D_%s_IntPools_Normalized_pTassoc%1.1fto%1.1f",fDmesonLabel.Data(),thrMin,thrMax),100,100,1200,700);
    TH1D* h1D_SignN = (TH1D*)h1D_Sign->Clone("h1D_Sign_n");
    TH1D* h1D_SidebN = (TH1D*)h1D_Sideb->Clone("h1D_Sideb_n");
    TH1D* h1D_SubtrN = (TH1D*)h1D_Subtr->Clone("h1D_Subtr_n");
    h1D_SignN->Scale(1./(Snorm+Bnorm));
    h1D_SidebN->Scale(1./Bnorm);
    h1D_SubtrN->Scale(1./Snorm);
    SetHistoCorrStyle(h1D_SignN);
    SetHistoCorrStyle(h1D_SidebN);
    SetHistoCorrStyle(h1D_SubtrN);
    h1D_SignN->SetLineColor(kBlue);
    h1D_SignN->SetMarkerColor(kBlue);
    h1D_SidebN->SetLineColor(kRed);
    h1D_SidebN->SetMarkerColor(kRed);
    h1D_SidebN->SetTitle("Comparison of signal region, sidebands and signal correlations, normalized");
    h1D_SidebN->Draw();
    h1D_SignN->Draw("same");
    h1D_SubtrN->Draw("same");
    c1D_Norm->SaveAs(Form("Output_png/h1D_%s_Canvas_Normalized_PtIntBins%dto%d_PoolInt_thr%1.1fto%1.1f.png",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin,thrMin,thrMax));
    c1D_Norm->SaveAs(Form("Output_Root/h1D_%s_Canvas_Normalized_PtIntBins%dto%d_PoolInt_thr%1.1fto%1.1f.root",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin,thrMin,thrMax));
  }

  //Apply normalization to number of triggers
  h1D_SubtrNorm = (TH1D*)h1D_Subtr->Clone("h1D_SubtrNorm");
  h1D_SubtrNorm->Scale(1./Snorm);
  h1D_SubtrNorm->SetTitle("Signal region after sideb. subt. corr. - Normalized to # of triggers");  

  //Draw 1D plots (Signal region, Sidebands, S-SB)
  TCanvas *cFinal = new TCanvas(Form("cFinal_%1.1fto%1.1f",thrMin,thrMax),Form("cFinal_%s_IntPools_pTassoc%1.1fto%1.1f",fDmesonLabel.Data(),thrMin,thrMax),100,100,1200,700);
  SetHistoCorrStyle(h1D_SubtrNorm);
  h1D_SubtrNorm->Draw(); 
  cFinal->SaveAs(Form("Output_png/AzimCorrDistr_%s_Canvas_PtIntBins%dto%d_PoolInt_thr%1.1fto%1.1f.png",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin,thrMin,thrMax));
  cFinal->SaveAs(Form("Output_Root/AzimCorrDistr_%s_Canvas_PtIntBins%dto%d_PoolInt_thr%1.1fto%1.1f.root",fDmesonLabel.Data(),fFirstpTbin,fLastpTbin,thrMin,thrMax));

  return kTRUE;
}


//___________________________________________________________________________________________
void AliDhCorrelationExtraction::MergeMassPlotsVsPt() {

  TH1F* massPlotMerged;

  if(fDmesonSpecies==kDplusKpipi) {
    THnSparse *h = (THnSparse*)fMassList->FindObject(Form("%s%d",fMassHistoName.Data(),fFirstpTbin));
    massPlotMerged = (TH1F*)h->Projection(0);
  }
  else massPlotMerged = (TH1F*)fMassList->FindObject(Form("%s%d",fMassHistoName.Data(),fFirstpTbin));

  for(int i=1; i<fNpTbins; i++) {
    if(fDmesonSpecies==kDplusKpipi) {
      THnSparse *h = (THnSparse*)fMassList->FindObject(Form("%s%d",fMassHistoName.Data(),i+fFirstpTbin));
      massPlotMerged->Add((TH1F*)h->Projection(0));
    }
    else massPlotMerged->Add((TH1F*)fMassList->FindObject(Form("%s%d",fMassHistoName.Data(),i+fFirstpTbin)));
  }

  fMassHisto[0]=(TH1F*)massPlotMerged->Clone(Form("%s_Merged",massPlotMerged->GetName()));
  return;
}


//___________________________________________________________________________________________
void AliDhCorrelationExtraction::MergeCorrelPlotsVsPt(THnSparse* &hsparse, Int_t SEorME, Int_t SorSB, Int_t pool) {

  switch(fDmesonSpecies) {

    case (kD0toKpi): { //take 1st pT bin
      if(SEorME==kSE) {
        if(fCorrectPoolsSeparately) hsparse = (THnSparse*)(fSECorrelationList->FindObject(Form("%s%d_p%d",fSECorrelHistoName.Data(),fFirstpTbin,pool)))->Clone();
        else hsparse = (THnSparse*)(fSECorrelationList->FindObject(Form("%s%d",fSECorrelHistoName.Data(),fFirstpTbin)))->Clone();
      } else if(SEorME==kME) {
        if(fCorrectPoolsSeparately) hsparse = (THnSparse*)(fMECorrelationList->FindObject(Form("%s%d_p%d%s",fSECorrelHistoName.Data(),fFirstpTbin,pool,fMEsuffix.Data())))->Clone();
        else hsparse = (THnSparse*)(fMECorrelationList->FindObject(Form("%s%d%s",fSECorrelHistoName.Data(),fFirstpTbin,fMEsuffix.Data())))->Clone();
      }
      for(int iBin=1; iBin<fNpTbins; iBin++) { //now add other bins
        if(SEorME==kSE) {
          if(fCorrectPoolsSeparately) hsparse->Add((THnSparse*)fSECorrelationList->FindObject(Form("%s%d_p%d",fSECorrelHistoName.Data(),iBin+fFirstpTbin,pool)));
          else hsparse->Add((THnSparse*)fSECorrelationList->FindObject(Form("%s%d",fSECorrelHistoName.Data(),iBin+fFirstpTbin)));
        } else if(SEorME==kME) {
          if(fCorrectPoolsSeparately) hsparse->Add((THnSparse*)fMECorrelationList->FindObject(Form("%s%d_p%d%s",fSECorrelHistoName.Data(),iBin+fFirstpTbin,pool,fMEsuffix.Data())));
          else hsparse->Add((THnSparse*)fMECorrelationList->FindObject(Form("%s%d%s",fSECorrelHistoName.Data(),iBin+fFirstpTbin,fMEsuffix.Data())));
        }    
      }  
      break;
    } //end case D0

    case (kDplusKpipi): {  //take 1st pT bin 
      if(SEorME==kSE) hsparse = (THnSparse*)(fSECorrelationList->FindObject(Form("%s%d",fSECorrelHistoName.Data(),fFirstpTbin)))->Clone();
      else if(SEorME==kME) hsparse = (THnSparse*)(fMECorrelationList->FindObject(Form("%s%d%s",fSECorrelHistoName.Data(),fFirstpTbin,fMEsuffix.Data())))->Clone();
      for(int iBin=1; iBin<fNpTbins; iBin++) { //now add other bins
        if(SEorME==kSE) hsparse->Add((THnSparse*)fSECorrelationList->FindObject(Form("%s%d",fSECorrelHistoName.Data(),iBin+fFirstpTbin)));
        else if(SEorME==kME) hsparse->Add((THnSparse*)fMECorrelationList->FindObject(Form("%s%d%s",fSECorrelHistoName.Data(),iBin+fFirstpTbin,fMEsuffix.Data())));
      }
      break;
    } //end case D+

    case (kDStarD0pi): {  //take 1st pT bin 
      TString histoname = fSECorrelHistoName;
      if(SorSB==kSideb) histoname = fSECorrelHistoName_DstarBkg;
      if(SEorME==kSE) hsparse = (THnSparse*)(fSECorrelationList->FindObject(Form("%s%d",histoname.Data(),fFirstpTbin)))->Clone();
      else if(SEorME==kME) hsparse = (THnSparse*)(fMECorrelationList->FindObject(Form("%s%d%s",histoname.Data(),fFirstpTbin,fMEsuffix.Data())))->Clone();
      for(int iBin=1; iBin<fNpTbins; iBin++) { //now add other bins
        if(SEorME==kSE) hsparse->Add((THnSparse*)fSECorrelationList->FindObject(Form("%s%d",histoname.Data(),iBin+fFirstpTbin)));
        else if(SEorME==kME) hsparse->Add((THnSparse*)fMECorrelationList->FindObject(Form("%s%d%s",histoname.Data(),iBin+fFirstpTbin,fMEsuffix.Data())));
      }
      break;
    } //end case D*

    default:    
      printf("Error! Wrong setting in the D meson specie - Returning...\n");
      return; 
  } 

  return;
}


//___________________________________________________________________________________________
TH2D* AliDhCorrelationExtraction::GetCorrelHisto(Int_t SEorME, Int_t SorSB, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax) {

  switch(fDmesonSpecies) {

    case (kD0toKpi): {
      return GetCorrelHistoDzero(SEorME,SorSB,pool,pTbin,thrMin,thrMax);
      break;
    } //end case D0

    case (kDplusKpipi): {
      return GetCorrelHistoDplus(SEorME,SorSB,pool,pTbin,thrMin,thrMax);
      break;
    } //end case D+

    case (kDStarD0pi): {
      return GetCorrelHistoDstar(SEorME,SorSB,pool,pTbin,thrMin,thrMax);
      break;
    } //end case D*

    default:    
      printf("Error! Wrong setting in the D meson specie - Returning...\n");
      return 0x0; 
  } 

}


//___________________________________________________________________________________________
TH2D* AliDhCorrelationExtraction::GetCorrelHistoDzero(Int_t SEorME, Int_t SorSB, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax) {

  TH2D* h2D = new TH2D(); //pointer to be returned

  TH3D* h3D; //for projecting the TH2Sparse onto
  if(SEorME==kSE) {
    THnSparse *hsparse = 0x0;
    if(!fIntegratePtBins) {
      if(fCorrectPoolsSeparately) hsparse = (THnSparse*)fSECorrelationList->FindObject(Form("%s%d_p%d",fSECorrelHistoName.Data(),pTbin+fFirstpTbin,pool)); //D0: one histo x pool
      else hsparse = (THnSparse*)fSECorrelationList->FindObject(Form("%s%d",fSECorrelHistoName.Data(),pTbin+fFirstpTbin));
    } else MergeCorrelPlotsVsPt(hsparse,SEorME,SorSB,pool);
    Int_t ptBinTrMin = (Int_t)(hsparse->GetAxis(2)->FindBin(thrMin+0.01)); //the 0.01 to avoid bin edges!
    Int_t ptBinTrMax = (Int_t)(hsparse->GetAxis(2)->FindBin(thrMax-0.01));
    Int_t etaLowBin = (Int_t)(hsparse->GetAxis(4)->FindBin(fDeltaEtaMin+0.01));
    Int_t etaHighBin = (Int_t)(hsparse->GetAxis(4)->FindBin(fDeltaEtaMax-0.01));
    if(ptBinTrMax > hsparse->GetAxis(2)->GetNbins()) ptBinTrMax = hsparse->GetAxis(2)->GetNbins();
    if(etaHighBin > hsparse->GetAxis(4)->GetNbins()) etaHighBin = hsparse->GetAxis(4)->GetNbins();
    if(etaLowBin < 1) etaLowBin = 1;
    hsparse->GetAxis(2)->SetRange(ptBinTrMin,ptBinTrMax); //Apply cut on pT
    hsparse->GetAxis(3)->SetRange(1,hsparse->GetAxis(3)->GetNbins()); //AApply cut on displacement
    hsparse->GetAxis(4)->SetRange(etaLowBin,etaHighBin); //Apply cut on deltaEta
    if(fDebug>=2 && pool==0) printf("Bin ranges - pT: %d-%d, dEta: %d-%d, displ: %d-%d\n",ptBinTrMin,ptBinTrMax,etaLowBin,etaHighBin,1,hsparse->GetAxis(3)->GetNbins());
    h3D = (TH3D*)hsparse->Projection(1,4,0);//x,y,z axes
  } else if(SEorME==kME) {
    THnSparse *hsparse = 0x0;
    if(!fIntegratePtBins) {
      if(fCorrectPoolsSeparately) hsparse = (THnSparse*)fMECorrelationList->FindObject(Form("%s%d_p%d%s",fSECorrelHistoName.Data(),pTbin+fFirstpTbin,pool,fMEsuffix.Data()));
      else hsparse = (THnSparse*)fMECorrelationList->FindObject(Form("%s%d%s",fSECorrelHistoName.Data(),pTbin+fFirstpTbin,fMEsuffix.Data()));
    } else MergeCorrelPlotsVsPt(hsparse,SEorME,SorSB,pool);
    Int_t ptBinTrMin = (Int_t)(hsparse->GetAxis(3)->FindBin(thrMin+0.01));
    Int_t ptBinTrMax = (Int_t)(hsparse->GetAxis(3)->FindBin(thrMax-0.01));
    Int_t etaLowBin = (Int_t)(hsparse->GetAxis(2)->FindBin(fDeltaEtaMin+0.01));
    Int_t etaHighBin = (Int_t)(hsparse->GetAxis(2)->FindBin(fDeltaEtaMax-0.01));
    if(ptBinTrMax > hsparse->GetAxis(3)->GetNbins()) ptBinTrMax = hsparse->GetAxis(3)->GetNbins();
    if(etaHighBin > hsparse->GetAxis(2)->GetNbins()) etaHighBin = hsparse->GetAxis(4)->GetNbins();
    if(etaLowBin < 1) etaLowBin = 1;
    hsparse->GetAxis(2)->SetRange(etaLowBin,etaHighBin); //Apply cut on deltaEta
    hsparse->GetAxis(3)->SetRange(ptBinTrMin,ptBinTrMax); //Apply cut on pT
    if(fDebug>=2 && pool==0) printf("Bin ranges - pT: %d-%d, dEta: %d-%d\n",ptBinTrMin,ptBinTrMax,etaLowBin,etaHighBin);
    h3D = (TH3D*)hsparse->Projection(1,2,0);//x,y,z axes
  }
  //now restrict to signal or to sidebands
  if(SorSB==kSign) {
    h3D->GetXaxis()->SetRange(h3D->GetXaxis()->FindBin(fRangesSignL[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSignR[pTbin]-0.00001));
    if(fDebug>=2 && pool==0) printf("Signal range bins: %d-%d\n",h3D->GetXaxis()->FindBin(fRangesSignL[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSignR[pTbin]-0.00001));
    h2D = (TH2D*)h3D->Project3D("yz");
  } else if (SorSB==kSideb) {
    TH3D* h3Da = (TH3D*)h3D->Clone(Form("%s_sb2",h3D->GetName()));
    if(!fSBSingleBin) {
      h3Da->GetXaxis()->SetRange(h3D->GetXaxis()->FindBin(fRangesSB2L[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSB2R[pTbin]-0.00001));
      h3D->GetXaxis()->SetRange(h3D->GetXaxis()->FindBin(fRangesSB1L[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSB1R[pTbin]-0.00001));
      if(fDebug>=2 && pool==0) printf("SB1 range bins: %d-%d\n",h3D->GetXaxis()->FindBin(fRangesSB1L[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSB1R[pTbin]-0.00001));
      if(fDebug>=2 && pool==0) printf("SB2 range bins: %d-%d\n",h3D->GetXaxis()->FindBin(fRangesSB2L[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSB2R[pTbin]-0.00001));
    } else {
      h3Da->GetXaxis()->SetRange(h3Da->GetXaxis()->GetNbins(),h3Da->GetXaxis()->GetNbins());
      h3D->GetXaxis()->SetRange(1,1);
    }
    TH2D* h2Da = (TH2D*)h3Da->Project3D("yz"); 
    h2D = (TH2D*)h3D->Project3D("yz"); 
    h2D->Add(h2Da);
  }
  h3D->SetName(Form("%s_SE-ME%d_S-SB%d_3D",h3D->GetName(),SEorME,SorSB));      
  h2D->SetName(Form("%s_SE-ME%d_S-SB%d_2D",h2D->GetName(),SEorME,SorSB));
  return h2D;
}


//___________________________________________________________________________________________
TH2D* AliDhCorrelationExtraction::GetCorrelHistoDplus(Int_t SEorME, Int_t SorSB, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax) {

  TH2D* h2D = new TH2D(); //pointer to be returned

  TH3D* h3D; //for projecting the TH2Sparse onto
  THnSparse *hsparse = 0x0;
  if(!fIntegratePtBins) {
    if(SEorME==kSE) hsparse = (THnSparse*)fSECorrelationList->FindObject(Form("%s%d",fSECorrelHistoName.Data(),pTbin+fFirstpTbin));
    else if(SEorME==kME) hsparse = (THnSparse*)fMECorrelationList->FindObject(Form("%s%d%s",fSECorrelHistoName.Data(),pTbin+fFirstpTbin,fMEsuffix.Data()));
  } else MergeCorrelPlotsVsPt(hsparse,SEorME);
  Int_t ptBinTrMin = (Int_t)(hsparse->GetAxis(3)->FindBin(thrMin+0.01)); //the 0.01 to avoid bin edges!
  Int_t ptBinTrMax = (Int_t)(hsparse->GetAxis(3)->FindBin(thrMax-0.01));
  Int_t etaLowBin = (Int_t)(hsparse->GetAxis(2)->FindBin(fDeltaEtaMin+0.01));
  Int_t etaHighBin = (Int_t)(hsparse->GetAxis(2)->FindBin(fDeltaEtaMax-0.01));
  Int_t poolBin = pool+1; //D+: the pools are bins in a dedicated axis
  if(ptBinTrMax > hsparse->GetAxis(3)->GetNbins()) ptBinTrMax = hsparse->GetAxis(3)->GetNbins();
  if(etaHighBin > hsparse->GetAxis(2)->GetNbins()) etaHighBin = hsparse->GetAxis(2)->GetNbins();
  if(etaLowBin < 1) etaLowBin = 1;
  hsparse->GetAxis(2)->SetRange(etaLowBin,etaHighBin); //Apply cut on deltaEta
  hsparse->GetAxis(3)->SetRange(ptBinTrMin,ptBinTrMax); //Apply cut on pT
  if(fCorrectPoolsSeparately) hsparse->GetAxis(4)->SetRange(poolBin,poolBin); //Apply cut on pool (if the pool axis is there!)
  if(fDebug>=2 && pool==0) printf("Bin ranges - pT: %d-%d, dEta: %d-%d\n",ptBinTrMin,ptBinTrMax,etaLowBin,etaHighBin);
  h3D = (TH3D*)hsparse->Projection(0,2,1);//x,y,z axes
  //now restrict to signal or to sidebands
  if(SorSB==kSign) {
    h3D->GetXaxis()->SetRange(h3D->GetXaxis()->FindBin(fRangesSignL[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSignR[pTbin]-0.00001));
    if(fDebug>=2 && pool==0) printf("Signal range bins: %d-%d\n",h3D->GetXaxis()->FindBin(fRangesSignL[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSignR[pTbin]-0.00001));
    h2D = (TH2D*)h3D->Project3D("yz");
  } else if (SorSB==kSideb) {
    TH3D* h3Da = (TH3D*)h3D->Clone(Form("%s_sb2",h3D->GetName()));
    if(!fSBSingleBin) {
      h3Da->GetXaxis()->SetRange(h3D->GetXaxis()->FindBin(fRangesSB2L[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSB2R[pTbin]-0.00001));
      h3D->GetXaxis()->SetRange(h3D->GetXaxis()->FindBin(fRangesSB1L[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSB1R[pTbin]-0.00001));
      if(fDebug>=2 && pool==0) printf("SB1 range bins: %d-%d\n",h3D->GetXaxis()->FindBin(fRangesSB1L[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSB1R[pTbin]-0.00001));
      if(fDebug>=2 && pool==0) printf("SB2 range bins: %d-%d\n",h3D->GetXaxis()->FindBin(fRangesSB2L[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSB2R[pTbin]-0.00001));
    } else {
      h3Da->GetXaxis()->SetRange(h3Da->GetXaxis()->GetNbins(),h3Da->GetXaxis()->GetNbins());
      h3D->GetXaxis()->SetRange(1,1);
    }
    TH2D* h2Da = (TH2D*)h3Da->Project3D("yz"); 
    h2D = (TH2D*)h3D->Project3D("yz"); 
    h2D->Add(h2Da);
  }
  h3D->SetName(Form("%s_SE-ME%d_S-SB%d_3D",h3D->GetName(),SEorME,SorSB));      
  h2D->SetName(Form("%s_SE-ME%d_S-SB%d_2D",h2D->GetName(),SEorME,SorSB));
  return h2D;
}


//___________________________________________________________________________________________
TH2D* AliDhCorrelationExtraction::GetCorrelHistoDstar(Int_t SEorME, Int_t SorSB, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax) {

  TH2D* h2D = new TH2D(); //pointer to be returned

  TH3D* h3D; //for projecting the TH2Sparse onto
  THnSparse *hsparse = 0x0;
 
  TString histoname = fSECorrelHistoName;
  if(SorSB==kSideb) histoname = fSECorrelHistoName_DstarBkg;

  if(!fIntegratePtBins) {
    if(SEorME==kSE) hsparse = (THnSparse*)fSECorrelationList->FindObject(Form("%s%d",histoname.Data(),pTbin+fFirstpTbin));
    else if(SEorME==kME) hsparse = (THnSparse*)fMECorrelationList->FindObject(Form("%s%d%s",histoname.Data(),pTbin+fFirstpTbin,fMEsuffix.Data()));
  } else MergeCorrelPlotsVsPt(hsparse,SEorME,SorSB);
  Int_t ptBinTrMin = (Int_t)(hsparse->GetAxis(3)->FindBin(thrMin+0.01)); //the 0.01 to avoid bin edges!
  Int_t ptBinTrMax = (Int_t)(hsparse->GetAxis(3)->FindBin(thrMax-0.01));
  Int_t etaLowBin = (Int_t)(hsparse->GetAxis(2)->FindBin(fDeltaEtaMin+0.01));
  Int_t etaHighBin = (Int_t)(hsparse->GetAxis(2)->FindBin(fDeltaEtaMax-0.01));
  Int_t poolBin = pool+1; //D*: the pools are bins in a dedicated axis
  if(ptBinTrMax > hsparse->GetAxis(3)->GetNbins()) ptBinTrMax = hsparse->GetAxis(3)->GetNbins();
  if(etaHighBin > hsparse->GetAxis(2)->GetNbins()) etaHighBin = hsparse->GetAxis(2)->GetNbins();
  if(etaLowBin < 1) etaLowBin = 1;
  hsparse->GetAxis(2)->SetRange(etaLowBin,etaHighBin); //Apply cut on deltaEta
  hsparse->GetAxis(3)->SetRange(ptBinTrMin,ptBinTrMax); //Apply cut on pT
  if(fCorrectPoolsSeparately) hsparse->GetAxis(4)->SetRange(poolBin,poolBin); //Apply cut on pool (if the pool axis is there!)
  if(fDebug>=2 && pool==0) printf("Bin ranges - pT: %d-%d, dEta: %d-%d\n",ptBinTrMin,ptBinTrMax,etaLowBin,etaHighBin);
  h3D = (TH3D*)hsparse->Projection(1,2,0);//x,y,z axes
  //now restrict to signal or to sidebands
  if(SorSB==kSign) {
    h3D->GetXaxis()->SetRange(h3D->GetXaxis()->FindBin(fRangesSignL[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSignR[pTbin]-0.00001));
    if(fDebug>=2 && pool==0) printf("Signal range bins: %d-%d\n",h3D->GetXaxis()->FindBin(fRangesSignL[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSignR[pTbin]-0.00001));
    h2D = (TH2D*)h3D->Project3D("yz");
  } else if (SorSB==kSideb) { //only one sideband for the D* meson
    if(!fSBSingleBin) {
      h3D->GetXaxis()->SetRange(h3D->GetXaxis()->FindBin(fRangesSB1L[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSB1R[pTbin]-0.00001));
      if(fDebug>=2 && pool==0) printf("SB range bins: %d-%d\n",h3D->GetXaxis()->FindBin(fRangesSB1L[pTbin]+0.00001),h3D->GetXaxis()->FindBin(fRangesSB1R[pTbin]-0.00001));
    } else h3D->GetXaxis()->SetRange(h3D->GetXaxis()->GetNbins(),h3D->GetXaxis()->GetNbins()); //last bin for the only D* sideband (not used, but for the future...)
    h2D = (TH2D*)h3D->Project3D("yz");
  }
  h3D->SetName(Form("%s_SE-ME%d_S-SB%d_3D",h3D->GetName(),SEorME,SorSB));      
  h2D->SetName(Form("%s_SE-ME%d_S-SB%d_2D",h2D->GetName(),SEorME,SorSB));
  return h2D;
}


//___________________________________________________________________________________________
void AliDhCorrelationExtraction::SetSBRanges(Double_t* rangesSB1L, Double_t* rangesSB1R, Double_t* rangesSB2L, Double_t* rangesSB2R) {

if(!fAutoSBRange && !fSBSingleBin) {
  printf("*** WARNING! You are passing external SB ranges to the framework, and in the THnSparse the sidebands are mass-binned ***\n");
  printf("*** This is perfectly fine, provided that you match mass bins of THnSparse and of invariant mass plots! ***\n");
  if(fRebinMassPlots!=1) {
    printf("*** You are rebinning the mass plots! external SB ranges WILL bias the SB normalization factor! ***\n");
    printf("*** To remove the bias, you shall adapt the SB ranges to the rebinned edges of the mass histos (they must also be THnSparse mass edges) ***\n");
    getchar();
  }
}

  fRangesSB1L = new Double_t[fNpTbins];
  fRangesSB1R = new Double_t[fNpTbins];
  if(fDmesonSpecies!=kDStarD0pi) {
    fRangesSB2L = new Double_t[fNpTbins];
    fRangesSB2R = new Double_t[fNpTbins];
  }
  for(int i=0;i<fNpTbins;i++) {
    fRangesSB1L[i]=rangesSB1L[i];
    fRangesSB1R[i]=rangesSB1R[i];
    if(fDmesonSpecies!=kDStarD0pi) {
      fRangesSB2L[i]=rangesSB2L[i];
      fRangesSB2R[i]=rangesSB2R[i];
    }
    if(fDmesonSpecies==kDStarD0pi && fRangesSB1L[i]<0.148) {
      printf("\n*** WARNING! Sidebands exceed lower limit in the bkg THnSparse (0.148)! Results will be biased if committed task was used!\n");
      getchar();
    } 
  }
}

//___________________________________________________________________________________________
void AliDhCorrelationExtraction::GetSignalAndBackgroundForNorm(Int_t i, TH1F* &histo) {

  switch(fSandBextraction) {

    case (kSandBFromFit): {
      Double_t SandB = fMassFit[i]->Integral(fRangesSignL[i],fRangesSignR[i])/histo->GetBinWidth(histo->FindBin(fRangesSignL[i]));
      fBackgrCorrel[i] = fBkgFit[i]->Integral(fRangesSignL[i],fRangesSignR[i])/histo->GetBinWidth(histo->FindBin(fRangesSignL[i]));
      fSignalCorrel[i] = SandB - fBackgrCorrel[i];
      if(fDebug>=2) printf("S+B = %1.2f, S = %1.2f, B = %1.2f\n",SandB,fSignalCorrel[i],fBackgrCorrel[i]);
      break;
    }

    case (kSfromBinCount): {
      Double_t SandB = histo->Integral(histo->FindBin(fRangesSignL[i]+0.00001),histo->FindBin(fRangesSignR[i]-0.00001)); //the 0.00001 to move a bit from the edges
      fBackgrCorrel[i] = fBkgFit[i]->Integral(fRangesSignL[i],fRangesSignR[i])/histo->GetBinWidth(histo->FindBin(fRangesSignL[i]));
      fSignalCorrel[i] = SandB - fBackgrCorrel[i];
      if(fDebug>=2) printf("S+B = %1.2f, S = %1.2f, B = %1.2f\n",SandB,fSignalCorrel[i],fBackgrCorrel[i]);
      break;
    }

    case (kBfromBinCount): {
      Double_t SandB = histo->Integral(histo->FindBin(fRangesSignL[i]+0.00001),histo->FindBin(fRangesSignR[i]-0.00001)); //the 0.00001 to move a bit from the edges
      Double_t overallFitInt =  fMassFit[i]->Integral(fRangesSignL[i],fRangesSignR[i])/histo->GetBinWidth(histo->FindBin(fRangesSignL[i]));
      fSignalCorrel[i] = overallFitInt - fBkgFit[i]->Integral(fRangesSignL[i],fRangesSignR[i])/histo->GetBinWidth(histo->FindBin(fRangesSignL[i])); //S(fit) = SandB(fit)-B(fit)
      fBackgrCorrel[i] = SandB - fSignalCorrel[i];
      if(fDebug>=2) printf("S+B = %1.2f, S = %1.2f, B = %1.2f\n",SandB,fSignalCorrel[i],fBackgrCorrel[i]);
      break;
     }    

    default:    
      printf("Error! S and B extraction method wrongly set! Returning without evaluating S and B...\n");
      return; 
  } 

}

//___________________________________________________________________________________________
void AliDhCorrelationExtraction::GetSBScalingFactor(Int_t i, TH1F* &histo) {

  switch(fSBscaling) {

    case (kFitScaling): {
      Double_t binwidth = histo->GetBinWidth(histo->FindBin(fRangesSignL[i]));
      Double_t valueSB = 0;
      if(fSBSingleBin) { //in this case do not extend the fit range to bin limits
        if(fDmesonSpecies!=kDStarD0pi) valueSB = (fBkgFit[i]->Integral(fRangesSB1L[i],fRangesSB1R[i])+fBkgFit[i]->Integral(fRangesSB2L[i],fRangesSB2R[i]))/binwidth;
        else valueSB = fBkgFit[i]->Integral(fRangesSB1L[i],fRangesSB1R[i])/binwidth;
      } else { //in this case, extend the fit range to mass bin edges
        if(fDmesonSpecies!=kDStarD0pi) valueSB = (fBkgFit[i]->Integral(histo->GetXaxis()->GetBinLowEdge(histo->FindBin(fRangesSB1L[i]+0.00001)),histo->GetXaxis()->GetBinUpEdge(histo->FindBin(fRangesSB1R[i]-0.00001)))+fBkgFit[i]->Integral(histo->GetXaxis()->GetBinLowEdge(histo->FindBin(fRangesSB2L[i]+0.00001)),histo->GetXaxis()->GetBinUpEdge(histo->FindBin(fRangesSB2R[i]-0.00001))))/binwidth;
        else valueSB = fBkgFit[i]->Integral(histo->GetXaxis()->GetBinLowEdge(histo->FindBin(fRangesSB1L[i]+0.00001)),histo->GetXaxis()->GetBinUpEdge(histo->FindBin(fRangesSB1R[i]-0.00001)))/binwidth;   
      }  
      fScaleFactor[i] = fBackgrCorrel[i]/valueSB;
      if(fDebug>=2) printf("SB value = %1.2f, B = %1.2f, Scale Factor = %1.4f\n",valueSB,fBackgrCorrel[i],fScaleFactor[i]);
      break;
    }

    case (kBinCountScaling): {
      Double_t valueSB = 0;
      if(fDmesonSpecies!=kDStarD0pi) valueSB = histo->Integral(histo->FindBin(fRangesSB1L[i]+0.00001),histo->FindBin(fRangesSB1R[i]-0.00001)) + histo->Integral(histo->FindBin(fRangesSB2L[i]+0.00001),histo->FindBin(fRangesSB2R[i]-0.00001));
      else valueSB = histo->Integral(histo->FindBin(fRangesSB1L[i]+0.00001),histo->FindBin(fRangesSB1R[i]-0.00001));
      fScaleFactor[i] = fBackgrCorrel[i]/valueSB;   
      if(!fAutoSBRange && fSBSingleBin) RescaleSidebandsInMassBinEdges(i); //match SB ranges to inv.mass bin edges (if SB range passed from outside and single SB bin used)
      if(fDebug>=2) printf("SB value = %1.2f, B = %1.2f, Scale Factor (*renorm) = %1.4f\n",valueSB,fBackgrCorrel[i],fScaleFactor[i]);
      break;
    }

    default:    
      printf("Error! SB rescaling method wrongly set! Returning without evaluating S and B...\n");
      return; 
  } 

}


//___________________________________________________________________________________________
void AliDhCorrelationExtraction::NormalizeMEplot(TH2D* &histoME) {

  Double_t bin0phi = histoME->GetXaxis()->FindBin(0.);
  Double_t bin0eta = histoME->GetYaxis()->FindBin(0.);

  Double_t factorAdd = 0;
  for(int in=-1;in<=0;in++) factorAdd+=histoME->GetBinContent(bin0phi+in,bin0eta);
  for(int in=-1;in<=0;in++) factorAdd+=histoME->GetBinContent(bin0phi+in,bin0eta-1);
  factorAdd/=4.;
  histoME->Scale(1./factorAdd);

}


//___________________________________________________________________________________________
void AliDhCorrelationExtraction::PrintRanges() {
  if((fDmesonSpecies != kDStarD0pi && (!fRangesSignL || !fRangesSignR || !fRangesSB1L || !fRangesSB1R || !fRangesSB2L || !fRangesSB2R)) || 
     (fDmesonSpecies == kDStarD0pi && (!fRangesSignL || !fRangesSignR || !fRangesSB1L || !fRangesSB1R))) {
    printf("Error! Ranges not correctly set! Exiting...\n");
    return;
  }
  printf("******************************************\n");
  printf("  RANGES OF SIGNAL REGION AND SIDEBANDS\n\n");  
  printf(" NOTE: limits rounded to mass bin edges!\n");
  printf(" Signal region sigmas = %1.1f\n",fSignalSigmas);
  printf(" SB sigmas = %1.1f-%1.1f (ext. range: %d)\n",fSBInnerSigmas,fSBOuterSigmas,!fAutoSBRange);
  printf("******************************************\n");

  for(int i=0; i<fNpTbins; i++) {
    if(fIntegratePtBins && i>0) continue;
    printf("****************** Bin %d *****************\n",i+fFirstpTbin);
    printf("  Signal region = %1.4f - %1.4f (bins %d-%d)\n",fRangesSignL[i],fRangesSignR[i],fMassHisto[i]->FindBin(fDmesonFitterMean[i] - fSignalSigmas*fDmesonFitterSigma[i]),fMassHisto[i]->FindBin(fDmesonFitterMean[i] + fSignalSigmas*fDmesonFitterSigma[i]));
    printf("  SB1 region: = %1.4f - %1.4f\n",fRangesSB1L[i],fRangesSB1R[i]);
    if(fDmesonSpecies!=kDStarD0pi) printf("  SB2 region: = %1.4f - %1.4f\n",fRangesSB2L[i],fRangesSB2R[i]);
  }
  printf("******************************************\n\n");
}

//___________________________________________________________________________________________
void AliDhCorrelationExtraction::PrintSandBForNormal() {
  if(!fSignalCorrel || !fBackgrCorrel) {
    printf("Error! S and/or B not correctly set! Exiting...\n");
    return;
  }
  printf("******************************************\n");
  printf("   VALUES OF S AND B FOR NORMALIZATION\n\n");  
  printf(" NOTE: limits rounded to mass bin edges!\n");
  printf(" Signal region sigmas = %1.1f\n",fSignalSigmas);
  printf("******************************************\n");

  for(int i=0; i<fNpTbins; i++) {
    if(fIntegratePtBins && i>0) continue;
    printf("****************** Bin %d *****************\n",i+fFirstpTbin);
    printf("  Signal = %1.1f\n",fSignalCorrel[i]);
    printf("  Background = %1.1f\n",fBackgrCorrel[i]);
    printf("  SB scaling factor = %1.4f\n",fScaleFactor[i]);
  }
  printf("******************************************\n\n");
}

//___________________________________________________________________________________________
void AliDhCorrelationExtraction::RescaleSidebandsInMassBinEdges(Int_t i) {
  //Done only if external sideband ranges are passed (!fAutoSB) and a single bin is used for the SB (fSBSingleBin), and if bin counting approach is used
  //Use the integral of bkg fit function in the SB range and in the same range but matching mass bins (for the bin counting)

  Double_t fitInt_invmassSBwidth, fitInt_corrSBwidth;
  if(fDmesonSpecies!=kDStarD0pi) {
    fitInt_invmassSBwidth = fBkgFit[i]->Integral(fMassHisto[i]->GetXaxis()->GetBinLowEdge(fMassHisto[i]->FindBin(fRangesSB1L[i]+0.00001)),fMassHisto[i]->GetXaxis()->GetBinUpEdge(fMassHisto[i]->FindBin(fRangesSB1R[i]-0.00001))) + fBkgFit[i]->Integral(fMassHisto[i]->GetXaxis()->GetBinLowEdge(fMassHisto[i]->FindBin(fRangesSB2L[i]+0.00001)),fMassHisto[i]->GetXaxis()->GetBinUpEdge(fMassHisto[i]->FindBin(fRangesSB2R[i]-0.00001)));
    fitInt_corrSBwidth = fBkgFit[i]->Integral(fRangesSB1L[i],fRangesSB1R[i])+fBkgFit[i]->Integral(fRangesSB2L[i],fRangesSB2R[i]);
   } else {
    fitInt_invmassSBwidth = fBkgFit[i]->Integral(fMassHisto[i]->GetXaxis()->GetBinLowEdge(fMassHisto[i]->FindBin(fRangesSB1L[i]+0.00001)),fMassHisto[i]->GetXaxis()->GetBinUpEdge(fMassHisto[i]->FindBin(fRangesSB1R[i]-0.00001)));
    fitInt_corrSBwidth = fBkgFit[i]->Integral(fRangesSB1L[i],fRangesSB1R[i]); 
  }
  Double_t renorm =  fitInt_invmassSBwidth/fitInt_corrSBwidth; //trick to report the bin counting to "correaltion" SB edges (which are not inv.mass bin edges)

  if(TMath::Abs(renorm-1)>0.001) printf("\n***** WARNING!! Bin %d - SB ranges not matching mass binning! Readjusting the SB scaling factor by %1.3f *****\n\n",fFirstpTbin+i,renorm);
  fScaleFactor[i]*=renorm;
  if(fDebug>=2) {
    printf("Bin %d - renormalization of SB scaling factor by %1.3f (from %1.3f to %1.3f)\n",fFirstpTbin+i,renorm,fScaleFactor[i]/renorm,fScaleFactor[i]);
    printf("--- integrals of bkg function are %1.4f and %1.4f, region from %1.4f to %1.4f, %1.4f to %1.4f \n",fitInt_invmassSBwidth,fitInt_corrSBwidth,fMassHisto[i]->GetXaxis()->GetBinLowEdge(fMassHisto[i]->FindBin(fRangesSB1L[i]+0.00001)),fMassHisto[i]->GetXaxis()->GetBinUpEdge(fMassHisto[i]->FindBin(fRangesSB1R[i]-0.00001)),fMassHisto[i]->GetXaxis()->GetBinLowEdge(fMassHisto[i]->FindBin(fRangesSB2L[i]+0.00001)),fMassHisto[i]->GetXaxis()->GetBinUpEdge(fMassHisto[i]->FindBin(fRangesSB2R[i]-0.00001)));
  }
}

//___________________________________________________________________________________________
void AliDhCorrelationExtraction::SetHistoStyle(TH1F* &histo, Int_t colour) {
  histo->SetMarkerColor(colour);
  histo->SetLineColor(colour);
  histo->SetMarkerStyle(20);
  histo->SetLineStyle(1);
  histo->SetStats(kFALSE);
}

//__________________________________________
void AliDhCorrelationExtraction::SetHistoCorrStyle(TH1D* &histo) {
  histo->SetMinimum(TMath::Min(0.,histo->GetBinContent(histo->GetMinimumBin())*1.3));
  histo->SetMaximum(histo->GetBinContent(histo->GetMaximumBin())*1.3);
  histo->SetEntries(histo->GetEntries());
  histo->SetMarkerSize(0.8);
  histo->SetMarkerStyle(8);
  histo->SetLineColor(kBlack);
}

//__________________________________________
void AliDhCorrelationExtraction::DefinePaveText(TPaveText* &paveText, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, TString name) {
    
  paveText = new TPaveText(xmin,ymin,xmax,ymax,"NDC");
  paveText->SetBorderSize(0);
  paveText->SetFillColor(0);
  paveText->SetName(name);
  paveText->AddText(" ");
  if(fDmesonSpecies == kD0toKpi) paveText->AddText("D^{0} #rightarrowK^{-}#pi^{+}");
  if(fDmesonSpecies == kDStarD0pi) paveText->AddText("D^{*+}#rightarrow D^{0}#pi^{+}, D^{0} #rightarrowK^{-}#pi^{+}");
  if(fDmesonSpecies == kDplusKpipi) paveText->AddText("D^{+} #rightarrowK^{-}#pi^{+}#pi^{+}");
  paveText->AddText(" ");
}

//__________________________________________
void AliDhCorrelationExtraction::ClearObjects() {
   
  if(fDmesonFitterSignal) delete[] fDmesonFitterSignal;
  if(fDmesonFitterSignalError) delete[] fDmesonFitterSignalError; 
  if(fDmesonFitterBackground) delete[] fDmesonFitterBackground; 
  if(fDmesonFitterBackgroundError) delete[] fDmesonFitterBackgroundError; 
  if(fDMesonFitterSBCand) delete[] fDMesonFitterSBCand; 
  if(fDMesonFitterSBCandErr) delete[] fDMesonFitterSBCandErr;  
  if(fDmesonFitterMean) delete[] fDmesonFitterMean; 
  if(fDmesonFitterMeanError) delete[] fDmesonFitterMeanError;  

  if(fDmesonFitterSigma) delete[] fDmesonFitterSigma; 
  if(fDmesonFitterSigmaError) delete[] fDmesonFitterSigmaError; 
  if(fDmesonFitterSignificance) delete[] fDmesonFitterSignificance;  
  if(fDmesonFitterSignificanceError) delete[] fDmesonFitterSignificanceError;  
  if(fDmesonFitterSOverB) delete[] fDmesonFitterSOverB;  

  if(fSignalCorrel) delete[] fSignalCorrel; 
  if(fBackgrCorrel) delete[] fBackgrCorrel;  
  if(fRangesSignL) delete[] fRangesSignL; 
  if(fRangesSignR) delete[] fRangesSignR;
  if(fRangesSB1L) delete[] fRangesSB1L;
  if(fRangesSB1R) delete[] fRangesSB1R;  
  if(fDmesonSpecies != kDStarD0pi && fRangesSB2L) delete[] fRangesSB2L; 
  if(fDmesonSpecies != kDStarD0pi && fRangesSB2R) delete[] fRangesSB2R;  
  if(fScaleFactor) delete[] fScaleFactor;

  if(fMassFit) delete[] fMassFit; 
  if(fBkgFit) delete[] fBkgFit;

  if(fMassHisto) delete[] fMassHisto;  

}
