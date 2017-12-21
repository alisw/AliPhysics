/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/*____________________________________________________________
 | Class for the standard calculation and visualization of 
 | systematic uncertainties on the fit of the azimuthal correlations 
 |
 |  Author: Sandro Bjelogrlic, sandro.bjelogrlic@cern.ch
 |_____________________________________________________________*/


#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include "AliHFCorrelationUtils.h"
#include "AliHFCorrFitter.h"
#include "AliHFCorrFitSystematics.h"
#include <Riostream.h>
#include <TBufferFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TVirtualPad.h>
#include <TMath.h>
#include <TLatex.h>
#include <TColor.h>
#include <TClass.h>
#include <iostream>
#include <string>
#include <sstream>

//Correlation histogram should be normalised per bin width and No. of trigger----------

using std::cout;
using std::endl;
using std::stringstream;


ClassImp(AliHFCorrFitSystematics)

//___________________________________________________________
AliHFCorrFitSystematics::AliHFCorrFitSystematics():
fIsFittingTemplate(kFALSE),
fVecHisto(0),
  fVecHistoMinVar(0),
  fVecHistoMaxVar(0),
  fVecMinCorrelatedSyst(0),
  fVecMaxCorrelatedSyst(0),
  fCorrelationCanvas(0),
  fCorrelationGraphAsymmErrors(0),
  fFitFunctions(0),
  fVecLowEdgeDpt(0),
  fVecUpEdgeDpt(0),
  fv2DmesonVsPt(0),
  fVecSystMode(0),
  fFitFuncType(0),
  fMinBaselineEstimationRange(0),
  fMaxBaselineEstimationRange(0),
  fNMinPointsBaselineEstimationRange(0),
  fVecIsReference(0),
  fEntries(0),
  fVecSize(-1),
  fVecSystModesSize(-1),
  fDim(-1),
  fReferenceIndex(-1),
  fIsReferenceAlreadySet(kFALSE),
  fUseCorrelatedSystematics(kFALSE),
  fUseMaximumVariation(kFALSE),
  fIspPb(kFALSE),
  fUseCorrelatedSystematicsForWidths(kFALSE),
  fPlotV2SystSeparately(kFALSE),
  fSaveDotC(kFALSE),
  fSaveRoot(kFALSE),
  fSavePng(kFALSE),
  fSaveEps(kFALSE),
  fSavePdf(kFALSE),
  fV2DvsPt(kFALSE),
  fMinFitRange(0),
  fMaxFitRange(-1),
  fv2had(0),
  fv2Dmeson(0),
  fAssocTrackPtMin(-1),
  fAssocTrackPtMax(100),
  fMinReferenceBaselineEstimationRange(0),
  fMaxReferenceBaselineEstimationRange(0),
  fMinNPointsReferenceBaseline(0),
  fValueNSYield(0x0),
  fRatioNSYield(0x0),
  fValueNSSigma(0x0),
  fRatioNSSigma(0x0),
  fValueASYield(0x0),
  fRatioASYield(0x0),
  fValueASSigma(0x0),
  fRatioASSigma(0x0),
  fValuePedestal(0x0),
  fRatioPedestal(0x0),
  fValuev2NSYield(0x0),
  fValuev2NSSigma(0x0),
  fValuev2ASYield(0x0),
  fValuev2ASSigma(0x0),
  fValuev2Pedestal(0x0),
  fSystValuev2NSYield(0x0),
  fSystValuev2NSSigma(0x0),
  fSystValuev2ASYield(0x0),
  fSystValuev2ASSigma(0x0),
  fSystValuev2Pedestal(0x0),
  fValueSystematicBaselineNSYield(0x0),
  fValueSystematicBaselineNSSigma(0x0),
  fValueSystematicBaselineASYield(0x0),
  fValueSystematicBaselineASSigma(0x0),
  fValueSystematicBaselinePedestal(0x0),
  fValueSystematicNSYieldUp(0x0),
  fValueSystematicNSSigmaUp(0x0),
  fValueSystematicASYieldUp(0x0),
  fValueSystematicASSigmaUp(0x0),
  fValueSystematicPedestalUp(0x0),
  fValueSystematicNSYieldLow(0x0),
  fValueSystematicNSSigmaLow(0x0),
  fValueSystematicASYieldLow(0x0),
  fValueSystematicASSigmaLow(0x0),
  fValueSystematicPedestalLow(0x0),
  fOutputFileName(""),
  fOutputDirectory(""),
  fSystMode(kTransverse),
  fCombineSystematics(kRMS),
  fDmeson(AliHFCorrelationUtils::kDaverage),
  fFitter(0x0),
  fOutputFile(0x0),
  fReferenceHistoNSYield(0x0),
  fReferenceHistoNSSigma(0x0),
  fReferenceHistoASYield(0x0),
  fReferenceHistoASSigma(0x0),
  fReferenceHistoPedestal(0x0),
  fValueHistoNSYield(0x0),
  fRatioHistoNSYield(0x0),
  fValueHistoNSSigma(0x0),
  fRatioHistoNSSigma(0x0),
  fValueHistoASYield(0x0),
  fRatioHistoASYield(0x0),
  fValueHistoASSigma(0x0),
  fRatioHistoASSigma(0x0),
  fValueHistoPedestal(0x0),
  fRatioHistoPedestal(0x0),
  fRMSHistoNSYield(0x0),
  fRMSHistoNSSigma(0x0),
  fRMSHistoASYield(0x0),
  fRMSHistoASSigma(0x0),
  fRMSHistoPedestal(0x0),
  fSystematicSourcesNSYield(0x0),
  fSystematicSourcesNSSigma(0x0),
  fSystematicSourcesASYield(0x0),
  fSystematicSourcesASSigma(0x0),
  fSystematicSourcesPedestal(0x0),
  fCanvasSystematicSourcesNSYield(0x0),
  fCanvasSystematicSourcesNSSigma(0x0),
  fCanvasSystematicSourcesASYield(0x0),
  fCanvasSystematicSourcesASSigma(0x0),
  fCanvasSystematicSourcesPedestal(0x0),
  fCanvasTotalSystematicSourcesNSYield(0x0),
  fCanvasTotalSystematicSourcesNSSigma(0x0),
  fCanvasTotalSystematicSourcesASYield(0x0),
  fCanvasTotalSystematicSourcesASSigma(0x0),
  fCanvasTotalSystematicSourcesPedestal(0x0),
  fCanvasFinalTrendNSYield(0x0),
  fCanvasFinalTrendNSSigma(0x0),
  fCanvasFinalTrendASYield(0x0),
  fCanvasFinalTrendASSigma(0x0),
  fCanvasFinalTrendPedestal(0x0),
  fCanvasVariationBaselineTrendPedestal(0x0),
  fFinalTrendNSYield(0x0),
  fFinalTrendNSSigma(0x0),
  fFinalTrendASYield(0x0),
  fFinalTrendASSigma(0x0),
  fFinalTrendPedestal(0x0),
  fFullSystematicsNSYield(0x0),
  fv2SystematicsNSYield(0x0),
  fFullSystematicsNSSigma(0x0),
  fv2SystematicsNSSigma(0x0),
  fFullSystematicsASYield(0x0),
  fv2SystematicsASYield(0x0),
  fFullSystematicsASSigma(0x0),
  fv2SystematicsASSigma(0x0),
  fFullSystematicsPedestal(0x0),
  fv2SystematicsPedestal(0x0),
  fBaselineVariationSystematicsPedestal(0x0),
  fCanvasRefernce(0x0),
  fCanvasFitting(0x0)
{
    // default constructor
  
}

//_________________________|Destructor|___________________
AliHFCorrFitSystematics::~AliHFCorrFitSystematics()
{
    Info("AliHFCorrFitSystematics.cxx","Destructor is calling");
    
    if(fFitter){delete fFitter;}
    if(fOutputFile){delete fOutputFile;}
    // add a destructro for fRatioHisto
    
    if(fValueNSYield){delete fValueNSYield;}
    if(fRatioNSYield){delete fRatioNSYield;}
    if(fValueNSSigma){delete fValueNSSigma;}
    if(fRatioNSSigma){delete fRatioNSSigma;}
    if(fValueASYield){delete fValueASYield;}
    if(fRatioASYield){delete fRatioASYield;}
    if(fValueASSigma){delete fValueASSigma;}
    if(fRatioASSigma){delete fRatioASSigma;}
    if(fValuePedestal){delete fValuePedestal;}
    if(fRatioPedestal){delete fRatioPedestal;}
   
    /* add destructors
    fOutputFile(0x0),
    fReferenceHistoNSYield(0x0),
    fReferenceHistoNSSigma(0x0),
    fReferenceHistoASYield(0x0),
    fReferenceHistoASSigma(0x0),
    fReferenceHistoPedestal(0x0),
 */
     if(fOutputFile){delete fOutputFile;}
    
    if(fReferenceHistoNSYield){delete fReferenceHistoNSYield;}
    if(fReferenceHistoNSSigma){delete fReferenceHistoNSSigma;}
    if(fReferenceHistoASYield){delete fReferenceHistoASYield;}
    if(fReferenceHistoASSigma){delete fReferenceHistoASSigma;}
    if(fReferenceHistoPedestal){delete fReferenceHistoPedestal;}
    
    if(fValueHistoNSYield){delete fValueHistoNSYield;}
    if(fRatioHistoNSYield){delete fRatioHistoNSYield;}
    if(fValueHistoNSSigma){delete fValueHistoNSSigma;}
    if(fRatioHistoNSSigma){delete fRatioHistoNSSigma;}
    if(fValueHistoASYield){delete fValueHistoASYield;}
    if(fRatioHistoASYield){delete fRatioHistoASYield;}
    if(fValueHistoASSigma){delete fValueHistoASSigma;}
    if(fRatioHistoASSigma){delete fRatioHistoASSigma;}
    if(fValueHistoPedestal){delete fValueHistoPedestal;}
    if(fRatioHistoPedestal){delete fRatioHistoPedestal;}
    
    if(fBaselineVariationSystematicsPedestal){delete fBaselineVariationSystematicsPedestal;}
    if(fCanvasFitting){delete fCanvasFitting;}
    if(fCanvasVariationBaselineTrendPedestal){delete fCanvasVariationBaselineTrendPedestal;}
}

//___________________________________________________________
Bool_t AliHFCorrFitSystematics::AddHistoToFit(TString path){
    
    //change me if needed
   return AddHistoToFit(path,"cDraw","fhDaverage","fhDaverageVariedMin","fhDaverageVariedMax");
    
}

//___________________________________________________________
Bool_t AliHFCorrFitSystematics::AddHistoToFit(TString path, TString cname, TString  havname, TString havmin, TString havmax)
{
    
    TFile * file = TFile::Open(path.Data());
    if(!file->IsOpen()){
        std::cout << "Cannot open file " << path << " - check your input! " << std::endl; return kFALSE;
    }
    TH1D * hist;
    TList * thelist;
    TCanvas * c = (TCanvas*)file->Get(cname.Data());
    if(c){
      thelist = (TList *)c->GetListOfPrimitives();
      hist = (TH1D*)c->FindObject(havname.Data());
      if(!hist){   std::cout << "Cannot extract histogram  " << havname << "  from canvas" << cname << " from file " << path << " - check your input! " << std::endl; return kFALSE;
      }
    }
    else {
      cout<<"Canvas with name: < "<<cname<<" > does not exist. Trying to take the histo directly from file "<<std::endl;
      c=0x0;
      hist=(TH1D*)file->Get(havname.Data());
      if(!hist){   
	std::cout << "Cannot extract histogram " << havname << "  from file " << path << " - check your input! " << std::endl; return kFALSE;
      }
    }

    

    if(!hist){   std::cout << "Cannot extract histogram of average " << havname << "  from file " << path << " - check your input! " << std::endl; return kFALSE;}
    
    TH1D * histminvar = (TH1D *)file->Get(havmin.Data());
    if(!histminvar){   
      if(!fIsFittingTemplate){
	std::cout << "Cannot extract histogram of min variation " << havmin << " from file " << path << " - check your input! " << std::endl; return kFALSE;
      }
      histminvar=0x0;
    }
    
    TH1D * histmaxvar = (TH1D *)file->Get(havmax.Data());
    if(!histmaxvar){         
      if(!fIsFittingTemplate){
	std::cout << "Cannot extract histogram of max variation " << havmax << " from file " << path << " - check your input! " << std::endl; return kFALSE;
      }
      histmaxvar=0x0;
    }
      
    hist->GetXaxis()->SetTitle("#Delta#varphi (rad)");
    hist->GetYaxis()->SetTitle("#frac{1}{N_{D}}#frac{dN}{d#Delta#varphi}(rad^{-1})");
    hist->GetXaxis()->SetTitleSize(0.03);
    hist->GetYaxis()->SetTitleSize(0.03);
    hist->GetXaxis()->SetTitleOffset(1.25);
    hist->GetYaxis()->SetTitleOffset(1.25);
    hist->GetXaxis()->SetLabelSize(0.03);
    hist->GetYaxis()->SetLabelSize(0.03);
    
    fVecHisto.push_back(hist);
    fVecHistoMinVar.push_back(histminvar);
    fVecHistoMaxVar.push_back(histmaxvar);
    fCorrelationCanvas.push_back(c);
    
    if(thelist && !fIsFittingTemplate){
      TGraphAsymmErrors * tgraphassocerror = (TGraphAsymmErrors *)thelist->At(2);
      tgraphassocerror->SetName(Form("grapherror%d",fEntries));
      fCorrelationGraphAsymmErrors.push_back(tgraphassocerror);     
      fEntries++;
    }

    /*delete file;
    delete c;
    delete list;
    delete hist;
    delete histminvar;
    delete histmaxvar;
    */
    
    return kTRUE;
    // find a fancy way of extracting the correalted systematics - now are set by hand
    
  //  TLatex * latex = (TLatex *)list->At(7);
    
  //  Int_t maxsyst, minsyst;
  //  string dummy_1, dummy_2, dummy_3;
    
    
    //_________________________________
  //  TString text = latex->GetTitle();
   
    
    //stringstream(text) >> maxsyst;
    
   // char *mystring = (char *)text.Data();
   // int size = strlen(mystring);
   // cout << "size " <<  size << endl;
    
    //for (char &c : text)
   // {
   //     if(c=="+")cout <<  c << endl;
   // }
    
   //  std:: cout << "latex = " << text << std::endl;
    
    
  //   cout << "Max syst = " << maxsyst << endl;
  //  return kTRUE;
    
}
//_______________________________________________________________________________
void AliHFCorrFitSystematics::AddSystematicMode(SystematicModes mode, Bool_t isReference, Double_t min, Double_t max, Int_t minNpoints,AliHFCorrFitter::FunctionType fitfunctype){
    if(fIsReferenceAlreadySet && isReference) {std::cout << " " << std::endl; std::cout << "You hava already set a reference histo!" << std::endl; return;}
    fVecSystMode.push_back(mode);
    fFitFuncType.push_back(fitfunctype);
    fVecIsReference.push_back(isReference);
    if(isReference) fIsReferenceAlreadySet = kTRUE;
    
    if(mode == kTransverse && !isReference){
      fMinBaselineEstimationRange.push_back(min);
      fMaxBaselineEstimationRange.push_back(max);
      //  cout << "Range set " << min << ","<< max << endl;
    }
    else if(mode == kTransverseUppStatUnc && !isReference){
      fMinBaselineEstimationRange.push_back(min);
      fMaxBaselineEstimationRange.push_back(max);
      //  cout << "Range set " << min << ","<< max << endl;
    }    
    else if(mode == kTransverseLowStatUnc && !isReference){
      fMinBaselineEstimationRange.push_back(min);
      fMaxBaselineEstimationRange.push_back(max);
      //  cout << "Range set " << min << ","<< max << endl;
    }   
    else if(mode == kTransverse && isReference){
      fMinBaselineEstimationRange.push_back(min);
      fMaxBaselineEstimationRange.push_back(max);
      fMinReferenceBaselineEstimationRange=min;
      fMaxReferenceBaselineEstimationRange=max;      
    }
    else{
      fMinBaselineEstimationRange.push_back(fMinReferenceBaselineEstimationRange);
      fMaxBaselineEstimationRange.push_back(fMaxReferenceBaselineEstimationRange);
      //cout << "Range set else" << fMinReferenceBaselineEstimationRange << ","<< fMaxReferenceBaselineEstimationRange << endl;
    }

    if(mode == kNLowest && !isReference){
      fNMinPointsBaselineEstimationRange.push_back(minNpoints);
    }
    else if(mode == kNLowest && isReference){
      fNMinPointsBaselineEstimationRange.push_back(minNpoints);
      fMinNPointsReferenceBaseline=minNpoints;
    }
    else{
      fNMinPointsBaselineEstimationRange.push_back(fMinNPointsReferenceBaseline);
    }
}
//_______________________________________________________________________________
void AliHFCorrFitSystematics::AddSystematicMode(SystematicModes mode, Bool_t isReference){
  AddSystematicMode(mode,isReference, fMinReferenceBaselineEstimationRange, fMaxReferenceBaselineEstimationRange,fMinNPointsReferenceBaseline);
}

/*
//_______________________________________________________________________________
void AliHFCorrFitSystematics::AddSystematicMode(SystematicModes mode, Bool_t isReference){
    
    
    if(fIsReferenceAlreadySet && isReference) {std::cout << " " << std::endl; std::cout << "You hava already set a reference histo!" << std::endl; return;}
     fVecSystMode.push_back(mode);
     fVecIsReference.push_back(isReference);
    if(isReference) fIsReferenceAlreadySet = kTRUE;
    
  //  if(!(mode == kTransverse && !isReference)){
  //      fMinBaselineEstimationRange.push_back(fMinReferenceBaselineEstimationRange);
  //      fMaxBaselineEstimationRange.push_back(fMaxReferenceBaselineEstimationRange);
  //  }
    
}
//_______________________________________________________________________________
void AliHFCorrFitSystematics::AddSystematicMode(SystematicModes mode, Bool_t isReference, Double_t min, Double_t max){
    
    if(mode == kTransverse && !isReference){
        fMinBaselineEstimationRange.push_back(min);
        fMaxBaselineEstimationRange.push_back(max);
      //  cout << "Range set " << min << ","<< max << endl;
    }
    else{
        fMinBaselineEstimationRange.push_back(fMinReferenceBaselineEstimationRange);
        fMaxBaselineEstimationRange.push_back(fMaxReferenceBaselineEstimationRange);
        //cout << "Range set else" << fMinReferenceBaselineEstimationRange << ","<< fMaxReferenceBaselineEstimationRange << endl;
    }
    AddSystematicMode(mode,isReference);
    
}
 */
//_______________________________________________________________________________
void AliHFCorrFitSystematics::CheckBaselineRanges(){
    std::cout << "======================== " <<fMinBaselineEstimationRange.size() <<  std::endl;
    for(Int_t k =0; k<(Int_t)fMinBaselineEstimationRange.size(); k++){
        std::cout << "Baseline in range (" << k << "): " << fMinBaselineEstimationRange[k]/TMath::Pi() << "*pi - " << fMaxBaselineEstimationRange[k]/TMath::Pi() << "*pi" << std::endl;
    }
    for(Int_t k=0; k<(Int_t)fNMinPointsBaselineEstimationRange.size();k++){
      std::cout << "Min Npoints for baseline (" << k <<"):" << fNMinPointsBaselineEstimationRange[k]<<std::endl;
    }
}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::CheckHisto(Int_t i){
    
    TCanvas * c = new TCanvas();
    c->cd();
    fVecHisto[i]->Draw();
    
}


//_______________________________________________________________________________
Bool_t AliHFCorrFitSystematics::SetUpSystematicsFitter(){
    
    Bool_t flag1 = CheckSize();
    Bool_t flag2 = CheckDiffSystematicsSize();
    
    if(!flag1) {
        std::cout << "AliHFCorrFitSystematics::SetUpSystematicsFitter() - something is wrong in setting up the D meson pt" << std::endl;
        return flag1;
    }
    if(!flag2) {
        std::cout << "AliHFCorrFitSystematics::SetUpSystematicsFitter() - something is wrong in setting up systematics modes" << std::endl;
         std::cout << "Did you use the function AddSystematicMode(SystematicModes mode)?" << std::endl;
        return flag2;
    }
    
    if(!fOutputFile){
        std::cout << " " << std::endl;
         std::cout << "No OutPut file created" << std::endl;
         std::cout << "Creating it with default name FitSystematicsOutput.root" << std::endl;
        CreateOutputFile("FitSystematicsOutput_1.root");
    }
    
    
    fDim = fVecSize * fVecSystModesSize;
    
    //initialize the ratop histos

    const int size = fVecSize+1;
    cout << "fVecSize = " << fVecSize << endl;
    Float_t ptbins[size];
    for(Int_t i = 0; i<fVecSize; i++){
        ptbins[i] = fVecLowEdgeDpt[i];
    }
    ptbins[fVecSize] = fVecUpEdgeDpt[fVecSize-1];
    
    Float_t *ptbinspointer = ptbins;
    
   // for(Int_t i = 0; i<size; i++){
     
     //   cout << i << " " << ptbins[i] << endl;
   // }
    
    fReferenceHistoNSYield = new TH1D("ReferenceHistoNSYield","ReferenceHistoNSYield",fVecSize,ptbinspointer);
    fReferenceHistoNSSigma = new TH1D("ReferenceHistoNSSigma","ReferenceHistoNSSigma",fVecSize,ptbinspointer);
    fReferenceHistoASYield= new TH1D("ReferenceHistoASYield","ReferenceHistoASYield",fVecSize,ptbinspointer);
    fReferenceHistoASSigma = new TH1D("ReferenceHistoASSigma","ReferenceHistoASSigma",fVecSize,ptbinspointer);
    fReferenceHistoPedestal = new TH1D("ReferenceHistoPedestal","ReferenceHistoPedestal",fVecSize,ptbinspointer);

  
    fValueHistoNSYield = new TH1D("ValueHistoNSYield","ValueHistoNSYield",fVecSize,ptbinspointer);
    fRatioHistoNSYield = new TH1D("RatioHistoNSYield","RatioHistoNSYield",fVecSize,ptbinspointer);
    fValueHistoNSSigma = new TH1D("ValueHistoNSSigma","ValueHistoNSSigma",fVecSize,ptbinspointer);
    fRatioHistoNSSigma = new TH1D("RatioHistoNSSigma","RatioHistoNSSigma",fVecSize,ptbinspointer);
    fValueHistoASYield = new TH1D("ValueHistoASYield","ValueHistoASYield",fVecSize,ptbinspointer);
    fRatioHistoASYield = new TH1D("RatioHistoASYield","RatioHistoASYield",fVecSize,ptbinspointer);
    fValueHistoASSigma = new TH1D("ValueHistoASSigma","ValueHistoASSigma",fVecSize,ptbinspointer);
    fRatioHistoASSigma = new TH1D("RatioHistoASSigma","RatioHistoASSigma",fVecSize,ptbinspointer);
    fValueHistoPedestal = new TH1D("ValueHistoPedestal","ValueHistoPedestal",fVecSize,ptbinspointer);
    fRatioHistoPedestal = new TH1D("RatioHistoPedestal","RatioHistoPedestal",fVecSize,ptbinspointer);
    
   
    if(fVecSize == 1)fCanvasRefernce = new TCanvas("CanvasRefernce","CanvasRefernce",0,0,1500,1000);
    if(fVecSize == 2){fCanvasRefernce = new TCanvas("CanvasRefernce","CanvasRefernce",0,0,1500,1000); fCanvasRefernce->Divide(2,1);}
    if(fVecSize == 3){fCanvasRefernce = new TCanvas("CanvasRefernce","CanvasRefernce",0,0,1500,1000); fCanvasRefernce->Divide(3,1);}
    if(fVecSize == 4){fCanvasRefernce = new TCanvas("CanvasRefernce","CanvasRefernce",0,0,1600,1000); fCanvasRefernce->Divide(2,2);}
    if(fVecSize == 5){fCanvasRefernce = new TCanvas("CanvasRefernce","CanvasRefernce",0,0,1600,1000); fCanvasRefernce->Divide(3,2);}
    if(fVecSize == 6){fCanvasRefernce = new TCanvas("CanvasRefernce","CanvasRefernce",0,0,1600,1000); fCanvasRefernce->Divide(3,2);}
    if(fVecSize == 7){fCanvasRefernce = new TCanvas("CanvasRefernce","CanvasRefernce",0,0,1600,1000); fCanvasRefernce->Divide(4,2);}
    if(fVecSize == 8){fCanvasRefernce = new TCanvas("CanvasRefernce","CanvasRefernce",0,0,1600,1000); fCanvasRefernce->Divide(4,2);}
    if(fVecSize == 9){fCanvasRefernce = new TCanvas("CanvasRefernce","CanvasRefernce",0,0,1600,1000); fCanvasRefernce->Divide(4,2);}
    
    
    fCanvasFitting = new TCanvas*[fVecSystModesSize];
    
    /*
     
    if(fVecSize == 1)fCanvasFitting = new TCanvas("cFitting","cFitting",0,0,1000,1000);
    if(fVecSize == 2){fCanvasFitting = new TCanvas("cFitting","cFitting",0,0,1500,1000); fCanvasFitting->Divide(2,1);}
    if(fVecSize == 3){fCanvasFitting = new TCanvas("cFitting","cFitting",0,0,1500,1000); fCanvasFitting->Divide(3,1);}
    if(fVecSize == 4){fCanvasFitting = new TCanvas("cFitting","cFitting",0,0,1600,1000); fCanvasFitting->Divide(2,2);}
    if(fVecSize == 5){fCanvasFitting = new TCanvas("cFitting","cFitting",0,0,1600,1000); fCanvasFitting->Divide(3,2);}
    if(fVecSize == 6){fCanvasFitting = new TCanvas("cFitting","cFitting",0,0,1600,1000); fCanvasFitting->Divide(3,2);}
    if(fVecSize == 7){fCanvasFitting = new TCanvas("cFitting","cFitting",0,0,1600,1000); fCanvasFitting->Divide(4,2);}
    if(fVecSize == 8){fCanvasFitting = new TCanvas("cFitting","cFitting",0,0,1600,1000); fCanvasFitting->Divide(4,2);}
    if(fVecSize == 9){fCanvasFitting = new TCanvas("cFitting","cFitting",0,0,1600,1000); fCanvasFitting->Divide(4,2);}*/
    if(fVecSize>9){std::cout << "Sorry I cannot handle more than 9 pT bins of the D meson - requires modifications!" << std::endl; return kFALSE;}
    
    
    fValueNSYield = new Double_t[fDim];
    fRatioNSYield = new Double_t[fDim];
    fValueNSSigma = new Double_t[fDim];
    fRatioNSSigma = new Double_t[fDim];
    fValueASYield = new Double_t[fDim];
    fRatioASYield = new Double_t[fDim];
    fValueASSigma = new Double_t[fDim];
    fRatioASSigma = new Double_t[fDim];
    fValuePedestal = new Double_t[fDim];
    fRatioPedestal = new Double_t[fDim];
    
    
   
    
    fValueSystematicBaselineNSYield = new Double_t[fVecSize];
    fValueSystematicBaselineNSSigma = new Double_t[fVecSize];
    fValueSystematicBaselineASYield = new Double_t[fVecSize];
    fValueSystematicBaselineASSigma = new Double_t[fVecSize];
    fValueSystematicBaselinePedestal = new Double_t[fVecSize];

    fValueSystematicBaseline_FromBaselStatUp_NSYield = new Double_t[fVecSize];
    fValueSystematicBaseline_FromBaselStatUp_NSSigma = new Double_t[fVecSize];
    fValueSystematicBaseline_FromBaselStatUp_ASYield = new Double_t[fVecSize];
    fValueSystematicBaseline_FromBaselStatUp_ASSigma = new Double_t[fVecSize];
    fValueSystematicBaseline_FromBaselStatUp_Pedestal = new Double_t[fVecSize];
    
    fValueSystematicBaseline_FromBaselStatLo_NSYield = new Double_t[fVecSize];
    fValueSystematicBaseline_FromBaselStatLo_NSSigma = new Double_t[fVecSize];
    fValueSystematicBaseline_FromBaselStatLo_ASYield = new Double_t[fVecSize];
    fValueSystematicBaseline_FromBaselStatLo_ASSigma = new Double_t[fVecSize];
    fValueSystematicBaseline_FromBaselStatLo_Pedestal = new Double_t[fVecSize];

    fRMSRelative_NSYield = new Double_t[fVecSize];
    fRMSRelative_NSSigma = new Double_t[fVecSize];
    fRMSRelative_ASYield = new Double_t[fVecSize];
    fRMSRelative_ASSigma = new Double_t[fVecSize];
    fRMSRelative_Pedestal = new Double_t[fVecSize];

    fValueSystematicNSYieldUp = new Double_t[fVecSize];
    fValueSystematicNSSigmaUp = new Double_t[fVecSize];
    fValueSystematicASYieldUp = new Double_t[fVecSize];
    fValueSystematicASSigmaUp = new Double_t[fVecSize];
    fValueSystematicPedestalUp = new Double_t[fVecSize];
    
    fValueSystematicNSYieldLow = new Double_t[fVecSize];
    fValueSystematicNSSigmaLow = new Double_t[fVecSize];
    fValueSystematicASYieldLow = new Double_t[fVecSize];
    fValueSystematicASSigmaLow = new Double_t[fVecSize];
    fValueSystematicPedestalLow = new Double_t[fVecSize];
    
    
    if(fCombineSystematics == kRMS || fCombineSystematics == kEnvelope_RMS_BaselStat){
    fRMSHistoNSYield = new TH1D *[fVecSize];
    fRMSHistoNSSigma = new TH1D *[fVecSize];
    fRMSHistoASYield = new TH1D *[fVecSize];
    fRMSHistoASSigma = new TH1D *[fVecSize];
        fRMSHistoPedestal = new TH1D *[fVecSize];
        
        for(Int_t k=0; k<fVecSize;k++){
            fRMSHistoNSYield[k] = new TH1D(Form("RMSHistoNSYieldBin%d",k),Form("RMSHistoNSYieldBin%d",k),1000,0,10);
            fRMSHistoNSSigma[k] =  new TH1D(Form("RMSHistoNSSigmaBin%d",k),Form("RMSHistoNSSigmaBin%d",k),1000,0,3);
            fRMSHistoASYield[k] = new TH1D(Form("RMSHistoASYieldBin%d",k),Form("RMSHistoASYieldBin%d",k),1000,0,10);
            fRMSHistoASSigma[k] =  new TH1D(Form("RMSHistoASSigmaBin%d",k),Form("RMSHistoASSigmaBin%d",k),1000,0,3);
            fRMSHistoPedestal[k] = new TH1D(Form("RMSHistoPedestalBin%d",k),Form("RMSHistoPedestalBin%d",k),1000,0,10);
        }
        
        
    }

    
    fSystematicSourcesNSYield = new TH1D("SystematicSourcesNSYield","SystematicSourcesNSYield",fVecSize,ptbinspointer);
    fSystematicSourcesNSSigma = new TH1D("SystematicSourcesNSSigma","SystematicSourcesNSSigma",fVecSize,ptbinspointer);
    fSystematicSourcesASYield= new TH1D("SystematicSourcesASYield","SystematicSourcesASYield",fVecSize,ptbinspointer);
    fSystematicSourcesASSigma = new TH1D("SystematicSourcesASSigma","SystematicSourcesASSigma",fVecSize,ptbinspointer);
    fSystematicSourcesPedestal = new TH1D("SystematicSourcesPedestal","SystematicSourcesPedestal",fVecSize,ptbinspointer);
    
    fSystematicSourcesNSYield->GetYaxis()->SetRangeUser(0,3);
    fSystematicSourcesNSSigma->GetYaxis()->SetRangeUser(0,3);
    fSystematicSourcesASYield->GetYaxis()->SetRangeUser(0,3);
    fSystematicSourcesASSigma->GetYaxis()->SetRangeUser(0,3);
    fSystematicSourcesPedestal->GetYaxis()->SetRangeUser(0,3);
    
    
    //if(fCombineSystematics == kRMS){
    fCanvasSystematicSourcesNSYield = new TCanvas("BaselineSystematicSourcesNSYield","BaselineSystematicSourcesNSYield",0,0,1000,800);
    fCanvasSystematicSourcesNSSigma = new TCanvas("BaselineSystematicSourcesNSSigma","BaselineSystematicSourcesNSSigma",0,0,1000,800);
    fCanvasSystematicSourcesASYield = new TCanvas("BaselineSystematicSourcesASYield","BaselineSystematicSourcesASYield",0,0,1000,800);
    fCanvasSystematicSourcesASSigma = new TCanvas("BaselineSystematicSourcesASSigma","BaselineSystematicSourcesASSigma",0,0,1000,800);
    fCanvasSystematicSourcesPedestal = new TCanvas("BaselineSystematicSourcesPedestal","BaselineSystematicSourcesPedestal",0,0,1000,800);
    
    fCanvasTotalSystematicSourcesNSYield = new TCanvas("TotalSystematicSourcesNSYield","TotalSystematicSourcesNSYield",0,0,1000,800);
    fCanvasTotalSystematicSourcesNSSigma = new TCanvas("TotalSystematicSourcesNSSigma","TotalSystematicSourcesNSSigma",0,0,1000,800);
    fCanvasTotalSystematicSourcesASYield = new TCanvas("TotalSystematicSourcesASYield","TotalSystematicSourcesASYield",0,0,1000,800);
    fCanvasTotalSystematicSourcesASSigma = new TCanvas("TotalSystematicSourcesASSigma","TotalSystematicSourcesASSigma",0,0,1000,800);
    fCanvasTotalSystematicSourcesPedestal = new TCanvas("TotalSystematicSourcesPedestal","TotalSystematicSourcesPedestal",0,0,1000,800);
    
    fCanvasSystematicSourcesNSYield->SetTicks();
    fCanvasSystematicSourcesNSSigma->SetTicks();
    fCanvasSystematicSourcesASYield->SetTicks();
    fCanvasSystematicSourcesASSigma->SetTicks();
    fCanvasSystematicSourcesPedestal->SetTicks();
    
    fCanvasTotalSystematicSourcesNSYield->SetTicks();
    fCanvasTotalSystematicSourcesNSSigma->SetTicks();
    fCanvasTotalSystematicSourcesASYield->SetTicks();
    fCanvasTotalSystematicSourcesASSigma->SetTicks();
    fCanvasSystematicSourcesPedestal->SetTicks();
    //}
    /*else{
        fCanvasSystematicSourcesNSYield = new TCanvas("SystematicSourcesNSYield","SystematicSourcesNSYield",0,0,1200,800); fCanvasSystematicSourcesNSYield->Divide(2,1);
        fCanvasSystematicSourcesNSSigma = new TCanvas("SystematicSourcesNSSigma","SystematicSourcesNSSigma",0,0,1200,800); fCanvasSystematicSourcesNSSigma->Divide(2,1);
        fCanvasSystematicSourcesASYield = new TCanvas("SystematicSourcesASYield","SystematicSourcesASYield",0,0,1200,800); fCanvasSystematicSourcesASYield->Divide(2,1);
        fCanvasSystematicSourcesASSigma = new TCanvas("SystematicSourcesASSigma","SystematicSourcesASSigma",0,0,1200,800); fCanvasSystematicSourcesASSigma->Divide(2,1);
        fCanvasSystematicSourcesPedestal = new TCanvas("SystematicSourcesPedestal","SystematicSourcesPedestal",0,0,1200,800); fCanvasSystematicSourcesPedestal->Divide(2,1);
    }*/

    
    /*fFullSystematics = new TBox*[fVecSize];
    fv2Systematics= new TBox*[fVecSize];
    fFullSystematicsNSSigma = new TBox*[fVecSize];
    fv2SystematicsNSSigma = new TBox*[fVecSize];
    fFullSystematicsASYield = new TBox*[fVecSize];
    fv2SystematicsASYield = new TBox*[fVecSize];
    fFullSystematicsASSigma = new TBox*[fVecSize];
    fv2SystematicsASSigma = new TBox*[fVecSize];
    fFullSystematicsPedestal = new TBox*[fVecSize];
    fv2SystematicsPedestal = new TBox*[fVecSize];
    */
 
    
    
    for(Int_t iSystMode =0; iSystMode<fVecSystModesSize; iSystMode++){
    if(fVecIsReference[iSystMode]) fReferenceIndex = iSystMode;
    }
    return kTRUE;
}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::Fitv2Systematics(){
    
    
    TPaveText * pave = NULL;
    
    fValuev2NSYield = new Double_t[fVecSize];
    fValuev2NSSigma = new Double_t[fVecSize];
    fValuev2ASYield = new Double_t[fVecSize];
    fValuev2ASSigma = new Double_t[fVecSize];
    fValuev2Pedestal = new Double_t[fVecSize];
    
    fSystValuev2NSYield = new Double_t[fVecSize];
    fSystValuev2NSSigma = new Double_t[fVecSize];
    fSystValuev2ASYield = new Double_t[fVecSize];
    fSystValuev2ASSigma = new Double_t[fVecSize];
    fSystValuev2Pedestal = new Double_t[fVecSize];
    
    TF1* v2=NULL;
    
    Double_t baseline = 0;

  
    TCanvas * cv2 = NULL;
    
    
    if(fVecSize == 1)cv2 = new TCanvas("cv2","cv2",0,0,1000,1000);
    if(fVecSize == 2){cv2 = new TCanvas("cv2","cv2",0,0,1500,1000); cv2->Divide(2,1);}
    if(fVecSize == 3){cv2 = new TCanvas("cv2","cv2",0,0,1500,1000); cv2->Divide(3,1);}
    if(fVecSize == 4){cv2 = new TCanvas("cv2","cv2",0,0,1600,1000); cv2->Divide(2,2);}
    if(fVecSize == 5){cv2 = new TCanvas("cv2","cv2",0,0,1600,1000); cv2->Divide(3,2);}
    if(fVecSize == 6){cv2 = new TCanvas("cv2","cv2",0,0,1600,1000); cv2->Divide(3,2);}
    if(fVecSize == 7){cv2 = new TCanvas("cv2","cv2",0,0,1600,1000); cv2->Divide(4,2);}
    if(fVecSize == 8){cv2 = new TCanvas("cv2","cv2",0,0,1600,1000); cv2->Divide(4,2);}
    if(fVecSize == 9){cv2 = new TCanvas("cv2","cv2",0,0,1600,1000); cv2->Divide(4,2);}
    
    TH1D ** v2subtractedhisto = new TH1D*[fVecSize];
    TH1D ** inputhisto = new TH1D*[fVecSize];
    Double_t value = 0;

    for(Int_t iPtBin =0; iPtBin<fVecSize; iPtBin++){
     
        pave = new TPaveText(0.15,0.6,0.45,0.8,"NDC");
        pave->SetName(Form("pave_%d",iPtBin));
        pave->AddText(Form("%.1f < p_{T}(D) < %.1f GeV/c",fVecLowEdgeDpt[iPtBin],fVecUpEdgeDpt[iPtBin]));
        if(fAssocTrackPtMax>50) pave->AddText(Form("p_{T}(hadron) > %.1f GeV/c",fAssocTrackPtMin));
        else pave->AddText(Form("%.1f <p_{T}(hadron) < %.1f GeV/c",fAssocTrackPtMin, fAssocTrackPtMin));
        pave->SetFillStyle(0);
        pave->SetBorderSize(0);
        pave->SetTextSize(0.03);
        
        baseline = fReferenceHistoPedestal->GetBinContent(iPtBin+1);

        v2subtractedhisto[iPtBin] = (TH1D*)fVecHisto[iPtBin]->Clone(Form("v2subtrac_%d",iPtBin));
        inputhisto[iPtBin] = (TH1D*)fVecHisto[iPtBin]->Clone(Form("input_%d",iPtBin));
        v2subtractedhisto[iPtBin]->Reset();
        inputhisto[iPtBin]->SetMarkerColor(4);
        inputhisto[iPtBin]->SetLineColor(4);
        
        
        v2 = new TF1(Form("v2Modulation_%d",iPtBin),"[0]*(1+2*[1]*[2]*TMath::Cos(2*x[0]))",fMinFitRange,fMaxFitRange);

        v2->SetParameter(0,baseline);
        v2->SetParameter(1,fv2had);
        if(!fV2DvsPt) v2->SetParameter(2,fv2Dmeson);
        else v2->SetParameter(2,fv2DmesonVsPt.at(iPtBin));
        v2->SetLineColor(4);
        v2->SetLineStyle(4);

        printf("v2 for D-meson = %1.3f, for hadron = %1.3f\n",fv2DmesonVsPt.at(iPtBin),fv2had);

        for(Int_t iBin = 1; iBin <= fVecHisto[iPtBin]->GetNbinsX();iBin++){
            value = inputhisto[iPtBin]->GetBinContent(iBin);
            value -= v2->Eval(inputhisto[iPtBin]->GetBinCenter(iBin));
            value += baseline;
            
            v2subtractedhisto[iPtBin]->SetBinContent(iBin,value);
            v2subtractedhisto[iPtBin]->SetBinError(iBin,inputhisto[iPtBin]->GetBinError(iBin));
        }
     
        cv2->cd(iPtBin+1);
       
        
        
        fFitter = new AliHFCorrFitter((TH1F*)v2subtractedhisto[iPtBin],fMinFitRange,fMaxFitRange);
	if(fVecSystMode[fReferenceIndex]==kNLowest){
	  fFitter->SetFixBasetype(-1*fNMinPointsBaselineEstimationRange[fReferenceIndex]);
	  cout << "Reference baseline for v2 calculation will be defined by " << fNMinPointsBaselineEstimationRange[fReferenceIndex] << "lowest points " << endl;
	}
	else if(fVecSystMode[fReferenceIndex]==kLowestPoint){
	  fFitter->SetFixBasetype(-1);
	  cout << "Reference baseline for v2 calculation will be defined by the minimum point " << endl;
	}
	else fFitter->SetFixBasetype(fVecSystMode[fReferenceIndex]);
	if(fVecSystMode[fReferenceIndex]==kTransverse)fFitter->SetBaselineEstimationRange(fMinReferenceBaselineEstimationRange,fMaxReferenceBaselineEstimationRange);

        fFitter->SetFixMeanType(3);
        fFitter->SetFuncType(fFitFuncType[fReferenceIndex]);
        fFitter->Fitting();
        fFitter->DrawLegendWithParameters();
        pave->Draw("same");
        v2subtractedhisto[iPtBin]->Draw("sameep");
        inputhisto[iPtBin]->Draw("samehist");
       
        //draw histos
        
      //  v2subtractedhisto[iPtBin]->Draw("sameep");
      //  fVecHisto[iPtBin]->Draw("samehist");
        v2->Draw("same");
        
        fValuev2NSYield[iPtBin] = fFitter->GetNSYield();
        fValuev2NSSigma[iPtBin] = fFitter->GetNSSigma();
        fValuev2ASYield[iPtBin] = fFitter->GetASYield();
        fValuev2ASSigma[iPtBin] = fFitter->GetASSigma();
        fValuev2Pedestal[iPtBin] = fFitter->GetPedestal();
        
    
        
        fSystValuev2NSYield[iPtBin] = (fValuev2NSYield[iPtBin]-fReferenceHistoNSYield->GetBinContent(iPtBin+1))/fReferenceHistoNSYield->GetBinContent(iPtBin+1);
        fSystValuev2NSSigma[iPtBin] = (fValuev2NSSigma[iPtBin]-fReferenceHistoNSSigma->GetBinContent(iPtBin+1))/fReferenceHistoNSSigma->GetBinContent(iPtBin+1);
        fSystValuev2ASYield[iPtBin] = (fValuev2ASYield[iPtBin]-fReferenceHistoASYield->GetBinContent(iPtBin+1))/fReferenceHistoASYield->GetBinContent(iPtBin+1);
        fSystValuev2ASSigma[iPtBin] = (fValuev2ASSigma[iPtBin]-fReferenceHistoASSigma->GetBinContent(iPtBin+1))/fReferenceHistoASSigma->GetBinContent(iPtBin+1);
        fSystValuev2Pedestal[iPtBin] = (fValuev2Pedestal[iPtBin]-fReferenceHistoPedestal->GetBinContent(iPtBin+1))/fReferenceHistoPedestal->GetBinContent(iPtBin+1);
        
    }

    cv2->SaveAs(Form("Canvas_v2SystemForFitting_AssocPt%1.1ftp%1.1f.png",fAssocTrackPtMin,fAssocTrackPtMax));
    cv2->SaveAs(Form("Canvas_v2SystemForFitting_AssocPt%1.1ftp%1.1f.root",fAssocTrackPtMin,fAssocTrackPtMax));
    
}

//_______________________________________________________________________________
Bool_t AliHFCorrFitSystematics::RunFits(){
    
    //AliHFCorrFitter * fitter = NULL;
    TPaveText * pave = NULL;
    TPaveText * paveRef = NULL;
    
    //TFile * outputfile = NULL;
    Int_t fixbase = -2;
    Int_t fixmean = -2;
    Int_t referenceIndex = -1;
    Int_t index = -1;
    // finish me
    if(fDim<0){
      std::cout << "Systematics fitter not set up! I am doing it for you" << std::endl;
      Bool_t setup = SetUpSystematicsFitter();
        if(!setup) return kFALSE;
    }
    CheckBaselineRanges();
    //fOutputFile->cd();
    
    //*************************************************************************************************************************
    // =================================================== calculating Reference values =======================================
    cout << "fReferenceIndex = " << fReferenceIndex << endl;
    for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){// loop on D meson pt
         
      
      
      fValueSystematicBaselineNSYield[iPtBin] =0;
      fValueSystematicBaselineNSSigma[iPtBin] =0;
      fValueSystematicBaselineASYield[iPtBin] =0;
      fValueSystematicBaselineASSigma[iPtBin] =0;
      fValueSystematicBaselinePedestal[iPtBin] =0;
      
      fValueSystematicNSYieldUp[iPtBin] =0;
      fValueSystematicNSSigmaUp[iPtBin] =0;
      fValueSystematicASYieldUp[iPtBin] =0;
      fValueSystematicASSigmaUp[iPtBin] =0;
      fValueSystematicPedestalUp[iPtBin] =0;
      
      fValueSystematicNSYieldLow[iPtBin] =0;
      fValueSystematicNSSigmaLow[iPtBin] =0;
      fValueSystematicASYieldLow[iPtBin] =0;
      fValueSystematicASSigmaLow[iPtBin] =0;
      fValueSystematicPedestalLow[iPtBin] =0;
      
      paveRef = new TPaveText(0.15,0.6,0.45,0.8,"NDC");
      paveRef->SetName(Form("pave_%d",iPtBin));
      paveRef->AddText(Form("%.1f < p_{T}(D) < %.1f GeV/c",fVecLowEdgeDpt[iPtBin],fVecUpEdgeDpt[iPtBin]));
      if(fAssocTrackPtMax>50) paveRef->AddText(Form("p_{T}(hadron) > %.1f GeV/c",fAssocTrackPtMin));
      else paveRef->AddText(Form("%.1f <p_{T}(hadron) < %.1f GeV/c",fAssocTrackPtMin, fAssocTrackPtMin));
      paveRef->SetFillStyle(0);
      paveRef->SetBorderSize(0);
      paveRef->SetTextSize(0.03);

      Printf("The reference will be calceulated with baseline mode %d and with the fit function mode %d",fVecSystMode[fReferenceIndex],fFitFuncType[fReferenceIndex]);
      fFitter = new AliHFCorrFitter((TH1F*)fVecHisto[iPtBin],fMinFitRange,fMaxFitRange);
      if(fVecSystMode[fReferenceIndex]==kNLowest){
	fFitter->SetFixBasetype(-1*fNMinPointsBaselineEstimationRange[fReferenceIndex]);
	cout << "Reference baseline will be defined by " << fNMinPointsBaselineEstimationRange[fReferenceIndex] << "lowest points " << endl;
      }
      else if(fVecSystMode[fReferenceIndex]==kLowestPoint){
	fFitter->SetFixBasetype(-1);
	cout << "Reference baseline will be defined by the minimum point " << endl;
      }
      else fFitter->SetFixBasetype(fVecSystMode[fReferenceIndex]);
      
      if(fVecSystMode[fReferenceIndex]==kTransverse)fFitter->SetBaselineEstimationRange(fMinReferenceBaselineEstimationRange,fMaxReferenceBaselineEstimationRange);
      fCanvasRefernce->cd(iPtBin+1);
      fFitter->SetFixMeanType(3);
      fFitter->SetFuncType(fFitFuncType[fReferenceIndex]);
      fFitter->Fitting();
      fFitter->DrawLegendWithParameters();
      paveRef->Draw("same");
      if(!fFitter->GetFitFunction()) std::cout << "cannot get fit function to store" << std::endl;
      else fFitFunctions.push_back(fFitter->GetFitFunction());
      // set them
      fReferenceHistoNSYield->SetBinContent(iPtBin+1,fFitter->GetNSYield());
      fReferenceHistoNSSigma->SetBinContent(iPtBin+1,fFitter->GetNSSigma());
      fReferenceHistoASYield->SetBinContent(iPtBin+1,fFitter->GetASYield());
      fReferenceHistoASSigma->SetBinContent(iPtBin+1,fFitter->GetASSigma());
      fReferenceHistoPedestal->SetBinContent(iPtBin+1,fFitter->GetPedestal());
      
      fReferenceHistoNSYield->SetBinError(iPtBin+1,fFitter->GetNSYieldError());
      fReferenceHistoNSSigma->SetBinError(iPtBin+1,fFitter->GetNSSigmaError());
      fReferenceHistoASYield->SetBinError(iPtBin+1,fFitter->GetASYieldError());
      fReferenceHistoASSigma->SetBinError(iPtBin+1,fFitter->GetASSigmaError());
      fReferenceHistoPedestal->SetBinError(iPtBin+1,fFitter->GetPedestalError());
    }
    fOutputFile->cd();
    fCanvasRefernce->Write();
    fReferenceHistoNSYield->Write();
    fReferenceHistoNSSigma->Write();
    fReferenceHistoASYield->Write();
    fReferenceHistoASSigma->Write();
    fReferenceHistoPedestal->Write();
    
    //*************************************************************************************************************************
    // =================================================== looping and fitting on different systematic sources for the baseline  =======================================
    //
    for(Int_t iSystMode = 0; iSystMode < fVecSystModesSize; iSystMode++){ // loop on systematic modes
      
      
        if(fVecSize == 1)fCanvasFitting[iSystMode] = new TCanvas(Form("cFitting_%d",iSystMode),"cFitting",0,0,1000,1000);
        if(fVecSize == 2){fCanvasFitting[iSystMode] = new TCanvas(Form("cFitting_%d",iSystMode),"cFitting",0,0,1500,1000); fCanvasFitting[iSystMode]->Divide(2,1);}
        if(fVecSize == 3){fCanvasFitting[iSystMode] = new TCanvas(Form("cFitting_%d",iSystMode),"cFitting",0,0,1500,1000); fCanvasFitting[iSystMode]->Divide(3,1);}
        if(fVecSize == 4){fCanvasFitting[iSystMode] = new TCanvas(Form("cFitting_%d",iSystMode),"cFitting",0,0,1600,1000); fCanvasFitting[iSystMode]->Divide(2,2);}
        if(fVecSize == 5){fCanvasFitting[iSystMode] = new TCanvas(Form("cFitting_%d",iSystMode),"cFitting",0,0,1600,1000); fCanvasFitting[iSystMode]->Divide(3,2);}
        if(fVecSize == 6){fCanvasFitting[iSystMode] = new TCanvas(Form("cFitting_%d",iSystMode),"cFitting",0,0,1600,1000); fCanvasFitting[iSystMode]->Divide(3,2);}
        if(fVecSize == 7){fCanvasFitting[iSystMode] = new TCanvas(Form("cFitting_%d",iSystMode),"cFitting",0,0,1600,1000); fCanvasFitting[iSystMode]->Divide(4,2);}
        if(fVecSize == 8){fCanvasFitting[iSystMode] = new TCanvas(Form("cFitting_%d",iSystMode),"cFitting",0,0,1600,1000); fCanvasFitting[iSystMode]->Divide(4,2);}
        if(fVecSize == 9){fCanvasFitting[iSystMode] = new TCanvas(Form("cFitting_%d",iSystMode),"cFitting",0,0,1600,1000); fCanvasFitting[iSystMode]->Divide(4,2);}
        
        //fCanvasFitting = (TCanvas*)fCanvasRefernce->Clone(Form("cFitting_%d",iSystMode));
        fCanvasFitting[iSystMode]->SetTitle(Form("Fit canvas %d",iSystMode));
        //cout << "check here = " << Form("ValueHistoNSYield_%d",iSystMode) << endl;
        fValueHistoNSYield->SetName(Form("ValueHistoNSYield_%d",iSystMode));
        fValueHistoNSSigma->SetName(Form("ValueHistoNSSigma_%d",iSystMode));
        fValueHistoASYield->SetName(Form("ValueHistoASYield_%d",iSystMode));
        fValueHistoASSigma->SetName(Form("ValueHistoASSigma_%d",iSystMode));
        fValueHistoPedestal->SetName(Form("ValueHistoPedestal_%d",iSystMode));
        
        
        
        for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
            // loop on D meson pt
            
            pave = new TPaveText(0.15,0.6,0.45,0.8,"NDC");
            pave->SetName(Form("pave_%d",iPtBin));
            pave->AddText(Form("%.1f < p_{T}(D) < %.1f GeV/c",fVecLowEdgeDpt[iPtBin],fVecUpEdgeDpt[iPtBin]));
            if(fAssocTrackPtMax>50) pave->AddText(Form("p_{T}(hadron) > %.1f GeV/c",fAssocTrackPtMin));
            else pave->AddText(Form("%.1f <p_{T}(hadron) < %.1f GeV/c",fAssocTrackPtMin, fAssocTrackPtMin));
            pave->SetFillStyle(0);
            pave->SetBorderSize(0);
            pave->SetTextSize(0.03);
            
            
            if(fVecSystMode[iSystMode]==kMinVar) fFitter = new AliHFCorrFitter((TH1F*)fVecHistoMinVar[iPtBin],fMinFitRange,fMaxFitRange);
            else if(fVecSystMode[iSystMode]==kMaxVar)fFitter = new AliHFCorrFitter((TH1F*)fVecHistoMaxVar[iPtBin],fMinFitRange,fMaxFitRange);
            else fFitter = new AliHFCorrFitter((TH1F*)fVecHisto[iPtBin],fMinFitRange,fMaxFitRange);
               
               
            fCanvasFitting[iSystMode]->cd(iPtBin+1);
            if(fVecSystMode[iSystMode]==kMinVar || fVecSystMode[iSystMode]==kMaxVar) {
	      if(fVecSystMode[fReferenceIndex]==kNLowest){
		fFitter->SetFixBasetype(-1*fNMinPointsBaselineEstimationRange[fReferenceIndex]);
		cout << "Reference baseline for v2 calculation will be defined by " << fNMinPointsBaselineEstimationRange[fReferenceIndex] << "lowest points " << endl;
	      }
	      else if(fVecSystMode[fReferenceIndex]==kLowestPoint){
		fFitter->SetFixBasetype(-1);
		cout << "Reference baseline for v2 calculation will be defined by the minimum point " << endl;
	      }
	      else fFitter->SetFixBasetype(fVecSystMode[fReferenceIndex]);
	      if(fVecSystMode[fReferenceIndex]==kTransverse)fFitter->SetBaselineEstimationRange(fMinReferenceBaselineEstimationRange,fMaxReferenceBaselineEstimationRange);
	    }
	    else if(fVecSystMode[iSystMode]==kNLowest){
	      fFitter->SetFixBasetype(-1*fNMinPointsBaselineEstimationRange[iSystMode]);
	      cout << "Baseline will be defined by " << fNMinPointsBaselineEstimationRange[iSystMode] << "lowest points " << endl;
	    }
	    else if(fVecSystMode[iSystMode]==kLowestPoint){
	      fFitter->SetFixBasetype(-1);
	      cout << "Baseline will be defined by the minimum point " << endl;
	    }
            else fFitter->SetFixBasetype(fVecSystMode[iSystMode]);
	    
            //if(fVecSystMode[iSystMode]==kTransverse) fFitter->SetBaselineEstimationRange(0.25*TMath::Pi(),0.5*TMath::Pi());
              cout << " " << index <<endl;
	      cout << "Baseline range for systmode ("<< iSystMode <<"): " << fMinBaselineEstimationRange[iSystMode] << "," << fMaxBaselineEstimationRange[iSystMode] << endl;
	      cout << "Nmin points for systmode ("<<iSystMode<<"): "<< fNMinPointsBaselineEstimationRange[iSystMode]<< "," <<fNMinPointsBaselineEstimationRange[iSystMode]<<endl;
             cout << " " << endl;
            if(fVecSystMode[iSystMode]==kTransverse)fFitter->SetBaselineEstimationRange(fMinBaselineEstimationRange[iSystMode],fMaxBaselineEstimationRange[iSystMode]);
            if(fVecSystMode[iSystMode]==kTransverseUppStatUnc)fFitter->SetBaselineEstimationRange(fMinBaselineEstimationRange[iSystMode],fMaxBaselineEstimationRange[iSystMode]);
            if(fVecSystMode[iSystMode]==kTransverseLowStatUnc)fFitter->SetBaselineEstimationRange(fMinBaselineEstimationRange[iSystMode],fMaxBaselineEstimationRange[iSystMode]);
	    if(fVecSystMode[iSystMode]==kMinVar || fVecSystMode[iSystMode]==kMaxVar) fFitter->SetBaselineEstimationRange(fMinBaselineEstimationRange[iSystMode],fMaxBaselineEstimationRange[iSystMode]);
            
            
                    //  fFitter->SetBaselineEstimationRange(0.25*TMath::Pi(),0.5*TMath::Pi());
            
            
            fFitter->SetFixMeanType(3);
            fFitter->SetFuncType(fFitFuncType[iSystMode]);
            fFitter->Fitting();
            fFitter->DrawLegendWithParameters();
            pave->Draw("same");
            
            fValueHistoNSYield->SetBinContent(iPtBin+1,fFitter->GetNSYield());
            fValueHistoNSSigma->SetBinContent(iPtBin+1,fFitter->GetNSSigma());
            fValueHistoASYield->SetBinContent(iPtBin+1,fFitter->GetASYield());
            fValueHistoASSigma->SetBinContent(iPtBin+1,fFitter->GetASSigma());
            fValueHistoPedestal->SetBinContent(iPtBin+1,fFitter->GetPedestal());
            
            fValueHistoNSYield->SetBinError(iPtBin+1,fFitter->GetNSYieldError());
            fValueHistoNSSigma->SetBinError(iPtBin+1,fFitter->GetNSSigmaError());
            fValueHistoASYield->SetBinError(iPtBin+1,fFitter->GetASYieldError());
            fValueHistoASSigma->SetBinError(iPtBin+1,fFitter->GetASSigmaError());
            fValueHistoPedestal->SetBinError(iPtBin+1,fFitter->GetPedestalError());
            
            
            if(!ComputeRatios(fReferenceHistoNSYield, fValueHistoNSYield, fRatioHistoNSYield)) return kFALSE;
            if(!ComputeRatios(fReferenceHistoNSSigma, fValueHistoNSSigma, fRatioHistoNSSigma)) return kFALSE;
            if(!ComputeRatios(fReferenceHistoASYield, fValueHistoASYield, fRatioHistoASYield)) return kFALSE;
            if(!ComputeRatios(fReferenceHistoASSigma, fValueHistoASSigma, fRatioHistoASSigma)) return kFALSE;
            if(!ComputeRatios(fReferenceHistoPedestal, fValueHistoPedestal, fRatioHistoPedestal)) return kFALSE;
            
        }// end loop on D meson pt
        
            //Int_t index = GetBinIndex(Int_t ptbinindex, Int_t systematicIndex)
      
        
      /*  fSystematicSourcesNSYield->SetName(Form("SystematicSourcesNSYield%d",iSystMode));
        fSystematicSourcesNSSigma->SetName(Form("SystematicSourcesNSSigma%d",iSystMode));
        fSystematicSourcesASYield->SetName(Form("SystematicSourcesASYield%d",iSystMode));
        fSystematicSourcesASSigma->SetName(Form("SystematicSourcesASSigma%d",iSystMode));
        fSystematicSourcesPedestal->SetName(Form("SystematicSourcesPedestal%d",iSystMode));
        
        fSystematicSourcesNSYield->SetTitle(Form("SystematicSourcesNSYield%d; D p_{T} GeV/c",iSystMode));
        fSystematicSourcesNSSigma->SetTitle(Form("SystematicSourcesNSSigma%d; D p_{T} GeV/c",iSystMode));
        fSystematicSourcesASYield->SetTitle(Form("SystematicSourcesASYield%d; D p_{T} GeV/c",iSystMode));
        fSystematicSourcesASSigma->SetTitle(Form("SystematicSourcesASSigma%d; D p_{T} GeV/c",iSystMode));
        fSystematicSourcesPedestal->SetTitle(Form("SystematicSourcesPedestal%d; D p_{T} GeV/c",iSystMode));
        
        fSystematicSourcesNSYield->GetYaxis()->SetRangeUser(0,3);
        fSystematicSourcesNSSigma->GetYaxis()->SetRangeUser(0,3);
        fSystematicSourcesASYield->GetYaxis()->SetRangeUser(0,3);
        fSystematicSourcesASSigma->GetYaxis()->SetRangeUser(0,3);
        fSystematicSourcesPedestal->GetYaxis()->SetRangeUser(0,3);*/
        
            for(Int_t iPtBin =0;iPtBin<fVecSize; iPtBin++){
            
                index = GetBinIndex(iPtBin, iSystMode);
               //std:: cout << "index = " << index <<  ";iPtBin = " << iPtBin <<  ";iSystMode = " << iSystMode << std::endl;
                if(fVecSystMode[iSystMode] == kBinCount){
                    fValueNSYield[index] = fValueHistoNSYield->GetBinContent(iPtBin+1);
                    fRatioNSYield[index] = TMath::Abs(fRatioHistoNSYield->GetBinContent(iPtBin+1)-1);
                    fValueNSSigma[index] = 0;
                    fRatioNSSigma[index] = 0;
                    fValueASYield[index] = 0;
                    fRatioASYield[index] = 0;
                    fValueASSigma[index] = 0;
                    fRatioASSigma[index] = 0;
                    fValuePedestal[index] = 0;
                    fRatioPedestal[index] = 0;
                }
                else{
                    fValueNSYield[index] = fValueHistoNSYield->GetBinContent(iPtBin+1);
                    fRatioNSYield[index] = TMath::Abs(fRatioHistoNSYield->GetBinContent(iPtBin+1)-1);
                    fValueNSSigma[index] = fValueHistoNSSigma->GetBinContent(iPtBin+1);
                    fRatioNSSigma[index] = TMath::Abs(fRatioHistoNSSigma->GetBinContent(iPtBin+1)-1);
                    fValueASYield[index] = fValueHistoASYield->GetBinContent(iPtBin+1);
                    fRatioASYield[index] = TMath::Abs(fRatioHistoASYield->GetBinContent(iPtBin+1)-1);
                    fValueASSigma[index] = fValueHistoASSigma->GetBinContent(iPtBin+1);
                    fRatioASSigma[index] = TMath::Abs(fRatioHistoASSigma->GetBinContent(iPtBin+1)-1);
                    fValuePedestal[index] = fValueHistoPedestal->GetBinContent(iPtBin+1);
                    fRatioPedestal[index] = TMath::Abs(fRatioHistoPedestal->GetBinContent(iPtBin+1)-1);
                }
                
               // cout << "Bin center = " << fValueHistoNSYield->GetBinCenter(iPtBin+1) << endl;
                
                std::cout << "------------------------------------" << std::endl;
                std::cout << "index = " << index <<  ";iPtBin = " << iPtBin <<  ";iSystMode = " << iSystMode << std::endl;
                std::cout << "Estimated Yield uncertainty = " << fRatioNSYield[index] << std::endl;
                // plot systematics
                
                /*
                fSystematicSourcesNSYield->SetBinContent(iPtBin+1,1);
                fSystematicSourcesNSSigma->SetBinContent(iPtBin+1,1);
                fSystematicSourcesASYield->SetBinContent(iPtBin+1,1);
                fSystematicSourcesASSigma->SetBinContent(iPtBin+1,1);
                fSystematicSourcesPedestal->SetBinContent(iPtBin+1,1);
                
                fSystematicSourcesNSYield->SetBinError(iPtBin+1,fRatioNSYield[index]);
                fSystematicSourcesNSSigma->SetBinError(iPtBin+1,fRatioNSSigma[index]);
                fSystematicSourcesASYield->SetBinError(iPtBin+1,fRatioASYield[index]);
                fSystematicSourcesASSigma->SetBinError(iPtBin+1,fRatioASSigma[index]);
                fSystematicSourcesPedestal->SetBinError(iPtBin+1,fRatioPedestal[index]);
                */
               
                
                
                // combine systematics
                if(fCombineSystematics == kSumQuadr){
                    
                
                    
                    fValueSystematicBaselineNSYield[iPtBin] += fRatioNSYield[index]*fRatioNSYield[index];
                    fValueSystematicBaselineNSSigma[iPtBin] += fRatioNSSigma[index]*fRatioNSSigma[index];
                    fValueSystematicBaselineASYield[iPtBin] += fRatioASYield[index]*fRatioASYield[index];
                    fValueSystematicBaselineASSigma[iPtBin] += fRatioASSigma[index]*fRatioASSigma[index];
                    fValueSystematicBaselinePedestal[iPtBin] += fRatioPedestal[index]*fRatioPedestal[index];
   
                    
                }
                if(fCombineSystematics == kMax){
                    // check me
                    if(fValueSystematicBaselineNSYield[iPtBin]< fRatioNSYield[index]) fValueSystematicBaselineNSYield[iPtBin] = fRatioNSYield[index];
                     if(fValueSystematicBaselineNSSigma[iPtBin]< fRatioNSSigma[index]) fValueSystematicBaselineNSSigma[iPtBin] = fRatioNSSigma[index];
                     if(fValueSystematicBaselineASYield[iPtBin]< fRatioASYield[index]) fValueSystematicBaselineASYield[iPtBin] = fRatioASYield[index];
                     if(fValueSystematicBaselineASSigma[iPtBin]< fRatioASSigma[index]) fValueSystematicBaselineASSigma[iPtBin] = fRatioASSigma[index];
                     if(fValueSystematicBaselinePedestal[iPtBin]< fRatioPedestal[index]) fValueSystematicBaselinePedestal[iPtBin] = fRatioPedestal[index];
                    
          
                    
                }
                 if(fCombineSystematics == kRMS){
                     fRMSHistoNSYield[iPtBin]->Fill(fValueNSYield[index]);
                     fRMSHistoNSSigma[iPtBin]->Fill(fValueNSSigma[index]);
                     fRMSHistoASYield[iPtBin]->Fill(fValueASYield[index]);
                     fRMSHistoASSigma[iPtBin]->Fill(fValueASSigma[index]);
                     fRMSHistoPedestal[iPtBin]->Fill(fValuePedestal[index]);



                 }
                 if(fCombineSystematics == kEnvelope_RMS_BaselStat){
                    if(fVecSystMode[iSystMode] == kTransverseUppStatUnc) {
                      fValueSystematicBaseline_FromBaselStatUp_NSYield[iPtBin] = fRatioNSYield[index];
                      fValueSystematicBaseline_FromBaselStatUp_NSSigma[iPtBin] = fRatioNSSigma[index];
                      fValueSystematicBaseline_FromBaselStatUp_ASYield[iPtBin] = fRatioASYield[index];
                      fValueSystematicBaseline_FromBaselStatUp_ASSigma[iPtBin] = fRatioASSigma[index];
                      fValueSystematicBaseline_FromBaselStatUp_Pedestal[iPtBin] = fRatioPedestal[index];
                    }
                    else if(fVecSystMode[iSystMode] == kTransverseLowStatUnc) {
                      fValueSystematicBaseline_FromBaselStatLo_NSYield[iPtBin] = fRatioNSYield[index];
                      fValueSystematicBaseline_FromBaselStatLo_NSSigma[iPtBin] = fRatioNSSigma[index];
                      fValueSystematicBaseline_FromBaselStatLo_ASYield[iPtBin] = fRatioASYield[index];
                      fValueSystematicBaseline_FromBaselStatLo_ASSigma[iPtBin] = fRatioASSigma[index];
                      fValueSystematicBaseline_FromBaselStatLo_Pedestal[iPtBin] = fRatioPedestal[index];
                    }
                    else {
                       fRMSHistoNSYield[iPtBin]->Fill(fValueNSYield[index]);
                       fRMSHistoNSSigma[iPtBin]->Fill(fValueNSSigma[index]);
                       fRMSHistoASYield[iPtBin]->Fill(fValueASYield[index]);
                       fRMSHistoASSigma[iPtBin]->Fill(fValueASSigma[index]);
                       fRMSHistoPedestal[iPtBin]->Fill(fValuePedestal[index]);                      
                    }
                 }                 
                
                
        } // end 2nd loop on D meson pt
        
        
        /*
        fSystematicSourcesNSYield->SetFillColor(iSystMode);
        fSystematicSourcesNSSigma->SetFillColor(iSystMode);
        fSystematicSourcesASYield->SetFillColor(iSystMode);
        fSystematicSourcesASSigma->SetFillColor(iSystMode);
        fSystematicSourcesPedestal->SetFillColor(iSystMode);
        
        fSystematicSourcesNSYield->SetFillStyle(3000+iSystMode);
        fSystematicSourcesNSSigma->SetFillStyle(3000+iSystMode);
        fSystematicSourcesASYield->SetFillStyle(3000+iSystMode);
        fSystematicSourcesASSigma->SetFillStyle(3000+iSystMode);
        fSystematicSourcesPedestal->SetFillStyle(3000+iSystMode);
        
        fCanvasSystematicSourcesNSYield->cd();
        if(!iSystMode)fSystematicSourcesNSYield->Draw("e2");
        fSystematicSourcesNSYield->Draw("e2same");
        
        fCanvasSystematicSourcesNSSigma->cd();
        if(!iSystMode)fSystematicSourcesNSSigma->Draw("e2");
        fSystematicSourcesNSSigma->Draw("e2same");
        
        fCanvasSystematicSourcesASYield->cd();
        if(!iSystMode)fSystematicSourcesASYield->Draw("e2");
        fSystematicSourcesASYield->Draw("e2same");
        
        fCanvasSystematicSourcesASSigma->cd();
        if(!iSystMode)fSystematicSourcesASSigma->Draw("e2");
        fSystematicSourcesASSigma->Draw("e2same");
        
        fCanvasSystematicSourcesPedestal->cd();
        if(!iSystMode)fSystematicSourcesPedestal->Draw("e2");
        fSystematicSourcesPedestal->Draw("e2same");
        */
        
       // outputfile = new TFile(fOutputFileName,"UPDATE");
     
        //outputfile = TFile::Open(fOutputFileName.Data());
        
        
        fCanvasFitting[iSystMode]->Write();
        fValueHistoNSYield->Write();
        fValueHistoNSSigma->Write();
        fValueHistoASYield->Write();
        fValueHistoASSigma->Write();
        fValueHistoPedestal->Write();
        
        fRatioHistoNSYield->Write();
        fRatioHistoNSSigma->Write();
        fRatioHistoASYield->Write();
        fRatioHistoASSigma->Write();
        fRatioHistoPedestal->Write();
      
        
       
     } // end loop on systematic modes
    
    
   
    
    //*************************************************************************************************************************
    // =================================================== Estimating total systematic from baseline determination =======================================
    // loop on pt bins
    for(Int_t iPtBin=0; iPtBin<fVecSize; iPtBin++){
        

        
        if(fCombineSystematics == kSumQuadr){
            fValueSystematicBaselineNSYield[iPtBin] = TMath::Sqrt(fValueSystematicBaselineNSYield[iPtBin]);
            fValueSystematicBaselineNSSigma[iPtBin] = TMath::Sqrt(fValueSystematicBaselineNSSigma[iPtBin]);
            fValueSystematicBaselineASYield[iPtBin] = TMath::Sqrt(fValueSystematicBaselineASYield[iPtBin]);
            fValueSystematicBaselineASSigma[iPtBin] = TMath::Sqrt(fValueSystematicBaselineASSigma[iPtBin]);
            fValueSystematicBaselinePedestal[iPtBin] = TMath::Sqrt(fValueSystematicBaselinePedestal[iPtBin]);
            

        }
        if(fCombineSystematics == kRMS){
	  std::cout << "Using RMS" << std::endl;
	  if(fReferenceHistoNSYield->GetBinContent(iPtBin+1)>1.e-6){
	    fValueSystematicBaselineNSYield[iPtBin] = fRMSHistoNSYield[iPtBin]->GetRMS()/fReferenceHistoNSYield->GetBinContent(iPtBin+1);
	    fValueSystematicBaselineNSSigma[iPtBin] = fRMSHistoNSSigma[iPtBin]->GetRMS()/fReferenceHistoNSSigma->GetBinContent(iPtBin+1);
	    fValueSystematicBaselineASYield[iPtBin] = fRMSHistoASYield[iPtBin]->GetRMS()/fReferenceHistoASYield->GetBinContent(iPtBin+1);
	    fValueSystematicBaselineASSigma[iPtBin] = fRMSHistoASSigma[iPtBin]->GetRMS()/fReferenceHistoASSigma->GetBinContent(iPtBin+1);
	    fValueSystematicBaselinePedestal[iPtBin] = fRMSHistoPedestal[iPtBin]->GetRMS()/fReferenceHistoPedestal->GetBinContent(iPtBin+1);
	  }
	  fRMSHistoNSYield[iPtBin]->Write();
	  fRMSHistoNSSigma[iPtBin]->Write();
	  fRMSHistoASYield[iPtBin]->Write();
	  fRMSHistoASSigma[iPtBin]->Write();
	  fRMSHistoPedestal[iPtBin]->Write();
	  
        }

        if(fCombineSystematics == kEnvelope_RMS_BaselStat){
          std::cout << "Using RMS" << std::endl;
          if(fReferenceHistoNSYield->GetBinContent(iPtBin+1)>1.e-6){
            fValueSystematicBaselineNSYield[iPtBin] = fRMSHistoNSYield[iPtBin]->GetRMS()/fReferenceHistoNSYield->GetBinContent(iPtBin+1);
            fValueSystematicBaselineNSSigma[iPtBin] = fRMSHistoNSSigma[iPtBin]->GetRMS()/fReferenceHistoNSSigma->GetBinContent(iPtBin+1);
            fValueSystematicBaselineASYield[iPtBin] = fRMSHistoASYield[iPtBin]->GetRMS()/fReferenceHistoASYield->GetBinContent(iPtBin+1);
            fValueSystematicBaselineASSigma[iPtBin] = fRMSHistoASSigma[iPtBin]->GetRMS()/fReferenceHistoASSigma->GetBinContent(iPtBin+1);
            fValueSystematicBaselinePedestal[iPtBin] = fRMSHistoPedestal[iPtBin]->GetRMS()/fReferenceHistoPedestal->GetBinContent(iPtBin+1);
            fRMSRelative_NSYield[iPtBin] = fValueSystematicBaselineNSYield[iPtBin];
            fRMSRelative_NSSigma[iPtBin] = fValueSystematicBaselineNSSigma[iPtBin];
            fRMSRelative_ASYield[iPtBin] = fValueSystematicBaselineASYield[iPtBin];
            fRMSRelative_ASSigma[iPtBin] = fValueSystematicBaselineASSigma[iPtBin];
            fRMSRelative_Pedestal[iPtBin] = fValueSystematicBaselinePedestal[iPtBin];
          }
          //now substitute RMS with lowstat/uppstat variations, if are larger
          if(fValueSystematicBaselineNSYield[iPtBin] < TMath::Abs(fValueSystematicBaseline_FromBaselStatUp_NSYield[iPtBin])) fValueSystematicBaselineNSYield[iPtBin] = TMath::Abs(fValueSystematicBaseline_FromBaselStatUp_NSYield[iPtBin]);
          if(fValueSystematicBaselineNSSigma[iPtBin] < TMath::Abs(fValueSystematicBaseline_FromBaselStatUp_NSSigma[iPtBin])) fValueSystematicBaselineNSSigma[iPtBin] = TMath::Abs(fValueSystematicBaseline_FromBaselStatUp_NSSigma[iPtBin]);
          if(fValueSystematicBaselineASYield[iPtBin] < TMath::Abs(fValueSystematicBaseline_FromBaselStatUp_ASYield[iPtBin])) fValueSystematicBaselineASYield[iPtBin] = TMath::Abs(fValueSystematicBaseline_FromBaselStatUp_ASYield[iPtBin]);
          if(fValueSystematicBaselineASSigma[iPtBin] < TMath::Abs(fValueSystematicBaseline_FromBaselStatUp_ASSigma[iPtBin])) fValueSystematicBaselineASSigma[iPtBin] = TMath::Abs(fValueSystematicBaseline_FromBaselStatUp_ASSigma[iPtBin]);
          if(fValueSystematicBaselinePedestal[iPtBin] < TMath::Abs(fValueSystematicBaseline_FromBaselStatUp_Pedestal[iPtBin])) fValueSystematicBaselinePedestal[iPtBin] = TMath::Abs(fValueSystematicBaseline_FromBaselStatUp_Pedestal[iPtBin]);

          if(fValueSystematicBaselineNSYield[iPtBin] < TMath::Abs(fValueSystematicBaseline_FromBaselStatLo_NSYield[iPtBin])) fValueSystematicBaselineNSYield[iPtBin] = TMath::Abs(fValueSystematicBaseline_FromBaselStatLo_NSYield[iPtBin]);
          if(fValueSystematicBaselineNSSigma[iPtBin] < TMath::Abs(fValueSystematicBaseline_FromBaselStatLo_NSSigma[iPtBin])) fValueSystematicBaselineNSSigma[iPtBin] = TMath::Abs(fValueSystematicBaseline_FromBaselStatLo_NSSigma[iPtBin]);
          if(fValueSystematicBaselineASYield[iPtBin] < TMath::Abs(fValueSystematicBaseline_FromBaselStatLo_ASYield[iPtBin])) fValueSystematicBaselineASYield[iPtBin] = TMath::Abs(fValueSystematicBaseline_FromBaselStatLo_ASYield[iPtBin]);
          if(fValueSystematicBaselineASSigma[iPtBin] < TMath::Abs(fValueSystematicBaseline_FromBaselStatLo_ASSigma[iPtBin])) fValueSystematicBaselineASSigma[iPtBin] = TMath::Abs(fValueSystematicBaseline_FromBaselStatLo_ASSigma[iPtBin]);
          if(fValueSystematicBaselinePedestal[iPtBin] < TMath::Abs(fValueSystematicBaseline_FromBaselStatLo_Pedestal[iPtBin])) fValueSystematicBaselinePedestal[iPtBin] = TMath::Abs(fValueSystematicBaseline_FromBaselStatLo_Pedestal[iPtBin]);

          fRMSHistoNSYield[iPtBin]->Write();
          fRMSHistoNSSigma[iPtBin]->Write();
          fRMSHistoASYield[iPtBin]->Write();
          fRMSHistoASSigma[iPtBin]->Write();
          fRMSHistoPedestal[iPtBin]->Write();
          
        }        
        
    }
    
    
    //*************************************************************************************************************************
    // =================================================== Plotting total systematic from baseline determination =======================================
    

    
        TLegend * SystematicsLegendNSYield = new TLegend(0.5,0.6,0.85,0.8);
        SystematicsLegendNSYield->SetFillColor(0);
        SystematicsLegendNSYield->SetTextSize(0.025);
        SystematicsLegendNSYield->SetBorderSize(0);
        fCanvasSystematicSourcesNSYield->cd();
        for(Int_t iSystMode = 0;iSystMode<fVecSystModesSize; iSystMode++ ){
            DrawBaselineSystematicsOnCanvas(fSystematicSourcesNSYield,iSystMode,fRatioNSYield,SystematicsLegendNSYield);
        }
        if(fCombineSystematics == kEnvelope_RMS_BaselStat) DrawRMSBaselineSystematicOnCanvas(fSystematicSourcesNSYield,fRMSRelative_NSYield,fRMSRelative_NSYield,SystematicsLegendNSYield,kTRUE);
        DrawTotalBaselineSystematicOnCanvas(fSystematicSourcesNSYield,fValueSystematicBaselineNSYield, fValueSystematicBaselineNSYield ,SystematicsLegendNSYield,kTRUE);
        SystematicsLegendNSYield->Draw("same");
       DefinePaveText();
        
        TLegend * SystematicsLegendNSSigma = new TLegend(0.5,0.6,0.85,0.8);
        SystematicsLegendNSSigma->SetFillColor(0);
        SystematicsLegendNSSigma->SetTextSize(0.025);
        SystematicsLegendNSSigma->SetBorderSize(0);
        fCanvasSystematicSourcesNSSigma->cd();
        for(Int_t iSystMode = 0;iSystMode<fVecSystModesSize; iSystMode++ ){
            DrawBaselineSystematicsOnCanvas(fSystematicSourcesNSSigma,iSystMode,fRatioNSSigma,SystematicsLegendNSSigma);
        }
        if(fCombineSystematics == kEnvelope_RMS_BaselStat) DrawRMSBaselineSystematicOnCanvas(fSystematicSourcesNSSigma,fRMSRelative_NSSigma,fRMSRelative_NSSigma,SystematicsLegendNSSigma,kTRUE);
        DrawTotalBaselineSystematicOnCanvas(fSystematicSourcesNSSigma,fValueSystematicBaselineNSSigma, fValueSystematicBaselineNSSigma ,SystematicsLegendNSSigma,kTRUE);
        SystematicsLegendNSSigma->Draw("same");
        DefinePaveText();
        
        TLegend * SystematicsLegendASYield = new TLegend(0.5,0.6,0.85,0.8);
        SystematicsLegendASYield->SetFillColor(0);
        SystematicsLegendASYield->SetTextSize(0.025);
        SystematicsLegendASYield->SetBorderSize(0);
        fCanvasSystematicSourcesASYield->cd();
        for(Int_t iSystMode = 0;iSystMode<fVecSystModesSize; iSystMode++ ){
            DrawBaselineSystematicsOnCanvas(fSystematicSourcesASYield,iSystMode,fRatioASYield,SystematicsLegendASYield);
        }
        if(fCombineSystematics == kEnvelope_RMS_BaselStat) DrawRMSBaselineSystematicOnCanvas(fSystematicSourcesASYield,fRMSRelative_ASYield,fRMSRelative_ASYield,SystematicsLegendASYield,kTRUE);
        DrawTotalBaselineSystematicOnCanvas(fSystematicSourcesASYield,fValueSystematicBaselineASYield, fValueSystematicBaselineASYield ,SystematicsLegendASYield,kTRUE);
        SystematicsLegendASYield->Draw("same");
        DefinePaveText();
        
        TLegend * SystematicsLegendASSigma = new TLegend(0.5,0.6,0.85,0.8);
        SystematicsLegendASSigma->SetFillColor(0);
        SystematicsLegendASSigma->SetTextSize(0.025);
        SystematicsLegendASSigma->SetBorderSize(0);
        
        fCanvasSystematicSourcesASSigma->cd();
        for(Int_t iSystMode = 0;iSystMode<fVecSystModesSize; iSystMode++ ){
            DrawBaselineSystematicsOnCanvas(fSystematicSourcesASSigma,iSystMode,fRatioASSigma,SystematicsLegendASSigma);
        }
        if(fCombineSystematics == kEnvelope_RMS_BaselStat) DrawRMSBaselineSystematicOnCanvas(fSystematicSourcesASSigma,fRMSRelative_ASSigma,fRMSRelative_ASSigma,SystematicsLegendASSigma,kTRUE);
        DrawTotalBaselineSystematicOnCanvas(fSystematicSourcesASSigma,fValueSystematicBaselineASSigma, fValueSystematicBaselineASSigma ,SystematicsLegendASSigma,kTRUE);
        SystematicsLegendASSigma->Draw("same");
    DefinePaveText();
    
        TLegend * SystematicsLegendPedestal = new TLegend(0.5,0.6,0.85,0.8);
        SystematicsLegendPedestal->SetFillColor(0);
        SystematicsLegendPedestal->SetTextSize(0.025);
        SystematicsLegendPedestal->SetBorderSize(0);
        fCanvasSystematicSourcesPedestal->cd();
        for(Int_t iSystMode = 0;iSystMode<fVecSystModesSize; iSystMode++ ){
            DrawBaselineSystematicsOnCanvas(fSystematicSourcesPedestal,iSystMode,fRatioPedestal,SystematicsLegendPedestal);
        }
        if(fCombineSystematics == kEnvelope_RMS_BaselStat) DrawRMSBaselineSystematicOnCanvas(fSystematicSourcesPedestal,fRMSRelative_Pedestal,fRMSRelative_Pedestal,SystematicsLegendPedestal,kTRUE);    
        DrawTotalBaselineSystematicOnCanvas(fSystematicSourcesPedestal,fValueSystematicBaselinePedestal, fValueSystematicBaselinePedestal ,SystematicsLegendPedestal,kTRUE);
        SystematicsLegendPedestal->Draw("same");
        DefinePaveText();
    
    cout << "about to fit v2 systematics" << endl;
    
    //*************************************************************************************************************************
    // =================================================== Computing v2 systematics =======================================
      if(fIspPb)  Fitv2Systematics();
    
    cout <<"drawing ns yield final " << endl;
    TLegend * TotalSystematicsLegendNSyield = new TLegend(0.5,0.6,0.85,0.8);
    TotalSystematicsLegendNSyield->SetFillColor(0);
    TotalSystematicsLegendNSyield->SetTextSize(0.02);
    TotalSystematicsLegendNSyield->SetBorderSize(0);
    fCanvasTotalSystematicSourcesNSYield->cd();
    DrawTotalSystematicsOnCanvas("NSYield",TotalSystematicsLegendNSyield);
    DefinePaveText();
    
    cout <<"drawing ns yield sigma " << endl;
    TLegend * TotalSystematicsLegendNSSigma = new TLegend(0.5,0.6,0.85,0.8);
    TotalSystematicsLegendNSSigma->SetFillColor(0);
    TotalSystematicsLegendNSSigma->SetTextSize(0.02);
    TotalSystematicsLegendNSSigma->SetBorderSize(0);
    fCanvasTotalSystematicSourcesNSSigma->cd();
    DrawTotalSystematicsOnCanvas("NSSigma",TotalSystematicsLegendNSSigma);
    DefinePaveText();
    
    cout <<"drawing ss yield final " << endl;
    TLegend * TotalSystematicsLegendASyield = new TLegend(0.5,0.6,0.85,0.8);
    TotalSystematicsLegendASyield->SetFillColor(0);
    TotalSystematicsLegendASyield->SetTextSize(0.02);
    TotalSystematicsLegendASyield->SetBorderSize(0);
    fCanvasTotalSystematicSourcesASYield->cd();
    DrawTotalSystematicsOnCanvas("ASYield",TotalSystematicsLegendASyield);
    DefinePaveText();
    
    cout <<"drawing as sigma final " << endl;
    TLegend * TotalSystematicsLegendASSigma = new TLegend(0.5,0.6,0.85,0.8);
    TotalSystematicsLegendASSigma->SetFillColor(0);
    TotalSystematicsLegendASSigma->SetTextSize(0.02);
    TotalSystematicsLegendASSigma->SetBorderSize(0);
    fCanvasTotalSystematicSourcesASSigma->cd();
    DrawTotalSystematicsOnCanvas("ASSigma",TotalSystematicsLegendASSigma);
    DefinePaveText();
    
    cout <<"drawing ns yield final " << endl;
    TLegend * TotalSystematicsLegendPedestal = new TLegend(0.5,0.6,0.85,0.8);
    TotalSystematicsLegendPedestal->SetFillColor(0);
    TotalSystematicsLegendPedestal->SetTextSize(0.02);
    TotalSystematicsLegendPedestal->SetBorderSize(0);
    fCanvasTotalSystematicSourcesPedestal->cd();
    DrawTotalSystematicsOnCanvas("Pedestal",TotalSystematicsLegendPedestal);
    DefinePaveText();
    
    //*************************************************************************************************************************
    // =================================================== Computing final systematics=======================================
    
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++ add correlated systematics + v2
        for(Int_t iPtBin=0; iPtBin<fVecSize; iPtBin++){
            
          //  compute v2 systematics
            
      
            
	  // summing the correalted syt
        fValueSystematicNSYieldUp[iPtBin] = fValueSystematicBaselineNSYield[iPtBin]*fValueSystematicBaselineNSYield[iPtBin];
        fValueSystematicNSSigmaUp[iPtBin] = fValueSystematicBaselineNSSigma[iPtBin]*fValueSystematicBaselineNSSigma[iPtBin];
        fValueSystematicASYieldUp[iPtBin] = fValueSystematicBaselineASYield[iPtBin]*fValueSystematicBaselineASYield[iPtBin];
        fValueSystematicASSigmaUp[iPtBin] = fValueSystematicBaselineASSigma[iPtBin]*fValueSystematicBaselineASSigma[iPtBin];
        fValueSystematicPedestalUp[iPtBin] = fValueSystematicBaselinePedestal[iPtBin]*fValueSystematicBaselinePedestal[iPtBin];
        
        fValueSystematicNSYieldLow[iPtBin] = fValueSystematicBaselineNSYield[iPtBin]*fValueSystematicBaselineNSYield[iPtBin];
        fValueSystematicNSSigmaLow[iPtBin]= fValueSystematicBaselineNSSigma[iPtBin]*fValueSystematicBaselineNSSigma[iPtBin];
        fValueSystematicASYieldLow[iPtBin] = fValueSystematicBaselineASYield[iPtBin]*fValueSystematicBaselineASYield[iPtBin];
        fValueSystematicASSigmaLow[iPtBin] = fValueSystematicBaselineASSigma[iPtBin]*fValueSystematicBaselineASSigma[iPtBin];
        fValueSystematicPedestalLow[iPtBin] = fValueSystematicBaselinePedestal[iPtBin]*fValueSystematicBaselinePedestal[iPtBin];
        
        
        
        if(fUseCorrelatedSystematics){
            fValueSystematicNSYieldUp[iPtBin] += fVecMaxCorrelatedSyst[iPtBin]*fVecMaxCorrelatedSyst[iPtBin];
            if(fUseCorrelatedSystematicsForWidths)fValueSystematicNSSigmaUp[iPtBin] += fVecMaxCorrelatedSyst[iPtBin]*fVecMaxCorrelatedSyst[iPtBin];
            fValueSystematicASYieldUp[iPtBin] += fVecMaxCorrelatedSyst[iPtBin]*fVecMaxCorrelatedSyst[iPtBin];
            if(fUseCorrelatedSystematicsForWidths)fValueSystematicASSigmaUp[iPtBin] += fVecMaxCorrelatedSyst[iPtBin]*fVecMaxCorrelatedSyst[iPtBin];
            fValueSystematicPedestalUp[iPtBin] += fVecMaxCorrelatedSyst[iPtBin]*fVecMaxCorrelatedSyst[iPtBin];
        
            fValueSystematicNSYieldLow[iPtBin] += fVecMinCorrelatedSyst[iPtBin]*fVecMinCorrelatedSyst[iPtBin];
            if(fUseCorrelatedSystematicsForWidths)fValueSystematicNSSigmaLow[iPtBin] += fVecMinCorrelatedSyst[iPtBin]*fVecMinCorrelatedSyst[iPtBin];
            fValueSystematicASYieldLow[iPtBin] += fVecMinCorrelatedSyst[iPtBin]*fVecMinCorrelatedSyst[iPtBin];
            if(fUseCorrelatedSystematicsForWidths)fValueSystematicASSigmaLow[iPtBin] += fVecMinCorrelatedSyst[iPtBin]*fVecMinCorrelatedSyst[iPtBin];
            fValueSystematicPedestalLow[iPtBin] += fVecMinCorrelatedSyst[iPtBin]*fVecMinCorrelatedSyst[iPtBin];
        }
            
            if(fIspPb){
               if(fSystValuev2NSYield[iPtBin]>=0) fValueSystematicNSYieldUp[iPtBin] += fSystValuev2NSYield[iPtBin]*fSystValuev2NSYield[iPtBin];
               if(fSystValuev2NSSigma[iPtBin]>=0)  fValueSystematicNSSigmaUp[iPtBin] += fSystValuev2NSSigma[iPtBin]*fSystValuev2NSSigma[iPtBin];
               if(fSystValuev2ASYield[iPtBin]>=0)  fValueSystematicASYieldUp[iPtBin] += fSystValuev2ASYield[iPtBin]*fSystValuev2ASYield[iPtBin];
               if(fSystValuev2ASSigma[iPtBin]>=0)  fValueSystematicASSigmaUp[iPtBin] += fSystValuev2ASSigma[iPtBin]*fSystValuev2ASSigma[iPtBin];
               if(fSystValuev2Pedestal[iPtBin]>=0)  fValueSystematicPedestalUp[iPtBin] += fSystValuev2Pedestal[iPtBin]*fSystValuev2Pedestal[iPtBin];
                
               if(fSystValuev2NSYield[iPtBin]<0) fValueSystematicNSYieldLow[iPtBin] +=  fSystValuev2NSYield[iPtBin]*fSystValuev2NSYield[iPtBin];
               if(fSystValuev2NSSigma[iPtBin]<0) fValueSystematicNSSigmaLow[iPtBin] += fSystValuev2NSSigma[iPtBin]*fSystValuev2NSSigma[iPtBin];
               if(fSystValuev2ASYield[iPtBin]<0) fValueSystematicASYieldLow[iPtBin] += fSystValuev2ASYield[iPtBin]*fSystValuev2ASYield[iPtBin];
               if(fSystValuev2ASSigma[iPtBin]<0) fValueSystematicASSigmaLow[iPtBin] += fSystValuev2ASSigma[iPtBin]*fSystValuev2ASSigma[iPtBin];
               if(fSystValuev2Pedestal[iPtBin]<0) fValueSystematicPedestalLow[iPtBin] += fSystValuev2Pedestal[iPtBin]*fSystValuev2Pedestal[iPtBin];
            }
            
        //cout << "Input yield = " << fValueSystematicNSYieldUp[iPtBin] << endl;
        fValueSystematicNSYieldUp[iPtBin] = TMath::Sqrt(fValueSystematicNSYieldUp[iPtBin]);
        fValueSystematicNSSigmaUp[iPtBin] = TMath::Sqrt(fValueSystematicNSSigmaUp[iPtBin]);
        fValueSystematicASYieldUp[iPtBin] = TMath::Sqrt(fValueSystematicASYieldUp[iPtBin]);
        fValueSystematicASSigmaUp[iPtBin] = TMath::Sqrt(fValueSystematicASSigmaUp[iPtBin]);
        fValueSystematicPedestalUp[iPtBin] = TMath::Sqrt(fValueSystematicPedestalUp[iPtBin]);
        
        fValueSystematicNSYieldLow[iPtBin] = TMath::Sqrt(fValueSystematicNSYieldLow[iPtBin]);
        fValueSystematicNSSigmaLow[iPtBin] = TMath::Sqrt(fValueSystematicNSSigmaLow[iPtBin]);
        fValueSystematicASYieldLow[iPtBin] = TMath::Sqrt(fValueSystematicASYieldLow[iPtBin]);
        fValueSystematicASSigmaLow[iPtBin] = TMath::Sqrt(fValueSystematicASSigmaLow[iPtBin]);
        fValueSystematicPedestalLow[iPtBin] = TMath::Sqrt(fValueSystematicPedestalLow[iPtBin]);
        //cout << "Output yield = " << fValueSystematicNSYieldUp[iPtBin] << endl;
        
      //  cout << "The systematic error on the yield (" << iPtBin << ") = " << fValueSystematicNSYield[iPtBin] << endl;
    }
    
    fCanvasTotalSystematicSourcesNSYield->cd();
    DrawTotalBaselineSystematicOnCanvas(fSystematicSourcesNSYield,fValueSystematicNSYieldLow, fValueSystematicNSYieldUp ,TotalSystematicsLegendNSyield,kFALSE);
    
    fCanvasTotalSystematicSourcesNSSigma->cd();
    DrawTotalBaselineSystematicOnCanvas(fSystematicSourcesNSSigma,fValueSystematicNSSigmaLow, fValueSystematicNSSigmaUp ,TotalSystematicsLegendNSSigma,kFALSE);
    
    fCanvasTotalSystematicSourcesASYield->cd();
    DrawTotalBaselineSystematicOnCanvas(fSystematicSourcesASYield,fValueSystematicASYieldLow, fValueSystematicASYieldUp ,TotalSystematicsLegendASyield,kFALSE);
    
    fCanvasTotalSystematicSourcesASSigma->cd();
    DrawTotalBaselineSystematicOnCanvas(fSystematicSourcesASSigma,fValueSystematicASSigmaLow, fValueSystematicASSigmaUp ,TotalSystematicsLegendASSigma,kFALSE);
    
    fCanvasTotalSystematicSourcesPedestal->cd();
    DrawTotalBaselineSystematicOnCanvas(fSystematicSourcesPedestal,fValueSystematicPedestalLow, fValueSystematicPedestalUp ,TotalSystematicsLegendPedestal,kFALSE);
    
    
    fCanvasSystematicSourcesNSYield->Write();
    fCanvasSystematicSourcesNSSigma->Write();
    fCanvasSystematicSourcesASYield->Write();
    fCanvasSystematicSourcesASSigma->Write();
    fCanvasSystematicSourcesPedestal->Write();
    
    fCanvasTotalSystematicSourcesNSYield->Write();
    fCanvasTotalSystematicSourcesNSSigma->Write();
    fCanvasTotalSystematicSourcesASYield->Write();
    fCanvasTotalSystematicSourcesASSigma->Write();
    fCanvasSystematicSourcesPedestal->Write();
   

    // add here part to compute the final systematics
    //
    
    
   // fOutputFile->Close();
    
    
    return kTRUE;
}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::DrawTotalSystematicsOnCanvas(TString variable, TLegend * legend){
    
    Float_t ptmean[fVecSize];
    Float_t ptwidth[fVecSize];
    
    Float_t ValueBaselineSyst[fVecSize];
    Float_t ValueCorrSyst[fVecSize];
    Float_t Valuev2Syst[fVecSize];
    
    Float_t ValueBaselineSystUp[fVecSize];
    Float_t ValueCorrSystUp[fVecSize];
    Float_t Valuev2SystUp[fVecSize];
    
    Float_t ValueBaselineSystLow[fVecSize];
    Float_t ValueCorrSystLow[fVecSize];
    Float_t Valuev2SystLow[fVecSize];
    
    TH1D * dummyhisto = (TH1D*)fValueHistoNSYield->Clone(Form("dummy%s",variable.Data()));
    dummyhisto->Reset();
    
    for(Int_t i = 0; i<fVecSize; i++){
        dummyhisto->SetBinContent(i+1,1);
        
        ptmean[i] = 0.5*(fVecLowEdgeDpt[i]+fVecUpEdgeDpt[i]);
        ptwidth[i] = 0.5*(fVecUpEdgeDpt[i]-fVecLowEdgeDpt[i]);
        
        ValueBaselineSyst[i]=1;
        ValueCorrSyst[i]=1;
        Valuev2Syst[i]=1;
        
        if(!fIspPb) {
         Valuev2SystLow[i] = 0;
            Valuev2SystUp[i] = 0;
        }
        if(variable == "NSYield"){
            ValueBaselineSystUp[i] = fValueSystematicBaselineNSYield[i];
            ValueCorrSystUp[i] = fVecMaxCorrelatedSyst[i];
            if(fIspPb){
            if(fSystValuev2NSYield[i]>=0){
                Valuev2SystUp[i] = TMath::Abs(fSystValuev2NSYield[i]);
                Valuev2SystLow[i] = 0;
                }
            }
            ValueBaselineSystLow[i] = fValueSystematicBaselineNSYield[i];
            ValueCorrSystLow[i] = fVecMinCorrelatedSyst[i];
            if(fIspPb){
                if(fSystValuev2NSYield[i]<0) {Valuev2SystLow[i] = TMath::Abs(fSystValuev2NSYield[i]);
                Valuev2SystUp[i] = 0;
                }
            }
        }
        if(variable == "NSSigma"){
            ValueBaselineSystUp[i] = fValueSystematicBaselineNSSigma[i];
            ValueCorrSystUp[i] = fVecMaxCorrelatedSyst[i];
            if(fIspPb){
            if(fSystValuev2NSSigma[i]>=0) {
                Valuev2SystUp[i] = TMath::Abs(fSystValuev2NSSigma[i]);
                Valuev2SystLow[i] = 0;
                }
            }
            ValueBaselineSystLow[i] = fValueSystematicBaselineNSSigma[i];
            ValueCorrSystLow[i] = fVecMinCorrelatedSyst[i];
            if(fIspPb){
            if(fSystValuev2NSSigma[i]<0){ Valuev2SystLow[i] = TMath::Abs(fSystValuev2NSSigma[i]);
                Valuev2SystUp[i] = 0;
                }
            }
        }
        if(variable == "ASYield"){
            ValueBaselineSystUp[i] = fValueSystematicBaselineASYield[i];
            ValueCorrSystUp[i] = fVecMaxCorrelatedSyst[i];
            if(fIspPb){
            if(fSystValuev2ASYield[i]>=0) {Valuev2SystUp[i] = TMath::Abs(fSystValuev2ASYield[i]);
                 Valuev2SystLow[i] = 0;
                }
            }
            ValueBaselineSystLow[i] = fValueSystematicBaselineASYield[i];
            ValueCorrSystLow[i] = fVecMinCorrelatedSyst[i];
                if(fIspPb){
                    if(fSystValuev2ASYield[i]<0) {Valuev2SystLow[i] = TMath::Abs(fSystValuev2ASYield[i]);
                    Valuev2SystUp[i] = 0;
                    }
                }
        }
        if(variable == "ASSigma"){
            ValueBaselineSystUp[i] = fValueSystematicBaselineASSigma[i];
            ValueCorrSystUp[i] = fVecMaxCorrelatedSyst[i];
            if(fIspPb){
            if(fSystValuev2ASSigma[i]>=0) {
                Valuev2SystUp[i] = TMath::Abs(fSystValuev2ASSigma[i]);
                Valuev2SystLow[i] = 0;
                }
            }
            ValueBaselineSystLow[i] = fValueSystematicBaselineASSigma[i];
            ValueCorrSystLow[i] = fVecMinCorrelatedSyst[i];
            if(fIspPb){
            if(fSystValuev2ASSigma[i]<0){ Valuev2SystLow[i] = TMath::Abs(fSystValuev2ASSigma[i]);
                Valuev2SystUp[i] = 0;
                }
            }
        }
        
        if(variable == "Pedestal"){
            ValueBaselineSystUp[i] = fValueSystematicBaselinePedestal[i];
            ValueCorrSystUp[i] = fVecMaxCorrelatedSyst[i];
            if(fIspPb){
            if(fSystValuev2Pedestal[i]>=0) {
                Valuev2SystUp[i] = TMath::Abs(fSystValuev2Pedestal[i]);
                 Valuev2SystLow[i] = 0;
                }
            }
            ValueBaselineSystLow[i] = fValueSystematicBaselinePedestal[i];
            ValueCorrSystLow[i] = fVecMinCorrelatedSyst[i];
            if(fIspPb){
            if(fSystValuev2Pedestal[i]<0) {Valuev2SystLow[i] = TMath::Abs(fSystValuev2Pedestal[i]);
                Valuev2SystUp[i] = 0;
                }
            }
        }
        
    }
    
  
    
    dummyhisto->GetYaxis()->SetRangeUser(0,3);
    dummyhisto->SetTitle("; p_{T} (D) GeV/c");
    
    cout << "defining tgraphassym errors" << endl;
    TGraphAsymmErrors * SystematicsBaseline = new TGraphAsymmErrors(fVecSize,ptmean,ValueBaselineSyst,ptwidth,ptwidth,ValueBaselineSystLow,ValueBaselineSystUp);
    TGraphAsymmErrors * SystematicsCorrelDPhi = new TGraphAsymmErrors(fVecSize,ptmean,ValueCorrSyst,ptwidth,ptwidth,ValueCorrSystLow,ValueCorrSystUp);
    TGraphAsymmErrors * Systematicsv2 = new TGraphAsymmErrors(fVecSize,ptmean,Valuev2Syst,ptwidth,ptwidth,Valuev2SystLow,Valuev2SystUp);
    SystematicsBaseline->SetName(Form("grSystBas%s",variable.Data()));
    SystematicsCorrelDPhi->SetName(Form("grCorrDphi%s",variable.Data()));
    Systematicsv2->SetName(Form("grv2%s",variable.Data()));
    
  
    legend->AddEntry(SystematicsBaseline,"Systematics from baseline estimation","f");
    if(fUseCorrelatedSystematics){
        if((variable== "NSSigma"||variable== "ASSigma")){
            if(fUseCorrelatedSystematicsForWidths)legend->AddEntry(SystematicsCorrelDPhi,"Systematics from #Delta#varphi scaling uncertainty","f");
        }
        else legend->AddEntry(SystematicsCorrelDPhi,"Systematics from #Delta#varphi scaling uncertainty","f");
    }
    if(fIspPb) legend->AddEntry(Systematicsv2,"Systematics from v_{2}(D) hypothesis","f");
    

    
    SystematicsBaseline->SetFillColor(2);
    SystematicsCorrelDPhi->SetFillColor(3);
    Systematicsv2->SetFillColor(4);
    
    SystematicsBaseline->SetLineColor(2);
    SystematicsCorrelDPhi->SetLineColor(3);
    Systematicsv2->SetLineColor(4);
    
    SystematicsBaseline->SetFillStyle(3004);
    SystematicsCorrelDPhi->SetFillStyle(3005);
    Systematicsv2->SetFillStyle(3003);
     cout << "drawing" << endl;
    dummyhisto->Draw();
    SystematicsBaseline->Draw("E2same");
    if(fUseCorrelatedSystematics){
        if((variable== "NSSigma"||variable== "ASSigma")){
            if(fUseCorrelatedSystematicsForWidths)SystematicsCorrelDPhi->Draw("E2same");
        }
        else SystematicsCorrelDPhi->Draw("E2same");
    }
    if(fIspPb) Systematicsv2->Draw("E2same");
    legend->Draw("same");
 
    
}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::DrawBaselineSystematicsOnCanvas(TH1D * histoinput, Int_t iSystMode, Double_t * array, TLegend * legend){
    
    TString name = histoinput->GetName();
    //name += Form("%d",iSystMode);
    
    TH1D * histo = (TH1D*)histoinput->Clone(Form("%s_%d",name.Data(),iSystMode));
    
    

    
   // cout << "Name = " << name << endl;
 //   histo->SetName(Form("culone_%d",iSystMode));
    //histo->SetName(name.Data());
    
    TString title = name;
    title += "; p_{T}(D) GeV/c";
    
    histo->SetTitle(title.Data());
    Int_t index = -1;
    for(Int_t iPtBin=0; iPtBin<fVecSize; iPtBin++){
        
     //   cout << "Name = " << name << "; iPtBin = " << iPtBin << "; iSystMode = " << iSystMode << ", value " << array[index] << endl;
        
        index = GetBinIndex(iPtBin, iSystMode);
        histo->SetBinContent(iPtBin+1,1);
        histo->SetBinError(iPtBin+1,array[index]);
    }
    histo->SetFillColor(iSystMode);
    histo->SetFillStyle(3000+iSystMode);

    TString suffix = "";
    if(fVecSystMode[iSystMode] == kLowestPoint) suffix = "Lowest Point";
    if(fVecSystMode[iSystMode] == kNLowest) suffix = Form("%dLowest Point",fNMinPointsBaselineEstimationRange[iSystMode]);
    if(fVecSystMode[iSystMode] == k2PointsAtPiHalf) suffix = "2 Points at pi half";
    if(fVecSystMode[iSystMode] == k4PointsAtPiHalf) suffix = "4 Points at pi half";
    if(fVecSystMode[iSystMode] == kTransverse) suffix = Form("Baseline from region %.3f#pi - %.3f#pi ",fMinBaselineEstimationRange[iSystMode]/TMath::Pi(),fMaxBaselineEstimationRange[iSystMode]/TMath::Pi());
    if(fVecSystMode[iSystMode] == kTransverseUppStatUnc) suffix = Form("Basel. from %.3f#pi - %.3f#pi + #sigma_{unc}",fMinBaselineEstimationRange[iSystMode]/TMath::Pi(),fMaxBaselineEstimationRange[iSystMode]/TMath::Pi());
    if(fVecSystMode[iSystMode] == kTransverseLowStatUnc) suffix = Form("Basel. from %.3f#pi - %.3f#pi - #sigma_{unc}",fMinBaselineEstimationRange[iSystMode]/TMath::Pi(),fMaxBaselineEstimationRange[iSystMode]/TMath::Pi());
    if(fVecSystMode[iSystMode] == kBinCount) suffix = "Bin Counting";
    if(fVecSystMode[iSystMode] == kMinVar) suffix = "Min Variation";
    if(fVecSystMode[iSystMode] == kMaxVar) suffix = "Max Variation";
    
    TString legendentry = "";
    legendentry += suffix;
    legend->AddEntry(histo,legendentry,"f");
    
    histo->Draw("e2same");
 

}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::DrawRMSBaselineSystematicOnCanvas(TH1D * histoinput, Double_t * arraymin, Double_t * arraymax , TLegend * legend, Bool_t isSystBaseline){
    
    TString name = histoinput->GetName();
    //name += Form("%d",iSystMode);
    
    TH1D * histoup = (TH1D*)histoinput->Clone(Form("%s_TotalUp",name.Data()));
    TH1D * histodown = (TH1D*)histoinput->Clone(Form("%s_TotalDown",name.Data()));
    histoup->Reset();
    histodown->Reset();
  //  cout << "Name = " << name << endl;
    //   histo->SetName(Form("culone_%d",iSystMode));
    //histo->SetName(name.Data());
    
    TString title = name;
    title += ";p_{T} (D) GeV/c";
    
    histoup->SetTitle(title.Data());
    Int_t index = -1;
    for(Int_t iPtBin=0; iPtBin<fVecSize; iPtBin++){
        
       // cout << "Name = " << name << "; iPtBin = " << iPtBin << "; iSystMode = " << iSystMode << ", value " << array[index] << endl;
        
        histoup->SetBinContent(iPtBin+1,1+arraymax[iPtBin]);
        histodown->SetBinContent(iPtBin+1,1-arraymin[iPtBin]);
       // cout << "iptbin = " << iPtBin << "value max syste = " <<  arraymax[iPtBin] << endl;
      //  cout << "iptbin = " << iPtBin << "value min syste = " <<  arraymin[iPtBin] << endl;
    }
    histoup->SetLineColor(3);
    histodown->SetLineColor(3);
    histoup->SetLineWidth(2);
    histodown->SetLineWidth(2);
    
   if(fCombineSystematics == kRMS && isSystBaseline) legend->AddEntry(histoup,"Total Syst on std vars (RMS)","l");
    
   // cout << "drawing histo" << endl;
    histoup->Draw("histsame");
    histodown->Draw("histsame");
    
    
}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::DrawTotalBaselineSystematicOnCanvas(TH1D * histoinput, Double_t * arraymin, Double_t * arraymax , TLegend * legend, Bool_t isSystBaseline){
    
    TString name = histoinput->GetName();
    //name += Form("%d",iSystMode);
    
    TH1D * histoup = (TH1D*)histoinput->Clone(Form("%s_TotalUp",name.Data()));
    TH1D * histodown = (TH1D*)histoinput->Clone(Form("%s_TotalDown",name.Data()));
    histoup->Reset();
    histodown->Reset();
  //  cout << "Name = " << name << endl;
    //   histo->SetName(Form("culone_%d",iSystMode));
    //histo->SetName(name.Data());
    
    TString title = name;
    title += ";p_{T} (D) GeV/c";
    
    histoup->SetTitle(title.Data());
    Int_t index = -1;
    for(Int_t iPtBin=0; iPtBin<fVecSize; iPtBin++){
        
       // cout << "Name = " << name << "; iPtBin = " << iPtBin << "; iSystMode = " << iSystMode << ", value " << array[index] << endl;
        
        histoup->SetBinContent(iPtBin+1,1+arraymax[iPtBin]);
        histodown->SetBinContent(iPtBin+1,1-arraymin[iPtBin]);
       // cout << "iptbin = " << iPtBin << "value max syste = " <<  arraymax[iPtBin] << endl;
      //  cout << "iptbin = " << iPtBin << "value min syste = " <<  arraymin[iPtBin] << endl;
    }
    histoup->SetLineColor(1);
    histodown->SetLineColor(1);
    histoup->SetLineWidth(2);
    histodown->SetLineWidth(2);
    
    

   if(fCombineSystematics == kMax && isSystBaseline) legend->AddEntry(histoup,"Total Syst max variation","l");
   if(fCombineSystematics == kRMS && isSystBaseline) legend->AddEntry(histoup,"Total Syst via RMS ","l");
   if(fCombineSystematics == kEnvelope_RMS_BaselStat && isSystBaseline) legend->AddEntry(histoup,"Total Syst via special envel.","l");
   if(fCombineSystematics == kSumQuadr || !isSystBaseline) legend->AddEntry(histoup,"Total Syst Sum quadrature ","l");
    
   // cout << "drawing histo" << endl;
    histoup->Draw("histsame");
    histodown->Draw("histsame");
    
    
}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::PrintAllSystematicsOnShell(){
 
    Int_t index = -1;
    
 //   SystematicModes{kLowestPoint = 1, k2PointsAtPiHalf = 2, k4PointsAtPiHalf=4,kTransverse = 5, kBinCount = -1, kMinVar=100, kMaxVar=200, kv2Mod = 300};
    for(Int_t iSystMode = 0; iSystMode < fVecSystModesSize; iSystMode++){
        std::cout << "  " << std::endl;
        std::cout << "================================================================= " << std::endl;
        std::cout << "  " << std::endl;

        if(fVecSystMode[iSystMode] == kFree) std::cout << "Systematic mode is free baseline " << std::endl;       
        if(fVecSystMode[iSystMode] == kLowestPoint) std::cout << "Systematic mode is Lowest point in Dphi " << std::endl;
        if(fVecSystMode[iSystMode] == kNLowest) std::cout << "Systematic mode is "<< fNMinPointsBaselineEstimationRange[iSystMode] <<" lowest points in Dphi " << std::endl;
        if(fVecSystMode[iSystMode] == k2PointsAtPiHalf) std::cout << "Systematic mode is two points at pi/2  " << std::endl;
        if(fVecSystMode[iSystMode] == k4PointsAtPiHalf) std::cout << "Systematic mode is 4 points at pi/2  " << std::endl;
        if(fVecSystMode[iSystMode] == kTransverse) std::cout << "Systematic mode is transverse region : (" << fMinBaselineEstimationRange[iSystMode]/TMath::Pi() << "*pi, "<< fMaxBaselineEstimationRange[iSystMode]/TMath::Pi() << "*pi)"<<std::endl;
        if(fVecSystMode[iSystMode] == kTransverseUppStatUnc) std::cout << "Systematic mode is transverse region : (" << fMinBaselineEstimationRange[iSystMode]/TMath::Pi() << "*pi, "<< fMaxBaselineEstimationRange[iSystMode]/TMath::Pi() << "*pi) plus its stat unc"<<std::endl;
        if(fVecSystMode[iSystMode] == kTransverseLowStatUnc) std::cout << "Systematic mode is transverse region : (" << fMinBaselineEstimationRange[iSystMode]/TMath::Pi() << "*pi, "<< fMaxBaselineEstimationRange[iSystMode]/TMath::Pi() << "*pi) minus its stat unc"<<std::endl;
        if(fVecSystMode[iSystMode] == kBinCount) std::cout << "Systematic mode is bin counting (yield only)  " << std::endl;
        if(fVecSystMode[iSystMode] == kMinVar) std::cout << "Systematic mode is fit min variation histo  " << std::endl;
        if(fVecSystMode[iSystMode] == kMaxVar) std::cout << "Systematic mode is fit max variation histo  " << std::endl;
        if(fVecSystMode[iSystMode] == kv2Mod) std::cout << "Systematic mode is baseline v2 modulation  " << std::endl;
        
        
         for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
              index = GetBinIndex(iPtBin, iSystMode);
             if(iPtBin != (fVecSize-1) )std:: cout << " pT D meson (" << fVecLowEdgeDpt[iPtBin] <<"," << fVecUpEdgeDpt[iPtBin] << "); "<< std::endl;
             if(iPtBin == (fVecSize-1) )std:: cout << " pT D meson (" << fVecLowEdgeDpt[iPtBin] <<"," << fVecUpEdgeDpt[iPtBin] << "); " << std::endl;
             std::cout << "Systematic on NS yield due to this source = " << fRatioNSYield[index] << std::endl;
             std::cout << "Systematic on NS sigma due to this source = " << fRatioNSSigma[index] << std::endl;
             std::cout << "Systematic on AS yield due to this source = " << fRatioASYield[index] << std::endl;
             std::cout << "Systematic on AS sigma due to this source = " << fRatioASSigma[index] << std::endl;
             std::cout << "Systematic on pedestal due to this source = " << fRatioPedestal[index] << std::endl;
             
            // fRatioNSSigma[index]
            //  fRatioASYield[index]
            // fRatioASSigma[index]
            // fRatioPedestal[index]
         }
    }
    
    
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << std::endl;
    std::cout << "SYSTEMATICS FROM BASELINE ESTIMATION VALUES " << std::endl;
    for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
        
        if(iPtBin != (fVecSize-1) )std:: cout << " pT D meson (" << fVecLowEdgeDpt[iPtBin] <<"," << fVecUpEdgeDpt[iPtBin] << "); "<< std::endl;
        if(iPtBin == (fVecSize-1) )std:: cout << " pT D meson (" << fVecLowEdgeDpt[iPtBin] <<"," << fVecUpEdgeDpt[iPtBin] << "); " << std::endl;
        std::cout << "Systematic on NS yield (baseline) = +/- " << fValueSystematicBaselineNSYield[iPtBin] << std::endl;
        std::cout << "Systematic on NS sigma (baseline) = +/- " << fValueSystematicBaselineNSSigma[iPtBin] << std::endl;
        std::cout << "Systematic on AS yield (baseline) = +/- " << fValueSystematicBaselineASYield[iPtBin] << std::endl;
        std::cout << "Systematic on AS sigma (baseline) = +/- " << fValueSystematicBaselineASSigma[iPtBin] << std::endl;
        std::cout << "Systematic on pedestal (baseline) = +/- " << fValueSystematicBaselinePedestal[iPtBin] << std::endl;
        
        // fRatioNSSigma[index]
        //  fRatioASYield[index]
        // fRatioASSigma[index]
        // fRatioPedestal[index]
    }
    
    if(fUseCorrelatedSystematics){
        std::cout << " " << std::endl;
        std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << std::endl;
        std::cout << "SYSTEMATICS FROM DPHI NORMALIZATION VALUES " << std::endl;
        for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
        
            if(iPtBin != (fVecSize-1) )std:: cout << " pT D meson (" << fVecLowEdgeDpt[iPtBin] <<"," << fVecUpEdgeDpt[iPtBin] << "); "<< std::endl;
            if(iPtBin == (fVecSize-1) )std:: cout << " pT D meson (" << fVecLowEdgeDpt[iPtBin] <<"," << fVecUpEdgeDpt[iPtBin] << "); " << std::endl;
            std::cout << "Systematic on NS yield (DPHI NORM)) = +" << fVecMaxCorrelatedSyst[iPtBin] << " , -" << fVecMinCorrelatedSyst[iPtBin] <<  std::endl;
           if(fUseCorrelatedSystematicsForWidths) std::cout << "Systematic on NS sigma (DPHI NORM) = +" << fVecMaxCorrelatedSyst[iPtBin] << " , -" << fVecMinCorrelatedSyst[iPtBin]<< std::endl;
            else std::cout << "Systematic on NS sigma (DPHI NORM)) = +" << 0<< " , -" << 0<< std::endl;
            std::cout << "Systematic on AS yield (DPHI NORM)) = +" << fVecMaxCorrelatedSyst[iPtBin] << " , -" << fVecMinCorrelatedSyst[iPtBin] << std::endl;
           if(fUseCorrelatedSystematicsForWidths) std::cout << "Systematic on AS sigma (DPHI NORM) = +" << fVecMaxCorrelatedSyst[iPtBin] << " , -" << fVecMinCorrelatedSyst[iPtBin]<< std::endl;
            else std::cout << "Systematic on AS sigma (DPHI NORM) = +" << 0 << " , -" << 0<< std::endl;
            std::cout << "Systematic on pedestal (DPHI NORM) = +" << fVecMaxCorrelatedSyst[iPtBin] << " , -" << fVecMinCorrelatedSyst[iPtBin]<< std::endl;
            
        // fRatioNSSigma[index]
        //  fRatioASYield[index]
        // fRatioASSigma[index]
        // fRatioPedestal[index]
        }

    }
    if(fIspPb){
        std::cout << " " << std::endl;
        std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << std::endl;
        std::cout << "SYSTEMATICS FROM v2 Hypothesis VALUES " << std::endl;
        for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
        
            if(iPtBin != (fVecSize-1) )std:: cout << " pT D meson (" << fVecLowEdgeDpt[iPtBin] <<"," << fVecUpEdgeDpt[iPtBin] << "); "<< std::endl;
            if(iPtBin == (fVecSize-1) )std:: cout << " pT D meson (" << fVecLowEdgeDpt[iPtBin] <<"," << fVecUpEdgeDpt[iPtBin] << "); " << std::endl;
            if(fSystValuev2NSYield[iPtBin]>=0)  std::cout << "Systematic on NS yield (v2 hyp.) = +" << fSystValuev2NSYield[iPtBin] << std::endl;
            if(fSystValuev2NSYield[iPtBin]<0)  std::cout << "Systematic on NS yield (v2 hyp.) = " << fSystValuev2NSYield[iPtBin] << std::endl;
            
            if(fSystValuev2NSSigma[iPtBin]>=0)  std::cout << "Systematic on NS Sigma (v2 hyp.) = +" << fSystValuev2NSSigma[iPtBin] << std::endl;
            if(fSystValuev2NSSigma[iPtBin]<0)  std::cout << "Systematic on NS Sigma (v2 hyp.) = " << fSystValuev2NSSigma[iPtBin] << std::endl;
            
            if(fSystValuev2ASYield[iPtBin]>=0)  std::cout << "Systematic on AS yield (v2 hyp.) = +" << fSystValuev2ASYield[iPtBin] << std::endl;
            if(fSystValuev2ASYield[iPtBin]<0)  std::cout << "Systematic on AS yield (v2 hyp.) = " << fSystValuev2ASYield[iPtBin] << std::endl;
            
            if(fSystValuev2ASSigma[iPtBin]>=0)  std::cout << "Systematic on AS Sigma (v2 hyp.) = +" << fSystValuev2ASSigma[iPtBin] << std::endl;
            if(fSystValuev2ASSigma[iPtBin]<0)  std::cout << "Systematic on AS Sigma (v2 hyp.) = " << fSystValuev2ASSigma[iPtBin] << std::endl;
            
            if(fSystValuev2Pedestal[iPtBin]>=0)  std::cout << "Systematic on Pedestal(v2 hyp.) = +" << fSystValuev2Pedestal[iPtBin] << std::endl;
            if(fSystValuev2Pedestal[iPtBin]<0)  std::cout << "Systematic on Pedestal (v2 hyp.) = " << fSystValuev2Pedestal[iPtBin] << std::endl;
        
        // fRatioNSSigma[index]
        //  fRatioASYield[index]
        // fRatioASSigma[index]
        // fRatioPedestal[index]
        }
    
    }
    std::cout << " " << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << std::endl;
    std::cout << "FINAL SYSTEMATICS VALUES " << std::endl;
    
    
    
    for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
       
        if(iPtBin != (fVecSize-1) )std:: cout << " pT D meson (" << fVecLowEdgeDpt[iPtBin] <<"," << fVecUpEdgeDpt[iPtBin] << "); "<< std::endl;
        if(iPtBin == (fVecSize-1) )std:: cout << " pT D meson (" << fVecLowEdgeDpt[iPtBin] <<"," << fVecUpEdgeDpt[iPtBin] << "); " << std::endl;
        std::cout << "Systematic on NS yield (total) = +" << fValueSystematicNSYieldUp[iPtBin] << " , -" << fValueSystematicNSYieldLow[iPtBin] <<  std::endl;
        std::cout << "Systematic on NS sigma (total) = +" << fValueSystematicNSSigmaUp[iPtBin] << " , -" << fValueSystematicNSSigmaLow[iPtBin]<< std::endl;
        std::cout << "Systematic on AS yield (total) = +" << fValueSystematicASYieldUp[iPtBin] << " , -" << fValueSystematicASYieldLow[iPtBin] << std::endl;
        std::cout << "Systematic on AS sigma (total) = +" << fValueSystematicASSigmaUp[iPtBin] << " , -" << fValueSystematicASSigmaLow[iPtBin]<< std::endl;
        std::cout << "Systematic on pedestal (total) = +" << fValueSystematicPedestalUp[iPtBin] << " , -" << fValueSystematicPedestalLow[iPtBin]<< std::endl;
        
        // fRatioNSSigma[index]
        //  fRatioASYield[index]
        // fRatioASSigma[index]
        // fRatioPedestal[index]
    }
    
    std::cout << " " << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< " << std::endl;
    
}
//_______________________________________________________________________________
Bool_t AliHFCorrFitSystematics::DrawFinalCorrelationPlot(Bool_t drawNSy, Bool_t drawNSsigma,Bool_t drawASy, Bool_t drawASsigma,Bool_t drawPed){
 
    TCanvas ** cFinalCorrelation = new TCanvas*[fVecSize];
    
    TPaveText * pave = NULL;
    
    for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
        pave = new TPaveText(0.5,0.6,0.9,0.8,"NDC");
        pave->SetName(Form("paveCorrelations_%d",iPtBin));
        pave->AddText(Form("#chi^{2}/NDF = %.1f/%d",fFitFunctions[iPtBin]->GetChisquare(),fFitFunctions[iPtBin]->GetNDF()));
        if(drawNSy){
            if(fValueSystematicNSYieldUp == fValueSystematicNSYieldLow) {pave->AddText(Form("NS Y = %.2f #pm %.2f (stat) #pm %.2f (syst) ",fReferenceHistoNSYield->GetBinContent(iPtBin+1), fReferenceHistoNSYield->GetBinError(iPtBin+1), fValueSystematicNSYieldUp[iPtBin]*fReferenceHistoNSYield->GetBinContent(iPtBin+1)));
            }
            else pave->AddText(Form("NS Yield = %.2f #pm %.2f (stat) #pm ^{%.2f}_{%.2f} (syst) ",fReferenceHistoNSYield->GetBinContent(iPtBin+1), fReferenceHistoNSYield->GetBinError(iPtBin+1), fValueSystematicNSYieldUp[iPtBin]*fReferenceHistoNSYield->GetBinContent(iPtBin+1),fValueSystematicNSYieldLow[iPtBin]*fReferenceHistoNSYield->GetBinContent(iPtBin+1)));
        }
        if(drawNSsigma){
            if(fValueSystematicNSSigmaUp == fValueSystematicNSSigmaLow) {pave->AddText(Form("NS Y = %.2f #pm %.2f (stat) #pm %.2f (syst) ",fReferenceHistoNSSigma->GetBinContent(iPtBin+1), fReferenceHistoNSSigma->GetBinError(iPtBin+1), fValueSystematicNSSigmaUp[iPtBin]*fReferenceHistoNSSigma->GetBinContent(iPtBin+1)));
            }
            else pave->AddText(Form("NS Sigma = %.2f #pm %.2f (stat) #pm ^{%.2f}_{%.2f} (syst) ",fReferenceHistoNSSigma->GetBinContent(iPtBin+1), fReferenceHistoNSSigma->GetBinError(iPtBin+1), fValueSystematicNSSigmaUp[iPtBin]*fReferenceHistoNSSigma->GetBinContent(iPtBin+1),fValueSystematicNSSigmaLow[iPtBin]*fReferenceHistoNSSigma->GetBinContent(iPtBin+1)));
        }
        
        if(drawASy){
            if(fValueSystematicASYieldUp == fValueSystematicASYieldLow) {pave->AddText(Form("AS Y = %.2f #pm %.2f (stat) #pm %.2f (syst) ",fReferenceHistoASYield->GetBinContent(iPtBin+1), fReferenceHistoASYield->GetBinError(iPtBin+1), fValueSystematicASYieldUp[iPtBin]*fReferenceHistoASYield->GetBinContent(iPtBin+1)));
            }
            else pave->AddText(Form("AS Yield = %.2f #pm %.2f (stat) #pm ^{%.2f}_{%.2f} (syst) ",fReferenceHistoASYield->GetBinContent(iPtBin+1), fReferenceHistoASYield->GetBinError(iPtBin+1), fValueSystematicASYieldUp[iPtBin]*fReferenceHistoASYield->GetBinContent(iPtBin+1),fValueSystematicASYieldLow[iPtBin]*fReferenceHistoASYield->GetBinContent(iPtBin+1)));
        }
        if(drawASsigma){
            if(fValueSystematicASSigmaUp == fValueSystematicASSigmaLow) {pave->AddText(Form("AS Y = %.2f #pm %.2f (stat) #pm %.2f (syst) ",fReferenceHistoASSigma->GetBinContent(iPtBin+1), fReferenceHistoASSigma->GetBinError(iPtBin+1), fValueSystematicASSigmaUp[iPtBin]*fReferenceHistoASSigma->GetBinContent(iPtBin+1)));
            }
            else pave->AddText(Form("AS Sigma = %.2f #pm %.2f (stat) #pm ^{%.2f}_{%.2f} (syst) ",fReferenceHistoASSigma->GetBinContent(iPtBin+1), fReferenceHistoASSigma->GetBinError(iPtBin+1), fValueSystematicASSigmaUp[iPtBin]*fReferenceHistoASSigma->GetBinContent(iPtBin+1),fValueSystematicASSigmaLow[iPtBin]*fReferenceHistoASSigma->GetBinContent(iPtBin+1)));
        }
        
        if(drawPed){
            if(fValueSystematicPedestalUp == fValueSystematicPedestalLow) {pave->AddText(Form("Pedestal = %.2f #pm %.2f (stat) #pm %.2f (syst) ",fReferenceHistoPedestal->GetBinContent(iPtBin+1), fReferenceHistoPedestal->GetBinError(iPtBin+1), fValueSystematicPedestalUp[iPtBin]*fReferenceHistoPedestal->GetBinContent(iPtBin+1)));
            }
            else pave->AddText(Form("Pedestal = %.2f #pm %.2f (stat) #pm ^{%.2f}_{%.2f} (syst) ",fReferenceHistoPedestal->GetBinContent(iPtBin+1), fReferenceHistoPedestal->GetBinError(iPtBin+1), fValueSystematicPedestalUp[iPtBin]*fReferenceHistoPedestal->GetBinContent(iPtBin+1),fValueSystematicPedestalLow[iPtBin]*fReferenceHistoPedestal->GetBinContent(iPtBin+1)));
        }
        
        pave->SetFillStyle(0);
        pave->SetBorderSize(0);
        pave->SetTextSize(0.023);
        
        
      
        
        TLatex *tlDmesonpt =new TLatex(0.17,0.6,Form("%.1f < p_{T}(D) < %.1f GeV/c",fVecLowEdgeDpt[iPtBin],fVecUpEdgeDpt[iPtBin]));
        tlDmesonpt->SetNDC();
        tlDmesonpt->SetTextFont(42);
        tlDmesonpt->SetTextSize(0.03);
        
        
        cFinalCorrelation[iPtBin] = new TCanvas(Form("cFinalCorrelation%d",iPtBin),Form("cFinalCorrelation%d",iPtBin),0,0,1000,1000);
        cFinalCorrelation[iPtBin]->cd();
        cFinalCorrelation[iPtBin]->SetTicks();
        fVecHisto[iPtBin]->Draw("ep");
        if(!fIsFittingTemplate && fCorrelationGraphAsymmErrors[iPtBin])fCorrelationGraphAsymmErrors[iPtBin]->Draw("E2same");
        DefinePaveText();
        fFitFunctions[iPtBin]->Draw("same");
        pave->Draw("same");
        tlDmesonpt->Draw("same");
        
        
        //cout << "saving " << iPtBin << endl;
        if(fSaveDotC) cFinalCorrelation[iPtBin]->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),cFinalCorrelation[iPtBin]->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
        if(fSaveEps) cFinalCorrelation[iPtBin]->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),cFinalCorrelation[iPtBin]->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
        if(fSavePdf) cFinalCorrelation[iPtBin]->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),cFinalCorrelation[iPtBin]->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
        if(fSavePng) cFinalCorrelation[iPtBin]->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),cFinalCorrelation[iPtBin]->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
        if(fSaveRoot) cFinalCorrelation[iPtBin]->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),cFinalCorrelation[iPtBin]->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
        fOutputFile->cd();
        cFinalCorrelation[iPtBin]->Write();

    }
   
    return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliHFCorrFitSystematics::DrawFinalPlots(Bool_t drawNSy, Bool_t drawNSsigma,Bool_t drawASy, Bool_t drawASsigma,Bool_t drawPed){
   
    cout << "==============================================================" << endl;
    const int size = fVecSize+1;
    Float_t ptbinsout[size+2];
    ptbinsout[0] = 0;
    for(Int_t i = 0; i<fVecSize; i++){
        ptbinsout[i+1] = fVecLowEdgeDpt[i];
    }
    ptbinsout[fVecSize+1] = fVecUpEdgeDpt[fVecSize-1];
    ptbinsout[fVecSize+2] = fVecUpEdgeDpt[fVecSize-1]+4;
    
    Float_t *ptbinspointerout = ptbinsout;
    
    cout << "check me" << endl;
    for(Int_t i = 0; i<size+2; i++){
        
        cout << i << " " << ptbinsout[i] << endl;
    }
    
    
   /* fFinalTrendNSYield = new TH1D("FinalTrendNSYield","; p_{T}(D) GeV/c; NS Yield",fVecSize+2,ptbinspointerout);
    fFinalTrendNSSigma = new TH1D("FinalTrendNSSigma","; p_{T}(D) GeV/c; NS#sigma",fVecSize+2,ptbinspointerout);
    fFinalTrendASYield = new TH1D("FinalTrendASYield"," p_{T}(D) GeV/c; AS Yield",fVecSize+2,ptbinspointerout);
    fFinalTrendASSigma = new TH1D("FinalTrendASSigma"," p_{T}(D) GeV/c; AS#sigma",fVecSize+2,ptbinspointerout);
    fFinalTrendPedestal = new TH1D("FinalTrendPedestal"," p_{T}(D) GeV/c; Pedestal",fVecSize+2,ptbinspointerout);
    */
    Bool_t isDrawn = kFALSE;
    
    TLegend * Legend = new TLegend(0.5,0.7,0.85,0.8);
    Legend->SetFillColor(0);
    Legend->SetTextSize(0.025);
    Legend->SetBorderSize(0);
   

    fFullSystematicsNSYield = new TGraphAsymmErrors();
    fFullSystematicsNSYield->SetName("fFullSystematicsNSYield");
    if(fPlotV2SystSeparately){
      fv2SystematicsNSYield= new TGraphAsymmErrors();
      fv2SystematicsNSYield->SetName("fv2SystematicsNSYield");
    }
    fFullSystematicsNSSigma = new TGraphAsymmErrors();
    fFullSystematicsNSSigma->SetName("fFullSystematicsNSSigma");
    if(fPlotV2SystSeparately){
      fv2SystematicsNSSigma = new TGraphAsymmErrors();
      fv2SystematicsNSSigma->SetName("fv2SystematicsNSSigma");
    }
    fFullSystematicsASYield = new TGraphAsymmErrors();
    fFullSystematicsASYield->SetName("fFullSystematicsASYield");
    if(fPlotV2SystSeparately){
      fv2SystematicsASYield = new TGraphAsymmErrors();
      fv2SystematicsASYield->SetName("fv2SystematicsASYield");
    }
    fFullSystematicsASSigma = new TGraphAsymmErrors();
    fFullSystematicsASSigma->SetName("fFullSystematicsASSigma");
    if(fPlotV2SystSeparately) {
      fv2SystematicsASSigma = new TGraphAsymmErrors();
      fv2SystematicsASSigma->SetName("fv2SystematicsASSigma");
    }
    fFullSystematicsPedestal = new TGraphAsymmErrors();
    fFullSystematicsPedestal->SetName("fFullSystematicsPedestal");
    if(fPlotV2SystSeparately){
      fv2SystematicsPedestal = new TGraphAsymmErrors();
      fv2SystematicsPedestal->SetName("fv2SystematicsPedestal");
    }
    fBaselineVariationSystematicsPedestal = new TGraphAsymmErrors();
    fBaselineVariationSystematicsPedestal->SetName("fBaselineVariationSystematicsPedestal");
    
    Double_t Avpt = 0;
    Color_t sysboxcolor;
    
    if(drawNSy){
        sysboxcolor = kGreen-2;
        fCanvasFinalTrendNSYield = new TCanvas("CanvasFinalTrendNSYield","CanvasFinalTrendNSYield",0,0,1000,800);
        fCanvasFinalTrendNSYield->SetTickx();
        fCanvasFinalTrendNSYield->SetTicky();
        fFinalTrendNSYield = new TH1D("FinalTrendNSYield","; p_{T}(D) GeV/c; NS Yield",fVecSize+2,ptbinspointerout);
        SetHisto(fFinalTrendNSYield,fReferenceHistoNSYield);
        fCanvasFinalTrendNSYield->cd();
        fFinalTrendNSYield->Draw("ep");
        DefinePaveText();
	fFullSystematicsNSYield->SetFillColor(sysboxcolor);
	fFullSystematicsNSYield->SetFillStyle(0);
	fFullSystematicsNSYield->SetLineColor(2);

	if(fPlotV2SystSeparately){
	  fv2SystematicsNSYield->SetFillColor(kGreen-2);
	  fv2SystematicsNSYield->SetFillStyle(3001);
	  fv2SystematicsNSYield->SetLineColor(kGreen-2);
	}
        for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
            
            Avpt = (fVecLowEdgeDpt[iPtBin]+fVecUpEdgeDpt[iPtBin])*0.5;
            fFullSystematicsNSYield->SetPoint(iPtBin, Avpt,fFinalTrendNSYield->GetBinContent(iPtBin+2)); 
            fFullSystematicsNSYield->SetPointError(iPtBin,0.9,0.9,fValueSystematicNSYieldLow[iPtBin]*fFinalTrendNSYield->GetBinContent(iPtBin+2),fValueSystematicNSYieldUp[iPtBin]*fFinalTrendNSYield->GetBinContent(iPtBin+2));
            
            
                if(fPlotV2SystSeparately){
		  fv2SystematicsNSYield->SetPoint(iPtBin,Avpt,fFinalTrendNSYield->GetBinContent(iPtBin+2));
		  if(fSystValuev2NSYield[iPtBin]>0.){
		    fv2SystematicsNSYield->SetPointError(iPtBin,0.9,0.9,0,fSystValuev2NSYield[iPtBin]*fFinalTrendNSYield->GetBinContent(iPtBin+2));
		  }
		  else{
		    fv2SystematicsNSYield->SetPointError(iPtBin,0.9,0.9,TMath::Abs(fSystValuev2NSYield[iPtBin]*fFinalTrendNSYield->GetBinContent(iPtBin+2)),0);
		  }
                   // if(fSystValuev2NSYield[iPtBin]>=0) fv2SystematicsNSYield[iPtBin] =new TBox(Avpt-0.9,fFinalTrendNSYield->GetBinContent(iPtBin+2),Avpt+0.9,fFinalTrendNSYield->GetBinContent(iPtBin+2)+fSystValuev2NSYield[iPtBin]*fFinalTrendNSYield->GetBinContent(iPtBin+2));
                   // else fv2SystematicsNSYield[iPtBin] =new TBox(Avpt-0.9,fFinalTrendNSYield->GetBinContent(iPtBin+2)-fSystValuev2NSYield[iPtBin]*fFinalTrendNSYield->GetBinContent(iPtBin+2),Avpt+0.9,fFinalTrendNSYield->GetBinContent(iPtBin+2));
                }
        
	}
	fFullSystematicsNSYield->Draw("E2");
	if(fPlotV2SystSeparately){
	  fv2SystematicsNSYield->Draw("E2");
          Legend->Clear();
	  Legend->AddEntry(fFullSystematicsNSYield,"Total Systematic Uncert.","l");
	  Legend->AddEntry(fv2SystematicsNSYield,"v_{2} Systematic Uncert.","f");
	  Legend->Draw("same");
	}
        if(fPlotV2SystSeparately) fFinalTrendNSYield->Draw("epsame");
        fCanvasFinalTrendNSYield->Write();
        
    }
    if(drawNSsigma){
        sysboxcolor = kBlue-4;
        fCanvasFinalTrendNSSigma = new TCanvas("CanvasFinalTrendNSSigma","CanvasFinalTrendNSSigma",0,0,1000,800);
        fCanvasFinalTrendNSSigma->SetTickx();
        fCanvasFinalTrendNSSigma->SetTicky();
	fFinalTrendNSSigma = new TH1D("FinalTrendNSSigma","; p_{T}(D) GeV/c; NS#sigma",fVecSize+2,ptbinspointerout);
	SetHisto(fFinalTrendNSSigma,fReferenceHistoNSSigma);
        fCanvasFinalTrendNSSigma->cd();
        fFinalTrendNSSigma->Draw("ep");
        DefinePaveText();
	fFullSystematicsNSSigma->SetFillColor(sysboxcolor);
	fFullSystematicsNSSigma->SetFillStyle(0);
	fFullSystematicsNSSigma->SetLineColor(2);
	if(fPlotV2SystSeparately){
	  fv2SystematicsNSSigma->SetFillColor(kGreen-2);
	  fv2SystematicsNSSigma->SetFillStyle(3001);
	  fv2SystematicsNSSigma->SetLineColor(kGreen-2);
        }
        for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
            
            Avpt = (fVecLowEdgeDpt[iPtBin]+fVecUpEdgeDpt[iPtBin])*0.5;
            fFullSystematicsNSSigma->SetPoint(iPtBin,Avpt,fFinalTrendNSSigma->GetBinContent(iPtBin+2));
            fFullSystematicsNSSigma->SetPointError(iPtBin,0.9,0.9,fValueSystematicNSSigmaLow[iPtBin]*fFinalTrendNSSigma->GetBinContent(iPtBin+2),fValueSystematicNSSigmaUp[iPtBin]*fFinalTrendNSSigma->GetBinContent(iPtBin+2));
  
            if(fPlotV2SystSeparately){
	      fv2SystematicsNSSigma->SetPoint(iPtBin,Avpt,fFinalTrendNSSigma->GetBinContent(iPtBin+2));
	      if(fSystValuev2NSSigma[iPtBin]>0.){
		fv2SystematicsNSSigma->SetPointError(iPtBin,0.9,0.9,0,fSystValuev2NSSigma[iPtBin]*fFinalTrendNSSigma->GetBinContent(iPtBin+2));	      
	      }
	      else{
		fv2SystematicsNSSigma->SetPointError(iPtBin,0.9,0.9,TMath::Abs(fSystValuev2NSSigma[iPtBin]*fFinalTrendNSSigma->GetBinContent(iPtBin+2)),0.);
	      }
               // if(fSystValuev2NSSigma[iPtBin]>=0) fv2SystematicsNSSigma[iPtBin] =new TBox(Avpt-0.9,fFinalTrendNSSigma->GetBinContent(iPtBin+2),Avpt+0.9,fFinalTrendNSSigma->GetBinContent(iPtBin+2)+fSystValuev2NSSigma[iPtBin]*fFinalTrendNSSigma->GetBinContent(iPtBin+2));
                //    else fv2SystematicsNSSigma[iPtBin] =new TBox(Avpt-0.9,fFinalTrendNSSigma->GetBinContent(iPtBin+2)-fSystValuev2NSSigma[iPtBin]*fFinalTrendNSSigma->GetBinContent(iPtBin+2),Avpt+0.9,fFinalTrendNSSigma->GetBinContent(iPtBin+2));
            }
            
        }
	fFullSystematicsNSSigma->Draw("E2");
	if(fPlotV2SystSeparately){
	  fv2SystematicsNSSigma->Draw("E2");
          Legend->Clear();
	  Legend->AddEntry(fFullSystematicsNSSigma,"Total Systematic Uncert.","l");
	  Legend->AddEntry(fv2SystematicsNSSigma,"v_{2} Systematic Uncert.","f");
	  Legend->Draw("same");
	}
        if(fPlotV2SystSeparately) fFinalTrendNSSigma->Draw("epsame");
        fCanvasFinalTrendNSSigma->Write();
    }
    if(drawASy){
      sysboxcolor = kGreen-2;
      fCanvasFinalTrendASYield = new TCanvas("CanvasFinalTrendASYield","CanvasFinalTrendASYield",0,0,1000,800);
      fCanvasFinalTrendASYield->SetTickx();
      fCanvasFinalTrendASYield->SetTicky();
      fFinalTrendASYield = new TH1D("FinalTrendASYield"," p_{T}(D) GeV/c; AS Yield",fVecSize+2,ptbinspointerout);
      SetHisto(fFinalTrendASYield,fReferenceHistoASYield);
      fCanvasFinalTrendASYield->cd();
      fFinalTrendASYield->Draw("ep");
      DefinePaveText();
      fFullSystematicsASYield->SetFillColor(sysboxcolor);
      fFullSystematicsASYield->SetFillStyle(0);
      fFullSystematicsASYield->SetLineColor(2);
      if(fPlotV2SystSeparately){
	fv2SystematicsASYield->SetFillColor(kGreen-2);
	fv2SystematicsASYield->SetFillStyle(3001);
	fv2SystematicsASYield->SetLineColor(kGreen-2);
      }
      for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
	
	Avpt = (fVecLowEdgeDpt[iPtBin]+fVecUpEdgeDpt[iPtBin])*0.5;
	fFullSystematicsASYield->SetPoint(iPtBin,Avpt,fFinalTrendASYield->GetBinContent(iPtBin+2));
	fFullSystematicsASYield->SetPointError(iPtBin,0.9,0.9,fValueSystematicASYieldLow[iPtBin]*fFinalTrendASYield->GetBinContent(iPtBin+2),fValueSystematicASYieldUp[iPtBin]*fFinalTrendASYield->GetBinContent(iPtBin+2)); 
        
	if(fPlotV2SystSeparately){
	  fv2SystematicsASYield->SetPoint(iPtBin,Avpt,fFinalTrendASYield->GetBinContent(iPtBin+2));
	  if(fSystValuev2ASYield[iPtBin]>0){
	      fv2SystematicsASYield->SetPointError(iPtBin,0.9,0.9,0.,fSystValuev2ASYield[iPtBin]*fFinalTrendASYield->GetBinContent(iPtBin+2));	  
	  }
	  else{
	    fv2SystematicsASYield->SetPointError(iPtBin,0.9,0.9,TMath::Abs(fSystValuev2ASYield[iPtBin]*fFinalTrendASYield->GetBinContent(iPtBin+2)),0.);  
	  }
	  
	  
	}
	  // if(fSystValuev2ASYield[iPtBin]>=0) fv2SystematicsNSSigma[iPtBin] =new TBox(Avpt-0.9,fFinalTrendASYield->GetBinContent(iPtBin+2),Avpt+0.9,fFinalTrendASYield->GetBinContent(iPtBin+2)+fSystValuev2ASYield[iPtBin]*fFinalTrendASYield->GetBinContent(iPtBin+2));
	  //     else fv2SystematicsASYield[iPtBin] =new TBox(Avpt-0.9,fFinalTrendASYield->GetBinContent(iPtBin+2)-fSystValuev2ASYield[iPtBin]*fFinalTrendASYield->GetBinContent(iPtBin+2),Avpt+0.9,fFinalTrendASYield->GetBinContent(iPtBin+2));
      }
      fFullSystematicsASYield->Draw("E2");
      if(fPlotV2SystSeparately){
	fv2SystematicsASYield->Draw("E2");
        Legend->Clear();
	Legend->AddEntry(fFullSystematicsASYield,"Total Systematic Uncert.","l");
	Legend->AddEntry(fv2SystematicsASYield,"v_{2} Systematic Uncert.","f");
	Legend->Draw("same");
      }
      if(fPlotV2SystSeparately) fFinalTrendASYield->Draw("epsame");
      fCanvasFinalTrendASYield->Write();
    }
    if(drawASsigma){
        sysboxcolor = kBlue-4;
        fCanvasFinalTrendASSigma = new TCanvas("CanvasFinalTrendASSigma","CanvasFinalTrendASSigma",0,0,1000,800);
        fCanvasFinalTrendASSigma->SetTickx();
        fCanvasFinalTrendASSigma->SetTicky();
        fFinalTrendASSigma = new TH1D("FinalTrendASSigma"," p_{T}(D) GeV/c; AS#sigma",fVecSize+2,ptbinspointerout);
        SetHisto(fFinalTrendASSigma,fReferenceHistoASSigma);
        fCanvasFinalTrendASSigma->cd();
        fFinalTrendASSigma->Draw("ep");
        DefinePaveText();
	fFullSystematicsASSigma->SetFillColor(sysboxcolor);
	fFullSystematicsASSigma->SetFillStyle(0);
	fFullSystematicsASSigma->SetLineColor(2);
	if(fPlotV2SystSeparately){
	  fv2SystematicsASSigma->SetFillColor(kGreen-2);
	  fv2SystematicsASSigma->SetFillStyle(3001);
	  fv2SystematicsASSigma->SetLineColor(kGreen-2);
	}	
        for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
            
            Avpt = (fVecLowEdgeDpt[iPtBin]+fVecUpEdgeDpt[iPtBin])*0.5;
	    fFullSystematicsASSigma->SetPoint(iPtBin,Avpt,fFinalTrendASSigma->GetBinContent(iPtBin+2));
	    fFullSystematicsASSigma->SetPointError(iPtBin,0.9,0.9,fValueSystematicASSigmaLow[iPtBin]*fFinalTrendASSigma->GetBinContent(iPtBin+2),fValueSystematicASSigmaUp[iPtBin]*fFinalTrendASSigma->GetBinContent(iPtBin+2));        
            
            if(fPlotV2SystSeparately){
	      fv2SystematicsASSigma->SetPoint(iPtBin,Avpt,fFinalTrendASSigma->GetBinContent(iPtBin+2));
	      if(fSystValuev2ASSigma[iPtBin]>0.){
		fv2SystematicsASSigma->SetPointError(iPtBin,0.9,0.9,0.,fSystValuev2ASSigma[iPtBin]*fFinalTrendASSigma->GetBinContent(iPtBin+2));
	      }
	      else{
		fv2SystematicsASSigma->SetPointError(iPtBin,0.9,0.9,TMath::Abs(fSystValuev2ASSigma[iPtBin]*fFinalTrendASSigma->GetBinContent(iPtBin+2)),0);
	      }
               // if(fSystValuev2ASSigma[iPtBin]>=0) fv2SystematicsNSSigma[iPtBin] =new TBox(Avpt-0.9,fFinalTrendASSigma->GetBinContent(iPtBin+2),Avpt+0.9,fFinalTrendASSigma->GetBinContent(iPtBin+2)+fSystValuev2ASSigma[iPtBin]*fFinalTrendASSigma->GetBinContent(iPtBin+2));
                 //   else fv2SystematicsASSigma[iPtBin] =new TBox(Avpt-0.9,fFinalTrendASSigma->GetBinContent(iPtBin+2)-fSystValuev2ASSigma[iPtBin]*fFinalTrendASSigma->GetBinContent(iPtBin+2),Avpt+0.9,fFinalTrendASSigma->GetBinContent(iPtBin+2));
                
                        
            }
            
        }
	fFullSystematicsASSigma->Draw("E2");
	if(fPlotV2SystSeparately){
          fv2SystematicsASSigma->Draw("E2");
          Legend->Clear();       
	  Legend->AddEntry(fFullSystematicsASSigma,"Total Systematic Uncert.","l");
	  Legend->AddEntry(fv2SystematicsASSigma,"v_{2} Systematic Uncert.","f");
	  Legend->Draw("same");
	}
        if(fPlotV2SystSeparately) fFinalTrendASSigma->Draw("epsame");
        fCanvasFinalTrendASSigma->Write();
    }Legend->Clear();
    if(drawPed){fCanvasFinalTrendPedestal = new TCanvas("CanvasFinalTrendPedestal","CanvasFinalTrendPedestal",0,0,1000,800);
        fCanvasFinalTrendPedestal->SetTickx();
        fCanvasFinalTrendPedestal->SetTicky();
        fFinalTrendPedestal = new TH1D("FinalTrendPedestal","; p_{T}(D) GeV/c; Pedestal",fVecSize+2,ptbinspointerout);
        SetHisto(fFinalTrendPedestal,fReferenceHistoPedestal);
        fCanvasFinalTrendPedestal->cd();
        fFinalTrendPedestal->Draw("ep");
        DefinePaveText();
	fFullSystematicsPedestal->SetFillColor(sysboxcolor);
	fFullSystematicsPedestal->SetFillStyle(0);
	fFullSystematicsPedestal->SetLineColor(2);
	if(fPlotV2SystSeparately){
	  fv2SystematicsPedestal->SetFillColor(kGreen-2);
	  fv2SystematicsPedestal->SetFillStyle(3001);
	  fv2SystematicsPedestal->SetLineColor(kGreen-2);
	}	
        for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
            
            Avpt = (fVecLowEdgeDpt[iPtBin]+fVecUpEdgeDpt[iPtBin])*0.5;
	    fFullSystematicsPedestal->SetPoint(iPtBin,Avpt,fFinalTrendPedestal->GetBinContent(iPtBin+2));
	    fFullSystematicsPedestal->SetPointError(iPtBin,0.9,0.9,fValueSystematicPedestalLow[iPtBin]*fFinalTrendPedestal->GetBinContent(iPtBin+2),fValueSystematicPedestalUp[iPtBin]*fFinalTrendPedestal->GetBinContent(iPtBin+2));
            if(fPlotV2SystSeparately){
	      fv2SystematicsPedestal->SetPoint(iPtBin,Avpt,fFinalTrendPedestal->GetBinContent(iPtBin+2));
	      if(fSystValuev2Pedestal[iPtBin]>0.){
		fv2SystematicsPedestal->SetPointError(iPtBin,0.9,0.9,0,fSystValuev2Pedestal[iPtBin]*fFinalTrendPedestal->GetBinContent(iPtBin+2));
	      }
	      else{
		fv2SystematicsPedestal->SetPointError(iPtBin,0.9,0.9,TMath::Abs(fSystValuev2Pedestal[iPtBin]*fFinalTrendPedestal->GetBinContent(iPtBin+2)),0.);
	      }
               // if(fSystValuev2Pedestal[iPtBin]>=0) fv2SystematicsNSSigma[iPtBin] =new TBox(Avpt-0.9,fFinalTrendPedestal->GetBinContent(iPtBin+2),Avpt+0.9,fFinalTrendPedestal->GetBinContent(iPtBin+2)+fSystValuev2Pedestal[iPtBin]*fFinalTrendPedestal->GetBinContent(iPtBin+2));
                //    else fv2SystematicsPedestal[iPtBin] =new TBox(Avpt-0.9,fFinalTrendPedestal->GetBinContent(iPtBin+2)-fSystValuev2Pedestal[iPtBin]*fFinalTrendPedestal->GetBinContent(iPtBin+2),Avpt+0.9,fFinalTrendPedestal->GetBinContent(iPtBin+2));
            }
        }
	if(fPlotV2SystSeparately){
          fv2SystematicsPedestal->Draw("E2");
          Legend->Clear();
	  Legend->AddEntry(fFullSystematicsPedestal,"Total Systematic Uncert.","l");
	  Legend->AddEntry(fv2SystematicsPedestal,"v_{2} Systematic Uncert.","f");
	  Legend->Draw("same");
	}
        if(fPlotV2SystSeparately) fFinalTrendPedestal->Draw("epsame");
        fFullSystematicsPedestal->Draw("E2");
        fCanvasFinalTrendPedestal->Write();
        
        
        
        
        /*NEW PART ADDED BY SANDRO ON 23 MAY 2015*/
        fCanvasVariationBaselineTrendPedestal = new TCanvas("CanvasBaselineVariationTrendPedestal","CanvasBaselineVariationTrendPedestal",0,0,1000,800);
        fCanvasVariationBaselineTrendPedestal->SetTickx();
        fCanvasVariationBaselineTrendPedestal->SetTicky();
        fCanvasVariationBaselineTrendPedestal->cd();
        fFinalTrendPedestal->Draw("ep");
        DefinePaveText();
        fBaselineVariationSystematicsPedestal->SetFillColor(sysboxcolor);
        fBaselineVariationSystematicsPedestal->SetFillStyle(0);
        fBaselineVariationSystematicsPedestal->SetLineColor(2);
        if(fPlotV2SystSeparately){
            fv2SystematicsPedestal->SetFillColor(kGreen-2);
            fv2SystematicsPedestal->SetFillStyle(3001);
            fv2SystematicsPedestal->SetLineColor(kGreen-2);
        }
        for(Int_t iPtBin = 0; iPtBin<fVecSize; iPtBin++){
            
            Avpt = (fVecLowEdgeDpt[iPtBin]+fVecUpEdgeDpt[iPtBin])*0.5;
            fBaselineVariationSystematicsPedestal->SetPoint(iPtBin,Avpt,fFinalTrendPedestal->GetBinContent(iPtBin+2));
            fBaselineVariationSystematicsPedestal->SetPointError(iPtBin,0.9,0.9,fValueSystematicBaselinePedestal[iPtBin]*fFinalTrendPedestal->GetBinContent(iPtBin+2),fValueSystematicBaselinePedestal[iPtBin]*fFinalTrendPedestal->GetBinContent(iPtBin+2));
            if(fPlotV2SystSeparately){
                fv2SystematicsPedestal->SetPoint(iPtBin,Avpt,fFinalTrendPedestal->GetBinContent(iPtBin+2));
                if(fSystValuev2Pedestal[iPtBin]>0.){
                    fv2SystematicsPedestal->SetPointError(iPtBin,0.9,0.9,0,fSystValuev2Pedestal[iPtBin]*fFinalTrendPedestal->GetBinContent(iPtBin+2));
                }
                else{
                    fv2SystematicsPedestal->SetPointError(iPtBin,0.9,0.9,TMath::Abs(fSystValuev2Pedestal[iPtBin]*fFinalTrendPedestal->GetBinContent(iPtBin+2)),0.);
                }
                // if(fSystValuev2Pedestal[iPtBin]>=0) fv2SystematicsNSSigma[iPtBin] =new TBox(Avpt-0.9,fFinalTrendPedestal->GetBinContent(iPtBin+2),Avpt+0.9,fFinalTrendPedestal->GetBinContent(iPtBin+2)+fSystValuev2Pedestal[iPtBin]*fFinalTrendPedestal->GetBinContent(iPtBin+2));
                //    else fv2SystematicsPedestal[iPtBin] =new TBox(Avpt-0.9,fFinalTrendPedestal->GetBinContent(iPtBin+2)-fSystValuev2Pedestal[iPtBin]*fFinalTrendPedestal->GetBinContent(iPtBin+2),Avpt+0.9,fFinalTrendPedestal->GetBinContent(iPtBin+2));
            }
        }
        fBaselineVariationSystematicsPedestal->Draw("E2");
        if(fPlotV2SystSeparately){
            Legend->Clear();
            Legend->AddEntry(fBaselineVariationSystematicsPedestal,"Baseline variation Systematic Uncert.","l");
            Legend->AddEntry(fv2SystematicsPedestal,"v_{2} Systematic Uncert.","f");
            fv2SystematicsPedestal->Draw("E2");
            Legend->Draw("same");
        }
        fCanvasVariationBaselineTrendPedestal->Write();

        
        
        
        
    } // end if draw pedestal
    
    
    
    return kTRUE;
    
}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::DefinePaveText(){
 
    TLatex *tlTitle=new TLatex(0.15,0.85,"D meson-charged hadron azimuthal correlations");
    tlTitle->SetNDC();
    tlTitle->SetTextFont(42);
    tlTitle->SetTextSize(0.03);
    tlTitle->Draw("same");
    TLatex *tlMesons= NULL;
    if(fDmeson == AliHFCorrelationUtils::kDaverage) tlMesons = new TLatex(0.17,0.8,"D^{0},D^{*+}, D^{+} average");
    if(fDmeson == AliHFCorrelationUtils::kDzero) tlMesons = new TLatex(0.17,0.8,"D^{0} #rightarrow K^{-}#pi^{+}");
    if(fDmeson == AliHFCorrelationUtils::kDstar)  tlMesons = new TLatex(0.17,0.8,"D^{*+} #rightarrow D^{0}#pi^{+} #rightarrow K^{-}#pi^{+}#pi^{+}");
    if(fDmeson == AliHFCorrelationUtils::kDplus)  tlMesons = new TLatex(0.17,0.8,"D^{+} #rightarrowK^{-} #pi^{+}#pi^{+}");
    tlMesons->SetNDC();
    tlMesons->SetTextFont(42);
    tlMesons->SetTextSize(0.03);
     tlMesons->Draw("same");
    TLatex *tlcoll;
    if(fIspPb)  tlcoll=new TLatex(0.17,0.76,"p-Pb, #sqrt{s_{NN}}=5.02 TeV");
   if(!fIspPb)  tlcoll=new TLatex(0.17,0.76,"pp, #sqrt{s}=7 TeV");
    tlcoll->SetNDC();
    tlcoll->SetTextFont(42);
    tlcoll->SetTextSize(0.03);
    tlcoll->Draw("same");
    
    
    
    TLatex *tlaspt =NULL;
    if(fAssocTrackPtMax>20 || fAssocTrackPtMax<0) tlaspt =new TLatex(0.17,0.7,Form("p^{assoc}_{T} > %.1f GeV/c",fAssocTrackPtMin));
    else tlaspt =new TLatex(0.17,0.7,Form("%.1f < p^{assoc}_{T} < %.1f GeV/c",fAssocTrackPtMin,fAssocTrackPtMax));
    tlaspt->SetNDC();
    tlaspt->SetTextFont(42);
    tlaspt->SetTextSize(0.03);
    tlaspt->Draw("same");
    
    TLatex *tleta =new TLatex(0.17,0.65,"|#Delta#eta|<1");
    tleta->SetNDC();
    tleta->SetTextFont(42);
    tleta->SetTextSize(0.03);
    tleta->Draw("same");
    
}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::SetHisto(TH1D * &outputhist, TH1D * inputhist){
    
    outputhist->SetBinContent(0,-1); outputhist->SetBinError(0,0);
   // outputhist->SetBinContent(1,-1); outputhist->SetBinError(1,0);

  //  cout << "Setting Histo " << endl;
    
    for(Int_t k = 0; k<=inputhist->GetNbinsX(); k++){
     
     //   cout << "k = " << k << "; " << inputhist->GetBinContent(k) << " set on bin center " << outputhist->GetBinCenter(k+1) <<  endl;
        
        outputhist->SetBinContent(k+1,inputhist->GetBinContent(k));
        outputhist->SetBinError(k+1,inputhist->GetBinError(k));
    }
    Int_t nbins = outputhist->GetNbinsX();
    
  
    // cout << "nbins-1 = " << nbins-1 << endl;
    outputhist->SetBinContent(nbins,-1); outputhist->SetBinError(nbins,0);
    
   // cout << "Binning check" << endl;
   // for(Int_t l =0 ; l <=nbins; l++){
   //     cout << "l = " << l << "; bin center = " << outputhist->GetBinCenter(l) << "; bin content = " << outputhist->GetBinContent(l) << "; bin error = " << outputhist->GetBinError(l) << endl;
   // }
    
   
    
    outputhist->GetYaxis()->SetRangeUser(0,2*outputhist->GetBinContent(outputhist->GetMaximumBin()));
    outputhist->SetLineColor(1);
    outputhist->SetMarkerColor(1);
    outputhist->SetMarkerStyle(20);
}



//_______________________________________________________________________________
Bool_t AliHFCorrFitSystematics::ComputeRatios(TH1D * &historef, TH1D * &histosys, TH1D * &outputhisto){
    TString r = "Ratio";
    TString histname = histosys->GetName();
    r+= histname;
    histname = r;
   // histname.Prepend("Ratio");
    outputhisto = (TH1D*)histosys->Clone(histname.Data());
    Bool_t divide = outputhisto->Divide(historef);
    if(!divide) return divide;
    return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliHFCorrFitSystematics::CheckSize(){
    
   Int_t size1 = fVecHisto.size();
   Int_t size2 = fVecHistoMinVar.size();
   Int_t size3 = fVecHistoMaxVar.size();
    

   Int_t size4 =  fVecMinCorrelatedSyst.size();
   Int_t size5 =  fVecMaxCorrelatedSyst.size();

   Int_t size6 =  fVecLowEdgeDpt.size();
   Int_t size7 =  fVecUpEdgeDpt.size();

   if(!fUseCorrelatedSystematics && fIsFittingTemplate){
     if(size4<size1){
       for(Int_t j=size4;j<size1;j++){
	 fVecMinCorrelatedSyst.push_back(0.);
       }
       size4=size1;
     }

     if(size5<size1){
       for(Int_t j=size5;j<size1;j++){
	 fVecMaxCorrelatedSyst.push_back(0.);
       }
       size5=size1;
     }
   }
    
    
    if((size1 == size2) && (size1 == size3) && (size1 == size4) && (size1 == size5) && (size1 == size6) && (size1 == size7)){
        fVecSize = size1; return kTRUE;
    }
    else{
        std::cout << "FIX ME: mismatch in sizes of one std::vectors" << std::endl;
        std::cout << "Your inputs: " << std::endl;
        std::cout << " " << std::endl;
        std::cout << "DeltaPhi mid point size = " << size1 << std::endl;
        std::cout << "DeltaPhi min variation size = " << size2 << std::endl;
        std::cout << "DeltaPhi max variation size = " << size3 << std::endl;
        std::cout << "DeltaPhi min corrsyst size = " << size4 << std::endl;
        std::cout << "DeltaPhi max corrsyst size = " << size5 << std::endl;
        std::cout << "DeltaPhi min D meson pt size = " << size6 << std::endl;
        std::cout << "DeltaPhi max D meson pt size = " << size7 << std::endl;
        
        fVecSize = -1;
        return kFALSE;
    }
}

//_______________________________________________________________________________
Bool_t AliHFCorrFitSystematics::CheckDiffSystematicsSize(){
    
    if((int)fVecSystMode.size() == 0) return kFALSE;
    
    if(fVecIsReference.size() != fVecSystMode.size()){std::cout << "Mismatch sizes" << std::endl; return kFALSE;}
    
    fVecSystModesSize = fVecSystMode.size();
    cout << "Systematics mode size = " << fVecSystModesSize << endl;
    cout << "fVecIsReference size = " << fVecIsReference.size() << endl;
    return kTRUE;
}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::SaveCanvasesDotC(){
    
    fSaveDotC = kTRUE;
    
   // TString command = "mkdir ";
   // command += outputdrectory;
   // command += "/dotC";
   // gSystem->Exec(command.Data());
        std::cout << "Saving in .C" <<std::endl;
    
    
    if(fCanvasSystematicSourcesNSYield) fCanvasSystematicSourcesNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasSystematicSourcesNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesNSSigma) fCanvasSystematicSourcesNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasSystematicSourcesNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesASYield) fCanvasSystematicSourcesASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasSystematicSourcesASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesASSigma) fCanvasSystematicSourcesASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasSystematicSourcesASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesPedestal) fCanvasSystematicSourcesPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasSystematicSourcesPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasTotalSystematicSourcesNSYield) fCanvasTotalSystematicSourcesNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesNSSigma) fCanvasTotalSystematicSourcesNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesASYield) fCanvasTotalSystematicSourcesASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesASSigma) fCanvasTotalSystematicSourcesASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesPedestal) fCanvasTotalSystematicSourcesPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasFinalTrendNSYield) fCanvasFinalTrendNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasFinalTrendNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendNSSigma) fCanvasFinalTrendNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasFinalTrendNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendASYield) fCanvasFinalTrendASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasFinalTrendASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendASSigma) fCanvasFinalTrendASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasFinalTrendASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendPedestal) fCanvasFinalTrendPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasFinalTrendPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasVariationBaselineTrendPedestal) fCanvasVariationBaselineTrendPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasVariationBaselineTrendPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    
    
    if(fCanvasRefernce) fCanvasRefernce->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasRefernce->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    for(Int_t iSystMode = 0; iSystMode<fVecSystModesSize; iSystMode++){
        
         if(fCanvasFitting[iSystMode]) fCanvasFitting[iSystMode]->SaveAs(Form("%s/%s_pthad%.1fto%.1f.C",fOutputDirectory.Data(),fCanvasFitting[iSystMode]->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    }

    cout << "Saved in .C" << endl;
 }

//_______________________________________________________________________________
void AliHFCorrFitSystematics::SaveCanvasesDotRoot(){
    
    fSaveRoot = kTRUE;
    
    // TString command = "mkdir ";
    // command += outputdrectory;
    // command += "/dotC";
    // gSystem->Exec(command.Data());
    std::cout << "Saving in .root" <<std::endl;
    
    
    if(fCanvasSystematicSourcesNSYield) fCanvasSystematicSourcesNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasSystematicSourcesNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesNSSigma) fCanvasSystematicSourcesNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasSystematicSourcesNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesASYield) fCanvasSystematicSourcesASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasSystematicSourcesASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesASSigma) fCanvasSystematicSourcesASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasSystematicSourcesASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesPedestal) fCanvasSystematicSourcesPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasSystematicSourcesPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasTotalSystematicSourcesNSYield) fCanvasTotalSystematicSourcesNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesNSSigma) fCanvasTotalSystematicSourcesNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesASYield) fCanvasTotalSystematicSourcesASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesASSigma) fCanvasTotalSystematicSourcesASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesPedestal) fCanvasTotalSystematicSourcesPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasFinalTrendNSYield) fCanvasFinalTrendNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasFinalTrendNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendNSSigma) fCanvasFinalTrendNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasFinalTrendNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendASYield) fCanvasFinalTrendASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasFinalTrendASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendASSigma) fCanvasFinalTrendASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasFinalTrendASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendPedestal) fCanvasFinalTrendPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasFinalTrendPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasVariationBaselineTrendPedestal) fCanvasVariationBaselineTrendPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasVariationBaselineTrendPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasRefernce) fCanvasRefernce->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasRefernce->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    for(Int_t iSystMode = 0; iSystMode<fVecSystModesSize; iSystMode++){
        
        if(fCanvasFitting[iSystMode]) fCanvasFitting[iSystMode]->SaveAs(Form("%s/%s_pthad%.1fto%.1f.root",fOutputDirectory.Data(),fCanvasFitting[iSystMode]->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    }
    
    cout << "Saved in .root" << endl;
}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::SaveCanvasesDotPng(){
    
    fSavePng = kTRUE;
    
    // TString command = "mkdir ";
    // command += outputdrectory;
    // command += "/dotC";
    // gSystem->Exec(command.Data());
    std::cout << "Saving in .png" <<std::endl;
    
    
    if(fCanvasSystematicSourcesNSYield) fCanvasSystematicSourcesNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasSystematicSourcesNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesNSSigma) fCanvasSystematicSourcesNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasSystematicSourcesNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesASYield) fCanvasSystematicSourcesASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasSystematicSourcesASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesASSigma) fCanvasSystematicSourcesASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasSystematicSourcesASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesPedestal) fCanvasSystematicSourcesPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasSystematicSourcesPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasTotalSystematicSourcesNSYield) fCanvasTotalSystematicSourcesNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesNSSigma) fCanvasTotalSystematicSourcesNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesASYield) fCanvasTotalSystematicSourcesASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesASSigma) fCanvasTotalSystematicSourcesASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesPedestal) fCanvasTotalSystematicSourcesPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasFinalTrendNSYield) fCanvasFinalTrendNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasFinalTrendNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendNSSigma) fCanvasFinalTrendNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasFinalTrendNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendASYield) fCanvasFinalTrendASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasFinalTrendASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendASSigma) fCanvasFinalTrendASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasFinalTrendASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendPedestal) fCanvasFinalTrendPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasFinalTrendPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasVariationBaselineTrendPedestal) fCanvasVariationBaselineTrendPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasVariationBaselineTrendPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasRefernce) fCanvasRefernce->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasRefernce->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    for(Int_t iSystMode = 0; iSystMode<fVecSystModesSize; iSystMode++){
        
        if(fCanvasFitting[iSystMode]) fCanvasFitting[iSystMode]->SaveAs(Form("%s/%s_pthad%.1fto%.1f.png",fOutputDirectory.Data(),fCanvasFitting[iSystMode]->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    }
    
    cout << "Saved in .png" << endl;
}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::SaveCanvasesDotPdf(){
    
    fSavePdf= kTRUE;
    
    // TString command = "mkdir ";
    // command += outputdrectory;
    // command += "/dotC";
    // gSystem->Exec(command.Data());
    std::cout << "Saving in .pdf" <<std::endl;
    
    
    if(fCanvasSystematicSourcesNSYield) fCanvasSystematicSourcesNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasSystematicSourcesNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesNSSigma) fCanvasSystematicSourcesNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasSystematicSourcesNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesASYield) fCanvasSystematicSourcesASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasSystematicSourcesASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesASSigma) fCanvasSystematicSourcesASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasSystematicSourcesASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesPedestal) fCanvasSystematicSourcesPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasSystematicSourcesPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasTotalSystematicSourcesNSYield) fCanvasTotalSystematicSourcesNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesNSSigma) fCanvasTotalSystematicSourcesNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesASYield) fCanvasTotalSystematicSourcesASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesASSigma) fCanvasTotalSystematicSourcesASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesPedestal) fCanvasTotalSystematicSourcesPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasFinalTrendNSYield) fCanvasFinalTrendNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasFinalTrendNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendNSSigma) fCanvasFinalTrendNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasFinalTrendNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendASYield) fCanvasFinalTrendASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasFinalTrendASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendASSigma) fCanvasFinalTrendASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasFinalTrendASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendPedestal) fCanvasFinalTrendPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasFinalTrendPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasVariationBaselineTrendPedestal) fCanvasVariationBaselineTrendPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasVariationBaselineTrendPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasRefernce) fCanvasRefernce->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasRefernce->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    for(Int_t iSystMode = 0; iSystMode<fVecSystModesSize; iSystMode++){
        
        if(fCanvasFitting[iSystMode]) fCanvasFitting[iSystMode]->SaveAs(Form("%s/%s_pthad%.1fto%.1f.pdf",fOutputDirectory.Data(),fCanvasFitting[iSystMode]->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    }
    
    cout << "Saved in .pdf" << endl;
}

//_______________________________________________________________________________
void AliHFCorrFitSystematics::SaveCanvasesDotEps(){
    
    fSaveEps = kTRUE;
    
    // TString command = "mkdir ";
    // command += outputdrectory;
    // command += "/dotC";
    // gSystem->Exec(command.Data());
    std::cout << "Saving in .eps" <<std::endl;
    
    
    if(fCanvasSystematicSourcesNSYield) fCanvasSystematicSourcesNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasSystematicSourcesNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesNSSigma) fCanvasSystematicSourcesNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasSystematicSourcesNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesASYield) fCanvasSystematicSourcesASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasSystematicSourcesASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesASSigma) fCanvasSystematicSourcesASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasSystematicSourcesASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasSystematicSourcesPedestal) fCanvasSystematicSourcesPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasSystematicSourcesPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasTotalSystematicSourcesNSYield) fCanvasTotalSystematicSourcesNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesNSSigma) fCanvasTotalSystematicSourcesNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesASYield) fCanvasTotalSystematicSourcesASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesASSigma) fCanvasTotalSystematicSourcesASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasTotalSystematicSourcesPedestal) fCanvasTotalSystematicSourcesPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasTotalSystematicSourcesPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasFinalTrendNSYield) fCanvasFinalTrendNSYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasFinalTrendNSYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendNSSigma) fCanvasFinalTrendNSSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasFinalTrendNSSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendASYield) fCanvasFinalTrendASYield->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasFinalTrendASYield->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendASSigma) fCanvasFinalTrendASSigma->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasFinalTrendASSigma->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasFinalTrendPedestal) fCanvasFinalTrendPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasFinalTrendPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    if(fCanvasVariationBaselineTrendPedestal) fCanvasVariationBaselineTrendPedestal->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasVariationBaselineTrendPedestal->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    if(fCanvasRefernce) fCanvasRefernce->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasRefernce->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    
    for(Int_t iSystMode = 0; iSystMode<fVecSystModesSize; iSystMode++){
        
        if(fCanvasFitting[iSystMode]) fCanvasFitting[iSystMode]->SaveAs(Form("%s/%s_pthad%.1fto%.1f.eps",fOutputDirectory.Data(),fCanvasFitting[iSystMode]->GetName(),fAssocTrackPtMin,fAssocTrackPtMax));
    }
    
    cout << "Saved in .eps" << endl;
}


