/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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


///////////////////////////////////////////////////////////////////////////
/*

Origin: marian.ivanov@cern.ch
Frequenlty used function for visualization 
marian.ivanov@cern.ch
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string.h>
//ROOT includes
#include "TROOT.h"
#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TString.h"
#include "TBenchmark.h"
#include "TStopwatch.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TTimer.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TView.h"
#include "TView3D.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TObjString.h"


//ALIROOT includes
#include "AliTrackPointArray.h"
#include "AliTreeDraw.h" 

#endif

//
//     Class for visualization and some statistacal analysis using tree
//     To be used in comparisons
//                and calib viewers based on tree    


ClassImp(AliTreeDraw)


AliTreeDraw::AliTreeDraw():
  fTree(0),
  fRes(0),
  fMean(0),
  fPoints(0){
  //
  // default constructor
  //
}

void  AliTreeDraw::ClearHisto(){
  //
  //
  delete fRes; 
  delete fMean;
  fRes=0;
  fMean=0;
}



TH1F * AliTreeDraw::DrawXY(const char * chx, const char *chy, const char* selection, 
		const char * quality, Int_t nbins, Float_t minx, Float_t maxx, Float_t miny, Float_t maxy, Int_t nBinsRes)
{
  //
  Double_t* bins = CreateLogBins(nbins, minx, maxx);
  TH2F* hRes2 = new TH2F("hRes2", "residuals", nbins, minx, maxx, nBinsRes, miny, maxy);
  char cut[1000];
  snprintf(cut,1000,"%s&&%s",selection,quality);
  char expression[1000];
  snprintf(expression,1000,"%s:%s>>hRes2",chy,chx);
  fTree->Draw(expression, cut, "groff");
  TH1F* hMean=0;
  TH1F* hRes = CreateResHisto(hRes2, &hMean);
  AliLabelAxes(hRes, chx, chy);
  //
  delete hRes2;
  delete[] bins;
  ClearHisto();
  fRes  = hRes;
  fMean = hMean;
  return hRes;
}



TH1F * AliTreeDraw::DrawLogXY(const char * chx, const char *chy, const char* selection, 
				    const char * quality, Int_t nbins, Float_t minx, Float_t maxx, Float_t miny, Float_t maxy, Int_t nBinsRes)
{
  //
  // 
  //
  Double_t* bins = CreateLogBins(nbins, minx, maxx);
  TH2F* hRes2 = new TH2F("hRes2", "residuals", nbins, bins, nBinsRes, miny, maxy);
  char cut[1000];
  snprintf(cut,1000,"%s&&%s",selection,quality);
  char expression[1000];
  snprintf(expression,1000,"%s:%s>>hRes2",chy,chx);
  fTree->Draw(expression, cut, "groff");
  TH1F* hMean=0;  
  TH1F* hRes = CreateResHisto(hRes2, &hMean);
  AliLabelAxes(hRes, chx, chy);
  //
  delete hRes2;
  delete[] bins;
  ClearHisto();
  fRes  = hRes;
  fMean = hMean;
  return hRes;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
TH1F * AliTreeDraw::Eff(const char *variable, const char* selection, const char * quality, 
			      Int_t nbins, Float_t min, Float_t max)
{
  //
  //
  TH1F* hGen = new TH1F("hGen", "gen. tracks", nbins, min, max);
  TH1F* hRec = new TH1F("hRec", "rec. tracks", nbins, min, max);
  char inputGen[1000];  
  snprintf(inputGen,1000,"%s>>hGen", variable);
  fTree->Draw(inputGen, selection, "groff");
  char selectionRec[256];
  snprintf(selectionRec,1000, "%s && %s", selection, quality);
  char inputRec[1000];  
  snprintf(inputRec,1000,"%s>>hRec", variable);
  fTree->Draw(inputRec, selectionRec, "groff");
  //
  TH1F* hEff = CreateEffHisto(hGen, hRec);
  AliLabelAxes(hEff, variable, "#epsilon [%]");
  fRes = hEff;
  delete hRec;
  delete hGen;
  return hEff;
}



///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
TH1F * AliTreeDraw::EffLog(const char *variable, const char* selection, const char * quality, 
			      Int_t nbins, Float_t min, Float_t max)
{
  //
  //
  Double_t* bins = CreateLogBins(nbins, min, max);
  TH1F* hGen = new TH1F("hGen", "gen. tracks", nbins, bins);
  TH1F* hRec = new TH1F("hRec", "rec. tracks", nbins, bins);
  char inputGen[1000];  
  snprintf(inputGen,1000,"%s>>hGen", variable);
  fTree->Draw(inputGen, selection, "groff");
  char selectionRec[256];
  snprintf(selectionRec,1000, "%s && %s", selection, quality);
  char inputRec[1000];  
  snprintf(inputRec,1000,"%s>>hRec", variable);
  fTree->Draw(inputRec, selectionRec, "groff");
  //
  TH1F* hEff = CreateEffHisto(hGen, hRec);
  AliLabelAxes(hEff, variable, "#epsilon [%]");
  fRes = hEff;
  delete hRec;
  delete hGen;
  delete[] bins;
  return hEff;
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

Double_t* AliTreeDraw::CreateLogBins(Int_t nBins, Double_t xMin, Double_t xMax)
{
  Double_t* bins = new Double_t[nBins+1];
  bins[0] = xMin;
  Double_t factor = pow(xMax/xMin, 1./nBins);
  for (Int_t i = 1; i <= nBins; i++)
    bins[i] = factor * bins[i-1];
  return bins;
}




void AliTreeDraw::AliLabelAxes(TH1* histo, const char* xAxisTitle, const char* yAxisTitle)
{
  //
  histo->GetXaxis()->SetTitle(xAxisTitle);
  histo->GetYaxis()->SetTitle(yAxisTitle);
}


TH1F* AliTreeDraw::CreateEffHisto(TH1F* hGen, TH1F* hRec)
{
  //
  Int_t nBins = hGen->GetNbinsX();
  TH1F* hEff = (TH1F*) hGen->Clone("hEff");
  hEff->SetTitle("");
  hEff->SetStats(kFALSE);
  hEff->SetMinimum(0.);
  hEff->SetMaximum(110.);
  //
  for (Int_t iBin = 0; iBin <= nBins; iBin++) {
    Double_t nGen = hGen->GetBinContent(iBin);
    Double_t nRec = hRec->GetBinContent(iBin);
    if (nGen > 0) {
      Double_t eff = nRec/nGen;
      hEff->SetBinContent(iBin, 100. * eff);
      Double_t error = sqrt((eff*(1.-eff)+0.01) / nGen);      
      //      if (error == 0) error = sqrt(0.1/nGen);
      //
      if (error == 0) error = 0.0001;
      hEff->SetBinError(iBin, 100. * error);
    } else {
      hEff->SetBinContent(iBin, 100. * 0.5);
      hEff->SetBinError(iBin, 100. * 0.5);
    }
  }
  return hEff;
}



TH1F* AliTreeDraw::CreateResHisto(TH2F* hRes2, TH1F **phMean,  Bool_t drawBinFits, 
		     Bool_t overflowBinFits)
{
  TVirtualPad* currentPad = gPad;
  TAxis* axis = hRes2->GetXaxis();
  Int_t nBins = axis->GetNbins();
  TH1F* hRes, *hMean;
  if (axis->GetXbins()->GetSize()){
    hRes = new TH1F("hRes", "", nBins, axis->GetXbins()->GetArray());
    hMean = new TH1F("hMean", "", nBins, axis->GetXbins()->GetArray());
  }
  else{
    hRes = new TH1F("hRes", "", nBins, axis->GetXmin(), axis->GetXmax());
    hMean = new TH1F("hMean", "", nBins, axis->GetXmin(), axis->GetXmax());

  }
  hRes->SetStats(false);
  hRes->SetOption("E");
  hRes->SetMinimum(0.);
  //
  hMean->SetStats(false);
  hMean->SetOption("E");
 
  // create the fit function
  TF1 * fitFunc = new TF1("G","[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",-3,3);
  
  fitFunc->SetLineWidth(2);
  fitFunc->SetFillStyle(0);
  // create canvas for fits
  TCanvas* canBinFits = NULL;
  Int_t nPads = (overflowBinFits) ? nBins+2 : nBins;
  Int_t nx = Int_t(sqrt(nPads-1.));// + 1;
  Int_t ny = (nPads-1) / nx + 1;
  if (drawBinFits) {
    canBinFits = (TCanvas*)gROOT->FindObject("canBinFits");
    if (canBinFits) delete canBinFits;
    canBinFits = new TCanvas("canBinFits", "fits of bins", 200, 100, 500, 700);
    canBinFits->Divide(nx, ny);
  }

  // loop over x bins and fit projection
  Int_t dBin = ((overflowBinFits) ? 1 : 0);
  for (Int_t bin = 1-dBin; bin <= nBins+dBin; bin++) {
    if (drawBinFits) canBinFits->cd(bin + dBin);
    TH1D* hBin = hRes2->ProjectionY("hBin", bin, bin);
    //    
    if (hBin->GetEntries() > 5) {
      fitFunc->SetParameters(hBin->GetMaximum(),hBin->GetMean(),hBin->GetRMS());
      hBin->Fit(fitFunc,"s");
      Double_t sigma = TMath::Abs(fitFunc->GetParameter(2));

      if (sigma > 0.){
	hRes->SetBinContent(bin, TMath::Abs(fitFunc->GetParameter(2)));
	hMean->SetBinContent(bin, fitFunc->GetParameter(1));	
      }
      else{
	hRes->SetBinContent(bin, 0.);
	hMean->SetBinContent(bin,0);
      }
      hRes->SetBinError(bin, fitFunc->GetParError(2));
      hMean->SetBinError(bin, fitFunc->GetParError(1));
      
      //
      //

    } else {
      hRes->SetBinContent(bin, 0.);
      hRes->SetBinError(bin, 0.);
      hMean->SetBinContent(bin, 0.);
      hMean->SetBinError(bin, 0.);
    }
    

    if (drawBinFits) {
      char name[256];
      if (bin == 0) {
	snprintf(name,256, "%s < %.4g", axis->GetTitle(), axis->GetBinUpEdge(bin));
      } else if (bin == nBins+1) {
	snprintf(name,256, "%.4g < %s", axis->GetBinLowEdge(bin), axis->GetTitle());
      } else {
	snprintf(name,256, "%.4g < %s < %.4g", axis->GetBinLowEdge(bin),
		axis->GetTitle(), axis->GetBinUpEdge(bin));
      }
      canBinFits->cd(bin + dBin);
      hBin->SetTitle(name);
      hBin->SetStats(kTRUE);
      hBin->DrawCopy("E");
      canBinFits->Update();
      canBinFits->Modified();
      canBinFits->Update();
    }
    
    delete hBin;
  }

  delete fitFunc;
  currentPad->cd();
  *phMean = hMean;
  return hRes;
}

TH1F* AliTreeDraw::CreateResHistoI(TH2F* hRes2, TH1F **phMean, Int_t integ,  Bool_t drawBinFits)
{
  TVirtualPad* currentPad = gPad;
  TAxis* axis = hRes2->GetXaxis();
  Int_t nBins = axis->GetNbins();
  //Bool_t overflowBinFits = kFALSE;
  
  TH1F* hRes, *hMean;
  if (axis->GetXbins()->GetSize()){
    hRes = new TH1F("hRes", "", nBins, axis->GetXbins()->GetArray());
    hMean = new TH1F("hMean", "", nBins, axis->GetXbins()->GetArray());
  }
  else{
    hRes = new TH1F("hRes", "", nBins, axis->GetXmin(), axis->GetXmax());
    hMean = new TH1F("hMean", "", nBins, axis->GetXmin(), axis->GetXmax());

  }
  hRes->SetStats(false);
  hRes->SetOption("E");
  hRes->SetMinimum(0.);
  //
  hMean->SetStats(false);
  hMean->SetOption("E");
 
  // create the fit function
  TF1 * fitFunc = new TF1("G","[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",-3,3);
  
  fitFunc->SetLineWidth(2);
  fitFunc->SetFillStyle(0);
  // create canvas for fits
  TCanvas* canBinFits = NULL;
  //Int_t nPads = (overflowBinFits) ? nBins+2 : nBins;
  Int_t nPads = nBins;
  Int_t nx = Int_t(sqrt(nPads-1.));// + 1;
  Int_t ny = (nPads-1) / nx + 1;
  if (drawBinFits) {
    canBinFits = (TCanvas*)gROOT->FindObject("canBinFits");
    if (canBinFits) delete canBinFits;
    canBinFits = new TCanvas("canBinFits", "fits of bins", 200, 100, 500, 700);
    canBinFits->Divide(nx, ny);
  }

  // loop over x bins and fit projection
  //Int_t dBin = ((overflowBinFits) ? 1 : 0);
  Int_t dBin =  0;
  for (Int_t bin = 1-dBin; bin <= nBins+dBin; bin++) {
    if (drawBinFits) canBinFits->cd(bin + dBin);
    Int_t bin0=TMath::Max(bin-integ,0);
    Int_t bin1=TMath::Min(bin+integ,nBins);
    TH1D* hBin = hRes2->ProjectionY("hBin", bin0, bin1);
    //    
    if (hBin->GetEntries() > 5) {
      fitFunc->SetParameters(hBin->GetMaximum(),hBin->GetMean(),hBin->GetRMS());
      hBin->Fit(fitFunc,"s");
      Double_t sigma = TMath::Abs(fitFunc->GetParameter(2));

      if (sigma > 0.){
	hRes->SetBinContent(bin, TMath::Abs(fitFunc->GetParameter(2)));
	hMean->SetBinContent(bin, fitFunc->GetParameter(1));	
      }
      else{
	hRes->SetBinContent(bin, 0.);
	hMean->SetBinContent(bin,0);
      }
      hRes->SetBinError(bin, fitFunc->GetParError(2));
      hMean->SetBinError(bin, fitFunc->GetParError(1));
      
      //
      //

    } else {
      hRes->SetBinContent(bin, 0.);
      hRes->SetBinError(bin, 0.);
      hMean->SetBinContent(bin, 0.);
      hMean->SetBinError(bin, 0.);
    }
    

    if (drawBinFits) {
      char name[256];
      if (bin == 0) {
	snprintf(name,256, "%s < %.4g", axis->GetTitle(), axis->GetBinUpEdge(bin));
      } else if (bin == nBins+1) {
	snprintf(name,256, "%.4g < %s", axis->GetBinLowEdge(bin), axis->GetTitle());
      } else {
	snprintf(name,256, "%.4g < %s < %.4g", axis->GetBinLowEdge(bin),
		axis->GetTitle(), axis->GetBinUpEdge(bin));
      }
      canBinFits->cd(bin + dBin);
      hBin->SetTitle(name);
      hBin->SetStats(kTRUE);
      hBin->DrawCopy("E");
      canBinFits->Update();
      canBinFits->Modified();
      canBinFits->Update();
    }
    
    delete hBin;
  }

  delete fitFunc;
  currentPad->cd();
  *phMean = hMean;
  return hRes;
}

TH1F* AliTreeDraw::CreateResHistoII(TH2F* hRes2, TH1F **phMean, Int_t integ,  Bool_t drawBinFits, Int_t cut)
{
  TVirtualPad* currentPad = gPad;
  TAxis* axis = hRes2->GetXaxis();
  Int_t nBins = axis->GetNbins();
  //Bool_t overflowBinFits = kFALSE;
  TH1F* hRes, *hMean;
  if (axis->GetXbins()->GetSize()){
    hRes = new TH1F("hRes", "", nBins, axis->GetXbins()->GetArray());
    hMean = new TH1F("hMean", "", nBins, axis->GetXbins()->GetArray());
  }
  else{
    hRes = new TH1F("hRes", "", nBins, axis->GetXmin(), axis->GetXmax());
    hMean = new TH1F("hMean", "", nBins, axis->GetXmin(), axis->GetXmax());

  }
  hRes->SetStats(false);
  hRes->SetOption("E");
  hRes->SetMinimum(0.);
  //
  hMean->SetStats(false);
  hMean->SetOption("E");
 
  // create the fit function
  TF1 * fitFunc = new TF1("G","[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",-3,3);
  
  fitFunc->SetLineWidth(2);
  fitFunc->SetFillStyle(0);
  // create canvas for fits
  TCanvas* canBinFits = NULL;
  //Int_t nPads = (overflowBinFits) ? nBins+2 : nBins;
  Int_t nPads = nBins;
  Int_t nx = Int_t(sqrt(nPads-1.));// + 1;
  Int_t ny = (nPads-1) / nx + 1;
  if (drawBinFits) {
    canBinFits = (TCanvas*)gROOT->FindObject("canBinFits");
    if (canBinFits) delete canBinFits;
    canBinFits = new TCanvas("canBinFits", "fits of bins", 200, 100, 500, 700);
    canBinFits->Divide(nx, ny);
  }

  // loop over x bins and fit projection
  //Int_t dBin = ((overflowBinFits) ? 1 : 0);
  Int_t dBin = 0;
  for (Int_t bin = 1-dBin; bin <= nBins+dBin; bin++) {
    if (drawBinFits) canBinFits->cd(bin + dBin);
    Int_t bin0=TMath::Max(bin-integ,0);
    Int_t bin1=TMath::Min(bin+integ,nBins);
    TH1D* hBin = hRes2->ProjectionY("hBin", bin0, bin1);
    //    
    if (hBin->GetEntries() > cut) {
      fitFunc->SetParameters(hBin->GetMaximum(),hBin->GetMean(),hBin->GetRMS());
      hBin->Fit(fitFunc,"s");
      Double_t sigma = TMath::Abs(fitFunc->GetParameter(2));

      if (sigma > 0.){
	hRes->SetBinContent(bin, TMath::Abs(fitFunc->GetParameter(2)));
	hMean->SetBinContent(bin, fitFunc->GetParameter(1));	
      }
      else{
	hRes->SetBinContent(bin, 0.);
	hMean->SetBinContent(bin,0);
      }
      hRes->SetBinError(bin, fitFunc->GetParError(2));
      hMean->SetBinError(bin, fitFunc->GetParError(1));
      
      //
      //

    } else {
      hRes->SetBinContent(bin, 0.);
      hRes->SetBinError(bin, 0.);
      hMean->SetBinContent(bin, 0.);
      hMean->SetBinError(bin, 0.);
    }
    

    if (drawBinFits) {
      char name[256];
      if (bin == 0) {
	snprintf(name,256, "%s < %.4g", axis->GetTitle(), axis->GetBinUpEdge(bin));
      } else if (bin == nBins+1) {
	snprintf(name,256, "%.4g < %s", axis->GetBinLowEdge(bin), axis->GetTitle());
      } else {
	snprintf(name,256, "%.4g < %s < %.4g", axis->GetBinLowEdge(bin),
		axis->GetTitle(), axis->GetBinUpEdge(bin));
      }
      canBinFits->cd(bin + dBin);
      hBin->SetTitle(name);
      hBin->SetStats(kTRUE);
      hBin->DrawCopy("E");
      canBinFits->Update();
      canBinFits->Modified();
      canBinFits->Update();
    }
    
    delete hBin;
  }

  delete fitFunc;
  currentPad->cd();
  *phMean = hMean;
  return hRes;
}




void   AliTreeDraw::GetPoints3D(const char * label, const char * chpoints, 
				const char* selection, TTree * tpoints, Int_t color,Float_t rmin){
  //
  // load selected points from tree
  //
   if (!fPoints) fPoints = new TObjArray;
   if (tpoints->GetIndex()==0) tpoints->BuildIndex("fLabel","Label");
   TBranch * br = tpoints->GetBranch(chpoints);
   if (!br) return;
   AliTrackPointArray * points = new AliTrackPointArray;
   br->SetAddress(&points);

   Int_t npoints = fTree->Draw(label,selection);
   Float_t xyz[30000];
   rmin*=rmin;
   for (Int_t ii=0;ii<npoints;ii++){
     Int_t index = (Int_t)fTree->GetV1()[ii];
     tpoints->GetEntryWithIndex(index,index);
     if (points->GetNPoints()<2) continue;
     Int_t goodpoints=0;
     for (Int_t i=0;i<points->GetNPoints(); i++){
       xyz[goodpoints*3]   = points->GetX()[i];
       xyz[goodpoints*3+1] = points->GetY()[i];
       xyz[goodpoints*3+2] = points->GetZ()[i];
       if ( points->GetX()[i]*points->GetX()[i]+points->GetY()[i]*points->GetY()[i]>rmin) goodpoints++;
     } 
     TPolyMarker3D * marker = new TPolyMarker3D(goodpoints,xyz); 
     marker->SetMarkerColor(color);
     marker->SetMarkerStyle(1);
     fPoints->AddLast(marker);
   }
}




TString* AliTreeDraw::FitPlane(const char* drawCommand, const char* formula, const char* cuts, Double_t & chi2, TVectorD &fitParam, TMatrixD &covMatrix, Int_t start, Int_t stop){
   //
   // fit an arbitrary function, specified by formula into the data, specified by drawCommand and cuts
   // returns chi2, fitParam and covMatrix
   // returns TString with fitted formula
   //
    
   TString formulaStr(formula); 
   TString drawStr(drawCommand);
   TString cutStr(cuts);
      
   formulaStr.ReplaceAll("++", "~");
   TObjArray* formulaTokens = formulaStr.Tokenize("~"); 
   Int_t dim = formulaTokens->GetEntriesFast();
   
   fitParam.ResizeTo(dim);
   covMatrix.ResizeTo(dim,dim);
   
   TLinearFitter* fitter = new TLinearFitter(dim+1, Form("hyp%d",dim));
   fitter->StoreData(kTRUE);   
   fitter->ClearPoints();
   
   Int_t entries = fTree->Draw(drawStr.Data(), cutStr.Data(), "goff",  stop-start, start);
   if (entries == -1) return new TString("An ERROR has occured during fitting!");
   Double_t **values = new Double_t*[dim+1] ; 
   
   for (Int_t i = 0; i < dim + 1; i++){
      Int_t centries = 0;
      if (i < dim) centries = fTree->Draw(((TObjString*)formulaTokens->At(i))->GetName(), cutStr.Data(), "goff", stop-start,start);
      else  centries = fTree->Draw(drawStr.Data(), cutStr.Data(), "goff", stop-start,start);
      
      if (entries != centries) { 
        return new TString("An ERROR has occured during fitting!");
      }
      values[i] = new Double_t[entries];
      memcpy(values[i],  fTree->GetV1(), entries*sizeof(Double_t)); 
   }
   
   // add points to the fitter
   for (Int_t i = 0; i < entries; i++){
      Double_t x[1000];
      for (Int_t j=0; j<dim;j++) x[j]=values[j][i];
      fitter->AddPoint(x, values[dim][i], 1);
   }

   fitter->Eval();
   fitter->GetParameters(fitParam);
   fitter->GetCovarianceMatrix(covMatrix);
   chi2 = fitter->GetChisquare();
   
   TString *preturnFormula = new TString(Form("( %f+",fitParam[0])), &returnFormula = *preturnFormula; 
   
   for (Int_t iparam = 0; iparam < dim; iparam++) {
     returnFormula.Append(Form("%s*(%f)",((TObjString*)formulaTokens->At(iparam))->GetName(),fitParam[iparam+1]));
     if (iparam < dim-1) returnFormula.Append("+");
   }
   returnFormula.Append(" )");
   delete formulaTokens;
   delete fitter;
   delete[] *values;
   delete[] values;
   return preturnFormula;
}

