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

/*

//0.  Libraries to lod
gSystem->Load("libANALYSIS");
gSystem->Load("libSTAT");
gSystem->Load("libTPCcalib");


//1. Do calibration ...
//
// compare reference

//
//2. Visulaize results
//
TFile fcalib("CalibObjects.root");
TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
AliTPCcalibTimeGain * gain = ( AliTPCcalibTimeGain *)array->FindObject("calibTimeGain");
TGraph * gr = AliTPCcalibBase::FitSlices(gain->GetHistGainTime(),0,1,2000,10);
gr->Draw("ALPsame")


//
// MakeSlineFit example
//
AliSplineFit fit;
fit.SetGraph(gr)
fit->SetMinPoints(gr->GetN()+1);
fit->InitKnots(gr,2,0,0.001)
fit.SplineFit(0)
fit.MakeDiffHisto(gr)->Draw();
TGraph * grfit = fit.MakeGraph(gr->GetX()[0],gr->GetX()[gr->GetN()-1],50000,0);

gr->SetMarkerStyle(25);
gr->Draw("alp");
grfit->SetLineColor(2);
grfit->Draw("lu");

*/


#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"

#include "AliTracker.h"
#include "AliMagF.h"
#include "AliTPCCalROC.h"

#include "AliLog.h"

#include "AliTPCcalibTimeGain.h"

#include "TTreeStream.h"
#include "AliTPCTracklet.h"
#include "TTimeStamp.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibLaser.h"
#include "AliDCSSensorArray.h"
#include "AliDCSSensor.h"

ClassImp(AliTPCcalibTimeGain)


AliTPCcalibTimeGain::AliTPCcalibTimeGain() 
  :AliTPCcalibBase(), 
   fHistGainTime(0),
   fGainVsTime(0),
   fHistDeDxTotal(0),
   fIntegrationTimeDeDx(0),
   fMIP(0),
   fLowerTrunc(0),
   fUpperTrunc(0),
   fUseShapeNorm(0),
   fUsePosNorm(0),
   fUsePadNorm(0),
   fIsCosmic(0)
{  
  AliInfo("Default Constructor");  
}


AliTPCcalibTimeGain::AliTPCcalibTimeGain(const Text_t *name, const Text_t *title, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeGain)
  :AliTPCcalibBase(),
   fHistGainTime(0),
   fGainVsTime(0),
   fHistDeDxTotal(0),
   fIntegrationTimeDeDx(0),
   fMIP(0),
   fLowerTrunc(0),
   fUpperTrunc(0),
   fUseShapeNorm(0),
   fUsePosNorm(0),
   fUsePadNorm(0),
   fIsCosmic(0)
 {
  
  SetName(name);
  SetTitle(title);

  AliInfo("Non Default Constructor");

  fIntegrationTimeDeDx = deltaIntegrationTimeGain;
 
  Double_t deltaTime = EndTime - StartTime;
  

  // main histogram for time dependence: dE/dx, time, type (0-all,1-muon cosmic,2-pion beam data), meanDriftlength, momenta (only filled if enough space is available)
  Int_t binsGainTime[5] = {100, TMath::Nint(deltaTime/deltaIntegrationTimeGain), 3, 25, 200};
  Double_t xminGainTime[5] = {0.5, StartTime, -0.5,   0, 0.1};
  Double_t xmaxGainTime[5] = {  4,   EndTime,  2.5, 250, 50};
  fHistGainTime = new THnSparseF("HistGainTime","dEdx l;time;dEdx",5,binsGainTime,xminGainTime,xmaxGainTime);
  BinLogX(fHistGainTime, 4);
  //
  fHistDeDxTotal = new TH2F("DeDx","dEdx; momentum p (GeV); TPC signal (a.u.)",250,0.01,100.,1000,0.,1000);
  BinLogX(fHistDeDxTotal);
  
  // default values for dE/dx
  fMIP = 50.;
  fLowerTrunc = 0.0;
  fUpperTrunc = 0.7;
  fUseShapeNorm = kTRUE;
  fUsePosNorm = kFALSE;
  fUsePadNorm = kFALSE;
  //
  fIsCosmic = kTRUE;
  //
  
 }



AliTPCcalibTimeGain::~AliTPCcalibTimeGain(){
  //
  //
  //
}


void AliTPCcalibTimeGain::Process(AliESDEvent *event) {
  //
  // main track loop
  //
  if (!event) {
    Printf("ERROR: ESD not available");
    return;
  }  

  Int_t ntracks=event->GetNumberOfTracks();
  
  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!ESDfriend) {
   Printf("ERROR: ESDfriend not available");
   return;
  }
  //
  UInt_t time = event->GetTimeStamp();
  if (time < 0.1) time = (UInt_t)(fHistGainTime->GetAxis(1)->GetXmin() + 1);
  //
  // track loop
  //
  for (Int_t i=0;i<ntracks;++i) {

    AliESDtrack *track = event->GetTrack(i);
    if (!track) continue;
        
    const AliExternalTrackParam * trackIn = track->GetInnerParam();
    const AliExternalTrackParam * trackOut = track->GetOuterParam();
    if (!trackIn) continue;
    if (!trackOut) continue;

    // calculate necessary track parameters
    //Double_t meanP = 0.5*(trackIn->GetP() + trackOut->GetP());
    Double_t meanP = trackIn->GetP();
    Double_t meanDrift = 250 - 0.5*TMath::Abs(trackIn->GetZ() + trackOut->GetZ());
    //Double_t d = trackIn->GetLinearD(0,0);
    Int_t NclsDeDx = track->GetTPCNcls();

    //if (meanP > 0.7 || meanP < 0.2) continue;
    if (fIsCosmic && meanP < 20) continue;
    if (NclsDeDx < 60) continue;     

    // exclude tracks which do not look like primaries or are simply too short or on wrong sectors

    //if (TMath::Abs(trackIn->GetSnp()) > 3*0.4) continue;
    //if (TMath::Abs(trackIn->GetZ()) > 150) continue;   
    //if (seed->CookShape(1) > 1) continue;
    //if (TMath::Abs(trackIn->GetY()) > 20) continue;
    //if (TMath::Abs(d)>20) continue;   // distance to the 0,0; select only tracks which cross chambers under proper angle
    if (TMath::Abs(trackIn->GetTgl()) > 1) continue;
    if (TMath::Abs(trackIn->GetSnp()) > 0.6) continue;
    
    // Get seeds
    AliESDfriendTrack *friendTrack = ESDfriend->GetTrack(i);
    if (!friendTrack) continue;
    TObject *calibObject;
    AliTPCseed *seed = 0;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }    

    if (seed) { 
      Double_t TPCsignalMax = (1/fMIP)*seed->CookdEdxNorm(fLowerTrunc,fUpperTrunc,1,0,159,fUseShapeNorm,fUsePosNorm,fUsePadNorm,0);
      //Double_t TPCsignalMax = (1/fMIP)*track->GetTPCsignal();
      fHistDeDxTotal->Fill(meanP, TPCsignalMax);

      //
      //dE/dx, time, type (0-all,1-muon cosmic,2-pion beam data), momenta
      Double_t vec[5] = {TPCsignalMax,time,0,meanDrift,meanP}; 
      fHistGainTime->Fill(vec); // avoid this filling if memory consumption is too high

      // only partial filling if memory consumption has to be kept low; for cosmic and beam data
      if (fIsCosmic) {
	Double_t vecCos[5] = {TPCsignalMax,time,1,meanDrift,20}; 
	if (meanP > 20) fHistGainTime->Fill(vecCos);
      } else {
	Double_t vecBeam[5] = {TPCsignalMax,time,2,meanDrift,0.5};
	if (meanP > 0.4 && meanP < 0.5) fHistGainTime->Fill(vecBeam);
      }
    

    } else {
      cout << "ERROR: TPC seed not found" << endl;
    }
    
  }
  
  
}


void AliTPCcalibTimeGain::Analyze() {
  //
  //
  //
  TObjArray arr;
  if (fIsCosmic) {
    fHistGainTime->GetAxis(2)->SetRangeUser(0.51,1.49);
  } else {
    fHistGainTime->GetAxis(2)->SetRangeUser(1.51,2.49);
  }
  fHistGainTime->Projection(0,1)->FitSlicesY(0,0,-1,0,"QNR",&arr);
  TH1D * fitMean = (TH1D*) arr.At(1);
  //
  fGainVsTime = new TGraph(fitMean);
  //
  return;
}


Long64_t AliTPCcalibTimeGain::Merge(TCollection *li) {

  TIterator* iter = li->MakeIterator();
  AliTPCcalibTimeGain* cal = 0;

  while ((cal = (AliTPCcalibTimeGain*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibTimeGain::Class())) {
      Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }

    // add histograms here...
    if (cal->GetHistGainTime()) fHistGainTime->Add(cal->GetHistGainTime());
    if (cal->GetHistDeDxTotal()) fHistDeDxTotal->Add(cal->GetHistDeDxTotal());

  }
  
  return 0;
  
}


TGraph * AliTPCcalibTimeGain::FitSlices(THnSparse *h, Int_t axisDim1, Int_t axisDim2, Int_t minEntries){
  //
  // Fitting slices of the projection(axisDim1,axisDim2) of a sparse histogram
  //
  TH2D * hist = h->Projection(axisDim1, axisDim2);

  Double_t xvec[10000];
  Double_t yvec[10000];
  Int_t counter = 0;

  for(Int_t i=1; i < hist->GetNbinsX(); i++) {
    Int_t interval = 0;
    if (hist->Integral(i,i,0,hist->GetNbinsY()) < minEntries) {
      if (hist->Integral(i,i+1,0,hist->GetNbinsY()) < minEntries) {
	if (hist->Integral(i,i+2,0,hist->GetNbinsY()) < minEntries) {
	  continue;
	} else {
	  interval = 2;
	}
      } else {
	interval = 1;
      }
    }
    counter++;
    i += interval;
    //
    Double_t x = hist->GetXaxis()->GetBinCenter(i); 
    TH1D * projectionHist = hist->ProjectionY("dummy",i,i + interval);
    TF1 funcGaus("funcGaus","gaus");
    projectionHist->Fit(&funcGaus,"QN");
    //  
    xvec[counter-1] = x;
    yvec[counter-1] = funcGaus.GetParameter(1);
    delete projectionHist;
  }
  
  TGraph * graph = new TGraph(counter, xvec, yvec);

  return graph;
}



void AliTPCcalibTimeGain::BinLogX(THnSparse *h, Int_t axisDim) {

  // Method for the correct logarithmic binning of histograms

  TAxis *axis = h->GetAxis(axisDim);
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *new_bins = new Double_t[bins + 1];
   
  new_bins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete new_bins;
  
}


void AliTPCcalibTimeGain::BinLogX(TH1 *h) {

  // Method for the correct logarithmic binning of histograms

  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *new_bins = new Double_t[bins + 1];
   
  new_bins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete new_bins;
  
}



void AliTPCcalibTimeGain::CalculateBetheAlephParams(TH2F *hist, Double_t * ini) {
  //{0.0762*MIP,10.632,1.34e-05,1.863,1.948}
  const Double_t sigma = 0.06;

  TH2F *histBG = new TH2F("histBG","dEdxBg; #beta #gamma; TPC signal (a.u.)",hist->GetNbinsX(),0.1,5000.,hist->GetNbinsY(),0.5,5.);
  BinLogX(histBG);

  TF1 *foElectron = new TF1("foElectron", "(1.7/1.6)*AliExternalTrackParam::BetheBlochAleph(x/0.000511,[0],[1],[2],[3],[4])",0.01,100);
  TF1 *foMuon = new TF1("foMuon", "AliExternalTrackParam::BetheBlochAleph(x/0.1056,[0],[1],[2],[3],[4])",0.01,100);
  TF1 *foPion = new TF1("foPion", "AliExternalTrackParam::BetheBlochAleph(x/0.138,[0],[1],[2],[3],[4])",0.01,100);
  TF1 *foKaon = new TF1("foKaon", "AliExternalTrackParam::BetheBlochAleph(x/0.498,[0],[1],[2],[3],[4])",0.01,100);
  TF1 *foProton = new TF1("foProton", "AliExternalTrackParam::BetheBlochAleph(x/0.938,[0],[1],[2],[3],[4])",0.01,100);
  foElectron->SetParameters(ini);
  foMuon->SetParameters(ini);
  foPion->SetParameters(ini);
  foKaon->SetParameters(ini);
  foProton->SetParameters(ini);
  
  TCanvas *CanvCheck1 = new TCanvas();
  hist->Draw("colz");
  foElectron->Draw("same");
  foMuon->Draw("same");
  foPion->Draw("same");
  foKaon->Draw("same");  
  foProton->Draw("same");
 
  // Loop over all points of the input histogram
  
  for(Int_t i=1; i < hist->GetNbinsX(); i++) {
    Double_t x = hist->GetXaxis()->GetBinCenter(i);   
    for(Int_t j=1; j < hist->GetNbinsY(); j++) {
      Long64_t n = hist->GetBin(i, j);
      Double_t y = hist->GetYaxis()->GetBinCenter(j);
      Double_t entries = hist->GetBinContent(n);
      Double_t mass = 0;

      // 1. identify protons
      mass = 0.938;
      if (TMath::Abs(y - foProton->Eval(x))/foProton->Eval(x) < 4*sigma && x < 0.65 ) {
	for(Int_t iEntries =0; iEntries < entries; iEntries++) histBG->Fill(x/mass, y);
      }

      // 2. identify electrons
      mass = 0.000511;
      if (fIsCosmic) {
	if (TMath::Abs(y - foElectron->Eval(x))/foElectron->Eval(x) < 4*sigma && (x>0.25&&x<0.7) && fIsCosmic) {
	  for(Int_t iEntries =0; iEntries < entries; iEntries++) histBG->Fill(x/mass, y);
	}
      } else {
	if (TMath::Abs(y - foElectron->Eval(x))/foElectron->Eval(x) < 3*sigma && ((x>0.25&&x<0.35) || (x>1.5&&x<1.8) || (x>0.65&&x<0.7))) {
	  for(Int_t iEntries =0; iEntries < entries; iEntries++) histBG->Fill(x/mass, y);
	}
      }
      
      // 3. identify either muons or pions depending on cosmic or not
      if (fIsCosmic) {
	mass = 0.1056;
	if (TMath::Abs(y - foMuon->Eval(x))/foMuon->Eval(x) < 4*sigma && x > 0.25 && x < 100 ) {
	  for(Int_t iEntries =0; iEntries < entries; iEntries++) histBG->Fill(x/mass, y);
	}
      } else {
	mass = 0.1396;
	if (TMath::Abs(y - foPion->Eval(x))/foPion->Eval(x) < 4*sigma && x > 0.15 && x < 2) {
	  for(Int_t iEntries =0; iEntries < entries; iEntries++) histBG->Fill(x/mass, y);
	}
      }
      
      // 4. for pp also Kaons must be included
      if (!fIsCosmic) {
	mass = 0.4936;
	if (TMath::Abs(y - foKaon->Eval(x))/foKaon->Eval(x) < 4*sigma && x > 0.2 && x < 0.3 ) {
	  for(Int_t iEntries =0; iEntries < entries; iEntries++) histBG->Fill(x/mass, y);
	}
      }
    }
  }

  // Fit Aleph-Parameters to the obtained profile
  TF1 * funcAlephD = new TF1("AlephParametrizationD", "AliExternalTrackParam::BetheBlochAleph(x,[0],[1],[2],[3],[4])",0.3,10000);
  funcAlephD->SetParameters(ini);

  TCanvas *CanvCheck2 = new TCanvas();
  histBG->Draw();
  
  //FitSlices
  TObjArray * arr = new TObjArray();
  histBG->FitSlicesY(0,0,-1,0,"QN",arr);
  TH1D * fitMean = (TH1D*) arr->At(1);
  fitMean->Draw("same");

  funcAlephD->SetParLimits(2,1e-3,1e-7);
  funcAlephD->SetParLimits(3,0.5,3.5);
  funcAlephD->SetParLimits(4,0.5,3.5);
  fitMean->Fit(funcAlephD, "QNR");
  funcAlephD->Draw("same");

  for(Int_t i=0;i<5;i++) ini[i] = funcAlephD->GetParameter(i);

  foElectron->SetParameters(ini);
  foMuon->SetParameters(ini);
  foPion->SetParameters(ini);
  foKaon->SetParameters(ini);
  foProton->SetParameters(ini);
  
  TCanvas *CanvCheck3 = new TCanvas();
  hist->Draw("colz");
  foElectron->Draw("same");
  foMuon->Draw("same");
  foPion->Draw("same");
  foKaon->Draw("same");  
  foProton->Draw("same");
  
  CanvCheck1->Print();
  CanvCheck2->Print();
  CanvCheck3->Print();

  return;


}
