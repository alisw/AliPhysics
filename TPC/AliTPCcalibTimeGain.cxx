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

This class provides the calibration of the time dependence of the TPC gain due to pressure and temperature changes etc.


//0.  Libraries to load
gSystem->Load("libANALYSIS");
gSystem->Load("libSTAT");
gSystem->Load("libTPCcalib");


//1. Do calibration ...
//
// compare reference

//
//2. Visualize results
//
TFile fcalib("CalibObjects.root");
TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
AliTPCcalibTimeGain * gain = ( AliTPCcalibTimeGain *)array->FindObject("calibTimeGain");
TGraphErrors * gr = gain->GetGraphGainVsTime(0,1000)

// gain->GetHistGainTime()->GetAxis(1)->SetRangeUser(1213.8e6,1214.3e6)
TH2D * GainTime = gain->GetHistGainTime()->Projection(0,1)
GainTime->GetXaxis()->SetTimeDisplay(kTRUE)
GainTime->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}")
GainTime->Draw("colz")

//
// MakeSlineFit example
//
AliSplineFit *fit = AliTPCcalibTimeGain::MakeSplineFit(gr)

TGraph * grfit = fit.MakeGraph(gr->GetX()[0],gr->GetX()[gr->GetN()-1],50000,0);

gr->SetMarkerStyle(25);
gr->Draw("lp");
grfit->SetLineColor(2);
grfit->Draw("lu");

//
// QA - dE/dx resoultion as a function of time
//TCa

TGraph * grSigma = AliTPCcalibBase::FitSlicesSigma(gain->GetHistGainTime(),0,1,1800,5)

TCanvas *c1 = new TCanvas("c1","transparent pad",200,10,700,500);
   TPad *pad1 = new TPad("pad1","",0,0,1,1);
   TPad *pad2 = new TPad("pad2","",0,0,1,1);
   pad2->SetFillStyle(4000); //will be transparent
   pad1->Draw();
   pad1->cd();

GainTime->Draw("colz")
gr->Draw("lp")



  c1->cd();
 Double_t ymin = -0.04;
 Double_t ymax = 0.12;
Double_t dy = (ymax-ymin)/0.8;
Double_t xmin = GainTime->GetXaxis()->GetXmin()
Double_t xmax = GainTime->GetXaxis()->GetXmax()
Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
   pad2->Draw();
   pad2->cd();
grSigma->SetLineColor(2);
grSigma->SetLineWidth(2);
grSigma->Draw("lp")
TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
   axis->SetLabelColor(kRed);
   axis->SetTitle("dE/dx resolution #sigma_{dE/dx}");
   axis->Draw();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 ----> make Attachment study

TFile fcalib("CalibObjects40366b.root");
TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
AliTPCcalibTimeGain * gain = ( AliTPCcalibTimeGain *)array->FindObject("calibTimeGain");
TGraphErrors * grAttach = gain->GetGraphAttachment(2000,4)

TCanvas *c1 = new TCanvas("c1","transparent pad",200,10,700,500);
   TPad *pad1 = new TPad("pad1","",0,0,1,1);
   TPad *pad2 = new TPad("pad2","",0,0,1,1);
   pad2->SetFillStyle(4000); //will be transparent
   pad1->Draw();
   pad1->cd();

gain->GetHistGainTime()->GetAxis(1)->SetRangeUser(1213.8e6,1214.3e6)
TH2D * GainTime = gain->GetHistGainTime()->Projection(0,1)
GainTime->GetXaxis()->SetTimeDisplay(kTRUE)
GainTime->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}")
GainTime->Draw("colz")
//gr->Draw("lp")

  c1->cd();
 Double_t ymin = -0.001;
 Double_t ymax = 0.001;
Double_t dy = (ymax-ymin)/0.8;
Double_t xmin = GainTime->GetXaxis()->GetXmin()
Double_t xmax = GainTime->GetXaxis()->GetXmax()
Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
   pad2->Draw();
   pad2->cd();
grAttach->SetLineColor(2);
grAttach->SetLineWidth(2);
grAttach->Draw("lp")
TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
   axis->SetLabelColor(kRed);
   axis->SetTitle("attachment coefficient b");
   axis->Draw();


*/


#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TGraphErrors.h"
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
#include "AliESDv0.h"

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
   fUseMax(0),
   fLowerTrunc(0),
   fUpperTrunc(0),
   fUseShapeNorm(0),
   fUsePosNorm(0),
   fUsePadNorm(0),
   fUseCookAnalytical(0),
   fIsCosmic(0),
   fLowMemoryConsumption(0)
{  
  //
  // Default constructor
  //
  AliInfo("Default Constructor");  
}


AliTPCcalibTimeGain::AliTPCcalibTimeGain(const Text_t *name, const Text_t *title, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeGain)
  :AliTPCcalibBase(),
   fHistGainTime(0),
   fGainVsTime(0),
   fHistDeDxTotal(0),
   fIntegrationTimeDeDx(0),
   fMIP(0),
   fUseMax(0),
   fLowerTrunc(0),
   fUpperTrunc(0),
   fUseShapeNorm(0),
   fUsePosNorm(0),
   fUsePadNorm(0),
   fUseCookAnalytical(0),
   fIsCosmic(0),
   fLowMemoryConsumption(0)
{
  //
  // No default constructor
  //
  SetName(name);
  SetTitle(title);
  
  AliInfo("Non Default Constructor");
  
  fIntegrationTimeDeDx = deltaIntegrationTimeGain;
  
  Double_t deltaTime = EndTime - StartTime;
  
  
  // main histogram for time dependence: dE/dx, time, type (1-muon cosmic,2-pion beam data, 3&4 - proton points at higher dE/dx), meanDriftlength, momenta (only filled if enough space is available), run number, eta
  Int_t timeBins = TMath::Nint(deltaTime/deltaIntegrationTimeGain);
  Int_t binsGainTime[7]    = {150,  timeBins,    4,  25, 200, 10000000, 20};
  Double_t xminGainTime[7] = {0.5, StartTime,  0.5,   0, 0.1,    -0.5,  -1};
  Double_t xmaxGainTime[7] = {  8,   EndTime,  4.5, 250,  50, 9999999.5, 1};
  fHistGainTime = new THnSparseF("HistGainTime","dEdx time dep.;dEdx,time,type,driftlength,momenta,run number, eta;dEdx",7,binsGainTime,xminGainTime,xmaxGainTime);
  BinLogX(fHistGainTime, 4);
  //
  fHistDeDxTotal = new TH2F("DeDx","dEdx; momentum p (GeV); TPC signal (a.u.)",250,0.01,100.,1000,0.,8);
  BinLogX(fHistDeDxTotal);
  
  // default values for dE/dx
  fMIP = 50.;
  fUseMax = kTRUE;
  fLowerTrunc = 0.02;
  fUpperTrunc = 0.6;
  fUseShapeNorm = kTRUE;
  fUsePosNorm = kFALSE;
  fUsePadNorm = kFALSE;
  fUseCookAnalytical = kFALSE;
  //
  fIsCosmic = kTRUE;
  fLowMemoryConsumption = kTRUE;
  //
  
}



AliTPCcalibTimeGain::~AliTPCcalibTimeGain(){
  //
  // Destructor
  //
  delete fHistGainTime;
  delete fGainVsTime;
  delete fHistDeDxTotal;
}


void AliTPCcalibTimeGain::Process(AliESDEvent *event) {
  //
  // main track loop
  //
  if (!event) {
    Printf("ERROR: ESD not available");
    return;
  }
  
  if (fIsCosmic) { // this should be removed at some point based on trigger mask !?
    ProcessCosmicEvent(event);
  } else {
    ProcessBeamEvent(event);
  }
  

  
  
}


void AliTPCcalibTimeGain::ProcessCosmicEvent(AliESDEvent *event) {
  //
  // Process in case of cosmic event
  //
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!esdFriend) {
   Printf("ERROR: ESDfriend not available");
   return;
  }
  //
  UInt_t time = event->GetTimeStamp();
  Int_t ntracks = event->GetNumberOfTracks();
  Int_t runNumber = event->GetRunNumber();
  //
  // track loop
  //
  for (Int_t i=0;i<ntracks;++i) {

    AliESDtrack *track = event->GetTrack(i);
    if (!track) continue;
    AliESDfriendTrack *friendTrack = esdFriend->GetTrack(i);
    if (!friendTrack) continue;        
    const AliExternalTrackParam * trackIn = track->GetInnerParam();
    const AliExternalTrackParam * trackOut = friendTrack->GetTPCOut();
    if (!trackIn) continue;
    if (!trackOut) continue;

    // calculate necessary track parameters
    Double_t meanP = trackIn->GetP();
    Double_t meanDrift = 250 - 0.5*TMath::Abs(trackIn->GetZ() + trackOut->GetZ());
    Int_t nclsDeDx = track->GetTPCNcls();

    // exclude tracks which do not look like primaries or are simply too short or on wrong sectors
    if (nclsDeDx < 60) continue;     
    if (TMath::Abs(trackIn->GetTgl()) > 1) continue;
    if (TMath::Abs(trackIn->GetSnp()) > 0.6) continue;
    
    // Get seeds
    TObject *calibObject;
    AliTPCseed *seed = 0;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }    

    if (seed) { 
      Double_t tpcSignal = GetTPCdEdx(seed);
      if (nclsDeDx > 100) fHistDeDxTotal->Fill(meanP, tpcSignal);
      //
      if (fLowMemoryConsumption) {
	if (meanP < 20) continue;
	meanP = 30; // set all momenta to one in order to save memory
      }
      //dE/dx, time, type (1-muon cosmic,2-pion beam data), momenta
      Double_t vec[6] = {tpcSignal,time,1,meanDrift,meanP,runNumber};
      fHistGainTime->Fill(vec);
    }
    
  }

}



void AliTPCcalibTimeGain::ProcessBeamEvent(AliESDEvent *event) {
  //
  // Process in case of beam event
  //
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!esdFriend) {
   Printf("ERROR: ESDfriend not available");
   return;
  }
  //
  UInt_t time = event->GetTimeStamp();
  Int_t ntracks = event->GetNumberOfTracks();
  Int_t runNumber = event->GetRunNumber();
  //
  // track loop
  //
  for (Int_t i=0;i<ntracks;++i) { // begin track loop

    AliESDtrack *track = event->GetTrack(i);
    if (!track) continue;
    AliESDfriendTrack *friendTrack = esdFriend->GetTrack(i);
    if (!friendTrack) continue;
        
    const AliExternalTrackParam * trackIn = track->GetInnerParam();
    const AliExternalTrackParam * trackOut = friendTrack->GetTPCOut();
    if (!trackIn) continue;
    if (!trackOut) continue;

    // calculate necessary track parameters
    Double_t meanP = trackIn->GetP();
    Double_t meanDrift = 250 - 0.5*TMath::Abs(trackIn->GetZ() + trackOut->GetZ());
    Int_t nclsDeDx = track->GetTPCNcls();

    // exclude tracks which do not look like primaries or are simply too short or on wrong sectors
    if (nclsDeDx < 60) continue;     
    //if (TMath::Abs(trackIn->GetTgl()) > 1) continue;
    //if (TMath::Abs(trackIn->GetSnp()) > 0.6) continue;
    if (TMath::Abs(trackIn->Eta()) > 1) continue;
    UInt_t status = track->GetStatus();
    if ((status&AliESDtrack::kTPCrefit)==0) continue;
    //if (track->GetNcls(0) < 3) continue; // ITS clusters
    Float_t dca[2], cov[3];
    track->GetImpactParameters(dca,cov);
    if (TMath::Abs(dca[0]) > 7 || TMath::Abs(dca[0]) < 0.0000001 || TMath::Abs(dca[1]) > 25 ) continue; // cut in xy
    Double_t eta = trackIn->Eta();
    
    // Get seeds
    TObject *calibObject;
    AliTPCseed *seed = 0;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }    

    if (seed) {
      Int_t particleCase = 0;
      if (meanP < 0.5  && meanP > 0.4)  particleCase = 2; // MIP pions
      if (meanP < 0.57 && meanP > 0.56) particleCase = 3; // protons 1
      if (meanP < 0.66 && meanP > 0.65) particleCase = 4; // protons 2
      //
      if (fLowMemoryConsumption && particleCase == 0) continue;
      //
      Double_t tpcSignal = GetTPCdEdx(seed);
      fHistDeDxTotal->Fill(meanP, tpcSignal);
      //
      //dE/dx, time, type (1-muon cosmic,2-pion beam data, 3&4 protons), momenta, runNumner, eta
      Double_t vec[7] = {tpcSignal,time,particleCase,meanDrift,meanP,runNumber, eta};
      fHistGainTime->Fill(vec);
    }
    
  } // end track loop
  //
  // V0 loop -- in beam events the cosmic part of the histogram is filled with GammaConversions
  //
  for(Int_t iv0 = 0; iv0 < event->GetNumberOfV0s(); iv0++) {
    AliESDv0 * v0 = event->GetV0(iv0);
    if (!v0->GetOnFlyStatus()) continue;
    if (v0->GetEffMass(0,0) > 0.02) continue; // select low inv. mass
    Double_t xyz[3];
    v0->GetXYZ(xyz[0], xyz[1], xyz[2]);
    if (TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) < 3 || TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) > 30) continue;
    //
    // "loop over daughters" 
    //
    for(Int_t idaughter = 0; idaughter < 2; idaughter++) { // daughter loop
      Int_t index = idaughter == 0 ? v0->GetPindex() : v0->GetNindex();
      AliESDtrack * trackP = event->GetTrack(index);
      AliESDfriendTrack *friendTrackP = esdFriend->GetTrack(index);
      if (!friendTrackP) continue;
      const AliExternalTrackParam * trackPIn = trackP->GetInnerParam();
      const AliExternalTrackParam * trackPOut = friendTrackP->GetTPCOut();
      if (!trackPIn) continue;
      if (!trackPOut) continue;
      // calculate necessary track parameters
      Double_t meanP = trackPIn->GetP();
      Double_t meanDrift = 250 - 0.5*TMath::Abs(trackPIn->GetZ() + trackPOut->GetZ());
      Int_t nclsDeDx = trackP->GetTPCNcls();
      // exclude tracks which do not look like primaries or are simply too short or on wrong sectors
      if (nclsDeDx < 60) continue;     
      if (TMath::Abs(trackPIn->GetTgl()) > 1) continue;
      //
      TObject *calibObject;
      AliTPCseed *seed = 0;
      for (Int_t l=0;(calibObject=friendTrackP->GetCalibObject(l));++l) {
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
      }    
      if (seed) { 
	if (fLowMemoryConsumption) {
	  if (meanP > 0.5 || meanP < 0.4) continue;
	  meanP = 0.45; // set all momenta to one in order to save memory
      }
	Double_t tpcSignal = GetTPCdEdx(seed);
	//dE/dx, time, type (1-muon cosmic,2-pion beam data), momenta
	Double_t vec[6] = {tpcSignal,time,1,meanDrift,meanP,runNumber};
	fHistGainTime->Fill(vec);
      }
    }
    
  }

}


Float_t AliTPCcalibTimeGain::GetTPCdEdx(AliTPCseed * seed) {
  //
  // calculate tpc dEdx
  //
  Double_t signal = 0;
  //
  if (!fUseCookAnalytical) {
    signal = (1/fMIP)*seed->CookdEdxNorm(fLowerTrunc,fUpperTrunc,fUseMax,0,159,fUseShapeNorm,fUsePosNorm,fUsePadNorm,0);
  } else {
    signal = (1/fMIP)*seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,fUseMax);
  }
  //
  return signal;
}


void AliTPCcalibTimeGain::AnalyzeRun(Int_t minEntries) {
  //
  // Analyze results of calibration
  //
  if (fIsCosmic) {
    fHistGainTime->GetAxis(2)->SetRangeUser(0.51,1.49); // only cosmics
    fHistGainTime->GetAxis(4)->SetRangeUser(20,100);    // only Fermi-Plateau muons
  } else {
    fHistGainTime->GetAxis(2)->SetRangeUser(1.51,2.49); // only beam data
    fHistGainTime->GetAxis(4)->SetRangeUser(0.39,0.51); // only MIP pions
  }
  //
  fGainVsTime = AliTPCcalibBase::FitSlices(fHistGainTime,0,1,minEntries,10);
  //
  return;
}


TGraphErrors * AliTPCcalibTimeGain::GetGraphGainVsTime(Int_t runNumber, Int_t minEntries) {
  //
  // Analyze results and get the graph 
  //
  if (runNumber == 0) {
    if (!fGainVsTime) {
      AnalyzeRun(minEntries);
    }
  } else {
    // 1st check if the current run was cosmic or beam event
    fHistGainTime->GetAxis(5)->SetRangeUser(runNumber,runNumber);
    AnalyzeRun(minEntries);
  }
  if (fGainVsTime->GetN() == 0) return 0;
  return fGainVsTime;
}

Long64_t AliTPCcalibTimeGain::Merge(TCollection *li) {
  //
  // merge component
  //
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


AliSplineFit * AliTPCcalibTimeGain::MakeSplineFit(TGraphErrors * graph) {
  //
  // make spline fit of gain
  //
  AliSplineFit *fit = new AliSplineFit();
  fit->SetGraph(graph);
  fit->SetMinPoints(graph->GetN()+1);
  fit->InitKnots(graph,2,0,0.001);
  fit->SplineFit(0);
  return fit;
  
}



TGraphErrors * AliTPCcalibTimeGain::GetGraphAttachment(Int_t minEntries, Int_t nmaxBin, Float_t /*fracLow*/, Float_t /*fracUp*/) {
  //
  // For each time bin the driftlength-dependence of the signal is fitted.
  //
  TH3D * hist = fHistGainTime->Projection(1, 0, 3);
  Double_t *xvec = new Double_t[hist->GetNbinsX()];
  Double_t *yvec = new Double_t[hist->GetNbinsX()];
  Double_t *xerr = new Double_t[hist->GetNbinsX()];
  Double_t *yerr = new Double_t[hist->GetNbinsX()];
  Int_t counter  = 0;
  TH2D * projectionHist = 0x0;
  //
  for(Int_t i=1; i < hist->GetNbinsX(); i++) {
    Int_t nsum=0;
    Int_t imin   =  i;
    Int_t imax   =  i;    
    for (Int_t idelta=0; idelta<nmaxBin; idelta++){
      //
      imin   =  TMath::Max(i-idelta,1);
      imax   =  TMath::Min(i+idelta,hist->GetNbinsX());
      nsum = TMath::Nint(hist->Integral(imin,imax,1,hist->GetNbinsY()-1,1,hist->GetNbinsZ()-1));
      if (nsum==0) break;
      if (nsum>minEntries) break;
    }
    if (nsum<minEntries) continue;
    //
    hist->GetXaxis()->SetRange(imin,imax);
    projectionHist = (TH2D*)hist->Project3D("yzNUFNOF");
    //
    TObjArray arr;
    projectionHist->FitSlicesY(0,2, projectionHist->GetNbinsX()-2,0,"QNR",&arr);
    TH1D * histAttach = (TH1D*)arr.At(1);
    TF1 pol1("polynom1","pol1",10,240);
    histAttach->Fit(&pol1,"QNR");
    xvec[counter] = 0.5*(hist->GetXaxis()->GetBinCenter(imin)+hist->GetXaxis()->GetBinCenter(imax));
    yvec[counter] = pol1.GetParameter(1)/pol1.GetParameter(0);
    xerr[counter] = 0;
    yerr[counter] = pol1.GetParError(1)/pol1.GetParameter(0);
    counter++;
    //
    delete projectionHist;
  }
  
  TGraphErrors * graphErrors = new TGraphErrors(counter, xvec, yvec, xerr, yerr);
  delete [] xvec;
  delete [] yvec;
  delete [] xerr;
  delete [] yerr;
  delete hist;
  return graphErrors;
  
}



void AliTPCcalibTimeGain::BinLogX(THnSparse *h, Int_t axisDim) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = h->GetAxis(axisDim);
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
   
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete[] newBins;
  
}


void AliTPCcalibTimeGain::BinLogX(TH1 *h) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
   
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete[] newBins;
  
}



void AliTPCcalibTimeGain::CalculateBetheAlephParams(TH2F *hist, Double_t * ini) {
  //
  // Fit the bethe bloch params
  //
  //{0.0762*MIP,10.632,1.34e-05,1.863,1.948}
  const Double_t sigma = 0.06;

  TH2F *histBG = new TH2F("histBG","dEdxBg; #beta #gamma; TPC signal (a.u.)",TMath::Nint(0.5*hist->GetNbinsX()),0.1,5000.,hist->GetNbinsY(),hist->GetYaxis()->GetXmin(),hist->GetYaxis()->GetXmax());
  BinLogX(histBG);

  TF1 *foElectron = new TF1("foElectron", "AliExternalTrackParam::BetheBlochAleph(x/0.000511,[0],[1],[2],[3],[4])",0.01,100);
  TF1 *foMuon = new TF1("foMuon", "AliExternalTrackParam::BetheBlochAleph(x/0.1056,[0],[1],[2],[3],[4])",0.01,100);
  TF1 *foPion = new TF1("foPion", "AliExternalTrackParam::BetheBlochAleph(x/0.138,[0],[1],[2],[3],[4])",0.01,100);
  TF1 *foKaon = new TF1("foKaon", "AliExternalTrackParam::BetheBlochAleph(x/0.498,[0],[1],[2],[3],[4])",0.01,100);
  TF1 *foProton = new TF1("foProton", "AliExternalTrackParam::BetheBlochAleph(x/0.938,[0],[1],[2],[3],[4])",0.01,100);
  foElectron->SetParameters(ini);
  foMuon->SetParameters(ini);
  foPion->SetParameters(ini);
  foKaon->SetParameters(ini);
  foProton->SetParameters(ini);
  
  TCanvas *canvCheck1 = new TCanvas();
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
	if (TMath::Abs(y - foElectron->Eval(x))/foElectron->Eval(x) < 2*sigma && ((x>0.25&&x<0.35) || (x>1.5&&x<1.8) || (x>0.65&&x<0.7))) {
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

  TCanvas *canvCheck2 = new TCanvas();
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
  
  TCanvas *canvCheck3 = new TCanvas();
  hist->Draw("colz");
  foElectron->Draw("same");
  foMuon->Draw("same");
  foPion->Draw("same");
  foKaon->Draw("same");  
  foProton->Draw("same");
  
  canvCheck1->Print();
  canvCheck2->Print();
  canvCheck3->Print();

  return;


}
