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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Base class for the calibration components using 
//  as input TPCseeds and ESDs
//  Event loop outside of the component
//
//
// Base functionality to be implemeneted by component 
/* 
   //In some cases only one of this function to be implemented
   virtual void     Process(AliESDEvent *event)
   virtual void     Process(AliTPCseed *track)
   //
   virtual Long64_t Merge(TCollection *li);
   virtual void     Analyze()
   void             Terminate();
*/
// Functionality provided by base class for Algorith debuging:
//  TTreeSRedirector * cstream =  GetDebugStreamer() - get debug streamer which can be use for numerical debugging
//                      



//  marian.ivanov@cern.ch
// 
#include "AliTPCcalibBase.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTreeStream.h"
#include "TTimeStamp.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "AliRecoParam.h"


#include "AliLog.h"
#include "AliVEvent.h"


ClassImp(AliTPCcalibBase)

AliTPCcalibBase::AliTPCcalibBase():
    TNamed(),
    fDebugStreamer(0),
    fStreamLevel(0),   
    fRun(0),                  //!  current Run number
    fEvent(0),                //!  current Event number
    fTime(0),                 //!  current Time
    fTrigger(0),              //! current trigger type
    fMagF(0),                 //! current magnetic field
    fTriggerMaskReject(-1),   //trigger mask - reject trigger
    fTriggerMaskAccept(-1),   //trigger mask - accept trigger
    fHasLaser(kFALSE),                    //flag the laser is overlayed with given event 
    fRejectLaser(kTRUE),                 //flag- reject laser
    fTriggerClass(),
    fCurrentEvent(0),           //! current event
    fCurrentTrack(0),           //! current esd track
    fCurrentFriendTrack(0),           //! current esd track
    fCurrentSeed(0),            //! current seed
    fDebugLevel(0)
{
  //
  // Constructor
  //
}

AliTPCcalibBase::AliTPCcalibBase(const char * name, const char * title):
  TNamed(name,title),
  fDebugStreamer(0),
  fStreamLevel(0),   
  fRun(0),                  //!  current Run number
  fEvent(0),                //!  current Event number
  fTime(0),                 //!  current Time
  fTrigger(0),              //! current trigger type
  fMagF(0),                 //! current magnetic field
  fTriggerMaskReject(-1),   //trigger mask - reject trigger
  fTriggerMaskAccept(-1),   //trigger mask - accept trigger
  fHasLaser(kFALSE),                    //flag the laser is overlayed with given event 
  fRejectLaser(kTRUE),                 //flag- reject laser
  fTriggerClass(),
  fCurrentEvent(0),           //! current event
  fCurrentTrack(0),           //! current esd track
  fCurrentFriendTrack(0),           //! current esd track
  fCurrentSeed(0),            //! current seed
  fDebugLevel(0)
{
  //
  // Constructor
  //
}

AliTPCcalibBase::AliTPCcalibBase(const AliTPCcalibBase&calib):
  TNamed(calib),
  fDebugStreamer(0),
  fStreamLevel(calib.fStreamLevel),
  fRun(0),                  //!  current Run number
  fEvent(0),                //!  current Event number
  fTime(0),                 //!  current Time
  fTrigger(0),              //! current trigger type
  fMagF(0),                 //! current magnetic field
  fTriggerMaskReject(calib.fTriggerMaskReject),   //trigger mask - reject trigger
  fTriggerMaskAccept(calib.fTriggerMaskAccept),   //trigger mask - accept trigger
  fHasLaser(calib.fHasLaser),                    //flag the laser is overlayed with given event
  fRejectLaser(calib.fRejectLaser),                 //flag- reject laser
  fTriggerClass(calib.fTriggerClass),
  fCurrentEvent(0),           //! current event
  fCurrentTrack(0),           //! current esd track
  fCurrentFriendTrack(0),           //! current esd track
  fCurrentSeed(0),            //! current seed
  fDebugLevel(calib.fDebugLevel)
{
  //
  // copy constructor
  //
}

AliTPCcalibBase &AliTPCcalibBase::operator=(const AliTPCcalibBase&calib){
  //
  // operator=
  //
  if (this == &calib) return (*this);

  ((TNamed *)this)->operator=(calib);
  fDebugStreamer=0;
  fStreamLevel=calib.fStreamLevel;
  fDebugLevel=calib.fDebugLevel;
  return *this;
}


AliTPCcalibBase::~AliTPCcalibBase() {
  //
  // destructor
  //
  if (fDebugLevel>0) printf("AliTPCcalibBase::~AliTPCcalibBase\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer=0;
}

void  AliTPCcalibBase::Terminate(){
  //
  //
  //
  if (fDebugLevel>0) printf("AliTPCcalibBase::Terminate\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer = 0;
  return;
}

TTreeSRedirector *AliTPCcalibBase::GetDebugStreamer(){
  //
  // Get Debug streamer
  // In case debug streamer not yet initialized and StreamLevel>0 create new one
  //
  if (fStreamLevel==0) return 0;
  if (fDebugStreamer) return fDebugStreamer;
  TString dsName;
  dsName=GetName();
  dsName+="Debug.root";
  dsName.ReplaceAll(" ",""); 
  fDebugStreamer = new TTreeSRedirector(dsName.Data());
  return fDebugStreamer;
}


void    AliTPCcalibBase::UpdateEventInfo(AliVEvent * event){
  //
  //
  //
  fRun     = event->GetRunNumber();
  fEvent   = event->GetEventNumberInFile();
  fTime    = event->GetTimeStamp();
  fTrigger = event->GetTriggerMask();
  fMagF    = event->GetMagneticField();
  fTriggerClass = event->GetFiredTriggerClasses().Data();
  fHasLaser = HasLaser(event); 
  
}


Bool_t AliTPCcalibBase::HasLaser(AliVEvent *event){
  //
  //
  //
  // Use event specie
  Bool_t isLaser=kFALSE;
  UInt_t specie = event->GetEventSpecie();  // select only cosmic events
  if (specie==AliRecoParam::kCalib) {
    isLaser = kTRUE;
  }
  return isLaser;
}



Bool_t AliTPCcalibBase::AcceptTrigger(){
  //
  // Apply trigger mask - Don't do calibration for non proper triggers
  // 
  if (fTriggerMaskReject==(Int_t)fTrigger) return kFALSE;
  if (fTriggerMaskAccept>0 && fTriggerMaskAccept!=(Int_t)fTrigger) return kFALSE;
  if (fHasLaser && fRejectLaser) return kFALSE;
  return kTRUE;
}


void AliTPCcalibBase::RegisterDebugOutput(const char *path){
  //
  // store  - copy debug output to the destination position
  // currently ONLY for local copy
  if (fDebugLevel>0) printf("AliTPCcalibBase::RegisterDebugOutput(%s)\n",path);
  if (fStreamLevel==0) return;
  TString dsName;
  dsName=GetName();
  dsName+="Debug.root";
  dsName.ReplaceAll(" ",""); 
  TString dsName2=path;
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+=gSystem->HostName();
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+="/";
  TTimeStamp s;
  dsName2+=Int_t(s.GetNanoSec());
  dsName2+="/";
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+=dsName;
  AliInfo(Form("copy %s\t%s\n",dsName.Data(),dsName2.Data()));
  printf("copy %s\t%s\n",dsName.Data(),dsName2.Data());
  TFile::Cp(dsName.Data(),dsName2.Data());
}

TGraphErrors * AliTPCcalibBase::FitSlices(THnSparse *h, Int_t axisDim1, Int_t axisDim2, Int_t minEntries, Int_t nmaxBin, Float_t fracLow, Float_t fracUp, Bool_t useMedian, TTreeSRedirector *cstream, Int_t ival){
  //
  // Fitting slices of the projection(axisDim1,axisDim2) of a sparse histogram
  //
  TH2D * hist = h->Projection(axisDim1, axisDim2);
  TGraphErrors *gr=FitSlices(hist, minEntries, nmaxBin, fracLow, fracUp, useMedian, cstream, ival);
  delete hist;
  return gr;
}

TGraphErrors * AliTPCcalibBase::FitSlices(TH2* hist, Int_t minEntries, Int_t nmaxBin, Float_t fracLow, Float_t fracUp, Bool_t useMedian, TTreeSRedirector *cstream, Int_t ival){
  //
  // Fitting slices of the projection(axisDim1,axisDim2) of a sparse histogram
  // 
  const Int_t kMinStat=100;
  TF1 funcGaus("funcGaus","gaus");
  Double_t *xvec = new Double_t[hist->GetNbinsX()];
  Double_t *yvec = new Double_t[hist->GetNbinsX()];
  Double_t *xerr = new Double_t[hist->GetNbinsX()];
  Double_t *yerr = new Double_t[hist->GetNbinsX()];
  Int_t counter  = 0;
  TH1D * projectionHist =0;
  //

  for(Int_t i=1; i <= hist->GetNbinsX(); i++) {
    Int_t nsum=0;
    Int_t imin   =  i;
    Int_t imax   =  i;

    for (Int_t idelta=0; idelta<nmaxBin; idelta++){
      //
      imin   =  TMath::Max(i-idelta,1);
      imax   =  TMath::Min(i+idelta,hist->GetNbinsX());
      nsum = TMath::Nint(hist->Integral(imin,imax,1,hist->GetNbinsY()));
      if (nsum==0) break;
      if (nsum>minEntries) break;
    }
    if (nsum<kMinStat) continue;
    //
    hist->GetXaxis()->SetRange(imin,imax);
    projectionHist = hist->ProjectionY("gain",imin,imax);
    // Determine Median:
    Float_t xMin = projectionHist->GetXaxis()->GetXmin();
    Float_t xMax = projectionHist->GetXaxis()->GetXmax();
    Float_t xMedian = (xMin+xMax)*0.5;
    Float_t integral = 0;
    for(Int_t jbin=1; jbin<projectionHist->GetNbinsX()-1; jbin++) {
      integral+=projectionHist->GetBinContent(jbin);
    }
    //printf("Integral %f\t%f\n",integral, projectionHist->GetSum());
    //
    //
    Float_t currentSum=0;
    for(Int_t jbin=1; jbin<projectionHist->GetNbinsX()-1; jbin++) {
      currentSum += projectionHist->GetBinContent(jbin);
      if (currentSum<fracLow*integral) xMin = projectionHist->GetBinCenter(jbin);
      if (currentSum<fracUp*integral)  xMax = projectionHist->GetBinCenter(jbin+1);      
      if (currentSum<0.5*integral && projectionHist->GetBinContent(jbin)>0){
	xMedian = (projectionHist->GetBinCenter(jbin)*projectionHist->GetBinContent(jbin)+
		   projectionHist->GetBinCenter(jbin+1)*projectionHist->GetBinContent(jbin+1))/
	  (projectionHist->GetBinContent(jbin)+projectionHist->GetBinContent(jbin+1));
      }
    }
    //
    Float_t rms  = projectionHist->GetRMS();
    //    i += interval;
    //
    Double_t xcenter =  hist->GetMean(); 
    Double_t xrms    =  hist->GetRMS()+hist->GetXaxis()->GetBinWidth(1)/TMath::Sqrt(12.); 
    Double_t binWidth = projectionHist->GetXaxis()->GetBinWidth(1);
    if (rms>0){
      // cut on +- 4 RMS
      projectionHist->Fit(&funcGaus,"QN","",xMin, xMax);
      Double_t chi2 = funcGaus.GetChisquare();
      //  
      xvec[counter] = xcenter;
      yvec[counter] = funcGaus.GetParameter(ival);
      xerr[counter] = xrms;
      yerr[counter] = funcGaus.GetParError(ival); 
      if (useMedian) yvec[counter] = xMedian;
      if (cstream){
	(*cstream)<<"fitDebug"<<
	  "xcenter="<<xcenter<<
	  "xMin="<<xMin<<
	  "xMax="<<xMax<<
	  "xMedian="<<xMedian<<
	  "xFitM"<<yvec[counter]<<
	  "xFitS"<<yerr[counter]<<
	  "chi2="<<chi2<<	  
	  "\n";
      }
      counter++;
    }else{
      xvec[counter] = xcenter;
      yvec[counter] = xMedian;
      xerr[counter] = xrms;
      yerr[counter] = binWidth/TMath::Sqrt(12.); 
      counter++;
    }
    delete projectionHist;
  }
  
  TGraphErrors * graphErrors = new TGraphErrors(counter, xvec, yvec, xerr, yerr);
  delete [] xvec;
  delete [] yvec;
  delete [] xerr;
  delete [] yerr;
  return graphErrors;
}

TH2* AliTPCcalibBase::NormalizedProjection(THnSparse *h, Int_t axisDim1, Int_t axisDim2, Int_t normDim, Float_t minStatFrac)
{
  // Create projection in axisDim1, axisDim2 in slices of normDim
  // Normalize the the entries in the dimension normDim and the sum up
  // Don't use bin for sum if statistics is less than minStatFrac of average number of entries

  if (!h) return 0x0;
  // ---| axis ranges of dimension in which to normalize |----------------------
  TAxis *axisNorm=h->GetAxis(normDim);
  const Int_t first=axisNorm->GetFirst();
  const Int_t last =axisNorm->GetLast();
  const Int_t nbins=last-first+1;

  // ---| final merged histogram |----------------------------------------------
  TH2 * hist = h->Projection(axisDim1, axisDim2);
  const Double_t statPerBin = hist->Integral()/Double_t(nbins);
  hist->Reset();
  hist->SetName(Form("%s_norm", hist->GetName()));
//   printf("statPerBin: %.2f\n", statPerBin);
//   printf("first, last, bins: %i, %i, %i\n", first, last, nbins);

  // ---| array with 2D projections per bin and statistics
  TObjArray arrProjections(nbins);
  arrProjections.SetOwner();
  Double_t *values=new Double_t[nbins];

  // ---| loop over normDim bins and do projections |---------------------------
  Int_t nused=0;
  for (Int_t ibin=first; ibin<=last; ++ibin) {
    axisNorm->SetRange(ibin,ibin);
    TH2 *proj=h->Projection(axisDim1, axisDim2);
    if (!proj) continue;
    proj->SetName(Form("%s_%d", h->GetName(), ibin));
//     proj->SetDirectory(0x0);

    proj->SetUniqueID(0); //steer nomalisation

    const Double_t stat=proj->Integral();
//     printf("bin: %i (%.2f)\n", ibin, stat);
    if (stat>minStatFrac*statPerBin) {
      proj->SetUniqueID(1);
      proj->Scale(1./stat);
      values[nused]=stat;
      ++nused;
//       printf("  used\n");
    }
    arrProjections.Add(proj);
  }

  // ---| Get Median of used bins, scale and sum up |---------------------------
  Double_t median=TMath::Median(nused, values);
//   printf("median: %.2f\n", median);
  for (Int_t ihist=0; ihist<nbins; ++ihist) {
    TH2 *proj=static_cast<TH2*>(arrProjections.At(ihist));
    if (!proj) continue;
    if (proj->GetUniqueID()) proj->Scale(median);
    hist->Add(proj);
  }

  axisNorm->SetRange(first,last);
  delete [] values;

  return hist;
}


void AliTPCcalibBase::BinLogX(THnSparse *h, Int_t axisDim) {

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
  delete [] new_bins;

}
void AliTPCcalibBase::BinLogX(TH1 *h) {

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
  delete [] new_bins;

}

void AliTPCcalibBase::BinLogX(TAxis *axis) {

  // Method for the correct logarithmic binning of histograms

  Int_t bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  if (from<0) return;
  Double_t *new_bins = new Double_t[bins + 1];

  new_bins[0] = from;
  Double_t factor = pow(to/from, 1./bins);

  for (int i = 1; i <= bins; i++) {
   new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete [] new_bins;
}
