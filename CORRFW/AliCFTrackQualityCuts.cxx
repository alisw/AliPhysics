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

// The class AliCFTrackQualityCuts is designed to select reconstructed tracks
// of high quality and to provide corresponding QA histograms.
// This class inherits from the Analysis' Framework abstract base class
// AliAnalysisCuts and is a part of the Correction Framework.
// This class acts on single, reconstructed tracks, it is applicable on
// ESD and AOD data.
// It mainly consists of a IsSelected function that returns a boolean.
// This function checks whether the considered track passes a set of cuts:
// - number of clusters in the TPC
// - number of clusters in the ITS
// - chi2 / cluster in the TPC
// - chi2 / cluster in the ITS
// - successful TPC refit
// - successful ITS refit
// - covariance matrix diagonal elements
//
// The cut values for these cuts are set with the corresponding set functions.
// All cut classes provided by the correction framework are supposed to be
// added in the Analysis Framwork's class AliAnalysisFilter and applied by
// the filter via a loop.
//
// author: I. Kraus (Ingrid.Kraus@cern.ch)
// idea taken form
// AliESDtrackCuts writte by Jan Fiete Grosse-Oetringhaus and
// AliRsnDaughterCut class written by A. Pulvirenti.


#include <TCanvas.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TBits.h>

#include <AliESDtrack.h>
#include <AliLog.h>
#include "AliCFTrackQualityCuts.h"
#include "AliAODTrack.h"

ClassImp(AliCFTrackQualityCuts)

//__________________________________________________________________________________
AliCFTrackQualityCuts::AliCFTrackQualityCuts() :
  AliCFCutBase(),
  fMinNClusterTPC(0),
  fMinNClusterITS(0),
  fMaxChi2PerClusterTPC(0),
  fMaxChi2PerClusterITS(0),
  fCovariance11Max(0),
  fCovariance22Max(0),
  fCovariance33Max(0),
  fCovariance44Max(0),
  fCovariance55Max(0),
  fStatus(0),
  fhCutStatistics(0),
  fhCutCorrelation(0),
  fBitmap(0x0),
  fhNBinsClusterTPC(0),
  fhNBinsClusterITS(0),
  fhNBinsChi2TPC(0),
  fhNBinsChi2ITS(0),
  fhNBinsCovariance11(0),
  fhNBinsCovariance22(0),
  fhNBinsCovariance33(0),
  fhNBinsCovariance44(0),
  fhNBinsCovariance55(0),
  fhBinLimClusterTPC(0x0),
  fhBinLimClusterITS(0x0),
  fhBinLimChi2TPC(0x0),
  fhBinLimChi2ITS(0x0),
  fhBinLimCovariance11(0x0),
  fhBinLimCovariance22(0x0),
  fhBinLimCovariance33(0x0),
  fhBinLimCovariance44(0x0),
  fhBinLimCovariance55(0x0)
{
  //
  // Default constructor
  //
  Initialise();
}
//__________________________________________________________________________________
AliCFTrackQualityCuts::AliCFTrackQualityCuts(Char_t* name, Char_t* title) :
  AliCFCutBase(name,title),
  fMinNClusterTPC(0),
  fMinNClusterITS(0),
  fMaxChi2PerClusterTPC(0),
  fMaxChi2PerClusterITS(0),
  fCovariance11Max(0),
  fCovariance22Max(0),
  fCovariance33Max(0),
  fCovariance44Max(0),
  fCovariance55Max(0),
  fStatus(0),
  fhCutStatistics(0),
  fhCutCorrelation(0),
  fBitmap(0x0),
  fhNBinsClusterTPC(0),
  fhNBinsClusterITS(0),
  fhNBinsChi2TPC(0),
  fhNBinsChi2ITS(0),
  fhNBinsCovariance11(0),
  fhNBinsCovariance22(0),
  fhNBinsCovariance33(0),
  fhNBinsCovariance44(0),
  fhNBinsCovariance55(0),
  fhBinLimClusterTPC(0x0),
  fhBinLimClusterITS(0x0),
  fhBinLimChi2TPC(0x0),
  fhBinLimChi2ITS(0x0),
  fhBinLimCovariance11(0x0),
  fhBinLimCovariance22(0x0),
  fhBinLimCovariance33(0x0),
  fhBinLimCovariance44(0x0),
  fhBinLimCovariance55(0x0)
{
  //
  // Constructor
  //
  Initialise();
}
//__________________________________________________________________________________
AliCFTrackQualityCuts::AliCFTrackQualityCuts(const AliCFTrackQualityCuts& c) :
  AliCFCutBase(c),
  fMinNClusterTPC(c.fMinNClusterTPC),
  fMinNClusterITS(c.fMinNClusterITS),
  fMaxChi2PerClusterTPC(c.fMaxChi2PerClusterTPC),
  fMaxChi2PerClusterITS(c.fMaxChi2PerClusterITS),
  fCovariance11Max(c.fCovariance11Max),
  fCovariance22Max(c.fCovariance22Max),
  fCovariance33Max(c.fCovariance33Max),
  fCovariance44Max(c.fCovariance44Max),
  fCovariance55Max(c.fCovariance55Max),
  fStatus(c.fStatus),
  fhCutStatistics(c.fhCutStatistics),
  fhCutCorrelation(c.fhCutCorrelation),
  fBitmap(c.fBitmap),
  fhNBinsClusterTPC(c.fhNBinsClusterTPC),
  fhNBinsClusterITS(c.fhNBinsClusterITS),
  fhNBinsChi2TPC(c.fhNBinsChi2TPC),
  fhNBinsChi2ITS(c.fhNBinsChi2ITS),
  fhNBinsCovariance11(c.fhNBinsCovariance11),
  fhNBinsCovariance22(c.fhNBinsCovariance22),
  fhNBinsCovariance33(c.fhNBinsCovariance33),
  fhNBinsCovariance44(c.fhNBinsCovariance44),
  fhNBinsCovariance55(c.fhNBinsCovariance55),
  fhBinLimClusterTPC(c.fhBinLimClusterTPC),
  fhBinLimClusterITS(c.fhBinLimClusterITS),
  fhBinLimChi2TPC(c.fhBinLimChi2TPC),
  fhBinLimChi2ITS(c.fhBinLimChi2ITS),
  fhBinLimCovariance11(c.fhBinLimCovariance11),
  fhBinLimCovariance22(c.fhBinLimCovariance22),
  fhBinLimCovariance33(c.fhBinLimCovariance33),
  fhBinLimCovariance44(c.fhBinLimCovariance44),
  fhBinLimCovariance55(c.fhBinLimCovariance55)
{
  //
  // copy constructor
  //
  ((AliCFTrackQualityCuts &) c).Copy(*this);
}
//__________________________________________________________________________________
AliCFTrackQualityCuts& AliCFTrackQualityCuts::operator=(const AliCFTrackQualityCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fMinNClusterTPC = c.fMinNClusterTPC ;
    fMinNClusterITS = c.fMinNClusterITS ;
    fMaxChi2PerClusterTPC = c.fMaxChi2PerClusterTPC ;
    fMaxChi2PerClusterITS = c.fMaxChi2PerClusterITS ;
    fCovariance11Max = c.fCovariance11Max ;
    fCovariance22Max = c.fCovariance22Max ;
    fCovariance33Max = c.fCovariance33Max ;
    fCovariance44Max = c.fCovariance44Max ;
    fCovariance55Max = c.fCovariance55Max ;
    fStatus = c.fStatus ;
    fhCutStatistics = c.fhCutStatistics ;
    fhCutCorrelation = c.fhCutCorrelation ;
    fBitmap =  c.fBitmap ;
    fhNBinsClusterTPC = c.fhNBinsClusterTPC ;
    fhNBinsClusterITS = c.fhNBinsClusterITS ;
    fhNBinsChi2TPC = c.fhNBinsChi2TPC ;
    fhNBinsChi2ITS = c.fhNBinsChi2ITS ;
    fhNBinsCovariance11 = c.fhNBinsCovariance11 ;
    fhNBinsCovariance22 = c.fhNBinsCovariance22 ;
    fhNBinsCovariance33 = c.fhNBinsCovariance33 ;
    fhNBinsCovariance44 = c.fhNBinsCovariance44 ;
    fhNBinsCovariance55 = c.fhNBinsCovariance55 ;
    fhBinLimClusterTPC = c.fhBinLimClusterTPC ;
    fhBinLimClusterITS = c.fhBinLimClusterITS ;
    fhBinLimChi2TPC = c.fhBinLimChi2TPC ;
    fhBinLimChi2ITS = c.fhBinLimChi2ITS ;
    fhBinLimCovariance11 = c.fhBinLimCovariance11 ;
    fhBinLimCovariance22 = c.fhBinLimCovariance22 ;
    fhBinLimCovariance33 = c.fhBinLimCovariance33 ;
    fhBinLimCovariance44 = c.fhBinLimCovariance44 ;
    fhBinLimCovariance55 = c.fhBinLimCovariance55 ;

    for (Int_t i=0; i<c.kNHist; i++){
      for (Int_t j=0; j<c.kNStepQA; j++){
	if(c.fhQA[i][j]) fhQA[i][j] = (TH1F*)c.fhQA[i][j]->Clone();
      }
    }
    ((AliCFTrackQualityCuts &) c).Copy(*this);
  }
  return *this;
}
//__________________________________________________________________________________
AliCFTrackQualityCuts::~AliCFTrackQualityCuts()
{
  //
  // destructor
  //
  if (fhCutStatistics)  delete fhCutStatistics;
  if (fhCutCorrelation) delete fhCutCorrelation;

  for (Int_t i=0; i<kNHist; i++){
    for (Int_t j=0; j<kNStepQA; j++){
      if(fhQA[i][j]) delete fhQA[i][j];
    }
  }
  if(fBitmap) delete fBitmap;
  if(fhBinLimClusterTPC) delete fhBinLimClusterTPC;
  if(fhBinLimClusterITS) delete fhBinLimClusterITS;
  if(fhBinLimChi2TPC) delete fhBinLimChi2TPC;
  if(fhBinLimChi2ITS) delete fhBinLimChi2ITS;
  if(fhBinLimCovariance11) delete fhBinLimCovariance11;
  if(fhBinLimCovariance22) delete fhBinLimCovariance22;
  if(fhBinLimCovariance33) delete fhBinLimCovariance33;
  if(fhBinLimCovariance44) delete fhBinLimCovariance44;
  if(fhBinLimCovariance55) delete fhBinLimCovariance55;
}
//__________________________________________________________________________________
void AliCFTrackQualityCuts::Initialise()
{
  //
  // sets everything to zero
  //
  fMinNClusterTPC = 0;
  fMinNClusterITS = 0;
  fMaxChi2PerClusterTPC = 0;
  fMaxChi2PerClusterITS = 0;
  fCovariance11Max = 0;
  fCovariance22Max = 0;
  fCovariance33Max = 0;
  fCovariance44Max = 0;
  fCovariance55Max = 0;
  fStatus = 0;

  SetMinNClusterTPC();
  SetMinNClusterITS();
  SetMaxChi2PerClusterTPC();
  SetMaxChi2PerClusterITS();
  SetMaxCovDiagonalElements();
  SetStatus();

  for (Int_t i=0; i<kNHist; i++){
    for (Int_t j=0; j<kNStepQA; j++)
      fhQA[i][j] = 0x0;
  }
  fhCutStatistics = 0;
  fhCutCorrelation = 0;
  fBitmap=new TBits(0);

  //set default bining for QA histograms
  SetHistogramBins(kCutClusterTPC,165,-0.5,164.5);
  SetHistogramBins(kCutClusterITS,8,-0.5,7.5);
  SetHistogramBins(kCutChi2TPC,500,0.,10.);
  SetHistogramBins(kCutChi2ITS,500,0.,10.);
  SetHistogramBins(kCutCovElement11,200,0.,20.);
  SetHistogramBins(kCutCovElement22,200,0.,20.);
  SetHistogramBins(kCutCovElement33,100,0.,1.);
  SetHistogramBins(kCutCovElement44,100,0.,5.);
  SetHistogramBins(kCutCovElement55,100,0.,5.);
}
//__________________________________________________________________________________
void AliCFTrackQualityCuts::Copy(TObject &c) const
{
  //
  // Copy function
  //
  AliCFTrackQualityCuts& target = (AliCFTrackQualityCuts &) c;

  target.Initialise();

  if (fhCutStatistics)  target.fhCutStatistics = (TH1F*) fhCutStatistics->Clone();
  if (fhCutCorrelation) target.fhCutCorrelation = (TH2F*) fhCutCorrelation->Clone();

  for (Int_t i=0; i<kNHist; i++){
    for (Int_t j=0; j<kNStepQA; j++){
      if(fhQA[i][j]) target.fhQA[i][j] = (TH1F*)fhQA[i][j]->Clone();
    }
  }
  TNamed::Copy(c);
}
//__________________________________________________________________________________
void AliCFTrackQualityCuts::SelectionBitMap(TObject* obj)
{
  //
  // test if the track passes the single cuts
  // and store the information in a bitmap
  //

  // bitmap stores the decision of each single cut
  for(Int_t i=0; i<kNCuts; i++)fBitmap->SetBitNumber(i,kFALSE);

  if (!obj) return;
  if (!obj->InheritsFrom("AliVParticle")) {
    AliError("object must derived from AliVParticle !");
    return;
  }

  Bool_t isESDTrack = strcmp(obj->ClassName(),"AliESDtrack") == 0 ? kTRUE : kFALSE ;
  Bool_t isAODTrack = strcmp(obj->ClassName(),"AliAODTrack") == 0 ? kTRUE : kFALSE ;

  AliESDtrack * esdTrack = 0x0 ;
  AliAODTrack * aodTrack = 0x0 ; 
  if (isESDTrack) esdTrack = dynamic_cast<AliESDtrack*>(obj);
  if (isAODTrack) aodTrack = dynamic_cast<AliAODTrack*>(obj);

  // get cut quantities
  Int_t    fIdxInt[200];
  Int_t    nClustersTPC = 0;
  Int_t    nClustersITS = 0 ;
  Float_t  chi2PerClusterTPC =  0 ;
  Float_t  chi2PerClusterITS = 0 ;
  Double_t extCov[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  if (isESDTrack) {
    nClustersTPC = esdTrack->GetTPCclusters(fIdxInt);
    nClustersITS = esdTrack->GetITSclusters(fIdxInt);
    if (nClustersTPC != 0) chi2PerClusterTPC = esdTrack->GetTPCchi2() / Float_t(nClustersTPC);
    if (nClustersITS != 0) chi2PerClusterITS = esdTrack->GetITSchi2() / Float_t(nClustersITS);
    esdTrack->GetExternalCovariance(extCov);
  }

  // fill the bitmap
  Int_t iCutBit = 0;

  if (nClustersTPC >= fMinNClusterTPC)
    fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;
  if (nClustersITS >= fMinNClusterITS)
    fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;
  if (chi2PerClusterTPC <= fMaxChi2PerClusterTPC)
    fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;
  if (chi2PerClusterITS <= fMaxChi2PerClusterITS)
    fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;
  if (extCov[0]  <= fCovariance11Max)
    fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;
  if (extCov[2]  <= fCovariance22Max)
    fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;
  if (extCov[5]  <= fCovariance33Max)
    fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;
  if (extCov[9]  <= fCovariance44Max)
    fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;
  if (extCov[14] <= fCovariance55Max)
    fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;
  
  if (isESDTrack) {
    if ( (esdTrack->GetStatus() & fStatus) == fStatus ) fBitmap->SetBitNumber(iCutBit,kTRUE);
  }
  else {
    if ( (aodTrack->GetStatus() & fStatus) == fStatus ) fBitmap->SetBitNumber(iCutBit,kTRUE);
  }

  return;
}
//__________________________________________________________________________________
Bool_t AliCFTrackQualityCuts::IsSelected(TObject* obj) {
  //
  // loops over decisions of single cuts and returns if the track is accepted
  //
  SelectionBitMap(obj);

  if (fIsQAOn) FillHistograms(obj,0);
  Bool_t isSelected = kTRUE;

  for (UInt_t icut=0; icut<fBitmap->GetNbits();icut++) {
    if(!fBitmap->TestBitNumber(icut)) {
	isSelected = kFALSE;
	break;
    }
  }
  if (!isSelected) return kFALSE ;
  if (fIsQAOn) FillHistograms(obj,1);
  return kTRUE;
}
//__________________________________________________________________________________
void AliCFTrackQualityCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins)
{
  //
  // variable bin size
  //
  if(!fIsQAOn) return;

  if (index<0 || index>=kNHist) {
    Error("SetHistogramBins","could not determine histogram from index %d",index);
    return;
  }

  switch(index){
  case kCutClusterTPC:
    fhNBinsClusterTPC=nbins+1;
    fhBinLimClusterTPC=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimClusterTPC[i]=bins[i];
    break;

  case kCutClusterITS:
    fhNBinsClusterITS=nbins+1;
    fhBinLimClusterITS=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimClusterITS[i]=bins[i];
    break;

  case kCutChi2TPC:
    fhNBinsChi2TPC=nbins+1;
    fhBinLimChi2TPC=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimChi2TPC[i]=bins[i];
    break;

  case kCutChi2ITS:
    fhNBinsChi2ITS=nbins+1;
    fhBinLimChi2ITS=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimChi2ITS[i]=bins[i];
    break;

  case kCutCovElement11:
    fhNBinsCovariance11=nbins+1;
    fhBinLimCovariance11=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimCovariance11[i]=bins[i];
    break;

  case kCutCovElement22:
    fhNBinsCovariance22=nbins+1;
    fhBinLimCovariance22=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimCovariance22[i]=bins[i];
    break;

  case kCutCovElement33:
    fhNBinsCovariance33=nbins+1;
    fhBinLimCovariance33=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimCovariance33[i]=bins[i];
    break;

  case kCutCovElement44:
    fhNBinsCovariance44=nbins+1;
    fhBinLimCovariance44=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimCovariance44[i]=bins[i];
    break;

  case kCutCovElement55:
    fhNBinsCovariance55=nbins+1;
    fhBinLimCovariance55=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimCovariance55[i]=bins[i];
    break;
 }
}
//__________________________________________________________________________________
void AliCFTrackQualityCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax)
{
  //
  // fixed bin size
  //
  if (index<0 || index>=kNHist) {
    Error("SetHistogramBins","could not determine histogram from index %d",index);
    return;
  }

  switch(index){
  case kCutClusterTPC:
    fhNBinsClusterTPC=nbins+1;
    fhBinLimClusterTPC=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimClusterTPC[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutClusterITS:
    fhNBinsClusterITS=nbins+1;
    fhBinLimClusterITS=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimClusterITS[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutChi2TPC:
    fhNBinsChi2TPC=nbins+1;
    fhBinLimChi2TPC=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimChi2TPC[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutChi2ITS:
    fhNBinsChi2ITS=nbins+1;
    fhBinLimChi2ITS=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimChi2ITS[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutCovElement11:
    fhNBinsCovariance11=nbins+1;
    fhBinLimCovariance11=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimCovariance11[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutCovElement22:
    fhNBinsCovariance22=nbins+1;
    fhBinLimCovariance22=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimCovariance22[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutCovElement33:
    fhNBinsCovariance33=nbins+1;
    fhBinLimCovariance33=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimCovariance33[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutCovElement44:
    fhNBinsCovariance44=nbins+1;
    fhBinLimCovariance44=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimCovariance44[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutCovElement55:
    fhNBinsCovariance55=nbins+1;
    fhBinLimCovariance55=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimCovariance55[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
 }
}
//__________________________________________________________________________________
 void AliCFTrackQualityCuts::DefineHistograms() {
  //
  // histograms for cut variables, cut statistics and cut correlations
  //
  Int_t color = 2;

  // book cut statistics and cut correlation histograms
  fhCutStatistics = new TH1F(Form("%s_cut_statistics",GetName()), Form("%s cut statistics",GetName()), kNCuts,0.5,kNCuts+0.5);
  fhCutStatistics->SetLineWidth(2);
  fhCutStatistics->GetXaxis()->SetBinLabel(1,"nClustersTPC");
  fhCutStatistics->GetXaxis()->SetBinLabel(2,"nClustersITS");
  fhCutStatistics->GetXaxis()->SetBinLabel(3,"chi2PerClusterTPC");
  fhCutStatistics->GetXaxis()->SetBinLabel(4,"chi2PerClusterITS");
  fhCutStatistics->GetXaxis()->SetBinLabel(5,"covElement11");
  fhCutStatistics->GetXaxis()->SetBinLabel(6,"covElement22");
  fhCutStatistics->GetXaxis()->SetBinLabel(7,"covElement33");
  fhCutStatistics->GetXaxis()->SetBinLabel(8,"covElement44");
  fhCutStatistics->GetXaxis()->SetBinLabel(9,"covElement55");

  fhCutCorrelation = new TH2F(Form("%s_cut_correlation",GetName()), Form("%s cut  correlation",GetName()), kNCuts,0.5,kNCuts+0.5,kNCuts,0.5,kNCuts+0.5);
  fhCutCorrelation->SetLineWidth(2);
  fhCutCorrelation->GetXaxis()->SetBinLabel(1,"nClustersTPC");
  fhCutCorrelation->GetXaxis()->SetBinLabel(2,"nClustersITS");
  fhCutCorrelation->GetXaxis()->SetBinLabel(3,"chi2PerClusterTPC");
  fhCutCorrelation->GetXaxis()->SetBinLabel(4,"chi2PerClusterITS");
  fhCutCorrelation->GetXaxis()->SetBinLabel(5,"covElement11");
  fhCutCorrelation->GetXaxis()->SetBinLabel(6,"covElement22");
  fhCutCorrelation->GetXaxis()->SetBinLabel(7,"covElement33");
  fhCutCorrelation->GetXaxis()->SetBinLabel(8,"covElement44");
  fhCutCorrelation->GetXaxis()->SetBinLabel(9,"covElement55");

  fhCutCorrelation->GetYaxis()->SetBinLabel(1,"nClustersTPC");
  fhCutCorrelation->GetYaxis()->SetBinLabel(2,"nClustersITS");
  fhCutCorrelation->GetYaxis()->SetBinLabel(3,"chi2PerClusterTPC");
  fhCutCorrelation->GetYaxis()->SetBinLabel(4,"chi2PerClusterITS");
  fhCutCorrelation->GetYaxis()->SetBinLabel(5,"covElement11");
  fhCutCorrelation->GetYaxis()->SetBinLabel(6,"covElement22");
  fhCutCorrelation->GetYaxis()->SetBinLabel(7,"covElement33");
  fhCutCorrelation->GetYaxis()->SetBinLabel(8,"covElement44");
  fhCutCorrelation->GetYaxis()->SetBinLabel(9,"covElement55");


  // book QA histograms
  Char_t str[256];
  for (Int_t i=0; i<kNStepQA; i++) {
    if (i==0) sprintf(str," ");
    else sprintf(str,"_cut");

    fhQA[kCutClusterTPC][i]	= new TH1F(Form("%s_nClustersTPC%s",GetName(),str)     ,"",fhNBinsClusterTPC-1,fhBinLimClusterTPC);
    fhQA[kCutClusterITS][i]	= new TH1F(Form("%s_nClustersITS%s",GetName(),str)     ,"",fhNBinsClusterITS-1,fhBinLimClusterITS);
    fhQA[kCutChi2TPC][i]	= new TH1F(Form("%s_chi2PerClusterTPC%s",GetName(),str),"",fhNBinsChi2TPC-1,fhBinLimChi2TPC);
    fhQA[kCutChi2ITS][i]	= new TH1F(Form("%s_chi2PerClusterITS%s",GetName(),str),"",fhNBinsChi2ITS-1,fhBinLimChi2ITS);
    fhQA[kCutCovElement11][i]	= new TH1F(Form("%s_covMatrixDiagonal11%s",GetName(),str),"",fhNBinsCovariance11-1,fhBinLimCovariance11);
    fhQA[kCutCovElement22][i]	= new TH1F(Form("%s_covMatrixDiagonal22%s",GetName(),str),"",fhNBinsCovariance22-1,fhBinLimCovariance22);
    fhQA[kCutCovElement33][i]	= new TH1F(Form("%s_covMatrixDiagonal33%s",GetName(),str),"",fhNBinsCovariance33-1,fhBinLimCovariance33);
    fhQA[kCutCovElement44][i]	= new TH1F(Form("%s_covMatrixDiagonal44%s",GetName(),str),"",fhNBinsCovariance44-1,fhBinLimCovariance44);
    fhQA[kCutCovElement55][i]	= new TH1F(Form("%s_covMatrixDiagonal55%s",GetName(),str),"",fhNBinsCovariance55-1,fhBinLimCovariance55);

    fhQA[kCutClusterTPC][i]	->SetXTitle("n TPC clusters");
    fhQA[kCutClusterITS][i]	->SetXTitle("n ITS clusters");
    fhQA[kCutChi2TPC][i]	->SetXTitle("#chi^{2} per TPC cluster");
    fhQA[kCutChi2ITS][i]	->SetXTitle("#chi^{2} per ITS cluster");
    fhQA[kCutCovElement11][i]	->SetXTitle("cov 11 : #sigma_{y}^{2} (cm^{2})");
    fhQA[kCutCovElement22][i]	->SetXTitle("cov 22 : #sigma_{z}^{2} (cm^{2})");
    fhQA[kCutCovElement33][i]	->SetXTitle("cov 33 : #sigma_{sin(#phi)}^{2}");
    fhQA[kCutCovElement44][i]	->SetXTitle("cov 44 : #sigma_{tan(#theta_{dip})}^{2}");
    fhQA[kCutCovElement55][i]	->SetXTitle("cov 55 : #sigma_{1/p_{T}}^{2} ((c/GeV)^{2})");
  }

  for(Int_t i=0; i<kNHist; i++) fhQA[i][1]->SetLineColor(color);
}
//__________________________________________________________________________________
void AliCFTrackQualityCuts::FillHistograms(TObject* obj, Bool_t b)
{
  //
  // fill the QA histograms
  //

  if (!obj) return;
  if (!obj->InheritsFrom("AliVParticle")) {
    AliError("object must derived from AliVParticle !");
    return;
  }

  Bool_t isESDTrack = strcmp(obj->ClassName(),"AliESDtrack") == 0 ? kTRUE : kFALSE ;
  Bool_t isAODTrack = strcmp(obj->ClassName(),"AliAODTrack") == 0 ? kTRUE : kFALSE ;

  AliESDtrack * esdTrack = 0x0 ;
  AliAODTrack * aodTrack = 0x0 ; 
  if (isESDTrack) esdTrack = dynamic_cast<AliESDtrack*>(obj);
  if (isAODTrack) aodTrack = dynamic_cast<AliAODTrack*>(obj);

  // b = 0: fill histograms before cuts
  // b = 1: fill histograms after cuts

  Int_t    fIdxInt[200];
  Int_t    nClustersTPC = 0;
  Int_t    nClustersITS = 0 ;
  Float_t  chi2PerClusterTPC =  0 ;
  Float_t  chi2PerClusterITS = 0 ;
  Double_t extCov[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  if (isESDTrack) {
    nClustersTPC = esdTrack->GetTPCclusters(fIdxInt);
    nClustersITS = esdTrack->GetITSclusters(fIdxInt);
    if (nClustersTPC != 0) chi2PerClusterTPC = esdTrack->GetTPCchi2() / Float_t(nClustersTPC);
    if (nClustersITS!=0) chi2PerClusterITS = esdTrack->GetITSchi2() / Float_t(nClustersITS);
    esdTrack->GetExternalCovariance(extCov);
  }
  fhQA[kCutClusterTPC][b]->Fill((float)nClustersTPC);
  fhQA[kCutChi2TPC][b]->Fill(chi2PerClusterTPC);
  fhQA[kCutClusterITS][b]->Fill((float)nClustersITS);
  fhQA[kCutChi2ITS][b]->Fill(chi2PerClusterITS);
  fhQA[kCutCovElement11][b]->Fill(extCov[0]);
  fhQA[kCutCovElement22][b]->Fill(extCov[2]);
  fhQA[kCutCovElement33][b]->Fill(extCov[5]);
  fhQA[kCutCovElement44][b]->Fill(extCov[9]);
  fhQA[kCutCovElement55][b]->Fill(extCov[14]);

  // fill cut statistics and cut correlation histograms with information from the bitmap
  if (b) return;

  // Get the bitmap of the single cuts
  if ( !obj ) return;
  SelectionBitMap(obj);

  // Number of single cuts in this class
  UInt_t ncuts = fBitmap->GetNbits();
  for(UInt_t bit=0; bit<ncuts;bit++) {
    if (!fBitmap->TestBitNumber(bit)) {
	fhCutStatistics->Fill(bit+1);
	for (UInt_t bit2=bit; bit2<ncuts;bit2++) {
	  if (!fBitmap->TestBitNumber(bit2)) 
	    fhCutCorrelation->Fill(bit+1,bit2+1);
	}
    }
  }
}
//__________________________________________________________________________________
void AliCFTrackQualityCuts::SaveHistograms(const Char_t* dir) {
  //
  // saves the histograms in a directory (dir)
  //
  if(!fIsQAOn) return;

  if (!dir)
    dir = GetName();

  gDirectory->mkdir(dir);
  gDirectory->cd(dir);

  gDirectory->mkdir("before_cuts");
  gDirectory->mkdir("after_cuts");

  fhCutStatistics->Write();
  fhCutCorrelation->Write();

  for (Int_t j=0; j<kNStepQA; j++) {
    if (j==0)
      gDirectory->cd("before_cuts");
    else
      gDirectory->cd("after_cuts");

    for(Int_t i=0; i<kNHist; i++) fhQA[i][j]->Write();

    gDirectory->cd("../");
  }
  gDirectory->cd("../");
}
//__________________________________________________________________________________
void AliCFTrackQualityCuts::DrawHistograms(Bool_t drawLogScale)
{
  //
  // draws some histograms
  //
  if(!fIsQAOn) return;

  // pad margins
  Float_t right = 0.03;
  Float_t left = 0.175;
  Float_t top = 0.03;
  Float_t bottom = 0.175;

  TCanvas* canvas1 = new TCanvas("Track_QA_Quality_1", "Track QA Quality 1", 800, 500);
  canvas1->Divide(2, 1);

  canvas1->cd(1);
  fhCutStatistics->SetStats(kFALSE);
  fhCutStatistics->LabelsOption("v");
  gPad->SetLeftMargin(left);
  gPad->SetBottomMargin(0.25);
  gPad->SetRightMargin(right);
  gPad->SetTopMargin(0.1);
  fhCutStatistics->Draw();

  canvas1->cd(2);
  fhCutCorrelation->SetStats(kFALSE);
  fhCutCorrelation->LabelsOption("v");
  gPad->SetLeftMargin(0.30);
  gPad->SetRightMargin(bottom);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.25);
  fhCutCorrelation->Draw("COLZ");

  canvas1->SaveAs(Form("%s.eps", canvas1->GetName()));
  canvas1->SaveAs(Form("%s.ps", canvas1->GetName()));

  // -----

  TCanvas* canvas2 = new TCanvas("Track_QA_Quality_2", "Track QA Quality 2", 1200, 800);
  canvas2->Divide(2, 2);

  canvas2->cd(1);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutClusterTPC][0]->SetStats(kFALSE);
  fhQA[kCutClusterTPC][0]->Draw();
  fhQA[kCutClusterTPC][1]->Draw("same");

  canvas2->cd(2);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutChi2TPC][0]->SetStats(kFALSE);
  fhQA[kCutChi2TPC][0]->Draw();
  fhQA[kCutChi2TPC][1]->Draw("same");

  canvas2->cd(3);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutClusterITS][0]->SetStats(kFALSE);
  fhQA[kCutClusterITS][0]->Draw();
  fhQA[kCutClusterITS][1]->Draw("same");

  canvas2->cd(4);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutChi2ITS][0]->SetStats(kFALSE);
  fhQA[kCutChi2ITS][0]->Draw();
  fhQA[kCutChi2ITS][1]->Draw("same");

  canvas2->SaveAs(Form("%s.eps", canvas2->GetName()));
  canvas2->SaveAs(Form("%s.ps", canvas2->GetName()));

  // -----

  TCanvas* canvas3 = new TCanvas("Track_QA_Quality_3", "Track QA Quality 3", 1200, 800);
  canvas3->Divide(3, 2);

  canvas3->cd(1);
  if(drawLogScale) gPad->SetLogy();
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutCovElement11][0]->SetStats(kFALSE);
  fhQA[kCutCovElement11][0]->Draw();
  fhQA[kCutCovElement11][1]->Draw("same");

  canvas3->cd(2);
  if(drawLogScale) gPad->SetLogy();
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutCovElement22][0]->SetStats(kFALSE);
  fhQA[kCutCovElement22][0]->Draw();
  fhQA[kCutCovElement22][1]->Draw("same");

  canvas3->cd(3);
  if(drawLogScale) gPad->SetLogy();
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutCovElement33][0]->SetStats(kFALSE);
  fhQA[kCutCovElement33][0]->Draw();
  fhQA[kCutCovElement33][1]->Draw("same");

  canvas3->cd(4);
  if(drawLogScale) gPad->SetLogy();
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutCovElement44][0]->SetStats(kFALSE);
  fhQA[kCutCovElement44][0]->Draw();
  fhQA[kCutCovElement44][1]->Draw("same");

  canvas3->cd(5);
  if(drawLogScale) gPad->SetLogy();
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutCovElement55][0]->SetStats(kFALSE);
  fhQA[kCutCovElement55][0]->Draw();
  fhQA[kCutCovElement55][1]->Draw("same");

  canvas3->SaveAs(Form("%s.eps", canvas3->GetName()));
  canvas3->SaveAs(Form("%s.ps", canvas3->GetName()));
}
//__________________________________________________________________________________
void AliCFTrackQualityCuts::AddQAHistograms(TList *qaList) {
  //
  // saves the histograms in a TList
  //
  DefineHistograms();

  qaList->Add(fhCutStatistics);
  qaList->Add(fhCutCorrelation);

  for (Int_t j=0; j<kNStepQA; j++) {
    for(Int_t i=0; i<kNHist; i++)
	qaList->Add(fhQA[i][j]);
  }
}
