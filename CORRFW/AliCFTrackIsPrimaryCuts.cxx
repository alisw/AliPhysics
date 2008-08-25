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

// The class AliCFTrackIsPrimaryCut is designed to select reconstructed tracks
// with a small impact parameter and tracks which are (not) daughters of kink
// decays and to provide corresponding QA histograms.
// This class inherits from the Analysis' Framework abstract base class
// AliAnalysisCuts and is a part of the Correction Framework.
// This class acts on single, reconstructed tracks, it is applicable on
// ESD and AOD data.
// It mainly consists of a IsSelected function that returns a boolean.
// This function checks whether the considered track passes a set of cuts:
// - distance to main vertex in units of sigma (resolution)
// - require that the dca calculation doesn't fail
// - accept or not accept daughter tracks of kink decays
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
#include "AliCFTrackIsPrimaryCuts.h"

ClassImp(AliCFTrackIsPrimaryCuts)

//__________________________________________________________________________________
AliCFTrackIsPrimaryCuts::AliCFTrackIsPrimaryCuts() :
  AliCFCutBase(),
  fNSigmaToVertex(0),
  fNSigmaToVertexMax(0),
  fRequireSigmaToVertex(0),
  fAODType(AliAODTrack::kUndef),
  fAcceptKinkDaughters(0),
  fhCutStatistics(0),
  fhCutCorrelation(0),
  fBitmap(0x0),
  fhNBinsNSigma(0),
  fhNBinsRequireSigma(0),
  fhNBinsAcceptKink(0),
  fhNBinsDcaXY(0),
  fhNBinsDcaZ(0),
  fhNBinsDcaXYnorm(0),
  fhNBinsDcaZnorm(0),
  fhBinLimNSigma(0x0),
  fhBinLimRequireSigma(0x0),
  fhBinLimAcceptKink(0x0),
  fhBinLimDcaXY(0x0),
  fhBinLimDcaZ(0x0),
  fhBinLimDcaXYnorm(0x0),
  fhBinLimDcaZnorm(0x0)
{
  //
  // Default constructor
  //
  Initialise();
}
//__________________________________________________________________________________
AliCFTrackIsPrimaryCuts::AliCFTrackIsPrimaryCuts(Char_t* name, Char_t* title) :
  AliCFCutBase(name,title),
  fNSigmaToVertex(0),
  fNSigmaToVertexMax(0),
  fRequireSigmaToVertex(0),
  fAODType(AliAODTrack::kUndef),
  fAcceptKinkDaughters(0),
  fhCutStatistics(0),
  fhCutCorrelation(0),
  fBitmap(0x0),
  fhNBinsNSigma(0),
  fhNBinsRequireSigma(0),
  fhNBinsAcceptKink(0),
  fhNBinsDcaXY(0),
  fhNBinsDcaZ(0),
  fhNBinsDcaXYnorm(0),
  fhNBinsDcaZnorm(0),
  fhBinLimNSigma(0x0),
  fhBinLimRequireSigma(0x0),
  fhBinLimAcceptKink(0x0),
  fhBinLimDcaXY(0x0),
  fhBinLimDcaZ(0x0),
  fhBinLimDcaXYnorm(0x0),
  fhBinLimDcaZnorm(0x0)
{
  //
  // Constructor
  //
  Initialise();
}
//__________________________________________________________________________________
AliCFTrackIsPrimaryCuts::AliCFTrackIsPrimaryCuts(const AliCFTrackIsPrimaryCuts& c) :
  AliCFCutBase(c),
  fNSigmaToVertex(c.fNSigmaToVertex),
  fNSigmaToVertexMax(c.fNSigmaToVertexMax),
  fRequireSigmaToVertex(c.fRequireSigmaToVertex),
  fAODType(c.fAODType),
  fAcceptKinkDaughters(c.fAcceptKinkDaughters),
  fhCutStatistics(c.fhCutStatistics),
  fhCutCorrelation(c.fhCutCorrelation),
  fBitmap(c.fBitmap),
  fhNBinsNSigma(c.fhNBinsNSigma),
  fhNBinsRequireSigma(c.fhNBinsRequireSigma),
  fhNBinsAcceptKink(c.fhNBinsAcceptKink),
  fhNBinsDcaXY(c.fhNBinsDcaXY),
  fhNBinsDcaZ(c.fhNBinsDcaZ),
  fhNBinsDcaXYnorm(c.fhNBinsDcaXYnorm),
  fhNBinsDcaZnorm(c.fhNBinsDcaZnorm),
  fhBinLimNSigma(c.fhBinLimNSigma),
  fhBinLimRequireSigma(c.fhBinLimRequireSigma),
  fhBinLimAcceptKink(c.fhBinLimAcceptKink),
  fhBinLimDcaXY(c.fhBinLimDcaXY),
  fhBinLimDcaZ(c.fhBinLimDcaZ),
  fhBinLimDcaXYnorm(c.fhBinLimDcaXYnorm),
  fhBinLimDcaZnorm(c.fhBinLimDcaZnorm)
{
  //
  // copy constructor
  //
  ((AliCFTrackIsPrimaryCuts &) c).Copy(*this);
}
//__________________________________________________________________________________
AliCFTrackIsPrimaryCuts& AliCFTrackIsPrimaryCuts::operator=(const AliCFTrackIsPrimaryCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fNSigmaToVertex = c.fNSigmaToVertex ;
    fNSigmaToVertexMax = c.fNSigmaToVertexMax ;
    fRequireSigmaToVertex = c.fRequireSigmaToVertex ;
    fAODType = c.fAODType ;
    fAcceptKinkDaughters = c.fAcceptKinkDaughters ;
    fhCutStatistics = c.fhCutStatistics ;
    fhCutCorrelation = c.fhCutCorrelation ;
    fBitmap =  c.fBitmap;
    fhNBinsNSigma = c.fhNBinsNSigma;
    fhNBinsRequireSigma = c.fhNBinsRequireSigma;
    fhNBinsAcceptKink = c.fhNBinsAcceptKink;
    fhNBinsDcaXY = c.fhNBinsDcaXY;
    fhNBinsDcaZ = c.fhNBinsDcaZ;
    fhNBinsDcaXYnorm = c.fhNBinsDcaXYnorm;
    fhNBinsDcaZnorm = c.fhNBinsDcaZnorm;
    fhBinLimNSigma = c.fhBinLimNSigma;
    fhBinLimRequireSigma = c.fhBinLimRequireSigma;
    fhBinLimAcceptKink = c.fhBinLimAcceptKink;
    fhBinLimDcaXY = c.fhBinLimDcaXY;
    fhBinLimDcaZ = c.fhBinLimDcaZ;
    fhBinLimDcaXYnorm = c.fhBinLimDcaXYnorm;
    fhBinLimDcaZnorm = c.fhBinLimDcaZnorm;

    for (Int_t j=0; j<c.kNStepQA; j++){
      if(c.fhDcaXYvsDcaZ[j]) fhDcaXYvsDcaZ[j] = (TH2F*)c.fhDcaXYvsDcaZ[j]->Clone();
      if(c.fhDcaXYvsDcaZnorm[j]) fhDcaXYvsDcaZnorm[j] = (TH2F*)c.fhDcaXYvsDcaZnorm[j]->Clone();
      for (Int_t i=0; i<c.kNHist; i++){
	if(c.fhQA[i][j]) fhQA[i][j] = (TH1F*)c.fhQA[i][j]->Clone();
      }
    }
    ((AliCFTrackIsPrimaryCuts &) c).Copy(*this);
 }
  return *this;
}
//__________________________________________________________________________________
AliCFTrackIsPrimaryCuts::~AliCFTrackIsPrimaryCuts()
{
  //
  // destructor
  //
  if (fhCutStatistics)			delete fhCutStatistics;
  if (fhCutCorrelation)			delete fhCutCorrelation;

  for (Int_t j=0; j<kNStepQA; j++){
    if(fhDcaXYvsDcaZ[j])	delete fhDcaXYvsDcaZ[j];
    if(fhDcaXYvsDcaZnorm[j])	delete fhDcaXYvsDcaZnorm[j];
    for (Int_t i=0; i<kNHist; i++)
      if(fhQA[i][j]) 		delete fhQA[i][j];
  }
  if(fBitmap) 			delete fBitmap;
  if(fhBinLimNSigma) 		delete fhBinLimNSigma;
  if(fhBinLimRequireSigma) 	delete fhBinLimRequireSigma;
  if(fhBinLimAcceptKink) 	delete fhBinLimAcceptKink;
  if(fhBinLimDcaXY) 		delete fhBinLimDcaXY;
  if(fhBinLimDcaZ) 		delete fhBinLimDcaZ;
  if(fhBinLimDcaXYnorm) 	delete fhBinLimDcaXYnorm;
  if(fhBinLimDcaZnorm) 		delete fhBinLimDcaZnorm;
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::Initialise()
{
  //
  // sets everything to zero
  //
  fNSigmaToVertex = 0;
  fNSigmaToVertexMax = 0;
  fRequireSigmaToVertex = 0;
  fAcceptKinkDaughters = 0;
  fAODType = AliAODTrack::kUndef;

  SetMaxNSigmaToVertex();
  SetRequireSigmaToVertex();
  SetAcceptKinkDaughters();
  SetAODType();

  for (Int_t j=0; j<kNStepQA; j++)  {
    fhDcaXYvsDcaZ[j] = 0x0;
    fhDcaXYvsDcaZnorm[j] = 0x0;
    for (Int_t i=0; i<kNHist; i++)
      fhQA[i][j] = 0x0;
  }
  fhCutStatistics = 0;
  fhCutCorrelation = 0;
  fBitmap=new TBits(0);

  //set default bining for QA histograms
  SetHistogramBins(kCutNSigmaToVertex,500,0.,50.);
  SetHistogramBins(kCutRequireSigmaToVertex,5,-0.75,1.75);
  SetHistogramBins(kCutAcceptKinkDaughters,5,-0.75,1.75);
  SetHistogramBins(kDcaXY,500,-10.,10.);
  SetHistogramBins(kDcaZ,500,-10.,10.);
  SetHistogramBins(kDcaXYnorm,500,-10.,10.);
  SetHistogramBins(kDcaZnorm,500,-10.,10.);
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::Copy(TObject &c) const
{
  //
  // Copy function
  //
  AliCFTrackIsPrimaryCuts& target = (AliCFTrackIsPrimaryCuts &) c;

  target.Initialise();

  if (fhCutStatistics)  target.fhCutStatistics = (TH1F*) fhCutStatistics->Clone();
  if (fhCutCorrelation) target.fhCutCorrelation = (TH2F*) fhCutCorrelation->Clone();

  for (Int_t j=0; j<kNStepQA; j++){
    if(fhDcaXYvsDcaZ[j]) target.fhDcaXYvsDcaZ[j] = (TH2F*)fhDcaXYvsDcaZ[j]->Clone();
    if(fhDcaXYvsDcaZnorm[j]) target.fhDcaXYvsDcaZnorm[j] = (TH2F*)fhDcaXYvsDcaZnorm[j]->Clone();
    for (Int_t i=0; i<kNHist; i++)
      if(fhQA[i][j]) target.fhQA[i][j] = (TH1F*)fhQA[i][j]->Clone();
  }
  TNamed::Copy(c);
}
//____________________________________________________________________
void AliCFTrackIsPrimaryCuts::GetSigmaToVertex(AliESDtrack* esdTrack)
{
  //
  // Calculates the number of sigma to the vertex.
  // Currently applicable for ESD tracks only
  //
  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  esdTrack->GetImpactParameters(b,bCov);
  if (bCov[0]<=0 || bCov[2]<=0) {
    AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);

  // -----------------------------------
  // How to get to a n-sigma cut?
  //
  // The accumulated statistics from 0 to d is
  //
  // ->  Erf(d/Sqrt(2)) for a 1-dim gauss (d = n_sigma)
  // ->  1 - Exp(-d**2) for a 2-dim gauss (d*d = dx*dx + dy*dy != n_sigma)
  //
  // It means that for a 2-dim gauss: n_sigma(d) = Sqrt(2)*ErfInv(1 - Exp((-x**2)/2)
  // Can this be expressed in a different way?

  if (bRes[0] == 0 || bRes[1] ==0) {
    fNSigmaToVertex = -1;
    return;
  }

  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));

  // stupid rounding problem screws up everything:
  // if d is too big, TMath::Exp(...) gets 0, and TMath::ErfInverse(1) that should be infinite, gets 0 :(
  if (TMath::Exp(-d * d / 2) < 1e-10) {
    fNSigmaToVertex = 1000;
    return;
  }

  d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  fNSigmaToVertex = d;
  return;
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::SelectionBitMap(TObject* obj)
{
  //
  // test if the track passes the single cuts
  // and store the information in a bitmap
  //

  // bitmap stores the decision of each single cut
  for(Int_t i=0; i<kNCuts; i++)fBitmap->SetBitNumber(i,kFALSE);

  // check TObject and cast into ESDtrack
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

  // get the track to vertex parameter for ESD track
  if (isESDTrack) GetSigmaToVertex(esdTrack);

  // fill the bitmap
  Int_t iCutBit = 0;

  if (isESDTrack) {
    if (fNSigmaToVertex <= fNSigmaToVertexMax) {
      fBitmap->SetBitNumber(iCutBit,kTRUE);
    }
  }
  else fBitmap->SetBitNumber(iCutBit,kTRUE);

  iCutBit++;

  if (isESDTrack) {
    if (!fRequireSigmaToVertex || (fNSigmaToVertex>=0 && fRequireSigmaToVertex)) {
      fBitmap->SetBitNumber(iCutBit,kTRUE);
    }
  }
  else fBitmap->SetBitNumber(iCutBit,kTRUE);

  iCutBit++;

  if (esdTrack) {
    if (fAcceptKinkDaughters || (!fAcceptKinkDaughters && esdTrack->GetKinkIndex(0)<=0)) {
      fBitmap->SetBitNumber(iCutBit,kTRUE);
    }
  }
  else fBitmap->SetBitNumber(iCutBit,kTRUE);
  
  iCutBit++;
  
  if (isAODTrack) {
    if (fAODType==AliAODTrack::kUndef || fAODType == aodTrack->GetType()) {
      fBitmap->SetBitNumber(iCutBit,kTRUE);
    }
  }
  else fBitmap->SetBitNumber(iCutBit,kTRUE);
  
  return;
}
//__________________________________________________________________________________
Bool_t AliCFTrackIsPrimaryCuts::IsSelected(TObject* obj) {
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
void AliCFTrackIsPrimaryCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins)
{
  //
  // variable bin size
  //
  if(!fIsQAOn) return;

  switch(index){
  case kCutNSigmaToVertex:
    fhNBinsNSigma=nbins+1;
    fhBinLimNSigma=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimNSigma[i]=bins[i];
    break;

  case kCutRequireSigmaToVertex:
    fhNBinsRequireSigma=nbins+1;
    fhBinLimRequireSigma=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimRequireSigma[i]=bins[i];
    break;

  case kCutAcceptKinkDaughters:
    fhNBinsAcceptKink=nbins+1;
    fhBinLimAcceptKink=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimAcceptKink[i]=bins[i];
    break;

  case kDcaXY:
    fhNBinsDcaXY=nbins+1;
    fhBinLimDcaXY=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaXY[i]=bins[i];
    break;

  case kDcaZ:
    fhNBinsDcaZ=nbins+1;
    fhBinLimDcaZ=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaZ[i]=bins[i];
    break;

  case kDcaXYnorm:
    fhNBinsDcaXYnorm=nbins+1;
    fhBinLimDcaXYnorm=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaXYnorm[i]=bins[i];
    break;

  case kDcaZnorm:
    fhNBinsDcaZnorm=nbins+1;
    fhBinLimDcaZnorm=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaZnorm[i]=bins[i];
    break;
  }
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax)
{
  //
  // fixed bin size
  //
  switch(index){
  case kCutNSigmaToVertex:
    fhNBinsNSigma=nbins+1;
    fhBinLimNSigma=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimNSigma[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutRequireSigmaToVertex:
    fhNBinsRequireSigma=nbins+1;
    fhBinLimRequireSigma=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimRequireSigma[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutAcceptKinkDaughters:
    fhNBinsAcceptKink=nbins+1;
    fhBinLimAcceptKink=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimAcceptKink[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kDcaXY:
    fhNBinsDcaXY=nbins+1;
    fhBinLimDcaXY=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaXY[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kDcaZ:
    fhNBinsDcaZ=nbins+1;
    fhBinLimDcaZ=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaZ[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kDcaXYnorm:
    fhNBinsDcaXYnorm=nbins+1;
    fhBinLimDcaXYnorm=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaXYnorm[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kDcaZnorm:
    fhNBinsDcaZnorm=nbins+1;
    fhBinLimDcaZnorm=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaZnorm[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
  }
}
//__________________________________________________________________________________
 void AliCFTrackIsPrimaryCuts::DefineHistograms() {
  //
  // histograms for cut variables, cut statistics and cut correlations
  //
  Int_t color = 2;

  // book cut statistics and cut correlation histograms
  fhCutStatistics = new TH1F(Form("%s_cut_statistics",GetName()), Form("%s cut statistics",GetName()), kNCuts,0.5,kNCuts+0.5);
  fhCutStatistics->SetLineWidth(2);
  fhCutStatistics->GetXaxis()->SetBinLabel(1,"n dca");
  fhCutStatistics->GetXaxis()->SetBinLabel(2,"require dca");
  fhCutStatistics->GetXaxis()->SetBinLabel(3,"kink daughter");
  fhCutStatistics->GetXaxis()->SetBinLabel(4,"AOD type");

  fhCutCorrelation = new TH2F(Form("%s_cut_correlation",GetName()), Form("%s cut  correlation",GetName()), kNCuts,0.5,kNCuts+0.5,kNCuts,0.5,kNCuts+0.5);
  fhCutCorrelation->SetLineWidth(2);
  fhCutCorrelation->GetXaxis()->SetBinLabel(1,"n dca");
  fhCutCorrelation->GetXaxis()->SetBinLabel(2,"require dca");
  fhCutCorrelation->GetXaxis()->SetBinLabel(3,"kink daughter");
  fhCutCorrelation->GetXaxis()->SetBinLabel(4,"AOD type");

  fhCutCorrelation->GetYaxis()->SetBinLabel(1,"n dca");
  fhCutCorrelation->GetYaxis()->SetBinLabel(2,"require dca");
  fhCutCorrelation->GetYaxis()->SetBinLabel(3,"kink daughter");
  fhCutCorrelation->GetYaxis()->SetBinLabel(4,"AOD type");

  // book QA histograms
  Char_t str[256];
  for (Int_t i=0; i<kNStepQA; i++) {
    if (i==0) sprintf(str," ");
    else sprintf(str,"_cut");

    fhDcaXYvsDcaZ[i]            = new  TH2F(Form("%s_dcaXYvsDcaZ%s",GetName(),str),"",200,-10,10,200,-10,10);
    fhDcaXYvsDcaZnorm[i]        = new  TH2F(Form("%s_dcaXYvsDcaZnorm%s",GetName(),str),"",200,-10,10,200,-10,10);
    fhQA[kCutNSigmaToVertex][i]	= new TH1F(Form("%s_nSigmaToVertex%s",GetName(),str),"",fhNBinsNSigma-1,fhBinLimNSigma);
    fhQA[kCutRequireSigmaToVertex][i] = new TH1F(Form("%s_requireSigmaToVertex%s",GetName(),str),"",fhNBinsRequireSigma-1,fhBinLimRequireSigma);
    fhQA[kCutAcceptKinkDaughters][i] = new TH1F(Form("%s_acceptKinkDaughters%s",GetName(),str),"",fhNBinsAcceptKink-1,fhBinLimAcceptKink);
    fhQA[kDcaXY][i]		= new TH1F(Form("%s_dcaXY%s",GetName(),str),"",fhNBinsDcaXY-1,fhBinLimDcaXY);
    fhQA[kDcaZ][i]		= new TH1F(Form("%s_dcaZ%s",GetName(),str),"",fhNBinsDcaZ-1,fhBinLimDcaZ);
    fhQA[kDcaXYnorm][i]		= new TH1F(Form("%s_dcaXYnorm%s",GetName(),str),"",fhNBinsDcaXYnorm-1,fhBinLimDcaXYnorm);
    fhQA[kDcaZnorm][i]		= new TH1F(Form("%s_dcaZnorm%s",GetName(),str),"",fhNBinsDcaZnorm-1,fhBinLimDcaZnorm);

    fhDcaXYvsDcaZ[i]->SetXTitle("impact par. d_{z}");
    fhDcaXYvsDcaZ[i]->SetYTitle("impact par. d_{xy}");
    fhDcaXYvsDcaZnorm[i]->SetXTitle("norm. impact par. d_{z} / #sigma_{z}");
    fhDcaXYvsDcaZnorm[i]->SetYTitle("norm. impact par. d_{xy} / #sigma_{xy}");

    fhQA[kCutNSigmaToVertex][i]->SetXTitle("n #sigma to vertex");
    fhQA[kCutRequireSigmaToVertex][i]->SetXTitle("require #sigma to vertex");
    fhQA[kCutAcceptKinkDaughters][i]->SetXTitle("accept kink daughters");
    fhQA[kDcaXY][i]->SetXTitle("impact par. d_{xy}");
    fhQA[kDcaZ][i]->SetXTitle("impact par. d_{z}");
    fhQA[kDcaXYnorm][i]->SetXTitle("norm. impact par. d_{xy} / #sigma_{xy}");
    fhQA[kDcaZnorm][i]->SetXTitle("norm. impact par. d_{z} / #sigma_{z}");
  }
  for(Int_t i=0; i<kNHist; i++) fhQA[i][1]->SetLineColor(color);
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::FillHistograms(TObject* obj, Bool_t f)
{
  //
  // fill the QA histograms
  //

  if (!obj) return;

  Bool_t isESDTrack = strcmp(obj->ClassName(),"AliESDtrack") == 0 ? kTRUE : kFALSE ;
  Bool_t isAODTrack = strcmp(obj->ClassName(),"AliAODTrack") == 0 ? kTRUE : kFALSE ;

  AliESDtrack * esdTrack = 0x0 ;
  AliAODTrack * aodTrack = 0x0 ; 
  if (isESDTrack) esdTrack = dynamic_cast<AliESDtrack*>(obj);
  if (isAODTrack) aodTrack = dynamic_cast<AliAODTrack*>(obj);

  // f = 0: fill histograms before cuts
  // f = 1: fill histograms after cuts

  if (esdTrack) {
    Float_t b[2];
    Float_t bRes[2];
    Float_t bCov[3];
    esdTrack->GetImpactParameters(b,bCov);
    if (bCov[0]<=0 || bCov[2]<=0) {
      AliDebug(1, "Estimated b resolution lower or equal zero!");
      bCov[0]=0; bCov[2]=0;
    }
    bRes[0] = TMath::Sqrt(bCov[0]);
    bRes[1] = TMath::Sqrt(bCov[2]);
  
    fhQA[kDcaZ][f]->Fill(b[1]);
    fhQA[kDcaXY][f]->Fill(b[0]);
    fhDcaXYvsDcaZ[f]->Fill(b[1],b[0]);
    
    if (bRes[0]!=0 && bRes[1]!=0) {
      fhQA[kDcaZnorm][f]->Fill(b[1]/bRes[1]);
      fhQA[kDcaXYnorm][f]->Fill(b[0]/bRes[0]);
      fhDcaXYvsDcaZnorm[f]->Fill(b[1]/bRes[1], b[0]/bRes[0]);
    }

    fhQA[kCutNSigmaToVertex][f]->Fill(fNSigmaToVertex);
    if (fNSigmaToVertex<0 && fRequireSigmaToVertex) fhQA[kCutRequireSigmaToVertex][f]->Fill(0.);
    if (!(fNSigmaToVertex<0 && fRequireSigmaToVertex)) fhQA[kCutRequireSigmaToVertex][f]->Fill(1.);
    
    if (!fAcceptKinkDaughters && esdTrack->GetKinkIndex(0)>0) fhQA[kCutAcceptKinkDaughters][f]->Fill(0.);
    if (!(!fAcceptKinkDaughters && esdTrack->GetKinkIndex(0)>0)) fhQA[kCutAcceptKinkDaughters][f]->Fill(0.);
  }
    
  // fill cut statistics and cut correlation histograms with information from the bitmap
  if (f) return;

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
void AliCFTrackIsPrimaryCuts::SaveHistograms(const Char_t* dir) {
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

    fhDcaXYvsDcaZ[j]    ->Write();
    fhDcaXYvsDcaZnorm[j]->Write();

    for(Int_t i=0; i<kNHist; i++) fhQA[i][j]->Write();

    gDirectory->cd("../");
  }
  gDirectory->cd("../");
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::DrawHistograms()
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

  TCanvas* canvas1 = new TCanvas("Track_QA_Primary_1", "Track QA Primary 1", 800, 500);
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

  TCanvas* canvas2 = new TCanvas("Track_QA_Primary_2", "Track QA Primary 2", 800, 800);
  canvas2->Divide(2, 2);

  canvas2->cd(1);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kDcaXY][0]->SetStats(kFALSE);
  fhQA[kDcaXY][0]->Draw();
  fhQA[kDcaXY][1]->Draw("same");

  canvas2->cd(2);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kDcaZ][0]->SetStats(kFALSE);
  fhQA[kDcaZ][0]->Draw();
  fhQA[kDcaZ][1]->Draw("same");

  canvas2->cd(3);
//   fhDXYvsDZ[0]->SetStats(kFALSE);
  gPad->SetLogz();
  gPad->SetLeftMargin(bottom);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(bottom);
  gPad->SetRightMargin(0.2);
  fhDcaXYvsDcaZ[0]->Draw("COLZ");

  canvas2->cd(4);
//   fhDXYvsDZ[1]->SetStats(kFALSE);
  gPad->SetLogz();
  gPad->SetLeftMargin(bottom);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(bottom);
  gPad->SetRightMargin(0.2);
  fhDcaXYvsDcaZ[1]->Draw("COLZ");

  canvas2->SaveAs(Form("%s.eps", canvas2->GetName()));
  canvas2->SaveAs(Form("%s.ps", canvas2->GetName()));

  // -----

  TCanvas* canvas3 = new TCanvas("Track_QA_Primary_3", "Track QA Primary 3", 800, 800);
  canvas3->Divide(2, 2);

  canvas3->cd(1);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kDcaXYnorm][0]->SetStats(kFALSE);
  fhQA[kDcaXYnorm][0]->Draw();
  fhQA[kDcaXYnorm][1]->Draw("same");

  canvas3->cd(2);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kDcaZnorm][0]->SetStats(kFALSE);
  fhQA[kDcaZnorm][0]->Draw();
  fhQA[kDcaZnorm][1]->Draw("same");

  canvas3->cd(3);
//   fhDXYvsDZ[0]->SetStats(kFALSE);
  gPad->SetLogz();
  gPad->SetLeftMargin(bottom);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(bottom);
  gPad->SetRightMargin(0.2);
  fhDcaXYvsDcaZnorm[0]->Draw("COLZ");

  canvas3->cd(4);
//   fhDXYvsDZ[1]->SetStats(kFALSE);
  gPad->SetLogz();
  gPad->SetLeftMargin(bottom);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(bottom);
  gPad->SetRightMargin(0.2);
  fhDcaXYvsDcaZnorm[1]->Draw("COLZ");

  canvas3->SaveAs(Form("%s.eps", canvas3->GetName()));
  canvas3->SaveAs(Form("%s.ps", canvas3->GetName()));

  // -----

  TCanvas* canvas4 = new TCanvas("Track_QA_Primary_4", "Track QA Primary 4", 1200, 500);
  canvas4->Divide(3, 1);

  canvas4->cd(1);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutNSigmaToVertex][0]->SetStats(kFALSE);
  fhQA[kCutNSigmaToVertex][0]->Draw();
  fhQA[kCutNSigmaToVertex][1]->Draw("same");

  canvas4->cd(2);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutRequireSigmaToVertex][0]->SetStats(kFALSE);
  fhQA[kCutRequireSigmaToVertex][0]->Draw();
  fhQA[kCutRequireSigmaToVertex][1]->Draw("same");

  canvas4->cd(3);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutAcceptKinkDaughters][0]->SetStats(kFALSE);
  fhQA[kCutAcceptKinkDaughters][0]->Draw();
  fhQA[kCutAcceptKinkDaughters][1]->Draw("same");

  canvas4->SaveAs(Form("%s.eps", canvas4->GetName()));
  canvas4->SaveAs(Form("%s.ps", canvas4->GetName()));
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::AddQAHistograms(TList *qaList) {
  //
  // saves the histograms in a TList
  //
  DefineHistograms();

  qaList->Add(fhCutStatistics);
  qaList->Add(fhCutCorrelation);

  for (Int_t j=0; j<kNStepQA; j++) {
    qaList->Add(fhDcaXYvsDcaZ[j]);
    qaList->Add(fhDcaXYvsDcaZnorm[j]);
    for(Int_t i=0; i<kNHist; i++)
	qaList->Add(fhQA[i][j]);
  }
}
