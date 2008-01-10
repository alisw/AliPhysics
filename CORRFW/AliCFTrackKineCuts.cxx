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

// The class AliCFTrackKineCuts is designed to select both generated 
// and reconstructed tracks of a given range in momentum space,
// electric charge and azimuthal emission angle phi
// and to provide corresponding QA histograms.
// This class inherits from the Analysis' Framework abstract base class
// AliAnalysisCuts and is a part of the Correction Framework.
// This class acts on single, generated and reconstructed tracks, it is 
// applicable on ESD and AOD data.
// It mainly consists of a IsSelected function that returns a boolean.
// This function checks whether the considered track passes a set of cuts:
// - total momentum
// - pt
// - px
// - py
// - pz
// - eta
// - rapidity
// - phi
// - charge
// - is charged
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
#include <TBits.h>
#include <TH2.h>

#include <AliVParticle.h>
#include <AliLog.h>
#include "AliCFTrackKineCuts.h"

ClassImp(AliCFTrackKineCuts)

//__________________________________________________________________________________
AliCFTrackKineCuts::AliCFTrackKineCuts() :
  AliCFCutBase(),
  fMomentumMin(0),
  fMomentumMax(0),
  fPtMin(0),
  fPtMax(0),
  fPxMin(0),
  fPxMax(0),
  fPyMin(0),
  fPyMax(0),
  fPzMin(0),
  fPzMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fRapidityMin(0),
  fRapidityMax(0),
  fPhiMin(0),
  fPhiMax(0),
  fCharge(0),
  fRequireIsCharged(0),
  fhCutStatistics(0),
  fhCutCorrelation(0),
  fBitmap(0x0),
  fhNBinsMomentum(0),
  fhNBinsPt(0),
  fhNBinsPx(0),
  fhNBinsPy(0),
  fhNBinsPz(0),
  fhNBinsEta(0),
  fhNBinsRapidity(0),
  fhNBinsPhi(0),
  fhNBinsCharge(0),
  fhBinLimMomentum(0x0),
  fhBinLimPt(0x0),
  fhBinLimPx(0x0),
  fhBinLimPy(0x0),
  fhBinLimPz(0x0),
  fhBinLimEta(0x0),
  fhBinLimRapidity(0x0),
  fhBinLimPhi(0x0),
  fhBinLimCharge(0x0)

{
  //
  // Default constructor
  //
  fBitmap=new TBits(0);
  Initialise();
}
//__________________________________________________________________________________
AliCFTrackKineCuts::AliCFTrackKineCuts(Char_t* name, Char_t* title) :
  AliCFCutBase(name,title),
  fMomentumMin(0),
  fMomentumMax(0),
  fPtMin(0),
  fPtMax(0),
  fPxMin(0),
  fPxMax(0),
  fPyMin(0),
  fPyMax(0),
  fPzMin(0),
  fPzMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fRapidityMin(0),
  fRapidityMax(0),
  fPhiMin(0),
  fPhiMax(0),
  fCharge(0),
  fRequireIsCharged(0),
  fhCutStatistics(0),
  fhCutCorrelation(0),
  fBitmap(0x0),
  fhNBinsMomentum(0),
  fhNBinsPt(0),
  fhNBinsPx(0),
  fhNBinsPy(0),
  fhNBinsPz(0),
  fhNBinsEta(0),
  fhNBinsRapidity(0),
  fhNBinsPhi(0),
  fhNBinsCharge(0),
  fhBinLimMomentum(0x0),
  fhBinLimPt(0x0),
  fhBinLimPx(0x0),
  fhBinLimPy(0x0),
  fhBinLimPz(0x0),
  fhBinLimEta(0x0),
  fhBinLimRapidity(0x0),
  fhBinLimPhi(0x0),
  fhBinLimCharge(0x0)

{
  //
  // Constructor
  //
  fBitmap=new TBits(0);
  Initialise();
}
//__________________________________________________________________________________
AliCFTrackKineCuts::AliCFTrackKineCuts(const AliCFTrackKineCuts& c) :
  AliCFCutBase(c),
  fMomentumMin(c.fMomentumMin),
  fMomentumMax(c.fMomentumMax),
  fPtMin(c.fPtMin),
  fPtMax(c.fPtMax),
  fPxMin(c.fPxMin),
  fPxMax(c.fPxMax),
  fPyMin(c.fPyMin),
  fPyMax(c.fPyMax),
  fPzMin(c.fPzMin),
  fPzMax(c.fPzMax),
  fEtaMin(c.fEtaMin),
  fEtaMax(c.fEtaMax),
  fRapidityMin(c.fRapidityMin),
  fRapidityMax(c.fRapidityMax),
  fPhiMin(c.fPhiMin),
  fPhiMax(c.fPhiMax),
  fCharge(c.fCharge),
  fRequireIsCharged(c.fRequireIsCharged),
  fhCutStatistics(c.fhCutStatistics),
  fhCutCorrelation(c.fhCutCorrelation),
  fBitmap(c.fBitmap),
  fhNBinsMomentum(c.fhNBinsMomentum),
  fhNBinsPt(c.fhNBinsPt),
  fhNBinsPx(c.fhNBinsPx),
  fhNBinsPy(c.fhNBinsPy),
  fhNBinsPz(c.fhNBinsPz),
  fhNBinsEta(c.fhNBinsEta),
  fhNBinsRapidity(c.fhNBinsRapidity),
  fhNBinsPhi(c.fhNBinsPhi),
  fhNBinsCharge(c.fhNBinsCharge),
  fhBinLimMomentum(c.fhBinLimMomentum),
  fhBinLimPt(c.fhBinLimPt),
  fhBinLimPx(c.fhBinLimPx),
  fhBinLimPy(c.fhBinLimPy),
  fhBinLimPz(c.fhBinLimPz),
  fhBinLimEta(c.fhBinLimEta),
  fhBinLimRapidity(c.fhBinLimRapidity),
  fhBinLimPhi(c.fhBinLimPhi),
  fhBinLimCharge(c.fhBinLimCharge)

{
  //
  // copy constructor
  //
  ((AliCFTrackKineCuts &) c).Copy(*this);
}
//__________________________________________________________________________________
AliCFTrackKineCuts& AliCFTrackKineCuts::operator=(const AliCFTrackKineCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fMomentumMin = c.fMomentumMin ;
    fMomentumMax = c.fMomentumMax ;
    fPtMin = c.fPtMin ;
    fPtMax = c.fPtMax ;
    fPxMin = c.fPxMin ;
    fPxMax = c.fPxMax ;
    fPyMin = c.fPyMin ;
    fPyMax = c.fPyMax ;
    fPzMin = c.fPzMin ;
    fPzMax = c.fPzMax ;
    fEtaMin = c.fEtaMin ;
    fEtaMax = c.fEtaMax ;
    fRapidityMin = c.fRapidityMin ;
    fRapidityMax = c.fRapidityMax ;
    fPhiMin = c.fPhiMin ;
    fPhiMax = c.fPhiMax ;
    fCharge = c.fCharge ;
    fRequireIsCharged = c.fRequireIsCharged ;
    fhCutStatistics = c.fhCutStatistics ;
    fhCutCorrelation = c.fhCutCorrelation ;
    fBitmap = c.fBitmap;
    fhNBinsMomentum = c.fhNBinsMomentum;
    fhNBinsPt = c.fhNBinsPt;
    fhNBinsPx = c.fhNBinsPx;
    fhNBinsPy = c.fhNBinsPy;
    fhNBinsPz = c.fhNBinsPz;
    fhNBinsEta = c.fhNBinsEta;
    fhNBinsRapidity = c.fhNBinsRapidity;
    fhNBinsPhi = c.fhNBinsPhi;
    fhNBinsCharge = c.fhNBinsCharge;
    fhBinLimMomentum = c.fhBinLimMomentum;
    fhBinLimPt = c.fhBinLimPt;
    fhBinLimPx = c.fhBinLimPx;
    fhBinLimPy = c.fhBinLimPy;
    fhBinLimPz = c.fhBinLimPz;
    fhBinLimEta = c.fhBinLimEta;
    fhBinLimRapidity = c.fhBinLimRapidity;
    fhBinLimPhi = c.fhBinLimPhi;
    fhBinLimCharge = c.fhBinLimCharge;
    
    for (Int_t i=0; i<c.kNHist; i++){
      for (Int_t j=0; j<c.kNStepQA; j++){
	if(c.fhQA[i][j]) fhQA[i][j] = (TH1F*)c.fhQA[i][j]->Clone();
      }
    }

    ((AliCFTrackKineCuts &) c).Copy(*this);
  }
  return *this ;
}
//__________________________________________________________________________________
AliCFTrackKineCuts::~AliCFTrackKineCuts()
{
  //
  // destructor
  //
  if (fhCutStatistics)			delete fhCutStatistics;
  if (fhCutCorrelation)			delete fhCutCorrelation;

  for (Int_t i=0; i<kNHist; i++){
    for (Int_t j=0; j<kNStepQA; j++){
      if(fhQA[i][j]) delete fhQA[i][j];
    }
  }

  if(fBitmap)	delete   fBitmap;

  if(fhBinLimMomentum) delete fhBinLimMomentum;
  if(fhBinLimPt) delete fhBinLimPt;
  if(fhBinLimPx) delete fhBinLimPx;
  if(fhBinLimPy) delete fhBinLimPy;
  if(fhBinLimPz) delete fhBinLimPz;
  if(fhBinLimEta) delete fhBinLimEta;
  if(fhBinLimRapidity) delete fhBinLimRapidity;
  if(fhBinLimPhi) delete fhBinLimPhi;
  if(fhBinLimCharge) delete fhBinLimCharge;
}
//__________________________________________________________________________________
void AliCFTrackKineCuts::Initialise()
{
  //
  // sets everything to zero
  //
  fMomentumMin = 0;
  fMomentumMax = 0;
  fPtMin = 0;
  fPtMax = 0;
  fPxMin = 0;
  fPxMax = 0;
  fPyMin = 0;
  fPyMax = 0;
  fPzMin = 0;
  fPzMax = 0;
  fEtaMin = 0;
  fEtaMax = 0;
  fRapidityMin = 0;
  fRapidityMax = 0;
  fPhiMin = 0;
  fPhiMax = 0;
  fCharge = 0;
  fRequireIsCharged = 0;

  SetMomentumRange();
  SetPtRange();
  SetPxRange();
  SetPyRange();
  SetPzRange();
  SetEtaRange();
  SetRapidityRange();
  SetPhiRange();
  SetChargeRec();
  SetChargeMC();
  SetRequireIsCharged();

  for (Int_t i=0; i<kNHist; i++){
    for (Int_t j=0; j<kNStepQA; j++)  {
      fhQA[i][j] = 0x0;
    }
  }

  fhCutStatistics = 0;
  fhCutCorrelation = 0;
    
    //set default bining for QA histograms
  SetHistogramBins(kCutP,200,0.,20.);
  SetHistogramBins(kCutPt,200,0.,20.);
  SetHistogramBins(kCutPx,400,-20.,20.);
  SetHistogramBins(kCutPy,400,-20.,20.);
  SetHistogramBins(kCutPz,400,-20.,20.);
  SetHistogramBins(kCutRapidity,200,-10.,10.);
  SetHistogramBins(kCutEta,200,-10.,10.);
  SetHistogramBins(kCutPhi,38,-0.6,7.);
  SetHistogramBins(kCutCharge,21,-10.5,10.5);

}
//__________________________________________________________________________________
void AliCFTrackKineCuts::Copy(TObject &c) const
{
  //
  // Copy function
  //
  AliCFTrackKineCuts& target = (AliCFTrackKineCuts &) c;

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
void AliCFTrackKineCuts::GetBitMap(TObject* obj, TBits *bitmap)  {
  //
  // retrieve the pointer to the bitmap
  //
  bitmap = SelectionBitMap(obj);
}
//__________________________________________________________________________________
TBits* AliCFTrackKineCuts::SelectionBitMap(TObject* obj) {
  //
  // test if the track passes the single cuts
  // and store the information in a bitmap
  //

  // bitmap stores the decision of each single cut
  for(Int_t i=0; i<kNCuts; i++)fBitmap->SetBitNumber(i,kFALSE);
  // cast TObject into VParticle
  AliVParticle* particle = dynamic_cast<AliVParticle *>(obj);
  if ( !particle ) return fBitmap ;


  for(Int_t i=0; i<kNCuts; i++)fBitmap->SetBitNumber(i,kTRUE);

  Int_t iCutBit = 0;
  if((particle->P() < fMomentumMin) || (particle->P() > fMomentumMax))
	fBitmap->SetBitNumber(iCutBit,kFALSE);
  iCutBit++;
  if ((particle->Pt() < fPtMin) || (particle->Pt() > fPtMax))
	fBitmap->SetBitNumber(iCutBit,kFALSE);
  iCutBit++;
  if ((particle->Px() < fPxMin) || (particle->Px() > fPxMax))
	fBitmap->SetBitNumber(iCutBit,kFALSE);
  iCutBit++;
  if ((particle->Py() < fPyMin) || (particle->Py() > fPyMax))
	fBitmap->SetBitNumber(iCutBit,kFALSE);
  iCutBit++;
  if ((particle->Pz() < fPzMin) || (particle->Pz() > fPzMax))
	fBitmap->SetBitNumber(iCutBit,kFALSE);
  iCutBit++;
  if ((particle->Eta() < fEtaMin) || (particle->Eta() > fEtaMax))
	fBitmap->SetBitNumber(iCutBit,kFALSE);
  iCutBit++;
  if ((particle->Y() < fRapidityMin) || (particle->Y() > fRapidityMax))
	fBitmap->SetBitNumber(iCutBit,kFALSE);
  iCutBit++;
  if ((particle->Phi() < fPhiMin) || (particle->Phi() > fPhiMax))
	fBitmap->SetBitNumber(iCutBit,kFALSE);
  iCutBit++;
  if (fCharge < 10 && particle->Charge() != fCharge)
	fBitmap->SetBitNumber(iCutBit,kFALSE);
  iCutBit++;
  if (fRequireIsCharged && particle->Charge()==0)
	fBitmap->SetBitNumber(iCutBit,kFALSE);

  return fBitmap;
}
//__________________________________________________________________________________
Bool_t AliCFTrackKineCuts::IsSelected(TObject* obj) {
  //
  // loops over decisions of single cuts and returns if the track is accepted
  //
  TBits* bitmap = SelectionBitMap(obj);

  Bool_t isSelected = kTRUE;

  for (UInt_t icut=0; icut<bitmap->GetNbits();icut++)
	if(!bitmap->TestBitNumber(icut)) isSelected = kFALSE;

  return isSelected;
}
//__________________________________________________________________________________
void AliCFTrackKineCuts::Init() {
  //
  // initialises all histograms and the TList which holds the histograms
  //
  if(fIsQAOn)
    DefineHistograms();
}
//__________________________________________________________________________________
void AliCFTrackKineCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins)
{
  //
  // variable bin size
  //
  if(!fIsQAOn) return;

  switch(index){
  case kCutP:
    fhNBinsMomentum=nbins;
    fhBinLimMomentum=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimMomentum[i]=bins[i];
    break;
    
  case kCutPt:
    fhNBinsPt=nbins;
    fhBinLimPt=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimPt[i]=bins[i];
    break;
    
  case kCutPx:
    fhNBinsPx=nbins;
    fhBinLimPx=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimPx[i]=bins[i];
    break;

  case kCutPy:
    fhNBinsPy=nbins;
    fhBinLimPy=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimPy[i]=bins[i];
    break;
        
  case kCutPz:
    fhNBinsPz=nbins;
    fhBinLimPz=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimPz[i]=bins[i];
    break;
    
  case kCutRapidity:
    fhNBinsRapidity=nbins;
    fhBinLimRapidity=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimRapidity[i]=bins[i];
    break;
    
  case kCutEta:
    fhNBinsEta=nbins;
    fhBinLimEta=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimEta[i]=bins[i];
    break;
    
  case kCutPhi:
    fhNBinsPhi=nbins;
    fhBinLimPhi=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimPhi[i]=bins[i];
    break;
    
  case kCutCharge:
    fhNBinsCharge=nbins;
    fhBinLimCharge=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimCharge[i]=bins[i];
    break;
  }
}
//__________________________________________________________________________________
void AliCFTrackKineCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax)
{
  //
  // fixed bin size
  //
  if(!fIsQAOn) return;

  switch(index){
  case kCutP:
    fhNBinsMomentum=nbins;
    fhBinLimMomentum=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimMomentum[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
    
  case kCutPt:
    fhNBinsPt=nbins;
    fhBinLimPt=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimPt[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
    
  case kCutPx:
    fhNBinsPx=nbins;
    fhBinLimPx=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimPx[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutPy:
    fhNBinsPy=nbins;
    fhBinLimPy=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimPy[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
        
  case kCutPz:
    fhNBinsPz=nbins;
    fhBinLimPz=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimPz[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
    
  case kCutRapidity:
    fhNBinsRapidity=nbins;
    fhBinLimRapidity=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimRapidity[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
    
  case kCutEta:
    fhNBinsEta=nbins;
    fhBinLimEta=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimEta[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
    
  case kCutPhi:
    fhNBinsPhi=nbins;
    fhBinLimPhi=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimPhi[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
    
  case kCutCharge:
    fhNBinsCharge=nbins;
    fhBinLimCharge=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimCharge[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
  }
}
//__________________________________________________________________________________
 void AliCFTrackKineCuts::DefineHistograms() {
  //
  // histograms for cut variables, cut statistics and cut correlations
  //

  Int_t color = 2;

  // book cut statistics and cut correlation histograms
  fhCutStatistics = new TH1F(Form("%s_cut_statistics",GetName()), Form("%s cut statistics",GetName()), kNCuts,0.5,kNCuts+0.5);
  fhCutStatistics->SetLineWidth(2);
  fhCutStatistics->GetXaxis()->SetBinLabel(1,"p");
  fhCutStatistics->GetXaxis()->SetBinLabel(2,"pt");
  fhCutStatistics->GetXaxis()->SetBinLabel(3,"px");
  fhCutStatistics->GetXaxis()->SetBinLabel(4,"py");
  fhCutStatistics->GetXaxis()->SetBinLabel(5,"pz");
  fhCutStatistics->GetXaxis()->SetBinLabel(6,"eta");
  fhCutStatistics->GetXaxis()->SetBinLabel(7,"y");
  fhCutStatistics->GetXaxis()->SetBinLabel(8,"phi");
  fhCutStatistics->GetXaxis()->SetBinLabel(9,"charge");
  fhCutStatistics->GetXaxis()->SetBinLabel(10,"is charged");

  fhCutCorrelation = new TH2F(Form("%s_cut_correlation",GetName()), Form("%s cut  correlation",GetName()), kNCuts,0.5,kNCuts+0.5,kNCuts,0.5,kNCuts+0.5);
  fhCutCorrelation->SetLineWidth(2);
  fhCutCorrelation->GetXaxis()->SetBinLabel(1,"p");
  fhCutCorrelation->GetXaxis()->SetBinLabel(2,"pt");
  fhCutCorrelation->GetXaxis()->SetBinLabel(3,"px");
  fhCutCorrelation->GetXaxis()->SetBinLabel(4,"py");
  fhCutCorrelation->GetXaxis()->SetBinLabel(5,"pz");
  fhCutCorrelation->GetXaxis()->SetBinLabel(6,"eta");
  fhCutCorrelation->GetXaxis()->SetBinLabel(7,"y");
  fhCutCorrelation->GetXaxis()->SetBinLabel(8,"phi");
  fhCutCorrelation->GetXaxis()->SetBinLabel(9,"charge");
  fhCutCorrelation->GetXaxis()->SetBinLabel(10,"is charged");

  fhCutCorrelation->GetYaxis()->SetBinLabel(1,"p");
  fhCutCorrelation->GetYaxis()->SetBinLabel(2,"pt");
  fhCutCorrelation->GetYaxis()->SetBinLabel(3,"px");
  fhCutCorrelation->GetYaxis()->SetBinLabel(4,"py");
  fhCutCorrelation->GetYaxis()->SetBinLabel(5,"pz");
  fhCutCorrelation->GetYaxis()->SetBinLabel(6,"eta");
  fhCutCorrelation->GetYaxis()->SetBinLabel(7,"y");
  fhCutCorrelation->GetYaxis()->SetBinLabel(8,"phi");
  fhCutCorrelation->GetYaxis()->SetBinLabel(9,"charge");
  fhCutCorrelation->GetYaxis()->SetBinLabel(10,"is charged");


  // book QA histograms
  Char_t str[256];
  for (Int_t i=0; i<kNStepQA; i++) {
    if (i==0) sprintf(str," ");
    else sprintf(str,"_cut");
  
    fhQA[kCutP][i]	= new  TH1F(Form("%s_momentum%s",GetName(),str),	"",fhNBinsMomentum,fhBinLimMomentum);
    fhQA[kCutPt][i]	= new  TH1F(Form("%s_transverse_momentum%s",GetName(),str),"",fhNBinsPt,fhBinLimPt);
    fhQA[kCutPx][i]	= new  TH1F(Form("%s_px%s",GetName(),str),		"",fhNBinsPx,fhBinLimPx);
    fhQA[kCutPy][i]	= new  TH1F(Form("%s_py%s",GetName(),str),		"",fhNBinsPy,fhBinLimPy);
    fhQA[kCutPz][i]	= new  TH1F(Form("%s_pz%s",GetName(),str),		"",fhNBinsPz,fhBinLimPz);
    fhQA[kCutRapidity][i]=new  TH1F(Form("%s_rapidity%s",GetName(),str),	"",fhNBinsRapidity,fhBinLimRapidity);
    fhQA[kCutEta][i]	= new  TH1F(Form("%s_eta%s",GetName(),str),		"",fhNBinsEta,fhBinLimEta);
    fhQA[kCutPhi][i]	= new  TH1F(Form("%s_phi%s",GetName(),str),		"",fhNBinsPhi,fhBinLimPhi);
    fhQA[kCutCharge][i]	= new  TH1F(Form("%s_charge%s",GetName(),str),		"",fhNBinsCharge,fhBinLimCharge);

    fhQA[kCutP][i]	->SetXTitle("momentum p (GeV/c)");
    fhQA[kCutPt][i]	->SetXTitle("p_{T} (GeV/c)");
    fhQA[kCutPx][i]	->SetXTitle("p_{x} (GeV/c)");
    fhQA[kCutPy][i]	->SetXTitle("p_{y} (GeV/c)");
    fhQA[kCutPz][i]	->SetXTitle("p_{z} (GeV/c)");
    fhQA[kCutRapidity][i]->SetXTitle("rapidity y");
    fhQA[kCutEta][i]	->SetXTitle("pseudo rapidity #eta");
    fhQA[kCutPhi][i]	->SetXTitle("azimuth #phi (rad)");
    fhQA[kCutCharge][i]	->SetXTitle("charge");
  }

  for(Int_t i=0; i<kNHist; i++) fhQA[i][1]->SetLineColor(color);

}
//__________________________________________________________________________________
void AliCFTrackKineCuts::FillHistograms(TObject* obj, Bool_t b)
{
  //
  // fill the QA histograms
  //
  if(!fIsQAOn) return;

  // cast TObject into VParticle
  AliVParticle* particle = dynamic_cast<AliVParticle *>(obj);
  if ( !particle ) return;

  // index = 0: fill histograms before cuts
  // index = 1: fill histograms after cuts
  Int_t index = -1;
  index = ((b) ? 1 : 0);

  fhQA[kCutP][index]->Fill(particle->P());
  fhQA[kCutPt][index]->Fill(particle->Pt());
  fhQA[kCutPx][index]->Fill(particle->Px());
  fhQA[kCutPy][index]->Fill(particle->Py());
  fhQA[kCutPz][index]->Fill(particle->Pz());
  fhQA[kCutRapidity][index]->Fill(particle->Y());
  fhQA[kCutEta][index]->Fill(particle->Eta());
  fhQA[kCutPhi][index]->Fill(particle->Phi());
  fhQA[kCutCharge][index]->Fill((float)particle->Charge());

  // fill cut statistics and cut correlation histograms with information from the bitmap
  if (b) return;

  if (!obj) return;
  TBits* bitmap = SelectionBitMap(obj);

  // Number of single cuts in this class
  UInt_t ncuts = bitmap->GetNbits();
  for(UInt_t bit=0; bit<ncuts;bit++) {
    if (!bitmap->TestBitNumber(bit)) {
	fhCutStatistics->Fill(bit+1);
	for (UInt_t bit2=bit; bit2<ncuts;bit2++) {
	  if (!bitmap->TestBitNumber(bit2)) 
	    fhCutCorrelation->Fill(bit+1,bit2+1);
	}
    }
  }
}
//__________________________________________________________________________________
void AliCFTrackKineCuts::SaveHistograms(const Char_t* dir) {
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
void AliCFTrackKineCuts::DrawHistograms(Bool_t drawLogScale)
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

  TCanvas* canvas1 = new TCanvas("Track_QA_Kinematics_1", "Track QA Kinematics 1", 800, 500);
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

  TCanvas* canvas2 = new TCanvas("Track_QA_Kinematics_2", "Track QA Kinematics 2", 1600, 800);
  canvas2->Divide(4, 2);

  canvas2->cd(1);
  fhQA[kCutP][0]->SetStats(kFALSE);
  if(drawLogScale) gPad->SetLogy();
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutP][0]->Draw();
  fhQA[kCutP][1]->Draw("same");

  canvas2->cd(2);
  fhQA[kCutPt][0]->SetStats(kFALSE);
  if(drawLogScale) gPad->SetLogy();
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutPt][0]->Draw();
  fhQA[kCutPt][1]->Draw("same");

  canvas2->cd(3);
  fhQA[kCutRapidity][0]->SetStats(kFALSE);
  if(drawLogScale) gPad->SetLogy();
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutRapidity][0]->Draw();
  fhQA[kCutRapidity][1]->Draw("same");

  canvas2->cd(4);
  fhQA[kCutEta][0]->SetStats(kFALSE);
  if(drawLogScale) gPad->SetLogy();
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutEta][0]->Draw();
  fhQA[kCutEta][1]->Draw("same");

  canvas2->cd(5);
  fhQA[kCutPx][0]->SetStats(kFALSE);
  if(drawLogScale) gPad->SetLogy();
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutPx][0]->Draw();
  fhQA[kCutPx][1]->Draw("same");

  canvas2->cd(6);
  fhQA[kCutPy][0]->SetStats(kFALSE);
  if(drawLogScale) gPad->SetLogy();
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutPy][0]->Draw();
  fhQA[kCutPy][1]->Draw("same");

  canvas2->cd(7);
  fhQA[kCutPz][0]->SetStats(kFALSE);
  if(drawLogScale) gPad->SetLogy();
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutPz][0]->Draw();
  fhQA[kCutPz][1]->Draw("same");

  canvas2->SaveAs(Form("%s.eps", canvas2->GetName()));
  canvas2->SaveAs(Form("%s.ps", canvas2->GetName()));

  // -----

  TCanvas* canvas3 = new TCanvas("Track_QA_Kinematics_3", "Track QA Kinematics 3", 800, 400);
  canvas3->Divide(2, 1);

  canvas3->cd(1);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutPhi][0]->SetStats(kFALSE);
  fhQA[kCutPhi][0]->Draw();
  fhQA[kCutPhi][1]->Draw("same");

  canvas3->cd(2);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutCharge][0]->SetStats(kFALSE);
  fhQA[kCutCharge][0]->Draw();
  fhQA[kCutCharge][1]->Draw("same");

  canvas3->SaveAs(Form("%s.eps", canvas3->GetName()));
  canvas3->SaveAs(Form("%s.ps", canvas3->GetName()));
}
//__________________________________________________________________________________
void AliCFTrackKineCuts::AddQAHistograms(TList *qaList) const {
  //
  // saves the histograms in a TList
  //
  if(!fIsQAOn) return;

  qaList->Add(fhCutStatistics);
  qaList->Add(fhCutCorrelation);

  for (Int_t j=0; j<kNStepQA; j++) {
    for(Int_t i=0; i<kNHist; i++)
	qaList->Add(fhQA[i][j]);
  }
}
