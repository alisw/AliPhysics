#include "AliHFCutVarFDsubMassFitter.h"

#include "TList.h"
#include "TH1F.h"
#include "AliHFMassFitter.h"

#include "AliHFCutVarFDsubAxis.h"
#include "AliHFCutVarFDsubCut.h"
#include "AliHFCutVarFDsubCutSet.h"

/// \cond CLASSIMP
ClassImp(AliHFCutVarFDsubMassFitter);
/// \endcond


AliHFCutVarFDsubMassFitter::AliHFCutVarFDsubMassFitter()
  : TObject()
  , fTHn(0x0)
  , fAxes(0x0)
  , fCutSet(0x0)
  , fMassHist(0x0)
  , fFitter(0x0)
  , fSig(-1.)
  , fSigErr(-1.)
  , fBkg(-1.)
  , fBkgErr(-1.)
  , fSigma(-1.)
  , fSigmaErr(-1.)
  , fMean(-1.)
  , fMeanErr(-1.)
  , fChiSquare(-1.)
  , fRedChiSquare(-1.)
{
  // Default constructor
}


AliHFCutVarFDsubMassFitter::AliHFCutVarFDsubMassFitter(THnSparseF* thn, TList* axes, AliHFCutVarFDsubCutSet* cutSet)
  : TObject()
  , fTHn(thn)
  , fAxes(axes)
  , fCutSet(cutSet)
  , fMassHist(0x0)
  , fFitter(0x0)
  , fSig(-1.)
  , fSigErr(-1.)
  , fBkg(-1.)
  , fBkgErr(-1.)
  , fSigma(-1.)
  , fSigmaErr(-1.)
  , fMean(-1.)
  , fMeanErr(-1.)
  , fChiSquare(-1.)
  , fRedChiSquare(-1.)
{
  // Constructor
  ObtainMassHistogram();
  FitMassHistogram();
}


AliHFCutVarFDsubMassFitter::~AliHFCutVarFDsubMassFitter() {
  delete fFitter;
  fFitter = 0x0;
  delete fMassHist;
  fMassHist = 0x0;
}


void AliHFCutVarFDsubMassFitter::ObtainMassHistogram() {
  Int_t SymmCutAxisCounter=0;
  for (Int_t iCut=0; iCut<fCutSet->GetEntries(); ++iCut) {
    AliHFCutVarFDsubCut* cut = fCutSet->GetCut(iCut);
    AliHFCutVarFDsubAxis* axis = (AliHFCutVarFDsubAxis*)fAxes->At(cut->fAxisId);
    if(axis->IsCutSymmetric() && cut->fHigh != -cut->fLow)
      SymmCutAxisCounter++;
  }
  if(SymmCutAxisCounter==0) { // if none of the axis has a two-region cut symmetric with respect to zero just project
    fMassHist = (TH1F*)ProjectDataSparse(kFALSE);
  }
  else { // else get the histograms corresponding to the two symmetric regions and sum them
    TH1F* hRegion1 = (TH1F*)ProjectDataSparse(kFALSE);
    TH1F* hRegion2 = (TH1F*)ProjectDataSparse(kTRUE);
    fMassHist=(TH1F*)hRegion1->Clone();
    fMassHist->Add(hRegion1,hRegion2,1.,1.);
    delete hRegion1;
    hRegion1=0x0;
    delete hRegion2;
    hRegion2=0x0;
  }
  fMassHist->Rebin(fCutSet->fFitRebin);
}


TH1F* AliHFCutVarFDsubMassFitter::ProjectDataSparse(Bool_t reflectedaxes) {
  for (Int_t iCut=0; iCut<fCutSet->GetEntries(); ++iCut) {
    AliHFCutVarFDsubCut* cut = fCutSet->GetCut(iCut);
    AliHFCutVarFDsubAxis* axis = (AliHFCutVarFDsubAxis*)fAxes->At(cut->fAxisId);
    UInt_t axisNo = axis->GetAxisNo(AliHFCutVarFDsubAxis::kData);
    Bool_t isCutSymm = axis->IsCutSymmetric();
    if (axisNo<(UInt_t)-1) {
      TAxis* ax = fTHn->GetAxis(axisNo);
      Int_t binMin;
      Int_t binMax;
      if(reflectedaxes && isCutSymm) {
        binMin = ax->FindBin(-cut->fHigh*1.0001);
        binMax = ax->FindBin(-cut->fLow*0.9999);
      }
      else {
        binMin = ax->FindBin(cut->fLow*1.0001);
        binMax = ax->FindBin(cut->fHigh*0.9999);
      }
      ax->SetRange(binMin, binMax);
    }
  }
  // axis 0 is assumed to be the invariant-mass axis
  Int_t binMin=fTHn->GetAxis(0)->FindBin(fCutSet->fFitLow*1.0001);
  Int_t binMax=fTHn->GetAxis(0)->FindBin(fCutSet->fFitHigh*0.9999);
  fTHn->GetAxis(0)->SetRange(binMin,binMax);
  TH1F* hProj = (TH1F*)fTHn->Projection(0);

  return hProj;
}


void AliHFCutVarFDsubMassFitter::FitMassHistogram() {
  fFitter = new AliHFMassFitter(fMassHist, fMassHist->GetBinLowEdge(1),
                                fMassHist->GetBinLowEdge(fMassHist->GetNbinsX()),
                                1, fCutSet->fFitTypeBkg, fCutSet->fFitTypeSig);

  if (fCutSet->fFitMean>0.) fFitter->SetFixGaussianMean(fCutSet->fFitMean);
  else fFitter->SetInitialGaussianMean(fCutSet->fFitMean*-1.);
  if (fCutSet->fFitSigma>0.) fFitter->SetFixGaussianSigma(fCutSet->fFitSigma);
  else fFitter->SetInitialGaussianSigma(fCutSet->fFitSigma*-1.);

  fFitter->SetReflectionSigmaFactor(fCutSet->fFitReflFactor);

  if (fFitter->MassFitter(kFALSE)) {
    fSig    = fFitter->GetRawYield();
    fSigErr = fFitter->GetRawYieldError();
    fSigma = fFitter->GetSigma();
    fSigmaErr = fFitter->GetSigmaUncertainty();
    fMean = fFitter->GetMean();
    fMeanErr = fFitter->GetMeanUncertainty();
    fChiSquare = fFitter->GetChiSquare();
    fRedChiSquare = fFitter->GetReducedChiSquare();
    fFitter->Background(3,fBkg,fBkgErr);
  }
}


AliHFCutVarFDsubMassFitter::AliHFCutVarFDsubMassFitter(const AliHFCutVarFDsubMassFitter& mf)
  : TObject()
  , fTHn(mf.fTHn)
  , fAxes(mf.fAxes)
  , fCutSet(mf.fCutSet)
  , fMassHist(mf.fMassHist)
  , fFitter(mf.fFitter)
  , fSig(mf.fSig)
  , fSigErr(mf.fSigErr)
  , fBkg(mf.fBkg)
  , fBkgErr(mf.fBkgErr)
  , fSigma(mf.fSigma)
  , fSigmaErr(mf.fSigmaErr)
  , fMean(mf.fMean)
  , fMeanErr(mf.fMeanErr)
  , fChiSquare(mf.fChiSquare)
{
  /// Copy constructor
}


AliHFCutVarFDsubMassFitter AliHFCutVarFDsubMassFitter::operator=(const AliHFCutVarFDsubMassFitter& mf)
{
  /// Assignment operator
  if (this != &mf) {
    fTHn = mf.fTHn;
    fAxes = mf.fAxes;
    fCutSet = mf.fCutSet;
    fMassHist = mf.fMassHist;
    fFitter = mf.fFitter;
    fSig = mf.fSig;
    fSigErr = mf.fSigErr;
    fSig = mf.fSig;
    fBkgErr = mf.fBkgErr;
    fSigma = mf.fSigma;
    fSigmaErr = mf.fSigmaErr;
    fMean = mf.fMean;
    fMeanErr = mf.fMeanErr;
    fChiSquare = mf.fChiSquare;
  }
  return *this;
}
