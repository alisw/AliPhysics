#include "AliHFCutVarFDsubCutSet.h"

#include "TList.h"
#include "TString.h"

/// \cond CLASSIMP
ClassImp(AliHFCutVarFDsubCutSet);
/// \endcond


AliHFCutVarFDsubCutSet::AliHFCutVarFDsubCutSet()
  : TObject()
  , fCuts(0x0)
  , fFitLow(-1.)
  , fFitHigh(-1.)
  , fFitMean(-1.)
  , fFitSigma(-1.)
  , fFitTypeBkg(0)
  , fFitTypeSig(0)
  , fFitReflFactor(0.)
  , fFitRebin(1.)
{
  /// Default constructor
}


AliHFCutVarFDsubCutSet::AliHFCutVarFDsubCutSet(Double_t fitLow, Double_t fitHigh,
                                               Double_t fitMean, Double_t fitSigma,
                                               Double_t fitRebin/*=1*/,
                                               Int_t fitTypeBkg/*=0*/, Int_t fitTypeSig/*=0*/,
                                               Double_t fitReflFactor/*=0.*/)
  : TObject()
  , fCuts(new TList())
  , fFitLow(fitLow)
  , fFitHigh(fitHigh)
  , fFitMean(fitMean)
  , fFitSigma(fitSigma)
  , fFitTypeBkg(fitTypeBkg)
  , fFitTypeSig(fitTypeSig)
  , fFitReflFactor(fitReflFactor)
  , fFitRebin(fitRebin)
{
  /// Constructor
  fCuts->SetOwner();
}


AliHFCutVarFDsubCutSet::~AliHFCutVarFDsubCutSet() {
  delete fCuts;
  fCuts = 0x0;
}


Int_t AliHFCutVarFDsubCutSet::GetEntries() {
  if (fCuts) return fCuts->GetEntries();
  else       return -1;
}


Int_t AliHFCutVarFDsubCutSet::AddCut(AliHFCutVarFDsubCut* cut) {
  if (fCuts) {
    fCuts->Add((TObject*)cut);
    return fCuts->GetEntries() - 1;
  }
  else return -1;
}


AliHFCutVarFDsubCut* AliHFCutVarFDsubCutSet::GetCut(Int_t iCut) {
  if (fCuts && fCuts->GetEntries()>iCut) {
    return (AliHFCutVarFDsubCut*)fCuts->At(iCut);
  }
  else return 0x0;
}


AliHFCutVarFDsubCutSet::AliHFCutVarFDsubCutSet(const AliHFCutVarFDsubCutSet& cutSet)
  : TObject()
  , fCuts((TList*)cutSet.fCuts->Clone())
  , fFitLow(cutSet.fFitLow)
  , fFitHigh(cutSet.fFitHigh)
  , fFitMean(cutSet.fFitMean)
  , fFitSigma(cutSet.fFitSigma)
  , fFitTypeBkg(cutSet.fFitTypeBkg)
  , fFitTypeSig(cutSet.fFitTypeSig)
  , fFitReflFactor(cutSet.fFitReflFactor)
  , fFitRebin(cutSet.fFitRebin)
{
    /// Copy constructor
}


AliHFCutVarFDsubCutSet AliHFCutVarFDsubCutSet::operator=(const AliHFCutVarFDsubCutSet& cutSet)
{
  /// Assignment operator
  if (this != &cutSet) {
    delete fCuts;
    fCuts = 0x0;
    fCuts = cutSet.fCuts;
    fFitLow = cutSet.fFitLow;
    fFitHigh = cutSet.fFitHigh;
    fFitMean = cutSet.fFitMean;
    fFitSigma = cutSet.fFitSigma;
    fFitTypeBkg = cutSet.fFitTypeBkg;
    fFitTypeSig = cutSet.fFitTypeSig;
    fFitReflFactor = cutSet.fFitReflFactor;
    fFitRebin = cutSet.fFitRebin;
  }
  return *this;
}
