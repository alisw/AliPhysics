#ifndef ALIHFCUTVARFDSUBCUTSET_H
#define ALIHFCUTVARFDSUBCUTSET_H
/// \class AliHFCutVarFDsubCutSet
/// \brief Cut set (bin) for the cut variation feed down method analysis
///
///
///
///
/// \author Felix Reidt <felix.reidt@cern.ch>, CERN
/// \author Fabrizio Grosa <grosa@to.infn.it>, INFN Torino
/// \date Aug 17, 2015

#include "TObject.h"

class TList;
class AliHFCutVarFDsubCut;

class AliHFCutVarFDsubCutSet : public TObject {
protected:
  TList* fCuts; ///!<! Cuts

  AliHFCutVarFDsubCutSet(const AliHFCutVarFDsubCutSet& cutSet); /// Copy constructor
  AliHFCutVarFDsubCutSet operator=(const AliHFCutVarFDsubCutSet& cutSet); // Assignment operator

public:
  Double_t fFitLow;        ///< Lower fit invariant mass boundary
  Double_t fFitHigh;       ///< Upper fit invariant mass boundary
  /// fit parameters lead to fixed values if positive and variable values if negative
  Double_t fFitMean;       ///< Mean (-1.*mean leaves it variable but initialises it to mean)
  Double_t fFitSigma;      ///< Sigma (-1.*sigma leaves it variable but initialises it to sigma)
  Int_t    fFitTypeBkg;    ///< Background fit type
  Int_t    fFitTypeSig;    ///< Signal fit type
  Double_t fFitReflFactor; ///< Reflection sigma factor
  Double_t fFitRebin;      ///< Reflection sigma factor

  AliHFCutVarFDsubCutSet(); /// Default constructor
  AliHFCutVarFDsubCutSet(Double_t fitLow, Double_t fitHigh, Double_t fitMean, Double_t fitSigma,
                         Double_t fitRebin=1., Int_t fitTypeBkg=0, Int_t fitTypeSig=0,
                         Double_t fitReflFactor=0.); /// Constructor
  ~ AliHFCutVarFDsubCutSet(); /// Destructor

  Int_t GetEntries();
  Int_t AddCut(AliHFCutVarFDsubCut* cut);
  AliHFCutVarFDsubCut* GetCut(Int_t iCut);

  /// \cond CLASSDEF
  ClassDef(AliHFCutVarFDsubCutSet, 1);
  /// \endcond
};
#endif // ALIHFCUTVARFDSUBCUTSET_H
