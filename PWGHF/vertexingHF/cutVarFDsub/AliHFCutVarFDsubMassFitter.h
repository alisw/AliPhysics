#ifndef ALIHFCUTVARFDSUBMASSFITTER_H
#define ALIHFCUTVARFDSUBMASSFITTER_H
/// \class AliHFCutVarFDsubMassFitter
/// \brief Cut for the cut variation feed down method analysis
///
///
///
///
/// \author Felix Reidt <felix.reidt@cern.ch>, CERN
/// \author Fabrizio Grosa <grosa@to.infn.it>, INFN Torino
/// \date Aug 17, 2015

#include "TObject.h"
#include "THnSparse.h"

class TList;
class TH1F;
class AliHFMassFitter;
class AliHFCutVarFDsubCutSet;


class AliHFCutVarFDsubMassFitter : public TObject {
protected:
  THnSparseF* fTHn;
  TList* fAxes;
  AliHFCutVarFDsubCutSet* fCutSet;

  TH1F* fMassHist;          //!<! Invariant mass histogram
  AliHFMassFitter* fFitter; //!<! Helper class used to fit the invariant mass distribution

  Double_t fSig;
  Double_t fSigErr;
  Double_t fBkg;
  Double_t fBkgErr;

  void ObtainMassHistogram();
  void FitMassHistogram();

  AliHFCutVarFDsubMassFitter(const AliHFCutVarFDsubMassFitter& mf); /// Copy constructor
  AliHFCutVarFDsubMassFitter operator=(const AliHFCutVarFDsubMassFitter& mf); /// Assignment operator

public:
  AliHFCutVarFDsubMassFitter(); /// Default constructor
  AliHFCutVarFDsubMassFitter(THnSparseF* thn, TList* axes, AliHFCutVarFDsubCutSet* cutSet); /// Constructor
  ~AliHFCutVarFDsubMassFitter(); /// Destructor

  Double_t GetSig() { return fSig; };
  Double_t GetSigErr() { return fSigErr; };
  Double_t GetBkg() { return fBkg; };
  Double_t GetBkgErr() { return fBkgErr; };

  AliHFMassFitter* GetFitter() { return fFitter; };

/// \cond CLASSDEF
  ClassDef(AliHFCutVarFDsubMassFitter, 1);
  /// \endcond
};
#endif // ALIHFCUTVARFDSUBMASSFITTER_H
