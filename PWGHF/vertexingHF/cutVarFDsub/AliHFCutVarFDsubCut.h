#ifndef ALIHFCUTVARFDSUBCUT_H
#define ALIHFCUTVARFDSUBCUT_H
/// \class AliHFCutVarFDsubCut
/// \brief Cut for the cut variation feed down method analysis
///
///
///
///
/// \author Felix Reidt <felix.reidt@cern.ch>, CERN
/// \author Fabrizio Grosa <grosa@to.infn.it>, INFN Torino
/// \date Aug 17, 2015

#include "TObject.h"
#include "TString.h"

class AliHFCutVarFDsubCut : public TObject {
public:
  UInt_t fAxisId;     ///< Number of the THnSparse Data axis to cut on
  Double_t fLow;      ///< lower range limit
  Double_t fHigh;     ///< upper range limit

  AliHFCutVarFDsubCut(); /// Default constructor
  AliHFCutVarFDsubCut(UInt_t axisId, Double_t low, Double_t high); /// Constructor

  /// \cond CLASSDEF
  ClassDef(AliHFCutVarFDsubCut, 1);
  /// \endcond
};
#endif // ALIHFCUTVARFDSUBCUT_H
