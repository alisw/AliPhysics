#ifndef ALIHFCUTVARFDSUBAXIS_H
#define ALIHFCUTVARFDSUBAXIS_H
/// \class AliHFCutVarFDsubAxis
/// \brief Axis description for the cut variation feed down method analysis
///
///
///
///
/// \author Felix Reidt <felix.reidt@cern.ch>, CERN
/// \author Fabrizio Grosa <grosa@to.infn.it>, INFN Torino
/// \date Aug 17, 2015

#include "TObject.h"
#include "TString.h"

class AliHFCutVarFDsubAxis : public TObject {
protected:
  UInt_t fAxisNoData;        ///< Number of the THnSparse Data axis to cut on
  UInt_t fAxisNoMCgenLevel;  ///< Number of the THnSparse MC axis to cut on
  UInt_t fAxisNoMCafterCuts; ///< Number of the THnSparse MC axis to cut on
  TString fAxisName;         ///< Name of the axis (optional)
  Bool_t fSymmCut; ///< Flag to activate a two-region cut symmetric with respect to zero 
public:

  AliHFCutVarFDsubAxis(); /// Default constructor
  AliHFCutVarFDsubAxis(UInt_t axisNoData, UInt_t axisNoMCgenLevel, UInt_t axisNoMCafterCuts, TString axisName, Bool_t iscutsymm=kFALSE); /// Constructor

  TString GetAxisName() { return fAxisName; } /// Get the name of the axis
  UInt_t GetAxisNo(UInt_t thnType); /// Get the axis number
  Bool_t IsCutSymmetric() { return fSymmCut; } /// Get if the cut on this axis is a two-region cut symmetric with respect to zero
  
  enum { kData=0, kMCgenLevel, kMCafterCuts };

  /// \cond CLASSDEF
  ClassDef(AliHFCutVarFDsubAxis, 2);
  /// \endcond
};
#endif // ALIHFCUTVARFDSUBAXIS_H
