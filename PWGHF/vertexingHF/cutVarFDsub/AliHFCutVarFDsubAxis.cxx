#include "AliHFCutVarFDsubAxis.h"

#include "TString.h"


/// \cond CLASSIMP
ClassImp(AliHFCutVarFDsubAxis);
/// \endcond


AliHFCutVarFDsubAxis::AliHFCutVarFDsubAxis()
  : TObject()
  , fAxisNoData((UInt_t)-1)
  , fAxisNoMCgenLevel((UInt_t)-1)
  , fAxisNoMCafterCuts((UInt_t)-1)
  , fAxisName("")
  , fSymmCut(kFALSE)
{
  // Default constructor
}


AliHFCutVarFDsubAxis::AliHFCutVarFDsubAxis(UInt_t axisNoData, UInt_t axisNoMCgenLevel, UInt_t axisNoMCafterCuts, TString axisName, Bool_t iscutsymm)
  : TObject()
  , fAxisNoData(axisNoData)
  , fAxisNoMCgenLevel(axisNoMCgenLevel)
  , fAxisNoMCafterCuts(axisNoMCafterCuts)
  , fAxisName(axisName)
  , fSymmCut(iscutsymm)
{
  // Constructor
}

UInt_t AliHFCutVarFDsubAxis::GetAxisNo(UInt_t thnType) {
  switch(thnType) {
  case kData:        return fAxisNoData;
  case kMCgenLevel:  return fAxisNoMCgenLevel;
  case kMCafterCuts: return fAxisNoMCafterCuts;
  default:           return (UInt_t)-1;
  }
  return (UInt_t)-1;
}

