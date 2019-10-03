#include "AliHFCutVarFDsubCut.h"

#include "TString.h"


/// \cond CLASSIMP
ClassImp(AliHFCutVarFDsubCut);
/// \endcond


AliHFCutVarFDsubCut::AliHFCutVarFDsubCut()
  : TObject()
  , fAxisId((UInt_t)-1)
  , fLow(1.)
  , fHigh(-1.)
{
  // Default constructor
}


AliHFCutVarFDsubCut::AliHFCutVarFDsubCut(UInt_t axisId, Double_t low, Double_t high)
  : TObject()
  , fAxisId(axisId)
  , fLow(low)
  , fHigh(high)
{
  // Constructor
}
