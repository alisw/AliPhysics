#include "AliUtilTOFParams.h"

////////////////////////////////////////////////////////////////////////////
///                                                                       //
///                                                                       //
/// Set of parameters and utilities for the Pi/Ka/Pr analysis with TOF    //
///                                                                       //
///                                                                       //
/// Authors:                                                              //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                        //
////////////////////////////////////////////////////////////////////////////

namespace AliUtilTOFParams {

const Double_t* CutValues[nCuts] = { CutValueTPCRows, CutValueMaxChi2, CutValueMaxDCAz, CutValueMaxDCAxy, CutValueGeo }; //Values of the cuts for each cut type

const Double_t* PtRange[2] = { PtRange_Hi, PtRange_pp };

//________________________________________________________________________
Bool_t DataInRange(const Double_t x, const Double_t min, const Double_t max)
{
  if (x >= max || x < min) {
    ::Info("DataInRange", "x (%f) is out of range [%f, %f]", x, min, max);
    return kFALSE;
  } else
    return kTRUE;
}

//________________________________________________________________________
Int_t BinData(const Double_t x, const Double_t min, const Double_t max, const Int_t base)
{
  const Double_t width = (max - min) / base;
  return TMath::Floor((x - min) / width);
}

//________________________________________________________________________
Double_t GetBinnedData(const Int_t x, const Double_t min, const Double_t max, const Int_t base)
{
  const Double_t width = (max - min) / base;
  return (0.5 + x) * width + min;
}
}
