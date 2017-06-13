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
  const Double_t *CutValues[nCuts] = {CutValueTPCRows, CutValueMaxChi2, CutValueMaxDCAz, CutValueMaxDCAxy, CutValueGeo};//Values of the cuts for each cut type
  
  const TString *MCProduction[2] = {MCProduction_Hi, MCProduction_pp};
  const Double_t *PtRange[2] = {PtRange_Hi, PtRange_pp};
}
