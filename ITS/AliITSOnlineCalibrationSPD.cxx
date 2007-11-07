#include "AliITSOnlineCalibrationSPD.h"

///////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                            //
// Implementation of the online container for dead and noisy pixels. //
//                                                                   //
///////////////////////////////////////////////////////////////////////

ClassImp(AliITSOnlineCalibrationSPD)	

AliITSOnlineCalibrationSPD::AliITSOnlineCalibrationSPD():
fEqNr(0),
fNrBad(0),
fBadChannels(0)
{}
//____________________________________________________________________________
Int_t AliITSOnlineCalibrationSPD::GetKeyAt(UInt_t index) const {
  // Get key of index-th bad pixel
  if (index<fNrBad) {
    return fBadChannels.At(index);
  }
  return -1;
}
