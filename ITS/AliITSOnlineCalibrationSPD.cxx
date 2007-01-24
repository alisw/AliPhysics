#include "AliITSOnlineCalibrationSPD.h"

///////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                            //
// Implementation of the online container for dead and noisy pixels. //
//                                                                   //
///////////////////////////////////////////////////////////////////////

ClassImp(AliITSOnlineCalibrationSPD)	

AliITSOnlineCalibrationSPD::AliITSOnlineCalibrationSPD():
fModuleNr(0),
fNrDead(0),
fDeadChannels(0),
fNrNoisy(0),
fNoisyChannels(0)
{}
//_________________________________________________________________________
void AliITSOnlineCalibrationSPD::AddDead(UInt_t col, UInt_t row) {
  //
  // Add a dead channel to fDeadChannel array
  //
  fDeadChannels.Set(fNrDead*2+2);
  fDeadChannels.AddAt(col,fNrDead*2);
  fDeadChannels.AddAt(row,fNrDead*2+1);
  fNrDead++;
}
//_________________________________________________________________________
Int_t AliITSOnlineCalibrationSPD::GetDeadColAt(UInt_t index) const {
  // 
  // Returns column of index-th dead channel
  //
  if (index<fNrDead) {
    return fDeadChannels.At(index*2);
  }
  return -1;
}
//_________________________________________________________________________
Int_t AliITSOnlineCalibrationSPD::GetDeadRowAt(UInt_t index) const {
  // 
  // Returns row of index-th dead channel
  //
  if (index<fNrDead) {
    return fDeadChannels.At(index*2+1);
  }
  return -1;
}
//_________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPD::IsPixelDead(Int_t col, Int_t row) const {
  //
  // Check if pixel (col,row) is dead
  //
  for (UInt_t i=0; i<fNrDead; i++) { 
    if (fDeadChannels.At(i*2)==col && fDeadChannels.At(i*2+1)==row) {
      return kTRUE;
    }
  }
  return kFALSE;
}
//____________________________________________________________________________
void AliITSOnlineCalibrationSPD::AddNoisy(UInt_t col, UInt_t row) {
  //
  // add noisy pixel 
  //
  fDeadChannels.Set(fNrNoisy*2+2);
  fNoisyChannels.AddAt(col,fNrNoisy*2);
  fNoisyChannels.AddAt(row,fNrNoisy*2+1);
  fNrNoisy++;
}
//____________________________________________________________________________
Int_t AliITSOnlineCalibrationSPD::GetNoisyColAt(UInt_t index) const {
  //
  // Get column of index-th noisy pixel
  //
  if (index<fNrNoisy) {
    return fNoisyChannels.At(index*2);
  }
  return -1;
}
//____________________________________________________________________________
Int_t AliITSOnlineCalibrationSPD::GetNoisyRowAt(UInt_t index) const {
  //
  // Get row of index-th noisy pixel
  //
  if (index<fNrNoisy) {
    return fNoisyChannels.At(index*2+1);
  }
  return -1;
}
//____________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPD::IsPixelNoisy(Int_t col, Int_t row) const {
  //
  // Check if pixel (col,row) is noisy
  //
  for (UInt_t i=0; i<fNrNoisy; i++) { 
    if (fNoisyChannels.At(i*2)==col && fNoisyChannels.At(i*2+1)==row) {
      return kTRUE;
    }
  }
  return kFALSE;
}
