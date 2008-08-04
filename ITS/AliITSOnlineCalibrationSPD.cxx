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
fBadChannels(0),
fActiveEq(kTRUE),
fDeadEq(kFALSE)
{
  ActivateALL();
  UnSetDeadALL();
}
//____________________________________________________________________________
Int_t AliITSOnlineCalibrationSPD::GetKeyAt(UInt_t index) const {
  // Get key of index-th bad pixel
  if (index<fNrBad) {
    return fBadChannels.At(index);
  }
  return -1;
}
//____________________________________________________________________________
void AliITSOnlineCalibrationSPD::ActivateALL() {
  // activate eq, all hs, all chips
  ActivateEq();
  for (UInt_t hs=0; hs<6; hs++) {
    ActivateHS(hs);
    for (UInt_t chip=0; chip<10; chip++) {
      ActivateChip(hs,chip);
    }
  }
}
//____________________________________________________________________________
void AliITSOnlineCalibrationSPD::ActivateEq(Bool_t setval) {
  // activate this eq
  fActiveEq = setval;
}
//____________________________________________________________________________
void AliITSOnlineCalibrationSPD::ActivateHS(UInt_t hs, Bool_t setval) {
  // activate hs on this eq
  if (hs>=6) {
    Error("AliITSOnlineCalibrationSPD::ActivateHS", "hs (%d) out of bounds.",hs);
    return;
  }
  fActiveHS[hs] = setval;
}
//____________________________________________________________________________
void AliITSOnlineCalibrationSPD::ActivateChip(UInt_t hs, UInt_t chip, Bool_t setval) {
  // activate chip on this eq
  if (hs>=6 || chip>=10) {
    Error("AliITSOnlineCalibrationSPD::ActivateChip", "hs,chip (%d,%d) out of bounds.",hs,chip);
    return;
  }
  fActiveChip[hs*10+chip] = setval;
}
//____________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPD::IsActiveEq() const {
  // is this eq active?
  return fActiveEq;
}
//____________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPD::IsActiveHS(UInt_t hs) const {
  // is this hs active?
  if (hs>=6) {
    Error("AliITSOnlineCalibrationSPD::IsActiveHS", "hs (%d) out of bounds.",hs);
    return kFALSE;
  }
  return fActiveHS[hs];
}
//____________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPD::IsActiveChip(UInt_t hs, UInt_t chip) const {
  // is this chip active?
  if (hs>=6 || chip>=10) {
    Error("AliITSOnlineCalibrationSPD::IsActiveChip", "hs,chip (%d,%d) out of bounds.",hs,chip);
    return kFALSE;
  }
  return fActiveChip[hs*10+chip];
}
//____________________________________________________________________________
void AliITSOnlineCalibrationSPD::UnSetDeadALL() {
  // activate eq, all hs, all chips
  SetDeadEq(kFALSE);
  for (UInt_t hs=0; hs<6; hs++) {
    SetDeadHS(hs,kFALSE);
    for (UInt_t chip=0; chip<10; chip++) {
      SetDeadChip(hs,chip,kFALSE);
    }
  }
}
//____________________________________________________________________________
void AliITSOnlineCalibrationSPD::SetDeadEq(Bool_t setval) {
  // set this eq dead
  fDeadEq = setval;
}
//____________________________________________________________________________
void AliITSOnlineCalibrationSPD::SetDeadHS(UInt_t hs, Bool_t setval) {
  // set dead hs on this eq
  if (hs>=6) {
    Error("AliITSOnlineCalibrationSPD::SetDeadHS", "hs (%d) out of bounds.",hs);
    return;
  }
  fDeadHS[hs] = setval;
}
//____________________________________________________________________________
void AliITSOnlineCalibrationSPD::SetDeadChip(UInt_t hs, UInt_t chip, Bool_t setval) {
  // set dead chip on this eq
  if (hs>=6 || chip>=10) {
    Error("AliITSOnlineCalibrationSPD::SetDeadChip", "hs,chip (%d,%d) out of bounds.",hs,chip);
    return;
  }
  fDeadChip[hs*10+chip] = setval;
}
//____________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPD::IsDeadEq() const {
  // is this eq dead?
  return fDeadEq;
}
//____________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPD::IsDeadHS(UInt_t hs) const {
  // is this hs dead?
  if (hs>=6) {
    Error("AliITSOnlineCalibrationSPD::IsDeadHS", "hs (%d) out of bounds.",hs);
    return kFALSE;
  }
  return fDeadHS[hs];
}
//____________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPD::IsDeadChip(UInt_t hs, UInt_t chip) const {
  // is this chip dead?
  if (hs>=6 || chip>=10) {
    Error("AliITSOnlineCalibrationSPD::IsDeadChip", "hs,chip (%d,%d) out of bounds.",hs,chip);
    return kFALSE;
  }
  return fDeadChip[hs*10+chip];
}
