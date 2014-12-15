/////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                          //
//                                                                 //
// This class is used to store Fast-OR noise values in OCDB.       //
// One value per pixel chip.                                       //
// The values are the probability that a pixel chip will generate  //
// a fast-OR signal independently (originating from noise).        //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliITSFONoiseSPD.h"

AliITSFONoiseSPD::AliITSFONoiseSPD() :
  TObject()
{
  // default constructor, sets all noise values to 0%
  ResetValues();
}
//______________________________________________________________________
AliITSFONoiseSPD::AliITSFONoiseSPD(const AliITSFONoiseSPD& foNoi) :
  TObject()
{
  // copy constructor, copy the array values from input object
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	fChipNoise[eq][hs][chip] = foNoi.fChipNoise[eq][hs][chip];
      }
    }
  }
}
//______________________________________________________________________
AliITSFONoiseSPD::~AliITSFONoiseSPD() {}
//______________________________________________________________________
void AliITSFONoiseSPD::ResetValues() {
  // set all noise values to 0%
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	fChipNoise[eq][hs][chip] = 0;
      }
    }
  }
}
//______________________________________________________________________
AliITSFONoiseSPD& AliITSFONoiseSPD::operator=(const AliITSFONoiseSPD& foNoi) {
  // assignment operator
  if (this!=&foNoi) {
    for (UInt_t eq=0; eq<20; eq++) {
      for (UInt_t hs=0; hs<6; hs++) {
	for (UInt_t chip=0; chip<10; chip++) {
	  fChipNoise[eq][hs][chip] = foNoi.fChipNoise[eq][hs][chip];
	}
      }
    }
  }
  return *this;
}
//______________________________________________________________________
void AliITSFONoiseSPD::SetChipNoise(UInt_t eq, UInt_t hs, UInt_t chip, Float_t value) {
  // set a chip noise value
  if (eq>=20) {
    Error("AliITSFONoiseSPD::SetChipNoise", "eq (%d) out of bounds.",eq);
    return;
  }
  if (hs>=6) {
    Error("AliITSFONoiseSPD::SetChipNoise", "hs (%d) out of bounds.",hs);
    return;
  }
  if (chip>=10) {
    Error("AliITSFONoiseSPD::SetChipNoise", "chip (%d) out of bounds.",chip);
    return;
  }
  
  fChipNoise[eq][hs][chip] = value;
}
//______________________________________________________________________
Float_t AliITSFONoiseSPD::GetChipNoise(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // get a chip noise value
  if (eq>=20) {
    Error("AliITSFONoiseSPD::GetChipNoise", "eq (%d) out of bounds.",eq);
    return 0;
  }
  if (hs>=6) {
    Error("AliITSFONoiseSPD::GetChipNoise", "hs (%d) out of bounds.",hs);
    return 0;
  }
  if (chip>=10) {
    Error("AliITSFONoiseSPD::GetChipNoise", "chip (%d) out of bounds.",chip);
    return 0;
  }

  return fChipNoise[eq][hs][chip];
}

