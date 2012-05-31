/////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                          //
//                                                                 //
// This class is used to store Fast-OR efficiency values in OCDB.  //
// One value per pixel chip in this base class (if per column      //
// accuracy is needed, use AliITSFOEfficiencySPDColumn class).     //
// The values are the probability that a pixel hit will generate a //
// fast-OR signal.                                                 //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliITSFOEfficiencySPD.h"

AliITSFOEfficiencySPD::AliITSFOEfficiencySPD() :
  TObject()
{
  // default constructor, puts all efficiency values to 100%
  ResetValues();
}
//______________________________________________________________________
AliITSFOEfficiencySPD::AliITSFOEfficiencySPD(const AliITSFOEfficiencySPD& foEff) :
  TObject() 
{  
  // copy constructor, copy the array values from input object
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	fChipEfficiency[eq][hs][chip] = foEff.fChipEfficiency[eq][hs][chip];
      }
    }
  }
}
//______________________________________________________________________
AliITSFOEfficiencySPD::~AliITSFOEfficiencySPD() {}
//______________________________________________________________________
void AliITSFOEfficiencySPD::ResetValues() {
  // Set all efficiency values to 100%
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	fChipEfficiency[eq][hs][chip] = 1;
      }
    }
  }
}
//______________________________________________________________________
AliITSFOEfficiencySPD& AliITSFOEfficiencySPD::operator=(const AliITSFOEfficiencySPD& foEff) {
  // assignment operator
  if (this!=&foEff) {
    for (UInt_t eq=0; eq<20; eq++) {
      for (UInt_t hs=0; hs<6; hs++) {
	for (UInt_t chip=0; chip<10; chip++) {
	  fChipEfficiency[eq][hs][chip] = foEff.fChipEfficiency[eq][hs][chip];
	}
      }
    }
  }
  return *this;
}
//______________________________________________________________________
void AliITSFOEfficiencySPD::SetChipEfficiency(UInt_t eq, UInt_t hs, UInt_t chip, Float_t value) {
  // set a chip efficiency value
  if (eq>=20) {
    Error("AliITSFOEfficiencySPD::SetChipEfficiency", "eq (%d) out of bounds.",eq);
    return;
  }
  if (hs>=6) {
    Error("AliITSFOEfficiencySPD::SetChipEfficiency", "hs (%d) out of bounds.",hs);
    return;
  }
  if (chip>=10) {
    Error("AliITSFOEfficiencySPD::SetChipEfficiency", "chip (%d) out of bounds.",chip);
    return;
  }
  
  fChipEfficiency[eq][hs][chip] = value;
}
//______________________________________________________________________
Float_t AliITSFOEfficiencySPD::GetChipEfficiency(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // get a chip efficiency value
  if (eq>=20) {
    Error("AliITSFOEfficiencySPD::GetChipEfficiency", "eq (%d) out of bounds.",eq);
    return 0;
  }
  if (hs>=6) {
    Error("AliITSFOEfficiencySPD::GetChipEfficiency", "hs (%d) out of bounds.",hs);
    return 0;
  }
  if (chip>=10) {
    Error("AliITSFOEfficiencySPD::GetChipEfficiency", "chip (%d) out of bounds.",chip);
    return 0;
  }

  return fChipEfficiency[eq][hs][chip];
}

