/////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                          //
//                                                                 //
// This class is used to store Fast-OR efficiency values in OCDB.  //
// One value per pixel chip column in this daughter class.         //
// The values are the probability that a pixel hit will generate a //
// fast-OR signal.                                                 //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliITSFOEfficiencySPDColumn.h"

AliITSFOEfficiencySPDColumn::AliITSFOEfficiencySPDColumn() : 
  AliITSFOEfficiencySPD() 
{
  // default constructor, sets all efficiency values to 100%
  ResetValues();
}
//______________________________________________________________________
AliITSFOEfficiencySPDColumn::AliITSFOEfficiencySPDColumn(const AliITSFOEfficiencySPDColumn& foEff) :
  AliITSFOEfficiencySPD()
{  
  // copy constructor, copy the array values from input object
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	for (UInt_t col=0; col<32; col++) {
	  fColumnEfficiency[eq][hs][chip][col] = foEff.fColumnEfficiency[eq][hs][chip][col];
	}
      }
    }
  }
}
//______________________________________________________________________
AliITSFOEfficiencySPDColumn::~AliITSFOEfficiencySPDColumn() {}
//______________________________________________________________________
AliITSFOEfficiencySPDColumn& AliITSFOEfficiencySPDColumn::operator=(const AliITSFOEfficiencySPDColumn& foEff) {
  // assignment operator
  if (this!=&foEff) {
    for (UInt_t eq=0; eq<20; eq++) {
      for (UInt_t hs=0; hs<6; hs++) {
	for (UInt_t chip=0; chip<10; chip++) {
	  for (UInt_t col=0; col<32; col++) {
	    fColumnEfficiency[eq][hs][chip][col] = foEff.fColumnEfficiency[eq][hs][chip][col];
	  }
	}
      }
    }
  }
  return *this;
}
//______________________________________________________________________
void AliITSFOEfficiencySPDColumn::ResetValues() {
  // set all efficiency values to 100%
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	for (UInt_t col=0; col<32; col++) {
	  fColumnEfficiency[eq][hs][chip][col] = 1;
	}
      }
    }
  }
}
//______________________________________________________________________
void AliITSFOEfficiencySPDColumn::SetColumnEfficiency(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, Float_t value) {
  // set a column efficiency value
  if (eq>=20) {
    Error("AliITSFOEfficiencySPDColumn::SetColumnEfficiency", "eq (%d) out of bounds.",eq);
    return;
  }
  if (hs>=6) {
    Error("AliITSFOEfficiencySPDColumn::SetColumnEfficiency", "hs (%d) out of bounds.",hs);
    return;
  }
  if (chip>=10) {
    Error("AliITSFOEfficiencySPDColumn::SetColumnEfficiency", "chip (%d) out of bounds.",chip);
    return;
  }
  if (col>=32) {
    Error("AliITSFOEfficiencySPDColumn::SetColumnEfficiency", "col (%d) out of bounds.",col);
    return;
  }

  fColumnEfficiency[eq][hs][chip][col] = value;
}
//______________________________________________________________________
Float_t AliITSFOEfficiencySPDColumn::GetColumnEfficiency(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col) const {
  // get a column efficiency value
  if (eq>=20) {
    Error("AliITSFOEfficiencySPDColumn::GetEfficiency", "eq (%d) out of bounds.",eq);
    return 0;
  }
  if (hs>=6) {
    Error("AliITSFOEfficiencySPDColumn::GetEfficiency", "hs (%d) out of bounds.",hs);
    return 0;
  }
  if (chip>=10) {
    Error("AliITSFOEfficiencySPDColumn::GetEfficiency", "chip (%d) out of bounds.",chip);
    return 0;
  }
  if (col>=32) {
    Error("AliITSFOEfficiencySPDColumn::GetEfficiency", "col (%d) out of bounds.",col);
    return 0;
  }

  return fColumnEfficiency[eq][hs][chip][col];
}

