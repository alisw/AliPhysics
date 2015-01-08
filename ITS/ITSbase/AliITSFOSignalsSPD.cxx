/////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                          //
//                                                                 //
// This class is used to store information on generated Fast-OR    //
// signals. 1200 bits, one per pixel chip.                         //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliITSFOSignalsSPD.h"

ClassImp(AliITSFOSignalsSPD)

AliITSFOSignalsSPD::AliITSFOSignalsSPD() :
  TObject(), fSignals(1200)
{
  // default constructor
}
//______________________________________________________________________
AliITSFOSignalsSPD::~AliITSFOSignalsSPD() {}
//______________________________________________________________________
AliITSFOSignalsSPD::AliITSFOSignalsSPD(const AliITSFOSignalsSPD& fo):
  TObject(), fSignals(fo.fSignals)
{
  // copy constructor
}
//______________________________________________________________________
AliITSFOSignalsSPD& AliITSFOSignalsSPD::operator=(const AliITSFOSignalsSPD& fo) {
  // assignment operator
  if (this!=&fo) {
    fSignals = fo.fSignals;
  }
  return *this;
}
//______________________________________________________________________
void AliITSFOSignalsSPD::SetSignal(UInt_t eq, UInt_t hs, UInt_t chip, Bool_t setVal) {
  // Set 0 or 1 for a specific chip
  fSignals.SetBitNumber(GetChipKey(eq,hs,chip),setVal);
}
//______________________________________________________________________
Bool_t AliITSFOSignalsSPD::GetSignal(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // check if a specific chip has a signal
  return fSignals.TestBitNumber(GetChipKey(eq,hs,chip));
}
//______________________________________________________________________
Bool_t AliITSFOSignalsSPD::GetNextSignal(Int_t& eq, Int_t& hs, Int_t& chip) const {
  // Returns true if a signal was found (start looking after the bit number
  // corresponding to the input parameters eq,hs,chip).
  // If either of eq,hs,chip < 0 , start from beginning of TBits array.
  // See example of usage in DumpSignals method.
  UInt_t searchIndex;
  if (eq<0 || hs<0 || chip<0) searchIndex = 0;
  else searchIndex = GetChipKey(eq, hs, chip) + 1;
  UInt_t nextIndex = fSignals.FirstSetBit(searchIndex);
  if (nextIndex==1200) return kFALSE;
  GetChipFromKey(nextIndex, eq, hs, chip);
  return kTRUE;
}
//__________________________________________________________________________________
void AliITSFOSignalsSPD::DumpSignals() {
  // print a list of the chips which have a signal
  printf("These chips (given in eq,hs,chip) have a signal:\n");
  UInt_t nrSignals=0;
  Int_t eq   = -1;
  Int_t hs   = -1;
  Int_t chip = -1;
  while (GetNextSignal(eq,hs,chip)) {
    printf("%d,%d,%d\n",eq,hs,chip);
    nrSignals++;
  }
  printf("In total %d signals.\n",nrSignals);
}
//______________________________________________________________________
UInt_t AliITSFOSignalsSPD::GetChipKey(Int_t eq, Int_t hs, Int_t chip) const {
  // translates eq,hs,chip numbers into one integer key (0-1199)
  if (eq>=20 || eq<0 || hs>=6 || hs<0 || chip>=10 || chip<0) {
    Error("AliITSFOSignalsSPD::GetChipKey", "eq,hs,chip = %d,%d,%d out of range",eq,hs,chip);
    return 0;
  }
  return eq*60 + hs*10 + chip;
}
//__________________________________________________________________________________
void AliITSFOSignalsSPD::GetChipFromKey(UInt_t key, Int_t& eq, Int_t& hs, Int_t& chip) const {
  // translates a chip key back into eq,hs,chip numbers
  if (key>=1200) {
    Error("AliITSFOSignalsSPD::GetChipFromKey", "key = %d out of range", key);
    return;
  }
  eq   = key/60;
  hs   = (key%60)/10;
  chip = key%10;
}
