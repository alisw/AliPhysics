////////////////////////////////////////////////
//  Digits classes for set:FMD                //
////////////////////////////////////////////////


#include "AliFMDdigit.h"

ClassImp(AliFMDdigit)

AliFMDdigit::AliFMDdigit(Int_t *digits) {
  //
  // Creates a real data digit object
  //
  fNumOfDet       = digits[0];
  fNumOfSector    = digits[1];
  fNumOfRing      = digits[2];
  fNelectrons     = digits[3];
  fADC            = digits[4];
}

