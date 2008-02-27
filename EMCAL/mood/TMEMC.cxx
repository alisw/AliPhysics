#include "TMEMC.h"

ClassImp(TMEMC)

//_____________________________________________________________________________
TMEMC::TMEMC(const TGWindow *p, UInt_t w, UInt_t h) : 
  TMCal(p, w, h, AliCaloCalibPedestal::kEmCal)
{
}

