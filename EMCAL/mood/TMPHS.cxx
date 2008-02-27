#include "TMPHS.h"

ClassImp(TMPHS)

//_____________________________________________________________________________
TMPHS::TMPHS(const TGWindow *p, UInt_t w, UInt_t h) : 
  TMCal(p, w, h, AliCaloCalibPedestal::kPhos)
{
}
