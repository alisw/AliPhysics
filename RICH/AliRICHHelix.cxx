#include "AliRICHHelix.h"
#include <AliLog.h>

ClassImp(AliRICHHelix)

void  AliRICHHelix::Print(Option_t *opt) const
{
// Debug printout
  fX0.Print(opt);fP0.Print(opt);
  AliInfo("Point of interest:");
  fX.Print(opt); fP.Print(opt);
}//Print()
