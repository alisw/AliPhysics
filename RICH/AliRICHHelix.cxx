#include "AliRICHHelix.h"
#include <AliLog.h>

ClassImp(AliRICHHelix)

void  AliRICHHelix::Print(Option_t *) const
{
//Debug printout
  AliInfo(Form("Q=%i, in B(0,0,%5.2f) x0=(%5.2f,%5.2f,%5.2f) p0=(%5.2f,%5.2f,%5.2f)",fQ,fBz,
                            fX0.X(),fX0.Y(),fX0.Z(),
                            fP0.X(),fP0.Y(),fP0.Z()  ));
  AliInfo(Form("At length %7.2f gives x=(%5.2f,%5.2f,%5.2f) p=(%5.2f,%5.2f,%5.2f)",fLen,
                            fX.X(),fX.Y(),fX.Z(),
                            fP.X(),fP.Y(),fP.Z()  ));

}
