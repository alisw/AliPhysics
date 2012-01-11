////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD cluster summary info for performance                              //
//                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TMath.h"

#include "AliLog.h"
#include "AliTRDcluster.h"

#include "AliTRDclusterInfo.h"

ClassImp(AliTRDclusterInfo)

//_________________________________________________
AliTRDclusterInfo::AliTRDclusterInfo()
  : TObject()
  ,fDet(0xffff)
  ,fCol(0xff)
  ,fRow(0xff)
  ,fNpad(0)
  ,fPdg(0)
  ,fLbl(-1)
  ,fLocalTime(-100)
  ,fQ(0.)
  ,fX(0.)
  ,fY(0.)
  ,fYd(0.)
  ,fZ(0.)
  ,fdydx(0.)
  ,fdzdx(0.)
  ,fXd(0.)
  ,fYt(0.)
  ,fZt(0.)
  ,fdy(0.)
  ,fD(0.)
  ,fTilt(0.)
{
//  Constructor. Resets all fields.
  fCov[0] = 1.; fCov[1] = 0.;
  fCov[2] = 1.;
  fCovCl[0] = 1.; fCovCl[1] = 0.;
  fCovCl[2] = 1.;
  memset(fSignal, 0, 7*sizeof(Short_t));
}

//_________________________________________________
void AliTRDclusterInfo::SetCluster(AliTRDcluster *c)
{
// Load rec cluster data
  if(!c) return;
  fDet = c->GetDetector();
  fCol = c->GetPadCol();
  fRow = c->GetPadRow();
  fNpad= c->GetNPads();
  fX   = c->GetX();
  fY   = c->GetY();
  fZ   = c->GetZ();
  fQ   = TMath::Abs(c->GetQ());
  fLocalTime = c->GetLocalTimeBin();
  fYd  = c->GetCenter();
  fCovCl[0] = c->GetSigmaY2();
  fCovCl[1] = 0.;
  fCovCl[2] = c->GetSigmaZ2();
  memcpy(fSignal, c->GetSignals(), 7*sizeof(Short_t));
}

//_________________________________________________
void AliTRDclusterInfo::Print(Option_t */*opt*/) const
{
// Dump info
  printf("Det[%3d] Col[%3d] Row[%2d] X[%7.2f] Y[%7.2f] Z[%7.2f] Q[%7.2f] N[%d]\n", (fDet==0xffff ? -1 : fDet), (fCol==0xff ? -1 : fCol), (fRow==0xff ? -1 : fRow), fX, fY, fZ, fQ, fNpad);
  printf("\tPdg[%d] Lbl[%d] Yt[%7.2f] Zt[%7.2f]\n", fPdg, fLbl, fYt, fZt);
}
