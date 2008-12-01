#include "TMath.h"

#include "AliLog.h"
#include "AliTRDcluster.h"

#include "AliTRDclusterInfo.h"

ClassImp(AliTRDclusterInfo)

//_________________________________________________
AliTRDclusterInfo::AliTRDclusterInfo()
  : TObject()
  ,fDet(0xffff)
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
{
  fCov[0] = 1.; fCov[1] = 0.;
  fCov[2] = 1.;
  fCovCl[0] = 1.; fCovCl[1] = 0.;
  fCovCl[2] = 1.;
}

//_________________________________________________
void AliTRDclusterInfo::SetCluster(const AliTRDcluster *c, Float_t *cov)
{
  if(!c) return;
  fDet = c->GetDetector();
  fX   = c->GetX();
  fY   = c->GetY();
  fZ   = c->GetZ();
  fQ   = TMath::Abs(c->GetQ());
  fLocalTime = c->GetLocalTimeBin();
  fYd  = c->GetCenter();

  if(cov) memcpy(cov, fCovCl, 3*sizeof(Float_t));
}

//_________________________________________________
void AliTRDclusterInfo::Print(Option_t */*opt*/) const
{
  printf("Det[%3d] X[%7.2f] Y[%7.2f] Z[%7.2f] Q[%7.2f]\n", fDet==0xffff ? -1 : fDet, fX, fY, fZ, fQ);
  printf("\tPdg[%d] Lbl[%d] Yt[%7.2f] Zt[%7.2f]\n", fPdg, fLbl, fYt, fZt);
}
