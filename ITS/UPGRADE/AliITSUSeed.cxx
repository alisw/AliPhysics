#include <TString.h>
#include <TMath.h>
#include "AliITSUSeed.h"
using namespace TMath;

ClassImp(AliITSUSeed)

//_________________________________________________________________________
AliITSUSeed::AliITSUSeed() 
: fHitsPattern(0)
  ,fClID(0)
  ,fChi2Glo(0)
  ,fChi2Cl(0)
  ,fParent(0)
{
  // def c-tor
  ResetFMatrix();
}

//_________________________________________________________________________
AliITSUSeed::~AliITSUSeed()
{
  // d-rot
}

//_________________________________________________________________________
AliITSUSeed::AliITSUSeed(const AliITSUSeed& src) 
  :AliExternalTrackParam(src)
  ,fHitsPattern(src.fHitsPattern)
  ,fClID(src.fClID)
  ,fChi2Glo(src.fChi2Glo)
  ,fChi2Cl(src.fChi2Cl)
  ,fParent(src.fParent) 
{
  // def c-tor
  for (int i=kNFElem;i--;) fFMatrix[i] = src.fFMatrix[i];
}

//_________________________________________________________________________
AliITSUSeed &AliITSUSeed::operator=(const AliITSUSeed& src) 
{
  // def c-tor
  if (this == &src) return *this;
  fClID        = src.fClID;
  fHitsPattern = src.fHitsPattern;
  fChi2Glo     = src.fChi2Glo;
  fChi2Cl      = src.fChi2Cl;
  fParent      = src.fParent;
  for (int i=kNFElem;i--;) fFMatrix[i] = src.fFMatrix[i];
  AliExternalTrackParam::operator=(src);
  return *this;
}

//_________________________________________________________________________
void AliITSUSeed::Print(Option_t* opt) const
{
  // print seed info
  int lr,cl = GetLrCluster(lr);
  printf("%cLr%d Cl:%4d Chi2Glo:%7.2f(%7.2f) Chi2Cl:",IsKilled() ? '-':' ',
	 lr,cl,GetChi2Glo(),GetChi2GloNrm());
  cl<0 ? printf("   NA  ") : printf("%7.2f",GetChi2Cl());
  printf(" |"); 
  for (int i=0;i<=12;i++) printf("%c",HasClusterOnLayer(i) ? '+':'-'); printf("|\n");
  TString opts = opt; opts.ToLower();
  if (opts.Contains("etp")) AliExternalTrackParam::Print();
  if (opts.Contains("parent") && GetParent()) GetParent()->Print(opt);
}

//______________________________________________________________________________
Float_t AliITSUSeed::GetChi2GloNrm() const
{
  int ndf = 2*GetNLayersHit() - 5;
  return ndf>0 ? fChi2Glo/ndf : fChi2Glo;
}


//______________________________________________________________________________
Int_t AliITSUSeed::Compare(const TObject* obj)  const
{
  // compare clusters accodring to specific mode
  const AliITSUSeed* sd = (const AliITSUSeed*)obj;
  const Float_t kTol = 1e-5;
  if (!IsKilled() && sd->IsKilled()) return -1;
  if ( IsKilled() &&!sd->IsKilled()) return  1;
  //
  if      (GetChi2Glo()+kTol<sd->GetChi2Glo()) return -1;
  else if (GetChi2Glo()-kTol>sd->GetChi2Glo()) return  1;
  return 0;
}

//______________________________________________________________________________
Bool_t AliITSUSeed::IsEqual(const TObject* obj)  const
{
  // compare clusters accodring to specific mode
  const AliITSUSeed* sd = (const AliITSUSeed*)obj;
  const Float_t kTol = 1e-5;
  if (IsKilled() != sd->IsKilled()) return kFALSE;
  return Abs(GetChi2Glo() - sd->GetChi2Glo())<kTol;
}

//______________________________________________________________________________
Bool_t AliITSUSeed::PropagateToX(Double_t xk, Double_t b) 
{
  // Propagate this track to the plane X=xk (cm) in the field "b" (kG)
  Double_t dx=xk-fX;
  if (TMath::Abs(dx)<=kAlmost0)  return kTRUE;

  Double_t crv=GetC(b);
  if (TMath::Abs(b) < kAlmost0Field) crv=0.;
  Double_t x2r = crv*dx;
  Double_t f1=fP[2], f2=f1 + x2r;
  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  if (TMath::Abs(f2) >= kAlmost1) return kFALSE;
  if (TMath::Abs(fP[4])< kAlmost0) return kFALSE;

  Double_t &fP0=fP[0], &fP1=fP[1], &fP2=fP[2], &fP3=fP[3], &fP4=fP[4];
  Double_t 
  &fC00=fC[0],
  &fC10=fC[1],   &fC11=fC[2],  
  &fC20=fC[3],   &fC21=fC[4],   &fC22=fC[5],
  &fC30=fC[6],   &fC31=fC[7],   &fC32=fC[8],   &fC33=fC[9],  
  &fC40=fC[10],  &fC41=fC[11],  &fC42=fC[12],  &fC43=fC[13], &fC44=fC[14];

  Double_t r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2));
  if (TMath::Abs(r1)<kAlmost0)  return kFALSE;
  if (TMath::Abs(r2)<kAlmost0)  return kFALSE;

  fX=xk;
  double dy2dx = (f1+f2)/(r1+r2);
  fP0 += dx*dy2dx;
  if (TMath::Abs(x2r)<0.05) {
    fP1 += dx*(r2 + f2*dy2dx)*fP3;  // Many thanks to P.Hristov !
    fP2 += x2r;
  }
  else { 
    // for small dx/R the linear apporximation of the arc by the segment is OK,
    // but at large dx/R the error is very large and leads to incorrect Z propagation
    // angle traversed delta = 2*asin(dist_start_end / R / 2), hence the arc is: R*deltaPhi
    // The dist_start_end is obtained from sqrt(dx^2+dy^2) = x/(r1+r2)*sqrt(2+f1*f2+r1*r2)
    // Similarly, the rotation angle in linear in dx only for dx<<R
    double chord = dx*TMath::Sqrt(1+dy2dx*dy2dx);   // distance from old position to new one
    double rot = 2*TMath::ASin(0.5*chord*crv); // angular difference seen from the circle center
    fP1 += rot/crv*fP3;
    fP2  = TMath::Sin(rot + TMath::ASin(fP2));
  }

  //f = F - 1
  double r1i = 1./r1;
  double r2i = 1./r2;
  double tg1 = f1*r1i;
  double tg2 = f2*r2i;
  double v0 = 1. + dy2dx*tg2;
  double v1 = (r1i+r2i)*(dy2dx*(tg1+tg2)+2);
  double v2 = (r1i+r2i)*v0;
  //
  double f24 = dx*crv/fP4;
  double f02 = dx*v1;
  double f04 = dx*v2*f24;
  double f12 = dx*fP3*    (f2*v1+dy2dx-tg2);
  double f13 = dx*r2*v0;
  double f14 = dx*f24*fP3*(f2*v2+dy2dx-tg2);
  //
  //b = C*ft
  Double_t b00=f02*fC20 + f04*fC40, b01=f12*fC20 + f14*fC40 + f13*fC30;
  Double_t b02=f24*fC40;
  Double_t b10=f02*fC21 + f04*fC41, b11=f12*fC21 + f14*fC41 + f13*fC31;
  Double_t b12=f24*fC41;
  Double_t b20=f02*fC22 + f04*fC42, b21=f12*fC22 + f14*fC42 + f13*fC32;
  Double_t b22=f24*fC42;
  Double_t b40=f02*fC42 + f04*fC44, b41=f12*fC42 + f14*fC44 + f13*fC43;
  Double_t b42=f24*fC44;
  Double_t b30=f02*fC32 + f04*fC43, b31=f12*fC32 + f14*fC43 + f13*fC33;
  Double_t b32=f24*fC43;
  
  //a = f*b = f*C*ft
  Double_t a00=f02*b20+f04*b40,a01=f02*b21+f04*b41,a02=f02*b22+f04*b42;
  Double_t a11=f12*b21+f14*b41+f13*b31,a12=f12*b22+f14*b42+f13*b32;
  Double_t a22=f24*b42;

  //F*C*Ft = C + (b + bt + a)
  fC00 += b00 + b00 + a00;
  fC10 += b10 + b01 + a01; 
  fC20 += b20 + b02 + a02;
  fC30 += b30;
  fC40 += b40;
  fC11 += b11 + b11 + a11;
  fC21 += b21 + b12 + a12;
  fC31 += b31; 
  fC41 += b41;
  fC22 += b22 + b22 + a22;
  fC32 += b32;
  fC42 += b42;
  //
  // update stored transformation matrix   F = Fnew*Fold
  fFMatrix[kF04] += f04 + f24*fFMatrix[kF02];
  fFMatrix[kF14] += f14 + f24*fFMatrix[kF12];
  fFMatrix[kF02] += f02;
  fFMatrix[kF12] += f12;
  fFMatrix[kF13] += f13;
  fFMatrix[kF24] += f24;
  //
  CheckCovariance();

  return kTRUE;
}
