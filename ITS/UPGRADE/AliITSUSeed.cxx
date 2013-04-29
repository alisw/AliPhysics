#include <TString.h>
#include <TMath.h>
#include "AliITSUSeed.h"
#include "AliLog.h"
#include "AliESDtrack.h"
using namespace TMath;

ClassImp(AliITSUSeed)

//_________________________________________________________________________
AliITSUSeed::AliITSUSeed() 
: fHitsPattern(0)
  ,fNChildren(0)
  ,fClID(0)
  ,fChi2Glo(0)
  ,fChi2Cl(0)
  ,fChi2Penalty(0)
  ,fChi2Match(0)
  ,fChi2ITSSA(0)
  ,fParent(0)
{
  // def c-tor
  ResetFMatrix();
}

//_________________________________________________________________________
AliITSUSeed::~AliITSUSeed()
{
  // d-tor
}

//_________________________________________________________________________
AliITSUSeed::AliITSUSeed(const AliITSUSeed& src) 
  :AliExternalTrackParam(src)
  ,fHitsPattern(src.fHitsPattern)
  ,fNChildren(src.fNChildren)
  ,fClID(src.fClID)
  ,fChi2Glo(src.fChi2Glo)
  ,fChi2Cl(src.fChi2Cl)
  ,fChi2Penalty(src.fChi2Penalty)
  ,fChi2Match(src.fChi2Match)
  ,fChi2ITSSA(src.fChi2ITSSA)
  ,fParent(src.fParent) 
{
  // def c-tor
  for (int i=kNFElem;i--;) fFMatrix[i] = src.fFMatrix[i];
  for (int i=kNKElem;i--;) fKMatrix[i] = src.fKMatrix[i];
  for (int i=kNRElem;i--;) fRMatrix[i] = src.fRMatrix[i];
  fResid[0]=src.fResid[0];
  fResid[1]=src.fResid[1];  
  fCovIYZ[0]=src.fCovIYZ[0];
  fCovIYZ[1]=src.fCovIYZ[1];
  fCovIYZ[2]=src.fCovIYZ[2];
  //
}

//_________________________________________________________________________
AliITSUSeed &AliITSUSeed::operator=(const AliITSUSeed& src) 
{
  // def c-tor
  if (this == &src) return *this;
  this->~AliITSUSeed();
  new(this) AliITSUSeed(src);
  return *this;
}

//_________________________________________________________________________
void AliITSUSeed::Print(Option_t* opt) const
{
  // print seed info
  int lr,cl = GetLrCluster(lr);
  printf("%cLr%d Nchild: %3d Cl:%4d Chi2Glo:%7.2f(%7.2f) Chi2Cl:%7.2f Penalty: %7.2f Mtc:%6.2f Bwd:%6.2f",IsKilled() ? '-':' ',
	 lr,GetNChildren(),cl,GetChi2Glo(),GetChi2GloNrm(),GetChi2Cl(), GetChi2Penalty(), GetChi2ITSTPC(), GetChi2ITSSA());
  printf(" |"); 
  int lrc=0;
  const AliITSUSeed *sdc = this;
  while(1) {
    if (lrc<lr) printf(".");
    else {
      sdc = sdc->GetParent(lrc);
      if (!sdc) break;
      printf("%c",sdc->GetClusterID()<0 ? '.': (sdc->IsFake() ? '-':'+')); 
    }
    lrc++;
  }
  printf("|\n");
  TString opts = opt; opts.ToLower();
  if (opts.Contains("etp")) AliExternalTrackParam::Print();
  if (opts.Contains("parent") && GetParent()) GetParent()->Print(opt);
}

//______________________________________________________________________________
void AliITSUSeed::InitFromSeed(const AliExternalTrackParam* seed)
{
  // init seed from ESD track
  TObject::Clear();
  AliExternalTrackParam::operator=(*seed);
  ResetFMatrix();
  fHitsPattern = 0;
  fClID = 0;
  fNChildren = 0;
  fChi2Glo = fChi2Cl = fChi2Penalty = 0;
  fParent = 0; //!!!
}


//______________________________________________________________________________
Float_t AliITSUSeed::GetChi2GloNrm() const
{
  int ndf = 2*GetNLayersHit() - 5;
  return (ndf>0 ? fChi2Glo/ndf : fChi2Glo) + fChi2Penalty;
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
  float chi2This  = GetChi2GloNrm();
  float chi2Other = sd->GetChi2GloNrm();

  if      (chi2This+kTol<chi2Other) return -1;
  else if (chi2This-kTol>chi2Other) return  1;
  return 0;
}

//______________________________________________________________________________
Bool_t AliITSUSeed::IsEqual(const TObject* obj)  const
{
  // compare clusters accodring to specific mode
  const AliITSUSeed* sd = (const AliITSUSeed*)obj;
  const Float_t kTol = 1e-5;
  if (IsKilled() != sd->IsKilled()) return kFALSE;
  return Abs(GetChi2GloNrm() - sd->GetChi2GloNrm())<kTol;
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

//__________________________________________________________________
Int_t AliITSUSeed::GetClusterIndex(Int_t ind) const
{
  // get ind-th cluster index
  int ncl = 0;
  const AliITSUSeed* seed = this;
  while(seed) {
    if ( seed->HasCluster() && (ncl++==ind) ) return seed->GetLrClusterID();//GetClusterID();
    seed = (AliITSUSeed*)seed->GetParent();
  }
  return -1;
  //
}

//______________________________________________________________________________
Bool_t AliITSUSeed::RotateToAlpha(Double_t alpha) 
{
  // Transform this track to the local coord. system rotated
  // by angle "alpha" (rad) with respect to the global coord. system. 
  //
  if (TMath::Abs(fP[2]) >= kAlmost1) {
     AliError(Form("Precondition is not satisfied: |sin(phi)|>1 ! %f",fP[2])); 
     return kFALSE;
  }
  //
  if      (alpha < -TMath::Pi()) alpha += 2*TMath::Pi();
  else if (alpha >= TMath::Pi()) alpha -= 2*TMath::Pi();
  //
  Double_t &fP0=fP[0];
  Double_t &fP2=fP[2];
  Double_t &fC00=fC[0];
  Double_t &fC10=fC[1];
  Double_t &fC20=fC[3];
  Double_t &fC21=fC[4];
  Double_t &fC22=fC[5];
  Double_t &fC30=fC[6];
  Double_t &fC32=fC[8];
  Double_t &fC40=fC[10];
  Double_t &fC42=fC[12];
  //
  Double_t x=fX;
  Double_t ca=TMath::Cos(alpha-fAlpha), sa=TMath::Sin(alpha-fAlpha);
  Double_t sf=fP2, cf=TMath::Sqrt((1.- fP2)*(1.+fP2)); // Improve precision
  // RS: check if rotation does no invalidate track model (cos(local_phi)>=0, i.e. particle
  // direction in local frame is along the X axis
  if ((cf*ca+sf*sa)<0) {
    AliDebug(1,Form("Rotation failed: local cos(phi) would become %.2f",cf*ca+sf*sa));
    return kFALSE;
  }
  //
  Double_t tmp=sf*ca - cf*sa;

  if (TMath::Abs(tmp) >= kAlmost1) {
     if (TMath::Abs(tmp) > 1.+ Double_t(FLT_EPSILON))  
        AliWarning(Form("Rotation failed ! %.10e",tmp));
     return kFALSE;
  }
  fAlpha = alpha;
  fX =  x*ca + fP0*sa;
  fP0= -x*sa + fP0*ca;
  fP2=  tmp;

  if (TMath::Abs(cf)<kAlmost0) {
    AliError(Form("Too small cosine value %f",cf)); 
    cf = kAlmost0;
  } 

  Double_t rr=(ca+sf/cf*sa);  

  fC00 *= (ca*ca);
  fC10 *= ca;
  fC20 *= ca*rr;
  fC21 *= rr;
  fC22 *= rr*rr;
  fC30 *= ca;
  fC32 *= rr;
  fC40 *= ca;
  fC42 *= rr;
  //
  fRMatrix[kR00] = ca;
  fRMatrix[kR22] = rr; 
  //
  CheckCovariance();

  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliITSUSeed::GetTrackingXAtXAlpha(double xOther, double alpOther, double bz, double &xdst)
{
  // calculate X and Y in the tracking frame of the track, corresponding to other X,Alpha tracking
  double ca=TMath::Cos(alpOther-fAlpha), sa=TMath::Sin(alpOther-fAlpha);
  double &y=fP[0], &sf=fP[2], cf=Sqrt((1.-sf)*(1.+sf));
  double eta = xOther - fX*ca - y*sa;
  double xi  = sf*ca - cf*sa;
  if (Abs(xi)>= kAlmost1) return kFALSE;
  double nu  = xi + GetC(bz)*eta;
  if (Abs(nu)>= kAlmost1) return kFALSE;
  xdst = xOther*ca - sa*( y*ca-fX*sa + eta*(xi+nu)/(Sqrt((1.-xi)*(1.+xi)) + Sqrt((1.-nu)*(1.+nu))) );
  return kTRUE;
}

//____________________________________________________________________
Double_t AliITSUSeed::GetPredictedChi2(Double_t p[2],Double_t cov[3]) 
{
  // Estimate the chi2 of the space point "p" with the cov. matrix "cov"
  // Store info needed for update and smoothing
  Double_t sdd = fC[0] + cov[0]; 
  Double_t sdz = fC[1] + cov[1];
  Double_t szz = fC[2] + cov[2];
  Double_t det = sdd*szz - sdz*sdz;
  if (TMath::Abs(det) < kAlmost0) return kVeryBig;
  det = 1./det;
  fCovIYZ[0] =  szz*det;
  fCovIYZ[1] = -sdz*det;
  fCovIYZ[2] =  sdd*det;
  double &dy = fResid[0] = p[0] - fP[0];
  double &dz = fResid[1] = p[1] - fP[1];
  //
  return dy*(dy*fCovIYZ[0]+dz*fCovIYZ[1]) + dz*(dy*fCovIYZ[1]+dz*(fCovIYZ[2]));
  //
}

//____________________________________________________________________
Bool_t AliITSUSeed::Update() 
{
  // Update the track parameters with the measurement stored during GetPredictedChi2
  //
  Double_t &fP0=fP[0], &fP1=fP[1], &fP2=fP[2], &fP3=fP[3], &fP4=fP[4],
    &fC00=fC[kS00],
    &fC10=fC[kS10],  &fC11=fC[kS11],  
    &fC20=fC[kS20],  &fC21=fC[kS21],  &fC22=fC[kS22],
    &fC30=fC[kS30],  &fC31=fC[kS31],  &fC32=fC[kS32],  &fC33=fC[kS33],  
    &fC40=fC[kS40],  &fC41=fC[kS41],  &fC42=fC[kS42],  &fC43=fC[kS43], &fC44=fC[kS44];
  //
  double &r00=fCovIYZ[0],&r01=fCovIYZ[1],&r11=fCovIYZ[2];
  double &dy=fResid[0], &dz=fResid[1];
  //
  // store info needed for smoothing in the fKMatrix
  double &k00 = fKMatrix[kK00] = fC00*r00+fC10*r01;
  double &k01 = fKMatrix[kK01] = fC00*r01+fC10*r11;
  double &k10 = fKMatrix[kK10] = fC10*r00+fC11*r01;
  double &k11 = fKMatrix[kK11] = fC10*r01+fC11*r11;  
  double &k20 = fKMatrix[kK20] = fC20*r00+fC21*r01;
  double &k21 = fKMatrix[kK21] = fC20*r01+fC21*r11;
  double &k30 = fKMatrix[kK30] = fC30*r00+fC31*r01;
  double &k31 = fKMatrix[kK31] = fC30*r01+fC31*r11;
  double &k40 = fKMatrix[kK40] = fC40*r00+fC41*r01;
  double &k41 = fKMatrix[kK41] = fC40*r01+fC41*r11;
  //
  Double_t sf=fP2 + k20*dy + k21*dz;
  if (TMath::Abs(sf) > kAlmost1) return kFALSE;  
  
  fP0 += k00*dy + k01*dz;
  fP1 += k10*dy + k11*dz;
  fP2  = sf;
  fP3 += k30*dy + k31*dz;
  fP4 += k40*dy + k41*dz;
  //
  Double_t c01=fC10, c02=fC20, c03=fC30, c04=fC40;
  Double_t c12=fC21, c13=fC31, c14=fC41;

  fC00-=k00*fC00+k01*fC10; fC10-=k00*c01+k01*fC11;
  fC20-=k00*c02+k01*c12;   fC30-=k00*c03+k01*c13;
  fC40-=k00*c04+k01*c14; 

  fC11-=k10*c01+k11*fC11;
  fC21-=k10*c02+k11*c12;   fC31-=k10*c03+k11*c13;
  fC41-=k10*c04+k11*c14; 

  fC22-=k20*c02+k21*c12;   fC32-=k20*c03+k21*c13;
  fC42-=k20*c04+k21*c14; 

  fC33-=k30*c03+k31*c13;
  fC43-=k30*c04+k31*c14; 
  
  fC44-=k40*c04+k41*c14; 
  //
  CheckCovariance();
  //
  return kTRUE;
}


//____________________________________________________________________
Bool_t AliITSUSeed::Smooth(Double_t vecL[5],Double_t matL[15]) 
{
  // Prepare MBF smoothing auxiliary params for smoothing at prev. point:
  // \hat{l_N} = 0
  // \hat{L_N} = 0
  // \tilde{l_j} = -H_j^T N_{j}^{-1} z_j + B_{j}^T \hat{l_j}
  // \tilde{L_j} =  H_j^T N_{j}^{-1} H_j + B_j^T \hat{L_j} B_j
  // \hat{l_j} = F_j^T \tilde{l_{j+1}}
  // \hat{L_j} = F_j^T \tilde{L_{j+1}} F_j
  //
  // P_{j/N} = P_{j/j} - P_{j/j} \hat{L_j} P_{j/j}
  // \hat{x_{j/N}} = \hat{x_{j/j}} - P_{j/j} \hat{l_j}
  //
  // N^-1 = fCovIYZ
  // z = fResid
  // B = I - K H
  // H = {{1,0,0,0,0},{0,1,0,0,0}}
  // 
  // calc. \tilde{l_j} 
  //
  if (GetClusterID()<0) return kTRUE;
  //

  double 
    &k00=fKMatrix[kK00],&k01=fKMatrix[kK01],
    &k10=fKMatrix[kK10],&k11=fKMatrix[kK11],
    &k20=fKMatrix[kK20],&k21=fKMatrix[kK21],
    &k30=fKMatrix[kK30],&k31=fKMatrix[kK31],
    &k40=fKMatrix[kK40],&k41=fKMatrix[kK41];
  double 
    &l00=matL[kS00],
    &l10=matL[kS10], &l11=matL[kS11],  
    &l20=matL[kS20], &l21=matL[kS21], &l22=matL[kS22],
    &l30=matL[kS30], &l31=matL[kS31], &l32=matL[kS32], &l33=matL[kS33],  
    &l40=matL[kS40], &l41=matL[kS41], &l42=matL[kS42], &l43=matL[kS43], &l44=matL[kS44];
  //
  // calculate correction
  double corrVec[5]={0},corrMat[15]={0};
  corrVec[0] = fC[kS00]*vecL[0] + fC[kS10]*vecL[1] + fC[kS20]*vecL[2] + fC[kS30]*vecL[3] + fC[kS40]*vecL[4]; 
  corrVec[1] = fC[kS10]*vecL[0] + fC[kS11]*vecL[1] + fC[kS21]*vecL[2] + fC[kS31]*vecL[3] + fC[kS41]*vecL[4]; 
  corrVec[2] = fC[kS20]*vecL[0] + fC[kS21]*vecL[1] + fC[kS22]*vecL[2] + fC[kS32]*vecL[3] + fC[kS42]*vecL[4]; 
  corrVec[3] = fC[kS30]*vecL[0] + fC[kS31]*vecL[1] + fC[kS32]*vecL[2] + fC[kS33]*vecL[3] + fC[kS43]*vecL[4]; 
  corrVec[4] = fC[kS40]*vecL[0] + fC[kS41]*vecL[1] + fC[kS42]*vecL[2] + fC[kS43]*vecL[3] + fC[kS44]*vecL[4]; 
  //
  double *crm = ProdABA(fC,matL);
  for (int i=0;i<15;i++) corrMat[i] = crm[i];

  double vcL0 = vecL[0], vcL1 = vecL[1];
  vecL[0] -= k00*vcL0+k10*vcL1+k20*vecL[2]+k30*vecL[3]+k40*vecL[4] + fCovIYZ[0]*fResid[0] + fCovIYZ[1]*fResid[1];
  vecL[1] -= k01*vcL0+k11*vcL1+k21*vecL[2]+k31*vecL[3]+k41*vecL[4] + fCovIYZ[1]*fResid[0] + fCovIYZ[2]*fResid[1];

  /*
  double vcL0 = vecL[0], vcL1 = vecL[1];
  vecL[0] -= k00*vcL0+k10*vcL1+fKMatrix[kK20]*vecL[2]+k30*vecL[3]+k40*vecL[4] + fCovIYZ[0]*fResid[0] + fCovIYZ[1]*fResid[1];
  vecL[1] -= k01*vcL0+fKMatrix[kK11]*vcL1+k21*vecL[2]+k31*vecL[3]+k41*vecL[4] + fCovIYZ[1]*fResid[0] + fCovIYZ[2]*fResid[1];
  vecL[3] += fFMatrix[kF13]*vecL[1]; 
  vecL[4]  = fFMatrix[kF04]*vecL[0] + fFMatrix[kF14]*vecL[1] + fFMatrix[kF24]*vecL[2] + fFMatrix[kF44]*vecL[4];
  vecL[2] += fFMatrix[kF02]*vecL[0] + fFMatrix[kF12]*vecL[1];
  //
  */
  // and \hat{l_j} in one go

  // L = H^T * sg * H + (I-KH)^T * L * (I - KH)
  double v00 =  k00*l00+k10*l10+k20*l20+k30*l30+k40*l40;
  double v10 =  k00*l10+k10*l11+k20*l21+k30*l31+k40*l41;
  double v20 =  k00*l20+k10*l21+k20*l22+k30*l32+k40*l42;
  double v30 =  k00*l30+k10*l31+k20*l32+k30*l33+k40*l43;
  double v40 =  k00*l40+k10*l41+k20*l42+k30*l43+k40*l44;
  //
  double v01 =  k01*l00+k11*l10+k21*l20+k31*l30+k41*l40;
  double v11 =  k01*l10+k11*l11+k21*l21+k31*l31+k41*l41;
  double v21 =  k01*l20+k11*l21+k21*l22+k31*l32+k41*l42;
  double v31 =  k01*l30+k11*l31+k21*l32+k31*l33+k41*l43;
  double v41 =  k01*l40+k11*l41+k21*l42+k31*l43+k41*l44;
  //
  // (H^T * K^T * L * K * H) - (L * K * H) - (H^T * K^T * L) + (H^T*N^-1*H)
  l00 += k00*v00 + k10*v10 + k20*v20 + k30*v30 + k40*v40 - v00 - v00 + fCovIYZ[0];
  l10 += k01*v00 + k11*v10 + k21*v20 + k31*v30 + k41*v40 - v01 - v10 + fCovIYZ[1];
  l11 += k01*v01 + k11*v11 + k21*v21 + k31*v31 + k41*v41 - v11 - v11 + fCovIYZ[2];
  //
  l20 -= v20;
  l21 -= v21;
  l30 -= v30;
  l31 -= v31;
  l40 -= v40;
  l41 -= v41;
  //
  printf("CorrMt:\n");
  printf("%+e\n%+e %+e\n%+e %+e %+e\n%+e %+e %+e %+e\n%+e %+e %+e %+e %+e\n",
	 corrMat[kS00],corrMat[kS10],corrMat[kS11],corrMat[kS20],corrMat[kS21],corrMat[kS22],
	 corrMat[kS30],corrMat[kS31],corrMat[kS32],corrMat[kS33],
	 corrMat[kS40],corrMat[kS41],corrMat[kS42],corrMat[kS43],corrMat[kS44]);
  
  printf("SMcorr: %+e %+e %+e %+e %+e\n",corrVec[0],corrVec[1],corrVec[2],corrVec[3],corrVec[4]);

  printf("State : "); this->AliExternalTrackParam::Print("");
  //
  printf("\nBefore transport back (RotElems: %+e %+e)\n",fRMatrix[kR00],fRMatrix[kR22]);
  printf("Res: %+e %+e | Err: %+e %+e %+e\n",fResid[0],fResid[1],fCovIYZ[0],fCovIYZ[1],fCovIYZ[2]);
  printf("Lr%d VecL: ",GetLayerID()); for (int i=0;i<5;i++) printf("%+e ",vecL[i]); printf("\n");
  //
  printf("%+e\n%+e %+e\n%+e %+e %+e\n%+e %+e %+e %+e\n%+e %+e %+e %+e %+e\n",
	 matL[kS00],matL[kS10],matL[kS11],matL[kS20],matL[kS21],matL[kS22],
	 matL[kS30],matL[kS31],matL[kS32],matL[kS33],matL[kS40],matL[kS41],matL[kS42],matL[kS43],matL[kS44]);
  //
  printf("F: "); for (int i=0;i<kNFElem;i++) printf("%+e ",fFMatrix[i]); printf("\n");  
  printf("K: "); for (int i=0;i<kNKElem;i++) printf("%+e ",fKMatrix[i]); printf("\n");  
  //
  // apply rotation matrix (diagonal)
  vecL[0] *= fRMatrix[kR00];
  vecL[2] *= fRMatrix[kR22];
  //
  l00 *= fRMatrix[kR00]*fRMatrix[kR00];
  l10 *= fRMatrix[kR00];
  l20 *= fRMatrix[kR22]*fRMatrix[kR00];
  l21 *= fRMatrix[kR22];
  l22 *= fRMatrix[kR22]*fRMatrix[kR22];
  l30 *= fRMatrix[kR00];
  l32 *= fRMatrix[kR22];
  l40 *= fRMatrix[kR00];
  l42 *= fRMatrix[kR22];
  //
  // Apply translation matrix F^T. Note, that fFMatrix keeps non-trivial elems of F-1 = f, except the e-loss coeff f44
  // We need F^T*L* F = L + (L*f) + (L*f)^T + f^T * (L*f)
  //
  double 
    &f02=fFMatrix[kF02],&f04=fFMatrix[kF04],
    &f12=fFMatrix[kF12],&f13=fFMatrix[kF13],&f14=fFMatrix[kF14],
    &f24=fFMatrix[kF24],
    f44 =fFMatrix[kF44];
  //
  vecL[4]  = f04*vecL[0]+f14*vecL[1]+f24*vecL[2]+f44*vecL[4];
  vecL[3] += f13*vecL[1];
  vecL[2] += f02*vecL[0]+f12*vecL[1];
  //
  f44 -= 1.0; // !!!!!
  //
  //b = L*f
  Double_t b02=l00*f02+l10*f12, b03=l10*f13, b04=l00*f04+l10*f14+l20*f24+l40*f44;
  Double_t b12=l10*f02+l11*f12, b13=l11*f13, b14=l10*f04+l11*f14+l21*f24+l41*f44;
  Double_t b22=l20*f02+l21*f12, b23=l21*f13, b24=l20*f04+l21*f14+l22*f24+l42*f44;
  Double_t b32=l30*f02+l31*f12, b33=l31*f13, b34=l30*f04+l31*f14+l32*f24+l43*f44;
  Double_t b42=l40*f02+l41*f12, b43=l41*f13, b44=l40*f04+l41*f14+l42*f24+l44*f44;
  //
  //a = f^T * b = f^T * L * f, profit from symmetry
  Double_t a22=f02*b02+f12*b12, a33=f13*b13, a44=f04*b04+f14*b14+f24*b24+f44*b44,
    a32=f13*b12, //= a23=f02*b03+f12*b13, 
    a42=f02*b04+f12*b14, //f04*b02+f14*b12+f24*b22+f44*b42 = a24
    a43=f13*b14;         //f04*b03+f14*b13+f24*b23+f44*b43 = a34
  //
  // F^T*L* F = L + (b + b^T + a)
  l44 += b44 + b44 + a44;
  l43 += b43 + b34 + a43;
  l42 += b42 + b24 + a42;
  l41 += b14;
  l40 += b04;
  l33 += b33 + b33 + a33;
  l32 += b32 + b23 + a32;
  l31 += b13;
  l30 += b03;
  l22 += b22 + b23 + a22;
  l21 += b12;
  l20 += b02;
  //
  printf("After transport back\n");
  printf("Lr%d VecL: ",GetLayerID()); for (int i=0;i<5;i++) printf("%+e ",vecL[i]); printf("\n");
  //
  printf("%+e\n%+e %+e\n%+e %+e %+e\n%+e %+e %+e %+e\n%+e %+e %+e %+e %+e\n",
	 matL[kS00],matL[kS10],matL[kS11],matL[kS20],matL[kS21],matL[kS22],
	 matL[kS30],matL[kS31],matL[kS32],matL[kS33],matL[kS40],matL[kS41],matL[kS42],matL[kS43],matL[kS44]);

  return kTRUE;
}

//____________________________________________________________________
Double_t* AliITSUSeed::ProdABA(const double a[15],const double b[15]) const
{
  // product of symmetric matrices A*B*A
  //
  const Short_t knd[5][5] = {
    {kS00,kS10,kS20,kS30,kS40},
    {kS10,kS11,kS21,kS31,kS41},
    {kS20,kS21,kS22,kS32,kS42},
    {kS30,kS31,kS32,kS33,kS43},
    {kS40,kS41,kS42,kS43,kS44}
  };
  //
  static double aba[15];
  // 1) ba = B*A
  double ba[5][5];
  for (int i=5;i--;) for (int j=5;j--;) {
      ba[i][j] = 0;
      for (int k=5;k--;) ba[i][j] += b[knd[i][k]]*a[knd[k][j]];
    }
  //
  // 2) A * ba, lower triangle only
  for (int i=5;i--;) for (int j=i+1;j--;) {
      aba[knd[i][j]] = 0;
      for (int k=5;k--;) aba[knd[i][j]] += a[knd[i][k]]*ba[k][j];
    }
  //
  return &aba[0];
}

//____________________________________________________________________
Bool_t AliITSUSeed::ContainsFake() const
{
  // check if the full branch containes a fake cluster
  const AliITSUSeed* seed = this;
  while(seed) {
    if ( seed->IsFake() ) return kTRUE;
    seed = (AliITSUSeed*)seed->GetParent();
  }  
  return kFALSE;
}

//__________________________________________________________________
Int_t AliITSUSeed::FetchClusterInfo(Int_t *clIDarr) const
{
  // fill cl.id's in the array. The clusters of layer L will be set at slots
  // clID[2L] (and clID[2L+1] if there is an extra cluster).
  Int_t lr,ncl=0;
  const AliITSUSeed* seed = this;
  do {
    int clID = seed->GetLrCluster(lr);
    if (clID>=0) {
      lr<<=1;
      clIDarr[ clIDarr[lr]<0 ? lr : lr+1 ] = clID;
      ncl++;
    }
  } while ((seed=(AliITSUSeed*)seed->GetParent()));
  return ncl;
}

/*
//____________________________________________________________________
Bool_t AliITSUSeed::Smooth(Double_t vecL[5],Double_t matL[15]) 
{
  // Prepare MBF smoothing auxiliary params for smoothing at prev. point:
  // \hat{l_N} = 0
  // \hat{L_N} = 0
  // \tilde{l_j} = -H_j^T N_{j}^{-1} z_j + B_{j}^T \hat{l_j}
  // \tilde{L_j} =  H_j^T N_{j}^{-1} H_j + B_j^T \hat{L_j} B_j
  // \hat{l_j} = F_j^T \tilde{l_{j+1}}
  // \hat{L_j} = F_j^T \tilde{L_{j+1}} F_j
  //
  // P_{j/N} = P_{j/j} - P_{j/j} \hat{L_j} P_{j/j}
  // \hat{x_{j/N}} = \hat{x_{j/j}} - P_{j/j} \hat{l_j}
  //
  // N^-1 = fCovIYZ
  // z = fResid
  // B = I - K H
  // H = {{1,0,0,0,0},{0,1,0,0,0}}
  // 
  // calc. \tilde{l_j} and \hat{l_j} in one go
  //
  if (GetClusterID()<0) return kTRUE;
  //
  double 
    &k00=fKMatrix[kK00],&k01=fKMatrix[kK01],
    &k10=fKMatrix[kK10],&k11=fKMatrix[kK11],
    &k20=fKMatrix[kK20],&k21=fKMatrix[kK21],
    &k30=fKMatrix[kK30],&k31=fKMatrix[kK31],
    &k40=fKMatrix[kK40],&k41=fKMatrix[kK41];
  double 
    &matL00=matL[kS00],
    &matL10=matL[kS01],  &matL11=matL[kS11],  
    &matL20=matL[kS20],  &matL21=matL[kS21],  &matL22=matL[kS22],
    &matL30=matL[kS30],  &matL31=matL[kS31],  &matL32=matL[kS32],  &matL33=matL[kS33],  
    &matL40=matL[kS40],  &matL41=matL[kS41],  &matL42=matL[kS42],  &matL43=matL[kS43], &matL44=matL[kS44];
  //
  double vcL0 = vecL[0], vcL1 = vecL[1];
  vecL[0] -= k00*vcL0+k10*vcL1+fKMatrix[kK20]*vecL[2]+k30*vecL[3]+k40*vecL[4] + fCovIYZ[0]*fResid[0] + fCovIYZ[1]*fResid[1];
  vecL[1] -= k01*vcL0+fKMatrix[kK11]*vcL1+k21*vecL[2]+k31*vecL[3]+k41*vecL[4] + fCovIYZ[1]*fResid[0] + fCovIYZ[2]*fResid[1];
  vecL[3] += fFMatrix[kF13]*vecL[1]; 
  vecL[4] =  fFMatrix[kF04]*vecL[0] + fFMatrix[kF14]*vecL[1] + fFMatrix[kF24]*vecL[2] + fFMatrix[kF44]*vecL[4];
  vecL[2] += fFMatrix[kF02]*vecL[0] + fFMatrix[kF12]*vecL[1];
  //

  // L = H^T * sg * H + (I-KH)^T * L * (I - KH)
  double v00 =  k00*matL00+k10*matL10+k20*matL20+k30*matL30+k40*matL40;
  double v10 =  k00*matL10+k10*matL11+k20*matL21+k30*matL31+k40*matL41;
  double v20 =  k00*matL20+k10*matL12+k20*matL22+k30*matL32+k40*matL42;
  double v30 =  k00*matL30+k10*matL13+k20*matL23+k30*matL33+k40*matL43;
  double v40 =  k00*matL40+k10*matL14+k20*matL24+k30*matL34+k40*matL44;
  //
  double v01 =  k01*matL00+k11*matL10+k21*matL20+k31*matL30+k41*matL40;
  double v11 =  k01*matL01+k11*matL11+k21*matL21+k31*matL31+k41*matL41;
  double v21 =  k01*matL02+k11*matL12+k21*matL22+k31*matL32+k41*matL42;
  double v31 =  k01*matL03+k11*matL13+k21*matL23+k31*matL33+k41*matL43;
  double v41 =  k01*matL04+k11*matL14+k21*matL24+k31*matL34+k41*matL44;
  //
  double t00 =  k00*matL00+k10*matL01+k20*matL02+k30*matL03+k40*matL04;
  double t10 =  k00*matL10+k10*matL11+k20*matL12+k30*matL13+k40*matL14;
  double t20 =  k00*matL20+k10*matL21+k20*matL22+k30*matL23+k40*matL24;
  double t30 =  k00*matL30+k10*matL31+k20*matL32+k30*matL33+k40*matL34;
  double t40 =  k00*matL40+k10*matL41+k20*matL42+k30*matL43+k40*matL44;
  //
  double t01 =  k01*matL00+k11*matL01+k21*matL02+k31*matL03+k41*matL04;
  double t11 =  k01*matL10+k11*matL11+k21*matL12+k31*matL13+k41*matL14;
  double t21 =  k01*matL20+k11*matL21+k21*matL22+k31*matL23+k41*matL24;
  double t31 =  k01*matL30+k11*matL31+k21*matL32+k31*matL33+k41*matL34;
  double t41 =  k01*matL40+k11*matL41+k21*matL42+k31*matL43+k41*matL44;
  //
  // (H^T * K^T * L * K * H) - (L * K * H) - (H^T * K^T * L) + (H^T*N^-1*H)
  matL00 += k00*v00+k10*v10+k20*v20*k30*v30+k40*v40 - t00 - v00 + fCovIYZ[0];
  matL01 += k01*v00+k11*v10+k21*v20*k31*v30+k41*v40 - t01 - v10 + fCovIYZ[1];
  matL10 += k00*v01+k10*v11+k20*v21*k30*v31+k40*v41 - t10 - v01 + fCovIYZ[1];
  matL11 += k01*v01+k11*v11+k21*v21*k31*v31+k41*v41 - t11 - v11 + fCovIYZ[2];
  //
  matL20 -= t20;
  matL21 -= t21;
  matL30 -= t30;
  matL31 -= t31;
  matL40 -= t40;
  matL41 -= t41;
  //
  matL02 -= v20;
  matL03 -= v30;
  matL04 -= v40;
  matL12 -= v21;
  matL13 -= v31;
  matL14 -= v41;
  //
  printf("Lr%d VecL: ",GetLayerID()); for (int i=0;i<5;i++) printf("%+e ",vecL[i]); printf("\n");
  printf("F: "); for (int i=0;i<kNFElem;i++) printf("%+e ",fFMatrix[i]); printf("\n");  
  printf("K: "); for (int i=0;i<kNKElem;i++) printf("%+e ",fKMatrix[i]); printf("\n");  
  //
  for (int j=0;j<5;j++) {
    for (int i=0;i<5;i++) printf("%+e ",matL[j][i]); printf("\n");  
  }
  //
  return kTRUE;
}

 */

