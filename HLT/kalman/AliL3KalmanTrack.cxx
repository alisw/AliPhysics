#include "TMath.h"
#include "AliL3KalmanTrack.h"
#include "AliL3SpacePointData.h"
#include "AliL3Logging.h"
ClassImp(AliL3KalmanTrack)

// Class for kalman tracks

AliL3KalmanTrack::AliL3KalmanTrack()
{
  // Constructor
}

AliL3KalmanTrack::~AliL3KalmanTrack()
{
  // Destructor
}

void AliL3KalmanTrack::Init()
{
  // Set state vector and covariance matrix to 0
  fP0 = 0; fP1 = 0; fP2 = 0; fP3 = 0; fP4 = 0;
  fC00 = 0; 
  fC10 = 0; fC11 = 0;
  fC20 = 0; fC21 = 0; fC22 = 0;
  fC30 = 0; fC31 = 0; fC32 = 0; fC33 = 0;
  fC40 = 0; fC41 = 0; fC42 = 0; fC43 = 0; fC44 = 0;

  fChisq = 0;

  fX = 0; 
}

Int_t AliL3KalmanTrack::MakeTrackSeed(AliL3SpacePointData *points1, UInt_t pos1, AliL3SpacePointData *points2, UInt_t pos2, AliL3SpacePointData *points3, UInt_t pos3)
{
  // Make track seed based on three clusters 
  fX = points1[pos1].fX;

  fP0 = points1[pos1].fY; 
  fP1 = points1[pos1].fZ; 

  Float_t X2 = points2[pos2].fX;
  Float_t Y2 = points2[pos2].fY;
  Float_t Z2 = points2[pos2].fZ;
  //Float_t alpha = TMath::ATan((fX2-fX)/(fY2-fP0));
  //Float_t X2 = XX2*TMath::Cos(alpha) + YY2*TMath::Sin(alpha); 

  Float_t X3 = points3[pos3].fX;
  Float_t Y3 = points3[pos3].fY;
  Float_t Z3 = points3[pos3].fZ;

  Float_t ZZ = fP1 - ((fP1 - Z3)/(fX-X3))*(fX-X2);
  if (TMath::Abs(ZZ - Z2) > 10) return 0; //What's this?? (fP1 - Z3)/(fX-X3)*(fX-X2) is an angle

  // It may make no difference. Check on a big event??.
  if ((X2-fX)*(0-Y2)-(0-X2)*(Y2-fP0) == 0) return 0; //Straight seed

  // Initial approximation of the state vector
  fP2 = f2(fX,fP0,X2,Y2,X3,Y3);
  fP3 = f3(fX,fP0,X2,Y2,fP1,Z2);
  fP4 = f4(fX,fP0,X2,Y2,X3,Y3);

  // Initial approximation of the covariance matrix
  Float_t sy1=points1[pos1].fSigmaY2;
  Float_t sz1=points1[pos1].fSigmaZ2;
  Float_t sy2=points2[pos2].fSigmaY2;
  Float_t sz2=points2[pos2].fSigmaZ2;
  Float_t sy3=25000*fP4*fP4+0.1;
  Float_t sy=0.1;
  Float_t sz=0.1;
  
  Float_t f40=(f4(fX,fP0+sy,X2,Y2,X3,Y3)-fP4)/sy;
  Float_t f42=(f4(fX,fP0,X2,Y2+sy,X3,Y3)-fP4)/sy;
  Float_t f43=(f4(fX,fP0,X2,Y2,X3,Y3+sy)-fP4)/sy;
  Float_t f20=(f2(fX,fP0+sy,X2,Y2,X3,Y3)-fP2)/sy;
  Float_t f22=(f2(fX,fP0,X2,Y2+sy,X3,Y3)-fP2)/sy;
  Float_t f23=(f2(fX,fP0,X2,Y2,X3,Y3+sy)-fP2)/sy;
  Float_t f30=(f3(fX,fP0+sy,X2,Y2,fP1,Z2)-fP3)/sy;
  Float_t f31=(f3(fX,fP0,X2,Y2,fP1+sz,Z2)-fP3)/sz;
  Float_t f32=(f3(fX,fP0,X2,Y2+sy,fP1,Z2)-fP3)/sy;
  Float_t f34=(f3(fX,fP0,X2,Y2,fP1,Z2+sz)-fP3)/sz;
  
  fC00 = sy1;
  fC10 = 0;  
  fC11 = sz1;
  fC20 = f20*sy1;  
  fC21 = 0;  
  fC22 = f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
  fC30 = f30*sy1;  
  fC31 = f31*sz1;  
  fC32 = f30*sy1*f20+f32*sy2*f22;
  fC33 = f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
  fC40 = f40*sy1; 
  fC41 = 0; 
  fC42 = f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
  fC43 = f30*sy1*f40+f32*sy2*f42; 
  fC44 = f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
  
  return 1;
}

Int_t AliL3KalmanTrack::Propagate(AliL3SpacePointData *points, UInt_t pos)
{
  // Propagetes track to the plane of the next found cluster

  Float_t Xold = fX; // X position for previous space point
  Float_t Xnew = points[pos].fX; // X position of current space point
  Float_t dx = Xnew - Xold;
  Float_t Yold = fP0; // Y position of old point
  Float_t Zold = fP1; // Z position of old point

  if (TMath::Abs(fP4*Xnew - fP2) >= 0.9) // What's this??
    {
      return 0;
    }

  // R must be approx. of the radius of the circle based on
  // the old and new spacepoint. What then is Cod and Cnew??
  // C has something to do with curvature at least.
  Float_t Cold = fP4*Xold - fP2; 
  Float_t Rold = TMath::Sqrt(1 - Cold*Cold);
  Float_t Cnew = fP4*Xnew - fP2; 
  Float_t Rnew = TMath::Sqrt(1 - Cnew*Cnew);

  // Prediction of the y- and z- coordinate in the next plane
  fP0 += dx*(Cold+Cnew)/(Rold+Rnew);
  fP1 += dx*(Cold+Cnew)/(Cold*Rnew + Cnew*Rold)*fP3; 

  // f = F - 1 //What is this??
  // Must be the f-matrix for the prediction, as in eq 1 in ALICE Kalman paper
  Float_t RR = Rold + Rnew;
  Float_t CC = Cold + Cnew;
  Float_t XX = Xold + Xnew;

  Float_t f02 = -dx*(2*RR + CC*(Cold/Rold + Cnew/Rnew))/(RR*RR);
  Float_t f04 = dx*(RR*XX + CC*(Cold*Xold/Rold + Cnew*Xnew/Rnew))/(RR*RR);
  Float_t CR = Cold*Rnew + Cnew*Rold;
  Float_t f12 = -dx*fP3*(2*CR + CC*(Cnew*(Cold/Rold)-Rold + Cold*(Cnew/Rnew)-Rnew))/(CR*CR);
  Float_t f13 = dx*CC/CR;
  Float_t f14 = dx*fP3*(CR*XX-CC*(Rold*Xnew-Cnew*Cold*Xold/Rold + Rnew*Xold-Cold*Cnew*Xnew/Rnew))/(CR*CR);

  // b = C*ft // This? 
  Float_t b00=f02*fC20 + f04*fC40;
  Float_t b01=f12*fC20 + f14*fC40 + f13*fC30;
  Float_t b10=f02*fC21 + f04*fC41;
  Float_t b11=f12*fC21 + f14*fC41 + f13*fC31;
  Float_t b20=f02*fC22 + f04*fC42;
  Float_t b21=f12*fC22 + f14*fC42 + f13*fC32;
  Float_t b30=f02*fC32 + f04*fC43;
  Float_t b31=f12*fC32 + f14*fC43 + f13*fC33;
  Float_t b40=f02*fC42 + f04*fC44;
  Float_t b41=f12*fC42 + f14*fC44 + f13*fC43;

  //a = f*b = f*C*ft
  Float_t a00 = f02*b20 + f04*b40;
  Float_t a01 = f02*b21 + f04*b41;
  Float_t a11 = f12*b21 + f14*b41+f13*b31;

  //F*C*Ft = C + (a + b + bt) /This is the covariance matrix, the samll t 
  // means transform. Then F must be df/dx 
  fC00 += a00 + 2*b00;
  fC10 += a01 + b01 + b10;
  fC20 += b20;
  fC30 += b30;
  fC40 += b40;
  fC11 += a11 + 2*b11;
  fC21 += b21;
  fC31 += b31;
  fC41 += b41;

  // Multiple scattering (from AliTPCtrack::PropagateTo)
  Float_t d = TMath::Sqrt((Xold-Xnew)*(Xold-Xnew)+(Yold-fP0)*(Yold-fP0)+(Zold-fP1)*(Zold-fP1));
  Float_t Pt = (1e-9*TMath::Abs(fP4)/fP4 + fP4 * (1000/0.299792458/4));
  if (TMath::Abs(Pt) > 100) {
    return 0;
  }
  if (TMath::Abs(Pt) < 0.01) return 0;

  Float_t p2 = (1+fP3*fP3)/(Pt*Pt);
  Float_t beta2 = p2/(p2 +0.14*0.14);
  Float_t theta2 = 1.0259e-6*10*10/20/(beta2*p2)*d*0.9e-3;

  Float_t ey=fP4*Xnew - fP2;
  Float_t ez=fP3;
  Float_t xz=fP4*ez;
  Float_t zz1=ez*ez+1;
  Float_t xy=fP2+ey;

  fC22 += (2*ey*ez*ez*fP2+1-ey*ey+ez*ez+fP2*fP2*ez*ez)*theta2;
  fC32 += ez*zz1*xy*theta2;
  fC33 += zz1*zz1*theta2;
  fC42 += xz*ez*xy*theta2;
  fC43 += xz*zz1*theta2;
  fC44 += xz*xz*theta2;
  if (TMath::Abs(beta2) >= 1) printf("%f %f\n",beta2,Pt);

  // Energy loss
  Float_t dE = 0.153e-3/beta2*(log(5940*beta2/(1-beta2))-beta2)*d*0.9e-3;
  if (Xold < Xnew) dE = -dE;
  CC = fP4;
  fP4 *= (1 - TMath::Sqrt(p2+0.14*0.14)/p2*dE);
  fP2 += Xnew*(fP4-CC);

  // Update the track parameters with the measured values of the new point
  UpdateTrack(points, pos);

  return 1;
}

Int_t AliL3KalmanTrack::UpdateTrack(AliL3SpacePointData *points, UInt_t pos)
{
  // Update the track parameters with the measured values
  fX = points[pos].fX; 

  // The errors from the measurement of the spacepoint
  Float_t sigmaY2 = points[pos].fSigmaY2;
  Float_t sigmaZ2 = points[pos].fSigmaZ2;
  Float_t sigmaYZ = 0;

  sigmaY2 += fC00;
  sigmaZ2 += fC11;
  sigmaYZ += fC10;

  Float_t det = sigmaY2*sigmaZ2 - sigmaYZ*sigmaYZ;
  Float_t tmp = sigmaY2;
  sigmaY2 = sigmaZ2/det;
  sigmaZ2 = tmp/det;
  sigmaYZ = -sigmaYZ/det;

  // What's this?? Must be the Kalman gain matrix
  Float_t k00 = fC00*sigmaY2 + fC10*sigmaYZ;
  Float_t k01 = fC00*sigmaYZ + fC10*sigmaZ2;
  Float_t k10 = fC10*sigmaY2 + fC11*sigmaYZ;
  Float_t k11 = fC10*sigmaYZ + fC11*sigmaZ2;
  Float_t k20 = fC20*sigmaY2 + fC21*sigmaYZ;
  Float_t k21 = fC20*sigmaYZ + fC21*sigmaZ2;
  Float_t k30 = fC30*sigmaY2 + fC31*sigmaYZ;
  Float_t k31 = fC30*sigmaYZ + fC31*sigmaZ2;
  Float_t k40 = fC40*sigmaY2 + fC41*sigmaYZ;
  Float_t k41 = fC40*sigmaYZ + fC41*sigmaZ2;

  // Deviation between the predicted and measured values of y and z 
  Float_t dy = points[pos].fY-fP0;
  Float_t dz = points[pos].fZ-fP1;
  //printf("%f,%f\n",fC00,fC10);

  // Prediction of fP2 and fP4
  Float_t cur = fP4 + k40*dy + k41*dz; 
  Float_t eta = fP2 + k20*dy + k21*dz;

  if (TMath::Abs(cur*fX - eta) >= 0.9) return 0;

  // Filtered state vector
  fP0 += k00*dy + k01*dz;
  fP1 += k10*dy + k11*dz;
  fP2 = eta;
  fP3 += k30*dy + k31*dz;
  fP4 = cur;

  Float_t c10 = fC10;
  Float_t c20 = fC20;
  Float_t c30 = fC30;
  Float_t c40 = fC40;
  Float_t c21 = fC21;
  Float_t c31 = fC31;
  Float_t c41 = fC41;

  // Filtered covariance matrix
  fC00 -= k00*fC00 + k01*fC10;
  fC10 -= k00*c10 + k01*fC11;
  fC11 -= k10*c10 + k11*fC11;
  fC20 -= k00*c20 + k01*c21;
  fC21 -= k10*c20 + k11*c21;
  fC22 -= k20*c20 + k21*c21;
  fC30 -= k00*c30 + k01*c31;
  fC31 -= k10*c30 + k11*c31;
  fC32 -= k20*c30 + k21*c31;
  fC33 -= k30*c30 + k31*c31;
  fC40 -= k00*c40 + k01*c41;
  fC41 -= k10*c40 + k11*c41;
  fC42 -= k20*c40 + k21*c41;
  fC43 -= k40*c30 + k41*c31;
  fC44 -= k40*c40 + k41*c41;

  sigmaY2 = sigmaY2*det;
  sigmaZ2 = sigmaZ2*det;
  sigmaYZ = sigmaYZ*det;

  // Calculate increase of chisquare
  fChisq = GetChisq() + (dy*sigmaY2*dy + 2*sigmaYZ*dy*dz + dz*sigmaZ2*dz) / (sigmaY2*sigmaZ2 - sigmaYZ*sigmaYZ);
  // Must at some point make an cut on chisq. Here?

  return 1;
} 

Float_t AliL3KalmanTrack::f2(Float_t x1,Float_t y1,
			     Float_t x2,Float_t y2,
			     Float_t x3,Float_t y3)
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature times center of curvature
  //-----------------------------------------------------------------
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);

  return -a/(d*y1-b)*xr/sqrt(xr*xr+yr*yr);
}

Float_t AliL3KalmanTrack::f3(Float_t x1,Float_t y1,
			     Float_t x2,Float_t y2,
			     Float_t z1,Float_t z2)
{
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

Float_t AliL3KalmanTrack::f4(Float_t x1,Float_t y1,
			     Float_t x2,Float_t y2,
			     Float_t x3,Float_t y3)
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);

  return -xr*yr/sqrt(xr*xr+yr*yr);
}

void AliL3KalmanTrack::Set(AliL3KalmanTrack *track)
{

  AliL3KalmanTrack *tpt = (AliL3KalmanTrack*)track;
  SetX0(tpt->GetX0());
  SetX1(tpt->GetX1());
  SetX2(tpt->GetX2());
  SetX3(tpt->GetX3());
  SetX4(tpt->GetX4());

  SetC0(tpt->GetC0());
  SetC1(tpt->GetC1());
  SetC2(tpt->GetC2());
  SetC3(tpt->GetC3());
  SetC4(tpt->GetC4());
  SetC5(tpt->GetC5());
  SetC6(tpt->GetC6());
  SetC7(tpt->GetC7());
  SetC8(tpt->GetC8());
  SetC9(tpt->GetC9());
  SetC10(tpt->GetC10());
  SetC11(tpt->GetC11());
  SetC12(tpt->GetC12());
  SetC13(tpt->GetC13());
  SetC14(tpt->GetC14());

  SetNHits(tpt->GetNHits());
}
