#include "TMath.h"
#include "AliL3KalmanTrack.h"
#include "AliL3SpacePointData.h"
#include "AliL3Logging.h"
#include "AliL3Transform.h"
#include "AliL3StandardIncludes.h"

// includes for offline comparison, will be removed
#include "AliTPCtrack.h"
// includes for offline comparison, will be removed

#include "Riostream.h" 


ClassImp(AliL3KalmanTrack)

// Class for kalman tracks

AliL3KalmanTrack::AliL3KalmanTrack()
{
  fX = 0;

  fMaxChi2 = 12;
  // Constructor
}

AliL3KalmanTrack::~AliL3KalmanTrack()
{
  // Destructor
}

Int_t AliL3KalmanTrack::Init(AliL3Track *track, AliL3SpacePointData *points, UInt_t pos, Int_t slice)
{

  Float_t xyz[3];

  xyz[0] = points[pos].fX;
  xyz[1] = points[pos].fY;

  // NB! I think boolean variable in the command under is true if single slice.
  // Better do a test if it is signle slice, and if so set it true.??.
  AliL3Transform::Global2Local(xyz,slice);
  
  fX = xyz[0];

  fP0 = xyz[1];
  fP1 = points[pos].fZ; 
  fP3 = track->GetTgl(); //cout << fP3; 
  if (TMath::Abs(fP3) > 1.2) return 0; //From AliTPCtrackerMI
  fP4 = 0.5*(-track->GetCharge()*1./track->GetPt())*0.0029980*AliL3Transform::GetBField(); // Calculation of the curvature 
  if (TMath::Abs(fP4) >= 0.0066) return 0; // From AliTPCtrackerMI


  Float_t firstXY[2];
  firstXY[0] = track->GetFirstPointX();
  firstXY[1] = track->GetFirstPointY();

  AliL3Transform::Global2Local(firstXY,slice);

  //cout << firstXY[0] << endl;

  //Float_t centerX = track->GetFirstPointX() - ((track->GetPt()/(0.0029980*AliL3Transform::GetBField())) * TMath::Cos(track->GetPsi() + track->GetCharge() * 0.5 * 3.14159265358979323846));
  Float_t centerX = firstXY[0] - ((track->GetPt()/(0.0029980*AliL3Transform::GetBField())) * TMath::Cos(track->GetPsi() + track->GetCharge() * 0.5 * 3.14159265358979323846));
  //Float_t centerY = track->GetFirstPointY() - ((track->GetPt()/(0.0029980*AliL3Transform::GetBField())) * TMath::Sin(track->GetPsi() + track->GetCharge() * 0.5 * 3.14159265358979323846));

  fP2 = fP4*centerX; // Curvature times center of curvature

  if (TMath::Abs(fP4*fX - fP2) >= 0.9) // What's this??
    {
      return 0;
    }
  //cout << track->GetCharge() << endl; 
  //cout << "AliL3KalmanTrack::Init, " << fP0 << " " << fP1 << " " << fP2 << " " << fP3 << " " << fP4 << endl;
  //cout << fP4 << endl;
  //Float_t num = 12;

  fC00 = points[pos].fSigmaY2; 
  fC10 = 0; fC11 = points[pos].fSigmaZ2;
  fC20 = 0; fC21 = 0; fC22 = 5e-05;
  fC30 = 0; fC31 = 0; fC32 = 0; fC33 = 5e-05;
  fC40 = 0; fC41 = 0; fC42 = 0; fC43 = 0; fC44 = 5e-09;
  /*fC00 = points[pos].fSigmaY2;
  fC10 = 0; fC11 = points[pos].fSigmaZ2;
  fC20 = 0; fC21 = 0; fC22 = (0.1/TMath::Sqrt(num))*(0.1/TMath::Sqrt(num));
  fC30 = 0; fC31 = -points[pos].fSigmaZ2; fC32 = 0; fC33 = 2*(points[pos].fSigmaZ2);
  fC40 = -points[pos].fSigmaY2; fC41 = 0; fC42 = 0; fC43 = 0; fC44 = 2*(points[pos].fSigmaY2);*/
  //cout << "Init: errors " << fC00 << " " << fC11 << " " << fC22 << " " << fC33 << " " << fC44 << endl;

  fChisq = 0;

  return 1;
}

Int_t AliL3KalmanTrack::Propagate(AliL3SpacePointData *points, UInt_t pos, Int_t slice)
{
  // Propagates track to the plane of the next found cluster
  Float_t Xold = fX; // X position for previous space point
  //Float_t Xnew = points[pos].fX; // X position of current space point
  //Float_t dx = Xnew - Xold; cout << Xnew << endl;
  Float_t Yold = fP0; // Y position of old point
  Float_t Zold = fP1; // Z position of old point
  Float_t par2old = fP2;  
  Float_t par3old = fP3;
  Float_t par4old = fP4;
  Float_t oldfC00 = fC00;
  Float_t oldfC10 = fC10;
  Float_t oldfC11 = fC11;
  Float_t oldfC20 = fC20;
  Float_t oldfC21 = fC21;
  Float_t oldfC22 = fC22;
  Float_t oldfC30 = fC30;
  Float_t oldfC31 = fC31;
  Float_t oldfC32 = fC32;
  Float_t oldfC33 = fC33;
  Float_t oldfC40 = fC40;
  Float_t oldfC41 = fC41;
  Float_t oldfC42 = fC42;
  Float_t oldfC43 = fC43;
  Float_t oldfC44 = fC44;

  Float_t xyz[3];

  xyz[0] = points[pos].fX;
  xyz[1] = points[pos].fY;

  AliL3Transform::Global2Local(xyz,slice);
  
  Float_t Xnew = xyz[0];
  Float_t dx = Xnew - Xold; //cout << Xnew << endl;
  
  if (TMath::Abs(fP4*Xnew - fP2) >= 0.9) // What's this??
    {
      //cout << "Propagation failed! Stiff track. " << pos << endl;
      return 0;
    }
  
  // R must be approx. of the radius of the circle based on
  // the old and new spacepoint. What then is Cod and Cnew??
  // C has something to do with curvature at least.

  Float_t Cold = fP4*Xold - fP2; 
  //if (TMath::Abs(Cold) >= 1.0) Cold = 0.9;
    /*{
    cout << "Cold " << endl << fP4*Xnew - fP2 << endl;
    return 0;
    }*/
  Float_t Rold = TMath::Sqrt(1 - Cold*Cold);
  Float_t Cnew = fP4*Xnew - fP2; 
  //if (TMath::Abs(Cnew) >= 1.0) Cnew = 0.9;
    /*{
    cout << "Cnew " << endl;
    return 0;
    }*/
  Float_t Rnew = TMath::Sqrt(1 - Cnew*Cnew);
  //if (Rold < 0.9) cout << "Cold = " << Cold << " , Rold = " << Rold << " , Cnew = " << Cnew << " , Rnew = " << Rnew << endl;

  // Prediction of the y- and z- coordinate in the next plane
  fP0 += dx*(Cold+Cnew)/(Rold+Rnew);
  fP1 += dx*(Cold+Cnew)/(Cold*Rnew + Cnew*Rold)*fP3; 
  //cout << "Old point " << Yold << " " << Zold << endl;
  //cout << "Propagate " << fP0 << " " << fP1 << endl;
  //cout << "Measured  " << points[pos].fY << " " << points[pos].fZ << endl;

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
  // Shold this be included??
  /*
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
  if (TMath::Abs(beta2) >= 1) cout << "here? " << beta2 << ", " << Pt << endl;
  
  // Energy loss assuming track is pion
  Float_t dE = 0.153e-3/beta2*(log(5940*beta2/(1-beta2))-beta2)*d*0.9e-3;
  if (Xold < Xnew) dE = -dE;
  CC = fP4;
  fP4 *= (1 - TMath::Sqrt(p2+0.14*0.14)/p2*dE);
  fP2 += Xnew*(fP4-CC);
  */
  // Update the track parameters with the measured values of the new point

  Int_t value = UpdateTrack(points, pos, slice);

  if (value == 0)
    {
      fP0 = Yold;
      fP1 = Zold;
      fP2 = par2old;
      fP3 = par3old;
      fP4 = par4old;
      fC00 = oldfC00;
      fC10 = oldfC10;
      fC11 = oldfC11;
      fC20 = oldfC20;
      fC21 = oldfC21;
      fC22 = oldfC22;
      fC30 = oldfC30;
      fC31 = oldfC31;
      fC32 = oldfC32;
      fC33 = oldfC33;
      fC40 = oldfC40;
      fC41 = oldfC41;
      fC42 = oldfC42;
      fC43 = oldfC43;
      fC44 = oldfC44;

      return value;
      }
  
  else 
    return value;
  //return UpdateTrack(points, pos);
}

Int_t AliL3KalmanTrack::UpdateTrack(AliL3SpacePointData *points, UInt_t pos, Int_t slice)
{
  // Update the track parameters with the measured values
  
  Float_t xyz[3];

  xyz[0] = points[pos].fX;
  xyz[1] = points[pos].fY;
  
  AliL3Transform::Global2Local(xyz,slice);
  
  //fX = points[pos].fX; 
  fX = xyz[0];

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

  // Kalman gain matrix
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
  //Float_t dy = points[pos].fY-fP0; //cout << "dy = " << dy;
  Float_t dy = xyz[1] - fP0; //cout << "dy = " << dy;
  Float_t dz = points[pos].fZ-fP1; //cout << ", dz = " << dz << endl; 
  //cout << "Measured " << xyz[2] << " " << points[pos].fZ << endl;

  // Prediction of fP2 and fP4
  Float_t cur = fP4 + k40*dy + k41*dz; 
  Float_t eta = fP2 + k20*dy + k21*dz;

  if (TMath::Abs(cur*fX - eta) >= 0.9)
    {
      //cout << "Update failed! Stiff track. " << pos << endl;
      return 0;
    }

  // Filtered state vector
  fP0 += k00*dy + k01*dz;
  fP1 += k10*dy + k11*dz;
  fP2 = eta;
  fP3 += k30*dy + k31*dz; //cout << "update " << fP3 << endl;
  fP4 = cur;
  //cout << "AliL3KalmanTrack::Update, " << fP0 << " " << fP1 << " " << fP2 << " " << fP3 << " " << fP4 << endl;
  //cout << "Measured, " << points[pos].fY << " " << points[pos].fZ << endl;

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
  //cout << "AliL3KalmanTrack::Update, error " << fC00 << " " << fC11 << " " << fC22 << " " << fC33 << " " << fC44 << endl;

  sigmaY2 = sigmaY2*det;
  sigmaZ2 = sigmaZ2*det;
  sigmaYZ = sigmaYZ*det;
  //cout << "AliL3KalmanTrack::Update, Chi2 = " << GetChisq() << endl;
  //cout << "AliKalmanTrack::Update, sigmaY2 = " << sigmaY2 << " sigmaZ2 = " << sigmaZ2 << " sigmaYZ = " << sigmaYZ << " dy = " << dy << " dz = " << dz << endl;
  // Calculate increase of chisquare
  fChisq = GetChisq() + GetChisqIncrement(xyz[1],sigmaY2,points[pos].fZ,sigmaZ2);
  //cout << "fChisq = " << fChisq << endl;
//(dy*sigmaY2*dy + 2*sigmaYZ*dy*dz + dz*sigmaZ2*dz) / (sigmaY2*sigmaZ2 - sigmaYZ*sigmaYZ);
  // Must at some point make an cut on chisq. Here?
  //if (fChisq > fMaxChi2) return 0;

  return 1;
} 

//Float_t AliL3KalmanTrack::GetChisqIncrement(AliL3SpacePointData *points, UInt_t pos)
Float_t AliL3KalmanTrack::GetChisqIncrement(Float_t y, Float_t error_y, Float_t z, Float_t error_z)
{
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  //Float_t r00=points[pos].fSigmaY2, r01=0., r11=points[pos].fSigmaZ2;
  Float_t r00=error_y, r01=0., r11=error_z;
  r00+=fC00; r01+=fC10; r11+=fC11;

  Double_t det=r00*r11 - r01*r01;
  /*if (TMath::Abs(det) < 1.e-10) {
    Int_t n=GetNumberOfClusters();
    if (n>4) cerr<<n<<" AliKalmanTrack warning: Singular matrix !\n";
    return 1e10;
    }*/
  Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;
  
  Double_t dy=y - fP0, dz=z - fP1;
  //cout << "AliTPCtrack::GetPredictedChi2, r00 = " << r00 << " r11 = " << r11 << " ro1 = " << r01 << " dy = " << dy << " dz = " << dz << endl;
  return (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det;
}

Float_t AliL3KalmanTrack::f2(Float_t x1,Float_t y1,
			     Float_t x2,Float_t y2,
			     Float_t x3,Float_t y3)
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature times center of curvature
  // Will be removed ??
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
  // Will be removed ??
  //-----------------------------------------------------------------
  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

Float_t AliL3KalmanTrack::f4(Float_t x1,Float_t y1,
			     Float_t x2,Float_t y2,
			     Float_t x3,Float_t y3)
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  // Will be removed??
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

Int_t AliL3KalmanTrack::PropagateOfflineTrack(Double_t x, Double_t y, Double_t z, Double_t ey, Double_t ez)
{
  // Propagates track to the plane of the next found cluster
  Float_t Xold = fX; // X position for previous space point
  Float_t Xnew = x; // X position of current space point
  Float_t dx = Xnew - Xold; //cout << dx << endl;
  Float_t Yold = fP0; // Y position of old point
  Float_t Zold = fP1; // Z position of old point
  Float_t par2old = fP2;
  Float_t par3old = fP3;
  Float_t par4old = fP4;
  Float_t oldfC00 = fC00;
  Float_t oldfC10 = fC10;
  Float_t oldfC11 = fC11;
  Float_t oldfC20 = fC20;
  Float_t oldfC21 = fC21;
  Float_t oldfC22 = fC22;
  Float_t oldfC30 = fC30;
  Float_t oldfC31 = fC31;
  Float_t oldfC32 = fC32;
  Float_t oldfC33 = fC33;
  Float_t oldfC40 = fC40;
  Float_t oldfC41 = fC41;
  Float_t oldfC42 = fC42;
  Float_t oldfC43 = fC43;
  Float_t oldfC44 = fC44;

  if (TMath::Abs(fP4*Xnew - fP2) >= 0.9) // What's this??
    {
      //cout << "Propagation failed! Stiff track. " << pos << endl;
      return 0;
    }

  if (TMath::Abs(dx) > 5) return 0;

  /*if (TMath::Abs(fP4*Xold - fP2) >= 0.9) // What's this??
    {
      return 0;
      }*/

  // R must be approx. of the radius of the circle based on
  // the old and new spacepoint. What then is Cod and Cnew??
  // C has something to do with curvature at least.
  
  Float_t Cold = fP4*Xold - fP2;

  //if (TMath::Abs(Cold) >= 1.0) return 0;
  /*{
    cout << "Cold " << endl << fP4*Xnew - fP2 << endl;
    return 0;
    }*/
  Float_t Rold = TMath::Sqrt(1 - Cold*Cold);
  Float_t Cnew = fP4*Xnew - fP2; 

  //if (TMath::Abs(Cnew) >= 1.0) Cnew = 0.9;
    /*{
    cout << "Cnew " << endl;
    return 0;
    }*/
  Float_t Rnew = TMath::Sqrt(1 - Cnew*Cnew);
  //if (Rold < 0.9) 
  //cout << "Cold = " << Cold << " , Rold = " << Rold << " , Cnew = " << Cnew << " , Rnew = " << Rnew << endl;  
  // Prediction of the y- and z- coordinate in the next plane
  fP0 += dx*(Cold+Cnew)/(Rold+Rnew);
  fP1 += dx*(Cold+Cnew)/(Cold*Rnew + Cnew*Rold)*fP3;
  //cout << "Old point " << Yold << " " << Zold << endl;
  //cout << "Propagate " << fP0 << " " << fP1 << endl;
  //cout << "Measured  " << points[pos].fY << " " << points[pos].fZ << endl;

  fX = Xnew;

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

  // Update the track parameters with the measured values of the new point
  Int_t value = UpdateOfflineTrack(x,y,z,ey,ez);
  
  if (value == 0)
    {
      fP0 = Yold;
      fP1 = Zold;
      fP2 = par2old;
      fP3 = par3old;
      fP4 = par4old;
      fC00 = oldfC00;
      fC10 = oldfC10;
      fC11 = oldfC11;
      fC20 = oldfC20;
      fC21 = oldfC21;
      fC22 = oldfC22;
      fC30 = oldfC30;
      fC31 = oldfC31;
      fC32 = oldfC32;
      fC33 = oldfC33;
      fC40 = oldfC40;
      fC41 = oldfC41;
      fC42 = oldfC42;
      fC43 = oldfC43;
      fC44 = oldfC44;
      
      return value;
      }
  
  else
    return value;
  //return 1;
  //return UpdateTrack(points, pos);

}

Int_t AliL3KalmanTrack::UpdateOfflineTrack(Double_t x, Double_t y, Double_t z, Double_t ey, Double_t ez)
{
  // Update the track parameters with the measured values
  //fX = x;

  // The errors from the measurement of the spacepoint
  Float_t sigmaY2 = ey;
  Float_t sigmaZ2 = ez;
  Float_t sigmaYZ = 0;
 
  sigmaY2 += fC00;
  sigmaZ2 += fC11;
  sigmaYZ += fC10;

  Float_t det = sigmaY2*sigmaZ2 - sigmaYZ*sigmaYZ;
  Float_t tmp = sigmaY2;
  sigmaY2 = sigmaZ2/det;
  sigmaZ2 = tmp/det;
  sigmaYZ = -sigmaYZ/det;
 
  // Kalman gain matrix
  Float_t k00 = fC00*sigmaY2 + fC10*sigmaYZ, k01 = fC00*sigmaYZ + fC10*sigmaZ2;
  Float_t k10 = fC10*sigmaY2 + fC11*sigmaYZ, k11 = fC10*sigmaYZ + fC11*sigmaZ2;
  Float_t k20 = fC20*sigmaY2 + fC21*sigmaYZ, k21 = fC20*sigmaYZ + fC21*sigmaZ2;
  Float_t k30 = fC30*sigmaY2 + fC31*sigmaYZ, k31 = fC30*sigmaYZ + fC31*sigmaZ2;
  Float_t k40 = fC40*sigmaY2 + fC41*sigmaYZ, k41 = fC40*sigmaYZ + fC41*sigmaZ2;
  //cout << "x = " << fX << endl;
  // Deviation between the predicted and measured values of y and z
  Float_t dy = y-fP0; //cout << "dy = " << dy;
  Float_t dz = z-fP1; //cout << ", dz = " << dz << endl;
  //cout << "Measured " << points[pos].fY << " " << points[pos].fZ << endl;
  
  // Prediction of fP2 and fP4
  Float_t cur = fP4 + k40*dy + k41*dz;
  Float_t eta = fP2 + k20*dy + k21*dz;
  
  if (TMath::Abs(cur*fX - eta) >= 0.9)
    {
      //cout << "Update failed! Stiff track. " << pos << endl;
      return 0;
    }

  // Filtered state vector
  fP0 += k00*dy + k01*dz;
  fP1 += k10*dy + k11*dz;
  fP2 = eta;
  fP3 += k30*dy + k31*dz; //cout << "update " << fP3 << endl;
  fP4 = cur;
  //cout << "AliL3KalmanTrack::Update, " << fP0 << " " << fP1 << " " << fP2 << " " << fP3 << " " << fP4 << endl;
  //cout << "Measured, " << points[pos].fY << " " << points[pos].fZ << endl;
  
  Float_t c01 = fC10;
  Float_t c02 = fC20;
  Float_t c03 = fC30;
  Float_t c04 = fC40;
  Float_t c12 = fC21;
  Float_t c13 = fC31;
  Float_t c14 = fC41;
  
  // Filtered covariance matrix
  fC00 -= k00*fC00 + k01*fC10; fC10 -= k00*c01 + k01*fC11;
  fC20 -= k00*c02 + k01*c12;   fC30 -= k00*c03 + k01*c13;
  fC40 -= k00*c04 + k01*c14;

  fC11 -= k10*c01 + k11*fC11;
  fC21 -= k10*c02 + k11*c12;   fC31 -= k10*c03 + k11*c13;
  fC41 -= k10*c04 + k11*c14;

  fC22 -= k20*c02 + k21*c12;   fC32 -= k20*c03 + k21*c13;
  fC42 -= k20*c04 + k21*c14;

  fC33 -= k30*c03 + k31*c13;
  fC43 -= k40*c03 + k41*c13;

  fC44 -= k40*c04 + k41*c14;

  //cout << "AliL3KalmanTrack::Update, error " << fC00 << " " << fC11 << " " << fC22 << " " << fC33 << " " << fC44 << endl;

  sigmaY2 = sigmaY2*det;
  sigmaZ2 = sigmaZ2*det;
  sigmaYZ = sigmaYZ*det;
  //cout << "AliL3KalmanTrack::Update, Chi2 = " << GetChisq() << endl;

  fChisq = GetChisq() + GetChisqIncrementOfflineTrack(y,z,ey,ez);
  //cout << "fChisq = " << fChisq << endl;

  // Must at some point make an cut on chisq. Here?
  //if (fChisq > fMaxChi2) return 0;
  
  return 1;
  
}

Float_t AliL3KalmanTrack::GetChisqIncrementOfflineTrack(Double_t y, Double_t z, Double_t ey, Double_t ez)
{
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  Float_t r00=ey, r01=0., r11=ez;
  r00+=fC00; r01+=fC10; r11+=fC11;

  Double_t det=r00*r11 - r01*r01;
  /*if (TMath::Abs(det) < 1.e-10) {
    Int_t n=GetNumberOfClusters();
    if (n>4) cerr<<n<<" AliKalmanTrack warning: Singular matrix !\n";
    return 1e10;
    }*/
  Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;
  
  Double_t dy=y - fP0, dz=z - fP1;
  //cout << "AliTPCtrack::GetPredictedChi2, r00 = " << r00 << " r11 = " << r11 << " ro1 = " << r01 << " dy = " << dy << " dz = " << dz << endl;
  return (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det;
}
