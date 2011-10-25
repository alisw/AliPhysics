/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliESDtrack.h" 
#include "AliTracker.h" 
#include "AliHMPIDtrack.h" 

ClassImp(AliHMPIDtrack)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDtrack::AliHMPIDtrack():AliKalmanTrack()
{
  //
  // def. ctor
  //
}                                
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDtrack::AliHMPIDtrack(const AliHMPIDtrack& t):AliKalmanTrack(t)
{
  //
  // cctor.
  //
 }                                
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDtrack::AliHMPIDtrack(const AliESDtrack& t):AliKalmanTrack()
{
  //
  // Constructor from AliESDtrack
  //
  SetLabel(t.GetLabel());
  SetChi2(0.);
  SetMass(t.GetMass());

  Set(t.GetX(),t.GetAlpha(),t.GetParameter(),t.GetCovariance());

  if ((t.GetStatus()&AliESDtrack::kTIME) == 0) return;
  StartTimeIntegral();
  Double_t times[10]; t.GetIntegratedTimes(times); SetIntegratedTimes(times);
  SetIntegratedLength(t.GetIntegratedLength());
}              
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDtrack& AliHMPIDtrack::operator=(const AliHMPIDtrack &/*source*/)
{
  // ass. op.
  
  return *this;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDtrack::PropagateTo(Double_t xk, Double_t xx0, Double_t xrho)
{
  //
  // Propagates this track to a reference plane defined by "xk" [cm] 
  // correcting for the mean crossed material.
  // Arguments:
  // "xx0"  - thickness/rad.length [units of the radiation length] 
  // "xrho" - thickness*density    [g/cm^2] 
  //  Returns: kTRUE if the track propagates to plane, else kFALSE
 
  if (xk == GetX()) {
    return kTRUE;
  }
  Double_t b[3]; GetBxByBz(b);
  if (!AliExternalTrackParam::PropagateToBxByBz(xk,b)) {
    return kFALSE;
  }
  if (!AliExternalTrackParam::CorrectForMeanMaterial(xx0,xrho,GetMass())) { 
    return kFALSE;
  }
  return kTRUE;            
}//PropagateTo()  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDtrack::Rotate(Double_t alpha, Bool_t absolute)
{
  //
  // Rotates track parameters in R*phi plane
  // if absolute rotation alpha is in global system
  // otherwise alpha rotation is relative to the current rotation angle
  //  

  if (absolute) {
    alpha -= GetAlpha();
  }
 

  return AliExternalTrackParam::Rotate(GetAlpha()+alpha);

}   
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
Int_t AliHMPIDtrack::GetProlongation(Double_t xk, Double_t &y, Double_t &z)
{
  //
  // Find a prolongation at given x
  // Return 0 if it does not exist
  //  

  Double_t bz = GetBz();

  if (!AliExternalTrackParam::GetYAt(xk,bz,y)) {
    return 0;
  }
  if (!AliExternalTrackParam::GetZAt(xk,bz,z)) {
    return 0;
  }

  return 1;  

}
 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t   AliHMPIDtrack::PropagateToR(Double_t r,Double_t step)
{
  //
  // Propagate track to the radial position
  // Rotation always connected to the last track position
  //

  Double_t xyz0[3];
  Double_t xyz1[3];
  Double_t y;
  Double_t z; 

  Double_t radius = TMath::Sqrt(GetX()*GetX() + GetY()*GetY());
  // Direction +-
  Double_t dir    = (radius > r) ? -1.0 : 1.0;   

  for (Double_t x = radius+dir*step; dir*x < dir*r; x += dir*step) {

    GetXYZ(xyz0);	
    Double_t alpha = TMath::ATan2(xyz0[1],xyz0[0]);
    if(!Rotate(alpha,kTRUE)) return kFALSE;
    GetXYZ(xyz0);	
    if (!GetProlongation(x,y,z)) return kFALSE;
    xyz1[0] = x * TMath::Cos(alpha) + y * TMath::Sin(alpha); 
    xyz1[1] = x * TMath::Sin(alpha) - y * TMath::Cos(alpha);
    xyz1[2] = z;
    Double_t param[7];
    AliTracker::MeanMaterialBudget(xyz0,xyz1,param);
    if (param[1] <= 0) {
      param[1] = 100000000;
    }
    PropagateTo(x,param[1],param[0]*param[4]);

  } 

  GetXYZ(xyz0);	
  Double_t alpha = TMath::ATan2(xyz0[1],xyz0[0]);
  if(!Rotate(alpha,kTRUE)) return kFALSE;
  GetXYZ(xyz0);	
  if (!GetProlongation(r,y,z)) return kFALSE;
  xyz1[0] = r * TMath::Cos(alpha) + y * TMath::Sin(alpha); 
  xyz1[1] = r * TMath::Sin(alpha) - y * TMath::Cos(alpha);
  xyz1[2] = z;
  Double_t param[7];
  AliTracker::MeanMaterialBudget(xyz0,xyz1,param);

  if (param[1] <= 0) {
    param[1] = 100000000;
  }
 
  return PropagateTo(r,param[1],param[0]*param[4]);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDtrack::GetPredictedChi2(const AliCluster3D *c) const {
  //
  // Arguments: AliCluster3D
  // Returns:   Chi2 of track for the cluster
  Double_t      p[3]={c->GetX(),       c->GetY(),       c->GetZ()};
  Double_t  covyz[3]={c->GetSigmaY2(), c->GetSigmaYZ(), c->GetSigmaZ2()};
  Double_t covxyz[3]={c->GetSigmaX2(), c->GetSigmaXY(), c->GetSigmaXZ()};
  return AliExternalTrackParam::GetPredictedChi2(p, covyz, covxyz);
}//GetPredictedChi2()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDtrack::PropagateTo(const AliCluster3D *c) {
  //
  // Arguments: AliCluster3D
  // Returns: kTRUE if the track propagates to the plane of the cluster
  Double_t      oldX=GetX(),   oldY=GetY(), oldZ=GetZ();
  Double_t      p[3]={c->GetX(), c->GetY(), c->GetZ()};
  Double_t  covyz[3]={c->GetSigmaY2(), c->GetSigmaYZ(), c->GetSigmaZ2()};
  Double_t covxyz[3]={c->GetSigmaX2(), c->GetSigmaXY(), c->GetSigmaXZ()};
  Double_t bz=-GetBz();
    
  if(!AliExternalTrackParam::PropagateTo(p, covyz, covxyz, bz)) return kFALSE;
  if(IsStartedTimeIntegral()) 
    {
      Double_t d = TMath::Sqrt((GetX()-oldX)*(GetX()-oldX) + (GetY()-oldY)*(GetY()-oldY) + (GetZ()-oldZ)*(GetZ()-oldZ));
      if (GetX()<oldX) d=-d;
      AddTimeStep(d);
    }
  return kTRUE;
}//PropagateTo()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDtrack::Intersect(Double_t pnt[3], Double_t norm[3]) const {
  //+++++++++++++++++++++++++++++++++++++++++    
  // Origin: K. Shileev (Kirill.Shileev@cern.ch)
  // Finds point of intersection (if exists) of the helix with the plane. 
  // Stores result in fX and fP.   
  // Arguments: planePoint,planeNorm - the plane defined by any plane's point 
  // and vector, normal to the plane
  // Returns: kTrue if helix intersects the plane, kFALSE otherwise.
  //+++++++++++++++++++++++++++++++++++++++++    
  Double_t x0[3]; GetXYZ(x0); //get track position in MARS
  
  //estimates initial helix length up to plane
  Double_t s=(pnt[0]-x0[0])*norm[0] + (pnt[1]-x0[1])*norm[1] + (pnt[2]-x0[2])*norm[2];
  Double_t dist=99999,distPrev=dist;
  Double_t p[3],x[3]; 
  while(TMath::Abs(dist)>0.00001){
    //calculates helix at the distance s from x0 ALONG the helix
    Propagate(s,x,p);
    //distance between current helix position and plane
    dist=(x[0]-pnt[0])*norm[0]+(x[1]-pnt[1])*norm[1]+(x[2]-pnt[2])*norm[2];  
    if(TMath::Abs(dist) >= TMath::Abs(distPrev)) {return kFALSE;}
    distPrev=dist;
    s-=dist;
  }
  //on exit pnt is intersection point,norm is track vector at that point, 
  //all in MARS
  for (Int_t i=0; i<3; i++) {pnt[i]=x[i]; norm[i]=p[i];}
  return kTRUE;
}//Intersect()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDtrack::Propagate(Double_t len, Double_t x[3],Double_t p[3]) const {
  //+++++++++++++++++++++++++++++++++++++++++    
  // Origin: K. Shileev (Kirill.Shileev@cern.ch)
  // Extrapolate track along simple helix in magnetic field
  // Arguments: len -distance alogn helix, [cm]
  //            bz  - mag field, [kGaus]   
  // Returns: x and p contain extrapolated positon and momentum  
  // The momentum returned for straight-line tracks is meaningless !
  //+++++++++++++++++++++++++++++++++++++++++    
  GetXYZ(x);    
  Double_t bField[3];
  TGeoGlobalMagField::Instance()->Field(x,bField);
  Double_t bz = -bField[2];
  if (OneOverPt() < kAlmost0 || TMath::Abs(bz) < kAlmost0Field ){ //straight-line tracks
     Double_t unit[3]; GetDirection(unit);
     x[0]+=unit[0]*len;   
     x[1]+=unit[1]*len;   
     x[2]+=unit[2]*len;

     p[0]=unit[0]/kAlmost0;   
     p[1]=unit[1]/kAlmost0;   
     p[2]=unit[2]/kAlmost0;   
  } else {
     GetPxPyPz(p);
     Double_t pp=GetP();
     Double_t a = -kB2C*bz*GetSign(); ////////// what is kB2C
     Double_t rho = a/pp;
     x[0] += p[0]*TMath::Sin(rho*len)/a - p[1]*(1-TMath::Cos(rho*len))/a;
     x[1] += p[1]*TMath::Sin(rho*len)/a + p[0]*(1-TMath::Cos(rho*len))/a;
     x[2] += p[2]*len/pp;
     Double_t p0=p[0];
     p[0] = p0  *TMath::Cos(rho*len) - p[1]*TMath::Sin(rho*len);
     p[1] = p[1]*TMath::Cos(rho*len) + p0  *TMath::Sin(rho*len);
  }
}//Propagate()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDtrack::Update(const AliHMPIDCluster *pClu, Double_t /*chisq*/, Int_t /*index*/)
{
  //
  // Arguments: AliCluster3D, chi sq, and clu index
  // Returns: kTRUE if the track parameters are successfully updated
  Double_t p[2]={pClu->GetY(), pClu->GetZ()};
  Double_t cov[3]={pClu->GetSigmaY2(), 0., pClu->GetSigmaZ2()};
  if (!AliExternalTrackParam::Update(p,cov)) return kFALSE;

  /*
  AliTracker::FillResiduals(this,p,cov,pClu->GetVolumeId());

  Int_t n=GetNumberOfClusters();
  fIndex[n]=index;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2()+chisq);
*/
  return kTRUE;

}//Update()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
