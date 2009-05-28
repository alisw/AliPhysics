// @(#) $Id$
// Original: AliHLTTrack.cxx,v 1.32 2005/06/14 10:55:21 cvetan 

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Anders Vestbo, Uli Frankenfeld, maintained by         *
 *                  Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCTrack.cxx
    @author Anders Vestbo, Uli Frankenfeld, maintained by Matthias Richter
    @date   
    @brief  HLT TPC track implementation (conformal mapping) */


#include "AliHLTTPCLogging.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCSpacePointData.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCTrack)


AliHLTTPCTrack::AliHLTTPCTrack()
  :
  fNHits(0),
  fMCid(-1),
  fKappa(0),
  fRadius(0),
  fCenterX(0),
  fCenterY(0),
  fFromMainVertex(0),
  fSector(0),
  fQ(0),

  fTanl(0),
  fPsi(0),
  fPt(0),
  fLength(0),

  fPterr(0),
  fPsierr(0),
  fZ0err(0),
  fY0err(0),
  fTanlerr(0),

  fPhi0(0),
  fR0(0),
  fZ0(0),

  //  fPoint({0,0,0}),
  fPointPsi(0),

  fIsPoint(0),
  fIsLocal(true),
  //  fRowRange({0,0}),

  fPID(0)
{
  //Constructor
  fRowRange[0]=0;
  fRowRange[1]=0;
  fPoint[0]=0;
  fPoint[1]=0;
  fPoint[2]=0;

  SetFirstPoint(0,0,0);
  SetLastPoint(0,0,0);
  memset(fHitNumbers,0,159*sizeof(UInt_t));
}

void AliHLTTPCTrack::Copy(AliHLTTPCTrack *tpt)
{
  //setter
  SetRowRange(tpt->GetFirstRow(),tpt->GetLastRow());
  SetPhi0(tpt->GetPhi0());
  SetKappa(tpt->GetKappa());
  SetNHits(tpt->GetNHits());
  SetFirstPoint(tpt->GetFirstPointX(),tpt->GetFirstPointY(),tpt->GetFirstPointZ());
  SetLastPoint(tpt->GetLastPointX(),tpt->GetLastPointY(),tpt->GetLastPointZ());
  SetPt(tpt->GetPt());
  SetPsi(tpt->GetPsi());
  SetTgl(tpt->GetTgl());
  SetPterr(tpt->GetPterr());
  SetPsierr(tpt->GetPsierr());
  SetTglerr(tpt->GetTglerr());
  SetZ0err(tpt->GetZ0err());
  SetY0err(tpt->GetY0err());
  SetCharge(tpt->GetCharge());
  SetHits(tpt->GetNHits(),(UInt_t *)tpt->GetHitNumbers());
#ifdef do_mc
  SetMCid(tpt->GetMCid());
#endif
  SetPID(tpt->GetPID());
  SetSector(tpt->GetSector());
}

Int_t AliHLTTPCTrack::Compare(const AliHLTTPCTrack *track) const
{
  // compare tracks
  if(track->GetNHits() < GetNHits()) return 1;
  if(track->GetNHits() > GetNHits()) return -1;
  return 0;
}

AliHLTTPCTrack::~AliHLTTPCTrack()
{
  //Nothing to do
}

Double_t AliHLTTPCTrack::GetP() const
{
  // Returns total momentum.  
  return fabs(GetPt())*sqrt(1. + GetTgl()*GetTgl());
}

Double_t AliHLTTPCTrack::GetPseudoRapidity() const
{ //get pseudo rap
  return 0.5 * log((GetP() + GetPz()) / (GetP() - GetPz()));
}

/*
Double_t AliHLTTPCTrack::GetEta() const
{
  return GetPseudoRapidity();
}
*/

Double_t AliHLTTPCTrack::GetRapidity() const
{ 
  //get rap
  const Double_t kmpi = 0.13957;
  return 0.5 * log((kmpi + GetPz()) / (kmpi - GetPz()));
}

void AliHLTTPCTrack::Rotate(Int_t slice,Bool_t tolocal)
{
  //Rotate track to global parameters
  //If flag tolocal is set, the track is rotated
  //to local coordinates.

  Float_t psi[1] = {GetPsi()};
  if(!tolocal)
    AliHLTTPCTransform::Local2GlobalAngle(psi,slice);
  else
    AliHLTTPCTransform::Global2LocalAngle(psi,slice);
  SetPsi(psi[0]);
  Float_t first[3];
  first[0] = GetFirstPointX();
  first[1] = GetFirstPointY();
  first[2] = GetFirstPointZ();
  if(!tolocal)
    AliHLTTPCTransform::Local2Global(first,slice);
  else
    AliHLTTPCTransform::Global2LocHLT(first,slice);
  //AliHLTTPCTransform::Global2Local(first,slice,kTRUE);
  
  SetFirstPoint(first[0],first[1],first[2]);
  Float_t last[3];
  last[0] = GetLastPointX();
  last[1] = GetLastPointY();
  last[2] = GetLastPointZ();
  if(!tolocal)
    AliHLTTPCTransform::Local2Global(last,slice);
  else
    AliHLTTPCTransform::Global2LocHLT(last,slice);    
  //AliHLTTPCTransform::Global2Local(last,slice,kTRUE);
  SetLastPoint(last[0],last[1],last[2]);
  
  Float_t center[3] = {GetCenterX(),GetCenterY(),0};
  if(!tolocal)
    AliHLTTPCTransform::Local2Global(center,slice);
  else
    AliHLTTPCTransform::Global2LocHLT(center,slice);
  //AliHLTTPCTransform::Global2Local(center,slice,kTRUE);
  SetCenterX(center[0]);
  SetCenterY(center[1]);
  
  SetPhi0(atan2(fFirstPoint[1],fFirstPoint[0]));
  SetR0(sqrt(fFirstPoint[0]*fFirstPoint[0]+fFirstPoint[1]*fFirstPoint[1]));
  
  if(!tolocal)
    fIsLocal=kFALSE;
  else
    fIsLocal=kTRUE;
}

void AliHLTTPCTrack::CalculateHelix()
{
  // fit assigned clusters to helix
  // for straight line fit
  if (AliHLTTPCTransform::GetBFieldValue() == 0.0 ){
    fRadius = 999999;  //just zero
    
    SetPhi0(atan2(fFirstPoint[1],fFirstPoint[0]));
    SetR0(sqrt(fFirstPoint[0]*fFirstPoint[0]+fFirstPoint[1]*fFirstPoint[1]));
  }
  // for helix fit
  else { 
    //Calculate Radius, CenterX and CenterY from Psi, X0, Y0
    fRadius = fPt / (AliHLTTPCTransform::GetBFieldValue());
    if(fRadius) fKappa = -fQ*1./fRadius;
    else fRadius = 999999;  //just zero
    Double_t trackPhi0 = fPsi + fQ * AliHLTTPCTransform::PiHalf();
    
    fCenterX = fFirstPoint[0] - fRadius *  cos(trackPhi0);
    fCenterY = fFirstPoint[1] - fRadius *  sin(trackPhi0);
    
    SetPhi0(atan2(fFirstPoint[1],fFirstPoint[0]));
    SetR0(sqrt(fFirstPoint[0]*fFirstPoint[0]+fFirstPoint[1]*fFirstPoint[1]));
  }
}

Double_t AliHLTTPCTrack::GetCrossingAngle(Int_t padrow,Int_t slice) 
{
  //Calculate the crossing angle between track and given padrow.
  //Take the dot product of the tangent vector of the track, and
  //vector perpendicular to the padrow.
  //In order to do this, we need the tangent vector to the track at the
  //point. This is done by rotating the radius vector by 90 degrees;
  //rotation matrix: (  0  1 )
  //                 ( -1  0 )

  Float_t angle=0;//Angle perpendicular to the padrow in local coordinates
  if(slice>=0)//Global coordinates
    {
      AliHLTTPCTransform::Local2GlobalAngle(&angle,slice);
      if(!CalculateReferencePoint(angle,AliHLTTPCTransform::Row2X(padrow)))
	cerr<<"AliHLTTPCTrack::GetCrossingAngle : Track does not cross line in slice "<<slice<<" row "<<padrow<<endl;
    }
  else //should be in local coordinates
    {
      Float_t xyz[3];
      GetCrossingPoint(padrow,xyz);
      fPoint[0] = xyz[0];
      fPoint[1] = xyz[1];
      fPoint[2] = xyz[2];
    }
    
  Double_t tangent[2];
  
  tangent[0] = (fPoint[1] - GetCenterY())/GetRadius();
  tangent[1] = -1.*(fPoint[0] - GetCenterX())/GetRadius();

  Double_t perppadrow[2] = {cos(angle),sin(angle)}; 
  Double_t cosbeta = fabs(tangent[0]*perppadrow[0] + tangent[1]*perppadrow[1]);
  if(cosbeta > 1) cosbeta=1;
  return acos(cosbeta);
}

Bool_t AliHLTTPCTrack::GetCrossingPoint(Int_t padrow,Float_t *xyz)
{
  //Assumes the track is given in local coordinates
  if(!IsLocal())
    {
      cerr<<"GetCrossingPoint: Track is given on global coordinates"<<endl;
      return false;
    }
  
  Double_t xHit = AliHLTTPCTransform::Row2X(padrow);

//if (xHit < xyz[0]){
//    LOG(AliHLTTPCLog::kError,"AliHLTTPCTRACK::GetCrossingPoint","")<< "Track doesn't cross padrow " 
//				<< padrow <<"(x=" << xHit << "). Smallest x=" << xyz[0] << ENDLOG;
//      return false;
//}

  // for straight line fit
  if (AliHLTTPCTransform::GetBFieldValue() == 0.0 ){
    
    Double_t yHit = GetFirstPointY() + (Double_t) tan( GetPsi() ) * (xHit - GetFirstPointX());   
    
    Double_t s = (xHit - GetFirstPointX())*(xHit - GetFirstPointX()) + (yHit - GetFirstPointY())*(yHit - GetFirstPointY()); 
    
    Double_t zHit = GetFirstPointZ() + s * GetTgl();
    
    xyz[0] = xHit;
    xyz[1] = yHit;
    xyz[2] = zHit;
  }
  // for helix fit
    else { 
      xyz[0] = xHit;
      Double_t aa = (xHit - GetCenterX())*(xHit - GetCenterX());
      Double_t r2 = GetRadius()*GetRadius();
      if(aa > r2)
	return false;
      
      Double_t aa2 = sqrt(r2 - aa);
      Double_t y1 = GetCenterY() + aa2;
      Double_t y2 = GetCenterY() - aa2;
      xyz[1] = y1;
      if(fabs(y2) < fabs(y1)) xyz[1] = y2;
      
      Double_t yHit = xyz[1];
      Double_t angle1 = atan2((yHit - GetCenterY()),(xHit - GetCenterX()));
      if(angle1 < 0) angle1 += 2.*AliHLTTPCTransform::Pi();
      Double_t angle2 = atan2((GetFirstPointY() - GetCenterY()),(GetFirstPointX() - GetCenterX()));
      if(angle2 < 0) angle2 += AliHLTTPCTransform::TwoPi();
      
      Double_t diffangle = angle1 - angle2;
      diffangle = fmod(diffangle,AliHLTTPCTransform::TwoPi());
      if((GetCharge()*diffangle) > 0) diffangle = diffangle - GetCharge()*AliHLTTPCTransform::TwoPi();
      
      Double_t stot = fabs(diffangle)*GetRadius();
      
      Double_t zHit = GetFirstPointZ() + stot*GetTgl();
      
      xyz[2] = zHit;
    }
  
  return true;
}

Bool_t AliHLTTPCTrack::CalculateReferencePoint(Double_t angle,Double_t radius)
{
  // Global coordinate: crossing point with y = ax+ b; 
  // a=tan(angle-AliHLTTPCTransform::PiHalf());
  //
  const Double_t krr=radius; //position of reference plane
  const Double_t kxr = cos(angle) * krr;
  const Double_t kyr = sin(angle) * krr;
  
  Double_t a = tan(angle-AliHLTTPCTransform::PiHalf());
  Double_t b = kyr - a * kxr;

  Double_t pp=(fCenterX+a*fCenterY-a*b)/(1+pow(a,2));
  Double_t qq=(pow(fCenterX,2)+pow(fCenterY,2)-2*fCenterY*b+pow(b,2)-pow(fRadius,2))/(1+pow(a,2));

  Double_t racine = pp*pp-qq;
  if(racine<0) return IsPoint(kFALSE);      //no Point

  Double_t rootRacine = sqrt(racine);
  Double_t x0 = pp+rootRacine;
  Double_t x1 = pp-rootRacine;
  Double_t y0 = a*x0 + b;
  Double_t y1 = a*x1 + b;

  Double_t diff0 = sqrt(pow(x0-kxr,2)+pow(y0-kyr,2));
  Double_t diff1 = sqrt(pow(x1-kxr,2)+pow(y1-kyr,2));
 
  if(diff0<diff1){
    fPoint[0]=x0;
    fPoint[1]=y0;
  }
  else{
    fPoint[0]=x1;
    fPoint[1]=y1;
  }

  Double_t pointPhi0  = atan2(fPoint[1]-fCenterY,fPoint[0]-fCenterX);
  Double_t trackPhi0  = atan2(fFirstPoint[1]-fCenterY,fFirstPoint[0]-fCenterX);
  if(fabs(trackPhi0-pointPhi0)>AliHLTTPCTransform::Pi()){
    if(trackPhi0<pointPhi0) trackPhi0 += AliHLTTPCTransform::TwoPi();
    else                    pointPhi0 += AliHLTTPCTransform::TwoPi();
  }
  Double_t stot = -fQ * (pointPhi0-trackPhi0) * fRadius ;
  fPoint[2]   = fFirstPoint[2] + stot * fTanl;

  fPointPsi = pointPhi0 - fQ * AliHLTTPCTransform::PiHalf();
  if(fPointPsi<0.)  fPointPsi+= AliHLTTPCTransform::TwoPi();
  fPointPsi = fmod(fPointPsi, AliHLTTPCTransform::TwoPi());

  return IsPoint(kTRUE);
}

Bool_t AliHLTTPCTrack::CalculateEdgePoint(Double_t angle)
{
  // Global coordinate: crossing point with y = ax; a=tan(angle);
  //
  Double_t rmin=AliHLTTPCTransform::Row2X(AliHLTTPCTransform::GetFirstRow(-1));  //min Radius of TPC
  Double_t rmax=AliHLTTPCTransform::Row2X(AliHLTTPCTransform::GetLastRow(-1)); //max Radius of TPC

  Double_t a = tan(angle);
  Double_t pp=(fCenterX+a*fCenterY)/(1+pow(a,2));
  Double_t qq=(pow(fCenterX,2)+pow(fCenterY,2)-pow(fRadius,2))/(1+pow(a,2));
  Double_t racine = pp*pp-qq;
  if(racine<0) return IsPoint(kFALSE);     //no Point
  Double_t rootRacine = sqrt(racine);
  Double_t x0 = pp+rootRacine;
  Double_t x1 = pp-rootRacine;
  Double_t y0 = a*x0;
  Double_t y1 = a*x1;

  Double_t r0 = sqrt(pow(x0,2)+pow(y0,2));
  Double_t r1 = sqrt(pow(x1,2)+pow(y1,2)); 
  //find the right crossing point:
  //inside the TPC modules
  Bool_t ok0 = kFALSE;
  Bool_t ok1 = kFALSE;

  if(r0>rmin&&r0<rmax){
    Double_t da=atan2(y0,x0);
    if(da<0) da+=AliHLTTPCTransform::TwoPi();
    if(fabs(da-angle)<0.5)
      ok0 = kTRUE;
  }
  if(r1>rmin&&r1<rmax){
    Double_t da=atan2(y1,x1);
    if(da<0) da+=AliHLTTPCTransform::TwoPi();
    if(fabs(da-angle)<0.5)
      ok1 = kTRUE;
  }
  if(!(ok0||ok1)) return IsPoint(kFALSE);   //no Point
  
  if(ok0&&ok1){
    Double_t diff0 = sqrt(pow(fFirstPoint[0]-x0,2)+pow(fFirstPoint[1]-y0,2));
    Double_t diff1 = sqrt(pow(fFirstPoint[0]-x1,2)+pow(fFirstPoint[1]-y1,2));
    if(diff0<diff1) ok1 = kFALSE; //use ok0
    else ok0 = kFALSE;            //use ok1
  }
  if(ok0){fPoint[0]=x0; fPoint[1]=y0;}
  else   {fPoint[0]=x1; fPoint[1]=y1;}

  Double_t pointPhi0  = atan2(fPoint[1]-fCenterY,fPoint[0]-fCenterX);
  Double_t trackPhi0  = atan2(fFirstPoint[1]-fCenterY,fFirstPoint[0]-fCenterX);
  if(fabs(trackPhi0-pointPhi0)>AliHLTTPCTransform::Pi()){
    if(trackPhi0<pointPhi0) trackPhi0 += AliHLTTPCTransform::TwoPi();
    else                    pointPhi0 += AliHLTTPCTransform::TwoPi();
  }
  Double_t stot = -fQ * (pointPhi0-trackPhi0) * fRadius ;
  fPoint[2]   = fFirstPoint[2] + stot * fTanl;

  fPointPsi = pointPhi0 - fQ * AliHLTTPCTransform::PiHalf();
  if(fPointPsi<0.)  fPointPsi+= AliHLTTPCTransform::TwoPi();
  fPointPsi = fmod(fPointPsi, AliHLTTPCTransform::TwoPi());

  return IsPoint(kTRUE);
}

Bool_t AliHLTTPCTrack::CalculatePoint(Double_t xplane)
{
  // Local coordinate: crossing point with x plane
  //
  Double_t racine = pow(fRadius,2)-pow(xplane-fCenterX,2);
  if(racine<0) return IsPoint(kFALSE);
  Double_t rootRacine = sqrt(racine);

  Double_t y0 = fCenterY + rootRacine;
  Double_t y1 = fCenterY - rootRacine;
  //Double_t diff0 = sqrt(pow(fFirstPoint[0]-xplane)+pow(fFirstPoint[1]-y0));
  //Double_t diff1 = sqrt(pow(fFirstPoint[0]-xplane)+pow(fFirstPoint[1]-y1));
  Double_t diff0 = fabs(y0-fFirstPoint[1]);
  Double_t diff1 = fabs(y1-fFirstPoint[1]);

  fPoint[0]=xplane;
  if(diff0<diff1) fPoint[1]=y0;
  else            fPoint[1]=y1;

  Double_t pointPhi0  = atan2(fPoint[1]-fCenterY,fPoint[0]-fCenterX);
  Double_t trackPhi0  = atan2(fFirstPoint[1]-fCenterY,fFirstPoint[0]-fCenterX);
  if(fabs(trackPhi0-pointPhi0)>AliHLTTPCTransform::Pi()){
    if(trackPhi0<pointPhi0) trackPhi0 += AliHLTTPCTransform::TwoPi();
    else                    pointPhi0 += AliHLTTPCTransform::TwoPi();
  }
  Double_t stot = -fQ * (pointPhi0-trackPhi0) * fRadius ;  
  fPoint[2]   = fFirstPoint[2] + stot * fTanl;

  fPointPsi = pointPhi0 - fQ * AliHLTTPCTransform::PiHalf();
  if(fPointPsi<0.)  fPointPsi+= AliHLTTPCTransform::TwoPi();
  fPointPsi = fmod(fPointPsi, AliHLTTPCTransform::TwoPi());

  return IsPoint(kTRUE);
}

void AliHLTTPCTrack::UpdateToFirstPoint()
{
  //Update track parameters to the innermost point on the track.
  //This means that the parameters of the track will be given in the point
  //of closest approach to the first innermost point, i.e. the point 
  //lying on the track fit (and not the coordinates of the innermost point itself).
  //This function assumes that fFirstPoint is already set to the coordinates of the innermost
  //assigned cluster.
  //
  //During the helix-fit, the first point on the track is set to the coordinates
  //of the innermost assigned cluster. This may be ok, if you just want a fast
  //estimate of the "global" track parameters; such as the momentum etc.
  //However, if you later on want to do more precise local calculations, such
  //as impact parameter, residuals etc, you need to give the track parameters
  //according to the actual fit.
  // for straight line fit
  if (AliHLTTPCTransform::GetBFieldValue() == 0.0 ){
    Double_t xc = GetCenterX() - GetFirstPointX();
    Double_t yc = GetCenterY() - GetFirstPointY();
    
    Double_t xn = (Double_t) sin( GetPsi() );
    Double_t yn = -1. * (Double_t) cos( GetPsi() );
    
    Double_t d = xc*xn + yc*yn;
    
    Double_t distx = d * xn;
    Double_t disty = d * yn;
    
    Double_t point[2];
    
    point[0] = distx + GetFirstPointX();
    point[1] = disty + GetFirstPointY();
    
    //Update the track parameters
    SetR0(sqrt(point[0]*point[0]+point[1]*point[1]));
    SetPhi0(atan2(point[1],point[0]));
    SetFirstPoint(point[0],point[1],GetZ0());
  }
  // for helix fit
  else { 
    Double_t xc = GetCenterX() - GetFirstPointX();
    Double_t yc = GetCenterY() - GetFirstPointY();
    
    Double_t distx1 = xc*(1 + GetRadius()/sqrt(xc*xc + yc*yc));
    Double_t disty1 = yc*(1 + GetRadius()/sqrt(xc*xc + yc*yc));
    Double_t distance1 = sqrt(distx1*distx1 + disty1*disty1);
    
    Double_t distx2 = xc*(1 - GetRadius()/sqrt(xc*xc + yc*yc));
    Double_t disty2 = yc*(1 - GetRadius()/sqrt(xc*xc + yc*yc));
    Double_t distance2 = sqrt(distx2*distx2 + disty2*disty2);
    
    //Choose the closest:
    Double_t point[2];
    if(distance1 < distance2)
      {
	point[0] = distx1 + GetFirstPointX();
	point[1] = disty1 + GetFirstPointY();
      }
    else
      {
	point[0] = distx2 + GetFirstPointX();
	point[1] = disty2 + GetFirstPointY();
      }
    
    Double_t pointpsi = atan2(point[1]-GetCenterY(),point[0]-GetCenterX());
    pointpsi -= GetCharge()*AliHLTTPCTransform::PiHalf();
    if(pointpsi < 0) pointpsi += AliHLTTPCTransform::TwoPi();
    
    //Update the track parameters
    SetR0(sqrt(point[0]*point[0]+point[1]*point[1]));
    SetPhi0(atan2(point[1],point[0]));
    SetFirstPoint(point[0],point[1],GetZ0());
    SetPsi(pointpsi);
  }
}

void AliHLTTPCTrack::GetClosestPoint(AliHLTTPCVertex *vertex,Double_t &closestX,Double_t &closestY,Double_t &closestZ)
{
  //Calculate the point of closest approach to the vertex
  //This function calculates the minimum distance from the helix to the vertex, and choose 
  //the corresponding point lying on the helix as the point of closest approach.
  
  Double_t xc = GetCenterX() - vertex->GetX();
  Double_t yc = GetCenterY() - vertex->GetY();
  
  Double_t distx1 = xc*(1 + GetRadius()/sqrt(xc*xc + yc*yc));
  Double_t disty1 = yc*(1 + GetRadius()/sqrt(xc*xc + yc*yc));
  Double_t distance1 = sqrt(distx1*distx1 + disty1*disty1);
  
  Double_t distx2 = xc*(1 - GetRadius()/sqrt(xc*xc + yc*yc));
  Double_t disty2 = yc*(1 - GetRadius()/sqrt(xc*xc + yc*yc));
  Double_t distance2 = sqrt(distx2*distx2 + disty2*disty2);
  
  //Choose the closest:
  if(distance1 < distance2)
    {
      closestX = distx1 + vertex->GetX();
      closestY = disty1 + vertex->GetY();
    }
  else
    {
      closestX = distx2 + vertex->GetX();
      closestY = disty2 + vertex->GetY();
    }
  
  //Get the z coordinate:
  Double_t angle1 = atan2((closestY-GetCenterY()),(closestX-GetCenterX()));
  if(angle1 < 0) angle1 = angle1 + AliHLTTPCTransform::TwoPi();
 
  Double_t angle2 = atan2((GetFirstPointY()-GetCenterY()),(GetFirstPointX()-GetCenterX()));
  if(angle2 < 0) angle2 = angle2 + AliHLTTPCTransform::TwoPi();
  
  Double_t diffAngle = angle1 - angle2;
  diffAngle = fmod(diffAngle,AliHLTTPCTransform::TwoPi());
  
  if((GetCharge()*diffAngle) < 0) diffAngle = diffAngle + GetCharge()*AliHLTTPCTransform::TwoPi();
  Double_t stot = fabs(diffAngle)*GetRadius();
  closestZ = GetFirstPointZ() - stot*GetTgl();
}

void AliHLTTPCTrack::Print(Option_t* /*option*/) const
{ 
//print out parameters of track

 LOG(AliHLTTPCLog::kInformational,"AliHLTTPCTrack::Print","Print values")
    <<"NH="<<fNHits<<" "<<fMCid<<" K="<<fKappa<<" R="<<fRadius<<" Cx="<<fCenterX<<" Cy="<<fCenterY<<" MVT="
    <<fFromMainVertex<<" Row0="<<fRowRange[0]<<" Row1="<<fRowRange[1]<<" Sector="<<fSector<<" Q="<<fQ<<" TgLam="
    <<fTanl<<" psi="<<fPsi<<" pt="<<fPt<<" L="<<fLength<<" "<<fPterr<<" "<<fPsierr<<" "<<fZ0err<<" "
    <<fTanlerr<<" phi0="<<fPhi0<<" R0="<<fR0<<" Z0="<<fZ0<<" X0="<<fFirstPoint[0]<<" Y0="<<fFirstPoint[1]<<" Z0="
    <<fFirstPoint[2]<<" XL="<<fLastPoint[0]<<" YL="<<fLastPoint[1]<<" ZL="<<fLastPoint[2]<<" "
    <<fPoint[0]<<" "<<fPoint[1]<<" "<<fPoint[2]<<" "<<fPointPsi<<" "<<fIsPoint<<" local="
    <<fIsLocal<<" "<<fPID<<ENDLOG; 

}

int AliHLTTPCTrack::Convert2AliKalmanTrack()
{
  // The method has been copied from AliHLTHoughKalmanTrack and adapted
  // to the TPC conformal mapping track parametrization
  int iResult=0;

  // sector A00 starts at 3 o'clock, sectors are counted counterclockwise
  // median of sector 00 is at 10 degrees, median of sector A04 at 90
  //
  Double_t xhit;
  Double_t charge=(double) GetCharge();
  Double_t xx[5];
  xx[1] = GetFirstPointZ();
  xx[3] = GetTgl();
  xx[4] = charge*(1.0/GetPt());

  Double_t alpha = 0;
  if(GetSector() == -1){

    const double kMaxPhi = TMath::Pi()/2 - 10./180.*TMath::Pi();
    alpha = TMath::ATan2(GetFirstPointY(),GetFirstPointX());
    double phi = GetPsi() - alpha;

    // normalize phi to [-Pi,+Pi]

    phi = phi - TMath::TwoPi() * TMath::Floor( phi /TMath::TwoPi()+.5);

    // extra rotation to keep phi in the range (-Pi/2,+Pi/2)

    double rotation = 0;
    if( phi>=kMaxPhi ) rotation = -TMath::Pi()/2;
    else if( phi<=-kMaxPhi ) rotation = TMath::Pi()/2;

    phi   += rotation;
    alpha -= rotation;
    
    xhit = GetFirstPointX()*TMath::Cos(alpha) + GetFirstPointY()*TMath::Sin(alpha);
    xx[0] = -(GetFirstPointX()*TMath::Sin(alpha)) + GetFirstPointY()*TMath::Cos(alpha);
    xx[2] = TMath::Sin(phi);
  }
  else{
    alpha = fmod((2*GetSector()+1)*(TMath::Pi()/18),2*TMath::Pi());
    if      (alpha < -TMath::Pi()) alpha += 2*TMath::Pi();
    else if (alpha >= TMath::Pi()) alpha -= 2*TMath::Pi();
    
    xhit = GetFirstPointX();
    xx[0] = GetFirstPointY();
    xx[2] = TMath::Sin(GetPsi());
  }
  
  //covariance matrix
  Double_t cov[15]={
    GetY0err(),                         //Error in Y (Y and X are the same)
    0.,  GetZ0err(),                    //Error in Z
    0.,  0.,  GetPsierr(),              //Error for Psi
    0.,  0.,  0.,  GetTglerr(),         //Error for Tgl
    0.,  0.,  0.,  0.,  GetPterr()      //Error for Pt
  };

  Int_t nCluster = GetNHits();
  fdEdx=0;

  // the Set function was not available in earlier versions, check done
  // during configure; for the AliRoot build, by default ON
#ifdef EXTERNALTRACKPARAM_V1
#warning track conversion to ESD format needs AliRoot version > v4-05-04
  //TODO (Feb 07): make this a real warning when logging system is adapted
  //HLTWarning("track conversion to ESD format needs AliRoot version > v4-05-04");
#else
  Set(xhit,alpha,xx,cov);
  SetNumberOfClusters(nCluster);
  SetChi2(0.);
  SetFakeRatio(0.);
  SetMass(0.13957);
#endif

  return iResult;
}

void AliHLTTPCTrack::SetHits(Int_t nhits,UInt_t *hits)
{
  // set hit array
  if (!hits) return;
  if (nhits>fgkHitArraySize) {
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCTrack::SetHits","too many hits")
      << "too many hits (" << nhits << ") for hit array of size " << fgkHitArraySize << ENDLOG; 
  }
  memcpy(fHitNumbers,hits,(nhits<=fgkHitArraySize?nhits:fgkHitArraySize)*sizeof(UInt_t));
}

Double_t AliHLTTPCTrack::GetLengthXY() const
{
  //calculates the length of the arc in XY-plane. This is the length of the track in XY-plane.
  //Using a^2 = b^2 + c^2 - 2bc * cosA for finding the angle between first and last point.
  //Length of arc is arc = r*A. Where A is the angle between first and last point.

  Double_t dx = GetLastPointX()-GetFirstPointX();
  Double_t dy = GetLastPointY()-GetFirstPointY();
  Double_t a = TMath::Sqrt((dx*dx)+(dy*dy)); 
  Double_t r = GetRadius();
  Double_t r2 = r*r;

  Double_t A = TMath::ACos((r2+r2-(a*a))/(2*r2));

  return r*A;
}

Double_t AliHLTTPCTrack::GetLengthTot() const
{
  //Calculates the length of the track in 3D





  return 100.0;

}

int AliHLTTPCTrack::CheckConsistency()
{
  // Check consistency of all members
  int iResult=0;
  if (CheckDoubleMember(&fPterr,   0., "fPterr")<0) iResult=-EDOM;
  if (CheckDoubleMember(&fPsierr,  0., "fPsierr")<0) iResult=-EDOM;
  if (CheckDoubleMember(&fZ0err,   0., "fZ0err")<0) iResult=-EDOM;
  if (CheckDoubleMember(&fY0err,   0., "fY0err")<0) iResult=-EDOM;  
  if (CheckDoubleMember(&fTanlerr, 0., "fTanlerr")<0) iResult=-EDOM;
  return iResult;
}

int AliHLTTPCTrack::CheckDoubleMember(double* pMember, double def, const char* name) const
{
  // Check consistency of a Double member
  if (!pMember) return -EINVAL;
  if (TMath::Abs(*pMember)>kVeryBig) {
    LOG(AliHLTTPCLog::kWarning,"AliHLTTPCTrack","member consistency")
      << "invalid Double number %f" << *pMember << " in member " << name << ENDLOG; 
    *pMember=def;
    return -EDOM;
  }
  return 0;
}
