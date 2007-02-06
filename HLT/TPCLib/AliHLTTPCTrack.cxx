// @(#) $Id$
// Original: AliHLTTrack.cxx,v 1.32 2005/06/14 10:55:21 cvetan 

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          for The ALICE Off-line Project.                               *
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
    @author Anders Vestbo, Uli Frankenfeld, Matthias Richter
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
{
  //Constructor
  fNHits = 0;
  fMCid = -1;
  fKappa=0;
  fRadius=0;
  fCenterX=0;
  fCenterY=0;
  ComesFromMainVertex(false);
  fQ = 0;
  fPhi0=0;
  fPsi=0;
  fR0=0;
  fTanl=0;
  fZ0=0;
  fPt=0;
  fLength=0;
  fIsLocal=true;
  fRowRange[0]=0;
  fRowRange[1]=0;
  SetFirstPoint(0,0,0);
  SetLastPoint(0,0,0);
  memset(fHitNumbers,0,159*sizeof(UInt_t));
  fPID = 0;

  fSector=0;
  fPterr=0;
  fPsierr=0;
  fZ0err=0;
  fTanlerr=0;
  fPoint[0]=fPoint[1]=fPoint[2]=0;
  fPointPsi=0;
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
// #### -B0-CHANGE-START == JMT
    // for straight line fit
    if (AliHLTTPCTransform::GetBFieldValue() == 0.0 ){
	fRadius = 999999;  //just zero

	SetPhi0(atan2(fFirstPoint[1],fFirstPoint[0]));
	SetR0(sqrt(fFirstPoint[0]*fFirstPoint[0]+fFirstPoint[1]*fFirstPoint[1]));
    }
    // for helix fit
    else { 
// #### -B0-UNCHANGED-START == JMT
	//Calculate Radius, CenterX and CenterY from Psi, X0, Y0
	fRadius = fPt / (AliHLTTPCTransform::GetBFieldValue());
	if(fRadius) fKappa = -fQ*1./fRadius;
	else fRadius = 999999;  //just zero
	Double_t trackPhi0 = fPsi + fQ * AliHLTTPCTransform::PiHalf();
	
	fCenterX = fFirstPoint[0] - fRadius *  cos(trackPhi0);
	fCenterY = fFirstPoint[1] - fRadius *  sin(trackPhi0);
	
	SetPhi0(atan2(fFirstPoint[1],fFirstPoint[0]));
	SetR0(sqrt(fFirstPoint[0]*fFirstPoint[0]+fFirstPoint[1]*fFirstPoint[1]));
// #### -B0-UNCHANGED-END == JMT
    }
// #### -B0-CHANGE-END == JMT
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

// BEGINN ############################################## MODIFIY JMT
//if (xHit < xyz[0]){
//    LOG(AliHLTTPCLog::kError,"AliHLTTPCTRACK::GetCrossingPoint","")<< "Track doesn't cross padrow " 
//				<< padrow <<"(x=" << xHit << "). Smallest x=" << xyz[0] << ENDLOG;
//      return false;
//}
// END ################################################# MODIFIY JMT

// #### -B0-CHANGE-START == JMT
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
// #### -B0-UNCHANGED-START == JMT
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
// #### -B0-UNCHANGED-END == JMT
    }
// #### -B0-CHANGE-END == JMT

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
// #### -B0-CHANGE-START == JMT
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
// #### -B0-UNCHANGED-START == JMT
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
// #### -B0-UNCHANGED-END == JMT
    }
// #### -B0-CHANGE-END == JMT
}

void AliHLTTPCTrack::GetClosestPoint(AliHLTTPCVertex *vertex,Double_t &closestx,Double_t &closesty,Double_t &closestz)
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
      closestx = distx1 + vertex->GetX();
      closesty = disty1 + vertex->GetY();
    }
  else
    {
      closestx = distx2 + vertex->GetX();
      closesty = disty2 + vertex->GetY();
    }
  
  //Get the z coordinate:
  Double_t angle1 = atan2((closesty-GetCenterY()),(closestx-GetCenterX()));
  if(angle1 < 0) angle1 = angle1 + AliHLTTPCTransform::TwoPi();
 
  Double_t angle2 = atan2((GetFirstPointY()-GetCenterY()),(GetFirstPointX()-GetCenterX()));
  if(angle2 < 0) angle2 = angle2 + AliHLTTPCTransform::TwoPi();
  
  Double_t diff_angle = angle1 - angle2;
  diff_angle = fmod(diff_angle,AliHLTTPCTransform::TwoPi());
  
  if((GetCharge()*diff_angle) < 0) diff_angle = diff_angle + GetCharge()*AliHLTTPCTransform::TwoPi();
  Double_t stot = fabs(diff_angle)*GetRadius();
  closestz = GetFirstPointZ() - stot*GetTgl();
}

void AliHLTTPCTrack::Print() const
{ //print out parameters of track
// BEGINN ############################################## MODIFIY JMT

#if 1
 LOG(AliHLTTPCLog::kInformational,"AliHLTTPCTrack::Print","Print values")
    <<"NH="<<fNHits<<" "<<fMCid<<" K="<<fKappa<<" R="<<fRadius<<" Cx="<<fCenterX<<" Cy="<<fCenterY<<" MVT="
    <<fFromMainVertex<<" Row0="<<fRowRange[0]<<" Row1="<<fRowRange[1]<<" Sector="<<fSector<<" Q="<<fQ<<" TgLam="
    <<fTanl<<" psi="<<fPsi<<" pt="<<fPt<<" L="<<fLength<<" "<<fPterr<<" "<<fPsierr<<" "<<fZ0err<<" "
    <<fTanlerr<<" phi0="<<fPhi0<<" R0="<<fR0<<" Z0"<<fZ0<<" X0"<<fFirstPoint[0]<<" Y0"<<fFirstPoint[1]<<" Z0"
    <<fFirstPoint[2]<<" XL"<<fLastPoint[0]<<" YL"<<fLastPoint[1]<<" ZL"<<fLastPoint[2]<<" "
    <<fPoint[0]<<" "<<fPoint[1]<<" "<<fPoint[2]<<" "<<fPointPsi<<" "<<fIsPoint<<" local="
    <<fIsLocal<<" "<<fPID<<ENDLOG; 



#else
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCTrack::Print","Print values")
    <<fNHits<<" "<<fMCid<<" "<<fKappa<<" "<<fRadius<<" "<<fCenterX<<" "<<fCenterY<<" "
    <<fFromMainVertex<<" "<<fRowRange[0]<<" "<<fRowRange[1]<<" "<<fSector<<" "<<fQ<<" "
    <<fTanl<<" "<<fPsi<<" "<<fPt<<" "<<fLength<<" "<<fPterr<<" "<<fPsierr<<" "<<fZ0err<<" "
    <<fTanlerr<<" "<<fPhi0<<" "<<fR0<<" "<<fZ0<<" "<<fFirstPoint[0]<<" "<<fFirstPoint[1]<<" "
    <<fFirstPoint[2]<<" "<<fLastPoint[0]<<" "<<fLastPoint[1]<<" "<<fLastPoint[2]<<" "
    <<fPoint[0]<<" "<<fPoint[1]<<" "<<fPoint[2]<<" "<<fPointPsi<<" "<<fIsPoint<<" "
    <<fIsLocal<<" "<<fPID<<ENDLOG; 
#endif

// END ################################################# MODIFIY JMT
}

int AliHLTTPCTrack::Convert2AliKalmanTrack()
{
  int iResult=0;
  // The method has been copied from AliHLTHoughKalmanTrack and adapted
  // to the TPC conformal mapping track parametrization

  SetChi2(0.);
  SetNumberOfClusters(GetLastRow()-GetFirstRow());
  SetLabel(GetMCid());
  SetFakeRatio(0.);
  SetMass(0.13957); // just a guess

  fdEdx=0;
  Double_t alpha = fmod((GetSector()+0.5)*(2*TMath::Pi()/18),2*TMath::Pi());
  if      (alpha < -TMath::Pi()) alpha += 2*TMath::Pi();
  else if (alpha >= TMath::Pi()) alpha -= 2*TMath::Pi();

  Double_t xhit=GetFirstPointX();
  Double_t yhit=GetFirstPointY();
  Double_t zhit=GetFirstPointZ();
  Double_t psi = GetPsi();
  Double_t kappa = GetKappa();
  Double_t radius = GetRadius();
  Double_t centerx = GetCenterX();

  Double_t tanl = GetTgl();

  Double_t cnv=1.;
  // TODO: think about how to get the magnetic field
  //Double_t cnv=1./(GetBz()*kB2C);

  //covariance matrix
  Double_t cov[15]={
    0.,
    0.,  0.,
    0.,  0.,  0.,
    0.,  0.,  0.,  0.,
    0.,  0.,  0.,  0.,  0.
  };

  Double_t xx[5];
  xx[0] = yhit;
  xx[1] = zhit;
  xx[2] = (xhit-centerx)/radius;
  xx[3] = tanl;
  xx[4] = kappa*cnv;
  // the Set function was not available in earlier versions, check required in
  // configure.ac
  //Set(xhit,alpha,xx,cov);

  return iResult;
}
