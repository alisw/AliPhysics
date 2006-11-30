// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>, Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"
#include "AliHLTRootTypes.h"
#include "AliHLTRootTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTTrack.h"
#include "AliHLTTransform.h"
#include "AliHLTVertex.h"
#include "AliHLTSpacePointData.h"

#if __GNUC__ >= 3
using namespace std;
#endif

/** \class AliHLTTrack
//<pre>
//_____________________________________________________________
// AliHLTTrack
//
// Track base class
//Begin_Html
//<img src="track_coordinates.gif">
//End_Html
</pre>
*/

ClassImp(AliHLTTrack)


AliHLTTrack::AliHLTTrack()
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

void AliHLTTrack::Set(AliHLTTrack *tpt)
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

Int_t AliHLTTrack::Compare(const AliHLTTrack *track) const
{
  // compare tracks
  if(track->GetNHits() < GetNHits()) return 1;
  if(track->GetNHits() > GetNHits()) return -1;
  return 0;
}

AliHLTTrack::~AliHLTTrack()
{
  //Nothing to do
}

Double_t AliHLTTrack::GetP() const
{
  // Returns total momentum.  
  return fabs(GetPt())*sqrt(1. + GetTgl()*GetTgl());
}

Double_t AliHLTTrack::GetPseudoRapidity() const
{ //get pseudo rap
  return 0.5 * log((GetP() + GetPz()) / (GetP() - GetPz()));
}

/*
Double_t AliHLTTrack::GetEta() const
{
  return GetPseudoRapidity();
}
*/

Double_t AliHLTTrack::GetRapidity() const
{ 
  //get rap
  const Double_t kmpi = 0.13957;
  return 0.5 * log((kmpi + GetPz()) / (kmpi - GetPz()));
}

void AliHLTTrack::Rotate(Int_t slice,Bool_t tolocal)
{
  //Rotate track to global parameters
  //If flag tolocal is set, the track is rotated
  //to local coordinates.

  Float_t psi[1] = {GetPsi()};
  if(!tolocal)
    AliHLTTransform::Local2GlobalAngle(psi,slice);
  else
    AliHLTTransform::Global2LocalAngle(psi,slice);
  SetPsi(psi[0]);
  Float_t first[3];
  first[0] = GetFirstPointX();
  first[1] = GetFirstPointY();
  first[2] = GetFirstPointZ();
  if(!tolocal)
    AliHLTTransform::Local2Global(first,slice);
  else
    AliHLTTransform::Global2LocHLT(first,slice);
  //AliHLTTransform::Global2Local(first,slice,kTRUE);
  
  SetFirstPoint(first[0],first[1],first[2]);
  Float_t last[3];
  last[0] = GetLastPointX();
  last[1] = GetLastPointY();
  last[2] = GetLastPointZ();
  if(!tolocal)
    AliHLTTransform::Local2Global(last,slice);
  else
    AliHLTTransform::Global2LocHLT(last,slice);    
  //AliHLTTransform::Global2Local(last,slice,kTRUE);
  SetLastPoint(last[0],last[1],last[2]);
  
  Float_t center[3] = {GetCenterX(),GetCenterY(),0};
  if(!tolocal)
    AliHLTTransform::Local2Global(center,slice);
  else
    AliHLTTransform::Global2LocHLT(center,slice);
  //AliHLTTransform::Global2Local(center,slice,kTRUE);
  SetCenterX(center[0]);
  SetCenterY(center[1]);
  
  SetPhi0(atan2(fFirstPoint[1],fFirstPoint[0]));
  SetR0(sqrt(fFirstPoint[0]*fFirstPoint[0]+fFirstPoint[1]*fFirstPoint[1]));
  
  if(!tolocal)
    fIsLocal=kFALSE;
  else
    fIsLocal=kTRUE;
}

void AliHLTTrack::CalculateHelix()
{
  //Calculate Radius, CenterX and CenterY from Psi, X0, Y0
  fRadius = fPt / (AliHLTTransform::GetBFieldValue());
  if(fRadius) fKappa = -fQ*1./fRadius;
  else fRadius = 999999;  //just zero
  Double_t trackPhi0 = fPsi + fQ * AliHLTTransform::PiHalf();

  fCenterX = fFirstPoint[0] - fRadius *  cos(trackPhi0);
  fCenterY = fFirstPoint[1] - fRadius *  sin(trackPhi0);
  
  SetPhi0(atan2(fFirstPoint[1],fFirstPoint[0]));
  SetR0(sqrt(fFirstPoint[0]*fFirstPoint[0]+fFirstPoint[1]*fFirstPoint[1]));
}

Double_t AliHLTTrack::GetCrossingAngle(Int_t padrow,Int_t slice) 
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
      AliHLTTransform::Local2GlobalAngle(&angle,slice);
      if(!CalculateReferencePoint(angle,AliHLTTransform::Row2X(padrow)))
	cerr<<"AliHLTTrack::GetCrossingAngle : Track does not cross line in slice "<<slice<<" row "<<padrow<<endl;
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

Bool_t AliHLTTrack::GetCrossingPoint(Int_t padrow,Float_t *xyz)
{
  //Assumes the track is given in local coordinates
  
  if(!IsLocal())
    {
      cerr<<"GetCrossingPoint: Track is given on global coordinates"<<endl;
      return false;
    }
  
  Double_t xHit = AliHLTTransform::Row2X(padrow);

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
  if(angle1 < 0) angle1 += 2.*AliHLTTransform::Pi();
  Double_t angle2 = atan2((GetFirstPointY() - GetCenterY()),(GetFirstPointX() - GetCenterX()));
  if(angle2 < 0) angle2 += AliHLTTransform::TwoPi();
  Double_t diffangle = angle1 - angle2;
  diffangle = fmod(diffangle,AliHLTTransform::TwoPi());
  if((GetCharge()*diffangle) > 0) diffangle = diffangle - GetCharge()*AliHLTTransform::TwoPi();
  Double_t stot = fabs(diffangle)*GetRadius();
  Double_t zHit = GetFirstPointZ() + stot*GetTgl();
  xyz[2] = zHit;
 
  return true;

}

Bool_t AliHLTTrack::CalculateReferencePoint(Double_t angle,Double_t radius)
{
  // Global coordinate: crossing point with y = ax+ b; 
  // a=tan(angle-AliHLTTransform::PiHalf());
  //
  const Double_t krr=radius; //position of reference plane
  const Double_t kxr = cos(angle) * krr;
  const Double_t kyr = sin(angle) * krr;
  
  Double_t a = tan(angle-AliHLTTransform::PiHalf());
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
  if(fabs(trackPhi0-pointPhi0)>AliHLTTransform::Pi()){
    if(trackPhi0<pointPhi0) trackPhi0 += AliHLTTransform::TwoPi();
    else                    pointPhi0 += AliHLTTransform::TwoPi();
  }
  Double_t stot = -fQ * (pointPhi0-trackPhi0) * fRadius ;
  fPoint[2]   = fFirstPoint[2] + stot * fTanl;

  fPointPsi = pointPhi0 - fQ * AliHLTTransform::PiHalf();
  if(fPointPsi<0.)  fPointPsi+= AliHLTTransform::TwoPi();
  fPointPsi = fmod(fPointPsi, AliHLTTransform::TwoPi());

  return IsPoint(kTRUE);
}

Bool_t AliHLTTrack::CalculateEdgePoint(Double_t angle)
{
  // Global coordinate: crossing point with y = ax; a=tan(angle);
  //
  Double_t rmin=AliHLTTransform::Row2X(AliHLTTransform::GetFirstRow(-1));  //min Radius of TPC
  Double_t rmax=AliHLTTransform::Row2X(AliHLTTransform::GetLastRow(-1)); //max Radius of TPC

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
    if(da<0) da+=AliHLTTransform::TwoPi();
    if(fabs(da-angle)<0.5)
      ok0 = kTRUE;
  }
  if(r1>rmin&&r1<rmax){
    Double_t da=atan2(y1,x1);
    if(da<0) da+=AliHLTTransform::TwoPi();
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
  if(fabs(trackPhi0-pointPhi0)>AliHLTTransform::Pi()){
    if(trackPhi0<pointPhi0) trackPhi0 += AliHLTTransform::TwoPi();
    else                    pointPhi0 += AliHLTTransform::TwoPi();
  }
  Double_t stot = -fQ * (pointPhi0-trackPhi0) * fRadius ;
  fPoint[2]   = fFirstPoint[2] + stot * fTanl;

  fPointPsi = pointPhi0 - fQ * AliHLTTransform::PiHalf();
  if(fPointPsi<0.)  fPointPsi+= AliHLTTransform::TwoPi();
  fPointPsi = fmod(fPointPsi, AliHLTTransform::TwoPi());

  return IsPoint(kTRUE);
}

Bool_t AliHLTTrack::CalculatePoint(Double_t xplane)
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
  if(fabs(trackPhi0-pointPhi0)>AliHLTTransform::Pi()){
    if(trackPhi0<pointPhi0) trackPhi0 += AliHLTTransform::TwoPi();
    else                    pointPhi0 += AliHLTTransform::TwoPi();
  }
  Double_t stot = -fQ * (pointPhi0-trackPhi0) * fRadius ;  
  fPoint[2]   = fFirstPoint[2] + stot * fTanl;

  fPointPsi = pointPhi0 - fQ * AliHLTTransform::PiHalf();
  if(fPointPsi<0.)  fPointPsi+= AliHLTTransform::TwoPi();
  fPointPsi = fmod(fPointPsi, AliHLTTransform::TwoPi());

  return IsPoint(kTRUE);
}

void AliHLTTrack::UpdateToFirstPoint()
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
  pointpsi -= GetCharge()*AliHLTTransform::PiHalf();
  if(pointpsi < 0) pointpsi += AliHLTTransform::TwoPi();
  
  //Update the track parameters
  SetR0(sqrt(point[0]*point[0]+point[1]*point[1]));
  SetPhi0(atan2(point[1],point[0]));
  SetFirstPoint(point[0],point[1],GetZ0());
  SetPsi(pointpsi);
  
}

void AliHLTTrack::GetClosestPoint(AliHLTVertex *vertex,Double_t &closestx,Double_t &closesty,Double_t &closestz)
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
  if(angle1 < 0) angle1 = angle1 + AliHLTTransform::TwoPi();
 
  Double_t angle2 = atan2((GetFirstPointY()-GetCenterY()),(GetFirstPointX()-GetCenterX()));
  if(angle2 < 0) angle2 = angle2 + AliHLTTransform::TwoPi();
  
  Double_t diff_angle = angle1 - angle2;
  diff_angle = fmod(diff_angle,AliHLTTransform::TwoPi());
  
  if((GetCharge()*diff_angle) < 0) diff_angle = diff_angle + GetCharge()*AliHLTTransform::TwoPi();
  Double_t stot = fabs(diff_angle)*GetRadius();
  closestz = GetFirstPointZ() - stot*GetTgl();
}

void AliHLTTrack::Print() const
{ //print out parameters of track
  LOG(AliHLTLog::kInformational,"AliHLTTrack::Print","Print values")
    <<fNHits<<" "<<fMCid<<" "<<fKappa<<" "<<fRadius<<" "<<fCenterX<<" "<<fCenterY<<" "
    <<fFromMainVertex<<" "<<fRowRange[0]<<" "<<fRowRange[1]<<" "<<fSector<<" "<<fQ<<" "
    <<fTanl<<" "<<fPsi<<" "<<fPt<<" "<<fLength<<" "<<fPterr<<" "<<fPsierr<<" "<<fZ0err<<" "
    <<fTanlerr<<" "<<fPhi0<<" "<<fR0<<" "<<fZ0<<" "<<fFirstPoint[0]<<" "<<fFirstPoint[1]<<" "
    <<fFirstPoint[2]<<" "<<fLastPoint[0]<<" "<<fLastPoint[1]<<" "<<fLastPoint[2]<<" "
    <<fPoint[0]<<" "<<fPoint[1]<<" "<<fPoint[2]<<" "<<fPointPsi<<" "<<fIsPoint<<" "
    <<fIsLocal<<" "<<fPID<<ENDLOG; 
}
