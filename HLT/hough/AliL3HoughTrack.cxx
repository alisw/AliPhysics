// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3HoughTrack.h"
#include "AliL3Transform.h"

#if __GNUC__ == 3
using namespace std;
#endif

/** \class AliL3HoughTrack
<pre>
//_____________________________________________________________
// AliL3HoughTrack
//
// Track class for Hough tracklets
//
</pre>
*/

ClassImp(AliL3HoughTrack)


AliL3HoughTrack::AliL3HoughTrack()
{
  //Constructor
  
  fWeight = 0;
  fMinDist=0;
  fDLine = 0;
  fPsiLine = 0;
  fIsHelix = true;
  fEtaIndex = -1;
  fEta = 0;
}

AliL3HoughTrack::~AliL3HoughTrack()
{
}

void AliL3HoughTrack::Set(AliL3Track *track)
{
  AliL3HoughTrack *tpt = (AliL3HoughTrack*)track;
  SetTrackParameters(tpt->GetKappa(),tpt->GetPsi(),tpt->GetWeight());
  SetEtaIndex(tpt->GetEtaIndex());
  SetEta(tpt->GetEta());
  SetTgl(tpt->GetTgl());
  SetPsi(tpt->GetPsi());
  SetCenterX(tpt->GetCenterX());
  SetCenterY(tpt->GetCenterY());
  SetFirstPoint(tpt->GetFirstPointX(),tpt->GetFirstPointY(),tpt->GetFirstPointZ());
  SetLastPoint(tpt->GetLastPointX(),tpt->GetLastPointY(),tpt->GetLastPointZ());
  SetCharge(tpt->GetCharge());
  SetRowRange(tpt->GetFirstRow(),tpt->GetLastRow());
  SetSlice(tpt->GetSlice());
  SetHits(tpt->GetNHits(),(UInt_t *)tpt->GetHitNumbers());
#ifdef ROWHOUGH
  SetMCid(tpt->GetMCid());
#endif
  /*
    fWeight = tpt->GetWeight();
    fDLine = tpt->GetDLine();
    fPsiLine = tpt->GetPsiLine();
    SetNHits(tpt->GetWeight());
    SetRowRange(tpt->GetFirstRow(),tpt->GetLastRow());
    fIsHelix = false;
  */
}

Int_t AliL3HoughTrack::Compare(const AliL3Track *tpt) const
{
  AliL3HoughTrack *track = (AliL3HoughTrack*)tpt;
  if(track->GetWeight() < GetWeight()) return 1;
  if(track->GetWeight() > GetWeight()) return -1;
  return 0;
}

void AliL3HoughTrack::SetEta(Double_t f)
{
  //Set eta, and calculate fTanl, which is the tan of dipangle

  fEta = f;
  Double_t theta = 2*atan(exp(-1.*fEta));
  Double_t dipangle = AliL3Transform::PiHalf() - theta;
  Double_t tgl = tan(dipangle);
  SetTgl(tgl);
}

void AliL3HoughTrack::UpdateToFirstRow()
{
  //Update the track parameters to the point where track cross
  //its first padrow.`
  
  //Get the crossing point with the first padrow:
  Float_t xyz[3];
  if(!GetCrossingPoint(GetFirstRow(),xyz))
    LOG(AliL3Log::kWarning,"AliL3HoughTrack::UpdateToFirstRow()","Track parameters")
      <<AliL3Log::kDec<<"Track does not cross padrow "<<GetFirstRow()<<" centerx "
      <<GetCenterX()<<" centery "<<GetCenterY()<<" Radius "<<GetRadius()<<" tgl "<<GetTgl()<<ENDLOG;
  
  //printf("Track with eta %f tgl %f crosses at x %f y %f z %f on padrow %d\n",GetEta(),GetTgl(),xyz[0],xyz[1],xyz[2],GetFirstRow());
  //printf("Before: first %f %f %f tgl %f center %f %f charge %d\n",GetFirstPointX(),GetFirstPointY(),GetFirstPointZ(),GetTgl(),GetCenterX(),GetCenterY(),GetCharge());
  
  Double_t radius = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);

  //Get the track parameters
  
  /*
    Double_t x0    = GetR0() * cos(GetPhi0()) ;
    Double_t y0    = GetR0() * sin(GetPhi0()) ;
  */
  Double_t rc    = GetRadius();//fabs(GetPt()) / AliL3Transform::GetBFieldValue();
  Double_t tPhi0 = GetPsi() + GetCharge() * AliL3Transform::PiHalf() / abs(GetCharge()) ;
  Double_t xc    = GetCenterX();//x0 - rc * cos(tPhi0) ;
  Double_t yc    = GetCenterY();//y0 - rc * sin(tPhi0) ;
  
  //Check helix and cylinder intersect
  Double_t fac1 = xc*xc + yc*yc ;
  Double_t sfac = sqrt( fac1 ) ;
  
  if ( fabs(sfac-rc) > radius || fabs(sfac+rc) < radius ) {
    LOG(AliL3Log::kError,"AliL3HoughTrack::UpdateToFirstRow","Tracks")<<AliL3Log::kDec<<
      "Track does not intersect"<<ENDLOG;
    return;
  }
  
  //Find intersection
  Double_t fac2   = (radius*radius + fac1 - rc*rc) / (2.00 * radius * sfac ) ;
  Double_t phi    = atan2(yc,xc) + GetCharge()*acos(fac2) ;
  Double_t td     = atan2(radius*sin(phi) - yc,radius*cos(phi) - xc) ;
  
  //Intersection in z
  if ( td < 0 ) td = td + AliL3Transform::TwoPi();
  Double_t deltat = fmod((-GetCharge()*td + GetCharge()*tPhi0),AliL3Transform::TwoPi());
  if ( deltat < 0. ) deltat += AliL3Transform::TwoPi();
  else if ( deltat > AliL3Transform::TwoPi() ) deltat -= AliL3Transform::TwoPi();
  Double_t z = GetZ0() + rc * GetTgl() * deltat ;
  
  Double_t xExtra = radius * cos(phi) ;
  Double_t yExtra = radius * sin(phi) ;
  
  Double_t tPhi = atan2(yExtra-yc,xExtra-xc);
  
  //if ( tPhi < 0 ) tPhi += 2. * M_PI ;
  Double_t tPsi = tPhi - GetCharge() * AliL3Transform::PiHalf() / abs(GetCharge()) ;
  if ( tPsi > AliL3Transform::TwoPi() ) tPsi -= AliL3Transform::TwoPi() ;
  else if ( tPsi < 0. ) tPsi += AliL3Transform::TwoPi();
  
  //And finally, update the track parameters
  SetR0(radius);
  SetPhi0(phi);
  SetZ0(z);
  SetPsi(tPsi);
  SetFirstPoint(xyz[0],xyz[1],z);
  //printf("After: first %f %f %f tgl %f center %f %f charge %d\n",GetFirstPointX(),GetFirstPointY(),GetFirstPointZ(),GetTgl(),GetCenterX(),GetCenterY(),GetCharge());
  
  //printf("First point set %f %f %f\n",xyz[0],xyz[1],z);
  
  //Also, set the coordinates of the point where track crosses last padrow:
  GetCrossingPoint(GetLastRow(),xyz);
  SetLastPoint(xyz[0],xyz[1],xyz[2]);
  //printf("last point %f %f %f\n",xyz[0],xyz[1],xyz[2]);
}

void AliL3HoughTrack::SetTrackParameters(Double_t kappa,Double_t eangle,Int_t weight)
{
  fWeight = weight;
  fMinDist = 100000;
  SetKappa(kappa);
  Double_t pt = fabs(AliL3Transform::GetBFieldValue()/kappa);
  SetPt(pt);
  Double_t radius = 1/fabs(kappa);
  SetRadius(radius);
  SetFirstPoint(0,0,0);
  SetPsi(eangle); //Psi = emission angle when first point is vertex
  SetPhi0(0);     //not defined for vertex reference point
  SetR0(0);
  Double_t charge = -1.*kappa;
  SetCharge((Int_t)copysign(1.,charge));
  Double_t trackPhi0 = GetPsi() + charge*0.5*AliL3Transform::Pi()/fabs(charge);
  Double_t xc = GetFirstPointX() - GetRadius() * cos(trackPhi0) ;
  Double_t yc = GetFirstPointY() - GetRadius() * sin(trackPhi0) ;
  SetCenterX(xc);
  SetCenterY(yc);
  SetNHits(1); //just for the trackarray IO
  fIsHelix = true;
}

void AliL3HoughTrack::SetLineParameters(Double_t psi,Double_t D,Int_t weight,Int_t *rowrange,Int_t /*ref_row*/)
{
  //Initialize a track piece, not yet a track
  //Used in case of straight line transformation

  //Transform line parameters to coordinate system of slice:
  
  //  D = D + fTransform->Row2X(ref_row)*cos(psi);

  fDLine = D;
  fPsiLine = psi;
  fWeight = weight;
  SetNHits(1);
  SetRowRange(rowrange[0],rowrange[1]);
  fIsHelix = false;
}

void AliL3HoughTrack::SetBestMCid(Int_t mcid,Double_t min_dist)
{
  if(min_dist < fMinDist)
    {
      fMinDist = min_dist;
      SetMCid(mcid);
    }
}

void AliL3HoughTrack::GetLineCrossingPoint(Int_t padrow,Float_t *xy)
{
  if(fIsHelix)
    {
      printf("AliL3HoughTrack::GetLineCrossingPoint : Track is not a line\n");
      return;
    }

  Float_t xhit = AliL3Transform::Row2X(padrow) - AliL3Transform::Row2X(GetFirstRow());
  Float_t a = -1/tan(fPsiLine);
  Float_t b = fDLine/sin(fPsiLine);
  Float_t yhit = a*xhit + b;
  xy[0] = xhit;
  xy[1] = yhit;
}
