//Author:        Anders Strand Vestbo 
//Author:        Uli Frankenfeld
//Last Modified: 06.03.2001

//____________________________________
// AliL3Track
//
// Base track class for L3

//Changes:

//14.03.01: Moved fHitNumbers from protected to private.-ASV
//          Set memory to zero in ctor.
//          Moved fNHits 2 private. Protected data members not a good idea after all.
//19.03.01: Made the method void Set(AliL3Track) virtual.

#include "AliL3RootTypes.h"

#include "AliL3Logging.h"
#include "AliL3Track.h"
#include <math.h>

ClassImp(AliL3Track)

Float_t AliL3Track::BFACT = 0.0029980;
Float_t AliL3Track::bField = 0.2;
Double_t AliL3Track::pi=3.14159265358979323846;

AliL3Track::AliL3Track()
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
  memset(fHitNumbers,0,174*sizeof(UInt_t));
}

void AliL3Track::Set(AliL3Track *tpt){
  
   SetRowRange(tpt->GetFirstRow(),tpt->GetLastRow());
   SetPhi0(tpt->GetPhi0());
   SetKappa(tpt->GetKappa());
   SetNHits(tpt->GetNHits());

   SetFirstPoint(tpt->GetFirstPointX(),tpt->GetFirstPointY(),tpt->GetFirstPointZ());

  SetLastPoint(tpt->GetLastPointX(),tpt->GetLastPointY(),tpt->GetLastPointZ());
    SetPt(tpt->GetPt());
   SetPsi(tpt->GetPsi());
     SetTgl(tpt->GetTgl());
     SetCharge(tpt->GetCharge());
      
    SetHits(tpt->GetNHits(),(UInt_t *)tpt->GetHitNumbers());

/*
  fPhi0 = track->GetPhi0();
  fKappa = track->GetKappa();
  
  fRowRange[0] = track->GetFirstRow();
  fRowRange[1] = track->GetLastRow();
  fQ = track->GetCharge();
  fFirstPoint[0] = track->GetFirstPointX();
  fFirstPoint[1] = track->GetFirstPointY();
  fFirstPoint[2] = track->GetFirstPointZ();
  fLastPoint[0] = track->GetLastPointX();
  fLastPoint[1] = track->GetLastPointY();
  fLastPoint[2] = track->GetLastPointZ();
  fPt = track->GetPt();
  fTanl = track->GetTgl();
  fPsi = track->GetPsi();
  fQ = track->GetCharge();
  fNHits = track->GetNHits();
  memcpy(fHitNumbers,track->GetHitNumbers(),fNHits*sizeof(UInt_t));
*/
}


AliL3Track::~AliL3Track()
{

}

Double_t AliL3Track::GetP() const
{
  // Returns total momentum.
  
  return fabs(GetPt())*sqrt(1. + GetTgl()*GetTgl());

}

Double_t AliL3Track::GetPseudoRapidity() const
{
  return 0.5 * log((GetP() + GetPz()) / (GetP() - GetPz()));
}

Double_t AliL3Track::GetEta() const
{
  return GetPseudoRapidity();
}

Double_t AliL3Track::GetRapidity() const
{
  Double_t m_pi = 0.13957;
  return 0.5 * log((m_pi + GetPz()) / (m_pi - GetPz()));
}

void AliL3Track::CalculateHelix(){
  //Calculate Radius, CenterX and Centery from Psi, X0, Y0
  //
  
  fRadius = fPt / (BFACT*bField);
  if(fRadius) fKappa = 1./fRadius;
  else fRadius = 999999;  //just zero
  Double_t trackPhi0 = fPsi + fQ *0.5 * pi;

  fCenterX = fFirstPoint[0] - fRadius *  cos(trackPhi0);
  fCenterY = fFirstPoint[1] - fRadius *  sin(trackPhi0);
}

Bool_t AliL3Track::CalculateReferencePoint(Double_t angle){
  // Global coordinate: crossing point with y = ax+ b; a=tan(angle-Pi/2);
  //
  const Double_t rr=132; //position of referece plane
  const Double_t xr = cos(angle) *rr;
  const Double_t yr = sin(angle) *rr;
    
  Double_t a = tan(angle-pi/2);
  Double_t b = yr - a * xr;

  Double_t pp=(fCenterX+a*fCenterY-a*b)/(1+pow(a,2));
  Double_t qq=(pow(fCenterX,2)+pow(fCenterY,2)-2*fCenterY*b+pow(b,2)-pow(fRadius,2))/(1+pow(a,2));

  Double_t racine = pp*pp-qq;
  if(racine<0) return IsPoint(kFALSE);      //no Point

  Double_t rootRacine = sqrt(racine);
  Double_t x0 = pp+rootRacine;
  Double_t x1 = pp-rootRacine;
  Double_t y0 = a*x0 + b;
  Double_t y1 = a*x1 + b;

  Double_t diff0 = sqrt(pow(x0-xr,2)+pow(y0-yr,2));
  Double_t diff1 = sqrt(pow(x1-xr,2)+pow(y1-yr,2));
 
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
  if(fabs(trackPhi0-pointPhi0)>pi){
    if(trackPhi0<pointPhi0) trackPhi0 += 2*pi;
    else                    pointPhi0 += 2*pi;
  }
  Double_t stot = -fQ * (pointPhi0-trackPhi0) * fRadius ;
  fPoint[2]   = fFirstPoint[2] + stot * fTanl;

  fPointPsi = pointPhi0 - fQ * 0.5 * pi;
  if(fPointPsi<0.)  fPointPsi+= 2*pi;
  fPointPsi = fmod(fPointPsi, 2*pi);

  return IsPoint(kTRUE);
}

Bool_t AliL3Track::CalculateEdgePoint(Double_t angle){
  // Global coordinate: crossing point with y = ax; a=tan(angle);
  //
  Double_t rmin=80;  //min Radius of TPC
  Double_t rmax=260; //max Radius of TPC

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
    if(da<0) da+=pi;
    if(fabs(da-angle)<0.5)
      ok0 = kTRUE;
  }
  if(r1>rmin&&r1<rmax){
    Double_t da=atan2(y1,y1);
    if(da<0) da+=pi;
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
  if(fabs(trackPhi0-pointPhi0)>pi){
    if(trackPhi0<pointPhi0) trackPhi0 += 2*pi;
    else                    pointPhi0 += 2*pi;
  }
  Double_t stot = -fQ * (pointPhi0-trackPhi0) * fRadius ;
  fPoint[2]   = fFirstPoint[2] + stot * fTanl;

  fPointPsi = pointPhi0 - fQ * 0.5 * pi;
  if(fPointPsi<0.)  fPointPsi+= 2*pi;
  fPointPsi = fmod(fPointPsi, 2*pi);

  return IsPoint(kTRUE);
}

Bool_t AliL3Track::CalculatePoint(Double_t xplane){
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
  if(fabs(trackPhi0-pointPhi0)>pi){
    if(trackPhi0<pointPhi0) trackPhi0 += 2*pi;
    else                    pointPhi0 += 2*pi;
  }
  Double_t stot = -fQ * (pointPhi0-trackPhi0) * fRadius ;  
  fPoint[2]   = fFirstPoint[2] + stot * fTanl;

  fPointPsi = pointPhi0 - fQ * 0.5 * pi;
  if(fPointPsi<0.)  fPointPsi+= 2*pi;
  fPointPsi = fmod(fPointPsi, 2*pi);

  return IsPoint(kTRUE);
}

