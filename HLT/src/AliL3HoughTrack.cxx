
#include <math.h>

#include "AliL3Transform.h"
#include "AliL3Defs.h"
#include "AliL3HoughTrack.h"

ClassImp(AliL3HoughTrack)


AliL3HoughTrack::AliL3HoughTrack()
{
  //Constructor
  
  fWeight = 0;
  fMinDist=0;
  fTransform = new AliL3Transform();
  fDLine = 0;
  fPsiLine = 0;
  fIsHelix = true;
}


AliL3HoughTrack::~AliL3HoughTrack()
{
  //Destructor
  if(fTransform)
    delete fTransform;
}

void AliL3HoughTrack::Set(AliL3Track *track)
{

  AliL3HoughTrack *tpt = (AliL3HoughTrack*)track;

  fWeight = tpt->GetWeight();
  fDLine = tpt->GetDLine();
  fPsiLine = tpt->GetPsiLine();
  SetNHits(tpt->GetWeight());
  SetRowRange(tpt->GetFirstRow(),tpt->GetLastRow());
  fIsHelix = false;


}

void AliL3HoughTrack::SetTrackParameters(Double_t kappa,Double_t phi,Int_t weight)
{

  fWeight = weight;
  fMinDist = 100000;
  SetKappa(kappa);
  SetPhi0(phi);
  Double_t pt = fabs(BFACT*bField/kappa);
  SetPt(pt);
  Double_t radius = 1/fabs(kappa);
  SetRadius(radius);
  
  //set nhits for sorting.
  SetNHits(weight);
    
  SetCharge(copysign(1,kappa));
  Double_t charge = -1.*kappa;
  Double_t trackPhi0 = GetPhi0() + charge*0.5*Pi/fabs(charge);

  //The first point on track is origo:
  Double_t x0=0;
  Double_t y0=0;

  Double_t xc = x0 - GetRadius() * cos(trackPhi0) ;
  Double_t yc = y0 - GetRadius() * sin(trackPhi0) ;
  SetCenterX(xc);
  SetCenterY(yc);
  fIsHelix = true;
}

void AliL3HoughTrack::SetLineParameters(Double_t psi,Double_t D,Int_t weight,Int_t *rowrange,Int_t ref_row)
{
  //Initialize a track piece, not yet a track
  //Used in case of straight line transformation

  //Transform line parameters to coordinate system of slice:
  
  D = D + fTransform->Row2X(ref_row)*cos(psi);

  fDLine = D;
  fPsiLine = psi;
  fWeight = weight;
  SetNHits(weight);
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

void AliL3HoughTrack::GetLineCrossingPoint(Int_t padrow,Double_t *xy)
{
  

  if(fIsHelix)
    {
      printf("AliL3HoughTrack::GetLineCrossingPoint : Track is not a line\n");
      return;
    }

  Double_t xhit = fTransform->Row2X(padrow);
  Double_t a = -1/tan(fPsiLine);
  Double_t b = fDLine/sin(fPsiLine);
  
  Double_t yhit = a*xhit + b;
  xy[0] = xhit;
  xy[1] = yhit;

}

Double_t AliL3HoughTrack::GetCrossingAngle(Int_t padrow)
{
  //Calculate the crossing angle between track and given padrow.

  if(!fIsHelix)
    {
      printf("AliL3HoughTrack::GetCrossingAngle : Track is not a helix\n");
      return 0;
    }

  if(!IsLocal())
    {
      printf("Track is not given in local coordinates\n");
      return 0;
    }

  Float_t xyz[3];
  if(!GetCrossingPoint(padrow,xyz))
    printf("AliL3HoughTrack::GetCrossingPoint : Track does not cross line!!\n");
  
  
  //Convert center of curvature to local coordinates:
  //Float_t xyz_coc[3] = {GetCenterX(),GetCenterY(),0};
  //fTransform->Global2Local(xyz_coc,slice);
  
  //Take the dot product of the tangent vector of the track, and
  //vector perpendicular to the padrow.
  
  Double_t tangent[2];
  //tangent[1] = (xyz[0] - xyz_coc[0])/GetRadius();
  //tangent[0] = -1.*(xyz[1] - xyz_coc[1])/GetRadius();
  tangent[1] = (xyz[0] - GetCenterX())/GetRadius();
  tangent[0] = -1.*(xyz[1] - GetCenterY())/GetRadius();

  Double_t perp_padrow[2] = {1,0}; //locally in slice

  Double_t cos_beta = fabs(tangent[0]*perp_padrow[0] + tangent[1]*perp_padrow[1]);
  return acos(cos_beta);
  
}

Bool_t AliL3HoughTrack::GetCrossingPoint(Int_t padrow,Float_t *xyz)
{
  //Assumes the track is given in local coordinates

  if(!fIsHelix)
    {
      printf("AliL3HoughTrack::GetCrossingPoint : Track is not a helix\n");
      return 0;
    }
    

  if(!IsLocal())
    {
      printf("GetCrossingPoint: Track is given on global coordinates\n");
      return false;
    }
  
  Double_t xHit = fTransform->Row2X(padrow);

  //xyz[0] = fTransform->Row2X(padrow);
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
  xyz[2] = 0; //only consider transverse plane
  
  return true;
}


Bool_t AliL3HoughTrack::GetCrossingPoint(Int_t slice,Int_t padrow,Float_t *xyz)
{
  //Calculate the crossing point in transverse plane of this track and given 
  //padrow (y = a*x + b). Point is given in local coordinates in slice.
  //Assumes the track is given in global coordinates.

  if(!fIsHelix)
    {
      printf("AliL3HoughTrack::GetCrossingPoint : Track is not a helix\n");
      return 0;
    }

  
  if(IsLocal())
    {
      printf("GetCrossingPoint: Track is given in local coordintes!!!\n");
      return false;
    }

  Double_t padrowradii = fTransform->Row2X(padrow);


  Float_t rotation_angle = (slice*20)*ToRad;
  
  Float_t cs,sn;
  cs = cos(rotation_angle);
  sn = sin(rotation_angle);

  Double_t a = -1.*cs/sn;
  Double_t b = padrowradii/sn;

  Double_t ycPrime = GetCenterY() - b ;
  Double_t aa = ( 1. + a * a ) ;
  Double_t bb = -2. * ( GetCenterX() + a * ycPrime ) ;
  Double_t cc = ( GetCenterX() * GetCenterX() + ycPrime * ycPrime - GetRadius() * GetRadius() ) ;

  Double_t racine = bb * bb - 4. * aa * cc ;
  if ( racine < 0 ) return false ;
  Double_t rootRacine = sqrt(racine) ;
  
  Double_t oneOverA = 1./aa;
//
//   First solution
//
   Double_t x1 = 0.5 * oneOverA * ( -1. * bb + rootRacine ) ; 
   Double_t y1 = a * x1 + b ;
   Double_t r1 = sqrt(x1*x1+y1*y1);
//
//   Second solution
//
   Double_t x2 = 0.5 * oneOverA * ( -1. * bb - rootRacine ) ; 
   Double_t y2 = a * x2 + b ;
   Double_t r2 = sqrt(x2*x2+y2*y2);
//
//    Choose close to (0,0) 
//
   Double_t xHit ;
   Double_t yHit ;
   if ( r1 < r2 ) {
      xHit = x1 ;
      yHit = y1 ;
   }
   else {
      xHit = x2 ;
      yHit = y2 ;
   }
  
   xyz[0] = xHit;
   xyz[1] = yHit;
   xyz[2] = 0;
   
   fTransform->Global2Local(xyz,slice);
   
   return true;
}
