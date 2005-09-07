// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCStandardIncludes.h"

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCConfMapPoint.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCConfMapTrack.h"

/**
<pre>
//_____________________________________________________________
// AliHLTTPCConfMapPoint
//
// Hit class for conformal mapper
</pre
*/

ClassImp(AliHLTTPCConfMapPoint)

Bool_t AliHLTTPCConfMapPoint::fDontMap=kFALSE;

AliHLTTPCConfMapPoint::AliHLTTPCConfMapPoint()
{
  //Constructor
  
  SetUsage(false);
  SetHitNumber(-1);
  SetX(0);
  SetY(0);
  SetZ(0);
  SetXerr(0);
  SetYerr(0);
  SetZerr(0);

  SetPhi(0.);
  SetEta(0.);
  
  SetXprime(0.);
  SetYprime(0.);
  SetXprimeerr(0.);
  SetYprimeerr(0.);
  SetIntPoint(0., 0., 0., 0., 0., 0.);
  SetShiftedCoord();
  SetMCTrackID(0,0,0);
}

AliHLTTPCConfMapPoint::~AliHLTTPCConfMapPoint()
{
  // Destructor.
  // Does nothing except destruct. 
}

Bool_t AliHLTTPCConfMapPoint::ReadHits(AliHLTTPCSpacePointData* hits ){
  
  SetHitNumber(hits->fID);
  SetPadRow(hits->fPadRow);
  Int_t slice = (hits->fID>>25) & 0x7f;
  SetSector(slice);
  SetX(hits->fX);
  SetY(hits->fY);
  SetZ(hits->fZ);
  SetXerr(sqrt(hits->fSigmaY2));
  SetYerr(sqrt(hits->fSigmaY2));
  SetZerr(sqrt(hits->fSigmaZ2));
  return kTRUE;
}

void AliHLTTPCConfMapPoint::Reset()
{
  //Reset this point.
  SetUsage(kFALSE);
  SetS(0);
  nextRowHit = 0;
  nextVolumeHit=0;
  nextTrackHit=0;
}

void AliHLTTPCConfMapPoint::Setup(AliHLTTPCVertex *vertex)
{
  //Setup. Sets the vertex, conformal coordinates, and phi and eta of each hit.
  
  SetIntPoint(vertex->GetX(),    vertex->GetY(),    vertex->GetZ(),
              vertex->GetXErr(), vertex->GetYErr(), vertex->GetZErr());
  SetShiftedCoord();
  SetConfCoord();
  // The angles are set properly if they are set after the interaction point and the shifted coordinates
  SetAngles();
  //SetDist(0., 0.);
  
  return;
}

void AliHLTTPCConfMapPoint::SetIntPoint(const Double_t in_x,const Double_t in_y,const Double_t in_z,
			            const Double_t in_x_err,const Double_t in_y_err,const Double_t in_z_err)
{
  // Defines a new interaction point. This point is needed to calculate
  // the conformal coordinates.

  SetXt(in_x);
  SetYt(in_y);
  SetZt(in_z);
  SetXterr(in_x_err);
  SetYterr(in_y_err);
  SetZterr(in_z_err);

  return;
}

void AliHLTTPCConfMapPoint::SetAllCoord(const AliHLTTPCConfMapPoint *preceding_hit)
{
  // Sets the interaction point, the shifted coordinates, and the conformal mapping coordinates.
  // These values are calculated from the interaction point of the given cluster which should be a
  // already found cluster on the same track.

  if (this == preceding_hit) {
    SetIntPoint(preceding_hit->GetX(),    preceding_hit->GetY(),    preceding_hit->GetZ(),
                preceding_hit->GetXerr(), preceding_hit->GetYerr(), preceding_hit->GetZerr());
  }

  else {
    SetIntPoint(preceding_hit->GetXt(),    preceding_hit->GetYt(),    preceding_hit->GetZt(),
                preceding_hit->GetXterr(), preceding_hit->GetYterr(), preceding_hit->GetZterr());
  }

  SetShiftedCoord();
  SetConfCoord();

  return;
}

void AliHLTTPCConfMapPoint::SetShiftedCoord()
{
  // Sets the coordinates with resepct to the given vertex point

  SetXv(GetX() - fXt);
  SetYv(GetY() - fYt);
  SetZv(GetZ() - fZt);
  /*
  SetXverr(TMath::Sqrt(GetXerr()*GetXerr() + fXterr*fXterr));
  SetYverr(TMath::Sqrt(GetYerr()*GetYerr() + fYterr*fYterr));
  SetZverr(TMath::Sqrt(GetZerr()*GetZerr() + fZterr*fZterr));
  */
  return;
}

void AliHLTTPCConfMapPoint::SetConfCoord()
{
  // Calculates the conformal coordinates of one cluster.
  // If the option "vertex_constraint" applies the interaction point is
  // assumed to be at (0, 0, 0). Otherwise the function will use the
  // interaction point specified by fXt and fYt.

  if(fDontMap){
    fXprime = x;
    fYprime = y;
    fWxy = 0;
    s = 0; //track trajectory
    fWz = 0;
    return;
  }

  Double_t r2;
  Double_t xyErrorScale = 1;
  Double_t szErrorScale = 1;

  if ((r2 = fXv*fXv + fYv*fYv)) 
    {
      fXprime =  fXv / r2;
      fYprime = -fYv / r2;
      
      //set weights:
      fWxy = r2*r2 / ((xyErrorScale*xyErrorScale)*((xerr*xerr)+(yerr*yerr)));
      s = 0; //track trajectory
      fWz = (Double_t)(1./(szErrorScale*zerr*zerr));
    } else {
    fXprime    = 0.;
    fYprime    = 0.;
    fXprimeerr = 0.;
    fYprimeerr = 0.;
    fWxy = 0;
    fWz = 0;
    s = 0;
  }

  return;
}

void AliHLTTPCConfMapPoint::SetAngles()
{
  // Calculates the angle phi and the pseudorapidity eta for each cluster.
  /*
  Double_t r = TMath::Sqrt(x*x + y*y);

  fPhi = TMath::ATan2(y,x);
  if(fPhi<0) fPhi = fPhi + 2*TMath::Pi();
  fEta = 3.*z/(TMath::Abs(z)+2.*r);
  return;
  */
  //  Double_t r3dim = TMath::Sqrt(fXv*fXv + fYv*fYv + fZv*fZv);
  Double_t r3dim = sqrt(fXv*fXv + fYv*fYv + fZv*fZv);
  //Double_t r2dim = TMath::Sqrt(fXv*fXv + fYv*fYv);

  /*if (r2dim == 0.) {
  // If r2dim == 0 the pseudorapidity eta cannot be calculated (division by zero)!
  // This can only happen if the point is lying on the z-axis and this should never be possible.
    cerr << "The pseudorapidity eta cannot be calculated (division by zero)! Set to 1.e-10." << endl;
    r2dim = 1.e-10;
  }

  if (fXv == 0.) {
    fPhi = (fYv > 0.) ? TMath::Pi() / 2. : - TMath::Pi() / 2.;
  }

  else {
    fPhi = (fXv > 0.) ? TMath::ASin(fYv/r2dim) : TMath::Pi() - TMath::ASin(fYv/r2dim);
  }

  if (fPhi < 0.) {
    fPhi += 2. * TMath::Pi();
  }
  */
  //fPhi = TMath::ATan2(y,x);
  fPhi = atan2(y,x);
  //if(fPhi<0) fPhi = fPhi + 2*TMath::Pi();
  
  //fEta = 0.5 * TMath::Log((r3dim + fZv)/(r3dim - fZv));
  fEta = 0.5 * log((r3dim + fZv)/(r3dim - fZv));
  return;
}
/*
AliHLTTPCConfMapTrack *AliHLTTPCConfMapPoint::GetTrack(TClonesArray *tracks) const
{
  // Returns the pointer to the track to which this hit belongs.
  
  return (AliHLTTPCConfMapTrack*)tracks->At(this->GetTrackNumber());
}
*/
