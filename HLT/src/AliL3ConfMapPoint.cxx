//Author: Anders Strand Vestbo

//________________________________
// AliL3ConfMapPoint
//
// Hit class for conformal mapper

#include <iostream.h>
#include <math.h>
#include "AliL3Logging.h"

#include "AliL3ConfMapPoint.h"
#include "AliL3SpacePointData.h"
#include "AliL3Vertex.h"
#include "AliL3ConfMapTrack.h"

//ClassImp(AliL3ConfMapPoint)

AliL3ConfMapPoint::AliL3ConfMapPoint()
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
  SetMCTrackID(0.,0.,0.);
}

AliL3ConfMapPoint::~AliL3ConfMapPoint()
{
  // Destructor.
  // Does nothing except destruct. 
}

Bool_t AliL3ConfMapPoint::ReadHits(AliL3SpacePointData* hits ){
  SetHitNumber(hits->fID);
  SetPadRow(hits->fPadRow);
  Int_t slice = (hits->fID>>25) & 0x7f;
  SetSector(slice);
  SetX(hits->fX);
  SetY(hits->fY);
  SetZ(hits->fZ);
  SetXerr(sqrt(hits->fXYErr));
  SetYerr(sqrt(hits->fXYErr));
  SetZerr(sqrt(hits->fZErr));
  return kTRUE;
}

void AliL3ConfMapPoint::Setup(AliL3Vertex *vertex)
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

void AliL3ConfMapPoint::SetIntPoint(const Double_t in_x,const Double_t in_y, 
			       const Double_t in_z,
			       const Double_t in_x_err, const Double_t in_y_err, 
			       const Double_t in_z_err)
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

void AliL3ConfMapPoint::SetAllCoord(const AliL3ConfMapPoint *preceding_hit)
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

void AliL3ConfMapPoint::SetShiftedCoord()
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

void AliL3ConfMapPoint::SetConfCoord()
{
  // Calculates the conformal coordinates of one cluster.
  // If the option "vertex_constraint" applies the interaction point is
  // assumed to be at (0, 0, 0). Otherwise the function will use the
  // interaction point specified by fXt and fYt.

  Double_t r2;
  Double_t xyErrorScale = 1;
  Double_t szErrorScale = 1;

  if ((r2 = fXv*fXv + fYv*fYv)) 
    {
      fXprime =  fXv / r2;
      fYprime = -fYv / r2;
      //  fXprimeerr = TMath::Sqrt(TMath::Power((-fXv * fXv +   fYv*fYv) * fXverr, 2) + TMath::Power( 2*fXv*fYv*fYverr, 2)) / TMath::Power(fXv*fXv + fYv*fYv, 2);
      // fXprimeerr = TMath::Sqrt(TMath::Power((-fXv * fXv - 3*fYv*fYv) * fYverr, 2) + TMath::Power(-2*fXv*fYv*fXverr, 2)) / TMath::Power(fXv*fXv + fYv*fYv, 2);
    
      
      //set weights:
      //fWxy = r2*r2 / (TMath::Power(xyErrorScale,2)*(TMath::Power(xerr,2)+TMath::Power(yerr,2)));
      fWxy = r2*r2 / ((xyErrorScale*xyErrorScale)*((xerr*xerr)+(yerr*yerr)));
      s = 0; //track trajectory
      //fWz = (Double_t)(1./TMath::Power(szErrorScale*zerr,2));
      fWz = (Double_t)(1./(szErrorScale*zerr*zerr));
    }
  
  else {
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

void AliL3ConfMapPoint::SetAngles()
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
AliL3ConfMapTrack *AliL3ConfMapPoint::GetTrack(TClonesArray *tracks) const
{
  // Returns the pointer to the track to which this hit belongs.
  
  return (AliL3ConfMapTrack*)tracks->At(this->GetTrackNumber());
}
*/
