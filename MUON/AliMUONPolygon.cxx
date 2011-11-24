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

// $Id$

/// \class AliMUONPolygon
///
/// A simple planar polygon, with a given orientation
///
/// \author Laurent Aphecetche, Subatech

#include "AliMUONPolygon.h"

#include "AliLog.h"
#include "Riostream.h"
#include "TMath.h"

///\cond CLASSIMP
ClassImp(AliMUONPolygon)
///\endcond

//_____________________________________________________________________________
AliMUONPolygon::AliMUONPolygon(Int_t nvertices) 
: TObject(),
fN(nvertices),
fX(new Double_t[fN]),
fY(new Double_t[fN])
{
  /// Ctor with a predefined number of vertices.
}

//_____________________________________________________________________________
AliMUONPolygon::AliMUONPolygon(Double_t xpos, Double_t ypos, Double_t halfsizex, Double_t halfsizey)
: TObject(),
fN(5),
fX(new Double_t[fN]),
fY(new Double_t[fN])
{
  /// Ctor. Polygon will be a rectangle.
  
  
  double xmin(xpos-halfsizex);
  double xmax(xpos+halfsizex);
  double ymin(ypos-halfsizey);
  double ymax(ypos+halfsizey);

  SetVertex(0,xmin,ymin);
  SetVertex(1,xmax,ymin);
  SetVertex(2,xmax,ymax);
  SetVertex(3,xmin,ymax);
  
  Close();
}


//_____________________________________________________________________________
AliMUONPolygon::~AliMUONPolygon()
{
  /// dtor
  delete[] fX;
  delete[] fY;
}

//______________________________________________________________________________
AliMUONPolygon::AliMUONPolygon(const AliMUONPolygon& rhs) 
: TObject(rhs), 
fN(0),
fX(0x0),
fY(0x0)
{
  /// Copy constructor.
  
  ((AliMUONPolygon&)rhs).Copy(*this);
}

//______________________________________________________________________________
AliMUONPolygon&
AliMUONPolygon::operator=(const AliMUONPolygon& rhs)
{
  /// Assignment operator
  if ( this != &rhs ) 
  {
    rhs.Copy(*this);
  }
  return *this;
}

//______________________________________________________________________________
Bool_t
AliMUONPolygon::Contains(Double_t x, Double_t y) const
{
  /// Whether the polygon contains point (x,y)
  
  // Note that the polygon must be a closed polygon (1st and last point
  // must be identical), which should be the case here.

  return TMath::IsInside(x,y,fN,fX,fY);
}

//______________________________________________________________________________
void AliMUONPolygon::Copy(TObject& obj) const
{
  /// Copy this to obj
  
  AliMUONPolygon& rhs = static_cast<AliMUONPolygon&>(obj);

  Double_t* x = new Double_t[fN];
  Double_t* y = new Double_t[fN];
  
  for ( Int_t i = 0; i < fN; ++i )
  {
    x[i] = fX[i];
    y[i] = fY[i];
  }
  
  delete [] rhs.fX;
  delete [] rhs.fY;
  
  rhs.fX = x;
  rhs.fY = y;
  rhs.fN = fN;  
}

//_____________________________________________________________________________
void
AliMUONPolygon::Close()
{
  /// Make that last point = first point
  
  SetVertex(fN-1,X(0),Y(0));
}

//_____________________________________________________________________________
void AliMUONPolygon::Print(Option_t*) const
{
  /// Printout
  cout << Form("AliMUONPolygon : %3d vertices. Signed Area=%e",NumberOfVertices(),SignedArea()) << endl;
  for ( Int_t i = 0; i < NumberOfVertices(); ++i )
  {
    cout << Form("%10.5f,%10.5f",X(i),Y(i)) << endl;
  }
}

//_____________________________________________________________________________
Double_t 
AliMUONPolygon::SignedArea() const
{
  /// Compute the signed area of this polygon
  /// Algorithm from F. Feito, J.C. Torres and A. Urena,
  /// Comput. & Graphics, Vol. 19, pp. 595-600, 1995
  
  Double_t area(0.0);
  
  for ( Int_t i = 0; i < NumberOfVertices()-1; ++i ) 
  {
    area += X(i)*Y(i+1) - X(i+1)*Y(i);
  }
 
  return area;
}

//_____________________________________________________________________________
void 
AliMUONPolygon::ReverseOrientation()
{
  /// Reverse the orientation of this polygon
  Double_t* x = new Double_t[fN];
  Double_t* y = new Double_t[fN];
  
  for ( Int_t i = fN-1; i >= 0; --i )
  {
    x[i] = X(fN-i-1);
    y[i] = Y(fN-i-1);
  }

  delete[] fX;
  delete[] fY;
  
  fX = x;
  fY = y;
}

//_____________________________________________________________________________
void 
AliMUONPolygon::SetVertex(Int_t i, Double_t x, Double_t y)
{
  /// Set one vertex
  if ( i >= fN ) 
  {
    AliFatal("Wrong index");
  }
  fX[i] = x;
  fY[i] = y;
}

