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

///\class AliMUONContour
///
/// A contour is a set of (closed and counter-clockwise-oriented) polygons
///
/// \author Laurent Aphecetche, Subatech

#include "AliMUONContour.h"

#include "AliLog.h"
#include "AliMUONPolygon.h"
#include "AliMpArea.h"
#include <Riostream.h>
#include <TGeoMatrix.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPolyLine.h>
#include <TString.h>
#include <TVector2.h>
#include <float.h>

///\cond CLASSIMP
ClassImp(AliMUONContour)
///\endcond

//_____________________________________________________________________________
AliMUONContour::AliMUONContour(const char* name) : TNamed(name,""), 
fPolygons(new TObjArray),
fXmin(FLT_MAX),
fXmax(-FLT_MAX),
fYmin(FLT_MAX),
fYmax(-FLT_MAX),
fNofVertices(0)
{
  /// ctor
  fPolygons->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUONContour::AliMUONContour(const char* name, const AliMpArea& area) 
: TNamed(name,""), 
fPolygons(new TObjArray),
fXmin(area.LeftBorder()),
fXmax(area.RightBorder()),
fYmin(area.DownBorder()),
fYmax(area.UpBorder()),
fNofVertices(0)
{
  /// ctor
  fPolygons->SetOwner(kTRUE);
  
  AliMUONPolygon* pol = new AliMUONPolygon(area.GetPositionX(),
                                           area.GetPositionY(),
                                           area.GetDimensionX(),
                                           area.GetDimensionY());
  
  fPolygons->AddLast(pol);
  
  fNofVertices = pol->NumberOfVertices();
}

//______________________________________________________________________________
AliMUONContour::AliMUONContour(const AliMUONContour& rhs) 
: TNamed(rhs), 
fPolygons(0x0),
fXmin(FLT_MAX),
fXmax(-FLT_MAX),
fYmin(FLT_MAX),
fYmax(-FLT_MAX),
fNofVertices(0)
{
  /// Copy constructor.
  
  ((AliMUONContour&)rhs).Copy(*this);
}

//______________________________________________________________________________
AliMUONContour&
AliMUONContour::operator=(const AliMUONContour& rhs)
{
  /// Assignment operator
  if ( this != &rhs ) 
  {
    delete fPolygons;
    fPolygons = 0;
    rhs.Copy(*this);
  }
  return *this;
}

//_____________________________________________________________________________
AliMUONContour::~AliMUONContour()
{
  /// dtor
  delete fPolygons;
}

//_____________________________________________________________________________
void 
AliMUONContour::Add(const AliMUONPolygon& polygon)
{
  /// Add points from the polygon
  
  for ( Int_t i = 0; i < polygon.NumberOfVertices(); ++i ) 
  {
    Double_t x = polygon.X(i);
    Double_t y = polygon.Y(i);
    fXmin = TMath::Min(fXmin,x);
    fXmax = TMath::Max(fXmax,x);
    fYmin = TMath::Min(fYmin,y);
    fYmax = TMath::Max(fYmax,y);
  }
  
  fPolygons->AddLast(new AliMUONPolygon(polygon));
  
  fNofVertices += polygon.NumberOfVertices();
}

//_____________________________________________________________________________
AliMpArea
AliMUONContour::Area() const
{
  /// Return the area covered by this contour (i.e. the area that
  /// contains all the poylines)
  
  return AliMpArea( (fXmax+fXmin)/2.0, (fYmax+fYmin)/2.0 ,
                    TMath::Abs(fXmax-fXmin)/2.0, TMath::Abs(fYmax-fYmin)/2.0 );
}

//______________________________________________________________________________
void 
AliMUONContour::AssertOrientation(Bool_t autoCorrect)
{
  /// Insure that all our polygons are counter-clockwise oriented
  /// If autoCorrect==kTRUE, we change the orientation if it is not 
  /// already correct.
  /// If autoCorrect==kFALSE and the orientation is not correct, we
  /// just issue an error message.
  
  for ( Int_t i = 0; i <= fPolygons->GetLast(); ++i )
  {
    AliMUONPolygon* pol = static_cast<AliMUONPolygon*>(fPolygons->UncheckedAt(i));
    if ( !pol->IsCounterClockwiseOriented() ) 
    {
      if ( autoCorrect ) 
      {
        pol->ReverseOrientation();
      }
      else
      {
        AliError("Got a polygon oriented the wrong way");
        StdoutToAliError(Print(););
        return;
      }
    }
  }
}

//______________________________________________________________________________
void AliMUONContour::Copy(TObject& obj) const
{
  /// Copy this to obj
  
  AliMUONContour& rhs = static_cast<AliMUONContour&>(obj);
  TNamed::Copy(rhs);
  delete rhs.fPolygons;
  rhs.fPolygons = new TObjArray(fPolygons->GetLast()+1);
  rhs.fPolygons->SetOwner(kTRUE);
  TIter next(fPolygons);
  AliMUONPolygon* pol;
  while ( ( pol = static_cast<AliMUONPolygon*>(next()) ) )
  {
    rhs.fPolygons->AddLast(pol->Clone());
  }
  rhs.fXmin = fXmin;
  rhs.fXmax = fXmax;
  rhs.fYmin = fYmin;
  rhs.fYmax = fYmax;
  rhs.fNofVertices = fNofVertices;
}

//_____________________________________________________________________________
Bool_t 
AliMUONContour::IsInside(Double_t x, Double_t y) const
{
  /// Whether the point (x,y) is inside one of ours polylines

  if ( x >= fXmin && x <= fXmax && y >= fYmin && y <= fYmax ) 
  {
    TIter next(fPolygons);
    AliMUONPolygon* pol;
    while ( ( pol = static_cast<AliMUONPolygon*>(next()) ) )
    {
      if ( pol->Contains(x,y) ) 
      {
        return kTRUE;
      }
    }      
  }
  
  return kFALSE;
}

//_____________________________________________________________________________
void 
AliMUONContour::Offset(Double_t x, Double_t y)
{
  /// Offset all lines by a given offset
  
  TIter next(fPolygons);
  AliMUONPolygon* pol;
  
  while ( ( pol = static_cast<AliMUONPolygon*>(next()) ) )
  {
    for ( Int_t i = 0; i < pol->NumberOfVertices(); ++i ) 
    {
      pol->SetVertex(i,pol->X(i)+x,pol->Y(i)+y);
    }
  }

  fXmin += x;
  fXmax += x;
  fYmin += y;
  fYmax += y;
}

//_____________________________________________________________________________
void 
AliMUONContour::Print(Option_t* opt) const
{
  /// Printout
  
  cout << GetName() << " NofVertices=" << NumberOfVertices() << " Ngroups=" << fPolygons->GetLast()+1 << endl;
  TString sopt(opt);
  sopt.ToUpper();
  if (sopt.Contains("B"))
  {
    Area().Print("B");
  }

  TIter next(fPolygons);
  AliMUONPolygon* pol;
  while ( ( pol = static_cast<AliMUONPolygon*>(next()) ) )
  {
    pol->Print(opt);
  }
  
  
  cout << endl;
}  

//_____________________________________________________________________________
void 
AliMUONContour::Transform(const TGeoHMatrix& matrix)
{
  /// Transform the polygons using the given transformation
  
  TIter next(fPolygons);
  AliMUONPolygon* pol;
  
  fXmin = fYmin = FLT_MAX;
  fXmax = fYmax = -FLT_MAX;
  
  while ( ( pol = static_cast<AliMUONPolygon*>(next()) ) )
  {
    for ( Int_t i = 0; i < pol->NumberOfVertices(); ++i ) 
    {
      Double_t pl[3] = { pol->X(i), pol->Y(i), 0 };
      Double_t pg[3] = { 0., 0., 0. };
      matrix.LocalToMaster(pl, pg);
      pol->SetVertex(i,pg[0],pg[1]);
      fXmin = TMath::Min(fXmin,pg[0]);
      fYmin = TMath::Min(fYmin,pg[1]);
      fXmax = TMath::Max(fXmax,pg[0]);
      fYmax = TMath::Max(fYmax,pg[1]);
    }
  }
  
  AssertOrientation(kTRUE);
}

//_____________________________________________________________________________
Bool_t 
AliMUONContour::IsValid() const
{
  /// A valid contour is one with a valid area and at least 3 vertices.
  return fNofVertices >= 3 && Area().IsValid();
}
