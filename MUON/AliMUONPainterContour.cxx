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

#include "AliMUONPainterContour.h"

#include "AliMpArea.h"
#include "AliLog.h"
#include <Riostream.h>
#include <TGeoMatrix.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPolyLine.h>
#include <TString.h>
#include <TVector2.h>
#include <TVirtualX.h>
#include <float.h>

///\class AliMUONPainterContour
///
/// Contour for one painter. A contour is a set of TPolyLine (one polyline
/// per closed shape).
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONPainterContour)
///\endcond

//_____________________________________________________________________________
AliMUONPainterContour::AliMUONPainterContour(const char* name) : TNamed(name,""), 
fPolyLines(new TObjArray),
fXmin(FLT_MAX),
fXmax(-FLT_MAX),
fYmin(FLT_MAX),
fYmax(-FLT_MAX)
{
  /// ctor
  fPolyLines->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUONPainterContour::AliMUONPainterContour(const char* name, const AliMpArea& area) 
: TNamed(name,""), 
fPolyLines(new TObjArray),
fXmin(area.LeftBorder()),
fXmax(area.RightBorder()),
fYmin(area.DownBorder()),
fYmax(area.UpBorder())
{
  /// ctor
  fPolyLines->SetOwner(kTRUE);
  TPolyLine* line = new TPolyLine(5);
  Double_t x, y, dx, dy;
  area.GetParameters(x, y, dx, dy);
  line->SetPoint(0,x-dx,y-dy);
  line->SetPoint(1,x-dx,y+dy);
  line->SetPoint(2,x+dx,y+dy);
  line->SetPoint(3,x+dx,y-dy);
  line->SetPoint(4,x-dx,y-dy);
  
  fPolyLines->AddLast(line);
}

//______________________________________________________________________________
AliMUONPainterContour::AliMUONPainterContour(const AliMUONPainterContour& rhs) 
: TNamed(rhs), 
fPolyLines(0x0),
fXmin(FLT_MAX),
fXmax(-FLT_MAX),
fYmin(FLT_MAX),
fYmax(-FLT_MAX)
{
  /// Copy constructor.
  
  ((AliMUONPainterContour&)rhs).Copy(*this);
}

//______________________________________________________________________________
AliMUONPainterContour&
AliMUONPainterContour::operator=(const AliMUONPainterContour& rhs)
{
  /// Assignment operator
  if ( this != &rhs ) 
  {
    delete fPolyLines;
    fPolyLines = 0;
    rhs.Copy(*this);
  }
  return *this;
}

//_____________________________________________________________________________
AliMUONPainterContour::~AliMUONPainterContour()
{
  /// dtor
  delete fPolyLines;
}

//_____________________________________________________________________________
void 
AliMUONPainterContour::AdoptPolyLine(TPolyLine* line)
{
  /// Adopt one polyline into our array of polylines
  fPolyLines->AddLast(line);
  for ( Int_t i = 0; i <= line->GetLastPoint(); ++i ) 
  {
    Double_t x = line->GetX()[i];
    Double_t y = line->GetY()[i];
    fXmin = TMath::Min(fXmin,x);
    fXmax = TMath::Max(fXmax,x);
    fYmin = TMath::Min(fYmin,y);
    fYmax = TMath::Max(fYmax,y);
  }
}

//_____________________________________________________________________________
AliMpArea
AliMUONPainterContour::Area() const
{
  /// Return the area covered by this contour (i.e. the area that
  /// contains all the poylines)
  
  return AliMpArea( ( fXmax+fXmin)/2.0, (fYmax+fYmin)/2.0 ,
                    TMath::Abs(fXmax-fXmin)/2.0, TMath::Abs(fYmax-fYmin)/2.0 );
}

//______________________________________________________________________________
void AliMUONPainterContour::Copy(TObject& obj) const
{
  /// Copy this to obj
  
  AliMUONPainterContour& rhs = static_cast<AliMUONPainterContour&>(obj);
  TNamed::Copy(rhs);
  rhs.fPolyLines = new TObjArray;
  rhs.fPolyLines->SetOwner(kTRUE);
  TIter next(fPolyLines);
  TPolyLine* line;
  while ( ( line = static_cast<TPolyLine*>(next()) ) )
  {
    rhs.fPolyLines->AddLast(line->Clone());
  }
  rhs.fXmin = fXmin;
  rhs.fXmax = fXmax;
  rhs.fYmin = fYmin;
  rhs.fYmax = fYmax;
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterContour::IsInside(Double_t x, Double_t y) const
{
  /// Whether the point (x,y) is inside one of ours polylines
//  if ( x >= fXmin && x <= fXmax && y >= fYmin && y <= fYmax ) 
  {
    TIter next(fPolyLines);
    TPolyLine* line;
    while ( ( line = static_cast<TPolyLine*>(next()) ) )
    {
      if ( TMath::IsInside(x,y,line->Size(),line->GetX(),line->GetY() ) )
      {
        return kTRUE;
      }
    }
  }
  return kFALSE;
}

//_____________________________________________________________________________
void 
AliMUONPainterContour::Offset(const TVector2& offset)
{
  /// Offset all lines by a given offset
  
  TIter next(fPolyLines);
  TPolyLine* line;
  
  while ( ( line = static_cast<TPolyLine*>(next()) ) )
  {
    for ( Int_t i = 0; i <= line->GetLastPoint(); ++i ) 
    {
      Double_t x = line->GetX()[i];
      Double_t y = line->GetY()[i];
      x += offset.X();
      y += offset.Y();
      line->SetPoint(i,x,y);
    }
  }

  fXmin += offset.X();
  fXmax += offset.X();
  fYmin += offset.Y();
  fYmax += offset.Y();
}

//_____________________________________________________________________________
void 
AliMUONPainterContour::PaintArea(Int_t fillColor, Int_t fillStyle)
{
  /// Paint a filled contour
  
  Int_t fc = gVirtualX->GetFillColor();
  Int_t fs = gVirtualX->GetFillStyle();
  
  TIter next(fPolyLines);
  TPolyLine* line;
  
  while ( ( line = static_cast<TPolyLine*>(next()) ) )
  {
    line->SetFillColor(fillColor);
    line->SetFillStyle(fillStyle);
    line->Paint("F");
  }
  
  gVirtualX->SetFillColor(fc);
  gVirtualX->SetFillStyle(fs);
}

//_____________________________________________________________________________
void 
AliMUONPainterContour::PaintOutline(Int_t lineColor, Int_t lineWidth)
{
  /// Paint the outline of this contour
  
  Int_t lc = gVirtualX->GetLineColor();
  Int_t lw = gVirtualX->GetLineWidth();
  
  TIter next(fPolyLines);
  TPolyLine* line;
  
  while ( ( line = static_cast<TPolyLine*>(next()) ) )
  {
    line->SetLineColor(lineColor);
    line->SetLineWidth(lineWidth);
    line->Paint();
  }
  
  gVirtualX->SetLineColor(lc);
  gVirtualX->SetLineWidth(lw);
}

//_____________________________________________________________________________
void 
AliMUONPainterContour::Print(Option_t* opt) const
{
  /// Printout
  
  cout << GetName() << " Ngroups=" << fPolyLines->GetLast()+1;
  TString sopt(opt);
  sopt.ToUpper();

  TIter next(fPolyLines);
  TPolyLine* line;
  while ( ( line = static_cast<TPolyLine*>(next()) ) )
  {
    cout << " (" << line->Size() << ")";
    if ( sopt.Contains("FULL") )
    {
      cout << endl;
      for ( Int_t i = 0; i < line->Size(); ++i ) 
      {
        Double_t x = line->GetX()[i];
        Double_t y = line->GetY()[i];
        cout << Form("Point %3d = %7.3f %7.3f",i,x,y) << endl;
      }
    }
  }
  cout << endl;
}  

//_____________________________________________________________________________
void 
AliMUONPainterContour::Transform(const TGeoHMatrix& matrix)
{
  /// Transform the polylines using the given transformation
  
  TIter next(fPolyLines);
  TPolyLine* line;
  while ( ( line = static_cast<TPolyLine*>(next()) ) )
  {
    for ( Int_t i = 0; i < line->Size(); ++i ) 
    {
      Double_t pl[3] = { line->GetX()[i], line->GetY()[i], 0 };
      Double_t pg[3] = { 0., 0., 0. };
      matrix.LocalToMaster(pl, pg);
      line->SetPoint(i,pg[0],pg[1]);
    }
  }
  

  Double_t pl[3] = { fXmin,fYmin, 0 };
  Double_t pg[3] = { 0., 0., 0. };
  matrix.LocalToMaster(pl, pg);
  
  fXmin = pg[0];
  fYmin = pg[1];
  
  pl[0] = fXmax;
  pl[1] = fYmax;

  matrix.LocalToMaster(pl, pg);
  fXmax = pg[0];
  fYmax= pg[1];
}
