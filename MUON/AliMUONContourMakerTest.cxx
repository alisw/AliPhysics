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

///
/// \class AliMUONContourMakerTest
/// 
/// Class used to test (and in particular time) the contour creation
/// algorithms.
///
/// \author Laurent Aphecetche, Subatech
///

#include "AliMUONContourMakerTest.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONContour.h"
#include "AliMUONPolygon.h"
#include "AliMUONSegment.h"
#include "AliMUONContourHandler.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpExMap.h"
#include <float.h>
#include "Riostream.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGeoMatrix.h"
#include "TLine.h"
#include "TMap.h"
#include "TObjArray.h"
#include "TPolyLine.h"
#include "TSystem.h"

///\cond CLASSIMP
ClassImp(AliMUONContourMakerTest)
///\endcond 

namespace
{
  //_____________________________________________________________________________
  void Plot(TPolyLine& line, Bool_t orientation)
  {
    if ( !orientation ) 
    {
      line.Draw();
    }
    else
    {
      Double_t* x = line.GetX();
      Double_t* y = line.GetY();
      
      for ( Int_t i = 0; i < line.GetLastPoint(); ++i ) 
      {
        Double_t x1 = x[i];
        Double_t y1 = y[i];
        Double_t x2 = x[i+1];
        Double_t y2 = y[i+1];
        
        Bool_t horizontal = AliMUONSegment::AreEqual(y1,y2);
        
        TLine* a = new TArrow(x1,y1,x2,y2,0.03,"->-");
        if (horizontal)
        {
          a->SetLineStyle(3);
        }
        a->Draw();
      }
    }
  }
}

//_____________________________________________________________________________
AliMUONContourMakerTest::AliMUONContourMakerTest()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONContourMakerTest::~AliMUONContourMakerTest()
{
  /// dtor
}


//_____________________________________________________________________________
void 
AliMUONContourMakerTest::Exec(const Option_t* opt)
{
  /// Main method
  /// Generate the geometry transformations, then
  /// contours for all manus, and then for all the elements
  /// (bus patches, detection elements, etc...)
  
  AliInfo("Resetting all timers before I start...");
  
  AliCodeTimer::Instance()->Reset();
  
  AliMpCDB::LoadDDLStore2();
  
  AliCodeTimer::Instance()->Print();
  
  AliInfo("Resetting all timers after loading the mapping...");
  
  AliCodeTimer::Instance()->Reset();
  
  AliCodeTimerAuto("",0);

  TString sopt(opt);
  
  Bool_t explodedView(kTRUE);
  
  if (sopt.Contains("REAL")) explodedView = kFALSE;
    
  AliMUONContourHandler ch(explodedView);
  
  if ( sopt.Contains("SAVE") )
  {
    TFile f2("AliMUONContourMakerTest.manuContours.root","RECREATE");
    ch.AllContourMap()->Write("ALL",TObject::kSingleKey);
    f2.Close();
  }

  AliCodeTimer::Instance()->Print();  
}


//_____________________________________________________________________________
void 
AliMUONContourMakerTest::GetBoundingBox(const TObjArray& array, 
                                        Double_t& xmin, Double_t& ymin, 
                                        Double_t& xmax, Double_t& ymax,
                                        Bool_t enlarge) const
{
  /// Get the bounding box of all the contours in array. 
  /// If enlarge = kTRUE, the bounding box is "enlarged" a bit
  /// (e.g. to leave some blank around a plot in a canvas)
  ///
  
  xmin=ymin=FLT_MAX;
  xmax=ymax=-FLT_MAX;
  TIter next(&array);
  AliMUONContour* contour;
  while ( ( contour = static_cast<AliMUONContour*>(next()) ) )
  {
    AliMpArea area(contour->Area());
    xmin = TMath::Min(xmin,area.LeftBorder());
    xmax = TMath::Max(xmax,area.RightBorder());
    ymin = TMath::Min(ymin,area.DownBorder());
    ymax = TMath::Max(ymax,area.UpBorder());
  }

  if (enlarge)
  {
    Double_t xsize = (xmax-xmin);
    Double_t ysize = (ymax-ymin);
    Double_t xshift = xsize*0.1;
    Double_t yshift = ysize*0.1;    
    xmin -= xshift;
    ymin -= yshift;
    xmax = xmin + xsize + xshift*2;
    ymax = ymin + ysize + yshift*2;
  }
}

//_____________________________________________________________________________
void
AliMUONContourMakerTest::PlotSegments(const TObjArray& segments, Int_t lineColor, Int_t lineWidth, Bool_t orientation) const
{
  /// Plot an array of segments 
  
  TIter next(&segments);
  AliMUONSegment* s;
  while ( ( s = static_cast<AliMUONSegment*>(next()) ) )
  {
    TPolyLine* line = new TPolyLine(2);
    line->SetPoint(0,s->StartX(),s->StartY());
    line->SetPoint(1,s->EndX(),s->EndY());
    line->SetLineColor(lineColor);
    line->SetLineWidth(lineWidth);
    ::Plot(*line,orientation);
  }
}

//_____________________________________________________________________________
void 
AliMUONContourMakerTest::Plot(const AliMUONPolygon& polygon, 
                              Int_t lineColor, Int_t lineWidth,
                              Bool_t orientation) const 
{
  /// Plot a polygon
  TPolyLine* line = new TPolyLine(polygon.NumberOfVertices());
  for ( Int_t i = 0; i < polygon.NumberOfVertices(); ++i )
  {
    line->SetPoint(i,polygon.X(i),polygon.Y(i));
  }
  
  line->SetLineColor(lineColor);
  line->SetLineWidth(lineWidth);
  ::Plot(*line,kFALSE);
  if ( orientation ) ::Plot(*line,kTRUE);
}

//_____________________________________________________________________________
void 
AliMUONContourMakerTest::Plot(const AliMUONContour& contour, Int_t lineColor, Int_t lineWidth,
                              Bool_t orientation) const 
{
  /// Plot a contour (i.e. a set of polygons)
  const TObjArray* polygons = contour.Polygons();
  TIter next(polygons);
  AliMUONPolygon* pol;
  while ( ( pol = static_cast<AliMUONPolygon*>(next()) ) )
  {
    Plot(*pol,lineColor,lineWidth,orientation);
  }
}

//_____________________________________________________________________________
void 
AliMUONContourMakerTest::PlotContours(const TObjArray& array, Bool_t orientations) const
{
  /// Plot an array of contours
  TIter next(&array);
  AliMUONContour* contour;
  while ( ( contour = static_cast<AliMUONContour*>(next()) ) )
  {
    Plot(*contour,5,4,orientations);
  }
}

//______________________________________________________________________________
void 
AliMUONContourMakerTest::PrintAsPNG(const char* basename, const TObjArray& contourArray,
                                    const TObjArray* verticals, const TObjArray* horizontals) const
{
  /// Output contours and segments into a PNG file.
  TCanvas* c = new TCanvas(basename,basename,0,0,600,600);
  double xmin,ymin,xmax,ymax;
  GetBoundingBox(contourArray,xmin,ymin,xmax,ymax,kTRUE);
  c->Range(xmin,ymin,xmax,ymax);
  PlotContours(contourArray,kTRUE);
  c->Modified();
  c->Update();
  TString name(Form("%s",basename));
  name.ReplaceAll("/","_");
  c->Print(Form("%s.png",name.Data()));
  if ( verticals || horizontals ) 
  {
    c->Clear();
    if ( verticals ) PlotSegments(*verticals,1);
    if ( horizontals) PlotSegments(*horizontals,2);
    c->Modified();
    c->Update();
    name = Form("%s",basename);
    name.ReplaceAll("/","_");
    c->Print(Form("%s-segments.png",name.Data()));
  }
  delete c;
}

