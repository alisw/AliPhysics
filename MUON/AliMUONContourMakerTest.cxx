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
#include "AliMUONContourMaker.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONManuContourMaker.h"
#include "AliMUONPolygon.h"
#include "AliMUONSegment.h"
#include "AliMpArea.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpExMap.h"
#include <float.h>
#include "Riostream.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGeoMatrix.h"
#include "TLine.h"
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


//______________________________________________________________________________
TObjArray*
AliMUONContourMakerTest::CreateContourList(const TObjArray& manuContours)
{    
  /// Create an array of maps of contour names
  ///
  /// Assyming that key is something like station#/chamber#/de#/buspatch#/manu#
  /// the idea here is to put one TMap for each level in mapArray :
  ///
  /// mapArray[0].key = station0
  /// mapArray[0].value = map of strings { station0/chamber0, station0/chamber1 }
  ///
  /// Then each entry in mapArray will be converted into a contour by
  /// merging its children (e.g. station0 contour will be made from the merging
  /// of station0/chamber0 and station0/chamber1 in the example above).
  ///
  
  AliCodeTimerAuto("");
  
  TIter next(&manuContours);
  AliMUONContour* contour;
  TObjArray* mapArray = new TObjArray;

  while ( ( contour = static_cast<AliMUONContour*>(next()) ) )
  {
    // Key is something like station#/chamber#/de#/buspatch#/manu#

    TString key(contour->GetName());
    TObjArray* s = key.Tokenize("/");
    for ( Int_t i = 0; i < s->GetLast(); ++i ) 
    {
      TMap* m = static_cast<TMap*>(mapArray->At(i));
      if (!m)
      {
        m = new TMap;
        if ( i > mapArray->GetSize() ) mapArray->Expand(i);
          mapArray->AddAt(m,i);
          }
      TString parent;
      for ( Int_t k = 0; k <= i; ++k )
      {
        TObjString* str = static_cast<TObjString*>(s->At(k));
        parent += str->String();
        if ( k < i ) parent += "/";
          }
      TString child(parent);
      child += "/";
      child += static_cast<TObjString*>(s->At(i+1))->String();
      
      TObjArray* ma = static_cast<TObjArray*>(m->GetValue(parent.Data()));
      if (!ma)
      {
        ma = new TObjArray;
        m->Add(new TObjString(parent.Data()),ma);
      }
      TPair* p = static_cast<TPair*>(ma->FindObject(child.Data()));
      if ( !p ) 
      {
        ma->Add(new TObjString(child.Data()));
      }
    }
    delete s;
  }
 
  return mapArray;
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
  
  AliCodeTimerAuto("");
  
  AliMpExMap* real(0x0);
  AliMpExMap* exploded(0x0);
  
  GenerateTransformations(real,exploded);
  
  TObjArray* manus(0x0);
  TObjArray* all(0x0);
  
  TString sopt(opt);
  
  if ( sopt.Contains("MANU") || sopt.Contains("ALL") ) 
  {
    AliMUONManuContourMaker manuMaker(exploded);
    manus = manuMaker.GenerateManuContours(kTRUE);
  }
  
  if ( sopt.Contains("ALL") && manus )
  {
    manus->SetOwner(kFALSE);
    all = GenerateAllContours(*manus);
    if ( sopt.Contains("SAVE") && all )
    {
      TFile f("AliMUONContourMakerTest.all.root","RECREATE");
      all->Write("ALL",TObject::kSingleKey);
      f.Close();
    }
    
  }
  
  AliCodeTimer::Instance()->Print();
  
  delete manus;
  delete all;
}

//______________________________________________________________________________
TObjArray* 
AliMUONContourMakerTest::GenerateAllContours(const TObjArray& manuContours)
{
  /// From a map of manu contours, generate the compound contours (bp, de, etc...)
  /// by merging them.
  /// Note that manuContours should NOT be the owner of its contours,
  /// as they are adopted by the array returned by this method.
  
  AliCodeTimerAuto("");
  
  // Get the list of contours to create
  TObjArray* mapArray = CreateContourList(manuContours);
    
  // Now loop over the mapArray to actually create the contours
  TIter next2(mapArray,kIterBackward);

  TMap allContourMap;
  allContourMap.SetOwnerKeyValue(kTRUE,kFALSE); // not owner of contours, as the returned array will be the owner
  TObjArray* allContourArray = new TObjArray;
  allContourArray->SetOwner(kTRUE);
  
  TIter nextContour(&manuContours);  
  AliMUONContour* contour(0x0);
  
  while ( ( contour = static_cast<AliMUONContour*>(nextContour()) ) )
  {
    allContourMap.Add(new TObjString(contour->GetName()),contour);
    allContourArray->Add(contour);
  }
  
  AliMUONContourMaker maker;
  
  for ( Int_t i = mapArray->GetLast(); i >= 1; --i ) 
    // end at 1 to avoid merging different cathodes together, which
    // would not work...
  {
    TMap* a = static_cast<TMap*>(mapArray->At(i));
    TIter next3(a);
    TObjString* str;
    while ( ( str = static_cast<TObjString*>(next3()) ) )
    {
      TObjArray* m = static_cast<TObjArray*>(a->GetValue(str->String().Data()));
      TIter next4(m);
      TObjString* k;
      TObjArray subcontours;
      subcontours.SetOwner(kFALSE);
      while ( ( k = static_cast<TObjString*>(next4()) ) )
      {
        contour = static_cast<AliMUONContour*>(allContourMap.GetValue(k->String().Data()));
        if ( contour ) 
        {
          subcontours.Add(contour);
        }
        else
        {
          AliError(Form("Did not find contour %s",k->String().Data()))
          return allContourArray;
        }
      }

      contour = maker.MergeContour(subcontours,str->String().Data());
        
      bool error(kFALSE);
      
      if (!contour)
      {
        error=kTRUE;
        AliError(Form("ERROR : could not merge into %s",str->String().Data()));
      }
      else
      {
        if ( contour->Area().IsValid() == kFALSE ) 
        {
          error=kTRUE;
          AliError(Form("ERROR : area of contour %s is invalid",str->String().Data()));
        }
      }
      
      if ( error ) 
      {
        // do it again, but get intermediate results to plot them
        PrintAsPNG(str->String().Data(),subcontours);
        if (contour ) 
        {
          StdoutToAliError(contour->Area().Print("B"););
        }
        AliError(Form("%d subcontours",subcontours.GetLast()+1));
        StdoutToAliError(subcontours.Print(););
        // check whether one of the subcontour itself is already invalid ?
        TIter next(&subcontours);
        AliMUONContour* cont;
        while ( ( cont = static_cast<AliMUONContour*>(next()) ) )
        {
          if (!cont->IsValid())
          {
            AliError(Form("subcontour %s is invalid",cont->GetName()));
          }
        }
        TFile f("subcontour.root","recreate");
        subcontours.Write("fault",TObject::kSingleKey);
        f.Close();
                
        return allContourArray;
      }
      
      allContourArray->Add(contour);
      allContourMap.Add(new TObjString(str->String().Data()),contour);
    }
  }
  
  return allContourArray;
}

//_____________________________________________________________________________
void 
AliMUONContourMakerTest::GenerateTransformations(AliMpExMap*& real, AliMpExMap*& exploded)
{
  /// Generate geometric transformations to be used to compute the contours
  /// (real are, as the name implies, real ones, while the other ones are 
  /// a bit tweaked to look fine on screen).
  
  AliCodeTimerAuto("");
  
  AliMUONGeometryTransformer transformer;
  Bool_t ok = transformer.LoadGeometryData("transform.dat");
  //  transformer.LoadGeometryData("geometry.root"); //FIXME: add a protection if geometry.root file does not exist
  if (!ok)
  {
    cout << "ERROR : cannot get geometry !" << endl;
    return;
  }
  real = new AliMpExMap;
  exploded = new AliMpExMap;
  AliMpDEIterator deIt;
  deIt.First();
  while ( !deIt.IsDone() )
  {
    Int_t detElemId = deIt.CurrentDEId();
    const AliMUONGeometryDetElement* de = transformer.GetDetElement(detElemId);
    
    real->Add(detElemId,de->GetGlobalTransformation()->Clone());
    
    TGeoHMatrix* matrix = static_cast<TGeoHMatrix*>(de->GetGlobalTransformation()->Clone());
    Double_t* translation = matrix->GetTranslation();
        
    if ( AliMpDEManager::GetStationType(detElemId) == AliMp::kStation345 ) 
    {
      translation[0] *= 1.0;
      translation[1] *= 1.5; 
    }
    else
    {
      Double_t shift = 5; // cm
      Double_t xshift[] = { shift, -shift, -shift, shift };
      Double_t yshift[] = { shift, shift, -shift, -shift };
      Int_t ishift = detElemId % 100;
      
      translation[0] += xshift[ishift];
      translation[1] += yshift[ishift];
    }
    matrix->SetTranslation(translation);
    exploded->Add(detElemId,matrix);
    deIt.Next();
  }
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

