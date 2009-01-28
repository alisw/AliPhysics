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

#include "AliMUONPainterContourMaker.h"

#include "AliMUONPainterContour.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDigit.h"
#include "AliMpConnection.h"
#include "AliMpConstants.h"
#include "AliMpDEManager.h"
#include "AliMpExMap.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifType.h"
#include "AliMpSector.h"
#include "AliMpSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpStationType.h"
#include "AliMpVMotif.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include <Riostream.h>
#include <TArrayI.h>
#include <TGeoMatrix.h>
#include <TLine.h>
#include <TMap.h>
#include <TMath.h>
#include <TMathBase.h>
#include <TObjArray.h>
#include <TPolyLine.h>
#include <cassert>
#include <float.h>

/// \class AliMUONPainterContourMaker
///
/// A class to build painter contours. 
///
/// The basics are to build one manu contour, and then to merge contours
/// to build higher order objects, like PCBS, DEs, etc...
///
/// \author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONPainterContourMaker)
ClassImp(AliMUONPainterContourMaker::AliMUONNeighbour)
///\endcond

//_____________________________________________________________________________
Int_t
AliMUONPainterContourMaker::AliMUONNeighbour::Compare(const TObject* obj) const
{
  /// Compare two neighbours objects
  
  const AliMUONNeighbour* n = static_cast<const AliMUONNeighbour*>(obj);
  
  if ( Position().X() < n->Position().X() )
  {
    return -1;
  }
  else if ( Position().X() > n->Position().X() )
  {
    return 1;
  }
  else
  {
    // same X
    if ( Position().Y() < n->Position().Y() )
    {
      return -1;
    }
    else if ( Position().Y() > n->Position().Y() )
    {
      return 1;
    }
  }
  return 0;
}

//_____________________________________________________________________________
void
AliMUONPainterContourMaker::AliMUONNeighbour::Print(Option_t*) const
{
  /// Printout
  cout << Form("ID %10d DE %4d Manu %4d Channel %2d "
               "(X,Y)=(%7.3f,%7.3f) L,R,T,B=%1d,%1d,%1d,%1d",
               ID(),
               AliMUONVDigit::DetElemId(ID()),
               AliMUONVDigit::ManuId(ID()),
               AliMUONVDigit::ManuChannel(ID()),
               Position().X(),Position().Y(),
               HasLeftNeighbour(),HasRightNeighbour(),
               HasTopNeighbour(),HasBottomNeighbour())
  << endl;
}

//_____________________________________________________________________________
AliMUONPainterContourMaker::AliMUONPainterContourMaker(AliMpExMap* globalTransformations)
: TObject(), 
  fGlobalTransformations(globalTransformations),
  fLocalManuContours(new TMap),
  fContours(new TMap)
{
    /// ctor
    fContours->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUONPainterContourMaker::~AliMUONPainterContourMaker()
{
  /// dtor
  fLocalManuContours->DeleteAll();
  delete fLocalManuContours;
  fContours->DeleteAll();
  delete fContours;
}

//_____________________________________________________________________________
void
AliMUONPainterContourMaker::Add(AliMUONPainterContour* contour)
{
  /// Add a contour to our store of contours
  fContours->Add(new TObjString(contour->GetName()),contour);
}

//_____________________________________________________________________________
void 
AliMUONPainterContourMaker::AddSegment(TObjArray& segments, Double_t x1, Double_t y1,
                                       Double_t x2, Double_t y2, Int_t id) const
{
  /// Add one segment defined by (x1,y1,x2,y2) to the array of segments
  AliCodeTimerAuto("")
  AliDebug(1,Form("AddSegment %7.3f,%7.3f -> %7.3f,%7.3f",x1,y1,x2,y2));
  TLine* line = new TLine(x1,y1,x2,y2);
  line->SetUniqueID(id);
  segments.Add(line);
}

//_____________________________________________________________________________
Bool_t
AliMUONPainterContourMaker::HasLine(const TObjArray& segments,
                                    const TLine& line) const
{
  /// Check whether line is already part of segments array
  
  TIter next(&segments);
  TLine* l;
  
  while ( ( l = static_cast<TLine*>(next()) ) )
  {
    if ( IsEqual(line,*l) ) return kTRUE;
  }
  
  return kFALSE;
}

//_____________________________________________________________________________
void 
AliMUONPainterContourMaker::AddSegments(TObjArray& segments, 
                                        const AliMUONPainterContour& contour) const

{
  /// Add all the segments (that are not already there) 
  /// of contour to the segments array
  
  AliCodeTimerAuto("")
  
  const TObjArray* pl = contour.AsPolyLines();

  TIter next(pl);

  Int_t n(0);
  
  TPolyLine* line;
  
  while ( ( line = static_cast<TPolyLine*>(next()) ) )
  {
    n += line->GetLastPoint();
  }  
  
  AliDebug(1,Form("Adding %d groups (%d lines) from contour %s ",pl->GetLast()+1,n,contour.GetName()));

  next.Reset();
  
  while ( ( line = static_cast<TPolyLine*>(next()) ) )    
  {
    AliDebug(1,"line=");
//    StdoutToAliDebug(1,line->Print(););
    for ( Int_t i = 0; i < line->GetLastPoint(); ++i ) 
    {
      Double_t x1 = line->GetX()[i];
      Double_t y1 = line->GetY()[i];
      Double_t x2 = line->GetX()[i+1];
      Double_t y2 = line->GetY()[i+1];
      
      TLine* l = new TLine(x1,y1,x2,y2);
      
      if ( !HasLine(segments,*l) )
      {
        segments.Add(l);      
        AliDebug(1,Form("Adding line %s",LineAsString(*l).Data()));
      }
      else
      {
        AliDebug(1,Form("Line %s is already there",LineAsString(*l).Data()));
      }      
    }
  }
}

//_____________________________________________________________________________
AliMUONPainterContour*
AliMUONPainterContourMaker::ConvertEdgePadsToContour(TObjArray& ePads, 
                                                     const char* name) const
{
  /// Convert an array of edge pads into a contour of a given name
  
  AliCodeTimerAuto("")
  
  ePads.Sort();
  
  AliDebug(1,Form("%d pads to convert:",ePads.GetEntries()));
//  StdoutToAliDebug(1,ePads.Print();)
    
  TObjArray segments;
  segments.SetOwner(kTRUE);
  
  TIter nextPad(&ePads);
  AliMUONNeighbour* ne;
  
  while ( ( ne = static_cast<AliMUONNeighbour*>(nextPad()) ) )
  {
    Int_t id = ne->ID();
    
    if ( ! ne->HasLeftNeighbour() )
    {
      AddSegment(segments,ne->LowerLeft().X(),ne->LowerLeft().Y(),
                 ne->LowerLeft().X(),ne->UpperRight().Y(),id);
    }
    if ( ! ne->HasRightNeighbour() )
    {
      AddSegment(segments,ne->UpperRight().X(),ne->LowerLeft().Y(),
                 ne->UpperRight().X(),ne->UpperRight().Y(),id);
    }
    if ( ! ne->HasTopNeighbour() )
    {
      AddSegment(segments,ne->LowerLeft().X(),ne->UpperRight().Y(),
                 ne->UpperRight().X(),ne->UpperRight().Y(),id);
    }
    if ( ! ne->HasBottomNeighbour() )
    {
      AddSegment(segments,ne->LowerLeft().X(),ne->LowerLeft().Y(),
                 ne->UpperRight().X(),ne->LowerLeft().Y(),id);
    }    
  }
  
  return ConvertSegmentsToContour(segments,name);
}

//_____________________________________________________________________________
void 
AliMUONPainterContourMaker::PrintLine(const TLine& line, const char* msg) const
{
  /// Printout of a line
  cout << Form("%10s %s",
               msg,LineAsString(line).Data()) << endl;   
}

//_____________________________________________________________________________
TString 
AliMUONPainterContourMaker::LineAsString(const TLine& line, Bool_t slope) const
{
  /// Return a string representation of the line
  
  TString rv(Form("%7.3f,%7.3f -> %7.3f,%7.3f",
                  line.GetX1(),line.GetY1(),
                  line.GetX2(),line.GetY2()));
             
  if ( slope ) 
  {
    if ( IsHorizontal(line) ) rv += " H";
    else if ( IsVertical(line) ) rv += " V";
    else rv += Form(" (slope %e)",Slope(line));
  }
  
  return rv;
}

//_____________________________________________________________________________
void 
AliMUONPainterContourMaker::PrintSegments(const TObjArray& segments) const
{
  /// Printout of segment arrays (debug)
  
  for ( Int_t i = 0; i <= segments.GetLast(); ++i )
  {
    TLine* l = static_cast<TLine*>(segments.UncheckedAt(i));
    
    cout << Form("***--- i %4d",i);
    if ( l ) 
    {
      PrintLine(*l);
    }
    else
    {
      cout << " line is null ?" << endl;
    }
  }
}

//_____________________________________________________________________________
TLine* 
AliMUONPainterContourMaker::AddToLine(TPolyLine& line, TObjArray& segments, Int_t i) const
{
  /// Add one segment (taken from position i in array) into polyline
  
  AliDebug(1,Form("i=%d",i));
  TLine* l = static_cast<TLine*>(segments.UncheckedAt(i));
  if (l)
  {
    line.SetNextPoint(l->GetX1(),l->GetY1());
    line.SetNextPoint(l->GetX2(),l->GetY2());
  }
  else
  {
    AliError(Form("Did not find the line at i=%d",i));
    PrintSegments(segments);
  }
  return l;
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterContourMaker::FindPoint(Double_t x, Double_t y, 
                                      TObjArray& segments) const
{
  /// Find if point (x,y) is in segments array, and return
  /// its index (=position within array)
  
  TIter next(&segments);
  TLine* l;
  
  while ( ( l = static_cast<TLine*>(next()) ) )
  {
    if ( IsEqual(l->GetX1(),x) && IsEqual(l->GetY1(),y) )
    {
      return segments.IndexOf(l);
    }
  }
  AliError(Form("Did not find point %7.3f %7.3f in those segments:",x,y));
//  StdoutToAliDebug(1,PrintSegments(segments););
  return -1;
}

//_____________________________________________________________________________
AliMUONPainterContour*
AliMUONPainterContourMaker::ConvertSegmentsToContour(TObjArray& segments, 
                                                     const char* name) const
{
  /// Convert an array of segments into a contour
  
  AliDebug(1,"");
  AliCodeTimerAuto("");
  
  AliMUONPainterContour* contour = new AliMUONPainterContour(name);

  Int_t n(0); // this is a protection against infinite loop (used for debug only)
  
  while ( segments.GetLast() >= 0 && n < 100 ) 
  {
    TPolyLine lines;
    TIter next(&segments);
    TLine* l;
    
    while ( ( l = static_cast<TLine*>(next() ) ) )
    {
      TLine* inserted = InsertSegment(lines,*l);
      if ( inserted ) 
      {
        segments.Remove(inserted);
        next.Reset();
      }
      
      // check for closure
      if ( IsLineClosed(lines) ) 
      {
        AliDebug(1,"Line closed. Starting a new one");
        break;
      }
    }
    
    TPolyLine* sl = Simplify(lines);
    
    contour->AdoptPolyLine(sl);
    ++n;
  }
  
  if ( segments.GetLast() >= 0 ) 
  {
    AliError("segment should be empty by now");
//    StdoutToAliError(PrintSegments(segments););
  }
  
  return contour;
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterContourMaker::FindPoint(const TPolyLine& lines, Double_t x, Double_t y) const
{
  /// Return position of (x,y) within the polyline
  
  AliCodeTimerAuto("")
  
  for ( Int_t i = 0; i < lines.Size(); ++i ) 
  {
    if ( IsEqual(lines.GetX()[i],x) && IsEqual(lines.GetY()[i],y) ) 
    {
      return i;
    }
  }
  return -1;
}

//_____________________________________________________________________________
void
AliMUONPainterContourMaker::CleanSegments(TObjArray& segments,
                                          const TArrayI& toBeRemoved) const
{
  /// Remove segments at indices stored in toBeRemoved array
  for ( Int_t i = 0; i < toBeRemoved.GetSize(); ++i ) 
  {
    if ( toBeRemoved[i] )
    {
      segments.RemoveAt(i);
    }
  }
  segments.Compress();
}

//_____________________________________________________________________________
Int_t
AliMUONPainterContourMaker::SplitSegments(TObjArray& segments) const
{
  /// Split segments that have partial overlap 
  
  AliCodeTimerAuto("")
  
  TArrayI toBeRemoved(segments.GetLast()+1);
  toBeRemoved.Reset(0);
  Bool_t added(kFALSE);
  
  for ( Int_t i = 0; i <= segments.GetLast() && !added; ++i ) 
  {
    if ( toBeRemoved[i] ) continue;
    
    TLine* li = static_cast<TLine*>(segments.UncheckedAt(i));
    
    for ( Int_t j = i+1; j <= segments.GetLast() && !added; ++j ) 
    {
      if ( toBeRemoved[j] ) continue;
      
      TLine* lj = static_cast<TLine*>(segments.UncheckedAt(j));
      
      Int_t o = Overlap(*li,*lj);

      if ( o ) 
      {
        toBeRemoved[i] = toBeRemoved[j] = 1;
        
        Double_t x[] = { li->GetX1(), lj->GetX1(), li->GetX2(), lj->GetX2() };
        Double_t y[] = { li->GetY1(), lj->GetY1(), li->GetY2(), lj->GetY2() };
        
        Double_t xmin(FLT_MAX), ymin(FLT_MAX);
        Double_t xmax(-FLT_MAX), ymax(-FLT_MAX);

        for ( Int_t k = 0; k < 4; ++k )
        {
          xmin = TMath::Min(x[k],xmin);
          ymin = TMath::Min(y[k],ymin);
          xmax = TMath::Max(x[k],xmax);
          ymax = TMath::Max(y[k],ymax);
        }
        
        TLine fullLine(xmin,ymin,xmax,ymax);
        
        for ( Int_t i1 = 0; i1 < 4; ++i1 )
        {
          for ( Int_t j1 = i1+1; j1 < 4; ++j1 )
          {
            if ( TMath::Abs(i1-j1) != 2 ) 
            {
              TLine test(x[i1],y[i1],x[j1],y[j1]);

              Bool_t isFullLine = IsEqual(test,fullLine);
              
              if ( !IsPoint(test) && !isFullLine ) 
              {
                segments.Add(new TLine(test));
                added = kTRUE;
              }              
            }
          }
        }
      }
    }
  }
  
  CleanSegments(segments,toBeRemoved);
  
  return added;
}

//_____________________________________________________________________________
Bool_t
AliMUONPainterContourMaker::ShouldBeRemoved(const TObjArray& contours,
                                            Double_t x, Double_t y) const
{
  /// Tells whether or not a point can be removed, because it lies
  /// inside the global contour
  
  const Double_t kPrecision(AliMpConstants::LengthTolerance());
  const Double_t kShiftX[] = { kPrecision,kPrecision,-kPrecision,-kPrecision };
  const Double_t kShiftY[] = { kPrecision,-kPrecision,kPrecision,-kPrecision };
  
  TIter next(&contours);
  AliMUONPainterContour* contour;

  Int_t n(0);

  while ( ( contour = static_cast<AliMUONPainterContour*>(next()) ) )
  {
    for ( Int_t i = 0; i < 4; ++i )
    {
      if ( contour->IsInside( x + kShiftX[i], y + kShiftY[i]) )
      {
        ++n;
      }
    }
  }
  
  return (n>=4);
}

//_____________________________________________________________________________
Int_t
AliMUONPainterContourMaker::RemoveInsideSegments(const TObjArray& contours,
                                                 TObjArray& segments) const
{
  /// Remove segments that have 2 triple points

  AliCodeTimerAuto("")
  
  TArrayI toBeRemoved(segments.GetLast()+1);
  toBeRemoved.Reset(0);
  
  for ( Int_t i = 0; i <= segments.GetLast(); ++i ) 
  {
    TLine* line = static_cast<TLine*>(segments.UncheckedAt(i)); 
    
    Double_t x = (line->GetX1() + line->GetX2())/2.0;
    Double_t y = (line->GetY1() + line->GetY2())/2.0;
    
    if ( ShouldBeRemoved(contours,x,y) )
    {
      toBeRemoved[i] = 1;
    }
  }
    
  Int_t before = segments.GetLast()+1;
  
  CleanSegments(segments,toBeRemoved);
  
  Int_t after = segments.GetLast()+1;
  
  AliDebug(1,Form("# of segments before = %d after = %d",before,after));  
  
  return after-before;
}

//_____________________________________________________________________________
AliMUONPainterContour* 
AliMUONPainterContourMaker::MergeContours(const TObjArray& contours,
                                          const char* contourName) const
{
  /// Merge an array of contours into a single contour, with a given name
  
  AliCodeTimerAuto("");
  
  AliDebug(1,Form("Merging %d contours into %s",contours.GetLast()+1,contourName));
  
  if ( contours.GetSize() == 0 ) return 0x0;
  
  TIter next(&contours);
  AliMUONPainterContour* contour;
  
  TObjArray segments;
  segments.SetOwner(kTRUE);
  
  while ( ( contour = static_cast<AliMUONPainterContour*>(next()) ) )
  {
    AddSegments(segments,*contour);
  }
  
//  AliDebug(1,"After AddSegments");
//  StdoutToAliDebug(1,PrintSegments(segments));
  
  while (SplitSegments(segments)) {}

//  AliDebug(1,"After SplitSegments");
//  StdoutToAliDebug(1,PrintSegments(segments));

//  if (!SanityCheck(contours,segments))
//  {
//    return 0x0;
//  }

  RemoveInsideSegments(contours,segments);
  
//  if (!SanityCheck(contours,segments))
//  {
//    return 0x0;
//  }
    
//  AliDebug(1,"After RemoveInsideSegments");
//  StdoutToAliDebug(1,PrintSegments(segments););
    
//  if (!SanityCheck(contours,segments))
//  {
//    return 0x0;
//  }
  
  return ConvertSegmentsToContour(segments,contourName);
}

//_____________________________________________________________________________
TString
AliMUONPainterContourMaker::NameIt(const AliMpMotifPosition& motifPosition) const
{
  /// Get the name of an AliMpMotifPosition
  
  AliMpVMotif* motif = motifPosition.GetMotif();
  TString name(Form("%s",motif->GetID().Data()));
  
  for ( Int_t i = 0; i < motif->GetNofPadDimensions(); ++i )
  {
    TVector2 padDim = motif->GetPadDimensions(i);
    name += Form("/%7.3f-%7.3f:",padDim.X(),padDim.Y());
  }
  return name;
}

//_____________________________________________________________________________
AliMUONPainterContour*
AliMUONPainterContourMaker::FindLocalManuContour(Int_t detElemId, Int_t manuId) const
{
  /// Get a pre-computed manu contour (in local coordinates)
  AliCodeTimerAuto("")
  
  AliMpMotifPosition* motifPos = FindMotifPosition(detElemId,manuId);
  
  TObject* o = fLocalManuContours->GetValue(NameIt(*motifPos));
  
  if (o) return static_cast<AliMUONPainterContour*>(o);
  return 0x0;
}

//_____________________________________________________________________________
AliMpMotifPosition*
AliMUONPainterContourMaker::FindMotifPosition(Int_t detElemId, Int_t manuId) const
{
  /// Find a given motifPosition object
  
  AliCodeTimerAuto("")
  
  AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
  
  if ( stationType == AliMp::kStation345 ) 
  {
    const AliMpSlat* kSlat 
      = AliMpSegmentation::Instance()->GetSlatByElectronics(detElemId,manuId);
    if ( ! kSlat ) {
      AliFatal(Form("Could not find motif for DE %d manu %d",detElemId,manuId));
    }
    return kSlat->FindMotifPosition(manuId);
  }
  else
  {
    const AliMpSector* kSector 
      = AliMpSegmentation::Instance()->GetSectorByElectronics(detElemId,manuId);
    if ( ! kSector ) {
      AliFatal(Form("Could not find motif for DE %d manu %d",detElemId,manuId));
    }
    return kSector->GetMotifMap()->FindMotifPosition(manuId);
  }
  return 0x0;
}

//_____________________________________________________________________________
AliMUONPainterContour* 
AliMUONPainterContourMaker::GenerateManuContour(const char* name, 
                                                Int_t detElemId, Int_t manuId,
                                                AliMUONAttPainter viewType) const
{  
  /// Generate the contour for a given manu
  
  AliDebug(3,Form("DE %04d ManuID %04d Name %s",detElemId,manuId,name));
  
  AliCodeTimerAuto("")
  
  AliMpMotifPosition* motifPosition = FindMotifPosition(detElemId,manuId);
  AliMpVMotif* motif = motifPosition->GetMotif();
  
  AliMUONPainterContour* contour = FindLocalManuContour(detElemId,manuId); 
  // do we already have the local contour for that manu ?
  
  // no : build it
  if  (!contour)
  {
    AliCodeTimerAuto("Generation of local contour");
    TObjArray ePads;
    ePads.SetOwner(kTRUE);
    AliMpMotifType* motifType = motif->GetMotifType();
    AliDebug(3,Form("motifType %s",motifType->GetID().Data()));
    
//    for ( Int_t i = 0; i <= motifType->GetNofPads(); ++i ) 
    for ( Int_t i = 0; i <= AliMpConstants::ManuNofChannels(); ++i ) 
    {
//      AliMpConnection* connection = motifType->FindConnectionByPadNum(i);
      AliMpConnection* connection = motifType->FindConnectionByGassiNum(i);

      AliDebug(3,Form("connection i =%d",i));
      
      if ( connection ) 
      {
        AliMpIntPair indices = connection->LocalIndices();
        Bool_t left(kTRUE);
        Bool_t right(kTRUE);
        Bool_t top(kTRUE);
        Bool_t bottom(kTRUE);
        
        if ( !motifType->FindConnectionByLocalIndices(indices+AliMpIntPair(1,0)) )
        {
          right = kFALSE;
        }
        if ( !motifType->FindConnectionByLocalIndices(indices+AliMpIntPair(-1,0)) )
        {
          left = kFALSE;
        }
        if ( !motifType->FindConnectionByLocalIndices(indices+AliMpIntPair(0,1)) )
        {
          top = kFALSE;
        }
        if ( !motifType->FindConnectionByLocalIndices(indices+AliMpIntPair(0,-1)) )
        {
          bottom = kFALSE;
        }
        
        AliDebug(3,Form("indices=(%3d,%3d) L %d R %d T %d B %d",
                        indices.GetFirst(),indices.GetSecond(),
                        left,right,top,bottom));
        
        TVector2 position = motif->PadPositionLocal(indices);
        TVector2 dimensions = motif->GetPadDimensions(indices);

        if ( !left  || !right || !top  || !bottom )
        {
          // the pad is on the edge
          Int_t id = AliMUONVDigit::BuildUniqueID(detElemId,manuId,
                                                  connection->GetGassiNum(),0);
          ePads.AddLast(new AliMUONNeighbour(id,position,dimensions,left,right,top,bottom));
        }
      }
    }
    
    contour = ConvertEdgePadsToContour(ePads,NameIt(*motifPosition));
    
    AliDebug(1,Form("localContour:"));
//    StdoutToAliDebug(1,contour->Print("full"));
    // register the local contour
    fLocalManuContours->Add(new TObjString(contour->GetName()),contour);
  }
  
  AliMUONPainterContour* globalContour = static_cast<AliMUONPainterContour*>(contour->Clone(name));
  
  // once we have the local contour, convert it to global
  
  TVector2 pos(motifPosition->Position());

  if ( AliMpDEManager::GetStationType(detElemId) == AliMp::kStation345 ) 
  {
    const AliMpSlat* slat = AliMUONPainterHelper::Instance()->GetSlat(detElemId,manuId);
    pos -= slat->Position();
  }
  globalContour->Offset(pos);
  TGeoHMatrix* matrix = static_cast<TGeoHMatrix*>(fGlobalTransformations->GetValue(detElemId));
  globalContour->Transform(*matrix);

  if ( viewType.IsBackView() )
  {
    AliWarning("Got a back view : will rotate ! This has not been really tested. Please do so now !");
    TGeoRotation rot;
    rot.RotateZ(180);
    globalContour->Transform(rot);
  }
  
  return globalContour;
}

//_____________________________________________________________________________
AliMUONPainterContour*
AliMUONPainterContourMaker::GetContour(const char* name) const
{
  /// Get contour by name
  
  TObject* o = fContours->GetValue(name);
  return static_cast<AliMUONPainterContour*>(o);
}

//_____________________________________________________________________________
Bool_t
AliMUONPainterContourMaker::HasContour(const char* name) const
{
  /// Whether contour named "name" exists
  TObject* o = fContours->GetValue(name);
  if (o) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________
TLine* 
AliMUONPainterContourMaker::InsertSegment(TPolyLine& lines, TLine& l) const
{
  /// Insert line into polyline, at the correct position
  
  AliCodeTimerAuto("")
//  AliDebug(2,Form("Trying to insert %7.3f,%7.3f -> %7.3f,%7.3f from "
//                  "(DE,manu,ch)=(%d,%d,%d) into",
//                  l.GetX1(),l.GetY1(),l.GetX2(),l.GetY2(),
//                  AliMUONVDigit::DetElemId(l.GetUniqueID()),
//                  AliMUONVDigit::ManuId(l.GetUniqueID()),
//                  AliMUONVDigit::ManuChannel(l.GetUniqueID())));
  
  if ( lines.Size()==0 ) 
  {
//    AliDebug(2,"Starting line");
//    
    lines.SetNextPoint(l.GetX1(),l.GetY1());
    lines.SetNextPoint(l.GetX2(),l.GetY2());
    return &l;
  }
  
  Int_t i1 = FindPoint(lines,l.GetX1(),l.GetY1());
  Int_t i2 = FindPoint(lines,l.GetX2(),l.GetY2());
  
  if ( i1 < 0 && i2 < 0 ) 
  {
//    AliDebug(2,"Not yet");
    return 0x0;
  }
  
  if ( i1 >= 0 && i2 >= 0 )
  {
    if ( i1==0 )
    {
      lines.SetNextPoint(l.GetX1(),l.GetY1());
    }
    else if ( i2==0 )
    {
      lines.SetNextPoint(l.GetX2(),l.GetY2());
    }
    else
    {
      AliError("Segment already there but does not correspond to ending the polyline !");
      AliError(Form("Segment is %7.3f,%7.3f -> %7.3f,%7.3f and existing points are : ",
                    l.GetX1(),l.GetY1(),l.GetX2(),l.GetY2()));
                      
      for ( Int_t i = 0; i < lines.Size(); ++i ) 
      {
        AliError(Form("Point %2d X %7.3f Y %7.3f",i,lines.GetX()[i],lines.GetY()[i]));
      }
//      TObject* o(0x0);
//      o->Print(); // to crash and throw gdb...
    }
    return &l;
  }
  
  Double_t x = (i1>=0) ? l.GetX2() : l.GetX1();
  Double_t y = (i1>=0) ? l.GetY2() : l.GetY1();
  
  Int_t iref = ( i1 >= 0 ? i1 : i2 ) ;
  
  Bool_t firstPoint = ( iref == 0 );
  
  if ( firstPoint ) 
  {
    // must insert segment before
    lines.SetPolyLine(lines.Size()+1);
//    AliDebug(2,Form("Inserting %7.3f,%7.3f",x,y));
    for ( Int_t i = lines.Size()-1; i > 0; --i ) 
    {
      lines.SetPoint(i,lines.GetX()[i-1],lines.GetY()[i-1]);
    }
    lines.SetPoint(0,x,y);
  }
  else
  {
//    AliDebug(2,Form("Appending %7.3f,%7.3f",x,y));
    lines.SetNextPoint(x,y);
  }
  
  return &l;
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterContourMaker::IsEqual(Double_t x, Double_t y) const
{
  /// Whether x==y
  
  if ( TMath::Abs(x-y) < AliMpConstants::LengthTolerance() ) return kTRUE;
  else return kFALSE;      
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterContourMaker::IsEqual(const TLine& line1,
                                    const TLine& line2) const
{
  /// Whether line1 == line2
  
  Bool_t check1 =  
    IsEqual(line1.GetX1(),line2.GetX1()) && 
    IsEqual(line1.GetY1(),line2.GetY1()) &&
    IsEqual(line1.GetX2(),line2.GetX2()) &&
    IsEqual(line1.GetY2(),line2.GetY2());
  
  Bool_t check2 =  
    IsEqual(line1.GetX1(),line2.GetX2()) && 
    IsEqual(line1.GetY1(),line2.GetY2()) &&
    IsEqual(line1.GetX2(),line2.GetX1()) &&
    IsEqual(line1.GetY2(),line2.GetY1());
  
  return (check1 || check2);
}

//_____________________________________________________________________________
Double_t 
AliMUONPainterContourMaker::Slope(const TLine& line) const
{
  /// Get the slope of line
  
  Double_t x = TMath::Abs(line.GetX2() - line.GetX1());
  
  if ( x  < AliMpConstants::LengthTolerance() ) return FLT_MAX;
  
  return TMath::Abs(line.GetY2() - line.GetY1())/x;
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterContourMaker::IsPoint(const TLine& line) const
{
  /// Whether the line is a point (sic ;-) )
  return 
  IsEqual(line.GetX1(),line.GetX2()) && 
  IsEqual(line.GetY1(),line.GetY2());
}

//_____________________________________________________________________________
TLine
AliMUONPainterContourMaker::Shift(const TLine& line, Double_t x, Double_t y) const
{
  /// Shift the line by a given offset
  
  return TLine(line.GetX1()-x,line.GetY1()-y,line.GetX2()-x,line.GetY2()-y);
}

//_____________________________________________________________________________
Bool_t
AliMUONPainterContourMaker::SameDirection(const TLine& line1, const TLine& line2) const
{
  /// Whether both lines have the same direction.
  
  TLine l1 = Shift(line1,line1.GetX1(),line1.GetY1());
  TLine l2 = Shift(line2,line2.GetX1(),line2.GetY1());
  
  Double_t v = l1.GetX2()*l2.GetX2() + l1.GetY2()*l2.GetY2();
  
  return v > 0 ;
}

//_____________________________________________________________________________
void 
AliMUONPainterContourMaker::Swap(TLine& line) const
{
  /// Swap both points of the line
  
  Double_t x = line.GetX1();
  Double_t y = line.GetY1();
  
  line.SetX1(line.GetX2());
  line.SetY1(line.GetY2());
  line.SetX2(x);
  line.SetY2(y);
}

//_____________________________________________________________________________
Int_t
AliMUONPainterContourMaker::IsInRange(Double_t x, Double_t a, Double_t b,
                                      Bool_t strict) const
{
  /// Whether w is in [a,b] (if strict=kFALSE) or in ]a,b[ (if strict=kTRUE)
  
  if ( a > b ) 
  {
    Double_t tmp(b);
    b = a;
    a = tmp;
  }
  
  Bool_t rv(kFALSE);
  
  if ( strict )  
  {
    rv = ( x > a && x < b );
  }
  else
  {
    rv = ( x >= a && x <= b);
  }

  AliDebug(4,Form("x = %7.3f a = %7.3f b = %7.3f strict = %d IsInRange = %d",x,a,b,strict,rv));
  
  return rv;
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterContourMaker::IsInLine(const TLine& line,
                                     Double_t x,
                                     Double_t y,
                                     Bool_t strict) const
{
  /// Check whether point (x,y) is belonging to the line segment
  /// by computing the distance point to line
  /// line1 must not be a single point.
  /// Returns the number of *coordinates* that matches, for a point
  /// that lies on line (if point is not on the line, returns 0 always).
  /// For instance, if (x,y) is on the line (and strict=kFALSE), 
  /// it will return 1 if x *or* y corresponds to line.GetX1() or X2 or Y1 or Y2,
  /// and 2 if the pair (x,y) corresponds to one of the line points.
  
  Double_t x1 = line.GetX1();
  Double_t x2 = line.GetX2();
  Double_t y1 = line.GetY1();
  Double_t y2 = line.GetY2();
  
  Double_t distance = TMath::Abs( (x2-x1)*(y1-y) - (x1-x)*(y2-y1) );
  
  distance /= TMath::Sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
  
  Bool_t online = ( distance < AliMpConstants::LengthTolerance() ) ;
  
  Int_t rv(0);
  
  if (online)
  {
    // point is on the line, 
    // check in addition that it's within the segment
    
    rv = IsInRange(x,x1,x2,strict) + IsInRange(y,y1,y2,strict);
  }
  else
  {
    rv = 0;
  }
  
  AliDebug(4,Form("Point (%7.3f,%7.3f) isinline=%d in line %s",
                  x,y,rv,LineAsString(line).Data()));
  
  return rv;
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterContourMaker::IsInside(const TLine& line1,
                                     const TLine& line2,
                                     Bool_t useEndPoints) const
{
  /// Check whether one or both points of line2 are within line1.
  /// Both line1 and line2 must have the same slope 
  /// and the same direction

  if (!IsEqual(Slope(line1),Slope(line2))) return 0;
  
  TLine l2(line2);
  
  if (!SameDirection(line1,line2)) 
  {
    Swap(l2);
  }
      
  Int_t rv = 
    IsInLine(line1,l2.GetX1(),l2.GetY1(),!useEndPoints) +
    IsInLine(line1,l2.GetX2(),l2.GetY2(),!useEndPoints);
  
  assert(rv<=4);
  
  return rv;
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterContourMaker::IsInside(const TObjArray& segments, 
                                     const TLine& line) const
{
  /// Whether the segment (line) is contained inside the contour defined
  /// by all the segments (i.e. is it on the boundary or not)
  /// Basic (and dirty) implementation only working with horizontal and vertical lines.
  /// I know there must be a better way to do it, but it took me way too long
  /// to get this stuff working, so I'm giving up on the optimisation/cleaning,
  /// at least for now...
  /// If you'd like to clean this (while keeping it working in all cases), be
  /// my guest and do it ;-) )
  
  Int_t p1 = CountPoint(segments,line.GetX1(),line.GetY1());
  Int_t p2 = CountPoint(segments,line.GetX2(),line.GetY2());
  
  Bool_t triplet = ( p1 >= 3 || p2 >= 3 );
  
  AliDebug(4,Form("IsInside(segments,%s) triplet=%d",
                  LineAsString(line).Data(),triplet));
  
  if (!triplet) return kFALSE;
    
  Bool_t top(kFALSE), bottom(kFALSE), left(kFALSE), right(kFALSE);
  
  Bool_t vertical(IsVertical(line));
  Bool_t horizontal(IsHorizontal(line));
  
  if (!vertical && !horizontal ) 
  {
    AliFatal("Only working with horizontal and vertical lines");
  }
  
  for ( Int_t i = 0; i <= segments.GetLast(); ++i ) 
  {
    TLine* l = static_cast<TLine*>(segments.UncheckedAt(i));
    
    if ( IsEqual(*l,line) ) continue;
    
    if ( vertical && IsVertical(*l) )
    {
      TLine tmpLine(l->GetX1(),line.GetY1(),
                    l->GetX1(),line.GetY2());
      
      AliDebug(4,Form("i=%2d VV\nIsInside(l=%s,%s)=%d\nIsInside(%s,l=%s)=%d",
                      i,
                      LineAsString(*l).Data(),LineAsString(tmpLine).Data(),
                      IsInside(*l,tmpLine,kTRUE),
                      LineAsString(tmpLine).Data(),LineAsString(*l).Data(),
                      IsInside(tmpLine,*l,kTRUE)));
                      
      if ( IsInside(*l,tmpLine,kTRUE) == 4 || IsInside(tmpLine,*l,kTRUE) == 4 ) 
      {
        if ( l->GetX1() > line.GetX1() ) 
        {
          right = kTRUE;
        }
        else
        {
          left = kTRUE;
        }
      }
    }
    
    if ( vertical && IsHorizontal(*l) )
    {
      if ( !IsEqual(l->GetY1(),line.GetX1()) && 
           !IsEqual(l->GetY1(),line.GetY2()) &&
           IsInLine(*l,line.GetX1(),l->GetY1(),kFALSE)==2 )
        {
        if ( line.GetY2() < l->GetY1() ) 
        {
          top = kTRUE;
        }
        else if ( line.GetY2() > l->GetY1() )
        {
          bottom = kTRUE;
        }
      }
    }
    
    if ( horizontal && IsHorizontal(*l) )
    {
      TLine tmpLine(line.GetX1(),l->GetY1(),
                    line.GetX2(),l->GetY1());
      
      AliDebug(4,Form("i=%2d HH\nIsInside(%s,%s)=%d\nIsInside(%s,%s)=%d",
                      i,
                      LineAsString(*l).Data(),LineAsString(tmpLine).Data(),
                      IsInside(*l,tmpLine),
                      LineAsString(tmpLine).Data(),LineAsString(*l).Data(),
                      IsInside(tmpLine,*l)));
      
      if ( IsInside(*l,tmpLine) == 4 || IsInside(tmpLine,*l) == 4 ) 
      {
        if ( l->GetY1() > line.GetY1() ) 
        {
          top = kTRUE;
        }
        else
        {
          bottom = kTRUE;
        }
      }
    }
    
    if ( horizontal && IsVertical(*l) )
    {
      if ( !IsEqual(l->GetX1(),line.GetX1()) && 
           !IsEqual(l->GetX1(),line.GetX2()) &&
           IsInLine(*l,l->GetX1(),line.GetY1(),kFALSE)==2 )
        {
        if ( line.GetX2() < l->GetX1() ) 
        {
          right = kTRUE;
        }
        else if ( line.GetX2() > l->GetX1() )
        {
          left = kTRUE;
        }
      }
    }
    
  }
  
  Bool_t rv(kFALSE);

  AliDebug(3,Form("%s %s R %d L %d T %d B% d IsInside %d",
                  IsVertical(line) ? 
                  "Vertical  " : 
                  "Horizontal",
                  LineAsString(line,kFALSE).Data(),right,left,top,bottom,rv));

  if ( vertical ) 
  {
    rv = (right && left) && ( top || bottom );
  }
  
  if ( horizontal ) 
  {
    rv = (top && bottom) && ( right || left );
  }
  
  return rv;
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterContourMaker::IsHorizontal(const TLine& line) const
{
  /// whether line is horizontal
  
  static Double_t l2 = AliMpConstants::LengthTolerance()*AliMpConstants::LengthTolerance();
  
  return ( Slope(line) < l2 );
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterContourMaker::IsVertical(const TLine& line) const
{
  /// whether line is vertical
  
  return ( TMath::Abs(Slope(line)) == FLT_MAX );
}

//_____________________________________________________________________________
Int_t
AliMUONPainterContourMaker::Overlap(const TLine& line1,
                                    const TLine& line2) const
{
  /// Whether line1 and line2 overlap
  
  Int_t rv(0);
  
  if ( IsEqual(line1,line2) ) 
  {
    // First things first. If both lines are the same one, 
    // they for sure overlap ;-)
    rv = 4;
  }
  else
  {
    rv = IsInside(line1,line2) + IsInside(line2,line1);
  }
  
  AliDebug(3,Form("%s and %s : overlap=%d",
                  LineAsString(line1).Data(),
                  LineAsString(line2).Data(),
                  rv));
  
  return rv;
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterContourMaker::IsLineClosed(const TPolyLine& line) const
{
  /// check if polyline is already closed (i.e. last point = first point)
  
  Double_t* x = line.GetX();
  Double_t* y = line.GetY();
  
  if ( IsEqual(x[line.GetLastPoint()],x[0]) &&
       IsEqual(y[line.GetLastPoint()],y[0]) )
  {
    return kTRUE;
  }
  else
  {
    return kFALSE;
  }
}  

//_____________________________________________________________________________
void 
AliMUONPainterContourMaker::Local2Global(Int_t detElemId, 
                                         Double_t xl, Double_t yl, Double_t zl,
                                         Double_t& xg, Double_t& yg, Double_t& zg) const
{
  /// Convert local coordinates to global ones
  TGeoHMatrix* matrix = static_cast<TGeoHMatrix*>(fGlobalTransformations->GetValue(detElemId));
  Double_t pl[3] = { xl, yl, zl };
  Double_t pg[3] = { 0., 0., 0. };
  matrix->LocalToMaster(pl, pg);
  xg = pg[0];
  yg = pg[1];
  zg = pg[2];
}

//_____________________________________________________________________________
void
AliMUONPainterContourMaker::Print(Option_t* opt) const
{
  /// Printout
  
  cout << "Local Contours" << endl;
  
  TIter next(fLocalManuContours);
  TObjString* key;
  
  while ( ( key = static_cast<TObjString*>(next()) ) )
  {
    cout << key->String().Data() << endl;
    AliMUONPainterContour* contour = static_cast<AliMUONPainterContour*>(fLocalManuContours->GetValue(key));
    contour->Print(opt);
  }

  cout << "Global Contours" << endl;
  
  TIter nextC(fContours);
  
  while ( ( key = static_cast<TObjString*>(nextC()) ) )
  {
    AliMUONPainterContour* contour = static_cast<AliMUONPainterContour*>(fContours->GetValue(key));
    contour->Print(opt);
  }
}

//_____________________________________________________________________________
Int_t
AliMUONPainterContourMaker::CountPoint(const TObjArray& segments, 
                                       Double_t x, Double_t y) const
{
  /// Count the number of times the point (x,y) appears in the segment array
  
  Int_t n(0);
  
  for ( Int_t i = 0; i <= segments.GetLast(); ++i ) 
  {
    TLine* line = static_cast<TLine*>(segments.UncheckedAt(i));
      
    if ( IsEqual(x,line->GetX1()) &&
         IsEqual(y,line->GetY1()) )
    {
      ++n;
    }
  
    if ( IsEqual(x,line->GetX2()) &&
         IsEqual(y,line->GetY2()) )
    {
      ++n;
    }
  }
  
  return n;
}

//_____________________________________________________________________________
Bool_t
AliMUONPainterContourMaker::SanityCheck(const TObjArray& contours,
                                        const TObjArray& segments, Bool_t check) const
{
  /// (debug) check 
  
  Bool_t ok(kTRUE);

  // cross-check that we have no more complete duplicates
  // and that we have no orphan point

  Double_t xmin(FLT_MAX), xmax(-FLT_MAX);
  Double_t ymin(FLT_MAX), ymax(-FLT_MAX);
  
  for ( Int_t i = 0; i <= segments.GetLast(); ++i ) 
  {
    TLine* li = static_cast<TLine*>(segments.UncheckedAt(i));
  
    if (!IsHorizontal(*li) && !IsVertical(*li))
    {
      AliError("Got an oblique line !");
      return kFALSE;
    }
    
    xmin = TMath::Min(xmin,li->GetX1());
    xmin = TMath::Min(xmin,li->GetX2());

    xmax = TMath::Max(xmax,li->GetX1());
    xmax = TMath::Max(xmax,li->GetX2());

    ymin = TMath::Min(ymin,li->GetY1());
    ymin = TMath::Min(ymin,li->GetY2());
    
    ymax = TMath::Max(ymax,li->GetY1());
    ymax = TMath::Max(ymax,li->GetY2());
    
  }
  
  AliDebug(1,Form("xmin=%7.3f ymin=%7.3f xmax=%7.3f ymax=%7.3f",
                  xmin,ymin,xmax,ymax));
  
  for ( Int_t i = 0; i <= segments.GetLast(); ++i ) 
  {
    TLine* li = static_cast<TLine*>(segments.UncheckedAt(i));

    if (!check) 
    {
      for ( Int_t j = 0; j <= segments.GetLast(); ++j ) 
      {
        TLine* lj = static_cast<TLine*>(segments.UncheckedAt(j));

        if ( i != j && IsEqual(*li,*lj) )
        {
          ok = kFALSE;
          PrintLine(*li);
          PrintLine(*lj);
          AliFatal("");
        }      
      }
    }
    

    Int_t rv(0);
    
    Double_t x = (li->GetX1()+li->GetX2())/2.0;
    Double_t y = (li->GetY1()+li->GetY2())/2.0;
    
    if ( ShouldBeRemoved(contours,x,y) ) rv = 1;
    
    AliDebug(1,Form("Line %4d %7.3f,%7.3f -> %7.3f,%7.3f [ %d ]",
                    i,
                    li->GetX1(),li->GetY1(),
                    li->GetX2(),li->GetY2(),
                    rv));
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
TPolyLine*
AliMUONPainterContourMaker::Simplify(const TPolyLine& lines) const
{
  /// try to simplify the polyline, by minimizing the number of points

  if ( lines.Size() < 3 ) 
  {
    AliError("Cannot simplify lines with less that 3 points !");
    return 0x0;
  }
  
  AliCodeTimerAuto("")
  
//  cout << "Before simplify" << endl;
//  
//  for ( Int_t i = 0; i < lines.Size(); ++i ) 
//  {
//    cout << Form("Point %3d %7.3f %7.3f",i,lines.GetX()[i],lines.GetY()[i]) << endl;
//  }
  
  TPolyLine* l = new TPolyLine;

  Double_t* x = lines.GetX();
  Double_t* y = lines.GetY();

  l->SetNextPoint(x[0],y[0]);

  Bool_t verticalCurrent = IsEqual(x[1],x[0]);
  Bool_t horizontalCurrent = IsEqual(y[1],y[0]);

  Int_t i(2);

  while ( i < lines.Size() )
  {
    Bool_t vertical = IsEqual(x[i],x[i-1]);
    Bool_t horizontal = IsEqual(y[i],y[i-1]);
    
//    cout << Form("i %3d %7.3f %7.3f vert %d horiz %d (current vert %d horiz %d)",
//                 i,x[i],y[i],vertical,horizontal,verticalCurrent,horizontalCurrent)
//      << endl;
    
    if ( ( vertical != verticalCurrent ) || 
         ( horizontal != horizontalCurrent ) )
    {
//      cout << Form("Changing direction : adding point %7.3f %7.3f",x[i-1],y[i-1]) << endl;
      l->SetNextPoint(x[i-1],y[i-1]);
      verticalCurrent = vertical;
      horizontalCurrent = horizontal;
    }
    ++i;
  }
  
  l->SetNextPoint(l->GetX()[0],l->GetY()[0]);
  
//  cout << "After simplify" << endl;
//  
//  for ( Int_t i = 0; i < l->Size(); ++i ) 
//  {
//    cout << Form("Point %3d %7.3f %7.3f",i,l->GetX()[i],l->GetY()[i]) << endl;
//  }
  
  return l;
}

//_____________________________________________________________________________
Int_t
AliMUONPainterContourMaker::Size() const
{
  /// Number of contours we have already
  
  return fContours->GetSize();
}
                                         
