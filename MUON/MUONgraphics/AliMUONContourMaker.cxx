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
/// Maker/merger of contours. Can create contour from a set polygons or
/// merger a set of contours into a single one.
/// 
/// This is based on (one of the) algorithm found in 
/// Diane L. Souvaine and Iliana Bjorling-Sachs,
/// Proceedings of the IEEE, Vol. 80, No. 9, September 1992, p. 1449
///
/// Note that besides the AliMUON prefix, nothing is really MUON specific
/// in this class...
///
/// \author Laurent Aphecetche, Subatech
///

#include "AliMUONContourMaker.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONContour.h"
#include "AliMUONPointWithRef.h"
#include "AliMUONPolygon.h"
#include "AliMUONSegment.h"
#include "AliMUONSegmentTree.h"
#include "Riostream.h"
#include "TArrayD.h"
#include "TMath.h"
#include <cassert>
#include "TArrayI.h"

/// \cond CLASSIMP
ClassImp(AliMUONContourMaker)
/// \endcond

//_____________________________________________________________________________
AliMUONContourMaker::AliMUONContourMaker() 
{
/// Default constructor
}


//_____________________________________________________________________________
AliMUONContourMaker::~AliMUONContourMaker()
{
/// Destructor
}

//_____________________________________________________________________________
AliMUONContour* 
AliMUONContourMaker::CreateContour(const TObjArray& polygons, const char* name) const
{
  /// Create the contour of the polygon array
  /// and get back the intermediate verticals and horizontal segments
  /// both arrays are arrays of AliMUONSegment objects.
  
  AliCodeTimerAuto("",0);
  
  if ( polygons.IsEmpty() ) return 0x0; // protection against user error...
  
  // Sanity check : insure that all polygons are oriented counter-clockwise
  TIter next(&polygons);
  AliMUONPolygon* pol;
  while ( ( pol = static_cast<AliMUONPolygon*>(next()) ) )
  {
    if ( !pol->IsCounterClockwiseOriented() )
    {
      AliError(Form("Got a clockwise oriented polygon in CreateContour(%s). That's not OK !",name));
      StdoutToAliError(polygons.Print());
      return 0x0;
    }
  }
  
  AliMUONContour* contour(0x0);
  
  if ( polygons.GetLast() == 0 ) 
  {
    AliCodeTimerAuto("Trivial case",1);
    contour = new AliMUONContour(name);
    pol = static_cast<AliMUONPolygon*>(polygons.First());
    contour->Add(*pol);
    contour->AssertOrientation();
    return contour;
  }
  
  TObjArray polygonVerticalEdges; // arrray of AliMUONSegment objects  
  polygonVerticalEdges.SetOwner(kTRUE);
  // get vertical edges of input polygons
  GetVerticalEdges(polygons,polygonVerticalEdges);
  
  // sort them in ascending x order
  // if same x, insure that left edges are before right edges
  // within same x, order by increasing bottommost y (see AliMUONSegment::Compare method)
  polygonVerticalEdges.Sort();
  
  if ( polygonVerticalEdges.GetLast()+1 < 2 ) 
  {
    polygons.Print();
    AliFatal(Form("Got too few edges here for createContour %s",name));
  }
  
  // Find the vertical edges of the merged contour. This is the meat of the algorithm...
  TObjArray contourVerticalEdges;
  contourVerticalEdges.SetOwner(kTRUE);
  Sweep(polygonVerticalEdges,contourVerticalEdges);
  
  TObjArray horizontals;
  horizontals.SetOwner(kTRUE);
  VerticalToHorizontal(contourVerticalEdges,horizontals);
  
  contour = FinalizeContour(contourVerticalEdges,horizontals);
  
  if ( contour && name ) contour->SetName(name);
  
  return contour;
}

//_____________________________________________________________________________
AliMUONContour* 
AliMUONContourMaker::FinalizeContour(const TObjArray& verticals,
                                     const TObjArray& horizontals) const
{  
  /// For a list of vertical and horizontal edges, we build the final
  /// contour object.
  
  AliCodeTimerAuto("",0);
  
  TObjArray all; // array of AliMUONSegment
  TObjArray inorder; // array of AliMUONSegment

  all.SetOwner(kFALSE);
  inorder.SetOwner(kFALSE);
    
  for ( Int_t i = 0; i <= verticals.GetLast(); ++i ) 
  {
    all.Add(verticals.UncheckedAt(i));
    all.Add(horizontals.UncheckedAt(i));
  }

  TArrayI alreadyAdded(all.GetLast()+1);
  alreadyAdded.Reset();
  
  Int_t i(0);
  
  AliMUONContour* contour = new AliMUONContour;
  
  int total(0);
  
  while ( !all.IsEmpty() )
  {
    total++;
    
    if ( total > 1000 ) 
    {
      AliError("Total 1000 reached !!!!");
      return 0x0;
    }
    
    AliMUONSegment* si = static_cast<AliMUONSegment*>(all.UncheckedAt(i));
    inorder.Add(si);
    alreadyAdded[i] = 1;
    const AliMUONSegment* all0 = static_cast<const AliMUONSegment*>(all.First());
    if ( i != 0 && AliMUONSegment::AreEqual(si->EndX(),all0->StartX()) && AliMUONSegment::AreEqual(si->EndY(),all0->StartY()) )
    {
      Int_t n(-1);
      
      AliMUONPolygon polygon(inorder.GetLast()+2);

      // we got a cycle. Add it to the contour
      for ( Int_t j = 0; j <= inorder.GetLast(); ++j ) 
      {
        AliMUONSegment* s = static_cast<AliMUONSegment*>(inorder.UncheckedAt(j));
        polygon.SetVertex(++n,s->StartX(),s->StartY());
        all.Remove(s);
      }
      
      all.Compress();

      polygon.Close();
      
      contour->Add(polygon);
      
      if ( ! all.IsEmpty() )
      {
        i = 0;
        inorder.Clear();
        alreadyAdded.Set(all.GetLast()+1);
        alreadyAdded.Reset();
      }
      continue;
    }
    
    for ( Int_t j = 0; j <= all.GetLast(); ++j) 
    {
      if ( j != i && alreadyAdded[j] == 0 ) 
      {        
        const AliMUONSegment* sj = static_cast<const AliMUONSegment*>(all.UncheckedAt(j));
        if ( AliMUONSegment::AreEqual(si->EndX(),sj->StartX()) && AliMUONSegment::AreEqual(si->EndY(),sj->StartY()))
        {
          i = j;
          break;
        }
      }
    }    
  }
  
  contour->AssertOrientation(kTRUE);
  return contour;
}


//_____________________________________________________________________________
void 
AliMUONContourMaker::GetVerticalEdges(const TObjArray& polygons, TObjArray& polygonVerticalEdges) const
{
  /// From an array of polygons, extract the list of vertical edges.
  /// Output array polygonVerticalEdges should be empty before calling.
  
  AliCodeTimerAuto("",0);
  
  for ( Int_t i = 0; i <= polygons.GetLast(); ++i ) 
  {
    const AliMUONPolygon* g = static_cast<const AliMUONPolygon*>(polygons.UncheckedAt(i));
    for ( Int_t j = 0; j < g->NumberOfVertices()-1; ++j ) 
    {
      if ( AliMUONSegment::AreEqual(g->X(j),g->X(j+1)) ) // segment is vertical
      {
        polygonVerticalEdges.Add(new AliMUONSegment(g->X(j),g->Y(j),g->X(j+1),g->Y(j+1)));
      }
    }
  }
}


//_____________________________________________________________________________
void
AliMUONContourMaker::GetYPositions(const TObjArray& polygonVerticalEdges,
                                   TArrayD& yPositions) const
{
  /// Fill the array yPositions with the different y positions found in 
  /// polygonVerticalEdges
  
  AliCodeTimerAuto("",0);
  
  Double_t* y = new Double_t[polygonVerticalEdges.GetSize()*2];
  Int_t n(0);
  
  for ( Int_t i = 0; i < polygonVerticalEdges.GetLast(); ++i ) 
  {
    AliMUONSegment* s = static_cast<AliMUONSegment*>(polygonVerticalEdges.UncheckedAt(i));
    y[n] = s->StartY();
    y[n+1] = s->EndY();
    n += 2;
  }
  Int_t* ix = new Int_t[n+1];
  
  TMath::Sort(n,y,ix,kFALSE);
  
  yPositions.Set(n+1);
  
  Int_t u(0);
  Double_t x(FLT_MAX);
  
  for ( Int_t i = 0; i < n; ++i ) 
  {
    if ( y[ix[i]] != x )
    {
      yPositions[u] = y[ix[i]];
      x = y[ix[i]];
      ++u;
    }
  }

  yPositions.Set(u);
  
  delete[] ix;
  delete[] y;
  
}

//_____________________________________________________________________________
AliMUONContour* 
AliMUONContourMaker::MergeContour(const TObjArray& contours, const char* name) const
{
  /// Merge all the polygons of all contours into a single contour
  
  AliCodeTimerAuto("",0);
  
  TObjArray polygons;
  polygons.SetOwner(kTRUE);
  
  TIter next(&contours);
  AliMUONContour* contour;
  while ( ( contour = static_cast<AliMUONContour*>(next()) ) )
  {
    const TObjArray* contourPolygons = contour->Polygons();
    TIter nextPol(contourPolygons);
    AliMUONPolygon* pol;
    while ( ( pol = static_cast<AliMUONPolygon*>(nextPol()) ) )
    {
      polygons.Add(new AliMUONPolygon(*pol));
    }
  }
  
  if ( polygons.IsEmpty() ) return 0x0;
  
  contour = CreateContour(polygons,name);
  
  return contour;
}

//_____________________________________________________________________________
void 
AliMUONContourMaker::SortPoints(const TObjArray& polygonVerticalEdges, 
                                TObjArray& sortedPoints) const
{
  /// Sort the point of the vertical edges in ascending order, first on ordinate, 
  /// then on abcissa, and put them in output vector sortedPoints.
  /// Output array sortedPoints should be empty before calling this method.
  
  AliCodeTimerAuto("",0);
  
  for ( Int_t i = 0; i <= polygonVerticalEdges.GetLast(); ++i )
  {
    const AliMUONSegment* e = static_cast<const AliMUONSegment*>(polygonVerticalEdges.UncheckedAt(i));
    sortedPoints.Add(new AliMUONPointWithRef(e->StartX(),e->StartY(),i));
    sortedPoints.Add(new AliMUONPointWithRef(e->EndX(),e->EndY(),i));
    // note that we keep track of the original edge, which is used
    // later on to deduce orientation of horizontal edges.
  }
  
  sortedPoints.Sort(); 
}

//_____________________________________________________________________________
void 
AliMUONContourMaker::Sweep(const TObjArray& polygonVerticalEdges, 
                           TObjArray& contourVerticalEdges) const
{
  /// This is the meat of the algorithm of the contour merging...
  
  AliCodeTimerAuto("",0);
  
  TArrayD yPositions;
  GetYPositions(polygonVerticalEdges,yPositions);
  
  AliMUONSegmentTree segmentTree(yPositions);
  
  for ( Int_t i = 0; i <= polygonVerticalEdges.GetLast(); ++i )
  {
    const AliMUONSegment* edge = static_cast<const AliMUONSegment*>(polygonVerticalEdges.UncheckedAt(i));
    
    assert(edge!=0x0);
    
    if ( edge->IsLeftEdge() ) 
    {
      segmentTree.Contribution(edge->Bottom(),edge->Top());
      segmentTree.InsertInterval(edge->Bottom(),edge->Top());
    }
    else
    {
      segmentTree.DeleteInterval(edge->Bottom(),edge->Top());
      segmentTree.Contribution(edge->Bottom(),edge->Top());
    }
    
    AliMUONSegment e1(*edge);
    
    if ( i < polygonVerticalEdges.GetLast() ) 
    {
      const AliMUONSegment* next = static_cast<const AliMUONSegment*>(polygonVerticalEdges.UncheckedAt(i+1));
      e1 = *next;
    }

    if ( ( edge->IsLeftEdge() != e1.IsLeftEdge() ) ||
        ( !AliMUONSegment::AreEqual(edge->StartX(),e1.StartX() ) ) ||
        ( i == polygonVerticalEdges.GetLast() ) )
    {
      const TObjArray& stack = segmentTree.Stack();
      
      double x = edge->StartX();
      
      for ( Int_t j = 0; j <= stack.GetLast(); ++j )
      {
        AliMUONSegment* sj = static_cast<AliMUONSegment*>(stack.UncheckedAt(j));
        AliMUONSegment* s = new AliMUONSegment(x,sj->StartY(),x,sj->EndY());
        
        if  (s->IsAPoint()) 
        {
          delete s;
          continue;
        }
        
        if ( edge->IsLeftEdge() != s->IsLeftEdge() ) 
        {
          s->Set(x,sj->EndY(),x,sj->StartY());
        }
        contourVerticalEdges.Add(s);
      }
      segmentTree.ResetStack();
    }
  }
}

//_____________________________________________________________________________
void 
AliMUONContourMaker::VerticalToHorizontal(const TObjArray& polygonVerticalEdges,
                                          TObjArray& horizontalEdges) const
{
  /// Deduce the set of horizontal edges from the vertical edges
  /// Output array horizontalEdges should be empty before calling this method
  
  AliCodeTimerAuto("",0);
  
  TObjArray points; // array of AliMUONPointWithRef
  points.SetOwner(kTRUE);
  
  SortPoints(polygonVerticalEdges,points);
  
  for ( Int_t k = 0; k < (points.GetLast()+1)/2; ++k )
  {
    const AliMUONPointWithRef* p1 = static_cast<AliMUONPointWithRef*>(points.UncheckedAt(k*2));
    const AliMUONPointWithRef* p2 = static_cast<AliMUONPointWithRef*>(points.UncheckedAt(k*2+1));
    
    const AliMUONSegment* refEdge = static_cast<const AliMUONSegment*>(polygonVerticalEdges.UncheckedAt(p1->Ref()));
    
    // (p1,p2) is the horizontal edge.
    // refEdge is used to deduce the orientation of (p1,p2)
    
    if ( AliMUONSegment::AreEqual(p1->X(),refEdge->EndX()) && AliMUONSegment::AreEqual(p1->Y(),refEdge->EndY()) )
//    if ( AreEqual(p1,refEdge->End()) )
    {
      horizontalEdges.Add(new AliMUONSegment(p1->X(),p1->Y(),p2->X(),p2->Y()));
    }
    else
    {
      horizontalEdges.Add(new AliMUONSegment(p2->X(),p2->Y(),p1->X(),p1->Y()));
    }
  }
}

