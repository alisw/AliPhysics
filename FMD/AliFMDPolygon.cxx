/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//////////////////////////////////////////////////////////////////////////////
//
// Utility class to help implement shape of the FMD modules
//
// Latest changes by Christian Holm Christensen
//
//////////////////////////////////////////////////////////////////////////////
#include "AliFMDPolygon.h"	// ALIFMDPOLYGON_H
#include "AliLog.h"		// ALILOG_H
#include "TString.h"		// ROOT_TString
#include "TVector2.h"		// ROOT_TVector2
#include "TCanvas.h"		// ROOT_TCanvas
#include "TText.h"		// ROOT_TText
#include "TGraph.h"		// ROOT_TGraph
#include "TError.h"		// ROOT_TError

//____________________________________________________________________
ClassImp(AliFMDPolygon)

//____________________________________________________________________
AliFMDPolygon::AliFMDPolygon()
  : fState(kUnknown)
{}

//____________________________________________________________________
AliFMDPolygon::~AliFMDPolygon() 
{
  fVerticies.Delete();
}

//____________________________________________________________________
bool
AliFMDPolygon::AddVertex(double x, double y) 
{
  // Add points in order to the polygon. 
  TVector2* c = new TVector2(x,y);
  return AddVertex(c);
}
//____________________________________________________________________
bool
AliFMDPolygon::AddVertex(TVector2* c) 
{
  // Add points in order to the polygon.
  //
  // Checks if the addition of the vertex makes the alipolygon
  // concave, and if it does, returns false.  Note, that the check
  // isn't performed until at least 3 verticies have already been
  // added to the alipolygon (a triangle is always convex).
  if (!c) return false;
  if (fVerticies.GetEntries() >= 3) {
    if (!IsOnLeftHand(c, fVerticies.GetEntries()-2,
		      fVerticies.GetEntries()-1)) {
      Error("AliFMDPolygon::AddVertex", 
	    "Adding vertex (%f,%f) makes alipolygon concave", 
	    c->X(), c->Y());
      return false;
    }
  }
  fState = kUnknown;
  fVerticies.AddAtAndExpand(c, fVerticies.GetEntries());
  return true;
}

//____________________________________________________________________
const TVector2&
AliFMDPolygon::GetVertex(size_t i) const
{
  // Get one of the verticies
  if (i > size_t(fVerticies.GetEntries()))
    Fatal("AliFMDPolygon::GetVertex", "Index %d out of range [0,%d]", 
	  i, fVerticies.GetEntries());
  return *(static_cast<TVector2*>(fVerticies.At(i)));
}

//____________________________________________________________________
bool
AliFMDPolygon::Contains(const TVector2* c) const
{
  /* This function checks if a point is inside the polygon.  It does
   * that by looping over the segments of the polygon, and decides
   * wether the point is on the left-hand-side (LHS) of the segment.
   * 
   * Suppose we had the polygon and point 
   *
   *     2 x----x 1
   *      /      \\
   *    3 x   P   x 0
   *      \\      /
   *     4 x----x 5
   * 
   * Then, P is on LHS of the segment 0->1, 1->2, ..., and 5->0, and
   * so inside the polygon.
   * 
   * Suppose instead the point was like
   * 
   *     2 x----x 1
   *      /      \\    P
   *    3 x       x 0
   *      \\      /
   *     4 x----x 5
   * 
   * Then it would still be on the LHS of the segments 1->2, 2->3,
   * 3->4, 4->5, but on the right-hand-side of 5->0, and 0->1 and
   * hence outside the polygon.
   */ 
  if (!c) return kFALSE;
  if (fState == kUnknown) fState = (ConvexCheck() ? kConvex : kConcave);
  if (fState == kConcave) return false;
  size_t n = fVerticies.GetEntries();
  for (size_t i = 0; i < n; ++i) 
    if (!IsOnLeftHand(c, i, (i + 1) % n))
      // If point is on the left hand side of the segment, it's
      // outside the polygon.
      return false;
  return true;
}

//____________________________________________________________________
bool
AliFMDPolygon::Contains(double x, double y) const
{
  TVector2 c(x,y);
  return Contains(&c);
}

//____________________________________________________________________
bool
AliFMDPolygon::ConvexCheck() const
{
  /* 
   * This member function loops through the verticies, and checks if
   * the next-to-next vertex is in between the current vertex and the
   * next vertex.
   * 
   * Suppose we have the polygon 
   *           
   *     2 x----x 1
   *      /      \\
   *    3 x       x 0
   *      \\      /
   *     4 x----x 5
   * 
   * Then the vertex 2 is on the left-hand-side (LHS) of the segment
   * 0->1, 3 is on the LHS of 1->2, ..., and 0 is on the RHS of 4->5,
   * and so the polygon is convex.
   * 
   * If we had a polygon like
   *
   *     2 x----x 1
   *        \\    \\
   *       3 x    x 0
   *        /     /
   *     4 x----x 5
   * 
   * Then, the vertex 4 is NOT on the LHS of the segment 2->3, and so
   * the polygon is NOT convex.
   */
  size_t n = fVerticies.GetEntries();
  // A triangle is always convex. 
  if (n <= 3) return true;
  for (size_t i = 0; i < n; i++) {
    // Get the next, and next-to-next indicies, taking care to go
    // around the polygon.
    size_t j = (i + 1) % n;
    size_t k = (i + 2) % n;
    // Check The next-to-next vertex is on the LHS of the current
    // vertex and the next vertex 
    if (!IsOnLeftHand(static_cast<TVector2*>(fVerticies.At(k)), i, j)) {
      Error("AliFMDPolygon::ConvexCheck", 
	    "AliFMDPolygon is concave at segment %d -> %d" , i, j);
      return false;
    }
  }
  return true;
}
    
//____________________________________________________________________
bool
AliFMDPolygon::IsOnLeftHand(const TVector2* c, size_t i1, size_t i2) const
{
  /* Check if a point is on the left-hand-side (LHS) or
   * right-hand-side (RHS) of the segment defined by the two indicies.
   * 
   * Suppose we have the segment and point 
   * 
   *     1 x
   *        \\    P
   *         x 0  
   *
   * Then, we define the vectors 
   *
   *   v: 0->P 
   *   u: perpendicular to 1->0 
   * 
   * The dot product of u and v is >= 0, meaning that the
   * angle between the two vectors is in [0,pi], and so P is on the
   * RHS of the segment
   *
   * Suppose we had 
   *
   *       1 x
   *    P     \\    
   *           x 0  
   *
   * And defining the vectors u,v as above.  In this case the dot
   * product is less than 0, meaning the angle between u,v is in
   * (pi,2pi], and so P is on the LHS of the segment
   */
  if (!c) return false;
  const TVector2& v1 = GetVertex(i1);
  const TVector2& v2 = GetVertex(i2);
  double dot = (  (c->X()  - v1.X())  * (v2.Y() - v1.Y()) 
		- (c->Y() - v1.Y()) * (v2.X()  - v1.X()));
  return (dot < 0 ? true : false);
}

//____________________________________________________________________
void
AliFMDPolygon::Clear(Option_t*) 
{
  fState = kUnknown;
  fVerticies.Delete();
  // fVerticies.Resize(0);
}

//____________________________________________________________________
void
AliFMDPolygon::Draw(const char* option, const char* name) const 
{
  size_t n = fVerticies.GetEntries();
  TGraph* g = new TGraph(n+1);
  if (name) g->SetName(name);
  g->SetMarkerStyle(20);
  for (size_t i = 0; i < n; i++) {
    const TVector2& v = GetVertex(i);
    g->SetPoint(i, v.X(), v.Y());
  }
  const TVector2& v = GetVertex(0);
  g->SetPoint(n, v.X(), v.Y());
  g->Draw(option);
  TString opt(option);
  if (!opt.Contains("p", TString::kIgnoreCase))
    return;
  for (size_t i = 0; i < n; i++) {
    const TVector2& v = GetVertex(i);
    TText* t = new TText(v.X(), v.Y(), Form("%d", i));
    t->Draw();
  }
}

//____________________________________________________________________
//
// EOf
//
