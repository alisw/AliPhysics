#ifndef ALIMUONCONTOURMAKER_H
#define ALIMUONCONTOURMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONContourMaker
/// \brief Creator/merger of AliMUONContour objects
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ROOT_TMap
#include "TMap.h"
#endif

class AliMUONContour;
class TObjArray;
class TArrayD;

class AliMUONContourMaker : public TObject
{
public:
  AliMUONContourMaker();
  virtual ~AliMUONContourMaker();
  
  AliMUONContour* CreateContour(const TObjArray& polygons, const char* name=0x0) const;

  AliMUONContour* MergeContour(const TObjArray& contours, const char* name=0x0) const;

private:
  
  AliMUONContour* FinalizeContour(const TObjArray& verticals, const TObjArray& horizontals) const;
  
  void GetYPositions(const TObjArray& polygonVerticalEdges, TArrayD& yPositions) const;
  
  void GetVerticalEdges(const TObjArray& polygons, TObjArray& polygonVerticalEdges) const;

  void SortPoints(const TObjArray& polygonVerticalEdges, TObjArray& sortedPoints) const;
  
  void Sweep(const TObjArray& polygonVerticalEdges, TObjArray& contourVerticalEdges) const;

  void VerticalToHorizontal(const TObjArray& verticalEdges, TObjArray& horizontalEdges) const;
  
  ClassDef(AliMUONContourMaker,1) // Maker/merger of AliMUONContour objects
};

#endif
