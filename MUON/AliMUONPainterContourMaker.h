#ifndef ALIMUONPAINTERCONTOURMAKER_H
#define ALIMUONPAINTERCONTOURMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterContourMaker
/// \brief Utility class to build painter contours
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TVector2
#  include "TVector2.h"
#endif

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ALIMUONVPAINTER_H
#  include "AliMUONVPainter.h"
#endif

class AliMpExMap;
class AliMpMotifPosition;
class AliMUONAttPainter;
class AliMUONPainterContour;
class AliMUONPainterPadStore;
class TArrayI;
class TMap;
#include <TLine.h>
class TPolyLine;
class TObjArray;

class AliMUONPainterContourMaker : public TObject
{
public:
  AliMUONPainterContourMaker(AliMpExMap* globalTransformations=0x0);
  virtual ~AliMUONPainterContourMaker();
  
  void Add(AliMUONPainterContour* contour);
  
  AliMUONPainterContour* FindLocalManuContour(Int_t detElemId, Int_t manuId) const;
  
  AliMUONPainterContour* GetContour(const char* name) const;

  AliMUONPainterContour* GenerateManuContour(const char* name,
                                             Int_t detElemId,
                                             Int_t manuId,
                                             AliMUONAttPainter viewType) const;

  Bool_t HasContour(const char* name) const;

  AliMUONPainterContour* MergeContours(const TObjArray& contours,
                                       const char* contourName) const;
  
  Int_t Size() const;
  
  void Print(Option_t* opt="") const;
  
public:
    
  /// \ingroup graphics
  /// \brief Store information about one pad's neighbours.
  /// \author Laurent Aphecetche, Subatech
  class AliMUONNeighbour : public TObject
  {
public:
    /// default ctor
    AliMUONNeighbour()
    : fID(-1), fPosition(), fDimensions(), 
    fLeft(kFALSE), fRight(kFALSE), 
    fTop(kFALSE), fBottom(kFALSE) {}
    
    /// normal ctor
    AliMUONNeighbour(Int_t absID, 
              const TVector2& position, 
              const TVector2& dimensions,
              Bool_t hasLeftNeighbour, Bool_t hasRightNeighbour,
              Bool_t hasTopNeighbour, Bool_t hasBottomNeighbour)
    : fID(absID), fPosition(position), fDimensions(dimensions), 
    fLeft(hasLeftNeighbour), fRight(hasRightNeighbour), 
      fTop(hasTopNeighbour), fBottom(hasBottomNeighbour) {}

    /// dtor
    virtual ~AliMUONNeighbour() {}
    
    /// we are sortable
    virtual Bool_t IsSortable() const { return kTRUE; }
    
    virtual Int_t Compare(const TObject* object) const;
    
    /// our id
    Int_t ID() const { return fID; }
    
    /// Whether we have a neighbour on our left
    Bool_t HasLeftNeighbour() const { return fLeft; }
    /// Whether we have a neighbour on our right
    Bool_t HasRightNeighbour() const { return fRight; }
    /// Whether we have a neighbour above
    Bool_t HasTopNeighbour() const { return fTop; }
    /// Whether we have a neighbour below
    Bool_t HasBottomNeighbour() const { return fBottom; }
    
    /// Our position
    TVector2 Position() const { return fPosition; }
    /// Our (half-)dimensions
    TVector2 Dimensions() const { return fDimensions; }
    
    /// Lower left corner
    TVector2 LowerLeft() const { return fPosition - fDimensions; }
    /// Upper right corner
    TVector2 UpperRight() const { return fPosition + fDimensions; }
    
    void Print(Option_t* opt="") const;
    
    /// Set our position
    void SetPosition(Double_t x, Double_t y) { fPosition.Set(x,y); }
    
private:
    Int_t fID; ///< id of the pad
    TVector2 fPosition; ///< position
    TVector2 fDimensions; ///< (half)dimension
    Bool_t fLeft; ///< do we have a neighbour on our left ?
    Bool_t fRight; ///< do we have a neighbour on our right ?
    Bool_t fTop; ///< do we have a neighbour on top of us ?
    Bool_t fBottom; ///< do we have a neighbour below us ?
    
    ClassDef(AliMUONNeighbour,1) // Neighbour internal class
  };
    

public:
    

  void AddSegments(TObjArray& segments, const AliMUONPainterContour& contour) const;

  void AddSegment(TObjArray& segments, Double_t x1, Double_t y1,
                    Double_t x2, Double_t y2, Int_t padID) const;
  
  TLine* AddToLine(TPolyLine& line, TObjArray& segments, Int_t i) const;

  AliMpMotifPosition* FindMotifPosition(Int_t detElemId, Int_t manuId) const;
  
  Int_t FindPoint(const TPolyLine& lines, Double_t x, Double_t y) const;
  
  Int_t FindPoint(Double_t x, Double_t y, TObjArray& segments) const;

  TLine* InsertSegment(TPolyLine& lines, TLine& l) const;

  using TObject::IsEqual;
  
  Bool_t IsEqual(Double_t x, Double_t y) const;

  Int_t Overlap(const TLine& line1, const TLine& line2) const;

  Bool_t IsLineClosed(const TPolyLine& line) const;
  
  AliMUONPainterContour* ConvertEdgePadsToContour(TObjArray& edgePads, const char* name) const;
  
  AliMUONPainterContour* ConvertSegmentsToContour(TObjArray& segments, const char* name) const;
  
  void Local2Global(Int_t detElemId, Double_t xl, Double_t yl, Double_t zl,
                    Double_t& xg, Double_t& yg, Double_t& zg) const;
    
  TString NameIt(const AliMpMotifPosition& motifPosition) const;

  TPolyLine* Simplify(const TPolyLine& lines) const;

  Double_t Slope(const TLine& line) const;

  Bool_t IsPoint(const TLine& line) const;

  void PrintLine(const TLine& line, const char* msg="") const;

  void PrintSegments(const TObjArray& segments) const;
    
  Bool_t SameDirection(const TLine& line1, const TLine& line2) const;

  void Swap(TLine& line) const;
  
  TLine Shift(const TLine& line, Double_t x, Double_t y) const;

  Int_t IsInside(const TLine& line1, const TLine& line2,
                 Bool_t useEndPoints=kFALSE) const;
    
  Bool_t IsInside(const TObjArray& segments, const TLine& line) const;
    
  Int_t IsInLine(const TLine& line, Double_t x, Double_t y, 
                 Bool_t strict=kTRUE) const;
  
  Bool_t IsEqual(const TLine& line1, const TLine& line2) const;
    
  Bool_t SanityCheck(const TObjArray& contours, const TObjArray& segments, Bool_t check=kTRUE) const;
  
  TString LineAsString(const TLine& line, Bool_t slope=kTRUE) const;
  
  Int_t IsInRange(Double_t x, Double_t a, Double_t b, 
                  Bool_t strict=kTRUE) const;
  
  Bool_t HasLine(const TObjArray& segments, const TLine& line) const;
    
  void CleanSegments(TObjArray& segments, const TArrayI& toBeRemoved) const;
    
  Int_t SplitSegments(TObjArray& segments) const;

  Int_t RemoveInsideSegments(const TObjArray& contours, TObjArray& segments) const;
      
  Bool_t IsHorizontal(const TLine& line) const;

  Bool_t IsVertical(const TLine& line) const;
    
  Int_t CountPoint(const TObjArray& segments, Double_t x, Double_t y) const;
    
  Bool_t ShouldBeRemoved(const TObjArray& contours, Double_t x, Double_t y) const;
    
private:
  /// not implemented
  AliMUONPainterContourMaker(const AliMUONPainterContourMaker& rhs);
  /// not implemented
  AliMUONPainterContourMaker& operator=(const AliMUONPainterContourMaker& rhs);

  AliMpExMap* fGlobalTransformations; ///< store of global transformations for DEs
  TMap* fLocalManuContours; ///< store for local contours of all manus
  TMap* fContours; ///< store for all our contours
  
  ClassDef(AliMUONPainterContourMaker,1) // Painter contour builder
};

#endif

