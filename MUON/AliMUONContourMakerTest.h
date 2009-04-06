#ifndef ALIMUONCONTOURMAKERTEST_H
#define ALIMUONCONTOURMAKERTEST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup evaluation
/// \class AliMUONContourMakerTest
/// \brief Test of ContourMaker classes
/// 
// author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMpExMap;
class TObjArray;
class TString;
class AliMpMotifPosition;
class AliMUONContour;
class AliMUONPolygon;

class AliMUONContourMakerTest : public TObject
{
public:
  AliMUONContourMakerTest();
  virtual ~AliMUONContourMakerTest();
  
  void Exec(const Option_t* opt="ALL");

  void GetBoundingBox(const TObjArray& array, 
                      Double_t& xmin, Double_t& ymin, 
                      Double_t& xmax, Double_t& ymax,
                      Bool_t enlarge=kFALSE) const;
    
  void Plot(const AliMUONContour& contour, Int_t lineColor=5, Int_t lineWidth=4, Bool_t orientation=kFALSE) const;

  void Plot(const AliMUONPolygon& polygon, Int_t lineColor=5, Int_t lineWidth=4, Bool_t orientation=kFALSE) const;

  void PlotContours(const TObjArray& array, Bool_t orientations=kFALSE) const;

  void PlotSegments(const TObjArray& segments, Int_t lineColor=1, Int_t lineWidth=2, Bool_t orientations=kFALSE) const;
  
  void PrintAsPNG(const char* basename, const TObjArray& contourArray,
                  const TObjArray* contourVerticalEdges=0x0, 
                  const TObjArray* horizontals=0x0) const;
  
private:

  TObjArray* CreateContourList(const TObjArray& manuContours);

  TObjArray* GenerateAllContours(const TObjArray& manuContours);
  
  void GenerateTransformations(AliMpExMap*& real, AliMpExMap*& exploded);

  TString NameIt(const AliMpMotifPosition& motifPosition) const;

  ClassDef(AliMUONContourMakerTest,1) // Test of AliMUONContourMaker
};

#endif
