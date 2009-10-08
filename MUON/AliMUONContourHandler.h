#ifndef ALIMUONCONTOURHANDLER_H
#define ALIMUONCONTOURHANDLER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONContourHandler
/// \brief Holder for MUON tracker contours
/// 
/// \author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TObjArray;
class AliMpExMap;
class TMap;
class AliMUONContour;

class AliMUONContourHandler : public TObject
{
public:
  AliMUONContourHandler(Bool_t explodedView=kTRUE);
  virtual ~AliMUONContourHandler();
  
  Bool_t Adopt(AliMUONContour* contour); 
  
  TMap* AllContourMap() const { return fAllContourMap; }
  
  TObjArray* AllContourArray() const { return fAllContourArray; }
  
  AliMUONContour* GetContour(const char* contourname) const;
  
  AliMpExMap* GetTransformations() const { return fTransformations; }
  
  void Print(Option_t* opt="") const;
  
private:
  
  AliMUONContourHandler(const AliMUONContourHandler& rhs);
  AliMUONContourHandler& operator=(const AliMUONContourHandler& rhs);
  
  AliMpExMap* GenerateTransformations(Bool_t exploded);
  
  TObjArray* CreateContourList(const TObjArray& manuContours);
  
  void GenerateAllContours(const TObjArray& manuContours);  
  
private:
  AliMpExMap* fTransformations; //< transformations used to go from local to global coordinates
  TMap* fAllContourMap; //< all (i.e. manus,  buspatches, detection elements, etc..) contours
  TObjArray* fAllContourArray; //< all contours, ordered by hierarchy level
  
  ClassDef(AliMUONContourHandler,1) // MUON tracker contour holder
};

#endif
