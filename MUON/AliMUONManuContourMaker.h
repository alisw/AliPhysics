#ifndef ALIMUONMANUCONTOURMAKER_H
#define ALIMUONMANUCONTOURMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONManuContourMaker
/// \brief Maker of AliMUONContour objects for all the tracker manus
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ROOT_TMap
#  include "TMap.h"
#endif

class AliMpExMap;
class AliMpMotifPosition;
class AliMUONContour;

class AliMUONManuContourMaker : public TObject
{
public:
  AliMUONManuContourMaker(AliMpExMap* deTransformations);
  virtual ~AliMUONManuContourMaker();

  AliMUONContour* CreateManuContour(Int_t detElemId, Int_t manuId, const char* name="") const;
  
  AliMUONContour* CreateMotifContour(const AliMpMotifPosition& motifPosition) const;

  TObjArray* GenerateManuContours(Bool_t stopAtError=kFALSE);
  
  static TString ManuPathName(Int_t detElemId, Int_t manu, Bool_t withCathodeName=kTRUE);
  
private:
  /// not implemented
  AliMUONManuContourMaker(const AliMUONManuContourMaker& rhs);
  /// not implemented
  AliMUONManuContourMaker& operator=(const AliMUONManuContourMaker& rhs);
  
  TString NameIt(const AliMpMotifPosition& motifPosition) const;

private:
  AliMpExMap* fDETransformations; //< map<int,TGeoHMatrix> of detElemId to matrix
  mutable TMap fLocalManuContours; //< map of local manu contours
  
  ClassDef(AliMUONManuContourMaker,1) // 
};

#endif
