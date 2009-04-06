#ifndef ALIMUONPAINTERHELPER_H
#define ALIMUONPAINTERHELPER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterHelper
/// \brief Utility class for the painters display
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ROOT_TString
#  include "TString.h"
#endif

#ifndef ROOT_TVector2
#  include "TVector2.h"
#endif

#ifndef ALI_MP_PLANE_TYPE_H
#  include "AliMpPlaneType.h"
#endif

#ifndef ALI_MP_CATHOD_TYPE_H
#  include "AliMpCathodType.h"
#endif

#ifndef ALIMUONVPAINTER_H
#  include "AliMUONVPainter.h"
#endif

#ifndef ALI_MP_PAD_H
#  include "AliMpPad.h"
#endif

#ifndef ROOT_TMap
#  include "TMap.h"
#endif

class AliMUONAttPainter;
class AliMUONContour;
class AliMUONPainterEnv;
class AliMUONVCalibParam;
class AliMUONVTrackerData;
class AliMpExMap;
class AliMpMotifPosition;
class AliMpPCB;
class AliMpSector;
class AliMpSlat;
class TArrayI;
class TList;
class TMap;
class TObjArray;

class AliMUONPainterHelper : public TObject
{
public:
  AliMUONPainterHelper();
  virtual ~AliMUONPainterHelper();
    
  TString ChamberName(Int_t chamberId) const;
  TString StationName(Int_t stationId) const;
  TString DEName(Int_t detElemId) const;
  TString ManuName(Int_t manuId) const;
  TString BusPatchName(Int_t busPatchId) const;
  TString PCBName(Int_t pcbNumber) const;

  TString ChamberPathName(Int_t chamberId) const;
  TString StationPathName(Int_t stationId) const;
  TString DEPathName(Int_t detElemId) const;
  TString ManuPathName(Int_t detElemId, Int_t manuId) const;
  TString BusPatchPathName(Int_t busPatchId) const;
  TString PCBPathName(Int_t detElemId, Int_t pcbNumber) const;
  
  Int_t ColorFromValue(Double_t value, Double_t min, Double_t max) const;
  
  AliMp::CathodType GetCathodeType(Int_t detElemId, Int_t manuId) const;

  AliMUONContour* GenerateManuContour(Int_t detElemId, Int_t manuId,
                                             AliMUONAttPainter viewType,
                                             const char* contourName);

  AliMUONContour* GetContour(const char* contourName) const;

  /// Return a contour by name
  AliMUONContour* GetContour(const TString& contourName) const { return GetContour(contourName.Data()); }

  AliMpMotifPosition* GetMotifPosition(Int_t detElemId, Int_t manuId) const;
  
  AliMpPCB* GetPCB(Int_t detElemId, AliMp::PlaneType planeType, 
                   Int_t pcbNumber) const;

  AliMpPCB* GetPCB(Int_t detElemId, AliMp::CathodType cathodType,
                   Int_t pcbNumber) const;

  AliMp::PlaneType GetPlaneType(Int_t manuId) const;
  
  const AliMpSector* GetSector(Int_t detElemId, AliMp::PlaneType planeType) const;
    
  const AliMpSlat* GetSlat(Int_t detElemId, AliMp::PlaneType planeType) const;

  const AliMpSlat* GetSlat(Int_t detElemId, AliMp::CathodType cathodeType) const;

  const AliMpSlat* GetSlat(Int_t detElemId, Int_t manuId) const;

  static AliMUONPainterHelper* Instance();

  AliMpPad PadByExplodedPosition(Int_t detElemId, Int_t manuId, Double_t x, Double_t y) const;
  
  void Exploded2Real(Int_t detElemId, Double_t xe, Double_t ye, Double_t ze,
                     Double_t& xr, Double_t& yr, Double_t& zr) const;
  
  void Local2Global(Int_t detElemId, Double_t xl, Double_t yl, Double_t zl,
                    Double_t& xg, Double_t& yg, Double_t& zg) const;

  void Local2GlobalReal(Int_t detElemId, Double_t xl, Double_t yl, Double_t zl,
                        Double_t& xg, Double_t& yg, Double_t& zg) const;

  void Global2Local(Int_t detElemId, Double_t xg, Double_t yg, Double_t zg,
                     Double_t& xl, Double_t& yl, Double_t& zl) const;

  void Global2LocalReal(Int_t detElemId, Double_t xg, Double_t yg, Double_t zg,
                        Double_t& xl, Double_t& yl, Double_t& zl) const;

  AliMUONContour* MergeContours(const TObjArray& contours, 
                                       const char* contourName);
  
  virtual void Print(Option_t* opt="") const;
  
  void RegisterContour(AliMUONContour* contour);
  
  TString FormatValue(const char* name, Double_t value) const;
  
  /// Return the environment
  AliMUONPainterEnv* Env() { return fEnv; }
  
private:
    
  /// Not implemented
  AliMUONPainterHelper(const AliMUONPainterHelper&);
  /// Not implemented
  AliMUONPainterHelper& operator=(const AliMUONPainterHelper&);
  
  void GenerateDefaultMatrices();
  void GenerateGeometry();

private:
  static AliMUONPainterHelper* fgInstance; ///< global instance
  
  Double_t fExplodeFactor[2]; ///< explosing factors for representation
  AliMpExMap* fExplodedGlobalTransformations; ///< global geometric transformations (exploded)
  AliMpExMap* fRealGlobalTransformations; ///< global geometric transformations (real)
  TObjArray* fPainterMatrices; ///< default matrices
  AliMUONPainterEnv* fEnv; ///< resources
  TMap fAllContours; ///< all contours
  
  ClassDef(AliMUONPainterHelper,2) // Helper class for painters
};

#endif
