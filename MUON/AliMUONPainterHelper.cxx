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

#include <cstdlib>
#include "AliMUONPainterHelper.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONContour.h"
#include "AliMUONContourMaker.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONManuContourMaker.h"
#include "AliMUONPainterEnv.h"
#include "AliMUONPainterMatrix.h"
#include "AliMUONPainterRegistry.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDigit.h"
#include "AliMUONVTrackerData.h"
#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpExMap.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpPCB.h"
#include "AliMpPad.h"
#include "AliMpSector.h"
#include "AliMpSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpStationType.h"
#include "AliMpVPadIterator.h"
#include <Riostream.h>
#include <TArrayI.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TCollection.h>
#include <TFile.h>
#include <TGLabel.h>
#include <TGMsgBox.h>
#include <TGeoMatrix.h>
#include <TLine.h>
#include <TList.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TVirtualPad.h>
#include <TVirtualX.h>

#include "AliMUONChamberPainter.h"

///\class AliMUONPainterHelper
///
/// Helper class for painters
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONPainterHelper)
///\endcond

AliMUONPainterHelper* AliMUONPainterHelper::fgInstance(0x0);

//_____________________________________________________________________________
AliMUONPainterHelper::AliMUONPainterHelper() : 
TObject(),
fExplodedGlobalTransformations(0x0),
fRealGlobalTransformations(0x0),
fPainterMatrices(0x0),
fEnv(0x0),
fAllContours(20000,1)
{
    /// ctor
  
  fAllContours.SetOwnerKeyValue(kTRUE,kTRUE);
  
    fExplodeFactor[0] = 1.00;
    fExplodeFactor[1] = 1.50;

    if ( ! AliMpCDB::LoadMpSegmentation() ) 
    {
      AliFatal("Could not access mapping from OCDB !");
    }
    
    // Load DDL store
    if ( ! AliMpCDB::LoadDDLStore() ) 
    {
      AliFatal("Could not access DDL Store from OCDB !");
    }        
}

//_____________________________________________________________________________
AliMUONPainterHelper::~AliMUONPainterHelper()
{
  /// dtor
  delete fExplodedGlobalTransformations;
  delete fRealGlobalTransformations;
  delete fPainterMatrices;
  fEnv->Save();
  fgInstance = 0;
}

//_____________________________________________________________________________
AliMUONContour*
AliMUONPainterHelper::GetContour(const char* contourName) const
{
  /// Get a contour by name  
  return static_cast<AliMUONContour*>(fAllContours.GetValue(contourName));
}

//_____________________________________________________________________________
void
AliMUONPainterHelper::GenerateDefaultMatrices()
{
  /// Kind of bootstrap method to trigger the generation of all contours
  
  AliCodeTimerAuto("");
  
  fPainterMatrices = new TObjArray;
  fPainterMatrices->SetOwner(kFALSE);
  
  TObjArray attributes;
  
  AliMUONAttPainter att;
  
  att.SetViewPoint(kTRUE,kFALSE);
  att.SetPlane(kFALSE,kFALSE);
  att.SetCathode(kTRUE,kFALSE);

  att.SetCathode(kFALSE,kFALSE);
  att.SetPlane(kTRUE,kFALSE);
  attributes.Add(new AliMUONAttPainter(att));

// commented the lines below, as it's now fast enough to be created on the fly,
// only when requested by the user
//
//  attributes.Add(new AliMUONAttPainter(att));  
//  att.SetCathode(kFALSE,kTRUE);
//  attributes.Add(new AliMUONAttPainter(att));
//  att.SetPlane(kFALSE,kTRUE);
//  attributes.Add(new AliMUONAttPainter(att));
  
  TIter next(&attributes);
  AliMUONAttPainter* a;
  
  while ( ( a = static_cast<AliMUONAttPainter*>(next()) ) )
  {
    AliMUONPainterMatrix* matrix = new AliMUONPainterMatrix("Tracker",5,2);
    
    for ( Int_t i = 0; i < 10; ++i )
    {
      AliMUONVPainter* painter = new AliMUONChamberPainter(*a,i);
      
      painter->SetResponder("Chamber");
      
      painter->SetOutlined("*",kFALSE);
      
      painter->SetOutlined("MANU",kTRUE);
      
      for ( Int_t j = 0; j < 3; ++j ) 
      {
        painter->SetLine(j,1,4-j);
      }
      
      matrix->Adopt(painter);    
    }
    AliMUONPainterRegistry::Instance()->Register(matrix);
    fPainterMatrices->Add(matrix);
  }
}

//_____________________________________________________________________________
void
AliMUONPainterHelper::GenerateGeometry()
{  
  /// Generate the geometry (FIXME: using transform.dat for the moment)
  /// The geometry is not the "normal" one as we "explode" it to avoid
  /// having overlapping detection elements as in the reality, which 
  /// would be inconvenient for a display ;-)
  
  AliCodeTimerAuto("");
  
  AliMUONGeometryTransformer transformer;
  transformer.LoadGeometryData("transform.dat");
  fExplodedGlobalTransformations = new AliMpExMap;
  fRealGlobalTransformations = new AliMpExMap;
  AliMpDEIterator deIt;
  deIt.First();
  while ( !deIt.IsDone() )
  {
    Int_t detElemId = deIt.CurrentDEId();
    const AliMUONGeometryDetElement* de = transformer.GetDetElement(detElemId);
    
    fRealGlobalTransformations->Add(detElemId,de->GetGlobalTransformation()->Clone());
                                    
    TGeoHMatrix* matrix = static_cast<TGeoHMatrix*>(de->GetGlobalTransformation()->Clone());
    Double_t* translation = matrix->GetTranslation();
    
    AliDebug(1,Form("Initial translation for DE %04d is %7.3f, %7.3f",
                    detElemId,translation[0],translation[1]));
    
    if ( AliMpDEManager::GetStationType(detElemId) == AliMp::kStation345 ) 
    {
      translation[0] *= fExplodeFactor[0];
      translation[1] *= fExplodeFactor[1];
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
    fExplodedGlobalTransformations->Add(detElemId,matrix);
    deIt.Next();
  }
}

//_____________________________________________________________________________
AliMUONContour* 
AliMUONPainterHelper::GenerateManuContour(Int_t detElemId,
                                          Int_t manuId,
                                          AliMUONAttPainter viewType,
                                          const char* contourName)
{
  /// Generate the contour of the list of pads
  
  static AliMUONManuContourMaker maker(fExplodedGlobalTransformations);
  
  if ( viewType.IsBackView() )
  {
    AliError("Backview not implemented yet (and will probably never be, after all...)");
    return 0x0;
  }
  
  AliMUONContour* contour = maker.CreateManuContour(detElemId,manuId,contourName);
  
  if ( !contour ) return 0x0;
  
  AliMUONContour* pContour = new AliMUONContour(*contour);
  
  if (pContour) 
  {
    RegisterContour(pContour);
  }
  
  return pContour;
}

//_____________________________________________________________________________
AliMp::CathodType
AliMUONPainterHelper::GetCathodeType(Int_t detElemId, Int_t manuId) const
{
  /// Get the cathode type of a given manu
  
  AliMp::PlaneType planeType(AliMp::kBendingPlane);
  if ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) )
  {
    planeType = AliMp::kNonBendingPlane;
  }
  return AliMpDEManager::GetCathod(detElemId,planeType);
}


//_____________________________________________________________________________
AliMpMotifPosition* 
AliMUONPainterHelper::GetMotifPosition(Int_t detElemId, Int_t manuId) const
{
  /// Get a given motif position
  AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
  if ( stationType == AliMp::kStation345 ) 
  {
    AliMp::PlaneType planeType(AliMp::kBendingPlane);
    if ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) )
    {
      planeType = AliMp::kNonBendingPlane;
    }
    const AliMpSlat* slat = GetSlat(detElemId,planeType);
    return slat->FindMotifPosition(manuId);
  }
  else if ( stationType != AliMp::kStationTrigger ) 
  {
    AliMp::PlaneType planeType(AliMp::kBendingPlane);
    if ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) )
    {
      planeType = AliMp::kNonBendingPlane;
    }
    const AliMpSector* sector = GetSector(detElemId,planeType);
    return sector->GetMotifMap()->FindMotifPosition(manuId);
  }
  AliFatalClass("Not supposed to work with trigger");
  return 0x0;
}


//_____________________________________________________________________________
AliMpPCB*
AliMUONPainterHelper::GetPCB(Int_t detElemId, AliMp::CathodType cathodeType, 
                             Int_t pcbNumber) const
{
  /// Get a given PCB
  const AliMpSlat* slat = GetSlat(detElemId,cathodeType);
  if ( slat ) return slat->GetPCB(pcbNumber);
  return 0x0;
}

//_____________________________________________________________________________
AliMpPCB*
AliMUONPainterHelper::GetPCB(Int_t detElemId, AliMp::PlaneType planeType, 
                             Int_t pcbNumber) const
{
  /// Get a given PCB
  AliMp::CathodType cathodeType = AliMpDEManager::GetCathod(detElemId,
                                                            planeType);
  return GetPCB(detElemId,cathodeType,pcbNumber);
}

//_____________________________________________________________________________
AliMp::PlaneType
AliMUONPainterHelper::GetPlaneType(Int_t manuId) const
{
  /// Get the planeType of a given manu
  
  if ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) )
  {
    return AliMp::kNonBendingPlane;
  }
  return AliMp::kBendingPlane;
}

//_____________________________________________________________________________
const AliMpSlat*
AliMUONPainterHelper::GetSlat(Int_t detElemId, AliMp::PlaneType planeType) const
{
  /// Get a given slat
  
  AliMp::CathodType cathodeType = AliMpDEManager::GetCathod(detElemId,
                                                          planeType);

  return GetSlat(detElemId,cathodeType);
}

//_____________________________________________________________________________
const AliMpSector*
AliMUONPainterHelper::GetSector(Int_t detElemId, AliMp::PlaneType planeType) const
{
  /// Get a given sector
  AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
  if ( stationType != AliMp::kStation12 ) return 0x0;
  
  AliMp::CathodType cathodeType = AliMpDEManager::GetCathod(detElemId,planeType);
  
  return AliMpSegmentation::Instance()->GetSector(detElemId,cathodeType);
}

//_____________________________________________________________________________
const AliMpSlat*
AliMUONPainterHelper::GetSlat(Int_t detElemId, AliMp::CathodType cathodeType) const
{
  /// Get a given slat
  AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
  if ( stationType != AliMp::kStation345 ) return 0x0;

  return AliMpSegmentation::Instance()->GetSlat(detElemId,cathodeType);
}

//_____________________________________________________________________________
const AliMpSlat*
AliMUONPainterHelper::GetSlat(Int_t detElemId, Int_t manuId) const
{
  /// Get a given slat
  AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
  if ( stationType != AliMp::kStation345 ) return 0x0;

  return AliMpSegmentation::Instance()->GetSlatByElectronics(detElemId,manuId);
}

//_____________________________________________________________________________
AliMUONPainterHelper*
AliMUONPainterHelper::Instance()
{
  /// Return the global and unique instance of this class
  
  if (fgInstance) return fgInstance;

  AliCodeTimerAutoClass("");

  fgInstance = new AliMUONPainterHelper;
  fgInstance->GenerateGeometry();
  fgInstance->GenerateDefaultMatrices();
  fgInstance->fEnv = new AliMUONPainterEnv;
  return fgInstance;
}

//_____________________________________________________________________________
void 
AliMUONPainterHelper::Global2Local(Int_t detElemId, 
                                    Double_t xg, Double_t yg, Double_t zg,
                                    Double_t& xl, Double_t& yl, Double_t& zl) const
{
  /// Local to global transformation of coordinates
  
  TGeoHMatrix* matrix = static_cast<TGeoHMatrix*>(fExplodedGlobalTransformations->GetValue(detElemId));
  Double_t pg[3] = { xg, yg, zg };
  Double_t pl[3] = { 0., 0., 0. };
  matrix->MasterToLocal(pg, pl);
  xl = pl[0];
  yl = pl[1];
  zl = pl[2];
}

//_____________________________________________________________________________
void 
AliMUONPainterHelper::Global2LocalReal(Int_t detElemId, 
                                       Double_t xg, Double_t yg, Double_t zg,
                                       Double_t& xl, Double_t& yl, Double_t& zl) const
{
  /// Local to global transformation of coordinates
  
  TGeoHMatrix* matrix = static_cast<TGeoHMatrix*>(fRealGlobalTransformations->GetValue(detElemId));
  Double_t pg[3] = { xg, yg, zg };
  Double_t pl[3] = { 0., 0., 0. };
  matrix->MasterToLocal(pg, pl);
  xl = pl[0];
  yl = pl[1];
  zl = pl[2];
}

//_____________________________________________________________________________
void 
AliMUONPainterHelper::Local2Global(Int_t detElemId, 
                                   Double_t xl, Double_t yl, Double_t zl,
                                   Double_t& xg, Double_t& yg, Double_t& zg) const
{
  /// Local to (exploded) global transformation of coordinates
  
  TGeoHMatrix* matrix = static_cast<TGeoHMatrix*>(fExplodedGlobalTransformations->GetValue(detElemId));
  Double_t pl[3] = { xl, yl, zl };
  Double_t pg[3] = { 0., 0., 0. };
  matrix->LocalToMaster(pl, pg);
  xg = pg[0];
  yg = pg[1];
  zg = pg[2];
}

//_____________________________________________________________________________
void 
AliMUONPainterHelper::Local2GlobalReal(Int_t detElemId, 
                                       Double_t xl, Double_t yl, Double_t zl,
                                       Double_t& xg, Double_t& yg, Double_t& zg) const
{
  /// Local to (real) global transformation of coordinates
  
  TGeoHMatrix* matrix = static_cast<TGeoHMatrix*>(fRealGlobalTransformations->GetValue(detElemId));
  Double_t pl[3] = { xl, yl, zl };
  Double_t pg[3] = { 0., 0., 0. };
  matrix->LocalToMaster(pl, pg);
  xg = pg[0];
  yg = pg[1];
  zg = pg[2];
}

//_____________________________________________________________________________
Int_t
AliMUONPainterHelper::ColorFromValue(Double_t value, Double_t min, Double_t max) const
{ 
  /// Convert a value into a color, fitting within a given range
  
  Int_t rv;
  
  if (value > max) rv = 1;
  else if (value <= min) rv = 0;
  else
  {
    if  ( max == min ) return gStyle->GetColorPalette(1);
    Double_t range = max - min;
    Double_t offset = value - min;
    rv = gStyle->GetColorPalette( 1 + int( offset*(gStyle->GetNumberOfColors()-2)/range - 0.5 ) );
  }
  return rv;
}

//_____________________________________________________________________________
AliMUONContour* 
AliMUONPainterHelper::MergeContours(const TObjArray& contours, 
                                    const char* contourName)
{
  /// Merge a set of contours (delegating to the contour maker)
  
  AliMUONContourMaker maker;
  
  AliMUONContour* contour = maker.MergeContour(contours,contourName);
  
  if (contour) 
  {
    RegisterContour(contour);
  }
  return contour;
}


//_____________________________________________________________________________
void
AliMUONPainterHelper::Print(Option_t* opt) const
{
  /// Printout
  TString sopt(opt);
  sopt.ToUpper();
  
  if ( sopt.Length() == 0 )
  {
    cout << Form("ExplodeFactor=%e,%e",fExplodeFactor[0],fExplodeFactor[1]) << endl;
    cout << endl;
    cout << Form("GlobalTransformations=%x",fExplodedGlobalTransformations);
    if ( fExplodedGlobalTransformations ) cout << Form(" with %d transformations",fExplodedGlobalTransformations->GetSize());
    cout << endl;
    cout << Form("Contour map : collisions = %5.3f size = %d capacity = %d", 
                 fAllContours.AverageCollisions(),
                 fAllContours.GetSize(),
                 fAllContours.Capacity()) << endl;
    cout << endl;
  }
  
  if ( sopt.Contains("MATRI") || sopt.Contains("FULL") )
  {
    fPainterMatrices->Print(opt);
  }
}

//_____________________________________________________________________________
void
AliMUONPainterHelper::RegisterContour(AliMUONContour* contour)
{
  /// contour is adopted by contourMaker
  AliCodeTimerAuto("")
  AliDebug(1,contour->GetName());
  if ( fAllContours.GetValue(contour->GetName()) ) 
  {
    AliError(Form("Contour with name %s is already there",contour->GetName()));
//    Print("CONTOUR");
    return;
  }
  fAllContours.Add(new TObjString(contour->GetName()),contour);
}

//_____________________________________________________________________________
AliMpPad 
AliMUONPainterHelper::PadByExplodedPosition(Int_t detElemId, Int_t manuId, 
                                            Double_t x, Double_t y) const
{
  /// Find a pad by exploded position. FIXME: not really used nor tested !
  
  Double_t xr, yr, zr;
  
//  Local2Global(detElemId,0.0,0.0,0.0,dummy,dummy,z); // to find z 

  AliDebug(1,Form("DE %04d ManuID %04d x %7.3f y %7.3f",detElemId,manuId,x,y));
  
  Exploded2Real(detElemId,x,y,0,xr,yr,zr);

  AliDebug(1,Form("xr %7.3f yr %7.3f zr %7.3f",xr,yr,zr));

  Double_t xl,yl,zl;

  Global2LocalReal(detElemId,xr,yr,zr,xl,yl,zl);

  AliDebug(1,Form("xl %7.3f yl %7.3f zl %7.3f",xl,yl,zl));

  const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
  
  AliDebug(1,Form("dx,dy=%7.3f,%7.3f",seg->GetDimensionX(),seg->GetDimensionY()));
  
  return seg->PadByPosition(xl,yl);
}

//_____________________________________________________________________________
void 
AliMUONPainterHelper::Exploded2Real(Int_t detElemId, 
                                    Double_t xe, Double_t ye, Double_t ze, 
                                    Double_t& xr, Double_t& yr, Double_t& zr) const
{
  /// Convert exploded coordinates into real ones. FIXME: not really used nor tested !
  
  // first go back to local
  
  Double_t xl,yl,zl;
  
  Global2Local(detElemId,xe,ye,ze,xl,yl,zl);
  
  // and then back to global but not exploded
  
  Local2GlobalReal(detElemId,xl,yl,zl,xr,yr,zr);
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::ChamberName(Int_t chamberId) const
{
  /// Build a name for one chamber
  return Form("Chamber%1d",chamberId);
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::StationName(Int_t stationId) const
{
  /// Build a name for one station
  return Form("Station%1d",stationId+1);
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::DEName(Int_t detElemId) const
{
  /// Build a name for one detection element
  return Form("DE%04d",detElemId);
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::ManuName(Int_t manuId) const
{
  /// Build a name for one manu
  return Form("MANU%04d",manuId);
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::BusPatchName(Int_t busPatchId) const
{
  /// Build a name for one buspatch
  return Form("BUSPATCH%04d",busPatchId);
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::PCBName(Int_t pcbNumber) const
{
  /// Build a name for one pcb
  return Form("PCB%1d",pcbNumber);
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::ChamberPathName(Int_t chamberId) const
{
  /// Build a name for one chamber
  return Form("%s/%s",StationName(chamberId/2).Data(),ChamberName(chamberId).Data());
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::StationPathName(Int_t stationId) const
{
  /// Build a name for one station
  return StationName(stationId);
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::DEPathName(Int_t detElemId) const
{
  /// Build a name for one detection element
  
  Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
  
  return Form("%s/%s/%s",
              StationName(chamberId/2).Data(),
              ChamberName(chamberId).Data(),
              DEName(detElemId).Data());
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::ManuPathName(Int_t detElemId, Int_t manuId) const
{
  /// Build a name for one manu
  return Form("%s/%s",DEPathName(detElemId).Data(),ManuName(manuId).Data());
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::BusPatchPathName(Int_t busPatchId) const
{
  /// Build a name for one buspatch
  Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(busPatchId);
  
  return Form("%s/%s",DEPathName(detElemId).Data(),BusPatchName(busPatchId).Data());
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::PCBPathName(Int_t detElemId, Int_t pcbNumber) const
{
  /// Build a name for one pcb
  return Form("%s/%s",DEPathName(detElemId).Data(),PCBName(pcbNumber).Data());
}

//_____________________________________________________________________________
TString
AliMUONPainterHelper::FormatValue(const char* name, Double_t value) const
{
  /// Format a double value to be displayed
  /// FIXME: should insure we have the right number of significant digits here...
  
  return Form("%s = %e",name,value);
}

