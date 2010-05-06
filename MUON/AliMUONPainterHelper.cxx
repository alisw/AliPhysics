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
#include "AliMUONContourHandler.h"
#include "AliMUONContourMaker.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONManuContourMaker.h"
#include "AliMUONPainterEnv.h"
#include "AliMUONPainterRegistry.h"
#include "AliMUONPadStatusMaker.h"
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
#include <TMath.h>
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
fEnv(0x0),
fReal(0x0),
fExploded(0x0)
{
  /// ctor
  
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
  delete fReal;
  delete fExploded;
  fEnv->Save();
  fgInstance = 0;
}

//_____________________________________________________________________________
AliMUONContourHandler*
AliMUONPainterHelper::Exploded() const
{
  /// Create exploded contour handler
  if (!fExploded) fExploded = new AliMUONContourHandler(kTRUE);
  return fExploded;
}

//_____________________________________________________________________________
AliMUONContourHandler*
AliMUONPainterHelper::Real() const
{
  /// Create real contour handler
  if (!fReal) fReal = new AliMUONContourHandler(kFALSE);
  return fReal;
}

//_____________________________________________________________________________
AliMUONContour*
AliMUONPainterHelper::GetContour(const char* contourName, Bool_t explodedView) const
{
  /// Get a contour by name  
  if (explodedView) 
  {
    return Exploded()->GetContour(contourName);
  }
  else
  {
    if ( fReal ) 
    {
      return fReal->GetContour(contourName);
    }
  }
  return 0x0;
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
  const AliMpVSegmentation* vseg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
  if (vseg)
  {
    return vseg->MotifPosition(manuId);
  }
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

  AliCodeTimerAutoClass("",0);

  fgInstance = new AliMUONPainterHelper;
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
  
  TGeoHMatrix* matrix = static_cast<TGeoHMatrix*>(Exploded()->GetTransformations()->GetValue(detElemId));
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
  
  TGeoHMatrix* matrix = static_cast<TGeoHMatrix*>(Exploded()->GetTransformations()->GetValue(detElemId));
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
    if  ( TMath::AreEqualRel(max,min,1E-6) ) return gStyle->GetColorPalette(1);
    Double_t range = max - min;
    Double_t offset = value - min;
    rv = gStyle->GetColorPalette( 1 + int( offset*(gStyle->GetNumberOfColors()-2)/range - 0.5 ) );
  }
  return rv;
}

//_____________________________________________________________________________
AliMUONContour* 
AliMUONPainterHelper::MergeContours(const TObjArray& contours, const char* contourName, Bool_t explodedGeometry)
{
  /// Merge a set of contours (delegating to the contour maker)
  
  AliMUONContourMaker maker;
  
  AliMUONContour* contour = maker.MergeContour(contours,contourName);
  
  if (contour) 
  {
    RegisterContour(contour,explodedGeometry);
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
    if ( fExploded ) fExploded->Print();
    if ( fReal ) fReal->Print();
  }  
}

//_____________________________________________________________________________
void
AliMUONPainterHelper::RegisterContour(AliMUONContour* contour, Bool_t explodedView)
{
  /// contour is adopted by contourMaker
  AliCodeTimerAuto("",0)
  AliDebug(1,contour->GetName());
  AliMUONContourHandler* ch = fReal;
  if ( explodedView ) 
  {
    ch = Exploded();
  }
  if (!ch)
  {
    AliError(Form("ContourHandler for %s view is not created yet !",explodedView ? "EXPLODED" : "REAL"));
  }
  else
  {
    if ( ch->GetContour(contour->GetName()) )
    {
      AliError(Form("Contour with name %s is already there",contour->GetName()));
      return;
    }
    ch->Adopt(contour);
  }
}

//_____________________________________________________________________________
TString 
AliMUONPainterHelper::ChamberName(Int_t chamberId) const
{
  /// Build a name for one chamber
  return Form("Chamber%1d",chamberId+1);
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
  
  TString sname(name);
  
  sname.ToUpper();
  if (sname.Contains("BIT"))
  {
    Int_t i = (Int_t)(value);
    TString rv = Form("%s = 0x%x",name,i);
    cout << rv << ":" << AliMUONPadStatusMaker::AsString(i) << endl;
    return rv;
  }
  else
  {
    return Form("%s = %e",name,value);
  }
}

//_____________________________________________________________________________
TObjArray*
AliMUONPainterHelper::GetAllContoursAsArray(Bool_t explodedView) const
{
  /// Get the contours in a specially arranged array (orderer by hierarchy level)
  
  if ( explodedView ) 
  {
    return Exploded()->AllContourArray();
  }
  else
  {
    return Real()->AllContourArray();
  }
}


