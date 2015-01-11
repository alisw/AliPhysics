/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Oystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************
*/

#include "AliHLTPHOSGeometry.h"
#include "AliPHOSGeoUtils.h"
#include "TGeoManager.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "TVector3.h"

AliHLTPHOSGeometry::AliHLTPHOSGeometry() :
AliHLTCaloGeometry("PHOS"),
fGeoUtils(0)
{
 // See header file for class documentation
 
}

AliHLTPHOSGeometry::~AliHLTPHOSGeometry()
{
// See header file for class documentation
}
// FR: PHOS doesn't use iParticle for now
void AliHLTPHOSGeometry::GetGlobalCoordinates ( AliHLTCaloRecPointDataStruct& recPoint, AliHLTCaloGlobalCoordinate& globalCoord, Int_t /* iParticle */ )
{
   // See header file for class documentation
   if(!fIsInitialised) { InitialiseGeometry(); }
   if(!fGeoUtils) 
   {
      Logging(kHLTLogError, "HLT", "PHOS", "AliHLTPHOSGeometry::GetGlobalCoordinates: no geometry initialised");
      return;
   }

   Float_t x = recPoint.fX;
   Float_t z = recPoint.fZ;
   

   ConvertRecPointCoordinates(x, z);
   
   TVector3 coord;
   fGeoUtils->Local2Global(fCaloConstants->GetNMODULES() - recPoint.fModule, x, z, coord);
   

   globalCoord.fX = coord[0];
   globalCoord.fY = coord[1];
   globalCoord.fZ = coord[2];
   
}

void AliHLTPHOSGeometry::ConvertRecPointCoordinates(Float_t &x, Float_t &z) const
{
   // See header file for class documentation
   x = (x - (float)(fCaloConstants->GetNXCOLUMNSMOD())/2)*fCaloConstants->GetCELLSTEP();
   z = (z - ((float)(fCaloConstants->GetNZROWSMOD()))/2)*fCaloConstants->GetCELLSTEP();
}

int AliHLTPHOSGeometry::GetGeometryFromCDB()
{
   // See header file for documentation

   AliCDBPath path("GRP","Geometry","Data");
   if(path.GetPath())
    {
      //      HLTInfo("configure from entry %s", path.GetPath());
      AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
      if (pEntry) 
	{
	  if(fGeoUtils) 
	    {
	      delete fGeoUtils;
	      fGeoUtils = 0;
	    }
	  
	  if(!gGeoManager) gGeoManager = (TGeoManager*) pEntry->GetObject();
	  
	  if(gGeoManager)
	    {
	      fGeoUtils = new AliPHOSGeoUtils("PHOS", "noCPV");
	      if(fGeoUtils) fIsInitialised = kTRUE;
	    }
	    else
	    {
	       HLTError("can not get gGeoManager from OCDB");
	    }
	}
      else
	{
	    HLTError("can not fetch object \"%s\" from OCDB", path.GetPath().Data());
	}
    }
    return 0;
}


void AliHLTPHOSGeometry::GetCellAbsId ( UInt_t module, UInt_t x, UInt_t z, Int_t& AbsId ) 
  {
      // See header file for class documentation
      if(!fGeoUtils)
      {
	 Logging(kHLTLogError, "HLT", "PHOS", "AliHLTPHOSGeometry::GetCellAbsId: no geometry initialised");
	 return;
      }
      fGeoUtils->RelPosToAbsId(module, x, z, AbsId);
  }

void AliHLTPHOSGeometry::GetLocalCoordinatesFromAbsId(Int_t absId, Int_t& module, Int_t& x, Int_t& z)
{
  Int_t rel[4];
  fGeoUtils->AbsToRelNumbering(absId, rel);
  module = rel[0]-1;
  z = rel[2];
  x = rel[3];
}
