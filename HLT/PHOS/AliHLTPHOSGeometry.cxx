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
 GetGeometryFromCDB();
}

AliHLTPHOSGeometry::~AliHLTPHOSGeometry()
{
// See header file for class documentation
}

void AliHLTPHOSGeometry::GetGlobalCoordinates ( AliHLTCaloRecPointDataStruct& recPoint, AliHLTCaloGlobalCoordinate& globalCoord )
{
   // See header file for class documentation
   Float_t x = recPoint.fX;
   Float_t z = recPoint.fZ;

   ConvertRecPointCoordinates(x, z);

   TVector3 coord;
   
   fGeoUtils->Local2Global(recPoint.fModule+1, x, z, coord);
   
   globalCoord.fX = coord[0];
   globalCoord.fZ = coord[1];
   globalCoord.fY = coord[2];
}

void AliHLTPHOSGeometry::ConvertRecPointCoordinates(Float_t &x, Float_t &z) const
{
   x = (x - fCaloConstants->GetNXCOLUMNSMOD())*fCaloConstants->GetCELLSTEP();
   z = (z - fCaloConstants->GetNZROWSMOD())*fCaloConstants->GetCELLSTEP();
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
	  if(!fGeoUtils) 
	    {
	      delete fGeoUtils;
	      fGeoUtils = 0;
	    }

	  gGeoManager = (TGeoManager*) pEntry->GetObject();
//	  HLTError("gGeoManager = 0x%x", gGeoManager);
	  if(gGeoManager)
	    {
	      fGeoUtils = new AliPHOSGeoUtils("PHOS", "noCPV");
	    }
	}
      else
	{
//	    HLTError("can not fetch object \"%s\" from OCDB", path);
	}
    }
}



