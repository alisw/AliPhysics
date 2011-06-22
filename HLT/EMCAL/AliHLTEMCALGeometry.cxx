/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: federico ronchetti         for the ALICE HLT Project.*
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTEMCALGeometry.h"
#include "AliHLTEMCALConstants.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBPath.h"


ClassImp(AliHLTEMCALGeometry);
TGeoManager *gGeoManager = 0;

AliHLTEMCALGeometry::AliHLTEMCALGeometry() :
	AliHLTCaloGeometry ("EMCAL"),
	fGeo(0),fReco(0)
{
  GetGeometryFromCDB();
}

Int_t AliHLTEMCALGeometry::InitialiseGeometry()
{
   
   return GetGeometryFromCDB();
}


AliHLTEMCALGeometry::~AliHLTEMCALGeometry()
{

}
  
void 
AliHLTEMCALGeometry::GetGlobalCoordinates(AliHLTCaloRecPointDataStruct &recPoint, AliHLTCaloGlobalCoordinate &globalCoord)
{

  Float_t fDepth;
  Float_t *fRot = fReco->GetMisalRotShiftArray();
  Float_t *fTrans = fReco->GetMisalTransShiftArray();
  Float_t glob[] = {0.,0.,0.};

  //assume photo for the moment
  fDepth = fReco->GetDepth(recPoint.fAmp,AliEMCALRecoUtils::kPhoton,recPoint.fModule);
  
  fGeo->RecalculateTowerPosition(recPoint.fX, recPoint.fZ,recPoint.fModule, fDepth, fTrans, fRot, glob);
  
  globalCoord.fX = glob[0];
  globalCoord.fY = glob[1];
  globalCoord.fZ = glob[2];
  

}
 
void 
AliHLTEMCALGeometry::GetCellAbsId(UInt_t module, UInt_t x, UInt_t z, Int_t& AbsId)
{

  if(!fGeo)
    {
      Logging(kHLTLogError, "HLT", "EMCAL", "AliHLTEMCALGeometry::GetCellAbsId: no geometry initialised");
      return;

    }
	AbsId = fGeo->GetAbsCellIdFromCellIndexes(module, (Int_t) x, (Int_t) z);


	
}


int
AliHLTEMCALGeometry::GetGeometryFromCDB()
{

  // AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  AliCDBPath path("GRP","Geometry","Data");
  if(path.GetPath())
    {
      //      HLTInfo("configure from entry %s", path.GetPath());
      AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
      if (pEntry)
	{
	  if(!fGeo)
	    {
	      delete fGeo;
	      fGeo = 0;
	    }

	  gGeoManager = (TGeoManager*) pEntry->GetObject();

	  if(gGeoManager)
	    {
	      fGeo = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
	      //fGeo = new AliEMCALGeoUtils("EMCAL_COMPLETE","EMCAL");
	      fReco = new AliEMCALRecoUtils;
	      // FIXME
	      // need to be parametrized
	      // misalignment corrections to be put into OCDB

	      fReco->SetMisalTransShift(0,1.08); 
	      fReco->SetMisalTransShift(1,8.35); 
	      fReco->SetMisalTransShift(2,1.12); //sector 0
	      fReco->SetMisalRotShift(3,-8.05); 
	      fReco->SetMisalRotShift(4,8.05); 
	      fReco->SetMisalTransShift(3,-0.42); 
	      fReco->SetMisalTransShift(5,1.55);//sector 1
	      
	    }

	}
      else
	{
    	  //HLTError("can not fetch object \"%s\" from OCDB", path);
    	  Logging(kHLTLogError, "HLT", "EMCAL", "can not fetch object from OCDB");

	}
    }
  return 0;
}
