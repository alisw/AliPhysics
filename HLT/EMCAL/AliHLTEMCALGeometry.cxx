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
#include "AliHLTEMCALConstant.h"
#include "assert.h"
#include "AliHLTCaloConstantsHandler.h"
#include "AliHLTEMCALSharedMemoryInterface.h" 
#include "TVector3.h"

using namespace EmcalHLTConst;

ClassImp(AliHLTEMCALGeometry);

AliHLTEMCALGeometry::AliHLTEMCALGeometry( TString det) : 
	AliHLTCaloGeometry (det), 
	fShmPtr(0),
	fGeo(0)
{
	//fGeo = new AliEMCALGeoUtils("EMCAL_COMPLETE","EMCAL");
	//fGeo = new AliEMCALGeometry("EMCAL_COMPLETE","EMCAL");
	//fGeo =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
//	TGeoManager::Import("/home/fedro/work/AliRoot/test/QA/geometry.root");
	fShmPtr = new AliHLTEMCALSharedMemoryInterface();
}


AliHLTEMCALGeometry::~AliHLTEMCALGeometry()
{

}
  
void 
AliHLTEMCALGeometry::GetGlobalCoordinates(AliHLTCaloRecPointDataStruct &/*recPoint*/, AliHLTCaloGlobalCoordinate &globalCoord)
{
	//TVector3 v1;
	Double_t glob[3] = {0, 0, 0};
 	//Int_t AbsId;

//	GetCellAbsId(recPoint.fModule, recPoint.fX, recPoint.fZ, AbsId);

	//fGeo->GetGlobalEMCAL(const AliEMCALRecPoint *rp, TVector3 &vglob);
//	fGeo->GetGlobal(AbsId, glob);

	globalCoord.fX = glob[0];
	globalCoord.fY = glob[1];
	globalCoord.fZ = glob[2];	
}
 
void 
AliHLTEMCALGeometry::GetCellAbsId(UInt_t /*module*/, UInt_t /*x*/, UInt_t /*z*/, Int_t& /*AbsId*/)  const
{
	
//	AbsId = fGeo->GetAbsCellIdFromCellIndexes(module, x, z); 
	
}


