#ifndef ALIHLTEMCALGEOMETRY_H
#define ALIHLTEMCALGEOMETRY_H
/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Federico Ronchetti for the ALICE HLT Project.*
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTCaloGeometry.h"
//#include "AliEMCALGeoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"

class TGeoManager;
class AliCDBManager;
class AliCDBPath;
class AliHLTEMCALRecPointDataStruct;

//class AliEMCALGeoUtils;
class AliEMCALGeometry;

class  AliHLTEMCALGeometry : public AliHLTCaloGeometry
{
 public:
	AliHLTEMCALGeometry();
	virtual ~AliHLTEMCALGeometry();
	void GetGlobalCoordinates(AliHLTCaloRecPointDataStruct &recPoint, AliHLTCaloGlobalCoordinate &globalCoord, Int_t iParticle );
	void GetCellAbsId(UInt_t module, UInt_t x, UInt_t z, Int_t& AbsId);
	virtual Int_t InitialiseGeometry();
	
	virtual void GetLocalCoordinatesFromAbsId(Int_t absId, Int_t& module, Int_t& x, Int_t& z);
	
	
protected:
	int GetGeometryFromCDB();

private:
	AliHLTEMCALGeometry(const AliHLTEMCALGeometry & );
	AliHLTEMCALGeometry & operator = (const AliHLTEMCALGeometry &);	

	// EMCal Geometry
	//AliEMCALGeoUtils *fGeo;
	
	AliEMCALGeometry *fGeo;
	AliEMCALRecoUtils *fReco;
};

#endif
