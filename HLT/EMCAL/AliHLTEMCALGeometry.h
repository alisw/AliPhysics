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


#include  "AliHLTCaloGeometry.h"
#include "AliHLTEMCALSharedMemoryInterface.h" 
#include "AliEMCALGeoUtils.h"
#include "AliEMCALGeometry.h"
#include "TGeoManager.h"

class  AliHLTEMCALGeometry : public AliHLTCaloGeometry
{
 public:
  AliHLTEMCALGeometry( TString det );
  virtual ~AliHLTEMCALGeometry();
  void GetGlobalCoordinates(AliHLTCaloRecPointDataStruct &recPoint, AliHLTCaloGlobalCoordinate &globalCoord );
  void GetCellAbsId(UInt_t module, UInt_t x, UInt_t z, Int_t& AbsId)const;
 
private:
  	AliHLTEMCALGeometry();
	AliHLTEMCALSharedMemoryInterface* fShmPtr;  
	//AliEMCALGeometry *fGeo;
	AliEMCALGeoUtils *fGeo;
      AliHLTEMCALGeometry(const AliHLTEMCALGeometry & );
      AliHLTEMCALGeometry & operator = (const AliHLTEMCALGeometry &);
  
};

#endif
