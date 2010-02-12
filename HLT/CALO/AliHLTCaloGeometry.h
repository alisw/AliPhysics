//-*- Mode: C++ -*-
#ifndef ALIHLTCALOGEOMETRY_H
#define ALIHLTCALOGEOMETRY_H
/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Federico Ronchetti 

#include "Rtypes.h"
#include "AliHLTCaloConstants.h"
#include "AliHLTCaloConstantsHandler.h"
#include "AliHLTCaloCoordinate.h"
#include "AliHLTCaloGlobalCoordinate.h"
#include "AliHLTCaloRecPointDataStruct.h"

class AliHLTCaloCoordinate;
class AliHLTCaloGlobalCoordinate;

class AliHLTCaloGeometry : public AliHLTCaloConstantsHandler
{

 public:
  AliHLTCaloGeometry (TString det);
  virtual ~AliHLTCaloGeometry();
  virtual void GetGlobalCoordinates(AliHLTCaloRecPointDataStruct &recPoint, AliHLTCaloGlobalCoordinate &globalCoord ) = 0;
  
  private:
     
    
   AliHLTCaloGeometry();


};

#endif

