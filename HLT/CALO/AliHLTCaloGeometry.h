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
   
   /** Constructor */
  AliHLTCaloGeometry (TString det);  //COMMENT
  
  /** Destructor */
  virtual ~AliHLTCaloGeometry(); //COMMENT
  
  /**
  * Get the global coordinates for a rec point. 
  * Pure virtual - must be implemented by child classes
  * @param recPoint is the rec point containing the local coordinates
  * @param globalCoord is a AliHLTCaloGlobalCoordinate to write to
  */
  virtual void GetGlobalCoordinates(AliHLTCaloRecPointDataStruct &recPoint, AliHLTCaloGlobalCoordinate &globalCoord ) = 0;  //COMMENT
  
  /**
  * Get the absolute ID from the relative position in the module
  * pure virtual - must be imlemented by child classes
  * @param module is the module id
  * @param x is the x position in the module
  * @param z is the z position in the module
  * @param AbsId is a the absolute id variable
  */
  virtual void GetCellAbsId(UInt_t module, UInt_t x, UInt_t z, Int_t& AbsId) const = 0; //COMMENT
  
  private:
     
    /** Default constructor, not implemented */
   AliHLTCaloGeometry();   //COMMENT


};

#endif

