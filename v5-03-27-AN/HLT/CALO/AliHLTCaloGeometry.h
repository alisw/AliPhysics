//-*- Mode: C++ -*-
#ifndef ALIHLTCALOGEOMETRY_H
#define ALIHLTCALOGEOMETRY_H
/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Federico Ronchetti 

#include "Rtypes.h"
#include "AliHLTCaloConstants.h"
#include "AliHLTLogging.h"
#include "AliHLTCaloConstantsHandler.h"
#include "AliHLTCaloCoordinate.h"
#include "AliHLTCaloGlobalCoordinate.h"
#include "AliHLTCaloRecPointDataStruct.h"

struct AliHLTCaloCoordinate;
struct AliHLTCaloGlobalCoordinate;

class AliHLTCaloGeometry : public AliHLTCaloConstantsHandler, public AliHLTLogging
{

 public:
  AliHLTCaloGeometry (TString det);
  virtual ~AliHLTCaloGeometry();
  // Particle: 0=photon, 1=electron, 2=hadron
  virtual void GetGlobalCoordinates(AliHLTCaloRecPointDataStruct &recPoint, AliHLTCaloGlobalCoordinate &globalCoord, Int_t iParticle ) = 0;

  virtual void GetCellAbsId(UInt_t module, UInt_t x, UInt_t z, Int_t& AbsId) = 0; //COMMENT
  
  virtual void GetLocalCoordinatesFromAbsId(Int_t AbsId, Int_t &module, Int_t &x, Int_t &z) = 0; //COMMENT
  
  virtual Int_t InitialiseGeometry() = 0;

protected:
  Bool_t fIsInitialised;

 private:

  /** Default constructor, not implemented */
   AliHLTCaloGeometry();   //COMMENT

   ClassDef(AliHLTCaloGeometry, 0);
   
};

#endif

