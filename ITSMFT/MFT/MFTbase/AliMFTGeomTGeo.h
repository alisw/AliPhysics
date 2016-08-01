#ifndef ALIITSUGEOMTGEO_H
#define ALIITSUGEOMTGEO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliMFTConstants.h"
#include "AliITSMFTGeomTGeo.h"

class AliMFTGeomTGeo : public AliITSMFTGeomTGeo {

public:

  AliMFTGeomTGeo();
  AliMFTGeomTGeo(const AliMFTGeomTGeo &src);
  ~AliMFTGeomTGeo();
  
  virtual TGeoPNEntry* GetPNEntry(Int_t index) const { return NULL; };
  virtual Bool_t GetOrigMatrix(Int_t index, TGeoHMatrix &m) const { 
    return kFALSE; };
  virtual TGeoHMatrix* ExtractMatrixSens(Int_t index) const { return NULL; };

  const Int_t GetNDisks() { return AliMFTConstants::kNDisks; }

private:

  AliMFTGeomTGeo& operator=(const AliMFTGeomTGeo &geom);

  ClassDef(AliMFTGeomTGeo, 1) // MFT geometry based on TGeo

};

#endif
