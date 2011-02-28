#ifndef ALIABSO_H
#define ALIABSO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: ABSO          //
////////////////////////////////////////////////
 
#include "AliModule.h"
 
 
class AliABSO : public AliModule {
    
 public:
    AliABSO();
    AliABSO(const char *name, const char *title);
    virtual      ~AliABSO() {}
    virtual void    CreateGeometry();
    virtual void    CreateMaterials();
    virtual void    Init();
    virtual Int_t   IsVersion() const {return 0;}
    virtual Int_t   GetMatId(Int_t imat) const;
    virtual Int_t   NumberOfLayers(Int_t i) const {return fNLayers[i];}
    virtual Float_t ZPositionOfLayer(Int_t i, Int_t il) const 
      {return fZLayers[i][il];}    
    virtual Int_t   MaterialOfLayer (Int_t i, Int_t il) const 
      {return fMLayers[i][il];}    	  
 protected:
    Int_t   fNLayers[2];        // Number of Material Layers in the tracking Region
    Float_t fZLayers[2][15];     // z-position of layers
    Int_t   fMLayers[2][15];     // Material type of layers
  ClassDef(AliABSO,1)  // Muon Absorber Class
};

#endif
