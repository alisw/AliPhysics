#ifndef ALIGRAYPARTICLEMODEL_H
#define ALIGRAYPARTICLEMODEL_H
/* Copyright(c) 198-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"
class AliCollisionGeometry;

class AliGrayParticleModel : public TObject
{
public:
    AliGrayParticleModel() {;}
    virtual ~AliGrayParticleModel(){;}
    virtual void GetNumberOfGrayNucleons(AliCollisionGeometry* geo, Int_t& np, Int_t& nn) {;}
    
 protected:
  ClassDef(AliGrayParticleModel,1) // Gray Particle Model
};
#endif






