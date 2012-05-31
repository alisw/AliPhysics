#ifndef ALISLOWNUCLEONMODEL_H
#define ALISLOWNUCLEONMODEL_H
/* Copyright(c) 198-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"
class AliCollisionGeometry;

class AliSlowNucleonModel : public TObject
{
public:
    AliSlowNucleonModel() {;}
    virtual ~AliSlowNucleonModel(){;}
    virtual void GetNumberOfSlowNucleons(AliCollisionGeometry* /*geo*/,
					 Int_t& /*ngp*/, Int_t& /*ngn*/,
					 Int_t& /*nbp*/, Int_t& /*nbn*/) const {;}
    
 protected:
  ClassDef(AliSlowNucleonModel,1) // Gray Particle Model
};
#endif






