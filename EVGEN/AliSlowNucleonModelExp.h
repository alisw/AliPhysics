#ifndef ALISLOWNUCLEONMODELEXP_H
#define ALISLOWNUCLEONMODELEXP_H
/* Copyright(c) 198-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Experimental data inspired Gray Particle Model for p-Pb collisions
// Fluctuations are calculated from a binomial distribution.
// Author: A.Morsch
//

#include "AliSlowNucleonModel.h"

class AliCollisionGeometry;

class AliSlowNucleonModelExp : public AliSlowNucleonModel
{
 public:
    AliSlowNucleonModelExp();
    virtual ~AliSlowNucleonModelExp(){;}
    virtual void GetNumberOfSlowNucleons(AliCollisionGeometry* geo,
					 Int_t& ngp, Int_t& ngn, Int_t& nbp, Int_t& nbn) const;
    virtual void SetParameters(Float_t alpha1, Float_t alpha2);
    
 protected:
    Float_t  fP;          // Number of protons  in the target 
    Float_t  fN;          // Number of neutrons in the target
    Float_t  fAlphaGray;  // Proportionality between gray   particles and number of collisions
    Float_t  fAlphaBlack; // Proportionality between black  particles and number of collisions
    
  ClassDef(AliSlowNucleonModelExp,1) // Gray Particle Model (Experiment inspired)
};
#endif






