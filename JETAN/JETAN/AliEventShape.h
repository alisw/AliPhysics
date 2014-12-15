#ifndef ALIEVENTSHAPE_H
#define ALIEVENTSHAPE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Event shape utility class
// Circularity, Thrust, ... 
// Authors: Antonio Ortiz Velasquez <Antonio.Ortiz.Velasquez@cern.ch>
//          
//---------------------------------------------------------------------

#include <TObject.h>


class AliMCEvent;
class TParticle;
class TArrayD;

class AliEventShape : public TObject
{
  public:
    static TArrayD * GetThrustParamMC(AliMCEvent* mcEvent, Int_t  NSTUDYMIN=3, Double_t ptcutoff=1, Double_t etacutoff=1, Bool_t chom=kFALSE);
    static Double_t GetCircularityMC(AliMCEvent* mcEvent, Int_t  NSTUDYMIN=3, Double_t ptcutoff=1, Double_t etacutoff=1, Bool_t chom=kFALSE);
  private:
    AliEventShape(const AliEventShape&);
    AliEventShape& operator=(const AliEventShape&);
    ClassDef(AliEventShape, 0)
};

#endif

