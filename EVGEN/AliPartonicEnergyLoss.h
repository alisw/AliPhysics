#ifndef ALIPARTONICENERGYLOSS_H
#define ALIPARTONICENERGYLOSS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Calculate parton energy loss
// in nucleus nucleus reactions
// Author: A.Morsch
//

#include <TObject.h>

class AliPartonicEnergyLoss : public TObject {
    
 public:
    AliPartonicEnergyLoss(){;}
    
    virtual ~AliPartonicEnergyLoss(){;}
    
    static void QuenchingWeight(Double_t r, Double_t x,
				Double_t& cont, Double_t& disc);
    static void RunTest();
 private:
    static void Init();
    
    ClassDef(AliPartonicEnergyLoss,1) // Library for partonic energy loss
};

#endif 



