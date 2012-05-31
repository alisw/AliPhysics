#ifndef ALIGENHIJINGPARABA_H
#define ALIGENHIJINGPARABA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Parameterisation of pi, K, n and p eta and pt distributions
// eta: according to HIJING (shadowing + quenching)              
// pT : according to CDF measurement at 1.8 TeV                  
// Author: andreas.morsch@cern.ch

#include "AliGenHIJINGpara.h"
class TF1;

class AliGenHIJINGparaBa : public AliGenHIJINGpara
{
 public:
    AliGenHIJINGparaBa();
    AliGenHIJINGparaBa(Int_t npart);
    virtual ~AliGenHIJINGparaBa();
    virtual void Generate();
    virtual void Init();
 protected:
    TF1* fPtba;          //! Parametrised pt distribution for baryons
    TF1* fETAba;         //! Parametrised eta distribution for baryons

 private:
    AliGenHIJINGparaBa(const AliGenHIJINGparaBa &HIJINGpara);
    AliGenHIJINGparaBa& operator = (const AliGenHIJINGparaBa &para);

    ClassDef(AliGenHIJINGparaBa,1) // Hijing parametrisation generator with baryons
};
#endif










