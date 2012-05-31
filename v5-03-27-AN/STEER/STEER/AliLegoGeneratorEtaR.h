#ifndef ALILEGOGENERATORETAR_H
#define ALILEGOGENERATORETAR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
//    Utility class to compute and draw Radiation Length Map                 //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliLegoGenerator.h"

class AliLegoGeneratorEtaR : public AliLegoGenerator
{

 public:
    AliLegoGeneratorEtaR() {}
    virtual ~AliLegoGeneratorEtaR() {}
    virtual void    Generate();
    ClassDef(AliLegoGeneratorEtaR, 1) // Lego GeneratorEtaR
};

#endif








