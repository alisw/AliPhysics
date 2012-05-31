#ifndef ALILEGOGENERATORETA_H
#define ALILEGOGENERATORETA_H
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

class AliLegoGeneratorEta : public AliLegoGenerator
{

 public:
    AliLegoGeneratorEta() {}
    virtual ~AliLegoGeneratorEta() {}
    virtual void    Generate();
    ClassDef(AliLegoGeneratorEta,1) //Lego GeneratorEta
};

#endif








