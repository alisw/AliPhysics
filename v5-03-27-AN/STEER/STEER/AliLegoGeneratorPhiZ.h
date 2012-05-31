#ifndef ALILEGOGENERATORPHIZ_H
#define ALILEGOGENERATORPHIZ_H
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

class AliLegoGeneratorPhiZ : public AliLegoGenerator
{

 public:
  AliLegoGeneratorPhiZ(){;}
  virtual ~AliLegoGeneratorPhiZ(){;}
    virtual void    Generate();
    ClassDef(AliLegoGeneratorPhiZ,1) //Lego GeneratorPhiZ
};
#endif








