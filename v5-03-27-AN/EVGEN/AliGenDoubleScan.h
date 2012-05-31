#ifndef ALIGENDOUBLESCAN_H
#define ALIGENDOUBLESCAN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// As AliGenScan,  generation of particles on a 3-dim grid
// but here double hits with a predefined distance are generated.
// Author: andreas.morsch@cern.ch

#include "AliGenScan.h"

class AliGenDoubleScan : public AliGenScan
{
 public:
    AliGenDoubleScan();
    AliGenDoubleScan(Int_t npart);
    virtual ~AliGenDoubleScan();
    virtual void Generate();
    virtual void SetDistance(Float_t d) {fDistance=d;}
 private:
    Float_t fDistance;           // Distance between particles
    ClassDef(AliGenDoubleScan,1) // Generation of particles (double hits) on a grid
};
#endif






