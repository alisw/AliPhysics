#ifndef AliGenDoubleScan_H
#define AliGenDoubleScan_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////
#include "AliGenScan.h"
#include "TNamed.h"
#include "TF1.h"
#include "TArrayF.h"
#include "TTree.h"


class AliGenDoubleScan : public AliGenScan
{
 private:
    Float_t fDistance;
 public:
    AliGenDoubleScan();
    AliGenDoubleScan(Int_t npart);
    virtual ~AliGenDoubleScan();
    // generate event
    virtual void Generate();
    virtual void SetDistance(Float_t d) {fDistance=d;}
    ClassDef(AliGenDoubleScan,1) // Generation of particles on a grid 
};
#endif






