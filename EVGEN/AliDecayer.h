#ifndef ALIDECAYER_H
#define ALIDECAYER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include "GenTypeDefs.h"
#include <TObject.h>
class TClonesArray;
class TLorentzVector;

class AliDecayer :
public TObject
{
 public:
//
    virtual void    Init()                                     =0;
    virtual void    Decay(Int_t idpart, TLorentzVector* p)     =0;
    virtual Int_t   ImportParticles(TClonesArray *particles)   =0;
    virtual void    SetForceDecay(Decay_t type)                =0;
    virtual void    ForceDecay()                               =0;
    virtual Float_t GetPartialBranchingRatio(Int_t ipart)      =0;
    ClassDef(AliDecayer,1) // Alice Decayer Base Class
};
#endif







