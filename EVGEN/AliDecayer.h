#ifndef ALIDECAYER_H
#define ALIDECAYER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Abstract base class for particle decays.
// Clients are the transport code and the primary particle generators
// Author: andreas.morsch@cern.ch

#include <TObject.h>
class TClonesArray;
class TLorentzVector;

typedef enum
{ kSemiElectronic, kDiElectron, kSemiMuonic, kDiMuon,
  kBJpsiDiMuon, kBJpsiDiElectron, 
  kBPsiPrimeDiMuon, kBPsiPrimeDiElectron, kPiToMu, kKaToMu, kNoDecay, kHadronicD, kAll}
Decay_t;

class AliDecayer :
public TObject
{
 public:
//
    virtual ~AliDecayer(){;}
    virtual void    Init()                                     =0;
    virtual void    Decay(Int_t idpart, TLorentzVector* p)     =0;
    virtual Int_t   ImportParticles(TClonesArray *particles)   =0;
    virtual void    SetForceDecay(Decay_t type)                =0;
    virtual void    ForceDecay()                               =0;
    virtual Float_t GetPartialBranchingRatio(Int_t ipart)      =0;
    ClassDef(AliDecayer,1) // Alice Decayer Base Class
};
#endif







