#ifndef ROOT_TMevSimConverter
#define ROOT_TMevSimConverter

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"

class TMevSimConverter : public TObject {

    enum {kMaxParticles = 35};

    Int_t  fNPDGCodes;                   // Number of PDG codes known by G3
    Int_t  fPDGCode [kMaxParticles];     // Translation table of PDG codes

    void DefineParticles();

  public:

    TMevSimConverter() {DefineParticles();}
    virtual ~TMevSimConverter() {}

    Int_t PDGFromId(Int_t gpid);
    Int_t IdFromPDG(Int_t pdg); 
    
    ClassDef(TMevSimConverter,1)
};

#endif




