#ifndef ALIFASTPARTICLE_H
#define ALIFASTPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TParticle.h>
class AliFastParticle : public TParticle {
 public:
    AliFastParticle(){;}
    virtual ~AliFastParticle(){;}
 protected:
    ClassDef(AliFastParticle,1) // Base class for fast particle
};

#endif 



