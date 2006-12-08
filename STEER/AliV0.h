#ifndef ALIV0_H
#define ALIV0_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
//                          V0 Vertex Class
//          This class is part of the reconstruction
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------

#include "AliESDv0.h"

class AliV0: public AliESDv0 {
public:
  AliV0(): AliESDv0() {} 
  void Update(Float_t vertex[3]);            //update
  ClassDef(AliV0,1)      // V0 MI vertex
};


#endif
