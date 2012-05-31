#ifndef ALIKINK_H
#define ALIKINK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          Kink Vertex Class
//          This class is part of the reconstruction
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------

#include "AliESDkink.h"

class AliKink : public AliESDkink {
public:
  AliKink(){;}             //constructor
  void Update();            //update
  ClassDef(AliKink,1)      //kink vertex
};

#endif


