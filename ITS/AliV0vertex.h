#ifndef ALIV0VERTEX_H
#define ALIV0VERTEX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          V0 Vertex Class
//
//    Origin: Iouri Belikov, IReS, Strasbourg, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include <AliESDv0.h>

class AliITStrackV2;

class AliV0vertex : public AliESDv0 {
public:
  AliV0vertex() : AliESDv0() {;}
  AliV0vertex(const AliESDv0 &v) : AliESDv0(v) {;}
  AliV0vertex(const AliITStrackV2 &neg, const AliITStrackV2 &pos);

  ClassDef(AliV0vertex,1) // reconstructed V0 vertex
};

#endif


