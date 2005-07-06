#ifndef ALICASCADEVERTEX_H
#define ALICASCADEVERTEX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                          Cascade Vertex Class                          //
//                                                                        //
// Origin: Christian Kuhn, IReS, Strasbourg, christian.kuhn@ires.in2p3.fr //
////////////////////////////////////////////////////////////////////////////

#include "AliESDcascade.h"

class AliITStrackV2;
class AliV0vertex;

#define kXiMinus       3312
#define kXiPlusBar    -3312
#define kOmegaMinus    3334
#define kOmegaPlusBar -3334

class AliCascadeVertex : public AliESDcascade {
public:
  AliCascadeVertex():AliESDcascade(){;}
  AliCascadeVertex(const AliV0vertex &vtx, const AliITStrackV2 &trk);

  ClassDef(AliCascadeVertex,1)   // reconstructed cascade vertex
};

#endif


