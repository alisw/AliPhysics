#ifndef ALIITSVERTEXERFAST_H
#define ALIITSVERTEXERFAST_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include <AliITSVertexer.h>

/////////////////////////////////////////////////////////////////////
//                                                                 //
// Fast vertexer - True (i.e. generated) vertex coordinates        //
// are smeared with gaussians of given width                       //
//                                                                 //
/////////////////////////////////////////////////////////////////////


class AliITSVertexerFast : public AliITSVertexer {

 public:
  AliITSVertexerFast();
  AliITSVertexerFast(Double_t *smear);
  virtual ~AliITSVertexerFast(); 
  virtual AliESDVertex* FindVertexForCurrentEvent(Int_t evnumb);
  virtual void FindVertices();
  virtual void PrintStatus() const;

 protected:

  // copy constructor (NO copy allowed: the constructor is protected
  // to avoid misuse)
  AliITSVertexerFast(const AliITSVertexerFast& vtxr);
  // assignment operator (NO assignment allowed)
  AliITSVertexerFast& operator=(const AliITSVertexerFast& /* vtxr */);

  Double_t *fSmear;         // rms of gaussians used for smearing


ClassDef(AliITSVertexerFast,1);
};

#endif
