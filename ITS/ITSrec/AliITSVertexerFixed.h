#ifndef ALIITSVERTEXERFIXED_H
#define ALIITSVERTEXERFIXED_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <AliITSVertexer.h>
class TString;

///////////////////////////////////////////////////////////////////////
//                                                                   //
// Create vertex in fixed position (e.g. TDI for injection tests)    //
// Origin: M. Masera, F. Prino (masera@to.infn.it, prino@to.infn.it) //
//                                                                   //
///////////////////////////////////////////////////////////////////////


class AliITSVertexerFixed : public AliITSVertexer {

 public:
  AliITSVertexerFixed();
  AliITSVertexerFixed(TString option);
  virtual AliESDVertex* FindVertexForCurrentEvent(TTree* /* itsClusterTree */ );
  void SetVtxPosition(Double_t pos[3]){
    for(Int_t k=0; k<3;k++) fVtxPos[k]=pos[k];
  }
  void SetVtxError(Double_t err[3]){
    for(Int_t k=0; k<3;k++) fVtxErr[k]=err[k];
  }

  virtual void PrintStatus() const;

 private:
  Double_t fVtxPos[3];    // vertex coordinates
  Double_t fVtxErr[3];    // vertex errors



ClassDef(AliITSVertexerFixed,0);
};

#endif
