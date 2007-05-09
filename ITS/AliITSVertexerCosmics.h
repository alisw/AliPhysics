#ifndef ALIITSVERTEXERCOSMICS_H
#define ALIITSVERTEXERCOSMICS_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#include "AliESDVertex.h"
#include "AliITSVertexer.h"

//-----------------------------------------------------------------------
//                                                                 
// Class for constructing a fake primary vertex for cosmics events 
//                                                                 
// Origin: A.Dainese andrea.dainese@lnl.infn.it                    
//-----------------------------------------------------------------------

class AliESDVertex;

class AliITSVertexerCosmics : public AliITSVertexer {

 public:

  AliITSVertexerCosmics();
  virtual ~AliITSVertexerCosmics() {}
  AliESDVertex* FindVertexForCurrentEvent(Int_t evnumb);
  void FindVertices();
  void PrintStatus() const;
  void SetSPD1Modules(Int_t m1=0,Int_t m2=79) {fFirstSPD1 = m1; fLastSPD1 = m2;}
  void SetSPD2Modules(Int_t m1=80, Int_t m2=239) {fFirstSPD2 = m1; fLastSPD2 = m2;}
  void SetMaxDistOnSPD2(Double_t max=0.1) {fMaxDistOnSPD2=max;}
  Double_t GetMaxDistOnSPD2() const {return fMaxDistOnSPD2;}
  void SetMaxVtxRadius(Double_t maxr=3.5) {fMaxVtxRadius=maxr;}
  Double_t GetMaVtxRadius() const {return fMaxVtxRadius;}
  void SetMinDist2Vtxs(Double_t mind=0.1) {fMinDist2Vtxs=mind;}
  Double_t GetMinDist2Vtxs() const {return fMinDist2Vtxs;}

 private:

  Int_t fFirstSPD1;          // first module of the first pixel layer used
  Int_t fLastSPD1;           // last module of the first pixel layer used
  Int_t fFirstSPD2;          // first module of the second pixel layer used
  Int_t fLastSPD2;           // last module of the second pixel layer used
  Double_t fMaxDistOnSPD2;   // max dca between SPD1 tracklet and SPD2 cls
  Double_t fMaxVtxRadius;    // maximum radial pos of vertex
  Double_t fMinDist2Vtxs;    // minimum distance between two vertices

  ClassDef(AliITSVertexerCosmics,1);
};

#endif
