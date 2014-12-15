#ifndef ALIITSVERTEXERCOSMICS_H
#define ALIITSVERTEXERCOSMICS_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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
  virtual AliESDVertex* FindVertexForCurrentEvent(TTree *itsClusterTree);
  virtual void PrintStatus() const;
  void SetFirstLastModules(Int_t ilayer=0,Int_t m1=0,Int_t m2=79) 
    {fFirst[ilayer] = m1; fLast[ilayer] = m2;}
  void SetMaxDistOnOuterLayer(Double_t max=1.0) {fMaxDistOnOuterLayer=max;}
  Double_t GetMaxDistOnOuterLayer() const {return fMaxDistOnOuterLayer;}
  void SetMaxVtxRadius(Int_t ilayer=0,Double_t maxr=3.5) {fMaxVtxRadius[ilayer]=maxr;}
  Double_t GetMaVtxRadius(Int_t ilayer=0) const {return fMaxVtxRadius[ilayer];}
  void SetMinDist2Vtxs(Double_t mind=0.1) {fMinDist2Vtxs=mind;}
  Double_t GetMinDist2Vtxs() const {return fMinDist2Vtxs;}

 private:

  Int_t fFirst[6];          // first module of each layer
  Int_t fLast[6];           // last module of each layer
  Double_t fMaxDistOnOuterLayer;  // max dca between tracklet & outer layer cls
  Double_t fMaxVtxRadius[6];    // maximum radial pos of vertex
  Double_t fMinDist2Vtxs;    // minimum distance between two vertices

  ClassDef(AliITSVertexerCosmics,4); // vertexer for cosmics
};

#endif
