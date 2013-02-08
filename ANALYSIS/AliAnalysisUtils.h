#ifndef ALIANALYSISUTILS_H
#define ALIANALYSISUTILS_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// Class with functions useful for different analyses                       //
// - vertex selection                                                       //
//    * 2013 pA default cuts                                                //
// - identification of the fist event of the chunk                          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliVEvent;

class AliAnalysisUtils : public TObject {

 public:

  AliAnalysisUtils();
  virtual ~AliAnalysisUtils(){};

  Bool_t IsVertexSelected2013pA(AliVEvent *event);
  Bool_t IsFirstEventInChunk(AliVEvent *event);

  void SetMinVtxContr(Int_t contr=1) {fMinVtxContr=contr;}
  void SetMaxVtxZ(Float_t z=1e6) {fMaxVtxZ=z;}
  void SetCutOnZVertexSPD(Bool_t iscut=true) { fCutOnZVertexSPD = iscut; }

 private:

  Bool_t fisAOD; // flag for AOD:1 or ESD:0

  Int_t fMinVtxContr; // minimum vertex contributors
  Float_t fMaxVtxZ;   // maximum |z| of primary vertex

  Bool_t fCutOnZVertexSPD; // 0: no cut, 1: |zvtx-SPD - zvtx-TPC|<0.5cm

  ClassDef(AliAnalysisUtils,0) // base helper class
};
#endif
 
