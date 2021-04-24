#ifndef EventStruct_cxx
#define EventStruct_cxx
#include "TTree.h"
#include "TObjString.h"
#include "Rtypes.h"


struct EventStruct {
  //V0 timing
  Float_t fTimeV0A;
  Float_t fTimeV0C;
  //Pileup
  Bool_t fIsPileup;
  Bool_t fIsPileupLowMult;
  //Physics selection
  Bool_t fIsPhysSel;
  //SPD vertex variables
  Int_t fNcontSPD;
  Float_t fVertexSPD_Z;
  Double_t fVtxSPDresZ;
  //Track vertex variables
  Int_t fNcontTrack;
  Float_t fVertexTrack_Z;
  Double_t fVtxTrackResZ;
  //Global vertex variables
  Int_t fNcontGlobal;
  Float_t fVertexGlobal_Z;
  Double_t fVtxGlobalResZ;
  //EventID
  ULong64_t fOrbit;
  unsigned int fBC;
  //Centrality(not used)
  Float_t fCentralityV0M;
  //T0 trigger
  UShort_t fTriggerT0;
  Bool_t fIsT0fired;
  //CTP trigger
  TObjString fTrigger;
  ULong64_t fTriggerMask;
  ULong64_t fTriggerMaskNext50;
  void clear();
  void makeTree(TTree *treeInput);
  void connectTree(TTree *treeInput);
};
#endif
