#ifndef EventStructT0_cxx
#define EventStructT0_cxx
#include "TTree.h"
#include "TObjString.h"
#include "Rtypes.h"

struct EventStructT0 {
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
  void clear() {
    fTimeV0A=-9999;
    fTimeV0C=-9999;
    fIsPileup=kFALSE;
    fIsPileupLowMult=kFALSE;
    fIsPhysSel=kFALSE;
    fNcontSPD=-1;
    fVertexSPD_Z=-9999;
    fVtxSPDresZ=-9999;
    fNcontTrack=-1;
    fVertexTrack_Z=-9999;
    fVtxTrackResZ=-9999;
    fNcontGlobal=-1;
    fVertexGlobal_Z=-9999;
    fVtxGlobalResZ=-9999;
    fOrbit=0;
    fBC=0;
    fCentralityV0M=-300;
    fTriggerT0=0;
    fIsT0fired=kFALSE;
    fTrigger="";
    fTriggerMask=0;
    fTriggerMaskNext50=0;
  }
  void makeTree(TTree *treeInput) {
    //V0 timing
    treeInput->Branch("fTimeV0A", &fTimeV0A);
    treeInput->Branch("fTimeV0C", &fTimeV0C);
    //Pileup
    treeInput->Branch("fIsPileup", &fIsPileup);
    treeInput->Branch("fIsPileupLowMult", &fIsPileupLowMult);
    //Physics selection
    treeInput->Branch("fIsPhysSel",&fIsPhysSel);
    //Vertex SPD
    treeInput->Branch("fNcontSPD",&fNcontSPD);
    treeInput->Branch("fVertexSPD_Z",&fVertexSPD_Z);
    treeInput->Branch("fVtxSPDresZ",&fVtxSPDresZ);
    //Vertex Track
    treeInput->Branch("fNcontTrack",&fNcontTrack);
    treeInput->Branch("fVertexTrack_Z",&fVertexTrack_Z);
    treeInput->Branch("fVtxTrackResZ",&fVtxTrackResZ);
    //Vertex Global
    treeInput->Branch("fNcontGlobal",&fNcontGlobal);
    treeInput->Branch("fVertexGlobal_Z",&fVertexGlobal_Z);
    treeInput->Branch("fVtxGlobalResZ",&fVtxGlobalResZ);
    //Event ID
    treeInput->Branch("fOrbit", &fOrbit);
    treeInput->Branch("fBC", &fBC);
    //Centrality (not used)
  //	treeInput->Branch("fCentralityV0M", &fCentralityV0M);
    //T0 trigger
    treeInput->Branch("fTriggerT0", &fTriggerT0);
    treeInput->Branch("fIsT0fired", &fIsT0fired);
    //CTP trigger
    treeInput->Branch("fTriggerMask", &fTriggerMask);
    treeInput->Branch("fTriggerMaskNext50", &fTriggerMaskNext50);
  }
  void connectTree(TTree *treeInput) {
    //V0 timing
    treeInput->SetBranchAddress("fTimeV0A", &fTimeV0A);
    treeInput->SetBranchAddress("fTimeV0C", &fTimeV0C);
    //Pileup
    treeInput->SetBranchAddress("fIsPileup", &fIsPileup);
    treeInput->SetBranchAddress("fIsPileupLowMult", &fIsPileupLowMult);
    //Physics selection
    treeInput->SetBranchAddress("fIsPhysSel",&fIsPhysSel);
    //Vertex SPD
    treeInput->SetBranchAddress("fNcontSPD",&fNcontSPD);
    treeInput->SetBranchAddress("fVertexSPD_Z",&fVertexSPD_Z);
    treeInput->SetBranchAddress("fVtxSPDresZ",&fVtxSPDresZ);
    //Vertex Track
    treeInput->SetBranchAddress("fNcontTrack",&fNcontTrack);
    treeInput->SetBranchAddress("fVertexTrack_Z",&fVertexTrack_Z);
    treeInput->SetBranchAddress("fVtxTrackResZ",&fVtxTrackResZ);
    //Vertex Global
    treeInput->SetBranchAddress("fNcontGlobal",&fNcontGlobal);
    treeInput->SetBranchAddress("fVertexGlobal_Z",&fVertexGlobal_Z);
    treeInput->SetBranchAddress("fVtxGlobalResZ",&fVtxGlobalResZ);
    //Event ID
    treeInput->SetBranchAddress("fOrbit", &fOrbit);
    treeInput->SetBranchAddress("fBC", &fBC);
    //Centrality (not used)
  //	treeInput->SetBranchAddress("fCentralityV0M", &fCentralityV0M);
    //T0 trigger
    treeInput->SetBranchAddress("fTriggerT0", &fTriggerT0);
    treeInput->SetBranchAddress("fIsT0fired", &fIsT0fired);
    //CTP trigger
    treeInput->SetBranchAddress("fTriggerMask", &fTriggerMask);
    treeInput->SetBranchAddress("fTriggerMaskNext50", &fTriggerMaskNext50);
  }

};
#endif
