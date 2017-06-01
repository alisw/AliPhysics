#ifndef AliAnalysisTaskTrigHMTF_h
#define AliAnalysisTaskTrigHMTF_h

class TTree;
#include "AliAnalysisTaskSE.h"
#include "TBits.h"
#include "TObjString.h"

class AliAnalysisTaskTrigHMTF : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskTrigHMTF(const char* name="AliAnalysisTaskTest");
  virtual ~AliAnalysisTaskTrigHMTF(){};
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 protected:
  TTree* fTree;       //! output tree
  Int_t fRunNumber;
  UShort_t fBC;
  TObjString fClassesFired; 
  UInt_t fInputsL0;
  Float_t fV0ATime;
  Float_t fV0CTime;
  Bool_t fBBFlag[64];
  Bool_t fBGFlag[64];
  Bool_t fIsFriend;
  ULong64_t fBBFlagPF[21];
  ULong64_t fBGFlagPF[21];
  Int_t fV0ADecision;
  Int_t fV0CDecision;
  Float_t fMTotV0A;
  Float_t fMTotV0C;
  Float_t fMRingV0A[4];
  Float_t fMRingV0C[4];
  UShort_t fTriggerChargeA;
  UShort_t fTriggerChargeC;
  TBits fIR1;
  TBits fIR2;
  UInt_t fNofTracklets;
  Int_t fNofITSClusters[6];
  TBits fFOmap;
  TBits fFiredChipMap;
  Int_t fVertexContributors;
  Int_t fPileupContributors;
  Bool_t fIsPileupSPD;
  Bool_t fIsIncomplete;
  Float_t fT0A[5];
  Float_t fT0C[5];
  Float_t fTVX[5];
  Float_t fADATime;
  Float_t fADCTime;
  Int_t fADADecision;
  Int_t fADCDecision;
  Float_t fMTotADA;
  Float_t fMTotADC;
  UShort_t fTriggerChargeADA;
  UShort_t fTriggerChargeADC;
  Float_t fVz;
  UInt_t fNITSsaTracks;
  UInt_t fNITSsaTracksHits[6]; // Number of tracks having hits in layer x

  ClassDef(AliAnalysisTaskTrigHMTF,1);
};
#endif

