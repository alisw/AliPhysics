// -*- C++ -*-
#ifndef ALIANALYSISTASKDG_H
#define ALIANALYSISTASKDG_H

class TH1;
class TTree;
class TList;

class AliVEvent;
class AliESDEvent;
class AliESDHeader;
class AliVTrack;
class AliESDtrackCuts;

#include <algorithm>

#include <TObject.h>
#include <TString.h>
#include <TBits.h>
#include <TClonesArray.h>

#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliAODVZERO.h"
#include "AliESDVZERO.h"
#include "AliAODAD.h"
#include "AliESDAD.h"
#include "AliAnalysisTaskSE.h"
#include "AliTOFHeader.h"
#include "AliTriggerAnalysis.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliAnalysisUtils.h"

class AliAnalysisTaskDG : public AliAnalysisTaskSE {
public:

  AliAnalysisTaskDG(const char *name="AliAnalysisTaskDG");
  virtual ~AliAnalysisTaskDG();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t *);
  virtual void NotifyRun();

  void SetIsMC(Bool_t b=kTRUE) { fIsMC = b; }
  void SetBranchNames(TString options) { fTreeBranchNames = options; }
  void SetTrackCutType(TString tc) { fTrackCutType = tc; }
  void SetTriggerSelection(TString ts) { fTriggerSelection = ts; }
  void SetTriggerSelectionSPD(TString ts) { fTriggerSelectionSPD = ts; }
  void SetMaxTracksSave(Int_t m);

  TString GetListName() const { return fTrackCutType+"_TL"; }
  TString GetTreeName() const { return fTrackCutType+"_TE"; }
  TString GetResultsFileName() const { return "results.root"; }

  struct EventInfo {
    EventInfo()
      : fClassMask(0)
      , fClassMaskNext50(0)
      , fBCID(0)
      , fPeriod(0)
      , fTimeStamp(0)
      , fL0Inputs(0)
      , fL1Inputs(0)
      , fRunNumber(0)
      , fnTrk(0)
      , fCharge(0)
      , fL2Inputs(0)
      , fOrbitID(0) {
      std::fill_n(fnTrklet, 4, 0);
      fnFO[0] = fnFO[1] = 0;
    }

    void Fill(const AliVEvent *);

    ULong64_t fClassMask;
    ULong64_t fClassMaskNext50;
    UInt_t    fBCID;
    UInt_t    fPeriod;
    UInt_t    fTimeStamp;
    UInt_t    fL0Inputs;
    UInt_t    fL1Inputs;
    Int_t     fRunNumber;
    UShort_t  fnTrk;
    UShort_t  fnTrklet[4]; // all, C,cent,A
    UShort_t  fnFO[2];     // inner,outer layer
    Short_t   fCharge;
    UShort_t  fL2Inputs;
    UShort_t  fOrbitID;

  } ;

  struct ADV0 {
    enum {
      kCside = 0,
      kAside = 1
    };

    ADV0() {
      for (Int_t i=0; i<2; ++i) {
	fTime[i] = -10240.0f;
	fBB[i] = fBG[i];
	fDecisionOnline[i] = fDecisionOffline[i] = -1;
      }
      std::fill_n(fMult,   8, -1);
      std::fill_n(fPFBBA, 21,  0);
      std::fill_n(fPFBBC, 21,  0);
      std::fill_n(fPFBGA, 21,  0);
      std::fill_n(fPFBGC, 21,  0);
    }

    void FillAD(const AliVEvent *, AliTriggerAnalysis &);
    void FillV0(const AliVEvent *, AliTriggerAnalysis &);

    void FillInvalid();

    Float_t    fTime[2];            //
    Short_t    fBB[2];              //
    Short_t    fBG[2];              //
    Double32_t fDecisionOnline[2];  //[-1,3,2]
    Double32_t fDecisionOffline[2]; //[-1,3,2]
    Double32_t fPFBBA[21];          //[0,32,5]
    Double32_t fPFBBC[21];          //[0,32,5]
    Double32_t fPFBGA[21];          //[0,32,5]
    Double32_t fPFBGC[21];          //[0,32,5]
    Float_t    fMult[8];            // multiplicity per ring
  } ;

  struct FMD {
    FMD()
      : fA(kFALSE)
      , fC(kFALSE) {}

    void Fill(const AliVEvent *, AliTriggerAnalysis &);

    Bool_t fA;
    Bool_t fC;
  };

  struct ZDC {
    ZDC()
      : fZNenergy()
      , fZPenergy()
      , fZEMenergy()
      , fZNtower0()
      , fZPtower0()
      , fZNTDC() {}

    void Fill(AliVZDC*);

    Float_t fZNenergy[2];
    Float_t fZPenergy[2];
    Float_t fZEMenergy[2];
    Float_t fZNtower0[2];
    Float_t fZPtower0[2];
    Float_t fZNTDC[2][4];
  } ;

  class TreeData : public TObject {
  public:
    TreeData()
      : TObject()
      , fEventInfo()
      , fV0Info()
      , fADInfo()
      , fFMDInfo()
      , fZDCInfo()
      , fIsIncompleteDAQ(kFALSE)
      , fIsSPDClusterVsTrackletBG(kFALSE)
      , fIskMB(kFALSE) {}

    EventInfo fEventInfo;
    ADV0      fV0Info;
    ADV0      fADInfo;
    FMD       fFMDInfo;
    ZDC       fZDCInfo;
    Bool_t    fIsIncompleteDAQ;
    Bool_t    fIsSPDClusterVsTrackletBG;
    Bool_t    fIskMB;
    ClassDef(TreeData, 9);
  } ;

  struct TrackData : public TObject {
    TrackData(AliVTrack *tr=NULL, AliPIDResponse *pidResponse=NULL)
      : TObject()
      , fSign(0)
      , fPx(0)
      , fPy(0)
      , fPz(0)
      , fLength(0)
      , fITSsignal(0)
      , fTPCsignal(0)
      , fTOFsignal(0)
      , fFilterMap(0)
      , fFlags(0) {
      fPIDStatus[0] = fPIDStatus[1] = fPIDStatus[2] = AliPIDResponse::kDetNoSignal;
      const Int_t nSpecies = AliPID::kSPECIES;
      std::fill_n(fNumSigmaITS, nSpecies, -32.0f);
      std::fill_n(fNumSigmaTPC, nSpecies, -32.0f);
      std::fill_n(fNumSigmaTOF, nSpecies, -32.0f);
      fChipKey[0] = fChipKey[1] = -1;
      fStatus[0]  = fStatus[1]  = -1;
      Fill(tr, pidResponse);
    }

    void Fill(AliVTrack *, AliPIDResponse*);

    Double32_t fSign;                          //[-1,1,2]
    Float_t    fPx,fPy,fPz;
    Float_t    fLength;
    Float_t    fITSsignal, fTPCsignal, fTOFsignal;
    Double32_t fNumSigmaITS[AliPID::kSPECIES]; //[-32,32,8]
    Double32_t fNumSigmaTPC[AliPID::kSPECIES]; //[-32,32,8]
    Double32_t fNumSigmaTOF[AliPID::kSPECIES]; //[-32,32,8]
    Double32_t fPIDStatus[3];                  //[0,4,2] ITS,TPC,TOF
    Short_t    fChipKey[2];                    // L0,L1 (SPD)
    Int_t      fStatus[2];                     // L0,L1 (SPD)
    UInt_t     fFilterMap;
    ULong_t    fFlags;
    ClassDef(TrackData, 6);
  } ;

  struct SPD_0STG : public TObject {
    SPD_0STG()
      : TObject()
      , fMinDeltaPhi(-1)
      , fMaxDeltaPhi(-1)
      , fNPseudoTracklets(-1) {}
    virtual ~SPD_0STG() {}

    const TBits& Fill(const TBits& );

    Short_t fMinDeltaPhi;      // minimum opening angle [1-10]
    Short_t fMaxDeltaPhi;      // maximum opening angle [1-10]
    Short_t fNPseudoTracklets; // number of pseudo-tracklets

    ClassDef(SPD_0STG, 1);
  } ;

protected:
  void SetBranches(TTree* t, Bool_t isAOD);
  void SetClassMask(TString triggerSel, ULong64_t &mask, ULong64_t &maskNext50);
  static void FindChipKeys(AliESDtrack *tr, Short_t chipKeys[2], Int_t status[2]);

  void FillSPDFOEffiencyHistograms(const AliESDEvent* );
  void FillTH3(Int_t idx, Double_t x, Double_t y, Double_t z, Double_t w=1);

  void FillTriggerIR(const AliESDHeader* );

  enum {
    kHistTrig,
    kHistSPDFiredTrk,
    kHistSPDFOTrk,
    kHistSPDFOFiredTrk,
    kHistSPDFiredTrkVsMult,
    kHistSPDFOTrkVsMult,
    kHistSPDFOFiredTrkVsMult,
    kHistSPDFiredVsMult,
    kHistSPDFOVsMult,
    kHistSPDFOFiredVsMult,
    kNHist
  };

private:
  AliAnalysisTaskDG(const AliAnalysisTaskDG&); // not implemented
  AliAnalysisTaskDG& operator=(const AliAnalysisTaskDG&); // not implemented

  Bool_t           fIsMC;                //
  TString          fTreeBranchNames;     //
  TString          fTrackCutType;        //
  UInt_t           fTrackFilterMask;     //
  TString          fTriggerSelection;    //
  TString          fTriggerSelectionSPD; //
  Int_t            fMaxTracksSave;       //

  AliTriggerAnalysis fTriggerAnalysis;   //!
  AliAnalysisUtils   fAnalysisUtils;     //!

  TList           *fList;                //!
  TH1             *fHist[kNHist];        //!
  TTree           *fTE;                  //!
  TBits            fIR1InteractionMap;   //!
  TBits            fIR2InteractionMap;   //!
  TBits            fFastOrMap;           //!
  TBits            fFiredChipMap;        //!
  typedef std::pair<AliESDVertex, AliAODVertex> VtxPairType;
  VtxPairType      fVertexSPD;           //!
  VtxPairType      fVertexTPC;           //!
  VtxPairType      fVertexTracks;        //!
  typedef std::pair<AliESDVZERO, AliAODVZERO> V0PairType;
  V0PairType       fV0;                  //!
  typedef std::pair<AliESDAD, AliAODAD> ADPairType;
  ADPairType       fAD;                  //!
  AliTOFHeader     fTOFHeader;           //!
  TClonesArray     fTriggerIRs;          //!
  TString          fFiredTriggerClasses; //!
  TreeData         fTreeData;            //!
  SPD_0STG         fSPD_0STG_Online;     //! using FastOrMap    (online)
  SPD_0STG         fSPD_0STG_Offline;    //! using FiredChipMap (offline)
  TClonesArray     fTrackData;           //!
  TClonesArray     fMCTracks;            //!
  AliESDtrackCuts *fTrackCuts;           //!

  ClassDef(AliAnalysisTaskDG, 17);
} ;

#endif // ALIANALYSISTASKDG_H
