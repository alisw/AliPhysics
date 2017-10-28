// -*- C++ indent-tabs-mode:nil; -*-
#ifndef ALIANALYSISTASKDIFFCROSSSECTIONS_H
#define ALIANALYSISTASKDIFFCROSSSECTIONS_H

class TH1;
class TTree;
class TList;

class AliESDFMD;
class AliMCEvemt;
class AliESDEvent;
class AliStack;

#include <TObject.h>
#include <TString.h>
#include <TBits.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TCutG.h>
#include <TMatrixD.h>

#include "AliESDVertex.h"
#include "AliESDVZERO.h"
#include "AliESDAD.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"
#include "AliAnalysisUtils.h"

class AliAnalysisTaskDiffCrossSections : public AliAnalysisTaskSE {
public:

  AliAnalysisTaskDiffCrossSections(const char *name="AliAnalysisTaskDiffCrossSections");
  virtual ~AliAnalysisTaskDiffCrossSections();

  virtual void NotifyRun();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t *);

  void SetIsMC(Bool_t b=kTRUE) { fIsMC = b; }
  void SetMCType(TString s) { fMCType = s; }
  void SetTriggerSelection(TString ts) { fTriggerSelection = ts; }

  void SetDetectorsUsed(TString det) { fDetectorsUsed = det; }
  void SetUseBranch(TString b)       { fUseBranch = b; }

  void SetFMDMultLowCut(Float_t cut) { fFMDMultLowCut = cut; }
  TString GetTreeName() const {
    TString s = "TE";
    if (!fIsMC && fTriggerSelection != "") {
      s += fTriggerSelection;
      s.ReplaceAll("|", "_");
    }
    return s;
  }
  TString GetResultsFileName() const { return "results.root"; }

  void SetTimeChargeCut(Int_t ch, TCutG* cut) {
    fTimeChargeCuts.AddAt(cut, ch);
  }

  Bool_t DoTimeChargeCut(Int_t ch, Float_t time, Float_t charge) const;

  struct EventInfo {
    EventInfo()
      : fClassMask(0)
      , fClassMaskNext50(0)
      , fBCID(0)
      , fPeriod(0)
      , fTimeStamp(0)
      , fL0Inputs(0)
      , fL1Inputs(0)
      , fEventNumberInFile(0)
      , fRunNumber(0)
      , fnTrklet(0)
      , fL2Inputs(0)
      , fOrbitID(0)
      , fInputFileName("") {
      fnSPDClusters[0] = fnSPDClusters[1] = 0;
    }

    void Fill(const AliESDEvent *, const char*);

    ULong64_t fClassMask;
    ULong64_t fClassMaskNext50;
    UInt_t    fBCID;
    UInt_t    fPeriod;
    UInt_t    fTimeStamp;
    UInt_t    fL0Inputs;
    UInt_t    fL1Inputs;
    UInt_t    fnSPDClusters[2]; // 0 -> inner layer, 1 -> outer layer
    Int_t     fEventNumberInFile;
    Int_t     fRunNumber;
    UShort_t  fnTrklet;
    UShort_t  fL2Inputs;
    UShort_t  fOrbitID;
    TString   fInputFileName;
  } ;

  struct ADV0 {
    enum {
      kCside = 0,
      kAside = 1
    };

    ADV0() {
      FillInvalid();
    }

    void FillAD(const AliESDEvent *, AliTriggerAnalysis &);
    void FillV0(const AliESDEvent *, AliTriggerAnalysis &);

    void FillInvalid();

    Float_t    fTime[2];            //
    Float_t    fCharge[2];          //
    UChar_t    fBBOnline[2];        // number of BB flags (for AD: coincidences) per side
    UChar_t    fBBOffline[2];       // number of BB flags (for AD: coincidences) per side
    UChar_t    fBGOnline[2];        // number of BG flags (for AD: coincidences) per side
    UChar_t    fBGOffline[2];       // number of BG flags (for AD: coincidences) per side
    Double32_t fDecisionOnline[2];  //[-1,3,2]
    Double32_t fDecisionOffline[2]; //[-1,3,2]
    UChar_t    fBBFlagsOnline[64];  //
    UChar_t    fBBFlagsOffline[64]; //
    UShort_t   fCharges[64];        //
    Float_t    fTimes[64];          //
  } ;

  struct VtxInfo {
    VtxInfo()
      : fX(0)
      , fY(0)
      , fZ(0)
      , fNcontr(-4) {}

    void Fill(const AliESDVertex*);

    TVector3 GetVtxPosition() const {
      return TVector3(fX, fY, fZ);
    }

    Float_t fX;
    Float_t fY;
    Float_t fZ;
    Char_t  fNcontr;
  } ;

  struct FMDInfo {
    FMDInfo() {
      for (Int_t i=0; i<5; ++i)
       	fMult[i] = 0;
    }
    void Fill(const AliESDEvent*, Float_t);

    Int_t   fMult[5];   // A-side: 1,2i,2o; C-side: 3i,3o
  } ;

  // single pseudo-track
  class PseudoTrack : public TObject {
  public:
    PseudoTrack(Float_t eta=0, Float_t phi=0, Float_t charge=0, Float_t charge2=0, UInt_t detFlags=0)
      : TObject()
      , fEta(eta)
      , fPhi(phi)
      , fCharge(charge)
      , fCharge2(charge2)
      , fDetFlags(detFlags) {}

    PseudoTrack(const PseudoTrack& t)
      : TObject(t)
      , fEta(t.fEta)
      , fPhi(t.fPhi)
      , fCharge(t.fCharge)
      , fCharge2(t.fCharge2)
      , fDetFlags(t.fDetFlags) {}

    virtual ~PseudoTrack() {}

    virtual void   Print(Option_t *option="") const;
    virtual Int_t  Compare(const TObject *obj) const;
    virtual Bool_t IsSortable() const { return kTRUE; }

    Float_t Eta() const { return fEta; }
    Float_t Phi() const { return fPhi; }
    UInt_t  DetFlags() const { return fDetFlags; }
    Float_t Charge() const { return fCharge; }
    Float_t Charge2() const { return fCharge2; }
    Bool_t  operator<(const PseudoTrack& t) const { return fEta < t.fEta; }
  protected:
  private:
    PseudoTrack& operator=(const PseudoTrack&); // not implemented

    Float_t fEta;      // pseudo-rapidity
    Float_t fPhi;      // phi (rad)
    Float_t fCharge;   // charge (=1 for SPD)
    Float_t fCharge2;  // charge (=0 for SPD,FMD,VZERO) for 2nd layer of AD
    UInt_t  fDetFlags; // flags indicating from which detector this pseudo-track is + online/offline information
    ClassDef(PseudoTrack, 1);
  } ;

  class PseudoTracks : public TObject {
  public:
    enum { // used in PseudoTrack::fDetFlags
      kSPD     = (1<< 0),
      kV0C     = (1<< 1),
      kV0A     = (1<< 2),
      kV0      = kV0C | kV0A,
      kADC     = (1<< 3),
      kADA     = (1<< 4),
      kAD      = kADC | kADA,
      kFMD1    = (1<< 5),
      kFMD2i   = (1<< 6),
      kFMD2o   = (1<< 7),
      kFMD3i   = (1<< 8),
      kFMD3o   = (1<< 9),
      kFMD     = kFMD1 | kFMD2i | kFMD2o | kFMD3i | kFMD3o,
      kFlagNotBB = (1<<29),
      kOnline  = (1<<30),
      kOffline = (1<<31)
    };

    PseudoTracks();
    virtual ~PseudoTracks();

    inline Int_t GetEntries() const { return fTracks.GetEntriesFast(); }
    void  AddTrack(const PseudoTrack&);
    virtual void  Clear(Option_t *opt="");
    virtual void  Delete(Option_t *opt="");
    const PseudoTrack& operator[](Int_t ) const;
    inline const PseudoTrack& GetTrackAt(Int_t i) const { return this->operator[](i); }

    Int_t RemoveTracks(UInt_t mask); // removes all tracks matching mask; returns number of removed pseudo-tracks

    void  FindAcceptance(UInt_t mask, Float_t &etaAccL, Float_t &etaAccR, const TVector3 &vertexPosition) const;
    void  FindAcceptance(UInt_t mask, Float_t &etaAccL, Float_t &etaAccR) const;

    template<typename F> // F is a function (object) of type Bool_t f(const PseudoTrack&)
    Int_t ClassifyEvent(Int_t &iEtaL, Int_t &iEtaR, Float_t &etaGap, Float_t &etaGapCenter,
			UInt_t mask, F& f) {
      static TBits bits(10000);
      SortIfNeeded();
      for (Int_t i=0, nt=fTracks.GetEntriesFast(); i<nt; ++i)
	bits.SetBitNumber(i, f(GetTrackAt(i)));
      return ClassifyEventBits(iEtaL, iEtaR, etaGap, etaGapCenter, mask, bits);
    }
  protected:
    void SortIfNeeded();

    Int_t ClassifyEventBits(Int_t &iEtaL, Int_t &iEtaR, Float_t &etaGap, Float_t &etaGapCenter,
			    UInt_t mask, const TBits& bits) const;
  private:
    PseudoTracks(const PseudoTracks& ); // not implemented
    PseudoTracks& operator=(const PseudoTracks& ); // not implemented

    mutable TClonesArray fTracks;
    ClassDef(PseudoTracks, 3);
  } ;

  static TVector3 GetADPseudoTrack(Int_t ch);
  static TVector3 GetV0PseudoTrack(Int_t ch);

  class TreeData : public TObject {
  public:
    TreeData()
      : TObject()
      , fEventInfo()
      , fVtxInfo()
      , fV0Info()
      , fFMDInfo()
      , fADInfo()
      , fPseudoTracks()
      , fPhysSelBits(0)
      , fIsIncompleteDAQ(kFALSE)
      , fIsSPDClusterVsTrackletBG(kFALSE) {}

    EventInfo    fEventInfo;
    VtxInfo      fVtxInfo;
    ADV0         fV0Info;
    FMDInfo      fFMDInfo;
    ADV0         fADInfo;
    PseudoTracks fPseudoTracks;
    UInt_t       fPhysSelBits;
    Bool_t       fIsIncompleteDAQ;
    Bool_t       fIsSPDClusterVsTrackletBG;
  private:
    TreeData(const TreeData&); // not implemented
    TreeData& operator=(const TreeData&); // not implemented

    ClassDef(TreeData, 9);
  } ;

  class MCInfo : public TObject {
  public:
    enum {
      kInvalid = -1,
      kSDL,
      kSDR,
      kDD,
      kCD,
      kND,
      kElastic
    };

    MCInfo()
      : TObject()
      , fEventType(kInvalid)
      , fEventTypeGen(kInvalid)
    {
      for (Int_t i=0; i<2; ++i)
	fDiffMass[i] = fDiffMassGen[i] = -1.0;
    }
    virtual ~MCInfo() {}

    void   Fill(const AliMCEvent *, TString);
    Bool_t FindSingleDiffraction(AliStack *stack, TString mcType,
				 Int_t &side, Double_t &mass) const;

    Float_t  fEventType;      //[-3,5,3]
    Float_t  fEventTypeGen;   //[-3,5,3]
    Double_t fDiffMass[2];    // 0 -> L, 1 -> R
    Double_t fDiffMassGen[2]; // 0 -> L, 1 -> R
    ClassDef(MCInfo, 2);
  } ;

protected:
  void SetBranches(TTree* t);
  TVector3 GetRandomVtxPosition() const;

private:
  AliAnalysisTaskDiffCrossSections(const AliAnalysisTaskDiffCrossSections&); // not implemented
  AliAnalysisTaskDiffCrossSections& operator=(const AliAnalysisTaskDiffCrossSections&); // not implemented

  Bool_t           fIsMC;                //
  TString          fMCType;              //
  TString          fTriggerSelection;    //
  TString          fDetectorsUsed;       //
  TString          fUseBranch;           //
  Float_t          fFMDMultLowCut;       //
  TObjArray        fTimeChargeCuts;      // TCutG (time, charge) -> Bool_t for each channel

  AliTriggerAnalysis fTriggerAnalysis;   //!
  AliAnalysisUtils   fAnalysisUtils;     //!

  TTree           *fTE;                  //!
  TBits            fFastOrMap;           //!
  TBits            fFiredChipMap;        //!
  AliESDVZERO      fESDVZERO;            //!
  AliESDAD         fESDAD;               //!
  TreeData         fTreeData;            //!
  MCInfo           fMCInfo;              //!

  TBits            fIR1InteractionMap;   //!
  TBits            fIR2InteractionMap;   //!

  TVectorD         fMeanVtxPos;          //!
  TMatrixD         fMeanVtxCov;          //!
  TMatrixD         fMeanVtxU;            //!

  Int_t            fEventType;    //! -1 -> 1L, +1 -> 1R, +2 -> 2A
  Int_t            fIdxEtaL;      //! index of etaL in fTreeData.fTracks
  Int_t            fIdxEtaR;      //! index of etaR in fTreeData.fTracks
  Float_t          fEtaL;         //!
  Float_t          fEtaR;         //!
  Float_t          fEtaGap;       //!
  Float_t          fEtaGapCenter; //!
  ClassDef(AliAnalysisTaskDiffCrossSections, 4);
} ;

#endif // ALIANALYSISTASKDIFFCROSSSECTIONS_H
