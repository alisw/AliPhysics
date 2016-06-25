// -*- C++ -*-
#ifndef ALIANALYSISTASKDIFFCROSSSECTIONS_H
#define ALIANALYSISTASKDIFFCROSSSECTIONS_H

class TH1;
class TTree;
class TList;

class AliESDHeader;
class AliESDAD;
class AliESDVZERO;
class AliMCEvemt;
class AliESDEvent;

#include <TObject.h>
#include <TString.h>
#include <TBits.h>
#include <TClonesArray.h>

#include "AliESDVertex.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"
#include "AliAnalysisUtils.h"

class AlianalysisTaskDiffCrossSections : public AliAnalysisTaskSE {
public:

  AlianalysisTaskDiffCrossSections(const char *name="AlianalysisTaskDiffCrossSections");
  virtual ~AlianalysisTaskDiffCrossSections();
  
  virtual void NotifyRun();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t *);

  void SetIsMC(Bool_t b=kTRUE) { fIsMC = b; }
  void SetMCType(TString s) { fMCType = s; }
  void SetTriggerSelection(TString ts) { fTriggerSelection = ts; }

  TString GetTreeName() const { return "TE"; }
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
      , fnTrklet(0)
      , fOrbitID(0) {
      fnSPDClusters[0] = fnSPDClusters[1] = 0;
    }

    void Fill(const AliESDEvent *);

    ULong64_t fClassMask;
    ULong64_t fClassMaskNext50;
    UInt_t    fBCID;
    UInt_t    fPeriod;
    UInt_t    fTimeStamp;
    UInt_t    fL0Inputs;
    UInt_t    fL1Inputs;
    UInt_t    fnSPDClusters[2]; // 0 -> inner layer, 1 -> outer layer
    Int_t     fRunNumber;
    UShort_t  fnTrklet;
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
	fCharge[i] = 0.0;
	fBB[i] = fBG[i] = -1;
	fDecisionOnline[i] = fDecisionOffline[i] = -1;
      }
    }

    void FillAD(const AliESDEvent *, AliTriggerAnalysis &);
    void FillV0(const AliESDEvent *, AliTriggerAnalysis &);

    void FillInvalid();

    Float_t    fTime[2];            //
    Float_t    fCharge[2];          //
    Char_t     fBB[2];              // 
    Char_t     fBG[2];              //
    Double32_t fDecisionOnline[2];  //[-1,3,2]
    Double32_t fDecisionOffline[2]; //[-1,3,2]
  } ;

  class TreeData : public TObject {
  public:
    TreeData() 
      : TObject()
      , fEventInfo()
      , fV0Info()
      , fADInfo()
      , fPhysSelBits(0)
      , fIsIncompleteDAQ(kFALSE)
      , fIsSPDClusterVsTrackletBG(kFALSE) {}

    EventInfo fEventInfo;
    ADV0      fV0Info;
    ADV0      fADInfo;
    UInt_t    fPhysSelBits;
    Bool_t    fIsIncompleteDAQ;
    Bool_t    fIsSPDClusterVsTrackletBG;
    ClassDef(TreeData, 1);
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
      : TObject() {}
    virtual ~MCInfo() {}

    void Fill(const AliMCEvent *, TString);

    Float_t fEventType;    //[-3,5,3]
    Float_t fDiffMass[2];  // 0 -> L, 1 -> R
    ClassDef(MCInfo, 1);
  } ;

protected:
  void SetBranches(TTree* t);

private:  
  AlianalysisTaskDiffCrossSections(const AlianalysisTaskDiffCrossSections&); // not implemented
  AlianalysisTaskDiffCrossSections& operator=(const AlianalysisTaskDiffCrossSections&); // not implemented

  Bool_t           fIsMC;                //
  TString          fMCType;              //
  TString          fTriggerSelection;    //

  AliTriggerAnalysis fTriggerAnalysis;   //!
  AliAnalysisUtils   fAnalysisUtils;     //!

  TTree           *fTE;                  //!
  TBits            fFastOrMap;           //!
  TBits            fFiredChipMap;        //!
  AliESDVertex     fVertexSPD;           //!
  TreeData         fTreeData;            //!
  MCInfo           fMCInfo;              //!
  
  ClassDef(AlianalysisTaskDiffCrossSections, 1);
} ;

#endif // ALIANALYSISTASKDIFFCROSSSECTIONS_H
