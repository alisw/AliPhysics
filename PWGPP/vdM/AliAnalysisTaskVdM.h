// -*- C++ -*-
#ifndef ALIANALYSISTASKVDM_H
#define ALIANALYSISTASKVDM_H

class TH1;
class TTree;
class TList;

class AliVEvent;
class AliESDEvent;
class AliESDHeader;

#include <TObject.h>
#include <TString.h>
#include <TBits.h>
#include <TClonesArray.h>

#include "AliESDVertex.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"


// task for VdM analysis producting a TTree with
//  * constrained and unconstrained vertex -> non-separation analysis
//  * timing information for V0 and for AD -> bkgd estimation
//
class AliAnalysisTaskVdM : public AliAnalysisTaskSE {
public:

  AliAnalysisTaskVdM(const char *name="AliAnalysisTaskVdM");
  virtual ~AliAnalysisTaskVdM();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t *);
  virtual void NotifyRun();

  void SetBranchNames(TString options) { fTreeBranchNames = options; }

  TString GetListName() const { return "TL"; }
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
      , fL2Inputs(0)
      , fOrbitID(0) {}

    void Fill(const AliVEvent *);

    ULong64_t fClassMask;
    ULong64_t fClassMaskNext50;
    UInt_t    fBCID;
    UInt_t    fPeriod;
    UInt_t    fTimeStamp;
    UInt_t    fL0Inputs;
    UInt_t    fL1Inputs;
    Int_t     fRunNumber;
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
        fBB[i] = fBG[i] = -1;
        fDecisionOnline[i] = fDecisionOffline[i] = -1;
      }
      for (Int_t bc=0; bc<21; ++bc)
        fPFBBA[bc] = fPFBBC[bc] = fPFBGA[bc] = fPFBGC[bc] = 0;
    }

    void FillAD(const AliVEvent *, AliTriggerAnalysis &);
    void FillV0(const AliVEvent *, AliTriggerAnalysis &);

    void FillInvalid();

    Float_t    fTime[2];            //
    Char_t     fBB[2];              //
    Char_t     fBG[2];              //
    Double32_t fDecisionOnline[2];  //[-1,3,2]
    Double32_t fDecisionOffline[2]; //[-1,3,2]
    Double32_t fPFBBA[21];          //[0,32,5]
    Double32_t fPFBBC[21];          //[0,32,5]
    Double32_t fPFBGA[21];          //[0,32,5]
    Double32_t fPFBGC[21];          //[0,32,5]
  } ;

  class TreeData : public TObject {
  public:
    TreeData()
      : TObject()
      , fEventInfo()
      , fV0Info()
      , fADInfo()
      , fIsIncompleteDAQ(kFALSE) {}

    EventInfo fEventInfo;
    ADV0      fV0Info;
    ADV0      fADInfo;
    Bool_t    fIsIncompleteDAQ;
    ClassDef(TreeData, 1);
  } ;

protected:
  void SetBranches(TTree* t);
  void FillTriggerIR(const AliESDHeader* );
  enum {
    kHistTrig,
    kNHist
  };

private:
  AliAnalysisTaskVdM(const AliAnalysisTaskVdM&); // not implemented
  AliAnalysisTaskVdM& operator=(const AliAnalysisTaskVdM&); // not implemented

  TString          fTreeBranchNames;     //

  AliTriggerAnalysis fTriggerAnalysis;   //!

  TList           *fList;                //!
  TH1             *fHist[kNHist];        //!
  TTree           *fTE;                  //!
  TBits            fIR1InteractionMap;   //!
  TBits            fIR2InteractionMap;   //!
  AliESDVertex     fVertexSPD;           //!
  AliESDVertex     fVertexTPC;           //!
  AliESDVertex     fVertexTracks;        //!
  AliESDVertex     fVertexTracksUnconstrained; //!
  TClonesArray     fTriggerIRs;          //!
  TString          fFiredTriggerClasses; //!
  TreeData         fTreeData;            //!

  ClassDef(AliAnalysisTaskVdM, 2);
} ;

#endif // ALIANALYSISTASKVDM_H
