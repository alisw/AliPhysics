#ifndef ALI_ANALYSIS_TASK_SE_LBTOLBPI4_H
#define ALI_ANALYSIS_TASK_SE_LBTOLBPI4_H

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH2F.h>
#include <TArrayD.h>
#include "TClonesArray.h"
#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCuts.h"
#include "AliAODMCHeader.h"
#include "AliVertexingHFUtils.h"

#include "AliAODRecoDecayHF3Prong.h"
#include <TH1F.h>

class TNtuple;
class TGraph;
class TList;
class AliAODTrack;
class TClonesArray;
class TObjArray;
class AliESDVertex;
class AliVVertex;

class AliAnalysisTaskSELbtoLcpi4:public AliAnalysisTaskSE {
public:
  AliAnalysisTaskSELbtoLcpi4();
  AliAnalysisTaskSELbtoLcpi4(const char *name,
                             Bool_t fillntuple,
                             AliRDHFCutsLctopKpi *lccutsanal,
                             AliRDHFCutsLctopKpi *lccutsprod,
                             Int_t ndebug);
  virtual ~AliAnalysisTaskSELbtoLcpi4();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  // virtual void UserCreateOutputHistos();
  virtual void Init();
  //  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  // options to fill Ntuple
  void SetFillNtupleSignal(Bool_t val = kTRUE) {fFillNtupleSignal = val;}
  void SetFillNtupleBackgroundRotated(Bool_t val = kTRUE) {fFillNtupleBackgroundRotated = val;}
  void SetFillNtupleBackgroundNonRotated(Bool_t val = kTRUE) {fFillNtupleBackgroundNonRotated = val;}

  //set parameters
  void SetCutsond0Lcdaughters(Bool_t val = kTRUE) {fCutsond0Lcdaughters = val; return;}

  void SetApplyFixesITS3AnalysisBit(Bool_t val = kTRUE){fApplyFixesITS3AnalysisBit = val;}
  void SetApplyFixesITS3AnalysiskAll(Bool_t val = kTRUE){fApplyFixesITS3AnalysiskAll = val;}
  void SetApplyFixesITS3AnalysisHijing(Bool_t val = kTRUE){fApplyFixesITS3AnalysisHijing = val;}

  void ApplyD0CutLcdaughters(Double_t d0cutd1, Double_t d0cutd2){fCutD0Daughter[0]=d0cutd1; fCutD0Daughter[1]=d0cutd2;}
  void SetPtConfiguration(Double_t ptbin, Double_t ptlcupper, Double_t ptlclower, Double_t ptpionupper, Double_t ptpionlower, Double_t ptlbupper, Double_t ptlblower){fCutsPerPt[0]=ptbin;fCutsPerPt[1]=ptlcupper;fCutsPerPt[2]=ptlclower;fCutsPerPt[3]=ptpionupper;fCutsPerPt[4]=ptpionlower;fCutsPerPt[5]=ptlbupper;fCutsPerPt[6]=ptlblower;}
  void SetNRotations(Double_t nrotations){fNRotations=nrotations;}
  void SetLctopKpiPreselection(Bool_t d){ fPreSelectLctopKpi = d; }

private:
  AliAnalysisTaskSELbtoLcpi4(const AliAnalysisTaskSELbtoLcpi4&);
  AliAnalysisTaskSELbtoLcpi4& operator=(const AliAnalysisTaskSELbtoLcpi4&);

  // Helper functions
  void AddDaughterRefs(AliAODVertex *v,const AliVEvent *event, const TObjArray *trkArray) const;

  Int_t CheckMCLc(AliAODRecoDecayHF3Prong *d, TClonesArray* arrayMC);
  Int_t CheckMCpartPIONaf(AliAODTrack *p, TClonesArray* arrayMC);
  void  FillLbHists(AliAODRecoDecayHF2Prong *part,Int_t lb,AliAODMCHeader *mcHeader,TClonesArray* arrayMC,AliAODTrack *p,AliAODRecoDecayHF3Prong *d, Int_t Lc,AliAODEvent *ev, Bool_t IsPromptLc);
  void  FillLbHistsnr(AliAODRecoDecayHF2Prong *part,Int_t lb,AliAODMCHeader *mcHeader,TClonesArray* arrayMC,AliAODTrack *p,AliAODRecoDecayHF3Prong *d,Int_t Lc,AliAODEvent *ev, Bool_t IsPromptLc);
  void FillHistos(AliAODRecoDecayHF3Prong* d,TClonesArray* arrayMC,AliAODEvent *ev,AliAODMCHeader *mcHeader);
  Bool_t CheckGenerator(AliAODTrack *p,AliAODRecoDecayHF3Prong *d,AliAODMCHeader *mcHeader,TClonesArray* arrayMC);
  Int_t IsSelectedLbMY(TObject* obj,Int_t selectionLevel,Int_t lb,Int_t isRot, Bool_t isHijing) const;
  Bool_t CountLc(AliAODRecoDecayHF3Prong* Lc, AliAODTrack* pion, TClonesArray* arrayMC,Int_t motherLabelLc,Int_t motherLabelpione);
  Bool_t IsCandidateInjected(AliAODRecoDecayHF *part, AliAODMCHeader *header,TClonesArray *arrayMC);
  Int_t IsTrackInjected(AliAODTrack *part,AliAODMCHeader *header,TClonesArray *arrayMC);
  AliAODVertex* RecalculateVertex(const AliVVertex *primary,TObjArray *tracks,Double_t bField);
  AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray,Double_t &dispersion,Bool_t useTRefArray=kTRUE) const; //Reconstruct vertex of B
  void DoRotations(AliAODEvent* ev, AliAODRecoDecayHF2Prong *decay, AliAODRecoDecayHF3Prong* lc3prong, AliAODTrack* piontrack, Int_t nRot, Bool_t isHijing, Int_t lb, TClonesArray* arrayMC, AliAODMCHeader *mcHeader);

  AliPIDResponse *fPIDResponse;
  TList   *fOutput; //! list send on output slot 3
  TH1F    *fHistNEvents; //!hist. for No. of events
  AliRDHFCutsLctopKpi *fRDCutsAnalysisLc; //Cuts for Analysis
  AliRDHFCutsLctopKpi *fRDCutsProductionLb; //Production Cuts
  TList *fListCuts; //list of cuts

  Double_t fBzkG;                // z component of magnetic field
  AliAODVertex *fvtx1;           // primary vertex
  TList    *fOutputList;         //! output slot 8
  TNtuple *fNtupleLambdabUPG; //! output ntuple

  Bool_t fFillNtupleSignal; /// flag to fill ntuple with signal candidates
  Bool_t fFillNtupleBackgroundRotated; /// flag to fill ntuple with background rotated candidates
  Bool_t fFillNtupleBackgroundNonRotated; /// flag to fill ntuple with background non rotated candidates
  Double_t fCutD0Daughter[2];
  Bool_t fCutsond0Lcdaughters;
  Double_t fCutsPerPt[7];
  Double_t fNRotations;
  Bool_t fIsPromptLc;
  Bool_t fPreSelectLctopKpi;

  //temporary to compare with old ITS2 analysis
  Bool_t fApplyFixesITS3AnalysisBit;
  Bool_t fApplyFixesITS3AnalysiskAll;
  Bool_t fApplyFixesITS3AnalysisHijing;

  ClassDef(AliAnalysisTaskSELbtoLcpi4,8);
};

#endif
