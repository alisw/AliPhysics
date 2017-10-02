/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

 /******************************** 
 * femtoscopy with multiparticle *
 *           technology          * 
 *                               * 
 * author: Ante Bilandzic        * 
 *        (abilandzic@gmail.com) *
 ********************************/ 

#ifndef ALIANALYSISTASKMULTIPARTICLEFEMTOSCOPY_H
#define ALIANALYSISTASKMULTIPARTICLEFEMTOSCOPY_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "TExMap.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TH1I.h"

//================================================================================================================

class AliAnalysisTaskMultiparticleFemtoscopy : public AliAnalysisTaskSE{
 public:
  AliAnalysisTaskMultiparticleFemtoscopy();
  AliAnalysisTaskMultiparticleFemtoscopy(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskMultiparticleFemtoscopy(); 
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

  // 0.) Methods called in the constructor:
  virtual void InitializeArrays(); // use this method temporarily for all objects not classified yet
  virtual void InitializeArraysForControlHistograms();
  virtual void InitializeArraysForEBEObjects();
  virtual void InitializeArraysForCorrelationFunctions();
  virtual void InitializeArraysForBackground();
  virtual void InitializeArraysForBuffers();
  virtual void InitializeArraysForQA();
  virtual void InitializeArraysForGlobalTrackCuts();
  virtual void InitializeArraysForCorrelationFunctionsTEST();
  virtual void InitializeArraysForBackgroundTEST();
  virtual void InitializeArraysForHybridApproach();

  // 1.) Methods called in UserCreateOutputObjects():
  //  2a) Directly:
  virtual void InsanityChecksUserCreateOutputObjects();
  virtual void BookAndNestAllLists();
  virtual void BookEverything(); // use this to book all un-classified objects
  virtual void BookEverythingForControlHistograms();
  virtual void BookEverythingForEBEObjects();
  virtual void BookEverythingForCorrelationFunctions();
  virtual void BookEverythingForBackground();
  virtual void BookEverythingForBuffers();
  virtual void BookEverythingForQA();
  virtual void BookEverythingForGlobalTrackCuts();
  virtual void BookEverythingForCorrelationFunctionsTEST();
  virtual void BookEverythingForBackgroundTEST();
  virtual void BookEverythingForHybridApproach();
  //  2b) Indirectly:
  Int_t InsanityChecksForGlobalTrackCuts(); // insanity checks for global track cuts

  // 2.) Methods called in UserExec(Option_t *):
  //  2a) Directly:
  virtual void InsanityChecksUserExec();
  virtual void QA(AliVEvent *ave);
  virtual void MC(AliMCEvent *aMC);
  virtual void ESD(AliESDEvent *aESD);
  virtual void AOD(AliAODEvent *aAOD);
  virtual void OnlineMonitoring();
  //  2b) Indirectly:
  virtual void EstimateBackground(AliVEvent *ave);
  virtual void EstimateBackgroundTEST(AliVEvent *ave);
  virtual void DoHybridApproach(AliVEvent *ave);
  virtual void FillControlHistogramsEvent(AliVEvent *ave);
  virtual void FillControlHistogramsParticle(AliVEvent *ave);
  virtual void FillControlHistogramsNonIdentifiedParticles(AliAODTrack *atrack);
  virtual void FillControlHistogramsNonIdentifiedParticles(AliAODMCParticle *amcparticle);
  virtual void FillControlHistogramsNonIdentifiedParticlesFTSF(AliAODTrack *atrack);
  virtual void FillControlHistogramsIdentifiedParticles(AliAODTrack *atrack, AliAODTrack *gtrack);
  virtual void FillControlHistogramsIdentifiedParticles(AliAODMCParticle *amcparticle);
  virtual void V0s(AliVEvent *ave);
  Int_t InsanityChecksForTracks(AliAODTrack *atrack); // insanity checks for each track ('atrack') in AOD
  Int_t InsanityChecksForGlobalTracks(AliAODTrack *gtrack); // insanity checks only for global tracks ('gtrack') in AOD
  Bool_t Pion(AliAODTrack *atrack, Int_t charge = 1, Bool_t bPrimary = kTRUE);
  Bool_t Kaon(AliAODTrack *atrack, Int_t charge = 1, Bool_t bPrimary = kTRUE);
  Bool_t Proton(AliAODTrack *atrack, Int_t charge = 1, Bool_t bPrimary = kTRUE);
  Bool_t PassesCommonEventCuts(AliVEvent *ave);
  Bool_t PassesMixedEventCuts(AliVEvent *ave);
  Bool_t PassesGlobalTrackCuts(AliAODTrack *gtrack); // common cuts for global tracks TBI make it uniform with MC
  Bool_t PassesCommonTrackCuts(AliAODTrack *atrack); // common cuts for analysis specific tracks (e.g. TPC-only) TBI make it uniform with MC
  Bool_t PassesCommonTrackCuts(AliAODMCParticle *amcparticle); // common cuts for analysis specific tracks TBI see above two lines
  virtual void GlobalTracksAOD(AliAODEvent *aAOD, Int_t index); // fill fGlobalTracksAOD in e-b-e . For the meaning of 'index', see declaration of fGlobalTracksAOD
  virtual void GlobalTracksAODTEST(AliAODEvent *aAOD, Int_t index); // fill fGlobalTracksAODTEST in e-b-e . For the meaning of 'index', see declaration of fGlobalTracksAODTEST
  virtual void GlobalTracksAODHA(AliAODEvent *aAOD, Int_t index); // fill fGlobalTracksAODHA in e-b-e . For the meaning of 'index', see declaration of fGlobalTracksAODHA
  virtual void GlobalTracksAOD(AliAODEvent *aAOD, Int_t indexX, Int_t indexY); // fill TExMap *fGlobalTracksAOD1[10][5];
  Double_t RelativeMomenta(AliAODTrack *agtrack1, AliAODTrack *agtrack2);
  Double_t RelativeMomentaComponent(AliAODTrack *agtrack1, AliAODTrack *agtrack2, const char *component);
  Double_t PairVectorComponent(AliAODTrack *agtrack1, AliAODTrack *agtrack2, const char *component);
  Double_t RelativeMomenta(AliAODMCParticle *amcparticle1, AliAODMCParticle *amcparticle2);
  Double_t Q2(AliAODTrack *agtrack1, AliAODTrack *agtrack2);
  Double_t Q3(AliAODTrack *agtrack1, AliAODTrack *agtrack2, AliAODTrack *agtrack3);
  Double_t Q4(AliAODTrack *agtrack1, AliAODTrack *agtrack2, AliAODTrack *agtrack3, AliAODTrack *agtrack4);

  virtual void ResetEBEObjects();
  Bool_t SpecifiedEvent(UInt_t run, UShort_t bunchCross, UInt_t orbit, UInt_t period);
  Int_t CurrentEventNumber();
  virtual void DoSomeDebugging(AliVEvent *ave);
  virtual void CalculateCorrelationFunctions(AliAODEvent *aAOD);
   virtual void Calculate3pCorrelationFunctions(AliAODEvent *aAOD);
   virtual void Calculate4pCorrelationFunctions(AliAODEvent *aAOD);
  virtual void CalculateCorrelationFunctions(AliMCEvent *aMC);
  virtual void CalculateCorrelationFunctionsTEST(AliAODEvent *aAOD);
  virtual void Calculate2pBackground(TClonesArray *ca1, TClonesArray *ca2); // TBI soon will become obsolete
  virtual void Calculate2pBackground(TClonesArray *ca1, TClonesArray *ca2, TExMap *em1, TExMap *em2);
  virtual void Calculate2pBackgroundTEST(TClonesArray *ca1, TClonesArray *ca2, TExMap *em1, TExMap *em2);
  virtual void Calculate3pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3); // TBI soon will become obsolete
  virtual void Calculate3pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TExMap *em1, TExMap *em2, TExMap *em3);
  virtual void Calculate3pBackgroundTEST(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TExMap *em1, TExMap *em2, TExMap *em3);
  virtual void Calculate4pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TClonesArray *ca4);
  virtual void Calculate2pBackground(TClonesArray *ca1, TClonesArray *ca2, Bool_t bMC); // TBI unify with the previous function
  virtual void HybridApproach1stTerm(AliAODEvent *aAOD);
  virtual void HybridApproach2ndTerm(AliAODEvent *aAOD, TClonesArray *ca3, TExMap *em3);
  //virtual void HybridApproach3rdTerm(AliAODEvent *aAOD); // TBI not needed for the time being
  //virtual void HybridApproach4thTerm(AliAODEvent *aAOD); // TBI not needed for the time being
  virtual void HybridApproach5thTerm(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TExMap *em1, TExMap *em2, TExMap *em3);

  // 3.) Methods called in Terminate(Option_t *):
  virtual void GetOutputHistograms(TList *histList);
   // TBI implement the rest as well
   virtual void GetPointersForCorrelationFunctions();
   virtual void GetPointersForBackground();
   virtual void GetPointersForBuffers();
  virtual void NormalizeCorrelationFunctions();
  // 4.) Utility:
  Int_t BinNoForSpecifiedValue(TH1F *hist, Double_t value);
  Int_t BinNoForSpecifiedValue(TProfile *pro, Double_t value);

  // Setters and getters:
  // 0.) Not classified yet;
  // 1.) Control histograms;
  // 2.) Event-by-event histograms;
  // 3.) Correlation functions;
  // 4.) Background;
  // 5.) Buffers;
  // 6.) QA;
  // 7.) ...
  // *.) Debugging

  // 0.) Not classified yet:
  void SetMaxNoGlobalTracksAOD(Int_t mngta) {this->fMaxNoGlobalTracksAOD = mngta;};
  Int_t GetMaxNoGlobalTracksAOD() const {return this->fMaxNoGlobalTracksAOD;};
  void SetProcessBothKineAndReco(Bool_t pbkar) {this->fProcessBothKineAndReco = pbkar;};
  Bool_t GetProcessBothKineAndReco() const {return this->fProcessBothKineAndReco;};
  void SetProcessOnlyKine(Bool_t pok) {this->fProcessOnlyKine = pok;};
  Bool_t GetProcessOnlyKine() const {return this->fProcessOnlyKine;};
  void SetProcessOnlyReco(Bool_t por) {this->fProcessOnlyReco = por;};
  Bool_t GetProcessOnlyReco() const {return this->fProcessOnlyReco;};
  void SetRejectFakeTracks(Bool_t rft) {this->fRejectFakeTracks = rft;};
  Bool_t GetRejectFakeTracks() const {return this->fRejectFakeTracks;};

  // 1.) Control histograms:
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;}
  void SetControlHistogramsFlagsPro(TProfile* const chfp) {this->fControlHistogramsFlagsPro = chfp;};
  TProfile* GetControlHistogramsFlagsPro() const {return this->fControlHistogramsFlagsPro;};
  //void SetFillControlHistograms(Bool_t fch) {this->fFillControlHistograms = fch;}; // TBI remove eventually
  //Bool_t GetFillControlHistograms() const {return this->fFillControlHistograms;}; // TBI remove eventually
  void SetFillControlHistogramsEvent(Bool_t fche) {this->fFillControlHistogramsEvent = fche;};
  Bool_t GetFillControlHistogramsEvent() const {return this->fFillControlHistogramsEvent;};
  void SetFillControlHistogramsNonIdentifiedParticles(Bool_t fchnip) {this->fFillControlHistogramsNonIdentifiedParticles = fchnip;};
  Bool_t GetFillControlHistogramsNonIdentifiedParticles() const {return this->fFillControlHistogramsNonIdentifiedParticles;};
  void SetFillControlHistogramsNonIdentifiedParticlesFTSF(Bool_t fchnipFTSF) {this->fFillControlHistogramsNonIdentifiedParticlesFTSF = fchnipFTSF;};
  Bool_t GetFillControlHistogramsNonIdentifiedParticlesFTSF() const {return this->fFillControlHistogramsNonIdentifiedParticlesFTSF;};
  void SetFilterBitFTSF(Int_t fbFTSF) {this->fFilterBitFTSF = fbFTSF;};
  Int_t GetFilterBitFTSF() const {return this->fFilterBitFTSF;};
  void SetFillControlHistogramsIdentifiedParticles(Bool_t fchip) {this->fFillControlHistogramsIdentifiedParticles = fchip;};
  Bool_t GetFillControlHistogramsIdentifiedParticles() const {return this->fFillControlHistogramsIdentifiedParticles;};
  void SetFillControlHistogramsWithGlobalTrackInfo(Bool_t fchwgti) {this->fFillControlHistogramsWithGlobalTrackInfo = fchwgti;};
  Bool_t GetFillControlHistogramsWithGlobalTrackInfo() const {return this->fFillControlHistogramsWithGlobalTrackInfo;};
  void SetFillControlHistogramsV0s(Bool_t fchv) {this->fFillControlHistogramsV0s = fchv;};
  Bool_t GetFillControlHistogramsV0s() const {return this->fFillControlHistogramsV0s;};
  // 1d) Identified particles:
  void SetInclusiveSigmaCuts(Int_t pidFunction, Double_t sigmaValue)
  {
   // pidFunction: [0=Electron(...),1=Muon(...),2=Pion(...),3=Kaon(...),4=Proton(...)]
   // Example: SetInclusiveSigmaCuts(2,3.); sets in function Pion() inclusive cuts for pions to 3.0 sigmas
   fUseDefaultInclusiveSigmaCuts = kFALSE;
   this->fInclusiveSigmaCuts[pidFunction] = sigmaValue;
  };
  void SetExclusiveSigmaCuts(Int_t pidFunction, Int_t pidExclusive, Double_t sigmaValue)
  {
   // pidFunction: [0=Electron(...),1=Muon(...),2=Pion(...),3=Kaon(...),4=Proton(...)]
   // Example: SetExclusiveSigmaCuts(3,4,4.); sets in function Kaon() exclusive cuts for protons to 4.0 sigmas
   fUseDefaultExclusiveSigmaCuts = kFALSE;
   this->fExclusiveSigmaCuts[pidFunction][pidExclusive] = sigmaValue;
  };

  // 2.) Event-by-event histograms:
  void SetEBEHistogramsList(TList* const ehl) {this->fEBEHistogramsList = ehl;};
  TList* GetEBEHistogramsList() const {return this->fEBEHistogramsList;} 
  void SetEBEObjectsFlagsPro(TProfile* const ehfp) {this->fEBEObjectsFlagsPro = ehfp;};
  TProfile* GetEBEObjectsFlagsPro() const {return this->fEBEObjectsFlagsPro;}; 
  //void SetFillEBEHistograms(Bool_t feh) {this->fFillEBEHistograms = feh;}; // TBI rethink
  //Bool_t GetFillEBEHistograms() const {return this->fFillEBEHistograms;};

  // 3.) Correlation functions:
  void SetCorrelationFunctionsList(TList* const cfl) {this->fCorrelationFunctionsList = cfl;};
  TList* GetCorrelationFunctionsList() const {return this->fCorrelationFunctionsList;}
  void SetCorrelationFunctionsFlagsPro(TProfile* const cffp) {this->fCorrelationFunctionsFlagsPro = cffp;};
  TProfile* GetCorrelationFunctionsFlagsPro() const {return this->fCorrelationFunctionsFlagsPro;};
  void Set2pCorrelationFunctionsFlagsPro(TProfile* const cffp2p) {this->f2pCorrelationFunctionsFlagsPro = cffp2p;};
  TProfile* Get2pCorrelationFunctionsFlagsPro() const {return this->f2pCorrelationFunctionsFlagsPro;};
  void Set3pCorrelationFunctionsFlagsPro(TProfile* const cffp3p) {this->f3pCorrelationFunctionsFlagsPro = cffp3p;};
  TProfile* Get3pCorrelationFunctionsFlagsPro() const {return this->f3pCorrelationFunctionsFlagsPro;};
  void Set4pCorrelationFunctionsFlagsPro(TProfile* const cffp4p) {this->f4pCorrelationFunctionsFlagsPro = cffp4p;};
  TProfile* Get4pCorrelationFunctionsFlagsPro() const {return this->f4pCorrelationFunctionsFlagsPro;};
  void SetFillCorrelationFunctions(Bool_t fcf) {this->fFillCorrelationFunctions = fcf;};
  Bool_t GetFillCorrelationFunctions() const {return this->fFillCorrelationFunctions;};
  void SetNormalizeCorrelationFunctions(Bool_t ncf) {this->fNormalizeCorrelationFunctions = ncf;};
  Bool_t GetNormalizeCorrelationFunctions() const {return this->fNormalizeCorrelationFunctions;};
  void SetFill3pCorrelationFunctions(Bool_t f3pcf) {this->fFill3pCorrelationFunctions = f3pcf;};
  Bool_t GetFill3pCorrelationFunctions() const {return this->fFill3pCorrelationFunctions;};
  void SetFill4pCorrelationFunctions(Bool_t f4pcf) {this->fFill4pCorrelationFunctions = f4pcf;};
  Bool_t GetFill4pCorrelationFunctions() const {return this->fFill4pCorrelationFunctions;};
  void SetNormalizationOption(Int_t fno) {this->fNormalizationOption = fno;};
  Int_t GetNormalizationOption() const {return this->fNormalizationOption;};
  void SetNormalizationInterval(Float_t min, Float_t max)
  {
   this->fNormalizeCorrelationFunctions = kTRUE;
   this->fNormalizationOption = 1;
   this->fNormalizationInterval[0] = min;
   this->fNormalizationInterval[1] = max;
  };
  void SetnMergedBins(Int_t fnmb) {this->fnMergedBins = fnmb;};
  Int_t GetnMergedBins() const {return this->fnMergedBins;};

  // 4.) Background:
  void SetBackgroundList(TList* const bl) {this->fBackgroundList = bl;};
  TList* GetBackgroundList() const {return this->fBackgroundList;}
  void SetBackgroundFlagsPro(TProfile* const bfp) {this->fBackgroundFlagsPro = bfp;};
  TProfile* GetBackgroundFlagsPro() const {return this->fBackgroundFlagsPro;};
  void Set2pBackgroundFlagsPro(TProfile* const bfp2p) {this->f2pBackgroundFlagsPro = bfp2p;};
  TProfile* Get2pBackgroundFlagsPro() const {return this->f2pBackgroundFlagsPro;};
  void Set3pBackgroundFlagsPro(TProfile* const bfp3p) {this->f3pBackgroundFlagsPro = bfp3p;};
  TProfile* Get3pBackgroundFlagsPro() const {return this->f3pBackgroundFlagsPro;};
  void Set4pBackgroundFlagsPro(TProfile* const bfp4p) {this->f4pBackgroundFlagsPro = bfp4p;};
  TProfile* Get4pBackgroundFlagsPro() const {return this->f4pBackgroundFlagsPro;};
  void SetBackgroundOption(Int_t bo) {this->fBackgroundOption = bo;};
  Int_t GetBackgroundOption() const {return this->fBackgroundOption;};
  void SetEstimate2pBackground(Bool_t fe2pb) {this->fEstimate2pBackground = fe2pb;};
  Bool_t GetEstimate2pBackground() const {return this->fEstimate2pBackground;};
  void SetEstimate3pBackground(Bool_t fe3pb) {this->fEstimate3pBackground = fe3pb;};
  Bool_t GetEstimate3pBackground() const {return this->fEstimate3pBackground;};
  void SetEstimate4pBackground(Bool_t fe4pb) {this->fEstimate4pBackground = fe4pb;};
  Bool_t GetEstimate4pBackground() const {return this->fEstimate4pBackground;};
  void SetMaxBufferSize1(Int_t mbs1) {this->fMaxBufferSize1 = mbs1;};
  Int_t GetMaxBufferSize1() const {return this->fMaxBufferSize1;};

  // 5.) Buffers:
  void SetBuffersList(TList* const bl) {this->fBuffersList = bl;};
  TList* GetBuffersList() const {return this->fBuffersList;}
  void SetBuffersFlagsPro(TProfile* const bfp) {this->fBuffersFlagsPro = bfp;};
  TProfile* GetBuffersFlagsPro() const {return this->fBuffersFlagsPro;};
  void SetFillBuffers(Int_t mb) {this->fFillBuffers = kTRUE; this->fMaxBuffer = mb;};

  // 6.) QA:
  void SetQAList(TList* const qal) {this->fQAList = qal;};
  TList* GetQAList() const {return this->fQAList;}
  void SetQAFlagsPro(TProfile* const qafp) {this->fQAFlagsPro = qafp;};
  TProfile* GetQAlagsPro() const {return this->fQAFlagsPro;};
  void SetBailOutAfterQA(Bool_t boaqa) {this->fBailOutAfterQA = boaqa;};
  Bool_t GetBailOutAfterQA() const {return this->fBailOutAfterQA;};
  void SetFillQAEvents(Bool_t ffqae) {this->fFillQAEvents = ffqae;};
  Bool_t GetFillQAEvents() const {return this->fFillQAEvents;};
  void SetFillQAParticles(Bool_t ffqap) {this->fFillQAParticles = ffqap;};
  Bool_t GetFillQAParticles() const {return this->fFillQAParticles;};
  void SetQAEventsList(TList* const qael) {this->fQAEventsList = qael;};
  TList* GetQAEventsList() const {return this->fQAEventsList;}
  void SetQAParticlesList(TList* const qapl) {this->fQAParticlesList = qapl;};
  TList* GetQAParticlesList() const {return this->fQAParticlesList;}

  // 7.) Common event cuts:
  void SetRejectEventsWithoutPrimaryVertex(Bool_t rewpv) {this->fRejectEventsWithoutPrimaryVertex = rewpv;};
  Bool_t GetRejectEventsWithoutPrimaryVertex() const {return this->fRejectEventsWithoutPrimaryVertex;};
  void SetMinMagneticField(Float_t mmf) {this->fMinMagneticField = mmf;};
  Float_t GetMinMagneticField() const {return this->fMinMagneticField;};
  void SetNumberOfTracks(Int_t minnoft, Int_t maxnoft)
  {
   fCutOnNumberOfTracks = kTRUE;
   this->fMinNumberOfTracks = minnoft;
   this->fMaxNumberOfTracks = maxnoft;
  };
  void SetNumberOfGlobalTracks(Int_t minnofgt, Int_t maxnofgt)
  {
   fCutOnNumberOfGlobalTracks = kTRUE;
   this->fMinNumberOfGlobalTracks = minnofgt;
   this->fMaxNumberOfGlobalTracks = maxnofgt;
  };
  void SetNumberOfV0s(Int_t minnofV0s, Int_t maxnofV0s)
  {
   fCutOnNumberOfV0s = kTRUE;
   this->fMinNumberOfV0s = minnofV0s;
   this->fMaxNumberOfV0s = maxnofV0s;
  };
  void SetNumberOfCascades(Int_t minnofc, Int_t maxnofc)
  {
   fCutOnNumberOfCascades = kTRUE;
   this->fMinNumberOfCascades = minnofc;
   this->fMaxNumberOfCascades = maxnofc;
  };
  void SetVertexX(Float_t minvX, Float_t maxvX)
  {
   fCutOnVertexX = kTRUE;
   this->fMinVertexX = minvX;
   this->fMaxVertexX = maxvX;
  };
  void SetVertexY(Float_t minvY, Float_t maxvY)
  {
   fCutOnVertexY = kTRUE;
   this->fMinVertexY = minvY;
   this->fMaxVertexY = maxvY;
  };
  void SetVertexZ(Float_t minvZ, Float_t maxvZ)
  {
   fCutOnVertexZ = kTRUE;
   this->fMinVertexZ = minvZ;
   this->fMaxVertexZ = maxvZ;
  };
  void SetNContributors(Int_t minNc, Int_t maxNc)
  {
   fCutOnNContributors = kTRUE;
   this->fMinNContributors = minNc;
   this->fMaxNContributors = maxNc;
  };

  // 8.) Common global track cuts: // TBI at the moment, they are applied both to 'atracks' and 'gtracks', decouple eventually
  void SetPtRange(Float_t ptMin, Float_t ptMax)
  {
   fApplyGlobalTrackCuts = kTRUE;
   this->fPtRange[0] = ptMin;
   this->fPtRange[1] = ptMax;
  };
  void SetEtaRange(Float_t etaMin, Float_t etaMax)
  {
   fApplyGlobalTrackCuts = kTRUE;
   this->fEtaRange[0] = etaMin;
   this->fEtaRange[1] = etaMax;
  };
  void SetPhiRange(Float_t phiMin, Float_t phiMax)
  {
   fApplyGlobalTrackCuts = kTRUE;
   this->fPhiRange[0] = phiMin;
   this->fPhiRange[1] = phiMax;
  };

  // *.) Testing new ways to calculate correlation functions:
  void SetCorrelationFunctionsTESTList(TList* const cfTl) {this->fCorrelationFunctionsTESTList = cfTl;};
  TList* GetCorrelationFunctionsTESTList() const {return this->fCorrelationFunctionsTESTList;}
  void SetCorrelationFunctionsTESTFlagsPro(TProfile* const cfTfp) {this->fCorrelationFunctionsTESTFlagsPro = cfTfp;};
  TProfile* GetCorrelationFunctionsTESTFlagsPro() const {return this->fCorrelationFunctionsTESTFlagsPro;};
  void SetBoostVelocity(const Double_t vx, const Double_t vy, const Double_t vz)
  {
   fBoost = kTRUE;
   fBoostVelocity = TVector3(vx,vy,vz);
  }
  
  void SetQ2binning(const Int_t nBins, const Double_t min, const Double_t max)
  {
   this->fnQ2bins = nBins;
   this->fnQ2min = min;
   this->fnQ2max = max;
  };
  void SetQ3binning(const Int_t nBins, const Double_t min, const Double_t max)
  {
   this->fnQ3bins = nBins;
   this->fnQ3min = min;
   this->fnQ3max = max;
  };
  void SetFillCorrelationFunctionsTEST(Int_t const testNO, Bool_t bDoTest) {this->fFillCorrelationFunctionsTEST[testNO] = bDoTest;};

  // *.) Testing new ways to calculate background:
  void SetBackgroundTESTList(TList* const bTl) {this->fBackgroundTESTList = bTl;};
  TList* GetBackgroundTESTList() const {return this->fBackgroundTESTList;}
  void SetBackgroundTESTFlagsPro(TProfile* const bTfp) {this->fBackgroundTESTFlagsPro = bTfp;};
  TProfile* GetBackgroundTESTFlagsPro() const {return this->fBackgroundTESTFlagsPro;};
  void SetFillBackgroundTEST(Int_t const testNO, Bool_t bDoTest) {this->fFillBackgroundTEST[testNO] = bDoTest;};

  // *) 'hybrid approach':
  void SetHybridApproachList(TList* const hal) {this->fHybridApproachList = hal;};
  TList* GetHybridApproachList() const {return this->fHybridApproachList;}
  void SetHybridApproachFlagsPro(TProfile* const hafp) {this->fHybridApproachFlagsPro = hafp;};
  TProfile* GetHybridApproachFlagsPro() const {return this->fHybridApproachFlagsPro;};
  void SetDoHybridApproach(Bool_t bhap) {this->fDoHybridApproach = bhap;};

  // *.) Online monitoring:
  void SetUpdateOutputFile(const Int_t uf, const char *uqof)
  {
   // Example usage: taskMPF->SetUpdateOutputFile(44,"AnalysisResults.root");
   this->fOnlineMonitoring = kTRUE;
   this->fUpdateOutputFile = kTRUE;
   this->fUpdateFrequency = uf;
   this->fUpdateWhichOutputFile = new TString(uqof);
  };
  void SetMaxNumberOfEvents(const Int_t mnof, const char *uqof)
  {
   // Example usage: taskMPF->SetMaxNumberOfEvents(44,"AnalysisResults.root");
   this->fOnlineMonitoring = kTRUE;
   this->fMaxNumberOfEvents = mnof;
   this->fUpdateWhichOutputFile = new TString(uqof);
  };

  // *.) Debugging:
  void SetWaitForSpecifiedEvent(UInt_t run, UShort_t bunchCross, UInt_t orbit, UInt_t period)
  {
   this->fDoSomeDebugging = kTRUE;
   this->fWaitForSpecifiedEvent = kTRUE; 
   this->fRun = run;
   this->fBunchCross = bunchCross;
   this->fOrbit = orbit;  
   this->fPeriod = period;
  }; // void SetWaitForSpecifiedEvent(UInt_t run, UShort_t bunchCross, UInt_t orbit, UInt_t period)

 private:
  AliAnalysisTaskMultiparticleFemtoscopy(const AliAnalysisTaskMultiparticleFemtoscopy& aatmpf);
  AliAnalysisTaskMultiparticleFemtoscopy& operator=(const AliAnalysisTaskMultiparticleFemtoscopy& aatmpf);

  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)
  TString *fAnalysisType; //! MC, AOD, ESD, MC_AOD or MC_ESD

  AliPIDResponse *fPIDResponse; //! PID response object

  Int_t fMaxNoGlobalTracksAOD; // maximum # of TExMap *fGlobalTracksAOD objects to be booked. Default is 3, one for default analysis, and two for event mixing
  TExMap *fGlobalTracksAOD[10]; //! global tracks in AOD. [0] is used in the default analysis, [1] and [2] for event mixing, etc.
  Bool_t fProcessBothKineAndReco; // process both aMC and aAOD, or aMC and aESD, typically to get purities (TBI: add support for ESD)
  Bool_t fProcessOnlyKine; // process only aMC
  Bool_t fProcessOnlyReco; // process only aAOD or aESD (i.e. disregard aMC even if available)
  Bool_t fRejectFakeTracks; // if set to kFALSE, and if fMC is available, get the corresponding MC particle by taking absolute value of label
  AliMCEvent *fMC; // placeholder for MC info

  // 1.) Control histograms:
  TList *fControlHistogramsList;        // list to hold all 'control histograms' objects
  TProfile *fControlHistogramsFlagsPro; // profile to hold all flags for control histograms
  Bool_t fFillControlHistograms;        // not set directly, but instead via: fFillControlHistogramsEvent || fFillControlHistogramsNonIdentifiedParticles || ...
  // 1a) Event (a.k.a. global event observables):
  TList *fControlHistogramsEventList;        // list to hold all 'control histograms' for events TBI
  TProfile *fControlHistogramsEventFlagsPro; // profile to hold all flags for control histograms for events TBI
  Bool_t fFillControlHistogramsEvent;        // fill or not control histograms for global event observables
  TH1I *fGetNumberOfTracksHist;              // a{AOD,MC}->GetNumberOfTracks()
  TH1I *fGetNumberOfGlobalTracksHist;        //! fGlobalTracksAOD[0]->GetSize() this is then my multiplicity...
  TH1I *fGetNumberOfV0sHist;                 // aAOD->GetNumberOfV0s()
  TH1I *fGetNumberOfCascadesHist;            // aAOD->GetNumberOfCascades()
  TH1D *fGetMagneticFieldHist;               // aAOD->GetMagneticField()
  TH1I *fGetEventTypeHist;                   // aAOD->GetEventType()
  TH1D *fGetCentralityHist;                  // aAOD->GetCentrality()
  TH1F *fVertexXYZ[3];                       //! [avtx->GetX(),avtx->GetY(),avtx->GetZ()]
  TH1I *fGetNContributorsHist;               // avtx->GetNContributors()
  TH1F *fGetChi2perNDFHist;                  // avtx->GetChi2perNDF();
  TH1I *fGetNDaughtersHist;                  // avtx->GetNDaughters();

  // 1b) Non-identified particles (for AOD these are "normal global" tracks, i.e. the ones which satisfy atrack->GetID()>=0 ):
  TList *fControlHistogramsNonIdentifiedParticlesList;        // list to hold all 'control histograms' for non-identified particles
  TProfile *fControlHistogramsNonIdentifiedParticlesFlagsPro; // profile to hold all flags for control histograms for non-identified particles
  Bool_t fFillControlHistogramsNonIdentifiedParticles;        // fill or not control histograms for non-identified particles
  TH1I *fChargeHist;                                          // atrack->Charge()
  TH1I *fGetTPCNclsHist;                                      // atrack->GetTPCNcls()
  TH1I *fGetTPCsignalNHist;                                   // atrack->GetTPCsignalN()
  TH1I *fGetITSNclsHist;                                      // atrack->GetITSNcls()
  TH2F *fdEdxVsPtHist;                                        // atrack->GetTPCmomentum(),atrack->GetTPCsignal()  
  TH1F *fPtHist;                                              // atrack->Pt()
  TH1F *fEtaHist;                                             // atrack->Eta()
  TH1F *fPhiHist;                                             // atrack->Phi()
  TH1F *fMassHist;                                            // atrack->M()
  TH1I *fGetFilterMap;                                        // atrack->GetFilterMap()
  TH1I *fGetPdgCode;                                          // atrack->GetPdgCode()

  // 1c) Non-identified particles for the specified filterbit (f.t.s.f.) (by default TPC-only):
  TList *fControlHistogramsNonIdentifiedParticlesFTSFList;        // list to hold all 'control histograms' for non-identified particles
  TProfile *fControlHistogramsNonIdentifiedParticlesFTSFFlagsPro; // profile to hold all flags for control histograms for non-identified particles
  Bool_t fFillControlHistogramsNonIdentifiedParticlesFTSF;        // fill or not control histograms for non-identified particles
  Int_t fFilterBitFTSF;                                           // filter bit, relevant only for these group of control histos. For the particle selection, there is another flag
  TH1I *fChargeFTSFHist;                                          // atrack->Charge()
  TH1I *fGetTPCNclsFTSFHist;                                      // atrack->GetTPCNcls()
  TH1I *fGetTPCsignalNFTSFHist;                                   // atrack->GetTPCsignalN()
  TH1I *fGetITSNclsFTSFHist;                                      // atrack->GetITSNcls()
  TH2F *fdEdxVsPtFTSFHist;                                        // atrack->GetTPCmomentum(),atrack->GetTPCsignal()
  TH1F *fPtFTSFHist;                                              // atrack->Pt()
  TH1F *fEtaFTSFHist;                                             // atrack->Eta()
  TH1F *fPhiFTSFHist;                                             // atrack->Phi()
  TH1F *fMassFTSFHist;                                            // atrack->M()
  TH1I *fGetFilterMapFTSF;                                        // atrack->GetFilterMap()
  TH1I *fGetPdgCodeFTSF;                                          // atrack->GetPdgCode()

  // 1d) Identified particles:
  TList *fControlHistogramsIdentifiedParticlesList;        // list to hold all 'control histograms' for identified particles
  TProfile *fControlHistogramsIdentifiedParticlesFlagsPro; // profile to hold all flags for control histograms for identified particles
  Bool_t fFillControlHistogramsIdentifiedParticles;        // fill or not control histograms for identified particles (by default they are not filled)
  Bool_t fFillControlHistogramsWithGlobalTrackInfo;        // by default, control histograms are filled with info from 'atrack'. If this flag is TRUE, then instead info from 'gtrack' is used. This then also applies to info used to get correlation functions and background as well
  TProfile *fInclusiveSigmaCutsPro;                        //! holds the values of fInclusiveSigmaCuts[5];
  TProfile2D *fExclusiveSigmaCutsPro;                      //! holds the values of fExclusiveSigmaCuts[5][5];
  TH1F *fMassPIDHist[5][2][2];                             //! [0=e,1=mu,2=pi,3=K,4=p][particle(+q)/antiparticle(-q)][kPrimary/kFromDecayVtx]
  TH1F *fPtPIDHist[5][2][2];                               //! [0=e,1=mu,2=pi,3=K,4=p][particle(+q)/antiparticle(-q)][kPrimary/kFromDecayVtx]
  TH1F *fPPIDHist[5][2][2][3];                             //! [0=e,1=mu,2=pi,3=K,4=p][particle(+q)/antiparticle(-q)][kPrimary/kFromDecayVtx][xyz]
  TH1F *fEtaPIDHist[5][2][2];                              //! [0=e,1=mu,2=pi,3=K,4=p][particle(+q)/antiparticle(-q)][kPrimary/kFromDecayVtx]
  TH1F *fPhiPIDHist[5][2][2];                              //! [0=e,1=mu,2=pi,3=K,4=p][particle(+q)/antiparticle(-q)][kPrimary/kFromDecayVtx]
  Bool_t fUseDefaultInclusiveSigmaCuts;                    // if the setter SetInclusiveSigmaCuts(...) (see above) is not called explicitly, the default hardwired values will be used
  Bool_t fUseDefaultExclusiveSigmaCuts;                    // if the setter SetExclusiveSigmaCuts(...) (see above) is not called explicitly, the default hardwired values will be used
  Double_t fInclusiveSigmaCuts[5];                         // [PID function] see .cxx for detailed documentation
  Double_t fExclusiveSigmaCuts[5][5];                      // [PID function][PID exclusive] see .cxx for detailed documentation

  // ...
  // 1e) V0s:
  TList *fControlHistogramsV0sList;        // list to hold all 'control histograms' for V0s
  TProfile *fControlHistogramsV0sFlagsPro; // profile to hold all flags for control histograms for V0s
  Bool_t fFillControlHistogramsV0s;        // fill or not control histograms for V0s (by default they are not filled)
  TH1I *fGetNProngsHist;                   // aAODv0->GetNProngs()
  TH1F *fMassK0ShortHist;                  // aAODv0->MassK0Short()
  TH1F *fMassLambdaHist;                   // aAODv0->MassLambda()
  TH1F *fMassAntiLambdaHist;               // aAODv0->MassAntiLambda()
  TH1F *fOpenAngleV0Hist;                  // aAODv0->OpenAngleV0() // same as aAODv0->ProngsRelAngle()
  TH1F *fRadiusV0Hist;                     // aAODv0->RadiusV0()
  TH1F *fDcaV0ToPrimVertexHist;            // aAODv0->DcaV0ToPrimVertex()
  TH1F *fMomV0XHist;                       // aAODv0->MomV0X()
  TH1F *fMomV0YHist;                       // aAODv0->MomV0Y()
  TH1F *fMomV0ZHist;                       // aAODv0->MomV0Z()
  TH1F *fPtV0Hist;                         // pow(aAODv0->Pt2V0(),0.5)
  TH1F *fPseudoRapV0Hist;                  // aAODv0->PseudoRapV0()
  TH2F *fPAHist;                           // Pod.-Arm.

  // ...
  // 1d) Cascades:
  // ...

  // 2.) Event-by-event objects:
  TList *fEBEHistogramsList;        // list to hold all stuff from e-b-e histograms
  TProfile *fEBEObjectsFlagsPro;    // profile to hold all flags for e-b-e histograms for V0s
  //Bool_t fFillEBEHistograms;      // fill or not e-b-e histograms TBI do I really need this?
  TH1I *fUniqueIDHistEBE;           // filled with aAODv0->GetPosID() and aAODv0->GetNegID(). If the bin corresponding to that ID is already filled, two V0s share the same daughter     
  TClonesArray *fPIDCA[5][2][2];    //! holds AliAODTrack candidates for each event [0=e,1=mu,2=pi,3=K,4=p][particle(+q)/antiparticle(-q)][kPrimary/kFromDecayVtx]
  TClonesArray *fPIDV0sCA[1];       //! holds AliAODv0 candidates for each event [0=Lambda,1=...]

  // 3.) Correlation functions:
  TList *fCorrelationFunctionsList;              // list to hold all correlation functions for primary particle
  TProfile *fCorrelationFunctionsFlagsPro;       // profile to hold all flags for correlation functions
  TProfile *f2pCorrelationFunctionsFlagsPro;     //! profile to hold all flags for 2p correlation functions (placed in fCorrelationFunctionsSublist[0])
  TProfile *f3pCorrelationFunctionsFlagsPro;     //! profile to hold all flags for 3p correlation functions (placed in fCorrelationFunctionsSublist[1])
  TProfile *f4pCorrelationFunctionsFlagsPro;     //! profile to hold all flags for 4p correlation functions (placed in fCorrelationFunctionsSublist[2])
  TList *fCorrelationFunctionsSublist[3];        //! lists to hold all correlation functions, for 2p [0], 3p [1], 4p [2], etc., separately
  Bool_t fFillCorrelationFunctions;              // fill or not correlation functions (by default they are not filled)
  Bool_t fNormalizeCorrelationFunctions;         // normalize correlation functions with the background
  TExMap *fCorrelationFunctionsIndices;          //! associates pdg code to index of correlation function
  TH1F *fCorrelationFunctions[10][10];           //! [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 5=e,6=mu,7=pi,8=K,9=p] x [same]. Booking only upper 1/2 of the matrix, diagonal included.
  Bool_t fFill3pCorrelationFunctions;            // fill 3-p correlation functions
  TH1F *f3pCorrelationFunctions[10][10][10];     //! [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 5=e,6=mu,7=pi,8=K,9=p] x [same] x [same].
  Bool_t fFill4pCorrelationFunctions;            // fill 4-p correlation functions
  TH1F *f4pCorrelationFunctions[10][10][10][10]; //! [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 5=e,6=mu,7=pi,8=K,9=p] x [same] x [same] x [same].
  Int_t fNormalizationOption;                    // set here how to normalize the correlation function: 0 = "just scale", 1 = "use concrete interval", 2 = ...
  Float_t fNormalizationInterval[2];             // concrete example: 0.15 < q < 0.175 GeV/c. Then, fNormalizationInterval[0] is the low edge, etc. See the relevant setter SetNormalizationInterval
  Int_t fnMergedBins;                            // before normalization, both signal and background will be rebinned with this value

  // 4.) Background:
  TList *fBackgroundList;              // list to hold all background objects primary particle
  TProfile *fBackgroundFlagsPro;       // profile to hold all flags for background
  TProfile *f2pBackgroundFlagsPro;     //! profile to hold all flags for 2p background (placed in fBackgroundSublist[0])
  TProfile *f3pBackgroundFlagsPro;     //! profile to hold all flags for 3p background (placed in fBackgroundSublist[1])
  TProfile *f4pBackgroundFlagsPro;     //! profile to hold all flags for 4p background (placed in fBackgroundSublist[2])
  TList *fBackgroundSublist[3];        //! lists to hold all background correlations, for 2p [0], 3p [1], 4p [2], etc., separately
  Int_t fBackgroundOption;             // set how to estimate background: 0 = "shifting", 1 = "permutations", etc. (see .cxx for further explanation). By default, it is "shifting"
  Bool_t fEstimate2pBackground;        // enable or not 2p background estimation
  Bool_t fEstimate3pBackground;        // enable or not 3p background estimation
  Bool_t fEstimate4pBackground;        // enable or not 4p background estimation
  TH1F *f2pBackground[10][10];         //! [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p] x [same]. Booking only upper 1/2 of the matrix, diagonal included.
  TH1F *f3pBackground[10][10][10];     //! [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p] x [same] x [same]
  TH1F *f4pBackground[10][10][10][10]; //! [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p] x [same] x [same]
  TClonesArray *fMixedEvents0[3];      //! tracks for mixed events (supporting up to 3-mixed events at the moment). Used only for fBackgroundOption = 0. Global tracks are in TExMap *fGlobalTracksAOD[10]; above, TBI make it uniform eventually, i.e. decouple global tracks for the default analysis from the background
  Int_t fMaxBufferSize1;                // the second index in fMixedEvents1[10][50]; and fGlobalTracksAOD1[10][50]; is booked only up to this number. When this number is reached, calculation is done, and buffer is cleaned. max = 50. defaulted to 10
  TClonesArray *fMixedEvents1[10][50];  //! tracks for mixed events. 10 vertex z-ranges. Keep at maximum 5 events in the buffer. Used only for fBackgroundOption = 1
  TExMap *fGlobalTracksAOD1[10][50];    //! global tracks in AOD. Used only for fBackgroundOption = 1. Indices must be the same as in TClonesArray *fMixedEvents1[10][5];

  // 5.) Buffers:
  TList *fBuffersList;                             // list to hold all objects for buffers
  TProfile *fBuffersFlagsPro;                      // profile to hold all flags for buffers
  Bool_t fFillBuffers;                             // hold some thingies for bunch of events in memories
  Int_t fMaxBuffer;                                // max buffer size (e.g. for 3-p correlations it is 3, etc.) TBI there is a problem apparently, re-think
  TClonesArray *fChargedParticlesCA[2][10][10000]; //! [0=AOD||ESD,1=MC][#events,max=10][particles]
  TExMap *fChargedParticlesEM[10];                 //! [#events,max=10,has to correspond to 2nd entry above] this is standard mapping, nothing more nor less than that...

  // 6.) QA:
  TList *fQAList;                   // list to holds all QA objects. It is nested in: a) "QA events"; b) "QA particles"; c) ...
  TProfile *fQAFlagsPro;            // list to holds all flags for QA objects
  Bool_t fFillQA;                   // fill QA objects. Not set directly, but instead via: fFillQA = fFillQAEvents || fFillQAParticles || ...
  Bool_t fBailOutAfterQA;           // fill QA objects and bail out, i.e. do not do the actual analysis
  Bool_t fFillQAEvents;             // fill two sets of histograms, before and after event cuts
  Bool_t fFillQAParticles;          // fill two sets of histograms, before and after particle cuts
  TList *fQAEventsList;             // list to holds all objects for "QA events"
  TList *fQAParticlesList;          // list to holds all objects for "QA particles"
  TH1I *fQAFilterBitScan;           // for each track in AOD, dump it's filterbits
  TH2I *fQAIDvsFilterBit;           // atrack->ID() vs. filterbit
  TH1F *fQAParticleHist[2][10][10]; //! [0="before rain",1="after rain"][distribution_index][cut_index]

  // 7.) Common event cuts (TBI validated only for AOD analysis, for the time being...):
  //  a) Cuts on AliAODEvent:
  Bool_t fRejectEventsWithoutPrimaryVertex; // as the name says it, by default set to kTRUE. has a setter
  Float_t fMinMagneticField;                // defaulted to 0.001, compared to aAOD->GetMagneticField()
  Bool_t fCutOnNumberOfTracks;              // cut on the total number of tracks in AOD, i.e. on aAOD->GetNumberOfTracks(). This is NOT multiplicity, since some tracks are stored multiple times in AOD
  Int_t fMinNumberOfTracks;                 // default values never in effect; if aAOD->GetNumberOfTracks() < fMinNumberOfTracks, event is rejected
  Int_t fMaxNumberOfTracks;                 // default values never in effect; if aAOD->GetNumberOfTracks() > fMaxNumberOfTracks, event is rejected
  Bool_t fCutOnNumberOfGlobalTracks;        // cut on the total number of 'normal' global tracks in AOD, i.e. on fGlobalTracksAOD[0]->GetSize()
  Int_t fMinNumberOfGlobalTracks;           // default values never in effect; if fGlobalTracksAOD[0]->GetSize() < fMinNumberOfGlobalTracks, event is rejected
  Int_t fMaxNumberOfGlobalTracks;           // default values never in effect; if fGlobalTracksAOD[0]->GetSize() > fMaxNumberOfGlobalTracks, event is rejected
  Bool_t fCutOnNumberOfV0s;                 // cut on the total number of V0s in AOD, i.e. on aAOD->GetNumberOfV0s()
  Int_t fMinNumberOfV0s;                    // default values never in effect; if aAOD->GetNumberOfV0s() < fMinNumberOfV0s, event is rejected
  Int_t fMaxNumberOfV0s;                    // default values never in effect; if aAOD->GetNumberOfV0s() > fMaxNumberOfV0s, event is rejected
  Bool_t fCutOnNumberOfCascades;            // cut on the total number of cascades in AOD, i.e. on aAOD->GetNumberOfCascades()
  Int_t fMinNumberOfCascades;               // default values never in effect; if aAOD->GetNumberOfCascades() < fMinNumberOfCascades, event is rejected
  Int_t fMaxNumberOfCascades;               // default values never in effect; if aAOD->GetNumberOfCascades() > fMaxNumberOfCascades, event is rejected

  //  b) Cuts on AliAODVertex:
  Bool_t fCutOnVertexX;       // cut on the x position of vertex, i.e. on avtx->GetX()
  Float_t fMinVertexX;        // default values never in effect; if avtx->GetX() < fMinVertexX, event is rejected
  Float_t fMaxVertexX;        // default values never in effect; if avtx->GetX() > fMaxVertexX, event is rejected
  Bool_t fCutOnVertexY;       // cut on the y position of vertex, i.e. on avtx->GetY()
  Float_t fMinVertexY;        // default values never in effect; if avtx->GetY() < fMinVertexY, event is rejected
  Float_t fMaxVertexY;        // default values never in effect; if avtx->GetY() > fMaxVertexY, event is rejected
  Bool_t fCutOnVertexZ;       // cut on the z position of vertex, i.e. on avtx->GetZ()
  Float_t fMinVertexZ;        // default values never in effect; if avtx->GetZ() < fMinVertexZ, event is rejected
  Float_t fMaxVertexZ;        // default values never in effect; if avtx->GetZ() > fMaxVertexZ, event is rejected
  Bool_t fCutOnNContributors; // cut on avtx->GetNContributors()
  Int_t fMinNContributors;    // default values never in effect; if avtx->GetNContributors() < fMinNContributors, event is rejected
  Int_t fMaxNContributors;    // default values never in effect; if avtx->GetNContributors() > fMaxNContributors, event is rejected

  // 8.) Common global track cuts (applied only on "normal" global tracks in AOD):
  TList *fGlobalTrackCutsList;        // list to hold all objects for common global track cuts
  TProfile *fGlobalTrackCutsFlagsPro; // profile to hold all flags
  Bool_t fApplyGlobalTrackCuts;       // if set to kFALSE, the default hardwired cuts will be used. TBI doesn't do anything at the moment in .cxx
  Float_t fPtRange[2];                // ptMin = fPtRange[0], ptMax = fPtRange[1]
  Float_t fEtaRange[2];               // etaMin = etaRange[0], etaMax = etaRange[1]
  Float_t fPhiRange[2];               // phiMin = phiRange[0], phiMax = phiRange[1]

  // *.) Testing new ways to calculate correlation functions and cumulants:
  TList *fCorrelationFunctionsTESTList;              // list to hold all TEST correlation functions for primary particle
  TProfile *fCorrelationFunctionsTESTFlagsPro;       // profile to hold all flags for TEST correlation functions
  Bool_t fBoost;                                     // boost or not   
  TVector3 fBoostVelocity;                           // boost everything in the system mocing with relative velocity fBoostVelocity    
  Int_t fnQ2bins;                                    // number of bins for all histos and profiles vs. Q2 (both for signal and background)
  Double_t fnQ2min;                                  // min bin for all histos and profiles vs. Q2 (both for signal and background)
  Double_t fnQ2max;                                  // max bin for all histos and profiles vs. Q2 (both for signal and background)
  Int_t fnQ3bins;                                    // number of bins for all histos and profiles vs. Q3 (both for signal and background)
  Double_t fnQ3min;                                  // min bin for all histos and profiles vs. Q3 (both for signal and background)
  Double_t fnQ3max;                                  // max bin for all histos and profiles vs. Q3 (both for signal and background)
  TList *fCorrelationFunctionsTESTSublist[10];       //! lists to hold all TEST correlation functions, they are enumerated in .cxx file
  Bool_t fFillCorrelationFunctionsTEST[10];          // fill or not particular TEST correlation functions, they are enumerated in .cxx file (by default all are set to FALSE)
  TProfile *fCorrelationFunctionsTEST[10][2][7][10]; //! [testNo][0=vs Q2, 1=vs Q3][example [0=<x1>][1=<x2>], ...,[6=<x1x2x3>]][differential index, e.g. for test 0 [0=Cx][1=Cy][2=Cz]]
  TProfile *fSignalCumulantsTEST[10][2][4][10];      //! [testNo][0=vs Q2, 1=vs Q3][[0=<x1x2>_c][1=<x1x3>_c][2=<x2x3>_c][3=<x1x2x3>_c]][differential index, e.g. for test 0 [0=Cx][1=Cy][2=Cz]]
  TH1F *fSignalYieldTEST[2];                         //! [0=for <X1X2> and Q2, 1=for <X1X2X3> and Q3]
  TH1F *fEab_TEST6[2];                               //! [0=for signal "Test 6", 1=for background "Test 6"]

  // *.) Testing new ways to calculate background functions:
  TList *fBackgroundTESTList;                       // list to hold all TEST background for primary particle
  TProfile *fBackgroundTESTFlagsPro;                // profile to hold all flags for TEST background
  TList *fBackgroundTESTSublist[10];                //! lists to hold all TEST background, they are enumerated in .cxx file
  Bool_t fFillBackgroundTEST[10];                   // fill or not particular TEST background, they are enumerated in .cxx file (by default all are set to FALSE)
  TProfile *fBackgroundTEST[10][2][7][10];          //! [testNo][0=vs Q2, 1=vs Q3][correlation][differential index, e.g. for test 0 [0=Cx][1=Cy][2=Cz]]
  TProfile *fBackgroundCumulantsTEST[10][2][4][10]; //! [testNo][0=vs Q2, 1=vs Q3][[0=<x1x2>_c][1=<x1x3>_c][2=<x2x3>_c][3=<x1x2x3>_c]][differential index, e.g. for test 0 [0=Cx][1=Cy][2=Cz]]
  TClonesArray *fMixedEventsTEST[3];                //! tracks for mixed events, using just 'shifting' for simplicity TBI make it it more sophisticated later
  TExMap *fGlobalTracksAODTEST[3];                  //! global tracks in AOD
  TH1F *fBackgroundYieldTEST[2];                    //! [0=for <X1X2> and Q2, 1=for <X1X2X3> and Q3]

  // *.) 'hybrid approach':
  TList *fHybridApproachList;        // list to hold all histos for 'hybrid approach' a la UH
  TProfile *fHybridApproachFlagsPro; // profile to hold all flags for 'hybrid approach' a la UH
  Bool_t fDoHybridApproach;          // do or not the correlations via the 'hybrid approach' a la UH
  TH1F *fDistributionsVsQ3[5];       //! five distinct distributions in the numerator of Eq. (18)
  TClonesArray *fMixedEventsHA[3];   //! tracks for mixed events, using just 'shifting' for simplicity TBI make it it more sophisticated later
  TExMap *fGlobalTracksAODHA[3];     //! global tracks in AOD

  // *.) Online monitoring:
  Bool_t fOnlineMonitoring;        // enable online monitoring (not set excplicitly!), the flags below just refine it
  Bool_t fUpdateOutputFile;        // update the output file after certain number of analysed events
  Int_t fUpdateFrequency;          // after how many events the output file will be updated
  TString *fUpdateWhichOutputFile; // which file will be regularly updated
  Int_t fMaxNumberOfEvents;        // if this number of events is reached, write to external file and bail out

  // *.) Debugging:
  Bool_t fDoSomeDebugging;        // enable call to function within which debugging is done. Set indirectly.
  Bool_t fWaitForSpecifiedEvent;  // do something only for the specified event
  UInt_t fRun;                    // do something only for the specified event
  UShort_t fBunchCross;           // do something only for the specified event
  UInt_t fOrbit;                  // do something only for the specified event
  UInt_t fPeriod;                 // do something only for the specified event

  ClassDef(AliAnalysisTaskMultiparticleFemtoscopy,20);

};

//================================================================================================================

#endif











