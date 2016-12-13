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

  // 1.) Methods called in UserCreateOutputObjects():
  virtual void InsanityChecksUserCreateOutputObjects();
  virtual void BookAndNestAllLists();
  virtual void BookEverything(); // use this to book all un-classified objects
  virtual void BookEverythingForControlHistograms();
  virtual void BookEverythingForEBEObjects();
  virtual void BookEverythingForCorrelationFunctions();
  virtual void BookEverythingForBackground();
  virtual void BookEverythingForBuffers();
  virtual void BookEverythingForQA();

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
  Bool_t PassesCommonGlobalTrackCuts(AliAODTrack *gtrack); // common cuts for global tracks TBI make it uniform with MC
  Bool_t PassesCommonTrackCuts(AliAODTrack *atrack); // common cuts for analysis specific tracks (e.g. TPC-only) TBI make it uniform with MC
  Bool_t PassesCommonTrackCuts(AliAODMCParticle *amcparticle); // common cuts for analysis specific tracks TBI see above two lines
  virtual void GlobalTracksAOD(AliAODEvent *aAOD, Int_t index); // fill fGlobalTracksAOD in e-b-e . For the meaning of 'index', see declaration of fGlobalTracksAOD
  Double_t RelativeMomenta(AliAODTrack *agtrack1, AliAODTrack *agtrack2);
  Double_t RelativeMomenta(AliAODMCParticle *amcparticle1, AliAODMCParticle *amcparticle2);
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
  virtual void CalculateBackground(TClonesArray *ca1, TClonesArray *ca2);
   virtual void Calculate3pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3);
   virtual void Calculate4pBackground(TClonesArray *ca1, TClonesArray *ca2, TClonesArray *ca3, TClonesArray *ca4);
  virtual void CalculateBackground(TClonesArray *ca1, TClonesArray *ca2, Bool_t bMC); // TBI unify with the previous function
  // 3.) Methods called in Terminate(Option_t *):
  virtual void GetOutputHistograms(TList *histList);
   // TBI implement the rest as well
   virtual void GetPointersForCorrelationFunctions();
   virtual void GetPointersForBackground();
   virtual void GetPointersForBuffers();
  virtual void NormalizeCorrelationFunctions();

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
  void SetFillControlHistogramsV0s(Bool_t fchv) {this->fFillControlHistogramsV0s = fchv;};
  Bool_t GetFillControlHistogramsV0s() const {return this->fFillControlHistogramsV0s;};
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
  void SetFillCorrelationFunctions(Bool_t fcf) {this->fFillCorrelationFunctions = fcf;};
  Bool_t GetFillCorrelationFunctions() const {return this->fFillCorrelationFunctions;};
  void SetNormalizeCorrelationFunctions(Bool_t ncf) {this->fNormalizeCorrelationFunctions = ncf;};
  Bool_t GetNormalizeCorrelationFunctions() const {return this->fNormalizeCorrelationFunctions;};
  void SetFill3pCorrelationFunctions(Bool_t f3pcf) {this->fFill3pCorrelationFunctions = f3pcf;};
  Bool_t GetFill3pCorrelationFunctions() const {return this->fFill3pCorrelationFunctions;};
  void SetFill4pCorrelationFunctions(Bool_t f4pcf) {this->fFill4pCorrelationFunctions = f4pcf;};
  Bool_t GetFill4pCorrelationFunctions() const {return this->fFill4pCorrelationFunctions;};

  // 4.) Background:
  void SetBackgroundList(TList* const bl) {this->fBackgroundList = bl;};
  TList* GetBackgroundList() const {return this->fBackgroundList;}
  void SetBackgroundFlagsPro(TProfile* const bfp) {this->fBackgroundFlagsPro = bfp;};
  TProfile* GetBackgroundFlagsPro() const {return this->fBackgroundFlagsPro;};
  void SetEstimate2pBackground(Bool_t fe2pb) {this->fEstimate2pBackground = fe2pb;};
  Bool_t GetEstimate2pBackground() const {return this->fEstimate2pBackground;};
  void SetEstimate3pBackground(Bool_t fe3pb) {this->fEstimate3pBackground = fe3pb;};
  Bool_t GetEstimate3pBackground() const {return this->fEstimate3pBackground;};
  void SetEstimate4pBackground(Bool_t fe4pb) {this->fEstimate4pBackground = fe4pb;};
  Bool_t GetEstimate4pBackground() const {return this->fEstimate4pBackground;};

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
  TH1I *fGetNumberOfGlobalTracksHist;        // fGlobalTracksAOD[0]->GetSize() this is then my multiplicity...
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
  TH1F *fMassPIDHist[5][2][2];                             //! [0=e,1=mu,2=pi,3=K,4=p][particle(+q)/antiparticle(-q)][kPrimary/kFromDecayVtx]
  TH1F *fPtPIDHist[5][2][2];                               //! [0=e,1=mu,2=pi,3=K,4=p][particle(+q)/antiparticle(-q)][kPrimary/kFromDecayVtx]
  TH1F *fEtaPIDHist[5][2][2];                              //! [0=e,1=mu,2=pi,3=K,4=p][particle(+q)/antiparticle(-q)][kPrimary/kFromDecayVtx]
  TH1F *fPhiPIDHist[5][2][2];                              //! [0=e,1=mu,2=pi,3=K,4=p][particle(+q)/antiparticle(-q)][kPrimary/kFromDecayVtx]
  Double_t fInclusiveSigmaCuts[5];                         //! [PID function] see .cxx for detailed documentation
  Double_t fExclusiveSigmaCuts[5][5];                      //! [PID function][PID exclusive] see .cxx for detailed documentation

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
  TProfile *fEBEObjectsFlagsPro; // profile to hold all flags for e-b-e histograms for V0s 
  //Bool_t fFillEBEHistograms;        // fill or not e-b-e histograms TBI do I really need this?
  TH1I *fUniqueIDHistEBE;           // filled with aAODv0->GetPosID() and aAODv0->GetNegID(). If the bin corresponding to that ID is already filled, two V0s share the same daughter     
  TClonesArray *fPIDCA[5][2][2];    //! holds AliAODTrack candidates for each event [0=e,1=mu,2=pi,3=K,4=p][particle(+q)/antiparticle(-q)][kPrimary/kFromDecayVtx]
  TClonesArray *fPIDV0sCA[1];       //! holds AliAODv0 candidates for each event [0=Lambda,1=...]

  // 3.) Correlation functions:
  TList *fCorrelationFunctionsList;              // list to hold all correlation functions for primary particle
  TProfile *fCorrelationFunctionsFlagsPro;       // profile to hold all flags for correlation functions
  Bool_t fFillCorrelationFunctions;              // fill or not correlation functions (by default they are not filled)
  Bool_t fNormalizeCorrelationFunctions;         // normalize correlation functions with the background
  TExMap *fCorrelationFunctionsIndices;          // associates pdg code to index of correlation function
  TH1F *fCorrelationFunctions[10][10];           //! [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 5=e,6=mu,7=pi,8=K,9=p] x [same]. Booking only upper 1/2 of the matrix, diagonal included.
  Bool_t fFill3pCorrelationFunctions;            // fill 3-p correlation functions
  TH1F *f3pCorrelationFunctions[10][10][10];     //! [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 5=e,6=mu,7=pi,8=K,9=p] x [same] x [same].
  Bool_t fFill4pCorrelationFunctions;            // fill 4-p correlation functions
  TH1F *f4pCorrelationFunctions[10][10][10][10]; //! [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 5=e,6=mu,7=pi,8=K,9=p] x [same] x [same] x [same].

  // 4.) Background:
  TList *fBackgroundList;              // list to hold all correlation functions for primary particle
  TProfile *fBackgroundFlagsPro;       // profile to hold all flags for correlation functions
  Bool_t fEstimate2pBackground;        // enable or not 2p background estimation
  Bool_t fEstimate3pBackground;        // enable or not 3p background estimation
  Bool_t fEstimate4pBackground;        // enable or not 4p background estimation
  TH1F *fBackground[10][10];           //! [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p] x [same]. Booking only upper 1/2 of the matrix, diagonal included.
  TH1F *f3pBackground[10][10][10];     //! [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p] x [same] x [same]
  TH1F *f4pBackground[10][10][10][10]; //! [particle(+q): 0=e,1=mu,2=pi,3=K,4=p, anti-particle(-q): 0=e,1=mu,2=pi,3=K,4=p] x [same] x [same]
  TClonesArray *fMixedEvents[4];       //! tracks for mixed events (supporting up to 4-mixed events at the moment)

  // 5.) Buffers:
  TList *fBuffersList;                             // list to hold all objects for buffers
  TProfile *fBuffersFlagsPro;                      // profile to hold all flags for buffers
  Bool_t fFillBuffers;                             // hold some thingies for bunch of events in memories
  Int_t fMaxBuffer;                                // max buffer size (e.g. for 3-p correlations it is 3, etc.)
  TClonesArray *fChargedParticlesCA[2][10][10000]; // [0=AOD||ESD,1=MC][#events,max=10][particles]
  TExMap *fChargedParticlesEM[10];                 // [#events,max=10,has to correspond to 2nd entry above] this is standard mapping, nothing more nor less than that...

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
  TH1F *fQAParticleHist[2][10][10]; // [0="before rain",1="after rain"][distribution_index][cut_index]

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

  // *.) Online monitoring:
  Bool_t fOnlineMonitoring;        // enable online monitoring (not set excplicitly!), the flags below just refine it
  Bool_t fUpdateOutputFile;        // update the output file after certain number of analysed events
  Int_t fUpdateFrequency;          // after how many events the output file will be updated
  TString *fUpdateWhichOutputFile; // which file will be regularly updated
  Int_t fMaxNumberOfEvents;        // if this number of events is reached, write to external file and bail out

  // *.) Debugging:
  Bool_t fDoSomeDebugging;        //! enable call to function within which debugging is done. Set indirectly.
  Bool_t fWaitForSpecifiedEvent;  //! do something only for the specified event
  UInt_t fRun;                    //! do something only for the specified event
  UShort_t fBunchCross;           //! do something only for the specified event
  UInt_t fOrbit;                  //! do something only for the specified event
  UInt_t fPeriod;                 //! do something only for the specified event

  ClassDef(AliAnalysisTaskMultiparticleFemtoscopy,8);

};

//================================================================================================================

#endif











