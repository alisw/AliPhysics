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
  virtual void InitializeArraysForControlHistograms(); 
  virtual void InitializeArraysForEBEObjects(); 

  // 1.) Methods called in UserCreateOutputObjects():
  virtual void BookAndNestAllLists();
  virtual void BookEverythingForControlHistograms();
  virtual void BookEverythingForEBEObjects();
  // 2.) Methods called in UserExec(Option_t *):
  virtual void InsanityChecksUserExec();
  Bool_t Pion(AliAODTrack *atrack, Int_t charge = 1, Bool_t bPrimary = kTRUE);
  Bool_t Kaon(AliAODTrack *atrack, Int_t charge = 1, Bool_t bPrimary = kTRUE);
  Bool_t Proton(AliAODTrack *atrack, Int_t charge = 1, Bool_t bPrimary = kTRUE);
  Bool_t PassesCommonEventCuts(AliAODEvent *aAOD);
  Bool_t PassesCommonTrackCuts(AliAODTrack *gtrack); // TBI I am applying them to global tracks... rethink
  virtual void GlobalTracksAOD(AliAODEvent *aAOD); // fill fGlobalTracksAOD in e-b-e
  virtual void ResetEBEObjects();
  Bool_t SpecifiedEvent(UInt_t run, UShort_t bunchCross, UInt_t orbit, UInt_t period);
  // 3.) Methods called in Terminate(Option_t *):
  //...


  // Setters and getters:
  // 1.) Control histograms;
  // 2.) Event-by-event histograms;
  // 3.) ... 
  // *.) Debugging

  // 1.) Control histograms:
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;} 
  void SetControlHistogramsFlagsPro(TProfile* const chfp) {this->fControlHistogramsFlagsPro = chfp;};
  TProfile* GetControlHistogramsFlagsPro() const {return this->fControlHistogramsFlagsPro;}; 
  void SetFillControlHistograms(Bool_t fch) {this->fFillControlHistograms = fch;};
  Bool_t GetFillControlHistograms() const {return this->fFillControlHistograms;};
  void SetFillControlHistogramsEvent(Bool_t fche) {this->fFillControlHistogramsEvent = fche;};
  Bool_t GetFillControlHistogramsEvent() const {return this->fFillControlHistogramsEvent;};
  void SetFillControlHistogramsNonIdentifiedParticles(Bool_t fchnip) {this->fFillControlHistogramsNonIdentifiedParticles = fchnip;};
  Bool_t GetFillControlHistogramsNonIdentifiedParticles() const {return this->fFillControlHistogramsNonIdentifiedParticles;};
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
  // 3.) TBI
  // *.) Debugging:
  void SetWaitForSpecifiedEvent(UInt_t run, UShort_t bunchCross, UInt_t orbit, UInt_t period)
  {
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

  AliPIDResponse *fPIDResponse; //! PID response object

  TExMap *fGlobalTracksAOD; //! global tracks in AOD  
  

  // 1.) Control histograms:  
  TList *fControlHistogramsList;        // list to hold all 'control histograms' objects
  TProfile *fControlHistogramsFlagsPro; // profile to hold all flags for control histograms
  Bool_t fFillControlHistograms;        // fill or not control histograms (by default they are all filled)
  // 1a) Event:  
  TList *fControlHistogramsEventList;        // list to hold all 'control histograms' for events TBI
  TProfile *fControlHistogramsEventFlagsPro; // profile to hold all flags for control histograms for events TBI
  Bool_t fFillControlHistogramsEvent;        // fill or not control histograms events TBI (by default they are not filled)
  TH1I *fGetNumberOfTracksHist;              // aAOD->GetNumberOfTracks()
  TH1I *fGetNumberOfV0sHist;                 // aAOD->GetNumberOfV0s()
  TH1F *fVertexXYZ[3];                       // [avtx->GetX(),avtx->GetY(),avtx->GetZ()]
  TH1I *fGetNContributorsHist;               // avtx->GetNContributors() 
  TH1F *fGetChi2perNDFHist;                  // avtx->GetChi2perNDF();
  TH1I *fGetNDaughtersHist;                  // avtx->GetNDaughters();
  // ...
  // 1b) Non-identified particles:
  TList *fControlHistogramsNonIdentifiedParticlesList;        // list to hold all 'control histograms' for non-identified particles
  TProfile *fControlHistogramsNonIdentifiedParticlesFlagsPro; // profile to hold all flags for control histograms for non-identified particles
  Bool_t fFillControlHistogramsNonIdentifiedParticles;        // fill or not control histograms for non-identified particles (by default they are not filled)
  TH1I *fChargeHist;                                          // atrack->Charge()
  TH1I *fGetTPCNclsHist;                                      // atrack->GetTPCNcls()
  TH1I *fGetTPCsignalNHist;                                   // atrack->GetTPCsignalN()
  TH1I *fGetITSNclsHist;                                      // atrack->GetITSNcls()
  TH2F *fdEdxVsPtHist;                                        // atrack->GetTPCmomentum(),atrack->GetTPCsignal()  
  TH1F *fPtHist;                                              // atrack->Pt()
  TH1F *fEtaHist;                                             // atrack->Eta()
  TH1F *fPhiHist;                                             // atrack->Phi()
  TH1F *fMassHist;                                            // atrack->M()

  // 1c) Identified particles:
  TList *fControlHistogramsIdentifiedParticlesList;        // list to hold all 'control histograms' for identified particles
  TProfile *fControlHistogramsIdentifiedParticlesFlagsPro; // profile to hold all flags for control histograms for identified particles
  Bool_t fFillControlHistogramsIdentifiedParticles;        // fill or not control histograms for identified particles (by default they are not filled)
  TH1F *fMassPIDHist[5][2][2];                             // [0=e,1=mu,2=pi,3=K,4=p][particle/antiparticle][kPrimary/kFromDecayVtx]   
  TH1F *fPtPIDHist[5][2][2];                               // [0=e,1=mu,2=pi,3=K,4=p][particle/antiparticle][kPrimary/kFromDecayVtx]   
  TH1F *fEtaPIDHist[5][2][2];                              // [0=e,1=mu,2=pi,3=K,4=p][particle/antiparticle][kPrimary/kFromDecayVtx]   
  TH1F *fPhiPIDHist[5][2][2];                              // [0=e,1=mu,2=pi,3=K,4=p][particle/antiparticle][kPrimary/kFromDecayVtx]   

  // ...
  // 1c) V0s:
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
  TClonesArray *fPIDCA[5][2][2];    //! holds AliAODTrack candidates for each event [0=e,1=mu,2=pi,3=K,4=p][particle/antiparticle][kPrimary/kFromDecayVtx]
  TClonesArray *fPIDV0sCA[1];       //! holds AliAODv0 candidates for each event [0=Lambda,1=...]

  // Internal flags:
  Bool_t fUseInternalFlags;  // use internal flags (automatically set if some internal flag is used)

  // *. Debugging:
  Bool_t fWaitForSpecifiedEvent; //! do something only for the specified event
  UInt_t fRun;                    //! do something only for the specified event
  UShort_t fBunchCross;           //! do something only for the specified event
  UInt_t fOrbit;                  //! do something only for the specified event
  UInt_t fPeriod;                 //! do something only for the specified event

  // Control histograms:
  //Bool_t fFillControlHistograms;     // fill or not control histograms (by default they are filled)

  ClassDef(AliAnalysisTaskMultiparticleFemtoscopy,1); 

};

//================================================================================================================

#endif











