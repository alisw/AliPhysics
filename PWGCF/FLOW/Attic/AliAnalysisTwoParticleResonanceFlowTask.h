/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliAnalysisTwoParticleResonanceFlowTask:
// origin: Redmer Alexander Bertens (rbertens@nikhef.nl)
// analyis task for Resonance-meson reconstruction and estimation of v_n

#ifndef ALIANALYSISTWOPARTICLERESONANCEFLOWTASK_H
#define ALIANALYSISTWOPARTICLERESONANCEFLOWTASK_H

class TH1F;
class TH2F;
class TProfile;
class AliFlowTrackCuts;
class AliFlowEvent;
class AliFlowTrackSimple;
class AliFlowEventSimple;
class TDirectoryFile;
class AliFlowCandidateTrack;
class AliFlowBayesianPID;
class AliEventPoolManager;
class AliResonanceFlowHelperTrack;

#include "AliAnalysisTaskSE.h"
class AliResonanceFlowHelperTrack : public TObject
{
public:
        AliResonanceFlowHelperTrack(Double_t eta, Double_t phi, Double_t p, Double_t px, Double_t py, Double_t pz, Double_t pt, Int_t charge, Double_t mass, Int_t id, Int_t species) : fEta(eta), fPhi(phi), fp(p), fpX(px), fpY(py), fpZ(pz), fpT(pt), fCharge(charge), fMass(mass), fID(id), fSpecies(species) {  }
    ~AliResonanceFlowHelperTrack() {}
    virtual Double_t P()                const { return fp; }
    virtual Double_t Px()               const { return fpX; }
    virtual Double_t Py()               const { return fpY; }
    virtual Double_t Pz()               const { return fpZ; }
    virtual Double_t Pt()               const { return fpT; }
    virtual Double_t Phi()              const { return fPhi; }
    virtual Double_t Eta()              const { return fEta; }
    virtual Int_t Charge()              const { return fCharge; }
    virtual Double_t Mass()             const { return fMass; }
    virtual Short_t ID()                const { return fID; }
    virtual Int_t Species()             const { return fSpecies; }
    void    InitializeHelperTrack(Double_t eta, Double_t phi, Double_t p, Double_t px, Double_t py, Double_t pz, Double_t pt, Int_t charge, Double_t mass, Int_t id, Int_t species) { fEta = eta; fPhi = phi; fp = p; fpX = px; fpY = py; fpZ = pz; fpT = pt; fCharge = charge; fMass = mass; fID = id; fSpecies = species; }
private:
    Double_t                             fEta;      // eta
    Double_t                             fPhi;      // phi
    Double_t                             fp;        // p
    Double_t                             fpX;       // pX
    Double_t                             fpY;       // pY
    Double_t                             fpZ;       // pZ
    Double_t                             fpT;       // pT
    Int_t                                fCharge;   // charge
    Double_t                             fMass;     // mass
    Short_t                              fID;       // id
    Int_t                                fSpecies;  // species
    ClassDef(AliResonanceFlowHelperTrack, 1); // lightweight helper track for phi reconstruction
};

class AliAnalysisTwoParticleResonanceFlowTask : public AliAnalysisTaskSE
{
public:
   AliAnalysisTwoParticleResonanceFlowTask();
   AliAnalysisTwoParticleResonanceFlowTask(const char *name);
   virtual ~AliAnalysisTwoParticleResonanceFlowTask();
   // technical aspects of the analysis
   void                                 ForceExit(Int_t type, const char* message);
   Bool_t                               SetIsMC(Bool_t ismc) {fIsMC = ismc; return fIsMC; }
   void                                 IsMC();
   Bool_t                               UseEventMixing(Bool_t mix) { fEventMixing = mix; return mix; }
   Bool_t                               UsePhiMinusPsiMethod(Bool_t p) {fPhiMinusPsiMethod = p; return p;}
   Bool_t                               SetVZEROSubEvents(Bool_t v0) { fV0 = v0; return v0; }
   // configure the output of the analysis
   TH1F*                                BookHistogram(const char * name);
   TH2F*                                BookPIDHistogram(const char * name, Bool_t TPC);
   TH1F*                                InitPtSpectraHistograms(Float_t nmin, Float_t nmax);
   TH1F*                                BookPtHistogram(const char* name);
   Bool_t                               InitializeAnalysis();
   //PK void                                 AddResonanceIdentificationOutputObjects();
   virtual void                         UserCreateOutputObjects();
   // setters
   void                                 SetPtBins(Float_t bin[19], Int_t n) { for(Int_t i = 0; i < n+1; i++) fPtBins[i] = bin[i]; fNPtBins = n; }
   void                                 SetdPhiBins(Float_t bin[19], Int_t n) { for(Int_t i = 0; i < n+1; i++) fdPhiBins[i] = bin[i]; fNdPhiBins = n;}
   void                                 SetCentralityParameters(Double_t min, Double_t max, const char* a, const char* b, Bool_t c, Bool_t d) { 
                                                                                          fCentralityMin = min; 
                                                                                          fCentralityMax = max; 
                                                                                          fkCentralityMethodA = a; 
                                                                                          fkCentralityMethodB = b;
                                                                                          fCentralityCut2010 = c; 
											  fCentralityCut2011 = d; }
   void                                 SetPOICuts(AliFlowTrackCuts *cutsPOI) { fPOICuts = cutsPOI; }
   void                                 SetRPCuts(AliFlowTrackCuts *cutsRP) { fCutsRP = cutsRP; }
   void                                 SetPIDConfiguration(Float_t prob[7]) { for(Int_t i = 0; i < 7; i++) fPIDConfig[i] = prob[i]; }
   Bool_t                               SetQA(Bool_t qa) {fQA = qa; return fQA;}
   void                                 SetAddTaskMacroSummary(Float_t m[12]) {for(Int_t i(0); i < 12; i++) fAddTaskMacroSummary[i] = m[i];}
   void                                 SetPOIDCAXYZ(Float_t dca[5]) { for(Int_t i = 0; i < 5; i++) fDCAConfig[i] = dca[i]; }
   void                                 SetMixingBins(Int_t c[20], Int_t v[20]) {for(Int_t i = 0; i < 20; i++) { fCentralityMixingBins[i] = c[i];
                                                                                                                 fVertexMixingBins[i] = v[i]; } }
   void                                 SetMixingParameters(Int_t p[3]) { for(Int_t i = 0; i < 3; i++) fMixingParameters[i] = p[i]; }
   void                                 SetupSpeciesA(Int_t species, Int_t charge, Float_t mass, Float_t minPtA, Float_t maxPtA) {fSpeciesA = species; fChargeA = charge; fMassA = mass; fMinPtA = minPtA; fMaxPtA = maxPtA;}
   void                                 SetupSpeciesB(Int_t species, Int_t charge, Float_t mass, Float_t minPtB, Float_t maxPtB) {fSpeciesB = species; fChargeB = charge; fMassB = mass; fMinPtB = minPtB; fMaxPtB = maxPtB;}   
   void                                 SetCandidateEtaAndPt(Float_t mineta, Float_t maxeta, Float_t minpt, Float_t maxpt) { fCandidateMinEta = mineta;
                                                                                                                                 fCandidateMaxEta = maxeta;
                                                                                                                                 fCandidateMinPt = minpt;
                                                                                                                                 fCandidateMaxPt = maxpt;
                                                                                                                                 fCandidateEtaPtCut = kTRUE;}
   void                                 SetCommonConstants(Int_t massBins, Float_t minMass, Float_t maxMass) { fMassBins = massBins;
                                                                                                                 fMinMass = minMass;
                                                                                                                 fMaxMass = maxMass; }
   void                                 SetVertexZ(Float_t z) { fVertexRange = z; }
   void                                 SetMaxDeltaDipAngleAndPt(Float_t a, Float_t pt) { fDeltaDipAngle = a;
                                                                                          fDeltaDipPt = pt;
                                                                                          fApplyDeltaDipCut = kTRUE; };
   //getters
   void                                 GetMixingParameters(Int_t p[3]) const { for(Int_t i = 0; i < 3; i++) p[i] = fMixingParameters[i]; } 
   Float_t                              GetCenMin() const {return fCentralityMin; }
   Float_t                              GetCenMax() const {return fCentralityMax; }
   const char*                          GetCentralityMethod() const {return fkCentralityMethodA; }
   Float_t                              GetVertexZ() const { return fVertexRange; }
   Float_t                              GetDeltaDipAngle() const {return fDeltaDipAngle; }
   Float_t                              GetDeltaDipPt() const {return fDeltaDipPt; }
   void                                 GetPIDConfiguration(Float_t prob[7]) const {for(Int_t i = 0; i < 7; i++) prob[i] = fPIDConfig[i]; }
   void                                 GetPOIDCZXYZ(Float_t dca[5]) const { for(Int_t i = 0; i < 5; i++) dca[i] = fDCAConfig[i]; }
   void                                 GetCandidateEtaAndPt(Float_t etapt[4]) const { etapt[0] = fCandidateMinEta;
                                                                                        etapt[1] = fCandidateMaxEta;
                                                                                        etapt[2] = fCandidateMinPt;
                                                                                        etapt[3] = fCandidateMaxPt; }
   AliFlowEvent*                        GetFlowEvent() const {return fFlowEvent;}
   // the analysis itself
   AliEventPoolManager*                 InitializeEventMixing();
   template <typename T> Float_t        InvariantMass(const T* track1, const T* track2) const;
   template <typename T> Float_t        DeltaDipAngle(const T* track1, const T* track2) const;
   template <typename T> Bool_t         CheckDeltaDipAngle(const T* track1, const T* track2) const;
   template <typename T> Bool_t         CheckCandidateEtaPtCut(const T* track1, const T* track2) const;
   template <typename T> Bool_t         EventCut(T* event);
   template <typename T> void           PlotMultiplcities(const T* event) const;
   template <typename T> Bool_t         CheckVertex(const T* event);
   template <typename T> Bool_t         CheckCentrality(T* event);
   void                                 InitializeBayesianPID(AliAODEvent* event);
   template <typename T> Bool_t         PassesTPCbayesianCut(T* track, Int_t species) const;
   Bool_t                               PassesDCACut(AliAODTrack* track) const;
   Bool_t                               DoOwnPID(AliAODTrack* track, Int_t species) const;
   Bool_t                               AcceptTrack(AliAODTrack* track, Int_t species) const;
   template <typename T> Float_t        PairPt(const T* track_1, const T* track_2, Bool_t phi = kFALSE) const;
   template <typename T> Float_t        PtSelector(Int_t _track_type, const T* track_1, const T* track_2, Float_t mass) const;
   template <typename T> Bool_t         QualityCheck(T* track) const;
   void                                 TrackQA(AliAODTrack* track, Int_t species, Bool_t allChargedParticles) const;
   template <typename T> void           SetNullCuts(T* esd);
   void                                 PrepareFlowEvent(Int_t iMulti);
   void                                 PhiMinusPsiMethod(TObjArray* MixingCandidates);
   void                                 PhiMinusPsiMethodWriteData(Bool_t signal, TObjArray* SpeciesA, TObjArray* SpeciesB, Float_t* abcPsi2);
   void                                 VZEROSubEventAnalysis();
   void                                 DoAnalysisOnTheFly(AliFlowEventSimple* event);
   void                                 DoAnalysisOnTheFly(TObjArray* MixingCandidates, TObjArray* SpeciesA, TObjArray* ocSpeciesA, TObjArray* SpeciesB, TObjArray* ocSpeciesB); 
   void                                 DoAnalysisOnTheFly(TDirectoryFile* outputFile);
   virtual void                         UserExec(Option_t *option);
   void                                 ResonanceSignal(TObjArray* SpeciesA, TObjArray* SpeciesB) const;
   void                                 ResonanceBackground(TObjArray* SpeciesA, TObjArray* SpeciesB, Bool_t checkAutoCorrelations = kTRUE) const;
   void                                 ReconstructionWithEventMixing(TObjArray* MixingCandidates) const;
   virtual void                         Terminate(Option_t *);

private:

   Int_t                fSpeciesA; // particle species a
   Int_t                fSpeciesB; // species b
   Int_t                fChargeA; // charge for species a
   Int_t                fChargeB; // charge for species b
   Float_t              fMassA; // mass for species a
   Float_t              fMassB; // mass  for species b
   Float_t              fMinPtA; // min pt for species a
   Float_t              fMaxPtA; // max pt for species a
   Float_t              fMinPtB; // min pt for species b
   Float_t              fMaxPtB; // max pt for species b
   Bool_t               fIsMC; // use mc mode
   Bool_t               fEventMixing; // use event mixing
   Bool_t               fPhiMinusPsiMethod; //  use phi minus psi method (default is invariant mass fit method)
   Bool_t               fQA; // make qa plots
   Bool_t               fV0; // use three subevents including vzero
   Int_t                fMassBins; // mass bins
   Float_t              fMinMass; // mass range
   Float_t              fMaxMass; // mass range
   AliFlowTrackCuts     *fCutsRP; // track cuts for reference particles
   AliFlowTrackCuts     *fNullCuts; // dummy cuts for flow event tracks
   AliPIDResponse       *fPIDResponse; //! pid response object
   AliFlowEvent         *fFlowEvent; //! flow events (one for each inv mass band)
   AliFlowBayesianPID   *fBayesianResponse; //!PID response object
   TObjArray            *fCandidates; // candidate array
   Bool_t               fCandidateEtaPtCut; // set eta and pt cut for candidate tracks and combinatorial background
   Float_t              fCandidateMinEta; // minimum eta for candidates
   Float_t              fCandidateMaxEta; // maximum eta for candidates
   Float_t              fCandidateMinPt; // minimum pt for candidates
   Float_t              fCandidateMaxPt; // maximum pt for candidates
   Float_t              fPIDConfig[7]; // configure pid routine
   Float_t              fDCAConfig[5]; // configure dca routine
   Int_t                fMixingParameters[3]; // mixing: poolsize, mixing tracks, pool buffer
   Int_t                fCentralityMixingBins[20]; // configure centrality bins for event mixing
   Int_t                fVertexMixingBins[20]; // configure vertex bins for event mixing
   Float_t              fPtBins[19]; // pt bin borders
   Float_t              fdPhiBins[19]; // dPhi bin borders
   Int_t                fNPtBins; // no of pt bins + 1
   Int_t                fNdPhiBins; // no of dphi bins + 1
   Float_t              fCentrality; // event centrality
   Float_t              fVertex; // event vertex z 
   AliAODEvent          *fAOD;    //! AOD oject
   AliEventPoolManager  *fPoolManager; //! event pool manager
   TList                *fOutputList; // ! Output list
   TH1F                 *fEventStats; // ! Histogram for event statistics
   TH1F                 *fCentralityPass; // ! QA histogram of events that pass centrality cut
   TH1F                 *fCentralityNoPass; //! QA histogram of events that do not pass centrality cut
   TH2F                 *fNOPID;//! QA histogram of TPC response of all charged particles
   TH2F                 *fPIDk;//! QA histogram of TPC response of kaons
   TH2F                 *fPIDp; //! QA histogram of TPC response of pions
   TH1F                 *fResonanceSignal[18]; //! signal histograms
   TH1F                 *fResonanceBackground[18]; //! like-sign kaon pairs
   TH1F                 *fPtSpectra[18]; //! pt spectra
   TH1F                 *fPtP; //! QA histogram of p_t distribution of positive particles
   TH1F                 *fPtN; //! QA histogram of p_t distribution of negative particles
   TH1F                 *fPtSpeciesA; //! QA histogram of p_t distribution of species A
   TH1F                 *fPtSpeciesB; //! QA histogram of p_t distribution of species B
   TH2F                 *fMultCorAfterCuts; //! QA profile global and tpc multiplicity after outlier cut
   TH2F                 *fMultvsCentr; //! QA profile of centralty vs multiplicity
   Float_t              fCentralityMin; // lower bound of cenrality bin
   Float_t              fCentralityMax; // upper bound of centrality bin
   const char           *fkCentralityMethodA; // centrality determiantion (primary method)
   const char           *fkCentralityMethodB; // centrality determination fallback
   Bool_t               fCentralityCut2010; // 3 sigma cut for multiplicity outliers 
   Bool_t               fCentralityCut2011; // 3 sigma cut for multiplicity outliers 
   AliFlowTrackCuts     *fPOICuts; // cuts for particles of interest (flow package)
   Float_t              fVertexRange; // absolute value of maximum distance of vertex along the z-axis
   TH1F                 *fPhi; //! QA plot of azimuthal distribution of POI daughters
   TH1F                 *fEta; //! QA plot of eta distribution of POI daughters
   TH1F                 *fVZEROA; //! QA plot vzeroa multiplicity (all tracks in event)
   TH1F                 *fVZEROC; //! QA plot vzeroc multiplicity (all tracks in event)
   TH1F                 *fTPCM; //! QA plot TPC multiplicity (tracks used for event plane estimation)
   Float_t              fDeltaDipAngle; // absolute value of delta dip angle to be excluded
   Float_t              fDeltaDipPt; // upper value of pt range in which delta dip angle must be applied
   Bool_t               fApplyDeltaDipCut; // enforce delta dip cut
   TH2F                 *fDCAAll;//! qa dca of all charged particles
   TH1F                 *fDCAXYQA; //! qa plot of dca xz
   TH1F                 *fDCAZQA; //!qa plot of dca z
   TH2F                 *fDCAPrim; //!dca of primaries (mc) or kaons (data)
   TH2F                 *fDCASecondaryWeak; //! dca of weak (mc)
   TH2F                 *fDCAMaterial; //!dca material (mc) all (data)
   TProfile             *fSubEventDPhiv2; //! subevent resolution info for v2
   TProfile             *fV0Data[18][2]; //! profiles for vzero vn(minv)
   TH1F                 *fPhiMinusPsiDataContainer[18][18][2]; //! histograms for phi minus psi results
   TH1F                 *fPhiMinusPsiBackgroundContainer[18][18][2]; //! histograms for phi minus psi background
   TH1F                 *fAnalysisSummary; //! plot analysis flags
   Float_t              fAddTaskMacroSummary[12]; // add task macro summary

   AliAnalysisTwoParticleResonanceFlowTask(const AliAnalysisTwoParticleResonanceFlowTask&); // Not implemented
   AliAnalysisTwoParticleResonanceFlowTask& operator=(const AliAnalysisTwoParticleResonanceFlowTask&); // Not implemented
   void                 MakeTrack(Float_t, Float_t, Float_t, Float_t, Int_t , Int_t[]) const;

   ClassDef(AliAnalysisTwoParticleResonanceFlowTask, 4);
};

#endif



