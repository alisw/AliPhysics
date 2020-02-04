/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKUNIFLOWMULTISTRANGE_H
#define ALIANALYSISTASKUNIFLOWMULTISTRANGE_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "THnSparse.h"

class TString;
class TComplex;
class TFile;
class AliDirList;
class TList;
class TClonesArray;
class TProfile;
class TH1D;
class TH2D;
class TH3D;

class AliPIDResponse;
class AliPIDCombined;
class AliVEvent;
class AliAODEvent;
class AliVTrack;
class AliAODTrack;
class AliPicoTrack;
class AliAODv0;
class AliAODcascade;
class AliAODMCParticle;

//_____________________________________________________________________________

class AliAnalysisTaskUniFlowMultiStrange : public AliAnalysisTaskSE
{
    public:
      enum    RunMode {kFull = 0, kTest, kSkipFlow}; // task running mode (NOT GRID MODE)
      enum    ColSystem {kPP = 0, kPPb, kPbPb}; // tag for collisional system
      enum    AnalType {kAOD = 0, kESD}; // tag for analysis type
      enum    CentEst {kRFP = 0, kV0A, kV0C, kV0M, kCL0, kCL1, kZNA, kZNC}; // multiplicity/centrality estimator as AliMultSelection
      enum    PartSpecies {kRefs = 0, kCharged, kPion, kKaon, kProton, kK0s, kLambda, kPhi,kXi, kOmega,kUnknown}; // list of all particle species of interest; NB: kUknown last as counter
      enum    SparseCand {kInvMass = 0, kCent, kPt, kEta, kSample,kDim}; // reconstructed candidates dist. dimensions

      class CorrTask
      {
        public:
                      CorrTask(); // default ctor
                      CorrTask(Bool_t refs, Bool_t pois, std::vector<Int_t> harm, std::vector<Double_t> gaps = std::vector<Double_t>()); // actual ctor
                      ~CorrTask() { fiHarm.clear(); fdGaps.clear(); }

          Bool_t      HasGap() const { return (Bool_t) fiNumGaps; }; // check if Gap
          void        Print() const; // print CorrTask properties

          Bool_t                fbDoRefs; // which particles are procesed (RFPs / POIs / both )
          Bool_t                fbDoPOIs; // which particles are procesed (RFPs / POIs / both )
          Int_t                 fiNumHarm; // correlation order <M>
          Int_t                 fiNumGaps; // number of subevents
          std::vector<Int_t>    fiHarm; // harmonics n1,n2,...,nM
          std::vector<Double_t> fdGaps; // gaps between subevents (standard GF notation)
          TString               fsName; // automatically generated name: see Init() for format
          TString               fsLabel; // automatically generated label see Init() for format
        protected:
        private:
      };

                              AliAnalysisTaskUniFlowMultiStrange(); // constructor
                              AliAnalysisTaskUniFlowMultiStrange(const char *name); // named (primary) constructor
                              AliAnalysisTaskUniFlowMultiStrange(const AliAnalysisTaskUniFlowMultiStrange&); // not implemented
                              AliAnalysisTaskUniFlowMultiStrange& operator=(const AliAnalysisTaskUniFlowMultiStrange&); // not implemented
      virtual                 ~AliAnalysisTaskUniFlowMultiStrange(); // destructor


      virtual void            UserCreateOutputObjects(); //
      virtual void            UserExec(Option_t* option); // main methond - called for each event
      virtual void            Terminate(Option_t* option); // called after all events are processed
      // analysis setters
      void                    SetRunMode(RunMode mode = kFull) { fRunMode = mode; }
      void                    Set2018Data(Bool_t use = kTRUE){Is2018Data = use;}
      void                    SetAdditional2018DataEventCut(){IsAdditional2018DataEventCut = kTRUE;}
      void                    SetPIDCorrection(const char* file){IsPIDorrection = kTRUE; fPIDCorrectionPath = file;}
      void                    SetNumEventsAnalyse(Int_t num) { fNumEventsAnalyse = num; }
      void                    SetDumpTObjectTable(Bool_t dump = kTRUE) { fDumpTObjectTable = dump; }
      void					          SetAnalysisType(AnalType type = kAOD) { fAnalType = type; }
      void                    SetMC(Bool_t mc = kTRUE) { fMC = mc; }
      void                    SetSampling(Bool_t sample = kTRUE, Int_t iNum = 10) { fSampling = sample; fNumSamples = iNum; }
      void                    SetFillQAhistos(Bool_t fill = kTRUE) { fFillQA = fill; }
      void                    SetProcessPID(Bool_t use = kTRUE) { fProcessSpec[kPion] = use; fProcessSpec[kKaon] = use; fProcessSpec[kProton] = use; }
      void                    SetProcessV0s(Bool_t use = kTRUE) { fProcessSpec[kK0s] = use; fProcessSpec[kLambda] = use; }
      void                    SetProcessPhi(Bool_t use = kTRUE) { fProcessSpec[kPhi] = use; }


      void       SetProcessCascades(Bool_t use = kTRUE) { fProcessSpec[kXi] = use; fProcessSpec[kOmega] = use; }




      // flow related setters
      void                    AddTwo(Int_t n1, Int_t n2, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { fVecCorrTask.push_back(new CorrTask(refs, pois, {n1,n2})); }
      void                    AddTwoGap(Int_t n1, Int_t n2, Double_t gap, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { fVecCorrTask.push_back(new CorrTask(refs, pois, {n1,n2}, {gap})); }
      void                    AddThree(Int_t n1, Int_t n2, Int_t n3, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { fVecCorrTask.push_back(new CorrTask(refs, pois, {n1,n2,n3})); }
      void                    AddThreeGap(Int_t n1, Int_t n2, Int_t n3, Double_t gap, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { fVecCorrTask.push_back(new CorrTask(refs, pois, {n1,n2,n3} ,{gap})); }
      void                    AddFour(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { fVecCorrTask.push_back(new CorrTask(refs, pois, {n1,n2,n3,n4})); }
      void                    AddFourGap(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Double_t gap, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { fVecCorrTask.push_back(new CorrTask(refs, pois, {n1,n2,n3,n4}, {gap})); }

      void                    SetFlowRFPsPt(Double_t min, Double_t max) { fFlowRFPsPtMin = min; fFlowRFPsPtMax = max; }
      void                    SetFlowPOIsPt(Double_t min, Double_t max, Int_t bins = 0) { fFlowPOIsPtMin = min; fFlowPOIsPtMax = max; fFlowPOIsPtBinNum = bins; }
      void                    SetFlowEta(Double_t max, Int_t bins = 0) { fFlowEtaMax = max; fFlowEtaBinNum = bins; }
      void                    SetFlowPhiBins(Int_t bins) { fFlowPhiBinNum = bins; }
      void                    SetFlowFillWeights(Bool_t weights = kTRUE) { fFlowFillWeights = weights; }
      void                    SetUseWeigthsFile(const char* file, Bool_t bRunByRun) { fFlowWeightsPath = file; fFlowRunByRunWeights = bRunByRun; fFlowUseWeights = kTRUE; } //! NOTE file has to include "alien:///" if the file is on grid
      void                    SetUseWeights3D(Bool_t use = kTRUE) { fFlowUse3Dweights = use; }
      // events setters
      void                    SetCollisionSystem(ColSystem colSystem = kPP) { fColSystem = colSystem; }
      void                    SetCentrality(CentEst est, Int_t min = 0, Int_t max = 0, Int_t bins = 0) { fCentEstimator = est; fCentMin = min; fCentMax = max; fCentBinNum = bins; }
      void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
      void					          SetPVtxZMax(Double_t z) { fPVtxCutZ = z; }
      void                    SetRejectAddPileUp(Bool_t use = kTRUE) { fEventRejectAddPileUp = use; }
      // track setters
      void                    SetChargedDCAzMax(Double_t dcaz) {  fCutChargedDCAzMax = dcaz; }
      void                    SetChargedDCAxyMax(Double_t dcaxy) {  fCutChargedDCAxyMax = dcaxy; }
      void                    SetChargedNumTPCclsMin(UShort_t tpcCls) { fCutChargedNumTPCclsMin = tpcCls; }
      void                    SetChargedTrackFilterBit(UInt_t filter) { fCutChargedTrackFilterBit = filter; }
      // PID (pi,K,p) setters
      void                    SetPIDUseAntiProtonOnly(Bool_t use = kTRUE) { fCutPIDUseAntiProtonOnly = use; }
      void                    SetPIDNumSigmasPionMax(Float_t numSigmas) { fCutPIDnSigmaMax[kPion] = numSigmas; }
      void                    SetPIDNumSigmasKaonMax(Float_t numSigmas) { fCutPIDnSigmaMax[kKaon] = numSigmas; }
      void                    SetPIDNumSigmasProtonMax(Float_t numSigmas) { fCutPIDnSigmaMax[kProton] = numSigmas; }
      void                    SetPIDNumSigmasTPCRejectElectron(Float_t numSigmas) { fCutPIDnSigmaTPCRejectElectron = numSigmas; }
      void                    SetPIDNumSigmasCombinedNoTOFrejection(Bool_t reject = kTRUE) { fCutPIDnSigmaCombinedTOFrejection = reject; }
      void                    SetUseBayesPID(Bool_t bayes = kTRUE) { fCutUseBayesPID = bayes; }
      void                    SetPIDBayesProbPionMin(Double_t probPi) { fCutPIDBayesMin[kPion] = probPi; }
      void                    SetPIDBayesProbKaonMin(Double_t probK) { fCutPIDBayesMin[kKaon] = probK; }
      void                    SetPIDBayesProbProtonMin(Double_t probP) { fCutPIDBayesMin[kProton] = probP; }
      // V0s setters
      void		      SetV0sOnFly(Bool_t onFly) { fCutV0sOnFly = onFly; }
      void		      SetV0sTPCRefit(Bool_t refit) { fCutV0srefitTPC = refit; }
      void		      SetV0sRejectKinks(Bool_t reject) { fCutV0srejectKinks = reject; }
      void                    SetV0sDaughterNumTPCClsMin(UShort_t cls) { fCutV0sDaughterNumTPCClsMin = cls; }
      void                    SetV0sDaughterNumTPCrossMin(UShort_t cls) { fCutV0sDaughterNumTPCCrossMin = cls; }
      void                    SetV0sDaughterNumTPCFindMin(UShort_t cls) { fCutV0sDaughterNumTPCFindMin = cls; }
      void                    SetV0sDaughterNumTPCClsPIDMin(UShort_t cls) { fCutV0sDaughterNumTPCClsPIDMin = cls; }
      void                    SetV0sDaughterRatioCrossFindMin(Double_t ratio) { fCutV0sDaughterRatioCrossFindMin = ratio; }
      void					          SetV0sUseCrossMassRejection(Bool_t reject) { fCutV0sCrossMassRejection = reject; }
      void					          SetV0sCrossMassCutK0s(Double_t mass) { fCutV0sCrossMassCutK0s = mass; }
      void					          SetV0sCrossMassCutLambda(Double_t mass) { fCutV0sCrossMassCutLambda = mass; }
      void					          SetV0sDCAPVMin(Double_t dca) { fCutV0sDCAtoPVMin = dca; }
      void					          SetV0sDCAPVMax(Double_t dca) { fCutV0sDCAtoPVMax = dca; }
      void					          SetV0sDCAPVzMax(Double_t dca) { fCutV0sDCAtoPVzMax = dca; }
      void                    SetV0sDaughtersFilterBit(UInt_t filter) { fCutV0sDaughterFilterBit = filter; }
      void					          SetV0sDCADaughtersMin(Double_t dca) { fCutV0sDCADaughtersMin = dca; }
      void					          SetV0sDCADaughtersMax(Double_t dca) { fCutV0sDCADaughtersMax = dca; }
      void					          SetV0sDecayRadiusMin(Double_t radius) { fCutV0sDecayRadiusMin = radius; }
      void					          SetV0sDecayRadiusMax(Double_t radius) { fCutV0sDecayRadiusMax = radius; }
      void					          SetV0sDaughterEtaMax(Double_t eta) { fCutV0sDaughterEtaMax = eta; }
      void					          SetV0sDaughterPtMin(Double_t pt) { fCutV0sDaughterPtMin = pt; }
      void					          SetV0sDaughterPtMax(Double_t pt) { fCutV0sDaughterPtMax = pt; }
      void                    SetV0sMotherRapMax(Double_t rap) { fCutV0sMotherRapMax = rap; }
      void					          SetV0sK0sInvMassMin(Double_t mass) { fCutV0sInvMassK0sMin = mass; }
      void					          SetV0sK0sInvMassMax(Double_t mass) { fCutV0sInvMassK0sMax = mass; }
      void					          SetV0sLambdaInvMassMin(Double_t mass) { fCutV0sInvMassLambdaMin = mass; }
      void					          SetV0sLambdaInvMassMax(Double_t mass) { fCutV0sInvMassLambdaMax = mass; }

     void      SetCascadesXiInvMassMin(Double_t mass) { fCutCascadesInvMassXiMin = mass; }
      void     SetCascadesXiInvMassMax(Double_t mass) { fCutCascadesInvMassXiMax = mass; }
   
      void      SetCascadesOmegaInvMassMin(Double_t mass) { fCutCascadesInvMassOmegaMin = mass; }
      void     SetCascadesOmegaInvMassMax(Double_t mass) { fCutCascadesInvMassOmegaMax = mass; }






      void					          SetV0sK0sCPAMin(Double_t cpa) { fCutV0sCPAK0sMin = cpa; }
      void					          SetV0sLambdaCPAMin(Double_t cpa) { fCutV0sCPALambdaMin = cpa; }
      void					          SetV0sK0sNumTauMax(Double_t nTau) { fCutV0sNumTauK0sMax = nTau; }
      void					          SetV0sLambdaNumTauMax(Double_t nTau) { fCutV0sNumTauLambdaMax = nTau; }
      void					          SetV0sK0sArmenterosAlphaMin(Double_t alpha) { fCutV0sArmenterosAlphaK0sMin = alpha; }
      void					          SetV0sLambdaArmenterosAlphaMax(Double_t alpha) { fCutV0sArmenterosAlphaLambdaMax = alpha; }
      void                    SetV0sK0sPionNumTPCSigmaMax(Float_t nSigma) { fCutV0sK0sPionNumTPCSigmaMax = nSigma; }
      void                    SetV0sLambdaPionNumTPCSigmaMax(Float_t nSigma) { fCutV0sLambdaPionNumTPCSigmaMax = nSigma; }
      void                    SetV0sLambdaProtonNumTPCSigmaMax(Float_t nSigma) { fCutV0sLambdaProtonNumTPCSigmaMax = nSigma; }
      // phi setters
      void					          SetPhiInvMassMin(Double_t mass) { fCutPhiInvMassMin = mass; }
      void					          SetPhiInvMassMax(Double_t mass) { fCutPhiInvMassMax = mass; }




      //--------------------------cascade----------------------------------------------------------------
      void  SetXiEtaMin(Double_t value){fXiPseMin = value;}
      void  SetXiEtaMax(Double_t value){fXiPseMax = value;}
      void  SetV0RadiusXiMin(Double_t value){fV0RadiusXiMin = value;}
      void  SetV0RadiusXiMax(Double_t value){fV0RadiusXiMax = value;}
      void  SetXiRadiusMin(Double_t value){fXiRadiusMin = value;}
      void  SetXiRadiusMax(Double_t value){fXiRadiusMax = value;}
      void  SetDCAXiDaughtersMax(Double_t value){fdcaXiDaughtersMax = value;}
      void  SetXiCosOfPointingAngleMin(Double_t value){fXiCosOfPointingAngleMin = value;}
      void  SetDCAV0ToPrimaryVtxXiMin(Double_t value){fdcaV0ToPrimaryVtxXiMin = value;}
      void  SetDCABachToPrimaryVtxXiMin(Double_t value){fdcaBachToPrimaryVtxXiMin = value;}
      void  SetLambdaMassWindow(Double_t value){fLambdaMassWind = value;}
      void  SetDCAV0DaughtersXi(Double_t value){fdcaV0DaughtersXi = value;}
      void  SetV0CosOfPointingAngleXiMin(Double_t value){fV0CosOfPointingAngleXiMin = value;}
      void  SetDCAPosToPrimaryVtxXiMin(Double_t value){fdcaPosToPrimaryVtxXiMin = value;}
      void  SetDCANegToPrimaryVtxXiMin(Double_t value){fdcaNegToPrimaryVtxXiMin = value;}
      void  SetCascadesRejectKinks(Bool_t reject) { fCutCascadesrejectKinks = reject; }
      void  SetXiPIDSigma(Double_t value){fXiPIDsigma = value;}
      void  SetXiMasswindow(Double_t value){fXiMasswindow = value;}




      void  SetRPTrackFromTPC(Bool_t value){fRPFromTPC = value;}
      void  SetTrackEta(Double_t value){fTrackEta = value;}
      void  SetTrackPtMin(Double_t value){fTrackPtMin = value;}
      void  SetTPCNcls(Double_t ncls = 70){fTPCNcls = ncls;}

 






      AliEventCuts            fEventCuts; //

    private:
      static const Int_t      fPIDNumSpecies = 5; // Number of considered species for PID
      static const Int_t      fFlowNumHarmonicsMax = 7; // maximum harmonics length of flow vector array
      static const Int_t      fFlowNumWeightPowersMax = 5; // maximum weight power length of flow vector array
      static const Int_t      fV0sNumBinsMass = 60; // number of InvMass bins for V0s distribution
      static const Int_t      fPhiNumBinsMass = 60; // number of InvMass bins for phi distribution
      static const Int_t      fCascadesNumBinsMass = 120;

      static const Int_t      fiNumIndexQA = 2; // QA indexes: 0: before cuts // 1: after cuts

      const char*             GetSpeciesName(PartSpecies species) const;
      const char*             GetSpeciesLabel(PartSpecies species) const;
      const char*             GetEtaGapName(Double_t dEtaGap) const { return Form("%02.2g",10.0*dEtaGap); }

      Bool_t                  InitializeTask(); // called once on beginning of task (within CreateUserObjects method)
      Bool_t                  LoadWeights(); // load weights histograms
      Bool_t                  FillFlowWeight(AliVTrack* track, PartSpecies species) const; // fill distribution for per-particle flow weight
      Double_t                GetFlowWeight(AliVTrack* track, PartSpecies species) const; // extract per-particle flow weight from input file
      Bool_t                  FillFlowWeightCascade(const AliAODcascade* xi, PartSpecies species) const; // fill distribution for per-particle flow weight
      Double_t                GetFlowWeightCascade(const AliAODcascade* xi, PartSpecies species) const; // extract per-particle flow weight from input file

      void                    ListParameters() const; // list all task parameters
      void                    ClearVectors(); // properly clear all particle vectors
      void                    DumpTObjTable(const char* note = 0x0, Option_t* opt = "") const; // add a printf statmenet given by note followed by gObjTable->Print() dump

      Bool_t                  IsEventSelected(); // event selection for Run 2 using AliEventCuts
      Bool_t                  IsEventRejectedAddPileUp() const; // additional pile-up rejection for Run2 Pb-Pb
      Int_t                   GetSamplingIndex() const; // returns sampling index based on sampling selection (number of samples)
      Int_t                   GetCentralityIndex() const; // returns centrality index based centrality estimator or number of selected tracks
      const char*             GetCentEstimatorLabel(CentEst est) const; // returns mult/cent estimator string with label or 'n/a' if not available

      void                    CalculateCorrelations(CorrTask* task, PartSpecies species, Double_t dPt = -1.0, Double_t dMass = -1.0) const; // wrapper for correlations methods
      Bool_t                  ProcessCorrTask(CorrTask* task); // procesisng of CorrTask
      Bool_t                  CalculateFlow(); // main (envelope) method for flow calculations in selected events

      void                    FilterCharged() const; // charged tracks filtering
      void                    FilterPID() const; // pi,K,p filtering
      void                    FilterV0s() const; // K0s, Lambda, ALambda filtering
      void                    FilterPhi() const; // reconstruction and filtering of Phi meson candidates

      AliAODMCParticle*       GetMCParticle(Int_t label) const; // find corresponding MC particle from fArrayMC depending of AOD track label
      Double_t                GetRapidity(Double_t mass, Double_t Pt, Double_t Eta) const; // calculate particle / track rapidity
      Bool_t                  HasMass(PartSpecies spec) const { return (spec == kK0s || spec == kLambda || spec == kPhi|| spec == kXi ||spec == kOmega); }
      Bool_t                  HasTrackPIDTPC(const AliAODTrack* track) const; // is TPC PID OK for this track ?
      Bool_t                  HasTrackPIDTOF(const AliAODTrack* track) const; // is TOF PID OK for this track ?
      Bool_t                  IsWithinRefs(const AliAODTrack* track) const; // check if track fulfill requirements for Refs (used for refs selection & autocorelations)
      Bool_t                  IsChargedSelected(const AliAODTrack* track = 0x0) const; // charged track selection
      PartSpecies             IsPIDSelected(const AliAODTrack* track) const; // PID tracks selections
      Bool_t                  IsV0Selected(const AliAODv0* v0 = 0x0) const; // general (common) V0 selection
      Bool_t                  IsV0aK0s(const AliAODv0* v0 = 0x0) const; // V0 selection: K0s specific
      Int_t                   IsV0aLambda(const AliAODv0* v0 = 0x0) const; // V0 selection: (A)Lambda specific
      AliPicoTrack*           MakeMother(const AliAODTrack* part1, const AliAODTrack* part2) const; // Combine two prongs into a mother particle stored in AliPicoTrack object
      void                    FillSparseCand(THnSparse* sparse, AliVTrack* track) const; // Fill sparse histogram for inv. mass distribution of candidates (V0s,Phi)


      void                    FillEventsQA(Int_t iQAindex) const; // filling QA plots related to event selection
      void                    FillQARefs(Int_t iQAindex, const AliAODTrack* track = 0x0) const; // filling QA plots for RFPs selection
      void                    FillQACharged(Int_t iQAindex, const AliAODTrack* track = 0x0) const; // filling QA plots for charged track selection
      void                    FillQAPID(Int_t iQAindex, const AliAODTrack* track = 0x0, PartSpecies species = kUnknown) const; // filling pi,K,p QA histograms
      void                    FillQAV0s(Int_t iQAindex, const AliAODv0* v0 = 0x0, Bool_t bIsK0s = kTRUE, Int_t bIsLambda = 2) const; // filling QA plots for V0s candidates
      void                    FillQAPhi(Int_t iQAindex, const AliPicoTrack* part = 0x0) const; // filling QA plots for V0s candidates

      // Flow related methods
      void                    FillRefsVectors(Double_t dGap); // fill flow vector Q with RFPs for reference flow
      Int_t                   FillPOIsVectors(Double_t dEtaGap, PartSpecies species, Double_t dPtLow, Double_t dPtHigh, Double_t dMassLow = 0.0, Double_t dMassHigh = 0.0); // fill flow vectors p,q and s with POIs (for given species) for differential flow calculations
      void                    ResetFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]); // set values to TComplex(0,0,0) for given array
      void                    ListFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]) const; // printf all values of given Flow vector array


//-----------------------------------------------------------------------------------------
   
      void                    FilterCascades() const;

      Bool_t                  IsCascadeSelected(const AliAODcascade *xi) const;

      Bool_t                  IsSelected(const AliAODTrack *t =0x0) const;


      void                    Propagate( Double_t vv[3],Double_t x[3],Double_t p[3],Double_t bz,Double_t sign) const;
      Double_t                PIDCorrectionHF(const AliAODTrack *track, const Int_t ispecies) const;

       Bool_t  Is2018Data;//
       Bool_t  IsPIDorrection;//
       Bool_t  IsAdditional2018DataEventCut;//
       Double_t fXiPseMin;//
       Double_t fXiPseMax;//
       Double_t fV0RadiusXiMin;//
       Double_t fV0RadiusXiMax;//
       Double_t fXiRadiusMin;//
       Double_t fXiRadiusMax;//
       Double_t fdcaXiDaughtersMax;//
       Double_t fXiCosOfPointingAngleMin;//
       Double_t fdcaV0ToPrimaryVtxXiMin;//
       Double_t fdcaBachToPrimaryVtxXiMin;//
       Double_t fLambdaMassWind;//
       Double_t fdcaV0DaughtersXi;//
       Double_t fV0CosOfPointingAngleXiMin;//
       Double_t fdcaPosToPrimaryVtxXiMin;//
       Double_t fdcaNegToPrimaryVtxXiMin;//
       Bool_t   fCutCascadesrejectKinks; // Reject Kink cascade daughter tracks ?
       Double_t fXiPIDsigma;//
       Double_t fXiMasswindow;//
       //----------------track-------------------------------------------------  
       Double_t fTPCNcls; // number of TPC clusters   
       Double_t fTrackEta;//
       Double_t fTrackPtMin;//
       Bool_t   fRPFromTPC;//

//----------------------------------------------------------------------------------------
      TComplex                Q(Int_t n, Int_t p) const;//
      TComplex                QGapPos(Int_t n, Int_t p) const;//
      TComplex                QGapNeg(Int_t n, Int_t p) const;//
      TComplex                QGapMid(Int_t n, Int_t p) const;//
      TComplex                P(Int_t n, Int_t p) const;//
      TComplex                PGapPos(Int_t n, Int_t p) const;//
      TComplex                PGapNeg(Int_t n, Int_t p) const;//
      TComplex                S(Int_t n, Int_t p) const;//
      TComplex                SGapPos(Int_t n, Int_t p) const;//
      TComplex                SGapNeg(Int_t n, Int_t p) const;//

      TComplex                Two(Int_t n1, Int_t n2) const; // Two particle reference correlation calculations (no eta gap)
      TComplex                TwoGap(Int_t n1, Int_t n2) const; // Two particle reference correlation calculations (with eta gap)
      TComplex                Three(Int_t n1, Int_t n2, Int_t n3) const; // Three particle reference correlation calculations (no eta gap)
      TComplex                Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle reference correlation calculations (no eta gap)
      TComplex                FourGap(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle reference correlation calculations (no eta gap)
      TComplex                Four3sub(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle reference correlation calculations (with 3 sub-events)

      TComplex                TwoDiff(Int_t n1, Int_t n2) const; // Two particle diff. correlation calculations (no eta gap)
      TComplex                TwoDiffGapPos(Int_t n1, Int_t n2) const; // Two particle diff. correlation calculations (with eta gap)
      TComplex                TwoDiffGapNeg(Int_t n1, Int_t n2) const; // Two particle diff. correlation calculations (with eta gap)
      TComplex                ThreeDiff(Int_t n1, Int_t n2, Int_t n3) const; // Three particle diff. correlation calculations (no eta gap)
      TComplex                ThreeDiffGapPos(Int_t n1, Int_t n2, Int_t n3) const; // Three particle diff. correlation calculations (with eta gap)
      TComplex                ThreeDiffGapNeg(Int_t n1, Int_t n2, Int_t n3) const; // Three particle diff. correlation calculations (with eta gap)
      TComplex                FourDiff(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle reference correlation calculations (no eta gap)
      TComplex                FourDiffGapPos(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle reference correlation calculations (with eta gap)
      TComplex                FourDiffGapNeg(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle reference correlation calculations (with eta gap)

      // array lenghts & constants
      const Double_t          fPDGMassPion; // [DPGMass] DPG mass of charged pion
      const Double_t          fPDGMassKaon; // [DPGMass] DPG mass of charged kaon
      const Double_t          fPDGMassProton; // [DPGMass] DPG mass of proton
      const Double_t          fPDGMassPhi; // [DPGMass] DPG mass of phi (333) meson
      const Double_t          fPDGMassK0s; // [DPGMass] DPG mass of K0s
      const Double_t          fPDGMassLambda; // [DPGMass] DPG mass of (Anti)Lambda

      AliAODEvent*            fEventAOD; //! AOD event countainer
      Double_t                fPVz; // PV z-coordinate used for weights
      AliPIDResponse*         fPIDResponse; //! AliPIDResponse container
      AliPIDCombined*         fPIDCombined; //! AliPIDCombined container
      TFile*                  fFlowWeightsFile; //! source file containing weights
      TFile*                  fPIDCorrectionFile; //! source file containing weights
      TClonesArray*           fArrayMC; //! input list of MC particles
      Bool_t                  fMC; // is running on mc?
      Bool_t                  fInit; // initialization check
      Int_t                   fIndexSampling; // sampling index (randomly generated)
      Int_t                   fIndexCentrality; // centrality bin index (based on centrality est. or number of selected tracks)
      Int_t                   fEventCounter; // event counter (used for local test runmode purpose)
      Int_t                   fNumEventsAnalyse; // [50] number of events to be analysed / after passing selection (only in test mode)
      Int_t                   fRunNumber; // [-1] run number of previous event (not the current one)

      TComplex                fFlowVecQpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecQneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecQmid[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecPpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecPneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecSpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecSneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation

      std::vector<CorrTask*>  fVecCorrTask; //
      std::vector<AliVTrack*>* fVector[kUnknown]; //! container for selected Refs charged particles

      //cuts & selection: analysis
      RunMode                 fRunMode; // running mode (not grid related)
      AnalType                fAnalType; // analysis type: AOD / ESD
      Bool_t                  fDumpTObjectTable; // [kFALSE] flag for dumping TObjectTable to the output stream
      Bool_t                  fSampling; // [kFALSE] Do random sampling ? (estimation of vn stat. uncertanity)
      Bool_t                  fFillQA; //[kTRUE] flag for filling the QA plots
      Bool_t                  fProcessSpec[kUnknown];  // [false] flag for processing species
      // cuts & selection: flow related
      Double_t                fFlowRFPsPtMin; // [0.2] (GeV/c) min pT treshold for RFPs particle for reference flow
      Double_t                fFlowRFPsPtMax; // [5.0] (GeV/c) max pT treshold for RFPs particle for reference flow
      Double_t                fFlowPOIsPtMin; // [0] (GeV/c) min pT treshold for POIs for differential flow
      Double_t                fFlowPOIsPtMax; // [10] (GeV/c) max pT treshold for POIs for differential flow
      Int_t                   fFlowPOIsPtBinNum; // [0] number of pt bins
      Double_t                fFlowEtaMax; // [0.8] max eta acceptance for flow particles (RFPs & POIs)
      Int_t                   fFlowEtaBinNum; // [0] number of eta bins
      Int_t                   fFlowPhiBinNum; // [100] number of phi bins
      Int_t                   fNumSamples; // [1] overall number of samples (from random sampling) used
      Bool_t                  fFlowFillWeights; //[kFALSE] flag for filling weights
      Bool_t                  fFlowUseWeights; //[kFALSE] flag for using the previously filled weights (NOTE: this is turned on only when path to file is applied via fFlowWeightsPath)
      Bool_t                  fFlowUse3Dweights; // [kFALSE] flag for using 3D GF weights, if kFALSE, 2D weights are expected
      Bool_t                  fFlowRunByRunWeights; // [kTRUE] flag for using rub-by-run weigths from weigths file; if false, only one set of histrograms is provided
      TString                 fFlowWeightsPath; //[] path to source root file with weigthts (if empty unit weights are applied) e.g. "alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/PhiWeight_LHC16kl.root"
      TString                 fPIDCorrectionPath;//
      //cuts & selection: events
      ColSystem               fColSystem; // collisional system
      AliVEvent::EOfflineTriggerTypes    fTrigger; // physics selection trigger
      CentEst                 fCentEstimator; // [kV0A] multiplicity/centrality estimator as in AliMultSelection
      Int_t                   fCentMin; // [0] min range for centrality/multiplicity histos
      Int_t                   fCentMax; // [0] max range for centrality/multiplicity histos
      Int_t                   fCentBinNum; // [0] number of centrality bins
      Double_t                fPVtxCutZ; // (cm) PV z cut
      Bool_t                  fEventRejectAddPileUp; // additional pile-up rejection for Pb-Pb collisions in Run2 (17n, 15o)
      //cuts & selection: tracks
      UInt_t                  fCutChargedTrackFilterBit; // (-) tracks filter bit
      UShort_t                fCutChargedNumTPCclsMin;  // (-) Minimal number of TPC clusters used for track reconstruction
      Double_t                fCutChargedDCAzMax; // (cm) Maximal DCA-z cuts for tracks (pile-up rejection suggested for LHC16)
      Double_t                fCutChargedDCAxyMax; // (cm) Maximal DCA-xy cuts for tracks (pile-up rejection suggested for LHC16)
      // cuts & selection: PID selection
      Bool_t                  fCutPIDUseAntiProtonOnly; // [kFALSE] check proton PID charge to select AntiProtons only
      Bool_t                  fCutPIDnSigmaCombinedTOFrejection; // [kTRUE] flag for rejection candidates in TPC+TOF pt region if TOF is not available (if true and no TOF track is skipped, otherwise only TPC is used)
      Bool_t                  fCutUseBayesPID; // [kFALSE] flag for using Bayes PID for pi,K,p instead nsigma cut
      Float_t                 fCutPIDnSigmaTPCRejectElectron; // [0] number of TPC nSigma for electron rejection
      Float_t                 fCutPIDnSigmaMax[fPIDNumSpecies]; // [0] maximum of nSigmas (TPC or TPC & TOF combined)
      Double_t                fCutPIDBayesMin[fPIDNumSpecies]; // [0.0] minimal value of Bayes PID probability for pion
      //cuts & selection: V0 reconstruction
      Bool_t                  fCutV0sOnFly;		// V0 reconstruction method: is On-the-fly? (or offline)
      Bool_t                  fCutV0srefitTPC; // Check TPC refit of V0 daughters ?
      Bool_t                  fCutV0srejectKinks; // Reject Kink V0 daughter tracks ?
      UShort_t                fCutV0sDaughterNumTPCClsMin; // min number of TPC clusters
      UShort_t                fCutV0sDaughterNumTPCCrossMin; // min number of crossed TPC rows
      UShort_t                fCutV0sDaughterNumTPCFindMin; // min number of findable TPC clusters
      UShort_t                fCutV0sDaughterNumTPCClsPIDMin; // min number of TPC clusters used for PID
      Double_t                fCutV0sDaughterRatioCrossFindMin; // min ratio of crossed / findable TPC clusters
      Bool_t				          fCutV0sCrossMassRejection; // competing V0 rejection based on InvMass
      Double_t                fCutV0sCrossMassCutK0s; // [0.005] (GeV/c2) restricted vicinity of Lambda/ALambda inv. mass peak for K0s candidates
      Double_t                fCutV0sCrossMassCutLambda; // [0.020] (GeV/c2) restricted vicinity of K0s inv. mass peak for Lambda/ALambda candidates
      Double_t                fCutV0sDCAtoPVMin;   // (cm) min DCA of V0 daughter to PV
      Double_t				        fCutV0sDCAtoPVMax;	// (cm) max DCA of V0 daughter to PV
      Double_t                fCutV0sDCAtoPVzMax; // (cm) max DCA-z coordinate of V0 daughters to PV
      Double_t				        fCutV0sDCADaughtersMin;	// (cm) min DCA of V0 daughters among themselves
      Double_t				        fCutV0sDCADaughtersMax;	// (cm) max DCA of V0 daughters among themselves
      Double_t                fCutV0sDecayRadiusMin; // (cm) min distance of secondary vertex from z-axis in transverse plane
      Double_t				        fCutV0sDecayRadiusMax; // (cm) max distance of secondary vertex from z-axis in transverse plane
      UInt_t                  fCutV0sDaughterFilterBit; // (-) V0 daughters filter bit
      Double_t                fCutV0sDaughterPtMin; // (GeV/c) min pT of V0 daughters
      Double_t                fCutV0sDaughterPtMax; // (GeV/c) max pT of V0 daughters
      Double_t                fCutV0sDaughterEtaMax; // (-) max value of Eta of V0 daughters
      Double_t                fCutV0sMotherRapMax; // (-) max rapidity value of V0 mother
      Double_t                fCutV0sCPAK0sMin;    // (-) min cosine of pointing angle of K0s candidate to PV
      Double_t                fCutV0sCPALambdaMin; // (-) min cosine of pointing angle of K0s candidate to PV
      Double_t                fCutV0sNumTauK0sMax; // (c*tau) max number of c*tau (K0s)
      Double_t                fCutV0sNumTauLambdaMax; // (c*tau) max number of c*tau ((A)Lambda)
      Double_t                fCutV0sInvMassK0sMin; // [0.4] (GeV/c2) min inv. mass window for selected K0s candidates
      Double_t                fCutV0sInvMassK0sMax; // [0.6] (GeV/c2) max inv. mass window for selected K0s candidates
      Double_t                fCutV0sInvMassLambdaMin; // [1.08] (GeV/c2) min inv. mass window for selected (Anti)Lambda candidates
      Double_t                fCutV0sInvMassLambdaMax; // [1.16] (GeV/c2) max inv. mass window for selected (Anti)Lambda candidates

      Double_t                fCutCascadesInvMassXiMin;//
      Double_t                fCutCascadesInvMassOmegaMin;//


      Double_t                fCutCascadesInvMassXiMax;//
      Double_t                fCutCascadesInvMassOmegaMax;//


      Double_t		      fCutV0sArmenterosAlphaK0sMin; // (alpha) min Armenteros alpha for K0s
      Double_t                fCutV0sArmenterosAlphaLambdaMax; // (alpha) max Armenteros alpha for (Anti)Lambda
      Float_t                 fCutV0sK0sPionNumTPCSigmaMax; // (sigmaTPC) max number of TPC sigmas for kaon PID (K0s candidates)
      Float_t                 fCutV0sLambdaPionNumTPCSigmaMax;    // (sigmaTPC) max number of TPC sigma for pion PID (Lambda candidates)
      Float_t                 fCutV0sLambdaProtonNumTPCSigmaMax;    // (sigmaTPC) max number of TPC sigma for proton PID (Lambda candidates)
      // cuts & selection: phi
      Double_t                fCutPhiInvMassMin; // [0.99] (GeV/c2) min inv. mass window for selected phi candidates
      Double_t                fCutPhiInvMassMax; // [1.07] (GeV/c2) min inv. mass window for selected phi candidates

      // output lists
      AliDirList*                  fQAEvents; //! events list
      TList*                       fQAEventCut; //! events list
      AliDirList*                  fQACharged; //! charged tracks list
      AliDirList*                  fQAPID; //! pi,K,p list
      AliDirList*                  fQAV0s; //! V0s candidates list
      AliDirList*                  fQAPhi; //! Phi candidates list
      AliDirList*                  fFlowWeights; //! list for flow weights
      AliDirList*                  fListFlow[kUnknown]; //! flow lists
      // histograms & profiles

      // Flow
      THnSparseD*             fhsCandK0s; //! distribution of K0s candidates
      THnSparseD*             fhsCandLambda; //!  distribution of Lambda candidates
      THnSparseD*             fhsCandXi; //! distribution of xi candidates
      THnSparseD*             fhsCandOmega; //!  distribution of omega candidates


      THnSparseD*             fhsCandPhi; //!  distribution of Phi candidates
      THnSparseD*             fhsCandPhiBg; //!  distribution of Phi background

      TH2D*                   fh2Weights[kUnknown]; //! container for GF weights (phi,eta,pt) (2D)
      TH3D*                   fh3Weights[kUnknown]; //! container for GF weights (phi,eta,pt)
      TH2D*                   fh2AfterWeights[kUnknown]; //! distribution after applying GF weights - lightweight QA (phi)
      TH3D*                   fh3AfterWeights[kUnknown]; //! distribution after applying GF weights - full QA (phi,eta,pt)

      // Events
      TH2D*                   fhEventSampling; //! distribution of sampled events (based on randomly generated numbers)
      TH1D*                   fhEventCentrality; //! distribution of event centrality
      TH2D*                   fh2EventCentralityNumRefs; //! distribution of event centrality vs number of selected charged tracks
      TH1D*                   fhEventCounter; //! counter following event selection
      // Charged
      TH1D*                   fhRefsMult; //!multiplicity distribution of selected RFPs
      TH1D*                   fhRefsPt; //! pt distribution of selected RFPs
      TH1D*                   fhRefsEta; //! pt distribution of selected RFPs
      TH1D*                   fhRefsPhi; //! pt distribution of selected RFPs
      TProfile*               fpRefsMult; //! <multiplicity>
      TH1D*                   fhChargedCounter; //! counter following charged track selection
      // PID
      TH1D*                   fhPIDCounter; //! counter for PID
      TH1D*                   fhPIDMult[3]; //! multiplicity distribution of selected pions
      TH1D*                   fhPIDPt[3]; //! pt distribution of selected pions
      TH1D*                   fhPIDPhi[3]; //! phi distribution of selected pions
      TH1D*                   fhPIDEta[3]; //! eta distribution of selected pions
      TH1D*                   fhPIDCharge[3]; //! charge distribution of selected pions
      TH2D*                   fh2PIDTPCdEdx[3]; //! TPC dEdx response of selected pions
      TH2D*                   fh2PIDTPCdEdxDelta[3]; //! TPC delta dEdx (measured - expected) response of selected pions
      TH2D*                   fh2PIDTOFbeta[3]; //! TOF beta of selected pions
      TH2D*                   fh2PIDTOFbetaDelta[3]; //! TOF delta beta (measured - expected) of selected pions
      TH2D*                   fh2PIDBayesElectron[3]; //! Bayesian PID probability vs pT for selected pions (pion hypothesis)
      TH2D*                   fh2PIDBayesMuon[3]; //! Bayesian PID probability vs pT for selected pions (pion hypothesis)
      TH2D*                   fh2PIDBayesPion[3]; //! Bayesian PID probability vs pT for selected pions (pion hypothesis)
      TH2D*                   fh2PIDBayesKaon[3]; //! Bayesian PID probability vs pT for selected pions (kaon hypothesis)
      TH2D*                   fh2PIDBayesProton[3]; //! Bayesian PID probability vs pT for selected pions (proton hypothesis)
      TH2D*                   fh2PIDTPCnSigmaElectron[3]; //! TPC nSigma vs pT for selected pions (pion hypothesis)
      TH2D*                   fh2PIDTOFnSigmaElectron[3]; //! TOF nSigma vs pT for selected pions (pion hypothesis)
      TH2D*                   fh2PIDTPCnSigmaMuon[3]; //! TPC nSigma vs pT for selected pions (pion hypothesis)
      TH2D*                   fh2PIDTOFnSigmaMuon[3]; //! TOF nSigma vs pT for selected pions (pion hypothesis)
      TH2D*                   fh2PIDTPCnSigmaPion[3]; //! TPC nSigma vs pT for selected pions (pion hypothesis)
      TH2D*                   fh2PIDTOFnSigmaPion[3]; //! TOF nSigma vs pT for selected pions (pion hypothesis)
      TH2D*                   fh2PIDTPCnSigmaKaon[3]; //! TPC nSigma vs pT for selected pions (kaon hypothesis)
      TH2D*                   fh2PIDTOFnSigmaKaon[3]; //! TOF nSigma vs pT for selected pions (kaon hypothesis)
      TH2D*                   fh2PIDTPCnSigmaProton[3]; //! TPC nSigma vs pT for selected pions (proton hypothesis)
      TH2D*                   fh2PIDTOFnSigmaProton[3]; //! TOF nSigma vs pT for selected pions (proton hypothesis)

      TH1D*                   fhMCRecoSelectedPionPt; //! pt dist of selected (MC reco) pions
      TH1D*                   fhMCRecoSelectedTruePionPt; //! pt dist of selected (MC reco) true (tagged in MC gen) pions
      TH1D*                   fhMCRecoAllPionPt; //! pt dist of all (MC reco) pions (i.e. selected charged tracks that are tagged in MC)
      TH1D*                   fhMCGenAllPionPt; //! pt dist of all (MC) generated pions
      TH1D*                   fhMCRecoSelectedKaonPt; //! pt dist of selected (MC reco) Kaons
      TH1D*                   fhMCRecoSelectedTrueKaonPt; //! pt dist of selected (MC reco) true (tagged in MC gen) Kaons
      TH1D*                   fhMCRecoAllKaonPt; //! pt dist of all (MC reco) Kaons (i.e. selected charged tracks that are tagged in MC)
      TH1D*                   fhMCGenAllKaonPt; //! pt dist of all (MC) generated Kaons
      TH1D*                   fhMCRecoSelectedProtonPt; //! pt dist of selected (MC reco) Protons
      TH1D*                   fhMCRecoSelectedTrueProtonPt; //! pt dist of selected (MC reco) true (tagged in MC gen) Protons
      TH1D*                   fhMCRecoAllProtonPt; //! pt dist of all (MC reco) Protons (i.e. selected charged tracks that are tagged in MC)
      TH1D*                   fhMCGenAllProtonPt; //! pt dist of all (MC) generated Protons
      // Phi
      TH1D*                   fhPhiCounter; //! counter following phi candidate selection
      TH1D*                   fhPhiMult; //! multiplicity distribution of selected phi candidates
      TH1D*                   fhPhiBGMult; //! multiplicity distribution of BG candidates
      TH1D*                   fhPhiInvMass; //! invariant mass distribution of phi candidates
      TH1D*                   fhPhiBGInvMass; //! invariant mass distribution of phi background candidates
      TH1D*                   fhPhiCharge; //! charge distribution of selected phi candidates
      TH1D*                   fhPhiBGCharge; //! charge distribution of phi BG candidates
      TH1D*                   fhPhiPt; //! pt distribution of selected phi candidates
      TH1D*                   fhPhiEta; //! eta distribution of selected phi candidates
      TH1D*                   fhPhiPhi; //! phi distribution of selected phi candidates
      // V0s
      TH1D*                   fhV0sCounter; //! counter following V0s selection
      TH1D*                   fhV0sCounterK0s; //! counter following K0s selection
      TH1D*                   fhV0sCounterLambda; //! counter following (Anti-)Lambda selection
      TH2D*                   fhV0sInvMassK0s; //! 2D inv. mass distiburion (K0s mass vs. Lambda/AntiLambda mass)
      TH2D*                   fhV0sInvMassLambda; //! 2D inv. mass distiburion (K0s mass vs. Lambda/AntiLambda mass)
      TH2D*                   fhV0sCompetingInvMassK0s; //! dist of InvMass of rejected K0s candidates in (Anti-)Lambda peak
      TH2D*                   fhV0sCompetingInvMassLambda; //! dist of InvMass of rejected (Anti-)Lambda candidates in K0s peak
      // PID correction
      TH3D*                   histMeanPionq;//!
      TH3D*                   histMeanKaonq;//!
      TH3D*                   histMeanProtonq;//!
      TH3D*                   histSigmaPionq;//!
      TH3D*                   histSigmaKaonq;//!
      TH3D*                   histSigmaProtonq;//!
      TH3D*                   histMeanPionr;//!
      TH3D*                   histMeanKaonr;//!
      TH3D*                   histMeanProtonr;//!
      TH3D*                   histSigmaPionr;//!
      TH3D*                   histSigmaKaonr;//!
      TH3D*                   histSigmaProtonr;//!

      // QA: events
      TH1D*                   fhQAEventsPVz[fiNumIndexQA]; //!
      TH1D*                   fhQAEventsNumContrPV[fiNumIndexQA]; //!
      TH1D*                   fhQAEventsNumSPDContrPV[fiNumIndexQA]; //!
      TH1D*                   fhQAEventsDistPVSPD[fiNumIndexQA]; //!
      TH1D*                   fhQAEventsSPDresol[fiNumIndexQA]; //!
      TH2D*                   fhQAEventsfMult32vsCentr;//!
      TH2D*                   fhQAEventsMult128vsCentr;//!
      TH2D*                   fhQAEventsfMultTPCvsTOF;//!
      TH2D*                   fhQAEventsfMultTPCvsESD;//!
      // QA: charged tracks
      TH1D*                   fhQAChargedMult[fiNumIndexQA];       //! number of AOD charged tracks distribution
      TH1D*                   fhQAChargedPt[fiNumIndexQA];         //! pT dist of charged tracks
      TH1D*                   fhQAChargedEta[fiNumIndexQA];        //! eta dist of charged tracks
      TH1D*                   fhQAChargedPhi[fiNumIndexQA];        //! phi dist of charged tracks
      TH1D*                   fhQAChargedCharge[fiNumIndexQA];     //! charge dist of charged tracks
      TH1D*                   fhQAChargedNumTPCcls[fiNumIndexQA];  //! dist of track number of TPC clusters
      TH1D*                   fhQAChargedDCAxy[fiNumIndexQA];      //! dist of Charged DCA in transverse plane
      TH1D*                   fhQAChargedDCAz[fiNumIndexQA];       //! dist of charged DCA in z coordinate
      // QA: PID tracks
      TH1D*                   fhQAPIDTPCstatus[fiNumIndexQA];  //! based on AliPIDResponse::CheckPIDStatus();
      TH1D*                   fhQAPIDTOFstatus[fiNumIndexQA];  //! based on AliPIDResponse::CheckPIDStatus();
      TH2D*                   fhQAPIDTPCdEdx[fiNumIndexQA];    //! TPC PID information
      TH2D*                   fhQAPIDTOFbeta[fiNumIndexQA];    //! TOF PID information
      TH3D*                   fh3QAPIDnSigmaTPCTOFPtPion[fiNumIndexQA]; //! nSigma TPC vs nSigma TOF vs pt
      TH3D*                   fh3QAPIDnSigmaTPCTOFPtKaon[fiNumIndexQA]; //! nSigma TPC vs nSigma TOF vs pt
      TH3D*                   fh3QAPIDnSigmaTPCTOFPtProton[fiNumIndexQA]; //! nSigma TPC vs nSigma TOF vs pt
      // QA: V0s candidates
      TH1D*			  		        fhQAV0sMultK0s[fiNumIndexQA];	//! number of K0s candidates
      TH1D*			  		        fhQAV0sMultLambda[fiNumIndexQA];	//! number of Lambda candidates
      TH1D*			  		        fhQAV0sMultALambda[fiNumIndexQA];	//! number of Anti-Lambda candidates
      TH1D*			  		        fhQAV0sRecoMethod[fiNumIndexQA];	//! offline/online V0 reconstruction method
      TH1D*			  		        fhQAV0sDaughterTPCRefit[fiNumIndexQA];	//! Daughters TPC refit true/false
      TH1D*			  		        fhQAV0sDaughterKinks[fiNumIndexQA];	//! Daughters kinks true/false
      TH1D*                   fhQAV0sDaughterNumTPCCls[fiNumIndexQA]; //! Daughter # of TPC findable clusters
      TH1D*                   fhQAV0sDaughterNumTPCFind[fiNumIndexQA]; //! Daughter # of TPC clusters
      TH1D*                   fhQAV0sDaughterNumTPCCrossRows[fiNumIndexQA]; //! Daughter # of TPC crossed rows
      TH1D*                   fhQAV0sDaughterTPCCrossFindRatio[fiNumIndexQA]; //! Daughter # of TPC cross / # of TPC findable cls ratio
      TH1D*                   fhQAV0sDaughterNumTPCClsPID[fiNumIndexQA]; //! Daughter # of TPC findable clusters used for PID
      TH1D*			  		        fhQAV0sDCAtoPV[fiNumIndexQA];	//! V0 DCA to PV
      TH1D*			  		        fhQAV0sDCADaughters[fiNumIndexQA];	//! DCA between V0 daughters
      TH1D*			  		        fhQAV0sDecayRadius[fiNumIndexQA];	//! Distance between PV and Secondary vertex in transverse plane
      TH1D*                   fhQAV0sInvMassK0s[fiNumIndexQA];    //! inv. mass dist of V0s (K0s mass hypothesis)
      TH1D*					          fhQAV0sInvMassLambda[fiNumIndexQA];	//! inv. mass dist of V0s ((A)Lambda mass hypothesis)
      TH1D*                   fhQAV0sMotherPt[fiNumIndexQA];  //! pT dist of V0s
      TH1D*					          fhQAV0sMotherPhi[fiNumIndexQA];	//! azimuthal dist of V0s
      TH1D*                   fhQAV0sMotherEta[fiNumIndexQA]; //! pseudorapidity dist of V0s
      TH1D*                   fhQAV0sMotherCharge[fiNumIndexQA]; //! charge distribution of mothers
      TH1D*                   fhQAV0sMotherRapK0s[fiNumIndexQA];  //! rapidity dist of V0s (K0s mass hypothesis)
      TH1D*                   fhQAV0sMotherRapLambda[fiNumIndexQA]; //! rapidity dist of V0s (Lambda mass hypothesis)
      TH1D*                   fhQAV0sDaughterPt[fiNumIndexQA];    //! pT dist of V0 daughters
      TH1D*					          fhQAV0sDaughterPhi[fiNumIndexQA];	//! pT dist of V0 daughters
      TH1D*                   fhQAV0sDaughterEta[fiNumIndexQA];   //! pseudorapidity dist of V0 daughters
      TH1D*                   fhQAV0sDaughterCharge[fiNumIndexQA]; //! charge distribution of daughters
      TH1D*					          fhQAV0sDaughterTPCstatus[fiNumIndexQA];	//! TPC dEdx vs p of K0s daughters
      TH1D*					          fhQAV0sDaughterTOFstatus[fiNumIndexQA];	//! TPC dEdx vs p of K0s daughters
      TH2D*					          fhQAV0sDaughterTPCdEdxK0s[fiNumIndexQA];	//! TPC dEdx vs p of K0s daughters
      TH2D*					          fhQAV0sDaughterNumSigmaPionK0s[fiNumIndexQA];	//! Number of TPC sigmas (pion) vs mother pT of K0s daughters
      TH2D*					          fhQAV0sDaughterTPCdEdxLambda[fiNumIndexQA];	//! TPC dEdx vs p of Lambda daughters
      TH2D*                   fhQAV0sDaughterNumSigmaPionLambda[fiNumIndexQA];  //! number of TPC sigmas vs mother pT of pion (Lambda candidates)
      TH2D*                   fhQAV0sDaughterNumSigmaProtonLambda[fiNumIndexQA];  //! number of TPC sigmas vs mother pT of proton (Lambda candidates)
      TH2D*                   fhQAV0sDaughterNumSigmaPionALambda[fiNumIndexQA];  //! number of TPC sigmas vs mother pT of pion (Anti-Lambda candidates)
      TH2D*                   fhQAV0sDaughterNumSigmaProtonALambda[fiNumIndexQA];  //! number of TPC sigmas vs mother pT of proton (Anti-Lambda candidates)
      TH1D*					          fhQAV0sCPAK0s[fiNumIndexQA];	//! cosine of pointing angle of K0s candidates
      TH1D*					          fhQAV0sCPALambda[fiNumIndexQA];	//! cosine of pointing angle of Lambda candidates
      TH1D*					          fhQAV0sNumTauK0s[fiNumIndexQA];	//! number of c*tau of K0s candidates
      TH1D*					          fhQAV0sNumTauLambda[fiNumIndexQA];	//! number of c*tau of Lambda candidates
      TH2D*				   	        fhQAV0sArmenterosK0s[fiNumIndexQA];	//! Armenteros-Podolanski plot for K0s candidates
      TH2D*			  		        fhQAV0sArmenterosLambda[fiNumIndexQA];	//! Armenteros-Podolanski plot for Lambda candidates
      TH2D*			  		        fhQAV0sArmenterosALambda[fiNumIndexQA];	//! Armenteros-Podolanski plot for ALambda candidates

      ClassDef(AliAnalysisTaskUniFlowMultiStrange, 13);
};

#endif
