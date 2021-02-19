/* Copyright(c) 2016, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKUNIFLOW_H
#define ALIANALYSISTASKUNIFLOW_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "THnSparse.h"

class TString;
class TComplex;
class TFile;
class TList;
class TProfile;
class TH1D;
class TH2D;
class TH3D;
class TObjArray;

class AliPIDResponse;
class AliPIDCombined;
class AliVEvent;
class AliAODEvent;
class AliMCEvent;
class AliVParticle;
class AliVTrack;
class AliAODTrack;
class AliPicoTrack;
class AliAODv0;
class AliAODMCParticle;
class AliEventPoolManager;
class AliEventPool;

class AliUniFlowCorrTask;

//_____________________________________________________________________________

class AliAnalysisTaskUniFlow : public AliAnalysisTaskSE
{
    public:
      enum    RunMode {kFull = 0, kTest, kSkipFlow}; // task running mode (NOT GRID MODE)
      enum    ColSystem {kPP = 0, kPPb, kPbPb}; // tag for collisional system
      enum    AnalType {kAOD = 0, kESD, kMC}; // tag for analysis type
      enum    CentEst {kRFP = 0, kV0A, kV0C, kV0M, kCL0, kCL1, kZNA, kZNC}; // multiplicity/centrality estimator as AliMultSelection
      enum    PartSpecies {kRefs = 0, kCharged, kPion, kKaon, kProton, kCharUnidentified, kK0s, kLambda, kPhi, kUnknown}; // list of all particle species of interest; NB: kUknown last as counter
      enum    SparseCand {kInvMass = 0, kCent, kPt, kEta, kSample, kDim}; // reconstructed candidates dist. dimensions
      enum    SparseWeights {wPhi = 0, wCent, wPt, wEta, wVz, wSpec, wDim}; // multidimensional weights sparse.. w as weights (to avoid redefinition from the previous one)
      enum    QAindex { kBefore = 0, kAfter, kNumQA}; // index for filling QA status

                              AliAnalysisTaskUniFlow(); // constructor
                              AliAnalysisTaskUniFlow(const char *name, ColSystem colSys, Bool_t bUseWeights = kFALSE, Bool_t bIsMC = kFALSE); // named (primary) constructor
                              AliAnalysisTaskUniFlow(const AliAnalysisTaskUniFlow&); // not implemented
                              AliAnalysisTaskUniFlow& operator=(const AliAnalysisTaskUniFlow&); // not implemented
      virtual                 ~AliAnalysisTaskUniFlow(); // destructor


      virtual void            UserCreateOutputObjects(); //
      virtual void            UserExec(Option_t* option); // main methond - called for each event
      virtual void            Terminate(Option_t* option); // called after all events are processed
      // analysis setters
      void                    SetRunMode(RunMode mode = kFull) { fRunMode = mode; }
      void                    SetNumEventsAnalyse(Int_t num) { fNumEventsAnalyse = num; }
      void                    SetDumpTObjectTable(Bool_t dump = kTRUE) { fDumpTObjectTable = dump; }
      void					          SetAnalysisType(AnalType type = kAOD) { fAnalType = type; }
      void					          SetNeedPIDCorrection(Bool_t pidCorr) { fNeedPIDCorrection = pidCorr; }
      void					          SetPIDCorrectionPhi(Bool_t pidCorr) { fPIDCorrectionPhi = pidCorr; }
      void					          SetIs2018data(Bool_t is2018data) { fIs2018data = is2018data; }
      void					          SetIsHMpp(Bool_t isHMpp = kTRUE) { fIsHMpp = isHMpp; }
      void                    SetSampling(Bool_t sample = kTRUE, Int_t iNum = 10) { fSampling = sample; fNumSamples = iNum; }
      void                    SetEtaCheckRFP(Bool_t check = kFALSE) { fEtaCheckRFP = check; }
      void                    SetFillQAhistos(Bool_t fill = kTRUE) { fFillQA = fill; }
      void                    SetFillMultiDimensionalWeights(Bool_t fill = kTRUE) { fFlowFillWeightsMultiD = fill; }
      void                    SetProcessPID(Bool_t use = kTRUE) { fProcessSpec[kPion] = use; fProcessSpec[kKaon] = use; fProcessSpec[kProton] = use; }
      void                    SetProcessV0s(Bool_t use = kTRUE) { fProcessSpec[kK0s] = use; fProcessSpec[kLambda] = use; }
      void                    SetProcessK0s(Bool_t use = kTRUE) { fProcessSpec[kK0s] = use; }
      void                    SetProcessLambda(Bool_t use = kTRUE) { fProcessSpec[kLambda] = use; }
      void                    SetProcessPhi(Bool_t use = kTRUE) { fProcessSpec[kPhi] = use; }
      void                    SetDoCorrelationsUsingGF(Bool_t use = kTRUE) { fCorrUsingGF = use;}
      void                    SetUseGeneralFormula(Bool_t use = kTRUE) { fUseGeneralFormula = use;}
      void                    SetUsePIDweights(Bool_t use = kTRUE) { fFlowUsePIDWeights = use;}
      // flow related setters
      void                    AddCorr(std::vector<Int_t> harms, std::vector<Double_t> gaps = std::vector<Double_t>(), Bool_t doRFPs = kTRUE, Bool_t doPOIs = kTRUE, std::vector<Int_t> maxPowVec = {});
      // void                    AddCorr(std::vector<Int_t> harms, std::vector<Double_t> gaps = std::vector<Double_t>(), Bool_t doRFPs = kTRUE, Bool_t doPOIs = kTRUE) { fVecCorrTask.push_back(new AliUniFlowCorrTask(doRFPs, doPOIs, harms, gaps)); }
      void                    AddTwo(Int_t n1, Int_t n2, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { AddCorr({n1,n2},{},refs,pois); }
      void                    AddTwoGap(Int_t n1, Int_t n2, Double_t gap, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { AddCorr({n1,n2},{gap},refs,pois); }
      void                    AddThree(Int_t n1, Int_t n2, Int_t n3, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { AddCorr({n1,n2,n3},{},refs,pois); }
      void                    AddThreeGap(Int_t n1, Int_t n2, Int_t n3, Double_t gap, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { AddCorr({n1,n2,n3},{gap},refs,pois); }
      void                    AddFour(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { AddCorr({n1,n2,n3,n4},{},refs,pois); }
      void                    AddFourGap(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Double_t gap, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { AddCorr({n1,n2,n3,n4},{gap},refs,pois); }

      void                    SetFlowRFPsPt(Double_t min, Double_t max) { fFlowRFPsPtMin = min; fFlowRFPsPtMax = max; }
      void                    SetFlowPOIsPt(Double_t min, Double_t max, Int_t bins = 0) { fFlowPOIsPtMin = min; fFlowPOIsPtMax = max; fFlowPOIsPtBinNum = bins; }
      void                    SetFlowPOIsPtBins(std::vector<Double_t> bins, PartSpecies species) { fFlowPOIsPtBinEdges[species] = bins; }
      void                    SetFlowEta(Double_t max, Int_t bins = 0) { fFlowEtaMax = max; fFlowEtaBinNum = bins; }
      void                    SetFlowPhiBins(Int_t bins) { fFlowPhiBinNum = bins; }
      void                    SetV0sMassBins(Int_t bins) { fV0sNumBinsMass = bins; }
      void                    SetPhiMassBins(Int_t bins) { fPhiNumBinsMass = bins; }
      void                    SetFlowFillWeights(Bool_t weights = kTRUE) { fFlowFillWeights = weights; }
      void                    SetFlowFillAfterWeights(Bool_t weights = kTRUE) { fFlowFillAfterWeights = weights; }
      void                    SetUseWeigthsRunByRun(Bool_t bRunByRun = kTRUE) { fFlowRunByRunWeights = bRunByRun; }
      void                    SetUsePeriodWeigths(Bool_t weight = kTRUE) { fFlowPeriodWeights = weight; }
      void                    SetWeightsTag(TString tag) { fFlowWeightsTag = tag; }
      void                    SetUseWeights3D(Bool_t use = kTRUE) { fFlowUse3Dweights = use; }
      void                    SetApplyWeightsForReco(Bool_t apply = kTRUE) { fFlowWeightsApplyForReco = apply; }
      // events setters
      void                    SetCentrality(CentEst est, Int_t min = 0, Int_t max = 0, Int_t bins = 0) { fCentEstimator = est; fCentMin = min; fCentMax = max; fCentBinNum = bins; }
      void                    SetAddCentCut(CentEst est, Int_t min, Int_t max) { fCentEstimatorAdd = est; fCentMinAdd = min; fCentMaxAdd = max; }
      void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
      void					          SetPVtxZMax(Double_t z) { fPVtxCutZ = z; }
      void                    SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) { fVxMax = vx; fVyMax = vy; fVzMax = vz; }
      void                    SetRejectAddPileUp(Bool_t use = kTRUE) { fEventRejectAddPileUp = use; }
      void                    SetRejectAddPileUpESDTPCCut(Int_t cut) { fPileUpCutESDTPC = cut; }
      void                    SetPileUpCutCentrality(Int_t cut) { fPileUpCutCentrality = cut; }
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
      void                    SetPIDonlyForRefs(Bool_t use = kTRUE) { fPIDonlyForRefs = use; }
      // V0s setters
      void					  SetV0sOnFly(Bool_t onFly) { fCutV0sOnFly = onFly; }
      void					  SetV0sTPCRefit(Bool_t refit) { fCutV0srefitTPC = refit; }
      void					  SetV0sRejectKinks(Bool_t reject) { fCutV0srejectKinks = reject; }
      void                    SetV0sDaughterNumTPCClsMin(UShort_t cls) { fCutV0sDaughterNumTPCClsMin = cls; }
      void                    SetV0sDaughterNumTPCrossMin(UShort_t cls) { fCutV0sDaughterNumTPCCrossMin = cls; }
      void                    SetV0sDaughterNumTPCFindMin(UShort_t cls) { fCutV0sDaughterNumTPCFindMin = cls; }
      void                    SetV0sDaughterNumTPCClsPIDMin(UShort_t cls) { fCutV0sDaughterNumTPCClsPIDMin = cls; }
      void                    SetV0sDaughterRatioCrossFindMin(Double_t ratio) { fCutV0sDaughterRatioCrossFindMin = ratio; }
      void					  SetV0sUseCrossMassRejection(Bool_t reject) { fCutV0sCrossMassRejection = reject; }
      void					  SetV0sCrossMassCutK0s(Double_t mass) { fCutV0sCrossMassCutK0s = mass; }
      void					  SetV0sCrossMassCutLambda(Double_t mass) { fCutV0sCrossMassCutLambda = mass; }
      void					  SetV0sDCAPVMin(Double_t dca) { fCutV0sDCAtoPVMin = dca; }
      void					  SetV0sDCAPVMax(Double_t dca) { fCutV0sDCAtoPVMax = dca; }
      void					  SetV0sDCAPVzMax(Double_t dca) { fCutV0sDCAtoPVzMax = dca; }
      void            SetV0sDaughtersFilterBit(UInt_t filter) { fCutV0sDaughterFilterBit = filter; }
      void					  SetV0sDCADaughtersMin(Double_t dca) { fCutV0sDCADaughtersMin = dca; }
      void					  SetV0sDCADaughtersMax(Double_t dca) { fCutV0sDCADaughtersMax = dca; }
      void					  SetV0sDecayRadiusMin(Double_t radius) { fCutV0sDecayRadiusMin = radius; }
      void					  SetV0sDecayRadiusMax(Double_t radius) { fCutV0sDecayRadiusMax = radius; }
      void					  SetV0sDaughterEtaMax(Double_t eta) { fCutV0sDaughterEtaMax = eta; }
      void					  SetV0sDaughterPtMin(Double_t pt) { fCutV0sDaughterPtMin = pt; }
      void					  SetV0sDaughterPtMax(Double_t pt) { fCutV0sDaughterPtMax = pt; }
      void            SetV0sMotherRapMax(Double_t rap) { fCutV0sMotherRapMax = rap; }
      void					  SetV0sK0sInvMassMin(Double_t mass) { fCutV0sInvMassK0sMin = mass; }
      void					  SetV0sK0sInvMassMax(Double_t mass) { fCutV0sInvMassK0sMax = mass; }
      void					  SetV0sLambdaInvMassMin(Double_t mass) { fCutV0sInvMassLambdaMin = mass; }
      void					  SetV0sLambdaInvMassMax(Double_t mass) { fCutV0sInvMassLambdaMax = mass; }
      void					  SetV0sK0sCPAMin(Double_t cpa) { fCutV0sCPAK0sMin = cpa; }
      void					  SetV0sLambdaCPAMin(Double_t cpa) { fCutV0sCPALambdaMin = cpa; }
      void					  SetV0sK0sNumTauMax(Double_t nTau) { fCutV0sNumTauK0sMax = nTau; }
      void					  SetV0sLambdaNumTauMax(Double_t nTau) { fCutV0sNumTauLambdaMax = nTau; }
      void					  SetV0sK0sArmenterosAlphaMin(Double_t alpha) { fCutV0sArmenterosAlphaK0sMin = alpha; }
      void					  SetV0sLambdaArmenterosAlphaMax(Double_t alpha) { fCutV0sArmenterosAlphaLambdaMax = alpha; }
      void            SetV0sK0sPionNumTPCSigmaMax(Float_t nSigma) { fCutV0sK0sPionNumTPCSigmaMax = nSigma; }
      void            SetV0sLambdaPionNumTPCSigmaMax(Float_t nSigma) { fCutV0sLambdaPionNumTPCSigmaMax = nSigma; }
      void            SetV0sLambdaProtonNumTPCSigmaMax(Float_t nSigma) { fCutV0sLambdaProtonNumTPCSigmaMax = nSigma; }
      // phi setters
      void					  SetPhiInvMassMin(Double_t mass) { fCutPhiInvMassMin = mass; }
      void					  SetPhiInvMassMax(Double_t mass) { fCutPhiInvMassMax = mass; }

      //correlations related setters
      void            SetDoCorrelations(Bool_t use = kTRUE) { fCorrFill = use;}
      void            SetDEta(Int_t nBins, Double_t min, Double_t max) { fCorrDEtaBinNum = nBins; fCorrdEtaMin = min; fCorrdEtaMax = max; }
      void            SetDPhi(Int_t nBins, Double_t min, Double_t max) { fCorrDPhiBinNum = nBins; fCorrdPhiMin = min; fCorrdPhiMax = max; }
      void            SetDisablePtBinnedPool(Bool_t dis = kFALSE) {fUsePtBinnedEventPool = dis; }
      void            SetDisableFillMixedCorrelations(Bool_t fill = kFALSE) { fFillMixed = fill;}
      void            SetMaxPoolSize(Int_t size){fPoolSize = size; }
      Bool_t          FillCorrelations();
      Double_t        RangePhi(Double_t dPhi);

      AliEventCuts            fEventCuts; //

    private:
      static const Int_t      fPIDNumSpecies = 5; // Number of considered species for PID
      static const Int_t      fFlowNumHarmonicsMax = 25; // maximum harmonics length of flow vector array
      static const Int_t      fFlowNumWeightPowersMax = 17; // maximum weight power length of flow vector array
      static const Int_t      fFlowBinNumberEtaSlices = 32; // number of eta bin slices (for correlation study)

      const char*             GetSpeciesName(PartSpecies species) const;
      const char*             GetSpeciesName(Int_t species) const { return GetSpeciesName(PartSpecies(species)); }
      const char*             GetSpeciesLabel(PartSpecies species) const;
      const char*             GetSpeciesLabel(Int_t species) const { return GetSpeciesLabel(PartSpecies(species)); }
      const char*             GetEtaGapName(Double_t dEtaGap) const { return Form("%02.2g",10.0*dEtaGap); }

      Bool_t                  sortPt(const AliVParticle* t1, const AliVParticle* t2) { return (t1->Pt() < t2->Pt()); } // function for std::sort

      Bool_t                  InitializeTask(); // called once on beginning of task (within CreateUserObjects method)
      Bool_t                  LoadWeights(); // load weights histograms
      Bool_t                  FillFlowWeight(const AliVParticle* track, PartSpecies species) const; // fill distribution for per-particle flow weight
      Double_t                GetFlowWeight(const AliVParticle* track, PartSpecies species) const; // extract per-particle flow weight from input file
      const char*             ReturnPPperiod(const Int_t runNumber) const;
      void                    ListParameters() const; // list all task parameters
      void                    ClearVectors(); // properly clear all particle vectors
      void                    DumpTObjTable(const char* note, Option_t* opt = "") const; // add a printf statmenet given by note followed by gObjTable->Print() dump
      std::vector<Double_t>   MakeBinsVector(Int_t num, Double_t min, Double_t max); // transform fixed sized bins into an array of Double_t

      Bool_t                  IsEventSelected(); // event selection for Run 2 using AliEventCuts
      Bool_t                  IsMCEventSelected(); // event selection for simulated events
      Bool_t                  IsEventRejectedAddPileUp() const; // additional pile-up rejection for Run2 Pb-Pb
      Int_t                   GetSamplingIndex() const; // returns sampling index based on sampling selection (number of samples)
      Int_t                   GetCentralityIndex(CentEst est) const; // returns centrality index based centrality estimator or number of selected tracks
      const char*             GetCentEstimatorLabel(CentEst est) const; // returns mult/cent estimator string with label or 'n/a' if not available

      void                    ProcessMC() const; // processing MC generated particles
      void                    FilterCharged() const; // charged tracks filtering
      void                    FilterChargedMC() const; // charged tracks filtering
      void                    FilterPID() const; // pi,K,p filtering
      void                    FilterV0s() const; // K0s, Lambda, ALambda filtering
      void                    FilterPhi() const; // reconstruction and filtering of Phi meson candidates

      void                    CalculateCorrelations(const AliUniFlowCorrTask* task, PartSpecies species, Double_t dPt = -1.0, Double_t dMass = -1.0) const; // wrapper for correlations methods
      Bool_t                  ProcessCorrTask(const AliUniFlowCorrTask* task, const Int_t iTask, Bool_t doLowerOrder); // procesisng of AliUniFlowCorrTask
      Bool_t                  CalculateFlow(); // main (envelope) method for flow calculations in selected events

      AliAODMCParticle*       GetMCParticle(Int_t label) const; // find corresponding MC particle from fArrayMC depending of AOD track label
      Bool_t                  CheckMCPDG(const AliVParticle* track, const Int_t iPDGCode) const; // check if track has an associated MC particle which is the same species
      Bool_t                  CheckMCPDG(const AliVParticle* track, const PartSpecies species) const; // check if track has an associated MC particle which is the same species
      Bool_t                  CheckMCTruthReco(const PartSpecies species, const AliVParticle* track, const AliVParticle* daughterPos = nullptr, const AliVParticle* daughterNeg = nullptr) const; // check if Reco track has an associated MC particle which is the same species
      Double_t                PIDCorrection(const AliAODTrack *track, const PartSpecies species) const; //PID correction for 2018 data
      Double_t                GetRapidity(Double_t mass, Double_t Pt, Double_t Eta) const; // calculate particle / track rapidity
      Bool_t                  HasMass(PartSpecies spec) const { return (spec == kK0s || spec == kLambda || spec == kPhi); }
      Bool_t                  HasTrackPIDTPC(const AliAODTrack* track) const; // is TPC PID OK for this track ?
      Bool_t                  HasTrackPIDTOF(const AliAODTrack* track) const; // is TOF PID OK for this track ?
      Bool_t                  IsWithinRefs(const AliVParticle* track) const; // check if track is in (pt,eta) acceptance for Refs (used for refs selection & autocorelations)
      Bool_t                  IsWithinPOIs(const AliVParticle* track) const; // check if track is in (pt,eta) acceptance for POIs
      Bool_t                  IsChargedSelected(const AliAODTrack* track) const; // charged track selection
      PartSpecies             IsPIDSelected(AliVParticle* track) const; // PID tracks selections
      PartSpecies             IsPIDSelectedMC(AliVParticle* track) const; // PID tracks selections
      Bool_t                  IsV0Selected(const AliAODv0* v0) const; // general (common) V0 selection
      Bool_t                  IsV0aK0s(const AliAODv0* v0) const; // V0 selection: K0s specific
      Int_t                   IsV0aLambda(const AliAODv0* v0) const; // V0 selection: (A)Lambda specific
      AliPicoTrack*           MakeMother(const AliAODTrack* part1, const AliAODTrack* part2) const; // Combine two prongs into a mother particle stored in AliPicoTrack object
      void                    FillSparseCand(THnSparse* sparse, const AliVTrack* track) const; // Fill sparse histogram for inv. mass distribution of candidates (V0s,Phi)
      void                    FillQAEvents(QAindex iQAindex) const; // filling QA plots related to event selection
      void                    FillQARefs(QAindex iQAindex, const AliVParticle* track) const; // filling QA plots for RFPs selection
      void                    FillQACharged(QAindex iQAindex, AliVParticle* track) const; // filling QA plots for charged track selection
      void                    FillQAPID(QAindex iQAindex, AliVParticle* track, PartSpecies species) const; // filling pi,K,p QA histograms
      void                    FillQAV0s(QAindex iQAindex, const AliAODv0* v0, Bool_t bIsK0s = kTRUE, Int_t bIsLambda = 2) const; // filling QA plots for V0s candidates
      void                    FillQAPhi(QAindex iQAindex, const AliPicoTrack* part) const; // filling QA plots for V0s candidates

      // Flow related methods
      void                    FillRefsVectors(const AliUniFlowCorrTask* task, Double_t dGap); // fill flow vector Q with RFPs for reference flow
      Int_t                   FillPOIsVectors(const AliUniFlowCorrTask* task, Double_t dEtaGap, PartSpecies species, Int_t& indStart, Int_t& tracksInBin, Double_t dPtLow, Double_t dPtHigh, Double_t dMassLow = 0.0, Double_t dMassHigh = 0.0); // fill flow vectors p,q and s with POIs (for given species) for differential flow calculations
      Int_t                   FillPOIsVectorsCharged(const AliUniFlowCorrTask* task, Double_t dEtaGap, Double_t dPtLow, Double_t dPtHigh, std::array<Int_t, 4> &indexStart); // fill flow vectors p,q and s with POIs (for given species) for differential flow calculations for charged species with GF weights fix
      void                    ResetFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], Int_t maxHarm = 8, Int_t maxWeightPower = 4, Bool_t usePow = kFALSE, std::vector<Int_t> maxPowVec = {}); // set values to TComplex(0,0,0) for given array
      void                    ListFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]) const; // printf all values of given Flow vector array

      //dihadron corr related method
      void                    ResetFlowVectorQdih(TComplex (&array)[fFlowBinNumberEtaSlices][6][3], Int_t harm); // set values to TComplex(0,0,0) for given array
      void                    FillFlowQVectorsForDih(Double_t dWeight, Double_t dPhi, Double_t dEta, Int_t harm); // fill flow vector Q with RFPs for dihadron correlation study
      void                    CalculateDihCorr(const AliUniFlowCorrTask* task) const;

      TComplex                Q(Int_t n, Int_t p) const;
      TComplex                QGapPos(Int_t n, Int_t p) const;
      TComplex                QGapNeg(Int_t n, Int_t p) const;
      TComplex                QGapMid(Int_t n, Int_t p) const;
      TComplex                P(Int_t n, Int_t p) const;
      TComplex                PGapPos(Int_t n, Int_t p) const;
      TComplex                PGapNeg(Int_t n, Int_t p) const;
      TComplex                PGapMid(Int_t n, Int_t p) const;
      TComplex                S(Int_t n, Int_t p) const;
      TComplex                SGapPos(Int_t n, Int_t p) const;
      TComplex                SGapNeg(Int_t n, Int_t p) const;
      TComplex                SGapMid(Int_t n, Int_t p) const;

      TComplex                Two(Int_t n1, Int_t n2) const; // Two particle reference correlation calculations (no eta gap)
      TComplex                TwoGap(Int_t n1, Int_t n2) const; // Two particle reference correlation calculations (with eta gap)
      TComplex                TwoGap3sub(Int_t n1, Int_t n2, Int_t rf1Pos, Int_t rf2Pos) const; // Two particle reference correlation calculations (with 3 subevents)
      TComplex                TwoPos(Int_t n1, Int_t n2) const; /// Two particle reference correlation calculations (just from the positive eta)
      TComplex                TwoNeg(Int_t n1, Int_t n2) const; /// Two particle reference correlation calculations (just from the negative eta)
      TComplex                TwoMid(Int_t n1, Int_t n2) const; /// Two particle reference correlation calculations (just from the mid eta)
      TComplex                Three(Int_t n1, Int_t n2, Int_t n3) const; // Three particle reference correlation calculations (no eta gap)
      TComplex                ThreePos(Int_t n1, Int_t n2, Int_t n3) const; // Three particle reference correlation calculations (just from the positive eta)
      TComplex                ThreeNeg(Int_t n1, Int_t n2, Int_t n3) const; // Three particle reference correlation calculations (just from the negative eta)
      TComplex                Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle reference correlation calculations (no eta gap)
      TComplex                FourGap(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle reference correlation calculations (with eta gap)
      TComplex                FourPos(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle reference correlation calculations (just from the positive eta)
      TComplex                FourNeg(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle reference correlation calculations (just from the negative eta)
      TComplex                Four3sub(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t twoParCorrPosition) const; // Four particle reference correlation calculations (with 3 sub-events)
      TComplex                Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5) const; // Five particle reference correlation calculations (no eta gap)
      TComplex                FivePos(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5) const; // Five particle reference correlation calculations (just from the positive eta)
      TComplex                FiveNeg(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5) const; // Five particle reference correlation calculations (just from the negative eta)
      TComplex                Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6) const; // Six particle reference correlation calculations (no eta gap)
      TComplex                SixPos(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6) const; // Six particle reference correlation calculations (just from the positive eta)
      TComplex                SixNeg(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6) const; // Six particle reference correlation calculations (just from the negative eta)
      TComplex                SixGap(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6) const; // Six particle reference correlation calculations (with eta gap)
      TComplex                Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7) const; // Seven particle reference correlation calculations (no eta gap)
      TComplex                SevenPos(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7) const; // Seven particle reference correlation calculations (just from the positive eta)
      TComplex                SevenNeg(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7) const; // Seven particle reference correlation calculations (just from the negative eta)
      TComplex                Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8) const; // Eight particle reference correlation calculations (no eta gap)
      TComplex                EightGap(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8) const; // Eight particle reference correlation calculations (with eta gap)
      TComplex                EightPos(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8) const; // Eight particle reference correlation calculations (just from the positive eta)
      TComplex                EightNeg(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8) const; // Eight particle reference correlation calculations (just from the negative eta)

      TComplex                TwoDiff(Int_t n1, Int_t n2) const; // Two particle diff. correlation calculations (no eta gap)
      TComplex                TwoDiffGapPos(Int_t n1, Int_t n2) const; // Two particle diff. correlation calculations (with eta gap)
      TComplex                TwoDiffGapNeg(Int_t n1, Int_t n2) const; // Two particle diff. correlation calculations (with eta gap)
      TComplex                TwoDiffPos(Int_t n1, Int_t n2) const; // Two particle diff. correlation calculations (just from the positive eta)
      TComplex                TwoDiffNeg(Int_t n1, Int_t n2) const; // Two particle diff. correlation calculations (just from the negative eta)
      TComplex                TwoDiffMid(Int_t n1, Int_t n2) const; // Two particle diff. correlation calculations (just from the mid eta)
      TComplex                TwoDiffGap3sub(Int_t n1, Int_t n2, Int_t poiPos, Int_t rfPos) const; // Two particle diff. correlation calculations (with 3 subevents)
      TComplex                ThreeDiff(Int_t n1, Int_t n2, Int_t n3) const; // Three particle diff. correlation calculations (no eta gap)
      TComplex                ThreeDiffGapPos(Int_t n1, Int_t n2, Int_t n3) const; // Three particle diff. correlation calculations (with eta gap)
      TComplex                ThreeDiffGapNeg(Int_t n1, Int_t n2, Int_t n3) const; // Three particle diff. correlation calculations (with eta gap)
      TComplex                ThreeDiffPos(Int_t n1, Int_t n2, Int_t n3) const; // Three particle diff. correlation calculations (just from the positive eta)
      TComplex                ThreeDiffNeg(Int_t n1, Int_t n2, Int_t n3) const; // Three particle diff. correlation calculations (just from the negative eta)
      TComplex                FourDiff(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle reference correlation calculations (no eta gap)
      TComplex                FourDiffGapPos(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle diff. correlation calculations (with eta gap)
      TComplex                FourDiffGapNeg(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle diff. correlation calculations (with eta gap)
      TComplex                FourDiffPos(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle diff. correlation calculations (just from the positive eta)
      TComplex                FourDiffNeg(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const; // Four particle diff. correlation calculations (just from the negative eta)
      TComplex                FourDiff3sub(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t poiPosition, Int_t twoParCorrPosition) const; // Four particle diff. correlation calculations (with 3 sub-events)
      TComplex                SixDiffGapPos(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6) const; // Six particle reference correlation calculations (with eta gap)
      TComplex                SixDiffGapNeg(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6) const; // Six particle reference correlation calculations (with eta gap)
      TComplex                EightDiffGapPos(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8) const; // Eight particle reference correlation calculations (with eta gap)
      TComplex                EightDiffGapNeg(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8) const; // Eight particle reference correlation calculations (with eta gap)
      TComplex                Correlator(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0) const; // general formula

      // array lenghts & constants
      AliAODEvent*            fEventAOD; //! AOD event countainer
      // AliMCEvent*             fEventMC; //! MC event countainer
      AliVEvent*              fEvent; //! V event countainer
      Double_t                fPVz; // PV z-coordinate used for weights
      AliPIDResponse*         fPIDResponse; //! AliPIDResponse container
      AliPIDCombined*         fPIDCombined; //! AliPIDCombined container
      TList*                  fFlowWeightsList; //! list of weights from input file
      Bool_t                  fMC; // is running on mc?
      Bool_t                  fNeedPIDCorrection; // does data need PID correction?
      Bool_t                  fPIDCorrectionPhi; // special case for phi
      Bool_t                  fIs2018data; // is 2018 data?
      Bool_t                  fIsHMpp; // is high multiplicity pp? (different event selection)
      Bool_t                  fInit; // initialization check
      Bool_t                  fUseGeneralFormula; // [kFALSE] using of new formula
      Bool_t                  fFlowUsePIDWeights; // [kFALSE] using PID weights for Q vectors filling
      Bool_t                  fPIDonlyForRefs; // [kFALSE] for modification of GF
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
      TComplex                fFlowVecPmid[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecSpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecSneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecSmid[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation

      TComplex                fFlowVecQ[fFlowBinNumberEtaSlices][6][3]; // flow vector array for flow calculation in very narrow delta eta slices (for the correlation study)


      std::vector<AliUniFlowCorrTask*>  fVecCorrTask; //
      std::vector<AliVParticle*>* fVector[kUnknown]; //! container for selected Refs charged particles

      //cuts & selection: analysis
      RunMode                 fRunMode; // running mode (not grid related)
      AnalType                fAnalType; // analysis type: AOD / ESD / MC
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
      std::vector<Double_t>   fFlowPOIsPtBinEdges[kUnknown]; // pt bin edges for fixed pt bins per species
      Double_t                fFlowEtaMax; // [0.8] max eta acceptance for flow particles (RFPs & POIs)
      Int_t                   fFlowEtaBinNum; // [0] number of eta bins
      Int_t                   fFlowPhiBinNum; // [100] number of phi bins
      Int_t                   fPhiNumBinsMass; // number of InvMass bins for phi distribution
      Int_t                   fV0sNumBinsMass; // number of InvMass bins for V0s distribution
      Int_t                   fNumSamples; // [1] overall number of samples (from random sampling) used
      Bool_t                  fEtaCheckRFP; // [kFALSE] flag for doing analysis of FMPs for positive and negative eta separately
      Bool_t                  fFlowFillWeights; //[kFALSE] flag for filling weights
      Bool_t                  fFlowFillWeightsMultiD; //[kFALSE] flag for filling weights - multidimensional (phi, centrality, pT, eta, vz)
      Bool_t                  fFlowFillAfterWeights; //[kTRUE] flag for filling weights after NUA (only if fUseWeights is on)
      Bool_t                  fFlowUseWeights; //[kFALSE] flag for using the previously filled weights
      Bool_t                  fFlowUse3Dweights; // [kFALSE] flag for using 3D GF weights, if kFALSE, 2D weights are expected
      Bool_t                  fFlowRunByRunWeights; // [kTRUE] flag for using rub-by-run weigths from weigths file; if false, only one set of histrograms is provided
      Bool_t                  fFlowPeriodWeights; // [kFALSE] flag for using period-averaged weigths from weigths file; if false, only one set of histrograms is provided (average); for pp only
      Bool_t                  fFlowWeightsApplyForReco; //[kFALSE] flag for applying weights for Reco particles
      TString                 fFlowWeightsTag; // [""] tag with TList name for weights (used for systematics)
      // cuts & selection: correlations related
      AliEventPoolManager*    fEventPoolMgr; // event pool manager
      AliEventPool*           fPool; //
      TObjArray*              fSelectedTracks; //! tracks for mixing
      Bool_t                  fCorrUsingGF; // [kFALSE] fill correlations flag (but using GF and Q-cumulants)
      Bool_t                  fCorrFill; // [kFALSE] fill correlations flag
      Bool_t		              fFillMixed;		// [kTRUE] enable event mixing
      Bool_t		              fUsePtBinnedEventPool;		// [kTRUE] enable filling mixed events based on pT dependence
      Int_t                   fPoolSize; // [-1] maximum number of events, -1 means no limit
      Int_t  		              fMixingTracks;	// [5000] size of track buffer for event mixing
      Int_t  		              fMinEventsToMix;	// [5] min number of events for event mixing
      Int_t                   fCorrDEtaBinNum; // [32] number of dEta bins for correlations
      Int_t                   fCorrDPhiBinNum; // [72] number of dPhi bins for correlations
      Double_t                fCorrdEtaMin; // [-1.6] min of dEta bins for correlations
      Double_t                fCorrdEtaMax; // [1.6] max of dEta bins for correlations
      Double_t                fCorrdPhiMin; // [-pi/2] min of dEta bins for correlations
      Double_t                fCorrdPhiMax ; // [3/2 pi] max of dEta bins for correlations
      Double_t                fEtaSlicesArr[fFlowBinNumberEtaSlices+1]; // array with lower

      //cuts & selection: events
      ColSystem               fColSystem; // collisional system
      AliVEvent::EOfflineTriggerTypes    fTrigger; // physics selection trigger
      CentEst                 fCentEstimator; // [kV0A] multiplicity/centrality estimator as in AliMultSelection
      Int_t                   fCentMin; // [0] min range for centrality/multiplicity histos
      Int_t                   fCentMax; // [0] max range for centrality/multiplicity histos
      Int_t                   fCentBinNum; // [0] number of centrality bins
      CentEst                 fCentEstimatorAdd; // [kRFP] multiplicity/centrality estimator as in AliMultSelection
      Int_t                   fCentMinAdd; // [0] min range for centrality/multiplicity histos
      Int_t                   fCentMaxAdd; // [0] max range for centrality/multiplicity histos
      Double_t                fPVtxCutZ; // (cm) PV z cut
      Double_t                fVxMax; // vx max - MC
      Double_t                fVyMax; // vy max - MC
      Double_t                fVzMax; // vz max - MC
      Double_t                fImpactParameterMC; // impact parameter MC
      Bool_t                  fEventRejectAddPileUp; // additional pile-up rejection for Pb-Pb collisions in Run2 (17n, 15o)
      Int_t                   fPileUpCutESDTPC; // [500] additional pile-up rejection for Pb-Pb collisions in Run2 (15o)
      Int_t                   fPileUpCutCentrality; // [10] additional pile-up rejection for Pb-Pb collisions in Run2 (15o), upper limit for centrality (as this cut usually matter only for central collisions)
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
      Bool_t				  fCutV0sCrossMassRejection; // competing V0 rejection based on InvMass
      Double_t                fCutV0sCrossMassCutK0s; // [0.005] (GeV/c2) restricted vicinity of Lambda/ALambda inv. mass peak for K0s candidates
      Double_t                fCutV0sCrossMassCutLambda; // [0.020] (GeV/c2) restricted vicinity of K0s inv. mass peak for Lambda/ALambda candidates
      Double_t                fCutV0sDCAtoPVMin;   // (cm) min DCA of V0 daughter to PV
      Double_t				  fCutV0sDCAtoPVMax;	// (cm) max DCA of V0 daughter to PV
      Double_t                fCutV0sDCAtoPVzMax; // (cm) max DCA-z coordinate of V0 daughters to PV
      Double_t				  fCutV0sDCADaughtersMin;	// (cm) min DCA of V0 daughters among themselves
      Double_t				  fCutV0sDCADaughtersMax;	// (cm) max DCA of V0 daughters among themselves
      Double_t                fCutV0sDecayRadiusMin; // (cm) min distance of secondary vertex from z-axis in transverse plane
      Double_t				  fCutV0sDecayRadiusMax; // (cm) max distance of secondary vertex from z-axis in transverse plane
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
      Double_t				  fCutV0sArmenterosAlphaK0sMin; // (alpha) min Armenteros alpha for K0s
      Double_t                fCutV0sArmenterosAlphaLambdaMax; // (alpha) max Armenteros alpha for (Anti)Lambda
      Float_t                 fCutV0sK0sPionNumTPCSigmaMax; // (sigmaTPC) max number of TPC sigmas for kaon PID (K0s candidates)
      Float_t                 fCutV0sLambdaPionNumTPCSigmaMax;    // (sigmaTPC) max number of TPC sigma for pion PID (Lambda candidates)
      Float_t                 fCutV0sLambdaProtonNumTPCSigmaMax;    // (sigmaTPC) max number of TPC sigma for proton PID (Lambda candidates)
      // cuts & selection: phi
      Double_t                fCutPhiInvMassMin; // [0.99] (GeV/c2) min inv. mass window for selected phi candidates
      Double_t                fCutPhiInvMassMax; // [1.07] (GeV/c2) min inv. mass window for selected phi candidates

      // output lists
      TList*                  fQAEvents; //! events list
      TList*                  fQACharged; //! charged tracks list
      TList*                  fQAPID; //! pi,K,p list
      TList*                  fQAV0s; //! V0s candidates list
      TList*                  fQAPhi; //! Phi candidates list
      TList*                  fFlowWeights; //! list for flow weights
      TList*                  fListFlow[kUnknown]; //! flow lists
      TList*                  fListMC; //! list for MC

      // histograms & profiles

      // Flow
      THnSparseD*             fhsCandK0s; //! distribution of K0s candidates
      THnSparseD*             fhsCandLambda; //!  distribution of Lambda candidates
      THnSparseD*             fhsCandPhi; //!  distribution of Phi candidates
      THnSparseD*             fhsCandPhiBg; //!  distribution of Phi background

      TH2D*                   fh2Weights[kUnknown]; //! container for GF weights (phi,eta,pt) (2D)
      TH3D*                   fh3Weights[kUnknown]; //! container for GF weights (phi,eta,pt)
      TH2D*                   fh2AfterWeights[kUnknown]; //! distribution after applying GF weights - lightweight QA (phi)
      TH3D*                   fh3AfterWeights[kUnknown]; //! distribution after applying GF weights - full QA (phi,eta,pt)
      THnSparseD*             fhWeightsMultiD; //!  distribution of Phi background

      // Events
      TH2D*                   fhEventSampling; //! distribution of sampled events (based on randomly generated numbers)
      TH1D*                   fhEventCentrality; //! distribution of event centrality
      TH2D*                   fh2EventCentralityNumRefs; //! distribution of event centrality vs number of selected charged tracks
      TH1D*                   fhEventCounter; //! counter following event selection
      TH1D*                   fhV0Mamplitude; //! V0M amplitude (imporant for HM pp data)
      TH1D*                   fhV0MamplitudeRatio; //! V0M amplitude / <V0M> amplitude (imporant for HM pp data)
      TH2D*                   fh2V0MnCharged; //! V0M amplitude / <V0M> amplitude (imporant for HM pp data) vs. N_charged
      TH2D*                   fh2MeanMultRFP[10]; //! counter following RFP multiplicity (pT vs. mult.)
      TH2D*                   fh2MCip; //! impact parameter vs. Nch (for on-the-fly)
      // Charged
      TH1D*                   fhRefsMult; //!multiplicity distribution of selected RFPs
      TH1D*                   fhRefsPt; //! pt distribution of selected RFPs
      TH1D*                   fhRefsEta; //! pt distribution of selected RFPs
      TH1D*                   fhRefsPhi; //! pt distribution of selected RFPs
      TProfile*               fpRefsMult; //! <multiplicity>
      TH1D*                   fhChargedCounter; //! counter following charged track selection
      THnSparseD*             fh4CorrelationsSE[kUnknown]; //! eta phi distributin of the same event
      THnSparseD*             fh4CorrelationsME[kUnknown]; //! eta phi distributin for mixed events
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
      // MC
      TH2D*                   fh2MCPtEtaGen[kUnknown]; //! (pt,eta) dist for generated particles
      TH2D*                   fh2MCPtEtaReco[kUnknown]; //! (pt,eta) dist for reconstructed particles
      TH2D*                   fh2MCPtEtaRecoTrue[kUnknown]; //! (pt,eta) dist for reconstructed particles with matching generating particle (true)
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


      // QA: events
      TH1D*                   fhQAEventsPVz[QAindex::kNumQA]; //!
      TH1D*                   fhQAEventsNumContrPV[QAindex::kNumQA]; //!
      TH1D*                   fhQAEventsNumSPDContrPV[QAindex::kNumQA]; //!
      TH1D*                   fhQAEventsDistPVSPD[QAindex::kNumQA]; //!
      TH1D*                   fhQAEventsSPDresol[QAindex::kNumQA]; //!
      TH2D*                   fhQAEventsfMult32vsCentr; //!
      TH2D*                   fhQAEventsMult128vsCentr; //!
      TH2D*                   fhQAEventsfMultTPCvsTOF; //!
      TH2D*                   fhQAEventsfMultTPCvsESD; //!
      // QA: charged tracks
      TH1D*                   fhQAChargedMult[QAindex::kNumQA];       //! number of AOD charged tracks distribution
      TH1D*                   fhQAChargedPt[QAindex::kNumQA];         //! pT dist of charged tracks
      TH1D*                   fhQAChargedEta[QAindex::kNumQA];        //! eta dist of charged tracks
      TH1D*                   fhQAChargedPhi[QAindex::kNumQA];        //! phi dist of charged tracks
      TH1D*                   fhQAChargedCharge[QAindex::kNumQA];     //! charge dist of charged tracks
      TH1D*                   fhQAChargedNumTPCcls[QAindex::kNumQA];  //! dist of track number of TPC clusters
      TH1D*                   fhQAChargedDCAxy[QAindex::kNumQA];      //! dist of Charged DCA in transverse plane
      TH1D*                   fhQAChargedDCAz[QAindex::kNumQA];       //! dist of charged DCA in z coordinate
      // QA: PID tracks
      TH1D*                   fhQAPIDTPCstatus[QAindex::kNumQA];  //! based on AliPIDResponse::CheckPIDStatus();
      TH1D*                   fhQAPIDTOFstatus[QAindex::kNumQA];  //! based on AliPIDResponse::CheckPIDStatus();
      TH2D*                   fhQAPIDTPCdEdx[QAindex::kNumQA];    //! TPC PID information
      TH2D*                   fhQAPIDTOFbeta[QAindex::kNumQA];    //! TOF PID information
      TH3D*                   fh3QAPIDnSigmaTPCTOFPtPion[QAindex::kNumQA]; //! nSigma TPC vs nSigma TOF vs pt
      TH3D*                   fh3QAPIDnSigmaTPCTOFPtKaon[QAindex::kNumQA]; //! nSigma TPC vs nSigma TOF vs pt
      TH3D*                   fh3QAPIDnSigmaTPCTOFPtProton[QAindex::kNumQA]; //! nSigma TPC vs nSigma TOF vs pt
      // QA: V0s candidates
      TH1D*			  		  fhQAV0sMultK0s[QAindex::kNumQA];	//! number of K0s candidates
      TH1D*			  		  fhQAV0sMultLambda[QAindex::kNumQA];	//! number of Lambda candidates
      TH1D*			  		  fhQAV0sMultALambda[QAindex::kNumQA];	//! number of Anti-Lambda candidates
      TH1D*			  		  fhQAV0sRecoMethod[QAindex::kNumQA];	//! offline/online V0 reconstruction method
      TH1D*			  		  fhQAV0sDaughterTPCRefit[QAindex::kNumQA];	//! Daughters TPC refit true/false
      TH1D*			  		  fhQAV0sDaughterKinks[QAindex::kNumQA];	//! Daughters kinks true/false
      TH1D*                   fhQAV0sDaughterNumTPCCls[QAindex::kNumQA]; //! Daughter # of TPC findable clusters
      TH1D*                   fhQAV0sDaughterNumTPCFind[QAindex::kNumQA]; //! Daughter # of TPC clusters
      TH1D*                   fhQAV0sDaughterNumTPCCrossRows[QAindex::kNumQA]; //! Daughter # of TPC crossed rows
      TH1D*                   fhQAV0sDaughterTPCCrossFindRatio[QAindex::kNumQA]; //! Daughter # of TPC cross / # of TPC findable cls ratio
      TH1D*                   fhQAV0sDaughterNumTPCClsPID[QAindex::kNumQA]; //! Daughter # of TPC findable clusters used for PID
      TH1D*			  		  fhQAV0sDCAtoPV[QAindex::kNumQA];	//! V0 DCA to PV
      TH1D*			  		  fhQAV0sDCADaughters[QAindex::kNumQA];	//! DCA between V0 daughters
      TH1D*			  		  fhQAV0sDecayRadius[QAindex::kNumQA];	//! Distance between PV and Secondary vertex in transverse plane
      TH1D*                   fhQAV0sInvMassK0s[QAindex::kNumQA];    //! inv. mass dist of V0s (K0s mass hypothesis)
      TH1D*					  fhQAV0sInvMassLambda[QAindex::kNumQA];	//! inv. mass dist of V0s ((A)Lambda mass hypothesis)
      TH1D*                   fhQAV0sMotherPt[QAindex::kNumQA];  //! pT dist of V0s
      TH1D*					  fhQAV0sMotherPhi[QAindex::kNumQA];	//! azimuthal dist of V0s
      TH1D*                   fhQAV0sMotherEta[QAindex::kNumQA]; //! pseudorapidity dist of V0s
      TH1D*                   fhQAV0sMotherCharge[QAindex::kNumQA]; //! charge distribution of mothers
      TH1D*                   fhQAV0sMotherRapK0s[QAindex::kNumQA];  //! rapidity dist of V0s (K0s mass hypothesis)
      TH1D*                   fhQAV0sMotherRapLambda[QAindex::kNumQA]; //! rapidity dist of V0s (Lambda mass hypothesis)
      TH1D*                   fhQAV0sDaughterPt[QAindex::kNumQA];    //! pT dist of V0 daughters
      TH1D*					  fhQAV0sDaughterPhi[QAindex::kNumQA];	//! pT dist of V0 daughters
      TH1D*                   fhQAV0sDaughterEta[QAindex::kNumQA];   //! pseudorapidity dist of V0 daughters
      TH1D*                   fhQAV0sDaughterCharge[QAindex::kNumQA]; //! charge distribution of daughters
      TH1D*					  fhQAV0sDaughterTPCstatus[QAindex::kNumQA];	//! TPC dEdx vs p of K0s daughters
      TH1D*					  fhQAV0sDaughterTOFstatus[QAindex::kNumQA];	//! TPC dEdx vs p of K0s daughters
      TH2D*					  fhQAV0sDaughterTPCdEdxK0s[QAindex::kNumQA];	//! TPC dEdx vs p of K0s daughters
      TH2D*					  fhQAV0sDaughterNumSigmaPionK0s[QAindex::kNumQA];	//! Number of TPC sigmas (pion) vs mother pT of K0s daughters
      TH2D*					  fhQAV0sDaughterTPCdEdxLambda[QAindex::kNumQA];	//! TPC dEdx vs p of Lambda daughters
      TH2D*                   fhQAV0sDaughterNumSigmaPionLambda[QAindex::kNumQA];  //! number of TPC sigmas vs mother pT of pion (Lambda candidates)
      TH2D*                   fhQAV0sDaughterNumSigmaProtonLambda[QAindex::kNumQA];  //! number of TPC sigmas vs mother pT of proton (Lambda candidates)
      TH2D*                   fhQAV0sDaughterNumSigmaPionALambda[QAindex::kNumQA];  //! number of TPC sigmas vs mother pT of pion (Anti-Lambda candidates)
      TH2D*                   fhQAV0sDaughterNumSigmaProtonALambda[QAindex::kNumQA];  //! number of TPC sigmas vs mother pT of proton (Anti-Lambda candidates)
      TH1D*					  fhQAV0sCPAK0s[QAindex::kNumQA];	//! cosine of pointing angle of K0s candidates
      TH1D*					  fhQAV0sCPALambda[QAindex::kNumQA];	//! cosine of pointing angle of Lambda candidates
      TH1D*					  fhQAV0sNumTauK0s[QAindex::kNumQA];	//! number of c*tau of K0s candidates
      TH1D*					  fhQAV0sNumTauLambda[QAindex::kNumQA];	//! number of c*tau of Lambda candidates
      TH2D*				   	  fhQAV0sArmenterosK0s[QAindex::kNumQA];	//! Armenteros-Podolanski plot for K0s candidates
      TH2D*			  		  fhQAV0sArmenterosLambda[QAindex::kNumQA];	//! Armenteros-Podolanski plot for Lambda candidates
      TH2D*			  		  fhQAV0sArmenterosALambda[QAindex::kNumQA];	//! Armenteros-Podolanski plot for ALambda candidates

      ClassDef(AliAnalysisTaskUniFlow, 26);
};

#endif
