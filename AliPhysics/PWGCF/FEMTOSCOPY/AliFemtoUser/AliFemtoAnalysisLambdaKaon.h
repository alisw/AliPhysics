// AliFemtoAnalysisLambdaKaon.h

#ifndef ALIFEMTOANALYSISLAMBDAKAON_H
#define ALIFEMTOANALYSISLAMBDAKAON_H

#include "AliFemtoVertexMultAnalysis.h"

#include "AliFemtoBasicEventCut.h"
#include "AliFemtoEventCutEstimators.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorV0.h"
#include "AliFemtoCutMonitorV0CosPointingAngle.h"
#include "AliFemtoCutMonitorXi.h"
#include "AliFemtoCutMonitorEventVertex.h"

#include "AliFemtoBasicTrackCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoAODTrackCut.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorParticlePID.h"

#include "AliFemtoXiTrackCut.h"
#include "AliFemtoXiTrackPairCut.h"
#include "AliFemtoXiV0PairCut.h"

#include "AliFemtoV0PairCut.h"
#include "AliFemtoV0TrackPairCut.h"
#include "AliFemtoPairOriginMonitor.h"

#include "AliFemtoCorrFctnKStar.h"
#include "AliFemtoAvgSepCorrFctn.h"

#include "AliFemtoDummyPairCut.h"
#include "AliFemtoV0PurityBgdEstimator.h"

#include "AliFemtoNSigmaFilter.h"
#include "AliFemtoV0TrackCutNSigmaFilter.h"
#include "AliFemtoXiTrackCutNSigmaFilter.h"
#include "AliFemtoESDTrackCutNSigmaFilter.h"

#include "AliFemtoCutMonitorEventPartCollSize.h"


#include "AliFemtoModelWeightGeneratorBasicLednicky.h"
#include "AliFemtoModelCorrFctnKStarFull.h"

#include "AliFemtoCorrFctnDirectYlm.h"

#include <string>
#include <iostream>
#include <stdio.h>
#include <typeinfo>


class AliFemtoAnalysisLambdaKaon : public AliFemtoVertexMultAnalysis {

public:
  enum AnalysisType {kLamK0=0, kALamK0=1, 
                     kLamKchP=2, kALamKchP=3, kLamKchM=4, kALamKchM=5, 
                     kLamLam=6, kALamALam=7, kLamALam=8, 
                     kLamPiP=9, kALamPiP=10, kLamPiM=11, kALamPiM=12, 
                     kXiKchP=13, kAXiKchP=14, kXiKchM=15, kAXiKchM=16,
                     kXiK0=17, kAXiK0=18,
                     kProtPiM=19, kAProtPiP=20, kPiPPiM=21};

  enum GeneralAnalysisType {kV0V0=0, kV0Track=1, kXiTrack=2, kXiV0=3, kTrackTrack=4};

  enum ParticlePDGType {kPDGProt   = 2212,  kPDGAntiProt = -2212, 
		        kPDGPiP    = 211,   kPDGPiM      = -211, 
                        kPDGK0     = 310,
                        kPDGKchP   = 321,   kPDGKchM     = -321,
		        kPDGLam    = 3122,  kPDGALam     = -3122,
		        kPDGSigma  = 3212,  kPDGASigma   = -3212,
		        kPDGXiC    = 3312,  kPDGAXiC     = -3312,
		        kPDGXi0    = 3322,  kPDGAXi0     = -3322,
		        kPDGOmega  = 3334,  kPDGAOmega   = -3334,
                        kPDGNull      = 0                        };

  enum GeneralParticleType {kV0=0, kTrack=1, kCascade=2};

struct AnalysisParams
{
  unsigned int nBinsVertex;
  double minVertex,
         maxVertex;

  unsigned int nBinsMult;
  double minMult, 
         maxMult;

  bool binEventsInRP;  //bin events in reaction plane angle (in addition to vertex z-position and multiplicity)
  int nBinsRP;

  AnalysisType analysisType;
  GeneralAnalysisType generalAnalysisType;

  unsigned int nEventsToMix;
  unsigned int minCollectionSize;

  bool verbose;

  bool implementAvgSepCuts;
  bool implementPairCutsOnlyInKStarCfs;

  bool writePairKinematics;
  bool isMCRun;
  bool isMBAnalysis;
  bool buildMultHist;
  bool buildmTBinned;
  bool implementVertexCorrections;
  bool removeMisidentifiedMCParticles;
  bool setV0SharedDaughterCut;

  bool addCutMonitors;
  bool monitorEvCutPassOnly;
  bool monitorPart1CutPassOnly;
  bool monitorPart2CutPassOnly;
  bool monitorPairCutPassOnly;
  bool useMCWeightGenerator;

  bool buildSphericalHarmonics;
  bool useLCMSforSH;

  bool buildIndmTBinnedCfs;
};

struct EventCutParams
{
  double minCentrality,
         maxCentrality;

  double minMult,
         maxMult;

  double minVertexZ,
         maxVertexZ;

  bool verboseMode;

  int centEstMethod;
  int multEstMethod;
};

struct V0CutParams
{
  ParticlePDGType particlePDGType;
  GeneralParticleType generalParticleType;

  int v0Type;  //0=kLambda, 1=kAntiLambda, 2=kK0s

  double mass;
  double minInvariantMass,
         maxInvariantMass;

  bool useLooseInvMassCut;
  double minLooseInvMass,
         maxLooseInvMass;

  int nBinsPurity;
  double minPurityMass,
         maxPurityMass;

  bool useCustomFilter;

  bool removeMisID;
  double minInvMassReject,
         maxInvMassReject;

  bool useSimpleMisID;
  bool buildMisIDHistograms;
  bool useCustomMisID;

  double eta;
  double minPt,
         maxPt;
  bool onFlyStatus;
  double maxDcaV0;
  double minCosPointingAngle;
  double maxV0DecayLength;

  double etaDaughters;
  double minPtPosDaughter,
         maxPtPosDaughter;
  double minPtNegDaughter,
         maxPtNegDaughter;
  int minTPCnclsDaughters;
  double maxDcaV0Daughters;
  double minPosDaughterToPrimVertex,
         minNegDaughterToPrimVertex;

  double radiusV0Min,
         radiusV0Max;

  bool ignoreOnFlyStatus;

  bool buildCosPointingAnglewParentInfo;
};

struct ESDCutParams
{
  ParticlePDGType particlePDGType;
  GeneralParticleType generalParticleType;

  double minPidProbPion,
         maxPidProbPion;
  double minPidProbMuon,
         maxPidProbMuon;
  double minPidProbKaon,
         maxPidProbKaon;
  double minPidProbProton,
         maxPidProbProton;
  int mostProbable;
  int charge;
  double mass;

  double minPt,
         maxPt;
  double eta;
  int minTPCncls;

  bool removeKinks;
  bool setLabel;
  double maxITSChiNdof;
  double maxTPCChiNdof;
  double maxSigmaToVertex;
  double minImpactXY;
  double maxImpactXY;
  double maxImpactZ;

  bool useCustomFilter;
  bool useCustomMisID;
  bool useElectronRejection;
  bool useCustomElectronRejection;
  bool useIsProbableElectronMethod;
  bool usePionRejection;
};

struct XiCutParams
{
  ParticlePDGType particlePDGType;
  GeneralParticleType generalParticleType;

  int charge;
  int xiType;
  double minPt,
         maxPt;
  double eta;
  double mass;
  double minInvariantMass,
         maxInvariantMass;

  double maxDecayLengthXi;
  double minCosPointingAngleXi;
  double minCosPointingAngleV0toXi;
  double maxDcaXi;
  double maxDcaXiDaughters;

  double minDcaXiBac;
  double etaBac;
  int minTPCnclsBac;
  double minPtBac,
         maxPtBac;

  double radiusXiMin,
         radiusXiMax;

  int v0Type;
  double minDcaV0;
  double minInvMassV0,
         maxInvMassV0;
  double minCosPointingAngleV0;  //this is V0 to primary vertex, not terribly useful for Xi analysis
  double etaV0;
  double minPtV0,
         maxPtV0;
  bool onFlyStatusV0;
  double maxV0DecayLength;
  double minV0PosDaughterToPrimVertex,
         minV0NegDaughterToPrimVertex;
  double maxDcaV0Daughters;
  double etaV0Daughters;
  double minPtPosV0Daughter,
         maxPtPosV0Daughter;
  double minPtNegV0Daughter,
         maxPtNegV0Daughter;

  double radiusV0Min,
         radiusV0Max;

  int minTPCnclsV0Daughters;

  bool useCustomV0Filter;
  bool useCustomV0MisID;
  bool useCustomBacPionFilter;
  bool useCustomBacPionMisID;

  bool ignoreOnFlyStatusV0;
};

struct PairCutParams
{
  bool removeSameLabel;
  double shareQualityMax;
  double shareFractionMax;
  bool tpcOnly;

  double tpcExitSepMinimum;
  double tpcEntranceSepMinimum;

  double minAvgSepPosPos,  // Set these for V0-V0 pair; V0-AntiV0 and AntiV0-AntiV0 cases
         minAvgSepPosNeg,  // will automatically be handled by AliFemtoAnalysisLambdaKaon::CreateV0PairCut
         minAvgSepNegPos,
         minAvgSepNegNeg;

  double minAvgSepTrackPos,  // Set these for V0-PosTrack; V0-NegTrack, AntiV0-Pos, etc. cases
         minAvgSepTrackNeg;  // will automatically be handled by AliFemtoAnalysisLambdaKaon::CreateV0TrackPairCut

  double minAvgSepTrackBacPion;

  double minAvgSepBacPos;
  double minAvgSepBacNeg;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------
//**********************************************************************************************************************************************************
//----------------------------------------------------------------------------------------------------------------------------------------------------------

  AliFemtoAnalysisLambdaKaon(AnalysisType aAnalysisType, unsigned int binsVertex, double minVertex, double maxVertex, unsigned int binsMult, double minMult, double maxMult, bool aIsMCRun, bool aImplementAvgSepCuts, bool aWritePairKinematics=false, TString aDirNameModifier="");

  AliFemtoAnalysisLambdaKaon(AnalysisParams &aAnParams, EventCutParams &aEvCutParams, TString aDirNameModifier="");
  AliFemtoAnalysisLambdaKaon(AnalysisParams &aAnParams, EventCutParams &aEvCutParams, PairCutParams &aPairCutParams, V0CutParams &aV0CutParams1, V0CutParams &aV0CutParams2, TString aDirNameModifier="");
  AliFemtoAnalysisLambdaKaon(AnalysisParams &aAnParams, EventCutParams &aEvCutParams, PairCutParams &aPairCutParams, V0CutParams &aV0CutParams1, ESDCutParams &aESDCutParams2, TString aDirNameModifier="");
  AliFemtoAnalysisLambdaKaon(AnalysisParams &aAnParams, EventCutParams &aEvCutParams, PairCutParams &aPairCutParams, XiCutParams &aXiCutParams1, ESDCutParams &aESDCutParams2, TString aDirNameModifier="");
  AliFemtoAnalysisLambdaKaon(AnalysisParams &aAnParams, EventCutParams &aEvCutParams, PairCutParams &aPairCutParams, XiCutParams &aXiCutParams1, V0CutParams &aV0CutParams1, TString aDirNameModifier="");
  AliFemtoAnalysisLambdaKaon(AnalysisParams &aAnParams, EventCutParams &aEvCutParams, PairCutParams &aPairCutParams, ESDCutParams &aESDCutParams1, ESDCutParams &aESDCutParams2, TString aDirNameModifier="");

    //Since I am using rdr->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality), 
      // in AliFemtoEventReaderAOD.cxx this causes tEvent->SetNormalizedMult(lrint(10*cent->GetCentralityPercentile("V0A"))), i.e. fNormalizedMult in [0,1000]
      // Therefore, since fNormalizedMult is presumably >= -1, in AliFemtoEvent.cxx the call UncorrectedNumberOfPrimaries returns fNormalizedMult
    //LONG STORY SHORT:  the inputs for multiplicity in the above are actually for 10*centrality (i.e. 0-100 for 0-10% centrality)
    //Note:  fNormalizedMult is typically in range [0,1000] (as can be seen in AliFemtoEventReaderAOD.cxx).  This appears true when SetUseMultiplicity is set to kCentrality, kCentralityV0A, kCentralityV0C, kCentralityZNA, kCentralityZNC, kCentralityCL1, kCentralityCL0, kCentralityTRK, kCentralityTKL, kCentralityCND, kCentralityNPA, kCentralityFMD.
      // fNormalizedMult WILL NOT be [0,1000] when SetUseMultiplicity is set to kGlobalCount, kReference, kTPCOnlyRef, and kVZERO 

  
  AliFemtoAnalysisLambdaKaon(const AliFemtoAnalysisLambdaKaon& TheOriginalAnalysis);
  AliFemtoAnalysisLambdaKaon& operator=(const AliFemtoAnalysisLambdaKaon& TheOriginalAnalysis);
  virtual ~AliFemtoAnalysisLambdaKaon();

  virtual void ProcessEvent(const AliFemtoEvent* ProcessThisEvent);  //will add fMultHist to the process event of AliFemtoVertexMultAnalysis
  virtual TList* GetOutputList();


  void SetParticleTypes(AnalysisType aAnType);

  AliFemtoBasicEventCut* CreateBasicEventCut(EventCutParams &aEvCutParams);
  AliFemtoEventCutEstimators* CreateEventCutEstimators(EventCutParams &aEvCutParams);

  void AddCustomV0SelectionFilters(ParticlePDGType aV0Type, AliFemtoV0TrackCutNSigmaFilter* aCut);
  void AddCustomV0RejectionFilters(ParticlePDGType aV0Type, AliFemtoV0TrackCutNSigmaFilter* aCut);
  AliFemtoV0TrackCutNSigmaFilter* CreateV0Cut(V0CutParams &aCutParams);

  void AddCustomESDSelectionFilters(ParticlePDGType aESDType, AliFemtoESDTrackCutNSigmaFilter* aCut);
  void AddCustomESDRejectionFilters(ParticlePDGType aESDType, AliFemtoESDTrackCutNSigmaFilter* aCut);
  AliFemtoESDTrackCutNSigmaFilter* CreateESDCut(ESDCutParams &aCutParams);

  void AddCustomXiSelectionFilters(ParticlePDGType aXiType, AliFemtoXiTrackCutNSigmaFilter* aCut);
  void AddCustomXiV0RejectionFilters(ParticlePDGType aXiType, AliFemtoXiTrackCutNSigmaFilter* aCut);
  AliFemtoXiTrackCutNSigmaFilter* CreateXiCut(XiCutParams &aCutParams);

  AliFemtoV0PairCut* CreateV0PairCut(PairCutParams &aPairCutParams);
  AliFemtoV0TrackPairCut* CreateV0TrackPairCut(PairCutParams &aPairCutParams);
  AliFemtoXiTrackPairCut* CreateXiTrackPairCut(PairCutParams &aPairCutParams);
  AliFemtoXiV0PairCut* CreateXiV0PairCut(PairCutParams &aPairCutParams);

  AliFemtoCorrFctnKStar* CreateCorrFctnKStar(const char* name, unsigned int bins, double min, double max);
  AliFemtoAvgSepCorrFctn* CreateAvgSepCorrFctn(const char* name, unsigned int bins, double min, double max);
  AliFemtoModelCorrFctnKStarFull* CreateModelCorrFctnKStarFull(const char* name, unsigned int bins, double min, double max);    //TODO check that enum to int is working
  AliFemtoV0PurityBgdEstimator* CreateV0PurityBgdEstimator();
  AliFemtoCorrFctnDirectYlm* CreateCorrFctnDirectYlm(const char* name, int maxl, unsigned int bins, double min, double max, int useLCMS);

  void AddCutMonitors(AliFemtoEventCut* aEventCut, AliFemtoParticleCut* aPartCut1, AliFemtoParticleCut* aPartCut2, AliFemtoPairCut* aPairCut);
  void SetAnalysis(AliFemtoEventCut* aEventCut, AliFemtoParticleCut* aPartCut1, AliFemtoParticleCut* aPartCut2, AliFemtoPairCut* aPairCut);
  void SetMultHist(const char* name, int aNbins=30, double aMin=0., double aMax=3000);

  //------Builders for default cut objects
  static AnalysisParams DefaultAnalysisParams();
  static EventCutParams DefaultEventCutParams();

  static V0CutParams DefaultLambdaCutParams();
  static V0CutParams DefaultAntiLambdaCutParams();
  static V0CutParams DefaultK0ShortCutParams();

  static ESDCutParams DefaultKchCutParams(int aCharge);
  static ESDCutParams DefaultPiCutParams(int aCharge);

  static XiCutParams DefaultXiCutParams();
  static XiCutParams DefaultAXiCutParams();

  static PairCutParams DefaultPairParams();

  static ESDCutParams LambdaPurityPiCutParams(int aCharge);
  static ESDCutParams K0ShortPurityPiCutParams(int aCharge);
  static ESDCutParams LambdaPurityProtonCutParams(int aCharge);


  //----------------------------------------






protected:
  static const char* const fAnalysisTags[];

  AnalysisParams fAnalysisParams;
  AnalysisType fAnalysisType;
  GeneralAnalysisType fGeneralAnalysisType;
  ParticlePDGType fParticlePDGType1, fParticlePDGType2;
  GeneralParticleType fGeneralParticleType1, fGeneralParticleType2;
  TString fOutputName;		      /* name given to output directory for specific analysis*/
  TH1F* fMultHist;			      //histogram of event multiplicities to ensure event cuts are properly implemented
  bool fImplementAvgSepCuts;		      //Self-explanatory, set to kTRUE when I want Avg Sep cuts implemented
  bool fImplementPairCutsOnlyInKStarCfs;  //This will allow me to have an unbiased sample for fNumerator_RotatePar2 in AliFemtoCorrFctnKStar
                                          //As implied by name, pair cut will not be implemented in AliFemtoSimpleAnalysis, but only in the
                                          //AliFemtoCorrFctnKStar objects
  bool fWritePairKinematics;
  bool fIsMCRun;
  bool fIsMBAnalysis;
  bool fBuildMultHist;
  bool fBuildmTBinned;
  bool fBuildCosPointingAnglewParentInfo;

  double fMinCent, fMaxCent;

  AliFemtoCorrFctnCollection* fCollectionOfCfs;

  //----------------------------------------

  AliFemtoCorrFctnKStar *KStarCf;
  AliFemtoAvgSepCorrFctn *AvgSepCf;

  AliFemtoModelCorrFctnKStarFull *KStarModelCfs;





#ifdef __ROOT__
  ClassDef(AliFemtoAnalysisLambdaKaon, 0)
#endif

};



#endif
