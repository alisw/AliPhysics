/*
 * AliAnalysisTaskThreeBodyFemtoAOD.h
 *
 *  Created on: May 13, 2019
 *      Author: Laura Serksnyte 
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYFEMTOAOD_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYFEMTOAOD_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamControlSample.h"

class AliAnalysisTaskThreeBodyFemtoAOD : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskThreeBodyFemtoAOD();
  AliAnalysisTaskThreeBodyFemtoAOD(const char* name, bool isMC, bool triggerOn);
  virtual ~AliAnalysisTaskThreeBodyFemtoAOD();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  
  void FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);
  void FillPairDistributionPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, TH1F* sameEventDistributionPL, int mult, TH2F* hist2d);
  void FillPairDistributionPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, TH1F* sameEventDistributionPP, TH2F **fEventTripletPhiThetaArray,  AliFemtoDreamCollConfig Config);
  void SetMixedEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer);
  void FillTripletDistributionMEPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);
  // test different mixing
  void SetMixedEventOnlyPLambdaTEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void FillTripletDistributionMEPPTEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);
  void FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray,  AliFemtoDreamCollConfig Config);
  // test different mixing 2
  void SetMixedEventOnlyPPLambdaTEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  // Create triplets like (pp)l (lp)p
  void FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);
  // Add the close pair cut
  bool DeltaEtaDeltaPhi(AliFemtoDreamBasePart &part1,AliFemtoDreamBasePart &part2, bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist, AliFemtoDreamCollConfig Config);
  bool MyLovely3BodyTrigger(AliAODEvent *evt ,  bool isMC, std::vector<int> PDGCodes);
  double CalculatePPLTriggerQ3Min(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, std::vector<int> PDGCodes );


  void SetRunTaskLightWeight(bool light) {
    fisLightWeight = light;
  }
  void SetEventCuts(AliFemtoDreamEventCuts* evtCuts) {
    fEventCuts = evtCuts;
  }
  void SetProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fProton = trkCuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fAntiProton = trkCuts;
  }
  void Setv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {
    fLambda = v0Cuts;
  }
  void SetAntiv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {
    fAntiLambda = v0Cuts;
  }

  void SetEventCutsTrigger(AliFemtoDreamEventCuts* evtCutsTrigger) {
    fEventCutsTrigger = evtCutsTrigger;
  }
  void SetProtonCutsTrigger(AliFemtoDreamTrackCuts* trkCutsTrigger) {
    fTrackCutsTrigger = trkCutsTrigger;
  }
  void SetAntiProtonCutsTrigger(AliFemtoDreamTrackCuts* trkCutsTrigger) {
    fAntiTrackCutTrigger = trkCutsTrigger;
  }
  void Setv0CutsTrigger(AliFemtoDreamv0Cuts* v0CutsTrigger) {
    fv0CutsTrigger = v0CutsTrigger;
  }
  void SetAntiv0CutsTrigger(AliFemtoDreamv0Cuts* v0CutsTrigger) {
    fAntiv0CutsTrigger = v0CutsTrigger;
  }

  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig=config;
  }  
  void SetRunThreeBodyHistograms(bool RunThreeBodyHistos) {
    fRunThreeBody=RunThreeBodyHistos;
  }
  void SetTriggerOn(bool triggerOn) {
    fTriggerOn=triggerOn;
  }  
  void SetIsMC(bool isMCLocal) {
    fIsMC=isMCLocal;
  }

  void SetTriggerOnSample(bool triggerOnSample) {
    fTriggerOnSample=triggerOnSample;
  } 

  void SetRunppp0ppL1(bool runpppppL) {
    fRunppp0ppL1=runpppppL;
  } 


  void SetQ3Limit(float Q3Limit) {
    fQ3Limit = Q3Limit;
  }  
  void SetQ3LimitSample(float Q3LimitSample) {
    fQ3LimitSample = Q3LimitSample;
  } 
  void SetQ3LimitSample2(float Q3LimitSample2) {
    fQ3LimitSample2 = Q3LimitSample2;
  } 
  void SetQ3LimitSampleFraction(float Q3LimitSampleFraction) {
    fQ3LimitSampleFraction = Q3LimitSampleFraction;
  } 
  void SetQ3LimitSampleFraction2(float Q3LimitSampleFraction2) {
    fQ3LimitSampleFraction2 = Q3LimitSampleFraction2;
  } 

  void SetQ3LimitFraction(float Q3LimitFraction) {
    fQ3LimitFraction = Q3LimitFraction;
  }
  
  static TLorentzVector RelativePairMomentum(TLorentzVector &PartOne, TLorentzVector &PartTwo);
 private:
  AliAnalysisTaskThreeBodyFemtoAOD(const AliAnalysisTaskThreeBodyFemtoAOD &task);
  AliAnalysisTaskThreeBodyFemtoAOD &operator=(const AliAnalysisTaskThreeBodyFemtoAOD &task);
  bool fisLightWeight;//
  AliFemtoDreamEvent* fEvent;//!
  AliFemtoDreamEventCuts* fEventCuts;//
  TList* fEvtList;//!
  AliFemtoDreamTrack* fTrack;//!
  AliFemtoDreamTrackCuts* fProton;//
  TList* fProtonList;//!
  TList* fProtonMCList;//!
  AliFemtoDreamTrackCuts* fAntiProton;//
  TList* fAntiProtonList;//!
  TList* fAntiProtonMCList;//!
  AliFemtoDreamv0* fv0;//!
  AliFemtoDreamv0Cuts* fLambda;//
  TList* fLambdaList;
  TList* fLambdaMCList;
  AliFemtoDreamv0Cuts* fAntiLambda;//
  TList* fAntiLambdaList;
  TList* fAntiLambdaMCList;
  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  TList *fResults;//!
  // Three particles same event
  TList *fResultsThreeBody;//!
  TList *fSameEvent;//! 
  TList *fMixedEvent;//!  
  TList *fSameEventMult;//! 
  TList *fMixedEventMult;//!  
  TList *fSameEventPhiTheta;//! 
  TList *fMixedEventPhiTheta;//!  
  TList *fOtherHistos;//! 

  TRandom3* fRandomGen;
  
  bool fRunThreeBody;
  TH1F **fSameEventTripletArray;
  TH2F **fSameEventTripletMultArray;
  TH2F **fSameEventTripletPhiThetaArray;
  // Three particles mixed events
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainer;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerTEST;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerTESTppL;
  TH1F **fMixedEventTripletArray;
  TH2F **fMixedEventTripletMultArray;
  TH2F **fMixedEventTripletPhiThetaArray;
  // Three particles trigger studies
  bool fTriggerOn;
  bool fTriggerOnSample;
  bool fIsMC;
  bool fRunppp0ppL1;

  float fQ3Limit;
  float fQ3LimitSample;
  float fQ3LimitSample2;
  float fQ3LimitSampleFraction;
  float fQ3LimitSampleFraction2;
  float fQ3LimitFraction;

  TH1F* fRejectedEvents;
  TH1F* fAcceptedEvents;
  TH1F* fAcceptedEventsButNoPPL;
  TH1F* fTripletsPerCollision;
  TH1F* fWhichSample;
  TH1F* fSameEventOnlyLowestQ3;
  TH1F* fSameEventOnlyLowestQ3Anti;
  int ftotalTripletCount;
  AliFemtoDreamEventCuts* fEventCutsTrigger;//
  TList* fEventCutsTriggerList;//!
  AliFemtoDreamTrackCuts* fTrackCutsTrigger;//
  TList* fTrackCutsTriggerList;//!
  AliFemtoDreamTrackCuts* fAntiTrackCutTrigger;//
  TList* fAntiTrackCutTriggerList;//!
  AliFemtoDreamv0Cuts* fv0CutsTrigger;//
  TList* fv0CutsTriggerList;
  AliFemtoDreamv0Cuts* fAntiv0CutsTrigger;//
  TList* fAntiv0CutsTriggerList;
  ////////////////////////77
  TList *fResultsQA;//!
  AliFemtoDreamControlSample *fSample;   //!
  TList *fResultsSample;//!
  TList *fResultsSampleQA;//!
  int fTrackBufferSize;//
  AliAODTrack **fGTI;  //!
  ClassDef(AliAnalysisTaskThreeBodyFemtoAOD,2)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskThreeBodyFemtoAOD_H_ */


