/*
 * AliAnalysisTaskThreeBodyFemtoAODPionProton.h
 *
 *  Created on: May 13, 2019
 *      Author: Laura Serksnyte 
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYFEMTOAODPIONPROTON_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYFEMTOAODPIONPROTON_H_
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

class AliAnalysisTaskThreeBodyFemtoAODPionProton : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskThreeBodyFemtoAODPionProton();
  AliAnalysisTaskThreeBodyFemtoAODPionProton(const char* name, bool isMC, bool triggerOn);
  virtual ~AliAnalysisTaskThreeBodyFemtoAODPionProton();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  
  void FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);
  void FillTripletDistributionPluspT(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH1F* histpTDist);
  void SetMixedEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer);
  void FillTripletDistributionMEPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);
    // Create triplets like (pp)l (lp)p
  void FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);
  // Add the close pair cut
  bool DeltaEtaDeltaPhi(AliFemtoDreamBasePart &part1,AliFemtoDreamBasePart &part2, bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist, AliFemtoDreamCollConfig Config);
  void MyLovely3BodyTrigger(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector ,  bool isMC, std::vector<int> PDGCodes,  TH1F** histAccepted, TH1F** histRejected, AliFemtoDreamCollConfig Config, float mult, TH2F** fEventTripletPhiThetaArray, int phiEtaHistNo);
  double CalculatePPLTriggerQ3Min(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, std::vector<int> PDGCodes, AliFemtoDreamCollConfig Config, TH2F** fEventTripletPhiThetaArray, int phiEtaHistNo );


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
  void SetPionCuts(AliFemtoDreamTrackCuts* pionCuts) {
    fPion = pionCuts;
  }
  void SetAntiPionCuts(AliFemtoDreamTrackCuts* antiptionCuts) {
    fAntiPion = antiptionCuts;
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

  void SetQ3Limit(float Q3Limit) {
    fQ3Limit = Q3Limit;
  }

  
  static TLorentzVector RelativePairMomentum(TLorentzVector &PartOne, TLorentzVector &PartTwo);
 private:
  AliAnalysisTaskThreeBodyFemtoAODPionProton(const AliAnalysisTaskThreeBodyFemtoAODPionProton &task);
  AliAnalysisTaskThreeBodyFemtoAODPionProton &operator=(const AliAnalysisTaskThreeBodyFemtoAODPionProton &task);
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
  AliFemtoDreamTrack* fTrackPion;//!
  AliFemtoDreamTrackCuts* fPion;//
  TList* fPionList;
  TList* fPionMCList;
  AliFemtoDreamTrackCuts* fAntiPion;//
  TList* fAntiPionList;
  TList* fAntiPionMCList;
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
  TList *fRejectedParticles;
  TList *fAcceptedParticles;
  TList *fPtOfSmallQ3Pion;
  TList *fOtherHistos;//!
  bool fRunThreeBody;
  TH1F **fPtOfSmallQ3PionArray;
  TH1F **fRejectedParticlesArray;
  TH1F **fAcceptedParticlesArray;
  TH1F **fSameEventTripletArray;
  TH2F **fSameEventTripletMultArray;
  TH2F **fSameEventTripletPhiThetaArray;
  // Three particles mixed events
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainer;
  TH1F **fMixedEventTripletArray;
  TH2F **fMixedEventTripletMultArray;
  TH2F **fMixedEventTripletPhiThetaArray;
  // Three particles trigger studies
  bool fTriggerOn;
  bool fIsMC;
  float fQ3Limit;
  TH1F* fAllEvents;
  ////////////////////////77
  TList *fResultsQA;//!
  AliFemtoDreamControlSample *fSample;   //!
  TList *fResultsSample;//!
  TList *fResultsSampleQA;//!
  int fTrackBufferSize;//
  AliAODTrack **fGTI;  //!
  ClassDef(AliAnalysisTaskThreeBodyFemtoAODPionProton,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskThreeBodyFemtoAODPionProton_H_ */


