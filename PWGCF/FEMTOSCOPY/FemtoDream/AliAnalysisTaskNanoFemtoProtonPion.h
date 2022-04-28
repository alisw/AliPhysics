/*
 * AliAnalysisTaskNanoFemtoProtonPion.h
 *
 *  Created on: 11 Mar 2022
 *  Author: Lesch Marcel
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOFEMTOPROTONPION_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOFEMTOPROTONPION_H_
#include "AliAnalysisTaskSE.h"
#include "AliConvEventCuts.h"
#include "AliEventCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamControlSample.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TChain.h"
#include "AliFemtoDreamBaseDump.h"
class AliVParticle;
class AliVTrack;

class AliAnalysisTaskNanoFemtoProtonPion : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoFemtoProtonPion();
  AliAnalysisTaskNanoFemtoProtonPion(const char *name, bool isMC);
  virtual ~AliAnalysisTaskNanoFemtoProtonPion();
  void InitHistograms(AliFemtoDreamTrackCuts *trkCuts, TString trkCutsName, TString MCName);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {};
  virtual void InitializeArrays();
  void FillPairDistributionSE(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies, std::vector<int> PDGCodes, int mult, TH1F* hist, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int CombinationNumber, AliFemtoDreamCollConfig Config); 
  void FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME, std::vector<int> PDGCodes, int mult, TH1F* hist, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int CombinationNumber, AliFemtoDreamCollConfig Config);

  bool DeltaEtaDeltaPhi(int species1, int species2, AliFemtoDreamBasePart &part1, AliFemtoDreamBasePart &part2, int part1PDGcode,int part2PDGcode, unsigned int PairDaughterIdentifier, TH2F* beforeHist,TH2F* afterHist, AliFemtoDreamCollConfig Config, double RelativeMomentum);
  void SetMixedEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *PartContainer);

  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts = evtCuts;};
  void SetTrackCutsPion(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsPion = trkCuts;};
  void SetTrackCutsAntiPion(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiPion = trkCuts;};
  void SetTrackCutsProton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsProton = trkCuts;};
  void SetTrackCutsAntiProton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiProton = trkCuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig = config;};
  void SetRunTaskLightWeight(bool light) {fisLightWeight = light;};
  void SetDoPairCleaning(bool DoPairCleaning){fDoPairCleaning = DoPairCleaning;};
  void SetCombinationInput(string CombinationInput){fCombinationInput = CombinationInput;};
  void SetNameTagInput(string NameTagInput){fNameTagInput = NameTagInput;};
  void SetDoOfficialFemto(bool DoOfficialFemto){fDoOfficialFemto = DoOfficialFemto;};
  void SetDoThreeDFemto(bool DoThreeDFemto){fDoThreeDFemto = DoThreeDFemto;};
  void SetRunPlotMult(bool RunPlotMult){fRunPlotMult = RunPlotMult;};
  void SetRunPlotPhiTheta(bool RunPlotPhiTheta){fRunPlotPhiTheta = RunPlotPhiTheta;};
  void SetDoClosePairRejection(bool DoClosePairRejection){fDoClosePairRejection = DoClosePairRejection;};

  private:
  AliAnalysisTaskNanoFemtoProtonPion(const AliAnalysisTaskNanoFemtoProtonPion &task);
  AliAnalysisTaskNanoFemtoProtonPion &operator=(const AliAnalysisTaskNanoFemtoProtonPion &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  bool fisLightWeight;                      //
  int fTrackBufferSize;                     //
  bool fIsMC;                               //
  bool fDoPairCleaning;			   //

  int fCombinations[10][2]; //
  string fCombinationInput; //
  string fNameTags[4]; //
  string fNameTagInput; //

  bool fDoOfficialFemto; // 
  bool fDoThreeDFemto; //  
  bool fRunPlotMult; //
  bool fRunPlotPhiTheta; // 
  bool fDoClosePairRejection; //

  //GANESHA add needed histograms!! 

  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //
  AliFemtoDreamTrackCuts *fTrackCutsPion;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiPion;  //
  AliFemtoDreamTrackCuts *fTrackCutsProton;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiProton;  //
  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliVTrack** fGTI;           //!
  TList *fEvtList;//!
  TList *fProtonList;//!
  TList* fProtonMCList;//!
  TList *fAntiProtonList;//!
  TList* fAntiProtonMCList;//!
  TList *fPionList;//!
  TList* fPionMCList;//!
  TList *fAntiPionList;//!
  TList* fAntiPionMCList;//!
  TList *fResults;                          //!
  TList *fResultsQA;                        //!
  TList *fResultsThreeDFemto; //!

  // Particles mixed events container
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainer;

  TList *fSameEvent_List_OneDimensional; 
  TList *fSameEventDeltaEtaDeltaPhi_List_OneDimensional;
  TH1F **fSameEvent_OneDimensional;
  TH2F **fSameEventMult_OneDimensional;
  TH2F **fSameEventPhiTheta_OneDimensional;

  TList *fMixedEvent_List_OneDimensional;
  TList *fMixedEventDeltaEtaDeltaPhi_List_OneDimensional;
  TH1F **fMixedEvent_OneDimensional;
  TH2F **fMixedEventMult_OneDimensional;
  TH2F **fMixedEventTripletPhiTheta_OneDimensional;

  ClassDef(AliAnalysisTaskNanoFemtoProtonPion, 2) //GANESHA 1->2
};
#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOFEMTOPROTONPION_H_ */
