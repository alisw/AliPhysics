/*
 * AliAnalysisTaskFemtoProtonPion.h
 *
 *  Created on: 11 Mar 2022
 *  Author: Lesch Marcel
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOPROTONPION_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOPROTONPION_H_
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEvent.h"
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
#include "AliVEvent.h"

class AliAnalysisTaskFemtoProtonPion : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFemtoProtonPion();
  AliAnalysisTaskFemtoProtonPion(const char *name, bool isMC);
  virtual ~AliAnalysisTaskFemtoProtonPion();
  void InitHistograms(AliFemtoDreamTrackCuts *trkCuts, TString trkCutsName, TString MCName);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {};
  virtual void InitializeArrays();
  void FillPairDistributionSE(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies, std::vector<int> PDGCodes, int mult, bool DoClosePairRejection, TH1F* hist, TH2F* hist2d, TH1F* HistInvMass, TH1F* HistInvMassMCResonance, TH2F **SameEventPhiTheta_OneDimensional, int CombinationNumber, AliFemtoDreamCollConfig Config); 
  void FillPairDistributionSEAncestors(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies, std::vector<int> PDGCodes, int mult, bool DoClosePairRejection, TH1F* hist, TH2F* hist2d, TH1F* HistInvMass, TH2F **SameEventPhiTheta_OneDimensional, TH1F **histAncestor, TH2F **hist2dAncestor, TH2F **SameEventPhiTheta_OneDimensionalAncestor, int CombinationNumber, AliFemtoDreamCollConfig Config);
  void FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME, std::vector<int> PDGCodes, int mult, bool DoClosePairRejection, TH1F* hist, TH2F* hist2d, TH1F* HistInvMass, TH2F **EventPhiThetaArray, int CombinationNumber, AliFemtoDreamCollConfig Config);

  bool DeltaEtaDeltaPhi(int species1, int species2, AliFemtoDreamBasePart &part1, AliFemtoDreamBasePart &part2, int part1PDGcode,int part2PDGcode, unsigned int PairDaughterIdentifier, TH2F* beforeHist,TH2F* afterHist, AliFemtoDreamCollConfig Config, double RelativeMomentum);
  void SetMixedEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *PartContainer);
  bool CommonAncestors(AliFemtoDreamBasePart& part1, AliFemtoDreamBasePart& part2); //Stolen from AliFemtoDreamHigherPairMath
  bool CommonMotherResonance(AliFemtoDreamBasePart& part1, AliFemtoDreamBasePart& part2); //check if two particles are from a certain resonance 
  bool IsResonance(int PDG); 
  bool PassedMCKineCuts(AliAODMCParticle *mcPart); 

  double GetQOutLCMS(const TLorentzVector Particle1, const TLorentzVector Particle2);
  double GetQSideLCMS(const TLorentzVector Particle1, const TLorentzVector Particle2);
  double GetQLongLCMS(const TLorentzVector Particle1, const TLorentzVector Particle2);

  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts = evtCuts;};
  void SetTrackCutsPion(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsPion = trkCuts;};
  void SetTrackCutsAntiPion(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiPion = trkCuts;};
  void SetTrackCutsProton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsProton = trkCuts;};
  void SetTrackCutsAntiProton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiProton = trkCuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig = config;};
  void SetRunTaskLightWeight(bool light) {fisLightWeight = light;};
  void SetDoPairCleaning(bool DoPairCleaning){fDoPairCleaning = DoPairCleaning;};
  void SetCombinationInput(string CombinationInput){fCombinationInput = CombinationInput;};
  void SetClosePairRejectionInput(string ClosePairRejectionInput){fClosePairRejectionInput = ClosePairRejectionInput;};
  void SetNameTagInput(string NameTagInput){fNameTagInput = NameTagInput;};
  void SetDoOfficialFemto(bool DoOfficialFemto){fDoOfficialFemto = DoOfficialFemto;};
  void SetDoOwnFemto(bool DoOwnFemto){fDoOwnFemto = DoOwnFemto;};
  void SetDoThreeDFemto(bool DoThreeDFemto){fDoThreeDFemto = DoThreeDFemto;};
  void SetRunPlotMult(bool RunPlotMult){fRunPlotMult = RunPlotMult;};
  void SetRunPlotPhiTheta(bool RunPlotPhiTheta){fRunPlotPhiTheta = RunPlotPhiTheta;};
  void SetDoAncestors(bool DoAncestors){fDoAncestors = DoAncestors;};
  void SetRemoveMCResonances(bool RemoveMCResonances, bool RemoveMCResonanceDaughters){fRemoveMCResonances = RemoveMCResonances; fRemoveMCResonanceDaughters = RemoveMCResonanceDaughters;};
  void SetDoInvMassPlot(bool DoInvMassPlot){fDoInvMassPlot = DoInvMassPlot;};
  void SetDoResonanceLorentzFactor(bool DoResonanceLorentzFactor){fDoResonanceLorentzFactor = DoResonanceLorentzFactor;};
  void SetDoReco(bool DoReco){fRecoDist = DoReco;}; 
  void SetDoKine(bool DoKine){fKineDist = DoKine;}; 

  private:
  AliAnalysisTaskFemtoProtonPion(const AliAnalysisTaskFemtoProtonPion &task);
  AliAnalysisTaskFemtoProtonPion &operator=(const AliAnalysisTaskFemtoProtonPion &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fisLightWeight;                      //
  int fTrackBufferSize;                     //
  bool fIsMC;                               //
  bool fDoPairCleaning;			   //

  int fCombinations[10][2]; //
  string fCombinationInput; //
  string fNameTags[4]; //
  string fNameTagInput; //
  bool fClosePairRejection[10]; //
  string fClosePairRejectionInput; //

  bool fDoOfficialFemto; // 
  bool fDoOwnFemto; // Do own looping and calculations
  bool fDoThreeDFemto; //  Three dimensional femtoscopy in own looping
  bool fRunPlotMult; //
  bool fRunPlotPhiTheta; // 
  bool fDoAncestors; //
  bool fRemoveMCResonances; //
  bool fRemoveMCResonanceDaughters; //
  bool fDoInvMassPlot; //
  bool fDoResonanceLorentzFactor; //

  //For k* dependent efficiency 
  bool fKineDist; //
  bool fRecoDist; //

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
  AliAODTrack** fGTI;           //!
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
  TH1F **fSameEvent_OneDimensional;
  TH2F **fSameEventMult_OneDimensional;
  TH1F **fSameEvent_OneDimensional_Ancestors;
  TH2F **fSameEventMult_OneDimensional_Ancestors;
  TH1F **fSameEvent_InvMass; 
  TH1F **fSameEvent_InvMass_MCResonance;

  TList *fMixedEvent_List_OneDimensional;
  TH1F **fMixedEvent_OneDimensional;
  TH2F **fMixedEventMult_OneDimensional;
  TH1F **fMixedEvent_InvMass; 

  TList *fSameEvent_List_ThreeDimensional;
  TH2F **fSameEvent_ThreeDimensional;
  TList *fMixedEvent_List_ThreeDimensional;
  TH2F **fMixedEvent_ThreeDimensional;

  TList *fSameEventDeltaEtaDeltaPhi_List;
  TH2F **fSameEventPhiTheta;
  TH2F **fSameEventPhiTheta_Ancestors;

  TList *fMixedEventDeltaEtaDeltaPhi_List;
  TH2F **fMixedEventPhiTheta;

  TH2F *fResonanceLorentzFactor; 
  TH2F *fInvMassResonancesMCTruth; 

  TH1F **fpTKineOrReco; 
  TH1F **fEtaKineOrReco;
  TH1F **fPhiKineOrReco;

  ClassDef(AliAnalysisTaskFemtoProtonPion, 4) 
};
#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOPROTONPION_H_ */
