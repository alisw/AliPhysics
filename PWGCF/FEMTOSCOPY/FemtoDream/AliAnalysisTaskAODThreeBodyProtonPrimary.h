/*
 * AliAnalysisTaskAODThreeBodyProtonPrimary.h
 *
 *  Created on: Oct 13, 2023
 *      Authors: Marcel Lesch
 *      Based on AliAnalysisTaskThreeBodyFemtoAOD.h from Laura Serksnyte 
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKAODTHREEBODYPROTONPRIMARY_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKAODTHREEBODYPROTONPRIMARY_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamControlSample.h"
#include "AliAODTrack.h"

class AliAnalysisTaskAODThreeBodyProtonPrimary : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskAODThreeBodyProtonPrimary();
  AliAnalysisTaskAODThreeBodyProtonPrimary(const char* name, bool isMC);
  virtual ~AliAnalysisTaskAODThreeBodyProtonPrimary();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  
  void FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* hist2d12, TH2F* hist2d23, TH2F* hist2d31, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* InvMass12, TH2F* InvMass23, TH2F* InvMass31, TH2F* histmTQ312, TH2F* histmTQ323, TH2F* histmTQ331, TH2F* InvMassVsmT12, TH2F* InvMassVsmT23, TH2F* InvMassVsmT31, TH2F* Res, TH2F* ResAll);
  
  void SetMixedEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer);

  void FillTripletDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* hist2d12, TH2F* hist2d23, TH2F* hist2d31, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* histmTQ312, TH2F* histmTQ323, TH2F* histmTQ331, TH2F* InvMassVsmT12, TH2F* InvMassVsmT23, TH2F* InvMassVsmT31, TH2F* Projector, TH2F* Res, TH2F* ResAll);

  // test different mixing
  void SetMixedEventPPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void SetMixedEventPPPrim(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void SetMixedEventPPAPrim(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void SetMixedEventPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void SetMixedEventPPrim(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void SetMixedEventPAPrim(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);

  // Create triplets like (pp)prim (prim p)p
  void FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* hist2d12, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* InvMass12, TH2F* InvMass23, TH2F* InvMass31, TH2F* histmTQ312, TH2F* histmTQ323, TH2F* histmTQ331, TH2F* InvMassVsmT12, TH2F* InvMassVsmT23, TH2F* InvMassVsmT31);
  
  // Add the close pair cut
  bool DeltaEtaDeltaPhi(int species1, int species2, AliFemtoDreamBasePart &part1,AliFemtoDreamBasePart &part2,  int part1PDGcode, int part2PDGcode, bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist, AliFemtoDreamCollConfig Config, double Q3);
  bool DeltaEtaDeltaPhi(int species1, int species2, AliFemtoDreamBasePart &part1,AliFemtoDreamBasePart &part2,  int part1PDGcode, int part2PDGcode, bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist, AliFemtoDreamCollConfig Config);
  void FillPairInvMass(AliFemtoDreamBasePart &part1, AliFemtoDreamBasePart &part2, AliFemtoDreamBasePart &part3, TH2F* hist, float Q3) ;
  void FillPDGPairInvMass(AliFemtoDreamBasePart &part1, float massPart1, AliFemtoDreamBasePart &part2, float massPart2, AliFemtoDreamBasePart &part3, float massPart3, TH2F* hist, float Q3);

  void MomentumResolution( TH2F* histAll, TH2F* hist, AliFemtoDreamBasePart& part1, int PDGPart1, float mass1, AliFemtoDreamBasePart& part2, int PDGPart2, float mass2, AliFemtoDreamBasePart& part3, int PDGPart3, float mass3, float Q3Reconstructed); 
  //Two Body--------------------------------------------------------------------
  void FillPairDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* histDeltaPhi, TH2F* histLowDeltaPhi, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);
  void FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* histDeltaPhi, TH2F* histLowDeltaPhi, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);
  void FillPairTransverseMass(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies, int secondSpecies, TH1F* hist1, std::vector<int> PDGCodes, TH2F* hist2); 

  bool CommonAncestors(AliFemtoDreamBasePart& part1, AliFemtoDreamBasePart& part2); //Stolen from AliFemtoDreamHigherPairMath
  bool CommonMotherResonance(AliFemtoDreamBasePart* part1, AliFemtoDreamBasePart* part2); //check if two particles are from a certain resonance 
  bool IsResonance(int PDG); 

  //For mT Calculation 
  float GetmT(TLorentzVector &Part1, float mass1, TLorentzVector &Part2, float mass2);
  //Setters --------------------------------------------------------------------
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
  void SetPrimaryCuts(AliFemtoDreamTrackCuts* primCuts) {
    fPrimary = primCuts;
  } 
  void SetAntiPrimaryCuts(AliFemtoDreamTrackCuts* primCuts) {
    fAntiPrimary = primCuts;
  } 
  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig=config;
  }


  void SetDoOnlyThreeBody(bool DoOnlyThreeBody) {
    fDoOnlyThreeBody=DoOnlyThreeBody;
  }
  void SetRunOfficialTwoBody(bool RunOfficialTwoBody){ 
    fRunOfficialTwoBody=RunOfficialTwoBody;
  }
  void SetStandardMixing(bool StandardMixing){
    fStandardMixing = StandardMixing;
  }

  void SetClosePairRejectionForAll(bool ClosePairRejectionForAll) {
    fClosePairRejectionForAll=ClosePairRejectionForAll;
  }
  void SetturnoffClosePairRejectionCompletely(bool turnoffClosePairRejectionCompletely) {
    fturnoffClosePairRejectionCompletely=turnoffClosePairRejectionCompletely;
  }
  void SetDeltaPhiMaxPP(double DeltaPhiMaxPP){
    fDeltaPhiMaxPP = DeltaPhiMaxPP;
  }      
  void SetDeltaEtaMaxPP(double DeltaEtaMaxPP){
    fDeltaEtaMaxPP = DeltaEtaMaxPP;
  }      
  void SetDeltaPhiMaxPPrim(double DeltaPhiMaxPPrim){
    fDeltaPhiMaxPPrim = DeltaPhiMaxPPrim;
  }      
  void SetDeltaEtaMaxPPrim(double DeltaEtaMaxPPrim){
    fDeltaEtaMaxPPrim = DeltaEtaMaxPPrim;
  }   
  void SetDeltaPhiMaxPAPrim(double DeltaPhiMaxPAPrim){
    fDeltaPhiMaxPAPrim = DeltaPhiMaxPAPrim;
  }      
  void SetDeltaEtaMaxPAPrim(double DeltaEtaMaxPAPrim){
    fDeltaEtaMaxPAPrim = DeltaEtaMaxPAPrim;
  }    

  void SetQ3LimitForDeltaPhiDeltaEta(double Q3LimitForDeltaPhiDeltaEta) {
    fQ3LimitForDeltaPhiDeltaEta=Q3LimitForDeltaPhiDeltaEta;
  }
  void SetQ3cutValue(double Q3cutValue){
    fQ3cutValue = Q3cutValue;
  }
  void SetQ3MinValue(double Q3MinValue){
    fQ3MinValue = Q3MinValue;
  }

  void SetRunThreeBodyHistograms(bool RunThreeBodyHistos) {
    fRunThreeBody=RunThreeBodyHistos;
  }
  void SetRunPlotInvMass(bool RunPlotInvMass) {
    fRunPlotInvMass=RunPlotInvMass;
  }
  void SetRunPlotPhiTheta(bool RunPlotPhiTheta) {
    fRunPlotPhiTheta=RunPlotPhiTheta;
  }
  void SetRunPlotMult(bool RunPlotMult) {
    fRunPlotMult=RunPlotMult;
  }
  void SetRunPairMultThreeBody(bool doThis) {
    fRunPairMultThreeBody=doThis;
  }
  
  void SetRunmTPlots(bool RunmTPlots){
    fRunmTPlots = RunmTPlots;
  }
  void SetRunPlotInvMassVSmT(bool RunPlotInvMassVSmT){
    fRunPlotInvMassVSmT = RunPlotInvMassVSmT;
  }

  void SetMCAndReso(bool isMC, bool removeReso){
    fIsMC = isMC;
    fRemoveMCResonances = removeReso; 
  }

  void SetRunProjector(bool runIt){
     fRunProjector = runIt; 
  }

  void SetGetMomentumResolution(bool doIt){
    fGetMomentumResolution = doIt; 
  }

  static TLorentzVector RelativePairMomentum(TLorentzVector &PartOne, TLorentzVector &PartTwo);
  static float BoostOneParticle(TLorentzVector &PartOne, TLorentzVector &PartTwo, TLorentzVector &PartThree);

  private:
  AliAnalysisTaskAODThreeBodyProtonPrimary(const AliAnalysisTaskAODThreeBodyProtonPrimary &task);
  AliAnalysisTaskAODThreeBodyProtonPrimary &operator=(const AliAnalysisTaskAODThreeBodyProtonPrimary &task);
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
  AliFemtoDreamTrackCuts* fPrimary;// 
  TList* fPrimaryList;//!
  TList* fPrimaryMCList;//!
  AliFemtoDreamTrackCuts* fAntiPrimary;//
  TList* fAntiPrimaryList;//!
  TList* fAntiPrimaryMCList;//!
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
  TList *fSameEventmT;//!
  TList *fMixedEventmT;//!
  TList *fSameEventPhiTheta_SamePair;//!
  TList *fMixedEventPhiTheta_SamePair;//!
  TList *fSameEventPhiTheta_DifferentPair;//!
  TList *fMixedEventPhiTheta_DifferentPair;//!
  TList *fInvMassList;//!
  TList *fP1Histos;//!
  TList *fKHistos;//!
  TList *fInvMassVSmTList;//!

  bool fRunThreeBody;
  bool fRunPlotInvMass;
  bool fRunPlotPhiTheta;
  bool fRunPlotMult;
  bool fRunPairMultThreeBody;
  bool fRunOfficialTwoBody;
  bool fRunPlotInvMassVSmT;

  bool fClosePairRejectionForAll; 
  bool fturnoffClosePairRejectionCompletely; 

  double fQ3LimitForDeltaPhiDeltaEta;
  double fDeltaPhiMaxPP;
  double fDeltaEtaMaxPP;
  double fDeltaPhiMaxPPrim;
  double fDeltaEtaMaxPPrim;
  double fDeltaPhiMaxPAPrim;
  double fDeltaEtaMaxPAPrim;
  double fQ3cutValue;
  double fQ3MinValue;

  bool fDoOnlyThreeBody; //if kTRUE: 3 Body analysis, else 2 Body
  bool fStandardMixing;
  bool fRunmTPlots;

  bool fIsMC;
  bool fRemoveMCResonances; 

  bool fRunProjector; 
  bool fGetMomentumResolution; 

  //0: Proton, 1: AntiProton, 2: Primary, 3: AntiPrimary
  // int fTripletCombinations[6][3] = {{0,0,0}, {1,1,1}, {0,0,2}, {1,1,3}, {0,0,3}, {1,1,2}}; //triplet combinations that go into same event (xxx) and mixed event (x)(x)(x)
  int fTripletCombinations[6][3] = {{0,0,0}, {1,1,1}, {2,0,0}, {3,1,1}, {3,0,0}, {2,1,1}}; //triplet combinations that go into same event (xxx) and mixed event (x)(x)(x)
  int fTripletContainerID[6] = {0, 0, 1, 1, 2, 2}; //0: PPP Container, 1: PPPrim Container, 2: PPAPrim Container
  
  int fSameMixedCominations[10][3] = {{0,0,0}, {1,1,1}, {0,0,2}, {0,2,0}, {1,1,3}, {1,3,1}, {0,0,3}, {0,3,0}, {1,1,2}, {1,2,1}}; //triplet combinations that go into same-mixed event (xx)(x) 
  int fSameMixedContainerID[10] = {0, 0, 1, 1, 1, 1, 2, 2, 2, 2}; //0: PPP Container, 1: PPPrim Container, 2: PPAPrim Container
  

  TString fParticleNames[4] = {"P", "AP", "Prim", "APrim"}; //name of particles
  int fPairCombinations[6][2] = {{0,0}, {1,1}, {0,2}, {1,3}, {0,3}, {1,2}}; //pair combinations that go into same event (xx) and mixed event (x)(x)
  int fPairContainerID[6] = {0, 0, 1, 1, 2, 2};

  TH1F **fSameEventTripletArray; //Q3 Distribution for (xxy) and (xx)(y) events
  TH2F **fSameEventTripletMultArray; //Q3 Distribution vs mult bin for (xxy) and (xx)(y) events

  TH2F **fSameEventTripletMultArray12; //k* of pair (ij) vs mult bin 
  TH2F **fSameEventTripletMultArray23;
  TH2F **fSameEventTripletMultArray31;

  TH2F **fSameEventTripletmTArray12; //Q3 vs mT of a pair for (xxy) and (xx)(y) events. In case of  (xx)(y): Array(ij), i,j = 1,2 same event, 3 mixed event
  TH2F **fSameEventTripletmTArray23; 
  TH2F **fSameEventTripletmTArray31; 

  TH2F **fSameEventTripletPhiThetaArray_SamePair;
  TH2F **fSameEventTripletPhiThetaArray_DifferentPair;

  TH1F **fSameEventPairArray_TwoBody; //k* Distribution for (xy)
  TH2F **fSameEventMultArray_TwoBody; //k* Distribution vs mult bin for (xy)
  TH2F **fSameEventDphiArray_TwoBody;
  TH2F **fSameEventPairMultArray_TwoBody;
  TH2F **fSameEventPairPhiThetaArray_TwoBody;
  TH1F **fPairTranverseMass_TwoBody; 
  TH2F **fPairTranverseMassVSkstar_TwoBody; 

  TH1F *fMixingConfig; 
  
  // Three particles mixed events
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainer;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerPPP;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerPPPrim;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerPPAPrim;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerPP;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerPPrim;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerPAPrim;
  std::vector<std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>>> fVectPartContainers;

  TH1F **fMixedEventTripletArray; //Q3 Distribution for (x)(x)(y)
  TH2F **fMixedEventTripletMultArray; //Q3 Distribution vs mult bin for (x)(x)(y)

  TH2F **fMixedEventTripletMultArray12; //k* of pair (ij) vs mult bin (mixed event)
  TH2F **fMixedEventTripletMultArray23;
  TH2F **fMixedEventTripletMultArray31;

  TH2F **fMixedEventTripletmTArray12; //Mixed Event Q3 vs mT Pair
  TH2F **fMixedEventTripletmTArray23; 
  TH2F **fMixedEventTripletmTArray31; 

  TH2F **fMixedEventTripletPhiThetaArray_SamePair;
  TH2F **fMixedEventTripletPhiThetaArray_DifferentPair;

  TH1F **fMixedEventPairArray_TwoBody; //k* Distribution for (xy)
  TH2F **fMixedEventPairMultArray_TwoBody; //k* Distribution vs mult bin for (xy)
  TH2F **fMixedEventMultArray_TwoBody;
  TH2F **fMixedEventDphiArray_TwoBody;
  TH2F **fMixedEventPairPhiThetaArray_TwoBody;

  // Inv Mass  
  TH2F **fInvMass12; //invariant Mass of pair (ij) vs Q3 of the triplet for (xxy) and (xx)(y) events. In case of  (xx)(y): Array(ij), i,j = 1,2 same event, 3 mixed event
  TH2F **fInvMass23;
  TH2F **fInvMass31;

  TH2F **fInvMassVSmT12; //invariant Mass of pair (ij) vs mT of pair (ij) for (xxy) and (xx)(y) events. In case of  (xx)(y): Array(ij), i,j = 1,2 same event, 3 mixed event
  TH2F **fInvMassVSmT23;
  TH2F **fInvMassVSmT31;

  TH2F **fInvMassVSmT12_MixedEvent; //invariant Mass of pair (ij) vs mT of pair (ij) for (x)(x)(y) events
  TH2F **fInvMassVSmT23_MixedEvent;
  TH2F **fInvMassVSmT31_MixedEvent;

  TH2F **fProjectorData; 

  TList *fSameEventTripletResolutionList; 
  TH2F **fSameEventTripletResolution;
  TH2F **fSameEventTripletResolutionAll;

  TList *fMixedEventTripletResolutionList; 
  TH2F **fMixedEventTripletResolution;
  TH2F **fMixedEventTripletResolutionAll;

  TList *fResultsQA;//!
  AliFemtoDreamControlSample *fSample;   //!
  TList *fResultsSample;//!
  TList *fResultsSampleQA;//!
  int fTrackBufferSize;//
  AliAODTrack **fGTI;  //!
  ClassDef(AliAnalysisTaskAODThreeBodyProtonPrimary,3)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKAODTHREEBODYPROTONPRIMARY_H_ */


