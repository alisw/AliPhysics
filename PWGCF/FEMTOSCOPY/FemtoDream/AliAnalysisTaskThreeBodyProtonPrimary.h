/*
 * AliAnalysisTaskThreeBodyProtonPrimary.h
 *
 *  Created on: May 13, 2019
 *      Authors: Raffaele Del Grande, Marcel Lesch
 *      Based on AliAnalysisTaskThreeBodyFemtoAOD.h from Laura Serksnyte 
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYPROTONPRIMARY_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYPROTONPRIMARY_H_
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

class AliAnalysisTaskThreeBodyProtonPrimary : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskThreeBodyProtonPrimary();
  AliAnalysisTaskThreeBodyProtonPrimary(const char* name, bool isMC);
  virtual ~AliAnalysisTaskThreeBodyProtonPrimary();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  
  void FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* hist2d12, TH2F* hist2d23, TH2F* hist2d31, TH1F* hPtPrimaries, TH1F* hPtProtons, TH1F* hPtPrimaries2, TH1F* hPtProtons2, TH2F* hPtvsQ3Primaries, TH2F* hPtvsQ3Protons, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* Q3VskDistribution12, TH2F* Q3VskDistribution23, TH2F* hKinematics, TH2F* hPrimAngles, TH2F* hDeta, TH2F* hDphi, TH2F* InvMass12, TH2F* InvMass23, TH2F* InvMass31, TH2F* histDeltaPhi12, TH2F* histDeltaPhi23, TH2F* histDeltaPhi31, TH2F* histDeltaPhiDeltaEta12, TH2F* histDeltaPhiDeltaEta23, TH2F* histDeltaPhiDeltaEta31, TH2F* histP, TH2F* histK, TH2F* histmTQ312, TH2F* histmTQ323, TH2F* histmTQ331, TH2F* InvMassVsmT12, TH2F* InvMassVsmT23, TH2F* InvMassVsmT31);

  void SetMixedEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer);

  void FillTripletDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* hist2d12, TH2F* hist2d23, TH2F* hist2d31, TH1F* hPtPrimaries, TH1F* hPtProtons, TH1F* hPtPrimaries2, TH1F* hPtProtons2, TH2F* hPtvsQ3Primaries, TH2F* hPtvsQ3Protons, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* Q3VskDistribution12Mixed, TH2F*  Q3VskDistribution23Mixed, TH2F* hKinematics, TH2F* hPrimAngles, TH2F* hDeta, TH2F* hDphi, TH2F* histDeltaPhi12, TH2F* histDeltaPhi23, TH2F* histDeltaPhi31, TH2F* histDeltaPhiDeltaEta12, TH2F* histDeltaPhiDeltaEta23, TH2F* histDeltaPhiDeltaEta31, TH2F* histP, TH2F* histK, TH2F* histmTQ312, TH2F* histmTQ323, TH2F* histmTQ331, TH2F* InvMassVsmT12, TH2F* InvMassVsmT23, TH2F* InvMassVsmT31); //, TH2F* InvMassMixed, TH2F* Q3VskDistribution12Mixed, TH2F*  Q3VskDistribution23Mixed);

  // test different mixing
  void SetMixedEventPPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void SetMixedEventPPPrim(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void SetMixedEventPPAPrim(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void SetMixedEventPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void SetMixedEventPPrim(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void SetMixedEventPAPrim(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);

  void FillTripletDistributionMETEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);

  // Create triplets like (pp)l (lp)p
  void FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH1F* hPtPrimaries, TH1F* hPtProtons, TH1F* hPtPrimaries2, TH1F* hPtProtons2, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config,TH2F* Q3VskDistribution12, TH2F* Q3VskDistribution23, TH2F* hKinematics, TH2F* hPrimAngles, TH2F* hDeta, TH2F* hDphi, TH2F* hDetaSEvsME, TH2F* hDphiSEvsME, TH2F* InvMass12, TH2F* InvMass23, TH2F* InvMass31, TH2F* histmTQ312,  TH2F* histmTQ323,  TH2F* histmTQ331, TH2F* InvMassVsmT12, TH2F* InvMassVsmT23, TH2F* InvMassVsmT31);

  // Add the close pair cut
  bool DeltaEtaDeltaPhi(int species1, int species2, AliFemtoDreamBasePart &part1,AliFemtoDreamBasePart &part2,  int part1PDGcode, int part2PDGcode, bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist, AliFemtoDreamCollConfig Config, double Q3);
  bool DeltaEtaDeltaPhi(int species1, int species2, AliFemtoDreamBasePart &part1,AliFemtoDreamBasePart &part2,  int part1PDGcode, int part2PDGcode, bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist, AliFemtoDreamCollConfig Config);
  void FillPairInvMass(AliFemtoDreamBasePart &part1, AliFemtoDreamBasePart &part2, AliFemtoDreamBasePart &part3, TH2F* hist, float Q3) ;
  void FillPDGPairInvMass(AliFemtoDreamBasePart &part1, float massPart1, AliFemtoDreamBasePart &part2, float massPart2, AliFemtoDreamBasePart &part3, float massPart3, TH2F* hist, float Q3);

  //Two Body--------------------------------------------------------------------
  void FillPairDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* histDeltaPhi, TH2F* histLowDeltaPhi, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);
  void FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* histDeltaPhi, TH2F* histLowDeltaPhi, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);
  void FillPairTransverseMass(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies, int secondSpecies, TH1F* hist1, std::vector<int> PDGCodes, TH2F* hist2); 


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
  void Setv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {
    fLambda = v0Cuts; 
  }
  void SetAntiv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {
    fAntiLambda = v0Cuts;
  }
  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig=config;
  }  
  void SetRunThreeBodyHistograms(bool RunThreeBodyHistos) {
    fRunThreeBody=RunThreeBodyHistos;
  }

  void SetRunPlotInvMass(bool RunPlotInvMass) {
    fRunPlotInvMass=RunPlotInvMass;
  }

  void SetRunPlotQ3Vsq(bool RunPlotQ3Vsq) {
    fRunPlotQ3Vsq=RunPlotQ3Vsq;
  }

  void SetRunPlotPhiTheta(bool RunPlotPhiTheta) {
    fRunPlotPhiTheta=RunPlotPhiTheta;
  }

  void SetRunPlotOtherHistos(bool RunPlotOtherHistos) {
    fRunPlotOtherHistos=RunPlotOtherHistos;
  }
  void SetRunPlotMult(bool RunPlotMult) {
    fRunPlotMult=RunPlotMult;
  }
  void SetRunCorrDeltaPhi(bool RunCorrDeltaPhi) {
    fRunCorrDeltaPhi=RunCorrDeltaPhi;
  }

  void SetClosePairRejectionForAll(bool ClosePairRejectionForAll) {
    fClosePairRejectionForAll=ClosePairRejectionForAll;
  }
  void SetturnoffClosePairRejectionCompletely(bool turnoffClosePairRejectionCompletely) {
    fturnoffClosePairRejectionCompletely=turnoffClosePairRejectionCompletely;
  }
  void SetClosePairRejectionPPPorPPL(bool ClosePairRejectionPPPorPPL) {
    fClosePairRejectionPPPorPPL=ClosePairRejectionPPPorPPL;
  }

  void SetQ3LimitForDeltaPhiDeltaEta(double Q3LimitForDeltaPhiDeltaEta) {
    fQ3LimitForDeltaPhiDeltaEta=Q3LimitForDeltaPhiDeltaEta;
  }
  
  void SetCleanWithLambdas(bool CleanWithLambdas) {
    fCleanWithLambdas=CleanWithLambdas;
  }
  
  void SetDoOnlyThreeBody(bool DoOnlyThreeBody) {
    fDoOnlyThreeBody=DoOnlyThreeBody;
  }
  void SetRunOfficialTwoBody(bool RunOfficialTwoBody){ // ADDED BY RAFFA
    fRunOfficialTwoBody=RunOfficialTwoBody;
  }

  void SetRunPlotPt(bool RunPlotPt){
    fRunPlotPt=RunPlotPt;
  }
  void SetCutElectrons(bool CutElectrons){
    fCutElectrons=CutElectrons;
  }


  void SetStandardMixing(bool StandardMixing){
    fStandardMixing = StandardMixing;
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

  void SetDoKinematicsPlots(double DoKinematicsPlots){
    fDoKinematicsPlots = DoKinematicsPlots;
  }   

  void SetPlotsMC(bool PlotsMC){
    fPlotsMC = PlotsMC;
  }   

  void SetPlotP1(bool PlotP1){
    fPlotP1 = PlotP1;
  }

  void SetQ3cutValue(double Q3cutValue){
    fQ3cutValue = Q3cutValue;
  }

  void SetQ3MinValue(double Q3MinValue){
    fQ3MinValue = Q3MinValue;
  }

  void SetRunmTPlots(bool RunmTPlots){
    fRunmTPlots = RunmTPlots;
  }

  void SetRunPlotInvMassVSmT(bool RunPlotInvMassVSmT){
    fRunPlotInvMassVSmT = RunPlotInvMassVSmT;
  }



  static TLorentzVector RelativePairMomentum(TLorentzVector &PartOne, TLorentzVector &PartTwo);
  static float BoostOneParticle(TLorentzVector &PartOne, TLorentzVector &PartTwo, TLorentzVector &PartThree);
 private:
  AliAnalysisTaskThreeBodyProtonPrimary(const AliAnalysisTaskThreeBodyProtonPrimary &task);
  AliAnalysisTaskThreeBodyProtonPrimary &operator=(const AliAnalysisTaskThreeBodyProtonPrimary &task);
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
  TList *fSameEventMultDphi;//!
  TList *fMixedEventMultDphi;//!
  TList *fSameEventmT;//!
  TList *fMixedEventmT;//!
  TList *fSameEventPtPrimaries;//!
  TList *fMixedEventPtPrimaries;//!
  TList *fSameEventPtPrimaries2;//!
  TList *fMixedEventPtPrimaries2;//!
  TList *fSameEventPtProtons;//!
  TList *fMixedEventPtProtons;//!
  TList *fSameEventPtProtons2;//!
  TList *fMixedEventPtProtons2;//!
  TList *fSameEventPhiTheta_SamePair;//!
  TList *fMixedEventPhiTheta_SamePair;//!
  TList *fSameEventPhiTheta_DifferentPair;//!
  TList *fMixedEventPhiTheta_DifferentPair;//!
  TList *fQ3Vskq12;//!
  TList *fQ3Vskq12Mixed;//!
  TList *fQ3Vskq23;//!
  TList *fQ3Vskq23Mixed;//!
  TList *fOtherHistos;//!
  TList *fInvMassList;//!
  TList *fKinematicsPlots;//!
  TList *fP1Histos;//!
  TList *fKHistos;//!
  TList *fInvMassVSmTList;//!


  bool fRunThreeBody;
  bool fRunPlotInvMass;
  bool fRunPlotQ3Vsq;
  bool fRunPlotPhiTheta;
  bool fRunPlotOtherHistos;
  bool fRunPlotMult;
  bool fRunCorrDeltaPhi;
  bool fRunPlotPt;
  bool fRunOfficialTwoBody; // ADDED BY RAFFA
  bool fDoKinematicsPlots;
  bool fPlotP1;
  bool fCutElectrons;
  bool fRunPlotInvMassVSmT;


  bool fClosePairRejectionForAll; 
  bool fturnoffClosePairRejectionCompletely; 
  bool fClosePairRejectionPPPorPPL; 

  double fQ3LimitForDeltaPhiDeltaEta;
  double fDeltaPhiMaxPP;
  double fDeltaEtaMaxPP;
  double fDeltaPhiMaxPPrim;
  double fDeltaEtaMaxPPrim;
  double fDeltaPhiMaxPAPrim;
  double fDeltaEtaMaxPAPrim;
  double fQ3cutValue;
  double fQ3MinValue;

  bool fCleanWithLambdas; //if kTRUE: reject Proton + Pi- / Antiproton + Pi+ from Lambda/AntiLambda
  bool fDoOnlyThreeBody; //if kTRUE: 3 Body analysis, else 2 Body
  bool fStandardMixing;
  bool fPlotsMC;
  bool fRunmTPlots;

  TH1F **fSameEventTripletArray;
  TH1F **fSameEventTripletPtPrimaries;
  TH1F **fSameEventTripletPtProtons;
  TH1F **fSameEventTripletPtPrimaries2;
  TH1F **fSameEventTripletPtProtons2;
  TH2F **fSameEventTripletMultArray;
  TH2F **fSameEventTripletMultArray12;
  TH2F **fSameEventTripletMultArray23;
  TH2F **fSameEventTripletMultArray31;
  TH2F **fSameEventTripletmTArray12; 
  TH2F **fSameEventTripletmTArray23; 
  TH2F **fSameEventTripletmTArray31; 
  TH2F **f2Same1MixedEventTripletmTArray12; 
  TH2F **f2Same1MixedEventTripletmTArray23; 
  TH2F **f2Same1MixedEventTripletmTArray31; 
  TH2F **fSameEventTripletPhiThetaArray_SamePair;
  TH2F **fSameEventTripletPhiThetaArray_DifferentPair;
  TH2F **fSameEventTripletPtvsQ3Primaries;
  TH2F **fSameEventTripletPtvsQ3Protons;

  TH1F **fSameEventTripletArray_TwoBody;
  TH2F **fSameEventMultArray_TwoBody;
  TH2F **fSameEventDphiArray_TwoBody;
  TH2F **fSameEventTripletMultArray_TwoBody;
  TH2F **fSameEventTripletPhiThetaArray_TwoBody;
  TH2F **fSameEventP1;
  TH2F **fSameEventK;


  TH1F **fPairTranverseMass_TwoBody; // ADDED BY RAFFA
  TH2F **fPairTranverseMassVSkstar_TwoBody; // ADDED BY RAFFA
  TH2F **fKinematics;
  TH2F **fPrimAngles;
  TH2F **fKinematicsME;
  TH2F **fPrimAnglesME;
  TH2F **fDeta;
  TH2F **fDphi;
  TH2F **fDetaSEvsME;
  TH2F **fDphiSEvsME;
  TH2F **fDetaME;
  TH2F **fDphiME;
  
  // Three particles mixed events
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainer;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerPPP;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerPPPrim;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerPPAPrim;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerPP;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerPPrim;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerPAPrim;
  std::vector<std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>>> fVectPartContainers;


  TH1F **fMixedEventTripletArray;
  TH1F **fMixedEventTripletPtPrimaries;
  TH1F **fMixedEventTripletPtProtons;
  TH1F **fMixedEventTripletPtPrimaries2;
  TH1F **fMixedEventTripletPtProtons2;
  TH2F **fMixedEventTripletMultArray;
  TH2F **fMixedEventTripletMultArray12;
  TH2F **fMixedEventTripletMultArray23;
  TH2F **fMixedEventTripletMultArray31;
  TH2F **fMixedEventTripletmTArray12; 
  TH2F **fMixedEventTripletmTArray23; 
  TH2F **fMixedEventTripletmTArray31; 
  TH2F **fMixedEventTripletPhiThetaArray_SamePair;
  TH2F **fMixedEventTripletPhiThetaArray_DifferentPair;
  TH2F **fMixedEventTripletPtvsQ3Primaries;
  TH2F **fMixedEventTripletPtvsQ3Protons;
  TH2F **fhistDeltaPhi12;
  TH2F **fhistDeltaPhi23;
  TH2F **fhistDeltaPhi31;
  TH2F **fhistDeltaPhiDeltaEta12;
  TH2F **fhistDeltaPhiDeltaEta23;
  TH2F **fhistDeltaPhiDeltaEta31;
  TH2F **fhistDeltaPhi12ME;
  TH2F **fhistDeltaPhi23ME;
  TH2F **fhistDeltaPhi31ME;
  TH2F **fhistDeltaPhiDeltaEta12ME;
  TH2F **fhistDeltaPhiDeltaEta23ME;
  TH2F **fhistDeltaPhiDeltaEta31ME;

  TH1F **fMixedEventTripletArray_TwoBody;
  TH2F **fMixedEventTripletMultArray_TwoBody;
  TH2F **fMixedEventMultArray_TwoBody;
  TH2F **fMixedEventDphiArray_TwoBody;
  TH2F **fMixedEventTripletPhiThetaArray_TwoBody;
  TH2F **fMixedEventP1;
  TH2F **fMixedEventK;

  // Q3 and k* dependense
  TH2F **fQ3VskDistributionsArrayq12;
  TH2F **fQ3VskDistributionsArrayq12Mixed;
  TH2F **fQ3VskDistributionsArrayq23;
  TH2F **fQ3VskDistributionsArrayq23Mixed;
  // doublet vs triplet !!! only for PPP+APAPAP
  TH1F *fDoubletVsTrippletPPP;
  // Inv Mass fits 
  TH2F **fInvMass12;
  TH2F **fInvMass23;
  TH2F **fInvMass31;
  TH2F **fInvMassVSmT12;
  TH2F **fInvMassVSmT23;
  TH2F **fInvMassVSmT31;
  TH2F **fInvMassVSmT12MixedEvent;
  TH2F **fInvMassVSmT23MixedEvent;
  TH2F **fInvMassVSmT31MixedEvent;


  TH2F *fpTvsEtaTrueKaons;
  TH2F *fpTvsEtaTrueAntiKaons;
  TH2F *fpTvsEtaTrueProtons;
  TH2F *fpTvsEtaTrueAntiProtons;
  TH2F *fpTvsEtaRecoKaons;
  TH2F *fpTvsEtaRecoAntiKaons;
  TH2F *fpTvsEtaRecoProtons;
  TH2F *fpTvsEtaRecoAntiProtons;

  TH2F *fTripletInvMassDet;//!
  TH2F *fTripletInvMassPDG;//!
  TH2F *fTripletInvMassDetMixed;//!
  TH2F *fTripletInvMassPDGMixed;//!

  TH2F *fTripletInvMassDetAnti;//!
  TH2F *fTripletInvMassPDGAnti;//!
  TH2F *fTripletInvMassDetMixedAnti;//!
  TH2F *fTripletInvMassPDGMixedAnti;//!

  TList *fResultsQA;//!
  AliFemtoDreamControlSample *fSample;   //!
  TList *fResultsSample;//!
  TList *fResultsSampleQA;//!
  int fTrackBufferSize;//
  AliVTrack **fGTI;  //!
  ClassDef(AliAnalysisTaskThreeBodyProtonPrimary,8)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYPROTONPRIMARY_H_ */


