/*
 * AliAnalysisTaskThreeBodyFemto.h
 *
 *  Created on: May 13, 2019
 *      Author: Laura Serksnyte 
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYFEMTO_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYFEMTO_H_
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

class AliAnalysisTaskThreeBodyFemto : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskThreeBodyFemto();
  AliAnalysisTaskThreeBodyFemto(const char* name, bool isMC);
  virtual ~AliAnalysisTaskThreeBodyFemto();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  
  void FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);//, TH2F* InvMassSame);
  void FillTripletDistributionPPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* InvMassSame, TH2F* InvMassDET,TH2F* InvMassPDG);

  void FillPairDistributionPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, TH1F* sameEventDistributionPL, int mult, TH2F* hist2d);
  void FillPairDistributionPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, TH1F* sameEventDistributionPP, TH2F **fEventTripletPhiThetaArray,  AliFemtoDreamCollConfig Config);
  void SetMixedEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer);
  void FillTripletDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config); //, TH2F* InvMassMixed, TH2F* Q3VskDistribution12Mixed, TH2F*  Q3VskDistribution23Mixed);
  void FillTripletDistributionMEPPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* InvMassMixed, TH2F* Q3VskDistribution12Mixed,TH2F* InvMassDET,TH2F* InvMassPDG);

// test different mixing
  void SetMixedEventOnlyPLambdaTEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void FillTripletDistributionMETEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config);
  void FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray,  AliFemtoDreamCollConfig Config);
  // test different mixing 2
  void SetMixedEventOnlyPPLambdaTEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void SetMixedEventOnlyPPPTEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  // Create triplets like (pp)l (lp)p
  void FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config); //,TH2F* Q3VskDistribution12, TH2F* Q3VskDistribution23);
  // Add the close pair cut
  bool DeltaEtaDeltaPhi(AliFemtoDreamBasePart &part1,AliFemtoDreamBasePart &part2, bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist, AliFemtoDreamCollConfig Config, double Q3);
  bool DeltaEtaDeltaPhi(AliFemtoDreamBasePart &part1,AliFemtoDreamBasePart &part2, bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist, AliFemtoDreamCollConfig Config);
  void FillPairInvMass(AliFemtoDreamBasePart &part1, AliFemtoDreamBasePart &part2, AliFemtoDreamBasePart &part3, TH2F* hist, float Q3) ;
  void FillPDGPairInvMass(AliFemtoDreamBasePart &part1, float massPart1, AliFemtoDreamBasePart &part2, float massPart2, AliFemtoDreamBasePart &part3, float massPart3, TH2F* hist, float Q3);


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
  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig=config;
  }  
  void SetRunThreeBodyHistograms(bool RunThreeBodyHistos) {
    fRunThreeBody=RunThreeBodyHistos;
  }

  void SetRunPlotInvMassTriplet(bool RunPlotInvMassTriplet) {
    fRunPlotInvMassTriplet=RunPlotInvMassTriplet;
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
  

  

  static TLorentzVector RelativePairMomentum(TLorentzVector &PartOne, TLorentzVector &PartTwo);
 private:
  AliAnalysisTaskThreeBodyFemto(const AliAnalysisTaskThreeBodyFemto &task);
  AliAnalysisTaskThreeBodyFemto &operator=(const AliAnalysisTaskThreeBodyFemto &task);
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
  TList *fQ3Vskq12;//!
  TList *fQ3Vskq12Mixed;//!
  TList *fQ3Vskq23;//!
  TList *fQ3Vskq23Mixed;//!
  TList *fOtherHistos;//!
  TList *fInvMassTripletSame;//!
  TList *fInvMassTripletMixed;//!

  bool fRunThreeBody;
  bool fRunPlotInvMassTriplet;
  bool fRunPlotQ3Vsq;
  bool fRunPlotPhiTheta;
  bool fRunPlotOtherHistos;
  bool fRunPlotMult;

  bool fClosePairRejectionForAll;
  bool fturnoffClosePairRejectionCompletely;
  bool fClosePairRejectionPPPorPPL;

  double fQ3LimitForDeltaPhiDeltaEta;

  TH1F **fSameEventTripletArray;
  TH2F **fSameEventTripletMultArray;
  TH2F **fSameEventTripletPhiThetaArray;
  // Three particles mixed events
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainer;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerTEST;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerTESTppL;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerTESTppp;

  TH1F **fMixedEventTripletArray;
  TH2F **fMixedEventTripletMultArray;
  TH2F **fMixedEventTripletPhiThetaArray;
  // Q3 and k* dependense
  TH2F **fQ3VskDistributionsArrayq12;
  TH2F **fQ3VskDistributionsArrayq12Mixed;
  TH2F **fQ3VskDistributionsArrayq23;
  TH2F **fQ3VskDistributionsArrayq23Mixed;
  // doublet vs triplet !!! only for PPP+APAPAP
  TH1F *fDoubletVsTrippletPPP;
  // Inv Mass fits 
  TH2F **fInvMassSame;
  TH2F **fInvMassMixed;

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
  ClassDef(AliAnalysisTaskThreeBodyFemto,3)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYFEMTO_H_ */


