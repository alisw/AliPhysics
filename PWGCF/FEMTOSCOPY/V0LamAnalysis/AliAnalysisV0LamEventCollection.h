#ifndef AliAnalysisV0LamEventCollection_cxx
#define AliAnalysisV0LamEventCollection_cxx

#include <iostream>
#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
using namespace std;

class AliReconstructedV0  // Reconstructed V0s
{
 public:
  AliReconstructedV0();
  ~AliReconstructedV0();
  enum MCV0Origin_t {kUnassigned, kFake, kFakeLambda, kPrimaryLambda, kPrimarySigmaZero, kExcitedSigma, kPrimaryCascadeZero, kPrimaryCascadeMinus, kPrimaryOmega, kOtherOriginLambda, kFakeAntiLambda, kPrimaryAntiLambda, kPrimaryAntiSigmaZero, kExcitedAntiSigma, kPrimaryAntiCascadeZero, kPrimaryAntiCascadePlus, kPrimaryAntiOmega, kOtherOriginAntiLambda, kKZeroShort, kOriginTypeMax = kKZeroShort};
  TVector3 v0Momentum; // reconstructed momentum
  TVector3 v0MomentumTruth; //true momentum, used in momentum smearing analysis
  double v0Pt;
  double v0Eta;
  double v0Phi;
  double massLam;
  double massALam;
  double massLamDifference;
  double massALamDifference;
  double massK0;
  double lorentzGammaLam;
  double decayLength;
  short daughter1ID;
  short daughter2ID;
  double v0DCA;
  double cosPointing;
  bool hasProtonDaughter;
  bool hasAntiProtonDaughter;
  bool hasPiPlusDaughter;
  bool hasPiMinusDaughter;
  bool isLamCenter[10];
  bool isALamCenter[10];
  bool isDeemedUnworthy[10];
  //vector<vector<bool> > hasPassedCut;
  bool hasPassedCut[10][10];
  bool isPassingAllReconstructionCuts;
  MCV0Origin_t mcOriginType;
  TVector3 daughterPosMomentum; //momentum of the positive daughter
  TVector3 daughterNegMomentum; //momentum of the negative daughter
  /* double daughterPosMomentumTruth[3]; //used in momentum smearing analysis */
  /* double doughterNegMomentumTruth[3]; */
  double daughterPosProtonE;     //energy of the positive daughter
  double daughterPosPionE;       //energy of the positive daughter
  double daughterNegProtonE;	 //energy of the negative daughter
  double daughterNegPionE;	 //energy of the negative daughter
  double daughtersDCA;
  double daughterPosDCAPrimaryVertex;
  double daughterNegDCAPrimaryVertex;
  double daughterPosPositionDCA[3]; //position of the positive daughter at DCA to primary
  double daughterNegPositionDCA[3]; //position of the negative daughter at DCA to primary
  double daughterPosMomentumDCA[3]; //momentum of the positive daughter at DCA to primary
  double daughterNegMomentumDCA[3]; //momentum of the negative daughter at DCA to primary
  double daughterPosCovariance[21]; //covariance matrix of the positive daughter
  double daughterNegCovariance[21]; //covariance matrix of the negative daughter
  vector<TVector3> daughterPosGlobalPositions;
  vector<TVector3> daughterNegGlobalPositions;
  vector<TVector3> daughterPosCorrectedGlobalPositions;
  vector<TVector3> daughterNegCorrectedGlobalPositions;
  /* TVector3 emissionPoint; */
};

class AliAnalysisV0LamEvent
{
 public:
  AliAnalysisV0LamEvent();
  AliAnalysisV0LamEvent(const AliAnalysisV0LamEvent &event);
  ~AliAnalysisV0LamEvent();
  AliAnalysisV0LamEvent& operator=(const AliAnalysisV0LamEvent &event);
  Int_t fMaxV0Mult;
  Int_t fNumberCandidateV0;
  TVector3 fPrimaryVertex; //Location of the primary vertex
  AliReconstructedV0 *fReconstructedV0;
};

class AliAnalysisV0LamEventCollection
{
 public:
  AliAnalysisV0LamEventCollection();
  AliAnalysisV0LamEventCollection(short eventBuffSize,int maxV0Mult);
  ~AliAnalysisV0LamEventCollection();
  AliAnalysisV0LamEventCollection(const AliAnalysisV0LamEventCollection &eventCollection);
  AliAnalysisV0LamEventCollection& operator=(const AliAnalysisV0LamEventCollection& eventCollection);
  void FifoShift();
  short GetFifoSize() const {return fFifoSize;};
  AliAnalysisV0LamEvent *fEvt;
 private:
  short fFifoSize; //Size of the Event Storage buffer
  void SetBuffSize(short eventBuffSize){fFifoSize = eventBuffSize;}
};
#endif

