#ifndef AliAnalysisV0LamEventCollection_cxx
#define AliAnalysisV0LamEventCollection_cxx

#include <iostream>
#include <string>
#include "TH1.h"
#include "TH2.h"
using namespace std;

class AliReconstructedV0  // Reconstructed V0s
{
 public:
  AliReconstructedV0();
  ~AliReconstructedV0();
  enum MCV0Origin_t {kUnassigned, kFake, kFakeLambda, kPrimaryLambda, kPrimarySigmaZero, kExcitedSigma, kPrimaryCascadeZero, kPrimaryCascadeMinus, kPrimaryOmega, kOtherOriginLambda, kFakeAntiLambda, kPrimaryAntiLambda, kPrimaryAntiSigmaZero, kExcitedAntiSigma, kPrimaryAntiCascadeZero, kPrimaryAntiCascadePlus, kPrimaryAntiOmega, kOtherOriginAntiLambda, kKZeroShort, kOriginTypeMax = kKZeroShort};
  double v0Momentum[3]; // reconstructed momentum
  double v0MomentumTruth[3]; //true momentum, used in momentum smearing analysis
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
  double decayVertexPosition[3]; //position of the V0 decay vertex
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
  double daughterPosMomentum[3]; //momentum of the positive daughter
  double daughterNegMomentum[3]; //momentum of the negative daughter
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
  Float_t daughterPosGlobalPositions[9][3];
  Float_t daughterNegGlobalPositions[9][3];
  Float_t daughterPosCorrectedGlobalPositions[9][3];
  Float_t daughterNegCorrectedGlobalPositions[9][3];
};

class AliAnalysisV0LamEvent
{
 public:
  int fNumberCandidateV0;
  double fPrimaryVertex[3]; //Location of the primary vertex
  AliReconstructedV0 *fReconstructedV0;
};

class AliAnalysisV0LamEventCollection
{
 public:
  AliAnalysisV0LamEventCollection();
  AliAnalysisV0LamEventCollection(short eventBuffSize,int maxV0Mult);
  ~AliAnalysisV0LamEventCollection();
  void FifoShift();
  AliAnalysisV0LamEvent *fEvt;
 private:
  short fifo; //Size of the Event Storage buffer
  void SetBuffSize(short eventBuffSize){fifo = eventBuffSize;}
};
#endif

