#ifndef AliAnalysishCascadeEventCollection_cxx
#define AliAnalysishCascadeEventCollection_cxx
#include <iostream>

#include <string>

#include "TH1.h"

#include "TH2.h"

// author: maria.nicassio@cern.ch (derived and adapted from D. Gangadharan PWGCF/FEMTOSCOPY/Chaoticity/AliChaoticityEventCollection
//                                 and J. Salzwedel PWGCF/FEMTOSCOPY/V0LamAnalysis/AliAnalysisV0LamEventCollection)
using namespace std;

 class AliReconstructedProton {


  public:
   AliReconstructedProton();

   ~AliReconstructedProton();

//   AliReconstructedProton(const AliReconstructedProton&);
//   AliReconstructedProton & operator=(const AliReconstructedProton&);
 
   enum MCProtonOrigin_t {kUnassigned, kPrimaryP, kSecondaryWeak, kMaterial};
   double pMomentum[3]; // 3 reconstructed momentum
   double pMomentumTruth[3]; // 3 true momentum, used in momentum smearing analysis
   double pPt;
   double pEta;
   double pPhi;
   double pRap;
   short pCharge;
   int nSigmaProtonTPC;
   int nSigmaProtonTOF;
   bool isTOFmismatch;
   bool isP;
   bool isaP;
   int  index;
   MCProtonOrigin_t mcProtonOriginType;
   double iptoPV[2]; 
   bool doSkipOver;

   Double_t pShiftedGlobalPosition[3];
   Double_t pEtaS;
   Double_t pPhiS;

   ClassDef(AliReconstructedProton, 1);   

};


 
 class AliReconstructedXi {

  public:
   AliReconstructedXi();
   ~AliReconstructedXi();
//   AliReconstructedXi(const AliReconstructedXi&);
//   AliReconstructedXi & operator=(const AliReconstructedXi&);

   enum MCXiOrigin_t {kUnassigned, kFake, kFakeXi, kPrimaryXi, kPrimaryOmega, kOtherOriginXi, kPrimaryAntiXi, kPrimaryAntiOmega, kOtherOriginAntiXi, kSecondaryXi, kSecondaryOmega,kSecondaryAntiXi, kSecondaryAntiOmega};
   double xiMomentum[3]; // reconstructed momentum
   double xiMomentumTruth[3]; //true momentum, used in momentum smearing analysis
   double xiPt;
   double xiEta;
   double xiPhi;
   double xiRap;
//   double* v0Momentum;//[3]; // reconstructed momentum
//   double* v0MomentumTruth;//[3]; //true momentum, used in momentum smearing analysis
//   double v0Pt;
//   double v0Eta;
//   double v0Phi;  // FIXME maybe needed also properties of bach
   double xiMass;
//   double massLam;
//   double massALam;
//   double massLamDifference;
//   double massALamDifference;
//   double massK0;
//   double lorentzGammaLam;
//   double decayLength;
//   short daughter1ID;
//   short daughter2ID;
//   double v0DCA;
//   double* decayVertexPosition;//[3]; //position of the V0 decay vertex
//   double cosPointing;
//   bool hasProtonDaughter;
//   bool hasAntiProtonDaughter;
//   bool hasPiPlusDaughter;
//   bool hasPiMinusDaughter;
   bool isXiMass; 
   bool isXi;
   bool isaXi;
   int  indexB;
   int  indexP;
   int  indexN;
   bool doSkipOver;
   bool doPickOne;

   MCXiOrigin_t mcXiOriginType;
//   double* daughterPosMomentum;//[3]; //momentum of the positive daughter
//   double* daughterNegMomentum;//[3]; //momentum of the negative daughter
   /* double daughterPosMomentumTruth[3]; //used in momentum smearing analysis */
   /* double doughterNegMomentumTruth[3]; */
//   double* daughterBacMomentum;//[3]; //momentum of the bachelor daughter
//   double daughterPosProtonE;     //energy of the positive daughter
//   double daughterPosPionE;       //energy of the positive daughter
//   double daughterNegProtonE;     //energy of the negative daughter
//   double daughterNegPionE;       //energy of the negative daughter
//   double daughtersDCA;
//   double daughterPosDCAPrimaryVertex;
//   double daughterNegDCAPrimaryVertex;
//   double daughterBacDCAPrimaryVertex;
//   double* daughterPosPositionDCA;//[3]; //position of the positive daughter at DCA to primary
//   double* daughterNegPositionDCA;//[3]; //position of the negative daughter at DCA to primary
//   double* daughterPosMomentumDCA;//[3]; //momentum of the positive daughter at DCA to primary
//   double* daughterNegMomentumDCA;//[3]; //momentum of the negative daughter at DCA to primary

//  double* daughterPosCovariance;//[21]; //covariance matrix of the positive daughter

//  double* daughterNegCovariance;//[21]; //covariance matrix of the negative daughter
  
//  double* daughterBacCovariance;//[21]; //covariance matrix of the negative daughter

//  Float_t* daughterPosGlobalPositions;//[9][3];

//  Float_t* daughterNegGlobalPositions;//[9][3];
  
//  Float_t* daughterBacGlobalPositions;//[9][3];

  Double_t daughterPosShiftedGlobalPosition[3]; 

  Double_t daughterNegShiftedGlobalPosition[3];

  Double_t daughterBacShiftedGlobalPosition[3];

  Double_t daughterBacEtaS;
  Double_t daughterPosEtaS;
  Double_t daughterNegEtaS;
  Double_t daugherBacPhiS;
  Double_t daugherPosPhiS;
  Double_t daugherNegPhiS;

  ClassDef(AliReconstructedXi, 1);

};


class AliAnalysishCascadeEvent {

   public:

    int fNumberCandidateXi;
    int fNumberCandidateProton; 
    double fPrimaryVertex[3]; //Location of the primary vertex

    AliReconstructedXi *fReconstructedXi;
    AliReconstructedProton *fReconstructedProton;

    ClassDef(AliAnalysishCascadeEvent, 1);

};


class AliAnalysishCascadeEventCollection  {

   public:

    AliAnalysishCascadeEventCollection();

    AliAnalysishCascadeEventCollection(short eventBuffSize, int maxXiMult, int maxProtonMult);
    AliAnalysishCascadeEventCollection(const AliAnalysishCascadeEventCollection&);
    AliAnalysishCascadeEventCollection & operator=(const AliAnalysishCascadeEventCollection&);

    ~AliAnalysishCascadeEventCollection();

    void FifoShift();

    AliAnalysishCascadeEvent *fEvt;

   private:

    short fifo; //Size of the Event Storage buffer

    void SetBuffSize(short eventBuffSize){fifo = eventBuffSize;}

   ClassDef(AliAnalysishCascadeEventCollection, 1);

};




#endif
