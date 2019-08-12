#ifndef AliAnalysisHDEventCollection_cxx
#define AliAnalysisHDEventCollection_cxx
#include <iostream>

#include <string>
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"

// author: ramona.lea@cern.ch from 
// maria.nicassio@cern.ch (derived and adapted from D. Gangadharan PWGCF/FEMTOSCOPY/Chaoticity/AliChaoticityEventCollection
//                                             and J. Salzwedel PWGCF/FEMTOSCOPY/V0LamAnalysis/AliAnalysisV0LamEventCollection)

using namespace std;

class AliReconstructed2pcFirst {
  
 public:
  AliReconstructed2pcFirst();
  virtual ~AliReconstructed2pcFirst();
  /* AliReconstructed2pcFirst(const AliReconstructed2pcFirst &obj); */
  /* AliReconstructed2pcFirst & operator=(const AliReconstructed2pcFirst &obj); */
  
  enum MCFirstOrigin_t {kUnassigned, kFake, kFakeP, kPrimaryP, kPrimaryL, kOtherOriginP, kPrimaryAntiP, kPrimaryAntiL, kOtherOriginAntiP};
  Double_t fMomentum[3]; // 3 reconstructed momentum
  Double_t fMomentumTruth[3]; // 3 true momentum, used in momentum smearing analysis
  Double_t fPt;
  Double_t fEta;
  Double_t fTheta;
  Double_t fPhi;
  Double_t fRap;
  Short_t  fCharge;
  Double_t fDCAxy;
  Double_t fDCAz;
  Double_t nSigmaFirstTPC[5];
  Double_t nSigmaFirstTOF[5];
  Bool_t   isTOFmismatch;
  Bool_t   isMCptc;
  Int_t     fMCcode;
  Int_t     fPDGcode;
  Int_t     fMCmumIdx;
  Int_t     fMCmumPDG;
  Int_t     fMCgrandmumIdx;
  Int_t     fMCgrandmumPDG;
  Int_t     index;
  MCFirstOrigin_t mcFirstOriginType;
  Bool_t   doSkipOver;
  Double_t iptoPV[2]; 
  
  Double_t fShiftedGlobalPosition[3];
  Double_t fEtaS;
  Double_t fPhiS;
  
  Bool_t   isP;
  Bool_t   isaP;
  
  ClassDef(AliReconstructed2pcFirst, 1);   

};

class AliReconstructed2pcSecond {

  
 public:
  AliReconstructed2pcSecond();
  virtual ~AliReconstructed2pcSecond();
  /* AliReconstructed2pcSecond(const AliReconstructed2pcSecond &obj); */
  /* AliReconstructed2pcSecond & operator=(const AliReconstructed2pcSecond &obj); */
  
  enum MCSecondOrigin_t {kUnassigned, kFake, kFakeP, kPrimaryP, kPrimaryL, kOtherOriginP, kPrimaryAntiP, kPrimaryAntiL, kOtherOriginAntiP};
  Double_t sMomentum[3]; // 3 reconstructed momentum
  Double_t sMomentumTruth[3]; // 3 true momentum, used in momentum smearing analysis
  Double_t sPt;
  Double_t sEta;
  Double_t sTheta;
  Double_t sPhi;
  Double_t sRap;
  Short_t  sCharge;
  Double_t sDCAxy;
  Double_t sDCAz;
  Double_t nSigmaSecondTPC[5];
  Double_t nSigmaSecondTOF[5];
  Bool_t   isTOFmismatch;
  Bool_t   isMCptc;
  Int_t     sMCcode;
  Int_t     sPDGcode;
  Int_t     sMCmumIdx;
  Int_t     sMCmumPDG;
  Int_t     sMCgrandmumIdx;
  Int_t     sMCgrandmumPDG;
  Int_t     index;
  MCSecondOrigin_t mcSecondOriginType;
  Bool_t   doSkipOver;
  Double_t iptoPV[2]; 
  
  Double_t sShiftedGlobalPosition[3];
  Double_t sEtaS;
  Double_t sPhiS;
  Double_t sTOFMass;
  Bool_t   isP;
  Bool_t   isaP;
  
  ClassDef(AliReconstructed2pcSecond, 1);   
  
};

class AliAnalysishDEvent {

 public:
  AliAnalysishDEvent();
  virtual ~AliAnalysishDEvent();
  AliAnalysishDEvent (const AliAnalysishDEvent &obj);
  AliAnalysishDEvent &operator=(const AliAnalysishDEvent &obj);

  Int_t  fNumberCandidateFirst;
  Int_t  fNumberCandidateSecond; 
  Double_t fPrimaryVertex[3]; //Location of the primary vertex
  
  AliReconstructed2pcFirst *fReconstructedFirst;
  AliReconstructed2pcSecond *fReconstructedSecond;
  
  ClassDef(AliAnalysishDEvent, 1);
  
};



class AliAnalysishDEventCollection  {

 public:
  
  AliAnalysishDEventCollection();
  AliAnalysishDEventCollection(Short_t eventBuffSize, Int_t  maxFirstMult, Int_t  maxSecondMult);
  AliAnalysishDEventCollection(const AliAnalysishDEventCollection &obj);
  AliAnalysishDEventCollection & operator=(const AliAnalysishDEventCollection &obj);
  virtual  ~AliAnalysishDEventCollection();

  void FifoShift();
  
  AliAnalysishDEvent *fEvt;
  
  // private:
  
  Short_t fifo; //Size of the Event Storage buffer
  void SetBuffSize(Short_t eventBuffSize){fifo = eventBuffSize;}
  
  ClassDef(AliAnalysishDEventCollection, 1);
  
};




#endif
