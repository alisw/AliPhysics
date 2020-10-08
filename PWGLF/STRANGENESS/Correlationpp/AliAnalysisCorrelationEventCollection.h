#ifndef AliAnalysisCorrelationEventCollection_cxx
#define AliAnalysisCorrelationEventCollection_cxx
#include <iostream>
#include <string>
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"

// author: chiara.de.martin@cern.ch from
// ramona.lea@cern.ch from 
// maria.nicassio@cern.ch (derived and adapted from D. Gangadharan PWGCF/FEMTOSCOPY/Chaoticity/AliChaoticityEventCollection
//                                             and J. Salzwedel PWGCF/FEMTOSCOPY/V0LamAnalysis/AliAnalysisV0LamEventCollection)

using namespace std;

class AliReconstructedFirstC {
  
 public:
  AliReconstructedFirstC();
  virtual ~AliReconstructedFirstC();
  /* AliReconstructedFirst(const AliReconstructedFirst &obj); */
  /* AliReconstructedFirst & operator=(const AliReconstructedFirst &obj); */
  
  enum MCFirstOrigin_t {kUnassigned, kFake, kFakeP, kPrimaryP, kPrimaryL, kOtherOriginP, kPrimaryAntiP, kPrimaryAntiL, kOtherOriginAntiP};
  Double_t fMomentum[3]; // 3 reconstructed momentum
  Double_t fMomentumTruth[3]; // 3 true momentum, used in momentum smearing analysis
  Double_t fPt;
  Double_t fEta;
  Double_t fTheta;
  Double_t fPhi;
  Double_t fRap;
  Int_t    fCharge;
  Double_t fDCAxy;
  Double_t fDCAz;
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
  
  Int_t   isP;
  Double_t fMultiplicity;
  Double_t  fZvertex;  
  ClassDef(AliReconstructedFirstC, 1);   

};

class AliReconstructedSecondC {

  
 public:
  AliReconstructedSecondC();
  virtual ~AliReconstructedSecondC();
  /* AliReconstructedSecond(const AliReconstructedSecond &obj); */
  /* AliReconstructedSecond & operator=(const AliReconstructedSecond &obj); */
  
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
  Int_t    isP;

 Double_t  sDcaPosV0;
 Double_t  sDcaNegV0;
 Double_t  sPtArmV0;
 Double_t  sAlphaV0;
 Double_t  sInvMassK0s;
 Double_t  sInvMassLambda;
 Double_t  sInvMassAntiLambda;
 Double_t  sCosPointingAngle;
 Double_t  sDcaV0ToPV;
 Double_t  sMultiplicity;
 Double_t  sZvertex;  
 Double_t  sctau;
 Int_t  sLabelMotherPos;
 Int_t  sLabelPos;
 Int_t  sLabelNeg;
 Bool_t     sAssocOrNot;
 Bool_t     sIsCommonParton;
 Int_t      sPdgCommonParton;
  
 Int_t      cLabelMotherBach;
 Int_t      cisPrimCasc;
 Double_t   cInvMassLambda;
 Double_t   cInvMassXi;
 Double_t   cInvMassOmega;
 Double_t   cCosPointingAngleXi;
 Double_t   cCosPointingAngleV0ToXi;
 Double_t   cDCAXiDaughters;
 Double_t   cRapCasc;
 Double_t   cPt;
 Double_t   cctau;
 Double_t   cEta;
 Double_t   cTheta;
 Double_t   cPhi;
 Int_t      cCharge;
 Bool_t     cAssocOrNot;
 Bool_t     cIsCommonParton;
 Int_t     cPdgCommonParton;
 
  ClassDef(AliReconstructedSecondC, 1);   
  
};

class AliAnalysisCorrelationEvent {

 public:
  AliAnalysisCorrelationEvent();
  virtual ~AliAnalysisCorrelationEvent();
  /* AliAnalysisCorrelationEvent (const AliAnalysisCorrelationEvent &obj); */
  /* AliAnalysisCorrelationEvent &operator=(const AliAnalysisCorrelationEvent &obj); */

  Int_t  fNumberCandidateFirst;
  Int_t  fNumberCandidateSecond; 
  Double_t fPrimaryVertex[3]; //Location of the primary vertex
  
  AliReconstructedFirstC *fReconstructedFirst;
  AliReconstructedSecondC *fReconstructedSecond;
  
  ClassDef(AliAnalysisCorrelationEvent, 1);
  
};



class AliAnalysisCorrelationEventCollection  {

 public:
  
  AliAnalysisCorrelationEventCollection();
  AliAnalysisCorrelationEventCollection(Short_t eventBuffSize, Int_t  maxFirstMult, Int_t  maxSecondMult);
  /* AliAnalysisCorrelationEventCollection(const AliAnalysisCorrelationEventCollection &obj); */
  /* AliAnalysisCorrelationEventCollection & operator=(const AliAnalysisCorrelationEventCollection &obj); */
  virtual  ~AliAnalysisCorrelationEventCollection();

  void FifoShift();
  void FifoClear();
  
  AliAnalysisCorrelationEvent *fEvt;
  
  // private:
  
  Short_t fifo; //Size of the Event Storage buffer
  void SetBuffSize(Short_t eventBuffSize){fifo = eventBuffSize;}
  
  ClassDef(AliAnalysisCorrelationEventCollection, 1);
  
};




#endif
