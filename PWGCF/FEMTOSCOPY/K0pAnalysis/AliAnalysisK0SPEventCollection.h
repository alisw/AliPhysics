#ifndef AliAnalysisK0SPEventCollection_cxx
#define AliAnalysisK0SPEventCollection_cxx
#include <iostream>
#include <string>
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"

// author: Marta Urioni  from
// ramona.lea@cern.ch from 
// maria.nicassio@cern.ch (derived and adapted from D. Gangadharan PWGCF/FEMTOSCOPY/Chaoticity/AliChaoticityEventCollection
//                                             and J. Salzwedel PWGCF/FEMTOSCOPY/V0LamAnalysis/AliAnalysisV0LamEventCollection)

using namespace std;

class AliReconstructedFirstK0SP {
  
 public:
  AliReconstructedFirstK0SP();
  virtual ~AliReconstructedFirstK0SP();
  /* AliReconstructedFirstK0SP(const AliReconstructedFirstK0SP &obj); */
  /* AliReconstructedFirstK0SP & operator=(const AliReconstructedFirstK0SP &obj); */
 
  enum MCFirstOrigin_t {kUnassigned, kFake, kFakeP, kPrimaryP, kPrimaryL, kOtherOriginP, kPrimaryAntiP, kPrimaryAntiL, kOtherOriginAntiP};
  Double_t fMomentum[3]; // 3 reconstructed momentum
  Double_t fMomentumTruth[3]; // 3 true momentum, used in momentum smearing analysis
  Double_t fPt;  ///
  Double_t fTheta; ///
  Double_t fPhi; ///
  Double_t fDcaPosV0; ///
  Double_t fDcaNegV0; ///
  Double_t fInvMassK0s; ///
  Double_t fInvMassLambda; ///
  Double_t fInvMassAntiLambda; ///
  Double_t fCosPointingAngle; ///
  Int_t    fPDGcode; ///

  Double_t fEta;
  Double_t fRap;
  Short_t  fCharge;
  Double_t fDCAxy;
  Double_t fDCAz;
  Int_t fMCcode;
  Bool_t  isMCptc;
  Int_t     fMCmumIdx;
  Int_t     fMCmumPDG;
  Int_t     fMCgrandmumIdx;
  Int_t     fMCgrandmumPDG;
  Int_t *fAncestorParticleLabel; //!
  Int_t *fAncestorPdg; //!
  Int_t     index;
  Int_t     indexPosdaughter;
  Int_t     indexNegdaughter;
  MCFirstOrigin_t mcFirstOriginType;
  Bool_t   doSkipOver;
  Double_t iptoPV[2]; 
  
  Double_t fShiftedGlobalPosition[3];
  Double_t fEtaS;
  Double_t fPhiS;
  Bool_t   isP;

  Double_t  fPtArmV0;
  Double_t  fAlphaV0;
  Double_t  fDcaV0ToPV;
  Double_t  fMultiplicity;
  Double_t  fZvertex;  
  Double_t  fctau;
  Int_t  fLabelMotherPos;
  Int_t  fLabelPos;
  Int_t  fLabelNeg;
  Bool_t     fAssocOrNot;
  
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
 
  ClassDef(AliReconstructedFirstK0SP, 1);   

};

class AliReconstructedSecondK0SP {

  
 public:
  AliReconstructedSecondK0SP();
  virtual ~AliReconstructedSecondK0SP();
  /* AliReconstructedSecondK0SP(const AliReconstructedSecondK0SP &obj); */
  /* AliReconstructedSecondK0SP & operator=(const AliReconstructedSecondK0SP &obj); */
  
  enum MCSecondOrigin_t {kUnassigned, kFake, kFakeP, kPrimaryP, kPrimaryL, kOtherOriginP, kPrimaryAntiP, kPrimaryAntiL, kOtherOriginAntiP};
  Double_t sMomentum[3]; // 3 reconstructed momentum
  Double_t sMomentumTruth[3]; // 3 true momentum, used in momentum smearing analysis
  Double_t sPt; ///
  Double_t sTheta; ///
  Double_t sPhi; ///
  Int_t    sCharge; ///
  Double_t sDCAxy; ///
  Double_t sDCAz; ///
  Double_t sMassTOF; /// 
  Int_t    sPDGcode; ///
  
  Double_t sRap;
  Double_t sEta;
  Int_t sMCcode;
  Bool_t  isMCptc;
  Int_t     sMCmumIdx;
  Int_t     sMCmumPDG;
  Int_t     sMCgrandmumIdx;
  Int_t     sMCgrandmumPDG;
  Int_t *sAncestorParticleLabel; //!
  Int_t *sAncestorPdg; //!
  Int_t     index;
  MCSecondOrigin_t mcSecondOriginType;
  Bool_t   doSkipOver;
  Double_t iptoPV[2]; 
  
  Double_t sShiftedGlobalPosition[3];
  Double_t sEtaS;
  Double_t sPhiS;
  
  Int_t   isP;
  Int_t   isaP;
  Double_t sMultiplicity;
  Double_t  sZvertex;  
 

  ClassDef(AliReconstructedSecondK0SP, 1);   
  
};

class AliAnalysisK0SPEvent {

 public:
  AliAnalysisK0SPEvent();
  virtual ~AliAnalysisK0SPEvent();
  /* AliAnalysisK0SPEvent (const AliAnalysisK0SPEvent &obj); */
  /* AliAnalysisK0SPEvent &operator=(const AliAnalysisK0SPEvent &obj); */

  Int_t  fNumberCandidateFirst;
  Int_t  fNumberCandidateSecond; 
  Double_t fPrimaryVertex[3]; //Location of the primary vertex
  
  AliReconstructedFirstK0SP *fReconstructedFirst;
  AliReconstructedSecondK0SP *fReconstructedSecond;
  
  ClassDef(AliAnalysisK0SPEvent, 1);
  
};



class AliAnalysisK0SPEventCollection  {

 public:
  
  AliAnalysisK0SPEventCollection();
  AliAnalysisK0SPEventCollection(Short_t eventBuffSize, Int_t  maxFirstMult, Int_t  maxSecondMult);
  /* AliAnalysisK0SPEventCollection(const AliAnalysisK0SPEventCollection &obj); */
  /* AliAnalysisK0SPEventCollection & operator=(const AliAnalysisK0SPEventCollection &obj); */
  virtual  ~AliAnalysisK0SPEventCollection();

  void FifoShift();
  void FifoClear();
  
  AliAnalysisK0SPEvent *fEvt;
  
  // private:
  
  Short_t fifo; //Size of the Event Storage buffer
  void SetBuffSize(Short_t eventBuffSize){fifo = eventBuffSize;}
  
  ClassDef(AliAnalysisK0SPEventCollection, 1);
  
};




#endif
