#ifndef ALIESDV0PARAMS_H
#define ALIESDV0PARAMS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          ESD V0 Vertex Class - parameterization 
//          This class is part of the Event Summary Data set of classes
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------

#include "TObject.h"

class AliESDV0Params : public TObject{
  friend class AliESDv0;
 public:
  AliESDV0Params();
 
  void SetMaxDist0(Float_t kMaxDist0=0.1) {fkMaxDist0=kMaxDist0;}
  void SetMaxDist1(Float_t kMaxDist1=0.1) {fkMaxDist1=kMaxDist1;}
  void SetMaxDist(Float_t kMaxDist=1.) {fkMaxDist=kMaxDist;}
  void SetMinPointAngle(Float_t kMinPointAngle=0.85){fkMinPointAngle=kMinPointAngle;}
  void SetMinPointAngle2(Float_t kMinPointAngle2=0.99){fkMinPointAngle2=kMinPointAngle2;}
  void SetMinR(Float_t kMinR= 0.5){fkMinR= kMinR;}
  void SetMaxR(Float_t kMaxR= 220.){fkMaxR= kMaxR;}
  void SetMinPABestConst(Float_t kMinPABestConst= 0.9999){fkMinPABestConst= kMinPABestConst;}
  void SetMaxRBestConst(Float_t kMaxRBestConst= 10.){fkMaxRBestConst= kMaxRBestConst;}

  void SetCausality0Cut(Float_t kCausality0Cut= 0.19){fkCausality0Cut= kCausality0Cut;}
  void SetLikelihood01Cut(Float_t kLikelihood01Cut = 0.45) {fkLikelihood01Cut=kLikelihood01Cut;}
  void SetLikelihood1Cut(Float_t kLikelihood1Cut = 0.5) {fkLikelihood1Cut=kLikelihood1Cut;}
  void SetCombinedCut(Float_t kCombinedCut= 0.55) {fkCombinedCut=kCombinedCut;}
  void SetMinClFullTrk(Float_t kMinClFullTrk=5.0){fkMinClFullTrk=kMinClFullTrk;}
  void SetMinTgl0(Float_t kMinTgl0=1.05){fkMinTgl0=kMinTgl0;}
  void SetMinRTgl0(Float_t kMinRTgl0=40.0){fkMinRTgl0=kMinRTgl0;}

  void SetMinNormDistForbTgl0(Float_t kMinNormDistForbTgl0=3.0){fkMinNormDistForbTgl0=kMinNormDistForbTgl0;}
  void SetMinClForb0(Float_t kMinClForb0=4.5){fkMinClForb0=kMinClForb0;}
  void SetMinNormDistForb1(Float_t kMinNormDistForb1=3.0){fkMinNormDistForb1=kMinNormDistForb1;}
  void SetMinNormDistForb2(Float_t kMinNormDistForb2=2.0){fkMinNormDistForb2=kMinNormDistForb2;}
  void SetMinNormDistForb3(Float_t kMinNormDistForb3=1.0){fkMinNormDistForb3=kMinNormDistForb3;}
  void SetMinNormDistForb4(Float_t kMinNormDistForb4=4.0){fkMinNormDistForb4=kMinNormDistForb4;}
  void SetMinNormDistForb5(Float_t kMinNormDistForb5=5.0){fkMinNormDistForb5=kMinNormDistForb5;}

  void SetMinRTPCdensity(Float_t kMinRTPCdensity=40.0){fkMinRTPCdensity=kMinRTPCdensity;}
  void SetMaxRTPCdensity0( Float_t kMaxRTPCdensity0=110.0){fkMaxRTPCdensity0=kMaxRTPCdensity0;}
  void SetMaxRTPCdensity10(Float_t kMaxRTPCdensity10=120.0){fkMaxRTPCdensity10=kMaxRTPCdensity10;}
  void SetMaxRTPCdensity20(Float_t kMaxRTPCdensity20=130.0){fkMaxRTPCdensity20=kMaxRTPCdensity20;}
  void SetMaxRTPCdensity30(Float_t kMaxRTPCdensity30=140.0){fkMaxRTPCdensity30=kMaxRTPCdensity30;}

  void SetMinTPCdensity(Float_t kMinTPCdensity=0.6){fkMinTPCdensity=kMinTPCdensity;}
  void SetMinTgl1(Float_t kMinTgl1=1.1){fkMinTgl1=kMinTgl1;}
  void SetMinTgl2(Float_t kMinTgl2=1.){fkMinTgl1=kMinTgl2;}
  void SetMinchi2before0(Float_t kMinchi2before0=16.){fkMinchi2before0=kMinchi2before0;}
  void SetMinchi2before1(Float_t kMinchi2before1=16.){fkMinchi2before1=kMinchi2before1;}
  void SetMinchi2after0(Float_t kMinchi2after0=16.){fkMinchi2after0=kMinchi2after0;}
  void SetMinchi2after1(Float_t kMinchi2after1=16.){fkMinchi2after1=kMinchi2after1;}
  void SetAddchi2SharedCl(Float_t kAddchi2SharedCl=18){fkAddchi2SharedCl=kAddchi2SharedCl;}
  void SetAddchi2NegCl0(Float_t kAddchi2NegCl0=25){fkAddchi2NegCl0=kAddchi2NegCl0;}
  void SetAddchi2NegCl1(Float_t kAddchi2NegCl1=30){fkAddchi2NegCl1=kAddchi2NegCl1;}
  void SetSigp0Par0(Float_t kSigp0Par0=0.0001){fkSigp0Par0=kSigp0Par0;}
  void SetSigp0Par1(Float_t kSigp0Par1=0.001){fkSigp0Par1=kSigp0Par1;}
  void SetSigp0Par2(Float_t kSigp0Par2=0.1){fkSigp0Par2=kSigp0Par2;}

  void SetSigpPar0(Float_t kSigpPar0=0.5){fkSigpPar0=kSigpPar0;}
  void SetSigpPar1(Float_t kSigpPar1=0.6){fkSigpPar1=kSigpPar1;}
  void SetSigpPar2(Float_t kSigpPar2=0.4){fkSigpPar2=kSigpPar2;}
  void SetMaxDcaLh0(Float_t kMaxDcaLh0=0.5) {fkMaxDcaLh0=kMaxDcaLh0;}

  void SetChi2KF(Float_t kChi2KF=100.){fkChi2KF=kChi2KF;}
  void SetRobustChi2KF(Float_t kRobustChi2KF=100.){fkRobustChi2KF=kRobustChi2KF;}


  Float_t GetMaxDist0() const {return fkMaxDist0;}
  Float_t GetMaxDist1() const {return fkMaxDist1;}
  Float_t GetMaxDist() const {return fkMaxDist;}
  Float_t GetMinPointAngle() const {return fkMinPointAngle;}
  Float_t GetMinPointAngle2() const {return fkMinPointAngle2;}
  Float_t GetMinR() const {return fkMinR;}
  Float_t GetMaxR() const {return fkMaxR;}
  Float_t GetMinPABestConst() const {return fkMinPABestConst;}
  Float_t GetMaxRBestConst() const {return fkMaxRBestConst;}
  Float_t GetCausality0Cut() const {return fkCausality0Cut;}
  Float_t GetLikelihood01Cut() const {return fkLikelihood01Cut;}
  Float_t GetLikelihood1Cut() const {return fkLikelihood1Cut;}
  Float_t GetCombinedCut() const {return fkCombinedCut;}
  Float_t GetMinClFullTrk() const {return fkMinClFullTrk;}
  Float_t GetMinTgl0() const {return fkMinTgl0;}
  Float_t GetMinRTgl0() const {return fkMinRTgl0;}

  Float_t GetMinNormDistForbTgl0() const {return fkMinNormDistForbTgl0;} 
  Float_t GetMinClForb0() const {return fkMinClForb0;}
  Float_t GetMinNormDistForb1() const {return fkMinNormDistForb1;}
  Float_t GetMinNormDistForb2() const {return fkMinNormDistForb2;}
  Float_t GetMinNormDistForb3() const {return fkMinNormDistForb3;}
  Float_t GetMinNormDistForb4() const {return fkMinNormDistForb4;}
  Float_t GetMinNormDistForb5() const {return fkMinNormDistForb5;}

  Float_t GetMinRTPCdensity() const  {return fkMinRTPCdensity;}
  Float_t GetMaxRTPCdensity0() const {return fkMaxRTPCdensity0;}
  Float_t GetMaxRTPCdensity10() const {return fkMaxRTPCdensity10;}
  Float_t GetMaxRTPCdensity20() const {return fkMaxRTPCdensity20;}
  Float_t GetMaxRTPCdensity30() const {return fkMaxRTPCdensity30;}

  Float_t GetMinTPCdensity() const {return fkMinTPCdensity;}
  Float_t GetMinTgl1() const {return fkMinTgl1;}
  Float_t GetMinTgl2() const {return fkMinTgl2;}
  Float_t GetMinchi2before0() const {return fkMinchi2before0;}
  Float_t GetMinchi2before1() const {return fkMinchi2before1;}
  Float_t GetMinchi2after0() const {return fkMinchi2after0;}
  Float_t GetMinchi2after1() const {return fkMinchi2after1;}
  Float_t GetAddchi2SharedCl() const{return fkAddchi2SharedCl;}
  Float_t GetAddchi2NegCl0() const{return fkAddchi2NegCl0;}
  Float_t GetAddchi2NegCl1() const{return fkAddchi2NegCl1;}
  Float_t GetSigp0Par0() const{return fkSigp0Par0;}
  Float_t GetSigp0Par1() const{return fkSigp0Par1;}
  Float_t GetSigp0Par2() const{return fkSigp0Par2;}
  Float_t GetSigpPar0() const{return fkSigpPar0;}
  Float_t GetSigpPar1() const{return fkSigpPar1;}
  Float_t GetSigpPar2() const{return fkSigpPar2;}
  Float_t GetMaxDcaLh0() const{return fkMaxDcaLh0;}

  Float_t GetChi2KF() const{return fkChi2KF;}
  Float_t GetRobustChi2KF() const{return fkRobustChi2KF;}

  Int_t StreamLevel()               { return fgStreamLevel;}
  void  SetStreamLevel(Int_t level=0) { fgStreamLevel = level;}
  

 private:
  Double_t  fPSigmaOffsetD0;        // sigma offset DCA
  Double_t  fPSigmaOffsetAP0;       // sigma offset AP
  // effective sigma DCA params    
  Double_t  fPSigmaMaxDE;           // maximal allowed sigma DCA
  Double_t  fPSigmaOffsetDE;        // offset sigma DCA
  Double_t  fPSigmaCoefDE;          // sigma coefiecient 
  Double_t  fPSigmaRminDE;          // max radius  - with momentum dependence 
  // effective sigma PA params
  Double_t  fPSigmaBase0APE;        // base sigma PA
  Double_t  fPSigmaMaxAPE;          // maximal sigma PA
  Double_t  fPSigmaR0APE;           // radial dependent part   - coeficient
  Double_t  fPSigmaR1APE;           // radial dependent part   - offset 
  Double_t  fPSigmaP0APE;           // momentum dependent part - coeficient
  Double_t  fPSigmaP1APE;           // momentum dependent part - offset
  // minimax parameters
  Double_t fPMinFractionAP0;        // minimal allowed fraction of effective params - PA
  Double_t fPMaxFractionAP0;        // maximal allowed fraction of effective params - PA
  Double_t fPMinAP0;                // minimal minimax - PA sigma 
  //
  Double_t fPMinFractionD0;         // minimal allowed fraction of effective params - DCA
  Double_t fPMaxFractionD0;         // maximal allowed fraction of effective params - DCA
  Double_t fPMinD0;                 // minimal minimax - DCA sigma
  //

  //Cuts for AliITSV0Finder


  Float_t fkMaxDist0;
  Float_t fkMaxDist1;
  Float_t fkMaxDist;
  Float_t fkMinPointAngle;
  Float_t fkMinPointAngle2;
  Float_t fkMinR;
  Float_t fkMaxR;
  Float_t fkMinPABestConst;
  Float_t fkMaxRBestConst;
  Float_t fkCausality0Cut;
  Float_t fkLikelihood01Cut;
  Float_t fkLikelihood1Cut;
  Float_t fkCombinedCut;
  Float_t fkMinClFullTrk;
  Float_t fkMinTgl0;

  Float_t fkMinClForb0;
  Float_t fkMinRTgl0; 
  Float_t fkMinNormDistForbTgl0;
  Float_t fkMinNormDistForb1;
  Float_t fkMinNormDistForb2;
  Float_t fkMinNormDistForb3;
  Float_t fkMinNormDistForb4;
  Float_t fkMinNormDistForb5;

  Float_t fkMinRTPCdensity;
  Float_t fkMaxRTPCdensity0;
  Float_t fkMaxRTPCdensity10;
  Float_t fkMaxRTPCdensity20;
  Float_t fkMaxRTPCdensity30;

  Float_t fkMinTPCdensity;
  Float_t fkMinTgl1;
  Float_t fkMinTgl2;
  Float_t fkMinchi2before0;
  Float_t fkMinchi2before1;
  Float_t fkMinchi2after0;
  Float_t fkMinchi2after1;
  Float_t fkAddchi2SharedCl;
  Float_t fkAddchi2NegCl0;
  Float_t fkAddchi2NegCl1;

  Float_t fkSigp0Par0;
  Float_t fkSigp0Par1;
  Float_t fkSigp0Par2;
  Float_t fkSigpPar0;
  Float_t fkSigpPar1;
  Float_t fkSigpPar2;
  Float_t fkMaxDcaLh0;

  Float_t fkChi2KF;
  Float_t fkRobustChi2KF;
  Int_t   fgStreamLevel; // flag for streaming - for ITS V0




  ClassDef(AliESDV0Params,2)      // ESD V0 vertex - error and likelihood parameterization constant
};



#endif
