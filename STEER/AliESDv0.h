#ifndef ALIESDV0_H
#define ALIESDV0_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          ESD V0 Vertex Class
//          This class is part of the Event Summary Data set of classes
//    Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//    Modified by: Marian Ivanov,  CERN, Marian.Ivanov@cern.ch
//            and  Boris Hippolyte,IPHC, hippolyt@in2p3.fr 
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TPDGCode.h>
#include "AliESDV0Params.h"
#include "AliExternalTrackParam.h"

class AliESDv0 : public TObject {
public:
  AliESDv0();
  AliESDv0(const AliExternalTrackParam &t1, Int_t i1,
           const AliExternalTrackParam &t2, Int_t i2);

  AliESDv0(const AliESDv0&);
  virtual ~AliESDv0();

  Double_t ChangeMassHypothesis(Int_t code=kK0Short); 

  Int_t    GetPdgCode() const {return fPdgCode;}
  Float_t  GetEffMass() const {return fEffMass;}
  Float_t  GetChi2V0()  const {return fChi2V0;}
  void     GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void     GetNPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void     GetPPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void     GetXYZ(Double_t &x, Double_t &y, Double_t &z) const;
  Float_t  GetD(Double_t x0=0.,Double_t y0=0.,Double_t z0=0.) const;
  Int_t    GetNindex() const {return fNidx;}
  Int_t    GetPindex() const {return fPidx;}
  void     SetDcaV0Daughters(Double_t rDcaV0Daughters=0.);
  Float_t  GetDcaV0Daughters() {return fDcaV0Daughters;}
  Float_t  GetV0CosineOfPointingAngle(Double_t&, Double_t&, Double_t&) const;
  Float_t  GetV0CosineOfPointingAngle() const {return fPointAngle;}
  void     SetV0CosineOfPointingAngle(Double_t cpa) {fPointAngle=cpa;}
  void     SetOnFlyStatus(Bool_t status){fOnFlyStatus=status;}
  Bool_t   GetOnFlyStatus() const {return fOnFlyStatus;}
  const AliExternalTrackParam *GetParamP() const {return &fParamP;}
  const AliExternalTrackParam *GetParamN() const {return &fParamN;}



  // **** The following member functions need to be revised ***

  Double_t GetSigmaY();     // sigma of y coordinate at vertex posistion
  Double_t GetSigmaZ();     // sigma of z coordinate at vertex posistion
  Double_t GetSigmaAP0();   // calculate sigma of Point angle resolution at vertex pos.
  Double_t GetSigmaD0();    // calculate sigma of position resolution at vertex pos.
  Double_t GetEffectiveSigmaAP0();   // calculate sigma of point angle resolution at vertex pos. effecive parameterization
  Double_t GetEffectiveSigmaD0();    // calculate sigma of position resolution at vertex pos.
  Double_t GetMinimaxSigmaAP0();    // calculate mini-max sigma of point angle resolution
  Double_t GetMinimaxSigmaD0();     // calculate mini-max sigma of dca resolution
  Double_t GetLikelihoodAP(Int_t mode0, Int_t mode1);   // get likelihood for point angle
  Double_t GetLikelihoodD(Int_t mode0, Int_t mode1);    // get likelihood for DCA
  Double_t GetLikelihoodC(Int_t mode0, Int_t mode1);    // get likelihood for Causality
  //
  //
  static const AliESDV0Params & GetParameterization(){return fgkParams;}
  void SetParamP(const AliExternalTrackParam & paramP) {fParamP = paramP;}
  void SetParamN(const AliExternalTrackParam & paramN) {fParamN = paramN;}
  void SetStatus(Int_t status){fStatus=status;}
  Int_t GetStatus() const {return fStatus;}
  Int_t GetIndex(Int_t i) const {return (i==0) ? fNidx : fPidx;}
  void SetIndex(Int_t i, Int_t ind) {(i==0) ? (fNidx=ind) : (fPidx=ind);}
  Double_t *GetAnglep() {return fAngle;}
  Double_t GetRr() const {return fRr;}
  Double_t GetDistSigma() const {return fDistSigma;}
  void SetDistSigma(Double_t ds) {fDistSigma=ds;}
  Float_t GetChi2Before() const {return fChi2Before;}
  void SetChi2Before(Float_t cb) {fChi2Before=cb;}
  Float_t GetChi2After() const {return fChi2After;}
  void SetChi2After(Float_t ca) {fChi2After=ca;}
  Float_t GetNAfter() const {return fNAfter;}
  void SetNAfter(Float_t na) {fNAfter=na;}
  Float_t GetNBefore() const {return fNBefore;}
  void SetNBefore(Float_t nb) {fNBefore=nb;}  
  void SetCausality(Float_t pb0, Float_t pb1, Float_t pa0, Float_t pa1);
  const Float_t * GetCausalityP() const {return fCausality;}
  void SetClusters(Int_t *clp, Int_t *clm);
  const Int_t * GetClusters(Int_t i) const {return fClusters[i];}
  void SetNormDCAPrim(Float_t nd0, Float_t nd1){fNormDCAPrim[0] = nd0; fNormDCAPrim[1]=nd1;}
  const Float_t  *GetNormDCAPrimP() const {return fNormDCAPrim;}

protected:
  Bool_t   fOnFlyStatus;    // if kTRUE, then this V0 is recontructed
                            // "on fly" during the tracking

  Int_t    fPdgCode;          // reconstructed V0's type (PDG code)
  Float_t  fEffMass;          // reconstructed V0's effective mass
  Float_t  fDcaV0Daughters;   // dca between V0's daughters
  Float_t  fPointAngle;       //cosine of the pointing angle
  Float_t  fChi2V0;           // V0's chi2 value

  Double32_t fPos[3];         // V0's position (global)
  Double32_t fPosCov[6];      // covariance matrix of the vertex position

  Int_t fNidx;                // index of the negative daughter
  Double32_t fNmom[3];        // momentum of the negative daughter (global)
  AliExternalTrackParam fParamN;  // external parameters of negative particle
  Int_t fPidx;                // index of the positive daughter
  Double32_t fPmom[3];        // momentum of the positive daughter (global)
  AliExternalTrackParam fParamP;  // external parameters of positive particle


  // **** The following data members need to be revised ***

  Int_t          fClusters[2][6]; //! its clusters 
  //
  Float_t        fNormDCAPrim[2];  // normalize distance to the priary vertex
  //
  Double32_t     fAngle[3];   //three angles
  Float_t        fRr;         //rec position of the vertex 
  Int_t          fStatus;       //status
  Float_t        fDistSigma; //sigma of distance
  Float_t        fCausality[4];  // causality information - see comments in SetCausality
  Float_t        fChi2Before;   //chi2 of the tracks before V0
  Float_t        fNBefore;      // number of possible points before V0
  Float_t        fChi2After;   // chi2 of the tracks after V0
  Float_t        fNAfter;      // number of possible points after V0
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  //
  // parameterization coefficients
  static AliESDV0Params fgkParams;  // resolution and likelihood parameterization  

private:
  AliESDv0& operator=(const AliESDv0&);

  ClassDef(AliESDv0,3)      // ESD V0 vertex
};

inline 
void AliESDv0::GetNPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
px=fNmom[0]; py=fNmom[1]; pz=fNmom[2];
}

inline 
void AliESDv0::GetPPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
px=fPmom[0]; py=fPmom[1]; pz=fPmom[2];
}

inline
void AliESDv0::SetDcaV0Daughters(Double_t rDcaV0Daughters){
  fDcaV0Daughters=rDcaV0Daughters;
}

#endif
