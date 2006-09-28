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
  Double_t GetEffMass() const {return fEffMass;}
  Double_t GetChi2V0() const {return fChi2V0;}
  void     GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void     GetNPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void     GetPPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  void     GetXYZ(Double_t &x, Double_t &y, Double_t &z) const;
  Double_t GetD(Double_t x0=0.,Double_t y0=0.,Double_t z0=0.) const;
  Int_t    GetNindex() const {return fNidx;}
  Int_t    GetPindex() const {return fPidx;}
  void     SetESDindexes(Int_t ip, Int_t im){fNidx=ip;fPidx=im;}
  void     SetDcaV0Daughters(Double_t rDcaV0Daughters=0.);
  Double_t GetDcaV0Daughters() {return fDcaV0Daughters;}
  Double_t GetV0CosineOfPointingAngle(Double_t&, Double_t&, Double_t&) const;
  void SetOnFlyStatus(Bool_t status){fOnFlyStatus=status;}
  Bool_t GetOnFlyStatus() const {return fOnFlyStatus;}

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
  const AliExternalTrackParam *GetParamP() const {return &fParamP;}
  const AliExternalTrackParam *GetParamM() const {return &fParamN;}
  static const AliESDV0Params & GetParameterization(){return fgkParams;}
  void SetP(const AliExternalTrackParam & paramp); 
  void SetM(const AliExternalTrackParam & paramd);
  void SetRp(const Double_t *rp);
  void SetRm(const Double_t *rm);
  void UpdatePID(Double_t pidp[5], Double_t pidm[5]);
  void SetStatus(Int_t status){fStatus=status;}
  Int_t GetStatus() const {return fStatus;}
  Float_t GetEffMass(UInt_t p1, UInt_t p2);
  Float_t GetProb(UInt_t p1, UInt_t p2);
  void SetID(Int_t id){fID =id;}
  Int_t GetID() const { return fID;}
  Int_t GetIndex(Int_t i) const {return fIndex[i];}
  void SetIndex(Int_t i, Int_t ind) {fIndex[i]=ind;}
  void SetDist1(Double_t d1) {fDist1=d1;}
  void SetDist2(Double_t d2) {fDist2=d2;}
  Double_t GetDist1() const {return fDist1;}
  Double_t GetDist2() const {return fDist2;}
  Double_t *GetAnglep() {return fAngle;}
  Double_t GetRr() const {return fRr;}
  void SetRr(Double_t rr) {fRr=rr;}
  Double_t *GetPMp() {return fPM;}
  Double_t *GetPPp() {return fPP;}
  Double_t *GetXrp() {return fXr;}
  Double_t GetXr(Int_t i) const {return fXr[i];}
  Double_t GetDistSigma() const {return fDistSigma;}
  void SetDistSigma(Double_t ds) {fDistSigma=ds;}
  Double_t GetDistNorm() const {return fDistNorm;}
  void SetDistNorm(Double_t ds) {fDistNorm=ds;}
  Float_t GetChi2Before() const {return fChi2Before;}
  void SetChi2Before(Float_t cb) {fChi2Before=cb;}
  Float_t GetChi2After() const {return fChi2After;}
  void SetChi2After(Float_t ca) {fChi2After=ca;}
  Float_t GetPointAngle() const {return fPointAngle;}
  void SetOrder(Int_t i, Int_t ord) {fOrder[i]=ord;}
  Float_t GetNAfter() const {return fNAfter;}
  void SetNAfter(Float_t na) {fNAfter=na;}
  Float_t GetNBefore() const {return fNBefore;}
  void SetNBefore(Float_t nb) {fNBefore=nb;}  
  Int_t GetLab(Int_t i) const {return fLab[i];}
  void  SetLab(Int_t i, Int_t lab) {fLab[i]=lab;}
  void SetCausality(Float_t pb0, Float_t pb1, Float_t pa0, Float_t pa1);
  const Float_t * GetCausalityP() const {return fCausality;}
  void SetClusters(Int_t *clp, Int_t *clm);
  const Int_t * GetClusters(Int_t i) const {return fClusters[i];}
  void SetNormDCAPrim(Float_t nd0, Float_t nd1){fNormDCAPrim[0] = nd0; fNormDCAPrim[1]=nd1;}
  const Float_t  *GetNormDCAPrimP() const {return fNormDCAPrim;}

protected:
  Bool_t   fOnFlyStatus;    // if kTRUE, then this V0 is recontructed
                            // "on fly" during the tracking

  Int_t    fPdgCode;        // reconstructed V0's type (PDG code)
  Double_t fEffMass;        // reconstructed V0's effective mass
  Double_t fDcaV0Daughters; // dca between V0's daughters
  Double_t fChi2V0;         // V0's chi2 value
  Double_t fPos[3];         // V0's position (global)
  Double_t fPosCov[6];      // covariance matrix of the vertex position

  Int_t    fNidx;           // index of the negative daughter
  Double_t fNmom[3];        // momentum of the negative daughter (global)
  Double_t fNmomCov[6];     // covariance matrix of the negative daughter mom.

  Int_t    fPidx;           // index of the positive daughter
  Double_t fPmom[3];        // momentum of the positive daughter (global)
  Double_t fPmomCov[6];     // covariance matrix of the positive daughter mom.

  // **** The following data members need to be revised ***

  AliExternalTrackParam fParamP;  // external parameters of positive particle
  AliExternalTrackParam fParamN;  // external parameters of negative particle
  Float_t        fRP[5];         // combined pid positive
  Float_t        fRM[5];         // combined pid positive
  Int_t          fID;            // ID number of the V0 in the ESDV0 container
  Int_t          fLab[2];         // MC label of the particle
  Int_t          fIndex[2];       // reconstructed labels of the tracks
  Int_t          fClusters[2][6]; //! its clusters 
  //
  //  
  Float_t        fNormDCAPrim[2];  // normalize distance to the priary vertex
  Double_t       fDist1;    //info about closest distance according closest MC - linear DCA
  Double_t       fDist2;    //info about closest distance parabolic DCA
  //
  Double_t       fPP[3];    //momentum  positive   - according approx at DCA
  Double_t       fPM[3];    //momentum negative
  //
  Double_t       fXr[3];      //rec. position according helix
  Double_t       fAngle[3];   //three angles
  Double_t       fRr;         //rec position of the vertex 
  Int_t          fStatus;       //status
  Int_t          fRow0;         // critical layer
  Int_t          fOrder[3]; //order of the vertex 
  //  quality information
  Double_t       fDistNorm; //normalized  DCA
  Double_t       fDistSigma; //sigma of distance
  Float_t        fCausality[4];  // causality information - see comments in SetCausality
  Float_t        fChi2Before;   //chi2 of the tracks before V0
  Float_t        fNBefore;      // number of possible points before V0
  Float_t        fChi2After;   // chi2 of the tracks after V0
  Float_t        fNAfter;      // number of possible points after V0
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Double_t       fPointAngle;   //point angle full
  //
  // parameterization coefficients
  static AliESDV0Params fgkParams;  // resolution and likelihood parameterization  

private:
  AliESDv0& operator=(const AliESDv0&);

  ClassDef(AliESDv0,2)      // ESD V0 vertex
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
