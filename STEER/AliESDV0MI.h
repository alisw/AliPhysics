#ifndef ALIESDV0MI_H
#define ALIESDV0MI_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          ESD V0 Vertex Class
//          This class is part of the Event Summary Data set of classes
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------

#include "AliESDv0.h"
#include "AliESDV0MIParams.h"
#include "AliExternalTrackParam.h"

class AliESDtrack;


class AliESDV0MI :  public AliESDv0 {
public:
  //  friend class AliITStrackerMI;
  AliESDV0MI();             //constructor
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
  const AliExternalTrackParam *GetParamM() const {return &fParamM;}
  static const AliESDV0MIParams & GetParameterization(){return fgkParams;}
  void SetP(const AliExternalTrackParam & paramp); 
  void SetM(const AliExternalTrackParam & paramd);
  void SetRp(const Double_t *rp);
  void SetRm(const Double_t *rm);
  void UpdatePID(Double_t pidp[5], Double_t pidm[5]);
  void SetStatus(Int_t status){fStatus=status;}
  Int_t  GetStatus() const {return fStatus;}
  Float_t GetEffMass(UInt_t p1, UInt_t p2);
  Float_t GetProb(UInt_t p1, UInt_t p2);
  void Update(Float_t vertex[3]);            //update
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
  void SetLab(Int_t i, Int_t lab) {fLab[i]=lab;}
  void SetCausality(Float_t pb0, Float_t pb1, Float_t pa0, Float_t pa1);
  const Float_t * GetCausalityP() const {return fCausality;}
  void SetClusters(Int_t *clp, Int_t *clm);
  const Int_t * GetClusters(Int_t i) const {return fClusters[i];}
  void SetNormDCAPrim(Float_t nd0, Float_t nd1){fNormDCAPrim[0] = nd0; fNormDCAPrim[1]=nd1;}
  const Float_t  *GetNormDCAPrimP() const {return fNormDCAPrim;}
private:
  AliExternalTrackParam fParamP;  // external parameters of positive particle
  AliExternalTrackParam fParamM;  // external parameters of negative particle
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
  static AliESDV0MIParams fgkParams;  // resolution and likelihood parameterization  
  ClassDef(AliESDV0MI,4)      // ESD V0 vertex
};


#endif
