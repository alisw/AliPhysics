#ifndef ALIUA1JETHEADERV1_H
#define ALIUA1JETHEADERV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//---------------------------------------------------------------------
// Jet Finder header class for algoritm using particles
// Stores the parameters of the particles algoritm
// Author: Rafael.Diaz.Valdes@cern.ch
//---------------------------------------------------------------------

#include "AliJetHeader.h"


class AliUA1JetHeaderV1 : public AliJetHeader
{
 public:

  AliUA1JetHeaderV1();
  virtual ~AliUA1JetHeaderV1() { }

  // Getters
  Float_t GetRadius()    const  {return fConeRadius;}
  Float_t GetMinMove()   const  {return fMinMove;}
  Float_t GetMaxMove()   const  {return fMaxMove;}
  Float_t GetEtSeed()    const  {return fEtSeed;}
  Float_t GetMinJetEt()  const  {return fMinJetEt;}
  Int_t   GetLegoNbinEta()   const  {return fLegoNbinEta;}
  Int_t   GetLegoNbinPhi()   const  {return fLegoNbinPhi;}
  Float_t GetLegoEtaMin()    const  {return fLegoEtaMin;}
  Float_t GetLegoEtaMax()    const  {return fLegoEtaMax;}
  Float_t GetLegoPhiMin()    const  {return fLegoPhiMin;}
  Float_t GetLegoPhiMax()    const  {return fLegoPhiMax;}
  Int_t   GetBackgMode()     const  {return fBackgMode;}
  Float_t GetPrecBg()    const  {return fPrecBg;}
  Float_t GetBackgStat()    const  {return fBackgStat;}
  Float_t GetBackgCutRatio()    const  {return fBackgCutRatio;}
  Int_t   GetNAcceptJets()     const  {return fNAcceptJets;}


  // Setters
  void SetRadius(Float_t f) {fConeRadius=f;}
  void SetMinMove(Float_t f) {fMinMove=f;}
  void SetMaxMove(Float_t f) {fMaxMove=f;}
  void SetEtSeed(Float_t f) {fEtSeed=f;}
  void SetMinJetEt(Float_t f) {fMinJetEt=f;}
  void SetLegoNbinEta(Int_t f) {fLegoNbinEta=f;}
  void SetLegoNbinPhi(Int_t f) {fLegoNbinPhi=f;}
  void SetLegoEtaMin(Float_t f) {fLegoEtaMin=f;}
  void SetLegoEtaMax(Float_t f) {fLegoEtaMax=f;}
  void SetLegoPhiMin(Float_t f) {fLegoPhiMin=f;}
  void SetLegoPhiMax(Float_t f) {fLegoPhiMax=f;}
  void BackgMode(Int_t mod ) {fBackgMode = mod;}
  void SetPrecBg(Float_t f) {fPrecBg=f;}
  void SetBackgStat(Float_t f) {fBackgStat=f;}
  void SetBackgCutRatio(Float_t f) {fBackgCutRatio=f;}
  void SetNAcceptJets(Int_t ajets)   {fNAcceptJets = ajets;}

  // others
  void PrintParameters() const;

protected:

  // parameters of algorithm
  Float_t fConeRadius;      //  Cone radius
  Float_t fEtSeed;          //  Min. Et for seed
  Float_t fMinJetEt;        //  Min Et of jet
  // parameters of backgound substraction
  Float_t fMinMove;         // min cone move
  Float_t fMaxMove;         // max cone move
  Int_t   fBackgMode;       // background subtraction mode
  Float_t fPrecBg;          // max value of change for BG (in %)
  Float_t fBackgStat;       // pre-calculated background used in statistic subtraction method
  Float_t fBackgCutRatio;   // pre-calculated pt-cut ratio used in ratio subtraction method
  Int_t   fNAcceptJets;     // number of accepted jets per events

  // parameters for legos
  Int_t   fLegoNbinEta;         // number of cells in eta
  Int_t   fLegoNbinPhi;         // number of cells in phi
  Float_t fLegoEtaMin;          // minimum eta
  Float_t fLegoEtaMax;          // maximum eta
  Float_t fLegoPhiMin;          // minimun phi
  Float_t fLegoPhiMax;          // maximum phi

  ClassDef(AliUA1JetHeaderV1,1)
};

#endif
