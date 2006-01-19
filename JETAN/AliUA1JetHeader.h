#ifndef ALIUA1JETHEADER_H
#define ALIUA1JETHEADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet header class for UA1 algorithm
// Stores the parameters of the UA1 jet algorithm
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------
 
#include "AliJetHeader.h"
 
class AliUA1JetHeader : public AliJetHeader
{
 public:
 
  AliUA1JetHeader();
  virtual ~AliUA1JetHeader() { }

  // Getters
  Float_t GetRadius()    const  {return fConeRadius;}
  Float_t GetEtSeed()    const  {return fEtSeed;}
  Float_t GetMinJetEt()  const  {return fMinJetEt;}
  Float_t GetMinCellEt() const  {return fMinCellEt;}
  Int_t   GetMode()      const  {return fMode;}
  Float_t GetMinMove()   const  {return fMinMove;}
  Float_t GetMaxMove()   const  {return fMaxMove;}
  Float_t GetPrecBg()    const  {return fPrecBg;}
  Int_t   GetLegoNbinEta()   const  {return fLegoNbinEta;}
  Int_t   GetLegoNbinPhi()   const  {return fLegoNbinPhi;}
  Float_t GetLegoEtaMin()    const  {return fLegoEtaMin;}
  Float_t GetLegoEtaMax()    const  {return fLegoEtaMax;}
  Float_t GetLegoPhiMin()    const  {return fLegoPhiMin;}
  Float_t GetLegoPhiMax()    const  {return fLegoPhiMax;}
  Bool_t GetOnlySignal()     const {return fOnlySignal;}
  Bool_t GetOnlyBkgd()     const {return fOnlyBkgd;}

  // Setters
  void SetRadius(Float_t f) {fConeRadius=f;}
  void SetEtSeed(Float_t f) {fEtSeed=f;}
  void SetMinJetEt(Float_t f) {fMinJetEt=f;}
  void SetMinCellEt(Float_t f) {fMinCellEt=f;}
  void SetMode(Int_t f) {fMode=f;}
  void SetMinMove(Float_t f) {fMinMove=f;}
  void SetMaxMove(Float_t f) {fMaxMove=f;}
  void SetPrecBg(Float_t f) {fPrecBg=f;}
  void SetLegoNbinEta(Int_t f) {fLegoNbinEta=f;}
  void SetLegoNbinPhi(Int_t f) {fLegoNbinPhi=f;}
  void SetLegoEtaMin(Float_t f) {fLegoEtaMin=f;}
  void SetLegoEtaMax(Float_t f) {fLegoEtaMax=f;}
  void SetLegoPhiMin(Float_t f) {fLegoPhiMin=f;}
  void SetLegoPhiMax(Float_t f) {fLegoPhiMax=f;}
  void SetOnlySignal(Bool_t b) {fOnlySignal= b;}
  void SetOnlyBkgd(Bool_t b) {fOnlyBkgd= b;}

  // others
  void PrintParameters() const; 

protected:

  // parameters of algorithm
  Float_t fConeRadius;      //  Cone radius
  Float_t fEtSeed;          //  Min. Et for seed
  Float_t fMinJetEt;        //  Min Et of jet
  Float_t fMinCellEt;       //  Min Et in one cell
  // parameters of backgound substraction
  Int_t   fMode;            // key for BG subtraction
  Float_t fMinMove;         // min cone move 
  Float_t fMaxMove;         // max cone move
  Float_t fPrecBg;          // max value of change for BG (in %)
  // parameters for legos
  Int_t   fLegoNbinEta;         //! number of cells in eta
  Int_t   fLegoNbinPhi;         //! number of cells in phi
  Float_t fLegoEtaMin;          //! minimum eta  
  Float_t fLegoEtaMax;          //! maximum eta
  Float_t fLegoPhiMin;          //! minimun phi
  Float_t fLegoPhiMax;          //! maximum phi
  // parameters for Jet search
  Bool_t fOnlySignal; // use only signal
  Bool_t fOnlyBkgd;   // use only background

  ClassDef(AliUA1JetHeader,1)
};
 
#endif
