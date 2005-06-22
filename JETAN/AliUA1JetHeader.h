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
  Int_t   GetNbinEta()   const  {return fNbinEta;}
  Int_t   GetNbinPhi()   const  {return fNbinPhi;}
  Float_t GetEtaMin()    const  {return fEtaMin;}
  Float_t GetEtaMax()    const  {return fEtaMax;}
  Float_t GetPhiMin()    const  {return fPhiMin;}
  Float_t GetPhiMax()    const  {return fPhiMax;}

  // Setters
  void SetRadius(Float_t f) {fConeRadius=f;}
  void SetEtSeed(Float_t f) {fEtSeed=f;}
  void SetMinJetEt(Float_t f) {fMinJetEt=f;}
  void SetMinCellEt(Float_t f) {fMinCellEt=f;}
  void SetMode(Int_t f) {fMode=f;}
  void SetMinMove(Float_t f) {fMinMove=f;}
  void SetMaxMove(Float_t f) {fMaxMove=f;}
  void SetPrecBg(Float_t f) {fPrecBg=f;}
  void SetNbinEta(Int_t f) {fNbinEta=f;}
  void SetNbinPhi(Int_t f) {fNbinPhi=f;}
  void SetEtaMin(Float_t f) {fEtaMin=f;}
  void SetEtaMax(Float_t f) {fEtaMax=f;}
  void SetPhiMin(Float_t f) {fPhiMin=f;}
  void SetPhiMax(Float_t f) {fPhiMax=f;}

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
  Int_t   fNbinEta;         //! number of cells in eta
  Int_t   fNbinPhi;         //! number of cells in phi
  Float_t fEtaMin;          //! minimum eta  
  Float_t fEtaMax;          //! maximum eta
  Float_t fPhiMin;          //! minimun phi
  Float_t fPhiMax;          //! maximum phi

  ClassDef(AliUA1JetHeader,1)
};
 
#endif
