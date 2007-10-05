#ifndef ALIFASTJETHEADER_H
#define ALIFASTJETHEADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// FastJet algorithm
//
// Author: Rafael.Diaz.Valdes@cern.ch
//---------------------------------------------------------------------
 
#include "AliJetHeader.h"

 
class AliFastJetHeader : public AliJetHeader
{
 public:
 
  AliFastJetHeader();
  virtual ~AliFastJetHeader() { }

  // Getters
  Int_t   GetLegoNbinEta()   const  {return fLegoNbinEta;}
  Int_t   GetLegoNbinPhi()   const  {return fLegoNbinPhi;}
  Float_t GetLegoEtaMin()    const  {return fLegoEtaMin;}
  Float_t GetLegoEtaMax()    const  {return fLegoEtaMax;}
  Float_t GetLegoPhiMin()    const  {return fLegoPhiMin;}
  Float_t GetLegoPhiMax()    const  {return fLegoPhiMax;}
  Float_t GetRadius()        const  {return fRadius;}
  Float_t GetMinJetEt()      const  {return fMinJetEt;}
  Bool_t  AddSoftBackg()     const  {return fSoftBackg;}
  Float_t GetPrecBg()        const  {return fPrecBg;}
  // Setters
  void SetLegoNbinEta(Int_t f) {fLegoNbinEta=f;}
  void SetLegoNbinPhi(Int_t f) {fLegoNbinPhi=f;}
  void SetLegoEtaMin(Float_t f) {fLegoEtaMin=f;}
  void SetLegoEtaMax(Float_t f) {fLegoEtaMax=f;}
  void SetLegoPhiMin(Float_t f) {fLegoPhiMin=f;}
  void SetLegoPhiMax(Float_t f) {fLegoPhiMax=f;}
  void SetRadius(Float_t f)     {fRadius = f;}
  void SetMinJetEt(Float_t f) {fMinJetEt=f;}
  void SetSoftBackg(Bool_t f) {fSoftBackg=f;}
  void SetPrecBg(Float_t f) {fPrecBg=f;}
  // others
  void PrintParameters() const;

 protected:
  // parameters for legos
  Int_t   fLegoNbinEta;         //! number of cells in eta
  Int_t   fLegoNbinPhi;         //! number of cells in phi
  Float_t fLegoEtaMin;          //! minimum eta
  Float_t fLegoEtaMax;          //! maximum eta
  Float_t fLegoPhiMin;          //! minimun phi
  Float_t fLegoPhiMax;          //! maximum phi
  //algorithm parameters
  Float_t fRadius;              // R parameter  (fastjet)
  Float_t fMinJetEt;            // Min Et of jet (geneal)
  Bool_t  fSoftBackg;           // add soft background
  Float_t fPrecBg;              // max value of change for BG (in %)


  ClassDef(AliFastJetHeader,1)
};
 
#endif
