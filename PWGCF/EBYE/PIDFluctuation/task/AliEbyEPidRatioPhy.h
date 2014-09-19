#ifndef ALIEBYEPIDRSTIOPHY_H
#define ALIEBYEPIDRSTIOPHY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    //
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//=========================================================================// 

#include "THnSparse.h"
#include "TList.h"

#include "AliEbyEPidRatioBase.h"

class AliEbyEPidRatioPhy : public AliEbyEPidRatioBase {

 public:

  AliEbyEPidRatioPhy();
  virtual ~AliEbyEPidRatioPhy();
  virtual void Process();
  void SetOutList(TList* l) {fOutList = l;}

 private:

  AliEbyEPidRatioPhy(const AliEbyEPidRatioPhy&);
  AliEbyEPidRatioPhy& operator=(const AliEbyEPidRatioPhy&); 
  virtual void Init();
  virtual void Reset();
  virtual void CreateHistograms();  
  Int_t ProcessTracks();
  Int_t ProcessParticles();  
  void ResetHistSet();
  void AddHistSetCent(const Char_t *name, const Char_t *title);
  void FillHistSetCent(const Char_t *name, Int_t idx, Bool_t isMC);

  void AddHistSetRatio(const Char_t *name, const Char_t *title);
  void FillHistSetRatio(const Char_t *name, Int_t idx, Bool_t isMC);
  
  
  TList                *fOutList;               //! Output data container
  Int_t                 fOrder;                 //  Max order of higher order distributions
  Int_t                 fNNp;                   //  N sets of arrays of particle/anti-particle counts
  Int_t              ***fNp;                    //  Array of particle/anti-particle counts
  Int_t                 fNMCNp;                 //  N sets of arrays of MC particle/anti-particle counts
  Int_t              ***fMCNp;                  //  Array of MC particle/anti-particle counts
  Double_t            **fRedFactp;              //  Array of particle/anti-particle reduced factorial

  ClassDef(AliEbyEPidRatioPhy, 1);
};

#endif
