#ifndef ALIUA1JETFINDER_H
#define ALIUA1JETFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//---------------------------------------------------------------------
// UA1 Cone Algorithm Finder 
// manages the search for jets
// Author: Rafael.Diaz.Valdes@cern.ch
// (version in c++)
// ** 2011
// Modified accordingly to reader/finder splitting and new handling of neutral information
// Versions V1 and V2 merged
//---------------------------------------------------------------------

#include "AliJetFinder.h"

class TH2F;
class AliJetBkg;

class AliUA1JetFinder : public AliJetFinder
{
 public:

  AliUA1JetFinder();
  ~AliUA1JetFinder();

  // others
  void FindJets();
  void RunAlgoritm(Float_t EtbgTotal, Double_t dEtTotal, Int_t& nJets,
		   Float_t* const etJet,Float_t* const etaJet, Float_t* const phiJet,
		   Float_t* const etallJet, Int_t* const ncellsJet);

  void Reset();
  void Init();
  void WriteJHeaderToFile() const;

  enum {kMaxJets = 60};

 protected:
  AliUA1JetFinder(const AliUA1JetFinder& rJetF1);
  AliUA1JetFinder& operator = (const AliUA1JetFinder& rhsf);

  TH2F*       fLego;          //  Lego Histo

  AliJetBkg*  fJetBkg;        //! pointer to bkg class

  ClassDef(AliUA1JetFinder,3) //  UA1 jet finder

};

#endif
