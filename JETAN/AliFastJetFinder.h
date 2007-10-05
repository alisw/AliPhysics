#ifndef ALIFASTJETFINDER_H
#define ALIFASTJETFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//---------------------------------------------------------------------
// Fast Jet finder algorithm interface
// manages the search for jets
// Author: Rafael.Diaz.Valdes@cern.ch
// kt algorithm using NlnN
//---------------------------------------------------------------------

#include "AliJetFinder.h"
class AliFastJetHeader;

class AliFastJetFinder : public AliJetFinder
{
 public:

  AliFastJetFinder();
  ~AliFastJetFinder();

  // getters

  // setters
  void SetJetHeader(AliFastJetHeader* h) {fHeader= h;}
  // others
  void FindJets();
  void RunAlgorithm(Int_t& nJets,Float_t* etJet,Float_t* etaJet,Float_t* phiJet,
                    Float_t* etsigJet, Float_t* etallJet, Int_t* ncellsJet,
                    Int_t& nCell,Float_t* etCell,Float_t* etaCell,Float_t* phiCell,
                    Float_t* etsigCell, Int_t* jetflagCell);
  void SubtractBackg(Int_t& nCell, Int_t* jetflagCell, Float_t* etCell,
                     Int_t& nJets, Float_t* etJet, Float_t* etallJet, Int_t* ncellsJet,
                     Float_t& meanptCell, Float_t& sqptCell, Float_t& etBackg);
  void SubtractBackgArea(Int_t& nCell, Int_t* jetflagCell, Float_t* etCell,
                     Int_t& nJets, Float_t* etJet, Float_t* etallJet);
  void Reset();
  void Init();
  void WriteJHeaderToFile();

 protected:

  AliFastJetHeader* fHeader;         // pointer to jet header
  TH2F*             fLego;           //! Lego Histo
  TH2F*             fLegoSignal;    // ! Lego histogram for signal

  ClassDef(AliFastJetFinder,1)
};

#endif
