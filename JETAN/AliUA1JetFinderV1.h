#ifndef ALIUA1JETFINDERV1_H
#define ALIUA1JETFINDERV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//---------------------------------------------------------------------
// UA1 Cone Algorithm Finder V1
// manages the search for jets
// Author: Rafael.Diaz.Valdes@cern.ch
// (version in c++)
//---------------------------------------------------------------------

#include "AliJetFinder.h"
class AliUA1JetHeaderV1;
class TH2F;

class AliUA1JetFinderV1 : public AliJetFinder
{
 public:

  AliUA1JetFinderV1();
  ~AliUA1JetFinderV1();

  // others
  void FindJets();
  void RunAlgoritm(Float_t EtbgTotal, Double_t dEtTotal, Int_t& nJets,
                   Float_t* etJet,Float_t* etaJet, Float_t* phiJet,
                   Float_t* etallJet, Int_t* ncellsJet);

  void SubtractBackg(Int_t& nIn, Int_t&nJ, Float_t&EtbgTotalN,
                      Float_t* ptT, Float_t* etaT, Float_t* phiT,
                      Float_t* etJet,Float_t* etaJet, Float_t* phiJet,
                      Float_t* etsigJet,Int_t* multJet, Int_t* injet);

  void SubtractBackgCone(Int_t& nIn, Int_t&nJ,Float_t& EtbgTotalN,
                      Float_t* ptT, Float_t* etaT, Float_t* phiT,
                      Float_t* etJet,Float_t* etaJet, Float_t* phiJet,
                      Float_t* etsigJet, Int_t* multJet, Int_t* injet);

  void SubtractBackgRatio(Int_t& nIn, Int_t&nJ,Float_t& EtbgTotalN,
                      Float_t* ptT, Float_t* etaT, Float_t* phiT,
                      Float_t* etJet,Float_t* etaJet, Float_t* phiJet,
                      Float_t* etsigJet, Int_t* multJet, Int_t* injet);

  void SubtractBackgStat(Int_t& nIn, Int_t&nJ,Float_t&EtbgTotalN,
                      Float_t* ptT, Float_t* etaT, Float_t* phiT,
                      Float_t* etJet,Float_t* etaJet, Float_t* phiJet,
                      Float_t* etsigJet, Int_t* multJet, Int_t* injet);
  void Reset();
  void Init();
  void WriteJHeaderToFile();

 protected:
  AliUA1JetFinderV1(const AliUA1JetFinderV1& rJetF1);
  AliUA1JetFinderV1& operator = (const AliUA1JetFinderV1& rhsf);
  TH2F           * fLego;           //! Lego Histo

  ClassDef(AliUA1JetFinderV1,1)
};

#endif
