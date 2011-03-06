#ifndef ALIUA1JETFINDERV2_H
#define ALIUA1JETFINDERV2_H

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
class TChain;

class AliUA1JetFinderV2 : public AliJetFinder
{
 public:

  AliUA1JetFinderV2();
  ~AliUA1JetFinderV2();

  // others
  void FindJetsC();
  void FindJets();
  void RunAlgoritmC(Float_t EtbgTotal, Double_t dEtTotal, Int_t& nJets,
		   Float_t* const etJet,Float_t* const etaJet, Float_t* const phiJet,
		   Float_t* const etallJet, Int_t* const ncellsJet);

  void RunAlgoritm(Int_t nIn, Float_t* etCell, Float_t* const etaCell, Float_t* phiCell, 
		   Int_t* const flagCell, const Float_t* etCell2, const Float_t* etaCell2, const Float_t* phiCell2, 
		   const Int_t* flagCell2, Float_t etbgTotal, Double_t dEtTotal, 
		   Int_t& nJets, Float_t* const etJet, Float_t* const etaJet, Float_t* const phiJet,
		   Float_t* const etallJet, Int_t* const ncellsJet);
   
  void SubtractBackgC(const Int_t& nIn, const Int_t&nJ, Float_t&EtbgTotalN,
                      const Float_t* ptT, const Float_t* etaT, const Float_t* phiT,
                      Float_t* const etJet, const Float_t* etaJet, const Float_t* phiJet,
                      Float_t* const etsigJet,Int_t* const multJet, Int_t* const injet);

  void SubtractBackg(const Int_t& nIn, const Int_t&nJ, Float_t&EtbgTotalN, const Float_t* ptT, const Int_t* vectT, 
		     const Float_t* etaT, const Float_t* phiT, const Int_t* cFlagT, const Int_t* cFlag2T, 
		     const Int_t* sFlagT, Float_t* const etJet, const Float_t* etaJet, const Float_t* phiJet, 
		     Float_t* const etsigJet, Int_t* const multJet, Int_t* const injet);

  void SubtractBackgCone(const Int_t& nIn, const Int_t&nJ,Float_t& EtbgTotalN,
                      const Float_t* ptT, const Float_t* etaT, const Float_t* phiT, const Int_t* cFlagT, const Int_t* sFlagT,
                      Float_t* const etJet, const Float_t* etaJet, const Float_t* phiJet,
                      Float_t* const etsigJet, Int_t* const multJet, Int_t* const injet);

  void SubtractBackgRatio(const Int_t& nIn, const Int_t&nJ,Float_t& EtbgTotalN,
                      const Float_t* ptT, const Float_t* etaT, const Float_t* phiT, const Int_t* cFlagT, const Int_t* sFlagT,
                      Float_t* const etJet, const Float_t* etaJet, const Float_t* phiJet,
                      Float_t* const etsigJet, Int_t* const multJet, Int_t* const injet);

  void SubtractBackgStat(const Int_t& nIn, const Int_t&nJ,Float_t&EtbgTotalN,
			 const Float_t* ptT, const Float_t* etaT, const Float_t* phiT, const Int_t* cFlagT, 
			 const Int_t* sFlagT, Float_t* const etJet, const Float_t* etaJet, const Float_t* phiJet,
			 Float_t* const etsigJet, Int_t* const multJet, Int_t* const injet);
  void Reset();
  void InitTask(TChain* tree);
  void WriteJHeaderToFile() const;
 protected:
  AliUA1JetFinderV2(const AliUA1JetFinderV2& rJetF1);
  AliUA1JetFinderV2& operator = (const AliUA1JetFinderV2& rhsf);
  TH2F           * fLego;           // Lego Histo
  Int_t fOpt;                       // Detector option (charged only or charged+neutral)

  ClassDef(AliUA1JetFinderV2,1)
};

#endif
