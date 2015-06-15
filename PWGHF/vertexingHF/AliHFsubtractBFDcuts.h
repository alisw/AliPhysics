#ifndef ALIHFSUBTRACTBFDCUTS_H
#define ALIHFSUBTRACTBFDCUTS_H
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Class for storing and handling D0 meson candidates properties         //
//  for estimating the feed-down fraction using several sets of cuts      //
//     Andrea Rossi <andrea.rossi@cern.ch>                                //
//     Felix Reidt  <felix.reidt@cern.ch>                                 //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
#include "TNamed.h"
#include "THnSparse.h"
#include "TH3F.h"

#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"

class AliHFsubtractBFDcuts : public TNamed{

public:
  AliHFsubtractBFDcuts();
  AliHFsubtractBFDcuts(const char* name,const char* title);

  void InitHistos();

  void FillGenStep(AliAODMCParticle *dzeroMC,Double_t pt=-1,Double_t weight=1.,TClonesArray* mcArray=0x0);
  THnSparseF* GetSparseData() const {return fCutsData;}
  THnSparseF* GetSparseMC()   const {return fCutsMC;}
  TH3F* GetHistoPtMCgen()     const {return fPtMCGenStep;}

  void SetHistoPtMCgen(TH3F *h){if(fPtMCGenStep)delete fPtMCGenStep;fPtMCGenStep=(TH3F*)h->Clone();return;}
  void SetSparseData(THnSparseF *h){if(fCutsData)delete fCutsData;fCutsData=(THnSparseF*)h->Clone();return;}
  void SetSparseMC(THnSparseF *h){if(fCutsMC)delete fCutsMC;fCutsMC=(THnSparseF*)h->Clone();return;}

  void SetFillMC (Bool_t fillMC = kTRUE) {fIsMC = fillMC;}

  void FillSparses(AliAODRecoDecayHF2Prong *dzeroPart,Int_t isSelected,Double_t pt=-1,Double_t massD0=-1,Double_t massD0bar=-1,Double_t weight=1.,TClonesArray* mcArray=0x0);

private:
  AliHFsubtractBFDcuts(const AliHFsubtractBFDcuts& c);
  AliHFsubtractBFDcuts operator=(const AliHFsubtractBFDcuts& c);

  void GetCandidateLabel(AliAODRecoDecayHF2Prong *dzerocand);
  Bool_t AnalyseDecay();                 // check in which decay process a particle was created
  void   CountProngs(Int_t labCurrMother, Int_t labCurrExcl);  // counting the prongs of labCurrMother,
                                                               // labCurrExcl is assumed to be a stable particle
  Bool_t IsStable(Int_t labProng) const;       // Is that prong a stable particle?
  Bool_t IsInAcceptance(Int_t labProng) const; // Is that prong within the fiducial acceptance

  Bool_t      fIsMC;              // flag for MC/Data
  Bool_t      fCheckAcceptance;   // flag for checking whether the decay prongs are within acceptance
  Bool_t      fResolveResonances; // flag resolve resonances in during the prong determination
  TH3F*       fPtMCGenStep;       //! histo with spectrum at generation level
  THnSparseF* fCutsData;          //! THnSparse for cut variables (data, with inv mass axis), first axis is always mass
  THnSparseF* fCutsMC;            //! THnSparse for cut variables (MC at PID level, w/o mass axis)

  // Event specific variables
  TClonesArray* fMCarray;    //! TClonesArray holding the particles of the event to be processed
  Int_t         fLabCand;    // Label of the candidate D0 (charmed hadron in case of a chained decay)
  Int_t         fLabMother;  // Label of the mother of the candidate D0 (or charmed hadron)
  UInt_t        fNprongs;    // Number of prongs, counting the first charmed hadron as one particle
  Bool_t        fDecayChain; // Chained decay of charmed hadrons
  Double_t      fMotherPt;   // Transverse momentum of the mother particle (B hadron in case of feed-down,
                             // the charmed hadron itsself in case of prompt production)

  ClassDef(AliHFsubtractBFDcuts,5);
};

#endif
