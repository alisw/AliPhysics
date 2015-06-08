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
#include <THnSparse.h>
#include <TH1F.h>
#include <TNamed.h>
#include <TClonesArray.h>
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"

class AliHFsubtractBFDcuts : public TNamed{

public:
  AliHFsubtractBFDcuts();
  AliHFsubtractBFDcuts(const char* name,const char* title);

  void InitHistos();

  void FillGenStep(AliAODMCParticle *dzeroMC,Double_t pt=-1,Double_t weight=1.);
  THnSparseF* GetSparseData() const {return fCutsData;}
  THnSparseF* GetSparseMC()   const {return fCutsMC;}
  TH1F* GetHistoPtMCgen()     const {return fPtMCGenStep;}

  void SetHistoPtMCgen(TH1F *h){if(fPtMCGenStep)delete fPtMCGenStep;fPtMCGenStep=(TH1F*)h->Clone();return;}
  void SetSparseData(THnSparseF *h){if(fCutsData)delete fCutsData;fCutsData=(THnSparseF*)h->Clone();return;}
  void SetSparseMC(THnSparseF *h){if(fCutsMC)delete fCutsMC;fCutsMC=(THnSparseF*)h->Clone();return;}

  void SetFillMC (Bool_t fillMC = kTRUE) {fIsMC = fillMC;}

  void FillSparses(AliAODRecoDecayHF2Prong *dzeroPart,Int_t isSelected,Double_t pt=-1,Double_t massD0=-1,Double_t massD0bar=-1,Double_t weight=1.,TClonesArray* mcArray=0x0);

private:
  Int_t GetCandidateLabel(AliAODRecoDecayHF2Prong *dzerocand,TClonesArray* mcArray) const;
  Bool_t DetermineDecayType(Int_t labCand,TClonesArray* mcArrray, UInt_t& nProngs, Bool_t& decayChain, Double_t& motherPt) const; // check in which decay process a particle was created

  Bool_t   fIsMC;        // flag for MC/Data
  TH1F *fPtMCGenStep;    //! histo with spectrum at generation level
  THnSparseF *fCutsData; //! THnSparse for cut variables (data, with inv mass axis), first axis is always mass
  THnSparseF *fCutsMC;   //! THnSparse for cut variables (MC at PID level, w/o mass axis)

  ClassDef(AliHFsubtractBFDcuts,2);
};

#endif
