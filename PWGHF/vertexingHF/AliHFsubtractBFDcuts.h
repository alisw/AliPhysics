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
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"
class AliHFsubtractBFDcuts : public TNamed{

public:
  AliHFsubtractBFDcuts();
  AliHFsubtractBFDcuts(const char* name,const char* title);

  void InitHistos();

  void FillGenStep(AliAODMCParticle *dzeroMC,Double_t pt=-1);
  THnSparseF* GetSparseData(){return fcutsData;}
  THnSparseF* GetSparseMC(){return fcutsMC;}
  TH1F* GetHistoPtMCgen(){return fPtMCGenStep;}

  void SetHistoPtMCgen(TH1F *h){if(fPtMCGenStep)delete fPtMCGenStep;fPtMCGenStep=(TH1F*)h->Clone();return;}
  void SetSparseData(THnSparseF *h){if(fcutsData)delete fcutsData;fcutsData=(THnSparseF*)h->Clone();return;}
  void SetSparseMC(THnSparseF *h){if(fcutsMC)delete fcutsMC;fcutsMC=(THnSparseF*)h->Clone();return;}

  void SetFillMC (Bool_t fillMC = kTRUE) {fisMC = fillMC;}

  void FillSparses(AliAODRecoDecayHF2Prong *dzeroPart,Int_t isSelected,Double_t pt=-1,Double_t massD0=-1,Double_t massD0bar=-1);

private:

  AliHFsubtractBFDcuts(const AliHFsubtractBFDcuts &source);
  AliHFsubtractBFDcuts& operator=(const AliHFsubtractBFDcuts& source); 
 

  Bool_t   fisMC;        // flag for MC/Data
  TH1F *fPtMCGenStep;    //! histo with spectrum at generation level
  THnSparseF *fcutsData; //! THnSparse for cut variables (data, with inv mass axis), first axis is always mass
  THnSparseF *fcutsMC;   //! THnSparse for cut variables (MC at PID level, w/o mass axis)



  ClassDef(AliHFsubtractBFDcuts,1);
};

#endif
