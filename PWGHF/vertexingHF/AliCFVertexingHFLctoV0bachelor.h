#ifndef ALICFVERTEXINGHFLCTOV0BACHELOR_H
#define ALICFVERTEXINGHFLCTOV0BACHELOR_H

/**************************************************************************
 * Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */ 

//-----------------------------------------------------------------------
// Class for HF corrections as a function of many variables and steps
// For Lc->V0 bachelor
// 
// Author : A. De Caro
//-----------------------------------------------------------------------


#include "AliCFVertexingHF.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"

class AliAODMCParticle;
class TClonesArray;
class AliCFVertexingHF;
class AliESDtrack;
class TDatabasePDG;

class AliCFVertexingHFLctoV0bachelor : public AliCFVertexingHF{
 public:

  enum ELctoV0Sel {
    kCountK0Sp=0,
    kCountLambdapi=1
  };

  AliCFVertexingHFLctoV0bachelor();
  AliCFVertexingHFLctoV0bachelor(TClonesArray *mcArray, UShort_t originDselection, Int_t lcDecay=0);

  virtual ~AliCFVertexingHFLctoV0bachelor(){};

  Bool_t GetGeneratedValuesFromMCParticle(Double_t* /*vectorMC*/);
  Bool_t GetRecoValuesFromCandidate(Double_t* /*vectorReco*/ ) const;
  Bool_t CheckMCChannelDecay()const;
  
  Bool_t SetRecoCandidateParam(AliAODRecoDecayHF *recoCand);

  Double_t GetEtaProng(Int_t iProng)const;
  Double_t GetPtProng(Int_t iProng) const;

  virtual Bool_t SetLabelArray();

 protected:

  Double_t Ctau(AliAODMCParticle *mcPartCandidate);
  Bool_t FillVectorFromMCarray(AliAODMCParticle *mcPartDaughterBachelor,
			       AliAODMCParticle *mcPartDaughterK0,
			       Double_t *vectorMC);

 private:	

  AliCFVertexingHFLctoV0bachelor(const AliCFVertexingHFLctoV0bachelor& c);
  AliCFVertexingHFLctoV0bachelor& operator= (const AliCFVertexingHFLctoV0bachelor& other);
  
  Int_t fGenLcOption;  // option for selection Lc (see enum)

  ClassDef(AliCFVertexingHFLctoV0bachelor, 1); // CF class for Lc->V0+bachelor and other cascades
  
};

#endif
