#ifndef ALICFVERTEXINGHFCASCADE_H
#define ALICFVERTEXINGHFCASCADE_H

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

//-----------------------------------------------------------------------
// Class for HF corrections as a function of many variables and steps
// For D* and other cascades
// 
// Author : A.GRELLI - a.grelli@uu.nl UTRECHT
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

class AliCFVertexingHFCascade : public AliCFVertexingHF{
 public:
		
  AliCFVertexingHFCascade(){};
  AliCFVertexingHFCascade(TClonesArray *mcArray, UShort_t originDselection);
	
  //virtual ~AliCFVertexingHFCascade(){};
  
  
  Bool_t GetGeneratedValuesFromMCParticle(Double_t* /*vectorMC*/);
  Bool_t GetRecoValuesFromCandidate(Double_t* /*vectorReco*/ ) const;
  Bool_t CheckMCChannelDecay()const;
  
  Bool_t SetRecoCandidateParam(AliAODRecoDecayHF *recoCand);
  Bool_t EvaluateIfD0toKpi(AliAODMCParticle* neutralDaugh, Double_t* VectorD0)const;

  void SetPtAccCut(Float_t* ptAccCut);
  void SetEtaAccCut(Float_t* etaAccCut);
  void SetAccCut(Float_t* ptAccCut, Float_t* etaAccCut);
  void SetAccCut();

 protected:
  
  
 private:	
  AliCFVertexingHFCascade(const AliCFVertexingHFCascade& c);
  AliCFVertexingHFCascade& operator= (const AliCFVertexingHFCascade& other);
  
  ClassDef(AliCFVertexingHFCascade, 1); // CF class for D* and other cascades
  
};

#endif
