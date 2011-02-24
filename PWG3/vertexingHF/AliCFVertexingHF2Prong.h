#ifndef ALICFVERTEXINGHF2PRONG_H
#define ALICFVERTEXINGHF2PRONG_H

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

/* $Id$ */ 

//-----------------------------------------------------------------------
// Class for HF corrections as a function of many variables and step 
// Author : C. Zampolli, CERN
// D. Caffarri, Univ & INFN Padova caffarri@pd.infn.it
// Base class for HF Unfolding - agrelli@uu.nl
//-----------------------------------------------------------------------


#include "AliCFVertexingHF.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"

class AliAODMCParticle;
class TClonesArray;
class AliCFVertexingHF;
class AliESDtrack;
class TDatabasePDG;

class AliCFVertexingHF2Prong : public AliCFVertexingHF{
	public:
		
	AliCFVertexingHF2Prong(){};
	AliCFVertexingHF2Prong(TClonesArray *mcArray, UShort_t originDselection);
	
	//  virtual ~AliCFVertexingHF2Prong(){};
	
 	Bool_t GetGeneratedValuesFromMCParticle(Double_t* /*vectorMC*/);
	Bool_t GetRecoValuesFromCandidate(Double_t* /*vectorReco*/ ) const;
	Bool_t CheckMCChannelDecay()const;
	
	Bool_t SetRecoCandidateParam(AliAODRecoDecayHF *recoCand);
	
 protected:
  
  
 private:	
	AliCFVertexingHF2Prong(const AliCFVertexingHF2Prong& c);
	AliCFVertexingHF2Prong& operator= (const AliCFVertexingHF2Prong& other);
	
	ClassDef(AliCFVertexingHF2Prong, 1);
  
};

#endif
