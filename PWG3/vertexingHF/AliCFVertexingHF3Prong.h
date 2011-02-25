#ifndef ALICFVERTEXINGHF3PRONG_H
#define ALICFVERTEXINGHF3PRONG_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to compute variables for correction framework           //  
// for 3-body decays of D mesons (D+, Ds, Lc)                    //
// in bins of cut variables                                      //
// Origin:       Francesco Prino (prino@to.infn.it)              //
//               Renu Bala       (bala@to.infn.it)               //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliCFVertexingHF.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF3Prong.h"

class AliAODMCParticle;
class TClonesArray;
class AliCFVertexingHF;
class AliESDtrack;
class TDatabasePDG;

class AliCFVertexingHF3Prong : public AliCFVertexingHF{
	public:
		
  AliCFVertexingHF3Prong(Int_t decay);
  AliCFVertexingHF3Prong(TClonesArray *mcArray, UShort_t originDselection, Int_t decay);  
  virtual ~AliCFVertexingHF3Prong(){};
  
  
  
  Bool_t GetGeneratedValuesFromMCParticle(Double_t* /*vectorMC*/);
  Bool_t GetRecoValuesFromCandidate(Double_t* /*vectorReco*/ ) const;
  Bool_t CheckMCChannelDecay()const;
  
  Bool_t SetRecoCandidateParam(AliAODRecoDecayHF *recoCand);
  virtual void SetDecay3Prong(Int_t decay){fDecay=decay;}
 
  
 protected:
  
  
  
 private:	
  AliCFVertexingHF3Prong(const AliCFVertexingHF3Prong& c);
  AliCFVertexingHF3Prong& operator= (const AliCFVertexingHF3Prong& other);
  
  Int_t fDecay;   // decay mode id

  ClassDef(AliCFVertexingHF3Prong, 1);
  
};

#endif
