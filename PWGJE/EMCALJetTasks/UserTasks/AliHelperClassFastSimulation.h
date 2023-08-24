// *************************************************************************
// * Class to do fast simulations
// *************************************************************************

#ifndef ALIHELPERCLASSFASTSIMULATION_H
#define ALIHELPERCLASSFASTSIMULATION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
class TF1;
class AliFJWrapper;
class AliAODMCParticle;
class TRandom3;
class TObjArray;
class TArrayI;
#include <fastjet/PseudoJet.hh>

class AliHelperClassFastSimulation {
public:
	AliHelperClassFastSimulation(TF1** effFunctions, Double_t simEffFactor, Double_t simResFactor, Double_t simRes, Double_t jetMinPt, Double_t jetMaxEta, Double_t jetMinEta, AliFJWrapper* wrapper = 0x0);
	virtual ~AliHelperClassFastSimulation();
	
	void 									AddParticle(AliAODMCParticle* part);
	void 									AddInputVector(Double_t pT, Double_t phi, Double_t eta, Double_t mass, Int_t index, Double_t charge, Int_t pdgCode);
	void 									Run();
	AliEmcalJet* 					GetJet(Int_t nOfJet);
  AliAODMCParticle* 		GetTrackOfJet(Int_t nOfTrack, Int_t nOfJet);
	AliAODMCParticle* 		GetParticleFromConstituent(fastjet::PseudoJet constituent);
	AliFJWrapper* 				GetStandardJetWrapper();
	
	UInt_t 								GetNJets();
	UInt_t                GetNParticlesOfJet(UInt_t jetNumber);
	
private:
  AliHelperClassFastSimulation(const AliHelperClassFastSimulation&)           ;  //Not implemented
  AliHelperClassFastSimulation &operator=(const AliHelperClassFastSimulation&);  //Not implemented
  
	Double_t 							fFastSimEffFactor;
	Double_t 							fFastSimResFactor;
	Double_t 							fFastSimRes;
	TF1** 	 							fEffFunctions;
	
	Double_t 							fJetMinPt;
	Double_t 							fJetMaxEta;
	Double_t 							fJetMinEta;
	
	TRandom3* 						fRandom;
	
	AliFJWrapper* 				fWrapper;
	
	TObjArray* 						fJetArray;
	TObjArray* 					  fParticlesArray;
	TArrayI*            	fPDGCodeArray;
	
	ClassDef(AliHelperClassFastSimulation, 1);
};

#endif
