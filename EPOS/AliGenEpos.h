/*
 * AliGenEpos.cpp
 *
 *  ALICE event generator based on EPOS model from Klaus Werner
 *
 *  Created on: Feb 28, 2009
 *      Author: Piotr Ostrowski, postrow@if.pw.edu.pl
 */

#ifndef ALIGENEPOS_H_
#define ALIGENEPOS_H_

#include "AliGenMC.h"
#include "TEpos.h"

class AliGenEpos: public AliGenMC {
public:
	AliGenEpos();
	AliGenEpos(Int_t npart);
	virtual void Init();
	virtual void Generate();

	virtual ~AliGenEpos();

	void SetImpactParameterRange(Float_t bmin, Float_t bmax) { fBmin = bmin; fBmax = bmax; }
	void SetReactionPlaneAngleRange(Float_t phimin, Float_t phimax) { fPhiMin = phimin; fPhiMax = phimax; }
	void AddNoDecay(Int_t nodecay) { GetTEpos()->AddNoDecay(nodecay); }
	void AddExtraInputLine(const char *line) { GetTEpos()->AddExtraInputLine(line); }
	Float_t GetPhiMin() { return fPhiMin; }
	Float_t GetPhiMax() { return fPhiMax; }
	Float_t GetBmin() { return fBmin; }
	Float_t GetBMax() { return fBmax; }

	void FilterModelOutput(Bool_t value) {fFilterModelOutput = value;}
	Bool_t IsModelOutputFiltered() { return fFilterModelOutput; }
protected:
	virtual TEpos* GetTEpos() { return (TEpos *)fMCEvGen; }

	Float_t fBmin; //minimum impact parameter
	Float_t fBmax; //maximum impact parameter

	Float_t fPhiMin; // reaction plane angle minimum
	Float_t fPhiMax; // reaction plane angle maximum

	Bool_t fFilterModelOutput; //if true it will filter out internal model entities from the stack
private:

	ClassDef(AliGenEpos,1)
};

#endif /* ALIGENEPOS_H_ */
