/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#ifndef ALIREDUCEDJETCONSTITUENT_H
#define ALIREDUCEDJETCONSTITUENT_H

#include <TObject.h>

class TLorentzVector;
class TParticlePDG;

namespace HighPtTracks {

class AliReducedJetConstituent : public TObject {
public:
	AliReducedJetConstituent();
	AliReducedJetConstituent(double px, double py, double pz, double e, int pdg);
	virtual ~AliReducedJetConstituent() {}

	void Set(double px, double py, double pz, double e){
		fPx = px;
		fPy = py;
		fPz = pz;
		fE = e;
	}
	void SetPdgCode(int pdg){ fPdgCode = pdg; }

	void FillLorentzVector(TLorentzVector &target) const;
	int GetPdgCode() const { return fPdgCode; }
	TParticlePDG *GetPDGParticle() const;

private:
	double 			fPx;
	double 			fPy;
	double 			fPz;
	double 			fE;
	int 			fPdgCode;

	ClassDef(AliReducedJetConstituent, 1);
};

} /* namespace HighPtTracks */

#endif /* PWGJE_EMCALJETTASKS_TRACKS_ALIREDUCEDJETCONSTITUENT_H_ */
