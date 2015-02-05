/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#ifndef ALIREDUCEDJETINFO_H
#define ALIREDUCEDJETINFO_H

#include <TObject.h>

class TLorentzVector;
class TObjArray;

namespace HighPtTracks {

class AliReducedJetParticle;
class AliReducedJetConstituent;

class AliReducedJetInfo  : public TObject {
public:
	AliReducedJetInfo();
	AliReducedJetInfo(double px, double py, double pz, double e);
	AliReducedJetInfo(const AliReducedJetInfo &ref);
	AliReducedJetInfo &operator=(const AliReducedJetInfo &ref);
	virtual ~AliReducedJetInfo();

	void Set(double px, double py, double pz, double e){
		fPx = px;
		fPy = py;
		fPz = pz;
		fE = e;
	}
	void AddParticleInCone(AliReducedJetParticle *part);
	void AddConstituent(AliReducedJetConstituent *con);

	void GetPxPyPxE(double &px, double &py, double &pz, double e){
		px = fPx;
		py = fPy;
		pz = fPz;
		e = fE;
	}
	void FillLorentzVector(TLorentzVector &vec ) const;
	int GetNumberOfMatchedParticles() const;
	TObjArray *GetListOfMatchedParticles() const { return fParticlesInCone; }
	TObjArray *GetListOfConstituents() const { return fConstituents; }
	AliReducedJetParticle *GetMatchedParticle(int ipart) const;

private:
	double					fPx;
	double 					fPy;
	double					fPz;
	double					fE;
	TObjArray				*fConstituents;
	TObjArray				*fParticlesInCone;

	ClassDef(AliReducedJetInfo, 1);
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDJETINFO_H */
