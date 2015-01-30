/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#ifndef ALIREDUCEDJETPARTICLE_H
#define ALIREDUCEDJETPARTICLE_H

#include <TObject.h>

class TLorentzVector;

namespace HighPtTracks {

class AliReducedJetParticle : public TObject {
public:
	AliReducedJetParticle();
	AliReducedJetParticle(double px, double py, double pz, double e, int pdgcode, bool rec);
	virtual ~AliReducedJetParticle() {}

	void FillLorentzVector(TLorentzVector &ref) const;
	double GetDistanceToJetMainAxis() const { return fDr; }
	int GetPdgCode() const { return fPdgCode; }
	bool IsReconstructed() const { return fIsReconstructed; }
	float GetDeltaPt() const { return fDeltaPt; }
	unsigned char GetNumberOfClustersTPC() const { return fNumberOfClustersTPC; }

	void SetKine(double px, double py, double pz, double e){
		fPx = px;
		fPy = py;
		fPz = pz;
		fE = e;
	}
	void SetDistanceToMainJetAxis(double dr) { fDr = dr; }
	void SetPdgCode(int pdg) { fPdgCode = pdg; }
	void SetReconstructed(bool rec) { fIsReconstructed = rec; }

	void SetDeltaPt(float deltaPt) { fDeltaPt = deltaPt; }
	void SetNumberOfClustersTPC(int ncls) { int fNumberOfClustersTPC = ncls; }

private:
	double 			fPx;
	double 			fPy;
	double			fPz;
	double			fE;
	double			fDr;
	int				fPdgCode;
	bool 			fIsReconstructed;
	float			fDeltaPt;
	unsigned char	fNumberOfClustersTPC;

	ClassDef(AliReducedJetParticle, 1);
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDJETPARTICLE_H */
