/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#ifndef ALIREDUCEDJETEVENT_H
#define ALIREDUCEDJETEVENT_H

#include <TObject.h>

class TObjArray;

namespace HighPtTracks {

class AliReducedJetInfo;

class AliReducedJetEvent : public TObject {
public:
	AliReducedJetEvent();
	AliReducedJetEvent(double crosssection, int ntrials, double pthard);
	AliReducedJetEvent(const AliReducedJetEvent &ref);
	AliReducedJetEvent &operator=(const AliReducedJetEvent &);
	virtual ~AliReducedJetEvent();

	void SetPythiaHardInfo(double crosssection, int ntrials, double pthard){
		fCrossSection = crosssection;
		fTrials = ntrials;
		fPtHard = pthard;
	}
	void AddReconstructedJet(AliReducedJetInfo *jet);

	double GetCrossSection() const { return fCrossSection; }
	int GetNumberOfTrials() const { return fTrials; }
	double GetPtHard() const { return fPtHard; }
	int GetNumberOfJets() const;
	AliReducedJetInfo *GetReconstructedJet(int ijet) const;
	TObjArray *GetListOfJets() const { return fReconstructedJets; }

private:
	double					fCrossSection;
	int 					fTrials;
	double					fPtHard;
	TObjArray				*fReconstructedJets;

	ClassDef(AliReducedJetEvent, 1);
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDJETEVENT_H */
