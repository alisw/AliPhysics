/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#ifndef ALIHIGHPTRECONSTRUCTIONEFFICIENCY_H
#define ALIHIGHPTRECONSTRUCTIONEFFICIENCY_H

#include "AliAnalysisTaskSE.h"

#include <map>
#include <vector>
#include <fastjet/PseudoJet.hh>
#include "AliParticleMap.h"

class AliESDtrackCuts;
class AliVParticle;
class AliVTrack;
class TList;
class TNtuple;

//namespace fastjet{
//	class PseudoJet;
//};
namespace HighPtTracks {

class AliReducedJetEvent;
class AliReducedJetInfo;

class AliHighPtReconstructionEfficiency : public AliAnalysisTaskSE {
public:
	AliHighPtReconstructionEfficiency();
	AliHighPtReconstructionEfficiency(const char *name);
	virtual ~AliHighPtReconstructionEfficiency();

	virtual void UserCreateOutputObjects();
	virtual bool UserNotify();
	virtual void UserExec(Option_t* /*option*/);
	void Terminate(Option_t* /*option*/) {}

	void SetMaxEtaParticles(double maxeta) { fMaxEtaParticles = maxeta; }
	void SetMaxEtaJets(double maxeta) { fMaxEtaJets = maxeta; }
	void SetMinPtTracks(double minpt) { fMinPtParticles = minpt; }
	void SetMaxDR(double maxdr) { fMaxDR = maxdr; }

	void SetTaskDebugMode() { fTaskDebugMode = true; }

protected:
	bool IsSelected(const AliVTrack * const track) const;
	bool IsTrueSelected(const AliVParticle *const track) const;
	void SelectParticlesForJetfinding(TList &particles) const;
	void CreateRectrackLookup();
	std::vector<AliReconstructedParticlePair> SelectParticles() const;
	double GetDR(const fastjet::PseudoJet &recjet, const AliVParticle *const inputtrack) const;
	AliVTrack *FindReconstructedParticle(int label) const;
	AliVTrack *FindReconstructedParticleFast(int label) const;
	void ProcessJet(AliReducedJetInfo *const jet, const std::vector<AliReconstructedParticlePair> &particles) const;
	void ConvertConstituents(AliReducedJetInfo * const recjet, const fastjet::PseudoJet &inputjet);
	bool IsPhysicalPrimary(const AliVParticle *const part) const;
	bool PythiaInfoFromFile(const char* currFile, double &fXsec, double &fTrials, int &pthard) const ;

private:
	AliHighPtReconstructionEfficiency(const AliHighPtReconstructionEfficiency &);
	AliHighPtReconstructionEfficiency &operator=(const AliHighPtReconstructionEfficiency &);

	AliESDtrackCuts					*fTrackCuts;				//!
	AliParticleMap					*fParticleMap;				//!

	// Output objects
	TTree							*fJetTree;
	AliReducedJetEvent				*fJetEvent;

	double 							fMaxEtaJets;
	double							fMaxEtaParticles;
	double 							fMinPtParticles;
	double							fMaxDR;

	double							fCrossSection;
	double							fNtrials;
	int								fPtHardBin;

	bool 							fTaskDebugMode;

	ClassDef(AliHighPtReconstructionEfficiency, 1);
};

}
#endif /* ALIHIGHPTRECONSTRUCTIONEFFICIENCY_H */
