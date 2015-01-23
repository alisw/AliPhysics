/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#ifndef ALIHIGHPTRECONSTRUCTIONEFFICIENCY_H
#define ALIHIGHPTRECONSTRUCTIONEFFICIENCY_H

#include "AliAnalysisTaskSE.h"

#include <map>
#include <vector>
#include <fastjet/PseudoJet.hh>

class AliESDtrackCuts;
class AliVParticle;
class AliVTrack;
class TList;
class TNtuple;

//namespace fastjet{
//	class PseudoJet;
//};

class AliReconstructedParticlePair{
public:
	AliReconstructedParticlePair():
		fTrueParticle(NULL),
		fRecParticle(NULL)
	{}
	AliReconstructedParticlePair(const AliReconstructedParticlePair &ref):
		fTrueParticle(ref.fTrueParticle),
		fRecParticle(ref.fRecParticle)
	{}
	AliReconstructedParticlePair &operator=(const AliReconstructedParticlePair &ref){
		if(this != &ref){
			fTrueParticle = ref.fTrueParticle;
			fRecParticle = ref.fRecParticle;
		}
		return *this;
	}
	~AliReconstructedParticlePair(){}

	const AliVParticle *GetMCTrueParticle() const { return fTrueParticle; }
	const AliVTrack *GetRecTrack() const { return fRecParticle; }

	void SetMCTrueParticle(const AliVParticle *const part) { fTrueParticle = part; }
	void SetRecParticle(const AliVTrack *const track) { fRecParticle = track; }

private:
	const AliVParticle 			*fTrueParticle;
	const AliVTrack 			*fRecParticle;
};

class AliParticleList{
public:
	AliParticleList():
		fParticles()
	{}
	~AliParticleList(){}

	void AddParticle(AliVTrack *track) { fParticles.push_back(track); }
	AliVTrack *GetParticle(int itrack) const { return fParticles[itrack]; }
	int GetNumberOfParticles() const { return fParticles.size(); }

private:
	std::vector<AliVTrack *> 			fParticles;
};

class AliParticleMap{
public:
	AliParticleMap():
		fParticles()
	{}
	~AliParticleMap();

	void AddParticle(AliVTrack *track);
	AliParticleList *GetParticles(int label) const;
	int GetNumberOfParticles() const { return fParticles.size(); }

	void Print() const;

private:
	std::map<int, AliParticleList *> 				fParticles;
};

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
	std::vector<fastjet::PseudoJet> FindJets(const TList &inputparticles) const;
	std::vector<AliReconstructedParticlePair> SelectParticles() const;
	double GetDR(const fastjet::PseudoJet &recjet, const AliVParticle *const inputtrack) const;
	AliVTrack *FindReconstructedParticle(int label) const;
	AliVTrack *FindReconstructedParticleFast(int label) const;
	void ProcessJet(const fastjet::PseudoJet &recjet, const std::vector<AliReconstructedParticlePair> &particles);
	bool IsPhysicalPrimary(const AliVParticle *const part) const;
	bool PythiaInfoFromFile(const char* currFile, double &fXsec, double &fTrials, int &pthard) const ;

private:
	AliHighPtReconstructionEfficiency(const AliHighPtReconstructionEfficiency &);
	AliHighPtReconstructionEfficiency &operator=(const AliHighPtReconstructionEfficiency &);

	AliESDtrackCuts					*fTrackCuts;				//!
	TNtuple							*fResults;					//!
	AliParticleMap					*fParticleMap;				//!

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

#endif /* ALIHIGHPTRECONSTRUCTIONEFFICIENCY_H */
