/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#ifndef ALIPARTICLEMAP_H
#define ALIPARTICLEMAP_H

#include <cstdlib>
#include <vector>
#include <map>

class AliVParticle;
class AliVTrack;

namespace HighPtTracks {

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

} /* namespace HighPtTracks */

#endif /* ALIPARTICLEMAP_H */
