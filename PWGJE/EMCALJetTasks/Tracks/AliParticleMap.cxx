/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/*
 * Helper structure: Particle map with correlation generated particle -> reconstructed particle,
 * sorted according to the MC label
 *
 * Author:
 *   Markus Fasel <markus.fasel@cern.ch>
 */
#include <iostream>

#include <TMath.h>

#include "AliVParticle.h"
#include "AliVTrack.h"

#include "AliParticleMap.h"

namespace HighPtTracks {

AliParticleMap::~AliParticleMap() {
	// Clean up all particle lists
	for(std::map<int, AliParticleList *>::iterator it = fParticles.begin(); it != fParticles.end(); ++it){
		delete it->second;
	}
}

void AliParticleMap::AddParticle(AliVTrack *track){
	int label = TMath::Abs(track->GetLabel());
	std::map<int, AliParticleList *>::iterator it = fParticles.find(label);
	if(it == fParticles.end()){			// not existing
		AliParticleList *nextparticle = new AliParticleList;
		nextparticle->AddParticle(track);
		fParticles.insert(std::pair<int, AliParticleList *>(label, nextparticle));
	} else {
		AliParticleList *mylist = it->second;
		mylist->AddParticle(track);
	}
}

AliParticleList* AliParticleMap::GetParticles(int label) const {
	AliParticleList *result = NULL, *content = NULL;
	std::map<int, AliParticleList *>::const_iterator found = fParticles.find(label);
	if(found != fParticles.end()){
		result = found->second;
	}
	return result;
}

void AliParticleMap::Print() const {
	for(std::map<int, AliParticleList *>::const_iterator it = fParticles.begin(); it != fParticles.end(); ++it){
		std::cout << "Particle with label " << it->first << ", number of reconstructed assigned: " << it->second->GetNumberOfParticles() << std::endl;
	}
}

} /* namespace HighPtTracks */
