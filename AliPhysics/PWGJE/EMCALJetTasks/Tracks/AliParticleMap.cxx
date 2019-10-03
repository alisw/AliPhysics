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
#include <iostream>

#include <TMath.h>

#include "AliVParticle.h"
#include "AliVTrack.h"

#include "AliParticleMap.h"

namespace HighPtTracks {

/**
 * Destructor. Clean up all particle lists.
 */
AliParticleMap::~AliParticleMap() {
	for(std::map<int, AliParticleList *>::iterator it = fParticles.begin(); it != fParticles.end(); ++it){
		delete it->second;
	}
}

/**
 * Add particle to the list. In case the same label is already existing in the list, the particle is just added to the list,
 * otherwise a new entry for the label is created.
 *
 * \param track particle to be added
 */
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

/**
 * Get list of all reconstructed particles associated with a given Monte-Carlo Label.
 *
 * \param label Label of the particle
 * \return List of reconstructed particles (NULL if not found)
 */
AliParticleList* AliParticleMap::GetParticles(int label) const {
	AliParticleList *result = NULL, *content = NULL;
	std::map<int, AliParticleList *>::const_iterator found = fParticles.find(label);
	if(found != fParticles.end()){
		result = found->second;
	}
	return result;
}

/**
 * Print status of the particle map.
 */
void AliParticleMap::Print() const {
	for(std::map<int, AliParticleList *>::const_iterator it = fParticles.begin(); it != fParticles.end(); ++it){
		std::cout << "Particle with label " << it->first << ", number of reconstructed assigned: " << it->second->GetNumberOfParticles() << std::endl;
	}
}

} /* namespace HighPtTracks */
