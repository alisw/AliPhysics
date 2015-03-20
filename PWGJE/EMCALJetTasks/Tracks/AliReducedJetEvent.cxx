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
#include <TObjArray.h>

#include "AliReducedJetInfo.h"
#include "AliReducedJetEvent.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedJetEvent)
/// \endcond

namespace HighPtTracks {

/**
 * Dummy (I/O) constructor, not to be uesd
 */
AliReducedJetEvent::AliReducedJetEvent():
	TObject(),
	fCrossSection(0),
	fTrials(0),
	fPtHard(0),
	fReconstructedJets(NULL)
{
}

/**
 * Main consturctor, initialising also PYTHIA information (cross section, number of trials, \f$ p_{t} \f$
 * hard)
 *
 * \param crosssection Cross section
 * \param ntrials Number of trials
 * \param pthard \f$ p_{t} \f$ of the hard interaction
 */
AliReducedJetEvent::AliReducedJetEvent(double crosssection, int ntrials, double pthard) :
	TObject(),
	fCrossSection(crosssection),
	fTrials(ntrials),
	fPtHard(pthard),
	fReconstructedJets(NULL)
{
	fReconstructedJets = new TObjArray;
	fReconstructedJets->SetOwner(true);
}

/**
 * Copy constructor, creating a new jet event copying information from a reference event
 * (reconstructed jets). Creates a deep copy - this event will take ownership over its objects.
 *
 * \param ref Reference jet event for the copy
 */
AliReducedJetEvent::AliReducedJetEvent(const AliReducedJetEvent& ref):
	TObject(ref),
	fCrossSection(ref.fCrossSection),
	fTrials(ref.fTrials),
	fPtHard(ref.fPtHard),
	fReconstructedJets(NULL)
{
	fReconstructedJets = new TObjArray;
	fReconstructedJets->SetOwner(true);

	AliReducedJetInfo *recjet(NULL);
	TIter jetiter(ref.fReconstructedJets);
	while((recjet = dynamic_cast<AliReducedJetInfo *>(jetiter()))){
		fReconstructedJets->Add(new AliReducedJetInfo(*recjet));
	}
}

/**
 * Assignment operator, copying information from a reference event (reconstructed jets) into
 * this event. Creates a deep copy - this event will take ownership over its objects.
 *
 * \param ref Reference jet event for the copy
 * \return This object
 */
AliReducedJetEvent& AliReducedJetEvent::operator=(const AliReducedJetEvent& ref) {
	TObject::operator=(ref);
	if(&ref != this){
		fCrossSection = ref.fCrossSection;
		fTrials = ref.fTrials;
		fPtHard = ref.fPtHard;

		fReconstructedJets->Clear();
		AliReducedJetInfo *recjet(NULL);
		TIter jetiter(ref.fReconstructedJets);
		while((recjet = dynamic_cast<AliReducedJetInfo *>(jetiter()))){
			fReconstructedJets->Add(new AliReducedJetInfo(*recjet));
		}
	}
	return *this;
}

/**
 * Destructor, deleting all associated jets
 */
AliReducedJetEvent::~AliReducedJetEvent() {
	delete fReconstructedJets;
}

/**
 * Adds a new reduced jet to this event. The jet event will take ownership over the reduced jet.
 *
 * \param jet
 */
void AliReducedJetEvent::AddReconstructedJet(AliReducedJetInfo* jet) {
	if(fReconstructedJets) fReconstructedJets->Add(jet);
}

/**
 * Get the number of jets stored in this jet event.
 *
 * \return The number of jets stored in this jet event
 */
int AliReducedJetEvent::GetNumberOfJets() const {
	if(!fReconstructedJets) return 0;
	return fReconstructedJets->GetEntries();
}

/**
 * Get jet stored in this jet event via its position in the container.
 *
 * \param ijet Index of the jet in the event.
 * \return The reduced jet (NULL if index is out of range)
 */
AliReducedJetInfo* AliReducedJetEvent::GetReconstructedJet(int ijet) const {
	if(!fReconstructedJets || ijet >= fReconstructedJets->GetEntries()) return NULL;
	return dynamic_cast<AliReducedJetInfo *>(fReconstructedJets->At(ijet));
}

} /* namespace HighPtTracks */
