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
 * Class with event base information and the reduced jets
 *
 * Author:
 *   Markus Fasel <markus.fasel@cern.ch>
 */
#include <TObjArray.h>

#include "AliReducedJetInfo.h"
#include "AliReducedJetEvent.h"

ClassImp(HighPtTracks::AliReducedJetEvent)

namespace HighPtTracks {

AliReducedJetEvent::AliReducedJetEvent():
	TObject(),
	fCrossSection(0),
	fTrials(0),
	fPtHard(0),
	fReconstructedJets(NULL)
{
}

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

AliReducedJetEvent::~AliReducedJetEvent() {
	delete fReconstructedJets;
}

void AliReducedJetEvent::AddReconstructedJet(AliReducedJetInfo* jet) {
	if(fReconstructedJets) fReconstructedJets->Add(jet);
}

int AliReducedJetEvent::GetNumberOfJets() const {
	if(!fReconstructedJets) return 0;
	return fReconstructedJets->GetEntries();
}

AliReducedJetInfo* AliReducedJetEvent::GetReconstructedJet(int ijet) const {
	if(!fReconstructedJets || ijet >= fReconstructedJets->GetEntries()) return NULL;
	return dynamic_cast<AliReducedJetInfo *>(fReconstructedJets->At(ijet));
}

} /* namespace HighPtTracks */
