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
#include <TLorentzVector.h>

#include "AliReducedJetParticle.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedJetParticle)
/// \endcond

namespace HighPtTracks {

/**
 * Dummy (I/O) constructor, not to be used
 */
AliReducedJetParticle::AliReducedJetParticle() :
	TObject(),
	fPx(0),
	fPy(0),
	fPz(0),
	fE(0),
	fDr(0),
	fPdgCode(0),
	fNTPCTrackReferences(0),
	fMatchedTrackContainer(NULL)
{
  fMatchedTrackContainer = new TObjArray();
  fMatchedTrackContainer->SetOwner();
}

/**
 * Main constructor, initializes the particle with the moment vector, the particle energy and
 * the PDG code. Also creates a container for tracks matched to the particle via the Monte-Carlo
 * label, which takes ownership over the matched tracks.
 *
 * \param px x-component of the momentum vector
 * \param py y-component of the momentum vector
 * \param pz z-component of the momentum vector
 * \param e Particle energy
 * \param pdgcode PDG code of the particle
 */
AliReducedJetParticle::AliReducedJetParticle(double px, double py, double pz, double e, int pdgcode) :
	TObject(),
	fPx(px),
	fPy(py),
	fPz(pz),
	fE(e),
	fDr(0),
	fPdgCode(pdgcode),
	fNTPCTrackReferences(0),
  fMatchedTrackContainer(NULL)
{
  fMatchedTrackContainer = new TObjArray();
  fMatchedTrackContainer->SetOwner();
}

/**
 * Copy constructor. For each matched track, a copy is created and stored in the
 * new object.
 *
 * \param ref Reference for the copy
 */
AliReducedJetParticle::AliReducedJetParticle(const AliReducedJetParticle& ref):
	TObject(ref),
	fPx(ref.fPx),
	fPy(ref.fPy),
	fPz(ref.fPz),
	fE(ref.fE),
	fDr(ref.fDr),
	fPdgCode(ref.fPdgCode),
	fNTPCTrackReferences(ref.fNTPCTrackReferences),
  fMatchedTrackContainer(NULL)
{
  fMatchedTrackContainer = new TObjArray;
  fMatchedTrackContainer->SetOwner();

  for(TIter trackIter = TIter(ref.fMatchedTrackContainer).Begin(); trackIter != TIter::End(); ++trackIter){
    AliReducedMatchedTrack *othermatched = static_cast<AliReducedMatchedTrack *>(*trackIter);
    fMatchedTrackContainer->Add(new AliReducedMatchedTrack(*othermatched));
  }
}

/**
 * Assignment operator. For each matched track, a copy is created and stored in the
 * new object.
 *
 * \param ref Reference for the copy
 * \return This object
 */
AliReducedJetParticle& AliReducedJetParticle::operator=(const AliReducedJetParticle& ref) {
  if(this != &ref){
    delete fMatchedTrackContainer;

    fMatchedTrackContainer = new TObjArray;
    fMatchedTrackContainer->SetOwner();

    for(TIter trackIter = TIter(ref.fMatchedTrackContainer).Begin(); trackIter != TIter::End(); ++trackIter){
      AliReducedMatchedTrack *othermatched = static_cast<AliReducedMatchedTrack *>(*trackIter);
      fMatchedTrackContainer->Add(new AliReducedMatchedTrack(*othermatched));
    }

  }
  return *this;
}

/**
 * Destructor. Deletes all matched tracks associated to the particle.
 */
AliReducedJetParticle::~AliReducedJetParticle() {
  delete fMatchedTrackContainer;
}

/**
 * Access to the particle kinematics information via a TLorentzVector which is filled
 * inside the function
 *
 * \param ref The TLorentzVector to be filled
 */
void AliReducedJetParticle::FillLorentzVector(TLorentzVector& ref) const {
	ref.SetPxPyPzE(fPx, fPy, fPz, fE);
}

/**
 * Adds a new matched reconstructed track to the jet particle. Tracks are matched via the
 * Monte-Carlo label.
 *
 * \param trk Reconstructed track matched to the particle
 */
void AliReducedJetParticle::AddMatchedTrack(AliReducedMatchedTrack *trk) {
  fMatchedTrackContainer->Add(trk);
}

/**
 * Calculate the difference between reconstructed and generated \f$ p_{t} \f$ relative to the
 * generated \f$ p_{t} \f$ . If no matching track is found for the index, returns -1000.
 *
 * \param itrk Number of the matched track to be checked (0 as default)
 * \return \f$(p_{t,gen} - p_{t,rec})/p_{t,rec}\f$ (-1000 if no track is available)
 */
double AliReducedJetParticle::GetDeltaPt(int itrk) const {
	if(itrk >=  fMatchedTrackContainer->GetEntries()) return -1000.;
	TLorentzVector partvec;
	partvec.SetPxPyPzE(fPx, fPy, fPz, fE);
	double particlePt = partvec.Pt();
	AliReducedMatchedTrack *mytrk = static_cast<AliReducedMatchedTrack *>(fMatchedTrackContainer->At(itrk));
	return (mytrk->Pt()-particlePt)/particlePt;
}

} /* namespace HighPtTracks */
