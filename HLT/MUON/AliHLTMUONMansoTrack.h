#ifndef ALIHLTMUONMANSOTRACK_H
#define ALIHLTMUONMANSOTRACK_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// @file   AliHLTMUONMansoTrack.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 Sep 2007
/// @brief  Declaration of the Manso track class used to store converted track data.
///

#include "TObject.h"
#include "TVector3.h"

class AliHLTMUONTriggerRecord;
class AliHLTMUONRecHit;

/**
 * AliHLTMUONMansoTrack stores converted dHLT raw track data as a ROOT object.
 * This class is mainly for testing or as a helper object for dHLT specific analysis,
 * since it is sometimes easier to store and handle ROOT objects.
 */
class AliHLTMUONMansoTrack : public TObject
{
	/**
	 * Stream operator for usage with std::ostream classes.
	 * Allows usage such as:
	 *   AliHLTMUONMansoTrack t; std::cout << t;
	 */
	friend std::ostream& operator << (
			std::ostream& stream,
			const AliHLTMUONMansoTrack& track
		);

public:

	/**
	 * Constructor for creating a track object with none, some or all 4 hits
	 * specified. Note: this class does not take ownership of the hit or trigger
	 * record objects and will not attempt to delete them.
	 * @param id       The track ID number which must be unique for any event.
	 * @param sign     The particle's sign: -1, 1 or 0 if unknown.
	 * @param px       X component of the particle's momentum (GeV/c).
	 * @param py       Y component of the particle's momentum (GeV/c).
	 * @param pz       Z component of the particle's momentum (GeV/c).
	 * @param chi2     The chi squared of the track fit.
	 * @param trigrec  Corresponding trigger record used as a seed to find
	 *                 this track.
	 * @param hit7     Hit on chamber 7, tracking station 4.
	 * @param hit8     Hit on chamber 8, tracking station 4.
	 * @param hit9     Hit on chamber 9, tracking station 5.
	 * @param hit10    Hit on chamber 10, tracking station 5.
	 * @param zf    The Z coordinate of the middle of the magnetic field assumed
	 *              during momentum calculation.
	 * @param qbl   The integrated magnetic field strength assumed during momentum
	 *              calculation.
	 */
	AliHLTMUONMansoTrack(
			Int_t id = -1, Int_t sign = 0,
			Float_t px = 0, Float_t py = 0, Float_t pz = 0,
			Float_t chi2 = -1,
			const AliHLTMUONTriggerRecord* trigrec = NULL,
			const AliHLTMUONRecHit* hit7 = NULL,
			const AliHLTMUONRecHit* hit8 = NULL,
			const AliHLTMUONRecHit* hit9 = NULL,
			const AliHLTMUONRecHit* hit10 = NULL,
			Float_t zf = 0, Float_t qbl = 0
		);
	
	/**
	 * Default destructor.
	 */
	virtual ~AliHLTMUONMansoTrack() {}

	/**
	 * Returns the track ID number, which is unique for an event.
	 */
	Int_t Id() const { return fId; }
	
	/**
	 * Returns the sign of the particle: -1, 1 or 0 if the sign is unknown.
	 */
	Int_t Sign() const { return fSign; }

	/**
	 * Returns the momentum vector with components in GeV/c.
	 */
	const TVector3& Momentum() const { return fMomentum; }

	/**
	 * Returns the X component of the particle's momentum in GeV/c.
	 */
	Double_t Px() const { return fMomentum.Px(); }

	/**
	 * Returns the Y component of the particle's momentum in GeV/c.
	 */
	Double_t Py() const { return fMomentum.Py(); }

	/**
	 * Returns the Z component of the particle's momentum in GeV/c.
	 */
	Double_t Pz() const { return fMomentum.Pz(); }

	/**
	 * Returns the momentum magnitude of the particle in GeV/c.
	 */
	Double_t P() const { return fMomentum.Mag(); }

	/**
	 * Returns the transverse momentum of the particle in GeV/c.
	 */
	Double_t Pt() const { return fMomentum.Pt(); }

	/**
	 * Returns the polar angle of the momentum vector in radians.
	 */
	Double_t Polar() const { return fMomentum.Theta(); }

	/**
	 * Returns the azimuthal angle of the transverse momentum in radians.
	 */
	Double_t Phi() const { return fMomentum.Phi(); }

	/**
	 * Returns the chi squared of the track fit, indicating the quality of
	 * the fit.
	 */
	Float_t Chi2() const { return fChi2; }

	/**
	 * Returns the trigger record corresponding to this track.
	 * If NULL is returned then no trigger record was found.
	 */
	const AliHLTMUONTriggerRecord* TriggerRecord() const { return fTrigRec; }

	/**
	 * Returns tje hit found on the specified tracking chamber.
	 * If NULL is returned then no hit was found or set.
	 * @param chamber  Specifies the chamber for which to return the hit.
	 *                 Valid values are in the range [7..10].
	 */
	const AliHLTMUONRecHit* Hit(Int_t chamber) const;
	
	/**
	 * Prints the details of the track.
	 * @param option  A case sensitive string that can contain one of the
	 *     following strings:
	 *       "compact" - Prints just the momentum, sign and ID of the track
	 *                   in a terse format.
	 *       "detail" - Prints also the hit information.
	 *       "all" - Prints all known information about this track.
	 *     If the string contains an empty option or NULL then the default is
	 *     to print compactly.
	 */
	virtual void Print(Option_t* option = NULL) const;
	
	// Methods inherited from TObject
	virtual Bool_t IsSortable() const { return kTRUE; }
	Int_t Compare(const TObject* obj) const;

	// Implement comparison operators.
	bool operator == (const AliHLTMUONMansoTrack& track) const;

	bool operator != (const AliHLTMUONMansoTrack& track) const
	{
		return not this->operator == (track);
	}

private:

	// Do not allow copying of this class.
	AliHLTMUONMansoTrack(const AliHLTMUONMansoTrack& track);
	AliHLTMUONMansoTrack& operator = (const AliHLTMUONMansoTrack& track);
	
	Int_t fId; ///< Track ID number which is unique for a particular event.
	Int_t fSign;  ///< The sign of the particle.
	TVector3 fMomentum; ///< Momentum vector of the particle in GeV/c.
	Float_t fChi2; ///< Chi squared of fit.
	const AliHLTMUONTriggerRecord* fTrigRec;  ///< Corresponding trigger record.
	const AliHLTMUONRecHit* fHit[4];   ///< Particle hits on tracking chambers 7 to 10.
	
	// The following is debugging information and may not be filled if the
	// dHLT components were not set to produce this information.
	
	// Parameters used in momentum estimation:
	Float_t fZmiddle; ///< Particle momentum X component in GeV/c.
	Float_t fQBL;     ///< The integrated magnetic field times charge in (T.m) tesla metres.

	ClassDef(AliHLTMUONMansoTrack, 3); // Manso track object containing data converted from a dHLT internal track structure.
};

#endif // ALIHLTMUONMANSOTRACK_H
