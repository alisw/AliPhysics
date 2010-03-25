#ifndef ALIHLTMUONTRACK_H
#define ALIHLTMUONTRACK_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

///
/// @file   AliHLTMUONTrack.h
/// @author Indranil Das <indra.ehep@gmail.com> and Artur Szostak <artursz@iafrica.com>
/// @date   10 March 2010
/// @brief  Declaration of the track class used to store converted dHLT track data.
///

#include "TObject.h"
#include "TVector3.h"

class AliHLTMUONTriggerRecord;
class AliHLTMUONRecHit;

/**
 * AliHLTMUONTrack stores converted dHLT raw track data from HLT raw data blocks
 * as a ROOT object. This class is mainly for testing or as a helper object for
 * dHLT specific analysis, since it is sometimes easier to store and handle ROOT
 * objects than the raw data blocks.
 */
class AliHLTMUONTrack : public TObject
{
	/**
	 * Stream operator for usage with std::ostream classes.
	 * Allows usage such as:
	 *   AliHLTMUONTrack t; std::cout << t;
	 */
	friend std::ostream& operator << (
			std::ostream& stream,
			const AliHLTMUONTrack& track
		);

public:

	/**
	 * Constructor for creating a new track object.
	 * \note this class does not take ownership of the hit or trigger record
	 * objects and will not attempt to delete them. This must be done by the
	 * caller.
	 * @param id       The track ID number which must be unique for any event.
	 * @param sign     The particle's sign: -1, 1 or 0 if unknown.
	 * @param px       X component of the particle's momentum (GeV/c).
	 * @param py       Y component of the particle's momentum (GeV/c).
	 * @param pz       Z component of the particle's momentum (GeV/c).
	 * @param invmom   Inverse bending momentum (GeV/c).
	 * @param thetaX   The non-bending plane slope of the fitted track.
	 * @param thetaY   The bending plane slope of the fitted track.
	 * @param x        X coordinate of the particle's distance of closest
	 *                 approach (DCA) position (cm).
	 * @param y        Y coordinate of the particle's DCA position (cm).
	 * @param z        Z coordinate of the particle's DCA position (cm).
	 * @param chi2     The chi squared of the track fit.
	 * @param trigrec  Corresponding trigger record used as a seed to find
	 *                 this track.
	 * @param hits     The array of 16 hit coordinates found for the track.
	 *                 If NULL then then all hit coordinates are set to empty.
	 */
	AliHLTMUONTrack(
			Int_t id = -1, Int_t sign = 0,
			Float_t px = 0, Float_t py = 0, Float_t pz = 0,
			Float_t invmom = 0, Float_t thetaX = 0, Float_t thetaY = 0,
			Float_t x = 0, Float_t y = 0, Float_t z = 0,
			Float_t chi2 = -1,
			const AliHLTMUONTriggerRecord* trigrec = NULL,
			const AliHLTMUONRecHit* hits[16] = NULL
		);
	
	/**
	 * Default destructor.
	 */
	virtual ~AliHLTMUONTrack() {}

	/**
	 * Returns the track ID number, which should be unique for an event.
	 */
	Int_t Id() const { return fId; }
	
	/**
	 * Returns the sign of the particle: -1, 1 or 0 if the sign is unknown.
	 */
	Int_t Sign() const { return fSign; }
	
	/**
	 * Returns the inverse momentum in the bending plane, i.e. 1/p_yz, in units c/GeV.
	 */
	Float_t InverseBendingMomentum() const { return fInverseBendingMomentum; }
	
	/**
	 * Returns the slope of the fitted track in the non-bending plane.
	 */
	Float_t ThetaX() const { return fThetaX; }
	
	/**
	 * Returns the slope of the fitted track in the bending plane.
	 */
	Float_t ThetaY() const { return fThetaY; }

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
	 * Returns the distance of closest approach (DCA) position in centimetres.
	 */
	const TVector3& VertexDCA() const { return fVertexDCA; }

	/**
	 * Returns the non-bending plane coordinate for the distance of closest approach (DCA) in cm.
	 */
	Double_t X() const { return fVertexDCA.X(); }

	/**
	 * Returns the bending plane coordinate for the distance of closest approach (DCA) in cm.
	 */
	Double_t Y() const { return fVertexDCA.Y(); }

	/**
	 * Returns the z coordinate for the distance of closest approach (DCA) in cm.
	 */
	Double_t Z() const { return fVertexDCA.Z(); }

	/**
	 * Returns the chi squared of the track fit, indicating the quality of
	 * the fit.
	 */
	Float_t Chi2() const { return fChi2; }

	/**
	 * Returns the trigger record corresponding to this track.
	 * If NULL is returned then no trigger record was found or associated
	 * with the track.
	 */
	const AliHLTMUONTriggerRecord* TriggerRecord() const { return fTrigRec; }

	/**
	 * Returns the i'th hit set for this track.
	 * If NULL is returned then no hit was found or set for that hit position.
	 * @param i  Specifies the number of the hit to return.
	 *           Valid values are in the range [0..15].
	 */
	const AliHLTMUONRecHit* Hit(Int_t i) const;

	/**
	 * Returns the first hit found on the specified tracking chamber.
	 * If NULL is returned then no hit was found or set for that chamber.
	 * @param chamber  Specifies the chamber for which to return the hit.
	 *                 Valid values are in the range [1..14].
	 */
	const AliHLTMUONRecHit* HitByChamber(Int_t chamber) const;
	
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
	bool operator == (const AliHLTMUONTrack& track) const;

	bool operator != (const AliHLTMUONTrack& track) const
	{
		return not this->operator == (track);
	}

private:

	// Do not allow copying of this class.
	AliHLTMUONTrack(const AliHLTMUONTrack& track);
	AliHLTMUONTrack& operator = (const AliHLTMUONTrack& track);
	
	Int_t fId; ///< Track ID number which is unique for a particular event.
	Int_t fSign;  ///< The sign of the particle.
	Float_t fInverseBendingMomentum; ///< One over the momentum of the fitted track in GeV/c.
	Float_t fThetaX; ///< The slope of the fitted track in the non-bending plane.
	Float_t fThetaY; ///< The slope of the fitted track in the bending plane.
	TVector3 fMomentum; ///< Momentum vector of the particle in GeV/c.
	TVector3 fVertexDCA; ///< The position for the distance of closest approach to the vertex in cm.
	Float_t fChi2; ///< Chi squared of fit.
	const AliHLTMUONTriggerRecord* fTrigRec;  ///< Corresponding trigger record.
	const AliHLTMUONRecHit* fHit[16];   ///< Particle hits on tracking chambers.

	ClassDef(AliHLTMUONTrack, 1); // Track object containing data converted from a dHLT internal track structure.
};

#endif // ALIHLTMUONTRACK_H
