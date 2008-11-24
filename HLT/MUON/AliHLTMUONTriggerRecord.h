#ifndef ALIHLTMUONTRIGGERRECORD_H
#define ALIHLTMUONTRIGGERRECORD_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// @file   AliHLTMUONTriggerRecord.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 Sep 2007
/// @brief  Declaration of the trigger record structure in ROOT object format for dHLT.
///

#include "TObject.h"
#include "TVector3.h"
#include <ostream>

/**
 * Trigger record class containing information about a dimuon L0 trigger
 * local board decision in a ROOT object.
 * This class is mainly for testing or as a helper object for dHLT specific analysis,
 * since it is sometimes easier to store and handle ROOT objects.
 */
class AliHLTMUONTriggerRecord : public TObject
{
	/**
	 * Stream operator for usage with std::ostream classes.
	 * Allows usage such as:
	 *   AliHLTMUONTriggerRecord tr; std::cout << tr;
	 */
	friend std::ostream& operator << (
			std::ostream& stream,
			const AliHLTMUONTriggerRecord& trigrec
		);

public:

	/**
	 * Constructor for creating a new trigger record.
	 * @param id    The trigger record ID number unique for an event.
	 * @param sign  The particle's sign. Must be -1, 1 or 0 if the sign is unknown.
	 * @param px    X component of the particle's momentum.
	 * @param py    Y component of the particle's momentum.
	 * @param pz    Z component of the particle's momentum.
	 * @param sourceDDL  The DDL from which this trigger record originates.
	 * @param zf    The Z coordinate of the middle of the magnetic field assumed
	 *              during momentum calculation.
	 * @param qbl   The integrated magnetic field strength assumed during momentum
	 *              calculation.
	 */
	AliHLTMUONTriggerRecord(
			Int_t id = -1,
			Int_t sign = 0,
			Float_t px = 0,
			Float_t py = 0,
			Float_t pz = 0,
			Int_t sourceDDL = -1,
			Float_t zf = 0,
			Float_t qbl = 0
		);
	
	/**
	 * Default destructor.
	 */
	virtual ~AliHLTMUONTriggerRecord() {}

	/**
	 * Returns the trigger record ID number, which is unique for an event.
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
	 * Returns the hit coordinate on the specified chamber in the AliRoot
	 * coordinate system.
	 * @param chamber  The chamber for which to fetch the hit. Valid values
	 *                 are in the range [11..14].
	 */
	const TVector3& Hit(Int_t chamber) const;

	/**
	 * Returns the X coordinate of the reconstructed hit in centimetres.
	 * @param chamber  The chamber for which to fetch the X coordinate.
	 *                 Valid values are in the range [11..14].
	 */
	Double_t X(Int_t chamber) const { return Hit(chamber).X(); }

	/**
	 * Returns the Y coordinate of the reconstructed hit in centimetres.
	 * @param chamber  The chamber for which to fetch the Y coordinate.
	 *                 Valid values are in the range [11..14].
	 */
	Double_t Y(Int_t chamber) const { return Hit(chamber).Y(); }

	/**
	 * Returns the Z coordinate of the reconstructed hit in centimetres.
	 * @param chamber  The chamber for which to fetch the Z coordinate.
	 *                 Valid values are in the range [11..14].
	 */
	Double_t Z(Int_t chamber) const { return Hit(chamber).Z(); }
	
	/**
	 * Returns the source DDL from which this trigger record originates.
	 * -1 is returned if this was not set.
	 */
	Int_t SourceDDL() const { return fSourceDDL; }
	
	/**
	 * Returns the detector element ID number for the hit on the specified
	 * chamber. -1 is returned if this information was not set.
	 * @param chamber  The chamber for which to fetch the detector element ID.
	 *                 Valid values are in the range [11..14].
	 */
	Int_t DetElemId(Int_t chamber) const;
	
	/**
	 * Returns the 16 bit X pattern from the local board.
	 * -1 is returned if this information was not set.
	 * @param chamber  The chamber for which to fetch the bit pattern.
	 *                 Valid values are in the range [11..14].
	 */
	Int_t PatternX(Int_t chamber) const;
	
	/**
	 * Returns the 16 bit Y pattern from the local board.
	 * -1 is returned if this information was not set.
	 * @param chamber  The chamber for which to fetch the bit pattern.
	 *                Valid values are in the range [11..14].
	 */
	Int_t PatternY(Int_t chamber) const;
	
	/**
	 * Returns the Z coordinate in the middle of the magnetic field used to
	 * calculate the momentum.
	 */
	Float_t Zmiddle() const { return fZmiddle; }
	
	/**
	 * Returns the integrated magnetic field strength times charge used in
	 * the calculation of the momentum. Value returned in (T.m) tesla metres.
	 */
	Float_t QBL() const { return fQBL; }
	
	/**
	 * Sets the hit coordinate (in AliRoot global coordinates) on the
	 * given chamber.
	 * @param chamber  The chamber for which to set the hit. Valid values
	 *                 are in the range [11..14].
	 * @param x  The X coordinate of the hit in centimetres.
	 * @param y  The Y coordinate of the hit in centimetres.
	 * @param z  The Z coordinate of the hit in centimetres.
	 */
	void SetHit(Int_t chamber, Float_t x, Float_t y, Float_t z);
	
	/**
	 * Sets the debugging information for the hit on the specified chamber.
	 * @param chamber  The chamber for which to set the debugging information.
	 *                Valid values are in the range [11..14].
	 * @param detElemId  The detector element ID.
	 * @param patterX    The X bit pattern from the local board.
	 * @param patterY    The Y bit pattern from the local board.
	 */
	void SetHitDebugInfo(
			Int_t chamber,
			Int_t detElemId, UShort_t patternX, UShort_t patternY
		);
	
	/**
	 * Prints the details of the trigger record.
	 * @param option  A case sensitive string that can contain one of the
	 *     following strings:
	 *       "compact" - Prints just the momentum, sign and ID of the trigger
	 *                   record in a terse format.
	 *       "detail" - Prints also the hit information.
	 *       "all" - Prints all known information about this trigger record.
	 *     If the string contains an empty option or NULL then the default is
	 *     to print compactly.
	 */
	virtual void Print(Option_t* option = NULL) const;
	
	// Methods inherited from TObject
	virtual Bool_t IsSortable() const { return kTRUE; }
	Int_t Compare(const TObject* obj) const;

	// Implement comparison operators.
	bool operator == (const AliHLTMUONTriggerRecord& trigrec) const;

	bool operator != (const AliHLTMUONTriggerRecord& trigrec) const
	{
		return not this->operator == (trigrec);
	}

private:

	Int_t fId; ///< Each trigger record should have an ID number unique for a given event.
	Int_t fSign;  ///< The sign of the particle: -1 or 1. 0 indicates unknown value.
	TVector3 fMomentum; ///< Momentum vector of the particle in GeV/c.
	TVector3 fHit[4];   ///< hit coordinates on trigger chambers 11 to 14.
	
	// The following is debugging information and may not be filled if the
	// dHLT components were not set to produce this information.
	Int_t fSourceDDL;  ///< The DDL from which this trigger record originates.
	Int_t fDetElemId[4]; ///< The detector element ID for the hit on each chamber 11 to 14.
	Int_t fPatternX[4];  ///< The X pattern from the local board structure for chambers 11 to 14. -1 if invalid.
	Int_t fPatternY[4];  ///< The Y pattern from the local board structure for chambers 11 to 14. -1 if invalid.
	
	// Parameters used in momentum estimation:
	Float_t fZmiddle; ///< Particle momentum X component in GeV/c.
	Float_t fQBL;     ///< The integrated magnetic field times charge in (T.m) tesla metres.
		
	ClassDef(AliHLTMUONTriggerRecord, 3);  // Trigger record object translated from dHLT internal raw data.
};

#endif // ALIHLTMUONTRIGGERRECORD_H
