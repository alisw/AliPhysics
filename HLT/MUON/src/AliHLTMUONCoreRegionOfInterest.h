#ifndef ALIHLTMUONCOREREGIONOFINTEREST_H
#define ALIHLTMUONCOREREGIONOFINTEREST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
// The region of interest object is used to encode/decode and work with boundary
// box type regions of interest. The 32 bit ROI codes are used to communicate
// regions of interest between different parts of the dHLT system. This is more
// efficient than sending 20 byte long region of interest objects.
//
////////////////////////////////////////////////////////////////////////////////

#include "AliHLTMUONBasicTypes.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONCoreCluster.h"


typedef UInt AliHLTMUONCoreROI;

enum
{
	kInvalidROI = 0xFFFFFFFF
};


/* Identification numbers specifying the tracking chambers of the dimuon
   spectrometer.
 */
enum AliHLTMUONCoreChamberID
{
	kUnknownChamber = -1,
	kChamber1 = 0,
	kChamber2 = 1,
	kChamber3 = 2,
	kChamber4 = 3,
	kChamber5 = 4,
	kChamber6 = 5,
	kChamber7 = 6,
	kChamber8 = 7,
	kChamber9 = 8,
	kChamber10 = 9
};


class AliHLTMUONCoreRegionOfInterest
{
public:

	AliHLTMUONCoreRegionOfInterest()
		: fChamber(kChamber1), fLeft(0.0), fRight(0.0), fBottom(0.0), fTop(0.0)
	{
		fChamber = kChamber1;
		fLeft = fRight = fBottom = fTop = 0.0;
	}
	
	/* This constructor decodes the ROI bit pattern into a region of
	   interest object.
	 */
	AliHLTMUONCoreRegionOfInterest(const AliHLTMUONCoreROI& code)
		: fChamber(kChamber1), fLeft(0.0), fRight(0.0), fBottom(0.0), fTop(0.0)
	{
		Decode(code);
	}

	/* Creates a region of interest around the given point for the
	   specified chamber.
	 */
	AliHLTMUONCoreRegionOfInterest(
			const AliHLTMUONCoreClusterPoint& point0,
			AliHLTMUONCoreChamberID chamber0
		)
		: fChamber(kChamber1), fLeft(0.0), fRight(0.0), fBottom(0.0), fTop(0.0)
	{
		CreateToContain(point0, chamber0);
	}

	/* Creates a region of interest around all the given points and for the
	   specified chamber.
	 */
	AliHLTMUONCoreRegionOfInterest(
			const AliHLTMUONCoreClusterPoint* points0, UInt count0,
			AliHLTMUONCoreChamberID chamber0
		)
		: fChamber(kChamber1), fLeft(0.0), fRight(0.0), fBottom(0.0), fTop(0.0)
	{
		CreateToContain(points0, count0, chamber0);
	}

	/* Creates a region of interest with the specified boundaries and for
	   the specified chamber.
	 */
	AliHLTMUONCoreRegionOfInterest(
			Float left, Float right, Float bottom, Float top,
			AliHLTMUONCoreChamberID chamber
		);

	/* Checks if the point is contained in this region of interest.
	 */
	bool Contains(const AliHLTMUONCoreClusterPoint& point) const;

	/* Checks if the point is contained in this region of interest and the
	   chamber number corresponds to this region object.
	 */
	bool Contains(
			const AliHLTMUONCoreClusterPoint& point,
			AliHLTMUONCoreChamberID chamber
		) const;

	/* Checks if the specified region of interest is contained in this
	   region of interest object.
	 */
	bool Contains(const AliHLTMUONCoreRegionOfInterest& roi) const;

	/* Creates a region of interest around the given point for the
	   specified chamber.
	 */
	void CreateToContain(
			const AliHLTMUONCoreClusterPoint& point,
			AliHLTMUONCoreChamberID chamber
		);

	/* Extends the region of interest to contain the specified point.
	 */
	void ExpandToContain(const AliHLTMUONCoreClusterPoint& point);

	/* Creates a region of interest around all the given points and for the
	   specified chamber.
	 */
	void CreateToContain(
			const AliHLTMUONCoreClusterPoint* points, UInt count,
			AliHLTMUONCoreChamberID chamber
		);

	/* Checks if the region of interest is within the boundaries imposed on
	   the specific chamber plane. This boundary is aproximately the square
	   box around the chamber's detection region.
	 */
	bool InBounds();

	/* Encodes the region of interest into a 32 bit code.
	 */
	AliHLTMUONCoreROI Encode() const;
	
	/* Encodes the region of interest into a 32 bit code, and returns the
	   hierarchal level the region was encoded at and the left and right
	   grid coordinate of the bottom left corner of the region boundary box.
	 */
	AliHLTMUONCoreROI Encode(UChar& level, UInt& l, UInt& b) const;

	/* Decodes a 32 bit region of interest code into this region of interest
	   object.
	 */
	void Decode(AliHLTMUONCoreROI code);

	/* Decodes the chamber number of the region of interest 32 bit code.
	 */
	static AliHLTMUONCoreChamberID DecodeChamber(AliHLTMUONCoreROI code);

	/* Decodes the 32 bit region of interest code into the chamber number,
	   hierarchal level, left and right grid coordinates of the region
	   boundary box. 
	 */
	static void Decode(
			AliHLTMUONCoreROI code, AliHLTMUONCoreChamberID& chamber,
			UChar& level, UInt& l, UInt& b
		);

	/* Returns the chamber number of the region of interest.
	 */
	AliHLTMUONCoreChamberID Chamber() const { return fChamber; };
	
	/* Returns the left hand boundary of the region of interest.
	 */
	Float Left() const    { return fLeft; };
	
	/* Returns the right hand boundary of the region of interest.
	 */
	Float Right() const   { return fRight; };
	
	/* Returns the bottom boundary of the region of interest.
	 */
	Float Bottom() const  { return fBottom; };
	
	/* Returns the top boundary of the region of interest.
	 */
	Float Top() const     { return fTop; };


	/* Typecast operator for implicit typecasing to 32 bit ROI codes.
	 */
	operator AliHLTMUONCoreROI () const { return Encode(); };


private:

	/* Converts the internal region of interest boundary box, which is
	   specified in as floats, into a regular integer grid.
	   l = left boundary, r = right boundary, b = bottom boundary, 
	   t = top boundary.
	 */
	inline void ConvertToGrid(UInt& l, UInt& r, UInt& b, UInt& t) const;
	
	/* Performs the inverse conversion of the method ConvertToGrid.
	   That is converts from a regular integer grid back to the internal
	   floating point boundary box specification.
	 */
	inline void ConvertBackFromGrid(register UInt l, register UInt r, register UInt b, register UInt t);

	/* Internal method for decoding 32 bit region codes. This method is
	   called by the Decode methods.
	 */
	static void DecodeBits(
			AliHLTMUONCoreROI code, AliHLTMUONCoreChamberID& chamber,
			UChar& colevel, UInt& l, UInt& b
		);

	enum {kNumberOfTrackingChambers = 10};  // Number of tracking chambers.

	// Boundary box scale numbers for each chamber. These are the boundary
	// boxes around the chambers detection surface.
	static Float fgPlaneScale[kNumberOfTrackingChambers];  // scale numbers.

	static UInt fgIndexOffsets[13];  // Offset numbers used in the encoding and decoding process.


	AliHLTMUONCoreChamberID fChamber; // Specifies the chamber the region of interest is on.
	Float fLeft;        // Left boundary of boundary box.
	Float fRight;       // Right boundary of boundary box.
	Float fBottom;      // Bottom boundary of boundary box.
	Float fTop;         // Top boundary of boundary box.
};


//-----------------------------------------------------------------------------
// Inline methods:


inline AliHLTMUONCoreRegionOfInterest::AliHLTMUONCoreRegionOfInterest(
		Float left, Float right, Float bottom, Float top,
		AliHLTMUONCoreChamberID chamber
	)
	: fChamber(kChamber1), fLeft(0.0), fRight(0.0), fBottom(0.0), fTop(0.0)
{
// Creates a region of interest with the specified boundaries and for
// the specified chamber.

	Assert( 0 <= chamber && chamber < (AliHLTMUONCoreChamberID)kNumberOfTrackingChambers );
	fChamber = chamber;
	fLeft = left;
	fRight = right;
	fBottom = bottom;
	fTop = top;
}


inline bool AliHLTMUONCoreRegionOfInterest::Contains(const AliHLTMUONCoreClusterPoint& point) const
{
// Checks if the point is contained in this region of interest.

	return fLeft <= point.X()
	  && point.X() <= fRight
	  && fBottom <= point.Y()
	  && point.Y() <= fTop;
}


inline bool AliHLTMUONCoreRegionOfInterest::Contains(
		const AliHLTMUONCoreClusterPoint& point,
		AliHLTMUONCoreChamberID chamber
	) const
{
// Checks if the point is contained in this region of interest and the
// chamber number corresponds to this region object.

	return fLeft <= point.X()
	  && point.X() <= fRight
	  && fBottom <= point.Y()
	  && point.Y() <= fTop
	  && fChamber == chamber;
}


inline bool AliHLTMUONCoreRegionOfInterest::Contains(
		const AliHLTMUONCoreRegionOfInterest& roi
	) const
{
// Checks if the specified region of interest is contained in this
// region of interest object.

	return fChamber == roi.fChamber
	  && fLeft <= roi.fLeft
	  && fRight >= roi.fRight
	  && fBottom <= roi.fBottom
	  && fTop >= roi.fTop;
}


#endif // ALIHLTMUONCOREREGIONOFINTEREST_H

