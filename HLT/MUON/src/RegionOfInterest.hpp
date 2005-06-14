////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_REGION_OF_INTEREST_HPP
#define dHLT_REGION_OF_INTEREST_HPP

#include "BasicTypes.hpp"
#include "Utils.hpp"
#include "Cluster.hpp"

namespace dHLT
{


const Int NUMBER_OF_TRACKING_CHAMBERS = 10;

typedef UInt ROI;

enum
{
	INVALID_ROI = 0xFFFFFFFF
};


/* Identification numbers specifying the tracking chambers of the dimuon
   spectrometer.
 */
enum ChamberID
{
	UnknownChamber = -1,
	Chamber1 = 0,
	Chamber2 = 1,
	Chamber3 = 2,
	Chamber4 = 3,
	Chamber5 = 4,
	Chamber6 = 5,
	Chamber7 = 6,
	Chamber8 = 7,
	Chamber9 = 8,
	Chamber10 = 9
};


/* The region of interest object is used to encode/decode and work with boundary
   box type regions of interest. The 32 bit ROI codes are used to communicate
   regions of interest between different parts of the dHLT system. This is more
   efficient than sending 20 byte long region of interest objects.
 */
class RegionOfInterest
{
public:

	RegionOfInterest()
	{
		chamber = Chamber1;
		left = right = top = bottom = 0.0;
	};
	
	/* This constructor decodes the ROI bit pattern into a region of
	   interest object.
	 */
	RegionOfInterest(const ROI& code)
	{
		Decode(code);
	};

	/* Creates a region of interest around the given point for the
	   specified chamber.
	 */
	RegionOfInterest(const ClusterPoint& point, const ChamberID chamber)
	{
		CreateToContain(point, chamber);
	};

	/* Creates a region of interest around all the given points and for the
	   specified chamber.
	 */
	RegionOfInterest(const ClusterPoint* points, const UInt count, const ChamberID chamber)
	{
		CreateToContain(points, count, chamber);
	};

	/* Creates a region of interest with the specified boundaries and for
	   the specified chamber.
	 */
	RegionOfInterest(const Float left, const Float right, const Float bottom, const Float top, const ChamberID chamber)
	{
		Assert( 0 <= chamber and chamber < NUMBER_OF_TRACKING_CHAMBERS );
		this->chamber = chamber;
		this->left = left;
		this->right = right;
		this->bottom = bottom;
		this->top = top;
	};


	/* Checks if the point is contained in this region of interest.
	 */
	bool Contains(const ClusterPoint& point) const
	{
		return left <= point.x and point.x <= right and
			bottom <= point.y and point.y <= top;
	};


	/* Checks if the point is contained in this region of interest and the
	   chamber number corresponds to this region object.
	 */
	bool Contains(const ClusterPoint& point, const ChamberID chamber) const
	{
		return left <= point.x and point.x <= right and
			bottom <= point.y and point.y <= top and
			this->chamber == chamber;
	};


	/* Checks if the specified region of interest is contained in this
	   region of interest object.
	 */
	bool Contains(const RegionOfInterest& roi) const
	{
		return chamber == roi.chamber and 
			left <= roi.left and right >= roi.right and
			bottom <= roi.bottom and top >= roi.top;
	};


	/* Creates a region of interest around the given point for the
	   specified chamber.
	 */
	void CreateToContain(const ClusterPoint& point, const ChamberID chamber);

	/* Extends the region of interest to contain the specified point.
	 */
	void ExpandToContain(const ClusterPoint& point);

	/* Creates a region of interest around all the given points and for the
	   specified chamber.
	 */
	void CreateToContain(const ClusterPoint* points, const UInt count, const ChamberID chamber);

	/* Checks if the region of interest is within the boundaries imposed on
	   the specific chamber plane. This boundary is aproximately the square
	   box around the chamber's detection region.
	 */
	bool InBounds();

	/* Encodes the region of interest into a 32 bit code.
	 */
	ROI Encode() const;
	
	/* Encodes the region of interest into a 32 bit code, and returns the
	   hierarchal level the region was encoded at and the left and right
	   grid coordinate of the bottom left corner of the region boundary box.
	 */
	ROI Encode(UChar& level, UInt& l, UInt& b) const;

	/* Decodes a 32 bit region of interest code into this region of interest
	   object.
	 */
	void Decode(const ROI code);

	/* Decodes the chamber number of the region of interest 32 bit code.
	 */
	static ChamberID DecodeChamber(const ROI code);

	/* Decodes the 32 bit region of interest code into the chamber number,
	   hierarchal level, left and right grid coordinates of the region
	   boundary box. 
	 */
	static void Decode(const ROI code, ChamberID& chamber, UChar& level, UInt& l, UInt& b);

	/* Returns the chamber number of the region of interest.
	 */
	ChamberID Chamber() const { return chamber; };
	
	/* Returns the left hand boundary of the region of interest.
	 */
	Float Left() const    { return left; };
	
	/* Returns the right hand boundary of the region of interest.
	 */
	Float Right() const   { return right; };
	
	/* Returns the bottom boundary of the region of interest.
	 */
	Float Bottom() const  { return bottom; };
	
	/* Returns the top boundary of the region of interest.
	 */
	Float Top() const     { return top; };


	/* Typecast operator for implicit typecasing to 32 bit ROI codes.
	 */
	operator ROI () const { return Encode(); };


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
	inline void ConvertBackFromGrid(const UInt l, const UInt r, const UInt b, const UInt t);

	/* Internal method for decoding 32 bit region codes. This method is
	   called by the Decode methods.
	 */
	static void DecodeBits(const ROI code, ChamberID& chamber, UChar& colevel, UInt& l, UInt& b);


	// Boundary box scale numbers for each chamber. These are the boundary
	// boxes around the chambers detection surface.
	static Float planescale[NUMBER_OF_TRACKING_CHAMBERS];

	static UInt indexoffsets[13];  // Offset numbers used in the encoding and decoding process.


	ChamberID chamber; // Specifies the chamber the region of interest is on.
	Float left;        // Left boundary of boundary box.
	Float right;       // Right boundary of boundary box.
	Float bottom;      // Bottom boundary of boundary box.
	Float top;         // Top boundary of boundary box.
};


} // dHLT

#endif // dHLT_REGION_OF_INTEREST_HPP
