////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "RegionOfInterest.hpp"
#include <math.h>


Float AliHLTMUONCoreRegionOfInterest::fgPlaneScale[gkNUMBER_OF_TRACKING_CHAMBERS]
///*
	= {	100.0f,
		100.0f,
		130.0f,
		130.0f,
		180.0f,
		180.0f,
		240.0f,
		240.0f,
		250.0f,
		250.0f
	};
//*/
/*
	= {	101.397718695768f,
		103.868515964831f,
		128.956611312237f,
		131.807531238079f,
		183.219120567424f,
		187.400469791991f,
		237.861752479389f,
		242.993408345904f,
		269.221871663647f,
		274.353527530162f
	};
*/


void AliHLTMUONCoreRegionOfInterest::CreateToContain(
		const AliHLTMUONCoreClusterPoint& point, AliHLTMUONCoreChamberID chamber
	)
{
	Assert( 0 <= chamber && chamber < gkNUMBER_OF_TRACKING_CHAMBERS );
	fChamber = chamber;
	fLeft = fRight = point.fX;
	fBottom = fTop = point.fY;
}


void AliHLTMUONCoreRegionOfInterest::ExpandToContain(const AliHLTMUONCoreClusterPoint& point)
{
	if (point.fX < fLeft)
		fLeft = point.fX;
	else
		if (point.fX > fRight) fRight = point.fX;
	if (point.fY < fBottom)
		fBottom = point.fY;
	else
		if (point.fY > fTop) fTop = point.fY;
}


void AliHLTMUONCoreRegionOfInterest::CreateToContain(
		const AliHLTMUONCoreClusterPoint* points, UInt count,
		AliHLTMUONCoreChamberID chamber
	)
{
	Assert( 0 <= chamber && chamber < gkNUMBER_OF_TRACKING_CHAMBERS );
	Assert( count > 0 );
	
	CreateToContain(points[0], chamber);
	for (UInt i = 1; i < count; i++)
		ExpandToContain(points[i]);
}


bool AliHLTMUONCoreRegionOfInterest::InBounds()
{
	Assert( 0 <= fChamber && fChamber < gkNUMBER_OF_TRACKING_CHAMBERS );
	register Float bound = fgPlaneScale[fChamber];
	return -bound <= fLeft
	  && fRight <= bound
	  && -bound <= fBottom
	  && fTop <= bound;
}


#define MAX_LEVELS           13
#define GRID_SIZE            (2 << MAX_LEVELS)

// The following numbers specifiy the total indices for lower levels.
// Computed with: x_level = sum from 0 to level ( ( 2^level - 1 )^2 )
#define MAX_INDICES          357848422
#define INDICES_TO_LEVEL_13  89445733
#define INDICES_TO_LEVEL_12  22353252
#define INDICES_TO_LEVEL_11  5584227
#define INDICES_TO_LEVEL_10  1394018
#define INDICES_TO_LEVEL_9   347489
#define INDICES_TO_LEVEL_8   86368
#define INDICES_TO_LEVEL_7   21343
#define INDICES_TO_LEVEL_6   5214
#define INDICES_TO_LEVEL_5   1245
#define INDICES_TO_LEVEL_4   284
#define INDICES_TO_LEVEL_3   59
#define INDICES_TO_LEVEL_2   10
#define INDICES_TO_LEVEL_1   1


UInt AliHLTMUONCoreRegionOfInterest::fgIndexOffsets[13]
	= {	INDICES_TO_LEVEL_1,
		INDICES_TO_LEVEL_2,
		INDICES_TO_LEVEL_3,
		INDICES_TO_LEVEL_4,
		INDICES_TO_LEVEL_5,
		INDICES_TO_LEVEL_6,
		INDICES_TO_LEVEL_7,
		INDICES_TO_LEVEL_8,
		INDICES_TO_LEVEL_9,
		INDICES_TO_LEVEL_10,
		INDICES_TO_LEVEL_11,
		INDICES_TO_LEVEL_12,
		INDICES_TO_LEVEL_13
	};


inline void AliHLTMUONCoreRegionOfInterest::ConvertToGrid(
		register UInt& l, register UInt& r, register UInt& b, register UInt& t
	) const
{
	Assert( 0 <= fChamber && fChamber < gkNUMBER_OF_TRACKING_CHAMBERS );

	register Float scale = fgPlaneScale[fChamber];
	l = (UInt) floor( (fLeft / scale + 1.0f) * 0.5f * GRID_SIZE );
	r = (UInt) ceil( (fRight / scale + 1.0f) * 0.5f * GRID_SIZE );
	b = (UInt) floor( (fBottom / scale + 1.0f) * 0.5f * GRID_SIZE );
	t = (UInt) ceil( (fTop / scale + 1.0f) * 0.5f * GRID_SIZE );

	/* 
	// Not required with proper coding. This case is handled properly by
	// encoding as a global boundary box.
	Assert( l <= GRID_SIZE );
	Assert( r <= GRID_SIZE );
	Assert( b <= GRID_SIZE );
	Assert( t <= GRID_SIZE );
	*/
}


inline void AliHLTMUONCoreRegionOfInterest::ConvertBackFromGrid(
		register UInt l, register UInt r, register UInt b, register UInt t
	)
{
	Assert( l <= GRID_SIZE );
	Assert( r <= GRID_SIZE );
	Assert( b <= GRID_SIZE );
	Assert( t <= GRID_SIZE );

	Assert( 0 <= fChamber && fChamber < gkNUMBER_OF_TRACKING_CHAMBERS );
	register Float scale = fgPlaneScale[fChamber];
	fLeft = ((Float)l / (Float)GRID_SIZE - 0.5f) * 2.0f * scale;
	fRight = ((Float)r / (Float)GRID_SIZE - 0.5f) * 2.0f * scale;
	fBottom = ((Float)b / (Float)GRID_SIZE - 0.5f) * 2.0f * scale;
	fTop = ((Float)t / (Float)GRID_SIZE - 0.5f) * 2.0f * scale;
}


AliHLTMUONCoreROI AliHLTMUONCoreRegionOfInterest::Encode() const
{
	UInt l, r, b, t, maxwidth, index;

	ConvertToGrid(l, r, b, t);

	// Work out which resolution level the location of the ROI
	// boundary box needs to be coded at. The higher the level
	// the smaller the box size: i.e. box_size = scale / 2^(level)
	// We use a binary search type method for finding the level.
	// More specificaly we search for which n,
	//   (l * 2^(-n) + 2) * 2^n >= r and (b * 2^(-n) + 2) * 2^n >= t
	
	if ( ((l >> 6) + 2) << 6 >= r && ((b >> 6) + 2) << 6 >= t )
	{
		if ( ((l >> 2) + 2) << 2 >= r && ((b >> 2) + 2) << 2 >= t )
		{
			if ( l + 2 >= r && b + 2 >= t )
			{
				index = INDICES_TO_LEVEL_13;
				maxwidth = (2 << 13) - 1;
			}
			else
			{
				if ( ((l >> 1) + 2) << 1 >= r && ((b >> 1) + 2) << 1 >= t )
				{
					index = INDICES_TO_LEVEL_12;
					maxwidth = (2 << 12) - 1;
					l >>= 1; b >>= 1;
				}
				else
				{
					index = INDICES_TO_LEVEL_11;
					maxwidth = (2 << 11) - 1;
					l >>= 2; b >>= 2;
				}
			}
		}
		else
		{
			if ( ((l >> 4) + 2) << 4 >= r && ((b >> 4) + 2) << 4 >= t )
			{
				if ( ((l >> 3) + 2) << 3 >= r && ((b >> 3) + 2) << 3 >= t )
				{
					index = INDICES_TO_LEVEL_10;
					maxwidth = (2 << 10) - 1;
					l >>= 3; b >>= 3;
				}
				else
				{
					index = INDICES_TO_LEVEL_9;
					maxwidth = (2 << 9) - 1;
					l >>= 4; b >>= 4;
				}
			}
			else
			{
				if ( ((l >> 5) + 2) << 5 >= r && ((b >> 5) + 2) << 5 >= t )
				{
					index = INDICES_TO_LEVEL_8;
					maxwidth = (2 << 8) - 1;
					l >>= 5; b >>= 5;
				}
				else
				{
					index = INDICES_TO_LEVEL_7;
					maxwidth = (2 << 7) - 1;
					l >>= 6; b >>= 6;
				}
			}
		}
	}
	else
	{
		if ( ((l >> 10) + 2) << 10 >= r && ((b >> 10) + 2) << 10 >= t )
		{
			if ( ((l >> 8) + 2) << 8 >= r && ((b >> 8) + 2) << 8 >= t )
			{
				if ( ((l >> 7) + 2) << 7 >= r && ((b >> 7) + 2) << 7 >= t )
				{
					index = INDICES_TO_LEVEL_6;
					maxwidth = (2 << 6) - 1;
					l >>= 7; b >>= 7;
				}
				else
				{
					index = INDICES_TO_LEVEL_5;
					maxwidth = (2 << 5) - 1;
					l >>= 8; b >>= 8;
				}
			}
			else
			{
				if ( ((l >> 9) + 2) << 9 >= r && ((b >> 9) + 2) << 9 >= t )
				{
					index = INDICES_TO_LEVEL_4;
					maxwidth = (2 << 4) - 1;
					l >>= 9; b >>= 9;
				}
				else
				{
					index = INDICES_TO_LEVEL_3;
					maxwidth = (2 << 3) - 1;
					l >>= 10; b >>= 10;
				}
			}
		}
		else
		{
			if ( ((l >> 12) + 2) << 12 >= r && ((b >> 12) + 2) << 12 >= t )
			{
				if ( ((l >> 11) + 2) << 11 >= r && ((b >> 11) + 2) << 11 >= t )
				{
					index = INDICES_TO_LEVEL_2;
					maxwidth = (2 << 2) - 1;
					l >>= 11; b >>= 11;
				}
				else
				{
					index = INDICES_TO_LEVEL_1;
					maxwidth = (2 << 1) - 1;
					l >>= 12; b >>= 12;
				}
			}
			else
			{
				index = 0;
				maxwidth = 0;
				l = 0; b = 0;
			}
		}
	}
	
	// Make sure the ROI boundary box does not go outside the
	// global region of interest.
	if ( l > maxwidth - 1 ) l = maxwidth - 1;
	if ( b > maxwidth - 1 ) b = maxwidth - 1;
	
	return MAX_INDICES * fChamber + b * maxwidth + l + index;
}


AliHLTMUONCoreROI AliHLTMUONCoreRegionOfInterest::Encode(UChar& level, UInt& l, UInt& b) const
{
	UInt r, t, maxwidth, index;

	ConvertToGrid(l, r, b, t);

	// Work out which resolution level the location of the ROI
	// boundary box needs to be coded at. The higher the level
	// the smaller the box size: i.e. box_size = scale / 2^(level)
	// We use a binary search type method for finding the level.
	// More specificaly we search for which n,
	//   (l * 2^(-n) + 2) * 2^n >= r and (b * 2^(-n) + 2) * 2^n >= t
	
	if ( ((l >> 6) + 2) << 6 >= r && ((b >> 6) + 2) << 6 >= t )
	{
		if ( ((l >> 2) + 2) << 2 >= r && ((b >> 2) + 2) << 2 >= t )
		{
			if ( l + 2 >= r && b + 2 >= t )
			{
				level = 13;
				index = INDICES_TO_LEVEL_13;
				maxwidth = (2 << 13) - 1;
			}
			else
			{
				if ( ((l >> 1) + 2) << 1 >= r && ((b >> 1) + 2) << 1 >= t )
				{
					level = 12;
					index = INDICES_TO_LEVEL_12;
					maxwidth = (2 << 12) - 1;
					l >>= 1; b >>= 1;
				}
				else
				{
					level = 11;
					index = INDICES_TO_LEVEL_11;
					maxwidth = (2 << 11) - 1;
					l >>= 2; b >>= 2;
				}
			}
		}
		else
		{
			if ( ((l >> 4) + 2) << 4 >= r && ((b >> 4) + 2) << 4 >= t )
			{
				if ( ((l >> 3) + 2) << 3 >= r && ((b >> 3) + 2) << 3 >= t )
				{
					level = 10;
					index = INDICES_TO_LEVEL_10;
					maxwidth = (2 << 10) - 1;
					l >>= 3; b >>= 3;
				}
				else
				{
					level = 9;
					index = INDICES_TO_LEVEL_9;
					maxwidth = (2 << 9) - 1;
					l >>= 4; b >>= 4;
				}
			}
			else
			{
				if ( ((l >> 5) + 2) << 5 >= r && ((b >> 5) + 2) << 5 >= t )
				{
					level = 8;
					index = INDICES_TO_LEVEL_8;
					maxwidth = (2 << 8) - 1;
					l >>= 5; b >>= 5;
				}
				else
				{
					level = 7;
					index = INDICES_TO_LEVEL_7;
					maxwidth = (2 << 7) - 1;
					l >>= 6; b >>= 6;
				}
			}
		}
	}
	else
	{
		if ( ((l >> 10) + 2) << 10 >= r && ((b >> 10) + 2) << 10 >= t )
		{
			if ( ((l >> 8) + 2) << 8 >= r && ((b >> 8) + 2) << 8 >= t )
			{
				if ( ((l >> 7) + 2) << 7 >= r && ((b >> 7) + 2) << 7 >= t )
				{
					level = 6;
					index = INDICES_TO_LEVEL_6;
					maxwidth = (2 << 6) - 1;
					l >>= 7; b >>= 7;
				}
				else
				{
					level = 5;
					index = INDICES_TO_LEVEL_5;
					maxwidth = (2 << 5) - 1;
					l >>= 8; b >>= 8;
				}
			}
			else
			{
				if ( ((l >> 9) + 2) << 9 >= r && ((b >> 9) + 2) << 9 >= t )
				{
					level = 4;
					index = INDICES_TO_LEVEL_4;
					maxwidth = (2 << 4) - 1;
					l >>= 9; b >>= 9;
				}
				else
				{
					level = 3;
					index = INDICES_TO_LEVEL_3;
					maxwidth = (2 << 3) - 1;
					l >>= 10; b >>= 10;
				}
			}
		}
		else
		{
			if ( ((l >> 12) + 2) << 12 >= r && ((b >> 12) + 2) << 12 >= t )
			{
				if ( ((l >> 11) + 2) << 11 >= r && ((b >> 11) + 2) << 11 >= t )
				{
					level = 2;
					index = INDICES_TO_LEVEL_2;
					maxwidth = (2 << 2) - 1;
					l >>= 11; b >>= 11;
				}
				else
				{
					level = 1;
					index = INDICES_TO_LEVEL_1;
					maxwidth = (2 << 1) - 1;
					l >>= 12; b >>= 12;
				}
			}
			else
			{
				level = 0;
				index = 0;
				maxwidth = 0;
				l = 0; b = 0;
			}
		}
	}
	
	// Make sure the ROI boundary box does not go outside the
	// global region of interest.
	if ( l > maxwidth - 1 ) l = maxwidth - 1;
	if ( b > maxwidth - 1 ) b = maxwidth - 1;
	
	return MAX_INDICES * fChamber + b * maxwidth + l + index;
}


void AliHLTMUONCoreRegionOfInterest::Decode(AliHLTMUONCoreROI code)
{
	UInt l, r, b, t;
	UChar colevel;
	DecodeBits(code, fChamber, colevel, l, b);

	// Complete decoding of bottom left corner.
	b = b << colevel;
	l = l << colevel;

	// Decode top right corner of boundry box.
	r = l + (2 << colevel);
	t = b + (2 << colevel);

	ConvertBackFromGrid(l, r, b, t);
}


void AliHLTMUONCoreRegionOfInterest::Decode(
		AliHLTMUONCoreROI code, AliHLTMUONCoreChamberID& chamber,
		UChar& level, UInt& l, UInt& b
	)
{
	UChar colevel;
	DecodeBits(code, chamber, colevel, l, b);
	level = MAX_LEVELS - colevel;
}


void AliHLTMUONCoreRegionOfInterest::DecodeBits(
		AliHLTMUONCoreROI code, AliHLTMUONCoreChamberID& chamber,
		UChar& colevel, UInt& l, UInt& b
	)
{
	// First decode the chamber number and the remainder 
	// contains the location index.
	chamber = (AliHLTMUONCoreChamberID)(code / MAX_INDICES);
	UInt index = code % MAX_INDICES;
	
	UInt width;

	// Binary search which colevel the index is for. 
	if (index >= INDICES_TO_LEVEL_7)
	{
		if (index >= INDICES_TO_LEVEL_11)
		{
			if (index >= INDICES_TO_LEVEL_13)
			{
				index -= INDICES_TO_LEVEL_13;
				width = (2 << 13) - 1;
				colevel = 0;
			}
			else
			{
				if (index >= INDICES_TO_LEVEL_12)
				{
					index -= INDICES_TO_LEVEL_12;
					width = (2 << 12) - 1;
					colevel = 1;
				}
				else
				{
					index -= INDICES_TO_LEVEL_11;
					width = (2 << 11) - 1;
					colevel = 2;
				}
			}
		}
		else
		{
			if (index >= INDICES_TO_LEVEL_9)
			{
				if (index >= INDICES_TO_LEVEL_10)
				{
					index -= INDICES_TO_LEVEL_10;
					width = (2 << 10) - 1;
					colevel = 3;
				}
				else
				{
					index -= INDICES_TO_LEVEL_9;
					width = (2 << 9) - 1;
					colevel = 4;
				}
			}
			else
			{
				if (index >= INDICES_TO_LEVEL_8)
				{
					index -= INDICES_TO_LEVEL_8;
					width = (2 << 8) - 1;
					colevel = 5;
				}
				else
				{
					index -= INDICES_TO_LEVEL_7;
					width = (2 << 7) - 1;
					colevel = 6;
				}
			}
		}
	}
	else
	{
		if (index >= INDICES_TO_LEVEL_3)
		{
			if (index >= INDICES_TO_LEVEL_5)
			{
				if (index >= INDICES_TO_LEVEL_6)
				{
					index -= INDICES_TO_LEVEL_6;
					width = (2 << 6) - 1;
					colevel = 7;
				}
				else
				{
					index -= INDICES_TO_LEVEL_5;
					width = (2 << 5) - 1;
					colevel = 8;
				}
			}
			else
			{
				if (index >= INDICES_TO_LEVEL_4)
				{
					index -= INDICES_TO_LEVEL_4;
					width = (2 << 4) - 1;
					colevel = 9;
				}
				else
				{
					index -= INDICES_TO_LEVEL_3;
					width = (2 << 3) - 1;
					colevel = 10;
				}
			}
		}
		else
		{
			if (index >= INDICES_TO_LEVEL_1)
			{
				if (index >= INDICES_TO_LEVEL_2)
				{
					index -= INDICES_TO_LEVEL_2;
					width = (2 << 2) - 1;
					colevel = 11;
				}
				else
				{
					index -= INDICES_TO_LEVEL_1;
					width = (2 << 1) - 1;
					colevel = 12;
				}
			}
			else
			{
				width = 1;
				colevel = 13;
			}
		}
	}

	// Can now decode the x, y position of the bottom left corner
	// of the ROI boundary box.
	b = (index / width);
	l = (index % width);
}


AliHLTMUONCoreChamberID AliHLTMUONCoreRegionOfInterest::DecodeChamber(AliHLTMUONCoreROI code)
{
	return (AliHLTMUONCoreChamberID)(code / MAX_INDICES);
}

