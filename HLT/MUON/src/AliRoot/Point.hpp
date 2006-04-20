////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONPOINT_H
#define ALIHLTMUONPOINT_H

#include <TObject.h>
#include <Riostream.h>


class AliHLTMUONPoint : public TObject
{
	// ostream operator usefull for text output.
	friend ostream& operator << (ostream& os, const AliHLTMUONPoint& p);

public:

	/* Default constructor initialises everything to zero.
	 */
	AliHLTMUONPoint();
	
	/* Constructs the point from the specified x and y coordinates.
	 */
	AliHLTMUONPoint(Float_t x, Float_t y);

	virtual ~AliHLTMUONPoint() {};
	

	Float_t fX;    // X coordinate of the 2D point.
	Float_t fY;    // Y coordinate of the 2D point.

	ClassDef(AliHLTMUONPoint, 1)  // A 2D space point.
};


#endif // ALIHLTMUONPOINT_H
