////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_POINT_HPP
#define dHLT_ALIROOT_POINT_HPP

#include <TObject.h>
#include <Riostream.h>


namespace AliMUONHLT
{


class Point : public TObject
{
public:

	/* Default constructor initialises everything to zero.
	 */
	Point();
	
	/* Constructs the point from the specified x and y coordinates.
	 */
	Point(Float_t x, Float_t y);

	virtual ~Point() {};
	
	// ostream operator usefull for text output.
	friend ostream& operator << (ostream& os, const Point& p);


	Float_t fX;    // X coordinate of the 2D point.
	Float_t fY;    // Y coordinate of the 2D point.

	ClassDef(Point, 1)  // A 2D space point.
};


} // AliMUONHLT

#endif // dHLT_ALIROOT_POINT_HPP
