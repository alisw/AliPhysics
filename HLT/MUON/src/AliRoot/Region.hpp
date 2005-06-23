////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_REGION_HPP
#define dHLT_ALIROOT_REGION_HPP

#include <TObject.h>
#include <Riostream.h>


namespace AliMUONHLT
{

class Point;


class Region : public TObject
{
public:

	/* Default constructor initialises everything to zero.
	 */
	Region();
	
	/* Constructs the boundary box from the specified boundary parameters.
	   Note: the following conditions must apply: left <= right and bottom <= top
	   If this is not the case then an error message is shown and everything
	   is set to zero.
	 */
	Region(const Float_t left, const Float_t right, const Float_t bottom, const Float_t top);

	virtual ~Region() {}
	
	/* Get/Set methods for the boundary lines.
	   Note: before assignment we check for the following parameter consistency:
	   fLeft <= fRight and fBottom <= fTop.
	   If this is condition cannot be achieved then an error message is shown
	   and nothing is changed internally.
	 */
	void Left(const Float_t value);
	Float_t Left() const { return fLeft; }
	void Right(const Float_t value);
	Float_t Right() const { return fRight; }
	void Bottom(const Float_t value);
	Float_t Bottom() const { return fBottom; }
	void Top(const Float_t value);
	Float_t Top() const { return fTop; }
	
	/* Checks if the point is within this region. If it is then kTRUE is returned
	   otherwise kFALSE is returned.
	 */
	Bool_t Contains(const Point& p) const;

	// ostream operator usefull for text output.
	friend ostream& operator << (ostream& os, const Region& r);

private:

	Float_t fLeft;    // Left boundary of boundary box.
	Float_t fRight;   // Right boundary of boundary box.
	Float_t fBottom;  // Bottom boundary of boundary box.
	Float_t fTop;     // Top boundary of boundary box.

	ClassDef(Region, 1)  // A boundary box region.
};


} // AliMUONHLT

#endif // dHLT_ALIROOT_REGION_HPP
