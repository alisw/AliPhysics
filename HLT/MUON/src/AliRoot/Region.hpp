////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONREGION_H
#define ALIHLTMUONREGION_H

#include <TObject.h>
#include <Riostream.h>

class AliHLTMUONPoint;


class AliHLTMUONRegion : public TObject
{
	// ostream operator usefull for text output.
	friend ostream& operator << (ostream& os, const AliHLTMUONRegion& r);

public:

	/* Default constructor initialises everything to zero.
	 */
	AliHLTMUONRegion();
	
	/* Constructs the boundary box from the specified boundary parameters.
	   Note: the following conditions must apply: left <= right and bottom <= top
	   If this is not the case then an error message is shown and everything
	   is set to zero.
	 */
	AliHLTMUONRegion(Float_t left, Float_t right, Float_t bottom, Float_t top);

	virtual ~AliHLTMUONRegion() {}
	
	/* Get/Set methods for the boundary lines.
	   Note: before assignment we check for the following parameter consistency:
	   fLeft <= fRight and fBottom <= fTop.
	   If this is condition cannot be achieved then an error message is shown
	   and nothing is changed internally.
	 */
	void Left(Float_t value);
	Float_t Left() const { return fLeft; }
	void Right(Float_t value);
	Float_t Right() const { return fRight; }
	void Bottom(Float_t value);
	Float_t Bottom() const { return fBottom; }
	void Top(Float_t value);
	Float_t Top() const { return fTop; }
	
	/* Checks if the point is within this region. If it is then kTRUE is returned
	   otherwise kFALSE is returned.
	 */
	Bool_t Contains(const AliHLTMUONPoint& p) const;

private:

	Float_t fLeft;    // Left boundary of boundary box.
	Float_t fRight;   // Right boundary of boundary box.
	Float_t fBottom;  // Bottom boundary of boundary box.
	Float_t fTop;     // Top boundary of boundary box.

	ClassDef(AliHLTMUONRegion, 1)  // A boundary box region.
};


#endif // ALIHLTMUONREGION_H
