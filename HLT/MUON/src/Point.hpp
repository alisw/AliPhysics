////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCOREPOINT_H
#define ALIHLTMUONCOREPOINT_H

#include "BasicTypes.hpp"


/* A 2D point structure using floats.
   These are used to store impact points on the trigger chambers and 
   cluster centroids.
 */
class AliHLTMUONCorePoint
{
public:

	Float fX, fY;

	AliHLTMUONCorePoint()
	{
		fX = 0.0;
		fY = 0.0;
	}

	AliHLTMUONCorePoint(Float x, Float y)
	{
		fX = x;
		fY = y;
	}
};


#endif // ALIHLTMUONCOREPOINT_H
