////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_POINT_HPP
#define dHLT_POINT_HPP

#include "BasicTypes.hpp"

namespace dHLT
{


/* A 2D point structure using floats.
   These are used to store impact points on the trigger chambers and 
   cluster centroids.
 */
class Point
{
public:

	Float x, y;

	Point()
	{
		x = 0.0;
		y = 0.0;
	};

	Point(const Float x, const Float y)
	{
		this->x = x;
		this->y = y;
	};
};


} // dHLT

#endif // dHLT_POINT_HPP
