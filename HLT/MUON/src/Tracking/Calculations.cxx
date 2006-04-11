////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Calculations.hpp"
#include <math.h>


Float AliHLTMUONCoreCalculateSignedPt(
		register Float x1,
		register Float y1, register Float y2,
		register Float z1, register Float z2,
		Float& p
	)
{
	register Float qBL = 3.0;
	register Float zf = 975.0;
	Float theta_times_zf = - (y1*z2 - y2*z1) / (z2-z1);
	Float xf = x1 * zf / z1;
	Float yf = y2 - ((y2-y1) * (z2-zf)) / (z2-z1);

	Float p_div_zf = (qBL * 0.3 / theta_times_zf);
	p = (Float) fabs( p_div_zf * zf );
	
	// Note: the 0.3 is a conversion factor to GeV. it is 1e9 / c where c is
	// the speed of light.
	Float pt = p_div_zf * sqrt(xf*xf+yf*yf);
	return pt;
};


Float AliHLTMUONCoreCalculateSignedPt(
		register Float x1,
		register Float y1, register Float y2,
		register Float z1, register Float z2,
		register Float zf, register Float qBL,
		Float& p
	)
{
	Float theta_times_zf = - (y1*z2 - y2*z1) / (z2-z1);
	Float xf = x1 * zf / z1;
	Float yf = y2 - ((y2-y1) * (z2-zf)) / (z2-z1);

	Float p_div_zf = (qBL * 0.3 / theta_times_zf);
	p = (Float) fabs( p_div_zf * zf );

	// Note: the 0.3 is a conversion factor to GeV. it is 1e9 / c where c is
	// the speed of light.
	Float pt = p_div_zf * sqrt(xf*xf+yf*yf);
	return pt;
};


Float AliHLTMUONCoreCalculateSignedPt(
		register Float x1,
		register Float y1, register Float y2,
		register Float z1, register Float z2
	)
{
	register Float qBL = 3.0;
	register Float zf = 975.0;
	Float theta_times_zf = - (y1*z2 - y2*z1) / (z2-z1);
	Float xf = x1 * zf / z1;
	Float yf = y2 - ((y2-y1) * (z2-zf)) / (z2-z1);

	// Note: the 0.3 is a conversion factor to GeV. it is 1e9 / c where c is
	// the speed of light.
	Float pt = (qBL * 0.3 / theta_times_zf) * sqrt(xf*xf+yf*yf);
	return pt;
};


Float AliHLTMUONCoreCalculateSignedPt(
		register Float x1,
		register Float y1, register Float y2,
		register Float z1, register Float z2,
		register Float zf, register Float qBL
	)
{
	Float theta_times_zf = - (y1*z2 - y2*z1) / (z2-z1);
	Float xf = x1 * zf / z1;
	Float yf = y2 - ((y2-y1) * (z2-zf)) / (z2-z1);

	// Note: the 0.3 is a conversion factor to GeV. it is 1e9 / c where c is
	// the speed of light.
	Float pt = (qBL * 0.3 / theta_times_zf) * sqrt(xf*xf+yf*yf);
	return pt;
};


Float AliHLTMUONCoreCalculatePt(
		register Float x1,
		register Float y1, register Float y2,
		register Float z1, register Float z2
	)
{
	return (Float) fabs(AliHLTMUONCoreCalculateSignedPt(x1, y1, y2, z1, z2));
};


Float AliHLTMUONCoreCalculatePt(
		register Float x1,
		register Float y1, register Float y2,
		register Float z1, register Float z2,
		register Float zf, register Float qBL
	)
{
	return (Float) fabs(AliHLTMUONCoreCalculateSignedPt(x1, y1, y2, z1, z2, zf, qBL));
};

