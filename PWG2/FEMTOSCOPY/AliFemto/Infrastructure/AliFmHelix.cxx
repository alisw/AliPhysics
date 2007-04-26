/***************************************************************************
*
* $Id$
*
* Author: Thomas Ullrich, Sep 1997
***************************************************************************
*
* Description: Parametrization of a helix
* 
***************************************************************************
*
* $Log$
* Revision 1.1.1.1  2007/04/25 15:38:41  panos
* Importing the HBT code dir
*
* Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
* First version on CVS
*
* Revision 1.26  2005/10/13 23:15:13  genevb
* Save a few calculations
*
* Revision 1.25  2005/10/13 22:25:35  genevb
* pathLength to plane now finds nearest approach to intersection regardless of # of loops
*
* Revision 1.24  2005/07/06 18:49:56  fisyak
* Replace AliFmHelixD, AliFmLorentzVectorD,AliFmLorentzVectorF,AliFmMatrixD,AliFmMatrixF,AliFmPhysicalHelixD,AliFmThreeVectorD,AliFmThreeVectorF by templated version
*
* Revision 1.23  2004/12/02 02:51:16  ullrich
* Added option to pathLenghth() and distance() to search for
* DCA only within one period. Default stays as it was.
*
* Revision 1.22  2004/05/03 23:35:31  perev
* Possible non init WarnOff
*
* Revision 1.21  2004/01/27 02:49:48  perev
* Big value appropriate for float
*
* Revision 1.20  2003/12/22 18:59:36  perev
* remove test for only +ve pathLeng added before
*
* Revision 1.19  2003/12/18 17:27:02  perev
* Small bug fix in number of iters. i++ ==> i--
*
* Revision 1.18  2003/10/30 20:06:46  perev
* Check of quality added
*
* Revision 1.17  2003/10/19 20:17:00  perev
* Protection agains overfloat added into pathLength(AliFmThreeVector,AliFmThreeVector)
*
* Revision 1.16  2003/10/06 23:39:21  perev
* sqrt(-ve) == no solution. infinity returns
*
* Revision 1.15  2003/09/02 17:59:34  perev
* gcc 3.2 updates + WarnOff
*
* Revision 1.14  2003/06/26 17:15:56  ullrich
* Changed local variable name in pathLenght.
*
* Revision 1.13  2002/06/21 17:49:25  genevb
* Some minor speed improvements
*
* Revision 1.12  2002/04/24 02:40:25  ullrich
* Restored old format and lost CVS statements.
*
* Revision 1.11  2002/02/12 19:37:51  jeromel
* fix for Linux 7.2 (float.h). Took oportunity to doxygenize.
*
* Revision 1.10  2000/07/17 21:44:19  ullrich
* Fixed problem in pathLength(cyl_rad).
*
* Revision 1.9  2000/05/22 21:38:28  ullrich
* Add parenthesis to make Linux compiler happy.
*
* Revision 1.8  2000/05/22 21:11:21  ullrich
* In pathLength(AliFmThreeVector&): Increased number of max iteration
* in Newton method from 10 to 100. Improved initial guess in case
* it is off by n period.
*
* Revision 1.7  2000/03/06 20:24:25  ullrich
* Parameter h for case B=0 correctly handled now.
*
* Revision 1.6  1999/12/22 15:14:39  ullrich
* Added analytical solution for dca between two helices
* in the case for B=0.
*
* Revision 1.5  1999/12/21 15:14:08  ullrich
* Modified to cope with new compiler version on Sun (CC5.0).
*
* Revision 1.4  1999/11/29 21:45:38  fisyak
* fix abs for HP
*
* Revision 1.3  1999/03/07 14:55:41  wenaus
* fix scope problem
*
* Revision 1.2  1999/03/02 19:47:35  ullrich
* Added method to find dca between two helices
*
* Revision 1.1  1999/01/30 03:59:02  fisyak
* Root Version of AliFmarClassLibrary
*
* Revision 1.1  1999/01/23 00:29:15  ullrich
* Initial Revision
*
**************************************************************************/
#if !defined(ST_NO_NUMERIC_LIMITS)
#    include <limits>
#    if !defined(ST_NO_NAMESPACES)
using std::numeric_limits;
#    endif
#endif
#define FOR_HELIX
#include <float.h>
#include <assert.h>

#include "AliFmHelix.h"
#include "Base/PhysicalConstants.h" 
#include "Base/SystemOfUnits.h"
#ifdef __ROOT__
ClassImpT(AliFmHelix,double);
#endif

#ifdef WIN32
#include "gcc2vs.h"
#endif

#include <iostream>
#include <fstream>
using namespace std;

const double AliFmHelix::NoSolution = 3.e+33;

AliFmHelix::AliFmHelix(){ /*noop*/ }

AliFmHelix::AliFmHelix(double c, double d, double phase,
				 const AliFmThreeVector<double>& o, int h)
{
	setParameters(c, d, phase, o, h);
}

AliFmHelix::~AliFmHelix() { /* noop */ };

void AliFmHelix::setParameters(double c, double dip, double phase,
							const AliFmThreeVector<double>& o, int h)
{
	//
	//  The order in which the parameters are set is important
	//  since setCurvature might have to adjust the others.
	//
	mH = (h>=0) ? 1 : -1;    // Default is: positive particle
	//             positive field
	mOrigin   = o;
	setDipAngle(dip);
	setPhase(phase);

	//
	// Check for singularity and correct for negative curvature.           
	// May change mH and mPhase. Must therefore be set last.
	//
	setCurvature(c);

	//
	// For the case B=0, h is ill defined. In the following we
	// always assume h = +1. Since phase = psi - h * pi/2
	// we have to correct the phase in case h = -1.
	// This assumes that the user uses the same h for phase
	// as the one he passed to the constructor.
	//
	if (mSingularity && mH == -1) {
		mH = +1;
		setPhase(mPhase-M_PI);
	}
}

void AliFmHelix::setCurvature(double val)
{
	if (val < 0) {
		mCurvature = -val;
		mH = -mH;
		setPhase(mPhase+M_PI);
	}
	else
		mCurvature = val;

#ifndef ST_NO_NUMERIC_LIMITS
	if (fabs(mCurvature) <= numeric_limits<double>::epsilon())
#else
	if (fabs(mCurvature) <= static_cast<double>(0))
#endif    
		mSingularity = true;			// straight line
	else
		mSingularity = false;            	// curved
}

void AliFmHelix::setPhase(double val)
{
	mPhase       = val;
	mCosPhase    = cos(mPhase);
	mSinPhase    = sin(mPhase);
	if (fabs(mPhase) > M_PI)
		mPhase = atan2(mSinPhase, mCosPhase);  // force range [-pi,pi]
}

void AliFmHelix::setDipAngle(double val)
{
	mDipAngle    = val;
	mCosDipAngle = cos(mDipAngle);
	mSinDipAngle = sin(mDipAngle);
}

double AliFmHelix::xcenter() const
{
	if (mSingularity)
		return 0;
	else
		return mOrigin.x()-mCosPhase/mCurvature;
}

double AliFmHelix::ycenter() const
{
	if (mSingularity)
		return 0;
	else
		return mOrigin.y()-mSinPhase/mCurvature;
}

double AliFmHelix::fudgePathLength(const AliFmThreeVector<double>& p) const
{
	double s;
	double dx = p.x()-mOrigin.x();
	double dy = p.y()-mOrigin.y();

	if (mSingularity) {
		s = (dy*mCosPhase - dx*mSinPhase)/mCosDipAngle;
	}
	else {
		s = atan2(dy*mCosPhase - dx*mSinPhase,
			1/mCurvature + dx*mCosPhase+dy*mSinPhase)/
			(mH*mCurvature*mCosDipAngle);
	}
	return s;
}

double AliFmHelix::distance(const AliFmThreeVector<double>& p, bool scanPeriods) const
{
	return abs(this->at(pathLength(p,scanPeriods))-p);
}

double AliFmHelix::pathLength(const AliFmThreeVector<double>& p, bool scanPeriods) const 
{
	//
	//  Returns the path length at the distance of closest 
	//  approach between the helix and point p. 
	//  For the case of B=0 (straight line) the path length
	//  can be calculated analytically. For B>0 there is
	//  unfortunately no easy solution to the problem.
	//  Here we use the Newton method to find the root of the
	//  referring equation. The 'fudgePathLength' serves
	//  as a starting value.
	//
	double s;
	double dx = p.x()-mOrigin.x();
	double dy = p.y()-mOrigin.y();
	double dz = p.z()-mOrigin.z();

	if (mSingularity) {
		s = mCosDipAngle*(mCosPhase*dy-mSinPhase*dx) +
			mSinDipAngle*dz;
	}
	else { //
#ifndef ST_NO_NAMESPACES
		{
			using namespace units;
#endif
			const double MaxPrecisionNeeded = micrometer;
			const int    MaxIterations      = 100;

			//
			// The math is taken from Maple with C(expr,optimized) and
			// some hand-editing. It is not very nice but efficient.
			//
			double t34 = mCurvature*mCosDipAngle*mCosDipAngle;
			double t41 = mSinDipAngle*mSinDipAngle;
			double t6, t7, t11, t12, t19;

			//
			// Get a first guess by using the dca in 2D. Since
			// in some extreme cases we might be off by n periods
			// we add (subtract) periods in case we get any closer.
			// 
			s = fudgePathLength(p);

			if (scanPeriods) {
				double ds = period();
				int    j, jmin = 0;
				double d, dmin = abs(at(s) - p);
				for(j=1; j<MaxIterations; j++) {
					if ((d = abs(at(s+j*ds) - p)) < dmin) {
						dmin = d;
						jmin = j;
					}
					else
						break;
				}
				for(j=-1; -j<MaxIterations; j--) {
					if ((d = abs(at(s+j*ds) - p)) < dmin) {
						dmin = d;
						jmin = j;
					}
					else
						break;
				}
				if (jmin) s += jmin*ds;
			}

			//
			// Newtons method:
			// Stops after MaxIterations iterations or if the required
			// precision is obtained. Whatever comes first.
			//
			double sOld = s;
			for (int i=0; i<MaxIterations; i++) {
				t6  = mPhase+s*mH*mCurvature*mCosDipAngle;
				t7  = cos(t6);
				t11 = dx-(1/mCurvature)*(t7-mCosPhase);
				t12 = sin(t6);
				t19 = dy-(1/mCurvature)*(t12-mSinPhase);
				s  -= (t11*t12*mH*mCosDipAngle-t19*t7*mH*mCosDipAngle -
					(dz-s*mSinDipAngle)*mSinDipAngle)/
					(t12*t12*mCosDipAngle*mCosDipAngle+t11*t7*t34 +
					t7*t7*mCosDipAngle*mCosDipAngle +
					t19*t12*t34+t41);
				if (fabs(sOld-s) < MaxPrecisionNeeded) break;
				sOld = s;
			}
#ifndef ST_NO_NAMESPACES
		}
#endif
	}
	return s;
}

double AliFmHelix::period() const
{
	if (mSingularity)
#ifndef ST_NO_NUMERIC_LIMITS
		return numeric_limits<double>::max();
#else
		return DBL_MAX;
#endif    
	else	
		return fabs(2*M_PI/(mH*mCurvature*mCosDipAngle)); 
}

pair<double, double> AliFmHelix::pathLength(double r) const
{
	pair<double,double> value;
	pair<double,double> VALUE(999999999.,999999999.);
	//
	// The math is taken from Maple with C(expr,optimized) and
	// some hand-editing. It is not very nice but efficient.
	// 'first' is the smallest of the two solutions (may be negative)
	// 'second' is the other.
	//
	if (mSingularity) {
		double t1 = mCosDipAngle*(mOrigin.x()*mSinPhase-mOrigin.y()*mCosPhase);
		double t12 = mOrigin.y()*mOrigin.y();
		double t13 = mCosPhase*mCosPhase;
		double t15 = r*r;
		double t16 = mOrigin.x()*mOrigin.x();
		double t20 = -mCosDipAngle*mCosDipAngle*(2.0*mOrigin.x()*mSinPhase*mOrigin.y()*mCosPhase +
			t12-t12*t13-t15+t13*t16);
		if (t20<0.) return VALUE;
		t20 = ::sqrt(t20);
		value.first  = (t1-t20)/(mCosDipAngle*mCosDipAngle);
		value.second = (t1+t20)/(mCosDipAngle*mCosDipAngle);
	}
	else {
		double t1 = mOrigin.y()*mCurvature;
		double t2 = mSinPhase;
		double t3 = mCurvature*mCurvature;
		double t4 = mOrigin.y()*t2;
		double t5 = mCosPhase;
		double t6 = mOrigin.x()*t5;
		double t8 = mOrigin.x()*mOrigin.x();
		double t11 = mOrigin.y()*mOrigin.y();
		double t14 = r*r;
		double t15 = t14*mCurvature;
		double t17 = t8*t8;
		double t19 = t11*t11;
		double t21 = t11*t3;
		double t23 = t5*t5;
		double t32 = t14*t14;
		double t35 = t14*t3;
		double t38 = 8.0*t4*t6 - 4.0*t1*t2*t8 - 4.0*t11*mCurvature*t6 +
			4.0*t15*t6 + t17*t3 + t19*t3 + 2.0*t21*t8 + 4.0*t8*t23 -
			4.0*t8*mOrigin.x()*mCurvature*t5 - 4.0*t11*t23 -
			4.0*t11*mOrigin.y()*mCurvature*t2 + 4.0*t11 - 4.0*t14 +
			t32*t3 + 4.0*t15*t4 - 2.0*t35*t11 - 2.0*t35*t8;
		double t40 = (-t3*t38);
		if (t40<0.) return VALUE;
		t40 = ::sqrt(t40);

		double t43 = mOrigin.x()*mCurvature;
		double t45 = 2.0*t5 - t35 + t21 + 2.0 - 2.0*t1*t2 -2.0*t43 - 2.0*t43*t5 + t8*t3;
		double t46 = mH*mCosDipAngle*mCurvature;

		value.first = (-mPhase + 2.0*atan((-2.0*t1 + 2.0*t2 + t40)/t45))/t46;
		value.second = -(mPhase + 2.0*atan((2.0*t1 - 2.0*t2 + t40)/t45))/t46;

		//
		//   Solution can be off by +/- one period, select smallest
		//
		double p = period();
		//	malisa deletes "isnan" check 22apr2006	if (!isnan(value.first)) {
		if (fabs(value.first-p) < fabs(value.first)) value.first = value.first-p;
		else if (fabs(value.first+p) < fabs(value.first)) value.first = value.first+p;
		//	malisa	}
		//	malisa deletes "isnan" check 22apr2006		if (!isnan(value.second)) {
		if (fabs(value.second-p) < fabs(value.second)) value.second = value.second-p;
		else if (fabs(value.second+p) < fabs(value.second)) value.second = value.second+p;
		//      malisa }
	}
	if (value.first > value.second)
		swap(value.first,value.second);
	return(value);
}

pair<double, double> AliFmHelix::pathLength(double r, double x, double y, bool scanPeriods)
{
	double x0 = mOrigin.x();
	double y0 = mOrigin.y();
	mOrigin.setX(x0-x);
	mOrigin.setY(y0-y);
	pair<double, double> result = this->pathLength(r);
	mOrigin.setX(x0);
	mOrigin.setY(y0);
	return result;  
}

double AliFmHelix::pathLength(const AliFmThreeVector<double>& r,
						   const AliFmThreeVector<double>& n) const
{
	//
	// Vector 'r' defines the position of the center and
	// vector 'n' the normal vector of the plane.
	// For a straight line there is a simple analytical
	// solution. For curvatures > 0 the root is determined
	// by Newton method. In case no valid s can be found
	// the max. largest value for s is returned.
	//
	double s;

	if (mSingularity) {
		double t = n.z()*mSinDipAngle +
			n.y()*mCosDipAngle*mCosPhase -
			n.x()*mCosDipAngle*mSinPhase;
		if (t == 0)
			s = NoSolution;
		else
			s = ((r - mOrigin)*n)/t;
	}
	else {
		const double MaxPrecisionNeeded = micrometer;
		const int    MaxIterations      = 20;

		double A = mCurvature*((mOrigin - r)*n) -
			n.x()*mCosPhase - 
			n.y()*mSinPhase;
		double t = mH*mCurvature*mCosDipAngle;
		double u = n.z()*mCurvature*mSinDipAngle;

		double a, f, fp;
		double sOld = s = 0;  
		double shiftOld = 0;
		double shift;
		//		(cos(angMax)-1)/angMax = 0.1
		const double angMax = 0.21;
		double deltas = fabs(angMax/(mCurvature*mCosDipAngle));
		//              dampingFactor = exp(-0.5);
		double dampingFactor = 0.60653;
		int i;

		for (i=0; i<MaxIterations; i++) {
			a  = t*s+mPhase;
			double sina = sin(a);
			double cosa = cos(a);
			f  = A +
				n.x()*cosa +
				n.y()*sina +
				u*s;
			fp = -n.x()*sina*t +
				n.y()*cosa*t +
				u;
			if ( fabs(fp)*deltas <= fabs(f) ) { //too big step
				int sgn = 1;
				if (fp<0.) sgn = -sgn;
				if (f <0.) sgn = -sgn;
				shift = sgn*deltas;
				if (shift == -shiftOld) { // don't get stuck shifting +/-deltas
					deltas *= dampingFactor; // dampen magnitude of shift
					shift = sgn*deltas;
					// allow iterations to run out
				} else {
					i--; // don't count against iterations
				}
			} else {
				shift = f/fp;
			}
			s -= shift;
			shiftOld = shift;
			if (fabs(sOld-s) < MaxPrecisionNeeded) break;
			sOld = s;
		}
		if (i == MaxIterations) return NoSolution;
	}
	return s;
}

pair<double, double>
AliFmHelix::pathLengths(const AliFmHelix& h, bool scanPeriods) const
{

	//
	//	Cannot handle case where one is a helix
	//  and the other one is a straight line.
	//
	if (mSingularity != h.mSingularity) 
		return pair<double, double>(NoSolution, NoSolution);

	double s1, s2;

	if (mSingularity) {
		//
		//  Analytic solution
		//
		AliFmThreeVector<double> dv = h.mOrigin - mOrigin;
		AliFmThreeVector<double> a(-mCosDipAngle*mSinPhase,
			mCosDipAngle*mCosPhase,
			mSinDipAngle);
		AliFmThreeVector<double> b(-h.mCosDipAngle*h.mSinPhase,
			h.mCosDipAngle*h.mCosPhase,
			h.mSinDipAngle);	
		double ab = a*b;
		double g  = dv*a;
		double k  = dv*b;
		s2 = (k-ab*g)/(ab*ab-1.);
		s1 = g+s2*ab;
		return pair<double, double>(s1, s2);
	}
	else {	
		//
		//  First step: get dca in the xy-plane as start value
		//
		double dx = h.xcenter() - xcenter();
		double dy = h.ycenter() - ycenter();
		double dd = ::sqrt(dx*dx + dy*dy);
		double r1 = 1/curvature();
		double r2 = 1/h.curvature();

		/* malisa 22apr2006 is commenting out the "isnan" check
		 *		if ( !finite(r1) || isnan(r1) || !finite(r2) || isnan(r2) ) {
		 *		 cerr << __FUNCTION__ << " *** error *** ";
		 *			cerr << "  r1=" << r1;
		 *			cerr << "  r2=" << r2 << endl;
		 *			return pair<double, double>(NoSolution, NoSolution);
		 *		}
		 */

		double cosAlpha = (r1*r1 + dd*dd - r2*r2)/(2*r1*dd);

		double s;
		double x, y;
		if (fabs(cosAlpha) < 1) {           // two solutions
			double sinAlpha = sin(acos(cosAlpha));
			x = xcenter() + r1*(cosAlpha*dx - sinAlpha*dy)/dd;
			y = ycenter() + r1*(sinAlpha*dx + cosAlpha*dy)/dd;
			s = pathLength(x, y);
			x = xcenter() + r1*(cosAlpha*dx + sinAlpha*dy)/dd;
			y = ycenter() + r1*(cosAlpha*dy - sinAlpha*dx)/dd;
			double a = pathLength(x, y);
			if (h.distance(at(a),scanPeriods) < h.distance(at(s),scanPeriods)) s = a;
		}
		else {                              // no intersection (or exactly one)
			x = xcenter() + r1*dx/dd;
			y = ycenter() + r1*dy/dd;
			s = pathLength(x, y);
		}

		//
		//   Second step: scan in decreasing intervals around seed 's'
		// 
		const double MinStepSize = 10*micrometer;
		const double MinRange    = 10*centimeter;    
		double dmin              = h.distance(at(s),scanPeriods);
		double range             = max(2*dmin, MinRange);

		/* malisa comments out the "isnan" check 22apr2006
		 *		if ( !finite(range) || isnan(range)) {
		 *			cerr << __FUNCTION__ << " *** error *** ";
		 *			cerr << "  range=" << range << endl;
		 *			return pair<double, double>(NoSolution, NoSolution);
		 *		}
		 */

		double ds = range/10.;
		double slast=-999999, d;
		s1 = s - range/2.;
		s2 = s + range/2.;

		double tmp1;
		double tmp2;

		while (ds > MinStepSize) {
			if ( !(s1<s2) ) {
				cerr << __FUNCTION__ << " *** error *** s1 = ";
				cerr << "    s2 = " << s2;
				cerr << "    range = " << range;
				cerr << "     ds = " << ds << endl;
				return pair<double, double>(NoSolution, NoSolution);
			}	
			tmp1 = (s2-s1);
			tmp2 = tmp1/10.0;
			ds = tmp2;
			if ( ds<0) {
			    cerr << __FUNCTION__ << " *** error *** ds = " << ds << endl;
			    return pair<double, double>(NoSolution, NoSolution);
			}
			int iterations = 0;
			for (double ss=s1; ss<(s2+ds); ss+=ds) {
			    iterations++;
			    if ( iterations > 100 ) {
				cerr << __FUNCTION__ << " *** error *** iterations = " << iterations << endl;
				return pair<double, double>(NoSolution, NoSolution);
			    }
			    d = h.distance(at(ss),scanPeriods);
			    if (d < dmin) {
				dmin = d;
				s = ss;
			    }
			    slast = ss;
			}
			//
			//  In the rare cases where the minimum is at the
			//  the border of the current range we shift the range
			//  and start all over, i.e we do not decrease 'ds'.
			//  Else we decrease the search intervall around the
			//  current minimum and redo the scan in smaller steps.
			//
			if (s == s1) {
				d = 0.8*(s2-s1);
				s1 -= d;
				s2 -= d;
			}
			else if (s == slast) {
				d = 0.8*(s2-s1);
				s1 += d;
				s2 += d;
			}
			else {           
				s1 = s-ds;
				s2 = s+ds;
			//	ds /= 10;
			}
		}
		return pair<double, double>(s, h.pathLength(at(s),scanPeriods));
	}
}


void AliFmHelix::moveOrigin(double s)
{
	if (mSingularity)
		mOrigin	= at(s);
	else {
		AliFmThreeVector<double> newOrigin = at(s);
		double newPhase = atan2(newOrigin.y() - ycenter(),
			newOrigin.x() - xcenter());
		mOrigin = newOrigin;
		setPhase(newPhase);	        
	}
}

int operator== (const AliFmHelix& a, const AliFmHelix& b)
{
	//
	// Checks for numerical identity only !
	//
	return (a.origin()    == b.origin()    &&
		a.dipAngle()  == b.dipAngle()  &&
		a.curvature() == b.curvature() &&
		a.phase()     == b.phase()     &&
		a.h()         == b.h());
}

int operator!= (const AliFmHelix& a, const AliFmHelix& b) {return !(a == b);}

ostream& operator<<(ostream& os, const AliFmHelix& h)
{
	return os << '('
		<< "curvature = "  << h.curvature() << ", " 
		<< "dip angle = "  << h.dipAngle()  << ", "
		<< "phase = "      << h.phase()     << ", "  
		<< "h = "          << h.h()         << ", "    
		<< "origin = "     << h.origin()    << ')';
}





