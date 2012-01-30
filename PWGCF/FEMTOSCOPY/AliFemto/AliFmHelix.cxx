///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFmHelix: a helper helix class                                      //
// Includes all the operations and specifications of the helix. Can be   //
// used to determine path lengths, distance of closest approach etc.     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
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
#include "PhysicalConstants.h" 
#include "SystemOfUnits.h"
#ifdef __ROOT__
ClassImpT(AliFmHelix,double);
#endif

#ifdef WIN32
#include "gcc2vs.h"
#endif

#include <iostream>
#include <fstream>
using namespace std;

const double AliFmHelix::fgkNoSolution = 3.e+33;

AliFmHelix::AliFmHelix() :
  fSingularity(0),
  fOrigin(0,0,0),
  fDipAngle(0),
  fCurvature(0),
  fPhase(0),
  fH(0),
  fCosDipAngle(0),
  fSinDipAngle(0),
  fCosPhase(0),
  fSinPhase(0)
{
  //Default constructor
/*noop*/ 
}

AliFmHelix::AliFmHelix(double c, double d, double phase,
		       const AliFmThreeVector<double>& o, int h) :
  fSingularity(0),
  fOrigin(0,0,0),
  fDipAngle(0),
  fCurvature(0),
  fPhase(0),
  fH(0),
  fCosDipAngle(0),
  fSinDipAngle(0),
  fCosPhase(0),
  fSinPhase(0)
{
  // Constructor with helix parameters
  SetParameters(c, d, phase, o, h);
}

AliFmHelix::~AliFmHelix() { 
  // Default destructor
/* noop */ 
}

void AliFmHelix::SetParameters(double c, double dip, double phase,
			       const AliFmThreeVector<double>& o, int h)
{
	//
	//  The order in which the parameters are set is important
	//  since setCurvature might have to adjust the others.
	//
	fH = (h>=0) ? 1 : -1;    // Default is: positive particle
	//             positive field
	fOrigin   = o;
	SetDipAngle(dip);
	SetPhase(phase);

	//
	// Check for singularity and correct for negative curvature.           
	// May change mH and mPhase. Must therefore be set last.
	//
	SetCurvature(c);

	//
	// For the case B=0, h is ill defined. In the following we
	// always assume h = +1. Since phase = psi - h * pi/2
	// we have to correct the phase in case h = -1.
	// This assumes that the user uses the same h for phase
	// as the one he passed to the constructor.
	//
	if (fSingularity && fH == -1) {
		fH = +1;
		SetPhase(fPhase-M_PI);
	}
}

void AliFmHelix::SetCurvature(double val)
{
  // Set helix curvature
	if (val < 0) {
		fCurvature = -val;
		fH = -fH;
		SetPhase(fPhase+M_PI);
	}
	else
		fCurvature = val;

#ifndef ST_NO_NUMERIC_LIMITS
	if (fabs(fCurvature) <= numeric_limits<double>::epsilon())
#else
	if (fabs(fCurvature) <= static_cast<double>(0))
#endif    
		fSingularity = true;			// straight line
	else
		fSingularity = false;            	// curved
}

void AliFmHelix::SetPhase(double val)
{
  // Set helix phase
	fPhase       = val;
	fCosPhase    = cos(fPhase);
	fSinPhase    = sin(fPhase);
	if (fabs(fPhase) > M_PI)
		fPhase = atan2(fSinPhase, fCosPhase);  // force range [-pi,pi]
}

void AliFmHelix::SetDipAngle(double val)
{
  // Set helix dip angle
	fDipAngle    = val;
	fCosDipAngle = cos(fDipAngle);
	fSinDipAngle = sin(fDipAngle);
}

double AliFmHelix::XCenter() const
{
  // Set helix center in X
	if (fSingularity)
		return 0;
	else
		return fOrigin.x()-fCosPhase/fCurvature;
}

double AliFmHelix::YCenter() const
{
  // Set helix center in Y
	if (fSingularity)
		return 0;
	else
		return fOrigin.y()-fSinPhase/fCurvature;
}

double AliFmHelix::FudgePathLength(const AliFmThreeVector<double>& p) const
{
  // Path length
	double s;
	double dx = p.x()-fOrigin.x();
	double dy = p.y()-fOrigin.y();

	if (fSingularity) {
		s = (dy*fCosPhase - dx*fSinPhase)/fCosDipAngle;
	}
	else {
		s = atan2(dy*fCosPhase - dx*fSinPhase,
			1/fCurvature + dx*fCosPhase+dy*fSinPhase)/
			(fH*fCurvature*fCosDipAngle);
	}
	return s;
}

double AliFmHelix::Distance(const AliFmThreeVector<double>& p, bool scanPeriods) const
{
  // calculate distance between origin an p along the helix
	return abs(this->At(PathLength(p,scanPeriods))-p);
}

double AliFmHelix::PathLength(const AliFmThreeVector<double>& p, bool scanPeriods) const 
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
	double dx = p.x()-fOrigin.x();
	double dy = p.y()-fOrigin.y();
	double dz = p.z()-fOrigin.z();

	if (fSingularity) {
		s = fCosDipAngle*(fCosPhase*dy-fSinPhase*dx) +
			fSinDipAngle*dz;
	}
	else { //
#ifndef ST_NO_NAMESPACES
		{
			using namespace units;
#endif
			const double ktMaxPrecisionNeeded = micrometer;
			const int    ktMaxIterations      = 100;

			//
			// The math is taken from Maple with C(expr,optimized) and
			// some hand-editing. It is not very nice but efficient.
			//
			double t34 = fCurvature*fCosDipAngle*fCosDipAngle;
			double t41 = fSinDipAngle*fSinDipAngle;
			double t6, t7, t11, t12, t19;

			//
			// Get a first guess by using the dca in 2D. Since
			// in some extreme cases we might be off by n periods
			// we add (subtract) periods in case we get any closer.
			// 
			s = FudgePathLength(p);

			if (scanPeriods) {
				double ds = Period();
				int    j, jmin = 0;
				double d, dmin = abs(At(s) - p);
				for(j=1; j<ktMaxIterations; j++) {
					if ((d = abs(At(s+j*ds) - p)) < dmin) {
						dmin = d;
						jmin = j;
					}
					else
						break;
				}
				for(j=-1; -j<ktMaxIterations; j--) {
					if ((d = abs(At(s+j*ds) - p)) < dmin) {
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
			// Stops after ktMaxIterations iterations or if the required
			// precision is obtained. Whatever comes first.
			//
			double sOld = s;
			for (int i=0; i<ktMaxIterations; i++) {
				t6  = fPhase+s*fH*fCurvature*fCosDipAngle;
				t7  = cos(t6);
				t11 = dx-(1/fCurvature)*(t7-fCosPhase);
				t12 = sin(t6);
				t19 = dy-(1/fCurvature)*(t12-fSinPhase);
				s  -= (t11*t12*fH*fCosDipAngle-t19*t7*fH*fCosDipAngle -
					(dz-s*fSinDipAngle)*fSinDipAngle)/
					(t12*t12*fCosDipAngle*fCosDipAngle+t11*t7*t34 +
					t7*t7*fCosDipAngle*fCosDipAngle +
					t19*t12*t34+t41);
				if (fabs(sOld-s) < ktMaxPrecisionNeeded) break;
				sOld = s;
			}
#ifndef ST_NO_NAMESPACES
		}
#endif
	}
	return s;
}

double AliFmHelix::Period() const
{
  // period
	if (fSingularity)
#ifndef ST_NO_NUMERIC_LIMITS
		return numeric_limits<double>::max();
#else
		return DBL_MAX;
#endif    
	else	
		return fabs(2*M_PI/(fH*fCurvature*fCosDipAngle)); 
}

pair<double, double> AliFmHelix::PathLength(double r) const
{
	//
	// The math is taken from Maple with C(expr,optimized) and
	// some hand-editing. It is not very nice but efficient.
	// 'first' is the smallest of the two solutions (may be negative)
	// 'second' is the other.
	//
	pair<double,double> tvalue;
	pair<double,double> tVALUE(999999999.,999999999.);
	if (fSingularity) {
		double t1 = fCosDipAngle*(fOrigin.x()*fSinPhase-fOrigin.y()*fCosPhase);
		double t12 = fOrigin.y()*fOrigin.y();
		double t13 = fCosPhase*fCosPhase;
		double t15 = r*r;
		double t16 = fOrigin.x()*fOrigin.x();
		double t20 = -fCosDipAngle*fCosDipAngle*(2.0*fOrigin.x()*fSinPhase*fOrigin.y()*fCosPhase +
			t12-t12*t13-t15+t13*t16);
		if (t20<0.) return tVALUE;
		t20 = ::sqrt(t20);
		tvalue.first  = (t1-t20)/(fCosDipAngle*fCosDipAngle);
		tvalue.second = (t1+t20)/(fCosDipAngle*fCosDipAngle);
	}
	else {
		double t1 = fOrigin.y()*fCurvature;
		double t2 = fSinPhase;
		double t3 = fCurvature*fCurvature;
		double t4 = fOrigin.y()*t2;
		double t5 = fCosPhase;
		double t6 = fOrigin.x()*t5;
		double t8 = fOrigin.x()*fOrigin.x();
		double t11 = fOrigin.y()*fOrigin.y();
		double t14 = r*r;
		double t15 = t14*fCurvature;
		double t17 = t8*t8;
		double t19 = t11*t11;
		double t21 = t11*t3;
		double t23 = t5*t5;
		double t32 = t14*t14;
		double t35 = t14*t3;
		double t38 = 8.0*t4*t6 - 4.0*t1*t2*t8 - 4.0*t11*fCurvature*t6 +
			4.0*t15*t6 + t17*t3 + t19*t3 + 2.0*t21*t8 + 4.0*t8*t23 -
			4.0*t8*fOrigin.x()*fCurvature*t5 - 4.0*t11*t23 -
			4.0*t11*fOrigin.y()*fCurvature*t2 + 4.0*t11 - 4.0*t14 +
			t32*t3 + 4.0*t15*t4 - 2.0*t35*t11 - 2.0*t35*t8;
		double t40 = (-t3*t38);
		if (t40<0.) return tVALUE;
		t40 = ::sqrt(t40);

		double t43 = fOrigin.x()*fCurvature;
		double t45 = 2.0*t5 - t35 + t21 + 2.0 - 2.0*t1*t2 -2.0*t43 - 2.0*t43*t5 + t8*t3;
		double t46 = fH*fCosDipAngle*fCurvature;

		tvalue.first = (-fPhase + 2.0*atan((-2.0*t1 + 2.0*t2 + t40)/t45))/t46;
		tvalue.second = -(fPhase + 2.0*atan((2.0*t1 - 2.0*t2 + t40)/t45))/t46;

		//
		//   Solution can be off by +/- one Period, select smallest
		//
		double p = Period();
		//	malisa deletes "isnan" check 22apr2006	if (!isnan(value.first)) {
		if (fabs(tvalue.first-p) < fabs(tvalue.first)) tvalue.first = tvalue.first-p;
		else if (fabs(tvalue.first+p) < fabs(tvalue.first)) tvalue.first = tvalue.first+p;
		//	malisa	}
		//	malisa deletes "isnan" check 22apr2006		if (!isnan(tvalue.second)) {
		if (fabs(tvalue.second-p) < fabs(tvalue.second)) tvalue.second = tvalue.second-p;
		else if (fabs(tvalue.second+p) < fabs(tvalue.second)) tvalue.second = tvalue.second+p;
		//      malisa }
	}
	if (tvalue.first > tvalue.second)
		swap(tvalue.first,tvalue.second);
	return(tvalue);
}

pair<double, double> AliFmHelix::PathLength(double r, double x, double y, bool /* scanPeriods */)
{
  // path length
	double x0 = fOrigin.x();
	double y0 = fOrigin.y();
	fOrigin.SetX(x0-x);
	fOrigin.SetY(y0-y);
	pair<double, double> result = this->PathLength(r);
	fOrigin.SetX(x0);
	fOrigin.SetY(y0);
	return result;  
}

double AliFmHelix::PathLength(const AliFmThreeVector<double>& r,
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

	if (fSingularity) {
		double t = n.z()*fSinDipAngle +
			n.y()*fCosDipAngle*fCosPhase -
			n.x()*fCosDipAngle*fSinPhase;
		if (t == 0)
			s = fgkNoSolution;
		else
			s = ((r - fOrigin)*n)/t;
	}
	else {
		const double ktMaxPrecisionNeeded = micrometer;
		const int    ktMaxIterations      = 20;

		double tA = fCurvature*((fOrigin - r)*n) -
			n.x()*fCosPhase - 
			n.y()*fSinPhase;
		double t = fH*fCurvature*fCosDipAngle;
		double u = n.z()*fCurvature*fSinDipAngle;

		double a, f, fp;
		double sOld = s = 0;  
		double shiftOld = 0;
		double shift;
		//		(cos(kangMax)-1)/kangMax = 0.1
		const double kangMax = 0.21;
		double deltas = fabs(kangMax/(fCurvature*fCosDipAngle));
		//              dampingFactor = exp(-0.5);
		double dampingFactor = 0.60653;
		int i;

		for (i=0; i<ktMaxIterations; i++) {
			a  = t*s+fPhase;
			double sina = sin(a);
			double cosa = cos(a);
			f  = tA +
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
			if (fabs(sOld-s) < ktMaxPrecisionNeeded) break;
			sOld = s;
		}
		if (i == ktMaxIterations) return fgkNoSolution;
	}
	return s;
}

pair<double, double>
AliFmHelix::PathLengths(const AliFmHelix& h, bool scanPeriods) const
{
	//
	//	Cannot handle case where one is a helix
	//  and the other one is a straight line.
	//
	if (fSingularity != h.fSingularity) 
		return pair<double, double>(fgkNoSolution, fgkNoSolution);

	double s1, s2;

	if (fSingularity) {
		//
		//  Analytic solution
		//
		AliFmThreeVector<double> dv = h.fOrigin - fOrigin;
		AliFmThreeVector<double> a(-fCosDipAngle*fSinPhase,
			fCosDipAngle*fCosPhase,
			fSinDipAngle);
		AliFmThreeVector<double> b(-h.fCosDipAngle*h.fSinPhase,
			h.fCosDipAngle*h.fCosPhase,
			h.fSinDipAngle);	
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
		double dx = h.XCenter() - XCenter();
		double dy = h.YCenter() - YCenter();
		double dd = ::sqrt(dx*dx + dy*dy);
		double r1 = 1/Curvature();
		double r2 = 1/h.Curvature();

		/* malisa 22apr2006 is commenting out the "isnan" check
		 *		if ( !finite(r1) || isnan(r1) || !finite(r2) || isnan(r2) ) {
		 *		 cerr << __FUNCTION__ << " *** error *** ";
		 *			cerr << "  r1=" << r1;
		 *			cerr << "  r2=" << r2 << endl;
		 *			return pair<double, double>(fgkNoSolution, fgkNoSolution);
		 *		}
		 */

		double cosAlpha = (r1*r1 + dd*dd - r2*r2)/(2*r1*dd);

		double s;
		double x, y;
		if (fabs(cosAlpha) < 1) {           // two solutions
			double sinAlpha = sin(acos(cosAlpha));
			x = XCenter() + r1*(cosAlpha*dx - sinAlpha*dy)/dd;
			y = YCenter() + r1*(sinAlpha*dx + cosAlpha*dy)/dd;
			s = PathLength(x, y);
			x = XCenter() + r1*(cosAlpha*dx + sinAlpha*dy)/dd;
			y = YCenter() + r1*(cosAlpha*dy - sinAlpha*dx)/dd;
			double a = PathLength(x, y);
			if (h.Distance(At(a),scanPeriods) < h.Distance(At(s),scanPeriods)) s = a;
		}
		else {                              // no intersection (or exactly one)
			x = XCenter() + r1*dx/dd;
			y = YCenter() + r1*dy/dd;
			s = PathLength(x, y);
		}

		//
		//   Second step: scan in decreasing intervals around seed 's'
		// 
		const double ktMinStepSize = 10*micrometer;
		const double ktMinRange    = 10*centimeter;    
		double dmin              = h.Distance(At(s),scanPeriods);
		double range             = max(2*dmin, ktMinRange);

		/* malisa comments out the "isnan" check 22apr2006
		 *		if ( !finite(range) || isnan(range)) {
		 *			cerr << __FUNCTION__ << " *** error *** ";
		 *			cerr << "  range=" << range << endl;
		 *			return pair<double, double>(fgkNoSolution, fgkNoSolution);
		 *		}
		 */

		double ds = range/10.;
		double slast=-999999, d;
		s1 = s - range/2.;
		s2 = s + range/2.;

		double tmp1;
		double tmp2;

		while (ds > ktMinStepSize) {
			if ( !(s1<s2) ) {
				cerr << __FUNCTION__ << " *** error *** s1 = ";
				cerr << "    s2 = " << s2;
				cerr << "    range = " << range;
				cerr << "     ds = " << ds << endl;
				return pair<double, double>(fgkNoSolution, fgkNoSolution);
			}	
			tmp1 = (s2-s1);
			tmp2 = tmp1/10.0;
			ds = tmp2;
			if ( ds<0) {
			    cerr << __FUNCTION__ << " *** error *** ds = " << ds << endl;
			    return pair<double, double>(fgkNoSolution, fgkNoSolution);
			}
			int iterations = 0;
			for (double ss=s1; ss<(s2+ds); ss+=ds) {
			    iterations++;
			    if ( iterations > 100 ) {
				cerr << __FUNCTION__ << " *** error *** iterations = " << iterations << endl;
				return pair<double, double>(fgkNoSolution, fgkNoSolution);
			    }
			    d = h.Distance(At(ss),scanPeriods);
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
		return pair<double, double>(s, h.PathLength(At(s),scanPeriods));
	}
}


void AliFmHelix::MoveOrigin(double s)
{
  // Move the helix origin
	if (fSingularity)
		fOrigin	= At(s);
	else {
		AliFmThreeVector<double> newOrigin = At(s);
		double newPhase = atan2(newOrigin.y() - YCenter(),
			newOrigin.x() - XCenter());
		fOrigin = newOrigin;
		SetPhase(newPhase);	        
	}
}

int operator== (const AliFmHelix& a, const AliFmHelix& b)
{
	//
	// Checks for numerical identity only !
	//
	return (a.Origin()    == b.Origin()    &&
		a.DipAngle()  == b.DipAngle()  &&
		a.Curvature() == b.Curvature() &&
		a.Phase()     == b.Phase()     &&
		a.H()         == b.H());
}

int operator!= (const AliFmHelix& a, const AliFmHelix& b) {return !(a == b);}

ostream& operator<<(ostream& os, const AliFmHelix& h)
{
	return os << '('
		<< "curvature = "  << h.Curvature() << ", " 
		<< "dip angle = "  << h.DipAngle()  << ", "
		<< "phase = "      << h.Phase()     << ", "  
		<< "h = "          << h.H()         << ", "    
		<< "origin = "     << h.Origin()    << ')';
}





