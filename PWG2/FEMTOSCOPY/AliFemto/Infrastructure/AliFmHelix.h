/**
 * \class AliFmHelix
 * \author Thomas Ullrich, Sep 26 1997
 * 
 * Parametrization of a helix. Can also cope with straight tracks, i.e.
 * with zero curvature. This represents only the mathematical model of 
 * a helix. See the SCL user guide for more. 
 */
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
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.11  2005/10/13 22:23:27  genevb
 * NoSolution is public
 *
 * Revision 1.10  2005/07/06 18:49:56  fisyak
 * Replace AliFmHelixD, AliFmLorentzVectorD,AliFmLorentzVectorF,AliFmMatrixD,AliFmMatrixF,AliFmPhysicalHelixD,AliFmThreeVectorD,AliFmThreeVectorF by templated version
 *
 * Revision 1.9  2004/12/02 02:51:16  ullrich
 * Added option to pathLenghth() and distance() to search for
 * DCA only within one period. Default stays as it was.
 *
 * Revision 1.8  2003/10/30 20:06:46  perev
 * Check of quality added
 *
 * Revision 1.7  2002/06/21 17:49:25  genevb
 * Some minor speed improvements
 *
 * Revision 1.6  2002/04/24 02:41:55  ullrich
 * Restored old format.
 *
 **************************************************************************/

#ifndef ST_HELIX_HH
#define ST_HELIX_HH

#include <math.h>
#include <utility>
#include <algorithm>
#include "AliFmThreeVector.h"
#if !defined(ST_NO_NAMESPACES)
using std::pair;
using std::swap;
using std::max;
#endif

#ifdef WIN32
#include "gcc2vs.h"
#endif

class AliFmHelix {
public:
    /// curvature, dip angle, phase, origin, h
    AliFmHelix(double c, double dip, double phase,
	    const AliFmThreeVector<double>& o, int h=-1);
    
    virtual ~AliFmHelix();
    // AliFmHelix(const AliFmHelix&);			// use default
    // AliFmHelix& operator=(const AliFmHelix&);	// use default

    double       dipAngle()   const;           
    double       curvature()  const;	/// 1/R in xy-plane
    double       phase()      const;	/// aziumth in xy-plane measured from ring center
    double       xcenter()    const;	/// x-center of circle in xy-plane
    double       ycenter()    const;	/// y-center of circle in xy-plane
    int          h()          const;	/// -sign(q*B);
    
    const AliFmThreeVector<double>& origin() const;	/// starting point

    void setParameters(double c, double dip, double phase, const AliFmThreeVector<double>& o, int h);
    
    double       x(double s)  const;
    double       y(double s)  const;
    double       z(double s)  const;

    AliFmThreeVector<double>  at(double s) const;

    /// returns period length of helix
    double       period()       const;
    
    /// path length at given r (cylindrical r)
    pair<double, double> pathLength(double r)   const;
    
    /// path length at given r (cylindrical r, cylinder axis at x,y)
    pair<double, double> pathLength(double r, double x, double y, bool scanPeriods = true);
    
    /// path length at distance of closest approach to a given point
    double       pathLength(const AliFmThreeVector<double>& p, bool scanPeriods = true) const;
    
    /// path length at intersection with plane
    double       pathLength(const AliFmThreeVector<double>& r,
			    const AliFmThreeVector<double>& n) const;

    /// path length at distance of closest approach in the xy-plane to a given point
    double       pathLength(double x, double y) const;

    /// path lengths at dca between two helices 
    pair<double, double> pathLengths(const AliFmHelix&, bool scanPeriods = true) const;
    
    /// minimal distance between point and helix
    double       distance(const AliFmThreeVector<double>& p, bool scanPeriods = true) const;    
    
    /// checks for valid parametrization
    bool         valid(double world = 1.e+5) const {return !bad(world);}
    int            bad(double world = 1.e+5) const;
    
    /// move the origin along the helix to s which becomes then s=0
    virtual void moveOrigin(double s);
    
    static const double NoSolution;
    
protected:
    AliFmHelix();
    
    void setCurvature(double);	/// performs also various checks   
    void setPhase(double);	        
    void setDipAngle(double);
    
    /// value of S where distance in x-y plane is minimal
    double fudgePathLength(const AliFmThreeVector<double>&) const;
    
protected:
    bool                   mSingularity;	// true for straight line case (B=0)
    AliFmThreeVector<double>  mOrigin;
    double                 mDipAngle;
    double                 mCurvature;
    double                 mPhase;
    int                    mH;			// -sign(q*B);

    double                 mCosDipAngle;
    double                 mSinDipAngle;
    double                 mCosPhase;
    double                 mSinPhase;
#ifdef __ROOT__
  ClassDef(AliFmHelix,1)
#endif
};

//
//     Non-member functions
//
int operator== (const AliFmHelix&, const AliFmHelix&);
int operator!= (const AliFmHelix&, const AliFmHelix&);
ostream& operator<<(ostream&, const AliFmHelix&);

//
//     Inline functions
//
inline int AliFmHelix::h() const {return mH;}

inline double AliFmHelix::dipAngle() const {return mDipAngle;}

inline double AliFmHelix::curvature() const {return mCurvature;}

inline double AliFmHelix::phase() const {return mPhase;}

inline double AliFmHelix::x(double s) const
{
    if (mSingularity)
	return mOrigin.x() - s*mCosDipAngle*mSinPhase;
    else
	return mOrigin.x() + (cos(mPhase + s*mH*mCurvature*mCosDipAngle)-mCosPhase)/mCurvature;
}
 
inline double AliFmHelix::y(double s) const
{
    if (mSingularity)
	return mOrigin.y() + s*mCosDipAngle*mCosPhase;
    else
	return mOrigin.y() + (sin(mPhase + s*mH*mCurvature*mCosDipAngle)-mSinPhase)/mCurvature;
}

inline double AliFmHelix::z(double s) const
{
    return mOrigin.z() + s*mSinDipAngle;
}

inline const AliFmThreeVector<double>& AliFmHelix::origin() const {return mOrigin;}

inline AliFmThreeVector<double> AliFmHelix::at(double s) const
{
    return AliFmThreeVector<double>(x(s), y(s), z(s));
}

inline double AliFmHelix::pathLength(double x, double y) const
{
    return fudgePathLength(AliFmThreeVector<double>(x, y, 0));
}
inline int AliFmHelix::bad(double WorldSize) const
{

    int ierr;
    if (!::finite(mDipAngle    )) 	return   11;
    if (!::finite(mCurvature   )) 	return   12;

    ierr = mOrigin.bad(WorldSize);
    if (ierr)                           return    3+ierr*100;

    if (::fabs(mDipAngle)  >1.58)	return   21;
    double qwe = ::fabs(::fabs(mDipAngle)-M_PI/2);
    if (qwe < 1./WorldSize      ) 	return   31; 

    if (::fabs(mCurvature) > WorldSize)	return   22;
    if (mCurvature < 0          )	return   32;

    if (abs(mH) != 1            )       return   24; 

    return 0;
}

#endif
