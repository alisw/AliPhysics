///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFmHelix: a helper helix class                                      //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFMHELIX_H
#define ALIFMHELIX_H

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

#ifdef __SUNPRO_CC
#include <ieeefp.h>
#define __FUNCTION__ "__FILE__:__LINE__"
#endif

class AliFmHelix {
public:
    /// curvature, dip angle, phase, origin, h
    AliFmHelix(double c, double dip, double phase,
	       const AliFmThreeVector<double>& o, int h=-1);
    
    virtual ~AliFmHelix();
    // AliFmHelix(const AliFmHelix&);			// use default
    // AliFmHelix& operator=(const AliFmHelix&);	// use default

    double       DipAngle()   const;           
    double       Curvature()  const;	/// 1/R in xy-plane
    double       Phase()      const;	/// aziumth in xy-plane measured from ring center
    double       XCenter()    const;	/// x-center of circle in xy-plane
    double       YCenter()    const;	/// y-center of circle in xy-plane
    int          H()          const;	/// -sign(q*B);
    
    const AliFmThreeVector<double>& Origin() const;	/// starting point

    void SetParameters(double c, double dip, double phase, const AliFmThreeVector<double>& o, int h);
    
    double       X(double s)  const;
    double       Y(double s)  const;
    double       Z(double s)  const;

    AliFmThreeVector<double>  At(double s) const;

    /// returns period length of helix
    double       Period()       const;
    
    /// path length at given r (cylindrical r)
    pair<double, double> PathLength(double r)   const;
    
    /// path length at given r (cylindrical r, cylinder axis at x,y)
    pair<double, double> PathLength(double r, double x, double y, bool scanPeriods = true);
    
    /// path length at distance of closest approach to a given point
    double       PathLength(const AliFmThreeVector<double>& p, bool scanPeriods = true) const;
    
    /// path length at intersection with plane
    double       PathLength(const AliFmThreeVector<double>& r,
			    const AliFmThreeVector<double>& n) const;

    /// path length at distance of closest approach in the xy-plane to a given point
    double       PathLength(double x, double y) const;

    /// path lengths at dca between two helices 
    pair<double, double> PathLengths(const AliFmHelix& h, bool scanPeriods = true) const;
    
    /// minimal distance between point and helix
    double       Distance(const AliFmThreeVector<double>& p, bool scanPeriods = true) const;    
    
    /// checks for valid parametrization
    bool         Valid(double world = 1.e+5) const {return !Bad(world);}
    int            Bad(double world = 1.e+5) const;
    
    /// move the origin along the helix to s which becomes then s=0
    virtual void MoveOrigin(double s);
    
    static const double fgkNoSolution;  // coinstant indicating lack of solution
    
protected:
    AliFmHelix();
    
    void SetCurvature(double d);	/// performs also various checks   
    void SetPhase(double d);	        
    void SetDipAngle(double d);
    
    /// value of S where distance in x-y plane is minimal
    double FudgePathLength(const AliFmThreeVector<double>& v) const;
    
protected:
    bool                   fSingularity;	// true for straight line case (B=0)
    AliFmThreeVector<double>  fOrigin;          // Helix origin
    double                 fDipAngle;           // Helix dip angle
    double                 fCurvature;          // curvature
    double                 fPhase;              // phase
    int                    fH;			// -sign(q*B);

    double                 fCosDipAngle;        // cosine of the dip angle
    double                 fSinDipAngle;        // sine of the dip angle
    double                 fCosPhase;           // cosine of the phase
    double                 fSinPhase;           // sine of the phase
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
inline int AliFmHelix::H() const {return fH;}

inline double AliFmHelix::DipAngle() const {return fDipAngle;}

inline double AliFmHelix::Curvature() const {return fCurvature;}

inline double AliFmHelix::Phase() const {return fPhase;}

inline double AliFmHelix::X(double s) const
{
    if (fSingularity)
	return fOrigin.x() - s*fCosDipAngle*fSinPhase;
    else
	return fOrigin.x() + (cos(fPhase + s*fH*fCurvature*fCosDipAngle)-fCosPhase)/fCurvature;
}
 
inline double AliFmHelix::Y(double s) const
{
    if (fSingularity)
	return fOrigin.y() + s*fCosDipAngle*fCosPhase;
    else
	return fOrigin.y() + (sin(fPhase + s*fH*fCurvature*fCosDipAngle)-fSinPhase)/fCurvature;
}

inline double AliFmHelix::Z(double s) const
{
    return fOrigin.z() + s*fSinDipAngle;
}

inline const AliFmThreeVector<double>& AliFmHelix::Origin() const {return fOrigin;}

inline AliFmThreeVector<double> AliFmHelix::At(double s) const
{
    return AliFmThreeVector<double>(X(s), Y(s), Z(s));
}

inline double AliFmHelix::PathLength(double x, double y) const
{
    return FudgePathLength(AliFmThreeVector<double>(x, y, 0));
}
inline int AliFmHelix::Bad(double WorldSize) const
{

    int ierr;
    if (!::finite(fDipAngle    )) 	return   11;
    if (!::finite(fCurvature   )) 	return   12;

    ierr = fOrigin.Bad(WorldSize);
    if (ierr)                           return    3+ierr*100;

    if (::fabs(fDipAngle)  >1.58)	return   21;
    double qwe = ::fabs(::fabs(fDipAngle)-M_PI/2);
    if (qwe < 1./WorldSize      ) 	return   31; 

    if (::fabs(fCurvature) > WorldSize)	return   22;
    if (fCurvature < 0          )	return   32;

    if (abs(fH) != 1            )       return   24; 

    return 0;
}

#endif
