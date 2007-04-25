/***************************************************************************
 *
 * $Id$
 *
 * Author: Brian Lasiuk, Sep 1997
 ***************************************************************************
 *
 * Description: 
 * Parametrization of a physical helix.
 * 
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.7  2005/07/06 18:49:56  fisyak
 * Replace StHelixD, StLorentzVectorD,StLorentzVectorF,StMatrixD,StMatrixF,AliFmPhysicalHelixD,StThreeVectorD,StThreeVectorF by templated version
 *
 * Revision 1.6  2003/09/02 17:59:35  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.5  2002/06/21 17:49:26  genevb
 * Some minor speed improvements
 *
 * Revision 1.4  2002/02/20 00:56:23  ullrich
 * Added methods to calculate signed DCA.
 *
 * Revision 1.3  1999/02/24 11:42:18  ullrich
 * Fixed bug in momentum().
 *
 * Revision 1.2  1999/02/12 01:01:04  wenaus
 * Fix bug in momentum calculation
 *
 * Revision 1.1  1999/01/30 03:59:04  fisyak
 * Root Version of StarClassLibrary
 *
 * Revision 1.1  1999/01/23 00:29:21  ullrich
 * Initial Revision
 *
 **************************************************************************/
#include <math.h>
#include "AliFmHelix.h"
#include "AliFmPhysicalHelix.h"
#include "PhysicalConstants.h" 
#include "SystemOfUnits.h"
#ifdef __ROOT__
ClassImpT(AliFmPhysicalHelix,double);
#endif
AliFmPhysicalHelix::AliFmPhysicalHelix(){}

AliFmPhysicalHelix::~AliFmPhysicalHelix() { /* nop */ }

AliFmPhysicalHelix::AliFmPhysicalHelix(const AliFmThreeVector<double>& p,
				       const AliFmThreeVector<double>& o,
				       double B, double q)
{
  mH = (q*B <= 0) ? 1 : -1;
  if(p.y() == 0 && p.x() == 0)
    setPhase((M_PI/4)*(1-2.*mH));
  else
    setPhase(atan2(p.y(),p.x())-mH*M_PI/2);
  setDipAngle(atan2(p.z(),p.perp()));
  mOrigin = o;
  
#ifndef ST_NO_NAMESPACES
  {
    using namespace units;
#endif
    setCurvature(fabs((c_light*nanosecond/meter*q*B/tesla)/
		      (abs(p.mag())/GeV*mCosDipAngle)/meter));   
#ifndef ST_NO_NAMESPACES
  }
#endif
}

AliFmPhysicalHelix::AliFmPhysicalHelix(double c, double d, double phase,
				 const AliFmThreeVector<double>& o, int h)
    : AliFmHelix(c, d, phase, o, h) { /* nop */}


AliFmThreeVector<double> AliFmPhysicalHelix::momentum(double B) const
{
    if (mSingularity)
	return(AliFmThreeVector<double>(0,0,0));
    else {
#ifndef ST_NO_NAMESPACES
	{
	    using namespace units;
#endif
	    double pt = GeV*fabs(c_light*nanosecond/meter*B/tesla)/(fabs(mCurvature)*meter);
    
	    return (AliFmThreeVector<double>(pt*cos(mPhase+mH*M_PI/2),   // pos part pos field
					  pt*sin(mPhase+mH*M_PI/2),
					  pt*tan(mDipAngle)));
#ifndef ST_NO_NAMESPACES
	}
#endif
    }
}

AliFmThreeVector<double> AliFmPhysicalHelix::momentumAt(double S, double B) const
{
    // Obtain phase-shifted momentum from phase-shift of origin
    double xc = this->xcenter();
    double yc = this->ycenter();
    double rx = (y(S)-yc)/(mOrigin.y()-yc);
    double ry = (x(S)-xc)/(mOrigin.x()-xc);
    return (this->momentum(B)).pseudoProduct(rx,ry,1.0);
}

int AliFmPhysicalHelix::charge(double B) const
{
    return (B > 0 ? -mH : mH);
}

double AliFmPhysicalHelix::geometricSignedDistance(double x, double y)  
{
    // Geometric signed distance
    double thePath = this->pathLength(x,y);
    AliFmThreeVector<double> DCA2dPosition = this->at(thePath);
    DCA2dPosition.setZ(0);
    AliFmThreeVector<double> position(x,y,0);
    AliFmThreeVector<double> DCAVec = (DCA2dPosition-position);
    AliFmThreeVector<double> momVec;
    // Deal with straight tracks
    if (this->mSingularity) {
	momVec = this->at(1)- this->at(0);
	momVec.setZ(0);
    }
    else {
	momVec = this->momentumAt(thePath,1./tesla); // Don't care about Bmag.  Helicity is what matters.
	momVec.setZ(0);
    }
    
    double cross = DCAVec.x()*momVec.y() - DCAVec.y()*momVec.x();
    double theSign = (cross>=0) ? 1. : -1.;
    return theSign*DCAVec.perp();
}

double AliFmPhysicalHelix::curvatureSignedDistance(double x, double y) 
{
    // Protect against mH = 0 or zero field
    if (this->mSingularity || abs(this->mH)<=0) {
	return (this->geometricSignedDistance(x,y));
    }
    else {
	return (this->geometricSignedDistance(x,y))/(this->mH);
    }
    
}

double AliFmPhysicalHelix::geometricSignedDistance(const AliFmThreeVector<double>& pos) 
{
    double sdca2d = this->geometricSignedDistance(pos.x(),pos.y());
    double theSign = (sdca2d>=0) ? 1. : -1.;
    return (this->distance(pos))*theSign;
}

double AliFmPhysicalHelix::curvatureSignedDistance(const AliFmThreeVector<double>& pos) 
{
    double sdca2d = this->curvatureSignedDistance(pos.x(),pos.y());
    double theSign = (sdca2d>=0) ? 1. : -1.;
    return (this->distance(pos))*theSign;
}
