///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFmHelix: a helper helix class                                      //
// Includes all the operations and specifications of the helix. Can be   //
// used to determine path lengths, distance of closest approach etc.     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
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
  // Constructor from given parameters
  fH = (q*B <= 0) ? 1 : -1;
  if(p.y() == 0 && p.x() == 0)
    SetPhase((M_PI/4)*(1-2.*fH));
  else
    SetPhase(atan2(p.y(),p.x())-fH*M_PI/2);
  SetDipAngle(atan2(p.z(),p.perp()));
  fOrigin = o;
  
#ifndef ST_NO_NAMESPACES
  {
    using namespace units;
#endif
    SetCurvature(fabs((c_light*nanosecond/meter*q*B/tesla)/
		      (abs(p.mag())/GeV*fCosDipAngle)/meter));   
#ifndef ST_NO_NAMESPACES
  }
#endif
}

AliFmPhysicalHelix::AliFmPhysicalHelix(double c, double d, double phase,
				       const AliFmThreeVector<double>& o, int h)
  : AliFmHelix(c, d, phase, o, h) { /* nop */}


AliFmThreeVector<double> AliFmPhysicalHelix::Momentum(double B) const
{
  // momentum for given magnetic field
  if (fSingularity)
    return(AliFmThreeVector<double>(0,0,0));
  else {
#ifndef ST_NO_NAMESPACES
    {
	    using namespace units;
#endif
	    double pt = GeV*fabs(c_light*nanosecond/meter*B/tesla)/(fabs(fCurvature)*meter);
	    
	    return (AliFmThreeVector<double>(pt*cos(fPhase+fH*M_PI/2),   // pos part pos field
					  pt*sin(fPhase+fH*M_PI/2),
					  pt*tan(fDipAngle)));
#ifndef ST_NO_NAMESPACES
	}
#endif
    }
}

AliFmThreeVector<double> AliFmPhysicalHelix::MomentumAt(double S, double B) const
{
    // Obtain phase-shifted momentum from phase-shift of origin
    double xc = this->XCenter();
    double yc = this->YCenter();
    double rx = (Y(S)-yc)/(fOrigin.y()-yc);
    double ry = (X(S)-xc)/(fOrigin.x()-xc);
    return (this->Momentum(B)).pseudoProduct(rx,ry,1.0);
}

int AliFmPhysicalHelix::Charge(double B) const
{
  // charge
    return (B > 0 ? -fH : fH);
}

double AliFmPhysicalHelix::GeometricSignedDistance(double x, double y)  
{
    // Geometric signed distance
    double thePath = this->PathLength(x,y);
    AliFmThreeVector<double> tDCA2dPosition = this->At(thePath);
    tDCA2dPosition.SetZ(0);
    AliFmThreeVector<double> position(x,y,0);
    AliFmThreeVector<double> tDCAVec = (tDCA2dPosition-position);
    AliFmThreeVector<double> momVec;
    // Deal with straight tracks
    if (this->fSingularity) {
	momVec = this->At(1)- this->At(0);
	momVec.SetZ(0);
    }
    else {
	momVec = this->MomentumAt(thePath,1./tesla); // Don't care about Bmag.  Helicity is what matters.
	momVec.SetZ(0);
    }
    
    double cross = tDCAVec.x()*momVec.y() - tDCAVec.y()*momVec.x();
    double theSign = (cross>=0) ? 1. : -1.;
    return theSign*tDCAVec.perp();
}

double AliFmPhysicalHelix::CurvatureSignedDistance(double x, double y) 
{
    // Protect against fH = 0 or zero field
    if (this->fSingularity || abs(this->fH)<=0) {
	return (this->GeometricSignedDistance(x,y));
    }
    else {
	return (this->GeometricSignedDistance(x,y))/(this->fH);
    }
    
}

double AliFmPhysicalHelix::GeometricSignedDistance(const AliFmThreeVector<double>& pos) 
{
  // Geometric distance
    double sdca2d = this->GeometricSignedDistance(pos.x(),pos.y());
    double theSign = (sdca2d>=0) ? 1. : -1.;
    return (this->Distance(pos))*theSign;
}

double AliFmPhysicalHelix::CurvatureSignedDistance(const AliFmThreeVector<double>& pos) 
{
  // Distance with the sign dependent on curvature sigm
    double sdca2d = this->CurvatureSignedDistance(pos.x(),pos.y());
    double theSign = (sdca2d>=0) ? 1. : -1.;
    return (this->Distance(pos))*theSign;
}
