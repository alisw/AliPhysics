///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFmHelix: a helper helix class                                      //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#ifndef ALIFMPHYSICALHELIX_H
#define ALIFMPHYSICALHELIX_H

#include "AliFmThreeVector.h"
#include "AliFmHelix.h"

class AliFmPhysicalHelix : public AliFmHelix {
 public:
  // Requires: momentum, origin, signed Magnetic Field
  //           and Charge of particle (+/- 1)
  AliFmPhysicalHelix(const AliFmThreeVector<double>& v,
		     const AliFmThreeVector<double>& v,
		     double x, double x);
    
  // curvature, dip angle, phase, origin, h
  AliFmPhysicalHelix(double curvature, double dipAngle, double phase,
		     const AliFmThreeVector<double>& origin, int h=-1);
  AliFmPhysicalHelix();
  
  ~AliFmPhysicalHelix();
  
  // Requires:  signed Magnetic Field
  AliFmThreeVector<double> Momentum(double x) const;     // returns the momentum at origin
  AliFmThreeVector<double> MomentumAt(double x, double x) const; // returns momemtum at S
  int                   Charge(double x)   const;     // returns charge of particle
  // 2d DCA to x,y point signed relative to curvature
  double CurvatureSignedDistance(double x, double y) ;
  // 2d DCA to x,y point signed relative to rotation 
  double GeometricSignedDistance(double x, double y) ;
  // 3d DCA to 3d point signed relative to curvature
  double CurvatureSignedDistance(const AliFmThreeVector<double>& v) ;
  // 3d DCA to 3d point signed relative to rotation
  double GeometricSignedDistance(const AliFmThreeVector<double>& v) ;
  
#ifdef __ROOT__
  ClassDef(AliFmPhysicalHelix,1)
#endif
};

#endif
