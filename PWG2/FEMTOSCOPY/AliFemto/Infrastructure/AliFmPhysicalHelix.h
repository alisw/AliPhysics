/***************************************************************************
 *
 * $Id$
 *
 * Author: Brian Lasiuk, Sep 1997
 ***************************************************************************
 *
 * Description: 
 * Parametrization of a physical helix. See the SCL user guide for more.
 * 
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.4  2005/07/06 18:49:56  fisyak
 * Replace AliFmHelixD, AliFmLorentzVectorD,AliFmLorentzVectorF,AliFmMatrixD,AliFmMatrixF,AliFmPhysicalHelixD,AliFmThreeVectorD,AliFmThreeVectorF by templated version
 *
 * Revision 1.3  2002/06/21 17:49:26  genevb
 * Some minor speed improvements
 *
 * Revision 1.2  2002/02/20 00:56:23  ullrich
 * Added methods to calculate signed DCA.
 *
 * Revision 1.1  1999/01/30 03:59:04  fisyak
 * Root Version of AliFmarClassLibrary
 *
 * Revision 1.1  1999/01/23 00:27:59  ullrich
 * Initial Revision
 *
 **************************************************************************/
#ifndef ST_PHYSICAL_HELIX_HH
#define ST_PHYSICAL_HELIX_HH

#include "AliFmThreeVector.h"
#include "AliFmHelix.h"

class AliFmPhysicalHelix : public AliFmHelix {
public:
    // Requires: momentum, origin, signed Magnetic Field
    //           and Charge of particle (+/- 1)
    AliFmPhysicalHelix(const AliFmThreeVector<double>&,
		    const AliFmThreeVector<double>&,
		    double, double);
    
    // curvature, dip angle, phase, origin, h
    AliFmPhysicalHelix(double, double, double,
		    const AliFmThreeVector<double>&, int h=-1);
    AliFmPhysicalHelix();
    
    ~AliFmPhysicalHelix();

    // Requires:  signed Magnetic Field
    AliFmThreeVector<double> momentum(double) const;     // returns the momentum at origin
    AliFmThreeVector<double> momentumAt(double, double) const; // returns momemtum at S
    int                   charge(double)   const;     // returns charge of particle
    // 2d DCA to x,y point signed relative to curvature
    double curvatureSignedDistance(double x, double y) ;
    // 2d DCA to x,y point signed relative to rotation 
    double geometricSignedDistance(double x, double y) ;
    // 3d DCA to 3d point signed relative to curvature
    double curvatureSignedDistance(const AliFmThreeVector<double>&) ;
    // 3d DCA to 3d point signed relative to rotation
    double geometricSignedDistance(const AliFmThreeVector<double>&) ;
    
#ifdef __ROOT__
  ClassDef(AliFmPhysicalHelix,1)
#endif
};

#endif
