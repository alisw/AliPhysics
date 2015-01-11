//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtKine.cc
//
// Description: routines to calculate decay angles. 
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <math.h>
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtdFunction.hh"
#include "EvtGenBase/EvtReport.hh"



double EvtDecayAngle(const EvtVector4R& p,const EvtVector4R& q,
		     const EvtVector4R& d) {

  double pd=p*d;
  double pq=p*q;
  double qd=q*d;
  double mp2=p.mass2();
  double mq2=q.mass2();
  double md2=d.mass2();

  double cost=(pd*mq2-pq*qd)/sqrt((pq*pq-mq2*mp2)*(qd*qd-mq2*md2));

  return cost;

}

double EvtDecayAngleChi(const EvtVector4R& p4_p,const EvtVector4R& p4_d1,
			const EvtVector4R& p4_d2,const EvtVector4R& p4_h1,
			const EvtVector4R& p4_h2 ) {

  EvtVector4R p4_d1p,p4_h1p,p4_h2p,p4_d2p;
  

  // boost all vectors parent restframe
  // This does not boost particle to parent rest frame !!!
  // It goes from parents rest frame to frame where parent has given momentum.
  p4_d1p=boostTo(p4_d1,p4_p,true);
  p4_d2p=boostTo(p4_d2,p4_p,true);
  p4_h1p=boostTo(p4_h1,p4_p,true);
  p4_h2p=boostTo(p4_h2,p4_p,true);
  

  EvtVector4R d1_perp,d1_prime,h1_perp;
  EvtVector4R D;
  
  D=p4_d1p+p4_d2p;

  d1_perp=p4_d1p-(D.dot(p4_d1p)/D.dot(D))*D;
  h1_perp=p4_h1p-(D.dot(p4_h1p)/D.dot(D))*D;
  
  // orthogonal to both D and d1_perp
  
  d1_prime=D.cross(d1_perp);
  
  d1_perp= d1_perp/d1_perp.d3mag();
  d1_prime= d1_prime/d1_prime.d3mag();
  
  double x,y;
  
  x=d1_perp.dot(h1_perp);
  y=d1_prime.dot(h1_perp);
  
  double chi=atan2(y,x);
  
  if (chi<0.0) chi+=EvtConst::twoPi;
  
  return chi;
  
}



double EvtDecayPlaneNormalAngle(const EvtVector4R& p,const EvtVector4R& q,
                          const EvtVector4R& d1,const EvtVector4R& d2){

  EvtVector4C lc=dual(EvtGenFunctions::directProd(d1,d2)).cont2(q);

  EvtVector4R l(real(lc.get(0)),real(lc.get(1)),
		real(lc.get(2)),real(lc.get(3)));

  double pq=p*q;

  return q.mass()*(p*l)/sqrt(-(pq*pq-p.mass2()*q.mass2())*l.mass2());


}


// Calculate phi using the given 4 vectors (all in the same frame)
double EvtDecayAnglePhi( const EvtVector4R& z, const EvtVector4R& p, const
        EvtVector4R& q, const EvtVector4R& d )
{
    double eq = (p * q) / p.mass();
    double ed = (p * d) / p.mass();
    double mq = q.mass();
    double q2 = p.mag2r3(q);
    double qd = p.dotr3(q,d);
    double zq = p.dotr3(z,q);
    double zd = p.dotr3(z,d);
    double alpha = (eq - mq)/(q2 * mq) * qd - ed/mq;

    double y = p.scalartripler3(z,q,d) + alpha * p.scalartripler3(z,q,q);
    double x = (zq * (qd + alpha * q2) - q2 * (zd + alpha * zq)) / sqrt(q2);

    double phi = atan2(y,x);

    return phi<0 ? (phi+EvtConst::twoPi) : phi;
}

EvtComplex wignerD( int j, int m1, int m2, double phi, 
		    double theta, double gamma )
{

    EvtComplex gp(0.0, -phi*m1);
    EvtComplex gm(0.0, -gamma*m2);

    return exp( gp ) * EvtdFunction::d(j, m1, m2, theta) * exp( gm );
}


