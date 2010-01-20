#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtBlattWeisskopf.cc,v 1.6 2004/12/21 19:58:41 ryd Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

#include <iostream>
#include <assert.h>
#include <math.h>
#include "EvtGenBase/EvtBlattWeisskopf.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;

EvtBlattWeisskopf::EvtBlattWeisskopf(int LL, double R, double p0)
  : _LL(LL), _radial(R), _p0(p0)
{
  if(R < 0) {

    report(INFO,"EvtGen") << "Radius " << R << " negative" << endl;
    assert(0);
  }

  _radial = R;

  // compute formula for nominal momentum

  _F0 = compute(_p0);
  if(_F0 <= 0) {
    
    report(INFO,"EvtGen") << "Invalid nominal form factor computed " << _F0 << endl;
    assert(0);
  } 
}

EvtBlattWeisskopf::EvtBlattWeisskopf(const EvtBlattWeisskopf& other)
  : _LL(other._LL), _radial(other._radial), _p0(other._p0), _F0(other._F0)
{}

EvtBlattWeisskopf::~EvtBlattWeisskopf()
{}

double EvtBlattWeisskopf::operator()(double p) const
{
  double ret = compute(p)/_F0;
  //  report(INFO,"EvtGen") << p << " " << _p0 << " " << _F0 << " " << _LL << " " << _radial << " " << ret << endl;
  return ret;
}

// Blatt-Weisskopf form factors
// see e.g. hep-ex/0011065
// Dalitz Analysis of the Decay D0->K-pi+pi0 (CLEO)
//
// p   - momentum of either daugher in the meson rest frame,
//       the mass of the meson is used
// pAB - momentum of either daughter in the candidate rest frame
//       the mass of the candidate is used
// R - meson radial parameter
// 
// In the CLEO paper R=5 GeV-1 for D0, R=1.5 for intermediate resonances

double EvtBlattWeisskopf::compute(double p) const
{
  if(p < 0) {
    
    report(INFO,"EvtGen") << "Momentum " << p << " negative in form factor calculation" << endl;
    assert(0);
  }
  else {
    
    double x = p*p*_radial*_radial;
    
    if(0 == _LL) return 1.;
    else
      if(1 == _LL) return sqrt(1.0/(1.0+x));
      else
	if(2 == _LL) return sqrt(1.0/(1.0+x/3.0+x*x/9.0));
	else {
	  report(INFO,"EvtGen") << "Angular momentum " << _LL << " not implemented" << endl;
	  assert(0);
	}
  }
}

