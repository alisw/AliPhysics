#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtBlattWeisskopf.cpp,v 1.3 2009-03-16 15:56:37 robbep Exp $
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

    report(Severity::Info,"EvtGen") << "Radius " << R << " negative" << endl;
    assert(0);
  }

  _radial = R;

  // compute formula for nominal momentum

  _F0 = compute(_p0);
  if(_F0 <= 0) {
    
    report(Severity::Info,"EvtGen") << "Invalid nominal form factor computed " << _F0 << endl;
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
  //  report(Severity::Info,"EvtGen") << p << " " << _p0 << " " << _F0 << " " << _LL << " " << _radial << " " << ret << endl;
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

  double value(1.0);

  double z = p*_radial;
  double zSq = z*z;

  if (_LL == 0) {
    value = 1.0;
  } else if (_LL == 1) {
    value = sqrt(1.0/(1.0 + zSq));
  } else if (_LL == 2) {
    value = sqrt(1.0/(zSq*(zSq + 3.0) + 9.0));
  } else if (_LL == 3) {
    double denom = zSq*(zSq*(zSq + 6.0) + 45.0) + 225.0;
    value = sqrt(1.0/denom);
  } else if (_LL == 4) {
    double denom = zSq*(zSq*(zSq*(zSq + 10.0) + 135.0) + 1575.0) + 11025.0;
    value = sqrt(1.0/denom);
  } else if (_LL == 5) {
    double denom = zSq*(zSq*(zSq*(zSq*(zSq + 15.0) + 315.0) + 6300.0) + 99225.0) + 893025.0;
    value = sqrt(1.0/denom);
  }

  return value;

}
