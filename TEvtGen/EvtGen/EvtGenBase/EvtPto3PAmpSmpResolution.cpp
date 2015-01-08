#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *  Author: D. Dujmic, ddujmic@slac.stanford.edu
 *
 * Copyright (C) 2005 SLAC
 *
 *******************************************************************************/

#include <assert.h>
#include <math.h>
#include <iostream>
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPto3PAmpSmpResolution.hh"
#include "EvtGenBase/EvtPto3PAmp.hh"
#include "EvtGenBase/EvtDalitzCoord.hh"
#include "EvtGenBase/EvtCyclic3.hh"
using std::cout;
using std::endl;
using EvtCyclic3::Index;
using EvtCyclic3::Pair;



EvtPto3PAmpSmpResolution::EvtPto3PAmpSmpResolution(EvtDalitzPlot dp, Pair pairAng, Pair pairRes, 
						   EvtSpinType::spintype spin, 
						   const EvtPropagator& prop, NumType typeN) 
  : EvtPto3PAmp(dp, pairAng, pairRes, spin, prop, typeN)
{}



EvtPto3PAmpSmpResolution::EvtPto3PAmpSmpResolution(const EvtPto3PAmp& other) 
  : EvtPto3PAmp(other)
{}


EvtPto3PAmpSmpResolution::~EvtPto3PAmpSmpResolution()
{}


EvtComplex 
EvtPto3PAmpSmpResolution::evalPropagator(double m) const
{
  EvtComplex prop(0,0);
  
  if (_sigma>0) { // convolved
    int nconv=20;
    double min=m+_bias-_sigma*2.5;
    double max=m+_bias+_sigma*2.5;
    double dm=(max-min)/nconv;
    static double  sqrt2pi = sqrt(2*3.14159);
    double ifact = 1./(sqrt2pi*_sigma);
    for (int i=0;i<nconv; i++) {
      double mprime = min+dm*(i+0.5);
      double t = (mprime-m)/_sigma;
      prop += ifact * exp(-0.5*t*t)*EvtPto3PAmp::evalPropagator(m)*dm;
    }
  } else {
    prop = EvtPto3PAmp::evalPropagator(m);
  }

  return prop;
}
