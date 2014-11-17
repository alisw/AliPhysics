/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtDalitzPoint.hh,v 1.2 2009-03-16 16:44:53 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// This class describes the complete kinematics of the Dalitz decay.
// It holds all the six invariant momentum products, three daughter
// particle masses and three invariant masses of pairs of particles.
// This description is completely symmetric with respect to particle 
// permutations.
//
// Another way to slice the six coordinate is to make a transformation 
// to the mass of the decaying particle. The four masses make up a 
// Dalitz plot. The other two are coordinates of a point in the plot.

#ifndef EVT_DALITZ_POINT_HH
#define EVT_DALITZ_POINT_HH

#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzCoord.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"

class EvtDalitzPoint {

public:

  EvtDalitzPoint();
  EvtDalitzPoint(double mA, double mB, double mC, 
		 double qAB, double qBC, double qCA);
  EvtDalitzPoint(double mA, double mB, double mC, 
		 EvtCyclic3::Pair i, double qres, double qhel, double qsum);
  EvtDalitzPoint(const EvtDalitzPlot&, const EvtDalitzCoord&);
  EvtDalitzPoint(const EvtDalitzPoint& other);
  ~EvtDalitzPoint();

  EvtDalitzCoord getDalitzPoint(EvtCyclic3::Pair i, EvtCyclic3::Pair j) const;
  EvtDalitzPlot  getDalitzPlot() const;

  double q(EvtCyclic3::Pair) const;
  double bigM() const;
  double m(EvtCyclic3::Index) const;

  // Zemach variables

  double qres(EvtCyclic3::Pair i) const;
  double qhel(EvtCyclic3::Pair i) const;
  double qsum() const;

  // Kinematic quantities
  //
  // pp  - 4 momentum product
  // e,p,cosTh - energy/moementum in rest-frame of j


  double qMin(EvtCyclic3::Pair i, EvtCyclic3::Pair j) const;
  double qMax(EvtCyclic3::Pair i, EvtCyclic3::Pair j) const;
  double pp(EvtCyclic3::Index i, EvtCyclic3::Index j) const; 
  double e(EvtCyclic3::Index i, EvtCyclic3::Pair j) const;
  double p(EvtCyclic3::Index i, EvtCyclic3::Pair j) const;
  double cosTh(EvtCyclic3::Pair pairAng, EvtCyclic3::Pair pairRes) const;

  bool isValid() const;

  void print() const;

private:

  double _mA, _mB, _mC;     // masses
  double _qAB, _qBC, _qCA;  // masses squared
};

#endif
