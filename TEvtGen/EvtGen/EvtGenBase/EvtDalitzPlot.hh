//-----------------------------------------------------------------------
// File and Version Information: 
//      $Id: EvtDalitzPlot.hh,v 1.2 2009-03-16 16:44:53 robbep Exp $
// 
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1998 Caltech, UCSB
//
// Module creator:
//      Alexei Dvoretskii, Caltech, 2001-2002.
//-----------------------------------------------------------------------

#ifndef EVT_DALITZ_PLOT_HH
#define EVT_DALITZ_PLOT_HH

#include <assert.h>
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtTwoBodyVertex.hh"
#include "EvtGenBase/EvtDecayMode.hh"

class EvtDalitzPlot {
public:

  EvtDalitzPlot();
  EvtDalitzPlot(double mA, double mB, double mC, double bigM, double ldel = 0., double rdel = 0.);
  EvtDalitzPlot(const EvtDecayMode& mode, double ldel = 0., double rdel = 0.);
  EvtDalitzPlot(const EvtDalitzPlot& other);
  ~EvtDalitzPlot();
  bool operator==(const EvtDalitzPlot& other) const;  
  const EvtDalitzPlot* clone() const;

  
  // Absolute limits for masses squared in the Dalitz plot
  // e.g. qAbsMin(0) is the lowest possible value
  // for m2 of particles {12}
  
  double qAbsMin(EvtCyclic3::Pair i) const;
  double qAbsMax(EvtCyclic3::Pair i) const;
  double mAbsMin(EvtCyclic3::Pair i) const;
  double mAbsMax(EvtCyclic3::Pair i) const;

  // Absolute limits for Zemach coordinate qres and qhel (approximate)
  // qHelAbsMin(BC,CA) means absolute minimum for (qCA-qAB)/2.

  double qResAbsMin(EvtCyclic3::Pair i) const;
  double qResAbsMax(EvtCyclic3::Pair i) const;
  double qHelAbsMin(EvtCyclic3::Pair i) const;
  double qHelAbsMax(EvtCyclic3::Pair i) const;
  inline double qSumMin() const { return sum() + _ldel; }
  inline double qSumMax() const { return sum() + _rdel; }
  inline bool fuzzy() const { return (_rdel - _ldel != 0.); }

  // Find the area of the Dalitz plot by numeric integration. (N bins for variable q(i) are used).
  // Very large numbers of N can result in a very long calculation. It should not 
  // matter which two pairs f variables are used. The integral should eventually 
  // converge to the same number

  double getArea(int N = 1000, EvtCyclic3::Pair i = EvtCyclic3::AB, EvtCyclic3::Pair j = EvtCyclic3::BC) const;

  // Limits for masses squared when one mass squared is known

  double qMin(EvtCyclic3::Pair i, EvtCyclic3::Pair j, double q) const;
  double qMax(EvtCyclic3::Pair i, EvtCyclic3::Pair j, double q) const;


  // Coordinate transformations

  double cosTh(EvtCyclic3::Pair i1, double q1, EvtCyclic3::Pair i2, double q2) const;
  double e(EvtCyclic3::Index i, EvtCyclic3::Pair j, double q) const;
  double p(EvtCyclic3::Index i, EvtCyclic3::Pair j, double q) const;

  double q(EvtCyclic3::Pair i1, double cosTh, EvtCyclic3::Pair i2, double q2) const;

  // |J| of transformation of qi to cosTh in the rest-frame of j

  double jacobian(EvtCyclic3::Pair i, double q) const;


  // Given resonance index and mass returns decay 
  // and birth vertices

  EvtTwoBodyVertex vD(EvtCyclic3::Pair iRes, double m0, int L) const;
  EvtTwoBodyVertex vB(EvtCyclic3::Pair iRes, double m0, int L) const;

  // Accessors

  double sum() const;
  inline double bigM() const { return _bigM; }
  inline double mA() const { return _mA; } 
  inline double mB() const { return _mB; } 
  inline double mC() const { return _mC; } 
  double m(EvtCyclic3::Index i) const;


  void print() const;

  void sanityCheck() const;

protected:

  // Defines two dimensional dalitz plot

  double _mA;
  double _mB;
  double _mC;
  double _bigM;  
  
  // Defines third dimension, or fuzziness. M^2 + ldel < M^2 < M^2 + rdel

  double _ldel;
  double _rdel;

}; 

#endif


