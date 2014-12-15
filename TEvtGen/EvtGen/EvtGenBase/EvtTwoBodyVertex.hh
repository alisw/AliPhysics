/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtTwoBodyVertex.hh,v 1.2 2009-03-16 16:34:38 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Two-body propagator vertex AB->A,B with an attached Blatt-Weisskopf form factor.

#ifndef EVT_TWO_BODY_VERTEX_HH
#define EVT_TWO_BODY_VERTEX_HH

#include <iostream>
#include "EvtGenBase/EvtTwoBodyKine.hh"
#include "EvtGenBase/EvtBlattWeisskopf.hh"

#include <iosfwd>

class EvtTwoBodyVertex {

public:

  EvtTwoBodyVertex();
  EvtTwoBodyVertex(double mA, double mB, double mAB, int L);
  EvtTwoBodyVertex(const EvtTwoBodyVertex& other);
  ~EvtTwoBodyVertex();

  double widthFactor(EvtTwoBodyKine x) const;
  double formFactor(EvtTwoBodyKine x) const;
  double phaseSpaceFactor(EvtTwoBodyKine x, EvtTwoBodyKine::Index) const;

  inline int L() const { return _LL; }
  inline double mA() const { return _kine.mA(); }
  inline double mB() const { return _kine.mB(); }
  inline double mAB() const { return _kine.mAB(); }
  inline double pD() const { return _p0; }
  void print(std::ostream& os) const; 

  void set_f(double R);

private:

  EvtTwoBodyKine _kine;
  int _LL;  
  double _p0; 
  EvtBlattWeisskopf* _f;  // optional Blatt-Weisskopf form factor

};

std::ostream& operator<<(std::ostream& os, const EvtTwoBodyVertex& v);

#endif
