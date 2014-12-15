/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtDalitzCoord.hh,v 1.2 2009-03-16 16:43:40 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Two dimensional coordinate of a point in a Dalitz plot

#ifndef EVT_DALITZ_COORD_HH
#define EVT_DALITZ_COORD_HH

#include "EvtGenBase/EvtCyclic3.hh"

#include <iostream>

class EvtDalitzCoord {
  
  
public:
  
  // ctor, dtor

  EvtDalitzCoord(); 
  EvtDalitzCoord(EvtCyclic3::Pair i1, double q1, EvtCyclic3::Pair i2, double q2);
  EvtDalitzCoord(const EvtDalitzCoord& other); 
  ~EvtDalitzCoord();

  inline EvtCyclic3::Pair pair1() const { return _i1; }
  inline EvtCyclic3::Pair pair2() const { return _i2; }
  inline double q1() const { return _q1; }
  inline double q2() const { return _q2; }


  // It's nice to have an equality operator for
  // a coordinate. However, beware effects of numerical precision
  
  bool operator==(const EvtDalitzCoord&) const;

  void print(std::ostream&) const;

private:

  // Two coordinates define the point

  EvtCyclic3::Pair _i1;
  EvtCyclic3::Pair _i2;

  double _q1;
  double _q2;
}; 

std::ostream& operator<<(std::ostream&,const EvtDalitzCoord&);

#endif






