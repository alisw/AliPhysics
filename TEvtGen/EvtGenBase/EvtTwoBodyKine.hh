/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtTwoBodyKine.hh,v 1.2 2009-03-16 16:34:38 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Descriptions of the kinematics of a two-body decay.

#ifndef EVT_TWO_BODY_KINE_HH
#define EVT_TWO_BODY_KINE_HH

#include <iostream>

class EvtTwoBodyKine {

public:
  
  enum Index {A,B,AB}; 

  EvtTwoBodyKine();
  EvtTwoBodyKine(double mA, double mB, double mAB);
  EvtTwoBodyKine(const EvtTwoBodyKine& other);
  ~EvtTwoBodyKine();

  // Accessors

  inline double mA()  const { return _mA; }
  inline double mB()  const { return _mB; }
  inline double mAB() const { return _mAB; } 
  double m(Index i) const;

  // Momentum of the other two particles in the 
  // rest-frame of particle i.

  double p(Index i = AB) const;

  // Energy of particle i in the rest frame of particle j

  double e(Index i, Index j) const;

  void print(std::ostream& os) const;
  
private:
  
  double _mA;
  double _mB;
  double _mAB;
};

std::ostream& operator<<(std::ostream& os, const EvtTwoBodyKine& p);

#endif
