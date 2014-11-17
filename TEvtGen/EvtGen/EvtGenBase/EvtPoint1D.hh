/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtPoint1D.hh,v 1.2 2009-03-16 16:40:16 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Point on a finite 1-D interval. isValid shows whether for a given specification,
// the coordinate _value is inside the interval defined by _min, _max.

#ifndef EVT_POINT_1D_HH
#define EVT_POINT_1D_HH

class EvtPoint1D {
public:

  EvtPoint1D();
  EvtPoint1D(double value);
  EvtPoint1D(double min, double max, double value);
  ~EvtPoint1D();

  bool isValid() const
  { 
    return _valid; 
  }

  double value() const 
  {
    return _value;
  }

  void print() const;

private:
  
  double _min;   // interval minimum
  double _max;   // interval maximum
  double _value;
  bool _valid;   // valid point inside the interval?

};

#endif

