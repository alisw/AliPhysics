/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtAbsBinning.hh,v 1.2 2009-03-16 16:43:39 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *****************************************************************************/

/*
 * Data point to bin value mapping
 */ 

#ifndef EVT_ABS_BINNING_HH
#define EVT_ABS_BINNING_HH
#define BIN_OUTSIDE -1 

#include <stdio.h>

template <class T> class EvtAbsBinning {
public:

  EvtAbsBinning()
  {}
  EvtAbsBinning(const EvtAbsBinning<T>& other)
  {}
  virtual ~EvtAbsBinning() 
  {}
  
  virtual EvtAbsBinning<T>* clone() const = 0;
  virtual int getBin(const T& point) const = 0;
  virtual T getBinPoint(int bin) const = 0;
  virtual double size(int bin) const = 0;

  virtual int nTypes() const = 0;

  virtual char* typeLabel(int i) const
  {
    char* a = new char[128];
    sprintf(a,"%d",i);
    return a;
  }

};

#endif
