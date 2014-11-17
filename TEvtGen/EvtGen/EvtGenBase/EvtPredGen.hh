/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtPredGen.hh,v 1.2 2009-03-16 16:40:16 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// A predicate is applied to a generator to get another generator.
// Accept-reject can be implemented in this way.
//
//           Predicate
// Generator    ->     Generator 

#ifndef EVT_PRED_GEN_HH
#define EVT_PRED_GEN_HH

#include <stdio.h>

template <class Generator, class Predicate> 
class EvtPredGen {

public:

  typedef typename Generator::result_type result_type;
  
  EvtPredGen()
    : itsTried(0), itsPassed(0)
  {}

  EvtPredGen(Generator gen, Predicate pred)
    : itsGen(gen), itsPred(pred), itsTried(0), itsPassed(0) 
  {}

  EvtPredGen(const EvtPredGen& other)
    : itsGen(other.itsGen), itsPred(other.itsPred), 
    itsTried(other.itsTried), itsPassed(other.itsPassed)
  {}

  ~EvtPredGen()
  {}
  
  result_type operator()() {

    int i = 0;
    int MAX = 10000;
    while(i++ < MAX) {

      itsTried++;
      result_type point = itsGen();
      if(itsPred(point)) {
	itsPassed++;
	return point;
      }
    }    
    
    printf("No random point generated after %d attempts\n",MAX);
    printf("Sharp peak? Consider using pole compensation.\n");
    printf("I will now pick a point at random to return.\n");
    return itsGen();
  }
  
  inline int getTried() const { return itsTried; }
  inline int getPassed() const { return itsPassed; }

protected:

  Generator itsGen;
  Predicate itsPred;
  int itsTried;
  int itsPassed;

};

#endif

