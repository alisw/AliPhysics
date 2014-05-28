/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtCyclic3.hh,v 1.2 2009-03-16 16:42:46 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Cyclic permutations of three indices A,B,C and their parings

#ifndef EVT_CYCLIC3_HH
#define EVT_CYCLIC3_HH

#include <iosfwd>

namespace EvtCyclic3 {


  enum Index {A=0,B=1,C=2};
  enum Pair  {BC=0,CB=BC,CA=1,AC=CA,AB=2,BA=AB}; 
  enum Perm {ABC=0,BCA=1,CAB=2,CBA=3,BAC=4,ACB=5};

  // Permutations (multiplication is not transitive)

  Index permute(Index i, Perm p);
  Perm  permutation(Index i1,Index i2,Index i3);
  Perm  permute(Perm i, Perm p);
  Pair permute(Pair i, Perm p); 

  Pair i2pair(int i);
  
  // Index-to-index

  Index prev(Index i); 
  Index next(Index i) ;
  Index other(Index i, Index j);
 
  // Index-to-pair

  Pair other(Index i); 
  Pair combine(Index i, Index j);


  // Pair-to-pair conversions

  Pair prev(Pair i);
  Pair next(Pair i);
  Pair other(Pair i, Pair j);
  
  // Pair-to-index conversions


  Index first(Pair i);
  Index second(Pair i); 
  Index other(Pair i) ;
  Index common(Pair i, Pair j);

  // String to Index, Pair

  Index strToIndex(const char* str);
  Pair   strToPair(const char* str);

  // To string conversions

  const char* c_str(Index i);
  const char* c_str(Pair i);
  const char* c_str(Perm i);

  // Useful name strings

  char* append(const char* str, EvtCyclic3::Index i);
  char* append(const char* str, EvtCyclic3::Pair i);

};

//where should these go?
//ostream& operator<<(ostream&, EvtCyclic3::Index);
//ostream& operator<<(ostream&, EvtCyclic3::Pair);

#endif
