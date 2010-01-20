//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtRandomEngine.cc
//
// Description: routines to generate random numbers
//              really trivial random number 
//              implementation.
//
// Modification history:
//
//    RYD      December 25, 1999           Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"


#include <stdio.h>
#include <math.h>
#include <iostream>
#include "EvtGenBase/EvtSimpleRandomEngine.hh"

double EvtSimpleRandomEngine::random(){
  
  _next=_next*1103515245+123345;
  unsigned temp=(unsigned)(_next/65536) % 32768;
  
  return ( temp + 1.0 ) / 32769.0;

}



