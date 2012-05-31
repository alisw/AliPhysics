// $Id: rly.cc,v 1.1 2009/02/18 03:34:11 ryd Exp $
//
// Define random number generators used by JETSET

#include "EvtGenBase/EvtRandom.hh"

extern "C" {
  extern float rly_();
}

float rly_(){
  return EvtRandom::Flat();
}
