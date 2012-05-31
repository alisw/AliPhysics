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
// Module: EvtGen/EvtStdHep.hh
//
// Description: Class produce the StdHep representation of the decay.
//
// Modification history:
//
//    RYD     March. 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSTDHEP_HH
#define EVTSTDHEP_HH

#include "EvtGenBase/EvtVector4R.hh"
#include <iosfwd>

const int EVTSTDHEPLENGTH =1000;

class EvtStdHep {

public:

  EvtStdHep(){}
  ~EvtStdHep(){}

  void init();

  int getFirstMother(int i) { return _prntfirst[i]; }
  int getLastMother(int i) { return _prntlast[i]; }
  int getFirstDaughter(int i) { return _daugfirst[i]; }
  int getLastDaughter(int i) { return _dauglast[i]; }
  
  int getStdHepID(int i) { return _id[i]; }
  int getIStat(int i) { return _istat[i]; }

  EvtVector4R getP4(int i) { return _p4[i]; }
  EvtVector4R getX4(int i) { return _x[i]; }

  void translate(EvtVector4R d);

  int getNPart();
  void createParticle(EvtVector4R p4,EvtVector4R x,int prntfirst,
		     int prntlast, int id);

  friend std::ostream& operator<<(std::ostream& s, const EvtStdHep& stdhep);  

private:

  int _npart;
  EvtVector4R _p4[EVTSTDHEPLENGTH];
  EvtVector4R _x[EVTSTDHEPLENGTH];
  int _prntfirst[EVTSTDHEPLENGTH];
  int _prntlast[EVTSTDHEPLENGTH];
  int _daugfirst[EVTSTDHEPLENGTH];
  int _dauglast[EVTSTDHEPLENGTH];
  int _id[EVTSTDHEPLENGTH];
  int _istat[EVTSTDHEPLENGTH];
  
}; 

#endif

