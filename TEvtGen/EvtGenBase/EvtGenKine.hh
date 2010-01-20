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
// Module: EvtGen/EvtGenKine.hh
//
// Description:Tools for generating phase space.
//
// Modification history:
//
//    RYD     March 24, 1998         Module created
//
//------------------------------------------------------------------------
#ifndef EVTGENKINE_HH
#define EVTGENKINE_HH

class EvtVector4R;

class EvtGenKine{

public:

static double PhaseSpace( int ndaug, double mass[30],
			  EvtVector4R p4[30], double mp );

static double PhaseSpacePole(double M, double m1, double m2, double m3, 
			     double a,EvtVector4R p4[10]);


};

#endif

