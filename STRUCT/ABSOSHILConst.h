#ifndef ABSOSHILCONST_H
#define ABSOSHILCONST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */

/* $Id$ */

//
// z-positions defining the absorber
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// start of inner opening cone
    const Float_t zOpen   = 300.;
// rear end of the absorber
    const Float_t zRear   = 503.;
// thickness of rear shield
    const Float_t dRear  =  35.;
//
// angles defining the absorber
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// min acceptance angle
    const Float_t accMin = 2. *kDegrad;
// angle of first opening cone
    Float_t thetaOpen1 = 0.78*kDegrad;    
//
// inner radius of heavy shield
    Float_t rAbs = 4.1;
// innner radius of beam tube
    Float_t rVacu=2.9;                     
#endif

