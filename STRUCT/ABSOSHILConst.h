#ifndef ABSOSHILCONST_H
#define ABSOSHILCONST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */

/* $Id$ */

//
// z-positions defining the absorber
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// start of inner opening cone
    const Float_t kZOpen   = 300.;
// rear end of the absorber
    const Float_t kZRear   = 503.;
// thickness of rear shield
    const Float_t kDRear  =  35.;
//
// angles defining the absorber
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// min acceptance angle
    const Float_t kAccMin = 2. *kDegrad;
// angle of first opening cone
    Float_t kThetaOpen1 = 0.75*kDegrad;    
//
// inner radius of heavy shield
    Float_t kRAbs = 4.5;
// innner radius of beam tube
    Float_t kRVacu=2.9;                     
#endif










