#ifndef ABSOCONST_H
#define ABSOCONST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */

/* $Id$ */

// start of 2deg cone
    const Float_t zTwoDeg     = 117.40864;           
// start of the absorber 
    const Float_t zAbsStart   = 90.;    
// end of the W-nose
    const Float_t zNose       = 102.;
// end of the 5deg line below the TPC field cage
    const Float_t zConeTPC    = 285.;
// start of concrete absorber
    const Float_t zAbsCc      = 315.; 

// max acceptance angle
    const Float_t accMax = 10.*kDegrad;
// angle of nose
    const Float_t theta1  = 24.*kDegrad;
// angle of second outer cone below field cage
    const Float_t theta2  = 5. *kDegrad;
// outer angler of W rear shield  
    const Float_t thetaR = 3. *kDegrad;
//
// thicknesses defining the absorber
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// steel envelope
    const Float_t dSteel = 1.;  
// poly-ethylene layer
    const Float_t dPoly  = 5.0;




#endif
