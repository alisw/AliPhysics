#ifndef ABSOCONST_H
#define ABSOCONST_H

#include "AliConst.h"

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */

/* $Id$ */

// start of 2deg cone
    const Float_t kZTwoDeg     = 128.863;           
// start of the absorber 
    const Float_t kZAbsStart   = 90.;    
// end of the W-nose
    const Float_t kZNose       = 102.;
// end of the 5deg line below the TPC field cage
    const Float_t kZConeTPC    = 285.;
// start of concrete absorber
    const Float_t kZAbsCc      = 315.; 

// max acceptance angle
    const Float_t kAccMax = 10.*kDegrad;
// angle of nose
    const Float_t kTheta1  = 24.*kDegrad;
// angle of second outer cone below field cage
    const Float_t kTheta2  = 5. *kDegrad;
// outer angler of W rear shield  
    const Float_t kThetaR = 3. *kDegrad;
//
// thicknesses defining the absorber
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// steel envelope
    const Float_t kDSteel = 1.;  
// poly-ethylene layer
    const Float_t kDPoly  = 5.0;




#endif

