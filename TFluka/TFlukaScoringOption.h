#ifndef TFLUKASCORINGOPTION
#define TFLUKASCORINGOPTION

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
// Class to store FLUKA specific scoring options
//                                                                           //
//
// Author: andreas.morsch@cern.ch 
//
//
///////////////////////////////////////////////////////////////////////////////


#include <TNamed.h>


class TFlukaScoringOption : public TNamed
{
public:
   // constructors
    TFlukaScoringOption();
    TFlukaScoringOption(const char* name, const char* sdum, Int_t npar,  Float_t what[12]);
 protected:
    Int_t         fNpar;        // Number of paramters
    Float_t       fWhat[12];    // WHAT()
    ClassDef(TFlukaScoringOption, 1)  // Fluka Scoring Option
};
	
#endif
	
