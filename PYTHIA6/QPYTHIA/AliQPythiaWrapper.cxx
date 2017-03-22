/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------------------------------
//   Class: AliQPythiaWrapper
//   Responsibilities: Interface to C++ Functionality needed by the QPythia Fortran Code 
//                     So far only AliFastGlauber is interfaced

#include <TMath.h>
#include <AliFastGlauber.h>

#include "AliQPythiaWrapper.h"

ClassImp(AliQPythiaWrapper)


//_______________________________________________________________________
#define getlength getlength_

extern "C" {
    void getlength(Double_t& ell, Double_t& b) 
    {
	// Wrapper to method GetLength from AliGlauber
	(AliFastGlauber::Instance())->GetLength(ell, b);
    }
}


//_______________________________________________________________________
#define getrandombhard getrandombhard_

extern "C" {
    void getrandombhard(Double_t& b) 
    {
	// Wrapper to method GetRandomBHard from AliGlauber
	(AliFastGlauber::Instance())->GetRandomBHard(b);
    }
}

//_______________________________________________________________________
#define getrandomxy getrandomxy_

extern "C" {
    void getrandomxy(Double_t& x,Double_t& y ) 
    {
	// Wrapper to method GetRandomXY from AliGlauber
	(AliFastGlauber::Instance())->GetRandomXY(x,y);
    }
}


//_______________________________________________________________________
#define calculatei0i1 calculatei0i1_

extern "C" {
  void calculatei0i1(Double_t& xintegral0,Double_t& xintegral1,Double_t& b,Double_t& x,Double_t& y,Double_t& phi,Double_t& ellCut)
    {
	// Wrapper to method GetRandomXY from AliGlauber
      (AliFastGlauber::Instance())->CalculateI0I1(xintegral0,xintegral1,b,x,y,phi,ellCut);
    }
}


//_______________________________________________________________________
#define savexy savexy_

extern "C" {
 void savexy(Double_t& x, Double_t& y)
    {
	// Wrapper to save current production point
	return (AliFastGlauber::Instance())->SaveXY(x,y);
    }
}

#define savei0i1 savei0i1_

extern "C" {
 void savei0i1(Double_t& i0, Double_t& i1)
    {
	// Wrapper to save current integrals I0 and I1
	return (AliFastGlauber::Instance())->SaveI0I1(i0, i1);
    }
}
