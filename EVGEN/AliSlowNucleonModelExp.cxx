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

//
// Experimental data inspired Gray Particle Model for p-Pb collisions
// The number of gray nucleons  is proportional to the number of collisions.
// The number of black nucleons is proportional to the number of collisions
// Fluctuations are calculated from a binomial distribution.
// Author: A.Morsch
//

#include "AliSlowNucleonModelExp.h"
#include "AliCollisionGeometry.h"
#include <TRandom.h>

ClassImp(AliSlowNucleonModelExp);


AliSlowNucleonModelExp::AliSlowNucleonModelExp():
    fP(82),
    fN (208 - 82),
    fAlphaGray(2.),
    fAlphaBlack(4.)
{
  //
  // Default constructor
  //
}


void AliSlowNucleonModelExp::GetNumberOfSlowNucleons(AliCollisionGeometry* geo, 
						      Int_t& ngp, Int_t& ngn, Int_t & nbp, Int_t & nbn) const
{
//
// Return the number of black and gray nucleons
//
// Number of collisions

    Int_t nu  = geo->NwN() +  geo->NNw();

// Mean number of gray nucleons 

    Float_t nGray         = fAlphaGray * nu;
    Float_t nGrayNeutrons = nGray * fN / (fN + fP);
    Float_t nGrayProtons  = nGray - nGrayNeutrons;

// Mean number of black nucleons 
    Float_t nBlack         = fAlphaBlack * nu;
    Float_t nBlackNeutrons = nBlack * fN / (fN + fP);
    Float_t nBlackProtons  = nBlack - nBlackNeutrons;

// Actual number (including fluctuations) from binomial distribution
    Double_t p;

//  gray neutrons
    p =  nGrayNeutrons/fN;
    ngn = gRandom->Binomial((Int_t) fN, p);
    
//  gray protons
    p =  nGrayProtons/fP;
    ngp = gRandom->Binomial((Int_t) fP, p);

//  black neutrons
    p =  nBlackNeutrons/fN;
    nbn = gRandom->Binomial((Int_t) fN, p);
    
//  black protons
    p =  nBlackProtons/fP;
    nbp = gRandom->Binomial((Int_t) fP, p);
	
}

void AliSlowNucleonModelExp::SetParameters(Float_t alpha1, Float_t alpha2)
{
    // Set the model parameters
    fAlphaGray  = alpha1;
    fAlphaBlack = alpha2;
}
