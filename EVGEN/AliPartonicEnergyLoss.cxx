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
// Calculate parton energy loss
// in nucleus nucleus
// reactions via calculation
// of quenching weights
// Author: A.Morsch
//

#include "AliPartonicEnergyLoss.h"
#include <TMath.h>


#ifndef WIN32
# define eloss_lookup   eloss_lookup_
# define eloss_read     eloss_read_
# define eloss_readca3  eloss_readca3_
# define type_of_call
#else
# define eloss_lookup   ELOSS_LOOKUP
# define eloss_read     ELOSS_READ
# define eloss_readca3  ELOSS_READCA3
# define type_of_call _stdcall
#endif


extern "C" void type_of_call eloss_lookup(Double_t &, Double_t &, 
				    Double_t &, Double_t &);
extern "C" void type_of_call eloss_read();
extern "C" void type_of_call eloss_readca3();


ClassImp(AliPartonicEnergyLoss)

void AliPartonicEnergyLoss::
QuenchingWeight(Double_t r, Double_t x, Double_t& cont, Double_t& disc)
{
//
//  Calculate quenching weight
//
    eloss_lookup(r,x,cont,disc);
}

void AliPartonicEnergyLoss::RunTest()
{
//
// Test routine
//
    Double_t cont, disc, rri, wwt;
    
    Init();

    rri = 400.;
    for (Int_t i = 1; i <= 40; i++) 
    {
	wwt = TMath::Power(10., (-3.+5.*(i-1)/49.));
	QuenchingWeight(rri,wwt,cont,disc);
	printf("%d %e %e\n", i, wwt, cont);
    }
}
 
void AliPartonicEnergyLoss::Init()
{
//
//  Read data
//
    eloss_read();
    eloss_readca3();
}

