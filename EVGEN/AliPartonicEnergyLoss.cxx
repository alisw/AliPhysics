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

/*
$Log$
Revision 1.1.2.1  2002/07/24 08:56:29  alibrary
Updating EVGEN on TVirtulaMC

Revision 1.1  2002/07/15 13:12:33  morsch
Library class for parton quenching calculations.

*/


#include "AliPartonicEnergyLoss.h"
#include <TMath.h>


#ifndef WIN32
# define lookup   lookup_
# define read     read_
# define readca3  readca3_
# define type_of_call
#else
# define lookup   LOOKUP
# define read     READ
# define readca3  READCA3
# define type_of_call _stdcall
#endif


extern "C" void type_of_call lookup(Double_t &, Double_t &, 
				    Double_t &, Double_t &);
extern "C" void type_of_call read();
extern "C" void type_of_call readca3();


ClassImp(AliPartonicEnergyLoss)

void AliPartonicEnergyLoss::
QuenchingWeight(Double_t r, Double_t x, Double_t& cont, Double_t& disc)
{
//
//  Calculate quenching weight
//
    lookup(r,x,cont,disc);
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
    read();
    readca3();
}

