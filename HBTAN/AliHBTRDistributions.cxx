
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

//__________________________________________________________________
////////////////////////////////////////////////////////////////////
//                                                                //
//
// Classes for Q's monitoring Vs Kt and Vs Qinv                   //
//                                                                //
// Author:                                                        //
// Zbigniew Chajecki <chajecki@if.pw.edu.pl>                      //
//                                                                //
////////////////////////////////////////////////////////////////////


#include "AliHBTRDistributions.h"

ClassImp(AliHBTRStarDistribution)


AliHBTRStarDistribution::AliHBTRStarDistribution(Int_t nXbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nXbins,maxXval,minXval)
{
//ctor
 Rename("RStarDistribution","R^{*} distribution");
}
/******************************************************************/

TH1* AliHBTRStarDistribution::GetResult()
{
 //returns the result histo
 return this->GetNumerator();
}


/******************************************************************/
/******************************************************************/

ClassImp(AliHBTRDistribution)

AliHBTRDistribution::AliHBTRDistribution(Int_t nXbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nXbins,maxXval,minXval)
{
//ctor
 Rename("RDistribution","R (distance between creation points) distribution ");
}

/******************************************************************/

TH1* AliHBTRDistribution::GetResult()
{
 //returns the result histo
 return this->GetNumerator();
}


/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/
