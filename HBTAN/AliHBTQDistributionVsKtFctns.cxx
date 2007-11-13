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
// class AliHBTQInvDistributionVsKtFctn                           //
// class AliHBTQOutDistributionVsKtFctn                           //
// class AliHBTQSideDistributionVsKtFctn                          //
// class AliHBTQLongDistributionVsKtFctn                          //
//                                                                //
// Classes for Q's monitoring Vs Kt and Vs Qinv                   //
//                                                                //
// Author:                                                        //
// Zbigniew Chajecki <chajecki@if.pw.edu.pl>                      //
//                                                                //
////////////////////////////////////////////////////////////////////

/******************************************************************/
/******************************************************************/

#include "AliHBTQDistributionVsKtFctns.h"

ClassImp( AliHBTQInvDistributionVsKtFctn )

AliHBTQInvDistributionVsKtFctn::AliHBTQInvDistributionVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                                                               Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 //ctor
 Rename("QInvDistributionVsKt","Q_{Inv} Distribution vs. K_{t}");
}
/******************************************************************/

TH1* AliHBTQInvDistributionVsKtFctn::GetResult()
{
 //returns the result histo
 return this->GetNumerator();
}

/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQOutDistributionVsKtFctn )

AliHBTQOutDistributionVsKtFctn::AliHBTQOutDistributionVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                                                               Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 //ctor
 Rename("QOutDistributionVsKt","Q_{Out} Distribution vs. K_{t}");
}
/******************************************************************/

TH1* AliHBTQOutDistributionVsKtFctn::GetResult()
{
 //returns the result histo
 return this->GetNumerator();
}

/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQSideDistributionVsKtFctn )

AliHBTQSideDistributionVsKtFctn::AliHBTQSideDistributionVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                                                                 Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 //ctor
 Rename("QSideDistributionVsKt","Q_{Side} Distribution vs. K_{t}");
}
/******************************************************************/

TH1* AliHBTQSideDistributionVsKtFctn::GetResult()
{
 //returns the result histo
 return this->GetNumerator();
}

/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQLongDistributionVsKtFctn )

AliHBTQLongDistributionVsKtFctn::AliHBTQLongDistributionVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                                                                 Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 //ctor
 Rename("QLongDistributionVsKt","Q_{Long} Distribution vs. K_{t}");
}
/******************************************************************/

TH1* AliHBTQLongDistributionVsKtFctn::GetResult()
{
 //returns the result histo
 return this->GetNumerator();
}

