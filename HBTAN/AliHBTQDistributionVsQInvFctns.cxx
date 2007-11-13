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
// class AliHBTQOutDistributionVsQInvFctn                         //
// class AliHBTQSideDistributionVsQInvFctn                        //
// class AliHBTQLongDistributionVsQInvFctn                        //
// class AliHBTPtDiffDistributionVsQInvFctn                       //
//                                                                //
// Classes for Q's monitoring Vs Kt and Vs Qinv                   //
//                                                                //
// Author:                                                        //
// Zbigniew Chajecki <chajecki@if.pw.edu.pl>                      //
//                                                                //
////////////////////////////////////////////////////////////////////



#include "AliHBTQDistributionVsQInvFctns.h"

ClassImp( AliHBTQOutDistributionVsQInvFctn )

AliHBTQOutDistributionVsQInvFctn::AliHBTQOutDistributionVsQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                                                                   Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 //ctor
 Rename("QOutDistributionVsQInv","Q_{Out} Distribution vs. Q_{inv}");
}
/******************************************************************/

TH1* AliHBTQOutDistributionVsQInvFctn::GetResult()
{
 //returns the result histo
 return this->GetNumerator();
}


/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQSideDistributionVsQInvFctn )

AliHBTQSideDistributionVsQInvFctn::AliHBTQSideDistributionVsQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                                                                     Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 //ctor
 Rename("QSideDistributionVsQInv","Q_{Side} Distribution vs. Q_{inv}");
}
/******************************************************************/

TH1* AliHBTQSideDistributionVsQInvFctn::GetResult()
{
 //returns the result histo
 return this->GetNumerator();
}


/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQLongDistributionVsQInvFctn )

AliHBTQLongDistributionVsQInvFctn::AliHBTQLongDistributionVsQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                                                                     Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 //ctor
 Rename("QLongDistributionVsQInv","Q_{Long} Distribution vs. Q_{inv}");
}
/******************************************************************/

TH1* AliHBTQLongDistributionVsQInvFctn::GetResult()
{
 //returns the result histo
 return this->GetNumerator();
}


/******************************************************************/
/******************************************************************/

ClassImp( AliHBTPtDiffDistributionVsQInvFctn )

AliHBTPtDiffDistributionVsQInvFctn::AliHBTPtDiffDistributionVsQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                                                                     Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 //ctor
 Rename("PtDiffDistributionVsQInv","P_{t} Difference Distribution vs. Q_{inv}");
}
/******************************************************************/

TH1* AliHBTPtDiffDistributionVsQInvFctn::GetResult()
{
 //returns the result histo
 return this->GetNumerator();
}

