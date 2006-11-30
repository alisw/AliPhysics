// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright & copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTTransBit.h"

#if __GNUC__ == 3
using namespace std;
#endif

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


/** \class AliHLTTransbit
<pre>
//_____________________________________________________________
// AliHLTTransbit (taken from the offline AliROOT code,
// original author: Marian Ivanov, GSI Darmstadt for AliROOT)
//
// Time Projection Chamber ADC bit compresion lookup table
//
//  Conversion equation:
//    For AliTransBitV1
//    dy/dx= Int_t(1+x/fX0)
//
//    For AliTransBitV2
//    y  =(2**bit0) ln(1+x/fX0) / ln(1+(2**bit1-1)/fX0)
//
//    where x0 is calculated in function GetOptimumX0()
//
//  Example session in aliroot:
//    Int_t b0=10;  // original number of bits
//    Int_t b1=8;   // compressed
//
//    AliTransBitV1 trans;
//    Int_t x0=TMath::Nint(TMath::Exp(b0*TMath::Log(2)));
//    Int_t x1=TMath::Nint(TMath::Exp(b1*TMath::Log(2)));
//    trans.SetBits(b0,b1);
//    trans.FindOptimumX0();
//    trans.Update();
//    cout<<trans.Get0to1(x0-2)<<"\n";
//    cout<<trans.Get1to0(x1-2)<<"\n";
//
//    // to produce table
//    for( Int_t i=0;i<x1;i++) cout<<i<<"\t"<<trans.Get1to0(i)<<"\n"; > table1to0.txt
//    for( Int_t i=0;i<x0;i++) cout<<i<<"\t"<<trans.Get0to1(i)<<"\n"; > table0to1.txt
//
//    for( Int_t i=0;i<x1-1;i++) cout<<i<<"\t"<<trans.Get1to0(i+1)-trans.Get1to0(i)<<"\n"; > tabled.txt
//
</pre>
*/

ClassImp(AliHLTTransBit)
ClassImp(AliHLTTransBitV1)
ClassImp(AliHLTTransBitV2)

AliHLTTransBit::AliHLTTransBit()
{
  // default constructor
  fTable0 = 0;
  fTable1 = 0;
  fBit0   = 10;
  fBit1   = 8;
  fX0     = 0;
}

AliHLTTransBit::~AliHLTTransBit() 
{
  // destructor
  if (fTable0!=0) delete [] fTable0;
  if (fTable1!=0) delete [] fTable1;
}

Double_t AliHLTTransBitV1::FindOptimumX0()
{
  //find x0 for which derivation at xder1 is equal 1
  Int_t x0=(Int_t)rint(exp(fBit0*log(2.)));
  Int_t x1=(Int_t)rint(exp(fBit1*log(2.)));

  fX0 = ((x1-2)*(x1-2)/2.)/(x0-x1-1);  //starting fX0
  Int_t digit=0;
  Int_t j;
  Double_t weight=1;
  for(j=0;(j<50)&&(digit!=(x0-1));j++){  
    Int_t olddigit=digit;
    digit=0;
    for (Int_t  i=0; i<x1-1;i++) digit+= Int_t(1+Double_t(i)/fX0);      
    fX0*= (1.+weight*Double_t(digit)/Double_t(x0-1))/(1.+weight);    
    if ( ((olddigit-(x0-1)) * (digit-(x0-1))) <0 ) 
      weight*=0.25;
    // cout<<j<<"\t"<<fX0<<"\t"<<digit<<"\n";
  }
  //  cout<<"End of iteration "<<j<<"\n";
  //cout<<"digit..."<<digit<<"\n";  
  return fX0;
}

void AliHLTTransBitV1::Update()
{
  //construct lookup tables for loosy compression from 
  if (fX0<1) fX0 = FindOptimumX0();  
  Int_t x0=(Int_t)rint(exp(fBit0*log(2.)));
  Int_t x1=(Int_t)rint(exp(fBit1*log(2.)));
  
  //fTable0 - conversion from bit0 coding to bit1 coding
  if (fTable0!=0) delete fTable0;
  fTable0= new Int_t[x0];
  //fTable1 - conversion table from bit1 to bit0 coding
  if (fTable1!=0) delete [] fTable1;
  fTable1= new Int_t[x1];

  Int_t digit=0;
  for (Int_t  i=0; i<x1-1;i++){    
    Int_t ddig=Int_t(1+Double_t(i)/fX0);
    for (Int_t j=0;j<ddig;j++) if ((digit+j)<x0) fTable0[digit+j] = i;    
    fTable1[i]= digit+ddig/2;
    digit+= ddig;
  }

  for (Int_t i=digit;i<x0-1;i++) fTable0[i]=x1-2;
  fTable1[x1-1]=x0-1; //overflow 
  fTable0[x0-1]=x1-1; //overflow    
  return;
}

Double_t AliHLTTransBitV2::FindOptimumX0()
{
  //find x0 for which derivation at xder1 is equal 1
  const Float_t kxder1=1;
  const Float_t kdx=0.1;

  Float_t x0=exp(fBit0*log(2.));
  Float_t x1=exp(fBit1*log(2.));
  Float_t deriv = 0;
  Float_t x;
  for (x=x1;( (x>1)&&(deriv<1)) ;x-=kdx)
    {
      deriv = (x1-1)/( log(1.+x0/x) *x *(1+kxder1/x));
    }
  x+=kdx/2.;
  fX0 = x;
  return fX0;
}

void AliHLTTransBitV2::Update()
{
  //construct lookup tables for loosy compresion from 
  if (fX0<1) fX0 = FindOptimumX0();  
  //Float_t x0=(Int_t)rint(exp(fBit0*log(2.)));
  //Float_t x1=(Int_t)rint(exp(fBit1*log(2.)));
  Int_t x0=(Int_t)rint(exp(fBit0*log(2.)));
  Int_t x1=(Int_t)rint(exp(fBit1*log(2.)));
  
  //fTable0 - conversion from bit0 coding to bit1 coding
  if (fTable0!=0) delete fTable0;
  fTable0= new Int_t[x0];

  //fTable1 - conversion table from bit1 to bit0 coding
  if (fTable1!=0) delete [] fTable1;
  fTable1= new Int_t[x1];

  Int_t i;

  for ( i=0; i<x0;i++)
      fTable0[i] =(Int_t)rint((x1-0.501)*log(1.+Float_t(i)/fX0)/
				 log(1.+(x0-1)/fX0));

  Int_t old0=-1;
  Int_t old1=-1;
  Int_t new1;
  for ( i=0; i<x0;i++){
      new1 = fTable0[i];
      if (new1!=old1){
	if (old1>=0) fTable1[old1] =(old0+i)/2;
	old0 = i;
	old1 = new1;
      }
  }
  fTable1[old1]=(Int_t)rint((Float_t)(old0+x0)/2);
  
  return;
}
