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
Revision 1.6  2002/10/14 14:57:32  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.4.12.1  2002/07/24 09:50:10  alibrary
Updating VirtualMC

Revision 1.5  2002/07/23 11:48:05  alla
new Digits structure

Revision 1.4  2000/10/13 13:14:08  hristov
Bug fixes and code cleaning

Revision 1.3  2000/07/13 16:41:29  fca
New START corrected for coding conventions

Revision 1.2  2000/03/24 17:40:35  alla
New AliSTART

*/ 

#include <TArrayI.h>
#include "AliSTARTdigit.h"
#include <Riostream.h>

ClassImp(AliSTARTdigit)

//------------------------------------
 AliSTARTdigit::AliSTARTdigit() : TObject()
{
  fTimeAverage=9999;
  fTimeDiff=9999;
  fTimeBestRight=9999;
  fTimeBestLeft=9999;

  fTimeRight = new TArrayI(12);  
  fTimeLeft  = new TArrayI(12);  
  fADCRight = new TArrayI(12);  
  fADCLeft  = new TArrayI(12);  
}
//-----------------------------------
void AliSTARTdigit::SetTimeRight (TArrayI &o)
{
  fTimeRight = new TArrayI(12);  

  Int_t i;
  for (i=0; i<12; i++)
    {
      Int_t buf=o.At(i);
      fTimeRight->AddAt(buf,i);
    }
}
//--------------------------------------------
void AliSTARTdigit::SetTimeLeft (TArrayI &o)
{

  fTimeLeft  = new TArrayI(12);  
  Int_t i;
  for (i=0; i<12; i++)
    {
      Int_t buf=o.At(i);
      fTimeLeft->AddAt(buf,i);
    }
}
//--------------------------------------------
void AliSTARTdigit::GetTimeLeft (TArrayI &o)
{

  Int_t i;
  for (i=0; i<12; i++)
    {
      o[i]=fTimeLeft->At(i);
    }
}
//--------------------------------------------
void AliSTARTdigit::GetTimeRight (TArrayI &o)
{

  Int_t i;
  for (i=0; i<12; i++)
    {
      o[i]=fTimeRight->At(i);
    }
}
//--------------------------------------------
void AliSTARTdigit::GetADCLeft (TArrayI &o)
{

  Int_t i;
  for (i=0; i<12; i++)
    {
      o[i]=fADCLeft->At(i);
    }
}
//--------------------------------------------
void AliSTARTdigit::GetADCRight (TArrayI &o)
{

 Int_t i;
  for (i=0; i<12; i++)
    {
      o[i]=fADCRight->At(i);
    }
}
//--------------------------------------------
void AliSTARTdigit::SetADCLeft (TArrayI &o)
{

  fADCLeft  = new TArrayI(12);  
  Int_t i;
  //  Float_t fProcessKoef=1; // for pb 0.001
  for (i=0; i<12; i++)
    {
      Int_t buf=(o.At(i));
      fADCLeft->AddAt(buf,i);
    }
}
//--------------------------------------------
void AliSTARTdigit::SetADCRight (TArrayI &o)
{

  //  Float_t fProcessKoef=1; // for pb 0.001
  fADCRight  = new TArrayI(12);  
  Int_t i;
  for (i=0; i<12; i++)
    {
      Int_t buf=(o.At(i));
      fADCRight->AddAt(buf,i);
    }
}
//------------------------------------------------------
void AliSTARTdigit::Print()
{
  printf("AliSTARTdigit: fTimeAverage=%d, fTimeDiff=%d\n",
	 fTimeAverage, fTimeDiff);
  cout<<" BestTimeRigh "<<fTimeBestRight<<
    " TimeBestLeft "<<fTimeBestLeft<<endl;



}
