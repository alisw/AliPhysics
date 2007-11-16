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
/// \class AliMUONRegionalTriggerBoard
/// Dimuon regional trigger implementation:
/// - entry are local board responses
/// - output is 12-bit word
/// - algorithm is similar to the global one
///
/// \author Rachid Guernane (LPCCFd)
/// Corrected by Christian Finck (Subatech)
//-----------------------------------------------------------------------------

#include "AliMUONRegionalTriggerBoard.h"

#include "AliLog.h"

#include "TBits.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONRegionalTriggerBoard)
/// \endcond

//___________________________________________
AliMUONRegionalTriggerBoard::AliMUONRegionalTriggerBoard()
  : AliMUONTriggerBoard(),
    fMask(0x0)   
{
/// Default constructor
   for (Int_t i=0; i<16; i++) fLocalResponse[i] = 0;
}

//___________________________________________
AliMUONRegionalTriggerBoard::AliMUONRegionalTriggerBoard(const char *name, Int_t a) 
  : AliMUONTriggerBoard(name, a),
    fMask(0x0)   
{
/// Standard constructor
   for (Int_t i=0; i<16; i++) fLocalResponse[i] = 0;
}

//___________________________________________
AliMUONRegionalTriggerBoard::~AliMUONRegionalTriggerBoard()
{
/// Destructor
}

//___________________________________________
void AliMUONRegionalTriggerBoard::Response()
{
/// response is given following the regional algorithm

  Int_t t[16];

   for (Int_t i = 0; i < 16; ++i)
   {
     if ((fMask >> i) & 0x1)
      t[i] = fLocalResponse[i];
     else
       t[i] = 0;
   }
   
   Int_t rank = 8;

   for (Int_t i = 0; i < 4; ++i)
   {
      Int_t ip = 0;
      
      for (Int_t j = 0; j < rank; ++j)
      {
         UShort_t lthres = Algo(t[2*j],t[2*j+1],"LPT",i);

         UShort_t hthres = Algo(t[2*j],t[2*j+1],"HPT",i); hthres <<= 4;

         t[ip] = lthres | hthres;

         ip++;
      }
      
      rank /= 2; 
   }

   fResponse = t[0]; // 8-bit [H4:L4]
}

//___________________________________________
UShort_t AliMUONRegionalTriggerBoard::Algo(UShort_t i, UShort_t j, char *thres, Int_t level)
{
/// implementation of the regional algorithm
/// similar to the global algorithm except for the
/// input layer

  TBits a(8), b(8); a.Set(8,&i); b.Set(8,&j);

   TBits trg1(2), trg2(2), trg(2);

   if (!strcmp(thres,"LPT"))
   {
      if (!level)
      {         
         trg1[0] = a[0]; trg1[1] = a[1]; 
         trg2[0] = b[0]; trg2[1] = b[1];
      }
      else
      {
         trg1[0] = a[2]; trg1[1] = a[3]; 
         trg2[0] = b[2]; trg2[1] = b[3];
      }
   }
   else
   {
      if (!level)
      {         
         trg1[0] = a[2]; trg1[1] = a[3]; 
         trg2[0] = b[2]; trg2[1] = b[3];
      }
      else
      {
         trg1[0] = a[6]; trg1[1] = a[7]; 
         trg2[0] = b[6]; trg2[1] = b[7];         
      }
   }
       
   TBits trgLS1(1), trgUS1(1), trgLS2(1), trgUS2(1), trgLS(1), trgUS(1);

   if (!level) 
   {
      trgLS1[0] = trgUS1[0] = trgLS2[0] = trgUS2[0] = 0;
   }
   else
   {
       if (!strcmp(thres,"LPT"))
      {
         trgLS1[0] = a[1]; trgUS1[0] = a[0]; 
         trgLS2[0] = b[1]; trgUS2[0] = b[0];
      }
      else
      {
         trgLS1[0] = a[5]; trgUS1[0] = a[4]; 
         trgLS2[0] = b[5]; trgUS2[0] = b[4];         
      }
   }

   trgLS[0] = ( trg1[0] & trg2[0] ) | ( trg1[1] & trg2[1] ) | trgLS1[0] | trgLS2[0];
   trgUS[0] = ( trg1[0] & trg2[1] ) | ( trg1[1] & trg2[0] ) | trgUS1[0] | trgUS2[0];
   
   trg[0] = trg1[0] | trg2[0];
   trg[1] = trg1[1] | trg2[1];
   
   TBits v(4);
   
   v[0] = trgUS[0];
   v[1] = trgLS[0];
   v[2] = trg[0];
   v[3] = trg[1];
   
   UShort_t rv = 0;
   v.Get(&rv);

   return rv;
}

//___________________________________________
void AliMUONRegionalTriggerBoard::Scan(Option_t*) const
{
/// scan local board entries 

  for (Int_t i=0; i<16; i++) 
   {
      TBits b;
      b.Set(6,&fLocalResponse[i]);
      
      cout << "Entry " << i << " is " << b << endl;
      
   }
   
}
//___________________________________________
void AliMUONRegionalTriggerBoard::Mask(UShort_t mask)
{
/// mask entry index

    fMask = mask;
}



