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
/// \class AliMUONGlobalTriggerBoard
/// Global trigger implementation:
/// - inputs are regional responses
/// - output is a 12-bit word
/// - 4 bits per trigger level
///
/// \author Rachid Guernane (LPCCFd), 
/// Corrected by Christian Finck (Subatech)
//-----------------------------------------------------------------------------

#include "AliMUONGlobalTriggerBoard.h"
#include "AliLog.h"
#include "TBits.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONGlobalTriggerBoard)
/// \endcond

//___________________________________________
AliMUONGlobalTriggerBoard::AliMUONGlobalTriggerBoard(): AliMUONTriggerBoard()
{
/// Default constructor

   for (Int_t i=0;i<16;i++) fRegionalResponse[i] = 0;
   for (Int_t i=0;i< 4;i++) fGlobalInput[i] = 0;
}

//___________________________________________
AliMUONGlobalTriggerBoard::AliMUONGlobalTriggerBoard(const char *name, Int_t a) : AliMUONTriggerBoard(name, a)
{
/// Standard constructor

   for (Int_t i=0;i<16;i++) fRegionalResponse[i] = 0;
   for (Int_t i=0;i< 4;i++) fGlobalInput[i] = 0;
}

//___________________________________________
AliMUONGlobalTriggerBoard::~AliMUONGlobalTriggerBoard()
{
/// Destructor
}

//___________________________________________
void AliMUONGlobalTriggerBoard::Mask(Int_t index, UInt_t mask)
{
  /// mask global trigger board input index with value mask
  if ( index >= 0 && index < 4 ) 
  {
    fMask[index]=mask;
  }
  else
  {
    AliError(Form("Index %d out of bounds (max %d)",index,3));
  }  
}

//___________________________________________
void AliMUONGlobalTriggerBoard::Response()
{
   /// compute the global trigger board
   /// response according to the algo() method
// output from global trigger algorithm
// [+, -, US, LS] * [Hpt, Lpt]
// transformed to [usHpt, usLpt, lsHpt, lsLpt, sHpt, sLpt] according
// to Global Trigger Unit user manual

   Int_t t[16];

   BuildGlobalInput();
   MaskGlobalInput();

   for (Int_t i = 0; i < 16; ++i) 
   {
     t[i] = fRegionalResponse[i];
   }
   
   
   Int_t rank = 8;

   for (Int_t i=0;i<4;i++)
   {
      Int_t ip = 0;
      
      for (Int_t j=0;j<rank;j++)
      {
         UShort_t lthres = Algo(t[2*j],t[2*j+1],"LPT");

         UShort_t hthres = Algo(t[2*j],t[2*j+1],"HPT"); hthres <<= 4;

         t[ip] = lthres | hthres;

         ip++;
      }
      
      rank /= 2; 
   }
   UChar_t sLpt, sHpt, lsLpt, lsHpt, usLpt, usHpt;
   sLpt  = ((t[0] & 0xC)  != 0);
   sHpt  = ((t[0] & 0xC0) != 0);
   lsLpt = ((t[0] & 0x1)  != 0);
   lsHpt = ((t[0] & 0x10) != 0);
   usLpt = ((t[0] & 0x2 ) != 0);
   usHpt = ((t[0] & 0x20) != 0);

   sHpt  <<= 1;
   lsLpt <<= 2;
   lsHpt <<= 3;
   usLpt <<= 4;
   usHpt <<= 5;

   fResponse = sLpt | sHpt | lsLpt | lsHpt | usLpt |usHpt;
   

}

//___________________________________________
UShort_t AliMUONGlobalTriggerBoard::Algo(UShort_t i, UShort_t j, const char *thres)
{
/// global trigger algorithm
///   a ,b = reg  response  =  Hpt (+|-|us|ls) |  Lpt (+|-|us|ls)  
                           
   TBits a(8), b(8); a.Set(8,&i); b.Set(8,&j);

   TBits trg1(2), trg2(2), trg(2);

   if (!strcmp(thres,"LPT"))
   {
      trg1[0] = a[2]; trg1[1] = a[3]; 
      trg2[0] = b[2]; trg2[1] = b[3];
   }
   else
   {
      trg1[0] = a[6]; trg1[1] = a[7]; 
      trg2[0] = b[6]; trg2[1] = b[7];         
   }
       
   TBits trgLS1(1), trgUS1(1), trgLS2(1), trgUS2(1), trgLS(1), trgUS(1);

   if (!strcmp(thres,"LPT"))
   {
      trgLS1[0] = a[0]; trgUS1[0] = a[1]; 
      trgLS2[0] = b[0]; trgUS2[0] = b[1];
   }
   else
   {
      trgLS1[0] = a[4]; trgUS1[0] = a[5]; 
      trgLS2[0] = b[4]; trgUS2[0] = b[5];         
   }

   trgLS[0] = ( trg1[0] & trg2[0] ) | ( trg1[1] & trg2[1] ) | trgLS1[0] | trgLS2[0];
   trgUS[0] = ( trg1[0] & trg2[1] ) | ( trg1[1] & trg2[0] ) | trgUS1[0] | trgUS2[0];
   
   trg[0] = trg1[0] | trg2[0];
   trg[1] = trg1[1] | trg2[1];
   
   TBits v(4);
   
   v[0] = trgLS[0];
   v[1] = trgUS[0];
   v[2] = trg[0];
   v[3] = trg[1];

   UShort_t rv = 0;
   v.Get(&rv);
   
   return rv;
}

//___________________________________________
void AliMUONGlobalTriggerBoard::BuildGlobalInput()
{
  /// build the 4 words (32bits) global input from the regional responses
  /// the order of regional responses is:
  /// 1R, 2R, 2-3R, 3R, 4R, 5R, 6R, 7R, 1L, 2L, 2-3L, 3L, 4L, 5L, 6L, 7L

  for (Int_t i=0;i< 4;i++) fGlobalInput[i] = 0;

  UShort_t regRespInv;
  TBits rs(8), rsi(8);
  for (Int_t iReg = 0; iReg < 16; iReg++) {

    // invert bit in regional response
    rs.Set(8,&fRegionalResponse[iReg]);
    for (Int_t i = 0; i < 4; i++) {
      rsi[2*i]   = rs[2*i+1];
      rsi[2*i+1] = rs[2*i];
    }
    regRespInv = 0;
    rsi.Get(&regRespInv);

    if (iReg < 8) {    // right
      // Lpt word
      fGlobalInput[0] |=  (regRespInv & 0x0F)       << (4*iReg);
      // Hpt word
      fGlobalInput[2] |= ((regRespInv & 0xF0) >> 4) << (4*iReg);
    } else {           // left
      // Lpt word
      fGlobalInput[1] |=  (regRespInv & 0x0F)       << (4*(iReg-8));
      // Hpt word
      fGlobalInput[3] |= ((regRespInv & 0xF0) >> 4) << (4*(iReg-8));
    }

  }

}

//___________________________________________
void AliMUONGlobalTriggerBoard::MaskGlobalInput()
{
  /// Apply masks to global input and recalculate regional inputs before
  /// applying the global response

  UShort_t regRespInv;
  TBits rs(8), rsi(8);

  // global input with masks applied
  UInt_t gitmp[4];

  for (Int_t i = 0; i < 4; i++) {
    gitmp[i] = fGlobalInput[i];
    gitmp[i] &= fMask[i];
  }

  for (Int_t iReg = 0; iReg < 16; iReg++) {
    fRegionalResponse[iReg] = 0;
    if (iReg < 8) {    // right
      // Lpt
      fRegionalResponse[iReg] |=  (gitmp[0] >> (4*iReg))     & 0xF;
      // Hpt
      fRegionalResponse[iReg] |= ((gitmp[2] >> (4*iReg))     & 0xF) << 4;
    } else {           // left
      // Lpt
      fRegionalResponse[iReg] |=  (gitmp[1] >> (4*(iReg-8))) & 0xF;
      // Hpt
      fRegionalResponse[iReg] |= ((gitmp[3] >> (4*(iReg-8))) & 0xF) << 4;
    }
    // invert bit in regional response
    rs.Set(8,&fRegionalResponse[iReg]);
    for (Int_t i = 0; i < 4; i++) {
      rsi[2*i]   = rs[2*i+1];
      rsi[2*i+1] = rs[2*i];
    }
    regRespInv = 0;
    rsi.Get(&regRespInv);
    fRegionalResponse[iReg] = regRespInv;
  }

}

//___________________________________________
void AliMUONGlobalTriggerBoard::Scan(Option_t*) const
{
  /// print global trigger output 
  TBits w(6); w.Set(6,&fResponse);

// TRG[1:0]
// 00 noth
// 01 negative track
// 10 positive track
// 11 undef

   Int_t iS[2] = {0,0};

   iS[0] = (Int_t)w.TestBitNumber(0);
   iS[1] = (Int_t)w.TestBitNumber(1);

   Int_t iPU[2] = {w[4],w[5]};
   Int_t iPL[2] = {w[2],w[3]};

   printf("============================================\n");
   printf(" Global Trigger output       Low pt  High pt\n");
   printf(" number of Single           :\t");
   for (Int_t i=0; i<2; i++) printf("%i\t",iS[i]);
   printf("\n");
   printf(" number of UnlikeSign pair  :\t"); 
   for (Int_t i=0; i<2; i++) printf("%i\t",iPU[i]);
   printf("\n");
   printf(" number of LikeSign pair    :\t");  
   for (Int_t i=0; i<2; i++) printf("%i\t",iPL[i]);
   printf("\n");
   printf("===================================================\n");
   printf("\n");
}

