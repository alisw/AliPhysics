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

//*-- Author: Rachid Guernane (LPCCFd)

#include "AliMUONGlobalTriggerBoard.h"

#include "TBits.h"

#include <Riostream.h>

ClassImp(AliMUONGlobalTriggerBoard)

//___________________________________________
AliMUONGlobalTriggerBoard::AliMUONGlobalTriggerBoard()
{
   for (Int_t i=0;i<16;i++) fRegionalResponse[i] = 0;
}

//___________________________________________
AliMUONGlobalTriggerBoard::AliMUONGlobalTriggerBoard(const char *name, Int_t a) : AliMUONTriggerBoard(name, a)
{
   for (Int_t i=0;i<16;i++) fRegionalResponse[i] = 0;
}

//___________________________________________
void AliMUONGlobalTriggerBoard::Response()
{
   Int_t t[16];

   for (Int_t i=0;i<16;i++) t[i] = fRegionalResponse[i];

   Int_t rank = 8;

   for (Int_t i=0;i<4;i++)
   {
      Int_t ip = 0;
      
      for (Int_t j=0;j<rank;j++)
      {
         UShort_t athres = Algo(t[2*j],t[2*j+1],"APT");

         UShort_t lthres = Algo(t[2*j],t[2*j+1],"LPT"); lthres <<= 4;

         UShort_t hthres = Algo(t[2*j],t[2*j+1],"HPT"); hthres <<= 8;

         t[ip] = athres | lthres | hthres;

         ip++;
      }
      
      rank /= 2; 
   }

   fResponse = t[0]; // 12-bit [H4:L4:A4]
}

//___________________________________________
UShort_t AliMUONGlobalTriggerBoard::Algo(UShort_t i, UShort_t j, char *thres)
{
   TBits a(12), b(12); a.Set(12,&i); b.Set(12,&j);

   TBits trg1(2), trg2(2), trg(2);

   if (!strcmp(thres,"APT"))
   {
      trg1[0] = a[2]; trg1[1] = a[3]; 
      trg2[0] = b[2]; trg2[1] = b[3];
   }
   else if (!strcmp(thres,"LPT"))
   {
      trg1[0] = a[6]; trg1[1] = a[7]; 
      trg2[0] = b[6]; trg2[1] = b[7];
   }
   else
   {
      trg1[0] = a[10]; trg1[1] = a[11]; 
      trg2[0] = b[10]; trg2[1] = b[11];         
   }
       
   TBits trgLS1(1), trgUS1(1), trgLS2(1), trgUS2(1), trgLS(1), trgUS(1);

   if (!strcmp(thres,"APT"))
   {
      trgLS1[0] = a[1]; trgUS1[0] = a[0]; 
      trgLS2[0] = b[1]; trgUS2[0] = b[0];
   }
   else if (!strcmp(thres,"LPT"))
   {
      trgLS1[0] = a[5]; trgUS1[0] = a[4]; 
      trgLS2[0] = b[5]; trgUS2[0] = b[4];
   }
   else
   {
      trgLS1[0] = a[9]; trgUS1[0] = a[8]; 
      trgLS2[0] = b[9]; trgUS2[0] = b[8];         
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
void AliMUONGlobalTriggerBoard::Scan(Option_t*)
{
   TBits w(12); w.Set(12,&fResponse);

// TRG[1:0]
// 00 noth
// 01 negative track
// 10 positive track
// 11 undef

   Int_t iSP[3] = {0,0,0}, iSM[3] = {0,0,0}, iSU[3] = {0,0,0};

   TBits a(2), n(2), p(2), u(2);
   
   UShort_t val;

   val = 1; n.Set(2,&val);
   val = 2; p.Set(2,&val);
   val = 3; u.Set(2,&val);
   
   a[0] = w[2];
   a[1] = w[3];
   
   if      (a==p) iSP[0] = 1;
   else if (a==n) iSM[0] = 1;
   else if (a==u) iSU[0] = 1;   

   a[0] = w[6];
   a[1] = w[7];

   if      (a==p) iSP[1] = 1;
   else if (a==n) iSM[1] = 1;
   else if (a==u) iSU[1] = 1;
   
   a[0] = w[10];
   a[1] = w[11];

   if      (a==p) iSP[2] = 1;
   else if (a==n) iSM[2] = 1;
   else if (a==u) iSU[2] = 1;

   Int_t iPU[3] = {w[0],w[4],w[8]};
   Int_t iPL[3] = {w[1],w[5],w[9]};

   printf("===================================================\n");
   printf(" Global Trigger output       Low pt  High pt   All\n");
   printf(" number of Single Plus      :\t");
   for (Int_t i=0; i<3; i++) printf("%i\t",iSP[i]);
   printf("\n");
   printf(" number of Single Minus     :\t");
   for (Int_t i=0; i<3; i++) printf("%i\t",iSM[i]);
   printf("\n");
   printf(" number of Single Undefined :\t"); 
   for (Int_t i=0; i<3; i++) printf("%i\t",iSU[i]);
   printf("\n");
   printf(" number of UnlikeSign pair  :\t"); 
   for (Int_t i=0; i<3; i++) printf("%i\t",iPU[i]);
   printf("\n");
   printf(" number of LikeSign pair    :\t");  
   for (Int_t i=0; i<3; i++) printf("%i\t",iPL[i]);
   printf("\n");
   printf("===================================================\n");
   printf("\n");
}

ClassImp(AliMUONGlobalTriggerBoard)
