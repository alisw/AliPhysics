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

//_________________________________________________________________________
// Implementation of local trigger board objects
// A local trigger board has as input a bit pattern and returns 
// the local trigger response after comparison w/ a LUT
//*-- Author: Rachid Guernane (LPCCFd)
//*
//*

#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONTriggerLut.h"

#include "AliLog.h"

#include <TBits.h>
#include <Riostream.h>

const Int_t AliMUONLocalTriggerBoard::fgkCircuitId[234] = 
{
  111,  121,  131,  141,  151,  161,  171,
  211,  212,  221,  222,  231,  232,  241,  242,  251,  252,  261,  262,  271,
  311,  312,  321,  322,  331,  332,  341,  342,  351,  352,  361,  362,  371,
  411,  412,  413,  421,  422,  423,  424,  431,  432,  433,  434,  441,  442,  451,  452,  461,  462,  471,
  521,  522,  523,  524,  531,  532,  533,  534,  541,  542,  551,  552,  561,  562,  571, 
  611,  612,  613,  621,  622,  623,  624,  631,  632,  633,  634,  641,  642,  651,  652,  661,  662,  671,
  711,  712,  721,  722,  731,  732,  741,  742,  751,  752,  761,  762,  771,
  811,  812,  821,  822,  831,  832,  841,  842,  851,  852,  861,  862,  871,
  911,  921,  931,  941,  951,  961,  971,
  -111, -121, -131, -141, -151, -161, -171,
  -211, -212, -221, -222, -231, -232, -241, -242, -251, -252, -261, -262, -271,
  -311, -312, -321, -322, -331, -332, -341, -342, -351, -352, -361, -362, -371,
  -411, -412, -413, -421, -422, -423, -424, -431, -432, -433, -434, -441, -442, -451, -452, -461, -462, -471,
  -521, -522, -523, -524, -531, -532, -533, -534, -541, -542, -551, -552, -561, -562, -571, 
  -611, -612, -613, -621, -622, -623, -624, -631, -632, -633, -634, -641, -642, -651, -652, -661, -662, -671,
  -711, -712, -721, -722, -731, -732, -741, -742, -751, -752, -761, -762, -771,
  -811, -812, -821, -822, -831, -832, -841, -842, -851, -852, -861, -862, -871,
  -911, -921, -931, -941, -951, -961, -971 
};

//___________________________________________
AliMUONLocalTriggerBoard::AliMUONLocalTriggerBoard()
    : AliMUONTriggerBoard(),
      fNumber(0),
      fCrate(0),
      fTC(kTRUE),
      fStripX11(0),
      fStripY11(0),
      fDev(0),
      fOutput(0),
      fLUT(0x0),
      fCoinc44(0)      
{
//* constructor
//*

   for (Int_t i=0; i<2; i++) 
      for (Int_t j=0; j<4; j++) 
      {
         fXY[i][j] = fXYU[i][j] = fXYD[i][j] = 0;

         fMask[i][j] = 0xFFFF;
      }

   for (Int_t i=0; i<10; i++) fSwitch[i] = 0;

   for (Int_t i=0; i<5; i++) fMinDevStrip[i] = fMinDev[i] = fCoordY[i] = 0;

   for (Int_t i=0; i<2; i++) fLutLpt[i] = fLutHpt[i] = 0;
}

//___________________________________________
AliMUONLocalTriggerBoard::AliMUONLocalTriggerBoard(const char *name, Int_t a,
                                                   AliMUONTriggerLut* lut) 
    : AliMUONTriggerBoard(name, a),
      fNumber(0),
      fCrate(0),
      fTC(kTRUE),
      fStripX11(0),
      fStripY11(0),
      fDev(0),
      fOutput(0),
      fLUT(lut),
      fCoinc44(0)
{
//* constructor
//*
   
   for (Int_t i=0; i<2; i++) 
      for (Int_t j=0; j<4; j++) 
      {
         fXY[i][j] = fXYU[i][j] = fXYD[i][j] = 0;

         fMask[i][j] = 0xFFFF;
      }

   for (Int_t i=0; i<10; i++) fSwitch[i] = 0;

   for (Int_t i=0; i<5; i++) fMinDevStrip[i] = fMinDev[i] = fCoordY[i] = 0;

   for (Int_t i=0; i<2; i++) fLutLpt[i] = fLutHpt[i] = 0;
}

//______________________________________________________________________________
AliMUONLocalTriggerBoard::AliMUONLocalTriggerBoard(const AliMUONLocalTriggerBoard& right) 
    : AliMUONTriggerBoard(right),
      fNumber(right.fNumber),
      fCrate(right.fCrate),
      fTC(right.fTC),
      fStripX11(right.fStripX11),
      fStripY11(right.fStripY11),
      fDev(right.fDev),
      fOutput(right.fOutput),
      fLUT(right.fLUT),
      fCoinc44(right.fCoinc44)
{  
/// Protected copy constructor (not implemented)

  AliFatal("Copy constructor not provided.");
}

//______________________________________________________________________________
AliMUONLocalTriggerBoard& 
AliMUONLocalTriggerBoard::operator=(const AliMUONLocalTriggerBoard& right)
{
/// Protected assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  AliFatal("Assignement operator not provided.");
    
  return *this;  
}    

//___________________________________________
void AliMUONLocalTriggerBoard::Reset()
{
//* reset board
//*
   for (Int_t i=0; i<2; i++) 
      for (Int_t j=0; j<4; j++) 
         fXY[i][j] = fXYU[i][j] = fXYD[i][j] = 0;

   fResponse = 0;

   for (Int_t i=0; i<5; i++) fMinDevStrip[i] = fMinDev[i] = fCoordY[i] = 0;

   fOutput = 0;
   
   fStripX11 = fStripY11 = fDev = 0;

   for (Int_t i=0; i<2; i++) fLutLpt[i] = fLutHpt[i] = 0;
}

//___________________________________________
void AliMUONLocalTriggerBoard::Setbit(Int_t strip, Int_t cathode, Int_t chamber)
{
// 0 .. LBS   :   N-1 .. MSB
   TBits w, m;

   UShort_t xy = fXY[cathode][chamber], mask = fMask[cathode][chamber];

   w.Set(16,&xy);
   m.Set(16,&mask);

   Int_t s = strip - int(strip / 16) * 16;

   w.SetBitNumber(s);
   
   w &= m;

   UShort_t value;

   w.Get(&value);

   fXY[cathode][chamber] = value;
}

//___________________________________________
void AliMUONLocalTriggerBoard::SetbitM(Int_t strip, Int_t cathode, Int_t chamber)
{
// 0 .. LBS   :   N-1 .. MSB
   TBits w, m;

   UShort_t xy = fXY[cathode][chamber], mask = fMask[cathode][chamber];

   w.Set(16,&xy);
   m.Set(16,&mask);

   w.SetBitNumber(strip);
   
   w &= m;

   UShort_t value;

   w.Get(&value);

   fXY[cathode][chamber] = value;
}

//___________________________________________
void AliMUONLocalTriggerBoard::Pattern(Option_t *option) const
{
//* print bit pattern
//*
   TString op = option;
   
   if (op.Contains("X")) BP("X");

   if (op.Contains("Y")) BP("Y");
}


//___________________________________________
void AliMUONLocalTriggerBoard::BP(Option_t *option) const
{
// RESPECT THE OLD PRINTOUT FORMAT
  
  const Int_t kModuleId[126] = 
  {11,12,13,14,15,16,17,         // right side of the chamber
  21,22,23,24,25,26,27,
  31,32,33,34,35,36,37,
  41,42,43,44,45,46,47,
  51,52,53,54,55,56,57,
  61,62,63,64,65,66,67,
  71,72,73,74,75,76,77,
  81,82,83,84,85,86,87,
  91,92,93,94,95,96,97,   
  -11,-12,-13,-14,-15,-16,-17,  // right side of chamber
  -21,-22,-23,-24,-25,-26,-27,
  -31,-32,-33,-34,-35,-36,-37,
  -41,-42,-43,-44,-45,-46,-47,
  -51,-52,-53,-54,-55,-56,-57,
  -61,-62,-63,-64,-65,-66,-67,
  -71,-72,-73,-74,-75,-76,-77,
  -81,-82,-83,-84,-85,-86,-87,
  -91,-92,-93,-94,-95,-96,-97};

  const Int_t kNstripY[126]=
  { 8, 8, 8, 8, 8, 8,16,  // right side of the chamber
  8, 8, 8, 8, 8, 8,16,
  16,16,16,16,16, 8,16,
  16,16,16,16,16, 8,16,
  0, 8,16,16,16, 8,16,
  16,16,16,16,16, 8,16,
  16,16,16,16,16, 8,16,
  8, 8, 8, 8, 8, 8,16,
  8, 8, 8, 8, 8, 8,16,  
  8, 8, 8, 8, 8, 8,16,  // left side of the chamber
  8, 8, 8, 8, 8, 8,16,
  16,16,16,16,16, 8,16,
  16,16,16,16,16, 8,16,
  0, 8,16,16,16, 8,16,
  16,16,16,16,16, 8,16,
  16,16,16,16,16, 8,16,
  8, 8, 8, 8, 8, 8,16,
  8, 8, 8, 8, 8, 8,16};

   TString op = option;
   	
   TString nn = GetName();

   if (op.Contains("X"))
   {
      printf("-------- TRIGGER INPUT ---------\n");
      printf("===============================================================\n");
      printf("                            5432109876543210");

      char *x[4] = {"XMC11","XMC12","XMC21","XMC22"};
      char *s[4] = {"                      ",
                    "                      ",
                    "              ",
                    "              "};
      
      for (Int_t ch=0; ch<4; ch++)
      { 
         printf("\n %s%s", x[ch], s[ch]);

         UShort_t xy = fXY[0][ch];

         TBits w(16); w.Set(16,&xy);

         if (ch<2) cout << w;
         else
         {
            UShort_t xyd = fXYD[0][ch], xyu = fXYU[0][ch];
            TBits dw(16), uw(16); dw.Set(16,&xyd); uw.Set(16,&xyu); 

            TBits ew(32);
         
            for (Int_t i=0;i<16;i++) ew[i+8] = w[i];

            for (Int_t i=0;i<8;i++) 
            {
               ew[i]    = dw[i+8]; // 8 MSB
               ew[i+24] = uw[i];   // 
            }

            cout << ew;
         }
      }
   
      printf("\n                    ");
      printf("10987654321098765432109876543210\n");
   }

   if (op.Contains("Y"))
   {
      printf("---------------------------------------------------------------\n");
      printf("                            ");

/*    OLD NUMBERING STYLE    */
/**/
      Int_t idCircuit = 0, absidModule = 0;

      if (!(nn.Contains("Int"))) 
      {	
	idCircuit   = fgkCircuitId[GetI()];
	absidModule = TMath::Abs(Int_t(idCircuit/10));
      }
		
      Int_t iModule=0;

      for (Int_t i=0; i<63; i++) 
      {
	if (kModuleId[i]==absidModule) 
         { 
            iModule=i;
            break;
         }
      }

      Int_t nStrip = kNstripY[iModule];
      for (Int_t istrip=nStrip-1; istrip>=0; istrip--) {
         if (istrip>9)  printf("%i",istrip-10*Int_t(istrip/10));
         if (istrip<10) printf("%i",istrip);
      }
/**/
/*                           */

      UShort_t xyval = 0;

      if (fSwitch[1])
      {
         xyval = fXY[1][0];
         TBits v11(8); v11.Set(8,&xyval);
         printf("\n YMC11                      ");
         cout << v11;         

         xyval = fXY[1][1];
         TBits v12(8); v12.Set(8,&xyval);
         printf("\n YMC12                      ");
         cout << v12;

         xyval = fXY[1][2];
         TBits v21(8); v21.Set(8,&xyval);
         printf("\n YMC21                      ");         
         cout << v21;

         xyval = fXY[1][3];
         TBits v22(8); v22.Set(8,&xyval);
         printf("\n YMC22                      ");
         cout << v22 << endl;
      }
      else
      {
         xyval = fXY[1][0];
         TBits v11(16); v11.Set(16,&xyval);
         printf("\n YMC11                      ");
         cout << v11;         

         xyval = fXY[1][1];
         TBits v12(16); v12.Set(16,&xyval);
         printf("\n YMC12                      ");
         cout << v12;

         xyval = fXY[1][2];
         TBits v21(16); v21.Set(16,&xyval);
         printf("\n YMC21                      ");         
         cout << v21;

         xyval = fXY[1][3];
         TBits v22(16); v22.Set(16,&xyval);
         printf("\n YMC22                      ");
         cout << v22 << endl;
      }

//    tmp
      printf("---------------------------------------------------------------");
      printf("\n upper part of circuit %i",idCircuit);
      printf("\n UMC21                      ");
      xyval = fXYU[1][2];
      TBits wu21(16); wu21.Set(16,&xyval);
      cout << wu21;
      printf("\n UMC22                      ");
      xyval = fXYU[1][3];
      TBits wu22(16); wu22.Set(16,&xyval);
      cout << wu22;
      printf("\n lower part of circuit %i",idCircuit);
      printf("\n LMC21                      ");
      xyval = fXYD[1][2];
      TBits wl21(16); wl21.Set(16,&xyval);
      cout << wl21;
      printf("\n LMC22                      ");
      xyval = fXYD[1][3];
      TBits wl22(16); wl22.Set(16,&xyval);
      cout << wl22;
      printf("\n");
      printf("===============================================================\n");
   }
}

//___________________________________________
void AliMUONLocalTriggerBoard::Conf() const
{
//* board switches
//*
   cout << "Switch(" << GetName() << ")" 
        << " x2d = "           << fSwitch[0] 
        << " x2m = "           << fSwitch[1] 
        << " x2u = "           << fSwitch[2] 
        << " OR[0] = "         << fSwitch[3] 
        << " OR[1] = "         << fSwitch[4] 
        << " EN-Y = "          << fSwitch[5] 
        << " ZERO-ALLY-LSB = " << fSwitch[6] 
        << " ZERO-down = "     << fSwitch[7] 
        << " ZERO-middle = "   << fSwitch[8] 
        << " ZERO-up = "       << fSwitch[9] 
        << " trans. conn. "    << fTC 
        << " Slot = "          << fSlot 
        << endl;
}

//___________________________________________
void AliMUONLocalTriggerBoard::Module(char *mod)
{
//* get module from name
//*
   const Int_t kMaxfields = 2; char **fields = new char*[kMaxfields];

   char s[100]; strcpy(s, GetName());

   Int_t numlines = 0;

   for (char *token = strtok(s, "B");
        token != NULL;
        token = strtok(NULL, " "))
   {
      fields[numlines] = new char[strlen(token)+1];
      strcpy(fields[numlines++],token);
   }
 
   strcpy(mod,fields[0]);
}

//___________________________________________
void AliMUONLocalTriggerBoard::TrigX(Int_t ch1q[16], Int_t ch2q[16], Int_t ch3q[32], Int_t ch4q[32])
{
// note : coinc44 = flag 0 or 1 (0 coincidence -> 3/4, 1 coincidence -> 4/4)
//---------------------------------------------------------
// step # 1 : declustering, reduction DS, calculate sgle & dble
//---------------------------------------------------------
   Int_t ch1e[19], ch2e[20], ch3e[35], ch4e[36]; 
   Int_t sgleHit1[31], sgleHit2[63];
   Int_t dbleHit1[31], dbleHit2[63];

   Int_t i;
   Int_t j;
   Int_t istrip;

   for (i=0; i<31; i++) {
      sgleHit1[i]=0;
      dbleHit1[i]=0;
   }
   for (i=0; i<63; i++) {
      sgleHit2[i]=0;
      dbleHit2[i]=0;
   }

//--- inititialize che using chq 
   for (i=0; i<19; i++) {
      if (i<1||i>16)  ch1e[i]=0; 
      else            ch1e[i]=ch1q[i-1]; 
   }
   for (i=0; i<20; i++) {
      if (i<2||i>17) ch2e[i]=0; 
      else           ch2e[i]=ch2q[i-2]; 
   }
   for (i=0; i<35; i++) {
      if (i<1||i>32) ch3e[i]=0; 
      else           ch3e[i]=ch3q[i-1];
   }
   for (i=0; i<36; i++) {
      if (i<2||i>33) ch4e[i]=0; 
      else           ch4e[i]=ch4q[i-2];
   }


//--- calculate dble & sgle first station
   for (i=0; i<=15; i++) {                   
      sgleHit1[2*i] = (!ch1e[i+1]|(ch1e[i]^ch1e[i+2])) & 
         (!ch2e[i+2] | (ch2e[i+1]^ch2e[i+3]));

      dbleHit1[2*i] = ch1e[i+1]&!(ch1e[i+2]^ch1e[i]) & 
         (ch2e[i+2] | (!ch2e[i]&ch2e[i+1]) | (ch2e[i+3]&!ch2e[i+4]));
   }

   for (i=0; i<=14; i++) {               
      sgleHit1[2*i+1] = (!ch1e[i+1]|!ch1e[i+2]|(ch1e[i]^ch1e[i+3])) & 
         (!ch2e[i+2] | !ch2e[i+3] | (ch2e[i+1]^ch2e[i+4]));
      dbleHit1[2*i+1] = ch1e[i+1]&ch1e[i+2]&!(ch1e[i]^ch1e[i+3]) & 
         (ch2e[i+2]&(!ch2e[i+1]|!ch2e[i]) | 
          ch2e[i+3]&(ch2e[i+2]|!ch2e[i+4]|!ch2e[i+5]));
   }

//--- calculate dble & sgle second station
   for (i=0; i<=31; i++) {               
      sgleHit2[2*i] = (!ch3e[i+1]|(ch3e[i]^ch3e[i+2])) & 
         (!ch4e[i+2] | (ch4e[i+1]^ch4e[i+3]));
      dbleHit2[2*i] = ch3e[i+1]&!(ch3e[i+2]^ch3e[i]) & 
         (ch4e[i+2] | (!ch4e[i]&ch4e[i+1]) | (ch4e[i+3]&!ch4e[i+4]));
   }
  
   for (i=0; i<=30; i++) {               
      sgleHit2[2*i+1] = (!ch3e[i+1]|!ch3e[i+2]|(ch3e[i]^ch3e[i+3])) & 
         (!ch4e[i+2] | !ch4e[i+3] | (ch4e[i+1]^ch4e[i+4]));
      dbleHit2[2*i+1] = ch3e[i+1]&ch3e[i+2]&!(ch3e[i]^ch3e[i+3]) & 
         (ch4e[i+2]&(!ch4e[i+1]|!ch4e[i]) | 
          ch4e[i+3]&(ch4e[i+2]|!ch4e[i+4]|!ch4e[i+5]));
   }

//--- 
   if(AliDebugLevel()==3||AliDebugLevel()==5) {
      printf("===============================================================\n");
      printf(" X plane after sgle and dble \n");
      printf("                       0987654321098765432109876543210");
      printf("\n SGLE1                 ");
      for (istrip=30; istrip>=0; istrip--) printf("%i",(!sgleHit1[istrip]));
      printf("\n DBLE1                 ");
      for (istrip=30; istrip>=0; istrip--) printf("%i",dbleHit1[istrip]);
      printf("\n SGLE2 ");
      for (istrip=62; istrip>=0; istrip--) printf("%i",(!sgleHit2[istrip]));
      printf("\n DBLE2 ");
      for (istrip=62; istrip>=0; istrip--) printf("%i",dbleHit2[istrip]);
      printf("\n       210987654321098765432109876543210987654321098765432109876543210\n");
   }
  
//---------------------------------------------------------
// step # 2 : coincidence 3/4
//---------------------------------------------------------
   Int_t rearImage[31][31];
   for (i=0; i<31; i++) {
      for (j=0; j<31; j++) {
         rearImage[i][j]=0;
      }
   }

   Int_t notOr1=!dbleHit1[30] & !dbleHit1[29] & !dbleHit1[28] & !dbleHit1[27] & 
      !dbleHit1[26] & !dbleHit1[25] & !dbleHit1[24] & !dbleHit1[23] &
      !dbleHit1[22] & !dbleHit1[21] & !dbleHit1[20] & !dbleHit1[19] & 
      !dbleHit1[18] & !dbleHit1[17] & !dbleHit1[16] & !dbleHit1[15] & 
      !dbleHit1[14] & !dbleHit1[13] & !dbleHit1[12] & !dbleHit1[11] & 
      !dbleHit1[10] & !dbleHit1[9]  & !dbleHit1[8]  & !dbleHit1[7]  & 
      !dbleHit1[6]  & !dbleHit1[5]  & !dbleHit1[4]  & !dbleHit1[3]  & 
      !dbleHit1[2]  & !dbleHit1[1]  & !dbleHit1[0]  & !fCoinc44;

   Int_t notOr2= !dbleHit2[62] & !dbleHit2[61] & !dbleHit2[60] & !dbleHit2[59] & 
      !dbleHit2[58] & !dbleHit2[57] & !dbleHit2[56] & !dbleHit2[55] & 
      !dbleHit2[54] & !dbleHit2[53] & !dbleHit2[52] & !dbleHit2[51] & 
      !dbleHit2[50] & !dbleHit2[49] & !dbleHit2[48] & !dbleHit2[47] & 
      !dbleHit2[46] & !dbleHit2[45] & !dbleHit2[44] & !dbleHit2[43] & 
      !dbleHit2[42] & !dbleHit2[41] & !dbleHit2[40] & !dbleHit2[39] & 
      !dbleHit2[38] & !dbleHit2[37] & !dbleHit2[36] & !dbleHit2[35] & 
      !dbleHit2[34] & !dbleHit2[33] & !dbleHit2[32] & !dbleHit2[31] &
      !dbleHit2[30] & !dbleHit2[29] & !dbleHit2[28] & !dbleHit2[27] & 
      !dbleHit2[26] & !dbleHit2[25] & !dbleHit2[24] & !dbleHit2[23] & 
      !dbleHit2[22] & !dbleHit2[21] & !dbleHit2[20] & !dbleHit2[19] & 
      !dbleHit2[18] & !dbleHit2[17] & !dbleHit2[16] & !dbleHit2[15] & 
      !dbleHit2[14] & !dbleHit2[13] & !dbleHit2[12] & !dbleHit2[11] & 
      !dbleHit2[10] & !dbleHit2[9]  & !dbleHit2[8]  & !dbleHit2[7]  & 
      !dbleHit2[6]  & !dbleHit2[5]  & !dbleHit2[4]  & !dbleHit2[3]  & 
      !dbleHit2[2]  & !dbleHit2[1]  & !dbleHit2[0]  & !fCoinc44;	

// DS reduction
   for (i=0; i<31; i++) {
      sgleHit1[i] = !sgleHit1[i]&notOr1;
   }
   for (i=0; i<63; i++) {
      sgleHit2[i] = !sgleHit2[i]&notOr2;
   }

// extract rearImage
   for (i=0; i<31; i++){
      Int_t tmpSgleHit2[31];
      Int_t tmpDbleHit2[31];
      for (j=0; j<31; j++){
         tmpSgleHit2[j] = sgleHit2[i+j+1];
         tmpDbleHit2[j] = dbleHit2[i+j+1];
      }

      for (Int_t k=0; k<31; k++) {
         rearImage[i][k]=(sgleHit1[i]&tmpDbleHit2[k])|
            (dbleHit1[i]&(tmpSgleHit2[k]|tmpDbleHit2[k]));
      }
   }

   //-----------
   if(AliDebugLevel()==3||AliDebugLevel()==5) {
      printf("===============================================================\n");
      for (i=30; i>=0; i--) {
         printf("%i \t",i);
         for (istrip=31; istrip>=0; istrip--) printf("%i",rearImage[i][istrip]);
         printf("\n");   
      }
   }

//---------------------------------------------------------
// step # 3 : calculate deviation
//--------------------------------------------------------- 
   Int_t dev[31][6];
   for (i=0; i<31; i++) {
      for (j=0; j<6; j++) {
         dev[i][j]=0;
      }
   }

   for (i=0; i<31; i++){
      Int_t leftDev[5], rightDev[5]; 
      Int_t orL1, andL1, andL2, orR1, orR2, andR1, andR2, andR3;

// calculate Left deviation
      orL1=rearImage[i][16]|rearImage[i][18]|rearImage[i][20]|rearImage[i][22];
      andL1=!rearImage[i][17]&!rearImage[i][19]&!rearImage[i][21] & !orL1; 
      andL2=!rearImage[i][23]&!rearImage[i][24]&!rearImage[i][25]&!rearImage[i][26];
 
      leftDev[0] = (rearImage[i][16]|!rearImage[i][17]) & 
         (rearImage[i][16]|rearImage[i][18]|!rearImage[i][19]&
          (rearImage[i][20]|!rearImage[i][21])) &
         (orL1|!rearImage[i][23]&(rearImage[i][24]|!rearImage[i][25])) & 
         (orL1|rearImage[i][24]|rearImage[i][26]|!rearImage[i][27]&
          (rearImage[i][28]|!rearImage[i][29]));
				
      leftDev[1] = !rearImage[i][16] & 
         !(!rearImage[i][17]&!rearImage[i][18]&!rearImage[i][21]&!rearImage[i][22] & 
           (!rearImage[i][25]&!rearImage[i][26]&(rearImage[i][27]|rearImage[i][28]))) &
         (rearImage[i][17]|rearImage[i][18] | !rearImage[i][19]&!rearImage[i][20]) &
         (rearImage[i][17]|rearImage[i][18]|rearImage[i][21]|rearImage[i][22] | 
          !rearImage[i][23]&!rearImage[i][24]);
				
      leftDev[2] = (!rearImage[i][16]&!rearImage[i][17]&!rearImage[i][18]) & 
         (rearImage[i][19]|rearImage[i][20]|rearImage[i][21]|rearImage[i][22] | andL2);
		
      leftDev[3] = andL1;
		
      leftDev[4] = 
         !rearImage[i][27]&!rearImage[i][28]&!rearImage[i][29]&!rearImage[i][30] & 
         andL1 & andL2;

      // calculate Right deviation
      orR1=rearImage[i][8]|rearImage[i][10]|rearImage[i][12]|rearImage[i][14];
      orR2=rearImage[i][8]|rearImage[i][9]|rearImage[i][10]|rearImage[i][11];
      andR1=!rearImage[i][12]&!rearImage[i][13]&!rearImage[i][14]&!rearImage[i][15];
      andR2=
         !rearImage[i][8]&!rearImage[i][9]&!rearImage[i][10]&!rearImage[i][11] & andR1;
      andR3=!rearImage[i][4]&!rearImage[i][5]&!rearImage[i][6]&!rearImage[i][7]; 
		
      rightDev[0] = !rearImage[i][15]&(rearImage[i][14]|!rearImage[i][13]) & 
         ((rearImage[i][12]|rearImage[i][14]|!rearImage[i][11]&
           (rearImage[i][10]|!rearImage[i][9])) &
          ((orR1|!rearImage[i][7]&(rearImage[i][6]|!rearImage[i][5])) & 
           (orR1|rearImage[i][4]|rearImage[i][6]|!rearImage[i][3]&(rearImage[i][2]|
                                                                   !rearImage[i][1]))));
				
      rightDev[1] = !rearImage[i][15]&!rearImage[i][14] & 
         !(!rearImage[i][4]&!rearImage[i][5]&!rearImage[i][8]&!rearImage[i][9] &
           (!rearImage[i][12]&!rearImage[i][13]&(rearImage[i][2]|rearImage[i][3]))) &
         (rearImage[i][12]|rearImage[i][13] | !rearImage[i][10]&!rearImage[i][11]) & 
         (rearImage[i][8]|rearImage[i][9]|rearImage[i][12]|rearImage[i][13] | 
          !rearImage[i][6]&!rearImage[i][7]);
		
      rightDev[2] = andR1 & (orR2 | andR3); 
      rightDev[3] = andR2;		
      rightDev[4] = 
         !rearImage[i][0]&!rearImage[i][1]&!rearImage[i][2]&!rearImage[i][3] & 
         andR2 & andR3 ;

      // compare Left & Right deviations
      Int_t tmpLeftDev=0, tmpRightDev=0;
      for (j=0; j<5; j++){
         tmpLeftDev  = tmpLeftDev + Int_t(leftDev[j]<<j); 
         tmpRightDev = tmpRightDev + Int_t(rightDev[j]<<j); 
      }

      // assign mimimum deviation do dev[][]
      if (tmpLeftDev < tmpRightDev ){
         for (j=0; j<5; j++){ dev[i][j]=leftDev[j];}
         dev[i][5]=1;
      } else {
         for (j=0; j<5; j++){ dev[i][j]=rightDev[j];}
         dev[i][5]=0;
      }
   }
  
//---
   if(AliDebugLevel()==3||AliDebugLevel()==5) {
      printf("===============================================================\n");
      for (i=30; i>=0; i--) {
         printf("%i \t",i);
         for (istrip=5; istrip>=0; istrip--) printf("%i",dev[i][istrip]);
         printf(" \n");
      }
   }

//---------------------------------------------------------
// step # 4 : sort deviation
//--------------------------------------------------------- 
   Int_t bga1[16], bga2[8], bga3[4], bga4[2], bga5;
   Int_t tmpbga1[16][6], tmpbga2[8][6], tmpbga3[4][6], tmpbga4[2][6], tmpbga5[6];
   Int_t tmpMax[6]={1,1,1,1,1,0};

   for (i=0; i<15; i++) {
      Sort2x5(dev[2*i],dev[2*i+1],tmpbga1[i],bga1[i]);
   }  
   Sort2x5(dev[30],tmpMax,tmpbga1[15],bga1[15]);

//--    
   if(AliDebugLevel()==3||AliDebugLevel()==5) {
      printf("===============================================================\n");
      printf(" sorting : 1st level \n");
      for (i=15; i>=0; i--) {
         printf("\t %i \t",bga1[i]); 	
         for (j=5; j>=0; j--) printf("%i",tmpbga1[i][j]); 
         printf(" \n");
      }
   }

   for (i=0; i<8; i++) {  
      Sort2x5(tmpbga1[2*i],tmpbga1[2*i+1],tmpbga2[i],bga2[i]);
   }

//--    
   if(AliDebugLevel()==3||AliDebugLevel()==5) {
      printf("===============================================================\n");
      printf(" sorting : 2nd level \n");
      for (i=7; i>=0; i--) {
         printf("\t %i \t",bga2[i]); 	
         for (j=5; j>=0; j--) printf("%i",tmpbga1[i][j]); 	
         printf(" \n");
      }
   }
  
   for (i=0; i<4; i++) {  
      Sort2x5(tmpbga2[2*i],tmpbga2[2*i+1],tmpbga3[i],bga3[i]);
   }

//--    
   if(AliDebugLevel()==3||AliDebugLevel()==5) {
      printf("===============================================================\n");
      printf(" sorting : 3rd level \n");
      for (i=3; i>=0; i--) {
         printf("\t %i \t",bga3[i]); 	
         for (j=5; j>=0; j--) printf("%i",tmpbga3[i][j]); 
         printf(" \n");
      }
   }

   for (i=0; i<2; i++) {  
      Sort2x5(tmpbga3[2*i],tmpbga3[2*i+1],tmpbga4[i],bga4[i]);
   }

//--    
   if(AliDebugLevel()==3||AliDebugLevel()==5) {
      printf("===============================================================\n");
      printf(" sorting : 4th level \n");
      for (i=1; i>=0; i--) {
         printf("\t %i \t",bga4[i]); 	
         for (j=5; j>=0; j--) printf("%i",tmpbga4[i][j]);
         printf(" \n");
      }
   }
  
   Sort2x5(tmpbga4[0],tmpbga4[1],tmpbga5,bga5);

   // coding from 6 to 5 bits 
   fMinDev[4] = tmpbga5[5] | tmpbga5[4];
   for (i=0; i<4; i++) { 
      fMinDev[i]=tmpbga5[i] & !tmpbga5[4];
   }

   // find address of strip with minimum deviation 
   fMinDevStrip[4]=bga5;
   if (bga5<=1) fMinDevStrip[3]=bga4[bga5];

   Int_t tmpAd=fMinDevStrip[3]+fMinDevStrip[4]*2;
   if (tmpAd<=3) fMinDevStrip[2]=bga3[tmpAd];

   tmpAd=fMinDevStrip[2]+fMinDevStrip[3]*2+fMinDevStrip[4]*4;
   if (tmpAd<=7) fMinDevStrip[1]=bga2[tmpAd];

   tmpAd=fMinDevStrip[1]+fMinDevStrip[2]*2+fMinDevStrip[3]*4+fMinDevStrip[4]*8;
   if (tmpAd<=15) fMinDevStrip[0]=bga1[tmpAd];

   if(AliDebugLevel()==3||AliDebugLevel()==5) {
      printf("===============================================================\n");
      printf("minDevStrip = ");
      for  (i=4; i>=0; i--) printf("%i",fMinDevStrip[i]);
      printf(" minDev = ");
      for  (i=4; i>=0; i--) printf("%i",fMinDev[i]); 
      printf(" \n");
      printf("===============================================================\n");
   }

}

//___________________________________________
void AliMUONLocalTriggerBoard::Sort2x5(Int_t dev1[6], Int_t dev2[6],
                                       Int_t minDev[6], Int_t &dev1GTdev2)
{ 
// returns minimun between dev1 and dev2
   Int_t tmpDev1=0, tmpDev2=0;

   for (Int_t j=0; j<5; j++)
   {
      tmpDev1 += Int_t(dev1[j]<<j); 
      tmpDev2 += Int_t(dev2[j]<<j); 
   }

   if (tmpDev1<=tmpDev2)
   {
      for (Int_t j=0; j<=5; j++) minDev[j]=dev1[j];
      dev1GTdev2=0;
   } 
   else 
   {
      for (Int_t j=0; j<=5; j++) minDev[j]=dev2[j];
      dev1GTdev2=1;   
   }
}

//___________________________________________
void AliMUONLocalTriggerBoard::TrigY(Int_t y1[16], Int_t y2[16], Int_t y3[16], Int_t y4[16],
                                     Int_t y3u[16], Int_t y3d[16], Int_t y4u[16], Int_t y4d[16])
{
// note : resMid = 1 -> cancel 
//---------------------------------------------------------
// step # 1 : prehandling Y
//--------------------------------------------------------- 
   Int_t i;
   Int_t istrip;

   for (i=0; i<16; i++)
   {
      y3[i]=y3[i]&!fSwitch[8];
      y4[i]=y4[i]&!fSwitch[8];
   }

// 10/29/04 fZeroAllYLSB added
//    for (i=0; i<8; i++)
//    {
//       y1[i] = y1[i]&!fSwitch[6];     
//       y2[i] = y2[i]&!fSwitch[6];      
//       y3[i] = y3[i]&!fSwitch[6];      
//       y4[i] = y4[i]&!fSwitch[6];
//    }

   Int_t ch1[16], ch2[16], ch3[16], ch4[16];

   Int_t tmpy3to16[16], tmpy4to16[16];
   Int_t tmpy3uto16[16], tmpy3dto16[16], tmpy4uto16[16], tmpy4dto16[16];
   for (i=0; i<8; i++){
      ch1[2*i]   = y1[i]&fSwitch[1] | y1[2*i]&!fSwitch[1];		
      ch1[2*i+1] = y1[i]&fSwitch[1] | y1[2*i+1]&!fSwitch[1];

      ch2[2*i]   = y2[i]&fSwitch[1] | y2[2*i]&!fSwitch[1];		
      ch2[2*i+1] = y2[i]&fSwitch[1] | y2[2*i+1]&!fSwitch[1];

      tmpy3to16[2*i  ] = y3[i]&fSwitch[1] | y3[2*i  ]&!fSwitch[1];		
      tmpy3to16[2*i+1] = y3[i]&fSwitch[1] | y3[2*i+1]&!fSwitch[1];

      tmpy4to16[2*i  ] = y4[i]&fSwitch[1] | y4[2*i  ]&!fSwitch[1];
      tmpy4to16[2*i+1] = y4[i]&fSwitch[1] | y4[2*i+1]&!fSwitch[1];

      tmpy3uto16[2*i  ] = y3u[i]&fSwitch[2] | y3u[2*i  ]&!fSwitch[2]; 
      tmpy3uto16[2*i+1] = y3u[i]&fSwitch[2] | y3u[2*i+1]&!fSwitch[2];

      tmpy4uto16[2*i  ] = y4u[i]&fSwitch[2] | y4u[2*i  ]&!fSwitch[2]; 
      tmpy4uto16[2*i+1] = y4u[i]&fSwitch[2] | y4u[2*i+1]&!fSwitch[2];

      tmpy3dto16[2*i  ] = y3d[i]&fSwitch[0] | y3d[2*i  ]&!fSwitch[0]; 
      tmpy3dto16[2*i+1] = y3d[i]&fSwitch[0] | y3d[2*i+1]&!fSwitch[0];
    
      tmpy4dto16[2*i  ] = y4d[i]&fSwitch[0] | y4d[2*i  ]&!fSwitch[0]; 
      tmpy4dto16[2*i+1] = y4d[i]&fSwitch[0] | y4d[2*i+1]&!fSwitch[0];
   }
  
   if (fSwitch[3]==0&&fSwitch[4]==0){
      for (i=0; i<16; i++){
         ch3[i] = tmpy3to16[i];
         ch4[i] = tmpy4to16[i];
      }
   }
   if (fSwitch[3]==0&&fSwitch[4]==1){
      for (i=0; i<16; i++){
         ch3[i] = tmpy3dto16[i]|tmpy3to16[i];
         ch4[i] = tmpy4dto16[i]|tmpy4to16[i];
      }
   }
   if (fSwitch[3]==1&&fSwitch[4]==0){
      for (i=0; i<16; i++){
         ch3[i] = tmpy3uto16[i]|tmpy3to16[i];
         ch4[i] = tmpy4uto16[i]|tmpy4to16[i];
      }
   }
   if (fSwitch[3]==1&&fSwitch[4]==1){
      for (i=0; i<16; i++){
         ch3[i] = tmpy3dto16[i]|tmpy3to16[i]|tmpy3uto16[i];
         ch4[i] = tmpy4dto16[i]|tmpy4to16[i]|tmpy4uto16[i];
      }
   }

// debug
   if(AliDebugLevel()==4||AliDebugLevel()==5) {
      printf("===============================================================\n");  
      printf(" Y plane after PreHandling x2m x2u x2d orMud %i %i %i %i %i \n",
             fSwitch[1],fSwitch[2], fSwitch[0],fSwitch[3],fSwitch[4]);
      printf("                            ");
      for (istrip=15; istrip>=0; istrip--) {
         if (istrip>9)  printf("%i",istrip-10*Int_t(istrip/10));
         if (istrip<10) printf("%i",istrip);
      }  
      printf("\n YMC11                      ");
      for (istrip=15; istrip>=0; istrip--) printf("%i",ch1[istrip]); 
      printf("\n YMC12                      ");
      for (istrip=15; istrip>=0; istrip--) printf("%i",ch2[istrip]); 
      printf("\n YMC21                      ");
      for (istrip=15; istrip>=0; istrip--) printf("%i",ch3[istrip]); 
      printf("\n YMC22                      ");
      for (istrip=15; istrip>=0; istrip--) printf("%i",ch4[istrip]); 
      printf(" \n"); 
   }
//debug
  
//---------------------------------------------------------
// step # 2 : calculate sgle and dble, apply DS reduction
//--------------------------------------------------------- 
   Int_t sgle1[16], dble1[16];
   Int_t sgle2[16], dble2[16];

   // Calculate simple and double hits
   for (i=0; i<16; i++) {
      dble1[i] = ch1[i] & ch2[i];
      dble2[i] = ch3[i] & ch4[i];
    
      sgle1[i] = (ch1[i]|ch2[i]);
      sgle2[i] = (ch3[i]|ch4[i]);
   }

   //debug
   if(AliDebugLevel()==4||AliDebugLevel()==5) {
      printf("===============================================================\n");
      printf(" Y plane after sgle dble \n"); 
      printf("                            ");
      for (istrip=15; istrip>=0; istrip--) {
         if (istrip>9)  printf("%i",istrip-10*Int_t(istrip/10));
         if (istrip<10) printf("%i",istrip);
      }  
      printf("\n SGLE1                      ");
      for (istrip=15; istrip>=0; istrip--) printf("%i",sgle1[istrip]); 
      printf("\n DBLE1                      ");
      for (istrip=15; istrip>=0; istrip--) printf("%i",dble1[istrip]); 
      printf("\n SGLE2                      ");
      for (istrip=15; istrip>=0; istrip--) printf("%i",sgle2[istrip]); 
      printf("\n DBLE2                      ");
      for (istrip=15; istrip>=0; istrip--) printf("%i",dble2[istrip]); 
      printf(" \n"); 
   }
   //debug

   // DS Reduction 
   Int_t notOr1, notOr2;

   notOr1=!dble1[15] & !dble1[14] & !dble1[13] & !dble1[12] & 
      !dble1[11] & !dble1[10] & !dble1[9]  & !dble1[8]  & 
      !dble1[7]  & !dble1[6]  & !dble1[5]  & !dble1[4]  & 
      !dble1[3]  & !dble1[2]  & !dble1[1]  & !dble1[0];

   notOr2=!dble2[15] & !dble2[14] & !dble2[13] & !dble2[12] & 
      !dble2[11] & !dble2[10] & !dble2[9]  & !dble2[8]  & 
      !dble2[7]  & !dble2[6]  & !dble2[5]  & !dble2[4]  & 
      !dble2[3]  & !dble2[2]  & !dble2[1]  & !dble2[0];

   for (i=0; i<16; i++) {
      sgle1[i] = sgle1[i] & notOr1 & !fCoinc44;
      sgle2[i] = sgle2[i] & notOr2 & !fCoinc44;
   }

//---------------------------------------------------------
// step # 3 : 3/4 coincidence 
//--------------------------------------------------------- 
   Int_t frontImage[16];

   for (i=1; i<15; i++) {
      frontImage[i] = (dble1[i] | sgle1[i]) & 
         (dble2[i+1] | dble2[i] | dble2[i-1]) |
         dble1[i] & (sgle2[i+1] | sgle2[i] | sgle2[i-1]);
   }
   frontImage[0] = (dble1[0] | sgle1[0]) & 
      (dble2[1] | dble2[0]) | dble1[0] & (sgle2[1] | sgle2[0]);

   frontImage[15] = (dble1[15] | sgle1[15]) & 
      (dble2[15] | dble2[14]) | dble1[15] & (sgle2[15] | sgle2[14]);


//debug
   if(AliDebugLevel()==4||AliDebugLevel()==5) {
      printf("===============================================================\n");
      printf(" Y plane frontImage\n");
      printf("                            ");
      for (istrip=15; istrip>=0; istrip--) {
         if (istrip>9)  printf("%i",istrip-10*Int_t(istrip/10));
         if (istrip<10) printf("%i",istrip);
      }
      printf("\n                            ");
      for (istrip=15; istrip>=0; istrip--) printf("%i",frontImage[istrip]); 
      printf("\n");
   }
//debug

//---------------------------------------------------------
// step # 4 : Y position 
//--------------------------------------------------------- 
   Int_t or1, or2, and1, and2, and3;

   or1  = frontImage[7]|frontImage[5]|frontImage[3]|frontImage[1];
   or2  = frontImage[7]|frontImage[6]|frontImage[5]|frontImage[4];
   and1 = !frontImage[3]&!frontImage[2]&!frontImage[1]&!frontImage[0];
   and2 = !frontImage[7]&!frontImage[6]&!frontImage[5]&!frontImage[4] & and1;
   and3 = !frontImage[11]&!frontImage[10]&!frontImage[9]&!frontImage[8]; 
 
   fCoordY[0] = !frontImage[0]&(frontImage[1]|!frontImage[2]) & 
      (frontImage[3]|frontImage[1]|!frontImage[4]&(frontImage[5]|!frontImage[6])) &
      (or1|!frontImage[8]&(frontImage[9]|!frontImage[10])) & 
      (or1|frontImage[11]|frontImage[9]|!frontImage[12]&(frontImage[13]|!frontImage[14]));
 
   fCoordY[1] = !frontImage[0]&!frontImage[1] & 
      !(!frontImage[11]&!frontImage[10]&!frontImage[7]&!frontImage[6] & 
        !frontImage[3]&!frontImage[2]&(frontImage[13]|frontImage[12])) &
      (frontImage[3]|frontImage[2] | !frontImage[5]&!frontImage[4]) & 
      (frontImage[7]|frontImage[6]|frontImage[3]|frontImage[2] | 
       !frontImage[9]&!frontImage[8]);
		
   fCoordY[2] = and1 & (or2 | and3);
		
   fCoordY[3] = and2;
		
   fCoordY[4] = !frontImage[15]&!frontImage[14]&!frontImage[13]&!frontImage[12] &
      and2 & and3 ;
}

//___________________________________________
void AliMUONLocalTriggerBoard::LocalTrigger()
{
//* L0 trigger after LUT
//*
   Int_t deviation=0, iStripY=0;

   for (Int_t i=0; i<4; i++) deviation += static_cast<int>( fMinDev[i] << i );
   for (Int_t i=0; i<4; i++) iStripY   += static_cast<int>( fCoordY[i] << i );

   if (fMinDev[4]==1 && !deviation) fOutput=0;     // No trigger
   else 
   {
      if (fCoordY[4]==1 && iStripY==15) fOutput=0; // No trigger
      else 
         fOutput=1;
   }
  
   if (fOutput) 
   { 
      for (Int_t i=0; i<5; i++) fStripX11 += static_cast<int>( fMinDevStrip[i] << i );

      fDev      = deviation;

      fStripY11 = iStripY;

      Int_t sign = 0;

      if ( !fMinDev[4] &&  deviation ) sign=-1;
      if ( !fMinDev[4] && !deviation ) sign= 0;
      if (  fMinDev[4] == 1 )          sign=+1;    

      fDev *= sign; 

//    calculate deviation in [0;+30]
      fDev += 15;

//    GET LUT OUTPUT FOR icirc/istripX1/deviation/istripY
      fLUT->GetLutOutput(fNumber, fStripX11, fDev, fStripY11, fLutLpt, fLutHpt);
   }  
}

//___________________________________________
Int_t AliMUONLocalTriggerBoard::GetI() const
{
//* old numbering
//*
   const Int_t kMaxfields = 2; char **fields = new char*[kMaxfields];

   char s[100]; strcpy(s, GetName());

   Int_t numlines = 0;

   for (char *token = strtok(s, "B");
        token != NULL;
        token = strtok(NULL, " "))
   {
      fields[numlines] = new char[strlen(token)+1];
      strcpy(fields[numlines++], token);
   }

   TString l(fields[0]);

   char copy = l[0];

   Int_t lL = atoi(&l[4]), cC = atoi(&l[2]), sS = (copy=='R') ? +1 : -1;

   char *b[4] = {"12", "34", "56", "78"};

   Int_t ib = 0;

   for (Int_t i=0; i<4; i++) if (!strcmp(fields[1],b[i])) {ib = i; break;} ib++;

// lL=1 ON TOP
   lL -= 9; lL = abs(lL); lL++;

   Int_t code = 100 * lL + 10 * cC + ib;

   code *= sS;

   Int_t ic = 0;

   for (Int_t i=0; i<234; i++) if (fgkCircuitId[i] == code) {ic = i; break;}

   return ic;
}

//___________________________________________
void AliMUONLocalTriggerBoard::Mask(Int_t index, UShort_t mask)
{
//* set mask
//*
  if ( index >= 0 && index < 2*4 )
  {
    Int_t i = index/4;
    Int_t j = index - i*4;
    fMask[i][j]=mask;
  }
  else
  {
    AliError(Form("Index %d out of bounds (max %d)",index,8));
  }
}

//___________________________________________
void AliMUONLocalTriggerBoard::Scan(Option_t *option) const
{
//* full dump
//*
   TString op = option;

   if (op.Contains("CONF")) Conf();

   if (op.Contains("BITP")) Pattern();

   if (op.Contains("RESPI")) Resp("I");

   if (op.Contains("RESPF")) Resp("F"); 

   if (op.Contains("ALL"))
   {
      Conf();
      Pattern();
      Resp("I");
      Resp("F");
   }
}

//___________________________________________
void AliMUONLocalTriggerBoard::Resp(Option_t *option) const
{
//* board I/O
//*
   TString op = option;

   if (op.Contains("I"))
   {
//    print Local trigger output before the LuT step
      printf("===============================================================\n");
      printf("-------- TRIGGER OUTPUT --------\n");
      printf("minDevStrip = ");
      for  (Int_t i=4; i>=0; i--) printf("%i",fMinDevStrip[i]);
      printf(" minDev = ");
      for  (Int_t i=4; i>=0; i--) printf("%i",fMinDev[i]);
      printf(" coordY = ");
      for  (Int_t i=4; i>=0; i--) printf("%i",fCoordY[i]); 
      printf(" \n");
   }

   if (op.Contains("F"))
   {
      Int_t icirc = GetI();
      Int_t idCircuit = fgkCircuitId[icirc];

      Int_t deviation = 0, iStripY = 0;

      for (Int_t i=0; i<4; i++) iStripY   += static_cast<int>( fCoordY[i] << i );

      for (Int_t i=0; i<4; i++) deviation += Int_t(fMinDev[i]<<i);   

      Float_t pt = 0.; //triggerCircuit->PtCal(fStripX11, fDev, fStripY11);
      printf("-------------------------------------\n");
      printf(" Local Trigger info for circuit Id %i (number %i ) \n", idCircuit, icirc);
      printf(" istripX1 signDev deviation istripY = %i %i %i %i \n", fStripX11, fMinDev[4], deviation, iStripY);
      printf(" pt = %f  (GeV/c) \n", pt);
      printf("-------------------------------------\n");
      printf(" Local Trigger Lut Output = Lpt : ");
      for (Int_t i=1; i>=0; i--) printf("%i", fLutLpt[i]);
      printf(" Hpt : ");
      for (Int_t i=1; i>=0; i--) printf("%i", fLutHpt[i]);
      printf("\n");
      printf("-------------------------------------\n");
   }      
}

//___________________________________________
void AliMUONLocalTriggerBoard::Response()
{
//* algo
//*
   Int_t xX1[16], xX2[16], xXX3[32], xXX4[32];

   TBits x1(16), x2(16), x3(16), x4(16);

   UShort_t xyv = 0;

   xyv = fXY[0][0]; x1.Set(16,&xyv);
   xyv = fXY[0][1]; x2.Set(16,&xyv);
   xyv = fXY[0][2]; x3.Set(16,&xyv);
   xyv = fXY[0][3]; x4.Set(16,&xyv);

   TBits x3u(16), x4u(16), x3d(16), x4d(16);

   xyv = fXYU[0][2]; x3u.Set(16,&xyv);
   xyv = fXYU[0][3]; x4u.Set(16,&xyv);

   xyv = fXYD[0][2]; x3d.Set(16,&xyv);
   xyv = fXYD[0][3]; x4d.Set(16,&xyv);

   for (Int_t i=0;i<16;i++)
   {
      xX1[i] = x1[i];
      xX2[i] = x2[i];
      
      xXX3[i+8] = x3[i];
      xXX4[i+8] = x4[i];  
   }

   for (Int_t i=0;i<8;i++)
   {
      xXX3[i] = x3d[i+8];
      xXX4[i] = x4d[i+8];

      xXX3[i+24] = x3u[i];
      xXX4[i+24] = x4u[i];
   }
   
//   Int_t coinc44 = 0;
   
   TrigX(xX1, xX2, xXX3, xXX4);   

   Int_t yY1[16], yY2[16], yY3[16], yY4[16];
   
   Int_t yY3U[16], yY3D[16], yY4U[16], yY4D[16];

   TBits y1(16), y2(16), y3(16), y4(16);

   TBits y3u(16), y3d(16), y4u(16), y4d(16);

   xyv = fXY[1][0]; y1.Set(16,&xyv);
   xyv = fXY[1][1]; y2.Set(16,&xyv);
   xyv = fXY[1][2]; y3.Set(16,&xyv);
   xyv = fXY[1][3]; y4.Set(16,&xyv);

   xyv = fXYU[1][2]; y3u.Set(16,&xyv);
   xyv = fXYD[1][2]; y3d.Set(16,&xyv);
   xyv = fXYU[1][3]; y4u.Set(16,&xyv);
   xyv = fXYD[1][3]; y4d.Set(16,&xyv);

   for (Int_t i=0;i<16;i++)
   {
      yY1[i] = y1[i];
      yY2[i] = y2[i];
      yY3[i] = y3[i];
      yY4[i] = y4[i];
      
      yY3U[i] = y3u[i];
      yY3D[i] = y3d[i];
      
      yY4U[i] = y4u[i];
      yY4D[i] = y4d[i];
   }

   TrigY(yY1, yY2, yY3, yY4, yY3U, yY3D, yY4U, yY4D);
   
// ASIGN fLutLpt, fLutHpt
   LocalTrigger(); 

   fResponse = fLutLpt[0]                      + 
       static_cast<int>(fLutLpt[1]<<1) + 
       static_cast<int>(fLutHpt[0]<<2) + 
       static_cast<int>(fLutHpt[1]<<3);
}

ClassImp(AliMUONLocalTriggerBoard)

