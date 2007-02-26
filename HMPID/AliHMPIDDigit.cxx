//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************

#include "AliHMPIDDigit.h" //class header
#include <TClonesArray.h>  //Hit2Sdi() 
#include <TBox.h>       //Draw() 
ClassImp(AliHMPIDDigit)

const Float_t AliHMPIDDigit::fMinPcX[]={0,SizePcX() + SizeDead(),0,SizePcX() + SizeDead(),       0,SizePcX() + SizeDead()};
const Float_t AliHMPIDDigit::fMinPcY[]={0,                     0,SizePcY()+SizeDead(),SizePcY() + SizeDead(),2*(SizePcY()+SizeDead()),2*(SizePcY()+SizeDead())};

const Float_t AliHMPIDDigit::fMaxPcX[]={SizePcX(),SizeAllX(),SizePcX(),SizeAllX(),SizePcX(),SizeAllX()};
const Float_t AliHMPIDDigit::fMaxPcY[]={SizePcY(),SizePcY(),SizeAllY() - SizePcY(),SizeAllY() - SizePcY(),SizeAllY(),SizeAllY()};


/*
Preface: all geometrical information (like left-right sides) is reported as seen from electronic side.


The DDL file starts with common header which size and structure is standartized and mandatory for all detectors. 
The header contains among other words, so called Equipment ID word. This unique value for each D-RORC is calculated as detector ID << 8 + DDL index. 
For HMPID the detector ID is 6 (reffered in the code as kRichRawId) while DDL indexes are from 0 to 13.

Common header might be followed by the private one although  HMPID has no any private header, just uses the common one.

Single HMPID D-RORC (with 2 channels) serves a single chamber so that channel 0 serves left half (PCs 0-2-4) 
                                                                              1 serves right half(PCs 1-3-5) 

So the LDC -chamber-ddl map is:
DDL index  0 -> ch 0 left -> DDL ID 0x600          DDL index  1 -> ch 1 right -> DDL ID 0x601 
DDL index  2 -> ch 1 left -> DDL ID 0x602          DDL index  3 -> ch 2 right -> DDL ID 0x603 
DDL index  4 -> ch 2 left -> DDL ID 0x604          DDL index  5 -> ch 3 right -> DDL ID 0x605 
DDL index  6 -> ch 3 left -> DDL ID 0x606          DDL index  7 -> ch 4 right -> DDL ID 0x607 
DDL index  8 -> ch 4 left -> DDL ID 0x608          DDL index  9 -> ch 5 right -> DDL ID 0x609 
DDL index 10 -> ch 5 left -> DDL ID 0x60a          DDL index 11 -> ch 6 right -> DDL ID 0x60b 
DDL index 12 -> ch 6 left -> DDL ID 0x60c          DDL index 13 -> ch 7 right -> DDL ID 0x60d 

HMPID FEE as seen by single D-RORC is composed from a number of DILOGIC chips organized in vertical stack of rows. 
Each DILOGIC chip serves 48 channels for the 8x6 pads Channels counted from 0 to 47.

The mapping inside DILOGIC chip has the following structure (see from electronics side):
pady

5  04 10 16 22 28 34 40 46                   due to repetition in column structure we may introduce per column map:   
4  02 08 14 20 26 32 38 44                   pady= 0 1 2 3 4 5          
3  00 06 12 18 24 30 36 42                   addr= 5 3 1 0 2 4
2  01 07 13 19 25 31 37 43                   or vice versa 
1  03 09 15 21 27 33 39 45                   addr= 0 1 2 3 4 5
0  05 11 17 23 29 35 41 47                   pady= 3 2 4 1 5 0  
    
    0  1  2  3  4  5  6  7  padx

10 DILOGIC chips composes so called "row" in horizontal direction (reffered in the code as kNdil), so the row is 80x6 pads structure. 
DILOGIC chips in the row are counted from right to left as seen from electronics side, from 1 to 10.
24 rows are piled up forming the whole FEE served by single D-RORC, so one DDL sees 80x144 pads separated in 3 photocathodes.
Rows are counted from 1 to 24 from top    to bottom for right  half of the chamber (PCs 1-3-5) as seen from electronics side, meaning even LDC number
                          and from bottom to top    for left   half of the chamber (PCs 0-2-4) as seen from electronics side, meaning odd LDC number.

HMPID raw word is 32 bits with the structure:   
   00000             rrrrr                      dddd                               aaaaaa                          qqqqqqqqqqqq 
 5 bits zero  5 bits row number (1..24)  4 bits DILOGIC chip number (1..10) 6 bits DILOGIC address (0..47)  12 bits QDC value (0..4095)
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigit::Draw(Option_t*)
{
//  TMarker *pMark=new TMarker(LorsX(),LorsY(),25); pMark->SetMarkerColor(kGreen);pMark->Draw();
  TBox *pad = new TBox(LorsX()-0.5*SizePadX(),LorsY()-0.5*SizePadY(),LorsX()+0.5*SizePadX(),LorsY()+0.5*SizePadY());
  pad->SetFillStyle(0);pad->SetLineColor(kGreen);
  pad->Draw();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigit::Hit2Sdi(AliHMPIDHit *pHit,TClonesArray *pSdiLst)
{
// Creates a list of sdigits out of provided hit
// Arguments: pHit- hit
//   Returns: none
  
  Int_t iSdiCnt=pSdiLst->GetEntries(); //list of sdigits contains sdigits from previous ivocations of Hit2Sdi, do not override them
  AliHMPIDDigit dig;
  for(Int_t i=0;i<9;i++){                                      //affected pads loop
    dig.Set(pHit,i); //c,q,tid,x,y   create tmp sdigit for pad i around hit position
    if(dig.PadPcX()==-1) continue;
    new((*pSdiLst)[iSdiCnt++]) AliHMPIDDigit(dig);
  }
}//Hit2Sdi()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigit::Print(Option_t*)const
{
// Print current digit  
// Arguments: option string not used
//   Returns: none    
  UInt_t w32; Raw(w32);
  Printf("DIG:(%7.3f,%7.3f) Q=%8.3f (ch=%1i,pc=%1i,x=%2i,y=%2i)  TID=(%5i,%5i,%5i) ddl=%i raw=0x%x (r=%2i,d=%2i,a=%2i) %s",
              LorsX(),LorsY(),Q(), A2C(fPad),A2P(fPad),A2X(fPad),A2Y(fPad),   
              fTracks[0],fTracks[1],fTracks[2],DdlIdx(),w32,Row(),Dilogic(),Addr(), 
              (IsOverTh(Q()))?"":"!!!");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigit::PrintSize()
{
// Print all segmentaion related sizes
// Arguments: none
//   Returns: none    
  Printf("-->pad    =(%6.2f,%6.2f) cm dead zone %.2f cm\n"
         "-->PC     =(%6.2f,%6.2f) cm (%3i,%3i) pads\n"
         "-->all PCs=(%6.2f,%6.2f) cm (%3i,%3i) pads",
               SizePadX(),SizePadY(),SizeDead(),
               SizePcX() ,SizePcY() ,kPadPcX ,kPadPcY,
               SizeAllX(),SizeAllY(),kPadAllX,kPadAllY);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigit::Test()
{
  AliHMPIDDigit d1,d2; Int_t ddl;UInt_t w32;
  for(Int_t i=0;i<10000;i++){
    Int_t ch=Int_t(gRandom->Rndm()*7);
    Int_t pc=Int_t(gRandom->Rndm()*6);
    Int_t px=Int_t(gRandom->Rndm()*80);
    Int_t py=Int_t(gRandom->Rndm()*48);
    d1.Manual2(ch,pc,px,py);                
    ddl=d1.Raw(w32);    d2.Raw(ddl,w32);
    if(d1.Compare(&d2)) Printf("Problem!!!");
  }
  Printf("OK");
}//Test()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

