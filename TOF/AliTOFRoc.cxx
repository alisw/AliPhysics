////////////////////////////////////////////////
//  Digitization class for set: TOF           //
//  AliTOFRoc  class                          //
//  Member variables                          //
//   fItems :  number of items		      //
//   fSize  :  size                           //
//   fNRoc  :  Roc number                     //
//   fHeader:  Roc header number              //
//   fChrgRow[1024]; // adc values            //
//   fTimeRow[1024]; // tdc values            //
//                                            //
//  Member functions implemented here         //
//                                            //
//*-- Authors: Pierella, Seganti, Vicinanza   //
//    (Bologna and Salerno University)        //
////////////////////////////////////////////////


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

#include <Riostream.h>
#include <assert.h>

#include "AliTOFRoc.h"


ClassImp(AliTOFRoc)

//______________________________________________________________________________
AliTOFRoc::AliTOFRoc()
{
//
// Constructor of AliTOFRoc class
// The class represents a ROC in the TARODA system
// here we make the initialization of the member variables
  fItems = 0;
  fSize  = 0;
  fNRoc  = 0;
  fHeader= 0;
// initialization of fChrgRow[1024] and fTimeRow[1024]
  for(Int_t i=0; i < 1024; i++){
     fChrgRow[i] = 0;
     fTimeRow[i] = 0;
  } // end loop
}

//______________________________________________________________________________
AliTOFRoc::AliTOFRoc(const AliTOFRoc& tofroc)
:TObject()
// : fItems(tofroc.fItems), fSize(tofroc.fSize), fNRoc(tofroc.fNRoc), fHeader(tofroc.fHeader) 
{ 
//
// copy ctor for AliTOFRoc class 
//
   assert(tofroc.fItems >= 0); // check for number of items
   assert(tofroc.fSize  >= 0); // check for roc size
// making a copy of adc and tdc vectors
   for(Int_t i=0; i < 1024; i++){
      fChrgRow[i] = tofroc.fChrgRow[i]; // coping adc values
      fTimeRow[i] = tofroc.fTimeRow[i]; // coping tdc values
   } // end loop
}

//______________________________________________________________________________
AliTOFRoc& AliTOFRoc::operator=(const AliTOFRoc& tofroc)
{
//
// Assignment operator for AliTOFRoc
// (used by copy ctor of AliTOFRawSector)
//
   if (this !=&tofroc) { // do nothing if assigned to self
     // setting head member data
     SetHeadVar(tofroc.fItems,tofroc.fSize,tofroc.fNRoc,tofroc.fHeader);
     // loop on adc and tdc values
     for(Int_t i=0; i < 1024; i++){
      fChrgRow[i] = tofroc.fChrgRow[i]; // coping adc values
      fTimeRow[i] = tofroc.fTimeRow[i]; // coping tdc values
     } // end loop
   }
   return *this;
}

//______________________________________________________________________________
AliTOFRoc::~AliTOFRoc(){}

//______________________________________________________________________________
Int_t AliTOFRoc::AddItem(Int_t Fec, Int_t Tdc, Int_t Error,
                         Float_t Charge, Float_t Time)
{
//
// Adds an item (i.e. the charge, the TOF and the 
// cohordinates of a hit pad) to the ROC class.
//
   fItems++;
   SetCharge(fItems,Fec,Tdc,Charge);
   SetTime  (fItems,Error,Time);
   return fItems; // return the number of current items   
}

//______________________________________________________________________________
void AliTOFRoc::SetHeadVar(Int_t items, Int_t size, Int_t nroc, UInt_t header)
{
//
// set header member variables for AliTOFRoc 
//
  fItems = items;
  fSize  = size ;
  fNRoc  = nroc ;
  fHeader= header ;
}

//______________________________________________________________________________
void AliTOFRoc::SetHeader()
{
//
// Calculate the header line of the ROC in the raw data file
//

   fHeader  = fNRoc<<28;
   fHeader += fSize;
}


//______________________________________________________________________________
void AliTOFRoc::SetTime(UInt_t Item, UInt_t Error, Float_t RealTime)
{
//
// Calculate the raw data line relative to the TDC
// output of a pad in the current ROC.
//

   UInt_t itime;
   itime = (UInt_t)(RealTime/50.);
   if (itime >= TMath::Power(2,24)) itime = 2^24-1;
   Error <<= 24;
   fTimeRow[Item]= Error+itime;
}

//______________________________________________________________________________
void AliTOFRoc::SetCharge(UInt_t Item, UInt_t Fec, UInt_t Tdc, Float_t RealCharge)
{
//
// Calculate the raw data line relative to the ADC 
// output of a pad in the current ROC.
//

   UInt_t iCharge;
   if (fNRoc>=TMath::Power(2,4)) fNRoc = 0;
   fNRoc <<= 28;
   if (Fec >=TMath::Power(2,6))  Fec = 0;
   Fec  <<= 22;
   if (Tdc >=TMath::Power(2,6))  Tdc = 0;
   Tdc  <<= 16;
   iCharge = (UInt_t)(RealCharge/50.); // 50 ps (TDC bin value)
   if(iCharge>=TMath::Power(2,16)) iCharge = (UInt_t)TMath::Power(2,16)-1;
   fChrgRow[Item] = iCharge+fNRoc+Fec+Tdc;
}

//______________________________________________________________________________
void AliTOFRoc::SetTime(UInt_t Item, UInt_t tir)
{
//
// Writes the raw data line relative to the TDC
//

   fChrgRow[Item]=tir;
}

//______________________________________________________________________________
void AliTOFRoc::SetCharge(UInt_t Item, UInt_t chr)
{
//
// Writes the raw data line relative to the ADC
//

   fChrgRow[Item]=chr;
}

//______________________________________________________________________________
Float_t AliTOFRoc::GetCharge(Int_t Item) const
{
//
// Reads the effective value of the charge starting
// from the line of the raw data
//

   UInt_t  icharge  = fChrgRow[Item]&0x0000ffff;
   Float_t charge = (Float_t)icharge*50.;
   return charge;
}

//______________________________________________________________________________
Float_t AliTOFRoc::GetTime(Int_t Item, UInt_t& Error) 
{
//
// Reads the effective value of the time of flight starting
// from the line of the raw data
//

   UInt_t  itime  = fTimeRow[Item]&0x00ffffff;
   Float_t time = (Float_t)itime*50.;
   Error = fTimeRow[Item]>>24; // the same as Error = fTimeRow[Item] / 24;
   return time; 
}

//______________________________________________________________________________
Int_t AliTOFRoc::GetTotPad(Int_t Item) const
{
//
// Reads the cohordinates of the pad starting
// from the line of the raw data
//

   UInt_t nRoc = (fChrgRow[Item]&0xf0000000)>>28; // >> the same as / (division by)
   UInt_t nFec = (fChrgRow[Item]&0x0fc00000)>>22;
   UInt_t nTdc = (fChrgRow[Item]&0x003f0000)>>16;
   UInt_t pad = nRoc*32*32+nFec*32+nTdc;
   return pad; 
}

//______________________________________________________________________________
UInt_t AliTOFRoc::GetCheckSum()
{
//
// Calculate the checksum word of the current ROC
// 

   UInt_t checkSum=0;
   for(Int_t i=0; i<fItems; i++){
      checkSum += BitCount(GetChrgRow(i));
      checkSum += BitCount(GetTimeRow(i));
   }
   return checkSum;
}

//______________________________________________________________________________
UInt_t AliTOFRoc::BitCount(UInt_t x) const
{
//
// Count the "1" bit in the current word
//

   UInt_t count=0;
   for (count=0; x!=0; x>>=1){
      if(x&0x00000001) count++;
   }
   return count;
}

//______________________________________________________________________________
UInt_t AliTOFRoc::SetSize()
{
//
// Reads the size of data from current ROC starting
// from the header line of the raw data
//

   fSize = fHeader&0x0000ffff;
   fItems = (fSize-4)/4;
   return fSize;
}



