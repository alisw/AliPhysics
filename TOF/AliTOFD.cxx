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
Revision 1.3  2000/10/02 21:28:17  fca
Removal of useless dependecies via forward declarations

Revision 1.2  2000/05/18 14:33:01  vicinanz
Modified to be full HP compliant

Revision 1.1  2000/05/10 16:52:18  vicinanz
New TOF version with holes for PHOS/RICH

*/

#include "AliTOF.h"
#include "AliTOFD.h"
#include "TObject.h"

//******************************************************************************

ClassImp(AliTOFRawDigit)

//______________________________________________________________________________
AliTOFRawDigit::AliTOFRawDigit()
//
// Constructor of AliTOFRawDigit class
//
{
  fTreeD     = 0;
  fRawDigits = 0;
}

//******************************************************************************

ClassImp(AliTOFRoc)

//______________________________________________________________________________
AliTOFRoc::AliTOFRoc()
//
// Constructor of AliTOFRoc class
// The class represents a ROC in the TARODA system
//
{
  fItems = 0;
  fSize  = 0;
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
   return fItems;   
}

//______________________________________________________________________________
void AliTOFRoc::SetHeader()
//
// Calculate the header line of the ROC in the raw data file
//
{
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
   iCharge = (UInt_t)(RealCharge/50.);
   if(iCharge>=TMath::Power(2,16)) iCharge = (UInt_t)TMath::Power(2,16)-1;
   fChrgRow[Item] = iCharge+fNRoc+Fec+Tdc;
}

//______________________________________________________________________________
void AliTOFRoc::SetTime(UInt_t Item, UInt_t tir)
//
// Writes the raw data line relative to the TDC
//
{
   fChrgRow[Item]=tir;
}

//______________________________________________________________________________
void AliTOFRoc::SetCharge(UInt_t Item, UInt_t chr)
//
// Writes the raw data line relative to the ADC
//
{
   fChrgRow[Item]=chr;
}

//______________________________________________________________________________
Float_t AliTOFRoc::GetCharge(Int_t Item)
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
   Error = fTimeRow[Item]>>24;
   return time; 
}

//______________________________________________________________________________
Int_t AliTOFRoc::GetTotPad(Int_t Item)
{
//
// Reads the cohordinates of the pad starting
// from the line of the raw data
//

   UInt_t nRoc = (fChrgRow[Item]&0xf0000000)>>28;
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
UInt_t AliTOFRoc::BitCount(UInt_t x)
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


//******************************************************************************

ClassImp(AliTOFRawSector)

//______________________________________________________________________________
AliTOFRawSector::AliTOFRawSector()
{
//
// Constructor of AliTOFRawSector class
// Each sector is in effect a 
// TClonesArray of 14 AliTOFRoc Objects
//
   fRocData = new TClonesArray("AliTOFRoc",14);   
}

//______________________________________________________________________________
AliTOFRawSector::~AliTOFRawSector()
{
// destructor of the AliTOFRawSector object
// Here we delete the 14 AliTOFRoc istances
//
   if (fRocData) {
       fRocData->Delete() ;
       delete fRocData;
       fRocData = 0;
   }
}

//______________________________________________________________________________
void AliTOFRawSector::WriteSector()
{
//
// Starting from the raw data objects writes a binary file
// similar to real raw data.
//

    FILE *rawfile;
    rawfile = fopen("rawdata.dat","w");
    
//    fprintf(rawfile,Header);
    
    Int_t nRoc;

// loop on all AliTOFRoc Headers to set them    
    for(nRoc=1; nRoc<=14; nRoc++){
       AliTOFRoc* currentRoc = (AliTOFRoc*)fRocData->UncheckedAt(nRoc);
       currentRoc->SetHeader();
       //       UInt_t RocHeader = currentRoc->fHeader;
//      fprintf(rawfile,RocHeader);
    }
    
    for(nRoc=1; nRoc<=14; nRoc++){
       AliTOFRoc* currentRoc = (AliTOFRoc*)fRocData->UncheckedAt(nRoc);
       Int_t rocItems = currentRoc->GetItems();

       for(Int_t nItem=1; nItem<=rocItems;nItem++){
	 //          UInt_t TimeRow = currentRoc->GetTimeRow(nItem);
//          fprintf(rawfile,TimeRow);
	  //          UInt_t ChrgRow = currentRoc->GetTimeRow(nItem);
//          fprintf(rawfile,ChrgRow);
       }
    }
    
    //    UInt_t EndOfSector = GlobalCheckSum;
//    fprintf(rawfile,EndOfSector);
}

//______________________________________________________________________________
void AliTOFRawSector::ReadSector()
{
//
// Starting from raw data initialize and write the 
// Raw Data objects 
//(i.e. a TClonesArray of 18 AliTOFRawSector)
//

    FILE *rawfile;
    rawfile = fopen("rawdata.dat","r");
    
//    fscanf(rawfile,Header);
//    fscanf(rawfile,Header);
    Int_t nRoc;
    
    for(nRoc=1; nRoc<=14; nRoc++){
       AliTOFRoc* currentRoc = (AliTOFRoc*)fRocData->UncheckedAt(nRoc);
       UInt_t rocHeader;
       fscanf(rawfile,"%u",&rocHeader);
       currentRoc->SetHeader(rocHeader);
    }
    
//      UInt_t SCMWord;
//      fscanf(rawfile,"%u",SCMWord);
    
    for(nRoc=1; nRoc<=14; nRoc++){
       AliTOFRoc* currentRoc = (AliTOFRoc*)fRocData->UncheckedAt(nRoc);
       //       Int_t Size = currentRoc->SetSize();
       Int_t nItems = currentRoc->GetItems();
       for(Int_t nrow=0; nrow<=nItems; nrow++){
          UInt_t charRow,timeRow;
	  fscanf(rawfile,"%u",&charRow);
	  currentRoc->SetTime(nrow, charRow);
          fscanf(rawfile,"%u",&timeRow);
	  currentRoc->SetTime(nrow, timeRow);
       }
         Int_t FinalWord;
         fscanf(rawfile,"%d",&FinalWord);              
    }
//    fscanf(rawfile,GlobalCheckSum);
}

