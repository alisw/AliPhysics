////////////////////////////////////////////////
//  Digitization class for set: TOF           //
//  AliTOFRawSector  class                    //
//  Member variables                          //
//                                            //
//  Member functions                          //
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

#include "TClonesArray.h"

#include "AliTOFRawSector.h"
#include "AliTOFRoc.h"

ClassImp(AliTOFRawSector)

//______________________________________________________________________________
AliTOFRawSector::AliTOFRawSector()
{
//
// Constructor of AliTOFRawSector class
// Each sector is in effect a 
// TClonesArray of 14 AliTOFRoc Objects
//
   fRocData = 0;   
}

//______________________________________________________________________________
AliTOFRawSector::AliTOFRawSector(const AliTOFRawSector& tofrawsector)
:TObject()
// : fHeader(tofrawsector.fHeader), fGlobalCheckSum(tofrawsector.fGlobalCheckSum) 
{
// copy ctor for AliTOFRawSector class
// (required also by RC10 Coding Convention)
//
// we make here a new istance of the TClonesArray containing Roc Data
  fRocData = new TClonesArray("AliTOFRoc",14);

// we make here a copy of Roc Data

  for(Int_t nroc=1; nroc<=14; nroc++){
// get the pointers to the current roc of new and to be copied TClonesArray
    AliTOFRoc* currentrocnew = (AliTOFRoc*)fRocData->UncheckedAt(nroc);
    AliTOFRoc* currentrocold = (AliTOFRoc*)tofrawsector.fRocData->UncheckedAt(nroc);

// we create here 2 references: one for the new current roc and another for the old one
    AliTOFRoc& newroc = *currentrocnew;
    AliTOFRoc& oldroc = *currentrocold;

    newroc = oldroc; // 'operator =' called for AliTOFRoc

  } // end loop on Roc Data
}

//______________________________________________________________________________
AliTOFRawSector& AliTOFRawSector::operator=(const  AliTOFRawSector& tofrawsector)
{
// Assignment operator for AliTOFRawSector 
// (required also by RC10 Coding Conventions)
//
    if (this !=&tofrawsector) { // do nothing if assigned to self
       fHeader=tofrawsector.fHeader;
       fGlobalCheckSum=tofrawsector.fGlobalCheckSum;
// loop on ROC data
         for(Int_t nroc=1; nroc<=14; nroc++){
// get the pointers to the current roc of new and to be copied TClonesArray
            AliTOFRoc* currentrocnew = (AliTOFRoc*)fRocData->UncheckedAt(nroc);
            AliTOFRoc* currentrocold = (AliTOFRoc*)tofrawsector.fRocData->UncheckedAt(nroc);
// we create here 2 references: one for the new current roc and another for the old one
            AliTOFRoc& newroc = *currentrocnew;   
            AliTOFRoc& oldroc = *currentrocold;
            newroc = oldroc; // 'operator =' called for AliTOFRoc
         } // end loop on Roc Data
    } // close if
    return *this;
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
         Int_t finalWord;
         fscanf(rawfile,"%d",&finalWord);              
    }
//    fscanf(rawfile,GlobalCheckSum);
}



