/***************************************************************************
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
#include <TClonesArray.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>

#include "AliRawDataHeader.h"
#include "AliPMDdigit.h"
#include "AliPMDRawStream.h"
#include "AliPMDDDLRawData.h"


ClassImp(AliPMDDDLRawData)

AliPMDDDLRawData::AliPMDDDLRawData():
  fDigits(new TClonesArray("AliPMDdigit", 1000))
{
  // Default Constructor
  //

}
//____________________________________________________________________________

AliPMDDDLRawData::~AliPMDDDLRawData()
{
  // Default Destructor
  //

}

//____________________________________________________________________________
void AliPMDDDLRawData::WritePMDRawData(TTree *treeD)
{
  // write digits into raw data format

  ofstream outfile;

  TBranch *branch = treeD->GetBranch("PMDDigit");
  if (!branch) return;
  branch->SetAddress(&fDigits);  
  
  //  Int_t   nmodules = (Int_t) treeD->GetEntries();
  //  cout << " nmodules = " << nmodules << endl;

  const Int_t kSize         = 4608;
  const Int_t kDDL          = 6;
  Int_t modulePerDDL        = 0;


  AliRawDataHeader header;
  UInt_t sizeRawData = 0;
  Int_t  totword     = 0;

  UInt_t buffer[kSize];

  Char_t filename[80];

  Int_t mmodule = 0;
  for(Int_t iddl = 0; iddl < kDDL; iddl++)
    {
      sprintf(filename,"PMD_%d.ddl",AliPMDRawStream::kDDLOffset+iddl);
#ifndef __DECCXX
      outfile.open(filename,ios::binary);
#else
      outfile.open(filename);
#endif
      
      if (iddl < 4)
	{
	  modulePerDDL = 6;
	  mmodule = 6*iddl;
	}
      else if (iddl == 4)
	{
	  modulePerDDL = 12;
	  mmodule = 24;
	}
      else if (iddl == 5)
	{
	  modulePerDDL = 12;
	  mmodule = 30;
	}


      // Write the Dummy Data Header into the file
      Int_t bHPosition = outfile.tellp();
      outfile.write((char*)(&header),sizeof(header));


      for(Int_t ium = 0; ium < modulePerDDL; ium++)
	{
	  
	  if (iddl == 4 && ium == 6) mmodule = 36;
	  if (iddl == 5 && ium == 6) mmodule = 42;

	  for (Int_t i = 0; i < kSize; i++)
	    {
	      buffer[i] = 0;
	    }
	  // Extract energy deposition per cell and pack it
	  // in a 32-bit word and returns all the total words
	  // per one unit-module
	  
	  GetUMDigitsData(treeD, mmodule, ium, iddl, totword, buffer);
	  
	  outfile.write((char*)buffer,totword*sizeof(UInt_t));

	  mmodule++;

	}

      // Write real data header
      // take the pointer to the beginning of the data header
      // write the total number of words per ddl and bring the
      // pointer to the current file position and close it
      UInt_t cFPosition = outfile.tellp();
      sizeRawData = cFPosition - bHPosition - sizeof(header);
      header.fSize = cFPosition - bHPosition;
      header.SetAttribute(0);  // valid data
      outfile.seekp(bHPosition);
      outfile.write((char*)(&header),sizeof(header));
      outfile.seekp(cFPosition);

      outfile.close();
    } // DDL Loop over


}
//____________________________________________________________________________
void AliPMDDDLRawData::GetUMDigitsData(TTree *treeD, Int_t imodule, Int_t ium,
				       Int_t ddlno, Int_t & totword,
				       UInt_t *buffer)
{

  UInt_t dataword, baseword;

  UInt_t mcmno, chno;
  UInt_t adc;
  Int_t  irownew = 0;
  Int_t  icolnew = 0;
  Int_t  irownew1 = 0;
  Int_t  icolnew1 = 0;
  Int_t  det, smn, irow, icol;


  treeD->GetEntry(imodule); 
  Int_t nentries = fDigits->GetLast();
  totword = nentries+1;

  for (Int_t ient = 0; ient < totword; ient++)
    {
      fPMDdigit = (AliPMDdigit*)fDigits->UncheckedAt(ient);
      
      det    = fPMDdigit->GetDetector();
      smn    = fPMDdigit->GetSMNumber();
      irow   = fPMDdigit->GetRow();
      icol   = fPMDdigit->GetColumn();
      adc    = (UInt_t) fPMDdigit->GetADC();

      //      cout << " imodule = " << imodule << " smn = " << smn << endl;

      if(smn < 12)
	{
	  irownew1 = icol;
	  icolnew1 = irow;
	}
      else if( smn >= 12 && smn < 24)
	{
	  irownew1 = irow;
	  icolnew1 = icol;
	}

      if(smn < 6)
	{
	  irownew = 95 - irownew1;
	  icolnew = icolnew1;
	}
      else if(smn >= 6 && smn < 12)
	{
	  irownew = irownew1;
	  icolnew = 47 - icolnew1;
	}
      else if(smn >= 12 && smn < 18)
	{
	  irownew = 47 - irownew1;
	  icolnew = icolnew1;
	}
      else if(smn >= 18 && smn < 24)
	{
	  irownew = irownew1;
	  icolnew = 95 - icolnew1;
	}
      


      GetMCMCh(ddlno, ium, irownew, icolnew, mcmno, chno);

      baseword = 0;
      dataword = adc;
      PackWord(0, 11,dataword,baseword);
      dataword = chno;
      PackWord(12,17,dataword,baseword);
      dataword = mcmno;
      PackWord(18,28,dataword,baseword);
      dataword = 0;
      PackWord(29,29,dataword,baseword);
      dataword = 0;
      PackWord(30,30,dataword,baseword);
      dataword = 1;
      PackWord(31,31,dataword,baseword);
      
      buffer[ient] = baseword;
      
    }


}
//____________________________________________________________________________
void AliPMDDDLRawData::GetMCMCh(Int_t ddlno, Int_t um,
				Int_t row, Int_t col,
				UInt_t &mcmno, UInt_t &chno)
{
  // This part will be modified once the final track layout on the PCB is
  // designed. This will only change the coordinate of the individual cell

  static const UInt_t kCh[16][4] = { {12, 13, 18, 19},
				     {11, 15, 17, 21},
				     {10, 14, 16, 22},
				     {8, 9, 20, 23},
				     {7, 4, 25, 24},
				     {6, 0, 30, 26},
				     {5, 1, 31, 27},
				     {3, 2, 29, 28},
				     {44, 45, 50, 51},
				     {43, 47, 49, 53},
				     {42, 46, 48, 54},
				     {40, 41, 52, 55},
				     {39, 36, 57, 56},
				     {38, 32, 62, 58},
				     {37, 33, 63, 59},
				     {35, 34, 61, 60} };
  
  if (ddlno == 0 || ddlno == 1)
    {
      // PRE plane, SU Mod = 1, 2
      Int_t irownew = row%16;
      Int_t icolnew = col%4;
      Int_t irowdiv = row/16;
      Int_t icoldiv = col/4;

      mcmno = 72*um + 12*irowdiv + icoldiv;
      chno  = kCh[irownew][icolnew];
    }
  else if (ddlno == 2 || ddlno == 3)
    {
      // PRE plane,  SU Mod = 3, 4
      Int_t irownew = row%16;
      Int_t icolnew = col%4;
      Int_t irowdiv = row/16;
      Int_t icoldiv = col/4;
      
      mcmno = 72*um + 24*irowdiv + icoldiv;
      chno  = kCh[irownew][icolnew];
    }
  else if (ddlno == 4 || ddlno == 5)
    {
      // CPV plane,  SU Mod = 1, 3 : ddl = 4
      // CPV plane,  SU Mod = 2, 4 : ddl = 5

      Int_t irownew = row%16;
      Int_t icolnew = col%4;
      Int_t irowdiv = row/16;
      Int_t icoldiv = col/4;

      if(um < 6)
	{
	  mcmno = 72*um + 12*irowdiv + icoldiv;
	}
      else if(um >= 6)
	{
	  mcmno = 72*um + 24*irowdiv + icoldiv;
	}
      chno  = kCh[irownew][icolnew];
    }

}
//____________________________________________________________________________

void AliPMDDDLRawData::PackWord(UInt_t startbit, UInt_t stopbit, 
				UInt_t dataword, UInt_t &packedword)
{
  UInt_t bitLength  = stopbit - startbit + 1;
  UInt_t bitContent = (UInt_t) (TMath::Power(2,bitLength) - 1);
  if(bitContent < dataword)
    {
      cout << " *** ERROR *** bitContent is less than the dataword" << endl;
      return;
    }
  UInt_t packedBits = 0;
  if (packedword != 0)
    packedBits = (UInt_t) (TMath::Log(packedword)/TMath::Log(2));

  UInt_t counter;
  if (packedBits <= stopbit)
    {
      counter   = 31 - stopbit;
    }
  else
    {
      counter   = 31 - packedBits;
    }
  UInt_t dummyword = 0xFFFFFFFF;
  dummyword >>= counter;
  UInt_t lword = dataword << startbit;
  UInt_t nword = lword | packedword;
  packedword = dummyword & nword;


}
//____________________________________________________________________________
void AliPMDDDLRawData::UnpackWord(UInt_t startbit, UInt_t stopbit, 
				UInt_t &dataword, UInt_t packedword)
{
  UInt_t bitLength  = stopbit - startbit + 1;
  UInt_t bitContent = (UInt_t) (TMath::Power(2,bitLength) - 1);
  bitContent <<= startbit;
  dataword = packedword & bitContent;
  dataword >>= startbit;

}
//____________________________________________________________________________

