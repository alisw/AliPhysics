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

#include "AliPMDdigit.h"
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
void AliPMDDDLRawData::WritePMDRawData(TTree *treeD, Int_t evtno)
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


  UInt_t sizeRawData = 0;
  UInt_t magicWord   = 0x123456;
  UInt_t detectorID  = 0;
  UInt_t ddlID       = 0;
  Int_t  flag        = 0;
  Int_t  version     = 1;
  UInt_t mHSize      = 3*sizeof(UInt_t);
  Int_t  totword     = 0;

  UInt_t miniHeader[3];
  UInt_t buffer[kSize];

  Char_t filename[80];

  Int_t mmodule = 0;
  for(Int_t iddl = 0; iddl < kDDL; iddl++)
    {
      sprintf(filename,"Ev%dPMDddl%d.dat",evtno,iddl);
#ifndef __DECCXX
      outfile.open(filename,ios::binary);
#else
      outfile.open(filename);
#endif
      
      if (iddl < 4)
	{
	  modulePerDDL = 6;
	  mmodule = 6*iddl;
	  detectorID = 0;
	}
      else if (iddl == 4)
	{
	  modulePerDDL = 12;
	  mmodule = 24;
	  detectorID = 1;
	}
      else if (iddl == 5)
	{
	  modulePerDDL = 12;
	  mmodule = 36;
	  detectorID = 1;
	}

      miniHeader[0] = sizeRawData;
      PackWord(0,23,magicWord,miniHeader[1]);
      PackWord(24,31,detectorID,miniHeader[1]);
      ddlID = iddl;
      PackWord(0,15,ddlID,miniHeader[2]);
      PackWord(16,23,flag,miniHeader[2]);
      PackWord(24,31,version,miniHeader[2]);


      // Write the Dummy Mini Header into the file
      Int_t bHPosition = outfile.tellp();
      outfile.write((char*)(miniHeader),mHSize);


      for(Int_t ium = 0; ium < modulePerDDL; ium++)
	{
	  
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

      // Write real mini header
      // take the pointer to the beginning of the mini header
      // write the total number of words per ddl and bring the
      // pointer to the current file position and close it
      UInt_t cFPosition = outfile.tellp();
      sizeRawData = cFPosition - bHPosition - mHSize;
      outfile.seekp(bHPosition);
      outfile.write((char*)(&sizeRawData),sizeof(UInt_t));
      outfile.seekp(cFPosition);

      outfile.close();
    } // DDL Loop over


}
//____________________________________________________________________________

void AliPMDDDLRawData::ReadPMDRawData(Int_t evtno)
{
  // reads the raw data
  ifstream infile;

  Char_t filename[80];

  const Int_t kDDL = 6;
  Int_t imodule = 0;

  for(Int_t iddl = 0; iddl < kDDL; iddl++)
    {
      sprintf(filename,"Ev%dPMDddl%d.dat",evtno,iddl);
#ifndef __DECCXX
      infile.open(filename,ios::binary);
#else
      infile.open(filename);
#endif
      
      Int_t  ium;
      Int_t  irow;
      Int_t  icol;
      UInt_t baseword;
      UInt_t miniHeader[3];
      UInt_t sizeRawData;
      UInt_t mHSize = 3*sizeof(UInt_t);

      infile.read((char*)(miniHeader), mHSize);
      
      sizeRawData = miniHeader[0];
      Int_t totword = sizeRawData/4;
      
      UInt_t adc;
      UInt_t chno;
      UInt_t mcmno;
      
      for (Int_t iword = 0; iword < totword; iword++)
	{
	  infile.read((char*)(&baseword),sizeof(UInt_t));
	  
	  UnpackWord(0,11,adc,baseword);
	  UnpackWord(12,17,chno,baseword);
	  UnpackWord(18,28,mcmno,baseword);
	  
	  GetRowCol(iddl, mcmno, chno, ium, irow, icol);
	  if (iddl < 4)
	    {
	      imodule = iddl*6 + ium;
	    }
	  else if (iddl == 4)
	    {
	      imodule = 24 + ium;
	    }
	  else if (iddl == 5)
	    {
	      imodule = 36 + ium;
	    }

	  
	}

    } // ddl loop
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

      if(smn < 12)
	{
	  irownew = icol;
	  icolnew = irow;
	}
      else if( smn >= 12 && smn < 24)
	{
	  irownew = irow;
	  icolnew = icol;
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
  UInt_t ch[16][4] = { {3, 2, 29, 28},
		       {5, 1, 31, 27},
		       {6, 0, 30, 26},
		       {7, 4, 25, 24},
		       {8, 9, 20, 23},
		       {10, 14, 16, 22},
		       {11, 15, 17, 21},
		       {12, 13, 18, 19},
		       {35, 34, 61, 60},
		       {37, 33, 63, 59},
		       {38, 32, 62, 58},
		       {39, 36, 57, 56},
		       {40, 41, 52, 55},
		       {42, 46, 48, 54},
		       {43, 47, 49, 53},
		       {44, 45, 50, 51} };
  
  if (ddlno == 0 || ddlno == 1)
    {
      // PRE plane, SU Mod = 1, 2
      Int_t irownew = row%16;
      Int_t icolnew = col%4;
      Int_t irowdiv = row/16;
      Int_t icoldiv = col/4;

      mcmno = 72*um + 6*icoldiv + irowdiv;
      chno  = ch[irownew][icolnew];
    }
  else if (ddlno == 2 || ddlno == 3)
    {
      // PRE plane,  SU Mod = 3, 4
      Int_t irownew = row%16;
      Int_t icolnew = col%4;
      Int_t irowdiv = row/16;
      Int_t icoldiv = col/4;
      
      mcmno = 72*um + 3*icoldiv + irowdiv;
      chno  = ch[irownew][icolnew];
    }
  else if (ddlno == 4)
    {
      // CPV plane,  SU Mod = 1, 2
      Int_t irownew = row%16;
      Int_t icolnew = col%4;
      Int_t irowdiv = row/16;
      Int_t icoldiv = col/4;

      mcmno = 72*um + 6*icoldiv + irowdiv;
      chno  = ch[irownew][icolnew];
    }
  else if (ddlno == 5)
    {
      // CPV plane,  SU Mod = 3, 4
      Int_t irownew = row%16;
      Int_t icolnew = col%4;
      Int_t irowdiv = row/16;
      Int_t icoldiv = col/4;
      
      mcmno = 72*um + 3*icoldiv + irowdiv;
      chno  = ch[irownew][icolnew];
    }

}
//____________________________________________________________________________

void AliPMDDDLRawData::GetRowCol(Int_t ddlno, UInt_t mcmno, UInt_t chno,
				 Int_t &um, Int_t &row, Int_t &col)
{
  UInt_t ch[64] = { 9, 5, 1, 0, 13, 4, 8, 12,
		    16, 17, 20, 24, 28, 29, 21, 25,
		    22, 26, 30, 31, 18, 27, 23, 19,
		    15, 14, 11, 7, 3, 2, 10, 6,
		    41, 37, 33, 32, 45, 36, 40, 44,
		    48, 49, 52, 56, 60, 61, 53, 57,
		    54, 58, 62, 63, 50, 59, 55, 51,
		    47, 46, 43, 39, 35, 34, 42, 38 };

  if (ddlno == 1 || ddlno == 2)
    {
      um  = mcmno/72;
      Int_t mcmnonew = mcmno - 72*um;
      Int_t rowcol  = ch[chno];
      Int_t irownew = rowcol/4;
      Int_t icolnew = rowcol%4;
      
      Int_t remmcm  = mcmnonew%6;
      Int_t divmcm  = mcmnonew/6;
      
      row = 16*remmcm + irownew;
      col =  4*divmcm + icolnew;
    }
  else   if (ddlno == 3 || ddlno == 4)
    {
      um  = mcmno/72;
      Int_t mcmnonew = mcmno - 72*um;
      Int_t rowcol  = ch[chno];
      Int_t irownew = rowcol/4;
      Int_t icolnew = rowcol%4;
      
      Int_t remmcm  = mcmnonew%3;
      Int_t divmcm  = mcmnonew/3;
      
      row = 16*remmcm + irownew;
      col =  4*divmcm + icolnew;
    }
  else if (ddlno == 4)
    {
      um  = mcmno/144;
      Int_t mcmnonew = mcmno - 72*um;
      Int_t rowcol  = ch[chno];
      Int_t irownew = rowcol/4;
      Int_t icolnew = rowcol%4;
      
      Int_t remmcm  = mcmnonew%6;
      Int_t divmcm  = mcmnonew/6;
      
      row = 16*remmcm + irownew;
      col =  4*divmcm + icolnew;
    }
  else if (ddlno == 5)
    {
      um  = mcmno/144;
      Int_t mcmnonew = mcmno - 72*um;
      Int_t rowcol  = ch[chno];
      Int_t irownew = rowcol/4;
      Int_t icolnew = rowcol%4;
      
      Int_t remmcm  = mcmnonew%3;
      Int_t divmcm  = mcmnonew/3;
      
      row = 16*remmcm + irownew;
      col =  4*divmcm + icolnew;
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

