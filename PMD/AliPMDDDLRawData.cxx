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
#include <TString.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliRawDataHeaderSim.h"
#include "AliBitPacking.h"
#include "AliPMDdigit.h"
#include "AliPMDBlockHeader.h"
#include "AliPMDDspHeader.h"
#include "AliPMDPatchBusHeader.h"
#include "AliPMDRawStream.h"
#include "AliPMDDDLRawData.h"
#include "AliDAQ.h"
#include "AliFstream.h"

ClassImp(AliPMDDDLRawData)

AliPMDDDLRawData::AliPMDDDLRawData():
  fDigits(new TClonesArray("AliPMDdigit", 1000))
{
  // Default Constructor
  //

}
//____________________________________________________________________________
AliPMDDDLRawData::AliPMDDDLRawData(const AliPMDDDLRawData& ddlraw):
  TObject(ddlraw),
  fDigits(ddlraw.fDigits)
{
  //Copy Constructor 
}
//____________________________________________________________________________
AliPMDDDLRawData & AliPMDDDLRawData::operator=(const AliPMDDDLRawData& ddlraw)
{
  //Assignment operator 
  if(this != &ddlraw)
    {
      fDigits = ddlraw.fDigits;
    }
  return *this;
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

  AliFstream *outfile;

  TBranch *branch = treeD->GetBranch("PMDDigit");
  if (!branch)
    {
      AliError("PMD Digit branch not found");
      return;
    }
  branch->SetAddress(&fDigits);  
  
  Int_t   nmodules = (Int_t) treeD->GetEntries();
  AliDebug(1,Form("Number of modules inside treeD = %d",nmodules));

  const Int_t kDDL          = AliDAQ::NumberOfDdls("PMD");
  Int_t modulePerDDL        = 0;


  AliRawDataHeaderSim header;
  UInt_t sizeRawData = 0;
  
  const Int_t kbusSize = 51;
  const Int_t kSize = 1536;
  UInt_t buffer[kSize];

  UInt_t busPatch[kbusSize][1536];

  Int_t contentsBus[kbusSize];

  Char_t filename[80];


  Int_t mmodule = 0;
  for(Int_t iddl = 0; iddl < kDDL; iddl++)
    {
      strcpy(filename,AliDAQ::DdlFileName("PMD",iddl));
      
      outfile = new AliFstream(filename);
      
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
      Int_t bHPosition = outfile->Tellp();
      outfile->WriteBuffer((char*)(&header),sizeof(header));

      for (Int_t ibus = 0; ibus < kbusSize; ibus++)
	{
	  contentsBus[ibus] = 0;
	  for (Int_t ich = 0; ich < kSize; ich++)
	    {
	      busPatch[ibus][ich] = 0;
	    }
	}

      for(Int_t ium = 0; ium < modulePerDDL; ium++)
	{
	  if (iddl == 4 && ium == 6) mmodule = 42;

	  // Extract energy deposition per cell and pack it
	  // in a 32-bit word and returns all the total words
	  // per one unit-module
	  
	  GetUMDigitsData(treeD, mmodule, iddl, contentsBus, busPatch);
	  mmodule++;
	}

      

      Int_t ij = 0;
      Int_t dsp[10];
      Int_t dspBus[10];
      for (Int_t i = 0; i < 10; i++)
	{
	  dsp[i] = 0;
	  dspBus[i] = 0;
	  for (Int_t ibus=0; ibus < 5; ibus++)
	    {
	      ij++;
	      if (contentsBus[ij] > 0)
		{
		  dsp[i] += contentsBus[ij];
		}
	      dspBus[i]++;
	    }
	  // Add the patch Bus header to the DSP contents
	  dsp[i] += 4*dspBus[i];
	}

      Int_t dspBlockARDL    = 0;
      Int_t dspBlockBRDL    = 0;
      Int_t remainder       = 0;


      for (Int_t i = 0; i < 5; i++)
	{
	  Int_t ieven = 2*i;
	  Int_t iodd  = 2*i + 1;
	  if (dsp[ieven] > 0)
	    {
	      dspBlockARDL += dsp[ieven];
	      remainder = dsp[ieven]%2;
	      if (remainder == 1)
		{
		  dspBlockARDL++;
		}
	    }
	  if (dsp[iodd] > 0)
	    {
	      dspBlockBRDL += dsp[iodd];
	      remainder = dsp[iodd]%2;
	      if (remainder == 1)
		{
		  dspBlockBRDL++;
		}
	    }
	}

      dspBlockARDL += 50;
      dspBlockBRDL += 50;

      // Start writing the DDL file

      AliPMDBlockHeader blHeader;
      AliPMDDspHeader   dspHeader;
      AliPMDPatchBusHeader pbusHeader;

      const Int_t kblHLen   = blHeader.GetHeaderLength();
      const Int_t kdspHLen  = dspHeader.GetHeaderLength();
      const Int_t kpbusHLen = pbusHeader.GetHeaderLength();

      UInt_t dspRDL = 0;
      UInt_t dspBlockHeaderWord[8];
      UInt_t dspHeaderWord[10];
      UInt_t patchBusHeaderWord[4];
      Int_t  iskip[5];
      UInt_t ddlEndWord[2] = {0xDEADFACE, 0xDEADFACE};

      Int_t bknJunk = 0;


      for (Int_t iblock = 0; iblock < 2; iblock++)
	{
	  // DSP Block Header
	  
	  for (Int_t i=0; i<kblHLen; i++)
	    {
	      dspBlockHeaderWord[i] = 0;
	    }
	  if (iblock == 0)
	    {
	      dspBlockHeaderWord[1] = (UInt_t) (dspBlockARDL + kblHLen);
	      dspBlockHeaderWord[2] = (UInt_t) dspBlockARDL;
	    }
	  else if (iblock == 1)
	    {
	      dspBlockHeaderWord[1] = (UInt_t) (dspBlockBRDL + kblHLen);
	      dspBlockHeaderWord[2] = (UInt_t) dspBlockBRDL;
	    }

	  outfile->WriteBuffer((char*)dspBlockHeaderWord,kblHLen*sizeof(UInt_t));

	  if (iblock == 0)
	    {
	      iskip[0] = 0;
	      iskip[1] = 10;
	      iskip[2] = 20;
	      iskip[3] = 30;
	      iskip[4] = 40;
	    }
	  else if (iblock == 1)
	    {
	      iskip[0] = 5;
	      iskip[1] = 15;
	      iskip[2] = 25;
	      iskip[3] = 35;
	      iskip[4] = 45;
	    }

	  for (Int_t idsp = 0; idsp < 5; idsp++)
	    {
	      // DSP Header
	      Int_t dspno = 0;
	      if (iblock == 0)
		{
		  dspno = 2*idsp;
		  dspRDL = (UInt_t) dsp[dspno];
		}
	      else if (iblock == 1)
		{
		  dspno = 2*idsp + 1;
		  dspRDL = (UInt_t) dsp[dspno];
		}

	      for (Int_t i=0; i<kdspHLen; i++)
		{
		  dspHeaderWord[i] = 0;
		}
	      remainder = dspRDL%2;
	      if (remainder == 1) dspRDL++;

	      dspHeaderWord[1] = dspRDL + kdspHLen;
	      dspHeaderWord[2] = dspRDL;
	      dspHeaderWord[3] = dspno;
	      if (remainder == 1) dspHeaderWord[8] = 1; // setting the padding word


	      outfile->WriteBuffer((char*)dspHeaderWord,kdspHLen*sizeof(UInt_t));

	      for (Int_t ibus = 0; ibus < 5; ibus++)
		{
		  // Patch Bus Header

		  Int_t busno = iskip[idsp] + ibus + 1;
		  Int_t patchbusRDL = contentsBus[busno];

		  if (patchbusRDL > 0)
		    {
		      patchBusHeaderWord[0] = 0;
		      patchBusHeaderWord[1] = (UInt_t) (patchbusRDL + kpbusHLen);
		      patchBusHeaderWord[2] = (UInt_t) patchbusRDL;
		      patchBusHeaderWord[3] = (UInt_t) busno;
		    }
		  else if (patchbusRDL == 0)
		    {
		      patchBusHeaderWord[0] = 0;
		      patchBusHeaderWord[1] = (UInt_t) kpbusHLen;
		      patchBusHeaderWord[2] = (UInt_t) 0;
		      patchBusHeaderWord[3] = (UInt_t) busno;
		    }


		  outfile->WriteBuffer((char*)patchBusHeaderWord,4*sizeof(UInt_t));

		  bknJunk += patchbusRDL;

		  for (Int_t iword = 0; iword < patchbusRDL; iword++)
		    {
		      buffer[iword] = busPatch[busno][iword];
		    }
		  
		  outfile->WriteBuffer((char*)buffer,patchbusRDL*sizeof(UInt_t));

		} // End of patch bus loop


	      // Adding a padding word if the total words odd
	      if (remainder == 1)
		{
		  UInt_t paddingWord = dspHeader.GetDefaultPaddingWord();
		  outfile->WriteBuffer((char*)(&paddingWord),sizeof(UInt_t));
		}
	    }
	}

      // Write two extra word at the end of each DDL file
      outfile->WriteBuffer((char*)ddlEndWord,2*sizeof(UInt_t));

      // Write real data header
      // take the pointer to the beginning of the data header
      // write the total number of words per ddl and bring the
      // pointer to the current file position and close it
      UInt_t cFPosition = outfile->Tellp();
      sizeRawData = cFPosition - bHPosition - sizeof(header);

      header.fSize = cFPosition - bHPosition;
      header.SetAttribute(0);  // valid data
      outfile->Seekp(bHPosition);
      outfile->WriteBuffer((char*)(&header),sizeof(header));
      outfile->Seekp(cFPosition);

      delete outfile;
    } // DDL Loop over


}
//____________________________________________________________________________
void AliPMDDDLRawData::GetUMDigitsData(TTree *treeD, Int_t imodule,
				       Int_t ddlno, Int_t *contentsBus,
				       UInt_t busPatch[][1536])
{
  // Retrives digits data UnitModule by UnitModule

  UInt_t baseword;
  UInt_t mcmno, chno;
  UInt_t adc;
  Int_t  det, smn, irow, icol;
  Int_t  parity;
  
  const Int_t kMaxBus = 51;   // BKN
  Int_t totPatchBus, bPatchBus, ePatchBus;
  Int_t ibus, totmcm, rows, cols, rowe, cole;
  Int_t moduleno;
  Int_t busno = 0;
  Int_t patchBusNo[kMaxBus], mcmperBus[kMaxBus];
  Int_t startRowBus[kMaxBus], startColBus[kMaxBus];
  Int_t endRowBus[kMaxBus], endColBus[kMaxBus];

  Int_t beginPatchBus = -1;
  Int_t endPatchBus   = -1;
  for(Int_t i = 0; i < kMaxBus; i++)
    {
      patchBusNo[i]  = -1;
      mcmperBus[i]   = -1;
      startRowBus[i] = -1;
      startColBus[i] = -1;
      endRowBus[i]   = -1;
      endColBus[i]   = -1;
    }
  Int_t modulePerDDL = 0;
  if (ddlno < 4)
    {
      modulePerDDL = 6;
    }
  else if (ddlno == 4 || ddlno == 5)
    {
      modulePerDDL = 12;
    }

  TString fileName(gSystem->Getenv("ALICE_ROOT"));

  if(ddlno == 0)
    {
      fileName += "/PMD/PMD_Mapping_ddl0.dat";
    }
  else if(ddlno == 1)
    {
      fileName += "/PMD/PMD_Mapping_ddl1.dat";
    }
  else if(ddlno == 2)
    {
      fileName += "/PMD/PMD_Mapping_ddl2.dat";
    }
  else if(ddlno == 3)
    {
      fileName += "/PMD/PMD_Mapping_ddl3.dat";
    }
  else if(ddlno == 4)
    {
      fileName += "/PMD/PMD_Mapping_ddl4.dat";
    }
  else if(ddlno == 5)
    {
      fileName += "/PMD/PMD_Mapping_ddl5.dat";
    }

  ifstream infile;
  infile.open(fileName.Data(), ios::in); // ascii file
  if(!infile)
    AliError(Form("Could not read the mapping file for DDL No = %d",ddlno));

  for (Int_t im = 0; im < modulePerDDL; im++)
    {
      infile >> moduleno;
      infile >> totPatchBus >> bPatchBus >> ePatchBus;

      if (moduleno == imodule)
	{
	  beginPatchBus = bPatchBus;
	  endPatchBus   = ePatchBus;
	}
      
      for(Int_t i=0; i<totPatchBus; i++)
	{
	  infile >> ibus >> totmcm >> rows >> rowe >> cols >> cole;

	  if (moduleno == imodule)
	    {
	      patchBusNo[ibus]   = ibus;
	      mcmperBus[ibus]    = totmcm;
	      startRowBus[ibus]  = rows;
	      startColBus[ibus]  = cols;
	      endRowBus[ibus]    = rowe;
	      endColBus[ibus]    = cole;
	    }
	}

    }

  infile.close();

  treeD->GetEntry(imodule); 
  Int_t nentries = fDigits->GetLast();
  Int_t totword = nentries+1;

  AliPMDdigit *pmddigit;

  for (Int_t ient = 0; ient < totword; ient++)
    {
      pmddigit = (AliPMDdigit*)fDigits->UncheckedAt(ient);
      
      det    = pmddigit->GetDetector();
      smn    = pmddigit->GetSMNumber();
      irow   = pmddigit->GetRow();
      icol   = pmddigit->GetColumn();
      adc    = (UInt_t) pmddigit->GetADC();

      TransformS2H(smn,irow,icol);

      GetMCMCh(ddlno, smn, irow, icol, beginPatchBus, endPatchBus,
	       mcmperBus, startRowBus, startColBus,
	       endRowBus, endColBus, busno, mcmno, chno);

      baseword = 0;
      AliBitPacking::PackWord(adc,baseword,0,11);
      AliBitPacking::PackWord(chno,baseword,12,17);
      AliBitPacking::PackWord(mcmno,baseword,18,28);
      AliBitPacking::PackWord(0,baseword,29,30);
      parity = ComputeParity(baseword);      // generate the parity bit
      AliBitPacking::PackWord(parity,baseword,31,31);

      Int_t jj = contentsBus[busno];
      busPatch[busno][jj] = baseword;

      contentsBus[busno]++;
    }

}

//____________________________________________________________________________
void AliPMDDDLRawData::TransformS2H(Int_t smn, Int_t &irow, Int_t &icol)
{
  // Does the Software to Hardware coordinate transformation
  //

  Int_t  irownew = 0;
  Int_t  icolnew = 0;

  // First in digits we have all dimension 48x96
  // Transform into the realistic one, i.e, For SM 0&1 96(row)x48(col)
  // and for SM 2&3 48(row)x96(col)
  // 
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

  irow = irownew;
  icol = icolnew;

}


//____________________________________________________________________________

void AliPMDDDLRawData::GetMCMCh(Int_t ddlno, Int_t smn, Int_t row, Int_t col,
				Int_t beginPatchBus, Int_t endPatchBus,
				Int_t *mcmperBus,
				Int_t *startRowBus, Int_t *startColBus,
				Int_t *endRowBus, Int_t *endColBus,
				Int_t & busno, UInt_t &mcmno, UInt_t &chno)
{
  // This converts row col to hardware channel number
  // This is the final version of mapping supplied by Mriganka

    UInt_t iCh[16][4];

    static const UInt_t kChDdl01[16][4] = { {6, 4, 5, 7},
					   {10, 2, 1, 9},
					   {12, 0, 3, 11},
					   {14, 8, 13, 15},
					   {16, 18, 23, 17},
					   {20, 28, 31, 19},
					   {22, 30, 29, 21},
					   {24, 26, 27, 25},
					   {38, 36, 37, 39},
					   {42, 34, 33, 41},
					   {44, 32, 35, 43},
					   {46, 40, 45, 47},
					   {48, 50, 55, 49},
					   {52, 60, 63, 51},
					   {54, 62, 61, 53},
					   {56, 58, 59, 57} };


    static const UInt_t kChDdl23[16][4] = { {57, 59, 58, 56},
					    {53, 61, 62, 54},
					    {51, 63, 60, 52},
					    {49, 55, 50, 48},
					    {47, 45, 40, 46},
					    {43, 35, 32, 44},
					    {41, 33, 34, 42},
					    {39, 37, 36, 38},
					    {25, 27, 26, 24},
					    {21, 29, 30, 22},
					    {19, 31, 28, 20},
					    {17, 23, 18, 16},
					    {15, 13, 8, 14},
					    {11, 3, 0, 12},
					    {9, 1, 2, 10},
					    {7, 5, 4, 6} };
    
    
    static const UInt_t kChDdl41[16][4] = { {56, 58, 59, 57},
					   {54, 62, 61, 53},
					   {52, 60, 63, 51},
					   {48, 50, 55, 49},
					   {46, 40, 45, 47},
					   {44, 32, 35, 43},
					   {42, 34, 33, 41},
					   {38, 36, 37, 39},
					   {24, 26, 27, 25},
					   {22, 30, 29, 21},
					   {20, 28, 31, 19},
					   {16, 18, 23, 17},
					   {14, 8, 13, 15},
					   {12, 0, 3, 11},
					   {10, 2, 1, 9},
					   {6, 4, 5, 7} };


    static const UInt_t kChDdl42[16][4] = { {7, 5, 4, 6},
					    {9, 1, 2, 10},
					    {11, 3, 0, 12},
					    {15, 13, 8, 14},
					    {17, 23, 18, 16},
					    {19, 31, 28, 20},
					    {21, 29, 30, 22},
					    {25, 27, 26, 24},
					    {39, 37, 36, 38},
					    {41, 33, 34, 42},
					    {43, 35, 32, 44},
					    {47, 45, 40, 46},
					    {49, 55, 50, 48},
					    {51, 63, 60, 52},
					    {53, 61, 62, 54},
					    {57, 59, 58, 56} };


    static const UInt_t kChDdl51[16][4] = { {7, 5, 4, 6},
					    {9, 1, 2, 10},
					    {11, 3, 0, 12},
					    {15, 13, 8, 14},
					    {17, 23, 18, 16},
					    {19, 31, 28, 20},
					    {21, 29, 30, 22},
					    {25, 27, 26, 24},
					    {39, 37, 36, 38},
					    {41, 33, 34, 42},
					    {43, 35, 32, 44},
					    {47, 45, 40, 46},
					    {49, 55, 50, 48},
					    {51, 63, 60, 52},
					    {53, 61, 62, 54},
					    {57, 59, 58, 56} };
    


    static const UInt_t kChDdl52[16][4] = { {56, 58, 59, 57},
					    {54, 62, 61, 53},
					    {52, 60, 63, 51},
					    {48, 50, 55, 49},
					    {46, 40, 45, 47},
					    {44, 32, 35, 43},
					    {42, 34, 33, 41},
					    {38, 36, 37, 39},
					    {24, 26, 27, 25},
					    {22, 30, 29, 21},
					    {20, 28, 31, 19},
					    {16, 18, 23, 17},
					    {14, 8, 13, 15},
					    {12, 0, 3, 11},
					    {10, 2, 1, 9},
					    {6, 4, 5, 7} };
    
    
    for (Int_t i = 0; i < 16; i++)
    {
	for (Int_t j = 0; j < 4; j++)
	{
	    if (ddlno == 0 || ddlno == 1) iCh[i][j] = kChDdl01[i][j];
	    if (ddlno == 2 || ddlno == 3) iCh[i][j] = kChDdl23[i][j];
	    
	    if (ddlno == 4 && smn < 6)                iCh[i][j] = kChDdl41[i][j];
	    if (ddlno == 4 && (smn >= 18 && smn < 24))iCh[i][j] = kChDdl42[i][j];
	    if (ddlno == 5 && (smn >= 12 && smn < 18))iCh[i][j] = kChDdl51[i][j];
	    if (ddlno == 5 && (smn >=  6 && smn < 12))iCh[i][j] = kChDdl52[i][j];
	}
    }


  Int_t irownew = row%16;
  Int_t icolnew = col%4;
  
  chno  = iCh[irownew][icolnew];
  
  
  for (Int_t ibus = beginPatchBus; ibus <= endPatchBus; ibus++)
    {
      Int_t srow = startRowBus[ibus];
      Int_t erow = endRowBus[ibus];
      Int_t scol = startColBus[ibus];
      Int_t ecol = endColBus[ibus];
      Int_t tmcm = mcmperBus[ibus];
      if ((row >= srow && row <= erow) && (col >= scol && col <= ecol))
	{
	  busno = ibus;
	  
	  // Find out the MCM Number
	  //
	  
	  if (ddlno == 0 || ddlno == 1)
	    {
	      // PRE plane, SU Mod = 0, 1
	      mcmno = (col-scol)/4 + 1;
	      
	    }
	  else if (ddlno == 2 || ddlno == 3)
	    {
	      // PRE plane,  SU Mod = 2, 3
	      icolnew = (col - scol)/4;
	      mcmno = tmcm - icolnew;
	    }
	  else if (ddlno == 4 )
	    {
	      // CPV plane,  SU Mod = 0, 3 : ddl = 4
	      
	      if(ibus <= 18)
		{
		  Int_t midrow = srow + 16;
		  if(row >= srow && row < midrow)
		    {
		      mcmno = 12 + (col-scol)/4 + 1;
		    }
		  else if(row >= midrow && row <= erow)
		  
		    {
		      mcmno = (col-scol)/4 + 1;
		    }
		}
	      
	      else if (ibus > 18)
		{
		  Int_t rowdiff = endRowBus[ibus] - startRowBus[ibus];
		  if(rowdiff > 16)
		    {
		      Int_t midrow = srow + 16;
		      if (row >= midrow && row <= erow)
			{
			  mcmno = 12 + (ecol -col)/4 + 1;
			}
		      else if (row >= srow && row < midrow)
			{
			  mcmno = (ecol - col)/4 + 1;
			}
		    }
		  else if (rowdiff < 16)
		    {
		      mcmno = (ecol - col)/4 + 1;
		    }
		}
	    }
	  else if ( ddlno == 5)
	    {
	      // CPV plane,  SU Mod = 1, 2 : ddl = 5
	      
	      if(ibus <= 18)
		{
		  Int_t midrow = srow + 16;
		  if(row >= srow && row < midrow)
		    {
		      mcmno = 12 + (col-scol)/4 + 1;
		    }
		  else if(row >= midrow && row <= erow)
		    {
		      mcmno = (col-scol)/4 + 1;
		    }
		}
	      
	      else if (ibus > 18)
		{
		  Int_t rowdiff = endRowBus[ibus] - startRowBus[ibus];
		  if(rowdiff > 16)
		    {
		      Int_t midrow = srow + 16;
		      if (row >= midrow && row <= erow)
			{
			  mcmno = 12 + (ecol - col)/4 + 1;
			}
		      else if (row >= srow && row < midrow)
			{
			  mcmno = (ecol - col)/4 + 1;
			}
		    }
		  else if (rowdiff < 16)
		    {
		      mcmno = (ecol - col)/4 + 1;
		    }
		}
	    }
	}
    }
} 

//____________________________________________________________________________

Int_t AliPMDDDLRawData::ComputeParity(UInt_t baseword)
{
  // Generate the parity bit

  Int_t count = 0;
  for(Int_t j=0; j<29; j++)
    {
      if (baseword & 0x01 ) count++;
      baseword >>= 1;
    }
  Int_t parity = count%2;
  return parity;
}

//____________________________________________________________________________
