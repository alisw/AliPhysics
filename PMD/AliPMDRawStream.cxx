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

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to PMD digits in raw data.
///
/// It loops over all PMD digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE.
/// Several getters provide information about the current digit.
///
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TObjArray.h>
#include <TString.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliPMDBlockHeader.h"
#include "AliPMDDspHeader.h"
#include "AliPMDPatchBusHeader.h"
#include "AliPMDddldata.h"
#include "AliPMDRawStream.h"
#include "AliRawReader.h"

ClassImp(AliPMDRawStream)


//_____________________________________________________________________________
AliPMDRawStream::AliPMDRawStream(AliRawReader* rawReader) :
  fRawReader(rawReader)
{
// create an object to read PMD raw digits

  fRawReader->Select("PMD");
}

//_____________________________________________________________________________
AliPMDRawStream::AliPMDRawStream(const AliPMDRawStream& stream) :
  TObject(stream),
  fRawReader(NULL)
{
// copy constructor

  AliFatal("Copy constructor not implemented");
}

//_____________________________________________________________________________
AliPMDRawStream& AliPMDRawStream::operator = (const AliPMDRawStream& 
					      /* stream */)
{
// assignment operator

  AliFatal("operator = assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliPMDRawStream::~AliPMDRawStream()
{
// destructor

}


//_____________________________________________________________________________

Bool_t AliPMDRawStream::DdlData(Int_t indexDDL, TObjArray *pmdddlcont)
{
// read the next raw digit
// returns kFALSE if there is no digit left

  AliPMDddldata *pmdddldata;

  if (!fRawReader->ReadHeader()) return kFALSE;
  Int_t  iddl  = fRawReader->GetDDLID();
  Int_t dataSize = fRawReader->GetDataSize();
  Int_t totaldataword = dataSize/4;
  Int_t equipId = fRawReader->GetEquipmentId();

  if (dataSize <= 0) return kFALSE;
  if (indexDDL != iddl)
    {
      AliError("Mismatch in the DDL index");
      return kFALSE;
    }

  UInt_t *buffer;
  buffer = new UInt_t[totaldataword];
  UInt_t data;
  for (Int_t i = 0; i < totaldataword; i++)
    {
      fRawReader->ReadNextInt(data);
      buffer[i] = data;
    }

  // --- Open the mapping file

  TString fileName(gSystem->Getenv("ALICE_ROOT"));
  if(iddl == 0)
    {
      fileName += "/PMD/PMD_Mapping_ddl0.dat";
    }
  else if(iddl == 1)
    {
      fileName += "/PMD/PMD_Mapping_ddl1.dat";
    }
  else if(iddl == 2)
    {
      fileName += "/PMD/PMD_Mapping_ddl2.dat";
    }
  else if(iddl == 3)
    {
      fileName += "/PMD/PMD_Mapping_ddl3.dat";
    }
  else if(iddl == 4)
    {
      fileName += "/PMD/PMD_Mapping_ddl4.dat";
    }
  else if(iddl == 5)
    {
      fileName += "/PMD/PMD_Mapping_ddl5.dat";
    }

  ifstream infile;
  infile.open(fileName.Data(), ios::in); // ascii file
  if(!infile)
    AliError(Form("Could not read the mapping file for DDL No = %d",iddl));
  
  Int_t modulePerDDL = 0;
  if (iddl < 4)
    {
      modulePerDDL = 6;
    }
  else if (iddl == 4 || iddl == 5)
    {
      modulePerDDL = 12;
    }

  const Int_t kNPatchBus = 50;
  
  Int_t modno, totPatchBus, bPatchBus, ePatchBus;
  Int_t ibus, totmcm, rows, rowe, cols, cole;
  Int_t moduleNo[kNPatchBus], mcmperBus[kNPatchBus];
  Int_t startRowBus[kNPatchBus], endRowBus[kNPatchBus];
  Int_t startColBus[kNPatchBus], endColBus[kNPatchBus];


  for (Int_t ibus = 0; ibus < kNPatchBus; ibus++)
    {
      mcmperBus[ibus]   = -1;
      startRowBus[ibus] = -1;
      endRowBus[ibus]   = -1;
      startColBus[ibus] = -1;
      endColBus[ibus]   = -1;
    }


  for (Int_t im = 0; im < modulePerDDL; im++)
    {
      infile >> modno;
      infile >> totPatchBus >> bPatchBus >> ePatchBus;
      
      for(Int_t i=0; i<totPatchBus; i++)
	{
	  infile >> ibus >> totmcm >> rows >> rowe >> cols >> cole;

	  moduleNo[ibus]     = modno;
	  mcmperBus[ibus]    = totmcm;
	  startRowBus[ibus]  = rows;
	  endRowBus[ibus]    = rowe;
	  startColBus[ibus]  = cols;
	  endColBus[ibus]    = cole;
	}
    }



  infile.close();


  AliPMDBlockHeader    blockHeader;
  AliPMDDspHeader      dspHeader;
  AliPMDPatchBusHeader pbusHeader;

  const Int_t kblHLen   = blockHeader.GetHeaderLength();
  const Int_t kdspHLen  = dspHeader.GetHeaderLength();
  const Int_t kpbusHLen = pbusHeader.GetHeaderLength();
  

  Int_t idet, ismn;
  Int_t irow = -1;
  Int_t icol = -1;

  Int_t blHeaderWord[8];
  Int_t dspHeaderWord[10];
  Int_t pbusHeaderWord[4];

  Int_t ilowLimit = 0;
  Int_t iuppLimit = 0;

  Int_t blRawDataLength = 0;
  Int_t iwordcount = 0;


  for (Int_t iblock = 0; iblock < 2; iblock++)
    {
      ilowLimit = iuppLimit;
      iuppLimit = ilowLimit + kblHLen;


      for (Int_t i = ilowLimit; i < iuppLimit; i++)
	{
	  blHeaderWord[i-ilowLimit] = (Int_t) buffer[i];
	}

      blockHeader.SetHeader(blHeaderWord);

      blRawDataLength = blockHeader.GetRawDataLength();

      for (Int_t idsp = 0; idsp < 5; idsp++)
	{
	  ilowLimit = iuppLimit;
	  iuppLimit = ilowLimit + kdspHLen;

	  for (Int_t i = ilowLimit; i < iuppLimit; i++)
	    {
	      iwordcount++;
	      dspHeaderWord[i-ilowLimit] = (Int_t) buffer[i];
	    }
	  dspHeader.SetHeader(dspHeaderWord);

	  for (Int_t ibus = 0; ibus < 5; ibus++)
	    {

	      ilowLimit = iuppLimit;
	      iuppLimit = ilowLimit + kpbusHLen;

	      for (Int_t i = ilowLimit; i < iuppLimit; i++)
		{
		  iwordcount++;
		  pbusHeaderWord[i-ilowLimit] = (Int_t) buffer[i];
		}
	      pbusHeader.SetHeader(pbusHeaderWord);
	      Int_t rawdatalength = pbusHeader.GetRawDataLength();
	      Int_t pbusid = pbusHeader.GetPatchBusId();

	      ilowLimit = iuppLimit;
	      iuppLimit = ilowLimit + rawdatalength;

	      Int_t imodule = moduleNo[pbusid];


	      for (Int_t iword = ilowLimit; iword < iuppLimit; iword++)
		{
		  iwordcount++;
		  data = buffer[iword];

		  Int_t isig =  data & 0x0FFF;
		  Int_t ich  = (data >> 12) & 0x003F;
		  Int_t imcm = (data >> 18) & 0x07FF;
		  Int_t ibit = (data >> 31) & 0x0001;

		  GetRowCol(iddl, pbusid, imcm, ich, 
			    startRowBus, endRowBus,
			    startColBus, endColBus,
			    irow, icol);

		  ConvertDDL2SMN(iddl, imodule, ismn, idet);
		  TransformH2S(ismn, irow, icol);

		  pmdddldata = new AliPMDddldata();

		  pmdddldata->SetDetector(idet);
		  pmdddldata->SetSMN(ismn);
		  pmdddldata->SetModule(imodule);
		  pmdddldata->SetPatchBusId(pbusid);
		  pmdddldata->SetMCM(imcm);
		  pmdddldata->SetChannel(ich);
		  pmdddldata->SetRow(irow);
		  pmdddldata->SetColumn(icol);
		  pmdddldata->SetSignal(isig);
		  pmdddldata->SetParityBit(ibit);
		  
		  pmdddlcont->Add(pmdddldata);
		  
		} // data word loop

	      if (iwordcount == blRawDataLength) break;

	    } // patch bus loop

	  if (dspHeader.GetPaddingWord() == 1) iuppLimit++;
	  if (iwordcount == blRawDataLength) break;

	} // end of DSP
      if (iwordcount == blRawDataLength) break;

    } // end of BLOCK

  
  delete [] buffer;

  return kTRUE;
}
//_____________________________________________________________________________
void AliPMDRawStream::GetRowCol(Int_t ddlno, Int_t pbusid,
				UInt_t mcmno, UInt_t chno,
				Int_t startRowBus[], Int_t endRowBus[],
				Int_t startColBus[], Int_t endColBus[],
				Int_t &row, Int_t &col) const
{
// decode: ddlno, patchbusid, mcmno, chno -> um, row, col


  static const UInt_t kCh[64] = { 53, 58, 57, 54, 61, 62, 60, 63,
				  49, 59, 56, 55, 52, 50, 48, 51,
				  44, 47, 45, 43, 40, 39, 36, 46,
				  32, 35, 33, 34, 41, 38, 37, 42,
				  21, 26, 25, 22, 29, 30, 28, 31,
				  17, 27, 24, 23, 20, 18, 16, 19,
				  12, 15, 13, 11,  8,  7,  4, 14,
				  0,   3,  1,  2,  9,  6,  5, 10 };


  Int_t rowcol  = kCh[chno];
  Int_t irownew = rowcol/4;
  Int_t icolnew = rowcol%4;

  if (ddlno == 0)
    {
      if (pbusid  <= 2)
	{
	  if (mcmno >= 12)
	    {
	      row = startRowBus[pbusid] + irownew;
	      col = startColBus[pbusid] + (mcmno-12)*4 + icolnew;
	    }
	  else
	    {
	      // Add 16 to skip the 1st 15 rows
	      row = startRowBus[pbusid] + irownew + 16;
	      col = startColBus[pbusid] + mcmno*4 + icolnew;
	    }
	}
      else if (pbusid > 2)
	{
	  row = startRowBus[pbusid] + irownew;
	  col = startColBus[pbusid] + mcmno*4 + icolnew;
	  
	}
    }
  else if (ddlno == 1)
    {
      if (pbusid  <= 2)
	{
	  if (mcmno >= 12)
	    {
	      row = endRowBus[pbusid] - (15 - irownew);
	      col = startColBus[pbusid] + (mcmno-12)*4 + icolnew;
	    }
	  else
	    {
	      // Subtract 16 to skip the 1st 15 rows
	      row = endRowBus[pbusid] - 16 - (15 - irownew) ;
	      col = startColBus[pbusid] + mcmno*4 + icolnew;
	    }
	}
      else if (pbusid > 2)
	{
	  row = endRowBus[pbusid] - (15 - irownew);
	  col = startColBus[pbusid] + mcmno*4 + icolnew;
	}
    }
  else if (ddlno == 2)
    {
      row = startRowBus[pbusid] + irownew;
      col = endColBus[pbusid] - mcmno*4 - (3 - icolnew);
    }
  else if (ddlno == 3)
    {
      row = endRowBus[pbusid] - (15 - irownew);
      col = endColBus[pbusid] - mcmno*4 - (3 - icolnew);
    }
  else if (ddlno == 4)
    {
      if (pbusid  <= 16)
	{
	  if (mcmno >= 12)
	    {
	      row = startRowBus[pbusid] + irownew;
	      col = startColBus[pbusid] + (mcmno-12)*4 + icolnew;
	    }
	  else
	    {
	      // Add 16 to skip the 1st 15 rows
	      row = startRowBus[pbusid] + irownew + 16;
	      col = startColBus[pbusid] + mcmno*4 + icolnew;
	    }
	}
      else if (pbusid > 16 && pbusid <= 20)
	{
	  row = startRowBus[pbusid] + irownew;
	  col = startColBus[pbusid] + mcmno*4 + icolnew;
	  
	}
      else if(pbusid > 20)
	{
	  row = endRowBus[pbusid] - (15 - irownew);
	  col = endColBus[pbusid] - mcmno*4 - (3 - icolnew);
	}
    }
  else if (ddlno == 5)
    {
      if (pbusid  <= 16)
	{
	  if (mcmno >= 12)
	    {
	      row = endRowBus[pbusid] - (15 - irownew);
	      col = startColBus[pbusid] + (mcmno-12)*4 + icolnew;
	    }
	  else
	    {
	      // Subtract 16 to skip the 1st 15 rows
	      row = endRowBus[pbusid] - 16 - (15 - irownew) ;
	      col = startColBus[pbusid] + mcmno*4 + icolnew;
	    }
	}
      else if (pbusid > 16 && pbusid <= 20)
	{
	  row = endRowBus[pbusid] - (15 - irownew);
	  col = startColBus[pbusid] + mcmno*4 + icolnew;
	}
      else if (pbusid > 20)
	{
	  row = startRowBus[pbusid] + irownew;
	  col = endColBus[pbusid] - mcmno*4 - (3 - icolnew);
	}
    }
}
//_____________________________________________________________________________
void AliPMDRawStream::ConvertDDL2SMN(Int_t iddl, Int_t imodule,
				     Int_t &smn, Int_t &detector) const
{
  // This converts the DDL number (0 to 5), Module Number (0-47)
  // to Serial module number in one detector (SMN : 0-23) and
  // detector number (0:PRE plane, 1:CPV plane)
  if (iddl < 4)
    {
      smn = imodule;
      detector = 0;
    }
  else
    {
      smn = imodule - 24;
      detector = 1;
    }
}
//_____________________________________________________________________________
void AliPMDRawStream::TransformH2S(Int_t smn, Int_t &row, Int_t &col) const
{
  // This does the transformation of the hardware coordinate to
  // software 
  // i.e., For SuperModule 0 &1, instead of 96x48(hardware),
  // it is 48x96 (software)
  // For Supermodule 3 & 4, 48x96

  Int_t irownew  = 0;
  Int_t icolnew  = 0;

  if(smn < 12)
    {
      irownew = col;
      icolnew = row;
    }
  else if(smn >= 12 && smn < 24)
    {
      irownew = row;
      icolnew = col;
    }

  row = irownew;
  col = icolnew;
}
//_____________________________________________________________________________
