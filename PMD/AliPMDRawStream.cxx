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
    fRawReader(rawReader),
    fData(NULL),
    fPosition(-1)
{
// create an object to read PMD raw digits

  fRawReader->Reset();
  fRawReader->Select("PMD");
}

//_____________________________________________________________________________
AliPMDRawStream::AliPMDRawStream(const AliPMDRawStream& stream) :
  TObject(stream),
  fRawReader(NULL),
  fData(NULL),
  fPosition(-1)
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

Int_t AliPMDRawStream::DdlData(TObjArray *pmdddlcont)
{
// read the next raw digit
// returns kFALSE if there is no digit left

    

    Int_t iddl = -1;

  AliPMDddldata *pmdddldata;

  if (!fRawReader->ReadHeader()) return iddl;

  iddl           = fRawReader->GetDDLID();
  Int_t dataSize = fRawReader->GetDataSize();
  Int_t totaldataword = dataSize/4;

  if (dataSize <= 0) return -1;

  UInt_t data;

  fRawReader->ReadNextData(fData);

  fPosition = 0;


  Int_t ibus;

  const Int_t kNPatchBus = 51;

  Int_t moduleNo[kNPatchBus], mcmperBus[kNPatchBus];
  Int_t startRowBus[kNPatchBus], endRowBus[kNPatchBus];
  Int_t startColBus[kNPatchBus], endColBus[kNPatchBus];

  for (ibus = 0; ibus < kNPatchBus; ibus++)
    {
      moduleNo[ibus]    = -1;
      mcmperBus[ibus]   = -1;
      startRowBus[ibus] = -1;
      endRowBus[ibus]   = -1;
      startColBus[ibus] = -1;
      endColBus[ibus]   = -1;
    }

  // Call the Mapping methods - hard coded


  if (iddl == 0)
  {
      Ddl0Mapping(moduleNo, mcmperBus, startRowBus, 
		  endRowBus, startColBus, endColBus);
  }
  else if (iddl == 1)
  {
      Ddl1Mapping(moduleNo, mcmperBus, startRowBus, 
		  endRowBus, startColBus, endColBus);
  }
  else if (iddl == 2)
  {
      Ddl2Mapping(moduleNo, mcmperBus, startRowBus, 
		  endRowBus, startColBus, endColBus);
  }
  else if (iddl == 3)
  {
      Ddl3Mapping(moduleNo, mcmperBus, startRowBus, 
		  endRowBus, startColBus, endColBus);
  }
  else if (iddl == 4)
  {
      Ddl4Mapping(moduleNo, mcmperBus, startRowBus, 
		  endRowBus, startColBus, endColBus);
  }
  else if (iddl == 5)
  {
      Ddl5Mapping(moduleNo, mcmperBus, startRowBus, 
		  endRowBus, startColBus, endColBus);
  }


  AliPMDBlockHeader    blockHeader;
  AliPMDDspHeader      dspHeader;
  AliPMDPatchBusHeader pbusHeader;

  const Int_t kblHLen   = blockHeader.GetHeaderLength();
  const Int_t kdspHLen  = dspHeader.GetHeaderLength();
  const Int_t kpbusHLen = pbusHeader.GetHeaderLength();
  
  Int_t parity;
  Int_t idet, ismn;
  Int_t irow = -1;
  Int_t icol = -1;

  Int_t blHeaderWord[8];
  Int_t dspHeaderWord[10];
  Int_t pbusHeaderWord[4];

  Int_t ilowLimit        = 0;
  Int_t iuppLimit        = 0;
  Int_t blRawDataLength  = 0;
  Int_t dspRawDataLength = 0;
  Int_t iwordddl         = 2;



  for (Int_t iblock = 0; iblock < 2; iblock++)
    {
      ilowLimit = iuppLimit;
      iuppLimit = ilowLimit + kblHLen;

      for (Int_t i = ilowLimit; i < iuppLimit; i++)
	{
	    iwordddl++;

	    blHeaderWord[i-ilowLimit] = (Int_t) GetNextWord();
	}

      blockHeader.SetHeader(blHeaderWord);
      blRawDataLength = blockHeader.GetRawDataLength();

      if (iwordddl == totaldataword) continue;

      Int_t iwordblk = 0;

      for (Int_t idsp = 0; idsp < 5; idsp++)
	{


	  ilowLimit = iuppLimit;
	  iuppLimit = ilowLimit + kdspHLen;

	  for (Int_t i = ilowLimit; i < iuppLimit; i++)
	    {
		iwordddl++;
		iwordblk++;
		dspHeaderWord[i-ilowLimit] = (Int_t) GetNextWord();
	    }
	  dspHeader.SetHeader(dspHeaderWord);
	  dspRawDataLength = dspHeader.GetRawDataLength();

	  if (iwordddl == totaldataword) continue;

	  Int_t iworddsp = 0;

	  for (ibus = 0; ibus < 5; ibus++)
	    {
	      ilowLimit = iuppLimit;
	      iuppLimit = ilowLimit + kpbusHLen;

	      for (Int_t i = ilowLimit; i < iuppLimit; i++)
		{
		    iwordddl++;
		    iwordblk++;
		    iworddsp++;
		  pbusHeaderWord[i-ilowLimit] = (Int_t) GetNextWord();
		}

	      pbusHeader.SetHeader(pbusHeaderWord);
	      Int_t rawdatalength = pbusHeader.GetRawDataLength();
	      Int_t pbusid = pbusHeader.GetPatchBusId();

	      if (pbusid < 0 || pbusid > 50) return -1;

	      ilowLimit = iuppLimit;
	      iuppLimit = ilowLimit + rawdatalength;

	      Int_t imodule = moduleNo[pbusid];

	      if (iwordddl == totaldataword) continue;

	      for (Int_t iword = ilowLimit; iword < iuppLimit; iword++)
		{
		    iwordddl++;
		    iwordblk++;
		    iworddsp++;
		    data = 0;
		    data = GetNextWord();

		  Int_t isig =  data & 0x0FFF;
		  Int_t ich  = (data >> 12) & 0x003F;
		  Int_t imcm = (data >> 18) & 0x07FF;
		  Int_t ibit = (data >> 31) & 0x0001;
		  parity = ComputeParity(data);
		  if (ibit != parity)
		    {
		      AliWarning("ComputeParity:: Parity Error");
		      fRawReader->AddMajorErrorLog(kParityError);
		      return -1;
		    }

		  ConvertDDL2SMN(iddl, imodule, ismn, idet);

		  GetRowCol(iddl, ismn, pbusid, imcm, ich, 
			    startRowBus, endRowBus,
			    startColBus, endColBus,
			    irow, icol);

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

	      if (iwordddl == totaldataword) break;

	      if (iworddsp == dspRawDataLength) break; // raw data

	    } // patch bus loop

	  if (dspHeader.GetPaddingWord() == 1)
	  {
	      iuppLimit++;
	      iwordddl++;
	      iwordblk++;
	      iworddsp++;
	      data = GetNextWord();
	  }

	  if (iwordblk == blRawDataLength) break; // for raw data

	} // end of DSP



    } // end of BLOCK
  
  return iddl;
}
//_____________________________________________________________________________
void AliPMDRawStream::GetRowCol(Int_t ddlno, Int_t smn, Int_t pbusid,
				UInt_t mcmno, UInt_t chno,
				Int_t startRowBus[], Int_t endRowBus[],
				Int_t startColBus[], Int_t endColBus[],
				Int_t &row, Int_t &col) const
{
// decode: ddlno, patchbusid, mcmno, chno -> um, row, col

  UInt_t iCh[64];

  static const UInt_t kChDdl01[64] = { 9, 6, 5, 10, 1, 2, 0, 3,
				       13, 7, 4, 11, 8, 14, 12, 15,
				       16, 19, 17, 23, 20, 27, 24, 18,
				       28, 31, 29, 30, 21, 26, 25, 22,
				       41, 38, 37, 42, 33, 34, 32, 35,
				       45, 39, 36, 43, 40, 46, 44, 47,
				       48, 51, 49, 55, 52, 59, 56, 50,
				       60, 63, 61, 62, 53, 58, 57, 54 };
    
  static const UInt_t kChDdl23[64] = { 54, 57, 58, 53, 62, 61, 63, 60,
				       50, 56, 59, 52, 55, 49, 51, 48,
				       47, 44, 46, 40, 43, 36, 39, 45,
				       35, 32, 34, 33, 42, 37, 38, 41,
				       22, 25, 26, 21, 30, 29, 31, 28,
				       18, 24, 27, 20, 23, 17, 19, 16,
				       15, 12, 14, 8, 11, 4, 7, 13,
				       3, 0, 2, 1, 10, 5, 6, 9 };
  
  static const UInt_t kChDdl41[64] = { 53, 58, 57, 54, 61, 62, 60, 63,
				       49, 59, 56, 55, 52, 50, 48, 51,
				       44, 47, 45, 43, 40, 39, 36, 46,
				       32, 35, 33, 34, 41, 38, 37, 42,
				       21, 26, 25, 22, 29, 30, 28, 31,
				       17, 27, 24, 23, 20, 18, 16, 19,
				       12, 15, 13, 11, 8, 7, 4, 14,
				       0, 3, 1, 2, 9, 6, 5, 10 };
  
  static const UInt_t kChDdl42[64] = { 10, 5, 6, 9, 2, 1, 3, 0,
				       14, 4, 7, 8, 11, 13, 15, 12,
				       19, 16, 18, 20, 23, 24, 27, 17,
				       31, 28, 30, 29, 22, 25, 26, 21,
				       42, 37, 38, 41, 34, 33, 35, 32,
				       46, 36, 39, 40, 43, 45, 47, 44,
				       51, 48, 50, 52, 55, 56, 59, 49,
				       63, 60, 62, 61, 54, 57, 58, 53 };
  
  static const UInt_t kChDdl51[64] = { 10, 5, 6, 9, 2, 1, 3, 0,
				       14, 4, 7, 8, 11, 13, 15, 12,
				       19, 16, 18, 20, 23, 24, 27, 17,
				       31, 28, 30, 29, 22, 25, 26, 21,
				       42, 37, 38, 41, 34, 33, 35, 32,
				       46, 36, 39, 40, 43, 45, 47, 44,
				       51, 48, 50, 52, 55, 56, 59, 49,
				       63, 60, 62, 61, 54, 57, 58, 53 };
  
  static const UInt_t kChDdl52[64] = { 53, 58, 57, 54, 61, 62, 60, 63,
				       49, 59, 56, 55, 52, 50, 48, 51,
				       44, 47, 45, 43, 40, 39, 36, 46,
				       32, 35, 33, 34, 41, 38, 37, 42,
				       21, 26, 25, 22, 29, 30, 28, 31,
				       17, 27, 24, 23, 20, 18, 16, 19,
				       12, 15, 13, 11, 8, 7, 4, 14,
				       0, 3, 1, 2, 9, 6, 5, 10 };
  
  for (Int_t i = 0; i < 64; i++)
  {
      if (ddlno == 0 || ddlno == 1) iCh[i] = kChDdl01[i];
      if (ddlno == 2 || ddlno == 3) iCh[i] = kChDdl23[i];
      
      if (ddlno == 4 && smn < 6)                iCh[i] = kChDdl41[i];
      if (ddlno == 4 && (smn >= 18 && smn < 24))iCh[i] = kChDdl42[i];
      if (ddlno == 5 && (smn >= 12 && smn < 18))iCh[i] = kChDdl51[i];
      if (ddlno == 5 && (smn >=  6 && smn < 12))iCh[i] = kChDdl52[i];
  }
  
  
  Int_t rowcol  = iCh[chno];
  Int_t irownew = rowcol/4;
  Int_t icolnew = rowcol%4;

  if (ddlno == 0 )
    {
      row = startRowBus[pbusid] + irownew;
      col = startColBus[pbusid] + (mcmno-1)*4 + icolnew;
    }
  else if (ddlno == 1)
    {
    row = endRowBus[pbusid] - (15 - irownew);
    col = startColBus[pbusid] + (mcmno-1)*4 + icolnew;
    
    }
  else if (ddlno == 2 )
    {
      row = startRowBus[pbusid] + irownew;
      col = endColBus[pbusid] - (mcmno-1)*4 - (3 - icolnew);
    }
  else if (ddlno == 3)
    {
    row = endRowBus[pbusid] - (15 - irownew);
    col = endColBus[pbusid] - (mcmno-1)*4 - (3 - icolnew);
    }
  else if (ddlno == 4 )
    {
      if (pbusid  < 19)
	{
	  if (mcmno <= 12)
	    {
	      // Add 16 to skip the 1st 15 rows
	      row = startRowBus[pbusid] + irownew + 16;
	      col = startColBus[pbusid] + (mcmno-1)*4 + icolnew;
	    }
	  else if(mcmno > 12)
	    {
	      row = startRowBus[pbusid] + irownew;
	      col = startColBus[pbusid] + (mcmno-12-1)*4 + icolnew;
	    }
	}
      else if(pbusid > 18)
	{
	  if (mcmno <= 12)
	    {
	      col = endColBus[pbusid] - (mcmno-1)*4 - (3 - icolnew); 

	      if(endRowBus[pbusid] - startRowBus[pbusid] > 16)
		row = endRowBus[pbusid] - (15 - irownew) - 16 ;
	      else
		row = endRowBus[pbusid] - (15 - irownew) ;
	    }
	  else if(mcmno > 12)
	    {
	      row = endRowBus[pbusid] - (15 - irownew)  ;
	      col = endColBus[pbusid] - (mcmno - 12 - 1)*4 - (3 - icolnew);
	    }
	}
    }
  else if (ddlno == 5)
    {
      if (pbusid  <= 18)
	{
	  if (mcmno > 12)
	    {
	      // Subtract 16 to skip the 1st 15 rows
	      row = endRowBus[pbusid] - 16 -(15 - irownew);
	      col = startColBus[pbusid] + (mcmno-12 -1)*4 + icolnew;
	    }
	  else
	    {
	      row = endRowBus[pbusid]  - (15 - irownew) ;
	      col = startColBus[pbusid] + (mcmno -1)*4 + icolnew;
	    }

	}
      
      else if (pbusid > 18)
	{
	  if(mcmno > 12)
	    {
	      // Add 16 to skip the 1st 15 rows
	      row = startRowBus[pbusid] + irownew + 16;
	      col = endColBus[pbusid] - (mcmno - 12 - 1)*4 - (3 - icolnew);
	    }
	  else 
	    {
	      row = startRowBus[pbusid] + irownew ;
	      col = endColBus[pbusid] - (mcmno - 1)*4 - (3 - icolnew); 
	    }
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
Int_t AliPMDRawStream::ComputeParity(UInt_t data)
{
// Calculate the parity bit

  Int_t count = 0;
  for(Int_t j = 0; j<29; j++)
    {
      if (data & 0x01 ) count++;
      data >>= 1;
    }
  
  Int_t parity = count%2;

  return parity;
}

//_____________________________________________________________________________
UInt_t AliPMDRawStream::GetNextWord()
{
    // Returns the next 32 bit word
    // inside the raw data payload.

    if (!fData || fPosition < 0) AliFatal("Raw data payload buffer is not yet initialized !");

    UInt_t word = 0;
    word |= fData[fPosition++];
    word |= fData[fPosition++] << 8;
    word |= fData[fPosition++] << 16;
    word |= fData[fPosition++] << 24;

    return word;
}

//_____________________________________________________________________________
void AliPMDRawStream::Ddl0Mapping(Int_t moduleNo[],    Int_t mcmperBus[],
				  Int_t startRowBus[], Int_t endRowBus[],
				  Int_t startColBus[], Int_t endColBus[])
{
// DDL0 Mapping 

    const Int_t ktotpbus = 36;
    static const Int_t modno0[ktotpbus] = {0, 0, 0, 0, 0, 0,
					   1, 1, 1, 1, 1, 1,
					   2, 2, 2, 2, 2, 2,
					   3, 3, 3, 3, 3, 3,
					   4, 4, 4, 4, 4, 4,
					   5, 5, 5, 5, 5, 5,};



    static const Int_t srbus0[ktotpbus] = {0, 16, 32, 48, 64, 80,
					   0, 16, 32, 48, 64, 80,
					   0, 16, 32, 48, 64, 80,
					   0, 16, 32, 48, 64, 80,
					   0, 16, 32, 48, 64, 80,
					   0, 16, 32, 48, 64, 80};


	
    static const Int_t erbus0[ktotpbus] = {15, 31, 47, 63, 79, 95,
					   15, 31, 47, 63, 79, 95,
					   15, 31, 47, 63, 79, 95,
					   15, 31, 47, 63, 79, 95,
					   15, 31, 47, 63, 79, 95,
					   15, 31, 47, 63, 79, 95};
    
	
    for (Int_t ibus = 1; ibus <= ktotpbus; ibus++)
    {
	moduleNo[ibus]      = modno0[ibus-1];
	mcmperBus[ibus]     = 12;
	startRowBus[ibus]   = srbus0[ibus-1];
	endRowBus[ibus]     = erbus0[ibus-1];
	startColBus[ibus]   = 0;
	endColBus[ibus]     = 47;
    }
	
}

//_____________________________________________________________________________
void AliPMDRawStream::Ddl1Mapping(Int_t moduleNo[],    Int_t mcmperBus[],
				  Int_t startRowBus[], Int_t endRowBus[],
				  Int_t startColBus[], Int_t endColBus[])
{
// DDL1 Mapping 

    const Int_t ktotpbus = 36;
    static const Int_t kmodno1[ktotpbus] = {6, 6, 6, 6, 6, 6,
					   7, 7, 7, 7, 7, 7,
					   8, 8, 8, 8, 8, 8,
					   9, 9, 9, 9, 9, 9,
					   10, 10, 10, 10, 10, 10,
					   11, 11, 11, 11, 11, 11};


    static const Int_t ksrbus1[ktotpbus] = {0, 16, 32, 48, 64, 80,
					   0, 16, 32, 48, 64, 80,
					   0, 16, 32, 48, 64, 80,
					   0, 16, 32, 48, 64, 80,
					   0, 16, 32, 48, 64, 80,
					   0, 16, 32, 48, 64, 80};


	
    static const Int_t kerbus1[ktotpbus] = {15, 31, 47, 63, 79, 95,
					   15, 31, 47, 63, 79, 95,
					   15, 31, 47, 63, 79, 95,
					   15, 31, 47, 63, 79, 95,
					   15, 31, 47, 63, 79, 95,
					   15, 31, 47, 63, 79, 95};
    
	
    for (Int_t ibus = 1; ibus <= ktotpbus; ibus++)
    {
	moduleNo[ibus]      = kmodno1[ibus-1];
	mcmperBus[ibus]     = 12;
	startRowBus[ibus]   = ksrbus1[ibus-1];
	endRowBus[ibus]     = kerbus1[ibus-1];
	startColBus[ibus]   = 0;
	endColBus[ibus]     = 47;
    }
	
}

//_____________________________________________________________________________
void AliPMDRawStream::Ddl2Mapping(Int_t moduleNo[],    Int_t mcmperBus[],
				  Int_t startRowBus[], Int_t endRowBus[],
				  Int_t startColBus[], Int_t endColBus[])
{
// DDL2 Mapping 

    const Int_t ktotpbus = 36;
    static const Int_t kmodno2[ktotpbus] = {12, 12, 12, 12, 12, 12,
					    13, 13, 13, 13, 13, 13,
					    14, 14, 14, 14, 14, 14,
					    15, 15, 15, 15, 15, 15,
					    16, 16, 16, 16, 16, 16,
					    17, 17, 17, 17, 17, 17};

    static const Int_t ksrbus2[ktotpbus] = {32, 32, 16, 16, 0, 0,
					    32, 32, 16, 16, 0, 0,
					    32, 32, 16, 16, 0, 0,
					    32, 32, 16, 16, 0, 0,
					    32, 32, 16, 16, 0, 0,
					    32, 32, 16, 16, 0, 0};

	
    static const Int_t kerbus2[ktotpbus] = {47, 47, 31, 31, 15, 15,
					    47, 47, 31, 31, 15, 15,
					    47, 47, 31, 31, 15, 15,
					    47, 47, 31, 31, 15, 15,
					    47, 47, 31, 31, 15, 15,
					    47, 47, 31, 31, 15, 15};

    static const Int_t kscbus2[ktotpbus] = {48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0};

    static const Int_t kecbus2[ktotpbus] = {95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47};


	
    for (Int_t ibus = 1; ibus <= ktotpbus; ibus++)
    {
	moduleNo[ibus]      = kmodno2[ibus-1];
	mcmperBus[ibus]     = 12;
	startRowBus[ibus]   = ksrbus2[ibus-1];
	endRowBus[ibus]     = kerbus2[ibus-1];
	startColBus[ibus]   = kscbus2[ibus-1];
	endColBus[ibus]     = kecbus2[ibus-1];
    }
	
}

//_____________________________________________________________________________
void AliPMDRawStream::Ddl3Mapping(Int_t moduleNo[],    Int_t mcmperBus[],
				  Int_t startRowBus[], Int_t endRowBus[],
				  Int_t startColBus[], Int_t endColBus[])
{
// DDL3 Mapping 

    const Int_t ktotpbus = 36;
    static const Int_t kmodno3[ktotpbus] = {18, 18, 18, 18, 18, 18,
					    19, 19, 19, 19, 19, 19,
					    20, 20, 20, 20, 20, 20,
					    21, 21, 21, 21, 21, 21,
					    22, 22, 22, 22, 22, 22,
					    23, 23, 23, 23, 23, 23};
    


    static const Int_t ksrbus3[ktotpbus] = {32, 32, 16, 16, 0, 0,
					    32, 32, 16, 16, 0, 0,
					    32, 32, 16, 16, 0, 0,
					    32, 32, 16, 16, 0, 0,
					    32, 32, 16, 16, 0, 0,
					    32, 32, 16, 16, 0, 0};
    
	
    static const Int_t kerbus3[ktotpbus] = {47, 47, 31, 31, 15, 15,
					    47, 47, 31, 31, 15, 15,
					    47, 47, 31, 31, 15, 15,
					    47, 47, 31, 31, 15, 15,
					    47, 47, 31, 31, 15, 15,
					    47, 47, 31, 31, 15, 15};
    
    static const Int_t kscbus3[ktotpbus] = {48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0};
	
    static const Int_t kecbus3[ktotpbus] = {95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47};
    
    for (Int_t ibus = 1; ibus <= ktotpbus; ibus++)
    {
	moduleNo[ibus]      = kmodno3[ibus-1];
	mcmperBus[ibus]     = 12;
	startRowBus[ibus]   = ksrbus3[ibus-1];
	endRowBus[ibus]     = kerbus3[ibus-1];
	startColBus[ibus]   = kscbus3[ibus-1];
	endColBus[ibus]     = kecbus3[ibus-1];
    }
	
}

//_____________________________________________________________________________
void AliPMDRawStream::Ddl4Mapping(Int_t moduleNo[],    Int_t mcmperBus[],
				  Int_t startRowBus[], Int_t endRowBus[],
				  Int_t startColBus[], Int_t endColBus[])
{
// DDL4 Mapping 

    const Int_t ktotpbus = 42;
    static const Int_t kmodno4[ktotpbus] = {24, 24, 24, 25, 25, 25,
					    26, 26, 26, 27, 27, 27,
					    28, 28, 28, 29, 29, 29,
					    42, 42, 42, 42, 43, 43, 43, 43,
					    44, 44, 44, 44, 45, 45, 45, 45,
					    46, 46, 46, 46, 47, 47, 47, 47};
    

    static const Int_t kmcbus4[ktotpbus] = {24, 24, 24, 24, 24, 24,
					    24, 24, 24, 24, 24, 24,
					    24, 24, 24, 24, 24, 24,
					    12, 12, 24, 24, 12, 12, 24, 24,
					    12, 12, 24, 24, 12, 12, 24, 24,
					    12, 12, 24, 24, 12, 12, 24, 24};


    static const Int_t ksrbus4[ktotpbus] = {0, 32, 64, 0, 32, 64,
					    0, 32, 64, 0, 32, 64,
					    0, 32, 64, 0, 32, 64,
					    32, 32, 0, 0, 32, 32, 0, 0,
					    32, 32, 0, 0, 32, 32, 0, 0,
					    32, 32, 0, 0, 32, 32, 0, 0};



	
    static const Int_t kerbus4[ktotpbus] = {31, 63, 95, 31, 63, 95,
					    31, 63, 95, 31, 63, 95,
					    31, 63, 95, 31, 63, 95,
					    47, 47, 31, 31, 47, 47, 31, 31,
					    47, 47, 31, 31, 47, 47, 31, 31,
					    47, 47, 31, 31, 47, 47, 31, 31};


	
    static const Int_t kscbus4[ktotpbus] = {0, 0, 0, 0, 0, 0,
					    0, 0, 0, 0, 0, 0,
					    0, 0, 0, 0, 0, 0,
					    48, 0, 48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0, 48, 0};


	
    static const Int_t kecbus4[ktotpbus] = {47, 47, 47, 47, 47, 47,
					    47, 47, 47, 47, 47, 47,
					    47, 47, 47, 47, 47, 47,
					    95, 47, 95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47, 95, 47};

	
    for (Int_t ibus = 1; ibus <= ktotpbus; ibus++)
    {
	moduleNo[ibus]      = kmodno4[ibus-1];
	mcmperBus[ibus]     = kmcbus4[ibus-1];
	startRowBus[ibus]   = ksrbus4[ibus-1];
	endRowBus[ibus]     = kerbus4[ibus-1];
	startColBus[ibus]   = kscbus4[ibus-1];
	endColBus[ibus]     = kecbus4[ibus-1];
    }
	
}

//_____________________________________________________________________________
void AliPMDRawStream::Ddl5Mapping(Int_t moduleNo[],    Int_t mcmperBus[],
				  Int_t startRowBus[], Int_t endRowBus[],
				  Int_t startColBus[], Int_t endColBus[])
{
// DDL5 Mapping 

    const Int_t ktotpbus = 42;
    static const Int_t kmodno5[ktotpbus] = {30, 30, 30, 31, 31, 31,
					    32, 32, 32, 33, 33, 33,
					    34, 34, 34, 35, 35, 35,
					    36, 36, 36, 36, 37, 37, 37, 37,
					    38, 38, 38, 38, 39, 39, 39, 39,
					    40, 40, 40, 40, 41, 41, 41, 41};


    static const Int_t kmcbus5[ktotpbus] = {24, 24, 24, 24, 24, 24,
					    24, 24, 24, 24, 24, 24,
					    24, 24, 24, 24, 24, 24,
					    12, 12, 24, 24, 12, 12, 24, 24,
					    12, 12, 24, 24, 12, 12, 24, 24,
					    12, 12, 24, 24, 12, 12, 24, 24};
					    


    static const Int_t ksrbus5[ktotpbus] = {0, 32, 64, 0, 32, 64,
					    0, 32, 64, 0, 32, 64,
					    0, 32, 64, 0, 32, 64,
					    32, 32, 0, 0, 32, 32, 0, 0,
					    32, 32, 0, 0, 32, 32, 0, 0,
					    32, 32, 0, 0, 32, 32, 0, 0};



	
    static const Int_t kerbus5[ktotpbus] = {31, 63, 95, 31, 63, 95,
					    31, 63, 95, 31, 63, 95,
					    31, 63, 95, 31, 63, 95,
					    47, 47, 31, 31, 47, 47, 31, 31,
					    47, 47, 31, 31, 47, 47, 31, 31,
					    47, 47, 31, 31, 47, 47, 31, 31};


	
    static const Int_t kscbus5[ktotpbus] = {0, 0, 0, 0, 0, 0,
					    0, 0, 0, 0, 0, 0,
					    0, 0, 0, 0, 0, 0,
					    48, 0, 48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0, 48, 0,
					    48, 0, 48, 0, 48, 0, 48, 0};


	
    static const Int_t kecbus5[ktotpbus] = {47, 47, 47, 47, 47, 47,
					    47, 47, 47, 47, 47, 47,
					    47, 47, 47, 47, 47, 47,
					    95, 47, 95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47, 95, 47,
					    95, 47, 95, 47, 95, 47, 95, 47};

	
    for (Int_t ibus = 1; ibus <= ktotpbus; ibus++)
    {
	moduleNo[ibus]      = kmodno5[ibus-1];
	mcmperBus[ibus]     = kmcbus5[ibus-1];
	startRowBus[ibus]   = ksrbus5[ibus-1];
	endRowBus[ibus]     = kerbus5[ibus-1];
	startColBus[ibus]   = kscbus5[ibus-1];
	endColBus[ibus]     = kecbus5[ibus-1];
    }
	
}


