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

//-----------------------------------------------------//
//                                                     //
//           Date   : October 09 2006                  //
//       converts raw to sdigit and digit              //
//                                                     //
//-----------------------------------------------------//

#include <Riostream.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliDAQ.h"
#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliRawReader.h"

#include "AliPMDsdigit.h"
#include "AliPMDdigit.h"
#include "AliPMDRawToSDigits.h"
#include "AliPMDRawStream.h"
#include "AliPMDddldata.h"


ClassImp(AliPMDRawToSDigits)

AliPMDRawToSDigits::AliPMDRawToSDigits():
  fSDigits(new TClonesArray("AliPMDsdigit", 1000)),
  fDigits(new TClonesArray("AliPMDdigit", 1000)),
  fNsdigit(0),
  fNdigit(0)
{
  //
// Constructor
//

}
// ------------------------------------------------------------------------- //
AliPMDRawToSDigits::AliPMDRawToSDigits(const AliPMDRawToSDigits & /*pmdr2sd*/):
  TObject(/* pmdr2sd */),
  fSDigits(NULL),
  fDigits(NULL),
  fNsdigit(0),
  fNdigit(0)
{
  // copy constructor
  AliFatal("Copy constructor not implemented");
}
// ------------------------------------------------------------------------- //


AliPMDRawToSDigits &AliPMDRawToSDigits::operator=(const AliPMDRawToSDigits &/* pmdr2sd */)
{
  // assignment operator
  AliFatal("Assignment operator not implemented");
  return *this;
}

// ------------------------------------------------------------------------- //

AliPMDRawToSDigits::~AliPMDRawToSDigits()
{
  // Destructor
  if (fSDigits)
    {
      fSDigits->Delete();
      delete fSDigits;
      fSDigits=0;
    }
  if (fDigits)
    {
      fDigits->Delete();
      delete fDigits;
      fDigits=0;
    }

}
// ------------------------------------------------------------------------- //

void AliPMDRawToSDigits::Raw2SDigits(AliRunLoader *runLoader, AliRawReader *rawReader)
{
  // Converts RAW data to sdigits
  //
  TObjArray pmdddlcont;
  AliLoader *pmdLoader = runLoader->GetLoader("PMDLoader");
  
  TTree* treeS = pmdLoader->TreeS();
  if (treeS == 0x0)
    {
      pmdLoader->MakeTree("S");
      treeS = pmdLoader->TreeS();
    }
  Int_t bufsize = 16000;
  if (!fSDigits) fSDigits = new TClonesArray("AliPMDsdigit", 1000);
  treeS->Branch("PMDSDigit", &fSDigits, bufsize); 

//  const Int_t kDDL = AliDAQ::NumberOfDdls("PMD");
  const Int_t kRow = 48;
  const Int_t kCol = 96;
  const Int_t kSMN = 48;

  Int_t idet = 0;
  Int_t indexsmn = 0;
  Int_t ismn = 0;
  
  Int_t ***precpvADC;
  precpvADC = new int **[kSMN];
  for (Int_t i=0; i<kSMN; i++) precpvADC[i] = new int *[kRow];
  for (Int_t i=0; i<kSMN;i++)
    {
      for (Int_t j=0; j<kRow; j++) precpvADC[i][j] = new int [kCol];
    }
  for (Int_t i = 0; i < kSMN; i++)
    {
      for (Int_t j = 0; j < kRow; j++)
	{
	  for (Int_t k = 0; k < kCol; k++)
	    {
	      precpvADC[i][j][k] = 0;
	    }
	}
    }
  
  AliPMDRawStream pmdinput(rawReader);
  Int_t indexDDL = -1;
  while ((indexDDL = pmdinput.DdlData(&pmdddlcont)) >=0)
  {
      Int_t ientries = pmdddlcont.GetEntries();
      for (Int_t ient = 0; ient < ientries; ient++)
	{
	  AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	  
	  Int_t det = pmdddl->GetDetector();
	  Int_t smn = pmdddl->GetSMN();
	  //Int_t mcm = pmdddl->GetMCM();
	  //Int_t chno = pmdddl->GetChannel();
	  Int_t row = pmdddl->GetRow();
	  Int_t col = pmdddl->GetColumn();
	  Int_t sig = pmdddl->GetSignal();
	  
	  if (indexDDL < 4)
	    {
	      if (det != 0)
		AliError(Form("*DDL %d and Detector NUMBER %d NOT MATCHING *",
			      indexDDL, det));
	      indexsmn = smn;
	    }
	  else if (indexDDL == 4 || indexDDL == 5)
	    {
	      if (det != 1)
		AliError(Form("*DDL %d and Detector NUMBER %d NOT MATCHING *",
			      indexDDL, det));
	      indexsmn = smn + 24;
	    }

	  precpvADC[indexsmn][row][col] = sig;
	}
      
      pmdddlcont.Clear();
      
    } // DDL Loop
      
  // Add the sdigits here
  
  for ( indexsmn = 0; indexsmn < kSMN; indexsmn++)
    {

      if (indexsmn < 23)
	{
	  idet = 0;
	  ismn = indexsmn;
	}
      else if (indexsmn > 23)
	{
	  idet = 0;
	  ismn = indexsmn - 24;
	}
      for (Int_t irow = 0; irow < kRow; irow++)
	{
	  for (Int_t icol = 0; icol < kCol; icol++)
	    {

	      Int_t trno  = -99999;   // when extracted from raw data
	      Int_t trpid = -99999;   // when extracted from raw data
	      Int_t sig1  = precpvADC[indexsmn][irow][icol];
	      
	      // plug in a function to convert to adc to MeV
	      Float_t edep = 0.;
	      if (sig1 > 0)
		{
		  AdcToMeV(sig1,edep);
		  AddSDigit(trno,trpid,idet,ismn,irow,icol,edep);
		}
	    } // row
	}     // col
      
      treeS->Fill();
      ResetSDigit();
    }

  pmdLoader->WriteSDigits("OVERWRITE");

  // Delete all the pointers
  
  for (Int_t i = 0; i < kSMN; i++)
    {
      for (Int_t j = 0; j < kRow; j++)
	{
	  delete [] precpvADC[i][j];
	}
    }
  for (Int_t j = 0; j < kSMN; j++)
    {
      delete [] precpvADC[j];
    }
  delete [] precpvADC;
  
}
// ------------------------------------------------------------------------- //
void AliPMDRawToSDigits::Raw2Digits(AliRunLoader *runLoader, AliRawReader *rawReader)
{
  // Converts RAW data to digits
  //
  TObjArray pmdddlcont;
  
  AliLoader *pmdLoader = runLoader->GetLoader("PMDLoader");
  
  TTree* treeD = pmdLoader->TreeD();
  if (treeD == 0x0)
    {
      pmdLoader->MakeTree("D");
      treeD = pmdLoader->TreeD();
    }
  Int_t bufsize = 16000;
  if (!fDigits) fDigits = new TClonesArray("AliPMDdigit", 1000);
  treeD->Branch("PMDDigit", &fDigits, bufsize); 
  
//  const Int_t kDDL = AliDAQ::NumberOfDdls("PMD");
  const Int_t kRow = 48;
  const Int_t kCol = 96;
  const Int_t kSMN = 48;
  
  Int_t idet = 0;
  Int_t ismn = 0;
  Int_t indexsmn = 0;
  
  Int_t ***precpvADC;
  precpvADC = new int **[kSMN];
  for (Int_t i=0; i<kSMN; i++) precpvADC[i] = new int *[kRow];
  for (Int_t i=0; i<kSMN;i++)
    {
      for (Int_t j=0; j<kRow; j++) precpvADC[i][j] = new int [kCol];
    }
  for (Int_t i = 0; i < kSMN; i++)
    {
      for (Int_t j = 0; j < kRow; j++)
	{
	  for (Int_t k = 0; k < kCol; k++)
	    {
	      precpvADC[i][j][k] = 0;
	    }
	}
    }

  AliPMDRawStream pmdinput(rawReader);
  Int_t indexDDL = -1;
  while ((indexDDL = pmdinput.DdlData(&pmdddlcont)) >=0)
  {
      Int_t ientries = pmdddlcont.GetEntries();
      for (Int_t ient = 0; ient < ientries; ient++)
	{
	  AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	  
	  Int_t det = pmdddl->GetDetector();
	  Int_t smn = pmdddl->GetSMN();
	  //Int_t mcm = pmdddl->GetMCM();
	  //Int_t chno = pmdddl->GetChannel();
	  Int_t row = pmdddl->GetRow();
	  Int_t col = pmdddl->GetColumn();
	  Int_t sig = pmdddl->GetSignal();
	  
	  
	  if (indexDDL < 4)
	    {
	      if (det != 0)
		AliError(Form("*DDL %d and Detector NUMBER %d NOT MATCHING *",
			      indexDDL, det));
	      indexsmn = smn ;
	    }
	  else if (indexDDL == 4 || indexDDL == 5)
	    {
	      if (det != 1)
		AliError(Form("*DDL %d and Detector NUMBER %d NOT MATCHING *",
			      indexDDL, det));
	      indexsmn = smn + 24;
	    }	      
	  precpvADC[indexsmn][row][col] = sig;
	}
      
      pmdddlcont.Clear();
    } // DDL Loop  
      
  // Add the digits here
  for (indexsmn = 0; indexsmn < kSMN; indexsmn++)
    {
      if (indexsmn < 23)
	{
	  ismn = indexsmn;
	  idet = 0;
	}
      else if (indexsmn > 23)
	{
	  ismn = indexsmn -24;
	  idet = 1;
	}
      for (Int_t irow = 0; irow < kRow; irow++)
	{
	  for (Int_t icol = 0; icol < kCol; icol++)
	    {
	      Int_t trno  = -99999;    // when extracted from raw
	      Int_t trpid = -99999;    // when extracted from raw
	      Int_t sig1  = precpvADC[indexsmn][irow][icol];
	      
	      // plug in a function to convert to adc to MeV
	      if (sig1 > 0)
		{
		  AddDigit(trno,trpid,idet,ismn,irow,icol,sig1);
		}
	    } // row
	}     // col
      treeD->Fill();
      ResetDigit();
    }	  
  
  pmdLoader->WriteDigits("OVERWRITE");

  // Delete all the pointers

  for (Int_t i = 0; i < kSMN; i++)
    {
      for (Int_t j = 0; j < kRow; j++)
	{
	  delete [] precpvADC[i][j];
	}
    }
  for (Int_t j = 0; j < kSMN; j++)
    {
      delete [] precpvADC[j];
    }
  delete [] precpvADC;
}
// ------------------------------------------------------------------------- //

void AliPMDRawToSDigits::AdcToMeV(Int_t adc, Float_t &edep)
{
  // To be implemented, this is just for the test

  const Float_t kConstant   = 7.181;
  //  const Float_t kErConstant = 0.6899;
  const Float_t kSlope      = 35.93;
  //  const Float_t kErSlope    = 0.306;



  Float_t adc10bit = (Float_t) adc/4;
  edep     = (1000.0/kSlope)*(adc10bit - kConstant);
}

// ------------------------------------------------------------------------- //

void AliPMDRawToSDigits::AddSDigit(Int_t trnumber, Int_t trpid, Int_t det,
				   Int_t smnumber, 
				   Int_t irow, Int_t icol, Float_t adc)
{
  // Add SDigit
  //
  if (!fSDigits) fSDigits = new TClonesArray("AliPMDsdigit", 1000);
  TClonesArray &lsdigits = *fSDigits;
  new(lsdigits[fNsdigit++]) AliPMDsdigit(trnumber,trpid,det,smnumber,irow,icol,adc);
}

// ------------------------------------------------------------------------- //
void AliPMDRawToSDigits::AddDigit(Int_t trnumber, Int_t trpid, Int_t det,
				  Int_t smnumber, 
				  Int_t irow, Int_t icol, Float_t adc)
{
  // Add Digit
  //
  if (!fDigits) fDigits = new TClonesArray("AliPMDdigit", 1000);
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigit++]) AliPMDdigit(trnumber,trpid,det,smnumber,irow,icol,adc);
}

// ------------------------------------------------------------------------- //
void AliPMDRawToSDigits::ResetSDigit()
{
  // Clears SDigits
  fNsdigit = 0;
  if (fSDigits) fSDigits->Clear();
}
// ------------------------------------------------------------------------- //
void AliPMDRawToSDigits::ResetDigit()
{
  // Clears SDigits
  fNdigit = 0;
  if (fDigits) fDigits->Clear();
}
