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


/**************************************************************************************

	ACORDE's RecoParams Version 1.0

	In this version we only consider:

	->) The ACORDE's Trigger Mode (Single Muon Trigger or Multi Muon Trigger)
	->) The ACORDE's Trigger Mask (Same in SMT and MMT)
	
	In Runs PbPb, pp, and cosmics by default we have the same RecoParams.

	From: 
		Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> @ CERN
		FCFM, BUAP. Puebla, Mexico

	Further comments:

		Arturo Fernandez <afernan@mail.cern.ch>	

	March 2nd. 2009	

	NOTE: Please suggest improvements if needed.

**************************************************************************************/

#include "AliLog.h"

#include "AliACORDERecoParam.h"

ClassImp(AliACORDERecoParam)

//_____________________________________________________________________________
AliACORDERecoParam::AliACORDERecoParam():
  AliDetectorRecoParam(),
  fAcordeSingleMuonTrigger(kFALSE),
  fAcordeMultiMuonTrigger(kFALSE),
  fAcordeWord0(0x00000000),
  fAcordeWord1(0x00000000),
  fAcordeWord2(0x00000000),
  fAcordeWord3(0x00000000)       
{
 	// ACORDE's RecoParam constructor
 
 	 SetNameTitle("ACORDE","ACORDE");
}


//____________________________________________________________________________
AliACORDERecoParam::AliACORDERecoParam(const AliACORDERecoParam &p):
  AliDetectorRecoParam(p),
  fAcordeSingleMuonTrigger(p.fAcordeSingleMuonTrigger),
  fAcordeMultiMuonTrigger(p.fAcordeMultiMuonTrigger),
  fAcordeWord0(p.fAcordeWord0),
  fAcordeWord1(p.fAcordeWord1),
  fAcordeWord2(p.fAcordeWord2),
  fAcordeWord3(p.fAcordeWord3)
{
	// Copy of constructor
}

//_____________________________________________________________________________
AliACORDERecoParam& AliACORDERecoParam::operator=(const AliACORDERecoParam &p)
{
 	// ACORDE's RecoParam Assign Operator

	if (this == &p)
   	return *this;
  
  	AliDetectorRecoParam::operator=(p);
  	fAcordeSingleMuonTrigger=p.fAcordeSingleMuonTrigger;
	fAcordeMultiMuonTrigger=p.fAcordeMultiMuonTrigger;      
   	fAcordeWord0=p.fAcordeWord0;
	fAcordeWord1=p.fAcordeWord1;
	fAcordeWord2=p.fAcordeWord2;
	fAcordeWord3=p.fAcordeWord3;
   	
	return *this;
}
//_____________________________________________________________________________
AliACORDERecoParam::~AliACORDERecoParam() 
{
  	// ACORDE's RecoParam destructor
}

//_____________________________________________________________________________
AliACORDERecoParam *AliACORDERecoParam::GetPbPbparam()
{
 
  	// Set ACORDE's default reconstruction parameters for PbPb 
  
  	AliACORDERecoParam *acordeRecoParam = new AliACORDERecoParam();

	// fAcordeTriggerMode = "ACORDE_SINGLE" for ACORDE's Single Muon Trigger
 	// fAcordeTriggerMode = "ACORDE_MULTI" for ACORDE's Multi Muon Trigger
	
	// By now we set fAcordeSingeMuonTrigger as default

	acordeRecoParam->fAcordeSingleMuonTrigger=kTRUE; // Enable ACORDE_SINGLE
	acordeRecoParam->fAcordeMultiMuonTrigger=kFALSE; // Disable ACORDE_MULTI
  	acordeRecoParam->fAcordeWord0 = 0x3FFFFFFF; // [1..30] ACORDE's modules in ACORDE_SINGLE
	acordeRecoParam->fAcordeWord1 = 0x7FFFFFFF; // [31..60] ACORDE's modules in ACORDE_SINGLE
	acordeRecoParam->fAcordeWord2 = 0xBFFFFFFF; // [1..30] ACORDE's modules in ACORDE_MULTI
	acordeRecoParam->fAcordeWord3 = 0xFFFFFFFF; // [31..60] ACORDE's modules in ACORDE_MULTI
	return acordeRecoParam;
}

//_____________________________________________________________________________
AliACORDERecoParam *AliACORDERecoParam::GetPPparam()
{
	// Set ACORDE's default reconstruction parameters for PbPb 
  
        AliACORDERecoParam *acordeRecoParam = new AliACORDERecoParam();

        // fAcordeTriggerMode = "ACORDE_SINGLE" for ACORDE's Single Muon Trigger
        // fAcordeTriggerMode = "ACORDE_MULTI" for ACORDE's Multi Muon Trigger
        
        // By now we set fAcordeSingeMuonTrigger as default

        acordeRecoParam->fAcordeSingleMuonTrigger=kTRUE; // Enable ACORDE_SINGLE
        acordeRecoParam->fAcordeMultiMuonTrigger=kFALSE; // Disable ACORDE_MULTI
        acordeRecoParam->fAcordeWord0 = 0x3FFFFFFF; // [1..30] ACORDE's modules in ACORDE_SINGLE
        acordeRecoParam->fAcordeWord1 = 0x7FFFFFFF; // [31..60] ACORDE's modules in ACORDE_SINGLE
        acordeRecoParam->fAcordeWord2 = 0xBFFFFFFF; // [1..30] ACORDE's modules in ACORDE_MULTI
        acordeRecoParam->fAcordeWord3 = 0xFFFFFFFF; // [31..60] ACORDE's modules in ACORDE_MULTI
        return acordeRecoParam;
}

//_____________________________________________________________________________
AliACORDERecoParam *AliACORDERecoParam::GetCosmicMuonParam()
{
	 // Set ACORDE's default reconstruction parameters for PbPb 
  
        AliACORDERecoParam *acordeRecoParam = new AliACORDERecoParam();

        // fAcordeTriggerMode = "ACORDE_SINGLE" for ACORDE's Single Muon Trigger
        // fAcordeTriggerMode = "ACORDE_MULTI" for ACORDE's Multi Muon Trigger
        
        // By now we set fAcordeSingeMuonTrigger as default

        acordeRecoParam->fAcordeSingleMuonTrigger=kTRUE; // Enable ACORDE_SINGLE
        acordeRecoParam->fAcordeMultiMuonTrigger=kFALSE; // Disable ACORDE_MULTI
        acordeRecoParam->fAcordeWord0 = 0x3FFFFFFF; // [1..30] ACORDE's modules in ACORDE_SINGLE
        acordeRecoParam->fAcordeWord1 = 0x7FFFFFFF; // [31..60] ACORDE's modules in ACORDE_SINGLE
        acordeRecoParam->fAcordeWord2 = 0xBFFFFFFF; // [1..30] ACORDE's modules in ACORDE_MULTI
        acordeRecoParam->fAcordeWord3 = 0xFFFFFFFF; // [31..60] ACORDE's modules in ACORDE_MULTI
        return acordeRecoParam;
}

//_____________________________________________________________________________
void AliACORDERecoParam::PrintParameters() const
{
  	// Printing of the used ACORDE reconstruction parameters

 	AliInfo(Form("ACORDE's Single Muon Trigger Mode: %b", fAcordeSingleMuonTrigger));
	AliInfo(Form("ACORDE's Multi Muon Trigger Mode: %b", fAcordeMultiMuonTrigger));
	if(fAcordeSingleMuonTrigger==kTRUE)
	{
		AliInfo(Form("ACORDE's Trigger Mask: 0x%08x 0x%08x",fAcordeWord0,fAcordeWord1));
	}
	if(fAcordeMultiMuonTrigger==kTRUE)
	{
                AliInfo(Form("ACORDE's Trigger Mask: 0x%08x 0x%08x",fAcordeWord2,fAcordeWord3));

	}
}
