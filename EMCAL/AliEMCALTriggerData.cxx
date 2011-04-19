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
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
 
 
EMCal trigger data container
for data (both raw & rec) persistency
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerData.h"
#include "AliEMCALTriggerPatch.h"
#include "AliLog.h"
#include "TIterator.h"
#include "Riostream.h"

ClassImp(AliEMCALTriggerData)

//_____________
AliEMCALTriggerData::AliEMCALTriggerData() : TObject(),
fMode(0),
fL0Patches(),
fL0Region(),
fL1GammaPatches(),
fL1JetPatches(),
fL1Region(),
fL1GammaThreshold(0),
fL1JetThreshold(0),
fL1V0(),
fL1FrameMask(0),
fL1TriggerType(),
fL1DataDecoded(0)
{  
	//
	for (Int_t i = 0; i < 2; i++)
	{
		       fL0Patches[i] = new TClonesArray("AliEMCALTriggerPatch");
		  fL1GammaPatches[i] = new TClonesArray("AliEMCALTriggerPatch");
		    fL1JetPatches[i] = new TClonesArray("AliEMCALTriggerPatch");
	}
	
	for (Int_t i = 0; i < 32; i++) for (Int_t j = 0; j < 24; j++) for (Int_t k = 0; k <  4; k++) fL0Region[i][j][k] = 0;
	for (Int_t i = 0; i <  2; i++) for (Int_t j = 0; j < 48; j++) for (Int_t k = 0; k < 64; k++) fL1Region[i][j][k] = 0;
	
	fL1V0[0] = fL1V0[1] = 0;
	for (Int_t i = 0; i < 8; i++) fL1TriggerType[i] = 0;	
}

//_____________
AliEMCALTriggerData::~AliEMCALTriggerData()
{
	//
	for (Int_t i = 0; i < 2; i++)
	{
		if (     fL0Patches[i])      fL0Patches[i]->Delete();
		if (fL1GammaPatches[i]) fL1GammaPatches[i]->Delete();
		if (  fL1JetPatches[i])   fL1JetPatches[i]->Delete();
	}
}

//_____________
void AliEMCALTriggerData::SetL0Region(Int_t i, const Int_t**& region)
{
	//
	if (i < 0 || i > 31) 
	{
		AliError("Bad index!");
		return;
	}
	
  	for (Int_t j=0;j<24;j++)
  		for (Int_t k=0;k<4;k++) fL0Region[i][j][k] = region[j][k];
}

//_____________
void AliEMCALTriggerData::GetPatches(TriggerType_t type, Int_t i, TClonesArray& patches) const
{
	//
	if (i < 0 || i > 1) 
	{
		AliError("Bad index!");
		return;
	}
	
	switch (type)
	{
		case kL0:
			patches =  *fL0Patches[i];
			break;
		case kL1Gamma:
			patches =  *fL1GammaPatches[i];
			break;
		case kL1Jet:
			patches =  *fL1JetPatches[i];
			break;
		default:
			AliError("Unknown trigger type!");
			break;
	}
}

//_____________
TClonesArray* AliEMCALTriggerData::GetPatches(TriggerType_t type, Int_t i) const
{
	//
	if (i < 0 || i > 1) 
	{
		AliError("Bad index!");
		return 0x0;
	}
	
	switch (type)
	{
		case kL0:
			return fL0Patches[i];
			break;
		case kL1Gamma:
			return fL1GammaPatches[i];
			break;
		case kL1Jet:
			return fL1JetPatches[i];
			break;
		default:
			AliError("Unknown trigger type!");
			break;
	}

	return 0x0;
}

//_____________
void AliEMCALTriggerData::SetPatches(TriggerType_t type, Int_t i, const TClonesArray& patches)
{
	//
	if (i < 0 || i > 1) 
	{
		AliError("Bad index!");
		return;
	}
	
	if (patches.GetEntriesFast())
	{
		TClonesArray* arr = 0x0;
		
		switch (type)
		{
			case kL0:
				arr = fL0Patches[i];
				break;
			case kL1Gamma:
				arr = fL1GammaPatches[i];
				break;
			case kL1Jet:
				arr = fL1JetPatches[i];
				break;
			default:
				AliError("Unknown trigger type!");
				return;
		}
		
		if (arr)
		{
			Int_t size = arr->GetSize() + patches.GetSize();
		
			arr->Expand(size);
		
			for (Int_t k = 0; k < patches.GetEntriesFast(); k++)
			{
				AliEMCALTriggerPatch* p = static_cast<AliEMCALTriggerPatch*>(patches.At(k));
				new((*arr)[arr->GetEntriesFast()]) AliEMCALTriggerPatch(*p);
			}
		}
		else
		{
			AliError("TClonesArray is NULL!");
		}
	}
}

//_____________
void AliEMCALTriggerData::SetL1Region(Int_t i, Int_t**& region)
{
	//
	if (i < 0 || i > 1) 
	{
		AliError("Bad index!");
		return;
	}
		
  	for (Int_t j = 0; j < 48; j++)
  		for (Int_t k = 0; k < 64; k++) fL1Region[i][j][k] = region[j][k];
}

//_____________
void AliEMCALTriggerData::GetL1Region(Int_t i, Int_t arr[][64]) const 
{ 
	//
	if (i < 0 || i > 1) 
	{
		AliError("Bad index!");
		return;
	}
	
	for (Int_t j = 0; j < 48; j++) for (Int_t k = 0; k < 64; k++) { arr[j][k] = fL1Region[i][j][k]; } 
}


//_____________
void AliEMCALTriggerData::Scan() const
{
	//
	TIterator* nP;

	printf("L0:\n");
	printf("\tFound (%2d,%2d) patches\n", fL0Patches[1]->GetEntriesFast(), fL0Patches[0]->GetEntriesFast());
	printf("\tRAW:\n");
	nP = fL0Patches[1]->MakeIterator();
	while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) {printf("\t"); p->Print("");}
	printf("\tREC:\n");
	nP = fL0Patches[0]->MakeIterator();
	while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) {printf("\t"); p->Print("");}
	printf("L1:\n");
	printf("\tFound (%4d,%4d) gamma patches\n",fL1GammaPatches[1]->GetEntriesFast(), fL1GammaPatches[0]->GetEntriesFast());
	printf("\tRAW:\n");
	nP = fL1GammaPatches[1]->MakeIterator();
	while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) {printf("\t"); p->Print("");}
	printf("\tREC:\n");
	nP = fL1GammaPatches[0]->MakeIterator();
	while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) {printf("\t"); p->Print("");}
	printf("\tFound (%4d,%4d) jet patches\n",fL1JetPatches[1]->GetEntriesFast(), fL1JetPatches[0]->GetEntriesFast());
	printf("\tRAW:\n");
	nP = fL1JetPatches[1]->MakeIterator();
	while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) {printf("\t"); p->Print("");}
	printf("\tREC:\n");
	nP = fL1JetPatches[0]->MakeIterator();
	while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) {printf("\t"); p->Print("");}
}

//_____________
void AliEMCALTriggerData::Reset()
{
	//
	for (Int_t i = 0; i < 2; i++)
	{
		if (     fL0Patches[i])      fL0Patches[i]->Delete();
		if (fL1GammaPatches[i]) fL1GammaPatches[i]->Delete();
		if (  fL1JetPatches[i])   fL1JetPatches[i]->Delete();	
	}
	 	
	for (Int_t i = 0; i < 2; i++) for (Int_t j = 0; j < 48; j++) for (Int_t k = 0; k < 64; k++) fL1Region[i][j][k] = 0;
	
	fL1DataDecoded = 0;
}



