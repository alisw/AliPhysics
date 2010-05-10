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
for persistency of produced data presently stored in TTreeD
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerData.h"
#include "AliEMCALTriggerPatch.h"

ClassImp(AliEMCALTriggerData)

//_____________
AliEMCALTriggerData::AliEMCALTriggerData() : TObject(),
fL0Patches( new TClonesArray("AliEMCALTriggerPatch") ),
fL0NPatches(),
fL0RegionSize(0,0),
fL0SubRegionSize(0,0),
fL0PatchSize(0,0),
fL1GammaPatches( new TClonesArray("AliEMCALTriggerPatch") ),
fL1JetPatches( new TClonesArray("AliEMCALTriggerPatch") ),
fL1RegionSize(0,0),
fL1GammaPatchSize(0,0),
fL1GammaSubRegionSize(0,0),
fL1JetPatchSize(0,0),
fL1JetSubRegionSize(0,0)
{  
	for (Int_t i=0;i<32;i++) fL0NPatches[i] = 0;
	for (Int_t i=0;i<48;i++) for (Int_t j=0;j<64;j++) fL1Region[i][j] = 0;
	
	fL1V0[0] = fL1V0[1] = 0;
}

//_____________
AliEMCALTriggerData::~AliEMCALTriggerData()
{
}

//_____________
void AliEMCALTriggerData::SetL0Patches(Int_t i, const TClonesArray& patches)
{
	Int_t new_size = patches.GetEntriesFast();
	Int_t old_size = fL0Patches->GetEntriesFast();

	fL0NPatches[i] = new_size;
	
	Int_t size = 0;
	for (Int_t j=0;j<=i;j++) size += fL0NPatches[j];
		
	fL0Patches->Expand( size );
			
	for (Int_t j=0;j<new_size;j++)
	{
		AliEMCALTriggerPatch* p = static_cast<AliEMCALTriggerPatch*>( patches.At(j) );
		new((*fL0Patches)[old_size+j]) AliEMCALTriggerPatch( *p );
	}
}

//_____________
void AliEMCALTriggerData::SetL1GammaPatches(const TClonesArray& patches)
{
	Int_t size = patches.GetEntriesFast();
	fL1GammaPatches->Expand( size );

	for (Int_t j=0;j<size;j++)
	{
		AliEMCALTriggerPatch* p = static_cast<AliEMCALTriggerPatch*>( patches.At(j) );
		new((*fL1GammaPatches)[j]) AliEMCALTriggerPatch( *p );
	}
}

//_____________
void AliEMCALTriggerData::SetL1JetPatches(const TClonesArray& patches)
{
	Int_t size = patches.GetEntriesFast();
	
	fL1JetPatches->Expand( size );
	
	for (Int_t j=0;j<size;j++)
	{
		AliEMCALTriggerPatch* p = static_cast<AliEMCALTriggerPatch*>( patches.At(j) );
		new((*fL1JetPatches)[j]) AliEMCALTriggerPatch( *p );
	}
}

//_____________
void AliEMCALTriggerData::SetL1Region(Int_t**& region)
{
	//
  	for (Int_t i=0;i<48;i++)
  		for (Int_t j=0;j<64;j++)
		{
			fL1Region[i][j] = region[i][j];
		}
}

//_____________
void AliEMCALTriggerData::SetL1V0(const Int_t*& arr)
{
	for (Int_t i=0;i<2;i++) fL1V0[i] = arr[i];
}

//_____________
void AliEMCALTriggerData::Scan() const
{
	//
	printf("L0:\n");
	for (Int_t i=0;i<32;i++) printf("\tFound %2d patches in TRU %2d\n",fL0NPatches[i],i);
	
	printf("L1:\n");
	printf("\tRegion of size.....................(%2d,%2d)\n",int(fL1RegionSize.X()),int(fL1RegionSize.Y()));
	printf("\tGamma sub-region size..............(%2d,%2d)\n",int(fL1GammaSubRegionSize.X()),int(fL1GammaSubRegionSize.Y()));
	printf("\tJet sub-region size................(%2d,%2d)\n",int(fL1JetSubRegionSize.X()),int(fL1JetSubRegionSize.Y()));
	printf("\tFound %4d gamma patches of size...(%2d,%2d)\n",fL1GammaPatches->GetEntriesFast(),int(fL1GammaPatchSize.X()),int(fL1GammaPatchSize.Y()));
	printf("\tFound %4d jet patches of size.....(%2d,%2d)\n",fL1JetPatches->GetEntriesFast(),int(fL1JetPatchSize.X()),int(fL1JetPatchSize.Y()));
}

//_____________
void AliEMCALTriggerData::Reset()
{
	//
    if (fL0Patches)           fL0Patches->Delete(); 
  	if (fL1GammaPatches) fL1GammaPatches->Delete();
	if (fL1JetPatches)     fL1JetPatches->Delete();	

	for (Int_t i=0;i<32;i++) fL0NPatches[i] = 0;
 	for (Int_t i=0;i<48;i++) for (Int_t j=0;j<64;j++) fL1Region[i][j] = 0;
 	fL1V0[0] = fL1V0[1] = 0;
}



