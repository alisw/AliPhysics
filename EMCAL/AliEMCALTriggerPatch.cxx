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

/*

 

Patch object implementation: one patch is made of subregions (Olivier's nomenclature)
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerPatch.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliEMCALGeometry.h"
#include "AliEMCAL.h"

#include "TArrayI.h"

ClassImp(AliEMCALTriggerPatch)

//____________
AliEMCALTriggerPatch::AliEMCALTriggerPatch() : TObject(),
fPosition(0x0),
fSum(0)
{
	// Default constructor
}

//____________
AliEMCALTriggerPatch::AliEMCALTriggerPatch( Int_t i, Int_t j,  Int_t k ) : TObject(),
fPosition(new TVector2( i , j )),
fSum(k)
{
}

//____________

//____________________________________________________________________
AliEMCALTriggerPatch::AliEMCALTriggerPatch(const AliEMCALTriggerPatch& other) : TObject(other), 
fPosition( new TVector2(*other.fPosition) ),
fSum( other.fSum )
{	
	// Copy ctor
}

//____________
AliEMCALTriggerPatch::~AliEMCALTriggerPatch()
{	
	if (fPosition) delete fPosition;
}

//____________
void AliEMCALTriggerPatch::Print(const Option_t*) const
{
	printf("]> Patch at (%2d , %2d) w/ sum %3d\n",
		   (int)fPosition->X(),(int)fPosition->Y(),fSum); 
}

//________________
void AliEMCALTriggerPatch::GetAbsCellIdsFromPatchPosition( TVector2& pSize, TVector2& sSize, TArrayI& absid )
{
	AliRunLoader*     runLoader = AliRunLoader::Instance();
	AliEMCALGeometry*      geom = dynamic_cast<AliEMCAL*>(runLoader->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
	
	Int_t nTowersinpatch = (Int_t) (pSize.X() * pSize.Y() * sSize.X() * sSize.Y() * 4);
	
	absid.Set( nTowersinpatch );
	
	// fPosition: patch position in the STU region
	Int_t ix = (Int_t)(( fPosition->X() + pSize.X() ) * sSize.X()); 
	Int_t iy = (Int_t)(( fPosition->Y() + pSize.Y() ) * sSize.Y());
	
	Int_t it = 0;
	
	for (Int_t i=(Int_t) (fPosition->X() * sSize.X()); i<ix; i++) // Loop over subregions FastOR
	{
	  for (Int_t j=(Int_t) (fPosition->Y() * sSize.Y()); j<iy; j++) 
		{
			Int_t nSupMod = int(i/24) + 2 * int(j/12);
			
			Int_t nModule = 0;
			
			if ( nSupMod<10 )
				nModule = ( i < 24 ) ? 12 * ( 23 - i ) + 11 - j%12 : 12 * ( i%24 ) + 11 - j%12;
			else
				nModule = ( i < 24 ) ?  6 * ( 23 - i ) +  5 - j%12 :  6 * ( i%24 ) +  5 - j;
			
			// Convert (TRU,eta,phi) to Id 
			for (Int_t k=0;k<2;k++)
			{
				for (Int_t l=0;l<2;l++)
				{
					Int_t iphi, ieta;
					geom->GetCellPhiEtaIndexInSModule(nSupMod, nModule, k, l, iphi, ieta);

					absid.SetAt( geom->GetAbsCellIdFromCellIndexes(nSupMod, iphi, ieta) , it++ );
				}
			}
		}
	}
}
