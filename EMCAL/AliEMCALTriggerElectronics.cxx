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
 
 
EMCal trigger electronics manager L0/L1
can handle both simulated digits and raw data
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerElectronics.h"
#include "AliEMCALTriggerTRU.h"
#include "AliEMCALTriggerSTU.h"
#include "AliEMCALGeometry.h"
#include "AliRunLoader.h"
#include "AliEMCAL.h" 
#include "AliRun.h" 
#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerData.h"
#include "AliEMCALDigit.h"
#include "AliCaloRawStreamV3.h"
#include "AliEMCALTriggerSTURawStream.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRawDigit.h"

#include <TVector2.h>
#include <TClonesArray.h>

namespace
{
	const Int_t kNTRU = 32;
}

ClassImp(AliEMCALTriggerElectronics)

//__________________
AliEMCALTriggerElectronics::AliEMCALTriggerElectronics(const AliEMCALTriggerDCSConfig *dcsConf) : TObject(),
fTRU(new TClonesArray("AliEMCALTriggerTRU",32)),
fSTU(0x0)
{
	TVector2 rSize;
	
	rSize.Set( 24.,  4. );

	// 32 TRUs
	for (Int_t i=0;i<kNTRU;i++) 
	{
		AliEMCALTriggerTRUDCSConfig* truConf = dcsConf->GetTRUDCSConfig(i);
		new ((*fTRU)[i]) AliEMCALTriggerTRU(truConf, rSize, int(i/3) % 2);
	}
	
	rSize.Set( 48., 64. );
	
	// 1 STU
	AliEMCALTriggerSTUDCSConfig* stuConf = dcsConf->GetSTUDCSConfig();
	fSTU = new AliEMCALTriggerSTU(stuConf, rSize);
	
	for (Int_t i=0;i<kNTRU;i++) fSTU->BuildMap( i, 
											  (static_cast<AliEMCALTriggerTRU*>(fTRU->At(i)))->Map(), 
											  (static_cast<AliEMCALTriggerTRU*>(fTRU->At(i)))->RegionSize() 
											  );
}

//________________
AliEMCALTriggerElectronics::~AliEMCALTriggerElectronics()
{
	//
	fTRU->Delete();
 	delete fSTU;
}

//__________________
void AliEMCALTriggerElectronics::Digits2Trigger(const TClonesArray* digits, const Int_t V0M[], AliEMCALTriggerData* data)
{
	//
	AliEMCALGeometry* geom = 0x0;
	
	AliRunLoader *rl = AliRunLoader::Instance();
	if (rl->GetAliRun() && rl->GetAliRun()->GetDetector("EMCAL"))
		geom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
	else 
		geom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());

	if (!geom) AliError("Cannot access geometry!");
	
	TIter NextDigit(digits);
	while (AliEMCALRawDigit* digit = (AliEMCALRawDigit*)NextDigit())
	{
		if ( digit )
		{
			Int_t id = digit->GetId();
			
//			digit->Print();
			
			Int_t iTRU, iADC;
			Bool_t isOK1 = geom->GetTRUFromAbsFastORIndex(id, iTRU, iADC);
			
			for (Int_t i = 0; i < digit->GetNSamples(); i++)
			{
				Int_t time, amp;
				Bool_t isOK2 = digit->GetTimeSample(i, time, amp);
				
				if (isOK1 && isOK2 && amp) (static_cast<AliEMCALTriggerTRU*>(fTRU->At(iTRU)))->SetADC(iADC, time, amp);
			}
		}
	}
	/*
	for (Int_t i=0; i<kNTRU; i++) 
	{
		printf("===========< TRU %2d >============\n",i);
		(static_cast<AliEMCALTriggerTRU*>(fTRU->At(i)))->Scan();
	}
	*/
	Int_t iL0 = 0;

	// At this point all FastOR are available for digitization
	// digitization is done in the TRU and produces time samples
	// Now run the trigger algo & consecutively write trigger outputs in TreeD dedicated branch

	for (Int_t i=0; i<kNTRU; i++) 
	{
		AliDebug(1,Form("===========< TRU %2d >============\n",i));
		
		AliEMCALTriggerTRU *iTRU = static_cast<AliEMCALTriggerTRU*>(fTRU->At(i));

	  	iL0 += iTRU->L0();

//		Int_t vL0Peaks[96][2]; iTRU->Peaks( vL0Peaks );
			
	   	data->SetL0Patches( i , iTRU->Patches() );			
//		data->SetL0Peaks(   i , vL0Peaks );

		if ( !i ) // do it once since identical for all TRU 
		{
			data->SetL0RegionSize(    *iTRU->RegionSize()    );
			data->SetL0SubRegionSize( *iTRU->SubRegionSize() );
	        data->SetL0PatchSize(     *iTRU->PatchSize()     );
		} 
		
		// 		if ( i == 31 ) i = 35;
		//
		// 		if ( ( i / 3 ) % 2 ) {
		// 			TRU->Print( 15 - 2 + ( i - int( i / 3 ) * 3 ) - 3 * ( (i / 3) / 2 ) , runLoader->GetEventNumber() );
		// 			printf("print data of TRU: from %2d to %2d\n",i,15 - 2 + ( i - int( i / 3 ) * 3 ) - 3 * ( (i / 3) / 2));
		// 		}
		// 		else
		// 		{
		// 			TRU->Print( 31 - i % 3 - 3 * ( (i / 3) / 2 ) , runLoader->GetEventNumber() );
		// 			printf("print data of TRU: from %2d to %2d\n",i,31 - i % 3 - 3 * ( (i / 3) / 2 ));
		// 		}		
	}

	// A L0 has been issued, run L1
	if ( iL0 ) 
	{
		for (Int_t i=0; i<kNTRU; i++) fSTU->FetchFOR( i, 
												    (static_cast<AliEMCALTriggerTRU*>(fTRU->At(i)))->Region(), 
												    (static_cast<AliEMCALTriggerTRU*>(fTRU->At(i)))->RegionSize() 
												    );
		
		fSTU->SetV0Multiplicity( V0M , 2 ); // C/A

		TVector2 size;

		size.Set( 1. , 1. );
		fSTU->SetSubRegionSize( size ); data->SetL1GammaSubRegionSize( size );
		
		size.Set( 2. , 2. );
		fSTU->SetPatchSize( size );     data->SetL1GammaPatchSize(     size );

		fSTU->L1( kGamma );

		data->SetL1GammaPatches(       fSTU->Patches()       );

		fSTU->Reset();

		size.Set( 4. , 4. ); 
		fSTU->SetSubRegionSize( size ); data->SetL1JetSubRegionSize( size );
		
		size.Set( 2. , 2. ); 
		fSTU->SetPatchSize( size );     data->SetL1JetPatchSize(     size );

		fSTU->L1( kJet );

		data->SetL1JetPatches(         fSTU->Patches()       );
		data->SetL1RegionSize(        *fSTU->RegionSize()    );

		Int_t** region = fSTU->Region();
		data->SetL1Region(      region      );
		const Int_t* mv0 = fSTU->V0(); 
		data->SetL1V0(          mv0         );
	}

	if ( AliDebugLevel() ) data->Scan();
	
	// Now reset the electronics for a fresh start with next event
	Reset();
}

//__________________
void AliEMCALTriggerElectronics::Reset()
{
	//
	TIter NextTRU(fTRU);
	while ( AliEMCALTriggerTRU *TRU = (AliEMCALTriggerTRU*)NextTRU() ) TRU->Reset();
	
	fSTU->Reset();
}
