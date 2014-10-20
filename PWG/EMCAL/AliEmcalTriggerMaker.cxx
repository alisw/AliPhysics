// $Id: AliEmcalTriggerMaker.cxx 64593 2013-10-18 10:23:58Z loizides $
//
// Class to make emcal particles in AOD/ESD events.
//
// Author: J.Kral

#include <TClonesArray.h>

#include "AliLog.h"
#include "AliEmcalTriggerPatchInfo.h"
#include "AliEmcalTriggerSetupInfo.h"

#include "AliEmcalTriggerMaker.h"
#include "AliEMCALTriggerTypes.h"
#include "AliEMCALGeometry.h"

#include "AliVCaloTrigger.h"
#include "AliAODCaloTrigger.h"
#include "AliVCaloCells.h"
#include "AliVVZERO.h"

ClassImp(AliEmcalTriggerMaker)

using namespace std;

//________________________________________________________________________
AliEmcalTriggerMaker::AliEmcalTriggerMaker() : 
  AliAnalysisTaskEmcal("AliEmcalTriggerMaker",kFALSE),
  fCaloTriggersOutName("EmcalTriggers"),
  fCaloTriggerSetupOutName("EmcalTriggersSetup"),
  fV0InName("AliAODVZERO"),
  fCaloTriggersOut(0),
  fCaloTriggerSetupOut(0),
  fSimpleOfflineTriggers(0),
  fV0(0),
  fITrigger(0)
{
  // Constructor.
  for( int i = 0; i < 4; i++ )
    for( int j = 0; j < 3; j++ )
      fThresholdConstants[i][j] = 0;
}

//________________________________________________________________________
AliEmcalTriggerMaker::AliEmcalTriggerMaker(const char *name) : 
  AliAnalysisTaskEmcal(name,kFALSE),
  fCaloTriggersOutName("EmcalTriggers"),
  fCaloTriggerSetupOutName("EmcalTriggersSetup"),
  fV0InName("AliAODVZERO"),
  fCaloTriggersOut(0),
  fCaloTriggerSetupOut(0),
  fSimpleOfflineTriggers(0),
  fV0(0),
  fITrigger(0)
{
  // Constructor.
  for( int i = 0; i < 4; i++ )
    for( int j = 0; j < 3; j++ )
      fThresholdConstants[i][j] = 0;
}

//________________________________________________________________________
AliEmcalTriggerMaker::~AliEmcalTriggerMaker()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalTriggerMaker::ExecOnce()
{
  // Init the analysis.

  AliAnalysisTaskEmcal::ExecOnce();

  if (!fInitialized)
    return;

  if (!fCaloTriggersOutName.IsNull()) {
    fCaloTriggersOut = new TClonesArray("AliEmcalTriggerPatchInfo");
    fCaloTriggersOut->SetName(fCaloTriggersOutName);

    if (!(InputEvent()->FindListObject(fCaloTriggersOutName))) {
      InputEvent()->AddObject(fCaloTriggersOut);
    }
    else {
      fInitialized = kFALSE;
      AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fCaloTriggersOutName.Data()));
      return;
    }
  }

  if (!fCaloTriggerSetupOutName.IsNull()) {
    fCaloTriggerSetupOut = new AliEmcalTriggerSetupInfo();
    fCaloTriggerSetupOut->SetName(fCaloTriggerSetupOutName);

    if (!(InputEvent()->FindListObject(fCaloTriggerSetupOutName))) {
      InputEvent()->AddObject(fCaloTriggerSetupOut);
    }
    else {
      fInitialized = kFALSE;
      AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fCaloTriggerSetupOutName.Data()));
      return;
    }
  }

  if ( ! fV0InName.IsNull()) {
    fV0 = (AliVVZERO*)InputEvent()->FindListObject( fV0InName );
  }

  // container for simple offline trigger processing
  fSimpleOfflineTriggers = new AliAODCaloTrigger();
  fSimpleOfflineTriggers->Allocate( 0 );
}

//________________________________________________________________________
Bool_t AliEmcalTriggerMaker::Run() 
{
  // Create the emcal particles

  Short_t cellId;
  Int_t globCol, globRow, v0[2];
  Int_t absId, adcAmp, tBits;
  Int_t i, j, nCell, iCell;
  Double_t amp;
  ULong64_t v0S, thresh;
  Bool_t isOfflineSimple;
  
  AliEmcalTriggerPatchInfo *trigger, *triggerMainJet, *triggerMainGamma, *triggerMainLevel0;
  AliEmcalTriggerPatchInfo *triggerMainJetSimple, *triggerMainGammaSimple;

  // delete patch array, clear setup object
  fCaloTriggersOut->Delete();
  fCaloTriggerSetupOut->Clean();

  if( !fCaloTriggers ){
    AliError(Form("Calo triggers container %s not available.", fCaloTriggersName.Data()));
    return kTRUE;
  }
  if( !fCaloCells ){
    AliError(Form("Calo cells container %s not available.", fCaloCellsName.Data()));
    return kTRUE;
  }
  if( !fCaloCells ){
    AliError(Form("V0 container %s not available.", fV0InName.Data()));
    return kTRUE;
  }
  
  // do not process, if sooner than 11h period
  // 160683 ??
  if( InputEvent()->GetRunNumber() < 167693 )
    return kTRUE;
 
//   // do not process any MC, since no MC was generated with correct
//   // EMCal trigger L1 jet trigger simulation, yet
//   // productions will be enabled, once some correct once are produced
//   if( MCEvent() != 0 )
//     return kTRUE;
  
  // must reset before usage, or the class will fail 
  fCaloTriggers->Reset();

	// first run over the patch array to compose a map of 2x2 patch energies
  // which is then needed to construct the full patch ADC energy
  // class is not empty
  if( fCaloTriggers->GetEntries() > 0 ){
		
		// zero the array
		for( i = 0; i < 48; i++ )
			for( j = 0; j < 64; j++ )
				fPatchADC[i][j] = 0;
		
    // go throuth the trigger channels
    while( fCaloTriggers->Next() ){
      // get position in global 2x2 tower coordinates
      // A0 left bottom (0,0)
      fCaloTriggers->GetPosition( globCol, globRow );

      // for some strange reason some ADC amps are initialized in reconstruction
      // as -1, neglect those :\\ wth
      fCaloTriggers->GetL1TimeSum( adcAmp );
			if( adcAmp > -1 )
				fPatchADC[globCol][globRow] = adcAmp;
      
		} // patches
	} // array not empty
  
  // fill the array for offline trigger processing
  // using calibrated cell energies
  for( i = 0; i < 48; i++ )
    for( j = 0; j < 64; j++ )
      fPatchADCSimple[i][j] = 0;

  // fill the patch ADCs from cells
  nCell = fCaloCells->GetNumberOfCells();
  
  for( iCell = 0; iCell < nCell; iCell++ ){
    // get the cell info, based in index in array
    cellId = fCaloCells->GetCellNumber( iCell );
    amp = fCaloCells->GetAmplitude( iCell );
    
    // get position
    fGeom->GetFastORIndexFromCellIndex( cellId, absId );
    fGeom->GetPositionInEMCALFromAbsFastORIndex( absId, globCol, globRow );
    
    // add
    fPatchADCSimple[globCol][globRow] += amp/kEMCL1ADCtoGeV;
  }

  // dig out common data (thresholds)
  // 0 - jet high, 1 - gamma high, 2 - jet low, 3 - gamma low
  fCaloTriggerSetupOut->SetThresholds( fCaloTriggers->GetL1Threshold( 0 ),
                                       fCaloTriggers->GetL1Threshold( 1 ),
                                       fCaloTriggers->GetL1Threshold( 2 ),
                                       fCaloTriggers->GetL1Threshold( 3 ));

  // get the V0 value and compute and set the offline thresholds
  // get V0, compute thresholds and save them as global parameters
  v0[0] = fV0->GetTriggerChargeA();
  v0[1] = fV0->GetTriggerChargeC();
  v0S = v0[0] + v0[1];
  
  fSimpleOfflineTriggers->SetL1V0( v0 );
  
  for( i = 0; i < 4; i++ ){
    // A*V0^2/2^32+B*V0/2^16+C
    thresh = ( ((ULong64_t)fThresholdConstants[i][0]) * v0S * v0S ) >> 32;
    thresh += ( ((ULong64_t)fThresholdConstants[i][1]) * v0S ) >> 16;
    thresh += ((ULong64_t)fThresholdConstants[i][2]);
    fSimpleOfflineTriggers->SetL1Threshold( i, thresh );
  }
  
  // save the thresholds in output object
  fCaloTriggerSetupOut->SetThresholdsSimple( fSimpleOfflineTriggers->GetL1Threshold( 0 ),
                                       fSimpleOfflineTriggers->GetL1Threshold( 1 ),
                                       fSimpleOfflineTriggers->GetL1Threshold( 2 ),
                                       fSimpleOfflineTriggers->GetL1Threshold( 3 ));
  
  // run the trigger
  RunSimpleOfflineTrigger();

  // reset for re-run
  fCaloTriggers->Reset();
  fSimpleOfflineTriggers->Reset();

  // class is not empty
  if( fCaloTriggers->GetEntries() > 0 ||  fSimpleOfflineTriggers->GetEntries() > 0 ){
    fITrigger = 0;
		triggerMainGamma = 0;
		triggerMainJet = 0;
		triggerMainGammaSimple = 0;
		triggerMainJetSimple = 0;
		triggerMainLevel0 = 0;

    // go throuth the trigger channels, real first, then offline
    while( NextTrigger( isOfflineSimple ) ){
			// process jet
			trigger = ProcessPatch( kTMEMCalJet, isOfflineSimple );
			
			// save main jet triggers in event
			if( trigger != 0 ){
				// check if more energetic than others for main patch marking
				if( ! isOfflineSimple ){
					if( triggerMainJet == 0 )
						triggerMainJet = trigger;
					else if( triggerMainJet->GetPatchE() < trigger->GetPatchE() )
						triggerMainJet = trigger;
				}
				else{
					if( triggerMainJetSimple == 0 )
						triggerMainJetSimple = trigger;
					else if( triggerMainJetSimple->GetPatchE() < trigger->GetPatchE() )
						triggerMainJetSimple = trigger;
				}
			}
      
			// process gamma
			trigger = ProcessPatch( kTMEMCalGamma, isOfflineSimple );
			
			// save main gamma triggers in event
			if( trigger != 0 ){
				// check if more energetic than others for main patch marking
				if( ! isOfflineSimple ){
					if( triggerMainGamma == 0 )
						triggerMainGamma = trigger;
					else if( triggerMainGamma->GetPatchE() < trigger->GetPatchE() )
						triggerMainGamma = trigger;
				}
				else{
					if( triggerMainGammaSimple == 0 )
						triggerMainGammaSimple = trigger;
					else if( triggerMainGammaSimple->GetPatchE() < trigger->GetPatchE() )
						triggerMainGammaSimple = trigger;
				}
			}

			// level 0 triggers
			trigger = ProcessPatch(kTMEMCalLevel0, isOfflineSimple);

			// save main level0 trigger in the event
			if(trigger){
				if(!triggerMainLevel0) triggerMainLevel0 = trigger;
				else if(triggerMainLevel0->GetPatchE() < trigger->GetPatchE()) triggerMainLevel0 = trigger;
			}

//       cout << " pi:" << trigger->GetPhiMin() << " px:" << trigger->GetPhiMax();
//       cout << " pg:" << trigger->GetPhiGeo() << " " << (trigger->GetPhiMin()+trigger->GetPhiMax()) / 2.;
//       cout << " pc:" << trigger->GetPhiCM();
//       cout << " ei:" << trigger->GetEtaMin() << " ex:" << trigger->GetEtaMax();
//       cout << " eg:" << trigger->GetEtaGeo() << " " << (trigger->GetEtaMin()+trigger->GetEtaMax()) / 2.;
//       cout << " ec:" << trigger->GetEtaCM();
//       cout << " e:" << trigger->GetPatchE();
//       cout << " jl:" << trigger->IsJetLow() << " jh:" << trigger->IsJetHigh() << endl;
      
    } // triggers
    
    // mark the most energetic patch as main
    // for real and also simple offline
    if( triggerMainJet != 0 ){
      tBits = triggerMainJet->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << 24 );
      triggerMainJet->SetTriggerBits( tBits );
    }
    if( triggerMainJetSimple != 0 ){
      tBits = triggerMainJetSimple->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << 24 );
      triggerMainJetSimple->SetTriggerBits( tBits );
    }
    if( triggerMainGamma != 0 ){
      tBits = triggerMainGamma->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << 24 );
      triggerMainGamma->SetTriggerBits( tBits );
    }
    if( triggerMainGammaSimple != 0 ){
      tBits = triggerMainGammaSimple->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << 24 );
      triggerMainGammaSimple->SetTriggerBits( tBits );
    }
    if(triggerMainLevel0){
    	tBits = triggerMainLevel0->GetTriggerBits();
    	// main trigger flag
    	tBits |= (1 << 24);
    	triggerMainLevel0->SetTriggerBits(tBits);
    }
    
  } // there are some triggers

  return kTRUE;
}

//________________________________________________________________________
AliEmcalTriggerPatchInfo* AliEmcalTriggerMaker::ProcessPatch( TriggerMakerTriggerType_t type, Bool_t isOfflineSimple ){

	
  Int_t globCol, globRow, tBits, cellAbsId[4], posOffset(0);
  Int_t absId, adcAmp;
  Int_t i, j, k, cmCol, cmRow, cmiCellCol, cmiCellRow;
  Double_t amp, ca, cmiCol, cmiRow;
  
  TVector3 centerGeo, center1, center2, centerMass, edge1, edge2, vertex, cellCoor;
  
  AliEmcalTriggerPatchInfo *trigger;

	// check if jet trigger low or high
	if( ! isOfflineSimple )
		fCaloTriggers->GetTriggerBits( tBits );
	else
		fSimpleOfflineTriggers->GetTriggerBits( tBits );
	
	if(( type == kTMEMCalJet && ! IsEJE( tBits )) || ( type == kTMEMCalGamma && ! IsEGA( tBits )) || (type == kTMEMCalLevel0 && !IsLevel0(tBits)))
		return 0;
	
	// save primary vertex in vector
	vertex.SetXYZ( fVertex[0], fVertex[1], fVertex[2] );

	// get position in global 2x2 tower coordinates
	// A0 left bottom (0,0)
	if( ! isOfflineSimple )
		fCaloTriggers->GetPosition( globCol, globRow );
	else
		fSimpleOfflineTriggers->GetPosition( globCol, globRow );

	// get the absolute trigger ID
	fGeom->GetAbsFastORIndexFromPositionInEMCAL( globCol, globRow, absId );
	// convert to the 4 absId of the cells composing the trigger channel
	fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );
	
	// get low left edge (eta max, phi min)
	fGeom->GetGlobal( cellAbsId[0], edge1 );
	
	// sum the available energy in the 32/32 window of cells
	// step over trigger channels and get all the corresponding cells
	// make CM
	amp = 0;
	cmiCol = 0;
	cmiRow = 0;
	adcAmp = 0;
	if( IsEJE( tBits ) && type == kTMEMCalJet ){
		for( i = 0; i < 16; i++ ){
			for( j = 0; j < 16; j++ ){
				// get the 4 cells composing the trigger channel
				fGeom->GetAbsFastORIndexFromPositionInEMCAL( globCol+i, globRow+j, absId );
				fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );
				// add amplitudes and find patch edges
				for( k = 0; k < 4; k++ ){
					ca = fCaloCells->GetCellAmplitude( cellAbsId[k] );
					fGeom->GetGlobal( cellAbsId[k], cellCoor );
					amp += ca;
					cmiCol += ca*(Double_t)i;
					cmiRow += ca*(Double_t)j;
				}
				// add the STU ADCs in the patch
				if( ! isOfflineSimple )
					adcAmp += fPatchADC[globCol+i][globRow+j];
				else
					adcAmp += fPatchADCSimple[globCol+i][globRow+j];
			}
		} // 32x32 cell window
	}
	
	if( (IsEGA( tBits ) && type == kTMEMCalGamma) || (IsLevel0( tBits ) && type == kTMEMCalLevel0)){
		for( i = 0; i < 2; i++ ){
			for( j = 0; j < 2; j++ ){
				// get the 4 cells composing the trigger channel
				fGeom->GetAbsFastORIndexFromPositionInEMCAL( globCol+i, globRow+j, absId );
				fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );
				// add amplitudes and find patch edges
				for( k = 0; k < 4; k++ ){
					ca = fCaloCells->GetCellAmplitude( cellAbsId[k] );
					fGeom->GetGlobal( cellAbsId[k], cellCoor );
					amp += ca;
					cmiCol += ca*(Double_t)i;
					cmiRow += ca*(Double_t)j;
				}
				// add the STU ADCs in the patch
				if( ! isOfflineSimple )
					adcAmp += fPatchADC[globCol+i][globRow+j];
				else
					adcAmp += fPatchADCSimple[globCol+i][globRow+j];
			}
		} // 4x4 cell window
	}

	if( amp == 0 ){
		AliDebug(2,"EMCal trigger patch with 0 energy.");
		return 0;
	}
	
	// get the CM and patch index
	cmiCol /= amp;
	cmiRow /= amp;
	cmCol = globCol + (Int_t)cmiCol;
	cmRow = globRow + (Int_t)cmiRow;

	// get the patch and corresponding cells
	fGeom->GetAbsFastORIndexFromPositionInEMCAL( cmCol, cmRow, absId );
	fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );

	// find which out of the 4 cells is closest to CM and get it's position
	cmiCellCol = TMath::Nint( cmiCol * 2. );
	cmiCellRow = TMath::Nint( cmiRow * 2. );
	fGeom->GetGlobal( cellAbsId[(cmiCellRow%2)*2 + cmiCellCol%2], centerMass );
	
	// get up right edge (eta min, phi max)
	// get the absolute trigger ID
	switch(type){
	case kTMEMCalJet:
		posOffset = 15;
		break;
	case kTMEMCalGamma:
		posOffset = 1;
		break;
	case kTMEMCalLevel0:
		posOffset = 1;
		break;
	default:
		posOffset = 0;
		break;
	};
	fGeom->GetAbsFastORIndexFromPositionInEMCAL( globCol+posOffset, globRow+posOffset, absId );
	// convert to the 4 absId of the cells composing the trigger channel
	fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );

	fGeom->GetGlobal( cellAbsId[3], edge2 );
	
	// get the geometrical center as an average of two diagonally
	// adjacent patches in the center
	// picking two diagonally closest cells from the patches
	switch(type){
	case kTMEMCalJet:
		posOffset = 7;
		break;
	case kTMEMCalGamma:
		posOffset = 0;
		break;
	case kTMEMCalLevel0:
		posOffset = 0;
		break;
	default:
		posOffset = 0;
		break;
	};
	fGeom->GetAbsFastORIndexFromPositionInEMCAL( globCol+posOffset, globRow+posOffset, absId );


	fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );
	fGeom->GetGlobal( cellAbsId[3], center1 );
	
	switch(type){
	case kTMEMCalJet:
		posOffset = 8;
		break;
	case kTMEMCalGamma:
		posOffset = 1;
		break;
	case kTMEMCalLevel0:
		posOffset = 1;
		break;
	};
	fGeom->GetAbsFastORIndexFromPositionInEMCAL( globCol+posOffset, globRow+posOffset, absId );

	fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );
	fGeom->GetGlobal( cellAbsId[0], center2 );
	
	centerGeo = center1;
	centerGeo += center2;
	centerGeo *= 0.5;
	
	// relate all to primary vertex
	centerGeo -= vertex;
	centerMass -= vertex;
	edge1 -= vertex;
	edge2 -= vertex;
	
	// fix tbits .. remove the unwanted type triggers
	// for Jet and Gamma triggers we remove also the level 0 bit since it will be stored in the level 0 patch
	// for level 0 we remove all gamma and jet trigger bits
	switch(type){
	case kTMEMCalJet:
		tBits = tBits & ~( 1 << (kTriggerTypeEnd + kL1GammaLow) | 1 << (kTriggerTypeEnd + kL1GammaHigh) | 1 << (kL1GammaLow) | 1 << (kL1GammaHigh) |
				1 << (kTriggerTypeEnd + kL0) | 1 << (kL0));
		break;
	case kTMEMCalGamma:
		tBits = tBits & ~( 1 << (kTriggerTypeEnd + kL1JetLow) | 1 << (kTriggerTypeEnd + kL1JetHigh) | 1 << (kL1JetLow) | 1 << (kL1JetHigh) |
				1 << (kTriggerTypeEnd + kL0) | 1 << (kL0));
		break;
	case kTMEMCalLevel0:
		tBits = tBits & ~( 1 << (kTriggerTypeEnd + kL1JetLow) | 1 << (kTriggerTypeEnd + kL1JetHigh) | 1 << (kL1JetLow) | 1 << (kL1JetHigh) |
				1 << (kTriggerTypeEnd + kL1GammaLow) | 1 << (kTriggerTypeEnd + kL1GammaHigh) | 1 << (kL1GammaLow) | 1 << (kL1GammaHigh));
		break;
	};

	// save the trigger object
	new ((*fCaloTriggersOut)[fITrigger])AliEmcalTriggerPatchInfo();
	trigger = (AliEmcalTriggerPatchInfo*)fCaloTriggersOut->At( fITrigger );
	fITrigger++;
	
	Int_t isMC = MCEvent() ? 1 : 0;
	Int_t offSet = (1 - isMC) * kTriggerTypeEnd;

	trigger->SetCenterGeo( centerGeo, amp );
	trigger->SetCenterMass( centerMass, amp );
	trigger->SetEdge1( edge1, amp );
	trigger->SetEdge2( edge2, amp );
	trigger->SetADCAmp( adcAmp );
	trigger->SetTriggerBits( tBits );
	trigger->SetOffSet(offSet);
	trigger->SetEdgeCell( globCol*2, globRow*2 ); // from triggers to cells
	
	return trigger;
	
}


//________________________________________________________________________
void AliEmcalTriggerMaker::RunSimpleOfflineTrigger() 
{
  // runs simple offline trigger algorithm
  
  // it creates separate patches for jet and gamma triggers
  // on the same positions (different from STU reconstruction behavior)
  // TODO:: change to merge

  Int_t i, j, k, l, tBits, tSum;
  TArrayI tBitsArray, rowArray, colArray;
  
  // 0 thresholds = no processing
  if( fCaloTriggerSetupOut->GetThresholdJetLowSimple() == 0 &&
    fCaloTriggerSetupOut->GetThresholdJetHighSimple() == 0 )
    return;
  
  // run the trigger algo, stepping by 8 towers (= 4 trigger channels)
  for( i = 0; i < 32; i += 4 ){
    for( j = 0; j < 48; j += 4 ){
      
      tSum = 0;
      tBits = 0;
      
      // window
      for( k = 0; k < 16; k++ )
        for( l = 0; l < 16; l++ )
          tSum += (ULong64_t)fPatchADCSimple[i+k][j+l];
      
      // check thresholds
      if( tSum > fCaloTriggerSetupOut->GetThresholdJetLowSimple() )
        tBits = tBits | ( 1 << ( kTriggerTypeEnd + kL1JetLow ));
      if( tSum > fCaloTriggerSetupOut->GetThresholdJetHighSimple() )
        tBits = tBits | ( 1 << ( kTriggerTypeEnd + kL1JetHigh ));
      
      // add trigger values
      if( tBits != 0 ){
        // add offline bit
        tBits = tBits | ( 1 << 25 );
        
        tBitsArray.Set( tBitsArray.GetSize() + 1 );
        colArray.Set( colArray.GetSize() + 1 );
        rowArray.Set( rowArray.GetSize() + 1 );
        
        tBitsArray[tBitsArray.GetSize()-1] = tBits;
        colArray[colArray.GetSize()-1] = i;
        rowArray[rowArray.GetSize()-1] = j;
      }
    }
  } // trigger algo
  
  // 4x4 trigger algo, stepping by 2 towers (= 1 trigger channel)
  for( i = 0; i < 46; i ++ ){
    for( j = 0; j < 62; j ++ ){
      
      tSum = 0;
      tBits = 0;
      
      // window
      for( k = 0; k < 2; k++ )
        for( l = 0; l < 2; l++ )
          tSum += (ULong64_t)fPatchADCSimple[i+k][j+l];
      
      // check thresholds
      if( tSum > fCaloTriggerSetupOut->GetThresholdGammaLowSimple() )
        tBits = tBits | ( 1 << ( kTriggerTypeEnd + kL1GammaLow ));
      if( tSum > fCaloTriggerSetupOut->GetThresholdGammaHighSimple() )
        tBits = tBits | ( 1 << ( kTriggerTypeEnd + kL1GammaHigh ));
      
      // add trigger values
      if( tBits != 0 ){
        // add offline bit
        tBits = tBits | ( 1 << 25 );
        
        tBitsArray.Set( tBitsArray.GetSize() + 1 );
        colArray.Set( colArray.GetSize() + 1 );
        rowArray.Set( rowArray.GetSize() + 1 );
        
        tBitsArray[tBitsArray.GetSize()-1] = tBits;
        colArray[colArray.GetSize()-1] = i;
        rowArray[rowArray.GetSize()-1] = j;
      }
    }
  } // trigger algo
  
  
  // save in object
  fSimpleOfflineTriggers->DeAllocate();
  fSimpleOfflineTriggers->Allocate( tBitsArray.GetSize() );

  for( i = 0; i < tBitsArray.GetSize(); i++ ){
    fSimpleOfflineTriggers->Add( colArray[i], rowArray[i], 0, 0, 0, 0, 0, tBitsArray[i] );
  }
  
}

//________________________________________________________________________
Bool_t AliEmcalTriggerMaker::NextTrigger( Bool_t &isOfflineSimple ) 
{
  // get next trigger
  
  Bool_t loopContinue;
  
  loopContinue = fCaloTriggers->Next();
  isOfflineSimple = kFALSE;

  if( ! loopContinue ){
    loopContinue = fSimpleOfflineTriggers->Next();
    isOfflineSimple = kTRUE;
  }
  
  return loopContinue;
}
