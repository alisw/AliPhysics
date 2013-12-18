// $Id$
//
// Class to make emcal particles in AOD/ESD events.
//
// Author: J.Kral

#include <TClonesArray.h>
#include <iostream>

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
  fV0(0)
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
  fV0(0)
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
  Int_t globCol, globRow, tBits, cellAbsId[4], v0[2];
  Int_t absId, adcAmp;
  Int_t i, j, k, iMain, iMainSimple, cmCol, cmRow, cmiCellCol, cmiCellRow, nCell, iCell;
  Int_t jetTrigger, iTriggers;
	Int_t patchADC[48][64];
  Double_t amp, ca, eMain, eMainSimple, cmiCol, cmiRow;
  ULong64_t v0S, thresh;
  Bool_t isOfflineSimple;
  
  TVector3 centerGeo, center1, center2, centerMass, edge1, edge2, vertex;
  
  AliEmcalTriggerPatchInfo *trigger;

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

  // must reset before usage, or the class will fail 
  fCaloTriggers->Reset();
  for (i=0; i<2; i++) {
    fEGA[i] = 0;
    fEJE[i] = 0;
  }
  // first run over the patch array to compose a map of 2x2 patch energies
  // which is then needed to construct the full patch ADC energy
  // class is not empty
  Int_t isMC = 0;
  if (MCEvent()) isMC = 1;
  Int_t offSet = (1 - isMC) * kTriggerTypeEnd;

  if( fCaloTriggers->GetEntries() > 0 ){
		
    // zero the array
    for( i = 0; i < 48; i++ )
      for( j = 0; j < 64; j++ )
	patchADC[i][j] = 0;
		
    // go throuth the trigger channels
    while( fCaloTriggers->Next() ){
      // get position in global 2x2 tower coordinates
      // A0 left bottom (0,0)
      fCaloTriggers->GetPosition( globCol, globRow );

      // for some strange reason some ADC amps are initialized in reconstruction
      // as -1, neglect those :\\ wth
      fCaloTriggers->GetL1TimeSum( adcAmp );
			if( adcAmp > -1 )
				patchADC[globCol][globRow] = adcAmp;
      
      fCaloTriggers->GetTriggerBits( tBits );
      
      if (tBits) {
        if ((tBits >> (offSet + kL1GammaHigh)) & 1 ) fEGA[0] = 1;
        if ((tBits >> (offSet + kL1GammaLow )) & 1 ) fEGA[1] = 1;
        if ((tBits >> (offSet + kL1JetHigh  )) & 1 ) fEJE[0] = 1;
        if ((tBits >> (offSet + kL1JetLow   )) & 1 ) fEJE[1] = 1;
      }
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
 
    iTriggers = 0;
    iMain = -1;
    eMain = -1;
    iMainSimple = -1;
    eMainSimple = -1;

    // save primary vertex in vector
    vertex.SetXYZ( fVertex[0], fVertex[1], fVertex[2] );

    // go throuth the trigger channels, real first, then offline
    while( NextTrigger( isOfflineSimple ) ){
      // check if jet trigger low or high
      if( ! isOfflineSimple )
	fCaloTriggers->GetTriggerBits( tBits );
      else
        fSimpleOfflineTriggers->GetTriggerBits( tBits );
      
      jetTrigger = 0;
      if(( tBits >> ( offSet + kL1JetLow )) & 1 )
        jetTrigger = 1;
      if(( tBits >> ( offSet + kL1JetHigh )) & 1)
        jetTrigger = jetTrigger | 2;
      
      if( jetTrigger == 0 )
        continue;
      
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
      for( i = 0; i < 16; i++ ){
        for( j = 0; j < 16; j++ ){
          // get the 4 cells composing the trigger channel
          fGeom->GetAbsFastORIndexFromPositionInEMCAL( globCol+i, globRow+j, absId );
          fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );
          // add amplitudes and find patch edges
          for( k = 0; k < 4; k++ ){
            ca = fCaloCells->GetCellAmplitude( cellAbsId[k] );
            amp += ca;
            cmiCol += ca*(Double_t)i;
            cmiRow += ca*(Double_t)j;
          }
          // add the STU ADCs in the patch
          if( ! isOfflineSimple )
          adcAmp += patchADC[globCol+i][globRow+j];
          else
            adcAmp += fPatchADCSimple[globCol+i][globRow+j];
        }
      } // 32x32 cell window
      if( amp == 0 ){
        AliDebug(2,"EMCal trigger patch with 0 energy.");
        continue;
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
      fGeom->GetAbsFastORIndexFromPositionInEMCAL( globCol+15, globRow+15, absId );
      // convert to the 4 absId of the cells composing the trigger channel
      fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );

      fGeom->GetGlobal( cellAbsId[3], edge2 );
      
      // get the geometrical center as an average of two diagonally
      // adjacent patches in the center
      // picking two diagonally closest cells from the patches
      fGeom->GetAbsFastORIndexFromPositionInEMCAL( globCol+7, globRow+7, absId );
      fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );
      fGeom->GetGlobal( cellAbsId[3], center1 );
      
      fGeom->GetAbsFastORIndexFromPositionInEMCAL( globCol+8, globRow+8, absId );
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
    
      // save the trigger object
      new ((*fCaloTriggersOut)[iTriggers])AliEmcalTriggerPatchInfo();
      trigger = (AliEmcalTriggerPatchInfo*)fCaloTriggersOut->At( iTriggers );
      iTriggers++;
      
      trigger->SetCenterGeo( centerGeo, amp );
      trigger->SetCenterMass( centerMass, amp );
      trigger->SetEdge1( edge1, amp );
      trigger->SetEdge2( edge2, amp );
      trigger->SetADCAmp( adcAmp );
      trigger->SetTriggerBits( tBits );
      trigger->SetEdgeCell( globCol*2, globRow*2 ); // from triggers to cells
      trigger->SetOffSet(offSet);

      // check if more energetic than others for main patch marking
      if( ! isOfflineSimple && eMain < amp ){
        eMain = amp;
        iMain = iTriggers - 1;
      }
      if( isOfflineSimple && eMainSimple < amp ){
        eMainSimple = amp;
        iMainSimple = iTriggers - 1;
      }
    } // triggers
    
    // mark the most energetic patch as main
    // for real and also simple offline
    if( iMain > -1 ){
      trigger = (AliEmcalTriggerPatchInfo*)fCaloTriggersOut->At( iMain );
      tBits = trigger->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << 24 );
      trigger->SetTriggerBits( tBits );
    }
    if( iMainSimple > -1 ){
      trigger = (AliEmcalTriggerPatchInfo*)fCaloTriggersOut->At( iMainSimple );
      tBits = trigger->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << 24 );
      trigger->SetTriggerBits( tBits );
    }
    
  } // there are some triggers

  return kTRUE;
}

//________________________________________________________________________
void AliEmcalTriggerMaker::RunSimpleOfflineTrigger() 
{
  // runs simple offline trigger algorithm

  Int_t i, j, k, l, tBits, tSum;
  TArrayI tBitsArray, rowArray, colArray;
  
  Int_t isMC = 0;
  if (MCEvent()) isMC = 1;
  Int_t offSet = (1 - isMC) * kTriggerTypeEnd;

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
        tBits = tBits | ( 1 << ( offSet + kL1JetLow ));
      if( tSum > fCaloTriggerSetupOut->GetThresholdJetHighSimple() )
        tBits = tBits | ( 1 << ( offSet + kL1JetHigh ));
      
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
