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
#include "AliVCaloCells.h"

ClassImp(AliEmcalTriggerMaker)

using namespace std;

//________________________________________________________________________
AliEmcalTriggerMaker::AliEmcalTriggerMaker() : 
  AliAnalysisTaskEmcal("AliEmcalTriggerMaker",kFALSE),
  fCaloTriggersOutName("EmcalTriggers"),
  fCaloTriggerSetupOutName("EmcalTriggersSetup"),
  fCaloTriggersOut(0),
  fCaloTriggerSetupOut(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalTriggerMaker::AliEmcalTriggerMaker(const char *name) : 
  AliAnalysisTaskEmcal(name,kFALSE),
  fCaloTriggersOutName("EmcalTriggers"),
  fCaloTriggerSetupOutName("EmcalTriggersSetup"),
  fCaloTriggersOut(0),
  fCaloTriggerSetupOut(0)
{
  // Constructor.

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

}

//________________________________________________________________________
Bool_t AliEmcalTriggerMaker::Run() 
{
  // Create the emcal particles

  Int_t globCol, globRow, tBits, cellAbsId[4];
  Int_t absId, adcAmp;
  Int_t i, j, k, iMain, cmCol, cmRow, cmiCellCol, cmiCellRow;
  Int_t jetTrigger, iTriggers;
  Double_t amp, ca, eMain, cmiCol, cmiRow;
  
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
  
  // must reset before usage, or the class will fail 
  fCaloTriggers->Reset();
  
  // dig out common data (thresholds)
  // 0 - jet high, 1 - gamma high, 2 - jet low, 3 - gamma low
  fCaloTriggerSetupOut->SetThresholds( fCaloTriggers->GetL1Threshold( 0 ),
                                       fCaloTriggers->GetL1Threshold( 1 ),
                                       fCaloTriggers->GetL1Threshold( 2 ),
                                       fCaloTriggers->GetL1Threshold( 3 ));

  // class is not empty
  if( fCaloTriggers->GetEntries() > 0 ){
 
    iTriggers = 0;
    iMain = -1;
    eMain = -1;

    // save primary vertex in vector
    vertex.SetXYZ( fVertex[0], fVertex[1], fVertex[2] );

    // go throuth the trigger channels
    while( fCaloTriggers->Next() ){
      
      // check if jet trigger low or high
      fCaloTriggers->GetTriggerBits( tBits );
      
      jetTrigger = 0;
      if(( tBits >> ( kTriggerTypeEnd + kL1JetLow )) & 1 )
        jetTrigger = 1;
      if(( tBits >> ( kTriggerTypeEnd + kL1JetHigh )) & 1)
        jetTrigger = jetTrigger | 2;
      
      if( jetTrigger == 0 )
        continue;
      
      // get position in global 2x2 tower coordinates
      // A0 left bottom (0,0)
      fCaloTriggers->GetPosition( globCol, globRow );

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
    
      // get ADC amplitude
      fCaloTriggers->GetL1TimeSum( adcAmp );

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
      
      // check if more energetic than others for main patch marking
      if( eMain < amp ){
        eMain = amp;
        iMain = iTriggers - 1;
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
    if( iMain > -1 ){
      trigger = (AliEmcalTriggerPatchInfo*)fCaloTriggersOut->At( iMain );
      tBits = trigger->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << 24 );
      trigger->SetTriggerBits( tBits );
    }
    
  } // there are some triggers

  return kTRUE;
}
