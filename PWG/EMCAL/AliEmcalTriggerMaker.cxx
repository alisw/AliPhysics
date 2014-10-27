// $Id: AliEmcalTriggerMaker.cxx 64593 2013-10-18 10:23:58Z loizides $
//
// Class to make array of trigger patch objects in AOD/ESD events.
//
// Author: J.Kral
#include <TClonesArray.h>
#include <TArrayI.h>
#include <THashList.h>
#include "AliAODCaloTrigger.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerTypes.h"
#include "AliEmcalTriggerPatchInfo.h"
#include "AliEmcalTriggerSetupInfo.h"
#include "AliLog.h"
#include "AliVCaloCells.h"
#include "AliVCaloTrigger.h"
#include "AliVVZERO.h"
#include "AliEmcalTriggerMaker.h"

#include "THistManager.h"
#include "TString.h"

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
  fITrigger(0),
  fDoQA(kFALSE),
  fQAHistos(NULL)
{
  // Constructor.
  memset(fThresholdConstants, 0, sizeof(Int_t) * 12);
  memset(fPatchADCSimple, 0, sizeof(Int_t) * kPatchCols * kPatchRows);
  memset(fPatchADC, 0, sizeof(Int_t) * kPatchCols * kPatchRows);
}

//________________________________________________________________________
AliEmcalTriggerMaker::AliEmcalTriggerMaker(const char *name, Bool_t doQA) :
  AliAnalysisTaskEmcal(name,doQA),
  fCaloTriggersOutName("EmcalTriggers"),
  fCaloTriggerSetupOutName("EmcalTriggersSetup"),
  fV0InName("AliAODVZERO"),
  fCaloTriggersOut(0),
  fCaloTriggerSetupOut(0),
  fSimpleOfflineTriggers(0),
  fV0(0),
  fITrigger(0),
  fDoQA(doQA),
  fQAHistos(NULL)
{
  // Constructor.
  memset(fThresholdConstants, 0, sizeof(Int_t) * 12);
  memset(fPatchADCSimple, 0, sizeof(Int_t) * kPatchCols * kPatchRows);
  memset(fPatchADC, 0, sizeof(Int_t) * kPatchCols * kPatchRows);
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
    fV0 = (AliVVZERO*)InputEvent()->FindListObject(fV0InName);
  }

  // container for simple offline trigger processing
  fSimpleOfflineTriggers = new AliAODCaloTrigger();
  fSimpleOfflineTriggers->Allocate(0);
}

//________________________________________________________________________
void AliEmcalTriggerMaker::UserCreateOutputObjects()
{
  // Do basic QA monitoring (if requested)
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if(fDoQA){
    fQAHistos = new THistManager("TriggerQA");

    TString trtypenames[3] = {"EJE", "EGA", "EL0"};
    for(int itype = 0; itype < 3; itype++){
      fQAHistos->CreateTH2(Form("RCPos%s", trtypenames[itype].Data()), Form("Lower edge position of %s patches (col-row)", trtypenames[itype].Data()), 48, -0.5, 47.5, 64, -0.5, 63.5);
      fQAHistos->CreateTH2(Form("EPCentPos%s", trtypenames[itype].Data()), Form("Center position of the %s trigger patches", trtypenames[itype].Data()), 20, -0.8, 0.8, 100., 1., 4.);
      fQAHistos->CreateTH2(Form("PatchADCvsE%s", trtypenames[itype].Data()), Form("Patch ADC value for trigger type %s", trtypenames[itype].Data()), 200, 0., 200, 200, 0., 200);
    }
    fQAHistos->CreateTH1("triggerBitsAll", "Trigger bits for all incoming patches", 64, -0.5, 63.5);
    fQAHistos->CreateTH1("triggerBitsSel", "Trigger bits for reconstructed patches", 64, -0.5, 63.5);
    fOutput->Add(fQAHistos->GetListOfHistograms());
    PostData(1, fOutput);
  }
}

//________________________________________________________________________
Bool_t AliEmcalTriggerMaker::Run() 
{
  // Create and fill the patch array.

  AliEmcalTriggerPatchInfo *trigger, *triggerMainJet, *triggerMainGamma, *triggerMainLevel0;
  AliEmcalTriggerPatchInfo *triggerMainJetSimple, *triggerMainGammaSimple;

  // delete patch array, clear setup object
  fCaloTriggersOut->Delete();
  fCaloTriggerSetupOut->Clean();

  if (!fCaloTriggers) {
    AliError(Form("Calo triggers container %s not available.", fCaloTriggersName.Data()));
    return kTRUE;
  }
  if (!fCaloCells) {
    AliError(Form("Calo cells container %s not available.", fCaloCellsName.Data()));
    return kTRUE;
  }
  if (!fCaloCells) {
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
  if (fCaloTriggers->GetEntries() > 0) {
    // zero the arrays
    memset(fPatchADC, 0, sizeof(Int_t) * kPatchCols * kPatchRows);

    // go throuth the trigger channels
    while (fCaloTriggers->Next()) {
      // get position in global 2x2 tower coordinates
      // A0 left bottom (0,0)
      Int_t globCol=-1, globRow=-1;
      fCaloTriggers->GetPosition(globCol, globRow);
      // for some strange reason some ADC amps are initialized in reconstruction
      // as -1, neglect those 
      Int_t adcAmp=-1;
      fCaloTriggers->GetL1TimeSum(adcAmp);
      if (adcAmp>-1)
	    fPatchADC[globCol][globRow] = adcAmp;
    } // patches
  } // array not empty
  
  // fill the array for offline trigger processing
  // using calibrated cell energies
  memset(fPatchADCSimple, 0, sizeof(Int_t) * kPatchRows * kPatchCols);

  // fill the patch ADCs from cells
  Int_t nCell = fCaloCells->GetNumberOfCells();
  for(Int_t iCell = 0; iCell < nCell; ++iCell) {
    // get the cell info, based in index in array
    Short_t cellId = fCaloCells->GetCellNumber(iCell);
    Double_t amp = fCaloCells->GetAmplitude(iCell);
    // get position
    Int_t absId=-1;
    fGeom->GetFastORIndexFromCellIndex(cellId, absId);
    Int_t globCol=-1, globRow=-1;
    fGeom->GetPositionInEMCALFromAbsFastORIndex(absId, globCol, globRow);
    // add
    fPatchADCSimple[globCol][globRow] += amp/kEMCL1ADCtoGeV;
  }

  // dig out common data (thresholds)
  // 0 - jet high, 1 - gamma high, 2 - jet low, 3 - gamma low
  fCaloTriggerSetupOut->SetThresholds(fCaloTriggers->GetL1Threshold(0),
                                      fCaloTriggers->GetL1Threshold(1),
                                      fCaloTriggers->GetL1Threshold(2),
                                      fCaloTriggers->GetL1Threshold(3));

  // get the V0 value and compute and set the offline thresholds
  // get V0, compute thresholds and save them as global parameters
  Int_t v0[2];
  v0[0] = fV0->GetTriggerChargeA();
  v0[1] = fV0->GetTriggerChargeC();
  ULong64_t v0S = v0[0] + v0[1];
  fSimpleOfflineTriggers->SetL1V0(v0);
  
  for (Int_t i = 0; i < 4; ++i) {
    // A*V0^2/2^32+B*V0/2^16+C
    ULong64_t thresh = ( ((ULong64_t)fThresholdConstants[i][0]) * v0S * v0S ) >> 32;
    thresh += ( ((ULong64_t)fThresholdConstants[i][1]) * v0S ) >> 16;
    thresh += ((ULong64_t)fThresholdConstants[i][2]);
    fSimpleOfflineTriggers->SetL1Threshold(i,thresh);
  }
  
  // save the thresholds in output object
  fCaloTriggerSetupOut->SetThresholdsSimple(fSimpleOfflineTriggers->GetL1Threshold(0),
					    fSimpleOfflineTriggers->GetL1Threshold(1),
					    fSimpleOfflineTriggers->GetL1Threshold(2),
					    fSimpleOfflineTriggers->GetL1Threshold(3));

  // run the trigger
  RunSimpleOfflineTrigger();

  // reset for re-run
  fCaloTriggers->Reset();
  fSimpleOfflineTriggers->Reset();

  // class is not empty
  if (fCaloTriggers->GetEntries() > 0 ||  fSimpleOfflineTriggers->GetEntries() > 0) {
    fITrigger = 0;
    triggerMainGamma = 0;
    triggerMainJet = 0;
    triggerMainGammaSimple = 0;
    triggerMainJetSimple = 0;
    triggerMainLevel0 = 0;

    // go throuth the trigger channels, real first, then offline
    Bool_t isOfflineSimple=0;
    while (NextTrigger(isOfflineSimple)) {
      // process jet
      trigger = ProcessPatch(kTMEMCalJet, isOfflineSimple);
      // save main jet triggers in event
      if (trigger != 0) {
        // check if more energetic than others for main patch marking
        if (!isOfflineSimple) {
          if (triggerMainJet == 0 || (triggerMainJet->GetPatchE() < trigger->GetPatchE()))
            triggerMainJet = trigger;
        } else {
          if (triggerMainJetSimple == 0 || (triggerMainJetSimple->GetPatchE() < trigger->GetPatchE()))
            triggerMainJetSimple = trigger;
        }
      }
      
      // process gamma
      trigger = ProcessPatch(kTMEMCalGamma, isOfflineSimple);
      // save main gamma triggers in event
      if (trigger != 0) {
        // check if more energetic than others for main patch marking
        if (!isOfflineSimple) {
          if (triggerMainGamma == 0 || (triggerMainGamma->GetPatchE() < trigger->GetPatchE()))
            triggerMainGamma = trigger;
        } else {
          if (triggerMainGammaSimple == 0 || (triggerMainGammaSimple->GetPatchE() < trigger->GetPatchE()))
            triggerMainGammaSimple = trigger;
        }
      }

      // level 0 triggers
      trigger = ProcessPatch(kTMEMCalLevel0, isOfflineSimple);
      // save main level0 trigger in the event
      if (trigger) {
        if (!triggerMainLevel0 || (triggerMainLevel0->GetPatchE() < trigger->GetPatchE()))
          triggerMainLevel0 = trigger;
      }
    } // triggers
    
    // mark the most energetic patch as main
    // for real and also simple offline
    if (triggerMainJet != 0) {
      Int_t tBits = triggerMainJet->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << 24 );
      triggerMainJet->SetTriggerBits( tBits );
    }
    if (triggerMainJetSimple != 0) {
      Int_t tBits = triggerMainJetSimple->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << 24 );
      triggerMainJetSimple->SetTriggerBits(tBits);
    }
    if (triggerMainGamma != 0) {
      Int_t tBits = triggerMainGamma->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << 24 );
      triggerMainGamma->SetTriggerBits( tBits );
    }
    if (triggerMainGammaSimple != 0) {
      Int_t tBits = triggerMainGammaSimple->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << 24 );
      triggerMainGammaSimple->SetTriggerBits( tBits );
    }
    if(triggerMainLevel0){
      Int_t tBits = triggerMainLevel0->GetTriggerBits();
      // main trigger flag
      tBits |= (1 << 24);
      triggerMainLevel0->SetTriggerBits(tBits);
    }
  } // there are some triggers

  return kTRUE;
}

//________________________________________________________________________
AliEmcalTriggerPatchInfo* AliEmcalTriggerMaker::ProcessPatch(TriggerMakerTriggerType_t type, Bool_t isOfflineSimple)
{
  // Process and fill trigger patch.
  // check if jet trigger low or high
  Int_t tBits=-1;
  if (!isOfflineSimple)
    fCaloTriggers->GetTriggerBits(tBits);
  else
    fSimpleOfflineTriggers->GetTriggerBits(tBits);

  Int_t nBitsFound = 0;
  Int_t bitsFound[64];
  if(fDoQA){
    for(int ibit = 0; ibit < sizeof(tBits)*8; ibit++) {
      if(tBits & (1 << ibit)){
        bitsFound[nBitsFound++] = ibit;
        fQAHistos->FillTH1("triggerBitsAll", ibit);
      }
    }
  }
	
  if ((type == kTMEMCalJet    && !IsEJE( tBits )) || 
      (type == kTMEMCalGamma  && !IsEGA( tBits )) || 
      (type == kTMEMCalLevel0 && !(CheckForL0(*fCaloTriggers))))
    return 0;
  TString trtypenames[3] = {"EJE", "EGA", "EL0"}; // For QA

  // save primary vertex in vector
  TVector3 vertex;
  vertex.SetXYZ(fVertex[0], fVertex[1], fVertex[2]);

  // get position in global 2x2 tower coordinates
  // A0 left bottom (0,0)
  Int_t globCol=-1, globRow=-1;
  if (!isOfflineSimple)
    fCaloTriggers->GetPosition(globCol,globRow);
  else
    fSimpleOfflineTriggers->GetPosition(globCol, globRow);

  // get the absolute trigger ID
  Int_t absId=-1;
  fGeom->GetAbsFastORIndexFromPositionInEMCAL(globCol, globRow, absId);
  // convert to the 4 absId of the cells composing the trigger channel
  Int_t cellAbsId[4]={-1,-1,-1,-1};
  fGeom->GetCellIndexFromFastORIndex(absId, cellAbsId);
	
  // get low left edge (eta max, phi min)
  TVector3 edge1;
  fGeom->GetGlobal(cellAbsId[0], edge1);
	
  // sum the available energy in the 32/32 window of cells
  // step over trigger channels and get all the corresponding cells
  // make CM
  Double_t amp = 0;
  Int_t cmiCol = 0;
  Int_t cmiRow = 0;
  Int_t adcAmp = 0;
  int nfastor = (type == kTMEMCalJet) ? 16 : 2; // 32x32 cell window for L1 Jet trigger, 4x4 for L1 Gamma or L0 trigger
  for (Int_t i = 0; i < nfastor; ++i) {
    for (Int_t j = 0; j < nfastor; ++j) {
	  // get the 4 cells composing the trigger channel
	  fGeom->GetAbsFastORIndexFromPositionInEMCAL(globCol+i, globRow+j, absId);
	  fGeom->GetCellIndexFromFastORIndex(absId, cellAbsId);
	  // add amplitudes and find patch edges
	  for (Int_t k = 0; k < 4; ++k) {
	    Double_t ca = fCaloCells->GetCellAmplitude(cellAbsId[k]);
	    //fGeom->GetGlobal(cellAbsId[k], cellCoor);
	    amp += ca;
	    cmiCol += ca*(Double_t)i;
	    cmiRow += ca*(Double_t)j;
	  }
	  // add the STU ADCs in the patch (in case of L1) or the TRU Amplitude (in case of L0)
	  if (!isOfflineSimple )
	      if(type == kTMEMCalLevel0){
	        adcAmp += fPatchADC[globCol+i][globRow+j] * 4; // precision loss in case of global integer field
	      } else
	        adcAmp += fPatchADC[globCol+i][globRow+j];
	  else
	    adcAmp += fPatchADCSimple[globCol+i][globRow+j];
    }
  }

  if (amp == 0) {
    AliDebug(2,"EMCal trigger patch with 0 energy.");
    return 0;
  }
  
  // get the CM and patch index
  cmiCol /= amp;
  cmiRow /= amp;
  Int_t cmCol = globCol + (Int_t)cmiCol;
  Int_t cmRow = globRow + (Int_t)cmiRow;

  // get the patch and corresponding cells
  fGeom->GetAbsFastORIndexFromPositionInEMCAL( cmCol, cmRow, absId );
  fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );

  // find which out of the 4 cells is closest to CM and get it's position
  Int_t cmiCellCol = TMath::Nint(cmiCol * 2.);
  Int_t cmiCellRow = TMath::Nint(cmiRow * 2.);
  TVector3 centerMass;
  fGeom->GetGlobal(cellAbsId[(cmiCellRow%2)*2 + cmiCellCol%2], centerMass);
	
  // get up right edge (eta min, phi max)
  // get the absolute trigger ID
  Int_t posOffset=-1;
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
  fGeom->GetAbsFastORIndexFromPositionInEMCAL(globCol+posOffset, globRow+posOffset, absId);
  fGeom->GetCellIndexFromFastORIndex(absId, cellAbsId);
  TVector3 edge2;
  fGeom->GetGlobal(cellAbsId[3], edge2);
	
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
  fGeom->GetAbsFastORIndexFromPositionInEMCAL(globCol+posOffset, globRow+posOffset, absId);
  fGeom->GetCellIndexFromFastORIndex(absId, cellAbsId);
  TVector3 center1;
  fGeom->GetGlobal(cellAbsId[3], center1);
	
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
  fGeom->GetAbsFastORIndexFromPositionInEMCAL(globCol+posOffset, globRow+posOffset, absId);
  fGeom->GetCellIndexFromFastORIndex(absId, cellAbsId);
  TVector3 center2;
  fGeom->GetGlobal(cellAbsId[0], center2);
	
  TVector3 centerGeo(center1);
  centerGeo += center2;
  centerGeo *= 0.5;
	
  // relate all to primary vertex
  centerGeo -= vertex;
  centerMass -= vertex;
  edge1 -= vertex;
  edge2 -= vertex;

  Int_t isMC = MCEvent() ? 1 : 0;
  Int_t offSet = (1 - isMC) * kTriggerTypeEnd;
	
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
    // Explicitly set the level 0 bit to overcome the masking out
    tBits |= 1 << (offSet + kL0);
    tBits = tBits & ~( 1 << (kTriggerTypeEnd + kL1JetLow) | 1 << (kTriggerTypeEnd + kL1JetHigh) | 1 << (kL1JetLow) | 1 << (kL1JetHigh) |
		       1 << (kTriggerTypeEnd + kL1GammaLow) | 1 << (kTriggerTypeEnd + kL1GammaHigh) | 1 << (kL1GammaLow) | 1 << (kL1GammaHigh));
    break;
  };

  // save the trigger object
  AliEmcalTriggerPatchInfo *trigger = 
    new ((*fCaloTriggersOut)[fITrigger]) AliEmcalTriggerPatchInfo();
  fITrigger++;
  trigger->SetCenterGeo(centerGeo, amp);
  trigger->SetCenterMass(centerMass, amp);
  trigger->SetEdge1(edge1, amp);
  trigger->SetEdge2(edge2, amp);
  trigger->SetADCAmp(adcAmp);
  trigger->SetTriggerBits(tBits);
  trigger->SetOffSet(offSet);
  trigger->SetEdgeCell(globCol*2, globRow*2); // from triggers to cells
  if(fDoQA){
    fQAHistos->FillTH2(Form("RCPos%s", trtypenames[type].Data()), globCol, globRow);
    fQAHistos->FillTH2(Form("EPCentPos%s", trtypenames[type].Data()), centerGeo.Eta(), centerGeo.Phi());
    fQAHistos->FillTH2(Form("PatchADCvsE%s", trtypenames[type].Data()), adcAmp, trigger->GetPatchE());
    if(nBitsFound){
      for(int ibit = 0; ibit < nBitsFound; ibit++)
        fQAHistos->FillTH1("triggerBitsSel", ibit);
    }
  }
  return trigger;
}

//________________________________________________________________________
void AliEmcalTriggerMaker::RunSimpleOfflineTrigger() 
{
  // Runs a simple offline trigger algorithm.
  // It creates separate patches for jet and gamma triggers
  // on the same positions (different from STU reconstruction behavior)
  // TODO:: change to merge

  TArrayI tBitsArray, rowArray, colArray;
  
  // 0 thresholds = no processing
  if (fCaloTriggerSetupOut->GetThresholdJetLowSimple() == 0 &&
      fCaloTriggerSetupOut->GetThresholdJetHighSimple() == 0 )
    return;
  
  // run the trigger algo, stepping by 8 towers (= 4 trigger channels)
  for (Int_t i = 0; i < 32; i += 4) {
    for (Int_t j = 0; j < 48; j += 4) {
      Int_t tSum  = 0;
      Int_t tBits = 0;
      // window
      for (Int_t k = 0; k < 16; ++k)
        for (Int_t l = 0; l < 16; ++l)
          tSum += (ULong64_t)fPatchADCSimple[i+k][j+l];
      
      // check thresholds
      if (tSum > fCaloTriggerSetupOut->GetThresholdJetLowSimple())
        tBits = tBits | ( 1 << ( kTriggerTypeEnd + kL1JetLow ));
      if (tSum > fCaloTriggerSetupOut->GetThresholdJetHighSimple())
        tBits = tBits | ( 1 << ( kTriggerTypeEnd + kL1JetHigh ));
      
      // add trigger values
      if (tBits != 0) {
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
  for (Int_t i = 0; i < 46; ++i) {
    for (Int_t j = 0; j < 62; ++j) {
      Int_t tSum = 0;
      Int_t tBits = 0;
      
      // window
      for (Int_t k = 0; k < 2; ++k)
        for (Int_t l = 0; l < 2; ++l)
          tSum += (ULong64_t)fPatchADCSimple[i+k][j+l];
      
      // check thresholds
      if (tSum > fCaloTriggerSetupOut->GetThresholdGammaLowSimple())
        tBits = tBits | ( 1 << ( kTriggerTypeEnd + kL1GammaLow ));
      if (tSum > fCaloTriggerSetupOut->GetThresholdGammaHighSimple())
        tBits = tBits | ( 1 << ( kTriggerTypeEnd + kL1GammaHigh ));
      
      // add trigger values
      if (tBits != 0) {
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
  fSimpleOfflineTriggers->Allocate(tBitsArray.GetSize());
  for (Int_t i = 0; i < tBitsArray.GetSize(); ++i){
    fSimpleOfflineTriggers->Add(colArray[i],rowArray[i], 0, 0, 0, 0, 0, tBitsArray[i]);
  }
}

//________________________________________________________________________
Bool_t AliEmcalTriggerMaker::NextTrigger(Bool_t &isOfflineSimple) 
{
  // Get next trigger
  
  isOfflineSimple = kFALSE;
  Bool_t loopContinue = fCaloTriggers->Next();
  if (!loopContinue) {
    loopContinue = fSimpleOfflineTriggers->Next();
    isOfflineSimple = kTRUE;
  }
  return loopContinue;
}

//________________________________________________________________________
Bool_t AliEmcalTriggerMaker::CheckForL0(const AliVCaloTrigger& trg) const {
  // Check whether the patch is a level0 patch
  if(MCEvent()){
    // For Monte-Carlo select
    Int_t tbits(-1);
    trg.GetTriggerBits(tbits);
    return tbits & (1 << kL0);
  } else {
    // For Data check from the level0 times if the trigger has fired at level0
    Int_t nl0times(0);
    Bool_t l0fired(kFALSE);
    trg.GetNL0Times(nl0times);
    if(nl0times){
      TArrayI l0times(nl0times);
      trg.GetL0Times(l0times.GetArray());
      // Apply timing cut to see if a L0 has fired
      for(Int_t *l0timeIter = l0times.GetArray(); l0timeIter < l0times.GetArray() + l0times.GetSize(); l0timeIter++){
        if(*l0timeIter > 7 && *l0timeIter < 10){
          l0fired = kTRUE;
          break;
        }
      }
    }
    return l0fired;
  }
}
