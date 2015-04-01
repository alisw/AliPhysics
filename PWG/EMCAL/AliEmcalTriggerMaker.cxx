/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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
#include <TClonesArray.h>
#include <TArrayI.h>
#include <THashList.h>
#include "AliAODCaloTrigger.h"
#include "AliEMCALGeometry.h"
//#include "AliEMCALTriggerTypes.h"
#include "AliEmcalTriggerPatchInfo.h"
#include "AliEmcalTriggerSetupInfo.h"
#include "AliLog.h"
#include "AliVCaloCells.h"
#include "AliVCaloTrigger.h"
#include "AliVVZERO.h"
#include "AliEmcalTriggerMaker.h"

#include "THistManager.h"
#include "TString.h"

#include <bitset>

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerMaker)
ClassImp(AliEmcalTriggerMaker::AliEmcalTriggerChannelPosition)
ClassImp(AliEmcalTriggerMaker::AliEmcalTriggerChannelContainer)
/// \endcond

using namespace std;

const TString AliEmcalTriggerMaker::fgkTriggerTypeNames[5] = {"EJE", "EGA", "EL0", "REJE", "REGA"};

/**
 * Dummy constructor
 */
AliEmcalTriggerMaker::AliEmcalTriggerMaker() : 
  AliAnalysisTaskEmcal("AliEmcalTriggerMaker",kFALSE),
  fCaloTriggersOutName("EmcalTriggers"),
  fCaloTriggerSetupOutName("EmcalTriggersSetup"),
  fV0InName("AliAODVZERO"),
  fUseTriggerBitConfig(kNewConfig),
  fTriggerBitConfig(NULL),
  fCaloTriggersOut(0),
  fCaloTriggerSetupOut(0),
  fSimpleOfflineTriggers(0),
  fV0(0),
  fITrigger(0),
  fRejectOffAcceptancePatches(kFALSE),
  fDoQA(kFALSE),
  fQAHistos(NULL)
{
  fRunTriggerType[kTMEMCalJet] = kTRUE;
  fRunTriggerType[kTMEMCalGamma] = kTRUE;
  fRunTriggerType[kTMEMCalLevel0] = kTRUE;
  fRunTriggerType[kTMEMCalRecalcJet] = kTRUE;
  fRunTriggerType[kTMEMCalRecalcGamma] = kTRUE;
  memset(fThresholdConstants, 0, sizeof(Int_t) * 12);
  memset(fPatchADCSimple, 0, sizeof(Double_t) * kPatchCols * kPatchRows);
  memset(fPatchADC, 0, sizeof(Int_t) * kPatchCols * kPatchRows);
  memset(fPatchAmplitudes, 0, sizeof(Float_t) * kPatchCols * kPatchRows);
  memset(fLevel0TimeMap, 0, sizeof(Char_t) * kPatchCols * kPatchRows);
}

/**
 * Named constructor.
 * \param name Name of the trigger maker task
 * \param doQA Switch on/off QA
 */
AliEmcalTriggerMaker::AliEmcalTriggerMaker(const char *name, Bool_t doQA) :
  AliAnalysisTaskEmcal(name,doQA),
  fCaloTriggersOutName("EmcalTriggers"),
  fCaloTriggerSetupOutName("EmcalTriggersSetup"),
  fV0InName("AliAODVZERO"),
  fUseTriggerBitConfig(kNewConfig),
  fTriggerBitConfig(NULL),
  fCaloTriggersOut(0),
  fCaloTriggerSetupOut(0),
  fSimpleOfflineTriggers(0),
  fV0(0),
  fITrigger(0),
  fRejectOffAcceptancePatches(kFALSE),
  fDoQA(doQA),
  fQAHistos(NULL)
{
  // Constructor.
  fRunTriggerType[kTMEMCalJet] = kTRUE;
  fRunTriggerType[kTMEMCalGamma] = kTRUE;
  fRunTriggerType[kTMEMCalLevel0] = kTRUE;
  fRunTriggerType[kTMEMCalRecalcJet] = kTRUE;
  fRunTriggerType[kTMEMCalRecalcGamma] = kTRUE;
  memset(fThresholdConstants, 0, sizeof(Int_t) * 12);
  memset(fPatchAmplitudes, 0, sizeof(Float_t) * kPatchCols * kPatchRows);
  memset(fPatchADCSimple, 0, sizeof(Double_t) * kPatchCols * kPatchRows);
  memset(fPatchADC, 0, sizeof(Int_t) * kPatchCols * kPatchRows);
  memset(fLevel0TimeMap, 0, sizeof(Char_t) * kPatchCols * kPatchRows);
}

/**
 * Destructor
 */
AliEmcalTriggerMaker::~AliEmcalTriggerMaker()
{
  if(fTriggerBitConfig) delete fTriggerBitConfig;
}

/**
 * Init the analysis.
 */
void AliEmcalTriggerMaker::ExecOnce()
{
  AliAnalysisTaskEmcal::ExecOnce();

  if (!fInitialized)
    return;

  if(!fTriggerBitConfig){
    switch(fUseTriggerBitConfig){
    case kNewConfig:
      fTriggerBitConfig = new AliEmcalTriggerBitConfigNew();
      break;
    case kOldConfig:
      fTriggerBitConfig = new AliEmcalTriggerBitConfigOld();
      break;
    }
  }

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

/**
 * Do basic QA monitoring (if requested)
 */
void AliEmcalTriggerMaker::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if(fDoQA && fOutput){
    fQAHistos = new THistManager("TriggerQA");
    const char *patchtypes[2] = {"Online", "Offline"};

    for(int itype = 0; itype < 5; itype++){
      for(const char **patchtype = patchtypes; patchtype < patchtypes + 2; ++patchtype){
        fQAHistos->CreateTH2(Form("RCPos%s%s", fgkTriggerTypeNames[itype].Data(), *patchtype), Form("Lower edge position of %s %s patches (col-row);iEta;iPhi", *patchtype, fgkTriggerTypeNames[itype].Data()), 48, -0.5, 47.5, 64, -0.5, 63.5);
        fQAHistos->CreateTH2(Form("EPCentPos%s%s", fgkTriggerTypeNames[itype].Data(), *patchtype), Form("Center position of the %s %s trigger patches;#eta;#phi", *patchtype, fgkTriggerTypeNames[itype].Data()), 20, -0.8, 0.8, 100, 1., 4.);
        fQAHistos->CreateTH2(Form("PatchADCvsE%s%s", fgkTriggerTypeNames[itype].Data(), *patchtype), Form("Patch ADC value for trigger type %s %s;Trigger ADC;FEE patch energy (GeV)", *patchtype, fgkTriggerTypeNames[itype].Data()), 2000, 0., 2000, 200, 0., 200);
      }
    }
    fQAHistos->CreateTH1("triggerBitsAll", "Trigger bits for all incoming patches;bit nr", 64, -0.5, 63.5);
    fQAHistos->CreateTH1("triggerBitsSel", "Trigger bits for reconstructed patches;bit nr", 64, -0.5, 63.5);
    fOutput->Add(fQAHistos->GetListOfHistograms());
    PostData(1, fOutput);
  }
}

/**
 * Create and fill the patch array.
 * \return Always true.
 */
Bool_t AliEmcalTriggerMaker::Run() 
{
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

  // zero the arrays
  memset(fPatchADC, 0, sizeof(Int_t) * kPatchCols * kPatchRows);
  memset(fPatchAmplitudes, 0, sizeof(Float_t) * kPatchCols * kPatchRows);
  memset(fLevel0TimeMap, 0, sizeof(Char_t) * kPatchCols * kPatchRows);

  // first run over the patch array to compose a map of 2x2 patch energies
  // which is then needed to construct the full patch ADC energy
  // class is not empty
  if (fCaloTriggers->GetEntries() > 0) {

    // go throuth the trigger channels
    while (fCaloTriggers->Next()) {
      // get position in global 2x2 tower coordinates
      // A0 left bottom (0,0)
      Int_t globCol=-1, globRow=-1;
      fCaloTriggers->GetPosition(globCol, globRow);
      // exclude channel completely if it is masked as hot channel
      if(fBadChannels.HasChannel(globCol, globRow)) continue;
      // for some strange reason some ADC amps are initialized in reconstruction
      // as -1, neglect those 
      Int_t adcAmp=-1;
      fCaloTriggers->GetL1TimeSum(adcAmp);
      if (adcAmp>-1)
	    fPatchADC[globCol][globRow] = adcAmp;

      // Handling for L0 triggers
      // For the ADC value we use fCaloTriggers->GetAmplitude()
      // In data, all patches which have 4 TRUs with proper level0 times are
      // valid trigger patches. Therefore we need to check all neighbors for
      // the level0 times, not only the bottom left. In order to obtain this
      // information, a lookup table with the L0 times for each TRU is created
      Float_t amplitude(0);
      fCaloTriggers->GetAmplitude(amplitude);
      if(amplitude < 0) amplitude = 0;
      fPatchAmplitudes[globCol][globRow] = amplitude;
      Int_t nl0times(0);
      fCaloTriggers->GetNL0Times(nl0times);
      if(nl0times){
        TArrayI l0times(nl0times);
        fCaloTriggers->GetL0Times(l0times.GetArray());
        for(int itime = 0; itime < nl0times; itime++){
          if(l0times[itime] >7 && l0times[itime] < 10){
            fLevel0TimeMap[globCol][globRow] = static_cast<Char_t>(l0times[itime]);
            break;
          }
        }
      }
    } // patches
  } // array not empty
  
  // fill the array for offline trigger processing
  // using calibrated cell energies
  memset(fPatchADCSimple, 0, sizeof(Double_t) * kPatchRows * kPatchCols); // This works, but in principle cannot assume that the representation for 0 in a double really consists of only zeros...

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
       if(fRunTriggerType[kTMEMCalJet]){
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
      }
      
      // process gamma
      if(fRunTriggerType[kTMEMCalGamma]){
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
      }

      // level 0 triggers
      if(fRunTriggerType[kTMEMCalLevel0]){
        trigger = ProcessPatch(kTMEMCalLevel0, isOfflineSimple);
        // save main level0 trigger in the event
        if (trigger) {
          if (!triggerMainLevel0 || (triggerMainLevel0->GetPatchE() < trigger->GetPatchE()))
            triggerMainLevel0 = trigger;
        }
      }

      // Recalculated triggers (max patches without threshold)
      if(fRunTriggerType[kTMEMCalRecalcJet])
	ProcessPatch(kTMEMCalRecalcJet, isOfflineSimple);
      if(fRunTriggerType[kTMEMCalRecalcGamma])
	ProcessPatch(kTMEMCalRecalcGamma, isOfflineSimple);
    } // triggers
    
    // mark the most energetic patch as main
    // for real and also simple offline
    if (triggerMainJet != 0) {
      Int_t tBits = triggerMainJet->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << AliEmcalTriggerPatchInfo::kMainTriggerBitNum );
      triggerMainJet->SetTriggerBits( tBits );
    }
    if (triggerMainJetSimple != 0) {
      Int_t tBits = triggerMainJetSimple->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << AliEmcalTriggerPatchInfo::kMainTriggerBitNum );
      triggerMainJetSimple->SetTriggerBits(tBits);
    }
    if (triggerMainGamma != 0) {
      Int_t tBits = triggerMainGamma->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << AliEmcalTriggerPatchInfo::kMainTriggerBitNum );
      triggerMainGamma->SetTriggerBits( tBits );
    }
    if (triggerMainGammaSimple != 0) {
      Int_t tBits = triggerMainGammaSimple->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << AliEmcalTriggerPatchInfo::kMainTriggerBitNum );
      triggerMainGammaSimple->SetTriggerBits( tBits );
    }
    if(triggerMainLevel0){
      Int_t tBits = triggerMainLevel0->GetTriggerBits();
      // main trigger flag
      tBits |= (1 << AliEmcalTriggerPatchInfo::kMainTriggerBitNum);
      triggerMainLevel0->SetTriggerBits(tBits);
    }
  } // there are some triggers

  // Diagnostics
  int npatchOnline = 0;
  for(TIter patchIter = TIter(fCaloTriggersOut).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEmcalTriggerPatchInfo *mypatch = static_cast<AliEmcalTriggerPatchInfo *>(*patchIter);
    if(mypatch->IsOfflineSimple()) continue;
    AliDebug(1,Form("Patch with bits: %s, types: JH[%s], JL[%s], GH[%s], GL[%s], L0[%s]",
        std::bitset<sizeof(int)*4>(mypatch->GetTriggerBits()).to_string().c_str(),
        (mypatch->IsJetHigh() ? "y" : "n"), (mypatch->IsJetLow() ? "y" : "n"),
        (mypatch->IsGammaHigh() ? "y" : "n"), (mypatch->IsGammaLow() ? "y" : "n"),(mypatch->IsLevel0() ? "y" : "n")));
    npatchOnline++;
  }
  AliDebug(1, Form("Number of online patches: %d", npatchOnline));

  return kTRUE;
}

/**
 * Process and fill trigger patch.
 * check if jet trigger low or high
 * \param type Type of the patch (Jet, gamma, Level0)
 * \param isOfflineSimple Switch between online and offline patches
 * \return The new patch (NULL in case of failure)
 */
AliEmcalTriggerPatchInfo* AliEmcalTriggerMaker::ProcessPatch(TriggerMakerTriggerType_t type, Bool_t isOfflineSimple)
{
  Int_t tBits=-1;
  if (!isOfflineSimple)
    fCaloTriggers->GetTriggerBits(tBits);
  else
    fSimpleOfflineTriggers->GetTriggerBits(tBits);

  if(fDoQA){
    for(unsigned int ibit = 0; ibit < sizeof(tBits)*8; ibit++) {
      if(tBits & (1 << ibit)){
        fQAHistos->FillTH1("triggerBitsAll", ibit);
      }
    }
  }
	
  if ((type == kTMEMCalJet    && !IsEJE( tBits )) || 
      (type == kTMEMCalGamma  && !IsEGA( tBits )) || 
      (type == kTMEMCalLevel0 && !(CheckForL0(*fCaloTriggers))) ||
      (type == kTMEMCalRecalcJet && (tBits & (1 << AliEmcalTriggerPatchInfo::kRecalcJetBitNum))==0) ||
      (type == kTMEMCalRecalcGamma && (tBits & (1 << AliEmcalTriggerPatchInfo::kRecalcGammaBitNum))==0) ) 
    return 0;

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

  if(globCol >= kPatchCols || globRow >= kPatchRows){
  }

  if(fRejectOffAcceptancePatches){
    int patchsize = 2;
    if(type == kTMEMCalJet || type == kTMEMCalRecalcJet) patchsize = 16;
    if((globCol + patchsize >= kPatchCols) || (globCol + patchsize >= kPatchRows)){
      AliError(Form("Invalid patch position for patch type %s: Col[%d], Row[%d] - patch rejected", fgkTriggerTypeNames[type].Data(), globCol, globRow));
      return NULL;
    }
  }

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
  Float_t amp = 0;
  Float_t cmiCol = 0;
  Float_t cmiRow = 0;
  Int_t adcAmp = 0;
  Double_t adcOfflineAmp = 0;
  int nfastor = (type == kTMEMCalJet || type == kTMEMCalRecalcJet) ? 16 : 2; // 32x32 cell window for L1 Jet trigger, 4x4 for L1 Gamma or L0 trigger
  for (Int_t i = 0; i < nfastor; ++i) {
    for (Int_t j = 0; j < nfastor; ++j) {
	  // get the 4 cells composing the trigger channel
	  fGeom->GetAbsFastORIndexFromPositionInEMCAL(globCol+i, globRow+j, absId);
	  fGeom->GetCellIndexFromFastORIndex(absId, cellAbsId);
	  // add amplitudes and find patch edges
	  for (Int_t k = 0; k < 4; ++k) {
	    Float_t ca = fCaloCells->GetCellAmplitude(cellAbsId[k]);
	    //fGeom->GetGlobal(cellAbsId[k], cellCoor);
	    amp += ca;
	    cmiCol += ca*(Float_t)i;
	    cmiRow += ca*(Float_t)j;
	  }
	  // add the STU ADCs in the patch (in case of L1) or the TRU Amplitude (in case of L0)
	  if(type == kTMEMCalLevel0){
	    adcAmp += static_cast<Int_t>(fPatchAmplitudes[globCol+i][globRow+j] * 4); // precision loss in case of global integer field
	  } else {
	    adcAmp += fPatchADC[globCol+i][globRow+j];
	  }

	  adcOfflineAmp += fPatchADCSimple[globCol+i][globRow+j];
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
  case kTMEMCalRecalcJet:
    posOffset = 15;
    break;
  case kTMEMCalGamma:
  case kTMEMCalRecalcGamma:
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
  case kTMEMCalRecalcJet:
    posOffset = 7;
    break;
  case kTMEMCalGamma:
  case kTMEMCalRecalcGamma:
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
  case kTMEMCalRecalcJet:
    posOffset = 8;
    break;
  case kTMEMCalGamma:
  case kTMEMCalRecalcGamma:
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
  Int_t offSet = (1 - isMC) * fTriggerBitConfig->GetTriggerTypesEnd();
	
  // fix tbits .. remove the unwanted type triggers
  // for Jet and Gamma triggers we remove also the level 0 bit since it will be stored in the level 0 patch
  // for level 0 we remove all gamma and jet trigger bits
  switch(type){
  case kTMEMCalJet:
    tBits = tBits & ~( 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetGammaLowBit()) | 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetGammaHighBit()) |
        1 << (fTriggerBitConfig->GetGammaLowBit()) | 1 << (fTriggerBitConfig->GetGammaHighBit()) |
		       1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetLevel0Bit()) | 1 << (fTriggerBitConfig->GetLevel0Bit()));
    break;
  case kTMEMCalGamma:
    tBits = tBits & ~( 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetJetLowBit()) | 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetJetHighBit()) |
        1 << (fTriggerBitConfig->GetJetLowBit()) | 1 << (fTriggerBitConfig->GetJetHighBit()) |
		       1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetLevel0Bit()) | 1 << (fTriggerBitConfig->GetLevel0Bit()));
    break;
  case kTMEMCalLevel0:
    // Explicitly set the level 0 bit to overcome the masking out
    tBits |= 1 << (offSet + fTriggerBitConfig->GetLevel0Bit());
    tBits = tBits & ~( 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetJetLowBit()) | 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetJetHighBit()) |
        1 << (fTriggerBitConfig->GetJetLowBit()) | 1 << (fTriggerBitConfig->GetJetHighBit()) | 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetGammaLowBit()) |
        1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetGammaHighBit()) | 1 << (fTriggerBitConfig->GetGammaLowBit()) | 1 << (fTriggerBitConfig->GetGammaHighBit()));
    break;
  default:  // recalculated patches don't need any action
    break;
  };

  // save the trigger object
  AliEmcalTriggerPatchInfo *trigger = 
    new ((*fCaloTriggersOut)[fITrigger]) AliEmcalTriggerPatchInfo();
  fITrigger++;
  trigger->SetTriggerBitConfig(fTriggerBitConfig);
  trigger->SetCenterGeo(centerGeo, amp);
  trigger->SetCenterMass(centerMass, amp);
  trigger->SetEdge1(edge1, amp);
  trigger->SetEdge2(edge2, amp);
  trigger->SetADCAmp(adcAmp);
  trigger->SetADCOfflineAmp(Int_t(adcOfflineAmp));
  trigger->SetTriggerBits(tBits);
  trigger->SetOffSet(offSet);
  trigger->SetEdgeCell(globCol*2, globRow*2); // from triggers to cells
  //if(isOfflineSimple)trigger->SetOfflineSimple();
  if(fDoQA){
    TString patchtype = isOfflineSimple ? "Offline" : "Online";
    fQAHistos->FillTH2(Form("RCPos%s%s", fgkTriggerTypeNames[type].Data(), patchtype.Data()), globCol, globRow);
    fQAHistos->FillTH2(Form("EPCentPos%s%s", fgkTriggerTypeNames[type].Data(), patchtype.Data()), centerGeo.Eta(), centerGeo.Phi());
    fQAHistos->FillTH2(Form("PatchADCvsE%s%s", fgkTriggerTypeNames[type].Data(), patchtype.Data()), isOfflineSimple ? adcOfflineAmp : adcAmp, trigger->GetPatchE());
    // Redo checking of found trigger bits after masking of unwanted triggers
    for(unsigned int ibit = 0; ibit < sizeof(tBits)*8; ibit++) {
      if(tBits & (1 << ibit)){
        fQAHistos->FillTH1("triggerBitsSel", ibit);
      }
    }
  }
  return trigger;
}

/**
 * Runs a simple algorithm to calculate patch energies based on
 * the offline/FEE ADC values (useOffline = kTRUE) or
 * the online/trigger values (useOffline = kFALSE.
 *
 *  It creates separate patches for jet and gamma triggers
 *  on the same positions (different from STU reconstruction behavior)
 */
void AliEmcalTriggerMaker::RunSimpleOfflineTrigger() 
{

  TArrayI tBitsArray(4), rowArray(4), colArray(4);
  
  // First entries are for recalculated patches

  tBitsArray[0] = 1 << AliEmcalTriggerPatchInfo::kRecalcJetBitNum;
  colArray[0] = -1;
  rowArray[0] = -1;

  tBitsArray[1] = 1 << AliEmcalTriggerPatchInfo::kRecalcJetBitNum | 1 << AliEmcalTriggerPatchInfo::kSimpleOfflineBitNum;
  colArray[1] = -1;
  rowArray[1] = -1;

  tBitsArray[2] = 1 << AliEmcalTriggerPatchInfo::kRecalcGammaBitNum;
  colArray[2] = -1;
  rowArray[2] = -1;

  tBitsArray[3] = 1 << AliEmcalTriggerPatchInfo::kRecalcGammaBitNum | 1 << AliEmcalTriggerPatchInfo::kSimpleOfflineBitNum;
  colArray[3] = -1;
  rowArray[3] = -1;

  Double_t maxPatchADCoffline = -1;
  Int_t maxPatchADC = -1;
  // run the trigger algo, stepping by 8 towers (= 4 trigger channels)
  Int_t maxCol = 48;
  Int_t maxRow = fGeom->GetNTotalTRU()*2;
  Int_t isMC = MCEvent() ? 1 : 0;
  Int_t bitOffSet = (1 - isMC) * fTriggerBitConfig->GetTriggerTypesEnd();
  for (Int_t i = 0; i <= (maxCol-16); i += 4) {
    for (Int_t j = 0; j <= (maxRow-16); j += 4) {
      Double_t tSumOffline  = 0;
      Int_t tSum  = 0;
      Int_t tBits = 0;
      // window
      for (Int_t k = 0; k < 16; ++k) {
        for (Int_t l = 0; l < 16; ++l) {
	  tSumOffline += fPatchADCSimple[i+k][j+l];	  
	  tSum += (ULong64_t)fPatchADC[i+k][j+l];
	}
      }

      if (tSum > maxPatchADC) { // Mark highest Jet patch
	maxPatchADC = tSum;
	colArray[0] = i;
	rowArray[0] = j;
      }

      if (tSumOffline > maxPatchADCoffline) { // Mark highest Jet patch
	maxPatchADCoffline = tSumOffline;
	colArray[1] = i;
	rowArray[1] = j;
      }

      // check thresholds
      if (tSumOffline > fCaloTriggerSetupOut->GetThresholdJetLowSimple())
	tBits = tBits | ( 1 << ( bitOffSet + fTriggerBitConfig->GetJetLowBit() ));
      if (tSumOffline > fCaloTriggerSetupOut->GetThresholdJetHighSimple())
	tBits = tBits | ( 1 << ( bitOffSet + fTriggerBitConfig->GetJetHighBit() ));
      
      // add trigger values
      if (tBits != 0) {
	// add offline bit
	tBits = tBits | ( 1 << AliEmcalTriggerPatchInfo::kSimpleOfflineBitNum );
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
  maxPatchADC = -1;
  maxPatchADCoffline = -1;

  for (Int_t i = 0; i <= (maxCol-2); ++i) {
    for (Int_t j = 0; j <= (maxRow-2); ++j) {
      Int_t tSum = 0;
      Double_t tSumOffline = 0;
      Int_t tBits = 0;
      
      // window
      for (Int_t k = 0; k < 2; ++k) {
        for (Int_t l = 0; l < 2; ++l) {
	  tSumOffline += fPatchADCSimple[i+k][j+l];
	  tSum += (ULong64_t)fPatchADC[i+k][j+l];
	}
      }

      if (tSum > maxPatchADC) { // Mark highest Gamma patch
	maxPatchADC = tSum;
	colArray[2] = i;
	rowArray[2] = j;
      }
      if (tSumOffline > maxPatchADCoffline) { // Mark highest Gamma patch
	maxPatchADCoffline = tSumOffline;
	colArray[3] = i;
	rowArray[3] = j;
      }

      // check thresholds
      if (tSumOffline > fCaloTriggerSetupOut->GetThresholdGammaLowSimple())
	tBits = tBits | ( 1 << ( bitOffSet + fTriggerBitConfig->GetGammaLowBit() ));
      if (tSumOffline > fCaloTriggerSetupOut->GetThresholdGammaHighSimple())
	tBits = tBits | ( 1 << ( bitOffSet + fTriggerBitConfig->GetGammaHighBit() ));
      
      // add trigger values
      if (tBits != 0) {
	// add offline bit
	tBits = tBits | ( 1 << AliEmcalTriggerPatchInfo::kSimpleOfflineBitNum );
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

/**
 * Get next trigger. Forwards the pointer of the trigger object inside the trigger maker
 * \param isOfflineSimple Switch between online and ofline patches
 * \return True if successful, false otherwise
 */
Bool_t AliEmcalTriggerMaker::NextTrigger(Bool_t &isOfflineSimple) 
{
  
  isOfflineSimple = kFALSE;
  Bool_t loopContinue = fCaloTriggers->Next();
  if (!loopContinue) {
    loopContinue = fSimpleOfflineTriggers->Next();
    isOfflineSimple = kTRUE;
  }
  return loopContinue;
}

/**
 * Accept trigger patch as Level0 patch. Level0 patches are identified as 2x2 FASTOR patches
 * in the same TRU
 * \param trg Triggers object with the pointer set to the patch to inspect
 * \return True if the patch is accepted, false otherwise.
 */
Bool_t AliEmcalTriggerMaker::CheckForL0(const AliVCaloTrigger& trg) const {
  Int_t row, col;trg.GetPosition(col, row);
  Int_t truref(-1), trumod(-1), absFastor(-1), adc(-1);
  fGeom->GetAbsFastORIndexFromPositionInEMCAL(col, row, absFastor);
  fGeom->GetTRUFromAbsFastORIndex(absFastor, truref, adc);
  int nvalid(0);
  for(int ipos = 0; ipos < 2; ipos++){
    if(row + ipos >= kPatchRows) continue;    // boundary check
    for(int jpos = 0; jpos < 2; jpos++){
      if(col + jpos >= kPatchCols) continue;  // boundary check
      // Check whether we are in the same TRU
      trumod = -1;
      fGeom->GetAbsFastORIndexFromPositionInEMCAL(col+jpos, row+ipos, absFastor);
      fGeom->GetTRUFromAbsFastORIndex(absFastor, trumod, adc);
      if(trumod != truref) continue;
      Char_t l0times = fLevel0TimeMap[col + jpos][row + ipos];
      if(l0times > 7 && l0times < 10) nvalid++;
    }
  }
  if (nvalid != 4) return false;
  return true;
}

/**
 * Add a new channel with the postion in column and row to the container, In case the channel
 * is already listed in the trigger channel container we don't add it again.
 * \param col Column of the channel
 * \param row Row of the channel
 */
void AliEmcalTriggerMaker::AliEmcalTriggerChannelContainer::AddChannel(int col, int row){
  if(HasChannel(col, row)) return;
  fChannels.Add(new AliEmcalTriggerChannelPosition(col, row));
}

/**
 * Check whether channel with the position (col, row) is listed in the trigger channel container
 * \param col Column of the channel
 * \param row Row of the channel
 * \return True if the channel is listed, false otherwise
 */
Bool_t AliEmcalTriggerMaker::AliEmcalTriggerChannelContainer::HasChannel(int col, int row){
  AliEmcalTriggerChannelPosition refChannel;
  if(fChannels.FindObject(&refChannel)) return true;
  return false;
}
