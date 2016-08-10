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
#include "AliEmcalTriggerMaker.h"
#include "AliEmcalTriggerSetupInfo.h"
#include "AliEMCALTriggerConstants.h"
#include "AliEMCALTriggerDataGrid.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliLog.h"
#include "AliVCaloCells.h"
#include "AliVCaloTrigger.h"
#include "AliVVZERO.h"
#include "THistManager.h"
#include "TString.h"

#include <array>
#include <bitset>
#include <iostream>
#include <string>

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerMaker)

/// \endcond

using namespace std;

const int AliEmcalTriggerMaker::kColsEta = 48;

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
  fPatchAmplitudes(NULL),
  fPatchADCSimple(NULL),
  fPatchADC(NULL),
  fLevel0TimeMap(NULL),
  fITrigger(0),
  fDoQA(kFALSE),
  fRejectOffAcceptancePatches(kFALSE),
  fQAHistos(NULL),
  fDebugLevel(0)
{
  memset(fThresholdConstants, 0, sizeof(Int_t) * 12);
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
  fPatchAmplitudes(NULL),
  fPatchADCSimple(NULL),
  fPatchADC(NULL),
  fLevel0TimeMap(NULL),
  fITrigger(0),
  fDoQA(doQA),
  fRejectOffAcceptancePatches(kFALSE),
  fQAHistos(NULL),
  fDebugLevel(0)
{
  // Constructor.
  memset(fThresholdConstants, 0, sizeof(Int_t) * 12);
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
      fTriggerBitConfig = new AliEMCALTriggerBitConfigNew();
      break;
    case kOldConfig:
      fTriggerBitConfig = new AliEMCALTriggerBitConfigOld();
      break;
    }
  }

  if (!fCaloTriggersOutName.IsNull()) {
    fCaloTriggersOut = new TClonesArray("AliEMCALTriggerPatchInfo");
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


  // Allocate containers for the ADC values
  int nrows = fGeom->GetNTotalTRU() * 2;
  std::cout << "Allocating channel grid with 48 columns in eta and " << nrows << " rows in phi" << std::endl;
  fPatchAmplitudes->Allocate(48, nrows);
  fPatchADC->Allocate(48, nrows);
  fPatchADCSimple->Allocate(48, nrows);
  fLevel0TimeMap->Allocate(48, nrows);

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

  // Create data structure for energy measurements;
  fPatchAmplitudes = new AliEMCALTriggerDataGrid<float>;
  fPatchADC =  new AliEMCALTriggerDataGrid<int>;
  fPatchADCSimple = new AliEMCALTriggerDataGrid<double>;
  fLevel0TimeMap = new AliEMCALTriggerDataGrid<char>;

  if(fDoQA && fOutput){
    fQAHistos = new THistManager("TriggerQA");
    std::array<std::string, 3> patchtypes = {"Online", "Offline", "Recalc"};

    for(int itype = 0; itype < 5; itype++){
      for(const auto & patchtype : patchtypes){
        fQAHistos->CreateTH2(Form("RCPos%s%s", fgkTriggerTypeNames[itype].Data(), patchtype.c_str()), Form("Lower edge position of %s %s patches (col-row);iEta;iPhi", patchtype.c_str(), fgkTriggerTypeNames[itype].Data()), 48, -0.5, 47.5, 104, -0.5, 103.5);
        fQAHistos->CreateTH2(Form("EPCentPos%s%s", fgkTriggerTypeNames[itype].Data(), patchtype.c_str()), Form("Center position of the %s %s trigger patches;#eta;#phi", patchtype.c_str(), fgkTriggerTypeNames[itype].Data()), 20, -0.8, 0.8, 700, 0., 7.);
        fQAHistos->CreateTH2(Form("PatchADCvsE%s%s", fgkTriggerTypeNames[itype].Data(), patchtype.c_str()), Form("Patch ADC value for trigger type %s %s;Trigger ADC;FEE patch energy (GeV)", patchtype.c_str(), fgkTriggerTypeNames[itype].Data()), 2000, 0., 2000, 200, 0., 200);
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
  AliEMCALTriggerPatchInfo  *trigger(nullptr), *triggerMainJet(nullptr),
                            *triggerMainGamma(nullptr), *triggerMainLevel0(nullptr),
                            *triggerMainJetSimple(nullptr), *triggerMainGammaSimple(nullptr),
                            *triggerMainJetRecalc(nullptr), *triggerMainGammaRecalc(nullptr);

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
  fPatchAmplitudes->Reset();
  fPatchADC->Reset();
  fPatchADCSimple->Reset();
  fLevel0TimeMap->Reset();

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
	    (*fPatchADC)(globCol,globRow) = adcAmp;

      // Handling for L0 triggers
      // For the ADC value we use fCaloTriggers->GetAmplitude()
      // In data, all patches which have 4 TRUs with proper level0 times are
      // valid trigger patches. Therefore we need to check all neighbors for
      // the level0 times, not only the bottom left. In order to obtain this
      // information, a lookup table with the L0 times for each TRU is created
      Float_t amplitude(0);
      fCaloTriggers->GetAmplitude(amplitude);
      if(amplitude < 0) amplitude = 0;
      (*fPatchAmplitudes)(globCol,globRow) = amplitude;
      Int_t nl0times(0);
      fCaloTriggers->GetNL0Times(nl0times);
      if(nl0times){
        TArrayI l0times(nl0times);
        fCaloTriggers->GetL0Times(l0times.GetArray());
        for(int itime = 0; itime < nl0times; itime++){
          if(l0times[itime] >7 && l0times[itime] < 10){
            (*fLevel0TimeMap)(globCol,globRow) = static_cast<Char_t>(l0times[itime]);
            break;
          }
        }
      }
    } // patches
  } // array not empty
  
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
    (*fPatchADCSimple)(globCol,globRow) += amp/EMCALTrigger::kEMCL1ADCtoGeV;
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

    // go throuth the trigger channels, real first, then offline
    Bool_t isOfflineSimple=0;
    while (NextTrigger(isOfflineSimple)) {
      if(isOfflineSimple){
        // process jet offline
        trigger = ProcessPatch(kTMEMCalJet, kTMOffline);
        // save main jet triggers in event
        if (trigger != 0) {
          // check if more energetic than others for main patch marking
          if (triggerMainJetSimple == 0 || (triggerMainJetSimple->GetPatchE() < trigger->GetPatchE()))
            triggerMainJetSimple = trigger;
        }

        // process jet recalc
        trigger = ProcessPatch(kTMEMCalJet, kTMRecalc);
        // save main jet triggers in event
        if (trigger != 0) {
          // check if more energetic than others for main patch marking
          if (triggerMainJetRecalc == 0 || (triggerMainJetRecalc->GetPatchE() < trigger->GetPatchE()))
            triggerMainJetRecalc = trigger;
        }

        // process gamma offline
        trigger = ProcessPatch(kTMEMCalGamma, kTMOffline);
        // save main gamma triggers in event
        if (trigger != 0) {
          // check if more energetic than others for main patch marking
          if (triggerMainGammaSimple == 0 || (triggerMainGammaSimple->GetPatchE() < trigger->GetPatchE()))
            triggerMainGammaSimple = trigger;
        }

        // process gamma recalc
        trigger = ProcessPatch(kTMEMCalGamma, kTMRecalc);
        // save main gamma triggers in event
        if (trigger != 0) {
          // check if more energetic than others for main patch marking
          if (triggerMainGammaRecalc == 0 || (triggerMainGammaRecalc->GetPatchE() < trigger->GetPatchE()))
            triggerMainGammaRecalc = trigger;
        }

      } else {
        // process jet
        trigger = ProcessPatch(kTMEMCalJet, kTMOnline);
        // save main jet triggers in event
        if (trigger != 0) {
          // check if more energetic than others for main patch marking
          if(triggerMainJet == 0 || (triggerMainJet->GetPatchE() < trigger->GetPatchE()))
            triggerMainJet = trigger;
        }

        trigger = ProcessPatch(kTMEMCalGamma, kTMOnline);
        // save main gamma triggers in event
        if (trigger != 0) {
          // check if more energetic than others for main patch marking
            if (triggerMainGamma == 0 || (triggerMainGamma->GetPatchE() < trigger->GetPatchE()))
              triggerMainGamma = trigger;
        }
        trigger = ProcessPatch(kTMEMCalLevel0, kTMOnline);
        // save main level0 trigger in the event
        if (trigger) {
          if (!triggerMainLevel0 || (triggerMainLevel0->GetPatchE() < trigger->GetPatchE()))
            triggerMainLevel0 = trigger;
        }
      }
    } // triggers
    
    // mark the most energetic patch as main
    // for real and also simple offline
    if (triggerMainJet != 0) {
      Int_t tBits = triggerMainJet->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << kMainTriggerBitNum );
      triggerMainJet->SetTriggerBits( tBits );
    }
    if (triggerMainJetSimple != 0) {
      Int_t tBits = triggerMainJetSimple->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << kMainTriggerBitNum );
      triggerMainJetSimple->SetTriggerBits(tBits);
    }
    if (triggerMainJetRecalc != 0) {
      Int_t tBits = triggerMainJetRecalc->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << kMainTriggerBitNum );
      triggerMainJetRecalc->SetTriggerBits(tBits);
    }
    if (triggerMainGamma != 0) {
      Int_t tBits = triggerMainGamma->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << kMainTriggerBitNum );
      triggerMainGamma->SetTriggerBits( tBits );
    }
    if (triggerMainGammaSimple != 0) {
      Int_t tBits = triggerMainGammaSimple->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << kMainTriggerBitNum );
      triggerMainGammaSimple->SetTriggerBits( tBits );
    }
    if (triggerMainGammaRecalc != 0) {
      Int_t tBits = triggerMainGammaRecalc->GetTriggerBits();
      // main trigger flag
      tBits = tBits | ( 1 << kMainTriggerBitNum );
      triggerMainGammaRecalc->SetTriggerBits( tBits );
    }
    if(triggerMainLevel0){
      Int_t tBits = triggerMainLevel0->GetTriggerBits();
      // main trigger flag
      tBits |= (1 << kMainTriggerBitNum);
      triggerMainLevel0->SetTriggerBits(tBits);
    }
  } // there are some triggers

  // Diagnostics
  int npatchOnline = 0;
  for(TIter patchIter = TIter(fCaloTriggersOut).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEMCALTriggerPatchInfo *mypatch = static_cast<AliEMCALTriggerPatchInfo *>(*patchIter);
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
AliEMCALTriggerPatchInfo* AliEmcalTriggerMaker::ProcessPatch(TriggerMakerTriggerType_t type, TriggerMakerPatchSource_t patchSource)
{
  Int_t tBits=-1;
  if (patchSource == kTMOnline)
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
      (type == kTMEMCalLevel0 && !(CheckForL0(*fCaloTriggers))))
    return 0;

  if((patchSource == kTMOffline && !IsOfflineSimple(tBits)) ||
     (patchSource == kTMRecalc && !IsRecalc(tBits)))
    return 0;

  // save primary vertex in vector
  TVector3 vertex;
  vertex.SetXYZ(fVertex[0], fVertex[1], fVertex[2]);

  // get position in global 2x2 tower coordinates
  // A0 left bottom (0,0)
  Int_t globCol=-1, globRow=-1;
  if (patchSource == kTMOnline)
    fCaloTriggers->GetPosition(globCol,globRow);
  else
    fSimpleOfflineTriggers->GetPosition(globCol, globRow);

  // Markus: For the moment reject jet patches with a row larger than 44 to overcome
  // an issue with patches containing inactive TRUs
  // (last 2 supermodules inactive but still included in the reconstruction)
  Int_t runno = InputEvent()->GetRunNumber();
  if(runno > 176000 && runno <= 197692){
    // Valid only for 2012 geometry
    if((type == kTMEMCalJet && IsEJE( tBits )) && (globRow > 44)) { // Hard coded number in order to be insensitive to changes in the geometry
      AliDebug(1, Form("Jet patch in inactive area: row[%d]", globRow));
      return NULL;
    }
  }

  if(fRejectOffAcceptancePatches){
    int patchsize = 2;
    const int kRowsPhi = fGeom->GetNTotalTRU() * 2;
    if(type == kTMEMCalJet) patchsize = 16;
    if((globCol + patchsize >= kColsEta) || (globCol + patchsize >= kRowsPhi)){
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
  Int_t colEdge1 = globCol, rowEdge1 = globRow, absIdEdge1 = absId, cellIdEdge1 = cellAbsId[0]; // Used in warning for invalid patch position
	
  // sum the available energy in the 32/32 window of cells
  // step over trigger channels and get all the corresponding cells
  // make CM
  Float_t amp = 0;
  Float_t cmiCol = 0;
  Float_t cmiRow = 0;
  Int_t adcAmp = 0;
  Double_t adcOfflineAmp = 0;
  int nfastor = (type == kTMEMCalJet) ? 16 : 2; // 32x32 cell window for L1 Jet trigger, 4x4 for L1 Gamma or L0 trigger
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
        try {
          adcAmp += static_cast<Int_t>((*fPatchAmplitudes)(globCol+i,globRow+j) * 4); // precision loss in case of global integer field
        } catch (const AliEMCALTriggerDataGrid<float>::OutOfBoundsException &e) {
          if(fDebugLevel){
            std::cerr << e.what() << std::endl;
          }
        }
      } else {
        try {
          adcAmp += (*fPatchADC)(globCol+i,globRow+j);
        } catch (AliEMCALTriggerDataGrid<int>::OutOfBoundsException &e){
          if(fDebugLevel){
            std::cerr << e.what() << std::endl;
          }
        }
      }

      try{
        adcOfflineAmp += (*fPatchADCSimple)(globCol+i,globRow+j);
      } catch (AliEMCALTriggerDataGrid<double>::OutOfBoundsException &e){
        if(fDebugLevel){
          std::cerr << e.what() << std::endl;
        }
      }
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
  Int_t colEdge2 = globCol+posOffset, rowEdge2 = globRow+posOffset, absIdEdge2 = absId, cellIdEdge2 = cellAbsId[3]; // Used in warning for invalid patch position
	
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
  TVector3 edge1tmp = edge1, edge2tmp = edge2; // Used in warning for invalid patch position
  centerGeo -= vertex;
  centerMass -= vertex;
  edge1 -= vertex;
  edge2 -= vertex;
  // Check for invalid patch positions
  if(!(edge1[0] || edge1[1] || edge1[2])){
    AliWarning(Form("Inconsistency in patch position for edge1: [%f|%f|%f]", edge1[0], edge1[1], edge1[2]));
    AliWarning("Original vectors:");
    AliWarning(Form("edge1: [%f|%f|%f]", edge1tmp[0], edge1tmp[1], edge1tmp[2]));
    AliWarning(Form("vertex: [%f|%f|%f]", vertex[0], vertex[1], vertex[2]));
    AliWarning(Form("Col: %d, Row: %d, FABSID: %d, Cell: %d", colEdge1, rowEdge1, absIdEdge1, cellIdEdge1));
    AliWarning(Form("Offline: %s", (patchSource == kTMOffline || patchSource == kTMRecalc) ? "yes" : "no"));
  }
  if(!(edge2[0] || edge2[1] || edge2[2])){
    AliWarning(Form("Inconsistency in patch position for edge2: [%f|%f|%f]", edge2[0], edge2[1], edge2[2]));
    AliWarning("Original vectors:");
    AliWarning(Form("edge2: [%f|%f|%f]", edge2tmp[0], edge2tmp[1], edge2tmp[2]));
    AliWarning(Form("vertex: [%f|%f|%f]", vertex[0], vertex[1], vertex[2]));
    AliWarning(Form("Col: %d, Row: %d, FABSID: %d, Cell: %d", colEdge2, rowEdge2, absIdEdge2, cellIdEdge2));
    AliWarning(Form("Offline: %s", (patchSource == kTMOffline || patchSource == kTMRecalc) ? "yes" : "no"));
  }

  Int_t isMC = MCEvent() ? 1 : 0;
  Int_t offSet = (1 - isMC) * fTriggerBitConfig->GetTriggerTypesEnd();
	
  // fix tbits .. remove the unwanted type triggers
  // for Jet and Gamma triggers we remove also the level 0 bit since it will be stored in the level 0 patch
  // for level 0 we remove all gamma and jet trigger bits
  if(patchSource == kTMOnline){
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
  }

  // Remove online bits from offline and recalc patches
  if(patchSource == kTMRecalc){
    tBits = tBits & kRecalcBitmask;
    // remove gamma bits from jet patches && vice versa
    if(type == kTMEMCalJet)
      tBits = tBits & (1 << (kRecalcOffset + fTriggerBitConfig->GetJetLowBit()) | 1 << (kRecalcOffset + fTriggerBitConfig->GetJetHighBit()));
    else
      tBits = tBits & (1 << (kRecalcOffset + fTriggerBitConfig->GetGammaLowBit()) | 1 << (kRecalcOffset + fTriggerBitConfig->GetGammaHighBit()));
  } else if(patchSource == kTMOffline){
    tBits = tBits & kOfflineBitmask;
    // remove gamma bits from jet patches && vice versa
    if(type == kTMEMCalJet)
      tBits = tBits & (1 << (kOfflineOffset + fTriggerBitConfig->GetJetLowBit()) | 1 << (kOfflineOffset + fTriggerBitConfig->GetJetHighBit()));
    else
      tBits = tBits & (1 << (kOfflineOffset + fTriggerBitConfig->GetGammaLowBit()) | 1 << (kOfflineOffset + fTriggerBitConfig->GetGammaHighBit()));
  }

  // save the trigger object
  AliEMCALTriggerPatchInfo *trigger =
    new ((*fCaloTriggersOut)[fITrigger]) AliEMCALTriggerPatchInfo();
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
  trigger->SetCol0(globCol);
  trigger->SetRowStart(globRow);
  trigger->SetEdgeCell(globCol*2, globRow*2); // from triggers to cells
  //if(isOfflineSimple)trigger->SetOfflineSimple();
  if(fDoQA){
    TString patchtype;
    switch(patchSource){
    case kTMOnline: patchtype =  "Online"; break;
    case kTMOffline: patchtype = "Offline"; break;
    case kTMRecalc: patchtype = "Recalc"; break;
    };
    fQAHistos->FillTH2(Form("RCPos%s%s", fgkTriggerTypeNames[type].Data(), patchtype.Data()), globCol, globRow);
    fQAHistos->FillTH2(Form("EPCentPos%s%s", fgkTriggerTypeNames[type].Data(), patchtype.Data()), centerGeo.Eta(), centerGeo.Phi());
    fQAHistos->FillTH2(Form("PatchADCvsE%s%s", fgkTriggerTypeNames[type].Data(), patchtype.Data()), (patchSource == kTMOffline) ? adcOfflineAmp : adcAmp, trigger->GetPatchE());
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

  // @TODO: inefficient implementation: Array needs to be
  // copied for every new patch
  TArrayI tBitsArray(4), rowArray(4), colArray(4);
  
  // First entries are for flagging the main patches
  // Logics:
  // 0. Recalc jet
  // 1. Offline jet
  // 2. Recalc gamma
  // 3. Offline gamma
  // Afterwards come all other entries. Attention:
  // In order to prevent multiple counting the the patch information
  // is written to the trigger stream only from index 4 on

  // Max Recalc jet bit
  tBitsArray[0] = 0;
  colArray[0] = -1;
  rowArray[0] = -1;

  // Max offline jet patch
  tBitsArray[1] = 0;
  colArray[1] = -1;
  rowArray[1] = -1;

  // Max recalc gamma patch
  tBitsArray[2] = 0;
  colArray[2] = -1;
  rowArray[2] = -1;

  // Max offline gamma patch
  tBitsArray[3] = 0;
  colArray[3] = -1;
  rowArray[3] = -1;

  Double_t maxPatchADCoffline = -1;
  Int_t maxPatchADC = -1;
  // run the trigger algo, stepping by 8 towers (= 4 trigger channels)
  Int_t maxCol = 48;
  Int_t maxRow = fGeom->GetNTotalTRU()*2;
  // Markus:
  // temp fix for the number of TRUs in the 2011 PbPb data to 30
  // @TODO: Fix in the geometry in the OCDB
  int runnumber = InputEvent()->GetRunNumber();
  if(runnumber > 139517 && runnumber <= 170593) maxRow = 60;
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
          tSumOffline += (*fPatchADCSimple)(i+k,j+l);
          tSum += static_cast<ULong64_t>((*fPatchADC)(i+k,j+l));
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
      // first we check the offline energy compared to thresholds and set the offline bits accordingly
      // second we check the recalc energy compared to thresholds and set the recalc bits accordingly
      if (tSumOffline > fCaloTriggerSetupOut->GetThresholdJetLowSimple()){
        // Set the trigger bit for jet low - it is needed by the ProcessPatch function
        // in order to handle jet patches - for offline and recalc patches these bits will
        // be removed later
        tBits = tBits | ( 1 << ( bitOffSet + fTriggerBitConfig->GetJetLowBit() ));
        // Add offline bit - it will be handled by the ProcessPatch function
        tBits = tBits | ( 1 << (kOfflineOffset + fTriggerBitConfig->GetJetLowBit()) );
      }
      if (tSumOffline > fCaloTriggerSetupOut->GetThresholdJetHighSimple()){
        // Set the trigger bit for jet high - it is needed by the ProcessPatch function
        // in order to handle jet patches - for offline and recalc patches these bits will
        // be removed later
        tBits = tBits | ( 1 << ( bitOffSet + fTriggerBitConfig->GetJetHighBit() ));
        // Add offline bit - it will be handled by the ProcessPatch function
        tBits = tBits | ( 1 << (kOfflineOffset + fTriggerBitConfig->GetJetHighBit()) );
      }
      if(tSum > fCaloTriggerSetupOut->GetThresholdJetLowSimple()){
        // Set the trigger bit for jet low - it is needed by the ProcessPatch function
        // in order to handle jet patches - for offline and recalc patches these bits will
        // be removed later
        tBits = tBits | ( 1 << ( bitOffSet + fTriggerBitConfig->GetJetLowBit() ));
        // Add recalc bit - it will be handled by the ProcessPatch function
        tBits = tBits | ( 1 << (kRecalcOffset + fTriggerBitConfig->GetJetLowBit()) );
      }
      if (tSum > fCaloTriggerSetupOut->GetThresholdJetHighSimple()){
        // Set the trigger bit for jet high - it is needed by the ProcessPatch function
        // in order to handle jet patches - for offline and recalc patches these bits will
        // be removed later
        tBits = tBits | ( 1 << ( bitOffSet + fTriggerBitConfig->GetJetHighBit() ));
        // Add recalc bit - it will be handled by the ProcessPatch function
        tBits = tBits | ( 1 << (kRecalcOffset + fTriggerBitConfig->GetJetHighBit()) );
      }
      
      // add trigger values
      if (tBits != 0) {
        tBitsArray.Set( tBitsArray.GetSize() + 1 );
        colArray.Set( colArray.GetSize() + 1 );
        rowArray.Set( rowArray.GetSize() + 1 );
        tBitsArray[tBitsArray.GetSize()-1] = tBits;
        colArray[colArray.GetSize()-1] = i;
        rowArray[rowArray.GetSize()-1] = j;
      }
    }
  } // trigger algo
  
  // Set trigger bits for the maximum patch
  if (maxPatchADC > fCaloTriggerSetupOut->GetThresholdJetLowSimple()){
    // Set the trigger bit for jet low - it is needed by the ProcessPatch function
    // in order to handle jet patches - for offline and recalc patches these bits will
    // be removed later
    tBitsArray[0] = tBitsArray[0] | ( 1 << ( bitOffSet + fTriggerBitConfig->GetJetLowBit() ));
    // Add recalc bit - it will be handled by the ProcessPatch function
    tBitsArray[0] = tBitsArray[0] | ( 1 << (kRecalcOffset + fTriggerBitConfig->GetJetLowBit()) );
  }
  if (maxPatchADC > fCaloTriggerSetupOut->GetThresholdJetHighSimple()){
    // Set the trigger bit for jet high - it is needed by the ProcessPatch function
    // in order to handle jet patches - for offline and recalc patches these bits will
    // be removed later
    tBitsArray[0] = tBitsArray[0] | ( 1 << ( bitOffSet + fTriggerBitConfig->GetJetHighBit() ));
    // Add recalc bit - it will be handled by the ProcessPatch function
    tBitsArray[0] = tBitsArray[0] | ( 1 << (kRecalcOffset + fTriggerBitConfig->GetJetHighBit()) );
  }
  if(maxPatchADCoffline > fCaloTriggerSetupOut->GetThresholdJetLowSimple()){
    // Set the trigger bit for jet low - it is needed by the ProcessPatch function
    // in order to handle jet patches - for offline and recalc patches these bits will
    // be removed later
    tBitsArray[1] = tBitsArray[1] | ( 1 << ( bitOffSet + fTriggerBitConfig->GetJetLowBit() ));
    // Add offline bit - it will be handled by the ProcessPatch function
    tBitsArray[1] = tBitsArray[1] | ( 1 << (kOfflineOffset + fTriggerBitConfig->GetJetLowBit()) );
  }
  if (maxPatchADCoffline > fCaloTriggerSetupOut->GetThresholdJetHighSimple()){
    // Set the trigger bit for jet high - it is needed by the ProcessPatch function
    // in order to handle jet patches - for offline and recalc patches these bits will
    // be removed later
    tBitsArray[1] = tBitsArray[1] | ( 1 << ( bitOffSet + fTriggerBitConfig->GetJetHighBit() ));
    // Add offline bit - it will be handled by the ProcessPatch function
    tBitsArray[1] = tBitsArray[1] | ( 1 << (kOfflineOffset + fTriggerBitConfig->GetJetHighBit()) );
  }



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
          tSumOffline += (*fPatchADCSimple)(i+k,j+l);
          tSum += static_cast<ULong64_t>((*fPatchADC)(i+k,j+l));
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
      if (tSumOffline > fCaloTriggerSetupOut->GetThresholdGammaLowSimple()){
        // Set the trigger bit for gamma low - it is needed by the ProcessPatch function
        // in order to handle gamma patches - for offline and recalc patches these bits will
        // be removed later
        tBits = tBits | ( 1 << ( bitOffSet + fTriggerBitConfig->GetGammaLowBit() ));
        // Add offline bit - it will be handled by the ProcessPatch function
        tBits = tBits | ( 1 << (kOfflineOffset + fTriggerBitConfig->GetGammaLowBit()) );
      }
      if (tSumOffline > fCaloTriggerSetupOut->GetThresholdGammaHighSimple()){
        // Set the trigger bit for gamma high - it is needed by the ProcessPatch function
        // in order to handle gamma patches - for offline and recalc patches these bits will
        // be removed later
        tBits = tBits | ( 1 << ( bitOffSet + fTriggerBitConfig->GetGammaHighBit() ));
        // Add offline bit - it will be handled by the ProcessPatch function
        tBits = tBits | ( 1 << (kOfflineOffset + fTriggerBitConfig->GetGammaHighBit()) );
      }
      if (tSum > fCaloTriggerSetupOut->GetThresholdGammaLowSimple()){
        // Set the trigger bit for gamma low - it is needed by the ProcessPatch function
        // in order to handle gamma patches - for offline and recalc patches these bits will
        // be removed later
        tBits = tBits | ( 1 << ( bitOffSet + fTriggerBitConfig->GetGammaLowBit() ));
        // Add recalc bit - it will be handled by the ProcessPatch function
        tBits = tBits | ( 1 << (kOfflineOffset + fTriggerBitConfig->GetGammaLowBit()) );
      }
      if (tSum > fCaloTriggerSetupOut->GetThresholdGammaHighSimple()){
        // Set the trigger bit for gamma high - it is needed by the ProcessPatch function
        // in order to handle gamma patches - for offline and recalc patches these bits will
        // be removed later
        tBits = tBits | ( 1 << ( bitOffSet + fTriggerBitConfig->GetGammaHighBit() ));
        // Add recalc bit - it will be handled by the ProcessPatch function
        tBits = tBits | ( 1 << (kOfflineOffset + fTriggerBitConfig->GetGammaHighBit()) );
      }
      
      // add trigger values
      if (tBits != 0) {
        tBitsArray.Set( tBitsArray.GetSize() + 1 );
        colArray.Set( colArray.GetSize() + 1 );
        rowArray.Set( rowArray.GetSize() + 1 );
        tBitsArray[tBitsArray.GetSize()-1] = tBits;
        colArray[colArray.GetSize()-1] = i;
        rowArray[rowArray.GetSize()-1] = j;
      }
    }
  } // trigger algo
  
  // Set bits for the maximum patch
  if (maxPatchADC > fCaloTriggerSetupOut->GetThresholdGammaLowSimple()){
    // Set the trigger bit for gamma low - it is needed by the ProcessPatch function
    // in order to handle gamma patches - for offline and recalc patches these bits will
    // be removed later
    tBitsArray[2] = tBitsArray[2] | ( 1 << ( bitOffSet + fTriggerBitConfig->GetGammaLowBit() ));
    // Add recalc bit - it will be handled by the ProcessPatch function
    tBitsArray[2] = tBitsArray[2] | ( 1 << (kRecalcOffset + fTriggerBitConfig->GetGammaLowBit()) );
  }
  if (maxPatchADC > fCaloTriggerSetupOut->GetThresholdGammaHighSimple()){
    // Set the trigger bit for gamma high - it is needed by the ProcessPatch function
    // in order to handle gamma patches - for offline and recalc patches these bits will
    // be removed later
    tBitsArray[2] = tBitsArray[2] | ( 1 << ( bitOffSet + fTriggerBitConfig->GetGammaHighBit() ));
    // Add recalc bit - it will be handled by the ProcessPatch function
    tBitsArray[2] = tBitsArray[2] | ( 1 << (kRecalcOffset + fTriggerBitConfig->GetGammaHighBit()) );
  }
  if (maxPatchADCoffline > fCaloTriggerSetupOut->GetThresholdGammaLowSimple()){
    // Set the trigger bit for gamma low - it is needed by the ProcessPatch function
    // in order to handle gamma patches - for offline and recalc patches these bits will
    // be removed later
    tBitsArray[3] = tBitsArray[3] | ( 1 << ( bitOffSet + fTriggerBitConfig->GetGammaLowBit() ));
    // Add offline bit - it will be handled by the ProcessPatch function
    tBitsArray[3] = tBitsArray[3] | ( 1 << (kOfflineOffset + fTriggerBitConfig->GetGammaLowBit()) );
  }
  if (maxPatchADCoffline > fCaloTriggerSetupOut->GetThresholdGammaHighSimple()){
    // Set the trigger bit for gamma high - it is needed by the ProcessPatch function
    // in order to handle gamma patches - for offline and recalc patches these bits will
    // be removed later
    tBitsArray[3] = tBitsArray[3] | ( 1 << ( bitOffSet + fTriggerBitConfig->GetGammaHighBit() ));
    // Add offline bit - it will be handled by the ProcessPatch function
    tBitsArray[3] = tBitsArray[3] | ( 1 << (kOfflineOffset + fTriggerBitConfig->GetGammaHighBit()) );
  }

  // save in object
  fSimpleOfflineTriggers->DeAllocate();
  if(tBitsArray.GetSize() - 4 > 0){
    fSimpleOfflineTriggers->Allocate(tBitsArray.GetSize() - 4);
    for (Int_t i = 4; i < tBitsArray.GetSize(); ++i){
      fSimpleOfflineTriggers->Add(colArray[i],rowArray[i], 0, 0, 0, 0, 0, tBitsArray[i]);
    }
  }

  // @TODO: Implement QA of the main patches
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
  Int_t row(-1), col(-1); trg.GetPosition(col, row);
  if(col < 0 || row < 0){
    AliError(Form("Patch outside range [col %d, row %d]", col, row));
    return kFALSE;
  }
  Int_t truref(-1), trumod(-1), absFastor(-1), adc(-1);
  fGeom->GetAbsFastORIndexFromPositionInEMCAL(col, row, absFastor);
  fGeom->GetTRUFromAbsFastORIndex(absFastor, truref, adc);
  int nvalid(0);
  const int kNRowsPhi = fGeom->GetNTotalTRU() * 2;
  for(int ipos = 0; ipos < 2; ipos++){
    if(row + ipos >= kNRowsPhi) continue;    // boundary check
    for(int jpos = 0; jpos < 2; jpos++){
      if(col + jpos >= kColsEta) continue;  // boundary check
      // Check whether we are in the same TRU
      trumod = -1;
      fGeom->GetAbsFastORIndexFromPositionInEMCAL(col+jpos, row+ipos, absFastor);
      fGeom->GetTRUFromAbsFastORIndex(absFastor, trumod, adc);
      if(trumod != truref) continue;
      if(col + jpos >= kColsEta) AliError(Form("Boundary error in col [%d, %d + %d]", col + jpos, col, jpos));
      if(row + ipos >= kNRowsPhi) AliError(Form("Boundary error in row [%d, %d + %d]", row + ipos, row, ipos));
      Char_t l0times = (*fLevel0TimeMap)(col + jpos,row + ipos);
      if(l0times > 7 && l0times < 10) nvalid++;
    }
  }
  if (nvalid != 4) return false;
  return true;
}
