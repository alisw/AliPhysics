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
#include <TGrid.h>
#include <THashList.h>
#include <THistManager.h>
#include <TObjArray.h>
#include <TParameter.h>

#include "AliEMCALTriggerBitConfig.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerMakerKernel.h"
#include "AliEmcalTriggerMakerTask.h"
#include "AliLog.h"
#include "AliOADBContainer.h"

#include <bitset>
#include <sstream>
#include <string>

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerMakerTask)
/// \endcond

AliEmcalTriggerMakerTask::AliEmcalTriggerMakerTask():
  AliAnalysisTaskEmcal(),
  fTriggerMaker(NULL),
  fV0(NULL),
  fCaloTriggersOutName("EmcalTriggers"),
  fV0InName("AliAODVZERO"),
  fBadFEEChannelOADB(""),
  fUseL0Amplitudes(kFALSE),
  fCaloTriggersOut(0),
  fDoQA(kFALSE),
  fQAHistos(NULL)
{

}

AliEmcalTriggerMakerTask::AliEmcalTriggerMakerTask(const char *name, Bool_t doQA):
  AliAnalysisTaskEmcal("AliEmcalTriggerMakerTask", doQA),
  fTriggerMaker(NULL),
  fV0(NULL),
  fCaloTriggersOutName("EmcalTriggers"),
  fV0InName("AliAODVZERO"),
  fBadFEEChannelOADB(""),
  fUseL0Amplitudes(kFALSE),
  fCaloTriggersOut(NULL),
  fDoQA(doQA),
  fQAHistos(NULL)
{
  fTriggerMaker = new AliEmcalTriggerMakerKernel;
}

AliEmcalTriggerMakerTask::~AliEmcalTriggerMakerTask() {
  if(fTriggerMaker) delete fTriggerMaker;
}

/**
 * Initialize output objets
 */
void AliEmcalTriggerMakerTask::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  const TString kTriggerTypeNames[3] = {"EJE", "EGA", "EL0"},
                kPatchTypes[3] = {"Online", "Offline", "Recalc"};

  if(fDoQA && fOutput){
    fQAHistos = new THistManager("TriggerQA");

    for(const TString *triggertype = kTriggerTypeNames; triggertype < kTriggerTypeNames + sizeof(kTriggerTypeNames)/sizeof(TString); triggertype++){
      for(const TString *patchtype = kPatchTypes; patchtype < kPatchTypes + sizeof(kPatchTypes)/sizeof(TString); ++patchtype){
        fQAHistos->CreateTH2(
            Form("RCPos%s%s", triggertype->Data(), patchtype->Data()),
            Form("Lower edge position of %s %s patches (col-row);iEta;iPhi", patchtype->Data(), triggertype->Data()),
            48, -0.5, 47.5, 104, -0.5, 103.5
            );
        fQAHistos->CreateTH2(
            Form("EPCentPos%s%s", triggertype->Data(), patchtype->Data()),
            Form("Center position of the %s %s trigger patches;#eta;#phi", patchtype->Data(), triggertype->Data()),
            20, -0.8, 0.8, 700, 0., 7.
            );
        fQAHistos->CreateTH2(
            Form("PatchADCvsE%s%s", triggertype->Data(), patchtype->Data()),
            Form("Patch ADC value for trigger type %s %s;Trigger ADC;FEE patch energy (GeV)", patchtype->Data(), triggertype->Data()),
            2000, 0., 2000, 200, 0., 200
            );
        fQAHistos->CreateTH2(
            Form("PatchADCOffvsE%s%s", triggertype->Data(), patchtype->Data()),
            Form("Patch offline ADC value for trigger type %s %s;Trigger ADC;FEE patch energy (GeV)", patchtype->Data(), triggertype->Data()),
            2000, 0., 2000, 200, 0., 200
            );
      }
    }
    fQAHistos->CreateTH1("triggerBitsAll", "Trigger bits for all incoming patches;bit nr", 64, -0.5, 63.5);
    fQAHistos->CreateTH1("triggerBitsSel", "Trigger bits for reconstructed patches;bit nr", 64, -0.5, 63.5);
    fOutput->Add(fQAHistos->GetListOfHistograms());
    PostData(1, fOutput);
  }
}

void AliEmcalTriggerMakerTask::SetUseTriggerBitConfig(TriggerMakerBitConfig_t bitConfig)
{
  AliEMCALTriggerBitConfig *triggerBitConfig(NULL);
  switch(bitConfig){
  case kNewConfig:
    triggerBitConfig = new AliEMCALTriggerBitConfigNew();
    break;
  case kOldConfig:
    triggerBitConfig = new AliEMCALTriggerBitConfigOld();
    break;
  }
  fTriggerMaker->SetTriggerBitConfig(triggerBitConfig);
}

/**
 * Initializes the trigger maker kernel
 */
void AliEmcalTriggerMakerTask::ExecOnce(){
  AliAnalysisTaskEmcal::ExecOnce();

  if (!fInitialized) return;

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

  if (!fV0InName.IsNull()) {
    fV0 = (AliVVZERO*)InputEvent()->FindListObject(fV0InName);
  }

  fTriggerMaker->SetGeometry(fGeom);
  fTriggerMaker->Init();
}

/**
 * Run the trigger maker
 * Move patches found by the trigger maker to the output clones array
 * Fill QA histograms if requested
 * @return True
 */
Bool_t AliEmcalTriggerMakerTask::Run(){
  fCaloTriggersOut->Clear();
  // prepare trigger maker
  fTriggerMaker->Reset();
  fTriggerMaker->ReadCellData(fCaloCells);
  fTriggerMaker->ReadTriggerData(fCaloTriggers);
  fTriggerMaker->BuildL1ThresholdsOffline(fV0);
  fTriggerMaker->SetIsMC(MCEvent());
  TObjArray *patches = fTriggerMaker->CreateTriggerPatches(InputEvent(), fUseL0Amplitudes);
  AliEMCALTriggerPatchInfo *recpatch = NULL;
  Int_t patchcounter = 0;
  TString triggerstring;
  AliDebug(2,Form("Trigger maker - Found %d patches\n", patches->GetEntries()));
  for(TIter patchIter = TIter(patches).Begin(); patchIter != TIter::End(); ++patchIter){
    recpatch = dynamic_cast<AliEMCALTriggerPatchInfo *>(*patchIter);
    if(fDoQA){
      //std::bitset<32> triggerbits = recpatch->GetTriggerBits();
      std::stringstream triggerbitstring;
      AliDebug(1, Form("Trigger maker - next patch: size %d, trigger bits %s", recpatch->GetPatchSize(), triggerbitstring.str().c_str()));
      // Handle types different - online - offline - re
      if(recpatch->IsJetHigh() || recpatch->IsJetLow())                   FillQAHistos("EJEOnline", *recpatch);
      if(recpatch->IsGammaHigh() || recpatch->IsGammaLow())               FillQAHistos("EGAOnline", *recpatch);
      if(recpatch->IsJetHighSimple() || recpatch->IsJetLowSimple())       FillQAHistos("EJEOffline", *recpatch);
      if(recpatch->IsGammaHighSimple() || recpatch->IsGammaLowSimple())   FillQAHistos("EGAOffline", *recpatch);
      if(recpatch->IsLevel0())                                            FillQAHistos("EL0Online", *recpatch);
      if(recpatch->IsRecalcJet())                                         FillQAHistos("EJERecalc", *recpatch);
      if(recpatch->IsRecalcGamma())                                       FillQAHistos("EGARecalc", *recpatch);
      // Redo checking of found trigger bits after masking of unwanted triggers
      int tBits = recpatch->GetTriggerBits();
      for(unsigned int ibit = 0; ibit < sizeof(tBits)*8; ibit++) {
        if(tBits & (1 << ibit)){
          fQAHistos->FillTH1("triggerBitsSel", ibit);
        }
      }
    }
    new((*fCaloTriggersOut)[patchcounter++]) AliEMCALTriggerPatchInfo(*recpatch);
  }
  if(patches) delete patches;
  return true;
}

void AliEmcalTriggerMakerTask::RunChanged(){
  if(fBadFEEChannelOADB.Length()) InitializeBadFEEChannels();
}

void AliEmcalTriggerMakerTask::InitializeBadFEEChannels(){
  fTriggerMaker->ClearOfflineBadChannels();
  if(fBadFEEChannelOADB.Contains("alien://") && !gGrid) TGrid::Connect("alien://");
  AliOADBContainer badchannelDB("EmcalBadChannelsAdditional");
  TObjArray *badchannelmap = static_cast<TObjArray *>(badchannelDB.GetObject(InputEvent()->GetRunNumber()));
  if(!badchannelmap || !badchannelmap->GetEntries()) return;
  for(TIter citer = TIter(badchannelmap).Begin(); citer != TIter::End(); ++citer){
    TParameter<int> *channelID = static_cast<TParameter<int> *>(*citer);
    fTriggerMaker->AddOfflineBadChannel(channelID->GetVal());
  }
}

void AliEmcalTriggerMakerTask::FillQAHistos(const TString &patchtype, const AliEMCALTriggerPatchInfo &recpatch){
  fQAHistos->FillTH2(Form("RCPos%s", patchtype.Data()), recpatch.GetColStart(), recpatch.GetRowStart());
  fQAHistos->FillTH2(Form("EPCentPos%s", patchtype.Data()), recpatch.GetEtaGeo(), recpatch.GetPhiGeo());
  fQAHistos->FillTH2(Form("PatchADCvsE%s", patchtype.Data()), recpatch.GetADCAmp(), recpatch.GetPatchE());
  fQAHistos->FillTH2(Form("PatchADCOffvsE%s", patchtype.Data()), recpatch.GetADCOfflineAmp(), recpatch.GetPatchE());
}
