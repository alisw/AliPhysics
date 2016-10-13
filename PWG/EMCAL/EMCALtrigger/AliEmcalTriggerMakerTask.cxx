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

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerBitConfig.h"
#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerTRUDCSConfig.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerMakerKernel.h"
#include "AliEmcalTriggerMakerTask.h"
#include "AliEMCALTriggerMapping.h"
#include "AliLog.h"
#include "AliOADBContainer.h"

#include <bitset>
#include <iostream>
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
  fMaskedFastorOADB(""),
  fUseL0Amplitudes(kFALSE),
  fLoadFastORMaskingFromOCDB(kFALSE),
  fCaloTriggersOut(0),
  fDoQA(kFALSE),
  fQAHistos(NULL)
{

}

AliEmcalTriggerMakerTask::AliEmcalTriggerMakerTask(const char *name, Bool_t doQA):
  AliAnalysisTaskEmcal("AliEmcalTriggerMakerTask", kTRUE),
  fTriggerMaker(NULL),
  fV0(NULL),
  fCaloTriggersOutName("EmcalTriggers"),
  fV0InName("AliAODVZERO"),
  fBadFEEChannelOADB(""),
  fMaskedFastorOADB(""),
  fUseL0Amplitudes(kFALSE),
  fLoadFastORMaskingFromOCDB(kFALSE),
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
 * Initialize output objects
 */
void AliEmcalTriggerMakerTask::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  const TString kTriggerTypeNames[3] = {"EJE", "EGA", "EL0"},
                kPatchTypes[3] = {"Online", "Offline", "Recalc"};

  if(fDoQA) std::cout << "Trigger maker - QA requested" << std::endl;
  else std::cout << "Trigger maker - no QA requested" << std::endl;
  if(!fOutput) std::cout << "No output container initialized" << std::endl;

  if(fDoQA){
    if(fOutput){
      AliInfoStream() << "Enabling QA for trigger maker" << std::endl;
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
      fQAHistos->CreateTH2("FastORMaskOnline", "Masked FastORs at online level; col; row", 48, -0.5, 47.5, 104, -0.5, 103.5);
      fOutput->Add(fQAHistos->GetListOfHistograms());
      PostData(1, fOutput);
    } else {
      AliWarningStream() << "QA requested but no output container initialized - QA needs to be disabled" << std::endl;
      fDoQA = kFALSE;
    }
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

  if (!fLocalInitialized) return;

  if (!fCaloTriggersOutName.IsNull()) {
    fCaloTriggersOut = new TClonesArray("AliEMCALTriggerPatchInfo");
    fCaloTriggersOut->SetName(fCaloTriggersOutName);

    if (!(InputEvent()->FindListObject(fCaloTriggersOutName))) {
      InputEvent()->AddObject(fCaloTriggersOut);
    }
    else {
      fLocalInitialized = kFALSE;
      AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fCaloTriggersOutName.Data()));
      return;
    }
  }

  if (!fV0InName.IsNull()) {
    fV0 = (AliVVZERO*)InputEvent()->FindListObject(fV0InName);
  }

  // Configure trigger maker
  if(!fTriggerMaker->IsConfigured()){
    AliInfoStream() << "Trigger maker not yet configure - automatically configuring ..." << std::endl;
    int runnumber = InputEvent()->GetRunNumber();
    std::string dataset = "";
    if(runnumber >= 145144 && runnumber <= 165746){
      fTriggerMaker->ConfigureForPP2011();
      dataset = "pp 2011";
    } else if(runnumber >= 167806 && runnumber <= 170593){
      fTriggerMaker->ConfigureForPbPb2011();
      dataset = "Pb-Pb 2011";
    } else if(runnumber >= 176326 && runnumber <= 193766){
      fTriggerMaker->ConfigureForPP2012();
      dataset = "pp 2012";
    } else if(runnumber >= 195344 && runnumber <= 197692){
      // dataset contains also the setup for pp 2.76 TeV
      fTriggerMaker->ConfigureForPPb2013();
      dataset = "p-Pb 2013";
    } else if((runnumber >= 224891 && runnumber <= 244628) || runnumber >= 253434){
      // Configuration starting with LHC15f
      fTriggerMaker->ConfigureForPbPb2015();
      dataset = "pp 2015-2016";
    } else if(runnumber >= 244824 && runnumber <= 246994){
      fTriggerMaker->ConfigureForPbPb2015();
      dataset = "Pb-Pb 2015";
    }

    if(fTriggerMaker->IsConfigured()){
      AliInfoStream() << "Applying configuration for " << dataset << std::endl;
    } else {
      AliErrorStream() << "No valid configuration found for the given dataset - trigger maker disabled" << std::endl;
      fLocalInitialized = false;
    }
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
  AliDebugStream(2) << GetName() << ": Found " << patches->GetEntries() << " patches" << std::endl;
  for(TIter patchIter = TIter(patches).Begin(); patchIter != TIter::End(); ++patchIter){
    recpatch = dynamic_cast<AliEMCALTriggerPatchInfo *>(*patchIter);
    if(fDoQA){
      AliDebugStream(3) <<  GetName() << ": Next patch: size " << recpatch->GetPatchSize() << " , trigger bits " << std::bitset<sizeof(Int_t)*8>(recpatch->GetTriggerBits()) << std::endl;
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

void AliEmcalTriggerMakerTask::RunChanged(Int_t newrun){
  fTriggerMaker->ClearOfflineBadChannels();
  if(fBadFEEChannelOADB.Length()) InitializeBadFEEChannels();
  fTriggerMaker->ClearFastORBadChannels();
  if(fLoadFastORMaskingFromOCDB) InitializeFastORMaskingFromOCDB();
  if(fMaskedFastorOADB.Length()) InitializeFastORMaskingFromOADB();
  // QA: Monitor all channels which are masked in the current run
  if(fDoQA && fQAHistos){
    Int_t globCol(-1), globRow(-1) ;
    for(const auto &ifastOrID : fTriggerMaker->GetListOfBadFastORAbsIDs()){
      fGeom->GetTriggerMapping()->GetPositionInEMCALFromAbsFastORIndex(ifastOrID, globCol, globRow);
      fQAHistos->FillTH2("FastORMaskOnline", globCol, globRow);
    }
  }
}

void AliEmcalTriggerMakerTask::InitializeBadFEEChannels(){
  AliInfoStream() << "Loading additional bad FEE channels from OADB container " << fBadFEEChannelOADB << std::endl;
  fTriggerMaker->ClearOfflineBadChannels();
  if(fBadFEEChannelOADB.Contains("alien://") && !gGrid) TGrid::Connect("alien://");
  AliOADBContainer badchannelDB("EmcalBadChannelsAdditional");
  badchannelDB.InitFromFile(fBadFEEChannelOADB.Data(), "EmcalBadChannelsAdditional");
  TObjArray *badchannelmap = static_cast<TObjArray *>(badchannelDB.GetObject(InputEvent()->GetRunNumber()));
  if(!badchannelmap || !badchannelmap->GetEntries()) return;
  for(TIter citer = TIter(badchannelmap).Begin(); citer != TIter::End(); ++citer){
    TParameter<int> *channelID = static_cast<TParameter<int> *>(*citer);
    fTriggerMaker->AddOfflineBadChannel(channelID->GetVal());
  }
}

void AliEmcalTriggerMakerTask::InitializeFastORMaskingFromOCDB(){
  AliInfoStream() << "Loading masked fastors from OCDB" << std::endl;
  AliCDBManager *cdb = AliCDBManager::Instance();

  AliCDBEntry *en = cdb->Get("EMCAL/Calib/Trigger");
  if(!en){
    AliErrorStream() << GetName() << ": FastOR masking from CDB required, but OCDB entry is not available. No masking will be applied." << std::endl;
    return;
  }

  AliEMCALTriggerDCSConfig *trgconf = dynamic_cast<AliEMCALTriggerDCSConfig *>(en->GetObject());
  if(!trgconf){
    AliErrorStream() << GetName() << ": Failed decoding OCDB entry: Object is not of type AliEMCALTriggerDCSConfig." << std::endl;
    return;
  }

  // In run 1 the small supermodules were not contributing to triggers.
  // Still the TRUs are counted. As access to the TRU config is not properly
  // protected the loop over NTRU from the geometry will produce a segfault.
  // As temporary workaround the loop limits are obtained from the DCS data itself.
  Int_t fastOrAbsID(-1), ic(-1);
  for(int itru = 0; itru < trgconf->GetTRUArr()->GetEntries(); itru++){
    AliEMCALTriggerTRUDCSConfig *truconf = trgconf->GetTRUDCSConfig(itru);
    // Test for each channel whether it is masked. The calculation is
    // done reversely as the channel mapping is different between run1
    // and run2: The loop is done over all masks and all bits inside the
    // mask, and a handler matching to the correct mapping converts them
    // into the channel ID. In case a masked channel is found, the absolute
    // ID is calculated. For this the function GetAbsFastORIndexFromTRU
    // is used - it is assumed that parameter 1 (iADC) corresponds to the
    // channel ID.
    for(unsigned int ifield = 0; ifield < 6; ifield++){
      for(unsigned int ibit = 0; ibit < 16; ibit ++){
        if((truconf->GetMaskReg(ifield) >> ibit) & 0x1){
          try{
            fGeom->GetTriggerMapping()->GetAbsFastORIndexFromTRU(itru, (ic =  GetMaskHandler()(ifield, ibit)), fastOrAbsID);
            AliDebugStream(1) << GetName() << "Channel " << ic  << " in TRU " << itru << " ( abs fastor " << fastOrAbsID << ") masked." << std::endl;
            fTriggerMaker->AddFastORBadChannel(fastOrAbsID);
          } catch (int exept){
            AliErrorStream() << GetName() << "Invalid mask: (" << ifield << "|" << ibit << "), exception " << exept << " thrown. Mask will not be recognized" << std::endl;
          }
        }
      }
    }
  }
}

void AliEmcalTriggerMakerTask::InitializeFastORMaskingFromOADB(){
  AliInfoStream() << "Initializing masked fastors from OADB container " << fMaskedFastorOADB.Data() << std::endl;
  if(fMaskedFastorOADB.Contains("alien://") && !gGrid) TGrid::Connect("alien://");
  AliOADBContainer badchannelDB("AliEmcalMaskedFastors");
  badchannelDB.InitFromFile(fMaskedFastorOADB.Data(), "AliEmcalMaskedFastors");
  TObjArray *badchannelmap = static_cast<TObjArray *>(badchannelDB.GetObject(InputEvent()->GetRunNumber()));
  if(!badchannelmap || !badchannelmap->GetEntries()) return;
  for(TIter citer = TIter(badchannelmap).Begin(); citer != TIter::End(); ++citer){
    TParameter<int> *channelID = static_cast<TParameter<int> *>(*citer);
    AliDebugStream(1) << GetName() << ": Found masked fastor channel " << channelID->GetVal() << std::endl;
    fTriggerMaker->AddFastORBadChannel(channelID->GetVal());
  }
}


std::function<int (unsigned int, unsigned int)> AliEmcalTriggerMakerTask::GetMaskHandler() const {
  if(fGeom->GetTriggerMappingVersion() == 2){
    // Run 2 - complicated TRU layout in 6 subregions
    return [] (unsigned int ifield, unsigned int ibit) -> int {
      if(ifield >= 6 || ibit >= 16) throw kInvalidChannelException;
      const int kChannelMap[6][16] = {{ 8, 9,10,11,20,21,22,23,32,33,34,35,44,45,46,47},   // Channels in mask0
                                      {56,57,58,59,68,69,70,71,80,81,82,83,92,93,94,95},   // Channels in mask1
                                      { 4, 5, 6, 7,16,17,18,19,28,29,30,31,40,41,42,43},   // Channels in mask2
                                      {52,53,54,55,64,65,66,67,76,77,78,79,88,89,90,91},   // Channels in mask3
                                      { 0, 1, 2, 3,12,13,14,15,24,25,26,27,36,37,38,39},   // Channels in mask4
                                      {48,49,50,51,60,61,62,63,72,73,74,75,84,85,86,87}};  // Channels in mask5
      return kChannelMap[ifield][ibit];
    };
  } else {
    // Run 1 - linear mapping was used
    return [] (int ifield, int ibit) -> int {
      if(ifield >= 6 || ibit >= 16) throw kInvalidChannelException;
      return ifield * 16 + ibit;
    };
  }
}

void AliEmcalTriggerMakerTask::FillQAHistos(const TString &patchtype, const AliEMCALTriggerPatchInfo &recpatch){
  fQAHistos->FillTH2(Form("RCPos%s", patchtype.Data()), recpatch.GetColStart(), recpatch.GetRowStart());
  fQAHistos->FillTH2(Form("EPCentPos%s", patchtype.Data()), recpatch.GetEtaGeo(), recpatch.GetPhiGeo());
  fQAHistos->FillTH2(Form("PatchADCvsE%s", patchtype.Data()), recpatch.GetADCAmp(), recpatch.GetPatchE());
  fQAHistos->FillTH2(Form("PatchADCOffvsE%s", patchtype.Data()), recpatch.GetADCOfflineAmp(), recpatch.GetPatchE());
}
