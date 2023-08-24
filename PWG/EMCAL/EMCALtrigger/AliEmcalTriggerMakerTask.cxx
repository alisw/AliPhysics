/**************************************************************************************
 * Copyright (C) 2014, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#include <TClonesArray.h>
#include <TGrid.h>
#include <THashList.h>
#include <THistManager.h>
#include <TObjArray.h>
#include <TParameter.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,8,0)
#include <ROOT/TSeq.hxx>
#endif

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliDataFile.h"
#include "AliEmcalFastorMaskContainer.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerBitConfig.h"
#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerSTUDCSConfig.h"
#include "AliEMCALTriggerTRUDCSConfig.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerMakerKernel.h"
#include "AliEmcalTriggerMakerTask.h"
#include "AliEMCALTriggerMapping.h"
#include "AliEmcalTriggerMaskHandlerOCDB.h"
#include "AliLog.h"
#include "AliOADBContainer.h"

#include <array>
#include <bitset>
#include <cfloat>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

ClassImp(AliEmcalTriggerMakerTask)

AliEmcalTriggerMakerTask::AliEmcalTriggerMakerTask():
  AliAnalysisTaskEmcal(),
  fTriggerMaker(NULL),
  fV0(NULL),
  fCaloTriggersOutName("EmcalTriggers"),
  fV0InName("AliAODVZERO"),
  fBadFEEChannelOADB(""),
  fMaskedFastorOADB(""),
  fUseL0Amplitudes(kFALSE),
  fLoadFastORMaskingFromOCDB(kTRUE),
  fUseDeadFastORsOADB(kFALSE),
  fUseBadFastORsOADB(kFALSE),
  fCaloTriggersOut(0),
  fRunSmearing(kTRUE),
  fSimulateNoise(kTRUE),
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
  fMaskedFastorOADB(""),
  fUseL0Amplitudes(kFALSE),
  fLoadFastORMaskingFromOCDB(kTRUE),
  fUseDeadFastORsOADB(kFALSE),
  fUseBadFastORsOADB(kFALSE),
  fCaloTriggersOut(NULL),
  fRunSmearing(kTRUE),
  fSimulateNoise(kTRUE),
  fDoQA(doQA),
  fQAHistos(NULL)
{
  fTriggerMaker = new AliEmcalTriggerMakerKernel;
}

AliEmcalTriggerMakerTask::~AliEmcalTriggerMakerTask() {
  if(fTriggerMaker) delete fTriggerMaker;
  if(fQAHistos) delete fQAHistos;
  if(fCaloTriggersOut) delete fCaloTriggersOut;
}

/**
 * Initialize output objects
 */
void AliEmcalTriggerMakerTask::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  const std::array<const TString, 3>  kTriggerTypeNames = {{"EJE", "EGA", "EL0"}},
                                      kPatchTypes = {{"Online", "Offline", "Recalc"}};

  if(fDoQA) AliInfoStream() << "Trigger maker - QA requested" << std::endl;
  else AliInfoStream() << "Trigger maker - no QA requested" << std::endl;
  if(!fOutput) AliErrorStream() << "No output container initialized" << std::endl;

  if(fDoQA){
    if(fOutput){
      AliInfoStream() << "Enabling QA for trigger maker" << std::endl;
      fQAHistos = new THistManager("TriggerQA");

      for(const auto &triggertype : kTriggerTypeNames){
        for(const auto &patchtype : kPatchTypes){
          fQAHistos->CreateTH2(
              "RCPos" + triggertype + patchtype,
              "Lower edge position of " + patchtype + " " + triggertype + " patches (col-row);iEta;iPhi",
              48, -0.5, 47.5, 104, -0.5, 103.5
              );
          fQAHistos->CreateTH2(
              "EPCentPos" + triggertype + patchtype,
              "Center position of the " + patchtype + "  " + triggertype + " trigger patches;#eta;#phi",
              20, -0.8, 0.8, 700, 0., 7.
              );
          fQAHistos->CreateTH2(
              "PatchADCvsE" + triggertype + patchtype,
              "Patch ADC value for trigger type " + patchtype + " " + triggertype + "; Trigger ADC; FEE patch energy (GeV)",
              2000, 0., 2000, 200, 0., 200
              );
          fQAHistos->CreateTH2(
              "PatchEvsEsmear" + triggertype + patchtype,
              "Patch energy vs. smeared energy for " + patchtype + " " + triggertype +"; FEE patch energy (GeV); smeared FEE energy (GeV)",
              200, 0., 200, 200, 0., 200
          );
          fQAHistos->CreateTH2(
              "PatchADCvsEsmear" + triggertype + patchtype,
              "Patch ADC vs. smeared energy for " + patchtype + " " + triggertype + "; Trigger ADC; smeared FEE energy (GeV)",
              2000, 0., 2000, 200, 0., 200
          );
          fQAHistos->CreateTH2(
              "PatchADCOffvsE" + triggertype + patchtype,
              "Patch offline ADC value for trigger type " + patchtype + " " + triggertype + "; Trigger ADC; FEE patch energy (GeV)",
              2000, 0., 2000, 200, 0., 200
              );
          fQAHistos->CreateTH2(
              "PatchEvsADCrough" + triggertype + patchtype,
              "Patch Energy vs. ADC rough for trigger type " + patchtype + " " + triggertype + "; FEE patch energy (GeV); ADC rough (GeV)",
              2000, 0., 2000, 200, 0., 200
              );
        }
      }
      fQAHistos->CreateTH1("triggerBitsAll", "Trigger bits for all incoming patches;bit nr", 64, -0.5, 63.5);
      fQAHistos->CreateTH1("triggerBitsSel", "Trigger bits for reconstructed patches;bit nr", 64, -0.5, 63.5);
      fQAHistos->CreateTH2("FastORMaskOnline", "Masked FastORs at online level; col; row", 48, -0.5, 47.5, 104, -0.5, 103.5);

      // Monitoring of FastOR raw signals:
      fQAHistos->CreateTH2("FastORCorrEnergyADCrough", "FastOR cell energy vs. ADC energy",  200, 0., 20., 200, 0., 20.);
      fQAHistos->CreateTH2("FastORCorrEnergyESmear", "FastOR cell energy vs. smeared energy",  200, 0., 20., 200, 0., 20.);
      fQAHistos->CreateTH2("FastORCorrADCroughEsmear", "FastOR ADC rough vs. smeared energy", 200, 0., 20., 200, 0., 20.);
      fQAHistos->CreateTH2("FastORDiffEnergyADCrough", "FastOR ADCrough - cell energy", 4994, -0.5, 4993.5, 200, -10., 10);
      fQAHistos->CreateTH2("FastORDiffEnergyEsmear", "FastOR smeared energy - cell energy", 4994, -0.5, 4993.5, 200, -10., 10);
      fQAHistos->CreateTH2("FastORDiffEsmearADCrough", "FastOR ADC rough - smeared energy", 4994, -0.5, 4993.5, 200, -10., 10);

      for(auto h : *(fQAHistos->GetListOfHistograms())) fOutput->Add(h);
      fQAHistos->GetListOfHistograms()->SetOwner(false);
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

void AliEmcalTriggerMakerTask::ExecOnce(){
  AliAnalysisTaskEmcal::ExecOnce();

  if (!fLocalInitialized) return;

  if (!fCaloTriggersOutName.IsNull()) {
    fCaloTriggersOut = new TClonesArray("AliEMCALTriggerPatchInfo");
    fCaloTriggersOut->SetName(fCaloTriggersOutName);

    if (!(InputEvent()->FindListObject(fCaloTriggersOutName))) {
      InputEvent()->AddObject(fCaloTriggersOut);
    } else {
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
    std::cout << "EMCAL trigger maker: Not yet configure - automatically configuring ..." << std::endl;
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
    } else if((runnumber >= 224891 && runnumber <= 244628) || (runnumber >= 252235 && runnumber <= 294960)){
      // Configuration starting with LHC15f
      fTriggerMaker->ConfigureForPP2015();
      dataset = "pp 2015-2018";
      // Load additional masked FastORs from OADB (of course only if not specified explicitly, to allow for new maskings to be tested)
      if(!fMaskedFastorOADB.Length()) fMaskedFastorOADB = "oadb";
    } else if((runnumber >= 244824 && runnumber <= 246994) || (runnumber >= 295581)){
      fTriggerMaker->ConfigureForPbPb2015();
      dataset = "Pb-Pb Run2";
    }

    if(fTriggerMaker->IsConfigured()){
      std::cout << "EMCAL trigger maker: Applying configuration for " << dataset << std::endl;
    } else {
      std::cout << "EMCAL trigger maker: No valid configuration found for the given dataset - trigger maker run loop disabled" << std::endl;
    }

    if(fRunSmearing || fTriggerMaker->HasSmearModel()) {
      if(!fTriggerMaker->HasSmearModel()){
        std::cout << "EMCAL trigger maker: Smear mode - Initialize standard smear model" << std::endl;
        InitializeSmearModel(); // Initialize smear model if not yet set from outside
      }
      fRunSmearing = true;
      std::cout << "EMCAL trigger maker: Smear mode - require online bad channel map for smeared signal" << std::endl;
      fTriggerMaker->SetApplyOnlineBadChannelMaskingToSmeared();
    } 

    if(MCEvent() && fSimulateNoise) {
      // In MC mode add noise 
      // Using noise sigma of 50 MeV/channel as found during the optimization of the trigger efficiency to run2 data
      std::cout << "EMCAL trigger maker: Initialize standard noise model" << std::endl;
      if(!fTriggerMaker->HasNoiseModel()) fTriggerMaker->SetGaussianNoiseFEESmear(0., 0.05);
    }
    std::cout << "EMCAL trigger maker: Configured ..." << std::endl;
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
  if(!fTriggerMaker->IsConfigured()){
    AliErrorStream() << "Trigger maker not configured" << std::endl;
    return false;    // only run trigger maker in case it is configured
  }
  AliDebugStream(2) << "Looking for trigger patches ..." << std::endl;
  fCaloTriggersOut->Delete(); // Needed to avoid memory leak
  // prepare trigger maker
  fTriggerMaker->Reset();
  fTriggerMaker->ReadCellData(fCaloCells);
  fTriggerMaker->ReadTriggerData(fCaloTriggers);
  fTriggerMaker->BuildL1ThresholdsOffline(fV0);
  fTriggerMaker->SetIsMC(MCEvent());

  // QA on FastOR level (if enabled)
  if(fDoQA){
    double ecell, eadc, esmear;
    int fastORAbsID;
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,8,0)
    // beautiful ROOT6 version - with generators
    for(auto icol : ROOT::TSeqI(0, 48)){
      for(auto irow : ROOT::TSeqI(0, fTriggerMaker->GetDataGridDimensionRows())){
#else
    // stupid ROOT5-compatible version
    for(int icol = 0; icol < 48; icol++){
      for(int irow = 0; irow < fTriggerMaker->GetDataGridDimensionRows(); irow++){
#endif
        ecell = fTriggerMaker->GetTriggerChannelEnergy(icol, irow);
        eadc = fTriggerMaker->GetTriggerChannelEnergyRough(icol, irow);
        esmear = fTriggerMaker->GetTriggerChannelEnergySmeared(icol, irow);
        fGeom->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(icol, irow, fastORAbsID);
        // suppress pairs at (0,0)
        if(TMath::Abs(ecell) > DBL_EPSILON || TMath::Abs(eadc) > DBL_EPSILON){
          fQAHistos->FillTH2("FastORCorrEnergyADCrough", ecell, eadc);
          fQAHistos->FillTH2("FastORDiffEnergyADCrough", fastORAbsID, eadc - ecell);
        }
        if(TMath::Abs(ecell) > DBL_EPSILON || TMath::Abs(esmear) > DBL_EPSILON){
          fQAHistos->FillTH2("FastORCorrEnergyESmear", ecell, esmear);
          fQAHistos->FillTH2("FastORDiffEnergyEsmear", fastORAbsID, esmear - ecell);
        }
        if(TMath::Abs(eadc) > DBL_EPSILON || TMath::Abs(esmear) >  DBL_EPSILON){
          fQAHistos->FillTH2("FastORCorrADCroughEsmear", eadc, esmear);
          fQAHistos->FillTH2("FastORDiffEsmearADCrough", fastORAbsID, esmear - eadc);
        }
      }
    }
  }

  std::vector<AliEMCALTriggerPatchInfo> patches;
  fTriggerMaker->CreateTriggerPatches(InputEvent(), patches, fUseL0Amplitudes);
  Int_t patchcounter = 0;
  TString triggerstring;
  AliDebugStream(2) << GetName() << ": Found " << patches.size() << " patches" << std::endl;
  for(const auto &patchIter : patches){
    if(fDoQA){
      AliDebugStream(3) <<  GetName() << ": Next patch: size " << patchIter.GetPatchSize() << " , trigger bits " << std::bitset<sizeof(Int_t)*8>(patchIter.GetTriggerBits()) << std::endl;
      // Handle types different - online - offline - re
      if(patchIter.IsJetHigh() || patchIter.IsJetLow())                   FillQAHistos("EJEOnline", patchIter);
      if(patchIter.IsGammaHigh() || patchIter.IsGammaLow())               FillQAHistos("EGAOnline", patchIter);
      if(patchIter.IsJetHighSimple() || patchIter.IsJetLowSimple())       FillQAHistos("EJEOffline", patchIter);
      if(patchIter.IsGammaHighSimple() || patchIter.IsGammaLowSimple())   FillQAHistos("EGAOffline", patchIter);
      if(patchIter.IsLevel0())                                            FillQAHistos("EL0Online", patchIter);
      if(patchIter.IsLevel0Simple())                                      FillQAHistos("EL0Offline", patchIter);
      if(patchIter.IsLevel0Recalc())                                      FillQAHistos("EL0Recalc", patchIter);
      if(patchIter.IsRecalcJet())                                         FillQAHistos("EJERecalc", patchIter);
      if(patchIter.IsRecalcGamma())                                       FillQAHistos("EGARecalc", patchIter);
      // Redo checking of found trigger bits after masking of unwanted triggers
      int tBits = patchIter.GetTriggerBits();
      for(unsigned int ibit = 0; ibit < sizeof(tBits)*8; ibit++) {
        if(tBits & (1 << ibit)){
          fQAHistos->FillTH1("triggerBitsSel", ibit);
        }
      }
    }
    new((*fCaloTriggersOut)[patchcounter++]) AliEMCALTriggerPatchInfo(patchIter);
  }
  return true;
}

void AliEmcalTriggerMakerTask::RunChanged(Int_t newrun){
  std::cout << "EMCAL trigger maker: Run changed, new run " << newrun << ", loading new maskings ..." <<  std::endl;
  fTriggerMaker->ClearOfflineBadChannels();
  if(fBadFEEChannelOADB.Length()) InitializeBadFEEChannels();
  fTriggerMaker->ClearFastORBadChannels();
  if(fLoadFastORMaskingFromOCDB) InitializeFastORMaskingFromOCDB(newrun);
  if(fMaskedFastorOADB.Length() && (fUseDeadFastORsOADB || fUseBadFastORsOADB)) InitializeFastORMaskingFromOADB();
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
  std::cout << "EMCAL trigger maker: Loading additional bad FEE channels from OADB container " << fBadFEEChannelOADB << std::endl;
  fTriggerMaker->ClearOfflineBadChannels();
  if(fBadFEEChannelOADB.Contains("alien://") && !gGrid) TGrid::Connect("alien");
  AliOADBContainer badchannelDB("EmcalBadChannelsAdditional");
  badchannelDB.InitFromFile(fBadFEEChannelOADB.Data(), "EmcalBadChannelsAdditional");
  TObjArray *badchannelmap = static_cast<TObjArray *>(badchannelDB.GetObject(InputEvent()->GetRunNumber()));
  if(!badchannelmap || !badchannelmap->GetEntries()) return;
  for(TIter citer = TIter(badchannelmap).Begin(); citer != TIter::End(); ++citer){
    TParameter<int> *channelID = static_cast<TParameter<int> *>(*citer);
    fTriggerMaker->AddOfflineBadChannel(channelID->GetVal());
  }
}

void AliEmcalTriggerMakerTask::InitializeFastORMaskingFromOCDB(int runnumber){
  std::cout << "EMCAL trigger maker: Loading masked fastors from OCDB" << std::endl;
  auto channels = PWG::EMCAL::AliEmcalTriggerMaskHandlerOCDB::Instance()->GetMaskedFastorIndicesL1(runnumber);
  std::cout << "Found " << channels.size() << " masked FastORs at Level1" << std::endl;
  for(auto fastORAbsId : channels) {
    AliDebugStream(1) << "Adding masked FastOR " << fastORAbsId << " at L1" << std::endl;
    fTriggerMaker->AddFastORBadChannel(fastORAbsId);
  }
}

void AliEmcalTriggerMakerTask::InitializeSmearModel(){
  TF1 *meanmodel = new TF1("meanmodel", "pol1", 0., 1000.);
  meanmodel->SetParameter(0, -0.0206247);
  meanmodel->SetParameter(1, 0.966160);
  // Power law smearing
  TF1 *widthmodel = new TF1("widthmodel", "[0] * TMath::Power(x, [1]) + [2]", 0., 1000.);
  widthmodel->SetParameter(0, 0.0273139);
  widthmodel->SetParameter(1, 1.36187);
  widthmodel->SetParameter(2, 0.0736051);
  fTriggerMaker->SetSmearModel(meanmodel, widthmodel);
}

void AliEmcalTriggerMakerTask::InitializeFastORMaskingFromOADB(){
  std::cout << "EMCAL trigger maker: Loading masked fastors from OADB" << std::endl;
  TString containername;
  if(fMaskedFastorOADB == "oadb") {
    containername = AliDataFile::GetFileNameOADB("EMCAL/MaskedFastors.root").data();
  } else {
    containername = fMaskedFastorOADB;
  }
  std::cout << "EMCAL trigger maker: Initializing masked fastors from OADB container " << containername << std::endl;
  if(containername.Contains("alien://") && !gGrid) TGrid::Connect("alien");
  AliOADBContainer badchannelDB("AliEmcalMaskedFastors");
  badchannelDB.InitFromFile(containername, "AliEmcalMaskedFastors");
  PWG::EMCAL::AliEmcalFastorMaskContainer *maskContainer = dynamic_cast<PWG::EMCAL::AliEmcalFastorMaskContainer *>(badchannelDB.GetObject(InputEvent()->GetRunNumber()));
  if(maskContainer) {
    std::string selectiontype;
    std::vector<int> maskedfastors;
    if(fUseDeadFastORsOADB) {
      if(fUseBadFastORsOADB) {
        maskedfastors = maskContainer->GetMaskAll();
        selectiontype = "all";
      } else{
        maskedfastors = maskContainer->GetMaskDead();
        selectiontype = "dead";
      }
    } else if(fUseBadFastORsOADB){
      maskedfastors = maskContainer->GetMaskBad();
      selectiontype = "bad";
    } else selectiontype = "no";
    std::cout << "EMCAL trigger maker: OADB Container is of new type AliEmcalTriggerMaskContainer - masking " << selectiontype << " FastORs" << std::endl;
    for(auto fastor : maskedfastors) {
      AliDebugStream(1) << GetName() << ": Found masked fastor channel " << fastor << std::endl;
      fTriggerMaker->AddFastORBadChannel(fastor);
    }
  } else {
    TObjArray *badchannelmap = dynamic_cast<TObjArray *>(badchannelDB.GetObject(InputEvent()->GetRunNumber()));
    if(badchannelmap && badchannelmap->GetEntries()) {
      std::cout << "EMCAL trigger maker: OADB Container is of old type (simple list) - no distinction between bad and dead channels possible" << std::endl;
      for(TIter citer = TIter(badchannelmap).Begin(); citer != TIter::End(); ++citer){
        TParameter<int> *channelID = static_cast<TParameter<int> *>(*citer);
        AliDebugStream(1) << GetName() << ": Found masked fastor channel " << channelID->GetVal() << std::endl;
        fTriggerMaker->AddFastORBadChannel(channelID->GetVal());
      }
    } else {
      std::cerr << "EMCAL trigger maker: Unsupported OADB container type - no FastORs will be loaded" << std::endl;
    }
  }
}

void AliEmcalTriggerMakerTask::FillQAHistos(const TString &patchtype, const AliEMCALTriggerPatchInfo &recpatch){
  fQAHistos->FillTH2("RCPos" + patchtype, recpatch.GetColStart(), recpatch.GetRowStart());
  fQAHistos->FillTH2("EPCentPos" + patchtype, recpatch.GetEtaGeo(), recpatch.GetPhiGeo());
  fQAHistos->FillTH2("PatchADCvsE" + patchtype, recpatch.GetADCAmp(), recpatch.GetPatchE());
  fQAHistos->FillTH2("PatchEvsEsmear" + patchtype, recpatch.GetPatchE(), recpatch.GetSmearedEnergy());
  fQAHistos->FillTH2("PatchADCvsEsmear" + patchtype, recpatch.GetADCAmp(), recpatch.GetSmearedEnergy());
  fQAHistos->FillTH2("PatchADCOffvsE" + patchtype, recpatch.GetADCOfflineAmp(), recpatch.GetPatchE());
  fQAHistos->FillTH2("PatchEvsADCrough" + patchtype, recpatch.GetPatchE(), recpatch.GetADCAmpGeVRough());
}
