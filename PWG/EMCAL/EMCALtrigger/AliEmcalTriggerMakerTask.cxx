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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,8,0)
#include <ROOT/TSeq.hxx>
#endif

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

#include <array>
#include <bitset>
#include <cfloat>
#include <iostream>
#include <memory>
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
  fRunSmearing(kTRUE),
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
  fLoadFastORMaskingFromOCDB(kFALSE),
  fCaloTriggersOut(NULL),
  fRunSmearing(kTRUE),
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
    } else if((runnumber >= 224891 && runnumber <= 244628) || (runnumber >= 252235 && runnumber <= 294960)){
      // Configuration starting with LHC15f
      fTriggerMaker->ConfigureForPP2015();
      dataset = "pp 2015-2016";
    } else if((runnumber >= 244824 && runnumber <= 246994) || (runnumber >= 295581)){
      fTriggerMaker->ConfigureForPbPb2015();
      dataset = "Pb-Pb 2015";
    }

    if(fTriggerMaker->IsConfigured()){
      AliInfoStream() << "Applying configuration for " << dataset << std::endl;
    } else {
      AliErrorStream() << "No valid configuration found for the given dataset - trigger maker run loop disabled" << std::endl;
    }

    if(fRunSmearing && !fTriggerMaker->HasSmearModel()){
      InitializeSmearModel(); // Initialize smear model if not yet set from outside
      fTriggerMaker->SetApplyOnlineBadChannelMaskingToSmeared();
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
  if(!fTriggerMaker->IsConfigured()){
    AliErrorStream() << "Trigger maker not configured" << std::endl;
    return false;    // only run trigger maker in case it is configured
  }
  AliDebugStream(1) << "Looking for trigger patches ..." << std::endl;
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
  AliDebugStream(1) << "Run changed, new run " << newrun << std::endl;
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
            fGeom->GetTriggerMapping()->GetAbsFastORIndexFromTRU(RemapTRUIndex(itru), (ic =  GetMaskHandler(itru)(ifield, ibit)), fastOrAbsID);
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

void AliEmcalTriggerMakerTask::InitializeSmearModel(){
  std::cout << "Initializing trigger maker with default smearing parameterization" << std::endl;
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
  AliInfoStream() << "Initializing masked fastors from OADB container " << fMaskedFastorOADB.Data() << std::endl;
  if(fMaskedFastorOADB.Contains("alien://") && !gGrid) TGrid::Connect("alien");
  AliOADBContainer badchannelDB("AliEmcalMaskedFastors");
  badchannelDB.InitFromFile(fMaskedFastorOADB, "AliEmcalMaskedFastors");
  TObjArray *badchannelmap = static_cast<TObjArray *>(badchannelDB.GetObject(InputEvent()->GetRunNumber()));
  if(!badchannelmap || !badchannelmap->GetEntries()) return;
  for(TIter citer = TIter(badchannelmap).Begin(); citer != TIter::End(); ++citer){
    TParameter<int> *channelID = static_cast<TParameter<int> *>(*citer);
    AliDebugStream(1) << GetName() << ": Found masked fastor channel " << channelID->GetVal() << std::endl;
    fTriggerMaker->AddFastORBadChannel(channelID->GetVal());
  }
}


std::function<int (unsigned int, unsigned int)> AliEmcalTriggerMakerTask::GetMaskHandler(int itru) const {
  bool isTRUsmallSM = ((itru >= 30 && itru < 31) || (itru >= 44 && itru < 45)) ;
  if(fGeom->GetTriggerMappingVersion() == 2 && !isTRUsmallSM){
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

int AliEmcalTriggerMakerTask::RemapTRUIndex(int itru) const {
  if(fGeom->GetTriggerMappingVersion() == 2){
    const int trumapping[46] = {0,1,2,5,4,3,6,7,8,11,10,9,12,13,14,17,16,15,18,19,20,23,22,21,24,25,26,29,28,27,30,31,32,33,37,36,38,39,43,42,44,45,49,48,50,51};
    return trumapping[itru];
  } else return itru;
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
