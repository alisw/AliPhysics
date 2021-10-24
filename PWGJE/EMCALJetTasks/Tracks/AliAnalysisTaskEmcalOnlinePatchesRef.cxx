/************************************************************************************
 * Copyright (C) 2015, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <TClonesArray.h>
#include <THashList.h>
#include <TLinearBinning.h>
#include <THistManager.h>
#include <TParameter.h>
#include <TString.h>

#include "AliAnalysisManager.h"
#include "AliDataFile.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliInputEventHandler.h"
#include "AliOADBContainer.h"
#include "AliVEvent.h"
#include "AliVCaloCells.h"
#include "AliVVertex.h"

#include "AliAnalysisTaskEmcalOnlinePatchesRef.h"

#include <array>
#include <sstream>

using namespace PWGJE::EMCALJetTasks;

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalOnlinePatchesRef)

AliAnalysisTaskEmcalOnlinePatchesRef::AliAnalysisTaskEmcalOnlinePatchesRef():
  AliAnalysisTaskEmcal(),
  fHistos(nullptr),
  fFastOREnergy(nullptr),
  fFEEnergy(nullptr),
  fInOnlinePatch(nullptr),
  fMaskedCellsFastor(nullptr),
  fMaskedFastors(),
  fCellTimeCut(-1., 1.),
  fNameMaskedFastorOADB(),
  fNameMaskedCellOADB(AliDataFile::GetFileNameOADB("EMCAL/EMCALBadChannels.root").data()),
  fMaskedFastorOADB(nullptr),
  fMaskedCellOADB(nullptr),
  fRecoUtils(nullptr)
{
}

AliAnalysisTaskEmcalOnlinePatchesRef::AliAnalysisTaskEmcalOnlinePatchesRef(const char *name):
  AliAnalysisTaskEmcal(name, kTRUE),
  fHistos(nullptr),
  fFastOREnergy(nullptr),
  fFEEnergy(nullptr),
  fInOnlinePatch(nullptr),
  fMaskedCellsFastor(nullptr),
  fMaskedFastors(),
  fCellTimeCut(-1., 1.),
  fNameMaskedFastorOADB(),
  fNameMaskedCellOADB(AliDataFile::GetFileNameOADB("EMCAL/EMCALBadChannels.root").data()),
  fMaskedFastorOADB(nullptr),
  fMaskedCellOADB(nullptr),
  fRecoUtils(nullptr)
{
  SetNeedEmcalGeom(true);
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

AliAnalysisTaskEmcalOnlinePatchesRef::~AliAnalysisTaskEmcalOnlinePatchesRef() {
  if(fFastOREnergy) delete fFastOREnergy;
  if(fFEEnergy) delete fFEEnergy;
  if(fInOnlinePatch) delete fInOnlinePatch;
  if(fMaskedCellsFastor) delete fMaskedCellsFastor;
  if(fRecoUtils) delete fRecoUtils;
}

void AliAnalysisTaskEmcalOnlinePatchesRef::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  TLinearBinning patchenergybinning(300, 0., 300.),
                 fastorenergybinning(100, 0., 100.),
                 patchadcbinning(5000, 0., 5000.),
                 fastoradcbinning(2500, 0., 2500),
                 etabinning(200, -1., 1.),
                 phibinning(200, 0., TMath::TwoPi()),
                 nfastorbinning(1001, -0.5, 1000.5),
                 ncellbinning(5, -0.5, 4.5),
                 fastorIDbinning(7001, -0.5, 7000.5),
                 patchresidualsbinning(1000, -50., 50.),
                 fastorresidualsbinning(400, -20., 20.),
                 residualsbinningNormalized(2000, -10., 10.),
                 onlinepatchbinning(2, -0.5, 1.5);

  fHistos = new THistManager("Ref");
  fHistos->CreateTH1("hEventCount", "Event counter", 1, 0.5, 1.5);
  const TBinning *binningPatchEnergy[7] = {&patchenergybinning, &patchadcbinning, &patchenergybinning, &etabinning, &phibinning, &nfastorbinning, &nfastorbinning},
                 *binningFastorEnergy[6] = {&fastorIDbinning, &fastorenergybinning, &fastoradcbinning, &fastorenergybinning, &ncellbinning, &onlinepatchbinning},
                 *binningPatchResiduals[4] = {&patchenergybinning, &patchresidualsbinning, &nfastorbinning, &nfastorbinning},
                 *binningPatchResidualsNormalized[4] = {&patchenergybinning, &residualsbinningNormalized, &nfastorbinning, &nfastorbinning},
                 *binningFastorResiduals[6] = {&fastorIDbinning, &fastorenergybinning, &fastorenergybinning, &fastorresidualsbinning, &ncellbinning, &onlinepatchbinning};
  fHistos->CreateTHnSparse("hPatchEnergy", "Patch Energy; energy; ADC; Energy from FastOR; #eta; #phi; nFastorOnline; nFastorOffline", 7, binningPatchEnergy);
  fHistos->CreateTHnSparse("hFastorEnergy", "FastOR Energy; ID; energy; ADC; Energy from FastOR; ncell; online patch status", 6, binningFastorEnergy);
  fHistos->CreateTHnSparse("hPatchResiduals", "Patch energy residual binning; energy; residuals; nFastorsOnline; nFastorsOffline", 4, binningPatchResiduals);
  fHistos->CreateTHnSparse("hPatchResidualsNormalized", "Patch energy residual binning; energy; residuals; nFastorsOnline; nFastorsOffline", 4, binningPatchResidualsNormalized);
  fHistos->CreateTHnSparse("hFastorResiduals", "FastOR energy residuals; ID; energy; energy from FastOR; residuals; ncells; online patch status", 6, binningFastorResiduals);
  // Helper histograms checking the mask status of cells and FastORs
  fHistos->CreateTH1("hMaskedCells", "Index of masked cell; Cell index; Counts", 20001, -0.5, 20000.5);
  fHistos->CreateTH1("hMaskedFastors", "Index of masked FastOR; FastOR index; Counts", 7001, -0.5, 7000.5);
  for(auto hist : *fHistos->GetListOfHistograms()) fOutput->Add(hist);
  PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalOnlinePatchesRef::IsTriggerSelected(){ 
  TString triggerstring = fInputEvent->GetFiredTriggerClasses();
  std::map<std::string, AliVEvent::EOfflineTriggerTypes> triggertypes = {
    {"EMC7", AliVEvent::kEMC7}, {"DMC7", AliVEvent::kEMC7}, {"EGA", AliVEvent::kEMCEGA}, {"EJE", AliVEvent::kEMCEJE},
    {"EG1", AliVEvent::kEMCEGA}, {"EG2", AliVEvent::kEMCEGA}, {"DG1", AliVEvent::kEMCEGA}, {"DG2", AliVEvent::kEMCEGA},
    {"EJ1", AliVEvent::kEMCEJE}, {"EJ2", AliVEvent::kEMCEJE}, {"DJ1", AliVEvent::kEMCEJE}, {"DJ2", AliVEvent::kEMCEJE}
  };
  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  auto selectionbits = static_cast<UInt_t>(triggertypes.find(fOnlineTriggerClass.Data())->second);
  if(!((selectionstatus & selectionbits) && (triggerstring.Contains(fOnlineTriggerClass)))) return false;
  return true;
}

bool AliAnalysisTaskEmcalOnlinePatchesRef::Run(){
  fHistos->FillTH1("hEventCount", 1);

  // Load fastor and cell data
  fInOnlinePatch->Reset();
  LoadFastorEnergies();
  LoadCellEnergies();

  AliEMCALTriggerPatchInfo *mypatch(nullptr);
  for(auto patchiter = TIter(fTriggerPatchInfo).Begin(); patchiter != TIter::End(); ++patchiter){
    mypatch = dynamic_cast<AliEMCALTriggerPatchInfo *>(*patchiter);
    if(!mypatch) continue;
    if(!SelectPatch(*mypatch)) continue;
    MarkFastorsContributing(*mypatch);
    auto nfastor = GetNumberOfFastors(*mypatch);
    auto residuals = mypatch->GetADCAmpGeVRough() - mypatch->GetPatchE();
    auto residualsNormalized = residuals / mypatch->GetPatchE();
    Double_t pointEnergy[7] = {mypatch->GetPatchE(), static_cast<double>(mypatch->GetADCAmp()), mypatch->GetADCAmpGeVRough(), mypatch->GetEtaGeo(), mypatch->GetPhiGeo(), static_cast<double>(nfastor.first), static_cast<double>(nfastor.second)},
             pointResidual[4] = {mypatch->GetPatchE(), residuals, static_cast<double>(nfastor.first), static_cast<double>(nfastor.second)};
    fHistos->FillTHnSparse("hPatchEnergy", pointEnergy);
    fHistos->FillTHnSparse("hPatchResiduals", pointResidual);
    pointResidual[1] = residualsNormalized;
    fHistos->FillTHnSparse("hPatchResidualsNormalized", pointResidual);
  }
  Int_t fastorID;
  for(int icol = 0; icol < fInOnlinePatch->GetNumberOfCols(); icol++) {
    for(int irow = 0; irow < fInOnlinePatch->GetNumberOfRows(); irow++) {
      // Exclude PHOS region
      if(IsPhosHole(icol, irow)) continue;
      int onlinestatus = 0;
      if((*fInOnlinePatch)(icol, irow) > 0) onlinestatus = 1;
      fGeom->GetAbsFastORIndexFromPositionInEMCAL(icol,irow,fastorID);
      auto fastorenergy = (*fFastOREnergy)(icol, irow),
           feeenergy = (*fFEEnergy)(icol, irow),
           ncells = static_cast<double>((*fMaskedCellsFastor)(icol, irow));
      auto residuals = fastorenergy * EMCALTrigger::kEMCL1ADCtoGeV - feeenergy;
      Double_t pointenergy[6] = {static_cast<double>(fastorID), fastorenergy * EMCALTrigger::kEMCL1ADCtoGeV, fastorenergy, feeenergy, ncells, static_cast<double>(onlinestatus)},
               pointresiduals[6] = {static_cast<double>(fastorID), fastorenergy * EMCALTrigger::kEMCL1ADCtoGeV, feeenergy, residuals, ncells, static_cast<double>(onlinestatus)};
      fHistos->FillTHnSparse("hFastorEnergy", pointenergy);
      fHistos->FillTHnSparse("hFastorResiduals", pointresiduals);
    }
  }
  return true;
}

void AliAnalysisTaskEmcalOnlinePatchesRef::UserExecOnce(){
  int nrow = fGeom->GetTriggerMappingVersion() == 2 ? 104 : 64;
  fFEEnergy = new AliEMCALTriggerDataGrid<double>;
  fFEEnergy->Allocate(48, nrow);
  fFastOREnergy = new AliEMCALTriggerDataGrid<double>;
  fFastOREnergy->Allocate(48, nrow);
  fInOnlinePatch = new AliEMCALTriggerDataGrid<int>;
  fInOnlinePatch->Allocate(48, nrow);
  fMaskedCellsFastor = new AliEMCALTriggerDataGrid<int>;
  fMaskedCellsFastor->Allocate(48, nrow);

  if(fNameMaskedCellOADB.Length()){
    std::cout << "Reading bad channels from file " << fNameMaskedCellOADB << std::endl;
    fMaskedCellOADB = new AliOADBContainer("AliEMCALBadChannels");
    fMaskedCellOADB->InitFromFile(fNameMaskedCellOADB, "AliEMCALBadChannels");
    fRecoUtils = new AliEMCALRecoUtils;
  } else {
    AliErrorStream() << "No bad channels found" << fNameMaskedCellOADB << std::endl;
  }

  if(fNameMaskedFastorOADB.Length()){
    fMaskedFastorOADB = new AliOADBContainer("AliEmcalMaskedFastors");
    fMaskedFastorOADB->InitFromFile(fNameMaskedFastorOADB, "AliEmcalMaskedFastors");
  }
}

void AliAnalysisTaskEmcalOnlinePatchesRef::RunChanged(int newrun){
    // Load masked FastOR data
  if(fMaskedFastorOADB){
    AliInfoStream() << "Loading masked fastors for run " << newrun << std::endl;
    fMaskedFastors.clear();
    TObjArray *maskedfastors = static_cast<TObjArray *>(fMaskedFastorOADB->GetObject(newrun));
    if(maskedfastors && maskedfastors->GetEntries()){
      AliDebugStream(1) << "Loading masked fastOR container" << std::endl;
      for(auto masked : *maskedfastors){
        TParameter<int> *fastOrAbsID = static_cast<TParameter<int> *>(masked);
        fMaskedFastors.push_back(fastOrAbsID->GetVal());
        AliDebugStream(1) << "Masking fastor " << fastOrAbsID->GetVal() << std::endl;
        fHistos->FillTH1("hMaskedFastors", fastOrAbsID->GetVal());
      }
      std::sort(fMaskedFastors.begin(), fMaskedFastors.end(), std::less<int>());
    } else {
      AliErrorStream() << "No masked fastor container found " << std::endl;
    }
  }

  // Load masked cell data
  if(fMaskedCellOADB){
    fMaskedCellsFastor->Reset();
    AliInfoStream() << "Loading masked cells for run " << newrun << std::endl;
    fRecoUtils->SetEMCALChannelStatusMap(static_cast<TObjArray *>(fMaskedCellOADB->GetObject(newrun)));
    for(int icell = 0; icell < fGeom->GetNCells(); icell++) {
      if(IsCellMasked(icell)) {
        AliDebugStream(1) << "Masking cell " << icell << std::endl;
        fHistos->FillTH1("hMaskedCells", icell);
      }
    }

    // Fill masked cell data grid
    for(int icol = 0; icol < fMaskedCellsFastor->GetNumberOfCols(); icol++) {
      for(int irow = 0; irow < fMaskedCellsFastor->GetNumberOfRows(); irow++) {
        if(IsPhosHole(icol, irow)) continue; 
        int absFastorID;
        int cellIndices[4];
        fGeom->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(icol, irow, absFastorID);
        fGeom->GetTriggerMapping()->GetCellIndexFromFastORIndex(absFastorID, cellIndices);
        int nmasked = 0;
        for(int icell = 0; icell < 4; icell++) {
          if(IsCellMasked(cellIndices[icell])) nmasked++;
        }
        (*fMaskedCellsFastor)(icol, irow) = nmasked;
      }
    }
  }
}

void AliAnalysisTaskEmcalOnlinePatchesRef::MarkFastorsContributing(const AliEMCALTriggerPatchInfo &patch) const {
  int fastorID;
  for(auto icol = patch.GetColStart(); icol < patch.GetColStart() + patch.GetPatchSize(); icol++) {
    for(auto irow = patch.GetRowStart(); irow < patch.GetRowStart() + patch.GetPatchSize(); irow++) {
      fGeom->GetAbsFastORIndexFromPositionInEMCAL(icol, irow, fastorID);
      if(IsFastORMasked(fastorID)) continue;    // Don't mark masked FastORs
      (*fInOnlinePatch)(icol, irow) = 1;
    }
  }
}

std::pair<int,int> AliAnalysisTaskEmcalOnlinePatchesRef::GetNumberOfFastors(const AliEMCALTriggerPatchInfo &patch) const {
  int nfastorOnline = 0, nfastorOffline = 0;
  for(auto icol = patch.GetColStart(); icol < patch.GetColStart() + patch.GetPatchSize(); icol++) {
    for(auto irow = patch.GetRowStart(); irow < patch.GetRowStart() + patch.GetPatchSize(); irow++) {
      if((*fFastOREnergy)(icol, irow) > 0.) nfastorOnline++;
      if((*fFEEnergy)(icol, irow) > 0.) nfastorOffline++;
    }
  }
  return {nfastorOnline, nfastorOffline};
}

bool AliAnalysisTaskEmcalOnlinePatchesRef::SelectPatch(const AliEMCALTriggerPatchInfo &patch) const {
  if(!patch.IsOnline()) return false;
  if(fOnlineTriggerClass[0] == 'E') {
    if(patch.IsDCalPHOS()) return false;
  } else {
    if(patch.IsEMCal()) return false;
  }
  if((fOnlineTriggerClass == "EG1" || fOnlineTriggerClass == "DG1" || fOnlineTriggerClass == "EGA") && patch.IsGammaHigh()) return true; 
  if((fOnlineTriggerClass == "EG2" || fOnlineTriggerClass == "DG2") && patch.IsGammaLow()) return true; 
  if((fOnlineTriggerClass == "EJ1" || fOnlineTriggerClass == "DJ1" || fOnlineTriggerClass == "EJE") && patch.IsJetHigh()) return true; 
  if((fOnlineTriggerClass == "EJ2" || fOnlineTriggerClass == "DJ2") && patch.IsJetLow()) return true; 
  return false;
}

void AliAnalysisTaskEmcalOnlinePatchesRef::LoadCellEnergies(){
  fFEEnergy->Reset();
  AliVCaloCells *emccells = InputEvent()->GetEMCALCells();
  for(int icell = 0; icell < emccells->GetNumberOfCells(); icell++){
    int position = emccells->GetCellNumber(icell);
    double amplitude = emccells->GetAmplitude(icell);
    if(amplitude > 0){
      if(IsCellMasked(position)) {
        AliDebugStream(1) << "Non-0 cell energy " << amplitude << " found for masked cell " << position << std::endl;
      }
      int absFastor, col, row;
      fGeom->GetTriggerMapping()->GetFastORIndexFromCellIndex(position, absFastor);
      fGeom->GetPositionInEMCALFromAbsFastORIndex(absFastor, col, row);
      (*fFEEnergy)(col, row) += amplitude;
    }
  }
}

void AliAnalysisTaskEmcalOnlinePatchesRef::LoadFastorEnergies(){
  fFastOREnergy->Reset();
  auto triggers = fInputEvent->GetCaloTrigger("EMCAL");
  triggers->Reset();
  Int_t l1timesum, fastOrID, globCol, globRow;   
  while(triggers->Next()) {
    triggers->GetPosition(globCol, globRow);
    triggers->GetL1TimeSum(l1timesum);
    if(l1timesum <= 0) continue;
    fGeom->GetAbsFastORIndexFromPositionInEMCAL(globCol, globRow, fastOrID);
    if(IsFastORMasked(fastOrID)) {
      AliDebugStream(1) << "Non-0 fastor L1 energy " << l1timesum <<  " found for masked FastOR " << fastOrID << " (" << globCol << ", " << globRow << ")" << std::endl;
      continue;
    }
    (*fFastOREnergy)(globCol, globRow) += l1timesum;
  }
}

bool AliAnalysisTaskEmcalOnlinePatchesRef::IsPhosHole(int col, int row) const {
  if(row >= 64 && row < 100) {
    if(col >= 16 && col < 32) return true;
  }
  return false;
}

bool AliAnalysisTaskEmcalOnlinePatchesRef::IsCellMasked(int absCellID) const {
  if(!fRecoUtils) return false; // In case bad cells are not initialized declare cell as good
  Int_t smcell, modcell, colcell, rowcell, colcellsm, rowcellsm, channelstatus;
  fGeom->GetCellIndex(absCellID, smcell, modcell, rowcell, colcell);
  fGeom->GetCellPhiEtaIndexInSModule(smcell, modcell, rowcell, colcell, rowcellsm, colcellsm);
  return fRecoUtils->GetEMCALChannelStatus(smcell, colcellsm, rowcellsm, channelstatus);
}

bool AliAnalysisTaskEmcalOnlinePatchesRef::IsFastORMasked(int absFastORID) const {
  if(std::find(fMaskedFastors.begin(),fMaskedFastors.end(), absFastORID) != fMaskedFastors.end()) return true;
  return false;
}

AliAnalysisTaskEmcalOnlinePatchesRef *AliAnalysisTaskEmcalOnlinePatchesRef::AddTaskOnlinePatchesRef(const char *name, const char *suffix) {
  auto mgr = AliAnalysisManager::GetAnalysisManager();
  
  std::stringstream taskname, dirname, containername;
  taskname << name;
  dirname << mgr->GetCommonFileName() << ":OnlinePatchQA";
  containername << "OnlinePatchResults";
  if(strlen(suffix)) {
    taskname << "_" << suffix;
    dirname << "_" << suffix;
    containername << "_" << suffix;
  }

  auto task = new AliAnalysisTaskEmcalOnlinePatchesRef(taskname.str().data());
  mgr->AddTask(task);

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(containername.str().data(), TList::Class(), AliAnalysisManager::kOutputContainer, dirname.str().data()));

  return task;
}
