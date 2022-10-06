/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
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
#include <algorithm>
#include <iostream>
#include <vector>
#include <THashList.h>
#include <TH2.h>
#include <THistManager.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParameter.h>
#include <TVector3.h>

#include "AliAnalysisManager.h"
#include "AliEmcalFastOrMonitorTask.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALTriggerConstants.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliOADBContainer.h"
#include "AliVCaloCells.h"
#include "AliVCaloTrigger.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliDataFile.h"

ClassImp(PWG::EMCAL::AliEmcalFastOrMonitorTask)

using namespace PWG::EMCAL;

AliEmcalFastOrMonitorTask::AliEmcalFastOrMonitorTask() :
  AliAnalysisTaskEmcal(),
  fHistosQA(nullptr),
  fLocalInitialized(false),
  fOldRun(-1),
  fRequestTrigger(AliVEvent::kAny),
  fTriggerPattern(""),
  fCellData(),
  fMaskedFastors(),
  fMinCellTimeNS(-1e9),
  fMaxCellTimeNS(1e9),
  fNameMaskedFastorOADB(),
  fNameMaskedCellOADB(AliDataFile::GetFileNameOADB("EMCAL/EMCALBadChannels_1D.root").data()),
  fMaskedFastorOADB(nullptr),
  fMaskedCellOADB(nullptr),
  fRecoUtils(nullptr)
{
  SetNeedEmcalGeom(kTRUE);
}

AliEmcalFastOrMonitorTask::AliEmcalFastOrMonitorTask(const char *name) :
  AliAnalysisTaskEmcal(name, kTRUE),
  fHistosQA(nullptr),
  fLocalInitialized(false),
  fOldRun(-1),
  fRequestTrigger(AliVEvent::kAny),
  fTriggerPattern(""),
  fMaskedFastors(),
  fMinCellTimeNS(-1e9),
  fMaxCellTimeNS(1e9),
  fNameMaskedFastorOADB(),
  fNameMaskedCellOADB(AliDataFile::GetFileNameOADB("EMCAL/EMCALBadChannels_1D.root").data()),
  fMaskedFastorOADB(nullptr),
  fMaskedCellOADB(nullptr),
  fRecoUtils(nullptr)
{
  SetNeedEmcalGeom(kTRUE);
  DefineOutput(1, TList::Class());
}

AliEmcalFastOrMonitorTask::~AliEmcalFastOrMonitorTask() {
  if(fMaskedFastorOADB) delete fMaskedFastorOADB;
  if(fMaskedCellOADB) delete fMaskedCellOADB;
  if(fRecoUtils) delete fRecoUtils;
}

void AliEmcalFastOrMonitorTask::UserCreateOutputObjects() {
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  TString cellsName;
  auto evhand = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (evhand->InheritsFrom("AliESDInputHandler")) {
    cellsName = "EMCALCells";
    AliInfoStream() << "ESD analysis, cellsName = \"" << cellsName << "\"" << std::endl;
  } else {
    cellsName = "emcalCells";
    AliInfoStream() << "AOD analysis, cellsName = \"" << cellsName << "\"" << std::endl;
  }
  SetCaloCellsName(cellsName.Data());

  fHistosQA = new THistManager("fastOrHistos");

  const int kMaxCol = 48, kMaxRow = 104, kMaxFastOr = kMaxRow * kMaxCol;

  fHistosQA->CreateTH1("hEvents", "Number of events", 1, 0.5, 1.5);
  fHistosQA->CreateTH1("hFastOrFrequencyL0", "FastOr frequency at Level0", kMaxFastOr, -0.5, kMaxFastOr - 0.5);
  fHistosQA->CreateTH1("hFastOrFrequencyL1", "FastOr frequency at Level1", kMaxFastOr, -0.5, kMaxFastOr - 0.5);
  fHistosQA->CreateTH2("hFastOrAmplitude", "FastOr amplitudes; FastOR Abs ID; Amplitude", kMaxFastOr, -0.5, kMaxFastOr - 0.5, 513, -0.5, 512.5);
  fHistosQA->CreateTH2("hFastOrTimeSum", "FastOr time sum; FastOR Abs ID; L0 time sum", kMaxFastOr, -0.5, kMaxFastOr - 0.5, 2049, -0.5, 2048.5);
  fHistosQA->CreateTH2("hFastOrTransverseTimeSum", "FastOr transverse time sum; FastOR Abs ID; L0 transverse time sum", kMaxFastOr, -0.5, kMaxFastOr - 0.5, 2049, -0.5, 2048.5);
  fHistosQA->CreateTH2("hFastOrNL0Times", "FastOr Number of L0 times; FastOR AbsID; Number of L0 times", kMaxFastOr, -0.5, kMaxFastOr - 0.5, 16, -0.5, 15.5);
  fHistosQA->CreateTH2("hFastOrColRowFrequencyL0", "FastOr Frequency (col-row) at Level0; col; row", kMaxCol, -0.5, kMaxCol - 0.5, kMaxRow, -0.5, kMaxRow - 0.5);
  fHistosQA->CreateTH2("hFastOrColRowFrequencyL1", "FastOr Frequency (col-row) at Level1; col; row", kMaxCol, -0.5, kMaxCol - 0.5, kMaxRow, -0.5, kMaxRow - 0.5);
  fHistosQA->CreateTH2("hEnergyFastorCell", "Sum of cell energy vs. fastor Energy", 1000, 0., 20., 1000 , 0., 20.);
  fHistosQA->CreateTH2("hEnergyFastorCellL0", "Sum of cell energy vs. fastor Energy (L0)", 1000, 0., 20., 1000 , 0., 20.);
  fHistosQA->CreateTH2("hEnergyFastorCellL0Amp", "Sum of cell energy vs. fastor Amplitude (L0)", 1000, 0., 20., 1000 , 0., 1000.);

  // Helper histograms checking the mask status of cells and FastORs
  fHistosQA->CreateTH1("hMaskedFastors", "Index of masked FastOR; FastOR index; Counts", 7001, -0.5, 7000.5);
  fHistosQA->CreateTH1("hMaskedCells", "Index of masked FastOR; FastOR index; Counts", 20001, -0.5, 20000.5);
  fHistosQA->CreateTH1("hCellEnergyCount", "Counts of non-0 cell entries; Cell index; Counts", 20001, -0.5, 20000.5);
  fHistosQA->CreateTH2("hCellTimeBefore", "Cell time before time cut", 2000, -1000., 1000., 200., 0., 200.);
  fHistosQA->CreateTH2("hCellTimeAfter", "Cell time after time cut", 2000, -1000., 1000., 200., 0., 200.);

  // THnSparse for fastor-by-fastor energy decalibration
  TAxis fastorIDAxis(4992, -0.5, 4991.5), offlineaxis(200, 0., 20.), onlineaxis(200, 0., 20.), residualaxis(1000, -10., 10.), cellmaskaxis(5, -0.5, 4.5);
  const TAxis *sparseaxis[5] = {&fastorIDAxis, &offlineaxis, &onlineaxis, &residualaxis, &cellmaskaxis};
  fastorIDAxis.SetNameTitle("fastorAbsID", "FastOR abs. ID");
  offlineaxis.SetNameTitle("offlinenergy", "E_{2x2 cells} (GeV)");
  onlineaxis.SetNameTitle("onlineenergy", "E_{FastOR} (GeV)");
  residualaxis.SetNameTitle("residuals", "E_{FastOR} - E_{2x2 cells} (GeV)");
  cellmaskaxis.SetNameTitle("maskedcells", "Number of masked cells");
  fHistosQA->CreateTHnSparse("hFastOrEnergyOfflineOnline", "FastOr Offline vs Online energy", 5, sparseaxis);
  fHistosQA->CreateTHnSparse("hFastOrEnergyOfflineOnlineL0", "FastOr Offline vs Online energy (L0)", 5, sparseaxis);

  TAxis adcAxisL0(2101, -0.5, 2100.5), adcAxisL1(2101, -0.5, 2100.5);
  const TAxis *adcaxes[3] = {&fastorIDAxis, &adcAxisL0, &adcAxisL1};
  adcAxisL0.SetNameTitle("adcL0", "L0 ADC amplitude");
  adcAxisL1.SetNameTitle("adcL1", "L1 ADC amplitude");
  fHistosQA->CreateTHnSparse("hFastOrADCL0L1", "FastOR ADC in L0 and L1", 3, adcaxes);


  for(auto h : *(fHistosQA->GetListOfHistograms())) fOutput->Add(h);
  PostData(1, fOutput);
}

void AliEmcalFastOrMonitorTask::UserExecOnce(){
  int nrow = fGeom->GetTriggerMappingVersion() == 2 ? 104 : 64;
  fCellData.Allocate(48, nrow);

  if(fNameMaskedCellOADB.Length()){
    fMaskedCellOADB = new AliOADBContainer("AliEMCALBadChannels");
    fMaskedCellOADB->InitFromFile(fNameMaskedCellOADB, "AliEMCALBadChannels");

    fRecoUtils = new AliEMCALRecoUtils;
  }

  if(fNameMaskedFastorOADB.Length()){
    fMaskedFastorOADB = new AliOADBContainer("AliEmcalMaskedFastors");
    fMaskedFastorOADB->InitFromFile(fNameMaskedFastorOADB, "AliEmcalMaskedFastors");
  }
}

void AliEmcalFastOrMonitorTask::RunChanged(Int_t newrun){
  // Load masked FastOR data
  if(fMaskedFastorOADB){
    AliInfoStream() << "Loading masked FastORs for run " << newrun << std::endl;
    fMaskedFastors.clear();
    TObjArray *maskedfastors = static_cast<TObjArray *>(fMaskedFastorOADB->GetObject(newrun));
    if(maskedfastors && maskedfastors->GetEntries()){
      for(auto masked : *maskedfastors){
        TParameter<int> *fastOrAbsID = static_cast<TParameter<int> *>(masked);
        fMaskedFastors.push_back(fastOrAbsID->GetVal());
        fHistosQA->FillTH1("hMaskedFastors", fastOrAbsID->GetVal());
      }
      std::sort(fMaskedFastors.begin(), fMaskedFastors.end(), std::less<int>());
    }
  }

  // Load masked cell data
  if(fMaskedCellOADB){
    AliInfoStream() << "Loading masked cells for run " << newrun << std::endl;
    TObjArray *maskedcells = static_cast<TObjArray *>(fMaskedCellOADB->GetObject(newrun));
    fRecoUtils->SetEMCALChannelStatusMap1D(static_cast<TH1C*>(maskedcells->FindObject("EMCALBadChannelMap")));
    Int_t cellstatus;
    for(int icell = 0; icell < fGeom->GetNCells(); icell++) {
      auto masked = fRecoUtils->GetEMCALChannelStatus1D(icell, cellstatus);
      if(masked) fHistosQA->FillTH1("hMaskedCells", icell);
    }
  }
}

bool AliEmcalFastOrMonitorTask::IsEventSelected() {
  // Check trigger
  if(!(fInputHandler->IsEventSelected() & fRequestTrigger)) return false;
  if(fTriggerPattern.Length()){
    if(!TString(InputEvent()->GetFiredTriggerClasses()).Contains(fTriggerPattern)) return false;
  }
  return true;
}

bool AliEmcalFastOrMonitorTask::Run() {
  const Double_t kEMCL0ADCtoGeV = 0.018970588*4;
  LoadEventCellData();
  fHistosQA->FillTH1("hEvents", 1);
  AliVCaloTrigger *triggerdata = InputEvent()->GetCaloTrigger("EMCAL");
  triggerdata->Reset();
  Int_t nl0times, l1timesum, fastOrID, globCol, globRow;
  Float_t amp;
  while(triggerdata->Next()){
    triggerdata->GetAmplitude(amp);
    triggerdata->GetNL0Times(nl0times);
    triggerdata->GetL1TimeSum(l1timesum);
    triggerdata->GetPosition(globCol, globRow);
    fGeom->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(globCol, globRow, fastOrID);
    if(amp > 1e-5){
      fHistosQA->FillTH2("hFastOrColRowFrequencyL0", globCol, globRow);
      fHistosQA->FillTH1("hFastOrFrequencyL0", fastOrID);
    }
    if(l1timesum){
      fHistosQA->FillTH2("hFastOrColRowFrequencyL1", globCol, globRow);
      fHistosQA->FillTH1("hFastOrFrequencyL1", fastOrID);
    }
    if(std::find(fMaskedFastors.begin(), fMaskedFastors.end(), fastOrID) == fMaskedFastors.end()){
      fHistosQA->FillTH2("hFastOrAmplitude", fastOrID, amp);
      fHistosQA->FillTH2("hFastOrTimeSum", fastOrID, l1timesum);
      fHistosQA->FillTH2("hFastOrNL0Times", fastOrID, nl0times);
      fHistosQA->FillTH2("hFastOrTransverseTimeSum", fastOrID, GetTransverseTimeSum(fastOrID, l1timesum, fVertex));
      fHistosQA->FillTH2("hEnergyFastorCell", fCellData(globCol, globRow), l1timesum * EMCALTrigger::kEMCL1ADCtoGeV);
      fHistosQA->FillTH2("hEnergyFastorCellL0", fCellData(globCol, globRow), amp * kEMCL0ADCtoGeV);
      fHistosQA->FillTH2("hEnergyFastorCellL0Amp", fCellData(globCol, globRow), amp);
      int ncellmasked = 0;
      if(fRecoUtils){
        int fastorCells[4];
        fGeom->GetTriggerMapping()->GetCellIndexFromFastORIndex(fastOrID, fastorCells);
        for(int icell = 0; icell < 4; icell++){
          if(IsCellMasked(fastorCells[icell])) ncellmasked++;
        }
      }
      double  feeenergy = fCellData(globCol, globRow),
              l0energy = amp * kEMCL0ADCtoGeV,
              l1energy = static_cast<double>(l1timesum) * EMCALTrigger::kEMCL1ADCtoGeV;
      double energydata[5] = {
            static_cast<double>(fastOrID),
            feeenergy,
            l1energy,
            l1energy - feeenergy,
            static_cast<double>(ncellmasked)
      };
      fHistosQA->FillTHnSparse("hFastOrEnergyOfflineOnline", energydata);
      energydata[2] = l0energy;
      energydata[3] = l0energy - feeenergy;
      fHistosQA->FillTHnSparse("hFastOrEnergyOfflineOnlineL0", energydata);
      double adcdata[3] = {
            static_cast<double>(fastOrID),
            amp * 4.,         // correction for the bit shift
            static_cast<double>(l1timesum)
      };
      fHistosQA->FillTHnSparse("hFastOrADCL0L1", adcdata);
    }
  }
  return true;
}

void AliEmcalFastOrMonitorTask::LoadEventCellData(){
  const Double_t kSecToNanoSec = 1e9;
  fCellData.Reset();
  for(int icell = 0; icell < fCaloCells->GetNumberOfCells(); icell++){
    int position = fCaloCells->GetCellNumber(icell);
    double amplitude = fCaloCells->GetAmplitude(icell),
           celltimeNS = fCaloCells->GetTime(icell) * kSecToNanoSec;
    if ( position < 0 ) {
      AliErrorStream() << "Skip: iCell = "<< icell << "Cell Id = " << position <<", E = "<< amplitude <<", time = "<<celltimeNS<<std::endl;
      continue;
    }
    if(amplitude > 0) fHistosQA->FillTH2("hCellTimeBefore", celltimeNS, amplitude);
    if(celltimeNS < fMinCellTimeNS || celltimeNS > fMaxCellTimeNS) continue;
    if(amplitude > 0){
      AliDebugStream(1) << "Found cell time " << celltimeNS << " nanosec" << std::endl;
      fHistosQA->FillTH2("hCellTimeAfter", celltimeNS, amplitude);
      fHistosQA->FillTH1("hCellEnergyCount", position);
      int absFastor, col, row;
      fGeom->GetTriggerMapping()->GetFastORIndexFromCellIndex(position, absFastor);
      fGeom->GetPositionInEMCALFromAbsFastORIndex(absFastor, col, row);
      fCellData(col, row) += amplitude;
    }
  }
}

Double_t AliEmcalFastOrMonitorTask::GetTransverseTimeSum(Int_t fastorAbsID, Double_t adc, const Double_t *vertex) const{
  Int_t cellIDs[4];
  fGeom->GetTriggerMapping()->GetCellIndexFromFastORIndex(fastorAbsID, cellIDs);
  std::vector<double> eta, phi;
  for(int i = 0l; i < 4; i++){
    double etatmp, phitmp;
    fGeom->EtaPhiFromIndex(cellIDs[i], etatmp, phitmp);
    eta.push_back(etatmp);
    phi.push_back(phitmp);
  }

  // Calculate FastOR position: for the approximation take mean eta and phi
  // Radius is taken from the geometry
  TVector3 fastorPos, vertexPos(vertex[0], vertex[1], vertex[2]);
  fastorPos.SetPtEtaPhi(fGeom->GetIPDistance(), TMath::Mean(eta.begin(), eta.end()), TMath::Mean(phi.begin(), phi.end()));
  fastorPos -= vertexPos;

  TLorentzVector evec(fastorPos, adc);
  return evec.Et();
}

bool AliEmcalFastOrMonitorTask::IsCellMasked(int absCellID) const {
  if(!fRecoUtils) return false; // In case bad cells are not initialized declare cell as good
  Int_t channelstatus;
  return fRecoUtils->GetEMCALChannelStatus1D(absCellID, channelstatus);
}
