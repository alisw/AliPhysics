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
#include <map>
#include <vector>
#include <TArrayD.h>
#include <TClonesArray.h>
#include <TGrid.h>
#include <THashList.h>
#include <THistManager.h>
#include <TLinearBinning.h>
#include <TObjArray.h>
#include <TParameter.h>

#include "AliEMCALGeometry.h"
#include "AliEmcalCellMonitorTask.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliOADBContainer.h"
#include "AliVCaloCells.h"
#include "AliVEvent.h"

/// \cond CLASSIMP
ClassImp(PWG::EMCAL::AliEmcalCellMonitorTask)
/// \endcond

using namespace PWG::EMCAL;

AliEmcalCellMonitorTask::AliEmcalCellMonitorTask() :
   AliAnalysisTaskSE(),
   fLocalInitialized(kFALSE),
   fHistManager(nullptr),
   fGeometry(nullptr),
   fMinCellAmplitude(0),
   fRequestTrigger(AliVEvent::kAnyINT),
   fBadChannelContainer(""),
   fTriggerString(""),
   fNumberOfCells(12288),
   fOldRun(-1),
   fMaskedCells()
{

}

AliEmcalCellMonitorTask::AliEmcalCellMonitorTask(const char *name) :
   AliAnalysisTaskSE(name),
   fLocalInitialized(kFALSE),
   fHistManager(nullptr),
   fGeometry(nullptr),
   fMinCellAmplitude(0),
   fRequestTrigger(AliVEvent::kAnyINT),
   fBadChannelContainer(""),
   fTriggerString(""),
   fNumberOfCells(12288),
   fOldRun(-1),
   fMaskedCells()
{
  DefineOutput(1, TList::Class());
}

AliEmcalCellMonitorTask::~AliEmcalCellMonitorTask() {
  if(fGeometry) delete fGeometry;
}

void AliEmcalCellMonitorTask::UserCreateOutputObjects(){
  fHistManager = new THistManager("EMCALCellMonitor");


  PostData(1, fHistManager->GetListOfHistograms());
}

void AliEmcalCellMonitorTask::ExecOnce(){
  if(!fGeometry)
    fGeometry = AliEMCALGeometry::GetInstanceFromRunNumber(fInputEvent->GetRunNumber());
  fNumberOfCells = fGeometry->GetNCells();
  CreateHistograms();
}

void AliEmcalCellMonitorTask::RunChanged(){
  if(fBadChannelContainer.Length()) LoadCellMasking();
  for(auto cellID : fMaskedCells) fHistManager->FillTH1("cellMasking", cellID);
}

void AliEmcalCellMonitorTask::UserExec(Option_t *){
  if(!fLocalInitialized) {
    ExecOnce();
    fLocalInitialized = kTRUE;
  }

  // Run change
  if(InputEvent()->GetRunNumber() != fOldRun){
    RunChanged();
    fOldRun = InputEvent()->GetRunNumber();
  }

  // Check trigger
  if(!(fInputHandler->IsEventSelected() & fRequestTrigger)) return;
  if(fTriggerString.Length()){
    if(!TString(InputEvent()->GetFiredTriggerClasses()).Contains(fTriggerString)) return;
  }

  fHistManager->FillTH1("events", 1);

  AliVCaloCells *emcalcells = fInputEvent->GetEMCALCells();

  // input data
  Short_t cellNumber;
  Double_t amplitude, celltime, efrac;
  Int_t mclabel;

  Int_t sm, mod, meta, mphi, ieta, iphi;
  for(int icell = 0; icell < emcalcells->GetNumberOfCells(); icell++){
    emcalcells->GetCell(icell, cellNumber, amplitude, celltime, mclabel, efrac);
    if(IsCellMasked(cellNumber)) continue;
    fHistManager->FillTH2("cellAmplitude", amplitude, cellNumber);
    if(amplitude < fMinCellAmplitude) continue;
    fHistManager->FillTH1("cellAmplitudeCut", amplitude, cellNumber);
    fHistManager->FillTH1("cellFrequency", cellNumber);
    fHistManager->FillTH2("cellTime", celltime, cellNumber);
    if(celltime >= 1e-6) fHistManager->FillTH2("cellTimeOutlier", celltime, cellNumber);
    if(celltime > -5e-8 && celltime < 1e-7) fHistManager->FillTH2("cellTimeMain", celltime, cellNumber);

    // Get Cell index in eta-phi of sm
    fGeometry->GetCellIndex(cellNumber, sm, mod, mphi, meta);
    fGeometry->GetCellPhiEtaIndexInSModule(sm, mod, mphi, meta, iphi, ieta);

    fHistManager->FillTH2(Form("cellCountSM%d", sm), ieta, iphi);
    fHistManager->FillTH2(Form("cellAmpSM%d", sm), ieta, iphi, amplitude);
    fHistManager->FillTH2(Form("cellAmpTimeCorrSM%d", sm), celltime, amplitude);
  }

  // Cluster loop
  if(fNameClusters.Length()){
    TClonesArray *clustercont = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject(fNameClusters.Data()));
    if(clustercont){
      const AliVCluster *myclust = nullptr;
      for(TIter clusteriter = TIter(clustercont).Begin(); clusteriter != TIter::End(); ++clusteriter){
        myclust = dynamic_cast<const AliVCluster *>(*clusteriter);
        if(!myclust) continue;
        for(int icell = 0; icell < myclust->GetNCells(); icell++){
          fHistManager->FillTH1("cellClusterOccurrency", myclust->GetCellAbsId(icell));
          fHistManager->FillTH2("cellAmplitudeFractionCluster", myclust->GetCellAbsId(icell), myclust->GetCellAmplitudeFraction(icell));
        }
      }
    } else {
      AliErrorStream() << GetName() << ": cluster container " << fNameClusters << " not found in the input event" << std::endl;
    }
  }
  PostData(1, fHistManager->GetListOfHistograms());
}

void AliEmcalCellMonitorTask::CreateHistograms(){
  fHistManager->CreateTH1("events", "Number of events", 1, 0.5, 1.5);
  fHistManager->CreateTH1("cellMasking", "Monitoring for masked cells", TLinearBinning(fNumberOfCells, -0.5, fNumberOfCells - 0.5));
  fHistManager->CreateTH1("cellFrequency", "Frequency of cell firing", TLinearBinning(fNumberOfCells, -0.5, fNumberOfCells - 0.5));
  fHistManager->CreateTH2("cellAmplitude", "Energy distribution per cell", AliEmcalCellMonitorAmplitudeBinning(), TLinearBinning(fNumberOfCells, -0.5, fNumberOfCells - 0.5));
  fHistManager->CreateTH2("cellAmplitudeCut", "Energy distribution per cell (after energy cut)", AliEmcalCellMonitorAmplitudeBinning(), TLinearBinning(fNumberOfCells, -0.5, fNumberOfCells - 0.5));
  fHistManager->CreateTH2("cellTime", "Time distribution per cell", 300, -3e-7, 1e-6, fNumberOfCells, -0.5, fNumberOfCells - 0.5);
  fHistManager->CreateTH2("cellTimeOutlier", "Outlier time distribution per cell", 100, 1e-6, 5e-5, fNumberOfCells, -0.5, fNumberOfCells - 0.5);
  fHistManager->CreateTH2("cellTimeMain", "Time distribution per cell for the main bunch", 150, -50e-9, 100e-9, fNumberOfCells, -0.5, fNumberOfCells - 0.5);
  fHistManager->CreateTH1("cellClusterOccurrency", "Occurrency of a cell in clusters", fNumberOfCells, -0.5, fNumberOfCells - 0.5);
  fHistManager->CreateTH2("cellAmplitudeFractionCluster", "Summed cell amplitude fraction in a cluster", fNumberOfCells, -0.5, fNumberOfCells - 0.5, 200, 0., 200.);
  for(int ism = 0; ism < 20; ++ism){
    fHistManager->CreateTH2(Form("cellAmpSM%d", ism), Form("Integrated cell amplitudes for SM %d; col; row", ism), 48, -0.5, 47.5, 24, -0.5, 23.5);
    fHistManager->CreateTH2(Form("cellCountSM%d", ism), Form("Count rate per cell for SM %d; col; row", ism), 48, -0.5, 47.5, 24, -0.5, 23.5);
  }

  for(int ism = 0; ism < 20; ++ism){
    fHistManager->CreateTH2(Form("cellAmpTimeCorrSM%d", ism), Form("Correlation between cell amplitude and time in Supermodule %d", ism), 1000, -5e-7, 5e-7, 1000, 0., 100.);
  }
}

void AliEmcalCellMonitorTask::LoadCellMasking(){
  if(!fBadChannelContainer.Length()) return;
  AliInfoStream() << GetName() << ": Loading bad channel map from " <<fBadChannelContainer << std::endl;
  fMaskedCells.clear();
  if(fBadChannelContainer.Contains("alien://") && ! gGrid) TGrid::Connect("alien://");    // Make sure alien connection is available as the AliOADBContainer doesn't handle it
  AliOADBContainer contreader("EmcalBadChannelsAdditional");
  contreader.InitFromFile(fBadChannelContainer.Data(), "EmcalBadChannelsAdditional");
  TObjArray *rundata = dynamic_cast<TObjArray *>(contreader.GetObject(InputEvent()->GetRunNumber()));
  if(!rundata) return;
  for(TIter channeliter = TIter(rundata).Begin(); channeliter != TIter::End(); ++channeliter){
    TParameter<int> *cellID = static_cast<TParameter<int> *>(*channeliter);
    if(cellID) SetBadCell(cellID->GetVal());
  }
}

void AliEmcalCellMonitorTask::SetBadCell(Int_t cellId){
  if(std::find(fMaskedCells.begin(), fMaskedCells.end(), cellId) != fMaskedCells.end()) return;
  fMaskedCells.push_back(cellId);
}

bool AliEmcalCellMonitorTask::IsCellMasked(Int_t cellId) const {
  return (std::find(fMaskedCells.begin(), fMaskedCells.end(), cellId) != fMaskedCells.end());
}

AliEmcalCellMonitorTask::AliEmcalCellMonitorAmplitudeBinning::AliEmcalCellMonitorAmplitudeBinning():
    TCustomBinning()
{
  SetMinimum(0);
  AddStep(2., 0.1);
  AddStep(5., 0.2);
  AddStep(10., 0.5);
  AddStep(20., 1.);
  AddStep(50., 2.);
  AddStep(100., 5.);
  AddStep(200., 10.);
}
