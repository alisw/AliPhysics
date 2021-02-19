/************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                 *
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
#include <sstream>
#include <string>

#include <TH2.h>
#include <TList.h>
#include <TNtuple.h>

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalResponseOutliers.h"
#include "AliClusterContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliTrackContainer.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalResponseOutliers)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalResponseOutliers::AliAnalysisTaskEmcalResponseOutliers():
    AliAnalysisTaskEmcalJet(),
    fOutlierData(nullptr),
    fHistNonOutliers(nullptr),
    fHistOutliers(nullptr)
{

}

AliAnalysisTaskEmcalResponseOutliers::AliAnalysisTaskEmcalResponseOutliers(const char *name):
    AliAnalysisTaskEmcalJet(name, true),
    fOutlierData(nullptr),
    fHistNonOutliers(nullptr),
    fHistOutliers(nullptr)
{
    SetIsPythia(true);
    SetUsePtHardBinScaling(true);
    SetMakeGeneralHistograms(true);
    ResetMCFilter();
}

AliAnalysisTaskEmcalResponseOutliers::~AliAnalysisTaskEmcalResponseOutliers() {
    if(fOutlierData) delete fOutlierData;
}

void AliAnalysisTaskEmcalResponseOutliers::UserCreateOutputObjects(){
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

    std::string varlist = "ptpart";
    varlist += ":ptdet";
    varlist += ":deltaR";
    varlist += ":etapart";
    varlist += ":phipart";
    varlist += ":etadet";
    varlist += ":phidet";
    varlist += ":nefpart";
    varlist += ":nefdet";
    varlist += ":npart";
    varlist += ":ncharged";
    varlist += ":nneutral";
    varlist += ":leadingpart";
    varlist += ":leadingcharged";
    varlist += ":leadingneutral";
    varlist += ":pdgmaxpart";
    fOutlierData = new TNtuple("fOutlierData", "Outlier information", varlist.data()); 

    fHistNonOutliers = new TH2D("fHistNonOutliers", "histogram for non-outlier jets; p_{t,det} (GeV/c); p_{t,part} (GeV/c)", 500, 0., 500., 500, 0., 500.);
    fHistOutliers = new TH2D("fHistOutliers", "histogram for outlier jets; p_{t,det} (GeV/c); p_{t,part} (GeV/c)", 500, 0., 500., 500, 0., 500.);

    fOutput->Add(fOutlierData);
    fOutput->Add(fHistNonOutliers);
    fOutput->Add(fHistOutliers);

    PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalResponseOutliers::Run(){
    auto detjets = GetJetContainer("detjets");
    auto partcont = GetParticleContainer("mcparticles");

    // Outlier selection based on pt-hard bin
    const double detptmax[21] = {0., 10., 20., 20., 30., 40., 50., 60., 80., 100., 120., 140., 150., 300., 300., 300., 300., 300., 300., 300., 300.},
                 partptmin[21] = {0., 20., 25., 30., 40., 50., 70., 80., 100., 120., 140., 180., 200., 250., 270., 300., 350., 380., 420., 450., 600.};
    AliDebugStream(1) << "Using outlier cuts for pt-hard bin " << fPtHardBinGlobal << ": pt max det: " << detptmax[fPtHardBinGlobal] 
                      << ", pt min part: " << partptmin[fPtHardBinGlobal] << std::endl;

    AliDebugStream(1) << "Found " << detjets->GetNEntries() << " jets at detector level" << std::endl;

    for(auto detjet : detjets->accepted()) {
        auto mcjet = detjet->ClosestJet();
        if(!mcjet) continue;
        if(!(detjet->Pt() < detptmax[fPtHardBinGlobal] && mcjet->Pt() > partptmin[fPtHardBinGlobal])){
            fHistNonOutliers->Fill(detjet->Pt(), mcjet->Pt());
            continue;
        } 
        // jet identified as outlier
        fHistOutliers->Fill(detjet->Pt(), mcjet->Pt());
        auto truepart = mcjet->GetLeadingTrack(partcont->GetArray());
        TVector3 detvec, mcvec;
        detvec.SetPtEtaPhi(detjet->Pt(), detjet->Eta(), detjet->Phi());
        mcvec.SetPtEtaPhi(mcjet->Pt(), mcjet->Eta(), mcjet->Phi());
        Float_t outlierDataBlock[16];
        outlierDataBlock[0] = mcjet->Pt();
        outlierDataBlock[1] = detjet->Pt();
        outlierDataBlock[2] = detvec.DeltaR(mcvec);
        outlierDataBlock[3] = mcjet->Eta();
        outlierDataBlock[4] = mcjet->Phi();
        outlierDataBlock[5] = detjet->Eta();
        outlierDataBlock[6] = detjet->Phi();
        outlierDataBlock[7] = mcjet->NEF();
        outlierDataBlock[8] = detjet->NEF();
        outlierDataBlock[9] = mcjet->GetNumberOfTracks() + mcjet->GetNumberOfClusters();
        outlierDataBlock[10] = detjet->GetNumberOfTracks();
        outlierDataBlock[11] = detjet->GetNumberOfClusters();
        outlierDataBlock[12] = mcjet->MaxPartPt();
        outlierDataBlock[13] = detjet->MaxChargedPt();
        outlierDataBlock[14] = detjet->MaxNeutralPt();
        outlierDataBlock[15] = static_cast<Float_t>(truepart->PdgCode());
        fOutlierData->Fill(outlierDataBlock);
    }
    return true;
}

AliAnalysisTaskEmcalResponseOutliers *AliAnalysisTaskEmcalResponseOutliers::AddTaskEmcalResponseOutliers(){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale", "No analysis manager available");
    return nullptr;
  } 

  auto inputhandler = mgr->GetInputEventHandler();
  auto isAOD = inputhandler->IsA() == AliAODInputHandler::Class();

  auto responsetask = new AliAnalysisTaskEmcalResponseOutliers("ResponseOutlierTask");
  mgr->AddTask(responsetask);


  auto partcont = responsetask->AddMCParticleContainer("mcparticles");
  partcont->SetMinPt(0.);
  auto clusters = responsetask->AddClusterContainer(AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
  clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  clusters->SetClusHadCorrEnergyCut(0.3);
  auto tracks = responsetask->AddTrackContainer(AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
  auto contdetjet = responsetask->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::E_scheme, 0.2,
                                                     AliJetContainer::kEMCALfid, tracks, clusters);
  contdetjet->SetName("detjets");

  std::stringstream outnamebuilder;
  outnamebuilder << mgr->GetCommonFileName() << ":OutlierResponse";

  mgr->ConnectInput(responsetask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(responsetask, 1, mgr->CreateContainer("OutlierResponseHists", TList::Class(), AliAnalysisManager::kOutputContainer, outnamebuilder.str().data()));

  return responsetask;
}