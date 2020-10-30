/************************************************************************************
 * Copyright (C) 2020, Copyright Holders of the ALICE Collaboration                 *
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
#include <THistManager.h>
#include <TLinearBinning.h>
#include <TCustomBinning.h>
#include <TString.h>
#include <TVector2.h>

#include "AliAnalysisTaskEmcalJetSpectrumSDPart.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetSpectrumSDPart)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalJetSpectrumSDPart::AliAnalysisTaskEmcalJetSpectrumSDPart():
    AliAnalysisTaskEmcalJet(),
    AliAnalysisEmcalSoftdropHelperImpl(),
    fHistos(nullptr),
    fBeta(0),
    fZcut(0.1),
    fUseChargedConstituents(true),
    fUseNeutralConstituents(true)
{
}

AliAnalysisTaskEmcalJetSpectrumSDPart::AliAnalysisTaskEmcalJetSpectrumSDPart(const char *name):
    AliAnalysisTaskEmcalJet(name, true),
    AliAnalysisEmcalSoftdropHelperImpl(),
    fHistos(nullptr),
    fBeta(0),
    fZcut(0.1),
    fUseChargedConstituents(true),
    fUseNeutralConstituents(true)
{
    SetMakeGeneralHistograms(true);
}

AliAnalysisTaskEmcalJetSpectrumSDPart::~AliAnalysisTaskEmcalJetSpectrumSDPart()
{
}


void AliAnalysisTaskEmcalJetSpectrumSDPart::UserCreateOutputObjects()
{
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
    fHistos = new THistManager(Form("Histos%s", GetName()));

    // Event property
    fHistos->CreateTH1("fNevents", "Number of events", 1, 0.5, 1.5);

    // Jet spectrum and QA
    fHistos->CreateTH1("hJetPt", "Jet Pt spectrum", 1000, 0., 1000.);
    fHistos->CreateTH2("hJetEtaPhi", "Jet #eta and #phi", 100, -1., 1, 100, 0., TMath::TwoPi());
    fHistos->CreateTH2("hJetNEFPt", "Neutral energy fraction vs. pt", 500, 0., 500., 100, 0., 1.);
    fHistos->CreateTH2("hJetNconstPt", "Number of jet constituents vs. pt", 500, 0., 500., 100, 0., 100.);
    fHistos->CreateTH2("hPtLeading", "Pt of the leading constituent vs. jet pt", 500, 0., 500., 500, 0., 500.);
    fHistos->CreateTH2("hDrLeading", "DeltaR of the leading constituent vs. jet pt", 500, 0., 500., 100, 0., 1.);
    fHistos->CreateTH2("hFailedSD", "Jets failing SoftDrop", 300, 0., 300., 101, -0.5, 100.5);

    // SoftDrop
    double R = double(int(GetJetContainer("partjets")->GetJetRadius() * 1000.))/1000.;  // Save cast from float to double truncating after 3rd decimal digit
    TLinearBinning ptbinning(500, 0., 500.),
                   nsdbinning(22, -1.5, 20.5),
                   thetagbinning(11, -0.5, 1.);
    std::unique_ptr<TBinning> zgbinning(GetZgBinning(fZcut)),
                              rgbinning(GetRgBinning(R));
    fHistos->CreateTH2("hSDZg", "Zg vs. pt", *zgbinning, ptbinning);
    fHistos->CreateTH2("hSDRg", "Rg vs. pt", *rgbinning, ptbinning);
    fHistos->CreateTH2("fSDNsd", "Nsd vs. pt", nsdbinning, ptbinning);
    fHistos->CreateTH2("fSDThetag", "Thetag vs. pt", thetagbinning, ptbinning);

    for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);
    PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalJetSpectrumSDPart::Run()
{
    fHistos->FillTH1("fNevents", 1.);
    auto jets = GetJetContainer("partjets");

    if(!jets) {
        AliErrorStream() << GetName() << ": Part. level jet container not found" << std::endl;
        return false;
    }

    Double_t vertex[3] = {0,0,0};
    AliAnalysisEmcalSoftdropHelperImpl::SoftdropParams sdsettings;
    sdsettings.fBeta = fBeta;
    sdsettings.fZcut = fZcut;
    sdsettings.fReclusterizer = AliAnalysisEmcalSoftdropHelperImpl::EReclusterizer_t::kCAAlgo;
    sdsettings.fUseChargedConstituents = fUseChargedConstituents;
    sdsettings.fUseNeutralConstituents = fUseNeutralConstituents;

    for(auto j : jets->accepted()) {
        fHistos->FillTH1("hJetPt", j->Pt());
        fHistos->FillTH2("hJetEtaPhi", j->Eta(), TVector2::Phi_0_2pi(j->Phi()));
        fHistos->FillTH2("hJetNEFPt", j->Pt(), j->NEF());
        fHistos->FillTH2("hJetNconstPt", j->Pt(), j->N());

        auto leading = j->GetLeadingTrack(jets->GetParticleContainer()->GetArray());
        if(leading) {
            fHistos->FillTH2("hPtLeading", j->Pt(), leading->Pt());
            TVector3 jetvec(j->Px(), j->Py(), j->Px()),
                     leadingvec(leading->Px(), leading->Py(), leading->Pz());
            fHistos->FillTH2("hDrLeading", j->Pt(), jetvec.DeltaR(leadingvec));
        }

        // SoftDrop
        try {
            auto sdparams = this->MakeSoftdrop(*j, jets->GetJetRadius(), true, sdsettings, AliVCluster::VCluUserDefEnergy_t::kNonLinCorr, vertex);
            auto splittings = this->IterativeDecluster(*j, jets->GetJetRadius(), true, sdsettings, AliVCluster::VCluUserDefEnergy_t::kNonLinCorr, vertex);

            fHistos->FillTH2("hSDZg", sdparams.fZg, j->Pt());
            fHistos->FillTH2("hSDRg", sdparams.fRg, j->Pt());
            fHistos->FillTH2("fSDNsd", splittings.size(), j->Pt());
            fHistos->FillTH2("fSDThetag", sdparams.fRg/jets->GetJetRadius(), j->Pt());
        } catch(...) {
            fHistos->FillTH2("hFailedSD", j->Pt(), j->N());
        }
    }

    return true;
}

AliAnalysisTaskEmcalJetSpectrumSDPart *AliAnalysisTaskEmcalJetSpectrumSDPart::AddTaskEmcalJetSpectrumSDPart(AliJetContainer::EJetType_t jettype, double R, const char *nameparticles) {
    auto mgr = AliAnalysisManager::GetAnalysisManager();
    if(!mgr) {
        AliErrorGeneralStream("AliAnalysisTaskEmcalJetSpectrumSDPart::AddTaskEmcalJetSpectrumSDPart") << "No analysis manager available" << std::endl;
        return nullptr;
    }

    TString rstring = Form("R%02d", int(R * 10.)),
            jtstring = "";

    switch (jettype) {
        case AliJetContainer::kChargedJet: jtstring = "ChargedJet"; break;
        case AliJetContainer::kFullJet:    jtstring = "FullJet"; break;
        case AliJetContainer::kNeutralJet: jtstring = "NeutralJet"; break;
        case AliJetContainer::kUndefinedJetType: break;
    };

    auto task = new AliAnalysisTaskEmcalJetSpectrumSDPart(Form("PartLevelJetTask%s%s", jtstring.Data(), rstring.Data()));
    mgr->AddTask(task);

    // Adding particle and jet container
    auto partcont = task->AddMCParticleContainer(nameparticles);
    partcont->SetMinPt(0.);
    auto jetcont = task->AddJetContainer(jettype, AliJetContainer::antikt_algorithm, AliJetContainer::E_scheme, R, AliJetContainer::kTPCfid, partcont, nullptr);
    jetcont->SetName("partjets");
    jetcont->SetMaxTrackPt(1000);
    jetcont->SetMinPt(0);
    jetcont->SetMaxPt(1000.);
    
    // Link input and output
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("PartLevelJetResults%s%s", jtstring.Data(), rstring.Data()), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()));

    return task;
}