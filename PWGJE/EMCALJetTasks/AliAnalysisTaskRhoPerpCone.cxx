/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
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
#include "AliAnalysisTaskRhoPerpCone.h"

#include <TClonesArray.h>
#include <TMath.h>

#include <AliLog.h>
#include <AliVEventHandler.h>
#include <AliAnalysisManager.h>

#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliAODTrack.h"

ClassImp(AliAnalysisTaskRhoPerpCone);

AliAnalysisTaskRhoPerpCone::AliAnalysisTaskRhoPerpCone() : AliAnalysisTaskRhoBase(),
                                                           fHistRhoVsCent(0),
                                                           fHistRhoVsLeadJetPt(),
                                                           fHistRhoVsLeadTrackPt(0),
                                                           fHistRhoVsNtrack(0)
{
}

AliAnalysisTaskRhoPerpCone::AliAnalysisTaskRhoPerpCone(const char *name, Bool_t histo) : AliAnalysisTaskRhoBase(name, histo),
                                                                                         fHistRhoVsCent(0),
                                                                                         fHistRhoVsLeadJetPt(),
                                                                                         fHistRhoVsLeadTrackPt(0),
                                                                                         fHistRhoVsNtrack(0)
{
}

void AliAnalysisTaskRhoPerpCone::UserCreateOutputObjects()
{
    if (!fCreateHisto)
        return;

    AliAnalysisTaskRhoBase::UserCreateOutputObjects();

    TString name;

    Int_t maxTracks = 6000;
    Double_t maxRho = 500;
    Int_t nRhoBins = 500;

    Double_t ptBinWidth = 0.5;
    Double_t maxPt = 250;

    if (fForceBeamType == kpp)
    {
        maxRho = 50;
        maxTracks = 200;
    }
    else if (fForceBeamType == kpA)
    {
        maxRho = 200;
        maxTracks = 500;
    }

    Int_t nPtBins = TMath::CeilNint(maxPt / ptBinWidth);

    fHistRhoVsCent = new TH2F("fHistRhoVsCent", "fHistRhoVsCent", 100, 0, 100, nRhoBins, 0, maxRho);
    fHistRhoVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistRhoVsCent->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoVsCent);

    if (fParticleCollArray.GetEntries() > 0)
    {
        fHistRhoVsNtrack = new TH2F("fHistRhoVsNtrack", "fHistRhoVsNtrack", 200, 0, maxTracks, nRhoBins, 0, maxRho);
        fHistRhoVsNtrack->GetXaxis()->SetTitle("No. of tracks");
        fHistRhoVsNtrack->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
        fOutput->Add(fHistRhoVsNtrack);

        fHistRhoVsLeadTrackPt = new TH2F("fHistRhoVsLeadTrackPt", "fHistRhoVsLeadTrackPt", nPtBins, 0, maxPt, nRhoBins, 0, maxRho);
        fHistRhoVsLeadTrackPt->GetXaxis()->SetTitle("#it{p}_{T,track} (GeV/c)");
        fHistRhoVsLeadTrackPt->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
        fOutput->Add(fHistRhoVsLeadTrackPt);
    }

    name = "fHistRhoVsLeadJetPt";
    fHistRhoVsLeadJetPt = new TH2F(name, name, nPtBins, 0, maxPt, nRhoBins, 0, maxRho);
    fHistRhoVsLeadJetPt->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
    fHistRhoVsLeadJetPt->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistRhoVsLeadJetPt);
}

Bool_t AliAnalysisTaskRhoPerpCone::Run()
{
    Double_t perpPtDensity1 = 0;
    Double_t perpPtDensity2 = 0;

    AliJetContainer *jetCont = static_cast<AliJetContainer *>(fJetCollArray[0]);
    AliParticleContainer *partCont = static_cast<AliParticleContainer *>(fParticleCollArray[0]);

    AliEmcalJet *LeadJet = jetCont->GetLeadingJet();
    if (!LeadJet)
        return kFALSE;

    Double_t dPhi1 = 999.;
    Double_t dPhi2 = 999.;
    Double_t dEta = 999.;
    Double_t PerpendicularConeAxisPhi1 = 999, PerpendicularConeAxisPhi2 = 999;
    PerpendicularConeAxisPhi1 = ((LeadJet->Phi() + (TMath::Pi() / 2.)) > TMath::TwoPi()) ? LeadJet->Phi() - ((3. / 2.) * TMath::Pi()) : LeadJet->Phi() + (TMath::Pi() / 2.);
    PerpendicularConeAxisPhi2 = ((LeadJet->Phi() - (TMath::Pi() / 2.)) < 0) ? LeadJet->Phi() + ((3. / 2.) * TMath::Pi()) : LeadJet->Phi() - (TMath::Pi() / 2.);
    partCont->ResetCurrentID();

    while (auto part = partCont->GetNextAcceptParticle())
    {
        if (!part)
            continue;

        dPhi1 = TMath::Abs(part->Phi() - PerpendicularConeAxisPhi1);
        dPhi1 = (dPhi1 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi1 : dPhi1;
        dPhi2 = TMath::Abs(part->Phi() - PerpendicularConeAxisPhi2);
        dPhi2 = (dPhi2 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi2 : dPhi2;
        dEta = LeadJet->Eta() - part->Eta();
        if (TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) <= (Double_t)jetCont->GetJetRadius())
            perpPtDensity1 += part->Pt();

        if (TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) <= (Double_t)jetCont->GetJetRadius())
            perpPtDensity2 += part->Pt();
    }

    Double_t perpPtDensity = (perpPtDensity1 + perpPtDensity2) / (2 * TMath::Pi() * TMath::Power(jetCont->GetJetRadius(), 2));

    fOutRho->SetVal(perpPtDensity);

    if (fCreateHisto)
    {
        fHistRhoVsCent->Fill(fCent, fOutRho->GetVal());

        if (partCont->GetLeadingParticle())
        {
            fHistRhoVsLeadTrackPt->Fill(partCont->GetLeadingParticle()->Pt(), fOutRho->GetVal());
        }

        fHistRhoVsNtrack->Fill(partCont->GetNAcceptedParticles(), fOutRho->GetVal());

        if (LeadJet)
            fHistRhoVsLeadJetPt->Fill(LeadJet->Pt(), fOutRho->GetVal());
    }

    return kTRUE;
}

AliAnalysisTaskRhoPerpCone *AliAnalysisTaskRhoPerpCone::AddTaskRhoPerpCone(TString trackName, Double_t trackPtCut, TString clusName, Double_t clusECut, TString nRho, Double_t jetradius, UInt_t acceptance, AliJetContainer::EJetType_t jetType, AliJetContainer::ERecoScheme_t rscheme, Bool_t histo, TString suffix)
{
    // Get the pointer to the existing analysis manager via the static access method.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AliAnalysisTaskRhoPerpCone::AddTaskRhoPerpCone", "No analysis manager to connect to.");
        return nullptr;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    AliVEventHandler *handler = mgr->GetInputEventHandler();
    if (!handler)
    {
        ::Error("AliAnalysisTaskRhoPerpCone::AddTaskRhoPerpCone", "This task requires an input event handler");
        return nullptr;
    }

    EDataType_t dataType = kUnknownDataType;

    if (handler->InheritsFrom("AliESDInputHandler"))
    {
        dataType = kESD;
    }
    else if (handler->InheritsFrom("AliAODInputHandler"))
    {
        dataType = kAOD;
    }

    // Init the task and do settings
    if (trackName == "usedefault")
    {
        if (dataType == kESD)
        {
            trackName = "Tracks";
        }
        else if (dataType == kAOD)
        {
            trackName = "tracks";
        }
        else
        {
            trackName = "";
        }
    }

    if (clusName == "usedefault")
    {
        if (dataType == kESD)
        {
            clusName = "CaloClusters";
        }
        else if (dataType == kAOD)
        {
            clusName = "caloClusters";
        }
        else
        {
            clusName = "";
        }
    }

    TString name(TString::Format("AliAnalysisTaskRhoPerpCone_%s", nRho.Data()));
    if (!suffix.IsNull())
    {
        name += "_";
        name += suffix;
    }

    AliAnalysisTaskRhoPerpCone *mgrTask = dynamic_cast<AliAnalysisTaskRhoPerpCone *>(mgr->GetTask(name.Data()));
    if (mgrTask)
    {
        ::Warning("AliAnalysisTaskRhoPerpCone::AddTaskRhoPerpCone", "Not adding the task again, since a task with the same name '%s' already exists", name.Data());
        return mgrTask;
    }

    AliAnalysisTaskRhoPerpCone *rhotask = new AliAnalysisTaskRhoPerpCone(name, histo);
    rhotask->SetOutRhoName(nRho);

    AliTrackContainer *trackCont;
    AliParticleContainer *partCont;
    if (trackName == "mcparticles")
    {
        partCont = rhotask->AddParticleContainer(trackName);
    }
    else if (trackName == "tracks" || trackName == "Tracks")
    {
        trackCont = rhotask->AddTrackContainer(trackName);
        partCont = rhotask->GetParticleContainer(0);
    }
    else
    {
        partCont = rhotask->AddParticleContainer(trackName);
    }
    partCont->SetMinPt(trackPtCut);

    AliClusterContainer *clusterCont = rhotask->AddClusterContainer(clusName.Data());
    if (clusterCont)
    {
        clusterCont->SetClusECut(0.);
        clusterCont->SetClusPtCut(0.);
        clusterCont->SetClusHadCorrEnergyCut(clusECut);
        clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    }

    AliJetContainer *jetCont = new AliJetContainer(jetType, AliJetContainer::antikt_algorithm, rscheme, jetradius, partCont, clusterCont);
    if (jetCont)
    {
        jetCont->SetJetPtCut(1);
        jetCont->SetJetAcceptanceType(acceptance);
        jetCont->SetName("Signal");
        rhotask->AdoptJetContainer(jetCont);
    }

    // Final settings, pass to manager and set the containers
    mgr->AddTask(rhotask);

    // Create containers for input/output
    mgr->ConnectInput(rhotask, 0, mgr->GetCommonInputContainer());
    if (histo)
    {
        TString contname(name);
        contname += "_histos";
        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                                                                  TList::Class(), AliAnalysisManager::kOutputContainer,
                                                                  Form("%s", AliAnalysisManager::GetCommonFileName()));
        mgr->ConnectOutput(rhotask, 1, coutput1);
    }

    return rhotask;
}
