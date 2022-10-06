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

#include <TChain.h>
#include <TH1D.h>
#include <TList.h>
#include <TMath.h>
#include <THnSparse.h>

#include <AliAnalysisManager.h>
#include "AliAODTrack.h"
#include <AliAODInputHandler.h>
#include "AliAnalysisTaskCorrelGen.h"
#include <AliEventPoolManager.h>
#include <AliESDtrackCuts.h>
#include "AliEventCuts.h"

class AliAnalysisTaskCorrelGen;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCorrelGen)
    /// \endcond

    /// Default constructor for ROOT I/O purposes
    AliAnalysisTaskCorrelGen::AliAnalysisTaskCorrelGen() : AliAnalysisTaskSE(),
                                                                 fAliEventCuts(nullptr),
                                                                 fmcEvent(nullptr),
                                                                 fOutputList(nullptr),
                                                                 fHistMCCorr(nullptr),
                                                                 fHistMCMixGen(nullptr),
                                                                 fHistNtrigGen(nullptr),
                                                                 fHistGenMult(nullptr),
                                                                 fMixTrk(5000),
                                                                 fPoolMgr(nullptr),
                                                                 fPoolMCGen(nullptr),
                                                                 fPtTrigMin(3.),
                                                                 fPtAsocMin(1.),
                                                                 fNMixEvts(10),
                                                                 fNbinsPtTrig(17),
                                                                 fNbinsPtAssoc(19),
                                                                 fV0hCorr(kTRUE),
                                                                 fhhCorr(kTRUE),
                                                                 fhV0Corr(kTRUE),
                                                                 fFilterBit(32),
                                                                 fIsAOD(kTRUE),
                                                                 fNbinsDPhi(72),
                                                                 fNbinsDEta(75),
                                                                 fNbinsEta(40),
                                                                 fNbinsPhi(72),
                                                                 fMixGen(kFALSE),
                                                                 fNbinsVz(15),
                                                                 fPrimVtxCut(10),
                                                                 fESDTrkCuts(nullptr),
                                                                 fTestPions(kFALSE),
                                                                 fPtAssocMax(20),
                                                                 fPtTrigMax(20),
                                                                 fMaxDCAtoVtxZ(2),
                                                                 fMaxDCAtoVtxXY(0.3),
                                                                 fMinNCrossedRowsTPCprimtracks(70),
                                                                 fMCTrkSel(0),
                                                                 fMCGenTrkMix(0),
                                                                 fMCTrkTrigSel(0),
                                                                 fMCTrkV0Sel(0),
                                                                 fMCV0AssocSel(0),
                                                                 fMCAssocSel(0),
                                                                 fselectedMCV0assoc(0),
                                                                 fMCTrigSel(0),
                                                                 fMCV0RecTrigSel(0),
                                                                 fTrkSel(0),
                                                                 fTrkAssocSel(0),
                                                                 fTrigTrkSel(0),
                                                                 fV0TrigSel(0),
                                                                 fV0AssocSel(0),
                                                                 fhighMult(kFALSE),
                                                                 fhighMultSPD(kFALSE),
                                                                 fEvtCutsQAPlots(kFALSE),
                                                                 fpp(kTRUE),
                                                                 fNbinsMult(10),
                                                                 fOnFlyMC(kFALSE)
{
    for (Int_t i = 0; i < 3; i++)
    {
        fPV[i] = -1.;
    }
}

/// Standard named constructor
///
/// \param name Name of the task
AliAnalysisTaskCorrelGen::AliAnalysisTaskCorrelGen(const char *name) : AliAnalysisTaskSE(name),
                                                                             fAliEventCuts(nullptr),
                                                                             fmcEvent(nullptr),
                                                                             fOutputList(nullptr),
                                                                             fHistMCCorr(nullptr),
                                                                             fHistMCMixGen(nullptr),
                                                                             fHistGenMult(nullptr),
                                                                             fMixTrk(5000),
                                                                             fPoolMgr(nullptr),
                                                                             fPoolMCGen(nullptr), fHistNtrigGen(nullptr),
                                                                             fPtTrigMin(3),
                                                                             fPtAsocMin(1.),
                                                                             fNMixEvts(10),
                                                                             fNbinsPtTrig(17),
                                                                             fNbinsPtAssoc(19),
                                                                             fV0hCorr(kTRUE),
                                                                             fhhCorr(kTRUE),
                                                                             fhV0Corr(kTRUE),
                                                                             fFilterBit(32),
                                                                             fIsAOD(kTRUE),
                                                                             fNbinsDPhi(72),
                                                                             fNbinsDEta(75),
                                                                             fNbinsEta(40),
                                                                             fNbinsPhi(72),
                                                                             fMixGen(kFALSE),
                                                                             fNbinsVz(15),
                                                                             fPrimVtxCut(10),
                                                                             fESDTrkCuts(nullptr),
                                                                             fTestPions(kFALSE),
                                                                             fPtAssocMax(20),
                                                                             fPtTrigMax(20),
                                                                             fMaxDCAtoVtxZ(2),
                                                                             fMaxDCAtoVtxXY(0.3),
                                                                             fMinNCrossedRowsTPCprimtracks(70),
                                                                             fMCTrkSel(0),
                                                                             fMCGenTrkMix(0),
                                                                             fMCTrkTrigSel(0),
                                                                             fMCTrkV0Sel(0),
                                                                             fMCV0AssocSel(0),
                                                                             fMCAssocSel(0),
                                                                             fselectedMCV0assoc(0),
                                                                             fMCTrigSel(0),
                                                                             fMCV0RecTrigSel(0),
                                                                             fTrkSel(0),
                                                                             fTrkAssocSel(0),
                                                                             fTrigTrkSel(0),
                                                                             fV0TrigSel(0),
                                                                             fV0AssocSel(0),
                                                                             fhighMult(kFALSE),
                                                                             fhighMultSPD(kFALSE),
                                                                             fEvtCutsQAPlots(kFALSE),
                                                                             fpp(kTRUE),
                                                                             fNbinsMult(10),
                                                                             fOnFlyMC(kFALSE)
{
    for (Int_t i = 0; i < 3; i++)
    {
        fPV[i] = -1.;
    }

    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

/// Destructor
AliAnalysisTaskCorrelGen::~AliAnalysisTaskCorrelGen()
{
    if (fESDTrkCuts)
    {
        delete fESDTrkCuts;
        fESDTrkCuts = nullptr;
    }
    if (fOutputList)
    {
        delete fOutputList;
    }
    if (fMCTrkSel)
        delete fMCTrkSel;
    if (fMCGenTrkMix)
        delete fMCGenTrkMix;
    if (fMCTrkTrigSel)
        delete fMCTrkTrigSel;
    if (fMCTrkV0Sel)
        delete fMCTrkV0Sel;
    if (fMCV0AssocSel)
        delete fMCV0AssocSel;
    if (fMCAssocSel)
        delete fMCAssocSel;
    if (fselectedMCV0assoc)
        delete fselectedMCV0assoc;
    if (fMCTrigSel)
        delete fMCTrigSel;
    if (fMCV0RecTrigSel)
        delete fMCV0RecTrigSel;
    if (fTrkSel)
        delete fTrkSel;
    if (fTrkAssocSel)
        delete fTrkAssocSel;
    if (fTrigTrkSel)
        delete fTrigTrkSel;
    if (fV0TrigSel)
        delete fV0TrigSel;
    if (fV0AssocSel)
        delete fV0AssocSel;
}

/// Overloads base class method. Creates output objects
void AliAnalysisTaskCorrelGen::UserCreateOutputObjects()
{
    Int_t NCentBins = 10;
    Double_t multBins[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100};

    const Int_t NbinsVtx = fNbinsVz;

    Double_t ZBins[NbinsVtx + 1];
    ZBins[0] = -1. * fPrimVtxCut;

    if (NbinsVtx == 9)
    {
        ZBins[1] = -7.0;
        for (Int_t i = 2; i < NbinsVtx; i++)
        {
            ZBins[i] = ZBins[i - 1] + 2;
        }
        ZBins[9] = fPrimVtxCut;
    }
    else
    {
        Double_t binstep = Double_t(2 * fPrimVtxCut) / NbinsVtx;
        for (Int_t i = 1; i < NbinsVtx + 1; i++)
        {
            ZBins[i] = ZBins[i - 1] + binstep;
        }
    }

    Int_t bins2d[5] = {fNbinsPtTrig, fNbinsVz, 7, 902, fNbinsMult};
    Double_t mis2d[5] = {fPtTrigMin, -10, 0., 0.44, 0};
    Double_t maxs2d[5] = {fPtTrigMax, 10, 7., 1.355, 100};

    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);

    Double_t *binsMass = new Double_t[903];
    binsMass[0] = 0.44;
    for (Int_t i = 0; i < 901; i++)
    {
        if (i < 301)
            binsMass[i + 1] = binsMass[i] + (0.56 - 0.44) / 300;
        if (i == 301)
            binsMass[i + 1] = 1.08;
        if (i > 301 && i < 602)
            binsMass[i + 1] = binsMass[i] + (1.15 - 1.08) / 300;
        if (i == 602)
            binsMass[i + 1] = 1.285;
        if (i > 602)
            binsMass[i + 1] = binsMass[i] + (1.355 - 1.285) / 300;
    }
    binsMass[902] = 1.355;

    Int_t bins[9] = {fNbinsPtTrig, fNbinsPtAssoc, fNbinsDPhi, fNbinsDEta, fNbinsVz, 12, 902, fNbinsMult, 2};
    Double_t min[9] = {fPtTrigMin, fPtAsocMin, -TMath::PiOver2(), -2., -10., 0., 0.44, 0, -2.};
    Double_t max[9] = {fPtTrigMax, fPtAssocMax, 3. * TMath::PiOver2(), 2., 10., 12., 1.355, 100, 2.};

    bins[7] = 500;
    max[7] = 500;

    bins[6] = 2;
    max[6] = 2;
    min[6] = -2;

    fHistMCCorr = new THnSparseF("fHistMCCorr", "fHistMCCorr", 9, bins, min, max);
    fHistMCCorr->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistMCCorr->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistMCCorr->GetAxis(2)->SetTitle("#Delta#varphi");
    fHistMCCorr->GetAxis(3)->SetTitle("#Delta#eta");
    fHistMCCorr->GetAxis(4)->SetTitle("p_{vz}");
    fHistMCCorr->GetAxis(5)->SetTitle("trigger");
    fHistMCCorr->GetAxis(6)->SetTitle("original parton");
    fHistMCCorr->GetAxis(7)->SetTitle("multiplicity");
    fHistMCCorr->GetAxis(8)->SetTitle("hadron charge");
    fHistMCCorr->Sumw2();
    fOutputList->Add(fHistMCCorr);
    fHistMCCorr->GetAxis(4)->Set(NbinsVtx, ZBins);

    fHistMCMixGen = new THnSparseF("fHistMCMixGen", "fHistMCMixGen", 9, bins, min, max);
    fHistMCMixGen->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistMCMixGen->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistMCMixGen->GetAxis(2)->SetTitle("#Delta#phi");
    fHistMCMixGen->GetAxis(3)->SetTitle("#Delta#eta");
    fHistMCMixGen->GetAxis(4)->SetTitle("p_{vz}");
    fHistMCMixGen->GetAxis(5)->SetTitle("trigger");
    fHistMCMixGen->GetAxis(6)->SetTitle("mass");
    fHistMCMixGen->GetAxis(7)->SetTitle("multiplicity percentile");
    fHistMCMixGen->GetAxis(8)->SetTitle("hadron charge");
    fOutputList->Add(fHistMCMixGen);
    fHistMCMixGen->Sumw2();
    fHistMCMixGen->GetAxis(4)->Set(NbinsVtx, ZBins);
    fHistMCMixGen->GetAxis(6)->Set(102, binsMass);

    bins2d[4] = 500;
    maxs2d[4] = 500;

    bins2d[3] = 2;
    maxs2d[3] = 2;
    mis2d[3] = -2;

    fHistNtrigGen = new THnSparseF("fHistNtrigGen", "fHistNtrigGen", 5, bins2d, mis2d, maxs2d);
    fHistNtrigGen->GetAxis(0)->SetTitle("p_{T}");
    fHistNtrigGen->GetAxis(1)->SetTitle("p_{vz}");
    fHistNtrigGen->GetAxis(2)->SetTitle("trigger");
    fHistNtrigGen->GetAxis(3)->SetTitle("original parton");
    fHistNtrigGen->GetAxis(4)->SetTitle("multiplicity");
    fOutputList->Add(fHistNtrigGen);
    fHistNtrigGen->Sumw2();
    fHistNtrigGen->GetAxis(1)->Set(NbinsVtx, ZBins);

    fHistGenMult = new TH1D("fHistGenMult", "fHistGenMult", 500, 0, 500);
    fOutputList->Add(fHistGenMult);

    if (!fIsAOD)
    {
        if (fFilterBit == 32)
            fESDTrkCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
        if (fFilterBit == 16 || fFilterBit == 256)
            fESDTrkCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
        if (fFilterBit == 128)
            fESDTrkCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
        if (fFilterBit == 1)
        {
            fESDTrkCuts = new AliESDtrackCuts;
            fESDTrkCuts->SetMaxChi2PerClusterTPC(4);
            fESDTrkCuts->SetAcceptKinkDaughters(kFALSE);
            fESDTrkCuts->SetDCAToVertex2D(kFALSE);

            fESDTrkCuts->SetMaxDCAToVertexZ(fMaxDCAtoVtxZ);
            fESDTrkCuts->SetMaxDCAToVertexXY(fMaxDCAtoVtxXY);
            fESDTrkCuts->SetRequireTPCRefit(kTRUE);
            fESDTrkCuts->SetRequireITSRefit(kTRUE);
            fESDTrkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kAny);
            fESDTrkCuts->SetRequireSigmaToVertex(kFALSE);

            fESDTrkCuts->SetMinNCrossedRowsTPC(fMinNCrossedRowsTPCprimtracks);
            fESDTrkCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
        }
    }
    fMCTrkSel = new TObjArray();
    fMCTrkSel->SetOwner(kTRUE);
    fMCGenTrkMix = new TObjArray();
    fMCGenTrkMix->SetOwner(kTRUE);
    fMCTrkTrigSel = new TObjArray();
    fMCTrkTrigSel->SetOwner(kTRUE);
    fMCTrkV0Sel = new TObjArray();
    fMCTrkV0Sel->SetOwner(kTRUE);
    fMCV0AssocSel = new TObjArray();
    fMCV0AssocSel->SetOwner(kTRUE);
    fMCAssocSel = new TObjArray();
    fMCAssocSel->SetOwner(kTRUE);
    fselectedMCV0assoc = new TObjArray();
    fselectedMCV0assoc->SetOwner(kTRUE);
    fMCTrigSel = new TObjArray();
    fMCTrigSel->SetOwner(kTRUE);
    fMCV0RecTrigSel = new TObjArray();
    fMCV0RecTrigSel->SetOwner(kTRUE);
    fTrkSel = new TObjArray();
    fTrkAssocSel = new TObjArray();
    fTrkAssocSel->SetOwner(kTRUE);
    fTrigTrkSel = new TObjArray();
    fTrigTrkSel->SetOwner(kTRUE);
    fV0TrigSel = new TObjArray();
    fV0TrigSel->SetOwner(kTRUE);
    fV0AssocSel = new TObjArray();
    fV0AssocSel->SetOwner(kTRUE);

    fAliEventCuts = new AliEventCuts();
    if (fEvtCutsQAPlots)
        fAliEventCuts->AddQAplotsToList(fOutputList, kTRUE);
    if (fpp)
        fAliEventCuts->SetupRun2pp();
    else
        fAliEventCuts->SetupRun2PbPb();
    if (fhighMult)
        fAliEventCuts->OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, kTRUE);
    if (fhighMultSPD)
        fAliEventCuts->OverrideAutomaticTriggerSelection(AliVEvent::kHighMultSPD, kTRUE);

    PostData(1, fOutputList); // Post data for ALL output slots > 0 here.

    // Settings for event mixing
    Int_t trackDepth = fMixTrk;
    Int_t poolSize = 100; // Maximum number of events, ignored in the present implemented of AliEventPoolManager

    fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, NCentBins, multBins, NbinsVtx, ZBins);

    fPoolMgr->SetTargetValues(trackDepth, 0.1, 5);
}

/// This function is called once for each event,
/// the manager will take care of reading the events from file,
/// and with the static function InputEvent() you have access to the current event.
/// Once you return from the UserExec function, the manager will retrieve the next event from the chain
void AliAnalysisTaskCorrelGen::UserExec(Option_t *)
{
    Double_t lPercentile = 302.;
    Int_t partCount = 0;
    AliMCParticle *mcTrack = nullptr;

    fmcEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fmcEvent)
    {
        AliFatal("No MC particle branch found");
        return;
    }

    AliVVertex *mcVertex = (AliVVertex *)fmcEvent->GetPrimaryVertex();
    fPV[2] = mcVertex->GetZ();
    if (TMath::Abs(fPV[2]) >= fPrimVtxCut)
        return;

    Int_t allTrksMC = fmcEvent->GetNumberOfTracks();

    AliVTrack *genTrkMix = nullptr;

    for (Int_t i = 0; i < allTrksMC; i++)
    {
        AliMCParticle *mcTrack = (AliMCParticle *)fmcEvent->GetTrack(i);
        if (!mcTrack)
        {
            ::Error("ReadEventAODMC", "Could not receive particle %d", i);
            continue;
        }

        Double_t trkEta = mcTrack->Eta();

        if (mcTrack->IsPhysicalPrimary() && mcTrack->Charge() != 0 && ((trkEta > -3.7 && trkEta < -1.7) || (trkEta > 2.8 && trkEta < 5.1)))
        {
            partCount++;
        }
    }

    fHistGenMult->Fill(partCount);

    AliMCParticle *mcMotherParticle = nullptr;
    AliMCParticle *daughter0 = nullptr;
    AliMCParticle *daughter1 = nullptr;

    for (Int_t i = 0; i < allTrksMC; i++)
    {
        mcTrack = (AliMCParticle *)fmcEvent->GetTrack(i);
        if (!mcTrack)
        {
            ::Error("ReadEventAODMC", "Could not receive particle %d", i);
            continue;
        }
        // track cuts for generated particles
        Double_t mcTrkEta = mcTrack->Eta();
        Double_t mcTrkPt = mcTrack->Pt();
        Bool_t IsPrimTrk = mcTrack->IsPhysicalPrimary();
        Bool_t IsMaxTrkEta = TMath::Abs(mcTrkEta) < 0.8;
        Bool_t IsMinTrkPt = mcTrkPt > fPtAsocMin;
        Bool_t IsTrkCh = (mcTrack->Charge()) != 0;
        Short_t ch;
        if (mcTrack->Charge() > 0)
            ch = 1.;
        else if (mcTrack->Charge() < 0)
            ch = -1.;
        else
            ch = 0;
        Int_t origPartonType = -30; // 0 - quark, 1 - gluon
        Int_t PDGparton = 10;

        if (IsPrimTrk && IsMinTrkPt && IsTrkCh && IsMaxTrkEta)
        {
            fMCTrkSel->Add(new AliV0ChParticleGen(mcTrack->Eta(), mcTrack->Phi(), mcTrack->Pt(), 4, mcTrack->GetLabel(), mcTrack->GetLabel(), ch));

            if (fMixGen)
            {
                genTrkMix = SetAliAODTrack(mcTrack->Theta(), mcTrack->Phi(), mcTrack->Pt(), mcTrack->Charge());
                fMCGenTrkMix->Add(genTrkMix);
                lPercentile = 10;
            }

            PDGparton = GetOrigParton(mcTrack);

            if (TMath::Abs(PDGparton) < 7)
                origPartonType = 0;
            else if (TMath::Abs(PDGparton) == 21)
                origPartonType = 1;

            if (mcTrkPt > fPtTrigMin)
            {
                fMCTrkTrigSel->Add(new AliV0ChParticleGen(mcTrack->Eta(), mcTrack->Phi(), mcTrack->Pt(), 4, mcTrack->GetLabel(), mcTrack->GetLabel(), origPartonType));
            }
        }

        Int_t mcPartPdg = mcTrack->PdgCode();
        Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
        Double_t V0genrapidity = mcTrack->Y();

        Bool_t IsK0, IsLambda, IsAntiLambda, IsPositiveXi, IsNegativeXi;
        Double_t etaDau0 = 0.;
        Double_t etaDau1 = 0.;

        Int_t labelPos = -1;
        Int_t labelNeg = -1;
        origPartonType = -30; // 0 - quark, 1 - gluon

        if (fTestPions)
        {
            if ((mcPartPdg != 111) && (mcPartPdg != 211) && (mcPartPdg != (-211)))
                continue; // keep only pions

            IsK0 = mcPartPdg == 111 && (isPhysPrim);
            IsLambda = mcPartPdg == 211 && (isPhysPrim);
            IsAntiLambda = mcPartPdg == -211 && (isPhysPrim);
        }
        else
        {
            if ((mcPartPdg != 310) && (mcPartPdg != 3122) && (mcPartPdg != (-3122)) && TMath::Abs(mcPartPdg) != 3312)
                continue; // keep only Lambdas and K0S and charged Xi (for the feedDownCorrection)
            Bool_t IsFromCascade = kFALSE;

            if (!fOnFlyMC)
            {
                Int_t mother = mcTrack->GetMother();
                mcMotherParticle = static_cast<AliMCParticle *>(fmcEvent->GetTrack(mother));
                Int_t motherPDG = 0;
                if (mother < 0)
                    motherPDG = 0;
                else
                    motherPDG = TMath::Abs(mcMotherParticle->PdgCode());

                Int_t dau0 = mcTrack->GetDaughterLabel(0);
                if (dau0 > 0)
                    daughter0 = (AliMCParticle *)fmcEvent->GetTrack(dau0);
                Int_t dau1 = mcTrack->GetDaughterLabel(1);
                if (dau1 > 0)
                    daughter1 = (AliMCParticle *)fmcEvent->GetTrack(dau1);

                if (!daughter0 || !daughter1)
                    continue;

                if (daughter0->Charge() < 0)
                {
                    labelPos = daughter1->GetLabel();
                    labelNeg = daughter0->GetLabel();
                }
                if (daughter0->Charge() > 0)
                {
                    labelPos = daughter0->GetLabel();
                    labelNeg = daughter1->GetLabel();
                }
                else
                {
                    labelPos = daughter0->GetLabel();
                    labelNeg = daughter1->GetLabel();
                }

                etaDau0 = daughter0->Eta();
                etaDau1 = daughter1->Eta();
            }
            PDGparton = GetOrigParton(mcTrack);

            if (TMath::Abs(PDGparton) < 7)
                origPartonType = 0;
            if (TMath::Abs(PDGparton) == 21)
                origPartonType = 1;

            IsK0 = mcPartPdg == 310 && (isPhysPrim);
            IsLambda = mcPartPdg == 3122 && (isPhysPrim || IsFromCascade);
            IsAntiLambda = mcPartPdg == -3122 && (isPhysPrim || IsFromCascade);
        }
        IsPositiveXi = mcPartPdg == -3312 && (isPhysPrim);
        IsNegativeXi = mcPartPdg == 3312 && (isPhysPrim);

        if (mcTrack->Pt() > fPtAsocMin && TMath::Abs(V0genrapidity) < 0.5 && TMath::Abs(etaDau0) < 0.8 && TMath::Abs(etaDau1) < 0.8)
        {
            if (IsK0)
            {
                fMCV0AssocSel->Add(new AliV0ChParticleGen(mcTrack->Eta(), mcTrack->Phi(), mcTrack->Pt(), 5, mcTrack->GetLabel(), labelPos, labelNeg, mcTrack->M()));
            }
            if (IsLambda)
            {
                fMCV0AssocSel->Add(new AliV0ChParticleGen(mcTrack->Eta(), mcTrack->Phi(), mcTrack->Pt(), 6, mcTrack->GetLabel(), labelPos, labelNeg, mcTrack->M()));
            }
            if (IsAntiLambda)
            {
                fMCV0AssocSel->Add(new AliV0ChParticleGen(mcTrack->Eta(), mcTrack->Phi(), mcTrack->Pt(), 7, mcTrack->GetLabel(), labelPos, labelNeg, mcTrack->M()));
            }
        }
        if (mcTrack->Pt() > fPtTrigMin && TMath::Abs(V0genrapidity) < 0.5 && TMath::Abs(etaDau0) < 0.8 && TMath::Abs(etaDau1) < 0.8)
        {
            if (IsK0)
                fMCTrkV0Sel->Add(new AliV0ChParticleGen(mcTrack->Eta(), mcTrack->Phi(), mcTrack->Pt(), 1, mcTrack->GetLabel(), labelPos, labelNeg, mcTrack->M(), origPartonType));
            if (IsLambda)
                fMCTrkV0Sel->Add(new AliV0ChParticleGen(mcTrack->Eta(), mcTrack->Phi(), mcTrack->Pt(), 2, mcTrack->GetLabel(), labelPos, labelNeg, mcTrack->M(), origPartonType));
            if (IsAntiLambda)
                fMCTrkV0Sel->Add(new AliV0ChParticleGen(mcTrack->Eta(), mcTrack->Phi(), mcTrack->Pt(), 3, mcTrack->GetLabel(), labelPos, labelNeg, mcTrack->M(), origPartonType));
        }

        if (mcTrack->Pt() > fPtAsocMin && TMath::Abs(V0genrapidity) < 0.5)
        {
            if (IsPositiveXi)
            {
                fMCV0AssocSel->Add(new AliV0ChParticleGen(mcTrack->Eta(), mcTrack->Phi(), mcTrack->Pt(), 11, mcTrack->M(), 1, 2, 3));
            }
            if (IsNegativeXi)
            {
                fMCV0AssocSel->Add(new AliV0ChParticleGen(mcTrack->Eta(), mcTrack->Phi(), mcTrack->Pt(), 10, mcTrack->M(), 1, 2, 3));
            }
        }
        if (mcTrack->Pt() > fPtTrigMin && TMath::Abs(V0genrapidity) < 0.5)
        {
            if (IsNegativeXi)
                fMCTrkV0Sel->Add(new AliV0ChParticleGen(mcTrack->Eta(), mcTrack->Phi(), mcTrack->Pt(), 8, mcTrack->M(), 1, 2, 3));
            if (IsPositiveXi)
                fMCTrkV0Sel->Add(new AliV0ChParticleGen(mcTrack->Eta(), mcTrack->Phi(), mcTrack->Pt(), 9, mcTrack->M(), 1, 2, 3));
        }
    }

    // V0-h
    if (fV0hCorr)
        Correlation(fMCTrkV0Sel, fMCTrkSel, kFALSE, kTRUE, partCount, kFALSE);

    // h-h
    if (fhhCorr)
        Correlation(fMCTrkTrigSel, fMCTrkSel, kFALSE, kFALSE, partCount, kFALSE);
    // h-V0
    if (fhV0Corr)
        Correlation(fMCTrkTrigSel, fMCV0AssocSel, kFALSE, kTRUE, partCount, kTRUE);

    // Doing mixing
    if (fMixGen)
    {
        fPoolMCGen = fPoolMgr->GetEventPool(lPercentile, fPV[2]);
        if (!fPoolMCGen)
        {
            AliWarning(Form("No pool MC Gen found for centrality = %f, zVtx = %f", lPercentile, fPV[2]));
            return;
        }
        Int_t nMixGen = fPoolMCGen->GetCurrentNEvents();
        if (fPoolMCGen->IsReady() || fPoolMCGen->NTracksInPool() > fMixTrk / 10 || nMixGen >= fNMixEvts)
        {
            for (Int_t jMix = 0; jMix < nMixGen; jMix++)
            {
                TObjArray *bgTracksGen = fPoolMCGen->GetEvent(jMix);
                if (fV0hCorr)
                    CorrelationMixing(fMCTrkV0Sel, bgTracksGen, partCount);
                if (fhhCorr)
                    CorrelationMixing(fMCTrkTrigSel, bgTracksGen, partCount);
                if (fhV0Corr)
                    CorrelationMixinghV0(bgTracksGen, fMCV0AssocSel, partCount);
            }
        }

        TObjArray *cloneArrayMCGen = (TObjArray *)fMCGenTrkMix->Clone();
        cloneArrayMCGen->SetOwner(kTRUE);
        fPoolMCGen->UpdatePool(cloneArrayMCGen);
    }

    fMCTrkSel->Clear();
    fMCTrkTrigSel->Clear();
    fMCTrkV0Sel->Clear();
    fMCAssocSel->Clear();
    fMCTrigSel->Clear();
    fMCV0RecTrigSel->Clear();
    fTrkAssocSel->Clear();
    fTrigTrkSel->Clear();
    fV0TrigSel->Clear();
    fV0AssocSel->Clear();
    fMCV0AssocSel->Clear();
    fselectedMCV0assoc->Clear();
    fMCGenTrkMix->Clear();
    fTrkSel->Clear();

    PostData(1, fOutputList); // stream the results the analysis of this event to
                              // the output manager which will take care of writing
                              // it to a file
}

/// This function is called once at the end of the analysis, after processing all events.
void AliAnalysisTaskCorrelGen::Terminate(Option_t *)
{
}

///
void AliAnalysisTaskCorrelGen::Correlation(TObjArray *triggers, TObjArray *associated, Bool_t hh, Bool_t V0h, Float_t perc, Bool_t hV0)
{
    Int_t nAssoc = associated->GetEntriesFast();
    Int_t nTrig = triggers->GetEntriesFast();
    Double_t trigPt = 0.;
    Double_t trigEta = 0.;
    Double_t trigPhi = 0;
    Double_t trigMass = 0.;
    Double_t assocPt = 0.;
    Double_t assocEta = 0.;
    Double_t assocPhi = 0.;
    Double_t assocCharge = 0.;

    for (Int_t i = 0; i < nTrig; i++)
    {
        AliV0ChParticleGen *trig = (AliV0ChParticleGen *)triggers->At(i);
        trigPt = trig->Pt();
        trigEta = trig->Eta();
        trigPhi = trig->Phi();
        if (trigPt < fPtTrigMin)
            continue;
        if (trig->WhichCandidate() > 7)
            continue;
        if (trig->WhichCandidate() < 4)
            trigMass = trig->M();

        if (!hV0)
        {
            if (trig->GetOrigalPartonType() == 0)
                trigMass = -1; // quark is the orinal parton
            else if (trig->GetOrigalPartonType() == 1)
                trigMass = 1; // gluon is the orinal parton
            else
                trigMass = -30;

            Double_t triggers[5] = {trigPt, fPV[2], trig->WhichCandidate() - 0.5, trigMass, perc};
            fHistNtrigGen->Fill(triggers);
        }
        else
        {
            Double_t triggers[5] = {trigPt, fPV[2], 4.5, trigMass, perc};
            fHistNtrigGen->Fill(triggers);
        }

        for (Int_t j = 0; j < nAssoc; j++)
        {
            AliV0ChParticleGen *assoc = (AliV0ChParticleGen *)associated->At(j);
            if (assoc->WhichCandidate() < 4 || assoc->WhichCandidate() > 7)
                continue;
            assocPt = assoc->Pt();
            assocEta = assoc->Eta();
            assocPhi = assoc->Phi();
            assocCharge = assoc->Charge();

            Double_t deltaEta = trigEta - assocEta;
            Double_t deltaPhi = GetDeltaPhi(trigPhi, assocPhi);

            if (trigPt <= assocPt)
                continue;

            // removing autocorrelations
            if (V0h)
            {
                Int_t negID = 0;
                Int_t posID = 0;
                Int_t atrID = 0;

                if (!hV0)
                {
                    negID = trig->GetIDNeg();
                    posID = trig->GetIDPos();
                    atrID = assoc->GetIDCh();
                }
                else
                {
                    negID = assoc->GetIDNeg();
                    posID = assoc->GetIDPos();
                    atrID = trig->GetIDCh();
                }

                if ((TMath::Abs(negID)) == (TMath::Abs(atrID)))
                    continue;
                if ((TMath::Abs(posID)) == (TMath::Abs(atrID)))
                    continue;
            }

            Int_t labelTrig = -2;
            Int_t labelAssoc = 0;
            if (hh)
            {
                labelTrig = trig->GetLabel();
                labelAssoc = assoc->GetLabel();
            }
            else if (hh)
            {
                labelTrig = trig->GetIDCh();
                labelAssoc = assoc->GetIDCh();
            }

            if (labelTrig == labelAssoc)
                continue;

            if (!hV0)
            {
                if (trig->GetOrigalPartonType() == 0)
                    trigMass = -1; // quark is the orinal parton
                else if (trig->GetOrigalPartonType() == 1)
                    trigMass = 1; // gluon is the orinal parton
                else
                    trigMass = -30;

                Double_t CorrelColl[9] = {trigPt, assocPt, deltaPhi, deltaEta, fPV[2], trig->WhichCandidate() - 0.5, trigMass, perc, (Double_t)assocCharge};
                fHistMCCorr->Fill(CorrelColl);
            }
            else
            {
                trigMass = assoc->M();
                Double_t CorrelColl[9] = {trigPt, assocPt, deltaPhi, deltaEta, fPV[2], assoc->WhichCandidate() - 0.5, trigMass, perc, (Double_t)trig->Charge()};
                fHistMCCorr->Fill(CorrelColl);
            }
        }
    }
}

///
void AliAnalysisTaskCorrelGen::CorrelationMixing(TObjArray *triggers, TObjArray *bgTracks, Float_t perc)
{
    Int_t nAssoc = bgTracks->GetEntriesFast();
    Int_t nTrig = triggers->GetEntriesFast();
    Double_t trigPt = 0.;
    Double_t trigEta = 0.;
    Double_t trigPhi = 0.;
    Double_t trigMass = 0.;
    Double_t assocPt = 0.;
    Double_t assocEta = 0.;
    Double_t assocPhi = 0.;
    Double_t assocCharge = 0.;

    AliV0ChParticleGen *trig = nullptr;
    AliVTrack *assoc = nullptr;

    for (Int_t i = 0; i < nTrig; i++)
    {
        trig = (AliV0ChParticleGen *)triggers->At(i);

        trigPt = trig->Pt();
        trigEta = trig->Eta();
        trigPhi = trig->Phi();

        if (trig->WhichCandidate() < 4)
            trigMass = trig->M();
        for (Int_t j = 0; j < nAssoc; j++)
        {
            assoc = (AliVTrack *)bgTracks->At(j);

            if (!assoc)
                continue;
            if (isnan(assoc->Eta()))
                continue;
            if (TMath::Abs(assoc->Eta()) > 0.8)
                continue;

            assocPt = assoc->Pt();
            assocEta = assoc->Eta();
            assocPhi = assoc->Phi();
            assocCharge = assoc->Charge();

            if ((assocPt >= trig->Pt()) || (assocPt < fPtAsocMin))
                continue;

            Double_t deltaEta = trigEta - assocEta;

            Double_t deltaPhi = GetDeltaPhi(trigPhi, assocPhi);

            Double_t CorrelColl[9] = {trigPt, assocPt, deltaPhi, deltaEta, fPV[2], trig->WhichCandidate() - 0.5, trigMass, perc, (Double_t)assocCharge};
            if (fMixGen)
                fHistMCMixGen->Fill(CorrelColl);
        }
    }
}

///
void AliAnalysisTaskCorrelGen::CorrelationMixinghV0(TObjArray *bgTracks, TObjArray *assocArray, Float_t perc)
{
    Int_t nAssoc = assocArray->GetEntriesFast();
    Int_t nTrig = bgTracks->GetEntriesFast();
    Double_t trigPt = 0.;
    Double_t trigEta = 0.;
    Double_t trigPhi = 0.;
    Double_t trigMass = 0.;
    Double_t assocPt = 0.;
    Double_t assocEta = 0.;
    Double_t assocPhi = 0.;

    for (Int_t i = 0; i < nTrig; i++)
    {
        AliVTrack *trig = (AliVTrack *)bgTracks->At(i);

        if (trig->Pt() < fPtTrigMin)
            continue;

        trigPt = trig->Pt();
        trigEta = trig->Eta();
        trigPhi = trig->Phi();

        if (TMath::Abs(trig->Eta()) > 0.8)
            continue;

        for (Int_t j = 0; j < nAssoc; j++)
        {
            AliV0ChParticleGen *assoc = (AliV0ChParticleGen *)assocArray->At(j);

            Double_t massAssoc = assoc->M();

            assocPt = assoc->Pt();
            assocEta = assoc->Eta();
            assocPhi = assoc->Phi();

            if ((assocPt >= trigPt) || (assocPt < fPtAsocMin))
                continue;

            Double_t deltaEta = trigEta - assocEta;

            Double_t deltaPhi = GetDeltaPhi(trigPhi, assocPhi);

            Double_t CorrelColl[9] = {trigPt, assocPt, deltaPhi, deltaEta, fPV[2], assoc->WhichCandidate() - 0.5, massAssoc, perc, (Double_t)trig->Charge()};
            if (fMixGen)
                fHistMCMixGen->Fill(CorrelColl);
        }
    }
}

///
AliAODTrack *AliAnalysisTaskCorrelGen::SetAliAODTrack(Double_t theta, Double_t phi, Double_t pt, Short_t charge)
{
    AliAODTrack *track = new AliAODTrack();
    track->SetPhi(phi);
    track->SetTheta(theta);
    track->SetPt(pt);
    track->SetCharge(charge);
    return track;
}

///
Int_t AliAnalysisTaskCorrelGen::GetOrigParton(const AliMCParticle *mcTrack) const
{
    Int_t label_mother = mcTrack->GetMother();
    Int_t pdg = mcTrack->PdgCode();

    while (label_mother > 4)
    {
        AliMCParticle *mother = static_cast<AliMCParticle *>(fmcEvent->GetTrack(label_mother));
        pdg = mother->PdgCode();

        label_mother = mother->GetMother();
    }
    return pdg;
}

/// Compute the angle between trigger and associated particles
Double_t AliAnalysisTaskCorrelGen::GetDeltaPhi(Double_t trigphi, Double_t assocphi) const
{
    Double_t dphi = trigphi - assocphi;
    if (dphi > 1.5 * TMath::Pi())
        dphi -= TMath::TwoPi();
    else if (dphi < -TMath::PiOver2())
        dphi += TMath::TwoPi();

    return dphi;
}


