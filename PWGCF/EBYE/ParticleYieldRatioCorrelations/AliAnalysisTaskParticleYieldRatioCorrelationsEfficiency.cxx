#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency.h"
#include "AliAODMCParticle.h"
#include "AliMultSelection.h"

using namespace std;

ClassImp(AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency)

    AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency::AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency() : AliAnalysisTaskSE(),
                                                                                                                         fAOD(0), fOutputList(0), fDeDx(0), fTOF(0), fHistEventsCut(0), fAliEventCuts(0), fHistTracksCut(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency::AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency(const char *name) : AliAnalysisTaskSE(name),
                                                                                                                                     fAOD(0), fOutputList(0), fDeDx(0), fTOF(0), fHistEventsCut(0), fAliEventCuts(0), fHistTracksCut(0)
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency::~AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency()
{
    // destructor
    if (fOutputList)
    {
        delete fOutputList;
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);

    fHistEventsCut = new TH1F("fHistEventsCut", ";NEvents", 6, 0, 6);
    fHistTracksCut = new TH1F("fHistTracksCut", ";Ntracks", 6, 0, 5);
    fHistEventsCut->GetXaxis()->SetBinLabel(1, "minBias");
    fHistEventsCut->GetXaxis()->SetBinLabel(2, "centrality");
    fHistEventsCut->GetXaxis()->SetBinLabel(3, "NContributors");
    fHistEventsCut->GetXaxis()->SetBinLabel(4, "vertex");
    fHistEventsCut->GetXaxis()->SetBinLabel(5, "AliEventCuts");
    fHistTracksCut->GetXaxis()->SetBinLabel(1, "minBias");
    fHistTracksCut->GetXaxis()->SetBinLabel(2, "FilterBit");
    fHistTracksCut->GetXaxis()->SetBinLabel(3, "Eta");
    fHistTracksCut->GetXaxis()->SetBinLabel(4, "pT");
    fHistTracksCut->GetXaxis()->SetBinLabel(5, "TOF");
    fDeDx = new TH2D("DeDx", ";p_{TPC};dE/dx (a.u.)", 200, 0.2, 2, 250, 0, 2000);
    fTOF = new TH2D("fTOF", ";p_{TPC};t", 200, 0.2, 2, 200, -4000, 4000);
    fOutputList->Add(fDeDx);
    fOutputList->Add(fTOF);
    for (int iSort = 0; iSort < 4; iSort++)
    {
        fDeDxSorts[iSort] = new TH2D(Form("DeDxSort%d", iSort), ";p_{TPC};dE/dx (a.u.)", 200, 0.2, 2, 250, 0, 2000);
        fTOFSorts[iSort] = new TH2D(Form("fTOFSort%d", iSort), ";p_{TPC};t", 200, 0.2, 2, 200, -4000, 4000);
        fOutputList->Add(fDeDxSorts[iSort]);
        fOutputList->Add(fTOFSorts[iSort]);
    }
    fOutputList->Add(fHistEventsCut);
    fOutputList->Add(fHistTracksCut);
nPBins = 5;
Double_t PBins[6] = {0.2, 0.3, 0.4, 0.6, 0.8, 2.0};
    Double_t PhiBins[nPhiBins + 1], VertexBins[nVertexBins + 1];
    for (int i = 0; i < nPhiBins + 1; i++)
        PhiBins[i] = i * TMath::TwoPi() / nPhiBins;
    for (int i = 0; i < nVertexBins + 1; i++)
        VertexBins[i] = Vertexmin + i * (Vertexmax - Vertexmin) / nVertexBins;

    for (int iCent = 0; iCent < nCentrClassesUsed; iCent++)
    {
        for (int iEta = 0; iEta < nEtaClasses; iEta++)
        {
            for (int iSort = 0; iSort < nSorts; iSort++)
            {
                NGenTracks[iCent][iEta][iSort] = new TH3D(Form("NGenTracksC%dEta%dSort%d", iCent, iEta, iSort), ";P_{T};Ntracks",
                                                          nPBins, PBins, nPhiBins, PhiBins, nVertexBins, VertexBins);
                EfficiencyTracking[iCent][iEta][iSort] = new TH3D(Form("EfficiencyTrackingC%dEta%dSort%d", iCent, iEta, iSort), ";P_{T};EfficiencyTracking",
                                                                  nPBins, PBins, nPhiBins, PhiBins, nVertexBins, VertexBins);
                ContaminationTracking[iCent][iEta][iSort] = new TH3D(Form("ContaminationTrackingC%dEta%dSort%d", iCent, iEta, iSort), ";P_{T};ContaminationTracking",
                                                                     nPBins, PBins, nPhiBins, PhiBins, nVertexBins, VertexBins);
                NReconstTracks[iCent][iEta][iSort] = new TH3D(Form("NReconstTracksC%dEta%dSort%d", iCent, iEta, iSort), ";P_{T};Ntracks",
                                                              nPBins, PBins, nPhiBins, PhiBins, nVertexBins, VertexBins);

                NTrueTracks[iCent][iEta][iSort] = new TH3D(Form("NTrueTracksC%dEta%dSort%d", iCent, iEta, iSort), ";P_{T};Ntracks",
                                                           nPBins, PBins, nPhiBins, PhiBins, nVertexBins, VertexBins);
                EfficiencyPID[iCent][iEta][iSort] = new TH3D(Form("EfficiencyPIDC%dEta%dSort%d", iCent, iEta, iSort), ";P_{T};EfficiencyPID",
                                                             nPBins, PBins, nPhiBins, PhiBins, nVertexBins, VertexBins);
                ContaminationPID[iCent][iEta][iSort] = new TH3D(Form("ContaminationPIDC%dEta%dSort%d", iCent, iEta, iSort), ";P_{T};ContaminationPID",
                                                                nPBins, PBins, nPhiBins, PhiBins, nVertexBins, VertexBins);
                NTracksInCut[iCent][iEta][iSort] = new TH3D(Form("NTracksInCutC%dEta%dSort%d", iCent, iEta, iSort), ";P_{T};Ntracks",
                                                            nPBins, PBins, nPhiBins, PhiBins, nVertexBins, VertexBins);

                fOutputList->Add(NGenTracks[iCent][iEta][iSort]);
                fOutputList->Add(EfficiencyTracking[iCent][iEta][iSort]);
                fOutputList->Add(ContaminationTracking[iCent][iEta][iSort]);
                fOutputList->Add(NReconstTracks[iCent][iEta][iSort]);

                fOutputList->Add(NTrueTracks[iCent][iEta][iSort]);
                fOutputList->Add(EfficiencyPID[iCent][iEta][iSort]);
                fOutputList->Add(ContaminationPID[iCent][iEta][iSort]);
                fOutputList->Add(NTracksInCut[iCent][iEta][iSort]);
            }
        }
    }

    for (int iSort = 0; iSort < nSorts; iSort++)
    {
        purityAll[iSort] = new TH1D(Form("pur%d", iSort), "", nPBins, PBins);
        fOutputList->Add(purityAll[iSort]);
        for (int jSort = 0; jSort < nSorts; jSort++)
        {
            purity[iSort][jSort] = new TH1D(Form("pur%din%d", jSort, iSort), "", nPBins, PBins);
            fOutputList->Add(purity[iSort][jSort]);
        }
    }
    for (Int_t i(0); i < 3; i++)
    {
        fHistQASPDTrackletsvsV0MCent[i] = new TH2D(Form("fHistQASPDTrackletsvsV0MCent%d", i), ";V0M Percentile;N Tracklets in SPD", 100, 0, 100, 400, 0, 1e4);
        fOutputList->Add(fHistQASPDTrackletsvsV0MCent[i]);
    }
    for (Int_t i(0); i < 2; i++)
    {
        fHistQAMultTPCvsESD[i] = new TH2D(Form("fHistQAMultTPCvsESD%d", i), ";MultTPC;MultESD", 400, 0, 1e4, 400, 0, 5e4);
        fOutputList->Add(fHistQAMultTPCvsESD[i]);
        fHistQAMultTPCvsV0[i] = new TH2D(Form("fHistQAMultTPCvsV0%d", i), ";MultTPC;V0", 400, 0, 1e4, 400, 0, 5e4);
        fOutputList->Add(fHistQAMultTPCvsV0[i]);
        fHistQAMultTrkvsMultTrkTOF[i] = new TH2D(Form("fHistQAMultTrkvsMultTrkTOF%d", i), ";MultTrk;MultTrkTOF", 400, 0, 5e3, 400, 0, 5e3);
        fOutputList->Add(fHistQAMultTrkvsMultTrkTOF[i]);
    }
    fAliEventCuts = new AliEventCuts();
    if (pbpb)
        fAliEventCuts->SetupPbPb2018();
    else
        fAliEventCuts->SetupRun2pp();
    fAliEventCuts->AddQAplotsToList(fOutputList);
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if (man)
    {
        AliInputEventHandler *inputHandler = dynamic_cast<AliInputEventHandler *>(man->GetInputEventHandler());
        if (inputHandler)
            fPIDResponse = inputHandler->GetPIDResponse();
        else
            AliFatal("Input handler needed");
    }
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());

    if (!fAOD)
    {
        PostData(1, fOutputList);
        return;
    }
    fHistEventsCut->Fill("minBias", 1);
    if (!fAliEventCuts->AcceptEvent(fAOD))
    {
        PostData(1, fOutputList);
        return;
    }

    Float_t centr;
    Int_t NTrackletsSPD;
    if (0)
    {
        AliCentrality *centrality = fAOD->GetCentrality();
        centr = centrality->GetCentralityPercentile("V0M");
    }
    else
    {
        AliMultSelection *MultSelection = 0x0;
        MultSelection = (AliMultSelection *)fAOD->FindListObject("MultSelection");
        if (MultSelection)
        {
            centr = MultSelection->GetMultiplicityPercentile("V0M");
            NTrackletsSPD = MultSelection->GetEstimator("SPDTracklets")->GetValue();
        }
        else
        {
            fHistEventsCut->Fill("AliEventCuts", 1);
            PostData(1, fOutputList);
            return;
        }
    }
    // centrality cut:
    if (centr < minCent || centr > maxCent)
    {
        PostData(1, fOutputList);
        return;
    }
    fHistEventsCut->Fill("centrality", 1);
    fHistQASPDTrackletsvsV0MCent[0]->Fill(centr, NTrackletsSPD);
    if (SPDvsV0MCut)
    {
        TF1 *fSPDvsV0M_DownLimit = new TF1("fSPDvsV0M_DownLimit", "exp(8.1456-0.0354*x-3.8e-04*x*x)", 0, 80);
        TF1 *fSPDvsV0M_UperLimit = new TF1("fSPDvsV0M_UperLimit", "exp(8.58-0.036*x-4.4e-05*x*x)", 0, 80);
        if (NTrackletsSPD < fSPDvsV0M_DownLimit->Eval(centr) || NTrackletsSPD > fSPDvsV0M_UperLimit->Eval(centr))
        {
            PostData(1, fOutputList);
            return;
        }
        fHistEventsCut->Fill("SPDvsV0M", 1);
    }
    fHistQASPDTrackletsvsV0MCent[1]->Fill(centr, NTrackletsSPD);
    // Events with large number of TPC clusters
    Int_t multEsd = ((AliAODHeader *)fAOD->GetHeader())->GetNumberOfESDTracks();
    Int_t multTPC = 0;
    const AliAODVertex *vtx = fAOD->GetPrimaryVertex();
    if (vtx->GetNContributors() < 1)
    {
        PostData(1, fOutputList);
        return;
    }
    fHistEventsCut->Fill("NContributors", 1);
    if (fabs((Float_t)vtx->GetZ()) > 8)
    {
        PostData(1, fOutputList);
        return;
    }
    fHistEventsCut->Fill("vertex", 1);
    Float_t vertex = (Float_t)vtx->GetZ();
    Int_t nTracks(fAOD->GetNumberOfTracks());
    Int_t EtaBin, CentrBin;
    int sort = 20;
    // CentrBin = (centr - minCent) / ((maxCent - minCent) / nCentrClasses);
    // CentrBin = centr<=5?0:(centr<=10?1: 2+(centr-10)/((maxCent-10)/(nCentrClasses-2)) ); //if first bins are 0-5 & 5-10
    for (Int_t i(0); i < nCentrClassesUsed; i++)
    {
        if (centr >= CentrPercentiles[i] && centr <= CentrPercentiles[i + 1])
        {
            CentrBin = i;
        }
    }
    if (CentrBin < 0 || CentrBin >= nCentrClassesUsed)
    {
        PostData(1, fOutputList);
        return;
    }
    int NAcceptedtracks = 0;
    fHistTracksCut->Fill("minBias", nTracks);

    AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    AliMCEvent *fMC = eventHandler->MCEvent();
    Int_t nMCTracks(fMC->GetNumberOfTracks());
    Float_t PtCut[3] = {0.2, 0.5, 0.5};
    Float_t nSigmaBoundary[3] = {0.5, 0.32, 0.7};
    Int_t multTrk = 0;
    Int_t multTrkTOF = 0;
    // Track loop:
    for (Int_t i(0); i < nTracks; i++)
    {
        AliAODTrack *track = static_cast<AliAODTrack *>(fAOD->GetTrack(i));
        if (track->TestFilterBit(32))
        {
            multTrk++;
            if (TMath::Abs(track->GetTOFsignalDz()) <= 10 && track->GetTOFsignal() >= 12000 && track->GetTOFsignal() <= 25000)
                multTrkTOF++;
        }
        if (track->TestFilterBit(128))
            multTPC++;
        if (!track || !track->TestFilterBit(filterBit))
            continue;
        fHistTracksCut->Fill("FilterBit", 1);
        if (fabs(track->Eta()) >= 0.8)
            continue;
        fHistTracksCut->Fill("Eta", 1);
        if (track->Pt() < minP || track->Pt() > maxP)
            continue;
        fHistTracksCut->Fill("pT", 1);
        EtaBin = (0.8 + track->Eta()) / (1.6 / nEtaClasses);

        Float_t nOfSigmasTPC_el = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        Float_t nOfSigmasTPC_pi = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        Float_t nOfSigmasTPC_K = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        Float_t nOfSigmasTPC_p = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        Float_t nOfSigmasTPC_d = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);

        Float_t nOfSigmasTOF_pi = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion, fPIDResponse->GetTOFResponse().GetTimeZero());
        Float_t nOfSigmasTOF_K = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon, fPIDResponse->GetTOFResponse().GetTimeZero());
        Float_t nOfSigmasTOF_p = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton, fPIDResponse->GetTOFResponse().GetTimeZero());

        if (fabs(nOfSigmasTOF_pi) > 900)
        {
            fHistTracksCut->Fill("TOF", 1);
            // continue;
        }
        /*if (fabs(nOfSigmasTOF_K) > 900)
        {
            continue;
        }
        if (fabs(nOfSigmasTOF_p) > 900)
        {
            continue;
        }*/
        if (track->GetTPCCrossedRows() <= nCrossedRows)
        {
            continue;
        }

        Float_t Pt = track->Pt();
        Float_t Px = track->Px();
        Float_t Py = track->Py();
        Float_t Pz = track->Pz();
        Float_t Moment = track->GetTPCmomentum();
        Float_t Phi = track->Phi();
        Float_t Eta = track->Eta();
        Float_t DeDx = track->GetTPCsignal();
        int Charge = track->Charge();
        fDeDx->Fill(Moment, DeDx);
        fTOF->Fill(Moment, track->GetTOFsignal());

        sort = 20;
        if (IsMC && fMC)
        {
            int label = track->GetLabel();
            AliAODMCParticle *part = (AliAODMCParticle *)fMC->GetTrack(abs(label));
            if (part)
            {
                if (fabs(part->GetPdgCode()) == 211)
                    sort = 0;
                if (fabs(part->GetPdgCode()) == 321)
                    sort = 1;
                if (fabs(part->GetPdgCode()) == 2212)
                    sort = 2;
                if (fabs(part->GetPdgCode()) == 11)
                    sort = 3;
                if (fabs(part->GetPdgCode()) == 1000010020)
                    sort = 4;
                if (sort > 3)
                    continue;
                if (Pt > PtCut[sort])
                {
                    if (part->IsPhysicalPrimary())
                        EfficiencyTracking[CentrBin][EtaBin][sort * 2 + (Charge < 0 ? 0 : 1)]->Fill(Pt, Phi, vertex);
                    else
                        ContaminationTracking[CentrBin][EtaBin][sort * 2 + (Charge < 0 ? 0 : 1)]->Fill(Pt, Phi, vertex);
                    NReconstTracks[CentrBin][EtaBin][sort * 2 + (Charge < 0 ? 0 : 1)]->Fill(Pt, Phi, vertex);
                }
            }
        }

        Float_t nSigma_comb_pi = sqrt(nOfSigmasTPC_pi * nOfSigmasTPC_pi + nOfSigmasTOF_pi * nOfSigmasTOF_pi);
        Float_t nSigma_comb_K = sqrt(nOfSigmasTPC_K * nOfSigmasTPC_K + nOfSigmasTOF_K * nOfSigmasTOF_K);
        Float_t nSigma_comb_p = sqrt(nOfSigmasTPC_p * nOfSigmasTPC_p + nOfSigmasTOF_p * nOfSigmasTOF_p);

        // Fill TH3D hists for maps
        bool selected_pi = false;
        if (Pt < nSigmaBoundary[0] && fabs(nOfSigmasTPC_pi) < nSigma && fabs(nOfSigmasTPC_K) > 3 && fabs(nOfSigmasTPC_p) > 3 && fabs(nOfSigmasTPC_el) > 1)
            selected_pi = true;
        if (Pt > nSigmaBoundary[0] && nSigma_comb_pi < nSigma && nSigma_comb_K > 3 && nSigma_comb_p > 3)
            selected_pi = true;
        if (selected_pi)
        {
            if (Charge < 0)
            {
                NTracksInCut[CentrBin][EtaBin][0]->Fill(Pt, Phi, vertex);
                if (sort == 0)
                    EfficiencyPID[CentrBin][EtaBin][0]->Fill(Pt, Phi, vertex);
                else
                    ContaminationPID[CentrBin][EtaBin][0]->Fill(Pt, Phi, vertex);
                purity[0][sort * 2 + (Charge < 0 ? 0 : 1)]->Fill(Pt);
                purityAll[0]->Fill(Pt);
            }
            if (Charge > 0)
            {
                NTracksInCut[CentrBin][EtaBin][1]->Fill(Pt, Phi, vertex);
                if (sort == 0)
                    EfficiencyPID[CentrBin][EtaBin][1]->Fill(Pt, Phi, vertex);
                else
                    ContaminationPID[CentrBin][EtaBin][1]->Fill(Pt, Phi, vertex);
                purity[1][sort * 2 + (Charge < 0 ? 0 : 1)]->Fill(Pt);
                purityAll[1]->Fill(Pt);
            }
            fDeDxSorts[0]->Fill(Moment, DeDx);
            fTOFSorts[0]->Fill(Moment, track->GetTOFsignal());
        }

        bool selected_K = false;
        if (Pt < nSigmaBoundary[1] && fabs(nOfSigmasTPC_K) < nSigma && fabs(nOfSigmasTPC_pi) > 3 && fabs(nOfSigmasTPC_p) > 3) // && fabs( nOfSigmasTPC_el)>1 )
            selected_K = true;
        if (Pt > nSigmaBoundary[1] && nSigma_comb_K < nSigma && nSigma_comb_pi > 3 && nSigma_comb_p > 3)
            selected_K = true;
        if (selected_K && Pt > PtCut[1])
        {
            if (Charge < 0)
            {
                NTracksInCut[CentrBin][EtaBin][2]->Fill(Pt, Phi, vertex);
                if (sort == 1)
                    EfficiencyPID[CentrBin][EtaBin][2]->Fill(Pt, Phi, vertex);
                else
                    ContaminationPID[CentrBin][EtaBin][2]->Fill(Pt, Phi, vertex);
                purity[2][sort * 2 + (Charge < 0 ? 0 : 1)]->Fill(Pt);
                purityAll[2]->Fill(Pt);
            }
            if (Charge > 0)
            {
                NTracksInCut[CentrBin][EtaBin][3]->Fill(Pt, Phi, vertex);
                if (sort == 1)
                    EfficiencyPID[CentrBin][EtaBin][3]->Fill(Pt, Phi, vertex);
                else
                    ContaminationPID[CentrBin][EtaBin][3]->Fill(Pt, Phi, vertex);
                purity[3][sort * 2 + (Charge < 0 ? 0 : 1)]->Fill(Pt);
                purityAll[3]->Fill(Pt);
            }
            fDeDxSorts[1]->Fill(Moment, DeDx);
            fTOFSorts[1]->Fill(Moment, track->GetTOFsignal());
        }

        bool selected_p = false;
        if (Pt < nSigmaBoundary[2] && fabs(nOfSigmasTPC_p) < nSigma && fabs(nOfSigmasTPC_pi) > 3 && fabs(nOfSigmasTPC_K) > 3 && fabs(nOfSigmasTPC_el) > 1 && fabs(nOfSigmasTPC_d) > 3)
            selected_p = true;
        if (Pt > nSigmaBoundary[2] && nSigma_comb_p < nSigma && nSigma_comb_pi > 3 && nSigma_comb_K > 3)
            selected_p = true;
        if (selected_p && Pt > PtCut[2])
        {
            if (Charge < 0)
            {
                NTracksInCut[CentrBin][EtaBin][4]->Fill(Pt, Phi, vertex);
                if (sort == 2)
                    EfficiencyPID[CentrBin][EtaBin][4]->Fill(Pt, Phi, vertex);
                else
                    ContaminationPID[CentrBin][EtaBin][4]->Fill(Pt, Phi, vertex);
                purity[4][sort * 2 + (Charge < 0 ? 0 : 1)]->Fill(Pt);
                purityAll[4]->Fill(Pt);
            }
            if (Charge > 0)
            {
                NTracksInCut[CentrBin][EtaBin][5]->Fill(Pt, Phi, vertex);
                if (sort == 2)
                    EfficiencyPID[CentrBin][EtaBin][5]->Fill(Pt, Phi, vertex);
                else
                    ContaminationPID[CentrBin][EtaBin][5]->Fill(Pt, Phi, vertex);
                purity[5][sort * 2 + (Charge < 0 ? 0 : 1)]->Fill(Pt);
                purityAll[5]->Fill(Pt);
            }
            fDeDxSorts[2]->Fill(Moment, DeDx);

            fTOFSorts[2]->Fill(Moment, track->GetTOFsignal());
        }

        bool selected_el = false;
        if (fabs(nOfSigmasTPC_p) > 3 && fabs(nOfSigmasTPC_pi) > 3 && fabs(nOfSigmasTPC_K) > 3 && fabs(nOfSigmasTPC_el) < 2)
            selected_el = true;
        if (selected_el)
        {
            if (Charge < 0)
            {
                NTracksInCut[CentrBin][EtaBin][6]->Fill(Pt, Phi, vertex);
                if (sort == 3)
                    EfficiencyPID[CentrBin][EtaBin][6]->Fill(Pt, Phi, vertex);
                else
                    ContaminationPID[CentrBin][EtaBin][6]->Fill(Pt, Phi, vertex);
                purity[6][sort * 2 + (Charge < 0 ? 0 : 1)]->Fill(Pt);
                purityAll[6]->Fill(Pt);
            }
            if (Charge > 0)
            {
                NTracksInCut[CentrBin][EtaBin][7]->Fill(Pt, Phi, vertex);
                if (sort == 3)
                    EfficiencyPID[CentrBin][EtaBin][7]->Fill(Pt, Phi, vertex);
                else
                    ContaminationPID[CentrBin][EtaBin][7]->Fill(Pt, Phi, vertex);
                purity[7][sort * 2 + (Charge < 0 ? 0 : 1)]->Fill(Pt);
                purityAll[7]->Fill(Pt);
            }
            fDeDxSorts[3]->Fill(Moment, DeDx);
            fTOFSorts[3]->Fill(Moment, track->GetTOFsignal());
        }
        if (sort > 3)
            continue;
        if (Pt > PtCut[sort])
            NTrueTracks[CentrBin][EtaBin][sort * 2 + (Charge < 0 ? 0 : 1)]->Fill(Pt, Phi, vertex);
        // end of filling TH3D hists for maps

        NAcceptedtracks++;
    }
    // end of track loop
    const AliAODVZERO *vzrData = fAOD->GetVZEROData();
    float sumV0ampl = vzrData->GetMTotV0A() + vzrData->GetMTotV0C();
    fHistQAMultTPCvsV0[0]->Fill(multTPC, sumV0ampl);
    fHistQAMultTPCvsESD[0]->Fill(multTPC, multEsd);
    fHistQAMultTrkvsMultTrkTOF[0]->Fill(multTrk, multTrkTOF);
    if (multEsd - 3.38 * multTPC >= 15000 && LargeTPCCut)
    {
        PostData(1, fOutputList);
        return;
    }
    fHistQASPDTrackletsvsV0MCent[2]->Fill(centr, NTrackletsSPD);
    fHistQAMultTPCvsESD[1]->Fill(multTPC, multEsd);
    fHistQAMultTPCvsV0[1]->Fill(multTPC, multEsd);
    fHistQAMultTrkvsMultTrkTOF[1]->Fill(multTrk, multTrkTOF);

    int nAcceptedGenTracks = 0;
    if (IsMC && fMC)
    {
        // Generated particles track loop
        for (Int_t i(0); i < nMCTracks; i++)
        {
            AliAODMCParticle *trackMC = (AliAODMCParticle *)fMC->GetTrack(i);
            if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC))
                continue;
            if (trackMC && trackMC->IsPhysicalPrimary() && fabs(trackMC->Eta()) < 0.8 && trackMC->Pt() >= minP && trackMC->Pt() <= maxP)
            {
                Float_t GenPt = trackMC->Pt();
                Float_t GenPhi = trackMC->Phi();
                Float_t GenEta = trackMC->Eta();
                int GenCharge = trackMC->Charge();

                sort = 20;
                if (fabs(trackMC->GetPdgCode()) == 211)
                    sort = 0;
                if (fabs(trackMC->GetPdgCode()) == 321)
                    sort = 1;
                if (fabs(trackMC->GetPdgCode()) == 2212)
                    sort = 2;
                if (fabs(trackMC->GetPdgCode()) == 11)
                    sort = 3;
                if (fabs(trackMC->GetPdgCode()) == 1000010020)
                    sort = 4;
                nAcceptedGenTracks++;

                if ((sort >= 0 && sort <= 3) && GenPt > PtCut[sort] && GenCharge != 0)
                {
                    EtaBin = (0.8 + GenEta) / (1.6 / nEtaClasses);
                    NGenTracks[CentrBin][EtaBin][sort * 2 + (GenCharge < 0 ? 0 : 1)]->Fill(GenPt, GenPhi, vertex);
                }
            }
        }
    }
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
