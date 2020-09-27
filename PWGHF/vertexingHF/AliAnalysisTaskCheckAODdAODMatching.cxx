/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. */

/////////////////////////////////////////////////////////////////////////////////////////
// \class AliAnalysisTaskCheckAODdAODMatching                                          //
// \brief analysis task for tagging the AOD files with AOD-dAOD mismatches             //
// \author: F. Grosa, fabrizio.grosa@cern.ch                                           //
// \author: F. Prino, francesco.prino@cern.ch                                          //
// \author: C. Terrevoli, cristina.terrevoli@cern.ch                                   //
/////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <TMath.h>
#include <TFile.h>
#include <TDatabasePDG.h>
#include <TClonesArray.h>
#include <TGrid.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliESDtrackCuts.h"
#include "AliRDHFCuts.h"
#include "AliESDtrack.h"

#include "AliAnalysisTaskCheckAODdAODMatching.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCheckAODdAODMatching);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskCheckAODdAODMatching::AliAnalysisTaskCheckAODdAODMatching() : AliAnalysisTaskSE()
    ,fOutput(nullptr)
    ,fHistFiles(nullptr)
    ,fTreeMismatch(nullptr)
    ,fAOD(nullptr)
    ,fAODMap{}
    ,fPrevInputFileName("")
    ,fStatus(0)
    ,fRunNumber(-1)
    ,fNevents(-1)
{
    /// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskCheckAODdAODMatching::AliAnalysisTaskCheckAODdAODMatching(const char *name) : AliAnalysisTaskSE(name)
    ,fOutput(nullptr)
    ,fHistFiles(nullptr)
    ,fTreeMismatch(nullptr)
    ,fAOD(nullptr)
    ,fAODMap{}
    ,fPrevInputFileName("")
    ,fStatus(0)
    ,fRunNumber(-1)
    ,fNevents(-1)
{
    /// Standard constructor

    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskCheckAODdAODMatching::~AliAnalysisTaskCheckAODdAODMatching()
{
    /// Destructor
    delete fOutput;
    delete fTreeMismatch;
}

//________________________________________________________________________
void AliAnalysisTaskCheckAODdAODMatching::UserCreateOutputObjects()
{
    /// Create the output container
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("fOutput");

    fHistFiles = new TH1F("fHistFiles", ";;number of files", 4, 0.5, 4.5);
    fHistFiles->GetXaxis()->SetBinLabel(1, "good files");
    fHistFiles->GetXaxis()->SetBinLabel(2, "AOD-dAOD mismatch num of event");
    fHistFiles->GetXaxis()->SetBinLabel(3, "AOD-dAOD mismatch TProcessID");
    fHistFiles->GetXaxis()->SetBinLabel(4, "AOD-dAOD mismatch from cand");
    fOutput->Add(fHistFiles);

    PostData(1, fOutput);

    fTreeMismatch = new TTree("fTreeMismatch", "fTreeMismatch");
    fTreeMismatch->Branch("file_name", &fPrevInputFileName);
    fTreeMismatch->Branch("mismatch_status", &fStatus);
    fTreeMismatch->Branch("n_events", &fNevents);
    fTreeMismatch->Branch("run_number", &fRunNumber);
    OpenFile(2);
        PostData(2, fTreeMismatch);
}

//________________________________________________________________________
void AliAnalysisTaskCheckAODdAODMatching::UserExec(Option_t * /*option*/)
{
    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());
    AliAODHandler *aodHandler = dynamic_cast<AliAODHandler *>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());

    TClonesArray *arrayD0toKpi = nullptr, *array3Prong = nullptr, *arrayDstar = nullptr, *arrayCasc = nullptr, *arrayLS2p = nullptr, *arrayLS3p = nullptr;
    if (!fAOD && AODEvent() && IsStandardAOD())
    {
        if (aodHandler->GetExtensions())
        {
            AliAODExtension *ext = (AliAODExtension *)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
            AliAODEvent *aodFromExt = ext->GetAOD();
            arrayD0toKpi = (TClonesArray *)aodFromExt->GetList()->FindObject("D0toKpi");
            array3Prong = (TClonesArray *)aodFromExt->GetList()->FindObject("Charm3Prong");
            arrayDstar = (TClonesArray *)aodFromExt->GetList()->FindObject("Dstar");
            arrayCasc = (TClonesArray *)aodFromExt->GetList()->FindObject("CascadesHF");
            arrayLS2p = (TClonesArray *)aodFromExt->GetList()->FindObject("LikeSign2Prong");
            arrayLS3p = (TClonesArray *)aodFromExt->GetList()->FindObject("LikeSign3Prong");
        }
    }
    else if (fAOD)
    {
        arrayD0toKpi = (TClonesArray *)fAOD->GetList()->FindObject("D0toKpi");
        array3Prong = (TClonesArray *)fAOD->GetList()->FindObject("Charm3Prong");
        arrayDstar = (TClonesArray *)fAOD->GetList()->FindObject("Dstar");
        arrayCasc = (TClonesArray *)fAOD->GetList()->FindObject("CascadesHF");
        arrayLS2p = (TClonesArray *)fAOD->GetList()->FindObject("LikeSign2Prong");
        arrayLS3p = (TClonesArray *)fAOD->GetList()->FindObject("LikeSign3Prong");
    }

    fStatus = 0;

    if(!fAOD || (!arrayD0toKpi && !array3Prong && !arrayDstar && !arrayCasc && !arrayLS2p && !arrayLS3p))
    {
        return;
    }

    fRunNumber = fAOD->GetRunNumber();

    // get name of analysed AliAOD.root
    AliAODInputHandler *aodInputHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    TTree *treeAOD = aodInputHandler->GetTree();
    TString currentFile = treeAOD->GetCurrentFile()->GetName();
    // analyse only one event per AliAOD.root (in any case loop over all the events of a given AliAOD.root in each UserExec loop)
    if (currentFile.EqualTo(fPrevInputFileName))
    {
        return;
    }

    bool okNevents = true;
    bool okTProcessID = true;
    int mismatchCheck = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if( mismatchCheck == -1)
    {
        okNevents = false;
        fHistFiles->Fill(2);
        fStatus |= kMismatchNevents;
    }
    else if(mismatchCheck == 0)
    {
        okTProcessID = false;
        fHistFiles->Fill(3);
        fStatus |= kMismatchTProcessID;
    }

    treeAOD = nullptr;
    if(!gGrid)
        TGrid::Connect("alien://");
    TFile* file = TFile::Open(currentFile.Data());
    treeAOD = (TTree*)file->Get("aodTree");
    fNevents = treeAOD->GetEntries();

    int nTotD0toKpi = 0, nTotDstar = 0, nTot3Prong = 0; //nTotHF = 0, 
    int nTotDplus = 0, nTotDs = 0, nTotLc = 0, nTotCasc = 0;
    int nTotLcV0 = 0, nTotDsV0 = 0, nTotDplusV0 = 0;
    int nTotLS2p = 0, nTotLS3p = 0;
    int nTotD0toKpiCuts = 0; //, nTotDplusCuts = 0, nTotDsCuts = 0, nTotLcCuts = 0;
    int nBad = 0;

    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "default");
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetMinNClustersTPC(50);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    // |d0|>25 micron for pt<2GeV, no cut above 2
    esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0025*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
    esdTrackCuts->SetMaxDCAToVertexXY(1.);
    esdTrackCuts->SetMaxDCAToVertexZ(1.);
    esdTrackCuts->SetPtRange(0.4, 1.e10);
    esdTrackCuts->SetEtaRange(-0.8, +0.8);

    double massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    double massK = TDatabasePDG::Instance()->GetParticle(321)->Mass();
    double massD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    bool okCand = true;
    for (int nEv = 0; nEv < fNevents; nEv++)
    {
        treeAOD->GetEvent(nEv);
        MapAODtracks();
        if (arrayD0toKpi)
        {
            int nD0toKpi = arrayD0toKpi->GetEntriesFast();
            nTotD0toKpi += nD0toKpi;
            for (int id0 = 0; id0 < nD0toKpi; id0++)
            {
                AliAODRecoDecayHF2Prong *cc = (AliAODRecoDecayHF2Prong *)arrayD0toKpi->UncheckedAt(id0);
                if (cc && cc->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts))
                {
                    nTotD0toKpiCuts++;
                    int ch[2];
                    double e2k[2], e2pi[2];
                    double sumPx = 0, sumPy = 0, sumPz = 0;
                    for (int iDau = 0; iDau < 2; iDau++)
                    {
                        int dau = cc->GetProngID(iDau);
                        AliAODTrack *track = (AliAODTrack *)fAOD->GetTrack(fAODMap[dau]);
                        if (!track)
                        {
                            okCand = false;
                            continue;
                        }
                        if (track->Pt() < 0.4)
                        {
                            okCand = false;
                            continue;
                        }
                        track->SetAODEvent(fAOD);
                        if (!esdTrackCuts->IsSelected(track))
                        {
                            okCand = false;
                            continue;
                        }
                        ch[iDau] = track->Charge();
                        double p2 = track->Px() * track->Px() + track->Py() * track->Py() + track->Pz() * track->Pz();
                        e2k[iDau] = massK * massK + p2;
                        e2pi[iDau] = massPi * massPi + p2;
                        sumPx += track->Px();
                        sumPy += track->Py();
                        sumPz += track->Pz();
                    }
                    if (ch[0] == ch[1])
                    {
                        okCand = false;
                        continue;
                    }
                    double esum1 = TMath::Sqrt(e2k[0]) + TMath::Sqrt(e2pi[1]);
                    double esum2 = TMath::Sqrt(e2k[1]) + TMath::Sqrt(e2pi[0]);
                    double p2 = sumPx * sumPx + sumPy * sumPy + sumPz * sumPz;
                    double ptD = TMath::Sqrt(sumPx * sumPx + sumPy * sumPy);
                    double im1 = TMath::Sqrt(esum1 * esum1 - p2);
                    double im2 = TMath::Sqrt(esum2 * esum2 - p2);
                    double massCut = 0.25;
                    if (ptD > 5)
                        massCut = 0.4;
                    massCut *= 1.01; // numerical tolerance...
                    if (TMath::Abs(im1 - massD0) > massCut && TMath::Abs(im2 - massD0) > massCut)
                    {
                        okCand = false;
                        continue;
                    }
                }
            }
        }
        if (array3Prong)
        {
            int n3Prong = array3Prong->GetEntriesFast();
            nTot3Prong += n3Prong;
            for (int i3p = 0; i3p < n3Prong; i3p++)
            {
                AliAODRecoDecayHF3Prong *ccc = (AliAODRecoDecayHF3Prong *)array3Prong->UncheckedAt(i3p);
                bool good = false;
                if (ccc->HasSelectionBit(AliRDHFCuts::kDplusCuts))
                {
                    nTotDplus++;
                    good = true;
                }
                if (ccc->HasSelectionBit(AliRDHFCuts::kDsCuts))
                {
                    nTotDs++;
                    good = true;
                }
                if (ccc->HasSelectionBit(AliRDHFCuts::kLcCuts))
                {
                    nTotLc++;
                    good = true;
                }
                if (good)
                {
                    int ch[3];
                    for (int iDau = 0; iDau < 3; iDau++)
                    {
                        int dau = ccc->GetProngID(iDau);
                        AliAODTrack *track = (AliAODTrack *)fAOD->GetTrack(fAODMap[dau]);
                        if (!track)
                        {
                            okCand = false;
                            continue;
                        }
                        if (track->Pt() < 0.4)
                        {
                            okCand = false;
                            continue;
                        }
                        track->SetAODEvent(fAOD);
                        if (!esdTrackCuts->IsSelected(track))
                        {
                            okCand = false;
                            continue;
                        }
                        ch[iDau] = track->Charge();
                    }
                    if (ch[1] == ch[0] || ch[1] == ch[2])
                    {
                        okCand = false;
                        continue;
                    }
                }
                else
                {
                    nBad++;
                }
            }
        }
        if (arrayDstar) // TODO: add checks on Dstar and HF cascades
        {
            int nDstar = arrayDstar->GetEntriesFast();
            nTotDstar += nDstar;
        }
        if (arrayCasc)
        {
            int nCasc = arrayCasc->GetEntriesFast();
            nTotCasc += nCasc;
            for (int ica = 0; ica < nCasc; ica++)
            {
                AliAODRecoCascadeHF *d = (AliAODRecoCascadeHF *)arrayCasc->UncheckedAt(ica);
                if (d->HasSelectionBit(AliRDHFCuts::kDstoK0sCuts))
                {
                    nTotDsV0++;
                }
                else if (d->HasSelectionBit(AliRDHFCuts::kLctoV0Cuts))
                {
                    nTotLcV0++;
                }
                else if (d->HasSelectionBit(AliRDHFCuts::kDplustoK0sCuts))
                {
                    nTotDplusV0++;
                }
            }
        }
        if (arrayLS2p)
        {
            int nLS2p = arrayLS2p->GetEntriesFast();
            nTotLS2p += nLS2p;
        }
        if (arrayLS3p)
        {
            int nLS3p = arrayLS3p->GetEntriesFast();
            nTotLS3p += nLS3p;
        }
    }

    if(!okCand)
    {
        fHistFiles->Fill(4);
        fStatus |= kMismatchCand;
    }

    if(okNevents && okTProcessID && okCand)
        fHistFiles->Fill(1);

    fPrevInputFileName = currentFile;
    fTreeMismatch->Fill();
    fStatus = 0;

    delete esdTrackCuts;
    esdTrackCuts = nullptr;
    treeAOD = nullptr;
    file->Close();

    PostData(1, fOutput);
    PostData(2, fTreeMismatch);
}

//________________________________________________________________________
void AliAnalysisTaskCheckAODdAODMatching::MapAODtracks()
{

    memset(fAODMap, 0, sizeof(int) * 100000);
    for (int iTr = 0; iTr < fAOD->GetNumberOfTracks(); iTr++)
    {
        AliAODTrack *track = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(iTr));
        if (track->GetStatus() & AliESDtrack::kITSpureSA)
            continue;
        if (!(track->GetStatus() & AliESDtrack::kITSin))
            continue;
        double covtest[21];
        if (!track->GetCovarianceXYZPxPyPz(covtest))
            continue;
        int ind = (int)track->GetID();
        if (ind > -1 && ind < 100000)
            fAODMap[ind] = iTr;
    }
}
