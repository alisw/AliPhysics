/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSENonPromptLc
// \brief Analysis task to produce trees of Lc candidates for ML analyses of non-prompt Lc
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include <TRandom3.h>

#include "AliAODMCHeader.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliHFMLVarHandlerLctopKpi.h"
#include "AliHFMLVarHandlerNonPromptLc2V0bachelor.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisUtils.h"
#include "AliAODHandler.h"
#include "AliAODExtension.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisVertexingHF.h"

#include "AliAnalysisTaskSENonPromptLc.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSENonPromptLc);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSENonPromptLc::AliAnalysisTaskSENonPromptLc() : AliAnalysisTaskSE()
{
    /// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSENonPromptLc::AliAnalysisTaskSENonPromptLc(const char *name, int decayChannel, AliRDHFCuts *analysisCuts, bool createMLtree) :
    AliAnalysisTaskSE(name),
    fDecChannel(decayChannel),
    fCreateMLtree(createMLtree)
{
    /// Standard constructor
    SetAnalysisCuts(analysisCuts);

    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    DefineOutput(3, AliNormalizationCounter::Class());
    if (fCreateMLtree)
        DefineOutput(4, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskSENonPromptLc::~AliAnalysisTaskSENonPromptLc()
{
    // Destructor
    delete fOutput;
    delete fCounter;
    delete fListCuts;
    delete fRDCuts;

    if (fCreateMLtree && fMLhandler)
        delete fMLhandler; // it also deletes the TTree
}

//________________________________________________________________________
void AliAnalysisTaskSENonPromptLc::LocalInit()
{
    // Initialization

    if (fDecChannel == kLctopKpi)
    {
        AliRDHFCutsLctopKpi *copycut = new AliRDHFCutsLctopKpi(*(static_cast<AliRDHFCutsLctopKpi *>(fRDCuts)));
        PostData(2, copycut);
    }
    else if (fDecChannel == kLctopK0s || fDecChannel == kLctopiL)
    {
        AliRDHFCutsLctoV0 *copycut = new AliRDHFCutsLctoV0(*(static_cast<AliRDHFCutsLctoV0 *>(fRDCuts)));
        PostData(2, copycut);
    }

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSENonPromptLc::UserCreateOutputObjects()
{
    /// Create the output container
    //

    // Several histograms are more conveniently managed in a TList
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos");

    fHistNEvents = new TH1F("hNEvents", "number of events ", 17, -0.5, 16.5);
    fHistNEvents->GetXaxis()->SetBinLabel(1, "nEventsRead");
    fHistNEvents->GetXaxis()->SetBinLabel(2, "nEvents Matched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(3, "nEvents Mismatched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(4, "nEventsAnal");
    fHistNEvents->GetXaxis()->SetBinLabel(5, "n. passing IsEvSelected");
    fHistNEvents->GetXaxis()->SetBinLabel(6, "n. rejected due to trigger");
    fHistNEvents->GetXaxis()->SetBinLabel(7, "n. rejected due to not reco vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(8, "n. rejected for contr vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(9, "n. rejected for vertex out of accept");
    fHistNEvents->GetXaxis()->SetBinLabel(10, "n. rejected for pileup events");
    fHistNEvents->GetXaxis()->SetBinLabel(11, "no. of out centrality events");
    fHistNEvents->GetXaxis()->SetBinLabel(12, "no. of Lc candidates");
    fHistNEvents->GetXaxis()->SetBinLabel(13, "no. of Lc after filtering cuts");
    fHistNEvents->GetXaxis()->SetBinLabel(14, "no. of Lc after selection cuts");
    fHistNEvents->GetXaxis()->SetBinLabel(15, "no. of not on-the-fly rec Lc");
    fHistNEvents->GetXaxis()->SetBinLabel(16, "no. of Lc rejected by preselect");
    fHistNEvents->GetXaxis()->SetBinLabel(17, "no. of Lc rejected because of on-the-fly V0");
    fHistNEvents->GetXaxis()->SetNdivisions(1, false);
    fHistNEvents->SetMinimum(0);
    fOutput->Add(fHistNEvents);

    // Sparses for efficiencies (only gen)
    if(fReadMC)
        CreateEffSparses();

    //Counter for Normalization
    fCounter = new AliNormalizationCounter("NormalizationCounter");
    fCounter->Init();
    PostData(3, fCounter);

    //Create ML tree
    if (fCreateMLtree)
    {
        OpenFile(4);
        if (fDecChannel == kLctopKpi)
            fMLhandler = new AliHFMLVarHandlerLctopKpi(fPIDopt);
        else if (fDecChannel == kLctopK0s || fDecChannel == kLctopiL)
        {
            fMLhandler = new AliHFMLVarHandlerNonPromptLc2V0bachelor(fPIDopt);
            if(fUseKFRecoForV0bachelor)
                (dynamic_cast<AliHFMLVarHandlerNonPromptLc2V0bachelor *>(fMLhandler))->SetUseKFParticleReco(); // currently for V0bachelor only
        }

        fMLhandler->SetAddSingleTrackVars(fAddSingleTrackVar);
        if (fReadMC)
        {
            if (fFillOnlySignal)
                fMLhandler->SetFillOnlySignal();
            fMLhandler->SetFillBeautyMotherPt();
        }
        fMLtree = fMLhandler->BuildTree("treeMLLc", "treeMLLc");
        fMLtree->SetMaxVirtualSize(1.e+8);
        PostData(4, fMLtree);
    }

    //Set seed of gRandom
    if (fCreateMLtree && fEnableEvtSampling)
        gRandom->SetSeed(fSeedSampling);

    PostData(1, fOutput);

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSENonPromptLc::UserExec(Option_t * /*option*/)
{
    if (fCreateMLtree && fEnableEvtSampling && ((fOptionSampling == 0 && gRandom->Rndm() > fFracEvtToKeep) || (fOptionSampling == 1 && gRandom->Rndm() < 1-fFracEvtToKeep)))
    {
        PostData(1, fOutput);
        return;
    }

    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());

    fHistNEvents->Fill(0); // all events
    if (fAODProtection >= 0)
    {
        //   Protection against different number of events in the AOD and deltaAOD
        //   In case of discrepancy the event is rejected.
        int matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
        if (matchingAODdeltaAODlevel < 0 || (matchingAODdeltaAODlevel == 0 && fAODProtection == 1))
        {
            // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
            fHistNEvents->Fill(2);
            PostData(1, fOutput);
            return;
        }
        fHistNEvents->Fill(1);
    }

    TClonesArray *arrayCand = nullptr;
    if (!fAOD && AODEvent() && IsStandardAOD())
    {
        // In case there is an AOD handler writing a standard AOD, use the AOD
        // event in memory rather than the input (ESD) event.
        fAOD = dynamic_cast<AliAODEvent *>(AODEvent());
        // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
        // have to taken from the AOD event hold by the AliAODExtension
        AliAODHandler *aodHandler = dynamic_cast<AliAODHandler *>((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
        if (aodHandler->GetExtensions())
        {
            AliAODExtension *ext = dynamic_cast<AliAODExtension *>(aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root"));
            AliAODEvent *aodFromExt = ext->GetAOD();
            if (fDecChannel == kLctopKpi)
                arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Charm3Prong"));
            else if (fDecChannel == kLctopK0s || fDecChannel == kLctopiL)
                arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("CascadesHF"));
        }
    }
    else if (fAOD)
    {
        if (fDecChannel == kLctopKpi)
            arrayCand = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Charm3Prong"));
        else if (fDecChannel == kLctopK0s || fDecChannel == kLctopiL)
            arrayCand = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("CascadesHF"));
    }

    if (!fAOD || !arrayCand)
    {
        AliWarning("Candidate branch not found!\n");
        PostData(1, fOutput);
        return;
    }

    // fix for temporary bug in ESDfilter
    // the AODs with null vertex pointer didn't pass the PhysSel
    if (!fAOD->GetPrimaryVertex() || TMath::Abs(fAOD->GetMagneticField()) < 0.001)
    {
        PostData(1, fOutput);
        return;
    }

    fHistNEvents->Fill(3); // count event

    fCounter->StoreEvent(fAOD, fRDCuts, fReadMC);

    bool isEvSel = fRDCuts->IsEventSelected(fAOD);

    if (fRDCuts->IsEventRejectedDueToTrigger())
        fHistNEvents->Fill(5);
    if (fRDCuts->IsEventRejectedDueToNotRecoVertex())
        fHistNEvents->Fill(6);
    if (fRDCuts->IsEventRejectedDueToVertexContributors())
        fHistNEvents->Fill(7);
    if (fRDCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion())
        fHistNEvents->Fill(8);
    if (fRDCuts->IsEventRejectedDueToPileup())
        fHistNEvents->Fill(9);
    if (fRDCuts->IsEventRejectedDueToCentrality())
        fHistNEvents->Fill(10);

    TClonesArray *arrayMC = nullptr;
    AliAODMCHeader *mcHeader = nullptr;

    // load MC particles
    if (fReadMC)
    {
        arrayMC = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
        if (!arrayMC)
        {
            AliWarning("MC particles branch not found!");
            PostData(1, fOutput);
            return;
        }

        // load MC header
        mcHeader = dynamic_cast<AliAODMCHeader *>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
        if (!mcHeader)
        {
            AliWarning("MC header branch not found!");
            PostData(1, fOutput);
            return;
        }

        // fill MC acceptance histos
        FillMCGenAccHistos(arrayMC, mcHeader);
    }

    if (!isEvSel)
    {
        PostData(1, fOutput);
        return;
    }

    fHistNEvents->Fill(4); // accepted event

    // vHF object is needed to call the method that refills the missing info of the candidates
    // if they have been deleted in dAOD reconstruction phase
    // in order to reduce the size of the file
    AliAnalysisVertexingHF vHF = AliAnalysisVertexingHF();

    for (int iCand = 0; iCand < arrayCand->GetEntriesFast(); iCand++)
    {
        AliAODRecoDecayHF *lc = dynamic_cast<AliAODRecoDecayHF *>(arrayCand->UncheckedAt(iCand));

        bool unsetVtx = false;
        bool recVtx = false;
        AliAODVertex *origOwnVtx = nullptr;

        int isSelected = IsCandidateSelected(lc, &vHF, unsetVtx, recVtx, origOwnVtx);
        if (!isSelected || (fDecChannel == kLctopK0s && isSelected % 2 == 0) || (fDecChannel == kLctopiL && isSelected < 2))
        {
            if (unsetVtx)
                lc->UnsetOwnPrimaryVtx();
            if (recVtx)
                fRDCuts->CleanOwnPrimaryVtx(lc, fAOD, origOwnVtx);
            continue;
        }

        fHistNEvents->Fill(13); // candidate selected

        // get MC truth
        AliAODMCParticle *partLc = nullptr;
        int labLc = -1;
        int pdgCode0 = -999;
        int orig = 0;
        bool isCandInjected = false;
        float ptB = -999.;
        int pdgLctopKpi[3] = {2212, 321, 211};
        int pdgLctopK0s[2] = {2212, 310};
        int pdgLctopiL[2] = {211, 3122};
        int pdgDgK0stoDaughters[2] = {211, 211};
        int pdgDgLtoDaughters[2] = {2212, 211};

        if (fReadMC)
        {
            switch (fDecChannel)
            {
                case kLctopKpi:
                    labLc = (dynamic_cast<AliAODRecoDecayHF3Prong *>(lc))->MatchToMC(4122, arrayMC, 3, pdgLctopKpi);
                break;
                case kLctopK0s:
                    labLc = (dynamic_cast<AliAODRecoCascadeHF *>(lc))->MatchToMC(4122, 310, pdgLctopK0s, pdgDgK0stoDaughters, arrayMC, true);
                break;
                case kLctopiL:
                    labLc = (dynamic_cast<AliAODRecoCascadeHF *>(lc))->MatchToMC(4122, 3122, pdgLctopiL, pdgDgLtoDaughters, arrayMC, true);
                break;
            }

            if (labLc >= 0)
            {
                partLc = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labLc));
                if (fDecChannel == kLctopKpi) // check if signal is reflected
                {
                    int labDau0 = dynamic_cast<AliAODTrack *>(lc->GetDaughter(0))->GetLabel();
                    AliAODMCParticle *dau0 = dynamic_cast<AliAODMCParticle *>(arrayMC->UncheckedAt(TMath::Abs(labDau0)));
                    pdgCode0 = TMath::Abs(dau0->GetPdgCode());
                }
            }
            else
            {
                if (fKeepOnlyBkgFromHIJING)
                    isCandInjected = AliVertexingHFUtils::IsCandidateInjected(lc, mcHeader, arrayMC);
            }

            if (partLc)
            {
                orig = AliVertexingHFUtils::CheckOrigin(arrayMC, partLc, true);
                ptB = AliVertexingHFUtils::GetBeautyMotherPt(arrayMC, partLc);
            }
        }
        // fill tree for ML
        AliAODPidHF *pidHF = fRDCuts->GetPidHF();
        if (fCreateMLtree)
        {
            if (fDecChannel == kLctopKpi) // Lc->pKpi
            {
                if (fReadMC)
                {
                    int labD[3] = {-1, -1, -1};
                    //check if resonant decay
                    int  decay = 0;
                    if(labLc && partLc)
                        decay = AliVertexingHFUtils::CheckLcpKpiDecay(arrayMC, partLc, labD);
                    (dynamic_cast<AliHFMLVarHandlerLctopKpi *>(fMLhandler))->SetIsLcpKpiRes(decay);
                }
                if (isSelected == 1 || isSelected == 3) // pKpi
                {
                    bool isSignal = false;
                    bool isBkg = false;
                    bool isPrompt = false;
                    bool isFD = false;
                    bool isRefl = false;

                    if (fReadMC)
                    {
                        if (labLc >= 0)
                        {
                            if (pdgCode0 == 211)
                                isRefl = true;
                            if (orig == 4)
                                isPrompt = true;
                            else if (orig == 5)
                                isFD = true;

                            if (orig >= 4)
                                isSignal = true;
                        }
                        else
                        {
                            if (!isCandInjected)
                                isBkg = true;
                        }
                        fMLhandler->SetBeautyMotherPt(ptB);
                    }

                    fMLhandler->SetCandidateType(isSignal, isBkg, isPrompt, isFD, isRefl);
                    bool okSetVar = fMLhandler->SetVariables(lc, fAOD->GetMagneticField(), AliHFMLVarHandlerLctopKpi::kpKpi, pidHF);
                    if (okSetVar && !(fReadMC && !isSignal && !isBkg && !isPrompt && !isFD && !isRefl))
                        fMLhandler->FillTree();
                }
                if (isSelected >= 2) // piKp
                {
                    bool isSignal = false;
                    bool isBkg = false;
                    bool isPrompt = false;
                    bool isFD = false;
                    bool isRefl = false;

                    if (fReadMC)
                    {
                        if (labLc >= 0)
                        {
                            if (pdgCode0 == 2212)
                                isRefl = true;
                            if (orig == 4)
                                isPrompt = true;
                            else if (orig == 5)
                                isFD = true;

                            if (orig >= 4)
                                isSignal = true;
                        }
                        else
                        {
                            if (!isCandInjected)
                                isBkg = true;
                        }
                        fMLhandler->SetBeautyMotherPt(ptB);
                    }

                    fMLhandler->SetCandidateType(isSignal, isBkg, isPrompt, isFD, isRefl);
                    bool okSetVar = fMLhandler->SetVariables(lc, fAOD->GetMagneticField(), AliHFMLVarHandlerLctopKpi::kpiKp, pidHF);
                    if (okSetVar && !(fReadMC && !isSignal && !isBkg && !isPrompt && !isFD && !isRefl)) // add tag in tree handler for signal from pileup events?
                        fMLhandler->FillTree();
                }
            }
            else // Lc->V0bachelor
            {
                bool isSignal = false;
                bool isBkg = false;
                bool isPrompt = false;
                bool isFD = false;
                bool isRefl = false;

                if (fReadMC)
                {
                    if (labLc >= 0)
                    {
                        if (orig == 4)
                            isPrompt = true;
                        else if (orig == 5)
                            isFD = true;

                        if (orig >= 4)
                            isSignal = true;
                    }
                    else
                    {
                        if (!isCandInjected)
                            isBkg = true;
                    }
                    fMLhandler->SetBeautyMotherPt(ptB);
                }

                int channel = (fDecChannel == kLctopK0s) ? AliHFMLVarHandlerNonPromptLc2V0bachelor::kpK0s : AliHFMLVarHandlerNonPromptLc2V0bachelor::kpiL;

                fMLhandler->SetCandidateType(isSignal, isBkg, isPrompt, isFD, isRefl);
                bool okSetVar = fMLhandler->SetVariables(lc, fAOD->GetMagneticField(), channel, pidHF);
                if (okSetVar && !(fReadMC && !isSignal && !isBkg && !isPrompt && !isFD && !isRefl)) // add tag in tree handler for signal from pileup events?
                    fMLhandler->FillTree();
            }
        }

        if (unsetVtx)
            lc->UnsetOwnPrimaryVtx();
        if (recVtx)
            fRDCuts->CleanOwnPrimaryVtx(lc, fAOD, origOwnVtx);
    }

    PostData(1, fOutput);
    PostData(3, fCounter);
    if(fCreateMLtree)
        PostData(4, fMLtree);
}

//________________________________________________________________________
int AliAnalysisTaskSENonPromptLc::IsCandidateSelected(AliAODRecoDecayHF *&lc, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx)
{

    if (!lc || !vHF)
        return 0;
    fHistNEvents->Fill(11);

    // Preselection to speed up task
    TObjArray arrDauTracks(3);
    int nDau = 3;
    if (fDecChannel == kLctopK0s)
        nDau = 2;

    for (int iDau = 0; iDau < nDau; iDau++)
    {
        AliAODTrack *track = vHF->GetProng(fAOD, lc, iDau);
        arrDauTracks.AddAt(track, iDau);
    }

    if (!fRDCuts->PreSelect(arrDauTracks))
    {
        fHistNEvents->Fill(15);
        return 0;
    }

    bool isSelBit = true;
    if (fDecChannel == kLctopKpi)
    {
        isSelBit = lc->HasSelectionBit(AliRDHFCuts::kLcCuts);
        if (!isSelBit || !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF3Prong *>(lc)))
        {
            fHistNEvents->Fill(14);
            return 0;
        }
    }
    else
    {
        isSelBit = lc->HasSelectionBit(AliRDHFCuts::kLctoV0Cuts);
        if (!isSelBit || !vHF->FillRecoCasc(fAOD, dynamic_cast<AliAODRecoCascadeHF *>(lc), false))
        {
            fHistNEvents->Fill(14);
            return 0;
        }

        AliAODv0 *v0part = dynamic_cast<AliAODv0*>(dynamic_cast<AliAODRecoCascadeHF *>(lc)->Getv0());
        if (v0part->GetOnFlyStatus()) {
            fHistNEvents->Fill(16);
            return 0;
        }
    }

    fHistNEvents->Fill(12);

    unsetVtx = false;
    if (!lc->GetOwnPrimaryVtx())
    {
        lc->SetOwnPrimaryVtx(dynamic_cast<AliAODVertex *>(fAOD->GetPrimaryVertex()));
        unsetVtx = true;
        // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
        // Pay attention if you use continue inside this loop!!!
    }

    double ptLc = lc->Pt();
    double yLc = lc->Y(4122);

    if (fCreateMLtree && fEnableCandSampling) // apply sampling in pt
    {
        double pseudoRand = ptLc * 1000. - (long)(ptLc * 1000);
        if (pseudoRand > fFracCandToKeep && ptLc < fMaxCandPtSampling)
        {
            if (unsetVtx)
                lc->UnsetOwnPrimaryVtx();
            return 0;
        }
    }

    int ptBin = fRDCuts->PtBin(ptLc);
    if (ptBin < 0)
    {
        if (unsetVtx)
            lc->UnsetOwnPrimaryVtx();
        return 0;
    }

    bool isFidAcc = fRDCuts->IsInFiducialAcceptance(ptLc, yLc);
    if (!isFidAcc)
    {
        if (unsetVtx)
            lc->UnsetOwnPrimaryVtx();
        return 0;
    }

    int isSelected = fRDCuts->IsSelected(lc, AliRDHFCuts::kAll, fAOD);
    if (!isSelected)
    {
        if (unsetVtx)
            lc->UnsetOwnPrimaryVtx();
        return 0;
    }

    recVtx = false;
    origOwnVtx = nullptr;

    if (fRDCuts->GetIsPrimaryWithoutDaughters())
    {
        if (lc->GetOwnPrimaryVtx())
            origOwnVtx = new AliAODVertex(*lc->GetOwnPrimaryVtx());
        if (fRDCuts->RecalcOwnPrimaryVtx(lc, fAOD))
            recVtx = true;
        else
            fRDCuts->CleanOwnPrimaryVtx(lc, fAOD, origOwnVtx);
    }

    // retvalue case for kLctopKpi
    // 1  Lc->pKpi
    // 2  Lc->piKp
    // 3  Lc->pKpi AND Lc->piKp

    // retvalue case for kLctopK0s and kLctopiL
    // 1  Lc->K0S + p
    // 2  Lc->LambdaBar + pi
    // 3  Lc->K0S + p AND Lc->LambdaBar + pi
    // 4  Lc->Lambda + pi
    // 5  Lc->K0S + p AND Lc->Lambda + pi
    // 6  Lc->LambdaBar + pi AND Lc->Lambda + pi
    // 7  Lc->K0S + p AND Lc->LambdaBar + pi AND Lc->Lambda + pi

    return isSelected;
}

//________________________________________________________________________
void AliAnalysisTaskSENonPromptLc::FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader)
{
    /// Fill MC histos for cuts study
    ///    - at GenLimAccStep and AccStep (if fFillAcceptanceLevel=false)
    ///    - at AccStep (if fFillAcceptanceLevel=true)

    double zMCVertex = mcHeader->GetVtxZ(); //vertex MC
    if (TMath::Abs(zMCVertex) <= fRDCuts->GetMaxVtxZ())
    {
        for (int iPart = 0; iPart < arrayMC->GetEntriesFast(); iPart++)
        {
            AliAODMCParticle *mcPart = dynamic_cast<AliAODMCParticle *>(arrayMC->At(iPart));

            if (TMath::Abs(mcPart->GetPdgCode()) == 4122)
            {
                int orig = AliVertexingHFUtils::CheckOrigin(arrayMC, mcPart, true); //Prompt = 4, FeedDown = 5
                bool isParticleFromOutOfBunchPileUpEvent = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iPart, mcHeader, arrayMC);

                int deca = 0;
                bool isGoodDecay = false;
                int labDau[3] = {-1, -1, -1};
                bool isFidAcc = false;
                bool isDaugInAcc = false;

                if(fDecChannel == kLctopKpi)
                    deca = AliVertexingHFUtils::CheckLcpKpiDecay(arrayMC, mcPart, labDau);
                else
                    deca = AliVertexingHFUtils::CheckLcV0bachelorDecay(arrayMC, mcPart, labDau);

                if (labDau[0] == -1)
                    continue; //protection against unfilled array of labels
                if ((fDecChannel == kLctopKpi && deca > 0) || (fDecChannel == kLctopK0s && deca == 1) || (fDecChannel == kLctopiL && deca == 2))
                    isGoodDecay = true;

                if (isGoodDecay)
                {
                    double pt = mcPart->Pt();
                    double rapid = mcPart->Y();
                    isFidAcc = fRDCuts->IsInFiducialAcceptance(pt, rapid);
                    isDaugInAcc = CheckDaugAcc(arrayMC, 3, labDau);

                    if ((fFillAcceptanceLevel && isFidAcc && isDaugInAcc) || (!fFillAcceptanceLevel && TMath::Abs(rapid) < 0.5))
                    {
                        if (orig == 4 && !isParticleFromOutOfBunchPileUpEvent)
                        {
                            double var4nSparseAcc[knVarForSparseAcc] = {pt, rapid,(double)deca};
                            fnSparseMC[0]->Fill(var4nSparseAcc);
                        }
                        else if (orig == 5 && !isParticleFromOutOfBunchPileUpEvent)
                        {
                            double ptB = AliVertexingHFUtils::GetBeautyMotherPt(arrayMC, mcPart);
                            double var4nSparseAcc[knVarForSparseAccFD] = {pt, rapid, ptB,(double)deca};
                            fnSparseMC[1]->Fill(var4nSparseAcc);
                        }
                    }
                }
            }
        }
    }
}

//________________________________________________________________________
bool AliAnalysisTaskSENonPromptLc::CheckDaugAcc(TClonesArray *arrayMC, int nProng, int *labDau)
{
    /// check if the decay products are in the good eta and pt range

    for (int iProng = 0; iProng < nProng; iProng++)
    {
        AliAODMCParticle *mcPartDaughter = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labDau[iProng]));
        if (!mcPartDaughter)
            return false;

        double eta = mcPartDaughter->Eta();
        double pt = mcPartDaughter->Pt();

        if (TMath::Abs(eta) > 0.9 || pt < 0.1)
            return false;
    }
    return true;
}

//_________________________________________________________________________
void AliAnalysisTaskSENonPromptLc::CreateEffSparses()
{
    /// use sparses to be able to add variables if needed (multiplicity, Zvtx, etc)

    int nPtBinsCutObj = fRDCuts->GetNPtBins();
    float *ptLims = fRDCuts->GetPtBinLimits();
    int nPtBins = (int)ptLims[nPtBinsCutObj];
    if (fUseFinPtBinsForSparse)
        nPtBins = nPtBins * 10;

    int nBinsAccFD[knVarForSparseAccFD] = {nPtBins, 20, nPtBins, 6};
    double xminAccFD[knVarForSparseAccFD] = {0., -1., 0., 0};
    double xmaxAccFD[knVarForSparseAccFD] = {ptLims[nPtBinsCutObj], 1., ptLims[nPtBinsCutObj], 6};

    int nBinsAcc[knVarForSparseAcc] = {nPtBins, 20, 6};
    double xminAcc[knVarForSparseAcc] = {0., -1., 0};
    double xmaxAcc[knVarForSparseAcc] = {ptLims[nPtBinsCutObj], 1., 6};

    TString label[2] = {"fromC", "fromB"};
    for (int iHist = 0; iHist < 2; iHist++)
    {
        TString titleSparse = Form("MC nSparse (%s)- %s", fFillAcceptanceLevel ? "Acc.Step" : "Gen.Acc.Step", label[iHist].Data());
        fnSparseMC[iHist] = new THnSparseF(Form("fnSparseAcc_%s", label[iHist].Data()), titleSparse.Data(),
                                           (iHist == 0) ? knVarForSparseAcc : knVarForSparseAccFD,
                                           (iHist == 0) ? nBinsAcc : nBinsAccFD,
                                           (iHist == 0) ? xminAcc : xminAccFD,
                                           (iHist == 0) ? xmaxAcc : xmaxAccFD);
        fnSparseMC[iHist]->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/c)");
        fnSparseMC[iHist]->GetAxis(1)->SetTitle("#it{y}");
        if (iHist == 0)
        {
            fnSparseMC[iHist]->GetAxis(2)->SetTitle("resonant channel");
        }
        else
        {
            fnSparseMC[iHist]->GetAxis(2)->SetTitle("#it{p}_{T}^{B} (GeV/c)");
            fnSparseMC[iHist]->GetAxis(3)->SetTitle("resonant channel");
        }
        fOutput->Add(fnSparseMC[iHist]);
    }
}
