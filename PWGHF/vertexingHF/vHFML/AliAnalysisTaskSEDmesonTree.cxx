/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSEDmesonTree
// \brief Analysis task to produce trees of Lc candidates for ML analyses of non-prompt Lc
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include <TRandom3.h>

#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliHFMLVarHandlerD0toKpi.h"
#include "AliHFMLVarHandlerDplustoKpipi.h"
#include "AliHFMLVarHandlerDstartoD0pi.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisUtils.h"
#include "AliAODHandler.h"
#include "AliAODExtension.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskSEDmesonTree.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEDmesonTree);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEDmesonTree::AliAnalysisTaskSEDmesonTree() : AliAnalysisTaskSE()
{
    /// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEDmesonTree::AliAnalysisTaskSEDmesonTree(const char *name, int decayChannel, AliRDHFCuts *analysisCuts, bool createMLtree) :
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
AliAnalysisTaskSEDmesonTree::~AliAnalysisTaskSEDmesonTree()
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
void AliAnalysisTaskSEDmesonTree::LocalInit()
{
    // Initialization

    switch (fDecChannel)
    {
        case kD0toKpi:
        {
            AliRDHFCutsD0toKpi *copycut = new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi *>(fRDCuts)));
            PostData(2, copycut);
            break;
        }
        case kDplustoKpipi:
        {
            AliRDHFCutsDplustoKpipi *copycut = new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi *>(fRDCuts)));
            PostData(2, copycut);
            break;
        }
        case kDstartoD0pi:
        {
            AliRDHFCutsDStartoKpipi *copycut = new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi *>(fRDCuts)));
            PostData(2, copycut);
            break;
        }
    }

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDmesonTree::UserCreateOutputObjects()
{
    /// Create the output container
    //

    // Several histograms are more conveniently managed in a TList
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos");

    fHistNEvents = new TH1F("hNEvents", "number of events ", 16, -0.5, 15.5);
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
    fHistNEvents->GetXaxis()->SetBinLabel(12, "no. of D candidates");
    fHistNEvents->GetXaxis()->SetBinLabel(13, "no. of D after filtering cuts");
    fHistNEvents->GetXaxis()->SetBinLabel(14, "no. of D after selection cuts");
    fHistNEvents->GetXaxis()->SetBinLabel(15, "no. of not on-the-fly rec D");
    fHistNEvents->GetXaxis()->SetBinLabel(16, "no. of D rejected by preselect");
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
        switch (fDecChannel)
        {
            case kD0toKpi:
                fMLhandler = new AliHFMLVarHandlerD0toKpi(fPIDopt);
                break;
            case kDplustoKpipi:
                fMLhandler = new AliHFMLVarHandlerDplustoKpipi(fPIDopt);
                break;
            case kDstartoD0pi:
                fMLhandler = new AliHFMLVarHandlerDstartoD0pi(fPIDopt);
                break;
        }

        fMLhandler->SetAddSingleTrackVars(fAddSingleTrackVar);
        fMLhandler->SetAddGlobalEventVariables(fAddNtrkl, fAddCentr, fCentEstimator);
        if (fReadMC)
        {
            if (fFillOnlySignal)
                fMLhandler->SetFillOnlySignal();
            fMLhandler->SetFillBeautyMotherPt();
        }

        fMLtree = fMLhandler->BuildTree("treeMLD", "treeMLD");
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
void AliAnalysisTaskSEDmesonTree::UserExec(Option_t * /*option*/)
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
            switch (fDecChannel)
            {
                case kD0toKpi:
                    arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Charm2Prong"));
                    break;
                case kDplustoKpipi:
                    arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Charm3Prong"));
                    break;
                case kDstartoD0pi:
                    arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("CascadesHF"));
                break;
            }
        }
    }
    else if (fAOD)
    {
        switch (fDecChannel)
        {
            case kD0toKpi:
                arrayCand = dynamic_cast<TClonesArray *>(fAOD->FindObject("Charm2Prong"));
                break;
            case kDplustoKpipi:
                arrayCand = dynamic_cast<TClonesArray *>(fAOD->FindObject("Charm3Prong"));
                break;
            case kDstartoD0pi:
                arrayCand = dynamic_cast<TClonesArray *>(fAOD->FindObject("CascadesHF"));
            break;
        }
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
        AliAODRecoDecayHF *dMeson = dynamic_cast<AliAODRecoDecayHF *>(arrayCand->UncheckedAt(iCand));

        bool unsetVtx = false;
        bool recVtx = false;
        AliAODVertex *origOwnVtx = nullptr;

        int isSelected = IsCandidateSelected(dMeson, &vHF, unsetVtx, recVtx, origOwnVtx);
        if (!isSelected)
        {
            if (unsetVtx)
                dMeson->UnsetOwnPrimaryVtx();
            if (recVtx)
                fRDCuts->CleanOwnPrimaryVtx(dMeson, fAOD, origOwnVtx);
            continue;
        }

        fHistNEvents->Fill(13); // candidate selected

        // get MC truth
        AliAODMCParticle *partD = nullptr;
        int labD = -1;
        int pdgCode0 = -999;
        int orig = 0;
        bool isCandInjected = false;
        float ptB = -999.;
        int pdgD0Dau[2] = {321, 211};
        int pdgDplusDau[3] = {321, 211, 211};
        int pdgDstarDau[2] = {421, 211};

        if (fReadMC)
        {
            switch (fDecChannel)
            {
                case kD0toKpi:
                    labD = (dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson))->MatchToMC(fPdgD, arrayMC, 2, pdgD0Dau);
                break;
                case kDplustoKpipi:
                    labD = (dynamic_cast<AliAODRecoDecayHF3Prong *>(dMeson))->MatchToMC(fPdgD, arrayMC, 3, pdgDplusDau);
                break;
                case kDstartoD0pi:
                    labD = (dynamic_cast<AliAODRecoCascadeHF *>(dMeson))->MatchToMC(fPdgD, 421, pdgDstarDau, pdgD0Dau, arrayMC, true);
                break;
            }

            if (labD >= 0)
            {
                partD = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labD));
                if (fDecChannel == kD0toKpi) // check if signal is reflected
                {
                    int labDau0 = dynamic_cast<AliAODTrack *>(dMeson->GetDaughter(0))->GetLabel();
                    AliAODMCParticle *dau0 = dynamic_cast<AliAODMCParticle *>(arrayMC->UncheckedAt(TMath::Abs(labDau0)));
                    pdgCode0 = TMath::Abs(dau0->GetPdgCode());
                }
            }
            else
            {
                if (fKeepOnlyBkgFromHIJING)
                    isCandInjected = AliVertexingHFUtils::IsCandidateInjected(dMeson, mcHeader, arrayMC);
            }

            if (partD)
            {
                orig = AliVertexingHFUtils::CheckOrigin(arrayMC, partD, true);
                ptB = AliVertexingHFUtils::GetBeautyMotherPt(arrayMC, partD);
            }
        }
        // fill tree for ML
        AliAODPidHF *pidHF = fRDCuts->GetPidHF();
        if (fCreateMLtree)
        {
            fMLhandler->SetGlobalEventVariables(fAOD);
            if (fDecChannel == kD0toKpi)
            {
                if (isSelected == 1 || isSelected == 3) // D0
                {
                    bool isSignal = false;
                    bool isBkg = false;
                    bool isPrompt = false;
                    bool isFD = false;
                    bool isRefl = false;

                    if (fReadMC)
                    {
                        if (labD >= 0)
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
                    bool okSetVar = fMLhandler->SetVariables(dMeson, fAOD->GetMagneticField(), AliHFMLVarHandlerD0toKpi::kD0, pidHF);
                    if (okSetVar && !(fReadMC && !isSignal && !isBkg && !isPrompt && !isFD && !isRefl))
                        fMLhandler->FillTree();
                }
                if (isSelected >= 2) // D0bar
                {
                    bool isSignal = false;
                    bool isBkg = false;
                    bool isPrompt = false;
                    bool isFD = false;
                    bool isRefl = false;

                    if (fReadMC)
                    {
                        if (labD >= 0)
                        {
                            if (pdgCode0 == 321)
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
                    bool okSetVar = fMLhandler->SetVariables(dMeson, fAOD->GetMagneticField(), AliHFMLVarHandlerD0toKpi::kD0bar, pidHF);
                    if (okSetVar && !(fReadMC && !isSignal && !isBkg && !isPrompt && !isFD && !isRefl)) // add tag in tree handler for signal from pileup events?
                        fMLhandler->FillTree();
                }
                break;
            }
            else
            {
                bool isSignal = false;
                bool isBkg = false;
                bool isPrompt = false;
                bool isFD = false;
                bool isRefl = false;

                if (fReadMC)
                {
                    if (labD >= 0)
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

                fMLhandler->SetCandidateType(isSignal, isBkg, isPrompt, isFD, isRefl);
                bool okSetVar = fMLhandler->SetVariables(dMeson, fAOD->GetMagneticField(), 0, pidHF);
                if (okSetVar && !(fReadMC && !isSignal && !isBkg && !isPrompt && !isFD && !isRefl)) // add tag in tree handler for signal from pileup events?
                    fMLhandler->FillTree();
                break;
            }
        }

        if (unsetVtx)
            dMeson->UnsetOwnPrimaryVtx();
        if (recVtx)
            fRDCuts->CleanOwnPrimaryVtx(dMeson, fAOD, origOwnVtx);
    }

    PostData(1, fOutput);
    PostData(3, fCounter);
    if(fCreateMLtree)
        PostData(4, fMLtree);
}

//________________________________________________________________________
int AliAnalysisTaskSEDmesonTree::IsCandidateSelected(AliAODRecoDecayHF *&dMeson, AliAnalysisVertexingHF *vHF, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx)
{

    if (!dMeson || !vHF)
        return 0;
    fHistNEvents->Fill(11);

    // Preselection to speed up task
    TObjArray arrDauTracks(3);
    int nDau = 3;
    if (fDecChannel == kD0toKpi)
        nDau = 2;

    for (int iDau = 0; iDau < nDau; iDau++)
    {
        AliAODTrack *track = vHF->GetProng(fAOD, dMeson, iDau);
        arrDauTracks.AddAt(track, iDau);
    }

    if (!fRDCuts->PreSelect(arrDauTracks))
    {
        fHistNEvents->Fill(15);
        return 0;
    }

    bool isSelBit = true;
    switch (fDecChannel)
    {
        case kD0toKpi:
            isSelBit = dMeson->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts);
            if (!isSelBit || !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF2Prong *>(dMeson)))
            {
                fHistNEvents->Fill(14);
                return 0;
            }
            break;
        case kDplustoKpipi:
            isSelBit = dMeson->HasSelectionBit(AliRDHFCuts::kDplusCuts);
            if (!isSelBit || !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF3Prong *>(dMeson)))
            {
                fHistNEvents->Fill(14);
                return 0;
            }
            break;
        case kDstartoD0pi:
            isSelBit = dMeson->HasSelectionBit(AliRDHFCuts::kDstarCuts);
            if (!isSelBit || !vHF->FillRecoCasc(fAOD, dynamic_cast<AliAODRecoCascadeHF *>(dMeson), false))
            {
                fHistNEvents->Fill(14);
                return 0;
            }
            break;
    }

    fHistNEvents->Fill(12);

    unsetVtx = false;
    if (!dMeson->GetOwnPrimaryVtx())
    {
        dMeson->SetOwnPrimaryVtx(dynamic_cast<AliAODVertex *>(fAOD->GetPrimaryVertex()));
        unsetVtx = true;
        // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
        // Pay attention if you use continue inside this loop!!!
    }

    double ptD = dMeson->Pt();
    double yD = dMeson->Y(fPdgD);

    if (fCreateMLtree && fEnableCandSampling) // apply sampling in pt
    {
        double pseudoRand = ptD * 1000. - (long)(ptD * 1000);
        if (pseudoRand > fFracCandToKeep && ptD < fMaxCandPtSampling)
        {
            if (unsetVtx)
                dMeson->UnsetOwnPrimaryVtx();
            return 0;
        }
    }

    int ptBin = fRDCuts->PtBin(ptD);
    if (ptBin < 0)
    {
        if (unsetVtx)
            dMeson->UnsetOwnPrimaryVtx();
        return 0;
    }

    bool isFidAcc = fRDCuts->IsInFiducialAcceptance(ptD, yD);
    if (!isFidAcc)
    {
        if (unsetVtx)
            dMeson->UnsetOwnPrimaryVtx();
        return 0;
    }

    int isSelected = fRDCuts->IsSelected(dMeson, AliRDHFCuts::kAll, fAOD);
    if (!isSelected)
    {
        if (unsetVtx)
            dMeson->UnsetOwnPrimaryVtx();
        return 0;
    }

    recVtx = false;
    origOwnVtx = nullptr;

    if (fRDCuts->GetIsPrimaryWithoutDaughters())
    {
        if (dMeson->GetOwnPrimaryVtx())
            origOwnVtx = new AliAODVertex(*dMeson->GetOwnPrimaryVtx());
        if (fRDCuts->RecalcOwnPrimaryVtx(dMeson, fAOD))
            recVtx = true;
        else
            fRDCuts->CleanOwnPrimaryVtx(dMeson, fAOD, origOwnVtx);
    }

    return isSelected;
}

//________________________________________________________________________
void AliAnalysisTaskSEDmesonTree::FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader)
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

            if (TMath::Abs(mcPart->GetPdgCode()) == fPdgD)
            {
                int orig = AliVertexingHFUtils::CheckOrigin(arrayMC, mcPart, true); //Prompt = 4, FeedDown = 5
                bool isParticleFromOutOfBunchPileUpEvent = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iPart, mcHeader, arrayMC);

                int deca = 0;
                bool isGoodDecay = false;
                int labDau[3] = {-1, -1, -1};
                bool isFidAcc = false;
                bool isDaugInAcc = false;
                int nDau = -1;

                switch(fDecChannel)
                {
                    case kD0toKpi:
                        nDau = 2;
                        deca = AliVertexingHFUtils::CheckD0Decay(arrayMC, mcPart, labDau);
                        break;
                    case kDplustoKpipi:
                        nDau = 3;
                        deca = AliVertexingHFUtils::CheckDplusDecay(arrayMC, mcPart, labDau);
                        break;
                    case kDstartoD0pi:
                        nDau = 3;
                        deca = AliVertexingHFUtils::CheckDstarDecay(arrayMC, mcPart, labDau);
                        break;
                }

                if (labDau[0] == -1)
                    continue; //protection against unfilled array of labels
                if (deca > 0)
                    isGoodDecay = true;

                if (isGoodDecay)
                {
                    double pt = mcPart->Pt();
                    double rapid = mcPart->Y();
                    isFidAcc = fRDCuts->IsInFiducialAcceptance(pt, rapid);
                    isDaugInAcc = CheckDaugAcc(arrayMC, nDau, labDau);

                    if ((fFillAcceptanceLevel && isFidAcc && isDaugInAcc) || (!fFillAcceptanceLevel && TMath::Abs(rapid) < 0.5))
                    {
                        if (orig == 4 && !isParticleFromOutOfBunchPileUpEvent)
                        {
                            double var4nSparseAcc[knVarForSparseAcc] = {pt, rapid};
                            fnSparseMC[0]->Fill(var4nSparseAcc);
                        }
                        else if (orig == 5 && !isParticleFromOutOfBunchPileUpEvent)
                        {
                            double ptB = AliVertexingHFUtils::GetBeautyMotherPt(arrayMC, mcPart);
                            double var4nSparseAcc[knVarForSparseAccFD] = {pt, rapid, ptB};
                            fnSparseMC[1]->Fill(var4nSparseAcc);
                        }
                    }
                }
            }
        }
    }
}

//________________________________________________________________________
bool AliAnalysisTaskSEDmesonTree::CheckDaugAcc(TClonesArray *arrayMC, int nProng, int *labDau)
{
    /// check if the decay products are in the good eta and pt range

    for (int iProng = 0; iProng < nProng; iProng++)
    {
        bool isSoftPion = false;
        AliAODMCParticle *mcPartDaughter = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labDau[iProng]));
        if (!mcPartDaughter)
            return false;

        if(fDecChannel == kDstartoD0pi)
        {
            AliAODMCParticle *mother = dynamic_cast<AliAODMCParticle *>(arrayMC->At(mcPartDaughter->GetMother()));
            if(TMath::Abs(mother->GetPdgCode()) == 413)
                isSoftPion = true;
        }

        double eta = mcPartDaughter->Eta();
        double pt = mcPartDaughter->Pt();
        double minPt = (!isSoftPion) ? 0.1 : 0.06;

        if (TMath::Abs(eta) > 0.9 || pt < minPt)
            return false;
    }
    return true;
}

//_________________________________________________________________________
void AliAnalysisTaskSEDmesonTree::CreateEffSparses()
{
    /// use sparses to be able to add variables if needed (multiplicity, Zvtx, etc)

    int nPtBinsCutObj = fRDCuts->GetNPtBins();
    float *ptLims = fRDCuts->GetPtBinLimits();
    int nPtBins = (int)ptLims[nPtBinsCutObj];
    if (fUseFinPtBinsForSparse)
        nPtBins = nPtBins * 10;

    int nBinsAcc[knVarForSparseAccFD] = {nPtBins, 20, nPtBins};
    double xminAcc[knVarForSparseAccFD] = {0., -1., 0.};
    double xmaxAcc[knVarForSparseAccFD] = {ptLims[nPtBinsCutObj], 1., ptLims[nPtBinsCutObj]};

    TString label[2] = {"fromC", "fromB"};
    for (int iHist = 0; iHist < 2; iHist++)
    {
        TString titleSparse = Form("MC nSparse (%s)- %s", fFillAcceptanceLevel ? "Acc.Step" : "Gen.Acc.Step", label[iHist].Data());
        fnSparseMC[iHist] = new THnSparseF(Form("fnSparseAcc_%s", label[iHist].Data()), titleSparse.Data(), (iHist == 0) ? knVarForSparseAcc : knVarForSparseAccFD, nBinsAcc, xminAcc, xmaxAcc);
        fnSparseMC[iHist]->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/c)");
        fnSparseMC[iHist]->GetAxis(1)->SetTitle("#it{y}");
        if (iHist == 1)
            fnSparseMC[iHist]->GetAxis(2)->SetTitle("#it{p}_{T}^{B} (GeV/c)");
        fOutput->Add(fnSparseMC[iHist]);
    }
}
