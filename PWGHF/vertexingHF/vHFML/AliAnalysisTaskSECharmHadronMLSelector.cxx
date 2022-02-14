/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSECharmHadronMLSelector
// \brief Analysis task to select charm-hadron candidates based on ML response for subsequent tasks
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch

#include <TDatabasePDG.h>

#include "AliAnalysisTaskSECharmHadronMLSelector.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliHFMLResponseDplustoKpipi.h"
#include "AliHFMLResponseDstoKKpi.h"
#include "AliHFMLResponseD0toKpi.h"
#include "AliHFMLResponseDstartoD0pi.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODHandler.h"
#include "AliAnalysisVertexingHF.h"

ClassImp(AliAnalysisTaskSECharmHadronMLSelector)

    //____________________________________________________________________________________________________
    AliAnalysisTaskSECharmHadronMLSelector::AliAnalysisTaskSECharmHadronMLSelector() : AliAnalysisTaskSE()
{
    /// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSECharmHadronMLSelector::AliAnalysisTaskSECharmHadronMLSelector(const char *name, int decayChannel, AliRDHFCuts *analysisCuts) : AliAnalysisTaskSE(name),
                                                                                                                                                fDecChannel(decayChannel)
{
    /// Standard constructor
    SetAnalysisCuts(analysisCuts);

    DefineOutput(1, TList::Class());
    DefineOutput(2, AliRDHFCuts::Class());
}

//________________________________________________________________________
AliAnalysisTaskSECharmHadronMLSelector::~AliAnalysisTaskSECharmHadronMLSelector()
{
    // Destructor
    delete fOutput;
    delete fListCuts;
    delete fRDCuts;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronMLSelector::LocalInit()
{
    // Initialization

    if (fDecChannel == kDplustoKpipi)
    {
        AliRDHFCutsDplustoKpipi *copycut = new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi *>(fRDCuts)));
        PostData(2, copycut);
    }
    else if (fDecChannel == kDstoKKpi)
    {
        AliRDHFCutsDstoKKpi *copycut = new AliRDHFCutsDstoKKpi(*(static_cast<AliRDHFCutsDstoKKpi *>(fRDCuts)));
        PostData(2, copycut);
    }
    else if (fDecChannel == kD0toKpi)
    {
        AliRDHFCutsD0toKpi *copycut = new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi *>(fRDCuts)));
        PostData(2, copycut);
    }
    else if (fDecChannel == kDstartoD0pi)
    {
        AliRDHFCutsDStartoKpipi *copycut = new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi *>(fRDCuts)));
        PostData(2, copycut);
    }
    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronMLSelector::UserCreateOutputObjects()
{
    /// Create the output container
    //

    // Several histograms are more conveniently managed in a TList
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos");

    fHistNEvents = new TH1F("hNEvents", "number of events ", 6, 0.5, 6.5);
    fHistNEvents->GetXaxis()->SetBinLabel(1, "nEventsRead");
    fHistNEvents->GetXaxis()->SetBinLabel(2, "nEvents Matched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(3, "nEvents Mismatched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(4, "nEventsAnal");
    fHistNEvents->GetXaxis()->SetBinLabel(5, "Rej. Phys. Sel. + Trig");
    fHistNEvents->GetXaxis()->SetBinLabel(6, "No vertex");
    fHistNEvents->GetXaxis()->SetNdivisions(1, false);
    fHistNEvents->SetMinimum(0);
    fOutput->Add(fHistNEvents);

    fHistNselCand = new TH1F("hNselCand", "number of selected candidates per event; number of candidates; entries", 1000, 0., 1000.);
    fHistNallCand = new TH1F("hNallCand", "number of all candidates per event; number of candidates; entries", 1000, 0., 1000.);
    fOutput->Add(fHistNselCand);
    fOutput->Add(fHistNallCand);

    for(int iHist = 0; iHist < 3; iHist++)
    {
        fHistBDTOutputVsPt[iHist] = new TH2F(Form("fHistBDTOutput%dVsPt", iHist),
                                             Form(";#it{p}_{T} (GeV/#it{c});BDT output score %d", iHist),
                                             500, 0., 50., 1000, 0., 1.); //assuming always probability
        fOutput->Add(fHistBDTOutputVsPt[iHist]);
    }

    //ML model
    double massD = -1.;
    switch(fDecChannel) {
        case kDplustoKpipi:
            fMLResponse = new AliHFMLResponseDplustoKpipi("DplustoKpipiMLResponse", "DplustoKpipiMLResponse", fConfigPath.Data());
            fMLResponse->MLResponseInit();
            massD = TDatabasePDG::Instance()->GetParticle(411)->Mass();
            break;
        case kDstoKKpi:
            fMLResponse = new AliHFMLResponseDstoKKpi("DstoKKpiMLResponse", "DstoKKpiMLResponse", fConfigPath.Data());
            fMLResponse->MLResponseInit();
            massD = TDatabasePDG::Instance()->GetParticle(431)->Mass();
            break;
        case kDstartoD0pi:
            fMLResponse = new AliHFMLResponseDstartoD0pi("DstartoD0piMLResponse", "DstartoD0piMLResponse", fConfigPath.Data());
            fMLResponse->MLResponseInit();
            massD = TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass();
            break; 
        case kD0toKpi:
            fMLResponse = new AliHFMLResponseD0toKpi("D0toKpiMLResponse", "D0toKpiMLResponse", fConfigPath.Data());
            fMLResponse->MLResponseInit();
            massD = TDatabasePDG::Instance()->GetParticle(421)->Mass();
            break; 
    }

    double minMass = massD-0.2;
    double maxMass = massD+0.2;
    if(fDecChannel == kDstartoD0pi)
    {
        minMass = 0.138;
        maxMass = 0.165;
    }

    fHistMassVsPt = new TH2F("fHistMassVsPt", ";#it{p}_{T} (GeV/#it{c});inv mass (GeV/#it{c}^{2})", 500, 0., 50., 200, minMass, maxMass);
    fOutput->Add(fHistMassVsPt);

    PostData(1, fOutput);

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronMLSelector::UserExec(Option_t * /*option*/)
{
    // main event loop
    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());
    if (!fAOD && AODEvent() && IsStandardAOD())
    {
        // In case there is an AOD handler writing a standard AOD, use the AOD
        // event in memory rather than the input (ESD) event.
        fAOD = dynamic_cast<AliAODEvent *>(AODEvent());
    }
    if (!fAOD)
    {
        AliWarning("AliAnalysisTaskSECharmHadronMLSelector::Exec(): bad AOD");
        PostData(1, fOutput);
        return;
    }
    fHistNEvents->Fill(1);

    fHistNEvents->Fill(0); // all events
    if (fAODProtection >= 0)
    {
        //   Protection against different number of events in the AOD and deltaAOD
        //   In case of discrepancy the event is rejected.
        int matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
        if (matchingAODdeltaAODlevel < 0 || (matchingAODdeltaAODlevel == 0 && fAODProtection == 1))
        {
            // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
            fHistNEvents->Fill(3);
            PostData(1, fOutput);
            return;
        }
        fHistNEvents->Fill(2);
    }

    // minimal event selection
    if (TMath::Abs(fAOD->GetMagneticField()) < 0.001)
    {
        PostData(1, fOutput);
        return;
    }
    TClonesArray *arrayCand = nullptr;
    TClonesArray *arrayCandDDau = nullptr;
    int absPdgMom = 0;
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
                    arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("D0toKpi"));
                    break;
                case kDplustoKpipi:
                    arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Charm3Prong"));
                    break;
                case kDstartoD0pi:
                    arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Dstar"));
                    arrayCandDDau = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("D0toKpi"));
                    break;
            }
        }
    }
    else if (fAOD)
    {
        switch (fDecChannel)
        {
            case kD0toKpi:
                arrayCand = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("D0toKpi"));
                break;
            case kDplustoKpipi:
                arrayCand = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Charm3Prong"));
                break;
            case kDstartoD0pi:
                arrayCand = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Dstar"));
                arrayCandDDau = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("D0toKpi"));
                break;
        }
    }

    if (!fAOD || !arrayCand || (fDecChannel == kDstartoD0pi && !arrayCandDDau))
    {
        AliWarning("Candidate branch not found!\n");
        PostData(1, fOutput);
        return;
    }

    switch(fDecChannel)
    {
        case kD0toKpi:
            absPdgMom = 421;
            break;
        case kDplustoKpipi:
            absPdgMom = 411;
            break;
        case kDstartoD0pi:
            absPdgMom = 413;
            break;
        case kDstoKKpi:
            absPdgMom = 431;
            break;
    }

    AliAODHandler *aodHandler = static_cast<AliAODHandler *>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    if (!aodHandler)
    {
        AliWarning("AliAnalysisTaskSECharmHadronMLSelector::Exec(): No AliInputEventHandler!");
        PostData(1, fOutput);
        return;
    }

    unsigned int maskPhysSel = aodHandler->IsEventSelected();
    TString firedTriggerClasses = fAOD->GetFiredTriggerClasses();
    if ((fAOD->GetRunNumber() < 136851 || fAOD->GetRunNumber() > 139517))
    {
        if (!(firedTriggerClasses.Contains(fTriggerClass.Data())))
        {
            fHistNEvents->Fill(5);
            PostData(1, fOutput);
            return;
        }
    }
    if ((maskPhysSel & fTriggerMask) == 0)
    {
        fHistNEvents->Fill(5);
        PostData(1, fOutput);
        return;
    }

    const AliAODVertex *vertex = fAOD->GetPrimaryVertex();
    if (!vertex || TMath::Abs(vertex->GetZ()) > 10. || vertex->GetNContributors() <= 0)
    {
        fHistNEvents->Fill(6);
        PostData(1, fOutput);
        return;
    }

    fHistNEvents->Fill(4);

    // needed to initialise PID response
    fRDCuts->IsEventSelected(fAOD);

    // select candidates
    fChHadIdx.clear();
    fMLScores.clear();
    fMLScoresSecond.clear();
    AliAnalysisVertexingHF vHF = AliAnalysisVertexingHF();

    for(int iCand = 0; iCand < arrayCand->GetEntriesFast(); iCand++)
    {
        AliAODRecoDecayHF *chHad = dynamic_cast<AliAODRecoDecayHF *>(arrayCand->UncheckedAt(iCand));
        AliAODRecoDecayHF *chHadWithVtx;
        if(fDecChannel == kDstartoD0pi)
        {
            if(chHad->GetIsFilled()<1)
                chHadWithVtx = dynamic_cast<AliAODRecoDecayHF *>(arrayCandDDau->UncheckedAt(chHad->GetProngID(1)));
            else
                chHadWithVtx = dynamic_cast<AliAODRecoDecayHF *>(((AliAODRecoCascadeHF *)chHad)->Get2Prong());
        }
        else
        {
            chHadWithVtx = chHad;
        }

        bool unsetVtx = false;
        bool recVtx = false;
        AliAODVertex *origOwnVtx = nullptr;

        std::vector<double> scores{}, scoresSecond{};
        int isSelected = IsCandidateSelected(chHad, chHadWithVtx, &vHF, absPdgMom, unsetVtx, recVtx, origOwnVtx, scores, scoresSecond);

        if (!isSelected || (fDecChannel == kDstoKKpi && !((isSelected & 4) || (isSelected & 8))))
        {
            if (unsetVtx)
                chHadWithVtx->UnsetOwnPrimaryVtx();
            if (recVtx)
                fRDCuts->CleanOwnPrimaryVtx(chHadWithVtx, fAOD, origOwnVtx);
            continue;
        }

        fChHadIdx.push_back(iCand);
        fMLScores.push_back(scores);
        fMLScoresSecond.push_back(scoresSecond);

        for(size_t iScore = 0; iScore < scores.size(); iScore++)
        {
            if(iScore > 2)
                break;
            if((fDecChannel == kD0toKpi && (isSelected == 1 || isSelected == 3)) || fDecChannel == kDplustoKpipi || fDecChannel == kDstartoD0pi || (fDecChannel == kDstoKKpi && (isSelected & 4)))
                fHistBDTOutputVsPt[iScore]->Fill(chHad->Pt(), scores[iScore]);
            if((fDecChannel == kD0toKpi && (isSelected >= 2)) || (fDecChannel == kDstoKKpi && (isSelected & 8)))
                fHistBDTOutputVsPt[iScore]->Fill(chHad->Pt(), scoresSecond[iScore]);
        }
        switch(fDecChannel)
        {
            case kD0toKpi:
                if(isSelected == 1 || isSelected == 3)
                    fHistMassVsPt->Fill(chHad->Pt(), dynamic_cast<AliAODRecoDecayHF2Prong*>(chHad)->InvMassD0());
                if(isSelected >= 2)
                    fHistMassVsPt->Fill(chHad->Pt(), dynamic_cast<AliAODRecoDecayHF2Prong*>(chHad)->InvMassD0bar());
                break;
            case kDplustoKpipi:
                fHistMassVsPt->Fill(chHad->Pt(), dynamic_cast<AliAODRecoDecayHF3Prong*>(chHad)->InvMassDplus());
                break;
            case kDstartoD0pi:
                fHistMassVsPt->Fill(chHad->Pt(), dynamic_cast<AliAODRecoCascadeHF*>(chHad)->DeltaInvMass());
                break;
            case kDstoKKpi:
                if(isSelected & 4)
                    fHistMassVsPt->Fill(chHad->Pt(), dynamic_cast<AliAODRecoDecayHF3Prong*>(chHad)->InvMassDsKKpi());
                if(isSelected & 8)
                    fHistMassVsPt->Fill(chHad->Pt(), dynamic_cast<AliAODRecoDecayHF3Prong*>(chHad)->InvMassDspiKK());
                break;
        }        

        if (unsetVtx)
            chHadWithVtx->UnsetOwnPrimaryVtx();
        if (recVtx)
            fRDCuts->CleanOwnPrimaryVtx(chHadWithVtx, fAOD, origOwnVtx);
    }

    fHistNselCand->Fill(fChHadIdx.size());
    fHistNallCand->Fill(arrayCand->GetEntriesFast());

    PostData(1, fOutput);
}

//________________________________________________________________________
int AliAnalysisTaskSECharmHadronMLSelector::IsCandidateSelected(AliAODRecoDecayHF *&chHad, AliAODRecoDecayHF *&chHadWithVtx, AliAnalysisVertexingHF *vHF,
                                                                int absPdgMom, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx,
                                                                std::vector<double> &modelPred, std::vector<double> &modelPredSecond)
{
    if(!chHad || !chHadWithVtx || !vHF )
        return 0;

    //Preselection to speed up task
    TObjArray arrDauTracks(3);
    int nDau = 3;
    if (fDecChannel == kD0toKpi)
        nDau = 2;

    for (int iDau = 0; iDau < nDau; iDau++)
    {
        AliAODTrack *track;
        if (fDecChannel != kDstartoD0pi || iDau == 0)
            track = vHF->GetProng(fAOD, chHad, iDau);
        else
            track = vHF->GetProng(fAOD, chHadWithVtx, iDau-1); //D0<-D* daughters
        arrDauTracks.AddAt(track, iDau);
    }

    if(!fRDCuts->PreSelect(arrDauTracks))
        return 0;

    bool isSelBit = true;
    switch (fDecChannel)
    {
        case kD0toKpi:
            isSelBit = chHad->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts);
            if (!isSelBit || !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF2Prong *>(chHad)))
                return 0;
            break;
        case kDplustoKpipi:
            isSelBit = chHad->HasSelectionBit(AliRDHFCuts::kDplusCuts);
            if (!isSelBit || !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF3Prong *>(chHad)))
                return 0;
            break;
        case kDstartoD0pi:
            if (!vHF->FillRecoCasc(fAOD, dynamic_cast<AliAODRecoCascadeHF *>(chHad), false))
                return 0;
            break;
        case kDstoKKpi:
            isSelBit = chHad->HasSelectionBit(AliRDHFCuts::kDsCuts);
            if (!isSelBit || !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF3Prong *>(chHad)))
                return 0;
            break;

    }

    unsetVtx = false;
    if (!chHadWithVtx->GetOwnPrimaryVtx())
    {
        chHadWithVtx->SetOwnPrimaryVtx(dynamic_cast<AliAODVertex *>(fAOD->GetPrimaryVertex()));
        unsetVtx = true;
        // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
        // Pay attention if you use continue inside this loop!!!
    }

    double ptD = chHad->Pt();
    int ptbin = fRDCuts->PtBin(ptD);
    if (ptbin < 0)
    {
        if (unsetVtx)
            chHadWithVtx->UnsetOwnPrimaryVtx();
        return 0;
    }

    int isSelected = fRDCuts->IsSelected(chHad, AliRDHFCuts::kAll, fAOD);
    if (!isSelected)
    {
        if (unsetVtx)
            chHadWithVtx->UnsetOwnPrimaryVtx();
        return 0;
    }

    recVtx = false;
    origOwnVtx = nullptr;

    if (fRDCuts->GetIsPrimaryWithoutDaughters())
    {
        if (chHadWithVtx->GetOwnPrimaryVtx())
            origOwnVtx = new AliAODVertex(*chHadWithVtx->GetOwnPrimaryVtx());
        if (fRDCuts->RecalcOwnPrimaryVtx(chHadWithVtx, fAOD))
            recVtx = true;
        else
            fRDCuts->CleanOwnPrimaryVtx(chHadWithVtx, fAOD, origOwnVtx);
    }
    
    // ML application
    AliAODPidHF *pidHF = fRDCuts->GetPidHF();
    int isMLsel = 0;
    modelPred = {};
    modelPredSecond = {};
    if((fDecChannel == kD0toKpi && (isSelected == 1 || isSelected == 3)) || fDecChannel == kDplustoKpipi || fDecChannel == kDstartoD0pi || (fDecChannel == kDstoKKpi && (isSelected & 4)))
    {
        if(fMLResponse->IsSelectedMultiClass(modelPred, chHad, fAOD->GetMagneticField(), pidHF, 0))
            isMLsel += (fDecChannel != kDstoKKpi) ? 1 : 4;
    }
    if((fDecChannel == kD0toKpi && (isSelected >= 2)) || (fDecChannel == kDstoKKpi && (isSelected & 8)))
    {
        if(fMLResponse->IsSelectedMultiClass(modelPredSecond, chHad, fAOD->GetMagneticField(), pidHF, 1))
            isMLsel += (fDecChannel != kDstoKKpi) ? 2 : 8;
    }

    if(modelPred.size() > modelPredSecond.size())
        for(int iScore=0; iScore<modelPred.size(); iScore++)
            modelPredSecond.push_back(-9999.);
    else if(modelPred.size() < modelPredSecond.size())
        for(int iScore=0; iScore<modelPredSecond.size(); iScore++)
            modelPred.push_back(-9999.);

    isSelected = isMLsel;

    return isSelected;
}
