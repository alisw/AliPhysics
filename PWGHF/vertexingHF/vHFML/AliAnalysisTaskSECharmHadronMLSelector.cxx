/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSECharmHadronMLSelector
// \brief Analysis task to select charm-hadron candidates based on ML response for subsequent tasks
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch

#include "AliAnalysisTaskSECharmHadronMLSelector.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliHFMLResponseDplustoKpipi.h"
#include "AliHFMLResponseDstoKKpi.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF3Prong.h"
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
    DefineOutput(2, TList::Class());
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
        fHistBDTOutput[iHist] = new TH1F(Form("fHistBDTOutput%d", iHist), Form(";BDT output score %d; entries", iHist), 1000, 0., 1.); //assuming always probability
        fOutput->Add(fHistBDTOutput[iHist]);
    }

    //ML model
    switch(fDecChannel) {
        case kDplustoKpipi:
            fMLResponse = new AliHFMLResponseDplustoKpipi("DplustoKpipiMLResponse", "DplustoKpipiMLResponse", fConfigPath.Data());
            fMLResponse->MLResponseInit();
            break;
        case kDstoKKpi:
            fMLResponse = new AliHFMLResponseDstoKKpi("DstoKKpiMLResponse", "DstoKKpiMLResponse", fConfigPath.Data());
            fMLResponse->MLResponseInit();
            break;
    }

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
            arrayCand = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Charm3Prong"));
        }
    }
    else if (fAOD)
    {
        arrayCand = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Charm3Prong"));
    }

    if (!fAOD || !arrayCand)
    {
        AliWarning("Candidate branch not found!\n");
        PostData(1, fOutput);
        return;
    }

    if(fDecChannel == kDplustoKpipi)
        absPdgMom = 411;
    else if(fDecChannel == kDstoKKpi)
        absPdgMom = 431;

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
            return;
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

    // select candidates
    fChHadIdx.clear();
    fMLScores.clear();
    AliAnalysisVertexingHF vHF = AliAnalysisVertexingHF();

    for(int iCand = 0; iCand < arrayCand->GetEntriesFast(); iCand++)
    {
        AliAODRecoDecayHF *chHad = dynamic_cast<AliAODRecoDecayHF *>(arrayCand->UncheckedAt(iCand));

        bool unsetVtx = false;
        bool recVtx = false;
        AliAODVertex *origOwnVtx = nullptr;

        std::vector<double> scores{};
        int isSelected = IsCandidateSelected(chHad, &vHF, absPdgMom, unsetVtx, recVtx, origOwnVtx, scores);

        if (!isSelected || (fDecChannel == kDstoKKpi && isSelected % 2 == 0))
        {
            if (unsetVtx)
                chHad->UnsetOwnPrimaryVtx();
            if (recVtx)
                fRDCuts->CleanOwnPrimaryVtx(chHad, fAOD, origOwnVtx);
            continue;
        }

        fChHadIdx.push_back(iCand);
        fMLScores.push_back(scores);

        for(size_t iScore = 0; iScore < scores.size(); iScore++)
        {
            if(iScore > 2)
                break;
            fHistBDTOutput[iScore]->Fill(scores[iScore]);
        }

        if (unsetVtx)
            chHad->UnsetOwnPrimaryVtx();
        if (recVtx)
            fRDCuts->CleanOwnPrimaryVtx(chHad, fAOD, origOwnVtx);
    }

    fHistNselCand->Fill(fChHadIdx.size());
    fHistNallCand->Fill(arrayCand->GetEntriesFast());

    PostData(1, fOutput);
}

//________________________________________________________________________
int AliAnalysisTaskSECharmHadronMLSelector::IsCandidateSelected(AliAODRecoDecayHF *&chHad, AliAnalysisVertexingHF *vHF,
                                                                int absPdgMom, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx, std::vector<double> modelPred)
{
    if(!chHad || !vHF )
        return 0;

    //Preselection to speed up task
    TObjArray arrDauTracks(3);
    for(int iDau=0; iDau < 3; iDau++){
        AliAODTrack *track = vHF->GetProng(fAOD, chHad, iDau);
        arrDauTracks.AddAt(track,iDau);
    }

    if(!fRDCuts->PreSelect(arrDauTracks))
        return 0;

    bool isSelBit = true;
    switch (fDecChannel)
    {
        case kDplustoKpipi:
            isSelBit = chHad->HasSelectionBit(AliRDHFCuts::kDplusCuts);
            if (!isSelBit || !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF3Prong *>(chHad)))
                return 0;
            break;
        case kDstoKKpi:
            isSelBit = chHad->HasSelectionBit(AliRDHFCuts::kDsCuts);
            if (!isSelBit || !vHF->FillRecoCand(fAOD, dynamic_cast<AliAODRecoDecayHF3Prong *>(chHad)))
                return 0;
            break;
    }

    unsetVtx = false;
    if (!chHad->GetOwnPrimaryVtx())
    {
        chHad->SetOwnPrimaryVtx(dynamic_cast<AliAODVertex *>(fAOD->GetPrimaryVertex()));
        unsetVtx = true;
        // NOTE: the own primary vertex should be unset, otherwise there is a memory leak
        // Pay attention if you use continue inside this loop!!!
    }

    double ptD = chHad->Pt();
    double yD = chHad->Y(absPdgMom);
    int ptbin = fRDCuts->PtBin(ptD);
    if (ptbin < 0)
    {
        if (unsetVtx)
            chHad->UnsetOwnPrimaryVtx();
        return 0;
    }

    bool isFidAcc = fRDCuts->IsInFiducialAcceptance(ptD, yD);
    if (!isFidAcc)
    {
        if (unsetVtx)
            chHad->UnsetOwnPrimaryVtx();
        return 0;
    }

    int isSelected = fRDCuts->IsSelected(chHad, AliRDHFCuts::kAll, fAOD);
    if (!isSelected)
    {
        if (unsetVtx)
            chHad->UnsetOwnPrimaryVtx();
        return 0;
    }

    // ML application
    AliAODPidHF *pidHF = fRDCuts->GetPidHF();
    bool isMLsel = true;
    isMLsel = fMLResponse->IsSelectedMultiClass(modelPred, chHad, fAOD->GetMagneticField(), pidHF);
    if (!isMLsel)
        isSelected = 0;

    recVtx = false;
    origOwnVtx = nullptr;

    if (fRDCuts->GetIsPrimaryWithoutDaughters())
    {
        if (chHad->GetOwnPrimaryVtx())
            origOwnVtx = new AliAODVertex(*chHad->GetOwnPrimaryVtx());
        if (fRDCuts->RecalcOwnPrimaryVtx(chHad, fAOD))
            recVtx = true;
        else
            fRDCuts->CleanOwnPrimaryVtx(chHad, fAOD, origOwnVtx);
    }

    return isSelected;
}
