// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliAnalysisTaskSECharmTriggerStudy
// \brief task that produces an output tree for the charm trigger studies
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TList.h>
#include <TString.h>
#include <TH1F.h>
#include <TDatabasePDG.h>
#include <TObjArray.h>

#include "AliAnalysisTaskSECharmTriggerStudy.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliRDHFCuts.h"
#include "AliLog.h"
#include "AliVertexerTracks.h"
#include "AliESDVertex.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSECharmTriggerStudy);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSECharmTriggerStudy::AliAnalysisTaskSECharmTriggerStudy() : AliAnalysisTaskSECharmTriggerStudy("", nullptr)
{
    /// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSECharmTriggerStudy::AliAnalysisTaskSECharmTriggerStudy(const char *name, TList *cutlist) : AliAnalysisTaskSE(name),
                                                                                                           fOutput(nullptr),
                                                                                                           fHistNEvents(nullptr),
                                                                                                           fRecoTree(nullptr),
                                                                                                           fGenTree(nullptr),
                                                                                                           fEventCuts{},
                                                                                                           fSystem(kpp),
                                                                                                           fAOD(nullptr),
                                                                                                           fAODProtection(1),
                                                                                                           fMCArray(nullptr),
                                                                                                           fRecoZvtx(-999.),
                                                                                                           fGenZvtx(-999.),
                                                                                                           fNtracklets(-1.),
                                                                                                           fCharm2Prong{},
                                                                                                           fCharm3Prong{},
                                                                                                           fDstar{},
                                                                                                           fCharmCascade{},
                                                                                                           fBeauty3Prong{},
                                                                                                           fBeauty4Prong{},
                                                                                                           fGenHadron{},
                                                                                                           fEnable2Prongs(true),
                                                                                                           fEnable3Prongs(0),
                                                                                                           fEnableDstars(false),
                                                                                                           fEnableCascades(false),
                                                                                                           fEnableBeauty3Prongs(false),
                                                                                                           fEnableBeauty4Prongs(false),
                                                                                                           fFillOnlySignal(false),
                                                                                                           fFillGenTree(true),
                                                                                                           fReadMC(true),
                                                                                                           fCutsD0toKpi(nullptr),
                                                                                                           fCutsDplustoKpipi(nullptr),
                                                                                                           fCutsDstartoKpipi(nullptr),
                                                                                                           fCutsDstoKKpi(nullptr),
                                                                                                           fCutsLctopKpi(nullptr),
                                                                                                           fCutsLctoV0bach(nullptr),
                                                                                                           fApplyCuts(false),
                                                                                                           fListCuts(cutlist)
{
    /// Default constructor
    AliRDHFCutsD0toKpi *cutsD0 = nullptr;
    AliRDHFCutsDplustoKpipi *cutsDplus = nullptr;
    AliRDHFCutsDstoKKpi *cutsDs = nullptr;
    AliRDHFCutsLctopKpi *cutsLc = nullptr;
    AliRDHFCutsDStartoKpipi *cutsDstar = nullptr;
    AliRDHFCutsLctoV0 *cutsLctoV0 = nullptr;

    if (fListCuts)
    {
        cutsD0 = static_cast<AliRDHFCutsD0toKpi *>(fListCuts->FindObject("D0toKpiCuts"));
        cutsDplus = static_cast<AliRDHFCutsDplustoKpipi *>(fListCuts->FindObject("DplustoKpipiCuts"));
        cutsDs = static_cast<AliRDHFCutsDstoKKpi *>(fListCuts->FindObject("DstoKKpiCuts"));
        cutsLc = static_cast<AliRDHFCutsLctopKpi *>(fListCuts->FindObject("LctopKpiCuts"));
        cutsDstar = static_cast<AliRDHFCutsDStartoKpipi *>(fListCuts->FindObject("DstartoKpipiCuts"));
        cutsLctoV0 = static_cast<AliRDHFCutsLctoV0 *>(fListCuts->FindObject("LctoV0bachCuts"));
    }

    if (cutsD0)
        fCutsD0toKpi = static_cast<AliRDHFCutsD0toKpi *>(cutsD0->Clone("fCutsD0toKpi"));
    if (cutsDplus)
        fCutsDplustoKpipi = static_cast<AliRDHFCutsDplustoKpipi *>(cutsDplus->Clone("fCutsDplustoKpipi"));
    if (cutsDs)
        fCutsDstoKKpi = static_cast<AliRDHFCutsDstoKKpi *>(cutsDs->Clone("fCutsDstoKKpi"));
    if (cutsLc)
        fCutsLctopKpi = static_cast<AliRDHFCutsLctopKpi *>(cutsLc->Clone("fCutsLctopKpi"));
    if (cutsDstar)
        fCutsDstartoKpipi = static_cast<AliRDHFCutsDStartoKpipi *>(cutsDstar->Clone("fCutsDstartoKpipi"));
    if (cutsLctoV0)
        fCutsLctoV0bach = static_cast<AliRDHFCutsLctoV0 *>(cutsLctoV0->Clone("fCutsLctoV0bach"));

    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
    DefineOutput(3, TTree::Class());
    DefineOutput(4, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskSECharmTriggerStudy::~AliAnalysisTaskSECharmTriggerStudy()
{
    // Destructor
    if (fOutput)
        delete fOutput;
    if (fRecoTree)
        delete fRecoTree;
    if (fGenTree)
        delete fGenTree;

    if (fCutsD0toKpi)
        delete fCutsD0toKpi;
    if (fCutsDplustoKpipi)
        delete fCutsDplustoKpipi;
    if (fCutsDstoKKpi)
        delete fCutsDstoKKpi;
    if (fCutsDstartoKpipi)
        delete fCutsDstartoKpipi;
    if (fCutsLctopKpi)
        delete fCutsLctopKpi;
    if (fCutsLctoV0bach)
        delete fCutsLctoV0bach;
    if (fListCuts)
        delete fListCuts;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::Init()
{
    /// Initialization

    PostData(4, fListCuts);

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::UserCreateOutputObjects()
{
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos");
    fEventCuts.AddQAplotsToList(fOutput);

    fHistNEvents = new TH1F("hNEvents", "number of events ", 10, 0.5, 10.5);
    fHistNEvents->GetXaxis()->SetBinLabel(1, "nEventsRead");
    fHistNEvents->GetXaxis()->SetBinLabel(2, "nEvents Matched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(3, "nEvents Mismatched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(4, "nEventsAnal");
    fHistNEvents->GetXaxis()->SetBinLabel(5, "n. rejected due to not reco vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(6, "n. passing IsEvSelected");
    fHistNEvents->GetXaxis()->SetBinLabel(7, "no. of 2 prong candidates");
    fHistNEvents->GetXaxis()->SetBinLabel(8, "no. of 3 prong candidates");
    fHistNEvents->GetXaxis()->SetBinLabel(9, "no. of Dstar candidates");
    fHistNEvents->GetXaxis()->SetBinLabel(10, "no. of cascade candidates");

    fHistNEvents->GetXaxis()->SetNdivisions(1, false);
    fHistNEvents->SetMinimum(0);

    fOutput->Add(fHistNEvents);

    fRecoTree = new TTree("fRecoTree", "Reconstructed charm hadron candidates");
    fRecoTree->SetMaxVirtualSize(1.e+8);
    fRecoTree->Branch("zVtxReco", &fRecoZvtx);
    fRecoTree->Branch("Ntracklets", &fNtracklets);
    if (fEnable2Prongs)
        fRecoTree->Branch("Charm2Prong", &fCharm2Prong);
    if (fEnable3Prongs)
        fRecoTree->Branch("Charm3Prong", &fCharm3Prong);
    if (fEnableDstars)
        fRecoTree->Branch("Dstar", &fDstar);
    if (fEnableCascades)
        fRecoTree->Branch("CharmCascade", &fCharmCascade);
    if (fEnableBeauty3Prongs)
        fRecoTree->Branch("Beauty3Prong", &fBeauty3Prong);
    if (fEnableBeauty4Prongs)
        fRecoTree->Branch("Beauty4Prong", &fBeauty4Prong);

    fGenTree = new TTree("fGenTree", "Generate charm hadrons");
    fGenTree->SetMaxVirtualSize(1.e+8);
    if(fReadMC && fFillGenTree)
    {
        fGenTree->Branch("zVtxGen", &fGenZvtx);
        fGenTree->Branch("GenHadron", &fGenHadron);
    }

    PostData(1, fOutput);
    PostData(2, fRecoTree);
    PostData(3, fGenTree);

    //cut objects
    float ptbins[2] = {0., 100.};
    if (fEnable2Prongs || fEnableBeauty3Prongs)
    {
        if (!fCutsD0toKpi)
        { //cut object for PID
            fCutsD0toKpi = new AliRDHFCutsD0toKpi("fCutsD0toKpi");
            fCutsD0toKpi->SetPtBins(2, ptbins);
            fCutsD0toKpi->SetUsePID(true);
            AliAODPidHF *pidObj = new AliAODPidHF();
            int mode = 1;
            const int nlims = 2;
            double plims[nlims] = {0.6, 0.8};
            bool compat = true;
            bool asym = true;
            double sigmas[5] = {2., 1., 0., 3., 0.};
            pidObj->SetAsym(asym);
            pidObj->SetMatch(mode);
            pidObj->SetPLimit(plims, nlims);
            pidObj->SetSigma(sigmas);
            pidObj->SetCompat(compat);
            pidObj->SetTPC(true);
            pidObj->SetTOF(true);
            fCutsD0toKpi->SetPidHF(pidObj);
            fCutsD0toKpi->SetUseDefaultPID(false);
        }
    }
    if (fEnable3Prongs || fEnableBeauty4Prongs)
    {
        if (!fCutsDplustoKpipi)
        { //cut object for PID
            fCutsDplustoKpipi = new AliRDHFCutsDplustoKpipi("fCutsDplustoKpipi");
            fCutsDplustoKpipi->SetPtBins(2, ptbins);
            fCutsDplustoKpipi->SetUsePID(true);
        }

        if (!fCutsDstoKKpi)
        { //cut object for PID
            fCutsDstoKKpi = new AliRDHFCutsDstoKKpi("fCutsDstoKKpi");
            fCutsDstoKKpi->SetPtBins(2, ptbins);
            fCutsDstoKKpi->SetUsePID(true);
        }

        if (!fCutsLctopKpi)
        { //cut object for PID
            fCutsLctopKpi = new AliRDHFCutsLctopKpi("fCutsLctopKpi");
            fCutsLctopKpi->SetPtBins(2, ptbins);
            fCutsLctopKpi->SetUsePID(true);
            AliAODPidHF *pidObjp = new AliAODPidHF();
            AliAODPidHF *pidObjK = new AliAODPidHF();
            AliAODPidHF *pidObjpi = new AliAODPidHF();
            pidObjp->SetMatch(1);
            pidObjK->SetMatch(1);
            pidObjpi->SetMatch(1);
            pidObjp->SetTPC(true);
            pidObjK->SetTPC(true);
            pidObjpi->SetTPC(true);
            pidObjp->SetTOF(true);
            pidObjK->SetTOF(true);
            pidObjpi->SetTOF(true);
            fCutsLctopKpi->SetPidprot(pidObjp);
            fCutsLctopKpi->SetPidHF(pidObjK);
            fCutsLctopKpi->SetPidpion(pidObjpi);
            fCutsLctopKpi->GetPidHF()->SetCombDetectors(AliAODPidHF::kTPCTOF);
            for (int iSpecies = 0; iSpecies < AliPID::kSPECIES; ++iSpecies)
                fCutsLctopKpi->SetPIDThreshold(static_cast<AliPID::EParticleType>(iSpecies), 0.);
            fCutsLctopKpi->GetPidHF()->SetUseCombined(true);
            fCutsLctopKpi->GetPidHF()->SetUpCombinedPID();
            fCutsLctopKpi->SetPIDStrategy(AliRDHFCutsLctopKpi::kCombinedpPb);
        }
    }
    if (fEnableDstars)
    {
        if (!fCutsDstartoKpipi)
        { //cut object for PID
            fCutsDstartoKpipi = new AliRDHFCutsDStartoKpipi("fCutsDstartoKpipi");
            fCutsDstartoKpipi->SetPtBins(2, ptbins);
            fCutsDstartoKpipi->SetUsePID(true);
            AliAODPidHF *pidObj = new AliAODPidHF();
            pidObj->SetMatch(1);
            pidObj->SetSigma(0, 2);
            pidObj->SetSigma(3, 3);
            pidObj->SetTPC(true);
            pidObj->SetTOF(true);
            fCutsDstartoKpipi->SetPidHF(pidObj);
        }
    }
    if (fEnableCascades)
    {
        if (!fCutsLctoV0bach)
        { //cut object for PID
            fCutsLctoV0bach = new AliRDHFCutsLctoV0("LctoV0bachCuts");
            fCutsLctoV0bach->SetPtBins(2, ptbins);
            fCutsLctoV0bach->SetUsePID(true);
        }
    }

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::UserExec(Option_t * /*option*/)
{

    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());

    fHistNEvents->Fill(1); // all events
    if (fAODProtection >= 0)
    {
        int matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
        if (matchingAODdeltaAODlevel < 0 || (matchingAODdeltaAODlevel == 0 && fAODProtection == 1))
        {
            fHistNEvents->Fill(3);
            PostData(1, fOutput);
            return;
        }
        fHistNEvents->Fill(2);
    }

    TClonesArray *array3Prong = nullptr, *array2Prong = nullptr, *arrayCasc = nullptr, *arrayDstar = nullptr;

    if (!fAOD && AODEvent() && IsStandardAOD())
    {
        fAOD = dynamic_cast<AliAODEvent *>(AODEvent());
        AliAODHandler *aodHandler = dynamic_cast<AliAODHandler *>((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
        if (aodHandler->GetExtensions())
        {
            AliAODExtension *ext = dynamic_cast<AliAODExtension *>(aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root"));
            AliAODEvent *aodFromExt = ext->GetAOD();
            array2Prong = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("D0toKpi"));
            array3Prong = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Charm3Prong"));
            arrayDstar = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Dstar"));
            arrayCasc = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("CascadesHF"));
        }
    }
    else if (fAOD)
    {
        array2Prong = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("D0toKpi"));
        array3Prong = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Charm3Prong"));
        arrayDstar = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Dstar"));
        arrayCasc = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("CascadesHF"));
    }

    if (!fAOD || ((fEnable3Prongs || fEnableBeauty4Prongs) && !array3Prong) || ((fEnable2Prongs || fEnableBeauty3Prongs) && !array2Prong) || (fEnableCascades && !arrayCasc) || (fEnableDstars && !arrayDstar))
    {
        AliWarning("Candidate branch not found!");
        return;
    }

    if (!fAOD->GetPrimaryVertex() || TMath::Abs(fAOD->GetMagneticField()) < 0.001)
        return;

    fHistNEvents->Fill(4); // count event

    fEventCuts.AcceptEvent(fAOD); // do not return yet (no physics selection applied for upgrade MC)

    if (!fEventCuts.PassedCut(AliEventCuts::kVertex) || !fEventCuts.PassedCut(AliEventCuts::kVertexPosition) || !fEventCuts.PassedCut(AliEventCuts::kVertexQuality))
    {
        fHistNEvents->Fill(5); // rejected for primary vtx
        PostData(1, fOutput);
        return;
    }

    fHistNEvents->Fill(6); // selected event

    AliAODMCHeader *mcHeader = nullptr;
    if(fReadMC)
    {
        fMCArray = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
        if (!fMCArray)
        {
            AliWarning("MC particles branch not found!");
            return;
        }

        // load MC header
        mcHeader = dynamic_cast<AliAODMCHeader *>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
        if (!mcHeader)
        {
            AliWarning("MC header branch not found!");
            return;
        }

        fGenZvtx = mcHeader->GetVtxZ();
    }

    AliAODVertex *primVtx = dynamic_cast<AliAODVertex *>(fAOD->GetPrimaryVertex());
    fRecoZvtx = primVtx->GetZ();
    fNtracklets = AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(fAOD,-1.,1.);

    if (fEnable2Prongs || fEnableBeauty3Prongs)
    {
        fCutsD0toKpi->SetupPID(fAOD);
    }
    if (fEnable3Prongs || fEnableBeauty4Prongs)
    {
        fCutsDplustoKpipi->SetupPID(fAOD);
        fCutsDstoKKpi->SetupPID(fAOD);
        fCutsLctopKpi->SetupPID(fAOD);
    }
    if (fEnableDstars)
    {
        fCutsDstartoKpipi->SetupPID(fAOD);
    }
    if (fEnableCascades)
    {
        fCutsLctoV0bach->SetupPID(fAOD);
    }

    //loop on generated particles
    if(fReadMC && fFillGenTree)
    {
        for (int iPart = 0; iPart < fMCArray->GetEntriesFast(); iPart++)
        {
            AliAODMCParticle *part = dynamic_cast<AliAODMCParticle *>(fMCArray->UncheckedAt(iPart));
            if (!part)
                continue;

            int pdgCode = TMath::Abs(part->GetPdgCode());
            if (pdgCode == 411 || pdgCode == 421 || pdgCode == 431 || pdgCode == 413 || pdgCode == 4122 || pdgCode == 521 || pdgCode == 511 || pdgCode == 531 || pdgCode == 5122)
            {
                int origin = 4; //beauty assumed to be always prompt
                if (pdgCode != 521 && pdgCode != 511 && pdgCode != 531 && pdgCode != 5122)
                    origin = AliVertexingHFUtils::CheckOrigin(fMCArray, part, true);
                if (origin != 4 && origin != 5)
                    continue; //keep only prompt or feed-down

                int labDau[4] = {-1, -1, -1, -1};
                int pdgCodeDau0 = -1;
                int decay = -1;
                bool dauInAcc = false;
                if (pdgCode == 421 && fEnable2Prongs) //Dzero
                {
                    if (part->GetNDaughters() == 2)
                        decay = AliVertexingHFUtils::CheckD0Decay(fMCArray, part, labDau);
                    if (decay != 1 || labDau[0] < 0 || labDau[1] < 0)
                        continue;
                    dauInAcc = AreDauInAcc(2, labDau);
                    pdgCodeDau0 = (dynamic_cast<AliAODMCParticle *>(fMCArray->UncheckedAt(labDau[0])))->GetPdgCode();
                    if (TMath::Abs(pdgCodeDau0) == 321)
                        FillGenerated(part, origin, kDzerotoKpi, dauInAcc);
                    else
                        FillGenerated(part, origin, kDzerotopiK, dauInAcc);
                }
                else if (pdgCode == 411 && fEnable3Prongs >> 0 & 1) //D+
                {
                    decay = AliVertexingHFUtils::CheckDplusDecay(fMCArray, part, labDau);
                    if (decay >= 1 && labDau[0] >= 0 && labDau[1] >= 0)
                    {
                        dauInAcc = AreDauInAcc(3, labDau);
                        FillGenerated(part, origin, kDplustoKpipi, dauInAcc);
                        continue;
                    }
                }
                else if ((pdgCode == 431 || pdgCode == 411) && fEnable3Prongs >> 1 & 1) //Ds (or D+ --> KKpi)
                {
                    decay = AliVertexingHFUtils::CheckDsDecay(fMCArray, part, labDau);
                    int decayDplus = AliVertexingHFUtils::CheckDplusKKpiDecay(fMCArray, part, labDau);

                    if (labDau[0] < 0 || labDau[1] < 0 || (decay != 1 && decayDplus!= 1)) //keep only Ds -> phipi --> KKpi (to be discussed)
                        continue;

                    dauInAcc = AreDauInAcc(3, labDau);
                    pdgCodeDau0 = (dynamic_cast<AliAODMCParticle *>(fMCArray->UncheckedAt(labDau[0])))->GetPdgCode();
                    if(decay == 1)
                    {
                        if (TMath::Abs(pdgCodeDau0) == 321)
                            FillGenerated(part, origin, kDstoKKpi, dauInAcc);
                        else
                            FillGenerated(part, origin, kDstopiKK, dauInAcc);
                        continue;
                    }
                    if(decayDplus == 1)
                    {
                        if (TMath::Abs(pdgCodeDau0) == 321)
                            FillGenerated(part, origin, kDplustoKKpi, dauInAcc);
                        else
                            FillGenerated(part, origin, kDplustopiKK, dauInAcc);
                    }
                }
                else if (pdgCode == 413 && fEnableDstars) //Dstar
                {
                    decay = AliVertexingHFUtils::CheckDstarDecay(fMCArray, part, labDau);
                    if (decay != 1 || labDau[0] < 0 || labDau[1] < 0)
                        continue;
                    dauInAcc = AreDauInAcc(3, labDau);
                    FillGenerated(part, origin, kDstartoKpipi, dauInAcc);
                }
                else if (pdgCode == 4122) //Lc
                {
                    if (fEnable3Prongs >> 2 & 1)
                    {
                        decay = AliVertexingHFUtils::CheckLcpKpiDecay(fMCArray, part, labDau);
                        dauInAcc = AreDauInAcc(3, labDau);
                        if (decay >= 1 && labDau[0] >= 0 && labDau[1] >= 0)
                        {
                            pdgCodeDau0 = (dynamic_cast<AliAODMCParticle *>(fMCArray->UncheckedAt(labDau[0])))->GetPdgCode();
                            if (TMath::Abs(pdgCodeDau0) == 2212)
                                FillGenerated(part, origin, kLctopKpi, dauInAcc);
                            else
                                FillGenerated(part, origin, kLctopiKp, dauInAcc);
                            continue;
                        }
                    }
                    if (fEnableCascades)
                    {
                        decay = AliVertexingHFUtils::CheckLcV0bachelorDecay(fMCArray, part, labDau);
                        dauInAcc = AreDauInAcc(3, labDau);
                        if (labDau[0] >= 0 && labDau[1] >= 0)
                        {
                            if (decay == 1)
                                FillGenerated(part, origin, kLctopiLambda, dauInAcc);
                            else if (decay == 2)
                                FillGenerated(part, origin, kLctopK0s, dauInAcc);
                        }
                    }
                }
                if (pdgCode == 521 && fEnableBeauty3Prongs) //Bplus
                {
                    decay = AliVertexingHFUtils::CheckBplusDecay(fMCArray, part, labDau);
                    if (decay != 1 || labDau[0] == -1 || labDau[1] < 0)
                        continue;

                    dauInAcc = AreDauInAcc(3, labDau);
                    FillGenerated(part, origin, kBplustoD0pi, dauInAcc);
                }
                if(pdgCode == 511 && fEnableBeauty4Prongs >> 0 & 1) //B0
                {
                    decay = AliVertexingHFUtils::CheckB0toDminuspiDecay(fMCArray, part, labDau);
                    if (decay < 1 || labDau[0] == -1 || labDau[1] < 0)
                        continue;

                    dauInAcc = AreDauInAcc(4, labDau);
                    FillGenerated(part, origin, kB0toDminuspi, dauInAcc);
                }
                if(pdgCode == 531 && fEnableBeauty4Prongs >> 1 & 1) //Bs
                {
                    decay = AliVertexingHFUtils::CheckBsDecay(fMCArray, part, labDau);
                    if (decay < 1 || labDau[0] == -1 || labDau[1] < 0)
                        continue;

                    dauInAcc = AreDauInAcc(4, labDau);
                    FillGenerated(part, origin, kBstoDsminuspi, dauInAcc);
                }
                if(pdgCode == 5122 && fEnableBeauty4Prongs >> 2 & 1) //Lb
                {
                    decay = AliVertexingHFUtils::CheckLbDecay(fMCArray, part, labDau);
                    if (decay < 1 || labDau[0] == -1 || labDau[1] < 0)
                        continue;

                    dauInAcc = AreDauInAcc(4, labDau);
                    FillGenerated(part, origin, kLbtoLcpluspi, dauInAcc);
                }
            }
        }
    }

    AliAnalysisVertexingHF vHF;

    //loop on 2 prongs
    if (fEnable2Prongs || fEnableBeauty3Prongs)
    {
        for (int i2Prong = 0; i2Prong < array2Prong->GetEntriesFast(); i2Prong++)
        {
            AliAODRecoDecayHF2Prong *d = dynamic_cast<AliAODRecoDecayHF2Prong *>(array2Prong->UncheckedAt(i2Prong));
            int issel = -1;
            if (d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts))
                issel = 3;

            if (!(vHF.FillRecoCand(fAOD, d)))
                continue;

            if (fApplyCuts && issel > 0)
                issel = fCutsD0toKpi->IsSelected(d, AliRDHFCuts::kAll, fAOD);
            if (!d || !issel)
                continue;

            fHistNEvents->Fill(7);

            //check if primary vtx is set
            bool unsetvtx = false;
            if (!d->GetOwnPrimaryVtx())
            {
                if (!d->GetOwnPrimaryVtx())
                {
                    d->SetOwnPrimaryVtx(primVtx);
                    unsetvtx = true;
                }
            }

            //if pp recompute primary vertex without daughters
            bool isvtxrecalc = false;
            AliAODVertex *origownvtx = nullptr;
            if (fSystem == kpp)
            {
                origownvtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
                isvtxrecalc = RecalcOwnPrimaryVertex(d);
                if (!isvtxrecalc)
                {
                    CleanOwnPrimaryVertex(d, origownvtx);
                }
                unsetvtx = true;                
            }

            //fill vector of 2prongs
            if (fEnable2Prongs)
                FillCharm2Prong(d, issel);

            if (fEnableBeauty3Prongs)
            {
                for (int iTrack = 0; iTrack < fAOD->GetNumberOfTracks(); iTrack++)
                {
                    AliAODTrack *track = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(iTrack));
                    if (!track)
                        continue;

                    if (!track->TestFilterBit(4) || TMath::Abs(track->Eta()) > 0.8 || track->Pt() < 0.3) //minimal track cuts
                        continue;
                    if ((track->Charge() > 0 && issel == 1) || (track->Charge() < 0 && issel == 2)) //keep only correct mass hypothesis
                        continue;

                    //we check if the IDs of the tracks are different
                    AliAODTrack *D0dau[2] = {dynamic_cast<AliAODTrack *>(d->GetDaughter(0)), dynamic_cast<AliAODTrack *>(d->GetDaughter(1))};
                    unsigned short idD0dau[2] = {static_cast<unsigned short>(D0dau[0]->GetID()), static_cast<unsigned short>(D0dau[1]->GetID())};
                    if (track->GetID() != idD0dau[0] && track->GetID() != idD0dau[1])
                    {
                        AliExternalTrackParam piTrackParams;
                        piTrackParams.CopyFromVTrack(track);
                        AliExternalTrackParam DTrackParams;
                        DTrackParams.CopyFromVTrack(d);

                        // we calculate the vertex of the mother candidate
                        TObjArray BplusdauTracks;
                        BplusdauTracks.Add(&piTrackParams); //first the pi
                        BplusdauTracks.Add(&DTrackParams); //then the D
                        double dispersion = 0;
                        AliAODVertex *vertexBplus = ReconstructDisplVertex(fAOD->GetPrimaryVertex(), &BplusdauTracks, fAOD->GetMagneticField(), dispersion);
                        if (vertexBplus)
                        {
                            //use the new vertex to create the Bplus candidate
                            double xdummy = 0., ydummy = 0.;
                            double d0z0[2], covd0z0[3], d0[2], d0err[2];
                            piTrackParams.PropagateToDCA(vertexBplus, fAOD->GetMagneticField(), 100., d0z0, covd0z0);
                            DTrackParams.PropagateToDCA(vertexBplus, fAOD->GetMagneticField(), 100., d0z0, covd0z0);

                            //we reconstruct the mother decay prong
                            double px[2], py[2], pz[2];
                            px[0] = piTrackParams.Px();
                            py[0] = piTrackParams.Py();
                            pz[0] = piTrackParams.Pz();
                            px[1] = DTrackParams.Px();
                            py[1] = DTrackParams.Py();
                            pz[1] = DTrackParams.Pz();
                            unsigned short id[2];
                            id[0] = piTrackParams.GetID();
                            id[1] = DTrackParams.GetID();
                            piTrackParams.PropagateToDCA((AliAODVertex *)fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 100., d0z0, covd0z0);
                            d0[0] = d0z0[0];
                            d0err[0] = TMath::Sqrt(covd0z0[0]);
                            DTrackParams.PropagateToDCA((AliAODVertex *)fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 100., d0z0, covd0z0);
                            d0[1] = d0z0[0];
                            d0err[1] = TMath::Sqrt(covd0z0[0]);

                            double dca = DTrackParams.GetDCA(&piTrackParams, fAOD->GetMagneticField(), xdummy, ydummy);
                            short chargeBplus = d->Charge() + track->Charge();

                            AliAODRecoDecayHF2Prong Bplus(vertexBplus, px, py, pz, d0, d0err, dca);
                            Bplus.SetCharge(chargeBplus);
                            Bplus.GetSecondaryVtx()->AddDaughter(track); //first the pi
                            Bplus.GetSecondaryVtx()->AddDaughter(d); //then the D
                            Bplus.SetPrimaryVtxRef((AliAODVertex *)fAOD->GetPrimaryVertex());
                            Bplus.SetProngIDs(2, id);

                            FillBeauty3Prong(&Bplus, d, true);

                            delete vertexBplus;
                            vertexBplus = nullptr;
                        }
                    }
                }
            }

            if (isvtxrecalc)
            {
                CleanOwnPrimaryVertex(d, origownvtx);
                unsetvtx = true;                
            }
            if (unsetvtx)
                d->UnsetOwnPrimaryVtx();
        }
    }

    //loop on 3 prongs
    if (fEnable3Prongs || fEnableBeauty4Prongs)
    {
        for (int i3Prong = 0; i3Prong < array3Prong->GetEntriesFast(); i3Prong++)
        {
            AliAODRecoDecayHF3Prong *d = dynamic_cast<AliAODRecoDecayHF3Prong *>(array3Prong->UncheckedAt(i3Prong));

            bool isselDplus = d->HasSelectionBit(AliRDHFCuts::kDplusCuts);
            int isselDs = -1;
            if (d->HasSelectionBit(AliRDHFCuts::kDsCuts))
                isselDs = 12;
            int isselLc = -1;
            if (d->HasSelectionBit(AliRDHFCuts::kLcCuts))
                isselLc = 3;

            if (!(vHF.FillRecoCand(fAOD, d)))
                continue;

            if (fApplyCuts)
            {
                if(isselDplus)
                    isselDplus = fCutsDplustoKpipi->IsSelected(d, AliRDHFCuts::kAll, fAOD);
                if(isselDs > 0)
                    isselDs = fCutsDstoKKpi->IsSelected(d, AliRDHFCuts::kAll, fAOD);
                if(isselLc > 0)
                    isselLc = fCutsLctopKpi->IsSelected(d, AliRDHFCuts::kAll, fAOD);
            }

            if (!d || (!isselDplus && !isselDs && !isselLc))
                continue;

            if (!(fEnable3Prongs >> 0 & 1) && !(fEnableBeauty4Prongs >> 0 & 1) && (!isselDs && !isselLc))
                continue;
            if (!(fEnable3Prongs >> 1 & 1) && !(fEnableBeauty4Prongs >> 1 & 1) && (!isselDplus && !isselLc))
                continue;
            if (!(fEnable3Prongs >> 2 & 1) && !(fEnableBeauty4Prongs >> 2 & 1) && (!isselDplus && !isselDs))
                continue;

            fHistNEvents->Fill(8);

            //check if primary vtx is set
            bool unsetvtx = false;
            if (!d->GetOwnPrimaryVtx())
            {
                if (!d->GetOwnPrimaryVtx())
                {
                    d->SetOwnPrimaryVtx(primVtx);
                    unsetvtx = true;
                }
            }

            //if pp recompute primary vertex without daughters
            bool isvtxrecalc = false;
            AliAODVertex *origownvtx = nullptr;
            if (fSystem == kpp)
            {
                origownvtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
                isvtxrecalc = RecalcOwnPrimaryVertex(d);
                if (!isvtxrecalc)
                {
                    CleanOwnPrimaryVertex(d, origownvtx);
                }
                unsetvtx = true;                
            }

            //fill vector of 3prongs
            if(fEnable3Prongs)
                FillCharm3Prong(d, isselDplus, isselDs, isselLc);

            if (fEnableBeauty4Prongs)
            {
                if( (!(fEnableBeauty4Prongs >> 0 & 1) || ((fEnableBeauty4Prongs >> 0 & 1) && !isselDplus)) &&
                (!(fEnableBeauty4Prongs >> 1 & 1) || ((fEnableBeauty4Prongs >> 1 & 1) && (!(isselDs & 4) && !(isselDs & 8))) ) && 
                (!(fEnableBeauty4Prongs >> 2 & 1) || ((fEnableBeauty4Prongs >> 2 & 1) && !isselLc)))
                {
                    if (isvtxrecalc)
                    {
                        CleanOwnPrimaryVertex(d, origownvtx);
                        unsetvtx = true;                
                    }
                    if (unsetvtx)
                        d->UnsetOwnPrimaryVtx();
                    continue;
                }

                for (int iTrack = 0; iTrack < fAOD->GetNumberOfTracks(); iTrack++)
                {

                    bool isselB0 = isselDplus;
                    bool isselBs = static_cast<bool>(isselDs);
                    bool isselLb = static_cast<bool>(isselLc);

                    AliAODTrack *track = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(iTrack));
                    if (!track)
                        continue;

                    if (!track->TestFilterBit(4) || TMath::Abs(track->Eta()) > 0.8 || track->Pt() < 0.3) //minimal track cuts
                        continue;

                    //keep only correct charge
                    if((track->Charge() > 0 && d->Charge() > 0) || (track->Charge() < 0 && d->Charge() < 0))
                        continue;

                    //we check if the IDs of the tracks are different
                    AliAODTrack *Ddau[3] = {dynamic_cast<AliAODTrack *>(d->GetDaughter(0)), dynamic_cast<AliAODTrack *>(d->GetDaughter(1)), dynamic_cast<AliAODTrack *>(d->GetDaughter(2))};
                    unsigned short idDdau[3] = {static_cast<unsigned short>(Ddau[0]->GetID()), static_cast<unsigned short>(Ddau[1]->GetID()), static_cast<unsigned short>(Ddau[2]->GetID())};
                    if (track->GetID() != idDdau[0] && track->GetID() != idDdau[1] && track->GetID() != idDdau[2])
                    {
                        AliExternalTrackParam piTrackParams;
                        piTrackParams.CopyFromVTrack(track);
                        AliExternalTrackParam DTrackParams;
                        DTrackParams.CopyFromVTrack(d);

                        // we calculate the vertex of the mother candidate
                        TObjArray BdauTracks;
                        BdauTracks.Add(&DTrackParams); //first the D
                        BdauTracks.Add(&piTrackParams); // then the pi
                        double dispersion = 0;
                        AliAODVertex *vertexB = ReconstructDisplVertex(fAOD->GetPrimaryVertex(), &BdauTracks, fAOD->GetMagneticField(), dispersion);
                        if (vertexB)
                        {
                            //use the new vertex to create the B candidate
                            double xdummy = 0., ydummy = 0.;
                            double d0z0[2], covd0z0[3], d0[2], d0err[2];
                            piTrackParams.PropagateToDCA(vertexB, fAOD->GetMagneticField(), 100., d0z0, covd0z0);
                            DTrackParams.PropagateToDCA(vertexB, fAOD->GetMagneticField(), 100., d0z0, covd0z0);

                            //we reconstruct the mother decay prong
                            double px[2], py[2], pz[2];
                            px[0] = piTrackParams.Px();
                            py[0] = piTrackParams.Py();
                            pz[0] = piTrackParams.Pz();
                            px[1] = DTrackParams.Px();
                            py[1] = DTrackParams.Py();
                            pz[1] = DTrackParams.Pz();
                            unsigned short id[2];
                            id[0] = piTrackParams.GetID();
                            id[1] = DTrackParams.GetID();
                            piTrackParams.PropagateToDCA((AliAODVertex *)fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 100., d0z0, covd0z0);
                            d0[0] = d0z0[0];
                            d0err[0] = TMath::Sqrt(covd0z0[0]);
                            DTrackParams.PropagateToDCA((AliAODVertex *)fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 100., d0z0, covd0z0);
                            d0[1] = d0z0[0];
                            d0err[1] = TMath::Sqrt(covd0z0[0]);

                            double dca = DTrackParams.GetDCA(&piTrackParams, fAOD->GetMagneticField(), xdummy, ydummy);
                            short chargeB = d->Charge() + track->Charge();

                            AliAODRecoDecayHF2Prong B(vertexB, px, py, pz, d0, d0err, dca);
                            B.SetCharge(chargeB);
                            B.GetSecondaryVtx()->AddDaughter(d); //first the D
                            B.GetSecondaryVtx()->AddDaughter(track); //then the pi
                            B.SetPrimaryVtxRef((AliAODVertex *)fAOD->GetPrimaryVertex());
                            B.SetProngIDs(2, id);

                            FillBeauty4Prong(&B, d, isselB0, isselBs, isselLb, isselDs, isselLc);

                            delete vertexB;
                            vertexB = nullptr;
                        }
                    }
                }
            }

            if (isvtxrecalc)
            {
                CleanOwnPrimaryVertex(d, origownvtx);
                unsetvtx = true;                
            }
            if (unsetvtx)
                d->UnsetOwnPrimaryVtx();
        }
    }

    //loop on Dstars
    if (fEnableDstars)
    {
        for (int iDstar = 0; iDstar < arrayDstar->GetEntriesFast(); iDstar++)
        {
            AliAODRecoCascadeHF *d = dynamic_cast<AliAODRecoCascadeHF *>(arrayDstar->UncheckedAt(iDstar));
            bool issel = d->HasSelectionBit(AliRDHFCuts::kDstarCuts);

            if (!(vHF.FillRecoCand(fAOD, d)))
                continue;

            if (fApplyCuts && issel)
                issel = fCutsDstartoKpipi->IsSelected(d, AliRDHFCuts::kAll, fAOD);

            if (!d || !issel)
                continue;

            AliAODRecoDecayHF2Prong *d0 = d->Get2Prong();
            if (!d0)
                continue;

            fHistNEvents->Fill(9);

            //check if primary vtx is set
            bool unsetvtx = false;
            if (!d->GetOwnPrimaryVtx())
            {
                if (!d->GetOwnPrimaryVtx())
                {
                    d->SetOwnPrimaryVtx(primVtx);
                    unsetvtx = true;
                }
            }

            //if pp recompute primary vertex without daughters
            bool isvtxrecalc = false;
            AliAODVertex *origownvtx = nullptr;
            if (fSystem == kpp)
            {
                origownvtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
                isvtxrecalc = RecalcOwnPrimaryVertex(d);
                if (!isvtxrecalc)
                {
                    CleanOwnPrimaryVertex(d, origownvtx);
                }
                unsetvtx = true;                
            }

            //fill vector of dstars
            FillDstar(d, d0, issel);

            if (isvtxrecalc)
            {
                CleanOwnPrimaryVertex(d, origownvtx);
                unsetvtx = true;                
            }
            if (unsetvtx)
                d->UnsetOwnPrimaryVtx();
        }
    }

    //loop on cascades
    if (fEnableCascades)
    {
        for (int iCasc = 0; iCasc < arrayCasc->GetEntriesFast(); iCasc++)
        {
            AliAODRecoCascadeHF *lc = dynamic_cast<AliAODRecoCascadeHF *>(arrayCasc->UncheckedAt(iCasc));
            if (!lc)
                continue;
            AliAODv0 *v0part = lc->Getv0();
            if (!v0part)
                continue;

            int issel = -1;
            if(lc->HasSelectionBit(AliRDHFCuts::kLctoV0Cuts))
                issel = 3;

            if (!(vHF.FillRecoCand(fAOD, lc)))
                continue;

            if (fApplyCuts && issel > 0)
                issel = fCutsLctoV0bach->IsSelected(lc, AliRDHFCuts::kAll, fAOD);

            if (issel)
                continue;

            fHistNEvents->Fill(10);

            //check if primary vtx is set
            bool unsetvtx = false;
            if (!lc->GetOwnPrimaryVtx())
            {
                if (!lc->GetOwnPrimaryVtx())
                {
                    lc->SetOwnPrimaryVtx(primVtx);
                    unsetvtx = true;
                }
            }

            //if pp recompute primary vertex without daughters
            bool isvtxrecalc = false;
            AliAODVertex *origownvtx = nullptr;
            if (fSystem == kpp)
            {
                origownvtx = new AliAODVertex(*lc->GetOwnPrimaryVtx());
                isvtxrecalc = RecalcOwnPrimaryVertex(lc);
                if (!isvtxrecalc)
                {
                    CleanOwnPrimaryVertex(lc, origownvtx);
                }
                unsetvtx = true;                
            }

            //fill vector of cascades
            FillCharmCascade(lc, v0part, issel);

            if (isvtxrecalc)
            {
                CleanOwnPrimaryVertex(lc, origownvtx);
                unsetvtx = true;                
            }
            if (unsetvtx)
                lc->UnsetOwnPrimaryVtx();
        }
    }

    fRecoTree->Fill();
    if(fReadMC && fFillGenTree)
        fGenTree->Fill();

    fCharm2Prong.clear();
    fCharm3Prong.clear();
    fDstar.clear();
    fCharmCascade.clear();
    fBeauty3Prong.clear();
    fBeauty4Prong.clear();
    fGenHadron.clear();

    PostData(1, fOutput);
    PostData(2, fRecoTree);
    PostData(3, fGenTree);

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::FillCharm2Prong(AliAODRecoDecayHF2Prong *cand, int issel)
{
    Charm2Prong ch2Prong;
    ch2Prong.fPt = cand->Pt();
    ch2Prong.fY = cand->Y(421);
    ch2Prong.fInvMassD0 = cand->InvMassD0();
    ch2Prong.fInvMassD0bar = cand->InvMassD0bar();
    ch2Prong.fCosP = cand->CosPointingAngle();
    ch2Prong.fCosPXY = cand->CosPointingAngleXY();
    ch2Prong.fDecayLength = cand->DecayLength();
    ch2Prong.fNormDecayLengthXY = cand->NormalizedDecayLengthXY();
    ch2Prong.fImpParProd = cand->Getd0Prong(0) * cand->Getd0Prong(1);
    double ptDau[2] = {dynamic_cast<AliAODTrack*>(cand->GetDaughter(0))->Pt(), dynamic_cast<AliAODTrack*>(cand->GetDaughter(1))->Pt()};
    double absd0Dau[2] = {TMath::Abs(cand->Getd0Prong(0)), TMath::Abs(cand->Getd0Prong(1))};
    ch2Prong.fd0MinDau = *std::min_element(absd0Dau, absd0Dau+2);
    ch2Prong.fPtMinDau = *std::min_element(ptDau, ptDau+2);
    ch2Prong.fProngIdx0 = dynamic_cast<AliAODTrack*>(cand->GetDaughter(0))->GetID();
    ch2Prong.fProngIdx1 = dynamic_cast<AliAODTrack*>(cand->GetDaughter(1))->GetID();

    ch2Prong.fSelBit = 0;
    if (issel == 1 || issel == 3)
        ch2Prong.fSelBit |= kDzerotoKpiCuts;
    if (issel >= 2)
        ch2Prong.fSelBit |= kDzerotopiKCuts;
    if (fCutsD0toKpi->IsSelectedPID(cand) == 1 || fCutsD0toKpi->IsSelectedPID(cand) == 3)
        ch2Prong.fSelBit |= kDzerotoKpiCutsPID;
    if (fCutsD0toKpi->IsSelectedPID(cand) >= 2)
        ch2Prong.fSelBit |= kDzerotopiKCutsPID;
    if (fCutsD0toKpi->IsInFiducialAcceptance(ch2Prong.fPt, ch2Prong.fY))
        ch2Prong.fSelBit |= kDzerotoKpiFidAcc;

    ch2Prong.fCandType = 0;
    ch2Prong.fGenLabel = -1;
    ch2Prong.fDecay = kNone;

    if(fReadMC)
    {
        int pdgDgD0toKpi[2] = {321, 211};
        ch2Prong.fGenLabel = cand->MatchToMC(421, fMCArray, 2, pdgDgD0toKpi);
        ch2Prong.fDecay = kNone;
        int origin = -1;
        if (ch2Prong.fGenLabel < 0)
        {
            ch2Prong.fCandType |= kBackground;
            if (fFillOnlySignal)
                return;
        }
        else
        {
            AliAODMCParticle *partD0 = dynamic_cast<AliAODMCParticle *>(fMCArray->At(ch2Prong.fGenLabel));
            origin = AliVertexingHFUtils::CheckOrigin(fMCArray, partD0, true);
            if (origin == 4)
            {
                ch2Prong.fCandType |= kSignal;
                ch2Prong.fCandType |= kPrompt;
            }
            else if (origin == 5)
            {
                ch2Prong.fCandType |= kSignal;
                ch2Prong.fCandType |= kFeedDown;
            }
            else
            {
                ch2Prong.fCandType |= kSignal; // no prompt, no feed-down --> weird stuff
            }

            int labDau0 = dynamic_cast<AliAODTrack *>(cand->GetDaughter(0))->GetLabel();
            AliAODMCParticle *dauPart0 = dynamic_cast<AliAODMCParticle *>(fMCArray->UncheckedAt(TMath::Abs(labDau0)));
            int pdgCode0 = TMath::Abs(dauPart0->GetPdgCode());
            if (pdgCode0 == 321)
                ch2Prong.fDecay = kDzerotoKpi;
            else if (pdgCode0 == 211)
                ch2Prong.fDecay = kDzerotopiK;
        }
    }

    fCharm2Prong.push_back(ch2Prong);
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::FillCharm3Prong(AliAODRecoDecayHF3Prong *cand, bool isselDplus, int isselDs, int isselLc)
{
    Charm3Prong ch3Prong;
    ch3Prong.fPt = cand->Pt();
    ch3Prong.fInvMassDplus = cand->InvMassDplus();
    ch3Prong.fInvMassDstoKKpi = cand->InvMassDsKKpi();
    ch3Prong.fInvMassDstopiKK = cand->InvMassDspiKK();
    ch3Prong.fInvMassLctopKpi = cand->InvMassLcpKpi();
    ch3Prong.fInvMassLctopiKp = cand->InvMassLcpiKp();
    ch3Prong.fYDplus = cand->YDplus();
    ch3Prong.fYDs = cand->YDs();
    ch3Prong.fYLc = cand->YLc();
    ch3Prong.fInvMassPhiKKpi = cand->InvMass2Prongs(0, 1, 321, 321);
    ch3Prong.fInvMassPhipiKK = cand->InvMass2Prongs(1, 2, 321, 321);
    ch3Prong.fCosP = cand->CosPointingAngle();
    ch3Prong.fCosPXY = cand->CosPointingAngleXY();
    ch3Prong.fDecayLength = cand->DecayLength();
    ch3Prong.fNormDecayLengthXY = cand->NormalizedDecayLengthXY();
    ch3Prong.fSigmaVtx = cand->GetSigmaVert();
    double ptDau[3] = {dynamic_cast<AliAODTrack*>(cand->GetDaughter(0))->Pt(), dynamic_cast<AliAODTrack*>(cand->GetDaughter(1))->Pt(), dynamic_cast<AliAODTrack*>(cand->GetDaughter(2))->Pt()};
    double absd0Dau[3] = {TMath::Abs(cand->Getd0Prong(0)), TMath::Abs(cand->Getd0Prong(1)), TMath::Abs(cand->Getd0Prong(2))};
    ch3Prong.fd0MinDau = *std::min_element(absd0Dau, absd0Dau+3);
    ch3Prong.fPtMinDau = *std::min_element(ptDau, ptDau+3);
    ch3Prong.fProngIdx0 = dynamic_cast<AliAODTrack*>(cand->GetDaughter(0))->GetID();
    ch3Prong.fProngIdx1 = dynamic_cast<AliAODTrack*>(cand->GetDaughter(1))->GetID();
    ch3Prong.fProngIdx2 = dynamic_cast<AliAODTrack*>(cand->GetDaughter(2))->GetID();

    ch3Prong.fSelBit = 0;
    if (isselDplus)
    {
        ch3Prong.fSelBit |= kDplustoKpipiCuts;
        if (fCutsDplustoKpipi->IsSelectedPID(cand))
            ch3Prong.fSelBit |= kDplustoKpipiCutsPID;
        if (fCutsDplustoKpipi->IsInFiducialAcceptance(ch3Prong.fPt, ch3Prong.fYDplus))
            ch3Prong.fSelBit |= kDplustoKpipiFidAcc;
    }
    if (isselDs)
    {
        if (isselDs & 4)
            ch3Prong.fSelBit |= kDstoKKpiCuts;
        if (isselDs & 8)
            ch3Prong.fSelBit |= kDstopiKKCuts;

        if (fCutsDstoKKpi->IsSelectedPID(cand) == 1 || fCutsDstoKKpi->IsSelectedPID(cand) == 3)
            ch3Prong.fSelBit |= kDstoKKpiCutsPID;
        if (fCutsDstoKKpi->IsSelectedPID(cand) >= 2)
            ch3Prong.fSelBit |= kDstopiKKCutsPID;
        if (fCutsDstoKKpi->IsInFiducialAcceptance(ch3Prong.fPt, ch3Prong.fYDs))
            ch3Prong.fSelBit |= kDstoKKpiFidAcc;
    }
    if (isselLc)
    {
        if (isselLc == 1 || isselLc == 3)
            ch3Prong.fSelBit |= kLctopKpiCuts;
        if (isselLc >= 2)
            ch3Prong.fSelBit |= kLctopiKpCuts;

        if (fCutsLctopKpi->IsSelectedPID(cand) == 1 || fCutsLctopKpi->IsSelectedPID(cand) == 1)
            ch3Prong.fSelBit |= kLctopKpiCutsPID;
        if (fCutsLctopKpi->IsSelectedPID(cand) >= 2)
            ch3Prong.fSelBit |= kLctopiKpCutsPID;
        if (fCutsLctopKpi->IsInFiducialAcceptance(ch3Prong.fPt, ch3Prong.fYLc))
            ch3Prong.fSelBit |= kLctopKpiFidAcc;
    }

    ch3Prong.fCandType = 0;
    ch3Prong.fGenLabel = -1;
    ch3Prong.fDecay = kNone;

    if(fReadMC)
    {
        int pdgDgDplustoKpipi[3] = {321, 211, 211};
        int pdgDgDstoKKpi[3] = {321, 321, 211};
        int pdgDgLctopKpi[3] = {2212, 321, 211};
        int origin = -1;

        int labDplus = cand->MatchToMC(411, fMCArray, 3, pdgDgDplustoKpipi);
        int labDs = -1;
        int labDplustoKKpi = -1;
        int labLc = -1;
        if (labDplus < 0)
        {
            labDs = cand->MatchToMC(431, fMCArray, 3, pdgDgDstoKKpi);
            if (labDs < 0)
            {
                labDplustoKKpi = cand->MatchToMC(411, fMCArray, 3, pdgDgDstoKKpi);
                if (labDplustoKKpi < 0)
                {
                    labLc = cand->MatchToMC(4122, fMCArray, 3, pdgDgLctopKpi);
                    if (labLc >= 0)
                    {
                        ch3Prong.fGenLabel = labLc;

                        int labDau0 = dynamic_cast<AliAODTrack *>(cand->GetDaughter(0))->GetLabel();
                        AliAODMCParticle *dauPart0 = dynamic_cast<AliAODMCParticle *>(fMCArray->UncheckedAt(TMath::Abs(labDau0)));
                        int pdgCode0 = TMath::Abs(dauPart0->GetPdgCode());
                        if (pdgCode0 == 2212)
                            ch3Prong.fDecay = kLctopKpi;
                        else if (pdgCode0 == 211)
                            ch3Prong.fDecay = kLctopiKp;
                    }
                }
                else
                {
                    ch3Prong.fGenLabel = labDplustoKKpi;
                    int labDau0 = dynamic_cast<AliAODTrack *>(cand->GetDaughter(0))->GetLabel();
                    AliAODMCParticle *dauPart0 = dynamic_cast<AliAODMCParticle *>(fMCArray->UncheckedAt(TMath::Abs(labDau0)));
                    int pdgCode0 = TMath::Abs(dauPart0->GetPdgCode());
                    if (pdgCode0 == 321)
                        ch3Prong.fDecay = kDplustoKKpi;
                    else if (pdgCode0 == 211)
                        ch3Prong.fDecay = kDplustopiKK;
                }
            }
            else
            {
                ch3Prong.fGenLabel = labDs;
                int labDau0 = dynamic_cast<AliAODTrack *>(cand->GetDaughter(0))->GetLabel();
                AliAODMCParticle *dauPart0 = dynamic_cast<AliAODMCParticle *>(fMCArray->UncheckedAt(TMath::Abs(labDau0)));
                int pdgCode0 = TMath::Abs(dauPart0->GetPdgCode());
                if (pdgCode0 == 321)
                    ch3Prong.fDecay = kDstoKKpi;
                else if (pdgCode0 == 211)
                    ch3Prong.fDecay = kDstopiKK;
            }
        }
        else
        {
            ch3Prong.fGenLabel = labDplus;
            ch3Prong.fDecay = kDplustoKpipi;
        }

        if (ch3Prong.fGenLabel < 0)
        {
            ch3Prong.fCandType |= kBackground;
            if (fFillOnlySignal)
                return;
        }
        else
        {
            AliAODMCParticle *part3prong = dynamic_cast<AliAODMCParticle *>(fMCArray->At(ch3Prong.fGenLabel));
            origin = AliVertexingHFUtils::CheckOrigin(fMCArray, part3prong, true);
            if (origin == 4)
            {
                ch3Prong.fCandType |= kSignal;
                ch3Prong.fCandType |= kPrompt;
            }
            else if (origin == 5)
            {
                ch3Prong.fCandType |= kSignal;
                ch3Prong.fCandType |= kFeedDown;
            }
            else
            {
                ch3Prong.fCandType |= kSignal; // no prompt, no feed-down --> weird stuff
            }
        }
    }

    fCharm3Prong.push_back(ch3Prong);
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::FillDstar(AliAODRecoCascadeHF *cand, AliAODRecoDecayHF2Prong *dau, bool issel)
{
    Dstar dstar;
    dstar.fInvMass = cand->InvMassDstarKpipi();
    dstar.fInvMassD0 = cand->InvMassD0();
    dstar.fPt = cand->Pt();
    dstar.fY = cand->Y(413);
    dstar.fCosPD0 = dau->CosPointingAngle();
    dstar.fCosPXYD0 = dau->CosPointingAngleXY();
    dstar.fDecayLengthD0 = dau->DecayLength();
    dstar.fNormDecayLengthXYD0 = dau->NormalizedDecayLengthXY();
    double ptDau[3] = {dynamic_cast<AliAODTrack*>(cand->GetBachelor())->Pt(), dynamic_cast<AliAODTrack*>(dau->GetDaughter(0))->Pt(), dynamic_cast<AliAODTrack*>(dau->GetDaughter(1))->Pt()};
    double absd0Dau[3] = {TMath::Abs(cand->Getd0Prong(0)), TMath::Abs(dau->Getd0Prong(0)), TMath::Abs(dau->Getd0Prong(1))};
    dstar.fd0MinDau = *std::min_element(absd0Dau, absd0Dau+3);
    dstar.fPtMinDau = *std::min_element(ptDau, ptDau+3);
    dstar.fProngIdx0 = dynamic_cast<AliAODTrack*>(cand->GetBachelor())->GetID();
    dstar.fProngIdx1 = dynamic_cast<AliAODTrack*>(dau->GetDaughter(0))->GetID();
    dstar.fProngIdx2 = dynamic_cast<AliAODTrack*>(dau->GetDaughter(1))->GetID();

    dstar.fSelBit = 0;
    if (issel)
        dstar.fSelBit |= kDstartoKpipiCuts;
    if (fCutsDstartoKpipi->IsSelectedPID(cand))
        dstar.fSelBit |= kDstartoKpipiCutsPID;
    if (fCutsDstartoKpipi->IsInFiducialAcceptance(dstar.fPt, dstar.fY))
        dstar.fSelBit |= kDstartoKpipiFidAcc;

    dstar.fDecay = kNone;
    dstar.fCandType = 0;
    dstar.fGenLabel = -1;

    if(fReadMC)
    {
        int pdgDgDstartoD0pi[3] = {421, 211};
        int pdgDgD0toKpi[2] = {321, 211};
        dstar.fGenLabel = cand->MatchToMC(413, 421, pdgDgDstartoD0pi, pdgDgD0toKpi, fMCArray);
        int origin = -1;
        if (dstar.fGenLabel < 0)
        {
            dstar.fCandType |= kBackground;
            if (fFillOnlySignal)
                return;
        }
        else
        {
            AliAODMCParticle *partDstar = dynamic_cast<AliAODMCParticle *>(fMCArray->At(dstar.fGenLabel));
            origin = AliVertexingHFUtils::CheckOrigin(fMCArray, partDstar, true);
            if (origin == 4)
            {
                dstar.fCandType |= kSignal;
                dstar.fCandType |= kPrompt;
            }
            else if (origin == 5)
            {
                dstar.fCandType |= kSignal;
                dstar.fCandType |= kFeedDown;
            }
            else
            {
                dstar.fCandType |= kSignal; // no prompt, no feed-down --> weird stuff
            }
            dstar.fDecay = kDstartoKpipi;
        }
    }

    fDstar.push_back(dstar);
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::FillCharmCascade(AliAODRecoCascadeHF *cand, AliAODv0 *dau, int issel)
{
    CharmCascade chCasc;

    chCasc.fInvMassLctopK0s = cand->InvMassLctoK0sP();
    chCasc.fInvMassLctopiLambda = cand->InvMassLctoLambdaPi();
    chCasc.fInvMassK0s = dau->MassK0Short();
    if (cand->Charge() > 0)
        chCasc.fInvMassLambda = dau->MassLambda();
    else
        chCasc.fInvMassLambda = dau->MassAntiLambda();
    chCasc.fPt = cand->Pt();
    chCasc.fY = cand->Y(4122);
    chCasc.fCosPV0 = cand->CosV0PointingAngle();
    chCasc.fCosPXYV0 = cand->CosV0PointingAngleXY();
    double ptDau[3] = {dynamic_cast<AliAODTrack*>(cand->GetBachelor())->Pt(), dynamic_cast<AliAODTrack*>(dau->GetDaughter(0))->Pt(), dynamic_cast<AliAODTrack*>(dau->GetDaughter(1))->Pt()};
    double absd0Dau[3] = {TMath::Abs(cand->Getd0Prong(0)), TMath::Abs(dau->Getd0Prong(0)), TMath::Abs(dau->Getd0Prong(1))};
    chCasc.fd0MinDau = *std::min_element(absd0Dau, absd0Dau+3);
    chCasc.fPtMinDau = *std::min_element(ptDau, ptDau+3);
    chCasc.fProngIdx0 = dynamic_cast<AliAODTrack*>(cand->GetBachelor())->GetID();
    chCasc.fProngIdx1 = dynamic_cast<AliAODTrack*>(dau->GetDaughter(0))->GetID();
    chCasc.fProngIdx2 = dynamic_cast<AliAODTrack*>(dau->GetDaughter(1))->GetID();

    chCasc.fSelBit = 0;
    if (issel)
        chCasc.fSelBit |= kLctoV0bachCuts;
    if (fCutsDstartoKpipi->IsSelectedPID(cand))
        chCasc.fSelBit |= kLctoV0bachCutsPID;
    if (fCutsDstartoKpipi->IsInFiducialAcceptance(chCasc.fPt, chCasc.fY))
        chCasc.fSelBit |= kLctoV0bachFidAcc;


    chCasc.fGenLabel = -1;
    chCasc.fDecay = kNone;
    chCasc.fCandType = 0;

    if(fReadMC)
    {
        int pdgDgLctopK0s[2] = {2212, 310};
        int pdgDgLctopiLambda[2] = {211, 3122};
        int pdgDgK0s[2] = {211, 211};
        int pdgDgLambda[2] = {2212, 211};

        int labtoK0s = cand->MatchToMC(4122, 310, pdgDgLctopK0s, pdgDgK0s, fMCArray, true);
        int labtoLambda = -1;
        if (labtoK0s < 0)
        {
            labtoLambda = cand->MatchToMC(4122, 3122, pdgDgLctopiLambda, pdgDgLambda, fMCArray, true);
            if (labtoLambda >= 0)
            {
                chCasc.fGenLabel = labtoLambda;
                chCasc.fDecay = kLctopK0s;
            }
        }
        else
        {
            chCasc.fGenLabel = labtoK0s;
            chCasc.fDecay = kLctopiLambda;
        }

        int origin = -1;
        if (chCasc.fGenLabel < 0)
        {
            chCasc.fCandType |= kBackground;
            if (fFillOnlySignal)
                return;
        }
        else
        {
            AliAODMCParticle *partCasc = dynamic_cast<AliAODMCParticle *>(fMCArray->At(chCasc.fGenLabel));
            origin = AliVertexingHFUtils::CheckOrigin(fMCArray, partCasc, true);
            if (origin == 4)
            {
                chCasc.fCandType |= kSignal;
                chCasc.fCandType |= kPrompt;
            }
            else if (origin == 5)
            {
                chCasc.fCandType |= kSignal;
                chCasc.fCandType |= kFeedDown;
            }
            else
            {
                chCasc.fCandType |= kSignal; // no prompt, no feed-down --> weird stuff
            }
        }
    }

    fCharmCascade.push_back(chCasc);
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::FillBeauty3Prong(AliAODRecoDecayHF2Prong *cand, AliAODRecoDecayHF2Prong *dau, bool issel)
{
    Beauty3Prong b3Prong;
    unsigned int pdgDgBplustoD0pi[2] = {211, 421};
    int pdgDgD0topiK[2] = {211, 321};
    double invmassBplus = cand->InvMass(2, pdgDgBplustoD0pi);
    double massBplusPDG = TDatabasePDG::Instance()->GetParticle(521)->Mass();
    if (TMath::Abs(massBplusPDG - invmassBplus) > 0.4) //check mass
        return;

    b3Prong.fInvMassBplustoD0pi = invmassBplus;
    b3Prong.fPt = cand->Pt();
    b3Prong.fY = cand->Y(521);
    b3Prong.fDecayLength = cand->DecayLength();
    b3Prong.fNormDecayLengthXY = cand->NormalizedDecayLengthXY();
    b3Prong.fCosP = cand->CosPointingAngle();
    b3Prong.fCosPXY = cand->CosPointingAngleXY();
    b3Prong.fImpParProd = cand->Getd0Prong(0) * cand->Getd0Prong(1);
    b3Prong.fPtD0 = dau->Pt();
    b3Prong.fCosPD0 = dau->CosPointingAngle();
    b3Prong.fCosPXYD0 = dau->CosPointingAngleXY();
    b3Prong.fDecayLengthD0 = dau->DecayLength();
    b3Prong.fNormDecayLengthXYD0 = dau->NormalizedDecayLengthXY();
    b3Prong.fImpParProdD0 = dau->Getd0Prong(0) * dau->Getd0Prong(1);
    if (cand->Charge() > 0)
        b3Prong.fInvMassD0 = dau->InvMassD0bar();
    else
        b3Prong.fInvMassD0 = dau->InvMassD0();

    double ptDau[2] = {dynamic_cast<AliAODTrack*>(cand->GetDaughter(0))->Pt(), dynamic_cast<AliAODTrack*>(cand->GetDaughter(1))->Pt()};
    double absd0Dau[2] = {TMath::Abs(dau->Getd0Prong(0)), TMath::Abs(dau->Getd0Prong(1))};
    b3Prong.fPtMinDauD0 = *std::min_element(ptDau, ptDau+2);
    b3Prong.fd0MinDauD0 = *std::min_element(absd0Dau, absd0Dau+2);

    b3Prong.fSelBit = 0;
    if (issel)
        b3Prong.fSelBit |= kBplustoD0piCuts;

    b3Prong.fDecay = kNone;
    b3Prong.fCandType = 0;
    b3Prong.fGenLabel = -1;

    if(fReadMC)
    {
        b3Prong.fGenLabel = cand->MatchToMCB2Prong(521, 421, reinterpret_cast<int *>(pdgDgBplustoD0pi), pdgDgD0topiK, fMCArray);
        if (b3Prong.fGenLabel < 0)
        {
            b3Prong.fCandType |= kBackground;
            if (fFillOnlySignal)
                return;
        }
        else
        {
            b3Prong.fCandType |= kSignal;
            b3Prong.fCandType |= kPrompt; //beauty always prompt
            b3Prong.fDecay = kBplustoD0pi;
        }

        if (IsInFiducialAcceptance(b3Prong.fPt, b3Prong.fY))
            b3Prong.fSelBit |= kBplustoD0piFidAcc;
    }

    fBeauty3Prong.push_back(b3Prong);
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::FillBeauty4Prong(AliAODRecoDecayHF2Prong *cand, AliAODRecoDecayHF3Prong *dau, bool isselB0, bool isselBs, bool isselLb, int isselDs, int isselLc)
{
    Beauty4Prong b4Prong;
    unsigned int pdgDgB0toDminuspi[2] = {211, 411};
    int pdgDgDplustopipiK[3] = {211, 211, 321};

    unsigned int pdgDgBstoDsminuspi[2] = {211, 431};
    int pdgDgDstopiKK[3] = {211, 321, 321};

    unsigned int pdgDgLbtoLcpi[2] = {211, 4122};
    int pdgDgLctopKpi[3] = {2212, 321, 211};

    double massB0PDG = TDatabasePDG::Instance()->GetParticle(511)->Mass();
    double massBsPDG = TDatabasePDG::Instance()->GetParticle(531)->Mass();
    double massLbPDG = TDatabasePDG::Instance()->GetParticle(5122)->Mass();

    if(isselB0)
    {
        double invmassB0 = cand->InvMass(2, pdgDgB0toDminuspi);
        b4Prong.fInvMassB0toDminuspi = invmassB0;
        if (TMath::Abs(massB0PDG - invmassB0) > 0.4) //check mass
            isselB0 = false;
    }
    if(isselBs)
    {
        double invmassBs = cand->InvMass(2, pdgDgBstoDsminuspi);
        b4Prong.fInvMassBstoDsminuspi = invmassBs;
        if (TMath::Abs(massBsPDG - invmassBs) > 0.4) //check mass
            isselBs = false;
    }
    if(isselLb)
    {
        double invmassLb = cand->InvMass(2, pdgDgLbtoLcpi);
        b4Prong.fInvMassLbtoLcpluspi = invmassLb;
        if (TMath::Abs(massLbPDG - invmassLb) > 0.4) //check mass
            isselLb = false;
    }

    b4Prong.fPt = cand->Pt();
    b4Prong.fYB0 = cand->Y(511);
    b4Prong.fYBs = cand->Y(531);
    b4Prong.fYLb = cand->Y(5122);
    b4Prong.fDecayLength = cand->DecayLength();
    b4Prong.fNormDecayLengthXY = cand->NormalizedDecayLengthXY();
    b4Prong.fCosP = cand->CosPointingAngle();
    b4Prong.fCosPXY = cand->CosPointingAngleXY();
    b4Prong.fImpParProd = cand->Getd0Prong(0) * cand->Getd0Prong(1);
    b4Prong.fPtD = dau->Pt();
    b4Prong.fCosPD = dau->CosPointingAngle();
    b4Prong.fCosPXYD = dau->CosPointingAngleXY();
    b4Prong.fDecayLengthD = dau->DecayLength();
    b4Prong.fNormDecayLengthXYD = dau->NormalizedDecayLengthXY();
    b4Prong.fSigmaVtxD = dau->GetSigmaVert();

    double ptDau[3] = {dynamic_cast<AliAODTrack*>(dau->GetDaughter(0))->Pt(), dynamic_cast<AliAODTrack*>(dau->GetDaughter(1))->Pt(), dynamic_cast<AliAODTrack*>(dau->GetDaughter(2))->Pt()};
    double absd0Dau[3] = {TMath::Abs(dau->Getd0Prong(0)), TMath::Abs(dau->Getd0Prong(1)), TMath::Abs(dau->Getd0Prong(2))};
    b4Prong.fd0MinDauD = *std::min_element(absd0Dau, absd0Dau+3);
    b4Prong.fPtMinDauD = *std::min_element(ptDau, ptDau+3);

    b4Prong.fInvMassDplus = dau->InvMassDplus();
    b4Prong.fInvMassDs = -1.;
    b4Prong.fInvMassLc = -1.;
    if(isselDs & 4)
        b4Prong.fInvMassDs = dau->InvMassDsKKpi();
    else if(isselDs & 8)
        b4Prong.fInvMassDs = dau->InvMassDspiKK();
    if(isselLc == 1 || isselLc == 3)
        b4Prong.fInvMassLc = dau->InvMassLcpKpi();
    else if(isselLc == 2)
        b4Prong.fInvMassLc = dau->InvMassLcpKpi();

    b4Prong.fSelBit = 0;
    if (isselB0)
        b4Prong.fSelBit |= kB0toDminuspiCuts;
    if (isselBs)
        b4Prong.fSelBit |= kBstoDminuspiCuts;
    if (isselLb)
        b4Prong.fSelBit |= kLbtoLcpluspiCuts;

    b4Prong.fDecay = kNone;
    b4Prong.fCandType = 0;
    b4Prong.fGenLabel = -1;

    if(fReadMC)
    {
        b4Prong.fGenLabel = cand->MatchToMCB3Prong(511, 411, reinterpret_cast<int *>(pdgDgB0toDminuspi), pdgDgDplustopipiK, fMCArray);
        if (b4Prong.fGenLabel < 0)
        {
            b4Prong.fGenLabel = cand->MatchToMCB3Prong(531, 431, reinterpret_cast<int *>(pdgDgBstoDsminuspi), pdgDgDstopiKK, fMCArray);
            if (b4Prong.fGenLabel < 0)
            {
                b4Prong.fGenLabel = cand->MatchToMCB3Prong(5122, 4122, reinterpret_cast<int *>(pdgDgLbtoLcpi), pdgDgLctopKpi, fMCArray);
                if (b4Prong.fGenLabel >= 0)
                {
                    b4Prong.fCandType |= kSignal;
                    b4Prong.fCandType |= kPrompt; //beauty always prompt
                    b4Prong.fDecay = kLbtoLcpluspi;
                }
            }
            else
            {
                b4Prong.fCandType |= kSignal;
                b4Prong.fCandType |= kPrompt; //beauty always prompt
                b4Prong.fDecay = kBstoDsminuspi;
            }
        }
        else
        {
            b4Prong.fCandType |= kSignal;
            b4Prong.fCandType |= kPrompt; //beauty always prompt
            b4Prong.fDecay = kB0toDminuspi;
        }

        if (b4Prong.fGenLabel < 0)
        {
            b4Prong.fCandType |= kBackground;
            if (fFillOnlySignal)
                return;
        }

        if (IsInFiducialAcceptance(b4Prong.fPt, b4Prong.fYB0))
            b4Prong.fSelBit |= kB0toDminuspiFidAcc;
        if (IsInFiducialAcceptance(b4Prong.fPt, b4Prong.fYBs))
            b4Prong.fSelBit |= kBstoDsminuspiFidAcc;
        if (IsInFiducialAcceptance(b4Prong.fPt, b4Prong.fYLb))
            b4Prong.fSelBit |= kLbtoLcpluspiFidAcc;
    }

    fBeauty4Prong.push_back(b4Prong);
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::FillGenerated(AliAODMCParticle *part, int origin, int decay, bool dauInAcc)
{
    GenHadron genCharm;
    genCharm.fPt = part->Pt();
    genCharm.fY = part->Y();
    genCharm.fGenLabel = part->GetLabel();
    genCharm.fCandType = 0;
    if (origin == 4)
        genCharm.fCandType |= kPrompt;
    else if (origin == 5)
        genCharm.fCandType |= kFeedDown;
    if (dauInAcc)
        genCharm.fCandType |= kHasDauInAcc;
    if (IsInFiducialAcceptance(genCharm.fPt, genCharm.fY))
        genCharm.fCandType |= kIsInFidAcc;
    genCharm.fDecay = decay;

    fGenHadron.push_back(genCharm);
}

//________________________________________________________________________
bool AliAnalysisTaskSECharmTriggerStudy::RecalcOwnPrimaryVertex(AliAODRecoDecayHF *cand)
{
    AliAODVertex *recvtx = cand->RemoveDaughtersFromPrimaryVtx(fAOD);
    if (!recvtx)
    {
        AliDebug(2, "Removal of daughter tracks failed");
        return false;
    }

    cand->SetOwnPrimaryVtx(recvtx);
    cand->RecalculateImpPars(recvtx, fAOD);
    delete recvtx;

    return true;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::CleanOwnPrimaryVertex(AliAODRecoDecayHF *cand, AliAODVertex *origvtx)
{
    cand->UnsetOwnPrimaryVtx();
    if (origvtx)
    {
        cand->SetOwnPrimaryVtx(origvtx);
        delete origvtx;
        origvtx = nullptr;
    }
    cand->RecalculateImpPars(cand->GetPrimaryVtx(), fAOD);
}

//________________________________________________________________________
bool AliAnalysisTaskSECharmTriggerStudy::AreDauInAcc(int nProng, int *labDau)
{
    for (int iProng = 0; iProng < nProng; iProng++)
    {
        AliAODMCParticle *mcPartDaughter = dynamic_cast<AliAODMCParticle *>(fMCArray->At(labDau[iProng]));
        if (!mcPartDaughter)
        {
            return false;
        }
        double eta = mcPartDaughter->Eta();
        double pt = mcPartDaughter->Pt();
        if (TMath::Abs(eta) > 0.9 || pt < 0.1)
        {
            return false;
        }
    }
    return true;
}

//________________________________________________________________________
bool AliAnalysisTaskSECharmTriggerStudy::IsInFiducialAcceptance(double pt, double y)
{
    if (pt > 5.)
    {
        if (TMath::Abs(y) > 0.8)
            return false;
    }
    else
    {
        double maxFiducialY = -0.2 / 15 * pt * pt + 1.9 / 15 * pt + 0.5;
        double minFiducialY = 0.2 / 15 * pt * pt - 1.9 / 15 * pt - 0.5;
        if (y < minFiducialY || y > maxFiducialY)
            return false;
    }

    return true;
}

//________________________________________________________________
AliAODVertex *AliAnalysisTaskSECharmTriggerStudy::ReconstructDisplVertex(const AliVVertex *primary, TObjArray *tracks, double bField, double dispersion)
{
    //
    // Helper function to recalculate a vertex.
    //

    AliESDVertex *vertexESD = nullptr;
    AliAODVertex *vertexAOD = nullptr;

    AliVertexerTracks vertexer;
    vertexer.SetFieldkG(bField);

    vertexer.SetVtxStart((AliESDVertex *)primary); //primary vertex
    vertexESD = dynamic_cast<AliESDVertex *>(vertexer.VertexForSelectedESDTracks(tracks));
    if (!vertexESD)
        return vertexAOD;

    if (vertexESD->GetNContributors() != tracks->GetEntriesFast())
    {
        delete vertexESD;
        vertexESD = nullptr;
        return vertexAOD;
    }

    // convert to AliAODVertex
    double pos[3], cov[6], chi2perNDF;
    for (int a = 0; a < 3; a++)
        pos[a] = 0.;
    for (int b = 0; b < 6; b++)
        cov[b] = 0.;
    chi2perNDF = 0;

    vertexESD->GetXYZ(pos);       // position
    vertexESD->GetCovMatrix(cov); //covariance matrix

    double vertRadius2 = pos[0] * pos[0] + pos[1] * pos[1];
    if (vertRadius2 > 8.)
    { //(2.82)^2 radius beam pipe
        delete vertexESD;
        vertexESD = nullptr;
        return vertexAOD;
    }

    chi2perNDF = vertexESD->GetChi2toNDF();
    dispersion = vertexESD->GetDispersion();
    delete vertexESD;
    vertexESD = nullptr;
    int nprongs = tracks->GetEntriesFast();
    vertexAOD = new AliAODVertex(pos, cov, chi2perNDF, nullptr, -1, AliAODVertex::kUndef, nprongs);

    return vertexAOD;
}
