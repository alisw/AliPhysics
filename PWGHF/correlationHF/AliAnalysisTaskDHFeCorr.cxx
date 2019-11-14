#include "AliAnalysisTaskDHFeCorr.h"

//ROOT Includes
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TMath.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TTree.h"

//AliRoot and AliPhysics includes
#include "AliVEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

//Electron
#include "AliKFParticle.h"

//Dmeson
#include "AliRDHFCuts.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliVertexingHFUtils.h"

//PID 
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVEventHandler.h"

//multiplicity
#include "AliMultSelection.h"

//STL includes
#include <iostream>
#include <exception>
//#include <algorithm>

ClassImp(AliAnalysisTaskDHFeCorr)

AliAnalysisTaskDHFeCorr::AliAnalysisTaskDHFeCorr() : AliAnalysisTaskSE() {}

AliAnalysisTaskDHFeCorr::AliAnalysisTaskDHFeCorr(const char *name) : AliAnalysisTaskSE(name) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    DefineOutput(3, TList::Class());
    DefineOutput(4, TTree::Class());
    DefineOutput(5, TTree::Class());
    DefineOutput(6, TTree::Class());
}

//_____________________________________________________________________________
void AliAnalysisTaskDHFeCorr::UserCreateOutputObjects() {
    fDMesonRequirements.fDMesonCuts->GetPidHF()->SetPidResponse(fInputHandler->GetPIDResponse());

    fOptEvent.SetOwner(kTRUE);
    fOptElectron.SetOwner(kTRUE);
    fOptDMeson.SetOwner(kTRUE);

    fEventTree = std::unique_ptr<TTree>(new TTree("event", "event"));
    fElectronTree = std::unique_ptr<TTree>(new TTree("electron", "electron"));
    fDmesonTree = std::unique_ptr<TTree>(new TTree("dmeson", "dmeson"));

    AddElectronVariables(fElectronTree);
    AddDMesonVariables(fDmesonTree, fDmesonSpecies);
    AddEventVariables(fEventTree);

    if (fIsMC) {
        AddDMesonMCVariables(fDmesonTree);
        AddElectronMCVariables(fElectronTree);
    }

    fEventQABeforeCuts = CreateQAEvents("before", fOptEvent);
    fEventQAAfterCuts = CreateQAEvents("after", fOptEvent);

    fElectronQABeforeCuts = CreateQAElectrons(fElectronOptConfig, "electron", "before", fOptElectron);
    fElectronQAAfterTrackCuts = CreateQAElectrons(fElectronOptConfig, "electron", "after_track_cuts", fOptElectron);
    fElectronQAAfterCuts = CreateQAElectrons(fElectronOptConfig, "electron", "after", fOptElectron);

    fDMesonQABeforeCuts = CreateQADMeson(fDMesonOptConfig, "dmeson", "before", fOptDMeson);
    fDMesonQAAfterCuts = CreateQADMeson(fDMesonOptConfig, "dmeson", "after", fOptDMeson);

    PostOutput();
}

void AliAnalysisTaskDHFeCorr::AddEventVariables(std::unique_ptr<TTree> &tree) {
    tree->Branch("RunNumber", &fRunNumber);
    tree->Branch("EventNumber", &fEventNumber);
    tree->Branch("VtxZ", &fVtxZ);
    tree->Branch("Centrality", &fCentrality);
}


void AliAnalysisTaskDHFeCorr::AddElectronVariables(std::unique_ptr<TTree> &tree) {
    tree->Branch("RunNumber", &fElectron.fRunNumber);
    tree->Branch("EventNumber", &fElectron.fEventNumber);

    tree->Branch("Pt", &fElectron.fPt);
    tree->Branch("Eta", &fElectron.fEta);
    tree->Branch("Phi", &fElectron.fPhi);
    tree->Branch("ID", &fElectron.fID);

    if (!fReducedElectronInfo) {
        tree->Branch("Charge", &fElectron.fCharge);
        tree->Branch("P", &fElectron.fP);
        tree->Branch("NClsTPC", &fElectron.fNClsTPC);
        tree->Branch("NClsTPCDeDx", &fElectron.fNClsTPCDeDx);
        tree->Branch("NITSCls", &fElectron.fNITSCls);
        tree->Branch("ITSHitFirstLayer", &fElectron.fITSHitFirstLayer);
        tree->Branch("ITSHitSecondLayer", &fElectron.fITSHitSecondLayer);
        tree->Branch("DCAxy", &fElectron.fDCAxy);
        tree->Branch("DCAz", &fElectron.fDCAz);
        tree->Branch("TPCNSigma", &fElectron.fTPCNSigma);
        tree->Branch("TOFNSigma", &fElectron.fTOFNSigma);
        tree->Branch("InvMassPartnersULS", &fElectron.fInvMassPartnersULS);
        tree->Branch("InvMassPartnersLS", &fElectron.fInvMassPartnersLS);
    }
}

void AliAnalysisTaskDHFeCorr::AddDMesonVariables(std::unique_ptr<TTree> &tree,
                                                 AliAnalysisTaskDHFeCorr::DMeson_t meson_species) {
    //Event information
    tree->Branch("RunNumber", &fDmeson.fRunNumber);
    tree->Branch("EventNumber", &fDmeson.fEventNumber);
    tree->Branch("ID", &fDmeson.fID);
    tree->Branch("IsParticleCandidate", &fDmeson.fIsParticleCandidate);

    //Basic information
    tree->Branch("Pt", &fDmeson.fPt);
    tree->Branch("Eta", &fDmeson.fEta);
    tree->Branch("Phi", &fDmeson.fPhi);
    tree->Branch("Y", &fDmeson.fY);
    tree->Branch("InvMass", &fDmeson.fInvMass);
    tree->Branch("ReducedChi2", &fDmeson.fReducedChi2);

    //Topological information
    tree->Branch("DecayLength", &fDmeson.fDecayLength);
    tree->Branch("DecayLengthXY", &fDmeson.fDecayLengthXY);

    tree->Branch("NormDecayLength", &fDmeson.fNormDecayLength);
    tree->Branch("NormDecayLengthXY", &fDmeson.fNormDecayLengthXY);

    tree->Branch("CosP", &fDmeson.fCosP);
    tree->Branch("CosPXY", &fDmeson.fCosPXY);

    tree->Branch("ImpParXY", &fDmeson.fImpParXY);
    tree->Branch("DCA", &fDmeson.fDCA);

    tree->Branch("Normd0MeasMinusExp", &fDmeson.fNormd0MeasMinusExp);
    tree->Branch("PtDaughter", &fDmeson.fPtDaughters);
    tree->Branch("D0Daughter", &fDmeson.fD0Daughters);

    tree->Branch("IDDaughters", &fDmeson.fIDDaughters);

    //PID
    tree->Branch("SelectionStatusDefaultPID", &fDmeson.fSelectionStatusDefaultPID);

    tree->Branch("NSigmaTPCDaughters0", &(fDmeson.fNSigmaTPCDaughters[0]));
    tree->Branch("NSigmaTPCDaughters1", &(fDmeson.fNSigmaTPCDaughters[1]));

    tree->Branch("NSigmaTOFDaughters0", &(fDmeson.fNSigmaTOFDaughters[0]));
    tree->Branch("NSigmaTOFDaughters1", &(fDmeson.fNSigmaTOFDaughters[1]));

    //Meson depended variables
    switch (meson_species) {
        case AliAnalysisTaskDHFeCorr::kD0: {
            tree->Branch("CosTs", &fDmeson.fCosTs);
        }
            break;
        case AliAnalysisTaskDHFeCorr::kDplus: {
            tree->Branch("SigmaVertex", &fDmeson.fSigmaVertex);
            tree->Branch("NSigmaTPCDaughters2", &(fDmeson.fNSigmaTPCDaughters[2]));
            tree->Branch("NSigmaTOFDaughters2", &(fDmeson.fNSigmaTOFDaughters[2]));
        }
            break;
        case AliAnalysisTaskDHFeCorr::kDstar: {
            tree->Branch("AngleD0dkpPisoft", &fDmeson.fAngleD0dkpPisoft);
            tree->Branch("CosTs", &fDmeson.fCosTs);
        }
            break;
    }
}

void AliAnalysisTaskDHFeCorr::AddElectronMCVariables(std::unique_ptr<TTree> &tree) {
    tree->Branch("PtMC", &fElectron.fPtMC);
    tree->Branch("PhiMC", &fElectron.fPhiMC);
    tree->Branch("EtaMC", &fElectron.fEtaMC);
    tree->Branch("Origin", &fElectron.fOrigin);
    tree->Branch("FirstMotherPDG", &fElectron.fFirstMotherPDG);
    tree->Branch("FirstMotherPt", &fElectron.fFirstMotherPt);
    tree->Branch("SecondMotherPDG", &fElectron.fSecondMotherPDG);
    tree->Branch("SecondMotherPt", &fElectron.fSecondMotherPt);
}

void AliAnalysisTaskDHFeCorr::AddDMesonMCVariables(std::unique_ptr<TTree> &tree) {
    tree->Branch("PtMC", &fDmeson.fPtMC);
    tree->Branch("IsD", &fDmeson.fIsD);
    tree->Branch("IsParticle", &fDmeson.fIsParticle);
    tree->Branch("IsPrompt", &fDmeson.fIsPrompt);
}

std::vector<AliDHFeCorr::AliElectron> AliAnalysisTaskDHFeCorr::ElectronAnalysis() {
    auto aod_event = dynamic_cast<AliAODEvent *>(InputEvent());

    std::vector<AliDHFeCorr::AliElectron> tracks;
    tracks.reserve(static_cast<unsigned long>(aod_event->GetNumberOfTracks()));

    for (int i(0); i < aod_event->GetNumberOfTracks(); i++) {
        auto particle = AliDHFeCorr::AliElectron();
        particle.fTrack = dynamic_cast<AliAODTrack *>(aod_event->GetTrack(i));
        tracks.push_back(particle);
    }

    FillElectronInformation(tracks);

    auto selected_tracks = FilterElectronsTracking(fElectronRequirements, tracks);
    auto selected_electrons = FilterElectronsPID(fElectronRequirements, selected_tracks);

    FillElectronMCInfo(selected_electrons);

    auto selected_partner_tracks = FilterElectronsTracking(fPartnerElectronRequirements, tracks);
    auto partner_electrons = FilterElectronsPID(fPartnerElectronRequirements, selected_partner_tracks);

    if (!fReducedElectronInfo) {
        for (auto &electron: selected_electrons)
            FindNonHFe(electron, partner_electrons);
    }

    //fill the electron (track) QA before applying the cuts
    FillElectronQA(tracks, fElectronQABeforeCuts);
    //fill the electron QA after track selection
    FillElectronQA(selected_tracks, fElectronQAAfterTrackCuts);
    //fill the electron QA after applying the cuts
    FillElectronQA(selected_electrons, fElectronQAAfterCuts);

    return selected_electrons;
}

std::vector<AliDHFeCorr::AliDMeson> AliAnalysisTaskDHFeCorr::DMesonAnalysis() {
    auto aod_event = dynamic_cast<AliAODEvent *>(InputEvent());

    //Move D meson candidates to vectors
    const TClonesArray *dmeson_candidates = dynamic_cast<TClonesArray *>(aod_event->GetList()->FindObject(
            fgkDMesonListName.at(fDmesonSpecies).c_str()));

    auto d_mesons = FillDMesonInfo(dmeson_candidates, fDmesonSpecies, fDMesonRequirements);
    //Select D mesons
    auto selected_d_mesons = FilterDmesons(d_mesons, fDMesonRequirements, aod_event, fgkDMesonPDG.at(fDmesonSpecies));

    FillDmesonQA(d_mesons, fDMesonQABeforeCuts, fDmesonSpecies);
    FillDmesonQA(selected_d_mesons, fDMesonQAAfterCuts, fDmesonSpecies);

    FillDmesonMCInfo(selected_d_mesons, fDmesonSpecies);

    return selected_d_mesons;
}

//_____________________________________________________________________________
void AliAnalysisTaskDHFeCorr::UserExec(Option_t *) {
    CheckConfiguration();
    SetRunAndEventNumber();

    //Set Global event variables
    fVtxZ = InputEvent()->GetPrimaryVertex()->GetZ();

    auto multiplicity_selection = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
    if (multiplicity_selection)
        fCentrality = multiplicity_selection->GetMultiplicityPercentile(fEventSelection.fMultiEstimator);

    auto aod_event = dynamic_cast<AliAODEvent *>(InputEvent());
    if (!aod_event)
        return;

    FillEventQA(fEventQABeforeCuts);

    if (!IsSelectedEvent(aod_event, fEventSelection)) {
        PostOutput();
        return;
    }

    auto selected_d_mesons = DMesonAnalysis();
    auto selected_electrons = ElectronAnalysis();

    //In case no D meson AND electron is present, the event is not saved at all
    if ((selected_d_mesons.empty()) && (selected_electrons.empty())) {
        PostOutput();
        return;
    }

    if (!fIsEffMode && !fKeepAllCandidates) {
        //Ignores event in case no electron and D meson is found, depending on the configuration
        if ((selected_d_mesons.empty()) || (selected_electrons.empty())) {
            PostOutput();
            return;
        }
    }

    //Fill the plots in case both D meson and electron is found the event
    FillEventQA(fEventQAAfterCuts);

    if (fIsEffMode && fIsMC) {
        selected_d_mesons = FilterTrueDMesons(selected_d_mesons);
    }

    FillElectronTree(selected_electrons);
    FillDMesonTree(selected_d_mesons);
    fEventTree->Fill();

    //Post information to the output
    PostOutput();

    //increase the number of the event used for the Trees
    fEventNumber++;
}

void AliAnalysisTaskDHFeCorr::CheckConfiguration() const {
    auto event_list = InputEvent()->GetList();
    auto mc_information = dynamic_cast<TClonesArray *>(event_list->FindObject(AliAODMCParticle::StdBranchName()));

    if (!mc_information) {
        if (fIsMC && fIsEffMode) {
            throw std::invalid_argument(
                    "The analysis has no MC information and it is set as MC analysis. Check the data/configuration.");
        }
        return;
    }

    if (fgkDMesonDaughterAliPID.count(fDmesonSpecies) < 1)
        throw std::invalid_argument("The daughter list (AliPID) is not defined. Check fgkDMesonDaughterAliPID.");

    if (fgkDMesonDaughterPDG.count(fDmesonSpecies) < 1)
        throw std::invalid_argument("The daughter list PDG is not defined. Check fgkDMesonDaughterPDG.");

    if (fgkDMesonListName.count(fDmesonSpecies) < 1)
        throw std::invalid_argument("The daughter AOD list name is not defined. Check fgkDMesonListName.");
}

void AliAnalysisTaskDHFeCorr::FillElectronInformation(std::vector<AliDHFeCorr::AliElectron> &electrons) {
    auto aod_event = dynamic_cast<AliAODEvent *>(InputEvent());
    const auto pid_response = fInputHandler->GetPIDResponse();

    for (auto &candidate: electrons) {
        const auto track = candidate.fTrack;
        candidate.fRunNumber = fRunNumber;
        candidate.fEventNumber = fEventNumber;
        candidate.fID = TMath::Abs(track->GetID());

        candidate.fCharge = track->Charge();
        candidate.fPt = track->Pt();
        candidate.fP = track->P();
        candidate.fEta = track->Eta();
        candidate.fPhi = track->Phi();

        candidate.fNClsTPC = track->GetTPCNcls();
        candidate.fNClsTPCDeDx = track->GetTPCsignalN();
        candidate.fNITSCls = track->GetITSNcls();

        candidate.fITSHitFirstLayer = track->HasPointOnITSLayer(0);
        candidate.fITSHitSecondLayer = track->HasPointOnITSLayer(1);

        Double_t d0z0[2] = {-999., -999.};
        Double_t cov[3] = {-999., -999., -999.};
        const AliVVertex *primaryVertex = aod_event->GetPrimaryVertex();
        AliAODTrack copyTrack = AliAODTrack(*track); //Copy track to not alter its parameters

        if (copyTrack.PropagateToDCA(primaryVertex, aod_event->GetMagneticField(), 20., d0z0, cov)) {
            candidate.fDCAxy = d0z0[0];
            candidate.fDCAz = d0z0[1];
        }

        candidate.fTPCNSigma = pid_response->NumberOfSigmasTPC(track, AliPID::kElectron);
        candidate.fTOFNSigma = pid_response->NumberOfSigmasTOF(track, AliPID::kElectron);
    }
}

std::vector<AliDHFeCorr::AliElectron>
AliAnalysisTaskDHFeCorr::FilterElectronsTracking(AliDHFeCorr::AliElectronSelection electronSelection,
                                                 const std::vector<AliDHFeCorr::AliElectron> &electrons) {
    std::vector<AliDHFeCorr::AliElectron> selected_electrons;
    selected_electrons.reserve(electrons.size());

    for (const auto &candidate: electrons) {
        if ((candidate.fPt < electronSelection.fPtMin) || (candidate.fPt > electronSelection.fPtMax))
            continue;

        if ((candidate.fEta < electronSelection.fEtaMin) || (candidate.fEta > electronSelection.fPtMax))
            continue;

        if (!candidate.fTrack->TestFilterBit(electronSelection.fFilterBit))
            continue;

        if (candidate.fNClsTPC < electronSelection.fTPCClsMin)
            continue;

        if (candidate.fNClsTPCDeDx < electronSelection.fTPCClsDeDxMin)
            continue;

        if (candidate.fNITSCls < electronSelection.fITSClsMin)
            continue;

        if (!FulfilPixelSelection(candidate, electronSelection.fITSPixel))
            continue;

        if (candidate.fDCAxy > electronSelection.fDCAxy || candidate.fDCAz > electronSelection.fDCAz)
            continue;

        auto electron = AliDHFeCorr::AliElectron(candidate);
        selected_electrons.push_back(electron);
    }

    selected_electrons.shrink_to_fit();
    return selected_electrons;
}

std::vector<AliDHFeCorr::AliElectron>
AliAnalysisTaskDHFeCorr::FilterElectronsPID(AliDHFeCorr::AliElectronSelection electronSelection,
                                            const std::vector<AliDHFeCorr::AliElectron> &electrons) {

    std::vector<AliDHFeCorr::AliElectron> selected_electrons;
    selected_electrons.reserve(electrons.size());

    for (const auto &candidate: electrons) {
        if ((candidate.fTPCNSigma < electronSelection.fTPCNSigmaMin) ||
            (candidate.fTPCNSigma > electronSelection.fTPCNSigmaMax))
            continue;

        if (electronSelection.fRequireTOF) {
            if ((candidate.fTOFNSigma < electronSelection.fTOFNSigmaMin) ||
                (candidate.fTOFNSigma > electronSelection.fTOFNSigmaMax))
                continue;
        }
        selected_electrons.push_back(AliDHFeCorr::AliElectron(candidate));
    }

    selected_electrons.shrink_to_fit();
    return selected_electrons;
}

void AliAnalysisTaskDHFeCorr::FillElectronMCInfo(std::vector<AliDHFeCorr::AliElectron> &electrons) {
    const auto mc_information = dynamic_cast<TClonesArray *>(InputEvent()->GetList()->FindObject(
            AliAODMCParticle::StdBranchName()));

    if (!mc_information)
        return;

    for (auto &candidate: electrons) {
        auto track = candidate.fTrack;
        auto label = TMath::Abs(track->GetLabel());

        auto mc_part = dynamic_cast<AliAODMCParticle *>(mc_information->At(label));

        if (!mc_part)
            continue;

        candidate.fPtMC = mc_part->Pt();
        candidate.fEtaMC = mc_part->Eta();
        candidate.fPhiMC = mc_part->Phi();
        candidate.fPDG = mc_part->PdgCode();
        candidate.fOrigin = AliVertexingHFUtils::CheckOrigin(mc_information, mc_part, true);

        if (mc_part->GetMother() > 0) {
            const auto mc_mother = dynamic_cast<AliAODMCParticle *>(mc_information->At(mc_part->GetMother()));
            candidate.fFirstMotherPDG = mc_mother->GetPdgCode();
            candidate.fFirstMotherPt = mc_mother->Pt();

            if (mc_mother->GetMother() > 0) {
                const auto mc_grandmother = dynamic_cast<AliAODMCParticle *>(mc_information->At(
                        mc_mother->GetMother()));
                candidate.fSecondMotherPDG = mc_grandmother->GetPdgCode();
                candidate.fSecondMotherPt = mc_grandmother->Pt();
            }
        }
    }
}


void AliAnalysisTaskDHFeCorr::FillDmesonMCInfo(std::vector<AliDHFeCorr::AliDMeson> &d_mesons,
                                               AliAnalysisTaskDHFeCorr::DMeson_t meson_species) const {
    if (!fIsMC) {
        return;
    }

    const auto mc_information = dynamic_cast<TClonesArray *>(InputEvent()->GetList()->FindObject(
            AliAODMCParticle::StdBranchName()));

    const auto pdg_daughters_int = fgkDMesonDaughterPDG.at(meson_species);
    const auto pdg_dmeson = fgkDMesonPDG.at(meson_species);

    for (auto &candidate: d_mesons) {
        const auto reco_cand = candidate.fRecoObj;
        const Int_t label = reco_cand->MatchToMC(pdg_dmeson, mc_information, pdg_daughters_int.size(),
                                                 &pdg_daughters_int[0]);

        if (label < 0)
            continue;

        candidate.fIsD = kTRUE;
        const auto mc_part = dynamic_cast<AliAODMCParticle *>(mc_information->At(label));
        candidate.fPtMC = mc_part->Pt();
        const auto origin = AliVertexingHFUtils::CheckOrigin(mc_information, mc_part, true); //Prompt = 4, FeedDown = 5
        candidate.fIsParticle = kFALSE;
        if (mc_part->PdgCode() > 0) {
            candidate.fIsParticle = kTRUE;
        }
        if (origin != 5) {
            candidate.fIsPrompt = kTRUE;
        }
    }
}

std::vector<AliDHFeCorr::AliDMeson>
AliAnalysisTaskDHFeCorr::FilterTrueDMesons(const std::vector<AliDHFeCorr::AliDMeson> &d_mesons) {
    std::vector<AliDHFeCorr::AliDMeson> true_d;
    true_d.reserve(d_mesons.size());

    for (const auto &candidate: d_mesons) {
        if (candidate.fIsD) {
            true_d.push_back(candidate);
        }
    }
    true_d.shrink_to_fit();
    return true_d;
}

void AliAnalysisTaskDHFeCorr::SetRunAndEventNumber() {
    //Based on AliAnalysisTaskSEHFTreeCreator::GetEvID

    fRunNumber = static_cast<UInt_t >(InputEvent()->GetRunNumber());

    std::string current_file_name = ((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetCurrentFile()))->GetName();

    if (fCurrentFile != current_file_name) {
        fEventNumber = 0;
        TObjArray *path = TString(current_file_name).Tokenize("/");
        TString s = (dynamic_cast<TObjString *>(path->At(((path->GetLast()) - 1))))->GetString();
        fDirNum = (unsigned int) s.Atoi();
        delete path;
    }

    Long64_t ev_number = Entry();

    fEventNumber = (unsigned int) ev_number + (unsigned int) (fDirNum << 17);
}


void AliAnalysisTaskDHFeCorr::PostOutput() {
    PostData(1, &fOptEvent);
    PostData(2, &fOptElectron);
    PostData(3, &fOptDMeson);
    PostData(4, fElectronTree.get());
    PostData(5, fDmesonTree.get());
    PostData(6, fEventTree.get());
}


void AliAnalysisTaskDHFeCorr::FillElectronQA(const std::vector<AliDHFeCorr::AliElectron> &tracks,
                                             AliDHFeCorr::AliElectronQAHistograms &histograms) {
    for (const auto &electron: tracks) {
        histograms.fPtEtaPhi->Fill(electron.fPt, electron.fEta, electron.fPhi);
        histograms.fTPCNCls->Fill(electron.fNClsTPC);
        histograms.fTPCNClsDeDx->Fill(electron.fNClsTPCDeDx);
        histograms.fITSCls->Fill(electron.fNITSCls);
        histograms.fDCAz->Fill(electron.fDCAz);
        histograms.fDCAxy->Fill(electron.fDCAxy);
        //PID
        histograms.fTPCNsigmaPt->Fill(electron.fPt, electron.fTPCNSigma);
        histograms.fTPCNsigmaP->Fill(electron.fP, electron.fTPCNSigma);
        histograms.fTOFNsigmaPt->Fill(electron.fPt, electron.fTOFNSigma);
        histograms.fTOFNsigmaP->Fill(electron.fP, electron.fTOFNSigma);
        histograms.fTPCTOFNSigma->Fill(electron.fTPCNSigma, electron.fTOFNSigma);
    }
}

void AliAnalysisTaskDHFeCorr::FillDmesonQA(const std::vector<AliDHFeCorr::AliDMeson> &d_mesons,
                                           AliDHFeCorr::AliDMesonQAHistos &histograms,
                                           AliAnalysisTaskDHFeCorr::DMeson_t meson_species) {
    for (const auto &candidate: d_mesons) {
        const auto reco_cand = candidate.fRecoObj;
        histograms.fPtEtaPhi->Fill(reco_cand->Pt(), reco_cand->Eta(), reco_cand->Phi());

        switch (meson_species) {
            case AliAnalysisTaskDHFeCorr::kD0: {
                const auto d0_cand = dynamic_cast<AliAODRecoDecayHF2Prong *>(reco_cand);
                if (candidate.fIsParticleCandidate)
                    histograms.fPtInvMass->Fill(reco_cand->Pt(), d0_cand->InvMassD0());
                else
                    histograms.fPtInvMass->Fill(reco_cand->Pt(), d0_cand->InvMassD0bar());
                break;
            }
            case AliAnalysisTaskDHFeCorr::kDplus: {
                const auto dplus_cand = dynamic_cast<AliAODRecoDecayHF3Prong *>(reco_cand);
                histograms.fPtInvMass->Fill(reco_cand->Pt(), dplus_cand->InvMassDplus());
                break;
            }
            case AliAnalysisTaskDHFeCorr::kDstar: {
                const auto dstar_cand = dynamic_cast<AliAODRecoCascadeHF *>(reco_cand);
                histograms.fPtInvMass->Fill(reco_cand->Pt(), dstar_cand->InvMassDstarKpipi());
                break;
            }
        }
    }
}

void AliAnalysisTaskDHFeCorr::FillEventQA(AliDHFeCorr::AliEventQAHistograms &histograms) {
    histograms.fNEvents->Fill(0);
    histograms.fVertexZ->Fill(fVtxZ);
    histograms.fCentrality->Fill(fCentrality);
}

void AliAnalysisTaskDHFeCorr::FindNonHFe(AliDHFeCorr::AliElectron &main_electron,
                                         const std::vector<AliDHFeCorr::AliElectron> &partners) const {

    //Reset all the vectors and reserve 10 (it should not have too many partners)
    std::vector<Float_t> inv_mass_uls, inv_mass_ls;
    std::vector<UInt_t> id_uls, id_ls;

    inv_mass_uls.reserve(10);
    inv_mass_ls.reserve(10);
    id_uls.reserve(10);
    id_ls.reserve(10);

    //Set the magnetic field
    AliKFParticle::SetField(InputEvent()->GetMagneticField());

    const auto track_main = main_electron.fTrack;
    const auto charge_main = track_main->Charge();

    int pdg_main = 11; //electron pdg is always 11
    if (charge_main > 0.) { //change sign in case it is a positron
        pdg_main = -1 * pdg_main;
    }

    const auto kfparticle_main = AliKFParticle(*static_cast<AliVTrack *>(track_main), pdg_main);

    for (const auto &partner: partners) {
        const auto track_partner = partner.fTrack;

        if (TMath::Abs(track_partner->GetID()) == TMath::Abs(track_main->GetID()))
            continue;

        const auto charge_partner = track_partner->Charge();

        int pdg_partner = 11; //electron pdg is always 11

        if (charge_partner > 0.) { //change sign in case it is a positron
            pdg_partner = -1 * pdg_partner;
        }

        const auto kfparticle_partner = AliKFParticle(*static_cast<AliVTrack *>(track_partner), pdg_partner);

        AliKFParticle photon(kfparticle_main, kfparticle_partner);

        //Reject in case no degree of freedom is smaller than 1
        if (photon.GetNDF() < 1)
            continue;

        const auto reduced_chi2 = photon.GetChi2() / photon.GetNDF();

        if (TMath::Abs(reduced_chi2) > fPhotonSelection.fReducedChi2Max)
            continue;

        const auto pair_mass = photon.GetMass();

        if (pair_mass > fPhotonSelection.fInvMassMax)
            continue;

        //If the pair arrived here, it fulfils all the requirements, so it will be saved
        if ((charge_main * charge_partner) > 0) { //++ or --, Like-sign pair
            inv_mass_ls.push_back(pair_mass);
            id_ls.push_back(static_cast<unsigned int &&>(TMath::Abs(track_partner->GetID())));
        } else {
            inv_mass_uls.push_back(pair_mass);
            id_uls.push_back(static_cast<unsigned int &&>(TMath::Abs(track_partner->GetID())));
        }
    }

    inv_mass_uls.shrink_to_fit();
    inv_mass_ls.shrink_to_fit();
    id_uls.shrink_to_fit();
    id_ls.shrink_to_fit();

    main_electron.fInvMassPartnersULS = inv_mass_uls;
    main_electron.fInvMassPartnersLS = inv_mass_ls;
    main_electron.fPartnersULSID = id_uls;
    main_electron.fPartnersLSID = id_ls;

}

void AliAnalysisTaskDHFeCorr::FillElectronTree(const std::vector<AliDHFeCorr::AliElectron> &tracks) {
    for (const auto &electron: tracks) {
        fElectron = AliDHFeCorr::AliElectron(electron);
        fElectronTree->Fill();
    }
}


void AliAnalysisTaskDHFeCorr::FillDMesonTree(const std::vector<AliDHFeCorr::AliDMeson> &d_mesons) {
    for (const auto &dmeson: d_mesons) {
        fDmeson = AliDHFeCorr::AliDMeson(dmeson);
        fDmesonTree->Fill();
    }
}

std::vector<Float_t> AliAnalysisTaskDHFeCorr::MakeBins(Float_t start, Float_t end, Int_t n_bins) {
    const auto step = (end - start) / n_bins;
    std::vector<Float_t> bins;
    bins.reserve(n_bins + 1);

    for (Int_t i(0); i < (n_bins + 1); i++)
        bins.push_back(start + i * step);

    return bins;
}


std::vector<AliDHFeCorr::AliDMeson> AliAnalysisTaskDHFeCorr::FillDMesonInfo(const TClonesArray *dmeson_candidates,
                                                                            AliAnalysisTaskDHFeCorr::DMeson_t meson_species,
                                                                            AliDHFeCorr::AliDMesonSelection &dmeson_selection) const {
    AliAnalysisVertexingHF vertexing_HF;
    std::vector<AliDHFeCorr::AliDMeson> d_mesons;
    d_mesons.reserve(2 * dmeson_candidates->GetEntriesFast());

    const float b_field = InputEvent()->GetMagneticField();

    const auto pdg_dmeson = fgkDMesonPDG.at(meson_species);
    const auto pdg_daughters = fgkDMesonDaughterPDG.at(meson_species);
    const auto id_daughter_alipid = fgkDMesonDaughterAliPID.at(fDmesonSpecies);

    const auto aod_event = dynamic_cast<AliAODEvent *>(InputEvent());
    const auto pid_response = fInputHandler->GetPIDResponse();


    for (int i(0); i < dmeson_candidates->GetEntriesFast(); i++) {
        auto cand = AliDHFeCorr::AliDMeson();
        cand.fRecoObj = dynamic_cast<AliAODRecoDecayHF *>(dmeson_candidates->At(i));
        cand.fID = static_cast<UInt_t>(i); // ID of the D meson in the event

        const auto reco_candidate = cand.fRecoObj;
        if (!reco_candidate)
            continue;
        //Fill missing info not saved in the AOD
        switch (meson_species) {
            case AliAnalysisTaskDHFeCorr::kD0: {
                if (!vertexing_HF.FillRecoCand(aod_event, dynamic_cast<AliAODRecoDecayHF2Prong *>(reco_candidate)))
                    continue;
                if (!reco_candidate->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts))
                    continue;
                break;
            }
            case AliAnalysisTaskDHFeCorr::kDplus: {
                if (!vertexing_HF.FillRecoCand(aod_event, dynamic_cast<AliAODRecoDecayHF3Prong *>(reco_candidate)))
                    continue;
                if (!reco_candidate->HasSelectionBit(AliRDHFCuts::kDplusCuts))
                    continue;
                break;
            }
            case AliAnalysisTaskDHFeCorr::kDstar: {
                if (!vertexing_HF.FillRecoCand(aod_event, dynamic_cast<AliAODRecoCascadeHF *>(reco_candidate)))
                    continue;
                if (!reco_candidate->HasSelectionBit(AliRDHFCuts::kDstarCuts))
                    continue;
                break;
            }
            default: {
                throw std::invalid_argument("The D meson requested is not defined in the task");
                break;
            }
        }
        //avoid a lot of computational time by rejecting minimum Pt here
        //Reject in case the pt is smaller than the minimum fPtMin
        if (reco_candidate->Pt() < dmeson_selection.fPtMin)
            continue;

        //Is it particle or antiparticle? Is positive -> True, Should work for D+
        //For D0, the Charge() should return 0, but the candidates will be saved twice (to handle reflections)
        cand.fIsParticleCandidate = reco_candidate->Charge() >= 0;
        cand.fRunNumber = fRunNumber;
        cand.fEventNumber = fEventNumber;

        const auto reco_cand = cand.fRecoObj;
        //Recalculate vertex w/o daughters
        AliAODVertex *original_own_vertex(nullptr);

        //GetIsPrimaryWithoutDaughters returns true if the user requested to recalculate the primary vertex
        //In Pb-Pb this option can be turned off
        if (dmeson_selection.fDMesonCuts->GetIsPrimaryWithoutDaughters()) {
            if (reco_cand->GetOwnPrimaryVtx())
                original_own_vertex = new AliAODVertex(*reco_cand->GetOwnPrimaryVtx());

            if (!dmeson_selection.fDMesonCuts->RecalcOwnPrimaryVtx(reco_cand, aod_event)) {
                dmeson_selection.fDMesonCuts->CleanOwnPrimaryVtx(reco_cand, aod_event, original_own_vertex);
                continue;
            }
        }

        cand.fPt = reco_cand->Pt();
        cand.fEta = reco_cand->Eta();
        cand.fPhi = reco_cand->Phi();
        cand.fY = reco_cand->Y(pdg_dmeson);
        cand.fReducedChi2 = reco_cand->GetReducedChi2();
        cand.fDecayLength = reco_cand->DecayLength();
        cand.fDecayLengthXY = reco_cand->DecayLengthXY();
        cand.fNormDecayLength = reco_cand->NormalizedDecayLength();
        cand.fNormDecayLengthXY = reco_cand->NormalizedDecayLengthXY();
        cand.fCosP = reco_cand->CosPointingAngle();
        cand.fCosPXY = reco_cand->CosPointingAngleXY();
        cand.fImpParXY = reco_cand->ImpParXY();
        cand.fDCA = reco_cand->GetDCA();
        cand.fNormd0MeasMinusExp = GetMaxd0MeasMinusExp(reco_cand, b_field);

        //Single particle information
        const auto n_prongs = reco_cand->GetNProngs();

        std::vector<Float_t> pt_daughters;
        pt_daughters.reserve(n_prongs);

        std::vector<Float_t> d0_daughters;
        d0_daughters.reserve(n_prongs);

        for (int j = 0; j < n_prongs; j++) {
            pt_daughters.push_back(reco_cand->PtProng(j));
            d0_daughters.push_back(reco_cand->Getd0Prong(j));
        }
        cand.fPtDaughters = pt_daughters;
        cand.fD0Daughters = d0_daughters;

        std::array<std::vector<Float_t>, 3> pid_info_tpc;
        std::array<std::vector<Float_t>, 3> pid_info_tof;
        std::vector<UInt_t> id_daughters;
        id_daughters.reserve(n_prongs);

        for (int j = 0; j < n_prongs; j++) {

            //TO DO: UNDERSTAND IT IS NOT SAVING THE pid_info_tpc_prong WITH 2X THE NUMBER OF id_daughter_alipid
            const auto track = dynamic_cast<AliAODTrack *>(reco_cand->GetDaughter(j));
            id_daughters.push_back(TMath::Abs(track->GetID()));

            std::vector<Float_t> pid_info_tpc_prong;
            std::vector<Float_t> pid_info_tof_prong;

            for (auto hypothesis: id_daughter_alipid) {
                pid_info_tpc_prong.push_back(pid_response->NumberOfSigmasTPC(track, hypothesis));
                pid_info_tof_prong.push_back(pid_response->NumberOfSigmasTOF(track, hypothesis));
            }

            pid_info_tpc[j] = pid_info_tpc_prong;
            pid_info_tof[j] = pid_info_tof_prong;
        }

        cand.fNSigmaTPCDaughters = pid_info_tpc;
        cand.fNSigmaTOFDaughters = pid_info_tof;
        cand.fIDDaughters = id_daughters;

        if (meson_species == kD0) { //The reflection will be saved later
            const auto d0_candidate = dynamic_cast<AliAODRecoDecayHF2Prong *>(reco_cand);
            cand.fInvMass = d0_candidate->InvMassD0();
            cand.fCosTs = d0_candidate->CosThetaStarD0();
        } else if (meson_species == kDplus) {
            const auto dplus_candidate = dynamic_cast<AliAODRecoDecayHF3Prong *>(reco_cand);
            cand.fSigmaVertex = dplus_candidate->GetSigmaVert();
            cand.fInvMass = dplus_candidate->InvMassDplus();
        } else if (meson_species == kDstar) {
            //const auto dstar_candidate = dynamic_cast<AliAODRecoCascadeHF *>(reco_cand);
        }

        d_mesons.push_back(cand);

        //Duplicate D0 candidates
        if (meson_species == kD0) {
            AliDHFeCorr::AliDMeson reflection(cand);
            reflection.fIsParticleCandidate = kFALSE; //should always be kFALSE, since the first one was kTRUE

            auto d0bar_candidate = dynamic_cast<AliAODRecoDecayHF2Prong *>(reflection.fRecoObj);
            reflection.fInvMass = d0bar_candidate->InvMassD0bar();
            reflection.fCosTs = d0bar_candidate->CosThetaStarD0bar();

            //Invert the vectors related to the reflections
            std::reverse(reflection.fPtDaughters.begin(), reflection.fPtDaughters.end());
            std::reverse(reflection.fD0Daughters.begin(), reflection.fD0Daughters.end());
            std::reverse(reflection.fIDDaughters.begin(), reflection.fIDDaughters.end());

            //Invert manually the NSigma responses
            std::swap(reflection.fNSigmaTPCDaughters[0], reflection.fNSigmaTPCDaughters[1]);
            std::swap(reflection.fNSigmaTOFDaughters[0], reflection.fNSigmaTOFDaughters[1]);

            d_mesons.push_back(reflection);
        }

        dmeson_selection.fDMesonCuts->CleanOwnPrimaryVtx(reco_cand, aod_event, original_own_vertex);
    }

    d_mesons.shrink_to_fit();
    return d_mesons;
}


std::vector<AliDHFeCorr::AliDMeson>
AliAnalysisTaskDHFeCorr::FilterDmesons(const std::vector<AliDHFeCorr::AliDMeson> &dmeson_candidates,
                                       const AliDHFeCorr::AliDMesonSelection &selectionDMeson, AliAODEvent *aod_event,
                                       int pdg_dmeson) {

    std::vector<AliDHFeCorr::AliDMeson> d_mesons;
    d_mesons.reserve(dmeson_candidates.size());

    for (const auto &candidate: dmeson_candidates) {
        const auto d_candidate = candidate.fRecoObj;

        //Check if it is in fiducial acceptance
        if (!selectionDMeson.fDMesonCuts->IsInFiducialAcceptance(d_candidate->Pt(), d_candidate->Y(pdg_dmeson)))
            continue;

        //Use PID only for values smaller than fPtMaxPID
        Bool_t pid_configuration = selectionDMeson.fDMesonCuts->GetIsUsePID();

        if (candidate.fPt > selectionDMeson.fPtMaxPID)
            selectionDMeson.fDMesonCuts->SetUsePID(kFALSE);

        auto selection_status = selectionDMeson.fDMesonCuts->IsSelected(d_candidate, AliRDHFCuts::kAll, aod_event);

        if (selection_status == AliAnalysisTaskDHFeCorr::kNotSelected)
            continue;

        if (candidate.fIsParticleCandidate) {
            if (selection_status == AliAnalysisTaskDHFeCorr::kAntiParticle)
                continue;
        } else {
            if (selection_status == AliAnalysisTaskDHFeCorr::kParticle)
                continue;
        }

        //Set PID status to True to obtain the default PID value
        selectionDMeson.fDMesonCuts->SetUsePID(kTRUE);
        Int_t pid_selection_status = selectionDMeson.fDMesonCuts->IsSelected(d_candidate, AliRDHFCuts::kPID, aod_event);
        //return the PID selection to the original value
        selectionDMeson.fDMesonCuts->SetUsePID(pid_configuration);

        AliDHFeCorr::AliDMeson cand(candidate);
        cand.fSelectionStatusDefaultPID = pid_selection_status;
        d_mesons.push_back(cand);
    }
    d_mesons.shrink_to_fit();

    return d_mesons;
}

bool AliAnalysisTaskDHFeCorr::IsSelectedEvent(AliAODEvent *event, const AliDHFeCorr::AliEventSelection &selection) {
    if (!event)
        return false;
    if ((fVtxZ < selection.fVertexZMin) || (fVtxZ > selection.fVertexZMax))
        return false;

    if (fCentrality < selection.fMultMin || fCentrality > selection.fMultMax)
        return false;

    return true;
}

bool AliAnalysisTaskDHFeCorr::FulfilPixelSelection(const AliDHFeCorr::AliElectron &track, Int_t requirement) {
    const bool hasHitFirstLayer = track.fITSHitFirstLayer;
    const bool hasHitSecondLayer = track.fITSHitSecondLayer;

    if ((hasHitFirstLayer || hasHitSecondLayer) && (requirement == kAny))
        return kTRUE;
    if (hasHitFirstLayer && hasHitSecondLayer && (requirement == kBoth))
        return kTRUE;
    if (hasHitFirstLayer && (requirement == kFirst))
        return kTRUE;
    if (hasHitSecondLayer && (requirement == kSecond))
        return kTRUE;
    if (!hasHitFirstLayer && !hasHitSecondLayer && (requirement == kNone))
        return kTRUE;
    if (!hasHitFirstLayer && hasHitSecondLayer && (requirement == kExclusiveSecond))
        return kTRUE;
    if (hasHitFirstLayer && !hasHitSecondLayer && (requirement == kExclusiveFirst))
        return kTRUE;

    return kFALSE;
}

AliAnalysisTaskDHFeCorr *AliAnalysisTaskDHFeCorr::AddTask(const std::string &name, const std::string &config_file) {
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        std::cout << "It was not possible to connect to the AnalysisManager. Returning null pointer." << std::endl;
        return nullptr;
    }

    if (!mgr->GetInputEventHandler()) {
        std::cout << "It was not possible to connect to the InputEventHandler. Returning null pointer." << std::endl;
        return nullptr;
    }

    auto *task = new AliAnalysisTaskDHFeCorr(name.c_str());
    if (!task) {
        std::cout << "Failed to create the task. Please check the results." << std::endl;
        return nullptr;
    }
    //Configure the task
    if (!task->Configure(config_file, name)) {
        std::cout << "Failed to configure the task. Please check the config file and error above." << std::endl;
        return nullptr;
    }

    task->SelectCollisionCandidates(AliVEvent::kINT7);

    mgr->AddTask(task);

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

    std::string fileName = static_cast<std::string>(AliAnalysisManager::GetCommonFileName());
    fileName += ":DHFeCorrelation_" + name;

    mgr->ConnectOutput(task, 1, mgr->CreateContainer("EventsQA", TList::Class(), AliAnalysisManager::kOutputContainer,
                                                     fileName.c_str()));
    mgr->ConnectOutput(task, 2, mgr->CreateContainer("ElectronQA", TList::Class(), AliAnalysisManager::kOutputContainer,
                                                     fileName.c_str()));
    mgr->ConnectOutput(task, 3, mgr->CreateContainer("DMesonQA", TList::Class(), AliAnalysisManager::kOutputContainer,
                                                     fileName.c_str()));
    mgr->ConnectOutput(task, 4,
                       mgr->CreateContainer("ElectronTree", TTree::Class(), AliAnalysisManager::kOutputContainer,
                                            fileName.c_str()));
    mgr->ConnectOutput(task, 5, mgr->CreateContainer("DMesonTree", TTree::Class(), AliAnalysisManager::kOutputContainer,
                                                     fileName.c_str()));

    mgr->ConnectOutput(task, 6, mgr->CreateContainer("EventTree", TTree::Class(), AliAnalysisManager::kOutputContainer,
                                                     fileName.c_str()));

    return task;
}

bool AliAnalysisTaskDHFeCorr::Configure(std::string config_file, std::string config_name, std::string default_file) {

    //restart the configuration every time the function configure is called (for wagons on LEGO-Trains)
    fYAMLConfig = PWG::Tools::AliYAMLConfiguration();

    fYAMLConfig.AddConfiguration(default_file, "default_configuration");
    // Will be checked first. This is so it can override values in the first added configuration.
    fYAMLConfig.AddConfiguration(config_file, config_name);

    //Read the global task settings
    fYAMLConfig.GetProperty("mc_mode", fIsMC, true);
    fYAMLConfig.GetProperty("calculate_only_efficiency", fIsEffMode, true);
    fYAMLConfig.GetProperty("keep_all_candidates", fKeepAllCandidates, true);
    fYAMLConfig.GetProperty("reduced_electron_info", fReducedElectronInfo, true);

    std::string meson_species;
    fYAMLConfig.GetProperty("d_meson_species", meson_species, true);
    try {
        fDmesonSpecies = fgkDMesonNames.at(meson_species);
    }
    catch (std::exception &exp) {
        std::cout
                << "The D meson chosen is not implemented. The following exception was raised:"
                << exp.what() << std::endl;
        return false;
    }

    if (!ConfigureEventSelection("event", fEventSelection))
        return false;

    if (!ConfigureElectrons("main_electron", fElectronRequirements))
        return false;

    if (!ConfigureElectrons("partner_electron", fPartnerElectronRequirements))
        return false;

    if (!ConfigureElectronOpt("electron_qa", fElectronOptConfig))
        return false;

    if (!ConfigureDMesons("D_meson", fDMesonRequirements))
        return false;

    if (!ConfigureDMesonOpt("d_meson_qa", fDMesonOptConfig))
        return false;

    fYAMLConfig.Initialize();

    return true;
}

bool
AliAnalysisTaskDHFeCorr::ConfigureEventSelection(const std::string &name,
                                                 AliDHFeCorr::AliEventSelection &event_selection) {
    std::cout << "Configuring event selection for: " << name << std::endl;

    AliDHFeCorr::Utils::GetPropertyRange<Float_t>(fYAMLConfig, {name, "properties", "zvxt_range"},
                                                  event_selection.fVertexZMin, event_selection.fVertexZMax, true);

    fYAMLConfig.GetProperty({name, "multiplicity", "estimator"}, event_selection.fMultiEstimator, true);

    AliDHFeCorr::Utils::GetPropertyRange<Float_t>(fYAMLConfig, {name, "multiplicity", "multiplicity_range"},
                                                  event_selection.fMultMin, event_selection.fMultMax, true);
    return true;

}

bool
AliAnalysisTaskDHFeCorr::ConfigureElectrons(const std::string &name,
                                            AliDHFeCorr::AliElectronSelection &electron_selection) {
    std::cout << "Configuring track and PID for: " << name << std::endl;

    AliDHFeCorr::Utils::GetPropertyRange<Float_t>(fYAMLConfig, {name, "track", "pt_range"}, electron_selection.fPtMin,
                                                  electron_selection.fPtMax, true);
    AliDHFeCorr::Utils::GetPropertyRange<Float_t>(fYAMLConfig, {name, "track", "eta_range"}, electron_selection.fEtaMin,
                                                  electron_selection.fEtaMax, true);

    //filter_bit
    std::string filter_bit;
    if (!fYAMLConfig.GetProperty({name, "track", "filter_bit"}, filter_bit, false)) {
        std::cout << "Check Config file. It is not possible to set filter_bit for electrons" << std::endl;
        return false;
    }

    try {
        electron_selection.fFilterBit = fgkAODFilterBitMap.at(filter_bit);
    }
    catch (std::exception &exp) {
        std::cout
                << "Problem to read 'filter_bit'. The following error was raised: " << exp.what() << std::endl;
        return false;
    }

    //min_TPC_cls
    fYAMLConfig.GetProperty({name, "track", "min_TPC_cls"}, electron_selection.fTPCClsMin, true);
    //min_TPC_cls_dedx
    fYAMLConfig.GetProperty({name, "track", "min_TPC_cls_dedx"}, electron_selection.fTPCClsDeDxMin, true);
    //min_ITS_hits
    fYAMLConfig.GetProperty({name, "track", "min_ITS_hits"}, electron_selection.fITSClsMin, true);

    std::string pixel_req;
    if (!fYAMLConfig.GetProperty({name, "track", "pixel_req"}, pixel_req, false)) {
        std::cout << "Check Config file. It is not possible to set pixel_req for electrons" << std::endl;
        return false;
    }

    try {
        electron_selection.fITSPixel = fgkITSPixelMap.at(pixel_req);
    }
    catch (std::exception &exp) {
        std::cout
                << "Problem to read 'pixel_req'. The following error was raised: "
                << exp.what() << std::endl;
        return false;
    }

    fYAMLConfig.GetProperty({name, "track", "dca_z"}, electron_selection.fDCAz, true);
    fYAMLConfig.GetProperty({name, "track", "dca_xy"}, electron_selection.fDCAxy, true);

    //TPC_selection
    AliDHFeCorr::Utils::GetPropertyRange<Float_t>(fYAMLConfig, {name, "PID", "TPC_selection"},
                                                  electron_selection.fTPCNSigmaMin, electron_selection.fTPCNSigmaMax,
                                                  true);

    //TOF Selection
    fYAMLConfig.GetProperty({name, "PID", "require_TOF"}, electron_selection.fRequireTOF, true);
    AliDHFeCorr::Utils::GetPropertyRange<Float_t>(fYAMLConfig, {name, "PID", "TPC_selection"},
                                                  electron_selection.fTOFNSigmaMin, electron_selection.fTOFNSigmaMax,
                                                  true);
    return true;
}

bool AliAnalysisTaskDHFeCorr::ConfigureDMesons(const std::string &name, AliDHFeCorr::AliDMesonSelection &opt_config) {
    std::cout << "Configuring selection for: " << name << std::endl;

    std::string fileName;
    if (!fYAMLConfig.GetProperty({name, "selection", "cut_file"}, fileName, false)) {
        std::cout << "Check Config file. It is not possible to set cut_file for D mesons." << std::endl;
        return false;
    }

    auto fileDMeson = std::unique_ptr<TFile>{TFile::Open(fileName.c_str())};

    std::string cut_name;
    if (!fYAMLConfig.GetProperty({name, "selection", "cut_name"}, cut_name, false)) {
        std::cout << "Check Config file. It is not possible to set cut_name for D mesons." << std::endl;
        return false;
    }

    if (fileDMeson->Get(cut_name.c_str())) {
        auto cuts = dynamic_cast<AliRDHFCuts *>((fileDMeson->Get(cut_name.c_str())));
        opt_config.fDMesonCuts = std::unique_ptr<AliRDHFCuts>(cuts);
        cuts->SaveAs(("Cuts_" + std::string(this->GetName()) + ".root").c_str());
    } else {
        std::cout << "Check Config file. It is not possible to set the cuts for D mesons." << std::endl;
        return false;
    }

    fYAMLConfig.GetProperty({name, "selection", "min_pt"}, opt_config.fPtMin, true);
    fYAMLConfig.GetProperty({name, "selection", "max_pt_pid"}, opt_config.fPtMaxPID, true);

    return true;
}

bool
AliAnalysisTaskDHFeCorr::ConfigureElectronOpt(const std::string &name, AliDHFeCorr::AliConfigureElectronOpt &opt_config,
                                              const std::string &base) {
    std::cout << "Configuring opt for: " << name << std::endl;
    fYAMLConfig.GetProperty({base, name, "pt_bins"}, opt_config.fPtBins, true);
    fYAMLConfig.GetProperty({base, name, "p_bins"}, opt_config.fPBins, true);
    fYAMLConfig.GetProperty({base, name, "n_eta_bins"}, opt_config.fNBinsEta, true);
    fYAMLConfig.GetProperty({base, name, "n_phi_bins"}, opt_config.fNBinsPhi, true);
    fYAMLConfig.GetProperty({base, name, "n_bins_TPC_cls"}, opt_config.fNBinsTPCCls, true);
    fYAMLConfig.GetProperty({base, name, "n_bins_n_sigma"}, opt_config.fNBinsNsigma, true);
    return true;
}

bool AliAnalysisTaskDHFeCorr::ConfigureDMesonOpt(const std::string &name,
                                                 AliDHFeCorr::AliConfigureDMesonOpt &opt_config,
                                                 const std::string &base) {
    std::cout << "Configuring opt for: " << name << std::endl;
    AliDHFeCorr::Utils::GetPropertyRange<Float_t>(fYAMLConfig, {base, name, "inv_mass_range"}, opt_config.fInvMassMin,
                                                  opt_config.fInvMassMax, true);
    fYAMLConfig.GetProperty({base, name, "pt_bins"}, opt_config.fPtBins, true);
    fYAMLConfig.GetProperty({base, name, "n_bins_inv_mass"}, opt_config.fNBinsInvMass, true);
    fYAMLConfig.GetProperty({base, name, "n_eta_bins"}, opt_config.fNBinsEta, true);
    fYAMLConfig.GetProperty({base, name, "n_phi_bins"}, opt_config.fNBinsPhi, true);
    return true;

}

void AliAnalysisTaskDHFeCorr::Terminate(Option_t *) {}

AliDHFeCorr::AliEventQAHistograms AliAnalysisTaskDHFeCorr::CreateQAEvents(const std::string &stage, TList &opt_list) {
    AliDHFeCorr::AliEventQAHistograms qa_hists;

    qa_hists.fNEvents = std::unique_ptr<TH1F>(new TH1F((std::string("NEvents_") + stage).c_str(),
                                                       (std::string("NEvents_") + stage).c_str(),
                                                       10, 0., 10.));
    qa_hists.fVertexZ = std::unique_ptr<TH1F>(new TH1F((std::string("VertexZ_") + stage).c_str(),
                                                       (std::string("VertexZ_") + stage).c_str(),
                                                       120, -15., 15.));
    qa_hists.fCentrality = std::unique_ptr<TH1F>(new TH1F((std::string("Centrality_") + stage).c_str(),
                                                          (std::string("Centrality_") + stage).c_str(),
                                                          120, -10, 110));

    opt_list.Add(qa_hists.fNEvents.get());
    opt_list.Add(qa_hists.fVertexZ.get());
    opt_list.Add(qa_hists.fCentrality.get());

    return qa_hists;
}

AliDHFeCorr::AliDMesonQAHistos
AliAnalysisTaskDHFeCorr::CreateQADMeson(AliDHFeCorr::AliConfigureDMesonOpt config, const std::string &type_particle,
                                        const std::string &stage, TList &opt_list) {
    AliDHFeCorr::AliDMesonQAHistos qa_hists;

    std::vector<Float_t> phi_bins = MakeBins(0.0, 2 * TMath::Pi(), config.fNBinsPhi);
    std::vector<Float_t> eta_bins = MakeBins(-0.8, 0.8, config.fNBinsEta);

    qa_hists.fPtEtaPhi = std::unique_ptr<TH3F>(new TH3F((type_particle + "_PtEtaPhi_" + stage).c_str(),
                                                        (type_particle + "_PtEtaPhi_" + stage +
                                                         ";p_{T};#eta;#varphi").c_str(),
                                                        config.fPtBins.size() - 1, &config.fPtBins[0], config.fNBinsEta,
                                                        &eta_bins[0],
                                                        config.fNBinsPhi, &phi_bins[0]));

    std::vector<Float_t> inv_mass_bins = MakeBins(config.fInvMassMin, config.fInvMassMax,
                                                  config.fNBinsInvMass);

    qa_hists.fPtInvMass = std::unique_ptr<TH2F>(new TH2F((type_particle + "_PtInvMass_" + stage).c_str(),
                                                         (type_particle + "_PtInvMass_" + stage +
                                                          ";p_{T};mass").c_str(),
                                                         config.fPtBins.size() - 1, &config.fPtBins[0],
                                                         config.fNBinsInvMass,
                                                         &inv_mass_bins[0]));

    opt_list.Add(qa_hists.fPtEtaPhi.get());
    opt_list.Add(qa_hists.fPtInvMass.get());

    return qa_hists;
}

AliDHFeCorr::AliElectronQAHistograms AliAnalysisTaskDHFeCorr::CreateQAElectrons(
        AliDHFeCorr::AliConfigureElectronOpt config, const std::string &type_particle, const std::string &stage,
        TList &opt_list) {

    AliDHFeCorr::AliElectronQAHistograms qa_hists;

    std::vector<Float_t> phi_bins = MakeBins(0.0, 2 * TMath::Pi(), config.fNBinsPhi);
    std::vector<Float_t> eta_bins = MakeBins(-0.8, 0.8, config.fNBinsEta);

    qa_hists.fPtEtaPhi = std::unique_ptr<TH3F>(new TH3F((type_particle + "_PtEtaPhi_" + stage).c_str(),
                                                        (type_particle + "_PtEtaPhi_" + stage +
                                                         ";p_{T};#eta;#varphi").c_str(),
                                                        config.fPtBins.size() - 1, &config.fPtBins[0], config.fNBinsEta,
                                                        &eta_bins[0],
                                                        config.fNBinsPhi, &phi_bins[0]));

    qa_hists.fTPCNCls = std::unique_ptr<TH1F>(new TH1F((type_particle + "_TPCNCls_" + stage).c_str(),
                                                       (type_particle + "_TPCNCls_" + stage + ";N cluster TPC").c_str(),
                                                       config.fNBinsTPCCls, 0, 160));

    qa_hists.fTPCNClsDeDx = std::unique_ptr<TH1F>(new TH1F((type_particle + "_TPCNClsDeDx_" + stage).c_str(),
                                                           (type_particle + "_TPCNClsDeDx_" + stage +
                                                            ";N cluster TPC dE/dx").c_str(),
                                                           config.fNBinsTPCCls, 0, 160));

    qa_hists.fITSCls = std::unique_ptr<TH1F>(new TH1F((type_particle + "_ITSCls_" + stage).c_str(),
                                                      (type_particle + "_ITSCls_" + stage + ";N ITS Cls").c_str(), 7, 0,
                                                      7));

    qa_hists.fDCAz = std::unique_ptr<TH1F>(new TH1F((type_particle + "_DCAz_" + stage).c_str(),
                                                    (type_particle + "_DCAz_" + stage + ";DCA z").c_str(), 25, 0, 5));


    qa_hists.fDCAxy = std::unique_ptr<TH1F>(new TH1F((type_particle + "_DCAxy_" + stage).c_str(),
                                                     (type_particle + "_DCAxy_" + stage = "; DCA xy").c_str(), 25, 0,
                                                     5));

    std::vector<Float_t> n_sigma_bins = MakeBins(-10, 10, config.fNBinsNsigma);

    qa_hists.fTPCNsigmaPt = std::unique_ptr<TH2F>(new TH2F((type_particle + "_TPCNsigmaPt_" + stage).c_str(),
                                                           (type_particle + "_fTPCNsigmaPt_" + stage +
                                                            ";p_{T};N#sigma TPC").c_str(),
                                                           config.fPtBins.size() - 1, &config.fPtBins[0],
                                                           config.fNBinsNsigma,
                                                           &n_sigma_bins[0]));

    qa_hists.fTPCNsigmaP = std::unique_ptr<TH2F>(new TH2F((type_particle + "_TPCNsigmaPt_" + stage).c_str(),
                                                          (type_particle + "_fTPCNsigmaPt_" + stage +
                                                           ";p;N#sigma TPC").c_str(),
                                                          config.fPBins.size() - 1, &config.fPBins[0],
                                                          config.fNBinsNsigma, &n_sigma_bins[0]));

    qa_hists.fTOFNsigmaPt = std::unique_ptr<TH2F>(new TH2F((type_particle + "_TOFNsigmaPt_" + stage).c_str(),
                                                           (type_particle + "_TOFNsigmaPt_" + stage +
                                                            ";p_{T};N#sigma TOF").c_str(),
                                                           config.fPBins.size() - 1, &config.fPBins[0],
                                                           config.fNBinsNsigma,
                                                           &n_sigma_bins[0]));

    qa_hists.fTOFNsigmaP = std::unique_ptr<TH2F>(new TH2F((type_particle + "_TOFNsigmaP_" + stage).c_str(),
                                                          (type_particle + "_TOFNsigmaP_" + stage +
                                                           ";p;N#sigma TOF").c_str(),
                                                          config.fPBins.size() - 1, &config.fPBins[0],
                                                          config.fNBinsNsigma, &n_sigma_bins[0]));

    qa_hists.fTPCTOFNSigma = std::unique_ptr<TH2F>(new TH2F((type_particle + "_TPCTOFNSigma_" + stage).c_str(),
                                                            (type_particle + "_TPCTOFNSigma_" + stage +
                                                             ";N#sigma TPC;N#sigma TOF").c_str(),
                                                            config.fNBinsNsigma, &n_sigma_bins[0], config.fNBinsNsigma,
                                                            &n_sigma_bins[0]));

    //Add all the histograms to the opt list
    opt_list.Add(qa_hists.fPtEtaPhi.get());
    opt_list.Add(qa_hists.fTPCNCls.get());
    opt_list.Add(qa_hists.fITSCls.get());
    opt_list.Add(qa_hists.fTPCNClsDeDx.get());
    opt_list.Add(qa_hists.fDCAz.get());
    opt_list.Add(qa_hists.fDCAxy.get());
    opt_list.Add(qa_hists.fTPCNsigmaPt.get());
    opt_list.Add(qa_hists.fTPCNsigmaP.get());
    opt_list.Add(qa_hists.fTOFNsigmaPt.get());
    opt_list.Add(qa_hists.fTOFNsigmaP.get());
    opt_list.Add(qa_hists.fTPCTOFNSigma.get());

    return qa_hists;
}

float AliAnalysisTaskDHFeCorr::GetMaxd0MeasMinusExp(AliAODRecoDecayHF *candidate, float b_field) {
    float d_d0_max = 0;
    auto n_prongs_candidate = static_cast<unsigned int>(candidate->GetNProngs());

    for (unsigned int i = 0; i < n_prongs_candidate; i++) {
        double d0_diff, error_d0_diff;
        candidate->Getd0MeasMinusExpProng(i, b_field, d0_diff, error_d0_diff);
        float norm_dd0 = d0_diff / error_d0_diff;
        if (i == 0)
            d_d0_max = norm_dd0;
        else if (TMath::Abs(norm_dd0) > TMath::Abs(d_d0_max))
            d_d0_max = norm_dd0;
    }
    return d_d0_max;
}

