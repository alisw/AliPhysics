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

    DefineOutput(1, TTree::Class());
    DefineOutput(2, TTree::Class());
    DefineOutput(3, TTree::Class());
    DefineOutput(4, TList::Class());
    DefineOutput(5, TList::Class());
    DefineOutput(6, TList::Class());
    DefineOutput(7, TTree::Class());
    DefineOutput(8, TTree::Class());

}

//_____________________________________________________________________________
void AliAnalysisTaskDHFeCorr::UserCreateOutputObjects() {
    //Additional setup for the task that needs to be done at the time of run time
    fDMesonRequirements.fDMesonCuts->GetPidHF()->SetPidResponse(fInputHandler->GetPIDResponse());

    //Remove the trigger mask from the automatic cuts
    fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAny);

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

        fElectronTreeMC = std::unique_ptr<TTree>(new TTree("electron_mc", "electron_mc"));
        fDmesonTreeMC = std::unique_ptr<TTree>(new TTree("dmeson_mc", "dmeson_mc"));

        AddMCTreeVariables(fDmesonTreeMC);
        AddMCTreeVariables(fElectronTreeMC);
    }

    fEventCuts.AddQAplotsToList(&fOptEvent);

    fElectronQABeforeCuts = CreateQAElectrons(fElectronOptConfig, "electron", "before", fOptElectron);
    fElectronQAAfterTrackCuts = CreateQAElectrons(fElectronOptConfig, "electron", "after_track_cuts", fOptElectron);
    fElectronQAAfterCuts = CreateQAElectrons(fElectronOptConfig, "electron", "after", fOptElectron);

    fDMesonQABeforeCuts = CreateQADMeson(fDMesonOptConfig, "dmeson", "before", fOptDMeson);
    fDMesonQAAfterCuts = CreateQADMeson(fDMesonOptConfig, "dmeson", "after", fOptDMeson);

    PostOutput();
}


void AliAnalysisTaskDHFeCorr::UserExec(Option_t *) {

    CheckConfiguration();
    SetRunAndEventNumber();

    if (!InputEvent()) {
        PostOutput();
        return;
    }

    if (!fEventCuts.AcceptEvent(InputEvent())) {
        PostOutput();
        return;
    }

    if (fSaveEvent) {
        fEventInfo = EventInfo();
        fEventTree->Fill();
    }

    if (fProcessDMeson)
        DMesonAnalysis();

    if (fProcessElectron)
        ElectronAnalysis();

    //Post information to the output
    PostOutput();

}

AliDHFeCorr::AliEvent AliAnalysisTaskDHFeCorr::EventInfo() {
    AliDHFeCorr::AliEvent event;

    event.fRunNumber = fRunNumber;
    event.fEventNumber = fEventNumber;

    event.fVtxZ = fEventCuts.GetPrimaryVertex()->GetZ();

    auto multiplicity_selection = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
    if (multiplicity_selection) {
        event.fMultV0MPercentile = multiplicity_selection->GetMultiplicityPercentile("V0M");
        event.fMultiRefMult08Percentile = multiplicity_selection->GetMultiplicityPercentile("RefMult08");
        event.fMultiSPDTrackletsPercentile = multiplicity_selection->GetMultiplicityPercentile("SPDTracklets");

        if (multiplicity_selection->GetEstimator("V0M"))
            event.fMultV0M = multiplicity_selection->GetEstimator("V0M")->GetValue();
        if (multiplicity_selection->GetEstimator("RefMult08"))
            event.fMultiRefMult08 = multiplicity_selection->GetEstimator("RefMult08")->GetValue();
        if (multiplicity_selection->GetEstimator("SPDTracklets"))
            event.fMultiSPDTracklets = multiplicity_selection->GetEstimator("SPDTracklets")->GetValue();
    }

    return event;
}


std::vector<AliDHFeCorr::AliElectron> AliAnalysisTaskDHFeCorr::ElectronAnalysis() {

    // Reconstructed candidates
    auto aod_event = dynamic_cast<AliAODEvent *>(InputEvent());

    std::vector<AliDHFeCorr::AliElectron> tracks;
    tracks.reserve(static_cast<unsigned long>(aod_event->GetNumberOfTracks()));

    for (int i(0); i < aod_event->GetNumberOfTracks(); i++) {
        tracks.emplace_back(dynamic_cast<AliAODTrack *>(aod_event->GetTrack(i)), fRunNumber,
                            fEventNumber, aod_event, fInputHandler->GetPIDResponse());
    }

    auto selected_tracks = FilterElectronsTracking(fElectronRequirements, tracks);
    auto selected_electrons = FilterElectronsPID(fElectronRequirements, selected_tracks);

    auto selected_partner_tracks = FilterElectronsTracking(fPartnerElectronRequirements, tracks);
    auto partner_electrons = FilterElectronsPID(fPartnerElectronRequirements, selected_partner_tracks);

    if (!fReducedElectronInfo) {
        for (auto &electron: selected_electrons)
            FindNonHFe(electron, partner_electrons);
    }

    if (fSaveHistograms) {
        //fill the electron (track) QA before applying the cuts
        FillElectronQA(tracks, fElectronQABeforeCuts);
        //fill the electron QA after track selection
        FillElectronQA(selected_tracks, fElectronQAAfterTrackCuts);
        //fill the electron QA after applying the cuts
        FillElectronQA(selected_electrons, fElectronQAAfterCuts);
    }

    if (fIsMC) {
        FillAllElectronsMCInfo(selected_electrons);

        auto mc_particles = FillMCParticleInfo();
        auto hfe_mc = FilterHFeInMCParticles(mc_particles);

        FillTreeFromStdContainer(hfe_mc, &fMCParticle, fElectronTreeMC);

        selected_electrons = selected_tracks; //Remove PID Cuts from MC
    }

    FillTreeFromStdContainer(selected_electrons, &fElectron, fElectronTree);

    return selected_electrons;
}

std::vector<AliDHFeCorr::AliParticleMC> AliAnalysisTaskDHFeCorr::FillMCParticleInfo() {

    const auto mc_information = dynamic_cast<TClonesArray *>(InputEvent()->GetList()->FindObject(
            AliAODMCParticle::StdBranchName()));

    std::vector<AliDHFeCorr::AliParticleMC> mc_particles;

    for (int i(0); i < mc_information->GetEntriesFast(); i++) {

        auto particle = dynamic_cast<AliAODMCParticle *>(mc_information->At(i));

        if (!particle)
            continue;

        auto origin = AliVertexingHFUtils::CheckOrigin(mc_information, particle, false);

        mc_particles.emplace_back(particle, fRunNumber, fEventNumber, i, origin);
    }

    return mc_particles;
}


std::vector<AliDHFeCorr::AliParticleMC> AliAnalysisTaskDHFeCorr::FindHFParticleInMC(int pdg,
                                                                                    vector<AliDHFeCorr::AliParticleMC> &mc_particles) {
    pdg = abs(pdg);

    std::vector<AliDHFeCorr::AliParticleMC> selected_particles;

    if (mc_particles.empty()) {
        mc_particles = FillMCParticleInfo();
    }

    for (auto &particle: mc_particles) {
        if (abs(particle.fPDGCode) == pdg)
            selected_particles.push_back(particle);
    }

    return selected_particles;
}

std::vector<AliDHFeCorr::AliParticleMC> AliAnalysisTaskDHFeCorr::FilterHFeInMCParticles(
        std::vector<AliDHFeCorr::AliParticleMC> &electrons) {

    std::vector<AliDHFeCorr::AliParticleMC> hfe;
    hfe.reserve(electrons.size());

    const auto mc_information = dynamic_cast<TClonesArray *>(InputEvent()->GetList()->FindObject(
            AliAODMCParticle::StdBranchName()));

    for (auto &mc_electron: electrons) {
        if (IsHFe(mc_electron.fMCParticle, mc_information))
            hfe.push_back(mc_electron);
    }

    hfe.shrink_to_fit();

    return hfe;
}

bool AliAnalysisTaskDHFeCorr::IsHFe(AliAODMCParticle *particle, TClonesArray *mc_information) {

    if (abs(particle->PdgCode()) != 11)
        return false;

    if (particle->GetMother() < 0)
        return false;

    if (!mc_information)
        throw std::runtime_error("The program is trying to access MC information, but the mc info is not present.");

    auto mother = dynamic_cast<AliAODMCParticle *>(mc_information->At(particle->GetMother()));

    if (!mother)
        return false;

    auto abs_pdg_mother = abs(mother->PdgCode());

    return (abs_pdg_mother > 400 && abs_pdg_mother < 600) || (abs_pdg_mother > 4000 && abs_pdg_mother < 6000);

}


std::vector<AliDHFeCorr::AliDMeson> AliAnalysisTaskDHFeCorr::DMesonAnalysis() {
    auto aod_event = dynamic_cast<AliAODEvent *>(InputEvent());

    //Move D meson candidates to vectors
    const TClonesArray *dmeson_candidates = dynamic_cast<TClonesArray *>(aod_event->GetList()->FindObject(
            fgkDMesonListName.at(fDmesonSpecies).c_str()));

    auto d_mesons = FillDMesonInfo(dmeson_candidates, fDmesonSpecies, fDMesonRequirements);
    //Select D mesons
    auto selected_d_mesons = FilterDmesons(d_mesons, fDMesonRequirements, aod_event, fgkDMesonPDG.at(fDmesonSpecies));

    if (fSaveHistograms) {
        FillDmesonQA(d_mesons, fDMesonQABeforeCuts, fDmesonSpecies);
        FillDmesonQA(selected_d_mesons, fDMesonQAAfterCuts, fDmesonSpecies);
    }

    FillDmesonMCInfo(selected_d_mesons, fDmesonSpecies);

    if (fIsMC) {
        auto mc_particles = vector<AliDHFeCorr::AliParticleMC>();

        auto d_mc = FindHFParticleInMC(fgkDMesonPDG.at(fDmesonSpecies), mc_particles);

        FillTreeFromStdContainer(d_mc, &fMCParticle, fDmesonTreeMC);

        if (fIsEffMode)
            selected_d_mesons = FilterTrueDMesons(selected_d_mesons);
    }

    FillTreeFromStdContainer(selected_d_mesons, &fDmeson, fDmesonTree);


    return selected_d_mesons;
}

//_____________________________________________________________________________


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

        if (candidate.fNCrossedRowsTPC < electronSelection.fNCrossedRowsTPCMin)
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
        selected_electrons.push_back(candidate);
    }

    selected_electrons.shrink_to_fit();
    return selected_electrons;
}

void AliAnalysisTaskDHFeCorr::FillAllElectronsMCInfo(std::vector<AliDHFeCorr::AliElectron> &electrons) {
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
        candidate.fPDGCode = mc_part->PdgCode();
        candidate.fOrigin = AliVertexingHFUtils::CheckOrigin(mc_information, mc_part, false);
        candidate.fLabel = label;

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
        candidate.fLabel = label;

        const auto mc_part = dynamic_cast<AliAODMCParticle *>(mc_information->At(label));
        candidate.fPtMC = mc_part->Pt();
        const auto origin = AliVertexingHFUtils::CheckOrigin(mc_information, mc_part, false); //Prompt = 4, FeedDown = 5
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

    //Long64_t ev_number = Entry();
    Long64_t ev_number = fInputHandler->GetReadEntry();

    fEventNumber = (unsigned int) ev_number + (unsigned int) (fDirNum << 17);
}


void AliAnalysisTaskDHFeCorr::PostOutput() {
    if (fProcessElectron) {
        PostData(1, fElectronTree.get());
    }
    if (fProcessDMeson) {
        PostData(2, fDmesonTree.get());
    }

    if (fSaveEvent) {
        PostData(3, fEventTree.get());
    }

    if (fSaveHistograms) {
        PostData(4, &fOptEvent);
        PostData(5, &fOptElectron);
        PostData(6, &fOptDMeson);
    }

    if (fIsMC) {
        if (fProcessElectron)
            PostData(7, fElectronTreeMC.get());
        if (fProcessDMeson)
            PostData(8, fDmesonTreeMC.get());
    }
}


void AliAnalysisTaskDHFeCorr::FillElectronQA(const std::vector<AliDHFeCorr::AliElectron> &tracks,
                                             AliDHFeCorr::AliElectronQAHistograms &histograms) {
    for (const auto &electron: tracks) {
        histograms.fPtEtaPhi->Fill(electron.fPt, electron.fEta, electron.fPhi);
        histograms.fTPCNCls->Fill(electron.fNCrossedRowsTPC);
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

void AliAnalysisTaskDHFeCorr::FindNonHFe(AliDHFeCorr::AliElectron &main_electron,
                                         const std::vector<AliDHFeCorr::AliElectron> &partners) const {

    std::vector<Float_t> inv_mass_uls, inv_mass_ls, pt_uls, pt_ls;
    std::vector<UShort_t> n_cross_row_uls, n_cross_row_ls;
    std::vector<UInt_t> id_uls, id_ls;

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
            pt_ls.push_back(track_partner->Pt());
            n_cross_row_ls.push_back(static_cast<UShort_t >(track_partner->GetTPCNCrossedRows()));
            inv_mass_ls.push_back(pair_mass);
            id_ls.push_back(static_cast<unsigned int &&>(TMath::Abs(track_partner->GetID())));

        } else {
            pt_uls.push_back(track_partner->Pt());
            n_cross_row_uls.push_back(static_cast<UShort_t >(track_partner->GetTPCNCrossedRows()));
            inv_mass_uls.push_back(pair_mass);
            id_uls.push_back(static_cast<unsigned int &&>(TMath::Abs(track_partner->GetID())));
        }
    }

    main_electron.fInvMassPartnersULS = inv_mass_uls;
    main_electron.fInvMassPartnersLS = inv_mass_ls;

    main_electron.fPartnersULSID = id_uls;
    main_electron.fPartnersLSID = id_ls;

    main_electron.fCrossedRowsTPCPartnersLS = n_cross_row_ls;
    main_electron.fCrossedRowsTPCPartnersULS = n_cross_row_uls;
    main_electron.fPtPartnersULS = pt_uls;
    main_electron.fPtPartnersLS = pt_ls;

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

        //Preselect to save  time
        TObjArray arrTracks(cand.fRecoObj->GetNProngs());
        for (Int_t itrack = 0; itrack < cand.fRecoObj->GetNProngs(); itrack++) {
            AliAODTrack *tr = vertexing_HF.GetProng(aod_event, cand.fRecoObj, itrack);
            arrTracks.AddAt(tr, itrack);
        }

        if (!dmeson_selection.fDMesonCuts->PreSelect(arrTracks)) {
            continue;
        }

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
            }
        }
        //avoid a lot of computational time by rejecting minimum Pt here
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
            AliDHFeCorr::AliDMeson reflection = BuildReflection(cand);
            d_mesons.push_back(reflection);
        }

        dmeson_selection.fDMesonCuts->CleanOwnPrimaryVtx(reco_cand, aod_event, original_own_vertex);
    }

    d_mesons.shrink_to_fit();
    return d_mesons;
}

AliDHFeCorr::AliDMeson AliAnalysisTaskDHFeCorr::BuildReflection(const AliDHFeCorr::AliDMeson &cand) const {
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
    return reflection;
}


std::vector<AliDHFeCorr::AliDMeson>
AliAnalysisTaskDHFeCorr::FilterDmesons(const std::vector<AliDHFeCorr::AliDMeson> &dmeson_candidates,
                                       const AliDHFeCorr::AliDMesonSelection &selectionDMeson, AliAODEvent *aod_event,
                                       int pdg_dmeson) {

    std::vector<AliDHFeCorr::AliDMeson> d_mesons;
    d_mesons.reserve(dmeson_candidates.size());

    for (const auto &candidate: dmeson_candidates) {
        const auto d_candidate = candidate.fRecoObj;
        // Enforce that the PID has the original value
        selectionDMeson.fDMesonCuts->SetUsePID(selectionDMeson.fUsePID);

        //Check if it is in fiducial acceptance
        if (!selectionDMeson.fDMesonCuts->IsInFiducialAcceptance(d_candidate->Pt(), d_candidate->Y(pdg_dmeson)))
            continue;

        //Use PID only for values smaller than fPtMaxPID, it will be restored to the original value after the
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
        AliDHFeCorr::AliDMeson cand(candidate);
        cand.fSelectionStatusDefaultPID = pid_selection_status;

        d_mesons.push_back(cand);
    }
    d_mesons.shrink_to_fit();

    return d_mesons;
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

    return hasHitFirstLayer && !hasHitSecondLayer && (requirement == kExclusiveFirst);
}

AliAnalysisTaskDHFeCorr *AliAnalysisTaskDHFeCorr::AddTask(const std::string &name,
                                                          const std::string &config_file, Int_t trigger) {
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
    task->SelectCollisionCandidates(trigger);

    mgr->AddTask(task);

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

    std::string fileName = static_cast<std::string>(AliAnalysisManager::GetCommonFileName());
    fileName += ":DHFeCorrelation_" + name;

    mgr->ConnectOutput(task, 1,
                       mgr->CreateContainer("ElectronTree", TTree::Class(), AliAnalysisManager::kOutputContainer,
                                            fileName.c_str()));

    mgr->ConnectOutput(task, 2, mgr->CreateContainer("DMesonTree", TTree::Class(), AliAnalysisManager::kOutputContainer,
                                                     fileName.c_str()));

    mgr->ConnectOutput(task, 3, mgr->CreateContainer("EventTree", TTree::Class(), AliAnalysisManager::kOutputContainer,
                                                     fileName.c_str()));

    mgr->ConnectOutput(task, 4, mgr->CreateContainer("EventsQA", TList::Class(), AliAnalysisManager::kOutputContainer,
                                                     fileName.c_str()));

    mgr->ConnectOutput(task, 5, mgr->CreateContainer("ElectronQA", TList::Class(), AliAnalysisManager::kOutputContainer,
                                                     fileName.c_str()));

    mgr->ConnectOutput(task, 6, mgr->CreateContainer("DMesonQA", TList::Class(), AliAnalysisManager::kOutputContainer,
                                                     fileName.c_str()));

    mgr->ConnectOutput(task, 7,
                       mgr->CreateContainer("ElectronTreeMC", TTree::Class(), AliAnalysisManager::kOutputContainer,
                                            fileName.c_str()));

    mgr->ConnectOutput(task, 8,
                       mgr->CreateContainer("DMesonTreeMC", TTree::Class(), AliAnalysisManager::kOutputContainer,
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
    fYAMLConfig.GetProperty("reduced_electron_info", fReducedElectronInfo, true);

    fYAMLConfig.GetProperty("process_electron", fProcessElectron, true);
    fYAMLConfig.GetProperty("process_dmeson", fProcessDMeson, true);
    fYAMLConfig.GetProperty("save_event", fSaveEvent, true);


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

    //min_TPC_crossedrows
    fYAMLConfig.GetProperty({name, "track", "min_TPC_crossedrows"}, electron_selection.fNCrossedRowsTPCMin, true);
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
    opt_config.fUsePID = opt_config.fDMesonCuts->GetIsUsePID();

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

void AliAnalysisTaskDHFeCorr::AddEventVariables(std::unique_ptr<TTree> &tree) {
    tree->Branch("RunNumber", &fEventInfo.fRunNumber);
    tree->Branch("EventNumber", &fEventInfo.fEventNumber);
    tree->Branch("VtxZ", &fEventInfo.fVtxZ);
    tree->Branch("MultV0M", &fEventInfo.fMultV0M);
    tree->Branch("MultiRefMult08", &fEventInfo.fMultiRefMult08);
    tree->Branch("MultiSPDTracklets", &fEventInfo.fMultiSPDTracklets);
    tree->Branch("MultV0MPercentile", &fEventInfo.fMultV0MPercentile);
    tree->Branch("MultiRefMult08Percentile", &fEventInfo.fMultiRefMult08Percentile);
    tree->Branch("MultiSPDTrackletsPercentile", &fEventInfo.fMultiSPDTrackletsPercentile);
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
        tree->Branch("NCrossedRowsTPC", &fElectron.fNCrossedRowsTPC);
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
        tree->Branch("PtPartnersULS", &fElectron.fPtPartnersULS);
        tree->Branch("PtPartnersLS", &fElectron.fPtPartnersLS);
        tree->Branch("CrossedRowsTPCPartnersULS", &fElectron.fCrossedRowsTPCPartnersULS);
        tree->Branch("CrossedRowsTPCPartnersLS", &fElectron.fCrossedRowsTPCPartnersLS);
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
    tree->Branch("Label", &fElectron.fLabel);
    tree->Branch("Origin", &fElectron.fOrigin);
    tree->Branch("PDGCode", &fElectron.fPDGCode);
    tree->Branch("FirstMotherPDG", &fElectron.fFirstMotherPDG);
    tree->Branch("FirstMotherPt", &fElectron.fFirstMotherPt);
    tree->Branch("SecondMotherPDG", &fElectron.fSecondMotherPDG);
    tree->Branch("SecondMotherPt", &fElectron.fSecondMotherPt);
}

void AliAnalysisTaskDHFeCorr::AddMCTreeVariables(std::unique_ptr<TTree> &tree) {
    tree->Branch("RunNumber", &fMCParticle.fRunNumber);
    tree->Branch("EventNumber", &fMCParticle.fEventNumber);
    tree->Branch("Label", &fMCParticle.fLabel);
    tree->Branch("E", &fMCParticle.fE);
    tree->Branch("Pt", &fMCParticle.fPt);
    tree->Branch("Eta", &fMCParticle.fEta);
    tree->Branch("Phi", &fMCParticle.fPhi);
    tree->Branch("Xv", &fMCParticle.fXv);
    tree->Branch("Yv", &fMCParticle.fYv);
    tree->Branch("Zv", &fMCParticle.fZv);
    tree->Branch("Tv", &fMCParticle.fTv);
    tree->Branch("Charge", &fMCParticle.fCharge);
    tree->Branch("PDGCode", &fMCParticle.fPDGCode);
    tree->Branch("Origin", &fMCParticle.fOrigin);
}

void AliAnalysisTaskDHFeCorr::AddDMesonMCVariables(std::unique_ptr<TTree> &tree) {
    tree->Branch("PtMC", &fDmeson.fPtMC);
    tree->Branch("Label", &fDmeson.fLabel);
    tree->Branch("IsD", &fDmeson.fIsD);
    tree->Branch("IsParticle", &fDmeson.fIsParticle);
    tree->Branch("IsPrompt", &fDmeson.fIsPrompt);
}

