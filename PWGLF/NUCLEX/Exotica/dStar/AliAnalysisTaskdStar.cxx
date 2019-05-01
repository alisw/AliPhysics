#include "AliAnalysisTaskdStar.h"


// ROOT includes
#include <TChain.h>
#include <TClonesArray.h>
#include <TH2F.h>
#include <TList.h>

// ALIROOT includes
#include "AliPDG.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

// std includes
#include <climits>

///\cond CLASSIMP
ClassImp(AliAnalysisTaskdStar);
///\endcond

using std::string;

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskdStar::AliAnalysisTaskdStar(const char* taskname) : AliAnalysisTaskSE(taskname),
fEventCut(false),
fList(),
fRequireYmin(-0.5f),
fRequireYmax(0.5f),
fDalitPlotMassCutMin(2.340),
fDalitPlotMassCutMax(2.420),
fPID(),
fProduction(),
fReconstructed(),
fTotal(),
fMCDalitzPlot(),
fTree(nullptr),
fDeuteronVector(),
fPiPlusVector(),
fPiMinusVector(),
fMCDeuteronVector(),
fMCPiPlusVector(),
fMCPiMinusVector()
{
  fFilterBit = BIT(5);
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

/// Standard destructor
///
AliAnalysisTaskdStar::~AliAnalysisTaskdStar() {
  if (fList)   delete fList;
  if (fPID)    delete fPID;
  if (fTree)   delete fTree;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskdStar::UserCreateOutputObjects() {

  OpenFile(1);

  fList = new TList();
  fList->SetOwner(true);

  char   letter[2] = {'a','m'};
  string tpctofMC[4] = {"TPC","TPC_TOF","TPC_(TOF)","ITS_sa"};
  float low_mass_limit = 2.2;
  float up_mass_limit = 2.7;

  for (int iC = 0; iC < 2; ++iC) {
    fProduction[iC] = new TH2F(Form("fProduction_dStar(2380)_%c",letter[iC]),";M (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});Counts",50,low_mass_limit,up_mass_limit,20,0,10);
    fList->Add(fProduction[iC]);
    fTotal[iC] = new TH2F(Form("fTotal_dStar(2380)_%c",letter[iC]),";M (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});Counts",50,low_mass_limit,up_mass_limit,20,0,10);
    fList->Add(fTotal[iC]);
    for (int iT = 0; iT < 4; ++iT) {
      fReconstructed[iC][iT] = new TH2F(Form("fRec_dStar(2380)_%c_ITS_%s",letter[iC],tpctofMC[iT].data()),";M (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});Counts",50,low_mass_limit,up_mass_limit,20,0,10);
      fList->Add(fReconstructed[iC][iT]);
    }
  }

  fMCDalitzPlot = new TH2F("dalitzplotMC_dStar", "Dalitz Plot d*(2380) #rightarrow d #pi^{+} #pi^{-} ;#it{M_{inv}}  #pi^{+} #pi^{-}; #it{M_{inv}} d #pi^{-}", 190, 0.05, 1.0, 800, 4.0, 8.0);
  fList->Add(fMCDalitzPlot);

  AliPDG::AddParticlesToPdgDataBase();
  fEventCut.AddQAplotsToList(fList);
  PostData(1,fList);

  OpenFile(2);
  fTree = new TTree("dStarTree", "Data for dStar background analysis");
  fTree->Branch("Deuteron", &fDeuteronVector);
  fTree->Branch("PiPlus", &fPiPlusVector);
  fTree->Branch("PiMinus", &fPiMinusVector);
  fTree->Branch("ZVertex", &fZvtx);
  fTree->Branch("MCDeuteron", &fMCDeuteronVector);
  fTree->Branch("MCPiPlus", &fMCPiPlusVector);
  fTree->Branch("MCPiMinus", &fMCPiMinusVector);
  fTree->SetAutoSave(100000000);
  PostData(2,fTree);

}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskdStar::UserExec(Option_t *) {
  /// The first check performed is on the particle type requested. If this type is AliPID::Unknown -
  /// the default one - the task does not process the information.
  AliVEvent *ev = InputEvent();
  if (!fEventCut.AcceptEvent(ev)) {
    PostData(1, fList);
    return;
  }

  /// To perform the majority of the analysis - and also this one - the standard PID handler is
  /// required.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  fPID = handl->GetPIDResponse();

  TClonesArray *stack = nullptr;
  stack = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (!stack)
  ::Fatal("AliAnalysisTaskdStar::UserExec","MC analysis requested on a sample without the MC particle array.");

  const AliVVertex *vertex = static_cast<const AliVVertex*>(fEventCut.GetPrimaryVertex());
  fZvtx = vertex->GetZ();

  fMCDeuteronVector.clear();
  fMCPiPlusVector.clear();
  fMCPiMinusVector.clear();

  // need comment
  for (int iMC = 0; iMC < stack->GetEntriesFast(); ++iMC) {

    AliAODMCParticle *part = (AliAODMCParticle*)stack->UncheckedAt(iMC);
    const int pdg = TMath::Abs(part->GetPdgCode());
    const int iC = part->Charge() > 0 ? 1 : 0;

    if (pdg == 211) {
      FourVector_t tmp_vec = {(float)part->Pt(),(float)part->Eta(),(float)part->Phi(),(float)part->M()};
      const int m_id = part->GetMother();
      AliAODMCParticle *m_part = (AliAODMCParticle*)stack->UncheckedAt(m_id);
      const int m_pdg = TMath::Abs(m_part->GetPdgCode());
      unsigned char prop = 0u;
      if (iC != 0)                          prop |= c;
      if (part->IsPhysicalPrimary())        prop |= p;
      if (part->IsSecondaryFromMaterial())  prop |= s;
      daughter_struct daug = {0};
      daug.mother_pdg = m_pdg;
      daug.mother_id  = m_id;
      daug.vec        = tmp_vec;
      daug.properties = prop;
      if (iC != 0) {
        prop |= c;
        fMCPiPlusVector.push_back(daug);
      } else {fMCPiMinusVector.push_back(daug);}
    }

    if (pdg == 1000010020) {
      FourVector_t tmp_vec = {(float)part->Pt(),(float)part->Eta(),(float)part->Phi(),(float)part->M()};
      const int m_id = part->GetMother();
      AliAODMCParticle *m_part = (AliAODMCParticle*)stack->UncheckedAt(m_id);
      const int m_pdg = TMath::Abs(m_part->GetPdgCode());
      unsigned char prop = 0u;
      if (iC != 0)                          prop |= c;
      if (part->IsPhysicalPrimary())        prop |= p;
      if (part->IsSecondaryFromMaterial())  prop |= s;
      daughter_struct daug = {0};
      daug.mother_pdg = m_pdg;
      daug.mother_id  = m_id;
      daug.vec        = tmp_vec;
      daug.properties = prop;
      fMCDeuteronVector.push_back(daug);
    }

    if (pdg == 900010020) {
      FourVector_t moth_vec = {0.f,0.f,0.f,0.f};
      for(int iD=0; iD<3; iD++){
        const int daughter_id = part->GetDaughterLabel(0)+iD;
        AliAODMCParticle *daughter_part = (AliAODMCParticle*)stack->At(TMath::Abs(daughter_id));
        FourVector_t tmp_vec = {(float)daughter_part->Pt(),(float)daughter_part->Eta(),(float)daughter_part->Phi(),(float)daughter_part->M()};
        moth_vec+=tmp_vec;
      }
      fProduction[iC]->Fill(part->M(),part->Pt());
      if (part->Y() < fRequireYmin || part->Y() > fRequireYmax) continue;
      fTotal[iC]->Fill(moth_vec.M(),moth_vec.Pt());
    }
  }

  // filling MCDalitzPlot (kinematics limits Mpp2 < 0.255    Mpd2 < 5.01)
  for (const auto& mcdeu : fMCDeuteronVector) {
    unsigned char pdeu = mcdeu.properties;
    if (mcdeu.mother_pdg != 900010020 || !(pdeu & c)) continue;
    for (const auto& mcpim : fMCPiMinusVector) {
      if (mcpim.mother_id != mcdeu.mother_id) continue;
      for (const auto& mcpip : fMCPiPlusVector) {
        if (mcpip.mother_id != mcdeu.mother_id) continue;
        FourVector_t deu_vec = mcdeu.vec;
        FourVector_t pim_vec = mcpim.vec;
        FourVector_t pip_vec = mcpip.vec;
        pip_vec += pim_vec;
        pim_vec += deu_vec;
        deu_vec += pip_vec;
        if (deu_vec.M() < fDalitPlotMassCutMin || deu_vec.M() > fDalitPlotMassCutMax) continue;
        fMCDalitzPlot->Fill(pip_vec.M2(), pim_vec.M2());
      }
    }
  }

  // Checking how many dstar in acceptance are reconstructed well
  std::vector<mother_struct> mothers;
  mothers.reserve(40);

  fDeuteronVector.clear();
  fPiPlusVector.clear();
  fPiMinusVector.clear();

  for (Int_t iT = 0; iT < (Int_t)ev->GetNumberOfTracks(); ++iT) {

    AliAODTrack *track = static_cast<AliAODTrack*>(ev->GetTrack(iT));
    if (track->GetID() < 0) continue;
    // if (!track->TestFilterBit(fFilterBit)) continue;

    bool fb = track->TestFilterBit(fFilterBit);
    bool itsSA = track->GetStatus() & AliVTrack::kITSrefit && (track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1));
    if (!fb && !itsSA) continue;

    AliAODMCParticle *part = (AliAODMCParticle*)stack->At(TMath::Abs(track->GetLabel()));
    if (!part) continue;
    const int pdg = TMath::Abs(part->GetPdgCode());
    const int mother_id = part->GetMother();
    AliAODMCParticle* mother = (mother_id >= 0) ? (AliAODMCParticle*)stack->At(mother_id) : nullptr;
    if (!mother) continue;
    const int mum_pdg = TMath::Abs(mother->GetPdgCode());   // IMPO!! for A/M analysis is necessary to take also the charge


    // add deuterons and pions to the Tree for background analysis (ITS TPC only)
    // if they are under 3 sigmas TPC response

    // PID for deuterons
    if (fb && TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kDeuteron)) < 3.) {
      FourVector_t tmp_deu = {(float)track->Pt(), (float)track->Eta(), (float)track->Phi(), (float)track->M(AliAODTrack::kDeuteron)};
      unsigned char prop = 0u;
      if (part->Charge() > 0)               prop |= c;
      if (part->IsPhysicalPrimary())        prop |= p;
      if (part->IsSecondaryFromMaterial())  prop |= s;
      if (AliAnalysisTaskdStar::HasTOF(track) && (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kDeuteron)) < 3.)) prop |= t;
      daughter_struct deu;
      deu.mother_pdg = mum_pdg;
      deu.mother_id  = mother_id;
      deu.mc_truth   = pdg;
      deu.vec        = tmp_deu;
      deu.properties = prop;
      fDeuteronVector.push_back(deu);
    }

    // PID for Pions
    if (fb) {    // TPC PID if available
      if (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kPion)) < 3.) {
        FourVector_t tmp_pi = {(float)track->Pt(), (float)track->Eta(), (float)track->Phi(), (float)track->M(AliAODTrack::kPion)};
        unsigned char prop = 0u;
        if (part->Charge() > 0)              prop |= c;
        if (part->IsPhysicalPrimary())       prop |= p;
        if (part->IsSecondaryFromMaterial()) prop |= s;
        if (AliAnalysisTaskdStar::HasTOF(track) && (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kPion)) < 3.)) prop |= t;
        daughter_struct pi;
        pi.mother_pdg = mum_pdg;
        pi.mother_id = mother_id;
        pi.mc_truth = pdg;
        pi.vec = tmp_pi;
        pi.properties = prop;
        track->Charge() > 0 ? fPiPlusVector.push_back(pi) : fPiMinusVector.push_back(pi);
      }
    } else if (itsSA) { // ITS standalone PID if TPC is not available
      if (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kITS, track, AliPID::kPion)) < 3.) {
        FourVector_t tmp_pi = {(float)track->Pt(), (float)track->Eta(), (float)track->Phi(), (float)track->M(AliAODTrack::kPion)};
        unsigned char prop = 0u;
        if (part->Charge() > 0)              prop |= c;
        if (part->IsPhysicalPrimary())       prop |= p;
        if (part->IsSecondaryFromMaterial()) prop |= s;
        prop |= i;
        daughter_struct pi;
        pi.mother_pdg = mum_pdg;
        pi.mother_id = mother_id;
        pi.mc_truth = pdg;
        pi.vec = tmp_pi;
        pi.properties = prop;
        track->Charge() > 0 ? fPiPlusVector.push_back(pi) : fPiMinusVector.push_back(pi);
     }
    }


    // Check wheter the track belongs to a deuteron
    if (pdg == 1000010020 && mum_pdg == 900010020 && fb) {
      if(TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kDeuteron)) > 3.) continue;
      FourVector_t tmp_vec = {(float)track->Pt(), (float)track->Eta(), (float)track->Phi(), (float)track->M(AliAODTrack::kDeuteron)};
      auto it = std::find(mothers.begin(), mothers.end(), mother_id);
      if (it == mothers.end()){
        mother_struct tmp_mum;
        tmp_mum.id = mother_id;
        tmp_mum.deuteron_tof = AliAnalysisTaskdStar::HasTOF(track) && (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kDeuteron)) < 3.);
        tmp_mum.n_daughters = 1;
        tmp_mum.vec = tmp_vec;
        mothers.push_back(tmp_mum);
      }
      else{
        it->n_daughters++;
        it->deuteron_tof *= AliAnalysisTaskdStar::HasTOF(track) && (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kDeuteron)) < 3.);
        it->vec+=tmp_vec;
      }
    }

    // Check wether the track belgons to a pion
    if (pdg == 211 && mum_pdg == 900010020) {
      if (fb) {
        if (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kPion)) > 3.) continue;
        FourVector_t tmp_vec = {(float)track->Pt(), (float)track->Eta(), (float)track->Phi(), (float)track->M(AliAODTrack::kPion)};
        auto it = std::find(mothers.begin(), mothers.end(), mother_id);
        if (it == mothers.end()) {
          mother_struct tmp_mum;
          tmp_mum.id = mother_id;
          tmp_mum.pi_tpc = true;
          tmp_mum.pi_tof = AliAnalysisTaskdStar::HasTOF(track) && (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kPion)) < 3.);
          tmp_mum.n_daughters = 1;
          tmp_mum.vec = tmp_vec;
          mothers.push_back(tmp_mum);
        } else {
          it->n_daughters++;
          it->pi_tof *= AliAnalysisTaskdStar::HasTOF(track) && (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kPion)) < 3.);
          it->pi_tpc = true;
          it->vec += tmp_vec;
        }
      } else if (itsSA) {
        if (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kITS, track, AliPID::kPion)) > 3.) continue;
        FourVector_t tmp_vec = {(float)track->Pt(), (float)track->Eta(), (float)track->Phi(), (float)track->M(AliAODTrack::kPion)};
        auto it = std::find(mothers.begin(), mothers.end(), mother_id);
        if (it == mothers.end()) {
          mother_struct tmp_mum;
          tmp_mum.id = mother_id;
          tmp_mum.pi_tof = false;
          tmp_mum.pi_tpc = false;
          tmp_mum.n_daughters = 1;
          tmp_mum.vec = tmp_vec;
          mothers.push_back(tmp_mum);
        } else {
          it->n_daughters++;
          it->pi_tof = false;
          it->pi_tpc = false;
          it->vec += tmp_vec;
        }
      }
    }

  } // End AOD track loop

  // Filling three
  fTree->Fill();

  // Filling histograms
  for (const auto& mum : mothers) {
    if (mum.n_daughters != 3) continue;
    if (mum.vec.Rapidity() < fRequireYmin || mum.vec.Rapidity() > fRequireYmax) continue;
    AliAODMCParticle *part = static_cast<AliAODMCParticle*>(stack->At(mum.id));
    const int iC = part->Charge() > 0 ? 1 : 0;
    const float pt_rec = mum.vec.Pt();
    const float mass_rec = mum.vec.M();
    if (mum.pi_tpc) fReconstructed[iC][0]->Fill(mass_rec, pt_rec);
    if (mum.pi_tof && mum.deuteron_tof && mum.pi_tpc) fReconstructed[iC][1]->Fill(mass_rec, pt_rec);
    if (mum.deuteron_tof && mum.pi_tpc) fReconstructed[iC][2]->Fill(mass_rec, pt_rec);
    fReconstructed[iC][3]->Fill(mass_rec, pt_rec);
  }

  //  Post output data.
  PostData(1,fList);
  PostData(2,fTree);
}

/// This function checks whether a track has or has not a prolongation in the TOF.
///
/// \param track Track that has to be checked
/// \return true if the track has a matching hit in the TOF.
///
bool AliAnalysisTaskdStar::HasTOF(AliVTrack *track) {
  bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  return hasTOFout && hasTOFtime && (track->GetIntegratedLength() > 350.);
}
