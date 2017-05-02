#include "AliAnalysisTaskPsEfficiency.h"

// ROOT includes
#include <TChain.h>
#include <TClonesArray.h>
#include <TH2F.h>
#include <TList.h>
#include <TMath.h>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
using namespace ROOT::Math;
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> FourVector_t ;
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
ClassImp(AliAnalysisTaskPsEfficiency);
///\endcond

struct mother_struct{
  int id;
  bool phi_tof;
  bool proton_tof;
  int n_daughters;
  FourVector_t vec;
  bool operator==(const int &id_comp) const {return id==id_comp;}
};

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskPsEfficiency::AliAnalysisTaskPsEfficiency(const char* taskname) : AliAnalysisTaskSE(taskname),
  fEventCut(false),
  fList(),
  fProduction(),
  fReconstructed(),
  fTotal(),
  fRequireYmin(-0.5f),
  fRequireYmax(0.5f),
  fPID()
{
  fFilterBit = BIT(8);
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

/// Standard destructor
///
AliAnalysisTaskPsEfficiency::~AliAnalysisTaskPsEfficiency() {
  if (fList) delete fList;
  if (fPID) delete fPID;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskPsEfficiency::UserCreateOutputObjects() {

  fList = new TList();
  fList->SetOwner(true);

  char   letter[2] = {'A','M'};
  string tpctof[2] = {"TPC","TOF"};
  string tpctofMC[3] = {"TPC","TPC_TOF","TPC_(TOF)"};
  string PsStates[2] = {"Ps(2100)","Ps(2500)"};
  float low_mass_limit[2] = {1.9,2.3};
  float up_mass_limit[2] = {2.3,2.7};

  for(int iS=0; iS<2; iS++){
    for (int iC = 0; iC < 2; ++iC) {

      fProduction[iS][iC] = new TH2F(Form("fProduction_%s_%c",PsStates[iS].data(),letter[iC]),";M (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});Counts",20,low_mass_limit[iS],up_mass_limit[iS],20,0,10);
      fList->Add(fProduction[iS][iC]);

      fTotal[iS][iC] = new TH2F(Form("fTotal_%s_%c",PsStates[iS].data(),letter[iC]),";M (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});Counts",20,low_mass_limit[iS],up_mass_limit[iS],20,0,10);
      fList->Add(fTotal[iS][iC]);

      for (int iT = 0; iT < 3; ++iT) {
        fReconstructed[iS][iC][iT] = new TH2F(Form("fRec_%s_%c_ITS_%s",PsStates[iS].data(),letter[iC],tpctofMC[iT].data()),";M (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});Counts",20,low_mass_limit[iS],up_mass_limit[iS],20,0,10);
        fList->Add(fReconstructed[iS][iC][iT]);
      }
    }
  }

  AliPDG::AddParticlesToPdgDataBase();
  fEventCut.AddQAplotsToList(fList);

  PostData(1,fList);
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskPsEfficiency::UserExec(Option_t *){
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
  // get branch "mcparticles"
  stack = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (!stack)
    ::Fatal("AliAnalysisTaskPsEfficiency::UserExec","MC analysis requested on a sample without the MC particle array.");

  /// Making the list of the pentaquarks we want to measure
  //std::cout << "LISTING ALL THE MOTHERS: "<<std::endl<<std::endl;
  int numero=0;
  for (int iMC = 0; iMC < stack->GetEntriesFast(); ++iMC) {
    AliAODMCParticle *part = (AliAODMCParticle*)stack->UncheckedAt(iMC);
    const int pdg = TMath::Abs(part->GetPdgCode());
    const int iC = part->Charge() > 0 ? 1 : 0;
    const int mult = -1 + 2 * iC;
    const int iS = pdg == 9322134 ? 0 : pdg == 9322136 ? 1 : -1;
    if (iS == -1) continue;
    FourVector_t moth_vec,tmp_vec;
    bool visible_dacay = false;
    moth_vec.SetPt(0.f);
    moth_vec.SetEta(0.f);
    moth_vec.SetPhi(0.f);
    moth_vec.SetM(0.f);
    const int phi_id = part->GetDaughter(0); //Phi meson
    AliAODMCParticle *phi_part = (AliAODMCParticle*)stack->At(TMath::Abs(phi_id));
    int phi_pdg = phi_part->GetPdgCode();
    int phi_dah_n = phi_part->GetNDaughters();
    const int kaon1_id = phi_part->GetDaughter(0);
    AliAODMCParticle *kaon1_part = (AliAODMCParticle*)stack->At(TMath::Abs(kaon1_id));
    int kaon1_pdg = TMath::Abs(kaon1_part->GetPdgCode());
    if(kaon1_pdg==321) visible_dacay = true;
    tmp_vec.SetPt(kaon1_part->Pt());
    tmp_vec.SetEta(kaon1_part->Eta());
    tmp_vec.SetPhi(kaon1_part->Phi());
    tmp_vec.SetM(kaon1_part->M());
    moth_vec+=tmp_vec;
    const int kaon2_id = phi_part->GetDaughter(1);
    AliAODMCParticle *kaon2_part = (AliAODMCParticle*)stack->At(TMath::Abs(kaon2_id));
    int kaon2_pdg = TMath::Abs(kaon2_part->GetPdgCode());
    tmp_vec.SetPt(kaon2_part->Pt());
    tmp_vec.SetEta(kaon2_part->Eta());
    tmp_vec.SetPhi(kaon2_part->Phi());
    tmp_vec.SetM(kaon2_part->M());
    moth_vec+=tmp_vec;
    const int proton_id = part->GetDaughter(1); // Proton
    AliAODMCParticle *proton_part = (AliAODMCParticle*)stack->At(TMath::Abs(proton_id));
    int proton_pdg = TMath::Abs(proton_part->GetPdgCode());
    tmp_vec.SetPt(proton_part->Pt());
    tmp_vec.SetEta(proton_part->Eta());
    tmp_vec.SetPhi(proton_part->Phi());
    tmp_vec.SetM(proton_part->M());
    moth_vec+=tmp_vec;

    // std::cout << "Mother #" << numero++ << ":" << std::endl;
    // std::cout << "  n daughters: " << part->GetNDaughters() << std::endl;
    // std::cout << "  n daughters phi meson: " << phi_dah_n << std::endl;
    // std::cout << "  first daughter PDG: " << phi_pdg << std::endl;
    // std::cout << "  last daughter PDG: " << proton_pdg << std::endl;
    // std::cout << "  Pt: " << moth_vec.Pt() << std::endl;
    // std::cout << "  Eta: " << moth_vec.Eta() << std::endl;
    // std::cout << "  Phi: " << moth_vec.Phi() << std::endl;
    // std::cout << "  M: " << moth_vec.M() << std::endl;
    // std::cout << "___________________________________________" <<std::endl;

    fProduction[iS][iC]->Fill(part->M(),part->Pt());
    if(visible_dacay){
      if ( (part->Y() < fRequireYmin || part->Y() > fRequireYmax) ) continue;
      fTotal[iS][iC]->Fill(moth_vec.M(),moth_vec.Pt());
    }
  }

  /// Checking how many pentaquarks in acceptance are reconstructed well
  std::vector<mother_struct> mothers;
  mothers.reserve(40);
  //int protoni = 0;
  //int kaoni = 0;
  for (Int_t iT = 0; iT < (Int_t)ev->GetNumberOfTracks(); ++iT) {
    AliAODTrack *track = static_cast<AliAODTrack*>(ev->GetTrack(iT));
    if (track->GetID() <= 0) continue;
    if (!track->TestFilterBit(fFilterBit)) continue;
    AliAODMCParticle *part = (AliAODMCParticle*)stack->At(TMath::Abs(track->GetLabel()));
    if (!part) continue;
    const int pdg = TMath::Abs(part->GetPdgCode());
    if (pdg != 321 && pdg != 2212) continue;
    //std::cout << "HasTOF : " << AliAnalysisTaskPsEfficiency::HasTOF(track) << std::endl;
    const int mother_id = part->GetMother();
    AliAODMCParticle* mother = (mother_id >= 0) ? (AliAODMCParticle*)stack->At(mother_id) : nullptr;
    if (!mother) continue;
    const int mother_pdg = TMath::Abs(mother->GetPdgCode());
    //Check wheter the track belongs to a proton
    if (pdg == 2212 && (mother_pdg == 9322134 || mother_pdg == 9322136)) {
      //protoni++;
      if(TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kProton))>3.) continue;
      FourVector_t tmp_vec;
      tmp_vec.SetPt((float)track->Pt());
      tmp_vec.SetEta((float)track->Eta());
      tmp_vec.SetPhi((float)track->Phi());
      tmp_vec.SetM((float)track->M(AliAODTrack::kProton));
      auto it = std::find(mothers.begin(),mothers.end(), mother_id);
      if (it == mothers.end()){
        mother_struct tmp_mum;
        tmp_mum.id=mother_id;
        tmp_mum.phi_tof=true;
        tmp_mum.proton_tof= AliAnalysisTaskPsEfficiency::HasTOF(track) && (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kProton))<3.);
        tmp_mum.n_daughters=1;
        tmp_mum.vec=tmp_vec;
        mothers.push_back(tmp_mum);
      }
      else{
        it->n_daughters++;
        it->proton_tof*=AliAnalysisTaskPsEfficiency::HasTOF(track) && (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kProton))<3.);
        it->vec+=tmp_vec;
      }
    }
    //Check wether the track belgons to a kaon from the decay of a phi
    if (pdg == 321 && mother_pdg == 333) {
      if(TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kKaon))>3.) continue;
      const int ancestor_id = mother->GetMother();
      AliAODMCParticle* ancestor = (ancestor_id >= 0) ? (AliAODMCParticle*)stack->At(ancestor_id) : nullptr;
      if (!ancestor) continue;
      const int ancestor_pdg = TMath::Abs(ancestor->GetPdgCode());
      if (ancestor_pdg == 9322134 || ancestor_pdg == 9322136) {
        //kaoni++;
        FourVector_t tmp_vec;
        tmp_vec.SetPt((float)track->Pt());
        tmp_vec.SetEta((float)track->Eta());
        tmp_vec.SetPhi((float)track->Phi());
        tmp_vec.SetM((float)track->M(AliAODTrack::kKaon));
        auto it = std::find(mothers.begin(),mothers.end(), ancestor_id);
        if (it == mothers.end()){
          mother_struct tmp_mum;
          tmp_mum.id=ancestor_id;
          tmp_mum.phi_tof=AliAnalysisTaskPsEfficiency::HasTOF(track) && (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kKaon))<3.);
          tmp_mum.proton_tof=true;
          tmp_mum.n_daughters=1;
          tmp_mum.vec=tmp_vec;
          mothers.push_back(tmp_mum);
        }
        else{
          it->n_daughters++;
          it->phi_tof*=AliAnalysisTaskPsEfficiency::HasTOF(track) && (TMath::Abs(fPID->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kKaon))<3.);
          it->vec+=tmp_vec;
        }
      }
    }

  } // End AOD track loop

  // int numero = 0;
  // int part_numero = 0;
  // std::cout << "LISTING ALL THE MOTHERS: "<<std::endl<<std::endl;
  for (const auto& mum : mothers) {
    // std::cout << "Mother #" << part_numero++ << ":" << std::endl;
    // std::cout << "  ID:" << mum.id << std::endl;
    // std::cout << "  Daughters: " << mum.n_daughters << std::endl;
    // std::cout << "  Phi tof: " << mum.phi_tof << std::endl;
    // std::cout << "  Proton tof: " << mum.proton_tof << std::endl;
    // std::cout << "  Pt: " << mum.vec.Pt() << std::endl;
    // std::cout << "  Eta: " << mum.vec.Eta() << std::endl;
    // std::cout << "  Phi: " << mum.vec.Phi() << std::endl;
    // std::cout << "  M: " << mum.vec.M() << std::endl;
    // std::cout << "___________________________________________" <<std::endl;
    if (mum.n_daughters != 3) continue;
    if (mum.vec.Rapidity() < fRequireYmin || mum.vec.Rapidity() > fRequireYmax) continue;
    //numero++;
    AliAODMCParticle *part = static_cast<AliAODMCParticle*>(stack->At(mum.id));
    const int iC = part->Charge() > 0 ? 1 : 0;
    const int pdg = TMath::Abs(part->GetPdgCode());
    const int iS = pdg == 9322134 ? 0 : pdg == 9322136 ? 1 : -1;
    const float pt_rec = mum.vec.Pt();
    const float mass_rec = mum.vec.M();
    fReconstructed[iS][iC][0]->Fill(mass_rec,pt_rec);
    if(mum.phi_tof && mum.proton_tof) fReconstructed[iS][iC][1]->Fill(mass_rec,pt_rec);
    if(mum.proton_tof) fReconstructed[iS][iC][2]->Fill(mass_rec,pt_rec);
  }
  // std::cout << "Number of selected protons: " << protoni << std::endl;
  // std::cout << "Number of selected kaons: " << kaoni << std::endl;
  // std::cout << "Number of reconstructed (at least partially): " << mothers.size() << std::endl;
  // std::cout << "Number of reconstructed: " << numero << std::endl;
  //  Post output data.
  PostData(1,fList);
}

/// This function checks whether a track has or has not a prolongation in the TOF.
///
/// \param track Track that has to be checked
/// \return true if the track has a matching hit in the TOF.
///
bool AliAnalysisTaskPsEfficiency::HasTOF(AliVTrack *track) {
  bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  return hasTOFout && hasTOFtime && (track->GetIntegratedLength() > 350.);
}
