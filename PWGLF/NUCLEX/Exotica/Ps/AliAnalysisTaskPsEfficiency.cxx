#include "AliAnalysisTaskPsEfficiency.h"

// ROOT includes
#include <TChain.h>
#include <TClonesArray.h>
#include <TH2F.h>
#include <TList.h>
#include <TMath.h>

// ALIROOT includes
#include "AliPDG.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"

using TMath::TwoPi;

///\cond CLASSIMP
ClassImp(AliAnalysisTaskPsEfficiency);
///\endcond

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskPsEfficiency::AliAnalysisTaskPsEfficiency(TString taskname) : AliAnalysisTaskSE(taskname.Data()),
  fEventCut{false},
  fFilterBit{BIT(8)},
  fPDG{9322134},
  fList{nullptr},
  fProduction{nullptr},
  fReconstructed{{nullptr}},
  fTotal{nullptr}
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

/// Standard destructor
///
AliAnalysisTaskPsEfficiency::~AliAnalysisTaskPsEfficiency() {
  if (fList) delete fList;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskPsEfficiency::UserCreateOutputObjects() {

  fList = new TList();
  fList->SetOwner(true);

  char   letter[2] = {'A','M'};
  string tpctof[2] = {"TPC","TOF"};
  string tpctofMC[2] = {"TPC","TPC_TOF"};

  fProduction = new TH2F("fProduction",";P_{s} state;#it{p} (GeV/#it{c});Entries",2,-0.5,1.5,100,-10,10);
  fList->Add(fProduction);

  for (int iC = 0; iC < 2; ++iC) {
    fTotal[iC] = new TH2F(Form("f%cTotal",letter[iC]),";P_{s} state;#it{p}_{T} (GeV/#it{c}); Counts",2,-0.5,1.5,20,0,10);
    fList->Add(fTotal[iC]);
    for (int iT = 0; iT < 2; ++iT) {
      fReconstructed[iT][iC] = new TH2F(Form("f%cITS_%s",letter[iC],tpctofMC[iT].data()),";P_{s} state;#it{p}_{T} (GeV/#it{c}); Counts",2,-0.5,1.5,20,0,10);
      fList->Add(fReconstructed[iT][iC]);
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

  TClonesArray *stack = nullptr;
  // get branch "mcparticles"
  stack = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (!stack)
    ::Fatal("AliAnalysisTaskPsEfficiency::UserExec","MC analysis requested on a sample without the MC particle array.");

  /// Making the list of the pentaquarks we want to measure
  for (int iMC = 0; iMC < stack->GetEntriesFast(); ++iMC) {
    AliAODMCParticle *part = (AliAODMCParticle*)stack->UncheckedAt(iMC);
    const int pdg = TMath::Abs(part->GetPdgCode());
    const int iC = part->Charge() > 0 ? 1 : 0;
    const int mult = -1 + 2 * iC;
    const int iS = pdg == 9322134 ? 0 : pdg == 9322136 ? 1 : -1;
    if (iS == -1) continue;
    fProduction->Fill(iS,mult * part->P());
    if (TMath::Abs(part->Y()) > 0.5) continue;
    if (part->IsPhysicalPrimary()) fTotal[iC]->Fill(iS,part->Pt());
  }

  /// Checking how many pentaquarks in acceptance are reconstructed well
  std::vector<AliAnalysisTaskPsEfficiency::mother_struct> mothers;
  mothers.reserve(40);
  for (Int_t iT = 0; iT < (Int_t)ev->GetNumberOfTracks(); ++iT) {
    AliAODTrack *track = static_cast<AliAODTrack*>(ev->GetTrack(iT));

    if (track->GetID() <= 0) continue;
    if (!track->TestBit(fFilterBit)) continue;

    AliAODMCParticle *part = (AliAODMCParticle*)stack->At(TMath::Abs(track->GetLabel()));
    if (!part) continue;
    const int pdg = TMath::Abs(part->GetPdgCode());
    if (pdg != 321 || pdg != 2212) continue;

    const int mother_id = part->GetMother();
    AliAODMCParticle* mother = (mother_id >= 0) ? (AliAODMCParticle*)stack->At(mother_id) : nullptr;
    if (!mother) continue;
    const int mother_pdg = TMath::Abs(mother->GetPdgCode());
    //Check wheter the track belongs to a proton
    if (pdg == 2212 && (mother_pdg == 9322134 || mother_pdg == 9322136)) {
      auto it = std::find(mothers.begin(),mothers.end(), mother_id);
      if (it == mothers.end())
      mothers.push_back({mother_id,AliAnalysisTaskPsEfficiency::HasTOF(track),1,(float)track->Px(),(float)track->Py()});
      else{
        it->n_daughters++;
        it->tof*=AliAnalysisTaskPsEfficiency::HasTOF(track);
        it->Px+=track->Px();
        it->Py+=track->Py();
      }
    }
    //Check wether the track belgons to a kaon from the decay of a phi
    if (pdg == 321 && mother_pdg == 333) {
      const int ancestor_id = mother->GetMother();
      AliAODMCParticle* ancestor = (ancestor_id >= 0) ? (AliAODMCParticle*)stack->At(ancestor_id) : nullptr;
      if (!ancestor) continue;
      const int ancestor_pdg = TMath::Abs(ancestor->GetPdgCode());
      if (ancestor_pdg == 9322134 || ancestor_pdg == 9322136) {
        auto it = std::find(mothers.begin(),mothers.end(), ancestor_id);
        if (it == mothers.end()) mothers.push_back({ancestor_id,AliAnalysisTaskPsEfficiency::HasTOF(track),1,(float)track->Px(),(float)track->Py()});
        else{
          it->n_daughters++;
          it->tof*=AliAnalysisTaskPsEfficiency::HasTOF(track);
          it->Px+=track->Px();
          it->Py+=track->Py();
        }
      }
    }

    //const int iTof = int(HasTOF(track));
    //float pT = track->Pt() * fCharge;
    //const int iC = part->Charge() > 0 ? 1 : 0;
    //for (int iR = iTof; iR >= 0; iR--) {
    //  fReconstructed[iR][iC]->Fill(centrality,part->Pt());
    //}

  } // End AOD track loop

  for (const auto& mum : mothers) {
    if (mum.n_daughters != 3) continue;
    AliAODMCParticle *part = static_cast<AliAODMCParticle*>(stack->At(mum.id));
    const int iC = part->Charge() > 0 ? 1 : 0;
    const int pdg = TMath::Abs(part->GetPdgCode());
    const int iS = pdg == 9322134 ? 0 : pdg == 9322136 ? 1 : -1;
    const int mult = -1 + 2 * iC;
    const float pt_rec = mult*TMath::Sqrt(mum.Px*mum.Px+mum.Py*mum.Py);
    fReconstructed[0][iC]->Fill(iS,pt_rec);
    if(mum.tof) fReconstructed[1][iC]->Fill(iS,pt_rec);
  }

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
