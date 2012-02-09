// $Id$

#include "AliEsdJetTask.h"
#include <TClonesArray.h>
#include <TParticle.h>
#include "AliESDJet.h"
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliFJWrapper.h"
#include "AliESDCaloCluster.h"

ClassImp(AliEsdJetTask)

#ifdef __HAVE_FJINTERFACE__
//________________________________________________________________________
AliEsdJetTask::AliEsdJetTask(const char *name) : 
  AliAnalysisTaskSE("AliEsdJetTask"),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fAlgo(1),
  fRadius(0.4),
  fType(0),
  fHadCorr(0),
  fJets(0)
{
  // Standard constructor.

  if (!name)
    return;

  SetName(name);
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";
}

//________________________________________________________________________
AliEsdJetTask::~AliEsdJetTask()
{
  // Destructor
}

//________________________________________________________________________
void AliEsdJetTask::UserCreateOutputObjects()
{
  // Create user objects.

  fJets = new TClonesArray("AliESDJet");
  fJets->SetName(fJetsName);
}

//________________________________________________________________________
void AliEsdJetTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.
  // Add jets to event if not yet there

  if (!(InputEvent()->FindListObject(fJetsName)))
    InputEvent()->AddObject(fJets);

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  TClonesArray *tracks = 0;
  TClonesArray *clus   = 0;
  TList *l = InputEvent()->GetList();
  if ((fType==0)||(fType==1)) {
    if (fTracksName == "Tracks")
      am->LoadBranch("Tracks");
    tracks = dynamic_cast<TClonesArray*>(l->FindObject(fTracksName));
    if (!tracks) {
      AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
      return;
    }
  }
  if ((fType==0)||(fType==2)) {
    if (fCaloName == "CaloClusters")
      am->LoadBranch("CaloClusters");
    clus = dynamic_cast<TClonesArray*>(l->FindObject(fCaloName));
    if (!clus) {
      AliError(Form("Pointer to clus %s == 0", fCaloName.Data() ));
      return;
    }
  }
      
  FindJets(tracks, clus, fAlgo, fRadius);
}

//________________________________________________________________________
void AliEsdJetTask::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.

}

//________________________________________________________________________
void AliEsdJetTask::FindJets(TObjArray *tracks, TObjArray *clus, Int_t algo, Double_t radius)
{
  // Find jets.

  TString name("kt");
  fastjet::JetAlgorithm jalgo(fastjet::kt_algorithm);
  if (algo>=1) {
    name  = "antikt";
    jalgo = fastjet::antikt_algorithm;
  }

  AliFJWrapper fjw(name, name);
  fjw.SetR(radius);
  fjw.SetAlgorithm(jalgo);
  fjw.SetMaxRap(0.9);
  fjw.Clear();

  if (tracks) {
    const Int_t Ntracks = tracks->GetEntries();
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliVTrack *t = static_cast<AliVTrack*>(tracks->At(iTracks));
      if (!t)
        continue;
      fjw.AddInputVector(t->Px(), t->Py(), t->Pz(), t->P(), iTracks);
    }
  }
  if (clus) {
    Double_t vertex[3] = {0, 0, 0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    const Int_t Nclus = clus->GetEntries();
    for (Int_t iClus = 0, iN = 0; iClus < Nclus; ++iClus) {
      AliVCluster *c = dynamic_cast<AliVCluster*>(clus->At(iClus));
      if (!c->IsEMCAL())
        continue;
      TLorentzVector nPart;
      c->GetMomentum(nPart, vertex);
      Double_t energy = nPart.P();
      if (fHadCorr>0) {
        Int_t imin = static_cast<Int_t>(c->GetEmcCpvDistance());
        if (imin>=0) {
          Double_t dPhiMin = c->GetTrackDx();
          Double_t dEtaMin = c->GetTrackDz();
          if (dPhiMin<0.05 && dEtaMin<0.025) { // pp cuts!!!
            AliVTrack *t = dynamic_cast<AliESDtrack*>(tracks->At(imin));
            if (t) {
              energy -= fHadCorr*t->P();
              if (energy<0)
                continue;
            }
          }
        }
      }
      fjw.AddInputVector(nPart.Px(), nPart.Py(), nPart.Pz(), energy, -iN-1);
    }
  }

  // run jet finder
  fjw.Run();

  std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();
  for(UInt_t ij=0, jetCount=0; ij<jets_incl.size(); ++ij) {
    if (jets_incl[ij].perp()<1) 
      continue;
#if 0
    AliAODJet *jet = new ((*fJets)[jetCount]) 
      AliAODJet(jets_incl[ij].px(), jets_incl[ij].py(), jets_incl[ij].pz(), jets_incl[ij].E());
    jet->SetEffArea(fjw.GetJetArea(ij),0);
    jet->GetRefTracks()->Clear();
    vector<fastjet::PseudoJet> constituents = fjw.GetJetConstituents(ij);
    Double_t neutralE = 0;
    for(UInt_t ic=0; ic<constituents.size(); ++ic) {
      Int_t uid = constituents[ic].user_index();
      if (uid>=0)
        jet->AddTrack(tracks->At(uid));
      else {
        TLorentzVector *nP = static_cast<TLorentzVector*>(fNeutrals->At(-(uid+1)));
        neutralE += nP->E();
        jet->AddTrack(nP);
      }
    }
    jet->SetNEF(neutralE/jet->E());
#endif
    //jet->Print("");
    jetCount++;
  }
}
#else
//________________________________________________________________________
AliEsdJetTask::AliEsdJetTask(const char *name) : 
  AliAnalysisTaskSE("AliEsdJetTask"),
  fPrimTrCuts(0),
  fPrimTracksName("Tracks"),
  fJetsName("Jets"),
  fCaloName("CaloClusters"),
  fAlgo(1),
  fRadius(0.4),
  fType(0),
  fHadCorr(0),
  fJets(0),
  fNeutrals(0)
{
  // Standard constructor.
  AliFatal("Compiled without FASTJET package. Task cannot be used!");
}

//________________________________________________________________________
AliEsdJetTask::~AliEsdJetTask()
{
  // Destructor
}

//________________________________________________________________________
void AliEsdJetTask::UserCreateOutputObjects()
{
  // Create user objects.
}

//________________________________________________________________________
void AliEsdJetTask::UserExec(Option_t *) 
{
}

//________________________________________________________________________
void AliEsdJetTask::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

//________________________________________________________________________
void AliEsdJetTask::FindJets(TObjArray *, TObjArray *, Int_t, Double_t)
{
  // Find jets.
}
#endif
