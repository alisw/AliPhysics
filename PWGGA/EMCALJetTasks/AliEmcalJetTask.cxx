// $Id$
//
// Emcal jet finder task.
//
//

#include "AliEmcalJetTask.h"
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliFJWrapper.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"

ClassImp(AliEmcalJetTask)

//________________________________________________________________________
AliEmcalJetTask::AliEmcalJetTask() : 
  AliAnalysisTaskSE("AliEmcalJetTask"),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fAlgo(1),
  fRadius(0.4),
  fType(0),
  fMinJetTrackPt(0.15),
  fMinJetClusPt(0.15),
  fJets(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliEmcalJetTask::AliEmcalJetTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fAlgo(1),
  fRadius(0.4),
  fType(0),
  fMinJetTrackPt(0.15),
  fMinJetClusPt(0.15),
  fJets(0)
{
  // Standard constructor.

  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";
}

//________________________________________________________________________
AliEmcalJetTask::~AliEmcalJetTask()
{
  // Destructor
}

//________________________________________________________________________
void AliEmcalJetTask::UserCreateOutputObjects()
{
  // Create user objects.

  fJets = new TClonesArray("AliEmcalJet");
  fJets->SetName(fJetsName);
}

//________________________________________________________________________
void AliEmcalJetTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  // add jets to event if not yet there
  if (!(InputEvent()->FindListObject(fJetsName)))
    InputEvent()->AddObject(fJets);

  // delete jet output
  fJets->Delete();


  // get input collections
  TClonesArray *tracks = 0;
  TClonesArray *clus   = 0;
  TList *l = InputEvent()->GetList();
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();

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

  // get centrality
  Float_t cent = -1; 
  AliCentrality *centrality = InputEvent()->GetCentrality() ;
  if (centrality)
    cent = centrality->GetCentralityPercentile("V0M");
  else
    cent=99; // probably pp data
  if (cent<0) {
    AliError(Form("Centrality negative: %f", cent));
    return;
  }
      
  FindJets(tracks, clus, fAlgo, fRadius, cent);
}

//________________________________________________________________________
void AliEmcalJetTask::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

//________________________________________________________________________
void AliEmcalJetTask::FindJets(TObjArray *tracks, TObjArray *clus, Int_t algo, Double_t radius, Float_t /*cent*/)
{
  // Find jets.

  if (!tracks && !clus)
    return;

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
      if (t->Pt()<fMinJetTrackPt) 
        continue;
      fjw.AddInputVector(t->Px(), t->Py(), t->Pz(), t->P(), iTracks);
    }
  }

  if (clus) {
    Double_t vertex[3] = {0, 0, 0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    const Int_t Nclus = clus->GetEntries();
    for (Int_t iClus = 0; iClus < Nclus; ++iClus) {
      AliVCluster *c = dynamic_cast<AliVCluster*>(clus->At(iClus));
      if (!c)
        continue;
      if (!c->IsEMCAL())
        continue;
      TLorentzVector nPart;
      c->GetMomentum(nPart, vertex);
      Double_t et = nPart.Pt();
      if (et<fMinJetClusPt) 
        continue;
      fjw.AddInputVector(nPart.Px(), nPart.Py(), nPart.Pz(), nPart.P(), -iClus-1);
    }
  }

  // run jet finder
  fjw.Run();
  
  std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();
  for(UInt_t ij=0, jetCount=0; ij<jets_incl.size(); ++ij) {
    if (jets_incl[ij].perp()<1) 
      continue;
    AliEmcalJet *jet = new ((*fJets)[jetCount]) 
      AliEmcalJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());
    jet->SetArea(fjw.GetJetArea(ij));
    Double_t vertex[3] = {0, 0, 0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    vector<fastjet::PseudoJet> constituents = fjw.GetJetConstituents(ij);
    Double_t neutralE = 0;Double_t maxTrack = 0;Double_t maxCluster=0;
    jet->SetNumberOfTracks(constituents.size());
    jet->SetNumberOfClusters(constituents.size());
    Int_t nt = 0;
    Int_t nc = 0;
    for(UInt_t ic=0; ic<constituents.size(); ++ic) {
      Int_t uid = constituents[ic].user_index();
      if (uid>=0){
	jet->AddTrackAt(uid, nt);
        AliVTrack *t = static_cast<AliVTrack*>(tracks->At(uid));
        if (t) {
          if (t->Pt()>maxTrack)
            maxTrack=t->Pt();
          nt++;
        }
      } else {
	jet->AddClusterAt(-(uid+1),nc);
        AliVCluster *c = static_cast<AliVCluster*>(clus->At(-(uid+1)));
        if (c) {
          TLorentzVector nP;
          c->GetMomentum(nP, vertex);
          neutralE += nP.P();
          if (nP.Pt()>maxCluster)
            maxCluster=nP.Pt();
        }
	nc++;
      }
    }
    jet->SetNumberOfTracks(nt);
    jet->SetNumberOfClusters(nc);
    jet->SetMaxTrackPt(maxTrack);
    jet->SetMaxClusterPt(maxCluster);
    jet->SetNEF(neutralE/jet->E());
    jetCount++;
  }
}
