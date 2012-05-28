// $Id$
//
// Emcal jet finder task.
//
// Authors: C.Loizides, S.Aiola

#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalJet.h"
#include "AliFJWrapper.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"

#include "AliEmcalJetTask.h"

ClassImp(AliEmcalJetTask)

//________________________________________________________________________
AliEmcalJetTask::AliEmcalJetTask() : 
  AliAnalysisTaskSE("AliEmcalJetTask"),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fMC(kFALSE),
  fAlgo(1),
  fRadius(0.4),
  fType(0),
  fMinJetTrackPt(0.15),
  fMinJetClusPt(0.15),
  fMinJetArea(0.01),
  fMinJetPt(1.0),
  fJets(0),
  fEvent(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliEmcalJetTask::AliEmcalJetTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fMC(kFALSE),
  fAlgo(1),
  fRadius(0.4),
  fType(0),
  fMinJetTrackPt(0.15),
  fMinJetClusPt(0.15),
  fMinJetArea(0.01),
  fMinJetPt(1.0),
  fJets(0),
  fEvent(0)
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

  // get the event
  if (fMC) 
    fEvent = MCEvent();
  else
    fEvent = InputEvent();

  // add jets to event if not yet there
  if (!(fEvent->FindListObject(fJetsName)))
    fEvent->AddObject(fJets);

  // delete jet output
  fJets->Delete();

  if (!fEvent) {
    AliError(Form("Could not get the event! fMC = %d", fMC));
    return;
  }

  // get input collections
  TClonesArray *tracks = 0;
  TClonesArray *clus   = 0;
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();

  if ((fType == 0) || (fType == 1)) {
    if (fTracksName == "Tracks")
      am->LoadBranch("Tracks");
    tracks = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fTracksName));
    if (!tracks) {
      AliError(Form("Pointer to tracks %s == 0", fTracksName.Data()));
      return;
    }
  }
  if ((fType == 0) || (fType == 2)) {
    if (fCaloName == "CaloClusters")
      am->LoadBranch("CaloClusters");
    clus = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fCaloName));
    if (!clus) {
      AliError(Form("Pointer to clus %s == 0", fCaloName.Data()));
      return;
    }
  }

  // get centrality
  Float_t cent = -1; 
  AliCentrality *centrality = fEvent->GetCentrality();
  if (centrality)
    cent = centrality->GetCentralityPercentile("V0M");
  else
    cent=99; // probably pp data
  if (cent < 0) {
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
  fjw.SetAreaType(fastjet::active_area_explicit_ghosts);
  fjw.SetR(radius);
  fjw.SetAlgorithm(jalgo);  
  fjw.SetMaxRap(0.9);
  fjw.Clear();

  if (tracks) {
    const Int_t Ntracks = tracks->GetEntries();
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliVParticle *t = static_cast<AliVParticle*>(tracks->At(iTracks));
      if (!t)
        continue;
      if (t->Pt()<fMinJetTrackPt) 
        continue;

      Int_t index = 1;
      if(fMC && t->Charge() == 0)
	index =- 1;

      fjw.AddInputVector(t->Px(), t->Py(), t->Pz(), t->P(), (iTracks + 123) * index);
    }
  }

  if (clus) {
    Double_t vertex[3] = {0, 0, 0};
    fEvent->GetPrimaryVertex()->GetXYZ(vertex);
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
      if (et < fMinJetClusPt) 
        continue;
      fjw.AddInputVector(nPart.Px(), nPart.Py(), nPart.Pz(), nPart.P(), -iClus - 123);
    }
  }

  // run jet finder
  fjw.Run();

  // get geometry
  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  if (!geom) {
    AliFatal("Can not create geometry");
    return;
  }
  
  std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();
  for (UInt_t ij=0, jetCount=0; ij<jets_incl.size(); ++ij) {
    if (jets_incl[ij].perp()<fMinJetPt) 
      continue;
    if (fjw.GetJetArea(ij)<fMinJetArea)
      continue;
    AliEmcalJet *jet = new ((*fJets)[jetCount]) 
      AliEmcalJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());
    Double_t vertex[3] = {0, 0, 0};
    fEvent->GetPrimaryVertex()->GetXYZ(vertex);
    vector<fastjet::PseudoJet> constituents = fjw.GetJetConstituents(ij);
    Int_t nt            = 0;
    Int_t nc            = 0;
    Double_t neutralE   = 0;
    Double_t maxTrack   = 0;
    Double_t maxCluster = 0;
    Int_t gall          = 0;
    Int_t gemc          = 0;

    jet->SetNumberOfTracks(constituents.size());
    jet->SetNumberOfClusters(constituents.size());

    for(UInt_t ic = 0; ic < constituents.size(); ++ic) {
      Int_t uid = constituents[ic].user_index();

      if ((uid == -1) && (constituents[ic].kt2() < 1e-25)) { //ghost particle
        ++gall;
        Double_t gphi = constituents[ic].phi() * TMath::RadToDeg();
        Double_t geta = constituents[ic].eta();
        if ((gphi > geom->GetArm1PhiMin()) && (gphi < geom->GetArm1PhiMax()) &&
            (geta > geom->GetArm1EtaMin()) && (geta < geom->GetArm1EtaMax()))
          ++gemc;
        continue;
      }	else if (fMC || uid >= 0) {
	Int_t tid = TMath::Abs(uid) - 123;
        AliVParticle *t = static_cast<AliVParticle*>(tracks->At(tid));
        if (t) {
	  if (uid >= 0)
	    neutralE += t->P();
          if (t->Pt() > maxTrack)
            maxTrack = t->Pt();
          jet->AddTrackAt(tid, nt);
          nt++;
        }
      } else {
	Int_t cid = TMath::Abs(uid) - 123;
	AliVCluster *c = static_cast<AliVCluster*>(clus->At(cid));
	if (c) {
	  TLorentzVector nP;
	  c->GetMomentum(nP, vertex);
	  neutralE += nP.P();
	  if (nP.Pt() > maxCluster)
	    maxCluster = nP.Pt();
	  jet->AddClusterAt(cid, nc);
	  nc++;
	}
      }
    }
    jet->SetNumberOfTracks(nt);
    jet->SetNumberOfClusters(nc);
    jet->SortConstituents();
    jet->SetMaxTrackPt(maxTrack);
    jet->SetMaxClusterPt(maxCluster);
    jet->SetNEF(neutralE / jet->E());
    jet->SetArea(fjw.GetJetArea(ij));
    if (gall > 0)
      jet->SetAreaEmc(fjw.GetJetArea(ij) * gemc / gall);
    else 
      jet->SetAreaEmc(-1);
    if ((jet->Phi() > geom->GetArm1PhiMin() * TMath::DegToRad()) && 
        (jet->Phi() < geom->GetArm1PhiMax() * TMath::DegToRad()) &&
        (jet->Eta() > geom->GetArm1EtaMin()) && 
        (jet->Eta() < geom->GetArm1EtaMax()))
      jet->SetAxisInEmcal(kTRUE);
    jetCount++;
  }
}
