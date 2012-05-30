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
#include "AliESDEvent.h"
#include "AliMCEvent.h"

#include "AliEmcalJetTask.h"

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

//_____________________________________________________
TString AliEmcalJetTask::GetBeamType()
{
  // Get beam type : pp-AA-pA
  // ESDs have it directly, AODs get it from hardcoded run number ranges

  TString beamType;

  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(fEvent);
  if (esd) {
    const AliESDRun *run = esd->GetESDRun();
    beamType = run->GetBeamType();
  }
  else
  {
    Int_t runNumber = fEvent->GetRunNumber();
    if ((runNumber >= 136851 && runNumber <= 139517) ||  // LHC10h
	(runNumber >= 166529 && runNumber <= 170593))    // LHC11h
    {
      beamType = "A-A";
    }
    else 
    {
      beamType = "p-p";
    }
  }

  return beamType;    
}


//________________________________________________________________________
void AliEmcalJetTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  // get the event
  fEvent = InputEvent();

  if (!fEvent) {
    AliError("Could not retrieve event! Returning");
    return;
  }

  // add jets to event if not yet there
  if (!(fEvent->FindListObject(fJetsName)))
    fEvent->AddObject(fJets);

  // delete jet output
  fJets->Delete();

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
  Double_t cent = 99; 
  
  if (GetBeamType() == "A-A") {
    AliCentrality *centrality = InputEvent()->GetCentrality();
    
    if (centrality)
      cent = centrality->GetCentralityPercentile("V0M");
    else
      cent = 99; // probably pp data
    
    if (cent < 0) {
      AliWarning(Form("Centrality negative: %f, assuming 99", cent));
      cent = 99;
    }
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
      AliVParticle *t = dynamic_cast<AliVParticle*>(tracks->At(iTracks));
      if (!t)
        continue;

      if (t->Pt() < fMinJetTrackPt) 
        continue;
      fjw.AddInputVector(t->Px(), t->Py(), t->Pz(), t->P(), iTracks + 100);  // offset of 100 for consistency with cluster ids
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
      if (nPart.Pt() < fMinJetClusPt) 
        continue;
      fjw.AddInputVector(nPart.Px(), nPart.Py(), nPart.Pz(), nPart.P(), -iClus - 100);  // offset of 100 to skip ghost particles uid = -1
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
    Double_t maxCh      = 0;
    Double_t maxNe      = 0;
    Int_t gall          = 0;
    Int_t gemc          = 0;
    Int_t ncharged      = 0;
    Int_t nneutral      = 0;
    Double_t MCpt       = 0;

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
      }	else if (uid > 0) {
	Int_t tid = uid - 100;
        AliVParticle *t = static_cast<AliVParticle*>(tracks->At(tid));
        if (t) {
	  if (t->Charge() == 0) {
	    neutralE += t->P();
	    ++nneutral;
	    if (t->Pt() > maxNe)
	      maxNe = t->Pt();
	  } else {
	    ++ncharged;
	    if (t->Pt() > maxCh)
	      maxCh = t->Pt();
	  }

	  if (t->InheritsFrom("AliMCParticle") || t->GetLabel() == 100) // MC particle
	    MCpt += t->Pt();

          jet->AddTrackAt(tid, nt);
          ++nt;
        }
      } else {
	Int_t cid = -uid - 100;
	AliVCluster *c = static_cast<AliVCluster*>(clus->At(cid));
	if (c) {
	  TLorentzVector nP;
	  c->GetMomentum(nP, vertex);
	  neutralE += nP.P();
	  if (nP.Pt() > maxNe)
	    maxNe = nP.Pt();

	  if (c->Chi2() == 100) // MC particle
	    MCpt += nP.Pt();

	  jet->AddClusterAt(cid, nc);
	  ++nc;
	  ++nneutral;
	}
      }
    }
    jet->SetNumberOfTracks(nt);
    jet->SetNumberOfClusters(nc);
    jet->SortConstituents();
    jet->SetMaxNeutralPt(maxNe);
    jet->SetMaxChargedPt(maxCh);
    jet->SetNEF(neutralE / jet->E());
    jet->SetArea(fjw.GetJetArea(ij));
    jet->SetNumberOfCharged(ncharged);
    jet->SetNumberOfNeutrals(nneutral);
    jet->SetMCPt(MCpt);
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
