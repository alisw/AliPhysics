// $Id$
//
// Emcal jet finder task.
//
// Authors: C.Loizides, S.Aiola

#include <vector>
#include "AliEmcalJetTask.h"

#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliFJWrapper.h"
#include "AliMCEvent.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVParticle.h"

ClassImp(AliEmcalJetTask)

//________________________________________________________________________
AliEmcalJetTask::AliEmcalJetTask() : 
  AliAnalysisTaskSE("AliEmcalJetTask"),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fJetType(kNone),
  fConstSel(kAllJets),
  fMCConstSel(kAllJets),
  fMarkConst(kFALSE),
  fRadius(0.4),
  fMinJetTrackPt(0.15),
  fMinJetClusPt(0.15),
  fPhiMin(0),
  fPhiMax(TMath::TwoPi()),
  fEtaMin(-1),
  fEtaMax(1),
  fMinJetArea(0.01),
  fMinJetPt(1.0),
  fGhostArea(0.01),
  fIsInit(0),
  fIsPSelSet(0),
  fIsMcPart(0),
  fIsEmcPart(0),
  fJets(0),
  fEvent(0),
  fTracks(0),
  fClus(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliEmcalJetTask::AliEmcalJetTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fJetType(kAKT|kFullJet|kRX1Jet),
  fConstSel(kAllJets),
  fMCConstSel(kAllJets),
  fMarkConst(kFALSE),
  fRadius(0.4),
  fMinJetTrackPt(0.15),
  fMinJetClusPt(0.15),
  fPhiMin(0),
  fPhiMax(TMath::TwoPi()),
  fEtaMin(-1),
  fEtaMax(1),
  fMinJetArea(0.01),
  fMinJetPt(1.0),
  fGhostArea(0.01),
  fIsInit(0),
  fIsPSelSet(0),
  fIsMcPart(0),
  fIsEmcPart(0),
  fJets(0),
  fEvent(0),
  fTracks(0),
  fClus(0)
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

  if (!fIsInit) {
    if (!DoInit())
      return;
    fIsInit = kTRUE;
  }

  FindJets();
}

//________________________________________________________________________
void AliEmcalJetTask::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

//________________________________________________________________________
void AliEmcalJetTask::FindJets()
{
  // Find jets.

  if (!fTracks && !fClus)
    return;

  TString name("kt");
  fastjet::JetAlgorithm jalgo(fastjet::kt_algorithm);
  if ((fJetType & kAKT) != 0) {
    name  = "antikt";
    jalgo = fastjet::antikt_algorithm;
  }

  if ((fJetType & kR020Jet) != 0)
    fRadius = 0.2;
  else if ((fJetType & kR030Jet) != 0)
    fRadius = 0.3; 
  else if ((fJetType & kR040Jet) != 0)
    fRadius = 0.4;

  // setup fj wrapper
  AliFJWrapper fjw(name, name);
  fjw.SetAreaType(fastjet::active_area_explicit_ghosts);
  fjw.SetGhostArea(fGhostArea);
  fjw.SetR(fRadius);
  fjw.SetAlgorithm(jalgo);  
  fjw.SetMaxRap(1);
  fjw.Clear();

  // get primary vertex
  Double_t vertex[3] = {0, 0, 0};
  fEvent->GetPrimaryVertex()->GetXYZ(vertex);

  AliDebug(2,Form("Jet type = %d)", fJetType));

  if ((fIsMcPart || ((fJetType & kFullJet) != 0) || ((fJetType & kChargedJet) != 0)) && fTracks) {
    const Int_t Ntracks = fTracks->GetEntries();
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliVParticle *t = static_cast<AliVParticle*>(fTracks->At(iTracks));
      if (!t)
        continue;
      if (fIsMcPart) {
	if (((fJetType & kChargedJet) != 0) && (t->Charge() == 0))
	  continue;
	if (((fJetType & kNeutralJet) != 0) && (t->Charge() != 0))
	  continue;
      }
      if (fIsMcPart || t->GetLabel() != 0) {
	if (fMCConstSel != kAllJets && t->TestBits(fMCConstSel) != (Int_t)fMCConstSel) {
	  AliDebug(2,Form("Skipping track %d because it does not match the bit mask (%d, %d)", iTracks, fMCConstSel, t->TestBits(TObject::kBitMask)));
	  continue;
	}
      }
      else {
	if (fConstSel != kAllJets && t->TestBits(fConstSel) != (Int_t)fConstSel) {
	  AliDebug(2,Form("Skipping track %d because it does not match the bit mask (%d, %d)", iTracks, fConstSel, t->TestBits(TObject::kBitMask)));
	  continue;
	}
      }
      if (t->Pt() < fMinJetTrackPt) 
        continue;
      Double_t eta = t->Eta();
      Double_t phi = t->Phi();
      if ((eta<fEtaMin) || (eta>fEtaMax) ||
          (phi<fPhiMin) || (phi>fPhiMax))
        continue;

      // offset of 100 for consistency with cluster ids
      fjw.AddInputVector(t->Px(), t->Py(), t->Pz(), t->P(), iTracks + 100);  
    }
  }

  if ((((fJetType & kFullJet) != 0) || ((fJetType & kNeutralJet) != 0)) && fClus) {
    const Int_t Nclus = fClus->GetEntries();
    for (Int_t iClus = 0; iClus < Nclus; ++iClus) {
      AliVCluster *c = 0;
      Double_t cEta=0,cPhi=0,cPt=0;
      Double_t cPx=0,cPy=0,cPz=0;
      if (fIsEmcPart) {
	AliEmcalParticle *ep = static_cast<AliEmcalParticle*>(fClus->At(iClus));
	if (!ep)
	  continue;

	c = ep->GetCluster();
	if (!c)
	  continue;

	if (c->GetLabel() > 0) {
	  if (fMCConstSel != kAllJets && ep->TestBits(fMCConstSel) != (Int_t)fMCConstSel)
	    continue;
	}
	else {
	  if (fConstSel != kAllJets && ep->TestBits(fConstSel) != (Int_t)fConstSel)
	    continue;
	}

	cEta = ep->Eta();
	cPhi = ep->Phi();
	cPt  = ep->Pt();
	cPx  = ep->Px();
	cPy  = ep->Py();
	cPz  = ep->Pz();
      } else {
	c = static_cast<AliVCluster*>(fClus->At(iClus));
	if (!c)
	  continue;

	if (c->GetLabel() > 0) {
	  if (fMCConstSel != kAllJets && c->TestBits(fMCConstSel) != (Int_t)fMCConstSel)
	    continue;
	}
	else {
	  if (fConstSel != kAllJets && c->TestBits(fConstSel) != (Int_t)fConstSel)
	    continue;
	}

	TLorentzVector nP;
	c->GetMomentum(nP, vertex);
	cEta = nP.Eta();
	cPhi = nP.Phi();
	cPt  = nP.Pt();
	cPx  = nP.Px();
	cPy  = nP.Py();
	cPz  = nP.Pz();
      }
      if (!c->IsEMCAL())
	continue;
      if (cPt < fMinJetClusPt) 
	continue;
      if ((cEta<fEtaMin) || (cEta>fEtaMax) ||
	  (cPhi<fPhiMin) || (cPhi>fPhiMax))
	continue;
      // offset of 100 to skip ghost particles uid = -1
      fjw.AddInputVector(cPx, cPy, cPz, TMath::Sqrt(cPx*cPx+cPy*cPy+cPz*cPz), -iClus - 100);
    }
  }
  
  // run jet finder
  fjw.Run();

  // get geometry
  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  if (!geom) {
    AliFatal(Form("%s: Can not create geometry", GetName()));
    return;
  }

  // loop over fastjet jets
  std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();
  for (UInt_t ij=0, jetCount=0; ij<jets_incl.size(); ++ij) {
    if (jets_incl[ij].perp()<fMinJetPt) 
      continue;
    if (fjw.GetJetArea(ij)<fMinJetArea)
      continue;

    AliEmcalJet *jet = new ((*fJets)[jetCount]) 
      AliEmcalJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());

    // loop over constituents
    std::vector<fastjet::PseudoJet> constituents(fjw.GetJetConstituents(ij));
    jet->SetNumberOfTracks(constituents.size());
    jet->SetNumberOfClusters(constituents.size());

    Int_t nt            = 0;
    Int_t nc            = 0;
    Double_t neutralE   = 0;
    Double_t maxCh      = 0;
    Double_t maxNe      = 0;
    Int_t gall          = 0;
    Int_t gemc          = 0;
    Int_t cemc          = 0;
    Int_t ncharged      = 0;
    Int_t nneutral      = 0;
    Double_t mcpt       = 0;
    Double_t emcpt      = 0;

    for(UInt_t ic = 0; ic < constituents.size(); ++ic) {
      Int_t uid = constituents[ic].user_index();

      if ((uid == -1) /*&& (constituents[ic].kt2() < 1e-25)*/) { //ghost particle
        ++gall;
        Double_t gphi = constituents[ic].phi();
        if (gphi<0) 
          gphi += TMath::TwoPi();
        gphi *= TMath::RadToDeg();
        Double_t geta = constituents[ic].eta();
        if ((gphi > geom->GetArm1PhiMin()) && (gphi < geom->GetArm1PhiMax()) &&
            (geta > geom->GetArm1EtaMin()) && (geta < geom->GetArm1EtaMax()))
          ++gemc; 
      }	else if ((uid > 0) && fTracks) { // track constituent
	Int_t tid = uid - 100;
        AliVParticle *t = static_cast<AliVParticle*>(fTracks->At(tid));
        if (!t)
          continue;
	if (fMarkConst)
	  t->SetBit(fJetType);
        Double_t cEta = t->Eta();
        Double_t cPhi = t->Phi();
        Double_t cPt  = t->Pt();
        Double_t cP   = t->P();
        if (t->Charge() == 0) {
          neutralE += cP;
          ++nneutral;
          if (cPt > maxNe)
            maxNe = cPt;
        } else {
          ++ncharged;
          if (cPt > maxCh)
            maxCh = cPt;
        }
        if (fIsMcPart || t->GetLabel() != 0) // check if MC particle
          mcpt += cPt;

        if (cPhi<0) 
          cPhi += TMath::TwoPi();
        cPhi *= TMath::RadToDeg();
        if ((cPhi > geom->GetArm1PhiMin()) && (cPhi < geom->GetArm1PhiMax()) &&
            (cEta > geom->GetArm1EtaMin()) && (cEta < geom->GetArm1EtaMax())) {
          emcpt += cPt;
          ++cemc;
        }

        jet->AddTrackAt(tid, nt);
        ++nt;
      } else if (fClus) { // cluster constituent
	Int_t cid = -uid - 100;
	AliVCluster *c = 0;
        Double_t cEta=0,cPhi=0,cPt=0,cP=0;
        if (fIsEmcPart) {
          AliEmcalParticle *ep = static_cast<AliEmcalParticle*>(fClus->At(cid));
          if (!ep)
            continue;
          c = ep->GetCluster();
          if (!c)
            continue;
	  if (fMarkConst)
	    ep->SetBit(fJetType);
          cEta = ep->Eta();
          cPhi = ep->Phi();
          cPt  = ep->Pt();
          cP   = ep->P();
        } else {
          c = static_cast<AliVCluster*>(fClus->At(cid));
          if (!c)
            continue;
	  if (fMarkConst)
	    c->SetBit(fJetType);
          TLorentzVector nP;
          c->GetMomentum(nP, vertex);
          cEta = nP.Eta();
          cPhi = nP.Phi();
          cPt  = nP.Pt();
          cP   = nP.P();
        }

        neutralE += cP;
        if (cPt > maxNe)
          maxNe = cPt;

        if (c->GetLabel() > 0) // MC particle
          mcpt += cPt;

        if (cPhi<0) 
          cPhi += TMath::TwoPi();
        cPhi *= TMath::RadToDeg();
        if ((cPhi > geom->GetArm1PhiMin()) && (cPhi < geom->GetArm1PhiMax()) &&
            (cEta > geom->GetArm1EtaMin()) && (cEta < geom->GetArm1EtaMax())) {
          emcpt += cPt;
          ++cemc;
        }

        jet->AddClusterAt(cid, nc);
        ++nc;
        ++nneutral;
      } else {
        AliError(Form("%s: No logical way to end up here.", GetName()));
        continue;
      }
    } /* end of constituent loop */

    jet->SetNumberOfTracks(nt);
    jet->SetNumberOfClusters(nc);
    jet->SortConstituents();
    jet->SetMaxNeutralPt(maxNe);
    jet->SetMaxChargedPt(maxCh);
    jet->SetNEF(neutralE / jet->E());
    fastjet::PseudoJet area(fjw.GetJetAreaVector(ij));
    jet->SetArea(area.perp());
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetNumberOfCharged(ncharged);
    jet->SetNumberOfNeutrals(nneutral);
    jet->SetMCPt(mcpt);
    jet->SetNEmc(cemc);
    jet->SetPtEmc(emcpt);

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
  fJets->Sort();
}

//________________________________________________________________________
Bool_t AliEmcalJetTask::DoInit()
{
  // Init. Return true if successful.

  // get input collections
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();

  // get the event
  fEvent = InputEvent();
  if (!fEvent) {
    AliError(Form("%s: Could not retrieve event! Returning", GetName()));
    return 0;
  }

  // add jets to event if not yet there
  fJets->Delete();
  if (!(fEvent->FindListObject(fJetsName)))
    fEvent->AddObject(fJets);
  else {
    AliError(Form("%s: Object with name %s already in event! Returning", GetName(), fJetsName.Data()));
    return 0;
  }

  if (fTracksName == "Tracks")
    am->LoadBranch("Tracks");
  if (!fTracks && !fTracksName.IsNull()) {
    fTracks = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fTracksName));
    if (!fTracks) {
      AliError(Form("%s: Pointer to tracks %s == 0", GetName(), fTracksName.Data()));
      return 0;
    }
  }
  if (fTracks) {
    TClass cls(fTracks->GetClass()->GetName());
    if (cls.InheritsFrom("AliMCParticle"))
      fIsMcPart = 1;
  }
  
  if (fCaloName == "CaloClusters")
    am->LoadBranch("CaloClusters");
  if (!fClus && !fCaloName.IsNull()) {
    fClus = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fCaloName));
    if (!fClus) {
      AliError(Form("%s: Pointer to clus %s == 0", GetName(), fCaloName.Data()));
      return 0;
    }
  }
  if (fClus) {
    TClass cls(fClus->GetClass()->GetName());
    if (cls.InheritsFrom("AliEmcalParticle"))
      fIsEmcPart = 1;
  }

  return 1;
}
