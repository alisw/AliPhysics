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
#include <TMath.h>

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
  fJetPhiMin(-10),
  fJetPhiMax(10),
  fJetEtaMin(-1),
  fJetEtaMax(1),
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
  fJetPhiMin(-10),
  fJetPhiMax(10),
  fJetEtaMin(-1),
  fJetEtaMax(1),
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
    AliDebug(1,"Using AKT algorithm");
  }
  else {
    AliDebug(1,"Using KT algorithm");
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

  AliDebug(2,Form("Jet type = %d", fJetType));

  if ((fIsMcPart || ((fJetType & kFullJet) != 0) || ((fJetType & kChargedJet) != 0)) && fTracks) {
    const Int_t Ntracks = fTracks->GetEntries();
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliVParticle *t = static_cast<AliVParticle*>(fTracks->At(iTracks));
      if (!t)
        continue;
      if (fIsMcPart) {
	if (((fJetType & kChargedJet) != 0) && (t->Charge() == 0)) {
	  AliDebug(2,Form("Skipping track %d because it is neutral.", iTracks));
	  continue;
	}
	if (((fJetType & kNeutralJet) != 0) && (t->Charge() != 0)) {
	  AliDebug(2,Form("Skipping track %d because it is charged.", iTracks));
	  continue;
	}
      }
      if (fIsMcPart || t->GetLabel() != 0) {
	if (fMCConstSel == kNone) {
	  AliDebug(2,Form("Skipping track %d because bit mask is 0.", iTracks));
	  continue;
	}
	if (fMCConstSel != kAllJets) {
	  if (t->TestBits(fMCConstSel) != (Int_t)fMCConstSel) {
	    AliDebug(2,Form("Skipping track %d because it does not match the bit mask (%d, %d)", iTracks, fMCConstSel, t->TestBits(TObject::kBitMask)));
	    continue;
	  }
	  else {
	    AliDebug(2,Form("Track %d matches the bit mask (%d, %d)", iTracks, fMCConstSel, t->TestBits(TObject::kBitMask)));
	  }
	}
      }
      else {
	if (fConstSel == kNone) {
	  AliDebug(2,Form("Skipping track %d because bit mask is 0.", iTracks));
	  continue;
	}
	if (fConstSel != kAllJets) {
	  if (t->TestBits(fConstSel) != (Int_t)fConstSel) {
	    AliDebug(2,Form("Skipping track %d because it does not match the bit mask (%d, %d)", iTracks, fConstSel, t->TestBits(TObject::kBitMask)));
	    continue;
	  }
	  else {
	    AliDebug(2,Form("Track %d matches the bit mask (%d, %d)", iTracks, fConstSel, t->TestBits(TObject::kBitMask)));	    
	  }
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
      AliDebug(2,Form("Track %d accepted (label = %d, pt = %f)", iTracks, t->GetLabel(), t->Pt()));
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
	  if (fMCConstSel == kNone) {
	    AliDebug(2,Form("Skipping cluster %d because bit mask is 0.", iClus));
	    continue;
	  }
	  if (fMCConstSel != kAllJets) {
	    if (ep->TestBits(fMCConstSel) != (Int_t)fMCConstSel) {
	      AliDebug(2,Form("Skipping cluster %d because it does not match the bit mask (%d, %d)", iClus, fMCConstSel, ep->TestBits(TObject::kBitMask)));
	      continue;
	    }
	    else {
	      AliDebug(2,Form("Cluster %d matches the bit mask (%d, %d)", iClus, fMCConstSel, ep->TestBits(TObject::kBitMask)));
	    }
	  }
	}
	else {
	  if (fConstSel == kNone) {
	    AliDebug(2,Form("Skipping cluster %d because bit mask is 0.", iClus));
	    continue;
	  }
	  if (fConstSel != kAllJets) {
	    if (ep->TestBits(fConstSel) != (Int_t)fConstSel) {
	      AliDebug(2,Form("Skipping cluster %d because it does not match the bit mask (%d, %d)", iClus, fConstSel, ep->TestBits(TObject::kBitMask)));
	      continue;
	    }
	    else {
	      AliDebug(2,Form("Cluster %d matches the bit mask (%d, %d)", iClus, fConstSel, ep->TestBits(TObject::kBitMask)));	    
	    }
	  }
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
	  if (fMCConstSel == kNone) {
	    AliDebug(2,Form("Skipping cluster %d because bit mask is 0.", iClus));
	    continue;
	  }
	  if (fMCConstSel != kAllJets) {
	    if (c->TestBits(fMCConstSel) != (Int_t)fMCConstSel) {
	      AliDebug(2,Form("Skipping cluster %d because it does not match the bit mask (%d, %d)", iClus, fMCConstSel, c->TestBits(TObject::kBitMask)));
	      continue;
	    }
	    else {
	      AliDebug(2,Form("Cluster %d matches the bit mask (%d, %d)", iClus, fMCConstSel, c->TestBits(TObject::kBitMask)));
	    }
	  }
	}
	else {
	  if (fConstSel == kNone) {
	    AliDebug(2,Form("Skipping cluster %d because bit mask is 0.", iClus));
	    continue;
	  }
	  if (fConstSel != kAllJets) {
	    if (c->TestBits(fConstSel) != (Int_t)fConstSel) {
	      AliDebug(2,Form("Skipping cluster %d because it does not match the bit mask (%d, %d)", iClus, fConstSel, c->TestBits(TObject::kBitMask)));
	      continue;
	    }
	    else {
	      AliDebug(2,Form("Cluster %d matches the bit mask (%d, %d)", iClus, fConstSel, c->TestBits(TObject::kBitMask)));	    
	    }
	  }
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
      AliDebug(2,Form("Cluster %d accepted (label = %d)", iClus, c->GetLabel()));
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
  // sort jets according to jet pt
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, jets_incl);

  AliDebug(1,Form("%d jets found", (Int_t)jets_incl.size()));
  for (UInt_t ijet=0, jetCount=0; ijet<jets_incl.size(); ++ijet) {
    Int_t ij = indexes[ijet];
    AliDebug(3,Form("Jet pt = %f, area = %f", jets_incl[ij].perp(), fjw.GetJetArea(ij)));

    if (jets_incl[ij].perp()<fMinJetPt) 
      continue;
    if (fjw.GetJetArea(ij)<fMinJetArea)
      continue;
    if ((jets_incl[ij].eta()<fJetEtaMin) || (jets_incl[ij].eta()>fJetEtaMax) ||
	(jets_incl[ij].phi()<fJetPhiMin) || (jets_incl[ij].phi()>fJetPhiMax))
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
      AliDebug(3,Form("Processing constituent %d", uid));
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
        if (!t) {
	  AliError(Form("Could not find track %d",tid));
          continue;
	}
	if (jetCount < fMarkConst) {
	  AliDebug(2,Form("Marking track %d with bit map %d", tid, fJetType));
	  t->SetBit(fJetType);
	}
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
	  if (jetCount < fMarkConst)
	    ep->SetBit(fJetType);
          cEta = ep->Eta();
          cPhi = ep->Phi();
          cPt  = ep->Pt();
          cP   = ep->P();
        } else {
          c = static_cast<AliVCluster*>(fClus->At(cid));
          if (!c)
            continue;
	  if (jetCount < fMarkConst)
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

    AliDebug(2,Form("Added jet n. %d, pt = %f, area = %f, constituents = %d", jetCount, jet->Pt(), jet->Area(), (Int_t)constituents.size()));
    jetCount++;
  }
  //fJets->Sort();
}

//________________________________________________________________________
Bool_t AliEmcalJetTask::GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const
{
  // Get the leading jets.

  static Float_t pt[9999] = {0};

  const Int_t n = (Int_t)array.size();

  if (n < 1)
    return kFALSE;
  
  for (Int_t i = 0; i < n; i++) 
    pt[i] = array[i].perp();

  TMath::Sort(n, pt, indexes);

  return kTRUE;
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
    if (cls.InheritsFrom("AliMCParticle") || cls.InheritsFrom("AliAODMCParticle"))
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
