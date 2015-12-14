//
// Emcal jet finder task.
//
// Authors: C.Loizides, S.Aiola, M. Verweij

#include <vector>
#include "AliEmcalJetTask.h"

#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TClass.h>

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
#include "AliEmcalJetUtility.h"
#include "AliAODTrack.h"

using std::cout;
using std::endl;
using std::cerr;

ClassImp(AliEmcalJetTask)

//________________________________________________________________________
AliEmcalJetTask::AliEmcalJetTask() :
  AliAnalysisTaskSE("AliEmcalJetTask"),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fJetType(kAKT|kFullJet|kRX1Jet),
  fMinLabelTracks(-kMaxInt),
  fMaxLabelTracks(kMaxInt),
  fMinLabelClusters(-kMaxInt),
  fMaxLabelClusters(kMaxInt),
  fMinMCLabel(0),
  fRadius(0.4),
  fMinJetTrackPt(0.15),
  fMinJetClusPt(0.),
  fMinJetClusE(0.30),
  fPhiMin(0),
  fPhiMax(TMath::TwoPi()),
  fEtaMin(-0.9),
  fEtaMax(+0.9),
  fMinJetArea(0.001),
  fMinJetPt(1.0),
  fJetPhiMin(-10),
  fJetPhiMax(+10),
  fJetEtaMin(-1),
  fJetEtaMax(+1),
  fGhostArea(0.005),
  fRecombScheme(fastjet::pt_scheme),
  fTrackEfficiency(1.),
  fMCFlag(0),
  fGeneratorIndex(-1),
  fUtilities(0),
  fFilterHybridTracks(kFALSE),
  fUseExchangeCont(0),
  fClusterEnergyType(-1),
  fLocked(0),
  fIsInit(0),
  fIsPSelSet(0),
  fIsEmcPart(0),
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fGeom(0),
  fJets(0),
  fEvent(0),
  fTracks(0),
  fClus(0),
  fFastJetWrapper("AliEmcalJetTask","AliEmcalJetTask")
{
  // Default constructor.

  fVertex[0] = 0.;
  fVertex[1] = 0.;
  fVertex[2] = 0.;
}

//________________________________________________________________________
AliEmcalJetTask::AliEmcalJetTask(const char *name, Int_t useExchangeCont) :
  AliAnalysisTaskSE(name),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fJetType(kAKT|kFullJet|kRX1Jet),
  fMinLabelTracks(-kMaxInt),
  fMaxLabelTracks(kMaxInt),
  fMinLabelClusters(-kMaxInt),
  fMaxLabelClusters(kMaxInt),
  fMinMCLabel(0),
  fRadius(0.4),
  fMinJetTrackPt(0.15),
  fMinJetClusPt(0.),
  fMinJetClusE(0.30),
  fPhiMin(0),
  fPhiMax(TMath::TwoPi()),
  fEtaMin(-0.9),
  fEtaMax(+0.9),
  fMinJetArea(0.001),
  fMinJetPt(1.0),
  fJetPhiMin(-10),
  fJetPhiMax(+10),
  fJetEtaMin(-1),
  fJetEtaMax(+1),
  fGhostArea(0.005),
  fRecombScheme(fastjet::pt_scheme),
  fTrackEfficiency(1.),
  fMCFlag(0),
  fGeneratorIndex(-1),
  fUtilities(0),
  fFilterHybridTracks(kFALSE),
  fUseExchangeCont(useExchangeCont),
  fClusterEnergyType(-1),
  fLocked(0),
  fIsInit(0),
  fIsPSelSet(0),
  fIsEmcPart(0),
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fGeom(0),
  fJets(0),
  fEvent(0),
  fTracks(0),
  fClus(0),
  fFastJetWrapper(name,name)
{
  // Standard constructor.

  for (Int_t i = 0; i < fUseExchangeCont; i++) {
    DefineInput(i+1, TClonesArray::Class());
  }

  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";

  fVertex[0] = 0.;
  fVertex[1] = 0.;
  fVertex[2] = 0.;
}

//________________________________________________________________________
AliEmcalJetTask::~AliEmcalJetTask()
{
  // Destructor
}

//________________________________________________________________________
AliEmcalJetUtility* AliEmcalJetTask::AddUtility(AliEmcalJetUtility* utility)
{
  // Addition of utilities.

  if (!fUtilities) fUtilities = new TObjArray();
  if (fUtilities->FindObject(utility)) {
    Error("AddSupply", "Jet utility %s already connected.", utility->GetName());
    return utility;
  }   
  fUtilities->Add(utility);
  utility->SetJetTask(this);

  return utility;
}

//________________________________________________________________________
void AliEmcalJetTask::InitUtilities()
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->Init();
}

//________________________________________________________________________
void AliEmcalJetTask::PrepareUtilities()
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->Prepare(fFastJetWrapper);
}

//________________________________________________________________________
void AliEmcalJetTask::ExecuteUtilities(AliEmcalJet* jet, Int_t ij)
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->ProcessJet(jet, ij, fFastJetWrapper);
}

//________________________________________________________________________
void AliEmcalJetTask::TerminateUtilities()
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->Terminate(fFastJetWrapper);
}

//________________________________________________________________________
void AliEmcalJetTask::UserCreateOutputObjects()
{
  // Create user objects.
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

  // clear the jet array (normally a null operation)
  fJets->Delete();

  // get primary vertex
  if(fEvent->GetPrimaryVertex()) fEvent->GetPrimaryVertex()->GetXYZ(fVertex);

  Int_t n = FindJets();

  if (n == 0) return;

  FillJetBranch();
}

//________________________________________________________________________
Int_t AliEmcalJetTask::FindJets()
{
  // Find jets.

  if (!fTracks && !fClus){
    AliError("No tracks or clusters, returning.");
    return 0;
  }

  fFastJetWrapper.Clear();

  AliDebug(2,Form("Jet type = %d", fJetType));

  if (fTracks) {
    const Int_t Ntracks = fTracks->GetEntries();
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliVParticle *t = static_cast<AliVParticle*>(fTracks->At(iTracks));

      if (!t) continue;

      if (fFilterHybridTracks) {  // the cast is safe because fFilterHybridTracks is reset in DoInit if the object type is not AliAODTrack
        AliAODTrack* aodTrack = static_cast<AliAODTrack*>(t);
        if (!aodTrack->IsHybridGlobalConstrainedGlobal()) continue;
      }

      if (t->Pt() < fMinJetTrackPt) continue;
      Double_t eta = t->Eta();
      Double_t phi = t->Phi();
      if ((eta<fEtaMin) || (eta>fEtaMax) ||
          (phi<fPhiMin) || (phi>fPhiMax))
        continue;

      if (((fJetType & kChargedJet) != 0) && (t->Charge() == 0)) {
        AliDebug(2,Form("Skipping track %d because it is neutral.", iTracks));
        continue;
      }

      if (((fJetType & kNeutralJet) != 0) && (t->Charge() != 0)) {
        AliDebug(2,Form("Skipping track %d because it is charged.", iTracks));
        continue;
      }

      Int_t lab = TMath::Abs(t->GetLabel());
      if (lab < fMinLabelTracks || lab > fMaxLabelTracks) {
        AliDebug(2,Form("Skipping track %d because label %d is not in range (%d, %d)", iTracks, lab, fMinLabelTracks, fMaxLabelTracks));
        continue;
      }

      if ((t->GetFlag() & fMCFlag) != fMCFlag) {
        AliDebug(2,Form("Skipping track %d because it does not match the MC flags", iTracks));
        continue;
      }

      if (fGeneratorIndex >= 0 && t->GetGeneratorIndex() != fGeneratorIndex) {
        AliDebug(2,Form("Skipping track %d because it does not match the MC generator index", iTracks));
        continue;
      }

      // artificial inefficiency
      if (fTrackEfficiency < 1.) {
        Double_t rnd = gRandom->Rndm();
        if (fTrackEfficiency < rnd) {
          AliDebug(2,Form("Track %d rejected due to artificial tracking inefficiency", iTracks));
          continue;
        }
      }

      // offset of 100 for consistency with cluster ids
      AliDebug(2,Form("Track %d accepted (label = %d, pt = %f)", iTracks, lab, t->Pt()));
      fFastJetWrapper.AddInputVector(t->Px(), t->Py(), t->Pz(), t->E(), iTracks + 100);
    }
  }

  if (fClus) {
    const Int_t Nclus = fClus->GetEntries();
    for (Int_t iClus = 0; iClus < Nclus; ++iClus) {
      AliVCluster *c = 0;
      Double_t cEta=0,cPhi=0,cPt=0;
      Double_t cPx=0,cPy=0,cPz=0;
      if (fIsEmcPart) {
        AliEmcalParticle *ep = static_cast<AliEmcalParticle*>(fClus->At(iClus));
        if (!ep) continue;
        c = ep->GetCluster();
        if (!c) continue;

        cEta = ep->Eta();
        cPhi = ep->Phi();
        cPt  = ep->Pt();
        cPx  = ep->Px();
        cPy  = ep->Py();
        cPz  = ep->Pz();
      } 
      else {
        c = static_cast<AliVCluster*>(fClus->At(iClus));
        if (!c) continue;

        TLorentzVector nP;
        if (fClusterEnergyType >= 0 &&  fClusterEnergyType <= AliVCluster::kLastUserDefEnergy) {
          c->GetMomentum(nP, fVertex, (AliVCluster::VCluUserDefEnergy_t)fClusterEnergyType);
        }
        else {
          c->GetMomentum(nP, fVertex);
        }
        cEta = nP.Eta();
        cPhi = TVector2::Phi_0_2pi(nP.Phi());
        cPt  = nP.Pt();
        cPx  = nP.Px();
        cPy  = nP.Py();
        cPz  = nP.Pz();
      }
      Double_t e = TMath::Sqrt(cPx*cPx+cPy*cPy+cPz*cPz);
      AliDebug(2,Form("Cluster %d (label = %d, energy = %.3f)", iClus, c->GetLabel(), e));

      if (!c->IsEMCAL()) continue;
      if (cPt < fMinJetClusPt) continue;
      if (e < fMinJetClusE) continue;
      if ((cEta<fEtaMin) || (cEta>fEtaMax) ||
          (cPhi<fPhiMin) || (cPhi>fPhiMax))
        continue;
      if (c->GetIsExotic()) continue;

      Int_t lab = TMath::Abs(c->GetLabel());
      if (lab < fMinLabelClusters || lab > fMaxLabelClusters) {
        AliDebug(2,Form("Skipping cluster %d because label %d is not in range (%d, %d)", iClus, lab, fMinLabelClusters, fMaxLabelClusters));
        continue;
      }

      // offset of 100 to skip ghost particles uid = -1
      AliDebug(2,Form("Cluster %d accepted (label = %d, energy = %.3f)", iClus, c->GetLabel(), e));
      fFastJetWrapper.AddInputVector(cPx, cPy, cPz, e, -iClus - 100);
    }
  }

  if (fFastJetWrapper.GetInputVectors().size() == 0) return 0;

  // run jet finder
  fFastJetWrapper.Run();

  return fFastJetWrapper.GetInclusiveJets().size();
}

//________________________________________________________________________
void AliEmcalJetTask::FillJetBranch()
{
  // Fill jet branch with fastjet jets.

  PrepareUtilities();

  // loop over fastjet jets
  std::vector<fastjet::PseudoJet> jets_incl = fFastJetWrapper.GetInclusiveJets();
  // sort jets according to jet pt
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, jets_incl);

  AliDebug(1,Form("%d jets found", (Int_t)jets_incl.size()));
  for (UInt_t ijet = 0, jetCount = 0; ijet < jets_incl.size(); ++ijet) {
    Int_t ij = indexes[ijet];
    AliDebug(3,Form("Jet pt = %f, area = %f", jets_incl[ij].perp(), fFastJetWrapper.GetJetArea(ij)));

    if (jets_incl[ij].perp() < fMinJetPt) continue;
    if (fFastJetWrapper.GetJetArea(ij) < fMinJetArea) continue;
    if ((jets_incl[ij].eta() < fJetEtaMin) || (jets_incl[ij].eta() > fJetEtaMax) ||
        (jets_incl[ij].phi() < fJetPhiMin) || (jets_incl[ij].phi() > fJetPhiMax))
      continue;

    AliEmcalJet *jet = new ((*fJets)[jetCount])
    		          AliEmcalJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());
    jet->SetLabel(ij);

    fastjet::PseudoJet area(fFastJetWrapper.GetJetAreaVector(ij));
    jet->SetArea(area.perp());
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());

    // Fill constituent info
    std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper.GetJetConstituents(ij));
    FillJetConstituents(jet, constituents, fTracks, fClus, constituents);

    if (fGeom) {
      if ((jet->Phi() > fGeom->GetArm1PhiMin() * TMath::DegToRad()) &&
          (jet->Phi() < fGeom->GetArm1PhiMax() * TMath::DegToRad()) &&
          (jet->Eta() > fGeom->GetArm1EtaMin()) &&
          (jet->Eta() < fGeom->GetArm1EtaMax()))
        jet->SetAxisInEmcal(kTRUE);
    }

    ExecuteUtilities(jet, ij);

    AliDebug(2,Form("Added jet n. %d, pt = %f, area = %f, constituents = %d", jetCount, jet->Pt(), jet->Area(), jet->GetNumberOfConstituents()));
    jetCount++;
  }

  TerminateUtilities();
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

  if (fTrackEfficiency < 1.) {
    if (gRandom) delete gRandom;
    gRandom = new TRandom3(0);
  }

  // get geometry
  fGeom = AliEMCALGeometry::GetInstance();
  if (!fGeom) {
    AliWarning(Form("%s: Can not create EMCal geometry, some features will not work...", GetName()));
  }

  // get input collections
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();

  // get the event
  fEvent = InputEvent();
  if (!fEvent) {
    AliError(Form("%s: Could not retrieve event! Returning", GetName()));
    return 0;
  }

  if (!fTracks) {
    if (fUseExchangeCont > 0) {
      fTracks = dynamic_cast<TClonesArray*>(GetInputData(1));
      if (!fTracks) {
        AliError(Form("%s: Could not get tracks form the input container n. 1", GetName()));
        return 0;
      }
    }
    else if (!fTracksName.IsNull()) {
      if (fTracksName == "Tracks") am->LoadBranch("Tracks");
      fTracks = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fTracksName));
      if (!fTracks) {
        AliError(Form("%s: Pointer to tracks %s == 0", GetName(), fTracksName.Data()));
        return 0;
      }
    }
  }

  if (fTracks) {
    TClass cls(fTracks->GetClass()->GetName());

    if (cls.InheritsFrom("AliVParticle")) {
      if (fFilterHybridTracks) {
        if (!cls.InheritsFrom("AliAODTrack")) {
          AliWarning(Form("Track collection contains objects of type '%s' that does not derive from 'AliAODTrack'. The hybrid track filter will not be applied.", cls.GetName()));
          fFilterHybridTracks = kFALSE;
        }
      }
    }
    else {
      AliError(Form("Track collection contains objects of type '%s' that does not derive from 'AliVParticle'. This input collection will be ignored!", cls.GetName()));
      fTracks = 0;
    }
  }

  if (!fClus) {
    if (fUseExchangeCont > 1) {
      fClus = dynamic_cast<TClonesArray*>(GetInputData(2));
      if (!fClus) {
        AliError(Form("%s: Could not get clusters form the input container n. 2", GetName()));
        return 0;
      }
    }
    else if (!fCaloName.IsNull()) {
      if (fCaloName == "CaloClusters") am->LoadBranch("CaloClusters");
      fClus = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fCaloName));
      if (!fClus) {
        AliError(Form("%s: Pointer to clus %s == 0", GetName(), fCaloName.Data()));
        return 0;
      }
    }
  }

  if (fClus) {
    TClass cls(fClus->GetClass()->GetName());
    if (cls.InheritsFrom("AliEmcalParticle")) {
      fIsEmcPart = 1;
    }
    else if (cls.InheritsFrom("AliVCluster")) {
      fIsEmcPart = 0;
    }
    else {
      AliError(Form("Cluster collection contains objects of type '%s' that does not derive from 'AliEmcalParticle' or 'AliVCluster'. This input collection will be ignored!", cls.GetName()));
      fTracks = 0;
    }
  }

  // add jets to event if not yet there
  if (!(fEvent->FindListObject(fJetsName))) {
    fJets = new TClonesArray("AliEmcalJet");
    fJets->SetName(fJetsName);
    fEvent->AddObject(fJets);
  }
  else {
    AliError(Form("%s: Object with name %s already in event! Returning", GetName(), fJetsName.Data()));
    return 0;
  }

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
  fFastJetWrapper.SetName(name);
  fFastJetWrapper.SetTitle(name);
  fFastJetWrapper.SetAreaType(fastjet::active_area_explicit_ghosts);
  fFastJetWrapper.SetGhostArea(fGhostArea);
  fFastJetWrapper.SetR(fRadius);
  fFastJetWrapper.SetAlgorithm(jalgo);
  fFastJetWrapper.SetRecombScheme(static_cast<fastjet::RecombinationScheme>(fRecombScheme));
  fFastJetWrapper.SetMaxRap(TMath::Max(TMath::Abs(fEtaMin),TMath::Abs(fEtaMax)));

  // setting legacy mode
  if (fLegacyMode) {
    fFastJetWrapper.SetLegacyMode(kTRUE);
  }

  InitUtilities();

  return kTRUE;
}

//___________________________________________________________________________________________________________________
void AliEmcalJetTask::FillJetConstituents(AliEmcalJet *jet, std::vector<fastjet::PseudoJet>& constituents, TClonesArray *tracks, TClonesArray *clusters,
    std::vector<fastjet::PseudoJet>& constituents_unsub, Int_t flag, TClonesArray *particles_sub)
{
  Int_t nt            = 0;
  Int_t nc            = 0;
  Double_t neutralE   = 0.;
  Double_t maxCh      = 0.;
  Double_t maxNe      = 0.;
  Int_t gall          = 0;
  Int_t gemc          = 0;
  Int_t cemc          = 0;
  Int_t ncharged      = 0;
  Int_t nneutral      = 0;
  Double_t mcpt       = 0.;
  Double_t emcpt      = 0.;

  Int_t uid   = -1;

  jet->SetNumberOfTracks(constituents.size());
  jet->SetNumberOfClusters(constituents.size());

  for (UInt_t ic = 0; ic < constituents.size(); ++ic) {

    if (flag == 0) {
      uid = constituents[ic].user_index();
    }
    else {
      if (constituents[ic].perp()<1.e-10) {
        uid=-1;
      }
      else {
        uid = GetIndexSub(constituents[ic].phi(), constituents_unsub);
      }
      if (uid==0) {
        AliError("correspondence between un/subtracted constituent not found");
        continue;
      }
    }

    AliDebug(3,Form("Processing constituent %d", uid));
    if (uid == -1) { //ghost particle
      ++gall;
      if (fGeom) {
        Double_t gphi = constituents[ic].phi();
        if (gphi < 0) gphi += TMath::TwoPi();
        gphi *= TMath::RadToDeg();
        Double_t geta = constituents[ic].eta();
        if ((gphi > fGeom->GetArm1PhiMin()) && (gphi < fGeom->GetArm1PhiMax()) &&
            (geta > fGeom->GetArm1EtaMin()) && (geta < fGeom->GetArm1EtaMax()))
          ++gemc;
      }

      if (fFillGhost) jet->AddGhost(constituents[ic].px(),
          constituents[ic].py(),
          constituents[ic].pz(),
          constituents[ic].e());
    }	
    else if ((uid >= 100) && tracks) { // track constituent
      Int_t tid = uid - 100;
      AliVParticle *t = static_cast<AliVParticle*>(tracks->At(tid));
      if (!t) {
        AliError(Form("Could not find track %d",tid));
        continue;
      }

      Double_t cEta = t->Eta();
      Double_t cPhi = t->Phi();
      Double_t cPt  = t->Pt();
      Double_t cP   = t->P();
      if (t->Charge() == 0) {
        neutralE += cP;
        ++nneutral;
        if (cPt > maxNe) maxNe = cPt;
      } else {
        ++ncharged;
        if (cPt > maxCh) maxCh = cPt;
      }

      // check if MC particle
      if (TMath::Abs(t->GetLabel()) > fMinMCLabel) mcpt += cPt;

      if (fGeom) {
        if (cPhi < 0) cPhi += TMath::TwoPi();
        cPhi *= TMath::RadToDeg();
        if ((cPhi > fGeom->GetArm1PhiMin()) && (cPhi < fGeom->GetArm1PhiMax()) &&
            (cEta > fGeom->GetArm1EtaMin()) && (cEta < fGeom->GetArm1EtaMax())) {
          emcpt += cPt;
          ++cemc;
        }
      }

      if (flag == 0 || particles_sub == 0) {
        jet->AddTrackAt(tid, nt);
      }
      else {
        Int_t part_sub_id = particles_sub->GetEntriesFast();
        AliEmcalParticle* part_sub = new ((*particles_sub)[part_sub_id]) AliEmcalParticle(dynamic_cast<AliVTrack*>(t));   // SA: probably need to be fixed!!
        part_sub->SetPtEtaPhiM(constituents[ic].perp(),constituents[ic].eta(),constituents[ic].phi(),constituents[ic].m());
        jet->AddTrackAt(part_sub_id, nt);
      }

      ++nt;
    } 
    else if ((uid <= -100) && clusters) { // cluster constituent
      Int_t cid = -uid - 100;
      AliVCluster *c = 0;
      Double_t cEta=0,cPhi=0,cPt=0,cP=0;
      if (fIsEmcPart) {
        AliEmcalParticle *ep = static_cast<AliEmcalParticle*>(fClus->At(cid));
        if (!ep) continue;
        c = ep->GetCluster();
        if (!c) continue;

        cEta = ep->Eta();
        cPhi = ep->Phi();
        cPt  = ep->Pt();
        cP   = ep->P();
      } 
      else {
        c = static_cast<AliVCluster*>(fClus->At(cid));
        if (!c) continue;

        TLorentzVector nP;
        c->GetMomentum(nP, fVertex);
        cEta = nP.Eta();
        cPhi = nP.Phi();
        cPt  = nP.Pt();
        cP   = nP.P();
      }

      neutralE += cP;
      if (cPt > maxNe) maxNe = cPt;

      // MC particle
      if (TMath::Abs(c->GetLabel()) > fMinMCLabel) mcpt += c->GetMCEnergyFraction() > 1e-6 ? cPt * c->GetMCEnergyFraction() : cPt;

      if (fGeom) {
        if (cPhi<0) cPhi += TMath::TwoPi();
        cPhi *= TMath::RadToDeg();
        if ((cPhi > fGeom->GetArm1PhiMin()) && (cPhi < fGeom->GetArm1PhiMax()) &&
            (cEta > fGeom->GetArm1EtaMin()) && (cEta < fGeom->GetArm1EtaMax())) {
          emcpt += cPt;
          ++cemc;
        }
      }

      if (flag == 0 || particles_sub == 0) {
        jet->AddClusterAt(cid, nc);
      }
      else {
        Int_t part_sub_id = particles_sub->GetEntriesFast();
        AliEmcalParticle* part_sub = new ((*particles_sub)[part_sub_id]) AliEmcalParticle(c);
        part_sub->SetPtEtaPhiM(constituents[ic].perp(),constituents[ic].eta(),constituents[ic].phi(),constituents[ic].m());
        jet->AddTrackAt(part_sub_id, nt);
      }

      ++nc;
      ++nneutral;
    } 
    else {
      AliError(Form("%s: No logical way to end up here.", GetName()));
      continue;
    }
  }

  jet->SetNumberOfTracks(nt);
  jet->SetNumberOfClusters(nc);
  jet->SetNEF(neutralE / jet->E());
  jet->SetMaxChargedPt(maxCh);
  jet->SetMaxNeutralPt(maxNe);
  if (gall > 0) jet->SetAreaEmc(jet->Area() * gemc / gall);
  else jet->SetAreaEmc(-1);
  jet->SetNEmc(cemc);
  jet->SetNumberOfCharged(ncharged);
  jet->SetNumberOfNeutrals(nneutral);
  jet->SetMCPt(mcpt);
  jet->SetPtEmc(emcpt);
  jet->SortConstituents();
}

//______________________________________________________________________________________
Int_t AliEmcalJetTask::GetIndexSub(Double_t phi_sub, std::vector<fastjet::PseudoJet>& constituents_unsub) 
{
  Double_t dphi=0;
  for (UInt_t ii = 0; ii < constituents_unsub.size(); ii++) {
    Double_t phi_unsub = constituents_unsub[ii].phi();
    dphi=phi_unsub-phi_sub;
    if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
    else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
    if (TMath::Abs(dphi)<0.1 && constituents_unsub[ii].user_index()!=-1)  return constituents_unsub[ii].user_index();
  }

  return 0;
}

//______________________________________________________________________________________
Bool_t AliEmcalJetTask::IsLocked() const
{
  if (fLocked) {
    AliFatal("Jet finder task is locked! Changing properties is not allowed."); 
    return kTRUE;
  } 
  else {
    return kFALSE;
  }
}

//______________________________________________________________________________________
void AliEmcalJetTask::SelectCollisionCandidates(UInt_t offlineTriggerMask)
{
  if(!fIsPSelSet) {
    fIsPSelSet = kTRUE;
    fOfflineTriggerMask = offlineTriggerMask;
  }
  else {
    AliWarning("Manually setting the event selection for jet finders is not allowed! Using trigger=old_trigger | your_trigger");  
    fOfflineTriggerMask = fOfflineTriggerMask | offlineTriggerMask;
  }
}

//______________________________________________________________________________________
void AliEmcalJetTask::SetRadius(Double_t r)            
{ 
  if (IsLocked()) return; 
  fRadius = r; 
  if ((fJetType & (kRX1Jet|kRX2Jet|kRX3Jet)) == 0) {
    AliWarning("Radius value will be ignored if jet type is not set to a user defined radius (kRX1Jet,kRX2Jet,kRX3Jet).");
  }
}

//______________________________________________________________________________________
void AliEmcalJetTask::SetType(Int_t t)                 
{ 
  if(IsLocked()) return; 
  if (t==0) fJetType |= kFullJet; 
  else if (t==1) fJetType |= kChargedJet; 
  else if (t==2) fJetType |= kNeutralJet; 
} // for backward compatibility only

//______________________________________________________________________________________
void AliEmcalJetTask::SelectPhysicalPrimaries(Bool_t s)    
{ 
  if(IsLocked()) return; 
  if (s) fMCFlag |=  AliAODMCParticle::kPhysicalPrim; 
  else   fMCFlag &= ~AliAODMCParticle::kPhysicalPrim; 
}
