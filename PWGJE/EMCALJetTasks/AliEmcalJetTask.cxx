//
// Emcal jet finder task.
//
// Authors: C.Loizides, S.Aiola, M. Verweij

#include <vector>
#include "AliEmcalJetTask.h"

#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TClass.h>

#include "AliTLorentzVector.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliFJWrapper.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliEmcalJetUtility.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

using std::cout;
using std::endl;
using std::cerr;

ClassImp(AliEmcalJetTask)

const Int_t AliEmcalJetTask::fgkConstIndexShift = 100000;

//________________________________________________________________________
AliEmcalJetTask::AliEmcalJetTask() :
  AliAnalysisTaskEmcal(),
  fJetsTag(),
  fJetAlgo(AliJetContainer::antikt_algorithm),
  fJetType(AliJetContainer::kFullJet),
  fRecombScheme(AliJetContainer::pt_scheme),
  fRadius(0.4),
  fMinJetArea(0.001),
  fMinJetPt(1.0),
  fJetPhiMin(-10),
  fJetPhiMax(+10),
  fJetEtaMin(-1),
  fJetEtaMax(+1),
  fGhostArea(0.005),
  fTrackEfficiency(1.),
  fUtilities(0),
  fLocked(0),
  fJetsName(),
  fIsInit(0),
  fIsPSelSet(0),
  fIsEmcPart(0),
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fJets(0),
  fFastJetWrapper("AliEmcalJetTask","AliEmcalJetTask")
{
  // Default constructor.
}

//________________________________________________________________________
AliEmcalJetTask::AliEmcalJetTask(const char *name) :
  AliAnalysisTaskEmcal(name),
  fJetsTag("Jets"),
  fJetAlgo(AliJetContainer::antikt_algorithm),
  fJetType(AliJetContainer::kFullJet),
  fRecombScheme(AliJetContainer::pt_scheme),
  fRadius(0.4),
  fMinJetArea(0.001),
  fMinJetPt(1.0),
  fJetPhiMin(-10),
  fJetPhiMax(+10),
  fJetEtaMin(-1),
  fJetEtaMax(+1),
  fGhostArea(0.005),
  fTrackEfficiency(1.),
  fUtilities(0),
  fLocked(0),
  fJetsName(),
  fIsInit(0),
  fIsPSelSet(0),
  fIsEmcPart(0),
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fJets(0),
  fFastJetWrapper(name,name)
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
Bool_t AliEmcalJetTask::Run()
{
  // Main loop, called for each event.

  // clear the jet array (normally a null operation)
  fJets->Delete();

  Int_t n = FindJets();

  if (n == 0) return kFALSE;

  FillJetBranch();

  return kTRUE;
}

//________________________________________________________________________
Int_t AliEmcalJetTask::FindJets()
{
  // Find jets.


  if (fParticleCollArray.GetEntriesFast() == 0 && fClusterCollArray.GetEntriesFast() == 0){
    AliError("No tracks or clusters, returning.");
    return 0;
  }

  fFastJetWrapper.Clear();

  AliDebug(2,Form("Jet type = %d", fJetType));

  AliTLorentzVector mom;

  Int_t iColl = 1;
  TIter nextPartColl(&fParticleCollArray);
  AliParticleContainer* tracks = 0;
  while ((tracks = static_cast<AliParticleContainer*>(nextPartColl()))) {
    AliDebug(2,Form("Tracks from collection %d: '%s'.", iColl-1, tracks->GetName()));
    tracks->ResetCurrentID();
    AliVParticle* t = 0;
    while ((t = tracks->GetNextAcceptParticle())) {
      tracks->GetMomentum(mom, tracks->GetCurrentID());
      if (((fJetType & AliJetContainer::kChargedJet) != 0) && (t->Charge() == 0)) {
        AliDebug(2,Form("Skipping track %d because it is neutral.", tracks->GetCurrentID()));
        continue;
      }

      if (((fJetType & AliJetContainer::kNeutralJet) != 0) && (t->Charge() != 0)) {
        AliDebug(2,Form("Skipping track %d because it is charged.", tracks->GetCurrentID()));
        continue;
      }

      // artificial inefficiency
      if (fTrackEfficiency < 1.) {
        Double_t rnd = gRandom->Rndm();
        if (fTrackEfficiency < rnd) {
          AliDebug(2,Form("Track %d rejected due to artificial tracking inefficiency", tracks->GetCurrentID()));
          continue;
        }
      }

      AliDebug(2,Form("Track %d accepted (label = %d, pt = %f)", tracks->GetCurrentID(), t->GetLabel(), t->Pt()));
      Int_t uid = tracks->GetCurrentID() + fgkConstIndexShift * iColl;
      fFastJetWrapper.AddInputVector(mom.Px(), mom.Py(), mom.Pz(), mom.E(), uid);
    }
    iColl++;
  }

  iColl = 1;
  TIter nextClusColl(&fClusterCollArray);
  AliClusterContainer* clusters = 0;
  while ((clusters = static_cast<AliClusterContainer*>(nextClusColl()))) {
    AliDebug(2,Form("Clusters from collection %d: '%s'.", iColl-1, clusters->GetName()));
    clusters->ResetCurrentID();
    AliVCluster* c = 0;
    while ((c = clusters->GetNextAcceptCluster())) {
      clusters->GetMomentum(mom, clusters->GetCurrentID());
      Double_t cEta = mom.Eta();
      Double_t cPhi = mom.Phi_0_2pi();
      Double_t cPt  = mom.Pt();
      Double_t cPx  = mom.Px();
      Double_t cPy  = mom.Py();
      Double_t cPz  = mom.Pz();

      Double_t e = TMath::Sqrt(cPx*cPx+cPy*cPy+cPz*cPz);

      AliDebug(2,Form("Cluster %d accepted (label = %d, energy = %.3f)", clusters->GetCurrentID(), c->GetLabel(), e));
      Int_t uid = -clusters->GetCurrentID() - fgkConstIndexShift * iColl;
      fFastJetWrapper.AddInputVector(cPx, cPy, cPz, e, uid);
    }
    iColl++;
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
    FillJetConstituents(jet, constituents, constituents);

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
void AliEmcalJetTask::ExecOnce()
{
  // Init the task.

  SetNeedEmcalGeom(kFALSE);

  if (fTrackEfficiency < 1.) {
    if (gRandom) delete gRandom;
    gRandom = new TRandom3(0);
  }

  fJetsName = AliJetContainer::GenerateJetName(fJetType, fJetAlgo, fRecombScheme, fRadius, GetParticleContainer(0), GetClusterContainer(0), fJetsTag);

  // add jets to event if not yet there
  if (!(InputEvent()->FindListObject(fJetsName))) {
    fJets = new TClonesArray("AliEmcalJet");
    fJets->SetName(fJetsName);
    ::Info("AliEmcalJetTask::ExecOnce", "Jet collection with name '%s' has been added to the event.", fJetsName.Data());
    InputEvent()->AddObject(fJets);
  }
  else {
    AliError(Form("%s: Object with name %s already in event! Returning", GetName(), fJetsName.Data()));
    return;
  }

  // setup fj wrapper
  fFastJetWrapper.SetAreaType(fastjet::active_area_explicit_ghosts);
  fFastJetWrapper.SetGhostArea(fGhostArea);
  fFastJetWrapper.SetR(fRadius);
  fFastJetWrapper.SetAlgorithm(ConvertToFJAlgo(fJetAlgo));
  fFastJetWrapper.SetRecombScheme(ConvertToFJRecoScheme(fRecombScheme));
  fFastJetWrapper.SetMaxRap(1);

  // setting legacy mode
  if (fLegacyMode) {
    fFastJetWrapper.SetLegacyMode(kTRUE);
  }

  InitUtilities();

  AliAnalysisTaskEmcal::ExecOnce();
}

//___________________________________________________________________________________________________________________
void AliEmcalJetTask::FillJetConstituents(AliEmcalJet *jet, std::vector<fastjet::PseudoJet>& constituents,
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
    else if (uid >= fgkConstIndexShift) { // track constituent
      Int_t iColl = uid / fgkConstIndexShift;
      Int_t tid = uid - iColl * fgkConstIndexShift;
      iColl--;
      AliDebug(3,Form("Constituent %d is a track from collection %d and with ID %d", uid, iColl, tid));
      AliParticleContainer* partCont = GetParticleContainer(iColl);
      if (!partCont) {
        AliError(Form("Could not find particle container %d",iColl));
        continue;
      }
      AliVParticle *t = partCont->GetParticle(tid);
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
    else if (uid <= -fgkConstIndexShift) { // cluster constituent
      Int_t iColl = -uid / fgkConstIndexShift;
      Int_t cid = -uid - iColl * fgkConstIndexShift;
      iColl--;
      AliClusterContainer* clusCont = GetClusterContainer(iColl);
      AliVCluster *c = clusCont->GetCluster(cid);

      if (!c) continue;

      AliTLorentzVector nP;
      clusCont->GetMomentum(nP, cid);

      Double_t cEta = nP.Eta();
      Double_t cPhi = nP.Phi_0_2pi();
      Double_t cPt  = nP.Pt();
      Double_t cP   = nP.P();

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
      AliError(Form("%s: No logical way to end up here (uid = %d).", GetName(), uid));
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
  Double_t phimin=100;
  Int_t index=0;
  for (UInt_t ii = 0; ii < constituents_unsub.size(); ii++) {
    dphi=0;
    Double_t phi_unsub = constituents_unsub[ii].phi();
    dphi=phi_unsub-phi_sub;
    if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
    else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
    if(TMath::Abs(dphi)<phimin){ phimin=TMath::Abs(dphi);
      index=ii;} }
    if (constituents_unsub[index].user_index()!=-1)  return constituents_unsub[index].user_index();

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
void AliEmcalJetTask::SetEtaRange(Double_t emi, Double_t ema)
{
  if (IsLocked()) return;

  TIter nextPartColl(&fParticleCollArray);
  AliParticleContainer* tracks = 0;
  while ((tracks = static_cast<AliParticleContainer*>(nextPartColl()))) {
    tracks->SetParticleEtaLimits(emi, ema);
  }
}

//______________________________________________________________________________________
void AliEmcalJetTask::SetMinJetClusPt(Double_t min)
{
  if (IsLocked()) return;

  TIter nextClusColl(&fClusterCollArray);
  AliClusterContainer* clusters = 0;
  while ((clusters = static_cast<AliClusterContainer*>(nextClusColl()))) {
    clusters->SetClusPtCut(min);
  }
}

//______________________________________________________________________________________
void AliEmcalJetTask::SetMinJetClusE(Double_t min)
{
  if (IsLocked()) return;

  TIter nextClusColl(&fClusterCollArray);
  AliClusterContainer* clusters = 0;
  while ((clusters = static_cast<AliClusterContainer*>(nextClusColl()))) {
    clusters->SetClusECut(min);
  }
}

//______________________________________________________________________________________
void AliEmcalJetTask::SetMinJetTrackPt(Double_t min)
{
  if (IsLocked()) return;

  TIter nextPartColl(&fParticleCollArray);
  AliParticleContainer* tracks = 0;
  while ((tracks = static_cast<AliParticleContainer*>(nextPartColl()))) {
    tracks->SetParticlePtCut(min);
  }
}

//______________________________________________________________________________________
void AliEmcalJetTask::SetPhiRange(Double_t pmi, Double_t pma)
{
  if (IsLocked()) return;

  TIter nextPartColl(&fParticleCollArray);
  AliParticleContainer* tracks = 0;
  while ((tracks = static_cast<AliParticleContainer*>(nextPartColl()))) {
    tracks->SetParticlePhiLimits(pmi, pma);
  }
}

//______________________________________________________________________________________
fastjet::JetAlgorithm AliEmcalJetTask::ConvertToFJAlgo(EJetAlgo_t algo)
{
  switch(algo) {
  case AliJetContainer::kt_algorithm:
    return fastjet::kt_algorithm;
  case AliJetContainer::antikt_algorithm:
    return fastjet::antikt_algorithm;
  case AliJetContainer::cambridge_algorithm:
    return fastjet::cambridge_algorithm;
  case AliJetContainer::genkt_algorithm:
    return fastjet::genkt_algorithm;
  case AliJetContainer::cambridge_for_passive_algorithm:
    return fastjet::cambridge_for_passive_algorithm;
  case AliJetContainer::genkt_for_passive_algorithm:
    return fastjet::genkt_for_passive_algorithm;
  case AliJetContainer::plugin_algorithm:
    return fastjet::plugin_algorithm;
  case AliJetContainer::undefined_jet_algorithm:
    return fastjet::undefined_jet_algorithm;

  default:
    ::Error("AliEmcalJetTask::ConvertToFJAlgo", "Jet algorithm %d not recognized!!!", algo);
    return fastjet::undefined_jet_algorithm;
  }
}

//______________________________________________________________________________________
fastjet::RecombinationScheme AliEmcalJetTask::ConvertToFJRecoScheme(ERecoScheme_t reco)
{
  switch(reco) {
  case AliJetContainer::E_scheme:
    return fastjet::E_scheme;

  case AliJetContainer::pt_scheme:
    return fastjet::pt_scheme;

  case AliJetContainer::pt2_scheme:
    return fastjet::pt2_scheme;

  case AliJetContainer::Et_scheme:
    return fastjet::Et_scheme;

  case AliJetContainer::Et2_scheme:
    return fastjet::Et2_scheme;

  case AliJetContainer::BIpt_scheme:
    return fastjet::BIpt_scheme;

  case AliJetContainer::BIpt2_scheme:
    return fastjet::BIpt2_scheme;

  case AliJetContainer::external_scheme:
    return fastjet::external_scheme;

  default:
    ::Error("AliEmcalJetTask::ConvertToFJRecoScheme", "Recombination scheme %d not recognized!!!", reco);
    return fastjet::external_scheme;
  }
}
