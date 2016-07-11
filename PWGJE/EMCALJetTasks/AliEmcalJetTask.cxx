/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <vector>

#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom3.h>

#include <AliVCluster.h>
#include <AliVEvent.h>
#include <AliVParticle.h>
#include <AliEMCALGeometry.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliFJWrapper.h"
#include "AliEmcalJetUtility.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliEmcalJetTask.h"

using std::cout;
using std::endl;
using std::cerr;

/// \cond CLASSIMP
ClassImp(AliEmcalJetTask);
/// \endcond

const Int_t AliEmcalJetTask::fgkConstIndexShift = 100000;

/**
 * Default constructor. This constructor is only for ROOT I/O and
 * not to be used by users.
 */
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
}

/**
 * Standard named constructor.
 * @param name Name of the task.
 */
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
}

/**
 * Destructor
 */
AliEmcalJetTask::~AliEmcalJetTask()
{
}

/**
 * Add a utility to the utility list. Utilities are instances of classes
 * derived from AliEmcalJetUtility that implements wrappers to FastJet contribs.
 */
AliEmcalJetUtility* AliEmcalJetTask::AddUtility(AliEmcalJetUtility* utility)
{
  if (!fUtilities) fUtilities = new TObjArray();
  if (fUtilities->FindObject(utility)) {
    Error("AddUtility", "Jet utility %s already connected.", utility->GetName());
    return utility;
  }   
  fUtilities->Add(utility);
  utility->SetJetTask(this);

  return utility;
}

/**
 * This method is called once before analyzing the first event. It executes
 * the Init() method of all utilities (if any).
 */
void AliEmcalJetTask::InitUtilities()
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->Init();
}

/**
 * This method is called in the event loop after jet finding but before filling
 * the output jet branch to prepare the utilities.
 * It executes the Prepare() method of all utilities (if any).
 */
void AliEmcalJetTask::PrepareUtilities()
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->Prepare(fFastJetWrapper);
}

/**
 * This method is called in the event loop for each jet found, while filling the output jet branch.
 * It executes the ProcessJet() method of all utilities (if any).
 */
void AliEmcalJetTask::ExecuteUtilities(AliEmcalJet* jet, Int_t ij)
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->ProcessJet(jet, ij, fFastJetWrapper);
}

/**
 * This method is called in the event loop after jet finding has been completed.
 * It executes the Terminate() method of all utilities (if any).
 */
void AliEmcalJetTask::TerminateUtilities()
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->Terminate(fFastJetWrapper);
}

/**
 * This method is called for each event.
 * @return Always kTRUE
 */
Bool_t AliEmcalJetTask::Run()
{
  // clear the jet array (normally a null operation)
  fJets->Delete();

  Int_t n = FindJets();

  if (n == 0) return kFALSE;

  FillJetBranch();

  return kTRUE;
}

/**
 * This method steers the jet finding. It first loops over all particle and cluster containers
 * that were provided when the task was initialized. All accepted objects (tracks, particle, clusters)
 * are added as input vectors to the FastJet wrapper. Then the jet finding is launched
 * in the wrapper.
 * @return Total number of jets found.
 */
Int_t AliEmcalJetTask::FindJets()
{
  if (fParticleCollArray.GetEntriesFast() == 0 && fClusterCollArray.GetEntriesFast() == 0){
    AliError("No tracks or clusters, returning.");
    return 0;
  }

  fFastJetWrapper.Clear();

  AliDebug(2,Form("Jet type = %d", fJetType));

  Int_t iColl = 1;
  TIter nextPartColl(&fParticleCollArray);
  AliParticleContainer* tracks = 0;
  while ((tracks = static_cast<AliParticleContainer*>(nextPartColl()))) {
    AliDebug(2,Form("Tracks from collection %d: '%s'.", iColl-1, tracks->GetName()));
    AliParticleIterableMomentumContainer itcont = tracks->accepted_momentum();
    for (AliParticleIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
      // artificial inefficiency
      if (fTrackEfficiency < 1.) {
        Double_t rnd = gRandom->Rndm();
        if (fTrackEfficiency < rnd) {
          AliDebug(2,Form("Track %d rejected due to artificial tracking inefficiency", it.current_index()));
          continue;
        }
      }

      AliDebug(2,Form("Track %d accepted (label = %d, pt = %f)", it.current_index(), it->second->GetLabel(), it->first.Pt()));
      Int_t uid = it.current_index() + fgkConstIndexShift * iColl;
      fFastJetWrapper.AddInputVector(it->first.Px(), it->first.Py(), it->first.Pz(), it->first.E(), uid);
    }
    iColl++;
  }

  iColl = 1;
  TIter nextClusColl(&fClusterCollArray);
  AliClusterContainer* clusters = 0;
  while ((clusters = static_cast<AliClusterContainer*>(nextClusColl()))) {
    AliDebug(2,Form("Clusters from collection %d: '%s'.", iColl-1, clusters->GetName()));
    AliClusterIterableMomentumContainer itcont = clusters->accepted_momentum();
    for (AliClusterIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
      AliDebug(2,Form("Cluster %d accepted (label = %d, energy = %.3f)", it.current_index(), it->second->GetLabel(), it->first.E()));
      Int_t uid = -it.current_index() - fgkConstIndexShift * iColl;
      fFastJetWrapper.AddInputVector(it->first.Px(), it->first.Py(), it->first.Pz(), it->first.E(), uid);
    }
    iColl++;
  }

  if (fFastJetWrapper.GetInputVectors().size() == 0) return 0;

  // run jet finder
  fFastJetWrapper.Run();

  return fFastJetWrapper.GetInclusiveJets().size();
}

/**
 * This method fills the jet output branch (TClonesArray) with the jet found by the FastJet
 * wrapper. Before filling the jet branch, the utilities are prepared. Then the utilities are
 * called for each jet and finally after jet finding the terminate method of all utilities is called.
 */
void AliEmcalJetTask::FillJetBranch()
{
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
    jet->SetJetAcceptanceType(FindJetAcceptanceType(jet->Eta(), jet->Phi_0_2pi(), fRadius));

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

/**
 * Sorts jets by pT (decreasing)
 * @param[out] indexes This array is used to return the indexes of the jets ordered by pT
 * @param[in] array Vector containing the list of jets obtained by the FastJet wrapper
 * @return kTRUE if at least one jet was found in array; kFALSE otherwise
 */
Bool_t AliEmcalJetTask::GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const
{
  static Float_t pt[9999] = {0};

  const Int_t n = (Int_t)array.size();

  if (n < 1)
    return kFALSE;

  for (Int_t i = 0; i < n; i++)
    pt[i] = array[i].perp();

  TMath::Sort(n, pt, indexes);

  return kTRUE;
}

/**
 * This method is called once before analzying the first event.
 * It generates the output jet branch name, initializes the FastJet wrapper
 * and the utilities (FJ contribs).
 */
void AliEmcalJetTask::ExecOnce()
{
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

/**
 * This method is called for each jet. It loops over the jet constituents and
 * adds them to the jet object.
 * @param jet Pointer to the AliEmcalJet object where the jet constituents will be added
 * @param constituents List of the jet constituents returned by the FastJet wrapper
 * @param constituents_unsub List of jet constituents before background subtraction
 * @param flag If kTRUE it means that the argument "constituents" is a list of subtracted constituents
 * @param particles_sub Array containing subtracted constituents
 */
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

/**
 * Search for the index of the unsubtracted constituents by comparing the azimuthal angles
 * @param phi_sub Azimuthal angle of the subtracted constituent
 * @param constituents_unsub Vector containing the list of unsubtracted constituents
 * @return Index of the subtracted constituent
 */
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

/**
 * An instance of this class can be "locked". Once locked, it cannot be unlocked.
 * If the instance is locked, attempting to change the configuration will throw a
 * fatal and stop the execution of the program. This method checks whether the instance
 * is locked and throw a fatal if it is locked.
 */
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

/**
 * This overloads the method of AliAnalysisTaskSE to set the trigger bits. Since
 * the output of this task is often a shared input of several consumer task
 * the event selection is not allwed.
 * @param offlineTriggerMask Trigger bit mask
 */
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

/**
 * Set the eta range of the track constituents.
 * @param emi Minimum eta
 * @param ema Maximum eta
 */
void AliEmcalJetTask::SetEtaRange(Double_t emi, Double_t ema)
{
  if (IsLocked()) return;

  TIter nextPartColl(&fParticleCollArray);
  AliParticleContainer* tracks = 0;
  while ((tracks = static_cast<AliParticleContainer*>(nextPartColl()))) {
    tracks->SetParticleEtaLimits(emi, ema);
  }
}

/**
 * Set the minimum pT of the cluster constituents.
 * @param min Minimum pT
 */
void AliEmcalJetTask::SetMinJetClusPt(Double_t min)
{
  if (IsLocked()) return;

  TIter nextClusColl(&fClusterCollArray);
  AliClusterContainer* clusters = 0;
  while ((clusters = static_cast<AliClusterContainer*>(nextClusColl()))) {
    clusters->SetClusPtCut(min);
  }
}

/**
 * Set the minimum energy of the cluster constituents.
 * @param min Minimum energy
 */
void AliEmcalJetTask::SetMinJetClusE(Double_t min)
{
  if (IsLocked()) return;

  TIter nextClusColl(&fClusterCollArray);
  AliClusterContainer* clusters = 0;
  while ((clusters = static_cast<AliClusterContainer*>(nextClusColl()))) {
    clusters->SetClusECut(min);
  }
}

/**
 * Set the minimum pT of the track constituents.
 * @param min Minimum pT
 */
void AliEmcalJetTask::SetMinJetTrackPt(Double_t min)
{
  if (IsLocked()) return;

  TIter nextPartColl(&fParticleCollArray);
  AliParticleContainer* tracks = 0;
  while ((tracks = static_cast<AliParticleContainer*>(nextPartColl()))) {
    tracks->SetParticlePtCut(min);
  }
}

/**
 * Set the phi range of the track constituents.
 * @param pmi Minimum phi
 * @param pma Maximum phi
 */
void AliEmcalJetTask::SetPhiRange(Double_t pmi, Double_t pma)
{
  if (IsLocked()) return;

  TIter nextPartColl(&fParticleCollArray);
  AliParticleContainer* tracks = 0;
  while ((tracks = static_cast<AliParticleContainer*>(nextPartColl()))) {
    tracks->SetParticlePhiLimits(pmi, pma);
  }
}

/**
 * Converts the internal enum values representing jet algorithms in
 * the corresponding values accepted by the FastJet wrapper.
 * @param algo Algorithm represented in the EJetAlgo_t enum
 * @return Algortithm represented in the fastjet::JetAlgorithm enum
 */
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

/**
 * Converts the internal enum values representing jet recombination schemes in
 * the corresponding values accepted by the FastJet wrapper.
 * @param reco Recombination scheme represented in the EJetAlgo_t enum
 * @return Recombination scheme represented in the fastjet::JetAlgorithm enum
 */
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

/**
 * Finds which geometrical acceptance types the jet satisfies.
 * @return bitwise jet acceptance type
 */
UInt_t AliEmcalJetTask::FindJetAcceptanceType(Double_t eta, Double_t phi, Double_t r) {
  
  //This method has to be called after the run number is known because it needs the EMCal geometry object.
  
  UInt_t jetAcceptanceType = AliEmcalJet::kUser; // all jets satify the "no acceptance cut" condition
  
  // Check if TPC
  if( eta < 0.9 && eta > -0.9 ) {
    jetAcceptanceType |= AliEmcalJet::kTPC;
    // Check if TPCfid
    if (eta < 0.9 - r && eta > -0.9 + r)
      jetAcceptanceType |= AliEmcalJet::kTPCfid;
  }
    
  // Check if EMCAL
  if( IsJetInEmcal(eta, phi, 0) ) {
    jetAcceptanceType |= AliEmcalJet::kEMCAL;
    // Check if EMCALfid
    if( IsJetInEmcal(eta, phi, r) )
      jetAcceptanceType |= AliEmcalJet::kEMCALfid;
  }
  
  // Check if DCAL
  if( IsJetInDcal(eta, phi, 0) ) {
    jetAcceptanceType |= AliEmcalJet::kDCAL;
    // Check if DCALfid
    if( IsJetInDcal(eta, phi, r) )
      jetAcceptanceType |= AliEmcalJet::kDCALfid;
  }
 
  return jetAcceptanceType;
}

/**
 * Returns whether or not jet with given eta, phi, R is in EMCal.
 */
Bool_t AliEmcalJetTask::IsJetInEmcal(Double_t eta, Double_t phi, Double_t r)
{
  if (!fGeom) SetEMCALGeometry();
  if (fGeom) {
    if ( eta < fGeom->GetArm1EtaMax() - r && eta > fGeom->GetArm1EtaMin() + r ) {
      if(fRunNumber >= 177295 && fRunNumber <= 197470) {//small SM masked in 2012 and 2013
        if ( phi < 3.135 - r && phi > 1.405 + r )
          return kTRUE;
      }
      else {
        if ( phi < fGeom->GetEMCALPhiMax() * TMath::DegToRad() - r && phi > fGeom->GetArm1PhiMin() * TMath::DegToRad() + r)
          return kTRUE;
      }
    }
  }
  else {
    AliWarning("Could not get instance of AliEMCALGeometry. Using manual settings for EMCAL year 2011!!");
    if (eta < 0.7 - r && eta > -0.7 + r ){
      if (phi < 3.135 - r && phi > 1.405 + r )
        return kTRUE;
    }
  }
  return kFALSE;
}

/**
 * Returns whether or not jet with given eta, phi, R is in DCal.
 */
Bool_t AliEmcalJetTask::IsJetInDcal(Double_t eta, Double_t phi, Double_t r)
{
  if (!fGeom) SetEMCALGeometry();
  if (fGeom) {
    if (eta < fGeom->GetArm1EtaMax() - r && eta > fGeom->GetArm1EtaMin() + r ) {
      if ( phi < fGeom->GetDCALPhiMax() * TMath::DegToRad() - r && phi > fGeom->GetDCALPhiMin() * TMath::DegToRad() + r)
        return kTRUE;
    }
  }
  else {
    AliWarning("Could not get instance of AliEMCALGeometry. Using manual settings for DCAL year 2015!!");
    if (eta < 0.7 - r && eta > -0.7 + r) {
      if ( phi < 5.727 - r && phi > 4.538 + r )
        return kTRUE;
    }
  }
  return kFALSE;
}

/**
 * Set a pointer to the EMCal geometry object.
 * It assumes that an instance of the object has been already
 * configured.
 */
void AliEmcalJetTask::SetEMCALGeometry()
{
  fGeom = AliEMCALGeometry::GetInstance();
  if (!fGeom) {
    AliError(Form("%s: Can not create geometry", GetName()));
    return;
  }
}
