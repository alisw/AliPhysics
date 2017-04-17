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

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliMultiInputEventHandler.h>

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
  fTrackEfficiencyOnlyForEmbedding(kFALSE),
  fUtilities(0),
  fLocked(0),
  fJetsName(),
  fIsInit(0),
  fIsPSelSet(0),
  fIsEmcPart(0),
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fJets(0),
  fClusterContainerIndexMap(),
  fParticleContainerIndexMap(),
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
  fTrackEfficiencyOnlyForEmbedding(kFALSE),
  fUtilities(0),
  fLocked(0),
  fJetsName(),
  fIsInit(0),
  fIsPSelSet(0),
  fIsEmcPart(0),
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fJets(0),
  fClusterContainerIndexMap(),
  fParticleContainerIndexMap(),
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
 * This method is called before analyzing each event. It executes
 * the InitEvent() method of all utilities (if any).
 */
void AliEmcalJetTask::InitEvent()
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->InitEvent(fFastJetWrapper);
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
  InitEvent();
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
    AliDebug(2,Form("Tracks from collection %d: '%s'. Embedded: %i, nTracks: %i", iColl-1, tracks->GetName(), tracks->GetIsEmbedding(), tracks->GetNParticles()));
    AliParticleIterableMomentumContainer itcont = tracks->accepted_momentum();
    for (AliParticleIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
      // artificial inefficiency
      if (fTrackEfficiency < 1.) {
        if (fTrackEfficiencyOnlyForEmbedding == kFALSE || (fTrackEfficiencyOnlyForEmbedding == kTRUE && tracks->GetIsEmbedding())) {
          Double_t rnd = gRandom->Rndm();
          if (fTrackEfficiency < rnd) {
            AliDebug(2,Form("Track %d rejected due to artificial tracking inefficiency", it.current_index()));
            continue;
          }
        }
      }

      AliDebug(2,Form("Track %d accepted (label = %d, pt = %f, eta = %f, phi = %f, E = %f, m = %f, px = %f, py = %f, pz = %f)", it.current_index(), it->second->GetLabel(), it->first.Pt(), it->first.Eta(), it->first.Phi(), it->first.E(), it->first.M(), it->first.Px(), it->first.Py(), it->first.Pz()));
      Int_t uid = it.current_index() + fgkConstIndexShift * iColl;
      fFastJetWrapper.AddInputVector(it->first.Px(), it->first.Py(), it->first.Pz(), it->first.E(), uid);
    }
    iColl++;
  }

  iColl = 1;
  TIter nextClusColl(&fClusterCollArray);
  AliClusterContainer* clusters = 0;
  while ((clusters = static_cast<AliClusterContainer*>(nextClusColl()))) {
    AliDebug(2,Form("Clusters from collection %d: '%s'. Embedded: %i, nClusters: %i", iColl-1, clusters->GetName(), clusters->GetIsEmbedding(), clusters->GetNClusters()));
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

  // Setup container utils. Must be called after AliAnalysisTaskEmcal::ExecOnce() so that the
  // containers' arrays are setup.
  fClusterContainerIndexMap.CopyMappingFrom(AliClusterContainer::GetEmcalContainerIndexMap(), fClusterCollArray);
  fParticleContainerIndexMap.CopyMappingFrom(AliParticleContainer::GetEmcalContainerIndexMap(), fParticleCollArray);
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
    std::vector<fastjet::PseudoJet>& constituents_unsub, Int_t flag, TString particlesSubName)
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
  TClonesArray * particles_sub = 0;

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
        uid = constituents[ic].user_index();
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

      if (flag == 0 || particlesSubName == "") {
        jet->AddTrackAt(fParticleContainerIndexMap.GlobalIndexFromLocalIndex(partCont, tid), nt);
      }
      else {
        // Get the particle container and array corresponding to the subtracted particles
        partCont = GetParticleContainer(particlesSubName);
        particles_sub = partCont->GetArray();
        // Create the new particle in the particles_sub array and add it to the jet
        Int_t part_sub_id = particles_sub->GetEntriesFast();
        AliEmcalParticle* part_sub = new ((*particles_sub)[part_sub_id]) AliEmcalParticle(dynamic_cast<AliVTrack*>(t));   // SA: probably need to be fixed!!
        part_sub->SetPtEtaPhiM(constituents[ic].perp(),constituents[ic].eta(),constituents[ic].phi(),constituents[ic].m());
        jet->AddTrackAt(fParticleContainerIndexMap.GlobalIndexFromLocalIndex(partCont, part_sub_id), nt);
      }

      ++nt;
    } 
    else if (uid <= -fgkConstIndexShift) { // cluster constituent
      Int_t iColl = -uid / fgkConstIndexShift;
      Int_t cid = -uid - iColl * fgkConstIndexShift;
      iColl--;
      AliDebug(3,Form("Constituent %d is a cluster from collection %d and with ID %d", uid, iColl, cid));
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

      if (flag == 0 || particlesSubName == "") {
        jet->AddClusterAt(fClusterContainerIndexMap.GlobalIndexFromLocalIndex(clusCont, cid), nc);
      }
      else {
        // Get the cluster container and array corresponding to the subtracted particles
        clusCont = GetClusterContainer(particlesSubName);
        particles_sub = clusCont->GetArray();
        // Create the new particle in the particles_sub array and add it to the jet
        Int_t part_sub_id = particles_sub->GetEntriesFast();
        AliEmcalParticle* part_sub = new ((*particles_sub)[part_sub_id]) AliEmcalParticle(c);
        part_sub->SetPtEtaPhiM(constituents[ic].perp(),constituents[ic].eta(),constituents[ic].phi(),constituents[ic].m());
        jet->AddClusterAt(fClusterContainerIndexMap.GlobalIndexFromLocalIndex(clusCont, part_sub_id), nc);
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
  
  // Check if DCAL (i.e. eta-phi rectangle spanning DCal, which includes most of PHOS)
  if( IsJetInDcal(eta, phi, 0) ) {
    jetAcceptanceType |= AliEmcalJet::kDCAL;
    // Check if DCALfid
    if( IsJetInDcal(eta, phi, r) )
      jetAcceptanceType |= AliEmcalJet::kDCALfid;
  }
  
  // Check if DCALonly (i.e. ONLY DCal, does not include any of PHOS region)
  if( IsJetInDcalOnly(eta, phi, 0) ) {
    jetAcceptanceType |= AliEmcalJet::kDCALonly;
    // Check if DCALonlyfid
    if( IsJetInDcalOnly(eta, phi, r) )
      jetAcceptanceType |= AliEmcalJet::kDCALonlyfid;
  }
  
  // Check if PHOS
  if( IsJetInPhos(eta, phi, 0) ) {
    jetAcceptanceType |= AliEmcalJet::kPHOS;
    // Check if PHOSfid
    if( IsJetInPhos(eta, phi, r) )
      jetAcceptanceType |= AliEmcalJet::kPHOSfid;
  }
 
  return jetAcceptanceType;
}

/**
 * Returns whether or not jet with given eta, phi, R is in EMCal.
 */
Bool_t AliEmcalJetTask::IsJetInEmcal(Double_t eta, Double_t phi, Double_t r)
{
  if (!fGeom) return kFALSE;

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

  return kFALSE;
}

/**
 * Returns whether or not jet with given eta, phi, R is in DCal region (note: spans most of PHOS as well).
 */
Bool_t AliEmcalJetTask::IsJetInDcal(Double_t eta, Double_t phi, Double_t r)
{
  if (!fGeom) return kFALSE;
  if (eta < fGeom->GetArm1EtaMax() - r && eta > fGeom->GetArm1EtaMin() + r ) {
    if ( phi < fGeom->GetDCALPhiMax() * TMath::DegToRad() - r && phi > fGeom->GetDCALPhiMin() * TMath::DegToRad() + r)
      return kTRUE;
  }
  return kFALSE;
}

/**
 * Returns whether or not jet with given eta, phi, R is in DCal (note: ONLY DCal -- none of PHOS included).
 * Assumes DCAL_8SM geometry.
 * For r=0, use the entire DCal acceptance, including both of the connecting 1/3 SMs.
 * For r>0, use only the two "disjoint" fiducial regions of the DCal (i.e. ignore the connecting portions of the 1/3 SMs)
 */
Bool_t AliEmcalJetTask::IsJetInDcalOnly(Double_t eta, Double_t phi, Double_t r)
{
  if (!fGeom) return kFALSE;
  
  if (eta < fGeom->GetArm1EtaMax() - r && eta > fGeom->GetArm1EtaMin() + r) {
    if ( phi < fGeom->GetDCALPhiMax() * TMath::DegToRad() - r && phi > fGeom->GetDCALPhiMin() * TMath::DegToRad() + r) {
      
      if (TMath::Abs(eta) > fGeom->GetDCALInnerExtandedEta() + r) {
        return kTRUE;
      }
      if (r < 1e-6) {
        if (phi > fGeom->GetEMCGeometry()->GetDCALStandardPhiMax() * TMath::DegToRad())
          return kTRUE;
      }
      
    }
  }
  
  return kFALSE;
}

/**
 * Returns whether or not jet with given eta, phi, R is in PHOS.
 */
Bool_t AliEmcalJetTask::IsJetInPhos(Double_t eta, Double_t phi, Double_t r)
{
  Double_t etaMax = 0.130;
  Double_t etaMin = -0.130;
  Double_t phiMax = 320;
  Double_t phiMin = 260; // Run 1
  if (fRunNumber > 209121)
    phiMin = 250; // Run 2
  
  if (eta < etaMax - r && eta > etaMin + r ) {
    if (phi < phiMax * TMath::DegToRad() - r && phi > phiMin * TMath::DegToRad() + r)
      return kTRUE;
  }
  return kFALSE;
}

/**
 * Add an instance of this class to the analysis manager
 * @param nTracks name of the track collection
 * @param nClusters name of the calorimeter cluster collection
 * @param jetAlgo jet finding algorithm (anti-kt, kt, etc.)
 * @param radius jet resolution parameter (0.2, 0.4, tyc.)
 * @param jetType full, charged or neutral
 * @param minTrPt cut on the minimum transverse momentum of tracks
 * @param minClPt cut on the minimum transverse momentum of calorimeter clusters
 * @param ghostArea area of ghost particles (determines the jet area resolution)
 * @param reco recombination scheme
 * @param tag addtional information to be appended at the end of the output jet collection name
 * @param minJetPt cut on the minimum jet pt
 * @param lockTask lock the task - no further changes are possible if kTRUE
 * @param bFillGhosts add ghosts particles among the jet constituents in the output
 * @return a pointer to the new AliEmcalJetTask instance
 */
AliEmcalJetTask* AliEmcalJetTask::AddTaskEmcalJet(
  const TString nTracks, const TString nClusters,
  const AliJetContainer::EJetAlgo_t jetAlgo, const Double_t radius, const AliJetContainer::EJetType_t jetType,
  const Double_t minTrPt, const Double_t minClPt,
  const Double_t ghostArea, const AliJetContainer::ERecoScheme_t reco,
  const TString tag, const Double_t minJetPt,
  const Bool_t lockTask, const Bool_t bFillGhosts
)
{
  // Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEmcalJet", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AddTaskEmcalJet", "This task requires an input event handler");
    return 0;
  }

  EDataType_t dataType = kUnknownDataType;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(nTracks);
  TString clusName(nClusters);

  if (trackName == "usedefault") {
    if (dataType == kESD) {
      trackName = "Tracks";
    }
    else if (dataType == kAOD) {
      trackName = "tracks";
    }
    else {
      trackName = "";
    }
  }

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    }
    else if (dataType == kAOD) {
      clusName = "caloClusters";
    }
    else {
      clusName = "";
    }
  }

  AliParticleContainer* partCont = 0;
  if (trackName == "mcparticles") {
    AliMCParticleContainer* mcpartCont = new AliMCParticleContainer(trackName);
    partCont = mcpartCont;
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    AliTrackContainer* trackCont = new AliTrackContainer(trackName);
    partCont = trackCont;
  }
  else if (!trackName.IsNull()) {
    partCont = new AliParticleContainer(trackName);
  }
  if (partCont) partCont->SetParticlePtCut(minTrPt);

  AliClusterContainer* clusCont = 0;
  if (!clusName.IsNull()) {
    clusCont = new AliClusterContainer(clusName);
    clusCont->SetClusECut(0.);
    clusCont->SetClusPtCut(0.);
    clusCont->SetClusHadCorrEnergyCut(minClPt);
    clusCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  switch (jetType) {
  case AliJetContainer::kChargedJet:
    if (partCont) partCont->SetCharge(AliParticleContainer::kCharged);
    break;
  case AliJetContainer::kNeutralJet:
    if (partCont) partCont->SetCharge(AliParticleContainer::kNeutral);
    break;
  default:
    break;
  }

  TString name = AliJetContainer::GenerateJetName(jetType, jetAlgo, reco, radius, partCont, clusCont, tag);

  Printf("Jet task name: %s", name.Data());

  AliEmcalJetTask* mgrTask = static_cast<AliEmcalJetTask *>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;

  AliEmcalJetTask* jetTask = new AliEmcalJetTask(name);
  jetTask->SetJetType(jetType);
  jetTask->SetJetAlgo(jetAlgo);
  jetTask->SetRecombScheme(reco);
  jetTask->SetRadius(radius);
  if (partCont) jetTask->AdoptParticleContainer(partCont);
  if (clusCont) jetTask->AdoptClusterContainer(clusCont);
  jetTask->SetJetsName(tag);
  jetTask->SetMinJetPt(minJetPt);
  jetTask->SetGhostArea(ghostArea);

  if (bFillGhosts) jetTask->SetFillGhost();
  if (lockTask) jetTask->SetLocked();

  // Final settings, pass to manager and set the containers

  mgr->AddTask(jetTask);

  // Create containers for input/output
  AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(jetTask, 0, cinput);

  TObjArray* cnt = mgr->GetContainers();

  return jetTask;
}
