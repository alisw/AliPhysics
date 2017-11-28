/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <vector>

#include <TClonesArray.h>
#include <TObjArray.h>
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
#include "AliEmcalClusterJetConstituent.h"
#include "AliEmcalParticleJetConstituent.h"

#include "AliEmcalJetFinderKernel.h"

using std::cout;
using std::endl;
using std::cerr;

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::AliEmcalJetFinderKernel);
/// \endcond

namespace PWGJE{

namespace EMCALJetTasks {

const Int_t AliEmcalJetFinderKernel::fgkConstIndexShift = 100000;

/**
 * Default constructor. This constructor is only for ROOT I/O and
 * not to be used by users.
 */
AliEmcalJetFinderKernel::AliEmcalJetFinderKernel() :
  TNamed(),
  fGeom(nullptr),
  fParticleCollArray(),
  fClusterCollArray(),
  fJetType(AliJetContainer::kFullJet),
  fJetAlgo(AliJetContainer::antikt_algorithm),
  fRecombScheme(AliJetContainer::pt_scheme),
  fRadius(0.4),
  fRunNumber(-1),
  fMinMCLabel(0),
  fMinJetArea(0.001),
  fMinJetPt(1.0),
  fJetPhiMin(-10),
  fJetPhiMax(+10),
  fJetEtaMin(-1),
  fJetEtaMax(+1),
  fGhostArea(0.005),
  fTrackEfficiency(1.),
  fUtilities(0),
  fTrackEfficiencyOnlyForEmbedding(kFALSE),
  fLocked(0),
  fFillConstituents(kTRUE),
  fIndexMapsPrepared(kFALSE),
  fIsInit(0),
  fIsEmcPart(0),
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fJets(0),
  fFastJetWrapper("AliEmcalJetFinderKernel","AliEmcalJetFinderKernel"),
  fClusterContainerIndexMap(),
  fParticleContainerIndexMap()
{
}

/**
 * Standard named constructor.
 * @param name Name of the task.
 */
AliEmcalJetFinderKernel::AliEmcalJetFinderKernel(const char *name) :
  TNamed(name,""),
  fGeom(nullptr),
  fParticleCollArray(),
  fClusterCollArray(),
  fJetType(AliJetContainer::kFullJet),
  fJetAlgo(AliJetContainer::antikt_algorithm),
  fRecombScheme(AliJetContainer::pt_scheme),
  fRadius(0.4),
  fRunNumber(-1),
  fMinMCLabel(0),
  fMinJetArea(0.001),
  fMinJetPt(1.0),
  fJetPhiMin(-10),
  fJetPhiMax(+10),
  fJetEtaMin(-1),
  fJetEtaMax(+1),
  fGhostArea(0.005),
  fTrackEfficiency(1.),
  fUtilities(0),
  fTrackEfficiencyOnlyForEmbedding(kFALSE),
  fLocked(0),
  fFillConstituents(kTRUE),
  fIndexMapsPrepared(kFALSE),
  fIsInit(0),
  fIsEmcPart(0),
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fJets(0),
  fFastJetWrapper(name,name),
  fClusterContainerIndexMap(),
  fParticleContainerIndexMap()
{
}

/**
 * Destructor
 */
AliEmcalJetFinderKernel::~AliEmcalJetFinderKernel()
{
  if(fUtilities) delete fUtilities;
  if(fJets) delete fJets;
}

/**
 * Add a utility to the utility list. Utilities are instances of classes
 * derived from AliEmcalJetUtility that implements wrappers to FastJet contribs.
 */
AliEmcalJetUtility* AliEmcalJetFinderKernel::AddUtility(AliEmcalJetUtility* utility)
{
  if (!fUtilities) fUtilities = new TObjArray();
  fUtilities->SetOwner();
  if (fUtilities->FindObject(utility)) {
    Error("AddUtility", "Jet utility %s already connected.", utility->GetName());
    return utility;
  }   
  fUtilities->Add(utility);

  return utility;
}

/**
 * This method is called once before analyzing the first event. It executes
 * the Init() method of all utilities (if any).
 */
void AliEmcalJetFinderKernel::InitUtilities()
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->Init();
}
/**
 * This method is called before analyzing each event. It executes
 * the InitEvent() method of all utilities (if any).
 */
void AliEmcalJetFinderKernel::InitEvent()
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
void AliEmcalJetFinderKernel::PrepareUtilities()
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->Prepare(fFastJetWrapper);
}

/**
 * This method is called in the event loop for each jet found, while filling the output jet branch.
 * It executes the ProcessJet() method of all utilities (if any).
 */
void AliEmcalJetFinderKernel::ExecuteUtilities(AliEmcalJet* jet, Int_t ij)
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->ProcessJet(jet, ij, fFastJetWrapper);
}

/**
 * This method is called in the event loop after jet finding has been completed.
 * It executes the Terminate() method of all utilities (if any).
 */
void AliEmcalJetFinderKernel::TerminateUtilities()
{
  TIter next(fUtilities);
  AliEmcalJetUtility *utility = 0;
  while ((utility=static_cast<AliEmcalJetUtility*>(next()))) utility->Terminate(fFastJetWrapper);
}

/**
 * This method is called for each event.
 * @return Always kTRUE
 */
void AliEmcalJetFinderKernel::RunJetFinder(TClonesArray *jetarray)
{
  InitEvent();
  // clear the jet array (normally a null operation)
  fJets->Delete();
  Int_t n = FindJets();

  if(n){
    TClonesArray *outputArray = jetarray;
    if(!outputArray){
      if(!fJets) {
        fJets = new TClonesArray("AliEmcalJet");
      }
      outputArray = fJets;
    }
    FillJetBranch(outputArray);
  }
}

/**
 * This method steers the jet finding. It first loops over all particle and cluster containers
 * that were provided when the task was initialized. All accepted objects (tracks, particle, clusters)
 * are added as input vectors to the FastJet wrapper. Then the jet finding is launched
 * in the wrapper.
 * @return Total number of jets found.
 */
Int_t AliEmcalJetFinderKernel::FindJets()
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
void AliEmcalJetFinderKernel::FillJetBranch(TClonesArray *outputcont)
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

    AliEmcalJet *jet = new ((*outputcont)[jetCount])
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
Bool_t AliEmcalJetFinderKernel::GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const
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
void AliEmcalJetFinderKernel::Init()
{
  if (fTrackEfficiency < 1.) {
    if (gRandom) delete gRandom;
    gRandom = new TRandom3(0);
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

}


void AliEmcalJetFinderKernel::InitIndexMaps(){
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
void AliEmcalJetFinderKernel::FillJetConstituents(AliEmcalJet *jet, std::vector<fastjet::PseudoJet>& constituents,
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
      if(fFillConstituents){
        jet->AddParticleConstituent(t, partCont->GetIsEmbedding(), fParticleContainerIndexMap.GlobalIndexFromLocalIndex(partCont, tid));
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
      Double_t pvec[3] = {nP.Px(), nP.Py(), nP.Pz()};
      if(fFillConstituents) jet->AddClusterConstituent(c, (AliVCluster::VCluUserDefEnergy_t)clusCont->GetDefaultClusterEnergy(), pvec, clusCont->GetIsEmbedding(), fClusterContainerIndexMap.GlobalIndexFromLocalIndex(clusCont, cid));

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
Bool_t AliEmcalJetFinderKernel::IsLocked() const
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
 * Set the eta range of the track constituents.
 * @param emi Minimum eta
 * @param ema Maximum eta
 */
void AliEmcalJetFinderKernel::SetEtaRange(Double_t emi, Double_t ema)
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
void AliEmcalJetFinderKernel::SetMinJetClusPt(Double_t min)
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
void AliEmcalJetFinderKernel::SetMinJetClusE(Double_t min)
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
void AliEmcalJetFinderKernel::SetMinJetTrackPt(Double_t min)
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
void AliEmcalJetFinderKernel::SetPhiRange(Double_t pmi, Double_t pma)
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
fastjet::JetAlgorithm AliEmcalJetFinderKernel::ConvertToFJAlgo(EJetAlgo_t algo)
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
    ::Error("AliEmcalJetFinderKernel::ConvertToFJAlgo", "Jet algorithm %d not recognized!!!", algo);
    return fastjet::undefined_jet_algorithm;
  }
}

/**
 * Converts the internal enum values representing jet recombination schemes in
 * the corresponding values accepted by the FastJet wrapper.
 * @param reco Recombination scheme represented in the EJetAlgo_t enum
 * @return Recombination scheme represented in the fastjet::JetAlgorithm enum
 */
fastjet::RecombinationScheme AliEmcalJetFinderKernel::ConvertToFJRecoScheme(ERecoScheme_t reco)
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
    ::Error("AliEmcalJetFinderKernel::ConvertToFJRecoScheme", "Recombination scheme %d not recognized!!!", reco);
    return fastjet::external_scheme;
  }
}

/**
 * Finds which geometrical acceptance types the jet satisfies.
 * @return bitwise jet acceptance type
 */
UInt_t AliEmcalJetFinderKernel::FindJetAcceptanceType(Double_t eta, Double_t phi, Double_t r) {
  
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
Bool_t AliEmcalJetFinderKernel::IsJetInEmcal(Double_t eta, Double_t phi, Double_t r)
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
Bool_t AliEmcalJetFinderKernel::IsJetInDcal(Double_t eta, Double_t phi, Double_t r)
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
Bool_t AliEmcalJetFinderKernel::IsJetInDcalOnly(Double_t eta, Double_t phi, Double_t r)
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
Bool_t AliEmcalJetFinderKernel::IsJetInPhos(Double_t eta, Double_t phi, Double_t r)
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

AliParticleContainer *AliEmcalJetFinderKernel::GetParticleContainer(Int_t ncontainer){
  if(ncontainer >= fParticleCollArray.GetEntriesFast()) return nullptr;
  return static_cast<AliParticleContainer *>(fParticleCollArray.At(ncontainer));
}

AliParticleContainer *AliEmcalJetFinderKernel::GetParticleContainer(const char * ncontainer){
  return static_cast<AliParticleContainer *>(fParticleCollArray.FindObject(ncontainer));
}


AliClusterContainer *AliEmcalJetFinderKernel::GetClusterContainer(Int_t ncontainer){
  if(ncontainer >= fParticleCollArray.GetEntriesFast()) return nullptr;
  return static_cast<AliClusterContainer *>(fParticleCollArray.At(ncontainer));
}

AliClusterContainer *AliEmcalJetFinderKernel::GetClusterContainer(const char * ncontainer){
  return static_cast<AliClusterContainer *>(fParticleCollArray.FindObject(ncontainer));
}

void AliEmcalJetFinderKernel::AppendParticleContainer(AliParticleContainer *cont) {
  fParticleCollArray.AddLast(cont);
}

void AliEmcalJetFinderKernel::AppendClusterContainer(AliClusterContainer *cont) {
  fClusterCollArray.AddLast(cont);
}

void AliEmcalJetFinderKernel::SetParticleContainer(AliParticleContainer *cont, UInt_t position) {
  fParticleCollArray.AddAt(cont, position);
}

void AliEmcalJetFinderKernel::SetClusterContainer(AliClusterContainer *cont, UInt_t position) {
  fClusterCollArray.AddAt(cont, position);
}

}

}
