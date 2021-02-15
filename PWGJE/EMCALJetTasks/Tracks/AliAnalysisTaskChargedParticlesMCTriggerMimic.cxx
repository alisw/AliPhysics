/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <array>
#include <bitset>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include <TClonesArray.h>
#include <TFile.h>
#include <THistManager.h>
#include <TLinearBinning.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TTree.h>

#include "AliAnalysisUtils.h"
#include "AliAODInputHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliClusterContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCalTriggerWeightHandler.h"
#include "AliESDtrack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVParticle.h"

#include "AliAnalysisTaskChargedParticlesMCTriggerMimic.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskChargedParticlesMCTriggerMimic)

using namespace PWGJE::EMCALJetTasks;

/**
 * Dummy constructor
 */
AliAnalysisTaskChargedParticlesMCTriggerMimic::AliAnalysisTaskChargedParticlesMCTriggerMimic():
            AliAnalysisTaskEmcal(),
            fTrackCuts(NULL),
            fHistos(NULL),
            fWeightHandler(NULL),
            fYshift(0.465),
            fEtaSign(1),
            fEtaLabCut(-0.5, 0.5),
            fEtaCmsCut(-0.13, 0.13),
            fPhiCut(0, TMath::TwoPi()),
            fPatchType(kUndef),
            fEnergyThreshold(0.),
            fObservables(),
            fNameClusters()
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetNeedEmcalGeom(true);
}

/**
 * Main constructor
 * @param name Name of the task
 */
AliAnalysisTaskChargedParticlesMCTriggerMimic::AliAnalysisTaskChargedParticlesMCTriggerMimic(const char *name):
            AliAnalysisTaskEmcal(name, true),
            fTrackCuts(nullptr),
            fHistos(nullptr),
            fWeightHandler(nullptr),
            fYshift(0.465),
            fEtaSign(1),
            fEtaLabCut(-0.6, 0.6),
            fEtaCmsCut(-0.13, 0.13),
            fPatchType(kUndef),
            fEnergyThreshold(0.),
            fObservables(),
            fNameClusters()
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetNeedEmcalGeom(true);
  SetMakeGeneralHistograms(true);
}

/**
 * Destructor - cleaning up
 */
AliAnalysisTaskChargedParticlesMCTriggerMimic::~AliAnalysisTaskChargedParticlesMCTriggerMimic() {
  //if(fTrackCuts) delete fTrackCuts;
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskChargedParticlesMCTriggerMimic::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if(fIsPythia){
    AliDebugStream(1) << GetName() << ": Running on PYTHIA Hard production" << std::endl;
  } else {
    AliDebugStream(1) << GetName() << ": Not running on PYTHIA Hard production" << std::endl;
  }

  if(!fAliAnalysisUtils) fAliAnalysisUtils = new AliAnalysisUtils;
  if(!fTrackCuts) InitializeTrackCuts("standard",fInputHandler->IsA() == AliAODInputHandler::Class());

  PtBinning newbinning;
  TLinearBinning smbinning(21, -0.5, 20.5), etabinning(100, -0.7, 0.7);
  fHistos = new THistManager("Ref");

  // The first two histograms should lead to a re-implementation inside AliAnalysisTaskEMCAL:
  // Users should be able to set event weights.
  fHistos->CreateTH1("hUserEventCount", "User event counter", 1, 0.5, 1.5);
  fHistos->CreateTH1("hUserVertexZ", "User vertex distribution after z-cut", 100, -10, 10);
  fHistos->CreateTH1("hUserPtHard", "User pt-hard distribution", 1000, 0., 300.);

  // Histograms for observable tracks
  if(HasObservable(kTracks)){
    const std::array<std::string, 2> kInputs = {"True", "Accept"};
    const std::array<std::string, 6> kSpecies = {"El", "Mu", "Pi", "Ka", "Pr", "Ot"};
    const std::array<double, 5> kPtCuts = {1., 2., 5., 10., 20.};
    for(const auto &input : kInputs){
      AliDebugStream(1) << GetName() << ": Creating histograms for case " << input << std::endl;
      fHistos->CreateTH1(Form("hTrackPtEtaAll%s", input.c_str()), Form("Charged particle p_{t} distribution all #eta %s", input.c_str()), newbinning, "s");
      fHistos->CreateTH1(Form("hTrackPtEtaCent%s", input.c_str()), Form("Charged particle p_{t} distribution central #eta  %s", input.c_str()), newbinning, "s");
      fHistos->CreateTH1(Form("hTrackPtEMCALEtaAll%s", input.c_str()), Form("Charged particle in EMCAL p_{t} distribution all #eta trigger %s", input.c_str()), newbinning);
      fHistos->CreateTH1(Form("hTrackPtEMCALEtaCent%s", input.c_str()), Form("Charged particle in EMCAL p_{t} distribution central eta trigger %s", input.c_str()), newbinning);
      for(const auto &piditer : kSpecies){
        fHistos->CreateTH1(Form("hTrackPtEtaAll%s%s", piditer.c_str(), input.c_str()), Form("Charged %s p_{t} distribution all #eta %s", piditer.c_str(), input.c_str()), newbinning);
        fHistos->CreateTH1(Form("hTrackPtEtaCent%s%s", piditer.c_str(), input.c_str()), Form("Charged %s p_{t} distribution central #eta %s", piditer.c_str(), input.c_str()), newbinning);
        fHistos->CreateTH1(Form("hTrackPtEMCALEtaAll%s%s", piditer.c_str(), input.c_str()), Form("Charged %s in EMCAL p_{t} distribution all #eta %s", piditer.c_str(), input.c_str()), newbinning);
        fHistos->CreateTH1(Form("hTrackPtEMCALEtaCent%s%s", piditer.c_str(), input.c_str()), Form("Charged %s in EMCAL p_{t} distribution central #eta %s", piditer.c_str(), input.c_str()), newbinning);
      }
      for(const auto &ptcut : kPtCuts){
        fHistos->CreateTH1(
            Form("hTrackEtaLabDistAllPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
            Form("#eta_{lab} distribution without #eta-cut for tracks with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
            100,
            -1.,
            1.
        );
        fHistos->CreateTH1(
            Form("hTrackEtaLabDistCutPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
            Form("#eta_{lab} distribution with #eta-cut for tracks with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
            100,
            -1.,
            1.
        );
        fHistos->CreateTH1(
            Form("hTrackEtaCentDistAllPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
            Form("#eta_{cent} distribution without #eta-cut for tracks with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
            160,
            -1.3,
            1.3
        );
        fHistos->CreateTH1(
            Form("hTrackEtaCentDistCutPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
            Form("#eta_{cent} distribution with #eta-cut for tracks with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
            160,
            -1.3,
            1.3
        );
        fHistos->CreateTH1(
            Form("hTrackEtaLabDistAllEMCALPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
            Form("#eta_{lab} distribution without #eta-cut for tracks in EMCAL with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
            100,
            -1.,
            1.
        );
        fHistos->CreateTH1(
            Form("hTrackEtaLabDistCutEMCALPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
            Form("#eta_{lab} distribution with #eta-cut for tracks in EMCAL with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
            100,
            -1.,
            1.
        );
        fHistos->CreateTH1(
            Form("hTrackEtaCentDistAllEMCALPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
            Form("#eta_{cent} distribution without #eta-cut for tracks in EMCAL with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
            160,
            -1.3,
            1.3
        );
        fHistos->CreateTH1(
            Form("hTrackEtaCentDistCutEMCALPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
            Form("Eta (cent) distribution with #eta-cut for tracks in EMCAL with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
            160,
            -1.3,
            1.3
        );
        fHistos->CreateTH1(
            Form("hTrackPhiDistAllPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
            Form("#phi distribution of particles with p_{t} above %.1f GeV/c trigger %s", ptcut, input.c_str()),
            300,
            0.,
            2*TMath::Pi()
        );
      }
    }
  }

  // Histograms for observable clusters
  std::array<Double_t, 5> kEnCuts = {1., 2., 5., 10., 20.};
  if(HasObservable(kClusters)){
    Int_t sectorsWithEMCAL[10] = {4, 5, 6, 7, 8, 9, 13, 14, 15, 16};
    fHistos->CreateTH1("hClusterEnergy", "Cluster energy", newbinning);
    fHistos->CreateTH1("hClusterET", "Cluster transverse energy", newbinning);
    fHistos->CreateTH2("hClusterEnergySM", "Cluster energy versus supermodule", smbinning, newbinning);
    fHistos->CreateTH2("hClusterETSM", "Cluster transverse energy versus supermodule", smbinning, newbinning);
    fHistos->CreateTH2("hClusterEtaEnergy", "Cluster energy vs. eta", etabinning, newbinning);
    fHistos->CreateTH2("hClusterEtaET", "Cluster transverse energy vs. eta", etabinning, newbinning);
    for(int ism = 0; ism < 20; ism++){
      fHistos->CreateTH2(Form("hClusterEtaEnergySM%d", ism), Form("Cluster energy vs. eta in Supermodule %d", ism), etabinning, newbinning);
      fHistos->CreateTH2(Form("hClusterEtaETSM%d", ism), Form("Cluster transverse energy vs. eta in Supermodule %d", ism), etabinning, newbinning);
    }
    for(int isec = 0; isec < 10; isec++){
      fHistos->CreateTH2(Form("hClusterEtaEnergySec%d", sectorsWithEMCAL[isec]), Form("Cluster energy vs.eta in tracking sector %d", sectorsWithEMCAL[isec]), etabinning, newbinning);
      fHistos->CreateTH2(Form("hClusterEtaETSec%d", sectorsWithEMCAL[isec]), Form("Cluster transverse energy vs.eta in tracking sector %d", sectorsWithEMCAL[isec]), etabinning, newbinning);
    }
    for(auto ien : kEnCuts){
      fHistos->CreateTH2(Form("hClusterEtaPhi%dG", static_cast<int>(ien)), Form("cluster #eta-#phi map for clusters with energy larger than %f GeV/c", ien), 100, -0.7, 0.7, 200, 0, 2*TMath::Pi());
    }
  }

  // Histograms for Observable EGA or EJE patches
  if(HasObservable(kEGAPatches) || HasObservable(kEJEPatches)){
    const int kNPatchTypes = 2;
    const std::array<Observable_t, kNPatchTypes> kPatchObservables = {kEGAPatches, kEJEPatches};
    const std::array<TString, kNPatchTypes> kPatchTypes = {"EGA","EJE"};
    for(int ipatch = 0; ipatch < kNPatchTypes; ipatch++){
      if(!HasObservable(kPatchObservables[ipatch])) continue;
      fHistos->CreateTH1(Form("h%sPatchEnergy", kPatchTypes[ipatch].Data()), Form("%s-patch energy", kPatchTypes[ipatch].Data()), newbinning);
      fHistos->CreateTH1(Form("h%sPatchET", kPatchTypes[ipatch].Data()), Form("%s-patch transverse energy", kPatchTypes[ipatch].Data()), newbinning);
      fHistos->CreateTH2(Form("h%sPatchEnergyEta", kPatchTypes[ipatch].Data()), Form("%s-patch energy", kPatchTypes[ipatch].Data()), newbinning, etabinning);
      fHistos->CreateTH2(Form("h%sPatchETEta", kPatchTypes[ipatch].Data()), Form("%s-patch transverse energy", kPatchTypes[ipatch].Data()), newbinning, etabinning);
      for(auto enc : kEnCuts){
        fHistos->CreateTH2(Form("h%sPatchEtaPhi%dG", kPatchTypes[ipatch].Data(), static_cast<int>(enc)), Form("%s-patch #eta-#phi map for patches with energy larger than %f GeV/c", kPatchTypes[ipatch].Data(), enc), 100, -0.7, 0.7, 200, 0, TMath::TwoPi());
        fHistos->CreateTH2(Form("h%sPatchColRow%dG", kPatchTypes[ipatch].Data(), static_cast<int>(enc)), Form("%s-patch col-row map for patches with energy larger than %f GeV/c", kPatchTypes[ipatch].Data(), enc), 48, -0.5, 47.5, 104, -0.5, 103.5);
        }
    }
  }

  for(auto hist : *(fHistos->GetListOfHistograms())) fOutput->Add(hist);
  PostData(1, fOutput);

  AliDebugStream(1) << GetName() << ": Output objects initialized" << std::endl;
}

/**
 * Perform event selection: This function overwrites the
 * event selection of the AliAnalysisTaskEmcal
 * @return Result of the event selection (true if event is selected)
 */
Bool_t AliAnalysisTaskChargedParticlesMCTriggerMimic::IsEventSelected(){
  AliDebugStream(2) << GetName() << ": Using custom event selection method" << std::endl;
  if(!fTriggerPatchInfo){
    AliErrorStream() << GetName() << ": Trigger patch container not found but required" << std::endl;
    return false;
  }
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return false;
  AliDebugStream(3) << GetName() << "Event is an INT7 event" << std::endl;

  // MC outlier cut
  if(fIsPythia){
    if(!CheckMCOutliers()){
      AliDebugStream(3) << GetName() << ": Event identified as outlier" << std::endl;
      return false;
    } else {
      AliDebugStream(3) << GetName() << ": Not an outlier event" << std::endl;
    }
  }

  // Generall event quality cuts
  // The vertex cut also contains cuts on the number
  // of contributors and the position in z
  AliDebugStream(3) << GetName() << ": Applying vertex selection" << std::endl;
  if(fAliAnalysisUtils){
    if(!fAliAnalysisUtils->IsVertexSelected2013pA(InputEvent())) return false;
    if(fAliAnalysisUtils->IsPileUpEvent(InputEvent())) return false;
    AliDebugStream(3) << GetName() << ": Vertex selection passed" << std::endl;
  }

  AliDebugStream(3) << GetName() << ": Applying EMCAL trigger selection" << std::endl;
  if(fPatchType != kUndef){
    if(!SelectEmcalTrigger(fTriggerPatchInfo)){
      AliDebugStream(3) << GetName() << ": Failed trigger selection" << std::endl;
      return false;
    }
  }

  AliDebugStream(2) << GetName() << "Event selected" << std::endl;
  return true;
}

/**
 * Run function
 * @return Always true
 */
Bool_t AliAnalysisTaskChargedParticlesMCTriggerMimic::Run(){
  AliDebugStream(1) << GetName() << ": Inspecting event" << std::endl;
  Double_t weight = 1.;
  if(fIsPythia){
    if(!fPythiaHeader){
      AliErrorStream() << GetName() << ": PYTHIA event header not found" << std::endl;
    } else {
      weight = fWeightHandler ? fWeightHandler->GetEventWeight(fPythiaHeader) : 1.;
    }
  }

  fHistos->FillTH1("hUserEventCount", 1, weight);
  const AliVVertex *vtx = InputEvent()->GetPrimaryVertex();
  fHistos->FillTH1("hUserVertexZ", vtx->GetZ(), weight);
  fHistos->FillTH1("hUserPtHard", fPtHard);

  if(HasObservable(kTracks)){
    AliDebugStream(3) << GetName() << ": eta-lab cut: " << fEtaLabCut << ", eta-cms cut: " << fEtaCmsCut << std::endl;

    AliVParticle *truepart(nullptr);
    Bool_t isEMCAL(kFALSE);
    if(MCEvent()){
      for(int ipart = 0; ipart < fMCEvent->GetNumberOfTracks(); ipart++){
        truepart = MCEvent()->GetTrack(ipart);

        // Select only particles within ALICE acceptance
        if(!fEtaLabCut.IsInRange(truepart->Eta())) continue;
        if(!fPhiCut.IsInRange(truepart->Phi())) continue;
        if(TMath::Abs(truepart->Pt()) < 0.1) continue;
        if(!truepart->Charge()) continue;

        if(!IsPhysicalPrimary(truepart, fMCEvent)) continue;
        isEMCAL = (truepart->Phi() > 1.5 && truepart->Phi() < 3.1) ? kTRUE : kFALSE;

        // Calculate eta in cms frame according
        // EPJC74 (2014) 3054:
        // eta_cms = - eta_lab - |yshift|
        Double_t etacent = -1. * truepart->Eta() - TMath::Abs(fYshift);
        etacent *= fEtaSign;

        Bool_t etacentcut = fEtaCmsCut.IsInRange(etacent);

        // Get PID
        TString pid = "";
        switch(TMath::Abs(truepart->PdgCode())){
        case kPiPlus: pid = "Pi"; break;
        case kMuonMinus: pid = "Mu"; break;
        case kElectron: pid = "El"; break;
        case kKPlus: pid = "Ka"; break;
        case kProton: pid = "Pr"; break;
        default: pid = "Ot"; break;
        };

        // Particle selected (do not filter TRD sectors for MC truth)
        FillTrackHistos("True", weight, TMath::Abs(truepart->Pt()), truepart->Eta() * fEtaSign, etacent, truepart->Phi(), etacentcut, isEMCAL, pid);
      }
    }

    // Loop over tracks, fill select particles
    // Histograms
    // - Full eta (-0.8, 0.8), new binning
    // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c without eta cut
    // - Central eta (-0.8, -0.2), new binning,
    // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c
    // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c with eta cut
    AliVTrack *checktrack(NULL);
    AliVParticle *assocMC(NULL);
    double ptparticle(-1.), etaparticle(-100.), etaEMCAL(0.), phiEMCAL(0.);
    for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); ++itrk){
      checktrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(itrk));
      if(!checktrack) continue;
      TString assocpid = "Ot";
      // Find associated particle
      if(MCEvent()){
        assocMC = MCEvent()->GetTrack(TMath::Abs(checktrack->GetLabel()));
        if(!assocMC) continue;        // Fake track
        if(!IsPhysicalPrimary(assocMC, fMCEvent)) continue;
        // Get PID
        switch(TMath::Abs(assocMC->PdgCode())){
        case kPiPlus: assocpid = "Pi"; break;
        case kMuonMinus: assocpid = "Mu"; break;
        case kElectron: assocpid = "El"; break;
        case kKPlus: assocpid = "Ka"; break;
        case kProton: assocpid = "Pr"; break;
        default: assocpid = "Ot"; break;
        };
      }

      // Select only particles within ALICE acceptance
      if(!fEtaLabCut.IsInRange(checktrack->Eta())) continue;
      if(!fPhiCut.IsInRange(checktrack->Phi())) continue;
      if(TMath::Abs(checktrack->Pt()) < 0.1) continue;
      if(checktrack->IsA() == AliESDtrack::Class()){
        AliESDtrack copytrack(*(static_cast<AliESDtrack *>(checktrack)));
        AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&copytrack);
        etaEMCAL = copytrack.GetTrackEtaOnEMCal();
        phiEMCAL = copytrack.GetTrackPhiOnEMCal();
      } else {
        AliAODTrack copytrack(*(static_cast<AliAODTrack *>(checktrack)));
        AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&copytrack);
        etaEMCAL = copytrack.GetTrackEtaOnEMCal();
        phiEMCAL = copytrack.GetTrackPhiOnEMCal();
      }
      Int_t supermoduleID = -1;
      isEMCAL = fGeom->SuperModuleNumberFromEtaPhi(etaEMCAL, phiEMCAL, supermoduleID);
      // Exclude supermodules 10 and 11 as they did not participate in the trigger
      isEMCAL = isEMCAL && supermoduleID < 10;

      if(!fTrackCuts->IsTrackAccepted(checktrack)) continue;

      // prefer true pt and eta, however in case of running on data take measured values
      ptparticle = assocMC ? TMath::Abs(assocMC->Pt()) : TMath::Abs(checktrack->Pt());
      etaparticle = assocMC ? assocMC->Eta() : checktrack->Eta();

      // Calculate eta in cms frame according
      // EPJC74 (2014) 3054:
      // eta_cms = - eta_lab - |yshift|
      Double_t etacent = -1. * checktrack->Eta() - TMath::Abs(fYshift);
      etacent *= fEtaSign;

      Bool_t etacentcut = fEtaCmsCut.IsInRange(etacent);

      FillTrackHistos("Accept", weight,  ptparticle, checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, assocpid);
    }
  }

  if(HasObservable(kClusters)){
    Double_t vertexpos[3];
    fInputEvent->GetPrimaryVertex()->GetXYZ(vertexpos);

    Double_t energy, et, eta, phi;
    AliClusterContainer *clustercont = this->GetClusterContainer(fNameClusters.Data());
    if(clustercont){
      for(auto clust : clustercont->all()){
        if(!clust->IsEMCAL()) continue;
        if(clust->GetIsExotic()) continue;

        TLorentzVector posvec;
        energy = clust->GetNonLinCorrEnergy();
        et = posvec.Et();
        clust->GetMomentum(posvec, vertexpos);
        eta = posvec.Eta();
        phi = posvec.Phi();
        FillClusterHistos(weight, energy, et, eta, phi);
      }
    }
  }

  if(HasObservable(kEGAPatches) || HasObservable(kEJEPatches)){
    AliEMCALTriggerPatchInfo *recpatch(nullptr);
    TString patchname;
    for(TIter patchiter = TIter(fTriggerPatchInfo).Begin(); patchiter != TIter::End(); ++patchiter){
      recpatch = static_cast<AliEMCALTriggerPatchInfo *>(*patchiter);
      if(!recpatch->IsOfflineSimple()) continue;
      if(recpatch->IsJetHighSimple() && HasObservable(kEJEPatches)){
        patchname = "EJE";
      } else if(recpatch->IsGammaHighSimple() && HasObservable(kEGAPatches)){
        patchname = "EGA";
      } else continue;

      FillPatchHistos(patchname.Data(), weight, recpatch->GetPatchE(), recpatch->GetPatchET(), recpatch->GetEtaGeo(), recpatch->GetPhiGeo(), recpatch->GetColStart(), recpatch->GetRowStart());
    }
  }

  return true;
}

/**
 * Fill track histograms
 * @param[in] eventclass Trigger class fired
 * @param[in] weight \f$ p_{t} \f$-hard dependent weight
 * @param[in] pt track \f$ p_{t} \f$
 * @param[in] etalab Track \f$ \eta \f$ in lab frame
 * @param[in] etacent Track \f$ \eta \f$ in cms frame
 * @param[in] phi Track \f$ \eta \f$ in lab frame
 * @param[in] etacut Track accepted by \f$ \eta \f$ cut
 * @param[in] inEmcal Track in EMCAL \f$ \phi \f$ acceptance
 */
void AliAnalysisTaskChargedParticlesMCTriggerMimic::FillTrackHistos(
    const char *eventclass,
    Double_t weight,
    Double_t pt,
    Double_t etalab,
    Double_t etacent,
    Double_t phi,
    Bool_t etacut,
    Bool_t inEmcal,
    const char *pid
    )
{
  fHistos->FillTH1(Form("hTrackPtEtaAll%s", eventclass), TMath::Abs(pt), weight);
  fHistos->FillTH1(Form("hTrackPtEtaAll%s%s", pid, eventclass), TMath::Abs(pt), weight);
  if(inEmcal){
    fHistos->FillTH1(Form("hTrackPtEMCALEtaAll%s", eventclass), TMath::Abs(pt), weight);
    fHistos->FillTH1(Form("hTrackPtEMCALEtaAll%s%s", pid, eventclass), TMath::Abs(pt), weight);
  }

  const std::array<int, 5> kPtMin = {1,2,5,10,20}; // for eta distributions
  for(const auto &ptmin : kPtMin){
    if(TMath::Abs(pt) > static_cast<double>(ptmin)){
      fHistos->FillTH1(Form("hTrackPhiDistAllPt%d%s", ptmin, eventclass), phi, weight);
      fHistos->FillTH1(Form("hTrackEtaLabDistAllPt%d%s", ptmin, eventclass), etalab, weight);
      fHistos->FillTH1(Form("hTrackEtaCentDistAllPt%d%s", ptmin, eventclass), etacent, weight);
      if(inEmcal){
        fHistos->FillTH1(Form("hTrackEtaLabDistAllEMCALPt%d%s", ptmin, eventclass), etalab, weight);
        fHistos->FillTH1(Form("hTrackEtaCentDistAllEMCALPt%d%s", ptmin, eventclass), etacent, weight);
      }
    }
  }

  if(etacut){
    fHistos->FillTH1(Form("hTrackPtEtaCent%s", eventclass), TMath::Abs(pt), weight);
    fHistos->FillTH1(Form("hTrackPtEtaCent%s%s", pid, eventclass), TMath::Abs(pt), weight);
    if(inEmcal){
      fHistos->FillTH1(Form("hTrackPtEMCALEtaCent%s", eventclass), TMath::Abs(pt), weight);
      fHistos->FillTH1(Form("hTrackPtEMCALEtaCent%s%s", pid, eventclass), TMath::Abs(pt), weight);
    }
    for(const auto &ptmin : kPtMin){
      if(TMath::Abs(pt) > static_cast<double>(ptmin)){
        fHistos->FillTH1(Form("hTrackEtaLabDistCutPt%d%s", ptmin, eventclass), etalab, weight);
        fHistos->FillTH1(Form("hTrackEtaCentDistCutPt%d%s", ptmin, eventclass), etacent, weight);
        if(inEmcal){
          fHistos->FillTH1(Form("hTrackEtaLabDistCutEMCALPt%d%s", ptmin, eventclass), etalab, weight);
          fHistos->FillTH1(Form("hTrackEtaCentDistCutEMCALPt%d%s", ptmin, eventclass), etacent, weight);
        }
      }
    }
  }
}

/**
 * Fill cluster histograms
 * @param[in] weight \f$ p_{t} \f$-hard dependent weight
 * @param[in] energy Cluster energy
 * @param[in] transverseenergy Cluster transverse energy
 * @param[in] eta \f$ \eta \f$ of the cluster with respect to the primary vertex
 * @param[in] phi \f$ \phi \f$ of the cluster with respect to the primary vertex
 */
void AliAnalysisTaskChargedParticlesMCTriggerMimic::FillClusterHistos(
    double weight,
    double energy,
    double transverseenergy,
    double eta,
    double phi
    )
{
  Int_t supermoduleID = -1, sector = -1;
  fGeom->SuperModuleNumberFromEtaPhi(eta, phi, supermoduleID);

  fHistos->FillTH1("hClusterEnergy", energy, weight);
  fHistos->FillTH1("hClusterET", transverseenergy, weight);
  fHistos->FillTH2("hClusterEtaEnergy", eta, energy, weight);
  fHistos->FillTH2("hClusterEtaET", eta, transverseenergy, weight);
  if(supermoduleID >= 0){
    fHistos->FillTH2("hClusterEnergySM", supermoduleID, energy, weight);
    fHistos->FillTH2("hClusterETSM", supermoduleID, transverseenergy, weight);
    fHistos->FillTH2(Form("hClusterEtaEnergySM%d", supermoduleID), eta, energy, weight);
    fHistos->FillTH2(Form("hClusterEtaETSM%d", supermoduleID), eta, transverseenergy, weight);
    if(supermoduleID < 12)
      sector = 4 + int(supermoduleID/2); // EMCAL
    else
      sector = 13 + int((supermoduleID-12)/2);  // DCAL
    fHistos->FillTH2(Form("hClusterEtaEnergySec%d", sector), eta, energy, weight);
    fHistos->FillTH2(Form("hClusterEtaETSec%d", sector), eta, transverseenergy, weight);
  }
  std::array<Double_t, 5> encuts = {1., 2., 5., 10., 20.};
  for(auto e : encuts){
    if(energy > e){
      fHistos->FillTH2(Form("hClusterEtaPhi%dG", static_cast<int>(e)), eta, phi, weight);
    }
  }
}

/**
 * Filling patch related histogram.
 * @param[in] patchname Name of the patchtype
 * @param[in] energy Calibrated energy of the patch
 * @param[in] eta Patch eta at the geometrical center
 * @param[in] phi Patch phi at the geometrical center
 */
void AliAnalysisTaskChargedParticlesMCTriggerMimic::FillPatchHistos(const char *patchname, double weight, double energy, double transverseenergy, double eta, double phi, int col, int row){
  fHistos->FillTH1(Form("h%sPatchEnergy", patchname), energy, weight);
  fHistos->FillTH1(Form("h%sPatchET", patchname), transverseenergy, weight);
  fHistos->FillTH2(Form("h%sPatchEnergyEta", patchname), energy, eta, weight);
  fHistos->FillTH2(Form("h%sPatchETEta", patchname), transverseenergy, eta, weight);
  const std::array<Double_t, 5> kEnCuts = {1., 2., 5., 10., 20.};
  for(auto e : kEnCuts){
    if(energy > e){
      fHistos->FillTH2(Form("h%sPatchEtaPhi%dG", patchname, static_cast<int>(e)), eta, phi, weight);
      fHistos->FillTH2(Form("h%sPatchColRow%dG", patchname, static_cast<int>(e)), col, row, weight);
    }
  }
}


/**
 * Check in a transparent way for ESDs and AODs whether the particle is physical primary or not
 * -# AOD: Information stored in the AliAODMCParticle
 * -# ESD: Information needs to be retrieved from the stack via the label of the MC particle
 * @param[in] part The particle to check
 * @param[in] mcevent The MC event containing the stack (ESD only)
 * @return True if particle is a physical primary particle, false otherwise
 */
Bool_t AliAnalysisTaskChargedParticlesMCTriggerMimic::IsPhysicalPrimary(const AliVParticle* const part, AliMCEvent* const mcevent) {
  Bool_t physprim = false;
  const AliAODMCParticle *aodmc = dynamic_cast<const AliAODMCParticle *>(part);
  if(aodmc){
    physprim = aodmc->IsPhysicalPrimary();
  } else {
    physprim = mcevent->IsPhysicalPrimary(part->GetLabel());
  }
  return physprim;
}

/**
 * Set the track selection
 * @param cutname Name of the track cuts
 * @param isAOD check whether we run on ESDs or AODs
 */
void AliAnalysisTaskChargedParticlesMCTriggerMimic::InitializeTrackCuts(TString cutname, bool isAOD){
  SetEmcalTrackSelection(AliEmcalAnalysisFactory::TrackCutsFactory(cutname, isAOD));
}

/**
 * Select EMCAL-triggered event based on the presence of a trigger patch above
 * (offline) energy threshold. The threshold is settable, as well as at the patch
 * type.
 * @param[in] triggerpatches Array of trigger patches used for the trigger decision
 * @return True if the event is selected as triggered, false otherwise
 */
bool AliAnalysisTaskChargedParticlesMCTriggerMimic::SelectEmcalTrigger(const TClonesArray *triggerpatches){
  bool selected = false;
  AliEMCALTriggerPatchInfo *patch(nullptr);
  AliDebugStream(2) << GetName() << ": Selecting EMCAL triggered event type (" << (fPatchType == kEMCEGA ? "EGA" : "EJE") << ") using patch energy above threshold" << std::endl;
  AliDebugStream(2) << GetName() << ": Energy threshold " << fEnergyThreshold << " GeV" << std::endl;
  AliDebugStream(2) << GetName() << ": Number of reconstructed patches " << triggerpatches->GetEntries() << std::endl;
  for(TIter patchiter = TIter(triggerpatches).Begin(); patchiter != TIter::End(); ++patchiter){
    patch = static_cast<AliEMCALTriggerPatchInfo *>(*patchiter);
    AliDebugStream(4) << GetName() << ": Next patch" << std::endl;
    if(!patch->IsOfflineSimple()) continue;
    AliDebugStream(4) << GetName() << "Patch is an offline simple patch" << std::endl;
    AliDebugStream(4) << GetName() << ": Trigger bits: " << std::bitset<32>(patch->GetTriggerBits()) << std::endl;
    AliDebugStream(4) << GetName() << ": J1(" << patch->IsJetHighSimple()  << "), J2(" << patch->IsJetLowSimple()
                                   << "), G1(" << patch->IsGammaHighSimple() << ") G2(" << patch->IsGammaLowSimple() << ")" << std::endl;
    if(fPatchType == kEMCEJE){
      if(!patch->IsJetHighSimple()) continue;
      AliDebugStream(3) << GetName() << ": Patch is jet high simple" << std::endl;
    } else if(fPatchType == kEMCEGA){
      if(!patch->IsGammaHighSimple()) continue;
      AliDebugStream(4) << GetName() << ": Patch is gamma high simple" << std::endl;
    }
    AliDebugStream(3) << GetName() << ": Found trigger patch of matching type, now cutting on energy ...." << std::endl;
    if(patch->GetPatchE() > fEnergyThreshold){
      // firing trigger patch found
      AliDebugStream(2) << GetName() << ": Firing trigger patch found at energy " << std::setprecision(1) << patch->GetPatchE() << std::endl;
      selected = true;
      break;
    }
  }
  return selected;
}

/**
 * Create \f$ p_{t} \f$ binning
 */
AliAnalysisTaskChargedParticlesMCTriggerMimic::PtBinning::PtBinning() :
    TCustomBinning()
{
  this->SetMinimum(0.);
  this->AddStep(1, 0.05);
  this->AddStep(2, 0.1);
  this->AddStep(4, 0.2);
  this->AddStep(7, 0.5);
  this->AddStep(16, 1);
  this->AddStep(36, 2);
  this->AddStep(40, 4);
  this->AddStep(50, 5);
  this->AddStep(100, 10);
  this->AddStep(200, 20);
}
