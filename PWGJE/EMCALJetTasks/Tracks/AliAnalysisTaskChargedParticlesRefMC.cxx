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
#include <vector>
#include <map>

#include <TArrayD.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <THashList.h>
#include <TH1.h>
#include <TKey.h>
#include <TList.h>
#include <TProfile.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliEmcalTriggerPatchInfo.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliVVertex.h"

#include "AliEMCalHistoContainer.h"
#include "AliAnalysisTaskChargedParticlesRefMC.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy constructor
 */
AliAnalysisTaskChargedParticlesRefMC::AliAnalysisTaskChargedParticlesRefMC():
        AliAnalysisTaskSE(),
        fTrackCuts(NULL),
        fAnalysisUtil(NULL),
        fHistos(NULL),
        fPtHard(0),
        fPtHardBin(0),
        fNTrials(0),
        fXsection(0)
{
}

/**
 * Main constructor
 * @param name Name of the task
 */
AliAnalysisTaskChargedParticlesRefMC::AliAnalysisTaskChargedParticlesRefMC(const char* name):
        AliAnalysisTaskSE(name),
        fTrackCuts(NULL),
        fAnalysisUtil(NULL),
        fHistos(NULL),
        fPtHard(0),
        fPtHardBin(0),
        fNTrials(0),
        fXsection(0)
{
  DefineOutput(1, TList::Class());
}

/**
 * Destuctor
 */
AliAnalysisTaskChargedParticlesRefMC::~AliAnalysisTaskChargedParticlesRefMC() {
  if(fTrackCuts) delete fTrackCuts;
  if(fAnalysisUtil) delete fAnalysisUtil;
  if(fHistos) delete fHistos;
}
/**
 * Create the output histograms
 */
void AliAnalysisTaskChargedParticlesRefMC::UserCreateOutputObjects() {
  fAnalysisUtil = new AliAnalysisUtils;

  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
  fTrackCuts->SetName("Standard Track cuts");
  fTrackCuts->SetMinNCrossedRowsTPC(120);
  fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

  TArrayD oldbinning, newbinning;
  CreateOldPtBinning(oldbinning);
  CreateNewPtBinning(newbinning);

  fHistos = new AliEMCalHistoContainer("Ref");
  fHistos->CreateTH1("hNtrials", "Number of trials", 1, 0.5, 1.5);
  fHistos->CreateTProfile("hCrossSection", "PYTHIA cross section", 1, 0.5, 1.5);
  fHistos->CreateTH1("hNtrialsEvent", "Number of trials (from header, after event selection)", 1, 0.5, 1.5);
  fHistos->CreateTProfile("hCrossSectionEvent", "PYTHIA cross section (from header, after event selection)", 1, 0.5, 1.5);
  fHistos->CreateTH1("hPtHard", "Pt of the hard interaction", 1000, 0., 500);
  TString triggers[6] = {"True", "MB", "EJ1", "EJ2", "EG1", "EG2"};
  for(TString *trg = triggers; trg < triggers + sizeof(triggers)/sizeof(TString); trg++){
    fHistos->CreateTH1(Form("hEventCount%s", trg->Data()), Form("Event Counter for trigger class %s", trg->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hVertexBefore%s", trg->Data()), Form("Vertex distribution before z-cut for trigger class %s", trg->Data()), 500, -50, 50);
    fHistos->CreateTH1(Form("hVertexAfter%s", trg->Data()), Form("Vertex distribution after z-cut for trigger class %s", trg->Data()), 100, -10, 10);
    fHistos->CreateTH1(Form("hPtEtaAllOldBinning%s", trg->Data()), Form("Charged particle pt distribution all eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEtaCentOldBinning%s", trg->Data()), Form("Charged particle pt distribution central eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEtaAllNewBinning%s", trg->Data()), Form("Charged particle pt distribution all eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEtaCentNewBinning%s", trg->Data()), Form("Charged particle pt distribution central eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hEtaDistAllPt1%s", trg->Data()), Form("Eta distribution without etacut for tracks with Pt above 1 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaDistAllPt2%s", trg->Data()), Form("Eta distribution without etacut for tracks with Pt above 2 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaDistAllPt5%s", trg->Data()), Form("Eta distribution without etacut for tracks with Pt above 5 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaDistAllPt10%s", trg->Data()), Form("Eta distribution without etacut for tracks with Pt above 10 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaDistAllPt20%s", trg->Data()), Form("Eta distribution without etacut for tracks with Pt above 20 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaDistCutPt1%s", trg->Data()), Form("Eta distribution with etacut for tracks with Pt above 1 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaDistCutPt2%s", trg->Data()), Form("Eta distribution with etacut for tracks with Pt above 2 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaDistCutPt5%s", trg->Data()), Form("Eta distribution with etacut for tracks with Pt above 5 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaDistCutPt10%s", trg->Data()), Form("Eta distribution with etacut for tracks with Pt above 10 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
    fHistos->CreateTH1(Form("hEtaDistCutPt20%s", trg->Data()), Form("Eta distribution with etacut for tracks with Pt above 20 GeV/c trigger %s", trg->Data()), 100, -1., 1.);
  }
  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskChargedParticlesRefMC::UserExec(Option_t*) {  // Select event
  if(!fMCEvent) return;
  TClonesArray *fTriggerPatches = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));
  if(!fTriggerPatches) return;

  TString triggerstring = GetFiredTriggerClasses(fTriggerPatches);
  Bool_t isMinBias = fInputHandler->IsEventSelected() & AliVEvent::kINT7,
      isEJ1 = triggerstring.Contains("EJ1"),
      isEJ2 = triggerstring.Contains("EJ2"),
      isEG1 = triggerstring.Contains("EG1"),
      isEG2 = triggerstring.Contains("EG2");
  if(!(isMinBias || isEG1 || isEG2 || isEJ1 || isEJ2)) return;
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1) return;
  // Fill reference distribution for the primary vertex before any z-cut
  fHistos->FillTH1("hVertexBeforeTrue", vtx->GetZ());
  if(isMinBias) fHistos->FillTH1("hVertexBeforeMB", vtx->GetZ());
  if(isEJ1) fHistos->FillTH1("hVertexBeforeEJ1", vtx->GetZ());
  if(isEJ2) fHistos->FillTH1("hVertexBeforeEJ2", vtx->GetZ());
  if(isEG1) fHistos->FillTH1("hVertexBeforeEG1", vtx->GetZ());
  if(isEG2) fHistos->FillTH1("hVertexBeforeEG2", vtx->GetZ());
  if(!fAnalysisUtil->IsVertexSelected2013pA(fInputEvent)) return;       // Apply new vertex cut
  if(fAnalysisUtil->IsPileUpEvent(fInputEvent)) return;       // Apply new vertex cut
  // Apply vertex z cut
  if(vtx->GetZ() < -10. || vtx->GetZ() > 10.) return;

  // Fill Event counter and reference vertex distributions for the different trigger classes
  fHistos->FillTH1("hEventCountTrue", 1);
  fHistos->FillTH1("hVertexAfterTrue", vtx->GetZ());
  if(isMinBias){
    fHistos->FillTH1("hEventCountMB", 1);
    fHistos->FillTH1("hVertexAfterMB", vtx->GetZ());
  }
  if(isEJ1){
    fHistos->FillTH1("hEventCountEJ1", 1);
    fHistos->FillTH1("hVertexAfterEJ1", vtx->GetZ());
  }
  if(isEJ2){
    fHistos->FillTH1("hEventCountEJ2", 1);
    fHistos->FillTH1("hVertexAfterEJ2", vtx->GetZ());
  }
  if(isEG1){
    fHistos->FillTH1("hEventCountEG1", 1);
    fHistos->FillTH1("hVertexAfterEG1", vtx->GetZ());
  }
  if(isEG2){
    fHistos->FillTH1("hEventCountEG2", 1);
    fHistos->FillTH1("hVertexAfterEG2", vtx->GetZ());
  }

  // Fill PYTHIA histograms from event header
  AliGenPythiaEventHeader *pyheader = GetPythiaHeader();
  if(pyheader){
   fHistos->FillTH1("hNtrialsEvent", 1., pyheader->Trials());
   fHistos->FillProfile("hCrossSectionEvent", 1., pyheader->GetXsection());
   fHistos->FillTH1("hPtHard", pyheader->GetPtHard());
  }

  // MonteCarlo Loop
  // Histograms
  // - Full eta (-0.8, 0.8), new binning
  // - Full eta (-0.8, 0.8), old binning
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c without eta cut
  // - Central eta (-0.8, -0.2), new binning,
  // - Central eta (-0.8, -0.2), old binning,
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c with eta cut
  int ptmin[5] = {1,2,5,10,20}; // for eta distributions
  AliVParticle *truepart = NULL;
  for(int ipart = 0; ipart < fMCEvent->GetNumberOfTracks(); ipart++){
    truepart = fMCEvent->GetTrack(ipart);

    // Select only particles within ALICE acceptance
    if(TMath::Abs(truepart->Eta()) > 0.8) continue;
    if(TMath::Abs(truepart->Pt()) < 0.1) continue;
    if(!truepart->Charge()) continue;

    if(!IsPhysicalPrimary(truepart, fMCEvent)) continue;

    // Particle selected
    fHistos->FillTH1("hPtEtaAllNewBinningTrue", TMath::Abs(truepart->Pt()));
    fHistos->FillTH1("hPtEtaAllOldBinningTrue", TMath::Abs(truepart->Pt()));

    for(int icut = 0; icut < 5; icut++){
      if(TMath::Abs(truepart->Pt()) > static_cast<double>(ptmin[icut])){
        fHistos->FillTH1(Form("hEtaDistAllPt%dTrue", ptmin[icut]), truepart->Eta());
      }
    }

    if(truepart->Eta() > -0.8 && truepart->Eta() < -0.2){
      fHistos->FillTH1("hPtEtaCentNewBinningTrue", TMath::Abs(truepart->Pt()));
      fHistos->FillTH1("hPtEtaCentOldBinningTrue", TMath::Abs(truepart->Pt()));
      for(int icut = 0; icut < 5; icut++){
        if(TMath::Abs(truepart->Pt()) > static_cast<double>(ptmin[icut])){
          fHistos->FillTH1(Form("hEtaDistCutPt%dTrue", ptmin[icut]), truepart->Eta());
        }
      }
    }
  }

  // Loop over tracks, fill select particles
  // Histograms
  // - Full eta (-0.8, 0.8), new binning
  // - Full eta (-0.8, 0.8), old binning
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c without eta cut
  // - Central eta (-0.8, -0.2), new binning,
  // - Central eta (-0.8, -0.2), old binning,
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c with eta cut
  AliVTrack *checktrack(NULL);
  AliVParticle *assocMC(NULL);
  double ptparticle(-1.), etaparticle(-100.);
  for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); ++itrk){
    checktrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(itrk));
    if(!checktrack) continue;
    // Find associated particle
    assocMC = fMCEvent->GetTrack(TMath::Abs(checktrack->GetLabel()));
    if(!assocMC) continue;        // Fake track
    if(!IsPhysicalPrimary(assocMC, fMCEvent)) continue;

    // Select only particles within ALICE acceptance
    if(TMath::Abs(checktrack->Eta()) > 0.8) continue;
    if(TMath::Abs(checktrack->Pt()) < 0.1) continue;

    // Distinguish track selection for ESD and AOD tracks
    AliESDtrack *esdtrack(NULL);
    AliAODTrack *aodtrack(NULL);
    if((esdtrack = dynamic_cast<AliESDtrack *>(checktrack))){
      if(!TrackSelectionESD(esdtrack)) continue;
    } else if((aodtrack = dynamic_cast<AliAODTrack *>(checktrack))){
      if(!TrackSelectionAOD(aodtrack)) continue;
    } else {
      continue;
    }

    ptparticle = TMath::Abs(assocMC->Pt());
    etaparticle = assocMC->Eta();

    // fill histograms allEta
    if(isMinBias){
      fHistos->FillTH1("hPtEtaAllNewBinningMB", ptparticle);
      fHistos->FillTH1("hPtEtaAllOldBinningMB", ptparticle);
      for(int icut = 0; icut < 5; icut++){
        if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
          fHistos->FillTH1(Form("hEtaDistAllPt%dMB", ptmin[icut]), etaparticle);
        }
      }
    }
    if(isEJ1){
      fHistos->FillTH1("hPtEtaAllNewBinningEJ1", ptparticle);
      fHistos->FillTH1("hPtEtaAllOldBinningEJ1", ptparticle);
      for(int icut = 0; icut < 5; icut++){
        if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
          fHistos->FillTH1(Form("hEtaDistAllPt%dEJ1", ptmin[icut]), etaparticle);
        }
      }
    }
    if(isEJ2){
      fHistos->FillTH1("hPtEtaAllNewBinningEJ2", ptparticle);
      fHistos->FillTH1("hPtEtaAllOldBinningEJ2", ptparticle);
      for(int icut = 0; icut < 5; icut++){
        if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
          fHistos->FillTH1(Form("hEtaDistAllPt%dEJ2", ptmin[icut]), etaparticle);
        }
      }
    }
    if(isEG1){
      fHistos->FillTH1("hPtEtaAllNewBinningEG1", ptparticle);
      fHistos->FillTH1("hPtEtaAllOldBinningEG1", ptparticle);
      for(int icut = 0; icut < 5; icut++){
        if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
          fHistos->FillTH1(Form("hEtaDistAllPt%dEG1", ptmin[icut]), etaparticle);
        }
      }
    }
    if(isEG2){
      fHistos->FillTH1("hPtEtaAllNewBinningEG2", ptparticle);
      fHistos->FillTH1("hPtEtaAllOldBinningEG2", ptparticle);
      for(int icut = 0; icut < 5; icut++){
        if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
          fHistos->FillTH1(Form("hEtaDistAllPt%dEG2", ptmin[icut]), etaparticle);
        }
      }
    }

    if(checktrack->Eta() > -0.8 && checktrack->Eta() < -0.2){
      // Fill Histograms in central eta
      if(isMinBias){
        fHistos->FillTH1("hPtEtaCentNewBinningMB", ptparticle);
        fHistos->FillTH1("hPtEtaCentOldBinningMB", ptparticle);
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaDistCutPt%dMB", ptmin[icut]), etaparticle);
          }
        }
      }
      if(isEJ1){
        fHistos->FillTH1("hPtEtaCentNewBinningEJ1", ptparticle);
        fHistos->FillTH1("hPtEtaCentOldBinningEJ1", ptparticle);
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaDistCutPt%dEJ1", ptmin[icut]), etaparticle);
          }
        }
      }
      if(isEJ2){
        fHistos->FillTH1("hPtEtaCentNewBinningEJ2", ptparticle);
        fHistos->FillTH1("hPtEtaCentOldBinningEJ2", ptparticle);
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaDistCutPt%dEJ2", ptmin[icut]), etaparticle);
          }
        }
      }
      if(isEG1){
        fHistos->FillTH1("hPtEtaCentNewBinningEG1", ptparticle);
        fHistos->FillTH1("hPtEtaCentOldBinningEG1", ptparticle);
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaDistCutPt%dEG1", ptmin[icut]), etaparticle);
          }
        }
      }
      if(isEG2){
        fHistos->FillTH1("hPtEtaCentNewBinningEG2", ptparticle);
        fHistos->FillTH1("hPtEtaCentOldBinningEG2", ptparticle);
        for(int icut = 0; icut < 5; icut++){
          if(TMath::Abs(checktrack->Pt()) > static_cast<double>(ptmin[icut])){
            fHistos->FillTH1(Form("hEtaDistCutPt%dEG2", ptmin[icut]), etaparticle);
          }
        }
      }
    }
  }
  PostData(1, fHistos->GetListOfHistograms());

}

/**
 * Perform actions when giles change
 * @return
 */
Bool_t AliAnalysisTaskChargedParticlesRefMC::UserNotify()
{
  // Called when file changes.

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliError(Form("%s - UserNotify: No current tree!",GetName()));
    return kFALSE;
  }

  Float_t xsection    = 0;
  Float_t trials      = 0;
  Int_t   pthardbin   = 0;

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliError(Form("%s - UserNotify: No current file!",GetName()));
    return kFALSE;
  }

  TChain *chain = dynamic_cast<TChain*>(tree);
  if (chain) tree = chain->GetTree();

  Int_t nevents = tree->GetEntriesFast();

  PythiaInfoFromFile(curfile->GetName(), xsection, trials, pthardbin);

  fHistos->FillTH1("hNtrials", 1., trials);
  fHistos->FillProfile("hCrossSection", 1., xsection);
  //fHistos->FillTH1(pthardbin, nevents);

  return kTRUE;
}

/**
 * Get the cross section and the trails either from pyxsec.root or from pysec_hists.root
 * Get the pt hard bin from the file path
 *  This is to called in Notify and should provide the path to the AOD/ESD file
 * (Partially copied from AliAnalysisHelperJetTasks)
 * From AliAnalysisTaskEmcal
 * @param currFile File name with PYTHIA hard cross section
 * @param fXsec Output storage for the cross section
 * @param fTrials Output storage for the number of trials
 * @param pthard Output storage of the pthardbin
 * @return True if reading was successful, false in case of errors
 */
Bool_t AliAnalysisTaskChargedParticlesRefMC::PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard) const {

  TString file(currFile);
  fXsec = 0;
  fTrials = 1;

  if (file.Contains(".zip#")) {
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  AliDebug(1,Form("File name: %s",file.Data()));

  // Get the pt hard bin
  TString strPthard(file);

  strPthard.Remove(strPthard.Last('/'));
  strPthard.Remove(strPthard.Last('/'));
  if (strPthard.Contains("AOD")) strPthard.Remove(strPthard.Last('/'));
  strPthard.Remove(0,strPthard.Last('/')+1);
  if (strPthard.IsDec())
    pthard = strPthard.Atoi();
  else
    AliWarning(Form("Could not extract file number from path %s", strPthard.Data()));

  // problem that we cannot really test the existance of a file in a archive so we have to live with open error message from root
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root"));

  if (!fxsec) {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if (!fxsec) {
  // not a severe condition but inciate that we have no information
      return kFALSE;
    } else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
      if (!key) {
        fxsec->Close();
        return kFALSE;
      }
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) {
        fxsec->Close();
        return kFALSE;
      }
      fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } else { // no tree pyxsec.root
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if (!xtree) {
      fxsec->Close();
      return kFALSE;
    }
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    fTrials = ntrials;
    fXsec = xsection;
    fxsec->Close();
  }
  return kTRUE;
}

/**
 * Access PYTHIA event header
 * @return pythia event header (if existing)
 */
AliGenPythiaEventHeader *AliAnalysisTaskChargedParticlesRefMC::GetPythiaHeader() const {
  AliGenPythiaEventHeader *pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
  if (!pythiaHeader) {
    // Check if AOD
    AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));
    if (aodMCH) {
      for (UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
        pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
        if (pythiaHeader) break;
      }
    }
  }
  return pythiaHeader;
}

/**
 * Create old pt binning
 * @param binning
 */
void AliAnalysisTaskChargedParticlesRefMC::CreateOldPtBinning(TArrayD& binning) const {
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double,double>(2.5, 0.1));
  definitions.insert(std::pair<double,double>(7., 0.25));
  definitions.insert(std::pair<double,double>(15., 0.5));
  definitions.insert(std::pair<double,double>(25., 1.));
  definitions.insert(std::pair<double,double>(40., 2.5));
  definitions.insert(std::pair<double,double>(50., 5.));
  definitions.insert(std::pair<double,double>(100., 10.));
  double currentval = 0;
  for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
    double limit = id->first, binwidth = id->second;
    while(currentval < limit){
      currentval += binwidth;
      mybinning.push_back(currentval);
     }
   }
   binning.Set(mybinning.size());
   int ib = 0;
   for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
     binning[ib++] = *it;
}

/**
 * Create new Pt binning
 * @param binning
 */
void AliAnalysisTaskChargedParticlesRefMC::CreateNewPtBinning(TArrayD& binning) const {
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double, double>(1, 0.05));
  definitions.insert(std::pair<double, double>(2, 0.1));
  definitions.insert(std::pair<double, double>(4, 0.2));
  definitions.insert(std::pair<double, double>(7, 0.5));
  definitions.insert(std::pair<double, double>(16, 1));
  definitions.insert(std::pair<double, double>(36, 2));
  definitions.insert(std::pair<double, double>(40, 4));
  definitions.insert(std::pair<double, double>(50, 5));
  definitions.insert(std::pair<double, double>(100, 10));
  definitions.insert(std::pair<double, double>(200, 20));
  double currentval = 0.;
  mybinning.push_back(currentval);
  for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
    double limit = id->first, binwidth = id->second;
    while(currentval < limit){
      currentval += binwidth;
      mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

/**
 * Run track selection for ESD tracks
 * @param track The track to check
 * @return True if the track is selected, false otherwise
 */
Bool_t AliAnalysisTaskChargedParticlesRefMC::TrackSelectionESD(AliESDtrack* track) {
  return fTrackCuts->AcceptTrack(track);
}

/**
 * Run track selection for AOD tracks
 * @param track The track to check
 * @return True if the track is selected, false otherwise
 */
Bool_t AliAnalysisTaskChargedParticlesRefMC::TrackSelectionAOD(AliAODTrack* track) {
  if(!track->TestFilterBit(AliAODTrack::kTrkGlobal)) return false;
  if(track->GetTPCNCrossedRows() < 120) return false;
  return true;
}

/**
 * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
 * @param triggerpatches Trigger patches found by the trigger maker
 * @return String with EMCAL trigger decision
 */
TString AliAnalysisTaskChargedParticlesRefMC::GetFiredTriggerClasses(const TClonesArray* triggerpatches) {
  TString triggerstring = "";
  Int_t nEJ1 = 0, nEJ2 = 0, nEG1 = 0, nEG2 = 0;
  double  minADC_EJ1 = 260.,
          minADC_EJ2 = 127.,
          minADC_EG1 = 140.,
          minADC_EG2 = 89.;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEmcalTriggerPatchInfo *patch = dynamic_cast<AliEmcalTriggerPatchInfo *>(*patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(patch->IsJetHighSimple() && patch->GetADCOfflineAmp() > minADC_EJ1) nEJ1++;
    if(patch->IsJetLowSimple() && patch->GetADCOfflineAmp() > minADC_EJ2) nEJ2++;
    if(patch->IsGammaHighSimple() && patch->GetADCOfflineAmp() > minADC_EG1) nEG1++;
    if(patch->IsGammaLowSimple() && patch->GetADCOfflineAmp() > minADC_EG2) nEG2++;
  }
  if(nEJ1) triggerstring += "EJ1";
  if(nEJ2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EJ2";
  }
  if(nEG1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EG1";
  }
  if(nEG2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EG2";
  }
  return triggerstring;
}

/**
 * Check in a transparent way for ESDs and AODs whether the particle is physical primary or not
 * -# AOD: Information stored in the AliAODMCParticle
 * -# ESD: Information needs to be retrieved from the stack via the label of the MC particle
 * @param part The particle to check
 * @param mcevent The MC event containing the stack (ESD only)
 * @return True if particle is a physical primary particle, false otherwise
 */
Bool_t AliAnalysisTaskChargedParticlesRefMC::IsPhysicalPrimary(const AliVParticle* const part, AliMCEvent* const mcevent) {
  Bool_t physprim = false;
  const AliAODMCParticle *aodmc = dynamic_cast<const AliAODMCParticle *>(part);
  if(aodmc){
    physprim = aodmc->IsPhysicalPrimary();
  } else {
    physprim = mcevent->IsPhysicalPrimary(part->GetLabel());
  }
  return physprim;
}

} /* namespace EMCalTriggerPtAnalysis */
