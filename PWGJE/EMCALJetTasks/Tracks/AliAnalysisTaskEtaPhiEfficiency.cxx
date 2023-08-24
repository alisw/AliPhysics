/*
 * AliAnalysisTaskEtaPhiEfficiency.cxx
 *
 *  Created on: Feb 4, 2016
 *      Author: markus
 */

#include <TChain.h>
#include <TFile.h>
#include <THashList.h>
#include <THistManager.h>
#include <TKey.h>
#include <TList.h>
#include <TMath.h>
#include <TProfile.h>
#include <TString.h>
#include <TSystem.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliVEvent.h"

#include "AliAnalysisTaskEtaPhiEfficiency.h"

using namespace PWGJE::EMCALJetTasks;

ClassImp(AliAnalysisTaskEtaPhiEfficiency)

AliAnalysisTaskEtaPhiEfficiency::AliAnalysisTaskEtaPhiEfficiency():
  AliAnalysisTaskSE(),
  fAnalysisUtil(NULL),
  fHistos(NULL),
  fTrackCuts(NULL)
{
}

AliAnalysisTaskEtaPhiEfficiency::AliAnalysisTaskEtaPhiEfficiency(const char *name):
  AliAnalysisTaskSE(name),
  fAnalysisUtil(NULL),
  fHistos(NULL),
  fTrackCuts(NULL)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskEtaPhiEfficiency::~AliAnalysisTaskEtaPhiEfficiency() {
  if(fAnalysisUtil) delete fAnalysisUtil;
  if(fTrackCuts) delete fTrackCuts;
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskEtaPhiEfficiency::UserCreateOutputObjects(){
  fAnalysisUtil = new AliAnalysisUtils;

  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
  fTrackCuts->SetName("Standard Track cuts");
  fTrackCuts->SetMinNCrossedRowsTPC(120);
  fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

  fHistos = new THistManager("EfficiencyMaps");
  fHistos->CreateTH1("hNtrials", "Number of trials", 1, 0.5, 1.5);
  fHistos->CreateTProfile("hCrossSection", "Cross section", 1, 0.5, 1.5);

  // make efficiency maps
  for(int i = 10; i < 100; i += 10){
    fHistos->CreateTH2(Form("hParticles_%d_%d", i, i + 10), Form("Efficiency Map between %d and %d GeV/c for true particles", i, i + 10), 100, -0.8, 0.8, 100, 0., 2 * TMath::Pi());
    fHistos->CreateTH2(Form("hTracks_%d_%d", i, i + 10), Form("Efficiency Map between %d and %d GeV/c for reconstructed tracks", i, i + 10), 100, -0.8, 0.8, 100, 0., 2 * TMath::Pi());
  }

  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEtaPhiEfficiency::UserExec(Option_t *){
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return;
  if(fAnalysisUtil->IsFirstEventInChunk(InputEvent())) return;
  if(fAnalysisUtil->IsPileUpEvent(InputEvent())) return;
  if(!fAnalysisUtil->IsVertexSelected2013pA(InputEvent())) return;

  AliESDtrack *track = NULL;
  AliMCParticle *part = 0;
  Int_t ptmin, ptmax;
  for(int ipart = 0; ipart < MCEvent()->GetNumberOfTracks(); ipart++){
    part = static_cast<AliMCParticle *>(MCEvent()->GetTrack(ipart));
    if(!part->Charge()) continue;
    if(!MCEvent()->IsPhysicalPrimary(ipart)) continue;
    if(TMath::Abs(part->Eta()) > 0.8) continue;
    if(!FindPtBin(TMath::Abs(part->Pt()), ptmin, ptmax)) continue;
    fHistos->FillTH2(Form("hParticles_%d_%d", ptmin, ptmax), part->Eta(), part->Phi());
  }
  for(int itrk = 0; itrk < InputEvent()->GetNumberOfTracks(); itrk++){
    track = static_cast<AliESDtrack *>(InputEvent()->GetTrack(itrk));
    part = static_cast<AliMCParticle *>(MCEvent()->GetTrack(TMath::Abs(track->GetLabel())));
    if(!MCEvent()->IsPhysicalPrimary(TMath::Abs(part->GetLabel()))) continue;
    if(TMath::Abs(part->Eta()) > 0.8) continue;
    if(!fTrackCuts->AcceptTrack(track)) continue;
    if(!FindPtBin(TMath::Abs(part->Pt()), ptmin, ptmax)) continue;
    fHistos->FillTH2(Form("hTracks_%d_%d", ptmin, ptmax), part->Eta(), part->Phi());

  }
}

Bool_t AliAnalysisTaskEtaPhiEfficiency::UserNotify(){
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

  PythiaInfoFromFile(curfile->GetName(), xsection, trials, pthardbin);

  fHistos->FillTH1("hNtrials", 1., trials);
  fHistos->FillProfile("hCrossSection", 1., xsection);
  //fHistos->FillTH1(pthardbin, nevents);

  return kTRUE;
}

Bool_t AliAnalysisTaskEtaPhiEfficiency::PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard) const {

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

Bool_t AliAnalysisTaskEtaPhiEfficiency::FindPtBin(Double_t ptin, Int_t &ptmin, Int_t &ptmax) const{
  if(ptin < 10 || ptin > 100) return false;
  for(double ptiter = 10; ptiter < 100; ptiter += 10){
    if(ptin >= ptiter && ptin < ptiter + 10){
      ptmin = static_cast<Int_t>(ptiter);
      ptmax = static_cast<Int_t>(ptiter + 10);
      break;
    }
  }
  return true;
}
