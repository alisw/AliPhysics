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
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <THashList.h>
#include <TKey.h>
#include <TList.h>
#include <TMath.h>
#include <TParticlePDG.h>
#include <TPDGCode.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TTree.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliAODHeader.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliEMCalHistoContainer.h"
#include "AliGenPythiaEventHeader.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliVEvent.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include "AliAnalysisTaskTracksInJet.h"

ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskTracksInJet)

namespace EMCalTriggerPtAnalysis {

/**
 * Constructor
 */
AliAnalysisTaskTracksInJet::AliAnalysisTaskTracksInJet() :
    AliAnalysisTaskSE(),
    fJetStructure(),
    fJetTree(NULL),
    fAnalysisUtils(NULL),
    fTrackCuts(NULL),
    fHybridCuts(NULL),
    fIsMC(kFALSE),
    fFracPtHard(-1),
    fHistosMC(NULL)
{
}

AliAnalysisTaskTracksInJet::AliAnalysisTaskTracksInJet(const char *taskname) :
    AliAnalysisTaskSE(taskname),
    fJetStructure(),
    fJetTree(NULL),
    fAnalysisUtils(NULL),
    fTrackCuts(NULL),
    fHybridCuts(NULL),
    fIsMC(kFALSE),
    fFracPtHard(-1),
    fHistosMC(NULL)
{
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

AliAnalysisTaskTracksInJet::~AliAnalysisTaskTracksInJet() {}

void AliAnalysisTaskTracksInJet::UserCreateOutputObjects(){
  fHistosMC = new AliEMCalHistoContainer("MCHistos");
  fHistosMC->CreateTH1("hNtrials", "Number of trials", 1., 0.5, 1.5);
  fHistosMC->CreateTH1("hNtrialsEvent", "Number of trials (from header, after event selection)", 1, 0.5, 1.5);
  fHistosMC->CreateTH1("hPtHard", "Pt of the hard interaction", 1000, 0., 500);
  fHistosMC->CreateTProfile("hCrossSection", "PYTHIA Cross section", 1, 0.5, 1.5);

  OpenFile(2);
  fJetTree = new TTree("JetTreeRec", "Jet Tree on reconstructed data");
  fJetTree->Branch("JetsRec", &fJetStructure, "pjet[3]/D:plead[3]/D:psub[3]/D:isData/I");

  fAnalysisUtils = new AliAnalysisUtils;

  // @TODO: Instance Hybrid track cuts for ESDs
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  fTrackCuts->SetMaxDCAToVertexXY(2.4);
  fTrackCuts->SetMaxDCAToVertexZ(3.2);
  fTrackCuts->SetDCAToVertex2D(kTRUE);
  fTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  fTrackCuts->SetMaxFractionSharedTPCClusters(0.4);

  // Second class of hybrid track cuts
  fHybridCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  fHybridCuts->SetMaxDCAToVertexXY(2.4);
  fHybridCuts->SetMaxDCAToVertexZ(3.2);
  fHybridCuts->SetDCAToVertex2D(kTRUE);
  fHybridCuts->SetMaxChi2TPCConstrainedGlobal(36);
  fHybridCuts->SetMaxFractionSharedTPCClusters(0.4);
  fHybridCuts->SetRequireITSRefit(kFALSE);
  fHybridCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);

  PostData(1, fHistosMC->GetListOfHistograms());
  PostData(2, fJetTree);
}

/**
 * Perform actions when giles change
 * @return
 */
Bool_t AliAnalysisTaskTracksInJet::UserNotify()
{
  // Called when file changes.

  if(fIsMC){
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

    fHistosMC->FillTH1("hNtrials", 1., trials);
    fHistosMC->FillProfile("hCrossSection", 1., xsection);
    //fHistos->FillTH1(pthardbin, nevents);
  }

  return kTRUE;
}

/**
 * Event loop
 * @param
 */
void AliAnalysisTaskTracksInJet::UserExec(Option_t *){
  // 1st Select event
  if(!fAnalysisUtils->IsVertexSelected2013pA(InputEvent())) return;
  if(fAnalysisUtils->IsPileUpEvent(InputEvent())) return;

  if(MCEvent()){
    // Apply outlier cut
    AliGenPythiaEventHeader *pyheader = GetPythiaHeader();
    if(fFracPtHard > 0){
      if(pyheader && IsOutlier(pyheader))
        return;
    }
    fHistosMC->FillTH1("hNtrialsEvent", pyheader->Trials());
    fHistosMC->FillTH1("hPtHard", pyheader->GetPtHard());
  }

  Double_t pvecjet[3], pveclead[3], pvecsublead[3];
  if(MCEvent()){
    // Prepare MC jet finding
    std::vector<fastjet::PseudoJet> mcpseudo;
    for(int itrk = 0; itrk < fMCEvent->GetNumberOfTracks(); itrk++){
      AliVParticle *mctrack = fMCEvent->GetTrack(itrk);
      if(!mctrack->Charge()) continue;
      if(TMath::Abs(mctrack->Pt()) < 0.15) continue;
      if(TMath::Abs(mctrack->Eta()) > 0.8) continue;
      if(!IsPhysicalPrimary(mctrack, MCEvent())) continue;

      fastjet::PseudoJet inputjet(mctrack->Px(), mctrack->Py(), mctrack->Pz(), mctrack->E());
      inputjet.set_user_index(itrk);
      mcpseudo.push_back(inputjet);
    }

    // Run jet finding
    fastjet::ClusterSequence jetfinder(mcpseudo, fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4));
    std::vector<fastjet::PseudoJet> genjets = sorted_by_pt(jetfinder.inclusive_jets());

    for(std::vector<fastjet::PseudoJet>::const_iterator jetit = genjets.begin(); jetit != genjets.end(); ++jetit){
      if(TMath::Abs(jetit->eta()) > 0.5) continue;       // expect jet to be within |eta| < +- 0.5
      std::vector<fastjet::PseudoJet> constituentssorted = sorted_by_pt(jetit->constituents());
      memset(pveclead, 0, sizeof(Double_t) * 3);
      memset(pvecsublead, 0, sizeof(Double_t) * 3);
      pvecjet[0] = jetit->px();
      pvecjet[1] = jetit->py();
      pvecjet[2] = jetit->pz();
      int icounter = 0;
      for(std::vector<fastjet::PseudoJet>::const_iterator cit = constituentssorted.begin(); cit != constituentssorted.end(); ++cit){
        AliVParticle *basepart = fMCEvent->GetTrack(cit->user_index());
        if(icounter == 0) {
          basepart->PxPyPz(pveclead);
        } else if(icounter == 1) {
          basepart->PxPyPz(pvecsublead);
        } else {
          break;
        }
        icounter++;
      }

      if(icounter >= 2){
        // At least 2 particles, accept jet
        fJetStructure.Reset();
        fJetStructure.fIsData = 0;
        memcpy(fJetStructure.fPvecJet, pvecjet, sizeof(Double_t) * 3);
        memcpy(fJetStructure.fPvecLead, pveclead, sizeof(Double_t) * 3);
        memcpy(fJetStructure.fPvecSubLead, pvecsublead, sizeof(Double_t) * 3);
        fJetTree->Fill();
      }
    }
  }

  // Select tracks
  TVector3 pvectrack;
  std::vector<fastjet::PseudoJet> datapseudo;
  AliESDtrack *esdtrack(NULL);
  AliAODTrack *aodtrack(NULL);
  double mpion = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
  for(int ipart = 0; ipart < fInputEvent->GetNumberOfTracks(); ipart++){
    AliVParticle *recpart = fInputEvent->GetTrack(ipart);
    if(MCEvent()){
      AliVParticle *assocMC = fMCEvent->GetTrack(TMath::Abs(recpart->GetLabel()));
      if(!assocMC) continue;
      if(!IsPhysicalPrimary(assocMC, fMCEvent)) continue;
    }

    // Select track
    if(TMath::Abs(recpart->Eta()) > 0.8) continue;
    if((esdtrack = dynamic_cast<AliESDtrack *>(recpart))){
      AliESDtrack copytrack(*esdtrack);
      if(!TrackSelectionESD(&copytrack)) continue;
      pvectrack.SetXYZ(copytrack.Px(), copytrack.Py(), copytrack.Pz());
    } else if((aodtrack = dynamic_cast<AliAODTrack *>(recpart))){
      if(!TrackSelectionAOD(aodtrack)) continue;
      pvectrack.SetXYZ(recpart->Px(), recpart->Py(), recpart->Pz());
    } else continue;

    // Create pseudojet under the assumption the particle is a pion
    fastjet::PseudoJet inputparticle(pvectrack.Px(), pvectrack.Py(), pvectrack.Pz(), TMath::Sqrt(pvectrack.Perp()*pvectrack.Perp() + mpion * mpion));
    inputparticle.set_user_index(ipart);
    datapseudo.push_back(inputparticle);
  }

  // Run jetfinder
  fastjet::ClusterSequence jetfinder(datapseudo, fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4));
  std::vector<fastjet::PseudoJet> recjets = sorted_by_pt(jetfinder.inclusive_jets());

  // Loop over jets, select jets with at least 2 tracks, store momentum vector of track, leading jet and subleading jet
  for(std::vector<fastjet::PseudoJet>::const_iterator jetit = recjets.begin(); jetit != recjets.end(); ++jetit){
    if(TMath::Abs(jetit->eta()) > 0.5) continue;       // expect jet to be within |eta| < +- 0.5
    std::vector<fastjet::PseudoJet> constituentssorted = sorted_by_pt(jetit->constituents());
    memset(pveclead, 0, sizeof(Double_t) * 3);
    memset(pvecsublead, 0, sizeof(Double_t) * 3);
    pvecjet[0] = jetit->px();
    pvecjet[1] = jetit->py();
    pvecjet[2] = jetit->pz();
    int icounter = 0;
    for(std::vector<fastjet::PseudoJet>::const_iterator cit = constituentssorted.begin(); cit != constituentssorted.end(); ++cit){
      AliVParticle *basepart = fInputEvent->GetTrack(cit->user_index());
      if(icounter == 0) {
        basepart->PxPyPz(pveclead);
      } else if(icounter == 1) {
        basepart->PxPyPz(pvecsublead);
      } else {
        break;
      }
      icounter++;
    }

    if(icounter >= 2){
      // At least 2 particles, accept jet
      fJetStructure.Reset();
      fJetStructure.fIsData = 1;
      memcpy(fJetStructure.fPvecJet, pvecjet, sizeof(Double_t) * 3);
      memcpy(fJetStructure.fPvecLead, pveclead, sizeof(Double_t) * 3);
      memcpy(fJetStructure.fPvecSubLead, pvecsublead, sizeof(Double_t) * 3);
      fJetTree->Fill();
    }
  }

  PostData(1, fHistosMC->GetListOfHistograms());
  PostData(2, fJetTree);
}

/**
 * Run track selection for ESD tracks. Does hybrid track selection. Remember to
 * do this on a copy
 * @param track The track to check
 * @return True if the track is selected, false otherwise
 */
Bool_t AliAnalysisTaskTracksInJet::TrackSelectionESD(AliESDtrack* track) const {
  if (fTrackCuts->AcceptTrack(track)) {
    track->SetBit(BIT(22),0);
    track->SetBit(BIT(23),0);
    return true;
  } else if (fHybridCuts->AcceptTrack(track)) {
    if (!track->GetConstrainedParam())
      return false;
    UInt_t status = track->GetStatus();
    if (((status&AliESDtrack::kITSrefit)==0))
      return false;
    const AliExternalTrackParam* constrainParam = track->GetConstrainedParam();
    track->Set(constrainParam->GetX(),
        constrainParam->GetAlpha(),
        constrainParam->GetParameter(),
        constrainParam->GetCovariance());
    if ((status&AliESDtrack::kITSrefit)==0) {
      track->SetBit(BIT(22),0); //type 2
      track->SetBit(BIT(23),1);
    } else {
      track->SetBit(BIT(22),1); //type 1
      track->SetBit(BIT(23),0);
    }
    return true;
  }

  return false;
}

/**
 * Run track selection for AOD tracks
 * @param track The track to check
 * @return True if the track is selected, false otherwise
 */
Bool_t AliAnalysisTaskTracksInJet::TrackSelectionAOD(AliAODTrack* track) const {
  // @TODO: Change to hybrid track cuts (256,512)
  if(!(track->TestFilterBit(256) || track->TestFilterBit(256))) return false;
  return true;
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
Bool_t AliAnalysisTaskTracksInJet::PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard) const {

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
 * Check in a transparent way for ESDs and AODs whether the particle is physical primary or not
 * -# AOD: Information stored in the AliAODMCParticle
 * -# ESD: Information needs to be retrieved from the stack via the label of the MC particle
 * @param part The particle to check
 * @param mcevent The MC event containing the stack (ESD only)
 * @return True if particle is a physical primary particle, false otherwise
 */
Bool_t AliAnalysisTaskTracksInJet::IsPhysicalPrimary(const AliVParticle* const part, AliMCEvent* const mcevent) const {
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
 * Access PYTHIA event header
 * @return pythia event header (if existing)
 */
AliGenPythiaEventHeader *AliAnalysisTaskTracksInJet::GetPythiaHeader() const {
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
 * Find outlier jets compared to the pt hard
 * @param header PYTHIA header with trigger jets and pt hard
 * @return True if event has at least one outlier, false otherwise
 */
Bool_t AliAnalysisTaskTracksInJet::IsOutlier(AliGenPythiaEventHeader * const header) const {
  Bool_t hasOutlier = kFALSE;
  Float_t pbuf[4];
  TLorentzVector jetvec;
  for(int ijet = 0; ijet < header->NTriggerJets(); ijet++){
    memset(pbuf, 0, sizeof(Float_t) * 4);
    header->TriggerJet(ijet, pbuf);
    jetvec.SetPxPyPzE(pbuf[0], pbuf[1], pbuf[2], pbuf[3]);
    if(TMath::Abs(jetvec.Pt()) >= fFracPtHard * header->GetPtHard()){
      hasOutlier = true;
      break;
    }
  }
  return hasOutlier;
}

} /* namespace EMCalTriggerPtAnalysis */
