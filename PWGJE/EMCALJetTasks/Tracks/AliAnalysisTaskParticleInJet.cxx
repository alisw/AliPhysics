/*
 * AliAnalysisTaskParticleInJet.cxx
 *
 *  Created on: Feb 17, 2016
 *      Author: markus
 */
#include <map>
#include <string>

#include <TArrayD.h>
#include <THashList.h>
#include <THistManager.h>
#include <TMath.h>
#include <TString.h>

#include "AliAODMCParticle.h"
#include "AliEmcalJet.h"
#include "AliEmcalTrackSelection.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliParticleContainer.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"

#include "AliAnalysisTaskParticleInJet.h"

ClassImp(AliAnalysisTaskParticleInJet)

AliAnalysisTaskParticleInJet::AliAnalysisTaskParticleInJet() :
AliAnalysisTaskEmcalJet(),
fHistMgr(NULL),
fTrackSelection(NULL),
fParticleContainerNameRec(""),
fParticleContainerNameMC(""),
fJetContainerNameRec(""),
fJetContainerNameMC("")
{
}

AliAnalysisTaskParticleInJet::AliAnalysisTaskParticleInJet(const char *name) :
AliAnalysisTaskEmcalJet(),
fHistMgr(NULL),
fTrackSelection(NULL),
fParticleContainerNameRec(""),
fParticleContainerNameMC(""),
fJetContainerNameRec(""),
fJetContainerNameMC("")
{
}

AliAnalysisTaskParticleInJet::~AliAnalysisTaskParticleInJet() {
  if(fTrackSelection) delete fTrackSelection;
}

void AliAnalysisTaskParticleInJet::UserCreateOutputObjects(){
  fHistMgr = new THistManager("histos");

  TArrayD ptbinning, jetptbinning;
  CreatePtBinning(ptbinning);
  CreateLinearBinning(jetptbinning, 200, 0., 200);

  std::map<std::string, std::string> histmap1D, histmap2D;
  histmap1D.insert(std::pair<std::string, std::string>("hMCall","MC true particles in full acceptance"));
  histmap1D.insert(std::pair<std::string, std::string>("hMCcont","Accepted true particles in particle containers"));
  histmap1D.insert(std::pair<std::string, std::string>("hTracksMB","All accepted tracks in the EMCAL acceptance in MB events"));
  histmap1D.insert(std::pair<std::string, std::string>("hTracksEG1","All accepted tracks in the EMCAL acceptance in EG1 events"));
  histmap1D.insert(std::pair<std::string, std::string>("hTracksEG2","All accepted tracks in the EMCAL acceptance in EG2 events"));
  histmap1D.insert(std::pair<std::string, std::string>("hTracksEJ1","All accepted tracks in the EMCAL acceptance in EJ1 events"));
  histmap1D.insert(std::pair<std::string, std::string>("hTracksEJ2","All accepted tracks in the EMCAL acceptance in EJ2 events"));
  histmap2D.insert(std::pair<std::string, std::string>("hMCjetTrack","Particles in jets with jet pt"));
  histmap2D.insert(std::pair<std::string, std::string>("hRecjetTrackMB","Particles in jets with jet pt in MB events"));
  histmap2D.insert(std::pair<std::string, std::string>("hRecjetTrackEG1","Particles in jets with jet pt in EG1 events"));
  histmap2D.insert(std::pair<std::string, std::string>("hRecjetTrackEG2","Particles in jets with jet pt in EG2 events"));
  histmap2D.insert(std::pair<std::string, std::string>("hRecjetTrackEJ1","Particles in jets with jet pt in EJ1 events"));
  histmap2D.insert(std::pair<std::string, std::string>("hRecjetTrackEJ2","Particles in jets with jet pt in EJ2 events"));

  for(std::map<std::string, std::string>::iterator it = histmap1D.begin(); it != histmap1D.end(); ++it){
    fHistMgr->CreateTH1(it->first.c_str(), it->second.c_str(), ptbinning);
  }
  for(std::map<std::string, std::string>::iterator it = histmap2D.begin(); it != histmap2D.end(); ++it){
    fHistMgr->CreateTH2(it->first.c_str(), it->second.c_str(), jetptbinning, ptbinning);
  }

  fHistMgr->GetListOfHistograms()->SetOwner(kFALSE);
  for(TIter histiter = TIter(fHistMgr->GetListOfHistograms()).Begin(); histiter != TIter::End(); ++histiter)
    fHistosQA->Add(*histiter);
}

Bool_t AliAnalysisTaskParticleInJet::Run(){
  TString triggerstring = fInputEvent->GetFiredTriggerClasses();
  Bool_t isMinBias = fInputHandler->IsEventSelected() & AliVEvent::kINT7,
      isEG1 = triggerstring.Contains("EG1"),
      isEG2 = triggerstring.Contains("EG2"),
      isEJ1 = triggerstring.Contains("EJ1"),
      isEJ2 = triggerstring.Contains("EJ2");
  if(!(isMinBias || isEG1 || isEG2 || isEJ1 || isEJ2)) return kFALSE;

  if(MCEvent()){
    AliVParticle *mctrack = 0;
    for(Int_t ipart = 0; ipart < MCEvent()->GetNumberOfTracks(); ipart++){
      mctrack = MCEvent()->GetTrack(ipart);
      if(!IsPhysicalPrimary(mctrack)) continue;
      if(!AcceptParticle(mctrack)) continue;
      fHistMgr->FillTH1("hMCall", mctrack->Pt());
    }
    std::vector<const AliVParticle *> particles = GetSelectedParticles(GetParticleContainer(fParticleContainerNameMC.Data()));
    for(std::vector<const AliVParticle *>::iterator it = particles.begin(); it != particles.end(); ++it) {
      const Double_t &pt = (*it)->Pt();
      fHistMgr->FillTH1("hMCcont", pt);
    }

    // Look at the MC Jet container
    GetJetContainer(fJetContainerNameMC.Data())->ResetCurrentID();
    AliEmcalJet *mcjet = GetJetContainer(fJetContainerNameMC.Data())->GetNextAcceptJet();
    do{
      for(int itrk = 0; itrk < mcjet->GetNumberOfTracks(); itrk++){
        mctrack = mcjet->TrackAt(itrk, GetParticleContainer(0)->GetArray());
        if(!AcceptParticle(mctrack)) continue;
        fHistMgr->FillTH2("hMCjetTrack", mcjet->Pt(), mctrack->Pt());
      }
    } while ((mcjet = GetJetContainer(fJetContainerNameMC.Data())->GetNextAcceptJet()));
  }

  // Loop over particles, select tracks, fill histograms of tracks in EMCAL
  std::vector<const AliVParticle *> tracks = GetSelectedParticles(GetParticleContainer(fParticleContainerNameRec.Data()));
  for(std::vector<const AliVParticle *>::iterator it = tracks.begin(); it != tracks.end(); ++it) {
    const Double_t &pt = (*it)->Pt();
    if(isMinBias) fHistMgr->FillTH1("hTracksMB", pt);
    if(isEJ1) fHistMgr->FillTH1("hTracksEJ1", pt);
    if(isEJ2) fHistMgr->FillTH1("hTracksEJ2", pt);
    if(isEG1) fHistMgr->FillTH1("hTracksEG1", pt);
    if(isEG2) fHistMgr->FillTH1("hTracksEG2", pt);
  }

  // Look at the MC Jet container
  AliVParticle *jettrack = NULL;
  GetJetContainer(fJetContainerNameRec.Data())->ResetCurrentID();
  AliEmcalJet *recjet = GetJetContainer(fJetContainerNameRec.Data())->GetNextAcceptJet();
  do{
    for(int itrk = 0; itrk < recjet->GetNumberOfTracks(); itrk++){
      jettrack = recjet->TrackAt(itrk, GetParticleContainer(1)->GetArray());
      if(!AcceptParticle(jettrack)) continue;
      if(!AcceptTrack(dynamic_cast<AliVTrack *>(jettrack))) continue;
      if(isMinBias) fHistMgr->FillTH2("hRecjetTrackMB", recjet->Pt(), jettrack->Pt());
      if(isEJ1) fHistMgr->FillTH2("hRecjetTrackEJ1", recjet->Pt(), jettrack->Pt());
      if(isEJ2) fHistMgr->FillTH2("hRecjetTrackEJ2", recjet->Pt(), jettrack->Pt());
      if(isEG1) fHistMgr->FillTH2("hRecjetTrackEG1", recjet->Pt(), jettrack->Pt());
      if(isEG2) fHistMgr->FillTH2("hRecjetTrackEG2", recjet->Pt(), jettrack->Pt());
    }
  } while ((recjet = GetJetContainer(fJetContainerNameRec.Data())->GetNextAcceptJet()));

  return true;
}

std::vector<const AliVParticle *> AliAnalysisTaskParticleInJet::GetSelectedParticles(AliParticleContainer *const cont) const {
  std::vector<const AliVParticle *> result;
  cont->ResetCurrentID();
  AliVParticle * test = cont->GetNextAcceptParticle();
  AliVTrack *track = NULL;
  do{
    if(!AcceptParticle(test)) continue;
    if((track = dynamic_cast<AliVTrack *>(test))){
      // apply extra track cuts
      if(!AcceptTrack(track)) continue;
    }
    result.push_back(test);
  } while ((test = cont->GetNextAcceptParticle()));

  return result;
}

Bool_t AliAnalysisTaskParticleInJet::AcceptParticle(const AliVParticle * const part) const{
  if(part->Phi() < 1.4 || part->Phi() > 3.1) return false;
  if(TMath::Abs(part->Eta()) > 0.5) return false;
  if(!part->Charge()) return false;
  return true;
}

Bool_t AliAnalysisTaskParticleInJet::AcceptTrack(AliVTrack * const track) const{
   if(!fTrackSelection->IsTrackAccepted(track)) return false;
   return true;
}

Bool_t AliAnalysisTaskParticleInJet::IsPhysicalPrimary(const AliVParticle *const part) const {
  const AliMCParticle *mcpart = dynamic_cast<const AliMCParticle *>(part);
  if(mcpart){
    return MCEvent()->IsPhysicalPrimary(part->GetLabel());
  } else {
    const AliAODMCParticle *aodpart = dynamic_cast<const AliAODMCParticle *>(part);
    if(aodpart) return aodpart->IsPhysicalPrimary();
  }
  return false;
}

void AliAnalysisTaskParticleInJet::CreatePtBinning(TArrayD& binning) const {
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

void AliAnalysisTaskParticleInJet::CreateLinearBinning(TArrayD& binning, int nbins, double min, double max) const {
  double binwidth = (max-min)/static_cast<double>(nbins);
  binning.Set(nbins+1);
  binning[0] = min;
  double currentlimit = min + binwidth;
  for(int ibin = 0; ibin < nbins; ibin++){
    binning[ibin+1] = currentlimit;
    currentlimit += binwidth;
  }
}
