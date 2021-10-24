/*
 * AliAnalysisTaskTrackDensityData.cxx
 *
 *  Created on: Mar 11, 2016
 *      Author: markus
 */
#include <THistManager.h>
#include <TMath.h>

#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliTrackContainer.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerBinningFactory.h"

#include "AliAnalysisTaskTrackDensityData.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskTrackDensityData)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskTrackDensityData::AliAnalysisTaskTrackDensityData() :
    AliAnalysisTaskEmcalJet(),
    fHistos(NULL),
    fTrackSelection(NULL),
    fBinHandler(NULL),
    fNameJetContainer(""),
    fNameTrackContainer("")
{
}

AliAnalysisTaskTrackDensityData::AliAnalysisTaskTrackDensityData(const char *name) :
    AliAnalysisTaskEmcalJet(name, true),
    fHistos(NULL),
    fTrackSelection(NULL),
    fBinHandler(new AliEMCalTriggerBinningComponent),
    fNameJetContainer(""),
    fNameTrackContainer("")
{
  SetMakeGeneralHistograms(true);
  double  defaultJetRadii[] = {0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4},
          defaultPtMinSteps[] = {0.5, 1, 2, 5, 10, 20},
          defaultJetPtBins[] = {20, 40, 60, 80, 100, 150, 200, 1000};

  fBinHandler->SetLinearBinning("jetpt", 100, 0, 100.);
  fBinHandler->SetLinearBinning("contributors", 101, -0.5, 100.5);
  fBinHandler->SetBinning("jetradii", sizeof(defaultJetRadii)/sizeof(double) - 1, defaultJetRadii);
  fBinHandler->SetBinning("ptmin", sizeof(defaultPtMinSteps)/sizeof(double) - 1, defaultPtMinSteps);
  fBinHandler->SetBinning("jetlarge", sizeof(defaultJetPtBins)/sizeof(double) - 1, defaultJetPtBins);
}


AliAnalysisTaskTrackDensityData::~AliAnalysisTaskTrackDensityData() {
  if(fHistos) delete fHistos;
  if(fBinHandler) delete fBinHandler;
}

void AliAnalysisTaskTrackDensityData::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  fHistos = new THistManager("histos");
  fHistos->ReleaseOwner();

  AliEMCalTriggerBinningFactory bininit;
  bininit.Create(fBinHandler);

  const TBinning *trackptbinning = fBinHandler->GetBinning("pt"),
          *jetptbinning = fBinHandler->GetBinning("jetpt"),
          *contributorbinning = fBinHandler->GetBinning("contributors"),
          *jetradii = fBinHandler->GetBinning("jetradii"),
          *ptminsteps = fBinHandler->GetBinning("ptmin"),
          *jetptlarge = fBinHandler->GetBinning("jetlarge");

  fHistos->CreateTH1("hTrackPtSel", "Pt spectrum of selected tracks", *trackptbinning);
  fHistos->CreateTH1("hTrackPtSelEvent", "Pt spectrum of selected tracks (directly from the input event)", *trackptbinning);
  fHistos->CreateTH2("hJetMultiplicity", "Multiplicity of particles in jets", *jetptbinning, *contributorbinning);
  fHistos->CreateTH2("hParticlePtJet", "Correlation between track pt and jet pt", *jetptbinning, *trackptbinning);

  TArrayD jrad, ptms, jptl;
  jetradii->CreateBinEdges(jrad);
  ptminsteps->CreateBinEdges(ptms);
  jetptlarge->CreateBinEdges(jptl);
  for(int irad  = 0; irad < jrad.GetSize()-1; irad++){
    for(int ptstep = 0; ptstep <ptms.GetSize(); ptstep++){
      fHistos->CreateTH2(Form("trackDensityJet_r%d_%d_minpt%d",
                              static_cast<int>(jrad[irad] * 100.),
                              static_cast<int>(jrad[irad+1] * 100.),
                              static_cast<int>(ptms[ptstep] * 10.)),
                         Form("Density of tracks with p_{t} > %f GeV/c in r [%.2f, %.2f]; p_{t, jet} (GeV/c); Number of tracks",
                             ptms[ptstep],
                             jrad[irad],
                             jrad[irad+1]),
                         200, 0., 200.,
                         102, -0.5, 100.5);
    }
    for(int jetptbin = 0 ; jetptbin < jptl.GetSize()-1; jetptbin++){
      fHistos->CreateTH2(Form("trackDensityParticle_r%d_%d_jetpt%d_%d",
                              static_cast<int>(jrad[irad] * 100.),
                              static_cast<int>(jrad[irad+1] * 100.),
                              static_cast<int>(jptl[jetptbin]),
                              static_cast<int>(jptl[jetptbin+1])),
                         Form("Density of tracks in jet with p_{t} [%.1f, %.1f] in r[%.2f,%2f]",
                             jptl[jetptbin],
                             jptl[jetptbin+1],
                             jrad[irad],
                             jrad[irad+1]),
                         *trackptbinning, *contributorbinning);
    }
  }

  for(THistManager::iterator it = fHistos->begin(); it != fHistos->end(); it++){
    fOutput->Add(*it);
  }

  /*
  for(auto it : *fHistos){
    fOutput->Add(it);
  }
  */
}

bool AliAnalysisTaskTrackDensityData::Run(){
  // Loop over jets
  AliJetContainer *jcont = GetJetContainer(fNameJetContainer.Data());
  AliTrackContainer *tcont = GetTrackContainer(fNameTrackContainer.Data());
  AliEmcalJet *myjet = NULL;
  AliVParticle *jetparticle = NULL;
  TArrayD particlePtBinning, jetradii, ptminsteps;
  fBinHandler->GetBinning("pt")->CreateBinEdges(particlePtBinning);
  fBinHandler->GetBinning("jetradii")->CreateBinEdges(jetradii);
  fBinHandler->GetBinning("ptmin")->CreateBinEdges(ptminsteps);
  for(int ipart = 0; ipart < InputEvent()->GetNumberOfTracks(); ipart++){
    jetparticle = InputEvent()->GetTrack(ipart);
    if(TMath::Abs(jetparticle->Eta()) > 0.8) continue;
    if(!fTrackSelection->IsTrackAccepted(static_cast<AliVTrack *>(jetparticle))) continue;
    fHistos->FillTH1("hTrackPtSelEvent", TMath::Abs(jetparticle->Pt()));
  }
  const AliTrackIterableContainer accepted_tracks = tcont->accepted();
  for(AliTrackIterableContainer::iterator trackiter = accepted_tracks.begin(); trackiter != accepted_tracks.end(); ++trackiter){
    jetparticle = *trackiter;
    if(TMath::Abs(jetparticle->Eta()) > 0.8) continue;
    if(!fTrackSelection->IsTrackAccepted(static_cast<AliVTrack *>(jetparticle))) continue;
    fHistos->FillTH1("hTrackPtSel", TMath::Abs(jetparticle->Pt()));
  }
  const AliJetIterableContainer accepted_jets = jcont->accepted();
  for(AliJetIterableContainer::iterator jetiter = accepted_jets.begin(); jetiter != accepted_jets.end(); ++jetiter){
    myjet = *jetiter;
    fHistos->FillTH2("hJetMultiplicity", myjet->Pt(), myjet->GetNumberOfConstituents());

    for(int iconst = 0; iconst < myjet->GetNumberOfTracks(); iconst++){
      jetparticle = myjet->TrackAt(iconst, tcont->GetArray());
      if(TMath::Abs(jetparticle->Eta()) > 0.8) continue;
      if(!fTrackSelection->IsTrackAccepted(static_cast<AliVTrack *>(jetparticle))) continue;
      fHistos->FillTH2("hParticlePtJet", myjet->Pt(), jetparticle->Pt());
    }

    for(int irad = 0 ; irad < jetradii.GetSize()-1; irad++){
      for(int ptstep = 0; ptstep < ptminsteps.GetSize(); ptstep++){
        fHistos->FillTH2(Form("trackDensityJet_r%d_%d_minpt%d", static_cast<int>(jetradii[irad] * 100.), static_cast<int>(jetradii[irad+1] * 100.), static_cast<int>(ptminsteps[ptstep] * 10.)),
            TMath::Abs(myjet->Pt()), GetParticleMultiplicity(*myjet, *tcont, ptminsteps[ptstep], 10000., jetradii[irad], jetradii[irad+1]));
      }
      double jetptmin, jetptmax;
      FindJetPtBin(myjet, jetptmin, jetptmax);
      if(jetptmin > 0 && jetptmax > 0){
        for(int ptstep =  0; ptstep < particlePtBinning.GetSize()-1; ptstep++){
          double mean = (particlePtBinning[ptstep] + particlePtBinning[ptstep+1])/2.;
          fHistos->FillTH2(Form("trackDensityParticle_r%d_%d_jetpt%d_%d", static_cast<int>(jetradii[irad] * 100.),static_cast<int>(jetradii[irad+1] * 100.),
              static_cast<int>(jetptmin), static_cast<int>(jetptmax)), mean, GetParticleMultiplicity(*myjet, *tcont, particlePtBinning[ptstep], particlePtBinning[+1], jetradii[irad], jetradii[irad+1]));
        }
      }
    }

  }
  return true;
}

int AliAnalysisTaskTrackDensityData::GetParticleMultiplicity(const AliEmcalJet &jet, const AliParticleContainer &partcont, double ptmin, double ptmax, double rmin, double rmax) const {
  AliDebug(1, Form("Next jet: %s\n", jet.toString().Data()));
  TLorentzVector jetaxis(jet.Px(), jet.Py(), jet.Pz(), jet.E());
  int nselected = 0;
  AliVParticle *jetparticle(NULL);
  for(int ipart = 0; ipart < jet.GetNumberOfTracks(); ipart++){
    jetparticle = static_cast<AliVParticle *>(jet.TrackAt(ipart, partcont.GetArray()));
    if(TMath::Abs(jetparticle->Eta()) > 0.8) continue;
    if(!fTrackSelection->IsTrackAccepted(static_cast<AliVTrack *>(jetparticle))) continue;
    double partpt = TMath::Abs(jetparticle->Pt());
    if(partpt >= ptmin && partpt < ptmax){
      TLorentzVector partvector(jetparticle->Px(), jetparticle->Py(), jetparticle->Pz(), jetparticle->E());
      double r = TMath::Abs(jetaxis.DeltaR(partvector));
      if(r >= rmin && r < rmax)  nselected++;
    }
  }
  return nselected;
}

void AliAnalysisTaskTrackDensityData::FindJetPtBin(const AliEmcalJet *const jet, double &ptmin, double &ptmax) const {
  TArrayD jetptlarge;
  const TBinning *jptlbin = fBinHandler->GetBinning("jetlarge");
  jptlbin->CreateBinEdges(jetptlarge);
  ptmin = ptmax = -1;
  double jetpt = TMath::Abs(jet->Pt());
  for(int ptstep = 0; ptstep < jetptlarge.GetSize() - 1; ptstep++){
    if(jetpt >= jetptlarge[ptstep] && jetpt < jetptlarge[ptstep+1]){
      ptmin = jetptlarge[ptstep];
      ptmax = jetptlarge[ptstep+1];
      break;
    }
  }
}
