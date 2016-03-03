/*
 * AliAnalysisTaskTrackDensity.cxx
 *
 *  Created on: Mar 2, 2016
 *      Author: markus
 */
#include <THashList.h>
#include <THistManager.h>
#include <TMath.h>

#include "AliEmcalJet.h"

#include "AliAnalysisTaskTrackDensity.h"

ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskTrackDensity)

namespace EMCalTriggerPtAnalysis {

AliAnalysisTaskTrackDensity::AliAnalysisTaskTrackDensity() :
    AliAnalysisTaskEmcalJet(),
    fHistos(NULL),
    fMCJetContainerName(""),
    fMCParticleContainerName(""),
    fJetRadii(),
    fJetPtBins(),
    fPtMinSteps(),
    fParticlePtBinning()
{
}

AliAnalysisTaskTrackDensity::AliAnalysisTaskTrackDensity(const char *name) :
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fHistos(NULL),
    fMCJetContainerName(""),
    fMCParticleContainerName(""),
    fJetRadii(),
    fJetPtBins(),
    fPtMinSteps(),
    fParticlePtBinning()
{
  double  defaultJetRadii[] = {0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4},
          defaultPtMinSteps[] = {0.5, 1, 2, 5, 10, 20},
          defaultJetPtBins[] = {20, 40, 60, 80, 100, 150, 200, 1000},
          defaultParticlePtBinning[] = {0., 0.5, 1., 1.5, 2., 3., 4., 5., 7.5, 10., 15., 20., 30., 40., 60., 80., 100.};
  fJetRadii.Set(sizeof(defaultJetRadii)/sizeof(double), defaultJetRadii);
  fJetPtBins.Set(sizeof(defaultJetPtBins)/sizeof(double), defaultJetPtBins);
  fPtMinSteps.Set(sizeof(defaultPtMinSteps)/sizeof(double), defaultPtMinSteps);
  fParticlePtBinning.Set(sizeof(defaultParticlePtBinning)/sizeof(double), defaultParticlePtBinning);
}

AliAnalysisTaskTrackDensity::~AliAnalysisTaskTrackDensity() {
}

void AliAnalysisTaskTrackDensity::UserCreateOutputObjects() {
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  fHistos = new THistManager("trackdensityhistos");

  TArrayD linearBinning(102);
  int counter = 0;
  for(double val = -0.5; val <= 100.5; val += 1.){
    linearBinning[counter++] = val;
  }

  for(int irad  = 0; irad < fJetRadii.GetSize()-1; irad++){
    for(int ptstep = 0; ptstep <fPtMinSteps.GetSize(); ptstep++){
      fHistos->CreateTH2(Form("trackDensityJet_r%d_%d_minpt%d",
                          static_cast<int>(fJetRadii[irad] * 100.),
                          static_cast<int>(fJetRadii[irad+1] * 100.),
                          static_cast<int>(fPtMinSteps[ptstep] * 10.)),
                         Form("Density of tracks with p_{t} > %f GeV/c in r [%.2f, %.2f]; p_{t, jet} (GeV/c); Number of tracks",
                             fPtMinSteps[ptstep],
                             fJetRadii[irad],
                             fJetRadii[irad+1]),
                         200, 0., 200.,
                         102, -0.5, 100.5);
    }
    for(int jetptbin = 0 ; jetptbin < fJetPtBins.GetSize()-1; jetptbin++)
      fHistos->CreateTH2(Form("trackDensityParticle_r%d_%d_jetpt%d_%d",
                          static_cast<int>(fJetRadii[irad] * 100.),
                          static_cast<int>(fJetRadii[irad+1] * 100.),
                          static_cast<int>(fJetPtBins[jetptbin]),
                          static_cast<int>(fJetPtBins[jetptbin+1])),
                         Form("Density of tracks in jet with p_{t} [%.1f, %.1f] in r[%.2f,%2f]",
                          fJetPtBins[jetptbin],
                          fJetPtBins[jetptbin+1],
                          fJetRadii[irad],
                          fJetRadii[irad+1]),
                         fParticlePtBinning, linearBinning);
  }

  for(TIter histiter = TIter(fHistos->GetListOfHistograms()).Begin(); histiter != TIter::End(); ++histiter){
    fOutput->Add(*histiter);
  }
  PostData(1, fOutput);
}

bool AliAnalysisTaskTrackDensity::Run(){
  AliEmcalJet *myjet = NULL;
  AliJetContainer *mcjetcontainer = GetJetContainer(fMCJetContainerName.Data());
  AliParticleContainer *mcparticleContainer = GetParticleContainer(fMCParticleContainerName.Data());
  if(!mcparticleContainer) printf("Error getting particle container with name %s\n", fMCParticleContainerName.Data());
  double jetptmin, jetptmax;
  mcjetcontainer->ResetCurrentID(-1);
  while((myjet = mcjetcontainer->GetNextAcceptJet())){
    for(int irad = 0 ; irad < fJetRadii.GetSize()-1; irad++){
      for(int ptstep = 0; ptstep < fPtMinSteps.GetSize(); ptstep++){
        fHistos->FillTH2(Form("trackDensityJet_r%d_%d_minpt%d", static_cast<int>(fJetRadii[irad] * 100.), static_cast<int>(fJetRadii[irad+1] * 100.), static_cast<int>(fPtMinSteps[ptstep] * 10.)),
            TMath::Abs(myjet->Pt()), GetParticleMultiplicity(*myjet, *mcparticleContainer, fPtMinSteps[ptstep], 10000., fJetRadii[irad], fJetRadii[irad+1]));
      }
      FindJetPtBin(myjet, jetptmin, jetptmax);
      if(jetptmin > 0 && jetptmax > 0){
        for(int ptstep =  0; ptstep < fParticlePtBinning.GetSize(); ptstep++){
          double mean = (fParticlePtBinning[ptstep], fParticlePtBinning[+1])/2.;
          fHistos->FillTH2(Form("trackDensityParticle_r%d_%d_jetpt%d_%d", static_cast<int>(fJetRadii[irad] * 100.),static_cast<int>(fJetRadii[irad+1] * 100.),
              static_cast<int>(jetptmin), static_cast<int>(jetptmax)), mean, GetParticleMultiplicity(*myjet, *mcparticleContainer, fParticlePtBinning[ptstep], fParticlePtBinning[+1], fJetRadii[irad], fJetRadii[irad+1]));
        }
      }
    }
  }

  return true;
}

int AliAnalysisTaskTrackDensity::GetParticleMultiplicity(const AliEmcalJet &jet, const AliParticleContainer &partcont, double ptmin, double ptmax, double rmin, double rmax) const {
  AliDebug(1, Form("Next jet: %s\n", jet.toString().Data()));
  TLorentzVector jetaxis(jet.Px(), jet.Py(), jet.Pz(), jet.E());
  int nselected = 0;
  const AliVParticle *jetparticle(NULL);
  for(int ipart = 0; ipart < jet.GetNumberOfTracks(); ipart++){
    jetparticle = jet.TrackAt(ipart, partcont.GetArray());
    if(!TMath::Abs(jetparticle->Charge())) continue;      // select charged particles only
    double partpt = TMath::Abs(jetparticle->Pt());
    if(partpt >= ptmin && partpt < ptmax){
      TLorentzVector partvector(jetparticle->Px(), jetparticle->Py(), jetparticle->Pz(), jetparticle->E());
      double r = TMath::Abs(jetaxis.DeltaR(partvector));
      if(r >= rmin && r < rmax)  nselected++;
    }
  }
  return nselected;
}

void AliAnalysisTaskTrackDensity::FindJetPtBin(const AliEmcalJet *const jet, double &ptmin, double &ptmax) const {
  ptmin = ptmax = -1;
  double jetpt = TMath::Abs(jet->Pt());
  for(int ptstep = 0; ptstep < fJetPtBins.GetSize() - 1; ptstep++){
    if(jetpt >= fJetPtBins[ptstep] && jetpt < fJetPtBins[ptstep+1]){
      ptmin = fJetPtBins[ptstep];
      ptmax = fJetPtBins[ptstep+1];
      break;
    }
  }
}

} /* namespace EMCalTriggerPtAnalysis */
