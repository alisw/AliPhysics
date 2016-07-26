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
          defaultJetPtBins[] = {0,5, 1, 2, 5, 10, 20},
          defaultPtMinSteps[] = {20, 40, 60, 80, 100, 150, 200, 1000},
          defaultParticlePtBinning[] = {0., 0.5, 1., 1.5, 2., 3., 4., 5., 7.5, 10., 15., 20., 30., 40., 60., 80., 100.};
  fJetRadii.Set(sizeof(defaultJetRadii)/sizeof(double), defaultJetRadii);
  fJetPtBins.Set(sizeof(defaultJetPtBins)/sizeof(double), defaultJetPtBins);
  fPtMinSteps.Set(sizeof(defaultPtMinSteps)/sizeof(double), defaultPtMinSteps);
  fParticlePtBinning.Set(sizeof(defaultParticlePtBinning)/sizeof(double), defaultParticlePtBinning);
}

AliAnalysisTaskTrackDensity::~AliAnalysisTaskTrackDensity() {
}

void AliAnalysisTaskTrackDensity::UserCreateOutputObjects() {
  fHistos = new THistManager("trackdensityhistos");

  const double *jetradii = fJetRadii.GetArray(), * ptsteps = fPtMinSteps.GetArray(), *jetptbins = fJetPtBins.GetArray();

  TArrayD linearBinning(101);
  int counter = 0;
  for(double val = -0.5; val <= 100.5; val += 1.) linearBinning[counter++] = val;

  for(const double *raditer = jetradii; raditer < jetradii + sizeof(jetradii)/sizeof(double) - 1; ++raditer){
    for(const double *trackptiter = ptsteps; trackptiter < ptsteps + sizeof(ptsteps)/sizeof(double); ++trackptiter)
      fHistos->CreateTH2(Form("trackDensityJet_r%d_%d_minpt%d",
                          static_cast<int>((*raditer) * 100.),
                          static_cast<int>((*(raditer+1)) * 100.),
                          static_cast<int>((*trackptiter) * 10)),
                         Form("Density of tracks with p_{t} > %f GeV/c in r [%.2f, %.2f]; p_{t, jet} (GeV/c); Number of tracks",
                             *trackptiter,
                             *raditer,
                             *(raditer+1)),
                         200, 0., 200.,
                         101, -0.5, 1.5);
    for(const double *jetptiter = jetptbins; jetptiter < jetptbins + sizeof(jetptbins)/sizeof(double); ++jetptiter)
      fHistos->CreateTH2(Form("trackDensityParticle_r%d_%d_jetpt%d_%d",
                          static_cast<int>((*raditer) * 100.),
                          static_cast<int>((*(raditer+1)) * 100.),
                          static_cast<int>(*jetptiter),
                          static_cast<int>(*(jetptiter+1))),
                         Form("Density of tracks in jet with p_{t} [%.1f, %.1f] in r[%.2f,%2f]",
                          *jetptiter,
                          *(jetptiter+1),
                          *raditer,
                          *(raditer+1)),
                         fParticlePtBinning, linearBinning);
  }

  for(TIter histiter = TIter(fHistos->GetListOfHistograms()).Begin(); histiter != TIter::End(); ++histiter){
    fOutput->Add(*histiter);
  }
}

bool AliAnalysisTaskTrackDensity::Run(){
  AliEmcalJet *myjet = NULL;
  AliJetContainer *mcjetcontainer = GetJetContainer(fMCJetContainerName.Data());
  AliParticleContainer *mcparticleContainer = GetParticleContainer(fMCParticleContainerName.Data());
  const double *jetradii = fJetRadii.GetArray(), *ptsteps = fPtMinSteps.GetArray(), *jetptbins = fJetPtBins.GetArray(), *trackptsteps = fParticlePtBinning.GetArray();
  double jetptmin, jetptmax;
  for(int ijet = 0; ijet < mcjetcontainer->GetNAcceptedJets(); ijet++){
    myjet = mcjetcontainer->GetAcceptJet(ijet);
    for(const double *raditer = jetradii; raditer < jetradii + sizeof(jetradii)/sizeof(double); raditer++){
      for(const double *trackptiter = ptsteps; trackptiter < ptsteps + sizeof(ptsteps)/sizeof(double); trackptiter++){
        fHistos->FillTH2(Form("trackDensityJet_r%d_%d_minpt%d", static_cast<int>((*raditer) * 100.), static_cast<int>((*(raditer+1)) * 100.),
            static_cast<int>((*trackptiter) * 10)),
            TMath::Abs(myjet->Pt()), GetParticleMultiplicity(myjet, mcparticleContainer, *trackptiter, 10000., *raditer, *(raditer+1)));
        FindJetPtBin(myjet, jetptmin, jetptmax);
        if(jetptmin > 0 && jetptmax > 0){
          for(const double *trackptiter = trackptsteps; trackptiter < trackptsteps + sizeof(trackptsteps)/sizeof(double) -1; trackptiter++){
            double mean = (*trackptiter + *(trackptiter+1))/2.;
            fHistos->FillTH2(Form("trackDensityParticle_r%d_%d_jetpt%d_%d", static_cast<int>((*raditer) * 100.),static_cast<int>((*(raditer+1)) * 100.),
                static_cast<int>(jetptmin), static_cast<int>(jetptmax)), mean, GetParticleMultiplicity(myjet, mcparticleContainer, *trackptiter, *(trackptiter+1), *raditer, *(raditer+1)));
          }
        }
      }
    }
  }

  return true;
}

int AliAnalysisTaskTrackDensity::GetParticleMultiplicity(const AliEmcalJet * const jet, AliParticleContainer *partcont, double ptmin, double ptmax, double rmin, double rmax) const {
  TLorentzVector jetaxis(jet->Px(), jet->Py(), jet->Pz(), jet->E());
  int nselected = 0;
  const AliVParticle *jetparticle(NULL);
  for(int ipart = 0; ipart < jet->GetNumberOfTracks(); ipart++){
    jetparticle = jet->TrackAt(ipart, partcont->GetArray());
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
  const double *jetptbins = fJetPtBins.GetArray();
  for(const double *ptiter = jetptbins; ptiter < jetptbins + sizeof(jetptbins)/sizeof(double); ptiter++){
    if(jetpt >= *ptiter && jetpt < *(ptiter + 1)){
      ptmin = *ptiter;
      ptmax = *(ptiter+1);
      break;
    }
  }
}

} /* namespace EMCalTriggerPtAnalysis */
