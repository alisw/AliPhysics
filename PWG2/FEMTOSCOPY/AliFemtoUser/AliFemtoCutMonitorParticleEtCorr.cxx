////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticleEtCorr - the cut monitor for particles           //
// which saves particles' et histogram and makes the bin-by-bin correlation   //
//                                                                            //
// Author: Adam.Kisiel@cern.ch                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticleEtCorr.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorParticleEtCorr::AliFemtoCutMonitorParticleEtCorr():
  AliFemtoCutMonitor(),
  fPhiBins(60),
  fPtPerPhi(0),     
  fPtCovPerPhi(0),
  fPtMultPerPhi(0),
  fNEventsProcessed(0)
{
  // Default constructor
}

AliFemtoCutMonitorParticleEtCorr::AliFemtoCutMonitorParticleEtCorr(const char *aName, int aPhiBins):
  AliFemtoCutMonitor(),
  fPhiBins(aPhiBins),
  fPtPerPhi(0),     
  fPtCovPerPhi(0),
  fPtMultPerPhi(0),
  fNEventsProcessed(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "EtCorrAvgPt%s", aName);
  fPtPerPhi = new TH1D(name, "Average Pt Per Phi", aPhiBins, -0.5, aPhiBins-0.5);
  snprintf(name, 200, "EtCorrMult%s", aName);
  fPtMultPerPhi = new TH2D(name, "Multiplicity Per Phi", aPhiBins, -0.5, aPhiBins-0.5, aPhiBins, -0.5, aPhiBins-0.5);
  snprintf(name, 200, "EtCorrAvgPtCov%s", aName);
  fPtCovPerPhi = new TH2D(name, "Covariance of Average Pt Per Phi", aPhiBins, -0.5, aPhiBins-0.5, aPhiBins, -0.5, aPhiBins-0.5);

  fPtPerPhi->Sumw2();
  fPtCovPerPhi->Sumw2();
  fPtMultPerPhi->Sumw2();
  fPhiBins = aPhiBins;
}

AliFemtoCutMonitorParticleEtCorr::AliFemtoCutMonitorParticleEtCorr(const AliFemtoCutMonitorParticleEtCorr &aCut):
  AliFemtoCutMonitor(),
  fPhiBins(0),
  fPtPerPhi(0),     
  fPtCovPerPhi(0),
  fPtMultPerPhi(0),
  fNEventsProcessed(0)
{
  // copy constructor
  if (fPtCovPerPhi) delete fPtCovPerPhi;
  fPtCovPerPhi = new TH2D(*aCut.fPtCovPerPhi);
  if (fPtPerPhi) delete fPtPerPhi;
  fPtPerPhi = new TH1D(*aCut.fPtPerPhi);
  if (fPtMultPerPhi) delete fPtMultPerPhi;
  fPtMultPerPhi = new TH2D(*aCut.fPtMultPerPhi);
  fPhiBins = aCut.fPhiBins;
  fNEventsProcessed = aCut.fNEventsProcessed;
}

AliFemtoCutMonitorParticleEtCorr::~AliFemtoCutMonitorParticleEtCorr()
{
  // Destructor
  delete fPtPerPhi;
  delete fPtMultPerPhi;
  delete fPtCovPerPhi;
}

AliFemtoCutMonitorParticleEtCorr& AliFemtoCutMonitorParticleEtCorr::operator=(const AliFemtoCutMonitorParticleEtCorr& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fPtCovPerPhi) delete fPtCovPerPhi;
  fPtCovPerPhi = new TH2D(*aCut.fPtCovPerPhi);
  if (fPtPerPhi) delete fPtPerPhi;
  fPtPerPhi = new TH1D(*aCut.fPtPerPhi);
  if (fPtMultPerPhi) delete fPtMultPerPhi;
  fPtMultPerPhi = new TH2D(*aCut.fPtMultPerPhi);
  fPhiBins = aCut.fPhiBins;
  fNEventsProcessed = aCut.fNEventsProcessed;
  
  return *this;
}

AliFemtoString AliFemtoCutMonitorParticleEtCorr::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorParticleEtCorr report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorParticleEtCorr::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track
  //  float tEnergy = ::sqrt(aTrack->P().mag2()+fMass*fMass);
  //  float tRapidity = 0.5*::log((tEnergy+aTrack->P().z())/(tEnergy-aTrack->P().z()));
  float tPt = ::sqrt((aTrack->P().x())*(aTrack->P().x())+(aTrack->P().y())*(aTrack->P().y()));
  //  float tEta = -TMath::Log(TMath::Tan(aTrack->P().theta()/2.0));
  float tPhi = aTrack->P().phi();
  Double_t tPiTwo = TMath::Pi()*2;

  while (tPhi > tPiTwo) tPhi -= tPiTwo;
  while (tPhi < 0) tPhi += tPiTwo;

  int nbin = (int) floor(tPhi * fPhiBins / tPiTwo);
  fPtSumEvent[nbin] += tPt;
  fMultSumEvent[nbin] += 1;
}

void AliFemtoCutMonitorParticleEtCorr::Write()
{
  // Write out the relevant histograms
  fPtPerPhi->Write();
  fPtCovPerPhi->Write();
  fPtMultPerPhi->Write();
}

TList *AliFemtoCutMonitorParticleEtCorr::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fPtPerPhi);
  tOutputList->Add(fPtCovPerPhi);
  tOutputList->Add(fPtMultPerPhi);

  return tOutputList;
}

void AliFemtoCutMonitorParticleEtCorr::EventBegin(const AliFemtoEvent* aEvent)
{
  if (aEvent)
    for (int iter=0; iter<fPhiBins; iter++) {
      fPtSumEvent[iter] = 0;
      fMultSumEvent[iter] = 0;
    }
}

void AliFemtoCutMonitorParticleEtCorr::EventEnd(const AliFemtoEvent* aEvent)
{
  if (aEvent) {
    for (int ispt=0; ispt<fPhiBins; ispt++) {
      fPtPerPhi->Fill(ispt, fPtSumEvent[ispt]);
      for (int ispt2=0; ispt2<fPhiBins; ispt2++) {
	fPtCovPerPhi->Fill(ispt, ispt2, fPtSumEvent[ispt]*fPtSumEvent[ispt2]);
	fPtMultPerPhi->Fill(ispt, ispt2, fMultSumEvent[ispt]*fMultSumEvent[ispt2]);
      }
    }
    fNEventsProcessed++;
  }
}

