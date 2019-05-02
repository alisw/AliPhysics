////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticleYPtWithWeights - the cut monitor for particles to study    //
// the difference between reconstructed and true momentum                     //
// author: Rafal Maselek rmaselek@cern.ch                                     //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticleYPtWithWeights.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorParticleYPtWithWeights::AliFemtoCutMonitorParticleYPtWithWeights(TH2D *filter, int calculateWeights = 1):
  filterHist(filter),
  fCalcWeights(calculateWeights),
  fYPt(0),
  fYPhi(0),
  fPtPhi(0),
  fEtaPhi(0),
  fEtaPt(0),
  fEtaPhiW(0),
  fEtaPtW(0),
  fDCARPt(0),
  fDCAZPt(0),
  fMass(0.13957)
{
  // Default constructor
  fYPt = new TH2D("YPt", "Rapidity vs Pt",              100, 0.0, 6.0, 100, 0.0, 3.0);
  fYPhi = new TH2D("YPhi", "Rapidity vs Phi",           100, 0.0, 6.0, 100, -TMath::Pi(), TMath::Pi());
  fPtPhi = new TH2D("PtPhi", "Pt vs Phi",               100,  0.0, 3.0, 100, -TMath::Pi(), TMath::Pi());
  fEtaPhi = new TH2D("EtaPhi", "Pseudorapidity vs Phi", 100, 0.0, 6.0, 100, -TMath::Pi(), TMath::Pi());
  fEtaPt = new TH2D("EtaPt", "Pseudorapidity vs Pt",    100, 0.0, 6.0, 100, 0.0, 3.0);
  // fEtaPhiW = new TH2D("EtaPhiW", "Pseudorapidity vs Phi chi2/N weighted", 140, -1.4, 1.4, 100, -TMath::Pi(), TMath::Pi());
  // fEtaPtW = new TH2D("EtaPtW", "Pseudorapidity vs Pt chi2/N weighted",    140, -1.4, 1.4, 100, 0.0, 5.0);
  fDCARPt = new TH2D("DCARPt", "DCA in XY vs. Pt", 400, -3.0, 3.0, 100,0.0,3.0);
  fDCAZPt = new TH2D("DCAZPt", "DCA in Z vs. Pt", 400, -3.0, 3.0, 100,0.0,3.0);


}

AliFemtoCutMonitorParticleYPtWithWeights::AliFemtoCutMonitorParticleYPtWithWeights(const char *aName, float aMass, TH2D *filter, int calculateWeights=1):
  filterHist(filter),
  fCalcWeights(calculateWeights),
  fYPt(0),
  fYPhi(0),
  fPtPhi(0),
  fEtaPhi(0),
  fEtaPt(0),
  fEtaPhiW(0),
  fEtaPtW(0),
  fDCARPt(0),
  fDCAZPt(0),
  fMass(aMass)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "YPt%s", aName);
  fYPt = new TH2D(name, "Rapdity vs Pt", 100, 0.0, 6.0, 100, 0.0, 3.0);
  snprintf(name, 200, "YPhi%s", aName);
  fYPhi = new TH2D(name, "Rapidity vs Phi",           100, 0.0, 6.0, 100, -TMath::Pi(), TMath::Pi());
  snprintf(name, 200, "PtPhi%s", aName);
  fPtPhi = new TH2D(name, "Pt vs Phi",               100,  0.0, 3.0, 100, -TMath::Pi(), TMath::Pi());
  snprintf(name, 200, "EtaPhi%s", aName);
  fEtaPhi = new TH2D(name, "Pseudorapidity vs Phi", 100, 0.0, 6.0, 100, -TMath::Pi(), TMath::Pi());
  snprintf(name, 200, "EtaPt%s", aName);
  fEtaPt = new TH2D(name, "Pseudorapidity vs Pt",    100, 0.0, 6.0, 100, 0.0, 3.0);
  // snprintf(name, 200, "EtaPhiW%s", aName);
  // fEtaPhiW = new TH2D(name, "Pseudorapidity vs Phi chi2/N weighted", 140, -1.4, 1.4, 100, -TMath::Pi(), TMath::Pi());
  // snprintf(name, 200, "EtaPtW%s", aName);
  // fEtaPtW = new TH2D(name, "Pseudorapidity vs Pt chi2/N weighted",    140, -1.4, 1.4, 100, 0.0, 2.0);
  snprintf(name, 200, "DCARPt%s", aName);
  fDCARPt = new TH2D(name, "DCA in XY vs. Pt", 400, -3.0, 3.0, 100,0.0,3.0);
  snprintf(name, 200, "DCAZPt%s", aName);
  fDCAZPt = new TH2D(name, "DCA in Z vs. Pt", 400, -3.0, 3.0, 100,0.0,3.0);
}

AliFemtoCutMonitorParticleYPtWithWeights::AliFemtoCutMonitorParticleYPtWithWeights(const AliFemtoCutMonitorParticleYPtWithWeights &aCut):
  filterHist(0),
  fCalcWeights(aCut.fCalcWeights),
  fYPt(0),
  fYPhi(0),
  fPtPhi(0),
  fEtaPhi(0),
  fEtaPt(0),
  fEtaPhiW(0),
  fEtaPtW(0),
  fDCARPt(0),
  fDCAZPt(0),
  fMass(aCut.fMass)
{
  // copy constructor
  fYPt = new TH2D(*aCut.fYPt);
  fYPhi = new TH2D(*aCut.fYPhi);
  fPtPhi = new TH2D(*aCut.fPtPhi);
  fEtaPhi = new TH2D(*aCut.fEtaPhi);
  fEtaPt = new TH2D(*aCut.fEtaPt);
  // fEtaPhiW = new TH2D(*aCut.fEtaPhiW);
  // fEtaPtW = new TH2D(*aCut.fEtaPtW);
  fDCARPt = new TH2D(*aCut.fDCARPt);
  fDCAZPt = new TH2D(*aCut.fDCAZPt);
  filterHist = new TH2D(*aCut.filterHist);
}

AliFemtoCutMonitorParticleYPtWithWeights::~AliFemtoCutMonitorParticleYPtWithWeights()
{
  // Destructor
  delete fYPt;
  delete fYPhi;
  delete fPtPhi;
  delete fEtaPhi;
  delete fEtaPt;
  // delete fEtaPhiW;
  // delete fEtaPtW;
  delete fDCARPt;
  delete fDCAZPt;
  if(filterHist) delete filterHist;
}

AliFemtoCutMonitorParticleYPtWithWeights& AliFemtoCutMonitorParticleYPtWithWeights::operator=(const AliFemtoCutMonitorParticleYPtWithWeights& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fYPt) delete fYPt;
  fYPt = new TH2D(*aCut.fYPt);
  if (fYPhi) delete fYPhi;
  fYPhi = new TH2D(*aCut.fYPhi);
  if (fPtPhi) delete fPtPhi;
  fPtPhi = new TH2D(*aCut.fPtPhi);
  if (fEtaPhi) delete fEtaPhi;
  fEtaPhi = new TH2D(*aCut.fEtaPhi);
  if (fEtaPt) delete fEtaPt;
  fEtaPt = new TH2D(*aCut.fEtaPt);
  // if (fEtaPhiW) delete fEtaPhiW;
  // fEtaPhiW = new TH2D(*aCut.fEtaPhiW);
  // if (fEtaPtW) delete fEtaPtW;
  // fEtaPtW = new TH2D(*aCut.fEtaPtW);
  if (fDCARPt) delete fDCARPt;
  fDCARPt = new TH2D(*aCut.fDCARPt);
  if (fDCAZPt) delete fDCAZPt;
  fDCAZPt = new TH2D(*aCut.fDCAZPt);

  if(filterHist) delete filterHist;
  filterHist = new TH2D(*aCut.filterHist);

  fCalcWeights = aCut.fCalcWeights;
  
  return *this;
}

AliFemtoString AliFemtoCutMonitorParticleYPtWithWeights::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorParticleYPtWithWeights report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorParticleYPtWithWeights::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track
  float tEnergy = ::sqrt(aTrack->P().Mag2()+fMass*fMass);
  if(tEnergy==abs(aTrack->P().z())) tEnergy+=0.001;

  float tRapidity = -1000;
  if((tEnergy+aTrack->P().z())/(tEnergy-aTrack->P().z())>0)
    tRapidity = 0.5*::log((tEnergy+aTrack->P().z())/(tEnergy-aTrack->P().z()));
  else
    tRapidity = -999;
  
  float tPt = -1000;
  if((aTrack->P().x())*(aTrack->P().x())+(aTrack->P().y())*(aTrack->P().y())>=0)
    tPt = ::sqrt((aTrack->P().x())*(aTrack->P().x())+(aTrack->P().y())*(aTrack->P().y()));
  else
    tPt = -999;

  double weight = filterHist->GetBinContent(filterHist->FindBin(tRapidity, tPt));
  if(fCalcWeights == 0)
    weight=1.0;

  float tEta;
  if(aTrack->P().Theta()==0)
    tEta=0;
  else
    tEta = -TMath::Log(TMath::Tan(aTrack->P().Theta()/2.0));
  float tPhi = aTrack->P().Phi();
  // float chi2w;
  float dcar = aTrack->ImpactD();
  float dcaz = aTrack->ImpactZ();
  // if (aTrack->TPCncls() > 0)
  //   chi2w = aTrack->TPCchi2()/aTrack->TPCncls();
  // else
  //   chi2w = 6.0;

  //  cout << " CMYPt: " << fYPt << " " << fYPt->GetEntries() << " " << tRapidity << " " << tPt << endl;

  fYPt->Fill(tRapidity, tPt, weight);
  fYPhi->Fill(tRapidity, tPhi, weight);
  fPtPhi->Fill(tPt, tPhi, weight);
  fEtaPhi->Fill(tEta, tPhi, weight);
  fEtaPt->Fill(tEta, tPt, weight);
  // fEtaPhiW->Fill(tEta, tPhi, chi2w);
  // fEtaPtW->Fill(tEta, tPt, chi2w);
  fDCARPt->Fill(dcar, tPt, weight);
  fDCAZPt->Fill(dcaz, tPt, weight);
}

void AliFemtoCutMonitorParticleYPtWithWeights::Write()
{
  // Write out the relevant histograms
  fYPt->Write();
  fYPhi->Write();
  fPtPhi->Write();
  fEtaPhi->Write();
  fEtaPt->Write();
  // fEtaPhiW->Write();
  // fEtaPtW->Write();
  fDCARPt->Write();
  fDCAZPt->Write();
}

TList *AliFemtoCutMonitorParticleYPtWithWeights::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fYPt);
  tOutputList->Add(fYPhi);
  tOutputList->Add(fPtPhi);
  tOutputList->Add(fEtaPhi);
  tOutputList->Add(fEtaPt);
  // tOutputList->Add(fEtaPhiW);
  // tOutputList->Add(fEtaPtW);
  tOutputList->Add(fDCARPt);
  tOutputList->Add(fDCAZPt);

  return tOutputList;
}
