////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticleYPt - the cut monitor for particles to study    //
// the difference between reconstructed and true momentum                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorParticleYPt::AliFemtoCutMonitorParticleYPt():
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
  fYPt = new TH2D("YPt", "Rapidity vs Pt",              140, -1.4, 1.4, 100, 0.0, 5.0);
  fYPhi = new TH2D("YPhi", "Rapidity vs Phi",           140, -1.4, 1.4, 100, -TMath::Pi(), TMath::Pi());
  fPtPhi = new TH2D("PtPhi", "Pt vs Phi",               100,  0.0, 5.0, 100, -TMath::Pi(), TMath::Pi());
  fEtaPhi = new TH2D("EtaPhi", "Pseudorapidity vs Phi", 140, -1.4, 1.4, 100, -TMath::Pi(), TMath::Pi());
  fEtaPt = new TH2D("EtaPt", "Pseudorapidity vs Pt",    140, -1.4, 1.4, 100, 0.0, 5.0);
  // fEtaPhiW = new TH2D("EtaPhiW", "Pseudorapidity vs Phi chi2/N weighted", 140, -1.4, 1.4, 100, -TMath::Pi(), TMath::Pi());
  // fEtaPtW = new TH2D("EtaPtW", "Pseudorapidity vs Pt chi2/N weighted",    140, -1.4, 1.4, 100, 0.0, 5.0)
  ;
  fDCARPt = new TH2D("DCARPt", "DCA in XY vs. Pt", 400, -3.0, 3.0, 100,0.0,5.0);
  fDCAZPt = new TH2D("DCAZPt", "DCA in Z vs. Pt", 400, -3.0, 3.0, 100,0.0,5.0);
}

AliFemtoCutMonitorParticleYPt::AliFemtoCutMonitorParticleYPt(const char *aName, float aMass):
  AliFemtoCutMonitor(),
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
    fYPt = new TH2D(name, "Rapdity vs Pt", 140, -1.4, 1.4, 100, 0.0, 5.0);
  snprintf(name, 200, "YPhi%s", aName);
  fYPhi = new TH2D(name, "Rapidity vs Phi",           140, -1.4, 1.4, 100, -TMath::Pi(), TMath::Pi());
  snprintf(name, 200, "PtPhi%s", aName);
    fPtPhi = new TH2D(name, "Pt vs Phi",               100,  0.0, 5.0, 100, -TMath::Pi(), TMath::Pi());
  snprintf(name, 200, "EtaPhi%s", aName);
  fEtaPhi = new TH2D(name, "Pseudorapidity vs Phi", 140, -1.4, 1.4, 100, -TMath::Pi(), TMath::Pi());
  snprintf(name, 200, "EtaPt%s", aName);
    fEtaPt = new TH2D(name, "Pseudorapidity vs Pt",    140, -1.4, 1.4, 100, 0.0, 5.0);
  // snprintf(name, 200, "EtaPhiW%s", aName);
  // fEtaPhiW = new TH2D(name, "Pseudorapidity vs Phi chi2/N weighted", 140, -1.4, 1.4, 100, -TMath::Pi(), TMath::Pi());
  // snprintf(name, 200, "EtaPtW%s", aName);
  // fEtaPtW = new TH2D(name, "Pseudorapidity vs Pt chi2/N weighted",    140, -1.4, 1.4, 100, 0.0, 2.0);
  snprintf(name, 200, "DCARPt%s", aName);
    fDCARPt = new TH2D(name, "DCA in XY vs. Pt", 400, -3.0, 3.0, 100,0.0,3.0);
  snprintf(name, 200, "DCAZPt%s", aName);
    fDCAZPt = new TH2D(name, "DCA in Z vs. Pt", 400, -3.0, 3.0, 100,0.0,3.0);
}

AliFemtoCutMonitorParticleYPt::AliFemtoCutMonitorParticleYPt(const AliFemtoCutMonitorParticleYPt &aCut):
  AliFemtoCutMonitor(),
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
  fMass = aCut.fMass; 
}

AliFemtoCutMonitorParticleYPt::~AliFemtoCutMonitorParticleYPt()
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
}

AliFemtoCutMonitorParticleYPt& AliFemtoCutMonitorParticleYPt::operator=(const AliFemtoCutMonitorParticleYPt& aCut)
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
  
  return *this;
}

AliFemtoString AliFemtoCutMonitorParticleYPt::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorParticleYPt report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorParticleYPt::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track
  float tEnergy = ::sqrt(aTrack->P().Mag2()+fMass*fMass);
  if(tEnergy==abs(aTrack->P().z())) tEnergy+=0.001;
  float tRapidity = 0.5*::log((tEnergy+aTrack->P().z())/(tEnergy-aTrack->P().z()));
  float tPt = ::sqrt((aTrack->P().x())*(aTrack->P().x())+(aTrack->P().y())*(aTrack->P().y()));
  
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

  fYPt->Fill(tRapidity, tPt);
  fYPhi->Fill(tRapidity, tPhi);
  fPtPhi->Fill(tPt, tPhi);
  fEtaPhi->Fill(tEta, tPhi);
  fEtaPt->Fill(tEta, tPt);
  // fEtaPhiW->Fill(tEta, tPhi, chi2w);
  // fEtaPtW->Fill(tEta, tPt, chi2w);
  fDCARPt->Fill(dcar, tPt);
  fDCAZPt->Fill(dcaz, tPt);
}

void AliFemtoCutMonitorParticleYPt::Write()
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

TList *AliFemtoCutMonitorParticleYPt::GetOutputList()
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
