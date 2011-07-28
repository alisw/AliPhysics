////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticlePID - the cut monitor for particles to study     //
// various aspects of the PID determination                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorParticlePID::AliFemtoCutMonitorParticlePID():
  fTPCdEdx(0),
  fTOFParticle(0),
  fTOFTime(0x0),
  ftofHist(0)
{
  // Default constructor
  fTPCdEdx =  new TH2D("TPCdEdx", "TPC dEdx vs. momentum", 200, 0.1, 4.0, 250, 0.0, 500.0);
  fTOFTime = new TH2D("TOFTime", "TOF Time vs. momentum", 190, 0.1, 2.0, 400, -4000.0, 4000.0);
  ftofHist=new TH2D("TOFHist","TOF momentum vs v",100,0.,1.1,100,0.,3.0);
}

AliFemtoCutMonitorParticlePID::AliFemtoCutMonitorParticlePID(const char *aName, Int_t aTOFParticle):
  AliFemtoCutMonitor(),
  fTPCdEdx(0),
  fTOFParticle(aTOFParticle),
  fTOFTime(0x0),
  ftofHist(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "TPCdEdx%s", aName);
  fTPCdEdx = new TH2D(name, "TPC dEdx vs. momentum", 200, 0.1, 4.0, 250, 0.0, 500.0);

  snprintf(name, 200, "TOFTime%s", aName);
  fTOFTime = new TH2D(name, "TOF Time vs. momentum", 190, 0.1, 2.0, 400, -4000.0, 4000.0);

  snprintf(name, 200, "TOFHist%s", aName);
  ftofHist=new TH2D(name,"TOF momentum vs v",100,0.,1.1,100,0.,3.0);
}

AliFemtoCutMonitorParticlePID::AliFemtoCutMonitorParticlePID(const AliFemtoCutMonitorParticlePID &aCut):
  AliFemtoCutMonitor(),
  fTPCdEdx(0),
  fTOFParticle(0),
  ftofHist(0),
  fTOFTime(0x0)
{
  // copy constructor
  if (fTPCdEdx) delete fTPCdEdx;
  fTPCdEdx = new TH2D(*aCut.fTPCdEdx);

  if (fTOFTime) delete fTOFTime;
  fTOFTime = new TH2D(*aCut.fTOFTime);

  if (ftofHist) delete ftofHist; 
  ftofHist= new TH2D(*aCut.ftofHist);
}

AliFemtoCutMonitorParticlePID::~AliFemtoCutMonitorParticlePID()
{
  // Destructor
  delete fTPCdEdx;
  delete fTOFTime;
  delete ftofHist;
}

AliFemtoCutMonitorParticlePID& AliFemtoCutMonitorParticlePID::operator=(const AliFemtoCutMonitorParticlePID& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fTPCdEdx) delete fTPCdEdx;
  fTPCdEdx = new TH2D(*aCut.fTPCdEdx);

  if (fTOFTime) delete fTOFTime;
  fTOFTime = new TH2D(*aCut.fTOFTime);
  
  if(ftofHist) delete ftofHist;
  ftofHist = new TH2D(*aCut.ftofHist);
  return *this;
}

AliFemtoString AliFemtoCutMonitorParticlePID::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorParticlePID report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorParticlePID::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track
  float tMom = aTrack->P().Mag();
  float tdEdx = aTrack->TPCsignal();
  float tTOF = 0.0;
  short tchg = aTrack->Charge();
  if (fTOFParticle == 0) tTOF = aTrack->TOFpionTime();
  if (fTOFParticle == 1) tTOF = aTrack->TOFkaonTime();
  if (fTOFParticle == 2) tTOF = aTrack->TOFprotonTime();

  fTPCdEdx->Fill(tMom, tdEdx);
  fTOFTime->Fill(tMom, tTOF);

  float vp= aTrack->VTOF();
  ftofHist->Fill(vp,tMom);
}

void AliFemtoCutMonitorParticlePID::Write()
{
  // Write out the relevant histograms
  fTPCdEdx->Write();
  fTOFTime->Write();
  ftofHist->Write();
}

TList *AliFemtoCutMonitorParticlePID::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fTPCdEdx);
  tOutputList->Add(fTOFTime);
  tOutputList->Add(ftofHist);
  return tOutputList;
}
