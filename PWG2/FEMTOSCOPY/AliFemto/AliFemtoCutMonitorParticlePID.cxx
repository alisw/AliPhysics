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
  fTPCdEdx(0)
{
  // Default constructor
  fTPCdEdx = new TH2D("TPCdEdx", "TPC dEdx vs momentum", 190, 0.1, 2.0, 500, 0.0, 500.0);
}

AliFemtoCutMonitorParticlePID::AliFemtoCutMonitorParticlePID(const char *aName):
  AliFemtoCutMonitor(),
  fTPCdEdx(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "TPCdEdx%s", aName);
  fTPCdEdx = new TH2D(name, "TPC dEdx vs. momentum", 190, 0.1, 2.0, 500, 0.0, 500.0);
}

AliFemtoCutMonitorParticlePID::AliFemtoCutMonitorParticlePID(const AliFemtoCutMonitorParticlePID &aCut):
  AliFemtoCutMonitor(),
  fTPCdEdx(0)
{
  // copy constructor
  if (fTPCdEdx) delete fTPCdEdx;
  fTPCdEdx = new TH2D(*aCut.fTPCdEdx);
}

AliFemtoCutMonitorParticlePID::~AliFemtoCutMonitorParticlePID()
{
  // Destructor
  delete fTPCdEdx;
}

AliFemtoCutMonitorParticlePID& AliFemtoCutMonitorParticlePID::operator=(const AliFemtoCutMonitorParticlePID& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fTPCdEdx) delete fTPCdEdx;
  fTPCdEdx = new TH2D(*aCut.fTPCdEdx);
  
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

  fTPCdEdx->Fill(tMom, tdEdx);
}

void AliFemtoCutMonitorParticlePID::Write()
{
  // Write out the relevant histograms
  fTPCdEdx->Write();
}

TList *AliFemtoCutMonitorParticlePID::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fTPCdEdx);

  return tOutputList;
}
