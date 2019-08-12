////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorTrackTPCncls - the cut monitor for tracks to study     ///
/// the number of TPC clusters distribution.                                 ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorTrackTPCncls.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoEvent.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>


AliFemtoCutMonitorTrackTPCncls::AliFemtoCutMonitorTrackTPCncls():
  fTrTPCncls(0)
{
  // Default constructor
  fTrTPCncls = new TH1D("TrTPCncls", "Track TPC Clusters", 5001, -0.5, 5000.5);
}

AliFemtoCutMonitorTrackTPCncls::AliFemtoCutMonitorTrackTPCncls(const char *aName):
  AliFemtoCutMonitor(),
  fTrTPCncls(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "TrTPCncls%s", aName);
  fTrTPCncls = new TH1D(name, "Track TPC Clusters", 5001, -0.5, 5000.5);
}

AliFemtoCutMonitorTrackTPCncls::AliFemtoCutMonitorTrackTPCncls(const AliFemtoCutMonitorTrackTPCncls &aCut):
  AliFemtoCutMonitor(),
  fTrTPCncls(0)
{
  // copy constructor
  if (fTrTPCncls) delete fTrTPCncls;
  fTrTPCncls = new TH1D(*aCut.fTrTPCncls);
}

AliFemtoCutMonitorTrackTPCncls::~AliFemtoCutMonitorTrackTPCncls()
{
  // Destructor
  delete fTrTPCncls;
}

AliFemtoCutMonitorTrackTPCncls& AliFemtoCutMonitorTrackTPCncls::operator=(const AliFemtoCutMonitorTrackTPCncls& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fTrTPCncls) delete fTrTPCncls;
  fTrTPCncls = new TH1D(*aCut.fTrTPCncls);
  
  return *this;
}

AliFemtoString AliFemtoCutMonitorTrackTPCncls::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorTrackTPCncls report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorTrackTPCncls::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track
  fTrTPCncls->Fill(aTrack->TPCncls());
}

void AliFemtoCutMonitorTrackTPCncls::Write()
{
  // Write out the relevant histograms
  fTrTPCncls->Write();
}

TList *AliFemtoCutMonitorTrackTPCncls::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fTrTPCncls);

  return tOutputList;
}


















