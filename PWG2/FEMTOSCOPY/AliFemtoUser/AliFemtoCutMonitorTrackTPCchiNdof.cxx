////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorTrackTPCchiNdof - the cut monitor for tracks to study  ///
/// the number of TPC clusters distribution.                                 ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorTrackTPCchiNdof.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoEvent.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>


AliFemtoCutMonitorTrackTPCchiNdof::AliFemtoCutMonitorTrackTPCchiNdof():
  fTrTPCchiNdof(0)
{
  // Default constructor
  fTrTPCchiNdof = new TH1D("TrTPCchiNdof", "Track TPC Clusters", 5001, -0.5, 5000.5);
}

AliFemtoCutMonitorTrackTPCchiNdof::AliFemtoCutMonitorTrackTPCchiNdof(const char *aName):
  AliFemtoCutMonitor(),
  fTrTPCchiNdof(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "TrTPCchiNdof%s", aName);
  fTrTPCchiNdof = new TH1D(name, "Track TPC Clusters", 5001, -0.5, 5000.5);
}

AliFemtoCutMonitorTrackTPCchiNdof::AliFemtoCutMonitorTrackTPCchiNdof(const AliFemtoCutMonitorTrackTPCchiNdof &aCut):
  AliFemtoCutMonitor(),
  fTrTPCchiNdof(0)
{
  // copy constructor
  if (fTrTPCchiNdof) delete fTrTPCchiNdof;
  fTrTPCchiNdof = new TH1D(*aCut.fTrTPCchiNdof);
}

AliFemtoCutMonitorTrackTPCchiNdof::~AliFemtoCutMonitorTrackTPCchiNdof()
{
  // Destructor
  delete fTrTPCchiNdof;
}

AliFemtoCutMonitorTrackTPCchiNdof& AliFemtoCutMonitorTrackTPCchiNdof::operator=(const AliFemtoCutMonitorTrackTPCchiNdof& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fTrTPCchiNdof) delete fTrTPCchiNdof;
  fTrTPCchiNdof = new TH1D(*aCut.fTrTPCchiNdof);
  
  return *this;
}

AliFemtoString AliFemtoCutMonitorTrackTPCchiNdof::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorTrackTPCchiNdof report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorTrackTPCchiNdof::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track
  if (aTrack->TPCncls() > 0) {
    fTrTPCchiNdof->Fill(aTrack->TPCchi2()/aTrack->TPCncls());
  }
}

void AliFemtoCutMonitorTrackTPCchiNdof::Write()
{
  // Write out the relevant histograms
  fTrTPCchiNdof->Write();
}

TList *AliFemtoCutMonitorTrackTPCchiNdof::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fTrTPCchiNdof);

  return tOutputList;
}


















