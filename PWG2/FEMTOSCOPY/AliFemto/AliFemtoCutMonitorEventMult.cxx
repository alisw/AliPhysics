////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorEventMult - the cut monitor for particles to study    //
// the difference between reconstructed and true momentum                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoEvent.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>

AliFemtoCutMonitorEventMult::AliFemtoCutMonitorEventMult():
  fEvMult(0)
{
  // Default constructor
  fEvMult = new TH1D("EvMult", "Event Multiplicity", 5001, -0.5, 5000.5);
}

AliFemtoCutMonitorEventMult::AliFemtoCutMonitorEventMult(const char *aName):
  AliFemtoCutMonitor(),
  fEvMult(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "EvMult%s", aName);
  fEvMult = new TH1D(name, "Event Multiplicity", 5001, -0.5, 5000.5);
}

AliFemtoCutMonitorEventMult::AliFemtoCutMonitorEventMult(const AliFemtoCutMonitorEventMult &aCut):
  AliFemtoCutMonitor(),
  fEvMult(0)
{
  // copy constructor
  if (fEvMult) delete fEvMult;
  fEvMult = new TH1D(*aCut.fEvMult);
}

AliFemtoCutMonitorEventMult::~AliFemtoCutMonitorEventMult()
{
  // Destructor
  delete fEvMult;
}

AliFemtoCutMonitorEventMult& AliFemtoCutMonitorEventMult::operator=(const AliFemtoCutMonitorEventMult& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fEvMult) delete fEvMult;
  fEvMult = new TH1D(*aCut.fEvMult);
  
  return *this;
}

AliFemtoString AliFemtoCutMonitorEventMult::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorEventMult report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorEventMult::Fill(const AliFemtoEvent* aEvent)
{
  // Fill in the monitor histograms with the values from the current track
  fEvMult->Fill(aEvent->NumberOfTracks());
}

void AliFemtoCutMonitorEventMult::Write()
{
  // Write out the relevant histograms
  fEvMult->Write();
}

TList *AliFemtoCutMonitorEventMult::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fEvMult);

  return tOutputList;
}
