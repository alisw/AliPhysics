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
  fEvMult(0),
  fNormEvMult(0),
  fSPDMult(0)
{
  // Default constructor
  fEvMult = new TH1D("EvMult", "Event Multiplicity", 5001, -0.5, 5000.5);
}

AliFemtoCutMonitorEventMult::AliFemtoCutMonitorEventMult(const char *aName):
  AliFemtoCutMonitor(),
  fEvMult(0),
  fNormEvMult(0),
  fSPDMult(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "EvMult%s", aName);
  fEvMult = new TH1D(name, "Event Multiplicity", 5001, -0.5, 5000.5);

  snprintf(name, 200, "NormEvMult%s", aName);
  fNormEvMult = new TH1D(name, "Normalized Event Multiplicity", 5001, -0.5, 5000.5);

  snprintf(name, 200, "SPDEvMult%s", aName);
  fSPDMult = new TH1D(name, "SPD Tracklet Multiplicity", 5001, -0.5, 5000.5);
}

AliFemtoCutMonitorEventMult::AliFemtoCutMonitorEventMult(const AliFemtoCutMonitorEventMult &aCut):
  AliFemtoCutMonitor(),
  fEvMult(0),
  fNormEvMult(0),
  fSPDMult(0)
{
  // copy constructor
  if (fEvMult) delete fEvMult;
  fEvMult = new TH1D(*aCut.fEvMult);

  if (fNormEvMult) delete fNormEvMult;
  fNormEvMult = new TH1D(*aCut.fNormEvMult);

  if (fSPDMult) delete fSPDMult;
  fSPDMult = new TH1D(*aCut.fSPDMult);
}

AliFemtoCutMonitorEventMult::~AliFemtoCutMonitorEventMult()
{
  // Destructor
  delete fEvMult;
  delete fNormEvMult;
  delete fSPDMult;
}

AliFemtoCutMonitorEventMult& AliFemtoCutMonitorEventMult::operator=(const AliFemtoCutMonitorEventMult& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fEvMult) delete fEvMult;
  fEvMult = new TH1D(*aCut.fEvMult);
  
  if (fNormEvMult) delete fNormEvMult;
  fNormEvMult = new TH1D(*aCut.fNormEvMult);
  
  if (fSPDMult) delete fSPDMult;
  fSPDMult = new TH1D(*aCut.fSPDMult);
  
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
  fNormEvMult->Fill(aEvent->UncorrectedNumberOfPrimaries());
  fSPDMult->Fill(aEvent->SPDMultiplicity());
}

void AliFemtoCutMonitorEventMult::Write()
{
  // Write out the relevant histograms
  fEvMult->Write();
  fNormEvMult->Write();
  fSPDMult->Write();
}

TList *AliFemtoCutMonitorEventMult::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fEvMult);
  tOutputList->Add(fNormEvMult);
  tOutputList->Add(fSPDMult);
  return tOutputList;
}
