////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorEventVertex - the cut monitor for events to study        //
// the distribution and error of the primary vertex                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoEvent.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorEventVertex::AliFemtoCutMonitorEventVertex():
  AliFemtoCutMonitor(),
  fEvVertRad(0),
  fEvVertXY(0),
  fEvVertSigXY(0),
  fEvVertZ(0),
  fEvVertSigZ(0)
{
  // Default constructor
  fEvVertRad   = new TH1D("EvVertRad",   "Vertex position radial", 200, 0.0, 0.2);
  fEvVertXY    = new TH2D("EvVertXY",    "Vertex position xy plane", 200, -0.2, 0.2, 200, -0.2, 0.2);
  fEvVertSigXY = new TH1D("EvVertSigXY", "Vertex error in xy plane", 200, 0.0, 0.2);
  fEvVertZ     = new TH1D("EvVertZ",     "Vertex position in z", 500, -50.0, 50.0);
  fEvVertSigZ  = new TH1D("EvVertSigZ",  "Vertex error in z", 100, 0.0, 0.2);
}

AliFemtoCutMonitorEventVertex::AliFemtoCutMonitorEventVertex(const char *aName):
  AliFemtoCutMonitor(),
  fEvVertRad(0),
  fEvVertXY(0),
  fEvVertSigXY(0),
  fEvVertZ(0),
  fEvVertSigZ(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "EvVertRad%s", aName);
  fEvVertRad   = new TH1D(name,   "Vertex position radial", 200, 0.0, 0.2);
  snprintf(name, 200, "EvVertXY%s", aName);
  fEvVertXY    = new TH2D(name,    "Vertex position xy plane", 200, -0.2, 0.2, 200, -0.2, 0.2);
  snprintf(name, 200, "EvVertSigXY%s", aName);
  fEvVertSigXY = new TH1D(name, "Vertex error in xy plane", 200, 0.0, 0.2);
  snprintf(name, 200, "EvVertZ%s", aName);
  fEvVertZ     = new TH1D(name,     "Vertex position in z", 500, -50.0, 50.0);
  snprintf(name, 200, "EvVertSigZ%s", aName);
  fEvVertSigZ  = new TH1D(name,  "Vertex error in z", 100, 0.0, 0.2);
}

AliFemtoCutMonitorEventVertex::AliFemtoCutMonitorEventVertex(const AliFemtoCutMonitorEventVertex &aCut):
  AliFemtoCutMonitor(),
  fEvVertRad(0),
  fEvVertXY(0),
  fEvVertSigXY(0),
  fEvVertZ(0),
  fEvVertSigZ(0)
{
  // copy constructor
  if (fEvVertRad) delete fEvVertRad;
  fEvVertRad = new TH1D(*aCut.fEvVertRad);
  if (fEvVertXY) delete fEvVertXY;
  fEvVertXY = new TH2D(*aCut.fEvVertXY);
  if (fEvVertSigXY) delete fEvVertSigXY;
  fEvVertSigXY = new TH1D(*aCut.fEvVertSigXY);
  if (fEvVertZ) delete fEvVertZ;
  fEvVertZ = new TH1D(*aCut.fEvVertZ);
  if (fEvVertSigZ) delete fEvVertSigZ;
  fEvVertSigZ = new TH1D(*aCut.fEvVertSigZ);
}

AliFemtoCutMonitorEventVertex::~AliFemtoCutMonitorEventVertex()
{
  // Destructor
  delete fEvVertRad;
  delete fEvVertXY;
  delete fEvVertSigXY;
  delete fEvVertZ;
  delete fEvVertSigZ;
}

AliFemtoCutMonitorEventVertex& AliFemtoCutMonitorEventVertex::operator=(const AliFemtoCutMonitorEventVertex& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fEvVertRad) delete fEvVertRad;
  fEvVertRad = new TH1D(*aCut.fEvVertRad);
  if (fEvVertXY) delete fEvVertXY;
  fEvVertXY = new TH2D(*aCut.fEvVertXY);
  if (fEvVertSigXY) delete fEvVertSigXY;
  fEvVertSigXY = new TH1D(*aCut.fEvVertSigXY);
  if (fEvVertZ) delete fEvVertZ;
  fEvVertZ = new TH1D(*aCut.fEvVertZ);
  if (fEvVertSigZ) delete fEvVertSigZ;
  fEvVertSigZ = new TH1D(*aCut.fEvVertSigZ);
  
  return *this;
}

AliFemtoString AliFemtoCutMonitorEventVertex::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorEventVertex report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorEventVertex::Fill(const AliFemtoEvent* aEvent)
{
  // Fill in the monitor histograms with the values from the current track
  fEvVertRad->Fill(TMath::Hypot(aEvent->PrimVertPos().x(), aEvent->PrimVertPos().y()));
  fEvVertXY->Fill(aEvent->PrimVertPos().x(), aEvent->PrimVertPos().y());
  fEvVertSigXY->Fill(TMath::Sqrt(aEvent->PrimVertCov()[0]+aEvent->PrimVertCov()[2]));
  fEvVertZ->Fill(aEvent->PrimVertPos().z());
  fEvVertSigZ->Fill(TMath::Sqrt(aEvent->PrimVertCov()[5]));
}

void AliFemtoCutMonitorEventVertex::Write()
{
  // Write out the relevant histograms
  fEvVertRad->Write();
  fEvVertXY->Write();
  fEvVertSigXY->Write();
  fEvVertZ->Write();
  fEvVertSigZ->Write();
}

TList *AliFemtoCutMonitorEventVertex::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fEvVertRad);
  tOutputList->Add(fEvVertXY);
  tOutputList->Add(fEvVertSigXY);
  tOutputList->Add(fEvVertZ);
  tOutputList->Add(fEvVertSigZ);

  return tOutputList;
}
