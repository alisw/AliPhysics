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

AliFemtoCutMonitorParticleYPt::AliFemtoCutMonitorParticleYPt():
  fYPt(0),
  fMass(0.13957)
{
  // Default constructor
  fYPt = new TH2D("YPt", "Rapidity vs Pt", 100, -1.0, 1.0, 100, 0.1, 2.0);
}

AliFemtoCutMonitorParticleYPt::AliFemtoCutMonitorParticleYPt(const char *aName, float aMass):
  fYPt(0),
  fMass(aMass)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "YPt%s", aName);
  fYPt = new TH2D(name, "Rapdity vs Pt", 100, -1.0, 1.0, 100, 0.1, 2.0);
}

AliFemtoCutMonitorParticleYPt::AliFemtoCutMonitorParticleYPt(const AliFemtoCutMonitorParticleYPt &aCut):
  fYPt(0),
  fMass(0.13957)
{
  // copy constructor
  if (fYPt) delete fYPt;
  fYPt = new TH2D(*aCut.fYPt);
  fMass = aCut.fMass; 
}

AliFemtoCutMonitorParticleYPt::~AliFemtoCutMonitorParticleYPt()
{
  // Destructor
  delete fYPt;
}

AliFemtoCutMonitorParticleYPt& AliFemtoCutMonitorParticleYPt::operator=(const AliFemtoCutMonitorParticleYPt& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fYPt) delete fYPt;
  fYPt = new TH2D(*aCut.fYPt);
  
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
  float tEnergy = ::sqrt(aTrack->P().mag2()+fMass*fMass);
  float tRapidity = 0.5*::log((tEnergy+aTrack->P().z())/(tEnergy-aTrack->P().z()));
  float tPt = ::sqrt((aTrack->P().x())*(aTrack->P().x())+(aTrack->P().y())*(aTrack->P().y()));
  fYPt->Fill(tRapidity, tPt);
}

void AliFemtoCutMonitorParticleYPt::Write()
{
  // Write out the relevant histograms
  fYPt->Write();
}

TList *AliFemtoCutMonitorParticleYPt::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fYPt);

  return tOutputList;
}
