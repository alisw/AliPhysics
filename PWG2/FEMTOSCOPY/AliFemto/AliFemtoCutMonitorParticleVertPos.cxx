////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticleVertPos - the cut monitor for particles to study    //
// the difference between reconstructed and true momentum                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorParticleVertPos::AliFemtoCutMonitorParticleVertPos():
  fVertPos(0),
  fEtaZ(0)
{
  // Default constructor
  fVertPos = new TH2D("VertPos", "Vertex position", 200, -20.0, 20.0, 200, -20.0, 20.0);
  fEtaZ    = new TH2D("EtaZPos", "Z vs. Eta", 200, -100.0, 100.0, 100, -1.5, 1.5);
}

AliFemtoCutMonitorParticleVertPos::AliFemtoCutMonitorParticleVertPos(const char *aName):
  AliFemtoCutMonitor(),
  fVertPos(0),
  fEtaZ(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "VertPos%s", aName);
  fVertPos = new TH2D(name, "Rapdity vs Pt", 200, -20.0, 20.0, 200, -20.0, 20.0);
  snprintf(name, 200, "EtaZPos%s", aName);
  fEtaZ    = new TH2D(name, "Z vs. Eta", 200, -100.0, 100.0, 100, -1.5, 1.5);
}

AliFemtoCutMonitorParticleVertPos::AliFemtoCutMonitorParticleVertPos(const AliFemtoCutMonitorParticleVertPos &aCut):
  AliFemtoCutMonitor(),
  fVertPos(0),
  fEtaZ(0)

{
  // copy constructor
  if (fVertPos) delete fVertPos;
  fVertPos = new TH2D(*aCut.fVertPos);
  if (fEtaZ) delete fEtaZ;
  fEtaZ = new TH2D(*aCut.fEtaZ);
}

AliFemtoCutMonitorParticleVertPos::~AliFemtoCutMonitorParticleVertPos()
{
  // Destructor
  delete fVertPos;
  delete fEtaZ;
}

AliFemtoCutMonitorParticleVertPos& AliFemtoCutMonitorParticleVertPos::operator=(const AliFemtoCutMonitorParticleVertPos& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fVertPos) delete fVertPos;
  fVertPos = new TH2D(*aCut.fVertPos);
  if (fEtaZ) delete fEtaZ;
  fEtaZ = new TH2D(*aCut.fEtaZ);
  
  return *this;
}

AliFemtoString AliFemtoCutMonitorParticleVertPos::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorParticleVertPos report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorParticleVertPos::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track
  AliFemtoModelHiddenInfo *hinfo = dynamic_cast<AliFemtoModelHiddenInfo *>(aTrack->GetHiddenInfo());
  if (hinfo) {
    float tEta = -TMath::Log(TMath::Tan(hinfo->GetTrueMomentum()->theta()/2.0));

    fVertPos->Fill(hinfo->GetEmissionPoint()->x(), hinfo->GetEmissionPoint()->y());
    fEtaZ->Fill(hinfo->GetEmissionPoint()->z(), tEta);
  }
}

void AliFemtoCutMonitorParticleVertPos::Write()
{
  // Write out the relevant histograms
  fVertPos->Write();
  fEtaZ->Write();
}

TList *AliFemtoCutMonitorParticleVertPos::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fVertPos);
  tOutputList->Add(fEtaZ);

  return tOutputList;
}
