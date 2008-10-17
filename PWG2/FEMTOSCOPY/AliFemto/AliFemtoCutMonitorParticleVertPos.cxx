////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticleVertPos - the cut monitor for particles to study    //
// the difference between reconstructed and true momentum                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelGlobalHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorParticleVertPos::AliFemtoCutMonitorParticleVertPos():
  fVertPos(0),
  fEtaZ(0),
  fRadPos(0)
{
  // Default constructor
  fVertPos = new TH2D("VertPos", "Vertex position", 200, -20.0, 20.0, 200, -20.0, 20.0);
  fEtaZ    = new TH2D("EtaZPos", "Z vs. Eta", 200, -100.0, 100.0, 100, -1.5, 1.5);
  fRadPos  = new TH1D("RadPos",  "Radial position", 200, 0.0, 1.0);
}

AliFemtoCutMonitorParticleVertPos::AliFemtoCutMonitorParticleVertPos(const char *aName):
  AliFemtoCutMonitor(),
  fVertPos(0),
  fEtaZ(0),
  fRadPos(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "VertPos%s", aName);
  fVertPos = new TH2D(name, "Rapdity vs Pt", 200, -20.0, 20.0, 200, -20.0, 20.0);
  snprintf(name, 200, "EtaZPos%s", aName);
  fEtaZ    = new TH2D(name, "Z vs. Eta", 200, -100.0, 100.0, 100, -1.5, 1.5);
  snprintf(name, 200, "RadPos%s", aName);
  fRadPos  = new TH1D(name,  "Radial position", 200, 0.0, 1.0);
}

AliFemtoCutMonitorParticleVertPos::AliFemtoCutMonitorParticleVertPos(const AliFemtoCutMonitorParticleVertPos &aCut):
  AliFemtoCutMonitor(),
  fVertPos(0),
  fEtaZ(0),
  fRadPos(0)
{
  // copy constructor
  if (fVertPos) delete fVertPos;
  fVertPos = new TH2D(*aCut.fVertPos);
  if (fEtaZ) delete fEtaZ;
  fEtaZ = new TH2D(*aCut.fEtaZ);
  if (fRadPos) delete fRadPos;
  fRadPos = new TH1D(*aCut.fRadPos);
}

AliFemtoCutMonitorParticleVertPos::~AliFemtoCutMonitorParticleVertPos()
{
  // Destructor
  delete fVertPos;
  delete fEtaZ;
  delete fRadPos;
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
  if (fRadPos) delete fRadPos;
  fRadPos = new TH1D(*aCut.fRadPos);
  
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
  AliFemtoModelGlobalHiddenInfo *hinfo = dynamic_cast<AliFemtoModelGlobalHiddenInfo *>(aTrack->GetHiddenInfo());
  if (hinfo) {
    float tEta = -TMath::Log(TMath::Tan(hinfo->GetTrueMomentum()->theta()/2.0));

    fVertPos->Fill(hinfo->GetGlobalEmissionPoint()->x(), hinfo->GetGlobalEmissionPoint()->y());
    fEtaZ->Fill(hinfo->GetGlobalEmissionPoint()->z(), tEta);
    fRadPos->Fill(hinfo->GetGlobalEmissionPoint()->perp());
  }
}

void AliFemtoCutMonitorParticleVertPos::Write()
{
  // Write out the relevant histograms
  fVertPos->Write();
  fEtaZ->Write();
  fRadPos->Write();
}

TList *AliFemtoCutMonitorParticleVertPos::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fVertPos);
  tOutputList->Add(fEtaZ);
  tOutputList->Add(fRadPos);

  return tOutputList;
}
