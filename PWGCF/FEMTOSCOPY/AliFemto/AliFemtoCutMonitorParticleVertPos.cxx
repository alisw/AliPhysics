///
/// \file AliFemtoCutMonitorParticleVertPos.cxx
///

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
  fRadPos(0),
  fEmPointX(0),
  fEmPointY(0),
  fEmPointZ(0),
  fEmPointT(0)
{
  // Default constructor
  fVertPos = new TH2D("VertPos", "Vertex position", 200, -20.0, 20.0, 200, -20.0, 20.0);
  fEtaZ    = new TH2D("EtaZPos", "Z vs. Eta", 200, -100.0, 100.0, 100, -1.5, 1.5);
  fRadPos  = new TH1D("RadPos",  "Radial position", 200, 0.0, 1.0);
  fEmPointX = new TH1D("EmPointX","Emission point x", 400, -200.0, 200.0);
  fEmPointY = new TH1D("EmPointY","Emission point y", 400, -200.0, 200.0);
  fEmPointZ = new TH1D("EmPointZ","Emission point z", 400, -200.0, 200.0);
  fEmPointT = new TH1D("EmPointT","Emission point t", 400, -200.0, 200.0);

}

AliFemtoCutMonitorParticleVertPos::AliFemtoCutMonitorParticleVertPos(const char *aName):
  AliFemtoCutMonitor(),
  fVertPos(0),
  fEtaZ(0),
  fRadPos(0),
  fEmPointX(0),
  fEmPointY(0),
  fEmPointZ(0),
  fEmPointT(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "VertPos%s", aName);
  fVertPos = new TH2D(name, "Rapdity vs Pt", 200, -20.0, 20.0, 200, -20.0, 20.0);
  snprintf(name, 200, "EtaZPos%s", aName);
  fEtaZ    = new TH2D(name, "Z vs. Eta", 200, -100.0, 100.0, 100, -1.5, 1.5);
  snprintf(name, 200, "RadPos%s", aName);
  fRadPos  = new TH1D(name,  "Radial position", 200, 0.0, 1.0);
  snprintf(name, 200, "EmPosX%s", aName);
  fEmPointX = new TH1D(name,"Emission point x", 400, -200.0, 200.0);
  snprintf(name, 200, "EmPosY%s", aName);
  fEmPointY = new TH1D(name,"Emission point y", 400, -200.0, 200.0);
  snprintf(name, 200, "EmPosZ%s", aName);
  fEmPointZ = new TH1D(name,"Emission point z", 400, -200.0, 200.0);
  snprintf(name, 200, "EmPosT%s", aName);
  fEmPointT = new TH1D(name,"Emission point t", 400, -200.0, 200.0);
}

AliFemtoCutMonitorParticleVertPos::AliFemtoCutMonitorParticleVertPos(const AliFemtoCutMonitorParticleVertPos &aCut):
  AliFemtoCutMonitor(),
  fVertPos(0),
  fEtaZ(0),
  fRadPos(0),
  fEmPointX(0),
  fEmPointY(0),
  fEmPointZ(0),
  fEmPointT(0)
{
  // copy constructor
  fVertPos = new TH2D(*aCut.fVertPos);
  fEtaZ = new TH2D(*aCut.fEtaZ);
  fRadPos = new TH1D(*aCut.fRadPos);
  fEmPointX = new TH1D(*aCut.fEmPointX);
  fEmPointY = new TH1D(*aCut.fEmPointY);
  fEmPointZ = new TH1D(*aCut.fEmPointZ);
  fEmPointT = new TH1D(*aCut.fEmPointT);
}

AliFemtoCutMonitorParticleVertPos::~AliFemtoCutMonitorParticleVertPos()
{
  // Destructor
  delete fVertPos;
  delete fEtaZ;
  delete fRadPos;
  delete fEmPointX;
  delete fEmPointY;
  delete fEmPointZ;
  delete fEmPointT;
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
  if (fEmPointX) delete fEmPointX;
  fEmPointX = new TH1D(*aCut.fEmPointX);
  if (fEmPointY) delete fEmPointY;
  fEmPointY = new TH1D(*aCut.fEmPointY);
  if (fEmPointZ) delete fEmPointZ;
  fEmPointZ = new TH1D(*aCut.fEmPointZ);
  if (fEmPointT) delete fEmPointT;
  fEmPointT = new TH1D(*aCut.fEmPointT);

  return *this;
}

AliFemtoString AliFemtoCutMonitorParticleVertPos::Report()
{
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
    float tEta = -TMath::Log(TMath::Tan(hinfo->GetTrueMomentum()->Theta()/2.0));

    fVertPos->Fill(hinfo->GetGlobalEmissionPoint()->x(), hinfo->GetGlobalEmissionPoint()->y());
    fEtaZ->Fill(hinfo->GetGlobalEmissionPoint()->z(), tEta);
    fRadPos->Fill(hinfo->GetGlobalEmissionPoint()->Perp());
  }

  AliFemtoModelHiddenInfo *hminfo = dynamic_cast<AliFemtoModelHiddenInfo *>(aTrack->GetHiddenInfo());
  if (hminfo) {
    fEmPointX->Fill(hminfo->GetEmissionPoint()->x());
    fEmPointY->Fill(hminfo->GetEmissionPoint()->y());
    fEmPointZ->Fill(hminfo->GetEmissionPoint()->z());
    fEmPointT->Fill(hminfo->GetEmissionPoint()->t());
  }
}

void AliFemtoCutMonitorParticleVertPos::Write()
{
  // Write out the relevant histograms
  fVertPos->Write();
  fEtaZ->Write();
  fRadPos->Write();
  fEmPointX->Write();
  fEmPointY->Write();
  fEmPointZ->Write();
  fEmPointT->Write();
}

TList *AliFemtoCutMonitorParticleVertPos::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fVertPos);
  tOutputList->Add(fEtaZ);
  tOutputList->Add(fRadPos);
  tOutputList->Add(fEmPointX);
  tOutputList->Add(fEmPointY);
  tOutputList->Add(fEmPointZ);
  tOutputList->Add(fEmPointT);

  return tOutputList;
}
