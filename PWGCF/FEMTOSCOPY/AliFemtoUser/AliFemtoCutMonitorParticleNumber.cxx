////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticleNumber - the cut monitor for particles to study    //
// the difference between reconstructed and true momentum                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticleNumber.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH3F.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorParticleNumber::AliFemtoCutMonitorParticleNumber():
  fNtrig(0)
{
  // Default constructor
}

AliFemtoCutMonitorParticleNumber::AliFemtoCutMonitorParticleNumber(const char *title, const int &pT1Bins=1, const double& pT1min=0, const double& pT1max=4, const int &zvtxBins=10, const double& zvtxmin=-10, const double& zvtxmax=10, const int &multBins=5, const int& multmin=0, const int& multmax=100):
  AliFemtoCutMonitor(),
  fNtrig(0)
{
  // Normal constructor
 
  // set up numerator
  char tTitNtrig[101] = "Ntrig";
  strncat(tTitNtrig,title, 100);
  fNtrig = new TH3F(tTitNtrig,title,pT1Bins,pT1min,pT1max,multBins,multmin,multmax,zvtxBins,zvtxmin,zvtxmax);
  
  // to enable error bar calculation...
  fNtrig->Sumw2();
}

AliFemtoCutMonitorParticleNumber::AliFemtoCutMonitorParticleNumber(const AliFemtoCutMonitorParticleNumber &aCut):
  AliFemtoCutMonitor(),
  fNtrig(0)
{
  // copy constructor
  fNtrig= new TH3F(*aCut.fNtrig);
  
}

AliFemtoCutMonitorParticleNumber::~AliFemtoCutMonitorParticleNumber()
{
  // Destructor
  delete fNtrig;
}

AliFemtoCutMonitorParticleNumber& AliFemtoCutMonitorParticleNumber::operator=(const AliFemtoCutMonitorParticleNumber& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fNtrig) delete fNtrig;
  fNtrig = new TH3F(*aCut.fNtrig);
    
  return *this;
}

AliFemtoString AliFemtoCutMonitorParticleNumber::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorParticleNumber report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorParticleNumber::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track

  double mult = -1;
  double zvtx = -11;

  double px1 = aTrack->P().x();
  double py1 = aTrack->P().y();
  double pt1 = TMath::Hypot(px1, py1);
  
  mult = aTrack->Multiplicity();
  zvtx = aTrack->Zvtx();
  fNtrig->Fill(pt1,mult,zvtx);  
 
}


void AliFemtoCutMonitorParticleNumber::Fill(const AliFemtoV0* aV0)
{
  // Fill in the monitor histograms with the values from the current V0

  double mult = -1;
  double zvtx = -11;
  
  double pt1 = aV0->PtV0();
    
  mult = aV0->Multiplicity();
  zvtx = aV0->Zvtx();
  fNtrig->Fill(pt1,mult,zvtx);  
 
}

void AliFemtoCutMonitorParticleNumber::Write()
{
  // Write out the relevant histograms
  fNtrig->Write();
}

TList *AliFemtoCutMonitorParticleNumber::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fNtrig);

  return tOutputList;
}
