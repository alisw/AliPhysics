////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticlePIDBeta - the cut monitor for particles to study     //
// various aspects of the PID determination                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticlePIDBeta.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorParticlePIDBeta::AliFemtoCutMonitorParticlePIDBeta():
  AliFemtoCutMonitorParticlePID()
  , fBeta(nullptr)
  , fMass(nullptr)
{
  // Default constructor
  fBeta   = new TH2D("Beta", "Beta vs. momentum", 100, 0.0,10.0, 250, 0.0, 1.0);
  fMass     = new TH2D("Mass", "m vs. p", 100, 0.0, 5.0, 250, -1.0, 10.0);

}


AliFemtoCutMonitorParticlePIDBeta::AliFemtoCutMonitorParticlePIDBeta(const char *aName, Int_t aTOFParticle, Double_t yTOFTimeMin, Double_t yTOFTimeMax):
  AliFemtoCutMonitorParticlePID(aName,aTOFParticle,yTOFTimeMin,yTOFTimeMax)
  , fBeta(nullptr)
  , fMass(nullptr)
{
  // Normal constructor
    fBeta     = new TH2D(TString::Format("Beta%s", aName), "Beta vs. momentum", 100, 0.0, 10.0, 250, 0.0, 1.0);
    fMass     = new TH2D(TString::Format("Mass%s", aName), "m2 vs. p", 100, 0.0, 5.0, 250, -1.0, 10.0);

}

AliFemtoCutMonitorParticlePIDBeta::AliFemtoCutMonitorParticlePIDBeta(const AliFemtoCutMonitorParticlePIDBeta &aCut):
  AliFemtoCutMonitorParticlePID(aCut)
  , fBeta(new TH2D(*aCut.fBeta))
  , fMass(new TH2D(*aCut.fMass))
{
  // copy constructor
}

AliFemtoCutMonitorParticlePIDBeta::~AliFemtoCutMonitorParticlePIDBeta()
{
  // Destructor
  delete fBeta;
  delete fMass;
}

AliFemtoCutMonitorParticlePIDBeta& AliFemtoCutMonitorParticlePIDBeta::operator=(const AliFemtoCutMonitorParticlePIDBeta& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  AliFemtoCutMonitorParticlePID::operator=(aCut);

  *fBeta = *aCut.fBeta;
  *fMass = *aCut.fMass;
  return *this;
}


void AliFemtoCutMonitorParticlePIDBeta::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track
  float tMom = (fIfUsePt) ? aTrack->Pt() : aTrack->P().Mag();
  float c=1;
  AliFemtoCutMonitorParticlePID::Fill(aTrack);
  double beta =aTrack->VTOF();
  if (fTOFParticle == 0) fBeta->Fill(tMom, beta);
  if (fTOFParticle == 1) fBeta->Fill(tMom, beta);
  if (fTOFParticle == 2) fBeta->Fill(tMom, beta);
  if (fTOFParticle == 3) fBeta->Fill(tMom, beta);  

  if (fTOFParticle == 0) fMass->Fill(tMom, tMom*tMom/c/c*(1/(beta*beta)-1));
  if (fTOFParticle == 1) fMass->Fill(tMom, tMom*tMom/c/c*(1/(beta*beta)-1));
  if (fTOFParticle == 2) fMass->Fill(tMom, tMom*tMom/c/c*(1/(beta*beta)-1));
  if (fTOFParticle == 3) fMass->Fill(tMom, tMom*tMom/c/c*(1/(beta*beta)-1));   
}

void AliFemtoCutMonitorParticlePIDBeta::Write()
{
  // Write out the relevant histograms
  AliFemtoCutMonitorParticlePID::Write();
  fBeta->Write();
  fMass->Write();
}

TList *AliFemtoCutMonitorParticlePIDBeta::GetOutputList()
{
  TList *tOutputList = AliFemtoCutMonitorParticlePID::GetOutputList();
  tOutputList->Add(fBeta);
  tOutputList->Add(fMass);
  tOutputList->Add(fParticleOrigin);

  return tOutputList;
}

