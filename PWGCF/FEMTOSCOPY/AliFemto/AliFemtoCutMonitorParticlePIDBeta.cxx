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
  , fDifference(nullptr)
{
  // Default constructor
  fBeta   = new TH2D("Beta", "Beta vs. momentum", 100, 0.0,10.0, 250, 0.0, 1.0);
  fMass     = new TH2D("Mass", "m vs. p", 100, 0.0, 5.0, 250, -1.0, 10.0);
  fDifference     = new TH1D("Difference", "counts vs. mTOF2-mPDG2", 500, -5.0, 5.0);
}


AliFemtoCutMonitorParticlePIDBeta::AliFemtoCutMonitorParticlePIDBeta(const char *aName, Int_t aTOFParticle, Double_t yTOFTimeMin, Double_t yTOFTimeMax):
  AliFemtoCutMonitorParticlePID(aName,aTOFParticle,yTOFTimeMin,yTOFTimeMax)
  , fBeta(nullptr)
  , fMass(nullptr)
  , fDifference(nullptr)
{
  // Normal constructor
    fBeta     = new TH2D(TString::Format("Beta%s", aName), "Beta vs. momentum", 100, 0.0, 10.0, 250, 0.0, 1.0);
    fMass     = new TH2D(TString::Format("Mass%s", aName), "m2 vs. p", 100, 0.0, 5.0, 250, -1.0, 10.0);
    fDifference     = new TH1D(TString::Format("Difference%s", aName), "counts vs. mTOF2-mPDG2", 500, -5.0, 5.0);

}

AliFemtoCutMonitorParticlePIDBeta::AliFemtoCutMonitorParticlePIDBeta(const AliFemtoCutMonitorParticlePIDBeta &aCut):
  AliFemtoCutMonitorParticlePID(aCut)
  , fBeta(new TH2D(*aCut.fBeta))
  , fMass(new TH2D(*aCut.fMass))
  , fDifference(new TH1D(*aCut.fDifference)) 
{
  // copy constructor
}

AliFemtoCutMonitorParticlePIDBeta::~AliFemtoCutMonitorParticlePIDBeta()
{
  // Destructor
  delete fBeta;
  delete fMass;
  delete fDifference;
}

AliFemtoCutMonitorParticlePIDBeta& AliFemtoCutMonitorParticlePIDBeta::operator=(const AliFemtoCutMonitorParticlePIDBeta& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  AliFemtoCutMonitorParticlePID::operator=(aCut);
 
  *fBeta = *aCut.fBeta;
  *fMass = *aCut.fMass;
  *fDifference = *aCut.fDifference;
  return *this;
}


void AliFemtoCutMonitorParticlePIDBeta::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track
  float tMom = ( AliFemtoCutMonitorParticlePID::fIfUsePt) ? aTrack->Pt() : aTrack->P().Mag();
  float c=1;
  AliFemtoCutMonitorParticlePID::Fill(aTrack);
  double beta =aTrack->VTOF();
  fBeta->Fill(tMom, beta);
  double massTOF;//Mass square
  
  double massPDGPi=0.13957018;
  double massPDGK=0.493677;
  double massPDGP=0.938272013;
  double massPDGD=1.8756;
  
  if(beta!=0){
    massTOF= tMom*tMom/c/c*(1/(beta*beta)-1);
 
    fMass->Fill(tMom,massTOF);
    if(AliFemtoCutMonitorParticlePID::fTOFParticle==0)
      fDifference->Fill(massTOF-massPDGPi*massPDGPi);
    else if(AliFemtoCutMonitorParticlePID::fTOFParticle==1)
      fDifference->Fill(massTOF-massPDGK*massPDGK);
    else if(AliFemtoCutMonitorParticlePID::fTOFParticle==2)
      fDifference->Fill(massTOF-massPDGP*massPDGP);
    else if(AliFemtoCutMonitorParticlePID::fTOFParticle==3)
      fDifference->Fill(massTOF-massPDGD*massPDGD);
    
    
  }

}

void AliFemtoCutMonitorParticlePIDBeta::Write()
{
  // Write out the relevant histograms
  AliFemtoCutMonitorParticlePID::Write();
  fBeta->Write();
  fMass->Write();
  fDifference->Write();
}

TList *AliFemtoCutMonitorParticlePIDBeta::GetOutputList()
{
  TList *tOutputList = AliFemtoCutMonitorParticlePID::GetOutputList();
  tOutputList->Add(fBeta);
  tOutputList->Add(fMass);
  tOutputList->Add(fDifference);

  return tOutputList;
}

