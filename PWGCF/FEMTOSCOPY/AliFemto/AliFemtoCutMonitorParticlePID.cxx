////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticlePID - the cut monitor for particles to study     //
// various aspects of the PID determination                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorParticlePID::AliFemtoCutMonitorParticlePID():
  fTPCdEdx(0),
  fTOFParticle(0),
  fTOFTime(0x0),
    fTOFNSigma(0),
    fTPCNSigma(0),
    fTPCTOFNSigma(0)
{
  // Default constructor
    fTPCdEdx =  new TH2D("TPCdEdx", "TPC dEdx vs. momentum", 100, 0.0, 5.0, 250, 0.0, 500.0);
    fTOFTime = new TH2D("TOFTime", "TOF Time vs. momentum", 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    fTOFNSigma = new TH2D("TOFNSigma","TOF NSigma vs. momentum", 100, 0.0, 5.0, 100, -5.0, 5.0);
    fTPCNSigma = new TH2D("TPCNSigma","TPC NSigma vs. momentum", 100, 0.0, 5.0, 100, -5.0, 5.0);
    fTPCTOFNSigma = new TH2D("TPCTOFNSigma","TPC & TOF NSigma vs. momentum", 100, 0.0, 5.0, 100, 0.0, 10.0);

}

AliFemtoCutMonitorParticlePID::AliFemtoCutMonitorParticlePID(const char *aName, Int_t aTOFParticle):
  AliFemtoCutMonitor(),
  fTPCdEdx(0),
  fTOFParticle(aTOFParticle),
  fTOFTime(0x0),
    fTOFNSigma(0),
    fTPCNSigma(0),
    fTPCTOFNSigma(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "TPCdEdx%s", aName);
    fTPCdEdx = new TH2D(name, "TPC dEdx vs. momentum", 100, 0.0, 6.0, 250, 0.0, 500.0);

  snprintf(name, 200, "TOFTime%s", aName);
    fTOFTime = new TH2D(name, "TOF Time vs. momentum", 100, 0.1, 5.0, 400, -4000.0, 4000.0);

    snprintf(name, 200, "TOFNSigma%s", aName);
    fTOFNSigma = new TH2D(name,"TOF NSigma vs. momentum", 100, 0.0, 5.0, 100, -5.0, 5.0);

    snprintf(name, 200, "TPCNSigma%s", aName);
    fTPCNSigma = new TH2D(name,"TPC NSigma vs. momentum", 100, 0.0, 5.0, 100, -5.0, 5.0);

    snprintf(name, 200, "TPCTOFNSigma%s", aName);
    fTPCTOFNSigma = new TH2D(name,"TPC & TOF NSigma vs. momentum", 100, 0.0, 5.0, 100, 0.0, 10.0);

}

AliFemtoCutMonitorParticlePID::AliFemtoCutMonitorParticlePID(const AliFemtoCutMonitorParticlePID &aCut):
  AliFemtoCutMonitor(),
  fTPCdEdx(0),
  fTOFParticle(0),
  fTOFTime(0x0),
    fTOFNSigma(0),
    fTPCNSigma(0),
    fTPCTOFNSigma(0)

{
  // copy constructor
  if (fTPCdEdx) delete fTPCdEdx;
  fTPCdEdx = new TH2D(*aCut.fTPCdEdx);

  if (fTOFTime) delete fTOFTime;
  fTOFTime = new TH2D(*aCut.fTOFTime);

    if (fTOFNSigma) delete fTOFNSigma;
    fTOFNSigma= new TH2D(*aCut.fTOFNSigma);

    if (fTPCNSigma) delete fTPCNSigma;
    fTPCNSigma= new TH2D(*aCut.fTPCNSigma);

    if (fTPCTOFNSigma) delete fTPCTOFNSigma;
    fTPCTOFNSigma= new TH2D(*aCut.fTPCTOFNSigma);
}

AliFemtoCutMonitorParticlePID::~AliFemtoCutMonitorParticlePID()
{
  // Destructor
  delete fTPCdEdx;
  delete fTOFTime;
    delete fTOFNSigma;
    delete fTPCNSigma;
    delete fTPCTOFNSigma;

}

AliFemtoCutMonitorParticlePID& AliFemtoCutMonitorParticlePID::operator=(const AliFemtoCutMonitorParticlePID& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fTPCdEdx) delete fTPCdEdx;
  fTPCdEdx = new TH2D(*aCut.fTPCdEdx);

  if (fTOFTime) delete fTOFTime;
  fTOFTime = new TH2D(*aCut.fTOFTime);
  
    if(fTOFNSigma) delete fTOFNSigma;
    fTOFNSigma = new TH2D(*aCut.fTOFNSigma);

    if(fTPCNSigma) delete fTPCNSigma;
    fTPCNSigma = new TH2D(*aCut.fTPCNSigma);

    if(fTPCTOFNSigma) delete fTPCTOFNSigma;
    fTPCTOFNSigma = new TH2D(*aCut.fTPCTOFNSigma);

  return *this;
}

AliFemtoString AliFemtoCutMonitorParticlePID::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorParticlePID report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorParticlePID::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track
  float tMom = aTrack->P().Mag();
  float tdEdx = aTrack->TPCsignal();
  float tTOF = 0.0;
  //  short tchg = aTrack->Charge();
  if (fTOFParticle == 0) tTOF = aTrack->TOFpionTime();
  if (fTOFParticle == 1) tTOF = aTrack->TOFkaonTime();
  if (fTOFParticle == 2) tTOF = aTrack->TOFprotonTime();

  fTPCdEdx->Fill(tMom, tdEdx);
  fTOFTime->Fill(tMom, tTOF);

  //  float vp= aTrack->VTOF();
  //     if (vp > 0.) {
  //         fTOFTime->Fill(tMom, tTOF);
  //         if (fTOFParticle == 0) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFPi());
  //         if (fTOFParticle == 1) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFK());
  //         if (fTOFParticle == 2) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFP());
  // }

    if (fTOFParticle == 0) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFPi());
    if (fTOFParticle == 1) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFK());
    if (fTOFParticle == 2) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFP());

    if (fTOFParticle == 0) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCPi());
    if (fTOFParticle == 1) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCK());
    if (fTOFParticle == 2) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCP());

    if (fTOFParticle == 0) fTPCTOFNSigma->Fill(tMom, TMath::Hypot( aTrack->NSigmaTPCPi(), aTrack->NSigmaTOFPi() )/TMath::Sqrt(2) );
    if (fTOFParticle == 1) fTPCTOFNSigma->Fill(tMom, TMath::Hypot( aTrack->NSigmaTPCK(), aTrack->NSigmaTOFK() )/TMath::Sqrt(2) );
    if (fTOFParticle == 2) fTPCTOFNSigma->Fill(tMom, TMath::Hypot( aTrack->NSigmaTPCP(), aTrack->NSigmaTOFP() )/TMath::Sqrt(2) );

}

void AliFemtoCutMonitorParticlePID::Write()
{
  // Write out the relevant histograms
  fTPCdEdx->Write();
  fTOFTime->Write();
    fTOFNSigma->Write();
    fTPCNSigma->Write();
    fTPCTOFNSigma->Write();

}

TList *AliFemtoCutMonitorParticlePID::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fTPCdEdx);
  tOutputList->Add(fTOFTime);
    tOutputList->Add(fTOFNSigma);
    tOutputList->Add(fTPCNSigma);
    tOutputList->Add(fTPCTOFNSigma);

  return tOutputList;
}
