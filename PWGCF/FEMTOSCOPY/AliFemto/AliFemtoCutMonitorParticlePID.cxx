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
  AliFemtoCutMonitor()
  , fTOFParticle(0)
  , fIfUsePt(false)
  , fTPCdEdx(nullptr)
  , fTOFTime(nullptr)
  , fTOFNSigma(nullptr)
  , fTPCNSigma(nullptr)
  , fTPCTOFNSigma(nullptr)
  , fTPCvsTOFNSigma(nullptr)
  , fParticleOrigin(nullptr)
  , fParticleId(nullptr)
{
  // Default constructor
  fTPCdEdx        = new TH2D("TPCdEdx", "TPC dEdx vs. momentum", 100, 0.0, 5.0, 250, 0.0, 500.0);
  fTOFTime        = new TH2D("TOFTime", "TOF Time vs. momentum", 100, 0.0, 5.0, 400, -4000.0, 4000.0);//to be consistent with dE/dx
  //fTOFTime        = new TH2D("TOFTime", "TOF Time vs. momentum", 100, 0.1, 5.0, 400, -4000.0, 4000.0);
  fTOFNSigma      = new TH2D("TOFNSigma","TOF NSigma vs. momentum", 100, 0.0, 5.0, 100, -5.0, 5.0);
  fTPCNSigma      = new TH2D("TPCNSigma","TPC NSigma vs. momentum", 100, 0.0, 5.0, 100, -5.0, 5.0);
  fTPCTOFNSigma   = new TH2D("TPCTOFNSigma","TPC & TOF NSigma vs. momentum", 100, 0.0, 5.0, 100, 0.0, 10.0);
  fTPCvsTOFNSigma = new TH2D("TPCvsTOFNSigma","TPC vs TOF Nsigma",100, -5.0, 5.0, 100, -5.0, 5.0);
  fParticleOrigin = new TH1D("POrigin", "Mothers PDG Codes", 6000, 0.0, 6000.0);
  fParticleId     = new TH1D("PId", "Particle PDG Codes", 6000, 0.0, 6000.0);

}

AliFemtoCutMonitorParticlePID::AliFemtoCutMonitorParticlePID(const char *aName, Int_t aTOFParticle, Double_t yTOFTimeMin, Double_t yTOFTimeMax):
  AliFemtoCutMonitor()
  , fTOFParticle(aTOFParticle)
  , fIfUsePt(false)
  , fTPCdEdx(nullptr)
  , fTOFTime(nullptr)
  , fTOFNSigma(nullptr)
  , fTPCNSigma(nullptr)
  , fTPCTOFNSigma(nullptr)
  , fTPCvsTOFNSigma(nullptr)
  , fParticleOrigin(nullptr)
  , fParticleId(nullptr)
{
  // Normal constructor
  fTPCdEdx        = new TH2D(TString::Format("TPCdEdx%s", aName), "TPC dEdx vs. momentum", 200, 0.0, 4.0, 250, 0.0, 500.0);
  fTOFTime        = new TH2D(TString::Format("TOFTime%s", aName), "TOF Time vs. momentum", 200, 0.0, 4.0, 400, yTOFTimeMin, yTOFTimeMax);//to be consistent with dE/dx
  //fTOFTime        = new TH2D(TString::Format("TOFTime%s", aName), "TOF Time vs. momentum", 100, 0.1, 5.0, 400, yTOFTimeMin, yTOFTimeMax);
  fTOFNSigma      = new TH2D(TString::Format("TOFNSigma%s", aName), "TOF NSigma vs. momentum", 100, 0.0, 5.0, 100, -5.0, 5.0);
  fTPCNSigma      = new TH2D(TString::Format("TPCNSigma%s", aName), "TPC NSigma vs. momentum", 100, 0.0, 5.0, 100, -5.0, 5.0);
  fTPCTOFNSigma   = new TH2D(TString::Format("TPCTOFNSigma%s", aName), "TPC & TOF NSigma vs. momentum", 100, 0.0, 5.0, 100, 0.0, 10.0);
  fTPCvsTOFNSigma = new TH2D(TString::Format("TPCvsTOFNSigma%s", aName), "TPC vs TOF Nsigma",100, -5.0, 5.0, 100, -5.0, 5.0);
  fParticleOrigin = new TH1D(TString::Format("POrigin%s", aName), "Mothers PDG Codes", 6000, 0.0, 6000.0);
  fParticleId     = new TH1D(TString::Format("PId%s", aName), "Particle PDG Codes", 6000, 0.0, 6000.0);
}

AliFemtoCutMonitorParticlePID::AliFemtoCutMonitorParticlePID(const AliFemtoCutMonitorParticlePID &aCut):
  AliFemtoCutMonitor(aCut)
  , fTOFParticle(aCut.fTOFParticle)
  , fIfUsePt(aCut.fIfUsePt)
  , fTPCdEdx(new TH2D(*aCut.fTPCdEdx))
  , fTOFTime(new TH2D(*aCut.fTOFTime))
  , fTOFNSigma(new TH2D(*aCut.fTOFNSigma))
  , fTPCNSigma(new TH2D(*aCut.fTPCNSigma))
  , fTPCTOFNSigma(new TH2D(*aCut.fTPCTOFNSigma))
  , fTPCvsTOFNSigma(new TH2D(*aCut.fTPCvsTOFNSigma))
  , fParticleOrigin(new TH1D(*aCut.fParticleOrigin))
  , fParticleId(new TH1D(*aCut.fParticleId))
{
  // copy constructor
}

AliFemtoCutMonitorParticlePID::~AliFemtoCutMonitorParticlePID()
{
  // Destructor
  delete fTPCdEdx;
  delete fTOFTime;
  delete fTOFNSigma;
  delete fTPCNSigma;
  delete fTPCTOFNSigma;
  delete fTPCvsTOFNSigma;
  delete fParticleOrigin;
  delete fParticleId;
}

AliFemtoCutMonitorParticlePID& AliFemtoCutMonitorParticlePID::operator=(const AliFemtoCutMonitorParticlePID& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  AliFemtoCutMonitor::operator=(aCut);

  fTOFParticle = aCut.fTOFParticle;
  fIfUsePt = aCut.fIfUsePt;

  *fTPCdEdx = *aCut.fTPCdEdx;
  *fTOFTime = *aCut.fTOFTime;
  *fTOFNSigma = *aCut.fTOFNSigma;
  *fTPCNSigma = *aCut.fTPCNSigma;
  *fTPCTOFNSigma = *aCut.fTPCTOFNSigma;
  *fTPCvsTOFNSigma = *aCut.fTPCvsTOFNSigma;
  *fParticleOrigin = *aCut.fParticleOrigin;
  *fParticleId = *aCut.fParticleId;

  return *this;
}

AliFemtoString AliFemtoCutMonitorParticlePID::Report()
{
  // Prepare report from the execution
  TString report = "*** AliFemtoCutMonitorParticlePID report";
  return AliFemtoString(report.Data());
}

void AliFemtoCutMonitorParticlePID::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track
  float tMom = (fIfUsePt) ? aTrack->Pt() : aTrack->P().Mag();
  float tdEdx = aTrack->TPCsignal();
  float tTOF = 0.0;
  //  short tchg = aTrack->Charge();
  if (fTOFParticle == 0) tTOF = aTrack->TOFpionTime();
  if (fTOFParticle == 1) tTOF = aTrack->TOFkaonTime();
  if (fTOFParticle == 2) tTOF = aTrack->TOFprotonTime();
  if (fTOFParticle == 3) {tTOF = aTrack->TOFdeuteronTime();
  }

  fTPCdEdx->Fill(tMom, tdEdx);
  fTOFTime->Fill(tMom, tTOF);

  AliFemtoModelHiddenInfo *tInfo = (AliFemtoModelHiddenInfo*)aTrack->GetHiddenInfo();
  if(tInfo!=NULL) {
    Int_t partID = TMath::Abs(tInfo->GetPDGPid());
    Int_t motherID = TMath::Abs(tInfo->GetMotherPdgCode());

    fParticleId->Fill(partID);
    fParticleOrigin->Fill(motherID);
  }
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
    if (fTOFParticle == 3) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFD());

    if (fTOFParticle == 0) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCPi());
    if (fTOFParticle == 1) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCK());
    if (fTOFParticle == 2) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCP());
    if (fTOFParticle == 3) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCD());

    if (fTOFParticle == 0) fTPCTOFNSigma->Fill(tMom, TMath::Hypot( aTrack->NSigmaTPCPi(), aTrack->NSigmaTOFPi() ) );
    if (fTOFParticle == 1) fTPCTOFNSigma->Fill(tMom, TMath::Hypot( aTrack->NSigmaTPCK(), aTrack->NSigmaTOFK() ) );
    if (fTOFParticle == 2) fTPCTOFNSigma->Fill(tMom, TMath::Hypot( aTrack->NSigmaTPCP(), aTrack->NSigmaTOFP() ) );
    if (fTOFParticle == 3) fTPCTOFNSigma->Fill(tMom, TMath::Hypot( aTrack->NSigmaTPCD(), aTrack->NSigmaTOFD() ) );
    
    if (fTOFParticle == 0) fTPCvsTOFNSigma->Fill(aTrack->NSigmaTPCPi(), aTrack->NSigmaTOFPi());
    if (fTOFParticle == 1) fTPCvsTOFNSigma->Fill(aTrack->NSigmaTPCK(), aTrack->NSigmaTOFK());
    if (fTOFParticle == 2) fTPCvsTOFNSigma->Fill(aTrack->NSigmaTPCP(), aTrack->NSigmaTOFP());
    if (fTOFParticle == 3) fTPCvsTOFNSigma->Fill(aTrack->NSigmaTPCD(), aTrack->NSigmaTOFD());


}

void AliFemtoCutMonitorParticlePID::Write()
{
  // Write out the relevant histograms
  fTPCdEdx->Write();
  fTOFTime->Write();
  fTOFNSigma->Write();
  fTPCNSigma->Write();
  fTPCTOFNSigma->Write();
  fTPCvsTOFNSigma->Write();
  fParticleId->Write();
  fParticleOrigin->Write();
}

TList *AliFemtoCutMonitorParticlePID::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fTPCdEdx);
  tOutputList->Add(fTOFTime);
  tOutputList->Add(fTOFNSigma);
  tOutputList->Add(fTPCNSigma);
  tOutputList->Add(fTPCTOFNSigma);
  tOutputList->Add(fTPCvsTOFNSigma);
  tOutputList->Add(fParticleId);
  tOutputList->Add(fParticleOrigin);

  return tOutputList;
}
