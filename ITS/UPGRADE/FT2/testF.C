#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TSystem.h>
#include <TMath.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TROOT.h>
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "FT2.h"

FT2* det=0;
#endif

TH1F* hdca=0,*hdcaN=0;
TH1F* hdcaZ=0,*hdcaZN=0;
TH1F* hFdca=0,*hFdcaN=0;
TH1F* hFdcaZ=0,*hFdcaZN=0;
TH1F* hPatternITS=0;
TH1F* hFPatternITS=0;

void testF(int ntrials=10000, double dndy=2100., Bool_t useKalmanOut=kTRUE)
{
#if defined(__CINT__) && !defined(__MAKECINT__)
  gSystem->Load("libITSUpgradeBase.so");
  gSystem->Load("libITSUpgradeSim.so");
  gSystem->Load("libITSUpgradeRec.so");
  gROOT->ProcessLine(".L FT2.cxx+");
  FT2* det=0;
#endif
  //
  det = new FT2();
  det->InitEnvLocal();
  det->InitDetector();
  det->SetSimMat(kTRUE);
  det->SetMaxStepTGeo(1.);
  det->SetdNdY(dndy);
  det->SetUseKalmanOut(useKalmanOut);
  //
  hdca  = new TH1F("hdca","dca",  100,-0.1,.1);
  hdcaN = new TH1F("hdcaN","dcaN",100,-10.,10.);
  hdcaZ  = new TH1F("hdcaZ","dcaZ",  100,-0.1,.1);
  hdcaZN = new TH1F("hdcaZN","dcaZN",100,-10.,10.);
  //
  hFdca  = new TH1F("hFdca","dca fake",  100,-0.1,.1);
  hFdcaN = new TH1F("hFdcaN","dcaN fake",100,-10.,10.);
  hFdcaZ  = new TH1F("hFdcaZ","dcaZ fake",  100,-0.1,.1);
  hFdcaZN = new TH1F("hFdcaZN","dcaZN fake",100,-10.,10.);
  //
  hPatternITS  = new TH1F("itsPattern","ITS hits pattern"      ,7,-0.5,6.5);
  hFPatternITS = new TH1F("itsFPattern","ITS fake hits pattern",7,-0.5,6.5);
  //
  AliESDVertex *vtx = new AliESDVertex();
  double vcov[6] = {1e-4, 0, 1e-4, 0, 0, 1e-4};
  vtx->SetCovarianceMatrix(vcov);
  vtx->Print();
  //
  TParticle prt;
  double p = 0.25;
  for (int ntr=0;ntr<ntrials;ntr++) {
    
    vtx->SetXv(gRandom->Gaus(-0.1, 50e-4));
    vtx->SetYv(gRandom->Gaus(0.2,  50e-4));
    vtx->SetZv(gRandom->Gaus(0.5, 5.0));
    //
    double phi = gRandom->Rndm()*TMath::TwoPi();
    double eta =  2*(gRandom->Rndm()-0.5)*0.8;
    double theta = 2*TMath::ATan(TMath::Exp(-eta));
    double pz = p*TMath::Cos(theta);
    double pt = p*TMath::Sin(theta);
    double pxyz[3]={pt*TMath::Cos(phi),pt*TMath::Sin(phi),pz};
    double en = TMath::Sqrt(p*p+0.14*0.14);
    prt.SetPdgCode(gRandom->Rndm()>0.5 ? 211 : -211);
    //prt.SetPdgCode(-211);
    prt.SetMomentum(pxyz[0],pxyz[1],pxyz[2],en);
    prt.SetProductionVertex(vtx->GetX(),vtx->GetY(),vtx->GetZ(),0);
    //
    if (det->ProcessTrack(&prt,vtx)) {
      //const AliExternalTrackParam& prob = det->GetProbe();
      //      printf("%d %d\n",det->GetNClITS(),det->GetNClTPC());
      //    prob.Print();
      //      /*
      const double* dca = det->GetDCA();
      const double* cov = det->GetDCACov();
      hdca->Fill(dca[0]);
      hdcaN->Fill(dca[0]/TMath::Sqrt(cov[0]));
      hdcaZ->Fill(dca[1]);
      hdcaZN->Fill(dca[1]/TMath::Sqrt(cov[2]));
      //
      if (det->GetNClITSFakes()) {
	hFdca->Fill(dca[0]);
	hFdcaN->Fill(dca[0]/TMath::Sqrt(cov[0]));
	hFdcaZ->Fill(dca[1]);
	hFdcaZN->Fill(dca[1]/TMath::Sqrt(cov[2]));
      }
      int hits  = det->GetITSPattern();
      int hitsF = det->GetITSPatternFakes();
      for (int j=7;j--;) {
	if (hits&(0x1<<j))  hPatternITS->Fill(j);
	if (hitsF&(0x1<<j)) hFPatternITS->Fill(j);	
      }
    }
    else {
      printf("Failed on track %d\n",ntr);
    }
  }
  //
}
