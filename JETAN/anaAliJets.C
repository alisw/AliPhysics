// $Id$

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <iostream.h>
#include <time.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TGraph.h>
#include <TNtuple.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TFile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TTimeStamp.h>
#include "AliTkConeJetEvent.h"
#include "AliTkConeJet.h"
#include "AliTkConeJetFinderV2.h"
#include <AliJetParticle.h>
#include <AliJetParticlesReader.h>
#include <AliJetParticlesReaderKine.h>
#include <AliJetParticlesReaderESD.h>
#include <AliJetParticlesReaderHLT.h>
#include <AliJetEventParticles.h>
#endif

Float_t relphi(Float_t phi1, Float_t phi2);
Float_t addphi(Float_t phi1, Float_t phi2);
Float_t diffphi(Float_t phi1, Float_t phi2);
Int_t eventindex(Float_t phi1, Float_t phi2);
void convert(Float_t pjet[4], Float_t &pt, Float_t &theta, Float_t &eta, Float_t &phi);

void anaAliJets(Char_t *filename,Char_t *savefilename,
		Int_t mEnergy=100,Int_t mBackEn=-1,Int_t nMaxEvents=-1,
		Char_t *evfilename=0,Char_t *sigevfilename=0,
		Char_t *signalfilename=0,Char_t *monteconefilename=0)
  /*
    filename       = cone finder event
    evfilename     = background event or jetevent only (if 0 take from jetevent)
    sigevfilename  = signal event if background event is used (otherwice == 0)
    signalfilename = original signal event (eg. for real tracking to get trigger jets)
    monteconefilename = reconstructed jets from signal event (esd=0,esd=10)
   */

{
  const Float_t minEtRecoCut=mEnergy*0.95;
  const Float_t maxEtRecoCut=mEnergy*1.05;
  Float_t minJetEt;
  if(mBackEn) minJetEt=mBackEn;
  else minJetEt=mEnergy/4.;
  const Float_t minPartPt=0.5;
  //const Char_t *figprefix=0;
  //const Char_t *figdirname=".";

  const Int_t nclasses=22;
  const Float_t cletmin[nclasses] = {0,10,20,30,40,50,60,70,80,90,100,
				     110,120,130,140,150,160,170,180,190,200,0};
  const Float_t cletmax[nclasses] = {10,20,30,40,50,60,70,80,90,100,
				     110,120,130,140,150,160,170,180,190,200,350,350};

  //differential shape
  const Float_t deltaR=0.1/2;

  Float_t corrfac=0.;
#ifdef APPLYCORRECTION
  if(!allparts) corrfac=2./3;
  Float_t conefluc=1.;
  if(Radius<=0.3) conefluc=0.8;
  else if(Radius<=0.5) conefluc=0.9;
  else if(Radius<=0.7) conefluc=0.9;
  corrfac*=conefluc;
#endif

  Char_t dummy[1024];
  Char_t name[1024];
  Char_t title[1024];

  TH1F *hJetEt = new TH1F("hJetEt","E_{T} (jet)",350,0,350);
  hJetEt->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEt->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");

  TH1F *hBackJetEt = new TH1F("hBackJetEt","E_{T} (jet)",350,0,350);
  hBackJetEt->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hBackJetEt->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");

  TH1F *hJetEtall = new TH1F("hAllJetEt","E_{T} (jet)",350,0,350);
  hJetEtall->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEtall->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");

  TProfile *hJetEttrue = new TProfile("hJetEttrue","E_{T} (jet)",350,0,350,0,1.5);
  hJetEttrue->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEttrue->GetYaxis()->SetTitle("E^{true}_{T}/E_{T}");

  TProfile *hBackJetEttrue = new TProfile("hBackJetEttrue","E_{T} (jet)",350,0,350,0,1.5);
  hBackJetEttrue->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hBackJetEttrue->GetYaxis()->SetTitle("E^{true}_{T}/E_{T}");

  TProfile *hJetEtalltrue = new TProfile("hAllJetEttrue","E_{T} (jet)",350,0,350,0,1.5);
  hJetEtalltrue->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEtalltrue->GetYaxis()->SetTitle("E^{true}_{T}/E_{T}");

  TH1F *hJetEtUQTrigger = new TH1F("hJetEtUQTrigger","E_{T} (trigger jet)",350,0,350);
  hJetEtUQTrigger->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEtUQTrigger->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");

  TH1F *hJetEtTrigger = new TH1F("hJetEtTrigger","E_{T} (trigger jet)",350,0,350);
  hJetEtTrigger->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEtTrigger->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");

  TH2F *hJetEtvsEll = new TH2F("hJetEtvsEll","E_{T} (jet) versus Length",350,0,350,100,0,10);
  hJetEtvsEll->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEtvsEll->GetYaxis()->SetTitle("L [fm]");

  TH2F *hJetEtallvsEll = new TH2F("hAllJetEtallvsEll","E_{T} (jet) versus Length",350,0,350,100,0,10);
  hJetEtallvsEll->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEtallvsEll->GetYaxis()->SetTitle("L [fm]");

  TH2F *hJetEtvsTrigger = new TH2F("hJetEtvsTrigger","",350,0,350,350,0,350);
  hJetEtvsTrigger->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEtvsTrigger->GetYaxis()->SetTitle("E_{T} (jettrigger) [GeV]");

  TH2F *hJetEtvsUQTrigger = new TH2F("hJetEtvsUQTrigger","",350,0,350,350,0,350);
  hJetEtvsUQTrigger->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEtvsUQTrigger->GetYaxis()->SetTitle("E_{T} (jettrigger) [GeV]");

  TH2F *hJetEtvsEt = new TH2F("hJetEtvsEt","",350,0,350,350,0,350);
  hJetEtvsEt->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEtvsEt->GetYaxis()->SetTitle("E_{T} (jet) [GeV]");

  TH1F *hJetZ     =  new TH1F("hjetZ","Z distribution",100,0,1);
  hJetZ->GetXaxis()->SetTitle("Z");
  hJetZ->GetYaxis()->SetTitle("dN/dZ");

  TH1F *hJet1    =  new TH1F("hjet1","E_{t} distribution",350,0,350);
  hJet1->GetXaxis()->SetTitle("E_{T} (parton) [GeV]");
  hJet1->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");

  TH1F *hJet2    =  new TH1F("hjet2","E_{t} distribution",350,0,350);
  hJet2->GetXaxis()->SetTitle("E_{T} (parton) [GeV]");
  hJet2->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");

  TH1F **hJettype=  new TH1F*[3];
  for(Int_t i=0;i<3;i++){
    Char_t t[100];
    sprintf(t,"hJettype%d",i);
    Char_t tit[100];
    if(i==0)
      sprintf(tit,"Unknown");
    else if(i==1)
      sprintf(tit,"Quark");
    else 
      sprintf(tit,"Gluon");
    hJettype[i]=new TH1F(t,tit,350,0,350);
    hJettype[i]->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
    hJettype[i]->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");
  }

  TH2F *hAxesDiff = new TH2F("hAxesDiff","",120,0,TMath::Pi(),40,0,2);
  hAxesDiff->GetXaxis()->SetTitle("#Delta #phi");
  hAxesDiff->GetYaxis()->SetTitle("#Delta #eta");

  //---
  TH1F *hJetEtres = new TH1F("hJetEtres","E_{T} (jet)",350,0,350);
  hJetEtres->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEtres->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");

  TH1F *hJetEtratio = new TH1F("hJetEtratio","E_{T} (jet)",200,0,2);
  hJetEtratio->GetXaxis()->SetTitle("E_{T} (jet)/E_{T} (trigger)");
  hJetEtratio->GetYaxis()->SetTitle("Number of jets");

  TProfile *hJetEtrestrue = new TProfile("hJetEttestrue","E_{T} (jet)",35,0,350,-2,2);
  hJetEtrestrue->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEtrestrue->GetYaxis()->SetTitle("E^{true}_{T}/E_{T}");

  TH1F *hJetEtresTrigger = new TH1F("hJetEtresTrigger","E_{T} (trigger jet)",350,0,350);
  hJetEtresTrigger->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEtresTrigger->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");

  TProfile2D *hAxesDiffres = new TProfile2D("hAxesDiffres","",240,0,TMath::Pi(),100,0,2);
  hAxesDiffres->GetXaxis()->SetTitle("#Delta #phi");
  hAxesDiffres->GetYaxis()->SetTitle("#Delta #eta");

  TH1F *hPhires = new TH1F("hPhires","",250,-1,1);
  hPhires->GetXaxis()->SetTitle("#Delta #phi [rad]");
  hPhires->GetYaxis()->SetTitle("Number of jets");

  TH1F *hEtares = new TH1F("hEtares","",250,-1,1);
  hEtares->GetXaxis()->SetTitle("#Delta #eta [rad]");
  hEtares->GetYaxis()->SetTitle("Number of jets");

  TH1F *hmJetEtres = new TH1F("hmJetEtres","E_{T} (jet)",350,0,350);
  hmJetEtres->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hmJetEtres->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");

  TH1F *hmJetEtratio = new TH1F("hmJetEtratio","E_{T} (jet)",200,0,2);
  hmJetEtratio->GetXaxis()->SetTitle("E_{T} (jet)/E_{T} (trigger)");
  hmJetEtratio->GetYaxis()->SetTitle("Number of jets");

  TProfile *hmJetEtrestrue = new TProfile("hmJetEttestrue","E_{T} (jet)",35,0,350,-2,2);
  hmJetEtrestrue->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hmJetEtrestrue->GetYaxis()->SetTitle("E^{true}_{T}/E_{T}");

  TProfile2D *hmAxesDiffres = new TProfile2D("hmAxesDiffres","",240,0,TMath::Pi(),100,0,2);
  hmAxesDiffres->GetXaxis()->SetTitle("#Delta #phi");
  hmAxesDiffres->GetYaxis()->SetTitle("#Delta #eta");

  TH1F *hmPhires = new TH1F("hmPhires","",250,-1,1);
  hmPhires->GetXaxis()->SetTitle("#Delta #phi [rad]");
  hmPhires->GetYaxis()->SetTitle("Number of jets");

  TH1F *hmEtares = new TH1F("hmEtares","",250,-1,1);
  hmEtares->GetXaxis()->SetTitle("#Delta #eta [rad]");
  hmEtares->GetYaxis()->SetTitle("Number of jets");

  TH1F *hPhiMonteres = new TH1F("hPhiMonteres","",250,-1,1);
  hPhiMonteres->GetXaxis()->SetTitle("#Delta #phi [rad]");
  hPhiMonteres->GetYaxis()->SetTitle("Number of jets");

  TH1F *hEtaMonteres = new TH1F("hEtaMonteres","",250,-1,1);
  hEtaMonteres->GetXaxis()->SetTitle("#Delta #eta [rad]");
  hEtaMonteres->GetYaxis()->SetTitle("Number of jets");

  TH1F *hEtMonteres = new TH1F("hEtMonteres","",250,-125,125);
  hEtMonteres->GetXaxis()->SetTitle("#Delta E_{T} [GeV]");
  hEtMonteres->GetYaxis()->SetTitle("Number of jets");

  TH1F *hEtMonteratio = new TH1F("hEtMonteratio","E_{T} (jet)",200,0,2);
  hEtMonteratio->GetXaxis()->SetTitle("E^{rec}_{T} / E^{monte}_{T}");
  hEtMonteratio->GetYaxis()->SetTitle("Number of jets");

  //---

  TH1F **hJetEttrigger=new TH1F*[9];
  TH1F **hJetEttrigger2=new TH1F*[9];
  for(Int_t i=0;i<9;i++){
    sprintf(dummy,"hJetEttrigger%d",i);
    hJetEttrigger[i] = new TH1F(dummy,"E_{T} (jet)",350,0,350);
    hJetEttrigger[i]->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
    hJetEttrigger[i]->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");
    sprintf(dummy,"hJetEttrigger2%d",i);
    hJetEttrigger2[i] = new TH1F(dummy,"E_{T} (jet)",350,0,350);
    hJetEttrigger2[i]->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
    hJetEttrigger2[i]->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");
  }
  sprintf(dummy,"hJetEttriggernorm");
  TH1F *hJetEttriggernorm = new TH1F(dummy,"E_{T} (jet)",350,0,350);
  hJetEttriggernorm->GetXaxis()->SetTitle("E_{T} (jet) [GeV]");
  hJetEttriggernorm->GetYaxis()->SetTitle("dN/dE_{T} [GeV^{-1}]");

  //---

  TH1F **hJetLeadingPt = new TH1F*[nclasses];
  for(Int_t i=0;i<nclasses;i++){
    sprintf(dummy,"hJetLeadingPt%d",i);
    sprintf(title,"Transverse Momentum Fraction of Leading Particle (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
    hJetLeadingPt[i]=new TH1F(dummy,title,100,0,1);
    hJetLeadingPt[i]->GetXaxis()->SetTitle("z=p_{T}^{lead.part.}/P_{T}^{jet}");
    hJetLeadingPt[i]->GetYaxis()->SetTitle("dN/dz");
  }

  TH1F **hJetFragLeadingPt = new TH1F*[nclasses];
  for(Int_t i=0;i<nclasses;i++){
    sprintf(dummy,"hJetFragLeadingPt%d",i);
    sprintf(title,"Fragmentation Function (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
    hJetFragLeadingPt[i]=new TH1F(dummy,title,100,0,1);
    hJetFragLeadingPt[i]->GetXaxis()->SetTitle("z=p_{T}/p_{T}^{lead.part.}");
    hJetFragLeadingPt[i]->GetYaxis()->SetTitle("dN/dz");
  }

  TH1F **hJetLeadingPtDist = new TH1F*[nclasses];
  for(Int_t i=0;i<nclasses;i++){
    sprintf(dummy,"hJetLeadingPtDist%d",i);
    sprintf(title,"Transverse Momentum Fraction of Leading Particle (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
    hJetLeadingPtDist[i]=new TH1F(dummy,title,600,0,150);
    hJetLeadingPtDist[i]->GetXaxis()->SetTitle("p_{T}^{lead.part.} [GeV]");
    hJetLeadingPtDist[i]->GetYaxis()->SetTitle("dN/dp_{T}^{lead.part.} [GeV^{-1}]");
  }

  TH1F **hJetFragL = new TH1F*[nclasses];
  for(Int_t i=0;i<nclasses;i++){
    sprintf(dummy,"hJetFragL%d",i);
    sprintf(title,"Longitudinal Fragmentation Function (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
    hJetFragL[i]=new TH1F(dummy,title,100,0,1);
    hJetFragL[i]->GetXaxis()->SetTitle("z=p_{L}/P^{jet}");
    hJetFragL[i]->GetYaxis()->SetTitle("dN/dz");
  }

  TH1F **hJetFragPL = new TH1F*[nclasses];
  for(Int_t i=0;i<nclasses;i++){
    sprintf(dummy,"hJetFragPL%d",i);
    sprintf(title,"Longitudinal Momentum (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
    hJetFragPL[i]=new TH1F(dummy,title,600,0,150);
    hJetFragPL[i]->GetXaxis()->SetTitle("p_{L}");
    hJetFragPL[i]->GetYaxis()->SetTitle("dN/dp_{L} [GeV^{-1}]");
  }

  TH1F **hJetFragT = new TH1F*[nclasses];
  for(Int_t i=0;i<nclasses;i++){
    sprintf(dummy,"hJetFragT%d",i);
    sprintf(title,"Transversal Fragmentation Function (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
    hJetFragT[i]=new TH1F(dummy,title,250,0,25);
    hJetFragT[i]->GetXaxis()->SetTitle("j_{T} [GeV]");
    hJetFragT[i]->GetYaxis()->SetTitle("dN/dj_{T}");
  }

  TH1F **hJetFragPt = new TH1F*[nclasses];
  for(Int_t i=0;i<nclasses;i++){
    sprintf(dummy,"hJetFragPt%d",i);
    sprintf(title,"Particle Transversal Momentum Distribution (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
    hJetFragPt[i]=new TH1F(dummy,title,600,0,150);
    hJetFragPt[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
    hJetFragPt[i]->GetYaxis()->SetTitle("dN/dp_{T}");
  }

  TH1F **hJetN = new TH1F*[nclasses];
  for(Int_t i=0;i<nclasses;i++){
    sprintf(dummy,"hJetN%d",i);
    sprintf(title,"Particle Multiplicity (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
    hJetN[i]=new TH1F(dummy,title,150+1,-0.5,150+0.5);
    hJetN[i]->GetXaxis()->SetTitle("n = Number of Particles within Jet");
    hJetN[i]->GetYaxis()->SetTitle("dN/dn");
  }

  TH1F **hJetMeanPt = new TH1F*[nclasses];
  for(Int_t i=0;i<nclasses;i++){
    sprintf(dummy,"hJetMeanPt%d",i);
    sprintf(title,"Mean Transverse Jet Energy per Particle (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
    hJetMeanPt[i]=new TH1F(dummy,title,125,0,25);
    hJetMeanPt[i]->GetXaxis()->SetTitle("#bar{p} = P_{T}^{jet}/n [GeV]");
    hJetMeanPt[i]->GetYaxis()->SetTitle("dN/d#bar{p}");
  }

  //intshape
  Float_t retall[nclasses][10];
  Float_t retlow[nclasses][10];
  Float_t retlow1[nclasses][10];
  Float_t retlow2[nclasses][10];
  Float_t retlow3[nclasses][10];
  Float_t retlow4[nclasses][10];
  Float_t rxet[10];
  for(Int_t i=0;i<10;i++) {
    for(Int_t j=0;j<nclasses;j++) {
      retall[j][i]=0.;
      retlow[j][i]=0.;
      retlow1[j][i]=0.;
      retlow2[j][i]=0.;
      retlow3[j][i]=0.;
      retlow4[j][i]=0.;
    }
    rxet[i]=(i+1.)/10;
  }

  //diffshape
  Float_t dretall[nclasses][11];
  Float_t dretlow[nclasses][11];
  Float_t dretlow1[nclasses][11];
  Float_t dretlow2[nclasses][11];
  Float_t dretlow3[nclasses][11];
  Float_t dretlow4[nclasses][11];
  Float_t drxet[11];
  for(Int_t i=0;i<11;i++) {
    for(Int_t j=0;j<nclasses;j++) {
      dretall[j][i]=0.;
      dretlow[j][i]=0.;
      dretlow1[j][i]=0.;
      dretlow2[j][i]=0.;
      dretlow3[j][i]=0.;
      dretlow4[j][i]=0.;
    }
    drxet[i]=1.*i/10.;
  }

  TH1F **hPhiCorr = new TH1F*[nclasses];
  for(Int_t i=0;i<nclasses;i++){
    sprintf(dummy,"hPartPhiCorr%d",i);
    sprintf(title,"Azimuthal Correlation (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
    hPhiCorr[i]= new TH1F(dummy,title,120,0,TMath::TwoPi());
    hPhiCorr[i]->GetXaxis()->SetTitle("#phi");
    hPhiCorr[i]->GetYaxis()->SetTitle("dN/d#phi");
  }

  //
  //global event properties (ue)
  // [0]=toward, [1]=away, [2]=transverse,
  //
  AliTkEtaPhiVector **centers=new AliTkEtaPhiVector*[4];

  TH1F ***hgJetFragLeadingPt = new TH1F**[3];
  TH1F ***hgJetFragL = new TH1F**[3];
  TH1F ***hgJetFragPL = new TH1F**[3];
  TH1F ***hgJetFragT = new TH1F**[3];
  TH1F ***hgJetFragPt = new TH1F**[3];
  TH1F ***hgDiJetFragLeadingPt = new TH1F**[3];
  TH1F ***hgDiJetFragL = new TH1F**[3];
  TH1F ***hgDiJetFragPL = new TH1F**[3];
  TH1F ***hgDiJetFragT = new TH1F**[3];
  TH1F ***hgDiJetFragPt = new TH1F**[3];
  for(Int_t k=0;k<3;k++){
    hgJetFragLeadingPt[k] = new TH1F*[nclasses];
    hgJetFragL[k] = new TH1F*[nclasses];
    hgJetFragPL[k] = new TH1F*[nclasses];
    hgJetFragT[k] = new TH1F*[nclasses];
    hgJetFragPt[k] = new TH1F*[nclasses];
    hgDiJetFragLeadingPt[k] = new TH1F*[nclasses];
    hgDiJetFragL[k] = new TH1F*[nclasses];
    hgDiJetFragPL[k] = new TH1F*[nclasses];
    hgDiJetFragT[k] = new TH1F*[nclasses];
    hgDiJetFragPt[k] = new TH1F*[nclasses];
    if(k==0){
      sprintf(name,"toward");
    } else if (k==1) {
      sprintf(name,"away");
    } else {
      sprintf(name,"transverse");
    }
    for(Int_t i=0;i<nclasses;i++){
      sprintf(dummy,"h%s-JetFragLeadingPt%d",name,i);
      sprintf(title,"Fragmentation Function (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
      hgJetFragLeadingPt[k][i]=new TH1F(dummy,title,100,0,1);
      hgJetFragLeadingPt[k][i]->GetXaxis()->SetTitle("z=p_{T}/p_{T}^{lead.part.}");
      hgJetFragLeadingPt[k][i]->GetYaxis()->SetTitle("dN/dz");

      sprintf(dummy,"h%s-JetFragL%d",name,i);
      sprintf(title,"Longitudinal Fragmentation Function (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
      hgJetFragL[k][i]=new TH1F(dummy,title,100,0,1);
      hgJetFragL[k][i]->GetXaxis()->SetTitle("z=p_{L}/P^{jet}");
      hgJetFragL[k][i]->GetYaxis()->SetTitle("dN/dz");

      sprintf(dummy,"h%s-JetFragPL%d",name,i);
      sprintf(title,"Longitudinal Momentum (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
      hgJetFragPL[k][i]=new TH1F(dummy,title,600,0,150);
      hgJetFragPL[k][i]->GetXaxis()->SetTitle("p_{L}");
      hgJetFragPL[k][i]->GetYaxis()->SetTitle("dN/dp_{L} [GeV^{-1}]");

      sprintf(dummy,"h%s-JetFragT%d",name,i);
      sprintf(title,"Transversal Fragmentation Function (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
      hgJetFragT[k][i]=new TH1F(dummy,title,250,0,25);
      hgJetFragT[k][i]->GetXaxis()->SetTitle("j_{T} [GeV]");
      hgJetFragT[k][i]->GetYaxis()->SetTitle("dN/dj_{T}");

      sprintf(dummy,"h%s-JetFragPt%d",name,i);
      sprintf(title,"Particle Transversal Momentum Distribution (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
      hgJetFragPt[k][i]=new TH1F(dummy,title,600,0,150);
      hgJetFragPt[k][i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      hgJetFragPt[k][i]->GetYaxis()->SetTitle("dN/dp_{T}");

      sprintf(dummy,"h%s-DiJetFragLeadingPt%d",name,i);
      sprintf(title,"Fragmentation Function (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
      hgDiJetFragLeadingPt[k][i]=new TH1F(dummy,title,100,0,1);
      hgDiJetFragLeadingPt[k][i]->GetXaxis()->SetTitle("z=p_{T}/p_{T}^{lead.part.}");
      hgDiJetFragLeadingPt[k][i]->GetYaxis()->SetTitle("dN/dz");

      sprintf(dummy,"h%s-DiJetFragL%d",name,i);
      sprintf(title,"Longitudinal Fragmentation Function (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
      hgDiJetFragL[k][i]=new TH1F(dummy,title,100,0,1);
      hgDiJetFragL[k][i]->GetXaxis()->SetTitle("z=p_{L}/P^{jet}");
      hgDiJetFragL[k][i]->GetYaxis()->SetTitle("dN/dz");

      sprintf(dummy,"h%s-DiJetFragPL%d",name,i);
      sprintf(title,"Longitudinal Momentum (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
      hgDiJetFragPL[k][i]=new TH1F(dummy,title,600,0,150);
      hgDiJetFragPL[k][i]->GetXaxis()->SetTitle("p_{L}");
      hgDiJetFragPL[k][i]->GetYaxis()->SetTitle("dN/dp_{L} [GeV^{-1}]");

      sprintf(dummy,"h%s-DiJetFragT%d",name,i);
      sprintf(title,"Transversal Fragmentation Function (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
      hgDiJetFragT[k][i]=new TH1F(dummy,title,250,0,25);
      hgDiJetFragT[k][i]->GetXaxis()->SetTitle("j_{T} [GeV]");
      hgDiJetFragT[k][i]->GetYaxis()->SetTitle("dN/dj_{T}");

      sprintf(dummy,"h%s-DiJetFragPt%d",name,i);
      sprintf(title,"Particle Transversal Momentum Distribution (%.0f GeV<E_{T}<%.0f GeV)",cletmin[i],cletmax[i]);
      hgDiJetFragPt[k][i]=new TH1F(dummy,title,600,0,150);
      hgDiJetFragPt[k][i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      hgDiJetFragPt[k][i]->GetYaxis()->SetTitle("dN/dp_{T}");
    }
  }

  //intshape 
  Float_t gretall[3][nclasses][10];
  Float_t gretlow[3][nclasses][10];
  Float_t gretlow1[3][nclasses][10];
  Float_t gretlow2[3][nclasses][10];
  Float_t gretlow3[3][nclasses][10];
  Float_t gretlow4[3][nclasses][10];
  Float_t grxet[10];
  for(Int_t i=0;i<10;i++) {
    for(Int_t j=0;j<nclasses;j++) {
      for(Int_t k=0;k<3;k++){
	gretall[k][j][i]=0.;
	gretlow[k][j][i]=0.;
	gretlow1[k][j][i]=0.;
	gretlow2[k][j][i]=0.;
	gretlow3[k][j][i]=0.;
	gretlow4[k][j][i]=0.;
      }
    }
    grxet[i]=(i+1.)/10;
  }

  //diffshape
  Float_t gdretall[3][nclasses][11];
  Float_t gdretlow[3][nclasses][11];
  Float_t gdretlow1[3][nclasses][11];
  Float_t gdretlow2[3][nclasses][11];
  Float_t gdretlow3[3][nclasses][11];
  Float_t gdretlow4[3][nclasses][11];
  Float_t gdrxet[11];
  for(Int_t i=0;i<11;i++) {
    for(Int_t j=0;j<nclasses;j++) {
      for(Int_t k=0;k<3;k++){
	gdretall[k][j][i]=0.;
	gdretlow[k][j][i]=0.;
	gdretlow1[k][j][i]=0.;
	gdretlow2[k][j][i]=0.;
	gdretlow3[k][j][i]=0.;
	gdretlow4[k][j][i]=0.;
      }
    }
    gdrxet[i]=1.*i/10.;
  }

  //intshape 
  Float_t gdiretall[3][nclasses][10];
  Float_t gdiretlow[3][nclasses][10];
  Float_t gdiretlow1[3][nclasses][10];
  Float_t gdiretlow2[3][nclasses][10];
  Float_t gdiretlow3[3][nclasses][10];
  Float_t gdiretlow4[3][nclasses][10];
  for(Int_t i=0;i<10;i++) {
    for(Int_t j=0;j<nclasses;j++) {
      for(Int_t k=0;k<3;k++){
	gdiretall[k][j][i]=0.;
	gdiretlow[k][j][i]=0.;
	gdiretlow1[k][j][i]=0.;
	gdiretlow2[k][j][i]=0.;
	gdiretlow3[k][j][i]=0.;
	gdiretlow4[k][j][i]=0.;
      }
    }
  }

  //diffshape
  Float_t gdidretall[3][nclasses][11];
  Float_t gdidretlow[3][nclasses][11];
  Float_t gdidretlow1[3][nclasses][11];
  Float_t gdidretlow2[3][nclasses][11];
  Float_t gdidretlow3[3][nclasses][11];
  Float_t gdidretlow4[3][nclasses][11];
  for(Int_t i=0;i<11;i++) {
    for(Int_t j=0;j<nclasses;j++) {
      for(Int_t k=0;k<3;k++){
	gdidretall[k][j][i]=0.;
	gdidretlow[k][j][i]=0.;
	gdidretlow1[k][j][i]=0.;
	gdidretlow2[k][j][i]=0.;
	gdidretlow3[k][j][i]=0.;
	gdidretlow4[k][j][i]=0.;
      }
    }
  }

  TH1F *hPartPtDist = new TH1F("hPartPtDist","Transverse Momentum Distribution",600,0,150);
  hPartPtDist->GetXaxis()->SetTitle("p_{T} [GeV]");
  hPartPtDist->GetYaxis()->SetTitle("dN/dp_{T} [GeV^{-1}]");

  TH1F *hPartPhiDist = new TH1F("hPartPhiDist","Azimuthal Distribution",120,0,TMath::TwoPi());
  hPartPhiDist->GetXaxis()->SetTitle("#phi");
  hPartPhiDist->GetYaxis()->SetTitle("dN/d#phi");

  TH1F *hPartEtaDist = new TH1F("hPartEtaDist","Pseudo-Rapidity Distribution",100,-1,1);
  hPartEtaDist->GetXaxis()->SetTitle("#eta");
  hPartEtaDist->GetYaxis()->SetTitle("dN/d#eta");

  TH1F *hPartPhiCorr = new TH1F("hPartPhiCorr","Azimuthal Correlation",120,0,TMath::TwoPi());
  hPartPhiCorr->GetXaxis()->SetTitle("#phi");
  hPartPhiCorr->GetYaxis()->SetTitle("dN/d#phi");

  TH1F *hPartDiPhiCorr = new TH1F("hPartDiPhiCorr","Azimuthal Correlation",120,0,TMath::TwoPi());
  hPartDiPhiCorr->GetXaxis()->SetTitle("#phi");
  hPartDiPhiCorr->GetYaxis()->SetTitle("dN/d#phi");

  TH1F *hPartACorr = new TH1F("hPartACorr","Alpha Correlation",100,-1.1,1.1);
  hPartACorr->GetXaxis()->SetTitle("#alpha");
  hPartACorr->GetYaxis()->SetTitle("dN/d#alpha");

  TH1F *hPartDiACorr = new TH1F("hPartDiACorr","Alpha Correlation",100,-1.1,1.1);
  hPartDiACorr->GetXaxis()->SetTitle("#alpha");
  hPartDiACorr->GetYaxis()->SetTitle("dN/d#alpha");
  
  Int_t nTotalEvents = 0;
  Int_t nclGoodEvents[nclasses];
  for(Int_t i=0;i<nclasses;i++) nclGoodEvents[i]=0;
  Int_t nclLeadEvents[nclasses];
  for(Int_t i=0;i<nclasses;i++) nclLeadEvents[i]=0;
  Int_t nclDiEvents[nclasses];
  for(Int_t i=0;i<nclasses;i++) nclDiEvents[i]=0;

  //connect to jets
  TChain *theTree = new TChain("jets");
  theTree->Add(filename);
  AliTkConeJetEvent *event = new AliTkConeJetEvent();
  theTree->SetBranchAddress("ConeFinder",&event);

  Int_t treeentries=(Int_t)theTree->GetEntries();
  if((nMaxEvents<0) || (nMaxEvents>treeentries))
    nMaxEvents=treeentries;

  cout << "Found " << nMaxEvents << " in " << filename << endl;

  TChain *theEvTree=0;
  AliJetEventParticles *jetevent=0;
  TChain *theEvTree_sig=0;
  AliJetEventParticles *jetevent_sig=0;
  Int_t backtreeentries=0;
  Int_t nPerBackground=0;
  if(nPerBackground==0) nPerBackground=1;

  if(evfilename){
    theEvTree = new TChain("AJEPtree");
    theEvTree->Add(evfilename);
    jetevent=new AliJetEventParticles();
    theEvTree->SetBranchAddress("particles",&jetevent);
    Int_t evtreeentries=(Int_t)theEvTree->GetEntries();

    if(sigevfilename){
      theEvTree_sig = new TChain("AJEPtree");
      theEvTree_sig->Add(sigevfilename);
      jetevent_sig=new AliJetEventParticles();
      theEvTree_sig->SetBranchAddress("particles",&jetevent_sig);
      evtreeentries=(Int_t)theEvTree_sig->GetEntries();
      backtreeentries=(Int_t)theEvTree->GetEntries();
      nPerBackground=nMaxEvents/backtreeentries;
      if(nPerBackground==0) nPerBackground=1;
    }

    if(evtreeentries!=treeentries){
      cerr << "WARNING: ";
      cerr << "Total jet event number not equals event number: " 
	   << evtreeentries << " " << treeentries <<endl; 
      exit(1);
    }
  }

  TChain *theSignalTree=0;
  AliJetEventParticles *jetsigevent=0;
  if(signalfilename){
    theSignalTree = new TChain("AJEPtree");
    theSignalTree->Add(signalfilename);
    jetsigevent=new AliJetEventParticles();
    theSignalTree->SetBranchAddress("particles",&jetsigevent);

    Int_t sigtreeentries=(Int_t)theSignalTree->GetEntries();
    if(sigtreeentries!=treeentries){
      cerr << "WARNING: ";
      cerr << "Total jet signal event number not equals event number: " 
	   << sigtreeentries << " " << treeentries <<endl; 
      //exit(1);
    }
  }

  //connect to monte jets (if wanted)
  TChain *theTreeMonte=0;
  AliTkConeJetEvent *monteevent = 0;
  if(monteconefilename){
    theTreeMonte = new TChain("jets");
    theTreeMonte->Add(monteconefilename);
    monteevent= new AliTkConeJetEvent();
    theTreeMonte->SetBranchAddress("ConeFinder",&monteevent);
    Int_t montetreeentries=(Int_t)theTreeMonte->GetEntries();
    if(montetreeentries!=treeentries){
      cerr << "WARNING: ";
      cerr << "Total monte jet event number not equals jet event number: " 
	   << montetreeentries << " " << treeentries <<endl; 
      //exit(1);
      if(theSignalTree){
	Int_t sigtreeentries=(Int_t)theSignalTree->GetEntries();
	if(montetreeentries!=sigtreeentries){
	  cerr << "WARNING: ";
	  cerr << "Total signal cone jet event number not equals signal event number: " 
	       << sigtreeentries << " " << montetreeentries <<endl; 
	  exit(1);
	}
      }
    }
  }

  //=========================================================================
  // start the event loop
  //=========================================================================
  Int_t nEvent = -1;
  Int_t nEventSig = -1;
  Int_t nEventHijing = -1;
  Int_t nEventHijingCounter = nPerBackground;
  while(nEvent<nMaxEvents-1){
    nEvent++;
    nEventSig++;

    if ((nEvent % 100) == 0) {
      cout << "Analysing event " << nEvent << endl;
    }

    //connect the cone event/jets
    event->Clear();
    theTree->GetEvent(nEvent);
    event->sortJets();
    //event->Print("");
    Float_t ptcutused=event->getPtCut();
    if(ptcutused<minPartPt) ptcutused=minPartPt;

    TClonesArray *jets=event->getJets();

    if(theEvTree_sig){ // need to mix
      if(nEventHijingCounter==nPerBackground){    
	jetevent->Reset();
	nEventHijing++;
	theEvTree->GetEvent(nEventHijing);
	//jetevent->Print();
	if(nEventHijing==backtreeentries) nEventHijing=0;
	nEventHijingCounter=0;
      }
      jetevent_sig->Reset();
      theEvTree_sig->GetEvent(nEvent);
      //jetevent_sig->Print();
      jetevent->AddSignal(*jetevent_sig);
      TString dummy="Counter: ";
      dummy+=nEvent;
      dummy+=" ";
      dummy+="(Pythia ";dummy+=nEvent;
      dummy+=" ";
      dummy+=", Hijing ";dummy+=nEventHijing;
      dummy+=")";
      jetevent->SetHeader(dummy);
      nEventHijingCounter++;
    } else if(theEvTree){
      theEvTree->GetEvent(nEvent);
    }
    else{
      jetevent=event->getJetParticles();
    }

    if(theSignalTree){ //get MC event containing signal
      jetsigevent->Reset();
      theSignalTree->GetEvent(nEventSig);
#if 1
      if(jetsigevent->GetEventNr()!=jetevent->GetEventNr()){
	Int_t js=jetsigevent->GetEventNr();
	Int_t es=jetevent->GetEventNr();
	cerr << "Need to skip event: " << nEvent << ": " <<  es  << " != " << js << endl;
	if(es>js) nEventSig++;
	else nEvent++;
	continue;
      }
#endif
    }    
    else {
      jetsigevent=jetevent;
    }

    //connect the monte cone jets
    TClonesArray *montejets=0;
    if(theTreeMonte){
      monteevent->Clear();
      theTreeMonte->GetEvent(nEventSig);
      monteevent->sortJets();
      //monteevent->Print("");
      montejets=monteevent->getJets();
    }

    //
    //jetevent->Print();
    //

    Int_t njets=jets->GetEntries();
    if(!jets || !njets){
      cerr << "No Cone jet found in event " << nEvent << ", continuing..." <<endl;
      continue;
    }

    Int_t nhard=0;    //get hard partons
    TClonesArray *chard_jets=new TClonesArray("AliTkConeJet",0);
    Float_t phard[4];
    Float_t type;
    //TString header=jetsigevent->getHeader();
    //Int_t ntrials=jetsigevent->Trials();
    jetsigevent->Hard(0,phard,type);
    if(type!=0){
      Float_t ptj,thj,etaj,phj;
      convert(phard,ptj,thj,etaj,phj);
      new ((*chard_jets)[nhard]) AliTkConeJet(ptj,etaj,phj,(Int_t)type); 
      nhard++;
    }
    jetsigevent->Hard(1,phard,type);
    if(type!=0){
      Float_t ptj,thj,etaj,phj;
      convert(phard,ptj,thj,etaj,phj);
      new ((*chard_jets)[nhard]) AliTkConeJet(ptj,etaj,phj,(Int_t)type); 
      nhard++;
    }
    chard_jets->Sort();
    //chard_jets->Print();

    Int_t ntr=jetsigevent->NTriggerJets();
    Float_t tr_jets[ntr][3];    // trigger jets
    TClonesArray *ctr_jets=new TClonesArray("AliTkConeJet",ntr);
    Float_t lead_tr_pt=-1.,lead_tr_phi=0,lead_tr_eta=0.;
    for(Int_t j=0;j<ntr;j++){
      Float_t pjet[4];
      Float_t ptj,thj,etaj,phj;
      jetsigevent->TriggerJet(j,pjet);
      convert(pjet,ptj,thj,etaj,phj);
      tr_jets[j][0]=ptj;
      tr_jets[j][1]=phj;
      tr_jets[j][2]=etaj;
      new ((*ctr_jets)[j]) AliTkConeJet(ptj,etaj,phj); 
      if(lead_tr_pt<ptj) {
	 lead_tr_pt=ptj;
	 lead_tr_phi=phj;
	 lead_tr_eta=etaj;
      }
      for(Int_t i=0;i<nhard;i++){
	Float_t diff,etdiff,etadiff,phidiff;
	AliTkConeJet *jet=(AliTkConeJet*)ctr_jets->At(j);
	diff=AliTkConeJet::Diff(jet,(AliTkConeJet*)chard_jets->At(i),etdiff,phidiff,etadiff);
	//cout << diff << " " << etdiff << " " << phidiff << " " << etadiff << endl;
	if(diff<0.25) {
	  jet->setType((Int_t)type); 
	  break;
	}
      }
    }
    ctr_jets->Sort();
    //ctr_jets->Print();
    Int_t ntruq=jetsigevent->NUQTriggerJets();
    Float_t uq_jets[ntruq][3];  // unquenched jets
    TClonesArray *cuq_jets=new TClonesArray("AliTkConeJet",ntruq);
    Float_t lead_uq_pt=-1.,lead_uq_phi=0,lead_uq_eta=0.;
    for(Int_t j=0;j<ntruq;j++){
      Float_t pjet[4];
      Float_t ptj,thj,etaj,phj;
      jetsigevent->UQJet(j,pjet);
      convert(pjet,ptj,thj,etaj,phj);
      uq_jets[j][0]=ptj;
      uq_jets[j][1]=phj;
      uq_jets[j][2]=etaj;
      new ((*cuq_jets)[j]) AliTkConeJet(ptj,etaj,phj); 
      if(lead_uq_pt<ptj) {
	lead_uq_pt=ptj;
	lead_uq_phi=phj;
	lead_uq_eta=etaj;
      }
    }
    cuq_jets->Sort();
    //cuq_jets->Print();

    Double_t x0=jetsigevent->GetXJet();
    Double_t y0=jetsigevent->GetYJet();
    Double_t prodlength=TMath::Sqrt(x0*x0+y0*y0);
    Double_t zquench[4];
    jetsigevent->GetZQuench(zquench);

    //fill trigger histos (without cuts)
    hJetEtUQTrigger->Fill(lead_uq_pt);   
    hJetEtTrigger->Fill(lead_tr_pt);   
    for(Int_t i=0;i<2;i++){
      hJetZ->Fill(zquench[i]);
    }
    //reconstruction efficiency plots
    if(ntr) for(Int_t j=0;j<1/*ntr*/;j++){ //loop over MC jets
      Float_t diff,etdiff,etadiff,phidiff;
      Float_t mindiff=100;
      Int_t index=-1;
      AliTkConeJet *jt=(AliTkConeJet*)ctr_jets->At(j);
      if(jt->getEt()<minEtRecoCut||jt->getEt()>maxEtRecoCut) continue;
      //if(TMath:Abs(jt->getEta())>0.1) continue,
      //jets in event
      Int_t njets=jets->GetEntries();
      AliTkConeJet *jet=0;
      for(Int_t i=0;i<njets;i++){
	jet=(AliTkConeJet*)jets->At(i);
	if(!jet) continue;
	jet->calculateValues();
	diff=AliTkConeJet::Diff(jt,jet,etdiff,phidiff,etadiff);
	//_cut_ here if wanted (prob. for mixed events)
	//if(phidiff>0.25||etadiff>0.25) continue;
	//if(TMath::Abs(etdiff)/jt->getEt()>0.1) continue;
	//cout << diff << " " << etdiff << " " << phidiff << " " << etadiff << endl;
	if(mindiff>diff){
	  mindiff=diff;
	  index=i;
	}
      }
      if(index>-1) {
	jet=(AliTkConeJet*)jets->At(index);
	hJetEtres->Fill(jet->getEt());
	hJetEtresTrigger->Fill(jt->getEt());
	hJetEtratio->Fill(jet->getEt()/jt->getEt());
	Float_t phidiff=(jt->getPhi()-jet->getPhi());
	if(phidiff>TMath::Pi()) phidiff=TMath::TwoPi()-phidiff;
	Float_t etadiff=(jt->getEta()-jet->getEta());
	Float_t etfract=jet->getEtMarkedFrac();
	hAxesDiffres->Fill(phidiff,etadiff,TMath::Abs(jet->getEt()-jt->getEt()));
	hJetEtrestrue->Fill(jet->getEt(),etfract);
	hPhires->Fill(phidiff);
	hEtares->Fill(etadiff);
	//cout << "test" << nEvent << " " << phidiff << " " << etadiff << endl;
      }
    }    
    if(theTreeMonte){
      //compare leading jets from Monte and Reconstruction
      AliTkConeJet *jt=(AliTkConeJet*)jets->At(0);
      AliTkConeJet *mjt=(AliTkConeJet*)montejets->At(0);
      if(jt && mjt){
	Float_t phidiff=(jt->getPhi()-mjt->getPhi());
	if(phidiff>TMath::Pi()) phidiff=TMath::TwoPi()-phidiff;
	Float_t etadiff=(jt->getEta()-mjt->getEta());
	Float_t etdiff=jt->getEt()-mjt->getEt();
	hPhiMonteres->Fill(phidiff);
	hEtaMonteres->Fill(etadiff);
	hEtMonteres->Fill(etdiff);
	hEtMonteratio->Fill(jt->getEt()/mjt->getEt());
      }
      for(Int_t j=0;j<1/*montejets->Entries()*/;j++){
	Float_t diff,etdiff,etadiff,phidiff;
	Float_t mindiff=100;
	Int_t index=-1;
	AliTkConeJet *jt=(AliTkConeJet*)montejets->At(j);
	if(!jt) continue;
	//if(jt->getEt()<minEtRecoCut||jt->getEt()>maxEtRecoCut) continue;
	//jets in event
	Int_t njets=jets->GetEntries();
	AliTkConeJet *jet=0;
	for(Int_t i=0;i<njets;i++){
	  jet=(AliTkConeJet*)jets->At(i);
	  if(!jet) continue;
	  jet->calculateValues();
	  diff=AliTkConeJet::Diff(jt,jet,etdiff,phidiff,etadiff);
	  //_cut_ here if wanted (prob. for mixed events)
	  //if(phidiff>0.25||etadiff>0.25) continue;
	  //if(TMath::Abs(etdiff)/jt->getEt()>0.1) continue;
	  //cout << diff << " " << etdiff << " " << phidiff << " " << etadiff << endl;
	  if(mindiff>diff){
	    mindiff=diff;
	    index=i;
	  }
	}
	if(index>-1) {
	  jet=(AliTkConeJet*)jets->At(index);
	  hmJetEtres->Fill(jet->getEt());
	  hmJetEtratio->Fill(jet->getEt()/jt->getEt());
	  Float_t phidiff=(jt->getPhi()-jet->getPhi());
	  if(phidiff>TMath::Pi()) phidiff=TMath::TwoPi()-phidiff;
	  Float_t etadiff=(jt->getEta()-jet->getEta());
	  Float_t etfract=jet->getEtMarkedFrac();
	  hmAxesDiffres->Fill(phidiff,etadiff,TMath::Abs(jet->getEt()-jt->getEt()));
	  hmJetEtrestrue->Fill(jet->getEt(),etfract);
	  hmPhires->Fill(phidiff);
	  hmEtares->Fill(etadiff);
	}
      }
    }

    //could _cut_ on event
#if 0
    //cout << lead_tr_pt << " " << lead_tr_eta << " " << lead_tr_phi << endl;
    if((lead_tr_eta>0.5) || (lead_tr_eta<-0.5)) continue;
    if(zquench[0]<0.1||zquench[1]<0.1) continue;
#endif

    //particles in event
    const TClonesArray *parts=jetevent->GetParticles();

    //jets in event
    AliTkConeJet *lead_jet=0;
    AliTkConeJet *back_jet=0;
    for(Int_t i=0;i<njets;i++){
      AliTkConeJet* jet=0;
      Int_t whichjet=0;
      jet=(AliTkConeJet*)jets->At(i);
      if(!jet) continue;
      jet->calculateValues();

      hJetEtall->Fill(jet->getEt()); //without any cuts
      Float_t etfract=jet->getEtMarkedFrac();
      hJetEtalltrue->Fill(jet->getEt(),etfract); //without any cuts

      //here could _cut_ on jet
      //-------
      Float_t et=jet->getEt();
      Float_t corret=et;
      if(corrfac>0) corret/=corrfac;
      if(corret<minJetEt) continue;

      Int_t clindex=Int_t(corret/10);
      if(clindex>nclasses-2) clindex=nclasses-2;
      nclGoodEvents[clindex]++;
      nclGoodEvents[nclasses-1]++;
      //-------

      //set jet type
      for(Int_t i=0;i<nhard;i++){
	Float_t diff,etdiff,etadiff,phidiff;
	AliTkConeJet *jh=(AliTkConeJet*)chard_jets->At(i);
	diff=AliTkConeJet::Diff(jet,jh,etdiff,phidiff,etadiff);
	//cout << diff << " " << etdiff << " " << phidiff << " " << etadiff << endl;
	if(diff<0.25) {
	  jet->setType(jh->getType()); 
	  break;
	}
      }

      //check leading jet
      if(!lead_jet){ 
	lead_jet=jet;
	whichjet=1;
	hJetEt->Fill(lead_jet->getEt());
	hJetEttrue->Fill(lead_jet->getEt(),etfract);
	Float_t mindiff=100;
	Float_t minetdiff=0.,minphidiff=0.,minetadiff=0.;
	Int_t index;
	for(Int_t j=0;j<ntr;j++){
	  Float_t diff,etdiff,etadiff,phidiff;
	  AliTkConeJet *jt=(AliTkConeJet*)ctr_jets->At(j);
	  diff=AliTkConeJet::Diff(lead_jet,jt,etdiff,phidiff,etadiff);
	  if(mindiff>diff){
	    mindiff=diff;
	    index=j;
	    minphidiff=phidiff;
	    minetadiff=etadiff;
	    minetdiff=etdiff;
	  }
	}
	if(TMath::Abs(minetdiff)/lead_jet->getEt()<0.15)
	  hAxesDiff->Fill(minphidiff,minetadiff);

	//trigger
	Float_t triget=lead_jet->getEt();
	Float_t leadet=((AliTkConeJet*)ctr_jets->At(0))->getEt();
	for(Int_t i=0;i<9;i++){
	  Float_t minet=i*10+10;
	  if(triget>=minet) hJetEttrigger[i]->Fill(leadet,1);
	  else hJetEttrigger2[i]->Fill(leadet,1);
	}
	hJetEttriggernorm->Fill(leadet,1);
      }// check the back-to-back jet
      else if ((jet->getEt()/lead_jet->getEt()>0.75) &&
	        (relphi(jet->getPhi(),lead_jet->getPhi())>5/6*TMath::Pi())){
	if(!back_jet) {
	  whichjet=2;
	  back_jet=jet;
	  hBackJetEt->Fill(back_jet->getEt());
	  hBackJetEttrue->Fill(back_jet->getEt(),etfract);
	  hJetEtvsEt->Fill(lead_jet->getEt(),back_jet->getEt());

	} else{
	  cerr << "Already found one back jet, disregarding the new one." << endl;
	}
      }
      
      //fill properties
      hJetEtvsTrigger->Fill(corret,lead_tr_pt);
      hJetEtvsUQTrigger->Fill(corret,lead_uq_pt);
      hJetEtallvsEll->Fill(corret,prodlength);
      hJetEtvsEll->Fill(corret,prodlength);
      hJettype[jet->getType()]->Fill(jet->getEt());

      Int_t njetparts=jet->getNParts();
      hJetN[clindex]->Fill(njetparts);
      hJetN[nclasses-1]->Fill(njetparts);
      Float_t leadPartPt=jet->getPtLead();
      Float_t ptRatio=0;
      if(corret>0) ptRatio=leadPartPt/corret;
      Float_t meanpt=0.;
      if(njetparts>0) meanpt=corret/njetparts;
      hJetMeanPt[clindex]->Fill(meanpt);
      hJetMeanPt[nclasses-1]->Fill(meanpt);
      hJetLeadingPt[clindex]->Fill(ptRatio);
      hJetLeadingPt[nclasses-1]->Fill(ptRatio);
      hJetLeadingPtDist[clindex]->Fill(leadPartPt);
      hJetLeadingPtDist[nclasses-1]->Fill(leadPartPt);

      Float_t jetAxisX=jet->getXAxis();
      Float_t jetAxisY=jet->getYAxis();
      Float_t jetAxisZ=jet->getZAxis();
      Float_t jetAxisLength=jet->getPLength();
      if(jetAxisLength) {
	TClonesArray *particles = jet->getParts();
	TIterator *partit = particles->MakeIterator();
	TParticle *particle = NULL;
	Int_t firstval=0;
	while ((particle = (TParticle *) partit->Next()) != NULL) {
	  Float_t al=particle->Px()*jetAxisX+particle->Py()*jetAxisY+particle->Pz()*jetAxisZ;
	  if(al<0){
	    //cout << "Should not happen!" << al << endl;
	    continue;
	  }
	  Float_t at=TMath::Sqrt(TMath::Abs(particle->P()*particle->P()-al*al));
	  hJetFragL[clindex]->Fill(al/jetAxisLength);
	  hJetFragT[clindex]->Fill(at);
	  hJetFragPL[clindex]->Fill(al);
	  hJetFragPt[clindex]->Fill(particle->Pt());
	  hJetFragL[nclasses-1]->Fill(al/jetAxisLength);
	  hJetFragT[nclasses-1]->Fill(at);
	  hJetFragPL[nclasses-1]->Fill(al);
	  hJetFragPt[nclasses-1]->Fill(particle->Pt());
	  if(leadPartPt>0){
	    Float_t z=particle->Pt()/leadPartPt;
	    if(firstval||z<1){
	      hJetFragLeadingPt[clindex]->Fill(z);
	      hJetFragLeadingPt[nclasses-1]->Fill(z);
	    } else if(z==1) firstval=1;
	  }
	}
	delete partit;
      }

      //calculate int/diff shapes
      Float_t jeteta=jet->getCEta();
      Float_t jetphi=jet->getCPhi();
      AliTkEtaPhiVector center(jeteta,jetphi);
      TIterator *partit = parts->MakeIterator();
      AliJetParticle *particle = NULL;
      while ((particle = (AliJetParticle *) partit->Next()) != NULL) {
	Float_t pt=particle->Pt();
	if(pt<ptcutused) continue;
	AliTkEtaPhiVector centerp;
	centerp.setVector(particle->Eta(),particle->Phi());
	Float_t diff=center.diff(centerp);
	for(Int_t loop=Int_t(diff*10)+1;loop<=10;loop++){
	  if(pt<1.) retlow1[clindex][loop-1]+=pt;
	  if(pt<2.) retlow[clindex][loop-1]+=pt;
	  if(pt<3.) retlow2[clindex][loop-1]+=pt;
	  if(pt<4.) retlow3[clindex][loop-1]+=pt;
	  if(pt<5.) retlow4[clindex][loop-1]+=pt;
	  retall[clindex][loop-1]+=pt;
	  if(pt<1.) retlow1[nclasses-1][loop-1]+=pt;
	  if(pt<2.) retlow[nclasses-1][loop-1]+=pt;
	  if(pt<3.) retlow2[nclasses-1][loop-1]+=pt;
	  if(pt<4.) retlow3[nclasses-1][loop-1]+=pt;
	  if(pt<5.) retlow4[nclasses-1][loop-1]+=pt;
	  retall[nclasses-1][loop-1]+=pt;
	}
	if(diff<=1){
	  Int_t index=0;
	  if(diff<=deltaR) index=0;
	  else index=Int_t((diff-deltaR)*10)+1;
	  if(pt<1.) dretlow1[clindex][index]+=pt;
	  if(pt<2.) dretlow[clindex][index]+=pt;
	  if(pt<3.) dretlow2[clindex][index]+=pt;
	  if(pt<4.) dretlow3[clindex][index]+=pt;
	  if(pt<5.) dretlow4[clindex][index]+=pt;
	  dretall[clindex][index]+=pt;
	  if(pt<1.) dretlow1[nclasses-1][index]+=pt;
	  if(pt<2.) dretlow[nclasses-1][index]+=pt;
	  if(pt<3.) dretlow2[nclasses-1][index]+=pt;
	  if(pt<4.) dretlow3[nclasses-1][index]+=pt;
	  if(pt<5.) dretlow4[nclasses-1][index]+=pt;
	  dretall[nclasses-1][index]+=pt;
	}
	if(pt<ptcutused) continue;
	Float_t phi=particle->Phi();

	Float_t dphi=diffphi(jet->getPhi(),phi);
	hPhiCorr[clindex]->Fill(dphi);
	hPhiCorr[nclasses-1]->Fill(dphi);
      }
      delete partit;

    } //loop over cone jets

    //global event studies
    TIterator *partit = parts->MakeIterator();
    AliJetParticle *particle = NULL;
    Float_t leta=0.;
    Float_t lphi=0.;
    Int_t clindex=-1;
    Float_t jetAxisX=-1;
    Float_t jetAxisY=-1;
    Float_t jetAxisZ=-1;
    Float_t jetAxisLength=-1;
    Float_t leadPartPt=-1;
    if(lead_jet){
      leta=lead_jet->getCEta();
      lphi=lead_jet->getCPhi();
      jetAxisX=lead_jet->getXAxis();
      jetAxisY=lead_jet->getYAxis();
      jetAxisZ=lead_jet->getZAxis();
      jetAxisLength=lead_jet->getPLength();
      leadPartPt=lead_jet->getPtLead();
      centers[0]=new AliTkEtaPhiVector(leta,lphi);
      centers[1]=new AliTkEtaPhiVector(leta,addphi(lphi,TMath::Pi()));
      centers[2]=new AliTkEtaPhiVector(leta,addphi(lphi,TMath::Pi()/2));
      centers[3]=new AliTkEtaPhiVector(leta,addphi(lphi,3*TMath::Pi()/2));
      Float_t et=lead_jet->getEt();
      Float_t corret=et;
      if(corrfac>0) corret/=corrfac;
      clindex=Int_t(corret/10);
      if(clindex>nclasses-2) clindex=nclasses-2;
      nclLeadEvents[clindex]++;
      nclLeadEvents[nclasses-1]++;
      if(back_jet){
	nclDiEvents[clindex]++;
	nclDiEvents[nclasses-1]++;
      }
    }

    //loop over particles in event
    Int_t firstval=0;
    Int_t difirstval=0;
    while ((particle = (AliJetParticle *) partit->Next()) != NULL) {
      Float_t pt=particle->Pt();
      if(pt<ptcutused) continue;
      Float_t eta=particle->Eta();
      Float_t phi=particle->Phi();
      hPartPtDist->Fill(pt);
      hPartPhiDist->Fill(phi);
      hPartEtaDist->Fill(eta);
      if(lead_jet){
	Float_t dphi=diffphi(lphi,phi);
	Float_t al=particle->Px()*jetAxisX+particle->Py()*jetAxisY+particle->Pz()*jetAxisZ;
	al/=particle->P();
	hPartPhiCorr->Fill(dphi);
	hPartACorr->Fill(al);
	if(back_jet){
	  hPartDiPhiCorr->Fill(dphi);
	  hPartDiACorr->Fill(al);
	}

	//ue studies
	for(Int_t ceindex=0;ceindex<4;ceindex++){
	  Int_t heindex=ceindex;
	  if(heindex==3) heindex=2; //store both sides of trans plane in one histo

	  AliTkEtaPhiVector centerp;
	  centerp.setVector(eta,phi);
	  Float_t diff=centers[ceindex]->diff(centerp);
	  for(Int_t loop=Int_t(diff*10)+1;loop<=10;loop++){
	    if(pt<1.) gretlow1[heindex][clindex][loop-1]+=pt;
	    if(pt<2.) gretlow[heindex][clindex][loop-1]+=pt;
	    if(pt<3.) gretlow2[heindex][clindex][loop-1]+=pt;
	    if(pt<4.) gretlow3[heindex][clindex][loop-1]+=pt;
	    if(pt<5.) gretlow4[heindex][clindex][loop-1]+=pt;
	    gretall[heindex][clindex][loop-1]+=pt;
	    if(pt<1.) gretlow1[heindex][nclasses-1][loop-1]+=pt;
	    if(pt<2.) gretlow[heindex][nclasses-1][loop-1]+=pt;
	    if(pt<3.) gretlow2[heindex][nclasses-1][loop-1]+=pt;
	    if(pt<4.) gretlow3[heindex][nclasses-1][loop-1]+=pt;
	    if(pt<5.) gretlow4[heindex][nclasses-1][loop-1]+=pt;
	    gretall[heindex][nclasses-1][loop-1]+=pt;
	  }
	  if(diff<=1){
	    Int_t index=0;
	    if(diff<=deltaR) index=0;
	    else index=Int_t((diff-deltaR)*10)+1;
	    if(pt<1.) gdretlow1[heindex][clindex][index]+=pt;
	    if(pt<2.) gdretlow[heindex][clindex][index]+=pt;
	    if(pt<3.) gdretlow2[heindex][clindex][index]+=pt;
	    if(pt<4.) gdretlow3[heindex][clindex][index]+=pt;
	    if(pt<5.) gdretlow4[heindex][clindex][index]+=pt;
	    gdretall[heindex][clindex][index]+=pt;
	    if(pt<1.) gdretlow1[heindex][nclasses-1][index]+=pt;
	    if(pt<2.) gdretlow[heindex][nclasses-1][index]+=pt;
	    if(pt<3.) gdretlow2[heindex][nclasses-1][index]+=pt;
	    if(pt<4.) gdretlow3[heindex][nclasses-1][index]+=pt;
	    if(pt<5.) gdretlow4[heindex][nclasses-1][index]+=pt;
	    gdretall[heindex][nclasses-1][index]+=pt;
	  }
	  if(back_jet){
	    for(Int_t loop=Int_t(diff*10)+1;loop<=10;loop++){
	      if(pt<1.) gdiretlow1[heindex][clindex][loop-1]+=pt;
	      if(pt<2.) gdiretlow[heindex][clindex][loop-1]+=pt;
	      if(pt<3.) gdiretlow2[heindex][clindex][loop-1]+=pt;
	      if(pt<4.) gdiretlow3[heindex][clindex][loop-1]+=pt;
	      if(pt<5.) gdiretlow4[heindex][clindex][loop-1]+=pt;
	      gdiretall[heindex][clindex][loop-1]+=pt;
	      if(pt<1.) gdiretlow1[heindex][nclasses-1][loop-1]+=pt;
	      if(pt<2.) gdiretlow[heindex][nclasses-1][loop-1]+=pt;
	      if(pt<3.) gdiretlow2[heindex][nclasses-1][loop-1]+=pt;
	      if(pt<4.) gdiretlow3[heindex][nclasses-1][loop-1]+=pt;
	      if(pt<5.) gdiretlow4[heindex][nclasses-1][loop-1]+=pt;
	      gdiretall[heindex][nclasses-1][loop-1]+=pt;
	    }
	    if(diff<=1){
	      Int_t index=0;
	      if(diff<=deltaR) index=0;
	      else index=Int_t((diff-deltaR)*10)+1;
	      if(pt<1.) gdidretlow1[heindex][clindex][index]+=pt;
	      if(pt<2.) gdidretlow[heindex][clindex][index]+=pt;
	      if(pt<3.) gdidretlow2[heindex][clindex][index]+=pt;
	      if(pt<4.) gdidretlow3[heindex][clindex][index]+=pt;
	      if(pt<5.) gdidretlow4[heindex][clindex][index]+=pt;
	      gdidretall[heindex][clindex][index]+=pt;
	      if(pt<1.) gdretlow1[heindex][nclasses-1][index]+=pt;
	      if(pt<2.) gdretlow[heindex][nclasses-1][index]+=pt;
	      if(pt<3.) gdidretlow2[heindex][nclasses-1][index]+=pt;
	      if(pt<4.) gdidretlow3[heindex][nclasses-1][index]+=pt;
	      if(pt<5.) gdidretlow4[heindex][nclasses-1][index]+=pt;
	      gdidretall[heindex][nclasses-1][index]+=pt;
	    }
	  }
	}//ceindex

	Int_t ceindex=eventindex(lphi,phi);
	Int_t heindex=ceindex;
	if(heindex==3) heindex=2; //store both sides of trans plane in one histo

	if(jetAxisLength) {
	  Float_t jetX=jetAxisX;
	  Float_t jetY=jetAxisY;
	  if(ceindex==1) {
	    jetX=-jetAxisX;
	    jetY=-jetAxisY;
	  } else if(ceindex==2) {
	    jetX=-jetAxisY;
	    jetY=jetAxisX;
	  }
	  else if(ceindex==3) {
	    jetX=jetAxisY;
	    jetY=-jetAxisX;
	  }
	  Float_t al=particle->Px()*jetX+particle->Py()*jetY+particle->Pz()*jetAxisZ;
	  Float_t at=TMath::Sqrt(TMath::Abs(particle->P()*particle->P()-al*al));
	  hgJetFragL[heindex][clindex]->Fill(al/jetAxisLength);
	  hgJetFragT[heindex][clindex]->Fill(at);
	  hgJetFragPL[heindex][clindex]->Fill(al);
	  hgJetFragPt[heindex][clindex]->Fill(particle->Pt());
	  hgJetFragL[heindex][nclasses-1]->Fill(al/jetAxisLength);
	  hgJetFragT[heindex][nclasses-1]->Fill(at);
	  hgJetFragPL[heindex][nclasses-1]->Fill(al);
	  hgJetFragPt[heindex][nclasses-1]->Fill(particle->Pt());
	  if(leadPartPt>0){
	    Float_t z=particle->Pt()/leadPartPt;
	    if(firstval||z<1){
	      hgJetFragLeadingPt[heindex][clindex]->Fill(z);
	      hgJetFragLeadingPt[heindex][nclasses-1]->Fill(z);
	    } else if(z==1) firstval=1;
	  }
	  if(back_jet){
	    hgDiJetFragL[heindex][clindex]->Fill(al/jetAxisLength);
	    hgDiJetFragT[heindex][clindex]->Fill(at);
	    hgDiJetFragPL[heindex][clindex]->Fill(al);
	    hgDiJetFragPt[heindex][clindex]->Fill(particle->Pt());
	    hgDiJetFragL[heindex][nclasses-1]->Fill(al/jetAxisLength);
	    hgDiJetFragT[heindex][nclasses-1]->Fill(at);
	    hgDiJetFragPL[heindex][nclasses-1]->Fill(al);
	    hgDiJetFragPt[heindex][nclasses-1]->Fill(particle->Pt());
	    if(leadPartPt>0){
	      Float_t z=particle->Pt()/leadPartPt;
	      if(difirstval||z<1){
		hgDiJetFragLeadingPt[heindex][clindex]->Fill(z);
		hgDiJetFragLeadingPt[heindex][nclasses-1]->Fill(z);
	      } else if(z==1) difirstval=1;
	    }
	  }
	}
      } //lead_jet

    }
    delete partit;
    if(lead_jet) 
      for(Int_t i=0;i<4;i++) delete centers[i];

    nTotalEvents++;
  } //end of nev loop

  cout << "Finished analysing events " << nTotalEvents << endl;

  delete centers;
  delete event;
  delete theTree;
  if(theEvTree){
    delete jetevent;
    delete theEvTree;
  }
  if(theSignalTree){
    delete jetsigevent;
    delete theSignalTree;
  }
  if(theTreeMonte){
    delete monteevent;
    delete theTreeMonte;
  }

  //========================================================================
  // draw histograms
  //========================================================================

  TTimeStamp tstamp;
  Char_t timestamp[255];
  sprintf(timestamp,"%d",tstamp.GetDate(0,0));  

  //store all objects in root file (you never know, when you need it)
  Char_t rootfilename[1024];
  if(savefilename!=NULL)
    sprintf(rootfilename,"%s",savefilename);
  else 
    sprintf(rootfilename,"%s/anaAliJets-%s.root",gSystem->Getenv("JF_PLOTSDIR"),timestamp);
  TFile *rootfile=new TFile(rootfilename,"RECREATE");

  // let's start with the drawing...
  Int_t nents=(Int_t)hJetEt->GetEntries();
  if(nents){
    cout << "hJetEt: " << nents << " entries in histogram" << endl;
    hJetEt->SetLineWidth(3);
    hJetEt->Write();
  }

  nents=(Int_t)hJetEttrue->GetEntries();
  if(nents){
    cout << "hJetEttrue: " << nents << " entries in histogram" << endl;
    hJetEttrue->SetLineWidth(3);
    hJetEttrue->Write();
  }

  nents=(Int_t)hBackJetEt->GetEntries();
  if(nents){
    cout << "hBackJetEt: " << nents << " entries in histogram" << endl;
    hBackJetEt->SetLineWidth(3);
    hBackJetEt->Write();
  }

  nents=(Int_t)hBackJetEttrue->GetEntries();
  if(nents){
    cout << "hBackJetEttrue: " << nents << " entries in histogram" << endl;
    hBackJetEttrue->SetLineWidth(3);
    hBackJetEttrue->Write();
  }

  nents=(Int_t)hJetEtall->GetEntries();
  if(nents){
    cout << "hJetEtall: " << nents << " entries in histogram" << endl;
    hJetEtall->SetLineWidth(3);
    hJetEtall->Write();
  }

  nents=(Int_t)hJetEtalltrue->GetEntries();
  if(nents){
    cout << "hJetEtalltrue: " << nents << " entries in histogram" << endl;
    hJetEtalltrue->SetLineWidth(3);
    hJetEtalltrue->Write();
  }

  nents=(Int_t)hJetEtUQTrigger->GetEntries();
  if(nents){
    cout << "hJetEtUQTrigger: " << nents << " entries in histogram" << endl;
    hJetEtUQTrigger->SetLineWidth(3);
    hJetEtUQTrigger->Write();
  }
  nents=(Int_t)hJetEtTrigger->GetEntries();
  if(nents){
    cout << "hJetEtTrigger: " << nents << " entries in histogram" << endl;
    hJetEtTrigger->SetLineWidth(3);
    hJetEtTrigger->Write();
  }

  rootfile->cd();
  rootfile->mkdir("LeadingPt");
  rootfile->cd("LeadingPt");
  for(Int_t i=0;i<nclasses;i++){
    nents=(Int_t)hJetLeadingPt[i]->GetEntries();
    if(nents){
      //cout << "hJetLeadingPt " << i << ": "  << nents << " entries in histogram" << endl;
      hJetLeadingPt[i]->SetLineWidth(3);
      hJetLeadingPt[i]->Write();
    }
  }

  rootfile->cd();
  rootfile->mkdir("FragLeadingPt");
  rootfile->cd("FragLeadingPt");
  for(Int_t i=0;i<nclasses;i++){
    nents=(Int_t)hJetFragLeadingPt[i]->GetEntries();
    if(nents){
      //cout << "hJetFragLeadingPt " << i << ": "  << nents << " entries in histogram" << endl;
      hJetFragLeadingPt[i]->SetLineWidth(3);
      hJetFragLeadingPt[i]->Write();
    }
  }

  rootfile->cd();
  rootfile->mkdir("LeadingPtDist");
  rootfile->cd("LeadingPtDist");
  for(Int_t i=0;i<nclasses;i++){
    nents=(Int_t)hJetLeadingPtDist[i]->GetEntries();
    if(nents){
      //cout << "hJetLeadingPtDist " << i << ": " << nents << " entries in histogram" << endl;
      hJetLeadingPtDist[i]->SetLineWidth(3);
      hJetLeadingPtDist[i]->Write();
    }
  }

  rootfile->cd();
  rootfile->mkdir("FragLong");
  rootfile->cd("FragLong");
  for(Int_t i=0;i<nclasses;i++){
    nents=(Int_t)hJetFragL[i]->GetEntries();
    if(nents){
      //cout << "hJetFragL " << i << ": " << nents << " entries in histogram " << endl;
      hJetFragL[i]->SetLineWidth(3);
      hJetFragL[i]->Write();
    }
  }

  rootfile->cd();
  rootfile->mkdir("FragPL");
  rootfile->cd("FragPL");
  for(Int_t i=0;i<nclasses;i++){
    nents=(Int_t)hJetFragPL[i]->GetEntries();
    if(nents){
      //cout << "hJetFragPL " << i << ": " << nents << " entries in histogram " << endl;
      hJetFragPL[i]->SetLineWidth(3);
      hJetFragPL[i]->Write();
    }
  }

  rootfile->cd();
  rootfile->mkdir("FragTrans");
  rootfile->cd("FragTrans");
  for(Int_t i=0;i<nclasses;i++){
    nents=(Int_t)hJetFragT[i]->GetEntries();
    if(nents){
      //cout << "hJetFragT " << i << ": " << nents << " entries in histogram" << endl;
      hJetFragT[i]->SetLineWidth(3);
      hJetFragT[i]->Write();
    }
  }

  rootfile->cd();
  rootfile->mkdir("Multiplicity");
  rootfile->cd("Multiplicity");
  for(Int_t i=0;i<nclasses;i++){
    nents=(Int_t)hJetN[i]->GetEntries();
    //cout <<"hJetN " << i << ": " << nents << " entries in histogram " << endl;  
    if(nents){
      hJetN[i]->SetLineWidth(3);
      hJetN[i]->Write();
    }
  }

  rootfile->cd();
  rootfile->mkdir("PtDistribution");
  rootfile->cd("PtDistribution");
  for(Int_t i=0;i<nclasses;i++){
    nents=(Int_t)hJetFragPt[i]->GetEntries();
    if(nents){
      //cout << "hJetFragPt " << i << ": " << nents << " entries in histogram " << endl;  
      hJetFragPt[i]->SetLineWidth(3);
      hJetFragPt[i]->Write();
    }
  }

  rootfile->cd();
  rootfile->mkdir("MeanPt");
  rootfile->cd("MeanPt");
  for(Int_t i=0;i<nclasses;i++){
    nents=(Int_t)hJetMeanPt[i]->GetEntries();
    if(nents){
      //cout << "hJetMeanPt " << i << ": " << nents << " entries in histogram " << endl;  
      hJetMeanPt[i]->SetLineWidth(3);
      hJetMeanPt[i]->Write();
    }
  }

  rootfile->cd();
  rootfile->mkdir("PhiCorr");
  rootfile->cd("PhiCorr");
  for(Int_t i=0;i<nclasses;i++){
    nents=(Int_t)hPhiCorr[i]->GetEntries();
    if(nents){
      //cout << "hPhiCorr " << i << ": " << nents << " entries in histogram " << endl;  
      hPhiCorr[i]->SetLineWidth(3);
      hPhiCorr[i]->Write();
    }
  }

  rootfile->cd();//intshape
  rootfile->mkdir("ConeFluc");
  rootfile->cd("ConeFluc");
  for(Int_t j=0;j<nclasses;j++){
    TGraph *graphall=new TGraph(10);
    TGraph *graphlow=new TGraph(10);
    TGraph *graphlow1=new TGraph(10);
    TGraph *graphlow2=new TGraph(10);
    TGraph *graphlow3=new TGraph(10);
    TGraph *graphlow4=new TGraph(10);
    for(Int_t i=0;i<10;i++) {
      graphall->SetPoint(i,rxet[i],retall[j][i]);
      graphlow->SetPoint(i,rxet[i],retlow[j][i]);
      graphlow1->SetPoint(i,rxet[i],retlow1[j][i]);
      graphlow2->SetPoint(i,rxet[i],retlow2[j][i]);
      graphlow3->SetPoint(i,rxet[i],retlow3[j][i]);
      graphlow4->SetPoint(i,rxet[i],retlow4[j][i]);
    }
    sprintf(dummy,"gretall%d",j);
    graphall->Write(dummy);
    sprintf(dummy,"gretlow%d",j);
    graphlow->Write(dummy);
    sprintf(dummy,"gret2low%d",j);
    graphlow->Write(dummy);
    sprintf(dummy,"gret1low%d",j);
    graphlow1->Write(dummy);
    sprintf(dummy,"gret3low%d",j);
    graphlow2->Write(dummy);
    sprintf(dummy,"gret4low%d",j);
    graphlow3->Write(dummy);
    sprintf(dummy,"gret5low%d",j);
    graphlow4->Write(dummy);
    delete graphall;
    delete graphlow;
    delete graphlow1;
    delete graphlow2;
    delete graphlow3;
  }

  rootfile->cd();
  rootfile->mkdir("Shape");
  rootfile->cd("Shape");
  for(Int_t j=0;j<nclasses;j++){
    TGraph *graphall=new TGraph(11);
    TGraph *graphlow=new TGraph(11);
    TGraph *graphlow1=new TGraph(11);
    TGraph *graphlow2=new TGraph(11);
    TGraph *graphlow3=new TGraph(11);
    TGraph *graphlow4=new TGraph(11);
    for(Int_t i=0;i<11;i++) {
      graphall->SetPoint(i,drxet[i],dretall[j][i]/deltaR);
      graphlow->SetPoint(i,drxet[i],dretlow[j][i]/deltaR);
      graphlow1->SetPoint(i,drxet[i],dretlow1[j][i]/deltaR);
      graphlow2->SetPoint(i,drxet[i],dretlow2[j][i]/deltaR);
      graphlow3->SetPoint(i,drxet[i],dretlow3[j][i]/deltaR);
      graphlow4->SetPoint(i,drxet[i],dretlow4[j][i]/deltaR);
    }
    sprintf(dummy,"gretall%d",j);
    graphall->Write(dummy);
    sprintf(dummy,"gretlow%d",j);
    graphlow->Write(dummy);
    sprintf(dummy,"gret2low%d",j);
    graphlow->Write(dummy);
    sprintf(dummy,"gret1low%d",j);
    graphlow1->Write(dummy);
    sprintf(dummy,"gret3low%d",j);
    graphlow2->Write(dummy);
    sprintf(dummy,"gret4low%d",j);
    graphlow3->Write(dummy);
    sprintf(dummy,"gret5low%d",j);
    graphlow4->Write(dummy);
    delete graphall;
    delete graphlow;
    delete graphlow1;
    delete graphlow2;
    delete graphlow3;
  }

  rootfile->cd();
  rootfile->mkdir("Extended");
  rootfile->cd("Extended");

  nents=(Int_t)hJetEtvsTrigger->GetEntries();
  if(nents){
    cout << "hJetEtvsTrigger" << nents << " entries in histogram" << endl;
    hJetEtvsTrigger->SetLineWidth(3);
    hJetEtvsTrigger->Write();
  }

  nents=(Int_t)hJetEtvsEt->GetEntries();
  if(nents){
    cout << "hJetEtvsEt" << nents << " entries in histogram" << endl;
    hJetEtvsEt->SetLineWidth(3);
    hJetEtvsEt->Write();
  }

  nents=(Int_t)hAxesDiff->GetEntries();
  if(nents){
    cout << "hAxesDiff" << nents << " entries in histogram" << endl;
    hAxesDiff->SetLineWidth(3);
    hAxesDiff->Write();
  }

  nents=(Int_t)hJet1->GetEntries();
  if(nents){
    cout << "hJet1" << nents << " entries in histogram" << endl;
    hJet1->SetLineWidth(3);
    hJet1->Write();
  }
  nents=(Int_t)hJet2->GetEntries();
  if(nents){
    cout << "hJet2" << nents << " entries in histogram" << endl;
    hJet2->SetLineWidth(3);
    hJet2->Write();
  }

  for(Int_t i=0;i<3;i++){
    nents=(Int_t)hJettype[i]->GetEntries();
    if(nents){
      hJettype[i]->SetLineWidth(3);
      hJettype[i]->Write();
    }
  }

  nents=(Int_t)hJetEtvsEll->GetEntries();
  if(nents){
    cout << "hJetEtvsEll" << nents << " entries in histogram " << endl;  
    hJetEtvsEll->SetLineWidth(3);
    hJetEtvsEll->Write();
  }

  nents=(Int_t)hJetEtallvsEll->GetEntries();
  if(nents){
    cout << "hJetEtallvsEll" << nents << " entries in histogram " << endl;  
    hJetEtallvsEll->SetLineWidth(3);
    hJetEtallvsEll->Write();
  }

  nents=(Int_t)hJetZ->GetEntries();
  if(nents){
    cout << "hJetZ: " << nents << " entries in histogram " << endl;  
    hJetZ->SetLineWidth(3);
    hJetZ->Write();
  }

  rootfile->cd();
  rootfile->mkdir("Global");
  rootfile->cd("Global");

  nents=(Int_t)hPartPtDist->GetEntries();
  if(nents){
    cout << "hPartPtDist: " << nents << " entries in histogram" << endl;
    hPartPtDist->SetLineWidth(3);
    hPartPtDist->Write();
  }

  nents=(Int_t)hPartEtaDist->GetEntries();
  if(nents){
    cout << "hPartEtaDist: " << nents << " entries in histogram" << endl;
    hPartEtaDist->SetLineWidth(3);
    hPartEtaDist->Write();
  }

  nents=(Int_t)hPartPhiDist->GetEntries();
  if(nents){
    cout << "hPartPhiDist: " << nents << " entries in histogram" << endl;
    hPartPhiDist->SetLineWidth(3);
    hPartPhiDist->Write();
  }

  nents=(Int_t)hPartPhiCorr->GetEntries();
  if(nents){
    cout << "hPartPhiCorr: " << nents << " entries in histogram" << endl;
    hPartPhiCorr->SetLineWidth(3);
    hPartPhiCorr->Write();
  }

  nents=(Int_t)hPartDiPhiCorr->GetEntries();
  if(nents){
    cout << "hPartDiPhiCorr: " << nents << " entries in histogram" << endl;
    hPartDiPhiCorr->SetLineWidth(3);
    hPartDiPhiCorr->Write();
  }

  nents=(Int_t)hPartACorr->GetEntries();
  if(nents){
    cout << "hPartACorr: " << nents << " entries in histogram" << endl;
    hPartACorr->SetLineWidth(3);
    hPartACorr->Write();
  }

  nents=(Int_t)hPartDiACorr->GetEntries();
  if(nents){
    cout << "hPartDiACorr: " << nents << " entries in histogram" << endl;
    hPartDiACorr->SetLineWidth(3);
    hPartDiACorr->Write();
  }

  rootfile->cd();
  rootfile->mkdir("gConeFluc");
  rootfile->cd("gConeFluc");
  for(Int_t k=0;k<3;k++){
    if(k==0){
      sprintf(name,"toward");
    } else if (k==1) {
      sprintf(name,"away");
    } else {
      sprintf(name,"transverse");
    }

    for(Int_t j=0;j<nclasses;j++){
      TGraph *graphall=new TGraph(10);
      TGraph *graphlow=new TGraph(10);
      TGraph *graphlow1=new TGraph(10);
      TGraph *graphlow2=new TGraph(10);
      TGraph *graphlow3=new TGraph(10);
      TGraph *graphlow4=new TGraph(10);
      for(Int_t i=0;i<10;i++) {
	graphall->SetPoint(i,grxet[i],gretall[k][j][i]);
	graphlow->SetPoint(i,grxet[i],gretlow[k][j][i]);
	graphlow1->SetPoint(i,grxet[i],gretlow1[k][j][i]);
	graphlow2->SetPoint(i,grxet[i],gretlow2[k][j][i]);
	graphlow3->SetPoint(i,grxet[i],gretlow3[k][j][i]);
	graphlow4->SetPoint(i,grxet[i],gretlow4[k][j][i]);
      }
      sprintf(dummy,"%s-gretall%d",name,j);
      graphall->Write(dummy);
      sprintf(dummy,"%s-gretlow%d",name,j);
      graphlow->Write(dummy);
      sprintf(dummy,"%s-gret2low%d",name,j);
      graphlow->Write(dummy);
      sprintf(dummy,"%s-gret1low%d",name,j);
      graphlow1->Write(dummy);
      sprintf(dummy,"%s-gret3low%d",name,j);
      graphlow2->Write(dummy);
      sprintf(dummy,"%s-gret4low%d",name,j);
      graphlow3->Write(dummy);
      sprintf(dummy,"%s-gret5low%d",name,j);
      graphlow4->Write(dummy);
      delete graphall;
      delete graphlow;
      delete graphlow1;
      delete graphlow2;
      delete graphlow3;
    }
  }

  rootfile->cd();
  rootfile->mkdir("gShape");
  rootfile->cd("gShape");
  for(Int_t k=0;k<3;k++){
    if(k==0){
      sprintf(name,"toward");
    } else if (k==1) {
      sprintf(name,"away");
    } else {
      sprintf(name,"transverse");
    }

    for(Int_t j=0;j<nclasses;j++){
      TGraph *graphall=new TGraph(10);
      TGraph *graphlow=new TGraph(10);
      TGraph *graphlow1=new TGraph(10);
      TGraph *graphlow2=new TGraph(10);
      TGraph *graphlow3=new TGraph(10);
      TGraph *graphlow4=new TGraph(10);
      for(Int_t i=0;i<10;i++) {
	graphall->SetPoint(i,gdrxet[i],gdretall[k][j][i]/deltaR);
	graphlow->SetPoint(i,gdrxet[i],gdretlow[k][j][i]/deltaR);
	graphlow1->SetPoint(i,gdrxet[i],gdretlow1[k][j][i]/deltaR);
	graphlow2->SetPoint(i,gdrxet[i],gdretlow2[k][j][i]/deltaR);
	graphlow3->SetPoint(i,gdrxet[i],gdretlow3[k][j][i]/deltaR);
	graphlow4->SetPoint(i,gdrxet[i],gdretlow4[k][j][i]/deltaR);
      }
      sprintf(dummy,"%s-gretall%d",name,j);
      graphall->Write(dummy);
      sprintf(dummy,"%s-gretlow%d",name,j);
      graphlow->Write(dummy);
      sprintf(dummy,"%s-gret2low%d",name,j);
      graphlow->Write(dummy);
      sprintf(dummy,"%s-gret1low%d",name,j);
      graphlow1->Write(dummy);
      sprintf(dummy,"%s-gret3low%d",name,j);
      graphlow2->Write(dummy);
      sprintf(dummy,"%s-gret4low%d",name,j);
      graphlow3->Write(dummy);
      sprintf(dummy,"%s-gret5low%d",name,j);
      graphlow4->Write(dummy);
      delete graphall;
      delete graphlow;
      delete graphlow1;
      delete graphlow2;
      delete graphlow3;
    }
  }

  rootfile->cd();
  rootfile->mkdir("gDiConeFluc");
  rootfile->cd("gDiConeFluc");
  for(Int_t k=0;k<3;k++){
    if(k==0){
      sprintf(name,"toward");
    } else if (k==1) {
      sprintf(name,"away");
    } else {
      sprintf(name,"transverse");
    }

    for(Int_t j=0;j<nclasses;j++){
      TGraph *graphall=new TGraph(10);
      TGraph *graphlow=new TGraph(10);
      TGraph *graphlow1=new TGraph(10);
      TGraph *graphlow2=new TGraph(10);
      TGraph *graphlow3=new TGraph(10);
      TGraph *graphlow4=new TGraph(10);
      for(Int_t i=0;i<10;i++) {
	graphall->SetPoint(i,grxet[i],gdiretall[k][j][i]);
	graphlow->SetPoint(i,grxet[i],gdiretlow[k][j][i]);
	graphlow1->SetPoint(i,grxet[i],gdiretlow1[k][j][i]);
	graphlow2->SetPoint(i,grxet[i],gdiretlow2[k][j][i]);
	graphlow3->SetPoint(i,grxet[i],gdiretlow3[k][j][i]);
	graphlow4->SetPoint(i,grxet[i],gdiretlow4[k][j][i]);
      }
      sprintf(dummy,"%s-gretall%d",name,j);
      graphall->Write(dummy);
      sprintf(dummy,"%s-gretlow%d",name,j);
      graphlow->Write(dummy);
      sprintf(dummy,"%s-gret2low%d",name,j);
      graphlow->Write(dummy);
      sprintf(dummy,"%s-gret1low%d",name,j);
      graphlow1->Write(dummy);
      sprintf(dummy,"%s-gret3low%d",name,j);
      graphlow2->Write(dummy);
      sprintf(dummy,"%s-gret4low%d",name,j);
      graphlow3->Write(dummy);
      sprintf(dummy,"%s-gret5low%d",name,j);
      graphlow4->Write(dummy);
      delete graphall;
      delete graphlow;
      delete graphlow1;
      delete graphlow2;
      delete graphlow3;
    }
  }

  rootfile->cd();
  rootfile->mkdir("gDiShape");
  rootfile->cd("gDiShape");
  for(Int_t k=0;k<3;k++){
    if(k==0){
      sprintf(name,"toward");
    } else if (k==1) {
      sprintf(name,"away");
    } else {
      sprintf(name,"transverse");
    }

    for(Int_t j=0;j<nclasses;j++){
      TGraph *graphall=new TGraph(10);
      TGraph *graphlow=new TGraph(10);
      TGraph *graphlow1=new TGraph(10);
      TGraph *graphlow2=new TGraph(10);
      TGraph *graphlow3=new TGraph(10);
      TGraph *graphlow4=new TGraph(10);
      for(Int_t i=0;i<10;i++) {
	graphall->SetPoint(i,gdrxet[i],gdidretall[k][j][i]/deltaR);
	graphlow->SetPoint(i,gdrxet[i],gdidretlow[k][j][i]/deltaR);
	graphlow1->SetPoint(i,gdrxet[i],gdidretlow1[k][j][i]/deltaR);
	graphlow2->SetPoint(i,gdrxet[i],gdidretlow2[k][j][i]/deltaR);
	graphlow3->SetPoint(i,gdrxet[i],gdidretlow3[k][j][i]/deltaR);
	graphlow4->SetPoint(i,gdrxet[i],gdidretlow4[k][j][i]/deltaR);
      }
      sprintf(dummy,"%s-gretall%d",name,j);
      graphall->Write(dummy);
      sprintf(dummy,"%s-gretlow%d",name,j);
      graphlow->Write(dummy);
      sprintf(dummy,"%s-gret2low%d",name,j);
      graphlow->Write(dummy);
      sprintf(dummy,"%s-gret1low%d",name,j);
      graphlow1->Write(dummy);
      sprintf(dummy,"%s-gret3low%d",name,j);
      graphlow2->Write(dummy);
      sprintf(dummy,"%s-gret4low%d",name,j);
      graphlow3->Write(dummy);
      sprintf(dummy,"%s-gret5low%d",name,j);
      graphlow4->Write(dummy);
      delete graphall;
      delete graphlow;
      delete graphlow1;
      delete graphlow2;
      delete graphlow3;
    }
  }

  rootfile->cd();
  rootfile->mkdir("gFragLeadingPt");
  rootfile->cd("gFragLeadingPt");
  for(Int_t k=0;k<3;k++){
    for(Int_t i=0;i<nclasses;i++){
      nents=(Int_t)hgJetFragLeadingPt[k][i]->GetEntries();
      if(nents){
	//cout << "hJetFragLeadingPt " << i << ": "  << nents << " entries in histogram" << endl;
	hgJetFragLeadingPt[k][i]->SetLineWidth(3);
	hgJetFragLeadingPt[k][i]->Write();
      }
    }
  }

  rootfile->cd();
  rootfile->mkdir("gFragLong");
  rootfile->cd("gFragLong");
  for(Int_t k=0;k<3;k++){
    for(Int_t i=0;i<nclasses;i++){
      nents=(Int_t)hgJetFragL[k][i]->GetEntries();
      if(nents){
	//cout << "hgJetFragL " << i << ": " << nents << " entries in histogram " << endl;
	hgJetFragL[k][i]->SetLineWidth(3);
	hgJetFragL[k][i]->Write();
      }
    }
  }

  rootfile->cd();
  rootfile->mkdir("gFragPL");
  rootfile->cd("gFragPL");
  for(Int_t k=0;k<3;k++){
    for(Int_t i=0;i<nclasses;i++){
      nents=(Int_t)hgJetFragPL[k][i]->GetEntries();
      if(nents){
	//cout << "hJetFragPL " << i << ": " << nents << " entries in histogram " << endl;
	hgJetFragPL[k][i]->SetLineWidth(3);
	hgJetFragPL[k][i]->Write();
      }
    }
  }

  rootfile->cd();
  rootfile->mkdir("gFragTrans");
  rootfile->cd("gFragTrans");
  for(Int_t k=0;k<3;k++){
    for(Int_t i=0;i<nclasses;i++){
      nents=(Int_t)hgJetFragT[k][i]->GetEntries();
      if(nents){
	//cout << "hJetFragT " << i << ": " << nents << " entries in histogram" << endl;
	hgJetFragT[k][i]->SetLineWidth(3);
	hgJetFragT[k][i]->Write();
      }
    }
  }

  rootfile->cd();
  rootfile->mkdir("gDiFragLeadingPt");
  rootfile->cd("gDiFragLeadingPt");
  for(Int_t k=0;k<3;k++){
    for(Int_t i=0;i<nclasses;i++){
      nents=(Int_t)hgDiJetFragLeadingPt[k][i]->GetEntries();
      if(nents){
	//cout << "hJetFragLeadingPt " << i << ": "  << nents << " entries in histogram" << endl;
	hgDiJetFragLeadingPt[k][i]->SetLineWidth(3);
	hgDiJetFragLeadingPt[k][i]->Write();
      }
    }
  }

  rootfile->cd();
  rootfile->mkdir("gDiFragLong");
  rootfile->cd("gDiFragLong");
  for(Int_t k=0;k<3;k++){
    for(Int_t i=0;i<nclasses;i++){
      nents=(Int_t)hgDiJetFragL[k][i]->GetEntries();
      if(nents){
	//cout << "hgDiJetFragL " << i << ": " << nents << " entries in histogram " << endl;
	hgDiJetFragL[k][i]->SetLineWidth(3);
	hgDiJetFragL[k][i]->Write();
      }
    }
  }

  rootfile->cd();
  rootfile->mkdir("gDiFragPL");
  rootfile->cd("gDiFragPL");
  for(Int_t k=0;k<3;k++){
    for(Int_t i=0;i<nclasses;i++){
      nents=(Int_t)hgDiJetFragPL[k][i]->GetEntries();
      if(nents){
	//cout << "hgDJetFragPL " << i << ": " << nents << " entries in histogram " << endl;
	hgDiJetFragPL[k][i]->SetLineWidth(3);
	hgDiJetFragPL[k][i]->Write();
      }
    }
  }

  rootfile->cd();
  rootfile->mkdir("gDiFragTrans");
  rootfile->cd("gDiFragTrans");
  for(Int_t k=0;k<3;k++){
    for(Int_t i=0;i<nclasses;i++){
      nents=(Int_t)hgDiJetFragT[k][i]->GetEntries();
      if(nents){
	//cout << "hJetFragT " << i << ": " << nents << " entries in histogram" << endl;
	hgDiJetFragT[k][i]->SetLineWidth(3);
	hgDiJetFragT[k][i]->Write();
      }
    }
  }

  //reconstruction
  rootfile->cd();
  rootfile->mkdir("reconstruction");
  rootfile->cd("reconstruction");

  nents=(Int_t)hJetEtres->GetEntries();
  if(nents){
    cout << "hJetEtres " << nents << " entries in histogram" << endl;
    hJetEtres->SetLineWidth(3);
    hJetEtres->Write();
  }

  nents=(Int_t)hJetEtratio->GetEntries();
  if(nents){
    cout << "hJetEtratio " << nents << " entries in histogram" << endl;
    hJetEtratio->SetLineWidth(3);
    hJetEtratio->Write();
  }

  nents=(Int_t)hJetEtrestrue->GetEntries();
  if(nents){
    cout << "hJetEtrestrue " << nents << " entries in histogram" << endl;
    hJetEtrestrue->SetLineWidth(3);
    hJetEtrestrue->Write();
  }

  nents=(Int_t)hJetEtresTrigger->GetEntries();
  if(nents){
    cout << "hJetEtresTrigger " << nents << " entries in histogram" << endl;
    hJetEtresTrigger->SetLineWidth(3);
    hJetEtresTrigger->Write();
  }

  nents=(Int_t)hAxesDiffres->GetEntries();
  if(nents){
    cout << "hAxesDiffres " << nents << " entries in histogram" << endl;
    hAxesDiffres->SetLineWidth(3);
    hAxesDiffres->Write();
  }

  nents=(Int_t)hPhires->GetEntries();
  if(nents){
    cout << "hPhires " << nents << " entries in histogram" << endl;
    hPhires->SetLineWidth(3);
    hPhires->Write();
  }

  nents=(Int_t)hEtares->GetEntries();
  if(nents){
    cout << "hEtares " << nents << " entries in histogram" << endl;
    hEtares->SetLineWidth(3);
    hEtares->Write();
  }

  nents=(Int_t)hEtMonteres->GetEntries();
  if(nents){
    cout << "hEtMonteres " << nents << " entries in histogram" << endl;
    hEtMonteres->SetLineWidth(3);
    hEtMonteres->Write();
  }

  nents=(Int_t)hEtaMonteres->GetEntries();
  if(nents){
    cout << "hEtaMonteres " << nents << " entries in histogram" << endl;
    hEtaMonteres->SetLineWidth(3);
    hEtaMonteres->Write();
  }

  nents=(Int_t)hPhiMonteres->GetEntries();
  if(nents){
    cout << "hPhiMonteres " << nents << " entries in histogram" << endl;
    hPhiMonteres->SetLineWidth(3);
    hPhiMonteres->Write();
  }

  nents=(Int_t)hEtMonteratio->GetEntries();
  if(nents){
    cout << "hEtMonteratio " << nents << " entries in histogram" << endl;
    hEtMonteratio->SetLineWidth(3);
    hEtMonteratio->Write();
  }

  nents=(Int_t)hmJetEtres->GetEntries();
  if(nents){
    cout << "hmJetEtres " << nents << " entries in histogram" << endl;
    hmJetEtres->SetLineWidth(3);
    hmJetEtres->Write();
  }

  nents=(Int_t)hmJetEtratio->GetEntries();
  if(nents){
    cout << "hmJetEtratio " << nents << " entries in histogram" << endl;
    hmJetEtratio->SetLineWidth(3);
    hmJetEtratio->Write();
  }

  nents=(Int_t)hmJetEtrestrue->GetEntries();
  if(nents){
    cout << "hmJetEtrestrue " << nents << " entries in histogram" << endl;
    hmJetEtrestrue->SetLineWidth(3);
    hmJetEtrestrue->Write();
  }

  nents=(Int_t)hmAxesDiffres->GetEntries();
  if(nents){
    cout << "hmAxesDiffres " << nents << " entries in histogram" << endl;
    hmAxesDiffres->SetLineWidth(3);
    hmAxesDiffres->Write();
  }

  nents=(Int_t)hmPhires->GetEntries();
  if(nents){
    cout << "hmPhires " << nents << " entries in histogram" << endl;
    hmPhires->SetLineWidth(3);
    hmPhires->Write();
  }

  nents=(Int_t)hmEtares->GetEntries();
  if(nents){
    cout << "hmEtares " << nents << " entries in histogram" << endl;
    hmEtares->SetLineWidth(3);
    hmEtares->Write();
  }

  rootfile->cd();
  rootfile->mkdir("trigger");
  rootfile->cd("trigger");
  for(Int_t i=0;i<9;i++){
    hJetEttrigger[i]->Write();
    hJetEttrigger2[i]->Write();
  }
  hJetEttriggernorm->Write();

  //store event info
  rootfile->cd();
  TGraph *graph=new TGraph(nclasses+1);
  cout << "Good jets counters" << endl;
  for(Int_t i=0;i<nclasses;i++){
    cout << i << " " << nclGoodEvents[i] << endl;
    graph->SetPoint(i,i,nclGoodEvents[i]);
  }
  graph->SetPoint(nclasses,nclasses,nTotalEvents);
  cout << nclasses << " " << nTotalEvents << endl;
  graph->Write("ggoodevents");
  delete graph;
  graph=new TGraph(nclasses+1);
  for(Int_t i=0;i<nclasses;i++){
    //cout << i << " " << nclLeadEvents[i] << endl;
    graph->SetPoint(i,i,nclLeadEvents[i]);
  }
  graph->SetPoint(nclasses,nclasses,nTotalEvents);
  graph->Write("ggoodlead");
  delete graph;
  graph=new TGraph(nclasses+1);
  for(Int_t i=0;i<nclasses;i++){
    //cout << i << " " << nclDiEvents[i] << endl;
    graph->SetPoint(i,i,nclDiEvents[i]);
  }
  graph->SetPoint(nclasses,nclasses,nTotalEvents);
  graph->Write("ggooddilead");
  delete graph;

  //close the result file
  rootfile->Close();
  //

  delete hJetEt;
  delete hBackJetEt;
  delete hJetEtall;
  delete hJetEttrue;
  delete hBackJetEttrue;
  delete hJetEtalltrue;
  delete hJetEtTrigger;
  delete hJetEtUQTrigger;
  delete hJetEtvsTrigger;
  delete hJetEtvsUQTrigger;
  delete hJetEtvsEt;
  delete hAxesDiff;
  delete hPartPtDist;
  delete hPartEtaDist;
  delete hPartPhiDist;
  delete hPartPhiCorr;
  delete hPartACorr;
  delete hPartDiACorr;
  delete hPartDiPhiCorr;
  delete hJet1;
  delete hJet2;
  for(Int_t i=0;i<3;i++) delete hJettype[i];
  delete[] hJettype;
  delete hJetEtvsEll;
  delete hJetEtallvsEll;
  delete hJetZ;

  delete hJetEtres;
  delete hJetEtratio;
  delete hJetEtrestrue;
  delete hJetEtresTrigger;
  delete hAxesDiffres;
  delete hPhires;
  delete hEtares;
  delete hmJetEtres;
  delete hmJetEtratio;
  delete hmJetEtrestrue;
  delete hmAxesDiffres;
  delete hmPhires;
  delete hmEtares;
  delete hPhiMonteres;
  delete hEtaMonteres;
  delete hEtMonteres;
  delete hEtMonteratio;

  for(Int_t i=0;i<nclasses;i++) delete hJetLeadingPt[i];
  delete[] hJetLeadingPt;
  for(Int_t i=0;i<nclasses;i++) delete hJetFragLeadingPt[i];
  delete[] hJetFragLeadingPt;
  for(Int_t i=0;i<nclasses;i++) delete hJetLeadingPtDist[i];
  delete[] hJetLeadingPtDist;
  for(Int_t i=0;i<nclasses;i++) delete hJetFragL[i];
  delete[] hJetFragL;
  for(Int_t i=0;i<nclasses;i++) delete hJetFragPL[i];
  delete[] hJetFragPL;
  for(Int_t i=0;i<nclasses;i++) delete hJetFragT[i];
  delete[] hJetFragT;
  for(Int_t i=0;i<nclasses;i++) delete hJetFragPt[i];
  delete[] hJetFragPt;
  for(Int_t i=0;i<nclasses;i++) delete hJetN[i];
  delete[] hJetN;
  for(Int_t i=0;i<nclasses;i++) delete hJetMeanPt[i];
  delete[] hJetMeanPt;
  for(Int_t i=0;i<nclasses;i++) delete hPhiCorr[i];
  delete[] hPhiCorr;

  for(Int_t k=0;k<3;k++){
    for(Int_t i=0;i<nclasses;i++){
      delete hgJetFragLeadingPt[k][i];
      delete hgJetFragL[k][i];
      delete hgJetFragPL[k][i];
      delete hgJetFragT[k][i];
      delete hgJetFragPt[k][i];
    }
    delete[] hgJetFragLeadingPt[k];
    delete[] hgJetFragL[k];
    delete[] hgJetFragPL[k];
    delete[] hgJetFragT[k];
    delete[] hgJetFragPt[k];
  }
  delete[] hgJetFragLeadingPt;
  delete[] hgJetFragL;
  delete[] hgJetFragPL;
  delete[] hgJetFragT;
  delete[] hgJetFragPt;

  for(Int_t k=0;k<3;k++){
    for(Int_t i=0;i<nclasses;i++){
      delete hgDiJetFragLeadingPt[k][i];
      delete hgDiJetFragL[k][i];
      delete hgDiJetFragPL[k][i];
      delete hgDiJetFragT[k][i];
      delete hgDiJetFragPt[k][i];
    }
    delete[] hgDiJetFragLeadingPt[k];
    delete[] hgDiJetFragL[k];
    delete[] hgDiJetFragPL[k];
    delete[] hgDiJetFragT[k];
    delete[] hgDiJetFragPt[k];
  }
  delete[] hgDiJetFragLeadingPt;
  delete[] hgDiJetFragL;
  delete[] hgDiJetFragPL;
  delete[] hgDiJetFragT;
  delete[] hgDiJetFragPt;

  for(Int_t i=0;i<9;i++){
    delete hJetEttrigger[i];
    delete hJetEttrigger2[i];
  }
  delete[] hJetEttrigger;
  delete[] hJetEttrigger2;
  delete hJetEttriggernorm;
  delete rootfile;
}

//------------------------------------------------------------------

Float_t relphi(Float_t phi1, Float_t phi2)
{ //rel to phi1
  Float_t ret=TMath::Abs(phi1-phi2);
  if(ret>TMath::Pi()) ret=TMath::TwoPi()-ret;
  return ret;
}

Float_t addphi(Float_t phi1, Float_t phi2)
{
  Float_t addphi=phi1+phi2;
  if(addphi>TMath::TwoPi()) addphi-=TMath::TwoPi();
  else if(addphi<0) addphi+=TMath::TwoPi();
  return addphi;
}

Float_t diffphi(Float_t phi1, Float_t phi2)
{ //and move correlation to pi/2
  Float_t diffphi=TMath::Pi()/2+phi1-phi2;
  if(diffphi>TMath::TwoPi()) diffphi-=TMath::TwoPi();
  else if(diffphi<0) diffphi+=TMath::TwoPi();
  return diffphi;
}

Int_t eventindex(Float_t phi1, Float_t phi2)
{ //rel to phi1
  Int_t ret=0; //toward (300 - 60)
  Float_t dphi=addphi(phi1,-phi2);
  const Float_t slice=TMath::Pi()/3.; 
  if (dphi>1*slice && dphi < 2*slice) 
    ret=2; //transverse (60-120)
  else if (dphi>4*slice && dphi < 5*slice) 
    ret=3; //transverse (240-300)
  else if(dphi>=2*slice && dphi<=4*slice)
    ret=1; //away side (120-200)

  return ret;
}

void convert(Float_t pjet[4], Float_t &pt, Float_t &theta, Float_t &eta, Float_t &phi)
{
  pt=TMath::Sqrt(pjet[0]*pjet[0]+pjet[1]*pjet[1]);
  theta=TMath::ATan2(pt,pjet[2]);
  eta=-TMath::Log(TMath::Tan(theta/2));
  phi=TMath::Pi()+TMath::ATan2(-pjet[1],-pjet[0]);
}

