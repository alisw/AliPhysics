#if !defined(__CINT__) || defined(__MAKECINT__)

#include "../ITS/UPGRADE/AliITSUTrackerGlo.h"
#include "../ITS/UPGRADE/AliITSURecoDet.h"
#include "../ITS/UPGRADE/AliITSURecoLayer.h"
#include "TProfile.h"
#include "TH2.h"
#include "TAxis.h"
#include "TRandom.h"
#include "TFile.h"
#endif

/*
  to run: 
  .x LoadLibs.C 
  .x chkTracker.C 
  .x matBudget.C(trk)
*/

TProfile* histo[20]={0};
TH2F *rphi=0, *rphic=0;

void matBudget(AliITSUTrackerGlo* tracker, int nTest=100000, double phiMin=0,double phiMax=1.5708, int nbins=1000)
{
  AliITSURecoDet* its = tracker->GetITSInterface();
  int nlr = its->GetNLayers();
  double pnt0[3],pnt1[3],par[10];

  for (int ilr=0;ilr<nlr;ilr++) {
    histo[ilr] = new TProfile(Form("lr%d",ilr),
			      Form("lr%d act%d",ilr,its->GetLayer(ilr)->GetActiveID()),
			      nbins,phiMin,phiMax
			      );
  }
  //
  rphi  = new TH2F("rphi" ,"",nbins,phiMin,phiMax,its->GetRMax()*10,0,its->GetRMax());
  rphic = new TH2F("rphic","",nbins,phiMin,phiMax,its->GetRMax()*10,0,its->GetRMax());
  //
  float frac = 0.1;
  int npr=nTest*frac;
  for (int it=0;it<nTest;it++) {
    if ( (it%npr)==0 ) printf("Done %.1f%%\n",100*float(it)/nTest);
    double phi = phiMin + gRandom->Rndm()*(phiMax-phiMin);
    pnt1[0]=pnt1[1]=0;
    pnt1[2]=gRandom->Rndm()*10;
    double rmin=0;
    for (int ilr=0;ilr<nlr;ilr++) {
      double rmax = ilr<nlr-1 ? (its->GetLayer(ilr+1)->GetRMin()+its->GetLayer(ilr)->GetRMax())/2. : its->GetLayer(ilr)->GetRMax();
      if (it==0) {
	printf("Test of Lr%d in %f < R < %f\n",ilr,rmin,rmax);
	histo[ilr]->SetTitle(Form("%s %.3f<R<%.3f",histo[ilr]->GetTitle(),rmin,rmax));
      }
      for (int i=3;i--;) pnt0[i]=pnt1[i];
      pnt1[0] = rmax*TMath::Cos(phi);
      pnt1[1] = rmax*TMath::Sin(phi);
      pnt1[2] += rmax*1e-3; //to avoid strictly normal tracks
      tracker->MeanMaterialBudget(pnt0,pnt1,par);
      histo[ilr]->Fill(phi,par[1]);
      rmin = rmax;
    }   
  }
  //
  for (int ilr=0;ilr<nlr;ilr++) {
    if (!histo[ilr]->GetEntries()) continue;
    histo[ilr]->Scale(100);
    histo[ilr]->Fit("pol0","l");
    double rmsy = histo[ilr]->GetRMS(2);
    if (rmsy<1e-6) {
      double cst = histo[ilr]->GetMean(2)*100;
      histo[ilr]->SetMinimum(cst*0.8);      
      histo[ilr]->SetMaximum(cst*1.2);
    }
  }
  //
  TAxis* xax = rphi->GetXaxis();
  TAxis* yax = rphi->GetYaxis();
  int nbr = yax->GetNbins();
  for (int iphi=1;iphi<=nbins;iphi++) {
    double phi = xax->GetBinCenter(iphi);
    double rmax = 0;
    pnt1[0]=pnt1[1]=0;
    double avCum = 0;
    for (int ir=2;ir<nbr;ir++) {
      rmax = yax->GetBinCenter(ir);
      pnt1[0] = rmax*TMath::Cos(phi);
      pnt1[1] = rmax*TMath::Sin(phi);
      for (int i=2;i--;) pnt0[i]=pnt1[i];
      double av = 0;
      for (int it=0;it<100;it++) {
	pnt1[2] = pnt0[2] = gRandom->Rndm()*10;
	pnt1[2] += rmax*1e-3; //to avoid strictly normal tracks
	tracker->MeanMaterialBudget(pnt0,pnt1,par);
	av += par[1];
      }
      av /= 10;
      avCum += av;
      rphi->SetBinContent(iphi,ir,av);
      rphic->SetBinContent(iphi,ir,avCum);
    }
  }
  //

}
