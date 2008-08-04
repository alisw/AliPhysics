/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
    Comments to be written here:
    1. What do we calibrate.
    2. How to interpret results
    3. Simple example
    4. Analysis using debug streamers.



    3.Simple example
    // To make cosmic scan the user interaction neccessary
    //
    .x ~/UliStyle.C
    gSystem->Load("libANALYSIS");
    gSystem->Load("libTPCcalib");
    TFile fcalib("CalibObjects.root");
    TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
    AliTPCcalibCosmic * cosmic = ( AliTPCcalibCosmic *)array->FindObject("cosmicTPC");
    


*/



#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"

#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"

#include "AliTracker.h"
#include "AliMagFMaps.h"
#include "AliTPCCalROC.h"

#include "AliLog.h"

#include "AliTPCcalibCosmic.h"

#include "TTreeStream.h"
#include "AliTPCTracklet.h"

ClassImp(AliTPCcalibCosmic)


AliTPCcalibCosmic::AliTPCcalibCosmic() 
  :AliTPCcalibBase(),
   fGainMap(0),
   fHistNTracks(0),
   fClusters(0),
   fModules(0),
   fHistPt(0),
   fDeDx(0),
   fDeDxMIP(0),
   fMIPvalue(1), 
   fCutMaxD(5),        // maximal distance in rfi ditection
   fCutTheta(0.03),    // maximal distan theta
   fCutMinDir(-0.99)   // direction vector products
{  
  AliInfo("Default Constructor");    
}


AliTPCcalibCosmic::AliTPCcalibCosmic(const Text_t *name, const Text_t *title) 
  :AliTPCcalibBase(),
   fGainMap(0),
   fHistNTracks(0),
   fClusters(0),
   fModules(0),
   fHistPt(0),
   fDeDx(0),
   fDeDxMIP(0),
   fMIPvalue(1),
   fCutMaxD(5),        // maximal distance in rfi ditection
   fCutTheta(0.03),    // maximal distan theta
   fCutMinDir(-0.99)   // direction vector products
{  
  SetName(name);
  SetTitle(title);
  AliMagFMaps * field = new AliMagFMaps("dummy1", "dummy2",0,5,0);
  AliTracker::SetFieldMap(field, kTRUE);  

  fHistNTracks = new TH1F("ntracks","Number of Tracks per Event; number of tracks per event; number of tracks",501,-0.5,500.5);
  fClusters = new TH1F("signal","Number of Clusters per track; number of clusters per track n_{cl}; counts",160,0,160);
  fModules = new TH2F("sector","Acorde hits; z (cm); x(cm)",1200,-650,650,600,-700,700);
  fHistPt = new TH1F("Pt","Pt distribution; p_{T} (GeV); counts",2000,0,50);
  fDeDx = new TH2F("DeDx","dEdx; momentum p (GeV); TPC signal (a.u.)",500,0.01,100.,500,2.,1000);
  BinLogX(fDeDx);
  fDeDxMIP =  new TH1F("DeDxMIP","MIP region; TPC signal (a.u.);counts ",500,2.,1000);

  AliInfo("Non Default Constructor");  
  //
}

AliTPCcalibCosmic::~AliTPCcalibCosmic(){
  //
  //
  //
}





void AliTPCcalibCosmic::Process(AliESDEvent *event) {
  //
  //
  //
  if (!event) {
    Printf("ERROR: ESD not available");
    return;
  }  
  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!ESDfriend) {
   Printf("ERROR: ESDfriend not available");
   return;
  }
   
  FindPairs(event); // nearly everything takes place in find pairs...

  if (GetDebugLevel()>20) printf("Hallo world: Im here and processing an event\n");
  Int_t ntracks=event->GetNumberOfTracks(); 
  fHistNTracks->Fill(ntracks);
  if (ntracks==0) return;

}



void AliTPCcalibCosmic::Analyze() {

  fMIPvalue = CalculateMIPvalue(fDeDxMIP);

  return;

}



void AliTPCcalibCosmic::FindPairs(AliESDEvent *event) {
  //
  // Find cosmic pairs
  // 
  // Track0 is choosen in upper TPC part
  // Track1 is choosen in lower TPC part
  //
  if (GetDebugLevel()>20) printf("Hallo world: Im here\n");
  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  Int_t ntracks=event->GetNumberOfTracks(); 
  TObjArray  tpcSeeds(ntracks);
  if (ntracks==0) return;
  Double_t vtxx[3]={0,0,0};
  Double_t svtxx[3]={0.000001,0.000001,100.};
  AliESDVertex vtx(vtxx,svtxx);
  //
  //track loop
  //
  for (Int_t i=0;i<ntracks;++i) {
   AliESDtrack *track = event->GetTrack(i);
   fClusters->Fill(track->GetTPCNcls()); 
  
   const AliExternalTrackParam * trackIn = track->GetInnerParam();
   const AliExternalTrackParam * trackOut = track->GetOuterParam();
   if (!trackIn) continue;
   if (!trackOut) continue;
   
   AliESDfriendTrack *friendTrack = ESDfriend->GetTrack(i);
   TObject *calibObject;
   AliTPCseed *seed = 0;
   for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
     if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
   }
   if (seed) tpcSeeds.AddAt(seed,i);

   Double_t meanP = 0.5*(trackIn->GetP() + trackOut->GetP());
   if (seed && track->GetTPCNcls() > 80 + 60/(1+TMath::Exp(-meanP+5))) {
     fDeDx->Fill(meanP, seed->CookdEdxNorm(0.0,0.45,0,0,159,fGainMap));
     //
     if (meanP > 0.4 && meanP < 0.45) fDeDxMIP->Fill(seed->CookdEdxNorm(0.0,0.45,0,0,159,fGainMap));
     //
     if (GetDebugLevel()>0&&meanP>0.2&&seed->CookdEdxNorm(0.0,0.45,0,0,159,fGainMap)>300) {
       TFile *curfile = AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile();
       if (curfile) printf(">>> p+ in file: %s \t event: %i \t Number of ESD tracks: %i \n", curfile->GetName(), (int)event->GetEventNumberInFile(), (int)ntracks);
       if (track->GetOuterParam()->GetAlpha()<0) cout << " Polartiy: " << track->GetSign() << endl;
     }

   }

  }

  if (ntracks<2) return;
  //
  // Find pairs
  //
  for (Int_t i=0;i<ntracks;++i) {
    AliESDtrack *track0 = event->GetTrack(i);     
    // track0 - choosen upper part
    if (!track0) continue;
    if (!track0->GetOuterParam()) continue;
    if (track0->GetOuterParam()->GetAlpha()<0) continue;
    Double_t d1[3];
    track0->GetDirection(d1);    
    for (Int_t j=0;j<ntracks;++j) {
      if (i==j) continue;
      AliESDtrack *track1 = event->GetTrack(j);   
      //track 1 lower part
      if (!track1) continue;
      if (!track1->GetOuterParam()) continue;
      if (track1->GetOuterParam()->GetAlpha()>0) continue;
      //
      Double_t d2[3];
      track1->GetDirection(d2);
      
      AliTPCseed * seed0 = (AliTPCseed*) tpcSeeds.At(i);
      AliTPCseed * seed1 = (AliTPCseed*) tpcSeeds.At(j);
      if (! seed0) continue;
      if (! seed1) continue;
      Float_t dedx0 = seed0->CookdEdxNorm(0.05,0.55,0,0,159,fGainMap);
      Float_t dedx1 = seed1->CookdEdxNorm(0.05,0.55,0,0,159,fGainMap);
      //
      Float_t dedx0I = seed0->CookdEdxNorm(0.05,0.55,0,0,63,fGainMap);
      Float_t dedx1I = seed1->CookdEdxNorm(0.05,0.55,0,0,63,fGainMap);
      //
      Float_t dedx0O = seed0->CookdEdxNorm(0.05,0.55,0,64,159,fGainMap);
      Float_t dedx1O = seed1->CookdEdxNorm(0.05,0.55,0,64,159,fGainMap);
      //
      Float_t dir = (d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]);
      Float_t d0  = track0->GetLinearD(0,0);
      Float_t d1  = track1->GetLinearD(0,0);
      //
      // conservative cuts - convergence to be guarantied
      // applying before track propagation
      if (TMath::Abs(d0+d1)>fCutMaxD) continue;   // distance to the 0,0
      if (dir>fCutMinDir) continue;               // direction vector product
      Float_t bz = AliTracker::GetBz();
      Float_t dvertex0[2];   //distance to 0,0
      Float_t dvertex1[2];   //distance to 0,0 
      track0->GetDZ(0,0,0,bz,dvertex0);
      track1->GetDZ(0,0,0,bz,dvertex1);
      if (TMath::Abs(dvertex0[1])>250) continue;
      if (TMath::Abs(dvertex1[1])>250) continue;
      //
      //
      //
      Float_t dmax = TMath::Max(TMath::Abs(d0),TMath::Abs(d1));
      AliExternalTrackParam param0(*track0);
      AliExternalTrackParam param1(*track1);
      //
      // Propagate using Magnetic field and correct fo material budget
      //
      AliTracker::PropagateTrackTo(&param0,dmax+1,0.0005,3,kTRUE);
      AliTracker::PropagateTrackTo(&param1,dmax+1,0.0005,3,kTRUE);
      //
      // Propagate rest to the 0,0 DCA - z should be ignored
      //
      Bool_t b0 = param0.PropagateToDCA(&vtx,bz,1000);
      Bool_t b1 = param1.PropagateToDCA(&vtx,bz,1000);
      //      
      param0.GetDZ(0,0,0,bz,dvertex0);
      param1.GetDZ(0,0,0,bz,dvertex1);
      //
      Double_t xyz0[3];//,pxyz0[3];
      Double_t xyz1[3];//,pxyz1[3];
      param0.GetXYZ(xyz0);
      param1.GetXYZ(xyz1);
      Bool_t isPair = IsPair(&param0,&param1);
      //
      if (isPair) FillAcordeHist(track0);
      //
      if (fStreamLevel>0){
	TTreeSRedirector * cstream =  GetDebugStreamer();
	printf("My stream=%p\n",(void*)cstream);
	if (cstream) {
	  (*cstream) << "Track0" <<
	    "dir="<<dir<<               //  direction
	    "OK="<<isPair<<             //  will be accepted
	    "b0="<<b0<<                 //  propagate status
	    "b1="<<b1<<                 //  propagate status
	    "Orig0.=" << track0 <<      //  original track  0
	    "Orig1.=" << track1 <<      //  original track  1
	    "Tr0.="<<&param0<<          //  track propagated to the DCA 0,0
	    "Tr1.="<<&param1<<          //  track propagated to the DCA 0,0	   
	    "v00="<<dvertex0[0]<<       //  distance using kalman
	    "v01="<<dvertex0[1]<<       // 
	    "v10="<<dvertex1[0]<<       //
	    "v11="<<dvertex1[1]<<       // 
	    "d0="<<d0<<                 //  linear distance to 0,0
	    "d1="<<d1<<                 //  linear distance to 0,0
	    //
	    "x00="<<xyz0[0]<<
	    "x01="<<xyz0[1]<<
	    "x02="<<xyz0[2]<<
	    //
	    "x10="<<xyz1[0]<<
	    "x11="<<xyz1[1]<<
	    "x12="<<xyz1[2]<<
	    //
	    "Seed0.=" << seed0 <<       //  original seed 0
	    "Seed1.=" << seed1 <<       //  original seed 1
	    //
	    "dedx0="<<dedx0<<           //  dedx0 - all
	    "dedx1="<<dedx1<<           //  dedx1 - all
	    //
	    "dedx0I="<<dedx0I<<         //  dedx0 - inner ROC
	    "dedx1I="<<dedx1I<<         //  dedx1 - inner ROC
	    //
	    "dedx0O="<<dedx0O<<         //  dedx0 - outer ROC
	    "dedx1O="<<dedx1O<<         //  dedx1 - outer ROC
	    "\n";
	}
      }      
    }
  }  
}    




void  AliTPCcalibCosmic::FillAcordeHist(AliESDtrack *upperTrack) {

  // Pt cut to select straight tracks which can be easily propagated to ACORDE which is outside the magnetic field
  if (upperTrack->Pt() < 10 || upperTrack->GetTPCNcls() < 80) return;
    
  const Double_t AcordePlane = 850.; // distance of the central Acorde detectors to the beam line at y =0
  const Double_t roof = 210.5;       // distance from x =0 to end of magnet roof

  Double_t r[3];
  upperTrack->GetXYZ(r);
  Double_t d[3];
  upperTrack->GetDirection(d);
  Double_t x,z;
  z = r[2] + (d[2]/d[1])*(AcordePlane - r[1]);
  x = r[0] + (d[0]/d[1])*(AcordePlane - r[1]);
  
  if (x > roof) {
    x = r[0] + (d[0]/(d[0]+d[1]))*(AcordePlane+roof-r[0]-r[1]);
    z = r[2] + (d[2]/(d[0]+d[1]))*(AcordePlane+roof-r[0]-r[1]);
  }
  if (x < -roof) {
    x = r[0] + (d[0]/(d[1]-d[0]))*(AcordePlane+roof+r[0]-r[1]);	      
    z = r[2] + (d[2]/(d[1]-d[0]))*(AcordePlane+roof+r[0]-r[1]);
  } 

  fModules->Fill(z, x);
 
}



Long64_t AliTPCcalibCosmic::Merge(TCollection *li) {

  TIterator* iter = li->MakeIterator();
  AliTPCcalibCosmic* cal = 0;

  while ((cal = (AliTPCcalibCosmic*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibCosmic::Class())) {
      Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }
    
    fHistNTracks->Add(cal->GetHistNTracks());
    fClusters->Add(cal-> GetHistClusters());
    fModules->Add(cal->GetHistAcorde());
    fHistPt->Add(cal->GetHistPt());
    fDeDx->Add(cal->GetHistDeDx());
    fDeDxMIP->Add(cal->GetHistMIP());
  
  }
  
  return 0;
  
}


Bool_t  AliTPCcalibCosmic::IsPair(AliExternalTrackParam *tr0, AliExternalTrackParam *tr1){
  //
  //
  /*
  // 0. Same direction - OPOSITE  - cutDir +cutT    
  TCut cutDir("cutDir","dir<-0.99")
  // 1. 
  TCut cutT("cutT","abs(Tr1.fP[3]+Tr0.fP[3])<0.03")
  //
  // 2. The same rphi 
  TCut cutD("cutD","abs(Tr0.fP[0]+Tr1.fP[0])<5")
  //
  //
  //
  TCut cutPt("cutPt","abs(Tr1.fP[4]+Tr0.fP[4])<1&&abs(Tr0.fP[4])+abs(Tr1.fP[4])<10");  
  // 1/Pt diff cut
  */
  const Double_t *p0 = tr0->GetParameter();
  const Double_t *p1 = tr1->GetParameter();
  if (TMath::Abs(p0[3]+p1[3])>fCutTheta) return kFALSE;
  if (TMath::Abs(p0[0]+p1[0])>fCutMaxD)  return kFALSE;
  Double_t d0[3], d1[3];
  tr0->GetDirection(d0);    
  tr1->GetDirection(d1);       
  if (d0[0]*d1[0] + d0[1]*d1[1] + d0[2]*d1[2] >fCutMinDir) return kFALSE;
  //
  return kTRUE;  
}
 


Double_t AliTPCcalibCosmic::CalculateMIPvalue(TH1F * hist) {

  TF1 * funcDoubleGaus = new TF1("funcDoubleGaus", "gaus(0)+gaus(3)",0,1000);
  funcDoubleGaus->SetParameters(hist->GetEntries()*0.75,hist->GetMean()/1.3,hist->GetMean()*0.10,
				hist->GetEntries()*0.25,hist->GetMean()*1.3,hist->GetMean()*0.10);
  hist->Fit(funcDoubleGaus);
  Double_t MIPvalue = TMath::Min(funcDoubleGaus->GetParameter(1),funcDoubleGaus->GetParameter(4));

  delete funcDoubleGaus;

  return MIPvalue;

}




void AliTPCcalibCosmic::CalculateBetheParams(TH2F *hist, Double_t * initialParam) {

  return;

}


void AliTPCcalibCosmic::BinLogX(TH1 *h) {

  // Method for the correct logarithmic binning of histograms

  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *new_bins = new Double_t[bins + 1];
   
  new_bins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete new_bins;
  
}

AliExternalTrackParam *AliTPCcalibCosmic::Invert(AliExternalTrackParam *input)
{
  //
  // Invert paramerameter  - not covariance yet
  //
  AliExternalTrackParam *output = new AliExternalTrackParam(*input);
  Double_t * param = (Double_t*)output->GetParameter();
  param[0]*=-1;
  param[3]*=-1;
  param[4]*=-1;
  //
  return output;
}


/*

void AliTPCcalibCosmic::dEdxCorrection(){
  TCut cutT("cutT","abs(Tr1.fP[3]+Tr0.fP[3])<0.03");  // OK
  TCut cutD("cutD","abs(Tr0.fP[0]+Tr1.fP[0])<5");     // OK
  TCut cutPt("cutPt","abs(Tr1.fP[4]+Tr0.fP[4])<0.2&&abs(Tr0.fP[4])+abs(Tr1.fP[4])<10");
  TCut cutN("cutN","min(Orig0.fTPCncls,Orig1.fTPCncls)>70");
  TCut cutA=cutT+cutD+cutPt+cutN;


 .x ~/UliStyle.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
 AliXRDPROOFtoolkit tool; 
  TChain * chain = tool.MakeChain("cosmic.txt","Track0",0,1000000);
  chain->Lookup();

  .x ~/rootlogon.C
   gSystem->Load("libSTAT.so");

  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParam;
  TMatrixD covMatrix;
  
  chain->Draw("Tr0.fP[4]+Tr1.fP[4]","OK"+cutA);
  
  TString strFit;
  strFit+="(Tr0.fP[1]/250)++";
  strFit+="(Tr0.fP[1]/250)^2++";
  strFit+="(Tr0.fP[3])++";
  strFit+="(Tr0.fP[3])^2++";

  TString * ptParam = TStatToolkit::FitPlane(chain,"Tr0.fP[4]+Tr1.fP[4]", strFit.Data(),cutA, chi2,npoints,fitParam,covMatrix) 

strFit+="(Tr0.fP[1]/250)++";
strFit+="(Tr0.fP[1]/250)^2++";
strFit+="(Tr0.fP[3])++";
strFit+="(Tr0.fP[3])^2++";
strFit+="(Tr0.fP[1]/250)^2*Tr0.fP[3]++";
strFit+="(Tr0.fP[1]/250)^2*Tr0.fP[3]^2++";
//

strFit+="sign(Tr0.fP[1])++"
strFit+="sign(Tr0.fP[1])*(1-abs(Tr0.fP[1]/250))"
					    
TString * thetaParam = TStatToolkit::FitPlane(chain,"Tr0.fP[3]+Tr1.fP[3]", strFit.Data(),cutA, chi2,npoints,fitParam,covMatrix)



(-0.009263+(Tr0.fP[1]/250)*(-0.009693)+(Tr0.fP[1]/250)^2*(0.009062)+(Tr0.fP[3])*(0.002256)+(Tr0.fP[3])^2*(-0.052775)+(Tr0.fP[1]/250)^2*Tr0.fP[3]*(-0.020627)+(Tr0.fP[1]/250)^2*Tr0.fP[3]^2*(0.049030))

*/


/*
gSystem->Load("libANALYSIS");
gSystem->Load("libSTAT");
gSystem->Load("libTPCcalib");

TStatToolkit toolkit;
Double_t chi2;
TVectorD fitParam;
TMatrixD covMatrix;
Int_t npoints;
//
TCut cutT("cutT","abs(Tr1.fP[3]+Tr0.fP[3])<0.03");  // OK
TCut cutD("cutD","abs(Tr0.fP[0]+Tr1.fP[0])<5");     // OK
TCut cutPt("cutPt","abs(Tr1.fP[4]+Tr0.fP[4])<0.2&&abs(Tr0.fP[4])+abs(Tr1.fP[4])<10");
TCut cutN("cutN","min(Orig0.fTPCncls,Orig1.fTPCncls)>110");
TCut cutA="abs(norm-1)<0.3"+cutT+cutD+cutPt+cutN;



TTree * chain = Track0;


chain->SetAlias("norm","signalTot0.fElements[3]/signalTot1.fElements[3]");
//
chain->SetAlias("dr1","(signalTot0.fElements[1]/signalTot0.fElements[3])");
chain->SetAlias("dr2","(signalTot0.fElements[2]/signalTot0.fElements[3])");
chain->SetAlias("dr4","(signalTot0.fElements[4]/signalTot0.fElements[3])");
chain->SetAlias("dr5","(signalTot0.fElements[5]/signalTot0.fElements[3])");

TString fstring="";  
fstring+="dr1++";
fstring+="dr2++";
fstring+="dr4++";
fstring+="dr5++";
//
fstring+="dr1*dr2++";
fstring+="dr1*dr4++";
fstring+="dr1*dr5++";
fstring+="dr2*dr4++";
fstring+="dr2*dr5++";
fstring+="dr4*dr5++";



TString *strqdedx = toolkit.FitPlane(chain,"norm",fstring->Data(), cutA, chi2,npoints,fitParam,covMatrix,-1,0,200000);
  
chain->SetAlias("corQT",strqdedx->Data());

*/






