/// \file MakeAlignCalPad.C
///
/// \author marian.ivanov@cern.ch
///
/// Macro to create  alignment/distortion maps
/// As a input output of AliTPCcalibAlign is used.
///
/// Algorithm:
///
/// In the setup without the magnetic field the tracks are fitted using the linear track model.
/// ( in y-rphi coordinate the primary vertex is also used as a constrain)
/// Residuals (deltas0  between the track and clusters in Y and in z direction are filled in the 4 D histogram:
/// Delta: phi(180 bins): localX(53 bins): tan(phi): tan(theta)(10 bins)
///
/// Distortion map are extracted form the residual histograms as a mean value at each bin.
/// Linear fits are then performed for each pad - delta as function of theta
///
/// ~~~{.cpp}
/// Delta Ymeas = offsetY+slopeY*tan(theta)  
/// Delta Zmeas = offsetZ+slopeZ*tan(theta)
/// ~~~
/// 
/// Resulting residuals exported into the OCDB are:
///
/// ~~~{.cpp}
/// DeltaY  = offsetY
/// DeltaZ  = offsetZ
/// DeltaR  = slopeZ;
/// ~~~
/// 
/// Example usage:
///
/// make calpad+ make report ps file:
///
/// ~~~
/// aliroot -b -q ~/NimStyle.C ../ConfigCalibTrain.C\(119037\) $ALICE_ROOT/TPC/CalibMacros/MakeAlignCalPad.C\(1\)
/// ~~~
///
/// make only report ps file:
///
/// ~~~
/// aliroot -b -q ~/NimStyle.C $ALICE_ROOT/TPC/CalibMacros/MakeAlignCalPad.C\(3\)
/// ~~~
///
/// Making fit - iterative procedure - see below:
///
/// ~~~{.cpp}
/// gROOT->Macro("~/rootlogon.C");
/// gSystem->Load("libANALYSIS");
/// gSystem->Load("libSTAT");
/// gSystem->Load("libTPCcalib");
/// gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros -I$ALICE_ROOT/TPC/TPC -I$ALICE_ROOT/STAT");
/// .L $ALICE_ROOT/TPC/CalibMacros/MakeAlignCalPad.C+
/// // load OCDB
///
/// gROOT->Macro("../ConfigCalibTrain.C(119037)");
/// //2
/// InitTPCalign(); 
/// MakeFits();      // this is logn proceure 30 minutes
///
/// //UpdateOCDB(0,AliCDBRunRange::Infinity());
/// //
/// gROOT->Macro("~/NimStyle.C")
/// ~~~

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1D.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TVectorD.h"
#include "TTreeStream.h"
#include "TFile.h"
#include "AliTPCcalibAlign.h"
#include "AliTPCROC.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TROOT.h"
#include "TProfile.h"
#include "AliTPCPreprocessorOnline.h"
#include "AliTPCcalibDB.h"
#include "AliTPCkalmanAlign.h"
#include "TPostScript.h"
#include "TLegend.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "AliTPCExBEffectiveSector.h"
#include "AliTPCComposedCorrection.h"
//
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "TStatToolkit.h"
#endif


TObject  *palign=0;
TTree * treeP=0;
TTree * treeM=0;
TTree * tree0=0;
TTree * treePZ=0;
TTree * treeMZ=0;
TTree * tree0Z=0;
AliTPCROC *proc = AliTPCROC::Instance();
AliTPCExBEffectiveSector * geffCorr=0;

Double_t GetCorr(Double_t sector, Double_t localX, Double_t kZ,Int_t type);

void MakeFits();
void InitTPCalign();


void MakeAlignCalPad(Int_t mode){
  /// Make AlignCalpad and make report

  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTAT");
  gSystem->Load("libTPCcalib");  
  printf("Make report mode\t%d\n", mode);
  if (mode<1) { 
    InitTPCalign();
    MakeFits();
  }
  printf("Make report\n");
}


void InitTPCalign(){
  /// read the TPC alignment

  TFile fcalib("CalibObjects.root");
  TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
  if (array){
    palign = ( AliTPCcalibAlign *)array->FindObject("alignTPC");
  }else{
    palign = ( AliTPCcalibAlign *)fcalib.Get("alignTPC");
  }
  if (!palign){
     TFile fcalib2("TPCAlignObjects.root");
     palign = ( AliTPCcalibAlign *)fcalib2.Get("alignTPC");
  }
}


void MakeFits(){
  AliTPCcalibAlign * align =  (AliTPCcalibAlign *)palign;
  THnSparse * hdY = align->GetClusterDelta(0);
  THnSparse * hdZ = align->GetClusterDelta(1);
  AliTPCExBEffectiveSector::MakeResidualMap(hdY,"clusterDY.root");
  AliTPCExBEffectiveSector::MakeResidualMap(hdZ,"clusterDZ.root");
}

void FitdY(TTree * tree){
  ///

  AliTPCROC * roc = AliTPCROC::Instance();
  Double_t xquadrant = roc->GetPadRowRadii(36,53); 
  Double_t xIO       = ( roc->GetPadRowRadii(0,roc->GetNRows(0)-1)+roc->GetPadRowRadii(36,0))*0.5;
  
  //
  tree->SetAlias("diffS","(sector-int(sector)-0.5)");
  tree->SetAlias("diffQ",Form("(localX-%f)",xquadrant));
  tree->SetAlias("diffIO",Form("(localX-%f)",xIO));
  tree->SetAlias("iroc",Form("(localX<%f)",xIO));   // bool IROC
  tree->SetAlias("dqLR0",Form("(localX<%f)*(-1+(diffS<0)*2)",xIO));  //sign LeftRight
  tree->SetAlias("dqLR2",Form("(localX>%f)*(-1+(diffS<0)*2)",xquadrant));  //sign LeftRight up
  //
  tree->SetAlias("diffIFC","abs(localX-80)");
  tree->SetAlias("diffOFC","abs(250-localX)");
  tree->SetAlias("drift","(1-abs(kZ*localX)/250.)");
  //
  tree->SetAlias("dIFC2",Form("drift/(1+abs(diffIFC^2/(%f^2)))",10.));
  tree->SetAlias("dIFC",Form("drift/(1+abs(diffIFC/(%f)))",5.));
  tree->SetAlias("dOFC2",Form("drift/(1+abs(diffOFC^2/(%f^2)))",10.));
  tree->SetAlias("dOFC",Form("drift/(1+abs(diffOFC/(%f)))",5.));

 

  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD param;
  TMatrixD covar;
  TString  fstringG="";		 
  fstringG+="(iroc)++";           //1 iroc shift
  fstringG+="(dqLR0)++";
  fstringG+="(dqLR2)++";
  fstringG+="(diffIO*iroc)++";     // iroc rotation
  fstringG+="(diffIO*dqLR0)++";
  fstringG+="(diffQ*dqLR2)++";
  //
  fstringG+="(dIFC+dIFC2)++";            //9
  fstringG+="(diffS*(dIFC+dIFC2))++";
  fstringG+="(diffS^2*(dIFC+dIFC2))++";
  fstringG+="(dIFC-dIFC2)++";            //9
  fstringG+="(diffS*(dIFC-dIFC2))++";
  fstringG+="(diffS^2*(dIFC-dIFC2))++";
  //
  //
  fstringG+="(dOFC+dOFC2)++";            //9
  fstringG+="(diffS*(dOFC+dOFC2))++";
  fstringG+="(diffS^2*(dOFC+dOFC2))++";
  fstringG+="(dOFC-dOFC2)++";            //9
  fstringG+="(diffS*(dOFC-dOFC2))++";
  fstringG+="(diffS^2*(dOFC-dOFC2))++";
  //
  //
  //
  TTreeSRedirector *pcstream=new TTreeSRedirector("dyFit.root");
  {for (Int_t isec=0; isec<18; isec++){
      {for (Int_t iside=0; iside<=1; iside++){
	  TCut cutS=Form("abs(sector-%d.5)<0.5&&kZ*%d>0&&entries>2000",isec,(iside>0)?-1:1);
	  TString *strFitG = TStatToolkit::FitPlane(tree,"mean", fstringG.Data(),cutS, chi2,npoints,param,covar,-1,0, 10000000, kTRUE);
	  tree->SetAlias("fitG",strFitG->Data());
	  tree->Draw("mean-fitG",cutS,"");
	  //
	  printf("%d\t%d\t%f\n",iside,isec,TMath::Sqrt(chi2/npoints));
	  (*pcstream)<<"fitDy"<<
	    "iside="<<iside<<
	    "sec="<<isec<<
	    "chi2="<<chi2<<
	    "npoints="<<npoints<<
	    "p.="<<&param<<
	    "\n";
	}
      }
    }
  }
  delete pcstream;
}



void DumpDerivative(TH3 * his){
  ///

  Int_t nbins0=his->GetXaxis()->GetNbins();
  Int_t nbins1=his->GetYaxis()->GetNbins();
  Int_t nbins2=his->GetZaxis()->GetNbins();
  TTreeSRedirector *pcstream = new TTreeSRedirector("aaa.root");
  for (Int_t ibin0=1; ibin0<nbins0; ibin0++)
    for (Int_t ibin1=1; ibin1<nbins1; ibin1++)
      for (Int_t ibin2=1; ibin2<nbins2; ibin2++){
	Float_t x= his->GetXaxis()->GetBinCenter(ibin0);
	Float_t y= his->GetYaxis()->GetBinCenter(ibin1);
	Float_t z= his->GetZaxis()->GetBinCenter(ibin2);
	Float_t c= his->GetBinContent(ibin0,ibin1,ibin2);
	Float_t c0M=his->GetBinContent(ibin0-1,ibin1,ibin2);
	Float_t c0P=his->GetBinContent(ibin0+1,ibin1,ibin2);
	Float_t c1M=his->GetBinContent(ibin0,ibin1-1,ibin2);
	Float_t c1P=his->GetBinContent(ibin0,ibin1+1,ibin2);
	Float_t c2M=his->GetBinContent(ibin0,ibin1,ibin2-1);
	Float_t c2P=his->GetBinContent(ibin0,ibin1,ibin2+1);
	printf("%f\t%f\t%f\t%f\n",x,y,z,c);
	(*pcstream)<<"aaa"<<
	  "x="<<x<<
	  "y="<<y<<
	  "z="<<z<<
	  "c="<<c<<
	  "c0M="<<c0M<<
	  "c0P="<<c0P<<
	  "c1M="<<c1M<<
	  "c1P="<<c1P<<
	  "c2M="<<c2M<<
	  "c2P="<<c2P<<
	  "\n";
      }						
  delete pcstream;
}

Double_t GetCorr(Double_t sector, Double_t localX, Double_t kZ, Int_t type){
  /// calculate the correction at given position - check the geffCorr

  Double_t phi=sector*TMath::Pi()/9.;
  Double_t gx = localX*TMath::Cos(phi);
  Double_t gy = localX*TMath::Sin(phi);
  Double_t gz = localX*kZ;
  Int_t nsector=(gz>0) ? 0:18; 
  //
  //
  //
  Float_t distPoint[3]={gx,gy,gz};
  geffCorr->DistortPoint(distPoint, nsector);
  Double_t r0=TMath::Sqrt(gx*gx+gy*gy);
  Double_t r1=TMath::Sqrt(distPoint[0]*distPoint[0]+distPoint[1]*distPoint[1]);
  Double_t phi0=TMath::ATan2(gy,gx);
  Double_t phi1=TMath::ATan2(distPoint[1],distPoint[0]);
  if (type==0) return r1-r0;
  if (type==1) return (phi1-phi0)*r0;
  return phi1-phi0;
}



void LoadDistortionTrees(){
  /// Load distortion tree

  TFile *fp = new TFile("clusterDYPlus.root");
  TFile *fm = new TFile("clusterDYMinus.root");
  TFile *f0 = new TFile("clusterDY0.root");
  TFile *fpz= new TFile("clusterDZPlus.root");
  TFile *fmz= new TFile("clusterDZMinus.root");
  TFile *f0z=new TFile("clusterDZ0.root");
  treeP=(TTree*)fp->Get("delta");
  treeM=(TTree*)fm->Get("delta");
  tree0=(TTree*)f0->Get("delta");
  treePZ=(TTree*)fpz->Get("delta");
  treeMZ=(TTree*)fmz->Get("delta");
  tree0Z=(TTree*)f0z->Get("delta");
  tree0->AddFriend(treeP,"P");
  tree0->AddFriend(treeM,"M");
  tree0->AddFriend(treePZ,"PZ");
  tree0->AddFriend(treeMZ,"MZ");
  tree0->AddFriend(tree0Z,"Z");
  tree0->SetMarkerStyle(25);
  tree0->SetMarkerSize(0.4);
}

void UpdateEffSectorOCDB(){
  /// Incremeantal update ot the correction maps
  /// corrections on top of previous corrections

  TFile fp("clusterDYPlus.root");
  TFile fm("clusterDYMinus.root");
  TFile f0("clusterDY0.root");
  TFile f0z("clusterDZ0.root");
  
  TH3F *his3D0=(TH3F*)f0.Get("his3D");
  TH3F *his3DP=(TH3F*)fp.Get("his3D");
  TH3F *his3DM=(TH3F*)fm.Get("his3D");
  TH3F *his3DZ=(TH3F*)f0z.Get("his3D");
  TH3F *hisR=0;   //half of difference (Plus-Minus) scaled by c1
  hisR=(TH3F*)his3DP->Clone();
  hisR->Add(his3DM,-1);
  hisR->Scale(0.5);         
  hisR->Scale(1/0.3); // c1 factor to be added there //is sign correct?        
  //
  //
  AliTPCExBEffectiveSector *effSector = new  AliTPCExBEffectiveSector;
  effSector->SetName("EffSector");
  effSector->SetTitle("EffSector");
  effSector->fCorrectionRPhi=(TH3F*)his3D0->Clone();
  effSector->fCorrectionR=(TH3F*)hisR->Clone();
  effSector->fCorrectionZ=(TH3F*)his3DZ->Clone();
  geffCorr=effSector;
  AliTPCCorrection::AddVisualCorrection(geffCorr,0);
  //
  //
  //
  AliTPCcalibDB * calib = AliTPCcalibDB::Instance();
  TObjArray *arr = calib->GetTPCComposedCorrectionArray();
  {for (Int_t i=0; i<arr->GetEntries(); i++){
    AliTPCComposedCorrection * corr=(AliTPCComposedCorrection*)arr->At(i);
    if (!corr) continue;
    AliTPCExBEffectiveSector * corrSec=( AliTPCExBEffectiveSector *)corr->GetCorrections()->FindObject("EffSector");
    if (i==0) geffCorr=(AliTPCExBEffectiveSector *)corrSec;
    if (!corrSec) corr->GetCorrections()->Add(effSector->Clone());
    if (corrSec){
      if (corrSec->fCorrectionRPhi) corrSec->fCorrectionRPhi->Add(effSector->fCorrectionRPhi);
      else corrSec->fCorrectionRPhi=(TH3F*)effSector->fCorrectionRPhi->Clone();
      if (corrSec->fCorrectionR) corrSec->fCorrectionR->Add(effSector->fCorrectionR);
      else corrSec->fCorrectionR=(TH3F*)effSector->fCorrectionR->Clone();
      if (corrSec->fCorrectionZ) corrSec->fCorrectionR->Add(effSector->fCorrectionZ);
      else corrSec->fCorrectionZ=(TH3F*)effSector->fCorrectionZ->Clone();
    }
    if (i==0) AliTPCCorrection::AddVisualCorrection(corrSec,1);
    }}
// make OCDB entries
  TString userName=gSystem->GetFromPipe("echo $USER");
  TString date=gSystem->GetFromPipe("date");
  TString ocdbStorage="local:////";
  ocdbStorage+=gSystem->GetFromPipe("pwd");
  ocdbStorage+="/OCDB";
  //
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-27-04"); //root version
  metaData->SetComment(Form("Correction calibration. User: %s\n Data%s",userName.Data(),date.Data()));
  AliCDBId* id1=NULL;
  id1=new AliCDBId("TPC/Calib/Correction", 0, AliCDBRunRange::Infinity());
  AliCDBStorage* gStorage = 0;

  gStorage=AliCDBManager::Instance()->GetStorage((ocdbStorage+"Update").Data());
  gStorage->Put(arr, (*id1), metaData);


}


void DrawDiff(){
  TFile f0("../mergeField0/mean.root");
  TFile fP("../mergePlus/mean.root");
  TFile fM("../mergeMinus/mean.root");
  //
  TTree *itsdy0=(TTree*)f0.Get("ITSdy");
  TTree *itsdyP=(TTree*)fP.Get("ITSdy");
  TTree *itsdyM=(TTree*)fM.Get("ITSdy");
  TTree *vdy0=(TTree*)f0.Get("Vertexdy");
  TTree *vdyP=(TTree*)fP.Get("Vertexdy");
  TTree *vdyM=(TTree*)fM.Get("Vertexdy");
  itsdy0->SetMarkerStyle(25);
  itsdy0->SetMarkerSize(0.3);
  itsdy0->AddFriend(itsdyP,"P");
  itsdy0->AddFriend(itsdyM,"M");
  itsdy0->AddFriend(vdy0,"V");
  itsdy0->AddFriend(vdyP,"VP");
  itsdy0->AddFriend(vdyM,"VM");
  itsdy0->SetMarkerStyle(25);
  itsdy0->SetMarkerSize(0.4);

}





void MakePlotDeltaZ(){
  ///

  TCut cut="entries>500&&PZ.entries>500&&MZ.entries>500";
  TCanvas *cA = new TCanvas("deltaZA","deltaZA",900,700);
  TCanvas *cC = new TCanvas("deltaZC","deltaZC",900,700);
  TCanvas *cARef = new TCanvas("deltaZARef","deltaZARef",1000,800);
  TCanvas *cCRef = new TCanvas("deltaZCRef","deltaZCRef",1000,800);
  //
  TH1::AddDirectory(0);
  TH1 * his=0;
  cA->Divide(4,5);
  for (Int_t isec=0; isec<18; isec++){
    cA->cd(isec+1);
    TCut cutS=Form("abs(sector-0.5-%d)<0.1",isec);
    tree0->Draw("Z.mean*10.:localX:abs(kZ)",cutS+"entries>4000&&P.entries>200&&localX<250&&kZ>0&&abs(kZ)<1","colz");
    his=(TH1*)tree0->GetHistogram()->Clone();
    his->GetXaxis()->SetTitle("r (cm)");
    his->GetYaxis()->SetTitle("#Delta_{z} (mm)");
    his->GetZaxis()->SetTitle("tan($theta)");
    his->Draw("colz");
  }
  cA->cd(19);

  cC->Divide(4,5);
  for (Int_t isec=0; isec<18; isec++){
    cC->cd(isec+1);
    TCut cutS=Form("abs(sector-0.5-%d)<0.1",isec);
    tree0->Draw("Z.mean*10.:localX:abs(kZ)",cutS+"entries>4000&&P.entries>200&&localX<250&&kZ<0&&abs(kZ)<1","colz"); 
    his=(TH1*)tree0->GetHistogram()->Clone();
    his->GetXaxis()->SetTitle("r (cm)");
    his->GetYaxis()->SetTitle("#Delta_{z} (mm)");
    his->GetZaxis()->SetTitle("tan($theta)");
    his->Draw("colz");

  }

  cARef->Divide(4,5);
  for (Int_t isec=0; isec<18; isec++){
    cARef->cd(isec+1);
    TCut cutS=Form("abs(sector-0.5-%d)<0.1",isec);
    tree0->Draw("(PZ.mean-MZ.mean)*10.:localX:abs(kZ)",cutS+cut+"kZ>0&&abs(kZ)<1","colz");
    his=(TH1*)tree0->GetHistogram()->Clone();
    his->GetXaxis()->SetTitle("r (cm)");
    his->GetYaxis()->SetTitle("#Delta_{z} (mm)");
    his->GetZaxis()->SetTitle("tan($theta)");
    his->Draw("colz");

  }
  
  cCRef->Divide(4,5);
  for (Int_t isec=0; isec<18; isec++){
    cCRef->cd(isec+1);
    TCut cutS=Form("abs(sector-0.5-%d)<0.1",isec);
    tree0->Draw("(PZ.mean-MZ.mean)*10.:localX:abs(kZ)",cutS+cut+"kZ<0&&abs(kZ)<1","colz");
    his=(TH1*)tree0->GetHistogram()->Clone();
    his->GetXaxis()->SetTitle("r (cm)");
    his->GetYaxis()->SetTitle("#Delta_{z} (mm)");
    his->GetZaxis()->SetTitle("tan($theta)");
    his->Draw("colz");
  }

  TPostScript *ps = new TPostScript("distortionZ.ps", 112);
  ps->NewPage();
  cA->Draw();
  ps->NewPage();
  cA->Draw();
  ps->NewPage();
  cC->Draw();
  ps->NewPage();
  cARef->Draw();
  ps->NewPage();
  cCRef->Draw();
  ps->Close();
}




void MakeAlign(){
  /// make  sector alignment - using Kalman filter method -AliTPCkalmanAlign
  /// Combined information is used, mean residuals are minimized:
  ///
  /// 1. TPC-TPC sector alignment
  /// 2. TPC-ITS alignment
  /// 3. TPC vertex alignment

  TFile fcalib("../mergeField0/TPCAlignObjects.root");
  AliTPCcalibAlign * align = ( AliTPCcalibAlign *)fcalib.Get("alignTPC");
  TFile f0("../mergeField0/mean.root");

  //
  TTree *itsdy=(TTree*)f0.Get("ITSdy");
  TTree *itsdp=(TTree*)f0.Get("ITSdsnp");
  TTree *itsdz=(TTree*)f0.Get("ITSdz");
  TTree *itsdt=(TTree*)f0.Get("ITSdtheta");
  //
  TTree *vdy=(TTree*)f0.Get("Vertexdy");
  TTree *vds=(TTree*)f0.Get("Vertexdsnp");
  TTree *vdz=(TTree*)f0.Get("Vertexdz");
  TTree *vdt=(TTree*)f0.Get("Vertexdtheta");

  itsdy->AddFriend(itsdp,"Snp");
  itsdy->AddFriend(itsdz,"Z");
  itsdy->AddFriend(itsdt,"T");
  //
  itsdy->AddFriend(vdy,"V");
  itsdy->AddFriend(vds,"VSnp");
  itsdy->AddFriend(vdz,"VZ");
  itsdy->AddFriend(vdt,"VT");
  itsdy->SetMarkerStyle(25);
  itsdy->SetMarkerSize(0.4);

  TCut cutQ="entries>500&&abs(theta)<0.8&&abs(snp)<0.2";
  TH1F his1("hdeltaY1","hdeltaY1",100,-0.5,0.5);
  TMatrixD vecAlign(72,1);
  TMatrixD covAlign(72,72);
  AliTPCkalmanAlign::BookAlign1D(vecAlign,covAlign,0,0.0005);
  TVectorD vecITSY(72);
  TVectorD vecITSS(72);
  TVectorD vecVS(72);
  TVectorD vecITSTan(72);
  TVectorD vecVTan(72);
  {for (Int_t isec0=0; isec0<36; isec0++){
      Double_t phi0=(isec0%18+0.5)*TMath::Pi()/9.;
      if (phi0>TMath::Pi()) phi0-=TMath::TwoPi();
      Int_t iside0=(isec0%36<18)? 0:1;
      TCut cutSector=Form("abs(%f-phi)<0.14",phi0);
      TCut cutSide = (iside0==0)? "theta>0":"theta<0";
      itsdy->Draw("mean",cutQ+cutSector+cutSide);
      Double_t meanITSY=itsdy->GetHistogram()->GetMean()/83.6;
      vecITSY[isec0]=meanITSY;
      vecITSY[isec0+36]=meanITSY;
      itsdy->Draw("Snp.mean",cutQ+cutSector+cutSide);
      Double_t meanITSS=itsdy->GetHistogram()->GetMean();
      vecITSS[isec0]=meanITSS;
      vecITSS[isec0+36]=meanITSS;
      itsdy->Draw("VSnp.mean",cutQ+cutSector+cutSide);
      Double_t meanVS=itsdy->GetHistogram()->GetMean();
      vecVS[isec0]=meanVS;
      vecVS[isec0+36]=meanVS;
    }
  }
  AliTPCROC * roc = AliTPCROC::Instance();
  Double_t fXIROC = (roc->GetPadRowRadii(0,0)+roc->GetPadRowRadii(0,roc->GetNRows(0)-1))*0.5;
  Double_t fXOROC = (roc->GetPadRowRadii(36,0)+roc->GetPadRowRadii(36,roc->GetNRows(36)-1))*0.5;
  Double_t fXmiddle   = ( roc->GetPadRowRadii(0,0)+roc->GetPadRowRadii(36,roc->GetNRows(36)-1))*0.5;
  Double_t fXIO       = ( roc->GetPadRowRadii(0,roc->GetNRows(0)-1)+roc->GetPadRowRadii(36,0))*0.5;

  TTreeSRedirector *pcstream=new TTreeSRedirector("combAlign.root");
  
  {
    for (Int_t iter=0; iter<2; iter++){
    for (Int_t isec0=0; isec0<72; isec0++){
    for (Int_t isec1=0; isec1<72; isec1++){
      TH1 * his = align->GetHisto(AliTPCcalibAlign::kY,isec0,isec1);
      TH1 * hisPhi = align->GetHisto(AliTPCcalibAlign::kPhi,isec0,isec1);
      if (!his) continue;
      if (his->GetEntries()<100) continue;
      Double_t xref=fXIO;
      if (isec0<36&&isec1<36) xref=fXIROC;
      if (isec0>=36&&isec1>=36) xref=fXOROC;
      Double_t meanTPC=his->GetMean()/xref;
      Double_t meanTPCPhi=hisPhi->GetMean();
      Double_t meanITS0=vecITSY[isec0];
      Double_t meanITS1=vecITSY[isec1];
      Double_t meanITSS0=vecITSS[isec0];
      Double_t meanITSS1=vecITSS[isec1];
      Double_t meanVS0=vecVS[isec0];
      Double_t meanVS1=vecVS[isec1];
      AliTPCkalmanAlign::UpdateAlign1D(meanTPC, 0.001, isec0,isec1,  vecAlign,covAlign);
      AliTPCkalmanAlign::UpdateAlign1D(meanTPCPhi, 0.001, isec0,isec1,  vecAlign,covAlign);
      AliTPCkalmanAlign::UpdateAlign1D(meanITS1-meanITS0, 0.001, isec0,isec1,  vecAlign,covAlign);
      AliTPCkalmanAlign::UpdateAlign1D(meanITSS1-meanITSS0, 0.001, isec0,isec1,  vecAlign,covAlign);
      AliTPCkalmanAlign::UpdateAlign1D(meanVS1-meanVS0, 0.001, isec0,isec1,  vecAlign,covAlign);
      //printf("isec0\t%d\tisec1\t%d\t%f\t%f\t%f\n",isec0,isec1, meanTPC, meanITS0,meanITS1);
      Double_t kalman0= vecAlign(isec0,0);
      Double_t kalman1= vecAlign(isec1,0);
      if (iter>0) (*pcstream)<<"align"<<
	"iter="<<iter<<
	"xref="<<xref<<
	"isec0="<<isec0<<
	"isec1="<<isec1<<
	"mTPC="<<meanTPC<<
	"mTPCPhi="<<meanTPCPhi<<
	"mITS0="<<meanITS0<<
	"mITS1="<<meanITS1<<
	"mITSS0="<<meanITSS0<<
	"mITSS1="<<meanITSS1<<
	"mVS0="<<meanVS0<<
	"mVS1="<<meanVS1<<
	"k0="<<kalman0<<
	"k1="<<kalman1<<
	"\n";
    }          
    }
    }
  }
  pcstream->GetFile()->cd();
  vecAlign.Write("alignPhiMean");
  covAlign.Write("alingPhiCovar");
  delete pcstream;
  TFile f("combAlign.root");
  TTree * treeA = (TTree*)f.Get("align"); 
  treeA->SetMarkerStyle(25);
  treeA->SetMarkerSize(0.5);
}



