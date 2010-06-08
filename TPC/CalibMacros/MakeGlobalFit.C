 /*
  gROOT->Macro("~/rootlogon.C")
  gSystem->AddIncludePath("-I$ALICE_ROOT/STAT")
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC")
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros")
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
 
  .L $ALICE_ROOT/TPC/CalibMacros/MakeGlobalFit.C+
  MakeGlobalFit();
 
  rm matching.ps
  aliroot -b -q $ALICE_ROOT/TPC/CalibMacros/MakeGlobalFit.C | tee globalFit.log
  
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "THnSparse.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TCut.h"
#include "TH3.h"
#include "TH2F.h"
#include "TProfile3D.h"
#include "TMath.h" 
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TStatToolkit.h"
#include "TTreeStream.h"
#include "AliExternalTrackParam.h"
#include "AliESDfriend.h"
#include "AliTPCcalibTime.h"
#include "TROOT.h"
#include "AliXRDPROOFtoolkit.h"
#include "AliTPCCorrection.h"
#include "AliTPCExBTwist.h"
#include "AliTPCGGVoltError.h"
#include "AliTPCComposedCorrection.h"
#include "AliTPCExBConical.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "AliTrackerBase.h"
#include "AliTPCCalibGlobalMisalignment.h"
#include "AliTPCExBEffective.h"
#include "TEntryList.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBManager.h" 
#include "AliCDBStorage.h"
#include "AliTPCcalibDB.h"
#endif

const char* chDetName[5]={"ITS","TRD", "Vertex", "TOF","Laser"};
const char* chParName[5]={"rphi","z", "snp", "tan","1/pt"};
const char* chSideName[2]={"A side","C side"};
Bool_t enableDet[5]={1,1,1,0,1};      // detector  for fit  ITS=0, TRD=1, Vertex=2, Laser=4
Bool_t enableParam[5]={1,0,1,0,1};    //
Bool_t enableSign=kFALSE;  
Bool_t useEff0=kFALSE;
Bool_t useEffD=kFALSE;
Bool_t useEffR=kFALSE;
//
TChain *chain    = 0;
TChain *chainRef = 0;
Bool_t printMatrix=kFALSE;
TString *fitString=0;
//
// cuts
//
TCut cut="1";
TCut cutD="1";

void MakeChain();
void MakeAliases();
void MakeCuts();
void MakeFit(TCut cutCustom);
void MakeOCDBEntry(Int_t refRun);
TCanvas* DrawFitITS(const char *name);
TCanvas* DrawFitVertex(const char *name);
TCanvas*  DrawFitLaser(const char *cname);
TCanvas* DrawCorrdY();
TCanvas* DrawCorrdSnp();
TCanvas * DrawFitdY(const char *name);
TCanvas * DrawFitdSnp(const char *name);
void PrintMatch();
TCanvas * MakeComposedCorrection(const char * name);

void MakeAliases(){
  chain->SetAlias("isITS","(dtype==0||dtype==2)");
  chain->SetAlias("isTRD","(dtype==1)");
  chain->SetAlias("isLaser","(dtype==4)");
  chain->SetAlias("r","sqrt(gx*gx+gy*gy)");
  chain->SetAlias("r250","(sqrt(gx*gx+gy*gy)/250.)");
  chain->SetAlias("weight","((dtype==4)*rms*10+rms)");  // downscale laser
  chain->SetAlias("side","(1+(theta>0)*2)");
  chain->SetAlias("mdelta","(mean-R.mean-isLaser*((dRrec-R.dRrec)*tan(asin(snp))))");
  chain->SetAlias("drift","(1-abs(gz/250))");
  chain->SetAlias("r0","(r-160)/80");
  chain->SetAlias("delta","(0*gx)");
}


void MakeGlobalFit(){
  //
  //
  //
  gROOT->Macro("~/rootlogon.C");
  //gROOT->Macro("NimStyle.C");
  gSystem->AddIncludePath("-I$ALICE_ROOT/STAT");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");

  MakeChain();  
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  TPostScript *ps = new TPostScript("exbFit.eps", 112); 
  ps->NewPage();
  TCanvas *c=0;
  useEff0=kFALSE; useEffD=kFALSE; useEffR=kFALSE; 
  ps->NewPage();
  MakeFit("1");  
  //
  ps->NewPage();
  c=DrawFitdY("dY-Physical"); 
  c->Update();
  ps->NewPage();
  c->Update();

  ps->NewPage();
  c=DrawFitdSnp("dSnp-Physical");
  c->Update();
  
  ps->NewPage();
  c=DrawFitITS("ITS Physical");
  c->Update();
  
  ps->NewPage();
  c=DrawFitVertex("Vertex Physical");
  c->Update();
  
  ps->NewPage();
  c=DrawFitLaser("Laser Physical");
  c->Update();
  
  ps->NewPage();
  c=MakeComposedCorrection("Correction physical");
  c->Update();  

  //
  //
  //
  //
  useEff0=kTRUE; useEffD=kTRUE; useEffR=kTRUE; enableSign=kFALSE;
  ps->NewPage();
  MakeFit("1");
  ps->NewPage();
  c=DrawFitdY("dY-Physical+Effective ");
  c->Update();  
  
  ps->NewPage();
  c=DrawFitdSnp("dSnp-Physical+Effective ");
  c->Update();  

  ps->NewPage();
  c=DrawFitITS("ITS Physical+Effective");
  c->Update();  

  ps->NewPage();
  c=DrawFitVertex("Vertex Physical+Effective");
  c->Update();  

  ps->NewPage();
  c=DrawFitLaser("Laser Physical +Effective ");
  c->Update();  

  ps->NewPage();
  c=MakeComposedCorrection("Correction physical+Effective");
  c->Update();  
  //
  useEff0=kTRUE; useEffD=kTRUE; useEffR=kTRUE; enableSign=kTRUE; 
  ps->NewPage();
  MakeFit("1");
  //
  ps->NewPage();
  c=DrawFitdY("dY-Physical+Effective Sign");
  c->Update();  

  ps->NewPage();
  c=DrawFitdSnp("dSnp-Physical+Effective Sign");
  c->Update();  

  ps->NewPage();
  c=DrawFitITS("ITS Physical+Effective Sign");
  c->Update();  

  ps->NewPage();
  c=DrawFitVertex("Vertex Physical+Effective Sign");
  c->Update();  
  
  ps->NewPage();
  c=DrawFitLaser("Laser Physical +Effective Sign");
  c->Update();  
 
  ps->NewPage();
  c=MakeComposedCorrection("Correction physical+Effective Sign");
  c->Update();  
  
  ps->Close();
  delete ps;

  //
}

void MakeChain(){
  //
  //
  TH1::AddDirectory(0);
  TFile * f0 =0;      // file 0 field
  TFile * fp =0;      // file plus
  TFile * fm =0;      // file minus
  TTree * tree0=0;
  TTree * treeP=0;
  TTree * treeM=0;
  //
  chain    = new TChain("fit","fit");
  chainRef = new TChain("fit","fit");
  //
  //
  //
  for (Int_t idet=0; idet<5; idet++){
    for (Int_t ipar=0; ipar<5; ipar++){
      f0= TFile::Open(Form("../mergeField0/distortion%d_%d.root",idet,ipar));
      fp= TFile::Open(Form("../mergePlus/distortion%d_%d.root",idet,ipar));
      fm= TFile::Open(Form("../mergeMinus/distortion%d_%d.root",idet,ipar));
      tree0 = (f0) ? (TTree*)f0->Get("fit"):0;
      treeP = (fp) ? (TTree*)fp->Get("fit"):0;
      treeM = (fm) ? (TTree*)fm->Get("fit"):0;
      //
      if ( ipar==0 || ipar==2){
	if (tree0 && treeP){
	  chain->Add(Form("../mergePlus/distortion%d_%d.root",idet,ipar));
	  chainRef->Add(Form("../mergeField0/distortion%d_%d.root",idet,ipar));	
	}
	if (tree0 && treeM){
	  chain->Add(Form("../mergeMinus/distortion%d_%d.root",idet,ipar));
	  chainRef->Add(Form("../mergeField0/distortion%d_%d.root",idet,ipar));	
	}
      }
      if ( ipar==1 || ipar==3 || ipar==4){
	if (treeP && treeM){
	  chain->Add(Form("../mergePlus/distortion%d_%d.root",idet,ipar));
	  chainRef->Add(Form("../mergeMinus/distortion%d_%d.root",idet,ipar));	
	}
      }
    }
  }  
  chain->AddFriend(chainRef,"R");
  MakeAliases();
  MakeCuts();
}


void MakeCuts(){
  //
  //
  //
  TCut cutS="((rms>0&&R.rms>0&&entries>0&&R.entries>0))";         // statistic cuts
  TCut cutType="((dtype==R.dtype)&&(ptype==R.ptype))";            // corresponding types
  TCut cutOut="(ptype==0)*abs(mdelta)<(0.3+rms)||(ptype==0&&abs(mdelta*85)<(0.3+rms*85))";            // corresponding types
  //
  
  chain->Draw(">>outList",cutS+cutType+cutOut+"abs(snp)<0.5","entryList");
  TEntryList *elistOut = (TEntryList*)gDirectory->Get("outList");
  chain->SetEntryList(elistOut);
 
  
}


TMatrixD * MakeCorrelation(TMatrixD &matrix){
  //
  //
  //
  Int_t nrows = matrix.GetNrows();
  TMatrixD * mat = new TMatrixD(nrows,nrows);
  for (Int_t irow=0; irow<nrows; irow++)
    for (Int_t icol=0; icol<nrows; icol++){
      (*mat)(irow,icol)= matrix(irow,icol)/TMath::Sqrt(matrix(irow,irow)*matrix(icol,icol));
    }
  return mat;
}

void MakeDetCut(){  
  cutD=Form("((dtype==%d)||(dtype==%d)||(dtype==%d)||(dtype==%d)||(dtype==%d))",enableDet[0]?0:0,enableDet[1]?1:0,enableDet[2]?2:0,enableDet[3]?3:0,enableDet[4]?4:0);
  cutD+=Form("((ptype==%d)||(ptype==%d)||(ptype==%d)||(ptype==%d)||(ptype==%d))",enableParam[0]?0:0,enableParam[1]?1:0,enableParam[2]?2:0,enableParam[3]?3:0,enableParam[4]?4:0);
}

void MakeFit(TCut cutCustom){
  MakeDetCut();
  Int_t  npointsMax=30000000;
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD param;
  TMatrixD covar;
  //
  TString  fstring="";                       // magnetic part
  fstring+="(tX1-R.tX1)++";                  // twist X
  fstring+="(tY1-R.tY1)++";                   // twist Y
  fstring+="(sign(bz)*(tX1-R.tX1))++";       // twist X
  fstring+="(sign(bz)*(tY1-R.tY1))++";       // twist Y

  {if (enableDet[0]){
    fstring+="(isITS*shiftX)++";             // shift X - ITS
    fstring+="(isITS*shiftY)++";             // shift Y - ITS
    fstring+="(isITS*sign(bz)*shiftX)++";    // shift X - ITS
    fstring+="(isITS*sign(bz)*shiftY)++";    // shift Y - ITS
    }}
  {if (enableDet[1]){
    fstring+="(isTRD*shiftX)++";             // shift X - TRD
    fstring+="(isTRD*shiftY)++";             // shift Y - TRD
    fstring+="(isTRD*sign(bz)*shiftX)++";    // shift X - TRD
    fstring+="(isTRD*sign(bz)*shiftY)++";    // shift Y - TRD
    }}
  //
  if (enableDet[4]){
  }
  TString fstringCustom="";
  if (useEff0){
    fstring+="(exbT1-exb11-(R.exbT1-R.exb11))++";                // T1 adjustment
    fstring+="(exbT2-exb11-(R.exbT2-R.exb11))++";                // T2 adjustment
    fstringCustom+="(eff0_0_0-R.eff0_0_0)++";                  // effective correction constant part 
    fstringCustom+="(eff1_0_0-R.eff1_0_0)++";                  // 
  }
  if (useEffD){
    fstringCustom+="(eff0_0_1-R.eff0_0_1)++";                  // effective correction drift part
    fstringCustom+="(eff1_0_1-R.eff1_0_1)++";                  // 
    fstringCustom+="(eff0_0_2-R.eff0_0_2)++";                  // effective correction drift part 
    fstringCustom+="(eff1_0_2-R.eff1_0_2)++";                  // 
  }
  if (useEffR){
    fstringCustom+="(eff0_1_0-R.eff0_1_0)++";                  // effective correction radial part
    fstringCustom+="(eff1_1_0-R.eff1_1_0)++";                  // 
    fstringCustom+="(eff0_1_1-R.eff0_1_1)++";                  // effective correction radial part
    fstringCustom+="(eff1_1_1-R.eff1_1_1)++";                  // 
    fstringCustom+="(eff0_1_2-R.eff0_1_2)++";                  // effective correction radial part
    fstringCustom+="(eff1_1_2-R.eff1_1_2)++";                  // 
    fstringCustom+="(eff0_2_0-R.eff0_2_0)++";                  // effective correction radial part
    fstringCustom+="(eff1_2_0-R.eff1_2_0)++";                  // 
    fstringCustom+="(eff0_2_1-R.eff0_2_1)++";                  // effective correction radial part
    fstringCustom+="(eff1_2_1-R.eff1_2_1)++";                  // 
    fstringCustom+="(eff0_2_2-R.eff0_2_2)++";                  // effective correction radial part
    fstringCustom+="(eff1_2_2-R.eff1_2_2)++";                  // 
  }
  if (enableSign){
    if (useEff0){
      fstringCustom+="sign(bz)*(eff0_0_0-R.eff0_0_0)++";                  // effective correction constant part 
      fstringCustom+="sign(bz)*(eff1_0_0-R.eff1_0_0)++";                  // 
    }
    if (useEffD){
      fstringCustom+="sign(bz)*(eff0_0_1-R.eff0_0_1)++";                  // effective correction drift part
      fstringCustom+="sign(bz)*(eff1_0_1-R.eff1_0_1)++";                  // 
      fstringCustom+="sign(bz)*(eff0_0_2-R.eff0_0_2)++";                  // effective correction drift part 
      fstringCustom+="sign(bz)*(eff1_0_2-R.eff1_0_2)++";                  // 
    } 
    if (useEffR){
      fstringCustom+="sign(bz)*(eff0_1_0-R.eff0_1_0)++";                  // effective correction radial part
      fstringCustom+="sign(bz)*(eff1_1_0-R.eff1_1_0)++";                  // 
      fstringCustom+="sign(bz)*(eff0_1_1-R.eff0_1_1)++";                  // effective correction radial part
      fstringCustom+="sign(bz)*(eff1_1_1-R.eff1_1_1)++";                  // 
      fstringCustom+="sign(bz)*(eff0_1_2-R.eff0_1_2)++";                  // effective correction radial part
      fstringCustom+="sign(bz)*(eff1_1_2-R.eff1_1_2)++";                  // 
      fstringCustom+="sign(bz)*(eff0_2_0-R.eff0_2_0)++";                  // effective correction radial part
      fstringCustom+="sign(bz)*(eff1_2_0-R.eff1_2_0)++";                  // 
      fstringCustom+="sign(bz)*(eff0_2_1-R.eff0_2_1)++";                  // effective correction radial part
      fstringCustom+="sign(bz)*(eff1_2_1-R.eff1_2_1)++";                  // 
      fstringCustom+="sign(bz)*(eff0_2_2-R.eff0_2_2)++";                  // effective correction radial part
      fstringCustom+="sign(bz)*(eff1_2_2-R.eff1_2_2)++";                  // 
    }
  }
  //
  fstring+=fstringCustom;
  //
  TString *strDelta = TStatToolkit::FitPlane(chain,"mdelta:weight", fstring.Data(),cut+cutD+cutCustom, chi2,npoints,param,covar,-1,0, npointsMax, kTRUE);
  if (printMatrix) MakeCorrelation(covar)->Print();
  chain->SetAlias("delta",strDelta->Data());
  strDelta->Tokenize("++")->Print();
  fitString = strDelta;
  PrintMatch();
}





void PrintMatch(){
  //
  // Print detector matching info
  //
  for (Int_t ipar=0; ipar<5; ipar++){      
    for (Int_t idet=0; idet<5; idet++){
      Double_t mean0,rms0,mean1,rms1;
      Int_t entries = chain->Draw("delta-mdelta>>rhis",cut+Form("dtype==%d&&ptype==%d",idet,ipar),"goff");
      if (entries==0) continue;
      TH1 * his = (TH1*)(chain->GetHistogram()->Clone());
      mean1=his->GetMean();
      rms1 =his->GetRMS();
      delete his;
      //
      entries = chain->Draw("mdelta>>rhis",cut+Form("dtype==%d&&ptype==%d",idet,ipar),"goff");
      if (entries==0) continue;
      his = (TH1*)(chain->GetHistogram()->Clone());
      //
      mean0=his->GetMean();
      rms0 =his->GetRMS();
      
      printf("ptype==%s\tdtype==%s\tMean=%f -> %f\tRMS=%f -> %f\n",chParName[ipar],chDetName[idet], mean0,mean1,rms0,rms1);
      delete his;
    }
  }
}




TCanvas* DrawFitITS(const char *name){
  //
  //
  //
  TLegend *legend=0;
  TCanvas *canvas = new TCanvas(name,name,800,800);
  canvas->Divide(1,2);
  Int_t entries=0;
  const char * grname[10]={"A side (B=0.5T)","A side (B=-0.5T)", "C side (B=0.5T)", "C side (B=-0.5T)",0};

  TGraph *graphsdyITS[10];
  TGraph *graphsdyITSc[10];
  entries=chain->Draw("mdelta:phi",cut+"ptype==0&&dtype==0&&bz>0&&theta>0","goff");
  graphsdyITS[0]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("mdelta:phi",cut+cut+"ptype==0&&dtype==0&&bz<0&&theta>0","goff");
  graphsdyITS[1]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("mdelta:phi",cut+cut+"ptype==0&&dtype==0&&bz>0&&theta<0","goff");
  graphsdyITS[2]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("mdelta:phi",cut+cut+"ptype==0&&dtype==0&&bz<0&&theta<0","goff");
  graphsdyITS[3]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  //
  entries=chain->Draw("mdelta-delta:phi",cut+cut+"ptype==0&&dtype==0&&bz>0&&theta>0","goff");
  graphsdyITSc[0]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("mdelta-delta:phi",cut+cut+"ptype==0&&dtype==0&&bz<0&&theta>0","goff");
  graphsdyITSc[1]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("mdelta-delta:phi",cut+cut+"ptype==0&&dtype==0&&bz>0&&theta<0","goff");
  graphsdyITSc[2]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("mdelta-delta:phi",cut+cut+"ptype==0&&dtype==0&&bz<0&&theta<0","goff");
  graphsdyITSc[3]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  
  canvas->cd(1);
  legend = new TLegend(0.6,0.6,1.0,1.0,"#Delta R#phi TPC-ITS");
  for (Int_t i=0; i<4; i++){
    graphsdyITS[i]->SetMaximum(0.5);
    graphsdyITS[i]->SetMinimum(-0.5);
    graphsdyITS[i]->SetMarkerColor(i+1);
    graphsdyITS[i]->SetMarkerStyle(25+i);
    //    graphsdyITS[i]->SetName(grname[i]);
    graphsdyITS[i]->GetXaxis()->SetTitle("#phi");
    graphsdyITS[i]->GetYaxis()->SetTitle("#Delta R#Phi (cm)");
    if (i==0) graphsdyITS[i]->Draw("ap");
    graphsdyITS[i]->Draw("p");
    legend->AddEntry(graphsdyITS[i],grname[i]);
  }
  legend->Draw();
  canvas->cd(2);
  legend = new TLegend(0.6,0.6,1.0,1.0,"#Delta R#phi TPC-ITS corrected");
  for (Int_t i=0; i<4; i++){
    graphsdyITSc[i]->SetMaximum(0.5);
    graphsdyITSc[i]->SetMinimum(-0.5);
    graphsdyITSc[i]->SetMarkerColor(i+1);
    graphsdyITSc[i]->SetMarkerStyle(25+i);
    //    graphsdyITS[i]->SetName(grname[i]);
    graphsdyITSc[i]->GetXaxis()->SetTitle("#phi");
    graphsdyITSc[i]->GetYaxis()->SetTitle("#Delta R#Phi (cm)");
    if (i==0) graphsdyITS[i]->Draw("ap");
    graphsdyITSc[i]->Draw("p");
    legend->AddEntry(graphsdyITSc[i],grname[i]);
  }
  legend->Draw();
  return canvas;
}


TCanvas*  DrawFitLaser(const char *cname){
  //
  //
  //
  TH1::AddDirectory(0);
  TCut cutLaser=cut+"isLaser&&bz<0";
  TCanvas *canvas= new TCanvas(cname, cname,800,800);
  canvas->Divide(2,2);
  canvas->Draw();
  TLegend*legend=0;
  //
  TH1 *his[16]={0};
  {for (Int_t icorr=0; icorr<2; icorr++){
    for (Int_t iside=0; iside<2; iside++){
      canvas->cd(iside+1);
      for (Int_t ib=0; ib<4; ib++){    
	TCut cutB="1";
	if (iside==0) cutB=Form("(id==%d&&gz>0)",ib);
	if (iside==1) cutB=Form("(id==%d&&gz<0)",ib);
	//cutB.Print();
	if (icorr==0) chain->Draw("10*mdelta:r",cutLaser+cutB,"prof");
	if (icorr==1) chain->Draw("10*(mdelta-delta):r",cutLaser+cutB,"prof");
	his[icorr*8+iside*4+ib] = (TH1*)(chain->GetHistogram()->Clone());
	his[icorr*8+iside*4+ib]->SetName(Form("B%d%d%d",icorr,iside,ib));
	his[icorr*8+iside*4+ib]->SetTitle(Form("Bundle %d",ib));
	his[icorr*8+iside*4+ib]->SetMarkerColor(ib+1);
	his[icorr*8+iside*4+ib]->SetMarkerStyle(ib+25);
	his[icorr*8+iside*4+ib]->SetMarkerSize(0.4);
	his[icorr*8+iside*4+ib]->SetMaximum(3);
	his[icorr*8+iside*4+ib]->SetMinimum(-3);
	his[icorr*8+iside*4+ib]->GetXaxis()->SetTitle("r (cm)");
	his[icorr*8+iside*4+ib]->GetYaxis()->SetTitle("#Delta r#phi (mm)");
      }
    }
    }}
    //
  for (Int_t icorr=0; icorr<2; icorr++){
    for (Int_t iside=0; iside<2; iside++){
      canvas->cd(icorr*2+iside+1);
      legend = new TLegend(0.6,0.6,1.0,1.0,Form("#Delta R#phi Laser-%s",chSideName[iside]));
      for (Int_t ib=0; ib<4; ib++){    
	if (ib==0) his[icorr*2+iside*4+ib]->Draw();
	his[icorr*8+iside*4+ib]->Draw("same");
	legend->AddEntry(his[icorr*8+iside*4+ib]);
      }
      legend->Draw();
    }
  }
  return canvas;
}











TCanvas* DrawFitVertex(const char *name){
  //
  //
  //
  TLegend *legend=0;
  TCanvas *canvas = new TCanvas(name,name,800,800);
  canvas->Divide(1,2);
  Int_t entries=0;
  const char * grname[10]={"A side (B=0.5T)","A side (B=-0.5T)", "C side (B=0.5T)", "C side (B=-0.5T)",0};

  TGraph *graphsdyVertex[10];
  TGraph *graphsdyVertexc[10];
  entries=chain->Draw("mdelta:phi",cut+cut+"ptype==0&&dtype==2&&bz>0&&theta>0","goff");
  graphsdyVertex[0]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("mdelta:phi",cut+cut+"ptype==0&&dtype==2&&bz<0&&theta>0","goff");
  graphsdyVertex[1]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("mdelta:phi",cut+cut+"ptype==0&&dtype==2&&bz>0&&theta<0","goff");
  graphsdyVertex[2]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("mdelta:phi",cut+cut+"ptype==0&&dtype==2&&bz<0&&theta<0","goff");
  graphsdyVertex[3]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  //
  entries=chain->Draw("mdelta-delta:phi",cut+cut+"ptype==0&&dtype==2&&bz>0&&theta>0","goff");
  graphsdyVertexc[0]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("mdelta-delta:phi",cut+cut+"ptype==0&&dtype==2&&bz<0&&theta>0","goff");
  graphsdyVertexc[1]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("mdelta-delta:phi",cut+cut+"ptype==0&&dtype==2&&bz>0&&theta<0","goff");
  graphsdyVertexc[2]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("mdelta-delta:phi",cut+cut+"ptype==0&&dtype==2&&bz<0&&theta<0","goff");
  graphsdyVertexc[3]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  
  canvas->cd(1);
  legend = new TLegend(0.6,0.6,1.0,1.0,"#Delta R#phi TPC-Vertex");
  for (Int_t i=0; i<4; i++){
    graphsdyVertex[i]->SetMaximum(1);
    graphsdyVertex[i]->SetMinimum(-1);
    graphsdyVertex[i]->SetMarkerColor(i+1);
    graphsdyVertex[i]->SetMarkerStyle(25+i);
    //    graphsdyVertex[i]->SetName(grname[i]);
    graphsdyVertex[i]->GetXaxis()->SetTitle("#phi");
    graphsdyVertex[i]->GetYaxis()->SetTitle("#Delta R#Phi (cm)");
    if (i==0) graphsdyVertex[i]->Draw("ap");
    graphsdyVertex[i]->Draw("p");
    legend->AddEntry(graphsdyVertex[i],grname[i]);
  }
  legend->Draw();
  canvas->cd(2);
  legend = new TLegend(0.6,0.6,1.0,1.0,"#Delta R#phi TPC-Vertex corrected");
  for (Int_t i=0; i<4; i++){
    graphsdyVertexc[i]->SetMaximum(1);
    graphsdyVertexc[i]->SetMinimum(-1);
    graphsdyVertexc[i]->SetMarkerColor(i+1);
    graphsdyVertexc[i]->SetMarkerStyle(25+i);
    //    graphsdyVertex[i]->SetName(grname[i]);
    graphsdyVertexc[i]->GetXaxis()->SetTitle("#phi");
    graphsdyVertexc[i]->GetYaxis()->SetTitle("#Delta R#Phi (cm)");
    if (i==0) graphsdyVertex[i]->Draw("ap");
    graphsdyVertexc[i]->Draw("p");
    legend->AddEntry(graphsdyVertexc[i],grname[i]);
  }
  legend->Draw();
  return canvas;
}






TCanvas* DrawCorrdY(){
  
  TLegend *legend=0;
  TCanvas *canvas = new TCanvas("Corrections","Corrections",800,800);
  canvas->Divide(1,3);

  TGraph *graphsITS[10];
  TGraph *graphsTRD[10];
  TGraph *graphsVertex[10];
  Int_t entries=0;
  const char * grname[10]={"X twist (1 mrad)","Y twist (1 mrad)", "X shift (1 mm)", "Y shift (1 mm)", "drPhi (1 mm)"};
  //
  entries=chain->Draw("tX1:phi",cut+cut+"ptype==0&&dtype==0&&bz>0","goff");
  graphsITS[0]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("tY1:phi",cut+cut+"ptype==0&&dtype==0&&bz>0","goff");
  graphsITS[1]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("shiftX:phi",cut+cut+"ptype==0&&dtype==0&&bz>0","goff");
  graphsITS[2]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("shiftY:phi",cut+cut+"ptype==0&&dtype==0&&bz>0","goff");
  graphsITS[3]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("drPhiA+drPhiC:phi",cut+cut+"ptype==0&&dtype==0&&bz>0","goff");
  graphsITS[4]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  //  
  entries=chain->Draw("tX1:phi",cut+cut+"ptype==0&&dtype==1&&bz>0","goff");
  graphsTRD[0]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("tY1:phi",cut+cut+"ptype==0&&dtype==1&&bz>0","goff");
  graphsTRD[1]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("shiftX:phi",cut+cut+"ptype==0&&dtype==1&&bz>0","goff");
  graphsTRD[2]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("shiftY:phi",cut+cut+"ptype==0&&dtype==1&&bz>0","goff");
  graphsTRD[3]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("drPhiA+drPhiC:phi",cut+cut+"ptype==0&&dtype==1&&bz>0","goff");
  graphsTRD[4]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  //
  entries=chain->Draw("tX1:phi",cut+cut+"ptype==0&&dtype==2&&bz>0","goff");
  graphsVertex[0]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("tY1:phi",cut+cut+"ptype==0&&dtype==2&&bz>0","goff");
  graphsVertex[1]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("shiftX:phi",cut+cut+"ptype==0&&dtype==2&&bz>0","goff");
  graphsVertex[2]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("shiftY:phi",cut+cut+"ptype==0&&dtype==2&&bz>0","goff");
  graphsVertex[3]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("drPhiA+drPhiC:phi",cut+cut+"ptype==0&&dtype==2&&bz>0","goff");
  graphsVertex[4]=new TGraph(entries,chain->GetV2(), chain->GetV1());

  canvas->cd(1);
  //
  legend = new TLegend(0.6,0.6,1.0,1.0,"#Delta R#Phi TPC-Vertex");
  for (Int_t i=0; i<5; i++){
    graphsVertex[i]->SetMaximum(0.2);
    graphsVertex[i]->SetMinimum(-0.2);
    graphsVertex[i]->SetMarkerColor(i+1);
    graphsVertex[i]->SetMarkerStyle(25+i);
    graphsVertex[i]->SetName(grname[i]);
    graphsVertex[i]->GetXaxis()->SetTitle("#phi");
    graphsVertex[i]->GetYaxis()->SetTitle("#DeltaR#Phi (cm)");
    if (i==0) graphsVertex[i]->Draw("ap");
    graphsVertex[i]->Draw("p");
    legend->AddEntry(graphsITS[i],grname[i]);
  }
  legend->Draw();

  canvas->cd(2);
  legend = new TLegend(0.6,0.6,1.0,1.0,"#Delta R#Phi ITS-TPC");
  for (Int_t i=0; i<5; i++){
    graphsITS[i]->SetMaximum(0.2);
    graphsITS[i]->SetMinimum(-0.2);
    graphsITS[i]->SetMarkerColor(i+1);
    graphsITS[i]->SetMarkerStyle(25+i);
    graphsITS[i]->SetName(grname[i]);
    graphsITS[i]->GetXaxis()->SetTitle("#phi");
    graphsITS[i]->GetYaxis()->SetTitle("#Delta R#Phi (cm)");
    if (i==0) graphsITS[i]->Draw("ap");
    graphsITS[i]->Draw("p");
    legend->AddEntry(graphsITS[i],grname[i]);
  }  
  legend->Draw();
  //
  canvas->cd(3);
  legend = new TLegend(0.6,0.6,1.0,1.0,"#Delta R#Phi TPC-TRD");
  for (Int_t i=0; i<5; i++){
    graphsTRD[i]->SetMaximum(0.2);
    graphsTRD[i]->SetMinimum(-0.2);
    graphsTRD[i]->SetMarkerColor(i+1);
    graphsTRD[i]->SetMarkerStyle(25+i);
    graphsTRD[i]->SetName(grname[i]);
    graphsTRD[i]->GetXaxis()->SetTitle("#phi");
    graphsTRD[i]->GetYaxis()->SetTitle("#DeltaR#Phi (cm)");
    if (i==0) graphsTRD[i]->Draw("ap");
    graphsTRD[i]->Draw("p");
    legend->AddEntry(graphsITS[i],grname[i]);
  }
  legend->Draw();
  return canvas;
}


TCanvas * DrawCorrdSnp(){
  
  TLegend *legend=0;
  TCanvas *canvas = new TCanvas("Corrections dSnp","Corrections dSnp",800,800);
  canvas->Divide(1,3);

  TGraph *graphsITS[10];
  TGraph *graphsTRD[10];
  TGraph *graphsVertex[10];
  Int_t entries=0;
  const char * grname[10]={"X twist (1 mrad)","Y twist (1 mrad)", "X shift (1 mm)", "Y shift (1 mm)",0};
  //
  entries=chain->Draw("1000*tX1:phi",cut+cut+"ptype==2&&dtype==0&&bz>0","goff");
  graphsITS[0]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("1000*tY1:phi",cut+cut+"ptype==2&&dtype==0&&bz>0","goff");
  graphsITS[1]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("1000*shiftX:phi",cut+cut+"ptype==2&&dtype==0&&bz>0","goff");
  graphsITS[2]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("1000*shiftY:phi",cut+cut+"ptype==2&&dtype==0&&bz>0","goff");
  graphsITS[3]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  //  
  entries=chain->Draw("1000*tX1:phi",cut+cut+"ptype==2&&dtype==1&&bz>0","goff");
  graphsTRD[0]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("1000*tY1:phi",cut+cut+"ptype==2&&dtype==1&&bz>0","goff");
  graphsTRD[1]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("1000*shiftX:phi",cut+cut+"ptype==2&&dtype==1&&bz>0","goff");
  graphsTRD[2]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("1000*shiftY:phi",cut+cut+"ptype==2&&dtype==1&&bz>0","goff");
  graphsTRD[3]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  //
  entries=chain->Draw("1000*tX1:phi",cut+cut+"ptype==2&&dtype==2&&bz>0","goff");
  graphsVertex[0]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("1000*tY1:phi",cut+cut+"ptype==2&&dtype==2&&bz>0","goff");
  graphsVertex[1]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("1000*shiftX:phi",cut+cut+"ptype==2&&dtype==2&&bz>0","goff");
  graphsVertex[2]=new TGraph(entries,chain->GetV2(), chain->GetV1());
  entries=chain->Draw("1000*shiftY:phi",cut+cut+"ptype==2&&dtype==2&&bz>0","goff");
  graphsVertex[3]=new TGraph(entries,chain->GetV2(), chain->GetV1());

  canvas->cd(1);
  //
  legend = new TLegend(0.6,0.6,1.0,1.0,"#Delta sin(#phi) TPC-Vertex");
  for (Int_t i=0; i<4; i++){
    graphsVertex[i]->SetMaximum(2);
    graphsVertex[i]->SetMinimum(2);
    graphsVertex[i]->SetMarkerColor(i+1);
    graphsVertex[i]->SetMarkerStyle(25+i);
    graphsVertex[i]->SetName(grname[i]);
    graphsVertex[i]->GetXaxis()->SetTitle("#phi");
    graphsVertex[i]->GetYaxis()->SetTitle("#Deltasin(#phi) (mrad)");
    if (i==0) graphsVertex[i]->Draw("ap");
    graphsVertex[i]->Draw("p");
    legend->AddEntry(graphsITS[i],grname[i]);
  }
  legend->Draw();

  canvas->cd(2);
  legend = new TLegend(0.6,0.6,1.0,1.0,"#Delta sin(#phi) ITS-TPC");
  for (Int_t i=0; i<4; i++){
    graphsITS[i]->SetMaximum(1);
    graphsITS[i]->SetMinimum(-1);
    graphsITS[i]->SetMarkerColor(i+1);
    graphsITS[i]->SetMarkerStyle(25+i);
    graphsITS[i]->SetName(grname[i]);
    graphsITS[i]->GetXaxis()->SetTitle("#phi");
    graphsITS[i]->GetYaxis()->SetTitle("#Delta sin(#phi) (mrad)");
    if (i==0) graphsITS[i]->Draw("ap");
    graphsITS[i]->Draw("p");
    legend->AddEntry(graphsITS[i],grname[i]);
  }  
  legend->Draw();
  //
  canvas->cd(3);
  legend = new TLegend(0.6,0.6,1.0,1.0,"#Delta sin(#phi) TPC-TRD");
  for (Int_t i=0; i<4; i++){
    graphsTRD[i]->SetMaximum(1);
    graphsTRD[i]->SetMinimum(-1);
    graphsTRD[i]->SetMarkerColor(i+1);
    graphsTRD[i]->SetMarkerStyle(25+i);
    graphsTRD[i]->SetName(grname[i]);
    graphsTRD[i]->GetXaxis()->SetTitle("#phi");
    graphsTRD[i]->GetYaxis()->SetTitle("#Deltasin(#phi) (mrad)");
    if (i==0) graphsTRD[i]->Draw("ap");
    graphsTRD[i]->Draw("p");
    legend->AddEntry(graphsITS[i],grname[i]);
  }
  legend->Draw();
  return canvas;
 }

 
TCanvas * DrawFitdY(const char *name){
  //
  //
  //
  TH1::AddDirectory(0);
  TCanvas *canvas = new TCanvas(name,name,800,800);
  canvas->Divide(3,5);
  for (Int_t idet=0; idet<5; idet++){
    chain->SetMarkerColor(4);
    chain->SetMarkerStyle(25);
    chain->SetMarkerSize(0.3);
    chain->SetLineColor(2);
    //
    canvas->cd(idet*3+1);
    chain->Draw("mdelta:delta",cut+Form("ptype==0&&dtype==%d",idet),"");
    //
    canvas->cd(idet*3+2);
    chain->SetLineColor(2);
    chain->Draw("mdelta-delta",cut+Form("ptype==0&&dtype==%d",idet),"");
    chain->SetLineColor(1);
    chain->Draw("mdelta",cut+Form("ptype==0&&dtype==%d",idet),"same");
    //
    canvas->cd(idet*3+3);
    chain->SetMarkerColor(1);
    chain->Draw("mdelta:phi",cut+Form("ptype==0&&dtype==%d&&theta>0&&bz>0",idet),"");
    chain->SetMarkerColor(2);
    chain->Draw("mdelta-delta:phi",cut+Form("ptype==0&&dtype==%d&&theta>0&&bz>0",idet),"same");

  }
  return canvas;
}

TCanvas * DrawFitdSnp(const char *name){
  //
  //
  //
  TH1::AddDirectory(0);
  TCanvas *canvas = new TCanvas(name,name,800,800);
  canvas->Divide(3,5);
  {for (Int_t idet=0; idet<5; idet++){
    chain->SetMarkerColor(4);
    chain->SetMarkerStyle(25);
    chain->SetMarkerSize(0.3);
    chain->SetLineColor(2);
    //
    canvas->cd(idet*3+1);
    chain->Draw("1000*mdelta:1000*delta",cut+Form("ptype==2&&dtype==%d",idet),"");
    //
    canvas->cd(idet*3+2);
    chain->SetLineColor(2);
    chain->Draw("1000*(mdelta-delta)",cut+Form("ptype==2&&dtype==%d",idet),"");
    chain->SetLineColor(1);
    chain->Draw("1000*mdelta",cut+Form("ptype==2&&dtype==%d",idet),"same");
    //
    canvas->cd(idet*3+3);    
    chain->SetMarkerColor(1);
    chain->Draw("1000*(mdelta):phi",cut+Form("ptype==2&&dtype==%d&&theta>0&&bz>0",idet),"");
    chain->SetMarkerColor(2);
    chain->Draw("1000*(mdelta-delta):phi",cut+Form("ptype==2&&dtype==%d&&theta>0&&bz>0",idet),"same");
    }}
  return canvas;
}






TCanvas * MakeComposedCorrection(const char *name){
  
  TString  fit = chain->GetAlias("delta");
  TObjArray * array = fit.Tokenize("++");
  Int_t nfun=array->GetEntries();
  Double_t wt = 0.3 ; 
  Double_t T1 = 1.0;
  Double_t T2 = 1.0;
  //
  // sign independent correction
  //
  AliTPCExBEffective            *eff      = new  AliTPCExBEffective;
  AliTPCCalibGlobalMisalignment *shiftITS = new  AliTPCCalibGlobalMisalignment;
  AliTPCExBTwist                *twist    = new  AliTPCExBTwist;
  //
  // sign dependent correction
  //
  AliTPCExBEffective            *effS      = new  AliTPCExBEffective;
  AliTPCCalibGlobalMisalignment *shiftITSS = new  AliTPCCalibGlobalMisalignment;
  AliTPCExBTwist                *twistS    = new  AliTPCExBTwist;
  
  TMatrixD polA(100,4);
  TMatrixD polC(100,4);
  TMatrixD valA(100,1);
  TMatrixD valC(100,1);
  TMatrixD valAS(100,1);
  TMatrixD valCS(100,1);
  Int_t counterA=0;
  Int_t counterC=0;
  Int_t counterAS=0;
  Int_t counterCS=0;

  {
    for (Int_t i=1; i<nfun;i++){
      TObjString fitName=array->At(i)->GetName();
      TObjString fitVal= &(fitName.String()[1+fitName.String().Last('(')]);
      Double_t value= fitVal.String().Atof();
      if (fitName.String().Contains("sign")) continue;
      if (fitName.String().Contains("tX1")){
	twist->SetXTwist(value*0.001);
      }
      if (fitName.String().Contains("tY1")){
	twist->SetYTwist(value*0.001);
      }
      
      if (fitName.String().Contains("shiftX")&&fitName.String().Contains("isITS")){
	shiftITS->SetXShift(value*0.1);
      }
      if (fitName.String().Contains("shiftY")&&fitName.String().Contains("isITS")){
	shiftITS->SetYShift(value*0.1);
      }
      
      if (fitName.String().Contains("eff")){
	Int_t index=fitName.String().First("_")-1;
	Int_t side=atoi(&(fitName.String()[index]));
	Int_t px  =atoi(&(fitName.String()[index+2]));
	Int_t pd  =atoi(&(fitName.String()[index+4]));
	Int_t pp  =0;
	//printf("%s\t%d\t%d\t%d\t%d\t%f\n",fitName.GetName(),side,px,pd,pp, value);
	if (side==0){
	  polA(counterA,0)=0;
	  polA(counterA,1)=px;
	  polA(counterA,2)=pd;
	  polA(counterA,3)=pp;
	  valA(counterA,0)=value*0.1;
	  counterA++;
	}
	if (side==1){
	  polC(counterC,0)=0;
	  polC(counterC,1)=px;
	  polC(counterC,2)=pd;
	  polC(counterC,3)=pp;
	  valC(counterC,0)=value*0.1;
	  counterC++;
	}
      }
    }
  }
  polA.ResizeTo(counterA,4);
  polC.ResizeTo(counterC,4);
  valA.ResizeTo(counterA,1);
  valC.ResizeTo(counterC,1);
  eff->SetPolynoms(&polA,&polC);
  eff->SetCoeficients(&valA,&valC);
  eff->SetOmegaTauT1T2(wt,T1,T2);
  shiftITS->SetOmegaTauT1T2(wt,T1,T2);
  twist->SetOmegaTauT1T2(wt,T1,T2);
  //
  //
  //
  {
    counterAS=0;
    counterCS=0;
    for (Int_t i=1; i<nfun;i++){
      TObjString fitName=array->At(i)->GetName();
      TObjString fitVal= &(fitName.String()[1+fitName.String().Last('(')]);
      if (!fitName.String().Contains("sign")) continue;
      Double_t value= fitVal.String().Atof();
      if (fitName.String().Contains("tX1")){
	twistS->SetXTwist(value*0.001);
      }
      if (fitName.String().Contains("tY1")){
	twistS->SetYTwist(value*0.001);
      }
      
      if (fitName.String().Contains("shiftX")&&fitName.String().Contains("isITS")){
	shiftITSS->SetXShift(value*0.1);
      }
      if (fitName.String().Contains("shiftY")&&fitName.String().Contains("isITS")){
	shiftITSS->SetYShift(value*0.1);
      }
      
      if (fitName.String().Contains("eff")){
	Int_t index=fitName.String().First("_")-1;
	Int_t side=atoi(&(fitName.String()[index]));
	Int_t px  =atoi(&(fitName.String()[index+2]));
	Int_t pd  =atoi(&(fitName.String()[index+4]));
	Int_t pp  =0;
	//printf("%s\t%d\t%d\t%d\t%d\t%f\n",fitName.GetName(),side,px,pd,pp, value);
	if (side==0){
	  polA(counterAS,0)=0;
	  polA(counterAS,1)=px;
	  polA(counterAS,2)=pd;
	  polA(counterAS,3)=pp;
	  valAS(counterAS,0)=value*0.1;
	  counterAS++;
	}
	if (side==1){
	  polC(counterCS,0)=0;
	  polC(counterCS,1)=px;
	  polC(counterCS,2)=pd;
	  polC(counterCS,3)=pp;
	  valCS(counterCS,0)=value*0.1;
	  counterCS++;
	}
      }
    }
  }
  polA.ResizeTo(counterA,4);
  polC.ResizeTo(counterC,4);
  valA.ResizeTo(counterA,1);
  valC.ResizeTo(counterC,1);
  effS->SetPolynoms(&polA,&polC);
  effS->SetCoeficients(&valAS,&valCS);
  effS->SetOmegaTauT1T2(wt,T1,T2);
  shiftITSS->SetOmegaTauT1T2(wt,T1,T2);
  twistS->SetOmegaTauT1T2(wt,T1,T2);
  //
  // Make combined correction
  //

  TObjArray * corr0 = new TObjArray;
  TObjArray * corrP = new TObjArray;
  TObjArray * corrM = new TObjArray;
  AliTPCExBEffective            *eff0      = new  AliTPCExBEffective;
  AliTPCCalibGlobalMisalignment *shiftITSP = new  AliTPCCalibGlobalMisalignment;
  AliTPCExBTwist                *twistP    = new  AliTPCExBTwist;
  AliTPCExBEffective            *effP      = new  AliTPCExBEffective;
  AliTPCCalibGlobalMisalignment *shiftITSM = new  AliTPCCalibGlobalMisalignment;
  AliTPCExBTwist                *twistM    = new  AliTPCExBTwist;
  AliTPCExBEffective            *effM      = new  AliTPCExBEffective;
  //
  shiftITSP->SetXShift(shiftITS->GetXShift()+shiftITSS->GetXShift());    // shift due to the B field
  shiftITSP->SetYShift(shiftITS->GetXShift()+shiftITSS->GetYShift());                  
  shiftITSM->SetXShift(shiftITS->GetXShift()-shiftITSS->GetXShift());
  shiftITSM->SetYShift(shiftITS->GetXShift()-shiftITSS->GetYShift());
  //
  twistP->SetXTwist(twist->GetXTwist()+twistS->GetXTwist());     // twist between field  - both used
  twistP->SetYTwist(twist->GetYTwist()+twistS->GetYTwist());     //
  twistM->SetXTwist(twist->GetXTwist()-twistS->GetXTwist());
  twistM->SetYTwist(twist->GetYTwist()-twistS->GetYTwist());
  //
  effP->SetPolynoms(&polA,&polC);                                // effective correction
  effP->SetCoeficients(&valA,&valC);
  effM->SetPolynoms(&polA,&polC);
  effM->SetCoeficients(&valA,&valC);
  //
  eff0->SetPolynoms((TMatrixD*)&polA,(TMatrixD*)&polC);
  eff0->SetCoeficients((TMatrixD*)&valA,(TMatrixD*)&valC);
  //
  //
  //
  
  //
  corrP->AddLast(shiftITSP);
  corrP->AddLast(twistP);
  corrP->AddLast(effP);
  corrM->AddLast(shiftITSM);
  corrM->AddLast(twistM);
  corrM->AddLast(effM);
  corr0->AddLast(eff0);
  //

  AliTPCComposedCorrection *comp0= new AliTPCComposedCorrection ;
  AliTPCComposedCorrection *compP= new AliTPCComposedCorrection ;
  AliTPCComposedCorrection *compM= new AliTPCComposedCorrection ;
  comp0->SetCorrections((TObjArray*)(corr0));
  compP->SetCorrections((TObjArray*)(corrP));
  compM->SetCorrections((TObjArray*)(corrM));
  //
  comp0->SetOmegaTauT1T2(0*wt,T1,T2);
  compP->SetOmegaTauT1T2(wt,T1,T2);
  compM->SetOmegaTauT1T2(-wt,T1,T2);

  TFile f("correctionsExB.root","update");
  f.mkdir(name);
  f.cd(name);  
  printf("\nDump correction B=+0.5T:\n");
  compP->Print("da");
  printf("\nDump correction B=-0.5T:\n");
  compM->Print("da");
  //
  comp0->Write("ExB-Field0");  
  compP->Write("ExB-Bplus");
  compM->Write("ExB-Bminus");
  //
  //
  TCanvas *c = new TCanvas(name,name,800,800);
  c->Divide(4,2);
  c->cd(1);
  shiftITS->CreateHistoDRPhiinXY()->Draw("surf2");
  c->cd(2);
  twist->CreateHistoDRPhiinXY()->Draw("surf2");
  c->cd(3);
  eff->CreateHistoDRPhiinZR()->Draw("surf2");
  c->cd(5);
  shiftITSS->CreateHistoDRPhiinXY()->Draw("surf2");
  c->cd(6);
  twistS->CreateHistoDRPhiinXY()->Draw("surf2");
  c->cd(7);
  effS->CreateHistoDRPhiinZR()->Draw("surf2");
  //
  c->cd(4);
  compP->CreateHistoDRPhiinZR()->Draw("surf2");
  c->cd(8);
  compM->CreateHistoDRPhiinZR()->Draw("surf2");


  return c;
}


void MakeOCDBEntry(Int_t refRun){
  //
  // make a Correction OCDB entry
  // take the fit values writen in config file
  //
  //
  // 1. Read previous value used in calibration
  //    OCDB has to be initialized before
  
  gROOT->Macro(Form("ConfigCalibTrain.C(%d)",refRun));  // configuring calib db
  gROOT->LoadMacro("AddTaskTPCCalib.C");
  gROOT->ProcessLine(Form("ConfigOCDB(%d);",refRun));
  AliTPCCorrection *corr  = AliTPCcalibDB::Instance()->GetTPCComposedCorrection();
  //
  TFile f("correctionsExB.root");
  AliTPCComposedCorrection * corrP0=(AliTPCComposedCorrection *)f.Get("/Correction physical/ExB-Field0");
  AliTPCComposedCorrection * corrPP=(AliTPCComposedCorrection *)f.Get("/Correction physical/ExB-Bplus");
  AliTPCComposedCorrection * corrPM=(AliTPCComposedCorrection *)f.Get("/Correction physical/ExB-Bminus");
  //
  AliTPCComposedCorrection * corrPE0=(AliTPCComposedCorrection *)f.Get("/Correction physical+Effective/ExB-Field0");
  AliTPCComposedCorrection * corrPEP=(AliTPCComposedCorrection *)f.Get("/Correction physical+Effective/ExB-Bplus");
  AliTPCComposedCorrection * corrPEM=(AliTPCComposedCorrection *)f.Get("/Correction physical+Effective/ExB-Bminus");
  //
  corrP0->SetName("Field0 correction");
  corrPP->SetName("FieldP correction");
  corrPM->SetName("FieldM correction");
  corrPE0->SetName("Field0 correction");
  corrPEP->SetName("FieldP correction");
  corrPEM->SetName("FieldM correction");

  TObjArray corrPhysical;
  TObjArray corrPhysicalEffective;
  //
  // add the base correction
  corrP0->GetCorrections()->Add(corr);
  corrPP->GetCorrections()->Add(corr);
  corrPM->GetCorrections()->Add(corr);
  corrPE0->GetCorrections()->Add(corr);
  corrPEP->GetCorrections()->Add(corr);
  corrPEM->GetCorrections()->Add(corr);
  //
  corrPhysical.AddLast(corrP0);  corrPhysical.AddLast(corrPP);  corrPhysical.AddLast(corrPM); 
  corrPhysicalEffective.AddLast(corrPE0);  corrPhysicalEffective.AddLast(corrPEP);  corrPhysicalEffective.AddLast(corrPEM); 
  //
  // make OCDB entries
  TString userName=gSystem->GetFromPipe("echo $USER");
  TString date=gSystem->GetFromPipe("date");
  TString ocdbStorage="local:////";
  ocdbStorage+=gSystem->GetFromPipe("pwd")+"/OCDB";
  //
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-25-01"); //root version
  metaData->SetComment(Form("Correction calibration. User: %s\n Data%s",userName.Data(),date.Data()));
  AliCDBId* id1=NULL;
  id1=new AliCDBId("TPC/Calib/Correction", 0, AliCDBRunRange::Infinity());
  AliCDBStorage* gStorage = 0;

  gStorage=AliCDBManager::Instance()->GetStorage((ocdbStorage+"Physical").Data());
  gStorage->Put(&corrPhysical, (*id1), metaData);  
  gStorage = AliCDBManager::Instance()->GetStorage((ocdbStorage+"PhysicalEffective").Data());
  gStorage->Put(&corrPhysicalEffective, (*id1), metaData);  
}

