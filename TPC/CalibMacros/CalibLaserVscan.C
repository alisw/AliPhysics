/// \file CalibLaserVscan.C
///
/// 0. Make a calibration
/// 1. Make a laser scan list
///    e.g in TPC workscape
///
/// ~~~
/// find `pwd`/*/laserMean.root >laserScan.txt
/// ~~~
///
/// ~~~{.cpp}
/// // 2. Define a reference data 
/// rrunA=84469/; for a in `cat laserScan.txt`; do echo `pwd`/$rrunA/laserMean.root; done >laserScanRefA.txt
/// rrunC=84469; for a in `cat laserScan.txt`; do echo `pwd`/$rrunC/laserMean.root; done >laserScanRefC.txt
/// rrun=84469; for a in `cat laserScan.txt`; do echo `pwd`/$rrun/laserMean.root; done >laserScanRef.txt
///   									  // 
/// .x ~/rootlogon.C
/// gSystem->Load("libANALYSIS");
/// gSystem->Load("libTPCcalib"); 
/// gSystem->Load("libSTAT");
///
/// gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
/// gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
/// .L $ALICE_ROOT/TPC/CalibMacros/CalibLaserVscan.C+
///
/// AliXRDPROOFtoolkit tool;
/// chain = tool.MakeChainRandom("laserScan.txt","Mean",0,10200);
/// chain->Lookup();
/// chainRef = tool.MakeChain("laserScanRef.txt","Mean",0,10200);
/// chain->AddFriend(chainRef,"R") 
/// chainRefA = tool.MakeChain("laserScanRefA.txt","Mean",0,10200);
/// chain->AddFriend(chainRefA,"RA") 
/// chainRefC = tool.MakeChain("laserScanRefC.txt","Mean",0,10200);
/// chain->AddFriend(chainRefC,"RC") 
/// //
/// // MakeMeanBundle();
/// // SaveResult();   
/// //
/// ReadRunSetup();
/// ReadResult();
/// 
/// MakeAnalysisBeam();
/// TFile fbundle("scanDeltaBeam.root");
/// chain->Draw("mphi:GetValueBundle(id,1)","isOK&&GetValueBundle(id,0)>3&&LTr.fSide==0","")
/// ~~~

#include <fstream>
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include "TChain.h"
#include "TCut.h"
#include "TH1F.h"
#include "TObjArray.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TTreeStream.h"
#include "TLegend.h"
#include "AliTPCLaserTrack.h"
#include "AliTPCcalibDB.h"
#include "TStatToolkit.h"


TMatrixD matrixP4Corr(7,2);       // P4 correction for P0 - y position and P2 -snp(phi)
TMatrixD matrixP2RMSIROC(7,5);    // rms of second derivative IROC 
TMatrixD matrixP2RMSOROC(7,5);    // rms of second derivative OROC
// indexes - beam rod
//           rod 5 means all rods
//
TMatrixD matrixMeanBundle(336,6); // mean position/delta per bundle
TMatrixD matrixRMSBundle(336,6);  // rms of distribution
// indexes:
// 0 - number of entries
// 1 - dy
// 2 - dphi
// 3 - dp4
// 4 - dy0
// 5 - p0:p4 slope
//
map<int,TVectorD*> mapRunVoltage;             // run to voltage
// run#  -0 ggA 1- ggC 2- coA 3-coC 4-skA 5-skC
map<int,int> mapRunVgg;           // run to gating grid voltage
map<int,int> mapRunVskirt;        // 
map<int,int> mapRunVcover;        // 
TArrayI runlist(10000);
//
TChain * chain=0;
TObjArray apic;
//
//
TCut tA="isOK&&RA.isOK";
TCut tC="isOK&&RC.isOK";
TCut cA="eY.fElements<0.01&&RA.eY.fElements<0.01&&X.fElements>10&&RA.X.fElements>10";
TCut cC="eY.fElements<0.01&&RA.eY.fElements<0.01&&X.fElements>10&&RA.X.fElements>10";




Double_t GetValueBundle(Int_t id, Int_t type){
  ///

  return matrixMeanBundle(id,type);
}

Double_t GetP4Corr(Int_t id, Int_t value){
  AliTPCLaserTrack *ltrack =(AliTPCLaserTrack *) AliTPCLaserTrack::GetTracks()->At(id);
  Int_t beam = ltrack->GetBeam();
  return matrixP4Corr(beam,value);
}

Double_t GetDyBundle(Int_t id){
  Double_t dy = matrixMeanBundle(id,4);
  Double_t dx = matrixMeanBundle(id,5);
  AliTPCLaserTrack *ltrack =(AliTPCLaserTrack *) AliTPCLaserTrack::GetTracks()->At(id);
  Double_t p2 = ltrack->GetParameter()[2];
  Double_t delta = dy+dx*p2;
  return delta;
}

Double_t GetRMSBundle(Int_t id, Int_t type){
  ///

  return matrixRMSBundle(id,type);
}
Int_t GetVoltage(Int_t run, Int_t type){
  /// Get the voltage
  /// run#  -0 ggA 1- ggC 2- coA 3-coC 4-skA 5-skC

  TVectorD *runVoltage = mapRunVoltage[run];
  if (!runVoltage) return -1;
  return (*runVoltage)[type];
}



void SaveResult(){
  TFile f("laserData.root","recreate");
  matrixMeanBundle.Write("matrixMeanBundle");
  matrixRMSBundle.Write("matrixRMSBundle");
  matrixP4Corr.Write("matrixP4Corr");
  matrixP2RMSIROC.Write("matrixP2RMSIROC");
  matrixP2RMSOROC.Write("matrixP2RMSOROC");
  f.Close();
}
void ReadResult(){
  AliTPCLaserTrack::LoadTracks();
  TFile f("laserData.root");
  matrixMeanBundle.Read("matrixMeanBundle");
  matrixRMSBundle.Read("matrixRMSBundle");
  matrixP4Corr.Read("matrixP4Corr");
  matrixP2RMSIROC.Read("matrixP2RMSIROC");
  matrixP2RMSOROC.Read("matrixP2RMSOROC");
  f.Close();
}

void ReadRunSetup(){
  ifstream in;
  in.open("runSetup.txt");
  TString objfile;
  string line;
  TObjArray *arr = 0;
  Int_t counter=0;
  
  while(in.good()) {
    //in >> objfile;
    getline(in,line);
    objfile=line;
    printf("%s\n",objfile.Data());
    arr = objfile.Tokenize(" ");
    if (!arr) continue;
    if (arr->GetEntries()>=7){
      Int_t run = atoi(arr->At(0)->GetName()); 
      Int_t vcover = atoi(arr->At(0)->GetName()); 
      Int_t vskirt = atoi(arr->At(1)->GetName()); 
      Int_t vgg = atoi(arr->At(3)->GetName()); 
      TVectorD * vsetup = new TVectorD(6);
      for (Int_t iv=0; iv<6; iv++){
	(*vsetup)[iv] = atoi(arr->At(iv+1)->GetName());
      }
      mapRunVoltage[run]=vsetup;
      mapRunVgg[run]=vgg;
      mapRunVskirt[run]=vskirt;
      mapRunVcover[run]=vcover;
      runlist[counter]=run;
      counter++;
    }
    delete arr;    
  }
  runlist[counter]=-1;
}





void MakeAnalysisBeam(){
  //
  // 
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("scanDeltaBeam.root");
  TH1* phisP0=new TH1F("hhisP0","hhisP0",100,-5.,5.);
  TH1* phisP0X=new TH1F("hhisP0X","hhisP0X",100,-5.,5.);
  TH1* phisP2=new TH1F("hhisP2","hhisP2",100,-6,6);
  TH1* phisP4=new TH1F("hhisP4","hhisP4",100,-0.05,0.05);
  // second derivative
  TH1* phisP2IROC=new TH1F("hhisP2IROC","hhisPIROC",100,-2.,2.);
  TH1* phisP2OROC=new TH1F("hhisP2OROC","hhisPOROC",100,-2.,2.);
  //
  //
  //
  for (Int_t iside=0;iside<2;iside++){
    for (Int_t ibeam=0; ibeam<8; ibeam++){
      for (Int_t irun=0; irun<100; irun++){
	Int_t run =runlist[irun]; 
	if (run==-1) break;
	if (run==0) continue;
	Int_t vskirt = GetVoltage(run,1); 
	Int_t vcover = GetVoltage(run,2); 
	Int_t vgg    = GetVoltage(run,0); 
	//
	TCut cut = Form("isOK&&GetValueBundle(id,0)>3&&LTr.fSide==%d&&LTr.fBeam==%d&&run==%d", iside, ibeam,run);
	if (ibeam==7){
	  cut = Form("isOK&&GetValueBundle(id,0)>3&&LTr.fSide==%d&&run==%d", iside,run);
	}
	Int_t entries = chain->Draw("gp41-GetValueBundle(id,3)>>hhisP4",cut,"goff");	
	if (entries==0) continue;
	chain->Draw("10*(mphi-GetValueBundle(id,1)-GetP4Corr(id,0)*gp41)>>hhisP0",cut,"goff");
	chain->Draw("10*(mphi-GetDyBundle(id)-GetP4Corr(id,0)*gp41)>>hhisP0X",cut,"goff");
	chain->Draw("1000*(mphiP-GetValueBundle(id,2)-GetP4Corr(id,1)*gp41)>>hhisP2",cut,"goff");
	//
	printf("%d\t%d\t%d\t%f\t%f\t%d\n", iside, ibeam,run, phisP0X->GetMean(), phisP0X->GetRMS(), Int_t(phisP0X->GetEntries()));
	//
	chain->Draw("10*(mphi-GetValueBundle(id,1)-GetP4Corr(id,0)*gp41)>>hhisP0",cut,"goff");
	//
	chain->Draw("mPy2vP2In*100^2>>hhisP2IROC",cut,"goff");
	chain->Draw("mPy2vP2Out*100^2>>hhisP2OROC",cut,"goff");
	//
	//
	Double_t mp0  = phisP0->GetMean();
	Double_t mp0X = phisP0X->GetMean();
	Double_t mp2  = phisP2->GetMean();
	Double_t mp4  = phisP4->GetMean();
	Double_t mp2I = phisP2IROC->GetMean();
	Double_t mp2O = phisP2OROC->GetMean();
	Double_t sp0  = phisP0->GetRMS();
	Double_t sp0X = phisP0X->GetRMS();
	Double_t sp2  = phisP2->GetRMS();
	Double_t sp4  = phisP4->GetRMS();	  
	Double_t sp2I = phisP2IROC->GetRMS();
	Double_t sp2O = phisP2OROC->GetRMS();
	(*pcstream)<<"vScanBeam"<<
	  "side="<<iside<<
	  "run="<<run<<
	  "ibeam="<<ibeam<<
	  "vgg="<<vgg<<
	  "vskirt="<<vskirt<<
	  "vcover="<<vcover<<
	  "entries="<<entries<<
	  "mp0="<<mp0<<
	  "mp0X="<<mp0X<<
	  "mp2="<<mp2<<
	  "mp4="<<mp4<<
	  "mp2I="<<mp2I<<
	  "mp2O="<<mp2O<<
	  "sp0="<<sp0<<
	  "sp0X="<<sp0X<<
	  "sp2="<<sp2<<
	  "sp4="<<sp4<<
	  "sp2I="<<sp2I<<
	  "sp2O="<<sp2O<<
	  "\n";	
      }
    }
  }
  delete pcstream;
}




void MakeMeanBundle(){  
  ///

  AliTPCLaserTrack::LoadTracks();
  AliTPCLaserTrack *ltrack;
  TF1 * fp1 = 0;
  TCut ccut;
  TH1F * phisP0 = 0;
  TH1F * phisP2 = 0;
  TH1F * phisP4 = 0;
  //
  // get p4 corr
  //  
  for (Int_t ibeam=0; ibeam<7; ibeam++){
    Int_t entries0=chain->Draw("mphi:gp41",Form("isOK&&LTr.fBeam==%d",ibeam));
    TGraph gr0(entries0,chain->GetV2(),chain->GetV1());
    gr0.Draw();
    gr0.Fit("pol1","Q","Q");
    fp1 = gr0.GetFunction("pol1");
    matrixP4Corr(ibeam,0)=fp1->GetParameter(1);
    Int_t entries2=chain->Draw("mphiP:gp41",Form("isOK&&LTr.fBeam==%d",ibeam));
    TGraph gr2(entries2,chain->GetV2(),chain->GetV1());
    gr2.Draw();
    gr2.Fit("pol1","Q","Q");
    fp1 = gr2.GetFunction("pol1");
    matrixP4Corr(ibeam,1)=fp1->GetParameter(1);
  }
  //
  // get RMS of second derivative
  //
  phisP0= new TH1F("hisP0","hisP0",500,-1.,1.);
  for (Int_t irod=0; irod<5;irod++)
    for (Int_t ibeam=0; ibeam<7; ibeam++){
      //
      TCut cut=Form("isOK&&LTr.fBeam==%d&&LTr.fRod==%d",ibeam,irod);
      if (irod==5) cut=Form("isOK&&LTr.fBeam==%d",ibeam);
      chain->Draw("mPy2vP2In*100^2>>hisP0",cut);
      matrixP2RMSIROC(ibeam,irod) = phisP0->GetRMS();
      chain->Draw("mPy2vP2Out*100^2>>hisP0",cut);
      matrixP2RMSOROC(ibeam,irod) = phisP0->GetRMS();
  }
  delete phisP0;

  //
  // Get mean bundle position
  //
  for (Int_t id=0; id<336; id+=7){
    ltrack =(AliTPCLaserTrack *) AliTPCLaserTrack::GetTracks()->At(id);
    Int_t side   = ltrack->GetSide();
    Int_t rod    = ltrack->GetRod();
    Int_t bundle = ltrack->GetBundle();
    //    Int_t beam   = ltrack->GetBeam();
    ccut = Form("isOK&&LTr.fSide==%d&&LTr.fBundle==%d&&LTr.fRod==%d",side,bundle,rod);
    phisP0= new TH1F("hisP0","hisP0",500,-0.5,0.5);
    phisP2= new TH1F("hisP2","hisP2",500,-0.02,0.02);
    phisP4= new TH1F("hisP4","hisP4",500,-0.2,0.2);
    Int_t entries =chain->Draw("mphi",Form("isOK&&id==%d",id),"goff");
    chain->Draw("mphi-GetP4Corr(id,0)*gp41>>hisP0",ccut,"");
    chain->Draw("mphiP-GetP4Corr(id,1)*gp41>>hisP2",ccut,"");
    chain->Draw("gp41>>hisP4",ccut,"");
    Int_t entriesG =chain->Draw("mphi-GetP4Corr(id,0)*gp41:tan(asin(LTr.fP[2]))",ccut,"");
    if (entriesG<3) continue;
    TGraph gr(entriesG,chain->GetV2(),chain->GetV1());
    gr.Draw();
    gr.Fit("pol1","Q","Q");
    fp1 = gr.GetFunction("pol1");
    for (Int_t id2=0;id2<6; id2++){
      Int_t jd = id+id2;
      matrixMeanBundle(jd,0) = entries;
      matrixMeanBundle(jd,1) =phisP0->GetMean();
      matrixRMSBundle(jd,1) =phisP0->GetRMS();
      matrixMeanBundle(jd,2) =phisP2->GetMean();
      matrixRMSBundle(jd,2) =phisP2->GetRMS();
      matrixMeanBundle(jd,3) =phisP4->GetMean();
      matrixRMSBundle(jd,3) =phisP4->GetRMS();
      matrixMeanBundle(jd,4) = fp1->GetParameter(0);  // delta y
      matrixMeanBundle(jd,5) = fp1->GetParameter(1);  // delta x
    }
    //
    printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\n",id,entries,matrixMeanBundle(id,1),matrixMeanBundle(id,2),matrixMeanBundle(id,3),matrixMeanBundle(id,4),matrixMeanBundle(id,5));
    printf("%d\t%d\t%f\t%f\t%f\n",id,entries,matrixRMSBundle(id,1),matrixRMSBundle(id,2),matrixRMSBundle(id,3));
    delete phisP0;
    delete phisP2;
    delete phisP4;
  }
}








void MakeGraphsdY(){
  /// Make delta Y pictures from voltage scan

  TObjArray *aprofY = new TObjArray(14);
  for (Int_t ib=0;ib<14;ib++){
    TProfile *profY = new TProfile("py","py",100,0,150);
    chain->Draw("10*(mphi-GetValueBundle(id,1)):bz>>py",Form("isOK&&GetValueBundle(id,0)>3&&LTr.fSide==%d&&LTr.fBeam==%d",ib/7,ib%7),"prof");
    profY->SetName(Form("#DeltaY side %d beam %d",ib/7,ib%7));
    profY->SetTitle(Form("#DeltaY side %d beam %d",ib/7,ib%7));
    aprofY->AddAt(profY,ib);
    profY->SetMaximum(2.5);
    profY->SetMinimum(-2.5);
    profY->SetMarkerColor((ib%7)+1);
    profY->SetMarkerStyle((ib%7)+22);
    profY->SetMarkerSize(1.2);
    profY->SetXTitle("U_{gg} (V)");
    profY->SetYTitle("#Delta_{y} (mm)");
    profY->Fit("pol1");
  }
  TCanvas *cY = new TCanvas("deltaY","deltaY",900,600);
  cY->Divide(5,3);
  cY->Draw();
  {
    for (Int_t ib=0;ib<14;ib++){
      cY->cd(ib+1);   
      aprofY->At(ib)->Draw("p"); 
    }
  }
  cY->SaveAs("pic/deltaYlaserVscna.eps");
  cY->SaveAs("pic/deltaYlaserVscna.gif");
  apic.AddLast(cY);
}

void MakeGraphsdP2(){
  /// Make delta Y pictures from voltage scan

  TObjArray *aprofP2 = new TObjArray(14);
  for (Int_t ib=0;ib<14;ib++){
    TProfile *profP2 = new TProfile("pyP","pyP",100,0,150);
    chain->Draw("(mphiP-GetValueBundle(id,2)):bz>>pyP",Form("isOK&&GetValueBundle(id,0)>3&&LTr.fSide==%d&&LTr.fBeam==%d",ib/7,ib%7),"prof");
    profP2->SetName(Form("#Delta_{#phi} side %d beam %d",ib/7,ib%7));
    profP2->SetTitle(Form("#Delta_{#phi} side %d beam %d",ib/7,ib%7));
    aprofP2->AddAt(profP2,ib);
    profP2->SetMaximum(0.005);
    profP2->SetMinimum(-0.005);
    profP2->SetMarkerColor((ib%7)+1);
    profP2->SetMarkerStyle((ib%7)+22);
    profP2->SetMarkerSize(1.2);
    profP2->SetXTitle("U_{gg} (V)");
    profP2->SetYTitle("#Delta_{#phi}");
    profP2->Fit("pol1");
  }
  TCanvas *cP2 = new TCanvas("deltaP2","deltaP2",900,600);
  cP2->Divide(5,3);
  cP2->Draw();
  {
    for (Int_t ib=0;ib<14;ib++){
      cP2->cd(ib+1);   
      aprofP2->At(ib)->Draw("p"); 
    }
  }
  cP2->SaveAs("pic/deltaP2laserVscna.eps");
  cP2->SaveAs("pic/deltaP2laserVscna.gif");
  apic.AddLast(cP2);
}




void MakeGraphsP4(){
  /// Make delta Y pictures from voltage scan

  TObjArray *aprofP4 = new TObjArray(14);
  for (Int_t ib=0;ib<14;ib++){
    TProfile *profP4 = new TProfile("pp4","pp4",100,0,150);
    chain->Draw("(gp41-GetValueBundle(id,3)):bz>>pp4",Form("isOK&&GetValueBundle(id,0)>3&&LTr.fSide==%d&&LTr.fBeam==%d",ib/7,ib%7),"prof");
    profP4->SetName(Form("1/p_{t} side %d beam %d",ib/7,ib%7));
    profP4->SetTitle(Form("1/p_{t} side %d beam %d",ib/7,ib%7));
    aprofP4->AddAt(profP4,ib);
    profP4->SetMaximum(0.025);
    profP4->SetMinimum(-0.025);
    profP4->SetMarkerColor((ib%7)+1);
    profP4->SetMarkerStyle((ib%7)+22);
    profP4->SetMarkerSize(1.2);
    profP4->SetXTitle("U_{gg} (V)");
    profP4->SetYTitle("1/p_{t} (GeV/c)");
    profP4->Fit("pol1");
  }
  TCanvas *cY = new TCanvas("P4","P4",900,600);
  cY->Divide(5,3);
  cY->Draw();
  {
    for (Int_t ib=0;ib<14;ib++){
      cY->cd(ib+1);   
      aprofP4->At(ib)->Draw("p"); 
    }
  }
  cY->SaveAs("pic/m1ptlaserVscna.eps");
  cY->SaveAs("pic/m1ptlaserVscna.gif");
  apic.AddLast(cY);
}


void MakePlotsP2GG(TCut ucut){
  ///

  TFile fbundle("scanDeltaBeam.root");
  TTree * treeScan = (TTree*)fbundle.Get("vScanBeam");
  TGraph *graph[4];
  Int_t  mstyle[4]={22,23,24,25};
  Int_t  mcolor[4]={2,3,4,6};
  Int_t entries=0;
  entries = treeScan->Draw("10*sp2I/4:vgg","ibeam==7&&side==0"+ucut,"");
  graph[0]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  entries = treeScan->Draw("10*sp2O/4:vgg","ibeam==7&&side==0"+ucut,"");
  graph[1]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  entries = treeScan->Draw("10*sp2I/4:vgg","ibeam==7&&side==1"+ucut,"");
  graph[2]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  entries = treeScan->Draw("10*sp2O/4:vgg","ibeam==7&&side==1"+ucut,"");
  graph[3]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  
  for (Int_t i=0; i<4; i++){
    graph[i]->GetXaxis()->SetTitle("U_{gg} (V)");
    graph[i]->GetYaxis()->SetTitle("Sagita (mm)");
    graph[i]->SetMinimum(0);
    graph[i]->SetMaximum(2.5);
    graph[i]->SetMarkerStyle(mstyle[i]);
    graph[i]->SetMarkerSize(1);
    graph[i]->SetMarkerColor(mcolor[i]);
    if (i==0) graph[i]->Draw("ap");
    graph[i]->Draw("p");
  }
  TLegend *legend = new TLegend(0.45,0.70,0.9,0.9, "Voltage scan - Sagita on 50 cm");
  legend->AddEntry(graph[0],"IROC A side");
  legend->AddEntry(graph[1],"OROC A side");
  legend->AddEntry(graph[2],"IROC C side");
  legend->AddEntry(graph[3],"OROC C side");
  legend->Draw();
  gPad->SaveAs("pic/scansagita_Vgg.eps");
  gPad->SaveAs("pic/scansagita_Vgg.gif");
}





void MakePlotsP2Cover(TCut ucut){
  ///

  TFile fbundle("scanDeltaBeam.root");
  TTree * treeScan = (TTree*)fbundle.Get("vScanBeam");
  TGraph *graph[4];
  Int_t  mstyle[4]={22,23,24,25};
  Int_t  mcolor[4]={2,3,4,6};
  Int_t entries=0;
  entries = treeScan->Draw("10*sp2I/4:vcover","ibeam==7&&side==0"+ucut,"");
  graph[0]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  entries = treeScan->Draw("10*sp2O/4:vcover","ibeam==7&&side==0"+ucut,"");
  graph[1]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  entries = treeScan->Draw("10*sp2I/4:vcover","ibeam==7&&side==1"+ucut,"");
  graph[2]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  entries = treeScan->Draw("10*sp2O/4:vcover","ibeam==7&&side==1"+ucut,"");
  graph[3]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  
  for (Int_t i=0; i<4; i++){
    graph[i]->GetXaxis()->SetTitle("U_{cover} (V)");
    graph[i]->GetYaxis()->SetTitle("Sagita (mm)");
    graph[i]->SetMinimum(0.25);
    graph[i]->SetMaximum(0.7);
    graph[i]->SetMarkerStyle(mstyle[i]);
    graph[i]->SetMarkerSize(1);
    graph[i]->SetMarkerColor(mcolor[i]);
    if (i==0) graph[i]->Draw("ap");
    graph[i]->Draw("p");
  }
  TLegend *legend = new TLegend(0.45,0.70,0.9,0.9, "Voltage scan - Sagita on 50 cm");
  legend->AddEntry(graph[0],"IROC A side");
  legend->AddEntry(graph[1],"OROC A side");
  legend->AddEntry(graph[2],"IROC C side");
  legend->AddEntry(graph[3],"OROC C side");
  legend->Draw();
  gPad->SaveAs("pic/scansagita_Cover.eps");
  gPad->SaveAs("pic/scansagita_Cover.gif");
}

void MakePlotsP2Skirt(TCut ucut){
  ///

  TFile fbundle("scanDeltaBeam.root");
  TTree * treeScan = (TTree*)fbundle.Get("vScanBeam");
  TGraph *graph[4];
  Int_t  mstyle[4]={22,23,24,25};
  Int_t  mcolor[4]={2,3,4,6};
  Int_t entries=0;
  entries = treeScan->Draw("10*sp2I/4:vskirt","ibeam==7&&side==0"+ucut,"");
  graph[0]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  entries = treeScan->Draw("10*sp2O/4:vskirt","ibeam==7&&side==0"+ucut,"");
  graph[1]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  entries = treeScan->Draw("10*sp2I/4:vskirt","ibeam==7&&side==1"+ucut,"");
  graph[2]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  entries = treeScan->Draw("10*sp2O/4:vskirt","ibeam==7&&side==1"+ucut,"");
  graph[3]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  
  for (Int_t i=0; i<4; i++){
    graph[i]->GetXaxis()->SetTitle("U_{skirt} (V)");
    graph[i]->GetYaxis()->SetTitle("Sagita (mm)");
    graph[i]->SetMinimum(0.25);
    graph[i]->SetMaximum(0.7);
    graph[i]->SetMarkerStyle(mstyle[i]);
    graph[i]->SetMarkerSize(1);
    graph[i]->SetMarkerColor(mcolor[i]);
    if (i==0) graph[i]->Draw("ap");
    graph[i]->Draw("p");
  }
  TLegend *legend = new TLegend(0.45,0.70,0.9,0.9, "Voltage scan - Sagita on 50 cm");
  legend->AddEntry(graph[0],"IROC A side");
  legend->AddEntry(graph[1],"OROC A side");
  legend->AddEntry(graph[2],"IROC C side");
  legend->AddEntry(graph[3],"OROC C side");
  legend->Draw();
  gPad->SaveAs("pic/scansagita_Skirt.eps");
  gPad->SaveAs("pic/scansagita_Skirt.gif");
}




void MakePlotsdYGG(){
  ///

  TFile fbundle("scanDeltaBeam.root");
  TTree * treeScan = (TTree*)fbundle.Get("vScanBeam");
  TGraph *graph[4];
  Int_t  mstyle[4]={22,23,24,25};
  Int_t  mcolor[4]={2,3,4,6};
  Int_t entries=0;
  entries = treeScan->Draw("sp0X:vgg","ibeam==7&&side==0","");
  graph[0]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  entries = treeScan->Draw("sp0X:vgg","ibeam==7&&side==1","");
  graph[1]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  
  for (Int_t i=0; i<2; i++){
    graph[i]->GetXaxis()->SetTitle("U_{gg} (V)");
    graph[i]->GetYaxis()->SetTitle("#sigma_{r#phi} (mm)");
    graph[i]->SetMinimum(0);
    graph[i]->SetMaximum(1);
    graph[i]->SetMarkerStyle(mstyle[i]);
    graph[i]->SetMarkerSize(1);
    graph[i]->SetMarkerColor(mcolor[i]);
    if (i==0) graph[i]->Draw("ap");
    graph[i]->Draw("p");
  }
  TLegend *legend = new TLegend(0.45,0.70,0.9,0.9, "Voltage scan-#sigma_{r#phi}");
  legend->AddEntry(graph[0],"A side");
  legend->AddEntry(graph[1],"C side");
  legend->Draw();
  gPad->SaveAs("pic/scandy_Vgg.eps");
  gPad->SaveAs("pic/scandy_Vgg.gif");
}


void MakePlotsdYCover(TCut ucut){
  ///

  TFile fbundle("scanDeltaBeam.root");
  TTree * treeScan = (TTree*)fbundle.Get("vScanBeam");
  TGraph *graph[4];
  Int_t  mstyle[4]={22,23,24,25};
  Int_t  mcolor[4]={2,3,4,6};
  Int_t entries=0;
  entries = treeScan->Draw("sp0X:vcover","ibeam==7&&side==0"+ucut,"");
  graph[0]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  entries = treeScan->Draw("sp0X:vcover","ibeam==7&&side==1"+ucut,"");
  graph[1]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  
  for (Int_t i=0; i<2; i++){
    graph[i]->GetXaxis()->SetTitle("U_{cover} (V)");
    graph[i]->GetYaxis()->SetTitle("#sigma_{r#phi} (mm)");
    graph[i]->SetMinimum(0);
    graph[i]->SetMaximum(0.75);
    graph[i]->SetMarkerStyle(mstyle[i]);
    graph[i]->SetMarkerSize(1);
    graph[i]->SetMarkerColor(mcolor[i]);
    if (i==0) graph[i]->Draw("ap");
    graph[i]->Draw("p");
  }
  TLegend *legend = new TLegend(0.45,0.70,0.9,0.9, "Voltage scan-#sigma_{r#phi}");
  legend->AddEntry(graph[0],"A side");
  legend->AddEntry(graph[1],"C side");
  legend->Draw();
  gPad->SaveAs("pic/scandy_Vcover.eps");
  gPad->SaveAs("pic/scandy_Vcover.gif");
}

void MakePlotsdYSkirt(TCut ucut){
  ///

  TFile fbundle("scanDeltaBeam.root");
  TTree * treeScan = (TTree*)fbundle.Get("vScanBeam");
  TGraph *graph[4];
  Int_t  mstyle[4]={22,23,24,25};
  Int_t  mcolor[4]={2,3,4,6};
  Int_t entries=0;
  entries = treeScan->Draw("sp0X:vskirt","ibeam==7&&side==0"+ucut,"");
  graph[0]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  entries = treeScan->Draw("sp0X:vskirt","ibeam==7&&side==1"+ucut,"");
  graph[1]= new TGraph(entries, treeScan->GetV2(), treeScan->GetV1());
  
  for (Int_t i=0; i<2; i++){
    graph[i]->GetXaxis()->SetTitle("U_{skirt} (V)");
    graph[i]->GetYaxis()->SetTitle("#sigma_{r#phi} (mm)");
    graph[i]->SetMinimum(0);
    graph[i]->SetMaximum(0.75);
    graph[i]->SetMarkerStyle(mstyle[i]);
    graph[i]->SetMarkerSize(1);
    graph[i]->SetMarkerColor(mcolor[i]);
    if (i==0) graph[i]->Draw("ap");
    graph[i]->Draw("p");
  }
  TLegend *legend = new TLegend(0.45,0.70,0.9,0.9, "Voltage scan-#sigma_{r#phi}");
  legend->AddEntry(graph[0],"A side");
  legend->AddEntry(graph[1],"C side");
  legend->Draw();
  gPad->SaveAs("pic/scandy_Vskirt.eps");
  gPad->SaveAs("pic/scandy_Vskirt.gif");
}



void GetOptimalSetting(){
  ///

  TFile fbundle("scanDeltaBeam.root");
  TTree * treeScan = (TTree*)fbundle.Get("vScanBeam");
  //
  Float_t meansp0X=0;
  Float_t meansp2I=0;
  Float_t meansp2O=0; 
  TH1F * his = new TH1F("hisDelta","hisDelta",300,0,1);
  treeScan->Draw("sp0X>>hisDelta","ibeam==7","");
  meansp0X = his->GetMean();
  treeScan->Draw("sp2I>>hisDelta","ibeam==7","");
  meansp2I = his->GetMean();
  treeScan->Draw("sp2O>>hisDelta","ibeam==7","");
  meansp2O = his->GetMean();
  //
  //
  //
  Int_t indexes[10000];
  for (Int_t iside=0;iside<2;iside++){
    printf("*************************************\n");
    printf("Side%d\tRun\tchi2\t\tVgg\tVskirt\tVcover\n",iside);

    Int_t entries  =  treeScan->Draw(Form("sp0X/%f+sp2I/%f+sp2O/%f:run",meansp0X,meansp2I,meansp2O),Form("ibeam==7&&side==%d",iside));
    TMath::Sort(entries,treeScan->GetV1(),indexes,kFALSE);
    for (Int_t i=0; i<entries; i++){
      Int_t index = indexes[i];
      Int_t run   = treeScan->GetV2()[index];
      Float_t chi2=treeScan->GetV1()[index];
      Float_t vgg = 
	printf("%d\t%d\t%f\t%d\t%d\t%d\n",iside,run, chi2,GetVoltage(run,0), GetVoltage(run,1),GetVoltage(run,2));
    }
  }

   for (Int_t iside=0;iside<2;iside++){
    printf("***********************\n");
    printf("Side%d\tRun\tchi2\t\tVgg\tVskirt\tVcover\n",iside);

    Int_t entries  =  treeScan->Draw(Form("sp2I/%f+sp2O/%f:run",meansp0X,meansp2I,meansp2O),Form("ibeam==7&&side==%d",iside));
    TMath::Sort(entries,treeScan->GetV1(),indexes,kFALSE);
    for (Int_t i=0; i<entries; i++){
      Int_t index = indexes[i];
      Int_t run   = treeScan->GetV2()[index];
      Float_t chi2=treeScan->GetV1()[index];
      Float_t vgg = 
	printf("%d\t%d\t%f\t%d\t%d\t%d\n",iside,run, chi2,GetVoltage(run,0), GetVoltage(run,1),GetVoltage(run,2));
    }
  }
}




void MakeAliases(){
  /// use table

  chain->SetAlias("VggA","GetVoltage(run,0)");
  chain->SetAlias("VggC","GetVoltage(run,1)");
  chain->SetAlias("VcoA","GetVoltage(run,2)");
  chain->SetAlias("VcoC","GetVoltage(run,3)");
  chain->SetAlias("VskA","GetVoltage(run,4)");
  chain->SetAlias("VskC","GetVoltage(run,5)");
  //
  // cuts
  //
  chain->SetAlias("TisOK","mdEdx>5&&entries>400");
  chain->SetAlias("CisOK","nCl.fElements>entries*0.5&&eY.fElements<0.01");
  chain->SetAlias("ATisOK","RA.mdEdx>5&&RA.entries>400");
  chain->SetAlias("ACisOK","RA.nCl.fElements>RA.entries*0.5&&abs(RA.dY.fElements-dY.fElements)<0.3&&RA.eY.fElements<0.01");
  //
  // shortcuts
  //
  chain->SetAlias("lX","X.fElements");   //
  chain->SetAlias("tY","kY.fElements");   //
  chain->SetAlias("dE","sign(Y.fElements)*(tan(10*pi/180.)*X.fElements-abs(Y.fElements))");
  chain->SetAlias("RdY","(dY.fElements-RA.dY.fElements)");
  
  chain->SetAlias("dIF","lX-85.2");
  chain->SetAlias("dOF","245.8-lX");
  chain->SetAlias("dUtY","((GetVoltage(run,0)-GetVoltage(RA.run,0))/400.*tY)");
  //
  chain->SetAlias("cIF","dUtY/(1+dIF/13.0)");
  chain->SetAlias("cOF","dUtY/(1+dOF/22.0)");
  chain->SetAlias("cIF2","(dUtY*dIF/13.0)/((1.+dIF/13.0)^2)");
  chain->SetAlias("cOF2","(dUtY*dOF/22.0)/((1.+dOF/22.0)^2)");

  chain->SetAlias("cE","(GetVoltage(run,0)-GetVoltage(RA.run,0))/400.*sign(dE)/(1+(abs(dE)-1.55)/1.4)");
  chain->SetAlias("cE2","(GetVoltage(run,0)-GetVoltage(RA.run,0))/400.*sign(dE)*((abs(dE)-1.55)/1.4)/((1+(abs(dE)-1.55)/1.4)^2)");
}



void MakeAliasesBoth(){
  /// cuts - slect good tracks

  chain->SetAlias("TisOK","mdEdx>5&&entries>400");
  chain->SetAlias("ATisOK","(LTr.fSide==0)*(RA.mdEdx>5&&RA.entries>500)");
  chain->SetAlias("CTisOK","(LTr.fSide==1)*(RC.mdEdx>5&&RC.entries>500)");
  //
  chain->Draw(">>run","TisOK&&(ATisOK||CTisOK)","entryList");
  TEntryList *elist = (TEntryList*)gDirectory->Get("run");
  chain->SetEntryList(elist);
  //
  //
  chain->SetAlias("CisOK","(nCl.fElements>entries*0.5&&eY.fElements<0.01)");
  chain->SetAlias("ACisOK","RA.nCl.fElements>RA.entries*0.5&&abs(RA.dY.fElements-dY.fElements)<0.3&&RA.eY.fElements<0.01");
  chain->SetAlias("CCisOK","RC.nCl.fElements>RC.entries*0.5&&abs(RC.dY.fElements-dY.fElements)<0.3&&abs(RC.dZ.fElements-dZ.fElements)<0.3&&RC.eY.fElements<0.01");


  //
  // voltage table
  //
  chain->SetAlias("Vgg","((LTr.fSide==0)*GetVoltage(run,0)+(LTr.fSide==1)*GetVoltage(run,1))");
  chain->SetAlias("VcoA","((LTr.fSide==0)*GetVoltage(run,2)+(LTr.fSide==1)*GetVoltage(run,3))");
  chain->SetAlias("VskA","((LTr.fSide==0)*GetVoltage(run,4)+(LTr.fSide==1)*GetVoltage(run,5))");
  
  chain->SetAlias("dVgg","(((LTr.fSide==0)*(GetVoltage(run,0)-GetVoltage(RA.run,0))+((LTr.fSide==0)*(GetVoltage(run,0)-GetVoltage(RA.run,0))))/400.)");
  //
  // shortcuts
  //
  chain->SetAlias("RdY","((LTr.fSide==0)*(dY.fElements-RA.dY.fElements)+(LTr.fSide==1)*(dY.fElements-RC.dY.fElements))");

  chain->SetAlias("drift","(1.-abs(LTr.fP[1]+0.)/250)");
  chain->SetAlias("lX","X.fElements");   //
  chain->SetAlias("tY","kY.fElements");   //
  chain->SetAlias("dE","sign(Y.fElements)*(tan(10*pi/180.)*X.fElements-abs(Y.fElements))");
  
  chain->SetAlias("dIF","lX-85.2");
  chain->SetAlias("dOF","245.8-lX");
  chain->SetAlias("dUtY","dVgg*tY");
  //
  //
  //
  chain->SetAlias("cIF","dUtY/(1+dIF/13.0)");
  chain->SetAlias("cOF","dUtY/(1+dOF/22.0)");
  chain->SetAlias("cIF2","(dUtY*dIF/13.0)/((1.+dIF/13.0)^2)");
  chain->SetAlias("cOF2","(dUtY*dOF/22.0)/((1.+dOF/22.0)^2)");

  chain->SetAlias("cE","dVgg*sign(dE)/(1+(abs(dE)-1.5)/1.3)");
  chain->SetAlias("cE2","dVgg*sign(dE)*((abs(dE)-1.5)/1.3)/((1+(abs(dE)-1.5)/1.3)^2)");

}
  



void MakeFit(){

  Int_t  ntracks=3000000;
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParam;
  TMatrixD covMatrix;
  TString fstring="";
  fstring+="cIF++";
  fstring+="cIF2++";
  fstring+="cOF++";
  fstring+="cOF2++";
  fstring+="cE++";
  fstring+="cE2++";
  TCut cutA="LTr.fBeam>-1&&CisOK&&(ACisOK||CCisOK)";

  TString * strA = TStatToolkit::FitPlane(chain,"RdY", fstring.Data(),cutA, chi2,npoints,fitParam,covMatrix,-1.,0, ntracks);
  chain->SetAlias("fdY",strA->Data());
  
  printf("sqrt(Chi2/npoints)=%f\n",TMath::Sqrt(chi2/npoints));
  chain->Draw("RdY-fdY",cutA);
  

  TF1 fe("fe","[0]/(1+(x-[2])/[1])",2,10);
  fe.SetParameters(0.4,1,1.5);
  

}


   
//Examples:
// chain->Draw("RdY:(GetVoltage(run,0)-GetVoltage(RA.run,0))*tY/(1+dIF/10.)","TisOK&&ATisOK&&CisOK&&ACisOK&&LTr.fBundle>0&&LTr.fBeam!=3&&dIF<30","",10000) 




