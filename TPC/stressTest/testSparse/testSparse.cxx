/// \class testSparse
/// \brief Test to check THnSparse - Layout example for TPC calibration
/// 
/// Tests categories:
/// 1. CPU and memory consumption - using AliSysInfo class
/// 2. Filling, reading and merging
/// 3. Correctness of filling, and merging
/// 
/// 
/// Functions:
/// 
/// 1. TestSparse(niter,nsize) - 
///          Fill the THnSparse in niter chunks
///          for each chunks nsize entries filled
///          The THnSparse is saved after each chunk
///          Current momory and CPU information dumped to the text file
/// 2. testRead(niter)
///          Read THnSparses from disk
///          Current momory and CPU information dumped to the text file
/// 
/// 3. testMerge(niter)
///          Read THnSpares from disk
///          Merge histograms
///          Current momory and CPU information dumped to the text file            




#include "THnSparse.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "AliSysInfo.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TMath.h"

void testSparse(Int_t niter, Int_t nsize){
  ///

  Double_t xminTrack[9], xmaxTrack[9];
  Int_t    binsTrack[9];
  TString  axisName[9],axisTitle[9];
  //
  // 0 - delta   of interest
  // 1 - global  phi in sector number  as float
  // 2 - local   x
  // 3 - local   ky
  // 4 - local   kz
  //
  axisName[0]="delta";   axisTitle[0]="#Delta (cm)";
  binsTrack[0]=60;       xminTrack[0]=-0.6;        xmaxTrack[0]=0.6;
  //
  axisName[1]="sector";   axisTitle[1]="Sector Number";
  binsTrack[1]=180;       xminTrack[1]=0;        xmaxTrack[1]=18;
  //
  axisName[2]="localX";   axisTitle[2]="x (cm)";
  binsTrack[2]=53;       xminTrack[2]=85.;        xmaxTrack[2]=245.;
  //
  //
  axisName[3]="kZ";      axisTitle[3]="dz/dx";
  binsTrack[3]=36;       xminTrack[3]=-1.8;        xmaxTrack[3]=1.8;
  //
  THnSparse * fClusterDelta[2]={0,0};
  fClusterDelta[0] = new THnSparseS("testFull","testFull", 4, binsTrack,xminTrack, xmaxTrack);
  for (Int_t iter=0; iter<niter; iter++){
    gRandom->SetSeed(0);
    fClusterDelta[1] = new THnSparseS("testM","testM", 4, binsTrack,xminTrack, xmaxTrack);

    for (Int_t ipoint=0; ipoint<nsize; ipoint++){
      Double_t x[4]={0,0,0,0};
      x[0]=gRandom->BreitWigner(0,0.2);
      x[1]=gRandom->Rndm()*18;
      x[2]=85+gRandom->Rndm()*(245.-85.);
      x[3]=-1.8+gRandom->Rndm()*(3.6);
      fClusterDelta[0]->Fill(x);
      fClusterDelta[1]->Fill(x);
    }
    TFile f(Form("testSparse_%d.root",iter),"recreate");
    fClusterDelta[0]->Write();
    f.Close();
    TFile f2(Form("testSparse2_%d.root",iter),"recreate");
    fClusterDelta[1]->Write();
    f2.Close();
    Int_t bins0=fClusterDelta[0]->GetNbins()/1000;
    Int_t bins1=fClusterDelta[1]->GetNbins()/1000;
    Int_t n0=Int_t(fClusterDelta[0]->GetEntries()/1000);
    Int_t n1=Int_t(fClusterDelta[1]->GetEntries()/1000);
    printf("iter=%d\t%d\t%d\t%d\t%d\n",iter, bins0, bins1, n0,n1);
    AliSysInfo::AddStamp(Form("iter%d",iter), iter,bins0,bins1);
    delete fClusterDelta[1];
  }
  TTree * tree =AliSysInfo::MakeTree("syswatch.log");
  gSystem->Exec("cp syswatch.log syswatchFill.log");
  TFile f("syswatchFill.root","recreate");
  tree->Write("sparse");
  f.Close();
}

void testRead(Int_t nmax=100000){
  /// test read of THnSparse

  AliSysInfo::AddStamp("start", 0,0,0);
  for (Int_t i=0; i<nmax; i++){
    TFile f(Form("testSparse_%d.root",i));
    AliSysInfo::AddStamp(Form("open%d",i), i,0,0);
    THnSparse * his0 = (THnSparse*)f.Get("testFull");    
    if (!his0) break;
    AliSysInfo::AddStamp(Form("read%d",i), i,1,0);
    Int_t bins0=his0->GetNbins()/1000;
    Int_t n0=Int_t(his0->GetEntries()/1000);
    printf("iter=%d\t%d\t%d\n",i, bins0, n0);
    AliSysInfo::AddStamp(Form("nbins%d",i), i,2,bins0);
    AliSysInfo::AddStamp(Form("entries%d",i), i,3,n0);
    delete his0;
    AliSysInfo::AddStamp(Form("delete%d",i), i,4,0);
  }
  //
  TTree * tree =AliSysInfo::MakeTree("syswatch.log");
  gSystem->Exec("cp syswatch.log syswatchRead.log");
  TFile f("syswatchRead.root","recreate");
  tree->Write("sparse");
  f.Close();
}

void testMerge(Int_t nmax,Int_t nmerge=1){
  /// test read of THnSparse

  THnSparse * hisM=0;
  THnSparse * hisA[nmerge];
  AliSysInfo::AddStamp("start", 0,0,0);
  for (Int_t i=0; i<nmax; i++){
    TFile f(Form("testSparse2_%d.root",i));
    AliSysInfo::AddStamp(Form("open%d",i), i,0,0);
    THnSparse * his0 = (THnSparse*)f.Get("testM");    
    if (!his0) break;
    AliSysInfo::AddStamp(Form("read%d",i), i,1,0);
    if (hisM==0) {
      hisM=his0;
      for (Int_t im=0; im<nmerge; im++) hisA[im]=(THnSparse*) hisM->Clone();
      continue;
    }
    hisM->Add(his0);
    for (Int_t im=0; im<nmerge; im++) hisA[im]->Add(hisM);
    AliSysInfo::AddStamp(Form("merge%d",i), i,1,1);  // stamp for merging - id1=1
    Int_t bins0=hisM->GetNbins()/1000;
    Int_t n0=Int_t(hisM->GetEntries()/1000);
    printf("iter=%d\t%d\t%d\n",i, bins0, n0);
    AliSysInfo::AddStamp(Form("nbins%d",i), i,2,bins0); // stamp for nbins id1=2
    AliSysInfo::AddStamp(Form("entries%d",i), i,3,n0);  // stamp for entries id1=3
    delete his0; 
    AliSysInfo::AddStamp(Form("delete%d",i), i,4,0);    // stamp for delete id1=4
  }
  //
  TTree * tree =AliSysInfo::MakeTree("syswatch.log");
  gSystem->Exec("cp syswatch.log syswatchMerge.log");
  TFile f("syswatchMerge.root","recreate");
  tree->Write("sparse");
  f.Close();
}


void DrawDiff(){
  ///

  TTree * treeNew =AliSysInfo::MakeTree("testSparseNew/syswatch.log");
  TTree * treeOld =AliSysInfo::MakeTree("testSparse/syswatch.log");
  treeNew->SetMarkerStyle(25);
  treeOld->SetMarkerStyle(27);
  treeNew->SetMarkerColor(2);
  treeOld->SetMarkerColor(4);

  TGraph *grTime[2];
  TGraph *grVM[2];
  Int_t entries=0; 
  Double_t maxT=0, maxM=0;
  entries=treeOld->Draw("deltaT:id0","","");
  grTime[0] = new TGraph(entries, treeOld->GetV2(), treeOld->GetV1());
  maxT=TMath::Max(maxT,TMath::MaxElement(entries, treeOld->GetV1()));
  entries=treeNew->Draw("deltaT:id0","","");
  grTime[1] = new TGraph(entries, treeNew->GetV2(), treeNew->GetV1());
  maxT=TMath::Max(maxT,TMath::MaxElement(entries, treeNew->GetV1()));

  entries=treeOld->Draw("VM:id0","","");
  grVM[0] = new TGraph(entries, treeOld->GetV2(), treeOld->GetV1());
  maxM=TMath::Max(maxM,TMath::MaxElement(entries, treeOld->GetV1()));
  entries=treeNew->Draw("VM:id0","","");
  grVM[1] = new TGraph(entries, treeNew->GetV2(), treeNew->GetV1());
  maxM=TMath::Max(maxM,TMath::MaxElement(entries, treeNew->GetV1()));


  grTime[0]->SetMarkerStyle(25);
  grTime[1]->SetMarkerStyle(27);
  grTime[0]->SetMarkerColor(2);
  grTime[1]->SetMarkerColor(4);
  grTime[0]->Draw("alp");  grTime[1]->Draw("alp");
  grTime[1]->SetMaximum(maxT);
  grTime[1]->Draw("alp");
  grTime[0]->Draw("lp");



}


/*
  
  
  aliroot -b -q /u/miranov/AliRoot/trunk/TPC/stressTest/testSparse/testSparse.cxx++\(150,1000000\)

  counter=0 
  while [ $counter -lt 5 ] ; do
   let counter=$counter+1
   echo $counter
   mkdir test$counter
   cd test$counter
   bsub -q alice-t3_8h -m batch_dgrid2 --oo out$counter.log  aliroot -b -q /u/miranov/AliRoot/trunk/TPC/tmp/testSparse.cxx+\(100,2000000\)
   cd ../
  done; 

*/
