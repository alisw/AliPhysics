//
// Process ESD tracks  - 
// Extract TPC tracks  - write them to tree
//
/*
  .L AnalyzeESDtracks.C+
  .L AliGenInfo.C+
  AnalyzeESDtracks(567);   // process tracks
  // Tracks are written to the file "TPCtracks.root"
  // Now yo can analyze it
  TFile fesd("AliESDs.root");
  TTree * treeE = (TTree*)fesd.Get("esdTree");
  TFile f("TPCtracks.root")
  TTree * tree =(TTree*)f.Get("Tracks");
  AliComparisonDraw comp;
  comp->fTree = tree;

*/

// 

/*
.L AnalyzeESDtracks.C+
.L AliGenInfo.C+

//AnalyzeESDtracks(567);


TFile fesd("AliESDs.root");
TTree * treeE = (TTree*)fesd.Get("esdTree");
TFile f("TPCtracks.root")
TTree * tree =(TTree*)f.Get("Tracks");
AliComparisonDraw comp;
comp->fTree = tree;

TFile fs("TPCsignal.root");
TTree *treeB =(TTree*)fs.Get("SignalB");
TTree *treeN =(TTree*)fs.Get("SignalN");
TTree *treeS =(TTree*)fs.Get("Signal");
TTree *treef =(TTree*)fs.Get("Fit");


FitSignals(treeB,"Max-Median>100&&RMS06<2.5")
TFile ffit("FitSignal.root");
TTree * treeF = (TTree*)ffit->Get("Fit");


TChain chaincl("TreeR","TreeR")
chaincl.Add("TPC.RecPoints.root/Event0/TreeR")
chaincl.Add("TPC.RecPoints1.root/Event1/TreeR")
chaincl.Add("TPC.RecPoints2.root/Event2/TreeR")
chaincl.Add("TPC.RecPoints3.root/Event3/TreeR")

*/



#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TTreeStream.h"
#include "TEventList.h"
#include "TCut.h"
#include "TFitter.h"
#include  "TMatrixD.h"

#include  "AliLog.h"
#include  "AliMagF.h"

#include  "AliESD.h"
#include  "AliESDfriend.h"
#include  "AliESDtrack.h"
#include  "AliTracker.h"
#include  "AliTPCseed.h"
#include  "AliTPCclusterMI.h"
#include  "TTreeStream.h"
#include  "TF1.h"
#include  "TGraph.h"
#include  "AliSignalProcesor.h"
#include  "TCanvas.h"


void FitSignals(TTree * treeB, TCut cut="Max-Median>150&&RMS06<2&&abs(Median-Mean09)<0.5");


void AnalyzeESDtracks(Int_t run){
  //
  // output redirect 
  TTreeSRedirector *  pcstream = new TTreeSRedirector("TPCtracks.root");
  TTreeSRedirector &cstream = *pcstream;
  //
  // dummy magnetic field
  AliMagF mag("aaa","aaa",1,1,10);
  AliTracker::SetFieldMap(&mag,kTRUE);
  TFile f("AliESDs.root");
  TTree * tree =(TTree*)f.Get("esdTree"); 
  AliESD * esd =0;
  tree->SetBranchAddress("ESD",&esd);
  AliESDfriend *evf=0;
  tree->AddFriend("esdFriendTree","AliESDfriends.root");
  tree->SetBranchAddress("ESDfriend",&evf);

  Int_t nevents = tree->GetEntries();
  TClonesArray *clusters = new TClonesArray("AliTPCclusterMI",160);
  for (Int_t irow=0; irow<160; irow++){
    new ((*clusters)[irow])  AliTPCclusterMI;   // iitial dummy clusters
  }
  
  for (Int_t ievent=0; ievent<nevents; ievent++){
    tree->GetEntry(ievent);
    if (!esd) continue;
    if (!evf) continue;
    esd->SetESDfriend(evf); //Attach the friend to the ESD      
    for (Int_t itrack =0; itrack<esd->GetNumberOfTracks(); itrack++){
      // Int_t itrack = 0;      
      if (esd->GetTrack(itrack)->GetFriendTrack()==0) continue;
      AliESDtrack * etrack = esd->GetTrack(itrack);
      AliESDfriendTrack * ftrack = (AliESDfriendTrack *)esd->GetTrack(itrack)->GetFriendTrack();
      AliTPCseed * seed =  (AliTPCseed*)(ftrack->GetCalibObject(0));
      if (!seed) continue;  
      for (Int_t irow=0; irow<160; irow++){
	if (seed->GetClusterFast(irow)){
	  AliTPCclusterMI * cl = new ((*clusters)[irow])  AliTPCclusterMI(*(seed->GetClusterFast(irow)));
	  cl->SetLabel(itrack,0);
	}
	else{
	  AliTPCclusterMI * cl = (AliTPCclusterMI*)clusters->At(irow);
	  cl->SetX(0); cl->SetY(0); cl->SetZ(0); cl->SetQ(0); cl->SetLabel(-1,0);
	}
      }
      Float_t dEdx = seed->GetdEdx();
      Float_t dEdxI = seed->CookdEdx(0.05,0.6,0,77);
      Float_t dEdxO = seed->CookdEdx(0.05,0.6,78,155);
      Int_t ncl = seed->GetNumberOfClusters();
      cstream<<"Tracks"<<
	"Run="<<run<<
	"Ncl="<<ncl<<
	"Event="<<ievent<<
	"dEdx="<<dEdx<<
	"dEdxI="<<dEdxI<<
	"dEdxO="<<dEdxO<<
	"Track.="<<seed<<
	"Etrack.="<<etrack<<
	"Cl.="<<clusters<<
	"\n";
    }  
  }
  delete pcstream;
  //
  // Fit signal part
  //
  TFile fs("TPCsignal.root");
  TTree *treeB =(TTree*)fs.Get("SignalB");
  //  FitSignals(treeB,"Max-Median>150&&RMS06<2.5");
}


void FitSignals(TTree * treeB, TCut cut){
  AliSignalProcesor proc;
  TF1 * f1 = proc.GetAsymGauss();
  TTreeSRedirector cstream("FitSignal.root");
  TFile *f = cstream.GetFile();

  char lname[100];
  sprintf(lname,"Fit%s", cut.GetTitle());
  TEventList *list = new TEventList(lname,lname);
  sprintf(lname,">>Fit%s", cut.GetTitle());
  treeB->Draw(lname,cut);
  treeB->SetEventList(list);
  for (Int_t ievent=0; ievent<list->GetN(); ievent++){
    char ename[100];
    sprintf(ename,"Fit%d", ievent);
    Double_t nsample = treeB->Draw("Graph.fY-Mean09:Graph.fX","","",1,ievent);
    Double_t * signal  = treeB->GetV1();
    Double_t * time  = treeB->GetV2();
    Double_t maxpos =0;
    Double_t max = 0;
    for (Int_t ipos = 0; ipos<nsample; ipos++){
      if (signal[ipos]>max){
	max    = signal[ipos];
	maxpos = ipos;
      }
    }

    Int_t first = TMath::Max(maxpos-10,0.);
    Int_t last  = TMath::Min(maxpos+60, nsample);
    //
    f->cd();
    TH1F his(ename,ename,last-first,first,last);
    for (Int_t ipos=0; ipos<last-first; ipos++){
      his.SetBinContent(ipos+1,signal[ipos+first]);
    }
    treeB->Draw("Sector:Row:Pad","","",1,ievent);
    Double_t sector = treeB->GetV1()[0];
    Double_t row    = treeB->GetV2()[0];
    Double_t pad    = treeB->GetV3()[0];
    //    TGraph  graph(last-first,&time[first],&signal[first]);
    f1->SetParameters(0.75*max,maxpos,1.1,0.8,0.25,0.2);
    //    TH1F * his = (TH1F*)graph.GetHistogram();
    his.Fit(f1,"q");
    his.Write(ename);
    gPad->Clear();
    his.Draw();
    gPad->Update();
    Double_t params[6];
    for (Int_t ipar=0; ipar<6; ipar++) params[ipar] = f1->GetParameters()[ipar];
    Double_t chi2 = TFitter::GetFitter()->Chisquare(6,params);
    TMatrixD cov(6,6);
    cov.SetMatrixArray(TFitter::GetFitter()->GetCovarianceMatrix());
    //
    // tail cancellation
    //
    Double_t x0[1000];
    Double_t x1[1000];
    Double_t x2[1000];
    for (Int_t ipos=0; ipos<last-first; ipos++){
      x0[ipos] = signal[ipos+first];
    }
    proc.TailCancelationALTRO1(x0,x1,0.85*0.339,0.09,last-first);
    proc.TailCancelationALTRO1(x1,x2,0.85,0.789,last-first);
    //
    sprintf(ename,"Cancel1_%d", ievent);
    TH1F his1(ename,ename,last-first,first,last);
    for (Int_t ipos=0; ipos<last-first; ipos++){
      his1.SetBinContent(ipos+1,x1[ipos]);
    }
    his1.Write(ename);
    sprintf(ename,"Cancel2_%d", ievent);
    TH1F his2(ename,ename,last-first,first,last);
    for (Int_t ipos=0; ipos<last-first; ipos++){
      his2.SetBinContent(ipos+1,x1[ipos]);
    }
    f1->SetParameters(0.75*max,maxpos,1.1,0.8,0.25,0.2);
    his2.Fit(f1,"q");
    his2.Write(ename);
    Double_t params2[6];
    for (Int_t ipar=0; ipar<6; ipar++) params2[ipar] = f1->GetParameters()[ipar];
    Double_t chi22 = TFitter::GetFitter()->Chisquare(6,params2);    
    TMatrixD cov2(6,6);
    cov2.SetMatrixArray(TFitter::GetFitter()->GetCovarianceMatrix());

    TGraph gr0(last-first, &time[first],x0);
    TGraph gr1(last-first, &time[first],x1);
    TGraph gr2(last-first, &time[first],x2);
    //
    cstream<<"Fit"<<
      "Sector="<<sector<<
      "Row="<<row<<
      "Pad="<<pad<<
      "First="<<first<<
      "Max="<<max<<
      "MaxPos="<<maxpos<<
      "chi2="<<chi2<<
      "chi22="<<chi22<<
      "Cov="<<&cov<<
      "Cov2="<<&cov2<<
      "gr0.="<<&gr0<<
      "gr1.="<<&gr1<<
      "gr2.="<<&gr2<<
      "p0="<<params[0]<<
      "p1="<<params[1]<<
      "p2="<<params[2]<<
      "p3="<<params[3]<<
      "p4="<<params[4]<<
      "p5="<<params[5]<<
      "p02="<<params2[0]<<
      "p12="<<params2[1]<<
      "p22="<<params2[2]<<
      "p32="<<params2[3]<<
      "p42="<<params2[4]<<
      "p52="<<params2[5]<<
      "\n";
    //    delete his;
  }

}



TChain *MakeChainCL(Int_t first, Int_t last){
  TChain *chaincl = new TChain("TreeR","TreeR");
  //
  char fname[100];
  for (Int_t i=first;i<last; i++){
    if (i>0) sprintf(fname,"TPC.RecPoints%d.root/Event%d/TreeR",i,i);
    if (i==0) sprintf(fname,"TPC.RecPoints.root/Event%d/TreeR",i);
    chaincl->Add(fname);
  }
  return chaincl;
}

TTree* GetTree(Int_t ievent){
  char fname[100];
  char tname[100];
  if (ievent>0) sprintf(fname,"TPC.RecPoints%d.root",ievent);
  if (ievent==0) sprintf(fname,"TPC.RecPoints.root");
  sprintf(tname,"Event%d/TreeR",ievent);
  TFile * f  = new TFile(fname);
  TTree * tree = (TTree*)f->Get(tname);
  return tree;

}
