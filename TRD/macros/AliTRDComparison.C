/*

//load mag field
AliRunLoader *loader = AliRunLoader::Open("galice.root");
loader->LoadgAlice();
gAlice = loader->GetAliRun();
AliMagF * field = gAlice->Field();
AliTracker::SetFieldMap(gAlice->Field(),1);
TFile fg("geometry.root");
gGeoManager = (TGeoManager*)fg.Get("Geo");
if (!gGeoManager) TGeoManager::Import("geometry.root");
.L AliGenInfo.C+
.L AliESDComparisonMI.C+
.L AliTRDComparison.C+g
//.L AliESDComparisonPIC.C+
.x TDR_style.C

MakeCompTr();
MakeTree();

 
 
MakeComp();
MakeAlias();
TFile fff("TRDdebug.root");
TTree * tree =(TTree*)fff.Get("tracklet");
TTree * tree2 =(TTree*)fff.Get("Tracks");
TTree * treefind =(TTree*)fff.Get("Find");
TTree * treebug =(TTree*)fff.Get("Bug");
TTree * treetofs =(TTree*)fff.Get("tofseed");
AliComparisonDraw compfind;
compfind.fTree = treefind;
AliComparisonDraw comptr;
comptr.fTree = tree;

TFile fff2("TOFdebug.root");
TTree * treetof =(TTree*)fff2.Get("Info");
TTree * treetoftrd =(TTree*)fff2.Get("TRD");

AliComparisonDraw comp2;
comp2->fTree =tree;
AliComparisonDraw comptof;
comptof->fTree =treetof;


*/

//ROOT includes
#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TString.h"
#include "TBenchmark.h"
#include "TStopwatch.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TTimer.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TView.h"
#include "TGeometry.h"
#include "TPolyLine3D.h"
#include "TLegend.h"
#include "TGraph.h"
//ALIROOT includes
#include "AliRun.h"
#include "AliStack.h"
#include "AliSimDigits.h"
#include "AliTPCParam.h"
#include "AliTPC.h"
#include "AliTPCLoader.h"
#include "AliDetector.h"
#include "AliTrackReference.h"
#include "AliTPCParamSR.h"
#include "AliTracker.h"
#include "AliMagF.h"
#include "AliHelix.h"
#include "AliESDtrack.h"
#include "AliTRDtrack.h"
#include "AliITStrackMI.h"
#include "AliESDV0MI.h"
#include "AliESDkink.h"
//#include "AliTRDparameter.h"
#include "AliTRDgeometryFull.h"
#include "AliTRDcluster.h"
#include "AliTRDtracker.h"
#include "AliTOFtrack.h"
#include "../TOF/AliTOFtrackerMI.h"
// COMPARISON includes
#include "../STEER/AliGenInfo.h"
#include "../STEER/AliESDComparisonMI.h"



AliComparisonDraw comp;
AliComparisonDraw comptr;
//AliTRDparameter   par("pica","vyjebana");
AliTRDgeometryFull    geom;
AliRunLoader  *fLoader=0; //= AliRunLoader::Open(fnGalice);

TCut cnesd("cnesd","abs(RC.fLabel)!=MC.fLabel");
TCut cprim("cprim","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<0.01&&abs(MC.fVDist[2])<0.01");
TCut citsin("cprim","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<0.01&&abs(MC.fVDist[2])<5");
TCut clength("clength","abs(TRD0.fIntegratedLength-MC.fTOFReferences[0].fLength)<30");
TCut c1("c1","TRD1.fStopped==0");    // track propagated to refernce
TCut c0("c0","TRD0.fStopped==0");    // track propagated to reference
TCut c2("c2","MC.fNTOFRef>0");       //track in tof
TCut c3("c3","MC.fNTOFRef>0&&MC.fNTRDRef>5");       //track in tof and TRD
TH1F hpull("hpull","hpull",50,-10,10);
TH1F hpull6("hpull6","hpull6",50,-6,6);


class AliTRDInfoMC: public TObject{
 public:
  void Update(AliMCInfo *info);  // update TRD info
  void Reset();                  // clear TRD info
  Short_t fLayer[6];         //indicates layer is on
  Short_t fLayerA[6];        // indicates layer after is ON
  Short_t fNLayers;          // number of crossed layers
  Short_t fLast;             // last crossed layer
  AliTrackReference fRef[6];    // trd reference
  AliTrackReference fRefLocal[6];    // trd local reference
  ClassDef(AliTRDInfoMC,1)   // TRD MC information 
};


class AliTRDInfoCl: public TObject{
  public:
  AliTRDInfoCl();
  void Update();
  AliTRDcluster fCl[6][25];      // clusters
  Float_t       fCount[6][25];   // number of clusters of given label
  Int_t   fNClusters;             // number of clusters
  Float_t fMeanQ[6];             // mean Q
  Int_t   fMaxBin[6];            // max time bin
  Int_t   fNcl[6];               // number of clusters 
  Float_t fMeanQS[6];             // mean Q
  Int_t   fMaxBinS[6];            // max time bin
  ClassDef(AliTRDInfoCl,1)      // TRD MC information 
};

ClassImp(AliTRDInfoCl)



class AliTRDInfoSeed: public TObject{
  public:
  AliTRDInfoSeed(){fLabel=-1; fLabel1=-1;fLabel2=-1;}
  void              Update(AliTRDseed *seeds);
  AliTRDseed        fSeeds[6];  // seeds from seeding
  Float_t           fFakeRatio; //
  Int_t             fLabel;     // seed label
  Int_t             fLabel1;    // full track label
  Int_t             fLabel2;    // second full track label
  Int_t             fNLayers;   //   
  Int_t             fNUsed;     //
  ClassDef(AliTRDInfoSeed,1)    // TRD MC information 
};


ClassImp(AliTRDInfoMC)
ClassImp(AliTRDInfoCl)
ClassImp(AliTRDInfoSeed)

void AliTRDInfoSeed::Update(AliTRDseed *seeds){
  //
  //
  //
  fNLayers=0;
  for (Int_t i=0;i<6;i++){
    fSeeds[i] = seeds[i];
    if (seeds[i].fN2>10) fNLayers++;
  }
}

void AliTRDInfoMC::Reset(){
  //
  //
  fNLayers=0;
  fLast =0;
  for (Int_t ilayer = 0; ilayer<6; ilayer++){
    fLayer[ilayer]=0;
    fLayerA[ilayer]=0;
  }
}
void  AliTRDInfoMC::Update(AliMCInfo *info)
{
  // 
  //
  Reset();
  for (Int_t i=0;i<info->fNTRDRef;i++){
    AliTrackReference * ref = (AliTrackReference*)info->fTRDReferences->At(i);
    Float_t localX = ref->LocalX();
    for (Int_t jplane = 0; jplane<6; jplane++){
      if (TMath::Abs(localX - (300+jplane*14.))<6){
	fLayer[jplane]++;
	for (Int_t kplane=jplane;kplane>=0;kplane--){
	  fLayerA[kplane]++;
	  fRef[jplane]=*ref;
	  //
	  fRefLocal[jplane].SetTrack(ref->GetTrack());
	  fRefLocal[jplane].SetPosition(ref->LocalX(), ref->LocalY(),ref->Z());
	  Float_t alpha  = ref->Alpha();
	  Float_t dir[3] = {ref->Px(), ref->Py(), ref->Pz()};
	  Float_t dirR[3];
	  dirR[0] = dir[0]*TMath::Cos(-alpha) - dir[1]*TMath::Sin(-alpha);
	  dirR[1] = dir[0]*TMath::Sin(-alpha) + dir[1]*TMath::Cos(-alpha);
	  fRefLocal[jplane].SetMomentum(1, dirR[1]/dirR[0],dir[2]/dir[0]);
	}
      }
    }
  }
  for (Int_t jplane = 0; jplane<6; jplane++){
    if (fLayer[jplane]>0) fNLayers++;
    if (fLayer[jplane]>0) fLast = jplane;
  }
}

AliTRDInfoCl::AliTRDInfoCl(){
  //
  //
  fNClusters =0;
  for (Int_t ilayer=0;ilayer<6;ilayer++){
    fMeanQ[ilayer]=0;
    fMaxBin[ilayer]=0;
    fNcl[ilayer]=0;
    fMeanQS[ilayer]=0;
    fMaxBinS[ilayer]=0;
    for (Int_t itime=0;itime<25;itime++){
      fCount[ilayer][itime]=0;
    }
  }
}


void  AliTRDInfoCl::Update(){
  fNClusters = 0;
  for (Int_t iplane=0;iplane<6;iplane++){
    fMeanQ[iplane]=-1;
    fNcl[iplane]  = 0;
    fMaxBin[iplane]=-1;
    //
    fMeanQS[iplane]=-1;
    fMaxBinS[iplane]=-1;

    Float_t maxQ  =-10;
    Float_t maxQS =-10;
    Int_t sums =0;
    for (Int_t itime=0;itime<25;itime++){
      if (fCl[iplane][itime].GetQ()>2){
	fNcl[iplane]++;
	fMeanQ[iplane]+=fCl[iplane][itime].GetQ();
	if (fCl[iplane][itime].GetQ()>=maxQ){
	  fMaxBin[iplane]= itime;
	  maxQ = fCl[iplane][itime].GetQ();
	}
      }
      if (fCl[iplane][itime].GetSumS()>2){
	sums++;
	fMeanQS[iplane]+=fCl[iplane][itime].GetSumS();
	if (fCl[iplane][itime].GetSumS()>=maxQS){
	  fMaxBinS[iplane]= itime;
	  maxQS = fCl[iplane][itime].GetSumS();
	}	
      }
    }
    if (fNcl[iplane]>0) fMeanQ[iplane]/=Float_t(fNcl[iplane]);    
    if (sums>0) fMeanQS[iplane]/=Float_t(sums);    
    fNClusters+=fNcl[iplane];
  }
}




void  MaxBinNorm(const char * expr = "TRD0.fTimBinPlane"){
  //
  // max time bin probability distribution
  //
  char  var[100];
  sprintf(var,"%s>>hispi",expr);
  TH1F * hispi= new TH1F("hispi","hispi",20,1,22);
  comp.fTree->Draw(var,"abs(MC.fPdg)==211&&TRD0.fN>90&&MC.fParticle.P()>3");
  hispi->Scale(1./(hispi->GetEntries()));
  hispi->Draw();
  sprintf(var,"%s>>hisel",expr);
  TH1F *hisel = new TH1F("hisel","hisel",20,1,22);
  comp.fTree->Draw(var,"abs(MC.fPdg)==11&&TRD0.fN>90&&MC.fParticle.P()>3");
  hisel->Scale(1./(hisel->GetEntries()));
  hispi->SetMarkerStyle(25); hispi->SetLineStyle(2);
  hisel->SetMarkerStyle(26); hisel->SetLineStyle(1);
  hisel->SetXTitle("Maximal Time Bin");
  hisel->SetYTitle("Probability");
  hisel->Draw();
  hispi->Draw("same");
  TLegend *legend = new TLegend(0.55,0.12,0.85,0.35, "Maximal Time Bin position");
  legend->SetBorderSize(1);
  legend->AddEntry(hispi ,"Pions");
  legend->AddEntry(hisel ,"Electrons");
  legend->Draw();

}


void  dedxNorm(const char * expr = "TRD0.fdEdxPlane/20"){
  //
  // max time bin probability distribution
  //
  char  var[100];
  sprintf(var,"%s>>hispi",expr);
  TH1F * hispi= new TH1F("hispi","hispi",100,1,120);
  comp.fTree->Draw(var,"abs(MC.fPdg)==211&&TRD0.fN>90&&MC.fParticle.P()>3");
  hispi->Scale(1./(hispi->GetEntries()));
  hispi->Draw();
  sprintf(var,"%s>>hisel",expr);
  TH1F *hisel = new TH1F("hisel","hisel",100,1,120);
  comp.fTree->Draw(var,"abs(MC.fPdg)==11&&TRD0.fN>90&&MC.fParticle.P()>3");
  hisel->Scale(1./(hisel->GetEntries()));
  hispi->SetMarkerStyle(25); hispi->SetLineStyle(2);
  hisel->SetMarkerStyle(26); hisel->SetLineStyle(1);
  hispi->SetXTitle("Normalized dEdx");
  hispi->SetYTitle("Probability");
  hispi->Draw();
  hisel->Draw("same");
  TLegend *legend = new TLegend(0.55,0.12,0.85,0.35, "Normalized dEdx");
  legend->SetBorderSize(1);
  legend->AddEntry(hispi ,"Pions");
  legend->AddEntry(hisel ,"Electrons");
  legend->Draw();
}


void EfficiencyPt()
{
  //
  //
  //
  comp.Eff("TR.Pt()","InfoMC.fNLayers>3&&abs(MC.fPdg)!=11"+c2+cprim,"TRD0.fN/InfoMC.fNLayers>10&&TRD0.fLab>0",10,0.3,2);
  TH1F * his = (TH1F*)comp.fRes->Clone();
  comp.Eff("TR.Pt()","InfoMC.fNLayers>3&&abs(MC.fPdg)!=11"+c2+cprim,"TRD0.fN/InfoMC.fNLayers>10&&TRD0.fLab<0",10,0.3,2);
  TH1F * hisf = (TH1F*)comp.fRes->Clone();
  his->SetMarkerStyle(25); his->SetLineStyle(2);
  hisf->SetMarkerStyle(26); hisf->SetLineStyle(1);
  his->SetXTitle("P_{t}(GeV/c)");
  his->SetYTitle("Tracking efficiency (%)");
  his->Draw();
  hisf->Draw("same");
  TLegend *legend = new TLegend(0.55,0.16,0.85,0.4, "TRD efficciency");
  legend->SetBorderSize(1);
  legend->AddEntry(his ,"Efficiency");
  legend->AddEntry(hisf ,"Fake tracks");
  legend->Draw();

}

void EfficiencyDip()
{
  //
  //
  //
  comp.Eff("atan(TR.fPz/TR.Pt())","InfoMC.fNLayers>3&&abs(MC.fPdg)!=11&&TR.Pt()>0.4"+c2+cprim,"TRD0.fN/InfoMC.fNLayers>10&&TRD0.fLab>0",10,-0.9,0.9);
  TH1F * his = (TH1F*)comp.fRes->Clone();
  comp.Eff("atan(TR.fPz/TR.Pt())","InfoMC.fNLayers>3&&abs(MC.fPdg)!=11&&TR.Pt()>0.4"+c2+cprim,"TRD0.fN/InfoMC.fNLayers>10&&TRD0.fLab<0",10,-0.9,0.9);
  TH1F * hisf = (TH1F*)comp.fRes->Clone();
  his->SetMarkerStyle(25); his->SetLineStyle(2);
  hisf->SetMarkerStyle(26); hisf->SetLineStyle(1);
  his->SetXTitle("#lambda(rad)");
  his->SetYTitle("Tracking efficiency (%)");
  his->Draw();
  hisf->Draw("same");
  TLegend *legend = new TLegend(0.55,0.16,0.85,0.4, "TRD efficciency");
  legend->SetBorderSize(1);
  legend->AddEntry(his ,"Efficiency");
  legend->AddEntry(hisf ,"Fake tracks");
  legend->Draw();
}

void ClusterRatio()
{
  TH1F *hisr = new TH1F("hisr","Percentage of found clusters",80,0,120);
  comp.fTree->Draw("100*TRD0.fN/(20*InfoMC.fNLayers)>>hisr","InfoMC.fNLayers>3&&TRD0.fLab>0&&abs(MC.fPdg)!=11&&TR.Pt()>0.4"+cprim+c2);
  hisr->SetXTitle("(%)");
  hisr->Draw();

}

void DEestimate(){
  TH1F * hderel = new TH1F("hderel","hderel",100,-150,150);
  comp.fTree->Draw("100*derel>>hderel",""+c0+c2+cprim);
  hderel->Fit("gaus");
  hderel->SetXTitle("Energy loss correction resolution (%)");
  hderel->Draw();
}

void DEestimate2(){
  //
  TH1F * hderel0 = new TH1F("hderel0","hderel0",100,-0,30);
  TH1F * hderel1 = new TH1F("hderel1","hderel1",100,-0,30);
  TH1F * hderela = new TH1F("hderela","hderela",100,-0,30);
  comp.fTree->Draw("100*(TPCE-TOFE)/TPCE>>hderela","Pt<1"+c0+c2+cprim);
  comp.fTree->Draw("100*(TPCE-TOFE)/TPCE>>hderel0","Pt<1&&TRD0.fBudget[0]<10"+c0+c2+cprim);
  comp.fTree->Draw("100*(TPCE-TOFE)/TPCE>>hderel1","Pt<1&&TRD0.fBudget[0]>10"+c0+c2+cprim);
  hderela->SetXTitle("(E_{TPC}-E_{TOF})/E_{TPC}(%)");
  hderela->SetLineStyle(1);
  hderel0->SetLineStyle(2);
  hderel1->SetLineStyle(3);
  hderela->Draw("");
  hderel0->Draw("same");
  hderel1->Draw("same");
  TLegend *legend = new TLegend(0.55,0.12,0.85,0.35, "Relative Energy Loss Between TPC and TOF");
  legend->SetBorderSize(1);
  legend->AddEntry(hderela ,"All tracks");
  legend->AddEntry(hderela ,"Non crossing of boundary");
  legend->AddEntry(hderela ,"Crossing boundary");
  legend->Draw();
}


void EffBudget(){
  //
  //
  //
  TF1 f1("f1","100*exp(-[0]*x)",0,100);
  comp.Eff("TRD.fBudget[0]","abs(ThetaM)<1&&abs(MC.fPdg)==211"+cprim,c2,10,1,100);
  TH1F * heffb = (TH1F*)comp.fRes->Clone();
  heffb->SetXTitle("Material budget (g/cm^{2})");
  heffb->SetYTitle("Efficiency to reach TOF (%)");
  heffb->Fit(&f1);
  heffb->Draw();

  
}

void TOFbackround(){
  TH1F * his = new TH1F("his","his",100,-100,10000);
  TCut cut = "InfoTOF.fNCluster>1&&tofm0";
  comp.fTree->Draw("InfoTOF.fTime1-33*sqrt(InfoTOF.fCl1.fZ**2+InfoTOF.fCl1.fR**2)>>his","InfoTOF.fTime1<100000"+cut);
  his->SetXTitle("TOF Time (ns)");
}

void TOFsignal(){
  TCut cut = "InfoTOF.fNCluster>0&&abs(MC.fPdg)==2212&&P<2";
  //
  TH1F * his0 = new TH1F("his0","his0",100,-1000,1000);
  comp.fTree->Draw("InfoTOF.fTime0-InfoTOF.fTimes0[4]>>his0","tofm0"+cut);
  TH1F * hisf = new TH1F("hisf","hisf",100,-1000,1000);
  comp.fTree->Draw("InfoTOF.fTime0-InfoTOF.fTimes0[4]>>hisf","!(tofm0||tofm0s)"+cut);
  his0->Draw();
  hisf->Draw("same");
  his0->SetXTitle("TOF Time (ns)");
}




void EffPIDK(Float_t xmin=0.4, Float_t xmax=2){
  //
  //
  //Float_t xmin=0.4; Float_t xmax=2;
 //  TCut cut= "1";
//   // TCut cut = c2;
//   // TCut cut = "tofm"
//   AliTOFtrackInfo::SetPIDMatchCuts();
//   comp.Eff("P","TRD0.fIntegratedLength>360&&abs(MC.fPdg)==321"+cut,"InfoTOF.GetCombPID(-1,0.1,0)==3",6,xmin,xmax);
//   TH1F *his0 = (TH1F*)comp.fRes->Clone();
//   comp.Eff("P","TRD0.fIntegratedLength>360&&abs(MC.fPdg)==321"+cut,"InfoTOF.GetCombPID(-1,0.1,1)==3",6,xmin,xmax);
//   TH1F *his1 = (TH1F*)comp.fRes->Clone();
//   //AliTOFtrackInfo::SetPIDMatchCuts(-1,-1,-1,100);  // disable pid match cuts  
//   //comp.Eff("P","TRD0.fIntegratedLength>360&&abs(MC.fPdg)==321"+cut,"InfoTOF.GetCombPID(-1,0.0,0)==3",6,xmin,xmax);
//   //TH1F *hisno = (TH1F*)comp.fRes->Clone();
//   //AliTOFtrackInfo::SetPIDMatchCuts(); //enable again
//   //
//   comp.Eff("P","TRD0.fIntegratedLength>360&&abs(MC.fPdg)==321"+cut,"InfoTOF.GetCombPID(-1,0.1,2)==3",6,xmin,xmax);
//   TH1F *hisold = (TH1F*)comp.fRes->Clone();
//   comp.Eff("P","TRD0.fIntegratedLength>360&&abs(MC.fPdg)==321"+cut,"InfoTOF.GetCombPID(-1,0.0,-1)==3",6,xmin,xmax);
//   TH1F *histpc = (TH1F*)comp.fRes->Clone();
//   his0->SetMarkerStyle(21);
//   his1->SetMarkerStyle(24);
//   //  hisno->SetMarkerStyle(26);
//   hisold->SetMarkerStyle(25);
//   histpc->SetMarkerStyle(23);
//   //
//   his0->Draw();
//   his1->Draw("same");
//   //hisno->Draw("same");
//   hisold->Draw("same");
//   histpc->Draw("same");
//   TLegend *legend = new TLegend(0.55,0.12,0.85,0.35, "PID efficiency ");
//   legend->SetBorderSize(1);
//   legend->AddEntry(his0 ,"TPC+New TOF v0");
//   legend->AddEntry(his1 ,"TPC+New TOF v1");
//   //legend->AddEntry(hisno ,"TPC+New TOF v0 - no PID match cuts");
//   legend->AddEntry(hisold ,"TPC+Old TOF");
//   legend->AddEntry(histpc ,"TPC only");
//   legend->Draw();
  //
}

void EffPIDP(Float_t xmin=0.4, Float_t xmax=2){
  //
  //
  //Float_t xmin=0.4; Float_t xmax=4;
  TCut cut= "1";
  // TCut cut= c2;
  //TCut cut = "tofm0"
  //AliTOFtrackInfo::SetPIDMatchCuts(-1,-1,-1,1000);   // default cuts
//   AliTOFtrackInfo::SetPIDMatchCuts();   // default cuts
//   comp.Eff("P","TRD0.fIntegratedLength>360&&abs(MC.fPdg)==2212"+cut,"InfoTOF.GetCombPID(-1,0.1,0)==4",6,xmin,xmax);
//   TH1F *his0 = (TH1F*)comp.fRes->Clone();
//   comp.Eff("P","TRD0.fIntegratedLength>360&&abs(MC.fPdg)==2212"+cut,"InfoTOF.GetCombPID(-1,0.1,1)==4",6,xmin,xmax);
//   TH1F *his1 = (TH1F*)comp.fRes->Clone();
//   comp.Eff("P","TRD0.fIntegratedLength>360&&abs(MC.fPdg)==2212"+cut,"InfoTOF.GetCombPID(-1,0.1,2)==4",6,xmin,xmax);
//   TH1F *hisold = (TH1F*)comp.fRes->Clone();
//   comp.Eff("P","TRD0.fIntegratedLength>360&&abs(MC.fPdg)==2212"+cut,"InfoTOF.GetCombPID(-1,0.1,-1)==4",6,xmin,xmax);
//   TH1F *histpc = (TH1F*)comp.fRes->Clone();
//   his0->SetMarkerStyle(21);
//   his1->SetMarkerStyle(24);
//   hisold->SetMarkerStyle(25);
//   histpc->SetMarkerStyle(23);
//   his0->Draw();
//   his1->Draw("same");
//   hisold->Draw("same");
//   histpc->Draw("same");
//   TLegend *legend = new TLegend(0.55,0.12,0.85,0.35, "PID efficiency ");
//   legend->SetBorderSize(1);
//   legend->AddEntry(his0 ,"TPC+New TOF v0");
//   legend->AddEntry(his1 ,"TPC+New TOF v1");
//   legend->AddEntry(hisold ,"TPC+Old TOF");
//   legend->AddEntry(histpc ,"TPC only");
//   legend->Draw();

  //
}





void PtRes(){
  
  
}


TH1F * MakeCumul(TH1F *his,Bool_t norm = kTRUE)
{
  TH1F *hcumul = (TH1F*)his->Clone();
  char name[1000];
  sprintf(name,"N%f",gRandom->Rndm());
  hcumul->SetName(name);
  hcumul->SetTitle(name);
  //
  Float_t all = (Float_t)his->GetSum();
  if (!norm) all =1;
  for (Int_t i=0;i<=his->GetNbinsX();i++){
    hcumul->SetBinContent(i,his->Integral(0,i)/all);
  } 
  hcumul->SetFillColor(0);
  return hcumul;
}

TGraph* His01(TH1F *his0,TH1F *his1){
  Int_t nbins = his0->GetNbinsX();
  Double_t *x = new Double_t[nbins];
  Double_t *y = new Double_t[nbins];
  for (Int_t i=0;i<nbins;i++) {
    x[i] = his0->GetBinContent(i+1);
    y[i] = his1->GetBinContent(i+1);
  }
  return new TGraph(nbins,x,y);
}


TH1F * PullCumul(const char *var, TCut cut, Int_t div=200, Float_t max=20){
  TH1F * hpullb = new TH1F("hpullb","hpullb",div,0,max);
  char v2[1000];
  sprintf(v2,"%s>>hpullb",var);
  comp.fTree->Draw(v2, cut);
  TH1F * res = MakeCumul(hpullb);
  delete hpullb;
  return res;
}



void MakeComp(){
  //
  // Set tree and set aliases
  //
  TFile *f = new TFile("trdComp1.root");
  TTree * tree= (TTree*)f->Get("Comparison");
  comp.fTree =tree;  
}

void MakeAlias(){

  comp.fTree->SetAlias("Pt","sqrt(MC.fParticle.fPx**2+MC.fParticle.fPy**2)");
  comp.fTree->SetAlias("P","sqrt(MC.fParticle.fPx**2+MC.fParticle.fPy**2+MC.fParticle.fPz**2)");
  comp.fTree->SetAlias("Ptpc","sqrt(fTrackRef.fPx**2+fTrackRef.fPy**2+fTrackRef.fPz**2)");
  comp.fTree->SetAlias("VRadius","sqrt(MC.fParticle.fVx**2+MC.fParticle.fVy**2)");
  comp.fTree->SetAlias("DRadius","sqrt(MC.fTRdecay.fX**2+MC.fTRdecay.fY**2)");
  comp.fTree->SetAlias("ThetaM","MC.fParticle.fPz/Pt");
  comp.fTree->SetAlias("rel0","(MC.fTOFReferences[0].P()-MC.fTPCReferences[0].P())/MC.fTOFReferences[0].P()");  //relative en. loss
  comp.fTree->SetAlias("rel1","(TR.P()-MC.fTPCReferences[0].P())/MC.fTPCReferences[0].P()");  //relative en. loss - in the last TRD
  comp.fTree->SetAlias("pt_tof","MC.fTOFReferences[0].Pt()");
  comp.fTree->SetAlias("dpt_tof","(abs(1/RC.fRp[4])-MC.fTOFReferences[0].Pt())/MC.fTOFReferences[0].Pt()");


  comp.fTree->SetAlias("d1pt_tof","(abs(RC.fRp[4])-1/MC.fTOFReferences[0].Pt())");
  comp.fTree->SetAlias("d1ptp_tof","(abs(RC.fRp[4])-1/MC.fTOFReferences[0].Pt())/sqrt(RC.fRc[14])");

  comp.fTree->SetAlias("dpt_trd","(abs(1/RC.fRp[4])-TR.Pt())/TR.Pt()");

  comp.fTree->SetAlias("tofth","(MC.fTOFReferences[0].fPz/MC.fTOFReferences[0].Pt())");
  comp.fTree->SetAlias("dth_tof","TRD0.fT[3]-tofth");
  comp.fTree->SetAlias("pullth0","(TRD0.fT[3]-tofth)/sqrt(TRD0.fCtt)");
  comp.fTree->SetAlias("phi_tof","atan2(MC.fTOFReferences[0].fPy,MC.fTOFReferences[0].fPx)");
  comp.fTree->SetAlias("phi_trd","asin(TRD0.GetSnp())+TRD0.fAlpha");
  comp.fTree->SetAlias("pullphi0","(sin(asin(TRD0.GetSnp())+TRD0.fAlpha)-sin(atan2(MC.fTOFReferences[0].fPy,MC.fTOFReferences[0].fPx)))/sqrt(TRD0.fCee)");
  comp.fTree->SetAlias("dphir0","(asin(TRD0.GetSnp())+TRD0.fAlpha-atan2(MC.fTOFReferences[0].fPy,MC.fTOFReferences[0].fPx))");
  comp.fTree->SetAlias("dz0","(TRD0.fZ-MC.fTOFReferences[0].fZ)");
  comp.fTree->SetAlias("pully0","TRD0.fY/sqrt(TRD0.fCyy)");
  comp.fTree->SetAlias("pullz0","(TRD0.fZ-MC.fTOFReferences[0].fZ)/sqrt(TRD0.fCzz)");
  //
  comp.fTree->SetAlias("sigmay","sqrt(TRD0.fCyy)*(1+1/(1+abs(4./(RC.fRp[4]))))");
  comp.fTree->SetAlias("pully3","TRD0.fY/sigmay");  
  comp.fTree->SetAlias("deltaz0","(TRD0.fZ-MC.fTOFReferences[0].fZ)");
  //
  comp.fTree->SetAlias("pullc1","(abs(RC.fRp[4])-1/TR.Pt())/sqrt(RC.fRc[14])");
  comp.fTree->SetAlias("deltac1","(abs(RC.fRp[4])-1/TR.Pt())*TR.Pt()");
  comp.fTree->SetAlias("dist","sqrt(deltaz0**2+TRD0.fY**2)");
  comp.fTree->SetAlias("pull0","pully0**2+pullz0**2");
  comp.fTree->SetAlias("pullc0","(abs(TRD0.fC)-(1/(TRD0.fLocalConvConst*MC.fTOFReferences[0].Pt())))/sqrt(TRD0.fCcc)");
  //
  comp.fTree->SetAlias("prob5","exp(-(TRD1.fTracklets[5].fChi2)/16.)*max((TRD1.fTracklets[5].fNFound-6)/9.,0)");
  comp.fTree->SetAlias("prob4","exp(-(TRD1.fTracklets[4].fChi2)/16.)*max((TRD1.fTracklets[4].fNFound-6)/9.,0)");
  //
  comp.fTree->SetAlias("TPCE","sqrt(MC.fTPCReferences[2].P()**2+MC.fMass**2)");
  comp.fTree->SetAlias("TOFE","sqrt(MC.fTOFReferences[0].P()**2+MC.fMass**2)");
  comp.fTree->SetAlias("TRDE","sqrt(TR.P()**2+MC.fMass**2)");
  comp.fTree->SetAlias("dtpi","TRD0.fIntegratedTime[2]-10^12*MC.fTOFReferences.fTime");
  comp.fTree->SetAlias("dtk","TRD0.fIntegratedTime[3]-10^12*MC.fTOFReferences.fTime");
  comp.fTree->SetAlias("dtp","TRD0.fIntegratedTime[4]-10^12*MC.fTOFReferences.fTime");
  comp.fTree->SetAlias("dt","(abs(MC.fPdg)==2212)*dtp+(abs(MC.fPdg)==211)*dtpi+(abs(MC.fPdg)==321)*dtk"); //delta time
  comp.fTree->SetAlias("Beta","sqrt(TR.P()**2/(TR.P()**2+MC.fMass**2))");
  comp.fTree->SetAlias("dtsigma","10+(1-Beta)*280");
  comp.fTree->SetAlias("dptrel","(abs(TRD0.GetSignedPt())-MC.fTOFReferences[0].Pt())/MC.fTOFReferences[0].Pt()");
  //
  comp.fTree->SetAlias("dphi","atan2(TR.fY,TR..fX)-atan2(TR.fPy,TR..fPx)");
  comp.fTree->SetAlias("dphi0","acos((TR.fX*TR.fPx+TR.fY*TR.fPy)/(TR.Pt()*TR.R()))");
  comp.fTree->SetAlias("derel","(TRD0.fDE-(TPCE-TOFE))/(TPCE-TOFE)");
  //
  // TOF stuff
  //
  comp.fTree->SetAlias("tofm0","InfoTOF.fCl0.fLab[0]==abs(MC.fLabel)");
  comp.fTree->SetAlias("tofm1","InfoTOF.fCl1.fLab[0]==abs(MC.fLabel)");
  comp.fTree->SetAlias("tofm0s","InfoTOF.fCl0.fLab[0]>=abs(MC.fParticle.fDaughter[0]) && InfoTOF.fCl0.fLab[0]<=abs(MC.fParticle.fDaughter[1])");
  comp.fTree->SetAlias("tofm1s","InfoTOF.fCl1.fLab[0]>=abs(MC.fParticle.fDaughter[0]) && InfoTOF.fCl1.fLab[0]<=abs(MC.fParticle.fDaughter[1])");
  comp.fTree->SetAlias("tofm","tofm0||tofm1||tofm0s||tofm1s");

  comp.fTree->SetAlias("Yref3","InfoMC.fSeeds[3].fYref[0]+InfoMC.fSeeds[3].fYref[1]*(InfoMC.fRef[3].LocalX()-InfoMC.fSeeds[3].fX0)");
  comp.fTree->SetAlias("dyref3","Yref3-InfoMC.fRef[3].LocalY()");

  comp.fTree->SetAlias("Yref5","InfoMC.fSeeds[5].fYref[0]+InfoMC.fSeeds[5].fYref[1]*(InfoMC.fRef[5].LocalX()-InfoMC.fSeeds[5].fX0)");
  comp.fTree->SetAlias("dyref5","Yref5-InfoMC.fRef[5].LocalY()");

}


void MakeCompTr(){
  //
  // Set tree and set aliases
  //
  TFile *ff = new TFile("TRDdebug.root");
  TTree * treetr= (TTree*)ff->Get("tracklet");
  comptr.fTree =treetr;  
  
}


void DrawDE(TCut cut){
  comp.DrawXY("(TRD0.fDE)","TPCE-TOFE",c2,cut,5,0.01,0.1,0,0.2);
  comp.fRes->Draw();
}




void DrawYResol1(TCut cut){
  comp.fTree->Draw("TRD1.fY","MC.fNTOFRef>0&&MC.fNTRDRef>5"+c1+cut,"");
}

void DrawYResol0(TCut cut){
  comp.fTree->Draw("TRD0.fY","MC.fNTOFRef>0&&MC.fNTRDRef>5"+c0+cut,"");
}


void DrawTOFResY(Int_t ndiv =10,Float_t ptmin=0.5, Float_t ptmax =1.5, Float_t dymin=-0.5, Float_t dymax=0.5){
  //
  //
  //Float_t ptmin=0.5, ptmax =1.5, dymin=-0.5, dymax=0.5;
  //Float_t ptmin=0.5, ptmax =5, dymin=-0.3, dymax=0.3;
  //  Int_t   ndiv =5;
  TCut cl("cl","InfoCl.fNClusters>60");
  comp.DrawXY("TR.Pt()","TRD0.fY",c0+c2+cl,"1",ndiv,ptmin, ptmax,dymin, dymax,60);
  TH1F * histof = (TH1F*)comp.fRes->Clone();
  comp.DrawXY("TR.Pt()","TRD1.fY",c0+c2+cl,"1",ndiv,ptmin, ptmax,dymin, dymax,60);
  TH1F * histrd = (TH1F*)comp.fRes->Clone();
  histof->SetMarkerStyle(24);
  histrd->SetMarkerStyle(25);
  histof->SetXTitle("P_{t}(GeV)");
  histof->SetYTitle("#Delta_{r-#phi}(cm)");
  histof->Draw();
  histrd->Draw("same"); 
  TLegend *legend = new TLegend(0.55,0.12,0.85,0.35, "TRD resolution r-#phi (cm)");
  legend->SetBorderSize(1);
  legend->AddEntry(histof ,"Resolution in TOF plane");
  legend->AddEntry(histrd ,"Resolution in last TRD plane");
  legend->Draw();
}

void DrawTOFResZ(Float_t ptmin=0.5, Float_t ptmax =1.5, Float_t dzmin=-2.5, Float_t dzmax=2.5){
  //
  comp.DrawXY("TR.Pt()","TRD0.fZ-MC.fTOFReferences.fZ",c0+c2,"1",10,ptmin,ptmax,dzmin,dzmax,30);
  TH1F * histof = (TH1F*)comp.fRes->Clone();
  comp.DrawXY("TR.Pt()","TRD1.fZ-TR.fZ",c0+c2,"1",10,ptmin,ptmax,dzmin,dzmax,30);
  TH1F * histrd = (TH1F*)comp.fRes->Clone();
  histof->SetMarkerStyle(24);
  histrd->SetMarkerStyle(25);
  histof->SetXTitle("P_{t}(GeV)");
  histof->SetYTitle("#Delta_{z}(cm)");
  histof->Draw();
  histrd->Draw("same"); 
  TLegend *legend = new TLegend(0.55,0.12,0.85,0.35, "TRD resolution z (cm)");
  legend->SetBorderSize(1);
  legend->AddEntry(histof ,"Resolution in TOF plane");
  legend->AddEntry(histrd ,"Resolution in last TRD plane");
  legend->Draw();
}


void MakeCumul(TCut & cut)
{
  //  TCut cut ="1&&TR.Pt()>0.4";
  TH1F * hisy = PullCumul("abs(pully0)",cut+cprim+c0+c2+c3,1000,20);
  TH1F * hisz = PullCumul("abs(pullz0)",cut+cprim+c0+c2+c3,1000,20);
  TH1F * hisyz = PullCumul("sqrt(pullz0**2+pully0**2)",cut+cprim+c0+c2+c3,1000,20);
  hisz->SetXTitle("Pull (unit)");
  hisz->SetYTitle("Cumulative function");
  hisy->SetLineStyle(1);hisz->SetLineStyle(2);hisyz->SetLineStyle(3);
  hisz->Draw(); hisy->Draw("same");
  //hisyz->Draw("same");  
  TLegend *legend = new TLegend(0.55,0.12,0.85,0.35, "TRD pull (unit)");
  legend->SetBorderSize(1);
  legend->AddEntry(hisy ,"Pull y");
  legend->AddEntry(hisz,"Pull z");
  //  legend->AddEntry(hisyz,"Pull yz");
  legend->Draw();
}

void MakeCumulAbs(TCut & cut)
{
  // TCut cut = "1";
  TH1F * hisy = PullCumul("abs(TRD0.fY)",cut+cprim+c0+c2+c3,1000,6);
  TH1F * hisz = PullCumul("abs(deltaz0)",cut+cprim+c0+c2+c3,1000,6);  
  TH1F * hisyz = PullCumul("sqrt(TRD0.fY**2+deltaz0**2)",cut+cprim+c0+c2+c3,1000,6);
  hisz->SetXTitle("Delta (cm)");
  hisz->SetYTitle("Cumulative function");
  //
  hisy->SetLineStyle(1);hisz->SetLineStyle(2);hisyz->SetLineStyle(3);
  hisz->Draw(); hisy->Draw("same");hisyz->Draw("same");  
  TLegend *legend = new TLegend(0.55,0.12,0.85,0.35, "TRD residulas (cm)");
  legend->SetBorderSize(1);
  legend->AddEntry(hisy ,"Delta y (cm)");
  legend->AddEntry(hisz,"Delta z (cm)");
  legend->AddEntry(hisyz,"Delta  (cm)");
  legend->Draw();
}

void LoadClusters(TObjArray *arr, Int_t event)
{
  //
  //
  // event=0
  printf("Load cluster - event%d\n",event);
  fLoader->SetEventNumber(event);
  AliLoader * trdl = (AliLoader*)fLoader->GetLoader("TRDLoader");
  trdl->LoadRecPoints();

  TTree * ClusterTree = trdl->TreeR();
  if (!ClusterTree) return;
  TBranch * branch = ClusterTree->GetBranch("TRDcluster"); 
  TObjArray *clusterArray = new TObjArray(1000); 
  branch->SetAddress(&clusterArray);
  Int_t nEntries = (Int_t) ClusterTree->GetEntries();
  //
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {    
    ClusterTree->GetEvent(iEntry);  
    for (Int_t icl =0; icl<clusterArray->GetEntriesFast();icl++){
      AliTRDcluster *cl = (AliTRDcluster*)clusterArray->UncheckedAt(icl);
      if (!cl) continue;
      Int_t detector=cl->GetDetector();
      Int_t localTimeBin=cl->GetLocalTimeBin();
      //      Int_t sector=geom.GetSector(detector);
      Int_t plane=geom.GetPlane(detector);
      for (Int_t ilab=0; ilab<3; ilab++){
	Int_t label = cl->GetLabel(ilab);
	Int_t pos = label;      
	if (label<0) continue;
	if (localTimeBin>21||localTimeBin<0) continue;
	if (plane<0||plane>5) continue;
	
	AliTRDInfoCl * info =( AliTRDInfoCl*)arr->At(pos); 
	if (!info) {
	  info = new AliTRDInfoCl;
	  arr->AddAt(info,pos);
	} 
	//
	if (info->fCount[plane][localTimeBin]>0){
	  if (TMath::Abs(info->fCl[plane][localTimeBin].GetY()-cl->GetY())<20 &&
	      TMath::Abs(info->fCl[plane][localTimeBin].GetZ()-cl->GetZ())<20){
	    info->fCount[plane][localTimeBin]++;
	  }
	}else{
	  info->fCl[plane][localTimeBin]=*cl;
	  info->fCount[plane][localTimeBin]++;
	}
      }
    }
  }  
  for (Int_t i=0;i<arr->GetEntriesFast();i++){
    AliTRDInfoCl * info =( AliTRDInfoCl*)arr->At(i); 
    if (info) info->Update();
  }
}



void ReadSeeds(Int_t eventNr, TObjArray *sarray)
{
  //
  // read seeds form the file
  //
  AliTRDseed seeds[6];
  AliTRDseed *pseeds[6]={&seeds[0],&seeds[1],&seeds[2],&seeds[3],&seeds[4],&seeds[5]};
  Int_t label, label1, label2, seventNr, nused;
  Float_t fakeRatio;
  TFile f("TRDdebug.root");
  TTree  *tree =(TTree*)f.Get("Seeds2");
  if (!tree) return;
  TBranch * br0 = tree->GetBranch("S0.");
  TBranch * br1 = tree->GetBranch("S1.");
  TBranch * br2 = tree->GetBranch("S2.");
  TBranch * br3 = tree->GetBranch("S3.");
  TBranch * br4 = tree->GetBranch("S4.");
  TBranch * br5 = tree->GetBranch("S5.");
  TBranch * brL = tree->GetBranch("Label");
  TBranch * brL1 = tree->GetBranch("Label1");
  TBranch * brL2 = tree->GetBranch("Label2");
  TBranch * brE = tree->GetBranch("EventNr");
  TBranch * brU = tree->GetBranch("NUsed");
  TBranch * brFR = tree->GetBranch("FakeRatio");
  brE->SetAddress(&seventNr);
  brL->SetAddress(&label);
  brL1->SetAddress(&label1);
  brL2->SetAddress(&label2);
  brU->SetAddress(&nused);
  brFR->SetAddress(&fakeRatio);
  br0->SetAddress(&pseeds[0]);
  br1->SetAddress(&pseeds[1]);
  br2->SetAddress(&pseeds[2]);
  br3->SetAddress(&pseeds[3]);
  br4->SetAddress(&pseeds[4]);
  br5->SetAddress(&pseeds[5]);
  //
  // delete all seeds
  for (Int_t ip=0;ip<sarray->GetEntriesFast();ip++){
    AliTRDseed * seed = (AliTRDseed*)sarray->UncheckedAt(ip);
    if (seed){
      delete seed;
    }
  }
  sarray->Clear();
  for (Int_t ip=0;ip<tree->GetEntries();ip++){
    tree->GetEntry(ip);
    if (seventNr!=eventNr) continue;    
    //if (nused>15) continue;
    if (!sarray->At(label)){
      AliTRDInfoSeed * seed =  new AliTRDInfoSeed;
      seed->Update(seeds);
      seed->fLabel       = label;
      seed->fLabel1      = label1;
      seed->fLabel2      = label2;
      seed->fFakeRatio = fakeRatio; 
      seed->fNUsed = nused; 
      sarray->AddAt(seed,label);
    }
  } 
}

void ReadTOFs(Int_t eventNr, TObjArray *tofs)
{
  //   
 //  tofs->Clear();
//   AliTOFtrackInfo *infotof = new AliTOFtrackInfo;
//   TFile ftof("TOFdebug.root");
//   TTree *   treetof =(TTree*)ftof.Get("Info");
//   if (!treetof) return;
//   TBranch * br30 = treetof->GetBranch("Info.");
//   br30->SetAddress(&infotof);
//   //
//   // read tof info
//   for (Int_t ip=0;ip<treetof->GetEntries();ip++){
//     treetof->GetEntry(ip);
//     if (eventNr!=infotof->fEventNr) continue;
//     //	if (eventNr>info->fEventNr) continue;
//     //
//     Int_t label = TMath::Abs(infotof->fLab);
//     if (label==0) continue;
//     //    if (label>=kmaxlabel) continue;
//     if (!(tofs->At(label))){
//       tofs->AddAt(new AliTOFtrackInfo(*infotof),label);	  
//     }
//     else{
//       if (infotof->fLength0>100)
// 	tofs->AddAt(new AliTOFtrackInfo(*infotof),label);
//     }
//   }
}

void ReadTracks(Int_t eventNr, TObjArray *esds, TObjArray *trds)
{
  //
  //  read debug tracks
  //
  AliESDtrack * esdp     = new AliESDtrack;
  AliTRDtrack * trdp     = new AliTRDtrack;
  Int_t teventNr;
  TFile f2("TRDdebug.root");
  TTree * tree2 =(TTree*)f2.Get("Tracks");
  if (!tree2) return;
  TBranch * br20 = tree2->GetBranch("ESD.");
  TBranch * br21 = tree2->GetBranch("trd.");
  TBranch * br22 = tree2->GetBranch("EventNr");
  br20->SetAddress(&esdp);
  br21->SetAddress(&trdp);  
  br22->SetAddress(&teventNr);
  //
  // delete all event
  for (Int_t ip=0;ip<esds->GetEntriesFast();ip++){
    if (esds->At(ip)){
      delete esds->At(ip);
      delete trds->At(ip);
    }
  }
  esds->Clear();
  trds->Clear(); 
  for (Int_t ip=0;ip<tree2->GetEntries();ip++){
    tree2->GetEntry(ip);
    if (eventNr!=teventNr) continue;
    //
    Int_t label = TMath::Abs(esdp->GetLabel());
    if (label==0) continue;
    if (!(esds->At(label))){
      esds->AddAt(new AliESDtrack(*esdp),label);
      trds->AddAt(new AliTRDtrack(*trdp),label);
    } else{
      AliTRDtrack * track =(AliTRDtrack*)trds->At(label);
      if (track->GetNumberOfClusters() <trdp->GetNumberOfClusters()){
	esds->AddAt(new AliESDtrack(*esdp),label);
	trds->AddAt(new AliTRDtrack(*trdp),label);
      }
    }
  }  
}



void MakeTree()
{
  fLoader = AliRunLoader::Open("galice.root");
  AliESDtrack * esd     = new AliESDtrack;
  AliTRDtrack * ptrd     = new AliTRDtrack;
  AliTRDtrack * ptrd0    = new AliTRDtrack;
  AliTRDtrack * ptrd1    = new AliTRDtrack;
  AliTRDtrack * ptrdb    = new AliTRDtrack;
  //
  AliTRDtrack * trd     = ptrd;
  AliTRDtrack * trd0    = ptrd0;
  AliTRDtrack * trd1    = ptrd1;
  AliTRDtrack * trdb    = ptrdb;
  //
  AliMCInfo* info       = new AliMCInfo;
  AliTRDInfoMC *trdinfo = new AliTRDInfoMC;
  AliTRDInfoCl *infocl  = new AliTRDInfoCl;
  AliTRDInfoSeed *infoSeed = new AliTRDInfoSeed;
  //  AliTOFtrackInfo *infotof = new AliTOFtrackInfo;
  // AliTOFtrackInfo dummytof;
  //  dummytof.fLab=-1;
  AliTRDInfoSeed  dummySeed;
  AliTRDtrack dummytrd;
  AliTRDtrack dummytrd0;
  AliTRDtrack dummytrd1;
  AliTRDtrack dummytrdb;
  AliESDtrack dummyesd;
  dummySeed.fLabel =-1;
  AliTrackReference *ref = new AliTrackReference;
  TObjArray         * clarray  = new TObjArray(6*1000000);
  const Int_t kmaxlabel=1500000;
  //
  // MC part
  //
  TFile f("cmpESDTracks.root");
  TTree * tree1 = (TTree*) f.Get("ESDcmpTracks");
  TBranch * brmc  = tree1->GetBranch("MC");
  brmc->SetAddress(&info);
  //
  // 
  TObjArray *esds  = new TObjArray(kmaxlabel);
  TObjArray *trds  = new TObjArray(kmaxlabel);
  TObjArray *tofs  = new TObjArray(kmaxlabel);  
  TObjArray *seeds = new TObjArray(kmaxlabel);
  //
  //
  //
  //
  // comparison part
  //
  TFile f3("trdComp1.root","new");
  TTree   * tree3   = new TTree("Comparison","Comparison");
  tree3->Branch("MC.","AliMCInfo",&info);
  tree3->Branch("RC.","AliESDtrack",&esd);
  tree3->Branch("TRD.","AliTRDtrack",&trd);
  tree3->Branch("TRD0.","AliTRDtrack",&trd0);  // track at tof ref
  tree3->Branch("TRD1.","AliTRDtrack",&trd1);  // track at last TRD ref
  tree3->Branch("TRDb.","AliTRDtrack",&trdb);  // track at tof ref
  tree3->Branch("TR.","AliTrackReference",&ref);  // track at last TRD ref
  tree3->Branch("InfoMC.","AliTRDInfoMC",&trdinfo);  // info about TRD track MC
  tree3->Branch("InfoCl.","AliTRDInfoCl",&infocl);  // info about TRD track MC
  tree3->Branch("InfoSeed.","AliTRDInfoSeed",&infoSeed);  // info about TRD track MC
  //  tree3->Branch("InfoTOF.","AliTOFtrackInfo",&infotof);  // info about TOF
  //
  //
  Int_t lastevent=-1;

  for (Int_t i=0;i<tree1->GetEntries();i++){
    brmc->GetEntry(i);
    //    if (info->fEventNr>5) break;
    if (lastevent!=info->fEventNr){
      printf("Read Event\t%d\n", info->fEventNr);
      ReadSeeds(info->fEventNr, seeds);
      ReadTOFs(info->fEventNr, tofs);
      ReadTracks(info->fEventNr,esds,trds);
      lastevent = info->fEventNr;
      clarray->Clear();
      LoadClusters(clarray,lastevent);
      //
    }
    //    if (info->fNTOFRef==0 && info->fNTRDRef<5) continue;
    Int_t label = info->fLabel;    
    trdinfo->Update(info);
    if (clarray->At(label)){
      infocl = (AliTRDInfoCl*)clarray->At(label);
      infocl->Update();
    }
    if (seeds->At(label)){
      infoSeed = (AliTRDInfoSeed*)seeds->At(label);
    }else{
      infoSeed = &dummySeed;
    }
    //    if (info->fNTPCRef<3) continue;
    trd0->SetStop(kTRUE);
    trdb->SetStop(kTRUE);
    trd1->SetStop(kTRUE);
 //    infotof = (AliTOFtrackInfo*)(tofs->At(label));
//     if (!infotof) {
//       infotof = &dummytof;
//       infotof->fNCluster=0; 
//       esd   = (AliESDtrack*)(esds->At(label));
//       if (esd){
// 	Double_t pid[5];
// 	esd->GetTPCpid(pid);
// 	Double_t suml=0;
// 	for (Int_t ip=0;ip<5;ip++) suml+=pid[ip];
// 	for (Int_t ip=0;ip<5;ip++) {
// 	  if (suml>0.0000000001) infotof->fRefPID[ip]= pid[ip]/suml;
// 	  else infotof->fRefPID[ip]=0.2;
// 	}
//       }
//    }
    trd   = &dummytrd;
    trd0  = &dummytrd0;
    trd1  = &dummytrd1;
    trdb  = &dummytrdb;
    esd   = &dummyesd;
    //
    if (esds->At(label)) {
      trd0  = ptrd0; 
      trd1  = ptrd1; 
      trdb  = ptrdb; 
      //      if (i%1000==0) printf("\n%d\n",i);
      trd   = (AliTRDtrack*)(trds->At(label));
      esd   = (AliESDtrack*)(esds->At(label));
      //
      *trd0 = *trd;
      if (trd->GetBackupTrack()) 
	*trdb = *trd->GetBackupTrack();
      else 
	 *trdb = *trd;
      *trd1 = *trd;
      if (info->fNTOFRef>0){
	AliTrackReference * ref = ((AliTrackReference*)(info->fTOFReferences->At(0)));
	Double_t x   = ((AliTrackReference*)(info->fTOFReferences->At(0)))->X();
	Double_t y   = ((AliTrackReference*)(info->fTOFReferences->At(0)))->Y();
	Double_t phi    = TMath::ATan2(y,x);	
	Double_t radius = TMath::Sqrt(x*x+y*y);
	//
	// rotate and propagate to local system
	//
	trd0->SetStop(kFALSE);
	if (!trd0->Rotate(phi-trd0->GetAlpha())) trd0->SetStop(kTRUE);
	else{
	  Double_t xyz0[3], param[7];
	  trd0->GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);
	  Double_t xyz1[3]={ref->X(),ref->Y(),ref->Z()};
	  AliTracker::MeanMaterialBudget(xyz0,xyz1,param);
	  if (!trd0->PropagateTo(radius,param[1],param[0]*param[4])) 
             trd0->SetStop(kTRUE);
	}
	trdb->SetStop(kFALSE);
	if (!trdb->Rotate(phi-trdb->GetAlpha())) trdb->SetStop(kTRUE);
	else
	  if (!trdb->PropagateTo(radius)) trdb->SetStop(kTRUE);
      }    

      if (info->fNTRDRef>0){
	Double_t cradius=0;
	Double_t cphi=0;
	for (Int_t itrd =0; itrd<info->fNTRDRef;itrd++){
	  Double_t x   = ((AliTrackReference*)(info->fTRDReferences->At(itrd)))->X();
	  Double_t y   = ((AliTrackReference*)(info->fTRDReferences->At(itrd)))->Y();
	  Double_t phi    = TMath::ATan2(y,x);	
	  //	  Double_t snp    =  TMath::ATan2(((AliTrackReference*)(info->fTRDReferences->At(itrd)))->Py(),
	  //				  ((AliTrackReference*)(info->fTRDReferences->At(itrd)))->Px());
	  //	  if (TMath::Abs(snp)>0.9 &&cradius>10) break; 
	  Double_t radius = TMath::Sqrt(x*x+y*y);
	  if (radius>cradius){
	    cradius = radius;
	    cphi    = phi;
	    ref     =  ((AliTrackReference*)(info->fTRDReferences->At(itrd)));
	  }else break;   
	}
	//
	// rotate and propagate to local system of the last track reference
	//
	trd1->SetStop(kFALSE);
	if (!trd1->Rotate(cphi-trd1->GetAlpha())) trd1->SetStop(kTRUE);
	else{
	  if (!trd1->PropagateTo(cradius)) trd1->SetStop(kTRUE);
	}
      }    
    }
    tree3->Fill();    
  }
  f3.cd();
  tree3->Write();
}


