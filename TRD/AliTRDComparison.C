/*

//load mag field
AliRunLoader *loader = AliRunLoader::Open("galice.root");
loader->LoadgAlice();
gAlice = loader->GetAliRun();
AliMagF * field = gAlice->Field();
AliKalmanTrack::SetFieldMap(gAlice->Field());
AliExternalTrackParam::SetFieldMap(gAlice->Field());
TFile fg("geometry.root");
gGeoManager = (TGeoManager*)fg.Get("Geo");

.L AliGenInfo.C+
.L AliESDComparisonMI.C+
.L AliTRDComparison.C+
//.L AliESDComparisonPIC.C+
.x TDR_style.C

MakeCompTr();
MakeTree();



MakeComp();
MakeAlias();
TFile fff("TRDdebug.root");
TTree * tree =(TTree*)fff.Get("tracklet");
TTree * tree2 =(TTree*)fff.Get("Tracks");

tree->Draw("mean0:tracklet.fZ-track.fZ","plane==0&&abs(angleb)<1","");



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
#include "AliPoints.h"
#include "AliESDtrack.h"
#include "AliTRDtrack.h"
#include "AliITStrackMI.h"
#include "AliESDV0MI.h"
#include "AliESDkink.h"
#include "AliTRDparameter.h"
#include "AliTRDgeometryFull.h"
#include "AliTRDcluster.h"

// COMPARISON includes
#include "../STEER/AliGenInfo.h"
#include "../STEER/AliESDComparisonMI.h"



AliComparisonDraw comp;
AliComparisonDraw comptr;
AliTRDparameter   par("pica","vyjebana");
AliTRDgeometryFull    geom;
AliRunLoader  *fLoader=0; //= AliRunLoader::Open(fnGalice);

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
  ClassDef(AliTRDInfoMC,1)   // TRD MC information 
};

ClassImp(AliTRDInfoMC)


class AliTRDInfoCl: public TObject{
  public:
  AliTRDcluster fCl[6][20];      // clusters
  Int_t  fNClusters;           // number of clusters
  ClassDef(AliTRDInfoCl,1)   // TRD MC information 
};

ClassImp(AliTRDInfoCl)



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
  //  Float_t radius[6]={305
  for (Int_t i=0;i<info->fNTRDRef;i++){
    AliTrackReference * ref = (AliTrackReference*)info->fTRDReferences->At(i);
    Float_t r = ref->R();
    for (Int_t jplane = 0; jplane<6; jplane++){
      if (TMath::Abs(r - par.GetTime0(jplane))<6){
      //      if (TMath::Abs(r - 305+jplane*14.)<6){
	fLayer[jplane]++;
	for (Int_t kplane=jplane;kplane>=0;kplane--){
	  fLayerA[kplane]++;
	  fRef[jplane]=*ref;
	}
      }
    }
  }
  for (Int_t jplane = 0; jplane<6; jplane++){
    if (fLayer[jplane]>0) fNLayers++;
    if (fLayer[jplane]>0) fLast = jplane;
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



void PtRes(){
  
  
}


TH1F * MakeCumul(TH1F *his)
{
  TH1F *hcumul = (TH1F*)his->Clone();
  Float_t all = (Float_t)his->GetSum();
  for (Int_t i=0;i<=his->GetNbinsX();i++){
    hcumul->SetBinContent(i,his->Integral(0,i)/all);
  } 
  hcumul->SetFillColor(0);
  return hcumul;
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
  comp.fTree->SetAlias("phi_tof","atan2(MC.fTOFReferences[0].fPy,MC.fTOFReferences[0].fPx)");
  comp.fTree->SetAlias("phi_trd","asin(TRD0.GetSnp())+TRD0.fAlpha");
  comp.fTree->SetAlias("pullphi0","(TRD0.GetSnp()+TRD0.fAlpha-atan2(MC.fTOFReferences[0].fPy,MC.fTOFReferences[0].fPx))/sqrt(TRD0.fCee)");
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
  comp.fTree->SetAlias("TPCE","sqrt(MC.fTPCReferences[4].P()**2+MC.fMass**2)");
  comp.fTree->SetAlias("TOFE","sqrt(MC.fTOFReferences[0].P()**2+MC.fMass**2)");
  comp.fTree->SetAlias("TRDE","sqrt(TR.P()**2+MC.fMass**2)");
  comp.fTree->SetAlias("dtpi","TRD0.fIntegratedTime[2]-10^12*MC.fTOFReferences.fTime");
  comp.fTree->SetAlias("dtk","TRD0.fIntegratedTime[3]-10^12*MC.fTOFReferences.fTime");
  comp.fTree->SetAlias("dtp","TRD0.fIntegratedTime[4]-10^12*MC.fTOFReferences.fTime");
  comp.fTree->SetAlias("dt","(abs(MC.fPdg)==2212)*dtp+(abs(MC.fPdg)==211)*dtpi+(abs(MC.fPdg)==321)*dtk"); //delta time
  comp.fTree->SetAlias("Beta","sqrt(TR.P()**2/(TR.P()**2+MC.fMass**2))");
  comp.fTree->SetAlias("dtsigma","10+(1-Beta)*280");
  comp.fTree->SetAlias("dptrel","(abs(TRD0.GetPt())-MC.fTOFReferences[0].Pt())/MC.fTOFReferences[0].Pt()");
  //
  comp.fTree->SetAlias("dphi","atan2(TR.fY,TR..fX)-atan2(TR.fPy,TR..fPx)");
  comp.fTree->SetAlias("dphi0","acos((TR.fX*TR.fPx+TR.fY*TR.fPy)/(TR.Pt()*TR.R()))");


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


void DrawTOFResY(){
  //
  comp.DrawXY("TR.Pt()","TRD0.fY",c2,"1",10,0.3,2.,-2,2,30);
  TH1F * histof = (TH1F*)comp.fRes->Clone();
  comp.DrawXY("TR.Pt()","TRD1.fY",c2,"1",10,0.3,2.,-2,2,30);
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

void DrawTOFResZ(){
  //
  comp.DrawXY("TR.Pt()","TRD0.fZ-MC.fTOFReferences.fZ",c2,"1",10,0.3,2.,-4,4,30);
  TH1F * histof = (TH1F*)comp.fRes->Clone();
  comp.DrawXY("TR.Pt()","TRD1.fZ-TR.fZ",c2,"1",10,0.3,2.,-4,4,30);
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
  TH1F * hisy = PullCumul("abs(pully0)",cut+cprim+c0+c2+c3,1000,20);
  TH1F * hisz = PullCumul("abs(pullz0)",cut+cprim+c0+c2+c3,1000,20);
  TH1F * hisyz = PullCumul("sqrt(pullz0**2+pully0**2)",cut+cprim+c0+c2+c3,1000,20);
  hisz->SetXTitle("Pull (unit)");
  hisz->SetYTitle("Cumulative function");
  hisy->SetLineColor(1);hisz->SetLineColor(2);hisyz->SetLineColor(3);
  hisz->Draw(); hisy->Draw("same");hisyz->Draw("same");  
  TLegend *legend = new TLegend(0.55,0.12,0.85,0.35, "TRD pull (unit)");
  legend->SetBorderSize(1);
  legend->AddEntry(hisy ,"Pull y");
  legend->AddEntry(hisz,"Pull z");
  legend->AddEntry(hisyz,"Pull yz");
  legend->Draw();
}

void MakeCumulAbs(TCut & cut)
{
  TH1F * hisy = PullCumul("abs(TRD0.fY)",cut+cprim+c0+c2+c3,1000,6);
  TH1F * hisz = PullCumul("abs(deltaz0)",cut+cprim+c0+c2+c3,1000,6);  
  TH1F * hisyz = PullCumul("sqrt(TRD0.fY**2+deltaz0**2)",cut+cprim+c0+c2+c3,1000,6);
  hisz->SetXTitle("Delta (cm)");
  hisz->SetYTitle("Cumulative function");
  //
  hisy->SetLineColor(1);hisz->SetLineColor(2);hisyz->SetLineColor(3);
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
      Int_t label = cl->GetLabel(0);
      Int_t detector=cl->GetDetector();
      Int_t localTimeBin=cl->GetLocalTimeBin();
      //      Int_t sector=geom.GetSector(detector);
      Int_t plane=geom.GetPlane(detector);
      Int_t pos = label;      
      if (label<0) continue;
      if (localTimeBin>19||localTimeBin<0) continue;
      if (plane<0||plane>5) continue;
      
      AliTRDInfoCl * info =( AliTRDInfoCl*)arr->At(pos); 
      if (!info) {
	info = new AliTRDInfoCl;
	arr->AddAt(info,pos);
      } 
     
      info->fCl[plane][localTimeBin]=*cl;
    }
  }  
}


void MakeTree()
{
  fLoader = AliRunLoader::Open("galice.root");
  AliESDtrack * esd     = new AliESDtrack;
  AliTRDtrack * trd     = new AliTRDtrack;
  AliTRDtrack * trd0    = new AliTRDtrack;
  AliTRDtrack * trd1    = new AliTRDtrack;
  AliTRDtrack * trdb    = new AliTRDtrack;
  AliMCInfo* info       = new AliMCInfo;
  AliTRDInfoMC *trdinfo = new AliTRDInfoMC;
  AliTRDInfoCl *infocl = new AliTRDInfoCl;
  AliTrackReference *ref = new AliTrackReference;
  TObjArray         * clarray = new TObjArray(6*1000000);
  Int_t       eventNr   =0;
  const Int_t kmaxlabel=1000000;
  //
  // MC part
  //
  TFile f("cmpESDTracks.root");
  TTree * tree1 = (TTree*) f.Get("ESDcmpTracks");
  TBranch * brmc  = tree1->GetBranch("MC");
  brmc->SetAddress(&info);
  //
  // 
  TObjArray *esds = new TObjArray(kmaxlabel);
  TObjArray *trds = new TObjArray(kmaxlabel);
  //
  //  rec part
  //
  TFile f2("TRDdebug.root");
  TTree * tree2 =(TTree*)f2.Get("Tracks");
  TBranch * br20 = tree2->GetBranch("ESD.");
  TBranch * br21 = tree2->GetBranch("trd.");
  TBranch * br22 = tree2->GetBranch("EventNr");
  br20->SetAddress(&esd);
  br21->SetAddress(&trd);  
  br22->SetAddress(&eventNr);
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
  //
  //
  Int_t lastevent=-1;

  for (Int_t i=0;i<tree1->GetEntries();i++){
    brmc->GetEntry(i);
    if (lastevent!=info->fEventNr){
      printf("Read Event\t%d\n", info->fEventNr);
      //
      // read new event
      //
      lastevent = info->fEventNr;
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
      clarray->Clear();
      LoadClusters(clarray,lastevent);
      //
      AliESDtrack * esdp     = new AliESDtrack;
      AliTRDtrack * trdp     = new AliTRDtrack;

      br20->SetAddress(&esdp);
      br21->SetAddress(&trdp);
      //
      // Load esd and trd tracks to memory
      Int_t ntracks0 =0;
      Int_t ntracks1 =0;
      Int_t ntracks2 =0;
      
      for (Int_t ip=0;ip<tree2->GetEntries();ip++){
	tree2->GetEntry(ip);
	if (eventNr!=info->fEventNr) continue;
	//	if (eventNr>info->fEventNr) continue;
	//
	Int_t label = TMath::Abs(esdp->GetTPCLabel());
	if (label==0) continue;
	if (label>=kmaxlabel) continue;
	ntracks0++;
	if (!(esds->At(label))){
	  esds->AddAt(new AliESDtrack(*esdp),label);
	  trds->AddAt(new AliTRDtrack(*trdp),label);
	  ntracks1++;
	} else{
	  if (esd->GetStatus() & AliESDtrack::kITSin){
	    esds->AddAt(new AliESDtrack(*esdp),label);
	    trds->AddAt(new AliTRDtrack(*trdp),label);
	    ntracks2++;
	  }
	}
      }
      printf(" readed\t%d\t%d\t%d\n",ntracks0, ntracks1,ntracks2);
    }

    //    if (info->fNTOFRef==0 && info->fNTRDRef<5) continue;
    Int_t label = info->fLabel;    
    trdinfo->Update(info);
    if (clarray->At(label)){
      infocl = (AliTRDInfoCl*)clarray->At(label);
    }
    if (info->fNTPCRef<3) continue;
    trd0->SetStop(kTRUE);
    trdb->SetStop(kTRUE);
    trd1->SetStop(kTRUE);
    //
    if (esds->At(label)) {
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
	  AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);
	  if (!trd0->PropagateTo(radius,param[1],param[0])) trd0->SetStop(kTRUE);
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
      tree3->Fill();
    }
  }
  f3.cd();
  tree3->Write();
}


