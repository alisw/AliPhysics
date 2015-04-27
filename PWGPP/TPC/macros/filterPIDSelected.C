/*
  Macro to make additional filter of the V0 to select clean sample of identified particles.
  As an input the   

  .x $HOME/rootlogon.C   
  .L $ALICE_PHYSICS/PWGPP/TPC/macros/filterPIDSelected.C+

 */

  



#include "TFile.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TH2.h"
#include "TF1.h"
#include "TTreeStream.h"
#include "AliMathBase.h"
#include "TSystem.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TRandom.h"
#include "AliTPCcalibBase.h"
#include "TCanvas.h"
#include "TLegend.h"
//
#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "TMath.h"
#include "AliXRDPROOFtoolkit.h"
#include "TStatToolkit.h"
#include "TCut.h"

TTree * tree  = 0;
TTreeSRedirector *pcstream = 0; //new TTreeSRedirector("trend.root");
//
void filterPIDSelectedTOF( const char * chfinput="highptAll.list");
void filterPIDSelectedV0( const char * chfinput="highptAll.list");

void filterPIDSelected( const char * chfinput="highptAll.list"){
  filterPIDSelectedTOF(chfinput);
  filterPIDSelectedV0(chfinput);
}

void filterPIDSelectedV0( const char * chfinput){
  //
  // Code to select identified V0 for the PID 
  // As an input chain of filter trees is used
  // Parameter:
  //   finput - name of the list file or the name of file itself
  // Oputput:
  //   file - V0Selected.root
  //
  //
  TTree * chain  = 0;  
  if (TString(chfinput).Contains(".list")) {
    chain = AliXRDPROOFtoolkit::MakeChainRandom(chfinput,"V0s",0,1000);
  }else{
    TFile * finput= TFile::Open(chfinput);
    if (!finput) finput= TFile::Open(TString::Format("%s#FilterEvents_Trees.root",finput));
    chain=(TTree*)finput->Get("V0s");
  }  
  chain->SetCacheSize(1000000000);
  //
  TDatabasePDG pdg;
  Double_t massLambda = pdg.GetParticle("Lambda0")->Mass();
  Double_t massK0 = pdg.GetParticle("K0")->Mass();
  Double_t massPion = pdg.GetParticle("pi+")->Mass();
  Double_t massProton = pdg.GetParticle("proton")->Mass();
  //
  //
  chain->SetAlias("massPion",Form("(%f+0)",massPion));
  chain->SetAlias("massProton",Form("(%f+0)",massProton));
  chain->SetAlias("massK0",Form("(%f+0)",massK0));
  chain->SetAlias("massLambda",Form("(%f+0)",massLambda));
  // delta of mass
  chain->SetAlias("K0Delta","(v0.GetEffMass(2,2)-massK0)");
  chain->SetAlias("LDelta","(v0.GetEffMass(4,2)-massLambda)");
  chain->SetAlias("ALDelta","(v0.GetEffMass(2,4)-massLambda)");
  chain->SetAlias("EDelta","(v0.GetEffMass(0,0))");
  // pull of the mass
  chain->SetAlias("K0Pull","(v0.GetEffMass(2,2)-massK0)/v0.GetKFInfo(2,2,1)");
  chain->SetAlias("LPull","(v0.GetEffMass(4,2)-massLambda)/v0.GetKFInfo(4,2,1)");
  chain->SetAlias("ALPull","(v0.GetEffMass(2,4)-massLambda)/v0.GetKFInfo(2,4,1)");
  chain->SetAlias("EPull","EDelta/v0.GetKFInfo(0,0,1)");
  // effective pull of the mass - (empirical values form fits)
  chain->SetAlias("K0PullEff","K0Delta/sqrt((3.63321e-03)**2+(5.68795e-04*v0.Pt())**2)");
  chain->SetAlias("LPullEff","LDelta/sqrt((1.5e-03)**2+(1.8e-04*v0.Pt())**2)");
  chain->SetAlias("ALPullEff","ALDelta/sqrt((1.5e-03)**2+(1.8e-04*v0.Pt())**2)");
  chain->SetAlias("EPullEff","v0.GetEffMass(0,0)/sqrt((5e-03)**2+(1.e-04*v0.Pt())**2)");
  //
  //
  chain->SetAlias("dEdx0DProton","AliMathBase::BetheBlochAleph(track0.fIp.P()/massProton)");
  chain->SetAlias("dEdx1DProton","AliMathBase::BetheBlochAleph(track1.fIp.P()/massProton)");
  chain->SetAlias("dEdx0DPion","AliMathBase::BetheBlochAleph(track0.fIp.P()/massPion)");
  chain->SetAlias("dEdx1DPion","AliMathBase::BetheBlochAleph(track1.fIp.P()/massPion)");
  //
  // V0 - cuts -PID, 
  //	
  chain->SetAlias("cutDist","sqrt((track0.fIp.fP[0]-track1.fIp.fP[0])**2+(track0.fIp.fP[1]-track1.fIp.fP[1])**2)>3");
  chain->SetAlias("cutLong","track0.GetTPCClusterInfo(3,1,0)+5*abs(track0.fP[4])>130&&track1.GetTPCClusterInfo(3,1,0)>130-5*abs(track1.fP[4])");
  chain->SetAlias("cutPID","track0.fTPCsignal>0&&track1.fTPCsignal>0");
  chain->SetAlias("cutResol","sqrt(track0.fC[14]/track0.fP[4])<0.15&&sqrt(track1.fC[14]/track1.fP[4])<0.15");
  chain->SetAlias("cutV0","cutPID&&cutDist&&cutLong&&cutResol");	
  //
  //  
  chain->SetAlias("K0Selected",      "abs(K0Pull)<3. &&abs(K0PullEff)<3.  && abs(LPull)>3  && abs(ALPull)>3  &&v0.PtArmV0()>0.11"); 
  chain->SetAlias("LambdaSelected",  "abs(LPull)<3.  &&abs(LPullEff)<3.   && abs(K0Pull)>3 && abs(EPull)>3 && abs(EDelta)>0.05");  
  chain->SetAlias("ALambdaSelected", "abs(ALPull)<3. &&abs(ALPullEff)<3   && abs(K0Pull)>3 && abs(EPull)>3 &&abs(EDelta)>0.05");
  //
  chain->SetAlias("GammaSelected", "abs(EPull)<3     && abs(K0Pull)>3 && abs(LPull)>3 && abs(ALPull)>3");
  //
  //
  TFile *fselected = TFile::Open("V0Selected.root","recreate");
  TTree * treeK0     =   chain->CopyTree("type==8&&cutV0&&K0Selected");
  TTree * treeLambda =   chain->CopyTree("type==4&&cutV0&&LambdaSelected");
  TTree * treeALambda =   chain->CopyTree("type==2&&cutV0&&ALambdaSelected");
  TTree * treeGamma =   chain->CopyTree("type==1&&cutV0&&GammaSelected");
  //
  TTree * trees[4]={treeK0,treeLambda, treeGamma,treeALambda};
  TList * aliases = chain->GetListOfAliases();
  Int_t nalias= aliases->GetEntries();

  for (Int_t i=0; i<4; i++){
    for (Int_t ialias=0; ialias<nalias; ialias++){
      TNamed *alias = (TNamed*)aliases->At(ialias);
      trees[i]->SetAlias(alias->GetName(),alias->GetTitle());
    }
  }  
  treeK0->Write("treeK0");
  treeLambda->Write("treeLambda");
  treeALambda->Write("treeALambda");
  treeGamma->Write("treeGamma");
  fselected->Close();
  //
}

void filterPIDSelectedTOF( const char * chfinput){
  //
  // Code to select identified V0 for the PID 
  // As an input chain of filter trees is used
   //
  TTree * chain  = 0;  
  if (TString(chfinput).Contains(".list")) {
    chain = AliXRDPROOFtoolkit::MakeChainRandom(chfinput,"highPt",0,1000);
  }else{
    TFile * finput= TFile::Open(chfinput);
    if (!finput) finput= TFile::Open(TString::Format("%s#FilterEvents_Trees.root",finput));
    chain=(TTree*)finput->Get("highPt");
  }  
  chain->SetCacheSize(1000000000);
  //
  TDatabasePDG pdg;
  Double_t massLambda = pdg.GetParticle("Lambda0")->Mass();
  Double_t massK0 = pdg.GetParticle("K0")->Mass();
  Double_t massPion = pdg.GetParticle("pi+")->Mass();
  Double_t massProton = pdg.GetParticle("proton")->Mass();
  chain->SetAlias("cutLong","esdTrack.GetTPCClusterInfo(3,1,0)+5*abs(esdTrack.fP[4])>130");
  TFile *fselected  = TFile::Open("TOFSelected.root","recreate");

  TCut cutDeltaProton="abs((esdTrack.fTrackTime[4]-esdTrack.fTrackTime[3]))>400&&abs(esdTrack.fTOFsignalDz<3)";
  TCut cutDeltaKaon="abs((esdTrack.fTrackTime[3]-esdTrack.fTrackTime[2]))>400&&abs(esdTrack.fTOFsignalDz<3)";
  TCut cutDeltaPion="abs((esdTrack.fTrackTime[2]-esdTrack.fTrackTime[0]))>200&&abs(esdTrack.fTOFsignalDz<3)";

  TTree * treeEl  = chain->CopyTree(cutDeltaPion+"cutLong&&esdTrack.fTOFr[0]>0.3+max(2*max(esdTrack.fTOFr[1],esdTrack.fTOFr[3]),esdTrack.fTOFr[4])");
  TTree * treePion  = chain->CopyTree(cutDeltaPion+"cutLong&&esdTrack.fTOFr[2]>0.3+max(max(esdTrack.fTOFr[0],esdTrack.fTOFr[3]),esdTrack.fTOFr[4])");
  TTree * treeKaon  = chain->CopyTree(cutDeltaKaon+"cutLong&&esdTrack.fTOFr[3]>0.3+max(max(esdTrack.fTOFr[0],2*esdTrack.fTOFr[2]),esdTrack.fTOFr[4])");
  TTree * treeProton  = chain->CopyTree(cutDeltaProton+"cutLong&&esdTrack.fTOFr[4]>0.3+max(max(esdTrack.fTOFr[0],2*esdTrack.fTOFr[2]),esdTrack.fTOFr[3])");
  treeEl->Write("treeEl");
  treePion->Write("treePion");
  treeKaon->Write("treeKaon");
  treeProton->Write("treeProton");
  //
  fselected->Close();
}

void FitPIDNCLSelected(){
  //
  // fit the probability to find the cluster
  //
  Int_t kmarkers[5]={20,21,25,24,22};
  Int_t kcolors[5]={1,2,4,5,7};
  const char *chname[5]={"Proton","Pion","Electron"};

  TFile f("V0Selected.root");
  TTree * treeLambda = (TTree*)f.Get("treeLambda"); 
  TTree * treeK0 = (TTree*)f.Get("treeK0"); 
  TTree * treeGamma = (TTree*)f.Get("treeGamma"); 
  //
  //
  TH2 *his3D=0;
  TObjArray  fitArray(3); 
  TH1D * hisNcldEdx[3]={0,0,0};
  TCut cutOut="abs(track0.fP[4])<2";
  treeLambda->Draw("(1-track0.GetTPCClusterInfo(2,0)):sqrt(1+track0.fP[3]**2)*AliMathBase::BetheBlochAleph(track0.fIp.P()/massProton)>>hisNclProton(50,1,2.)",cutOut,"prof");
  hisNcldEdx[0]=(TH1D*)treeLambda->GetHistogram()->Clone();
  treeK0->Draw("(1-track0.GetTPCClusterInfo(2,0)):sqrt(1+track0.fP[3]**2)*AliMathBase::BetheBlochAleph(track0.fIp.P()/massPion)>>hisNclPion(50,1,2.)",cutOut,"prof");
  hisNcldEdx[1]=(TH1D*)treeK0->GetHistogram()->Clone();
  treeGamma->Draw("(1-track0.GetTPCClusterInfo(2,0)):sqrt(1+track0.fP[3]**2)*AliMathBase::BetheBlochAleph(track0.fIp.P()/0.0005)>>hisNclPion(50,1,2.5)",cutOut,"prof");
  hisNcldEdx[2]=(TH1D*)treeGamma->GetHistogram()->Clone();

  TF1 fnclQ("fnclQ","[0]+[1]*exp(-[2]*abs(x))");
  fnclQ.SetParameters(0.02,5,3);
  hisNcldEdx[0]->Fit(&fnclQ);
  hisNcldEdx[1]->Fit(&fnclQ);
  //hisNcldEdx[2]->Fit(&fnclQ);

  TCanvas * canvas = new TCanvas("canvasNCL0","canvasNCL0",600,500);
  TLegend *legend = new TLegend(0.5,0.6,0.89,0.89,"Cluster finder eff.");
  legend->SetBorderSize(0);
  for (Int_t i=0; i<3; i++){
    hisNcldEdx[i]->SetMinimum(0);
    hisNcldEdx[i]->SetMaximum(0.3);
    hisNcldEdx[i]->SetMarkerColor(kcolors[i]);
    hisNcldEdx[i]->SetMarkerStyle(kmarkers[i]);
    hisNcldEdx[i]->GetXaxis()->SetTitle("Q ~ #sqrt{1+tan(#theta)^2}dE/dx");
    hisNcldEdx[i]->GetYaxis()->SetTitle("p_{cl}");
    if (i==0)hisNcldEdx[i]->Draw("");
    hisNcldEdx[i]->Draw("same");
    legend->AddEntry(hisNcldEdx[i],chname[i]);
  }
  legend->Draw();
  
}
