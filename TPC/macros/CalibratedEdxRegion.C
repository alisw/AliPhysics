/*
  Macro to make an per region dEdx calibration.
  Here we assue we have suufiecint ammount of V) PID selected tracks for the calibration

  .x $HOME/rootlogon.C
  .L $ALICE_ROOT/TPC/macros/CalibratedEdxRegion.C+
  
 */

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h" 
#include "TCut.h"
#include "TStatToolkit.h"
#include "TGraphErrors.h"

#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliExternalTrackParam.h"
#include "AliTPCdEdxInfo.h"
#include "AliTPCParamSR.h"

//TTree * trees[4]={treeK0, treeLambda,treeALambda,treeGamma};
TTree * trees[4]={0};
TFile * calibrationFile=0;

void InitTrees(const char * v0file= "/hera/alice/miranov/ExpertQA/data/LHC13c/pass1/V0Selected.root"){
  //
  //
  //
  //  const char * v0file = "/hera/alice/miranov/ExpertQA/data/LHC13c/pass1/V0Selected.root";
  TFile * f = TFile::Open(v0file);
  TTree * treeLambda=(TTree*)f->Get("treeLambda");
  TTree * treeALambda=(TTree*)f->Get("treeALambda");
  TTree * treeK0=(TTree*)f->Get("treeK0");
  TTree * treeGamma=(TTree*)f->Get("treeGamma");
  calibrationFile=TFile::Open("dEdxCalibration.root","update");
  //
  // register BetheBloch parameterization for each type of identified particles
  //
  AliTPCParam param;
  param.RegisterBBParam(param.GetBetheBlochParamAlice(),1);
  treeK0->SetAlias("dEdxExp0","AliTPCParam::BetheBlochAleph(track0.fIp.P()/massPion,1)");
  treeK0->SetAlias("dEdxExp1","AliTPCParam::BetheBlochAleph(track1.fIp.P()/massPion,1)");
  treeLambda->SetAlias("dEdxExp0","AliTPCParam::BetheBlochAleph(track0.fIp.P()/massProton,1)");
  treeLambda->SetAlias("dEdxExp1","AliTPCParam::BetheBlochAleph(track1.fIp.P()/massPion,1)");
  treeALambda->SetAlias("dEdxExp0","AliTPCParam::BetheBlochAleph(track0.fIp.P()/massPion,1)");
  treeALambda->SetAlias("dEdxExp1","AliTPCParam::BetheBlochAleph(track1.fIp.P()/massProton,1)");
  treeGamma->SetAlias("dEdxExp0","AliTPCParam::BetheBlochAleph(track0.fIp.P()/0.005,1)");
  treeGamma->SetAlias("dEdxExp1","AliTPCParam::BetheBlochAleph(track1.fIp.P()/0.005,1)");
  //
  //
  //
  AliTPCParamSR paramSR;
  TString rROC0=TString::Format("%2.2f", 0.5*(paramSR.GetInnerRadiusLow()+paramSR.GetInnerRadiusUp()));
  TString rROC1=TString::Format("%2.2f", 0.5*(paramSR.GetPadRowRadii(36,0)+paramSR.GetPadRowRadii(36,paramSR.GetNRowUp1()-1)));
  TString rROC2=TString::Format("%2.2f", 0.5*(paramSR.GetPadRowRadii(36,0)+paramSR.GetPadRowRadii(36,paramSR.GetNRowUp()-1)));
  //
  //
  // 
  Double_t bz=-5;
  trees={treeK0, treeLambda,treeALambda,treeGamma};
  for (Int_t itree=0; itree<4;itree++){
    trees[itree]->SetAlias("bz","-5");
    trees[itree]->SetAlias("phi0ROC0",TString::Format("track0.fIp.GetParameterAtRadius(%s+0,%2.2f,7)",rROC0.Data(),bz));
    trees[itree]->SetAlias("phi0ROC1",TString::Format("track0.fIp.GetParameterAtRadius(%s+0,%2.2f,7)",rROC1.Data(),bz));
    trees[itree]->SetAlias("phi0ROC2",TString::Format("track0.fIp.GetParameterAtRadius(%s+0,%2.2f,7)",rROC2.Data(),bz));
    trees[itree]->SetAlias("phi1ROC0",TString::Format("track1.fIp.GetParameterAtRadius(%s+0,%2.2f,7)",rROC0.Data(),bz));
    trees[itree]->SetAlias("phi1ROC1",TString::Format("track1.fIp.GetParameterAtRadius(%s+0,%2.2f,7)",rROC1.Data(),bz));
    trees[itree]->SetAlias("phi1ROC2",TString::Format("track1.fIp.GetParameterAtRadius(%s+0,%2.2f,7)",rROC2.Data(),bz));
    //
    //
    trees[itree]->SetAlias("side0","(track0.fIp.fP[3]>0&&track0.fIp.fP[1]>0)+2*(track0.fIp.fP[3]<0&&track0.fIp.fP[1]<0)");  // 0 - both sides, 1- A side, 2- C side
    trees[itree]->SetAlias("side1","(track1.fIp.fP[3]>0&&track1.fIp.fP[1]>0)+2*(track1.fIp.fP[3]<0&&track1.fIp.fP[1]<0)");
  }
  /*
    consistency check: 
    // Phi position extrapolated to the IFC OK
    treeK0->Draw("track0.GetParameterAtRadius(83.8,-5,7)-track0.fIp.fAlpha-track0.fIp.fP[0]/track0.fIp.fX>>his(100,-0.02,0.02)","abs(track0.GetParameterAtRadius(83.8,-5,7)-track0.fIp.fAlpha)<0.5","",100000);
    // Phi angle extrapolated to the IFC
    treeK0->Draw("track0.GetParameterAtRadius(83.8,-5,8)-track0.fIp.Phi()>>his(100,-0.02,0.02)","abs(track0.GetParameterAtRadius(83.8,-5,7)-track0.fIp.fAlpha)<0.5","",100000);
    //
    treeK0->Draw("track0.GetParameterAtRadius(103.8,-5,9)-(track0.fIp.GetParameterAtRadius(103.8,-5,8)-track0.fIp.GetParameterAtRadius(103.8,-5,7))>>his(100,-0.05,0.05)","abs(track0.GetParameterAtRadius(103.8,-5,8)-track0.GetParameterAtRadius(103.8,-5,7))<0.5","",10000);

   */
}

void MakeSectorCalibration(){
  //
  // Get sector phi correction - separate for 
  //
  TCut cutASide0="side0>0";
  TCut cutPt0="abs(track0.fIp.fP[4])<3";
  TCut cutASide1="side1>0";
  TCut cutPt1="abs(track0.fIp.fP[4])<3";
  TObjArray * histoArray = new TObjArray(100);
  TObjArray * graphArray = new TObjArray(100);
  //
  //
  //
  for (Int_t iregion=0; iregion<3; iregion++){
    //    
    TH2F *hisSectorROC=0;
    TH2F *hisSectorROCP=0;
    TH2F *hisSectorROCM=0;
    TGraphErrors * graphs[3]={0};

    for (Int_t itree=0; itree<3; itree++){
      TString var0=TString::Format("track0.fTPCdEdxInfo.fTPCsignalRegion[%d]/dEdxExp0:18*(phi0ROC%d/pi+(side0-1))>>hisSectorROCPlus%d(36,0,36,60,20,80)",iregion,iregion,iregion);
      TString var1=TString::Format("track1.fTPCdEdxInfo.fTPCsignalRegion[%d]/dEdxExp1:18*(phi1ROC%d/pi+(side1-1))>>hisSectorROCMinus%d(36,0,36,60,20,80)",iregion,iregion,iregion);
      
      trees[itree]->Draw(var0.Data(),cutASide0+cutPt0,"colzgoff");
      if (hisSectorROCP==0) {
	hisSectorROCP=(TH2F*)trees[itree]->GetHistogram()->Clone();
      }
      if (hisSectorROCP) hisSectorROCP->Add((TH2F*)trees[itree]->GetHistogram());
      trees[itree]->Draw(var1.Data(),cutASide1+cutPt1,"colzgoff");
      if (hisSectorROCM==0) hisSectorROCM=(TH2F*)trees[itree]->GetHistogram()->Clone();
      if (hisSectorROCM) hisSectorROCM->Add((TH2F*)trees[itree]->GetHistogram());      
    }
    hisSectorROC=(TH2F*)hisSectorROCP->Clone();
    hisSectorROC->SetName(TString::Format("hisSectorROCBoth%d(36,0,36,60,20,80)",iregion));
    hisSectorROC->Add(hisSectorROCM);
    //
    graphs[0]=TStatToolkit::MakeStat1D(hisSectorROC, 0, 0.85,4,20,1);    
    graphs[1]=TStatToolkit::MakeStat1D(hisSectorROCP, 0, 0.85,4,24,2);
    graphs[2]=TStatToolkit::MakeStat1D(hisSectorROCM, 0, 0.85,4,25,4);    
    graphs[0]->SetName(TString::Format("graphSectorROCBoth%d(36,0,36,60,20,80)",iregion));
    graphs[1]->SetName(TString::Format("graphSectorROCPlus%d(36,0,36,60,20,80)",iregion));
    graphs[2]->SetName(TString::Format("graphSectorROCMinus%d(36,0,36,60,20,80)",iregion));
    graphs[0]->Draw("alp");
    TStatToolkit::MakeStat1D(hisSectorROCP, 0, 0.85,4,24,2)->Draw("lp");
    TStatToolkit::MakeStat1D(hisSectorROCM, 0, 0.85,4,25,4)->Draw("lp");
    histoArray->AddLast(hisSectorROC);
    histoArray->AddLast(hisSectorROCP);
    histoArray->AddLast(hisSectorROCM); 
    for (Int_t i=0; i<3; i++){
      graphArray->AddLast(graphs[i]);
    }
  }
  calibrationFile->mkdir("histos");
  calibrationFile->cd("histos");
  histoArray->Write();
  calibrationFile->mkdir("graphs");
  calibrationFile->cd("graphs");
  graphArray->Write();
  
}
