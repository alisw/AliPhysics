/*
 gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TPC -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TOF -I$ALICE_ROOT/RAW -I$ALICE_ROOT/STAT");

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
#include "TGraph.h"
#include "AliCDBManager.h"
#include "AliTPCExBBShape.h"
#include "TRandom.h"
#include "AliGeomManager.h"
#include "AliESDVertex.h"
#include "AliTPCcalibDB.h"
#endif


void AliTPCCorrectionDemo() {

  //
  // This is a Demo function of the general class AliTPCCorrection, which is used for 
  // general space point correction due to different effects.
  // The effects used in this Demo are:
  //   1. ExB twist - general offset of the TPC axis in comparison to the B field axis
  //   2. GG error (Gating Grid volt. error) - not perfectly aligned GG voltage (in terms of voltage)
  //   3. ExBBShape - B field shape correction of the secound order
  //
  // See class descriptions for further details 
  //
  // Authors: Magnus Mager, Stefan Rossegger, Jim Thomas
  //
  //
  // omegaTau (wt) of the langevin equation
  // This is a function of the drift vel., the magnetic and electric field
  // e.g. vd=2.6 cm/usc; Ez=400 V/cm; Bz=0.5 T
  // wt =  -10.0*(Bz*10)*vd/Ez = -0.325 

  Double_t vdrift = 2.6; // [cm/us]   // to be updated: per second (ideally)
  Double_t bzField = -0.5; // [Tesla] // to be updated: per run
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 

  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  Double_t T1 = 0.9;
  Double_t T2 = 1.5;

  AliMagF mag("mag","mag");

  AliTPCExBTwist twist;
  twist.SetXTwist(0.001);
  
  AliTPCGGVoltError GGerror;
  GGerror.SetDeltaVGGA(50.);
  GGerror.SetDeltaVGGC(50.);
  GGerror.InitGGVoltErrorDistortion();

  AliTPCExBBShape exb;
  exb.SetBField(&mag);

  TObjArray cs;
  cs.Add(&twist);
  cs.Add(&GGerror);
  cs.Add(&exb);

  AliTPCComposedCorrection cc;
  cc.SetCorrections(&cs);
  cc.SetOmegaTauT1T2(wt,T1,T2);
  //cc.SetMode(1);

  cc.Print("DA"); // Print used correction classes

  TCanvas *c=new TCanvas;  // Plots
  c->Divide(2,2);
  c->cd(1);twist.CreateHistoDRPhiinZR(1.,100,100)->Draw("surf2");
  c->cd(2);GGerror.CreateHistoDRPhiinZR(1.,100,100)->Draw("surf2");
  c->cd(3);exb.CreateHistoDRPhiinZR(1.)->Draw("surf2");
  c->cd(4);cc.CreateHistoDRPhiinZR(1.)->Draw("surf2");
}



void MakeDistortionMap(){
  //
  // make distortiona map example for specific transformation 
  //
  Int_t run=0;
  TTreeSRedirector *pcstream =  new TTreeSRedirector("distort.root");
  Double_t vdrift = 2.6; // [cm/us]   // to be updated: per second (ideally)
  Double_t bzField = -0.5; // [Tesla] // to be updated: per run
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 

  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  Double_t T1 = 0.9;
  Double_t T2 = 1.5;
  AliMagF* magF= new AliMagF("Maps","Maps", bzField/0.5, 1.);
  TGeoGlobalMagField::Instance()->SetField(magF);

  AliCDBManager::Instance()->SetDefaultStorage("local:////$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(run);
  AliGeomManager::LoadGeometry();
  AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC");


  AliTPCExBTwist twist;
  twist.SetName("twistX001");
  twist.SetXTwist(0.001);
  
  AliTPCGGVoltError GGerror;
  GGerror.SetDeltaVGGA(50.);
  GGerror.SetDeltaVGGC(50.);
  GGerror.InitGGVoltErrorDistortion();
  GGerror.SetName("ggoffsetA50C50");
  AliTPCExBBShape exb;
  exb.SetBField(magF);
  exb.SetName("ExB");
  TObjArray cs;
  cs.Add(&twist);
  cs.Add(&GGerror);
  cs.Add(&exb);
  
  AliTPCComposedCorrection *cc = new AliTPCComposedCorrection;
  cc->SetCorrections(&cs);
  cc->SetOmegaTauT1T2(wt,T1,T2);
  //cc->SetMode(1);
  cc->SetName("composed");
  Double_t par[5]={0,0,0,0,0};
  Double_t cov[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for (Double_t theta=-1.; theta<1; theta+=0.1){
    for (Double_t alpha=-3.14; alpha<3.14; alpha+=0.2){
      par[0]=0;
      par[1]=theta*85;
      par[2]=0;
      par[3]=theta;
      par[4]=(gRandom->Rndm()-0.5);
      AliExternalTrackParam trackIn(85,alpha,par,cov);
      //FitDistortedTrack(&trackIn, cc, 85, -1,pcstream);
      twist.FitDistortedTrack(trackIn, 85, -1,pcstream);
      GGerror.FitDistortedTrack(trackIn, 85, -1,pcstream);
      exb.FitDistortedTrack(trackIn, 85, -1,pcstream);
      cc->FitDistortedTrack(trackIn, 85, -1,pcstream);      
    }
  }
  delete pcstream;
}

void DrawDistortionMap(){
  //
  // Example -drawing of distortion maps
  //
  TFile f("distort.root");
  TTree *   fitDistorttwistX001= (TTree*)f.Get("fitDistorttwistX001");
  TTree *   fitDistortggoffsetA50C50=(TTree*)f.Get("fitDistortggoffsetA50C50");
  TTree *   fitDistortExB=(TTree*)f.Get("fitDistortExB");
  TTree *   fitDistortcomposed=(TTree*)f.Get("fitDistortcomposed");
  fitDistortExB->SetMarkerColor(1);
  fitDistorttwistX001->SetMarkerColor(2);
  fitDistortggoffsetA50C50->SetMarkerColor(4);
  fitDistortExB->SetMarkerStyle(20);
  fitDistorttwistX001->SetMarkerStyle(21);
  fitDistortggoffsetA50C50->SetMarkerStyle(22);
  //
  // example draw delta local y/r-phi as function of fi
  Int_t entries=0;

  entries = fitDistortExB->Draw("track1.fP[0]-track0.fP[0]:track0.fAlpha","","");
  TGraph * grY0= new TGraph(entries, fitDistortExB->GetV2(),fitDistortExB->GetV1());
  entries = fitDistortggoffsetA50C50->Draw("track1.fP[0]-track0.fP[0]:track0.fAlpha","","");
  TGraph * grY1= new TGraph(entries, fitDistortggoffsetA50C50->GetV2(),fitDistortggoffsetA50C50->GetV1());
  entries = fitDistorttwistX001->Draw("track1.fP[0]-track0.fP[0]:track0.fAlpha","","");
  TGraph * grY2= new TGraph(entries, fitDistorttwistX001->GetV2(),fitDistorttwistX001->GetV1());
  entries = fitDistortcomposed->Draw("track1.fP[0]-track0.fP[0]:track0.fAlpha","","");
  TGraph * grY3= new TGraph(entries, fitDistortcomposed->GetV2(),fitDistortcomposed->GetV1());
  grY0->SetMinimum(-0.4);
  grY0->SetMaximum(0.4);
  grY0->SetMarkerStyle(20), grY0->SetMarkerColor(1);
  grY1->SetMarkerStyle(21), grY1->SetMarkerColor(2);
  grY2->SetMarkerStyle(22), grY2->SetMarkerColor(4);
  grY3->SetMarkerStyle(24), grY3->SetMarkerColor(5); 
  grY0->Draw("ap");
  grY1->Draw("p");
  grY2->Draw("p");
  grY3->Draw("p");
}


void TestVertex(){
  //
  //
  //  .x ConfigCalibTrain.C(120829)

  AliTPCComposedCorrection * corrC = ( AliTPCComposedCorrection *)AliTPCcalibDB::Instance()->GetTPCComposedCorrection(0.5);
  AliTPCCorrection * corrT = (AliTPCCorrection *)corrC->GetCorrections()->FindObject("exb_twist");
  AliTPCCorrection * corrS = (AliTPCCorrection *)corrC->GetCorrections()->FindObject("ExB");
  
  TTreeSRedirector *pcstream = new TTreeSRedirector("vertexDistort.root");
  TTreeSRedirector *pcstreamS = new TTreeSRedirector("vertexDistortS.root");
  TTreeSRedirector *pcstreamT = new TTreeSRedirector("vertexDistortT.root");
  Double_t orgVertex[3] ;
  AliESDVertex aV,cV,aVO,cVO;
  for (Int_t iv=0; iv<100; iv++){
    printf("%d\n",iv);
    orgVertex[0]=gRandom->Gaus()*0.01;
    orgVertex[1]=gRandom->Gaus()*0.01;
    orgVertex[2]=gRandom->Gaus()*3;
    corrC->FastSimDistortedVertex(orgVertex,100, aV,aVO,cV,cVO,pcstream);
    corrS->FastSimDistortedVertex(orgVertex,100, aV,aVO,cV,cVO,pcstreamS);
    corrT->FastSimDistortedVertex(orgVertex,100, aV,aVO,cV,cVO,pcstreamT);
  }
  delete pcstream;
  delete pcstreamS;
  delete pcstreamT;
  TFile fC("vertexDistortC.root");
  TFile fT("vertexDistortT.root");
  TFile fS("vertexDistortS.root");
  TTree * treeT=(TTree*)fT.Get("vertex"); // twist distrotions
  TTree * treeS=(TTree*)fS.Get("vertex"); // shape distortions
  TCanvas *canvasD=new TCanvas("canvasD","canvasD");
  canvasD->Divide(2,2);
  canvasD->cd(1);
  treeT->SetLineColor(2);
  treeT->Draw("av.fPosition[0]-x>>hisATX(100,-0.2,0.2)","","");
  treeT->SetLineColor(4);
  treeT->Draw("cv.fPosition[0]-x>>hisCTX(100,-0.2,0.2)","","same");
  canvasD->cd(2);
  treeT->SetLineColor(2);
  treeT->Draw("av.fPosition[1]-y>>hisATY(100,-0.2,0.2)","","");
  treeT->SetLineColor(4);
  treeT->Draw("cv.fPosition[1]-y>>hisCTY(100,-0.2,0.2)","","same");
  //
  canvasD->cd(3);
  treeS->SetLineColor(2);
  treeS->Draw("av.fPosition[0]-x>>hisASX(100,-0.2,0.2)","","");
  treeS->SetLineColor(4);
  treeS->Draw("cv.fPosition[0]-x>>hisCSX(100,-0.2,0.2)","","same");
  canvasD->cd(4);
  treeS->SetLineColor(2);
  treeS->Draw("av.fPosition[1]-y>>hisASY(100,-0.2,0.2)","","");
  treeS->SetLineColor(4);
  treeS->Draw("cv.fPosition[1]-y>>hisCSY(100,-0.2,0.2)","","same");
  canvasD->SaveAs("vertexShift.ps");

}
