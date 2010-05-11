//
// Make lookup of distortion TPC distortion +(ITS and TRD) 
// Input: Residual histograms obtained in the AliTPCcalibTime

// Residual histograms:  
//   1. TPC-ITS  - entrance of the TPC  
//   2. TPC-ITS  - at the vertex
//   3. TPC-TRD  - outer wall of the TPC 

// Histogram binning:
//   1. Theta    - fP3
//   2. Phi      - global phi at the entrance (case 1,2) and at the outer wall of TPC (case 3)
//   3. snp(phi) - fP2 - local inclination angle at reference X 
              
// Output value:
//   Mean residuals, rms and number of entries in each bing
// 
// 

 
/* 
  Example usage:
 
  gSystem->AddIncludePath("-I$ALICE_ROOT/STAT")
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC")
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros")
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  .L  $ALICE_ROOT/TPC/CalibMacros/MakeLookup.C+
  Int_t run=11111
  MakeLookup(run,0);

  //Local check of procedure
  TTreeSRedirector * pcstream  = new TTreeSRedirector("mean.root");
  TFile f("CalibObjects.root");
  AliTPCcalibTime  *calibTime= f.Get("calibTime");
  THnSparse * his = calibTime->GetResHistoTPCITS(0);
  his->GetAxis(1)->SetRangeUser(-1,1);
  his->GetAxis(2)->SetRange(0,1000000);
  his->GetAxis(3)->SetRangeUser(-0.3,0.3);
  //

  

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
#include "AliTrackerBase.h"
#endif

TChain * chain=0;
void MakeLookup(Int_t run, Int_t mode=0);
//void MakeLookupHisto(THnSparse * his, TTreeSRedirector *pcstream, const char* hname, Int_t run);
void FitLookup(TChain *chainIn, const char *prefix, TVectorD &vecA, TVectorD &vecC, TVectorD& vecStatA, TVectorD &vecStatD, TCut cut, TObjArray *picArray);
void MakeFits(Int_t run);
void MakeFitTree();
void DrawDistortionDy(TCut cutUser="", Double_t ymin=-0.6, Double_t ymax=0.6);
void DrawDistortionMaps(const char *fname="mean.root");


void MakeLookup(Int_t run, Int_t mode){
  //
  // make a lookup tree with mean values
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/STAT");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  if (mode==1){
    DrawDistortionDy("abs(snp)<0.25"); 
    DrawDistortionMaps();
    return;
  }
  if (mode==2) return  MakeFitTree();

  TFile f("CalibObjects.root");
  AliTPCcalibTime  *calibTime= (AliTPCcalibTime*)f.Get("calibTime");
  //
  TTreeSRedirector * pcstream  = new TTreeSRedirector("mean.root");
  THnSparse * his = 0;
  const char * hname[5]={"dy","dz","dsnp","dtheta","d1pt"};
  if (calibTime){
    for (Int_t ihis=0; ihis<5;ihis++){
      his = calibTime->GetResHistoTPCITS(ihis);
      his->GetAxis(1)->SetRangeUser(-1.1,1.1);
      his->GetAxis(2)->SetRange(0,1000000);
      his->GetAxis(3)->SetRangeUser(-0.5,0.5);
      AliTPCCorrection::MakeDistortionMap(his,pcstream, Form("ITS%s",hname[ihis]),run);
      //
      his = calibTime->GetResHistoTPCvertex(ihis);
      his->GetAxis(1)->SetRangeUser(-1.1,1.1);
      his->GetAxis(2)->SetRange(0,1000000);
      his->GetAxis(3)->SetRangeUser(-0.5,0.5);
      AliTPCCorrection::MakeDistortionMap(his,pcstream, Form("Vertex%s",hname[ihis]),run);
      //
      his = calibTime->GetResHistoTPCTRD(ihis);
      his->GetAxis(1)->SetRangeUser(-1.1,1.1);
      his->GetAxis(2)->SetRange(0,1000000);
      his->GetAxis(3)->SetRangeUser(-0.5,0.5);
      AliTPCCorrection::MakeDistortionMap(his,pcstream, Form("TRD%s",hname[ihis]),run);
    }
  }
  delete pcstream;
  gSystem->Exec("echo `pwd`/mean.root > distort.txt");
  MakeFits(run);
}


void MakeFits(Int_t run){
  //
  // Make the fits of distortion
  //   store fit results and QA pictures in the file distortFit.root 
  //
  TCut cut="entries>50&&rms>0";
  TTreeSRedirector *pcstream = new TTreeSRedirector("distortFit.root");
  AliXRDPROOFtoolkit tool;
  TVectorD vecA[100],vecC[100], vecStatA[100],vecStatC[100];
  TObjArray * picArray= new TObjArray();
  //
  TChain *chainITSdy      = tool.MakeChain("distort.txt","ITSdy",0,100000);
  TChain *chainITSdsnp     = tool.MakeChain("distort.txt","ITSdsnp",0,100000);
  TChain *chainTRDdy      = tool.MakeChain("distort.txt","TRDdy",0,100000);
  TChain *chainTRDdsnp     = tool.MakeChain("distort.txt","TRDdsnp",0,100000);
  TChain *chainVertexdy   = tool.MakeChain("distort.txt","Vertexdy",0,100000);
  TChain *chainVertexdsnp  = tool.MakeChain("distort.txt","Vertexdsnp",0,100000);

  (*pcstream)<<"fits"<<"run="<<run;

  FitLookup(chainITSdy, "yits-tpc",vecA[0],vecC[0],vecStatA[0],vecStatC[0],cut, picArray);
  (*pcstream)<<"fits"<<"itsdyA.="<<&vecA[0]<<"itsdyC.="<<&vecC[0];
  (*pcstream)<<"fits"<<"itsdyAS.="<<&vecStatA[0]<<"itsdyCS.="<<&vecStatC[0];
  FitLookup(chainITSdsnp, "snpits-tpc",vecA[1],vecC[1],vecStatA[1],vecStatC[1],cut, picArray);
  (*pcstream)<<"fits"<<"itsdsnpA.="<<&vecA[1]<<"itsdsnpC.="<<&vecC[1];
  (*pcstream)<<"fits"<<"itsdsnpAS.="<<&vecStatA[1]<<"itsdsnpCS.="<<&vecStatC[1];
  //
  FitLookup(chainTRDdy, "ytrd-tpc",vecA[5],vecC[5],vecStatA[5],vecStatC[5],cut, picArray);
  (*pcstream)<<"fits"<<"trddyA.="<<&vecA[5]<<"trddyC="<<&vecC[5];
  (*pcstream)<<"fits"<<"trddyAS.="<<&vecStatA[5]<<"trddyCS.="<<&vecStatC[5];
  FitLookup(chainTRDdsnp, "snptrd-tpc",vecA[6],vecC[6],vecStatA[6],vecStatC[6],cut, picArray);
  (*pcstream)<<"fits"<<"trddsnpA.="<<&vecA[6]<<"trddsnpC.="<<&vecC[6];
  (*pcstream)<<"fits"<<"trddsnpAS.="<<&vecStatA[6]<<"trddsnpCS.="<<&vecStatC[6];
  //
  FitLookup(chainVertexdy, "yvertex-tpc",vecA[10],vecC[10],vecStatA[10],vecStatC[10],cut, picArray);
  (*pcstream)<<"fits"<<"vertexdyA.="<<&vecA[10]<<"vertexdyC.="<<&vecC[10];
  (*pcstream)<<"fits"<<"vertexdyAS.="<<&vecStatA[10]<<"vertexdyCS.="<<&vecStatC[10];
  FitLookup(chainVertexdsnp, "snpvertex-tpc",vecA[11],vecC[11],vecStatA[11],vecStatC[11],cut, picArray);
  (*pcstream)<<"fits"<<"vertexdsnpA.="<<&vecA[11]<<"vertexdsnpC.="<<&vecC[11];
  (*pcstream)<<"fits"<<"vertexdsnpAS.="<<&vecStatA[11]<<"vertexdsnpCS.="<<&vecStatC[11];
  //
  (*pcstream)<<"fits"<<"\n";

  pcstream->GetFile()->cd();
  for (Int_t i=0;i<picArray->GetEntries(); i++){
    TObject * obj = picArray->At(i);
    if (obj) obj->Write(obj->GetName());
  }
  //
  //
  //
  chainITSdy->AddFriend(chainITSdy,"ITSdy");
  chainITSdy->AddFriend(chainITSdsnp,"ITSdsnp");
  chainITSdy->AddFriend(chainTRDdy,"TRDdy");
  chainITSdy->AddFriend(chainTRDdsnp,"TRDdsnp");
  chainITSdy->AddFriend(chainVertexdy,"Vertexdy");
  chainITSdy->AddFriend(chainVertexdsnp,"Vertexdsnp");
  TTree * tree = chainITSdy->CloneTree();
  tree->Write("distortionTree");
  delete pcstream;

  //
  //
}





void FitLookup(TChain *chainIn, const char *prefix, TVectorD &vecA, TVectorD &vecC, TVectorD& vecStatA, TVectorD &vecStatC,  TCut cut, TObjArray *picArray){ 
  //  TCut cut="entries>100&&rms>0";
  vecStatA.ResizeTo(6);
  vecStatC.ResizeTo(6);
  vecA.ResizeTo(10);
  vecC.ResizeTo(10);
  Int_t  npointsMax=30000000;
  TStatToolkit toolkit;
  Double_t chi2A=0;
  Double_t chi2C=0;
  Int_t    npointsA=0;
  Int_t    npointsC=0;
  TVectorD paramA;
  TVectorD paramC;
  TMatrixD covar;
  TString *strFitYA;
  TString *strFitYC;

  TString  fstring="";   // magnetic part
  //fstring+="(1)++";             //0    - offset value
  fstring+="(cos(phi))++";        //1    - cos part 
  fstring+="(sin(phi))++";        //2    - sin part
  //
  fstring+="(theta)++";           //3
  fstring+="(theta*cos(phi))++";  //4
  fstring+="(theta*sin(phi))++";  //5
  //
  fstring+="(snp)++";             //6   - delta X(radial)  coeficients - offset
  fstring+="(snp*cos(phi))++";    //7   - delta X(radial) coeficient  
  fstring+="(snp*sin(phi))++";    //8   
  //
  //
  strFitYA = TStatToolkit::FitPlane(chainIn,"mean", fstring.Data(),cut+"theta>0", chi2A,npointsA,paramA,covar,-1,0, npointsMax, kFALSE);
  strFitYC = TStatToolkit::FitPlane(chainIn,"mean", fstring.Data(),cut+"theta<0", chi2C,npointsC,paramC,covar,-1,0, npointsMax,kFALSE);
  strFitYA->Tokenize("++")->Print();
  strFitYC->Tokenize("++")->Print();
  chainIn->SetAlias("dyA",strFitYA->Data());
  chainIn->SetAlias("dyC",strFitYC->Data());
  //
  TH1 * his=0;
  vecA.ResizeTo(paramA.GetNrows());
  vecC.ResizeTo(paramC.GetNrows());
  vecA=paramA;
  vecC=paramC;
  // A side stat
  vecStatA[0]=npointsA;
  vecStatA[1]=TMath::Sqrt(chi2A/npointsA);
  chainIn->Draw("mean",cut+"theta>0");
  his=chainIn->GetHistogram();
  his->SetName(Form("Orig #Delta%sA",prefix));
  his->SetTitle(Form("Orig #Delta%sA",prefix));
  picArray->AddLast(his->Clone());
  vecStatA[2] = his->GetMean();
  vecStatA[3] = his->GetRMS();
  chainIn->Draw("mean-dyA",cut+"theta>0");
  his=chainIn->GetHistogram();
  his->SetName(Form("Corr. #Delta%sA",prefix));
  his->SetTitle(Form("Corr. #Delta%sA",prefix));
  picArray->AddLast(his->Clone());
  vecStatA[4] = his->GetMean();
  vecStatA[5] = his->GetRMS();
  //
  // C side stat
  vecStatC[0]=npointsC;
  vecStatC[1]=TMath::Sqrt(chi2C/npointsC);
  chainIn->Draw("mean",cut+"theta<0");
  his=chainIn->GetHistogram();
  his->SetName(Form("Orig #Delta%sC",prefix));
  his->SetTitle(Form("Orig #Delta%sC",prefix));
  picArray->AddLast(his->Clone());
  vecStatC[2] = his->GetMean();
  vecStatC[3] = his->GetRMS();
  chainIn->Draw("mean-dyC",cut+"theta<0");
  his=chainIn->GetHistogram();
  his->SetName(Form("Corr. #Delta%sC",prefix));
  his->SetTitle(Form("Corr. #Delta%sC",prefix));
  picArray->AddLast(his->Clone());
  vecStatC[4] = his->GetMean();
  vecStatC[5] = his->GetRMS();
  
  vecStatA.Print();
  vecStatC.Print();
}



void DrawDistortionDy(TCut cutUser, Double_t ymin, Double_t ymax){
  //
  //
  //
  TFile fplus("meanBplus.root");
  TFile fminus("meanBminus.root");
  TTree * titsDyPlus=  (TTree*)fplus.Get("ITSdy");
  TTree * titsDyMinus= (TTree*)fminus.Get("ITSdy");
  TTree * ttrdDyPlus=  (TTree*)fplus.Get("TRDdy");
  TTree * ttrdDyMinus= (TTree*)fminus.Get("TRDdy");
  //
  TTree * titsDsnpPlus=  (TTree*)fplus.Get("ITSdsnp");
  TTree * titsDsnpMinus= (TTree*)fminus.Get("ITSdsnp");
  TTree * ttrdDsnpPlus=  (TTree*)fplus.Get("TRDdsnp");
  TTree * ttrdDsnpMinus= (TTree*)fminus.Get("TRDdsnp");
  TTree * titsDthetaPlus=  (TTree*)fplus.Get("ITSdtheta");
  TTree * titsDthetaMinus= (TTree*)fminus.Get("ITSdtheta");
  TTree * ttrdDthetaPlus=  (TTree*)fplus.Get("TRDdtheta");
  TTree * ttrdDthetaMinus= (TTree*)fminus.Get("TRDdtheta");
  TCut cut="rms>0&&entries>100&&abs(theta)<1&&abs(snp)<0.3";
  cut+=cutUser;
  //
  titsDyPlus->AddFriend(titsDyMinus,"TM.");
  titsDyPlus->AddFriend(titsDsnpMinus,"Tsnp.");
  titsDyPlus->AddFriend(titsDthetaMinus,"Ttheta.");
  ttrdDyPlus->AddFriend(ttrdDyMinus,"TM.");
  titsDsnpPlus->AddFriend(titsDsnpMinus,"TM.");
  ttrdDsnpPlus->AddFriend(ttrdDsnpMinus,"TM.");
 
  //
  //
  //
  Int_t entries=0;
  TGraph * graphITSYC[3];
  TGraph * graphITSYA[3];
  TGraph * graphTRDYC[3];
  TGraph * graphTRDYA[3];
  //
  //
  entries=titsDyPlus->Draw("mean-TM..mean:phi",cut+"theta<-0.1","goff");
  graphITSYC[0] = new TGraph(entries, titsDyPlus->GetV2(),titsDyPlus->GetV1());
  titsDyMinus->Draw("mean:phi",cut+"theta<-0.1","goff");
  graphITSYC[2] = new TGraph(entries, titsDyMinus->GetV2(),titsDyMinus->GetV1());
  titsDyPlus->Draw("mean:phi",cut+"theta<-0.1","goff");
  graphITSYC[1] = new TGraph(entries, titsDyPlus->GetV2(),titsDyPlus->GetV1());
  //
  entries=titsDyPlus->Draw("mean-TM..mean:phi",cut+"theta>0.1","goff");
  graphITSYA[0] = new TGraph(entries, titsDyPlus->GetV2(),titsDyPlus->GetV1());
  titsDyMinus->Draw("mean:phi",cut+"theta>0.1","goff");
  graphITSYA[2] = new TGraph(entries, titsDyMinus->GetV2(),titsDyMinus->GetV1());
  titsDyPlus->Draw("mean:phi",cut+"theta>0.1","goff");
  graphITSYA[1] = new TGraph(entries, titsDyPlus->GetV2(),titsDyPlus->GetV1());
  //
  entries=ttrdDyPlus->Draw("mean-TM..mean:phi",cut+"theta<-0.1&&TM..entries>100","goff");
  graphTRDYC[0] = new TGraph(entries, ttrdDyPlus->GetV2(),ttrdDyPlus->GetV1());
  ttrdDyMinus->Draw("mean:phi",cut+"theta<-0.1","goff");
  graphTRDYC[2] = new TGraph(entries, ttrdDyMinus->GetV2(),ttrdDyMinus->GetV1());
  ttrdDyPlus->Draw("mean:phi",cut+"theta<-0.1","goff");
  graphTRDYC[1] = new TGraph(entries, ttrdDyPlus->GetV2(),ttrdDyPlus->GetV1());
  //
  entries=ttrdDyPlus->Draw("mean-TM..mean:phi",cut+"theta>0.1&&TM..entries>100","goff");
  graphTRDYA[0] = new TGraph(entries, ttrdDyPlus->GetV2(),ttrdDyPlus->GetV1());
  ttrdDyMinus->Draw("mean:phi",cut+"theta>0.1","goff");
  graphTRDYA[2] = new TGraph(entries, ttrdDyMinus->GetV2(),ttrdDyMinus->GetV1());
  ttrdDyPlus->Draw("mean:phi",cut+"theta>0.1","goff");
  graphTRDYA[1] = new TGraph(entries, ttrdDyPlus->GetV2(),ttrdDyPlus->GetV1());
  //


  //
  {for (Int_t i=0; i<3; i++){
    graphITSYC[i]->SetMaximum(ymax+0.2); graphITSYC[i]->SetMinimum(ymin);
    graphITSYC[i]->GetXaxis()->SetTitle("#phi"); graphITSYC[i]->GetYaxis()->SetTitle("#Delta r#phi (cm)");
    graphITSYC[i]->SetMarkerColor(i+1); graphITSYC[i]->SetMarkerStyle(20+i);
    graphITSYA[i]->SetMaximum(ymax+0.2); graphITSYA[i]->SetMinimum(ymin);
    graphITSYA[i]->GetXaxis()->SetTitle("#phi"); graphITSYA[i]->GetYaxis()->SetTitle("#Delta r#phi (cm)");
    graphITSYA[i]->SetMarkerColor(i+1); graphITSYA[i]->SetMarkerStyle(20+i);
    graphTRDYC[i]->SetMaximum(ymax+0.2); graphTRDYC[i]->SetMinimum(ymin);
    graphTRDYC[i]->GetXaxis()->SetTitle("#phi"); graphTRDYC[i]->GetYaxis()->SetTitle("#Delta r#phi (cm)");
    graphTRDYC[i]->SetMarkerColor(i+1); graphTRDYC[i]->SetMarkerStyle(20+i);
    graphTRDYA[i]->SetMaximum(ymax+0.2); graphTRDYA[i]->SetMinimum(ymin);
    graphTRDYA[i]->GetXaxis()->SetTitle("#phi"); graphTRDYA[i]->GetYaxis()->SetTitle("#Delta r#phi (cm)");
    graphTRDYA[i]->SetMarkerColor(i+1); graphTRDYA[i]->SetMarkerStyle(20+i);
    }
  }
  TLegend *legend=0;
  TString cname="cdeltaY";
  TString dirName=gSystem->GetFromPipe("dirname `pwd` | xargs basename ");
  TString splus=gSystem->GetFromPipe("ls -al meanBplus.root | gawk '{print \"File: \"$8\" \"$6\" \"$7\" \";}'");
  TString sminus=gSystem->GetFromPipe("ls -al meanBminus.root | gawk '{print \"File: \"$8\" \"$6\" \"$7\" \";}'");
  
  TCanvas *cdeltaY = new TCanvas("cdeltaY","cdeltaY",1300,800);
  cdeltaY->Divide(2,2);

  cdeltaY->cd(1);
  graphITSYC[0]->Draw("ap");
  graphITSYC[1]->Draw("p");
  graphITSYC[2]->Draw("p");
  legend = new TLegend(0.6,0.6,1.0,1.0,"ITS-TPC C side #Delta r-#phi");
  legend->AddEntry("",dirName.Data());
  legend->AddEntry("",splus.Data());
  legend->AddEntry("",sminus.Data());
  legend->AddEntry(graphITSYC[0],"B_{plus}-B_{minus}");
  legend->AddEntry(graphITSYC[1],"B_{plus}");
  legend->AddEntry(graphITSYC[2],"B_{minus}");
  legend->Draw();

  cdeltaY->cd(2);
  graphITSYA[0]->Draw("ap");
  graphITSYA[1]->Draw("p");
  graphITSYA[2]->Draw("p");
  legend = new TLegend(0.6,0.6,1.0,1.0,"ITS-TPC A side #Delta r-#phi");
  legend->AddEntry(graphITSYA[0],"B_{plus}-B_{minus}");
  legend->AddEntry(graphITSYA[1],"B_{plus}");
  legend->AddEntry(graphITSYA[2],"B_{minus}");
  legend->Draw();

  //
  cdeltaY->cd(3);
  graphTRDYC[0]->Draw("ap");
  graphTRDYC[1]->Draw("p");
  graphTRDYC[2]->Draw("p");
  legend = new TLegend(0.6,0.6,1.0,1.0,"TRD-TPC C side #Delta r-#phi");
  legend->AddEntry(graphTRDYC[0],"B_{plus}-B_{minus}");
  legend->AddEntry(graphTRDYC[1],"B_{plus}");
  legend->AddEntry(graphTRDYC[2],"B_{minus}");
  legend->Draw();

  cdeltaY->cd(4);
  graphTRDYA[0]->Draw("ap");
  graphTRDYA[1]->Draw("p");
  graphTRDYA[2]->Draw("p");
  legend = new TLegend(0.6,0.6,1.0,1.0,"TRD-TPC A side #Delta r-#phi");
  legend->AddEntry(graphTRDYA[0],"B_{plus}-B_{minus}");
  legend->AddEntry(graphTRDYA[1],"B_{plus}");
  legend->AddEntry(graphTRDYA[2],"B_{minus}");
  legend->Draw();
  cdeltaY->SaveAs(dirName+"_deltaY.pdf");

}


void DrawDistortionMaps(const char *fname){
  //
  //
  //
  TFile f(fname);
  TTree *ITSdy=(TTree*)f.Get("ITSdy");
  TTree *ITSdsnp=(TTree*)f.Get("ITSdsnp");
  TTree *ITSd1pt=(TTree*)f.Get("ITSd1pt");
  TCanvas *cdist = new TCanvas("distITS","distITS",1200,800); 
  cdist->Divide(3,2);
  cdist->Draw("");
  TH1 * his=0;

  cdist->cd(1);
  ITSdy->Draw("10*mean","rms>0&&abs(snp)<0.25&&entries>100","");
  his = (TH1*)ITSdy->GetHistogram()->Clone();
  his->SetName("dY"); his->SetXTitle("#Delta_{r#phi} (cm)");
  his->Draw();
  //
  cdist->cd(2);
  ITSdsnp->Draw("1000*mean","rms>0&&abs(snp)<0.25&&entries>200","");
  his = (TH1*)ITSdsnp->GetHistogram()->Clone();
  his->SetName("dsnp"); his->SetXTitle("#Delta_{sin(#phi)}");
  his->Draw();
  //
  cdist->cd(3);
  ITSd1pt->Draw("mean","rms>0&&abs(snp)<0.15&&entries>400","");
  his = (TH1*)ITSd1pt->GetHistogram()->Clone();
  his->SetName("d1pt"); his->SetXTitle("#Delta_{1/pt} (1/GeV)");
  his->Draw();
  //
  ITSdy->SetMarkerSize(0.3);
  ITSdsnp->SetMarkerSize(0.3);
  ITSd1pt->SetMarkerSize(0.3);
  {for (Int_t type=0; type<3; type++){
      cdist->cd(4+type);
      TTree * tree =ITSdy;
      if (type==1) tree=ITSdsnp;
      if (type==2) tree=ITSd1pt;
      Int_t counter=0;
      for (Double_t theta=-1; theta<=1; theta+=0.2){
	TCut cut=Form("rms>0&&abs(snp)<0.25&&entries>20&&abs(theta-%f)<0.05",theta);
	tree->SetMarkerStyle(20+counter);
	Option_t *option=(counter==0)? "prof": "profsame";
	if (theta>0) tree->SetMarkerColor(2);
	if (theta<0) tree->SetMarkerColor(4);
	if (type==0) tree->Draw("10*mean:phi>>his(45,-pi,pi)",cut,option);      
	if (type==1) tree->Draw("1000*mean:phi>>his(45,-pi,pi)",cut,option);      
	if (type==2) tree->Draw("mean:phi>>his(45,-pi,pi)",cut,option);      
	his = (TH1*)tree->GetHistogram()->Clone();
	his->SetName(Form("%d_%d",counter,type));
	his->SetXTitle("#phi");
	his->SetMaximum(4);
	his->SetMinimum(-4);
	if (type==2){
	  his->SetMaximum(0.06);
	  his->SetMinimum(-0.06);	  
	}
	his->Draw(option);
	counter++;
      }
    }
  }
  cdist->SaveAs("itstpcdistrotion.pdf");
}


void MakeFitTree(){
  //
  // 1. Initialize ocdb e.g 
  //     .x $ALICE_ROOT/ANALYSIS/CalibMacros/Pass0/ConfigCalibTrain.C(114972,)  
  //     .x $ALICE_ROOT/ANALYSIS/CalibMacros/Pass0/ConfigCalibTrain.C(114972,"local:///lustre/alice/alien/alice/data/2010/OCDB")
  //
  //
  Double_t bzField=AliTrackerBase::GetBz();  
  Double_t vdrift = 2.6; // [cm/us]   // to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 
  Double_t T1 = 0.9;
  Double_t T2 = 1.5;
  //
  AliTPCExBTwist *twistX1= new AliTPCExBTwist;
  twistX1->SetXTwist(0.001);  // 1 mrad twist in x
  twistX1->SetName("tX1");
  //
  AliTPCExBTwist *twistY1= new AliTPCExBTwist;
  twistY1->SetYTwist(0.001);   // 1 mrad twist in Y
  twistY1->SetName("tY1");
  //
  AliTPCGGVoltError* ggErrorA = new AliTPCGGVoltError;
  ggErrorA->SetDeltaVGGA(40.);  // delta GG - 1 mm
  ggErrorA->SetName("ggA");
  //
  AliTPCGGVoltError* ggErrorC = new AliTPCGGVoltError;
  ggErrorC->SetDeltaVGGC(40.);  // delta GG - 1 mm
  ggErrorC->SetName("ggC");
  // conical free param
  //
  AliTPCExBConical *exbConA = new AliTPCExBConical;
  Float_t conicalAA[3]={1., 0.0,0.00};
  Float_t conicalCA[3]={0., 0.0,0.00};
  exbConA->SetConicalA(conicalAA);
  exbConA->SetConicalC(conicalCA);
  exbConA->SetConicalFactor(0.25);
  exbConA->SetName("ExBConA");          // conical factor - A side

  //conical free param
  AliTPCExBConical *exbConC = new AliTPCExBConical;
  Float_t conicalAC[3]={0., 0.0,0.00};
  Float_t conicalCC[3]={1., 0.0,0.00};
  exbConC->SetConicalA(conicalAC);
  exbConC->SetConicalC(conicalCC);
  exbConC->SetConicalFactor(0.25);
  exbConC->SetName("ExBConC");          // conical factor - c side
  // conical as measured
  AliTPCExBConical *exbCon = new AliTPCExBConical;
  Float_t conicalC[3]={1., 0.6,-0.08};
  Float_t conicalA[3]={0.9, 0.3,0.04};
  exbCon->SetConicalA(conicalA);
  exbCon->SetConicalC(conicalC);
  exbCon->SetConicalFactor(0.25);
  exbCon->SetName("ExBCon");          // conical factor


  TObjArray * corr = new TObjArray;
  corr->AddLast(twistX1);
  corr->AddLast(twistY1);
  corr->AddLast(ggErrorA);
  corr->AddLast(ggErrorC);
  corr->AddLast(exbConA);
  corr->AddLast(exbConC);
  corr->AddLast(exbCon);
  //
  AliTPCComposedCorrection *cc= new AliTPCComposedCorrection ;
  cc->SetCorrections(corr);
  cc->SetOmegaTauT1T2(wt,T1,T2);
  cc->Init();
  //  cc->SetMode(1);
  cc->Print("DA"); // Print used correction classes
  cc->SetName("Comp");
  //corr->AddLast(cc);
  //
  TFile f("mean.root");
  TTree * tree= 0;
  tree = (TTree*)f.Get("ITSdy");
  AliTPCCorrection::MakeTrackDistortionTree(tree,0,0,corr,7);
  tree = (TTree*)f.Get("TRDdy");
  AliTPCCorrection::MakeTrackDistortionTree(tree,1,0,corr,7);
  tree = (TTree*)f.Get("Vertexdy");
  AliTPCCorrection::MakeTrackDistortionTree(tree,2,0,corr,7);
  //
  tree = (TTree*)f.Get("ITSdsnp");
  AliTPCCorrection::MakeTrackDistortionTree(tree,0,2,corr,7);
  tree = (TTree*)f.Get("TRDdsnp");
  AliTPCCorrection::MakeTrackDistortionTree(tree,1,2,corr,7);
  tree = (TTree*)f.Get("Vertexdsnp");
  AliTPCCorrection::MakeTrackDistortionTree(tree,2,2,corr,7);
  //
  tree = (TTree*)f.Get("ITSdz");
  AliTPCCorrection::MakeTrackDistortionTree(tree,0,1,corr,7);
  //tree = (TTree*)f.Get("TRDdz");
  //AliTPCCorrection::MakeTrackDistortionTree(tree,1,1,corr,3);
  tree = (TTree*)f.Get("Vertexdz");
  AliTPCCorrection::MakeTrackDistortionTree(tree,2,1,corr,7);

  //
  tree = (TTree*)f.Get("ITSdtheta");
  AliTPCCorrection::MakeTrackDistortionTree(tree,0,3,corr,7);
  //tree = (TTree*)f.Get("TRDdtheta");
  //  AliTPCCorrection::MakeTrackDistortionTree(tree,1,3,corr,3);
  //
  tree = (TTree*)f.Get("Vertexdtheta");
  AliTPCCorrection::MakeTrackDistortionTree(tree,2,3,corr,7);

  tree = (TTree*)f.Get("ITSd1pt");
  AliTPCCorrection::MakeTrackDistortionTree(tree,0,4,corr,7);
  tree = (TTree*)f.Get("TRDd1pt");
  AliTPCCorrection::MakeTrackDistortionTree(tree,1,4,corr,7);
  tree = (TTree*)f.Get("Vertexd1pt");
  AliTPCCorrection::MakeTrackDistortionTree(tree,2,4,corr,7);


}



void MakeGlobalFit(){
  AliXRDPROOFtoolkit tool;
  TChain *chain      = tool.MakeChain("distortion.txt","fit",0,100000);

  TCut cutS="rms>0&&entries>100";
  TCut cut=cutS+"(ptype==0||ptype==2||ptype==3||(ptype==4&&dtype==0))";
  Int_t  npointsMax=30000000;
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD param;
  TMatrixD covar;
  chain->SetAlias("side","(1+(theta>0)*2)");
  TString  fstring="";     // magnetic part
  fstring+="tX1++";        //1    - twist x in mrad
  fstring+="tY1++";        //2    - twist y in mrad 
  fstring+="ExBCon++";     //3    - custom
  fstring+="ExBCon*cos(phi)++";     //3    - custom
  fstring+="ExBCon*sin(phi)++";     //3    - custom
  fstring+="ggA++";        //4    - gating grid
  fstring+="ggA*cos(phi)++";        //4    - gating grid
  fstring+="ggA*sin(phi)++";        //4    - gating grid
  fstring+="ggC++";        //5    - gating grid
  fstring+="ggC*cos(phi)++";        //5    - gating grid
  fstring+="ggC*sin(phi)++";        //5    - gating grid
  fstring+="(ptype==2)++";
  fstring+="(ptype==3)++";
  fstring+="(ptype==4)++";
  fstring+="(ptype==4)*bz++";

  TString *strDelta = TStatToolkit::FitPlane(chain,"mean:rms", fstring.Data(),cut+"(Entry$%13==0)", chi2,npoints,param,covar,-1,0, npointsMax, kTRUE);
  chain->SetAlias("delta",strDelta->Data());
  strDelta->Tokenize("++")->Print();

}


void MakeGlobalFitRelative(Int_t highFrequency=0){
  //
  // Make a global fit of ExB
  // To get rid of the misalignment errors -
  //   Use relative change of deltas for 2 different filed settings
  //    

  AliXRDPROOFtoolkit tool;
  TChain *chain         = tool.MakeChain("distortion.txt","fit",0,100000);
  TChain *chainPlus     = tool.MakeChain("distortionPlus.txt","fit",0,100000);
  TChain *chainMinus    = tool.MakeChain("distortionMinus.txt","fit",0,100000);
  chainPlus->AddFriend(chainMinus,"M");

  TCut cut="rms>0&&M.rms>0&&entries>100&&dtype==0&&(ptype==0||ptype==2||ptype==4)";
  chainPlus->SetAlias("deltaM","mean-M.mean");
  Int_t  npointsMax=30000000;
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD param;
  TMatrixD covar;
  //
  TString  fstring="";     // fit string
  fstring+="(ptype==0)*cos(phi)++";   //     - primitive gx movement
  fstring+="(ptype==0)*sin(phi)++";   //     - primitive gy movement

  fstring+="(tX1-M.tX1)++";           // 1    - twist x in mrad
  fstring+="(tY1-M.tY1)++";           // 2    - twist y in mrad 
  //
  fstring+="(ExBConA-M.ExBConA)++";   // 3 conical shape A side 
  fstring+="(ExBConC-M.ExBConC)++";   // 4 conical shape C side 
  //
  fstring+="(ggA-M.ggA)++";   // 3 conical shape A side 
  fstring+="(ggC-M.ggC)++";   // 4 conical shape C side 
  
  if (highFrequency) {for (Int_t i=1; i<highFrequency; i++){
      fstring+=Form("((ggA-M.ggA)*cos(%d*phi))++",i);
      fstring+=Form("((ggA-M.ggA)*sin(%d*phi))++",i);
      //
      fstring+=Form("((ggC-M.ggC)*cos(%d*phi))++",i);
      fstring+=Form("((ggC-M.ggC)*sin(%d*phi))++",i);
      if (i>1){
	fstring+=Form("((ExBConA-M.ExBConA)*cos(%d*phi))++",i);
	fstring+=Form("((ExBConA-M.ExBConA)*sin(%d*phi))++",i);
	fstring+=Form("((ExBConC-M.ExBConC)*cos(%d*phi))++",i);
	fstring+=Form("((ExBConC-M.ExBConC)*sin(%d*phi))++",i);
      }
    }}

  TString *strDelta = TStatToolkit::FitPlane(chainPlus,"deltaM:rms", fstring.Data(),cut+"(Entry$%7==0)", chi2,npoints,param,covar,-1,0, npointsMax, kTRUE);
  chainPlus->SetAlias("delta",strDelta->Data());
  TObjArray *fitArray = strDelta->Tokenize("++");
  fitArray->Print();

  //
  // Draw results
  //  
  TH1::AddDirectory(0);
  TLatex *lfit=new TLatex;


  TCanvas * canvasdY= new TCanvas(Form("deltaY%d",highFrequency),Form("deltaY%d",highFrequency),1200,800);
  TCanvas * canvasdSnp= new TCanvas(Form("deltaSnp%d",highFrequency),Form("deltaSnp%d",highFrequency),1200,800);
  TCanvas * canvasd1pt= new TCanvas(Form("delta1pt%d",highFrequency),Form("delta1pt%d",highFrequency),1200,800);

  canvasdY->Divide(3,2);
  chainPlus->SetMarkerStyle(25);
  chainPlus->SetMarkerSize(0.5);
  chainPlus->SetMarkerColor(1);
  canvasdY->cd(1);
  chainPlus->Draw("deltaM:delta",cut+"theta>0&&abs(snp)<0.2&&ptype==0","");
  canvasdY->cd(2);
  chainPlus->SetMarkerColor(1);
  chainPlus->Draw("deltaM-delta>>deltaCorr(50,-0.4,0.4)",cut+"theta>0&&abs(snp)<0.2&&ptype==0","err");
  chainPlus->SetMarkerColor(3);
  chainPlus->Draw("deltaM",cut+"theta>0&&abs(snp)<0.2&&ptype==0","errsame");
  canvasdY->cd(3);
  chainPlus->SetMarkerColor(1);
  chainPlus->Draw("deltaM:phi",cut+"theta>0&&abs(snp)<0.2&&ptype==0","");
  chainPlus->SetMarkerColor(2);
  chainPlus->Draw("deltaM-delta:phi",cut+"theta>0&&abs(snp)<0.2&&ptype==0","same");
  canvasdY->cd(6);
  chainPlus->SetMarkerStyle(26);
  chainPlus->SetMarkerColor(1);  
  chainPlus->Draw("deltaM:phi",cut+"theta<0&&abs(snp)<0.2&&ptype==0","");
  chainPlus->SetMarkerColor(4);  
  chainPlus->Draw("deltaM-delta:phi",cut+"theta<0&&abs(snp)<0.2&&ptype==0","same");
  
  

  canvasdSnp->Divide(3,2);
  chainPlus->SetMarkerStyle(25);
  chainPlus->SetMarkerSize(0.5);
  chainPlus->SetMarkerColor(1);
  canvasdSnp->cd(1);
  chainPlus->Draw("1000*deltaM:1000*delta",cut+"theta>0&&abs(snp)<0.2&&ptype==2","");
  canvasdSnp->cd(2);
  chainPlus->SetMarkerColor(1);
  chainPlus->Draw("1000*(deltaM-delta)>>deltaCorr(50,-4.,4)",cut+"theta>0&&abs(snp)<0.2&&ptype==2","err");
  chainPlus->SetMarkerColor(3);
  chainPlus->Draw("1000*deltaM",cut+"theta>0&&abs(snp)<0.2&&ptype==2","errsame");
  canvasdSnp->cd(3);
  chainPlus->SetMarkerColor(1);
  chainPlus->Draw("1000*(deltaM):phi",cut+"theta>0&&abs(snp)<0.2&&ptype==2","");
  chainPlus->SetMarkerColor(2);
  chainPlus->Draw("1000*(deltaM-delta):phi",cut+"theta>0&&abs(snp)<0.2&&ptype==2","same");
  canvasdSnp->cd(6);
  chainPlus->SetMarkerStyle(26);
  chainPlus->SetMarkerColor(1);  
  chainPlus->Draw("1000*(deltaM):phi",cut+"theta<0&&abs(snp)<0.2&&ptype==2","");
  chainPlus->SetMarkerColor(4);  
  chainPlus->Draw("1000*(deltaM-delta):phi",cut+"theta<0&&abs(snp)<0.2&&ptype==2","same");



  canvasd1pt->Divide(3,2);
  chainPlus->SetMarkerStyle(25);
  chainPlus->SetMarkerSize(0.5);
  chainPlus->SetMarkerColor(1);
  canvasd1pt->cd(1);
  chainPlus->Draw("deltaM:delta",cut+"theta>0&&abs(snp)<0.2&&ptype==4","");
  canvasd1pt->cd(2);
  chainPlus->SetMarkerColor(1);
  chainPlus->Draw("deltaM-delta>>deltaCorr(50,-0.15,0.15)",cut+"theta>0&&abs(snp)<0.2&&ptype==4","err");
  chainPlus->GetHistogram()->Fit("gaus");
  chainPlus->SetMarkerColor(3);
  chainPlus->Draw("deltaM",cut+"theta>0&&abs(snp)<0.2&&ptype==4","errsame");
  canvasd1pt->cd(3);
  chainPlus->SetMarkerColor(1);
  chainPlus->Draw("deltaM:phi",cut+"theta>0&&abs(snp)<0.2&&ptype==4","");
  chainPlus->SetMarkerColor(2);
  chainPlus->Draw("deltaM-delta:phi",cut+"theta>0&&abs(snp)<0.2&&ptype==4","same");
  canvasd1pt->cd(6);
  chainPlus->SetMarkerStyle(26);
  chainPlus->SetMarkerColor(1);  
  chainPlus->Draw("deltaM:phi",cut+"theta<0&&abs(snp)<0.2&&ptype==4","");
  chainPlus->SetMarkerColor(4);  
  chainPlus->Draw("deltaM-delta:phi",cut+"theta<0&&abs(snp)<0.2&&ptype==4","same");
  
  {for (Int_t ipad=0; ipad<3;ipad++){
    if (ipad==0) canvasdY->cd(5);
    if (ipad==1) canvasdSnp->cd(5);
    if (ipad==2) canvasd1pt->cd(5);
    lfit->SetTextAlign(12); lfit->SetTextSize(0.04);  
    {for (Int_t i=1; i<fitArray->GetEntries(); i++){
	lfit->DrawLatex(0.1,1-i/9.,fitArray->At(i)->GetName());
      }}
    if (ipad==0) canvasdY->cd(4);
    if (ipad==1) canvasdSnp->cd(4);
    if (ipad==2) canvasd1pt->cd(4);
    lfit->DrawLatex(0.1,0.9,"Global fit - TPC ITS matching in r-#phi");
    lfit->DrawLatex(0.1,0.8,"Residual minimization:");
    lfit->DrawLatex(0.1,0.7,"r#phi_{0.5T}-r#phi_{-0.5T} (cm)");
    lfit->DrawLatex(0.1,0.6,"sin(r#phi)_{0.5T}-sin(r#phi)_{-0.5T} (mrad)");
    lfit->DrawLatex(0.1,0.5,"1/pt_{0.5T}-1/pt_{-0.5T} (1/GeV)");
    }}


  TFile *fpic = new TFile(Form("fitPictures%d.root",highFrequency),"update");
  canvasd1pt->Write();
  canvasdY->Write();
  canvasdSnp->Write();
  //
  canvasd1pt->SaveAs(Form("fitd1pt%d.pdf",highFrequency));
  canvasdY->SaveAs(Form("fitdy%d.pdf",highFrequency));
  canvasdSnp->SaveAs(Form("fitdsnp%d.pdf",highFrequency));
  delete fpic;
}
