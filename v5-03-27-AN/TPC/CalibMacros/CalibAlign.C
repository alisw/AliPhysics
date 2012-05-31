/*
  .x ~/UliStyle.C
  .x ~/rootlogon.C
  gSystem->Load("libSTAT.so");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib"); 
  gSystem->Load("libSTAT.so");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  


  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool; 
  TChain * chainCosmic = tool.MakeChain("cosmic.txt","Track0",0,1000000);
  chainCosmic->Lookup();
  chainCosmic->SetProof(kTRUE);
  TChain * chainTr = tool.MakeChain("align.txt","Tracklet",0,10200);
  chainTr->Lookup();
  //chainTr->SetProof(kTRUE);
  //
  //
  //
  .L $ALICE_ROOT/TPC/CalibMacros/CalibAlign.C
  SetAlias();
  InitCutsAlign();
  MakeAlign();
  
*/
/*
#include "TMath.h"
#include "TFile.h"
#include "TLinearFitter.h"
#include "TChain.h"
#include "TTreeStream.h"
#include "TStatToolkit.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TCut.h"
#include "TPad.h"
#include "AliTPCcalibAlign.h"
*/



AliTPCcalibAlign align;
TChain * chainTr;

void SetAlias(){
  //
  //
  //
  chainTr->SetAlias("dP0","tp1.fP[0]-tp2.fP[0]");
  chainTr->SetAlias("dP1","tp1.fP[1]-tp2.fP[1]");
  chainTr->SetAlias("dP2","tp1.fP[2]-tp2.fP[2]");
  chainTr->SetAlias("dP3","tp1.fP[3]-tp2.fP[3]");
  chainTr->SetAlias("dP4","tp1.fP[4]-tp2.fP[4]");
  //
  chainTr->SetAlias("sP0","sqrt(tp1.fC[0]+tp2.fC[0])");
  chainTr->SetAlias("sP1","sqrt(tp1.fC[2]+tp2.fC[2])");
  chainTr->SetAlias("sP2","sqrt(tp1.fC[5]+tp2.fC[5])");
  chainTr->SetAlias("sP3","sqrt(tp1.fC[9]+tp2.fC[9])");
  chainTr->SetAlias("sP4","sqrt(tp1.fC[14]+tp2.fC[14])");
  //
  chainTr->SetAlias("pP0","dP0/sP0");
  chainTr->SetAlias("pP1","dP1/sP1");
  chainTr->SetAlias("pP2","dP2/sP2");
  chainTr->SetAlias("pP3","dP3/sP3");
  //
  chainTr->SetAlias("side","(sign(tp1.fP[1])+0)");
  chainTr->SetAlias("dR","(1-abs(tp1.fP[1]/250.))");
  chainTr->SetAlias("ta0","(tp1.fP[2]+tp2.fP[2])*0.5");
  chainTr->SetAlias("ta1","(tp1.fP[3]+tp2.fP[3])*0.5");
  chainTr->SetAlias("ca","cos(tp1.fAlpha+0)");
  chainTr->SetAlias("sa","sin(tp1.fAlpha+0)");
  //
  chainTr->SetAlias("meanZ","(tp1.fP[1]+tp2.fP[1])*0.5");
  chainTr->SetAlias("vx1","(v1.fElements[0]+0)");
  chainTr->SetAlias("vy1","(v1.fElements[1]+0)");
  chainTr->SetAlias("vz1","(v1.fElements[2]+0)");
  chainTr->SetAlias("vdydx1","(v1.fElements[3]+0)");
  chainTr->SetAlias("vdzdx1","(v1.fElements[4]+0)");
  chainTr->SetAlias("vx2","(v2.fElements[0]+0)");
  chainTr->SetAlias("vy2","(v2.fElements[1]+0)");
  chainTr->SetAlias("vz2","(v2.fElements[2]+0)");
  chainTr->SetAlias("vdydx2","(v2.fElements[3]+0)");
  chainTr->SetAlias("vdzdx2","(v2.fElements[4]+0)");
  //
  //
  //
  chainTr->SetAlias("dx1","(AliTPCcalibAlign::SCorrect(0,0,s1,s2,vx1,vy1,vz1,vdydx1,vdzdx1)-vx1)");
  chainTr->SetAlias("dy1","(AliTPCcalibAlign::SCorrect(0,1,s1,s2,vx1,vy1,vz1,vdydx1,vdzdx1)-vy1)");
  chainTr->SetAlias("dz1","(AliTPCcalibAlign::SCorrect(0,2,s1,s2,vx1,vy1,vz1,vdydx1,vdzdx1)-vz1)");
  chainTr->SetAlias("ddy1","(AliTPCcalibAlign::SCorrect(0,3,s1,s2,vx1,vy1,vz1,vdydx1,vdzdx1)-vdydx1)");
  chainTr->SetAlias("ddz1","(AliTPCcalibAlign::SCorrect(0,4,s1,s2,vx1,vy1,vz1,vdydx1,vdzdx1)-vdzdx1)");
  chainTr->SetAlias("ddy01","(AliTPCcalibAlign::SCorrect(0,5,s1,s2,vx1,vy1,vz1,vdydx1,vdzdx1)-vdydx1)");
  chainTr->SetAlias("ddz01","(AliTPCcalibAlign::SCorrect(0,6,s1,s2,vx1,vy1,vz1,vdydx1,vdzdx1)-vdzdx1)");


}


void InitCutsAlign(){
  // resolution cuts
  TCut cutS0("sqrt(tp2.fC[0]+tp2.fC[0])<0.2");
  TCut cutS1("sqrt(tp2.fC[2]+tp2.fC[2])<0.2");
  TCut cutS2("sqrt(tp2.fC[5]+tp2.fC[5])<0.01");
  TCut cutS3("sqrt(tp2.fC[9]+tp2.fC[9])<0.01");
  TCut cutS4("sqrt(tp2.fC[14]+tp2.fC[14])<0.5");
  TCut cutS=cutS0+cutS1+cutS2+cutS3+cutS4;
  //
  // parameters matching cuts
  TCut cutP0("abs(tp1.fP[0]-tp2.fP[0])<0.6");
  TCut cutP1("abs(tp1.fP[1]-tp2.fP[1])<0.6");
  TCut cutP2("abs(tp1.fP[2]-tp2.fP[2])<0.03");
  TCut cutP3("abs(tp1.fP[3]-tp2.fP[3])<0.03");
  TCut cutP=cutP0+cutP1+cutP2+cutP3;
  //
  TCut cutAll = cutS+cutP;
  chainTr->Draw(">>listELP",cutAll,"entryList");
  TEntryList *elist = (TEntryList*)gDirectory->Get("listELP");
  chainTr->SetEntryList(elist);
  
  TCut cutRun("1");
  TCut cutN120("1");

}

void MakeAlign(){
  
  align.ProcessTree(chainTr);
  align.EvalFitters();
  align.MakeTree("alignTree.root");
  align.SetInstance(&align);
}
 
void MakeCompareAlign(){
  //
  //
  //
  TFile falignTreeNoMag("/lustre_alpha/alice/miranov/rec/LHC08d/nomag/alignTree.root");
  TTree * treeAlignNoMag = (TTree*)falignTreeNoMag.Get("Align");
  TFile falignTree("alignTree.root");
  TTree * treeAlign = (TTree*)falignTree.Get("Align");
  treeAlign->AddFriend(treeAlignNoMag,"NoMag");
  treeAlignNoMag->SetMarkerStyle(26);
  treeAlign->SetMarkerStyle(25);
}

TMatrixD * arrayAlign[72];
TMatrixD * arrayAlignTmp[72];

void ClearMatrices(){
  //
  //
  for (Int_t i=0;i<72; i++) {
    TMatrixD * mat = new TMatrixD(4,4);
    mat->UnitMatrix();
    arrayAlign[i]=mat;
    arrayAlignTmp[i]=(TMatrixD*)(mat->Clone());
  }

}

void GlobalAlign(){
  //
  // Global Align
  //
  TTreeSRedirector *cstream = new TTreeSRedirector("galign.root");

  for (Int_t iter=0; iter<10;iter++){
    printf("Iter=\t%d\n",iter);
    for (Int_t is0=0;is0<72; is0++) {
      //
      TMatrixD  *mati0 = arrayAlign[is0];
      TMatrixD matDiff(4,4);
      Double_t sumw=0;      
      for (Int_t is1=0;is1<72; is1++) {
	Bool_t invers=kFALSE;
	const TMatrixD *mat = align.GetTransformation(is0,is1,0); 
	if (!mat){
	  invers=kTRUE;
	  mat = align.GetTransformation(is1,is0,0); 
	}
	if (!mat) continue;
	Double_t weight=1;
	if  ( (is1%18-is0%18)!=0) weight*=0.3;
	if (is1/36>is0/36) weight*=2./3.; 
	if (is1/36<is0/36) weight*=1./3.;
	//
	//
	TMatrixD matT = *mat;	
	if (invers)  matT.Invert();
	matDiff+=weight*(*(arrayAlign[is1]))*matT;
	sumw+=weight;
	(*cstream)<<"LAlign"<<
	  "iter="<<iter<<
	  "s0="<<is0<<
	  "m6.="<<arrayAlign[is0]<<
	  "\n";

      }
      if (sumw>0){
	matDiff*=1/sumw;
	(*arrayAlignTmp[is0]) = matDiff;
      }
    }
    for (Int_t is0=0;is0<72; is0++) {
      (*arrayAlign[is0]) = (*arrayAlignTmp[is0]);
      //      TMatrixD * matM1=  align.GetTransformation(is0,36+(is0+35)%36,0);
      //TMatrixD * mat  =  align.GetTransformation(is0,36+(is0+36)%36,0);
      //TMatrixD * matP1=  align.GetTransformation(is0,36+(is0+37)%36,0);
      //
      (*cstream)<<"GAlign"<<
	"iter="<<iter<<
	"s0="<<is0<<
	"m6.="<<arrayAlign[is0]<<
	"\n";
    }    
  }
  delete cstream;
}





void MakeGlobalCorr(){
  //
  //
  //
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParam;
  TMatrixD covMatrix;
  TVectorD chi2V(5);
  //
  TString fstring="";
  fstring+="side++";
  //
  fstring+="dR++";
  fstring+="dR*dR++";
  fstring+="dR*sa++";
  fstring+="dR*ca++";
  fstring+="ta0++";
  fstring+="ta1++";
  //
  fstring+="dR*side++";
  fstring+="dR*dR*side++";
  fstring+="dR*sa*side++";
  fstring+="dR*ca*side++";
  fstring+="ta0*side++";
  fstring+="ta1*side++";

  TString * strP0 = TStatToolkit::FitPlane(chainTr,"dP0:sP0", fstring.Data(),cutS, chi2,npoints,fitParam,covMatrix);
  chi2V[0]=TMath::Sqrt(chi2/npoints);
  chainTr->SetAlias("corrP0",strP0->Data());

  TString * strP1 = TStatToolkit::FitPlane(chainTr,"dP1:sP1", fstring.Data(),cutS, chi2,npoints,fitParam,covMatrix);
  chi2V[1]=TMath::Sqrt(chi2/npoints);
  chainTr->SetAlias("corrP1",strP1->Data());

  TString * strP2 = TStatToolkit::FitPlane(chainTr,"dP2:sP2", fstring.Data(),cutS, chi2,npoints,fitParam,covMatrix);
  chi2V[2]=TMath::Sqrt(chi2/npoints);
  chainTr->SetAlias("corrP2",strP2->Data());

  TString * strP3 = TStatToolkit::FitPlane(chainTr,"dP3:sP3", fstring.Data(),cutS, chi2,npoints,fitParam,covMatrix);
  chi2V[3]=TMath::Sqrt(chi2/npoints);
  chainTr->SetAlias("corrP3",strP3->Data());

  TString * strP4 = TStatToolkit::FitPlane(chainTr,"dP4:sP4", fstring.Data(),cutS, chi2,npoints,fitParam,covMatrix);
  chi2V[4]=TMath::Sqrt(chi2/npoints);
  chainTr->SetAlias("corrP4",strP4->Data());
}

void P0resolZ(){
  //
  //
  //
  TH2F * hdP0Z = new TH2F("hdP0Z","hdP0Z",10,-250,250,100,-0.5,0.5);
  TH2F * hdP0ZNoCor = new TH2F("hdP0ZNoCor","hdP0ZNoCor",10,-250,250,100,-0.5,0.5);
  chainTr->Draw("(dP0-corrP0)/sqrt(2.):meanZ>>hdP0Z",""+cutRun+cutS+cutN120,"");
  chainTr->Draw("(dP0-0)/sqrt(2.):meanZ>>hdP0ZNoCor",""+cutRun+cutS+cutN120,"");
  
  hdP0Z->FitSlicesY();
  hdP0ZNoCor->FitSlicesY();
  hdP0Z_2->SetMinimum(0);
  hdP0Z_2->SetXTitle("Z position (cm)");
  hdP0Z_2->SetYTitle("#sigma_{y} (cm)");
  hdP0Z_2->SetMarkerStyle(25);
  hdP0ZNoCor_2->SetMarkerStyle(26);
  hdP0Z_2->Draw();
  hdP0ZNoCor_2->Draw("same");  
  gPad->SaveAs("picAlign/SigmaY_z.gif");
  gPad->SaveAs("picAlign/SigmaY_z.eps");
  //
  hdP0ZNoCor_1->SetXTitle("Z position (cm)");
  hdP0ZNoCor_1->SetYTitle("#Delta{y} (cm)");
  hdP0Z_1->SetMarkerStyle(25);
  hdP0ZNoCor_1->SetMarkerStyle(26);
  hdP0ZNoCor_1->Draw("");  
  hdP0Z_1->Draw("same");
  gPad->SaveAs("picAlign/DeltaY_z.gif");
  gPad->SaveAs("picAlign/DeltaY_z.eps");
  //
  //
  TH2F * hdPP0Z = new TH2F("hdPP0Z","hdPP0Z",8,-200,200,50,-5.05,5.05);
  TH2F * hncdPP0Z = new TH2F("hncdPP0Z","hncdPP0Z",8,-200,200,50,-5.05,5.05);
  chainTr->Draw("(dP0-corrP0)/sP0:meanZ>>hdPP0Z",""+cutRun+cutS+cutN120,"");
  chainTr->Draw("(dP0-0)/sP0:meanZ>>hncdPP0Z",""+cutRun+cutS+cutN120,"");
  hdPP0Z->FitSlicesY();
  hncdPP0Z->FitSlicesY();
  hdPP0Z_2->SetMarkerStyle(25);
  hncdPP0Z_2->SetMarkerStyle(26);
  hdPP0Z_2->SetMinimum(0);
  hdPP0Z_2->SetXTitle("Z position (cm)");
  hdPP0Z_2->SetYTitle("#sigma y (Unit) ");
  hdPP0Z_2->Draw();
  hncdPP0Z_2->Draw("same");
  gPad->SaveAs("picAlign/PullY_z.gif");
  gPad->SaveAs("picAlign/PullY_z.eps");
}

void P1resolZ(){
  //
  //
  //
  TH2F * hdP1Z = new TH2F("hdP1Z","hdP1Z",10,-250,250,100,-0.2,0.2);
  TH2F * hdP1ZNoCor = new TH2F("hdP1ZNoCor","hdP1ZNoCor",10,-250,250,100,-0.2,0.2);
  chainTr->Draw("(dP1-corrP1)/sqrt(2.):meanZ>>hdP1Z",""+cutRun+cutS+cutN120,"");
  chainTr->Draw("(dP1-0)/sqrt(2.):meanZ>>hdP1ZNoCor",""+cutRun+cutS+cutN120,"");
  
  hdP1Z->FitSlicesY();
  hdP1ZNoCor->FitSlicesY();
  hdP1Z_2->SetMinimum(0);
  hdP1Z_2->SetXTitle("Z position (cm)");
  hdP1Z_2->SetYTitle("#sigma_{z} (cm)");
  hdP1Z_2->SetMarkerStyle(25);
  hdP1ZNoCor_2->SetMarkerStyle(26);
  hdP1Z_2->Draw();
  hdP1ZNoCor_2->Draw("same");  
  gPad->SaveAs("picAlign/SigmaZ_z.gif");
  gPad->SaveAs("picAlign/SigmaZ_z.eps");
  //
  hdP1ZNoCor_1->SetXTitle("Z position (cm)");
  hdP1ZNoCor_1->SetYTitle("#Delta{z} (cm)");
  hdP1Z_1->SetMarkerStyle(25);
  hdP1ZNoCor_1->SetMarkerStyle(26);
  hdP1ZNoCor_1->Draw("");  
  hdP1Z_1->Draw("same");
  gPad->SaveAs("picAlign/DeltaZ_z.gif");
  gPad->SaveAs("picAlign/DeltaZ_z.eps");
  //
  //
  TH2F * hdPP1Z = new TH2F("hdPP1Z","hdPP1Z",8,-200,200,50,-5.05,5.05);
  TH2F * hncdPP1Z = new TH2F("hncdPP1Z","hncdPP1Z",8,-200,200,50,-5.05,5.05);
  chainTr->Draw("(dP1-corrP1)/sP1:meanZ>>hdPP1Z",""+cutRun+cutS+cutN120,"");
  chainTr->Draw("(dP1-0)/sP1:meanZ>>hncdPP1Z",""+cutRun+cutS+cutN120,"");
  hdPP1Z->FitSlicesY();
  hncdPP1Z->FitSlicesY();
  hdPP1Z_2->SetMarkerStyle(25);
  hncdPP1Z_2->SetMarkerStyle(26);
  hdPP1Z_2->SetMinimum(0);
  hdPP1Z_2->SetXTitle("Z position (cm)");
  hdPP1Z_2->SetYTitle("#sigma z (Unit) ");
  hdPP1Z_2->Draw();
  hncdPP1Z_2->Draw("same");
  gPad->SaveAs("picAlign/PullZ_z.gif");
  gPad->SaveAs("picAlign/PullZ_z.eps");
}


void P4resolZ(){
  //
  //
  //
  TH2F * hdP4Z = new TH2F("hdP4Z","hdP4Z",10,-250,250,100,-0.4,0.4);
  TH2F * hdP4ZNoCor = new TH2F("hdP4ZNoCor","hdP4ZNoCor",10,-250,250,100,-0.4,0.4);
  chainTr->Draw("(dP4-corrP4)/sqrt(2.):meanZ>>hdP4Z",""+cutRun+cutS+cutN120,"");
  chainTr->Draw("(dP4-0)/sqrt(2.):meanZ>>hdP4ZNoCor",""+cutRun+cutS+cutN120,"");
  
  hdP4Z->FitSlicesY();
  hdP4ZNoCor->FitSlicesY();
  hdP4Z_2->SetMinimum(0);
  hdP4Z_2->SetXTitle("Z position (cm)");
  hdP4Z_2->SetYTitle("#sigma_{1/pt} (1/GeV)");
  hdP4Z_2->SetMarkerStyle(25);
  hdP4ZNoCor_2->SetMarkerStyle(26);
  hdP4Z_2->Draw();
  hdP4ZNoCor_2->Draw("same");  
  gPad->SaveAs("picAlign/SigmaP4_z.gif");
  gPad->SaveAs("picAlign/SigmaP4_z.eps");
  //
  hdP4ZNoCor_1->SetXTitle("Z position (cm)");
  hdP4ZNoCor_1->SetYTitle("#Delta_{1/p_{t}} (1/GeV)");
  hdP4Z_1->SetMarkerStyle(25);
  hdP4ZNoCor_1->SetMarkerStyle(26);
  hdP4ZNoCor_1->Draw("");  
  hdP4Z_1->Draw("same");
  gPad->SaveAs("picAlign/DeltaP4_z.gif");
  gPad->SaveAs("picAlign/DeltaP4_z.eps");
  //
  //
  TH2F * hdPP4Z = new TH2F("hdPP4Z","hdPP4Z",8,-200,200,50,-5.05,5.05);
  TH2F * hncdPP4Z = new TH2F("hncdPP4Z","hncdPP4Z",8,-200,200,50,-5.05,5.05);
  chainTr->Draw("(dP4-corrP4)/sP4:meanZ>>hdPP4Z",""+cutRun+cutS+cutN120,"");
  chainTr->Draw("(dP4-0)/sP4:meanZ>>hncdPP4Z",""+cutRun+cutS+cutN120,"");
  hdPP4Z->FitSlicesY();
  hncdPP4Z->FitSlicesY();
  hdPP4Z_2->SetMarkerStyle(25);
  hncdPP4Z_2->SetMarkerStyle(26);
  hdPP4Z_2->SetMinimum(0);
  hdPP4Z_2->SetXTitle("Z position (cm)");
  hdPP4Z_2->SetYTitle("#sigma 1/p_{t} (Unit) ");
  hdPP4Z_2->Draw();
  hncdPP4Z_2->Draw("same");
  gPad->SaveAs("picAlign/PullP4_z.gif");
  gPad->SaveAs("picAlign/PullP4_z.eps");
}





void MakeFit(){
  //
  //
  //
  TChain *chainTracklet=chainTr;
  AliTPCcalibAlign align;
  //
  TVectorD * vec1 = 0;
  TVectorD * vec2 = 0;
  AliExternalTrackParam * tp1 = 0;
  AliExternalTrackParam * tp2 = 0;  
  Int_t      s1 = 0;
  Int_t      s2 = 0;
  {
    for (Int_t i=0; i< elist->GetN(); i++){
      //for (Int_t i=0; i<100000; i++){
      chainTracklet->GetBranch("tp1.")->SetAddress(&tp1);
      chainTracklet->GetBranch("tp2.")->SetAddress(&tp2);
      chainTracklet->GetBranch("v1.")->SetAddress(&vec1);
      chainTracklet->GetBranch("v2.")->SetAddress(&vec2);
      chainTracklet->GetBranch("s1")->SetAddress(&s1);
      chainTracklet->GetBranch("s2")->SetAddress(&s2);      
      chainTracklet->GetEntry(i);
      if (i%100==0) printf("%d\t%d\t%d\t\n",i, s1,s2);
      //vec1.Print();
      TLinearFitter * fitter = align.GetOrMakeFitter6(s1,s2);
      if (fitter) align.Process6(vec1->GetMatrixArray(),vec2->GetMatrixArray(),fitter);
    }
  }
}


// chainTr->Scan("vy1:AliTPCcalibAlign::SCorrect(0,1,s1,s2,vx1,vy1,vz1,vdydx1,vdzdx1)","","",100); 

void MakePlot(){
  //
  chainTr->Draw("vy2-vy1:s1>>hisdYNoCor(36,0,36,100,-0.3,0.3)","abs(s2-s1-36)<1","",40000);
  chainTr->Draw("vy2-vy1-dy1:s1>>hisdYCor(36,0,36,100,-0.3,0.3)","abs(s2-s1-36)<1","",40000);
  hisdYCor->FitSlicesY();
  hisdYNoCor->FitSlicesY();
  hisdYCor_1->Draw("");

}


void dPhi(){
  //
  //
  //  
  treeAlign->Draw("1000*m6.fElements[4]","s2==s1+36&&nphi>100");
  htemp->SetXTitle("#Delta_{#phi} (mrad)");
  gPad->SaveAs("picAlign/mag5dPhi.eps");
  gPad->SaveAs("picAlign/mag5dPhi.gif");
  //
  treeAlignNoMag->Draw("1000*m6.fElements[4]","s2==s1+36&&nphi>100");
  htemp->SetXTitle("#Delta_{#phi} (mrad)");
  gPad->SaveAs("picAlign/nomagdPhi.eps");
  gPad->SaveAs("picAlign/nomagdPhi.gif");
  //
  //
  //
  treeAlign->Draw("1000*(m6.fElements[4]-NoMag.m6.fElements[4])","s2==s1+36&&nphi>100");
  htemp->SetXTitle("#Delta_{#phi} (mrad)");
  gPad->SaveAs("picAlign/diffnomagmag5dPhi.eps");
  gPad->SaveAs("picAlign/diffnomagmag5dPhi.gif");
  //
}


void dTheta(){
  //
  //
  //  
  treeAlign->Draw("1000*m6.fElements[8]","s2==s1+36&&nphi>100");
  htemp->SetXTitle("#Delta_{#theta} (mrad)");
  gPad->SaveAs("picAlign/mag5dTheta.eps");
  gPad->SaveAs("picAlign/mag5dTheta.gif");
  //
  treeAlignNoMag->Draw("1000*m6.fElements[8]","s2==s1+36&&nphi>100");
  htemp->SetXTitle("#Delta_{#theta} (mrad)");
  gPad->SaveAs("picAlign/nomagdTheta.eps");
  gPad->SaveAs("picAlign/nomagdTheta.gif");
  //
  //
  //
  treeAlign->Draw("1000*(m6.fElements[8]-NoMag.m6.fElements[8])","s2==s1+36&&nphi>100");
  htemp->SetXTitle("#Delta_{#theta} (mrad)");
  gPad->SaveAs("picAlign/diffnomagmag5dTheta.eps");
  gPad->SaveAs("picAlign/diffnomagmag5dTheta.gif");
  //
  //
  treeAlign->Draw("1000*(m6.fElements[8]):s1","s2==s1+36&&nphi>100");
  htemp->SetYTitle("#Delta_{#theta} (mrad)");
  htemp->SetXTitle("Sector number");
  treeAlignNoMag->Draw("1000*(m6.fElements[8]):s1","s2==s1+36&&nphi>100","same");
  gPad->SaveAs("picAlign/diffnomagmag5dTheta_s1.eps");
  gPad->SaveAs("picAlign/diffnomagmag5dTheta_s1.gif");
  //
}



void dZ(){
  //
  //
  //  
  treeAlign->Draw("dz","s2==s1+36&&nphi>100");
  htemp->SetXTitle("#Delta_{Z} (cm)");
  gPad->SaveAs("picAlign/mag5dZ.eps");
  gPad->SaveAs("picAlign/mag5dZ.gif");
  //
  treeAlignNoMag->Draw("dz","s2==s1+36&&nphi>100");
  htemp->SetXTitle("#Delta_{Z} (cm)");
  gPad->SaveAs("picAlign/nomagdZ.eps");
  gPad->SaveAs("picAlign/nomagdZ.gif");
  //
  //
  //
  treeAlign->Draw("(dz-NoMag.dz)","s2==s1+36&&nphi>100");
  htemp->SetXTitle("#Delta_{Z} (cm)");
  gPad->SaveAs("picAlign/diffnomagmag5dZ.eps");
  gPad->SaveAs("picAlign/diffnomagmag5dZ.gif");
  //
  //
  treeAlign->Draw("dz:s1","s2==s1+36&&nphi>100");
  htemp->SetYTitle("#Delta_{Z} (cm)");
  htemp->SetXTitle("Sector number");
  treeAlignNoMag->Draw("dz:s1","s2==s1+36&&nphi>100","same");
  gPad->SaveAs("picAlign/diffnomagmag5dZ_s1.eps");
  gPad->SaveAs("picAlign/diffnomagmag5dZ_s1.gif");
  //
}

void dY(){
  //
  //
  //  
  treeAlign->Draw("dy","s2==s1+36&&nphi>100");
  htemp->SetXTitle("#Delta_{Y} (cm)");
  gPad->SaveAs("picAlign/mag5dY.eps");
  gPad->SaveAs("picAlign/mag5dY.gif");
  //
  treeAlignNoMag->Draw("dy","s2==s1+36&&nphi>100");
  htemp->SetXTitle("#Delta_{Y} (cm)");
  gPad->SaveAs("picAlign/nomagdY.eps");
  gPad->SaveAs("picAlign/nomagdY.gif");
  //
  //
  //
  treeAlign->Draw("(dy-NoMag.dy)","s2==s1+36&&nphi>100");
  htemp->SetXTitle("#Delta_{Y} (cm)");
  gPad->SaveAs("picAlign/diffnomagmag5dY.eps");
  gPad->SaveAs("picAlign/diffnomagmag5dY.gif");
  //
  //
  treeAlign->Draw("dy:s1","s2==s1+36&&nphi>100");
  htemp->SetYTitle("#Delta_{Y} (cm)");
  htemp->SetXTitle("Sector number");
  treeAlignNoMag->Draw("dy:s1","s2==s1+36&&nphi>100","same");
  gPad->SaveAs("picAlign/diffnomagmag5dY_s1.eps");
  gPad->SaveAs("picAlign/diffnomagmag5dY_s1.gif");
  //
}









