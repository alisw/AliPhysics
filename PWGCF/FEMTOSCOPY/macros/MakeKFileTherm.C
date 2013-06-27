#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>
#include <complex>

#include "TVector2.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TText.h"
#include "TRandom3.h"
#include "TArray.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMinuit.h"

using namespace std;

void MakeKFileTherm(){
  // 2-particles

  const int kbVALUES=6;// 6 b values (2,3,5,7,8,9)

  TFile *probeFile=new TFile("Therm_FSI_b2.root","READ");// all files have same binning.
  TH2D *probeHisto = (TH2D*)probeFile->Get("K2_ss");// kt x qinv
  int binsY= probeHisto->GetNbinsY();
  double lowY = probeHisto->GetYaxis()->GetBinLowEdge(1);
  double highY = probeHisto->GetYaxis()->GetBinUpEdge(binsY);
  probeFile->Close();

  TFile *OutFile=new TFile("KFile_temp.root","RECREATE");
  TH2D *K2ssT = new TH2D("K2ssT","",kbVALUES,0.5,kbVALUES+0.5, binsY,lowY,highY);// kt integrated
  TH2D *K2osT = new TH2D("K2osT","",kbVALUES,0.5,kbVALUES+0.5, binsY,lowY,highY);// kt integrated
  TH3D *K2ssT_kt = new TH3D("K2ssT_kt","",kbVALUES,0.5,kbVALUES+0.5, 6,0.5,6+0.5, binsY,lowY,highY);// kt differential
  TH3D *K2osT_kt = new TH3D("K2osT_kt","",kbVALUES,0.5,kbVALUES+0.5, 6,0.5,6+0.5, binsY,lowY,highY);// kt differential
  K2ssT_kt->GetXaxis()->SetTitle("b bin"); K2osT_kt->GetXaxis()->SetTitle("b bin");
  K2ssT_kt->GetYaxis()->SetTitle("kt bin"); K2osT_kt->GetYaxis()->SetTitle("kt bin");
  K2ssT_kt->GetZaxis()->SetTitle("qinv"); K2osT_kt->GetZaxis()->SetTitle("qinv");

  for(int r=0; r<kbVALUES; r++){ 
    
    TString *nameT=new TString("Therm_FSI_b");
        
     
    if(r==0) {*nameT += 2;}
    if(r==1) {*nameT += 3;}
    if(r==2) {*nameT += 5;}
    if(r==3) {*nameT += 7;}
    if(r==4) {*nameT += 8;}
    if(r==5) {*nameT += 9;}
    
    nameT->Append(".root");
        
   

    //
    TFile *file_T=new TFile(nameT->Data(),"READ");
    TH2D *Num2_ss = (TH2D*)file_T->Get("K2_ss");
    TH2D *Den2_ss = (TH2D*)file_T->Get("PlaneWF_ss");
    TH1D *Num_ss = (TH1D*)Num2_ss->ProjectionY();// kt integrated
    TH1D *Den_ss = (TH1D*)Den2_ss->ProjectionY();// kt integrated
    TH2D *Num2_os = (TH2D*)file_T->Get("K2_os");
    TH2D *Den2_os = (TH2D*)file_T->Get("PlaneWF_os");
    TH1D *Num_os = (TH1D*)Num2_os->ProjectionY();// kt integrated
    TH1D *Den_os = (TH1D*)Den2_os->ProjectionY();// kt integrated
    Num_ss->Divide(Den_ss);
    Num_os->Divide(Den_os);
    for(int i=1; i<=binsY; i++){
      K2ssT->SetBinContent(r+1, i, Num_ss->GetBinContent(i));
      K2osT->SetBinContent(r+1, i, Num_os->GetBinContent(i));
    }
   
    // kt differential
    for(int ktbin=1; ktbin<=6; ktbin++){
      TString *name=new TString("pro_");
      *name += ktbin;
      *name += 1;
      TH1D *Num_ss = (TH1D*)Num2_ss->ProjectionY(name->Data(), ktbin*2+3, ktbin*2+4);
      *name += 2;
      TH1D *Den_ss = (TH1D*)Den2_ss->ProjectionY(name->Data(), ktbin*2+3, ktbin*2+4);
      *name += 3;
      TH1D *Num_os = (TH1D*)Num2_os->ProjectionY(name->Data(), ktbin*2+3, ktbin*2+4);
      *name += 4;
      TH1D *Den_os = (TH1D*)Den2_os->ProjectionY(name->Data(), ktbin*2+3, ktbin*2+4);
      Num_ss->Divide(Den_ss);
      Num_os->Divide(Den_os);
       for(int i=1; i<=binsY; i++){
	 K2ssT_kt->SetBinContent(r+1, ktbin, i, Num_ss->GetBinContent(i));
	 K2osT_kt->SetBinContent(r+1, ktbin, i, Num_os->GetBinContent(i));
       }
    }  
    
    // 3-particles

    TH3D *Num3_ss=(TH3D*)file_T->Get("K3ss_3D");
    TH3D *Den3_ss=(TH3D*)file_T->Get("PlaneWF3ss_3D");
    Num3_ss->Divide(Den3_ss);
    TString *OutNameSS = new TString("K3ss_");
    *OutNameSS += r;
    OutFile->cd();
    TH3D *K3ss = (TH3D*)Num3_ss->Clone();
    TH3D *temp3ss = (TH3D*)Num3_ss->Clone();
    for(int i=1; i<=K3ss->GetNbinsX(); i++){
      for(int j=1; j<=K3ss->GetNbinsY(); j++){
	for(int k=1; k<=K3ss->GetNbinsZ(); k++){

	  // GRS
	  //double GRS = K2ssT->GetBinContent(r+1, i) * K2ssT->GetBinContent(r+1, j) * K2ssT->GetBinContent(r+1, k);
	  //K3ss->SetBinContent(i,j,k, GRS);

	  // Omega0
	  if(temp3ss->GetBinContent(i,j,k) > 1.0) K3ss->SetBinContent(i,j,k, 1.0);
	  else{
	    double mean=0;
	    double terms=0;
	    if(temp3ss->GetBinContent(i,j,k) < 1.0) {mean += temp3ss->GetBinContent(i,j,k); terms++;}
	    if(temp3ss->GetBinContent(i,k,j) < 1.0) {mean += temp3ss->GetBinContent(i,k,j); terms++;}
	    if(temp3ss->GetBinContent(j,i,k) < 1.0) {mean += temp3ss->GetBinContent(j,i,k); terms++;}
	    if(temp3ss->GetBinContent(j,k,i) < 1.0) {mean += temp3ss->GetBinContent(j,k,i); terms++;}
	    if(temp3ss->GetBinContent(k,i,j) < 1.0) {mean += temp3ss->GetBinContent(k,i,j); terms++;}
	    if(temp3ss->GetBinContent(k,j,i) < 1.0) {mean += temp3ss->GetBinContent(k,j,i); terms++;}
	    
	    if(terms > 0) {mean /= terms; K3ss->SetBinContent(i,j,k, mean);}
	    else K3ss->SetBinContent(i,j,k, 0);
	    }

	}
      }
    }
    K3ss->Write(OutNameSS->Data());
    //
    TH3D *Num3_os=(TH3D*)file_T->Get("K3os_3D");
    TH3D *Den3_os=(TH3D*)file_T->Get("PlaneWF3os_3D");
    Num3_os->Divide(Den3_os);
    //
    TString *OutNameOS = new TString("K3os_");
    *OutNameOS += r;
    OutFile->cd();
    TH3D *K3os = (TH3D*)Num3_os->Clone();
    TH3D *temp3os = (TH3D*)Num3_os->Clone();
    for(int i=1; i<=K3os->GetNbinsX(); i++){
      for(int j=1; j<=K3os->GetNbinsY(); j++){
	for(int k=1; k<=K3os->GetNbinsZ(); k++){

	  // GRS
	  //double GRS = K2ssT->GetBinContent(r+1, i) * K2osT->GetBinContent(r+1, j) * K2osT->GetBinContent(r+1, k);
	  //K3os->SetBinContent(i,j,k, GRS);
	  
	  // Omega0
	  if(temp3os->GetBinContent(i,j,k) > 3.0) K3os->SetBinContent(i,j,k, 1.0);
	  else {
	    double mean=0;
	    double terms=0;
	    if(temp3os->GetBinContent(i,j,k) < 3.0) {mean += temp3os->GetBinContent(i,j,k); terms++;}
	    if(temp3os->GetBinContent(i,k,j) < 3.0) {mean += temp3os->GetBinContent(i,k,j); terms++;}
	    
	    if(terms > 0) {mean /= terms; K3os->SetBinContent(i,j,k, mean);}
	    else K3os->SetBinContent(i,j,k, 0);
	    }
	  
	}
      }
    }
    
    K3os->Write(OutNameOS->Data());
    
    
    file_T->Close();
    
  }// r loop
  
  //
  OutFile->cd();
  //
  K2ssT->Write("K2ssT");
  K2osT->Write("K2osT");
  K2ssT_kt->Write("K2ssT_kt");
  K2osT_kt->Write("K2osT_kt");
  //
  
  OutFile->Close();

}
