#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>

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


void MakeWeightFile()
{
  
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  

  //TFile *InputFile = new TFile("Results/RawWeightFile_11h.root","READ");
  //TFile *InputFile = new TFile("Results/RawWeightFile_genSignal_Rinv11.root","READ");
  //TFile *InputFile = new TFile("Results/RawWeightFile_genSignal_Rinv11_Smeared.root","READ");
  //TFile *InputFile = new TFile("Results/RawWeightFile_11h_PID1p5.root","READ");
  //TFile *InputFile = new TFile("Results/RawWeightFile_FB5.root","READ");
  //TFile *InputFile = new TFile("Results/RawWeightFile_11h_Chi3p1.root","READ");
  //TFile *InputFile = new TFile("Results/PDC_11h_standard_and_Raw0p04TTC.root","READ");// standard weights
  //TFile *InputFile = new TFile("Results/RawWeightFile_11h_4ktbins_new.root","READ");
  TFile *InputFile = new TFile("Results/RawWeightFile_11h_4ktbins_2sigma_3sigmaTTC.root","READ");// _1(2sigmaTTC), _2(3sigmaTTC)

  //TFile *InputFile = new TFile("MyOutput.root","READ");
  TDirectoryFile *tdir = (TDirectoryFile*)InputFile->Get("PWGCF.outputChaoticityAnalysis.root");
  TList *MyList=(TList*)tdir->Get("ChaoticityOutput_2");
  //TList *MyList=(TList*)InputFile->Get("MyList");
  InputFile->Close();

  const int KtBins=4;
  const int KyBins=1;
  const int EDBins=1;
  const int MBins=10;//10, make sure MB assignment below is correct!!!!!!!!!!!!!!!
  
  TFile *OutFile = new TFile("WeightFile_temp.root","RECREATE");
  TH3F *WeightHistos[KtBins][MBins];
  
  for(int ktB=1; ktB<=KtBins; ktB++){
    for(int MB=1; MB<=MBins; MB++){
      TString *InNameNum = new TString("TwoPart_num_Kt_");
      *InNameNum += ktB-1;
      InNameNum->Append("_Ky_0_M_");
      if(MB<=10) *InNameNum += MB-1;// make sure MB assignment is correct!!!!!!!!!!!!!!!
      else *InNameNum += 1;
      InNameNum->Append("_ED_0");
      //
      TString *InNameDen = new TString("TwoPart_den_Kt_");
      *InNameDen += ktB-1;
      InNameDen->Append("_Ky_0_M_");
      if(MB<=10) *InNameDen += MB-1;// make sure MB assignment is correct!!!!!!!!!!!!!!!
      else *InNameDen += 1;
      InNameDen->Append("_ED_0");
      //
      TH3D *tempNum = (TH3D*)MyList->FindObject(InNameNum->Data());
      TH3D *tempDen = (TH3D*)MyList->FindObject(InNameDen->Data());
      //
      int NormBinStart = tempNum->GetXaxis()->FindBin(0.135);//0.135
      int NormBinEnd = tempNum->GetXaxis()->FindBin(0.2);//0.2
      double Norm = tempNum->Integral(NormBinStart,NormBinEnd, NormBinStart,NormBinEnd, NormBinStart,NormBinEnd);
      Norm /= tempDen->Integral(NormBinStart,NormBinEnd, NormBinStart,NormBinEnd, NormBinStart,NormBinEnd);
      cout<<"Normalization = "<<Norm<<endl;
      cout<<tempNum->Integral(NormBinStart,NormBinEnd, NormBinStart,NormBinEnd, NormBinStart,NormBinEnd)<<"  "<<tempDen->Integral(NormBinStart,NormBinEnd, NormBinStart,NormBinEnd, NormBinStart,NormBinEnd)<<endl;
      //
      TString *OutNameWeight = new TString("Weight_Kt_");
      *OutNameWeight += ktB-1;
      OutNameWeight->Append("_Ky_0_M_");
      *OutNameWeight += MB-1;
      OutNameWeight->Append("_ED_0");
      
      int Nbins=tempNum->GetNbinsX();
      double QLimit = tempNum->GetXaxis()->GetBinUpEdge(Nbins);
      WeightHistos[ktB-1][MB-1] = new TH3F(OutNameWeight->Data(),"r3 Weights", Nbins,0,QLimit, Nbins,0,QLimit, Nbins,0,QLimit);
      WeightHistos[ktB-1][MB-1]->GetXaxis()->SetTitle("out");
      WeightHistos[ktB-1][MB-1]->GetYaxis()->SetTitle("side");
      WeightHistos[ktB-1][MB-1]->GetZaxis()->SetTitle("long");
      WeightHistos[ktB-1][MB-1]->GetXaxis()->SetTitleOffset(1.8);
      WeightHistos[ktB-1][MB-1]->GetYaxis()->SetTitleOffset(1.8);
      WeightHistos[ktB-1][MB-1]->GetZaxis()->SetTitleOffset(1.8);
      for(int outB=1; outB<=Nbins; outB++){
	for(int sideB=1; sideB<=Nbins; sideB++){
	  for(int longB=1; longB<=Nbins; longB++){
	    double weight=1, weight_e=0;
	    if(tempDen->GetBinContent(outB,sideB,longB) > 0 && tempNum->GetBinContent(outB,sideB,longB) > 0) {
	      //if(outB==sideB && outB==longB) cout<<outB<<"  "<<tempNum->GetBinContent(outB,sideB,longB)<<"  "<<tempDen->GetBinContent(outB,sideB,longB)<<endl;
	      weight = double(tempNum->GetBinContent(outB,sideB,longB))/double(tempDen->GetBinContent(outB,sideB,longB)) / Norm;
	      weight_e = pow(sqrt(double(tempNum->GetBinContent(outB,sideB,longB)))/double(tempDen->GetBinContent(outB,sideB,longB)) / Norm,2);
	      weight_e += pow(sqrt(double(tempDen->GetBinContent(outB,sideB,longB)))*double(tempNum->GetBinContent(outB,sideB,longB))/pow(double(tempDen->GetBinContent(outB,sideB,longB)),2) / Norm,2);
	      weight_e = sqrt(weight_e);
	      //if(weight < 0.8 && weight > 0.1) cout<<outB<<"  "<<sideB<<"  "<<longB<<"  "<<tempNum->GetBinContent(outB,sideB,longB)<<"  "<<tempDen->GetBinContent(outB,sideB,longB)<<endl;
	    }
	    if(weight==0){
	      WeightHistos[ktB-1][MB-1]->SetBinContent(outB,sideB,longB, 0);
	      WeightHistos[ktB-1][MB-1]->SetBinError(outB,sideB,longB, 0);
	    }else {
	      WeightHistos[ktB-1][MB-1]->SetBinContent(outB,sideB,longB, weight-1.0);// difference from unity
	      WeightHistos[ktB-1][MB-1]->SetBinError(outB,sideB,longB, weight_e);
	    }

	  }
	}
      }
      
      WeightHistos[ktB-1][MB-1]->Write();

    }// MB
  }// ktB

  int BOI1=2;
  int BOI2=2;
  TH3D *histoOI = WeightHistos[0][0]->Clone();
  histoOI->SetDirectory(0);
  //OutFile->Close();
  
  TH1D *pro = WeightHistos[0][0]->ProjectionY("pro",BOI1,BOI1,BOI2,BOI2);
  pro->Draw();

  /*TFile *projections=new TFile("projections.root","RECREATE");
  histoOI->GetXaxis()->SetRange(BOI1,BOI1);
  TH2D *proYZ = (TH2D*)histoOI->Project3D("yz");
  histoOI->GetXaxis()->SetRange(1,40);
  histoOI->GetYaxis()->SetRange(BOI1,BOI1);
  TH2D *proXZ = (TH2D*)histoOI->Project3D("xz");
  histoOI->GetYaxis()->SetRange(1,40);
  histoOI->GetZaxis()->SetRange(BOI1,BOI1);
  TH2D *proXY = (TH2D*)histoOI->Project3D("xy");
  proYZ->GetXaxis()->SetRangeUser(0,0.1);
  proYZ->GetYaxis()->SetRangeUser(0,0.1);
  proXZ->GetXaxis()->SetRangeUser(0,0.1);
  proXZ->GetYaxis()->SetRangeUser(0,0.1);
  proXY->GetXaxis()->SetRangeUser(0,0.1);
  proXY->GetYaxis()->SetRangeUser(0,0.1);
  proXY->Draw("lego");
  
  proYZ->Write();
  proXZ->Write();
  proXY->Write();
  */

}
