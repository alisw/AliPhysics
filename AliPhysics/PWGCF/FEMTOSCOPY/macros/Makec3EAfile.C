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
#include "TCanvas.h"
#include "TPad.h"

#define BohrR 1963.6885
#define FmToGeV 0.19733 // conversion to fm
#define PI 3.1415926
#define masspiC 0.1395702 // pi+ mass (GeV/c^2)

using namespace std;



void Makec3EAfile(){

  int ExpanType=0;// 0(EW) or 1(LG)
  int RadiiCount=7;

  TFile *infile;
  TMinuit *fit;
  //
  TH3D *PbPbEA[2];
  TH3D *pPbEA[2];
  TH3D *ppEA[2];
  //
  PbPbEA[0] = new TH3D("PbPbEA_c3","",RadiiCount,0.5,RadiiCount+0.5, 6,0.5,6.5, 50,-0.5,49.5);// Rcoh type, parNum, Gindex
  PbPbEA[0]->SetDirectory(0);
  //
  pPbEA[0] = new TH3D("pPbEA_c3","",RadiiCount,0.5,RadiiCount+0.5, 6,0.5,6.5, 50,-0.5,49.5);
  pPbEA[0]->SetDirectory(0);
  //
  ppEA[0] = new TH3D("ppEA_c3","",RadiiCount,0.5,RadiiCount+0.5, 6,0.5,6.5, 50,-0.5,49.5);
  ppEA[0]->SetDirectory(0);
  //
  //
  PbPbEA[1] = new TH3D("PbPbEA_C3","",RadiiCount,0.5,RadiiCount+0.5, 6,0.5,6.5, 50,-0.5,49.5);// Rcoh type, parNum, Gindex
  PbPbEA[1]->SetDirectory(0);
  //
  pPbEA[1] = new TH3D("pPbEA_C3","",RadiiCount,0.5,RadiiCount+0.5, 6,0.5,6.5, 50,-0.5,49.5);
  pPbEA[1]->SetDirectory(0);
  //
  ppEA[1] = new TH3D("ppEA_C3","",RadiiCount,0.5,RadiiCount+0.5, 6,0.5,6.5, 50,-0.5,49.5);
  ppEA[1]->SetDirectory(0);
  //

  
  
  //
  //////////////////////////////
  double value=0, value_e=0;
  // 
  for(int FT=0; FT<2; FT++){// c3 or C3
    //
    for(int Gindex=0; Gindex<=25; Gindex++){
      // PbPb
      for(int RT=0; RT<RadiiCount; RT++){// Rcoh type
	
	TString *name1 = new TString("FitFiles/");
	//if(FT==0) name1->Append("c3/");
	//else name1->Append("C3/");
	name1->Append("FitFile_CT0_FT");
	*name1 += ExpanType;
	name1->Append("_R");
	if(RT<RadiiCount-1) *name1 += RT;
	else *name1 += 100;
	name1->Append("_G");
	*name1 += Gindex;
	name1->Append(".root");
	infile = new TFile(name1->Data(),"READ");
	if(FT==0) fit = (TMinuit*)infile->Get("MyMinuit_c3");
	else fit = (TMinuit*)infile->Get("MyMinuit_C3");
	for(int parNum=0; parNum<6; parNum++){
	  fit->GetParameter(parNum+1, value,value_e);
	  PbPbEA[FT]->SetBinContent(RT+1, parNum+1, Gindex+1, value);
	}
	infile->Close();
   
	//
	// pPb
	TString *name1 = new TString("FitFiles/");
	//if(FT==0) name1->Append("c3/");
	//else name1->Append("C3/");
	name1->Append("FitFile_CT1_FT");
	*name1 += ExpanType;
	name1->Append("_R");
	if(RT<RadiiCount-1) *name1 += RT;
	else *name1 += 100;
	name1->Append("_G");
	*name1 += Gindex;
	name1->Append(".root");
	infile = new TFile(name1->Data(),"READ");
	if(FT==0) fit = (TMinuit*)infile->Get("MyMinuit_c3");
	else fit = (TMinuit*)infile->Get("MyMinuit_C3");
	for(int parNum=0; parNum<6; parNum++){
	  fit->GetParameter(parNum+1, value,value_e);
	  pPbEA[FT]->SetBinContent(RT+1, parNum+1, Gindex+1, value);
	}
	infile->Close();
	//
	//
	TString *name1 = new TString("FitFiles/");
	//if(FT==0) name1->Append("c3/");
	//else name1->Append("C3/");
	name1->Append("FitFile_CT2_FT");
	*name1 += ExpanType;
	name1->Append("_R");
	if(RT<RadiiCount-1) *name1 += RT;
	else *name1 += 100;
	name1->Append("_G");
	*name1 += Gindex;
	name1->Append(".root");
	infile = new TFile(name1->Data(),"READ");
	if(FT==0) fit = (TMinuit*)infile->Get("MyMinuit_c3");
	else fit = (TMinuit*)infile->Get("MyMinuit_C3");
	for(int parNum=0; parNum<6; parNum++){
	  fit->GetParameter(parNum+1, value,value_e);
	  ppEA[FT]->SetBinContent(RT+1, parNum+1, Gindex+1, value);
	}
	infile->Close();
      }
    }
  }
  // blank for the rest
  for(int FT=0; FT<2; FT++){// c3 or C3
    for(int Gindex=26; Gindex<50; Gindex++){
      for(int RT=0; RT<RadiiCount; RT++){// EW or LG
	for(int parNum=0; parNum<6; parNum++){
	  PbPbEA[FT]->SetBinContent(RT+1, parNum+1, Gindex+1, PbPbEA[FT]->GetBinContent(RT+1, parNum+1, 26));
	  pPbEA[FT]->SetBinContent(RT+1, parNum+1, Gindex+1, pPbEA[FT]->GetBinContent(RT+1, parNum+1, 26));
	  ppEA[FT]->SetBinContent(RT+1, parNum+1, Gindex+1, ppEA[FT]->GetBinContent(RT+1, parNum+1, 26));
	}
      }
    }
  }
  // Old way to convert Lam_3 to proper EA normalization
  /*for(int Gindex=0; Gindex<50; Gindex++){
    for(int RT=0; RT<RadiiCount; RT++){// EW or LG
      PbPbEA->SetBinContent(RT+1, 1, Gindex+1, pow(PbPbEA->GetBinContent(RT+1, 1, Gindex+1)/2., 1/3.));
      pPbEA->SetBinContent(RT+1, 1, Gindex+1, pow(pPbEA->GetBinContent(RT+1, 1, Gindex+1)/2., 1/3.));
      ppEA->SetBinContent(RT+1, 1, Gindex+1, pow(ppEA->GetBinContent(RT+1, 1, Gindex+1)/2., 1/3.));
    }
  }
  */

  TFile *outfile=new TFile("c3EAfile_temp.root","RECREATE");
  PbPbEA[0]->Write();
  pPbEA[0]->Write();
  ppEA[0]->Write();
  PbPbEA[1]->Write();
  pPbEA[1]->Write();
  ppEA[1]->Write();
  //
  outfile->Close();
  

}
