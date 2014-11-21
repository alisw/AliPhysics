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

  int FT=0;// 0(EW) or 1(LG)

  TFile *infile;
  TMinuit *fit;
  //
  TH3D *PbPbEA = new TH3D("PbPbEA","",6,0.5,6.5, 4,0.5,4.5, 50,-0.5,49.5);// Rcoh type, parNum, Gindex
  PbPbEA->SetDirectory(0);
  //
  TH3D *pPbEA = new TH3D("pPbEA","",6,0.5,6.5, 4,0.5,4.5, 50,-0.5,49.5);
  pPbEA->SetDirectory(0);
  //
  TH3D *ppEA = new TH3D("ppEA","",6,0.5,6.5, 4,0.5,4.5, 50,-0.5,49.5);
  ppEA->SetDirectory(0);
  //
 

  
  
  //
  //////////////////////////////
  double value=0, value_e=0;
  // 
  for(int Gindex=0; Gindex<=25; Gindex++){
    // PbPb
    for(int RT=0; RT<6; RT++){// Rcoh type
      
      TString *name1 = new TString("FitFiles/FitFile_CT0_FT");
      *name1 += FT;
      name1->Append("_R");
      *name1 += RT;
      name1->Append("_G");
      *name1 += Gindex;
      name1->Append(".root");
      infile = new TFile(name1->Data(),"READ");
      fit = (TMinuit*)infile->Get("MyMinuit_c3");
      for(int parNum=0; parNum<4; parNum++){
	fit->GetParameter(parNum+1, value,value_e);
	PbPbEA->SetBinContent(RT+1, parNum+1, Gindex+1, value);
      }
      infile->Close();
   
      //
      // pPb
      TString *name1 = new TString("FitFiles/FitFile_CT1_FT");
      *name1 += FT;
      name1->Append("_R");
      *name1 += RT;
      name1->Append("_G");
      *name1 += Gindex;
      name1->Append(".root");
      infile = new TFile(name1->Data(),"READ");
      fit = (TMinuit*)infile->Get("MyMinuit_c3");
      for(int parNum=0; parNum<4; parNum++){
	fit->GetParameter(parNum+1, value,value_e);
	pPbEA->SetBinContent(RT+1, parNum+1, Gindex+1, value);
      }
      infile->Close();
      //
      //
      TString *name1 = new TString("FitFiles/FitFile_CT2_FT");
      *name1 += FT;
      name1->Append("_R");
      *name1 += RT;
      name1->Append("_G");
      *name1 += Gindex;
      name1->Append(".root");
      infile = new TFile(name1->Data(),"READ");
      fit = (TMinuit*)infile->Get("MyMinuit_c3");
      for(int parNum=0; parNum<4; parNum++){
	fit->GetParameter(parNum+1, value,value_e);
	ppEA->SetBinContent(RT+1, parNum+1, Gindex+1, value);
      }
      infile->Close();
    }
  }
  // blank for the rest
  for(int Gindex=26; Gindex<50; Gindex++){
    for(int RT=0; RT<6; RT++){// EW or LG
      for(int parNum=0; parNum<4; parNum++){
	PbPbEA->SetBinContent(RT+1, parNum+1, Gindex+1, PbPbEA->GetBinContent(RT+1, parNum+1, 26));
	pPbEA->SetBinContent(RT+1, parNum+1, Gindex+1, pPbEA->GetBinContent(RT+1, parNum+1, 26));
	ppEA->SetBinContent(RT+1, parNum+1, Gindex+1, ppEA->GetBinContent(RT+1, parNum+1, 26));
      }
    }
  }

  // Convert Lam_3 to proper EA normalization
  for(int Gindex=0; Gindex<50; Gindex++){
    for(int RT=0; RT<6; RT++){// EW or LG
      PbPbEA->SetBinContent(RT+1, 1, Gindex+1, pow(PbPbEA->GetBinContent(RT+1, 1, Gindex+1), 1/3.));
      pPbEA->SetBinContent(RT+1, 1, Gindex+1, pow(pPbEA->GetBinContent(RT+1, 1, Gindex+1), 1/3.));
      ppEA->SetBinContent(RT+1, 1, Gindex+1, pow(ppEA->GetBinContent(RT+1, 1, Gindex+1), 1/3.));
    }
  }
  

  TFile *outfile=new TFile("c3EAfile_temp.root","RECREATE");
  PbPbEA->Write();
  pPbEA->Write();
  ppEA->Write();
  //
  outfile->Close();
  

}
