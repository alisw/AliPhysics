#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <functional>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TTree.h>
#include <TSpectrum.h>

#include <TFile.h>
#include <TString.h>
#include <TList.h>

using namespace std;

void GetCorrectL1Phase(const Int_t runnumber, const TString trigger, Int_t *CorrectL1phase)
{
  TString infile = Form("AnalysisResults_%d.root",runnumber);
  cout << "infile = " << infile << endl;

  TFile *rootfile = TFile::Open(infile,"READ");
  if(!rootfile){
    printf("rootfile : %s does not exist. retrun.\n",infile.Data());
    return;
  }

  TString listname = Form("hist_%s",trigger.Data());
  TList *list = (TList*)rootfile->Get(listname);
  if(!list){
    printf("list : %s does not exist. retrun.\n",listname.Data());
    return;
  }

  TFile *outfile = new TFile(Form("Timing_%d.root",runnumber),"RECREATE");

  for(Int_t iddl=6;iddl<20;iddl++){
    TString histname = Form("hBC4vsTimeDDL%d_HG_afterCorr",iddl);
    //TString histname = Form("hBC4vsTimeDDL%d_HG",iddl);
    TH2F *h2 = (TH2F*)list->FindObject(histname);
    h2->SetXTitle("Time_{cell} (ns)");
    h2->SetYTitle("BC%4");
    outfile->WriteTObject(h2);

    Double_t mean[4]={};

    for(Int_t bc=0;bc<4;bc++){
      TH1F *h1 = (TH1F*)h2->ProjectionX(Form("h%d_BC%d",iddl,bc),bc+1,bc+1,"");
      outfile->WriteTObject(h1);
      mean[bc] = h1->GetMean();
    }

    Int_t index=-1;
    Double_t min=999;
    for(Int_t bc=0;bc<4;bc++){//find minimum value and index
      if(min>mean[bc]){
        min = mean[bc];
        index = bc;
      }
    }

    //cout << "DDL = " << iddl << " , index = " << index << " , min = " << min << endl;

    CorrectL1phase[iddl] = index;

  }


  rootfile->Close();
  outfile->Close();

}

