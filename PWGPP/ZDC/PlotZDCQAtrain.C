#if !defined(__CINT__) || defined(__MAKECINT__)

#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TClassTable.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TLine.h>
#include <TGrid.h>
#include <TBits.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFileMerger.h>
#include <TGridResult.h>
#include <TSystem.h>

#endif

void PlotZDCQAtrain(TString option="local",
		    Int_t nRun=0,
		    TString period="LHC11h",
		    TString recoPass="pass1_HLT",
		    TString qaTrain="QA",
		    TString fileName="QAresults.root"){

  gStyle->SetOptStat(111);
  //  gStyle->SetOptTitle(0);
  
  TFile *f;
  TString path;
  Int_t year=2011;
  if(period.Contains("LHC10")) year=2010;
  else if(period.Contains("LHC11")) year=2011;

  if(option.Contains("local")){
    char fnamerun[20];
    sprintf(fnamerun,"QAresults%d.root",nRun);
    f = new TFile(fnamerun);
    printf("Opened file %s\n",f->GetName());
  }
  else{
    TGrid::Connect("alien:");
    char filenameAlien[250];
    sprintf(filenameAlien,"alien:///alice/data/%d/%s/%09d/ESDs/%s/%s",year,period.Data(),nRun,recoPass.Data(),fileName.Data());    
    f = TFile::Open(filenameAlien);
    printf("Opened file %s\n",f->GetName());
  }

  TDirectoryFile* df=(TDirectoryFile*)f->Get("ZDC_Performance");
  if(!df){
    printf("ZDC_Performance MISSING -> Exit\n");
    return;
  }
  TList* l=(TList*)df->Get("QAZDCHists");
  if(!df){
    printf("QAZDCHists TList MISSING -> Exit\n");
    return;
  }

  TH2F *fDebunch       = (TH2F*)l->FindObject("fDebunch");
  TH1F *fhTDCZNSum     = (TH1F*)l->FindObject("fhTDCZNSum");	 
  TH1F *fhTDCZNDiff    = (TH1F*)l->FindObject("fhTDCZNDiff");	 
  TH1F *fhZNCSpectrum  = (TH1F*)l->FindObject("fhZNCSpectrum");  
  TH1F *fhZNASpectrum  = (TH1F*)l->FindObject("fhZNASpectrum");  
  TH1F *fhZPCSpectrum  = (TH1F*)l->FindObject("fhZPCSpectrum");  
  TH1F *fhZPASpectrum  = (TH1F*)l->FindObject("fhZPASpectrum");  
  TH1F *fhZEM1Spectrum = (TH1F*)l->FindObject("fhZEM1Spectrum"); 
  TH1F *fhZEM2Spectrum = (TH1F*)l->FindObject("fhZEM2Spectrum"); 
  TH1F *fhZNCpmc   = (TH1F*)l->FindObject("fhZNCpmc");   
  TH1F *fhZNApmc   = (TH1F*)l->FindObject("fhZNApmc");   
  TH1F *fhZPCpmc   = (TH1F*)l->FindObject("fhZPCpmc");   
  TH1F *fhZPApmc   = (TH1F*)l->FindObject("fhZPApmc");  
  TH2F *fhZNCCentroid  = (TH2F*)l->FindObject("fhZNCCentroid");  
  TH2F *fhZNACentroid  = (TH2F*)l->FindObject("fhZNACentroid");   
  TH1F *fhZNCemd   = (TH1F*)l->FindObject("fhZNCemd");   
  TH1F *fhZNAemd   = (TH1F*)l->FindObject("fhZNAemd");   
  TH1F *fhPMCZNCemd    = (TH1F*)l->FindObject("fhPMCZNCemd");	 
  TH1F *fhPMCZNAemd    = (TH1F*)l->FindObject("fhPMCZNAemd");	 
 
  /*TH1D *hxZNC = fhZNCCentroid->ProjectionX("hxZNC");
  TH1D *hyZNC = fhZNCCentroid->ProjectionY("hyZNC");
  TH1D *hxZNA = fhZNACentroid->ProjectionX("hxZNA");
  TH1D *hyZNA = fhZNACentroid->ProjectionY("hyZNA");*/


  gStyle->SetPalette(1);
  
  char nam[25];
  sprintf(nam,"QAhistos%d.root",nRun);
  TFile *fileout = TFile::Open(nam,"RECREATE");
        fDebunch->Write();
  	fhTDCZNSum->Write();    
	fhTDCZNDiff ->Write();   
	fhZNCSpectrum->Write();  
	fhZNASpectrum ->Write(); 
	fhZPCSpectrum ->Write(); 
	fhZPASpectrum ->Write(); 
	fhZEM1Spectrum->Write(); 
	fhZEM2Spectrum ->Write();
	fhZNCpmc->Write();
	fhZNApmc->Write();
	fhZPCpmc->Write();
	fhZPApmc->Write();
	fhZNCCentroid->Write();  
	fhZNACentroid->Write();  
	fhZNCemd->Write();
	fhZNAemd->Write();
	fhPMCZNCemd->Write();    
	fhPMCZNAemd->Write();
	fhZNCCentroid->Write();
  
  fileout->Close();
}

