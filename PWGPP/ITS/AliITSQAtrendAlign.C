#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMath.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <Riostream.h>
#include "AliITSgeomTGeo.h"
#include <TLine.h>
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliITSCalibration.h"
#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliITSBadChannelsSSD.h"
#include "AliITSBadChannelsSSDv2.h"
#include "AliITSgeomTGeo.h"
#include <TObjArray.h>
#endif

TString pdfFileNames="";
void MakePlots(Int_t run1=-1,Int_t run2=99999999,TString ntupleFileName="TrendingITSAlign_norm.root");
void ReadOldSSDBadChannels(TObjArray *array, AliITSBadChannelsSSDv2 *badChannelsSSD);
void AliITSQAtrendAlign(TString runListFile="QA_Trend_list.txt",TString ntupleFileName="TrendingITSAlign_norm.root");

////////////////////////////////////////////////////////////////
//   Please, read this comment before using this macro 
//
//   INPUT FILE: a text file (by default QA_Trend_list.txt) which contains
//   a list of the complete path+file name of the QAresults root files 
//   (without the alien:// prefix).
//   One file per line. The order is irrelevant.
//
//   USAGE:
//
//   Function AliITSQAtrendAlign():  
//   it looks for a local root file named TrendingITSAlign_norm.root,
//   where the ntuples used to build the trending plots are
//   stored. This file is used to generate incrementally the ntuple contents:
//   when you add a new entry in the QA_Trend_list.txt file and you have a
//   local copy of the TrendingITSAlign_norm.root file, only the additional run will
//   be processed. The whole list is processed only the first time you use the
//   macro. 
//   For each run, the macro performs a gaussian fit of the xlocal and zlocal 
//   residuals from QA train of each module of each ITS layer; 
//   then evaluates the fractions of modules  
//   whose xlocal and zlocal residual mean and sigma are smaller than chosen 
//   threshold values; the fractions are normalized to the number of good modules 
//   directly read from OCDB fpor each run. 
//   Please, bear in mind that this macro is RAM-intensive: all the
//   ntuples are kept in memory. It is better to add few runs to the list at 
//   each time, according to the RAM of your computer. 
//   The function AliITSQAtrendAlign does not produce any plot.
//
//   Function MakePlots(run1,run2):
//   it produces the plots. For each canvas a PDF file is created.
//   A PDF file with all the canvases merged is also produced, named by default TrendingITSAlign_norm.root
//   The first two arguments define a range for the runs to be displayed
//   These two arguments are optional: by default, all the runs
//   found in the ntuples are displayed
////////////////////////////////////////////////////////////////

/* $Id$ */

//_____________________________________________________________________//
void ReadOldSSDBadChannels(TObjArray *array, 
						   AliITSBadChannelsSSDv2 *badChannelsSSD) {
	Int_t fNMod = array->GetEntries();
	cout<<"Converting old calibration object for bad channels..."<<endl;
	for (Int_t iModule = 0; iModule < fNMod; iModule++) {
		//for (Int_t iModule = 0; iModule < 1; iModule++) {
		AliITSBadChannelsSSD *bad = (AliITSBadChannelsSSD*) (array->At(iModule));
		TArrayI arrayPSide = bad->GetBadPChannelsList();
		for(Int_t iPCounter = 0; iPCounter < arrayPSide.GetSize(); iPCounter++) 
			badChannelsSSD->AddBadChannelP(iModule,
										   iPCounter,
										   (Char_t)arrayPSide.At(iPCounter));
        
		TArrayI arrayNSide = bad->GetBadNChannelsList();
		for(Int_t iNCounter = 0; iNCounter < arrayNSide.GetSize(); iNCounter++) 
			badChannelsSSD->AddBadChannelN(iModule,
										   iNCounter,
										   (Char_t)arrayNSide.At(iNCounter));
		
	}//loop over modules      
}

//_____________________________________________________________________//


void AliITSQAtrendAlign(TString runListFile,TString ntupleFileName){

  TGrid::Connect("alien:");

  const Int_t nVariables=169;
  TNtuple* ntits=new TNtuple("ntitsalign","ITS ALIGN trending","nrun:SDDfracModmux20_3:errSDDfracModmux20_3:SDDfracModmux50_3:errSDDfracModmux50_3:SDDfracModmux100_3:errSDDfracModmux100_3:SDDfracModsigx100_3:errSDDfracModsigx100_3:SDDfracModsigx200_3:errSDDfracModsigx200_3:SDDfracModsigx300_3:errSDDfracModsigx300_3:SDDfracModmuxlr20_3:errSDDfracModmuxlr20_3:SDDfracModmuxlr50_3:errSDDfracModmuxlr50_3:SDDfracModmuxlr100_3:errSDDfracModmuxlr100_3:SDDfracModexcx100_3:errSDDfracModexcx100_3:SDDfracModexcx200_3:errSDDfracModexcx200_3:SDDfracModexcx300_3:errSDDfracModexcx300_3:SDDfracModmuz50_3:errSDDfracModmuz50_3:SDDfracModmuz100_3:errSDDfracModmuz100_3:SDDfracModmuz300_3:errSDDfracModmuz300_3:SDDfracModsigz100_3:errSDDfracModsigz100_3:SDDfracModsigz300_3:errSDDfracModsigz300_3:SDDfracModsigz500_3:errSDDfracModsigz500_3:SDDfracModmux20_4:errSDDfracModmux20_4:SDDfracModmux50_4:errSDDfracModmux50_4:SDDfracModmux100_4:errSDDfracModmux100_4:SDDfracModsigx100_4:errSDDfracModsigx100_4:SDDfracModsigx200_4:errSDDfracModsigx200_4:SDDfracModsigx300_4:errSDDfracModsigx300_4:SDDfracModmuxlr20_4:errSDDfracModmuxlr20_4:SDDfracModmuxlr50_4:errSDDfracModmuxlr50_4:SDDfracModmuxlr100_4:errSDDfracModmuxlr100_4:SDDfracModexcx100_4:errSDDfracModexcx100_4:SDDfracModexcx200_4:errSDDfracModexcx200_4:SDDfracModexcx300_4:errSDDfracModexcx300_4:SDDfracModmuz50_4:errSDDfracModmuz50_4:SDDfracModmuz100_4:errSDDfracModmuz100_4:SDDfracModmuz300_4:errSDDfracModmuz300_4:SDDfracModsigz100_4:errSDDfracModsigz100_4:SDDfracModsigz300_4:errSDDfracModsigz300_4:SDDfracModsigz500_4:errSDDfracModsigz500_4:SPDfracModmux20_1:errSPDfracModmux20_1:SPDfracModmux50_1:errSPDfracModmux50_1:SPDfracModmux100_1:errSPDfracModmux100_1:SPDfracModsigx100_1:errSPDfracModsigx100_1:SPDfracModsigx200_1:errSPDfracModsigx200_1:SPDfracModsigx300_1:errSPDfracModsigx300_1:SPDfracModmuz50_1:errSPDfracModmuz50_1:SPDfracModmuz100_1:errSPDfracModmuz100_1:SPDfracModmuz300_1:errSPDfracModmuz300_1:SPDfracModsigz100_1:errSPDfracModsigz100_1:SPDfracModsigz300_1:errSPDfracModsigz300_1:SPDfracModsigz500_1:errSPDfracModsigz500_1:SPDfracModmux20_2:errSPDfracModmux20_2:SPDfracModmux50_2:errSPDfracModmux50_2:SPDfracModmux100_2:errSPDfracModmux100_2:SPDfracModsigx100_2:errSPDfracModsigx100_2:SPDfracModsigx200_2:errSPDfracModsigx200_2:SPDfracModsigx300_2:errSPDfracModsigx300_2:SPDfracModmuz50_2:errSPDfracModmuz50_2:SPDfracModmuz100_2:errSPDfracModmuz100_2:SPDfracModmuz300_2:errSPDfracModmuz300_2:SPDfracModsigz100_2:errSPDfracModsigz100_2:SPDfracModsigz300_2:errSPDfracModsigz300_2:SPDfracModsigz500_2:errSPDfracModsigz500_2:SSDfracModmux20_5:errSSDfracModmux20_5:SSDfracModmux50_5:errSSDfracModmux50_5:SSDfracModmux100_5:errSSDfracModmux100_5:SSDfracModsigx100_5:errSSDfracModsigx100_5:SSDfracModsigx200_5:errSSDfracModsigx200_5:SSDfracModsigx300_5:errSSDfracModsigx300_5:SSDfracModmuz50_5:errSSDfracModmuz50_5:SSDfracModmuz100_5:errSSDfracModmuz100_5:SSDfracModmuz300_5:errSSDfracModmuz300_5:SSDfracModsigz100_5:errSSDfracModsigz100_5:SSDfracModsigz300_5:errSSDfracModsigz300_5:SSDfracModsigz500_5:errSSDfracModsigz500_5:SSDfracModmux20_6:errSSDfracModmux20_6:SSDfracModmux50_6:errSSDfracModmux50_6:SSDfracModmux100_6:errSSDfracModmux100_6:SSDfracModsigx100_6:errSSDfracModsigx100_6:SSDfracModsigx200_6:errSSDfracModsigx200_6:SSDfracModsigx300_6:errSSDfracModsigx300_6:SSDfracModmuz50_6:errSSDfracModmuz50_6:SSDfracModmuz100_6:errSSDfracModmuz100_6:SSDfracModmuz300_6:errSSDfracModmuz300_6:SSDfracModsigz100_6:errSSDfracModsigz100_6:SSDfracModsigz300_6:errSSDfracModsigz300_6:SSDfracModsigz500_6:errSSDfracModsigz500_6");
  Float_t xnt[nVariables];

  // inizializzazione variabili OCDB
	
	AliCDBManager* man = AliCDBManager::Instance();
	man->SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
	AliCDBEntry *entrySDD;
	AliITSCalibration *cal = NULL;
	TObjArray *calSDD;
//
	
  TBits* readRun=new TBits(999999);
  readRun->ResetAllBits();

    if(!gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",ntupleFileName.Data()))){  // esiste gia' un file con l'ntupla ...
      TFile* oldfil=new TFile(ntupleFileName.Data());                             // che viene letta per aggiungere altre entries ....
      TNtuple* ntmp=(TNtuple*)oldfil->Get("ntitsalign");
		
      Bool_t isOK=kFALSE;
      if(ntmp){
	//	cout << "old ntu = " << ntmp->GetNvar() << ", new ntu  " <<  ntits->GetNvar() << endl;
	if(ntmp->GetNvar()==ntits->GetNvar()){
	  isOK=kTRUE;
	  TObjArray* arr1=(TObjArray*)ntits->GetListOfBranches();
	  TObjArray* arr2=(TObjArray*)ntmp->GetListOfBranches();
	  for(Int_t iV=0; iV<ntmp->GetNvar(); iV++){
	    TString vnam1=arr1->At(iV)->GetName();
	    TString vnam2=arr2->At(iV)->GetName();
	    if(vnam1!=vnam2) isOK=kFALSE;
	    ntmp->SetBranchAddress(vnam2.Data(),&xnt[iV]);
	  }
	  if(isOK){
	    for(Int_t nE=0; nE<ntmp->GetEntries(); nE++){
	      ntmp->GetEvent(nE);
	      Int_t theRun=(Int_t)(xnt[0]+0.0001);
	      readRun->SetBitNumber(theRun);
	      ntits->Fill(xnt);
	    }
	  }
	}
      }
      if(!isOK){
	printf("Ntuple in local file not OK -> will be recreated\n");
      }
      oldfil->Close();
      delete oldfil;
    }
	
#define MAX_LINES 200
#define MAX_LINE_LEN 255
	
	char strings[MAX_LINES][MAX_LINE_LEN];
	ifstream in(runListFile.Data());
	int j = 0;
	Int_t nrun=0;
	Int_t runNumb[MAX_LINES];
	Bool_t goout = kFALSE;
	while ( in ) {
		in.getline(strings[j], MAX_LINE_LEN);
		TString aux(strings[j]);
		Int_t lentrail=0;
		if(aux.Contains("LHC11h/") || aux.Contains("LHC11e/")){
			lentrail = 27;
		}
		else if(aux.Contains("LHC11h_2/")){
			lentrail = 29;
		}
		else {
			if(!aux.IsNull())printf("Unrecognised path name %s \n",aux.Data());
			goout = kTRUE;
		}
		if(goout)break;
		if(aux.Length()<lentrail)continue;
		aux=aux.Remove(0,lentrail);
		aux=aux.Remove(6,aux.Length());  
		runNumb[j]=atoi(aux.Data());
		printf("%d ) - path %s \n",runNumb[j],strings[j]);
		j++;
		nrun++;
	}
	
	// loop on runs
	
	printf("\n *******************   Loop on runs *********** \n");

	Int_t twice=0;

	for(Int_t jru=0;jru<nrun;jru++) {
		printf("jru=%d - run number= %d \n",jru,runNumb[jru]);
		Int_t iRun=runNumb[jru];
		if(readRun->TestBitNumber(iRun))printf("Run %d - already processed\n",iRun);
		if(readRun->TestBitNumber(iRun))continue;
		//cout << "Value from file is " <<t << endl;

		if(jru>0){
			for(Int_t jru2=0;jru2<jru;jru2++) {
				if(runNumb[jru2]==runNumb[jru]){
					printf("Run %d - twice in the list - skipped\n",runNumb[jru]);
					twice=1;
				}
			}
		}
		if(twice==1){twice=0;continue;}
		
		printf("%s\n",strings[jru]);
		
		if(!gGrid||!gGrid->IsConnected()) {
    printf("gGrid not found! exit macro\n");
    return;
		}
		TFile *f=TFile::Open(Form("alien://%s",strings[jru])); 

      TDirectoryFile *df=(TDirectoryFile*)f->Get("ITSAlignQA");
      if(!df){
	printf("Run %d ITSAlignQA MISSING -> Exit\n",iRun);
	continue;
      }
      TList *l = (TList*)df->Get("clistITSAlignQA");
      if(!df){
	printf("Run %d clistITSAlignQA TList MISSING -> Exit\n",iRun);
	continue;
      }

      TH1F* hev=(TH1F*)l->FindObject("hNEvents");  // histo 1d numero di eventi processati
      Int_t nTotEvents=hev->GetBinContent(2);
      Int_t nTrigEvents=hev->GetBinContent(3);
       printf("Run %d Number of Events = %d Triggered=%d\n",iRun,nTotEvents,nTrigEvents);
      if(nTotEvents==0) continue;

      // working variables 
      Int_t nModmux20_1=0; 
      Int_t nModmux50_1=0;
      Int_t nModmux100_1=0;
      Int_t nModsigx100_1=0; 
      Int_t nModsigx200_1=0;
      Int_t nModsigx300_1=0;
      Int_t nModmuz50_1=0;
      Int_t nModmuz100_1=0;
      Int_t nModmuz300_1=0;
      Int_t nModsigz100_1=0;
      Int_t nModsigz300_1=0;
      Int_t nModsigz500_1=0;

      Int_t nModmux20_2=0; 
      Int_t nModmux50_2=0;
      Int_t nModmux100_2=0;
      Int_t nModsigx100_2=0; 
      Int_t nModsigx200_2=0;
      Int_t nModsigx300_2=0;
      Int_t nModmuz50_2=0;
      Int_t nModmuz100_2=0;
      Int_t nModmuz300_2=0;
      Int_t nModsigz100_2=0;
      Int_t nModsigz300_2=0;
      Int_t nModsigz500_2=0;

      Int_t nModmux20_3=0; 
      Int_t nModmux50_3=0;
      Int_t nModmux100_3=0;
      Int_t nModsigx100_3=0; 
      Int_t nModsigx200_3=0;
      Int_t nModsigx300_3=0;
      Int_t nModmuxlr20_3=0;
      Int_t nModmuxlr50_3=0;
      Int_t nModmuxlr100_3=0;
      Int_t nModexcx100_3=0;
      Int_t nModexcx200_3=0;
      Int_t nModexcx300_3=0;
      Int_t nModmuz50_3=0;
      Int_t nModmuz100_3=0;
      Int_t nModmuz300_3=0;
      Int_t nModsigz100_3=0;
      Int_t nModsigz300_3=0;
      Int_t nModsigz500_3=0;

      Int_t nModmux20_4=0; 
      Int_t nModmux50_4=0;
      Int_t nModmux100_4=0;
      Int_t nModsigx100_4=0;
      Int_t nModsigx200_4=0;
      Int_t nModsigx300_4=0;
      Int_t nModmuxlr20_4=0;
      Int_t nModmuxlr50_4=0;
      Int_t nModmuxlr100_4=0;
      Int_t nModexcx100_4=0;
      Int_t nModexcx200_4=0;
      Int_t nModexcx300_4=0;
      Int_t nModmuz50_4=0;
      Int_t nModmuz100_4=0;
      Int_t nModmuz300_4=0;
      Int_t nModsigz100_4=0;
      Int_t nModsigz300_4=0;
      Int_t nModsigz500_4=0;

      Int_t nModmux20_5=0; 
      Int_t nModmux50_5=0;
      Int_t nModmux100_5=0;
      Int_t nModsigx100_5=0;
      Int_t nModsigx200_5=0;
      Int_t nModsigx300_5=0;
      Int_t nModmuz50_5=0;
      Int_t nModmuz100_5=0;
      Int_t nModmuz300_5=0;
      Int_t nModsigz100_5=0;
      Int_t nModsigz300_5=0;
      Int_t nModsigz500_5=0;

      Int_t nModmux20_6=0; 
      Int_t nModmux50_6=0;
      Int_t nModmux100_6=0;
      Int_t nModsigx100_6=0;
      Int_t nModsigx200_6=0;
      Int_t nModsigx300_6=0;
      Int_t nModmuz50_6=0;
      Int_t nModmuz100_6=0;
      Int_t nModmuz300_6=0;
      Int_t nModsigz100_6=0;
      Int_t nModsigz300_6=0;
      Int_t nModsigz500_6=0;
   
      // mio pezzo lettura e fit

    for(Int_t nmod=0; nmod<240; nmod++){  // loop on SPD modules
      
      TString nameResXvsX=Form("hSPDResidX%d",nmod);
      TH2F *histoResXvsX=(TH2F*)l->FindObject(nameResXvsX);
      TH1D *histoResX=histoResXvsX->ProjectionY();

      TString nameResZvsZ=Form("hSPDResidZ%d",nmod);
      TH2F *histoResZvsZ=(TH2F*)l->FindObject(nameResZvsZ);
      TH1D *histoResZ=histoResZvsZ->ProjectionY();

      TF1 *gg_x = new TF1("gg_x","gaus",-0.050, 0.050);
      TF1 *gg_z = new TF1("gg_z","gaus",-0.15,0.15);

      if(histoResX->GetEntries()){
	histoResX->Fit("gg_x","RQ0");
	if(nmod<80){  // layer 1
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<5.) nModmux20_1+=1; // micron
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<10.) nModmux50_1+=1; 
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<50.) nModmux100_1+=1; 
	  if(gg_x->GetParameter(2)*10000<50.) nModsigx100_1+=1;  
	  if(gg_x->GetParameter(2)*10000<70.) nModsigx200_1+=1;  
	  if(gg_x->GetParameter(2)*10000<200.) nModsigx300_1+=1;  
	}
	else { // layer 2
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<10.) nModmux20_2+=1; // micron
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<20.) nModmux50_2+=1; 
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<50.) nModmux100_2+=1; 
	  if(gg_x->GetParameter(2)*10000<130.) nModsigx100_2+=1;  
	  if(gg_x->GetParameter(2)*10000<160.) nModsigx200_2+=1;  
	  if(gg_x->GetParameter(2)*10000<200.) nModsigx300_2+=1;  
	}
      }
  
      
      if(histoResZ->GetEntries()){
	histoResZ->Fit("gg_z","RQ0");
	if(nmod<80){  // layer 1
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<10.) nModmuz50_1+=1;
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<20.) nModmuz100_1+=1;
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<50.) nModmuz300_1+=1;
	  if(gg_z->GetParameter(2)*10000<150.) nModsigz100_1+=1;  
	  if(gg_z->GetParameter(2)*10000<160.) nModsigz300_1+=1;  
	  if(gg_z->GetParameter(2)*10000<200.) nModsigz500_1+=1;  
	}

	else { // layer 2
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<10.) nModmuz50_2+=1;
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<20.) nModmuz100_2+=1;
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<50.) nModmuz300_2+=1;
	  if(gg_z->GetParameter(2)*10000<180.) nModsigz100_2+=1;  
	  if(gg_z->GetParameter(2)*10000<190.) nModsigz300_2+=1;  
	  if(gg_z->GetParameter(2)*10000<250.) nModsigz500_2+=1;  
	}
      }

    } // loop on SPD modules -- end 

    for(Int_t nmod=240; nmod<500; nmod++){  // loop on SDD modules
      
      TString nameResXvsX=Form("hSDDResidXvsX%d",nmod);
      TH2F *histoResXvsX=(TH2F*)l->FindObject(nameResXvsX);
      TH1D *histoResX=histoResXvsX->ProjectionY();
      TH1D *histoResX_l=histoResXvsX->ProjectionY("residual xloc<0",1,20);
      TH1D *histoResX_r=histoResXvsX->ProjectionY("residual xloc>0",21,40);

      TString nameResZvsZ=Form("hSDDResidZvsZ%d",nmod);
      TH2F *histoResZvsZ=(TH2F*)l->FindObject(nameResZvsZ);
      TH1D *histoResZ=histoResZvsZ->ProjectionY();

      //      TString nameResZvsX=Form("hSDDResidZvsX%d",nmod);
      //      TH2F *histoResZvsX=(TH2F*)l->FindObject(nameResZvsX);

      TF1 *gg_x = new TF1("gg_x","gaus",-0.050, 0.050);
      TF1 *gg_l = new TF1("gg_l","gaus",-0.050, 0.050);
      TF1 *gg_r = new TF1("gg_r","gaus",-0.050, 0.050);
      TF1 *gg_z = new TF1("gg_z","gaus",-0.15,0.15);

      if(histoResX->GetEntries()){
	histoResX->Fit("gg_x","RQ0");
	if(nmod<324){  // layer 3
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<20.) nModmux20_3+=1; // micron
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<50.) nModmux50_3+=1; 
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<100.) nModmux100_3+=1; 
	  if(gg_x->GetParameter(2)*10000<200.) nModsigx100_3+=1;  
	  if(gg_x->GetParameter(2)*10000<250.) nModsigx200_3+=1;  
	  if(gg_x->GetParameter(2)*10000<400.) nModsigx300_3+=1;  
	}
	else { // layer 4
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<20.) nModmux20_4+=1; // micron
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<50.) nModmux50_4+=1; 
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<100.) nModmux100_4+=1; 
	  if(gg_x->GetParameter(2)*10000<200.) nModsigx100_4+=1;  
	  if(gg_x->GetParameter(2)*10000<250.) nModsigx200_4+=1;  
	  if(gg_x->GetParameter(2)*10000<400.) nModsigx300_4+=1;  
	}
      }
  
      if(histoResX_l->GetEntries()&&histoResX_r->GetEntries()){
	histoResX_l->Fit("gg_l","RQ0");
	histoResX_r->Fit("gg_r","RQ0");
	if(nmod<324){  // layer 3
	  if(TMath::Abs(gg_l->GetParameter(1)*10000-gg_r->GetParameter(1)*10000)<20.) nModmuxlr20_3+=1;
	  if(TMath::Abs(gg_l->GetParameter(1)*10000-gg_r->GetParameter(1)*10000)<50.) nModmuxlr50_3+=1;
	  if(TMath::Abs(gg_l->GetParameter(1)*10000-gg_r->GetParameter(1)*10000)<100.) nModmuxlr100_3+=1;
	}
	else { // layer 4
	  if(TMath::Abs(gg_l->GetParameter(1)*10000-gg_r->GetParameter(1)*10000)<20.) nModmuxlr20_4+=1;
	  if(TMath::Abs(gg_l->GetParameter(1)*10000-gg_r->GetParameter(1)*10000)<50.) nModmuxlr50_4+=1;
	  if(TMath::Abs(gg_l->GetParameter(1)*10000-gg_r->GetParameter(1)*10000)<100.) nModmuxlr100_4+=1;
	}
      }
      
      if(histoResZ->GetEntries()){
	histoResZ->Fit("gg_z","RQ0");
	if(nmod<324){  // layer 3
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<30.) nModmuz50_3+=1;
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<50.) nModmuz100_3+=1;
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<100.) nModmuz300_3+=1;
	  if(gg_z->GetParameter(2)*10000<150.) nModsigz100_3+=1;  
	  if(gg_z->GetParameter(2)*10000<200.) nModsigz300_3+=1;  
	  if(gg_z->GetParameter(2)*10000<350.) nModsigz500_3+=1;  
	}

	else { // layer 4
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<30.) nModmuz50_4+=1;
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<50.) nModmuz100_4+=1;
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<100.) nModmuz300_4+=1;
	  if(gg_z->GetParameter(2)*10000<250.) nModsigz100_4+=1;  
	  if(gg_z->GetParameter(2)*10000<300.) nModsigz300_4+=1;  
	  if(gg_z->GetParameter(2)*10000<450.) nModsigz500_4+=1;  
	}
      }


      Double_t min=-10000.0;
      Double_t max=10000.0;
      
	TProfile *profResXvsX = histoResXvsX->ProfileX();
	if(profResXvsX->GetEntries()>0){
	  max=profResXvsX->GetBinContent(1);
	  min=profResXvsX->GetBinContent(1);
	  for(Int_t jj=0;jj<39;jj++){
	    if(profResXvsX->GetBinContent(jj+1)>max) max = profResXvsX->GetBinContent(jj+1);
	    if(profResXvsX->GetBinContent(jj+1)<min) min = profResXvsX->GetBinContent(jj+1);
	  }	
	  
	  if(nmod<324){  // layer 3
	    if(TMath::Abs(max-min)*10000.<100) nModexcx100_3+=1;
	    if(TMath::Abs(max-min)*10000.<200) nModexcx200_3+=1;
	    if(TMath::Abs(max-min)*10000.<300) nModexcx300_3+=1;
	  }
	  
	  else { // layer 4
	    if(TMath::Abs(max-min)*10000.<100) nModexcx100_4+=1;
	    if(TMath::Abs(max-min)*10000.<200) nModexcx200_4+=1;
	    if(TMath::Abs(max-min)*10000.<300) nModexcx300_4+=1;
	  }
	}

    } // loop on SDD modules -- end

    for(Int_t nmod=500; nmod<2198; nmod++){  // loop on SSD modules
      
      TString nameResXvsX=Form("hSSDResidX%d",nmod);
      TH2F *histoResXvsX=(TH2F*)l->FindObject(nameResXvsX);
      TH1D *histoResX=histoResXvsX->ProjectionY();

      TString nameResZvsZ=Form("hSSDResidZ%d",nmod);
      TH2F *histoResZvsZ=(TH2F*)l->FindObject(nameResZvsZ);
      TH1D *histoResZ=histoResZvsZ->ProjectionY();

      TF1 *gg_x = new TF1("gg_x","gaus",-0.050, 0.050);
      TF1 *gg_z = new TF1("gg_z","gaus",-0.15,0.15);

      if(histoResX->GetEntries()){
	histoResX->Fit("gg_x","RQ0");
	if(nmod<1248){  // layer 5
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<20.) nModmux20_5+=1; // micron
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<40.) nModmux50_5+=1; 
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<100.) nModmux100_5+=1; 
	  if(gg_x->GetParameter(2)*10000<70.) nModsigx100_5+=1;  
	  if(gg_x->GetParameter(2)*10000<100.) nModsigx200_5+=1;  
	  if(gg_x->GetParameter(2)*10000<300.) nModsigx300_5+=1;  
	}
	else { // layer 6
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<20.) nModmux20_6+=1; // micron
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<40.) nModmux50_6+=1; 
	  if(TMath::Abs(gg_x->GetParameter(1))*10000<100.) nModmux100_6+=1; 
	  if(gg_x->GetParameter(2)*10000<200.) nModsigx100_6+=1;  
	  if(gg_x->GetParameter(2)*10000<300.) nModsigx200_6+=1;  
	  if(gg_x->GetParameter(2)*10000<400.) nModsigx300_6+=1;  
	}
      }
  
      
      if(histoResZ->GetEntries()){
	histoResZ->Fit("gg_z","RQ0");
	if(nmod<1248){  // layer 5
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<50.) nModmuz50_5+=1;
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<100.) nModmuz100_5+=1;
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<300.) nModmuz300_5+=1;
	  if(gg_z->GetParameter(2)*10000<1000.) nModsigz100_5+=1;  
	  if(gg_z->GetParameter(2)*10000<1100.) nModsigz300_5+=1;  
	  if(gg_z->GetParameter(2)*10000<1500.) nModsigz500_5+=1;  
	}

	else { // layer 6
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<50.) nModmuz50_6+=1;
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<100.) nModmuz100_6+=1;
	  if(TMath::Abs(gg_z->GetParameter(1))*10000<300.) nModmuz300_6+=1;
	  if(gg_z->GetParameter(2)*10000<1000.) nModsigz100_6+=1;  
	  if(gg_z->GetParameter(2)*10000<1100.) nModsigz300_6+=1;  
	  if(gg_z->GetParameter(2)*10000<1500.) nModsigz500_6+=1;  
	}
      }

    } // loop on SSD modules -- end

		Int_t count1 = 0;
		Int_t count2 = 0;
		Int_t count3 = 0;
		Int_t count4 = 0;
		Int_t count5 = 0;
		Int_t count6 = 0;
		
//		Int_t nPSideChannelsTotal = 0, nNSideChannelsTotal = 0;
//		Int_t nBadPSideChannelsTotal = 0, nBadNSideChannelsTotal = 0;
		Int_t nBadPSideChannels = 0, nBadNSideChannels = 0;
		Int_t layer = 0, ladder = 0, module = 0;
//		Int_t nPSideChannelsLayer5 = 0, nNSideChannelsLayer5 = 0;
//		Int_t nPSideChannelsLayer6 = 0, nNSideChannelsLayer6 = 0;
//		Int_t nPSideChannelsLayer5Total = 0, nNSideChannelsLayer5Total = 0;
//		Int_t nPSideChannelsLayer6Total = 0, nNSideChannelsLayer6Total = 0;
//		Int_t nBadPSideChannelsLayer5Total = 0,  nBadNSideChannelsLayer5Total = 0;
//		Int_t nBadPSideChannelsLayer6Total = 0,  nBadNSideChannelsLayer6Total = 0;
//		Int_t badPmodule5 = 0, badNmodule5 = 0, badPmodule6 = 0, badNmodule6 = 0;
		
		// lettura moduli SDD Ok da OCDB
		man->SetRun(iRun);
		entrySDD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSDD");
		calSDD = (TObjArray *)entrySDD->GetObject();
		for(Int_t i=0;i<84;i++){cal = (AliITSCalibration*) 
			calSDD->At(i); if(cal->IsBad())count3++;}
		for(Int_t i=84;i<calSDD->GetEntries();i++){cal = (AliITSCalibration*) 
			calSDD->At(i); if(cal->IsBad())count4++;}
//		cout << endl << endl << count3 << " " << count4 << endl;
		
    // ntuple filling
    Int_t index=0;
    // 0
    xnt[index++]=(Float_t)iRun;
    // 1 SDD L3
    xnt[index++]=nModmux20_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModmux20_3)/(84.-count3);
    xnt[index++]=nModmux50_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModmux50_3)/(84.-count3);
    xnt[index++]=nModmux100_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModmux100_3)/(84.-count3);
    xnt[index++]=nModsigx100_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModsigx100_3)/(84.-count3);
    xnt[index++]=nModsigx200_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModsigx200_3)/(84.-count3);
    xnt[index++]=nModsigx300_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModsigx300_3)/(84.-count3);
    xnt[index++]=nModmuxlr20_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModmuxlr20_3)/(84.-count3);
    xnt[index++]=nModmuxlr50_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModmuxlr50_3)/(84.-count3);
    xnt[index++]=nModmuxlr100_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModmuxlr100_3)/(84.-count3);
    xnt[index++]=nModexcx100_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModexcx100_3)/(84.-count3);
    xnt[index++]=nModexcx200_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModexcx200_3)/(84.-count3);
    xnt[index++]=nModexcx300_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModexcx300_3)/(84.-count3);
    xnt[index++]=nModmuz50_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModmuz50_3)/(84.-count3);
    xnt[index++]=nModmuz100_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModmuz100_3)/(84.-count3);
    xnt[index++]=nModmuz300_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModmuz300_3)/(84.-count3);
    xnt[index++]=nModsigz100_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModsigz100_3)/(84.-count3);
    xnt[index++]=nModsigz300_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModsigz300_3)/(84.-count3);
    xnt[index++]=nModsigz500_3/(84.-count3);
    xnt[index++]=TMath::Sqrt(nModsigz500_3)/(84.-count3);
    // 37 SDD L4
    xnt[index++]=nModmux20_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModmux20_4)/(176.-count4);
    xnt[index++]=nModmux50_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModmux50_4)/(176.-count4);
    xnt[index++]=nModmux100_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModmux100_4)/(176.-count4);
    xnt[index++]=nModsigx100_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModsigx100_4)/(176.-count4);
    xnt[index++]=nModsigx200_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModsigx200_4)/(176.-count4);
    xnt[index++]=nModsigx300_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModsigx300_4)/(176.-count4);
    xnt[index++]=nModmuxlr20_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModmuxlr20_4)/(176.-count4);
    xnt[index++]=nModmuxlr50_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModmuxlr50_4)/(176.-count4);
    xnt[index++]=nModmuxlr100_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModmuxlr100_4)/(176.-count4);
    xnt[index++]=nModexcx100_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModexcx100_4)/(176.-count4);
    xnt[index++]=nModexcx200_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModexcx200_4)/(176.-count4);
    xnt[index++]=nModexcx300_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModexcx300_4)/(176.-count4);
    xnt[index++]=nModmuz50_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModmuz50_4)/(176.-count4);
    xnt[index++]=nModmuz100_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModmuz100_4)/(176.-count4);
    xnt[index++]=nModmuz300_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModmuz300_4)/(176.-count4);
    xnt[index++]=nModsigz100_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModsigz100_4)/(176.-count4);
    xnt[index++]=nModsigz300_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModsigz300_4)/(176.-count4);
    xnt[index++]=nModsigz500_4/(176.-count4);
    xnt[index++]=TMath::Sqrt(nModsigz500_4)/(176.-count4);
    // 73 SPD L1
		// lettura moduli SPD bad da OCDB		
		AliITSOnlineCalibrationSPDhandler *h = new AliITSOnlineCalibrationSPDhandler();
		h->ReadDeadFromDB(iRun,"alien://folder=/alice/data/2011/OCDB");
//		Int_t nDeadModules[2]={0,0}; // n dead inner and outer
//		Double_t fractions[2]={0.,0.}; // fraction of dead modules
		// # dead modules inner layer
		for(Int_t i=0; i<80; i++) count1+=h->GetNrBad(i)/(5*8192);
		// # dead modules outer layer
		for(Int_t jj=80; jj<240; jj++) count2+=h->GetNrBad(jj)/(5*8192);
		
    xnt[index++]=nModmux20_1/(80.-count1);
    xnt[index++]=TMath::Sqrt(nModmux20_1)/(80.-count1);
    xnt[index++]=nModmux50_1/(80.-count1);
    xnt[index++]=TMath::Sqrt(nModmux50_1)/(80.-count1);
    xnt[index++]=nModmux100_1/(80.-count1);
    xnt[index++]=TMath::Sqrt(nModmux100_1)/(80.-count1);
    xnt[index++]=nModsigx100_1/(80.-count1);
    xnt[index++]=TMath::Sqrt(nModsigx100_1)/(80.-count1);
    xnt[index++]=nModsigx200_1/(80.-count1);
    xnt[index++]=TMath::Sqrt(nModsigx200_1)/(80.-count1);
    xnt[index++]=nModsigx300_1/(80.-count1);
    xnt[index++]=TMath::Sqrt(nModsigx300_1)/(80.-count1);
    xnt[index++]=nModmuz50_1/(80.-count1);
    xnt[index++]=TMath::Sqrt(nModmuz50_1)/(80.-count1);
    xnt[index++]=nModmuz100_1/(80.-count1);
    xnt[index++]=TMath::Sqrt(nModmuz100_1)/(80.-count1);
    xnt[index++]=nModmuz300_1/(80.-count1);
    xnt[index++]=TMath::Sqrt(nModmuz300_1)/(80.-count1);
    xnt[index++]=nModsigz100_1/(80.-count1);
    xnt[index++]=TMath::Sqrt(nModsigz100_1)/(80.-count1);
    xnt[index++]=nModsigz300_1/(80.-count1);
    xnt[index++]=TMath::Sqrt(nModsigz300_1)/(80.-count1);
    xnt[index++]=nModsigz500_1/(80.-count1);
    xnt[index++]=TMath::Sqrt(nModsigz500_1)/(80.-count1);
    // 97 SPD L2
    xnt[index++]=nModmux20_2/(160.-count2);
    xnt[index++]=TMath::Sqrt(nModmux20_2)/(160.-count2);
    xnt[index++]=nModmux50_2/(160.-count2);
    xnt[index++]=TMath::Sqrt(nModmux50_2)/(160.-count2);
    xnt[index++]=nModmux100_2/(160.-count2);
    xnt[index++]=TMath::Sqrt(nModmux100_2)/(160.-count2);
    xnt[index++]=nModsigx100_2/(160.-count2);
    xnt[index++]=TMath::Sqrt(nModsigx100_2)/(160.-count2);
    xnt[index++]=nModsigx200_2/(160.-count2);
    xnt[index++]=TMath::Sqrt(nModsigx200_2)/(160.-count2);
    xnt[index++]=nModsigx300_2/(160.-count2);
    xnt[index++]=TMath::Sqrt(nModsigx300_2)/(160.-count2);
    xnt[index++]=nModmuz50_2/(160.-count2);
    xnt[index++]=TMath::Sqrt(nModmuz50_2)/(160.-count2);
    xnt[index++]=nModmuz100_2/(160.-count2);
    xnt[index++]=TMath::Sqrt(nModmuz100_2)/(160.-count2);
    xnt[index++]=nModmuz300_2/(160.-count2);
    xnt[index++]=TMath::Sqrt(nModmuz300_2)/(160.-count2);
    xnt[index++]=nModsigz100_2/(160.-count2);
    xnt[index++]=TMath::Sqrt(nModsigz100_2)/(160.-count2);
    xnt[index++]=nModsigz300_2/(160.-count2);
    xnt[index++]=TMath::Sqrt(nModsigz300_2)/(160.-count2);
    xnt[index++]=nModsigz500_2/(160.-count2);
    xnt[index++]=TMath::Sqrt(nModsigz500_2)/(160.-count2);
    // 121 SSD L5
		// lettura moduli SSD bad da OCDB
		const Int_t fgkSSDMODULES = 1698;
		static const Int_t fgkDefaultNStripsSSD = 768;
		AliITSBadChannelsSSDv2 *badChannelsSSD = new AliITSBadChannelsSSDv2();
		AliCDBEntry *entryBadChannelsSSD = man->Get("ITS/Calib/BadChannelsSSD");
		TObject *empty = (TObject *)entryBadChannelsSSD->GetObject();
		TString objectname = empty->GetName();
		if(objectname=="TObjArray") {
			TObjArray *badChannelsSSDOld = (TObjArray *)entryBadChannelsSSD->GetObject();
			ReadOldSSDBadChannels(badChannelsSSDOld, badChannelsSSD);
		}	
		else if(objectname=="AliITSBadChannelsSSDv2") {
			cout<<"Reading the new format of the calibration file..."<<endl;
			badChannelsSSD = (AliITSBadChannelsSSDv2 *)entryBadChannelsSSD->GetObject();
		}
		nBadPSideChannels = 0; 
		nBadNSideChannels = 0;
		layer = 0; 
		ladder = 0; 
		module = 0;
		for(Int_t i = 0; i < fgkSSDMODULES; i++) { // loop sui moduli
			AliITSgeomTGeo::GetModuleId(i+500,layer,ladder,module);
			nBadPSideChannels = 0, nBadNSideChannels = 0;
			
			Int_t badChannel = 0;
			for(Int_t ja = 0; ja < fgkDefaultNStripsSSD; ja++) {   // loop sulle strips del modulo 
				badChannel = (Int_t)(badChannelsSSD->GetBadChannelP(i,ja));   // bad channels P side
				//cout<<"Module: "<<i+500<< " Strip: "<<j<<" - "<<badChannel<<endl;
				if(badChannel != 0) nBadPSideChannels += 1;
				
				badChannel = (Int_t)(badChannelsSSD->GetBadChannelN(i,ja));  // bad channels N side
				//cout<<"Module: "<<i+500<< " Strip: "<<fgkDefaultNStripsSSD+j+1<<" - "<<badChannel<<endl;
				if(badChannel != 0) nBadNSideChannels += 1;
			} // end loop sulle strip del modulo
			
			// determinazione numero moduli morti sui due lati per ogni layer 	
			if(nBadNSideChannels == 768  && nBadPSideChannels == 768){
				if(layer == 5){
					count5++;
				}
				else {
					count6++;
				}
			}
		} // endl loop sui moduli
			
    xnt[index++]=nModmux20_5/(748.-count5);
    xnt[index++]=TMath::Sqrt(nModmux20_5)/(748.-count5);
    xnt[index++]=nModmux50_5/(748.-count5);
    xnt[index++]=TMath::Sqrt(nModmux50_5)/(748.-count5);
    xnt[index++]=nModmux100_5/(748.-count5);
    xnt[index++]=TMath::Sqrt(nModmux100_5)/(748.-count5);
    xnt[index++]=nModsigx100_5/(748.-count5);
    xnt[index++]=TMath::Sqrt(nModsigx100_5)/(748.-count5);
    xnt[index++]=nModsigx200_5/(748.-count5);
    xnt[index++]=TMath::Sqrt(nModsigx200_5)/(748.-count5);
    xnt[index++]=nModsigx300_5/(748.-count5);
    xnt[index++]=TMath::Sqrt(nModsigx300_5)/(748.-count5);
    xnt[index++]=nModmuz50_5/(748.-count5);
    xnt[index++]=TMath::Sqrt(nModmuz50_5)/(748.-count5);
    xnt[index++]=nModmuz100_5/(748.-count5);
    xnt[index++]=TMath::Sqrt(nModmuz100_5)/(748.-count5);
    xnt[index++]=nModmuz300_5/(748.-count5);
    xnt[index++]=TMath::Sqrt(nModmuz300_5)/(748.-count5);
    xnt[index++]=nModsigz100_5/(748.-count5);
    xnt[index++]=TMath::Sqrt(nModsigz100_5)/(748.-count5);
    xnt[index++]=nModsigz300_5/(748.-count5);
    xnt[index++]=TMath::Sqrt(nModsigz300_5)/(748.-count5);
    xnt[index++]=nModsigz500_5/(748.-count5);
    xnt[index++]=TMath::Sqrt(nModsigz500_5)/(748.-count5);
    // 145 SSD L6
    xnt[index++]=nModmux20_6/(950.-count6);
    xnt[index++]=TMath::Sqrt(nModmux20_6)/(950.-count6);
    xnt[index++]=nModmux50_6/(950.-count6);
    xnt[index++]=TMath::Sqrt(nModmux50_6)/(950.-count6);
    xnt[index++]=nModmux100_6/(950.-count6);
    xnt[index++]=TMath::Sqrt(nModmux100_6)/(950.-count6);
    xnt[index++]=nModsigx100_6/(950.-count6);
    xnt[index++]=TMath::Sqrt(nModsigx100_6)/(950.-count6);
    xnt[index++]=nModsigx200_6/(950.-count6);
    xnt[index++]=TMath::Sqrt(nModsigx200_6)/(950.-count6);
    xnt[index++]=nModsigx300_6/(950.-count6);
    xnt[index++]=TMath::Sqrt(nModsigx300_6)/(950.-count6);
    xnt[index++]=nModmuz50_6/(950.-count6);
    xnt[index++]=TMath::Sqrt(nModmuz50_6)/(950.-count6);
    xnt[index++]=nModmuz100_6/(950.-count6);
    xnt[index++]=TMath::Sqrt(nModmuz100_6)/(950.-count6);
    xnt[index++]=nModmuz300_6/(950.-count6);
    xnt[index++]=TMath::Sqrt(nModmuz300_6)/(950.-count6);
    xnt[index++]=nModsigz100_6/(950.-count6);
    xnt[index++]=TMath::Sqrt(nModsigz100_6)/(950.-count6);
    xnt[index++]=nModsigz300_6/(950.-count6);
    xnt[index++]=TMath::Sqrt(nModsigz300_6)/(950.-count6);
    xnt[index++]=nModsigz500_6/(950.-count6);
    xnt[index++]=TMath::Sqrt(nModsigz500_6)/(950.-count6);
    // 168 - fine SSD L6
    ntits->Fill(xnt);
    //		cout<<"\n\nirun sDD"<<iRun<<endl<<endl;

	} // loop on runs

    TFile* outfil=new TFile(ntupleFileName.Data(),"recreate");
    outfil->cd();
    ntits->Write();
    outfil->Close();
    delete outfil;

    //    MakePlots(ntupleFileName);
}



void MakePlots(Int_t run1,Int_t run2,TString ntupleFileName){
  if(run1>=run2){
    printf("******   ERROR: invalid run range %d - %d\n",run1,run2);
    return;
  }
  TFile* fil=new TFile(ntupleFileName.Data(),"read");
  if(!fil){
    printf("File with ntuple does not exist\n");
    return;
  }
  TNtuple* ntits=(TNtuple*)fil->Get("ntitsalign");

  Float_t nrun;

  Float_t SDDfracModmux20_3,SDDfracModmux50_3,SDDfracModmux100_3,SDDfracModsigx100_3,SDDfracModsigx200_3,SDDfracModsigx300_3;
  Float_t SDDfracModmuxlr20_3,SDDfracModmuxlr50_3,SDDfracModmuxlr100_3,SDDfracModexcx100_3,SDDfracModexcx200_3,SDDfracModexcx300_3;
  Float_t SDDfracModmuz50_3,SDDfracModmuz100_3,SDDfracModmuz300_3,SDDfracModsigz100_3,SDDfracModsigz300_3,SDDfracModsigz500_3;
  Float_t errSDDfracModmux20_3,errSDDfracModmux50_3,errSDDfracModmux100_3,errSDDfracModsigx100_3,errSDDfracModsigx200_3,errSDDfracModsigx300_3;
  Float_t errSDDfracModmuxlr20_3,errSDDfracModmuxlr50_3,errSDDfracModmuxlr100_3,errSDDfracModexcx100_3,errSDDfracModexcx200_3,errSDDfracModexcx300_3;
  Float_t errSDDfracModmuz50_3,errSDDfracModmuz100_3,errSDDfracModmuz300_3,errSDDfracModsigz100_3,errSDDfracModsigz300_3,errSDDfracModsigz500_3;

  Float_t SDDfracModmux20_4,SDDfracModmux50_4,SDDfracModmux100_4,SDDfracModsigx100_4,SDDfracModsigx200_4,SDDfracModsigx300_4;
  Float_t SDDfracModmuxlr20_4,SDDfracModmuxlr50_4,SDDfracModmuxlr100_4,SDDfracModexcx100_4,SDDfracModexcx200_4,SDDfracModexcx300_4;
  Float_t SDDfracModmuz50_4,SDDfracModmuz100_4,SDDfracModmuz300_4,SDDfracModsigz100_4,SDDfracModsigz300_4,SDDfracModsigz500_4;
  Float_t errSDDfracModmux20_4,errSDDfracModmux50_4,errSDDfracModmux100_4,errSDDfracModsigx100_4,errSDDfracModsigx200_4,errSDDfracModsigx300_4;
  Float_t errSDDfracModmuxlr20_4,errSDDfracModmuxlr50_4,errSDDfracModmuxlr100_4,errSDDfracModexcx100_4,errSDDfracModexcx200_4,errSDDfracModexcx300_4;
  Float_t errSDDfracModmuz50_4,errSDDfracModmuz100_4,errSDDfracModmuz300_4,errSDDfracModsigz100_4,errSDDfracModsigz300_4,errSDDfracModsigz500_4;

  Float_t SPDfracModmux20_1,SPDfracModmux50_1,SPDfracModmux100_1,SPDfracModsigx100_1,SPDfracModsigx200_1,SPDfracModsigx300_1;
  Float_t SPDfracModmuz50_1,SPDfracModmuz100_1,SPDfracModmuz300_1,SPDfracModsigz100_1,SPDfracModsigz300_1,SPDfracModsigz500_1;
  Float_t errSPDfracModmux20_1,errSPDfracModmux50_1,errSPDfracModmux100_1,errSPDfracModsigx100_1,errSPDfracModsigx200_1,errSPDfracModsigx300_1;
  Float_t errSPDfracModmuz50_1,errSPDfracModmuz100_1,errSPDfracModmuz300_1,errSPDfracModsigz100_1,errSPDfracModsigz300_1,errSPDfracModsigz500_1;
  Float_t SPDfracModmux20_2,SPDfracModmux50_2,SPDfracModmux100_2,SPDfracModsigx100_2,SPDfracModsigx200_2,SPDfracModsigx300_2;
  Float_t SPDfracModmuz50_2,SPDfracModmuz100_2,SPDfracModmuz300_2,SPDfracModsigz100_2,SPDfracModsigz300_2,SPDfracModsigz500_2;
  Float_t errSPDfracModmux20_2,errSPDfracModmux50_2,errSPDfracModmux100_2,errSPDfracModsigx100_2,errSPDfracModsigx200_2,errSPDfracModsigx300_2;
  Float_t errSPDfracModmuz50_2,errSPDfracModmuz100_2,errSPDfracModmuz300_2,errSPDfracModsigz100_2,errSPDfracModsigz300_2,errSPDfracModsigz500_2;

  Float_t SSDfracModmux20_5,SSDfracModmux50_5,SSDfracModmux100_5,SSDfracModsigx100_5,SSDfracModsigx200_5,SSDfracModsigx300_5;
  Float_t SSDfracModmuz50_5,SSDfracModmuz100_5,SSDfracModmuz300_5,SSDfracModsigz100_5,SSDfracModsigz300_5,SSDfracModsigz500_5;
  Float_t errSSDfracModmux20_5,errSSDfracModmux50_5,errSSDfracModmux100_5,errSSDfracModsigx100_5,errSSDfracModsigx200_5,errSSDfracModsigx300_5;
  Float_t errSSDfracModmuz50_5,errSSDfracModmuz100_5,errSSDfracModmuz300_5,errSSDfracModsigz100_5,errSSDfracModsigz300_5,errSSDfracModsigz500_5;
  Float_t SSDfracModmux20_6,SSDfracModmux50_6,SSDfracModmux100_6,SSDfracModsigx100_6,SSDfracModsigx200_6,SSDfracModsigx300_6;
  Float_t SSDfracModmuz50_6,SSDfracModmuz100_6,SSDfracModmuz300_6,SSDfracModsigz100_6,SSDfracModsigz300_6,SSDfracModsigz500_6;
  Float_t errSSDfracModmux20_6,errSSDfracModmux50_6,errSSDfracModmux100_6,errSSDfracModsigx100_6,errSSDfracModsigx200_6,errSSDfracModsigx300_6;
  Float_t errSSDfracModmuz50_6,errSSDfracModmuz100_6,errSSDfracModmuz300_6,errSSDfracModsigz100_6,errSSDfracModsigz300_6,errSSDfracModsigz500_6;

  ntits->SetBranchAddress("nrun",&nrun);                                         // SDD
  ntits->SetBranchAddress("SDDfracModmux20_3",&SDDfracModmux20_3);
  ntits->SetBranchAddress("errSDDfracModmux20_3",&errSDDfracModmux20_3);
  ntits->SetBranchAddress("SDDfracModmux50_3",&SDDfracModmux50_3);
  ntits->SetBranchAddress("errSDDfracModmux50_3",&errSDDfracModmux50_3);
  ntits->SetBranchAddress("SDDfracModmux100_3",&SDDfracModmux100_3);
  ntits->SetBranchAddress("errSDDfracModmux100_3",&errSDDfracModmux100_3);
  ntits->SetBranchAddress("SDDfracModsigx100_3",&SDDfracModsigx100_3);
  ntits->SetBranchAddress("errSDDfracModsigx100_3",&errSDDfracModsigx100_3);
  ntits->SetBranchAddress("SDDfracModsigx200_3",&SDDfracModsigx200_3);
  ntits->SetBranchAddress("errSDDfracModsigx200_3",&errSDDfracModsigx200_3);
  ntits->SetBranchAddress("SDDfracModsigx300_3",&SDDfracModsigx300_3);
  ntits->SetBranchAddress("errSDDfracModsigx300_3",&errSDDfracModsigx300_3);
  ntits->SetBranchAddress("SDDfracModmuxlr20_3",&SDDfracModmuxlr20_3);
  ntits->SetBranchAddress("errSDDfracModmuxlr20_3",&errSDDfracModmuxlr20_3);
  ntits->SetBranchAddress("SDDfracModmuxlr50_3",&SDDfracModmuxlr50_3);
  ntits->SetBranchAddress("errSDDfracModmuxlr50_3",&errSDDfracModmuxlr50_3);
  ntits->SetBranchAddress("SDDfracModmuxlr100_3",&SDDfracModmuxlr100_3);
  ntits->SetBranchAddress("errSDDfracModmuxlr100_3",&errSDDfracModmuxlr100_3);
  ntits->SetBranchAddress("SDDfracModexcx100_3",&SDDfracModexcx100_3);
  ntits->SetBranchAddress("errSDDfracModexcx100_3",&errSDDfracModexcx100_3);
  ntits->SetBranchAddress("SDDfracModexcx200_3",&SDDfracModexcx200_3);
  ntits->SetBranchAddress("errSDDfracModexcx200_3",&errSDDfracModexcx200_3);
  ntits->SetBranchAddress("SDDfracModexcx300_3",&SDDfracModexcx300_3);
  ntits->SetBranchAddress("errSDDfracModexcx300_3",&errSDDfracModexcx300_3);
  ntits->SetBranchAddress("SDDfracModmuz50_3",&SDDfracModmuz50_3);
  ntits->SetBranchAddress("errSDDfracModmuz50_3",&errSDDfracModmuz50_3);
  ntits->SetBranchAddress("SDDfracModmuz100_3",&SDDfracModmuz100_3);
  ntits->SetBranchAddress("errSDDfracModmuz100_3",&errSDDfracModmuz100_3);
  ntits->SetBranchAddress("SDDfracModmuz300_3",&SDDfracModmuz300_3);
  ntits->SetBranchAddress("errSDDfracModmuz300_3",&errSDDfracModmuz300_3);
  ntits->SetBranchAddress("SDDfracModsigz100_3",&SDDfracModsigz100_3);
  ntits->SetBranchAddress("errSDDfracModsigz100_3",&errSDDfracModsigz100_3);
  ntits->SetBranchAddress("SDDfracModsigz300_3",&SDDfracModsigz300_3);
  ntits->SetBranchAddress("errSDDfracModsigz300_3",&errSDDfracModsigz300_3);
  ntits->SetBranchAddress("SDDfracModsigz500_3",&SDDfracModsigz500_3);
  ntits->SetBranchAddress("errSDDfracModsigz500_3",&errSDDfracModsigz500_3);
  ntits->SetBranchAddress("SDDfracModmux20_4",&SDDfracModmux20_4);
  ntits->SetBranchAddress("errSDDfracModmux20_4",&errSDDfracModmux20_4);
  ntits->SetBranchAddress("SDDfracModmux50_4",&SDDfracModmux50_4);
  ntits->SetBranchAddress("errSDDfracModmux50_4",&errSDDfracModmux50_4);
  ntits->SetBranchAddress("SDDfracModmux100_4",&SDDfracModmux100_4);
  ntits->SetBranchAddress("errSDDfracModmux100_4",&errSDDfracModmux100_4);
  ntits->SetBranchAddress("SDDfracModsigx100_4",&SDDfracModsigx100_4);
  ntits->SetBranchAddress("errSDDfracModsigx100_4",&errSDDfracModsigx100_4);
  ntits->SetBranchAddress("SDDfracModsigx200_4",&SDDfracModsigx200_4);
  ntits->SetBranchAddress("errSDDfracModsigx200_4",&errSDDfracModsigx200_4);
  ntits->SetBranchAddress("SDDfracModsigx300_4",&SDDfracModsigx300_4);
  ntits->SetBranchAddress("errSDDfracModsigx300_4",&errSDDfracModsigx300_4);
  ntits->SetBranchAddress("SDDfracModmuxlr20_4",&SDDfracModmuxlr20_4);
  ntits->SetBranchAddress("errSDDfracModmuxlr20_4",&errSDDfracModmuxlr20_4);
  ntits->SetBranchAddress("SDDfracModmuxlr50_4",&SDDfracModmuxlr50_4);
  ntits->SetBranchAddress("errSDDfracModmuxlr50_4",&errSDDfracModmuxlr50_4);
  ntits->SetBranchAddress("SDDfracModmuxlr100_4",&SDDfracModmuxlr100_4);
  ntits->SetBranchAddress("errSDDfracModmuxlr100_4",&errSDDfracModmuxlr100_4);
  ntits->SetBranchAddress("SDDfracModexcx100_4",&SDDfracModexcx100_4);
  ntits->SetBranchAddress("errSDDfracModexcx100_4",&errSDDfracModexcx100_4);
  ntits->SetBranchAddress("SDDfracModexcx200_4",&SDDfracModexcx200_4);
  ntits->SetBranchAddress("errSDDfracModexcx200_4",&errSDDfracModexcx200_4);
  ntits->SetBranchAddress("SDDfracModexcx300_4",&SDDfracModexcx300_4);
  ntits->SetBranchAddress("errSDDfracModexcx300_4",&errSDDfracModexcx300_4);
  ntits->SetBranchAddress("SDDfracModmuz50_4",&SDDfracModmuz50_4);
  ntits->SetBranchAddress("errSDDfracModmuz50_4",&errSDDfracModmuz50_4);
  ntits->SetBranchAddress("SDDfracModmuz100_4",&SDDfracModmuz100_4);
  ntits->SetBranchAddress("errSDDfracModmuz100_4",&errSDDfracModmuz100_4);
  ntits->SetBranchAddress("SDDfracModmuz300_4",&SDDfracModmuz300_4);
  ntits->SetBranchAddress("errSDDfracModmuz300_4",&errSDDfracModmuz300_4);
  ntits->SetBranchAddress("SDDfracModsigz100_4",&SDDfracModsigz100_4);
  ntits->SetBranchAddress("errSDDfracModsigz100_4",&errSDDfracModsigz100_4);
  ntits->SetBranchAddress("SDDfracModsigz300_4",&SDDfracModsigz300_4);
  ntits->SetBranchAddress("errSDDfracModsigz300_4",&errSDDfracModsigz300_4);
  ntits->SetBranchAddress("SDDfracModsigz500_4",&SDDfracModsigz500_4);
  ntits->SetBranchAddress("errSDDfracModsigz500_4",&errSDDfracModsigz500_4);

  ntits->SetBranchAddress("SPDfracModmux20_1",&SPDfracModmux20_1);                                  // SPD
  ntits->SetBranchAddress("errSPDfracModmux20_1",&errSPDfracModmux20_1);
  ntits->SetBranchAddress("SPDfracModmux50_1",&SPDfracModmux50_1);
  ntits->SetBranchAddress("errSPDfracModmux50_1",&errSPDfracModmux50_1);
  ntits->SetBranchAddress("SPDfracModmux100_1",&SPDfracModmux100_1);
  ntits->SetBranchAddress("errSPDfracModmux100_1",&errSPDfracModmux100_1);
  ntits->SetBranchAddress("SPDfracModsigx100_1",&SPDfracModsigx100_1);
  ntits->SetBranchAddress("errSPDfracModsigx100_1",&errSPDfracModsigx100_1);
  ntits->SetBranchAddress("SPDfracModsigx200_1",&SPDfracModsigx200_1);
  ntits->SetBranchAddress("errSPDfracModsigx200_1",&errSPDfracModsigx200_1);
  ntits->SetBranchAddress("SPDfracModsigx300_1",&SPDfracModsigx300_1);
  ntits->SetBranchAddress("errSPDfracModsigx300_1",&errSPDfracModsigx300_1);
  ntits->SetBranchAddress("SPDfracModmuz50_1",&SPDfracModmuz50_1);
  ntits->SetBranchAddress("errSPDfracModmuz50_1",&errSPDfracModmuz50_1);
  ntits->SetBranchAddress("SPDfracModmuz100_1",&SPDfracModmuz100_1);
  ntits->SetBranchAddress("errSPDfracModmuz100_1",&errSPDfracModmuz100_1);
  ntits->SetBranchAddress("SPDfracModmuz300_1",&SPDfracModmuz300_1);
  ntits->SetBranchAddress("errSPDfracModmuz300_1",&errSPDfracModmuz300_1);
  ntits->SetBranchAddress("SPDfracModsigz100_1",&SPDfracModsigz100_1);
  ntits->SetBranchAddress("errSPDfracModsigz100_1",&errSPDfracModsigz100_1);
  ntits->SetBranchAddress("SPDfracModsigz300_1",&SPDfracModsigz300_1);
  ntits->SetBranchAddress("errSPDfracModsigz300_1",&errSPDfracModsigz300_1);
  ntits->SetBranchAddress("SPDfracModsigz500_1",&SPDfracModsigz500_1);
  ntits->SetBranchAddress("errSPDfracModsigz500_1",&errSPDfracModsigz500_1);
  ntits->SetBranchAddress("SPDfracModmux20_2",&SPDfracModmux20_2);
  ntits->SetBranchAddress("errSPDfracModmux20_2",&errSPDfracModmux20_2);
  ntits->SetBranchAddress("SPDfracModmux50_2",&SPDfracModmux50_2);
  ntits->SetBranchAddress("errSPDfracModmux50_2",&errSPDfracModmux50_2);
  ntits->SetBranchAddress("SPDfracModmux100_2",&SPDfracModmux100_2);
  ntits->SetBranchAddress("errSPDfracModmux100_2",&errSPDfracModmux100_2);
  ntits->SetBranchAddress("SPDfracModsigx100_2",&SPDfracModsigx100_2);
  ntits->SetBranchAddress("errSPDfracModsigx100_2",&errSPDfracModsigx100_2);
  ntits->SetBranchAddress("SPDfracModsigx200_2",&SPDfracModsigx200_2);
  ntits->SetBranchAddress("errSPDfracModsigx200_2",&errSPDfracModsigx200_2);
  ntits->SetBranchAddress("SPDfracModsigx300_2",&SPDfracModsigx300_2);
  ntits->SetBranchAddress("errSPDfracModsigx300_2",&errSPDfracModsigx300_2);
  ntits->SetBranchAddress("SPDfracModmuz50_2",&SPDfracModmuz50_2);
  ntits->SetBranchAddress("errSPDfracModmuz50_2",&errSPDfracModmuz50_2);
  ntits->SetBranchAddress("SPDfracModmuz100_2",&SPDfracModmuz100_2);
  ntits->SetBranchAddress("errSPDfracModmuz100_2",&errSPDfracModmuz100_2);
  ntits->SetBranchAddress("SPDfracModmuz300_2",&SPDfracModmuz300_2);
  ntits->SetBranchAddress("errSPDfracModmuz300_2",&errSPDfracModmuz300_2);
  ntits->SetBranchAddress("SPDfracModsigz100_2",&SPDfracModsigz100_2);
  ntits->SetBranchAddress("errSPDfracModsigz100_2",&errSPDfracModsigz100_2);
  ntits->SetBranchAddress("SPDfracModsigz300_2",&SPDfracModsigz300_2);
  ntits->SetBranchAddress("errSPDfracModsigz300_2",&errSPDfracModsigz300_2);
  ntits->SetBranchAddress("SPDfracModsigz500_2",&SPDfracModsigz500_2);
  ntits->SetBranchAddress("errSPDfracModsigz500_2",&errSPDfracModsigz500_2);

  ntits->SetBranchAddress("SSDfracModmux20_5",&SSDfracModmux20_5);                                   // SSD
  ntits->SetBranchAddress("errSSDfracModmux20_5",&errSSDfracModmux20_5);
  ntits->SetBranchAddress("SSDfracModmux50_5",&SSDfracModmux50_5);
  ntits->SetBranchAddress("errSSDfracModmux50_5",&errSSDfracModmux50_5);
  ntits->SetBranchAddress("SSDfracModmux100_5",&SSDfracModmux100_5);
  ntits->SetBranchAddress("errSSDfracModmux100_5",&errSSDfracModmux100_5);
  ntits->SetBranchAddress("SSDfracModsigx100_5",&SSDfracModsigx100_5);
  ntits->SetBranchAddress("errSSDfracModsigx100_5",&errSSDfracModsigx100_5);
  ntits->SetBranchAddress("SSDfracModsigx200_5",&SSDfracModsigx200_5);
  ntits->SetBranchAddress("errSSDfracModsigx200_5",&errSSDfracModsigx200_5);
  ntits->SetBranchAddress("SSDfracModsigx300_5",&SSDfracModsigx300_5);
  ntits->SetBranchAddress("errSSDfracModsigx300_5",&errSSDfracModsigx300_5);
  ntits->SetBranchAddress("SSDfracModmuz50_5",&SSDfracModmuz50_5);
  ntits->SetBranchAddress("errSSDfracModmuz50_5",&errSSDfracModmuz50_5);
  ntits->SetBranchAddress("SSDfracModmuz100_5",&SSDfracModmuz100_5);
  ntits->SetBranchAddress("errSSDfracModmuz100_5",&errSSDfracModmuz100_5);
  ntits->SetBranchAddress("SSDfracModmuz300_5",&SSDfracModmuz300_5);
  ntits->SetBranchAddress("errSSDfracModmuz300_5",&errSSDfracModmuz300_5);
  ntits->SetBranchAddress("SSDfracModsigz100_5",&SSDfracModsigz100_5);
  ntits->SetBranchAddress("errSSDfracModsigz100_5",&errSSDfracModsigz100_5);
  ntits->SetBranchAddress("SSDfracModsigz300_5",&SSDfracModsigz300_5);
  ntits->SetBranchAddress("errSSDfracModsigz300_5",&errSSDfracModsigz300_5);
  ntits->SetBranchAddress("SSDfracModsigz500_5",&SSDfracModsigz500_5);
  ntits->SetBranchAddress("errSSDfracModsigz500_5",&errSSDfracModsigz500_5);
  ntits->SetBranchAddress("SSDfracModmux20_6",&SSDfracModmux20_6);
  ntits->SetBranchAddress("errSSDfracModmux20_6",&errSSDfracModmux20_6);
  ntits->SetBranchAddress("SSDfracModmux50_6",&SSDfracModmux50_6);
  ntits->SetBranchAddress("errSSDfracModmux50_6",&errSSDfracModmux50_6);
  ntits->SetBranchAddress("SSDfracModmux100_6",&SSDfracModmux100_6);
  ntits->SetBranchAddress("errSSDfracModmux100_6",&errSSDfracModmux100_6);
  ntits->SetBranchAddress("SSDfracModsigx100_6",&SSDfracModsigx100_6);
  ntits->SetBranchAddress("errSSDfracModsigx100_6",&errSSDfracModsigx100_6);
  ntits->SetBranchAddress("SSDfracModsigx200_6",&SSDfracModsigx200_6);
  ntits->SetBranchAddress("errSSDfracModsigx200_6",&errSSDfracModsigx200_6);
  ntits->SetBranchAddress("SSDfracModsigx300_6",&SSDfracModsigx300_6);
  ntits->SetBranchAddress("errSSDfracModsigx300_6",&errSSDfracModsigx300_6);
  ntits->SetBranchAddress("SSDfracModmuz50_6",&SSDfracModmuz50_6);
  ntits->SetBranchAddress("errSSDfracModmuz50_6",&errSSDfracModmuz50_6);
  ntits->SetBranchAddress("SSDfracModmuz100_6",&SSDfracModmuz100_6);
  ntits->SetBranchAddress("errSSDfracModmuz100_6",&errSSDfracModmuz100_6);
  ntits->SetBranchAddress("SSDfracModmuz300_6",&SSDfracModmuz300_6);
  ntits->SetBranchAddress("errSSDfracModmuz300_6",&errSSDfracModmuz300_6);
  ntits->SetBranchAddress("SSDfracModsigz100_6",&SSDfracModsigz100_6);
  ntits->SetBranchAddress("errSSDfracModsigz100_6",&errSSDfracModsigz100_6);
  ntits->SetBranchAddress("SSDfracModsigz300_6",&SSDfracModsigz300_6);
  ntits->SetBranchAddress("errSSDfracModsigz300_6",&errSSDfracModsigz300_6);
  ntits->SetBranchAddress("SSDfracModsigz500_6",&SSDfracModsigz500_6);
  ntits->SetBranchAddress("errSSDfracModsigz500_6",&errSSDfracModsigz500_6);

  // Sort entries according to run number in the chosen range
  // Same order is assumed for all the subsequent ntuples
  Int_t nr=ntits->GetEntries();
  Int_t *myIndex = new Int_t [nr];
  Int_t *noRuns = new Int_t [nr];
  Int_t kRunsToPlot=0;
  printf("Processing runs from %d up to %d\n",run1,run2);
  for(Int_t i=0; i<nr;i++){
    ntits->GetEvent(i);
    Int_t intrun = static_cast<Int_t>(nrun+0.01);
    if(intrun>=run1 && intrun<=run2){
      printf("Accepting run number %d in position %d\n",intrun,kRunsToPlot);
      noRuns[i]=intrun;
      kRunsToPlot++;
    }
    else {
      noRuns[i]=run2+10;  
      printf("Rejecting run number %d - out of range\n",intrun);
    }
  }
  TMath::Sort(nr,noRuns,myIndex,kFALSE);
  printf("Total number of runs accepted for display %d\n",kRunsToPlot);
  if(kRunsToPlot==0)return;
  for(Int_t i=0;i<kRunsToPlot;i++)printf("Position %d ) Run: %d\n",i,noRuns[myIndex[i]]);

  TH1F* hfmux20_3=new TH1F("hfmux20_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmux50_3=new TH1F("hfmux50_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmux100_3=new TH1F("hfmux100_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx100_3=new TH1F("hfsigx100_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx200_3=new TH1F("hfsigx200_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx300_3=new TH1F("hfsigx300_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuxlr20_3=new TH1F("hfmuxlr20_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuxlr50_3=new TH1F("hfmuxlr50_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuxlr100_3=new TH1F("hfmuxlr100_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfexcx100_3=new TH1F("hfexcx100_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfexcx200_3=new TH1F("hfexcx200_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfexcx300_3=new TH1F("hfexcx300_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz50_3=new TH1F("hfmuz50_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz100_3=new TH1F("hfmuz100_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz300_3=new TH1F("hfmuz300_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz100_3=new TH1F("hfsigz100_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz300_3=new TH1F("hfsigz300_3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz500_3=new TH1F("hfsigz500_3","",kRunsToPlot,0.,kRunsToPlot);

  TH1F* hfmux20_4=new TH1F("hfmux20_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmux50_4=new TH1F("hfmux50_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmux100_4=new TH1F("hfmux100_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx100_4=new TH1F("hfsigx100_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx200_4=new TH1F("hfsigx200_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx300_4=new TH1F("hfsigx300_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuxlr20_4=new TH1F("hfmuxlr20_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuxlr50_4=new TH1F("hfmuxlr50_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuxlr100_4=new TH1F("hfmuxlr100_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfexcx100_4=new TH1F("hfexcx100_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfexcx200_4=new TH1F("hfexcx200_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfexcx300_4=new TH1F("hfexcx300_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz50_4=new TH1F("hfmuz50_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz100_4=new TH1F("hfmuz100_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz300_4=new TH1F("hfmuz300_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz100_4=new TH1F("hfsigz100_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz300_4=new TH1F("hfsigz300_4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz500_4=new TH1F("hfsigz500_4","",kRunsToPlot,0.,kRunsToPlot);

  TH1F* hfmux20_1=new TH1F("hfmux20_1","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmux50_1=new TH1F("hfmux50_1","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmux100_1=new TH1F("hfmux100_1","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx100_1=new TH1F("hfsigx100_1","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx200_1=new TH1F("hfsigx200_1","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx300_1=new TH1F("hfsigx300_1","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz50_1=new TH1F("hfmuz50_1","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz100_1=new TH1F("hfmuz100_1","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz300_1=new TH1F("hfmuz300_1","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz100_1=new TH1F("hfsigz100_1","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz300_1=new TH1F("hfsigz300_1","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz500_1=new TH1F("hfsigz500_1","",kRunsToPlot,0.,kRunsToPlot);

  TH1F* hfmux20_2=new TH1F("hfmux20_2","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmux50_2=new TH1F("hfmux50_2","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmux100_2=new TH1F("hfmux100_2","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx100_2=new TH1F("hfsigx100_2","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx200_2=new TH1F("hfsigx200_2","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx300_2=new TH1F("hfsigx300_2","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz50_2=new TH1F("hfmuz50_2","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz100_2=new TH1F("hfmuz100_2","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz300_2=new TH1F("hfmuz300_2","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz100_2=new TH1F("hfsigz100_2","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz300_2=new TH1F("hfsigz300_2","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz500_2=new TH1F("hfsigz500_2","",kRunsToPlot,0.,kRunsToPlot);

  TH1F* hfmux20_5=new TH1F("hfmux20_5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmux50_5=new TH1F("hfmux50_5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmux100_5=new TH1F("hfmux100_5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx100_5=new TH1F("hfsigx100_5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx200_5=new TH1F("hfsigx200_5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx300_5=new TH1F("hfsigx300_5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz50_5=new TH1F("hfmuz50_5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz100_5=new TH1F("hfmuz100_5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz300_5=new TH1F("hfmuz300_5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz100_5=new TH1F("hfsigz100_5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz300_5=new TH1F("hfsigz300_5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz500_5=new TH1F("hfsigz500_5","",kRunsToPlot,0.,kRunsToPlot);

  TH1F* hfmux20_6=new TH1F("hfmux20_6","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmux50_6=new TH1F("hfmux50_6","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmux100_6=new TH1F("hfmux100_6","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx100_6=new TH1F("hfsigx100_6","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx200_6=new TH1F("hfsigx200_6","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigx300_6=new TH1F("hfsigx300_6","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz50_6=new TH1F("hfmuz50_6","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz100_6=new TH1F("hfmuz100_6","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfmuz300_6=new TH1F("hfmuz300_6","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz100_6=new TH1F("hfsigz100_6","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz300_6=new TH1F("hfsigz300_6","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* hfsigz500_6=new TH1F("hfsigz500_6","",kRunsToPlot,0.,kRunsToPlot);


  for(Int_t i=0; i<kRunsToPlot;i++){
    ntits->GetEvent(myIndex[i]);

    hfmux20_3->SetBinContent(i+1,SDDfracModmux20_3);
//    hfmux20_3->SetBinError(i+1,errSDDfracModmux20_3);
    hfmux20_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmux50_3->SetBinContent(i+1,SDDfracModmux50_3);
//    hfmux50_3->SetBinError(i+1,errSDDfracModmux50_3);
    hfmux50_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmux100_3->SetBinContent(i+1,SDDfracModmux100_3);
//    hfmux100_3->SetBinError(i+1,errSDDfracModmux100_3);
    hfmux100_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx100_3->SetBinContent(i+1,SDDfracModsigx100_3);
//    hfsigx100_3->SetBinError(i+1,errSDDfracModsigx100_3);
    hfsigx100_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx200_3->SetBinContent(i+1,SDDfracModsigx200_3);
//    hfsigx200_3->SetBinError(i+1,errSDDfracModsigx200_3);
    hfsigx200_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx300_3->SetBinContent(i+1,SDDfracModsigx300_3);
//    hfsigx300_3->SetBinError(i+1,errSDDfracModsigx300_3);
    hfsigx300_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuxlr20_3->SetBinContent(i+1,SDDfracModmuxlr20_3);
//    hfmuxlr20_3->SetBinError(i+1,errSDDfracModmuxlr20_3);
    hfmuxlr20_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuxlr50_3->SetBinContent(i+1,SDDfracModmuxlr50_3);
//    hfmuxlr50_3->SetBinError(i+1,errSDDfracModmux50_3);
    hfmuxlr50_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuxlr100_3->SetBinContent(i+1,SDDfracModmux100_3);
//    hfmuxlr100_3->SetBinError(i+1,errSDDfracModmux100_3);
    hfmuxlr100_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfexcx100_3->SetBinContent(i+1,SDDfracModexcx100_3);
//    hfexcx100_3->SetBinError(i+1,errSDDfracModexcx100_3);
    hfexcx100_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfexcx200_3->SetBinContent(i+1,SDDfracModexcx200_3);
//    hfexcx200_3->SetBinError(i+1,errSDDfracModexcx200_3);
    hfexcx200_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfexcx300_3->SetBinContent(i+1,SDDfracModexcx300_3);
//    hfexcx300_3->SetBinError(i+1,errSDDfracModexcx300_3);
    hfexcx300_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz50_3->SetBinContent(i+1,SDDfracModmuz50_3);
//    hfmuz50_3->SetBinError(i+1,errSDDfracModmuz50_3);
    hfmuz50_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz100_3->SetBinContent(i+1,SDDfracModmuz100_3);
//    hfmuz100_3->SetBinError(i+1,errSDDfracModmuz100_3);
    hfmuz100_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz300_3->SetBinContent(i+1,SDDfracModmuz300_3);
//    hfmuz300_3->SetBinError(i+1,errSDDfracModmuz300_3);
    hfmuz300_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz100_3->SetBinContent(i+1,SDDfracModsigz100_3);
//    hfsigz100_3->SetBinError(i+1,errSDDfracModsigz100_3);
    hfsigz100_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz300_3->SetBinContent(i+1,SDDfracModsigz300_3);
//    hfsigz300_3->SetBinError(i+1,errSDDfracModsigz300_3);
    hfsigz300_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz500_3->SetBinContent(i+1,SDDfracModsigz500_3);
//    hfsigz500_3->SetBinError(i+1,errSDDfracModsigz500_3);
    hfsigz500_3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

    hfmux20_4->SetBinContent(i+1,SDDfracModmux20_4);
//    hfmux20_4->SetBinError(i+1,errSDDfracModmux20_4);
    hfmux20_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmux50_4->SetBinContent(i+1,SDDfracModmux50_4);
//    hfmux50_4->SetBinError(i+1,errSDDfracModmux50_4);
    hfmux50_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmux100_4->SetBinContent(i+1,SDDfracModmux100_4);
//    hfmux100_4->SetBinError(i+1,errSDDfracModmux100_4);
    hfmux100_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx100_4->SetBinContent(i+1,SDDfracModsigx100_4);
//    hfsigx100_4->SetBinError(i+1,errSDDfracModsigx100_4);
    hfsigx100_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx200_4->SetBinContent(i+1,SDDfracModsigx200_4);
//    hfsigx200_4->SetBinError(i+1,errSDDfracModsigx200_4);
    hfsigx200_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx300_4->SetBinContent(i+1,SDDfracModsigx300_4);
//    hfsigx300_4->SetBinError(i+1,errSDDfracModsigx300_4);
    hfsigx300_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuxlr20_4->SetBinContent(i+1,SDDfracModmuxlr20_4);
//    hfmuxlr20_4->SetBinError(i+1,errSDDfracModmuxlr20_4);
    hfmuxlr20_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuxlr50_4->SetBinContent(i+1,SDDfracModmuxlr50_4);
//    hfmuxlr50_4->SetBinError(i+1,errSDDfracModmux50_4);
    hfmuxlr50_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuxlr100_4->SetBinContent(i+1,SDDfracModmux100_4);
//    hfmuxlr100_4->SetBinError(i+1,errSDDfracModmux100_4);
    hfmuxlr100_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfexcx100_4->SetBinContent(i+1,SDDfracModexcx100_4);
//    hfexcx100_4->SetBinError(i+1,errSDDfracModexcx100_4);
    hfexcx100_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfexcx200_4->SetBinContent(i+1,SDDfracModexcx200_4);
//    hfexcx200_4->SetBinError(i+1,errSDDfracModexcx200_4);
    hfexcx200_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfexcx300_4->SetBinContent(i+1,SDDfracModexcx300_4);
//    hfexcx300_4->SetBinError(i+1,errSDDfracModexcx300_4);
    hfexcx300_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz50_4->SetBinContent(i+1,SDDfracModmuz50_4);
//    hfmuz50_4->SetBinError(i+1,errSDDfracModmuz50_4);
    hfmuz50_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz100_4->SetBinContent(i+1,SDDfracModmuz100_4);
//    hfmuz100_4->SetBinError(i+1,errSDDfracModmuz100_4);
    hfmuz100_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz300_4->SetBinContent(i+1,SDDfracModmuz300_4);
//    hfmuz300_4->SetBinError(i+1,errSDDfracModmuz300_4);
    hfmuz300_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz100_4->SetBinContent(i+1,SDDfracModsigz100_4);
//    hfsigz100_4->SetBinError(i+1,errSDDfracModsigz100_4);
    hfsigz100_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz300_4->SetBinContent(i+1,SDDfracModsigz300_4);
//    hfsigz300_4->SetBinError(i+1,errSDDfracModsigz300_4);
    hfsigz300_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz500_4->SetBinContent(i+1,SDDfracModsigz500_4);
//    hfsigz500_4->SetBinError(i+1,errSDDfracModsigz500_4);
    hfsigz500_4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

    hfmux20_1->SetBinContent(i+1,SPDfracModmux20_1);
//    hfmux20_1->SetBinError(i+1,errSPDfracModmux20_1);
    hfmux20_1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmux50_1->SetBinContent(i+1,SPDfracModmux50_1);
//    hfmux50_1->SetBinError(i+1,errSPDfracModmux50_1);
    hfmux50_1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmux100_1->SetBinContent(i+1,SPDfracModmux100_1);
//    hfmux100_1->SetBinError(i+1,errSPDfracModmux100_1);
    hfmux100_1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx100_1->SetBinContent(i+1,SPDfracModsigx100_1);
//    hfsigx100_1->SetBinError(i+1,errSPDfracModsigx100_1);
    hfsigx100_1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx200_1->SetBinContent(i+1,SPDfracModsigx200_1);
//    hfsigx200_1->SetBinError(i+1,errSPDfracModsigx200_1);
    hfsigx200_1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx300_1->SetBinContent(i+1,SPDfracModsigx300_1);
//    hfsigx300_1->SetBinError(i+1,errSPDfracModsigx300_1);
    hfsigx300_1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz50_1->SetBinContent(i+1,SPDfracModmuz50_1);
//    hfmuz50_1->SetBinError(i+1,errSPDfracModmuz50_1);
    hfmuz50_1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz100_1->SetBinContent(i+1,SPDfracModmuz100_1);
//    hfmuz100_1->SetBinError(i+1,errSPDfracModmuz100_1);
    hfmuz100_1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz300_1->SetBinContent(i+1,SPDfracModmuz300_1);
//    hfmuz300_1->SetBinError(i+1,errSPDfracModmuz300_1);
    hfmuz300_1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz100_1->SetBinContent(i+1,SPDfracModsigz100_1);
//    hfsigz100_1->SetBinError(i+1,errSPDfracModsigz100_1);
    hfsigz100_1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz300_1->SetBinContent(i+1,SPDfracModsigz300_1);
//    hfsigz300_1->SetBinError(i+1,errSPDfracModsigz300_1);
    hfsigz300_1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz500_1->SetBinContent(i+1,SPDfracModsigz500_1);
//    hfsigz500_1->SetBinError(i+1,errSPDfracModsigz500_1);
    hfsigz500_1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

    hfmux20_2->SetBinContent(i+1,SPDfracModmux20_2);
//    hfmux20_2->SetBinError(i+1,errSPDfracModmux20_2);
    hfmux20_2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmux50_2->SetBinContent(i+1,SPDfracModmux50_2);
//    hfmux50_2->SetBinError(i+1,errSPDfracModmux50_2);
    hfmux50_2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmux100_2->SetBinContent(i+1,SPDfracModmux100_2);
//    hfmux100_2->SetBinError(i+1,errSPDfracModmux100_2);
    hfmux100_2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx100_2->SetBinContent(i+1,SPDfracModsigx100_2);
//    hfsigx100_2->SetBinError(i+1,errSPDfracModsigx100_2);
    hfsigx100_2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx200_2->SetBinContent(i+1,SPDfracModsigx200_2);
//    hfsigx200_2->SetBinError(i+1,errSPDfracModsigx200_2);
    hfsigx200_2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx300_2->SetBinContent(i+1,SPDfracModsigx300_2);
//    hfsigx300_2->SetBinError(i+1,errSPDfracModsigx300_2);
    hfsigx300_2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz50_2->SetBinContent(i+1,SPDfracModmuz50_2);
//    hfmuz50_2->SetBinError(i+1,errSPDfracModmuz50_2);
    hfmuz50_2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz100_2->SetBinContent(i+1,SPDfracModmuz100_2);
//    hfmuz100_2->SetBinError(i+1,errSPDfracModmuz100_2);
    hfmuz100_2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz300_2->SetBinContent(i+1,SPDfracModmuz300_2);
//    hfmuz300_2->SetBinError(i+1,errSPDfracModmuz300_2);
    hfmuz300_2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz100_2->SetBinContent(i+1,SPDfracModsigz100_2);
//    hfsigz100_2->SetBinError(i+1,errSPDfracModsigz100_2);
    hfsigz100_2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz300_2->SetBinContent(i+1,SPDfracModsigz300_2);
//    hfsigz300_2->SetBinError(i+1,errSPDfracModsigz300_2);
    hfsigz300_2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz500_2->SetBinContent(i+1,SPDfracModsigz500_2);
//    hfsigz500_2->SetBinError(i+1,errSPDfracModsigz500_2);
    hfsigz500_2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

    hfmux20_5->SetBinContent(i+1,SSDfracModmux20_5);
//    hfmux20_5->SetBinError(i+1,errSSDfracModmux20_5);
    hfmux20_5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmux50_5->SetBinContent(i+1,SSDfracModmux50_5);
//    hfmux50_5->SetBinError(i+1,errSSDfracModmux50_5);
    hfmux50_5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmux100_5->SetBinContent(i+1,SSDfracModmux100_5);
//    hfmux100_5->SetBinError(i+1,errSSDfracModmux100_5);
    hfmux100_5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx100_5->SetBinContent(i+1,SSDfracModsigx100_5);
//    hfsigx100_5->SetBinError(i+1,errSSDfracModsigx100_5);
    hfsigx100_5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx200_5->SetBinContent(i+1,SSDfracModsigx200_5);
//    hfsigx200_5->SetBinError(i+1,errSSDfracModsigx200_5);
    hfsigx200_5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx300_5->SetBinContent(i+1,SSDfracModsigx300_5);
//    hfsigx300_5->SetBinError(i+1,errSSDfracModsigx300_5);
    hfsigx300_5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz50_5->SetBinContent(i+1,SSDfracModmuz50_5);
//    hfmuz50_5->SetBinError(i+1,errSSDfracModmuz50_5);
    hfmuz50_5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz100_5->SetBinContent(i+1,SSDfracModmuz100_5);
//    hfmuz100_5->SetBinError(i+1,errSSDfracModmuz100_5);
    hfmuz100_5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz300_5->SetBinContent(i+1,SSDfracModmuz300_5);
//    hfmuz300_5->SetBinError(i+1,errSSDfracModmuz300_5);
    hfmuz300_5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz100_5->SetBinContent(i+1,SSDfracModsigz100_5);
//    hfsigz100_5->SetBinError(i+1,errSSDfracModsigz100_5);
    hfsigz100_5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz300_5->SetBinContent(i+1,SSDfracModsigz300_5);
//    hfsigz300_5->SetBinError(i+1,errSSDfracModsigz300_5);
    hfsigz300_5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz500_5->SetBinContent(i+1,SSDfracModsigz500_5);
//    hfsigz500_5->SetBinError(i+1,errSSDfracModsigz500_5);
    hfsigz500_5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

    hfmux20_6->SetBinContent(i+1,SSDfracModmux20_6);
//    hfmux20_6->SetBinError(i+1,errSSDfracModmux20_6);
    hfmux20_6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmux50_6->SetBinContent(i+1,SSDfracModmux50_6);
//    hfmux50_6->SetBinError(i+1,errSSDfracModmux50_6);
    hfmux50_6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmux100_6->SetBinContent(i+1,SSDfracModmux100_6);
//    hfmux100_6->SetBinError(i+1,errSSDfracModmux100_6);
    hfmux100_6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx100_6->SetBinContent(i+1,SSDfracModsigx100_6);
//    hfsigx100_6->SetBinError(i+1,errSSDfracModsigx100_6);
    hfsigx100_6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx200_6->SetBinContent(i+1,SSDfracModsigx200_6);
//    hfsigx200_6->SetBinError(i+1,errSSDfracModsigx200_6);
    hfsigx200_6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigx300_6->SetBinContent(i+1,SSDfracModsigx300_6);
//    hfsigx300_6->SetBinError(i+1,errSSDfracModsigx300_6);
    hfsigx300_6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz50_6->SetBinContent(i+1,SSDfracModmuz50_6);
//    hfmuz50_6->SetBinError(i+1,errSSDfracModmuz50_6);
    hfmuz50_6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz100_6->SetBinContent(i+1,SSDfracModmuz100_6);
//    hfmuz100_6->SetBinError(i+1,errSSDfracModmuz100_6);
    hfmuz100_6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfmuz300_6->SetBinContent(i+1,SSDfracModmuz300_6);
//    hfmuz300_6->SetBinError(i+1,errSSDfracModmuz300_6);
    hfmuz300_6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz100_6->SetBinContent(i+1,SSDfracModsigz100_6);
//    hfsigz100_6->SetBinError(i+1,errSSDfracModsigz100_6);
    hfsigz100_6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz300_6->SetBinContent(i+1,SSDfracModsigz300_6);
//    hfsigz300_6->SetBinError(i+1,errSSDfracModsigz300_6);
    hfsigz300_6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hfsigz500_6->SetBinContent(i+1,SSDfracModsigz500_6);
//    hfsigz500_6->SetBinError(i+1,errSSDfracModsigz500_6);
    hfsigz500_6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));


  }

  gROOT->SetStyle("Plain");
  //  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetFillColor(0);
  gStyle->SetTextFont(32);

  

  TCanvas* c1=new TCanvas("c1","SDD Layer 3", 1200,800);
  c1->Divide(1,2);
  c1->cd(1);
  hfmux20_3->SetTitle("fraction(| #mu x-residual | < threshold) - layer 3");
  hfmux20_3->SetLineColor(2);
  hfmux20_3->SetMarkerStyle(31);
  hfmux20_3->SetMarkerSize(2);
  hfmux20_3->SetMarkerColor(2);
  hfmux20_3->SetMinimum(0.);
  hfmux20_3->SetMaximum(1.2);
  hfmux20_3->GetXaxis()->SetTitle("run number");
  hfmux20_3->GetYaxis()->SetTitleOffset(1.2);
  //  hfmux20_3->GetYaxis()->SetTitle("%");
  hfmux20_3->Draw("P");
  hfmux50_3->SetLineColor(4);
  hfmux50_3->SetMarkerStyle(21);
  hfmux50_3->SetMarkerColor(4);
  hfmux50_3->Draw("same P");
  hfmux100_3->SetLineColor(6);
  hfmux100_3->SetMarkerStyle(22);
  hfmux100_3->SetMarkerColor(6);
  hfmux100_3->Draw("same P");

	// quiqui
	TLine *l3 = new TLine(0.,0.877,kRunsToPlot,0.877);
	l3->SetLineStyle(2);
	l3->SetLineWidth(2);	
//	l3->Draw();
  TLegend* leg=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
  TLegendEntry* enta=leg->AddEntry(hfmux20_3,"|#mu res x | < 20 #mum","PL"); // P: marker, L; line
  enta->SetTextColor(hfmux20_3->GetMarkerColor());
  TLegendEntry* entb=leg->AddEntry(hfmux50_3,"|#mu res x | < 50 #mum","PL"); // P: marker, L; line
  entb->SetTextColor(hfmux50_3->GetMarkerColor());
  TLegendEntry* entc=leg->AddEntry(hfmux100_3,"|#mu res x | < 100 #mum","PL"); // P: marker, L; line
  entc->SetTextColor(hfmux100_3->GetMarkerColor());
  leg->SetFillStyle(0);
  leg->Draw();
  //  c1->Update();

  c1->cd(2);
  hfmuxlr20_3->SetTitle("fraction(| #muL - #muR | < threshold) - layer 3");
  hfmuxlr20_3->SetLineColor(2);
  hfmuxlr20_3->SetMarkerStyle(31);
  hfmuxlr20_3->SetMarkerSize(2);
  hfmuxlr20_3->SetMarkerColor(2);
  hfmuxlr20_3->SetMinimum(0.);
  hfmuxlr20_3->SetMaximum(1.2);
  hfmuxlr20_3->GetXaxis()->SetTitle("run number");
  hfmuxlr20_3->GetYaxis()->SetTitleOffset(1.2);
  //  hfmuxlr20_3->GetYaxis()->SetTitle("%");
  hfmuxlr20_3->Draw("P");
  hfmuxlr50_3->SetLineColor(4);
  hfmuxlr50_3->SetMarkerStyle(21);
  hfmuxlr50_3->SetMarkerColor(4);
  hfmuxlr50_3->Draw("same P");
  hfmuxlr100_3->SetLineColor(6);
  hfmuxlr100_3->SetMarkerStyle(22);
  hfmuxlr100_3->SetMarkerColor(6);
  hfmuxlr100_3->Draw("same P");
	
//	l3->Draw();

  TLegend* leg2=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
  TLegendEntry* ent2a=leg2->AddEntry(hfmuxlr20_3,"|#muL - #muR| < 20 #mum","PL"); // P: marker, L; line
  ent2a->SetTextColor(hfmuxlr20_3->GetMarkerColor());
  TLegendEntry* ent2b=leg2->AddEntry(hfmuxlr50_3,"|#muL - #muR| < 50 #mum","PL"); // P: marker, L; line
  ent2b->SetTextColor(hfmuxlr50_3->GetMarkerColor());
  TLegendEntry* ent2c=leg2->AddEntry(hfmuxlr100_3,"|#muL - #muR| < 100 #mum","PL"); // P: marker, L; line
  ent2c->SetTextColor(hfmuxlr100_3->GetMarkerColor());
  leg2->SetFillStyle(0);
  leg2->Draw();

	c1->SaveAs("SDD_resx_L3");
//	pdfFileNames+=" SDD_resx_L3";
	c1->Update();

  TCanvas* c2=new TCanvas("c2","SDD Layer 4", 1200,800);
  c2->Divide(1,2);
  c2->cd(1);
  hfmux20_4->SetTitle("fraction(| #mu x-residual | < threshold) - layer 4");
  hfmux20_4->SetLineColor(2);
  hfmux20_4->SetMarkerStyle(31);
  hfmux20_4->SetMarkerSize(2);
  hfmux20_4->SetMarkerColor(2);
  hfmux20_4->SetMinimum(0.);
  hfmux20_4->SetMaximum(1.2);
  hfmux20_4->GetXaxis()->SetTitle("run number");
  hfmux20_4->GetYaxis()->SetTitleOffset(1.2);
  //  hfmux20_4->GetYaxis()->SetTitle("%");
  hfmux20_4->Draw("P");
  hfmux50_4->SetLineColor(4);
  hfmux50_4->SetMarkerStyle(21);
  hfmux50_4->SetMarkerColor(4);
  hfmux50_4->Draw("same P");
  hfmux100_4->SetLineColor(6);
  hfmux100_4->SetMarkerStyle(22);
  hfmux100_4->SetMarkerColor(6);
  hfmux100_4->Draw("same P");
	TLine *l4 = new TLine(0.,0.839,kRunsToPlot,0.839);
	l4->SetLineStyle(2);
	l4->SetLineWidth(2);	
//	l4->Draw();
	
  TLegend* legg=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
  TLegendEntry* entta=legg->AddEntry(hfmux20_4,"|#mu res x | < 20 #mum","PL"); // P: marker, L; line
  entta->SetTextColor(hfmux20_4->GetMarkerColor());
  TLegendEntry* enttb=legg->AddEntry(hfmux50_4,"|#mu res x | < 50 #mum","PL"); // P: marker, L; line
  enttb->SetTextColor(hfmux50_4->GetMarkerColor());
  TLegendEntry* enttc=legg->AddEntry(hfmux100_4,"|#mu res x | < 100 #mum","PL"); // P: marker, L; line
  enttc->SetTextColor(hfmux100_4->GetMarkerColor());
  legg->SetFillStyle(0);
  legg->Draw();
  //  c2->Update();

  c2->cd(2);
  hfmuxlr20_4->SetTitle("fraction(| #muL - #muR | < threshold) - layer 4");
  hfmuxlr20_4->SetLineColor(2);
  hfmuxlr20_4->SetMarkerStyle(31);
  hfmuxlr20_4->SetMarkerSize(2);
  hfmuxlr20_4->SetMarkerColor(2);
  hfmuxlr20_4->SetMinimum(0.);
  hfmuxlr20_4->SetMaximum(1.2);
  hfmuxlr20_4->GetXaxis()->SetTitle("run number");
  hfmuxlr20_4->GetYaxis()->SetTitleOffset(1.2);
  //  hfmuxlr20_4->GetYaxis()->SetTitle("%");
  hfmuxlr20_4->Draw("P");
  hfmuxlr50_4->SetLineColor(4);
  hfmuxlr50_4->SetMarkerStyle(21);
  hfmuxlr50_4->SetMarkerColor(4);
  hfmuxlr50_4->Draw("same P");
  hfmuxlr100_4->SetLineColor(6);
  hfmuxlr100_4->SetMarkerStyle(22);
  hfmuxlr100_4->SetMarkerColor(6);
  hfmuxlr100_4->Draw("same P");

//	l4->Draw();
	
  TLegend* legg2=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
  TLegendEntry* entt2a=legg2->AddEntry(hfmuxlr20_4,"|#muL - #muR| < 20 #mum","PL"); // P: marker, L; line
  entt2a->SetTextColor(hfmuxlr20_4->GetMarkerColor());
  TLegendEntry* entt2b=legg2->AddEntry(hfmuxlr50_4,"|#muL - #muR| < 50 #mum","PL"); // P: marker, L; line
  entt2b->SetTextColor(hfmuxlr50_4->GetMarkerColor());
  TLegendEntry* entt2c=legg2->AddEntry(hfmuxlr100_4,"|#muL - #muR| < 100 #mum","PL"); // P: marker, L; line
  entt2c->SetTextColor(hfmuxlr100_4->GetMarkerColor());
  legg2->SetFillStyle(0);
  legg2->Draw();

	c2->SaveAs("SDD_resx_L4");
//	pdfFileNames+=" SDD_resx_L4";
    c2->Update();
	
  //////////qui

 TCanvas* c3=new TCanvas("c3","SDD z residuals", 1200,800);
  c3->Divide(1,2);
  c3->cd(1);
  hfmuz50_3->SetTitle("fraction(| #mu z-residual | < threshold) - layer 3");
  hfmuz50_3->SetLineColor(2);
  hfmuz50_3->SetMarkerStyle(31);
  hfmuz50_3->SetMarkerSize(2);
  hfmuz50_3->SetMarkerColor(2);
  hfmuz50_3->SetMinimum(0.);
  hfmuz50_3->SetMaximum(1.2);
  hfmuz50_3->GetXaxis()->SetTitle("run number");
  hfmuz50_3->GetYaxis()->SetTitleOffset(1.2);
  //  hfmuz50_3->GetYaxis()->SetTitle("%");
  hfmuz50_3->Draw("P");
  hfmuz100_3->SetLineColor(4);
  hfmuz100_3->SetMarkerStyle(21);
  hfmuz100_3->SetMarkerColor(4);
  hfmuz100_3->Draw("same P");
  hfmuz300_3->SetLineColor(6);
  hfmuz300_3->SetMarkerStyle(22);
  hfmuz300_3->SetMarkerColor(6);
  hfmuz300_3->Draw("same P");
//	l3->Draw();
  TLegend* lege=new TLegend(0.5,0.15,0.88,0.35); // x1, y1, x2, y2
  TLegendEntry* ent_a=lege->AddEntry(hfmuz50_3,"|#mu res z| < 30 #mum","PL"); // P: marker, L; line
  ent_a->SetTextColor(hfmuz50_3->GetMarkerColor());
  TLegendEntry* ent_b=lege->AddEntry(hfmuz100_3,"|#mu res z| < 50 #mum","PL"); // P: marker, L; line
  ent_b->SetTextColor(hfmuz100_3->GetMarkerColor());
  TLegendEntry* ent_c=lege->AddEntry(hfmuz300_3,"|#mu res z| < 100 #mum","PL"); // P: marker, L; line
  ent_c->SetTextColor(hfmuz300_3->GetMarkerColor());
  lege->SetFillStyle(0);
  lege->Draw();
  //  c3->Update();

  c3->cd(2);
  hfmuz50_4->SetTitle("fraction(| #mu z-residual | < threshold) - layer 4");
  hfmuz50_4->SetLineColor(2);
  hfmuz50_4->SetMarkerStyle(31);
  hfmuz50_4->SetMarkerSize(2);
  hfmuz50_4->SetMarkerColor(2);
  hfmuz50_4->SetMinimum(0.);
  hfmuz50_4->SetMaximum(1.2);
  hfmuz50_4->GetXaxis()->SetTitle("run number");
  hfmuz50_4->GetYaxis()->SetTitleOffset(1.2);
  //  hfmuz50_4->GetYaxis()->SetTitle("%");
  hfmuz50_4->Draw("P");
  hfmuz100_4->SetLineColor(4);
  hfmuz100_4->SetMarkerStyle(21);
  hfmuz100_4->SetMarkerColor(4);
  hfmuz100_4->Draw("same P");
  hfmuz300_4->SetLineColor(6);
  hfmuz300_4->SetMarkerStyle(22);
  hfmuz300_4->SetMarkerColor(6);
  hfmuz300_4->Draw("same P");
//	l4->Draw();
  TLegend* leggi=new TLegend(0.5,0.15,0.88,0.35); // x1, y1, x2, y2
  TLegendEntry* en_ta=leggi->AddEntry(hfmuz50_4,"|#mu res z| < 30 #mum","PL"); // P: marker, L; line
  en_ta->SetTextColor(hfmuz50_4->GetMarkerColor());
  TLegendEntry* en_tb=leggi->AddEntry(hfmuz100_4,"|#mu res z| < 50 #mum","PL"); // P: marker, L; line
  en_tb->SetTextColor(hfmuz100_4->GetMarkerColor());
  TLegendEntry* en_tc=leggi->AddEntry(hfmuz300_4,"|#mu res z| < 100 #mum","PL"); // P: marker, L; line
  en_tc->SetTextColor(hfmuz300_4->GetMarkerColor());
  leggi->SetFillStyle(0);
  leggi->Draw();
	c3->SaveAs("SDD_resz");
//	pdfFileNames+=" SDD_resz";
  c3->Update();

  //////////////////////////////////////////////// SPD ////////////////////////////////////////////////////

 TCanvas* c4=new TCanvas("c4","SPD x residuals", 1200,800);
  c4->Divide(1,2);
  c4->cd(1);
  hfmux20_1->SetTitle("fraction(| #mu x-residual | < threshold) - layer 1");
  hfmux20_1->SetLineColor(2);
  hfmux20_1->SetMarkerStyle(31);
  hfmux20_1->SetMarkerSize(2);
  hfmux20_1->SetMarkerColor(2);
  hfmux20_1->SetMinimum(0.);
  hfmux20_1->SetMaximum(1.2);
  hfmux20_1->GetXaxis()->SetTitle("run number");
  hfmux20_1->GetYaxis()->SetTitleOffset(1.2);
  //  hfmux20_1->GetYaxis()->SetTitle("%");
  hfmux20_1->Draw("P");
  hfmux50_1->SetLineColor(4);
  hfmux50_1->SetMarkerStyle(21);
  hfmux50_1->SetMarkerColor(4);
  hfmux50_1->Draw("same P");
  hfmux100_1->SetLineColor(6);
  hfmux100_1->SetMarkerStyle(22);
  hfmux100_1->SetMarkerColor(6);
  hfmux100_1->Draw("same P");
	TLine *l1 = new TLine(0.,0.625,kRunsToPlot,0.625);
	l1->SetLineStyle(2);
	l1->SetLineWidth(2);	
//	l1->Draw();
	
  TLegend* lege1=new TLegend(0.5,0.15,0.88,0.35); // x1, y1, x2, y2
  TLegendEntry* ent_a1=lege1->AddEntry(hfmux20_1,"|#mu res x| < 5 #mum","PL"); // P: marker, L; line
  ent_a1->SetTextColor(hfmux20_1->GetMarkerColor());
  TLegendEntry* ent_b1=lege1->AddEntry(hfmux50_1,"|#mu res x| < 10 #mum","PL"); // P: marker, L; line
  ent_b1->SetTextColor(hfmux50_1->GetMarkerColor());
  TLegendEntry* ent_c1=lege1->AddEntry(hfmux100_1,"|#mu res x| < 50 #mum","PL"); // P: marker, L; line
  ent_c1->SetTextColor(hfmux100_1->GetMarkerColor());
  lege1->SetFillStyle(0);
  lege1->Draw();
  //  c4->Update();

  c4->cd(2);
  hfmux20_2->SetTitle("fraction(| #mu x-residual | < threshold) - layer 2");
  hfmux20_2->SetLineColor(2);
  hfmux20_2->SetMarkerStyle(31);
  hfmux20_2->SetMarkerSize(2);
  hfmux20_2->SetMarkerColor(2);
  hfmux20_2->SetMinimum(0.);
  hfmux20_2->SetMaximum(1.2);
  hfmux20_2->GetXaxis()->SetTitle("run number");
  hfmux20_2->GetYaxis()->SetTitleOffset(1.2);
  //  hfmux20_2->GetYaxis()->SetTitle("%");
  hfmux20_2->Draw("P");
  hfmux50_2->SetLineColor(4);
  hfmux50_2->SetMarkerStyle(21);
  hfmux50_2->SetMarkerColor(4);
  hfmux50_2->Draw("same P");
  hfmux100_2->SetLineColor(6);
  hfmux100_2->SetMarkerStyle(22);
  hfmux100_2->SetMarkerColor(6);
  hfmux100_2->Draw("same P");
	TLine *l2 = new TLine(0.,0.631,kRunsToPlot,0.631);
	l2->SetLineStyle(2);
	l2->SetLineWidth(2);	
//	l2->Draw();
	
  TLegend* leggi1=new TLegend(0.5,0.15,0.88,0.35); // x1, y1, x2, y2
  TLegendEntry* en_ta1=leggi1->AddEntry(hfmux20_2,"|#mu res x| < 10 #mum","PL"); // P: marker, L; line
  en_ta1->SetTextColor(hfmux20_2->GetMarkerColor());
  TLegendEntry* en_tb1=leggi1->AddEntry(hfmux50_2,"|#mu res x| < 20 #mum","PL"); // P: marker, L; line
  en_tb1->SetTextColor(hfmux50_2->GetMarkerColor());
  TLegendEntry* en_tc1=leggi1->AddEntry(hfmux100_2,"|#mu res x| < 50 #mum","PL"); // P: marker, L; line
  en_tc1->SetTextColor(hfmux100_2->GetMarkerColor());
  leggi1->SetFillStyle(0);
  leggi1->Draw();
	c4->SaveAs("SPD_resx");
//	pdfFileNames+=" SPD_resx";
  c4->Update();

 TCanvas* c5=new TCanvas("c5","SPD z residuals", 1200,800);
  c5->Divide(1,2);
  c5->cd(1);
  hfmuz50_1->SetTitle("fraction(| #mu z-residual | < threshold) - layer 1");
  hfmuz50_1->SetLineColor(2);
  hfmuz50_1->SetMarkerStyle(31);
  hfmuz50_1->SetMarkerColor(2);
  hfmuz50_1->SetMinimum(0.);
  hfmuz50_1->SetMaximum(1.2);
  hfmuz50_1->GetXaxis()->SetTitle("run number");
  hfmuz50_1->GetYaxis()->SetTitleOffset(1.2);
  //  hfmuz50_1->GetYaxis()->SetTitle("%");
  hfmuz50_1->Draw("P");
  hfmuz100_1->SetLineColor(4);
  hfmuz100_1->SetMarkerStyle(21);
  hfmuz100_1->SetMarkerColor(4);
  hfmuz100_1->Draw("same P");
  hfmuz300_1->SetLineColor(6);
  hfmuz300_1->SetMarkerStyle(22);
  hfmuz300_1->SetMarkerColor(6);
  hfmuz300_1->Draw("same P");
//	l1->Draw();

  TLegend* lege2=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
  TLegendEntry* ent_a2=lege2->AddEntry(hfmuz50_1,"|#mu res z| < 10 #mum","PL"); // P: marker, L; line
  ent_a2->SetTextColor(hfmuz50_1->GetMarkerColor());
  TLegendEntry* ent_b2=lege2->AddEntry(hfmuz100_1,"|#mu res z| < 20 #mum","PL"); // P: marker, L; line
  ent_b2->SetTextColor(hfmuz100_1->GetMarkerColor());
  TLegendEntry* ent_c2=lege2->AddEntry(hfmuz300_1,"|#mu res z| < 50 #mum","PL"); // P: marker, L; line
  ent_c2->SetTextColor(hfmuz300_1->GetMarkerColor());
  lege2->SetFillStyle(0);
  lege2->Draw();
  //  c3->Update();

  c5->cd(2);
  hfmuz50_2->SetTitle("fraction(| #mu z-residual | < threshold) - layer 2");
  hfmuz50_2->SetLineColor(2);
  hfmuz50_2->SetMarkerStyle(31);
  hfmuz50_2->SetMarkerSize(2);
  hfmuz50_2->SetMarkerColor(2);
  hfmuz50_2->SetMinimum(0.);
  hfmuz50_2->SetMaximum(1.2);
  hfmuz50_2->GetXaxis()->SetTitle("run number");
  hfmuz50_2->GetYaxis()->SetTitleOffset(1.2);
  //  hfmuz50_2->GetYaxis()->SetTitle("%");
  hfmuz50_2->Draw("P");
  hfmuz100_2->SetLineColor(4);
  hfmuz100_2->SetMarkerStyle(21);
  hfmuz100_2->SetMarkerColor(4);
  hfmuz100_2->Draw("same P");
  hfmuz300_2->SetLineColor(6);
  hfmuz300_2->SetMarkerStyle(22);
  hfmuz300_2->SetMarkerColor(6);
  hfmuz300_2->Draw("same P");
//	l2->Draw();

  TLegend* leggi2=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
  TLegendEntry* en_ta2=leggi2->AddEntry(hfmuz50_2,"|#mu res z| < 10 #mum","PL"); // P: marker, L; line
  en_ta2->SetTextColor(hfmuz50_2->GetMarkerColor());
  TLegendEntry* en_tb2=leggi2->AddEntry(hfmuz100_2,"|#mu res z| < 20 #mum","PL"); // P: marker, L; line
  en_tb2->SetTextColor(hfmuz100_2->GetMarkerColor());
  TLegendEntry* en_tc2=leggi2->AddEntry(hfmuz300_2,"|#mu res z| < 50 #mum","PL"); // P: marker, L; line
  en_tc2->SetTextColor(hfmuz300_2->GetMarkerColor());
  leggi2->SetFillStyle(0);
  leggi2->Draw();
	c5->SaveAs("SPD_resz");
//	pdfFileNames+=" SPD_resz";
  c5->Update();

  //////////////////////////////////////////////// SSD ////////////////////////////////////////////////////

 TCanvas* c6=new TCanvas("c6","SSD x residuals", 1200,800);
  c6->Divide(1,2);
  c6->cd(1);
  hfmux20_5->SetTitle("fraction(| #mu x-residual | < threshold) - layer 5");
  hfmux20_5->SetLineColor(2);
  hfmux20_5->SetMarkerStyle(31);
  hfmux20_5->SetMarkerSize(2);
  hfmux20_5->SetMarkerColor(2);
  hfmux20_5->SetMinimum(0.);
  hfmux20_5->SetMaximum(1.2);
  hfmux20_5->GetXaxis()->SetTitle("run number");
  hfmux20_5->GetYaxis()->SetTitleOffset(1.2);
  //  hfmux20_5->GetYaxis()->SetTitle("%");
  hfmux20_5->Draw("P");
  hfmux50_5->SetLineColor(4);
  hfmux50_5->SetMarkerStyle(21);
  hfmux50_5->SetMarkerColor(4);
  hfmux50_5->Draw("same P");
  hfmux100_5->SetLineColor(6);
  hfmux100_5->SetMarkerStyle(22);
  hfmux100_5->SetMarkerColor(6);
  hfmux100_5->Draw("same P");
	TLine *l5 = new TLine(0.,0.915,kRunsToPlot,0.915);
	l5->SetLineStyle(2);
	l5->SetLineWidth(2);	
//	l5->Draw();
	
  TLegend* lege3=new TLegend(0.5,0.15,0.88,0.35); // x1, y1, x2, y2
  TLegendEntry* ent_a3=lege3->AddEntry(hfmux20_5,"|#mu res x| < 20 #mum","PL"); // P: marker, L; line
  ent_a3->SetTextColor(hfmux20_5->GetMarkerColor());
  TLegendEntry* ent_b3=lege3->AddEntry(hfmux50_5,"|#mu res x| < 40 #mum","PL"); // P: marker, L; line
  ent_b3->SetTextColor(hfmux50_5->GetMarkerColor());
  TLegendEntry* ent_c3=lege3->AddEntry(hfmux100_5,"|#mu res x| < 100 #mum","PL"); // P: marker, L; line
  ent_c3->SetTextColor(hfmux100_5->GetMarkerColor());
  lege3->SetFillStyle(0);
  lege3->Draw();
  //  c6->Update();

  c6->cd(2);
  hfmux20_6->SetTitle("fraction(| #mu x-residual | < threshold) - layer 6");
  hfmux20_6->SetLineColor(2);
  hfmux20_6->SetMarkerStyle(31);
  hfmux20_6->SetMarkerSize(2);
  hfmux20_6->SetMarkerColor(2);
  hfmux20_6->SetMinimum(0.);
  hfmux20_6->SetMaximum(1.2);
  hfmux20_6->GetXaxis()->SetTitle("run number");
  hfmux20_6->GetYaxis()->SetTitleOffset(1.2);
  //  hfmux20_6->GetYaxis()->SetTitle("%");
  hfmux20_6->Draw("P");
  hfmux50_6->SetLineColor(4);
  hfmux50_6->SetMarkerStyle(21);
  hfmux50_6->SetMarkerColor(4);
  hfmux50_6->Draw("same P");
  hfmux100_6->SetLineColor(6);
  hfmux100_6->SetMarkerStyle(22);
  hfmux100_6->SetMarkerColor(6);
  hfmux100_6->Draw("same P");
	TLine *l6 = new TLine(0.,0.881,kRunsToPlot,0.881);
	l6->SetLineStyle(2);
	l6->SetLineWidth(2);	
//	l6->Draw();
	
  TLegend* leggi3=new TLegend(0.5,0.15,0.88,0.35); // x1, y1, x2, y2
  TLegendEntry* en_ta3=leggi3->AddEntry(hfmux20_6,"|#mu res x| < 20 #mum","PL"); // P: marker, L; line
  en_ta3->SetTextColor(hfmux20_6->GetMarkerColor());
  TLegendEntry* en_tb3=leggi3->AddEntry(hfmux50_6,"|#mu res x| < 40 #mum","PL"); // P: marker, L; line
  en_tb3->SetTextColor(hfmux50_6->GetMarkerColor());
  TLegendEntry* en_tc3=leggi3->AddEntry(hfmux100_6,"|#mu res x| < 100 #mum","PL"); // P: marker, L; line
  en_tc3->SetTextColor(hfmux100_6->GetMarkerColor());
  leggi3->SetFillStyle(0);
  leggi3->Draw();
	c6->SaveAs("SSD_resx");
//	pdfFileNames+=" SSD_resx";
  c6->Update();

 TCanvas* c7=new TCanvas("c7","SSD z residuals", 1200,800);
  c7->Divide(1,2);
  c7->cd(1);
  hfmuz50_5->SetTitle("fraction(| #mu z-residual | < threshold) - layer 5");
  hfmuz50_5->SetLineColor(2);
  hfmuz50_5->SetMarkerStyle(31);
  hfmuz50_5->SetMarkerSize(2);
  hfmuz50_5->SetMarkerColor(2);
  hfmuz50_5->SetMinimum(0.);
  hfmuz50_5->SetMaximum(1.2);
  hfmuz50_5->GetXaxis()->SetTitle("run number");
  hfmuz50_5->GetYaxis()->SetTitleOffset(1.2);
  //  hfmuz50_5->GetYaxis()->SetTitle("%");
  hfmuz50_5->Draw("P");
  hfmuz100_5->SetLineColor(4);
  hfmuz100_5->SetMarkerStyle(21);
  hfmuz100_5->SetMarkerColor(4);
  hfmuz100_5->Draw("same P");
  hfmuz300_5->SetLineColor(6);
  hfmuz300_5->SetMarkerStyle(22);
  hfmuz300_5->SetMarkerColor(6);
  hfmuz300_5->Draw("same P");
//	l5->Draw();

  TLegend* lege4=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
  TLegendEntry* ent_a4=lege4->AddEntry(hfmuz50_5,"|#mu res z| < 50 #mum","PL"); // P: marker, L; line
  ent_a4->SetTextColor(hfmuz50_5->GetMarkerColor());
  TLegendEntry* ent_b4=lege4->AddEntry(hfmuz100_5,"|#mu res z| < 100 #mum","PL"); // P: marker, L; line
  ent_b4->SetTextColor(hfmuz100_5->GetMarkerColor());
  TLegendEntry* ent_c4=lege4->AddEntry(hfmuz300_5,"|#mu res z| < 300 #mum","PL"); // P: marker, L; line
  ent_c4->SetTextColor(hfmuz300_5->GetMarkerColor());
  lege4->SetFillStyle(0);
  lege4->Draw();
  //  c7->Update();

  c7->cd(2);
  hfmuz50_6->SetTitle("fraction(| #mu z-residual | < threshold) - layer 6");
  hfmuz50_6->SetLineColor(2);
  hfmuz50_6->SetMarkerStyle(31);
  hfmuz50_6->SetMarkerSize(2);
  hfmuz50_6->SetMarkerColor(2);
  hfmuz50_6->SetMinimum(0.);
  hfmuz50_6->SetMaximum(1.2);
  hfmuz50_6->GetXaxis()->SetTitle("run number");
  hfmuz50_6->GetYaxis()->SetTitleOffset(1.2);
  //  hfmuz50_6->GetYaxis()->SetTitle("%");
  hfmuz50_6->Draw("P");
  hfmuz100_6->SetLineColor(4);
  hfmuz100_6->SetMarkerStyle(21);
  hfmuz100_6->SetMarkerColor(4);
  hfmuz100_6->Draw("same P");
  hfmuz300_6->SetLineColor(6);
  hfmuz300_6->SetMarkerStyle(22);
  hfmuz300_6->SetMarkerColor(6);
  hfmuz300_6->Draw("same P");
//	l6->Draw();
  TLegend* leggi4=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
  TLegendEntry* en_ta4=leggi4->AddEntry(hfmuz50_6,"|#mu res z| < 50 #mum","PL"); // P: marker, L; line
  en_ta4->SetTextColor(hfmuz50_6->GetMarkerColor());
  TLegendEntry* en_tb4=leggi4->AddEntry(hfmuz100_6,"|#mu res z| < 100 #mum","PL"); // P: marker, L; line
  en_tb4->SetTextColor(hfmuz100_6->GetMarkerColor());
  TLegendEntry* en_tc4=leggi4->AddEntry(hfmuz300_6,"|#mu res z| < 300 #mum","PL"); // P: marker, L; line
  en_tc4->SetTextColor(hfmuz300_6->GetMarkerColor());
  leggi4->SetFillStyle(0);
  leggi4->Draw();
	c7->SaveAs("SSD_resz");
//	pdfFileNames+=" SSD_resz";
  c7->Update();
  

	// -------------------------- plot sigma --------------------------------------
	// -------------------------- SPD  x -------------------------------------
	
	TCanvas* c8=new TCanvas("c8","SPD x res width", 1200,800);
	c8->Divide(1,2);
	c8->cd(1);
	hfsigx100_1->SetTitle("fraction(#sigma x-residual < threshold) - layer 1");
	hfsigx100_1->SetLineColor(2);
	hfsigx100_1->SetMarkerStyle(31);
	hfsigx100_1->SetMarkerSize(2);
	hfsigx100_1->SetMarkerColor(2);
	hfsigx100_1->SetMinimum(0.);
	hfsigx100_1->SetMaximum(1.2);
	hfsigx100_1->GetXaxis()->SetTitle("run number");
	hfsigx100_1->GetYaxis()->SetTitleOffset(1.2);
	//  hfsigx100_1->GetYaxis()->SetTitle("%");
	hfsigx100_1->Draw("P");
	hfsigx200_1->SetLineColor(4);
	hfsigx200_1->SetMarkerStyle(21);
	hfsigx200_1->SetMarkerColor(4);
	hfsigx200_1->Draw("same P");
	hfsigx300_1->SetLineColor(6);
	hfsigx300_1->SetMarkerStyle(22);
	hfsigx300_1->SetMarkerColor(6);
	hfsigx300_1->Draw("same P");
//	l1->Draw();
	
	TLegend* lege2s=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
	TLegendEntry* ent_a2s=lege2s->AddEntry(hfsigx100_1,"|#sigma res x| < 50 #mum","PL"); // P: marker, L; line
	ent_a2s->SetTextColor(hfsigx100_1->GetMarkerColor());
	TLegendEntry* ent_b2s=lege2s->AddEntry(hfsigx200_1,"|#sigma res x| < 70 #mum","PL"); // P: marker, L; line
	ent_b2s->SetTextColor(hfsigx200_1->GetMarkerColor());
	TLegendEntry* ent_c2s=lege2s->AddEntry(hfsigx300_1,"|#sigma res x| < 200 #mum","PL"); // P: marker, L; line
	ent_c2s->SetTextColor(hfsigx300_1->GetMarkerColor());
	lege2s->SetFillStyle(0);
	lege2s->Draw();
	
	c8->cd(2);
	hfsigx100_2->SetTitle("fraction(#sigma x-residual < threshold) - layer 2");
	hfsigx100_2->SetLineColor(2);
	hfsigx100_2->SetMarkerStyle(31);
	hfsigx100_2->SetMarkerSize(2);
	hfsigx100_2->SetMarkerColor(2);
	hfsigx100_2->SetMinimum(0.);
	hfsigx100_2->SetMaximum(1.2);
	hfsigx100_2->GetXaxis()->SetTitle("run number");
	hfsigx100_2->GetYaxis()->SetTitleOffset(1.2);
	//  hfsigx100_2->GetYaxis()->SetTitle("%");
	hfsigx100_2->Draw("P");
	hfsigx200_2->SetLineColor(4);
	hfsigx200_2->SetMarkerStyle(21);
	hfsigx200_2->SetMarkerColor(4);
       	hfsigx200_2->Draw("same P");
	hfsigx300_2->SetLineColor(6);
	hfsigx300_2->SetMarkerStyle(22);
	hfsigx300_2->SetMarkerColor(6);
	hfsigx300_2->Draw("same P");
//	l2->Draw();
	
	TLegend* leggi2s=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
	TLegendEntry* en_ta2s=leggi2s->AddEntry(hfsigx100_2,"#sigma res x < 130 #mum","PL"); // P: marker, L; line
	en_ta2s->SetTextColor(hfsigx100_2->GetMarkerColor());
	TLegendEntry* en_tb2s=leggi2s->AddEntry(hfsigx200_2,"#sigma res x < 160 #mum","PL"); // P: marker, L; line
	en_tb2s->SetTextColor(hfsigx200_2->GetMarkerColor());
	TLegendEntry* en_tc2s=leggi2s->AddEntry(hfsigx300_2,"#sigma res x < 200 #mum","PL"); // P: marker, L; line
	en_tc2s->SetTextColor(hfsigx300_2->GetMarkerColor());
	leggi2s->SetFillStyle(0);
	leggi2s->Draw();
	c8->SaveAs("SPD_sigx");
//	pdfFileNames+=" SPD_sigx";
	c8->Update();
	
	// -------------------------- SPD  z -------------------------------------
	
	TCanvas* c9=new TCanvas("c9","SPD z res width", 1200,800);
	c9->Divide(1,2);
	c9->cd(1);
	hfsigz100_1->SetTitle("fraction(#sigma z-residual < threshold) - layer 1");
	hfsigz100_1->SetLineColor(2);
	hfsigz100_1->SetMarkerStyle(31);
	hfsigz100_1->SetMarkerSize(2);
	hfsigz100_1->SetMarkerColor(2);
	hfsigz100_1->SetMinimum(0.);
	hfsigz100_1->SetMaximum(1.2);
	hfsigz100_1->GetXaxis()->SetTitle("run number");
	hfsigz100_1->GetYaxis()->SetTitleOffset(1.2);
	//  hfsigz100_1->GetYaxis()->SetTitle("%");
	hfsigz100_1->Draw("P");
	hfsigz300_1->SetLineColor(4);
	hfsigz300_1->SetMarkerStyle(21);
	hfsigz300_1->SetMarkerColor(4);
	hfsigz300_1->Draw("same P");
	hfsigz500_1->SetLineColor(6);
	hfsigz500_1->SetMarkerStyle(22);
	hfsigz500_1->SetMarkerColor(6);
	hfsigz500_1->Draw("same P");
//	l1->Draw();
	TLegend* lege20s=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
	TLegendEntry* ent_a20s=lege20s->AddEntry(hfsigz100_1,"#sigma res z < 130 #mum","PL"); // P: marker, L; line
	ent_a20s->SetTextColor(hfsigz100_1->GetMarkerColor());
	TLegendEntry* ent_b20s=lege20s->AddEntry(hfsigz300_1,"#sigma res z < 160 #mum","PL"); // P: marker, L; line
	ent_b20s->SetTextColor(hfsigz300_1->GetMarkerColor());
	TLegendEntry* ent_c20s=lege20s->AddEntry(hfsigz500_1,"#sigma res z < 200 #mum","PL"); // P: marker, L; line
	ent_c20s->SetTextColor(hfsigz500_1->GetMarkerColor());
	lege20s->SetFillStyle(0);
	lege20s->Draw();
	
	c9->cd(2);
	hfsigz100_2->SetTitle("fraction(#sigma z-residual < threshold) - layer 2");
	hfsigz100_2->SetLineColor(2);
	hfsigz100_2->SetMarkerStyle(31);
	hfsigz100_2->SetMarkerSize(2);
	hfsigz100_2->SetMarkerColor(2);
	hfsigz100_2->SetMinimum(0.);
	hfsigz100_2->SetMaximum(1.2);
	hfsigz100_2->GetXaxis()->SetTitle("run number");
	hfsigz100_2->GetYaxis()->SetTitleOffset(1.2);
	//  hfsigz100_2->GetYaxis()->SetTitle("%");
	hfsigz100_2->Draw("P");
	hfsigz300_2->SetLineColor(4);
	hfsigz300_2->SetMarkerStyle(21);
	hfsigz300_2->SetMarkerColor(4);
	hfsigz300_2->Draw("same P");
	hfsigz500_2->SetLineColor(6);
	hfsigz500_2->SetMarkerStyle(22);
	hfsigz500_2->SetMarkerColor(6);
	hfsigz500_2->Draw("same P");
//	l2->Draw();
	TLegend* leggi21s=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
	TLegendEntry* en_ta21s=leggi21s->AddEntry(hfsigz100_2,"#sigma res z < 180 #mum","PL"); // P: marker, L; line
	en_ta21s->SetTextColor(hfsigz100_2->GetMarkerColor());
	TLegendEntry* en_tb21s=leggi21s->AddEntry(hfsigz300_2,"#sigma res z < 190 #mum","PL"); // P: marker, L; line
	en_tb21s->SetTextColor(hfsigz300_2->GetMarkerColor());
	TLegendEntry* en_tc21s=leggi21s->AddEntry(hfsigz500_2,"#sigma res z < 250 #mum","PL"); // P: marker, L; line
	en_tc21s->SetTextColor(hfsigz500_2->GetMarkerColor());
	leggi21s->SetFillStyle(0);
	leggi21s->Draw();
	c9->SaveAs("SPD_sigz");
//	pdfFileNames+=" SPD_sigz";
	c9->Update();

	//---------------------------- SDD x ----------------------------------------

	TCanvas* c10=new TCanvas("c10","SDD x res width", 1200,800);
	c10->Divide(1,2);
	c10->cd(1);
	hfsigx100_3->SetTitle("fraction(#sigma x-residual < threshold) - layer 3");
	hfsigx100_3->SetLineColor(2);
	hfsigx100_3->SetMarkerStyle(31);
	hfsigx100_3->SetMarkerSize(2);
	hfsigx100_3->SetMarkerColor(2);
	hfsigx100_3->SetMinimum(0.);
	hfsigx100_3->SetMaximum(1.2);
	hfsigx100_3->GetXaxis()->SetTitle("run number");
	hfsigx100_3->GetYaxis()->SetTitleOffset(1.2);
	//  hfsigx100_3->GetYaxis()->SetTitle("%");
	hfsigx100_3->Draw("P");
	hfsigx200_3->SetLineColor(4);
	hfsigx200_3->SetMarkerStyle(21);
	hfsigx200_3->SetMarkerColor(4);
	hfsigx200_3->Draw("same P");
	hfsigx300_3->SetLineColor(6);
	hfsigx300_3->SetMarkerStyle(22);
	hfsigx300_3->SetMarkerColor(6);
	hfsigx300_3->Draw("same P");
//	l3->Draw();
	TLegend* lege22s=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
	TLegendEntry* ent_a22s=lege22s->AddEntry(hfsigx100_3,"#sigma res x < 200 #mum","PL"); // P: marker, L; line
	ent_a22s->SetTextColor(hfsigx100_3->GetMarkerColor());
	TLegendEntry* ent_b22s=lege22s->AddEntry(hfsigx200_3,"#sigma res x < 250 #mum","PL"); // P: marker, L; line
	ent_b22s->SetTextColor(hfsigx200_3->GetMarkerColor());
	TLegendEntry* ent_c22s=lege22s->AddEntry(hfsigx300_3,"#sigma res x < 400 #mum","PL"); // P: marker, L; line
	ent_c22s->SetTextColor(hfsigx300_3->GetMarkerColor());
	lege22s->SetFillStyle(0);
	lege22s->Draw();
	
	c10->cd(2);
	hfsigx100_4->SetTitle("fraction(#sigma x-residual < threshold) - layer 4");
	hfsigx100_4->SetLineColor(2);
	hfsigx100_4->SetMarkerStyle(31);
	hfsigx100_4->SetMarkerSize(2);
	hfsigx100_4->SetMarkerColor(2);
	hfsigx100_4->SetMinimum(0.);
	hfsigx100_4->SetMaximum(1.2);
	hfsigx100_4->GetXaxis()->SetTitle("run number");
	hfsigx100_4->GetYaxis()->SetTitleOffset(1.2);
	//  hfsigx100_4->GetYaxis()->SetTitle("%");
	hfsigx100_4->Draw("P");
	hfsigx200_4->SetLineColor(4);
	hfsigx200_4->SetMarkerStyle(21);
	hfsigx200_4->SetMarkerColor(4);
	hfsigx200_4->Draw("same P");
	hfsigx300_4->SetLineColor(6);
	hfsigx300_4->SetMarkerStyle(22);
	hfsigx300_4->SetMarkerColor(6);
	hfsigx300_4->Draw("same P");
//	l4->Draw();
	TLegend* leggi23s=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
	TLegendEntry* en_ta23s=leggi23s->AddEntry(hfsigx100_4,"#sigma res x < 200 #mum","PL"); // P: marker, L; line
	en_ta23s->SetTextColor(hfsigx100_4->GetMarkerColor());
	TLegendEntry* en_tb23s=leggi23s->AddEntry(hfsigx200_4,"#sigma res x < 250 #mum","PL"); // P: marker, L; line
	en_tb23s->SetTextColor(hfsigx200_4->GetMarkerColor());
	TLegendEntry* en_tc23s=leggi23s->AddEntry(hfsigx300_4,"#sigma res x < 400 #mum","PL"); // P: marker, L; line
	en_tc23s->SetTextColor(hfsigx300_4->GetMarkerColor());
	leggi23s->SetFillStyle(0);
	leggi23s->Draw();
	c10->SaveAs("SDD_sigx");
//	pdfFileNames+=" SDD_sigx";
	c10->Update();
	
	// -------------------------- SDD  z -------------------------------------
	
	TCanvas* c11=new TCanvas("c11","SDD z res width", 1200,800);
	c11->Divide(1,2);
	c11->cd(1);
	hfsigz100_3->SetTitle("fraction(#sigma z-residual < threshold) - layer 3");
	hfsigz100_3->SetLineColor(2);
	hfsigz100_3->SetMarkerStyle(31);
	hfsigz100_3->SetMarkerSize(2);
	hfsigz100_3->SetMarkerColor(2);
	hfsigz100_3->SetMinimum(0.);
	hfsigz100_3->SetMaximum(1.2);
	hfsigz100_3->GetXaxis()->SetTitle("run number");
	hfsigz100_3->GetYaxis()->SetTitleOffset(1.2);
	//  hfsigz100_3->GetYaxis()->SetTitle("%");
	hfsigz100_3->Draw("P");
	hfsigz300_3->SetLineColor(4);
	hfsigz300_3->SetMarkerStyle(21);
	hfsigz300_3->SetMarkerColor(4);
	hfsigz300_3->Draw("same P");
	hfsigz500_3->SetLineColor(6);
	hfsigz500_3->SetMarkerStyle(22);
	hfsigz500_3->SetMarkerColor(6);
	hfsigz500_3->Draw("same P");
//	l3->Draw();
	TLegend* lege24s=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
	TLegendEntry* ent_a24s=lege24s->AddEntry(hfsigz100_3,"#sigma res z < 150 #mum","PL"); // P: marker, L; line
	ent_a24s->SetTextColor(hfsigz100_3->GetMarkerColor());
	TLegendEntry* ent_b24s=lege24s->AddEntry(hfsigz300_3,"#sigma res z < 200 #mum","PL"); // P: marker, L; line
	ent_b24s->SetTextColor(hfsigz300_3->GetMarkerColor());
	TLegendEntry* ent_c24s=lege24s->AddEntry(hfsigz500_3,"#sigma res z < 350 #mum","PL"); // P: marker, L; line
	ent_c24s->SetTextColor(hfsigz500_3->GetMarkerColor());
	lege24s->SetFillStyle(0);
	lege24s->Draw();
	
	c11->cd(2);
	hfsigz100_4->SetTitle("fraction(#sigma z-residual < threshold) - layer 4");
	hfsigz100_4->SetLineColor(2);
	hfsigz100_4->SetMarkerStyle(31);
	hfsigz100_4->SetMarkerSize(2);
	hfsigz100_4->SetMarkerColor(2);
	hfsigz100_4->SetMinimum(0.);
	hfsigz100_4->SetMaximum(1.2);
	hfsigz100_4->GetXaxis()->SetTitle("run number");
	hfsigz100_4->GetYaxis()->SetTitleOffset(1.2);
	//  hfsigz100_4->GetYaxis()->SetTitle("%");
	hfsigz100_4->Draw("P");
	hfsigz300_4->SetLineColor(4);
	hfsigz300_4->SetMarkerStyle(21);
	hfsigz300_4->SetMarkerColor(4);
	hfsigz300_4->Draw("same P");
	hfsigz500_4->SetLineColor(6);
	hfsigz500_4->SetMarkerStyle(22);
	hfsigz500_4->SetMarkerColor(6);
	hfsigz500_4->Draw("same P");
//	l4->Draw();
	TLegend* leggi25s=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
	TLegendEntry* en_ta25s=leggi25s->AddEntry(hfsigz100_4,"#sigma res z < 250 #mum","PL"); // P: marker, L; line
	en_ta25s->SetTextColor(hfsigz100_4->GetMarkerColor());
	TLegendEntry* en_tb25s=leggi25s->AddEntry(hfsigz300_4,"#sigma res z < 300 #mum","PL"); // P: marker, L; line
	en_tb25s->SetTextColor(hfsigz300_4->GetMarkerColor());
	TLegendEntry* en_tc25s=leggi25s->AddEntry(hfsigz500_4,"#sigma res | < 450 #mum","PL"); // P: marker, L; line
	en_tc25s->SetTextColor(hfsigz500_4->GetMarkerColor());
	leggi25s->SetFillStyle(0);
	leggi25s->Draw();
	c11->SaveAs("SDD_sigz");
//	pdfFileNames+=" SDD_sigz";
	c11->Update();

	//---------------------------- SSD x ----------------------------------------

	TCanvas* c12=new TCanvas("c12","SSD x res width", 1200,800);
	c12->Divide(1,2);
	c12->cd(1);
	hfsigx100_5->SetTitle("fraction(#sigma x-residual < threshold) - layer 5");
	hfsigx100_5->SetLineColor(2);
	hfsigx100_5->SetMarkerStyle(31);
	hfsigx100_5->SetMarkerSize(2);
	hfsigx100_5->SetMarkerColor(2);
	hfsigx100_5->SetMinimum(0.);
	hfsigx100_5->SetMaximum(1.2);
	hfsigx100_5->GetXaxis()->SetTitle("run number");
	hfsigx100_5->GetYaxis()->SetTitleOffset(1.2);
	//  hfsigx100_5->GetYaxis()->SetTitle("%");
	hfsigx100_5->Draw("P");
	hfsigx200_5->SetLineColor(4);
	hfsigx200_5->SetMarkerStyle(21);
	hfsigx200_5->SetMarkerColor(4);
	hfsigx200_5->Draw("same P");
	hfsigx300_5->SetLineColor(6);
	hfsigx300_5->SetMarkerStyle(22);
	hfsigx300_5->SetMarkerColor(6);
	hfsigx300_5->Draw("same P");
//	l5->Draw();
	TLegend* lege26s=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
	TLegendEntry* ent_a26s=lege26s->AddEntry(hfsigx100_5,"#sigma res x < 70 #mum","PL"); // P: marker, L; line
	ent_a26s->SetTextColor(hfsigx100_5->GetMarkerColor());
	TLegendEntry* ent_b26s=lege26s->AddEntry(hfsigx200_5,"#sigma res x < 100 #mum","PL"); // P: marker, L; line
	ent_b26s->SetTextColor(hfsigx200_5->GetMarkerColor());
	TLegendEntry* ent_c26s=lege26s->AddEntry(hfsigx300_5,"#sigma res x < 300 #mum","PL"); // P: marker, L; line
	ent_c26s->SetTextColor(hfsigx300_5->GetMarkerColor());
	lege26s->SetFillStyle(0);
	lege26s->Draw();
	
	c12->cd(2);
	hfsigx100_6->SetTitle("fraction(#sigma x-residual < threshold) - layer 6");
	hfsigx100_6->SetLineColor(2);
	hfsigx100_6->SetMarkerStyle(31);
	hfsigx100_6->SetMarkerSize(2);
	hfsigx100_6->SetMarkerColor(2);
	hfsigx100_6->SetMinimum(0.);
	hfsigx100_6->SetMaximum(1.2);
	hfsigx100_6->GetXaxis()->SetTitle("run number");
	hfsigx100_6->GetYaxis()->SetTitleOffset(1.2);
	//  hfsigx100_6->GetYaxis()->SetTitle("%");
	hfsigx100_6->Draw("P");
	hfsigx200_6->SetLineColor(4);
	hfsigx200_6->SetMarkerStyle(21);
	hfsigx200_6->SetMarkerColor(4);
	hfsigx200_6->Draw("same P");
	hfsigx300_6->SetLineColor(6);
	hfsigx300_6->SetMarkerStyle(22);
	hfsigx300_6->SetMarkerColor(6);
	hfsigx300_6->Draw("same P");
//	l6->Draw();
	TLegend* leggi27s=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
	TLegendEntry* en_ta27s=leggi27s->AddEntry(hfsigx100_6,"#sigma res x < 200 #mum","PL"); // P: marker, L; line
	en_ta27s->SetTextColor(hfsigx100_6->GetMarkerColor());
	TLegendEntry* en_tb27s=leggi27s->AddEntry(hfsigx200_6,"#sigma res x < 300 #mum","PL"); // P: marker, L; line
	en_tb27s->SetTextColor(hfsigx200_6->GetMarkerColor());
	TLegendEntry* en_tc27s=leggi27s->AddEntry(hfsigx300_6,"#sigma res x < 400 #mum","PL"); // P: marker, L; line
	en_tc27s->SetTextColor(hfsigx300_6->GetMarkerColor());
	leggi27s->SetFillStyle(0);
	leggi27s->Draw();
	c12->SaveAs("SSD_sigx");
//	pdfFileNames+=" SSD_sigx";
	c12->Update();
	
	// -------------------------- SSD  z -------------------------------------
	
	TCanvas* c13=new TCanvas("c13","SSD z res width", 1200,800);
	c13->Divide(1,2);
	c13->cd(1);
	hfsigz100_5->SetTitle("fraction(#sigma z-residual < threshold) - layer 5");
	hfsigz100_5->SetLineColor(2);
	hfsigz100_5->SetMarkerStyle(31);
	hfsigz100_5->SetMarkerSize(2);
	hfsigz100_5->SetMarkerColor(2);
	hfsigz100_5->SetMinimum(0.);
	hfsigz100_5->SetMaximum(1.2);
	hfsigz100_5->GetXaxis()->SetTitle("run number");
	hfsigz100_5->GetYaxis()->SetTitleOffset(1.2);
	//  hfsigz100_5->GetYaxis()->SetTitle("%");
	hfsigz100_5->Draw("P");
	hfsigz300_5->SetLineColor(4);
	hfsigz300_5->SetMarkerStyle(21);
	hfsigz300_5->SetMarkerColor(4);
	hfsigz300_5->Draw("same P");
	hfsigz500_5->SetLineColor(6);
	hfsigz500_5->SetMarkerStyle(22);
	hfsigz500_5->SetMarkerColor(6);
	hfsigz500_5->Draw("same P");
//	l5->Draw();
	TLegend* lege28s=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
	TLegendEntry* ent_a28s=lege28s->AddEntry(hfsigz100_5,"#sigma res z < 1000 #mum","PL"); // P: marker, L; line
	ent_a28s->SetTextColor(hfsigz100_5->GetMarkerColor());
	TLegendEntry* ent_b28s=lege28s->AddEntry(hfsigz300_5,"#sigma res z < 1100 #mum","PL"); // P: marker, L; line
	ent_b28s->SetTextColor(hfsigz300_5->GetMarkerColor());
	TLegendEntry* ent_c28s=lege28s->AddEntry(hfsigz500_5,"#sigma res z < 1500 #mum","PL"); // P: marker, L; line
	ent_c28s->SetTextColor(hfsigz500_5->GetMarkerColor());
	lege28s->SetFillStyle(0);
	lege28s->Draw();
	
	c13->cd(2);
	hfsigz100_6->SetTitle("fraction(#sigma z-residual < threshold) - layer 6");
	hfsigz100_6->SetLineColor(2);
	hfsigz100_6->SetMarkerStyle(31);
	hfsigz100_6->SetMarkerSize(2);
	hfsigz100_6->SetMarkerColor(2);
	hfsigz100_6->SetMinimum(0.);
	hfsigz100_6->SetMaximum(1.2);
	hfsigz100_6->GetXaxis()->SetTitle("run number");
	hfsigz100_6->GetYaxis()->SetTitleOffset(1.2);
	//  hfsigz100_6->GetYaxis()->SetTitle("%");
	hfsigz100_6->Draw("P");
	hfsigz300_6->SetLineColor(4);
	hfsigz300_6->SetMarkerStyle(21);
	hfsigz300_6->SetMarkerColor(4);
	hfsigz300_6->Draw("same P");
	hfsigz500_6->SetLineColor(6);
	hfsigz500_6->SetMarkerStyle(22);
	hfsigz500_6->SetMarkerColor(6);
	hfsigz500_6->Draw("same P");
//	l6->Draw();
	TLegend* leggi29s=new TLegend(0.5,0.80,0.88,1.00); // x1, y1, x2, y2
	TLegendEntry* en_ta29s=leggi29s->AddEntry(hfsigz100_6,"#sigma res z < 1000 #mum","PL"); // P: marker, L; line
	en_ta29s->SetTextColor(hfsigz100_6->GetMarkerColor());
	TLegendEntry* en_tb29s=leggi29s->AddEntry(hfsigz300_6,"#sigma res z < 1100 #mum","PL"); // P: marker, L; line
	en_tb29s->SetTextColor(hfsigz300_6->GetMarkerColor());
	TLegendEntry* en_tc29s=leggi29s->AddEntry(hfsigz500_6,"#sigma res z  < 1500 #mum","PL"); // P: marker, L; line
	en_tc29s->SetTextColor(hfsigz500_6->GetMarkerColor());
	leggi29s->SetFillStyle(0);
	leggi29s->Draw();
	c13->SaveAs("SSD_sigz");
//	pdfFileNames+=" SSD_sigz";
	c13->Update();
  
   // order pdf files
	pdfFileNames+=" SPD_resx";
	pdfFileNames+=" SPD_sigx";
	pdfFileNames+=" SPD_resz";
	pdfFileNames+=" SPD_sigz";
	pdfFileNames+=" SDD_resx_L3";
	pdfFileNames+=" SDD_resx_L4";
	pdfFileNames+=" SDD_sigx";
	pdfFileNames+=" SDD_resz";
	pdfFileNames+=" SDD_sigz";
	pdfFileNames+=" SSD_resx";
	pdfFileNames+=" SSD_sigx";
	pdfFileNames+=" SSD_resz";
	pdfFileNames+=" SSD_sigz";
	
   // merge the pdf files
  TString command("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=merged");
  command=command+"ITS_Align_trend_norm.pdf "+pdfFileNames;
  gSystem->Exec(command.Data());
  printf(" Merging the pdf file:  %s \n",command.Data());
  delete [] myIndex;
  delete [] noRuns;

}
