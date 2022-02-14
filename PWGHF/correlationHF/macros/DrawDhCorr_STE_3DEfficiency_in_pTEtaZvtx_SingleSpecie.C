

//Macro to create tracking 3D (pt, eta, Zvtx) efficiency file for Dh Corr analysis
//Jitendra

#include <iostream>
#include "TSystem.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "THnSparse.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH3F.h"
#include "TH2F.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCFContainer.h"
#include "AliCFGridSparse.h"
#endif

using namespace std;


void CheckFileQuality(TString specie="pi");
void DrawDhCorr_STE_3DEfficiency_in_pTEtaZvtx(TString filename="AnalysisResults.root", TString specie="All") { //Species: "All", "pi", "K", "p", "e", "mu", "5Species", "SigP", "SigM", "Rest"
    
    //New bins for final rebinning
    const Int_t newLimitsBins = 22;
    Double_t* newLimits = new Double_t[newLimitsBins+1];
    newLimits[0]   = 0.30;
    newLimits[1]   = 0.40;
    newLimits[2]   = 0.50;
    newLimits[3]   = 0.60;
    newLimits[4]   = 0.70;
    newLimits[5]   = 0.80;
    newLimits[6]   = 0.90;
    newLimits[7]   = 1.00;
    newLimits[8]   = 1.20;
    newLimits[9]   = 1.40;
    newLimits[10]  = 1.60;
    newLimits[11]  = 1.80;
    newLimits[12]  = 2.00;
    newLimits[13]  = 2.40;
    newLimits[14]  = 2.80;
    newLimits[15]  = 3.40;
    newLimits[16]  = 4.00;
    newLimits[17]  = 4.80;
    newLimits[18]  = 6.00;
    newLimits[19]  = 8.00;
    newLimits[20]  = 10.00;
    newLimits[21]  = 16.00;
    newLimits[22]  = 24.00;

    TFile* f = new TFile(filename.Data(),"read");
    TCanvas *c = new TCanvas("c", "pT, Eta, Zvtx Efficiency distrubution", 400,900);
    c->Divide(1,2);
    TH3D* hnum; TH3D* hden;

    if(specie=="All"||specie=="pi"||specie=="K"||specie=="p"||specie=="e"||specie=="mu"||specie=="SigmaP"||specie=="SigmaM") {

      TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
      AliCFContainer *data = (AliCFContainer*) (d->Get(Form("containerpp13TeV_%s",specie.Data())));
    
      AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
      THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
      AliCFGridSparse* gridSparseden = (AliCFGridSparse*)data->GetGrid(6); // Reco
      THnSparse* denData = (THnSparse*)gridSparseden->GetGrid();    
    
      hnum = (TH3D*)numData->Projection(0,1,4);
      hden = (TH3D*)denData->Projection(0,1,4);

    } else if(specie=="5Species") {  //*** Sum of p,K,pi,e,mu ***

      //pi
      TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");

      AliCFContainer *data1 = (AliCFContainer*) (d->Get("containerpp13TeV_pi"));
    
      AliCFGridSparse* gridSparsenum1 = (AliCFGridSparse*)data1->GetGrid(1); // GenAcc
      THnSparse* numData1 = (THnSparse*)gridSparsenum1->GetGrid();
      AliCFGridSparse* gridSparseden1 = (AliCFGridSparse*)data1->GetGrid(6); // Reco
      THnSparse* denData1 = (THnSparse*)gridSparseden1->GetGrid();

      //K
      AliCFContainer *data2 = (AliCFContainer*) (d->Get("containerpp13TeV_K"));
    
      AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
      THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();
      AliCFGridSparse* gridSparseden2 = (AliCFGridSparse*)data2->GetGrid(6); // Reco
      THnSparse* denData2 = (THnSparse*)gridSparseden2->GetGrid();
    
      //p
      AliCFContainer *data3 = (AliCFContainer*) (d->Get("containerpp13TeV_p"));
    
      AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
      THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();
      AliCFGridSparse* gridSparseden3 = (AliCFGridSparse*)data3->GetGrid(6); // Reco
      THnSparse* denData3 = (THnSparse*)gridSparseden3->GetGrid();
    
      //e
      AliCFContainer *data4 = (AliCFContainer*) (d->Get("containerpp13TeV_e"));
    
      AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
      THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();
      AliCFGridSparse* gridSparseden4 = (AliCFGridSparse*)data4->GetGrid(6); // Reco
      THnSparse* denData4 = (THnSparse*)gridSparseden4->GetGrid();
    
      //mu
      AliCFContainer *data5 = (AliCFContainer*) (d->Get("containerpp13TeV_mu"));
    
      AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
      THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();
      AliCFGridSparse* gridSparseden5 = (AliCFGridSparse*)data5->GetGrid(6); // Reco
      THnSparse* denData5 = (THnSparse*)gridSparseden5->GetGrid();
             
      TH3D* hnum1 = (TH3D*)numData1->Projection(0,1,4);
      TH3D* hden1 = (TH3D*)denData1->Projection(0,1,4);
    
      TH3D* hnum2 = (TH3D*)numData2->Projection(0,1,4);
      TH3D* hden2 = (TH3D*)denData2->Projection(0,1,4);
    
      TH3D* hnum3 = (TH3D*)numData3->Projection(0,1,4);
      TH3D* hden3 = (TH3D*)denData3->Projection(0,1,4);
    
      TH3D* hnum4 = (TH3D*)numData4->Projection(0,1,4);
      TH3D* hden4 = (TH3D*)denData4->Projection(0,1,4);
    
      TH3D* hnum5 = (TH3D*)numData5->Projection(0,1,4);
      TH3D* hden5 = (TH3D*)denData5->Projection(0,1,4);    

      hnum = (TH3D*)hnum1->Clone("hnumtot");
      hnum->Add(hnum2);
      hnum->Add(hnum3);
      hnum->Add(hnum4);
      hnum->Add(hnum5);
      hden = (TH3D*)hden1->Clone("hdentot");
      hden->Add(hden2);
      hden->Add(hden3);
      hden->Add(hden4);
      hden->Add(hden5);

    } else if(specie=="Rest") {  //***"Rest" means ALL - p,K,pi,SigmaP,SigmaM, but includes electrons and muons! Is used for LF-style reweightings, not for our 5-specie procedure!!

      //All
      TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");

      AliCFContainer *dataAll = (AliCFContainer*) (d->Get("containerpp13TeV_All"));
    
      AliCFGridSparse* gridSparsenumAll = (AliCFGridSparse*)dataAll->GetGrid(1); // GenAcc
      THnSparse* numDataAll = (THnSparse*)gridSparsenumAll->GetGrid();
      AliCFGridSparse* gridSparsedenAll = (AliCFGridSparse*)dataAll->GetGrid(6); // Reco
      THnSparse* denDataAll = (THnSparse*)gridSparsedenAll->GetGrid();

      //pi
      AliCFContainer *data1 = (AliCFContainer*) (d->Get("containerpp13TeV_pi"));
    
      AliCFGridSparse* gridSparsenum1 = (AliCFGridSparse*)data1->GetGrid(1); // GenAcc
      THnSparse* numData1 = (THnSparse*)gridSparsenum1->GetGrid();
      AliCFGridSparse* gridSparseden1 = (AliCFGridSparse*)data1->GetGrid(6); // Reco
      THnSparse* denData1 = (THnSparse*)gridSparseden1->GetGrid();

      //K
      AliCFContainer *data2 = (AliCFContainer*) (d->Get("containerpp13TeV_K"));
    
      AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
      THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();
      AliCFGridSparse* gridSparseden2 = (AliCFGridSparse*)data2->GetGrid(6); // Reco
      THnSparse* denData2 = (THnSparse*)gridSparseden2->GetGrid();
    
      //p
      AliCFContainer *data3 = (AliCFContainer*) (d->Get("containerpp13TeV_p"));
    
      AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
      THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();
      AliCFGridSparse* gridSparseden3 = (AliCFGridSparse*)data3->GetGrid(6); // Reco
      THnSparse* denData3 = (THnSparse*)gridSparseden3->GetGrid();
    
      //Sigma+
      AliCFContainer *data4 = (AliCFContainer*) (d->Get("containerpp13TeV_SigmaP"));
    
      AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
      THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();
      AliCFGridSparse* gridSparseden4 = (AliCFGridSparse*)data4->GetGrid(6); // Reco
      THnSparse* denData4 = (THnSparse*)gridSparseden4->GetGrid();
    
      //sigma-
      AliCFContainer *data5 = (AliCFContainer*) (d->Get("containerpp13TeV_SigmaM"));
    
      AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
      THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();
      AliCFGridSparse* gridSparseden5 = (AliCFGridSparse*)data5->GetGrid(6); // Reco
      THnSparse* denData5 = (THnSparse*)gridSparseden5->GetGrid();
    
      TH3D* hnumAll = (TH3D*)numDataAll->Projection(0,1,4);
      TH3D* hdenAll = (TH3D*)denDataAll->Projection(0,1,4);
           
      TH3D* hnum1 = (TH3D*)numData1->Projection(0,1,4);
      TH3D* hden1 = (TH3D*)denData1->Projection(0,1,4);
    
      TH3D* hnum2 = (TH3D*)numData2->Projection(0,1,4);
      TH3D* hden2 = (TH3D*)denData2->Projection(0,1,4);
    
      TH3D* hnum3 = (TH3D*)numData3->Projection(0,1,4);
      TH3D* hden3 = (TH3D*)denData3->Projection(0,1,4);
    
      TH3D* hnum4 = (TH3D*)numData4->Projection(0,1,4);
      TH3D* hden4 = (TH3D*)denData4->Projection(0,1,4);
    
      TH3D* hnum5 = (TH3D*)numData5->Projection(0,1,4);
      TH3D* hden5 = (TH3D*)denData5->Projection(0,1,4);    

      hnum = (TH3D*)hnumAll->Clone("hnumtot");
      hnum->Add(hnum1,-1);
      hnum->Add(hnum2,-1);
      hnum->Add(hnum3,-1);
      hnum->Add(hnum4,-1);
      hnum->Add(hnum5,-1);
      hden = (TH3D*)hdenAll->Clone("hdentot");
      hden->Add(hden1,-1);
      hden->Add(hden2,-1);
      hden->Add(hden3,-1);
      hden->Add(hden4,-1);
      hden->Add(hden5,-1);

    } else {
    	printf("Error! Specie not supported!\n");
    	return;
    }

    TH3D* heff = (TH3D*)hden->Clone("heff");
    heff->Divide(hden,hnum,1,1,"B"); //ok, names of den and num are inverted, because they are inverted also in the definition :/ It works fine...
    
    c->cd(1);
    heff->Draw();
    
    //Now we do the rebinned efficiency
    TH3D* Rbhnum; TH3D* Rbhden;

    if(specie=="All"||specie=="pi"||specie=="K"||specie=="p"||specie=="e"||specie=="mu"||specie=="SigmaP"||specie=="SigmaM") {   

      TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
      AliCFContainer *data = (AliCFContainer*) (d->Get(Form("containerpp13TeV_%s",specie.Data())));

      AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
      THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
      THnSparse* RbnumData = (THnSparse*)numData->Clone("numNew");
      RbnumData->Reset();
      TAxis* axis = (TAxis*)RbnumData->GetAxis(0);
      axis->Set(newLimitsBins,newLimits);
      RbnumData->SetBinEdges(0,newLimits);
      RbnumData->RebinnedAdd(numData, 1);
      TH3D* Rbhnum = (TH3D*)RbnumData->Projection(0,1,4);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparseden = (AliCFGridSparse*)data->GetGrid(6); // Reco
      THnSparse* denData = (THnSparse*)gridSparseden->GetGrid();       
      THnSparse* RbdenData = (THnSparse*)denData->Clone("denNew");
      RbdenData->Reset();
      TAxis* axis2 = (TAxis*)RbdenData->GetAxis(0);
      axis2->Set(newLimitsBins,newLimits);
      RbdenData->SetBinEdges(0,newLimits);
      RbdenData->RebinnedAdd(denData, 1);
      TH3D* Rbhden = (TH3D*)RbdenData->Projection(0,1,4);
      gStyle->SetOptStat(kTRUE);

    } else if(specie=="5Species") {  //*** Sum of pi,K,p,e,mu ***

	//pi
	TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");

	AliCFContainer *data1 = (AliCFContainer*) (d->Get("containerpp13TeV_pi"));

	AliCFGridSparse* gridSparsenum1 = (AliCFGridSparse*)data1->GetGrid(1); // GenAcc
	THnSparse* numData1 = (THnSparse*)gridSparsenum1->GetGrid();	
	THnSparse* RbnumData1 = (THnSparse*)numData1->Clone("numNew1");
	RbnumData1->Reset();
	TAxis* axis3 = (TAxis*)RbnumData1->GetAxis(0);
	axis3->Set(newLimitsBins,newLimits);
	RbnumData1->SetBinEdges(0,newLimits);
	RbnumData1->RebinnedAdd(numData1, 1);
	TH3D* Rbhnum1 = (TH3D*)RbnumData1->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	AliCFGridSparse* gridSparseden1 = (AliCFGridSparse*)data1->GetGrid(6); // Reco
	THnSparse* denData1 = (THnSparse*)gridSparseden1->GetGrid();   		
	THnSparse* RbdenData1 = (THnSparse*)denData1->Clone("denNew1");
	RbdenData1->Reset();
	TAxis* axis4 = (TAxis*)RbdenData1->GetAxis(0);
	axis4->Set(newLimitsBins,newLimits);
	RbdenData1->SetBinEdges(0,newLimits);
	RbdenData1->RebinnedAdd(denData1, 1);
	TH3D* Rbhden1 = (TH3D*)RbdenData1->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	//K
	AliCFContainer *data2 = (AliCFContainer*) (d->Get("containerpp13TeV_K"));

	AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
	THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();	
	THnSparse* RbnumData2 = (THnSparse*)numData2->Clone("numNew2");
	RbnumData2->Reset();
	TAxis* axis5 = (TAxis*)RbnumData2->GetAxis(0);
	axis5->Set(newLimitsBins,newLimits);
	RbnumData2->SetBinEdges(0,newLimits);
	RbnumData2->RebinnedAdd(numData2, 1);
	TH3D* Rbhnum2 = (TH3D*)RbnumData2->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	AliCFGridSparse* gridSparseden2 = (AliCFGridSparse*)data2->GetGrid(6); // Reco
	THnSparse* denData2 = (THnSparse*)gridSparseden2->GetGrid();   	
	THnSparse* RbdenData2 = (THnSparse*)denData2->Clone("denNew2");
	RbdenData2->Reset();
	TAxis* axis6 = (TAxis*)RbdenData2->GetAxis(0);
	axis6->Set(newLimitsBins,newLimits);
	RbdenData2->SetBinEdges(0,newLimits);
	RbdenData2->RebinnedAdd(denData2, 1);
	TH3D* Rbhden2 = (TH3D*)RbdenData2->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);
	
	//p
	AliCFContainer *data3 = (AliCFContainer*) (d->Get("containerpp13TeV_p"));

	AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
	THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();		
	THnSparse* RbnumData3 = (THnSparse*)numData3->Clone("numNew3");
	RbnumData3->Reset();
	TAxis* axis7 = (TAxis*)RbnumData3->GetAxis(0);
	axis7->Set(newLimitsBins,newLimits);
	RbnumData3->SetBinEdges(0,newLimits);
	RbnumData3->RebinnedAdd(numData3, 1);
	TH3D* Rbhnum3 = (TH3D*)RbnumData3->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	AliCFGridSparse* gridSparseden3 = (AliCFGridSparse*)data3->GetGrid(6); // Reco
	THnSparse* denData3 = (THnSparse*)gridSparseden3->GetGrid();   		
	THnSparse* RbdenData3 = (THnSparse*)denData3->Clone("denNew3");
	RbdenData3->Reset();
	TAxis* axis8 = (TAxis*)RbdenData3->GetAxis(0);
	axis8->Set(newLimitsBins,newLimits);
	RbdenData3->SetBinEdges(0,newLimits);
	RbdenData3->RebinnedAdd(denData3, 1);
	TH3D* Rbhden3 = (TH3D*)RbdenData3->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	//e
	AliCFContainer *data4 = (AliCFContainer*) (d->Get("containerpp13TeV_e"));

	AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
	THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();		
	THnSparse* RbnumData4 = (THnSparse*)numData4->Clone("numNew4");
	RbnumData4->Reset();
	TAxis* axis9 = (TAxis*)RbnumData4->GetAxis(0);
	axis9->Set(newLimitsBins,newLimits);
	RbnumData4->SetBinEdges(0,newLimits);
	RbnumData4->RebinnedAdd(numData4, 1);
	TH3D* Rbhnum4 = (TH3D*)RbnumData4->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	AliCFGridSparse* gridSparseden4 = (AliCFGridSparse*)data4->GetGrid(6); // Reco
	THnSparse* denData4 = (THnSparse*)gridSparseden4->GetGrid();   		
	THnSparse* RbdenData4 = (THnSparse*)denData4->Clone("denNew4");
	RbdenData4->Reset();
	TAxis* axis10 = (TAxis*)RbdenData4->GetAxis(0);
	axis10->Set(newLimitsBins,newLimits);
	RbdenData4->SetBinEdges(0,newLimits);
	RbdenData4->RebinnedAdd(denData4, 1);
	TH3D* Rbhden4 = (TH3D*)RbdenData4->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	//mu
	AliCFContainer *data5 = (AliCFContainer*) (d->Get("containerpp13TeV_mu"));

	AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
	THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();		
	THnSparse* RbnumData5 = (THnSparse*)numData5->Clone("numNew5");
	RbnumData5->Reset();
	TAxis* axis11 = (TAxis*)RbnumData5->GetAxis(0);
	axis11->Set(newLimitsBins,newLimits);
	RbnumData5->SetBinEdges(0,newLimits);
	RbnumData5->RebinnedAdd(numData5, 1);
	TH3D* Rbhnum5 = (TH3D*)RbnumData5->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	AliCFGridSparse* gridSparseden5 = (AliCFGridSparse*)data5->GetGrid(6); // Reco
	THnSparse* denData5 = (THnSparse*)gridSparseden5->GetGrid();   	
	THnSparse* RbdenData5 = (THnSparse*)denData5->Clone("denNew5");
	RbdenData5->Reset();
	TAxis* axis12 = (TAxis*)RbdenData5->GetAxis(0);
	axis12->Set(newLimitsBins,newLimits);
	RbdenData5->SetBinEdges(0,newLimits);
	RbdenData5->RebinnedAdd(denData5, 1);
	TH3D* Rbhden5 = (TH3D*)RbdenData5->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);
    
	Rbhnum = (TH3D*)Rbhnum1->Clone("Rbhnumtot");
	Rbhnum->Add(Rbhnum2);
	Rbhnum->Add(Rbhnum3);
	Rbhnum->Add(Rbhnum4);
	Rbhnum->Add(Rbhnum5);
	Rbhden = (TH3D*)Rbhden1->Clone("Rbhdentot");
	Rbhden->Add(Rbhden2);
	Rbhden->Add(Rbhden3);
	Rbhden->Add(Rbhden4);
	Rbhden->Add(Rbhden5);

    } else if(specie=="Rest") {  //***"Rest" means ALL - p,K,pi,SigmaP,SigmaM, but includes electrons and muons! Is used for LF-style reweightings, not for our 5-specie procedure!!
                                 //one takes "All" num and den, and subtracts the identified species's nums and dens

      //All
	TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");

	AliCFContainer *dataAll = (AliCFContainer*) (d->Get("containerpp13TeV_All"));

	AliCFGridSparse* gridSparsenumAll = (AliCFGridSparse*)dataAll->GetGrid(1); // GenAcc
	THnSparse* numDataAll = (THnSparse*)gridSparsenumAll->GetGrid();
	THnSparse* RbnumDataAll = (THnSparse*)numDataAll->Clone("numNewAll");
	RbnumDataAll->Reset();
	TAxis* axis = (TAxis*)RbnumDataAll->GetAxis(0);
	axis->Set(newLimitsBins,newLimits);
	RbnumDataAll->SetBinEdges(0,newLimits);
	RbnumDataAll->RebinnedAdd(numDataAll, 1);
	TH3D* RbhnumAll = (TH3D*)RbnumDataAll->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	AliCFGridSparse* gridSparsedenAll = (AliCFGridSparse*)dataAll->GetGrid(6); // Reco
	THnSparse* denDataAll = (THnSparse*)gridSparsedenAll->GetGrid();   	
	THnSparse* RbdenDataAll = (THnSparse*)denDataAll->Clone("denNew");
	RbdenDataAll->Reset();
	TAxis* axis2 = (TAxis*)RbdenDataAll->GetAxis(0);
	axis2->Set(newLimitsBins,newLimits);
	RbdenDataAll->SetBinEdges(0,newLimits);
	RbdenDataAll->RebinnedAdd(denDataAll, 1);
	TH3D* RbhdenAll = (TH3D*)RbdenDataAll->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	//pi
	AliCFContainer *data1 = (AliCFContainer*) (d->Get("containerpp13TeV_pi"));

	AliCFGridSparse* gridSparsenum1 = (AliCFGridSparse*)data1->GetGrid(1); // GenAcc
	THnSparse* numData1 = (THnSparse*)gridSparsenum1->GetGrid();	
	THnSparse* RbnumData1 = (THnSparse*)numData1->Clone("numNew1");
	RbnumData1->Reset();
	TAxis* axis3 = (TAxis*)RbnumData1->GetAxis(0);
	axis3->Set(newLimitsBins,newLimits);
	RbnumData1->SetBinEdges(0,newLimits);
	RbnumData1->RebinnedAdd(numData1, 1);
	TH3D* Rbhnum1 = (TH3D*)RbnumData1->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	AliCFGridSparse* gridSparseden1 = (AliCFGridSparse*)data1->GetGrid(6); // Reco
	THnSparse* denData1 = (THnSparse*)gridSparseden1->GetGrid();   		
	THnSparse* RbdenData1 = (THnSparse*)denData1->Clone("denNew1");
	RbdenData1->Reset();
	TAxis* axis4 = (TAxis*)RbdenData1->GetAxis(0);
	axis4->Set(newLimitsBins,newLimits);
	RbdenData1->SetBinEdges(0,newLimits);
	RbdenData1->RebinnedAdd(denData1, 1);
	TH3D* Rbhden1 = (TH3D*)RbdenData1->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	//K
	AliCFContainer *data2 = (AliCFContainer*) (d->Get("containerpp13TeV_K"));

	AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
	THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();	
	THnSparse* RbnumData2 = (THnSparse*)numData2->Clone("numNew2");
	RbnumData2->Reset();
	TAxis* axis5 = (TAxis*)RbnumData2->GetAxis(0);
	axis5->Set(newLimitsBins,newLimits);
	RbnumData2->SetBinEdges(0,newLimits);
	RbnumData2->RebinnedAdd(numData2, 1);
	TH3D* Rbhnum2 = (TH3D*)RbnumData2->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	AliCFGridSparse* gridSparseden2 = (AliCFGridSparse*)data2->GetGrid(6); // Reco
	THnSparse* denData2 = (THnSparse*)gridSparseden2->GetGrid();   	
	THnSparse* RbdenData2 = (THnSparse*)denData2->Clone("denNew2");
	RbdenData2->Reset();
	TAxis* axis6 = (TAxis*)RbdenData2->GetAxis(0);
	axis6->Set(newLimitsBins,newLimits);
	RbdenData2->SetBinEdges(0,newLimits);
	RbdenData2->RebinnedAdd(denData2, 1);
	TH3D* Rbhden2 = (TH3D*)RbdenData2->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);
	
	//p
	AliCFContainer *data3 = (AliCFContainer*) (d->Get("containerpp13TeV_p"));

	AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
	THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();		
	THnSparse* RbnumData3 = (THnSparse*)numData3->Clone("numNew3");
	RbnumData3->Reset();
	TAxis* axis7 = (TAxis*)RbnumData3->GetAxis(0);
	axis7->Set(newLimitsBins,newLimits);
	RbnumData3->SetBinEdges(0,newLimits);
	RbnumData3->RebinnedAdd(numData3, 1);
	TH3D* Rbhnum3 = (TH3D*)RbnumData3->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	AliCFGridSparse* gridSparseden3 = (AliCFGridSparse*)data3->GetGrid(6); // Reco
	THnSparse* denData3 = (THnSparse*)gridSparseden3->GetGrid();   		
	THnSparse* RbdenData3 = (THnSparse*)denData3->Clone("denNew3");
	RbdenData3->Reset();
	TAxis* axis8 = (TAxis*)RbdenData3->GetAxis(0);
	axis8->Set(newLimitsBins,newLimits);
	RbdenData3->SetBinEdges(0,newLimits);
	RbdenData3->RebinnedAdd(denData3, 1);
	TH3D* Rbhden3 = (TH3D*)RbdenData3->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	//Sigma+
	AliCFContainer *data4 = (AliCFContainer*) (d->Get("containerpp13TeV_SigmaP"));

	AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
	THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();		
	THnSparse* RbnumData4 = (THnSparse*)numData4->Clone("numNew4");
	RbnumData4->Reset();
	TAxis* axis9 = (TAxis*)RbnumData4->GetAxis(0);
	axis9->Set(newLimitsBins,newLimits);
	RbnumData4->SetBinEdges(0,newLimits);
	RbnumData4->RebinnedAdd(numData4, 1);
	TH3D* Rbhnum4 = (TH3D*)RbnumData4->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	AliCFGridSparse* gridSparseden4 = (AliCFGridSparse*)data4->GetGrid(6); // Reco
	THnSparse* denData4 = (THnSparse*)gridSparseden4->GetGrid();   		
	THnSparse* RbdenData4 = (THnSparse*)denData4->Clone("denNew4");
	RbdenData4->Reset();
	TAxis* axis10 = (TAxis*)RbdenData4->GetAxis(0);
	axis10->Set(newLimitsBins,newLimits);
	RbdenData4->SetBinEdges(0,newLimits);
	RbdenData4->RebinnedAdd(denData4, 1);
	TH3D* Rbhden4 = (TH3D*)RbdenData4->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	//Sigma-
	AliCFContainer *data5 = (AliCFContainer*) (d->Get("containerpp13TeV_SigmaM"));

	AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
	THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();		
	THnSparse* RbnumData5 = (THnSparse*)numData5->Clone("numNew5");
	RbnumData5->Reset();
	TAxis* axis11 = (TAxis*)RbnumData5->GetAxis(0);
	axis11->Set(newLimitsBins,newLimits);
	RbnumData5->SetBinEdges(0,newLimits);
	RbnumData5->RebinnedAdd(numData5, 1);
	TH3D* Rbhnum5 = (TH3D*)RbnumData5->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);

	AliCFGridSparse* gridSparseden5 = (AliCFGridSparse*)data5->GetGrid(6); // Reco
	THnSparse* denData5 = (THnSparse*)gridSparseden5->GetGrid();   	
	THnSparse* RbdenData5 = (THnSparse*)denData5->Clone("denNew5");
	RbdenData5->Reset();
	TAxis* axis12 = (TAxis*)RbdenData5->GetAxis(0);
	axis12->Set(newLimitsBins,newLimits);
	RbdenData5->SetBinEdges(0,newLimits);
	RbdenData5->RebinnedAdd(denData5, 1);
	TH3D* Rbhden5 = (TH3D*)RbdenData5->Projection(0,1,4);
	gStyle->SetOptStat(kTRUE);
    
	Rbhden = (TH3D*)RbhdenAll->Clone("Rbhdentot");
	Rbhden->Add(Rbhden1,-1);
	Rbhden->Add(Rbhden2,-1);
	Rbhden->Add(Rbhden3,-1);
	Rbhden->Add(Rbhden4,-1);
	Rbhden->Add(Rbhden5,-1);
	Rbhnum = (TH3D*)RbhnumAll->Clone("Rbhnumtot");
	Rbhnum->Add(Rbhnum1,-1);
	Rbhnum->Add(Rbhnum2,-1);
	Rbhnum->Add(Rbhnum3,-1);
	Rbhnum->Add(Rbhnum4,-1);
	Rbhnum->Add(Rbhnum5,-1);

    } else {
    	printf("Error! Specie not supported!\n");
    	return;
    }    
    
    TH3D* Rbheff = (TH3D*)Rbhden->Clone("heff_rebin");
    Rbheff->Divide(Rbhden,Rbhnum,1,1,"B");
    gStyle->SetOptStat(kTRUE);
    c->cd(2);
    Rbheff->Draw();
    
    // Num bins checking
    Bool_t Rebintest = kTRUE;
    Float_t sum = 0, sumnew = 0;
    Float_t content = 0, content2 = 0;
    if(Rebintest){
        Printf("Original THnSparse <-------------");
        for (Int_t ibin = 1; ibin<=hnum->GetNbinsX(); ibin++){
            content=0;
            for (Int_t ibin2 = 1; ibin2<=hnum->GetNbinsY(); ibin2++){
                for (Int_t ibin3 = 1; ibin3<=hnum->GetNbinsZ(); ibin3++){
                    content+=hnum->GetBinContent(ibin,ibin2,ibin3);
                }
            }
            Printf("ibin = %d, low edge = %0.2f, content = %0.2f",ibin,hnum->GetXaxis()->GetBinLowEdge(ibin), content);
            sum+=content;
            
        }
        Printf("Rebinned THnSparse <-------------");
        for (Int_t ibin = 1; ibin<=Rbhnum->GetNbinsX(); ibin++){
            content2=0;
            for (Int_t ibin2 = 1; ibin2<=Rbhnum->GetNbinsY(); ibin2++){
                for (Int_t ibin3 = 1; ibin3<=Rbhnum->GetNbinsZ(); ibin3++){
                    content2+=Rbhnum->GetBinContent(ibin,ibin2,ibin3);
                }
            }
            Printf("ibin = %d, low edge = %0.2f, content = %0.2f",ibin,Rbhnum->GetXaxis()->GetBinLowEdge(ibin), content2);
            sumnew+=content2;
        }
        Printf("Sum_before = %f, Sum_After = %f, Checking difference.. =%f%s", sum, sumnew, (sum-sumnew)*100/sum,"%");
    }
    
    
    if(((sum-sumnew)*100/sum)<.01){
        TFile* fout = new TFile(Form("3D_TrackingEffMap_%s.root",specie.Data()),"RECREATE");
        c->Write();
        fout->Close();
        printf( "1>>> rebin OKAY, efficiency map file created");
        c->SaveAs(Form("3D_TrackingEffMap_%s.png",specie.Data()));
        CheckFileQuality(specie);
    }
    else Printf(" ======> rebin Failed");
}


//function to check 3D file Quality (e.g. Empty bins or low/high efficiency)
void CheckFileQuality(TString specie){

    Printf( "\n2>>> Checking file Quality..");
    TFile fC(Form("3D_TrackingEffMap_%s.root",specie.Data()));
    TCanvas *cPass2 = (TCanvas*)fC.Get("c");
    TH3D *hcEffPass2=(TH3D*)cPass2->FindObject("heff_rebin");
    
    cout << "Number of Pt bin = " <<hcEffPass2->GetNbinsX() << endl;
    cout << "Number of Eta bin = " <<hcEffPass2->GetNbinsY() << endl;
    cout << "Number of Zvtx bin = " <<hcEffPass2->GetNbinsZ() << endl;

    Printf( "\n**Suspected bins but can be due to less stats in bin");
    for(int i=1; i<=hcEffPass2->GetNbinsX(); i++){
        for(int j=1; j<=hcEffPass2->GetNbinsY(); j++){
            for(int k=1; k<=hcEffPass2->GetNbinsZ(); k++){
                if(hcEffPass2->GetBinContent(i,j,k) < 0.50 || hcEffPass2->GetBinContent(i,j,k) > 1.00){
                    cout <<"(Xbin = " << i <<", Ybin = "<<j<< ", Zbin = "<< k <<") -----> " <<  hcEffPass2->GetBinContent(i,j,k) << endl;
                }
            }
        }
    }
}


void DrawDhCorr_STE_1DEfficiency_in_pT(TString filename="AnalysisResults.root", TString specie="All") { //Species: "All", "pi", "K", "p", "e", "mu", "5Species", "SigP", "SigM", "Rest"
    
    //New bins for final rebinning
    const Int_t newLimitsBins = 22;
    Double_t* newLimits = new Double_t[newLimitsBins+1];
    newLimits[0]   = 0.30;
    newLimits[1]   = 0.40;
    newLimits[2]   = 0.50;
    newLimits[3]   = 0.60;
    newLimits[4]   = 0.70;
    newLimits[5]   = 0.80;
    newLimits[6]   = 0.90;
    newLimits[7]   = 1.00;
    newLimits[8]   = 1.20;
    newLimits[9]   = 1.40;
    newLimits[10]  = 1.60;
    newLimits[11]  = 1.80;
    newLimits[12]  = 2.00;
    newLimits[13]  = 2.40;
    newLimits[14]  = 2.80;
    newLimits[15]  = 3.40;
    newLimits[16]  = 4.00;
    newLimits[17]  = 4.80;
    newLimits[18]  = 6.00;
    newLimits[19]  = 8.00;
    newLimits[20]  = 10.00;
    newLimits[21]  = 16.00;
    newLimits[22]  = 24.00;

    TFile* f = new TFile(filename.Data(),"read");
    TCanvas *c = new TCanvas("c", "pT, Eta, Zvtx Efficiency distrubution", 400,900);
    c->Divide(1,2);
    TH1D* hnum; TH1D* hden;

    if(specie=="All"||specie=="pi"||specie=="K"||specie=="p"||specie=="e"||specie=="mu"||specie=="SigmaP"||specie=="SigmaM") {

      TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
      AliCFContainer *data = (AliCFContainer*) (d->Get(Form("containerpp13TeV_%s",specie.Data())));
    
      AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
      THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
      AliCFGridSparse* gridSparseden = (AliCFGridSparse*)data->GetGrid(6); // Reco
      THnSparse* denData = (THnSparse*)gridSparseden->GetGrid();    
    
      hnum = (TH1D*)numData->Projection(0);
      hden = (TH1D*)denData->Projection(0);

    } else if(specie=="5Species") {  //***"Rest" means ALL - p,K,pi,SigmaP,SigmaM, but includes electrons and muons! Is used for LF-style reweightings, not for our 5-specie procedure!!

      //pi
      TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");

      AliCFContainer *data1 = (AliCFContainer*) (d->Get("containerpp13TeV_pi"));
    
      AliCFGridSparse* gridSparsenum1 = (AliCFGridSparse*)data1->GetGrid(1); // GenAcc
      THnSparse* numData1 = (THnSparse*)gridSparsenum1->GetGrid();
      AliCFGridSparse* gridSparseden1 = (AliCFGridSparse*)data1->GetGrid(6); // Reco
      THnSparse* denData1 = (THnSparse*)gridSparseden1->GetGrid();

      //K
      AliCFContainer *data2 = (AliCFContainer*) (d->Get("containerpp13TeV_K"));
    
      AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
      THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();
      AliCFGridSparse* gridSparseden2 = (AliCFGridSparse*)data2->GetGrid(6); // Reco
      THnSparse* denData2 = (THnSparse*)gridSparseden2->GetGrid();
    
      //p
      AliCFContainer *data3 = (AliCFContainer*) (d->Get("containerpp13TeV_p"));
    
      AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
      THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();
      AliCFGridSparse* gridSparseden3 = (AliCFGridSparse*)data3->GetGrid(6); // Reco
      THnSparse* denData3 = (THnSparse*)gridSparseden3->GetGrid();
    
      //e
      AliCFContainer *data4 = (AliCFContainer*) (d->Get("containerpp13TeV_e"));
    
      AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
      THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();
      AliCFGridSparse* gridSparseden4 = (AliCFGridSparse*)data4->GetGrid(6); // Reco
      THnSparse* denData4 = (THnSparse*)gridSparseden4->GetGrid();
    
      //mu
      AliCFContainer *data5 = (AliCFContainer*) (d->Get("containerpp13TeV_mu"));
    
      AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
      THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();
      AliCFGridSparse* gridSparseden5 = (AliCFGridSparse*)data5->GetGrid(6); // Reco
      THnSparse* denData5 = (THnSparse*)gridSparseden5->GetGrid();
           
      TH1D* hnum1 = (TH1D*)numData1->Projection(0);
      TH1D* hden1 = (TH1D*)denData1->Projection(0);
    
      TH1D* hnum2 = (TH1D*)numData2->Projection(0);
      TH1D* hden2 = (TH1D*)denData2->Projection(0);
    
      TH1D* hnum3 = (TH1D*)numData3->Projection(0);
      TH1D* hden3 = (TH1D*)denData3->Projection(0);
    
      TH1D* hnum4 = (TH1D*)numData4->Projection(0);
      TH1D* hden4 = (TH1D*)denData4->Projection(0);
    
      TH1D* hnum5 = (TH1D*)numData5->Projection(0);
      TH1D* hden5 = (TH1D*)denData5->Projection(0);    

      hnum = (TH1D*)hnum1->Clone("hnumtot");
      hnum->Add(hnum2);
      hnum->Add(hnum3);
      hnum->Add(hnum4);
      hnum->Add(hnum5);
      hden = (TH1D*)hden1->Clone("hdentot");
      hden->Add(hden2);
      hden->Add(hden3);
      hden->Add(hden4);
      hden->Add(hden5);

    } else if(specie=="Rest") {  //***"Rest" means ALL - p,K,pi,SigmaP,SigmaM, but includes electrons and muons! Is used for LF-style reweightings, not for our 5-specie procedure!!

      //All
      TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");

      AliCFContainer *dataAll = (AliCFContainer*) (d->Get("containerpp13TeV_All"));
    
      AliCFGridSparse* gridSparsenumAll = (AliCFGridSparse*)dataAll->GetGrid(1); // GenAcc
      THnSparse* numDataAll = (THnSparse*)gridSparsenumAll->GetGrid();
      AliCFGridSparse* gridSparsedenAll = (AliCFGridSparse*)dataAll->GetGrid(6); // Reco
      THnSparse* denDataAll = (THnSparse*)gridSparsedenAll->GetGrid();

      //pi
      AliCFContainer *data1 = (AliCFContainer*) (d->Get("containerpp13TeV_pi"));
    
      AliCFGridSparse* gridSparsenum1 = (AliCFGridSparse*)data1->GetGrid(1); // GenAcc
      THnSparse* numData1 = (THnSparse*)gridSparsenum1->GetGrid();
      AliCFGridSparse* gridSparseden1 = (AliCFGridSparse*)data1->GetGrid(6); // Reco
      THnSparse* denData1 = (THnSparse*)gridSparseden1->GetGrid();

      //K
      AliCFContainer *data2 = (AliCFContainer*) (d->Get("containerpp13TeV_K"));
    
      AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
      THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();
      AliCFGridSparse* gridSparseden2 = (AliCFGridSparse*)data2->GetGrid(6); // Reco
      THnSparse* denData2 = (THnSparse*)gridSparseden2->GetGrid();
    
      //p
      AliCFContainer *data3 = (AliCFContainer*) (d->Get("containerpp13TeV_p"));
    
      AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
      THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();
      AliCFGridSparse* gridSparseden3 = (AliCFGridSparse*)data3->GetGrid(6); // Reco
      THnSparse* denData3 = (THnSparse*)gridSparseden3->GetGrid();
    
      //Sigma+
      AliCFContainer *data4 = (AliCFContainer*) (d->Get("containerpp13TeV_SigmaP"));
    
      AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
      THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();
      AliCFGridSparse* gridSparseden4 = (AliCFGridSparse*)data4->GetGrid(6); // Reco
      THnSparse* denData4 = (THnSparse*)gridSparseden4->GetGrid();
    
      //sigma-
      AliCFContainer *data5 = (AliCFContainer*) (d->Get("containerpp13TeV_SigmaM"));
    
      AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
      THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();
      AliCFGridSparse* gridSparseden5 = (AliCFGridSparse*)data5->GetGrid(6); // Reco
      THnSparse* denData5 = (THnSparse*)gridSparseden5->GetGrid();
    
      TH1D* hnumAll = (TH1D*)numDataAll->Projection(0);
      TH1D* hdenAll = (TH1D*)denDataAll->Projection(0);
           
      TH1D* hnum1 = (TH1D*)numData1->Projection(0);
      TH1D* hden1 = (TH1D*)denData1->Projection(0);
    
      TH1D* hnum2 = (TH1D*)numData2->Projection(0);
      TH1D* hden2 = (TH1D*)denData2->Projection(0);
    
      TH1D* hnum3 = (TH1D*)numData3->Projection(0);
      TH1D* hden3 = (TH1D*)denData3->Projection(0);
    
      TH1D* hnum4 = (TH1D*)numData4->Projection(0);
      TH1D* hden4 = (TH1D*)denData4->Projection(0);
    
      TH1D* hnum5 = (TH1D*)numData5->Projection(0);
      TH1D* hden5 = (TH1D*)denData5->Projection(0);    

      hnum = (TH1D*)hnumAll->Clone("hnumtot");
      hnum->Add(hnum1,-1);
      hnum->Add(hnum2,-1);
      hnum->Add(hnum3,-1);
      hnum->Add(hnum4,-1);
      hnum->Add(hnum5,-1);
      hden = (TH1D*)hdenAll->Clone("hdentot");
      hden->Add(hden1,-1);
      hden->Add(hden2,-1);
      hden->Add(hden3,-1);
      hden->Add(hden4,-1);
      hden->Add(hden5,-1);

    } else {
      printf("Error! Specie not supported!\n");
      return;
    }

    TH1D* heff = (TH1D*)hden->Clone("heff");
    heff->Divide(hden,hnum,1,1,"B"); //ok, names of den and num are inverted, because they are inverted also in the definition :/ It works fine...
    
    c->cd(1);
    heff->SetMarkerStyle(21);
    heff->SetMarkerSize(0.8);
    heff->Draw();
    
    //Now we do the rebinned efficiency
    TH1D* Rbhnum; TH1D* Rbhden;

    if(specie=="All"||specie=="pi"||specie=="K"||specie=="p"||specie=="e"||specie=="mu"||specie=="SigmaP"||specie=="SigmaM") {   

      TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
      AliCFContainer *data = (AliCFContainer*) (d->Get(Form("containerpp13TeV_%s",specie.Data())));

      AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
      THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
      THnSparse* RbnumData = (THnSparse*)numData->Clone("numNew");
      RbnumData->Reset();
      TAxis* axis = (TAxis*)RbnumData->GetAxis(0);
      axis->Set(newLimitsBins,newLimits);
      RbnumData->SetBinEdges(0,newLimits);
      RbnumData->RebinnedAdd(numData, 1);
      TH1D* Rbhnum = (TH1D*)RbnumData->Projection(0);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparseden = (AliCFGridSparse*)data->GetGrid(6); // Reco
      THnSparse* denData = (THnSparse*)gridSparseden->GetGrid();       
      THnSparse* RbdenData = (THnSparse*)denData->Clone("denNew");
      RbdenData->Reset();
      TAxis* axis2 = (TAxis*)RbdenData->GetAxis(0);
      axis2->Set(newLimitsBins,newLimits);
      RbdenData->SetBinEdges(0,newLimits);
      RbdenData->RebinnedAdd(denData, 1);
      TH1D* Rbhden = (TH1D*)RbdenData->Projection(0);
      gStyle->SetOptStat(kTRUE);

	} else if(specie=="5Species") {  //Sum of pi,K,p,e,mu

      //pi
      TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");

      AliCFContainer *data1 = (AliCFContainer*) (d->Get("containerpp13TeV_pi"));

      AliCFGridSparse* gridSparsenum1 = (AliCFGridSparse*)data1->GetGrid(1); // GenAcc
      THnSparse* numData1 = (THnSparse*)gridSparsenum1->GetGrid();      
      THnSparse* RbnumData1 = (THnSparse*)numData1->Clone("numNew1");
      RbnumData1->Reset();
      TAxis* axis3 = (TAxis*)RbnumData1->GetAxis(0);
      axis3->Set(newLimitsBins,newLimits);
      RbnumData1->SetBinEdges(0,newLimits);
      RbnumData1->RebinnedAdd(numData1, 1);
      TH1D* Rbhnum1 = (TH1D*)RbnumData1->Projection(0);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparseden1 = (AliCFGridSparse*)data1->GetGrid(6); // Reco
      THnSparse* denData1 = (THnSparse*)gridSparseden1->GetGrid();              
      THnSparse* RbdenData1 = (THnSparse*)denData1->Clone("denNew1");
      RbdenData1->Reset();
      TAxis* axis4 = (TAxis*)RbdenData1->GetAxis(0);
      axis4->Set(newLimitsBins,newLimits);
      RbdenData1->SetBinEdges(0,newLimits);
      RbdenData1->RebinnedAdd(denData1, 1);
      TH1D* Rbhden1 = (TH1D*)RbdenData1->Projection(0);
      gStyle->SetOptStat(kTRUE);

      //K
      AliCFContainer *data2 = (AliCFContainer*) (d->Get("containerpp13TeV_K"));

      AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
      THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();  
      THnSparse* RbnumData2 = (THnSparse*)numData2->Clone("numNew2");
      RbnumData2->Reset();
      TAxis* axis5 = (TAxis*)RbnumData2->GetAxis(0);
      axis5->Set(newLimitsBins,newLimits);
      RbnumData2->SetBinEdges(0,newLimits);
      RbnumData2->RebinnedAdd(numData2, 1);
      TH1D* Rbhnum2 = (TH1D*)RbnumData2->Projection(0);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparseden2 = (AliCFGridSparse*)data2->GetGrid(6); // Reco
      THnSparse* denData2 = (THnSparse*)gridSparseden2->GetGrid();        
      THnSparse* RbdenData2 = (THnSparse*)denData2->Clone("denNew2");
      RbdenData2->Reset();
      TAxis* axis6 = (TAxis*)RbdenData2->GetAxis(0);
      axis6->Set(newLimitsBins,newLimits);
      RbdenData2->SetBinEdges(0,newLimits);
      RbdenData2->RebinnedAdd(denData2, 1);
      TH1D* Rbhden2 = (TH1D*)RbdenData2->Projection(0);
      gStyle->SetOptStat(kTRUE);
      
      //p
      AliCFContainer *data3 = (AliCFContainer*) (d->Get("containerpp13TeV_p"));

      AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
      THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();        
      THnSparse* RbnumData3 = (THnSparse*)numData3->Clone("numNew3");
      RbnumData3->Reset();
      TAxis* axis7 = (TAxis*)RbnumData3->GetAxis(0);
      axis7->Set(newLimitsBins,newLimits);
      RbnumData3->SetBinEdges(0,newLimits);
      RbnumData3->RebinnedAdd(numData3, 1);
      TH1D* Rbhnum3 = (TH1D*)RbnumData3->Projection(0);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparseden3 = (AliCFGridSparse*)data3->GetGrid(6); // Reco
      THnSparse* denData3 = (THnSparse*)gridSparseden3->GetGrid();          
      THnSparse* RbdenData3 = (THnSparse*)denData3->Clone("denNew3");
      RbdenData3->Reset();
      TAxis* axis8 = (TAxis*)RbdenData3->GetAxis(0);
      axis8->Set(newLimitsBins,newLimits);
      RbdenData3->SetBinEdges(0,newLimits);
      RbdenData3->RebinnedAdd(denData3, 1);
      TH1D* Rbhden3 = (TH1D*)RbdenData3->Projection(0);
      gStyle->SetOptStat(kTRUE);

      //e
      AliCFContainer *data4 = (AliCFContainer*) (d->Get("containerpp13TeV_e"));

      AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
      THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();        
      THnSparse* RbnumData4 = (THnSparse*)numData4->Clone("numNew4");
      RbnumData4->Reset();
      TAxis* axis9 = (TAxis*)RbnumData4->GetAxis(0);
      axis9->Set(newLimitsBins,newLimits);
      RbnumData4->SetBinEdges(0,newLimits);
      RbnumData4->RebinnedAdd(numData4, 1);
      TH1D* Rbhnum4 = (TH1D*)RbnumData4->Projection(0);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparseden4 = (AliCFGridSparse*)data4->GetGrid(6); // Reco
      THnSparse* denData4 = (THnSparse*)gridSparseden4->GetGrid();          
      THnSparse* RbdenData4 = (THnSparse*)denData4->Clone("denNew4");
      RbdenData4->Reset();
      TAxis* axis10 = (TAxis*)RbdenData4->GetAxis(0);
      axis10->Set(newLimitsBins,newLimits);
      RbdenData4->SetBinEdges(0,newLimits);
      RbdenData4->RebinnedAdd(denData4, 1);
      TH1D* Rbhden4 = (TH1D*)RbdenData4->Projection(0);
      gStyle->SetOptStat(kTRUE);

      //mu
      AliCFContainer *data5 = (AliCFContainer*) (d->Get("containerpp13TeV_mu"));

      AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
      THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();        
      THnSparse* RbnumData5 = (THnSparse*)numData5->Clone("numNew5");
      RbnumData5->Reset();
      TAxis* axis11 = (TAxis*)RbnumData5->GetAxis(0);
      axis11->Set(newLimitsBins,newLimits);
      RbnumData5->SetBinEdges(0,newLimits);
      RbnumData5->RebinnedAdd(numData5, 1);
      TH1D* Rbhnum5 = (TH1D*)RbnumData5->Projection(0);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparseden5 = (AliCFGridSparse*)data5->GetGrid(6); // Reco
      THnSparse* denData5 = (THnSparse*)gridSparseden5->GetGrid();        
      THnSparse* RbdenData5 = (THnSparse*)denData5->Clone("denNew5");
      RbdenData5->Reset();
      TAxis* axis12 = (TAxis*)RbdenData5->GetAxis(0);
      axis12->Set(newLimitsBins,newLimits);
      RbdenData5->SetBinEdges(0,newLimits);
      RbdenData5->RebinnedAdd(denData5, 1);
      TH1D* Rbhden5 = (TH1D*)RbdenData5->Projection(0);
      gStyle->SetOptStat(kTRUE);
    
      Rbhden = (TH1D*)Rbhden1->Clone("Rbhdentot");
      Rbhden->Add(Rbhden2);
      Rbhden->Add(Rbhden3);
      Rbhden->Add(Rbhden4);
      Rbhden->Add(Rbhden5);
      Rbhnum = (TH1D*)Rbhnum1->Clone("Rbhnumtot");
      Rbhnum->Add(Rbhnum2);
      Rbhnum->Add(Rbhnum3);
      Rbhnum->Add(Rbhnum4);
      Rbhnum->Add(Rbhnum5);

    } else if(specie=="Rest") {  //***"Rest" means ALL - p,K,pi,SigmaP,SigmaM, but includes electrons and muons! Is used for LF-style reweightings, not for our 5-specie procedure!!
                                 //one takes "All" num and den, and subtracts the identified species's nums and dens

      //All
      TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");

      AliCFContainer *dataAll = (AliCFContainer*) (d->Get("containerpp13TeV_All"));

      AliCFGridSparse* gridSparsenumAll = (AliCFGridSparse*)dataAll->GetGrid(1); // GenAcc
      THnSparse* numDataAll = (THnSparse*)gridSparsenumAll->GetGrid();
      THnSparse* RbnumDataAll = (THnSparse*)numDataAll->Clone("numNewAll");
      RbnumDataAll->Reset();
      TAxis* axis = (TAxis*)RbnumDataAll->GetAxis(0);
      axis->Set(newLimitsBins,newLimits);
      RbnumDataAll->SetBinEdges(0,newLimits);
      RbnumDataAll->RebinnedAdd(numDataAll, 1);
      TH1D* RbhnumAll = (TH1D*)RbnumDataAll->Projection(0);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparsedenAll = (AliCFGridSparse*)dataAll->GetGrid(6); // Reco
      THnSparse* denDataAll = (THnSparse*)gridSparsedenAll->GetGrid();        
      THnSparse* RbdenDataAll = (THnSparse*)denDataAll->Clone("denNew");
      RbdenDataAll->Reset();
      TAxis* axis2 = (TAxis*)RbdenDataAll->GetAxis(0);
      axis2->Set(newLimitsBins,newLimits);
      RbdenDataAll->SetBinEdges(0,newLimits);
      RbdenDataAll->RebinnedAdd(denDataAll, 1);
      TH1D* RbhdenAll = (TH1D*)RbdenDataAll->Projection(0);
      gStyle->SetOptStat(kTRUE);

      //pi
      AliCFContainer *data1 = (AliCFContainer*) (d->Get("containerpp13TeV_pi"));

      AliCFGridSparse* gridSparsenum1 = (AliCFGridSparse*)data1->GetGrid(1); // GenAcc
      THnSparse* numData1 = (THnSparse*)gridSparsenum1->GetGrid();      
      THnSparse* RbnumData1 = (THnSparse*)numData1->Clone("numNew1");
      RbnumData1->Reset();
      TAxis* axis3 = (TAxis*)RbnumData1->GetAxis(0);
      axis3->Set(newLimitsBins,newLimits);
      RbnumData1->SetBinEdges(0,newLimits);
      RbnumData1->RebinnedAdd(numData1, 1);
      TH1D* Rbhnum1 = (TH1D*)RbnumData1->Projection(0);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparseden1 = (AliCFGridSparse*)data1->GetGrid(6); // Reco
      THnSparse* denData1 = (THnSparse*)gridSparseden1->GetGrid();              
      THnSparse* RbdenData1 = (THnSparse*)denData1->Clone("denNew1");
      RbdenData1->Reset();
      TAxis* axis4 = (TAxis*)RbdenData1->GetAxis(0);
      axis4->Set(newLimitsBins,newLimits);
      RbdenData1->SetBinEdges(0,newLimits);
      RbdenData1->RebinnedAdd(denData1, 1);
      TH1D* Rbhden1 = (TH1D*)RbdenData1->Projection(0);
      gStyle->SetOptStat(kTRUE);

      //K
      AliCFContainer *data2 = (AliCFContainer*) (d->Get("containerpp13TeV_K"));

      AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
      THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();  
      THnSparse* RbnumData2 = (THnSparse*)numData2->Clone("numNew2");
      RbnumData2->Reset();
      TAxis* axis5 = (TAxis*)RbnumData2->GetAxis(0);
      axis5->Set(newLimitsBins,newLimits);
      RbnumData2->SetBinEdges(0,newLimits);
      RbnumData2->RebinnedAdd(numData2, 1);
      TH1D* Rbhnum2 = (TH1D*)RbnumData2->Projection(0);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparseden2 = (AliCFGridSparse*)data2->GetGrid(6); // Reco
      THnSparse* denData2 = (THnSparse*)gridSparseden2->GetGrid();        
      THnSparse* RbdenData2 = (THnSparse*)denData2->Clone("denNew2");
      RbdenData2->Reset();
      TAxis* axis6 = (TAxis*)RbdenData2->GetAxis(0);
      axis6->Set(newLimitsBins,newLimits);
      RbdenData2->SetBinEdges(0,newLimits);
      RbdenData2->RebinnedAdd(denData2, 1);
      TH1D* Rbhden2 = (TH1D*)RbdenData2->Projection(0);
      gStyle->SetOptStat(kTRUE);
      
      //p
      AliCFContainer *data3 = (AliCFContainer*) (d->Get("containerpp13TeV_p"));

      AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
      THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();        
      THnSparse* RbnumData3 = (THnSparse*)numData3->Clone("numNew3");
      RbnumData3->Reset();
      TAxis* axis7 = (TAxis*)RbnumData3->GetAxis(0);
      axis7->Set(newLimitsBins,newLimits);
      RbnumData3->SetBinEdges(0,newLimits);
      RbnumData3->RebinnedAdd(numData3, 1);
      TH1D* Rbhnum3 = (TH1D*)RbnumData3->Projection(0);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparseden3 = (AliCFGridSparse*)data3->GetGrid(6); // Reco
      THnSparse* denData3 = (THnSparse*)gridSparseden3->GetGrid();          
      THnSparse* RbdenData3 = (THnSparse*)denData3->Clone("denNew3");
      RbdenData3->Reset();
      TAxis* axis8 = (TAxis*)RbdenData3->GetAxis(0);
      axis8->Set(newLimitsBins,newLimits);
      RbdenData3->SetBinEdges(0,newLimits);
      RbdenData3->RebinnedAdd(denData3, 1);
      TH1D* Rbhden3 = (TH1D*)RbdenData3->Projection(0);
      gStyle->SetOptStat(kTRUE);

      //Sigma+
      AliCFContainer *data4 = (AliCFContainer*) (d->Get("containerpp13TeV_SigmaP"));

      AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
      THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();        
      THnSparse* RbnumData4 = (THnSparse*)numData4->Clone("numNew4");
      RbnumData4->Reset();
      TAxis* axis9 = (TAxis*)RbnumData4->GetAxis(0);
      axis9->Set(newLimitsBins,newLimits);
      RbnumData4->SetBinEdges(0,newLimits);
      RbnumData4->RebinnedAdd(numData4, 1);
      TH1D* Rbhnum4 = (TH1D*)RbnumData4->Projection(0);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparseden4 = (AliCFGridSparse*)data4->GetGrid(6); // Reco
      THnSparse* denData4 = (THnSparse*)gridSparseden4->GetGrid();          
      THnSparse* RbdenData4 = (THnSparse*)denData4->Clone("denNew4");
      RbdenData4->Reset();
      TAxis* axis10 = (TAxis*)RbdenData4->GetAxis(0);
      axis10->Set(newLimitsBins,newLimits);
      RbdenData4->SetBinEdges(0,newLimits);
      RbdenData4->RebinnedAdd(denData4, 1);
      TH1D* Rbhden4 = (TH1D*)RbdenData4->Projection(0);
      gStyle->SetOptStat(kTRUE);

      //Sigma-
      AliCFContainer *data5 = (AliCFContainer*) (d->Get("containerpp13TeV_SigmaM"));

      AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
      THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();        
      THnSparse* RbnumData5 = (THnSparse*)numData5->Clone("numNew5");
      RbnumData5->Reset();
      TAxis* axis11 = (TAxis*)RbnumData5->GetAxis(0);
      axis11->Set(newLimitsBins,newLimits);
      RbnumData5->SetBinEdges(0,newLimits);
      RbnumData5->RebinnedAdd(numData5, 1);
      TH1D* Rbhnum5 = (TH1D*)RbnumData5->Projection(0);
      gStyle->SetOptStat(kTRUE);

      AliCFGridSparse* gridSparseden5 = (AliCFGridSparse*)data5->GetGrid(6); // Reco
      THnSparse* denData5 = (THnSparse*)gridSparseden5->GetGrid();        
      THnSparse* RbdenData5 = (THnSparse*)denData5->Clone("denNew5");
      RbdenData5->Reset();
      TAxis* axis12 = (TAxis*)RbdenData5->GetAxis(0);
      axis12->Set(newLimitsBins,newLimits);
      RbdenData5->SetBinEdges(0,newLimits);
      RbdenData5->RebinnedAdd(denData5, 1);
      TH1D* Rbhden5 = (TH1D*)RbdenData5->Projection(0);
      gStyle->SetOptStat(kTRUE);
    
      Rbhden = (TH1D*)RbhdenAll->Clone("Rbhdentot");
      Rbhden->Add(Rbhden1,-1);
      Rbhden->Add(Rbhden2,-1);
      Rbhden->Add(Rbhden3,-1);
      Rbhden->Add(Rbhden4,-1);
      Rbhden->Add(Rbhden5,-1);
      Rbhnum = (TH1D*)RbhnumAll->Clone("Rbhnumtot");
      Rbhnum->Add(Rbhnum1,-1);
      Rbhnum->Add(Rbhnum2,-1);
      Rbhnum->Add(Rbhnum3,-1);
      Rbhnum->Add(Rbhnum4,-1);
      Rbhnum->Add(Rbhnum5,-1);

    } else {
      printf("Error! Specie not supported!\n");
      return;
    }    
    
    TH1D* Rbheff = (TH1D*)Rbhden->Clone("heff_rebin");
    Rbheff->Divide(Rbhden,Rbhnum,1,1,"B");
    gStyle->SetOptStat(kTRUE);
    c->cd(2);
    Rbheff->SetMarkerStyle(21);
    Rbheff->SetMarkerSize(0.8);    
    Rbheff->Draw();
    
    // Num bins checking
    Bool_t Rebintest = kTRUE;
    Float_t sum = 0, sumnew = 0;
    Float_t content = 0, content2 = 0;
    if(Rebintest){
        Printf("Original THnSparse <-------------");
        for (Int_t ibin = 1; ibin<=hnum->GetNbinsX(); ibin++){
            content=0;
            for (Int_t ibin2 = 1; ibin2<=hnum->GetNbinsY(); ibin2++){
                for (Int_t ibin3 = 1; ibin3<=hnum->GetNbinsZ(); ibin3++){
                    content+=hnum->GetBinContent(ibin,ibin2,ibin3);
                }
            }
            Printf("ibin = %d, low edge = %0.2f, content = %0.2f",ibin,hnum->GetXaxis()->GetBinLowEdge(ibin), content);
            sum+=content;
            
        }
        Printf("Rebinned THnSparse <-------------");
        for (Int_t ibin = 1; ibin<=Rbhnum->GetNbinsX(); ibin++){
            content2=0;
            for (Int_t ibin2 = 1; ibin2<=Rbhnum->GetNbinsY(); ibin2++){
                for (Int_t ibin3 = 1; ibin3<=Rbhnum->GetNbinsZ(); ibin3++){
                    content2+=Rbhnum->GetBinContent(ibin,ibin2,ibin3);
                }
            }
            Printf("ibin = %d, low edge = %0.2f, content = %0.2f",ibin,Rbhnum->GetXaxis()->GetBinLowEdge(ibin), content2);
            sumnew+=content2;
        }
        Printf("Sum_before = %f, Sum_After = %f, Checking difference.. =%f%s", sum, sumnew, (sum-sumnew)*100/sum,"%");
    }
    
    
    if(((sum-sumnew)*100/sum)<.01){
        TFile* fout = new TFile(Form("1D_TrackingEffMap_%s.root",specie.Data()),"RECREATE");
        c->Write();
        fout->Close();
        printf( "1>>> rebin OKAY, efficiency map file created");
        c->SaveAs(Form("1D_TrackingEffMap_%s.png",specie.Data()));
        //CheckFileQuality(specie);
    }
    else Printf(" ======> rebin Failed");
}

void DrawRelativeAbundances_VsPt_7Species(TString filename="AnalysisResults.root") {
    
    TFile* f = new TFile(filename.Data(),"read");
    TDirectoryFile* dtot = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *datatot = (AliCFContainer*) (dtot->Get("containerpp13TeV_All"));
    AliCFGridSparse* gridSparsenum_tot = (AliCFGridSparse*)datatot->GetGrid(1); // GenAcc
    THnSparse* numData_tot = (THnSparse*)gridSparsenum_tot->GetGrid();

    TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data = (AliCFContainer*) (d->Get("containerpp13TeV_pi"));
    AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
    THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
    
    TDirectoryFile* d2 = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data2 = (AliCFContainer*) (d2->Get("containerpp13TeV_K"));
    AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
    THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();
    
    TDirectoryFile* d3 = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data3 = (AliCFContainer*) (d3->Get("containerpp13TeV_p"));
    AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
    THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();
    
    TDirectoryFile* d4 = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data4 = (AliCFContainer*) (d4->Get("containerpp13TeV_e"));
    AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
    THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();
    
    TDirectoryFile* d5 = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data5 = (AliCFContainer*) (d5->Get("containerpp13TeV_mu"));
    AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
    THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();

    TDirectoryFile* d6 = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data6 = (AliCFContainer*) (d6->Get("containerpp13TeV_SigmaP"));
    AliCFGridSparse* gridSparsenum6 = (AliCFGridSparse*)data6->GetGrid(1); // GenAcc
    THnSparse* numData6 = (THnSparse*)gridSparsenum6->GetGrid();

    TDirectoryFile* d7 = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data7 = (AliCFContainer*) (d7->Get("containerpp13TeV_SigmaM"));
    AliCFGridSparse* gridSparsenum7 = (AliCFGridSparse*)data7->GetGrid(1); // GenAcc
    THnSparse* numData7 = (THnSparse*)gridSparsenum7->GetGrid();
            
    TH1D* hnum_total = (TH1D*)numData_tot->Projection(0);
    TH1D* hnum_pi = (TH1D*)numData->Projection(0); 
    TH1D* hnum_k  = (TH1D*)numData2->Projection(0);    
    TH1D* hnum_p  = (TH1D*)numData3->Projection(0);    
    TH1D* hnum_e  = (TH1D*)numData4->Projection(0);    
    TH1D* hnum_mu = (TH1D*)numData5->Projection(0);
    TH1D* hnum_Sp = (TH1D*)numData6->Projection(0);
    TH1D* hnum_Sm = (TH1D*)numData7->Projection(0);

    TCanvas *c = new TCanvas();
    hnum_pi->Divide(hnum_pi,hnum_total,1,1,"B");
    hnum_pi->SetMarkerColor(kRed);
    hnum_pi->SetLineColor(kRed);
    hnum_pi->SetMarkerStyle(21);
    hnum_pi->SetMarkerSize(0.8);
    hnum_pi->SetMinimum(0.0001);
    hnum_pi->SetMaximum(1.8);
    hnum_pi->Draw();
    hnum_k->Divide(hnum_k,hnum_total,1,1,"B");
    hnum_k->SetMarkerColor(kBlue);
    hnum_k->SetLineColor(kBlue);
    hnum_k->SetMarkerStyle(21);
    hnum_k->SetMarkerSize(0.8);    
    hnum_k->Draw("same");
    hnum_p->Divide(hnum_p,hnum_total,1,1,"B");
    hnum_p->SetMarkerColor(kGreen+1);
    hnum_p->SetLineColor(kGreen+1);
    hnum_p->SetMarkerStyle(21);
    hnum_p->SetMarkerSize(0.8);    
    hnum_p->Draw("same");
    hnum_e->Divide(hnum_e,hnum_total,1,1,"B");
    hnum_e->SetMarkerColor(kMagenta);
    hnum_e->SetLineColor(kMagenta);
    hnum_e->SetMarkerStyle(21);
    hnum_e->SetMarkerSize(0.8);    
    hnum_e->Draw("same");
    hnum_mu->Divide(hnum_mu,hnum_total,1,1,"B");
    hnum_mu->SetMarkerColor(kCyan+1);
    hnum_mu->SetLineColor(kCyan+1);
    hnum_mu->SetMarkerStyle(21);
    hnum_mu->SetMarkerSize(0.8);    
    hnum_mu->Draw("same");
    hnum_Sp->Divide(hnum_mu,hnum_total,1,1,"B");
    hnum_Sp->SetMarkerColor(11);
    hnum_Sp->SetLineColor(11);
    hnum_Sp->SetMarkerStyle(21);
    hnum_Sp->SetMarkerSize(0.8);    
    hnum_Sp->Draw("same");
    hnum_Sm->Divide(hnum_mu,hnum_total,1,1,"B");
    hnum_Sm->SetMarkerColor(12);
    hnum_Sm->SetLineColor(12);
    hnum_Sm->SetMarkerStyle(21);
    hnum_Sm->SetMarkerSize(0.8);    
    hnum_Sm->Draw("same");

    TLegend* leg=new TLegend(0.65,0.15,0.95,0.40);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hnum_pi,"pions","epl");
    leg->AddEntry(hnum_k,"kaons","epl");
    leg->AddEntry(hnum_p,"protons","epl");
    leg->AddEntry(hnum_e,"electrons","epl");
    leg->AddEntry(hnum_mu,"muons","epl");
    leg->AddEntry(hnum_Sp,"sigmaP","epl");
    leg->AddEntry(hnum_Sm,"sigmaM","epl");        
    leg->Draw(); 

    c->SetLogy();
    c->Draw();

    c->SaveAs(Form("RelAbundances_7Species.png"));
    c->SaveAs(Form("RelAbundances_7Species.root"));
}

void DrawRelativeAbundances_VsPt_5Species(TString filename="AnalysisResults.root") {
    
    TFile* f = new TFile(filename.Data(),"read");
    TDirectoryFile* dtot = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *datatot = (AliCFContainer*) (dtot->Get("containerpp13TeV_All"));
    AliCFGridSparse* gridSparsenum_tot = (AliCFGridSparse*)datatot->GetGrid(1); // GenAcc
    THnSparse* numData_tot = (THnSparse*)gridSparsenum_tot->GetGrid();

    TDirectoryFile* d = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data = (AliCFContainer*) (d->Get("containerpp13TeV_pi"));
    AliCFGridSparse* gridSparsenum = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
    THnSparse* numData = (THnSparse*)gridSparsenum->GetGrid();
    
    TDirectoryFile* d2 = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data2 = (AliCFContainer*) (d2->Get("containerpp13TeV_K"));
    AliCFGridSparse* gridSparsenum2 = (AliCFGridSparse*)data2->GetGrid(1); // GenAcc
    THnSparse* numData2 = (THnSparse*)gridSparsenum2->GetGrid();
    
    TDirectoryFile* d3 = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data3 = (AliCFContainer*) (d3->Get("containerpp13TeV_p"));
    AliCFGridSparse* gridSparsenum3 = (AliCFGridSparse*)data3->GetGrid(1); // GenAcc
    THnSparse* numData3 = (THnSparse*)gridSparsenum3->GetGrid();
    
    TDirectoryFile* d4 = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data4 = (AliCFContainer*) (d4->Get("containerpp13TeV_e"));
    AliCFGridSparse* gridSparsenum4 = (AliCFGridSparse*)data4->GetGrid(1); // GenAcc
    THnSparse* numData4 = (THnSparse*)gridSparsenum4->GetGrid();
    
    TDirectoryFile* d5 = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data5 = (AliCFContainer*) (d5->Get("containerpp13TeV_mu"));
    AliCFGridSparse* gridSparsenum5 = (AliCFGridSparse*)data5->GetGrid(1); // GenAcc
    THnSparse* numData5 = (THnSparse*)gridSparsenum5->GetGrid();

    TDirectoryFile* d6 = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data6 = (AliCFContainer*) (d6->Get("containerpp13TeV_SigmaP"));
    AliCFGridSparse* gridSparsenum6 = (AliCFGridSparse*)data6->GetGrid(1); // GenAcc
    THnSparse* numData6 = (THnSparse*)gridSparsenum6->GetGrid();

    TDirectoryFile* d7 = (TDirectoryFile*)f->Get("PWGPP_CFSingleTrack");
    AliCFContainer *data7 = (AliCFContainer*) (d7->Get("containerpp13TeV_SigmaM"));
    AliCFGridSparse* gridSparsenum7 = (AliCFGridSparse*)data7->GetGrid(1); // GenAcc
    THnSparse* numData7 = (THnSparse*)gridSparsenum7->GetGrid();
            
    TH1D* hnum_total_7sp = (TH1D*)numData_tot->Projection(0);
    TH1D* hnum_pi = (TH1D*)numData->Projection(0); 
    TH1D* hnum_k  = (TH1D*)numData2->Projection(0);    
    TH1D* hnum_p  = (TH1D*)numData3->Projection(0);    
    TH1D* hnum_e  = (TH1D*)numData4->Projection(0);    
    TH1D* hnum_mu = (TH1D*)numData5->Projection(0);
    TH1D* hnum_Sp = (TH1D*)numData6->Projection(0);
    TH1D* hnum_Sm = (TH1D*)numData7->Projection(0);

    TH1D* hnum_total = (TH1D*)hnum_total_7sp->Clone("htotal_5sp");
    hnum_total->Add(hnum_Sp,-1);
    hnum_total->Add(hnum_Sm,-1);

    TCanvas *c = new TCanvas();
    hnum_pi->Divide(hnum_pi,hnum_total,1,1,"B");
    hnum_pi->SetMarkerColor(kRed);
    hnum_pi->SetLineColor(kRed);
    hnum_pi->SetMarkerStyle(21);
    hnum_pi->SetMarkerSize(0.8);
    hnum_pi->SetMinimum(0.0001);
    hnum_pi->SetMaximum(1.8);
    hnum_pi->Draw();
    hnum_k->Divide(hnum_k,hnum_total,1,1,"B");
    hnum_k->SetMarkerColor(kBlue);
    hnum_k->SetLineColor(kBlue);
    hnum_k->SetMarkerStyle(21);
    hnum_k->SetMarkerSize(0.8);    
    hnum_k->Draw("same");
    hnum_p->Divide(hnum_p,hnum_total,1,1,"B");
    hnum_p->SetMarkerColor(kGreen+1);
    hnum_p->SetLineColor(kGreen+1);
    hnum_p->SetMarkerStyle(21);
    hnum_p->SetMarkerSize(0.8);    
    hnum_p->Draw("same");
    hnum_e->Divide(hnum_e,hnum_total,1,1,"B");
    hnum_e->SetMarkerColor(kMagenta);
    hnum_e->SetLineColor(kMagenta);
    hnum_e->SetMarkerStyle(21);
    hnum_e->SetMarkerSize(0.8);    
    hnum_e->Draw("same");
    hnum_mu->Divide(hnum_mu,hnum_total,1,1,"B");
    hnum_mu->SetMarkerColor(kCyan+1);
    hnum_mu->SetLineColor(kCyan+1);
    hnum_mu->SetMarkerStyle(21);
    hnum_mu->SetMarkerSize(0.8);    
    hnum_mu->Draw("same");

    TLegend* leg=new TLegend(0.65,0.15,0.95,0.40);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hnum_pi,"pions","epl");
    leg->AddEntry(hnum_k,"kaons","epl");
    leg->AddEntry(hnum_p,"protons","epl");
    leg->AddEntry(hnum_e,"electrons","epl");
    leg->AddEntry(hnum_mu,"muons","epl");    
    leg->Draw(); 

    c->SetLogy();
    c->Draw();

    c->SaveAs(Form("RelAbundances_5Species.png"));
    c->SaveAs(Form("RelAbundances_5Species.root"));
}

/*
//additional lines 
 1. macro compliation
 root [0] gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_PHYSICS/macros -I$ALICE_PHYSICS/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF/ -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWGHF/correlationHF -I$ALICE_PHYSICS/PWGHF/correlationHF/macros  -g");
 root [1] gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/LoadLibraries.C")
 root [2] .L DrawDhCorr_STE_3DEfficiency_in_pTEtaZvtx.C++
 
 
 2. Getting histogram from file
 TFile* fAssoTracksEffMap = TFile::Open("STE3D_DhCorrPP_Pass4_03Dec15.root");
 TCanvas *cTrackEffMap = (TCanvas*)fAssoTracksEffMap->Get("c");
 TCanvas *c = new TCanvas ("c", "Eff histogram from root file", 500, 600);
 TH3D *hEffTrackMap = (TH3D*)cTrackEffMap->FindObject("heff_rebin");
 hEffTrackMap->Draw();
 */

