// provided by Gamma Conversion Group, PWG4, Kathrin Koch, kkoch@physi.uni-heidelberg.de and Friederike Bock, fbock@physi.uni-heidelberg.de

#include <Riostream.h>
#include <fstream>
#include "PlottingGammaHistos.h"
using namespace std;

extern TRandom *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem *gSystem;

void  Plot_Mapping_Histos_Events(const char *input1 = "myOutput", const char *secondinput = "", const char *cutsel = "", const char *path = "", const char *output = "Mapping", const char *Plots = "kFALSE", const char *suffix = "gif"){	
	
	gROOT->Reset();
	gROOT->SetStyle("Plain");
	StyleSettings();	

		// Which file you want to analyse
	Char_t filename_input1[200] = (Form("%s%s",path,input1));
	
	Bool_t SinglePlots = kFALSE;
	if(Plots == "kTRUE") SinglePlots = kTRUE;
	
	// the outputfile
	TPostScript *ps_mapping = 0x0;
	if(!SinglePlots)ps_mapping = new TPostScript(Form("%s%s.ps",path,output),111);
	
	//rebinning
	const int rebin = 4;

	//Which Minimum do you want to set for XY-2Dplot
	Int_t minimumXY = 10;

	//How big should the right margin in 2D plots be? Others are set by default setting.
	Float_t RightMargin = 0.17;

	//How many slices do you have?
//	const int NbinsPhi = 8;  
//	const int NbinsR = 12;
	const int NbinsZ = 6;
	const int NbinsR = 13;	
 	const int NbinsZ = 12;

	// how many raws and colums you want to have in your output?
	Int_t column = 2;
	Int_t raw = 7;
	
	//Array for ZinR-Ranges
	Float_t *RangeZinR[NbinsR] = {30,30,40,50,60,70,80,100,120,140,160,180,240}; 
	//Float_t *RangeZinR[NbinsR] = {30,40,50,60,70,80,100,120,140,160,180,240};
	
	//Array of Rbins
	Float_t ArrayRbins[12] = {3.5,5.75,9.5,13.,21.,27.5,35.,42.,55.,81.5,90,72.};
	//Array of Zbins
	Float_t ArrayZbins[11] = {0,15,30,50,100,200,-15,-30,-50,-100,-200};

	//Array defintion for printing Logo in right upper corner
	Float_t right_up[4]={0.7,0.63,0.15, 0.02};
	Float_t right_down[4]={0.7,0.23,0.15, 0.02};
	Float_t right_up2D[4]={0.68,0.775,0.11, 0.02};
	//Array defintion for printing Logo in left upper corner
	Float_t left_up[4]={0.17,0.73, 0.15, 0.02};
	Float_t left_down[4] = {0.15,0.17, 0.15, 0.02};
	//Array defintion for printing text in right upper corner
	Float_t right_up_text[4]={0.6,0.8,0.15, 0.04};

	// get the histos
	TFile f1(filename_input1);  


	// choice of dateset
	if(cutsel != ""){
		char *GammaDirectory = Form("PWG4_GammaConversion_%s",  cutsel);
		cout << GammaDirectory << endl;
		char *GammaList = Form("histogramsAliGammaConversion_%s", cutsel);
		cout << GammaList << endl;
	}else{
		char *GammaDirectory = "PWG4_GammaConversion";
		cout << GammaDirectory << endl;
		char *GammaList = "histogramsAliGammaConversion";
		cout << GammaList << endl;
	}	

	// labeling
	char *StandardYAxis = "#gamma/ event scaled by multiplicity";
	char *Date = "3rd June 2010";
	TLatex *EtaRange = new TLatex(0.15,0.845,"|#eta| < 0.9 "); // Bo: this was modified
	            EtaRange->SetNDC();
	            EtaRange->SetTextFont(62);
	            EtaRange->SetTextSize(0.04);
	            EtaRange->SetLineWidth(6);	

//------------------------------- Reading FILES	----------------------------------------------------------------------
	 TDirectory *fPWG4GammaConversion_input1 = new TDirectory(); // definition of first folder / list	
	 TList *fHistosGammaConversion_input1 = new TList(); // definition of first folder / list
	 TList *fESDContainer_input1 = new TList();  // definition of following folder / list
	 TList *fMappingContainer_input1 = new TList();

	fPWG4GammaConversion_input1 = (TDirectory*)f1.Get(GammaDirectory); 
	fHistosGammaConversion_input1 = (TList*)fPWG4GammaConversion_input1->Get(GammaList); 
	fMappingContainer_input1 = (TList*)fHistosGammaConversion_input1->FindObject("Mapping histograms"); 
	fESDContainer_input1 = (TList*)fHistosGammaConversion_input1->FindObject("ESD histograms"); 
	
	TH1F *ESD_Conversion_R_input1=fESDContainer_input1->FindObject("ESD_Conversion_R");        
	TH2F *ESD_Conversion_ZR_input1=fESDContainer_input1->FindObject("ESD_Conversion_ZR");
	TH2F *ESD_Conversion_XY_input1=fESDContainer_input1->FindObject("ESD_Conversion_XY");
	TH1F *ESD_Conversion_OpeningAngle_input1=fESDContainer_input1->FindObject("ESD_Conversion_OpeningAngle");
	TH1D *ESD_Conversion_Z_input1=ESD_Conversion_ZR_input1->ProjectionX("ESD_Conversion_Z_input1");
	TH1F *ESD_NumberOfContributorsVtx_input1=fESDContainer_input1->FindObject("ESD_NumberOfContributorsVtx");
	TH1D *ESD_Conversion_Z_input1=ESD_Conversion_ZR_input1->ProjectionX("ESD_Cocdnversion_Z_input1");
	TH1F *ScalingDiagramm_input1= fMappingContainer_input1->FindObject("ESD_Conversion_Mapping_Phi_in_R_11");
	TH1F *ESD_NumberOfGoodESDTracks_input1=fESDContainer_input1->FindObject("ESD_NumberOfGoodESDTracksVtx");

	ESD_NumberOfContributorsVtx_input1->SetAxisRange(1.,100.);	
//	ESD_NumberOfGoodESDTracks_input1->SetAxisRange(1.,100.);

	Float_t Scaling_input1 = ScalingDiagramm_input1->Integral();
	Float_t nGoodEvents_input1 = ESD_NumberOfContributorsVtx_input1->Integral();
	Float_t nGoodTrig_input1 = ESD_NumberOfContributorsVtx_input1->GetEntries();
	Float_t nRecGamma_input1 = ESD_Conversion_R_input1->GetEntries();
	cout<< input1 << "    Number of events::   " << nGoodEvents_input1 << "    Number of triggers::   " << nGoodTrig_input1 << "    Number reconstructed gammas::    "<< nRecGamma_input1 <<endl;

	Double_t mean_input1 = ESD_NumberOfGoodESDTracks_input1->GetMean();

	Float_t normFacRec_input1=1./nGoodEvents_input1;
//	Float_t normFacRec_input1=1./Scaling_input1;
	//Scaling reconstr.

	GammaScalingHistogramm(ESD_Conversion_R_input1,normFacRec_input1);
//	GammaScalingHistogramm(ESD_Conversion_ZR_input1,normFacRec_input1);
//	GammaScalingHistogramm(ESD_Conversion_XY_input1,normFacRec_input1);
	GammaScalingHistogramm(ESD_Conversion_OpeningAngle_input1,normFacRec_input1);
	GammaScalingHistogramm(ESD_Conversion_Z_input1,normFacRec_input1);

	TH1F *ESD_Conversion_R_summed_input1 = (TH1F*) ESD_Conversion_R_input1->Clone("ESD_Conversion_R_summed_input1");
	ESD_Conversion_R_summed_input1->Reset();
	Int_t ConvBinsR = ESD_Conversion_R_input1->GetNbinsX();
	for(Int_t ConvR = 1; ConvR < ConvBinsR +1; ConvR++){
		if (ConvR == 1){
			ESD_Conversion_R_summed_input1->AddBinContent(ConvR, ESD_Conversion_R_input1-> GetBinContent(ConvR));
		} else {
			ESD_Conversion_R_summed_input1->AddBinContent(ConvR,(ESD_Conversion_R_summed_input1->GetBinContent(ConvR -1) + ESD_Conversion_R_input1->GetBinContent(ConvR)));		
		}
	}

	TH1F *ESD_Conversion_Mapping_Phi_in_R_input1[NbinsR];
	Char_t histoname_Phi_in_R_input1[100];
	for(Int_t iR = 0; iR < NbinsR; iR++){
		sprintf(histoname_Phi_in_R_input1,"ESD_Conversion_Mapping_Phi_in_R_%02d",iR);
		ESD_Conversion_Mapping_Phi_in_R_input1[iR] = dynamic_cast<TH1F*>(fMappingContainer_input1->FindObject(histoname_Phi_in_R_input1));
		GammaScalingHistogramm(ESD_Conversion_Mapping_Phi_in_R_input1[iR],normFacRec_input1);
		ESD_Conversion_Mapping_Phi_in_R_input1[iR]->Rebin(rebin);	
	}
	
	TH1F *ESD_Conversion_Mapping_Z_in_R_input1[NbinsR];
	Char_t histoname_Z_in_R_input1[100];
	for(Int_t iR = 0; iR < NbinsR; iR++){
		sprintf(histoname_Z_in_R_input1,"ESD_Conversion_Mapping_Z_in_R_%02d",iR);
		ESD_Conversion_Mapping_Z_in_R_input1[iR] = dynamic_cast<TH1F*>(fMappingContainer_input1->FindObject(histoname_Z_in_R_input1));
		GammaScalingHistogramm(ESD_Conversion_Mapping_Z_in_R_input1[iR],normFacRec_input1);
		ESD_Conversion_Mapping_Z_in_R_input1[iR]->Rebin(rebin);
	}
	
	TH1F *ESD_Conversion_Mapping_Phi_in_Z_input1[NbinsZ];
	Char_t histoname_Phi_in_Z_input1[100];
	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		sprintf(histoname_Phi_in_Z_input1,"ESD_Conversion_Mapping_Phi_in_Z_%02d",iZ);
		ESD_Conversion_Mapping_Phi_in_Z_input1[iZ] = dynamic_cast<TH1F*>(fMappingContainer_input1->FindObject(histoname_Phi_in_Z_input1));
		GammaScalingHistogramm(ESD_Conversion_Mapping_Phi_in_Z_input1[iZ],normFacRec_input1);
		ESD_Conversion_Mapping_Phi_in_Z_input1[iZ]->Rebin(rebin);
	}

	TH1F *ESD_Conversion_Mapping_R_in_Z_input1[NbinsZ];
	Char_t histoname_R_in_Z_input1[100];
	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		sprintf(histoname_R_in_Z_input1,"ESD_Conversion_Mapping_R_in_Z_%02d",iZ);
		ESD_Conversion_Mapping_R_in_Z_input1[iZ] = dynamic_cast<TH1F*>(fMappingContainer_input1->FindObject(histoname_R_in_Z_input1));
		GammaScalingHistogramm(ESD_Conversion_Mapping_R_in_Z_input1[iZ],normFacRec_input1);
		ESD_Conversion_Mapping_R_in_Z_input1[iZ]->Rebin(rebin);
	}
	
	// Middle Pt	
	TH1F *ESD_Conversion_Mapping_MidPt_Phi_in_R_input1[NbinsR];
	Char_t histoname_MidPt_Phi_in_R_input1[100];
	for(Int_t iR = 0; iR < NbinsR; iR++){
		sprintf(histoname_MidPt_Phi_in_R_input1,"ESD_Conversion_Mapping_MidPt_Phi_in_R_%02d",iR);
		ESD_Conversion_Mapping_MidPt_Phi_in_R_input1[iR] = dynamic_cast<TH1F*>(fMappingContainer_input1->FindObject(histoname_MidPt_Phi_in_R_input1));
		GammaScalingHistogramm(ESD_Conversion_Mapping_MidPt_Phi_in_R_input1[iR],normFacRec_input1);
		ESD_Conversion_Mapping_MidPt_Phi_in_R_input1[iR] ->Rebin(rebin);
	}
	
	TH1F *ESD_Conversion_Mapping_MidPt_Z_in_R_input1[NbinsR];
	Char_t histoname_MidPt_Z_in_R_input1[100];
	for(Int_t iR = 0; iR < NbinsR; iR++){
		sprintf(histoname_MidPt_Z_in_R_input1,"ESD_Conversion_Mapping_MidPt_Z_in_R_%02d",iR);
		ESD_Conversion_Mapping_MidPt_Z_in_R_input1[iR] = dynamic_cast<TH1F*>(fMappingContainer_input1->FindObject(histoname_MidPt_Z_in_R_input1));
		GammaScalingHistogramm(ESD_Conversion_Mapping_MidPt_Z_in_R_input1[iR],normFacRec_input1);
		ESD_Conversion_Mapping_MidPt_Z_in_R_input1[iR] ->Rebin(rebin);
	}
	
	TH1F *ESD_Conversion_Mapping_MidPt_Phi_in_Z_input1[NbinsZ];
	Char_t histoname_MidPt_Phi_in_Z_input1[100];
	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		sprintf(histoname_MidPt_Phi_in_Z_input1,"ESD_Conversion_Mapping_MidPt_Phi_in_Z_%02d",iZ);
		ESD_Conversion_Mapping_MidPt_Phi_in_Z_input1[iZ] = dynamic_cast<TH1F*>(fMappingContainer_input1->FindObject(histoname_MidPt_Phi_in_Z_input1));
		GammaScalingHistogramm(ESD_Conversion_Mapping_MidPt_Phi_in_Z_input1[iZ],normFacRec_input1);
		ESD_Conversion_Mapping_MidPt_Phi_in_Z_input1[iZ] ->Rebin(rebin);
	}
	
	TH1F *ESD_Conversion_Mapping_MidPt_R_in_Z_input1[NbinsZ];
	Char_t histoname_MidPt_R_in_Z_input1[100];
	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		sprintf(histoname_MidPt_R_in_Z_input1,"ESD_Conversion_Mapping_MidPt_R_in_Z_%02d",iZ);
		ESD_Conversion_Mapping_MidPt_R_in_Z_input1[iZ] = dynamic_cast<TH1F*>(fMappingContainer_input1->FindObject(histoname_MidPt_R_in_Z_input1));
		GammaScalingHistogramm(ESD_Conversion_Mapping_MidPt_R_in_Z_input1[iZ],normFacRec_input1);
		ESD_Conversion_Mapping_MidPt_R_in_Z_input1[iZ] ->Rebin(rebin);
	}

	// ------------------------------------- second file-----------------------------------------------------------
	TFile *input2 = 0x0 ;
	
	TH1F *ESD_Conversion_Mapping_Phi_in_R_input2[NbinsR];
	TH1F *ESD_Conversion_Mapping_Z_in_R_input2[NbinsR];
	TH1F *ESD_Conversion_Mapping_Phi_in_Z_input2[NbinsZ];
	TH1F *ESD_Conversion_Mapping_R_in_Z_input2[NbinsZ];
	TH1F *ESD_Conversion_Mapping_MidPt_Phi_in_R_input2[NbinsR];
	TH1F *ESD_Conversion_Mapping_MidPt_Z_in_R_input2[NbinsR];
	TH1F *ESD_Conversion_Mapping_MidPt_Phi_in_Z_input2[NbinsZ];
	TH1F *ESD_Conversion_Mapping_MidPt_R_in_Z_input2[NbinsZ];
	
	if(secondinput != ""){
		
		input2 = new TFile(Form("%s%s",path, secondinput));	

		 TDirectory *fPWG4GammaConversion_input2 = new TDirectory(); // definition of first folder / list	
		 TList *fHistosGammaConversion_input2 = new TList(); // definition of first folder / list
		 TList *fMappingContainer_input2 = new TList();  // definition of following folder / list
		 TList *fESDContainer_input2 = new TList();  // definition of following folder / list
		 
		fPWG4GammaConversion_input2 = (TDirectory*)input2->Get(GammaDirectory); 
		fHistosGammaConversion_input2 = (TList*)fPWG4GammaConversion_input2->Get(GammaList); 
		fMappingContainer_input2 = (TList*)fHistosGammaConversion_input2->FindObject("Mapping histograms");
		fESDContainer_input2 = (TList*)fHistosGammaConversion_input2->FindObject("ESD histograms"); 

		TH1F * ESD_NumberOfGoodESDTracks_input2=fESDContainer_input2->FindObject("ESD_NumberOfGoodESDTracksVtx");			
		TH1F *ESD_Conversion_R_input2=fESDContainer_input2->FindObject("ESD_Conversion_R");        
		TH2F *ESD_Conversion_ZR_input2=fESDContainer_input2->FindObject("ESD_Conversion_ZR");
		TH2F *ESD_Conversion_XY_input2=fESDContainer_input2->FindObject("ESD_Conversion_XY");
		TH1F *ESD_Conversion_OpeningAngle_input2=fESDContainer_input2->FindObject("ESD_Conversion_OpeningAngle");
		TH1D *ESD_Conversion_Z_input2=ESD_Conversion_ZR_input2->ProjectionX("ESD_Conversion_Z_input2");
		TH1F * ESD_NumberOfContributorsVtx_input2=fESDContainer_input2->FindObject("ESD_NumberOfContributorsVtx");
		TH1F *ScalingDiagramm_input2= fMappingContainer_input2->FindObject("ESD_Conversion_Mapping_Phi_in_R_11");

		ESD_NumberOfContributorsVtx_input2->SetAxisRange(1.,100.);
	//	ESD_NumberOfGoodESDTracks_input2->SetAxisRange(1.,100.);
	
		Float_t Scaling_input2 = ScalingDiagramm_input2->Integral();
		Float_t nGoodEvents_input2 = ESD_NumberOfContributorsVtx_input2->Integral();
		Float_t nGoodTrig_input2 = ESD_NumberOfContributorsVtx_input2->GetEntries();
		Float_t nRecGamma_input2 = ESD_Conversion_R_input2->GetEntries();
		cout<< secondinput << "    Number of events::   " << nGoodEvents_input2 << "    Number of triggers::   " << nGoodTrig_input2 << "    Number reconstructed gammas::    "<< nRecGamma_input2 <<endl;

		Double_t mean_input2 = ESD_NumberOfGoodESDTracks_input2->GetMean();

	//	Float_t normFacRec_input2=1./nGoodEvents_input2;
		Float_t normFacRec_input2=1./nGoodEvents_input2 * mean_input1/mean_input2;
		
		//Scaling reconstr.
		GammaScalingHistogramm(ESD_Conversion_R_input2,normFacRec_input2);
//		GammaScalingHistogramm(ESD_Conversion_ZR_input2,normFacRec_input2);
//		GammaScalingHistogramm(ESD_Conversion_XY_input2,normFacRec_input2);
		GammaScalingHistogramm(ESD_Conversion_OpeningAngle_input2,normFacRec_input2);
		GammaScalingHistogramm(ESD_Conversion_Z_input2,normFacRec_input2);

		TH1F *ESD_Conversion_R_summed_input2 = (TH1F*) ESD_Conversion_R_input2->Clone("ESD_Conversion_R_summed_input2");
		ESD_Conversion_R_summed_input2->Reset();
		Int_t ConvBinsR = ESD_Conversion_R_input2->GetNbinsX();
		for(Int_t ConvR = 1; ConvR < ConvBinsR +1; ConvR++){
			if (ConvR == 1){
				ESD_Conversion_R_summed_input2->AddBinContent(ConvR, ESD_Conversion_R_input2-> GetBinContent(ConvR));
			} else {
				ESD_Conversion_R_summed_input2->AddBinContent(ConvR,(ESD_Conversion_R_summed_input2->GetBinContent(ConvR -1) + ESD_Conversion_R_input2->GetBinContent(ConvR)));		
			}
		}

		Char_t histoname_Phi_in_R_input2[100];
		Char_t histoname_Ratio_Phi_in_R_input2[100];
		TH1F *Ratio_Mapping_Phi_in_R[NbinsR];
		for(Int_t iR = 0; iR < NbinsR; iR++){
			sprintf(histoname_Phi_in_R_input2,"ESD_Conversion_Mapping_Phi_in_R_%02d",iR);
			ESD_Conversion_Mapping_Phi_in_R_input2[iR] = dynamic_cast<TH1F*>(fMappingContainer_input2->FindObject(histoname_Phi_in_R_input2));
			GammaScalingHistogramm(ESD_Conversion_Mapping_Phi_in_R_input2[iR],normFacRec_input2);
			ESD_Conversion_Mapping_Phi_in_R_input2[iR]->Rebin(rebin);
	
			sprintf(histoname_Ratio_Phi_in_R_input2,"Ratio_Mapping_Phi_in_R_%02d",iR);
			Ratio_Mapping_Phi_in_R[iR]= (TH1F*)ESD_Conversion_Mapping_Phi_in_R_input1[iR]->Clone();
			Ratio_Mapping_Phi_in_R[iR]->SetName(histoname_Ratio_Phi_in_R_input2); 
			Ratio_Mapping_Phi_in_R[iR]->Divide(Ratio_Mapping_Phi_in_R[iR],ESD_Conversion_Mapping_Phi_in_R_input2[iR]);
		}
		
		Char_t histoname_Z_in_R_input2[100];
		Char_t histoname_Ratio_Z_in_R_input2[100];
		TH1F *Ratio_Mapping_Z_in_R[NbinsR];
		for(Int_t iR = 0; iR < NbinsR; iR++){
			sprintf(histoname_Z_in_R_input2,"ESD_Conversion_Mapping_Z_in_R_%02d",iR);
			ESD_Conversion_Mapping_Z_in_R_input2[iR] = dynamic_cast<TH1F*>(fMappingContainer_input2->FindObject(histoname_Z_in_R_input2));
			GammaScalingHistogramm(ESD_Conversion_Mapping_Z_in_R_input2[iR],normFacRec_input2);
			ESD_Conversion_Mapping_Z_in_R_input2[iR]->Rebin(rebin);
			sprintf(histoname_Ratio_Z_in_R_input2,"Ratio_Mapping_Z_in_R_%02d",iR);
			Ratio_Mapping_Z_in_R[iR]= (TH1F*)ESD_Conversion_Mapping_Z_in_R_input1[iR]->Clone();
			Ratio_Mapping_Z_in_R[iR]->SetName(histoname_Ratio_Z_in_R_input2); 
			Ratio_Mapping_Z_in_R[iR]->Divide(Ratio_Mapping_Z_in_R[iR],ESD_Conversion_Mapping_Z_in_R_input2[iR]);
		}
		
		Char_t histoname_Phi_in_Z_input2[100];
		Char_t histoname_Ratio_Phi_in_Z_input2[100];
		TH1F *Ratio_Mapping_Phi_in_Z[NbinsZ];
		for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
			sprintf(histoname_Phi_in_Z_input2,"ESD_Conversion_Mapping_Phi_in_Z_%02d",iZ);
			ESD_Conversion_Mapping_Phi_in_Z_input2[iZ] = dynamic_cast<TH1F*>(fMappingContainer_input2->FindObject(histoname_Phi_in_Z_input2));
			GammaScalingHistogramm(ESD_Conversion_Mapping_Phi_in_Z_input2[iZ],normFacRec_input2);
			ESD_Conversion_Mapping_Phi_in_Z_input2[iZ] ->Rebin(rebin);			
			sprintf(histoname_Ratio_Phi_in_Z_input2,"Ratio_Mapping_Phi_in_Z_%02d",iZ);
			Ratio_Mapping_Phi_in_Z[iZ]= (TH1F*)ESD_Conversion_Mapping_Phi_in_Z_input1[iZ]->Clone();
			Ratio_Mapping_Phi_in_Z[iZ]->SetName(histoname_Ratio_Phi_in_Z_input2); 
			Ratio_Mapping_Phi_in_Z[iZ]->Divide(Ratio_Mapping_Phi_in_Z[iZ],ESD_Conversion_Mapping_Phi_in_Z_input2[iZ]);
		}
		
		Char_t histoname_R_in_Z_input2[100];
		Char_t histoname_Ratio_R_in_Z_input2[100];
		TH1F *Ratio_Mapping_R_in_Z[NbinsZ];
		for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
			sprintf(histoname_R_in_Z_input2,"ESD_Conversion_Mapping_R_in_Z_%02d",iZ);
			ESD_Conversion_Mapping_R_in_Z_input2[iZ] = dynamic_cast<TH1F*>(fMappingContainer_input2->FindObject(histoname_R_in_Z_input2));
			GammaScalingHistogramm(ESD_Conversion_Mapping_R_in_Z_input2[iZ],normFacRec_input2);
			ESD_Conversion_Mapping_R_in_Z_input2[iZ] ->Rebin(rebin);
			sprintf(histoname_Ratio_Phi_in_Z_input2,"Ratio_Mapping_R_in_Z_%02d",iZ);
			Ratio_Mapping_R_in_Z[iZ]= (TH1F*)ESD_Conversion_Mapping_R_in_Z_input1[iZ]->Clone();
			Ratio_Mapping_R_in_Z[iZ]->SetName(histoname_Ratio_R_in_Z_input2); 
			Ratio_Mapping_R_in_Z[iZ]->Divide(Ratio_Mapping_R_in_Z[iZ],ESD_Conversion_Mapping_R_in_Z_input2[iZ]);
		}
		
		// Middle Pt	
		Char_t histoname_MidPt_Phi_in_R_input2[100];
		for(Int_t iR = 0; iR < NbinsR; iR++){
			sprintf(histoname_MidPt_Phi_in_R_input2,"ESD_Conversion_Mapping_MidPt_Phi_in_R_%02d",iR);
			ESD_Conversion_Mapping_MidPt_Phi_in_R_input2[iR] = dynamic_cast<TH1F*>(fMappingContainer_input2->FindObject(histoname_MidPt_Phi_in_R_input2));
			GammaScalingHistogramm(ESD_Conversion_Mapping_MidPt_Phi_in_R_input2[iR],normFacRec_input2);
			ESD_Conversion_Mapping_MidPt_Phi_in_R_input2[iR] ->Rebin(rebin);
		}
		
		Char_t histoname_MidPt_Z_in_R_input2[100];
		for(Int_t iR = 0; iR < NbinsR; iR++){
			sprintf(histoname_MidPt_Z_in_R_input2,"ESD_Conversion_Mapping_MidPt_Z_in_R_%02d",iR);
			ESD_Conversion_Mapping_MidPt_Z_in_R_input2[iR] = dynamic_cast<TH1F*>(fMappingContainer_input2->FindObject(histoname_MidPt_Z_in_R_input2));
			GammaScalingHistogramm(ESD_Conversion_Mapping_MidPt_Z_in_R_input2[iR],normFacRec_input2);
			ESD_Conversion_Mapping_MidPt_Z_in_R_input2[iR] ->Rebin(rebin);
		}
		
		Char_t histoname_MidPt_Phi_in_Z_input2[100];
		for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
			sprintf(histoname_MidPt_Phi_in_Z_input2,"ESD_Conversion_Mapping_MidPt_Phi_in_Z_%02d",iZ);
			ESD_Conversion_Mapping_MidPt_Phi_in_Z_input2[iZ] = dynamic_cast<TH1F*>(fMappingContainer_input2->FindObject(histoname_MidPt_Phi_in_Z_input2));
			GammaScalingHistogramm(ESD_Conversion_Mapping_MidPt_Phi_in_Z_input2[iZ],normFacRec_input2);
			ESD_Conversion_Mapping_MidPt_Phi_in_Z_input2[iZ] ->Rebin(rebin);
		}
		
		Char_t histoname_MidPt_R_in_Z_input2[100];
		for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
			sprintf(histoname_MidPt_R_in_Z_input2,"ESD_Conversion_Mapping_MidPt_R_in_Z_%02d",iZ);
			ESD_Conversion_Mapping_MidPt_R_in_Z_input2[iZ] = dynamic_cast<TH1F*>(fMappingContainer_input2->FindObject(histoname_MidPt_R_in_Z_input2));
			GammaScalingHistogramm(ESD_Conversion_Mapping_MidPt_R_in_Z_input2[iZ],normFacRec_input2);
			ESD_Conversion_Mapping_MidPt_R_in_Z_input2[iZ] ->Rebin(rebin);
		}	
	}
	// end second file
	
	
	// -----------------page 1---------------------------------------------------------------------------
	
	// plot the ESD histos
		
	TLine * linePhi = new TLine (-3.2,1,3.2,1);
	TLine * lineZ = new TLine (-300,1,300,1);
	TLine * lineR = new TLine (0,1,200,1);
	linePhi->SetLineColor(2);
	lineZ->SetLineColor(2);
	lineR->SetLineColor(2);

	leg1 = new TLegend(0.6,0.82,0.9,0.9);
	leg1->AddEntry(ESD_Conversion_R_input1,("Data"),"l");
	if(secondinput != ""){
	leg1->AddEntry(ESD_Conversion_R_input2,("MC"),"l");}
	leg1->SetFillColor(0);
	leg1->SetTextSize(0.04);


if (!SinglePlots) {	
	ps_mapping->NewPage();
//	c_0 = new TCanvas("c_0","",200,10,700,1000);  // gives the page size
	
	title0 = new TPaveLabel(0.05,0.92,0.95,0.96,(Form("Input1: %s",input1)));	
	title0->SetFillColor(16);
	title0->SetTextSize(0.25);
//	title0->Draw();
	
	if(secondinput != ""){
		title1 = new TPaveLabel(0.05,0.87,0.95,0.91,(Form("Input2: %s",secondinput)));	
		title1->SetFillColor(16);
		title1->SetTextSize(0.25);
//		title1->Draw();
	}
//	c_0->Update();
	
	// --------------------------------page 2 - R-distributions -----------------------------------		
	
	ps_mapping->NewPage();
	
	c_R = new TCanvas("c_R","",200,10,700,1000);  // gives the page size
	
	pad_R = new TPad("pad_R","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas  	
	pad_R->SetFillColor(0);
	pad_R->GetFrame()->SetFillColor(0);
	pad_R->SetBorderMode(0);
	pad_R->Divide(2,2);
	pad_R->Draw();
	
	title0->Draw(); if(secondinput != ""){title1->Draw();}
	
	pad_R->cd(1);
	pad_R->cd(1)->SetLogy(1);
		if(secondinput == ""){
				DrawAutoGammaHisto( ESD_Conversion_R_input1, 
								 "Conversions R distribution","R [cm]",StandardYAxis,
								 kTRUE, 1.5,0.00001,
								 kFALSE,0. ,0,
								 kTRUE, 0.,180.);
		}else{
				DrawAutoGammaHistos( ESD_Conversion_R_input1, 
								 ESD_Conversion_R_input2, 
								 "Conversions R distribution","R [cm]",StandardYAxis,
								 kTRUE, 1.5,0.00001,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,180.);
		}
		for(Int_t i=0; i < 12 ; i++){
			DrawGammaLines(ArrayRbins[i], ArrayRbins[i], 0.00001,ESD_Conversion_R_input1->GetMaximum());
		}		
		DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);	
		
	pad_R->cd(2);	
		if(secondinput == ""){
				DrawAutoGammaHisto( ESD_Conversion_R_input1, 
								 "Conversions R distribution","R [cm]",StandardYAxis,
								 kFALSE, 1.5,0,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,180.);
		}else{
				DrawAutoGammaHistos( ESD_Conversion_R_input1, 
								 ESD_Conversion_R_input2, 
								 "Conversions R distribution","R [cm] ",StandardYAxis,
								 kFALSE, 1.5,0,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,180.);
		}
		for(Int_t i=0; i < 12 ; i++){
			DrawGammaLines(ArrayRbins[i], ArrayRbins[i], 0.,ESD_Conversion_R_input1->GetMaximum());
		}

		DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);	

	
	//----------------------------- Integrated Radius -----------------------
	pad_R->cd(3);
	pad_R->cd(3)->SetLogy(1);
		if(secondinput == ""){
				DrawAutoGammaHisto( ESD_Conversion_R_summed_input1, 
								 "Conversions R distribution","R [cm]","Integrated (0-to-R) #gamma candidates/event scaled by multiplicity",
								 kTRUE, 1.5,0.000001,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,180.);
		}else{
				DrawAutoGammaHistos( ESD_Conversion_R_summed_input1, 
								 ESD_Conversion_R_summed_input2, 
								 "Integrated Radius of Conversion","R [cm] ","Integrated (0-to-R) #gamma candidates/event scaled by multiplicity",
								 kTRUE, 1.5,0.000001,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,180.);
		}
		DrawAliceLogo(right_down [0],right_down[1],right_down[2],right_down[3]);	
		
	pad_R->cd(4);
	pad_R->cd(4)->SetLogy(0);
		if(secondinput == ""){
				DrawAutoGammaHisto( ESD_Conversion_R_summed_input1, 
								 "Conversions R distribution","R [cm]","Integrated (0-to-R) #gamma candidates/event scaled by multiplicity",
								 kTRUE, 1.5,0,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,180.);
		}else{
				DrawAutoGammaHistos( ESD_Conversion_R_summed_input1, 
								 ESD_Conversion_R_summed_input2, 
								 "Integrated Radius of Conversion","R [cm] ","Integrated (0-to-R) #gamma candidates/event scaled by multiplicity",
								 kTRUE, 1.5,0,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,180.);
		}
		DrawAliceLogo(left_up [0],left_up[1],left_up[2],left_up[3]);	

	pad_R->Update();

	//--------------------------- Page 3  - Z-Distributions -----------------------------------------------
	ps_mapping->NewPage();
	
	c_Z = new TCanvas("c_Z","",200,10,700,1000);  // gives the page size
	
	pad_Z = new TPad("pad_Z","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_Z->SetFillColor(0);
	pad_Z->GetFrame()->SetFillColor(0);
	pad_Z->SetBorderMode(0);
	pad_Z->Divide(2,2);
	pad_Z->Draw();
	
	title0->Draw(); if(secondinput != ""){title1->Draw();}
	
	pad_Z->cd(1);
	pad_Z->cd(1)->SetLogy(1);
		if(secondinput == ""){
				DrawAutoGammaHisto( ESD_Conversion_Z_input1
								 "Conversions Z distribution","Z [cm]",StandardYAxis,
								 kTRUE, 1.5,0.00001,
								 kFALSE,0. ,0,
								 kTRUE, -201.,201.);
		}else{
				DrawAutoGammaHistos( ESD_Conversion_Z_input1, 
								 ESD_Conversion_Z_input2, 
								 "Conversions Z distribution","Z [cm] ",StandardYAxis,
								 kTRUE, 1.5,0.00001,
								 kFALSE,0. ,0.,
								 kTRUE, -201.,201.);
		}
		for(Int_t i=0; i < 12 ; i++){
			DrawGammaLines(ArrayZbins[i], ArrayZbins[i], 0.00001,ESD_Conversion_Z_input1->GetMaximum());
		}
		DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);	
	pad_Z->cd(2);
	pad_Z->cd(2)->SetLogy(0);
		if(secondinput == ""){
				DrawAutoGammaHisto( ESD_Conversion_Z_input1
								 "Conversions Z distribution","Z [cm]",StandardYAxis,
								 kTRUE, 1.5,0.00001,
								 kFALSE,0. ,0,
								 kTRUE, -201.,201.);
		}else{
				DrawAutoGammaHistos( ESD_Conversion_Z_input1, 
								 ESD_Conversion_Z_input2, 
								 "Conversions Z distribution","Z [cm] ",StandardYAxis,
								 kTRUE, 1.5,0.00001,
								 kFALSE,0. ,0.,
								 kTRUE, -201.,201.);
		}
		for(Int_t i=0; i < 12 ; i++){
			DrawGammaLines(ArrayZbins[i], ArrayZbins[i], 0.00001,ESD_Conversion_Z_input1->GetMaximum());
		}		
		DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);	

	pad_Z->Update();

	
	// --------------------------------page 4 - 2 Dimensional Plots-----------------------------------		
	
	//hier
	
	ps_mapping->NewPage();
	
	c_2dims = new TCanvas("c_2dims","",200,10,700,1000);  // gives the page size

	pad_2dims = new TPad("pad_2dims","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_2dims->SetFillColor(0);
	pad_2dims->GetFrame()->SetFillColor(0);
	pad_2dims->SetBorderMode(0);
	pad_2dims->Divide(2,2);
	pad_2dims->Draw();
	
	title0->Draw(); if(secondinput != ""){title1->Draw();}

	pad_2dims->cd(1)->SetRightMargin(RightMargin); 				
	pad_2dims->cd(1);
	pad_2dims->cd(1)->SetLogz(1);		
			DrawAutoGammaHisto2D(	ESD_Conversion_ZR_input1,
								"Conversion in ZR- input 1", "Z [cm]", "R [cm]", "Data",
								kTRUE, 0., 200.,
								kFALSE, 0., 20.);
			DrawAliceLogo(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3]);
	pad_2dims->cd(2)->SetRightMargin(RightMargin); 				           
	pad_2dims->cd(2);
	pad_2dims->cd(2)->SetLogz(1);	
			DrawAutoGammaHisto2D(	ESD_Conversion_XY_input1,
								"Conversion in XY- input 1", "X [cm]", "Y [cm]", "Data",
								kTRUE, -180., 180.,
								kTRUE, -180., 180.);
			DrawAliceLogo(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3]);
	
	if(secondinput != ""){
		pad_2dims->cd(3)->SetRightMargin(RightMargin); 				
		pad_2dims->cd(3);
		pad_2dims->cd(3)->SetLogz(1);	
			DrawAutoGammaHisto2D(	ESD_Conversion_ZR_input2,
								"Conversion in ZR- input 2", "Z [cm]", "R [cm]", "MC",
								kTRUE, 0., 200.,
								kFALSE, 0., 20.);
			DrawAliceLogo(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3]);
		pad_2dims->cd(4)->SetRightMargin(RightMargin); 				
		pad_2dims->cd(4);
		pad_2dims->cd(4)->SetLogz(1);	
			DrawAutoGammaHisto2D(	ESD_Conversion_XY_input2,
								"Conversion in XY- input 2", "X [cm]", "Y [cm]", "MC",
								kTRUE, -180., 180.,
								kTRUE, -180., 180.);
			DrawAliceLogo(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3]);
	}	
	
	pad_2dims->Update();


	// ----------------------------------- Page 5 - Phi in R -------------------------------------
	ps_mapping->NewPage();

	c_Phi_in_R = new TCanvas("c_Phi_in_R","",200,10,700,1000);  // gives the page size
	
	pad_Phi_in_R = new TPad("pad_Phi_in_R","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_Phi_in_R->SetFillColor(0);
	pad_Phi_in_R->GetFrame()->SetFillColor(0);
	pad_Phi_in_R->SetBorderMode(0);
	pad_Phi_in_R->Divide(column,raw);
	pad_Phi_in_R->Draw();
	
	title0->Draw();
	if(secondinput != ""){
		title1->Draw();
	}
	
	Double_t Integ_R_input1[NbinsR];	
	Double_t Integ_R_input2[NbinsR];			

	for(Int_t iR = 0; iR < NbinsR; iR++){
		Int_t place = iR + 1;
		
		pad_Phi_in_R->cd(place);
		pad_Phi_in_R->cd(place)->SetLogy(1);
		sprintf(histoname_Phi_in_R_input1,"ESD_Conversion_Mapping_Phi_in_R_iR%02d",iR);
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_Phi_in_R_input1[iR], 
								 histoname_Phi_in_R_input1,"#Phi",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos( ESD_Conversion_Mapping_Phi_in_R_input1[iR], 
								 ESD_Conversion_Mapping_Phi_in_R_input2[iR], 
								 histoname_Phi_in_R_input1,"#Phi",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
			Integ_R_input1[iR]=ESD_Conversion_Mapping_Phi_in_R_input1[iR]->Integral() ;			
			Integ_R_input2[iR]=ESD_Conversion_Mapping_Phi_in_R_input2[iR]->Integral() ;	
		}	
	}
	pad_Phi_in_R->Update();
	
	//------------------------------ Page 6 - Ratio Phi in R --------------------------------------
	if(secondinput != ""){	ps_mapping->NewPage();

	c_Ratio_Phi_in_R = new TCanvas("c_Ratio_Phi_in_R","",200,10,700,1000);  // gives the page size
	
	pad_Ratio_Phi_in_R = new TPad("pad_Ratio_Phi_in_R","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_Ratio_Phi_in_R->SetFillColor(0);
	pad_Ratio_Phi_in_R->GetFrame()->SetFillColor(0);
	pad_Ratio_Phi_in_R->SetBorderMode(0);
	pad_Ratio_Phi_in_R->Divide(column,raw);
	pad_Ratio_Phi_in_R->Draw();
	
	title0->Draw();
	if(secondinput != ""){
		title1->Draw();
	}

	for(Int_t iR = 0; iR < NbinsR; iR++){
		Int_t place = iR + 1;
		
		pad_Ratio_Phi_in_R->cd(place);
		pad_Ratio_Phi_in_R->cd(place)->SetLogy(1);
		sprintf(histoname_Ratio_Phi_in_R_input2,"Ratio_Phi_in_R_iR%02d",iR);
		if(secondinput != ""){		
			DrawRatioGammaHisto( 	Ratio_Mapping_Phi_in_R[iR], 
								 histoname_Ratio_Phi_in_R_input2,"#Phi","norm Data/norm MC",
								 kFALSE,3 ,0.000001,
								 kTRUE,0.1 ,5.,
								 kFALSE, 0.,180.);
			linePhi->Draw("same");
		}	
	}
	pad_Ratio_Phi_in_R->Update();
}


	//--------------- page 7 - Z in R ---------------------------------------------------
	
	ps_mapping->NewPage();
	
	c_Z_in_R = new TCanvas("c_Z_in_R","",200,10,700,1000);  // gives the page size
	
	pad_Z_in_R = new TPad("pad_Z_in_R","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_Z_in_R->SetFillColor(0);
	pad_Z_in_R->GetFrame()->SetFillColor(0);
	pad_Z_in_R->SetBorderMode(0);
	pad_Z_in_R->Divide(column,raw);
	pad_Z_in_R->Draw();
	
	title0->Draw(); 
	if(secondinput != ""){title1->Draw();}
	
	Double_t Integ_ZinR_input1[NbinsR];	
	Double_t Integ_ZinR_input2[NbinsR];			
	
	
	for(Int_t iR = 0; iR < NbinsR; iR++){
		Int_t place = iR + 1;
		pad_Z_in_R->cd(place);
		pad_Z_in_R->cd(place)->SetLogy(1);
		sprintf(histoname_Z_in_R_input1,"ESD_Conversion_Mapping_Z_in_R_iR%02d",iR);
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_Z_in_R_input1[iR], 
								 histoname_Z_in_R_input1,"Z [cm]",StandardYAxis,
								 kTRUE,3 ,0.00001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos( ESD_Conversion_Mapping_Z_in_R_input1[iR], 
								 ESD_Conversion_Mapping_Z_in_R_input2[iR], 
								 histoname_Z_in_R_input1,"Z [cm]",StandardYAxis,
								 kTRUE,3 ,0.00001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
			Integ_ZinR_input1[iR]=ESD_Conversion_Mapping_Z_in_R_input1[iR]->Integral() ;			
			Integ_ZinR_input2[iR]=ESD_Conversion_Mapping_Z_in_R_input2[iR]->Integral() ;	
		}	
	}
	
	pad_Z_in_R->Update();
	
	//------------------------------ Page 8 - Ratio Z in R --------------------------------------
	if(secondinput != ""){	ps_mapping->NewPage();

	c_Ratio_Z_in_R = new TCanvas("c_Ratio_Z_in_R","",200,10,700,1000);  // gives the page size
	
	pad_Ratio_Z_in_R = new TPad("pad_Ratio_Z_in_R","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_Ratio_Z_in_R->SetFillColor(0);
	pad_Ratio_Z_in_R->GetFrame()->SetFillColor(0);
	pad_Ratio_Z_in_R->SetBorderMode(0);
	pad_Ratio_Z_in_R->Divide(column,raw);
	pad_Ratio_Z_in_R->Draw();
	
	title0->Draw();
	if(secondinput != ""){
		title1->Draw();
	}

	for(Int_t iR = 0; iR < NbinsR; iR++){
		Int_t place = iR + 1;
		
		pad_Ratio_Z_in_R->cd(place);
		pad_Ratio_Z_in_R->cd(place)->SetLogy(1);
		sprintf(histoname_Ratio_Z_in_R_input2,"Ratio_Z_in_R_iR%02d",iR);
		if(secondinput != ""){		
			Float_t ZRange = RangeZinR[iR];
			DrawRatioGammaHisto( Ratio_Mapping_Z_in_R[iR], 
								 histoname_Ratio_Z_in_R_input2,"Z","norm Data/norm MC",
								 kFALSE,3 ,0.000001,
								 kTRUE,0.1 ,5.,
								 kTRUE, -ZRange,ZRange);
			DrawGammaLines(-ZRange,ZRange,1,1);
		}	
	}

	pad_Ratio_Z_in_R->Update();
}



	//---------page 9 - Phi in Z --------------------------------	
	
	ps_mapping->NewPage();
	
	
	c_Phi_in_Z = new TCanvas("c_Phi_in_Z","",200,10,700,1000);  // gives the page size
	
	pad_Phi_in_Z = new TPad("pad_Phi_in_Z","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_Phi_in_Z->SetFillColor(0);
	pad_Phi_in_Z->GetFrame()->SetFillColor(0);
	pad_Phi_in_Z->SetBorderMode(0);
	pad_Phi_in_Z->Divide(column,raw);
	pad_Phi_in_Z->Draw();
	
	title0->Draw(); if(secondinput != ""){title1->Draw();}
	
	Double_t Integ_Z_input1[NbinsZ];	
	Double_t Integ_Z_input2[NbinsZ];			
	
	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		Int_t place = iZ + 1;
		pad_Phi_in_Z->cd(place);
		pad_Phi_in_Z->cd(place)->SetLogy(1);
		sprintf(histoname_Phi_in_Z_input1,"ESD_Conversion_Mapping_Phi_in_Z_iZ%02d",iZ);
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_Phi_in_Z_input1[iZ], 
								 histoname_Phi_in_Z_input1,"#Phi",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos( ESD_Conversion_Mapping_Phi_in_Z_input1[iZ], 
								 ESD_Conversion_Mapping_Phi_in_Z_input2[iZ], 
								 histoname_Phi_in_Z_input1,"#Phi",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
			Integ_Z_input1[iZ]=ESD_Conversion_Mapping_Phi_in_Z_input1[iZ]->Integral() ;			
			Integ_Z_input2[iZ]=ESD_Conversion_Mapping_Phi_in_Z_input2[iZ]->Integral() ;	
		}	
	}
	
	pad_Phi_in_Z->Update();

	//------------------------------------- Page 10 - Ratio Phi in Z -----------------------------------	

	if(secondinput != ""){	ps_mapping->NewPage();

	c_Ratio_Phi_in_Z = new TCanvas("c_Ratio_Phi_in_Z","",200,10,700,1000);  // gives the page size
	
	pad_Ratio_Phi_in_Z = new TPad("pad_Ratio_Phi_in_Z","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_Ratio_Phi_in_Z->SetFillColor(0);
	pad_Ratio_Phi_in_Z->GetFrame()->SetFillColor(0);
	pad_Ratio_Phi_in_Z->SetBorderMode(0);
	pad_Ratio_Phi_in_Z->Divide(column,raw);
	pad_Ratio_Phi_in_Z->Draw();
	
	title0->Draw();
	if(secondinput != ""){
		title1->Draw();
	}

	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		Int_t place = iZ + 1;
		
		pad_Ratio_Phi_in_Z->cd(place);
		pad_Ratio_Phi_in_Z->cd(place)->SetLogy(1);
		sprintf(histoname_Ratio_Phi_in_Z_input2,"Ratio_Phi_in_Z_iZ%02d",iZ);
		if(secondinput != ""){		
			DrawRatioGammaHisto( Ratio_Mapping_Phi_in_Z[iZ], 
								 histoname_Ratio_Phi_in_Z_input2,"#Phi","norm Data/norm MC",
								 kFALSE,3 ,0.000001,
								 kTRUE,0.1 ,5.,
								 kFALSE, 0,0);
			linePhi->Draw("same");
		}	
	}

	pad_Ratio_Phi_in_Z->Update();
}

	//---------------page 11 - R in Z ---------------------------------------------------
	
	ps_mapping->NewPage();
	
	c_R_in_Z = new TCanvas("c_R_in_Z","",200,10,700,1000);  // gives the page size
	
	pad_R_in_Z = new TPad("pad_R_in_Z","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_R_in_Z->SetFillColor(0);
	pad_R_in_Z->GetFrame()->SetFillColor(0);
	pad_R_in_Z->SetBorderMode(0);
	pad_R_in_Z->Divide(column,raw);
	pad_R_in_Z->Draw();
	
	title0->Draw(); if(secondinput != ""){title1->Draw();}
	
	Double_t Integ_RinZ_input1[NbinsZ];	
	Double_t Integ_RinZ_input2[NbinsZ];			

	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		Int_t place = iZ + 1;
		pad_R_in_Z->cd(place);
		pad_R_in_Z->cd(place)->SetLogy(1);
		sprintf(histoname_R_in_Z_input1,"ESD_Conversion_Mapping_R_in_Z_iZ%02d",iZ);
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_R_in_Z_input1[iZ], 
								 histoname_R_in_Z_input1,"R [cm]",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos( ESD_Conversion_Mapping_R_in_Z_input1[iZ], 
								 ESD_Conversion_Mapping_R_in_Z_input2[iZ], 
								 histoname_R_in_Z_input1,"R [cm]",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
			Integ_RinZ_input1[iZ]=ESD_Conversion_Mapping_R_in_Z_input1[iZ]->Integral() ;			
			Integ_RinZ_input2[iZ]=ESD_Conversion_Mapping_R_in_Z_input2[iZ]->Integral() ;	
		}	
	}
	pad_R_in_Z->Update();
	
	//------------------------------ Page 12 - Ratio R in Z --------------------------------------
	if(secondinput != ""){	ps_mapping->NewPage();

	c_Ratio_R_in_Z = new TCanvas("c_Ratio_R_in_Z","",200,10,700,1000);  // gives the page size
	
	pad_Ratio_R_in_Z = new TPad("pad_Ratio_R_in_Z","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_Ratio_R_in_Z->SetFillColor(0);
	pad_Ratio_R_in_Z->GetFrame()->SetFillColor(0);
	pad_Ratio_R_in_Z->SetBorderMode(0);
	pad_Ratio_R_in_Z->Divide(column,raw);
	pad_Ratio_R_in_Z->Draw();
	
	title0->Draw();
	if(secondinput != ""){
		title1->Draw();
	}

	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		Int_t place = iZ + 1;
		
		pad_Ratio_R_in_Z->cd(place);
		pad_Ratio_R_in_Z->cd(place)->SetLogy(1);
		sprintf(histoname_Ratio_R_in_Z_input2,"Ratio_R_in_Z_iZ%02d",iZ);
		if(secondinput != ""){		
			DrawRatioGammaHisto( Ratio_Mapping_R_in_Z[iZ], 
								 histoname_Ratio_R_in_Z_input2,"R [cm]","norm Data/norm MC",
								 kFALSE,3 ,0.000001,
								 kTRUE,0.1 ,5.,
								 kFALSE, 0,0);
			lineR->Draw("same");
		}	
	}
	
	pad_Ratio_R_in_Z->Update();
}	

	// --------------------page 13 - MidPt Phi in R-------------------------------------------------------------
	
	
	ps_mapping->NewPage();
	
	
	c_MidPt_Phi_in_R = new TCanvas("c_MidPt_Phi_in_R","",200,10,700,1000);  // gives the page size

	pad_MidPt_Phi_in_R = new TPad("pad_MidPt_Phi_in_R","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_MidPt_Phi_in_R->SetFillColor(0);
	pad_MidPt_Phi_in_R->GetFrame()->SetFillColor(0);
	pad_MidPt_Phi_in_R->SetBorderMode(0);
	pad_MidPt_Phi_in_R->Divide(column,raw);
	pad_MidPt_Phi_in_R->Draw();
	
	title0->Draw(); if(secondinput != ""){title1->Draw();}

	Double_t Integ_midpt_R_input1[NbinsR];	
	Double_t Integ_midpt_R_input2[NbinsR];			
	
	for(Int_t iR = 0; iR < NbinsR; iR++){
		Int_t place = iR + 1;
		pad_MidPt_Phi_in_R->cd(place);
		pad_MidPt_Phi_in_R->cd(place)->SetLogy(1);
		sprintf(histoname_MidPt_Phi_in_R_input1,"ESD_Conversion_Mapping_MidPt_Phi_in_R_iR%02d",iR);
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_MidPt_Phi_in_R_input1[iR], 
								 histoname_MidPt_Phi_in_R_input1,"#Phi",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos(ESD_Conversion_Mapping_MidPt_Phi_in_R_input1[iR], 
							ESD_Conversion_Mapping_MidPt_Phi_in_R_input2[iR], 
								 histoname_MidPt_Phi_in_Z_input1,"#Phi",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
			Integ_midpt_R_input1[iR]=ESD_Conversion_Mapping_MidPt_Phi_in_R_input1[iR]->Integral() ;			
			Integ_midpt_R_input2[iR]=ESD_Conversion_Mapping_MidPt_Phi_in_R_input2[iR]->Integral() ;	
		}	
	}
	
	pad_MidPt_Phi_in_R->Update();
	
	//---------------------------page 14 - MidPt Z in R ---------------------------------------
	
	ps_mapping->NewPage();
	
	c_midpt_Z_in_R = new TCanvas("c_midpt_Z_in_R","",200,10,700,1000);  // gives the page size
	
	pad_midpt_Z_in_R = new TPad("pad_midpt_Z_in_R","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_midpt_Z_in_R->SetFillColor(0);
	pad_midpt_Z_in_R->GetFrame()->SetFillColor(0);
	pad_midpt_Z_in_R->SetBorderMode(0);
	pad_midpt_Z_in_R->Divide(column,raw);
	pad_midpt_Z_in_R->Draw();
	
	title0->Draw(); if(secondinput != ""){title1->Draw();}
	
	Double_t Integ_midpt_ZinR_input1[NbinsR];	
	Double_t Integ_midpt_ZinR_input2[NbinsR];			
	
	

	for(Int_t iR = 0; iR < NbinsR; iR++){
		Int_t place = iR + 1;
		pad_midpt_Z_in_R->cd(place);
		pad_midpt_Z_in_R->cd(place)->SetLogy(1);
		sprintf(histoname_MidPt_Z_in_R_input1,"ESD_Conversion_Mapping_MidPt_Z_in_R_iR%02d",iR);
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_MidPt_Z_in_R_input1[iR], 
								 histoname_MidPt_Z_in_R_input1,"Z [cm]",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos(ESD_Conversion_Mapping_MidPt_Z_in_R_input1[iR], 
							ESD_Conversion_Mapping_MidPt_Z_in_R_input2[iR], 
								 histoname_MidPt_Z_in_R_input1,"Z [cm]",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
			Integ_midpt_ZinR_input1[iR]=ESD_Conversion_Mapping_MidPt_Z_in_R_input1[iR]->Integral() ;			
			Integ_midpt_ZinR_input2[iR]=ESD_Conversion_Mapping_MidPt_Z_in_R_input2[iR]->Integral() ;	
		}	
	}

	pad_midpt_Z_in_R->Update();

	
	//---------------------page 15 - MidPt Phi in Z--------------------	
	
	if(!SinglePlots)ps_mapping->NewPage();
	
	
	c_MidPt_Phi_in_Z = new TCanvas("c_MidPt_Phi_in_Z","",200,10,700,1000);  // gives the page size
	
	pad_MidPt_Phi_in_Z = new TPad("pad_MidPt_Phi_in_Z","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_MidPt_Phi_in_Z->SetFillColor(0);
	pad_MidPt_Phi_in_Z->GetFrame()->SetFillColor(0);
	pad_MidPt_Phi_in_Z->SetBorderMode(0);
	pad_MidPt_Phi_in_Z->Divide(column,raw);
	pad_MidPt_Phi_in_Z->Draw();
	
	title0->Draw(); if(secondinput != ""){title1->Draw();}
	
	Double_t Integ_midpt_Z_input1[NbinsZ];	
	Double_t Integ_midpt_Z_input2[NbinsZ];			
	
	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		Int_t place = iZ + 1;
		pad_MidPt_Phi_in_Z->cd(place);
		pad_MidPt_Phi_in_Z->cd(place)->SetLogy(1);
		sprintf(histoname_MidPt_Phi_in_Z_input1,"ESD_Conversion_Mapping_MidPt_Phi_in_Z_iR%02d",iR);
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_MidPt_Phi_in_Z_input1[iZ], 
								 histoname_MidPt_Phi_in_Z_input1,"#Phi",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos(ESD_Conversion_Mapping_MidPt_Phi_in_Z_input1[iZ], 
							ESD_Conversion_Mapping_MidPt_Phi_in_Z_input2[iZ], 
								 histoname_MidPt_Phi_in_Z_input1,"#Phi",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
			Integ_midpt_Z_input1[iZ]=ESD_Conversion_Mapping_MidPt_Phi_in_Z_input1[iZ]->Integral() ;			
			Integ_midpt_Z_input2[iZ]=ESD_Conversion_Mapping_MidPt_Phi_in_Z_input2[iZ]->Integral() ;	
		}	
	}
	pad_MidPt_Phi_in_Z->Update();
	
	//---------------------------page 16 - MidPt R in Z ---------------------------------------
	
	if(!SinglePlots)ps_mapping->NewPage();
	
	c_MidPt_R_in_Z = new TCanvas("c_MidPt_R_in_Z","",200,10,700,1000);  // gives the page size
	
	pad_MidPt_R_in_Z = new TPad("pad_MidPt_R_in_Z","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
	pad_MidPt_R_in_Z->SetFillColor(0);
	pad_MidPt_R_in_Z->GetFrame()->SetFillColor(0);
	pad_MidPt_R_in_Z->SetBorderMode(0);
	pad_MidPt_R_in_Z->Divide(column,raw);
	pad_MidPt_R_in_Z->Draw();
	
	title0->Draw(); if(secondinput != ""){title1->Draw();}
	
	Double_t Integ_midpt_RinZ_input1[NbinsZ];	
	Double_t Integ_midpt_RinZ_input2[NbinsZ];			

	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		Int_t place = iZ + 1;
		pad_MidPt_R_in_Z->cd(place);
		pad_MidPt_R_in_Z->cd(place)->SetLogy(1);
		sprintf(histoname_MidPt_R_in_Z_input1,"ESD_Conversion_Mapping_MidPt_R_in_Z_iR%02d",iR);
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_MidPt_R_in_Z_input1[iZ], 
								 histoname_MidPt_R_in_Z_input1,"R [cm]",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos(ESD_Conversion_Mapping_MidPt_R_in_Z_input1[iZ], 
							ESD_Conversion_Mapping_MidPt_R_in_Z_input2[iZ], 
								 histoname_MidPt_R_in_Z_input1,"R [cm]",StandardYAxis,
								 kTRUE,3 ,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
			Integ_midpt_RinZ_input1[iZ]=ESD_Conversion_Mapping_MidPt_R_in_Z_input1[iZ]->Integral() ;			
			Integ_midpt_RinZ_input2[iZ]=ESD_Conversion_Mapping_MidPt_R_in_Z_input2[iZ]->Integral() ;	
		}	
	}
	pad_MidPt_R_in_Z->Update();
	
	   if(secondinput != ""){
		Double_t Sum_phi_in_R_input1 = 0., Sum_phi_in_Z_input1 =0., Sum_midpt_phi_in_R_input1 = 0., Sum_midpt_phi_in_Z_input1= 0.;
		Double_t Sum_phi_in_R_input2 = 0., Sum_phi_in_Z_input2 =0., Sum_midpt_phi_in_R_input2 = 0., Sum_midpt_phi_in_Z_input2= 0.;
		Double_t Sum_Z_in_R_input1 = 0., Sum_R_in_Z_input1 =0., Sum_midpt_Z_in_R_input1 = 0., Sum_midpt_R_in_Z_input1= 0.;
		Double_t Sum_Z_in_R_input2 = 0., Sum_R_in_Z_input2 =0., Sum_midpt_Z_in_R_input2 = 0., Sum_midpt_R_in_Z_input2= 0.;
		
		for(Int_t iR = 0; iR < NbinsR; iR++){
			Sum_phi_in_R_input1 = Sum_phi_in_R_input1 + Integ_R_input1[iR];
			Sum_phi_in_R_input2 = Sum_phi_in_R_input2 + Integ_R_input2[iR];
			Sum_midpt_phi_in_R_input1 = Sum_midpt_phi_in_R_input1 + Integ_midpt_R_input1[iR];
			Sum_midpt_phi_in_R_input2 = Sum_midpt_phi_in_R_input2 + Integ_midpt_R_input2[iR];
			Sum_Z_in_R_input1 = Sum_Z_in_R_input1 + Integ_ZinR_input1[iR];
			Sum_Z_in_R_input2 = Sum_Z_in_R_input2 + Integ_ZinR_input2[iR];
			Sum_midpt_Z_in_R_input1 = Sum_midpt_Z_in_R_input1 + Integ_midpt_ZinR_input1[iR];
			Sum_midpt_Z_in_R_input2 = Sum_midpt_Z_in_R_input2 + Integ_midpt_ZinR_input2[iR];
				
		}
		for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
			Sum_phi_in_Z_input1 = Sum_phi_in_Z_input1 + Integ_Z_input1[iZ];
			Sum_phi_in_Z_input2 = Sum_phi_in_Z_input2 + Integ_Z_input2[iZ];
			Sum_midpt_phi_in_Z_input1 = Sum_midpt_phi_in_Z_input1 + Integ_midpt_Z_input1[iZ];
			Sum_midpt_phi_in_Z_input2 = Sum_midpt_phi_in_Z_input2 + Integ_midpt_Z_input2[iZ];
			Sum_R_in_Z_input1 = Sum_R_in_Z_input1 + Integ_RinZ_input1[iZ];
			Sum_R_in_Z_input2 = Sum_R_in_Z_input2 + Integ_RinZ_input2[iZ];
			Sum_midpt_R_in_Z_input1 = Sum_midpt_R_in_Z_input1 + Integ_midpt_RinZ_input1[iZ];
			Sum_midpt_R_in_Z_input2 = Sum_midpt_R_in_Z_input2 + Integ_midpt_RinZ_input2[iZ];
		}
					
		const char *outlname = "Export_Mapping.dat";
		fstream outl;
    		outl.open(outlname, ios::out);
		outl << "#Calculating Integrals" << endl;	
		outl << "------------------------------------------------------------------------------------------" << endl;
		outl << "# This file is created to display the Integrals of the different bins of the different diagrams. The forth column displays the Integral over the input1 in that bin divided by the sum over all bins of input1, the same does the sixth column for input2. The eigth column displays the difference of input1 - input2 divided by input2 of that bin. Therefore you should put the data in the first input and the Montecarlo in the second." << endl;
		outl << "------------------------------------------------------------------------------------------" << endl;

		outl << "input1 :\t" << input1  << endl;
		outl << "input2 :\t" << secondinput << endl;
		outl << "------------------------------------------------------------------------------------------" << endl;
		outl << endl;
		outl << "\t input1 \t input2 " << endl;
		outl << "Number of events" << "\t" << nGoodEvents_input1 << "\t" << nGoodEvents_input2 << endl;
		outl << "Number of triggers" << "\t" << nGoodTrig_input1 << "\t" << nGoodTrig_input2 << endl;
		outl << "Number reconstructed gammas"<< "\t" << nRecGamma_input1 << "\t" << nRecGamma_input2 <<endl;
		outl << "Mean Multiplicity" << "\t" <<  mean_input1 << "\t" << mean_input2 << endl;
		outl << endl;
		outl << endl;

		outl << "------------------------------------------------------------------------------------------" << endl;		
		outl << "graph \t bin \t Input1 \t % \t Input2 \t % \t Input1-Input2 \t %" <<endl;


		outl << "Phi in R" <<endl;
		outl <<"\t0\t" << Integ_R_input1[0] << "\t" << Integ_R_input1[0]/Sum_phi_in_R_input1 * 100 <<"\t" << Integ_R_input2[0] << "\t" << Integ_R_input2[0]/Sum_phi_in_R_input2 * 100 << "\t" << Integ_R_input1[0]- Integ_R_input2[0] << "\t" << (Integ_R_input1[0]- Integ_R_input2[0])/Integ_R_input2[0] *100 << endl; 	
		for(Int_t iR = 1; iR < NbinsR; iR++){
			outl << "\t" << iR << "\t" << Integ_R_input1[iR] << "\t" << Integ_R_input1[iR]/Sum_phi_in_R_input1 * 100 <<"\t" << Integ_R_input2[iR] << "\t" << Integ_R_input2[iR]/Sum_phi_in_R_input2 * 100 << "\t" << Integ_R_input1[iR]- Integ_R_input2[iR] << "\t" << (Integ_R_input1[iR]- Integ_R_input2[iR])/Integ_R_input2[iR]* 100 << endl;
		}
		outl << "\t" << "sum" << "\t" << ESD_Conversion_R_input1->Integral() << "\t" << "100" << "\t" << ESD_Conversion_R_input2->Integral() << "\t" << "100" <<"\t" << ESD_Conversion_R_input1->Integral() - ESD_Conversion_R_input2->Integral() << "\t" << (ESD_Conversion_R_input1->Integral()-ESD_Conversion_R_input2->Integral())/ESD_Conversion_R_input2->Integral() *100 << endl;
		outl << endl;

		outl << "Phi in Z" <<endl;
//		outl <<"\t0\t" << Integ_Z_input1[0] << "\t"<< Integ_Z_input1[0]/Sum_phi_in_Z_input1 * 100 <<"\t" << Integ_Z_input2[0] << "\t" << Integ_Z_input2[0]/Sum_phi_in_Z_input2 * 100 << "\t" << Integ_Z_input1[0]- Integ_Z_input2[0] << "\t" << (Integ_Z_input1[0]- Integ_Z_input2[0])/Integ_Z_input2[0] * 100 << endl;
		for(Int_t iZ = 1; iZ < NbinsZ; iZ++){
			outl << "\t" << iZ<< "\t" << Integ_Z_input1[iZ] << "\t" << Integ_Z_input1[iZ]/Sum_phi_in_Z_input1 * 100 <<"\t" << Integ_Z_input2[iZ] << "\t" << Integ_Z_input2[iZ]/Sum_phi_in_Z_input2 * 100 << "\t" << Integ_Z_input1[iZ]- Integ_Z_input2[iZ] << "\t" << (Integ_Z_input1[iZ]- Integ_Z_input2[iZ])/Integ_Z_input2[iZ]*100 << endl;
		}
		outl << "\t" << "sum" << "\t" << ESD_Conversion_Z_input1->Integral() << "\t" << "100" << "\t" << ESD_Conversion_Z_input2->Integral() << "\t" << "100" <<"\t" << ESD_Conversion_Z_input1->Integral() - ESD_Conversion_Z_input2->Integral() << "\t" << (ESD_Conversion_Z_input1->Integral()-ESD_Conversion_Z_input2->Integral())/ESD_Conversion_Z_input2->Integral() *100 << endl;
		outl << endl;
	
		outl << "MidPt Phi in R" <<endl;
		outl <<" \t0\t" << Integ_midpt_R_input1[0] << "\t" << Integ_midpt_R_input1[0]/Sum_midpt_phi_in_R_input1 * 100 <<"\t" << Integ_midpt_R_input2[0] << "\t" << Integ_midpt_R_input2[0]/Sum_midpt_phi_in_R_input2 * 100 << "\t" << Integ_midpt_R_input1[0]- Integ_midpt_R_input2[0] << "\t" << (Integ_midpt_R_input1[0]- Integ_midpt_R_input2[0])/Integ_midpt_R_input2[0]*100 << endl; 	
		for(Int_t iR = 1; iR < NbinsR; iR++){
			outl << "\t" << iR << "\t" << Integ_midpt_R_input1[iR] << "\t" << Integ_midpt_R_input1[iR]/Sum_midpt_phi_in_R_input1 * 100 <<"\t" << Integ_midpt_R_input2[iR] << "\t" << Integ_midpt_R_input2[iR]/Sum_midpt_phi_in_R_input2 * 100 << "\t" << Integ_midpt_R_input1[iR]- Integ_midpt_R_input2[iR] << "\t" << (Integ_midpt_R_input1[iR]- Integ_midpt_R_input2[iR])/Integ_midpt_R_input2[iR] *100 << endl;
		}
		outl << "\t" << "sum" << "\t" << Sum_midpt_phi_in_R_input1 << "\t" << "100" << "\t" << Sum_midpt_phi_in_R_input2 << "\t" << "100" <<"\t" << Sum_midpt_phi_in_R_input1 - Sum_midpt_phi_in_R_input2 << "\t" << (Sum_midpt_phi_in_R_input1-Sum_midpt_phi_in_R_input2)/Sum_midpt_phi_in_R_input2 *100 << endl;
	
		outl << endl;

		outl << "MidPt Phi in Z" <<endl;
//		outl <<" \t0\t" << Integ_midpt_Z_input1[0] << "\t"<< Integ_midpt_Z_input1[0]/Sum_midpt_phi_in_Z_input1 * 100 <<"\t" << Integ_midpt_Z_input2[0] << "\t" << Integ_midpt_Z_input2[0]/Sum_midpt_phi_in_Z_input2 * 100 << "\t" << Integ_midpt_Z_input1[0]- Integ_midpt_Z_input2[0] << "\t" << (Integ_midpt_Z_input1[0]- Integ_midpt_Z_input2[0])/Integ_midpt_Z_input2[0] *100<< endl;
		for(Int_t iZ = 1; iZ < NbinsZ; iZ++){
			outl << "\t" << iZ << "\t" << Integ_midpt_Z_input1[iZ] << "\t" << Integ_midpt_Z_input1[iZ]/Sum_midpt_phi_in_Z_input1 * 100 <<"\t" << Integ_midpt_Z_input2[iZ] << "\t" << Integ_midpt_Z_input2[iZ]/Sum_midpt_phi_in_Z_input2 * 100 << "\t" << Integ_midpt_Z_input1[iZ]- Integ_midpt_Z_input2[iZ] << "\t" << (Integ_midpt_Z_input1[iZ]- Integ_midpt_Z_input2[iZ])/Integ_midpt_Z_input2[iZ] *100 << endl;
		}
		outl << "\t" << "sum" << "\t" << Sum_midpt_phi_in_Z_input1 << "\t" << "100" << "\t" << Sum_midpt_phi_in_Z_input2 << "\t" << "100" <<"\t" << Sum_midpt_phi_in_Z_input1 - Sum_midpt_phi_in_Z_input2 << "\t" << (Sum_midpt_phi_in_Z_input1-Sum_midpt_phi_in_Z_input2)/Sum_midpt_phi_in_Z_input2 *100 << endl;
		outl << endl;

		outl << "Z in R" <<endl;
		outl <<" \t0\t" << Integ_ZinR_input1[0] << "\t" << Integ_R_input1[0]/Sum_Z_in_R_input1 * 100 <<"\t" << Integ_ZinR_input2[0] << "\t" << Integ_ZinR_input2[0]/Sum_Z_in_R_input2 * 100 << "\t" << Integ_ZinR_input1[0]- Integ_ZinR_input2[0] << "\t" << (Integ_ZinR_input1[0]- Integ_ZinR_input2[0])/Integ_ZinR_input2[0] *100 << endl; 	
		for(Int_t iR = 1; iR < NbinsR; iR++){
			outl << "\t" << iR << "\t" << Integ_ZinR_input1[iR] << "\t" << Integ_ZinR_input1[iR]/Sum_Z_in_R_input1 * 100 <<"\t" << Integ_ZinR_input2[iR] << "\t" << Integ_ZinR_input2[iR]/Sum_Z_in_R_input2 * 100 << "\t" << Integ_ZinR_input1[iR]- Integ_ZinR_input2[iR] << "\t" << (Integ_ZinR_input1[iR]- Integ_ZinR_input2[iR])/Integ_ZinR_input2[iR]*100 << endl;
		}
		outl << "\t" << "sum" << "\t" << ESD_Conversion_R_input1->Integral() << "\t" << "100" << "\t" << ESD_Conversion_R_input2->Integral() << "\t" << "100" <<"\t" << ESD_Conversion_R_input1->Integral() - ESD_Conversion_R_input2->Integral() << "\t" << (ESD_Conversion_R_input1->Integral()-ESD_Conversion_R_input2->Integral())/ESD_Conversion_R_input2->Integral() *100 << endl;
		outl << endl;

		outl << "R in Z" <<endl;
//		outl <<" \t0\t" << Integ_RinZ_input1[0] << "\t"<< Integ_RinZ_input1[0]/Sum_R_in_Z_input1 * 100 <<"\t" << Integ_RinZ_input2[0] << "\t" << Integ_RinZ_input2[0]/Sum_R_in_Z_input2 * 100 << "\t" << Integ_RinZ_input1[0]- Integ_RinZ_input2[0] << "\t" << (Integ_RinZ_input1[0]- Integ_RinZ_input2[0])/Integ_RinZ_input2[0] * 100 << endl;
		for(Int_t iZ = 1; iZ < NbinsZ; iZ++){
			outl << "\t" << iZ<< "\t" << Integ_RinZ_input1[iZ] << "\t" << Integ_RinZ_input1[iZ]/Sum_R_in_Z_input1 * 100 <<"\t" << Integ_RinZ_input2[iZ] << "\t" << Integ_RinZ_input2[iZ]/Sum_R_in_Z_input2 * 100 << "\t" << Integ_RinZ_input1[iZ]- Integ_RinZ_input2[iZ] << "\t" << (Integ_RinZ_input1[iZ]- Integ_RinZ_input2[iZ])/Integ_RinZ_input2[iZ]*100 << endl;
		}
		outl << "\t" << "sum" << "\t" << ESD_Conversion_Z_input1->Integral() << "\t" << "100" << "\t" << ESD_Conversion_Z_input2->Integral() << "\t" << "100" <<"\t" << ESD_Conversion_Z_input1->Integral() - ESD_Conversion_Z_input2->Integral() << "\t" << (ESD_Conversion_Z_input1->Integral()-ESD_Conversion_Z_input2->Integral())/ESD_Conversion_Z_input2->Integral() *100 << endl;
		outl << endl;

		outl << "MidPt Z in R" <<endl;
		outl <<" \t0\t" << Integ_midpt_ZinR_input1[0] << "\t" << Integ_midpt_ZinR_input1[0]/Sum_midpt_Z_in_R_input1 * 100 <<"\t" << Integ_midpt_ZinR_input2[0] << "\t" << Integ_midpt_ZinR_input2[0]/Sum_midpt_Z_in_R_input2 * 100 << "\t" << Integ_midpt_ZinR_input1[0]- Integ_midpt_ZinR_input2[0] << "\t" << (Integ_midpt_ZinR_input1[0]- Integ_midpt_ZinR_input2[0])/Integ_midpt_ZinR_input2[0]*100 << endl; 	
		for(Int_t iR = 1; iR < NbinsR; iR++){
			outl << "\t" << iR << "\t" << Integ_midpt_ZinR_input1[iR] << "\t" << Integ_midpt_ZinR_input1[iR]/Sum_midpt_Z_in_R_input1 * 100 <<"\t" << Integ_midpt_ZinR_input2[iR] << "\t" << Integ_midpt_ZinR_input2[iR]/Sum_midpt_Z_in_R_input2 * 100 << "\t" << Integ_midpt_ZinR_input1[iR]- Integ_midpt_ZinR_input2[iR] << "\t" << (Integ_midpt_ZinR_input1[iR]- Integ_midpt_ZinR_input2[iR])/Integ_midpt_ZinR_input2[iR] *100 << endl;
		}
		outl << "\t" << "sum" << "\t" << Sum_midpt_Z_in_R_input1 << "\t" << "100" << "\t" << Sum_midpt_Z_in_R_input2 << "\t" << "100" <<"\t" << Sum_midpt_Z_in_R_input1 - Sum_midpt_Z_in_R_input2 << "\t" << (Sum_midpt_Z_in_R_input1-Sum_midpt_Z_in_R_input2)/Sum_midpt_Z_in_R_input2 *100 << endl;
		outl << endl;
		
		outl << "MidPt R in Z" <<endl;
//		outl <<" \t0\t" << Integ_midpt_RinZ_input1[0] << "\t"<< Integ_midpt_RinZ_input1[0]/Sum_midpt_R_in_Z_input1 * 100 <<"\t" << Integ_midpt_RinZ_input2[0] << "\t" << Integ_midpt_RinZ_input2[0]/Sum_midpt_R_in_Z_input2 * 100 << "\t" << Integ_midpt_RinZ_input1[0]- Integ_midpt_RinZ_input2[0] << "\t" << (Integ_midpt_RinZ_input1[0]- Integ_midpt_RinZ_input2[0])/Integ_midpt_RinZ_input2[0] *100<< endl;
		for(Int_t iZ = 1; iZ < NbinsZ; iZ++){
			outl << "\t" << iZ << "\t" << Integ_midpt_RinZ_input1[iZ] << "\t" << Integ_midpt_RinZ_input1[iZ]/Sum_midpt_R_in_Z_input1 * 100 <<"\t" << Integ_midpt_RinZ_input2[iZ] << "\t" << Integ_midpt_RinZ_input2[iZ]/Sum_midpt_R_in_Z_input2 * 100 << "\t" << Integ_midpt_RinZ_input1[iZ]- Integ_midpt_RinZ_input2[iZ] << "\t" << (Integ_midpt_RinZ_input1[iZ]- Integ_midpt_RinZ_input2[iZ])/Integ_midpt_RinZ_input2[iZ] *100 << endl;
		}
		outl << "\t" << "sum" << "\t" << Sum_midpt_R_in_Z_input1 << "\t" << "100" << "\t" << Sum_midpt_R_in_Z_input2 << "\t" << "100" <<"\t" << Sum_midpt_R_in_Z_input1 - Sum_midpt_R_in_Z_input2 << "\t" << (Sum_midpt_R_in_Z_input1-Sum_midpt_R_in_Z_input2)/Sum_midpt_R_in_Z_input2 *100 << endl;

		outl << "------------------------------------------------------------------------------------------" << endl;		
	   }
	ps_mapping->Close();
	// hier ende
	}	
	
	// --------------------------- single plots  ----------------------------
	if(SinglePlots){
	

		// ---------------------------- R- Distribution -----------------------
		TCanvas * c_SinglePlot_2 = new TCanvas("c_SinglePlot_2","",1000,1000);  // gives the page size	
		c_SinglePlot_2->SetLogy(0);		
		c_SinglePlot_2->cd();

		if(secondinput == ""){
				DrawAutoGammaHisto( ESD_Conversion_R_input1, 
								 "","R [cm]",StandardYAxis,
								 kTRUE, 1.1,0.00001,
								 kFALSE,0. ,0,
								 kTRUE, 0.,180.);
		}else{
				DrawAutoGammaHistos( ESD_Conversion_R_input1, 
								 ESD_Conversion_R_input2, 
								 "","R [cm]",StandardYAxis,
								 kTRUE, 1.1,0.00001,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,180.);
		}
		for(Int_t i=0; i < 12 ; i++){
			DrawGammaLines(ArrayRbins[i], ArrayRbins[i], 0.00001,ESD_Conversion_R_input1->GetMaximum());
		}		
		DrawAliceLogoPerformance(right_up [0],right_up[1],right_up[2],right_up[3],0.03, Date);	
		c_SinglePlot_2->Update();
		c_SinglePlot_2->SaveAs(Form("%sR_distribution_lin.%s",path,suffix));
		delete c_SinglePlot_2;		


		TCanvas * c_SinglePlot_1 = new TCanvas("c_SinglePlot_1","",1000,1000);  // gives the page size
		c_SinglePlot_1->cd();
		c_SinglePlot_1->SetLogy(1);

		if(secondinput == ""){
				DrawAutoGammaHisto( ESD_Conversion_R_input1, 
								 "","R [cm]",StandardYAxis,
								 kTRUE, 1.5,0.00001,
								 kFALSE,0. ,0,
								 kTRUE, 0.,180.);
		}else{
				DrawAutoGammaHistos( ESD_Conversion_R_input1, 
								 ESD_Conversion_R_input2, 
								 "","R [cm]",StandardYAxis,
								 kTRUE, 1.5,0.00001,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,180.);
		}
		for(Int_t i=0; i < 12 ; i++){
			DrawGammaLines(ArrayRbins[i], ArrayRbins[i], 0.00001,ESD_Conversion_R_input1->GetMaximum());
		}		
		DrawAliceLogoPerformance(right_up [0],right_up[1],right_up[2],right_up[3],0.03, Date);	

		c_SinglePlot_1->Update();
		c_SinglePlot_1->SaveAs(Form("%sR_distribution_log.%s",path,suffix));
		delete c_SinglePlot_1;


		//----------------------------- Integrated Radius -----------------------
		TCanvas * c_SinglePlot_18 = new TCanvas("c_SinglePlot_18","",1000,1000);  // gives the page size		
		c_SinglePlot_18->SetLogy(0);
		c_SinglePlot_18->cd();
		if(secondinput == ""){
				DrawAutoGammaHisto( ESD_Conversion_R_summed_input1, 
								 "Conversions R distribution","R [cm]","Integrated (0-to-R) #gamma candidates/event scaled by multiplicity",
								 kTRUE, 1.2,0,
								 kTRUE,0. ,0.,
								 kTRUE, 0.,180.);
		}else{
				DrawAutoGammaHistos( ESD_Conversion_R_summed_input1, 
								 ESD_Conversion_R_summed_input2, 
								 "Integrated Radius of Conversion","R [cm] ","Integrated (0-to-R) #gamma candidates/event scaled by multiplicity",
								 kTRUE, 1.2,0,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,180.);
		}
		DrawAliceLogo(left_up [0],left_up[1],left_up[2],left_up[3]);	

		c_SinglePlot_18->Update();
		c_SinglePlot_18->SaveAs(Form("%sInteg_Radius.%s",path,suffix));
		delete c_SinglePlot_18;		


		TCanvas * c_SinglePlot_17 = new TCanvas("c_SinglePlot_17","",1000,1000);  // gives the page size
		c_SinglePlot_17->SetLogy(1);
		c_SinglePlot_17->cd();
		if(secondinput == ""){
				DrawAutoGammaHisto( ESD_Conversion_R_summed_input1, 
								 "Conversions R distribution","R [cm]","Integrated (0-to-R) #gamma candidates/event scaled by multiplicity",
								 kTRUE,2.5,0,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,180.);
		}else{
				DrawAutoGammaHistos( ESD_Conversion_R_summed_input1, 
								 ESD_Conversion_R_summed_input2, 
								 "Integrated Radius of Conversion","R [cm] ","Integrated (0-to-R) #gamma candidates/event scaled by multiplicity",
								 kTRUE, 2.5,0,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,180.);
		}
		DrawAliceLogo(right_down [0],right_down[1],right_down[2],right_down[3]);	

			
		c_SinglePlot_17->Update();
		c_SinglePlot_17->SaveAs(Form("%sInteg_Radius_log.%s",path,suffix));
		delete c_SinglePlot_17;



		// -----------------------  2 dim Plots -----------------------------------
		TCanvas * c_SinglePlot_3 = new TCanvas("c_SinglePlot_3","",1000,1000);  // gives the page size
		c_SinglePlot_3->SetLogz(1);	
		c_SinglePlot_3->SetRightMargin(RightMargin); 						
		c_SinglePlot_3->cd();
			DrawAutoGammaHisto2D(	ESD_Conversion_ZR_input1,
								"", "Z [cm]", "R [cm]", "",
								kTRUE, 0., 200.,
								kFALSE, 0., 20.);
			DrawAliceLogoPerformance(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3],0.03,Date);
//			DrawAliceText(left_down[0], left_down[1 ], left_down[3]);
			EtaRange->Draw();
		c_SinglePlot_3->Update();
		c_SinglePlot_3->SaveAs(Form("%sZR_distribution.%s",path,suffix));
		delete c_SinglePlot_3;

		TCanvas * c_SinglePlot_4 = new TCanvas("c_SinglePlot_4","",1000,1000);  // gives the page size
		c_SinglePlot_4->SetLogz(1);
		c_SinglePlot_4->SetRightMargin(RightMargin); 						
		c_SinglePlot_4->cd();
			DrawAutoGammaHisto2D(	ESD_Conversion_XY_input1,
								"", "X [cm]", "Y [cm]", "",
								kTRUE, -180., 180.,
								kTRUE, -180., 180.);
			DrawStructure();
			EtaRange->Draw();
			DrawAliceLogoPerformance(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3],0.03,Date);
//			DrawAliceText(left_down[0], left_down[1 ], left_down[3]);
		c_SinglePlot_4->Update();
		c_SinglePlot_4->SaveAs(Form("%sXY_distribution.%s",path,suffix));
		delete c_SinglePlot_4;

		TCanvas * c_SinglePlot_4 = new TCanvas("c_SinglePlot_4","",1000,1000);  // gives the page size
		c_SinglePlot_4->SetLogz(1);
		c_SinglePlot_4->SetRightMargin(RightMargin); 						
		c_SinglePlot_4->cd();
			ESD_Conversion_XY_input1->SetMinimum(minimumXY);
			DrawAutoGammaHisto2D(	ESD_Conversion_XY_input1,
								"", "X [cm]", "Y [cm]", "",
								kTRUE, -180., 180.,
								kTRUE, -180., 180.);
			DrawStructure();
			EtaRange->Draw();
			DrawAliceLogoPerformance(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3],0.03,Date);
//			DrawAliceText(left_down[0], left_down[1 ], left_down[3]);
		c_SinglePlot_4->Update();
		c_SinglePlot_4->SaveAs(Form("%sXY_distribution_minimum.%s",path,suffix));
		delete c_SinglePlot_4;


		
		if(secondinput != ""){
			TCanvas * c_SinglePlot_5 = new TCanvas("c_SinglePlot_5","",1000,1000);  // gives the page size
			c_SinglePlot_5->SetLogz(1);
			c_SinglePlot_5->SetRightMargin(RightMargin); 						
			c_SinglePlot_5->cd();
			Diff_XY_distribution = (TH1D*)ESD_Conversion_XY_input1->Clone();
			Diff_XY_distribution->Divide(ESD_Conversion_XY_input1,ESD_Conversion_XY_input2,1.,1.,"B");
			Diff_XY_distribution->Draw("col2");	
			DrawAliceLogo(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3]);
			c_SinglePlot_5->Update();
			c_SinglePlot_5->SaveAs(Form("%sDiff_XY_distributions.%s",path,suffix));
			delete c_SinglePlot_5;
		}
		

		//-------------------- Z - Distribution ------------------------------

		TCanvas * c_SinglePlot_19 = new TCanvas("c_SinglePlot_19","",1000,1000);  // gives the page size
		c_SinglePlot_19->SetLogy(0);
		c_SinglePlot_19->cd();
		if(secondinput == ""){
				DrawAutoGammaHisto( ESD_Conversion_Z_input1
								 "Conversions Z distribution","Z [cm]",StandardYAxis,
								 kTRUE, 1.1,0.00001,
								 kFALSE,0. ,0,
								 kTRUE, -201.,201.);
		}else{
				DrawAutoGammaHistos( ESD_Conversion_Z_input1, 
								 ESD_Conversion_Z_input2, 
								 "Conversions Z distribution","Z [cm] ",StandardYAxis,
								 kTRUE, 1.1,0.00001,
								 kFALSE,0. ,0.,
								 kFALSE, -201.,201.);
		}
		for(Int_t i=0; i < 12 ; i++){
			DrawGammaLines(ArrayZbins[i], ArrayZbins[i], 0.00001,ESD_Conversion_Z_input1->GetMaximum());
		}		
		
		leg1 = new TLegend( 0.6,0.82,0.9,0.9);
		leg1->SetTextSize(0.04);			
		leg1->SetFillColor(0);
		leg1->AddEntry(ESD_Conversion_Z_input1,("Data"));
		if(secondinput != ""){leg1->AddEntry(ESD_Conversion_Z_input2,("MC"));}
		leg1->Draw();

		DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);	

		c_SinglePlot_19->Update();
		c_SinglePlot_19->SaveAs(Form("%sZ_distribution_lin.%s",path,suffix));
		delete c_SinglePlot_19;

		TCanvas * c_SinglePlot_6 = new TCanvas("c_SinglePlot_6","",1000,1000);  // gives the page size
		c_SinglePlot_6->SetLogy(1);
		c_SinglePlot_6->cd();
		if(secondinput == ""){
				DrawAutoGammaHisto( ESD_Conversion_Z_input1
								 "Conversions Z distribution","Z [cm]",StandardYAxis,
								 kTRUE, 1.5,0.00001,
								 kFALSE,0. ,0,
								 kTRUE, -250.,250.);
		}else{
				DrawAutoGammaHistos( ESD_Conversion_Z_input1, 
								 ESD_Conversion_Z_input2, 
								 "Conversions Z distribution","Z [cm] ",StandardYAxis,
								 kTRUE, 1.5,0.0000001,
								 kFALSE,0. ,0.,
								 kFALSE, -201.,201.);
		}
		for(Int_t i=0; i < 12 ; i++){
			DrawGammaLines(ArrayZbins[i], ArrayZbins[i], 0.0000001,ESD_Conversion_Z_input1->GetMaximum());
		}		
		leg1->Draw();
		DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);	

		c_SinglePlot_6->Update();
		c_SinglePlot_6->SaveAs(Form("%sZ_distribution_log.%s",path,suffix));
		delete c_SinglePlot_6;

		
		//------------------- Giving the Pad Phi in R in better resolution

		c_Single_Phi_in_R = new TCanvas("c_Single_Phi_in_R","",400,20,1400,2000);  // gives the page size
		pad_Single_Phi_in_R = new TPad("pad_Single_Phi_in_R","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
		pad_Single_Phi_in_R->SetFillColor(0);
		pad_Single_Phi_in_R->GetFrame()->SetFillColor(0);
		pad_Single_Phi_in_R->SetBorderMode(0);
		pad_Single_Phi_in_R->Divide(column,raw);
		pad_Single_Phi_in_R->Draw();

		for(Int_t iR = 0; iR < NbinsR; iR++){
			Int_t place = iR + 1;
			
			pad_Single_Phi_in_R->cd(place);
			pad_Single_Phi_in_R->cd(place)->SetLogy(1);
			sprintf(histoname_Phi_in_R_input1,"ESD_Conversion_Mapping_Phi_in_R_iR%02d",iR);
			if(secondinput == ""){		
				DrawAutoGammaHisto( ESD_Conversion_Mapping_Phi_in_R_input1[iR], 
									 histoname_Phi_in_R_input1,"#Phi",StandardYAxis,
									 kTRUE,2 ,0.0000002,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}else{
				DrawAutoGammaHistos( ESD_Conversion_Mapping_Phi_in_R_input1[iR], 
									 ESD_Conversion_Mapping_Phi_in_R_input2[iR], 
									 histoname_Phi_in_R_input1,"#Phi",StandardYAxis,
									 kTRUE,2 ,0.0000002,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
		}
			pad_Single_Phi_in_R->Update();
			c_Single_Phi_in_R->Print(Form("%spad_Phi_in_R.%s",path,suffix));	
			delete	pad_Single_Phi_in_R;
			delete c_Single_Phi_in_R;		


		// --------------- Giving Phi in R for several bins -------------------------------------------------------------
		TCanvas * c_SinglePlot_7 = new TCanvas("c_SinglePlot_7","",1000,700);  // gives the page size
		c_SinglePlot_7->SetLogy(1);
		c_SinglePlot_7->cd();
			if(secondinput == ""){		
				DrawAutoGammaHisto( ESD_Conversion_Mapping_Phi_in_R_input1[4],
									"ESD_Conversion_Mapping_Phi_in_R_R04","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}else{
				DrawAutoGammaHistos( ESD_Conversion_Mapping_Phi_in_R_input1[4], 
									 ESD_Conversion_Mapping_Phi_in_R_input2[4], 
									 "ESD_Conversion_Mapping_Phi_in_R_R04","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		

		c_SinglePlot_7->Update();
		c_SinglePlot_7->SaveAs(Form("%sPhi_in_R_04.%s",path,suffix));
		delete c_SinglePlot_7;
		
		TCanvas * c_SinglePlot_8 = new TCanvas("c_SinglePlot_8","",1000,700);  // gives the page size		
		c_SinglePlot_8->SetLogy(1);		
		c_SinglePlot_8->cd();
			if(secondinput == ""){		
				DrawAutoGammaHisto( ESD_Conversion_Mapping_Phi_in_R_input1[5],
									"ESD_Conversion_Mapping_Phi_in_R_R05","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}else{
				DrawAutoGammaHistos( ESD_Conversion_Mapping_Phi_in_R_input1[5], 
									 ESD_Conversion_Mapping_Phi_in_R_input2[5], 
									 "ESD_Conversion_Mapping_Phi_in_R_R05","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
		c_SinglePlot_8->Update();
		c_SinglePlot_8->SaveAs(Form("%sPhi_in_R_05.%s",path,suffix));
		delete c_SinglePlot_8;

		TCanvas * c_SinglePlot_12 = new TCanvas("c_SinglePlot_12","",1000,700);  // gives the page size		
		c_SinglePlot_12->SetLogy(1);
		c_SinglePlot_12->cd();
			if(secondinput == ""){		
				DrawAutoGammaHisto( ESD_Conversion_Mapping_Phi_in_R_input1[2],
									"ESD_Conversion_Mapping_Phi_in_R_R02","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}else{
				DrawAutoGammaHistos( ESD_Conversion_Mapping_Phi_in_R_input1[2], 
									 ESD_Conversion_Mapping_Phi_in_R_input2[2], 
									 "ESD_Conversion_Mapping_Phi_in_R_R02","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		

		c_SinglePlot_12->Update();
		c_SinglePlot_12->SaveAs(Form("%sPhi_in_R_02.%s",path,suffix));
		delete c_SinglePlot_12;				

		TCanvas * c_SinglePlot_13 = new TCanvas("c_SinglePlot_13","",1000,700);  // gives the page size
		c_SinglePlot_13->SetLogy(1);
		c_SinglePlot_13->cd();
			if(secondinput == ""){		
				DrawAutoGammaHisto( ESD_Conversion_Mapping_Phi_in_R_input1[6],
									"ESD_Conversion_Mapping_Phi_in_R_R06","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}else{
				DrawAutoGammaHistos( ESD_Conversion_Mapping_Phi_in_R_input1[6], 
									 ESD_Conversion_Mapping_Phi_in_R_input2[6], 
									 "ESD_Conversion_Mapping_Phi_in_R_R06","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
		c_SinglePlot_13->Update();
		c_SinglePlot_13->SaveAs(Form("%sPhi_in_R_06.%s",path,suffix));
		delete c_SinglePlot_13;		

		TCanvas * c_SinglePlot_14 = new TCanvas("c_SinglePlot_14","",1000,700);  // gives the page size
		c_SinglePlot_14->SetLogy(1);
		c_SinglePlot_14->cd();		
		if(secondinput == ""){		
				DrawAutoGammaHisto( ESD_Conversion_Mapping_Phi_in_R_input1[7],
									"ESD_Conversion_Mapping_Phi_in_R_R07","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}else{
				DrawAutoGammaHistos( ESD_Conversion_Mapping_Phi_in_R_input1[7], 
									 ESD_Conversion_Mapping_Phi_in_R_input2[7], 
									 "ESD_Conversion_Mapping_Phi_in_R_R07","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		

		c_SinglePlot_14->Update();
		c_SinglePlot_14->SaveAs(Form("%sPhi_in_R_07.%s",path,suffix));
		delete c_SinglePlot_14;				

		TCanvas * c_SinglePlot_15 = new TCanvas("c_SinglePlot_15","",1000,700);  // gives the page size
		c_SinglePlot_15->SetLogy(1);
		c_SinglePlot_15->cd();
			if(secondinput == ""){		
				DrawAutoGammaHisto( ESD_Conversion_Mapping_Phi_in_R_input1[9 ],
									"ESD_Conversion_Mapping_Phi_in_R_R09","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}else{
				DrawAutoGammaHistos( ESD_Conversion_Mapping_Phi_in_R_input1[9], 
									 ESD_Conversion_Mapping_Phi_in_R_input2[9], 
									 "ESD_Conversion_Mapping_Phi_in_R_R09","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
		c_SinglePlot_15->Update();
		c_SinglePlot_15->SaveAs(Form("%sPhi_in_R_09.%s",path,suffix));
		delete c_SinglePlot_15;

		TCanvas * c_SinglePlot_16 = new TCanvas("c_SinglePlot_16","",1000,700);  // gives the page size		
		c_SinglePlot_16->SetLogy(1);
		c_SinglePlot_16->cd();
			if(secondinput == ""){		
				DrawAutoGammaHisto( ESD_Conversion_Mapping_Phi_in_R_input1[10],
									"ESD_Conversion_Mapping_Phi_in_R_R10","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}else{
				DrawAutoGammaHistos( ESD_Conversion_Mapping_Phi_in_R_input1[10], 
									 ESD_Conversion_Mapping_Phi_in_R_input2[10], 
									 "ESD_Conversion_Mapping_Phi_in_R_R10","#Phi",StandardYAxis,
									 kFALSE,3 ,0.000001,
									 kFALSE,0. ,0.,
									 kFALSE, 0.,180.);
			}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
		c_SinglePlot_16->Update();
		c_SinglePlot_16->SaveAs(Form("%sPhi_in_R_10.%s",path,suffix));
		delete c_SinglePlot_16;
		
	//----------- Giving Z in R as single Pad ----------------
		c_Single_Z_in_R = new TCanvas("c_Single_Z_in_R","",400,20,1400,2000);  // gives the page size
	
		pad_Single_Z_in_R = new TPad("pad_Single_Z_in_R","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
		pad_Single_Z_in_R->SetFillColor(0);
		pad_Single_Z_in_R->GetFrame()->SetFillColor(0);
		pad_Single_Z_in_R->SetBorderMode(0);
		pad_Single_Z_in_R->Divide(column,raw);
		pad_Single_Z_in_R->Draw();
	
/*		title0->Draw(); 
		if(secondinput != ""){title1->Draw();}	
*/		
	for(Int_t iR = 0; iR < NbinsR; iR++){
		Int_t place = iR + 1;
		pad_Single_Z_in_R->cd(place);
		pad_Single_Z_in_R->cd(place)->SetLogy(1);
		sprintf(histoname_Z_in_R_input1,"ESD_Conversion_Mapping_Z_in_R_iR%02d",iR);
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_Z_in_R_input1[iR], 
								 histoname_Z_in_R_input1,"Z [cm]",StandardYAxis,
								 kTRUE,2 ,0.00001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos( ESD_Conversion_Mapping_Z_in_R_input1[iR], 
								 ESD_Conversion_Mapping_Z_in_R_input2[iR], 
								 histoname_Z_in_R_input1,"Z [cm]",StandardYAxis,
								 kTRUE,2 ,0.00001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
	}
	
		pad_Single_Z_in_R->Update();
		c_Single_Z_in_R->SaveAs(Form("%spad_Z_in_R.%s",path,suffix));
		delete pad_Single_Z_in_R;
		delete c_Single_Z_in_R;

		// -------------- Giving Z in R for 4th and 5th bin
		TCanvas * c_SinglePlot_9 = new TCanvas("c_SinglePlot_9","",1000,1000);  // gives the page size
		c_SinglePlot_9->SetLogy(1);		
		c_SinglePlot_9->cd();
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_Z_in_R_input1[9], 
								 "ESD_Conversion_Mapping_Z_in_R_iR09","Z [cm]",StandardYAxis,
								 kFALSE,3 ,0.0001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos( ESD_Conversion_Mapping_Z_in_R_input1[9], 
								 ESD_Conversion_Mapping_Z_in_R_input2[9], 
								 "ESD_Conversion_Mapping_Z_in_R_iR09","Z [cm]",StandardYAxis,
								 kFALSE,3 ,0.0001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
		c_SinglePlot_9->Update();
		c_SinglePlot_9->SaveAs(Form("%sZ_in_R_09.%s",path,suffix));
		delete c_SinglePlot_9;
				
		TCanvas * c_SinglePlot_10 = new TCanvas("c_SinglePlot_10","",1000,1000);  // gives the page size
		c_SinglePlot_10->SetLogy(1);		
		c_SinglePlot_10->cd();
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_Z_in_R_input1[10], 
								 "ESD_Conversion_Mapping_Z_in_R_iR10","Z [cm]",StandardYAxis,
								 kFALSE,3 ,0.0001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos( ESD_Conversion_Mapping_Z_in_R_input1[10], 
								 ESD_Conversion_Mapping_Z_in_R_input2[10], 
								 "ESD_Conversion_Mapping_Z_in_R_iR10","Z [cm]",StandardYAxis,
								 kFALSE,3 ,0.0001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
		c_SinglePlot_10->Update();
		c_SinglePlot_10->SaveAs(Form("%sZ_in_R_10.%s",path,suffix));
		delete c_SinglePlot_10;

		// ------------- Giving Z in phi in singleplot
		c_Single_Phi_in_Z = new TCanvas("c_Single_Phi_in_Z","",200,10,1400,2000);  // gives the page size
	
		pad_Single_Phi_in_Z = new TPad("pad_Single_Phi_in_Z","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
		pad_Single_Phi_in_Z->SetFillColor(0);
		pad_Single_Phi_in_Z->GetFrame()->SetFillColor(0);
		pad_Single_Phi_in_Z->SetBorderMode(0);
		pad_Single_Phi_in_Z->Divide(column,raw);
		pad_Single_Phi_in_Z->Draw();
	

	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		Int_t place = iZ + 1;
		pad_Single_Phi_in_Z->cd(place);
		pad_Single_Phi_in_Z->cd(place)->SetLogy(1);
		sprintf(histoname_Phi_in_Z_input1,"ESD_Conversion_Mapping_Phi_in_Z_iZ%02d",iZ);
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_Phi_in_Z_input1[iZ], 
								 histoname_Phi_in_Z_input1,"#Phi",StandardYAxis,
								 kTRUE,2 ,0.0000002,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos( ESD_Conversion_Mapping_Phi_in_Z_input1[iZ], 
								 ESD_Conversion_Mapping_Phi_in_Z_input2[iZ], 
								 histoname_Phi_in_Z_input1,"#Phi",StandardYAxis,
								 kTRUE,2 ,0.0000002,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
	}

		pad_Single_Phi_in_Z->Update();
		c_Single_Phi_in_Z->SaveAs(Form("%spad_Phi_in_Z.%s",path,suffix));
		delete pad_Single_Phi_in_Z ;
		delete c_Single_Phi_in_Z;

		// ----- Giving R in Z in SinglePlot	
		c_Single_R_in_Z = new TCanvas("c_Single_R_in_Z","",400,20,1400,2000);  // gives the page size
	
		pad_Single_R_in_Z = new TPad("pad_Single_R_in_Z","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
		pad_Single_R_in_Z->SetFillColor(0);
		pad_Single_R_in_Z->GetFrame()->SetFillColor(0);
		pad_Single_R_in_Z->SetBorderMode(0);
		pad_Single_R_in_Z->Divide(column,raw);
		pad_Single_R_in_Z->Draw();
	

	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		Int_t place = iZ + 1;
		pad_Single_R_in_Z->cd(place);
		pad_Single_R_in_Z->cd(place)->SetLogy(1);
		sprintf(histoname_R_in_Z_input1,"ESD_Conversion_Mapping_R_in_Z_iZ%02d",iZ);
		if(secondinput == ""){		
			DrawAutoGammaHisto( ESD_Conversion_Mapping_R_in_Z_input1[iZ], 
								 histoname_R_in_Z_input1,"R [cm]",StandardYAxis,
								 kTRUE,2 ,0.0000002,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}else{
			DrawAutoGammaHistos( ESD_Conversion_Mapping_R_in_Z_input1[iZ], 
								 ESD_Conversion_Mapping_R_in_Z_input2[iZ], 
								 histoname_R_in_Z_input1,"R [cm]",StandardYAxis,
								 kTRUE,2 ,0.0000002,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,180.);
		}	
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
	}

		pad_Single_R_in_Z->Update();
		c_Single_R_in_Z->SaveAs(Form("%spad_R_in_Z.%s",path,suffix));
		delete pad_Single_R_in_Z ;
		delete c_Single_R_in_Z;
		//---------------------------- Ratio Pads ----------------------------------

if(secondinput != ""){

	c_Ratio_Phi_in_R = new TCanvas("c_Ratio_Phi_in_R","",400,20,1400,2000);  // gives the page size
	
	pad_Ratio_Phi_in_R = new TPad("pad_Ratio_Phi_in_R","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	pad_Ratio_Phi_in_R->SetFillColor(0);
	pad_Ratio_Phi_in_R->GetFrame()->SetFillColor(0);
	pad_Ratio_Phi_in_R->SetBorderMode(0);
	pad_Ratio_Phi_in_R->Divide(column,raw);
	pad_Ratio_Phi_in_R->Draw();
	
	for(Int_t iR = 0; iR < NbinsR; iR++){
		Int_t place = iR + 1;
		
		pad_Ratio_Phi_in_R->cd(place);
		pad_Ratio_Phi_in_R->cd(place)->SetLogy(1);
		sprintf(histoname_Ratio_Phi_in_R_input2,"Ratio_Phi_in_R_iR%02d",iR);
			DrawRatioGammaHisto( 	Ratio_Mapping_Phi_in_R[iR], 
								 histoname_Ratio_Phi_in_R_input2,"#Phi","norm Data/norm MC",
								 kFALSE,3 ,0.000001,
								 kTRUE,0.1 ,5.,
								 kFALSE, 0.,180.);
			linePhi->Draw("same");
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
	}
	
	pad_Ratio_Phi_in_R->Update();
	c_Ratio_Phi_in_R->SaveAs(Form("%sRatio_Phi_in_R.%s",path,suffix));
	delete pad_Ratio_Phi_in_R ;
	delete	c_Ratio_Phi_in_R;
}

if(secondinput != ""){

	c_Ratio_Z_in_R = new TCanvas("c_Ratio_Z_in_R","",400,20,1400,2000);  // gives the page size
	
	pad_Ratio_Z_in_R = new TPad("pad_Ratio_Z_in_R","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	pad_Ratio_Z_in_R->SetFillColor(0);
	pad_Ratio_Z_in_R->GetFrame()->SetFillColor(0);
	pad_Ratio_Z_in_R->SetBorderMode(0);
	pad_Ratio_Z_in_R->Divide(column,raw);
	pad_Ratio_Z_in_R->Draw();

	for(Int_t iR = 0; iR < NbinsR; iR++){
		Int_t place = iR + 1;
		
		pad_Ratio_Z_in_R->cd(place);
		pad_Ratio_Z_in_R->cd(place)->SetLogy(1);
		sprintf(histoname_Ratio_Z_in_R_input2,"Ratio_Z_in_R_iR%02d",iR);	
			Float_t ZRange = RangeZinR[iR];
			DrawRatioGammaHisto( Ratio_Mapping_Z_in_R[iR], 
								 histoname_Ratio_Z_in_R_input2,"Z","norm Data/norm MC",
								 kFALSE,3 ,0.000001,
								 kTRUE,0.1 ,5.,
								 kTRUE, -ZRange,ZRange);
			DrawGammaLines(-ZRange,ZRange,1,1);
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
	}
	pad_Ratio_Z_in_R->Update();
	c_Ratio_Z_in_R->SaveAs(Form("%sRatio_Z_in_R.%s",path,suffix));
	delete pad_Ratio_Z_in_R ;
	delete c_Ratio_Z_in_R;
}

if(secondinput != ""){

	c_Ratio_Phi_in_Z = new TCanvas("c_Ratio_Phi_in_Z","",400,20,1400,2000);  // gives the page size
	
	pad_Ratio_Phi_in_Z = new TPad("pad_Ratio_Phi_in_Z","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	pad_Ratio_Phi_in_Z->SetFillColor(0);
	pad_Ratio_Phi_in_Z->GetFrame()->SetFillColor(0);
	pad_Ratio_Phi_in_Z->SetBorderMode(0);
	pad_Ratio_Phi_in_Z->Divide(column,raw);
	pad_Ratio_Phi_in_Z->Draw();

	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		Int_t place = iZ + 1;
		
		pad_Ratio_Phi_in_Z->cd(place);
		pad_Ratio_Phi_in_Z->cd(place)->SetLogy(1);
		sprintf(histoname_Ratio_Phi_in_Z_input2,"Ratio_Phi_in_Z_iZ%02d",iZ);
			DrawRatioGammaHisto( Ratio_Mapping_Phi_in_Z[iZ], 
								 histoname_Ratio_Phi_in_Z_input2,"#Phi","norm Data/norm MC",
								 kFALSE,3 ,0.000001,
								 kTRUE,0.1 ,5.,
								 kFALSE, 0,0);
			linePhi->Draw("same");
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
	}
	
	pad_Ratio_Phi_in_Z->Update();
	c_Ratio_Phi_in_Z->SaveAs(Form("%sRatio_Phi_in_Z.%s",path,suffix));
	delete pad_Ratio_Phi_in_Z ;
	delete c_Ratio_Phi_in_Z;
}

if(secondinput != ""){

	c_Ratio_R_in_Z = new TCanvas("c_Ratio_R_in_Z","",400,20,1400,2000);  // gives the page size
	
	pad_Ratio_R_in_Z = new TPad("pad_Ratio_R_in_Z","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	pad_Ratio_R_in_Z->SetFillColor(0);
	pad_Ratio_R_in_Z->GetFrame()->SetFillColor(0);
	pad_Ratio_R_in_Z->SetBorderMode(0);
	pad_Ratio_R_in_Z->Divide(column,raw);
	pad_Ratio_R_in_Z->Draw();
	
	for(Int_t iZ = 0; iZ < NbinsZ; iZ++){
		Int_t place = iZ + 1;
		
		pad_Ratio_R_in_Z->cd(place);
		pad_Ratio_R_in_Z->cd(place)->SetLogy(1);
		sprintf(histoname_Ratio_R_in_Z_input2,"Ratio_R_in_Z_iZ%02d",iZ);
			DrawRatioGammaHisto( Ratio_Mapping_R_in_Z[iZ], 
								 histoname_Ratio_R_in_Z_input2,"R [cm]","norm Data/norm MC",
								 kFALSE,3 ,0.000001,
								 kTRUE,0.1 ,5.,
								 kFALSE, 0,0);
			lineR->Draw("same");
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);		
	}
	
	pad_Ratio_R_in_Z->Update();
	c_Ratio_R_in_Z->SaveAs(Form("%sRatio_R_in_Z.%s",path,suffix));
	delete pad_Ratio_R_in_Z ;
	delete c_Ratio_R_in_Z;
}	
	
	}	

//Dealocating all reserved resources

}
