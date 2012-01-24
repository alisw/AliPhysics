/***********************************************************************************************
*** provided by Gamma Conversion Group, PWGGA, 									******
***	    Friederike Bock, fbock@physi.uni-heidelberg.de ***							******
************************************************************************************************
************************************************************************************************
*** With this Macro Resolution Studies for the Conversion Method can be carried out 		******
*** For this a run with the flag "resolution" of the Gamma Conversion Software needs 	******
*** to be performed, this study can only be carried out on Montecarlo Productions 		******
************************************************************************************************/

#include <Riostream.h>
#include <TMath.h>
#include "PlottingGammaConversionHistos.h"
#include "PlottingGammaConversionAdditional.h"
extern TRandom *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem *gSystem;

void  Photon_Resolution(const char *MCfile = "", const char *cutsel = "",const char *path = "", const char *output = "Photon-Resolution", const char *Plots = "kTRUE", const char *suffix = "eps"){	
	
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettings();	
	set_plot_style();

	Bool_t SinglePlots = kFALSE;
	if(Plots == "kTRUE") SinglePlots = kTRUE;


//**********************************************************************************************************************
//******************************** Definition of some plotting variables ***********************************************
//**********************************************************************************************************************

	//Array defintion for printing Logo 
	Float_t right_up[4]={0.7,0.63,0.15, 0.02};
	Float_t right_up2D[4]={0.65,0.63,0.11, 0.02};
	Float_t right_down[4]={0.7,0.23,0.15, 0.02};
	Float_t left_up[4]={0.17,0.73, 0.15, 0.02};

	//Defintion of Rebin
	Float_t rebin = 4;	

//**********************************************************************************************************************
//******************************** Defintion of arrays for rebinning ***************************************************
//**********************************************************************************************************************

	//pt rebinning array
	Float_t ptcalcbinning1[5] =  {0.5 ,1 ,2 ,5,10};
	Float_t ptcalcbinning2[4] = {14,40,100,200};
	Float_t ptcalcbinning3[6] = {1,5,13,21,24,31};
	Double_t ptbinning[31];
	ptbinning[0] = 0;
	for ( int j = 0; j < 5 ; j ++ ){
		for ( int i = ptcalcbinning3[j]; i<ptcalcbinning3[j+1]; i++){ 
		ptbinning[i] = ptbinning[i-1] +	 ptcalcbinning1[j];
		}
	}
	
	//R rebinning array
	Float_t Rcalcbinning1[18] =  {5, 0.5, 2, 0.5, 4,0.5,1,0.5,3.5,0.5,4,0.5,5,0.5,4,1,11,15};
	Float_t Rcalcbinning3[19] = {1,2,16,17,27,28,38,39,41,43,65,68,71,74,79,80,81,82,89};
	Double_t Rbinning[89];
	Rbinning[0] = 0;
	for ( int j = 0; j < 18 ; j ++ ){
		for ( int i = Rcalcbinning3[j]; i<Rcalcbinning3[j+1]; i++){ 
		Rbinning[i] = Rbinning[i-1] +	 Rcalcbinning1[j];
		}
	}
//**********************************************************************************************************************
//******************************************* Defintion of file name ***************************************************
//**********************************************************************************************************************

	Char_t filename_MC[100] = (Form("%s%s",path,MCfile));

	if(cutsel != ""){
		char *GammaDirectory = Form("PWGGA_GammaConversion_%s",  cutsel);
		cout << GammaDirectory << endl;
		char *GammaList = Form("histogramsAliGammaConversion_%s", cutsel);
		cout << GammaList << endl;
	}else{
		char *GammaDirectory = "PWGGA_GammaConversion";
		cout << GammaDirectory << endl;
		char *GammaList = "histogramsAliGammaConversion";
		cout << GammaList << endl;
	}

//**********************************************************************************************************************
//****************************************** Loading of Histograms *****************************************************
//**********************************************************************************************************************

	TFile f(filename_MC);  
	
	//************************** Container Loading ********************************************************************
	TDirectory *fPWGGAGammaConversion_montecarlo = new TDirectory(); 
	TList *fHistosGammaConversion_montecarlo = new TList(); 
	TList *fResolutionContainer_montecarlo = new TList();  
     TList *fESDContainer_montecarlo = new TList();
	if(!(fPWGGAGammaConversion_montecarlo = (TDirectory*)f.Get(GammaDirectory))) cout <<"PWGGAGammConversion TList NOT loaded correctly"<<endl; 
	if(!(fHistosGammaConversion_montecarlo = (TList*)fPWGGAGammaConversion_montecarlo->Get(GammaList))) cout<<"histogramsAliGammaConversion NOT loaded correctly!"<<endl; 
	if(fResolutionContainer_montecarlo = (TList*)fHistosGammaConversion_montecarlo->FindObject("Resolution histograms")) cout<<"list loaded correctly!"<<endl; 
	if(!(fESDContainer_montecarlo = (TList*)fHistosGammaConversion_montecarlo->FindObject("ESD histograms"))) cout<<"ESD histograms NOT loaded correctly!"<<endl; 
	TH1F * ESD_ConvGamma_Pt=fESDContainer_montecarlo->FindObject("ESD_ConvGamma_Pt");
	TH1F * ESD_E_Pt=fESDContainer_montecarlo->FindObject("ESD_E_Pt");
	TH1F * ESD_P_Pt=fESDContainer_montecarlo->FindObject("ESD_P_Pt");

	//************************** Histo loading *************************************************************************

	//****************************** Resolution dPt vs Pt **********************************************
	TH1F * Resolution_MC_Pt = fResolutionContainer_montecarlo->FindObject("Resolution_MC_Pt" ); //
	TH1F * Resolution_ESD_Pt = fResolutionContainer_montecarlo->FindObject("Resolution_ESD_Pt" ); //
	TH2F * Resolution_Gamma_dPt_Pt = fResolutionContainer_montecarlo->FindObject("Resolution_Gamma_dPt_Pt" ); //
	TH2F * Resolution_E_dPt_Pt_ITS0 = fResolutionContainer_montecarlo->FindObject("Resolution_E_dPt_Pt_ITS0" ); //
	TH2F * Resolution_E_dPt_Pt_ITS1 = fResolutionContainer_montecarlo->FindObject("Resolution_E_dPt_Pt_ITS1" ); //
	TH2F * Resolution_E_dPt_Pt_ITS2 = fResolutionContainer_montecarlo->FindObject("Resolution_E_dPt_Pt_ITS2" ); //
	TH2F * Resolution_E_dPt_Pt_ITS3 = fResolutionContainer_montecarlo->FindObject("Resolution_E_dPt_Pt_ITS3" ); //
	TH2F * Resolution_E_dPt_Pt_ITS4 = fResolutionContainer_montecarlo->FindObject("Resolution_E_dPt_Pt_ITS4" ); //
	TH2F * Resolution_E_dPt_Pt_ITS5 = fResolutionContainer_montecarlo->FindObject("Resolution_E_dPt_Pt_ITS5" ); //
	TH2F * Resolution_E_dPt_Pt_ITS6 = fResolutionContainer_montecarlo->FindObject("Resolution_E_dPt_Pt_ITS6" ); //
	TH2F * Resolution_E_dPt_Pt = fResolutionContainer_montecarlo->FindObject("Resolution_E_dPt_Pt" ); //
	TH2F * Resolution_P_dPt_Pt_ITS0 = fResolutionContainer_montecarlo->FindObject("Resolution_P_dPt_Pt_ITS0" ); //
	TH2F * Resolution_P_dPt_Pt_ITS1 = fResolutionContainer_montecarlo->FindObject("Resolution_P_dPt_Pt_ITS1" ); //
	TH2F * Resolution_P_dPt_Pt_ITS2 = fResolutionContainer_montecarlo->FindObject("Resolution_P_dPt_Pt_ITS2" ); //
	TH2F * Resolution_P_dPt_Pt_ITS3 = fResolutionContainer_montecarlo->FindObject("Resolution_P_dPt_Pt_ITS3" ); //
	TH2F * Resolution_P_dPt_Pt_ITS4 = fResolutionContainer_montecarlo->FindObject("Resolution_P_dPt_Pt_ITS4" ); //
	TH2F * Resolution_P_dPt_Pt_ITS5 = fResolutionContainer_montecarlo->FindObject("Resolution_P_dPt_Pt_ITS5" ); //
	TH2F * Resolution_P_dPt_Pt_ITS6 = fResolutionContainer_montecarlo->FindObject("Resolution_P_dPt_Pt_ITS6" ); //
	TH2F * Resolution_P_dPt_Pt = fResolutionContainer_montecarlo->FindObject("Resolution_P_dPt_Pt" ); //

	//****************************** TRD performance ***************************************************
	TH2F * Resolution_E_nTRDtracklets_ESDPt = fResolutionContainer_montecarlo->FindObject("Resolution_E_nTRDtracklets_ESDPt" );
//	TH2F * Resolution_E_nTRDtracklets_MCPt = fResolutionContainer_montecarlo->FindObject("Resolution_E_nTRDtracklets_MCPt" );
	TH2F * Resolution_E_nTRDclusters_ESDPt = fResolutionContainer_montecarlo->FindObject("Resolution_E_nTRDclusters_ESDPt" );
//	TH2F * Resolution_E_nTRDclusters_MCPt = fResolutionContainer_montecarlo->FindObject("Resolution_E_nTRDclusters_MCPt" );
	TH2F * Resolution_P_nTRDtracklets_ESDPt = fResolutionContainer_montecarlo->FindObject("Resolution_P_nTRDtracklets_ESDPt" );
//	TH2F * Resolution_P_nTRDtracklets_MCPt = fResolutionContainer_montecarlo->FindObject("Resolution_P_nTRDtracklets_MCPt" );
	TH2F * Resolution_P_nTRDclusters_ESDPt = fResolutionContainer_montecarlo->FindObject("Resolution_P_nTRDclusters_ESDPt" );
//	TH2F * Resolution_P_nTRDclusters_MCPt = fResolutionContainer_montecarlo->FindObject("Resolution_P_nTRDclusters_MCPt" );

	//******************************** Resolution in Absolute ********************************************
	TH2F * Resolution_dZAbs_VS_R = fResolutionContainer_montecarlo->FindObject("Resolution_dZAbs_VS_R" ); //
	TH2F * Resolution_dPhiAbs_VS_R = fResolutionContainer_montecarlo->FindObject("Resolution_dPhiAbs_VS_R" ); //
	TH2F * Resolution_dRAbs_VS_R = fResolutionContainer_montecarlo->FindObject("Resolution_dRAbs_VS_R" ); //

	// ******************************* Resolution in R ***************************************************
	TH1F * Resolution_MC_R = fResolutionContainer_montecarlo->FindObject("Resolution_MC_R" ); //
	TH1F * Resolution_ESD_R = fResolutionContainer_montecarlo->FindObject("Resolution_ESD_R" ); //
//	TH2F * Resolution_dR = fResolutionContainer_montecarlo->FindObject("Resolution_dR" );	
	TH2F * Resolution_R_dPt = fResolutionContainer_montecarlo->FindObject("Resolution_R_dPt" );

	//******************************** Resolution in Z ***************************************************
	TH1F * Resolution_MC_Z = fResolutionContainer_montecarlo->FindObject("Resolution_MC_Z" ); //
	TH1F * Resolution_ESD_Z = fResolutionContainer_montecarlo->FindObject("Resolution_ESD_Z" ); //
//	TH2F * Resolution_dZ = fResolutionContainer_montecarlo->FindObject("Resolution_dZ" ); 


	//******************************** Resolution dPt vs Phi *********************************************
	TH2F * Resolution_Gamma_dPt_Phi = fResolutionContainer_montecarlo->FindObject("Resolution_Gamma_dPt_Phi" );
	TH2F * Resolution_E_dPt_Phi = fResolutionContainer_montecarlo->FindObject("Resolution_E_dPt_Phi" );
	TH2F * Resolution_P_dPt_Phi = fResolutionContainer_montecarlo->FindObject("Resolution_P_dPt_Phi" ); 

//**********************************************************************************************************************************************
//*********************************** Rebinning of 2D histograms and fitslices *****************************************************************
//**********************************************************************************************************************************************
	
		//********************************* Rebinning Gamma dPt Pt *****************************************************************
		TH2F *Resolution_Gamma_dPt_Ptrebin = new TH2F("Resolution_Gamma_dPt_Ptrebin", "Photon Resolution dPt vs Pt ", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_Gamma_dPt_Pt->Sumw2();
		Int_t maxYbin = Resolution_Gamma_dPt_Pt->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_Gamma_dPt_Pt->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_Gamma_dPt_Pt->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_Gamma_dPt_Pt->GetBinContent(i, j);
					Resolution_Gamma_dPt_Ptrebin->Fill(binoldpoint, biny, value);
				}
		}

		//********************************* Rebinning of E dPt Pt all **************************************************************
		TH2F *Resolution_E_dPt_Ptrebin = new TH2F("Resolution_E_dPt_Ptrebin", "Electron Resolution dPt vs Pt ", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_E_dPt_Pt->Sumw2();
		Int_t maxYbin = Resolution_E_dPt_Pt->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_E_dPt_Pt->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_E_dPt_Pt->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_E_dPt_Pt->GetBinContent(i, j);
					Resolution_E_dPt_Ptrebin->Fill(binoldpoint, biny, value);
				}
		}
		
		//**********************************Rebinning of E dPt Pt for different Nr of ITS clusters *********************************
		TH2F *Resolution_E_dPt_Pt_ITS0rebin = new TH2F("Resolution_E_dPt_Pt_ITS0rebin", "Electron Resolution dPt vs Pt 0 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_E_dPt_Pt_ITS0->Sumw2();
		Int_t maxYbin = Resolution_E_dPt_Pt_ITS0->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_E_dPt_Pt_ITS0->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_E_dPt_Pt_ITS0->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_E_dPt_Pt_ITS0->GetBinContent(i, j);
					Resolution_E_dPt_Pt_ITS0rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_E_dPt_Pt_ITS1rebin = new TH2F("Resolution_E_dPt_Pt_ITS1rebin", "Electron Resolution dPt vs Pt 1 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_E_dPt_Pt_ITS1->Sumw2();
		Int_t maxYbin = Resolution_E_dPt_Pt_ITS1->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_E_dPt_Pt_ITS1->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_E_dPt_Pt_ITS1->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_E_dPt_Pt_ITS1->GetBinContent(i, j);
					Resolution_E_dPt_Pt_ITS1rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_E_dPt_Pt_ITS2rebin = new TH2F("Resolution_E_dPt_Pt_ITS2rebin", "Electron Resolution dPt vs Pt 2 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_E_dPt_Pt_ITS2->Sumw2();
		Int_t maxYbin = Resolution_E_dPt_Pt_ITS2->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_E_dPt_Pt_ITS2->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_E_dPt_Pt_ITS2->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_E_dPt_Pt_ITS2->GetBinContent(i, j);
					Resolution_E_dPt_Pt_ITS2rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_E_dPt_Pt_ITS3rebin = new TH2F("Resolution_E_dPt_Pt_ITS3rebin", "Electron Resolution dPt vs Pt 3 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_E_dPt_Pt_ITS3->Sumw2();
		Int_t maxYbin = Resolution_E_dPt_Pt_ITS3->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_E_dPt_Pt_ITS3->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_E_dPt_Pt_ITS3->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_E_dPt_Pt_ITS3->GetBinContent(i, j);
					Resolution_E_dPt_Pt_ITS3rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_E_dPt_Pt_ITS4rebin = new TH2F("Resolution_E_dPt_Pt_ITS4rebin", "Electron Resolution dPt vs Pt 4 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_E_dPt_Pt_ITS4->Sumw2();
		Int_t maxYbin = Resolution_E_dPt_Pt_ITS4->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_E_dPt_Pt_ITS4->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_E_dPt_Pt_ITS4->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_E_dPt_Pt_ITS4->GetBinContent(i, j);
					Resolution_E_dPt_Pt_ITS4rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_E_dPt_Pt_ITS5rebin = new TH2F("Resolution_E_dPt_Pt_ITS5rebin", "Electron Resolution dPt vs Pt 5 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_E_dPt_Pt_ITS5->Sumw2();
		Int_t maxYbin = Resolution_E_dPt_Pt_ITS5->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_E_dPt_Pt_ITS5->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_E_dPt_Pt_ITS5->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_E_dPt_Pt_ITS5->GetBinContent(i, j);
					Resolution_E_dPt_Pt_ITS5rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_E_dPt_Pt_ITS6rebin = new TH2F("Resolution_E_dPt_Pt_ITS6rebin", "Electron Resolution dPt vs Pt 6 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_E_dPt_Pt_ITS6->Sumw2();
		Int_t maxYbin = Resolution_E_dPt_Pt_ITS6->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_E_dPt_Pt_ITS6->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_E_dPt_Pt_ITS6->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_E_dPt_Pt_ITS6->GetBinContent(i, j);
					Resolution_E_dPt_Pt_ITS6rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_P_dPt_Ptrebin = new TH2F("Resolution_P_dPt_Ptrebin", "Positron Resolution dPt vs Pt ", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_P_dPt_Pt->Sumw2();
		Int_t maxYbin = Resolution_P_dPt_Pt->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_P_dPt_Pt->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_P_dPt_Pt->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_P_dPt_Pt->GetBinContent(i, j);
					Resolution_P_dPt_Ptrebin->Fill(binoldpoint, biny, value);
				}
		}

		//**********************************Rebinning of P dPt Pt for different Nr of ITS clusters *********************************
		TH2F *Resolution_P_dPt_Pt_ITS0rebin = new TH2F("Resolution_P_dPt_Pt_ITS0rebin", "Positron Resolution dPt vs Pt 0 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_P_dPt_Pt_ITS0->Sumw2();
		Int_t maxYbin = Resolution_P_dPt_Pt_ITS0->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_P_dPt_Pt_ITS0->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_P_dPt_Pt_ITS0->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_P_dPt_Pt_ITS0->GetBinContent(i, j);
					Resolution_P_dPt_Pt_ITS0rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_P_dPt_Pt_ITS1rebin = new TH2F("Resolution_P_dPt_Pt_ITS1rebin", "Positron Resolution dPt vs Pt 1 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_P_dPt_Pt_ITS1->Sumw2();
		Int_t maxYbin = Resolution_P_dPt_Pt_ITS1->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_P_dPt_Pt_ITS1->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_P_dPt_Pt_ITS1->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_P_dPt_Pt_ITS1->GetBinContent(i, j);
					Resolution_P_dPt_Pt_ITS1rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_P_dPt_Pt_ITS2rebin = new TH2F("Resolution_P_dPt_Pt_ITS2rebin", "Positron Resolution dPt vs Pt 2 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_P_dPt_Pt_ITS2->Sumw2();
		Int_t maxYbin = Resolution_P_dPt_Pt_ITS2->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_P_dPt_Pt_ITS2->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_P_dPt_Pt_ITS2->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_P_dPt_Pt_ITS2->GetBinContent(i, j);
					Resolution_P_dPt_Pt_ITS2rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_P_dPt_Pt_ITS3rebin = new TH2F("Resolution_P_dPt_Pt_ITS3rebin", "Positron Resolution dPt vs Pt 3 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_P_dPt_Pt_ITS3->Sumw2();
		Int_t maxYbin = Resolution_P_dPt_Pt_ITS3->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_P_dPt_Pt_ITS3->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_P_dPt_Pt_ITS3->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_P_dPt_Pt_ITS3->GetBinContent(i, j);
					Resolution_P_dPt_Pt_ITS3rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_P_dPt_Pt_ITS4rebin = new TH2F("Resolution_P_dPt_Pt_ITS4rebin", "Positron Resolution dPt vs Pt 4 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_P_dPt_Pt_ITS4->Sumw2();
		Int_t maxYbin = Resolution_P_dPt_Pt_ITS4->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_P_dPt_Pt_ITS4->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_P_dPt_Pt_ITS4->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_P_dPt_Pt_ITS4->GetBinContent(i, j);
					Resolution_P_dPt_Pt_ITS4rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_P_dPt_Pt_ITS5rebin = new TH2F("Resolution_P_dPt_Pt_ITS5rebin", "Positron Resolution dPt vs Pt 5 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_P_dPt_Pt_ITS5->Sumw2();
		Int_t maxYbin = Resolution_P_dPt_Pt_ITS5->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_P_dPt_Pt_ITS5->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_P_dPt_Pt_ITS5->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_P_dPt_Pt_ITS5->GetBinContent(i, j);
					Resolution_P_dPt_Pt_ITS5rebin->Fill(binoldpoint, biny, value);
				}
		}

		TH2F *Resolution_P_dPt_Pt_ITS6rebin = new TH2F("Resolution_P_dPt_Pt_ITS6rebin", "Positron Resolution dPt vs Pt 6 ITS clusters", 30, ptbinning , 100, -10., 10.) ;	
		Resolution_P_dPt_Pt_ITS6->Sumw2();
		Int_t maxYbin = Resolution_P_dPt_Pt_ITS6->GetNbinsY();
		for (Int_t i = 1; i <201; i++) {
				TAxis *xaxis_old = 	Resolution_P_dPt_Pt_ITS6->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_P_dPt_Pt_ITS6->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_P_dPt_Pt_ITS6->GetBinContent(i, j);
					Resolution_P_dPt_Pt_ITS6rebin->Fill(binoldpoint, biny, value);
				}
		}

		//************************************** Rebinning of dZAbs R **********************************************
		TH2F *Resolution_dZAbs_VS_Rrebin = new TH2F("Resolution_dZAbs_VS_Rrebin", "Photon Resolution dZabs vs R ", 88, Rbinning , 50, -25., 25.) ;	
		Resolution_dZAbs_VS_R->Sumw2();
		Int_t maxXbin = Resolution_dZAbs_VS_R->GetNbinsX();		
		Int_t maxYbin = Resolution_dZAbs_VS_R->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_dZAbs_VS_R->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_dZAbs_VS_R->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_dZAbs_VS_R->GetBinContent(i, j);
					Resolution_dZAbs_VS_Rrebin->Fill(binoldpoint, biny, value);
				}
		}

		//************************************** Rebinning of dPhiAbs R **********************************************
		TH2F *Resolution_dPhiAbs_VS_Rrebin = new TH2F("Resolution_dPhiAbs_VS_Rrebin", "Photon Resolution dPhiabs vs R ", 88, Rbinning , 100, -TMath::Pi()/30, TMath::Pi()/30) ;	
		Resolution_dPhiAbs_VS_R->Sumw2();
		Int_t maxXbin = Resolution_dPhiAbs_VS_R->GetNbinsX();		
		Int_t maxYbin = Resolution_dPhiAbs_VS_R->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_dPhiAbs_VS_R->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_dPhiAbs_VS_R->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_dPhiAbs_VS_R->GetBinContent(i, j);
					Resolution_dPhiAbs_VS_Rrebin->Fill(binoldpoint, biny, value);
				}
		}

		//************************************** Rebinning of dRAbs R **********************************************
		TH2F *Resolution_dRAbs_VS_Rrebin = new TH2F("Resolution_dRAbs_VS_Rrebin", "Photon Resolution dRabs vs R ", 88, Rbinning , 50, -25,25) ;	
		Resolution_dRAbs_VS_R->Sumw2();
		Int_t maxXbin = Resolution_dRAbs_VS_R->GetNbinsX();		
		Int_t maxYbin = Resolution_dRAbs_VS_R->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_dRAbs_VS_R->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_dRAbs_VS_R->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_dRAbs_VS_R->GetBinContent(i, j);
					Resolution_dRAbs_VS_Rrebin->Fill(binoldpoint, biny, value);
				}
		}

		//************************************** Rebinning of dRAbs R **********************************************
		TH2F *Resolution_R_dPtrebin = new TH2F("Resolution_R_dPtrebin", "Photon Resolution dp_{t} vs R ", 88, Rbinning , 50, -10,10) ;	
		Resolution_R_dPt->Sumw2();
		Int_t maxXbin = Resolution_R_dPt->GetNbinsX();		
		Int_t maxYbin = Resolution_R_dPt->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_R_dPt->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_R_dPt->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_R_dPt->GetBinContent(i, j);
					Resolution_R_dPtrebin->Fill(binoldpoint, biny, value);
				}
		}




	//**************************** Writing all rebined histos in a file ********************************************
	Resolution = new TFile("ResolutionGamma.root","RECREATE");
		Resolution_Gamma_dPt_Ptrebin->Write();
		Resolution_E_dPt_Ptrebin->Write();
		Resolution_E_dPt_Pt_ITS0rebin->Write();
		Resolution_E_dPt_Pt_ITS1rebin->Write();
		Resolution_E_dPt_Pt_ITS2rebin->Write();
		Resolution_E_dPt_Pt_ITS3rebin->Write();
		Resolution_E_dPt_Pt_ITS4rebin->Write();
		Resolution_E_dPt_Pt_ITS5rebin->Write();
		Resolution_E_dPt_Pt_ITS6rebin->Write();
		Resolution_P_dPt_Ptrebin->Write();
		Resolution_P_dPt_Pt_ITS0rebin->Write();
		Resolution_P_dPt_Pt_ITS1rebin->Write();
		Resolution_P_dPt_Pt_ITS2rebin->Write();
		Resolution_P_dPt_Pt_ITS3rebin->Write();
		Resolution_P_dPt_Pt_ITS4rebin->Write();
		Resolution_P_dPt_Pt_ITS5rebin->Write();
		Resolution_P_dPt_Pt_ITS6rebin->Write();
		Resolution_dZAbs_VS_Rrebin->Write();
		Resolution_dRAbs_VS_Rrebin->Write();
		Resolution_dPhiAbs_VS_Rrebin->Write();
		Resolution_R_dPtrebin->Write();
	Resolution->Write();
	Resolution->Close();

		Double_t precision = 10E-5;

//*********************************** Fitting for dRabs vs R **********************************************
		TH1F *Resolution_dRAbs_VS_Rrebin_1 = new TH1F("Resolution_dRAbs_VS_Rrebin_1", "Mean Photon Resolution dRabs vs R ", 88, Rbinning ) ;	
		TH1F *Resolution_dRAbs_VS_Rrebin_2 = new TH1F("Resolution_dRAbs_VS_Rrebin_2", "Sigma Photon Resolution dRabs vs R ", 88, Rbinning ) ;
		for (Int_t i = 1; i < 88 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -25,25);
			TH1D *hy = Resolution_dRAbs_VS_Rrebin->ProjectionY("hy", i, i);
			hy->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2);
			Double_t ymax = mp + (rp * 2);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				hy->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp + (rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_dRAbs_VS_Rrebin_1->SetBinContent(i, mean);
			Resolution_dRAbs_VS_Rrebin_1->SetBinError(i, meanE);
			Resolution_dRAbs_VS_Rrebin_2->SetBinContent(i, sigma);
			Resolution_dRAbs_VS_Rrebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

//********************************************* fitting for dZAbs vs R *************************************************************
		TH1F *Resolution_dZAbs_VS_Rrebin_1 = new TH1F("Resolution_dZAbs_VS_Rrebin_1", "Mean Photon Resolution dZabs vs R ", 88, Rbinning ) ;	
		TH1F *Resolution_dZAbs_VS_Rrebin_2 = new TH1F("Resolution_dZAbs_VS_Rrebin_2", "Sigma Photon Resolution dZabs vs R ", 88, Rbinning ) ;
		for (Int_t i = 1; i < 88 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -25,25);
			TH1D *hx = Resolution_dZAbs_VS_Rrebin->ProjectionY("hx", i, i);
			hx->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2);
			Double_t ymax = mp +(rp * 2);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				hx->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_dZAbs_VS_Rrebin_1->SetBinContent(i, mean);
			Resolution_dZAbs_VS_Rrebin_1->SetBinError(i, meanE);
			Resolution_dZAbs_VS_Rrebin_2->SetBinContent(i, sigma);
			Resolution_dZAbs_VS_Rrebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

//********************************************* fitting for dPt vs R *************************************************************
		TH1F *Resolution_R_dPtrebin_1 = new TH1F("Resolution_R_dPtrebin_1", "Mean Photon Resolution dPt vs R ", 88, Rbinning ) ;	
		TH1F *Resolution_R_dPtrebin_2 = new TH1F("Resolution_R_dPtrebin_2", "Sigma Photon Resolution dPt vs R ", 88, Rbinning ) ;
		for (Int_t i = 1; i < 88 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -25,25);
			TH1D *hu = Resolution_R_dPtrebin->ProjectionY("hu", i, i);
			hu->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2);
			Double_t ymax = mp +(rp * 2);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				hu->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_R_dPtrebin_1->SetBinContent(i, mean);
			Resolution_R_dPtrebin_1->SetBinError(i, meanE);
			Resolution_R_dPtrebin_2->SetBinContent(i, sigma);
			Resolution_R_dPtrebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 


		Double_t precisionPhi = 10E-10;
//********************************************* fitting for dPhiAbs vs R *************************************************************
		TH1F *Resolution_dPhiAbs_VS_Rrebin_1 = new TH1F("Resolution_dPhiAbs_VS_Rrebin_1", "Mean Photon Resolution dPhiabs vs R ", 88, Rbinning ) ;	
		TH1F *Resolution_dPhiAbs_VS_Rrebin_2 = new TH1F("Resolution_dPhiAbs_VS_Rrebin_2", "Sigma Photon Resolution dPhiabs vs R ", 88, Rbinning ) ;
		for (Int_t i = 1; i < 88 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -0.5,0.5);
			TH1D *hp = Resolution_dPhiAbs_VS_Rrebin->ProjectionY("hp", i, i);
			hp->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2);
			Double_t ymax = mp +(rp * 2);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precisionPhi && counter < 100){
				f1->SetRange(ymin,ymax);
				hp->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_dPhiAbs_VS_Rrebin_1->SetBinContent(i, mean);
			Resolution_dPhiAbs_VS_Rrebin_1->SetBinError(i, meanE);
			Resolution_dPhiAbs_VS_Rrebin_2->SetBinContent(i, sigma);
			Resolution_dPhiAbs_VS_Rrebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

//************************************************* dPt fitting for ITS bins Positron ************************************************************
		TH1F *Resolution_P_dPt_Pt_ITS0rebin_1 = new TH1F("Resolution_P_dPt_Pt_ITS0rebin_1", "mean Positron Resolution dPt vs Pt 0 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_P_dPt_Pt_ITS0rebin_2 = new TH1F("Resolution_P_dPt_Pt_ITS0rebin_2", "sigma Positron Resolution dPt vs Pt 0 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PITS0 = Resolution_P_dPt_Pt_ITS0rebin->ProjectionY("PITS0", i, i);
			PITS0->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PITS0->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_P_dPt_Pt_ITS0rebin_1->SetBinContent(i, mean);
			Resolution_P_dPt_Pt_ITS0rebin_1->SetBinError(i, meanE);
			Resolution_P_dPt_Pt_ITS0rebin_2->SetBinContent(i, sigma);
			Resolution_P_dPt_Pt_ITS0rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_P_dPt_Pt_ITS1rebin_1 = new TH1F("Resolution_P_dPt_Pt_ITS1rebin_1", "mean Positron Resolution dPt vs Pt 1 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_P_dPt_Pt_ITS1rebin_2 = new TH1F("Resolution_P_dPt_Pt_ITS1rebin_2", "sigma Positron Resolution dPt vs Pt 1 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PITS1 = Resolution_P_dPt_Pt_ITS1rebin->ProjectionY("PITS1", i, i);
			PITS1->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PITS1->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_P_dPt_Pt_ITS1rebin_1->SetBinContent(i, mean);
			Resolution_P_dPt_Pt_ITS1rebin_1->SetBinError(i, meanE);
			Resolution_P_dPt_Pt_ITS1rebin_2->SetBinContent(i, sigma);
			Resolution_P_dPt_Pt_ITS1rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_P_dPt_Pt_ITS2rebin_1 = new TH1F("Resolution_P_dPt_Pt_ITS2rebin_1", "mean Positron Resolution dPt vs Pt 2 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_P_dPt_Pt_ITS2rebin_2 = new TH1F("Resolution_P_dPt_Pt_ITS2rebin_2", "sigma Positron Resolution dPt vs Pt 2 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PITS2 = Resolution_P_dPt_Pt_ITS2rebin->ProjectionY("PITS2", i, i);
			PITS2->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PITS2->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_P_dPt_Pt_ITS2rebin_1->SetBinContent(i, mean);
			Resolution_P_dPt_Pt_ITS2rebin_1->SetBinError(i, meanE);
			Resolution_P_dPt_Pt_ITS2rebin_2->SetBinContent(i, sigma);
			Resolution_P_dPt_Pt_ITS2rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_P_dPt_Pt_ITS3rebin_1 = new TH1F("Resolution_P_dPt_Pt_ITS3rebin_1", "mean Positron Resolution dPt vs Pt 3 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_P_dPt_Pt_ITS3rebin_2 = new TH1F("Resolution_P_dPt_Pt_ITS3rebin_2", "sigma Positron Resolution dPt vs Pt 3 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PITS3 = Resolution_P_dPt_Pt_ITS3rebin->ProjectionY("PITS3", i, i);
			PITS3->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PITS3->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_P_dPt_Pt_ITS3rebin_1->SetBinContent(i, mean);
			Resolution_P_dPt_Pt_ITS3rebin_1->SetBinError(i, meanE);
			Resolution_P_dPt_Pt_ITS3rebin_2->SetBinContent(i, sigma);
			Resolution_P_dPt_Pt_ITS3rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_P_dPt_Pt_ITS4rebin_1 = new TH1F("Resolution_P_dPt_Pt_ITS4rebin_1", "mean Positron Resolution dPt vs Pt 4 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_P_dPt_Pt_ITS4rebin_2 = new TH1F("Resolution_P_dPt_Pt_ITS4rebin_2", "sigma Positron Resolution dPt vs Pt 4 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PITS4 = Resolution_P_dPt_Pt_ITS4rebin->ProjectionY("PITS4", i, i);
			PITS4->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PITS4->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_P_dPt_Pt_ITS4rebin_1->SetBinContent(i, mean);
			Resolution_P_dPt_Pt_ITS4rebin_1->SetBinError(i, meanE);
			Resolution_P_dPt_Pt_ITS4rebin_2->SetBinContent(i, sigma);
			Resolution_P_dPt_Pt_ITS4rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_P_dPt_Pt_ITS5rebin_1 = new TH1F("Resolution_P_dPt_Pt_ITS5rebin_1", "mean Positron Resolution dPt vs Pt 5 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_P_dPt_Pt_ITS5rebin_2 = new TH1F("Resolution_P_dPt_Pt_ITS5rebin_2", "mean Positron Resolution dPt vs Pt 5 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PITS5 = Resolution_P_dPt_Pt_ITS5rebin->ProjectionY("PITS5", i, i);
			PITS5->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PITS5->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_P_dPt_Pt_ITS5rebin_1->SetBinContent(i, mean);
			Resolution_P_dPt_Pt_ITS5rebin_1->SetBinError(i, meanE);
			Resolution_P_dPt_Pt_ITS5rebin_2->SetBinContent(i, sigma);
			Resolution_P_dPt_Pt_ITS5rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_P_dPt_Pt_ITS6rebin_1 = new TH1F("Resolution_P_dPt_Pt_ITS6rebin_1", "mean Positron Resolution dPt vs Pt 6 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_P_dPt_Pt_ITS6rebin_2 = new TH1F("Resolution_P_dPt_Pt_ITS6rebin_2", "sigma Positron Resolution dPt vs Pt 6 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PITS6 = Resolution_P_dPt_Pt_ITS6rebin->ProjectionY("PITS6", i, i);
			PITS6->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PITS6->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_P_dPt_Pt_ITS6rebin_1->SetBinContent(i, mean);
			Resolution_P_dPt_Pt_ITS6rebin_1->SetBinError(i, meanE);
			Resolution_P_dPt_Pt_ITS6rebin_2->SetBinContent(i, sigma);
			Resolution_P_dPt_Pt_ITS6rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

//*************************************** fitting dPt for ITS bins for Electron ***************************************
		TH1F *Resolution_E_dPt_Pt_ITS0rebin_1 = new TH1F("Resolution_E_dPt_Pt_ITS0rebin_1", "mean Electron Resolution dPt vs Pt 0 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_E_dPt_Pt_ITS0rebin_2 = new TH1F("Resolution_E_dPt_Pt_ITS0rebin_2", "sigma Electron Resolution dPt vs Pt 0 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *EITS0 = Resolution_E_dPt_Pt_ITS0rebin->ProjectionY("EITS0", i, i);
			EITS0->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				EITS0->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_E_dPt_Pt_ITS0rebin_1->SetBinContent(i, mean);
			Resolution_E_dPt_Pt_ITS0rebin_1->SetBinError(i, meanE);
			Resolution_E_dPt_Pt_ITS0rebin_2->SetBinContent(i, sigma);
			Resolution_E_dPt_Pt_ITS0rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_E_dPt_Pt_ITS1rebin_1 = new TH1F("Resolution_E_dPt_Pt_ITS1rebin_1", "mean Electron Resolution dPt vs Pt 1 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_E_dPt_Pt_ITS1rebin_2 = new TH1F("Resolution_E_dPt_Pt_ITS1rebin_2", "sigma Electron Resolution dPt vs Pt 1 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *EITS1 = Resolution_E_dPt_Pt_ITS1rebin->ProjectionY("EITS1", i, i);
			EITS1->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				EITS1->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_E_dPt_Pt_ITS1rebin_1->SetBinContent(i, mean);
			Resolution_E_dPt_Pt_ITS1rebin_1->SetBinError(i, meanE);
			Resolution_E_dPt_Pt_ITS1rebin_2->SetBinContent(i, sigma);
			Resolution_E_dPt_Pt_ITS1rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_E_dPt_Pt_ITS2rebin_1 = new TH1F("Resolution_E_dPt_Pt_ITS2rebin_1", "mean Electron Resolution dPt vs Pt 2 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_E_dPt_Pt_ITS2rebin_2 = new TH1F("Resolution_E_dPt_Pt_ITS2rebin_2", "sigma Electron Resolution dPt vs Pt 2 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *EITS2 = Resolution_E_dPt_Pt_ITS2rebin->ProjectionY("EITS2", i, i);
			EITS2->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				EITS2->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_E_dPt_Pt_ITS2rebin_1->SetBinContent(i, mean);
			Resolution_E_dPt_Pt_ITS2rebin_1->SetBinError(i, meanE);
			Resolution_E_dPt_Pt_ITS2rebin_2->SetBinContent(i, sigma);
			Resolution_E_dPt_Pt_ITS2rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_E_dPt_Pt_ITS3rebin_1 = new TH1F("Resolution_E_dPt_Pt_ITS3rebin_1", "mean Electron Resolution dPt vs Pt 3 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_E_dPt_Pt_ITS3rebin_2 = new TH1F("Resolution_E_dPt_Pt_ITS3rebin_2", "sigma Electron Resolution dPt vs Pt 3 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *EITS3 = Resolution_E_dPt_Pt_ITS3rebin->ProjectionY("EITS3", i, i);
			EITS3->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				EITS3->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_E_dPt_Pt_ITS3rebin_1->SetBinContent(i, mean);
			Resolution_E_dPt_Pt_ITS3rebin_1->SetBinError(i, meanE);
			Resolution_E_dPt_Pt_ITS3rebin_2->SetBinContent(i, sigma);
			Resolution_E_dPt_Pt_ITS3rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_E_dPt_Pt_ITS4rebin_1 = new TH1F("Resolution_E_dPt_Pt_ITS4rebin_1", "mean Electron Resolution dPt vs Pt 4 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_E_dPt_Pt_ITS4rebin_2 = new TH1F("Resolution_E_dPt_Pt_ITS4rebin_2", "sigma Electron Resolution dPt vs Pt 4 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *EITS4 = Resolution_E_dPt_Pt_ITS4rebin->ProjectionY("EITS4", i, i);
			EITS4->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				EITS4->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_E_dPt_Pt_ITS4rebin_1->SetBinContent(i, mean);
			Resolution_E_dPt_Pt_ITS4rebin_1->SetBinError(i, meanE);
			Resolution_E_dPt_Pt_ITS4rebin_2->SetBinContent(i, sigma);
			Resolution_E_dPt_Pt_ITS4rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_E_dPt_Pt_ITS5rebin_1 = new TH1F("Resolution_E_dPt_Pt_ITS5rebin_1", "mean Electron Resolution dPt vs Pt 5 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_E_dPt_Pt_ITS5rebin_2 = new TH1F("Resolution_E_dPt_Pt_ITS5rebin_2", "mean Electron Resolution dPt vs Pt 5 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *EITS5 = Resolution_E_dPt_Pt_ITS5rebin->ProjectionY("EITS5", i, i);
			EITS5->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				EITS5->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_E_dPt_Pt_ITS5rebin_1->SetBinContent(i, mean);
			Resolution_E_dPt_Pt_ITS5rebin_1->SetBinError(i, meanE);
			Resolution_E_dPt_Pt_ITS5rebin_2->SetBinContent(i, sigma);
			Resolution_E_dPt_Pt_ITS5rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_E_dPt_Pt_ITS6rebin_1 = new TH1F("Resolution_E_dPt_Pt_ITS6rebin_1", "mean Electron Resolution dPt vs Pt 6 ITS clusters", 30, ptbinning) ;	
		TH1F *Resolution_E_dPt_Pt_ITS6rebin_2 = new TH1F("Resolution_E_dPt_Pt_ITS6rebin_2", "sigma Electron Resolution dPt vs Pt 6 ITS clusters", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *EITS6 = Resolution_E_dPt_Pt_ITS6rebin->ProjectionY("EITS6", i, i);
			EITS6->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				EITS6->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_E_dPt_Pt_ITS6rebin_1->SetBinContent(i, mean);
			Resolution_E_dPt_Pt_ITS6rebin_1->SetBinError(i, meanE);
			Resolution_E_dPt_Pt_ITS6rebin_2->SetBinContent(i, sigma);
			Resolution_E_dPt_Pt_ITS6rebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 


//********************************** Fitting with TPC bins ****************************************

		TH1F *Resolution_E_dPt_Ptrebin_1 = new TH1F("Resolution_E_dPt_Ptrebin_1", "mean Electron Resolution dPt vs Pt", 30, ptbinning) ;	
		TH1F *Resolution_E_dPt_Ptrebin_2 = new TH1F("Resolution_E_dPt_Ptrebin_2", "sigma Electron Resolution dPt vs Pt", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PtbinE = Resolution_E_dPt_Ptrebin->ProjectionY("PtbinE", i, i);
			PtbinE->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PtbinE->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_E_dPt_Ptrebin_1->SetBinContent(i, mean);
			Resolution_E_dPt_Ptrebin_1->SetBinError(i, meanE);
			Resolution_E_dPt_Ptrebin_2->SetBinContent(i, sigma);
			Resolution_E_dPt_Ptrebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_P_dPt_Ptrebin_1 = new TH1F("Resolution_P_dPt_Ptrebin_1", "mean Positron Resolution dPt vs Pt", 30, ptbinning) ;	
		TH1F *Resolution_P_dPt_Ptrebin_2 = new TH1F("Resolution_P_dPt_Ptrebin_2", "sigma Positron Resolution dPt vs Pt", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PtbinP = Resolution_P_dPt_Ptrebin->ProjectionY("PtbinP", i, i);
			PtbinP->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PtbinP->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2.);
				ymax = mp +(rp * 2.);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_P_dPt_Ptrebin_1->SetBinContent(i, mean);
			Resolution_P_dPt_Ptrebin_1->SetBinError(i, meanE);
			Resolution_P_dPt_Ptrebin_2->SetBinContent(i, sigma);
			Resolution_P_dPt_Ptrebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_Gamma_dPt_Ptrebin_1 = new TH1F("Resolution_Gamma_dPt_Ptrebin_1", "mean Gamma Resolution dPt vs Pt", 30, ptbinning) ;	
		TH1F *Resolution_Gamma_dPt_Ptrebin_2 = new TH1F("Resolution_Gamma_dPt_Ptrebin_2", "sigma Gamma Resolution dPt vs Pt", 30, ptbinning) ;	
		for (Int_t i = 1; i < 30 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PtbinGamma = Resolution_Gamma_dPt_Ptrebin->ProjectionY("PtbinGamma", i, i);
			PtbinGamma->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2.);
			Double_t ymax = mp + (rp * 2.);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PtbinGamma->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_Gamma_dPt_Ptrebin_1->SetBinContent(i, mean);
			Resolution_Gamma_dPt_Ptrebin_1->SetBinError(i, meanE);
			Resolution_Gamma_dPt_Ptrebin_2->SetBinContent(i, sigma);
			Resolution_Gamma_dPt_Ptrebin_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 


//******************************** Fitting for Phi binning *******************************************************************

		TH1F *Resolution_Gamma_dPt_Phi_1 = new TH1F("Resolution_Gamma_dPt_Phi_1", "mean Gamma Resolution dPt vs Phi", 100, -TMath::Pi(), TMath::Pi()) ;	
		TH1F *Resolution_Gamma_dPt_Phi_2 = new TH1F("Resolution_Gamma_dPt_Phi_2", "sigma Gamma Resolution dPt vs Phi", 100, -TMath::Pi(), TMath::Pi()) ;	
		for (Int_t i = 1; i < 100 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PhibinGamma = Resolution_Gamma_dPt_Phi->ProjectionY("PhibinGamma", i, i);
			PhibinGamma->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2);
			Double_t ymax = mp +(rp * 2);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PhibinGamma->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_Gamma_dPt_Phi_1->SetBinContent(i, mean);
			Resolution_Gamma_dPt_Phi_1->SetBinError(i, meanE);
			Resolution_Gamma_dPt_Phi_2->SetBinContent(i, sigma);
			Resolution_Gamma_dPt_Phi_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_E_dPt_Phi_1 = new TH1F("Resolution_E_dPt_Phi_1", "mean Electron Resolution dPt vs Phi", 100, -TMath::Pi(), TMath::Pi()) ;	
		TH1F *Resolution_E_dPt_Phi_2 = new TH1F("Resolution_E_dPt_Phi_2", "sigma Electron Resolution dPt vs Phi", 100, -TMath::Pi(), TMath::Pi()) ;	
		for (Int_t i = 1; i < 100 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PhibinE = Resolution_E_dPt_Phi->ProjectionY("PhibinE", i, i);
			PhibinE->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2);
			Double_t ymax = mp +(rp * 2);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PhibinE->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_E_dPt_Phi_1->SetBinContent(i, mean);
			Resolution_E_dPt_Phi_1->SetBinError(i, meanE);
			Resolution_E_dPt_Phi_2->SetBinContent(i, sigma);
			Resolution_E_dPt_Phi_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 

		TH1F *Resolution_P_dPt_Phi_1 = new TH1F("Resolution_P_dPt_Phi_1", "mean Positron Resolution dPt vs Phi", 100, -TMath::Pi(), TMath::Pi()) ;	
		TH1F *Resolution_P_dPt_Phi_2 = new TH1F("Resolution_P_dPt_Phi_2", "sigma Positron Resolution dPt vs Phi", 100, -TMath::Pi(), TMath::Pi()) ;	
		for (Int_t i = 1; i < 100 + 1 ; i ++){
			TF1 *f0 = new TF1("f0", "gaus", -10.,10.);
			TH1D *PhibinP = Resolution_P_dPt_Phi->ProjectionY("PhibinP", i, i);
			PhibinP->Fit(f0,"0RME");
			Double_t rp = f0->GetParameter(2);
			Double_t mp = f0->GetParameter(1);
			Double_t ymin = mp -(rp * 2);
			Double_t ymax = mp +(rp * 2);
			Double_t deviation = 100;
			Int_t counter = 0;
			TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
			while(deviation > precision && counter < 100){
				f1->SetRange(ymin,ymax);
				PhibinP->Fit(f1,"0RME");
				Double_t rp2 = f1->GetParameter(2);
				if (rp2>rp){ deviation = rp2-rp;} 
					else {deviation = rp -rp2 ;}
				rp = rp2 ;
				mp = f1->GetParameter(1);
				ymin = mp -(rp * 2);
				ymax = mp +(rp * 2);
				counter++;
			}
			Double_t mean = f1->GetParameter(1);
			Double_t meanE = f1->GetParError(1);
			Double_t sigma = f1->GetParameter(2);
			Double_t sigmaE = f1->GetParError(2);
			
			Resolution_P_dPt_Phi_1->SetBinContent(i, mean);
			Resolution_P_dPt_Phi_1->SetBinError(i, meanE);
			Resolution_P_dPt_Phi_2->SetBinContent(i, sigma);
			Resolution_P_dPt_Phi_2->SetBinError(i, sigmaE);
			delete f1;
			delete f0;
		} //end of fitting 



//********************************************************************************************************
//*************************************** creating ps-file ***********************************************
//********************************************************************************************************
if(!SinglePlots) { 	
		TPostScript *ps_characteristics;
		ps_characteristics = new TPostScript(Form("%s%s.ps",path,output),111);	

		//******************************* Performance *********************************************

		ps_characteristics->NewPage();
		
		TCanvas * c5 = new TCanvas("c5","",700,1000);  // gives the page size
		pad_c5 = new TPad("pad_c4","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c5->SetFillColor(0);
		pad_c5->GetFrame()->SetFillColor(0);
		pad_c5->SetBorderMode(0);
		pad_c5->Divide(2,2);
		pad_c5->Draw();


		leg1 = new TLegend( 0.6,0.82,0.92,0.9);
		leg1->SetTextSize(0.03);			
		leg1->SetFillColor(0);
		leg1->AddEntry(Resolution_MC_Pt,("MC truth"));
		leg1->AddEntry(Resolution_ESD_Pt,("MC reconstructed"));
			pad_c5->cd(1)->SetLogy(1);
			pad_c5->cd(1);
				DrawAutoGammaHistos( Resolution_MC_Pt, 
								 Resolution_ESD_Pt, 
								 "p_{t} distribution of the #gamma","p_{t} ","N_{#gamma}",
								 kTRUE, 1.2,1,
								 kFALSE,0. ,0.,
								 kFALSE, 50.,100.);
				leg1->Draw();	

			pad_c5->cd(1)->Update();	

			pad_c5->cd(3);
				DrawAutoGammaHistos( Resolution_MC_Z, 
								 Resolution_ESD_Z, 
								 "Z distribution of conversion point","Z [cm] ","N_{#gamma}",
								 kTRUE, 1.2,0,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,200.);
				leg1->Draw();	

			pad_c5->cd(2)->Update();	
		
			pad_c5->cd(4);
				DrawAutoGammaHistos( Resolution_MC_R, 
								 Resolution_ESD_R, 
								 "R distribution of conversion point","R [cm] ","N_{#gamma}",
								 kTRUE, 1.2,0,
								 kFALSE,0. ,0.,
								 kTRUE, 0.,200.);
				leg1->Draw();	

			pad_c5->cd(3)->Update();			
		c5->Update();

		//************************************ Resolution Gamma dPt vs Pt **********************************
		ps_characteristics->NewPage();
		TCanvas * c0 = new TCanvas("c0","",700,1000);  // gives the page size
		pad_c0 = new TPad("pad_c0","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c0->SetFillColor(0);
		pad_c0->GetFrame()->SetFillColor(0);
		pad_c0->SetBorderMode(0);
		pad_c0->Divide(2,2);
		pad_c0->Draw();

		pad_c0->cd(1)->SetLogz(1);
		pad_c0->cd(1)->SetLogx(1);
		pad_c0->cd(2)->SetLogz(1);
		pad_c0->cd(2)->SetLogx(1);

		pad_c0->cd(1)->SetRightMargin(0.17); 		
		pad_c0->cd(1);
			DrawAutoGammaHisto2D(	Resolution_Gamma_dPt_Pt,
								"", "p_{t} [GeV/c]", "dp_{t}/p_{t,MC} [%]", "",
								kFALSE, 10., 140.,
								kFALSE, 0., 100.);
		pad_c0->cd(1)->Update();	

		pad_c0->cd(2)->SetRightMargin(0.17); 		
		pad_c0->cd(2);
			DrawAutoGammaHisto2D(	Resolution_Gamma_dPt_Ptrebin,
								"", "p_{t} [GeV/c]", "dp_{t}/p_{t,MC} [%]", "",
								kFALSE, 10., 140.,
								kTRUE, 0., 100.);
		pad_c0->cd(2)->Update();	


		pad_c0->cd(3);
			StylingSliceHistos(Resolution_Gamma_dPt_Ptrebin_1,0.8);
			Resolution_Gamma_dPt_Ptrebin_1->SetMarkerColor(kRed+2);
			Resolution_Gamma_dPt_Ptrebin_1->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_Gamma_dPt_Ptrebin_1, 
								"", "p_{t} [GeV/c]", "Peak pos. dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, -2, 4.,
								kTRUE, 0.8, 100.);
		pad_c0->cd(3)->Update();

		pad_c0->cd(4);
			StylingSliceHistos(Resolution_Gamma_dPt_Ptrebin_2,0.8);
			Resolution_Gamma_dPt_Ptrebin_2->SetMarkerColor(kRed+2);
			Resolution_Gamma_dPt_Ptrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_Gamma_dPt_Ptrebin_2, 
							"", "p_{t} [GeV/c]", " #sigma dp_{t}/p_{t,MC} [%] ", 
								kFALSE, 10., 140.,
								kTRUE, 0, 16.,
								kTRUE, 0.8, 100.);
		pad_c0->cd(4)->Update();

		c0->Update();
		

		//************************************ Resolution E dPt vs Pt **********************************
		ps_characteristics->NewPage();
		
		TCanvas * c1 = new TCanvas("c1","",700,1000);  // gives the page size
		pad_c1 = new TPad("pad_c1","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c1->SetFillColor(0);
		pad_c1->GetFrame()->SetFillColor(0);
		pad_c1->SetBorderMode(0);
		pad_c1->Divide(2,2);
		pad_c1->Draw();

		pad_c1->cd(1)->SetLogz(1);
		pad_c1->cd(1)->SetLogx(1);
		pad_c1->cd(2)->SetLogz(1);
		pad_c1->cd(2)->SetLogx(1);

		pad_c1->cd(1)->SetRightMargin(0.17); 		
		pad_c1->cd(1);
			DrawAutoGammaHisto2D(	Resolution_E_dPt_Pt,
								"", "p_{t} [GeV/c]", "dp_{t}/p_{t,MC} [%]", "",
								kFALSE, 10., 140.,
								kFALSE, 0., 100.);
		pad_c1->cd(1)->Update();	

		pad_c1->cd(2)->SetRightMargin(0.17); 		
		pad_c1->cd(2);
			DrawAutoGammaHisto2D(	Resolution_E_dPt_Ptrebin,
								"", "p_{t} [GeV/c]", "dp_{t}/p_{t,MC} [%]", "",
								kFALSE, 10., 140.,
								kTRUE, 0., 100.);
		pad_c1->cd(2)->Update();	


		pad_c1->cd(3);
			StylingSliceHistos(Resolution_E_dPt_Ptrebin_1,0.8);
			Resolution_E_dPt_Ptrebin_1->SetMarkerColor(kRed+2);
			Resolution_E_dPt_Ptrebin_1->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_E_dPt_Ptrebin_1, 
								"", "p_{t} [GeV/c]", "Peak pos. dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, -2, 2.,
								kTRUE, 0.8, 100.);
		pad_c1->cd(3)->Update();

		pad_c1->cd(4);
			StylingSliceHistos(Resolution_E_dPt_Ptrebin_2,0.8);
			Resolution_E_dPt_Ptrebin_2->SetMarkerColor(kRed+2);
			Resolution_E_dPt_Ptrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_E_dPt_Ptrebin_2, 
							"", "p_{t} [GeV/c]", " #sigma dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, 0, 15.,
								kTRUE, 0.8, 100.);
		pad_c1->cd(4)->Update();

		c1->Update();

		//******************************* Resoltution E dPt vs Pt ITS dependent ***************************

		ps_characteristics->NewPage();
		
		TCanvas * c2 = new TCanvas("c2","",700,1000);  // gives the page size
		pad_c2 = new TPad("pad_c2","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c2->SetFillColor(0);
		pad_c2->GetFrame()->SetFillColor(0);
		pad_c2->SetBorderMode(0);
		pad_c2->Divide(1,2);
		pad_c2->Draw();


		leg2 = new TLegend( 0.6,0.72,0.92,0.9);
		leg2->SetTextSize(0.03);			
		leg2->SetFillColor(0);			pad_c5->cd(1);
		leg2->AddEntry(Resolution_E_dPt_Pt_ITS0rebin_1,("only TPC"));
//		leg2->AddEntry(Resolution_E_dPt_Pt_ITS1rebin_1,("1 ITS cluster"));
		leg2->AddEntry(Resolution_E_dPt_Pt_ITS2rebin_1,("2 ITS cluster"));
//		leg2->AddEntry(Resolution_E_dPt_Pt_ITS3rebin_1,("3 ITS cluster"));
		leg2->AddEntry(Resolution_E_dPt_Pt_ITS4rebin_1,("4 ITS cluster"));
		leg2->AddEntry(Resolution_E_dPt_Pt_ITS5rebin_1,("5 ITS cluster"));
		leg2->AddEntry(Resolution_E_dPt_Pt_ITS6rebin_1,("6 ITS cluster"));		

		pad_c2->cd(1);
			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS0rebin_1,0.8);
			Resolution_E_dPt_Pt_ITS0rebin_1->SetMarkerColor(kRed+2);
			Resolution_E_dPt_Pt_ITS0rebin_1->SetLineColor(kRed+2);
			DrawResolutionGammaHisto( Resolution_E_dPt_Pt_ITS0rebin_1, 
								"", "p_{t} [GeV/c]", "Peak pos. dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, -20, 20.,
								kTRUE, 0, 18.);
/*			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS1rebin_1,1);
			Resolution_E_dPt_Pt_ITS1rebin_1->SetMarkerColor(kBlue+2);
			Resolution_E_dPt_Pt_ITS1rebin_1->SetLineColor(kBlue+2);
			Resolution_E_dPt_Pt_ITS1rebin_1->Draw("same");
*/			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS2rebin_1,0.8);
			Resolution_E_dPt_Pt_ITS2rebin_1->SetMarkerColor(kViolet+2);
			Resolution_E_dPt_Pt_ITS2rebin_1->SetLineColor(kViolet+2);
			Resolution_E_dPt_Pt_ITS2rebin_1->Draw("same");
/*			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS3rebin_1,1);
			Resolution_E_dPt_Pt_ITS3rebin_1->SetMarkerColor(kOrange+2);
			Resolution_E_dPt_Pt_ITS3rebin_1->SetLineColor(kOrange+2);
			Resolution_E_dPt_Pt_ITS3rebin_1->Draw("same");
*/			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS4rebin_1,0.8);
			Resolution_E_dPt_Pt_ITS4rebin_1->SetMarkerColor(kCyan+2);
			Resolution_E_dPt_Pt_ITS4rebin_1->SetLineColor(kCyan+2);
			Resolution_E_dPt_Pt_ITS4rebin_1->Draw("same");
			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS5rebin_1,1);
			Resolution_E_dPt_Pt_ITS5rebin_1->SetMarkerColor(kGreen+2);
			Resolution_E_dPt_Pt_ITS5rebin_1->SetLineColor(kGreen+2);
			Resolution_E_dPt_Pt_ITS5rebin_1->Draw("same");
			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS6rebin_1,0.8);
			Resolution_E_dPt_Pt_ITS6rebin_1->SetMarkerColor(kMagenta+2);
			Resolution_E_dPt_Pt_ITS6rebin_1->SetLineColor(kMagenta+2);
			Resolution_E_dPt_Pt_ITS6rebin_1->Draw("same");
			leg2->Draw("same");
		pad_c2->cd(1)->Update();	

		pad_c2->cd(2);
			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS0rebin_2,0.8);
			Resolution_E_dPt_Pt_ITS0rebin_2->SetMarkerColor(kRed+2);
			Resolution_E_dPt_Pt_ITS0rebin_2->SetLineColor(kRed+2);
			DrawResolutionGammaHisto( Resolution_E_dPt_Pt_ITS0rebin_2, 
								"", "p_{t} [GeV/c]", "#sigma dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, -1, 30.,
								kTRUE, 0, 18.);
/*			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS1rebin_2,1);
			Resolution_E_dPt_Pt_ITS1rebin_2->SetMarkerColor(kBlue+2);
			Resolution_E_dPt_Pt_ITS1rebin_2->SetLineColor(kBlue+2);
			Resolution_E_dPt_Pt_ITS1rebin_2->Draw("same");
*/			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS2rebin_2,0.8);
			Resolution_E_dPt_Pt_ITS2rebin_2->SetMarkerColor(kViolet+2);
			Resolution_E_dPt_Pt_ITS2rebin_2->SetLineColor(kViolet+2);
			Resolution_E_dPt_Pt_ITS2rebin_2->Draw("same");
/*			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS3rebin_2,1);
			Resolution_E_dPt_Pt_ITS3rebin_2->SetMarkerColor(kOrange+2);
			Resolution_E_dPt_Pt_ITS3rebin_2->SetLineColor(kOrange+2);
			Resolution_E_dPt_Pt_ITS3rebin_2->Draw("same");
*/			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS4rebin_2,0.8);
			Resolution_E_dPt_Pt_ITS4rebin_2->SetMarkerColor(kCyan+2);
			Resolution_E_dPt_Pt_ITS4rebin_2->SetLineColor(kCyan+2);
			Resolution_E_dPt_Pt_ITS4rebin_2->Draw("same");
			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS5rebin_2,0.8);
			Resolution_E_dPt_Pt_ITS5rebin_2->SetMarkerColor(kGreen+2);
			Resolution_E_dPt_Pt_ITS5rebin_2->SetLineColor(kGreen+2);
			Resolution_E_dPt_Pt_ITS5rebin_2->Draw("same");
			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS6rebin_2,0.8);
			Resolution_E_dPt_Pt_ITS6rebin_2->SetMarkerColor(kMagenta+2);
			Resolution_E_dPt_Pt_ITS6rebin_2->SetLineColor(kMagenta+2);
			Resolution_E_dPt_Pt_ITS6rebin_2->Draw("same");
			leg2->Draw("same");
		pad_c2->cd(2)->Update();	


		c2->Update();

		//************************************ Resolution P dPt vs Pt **********************************
		
		ps_characteristics->NewPage();
		
		TCanvas * c3 = new TCanvas("c3","",700,1000);  // gives the page size
		pad_c3 = new TPad("pad_c3","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c3->SetFillColor(0);
		pad_c3->GetFrame()->SetFillColor(0);
		pad_c3->SetBorderMode(0);
		pad_c3->Divide(2,2);
		pad_c3->Draw();

		pad_c3->cd(1)->SetLogz(1);
		pad_c3->cd(1)->SetLogx(1);
		pad_c3->cd(2)->SetLogz(1);
		pad_c3->cd(2)->SetLogx(1);

		pad_c3->cd(1)->SetRightMargin(0.17); 		
		pad_c3->cd(1);
			DrawAutoGammaHisto2D(	Resolution_P_dPt_Pt,
								"", "p_{t} [GeV/c]", "dp_{t}/p_{t,MC} [%]", "",
								kFALSE, 10., 140.,
								kFALSE, 0., 100.);
		pad_c3->cd(1)->Update();	

		pad_c3->cd(2)->SetRightMargin(0.17); 		
		pad_c3->cd(2);
			DrawAutoGammaHisto2D(	Resolution_P_dPt_Ptrebin,
								"", "p_{t} [GeV/c]", "dp_{t}/p_{t,MC} [%]", "",
								kFALSE, 10., 140.,
								kTRUE, 0., 100.);
		pad_c3->cd(2)->Update();	


		pad_c3->cd(3);
			StylingSliceHistos(Resolution_P_dPt_Ptrebin_1,0.8);
			Resolution_P_dPt_Ptrebin_1->SetMarkerColor(kRed+2);
			Resolution_P_dPt_Ptrebin_1->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_P_dPt_Ptrebin_1, 
								"", "p_{t} [GeV/c]", "Peak pos. dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, -2, 4.,
								kTRUE, 0.8, 100.);
		pad_c3->cd(3)->Update();

		pad_c3->cd(4);
			StylingSliceHistos(Resolution_P_dPt_Ptrebin_2,0.8);
			Resolution_P_dPt_Ptrebin_2->SetMarkerColor(kRed+2);
			Resolution_P_dPt_Ptrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_P_dPt_Ptrebin_2, 
							"", "p_{t} [GeV/c]", " #sigma dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, 0, 16.,
								kTRUE, 0.8, 100.);
		pad_c3->cd(4)->Update();

		c3->Update();

		//******************************* Resoltution P dPt vs Pt ITS dependent ***************************

		ps_characteristics->NewPage();
		
		TCanvas * c4 = new TCanvas("c4","",700,1000);  // gives the page size
		pad_c4 = new TPad("pad_c4","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c4->SetFillColor(0);
		pad_c4->GetFrame()->SetFillColor(0);
		pad_c4->SetBorderMode(0);
		pad_c4->Divide(1,2);
		pad_c4->Draw();

		pad_c4->cd(1);
			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS0rebin_1,0.8);
			Resolution_P_dPt_Pt_ITS0rebin_1->SetMarkerColor(kRed+2);
			Resolution_P_dPt_Pt_ITS0rebin_1->SetLineColor(kRed+2);
			DrawResolutionGammaHisto( Resolution_P_dPt_Pt_ITS0rebin_1, 
								"", "p_{t} [GeV/c]", "Peak pos. dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, -20, 20.,
								kTRUE, 0., 18.);
/*			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS1rebin_1,1);
			Resolution_P_dPt_Pt_ITS1rebin_1->SetMarkerColor(kBlue+2);
			Resolution_P_dPt_Pt_ITS1rebin_1->SetLineColor(kBlue+2);
			Resolution_P_dPt_Pt_ITS1rebin_1->Draw("same");
*/			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS2rebin_1,0.8);
			Resolution_P_dPt_Pt_ITS2rebin_1->SetMarkerColor(kViolet+2);
			Resolution_P_dPt_Pt_ITS2rebin_1->SetLineColor(kViolet+2);
			Resolution_P_dPt_Pt_ITS2rebin_1->Draw("same");
/*			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS3rebin_1,1);
			Resolution_P_dPt_Pt_ITS3rebin_1->SetMarkerColor(kOrange+2);
			Resolution_P_dPt_Pt_ITS3rebin_1->SetLineColor(kOrange+2);
			Resolution_P_dPt_Pt_ITS3rebin_1->Draw("same");
*/			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS4rebin_1,0.8);
			Resolution_P_dPt_Pt_ITS4rebin_1->SetMarkerColor(kCyan+2);
			Resolution_P_dPt_Pt_ITS4rebin_1->SetLineColor(kCyan+2);
			Resolution_P_dPt_Pt_ITS4rebin_1->Draw("same");
			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS5rebin_1,0.8);
			Resolution_P_dPt_Pt_ITS5rebin_1->SetMarkerColor(kGreen+2);
			Resolution_P_dPt_Pt_ITS5rebin_1->SetLineColor(kGreen+2);
			Resolution_P_dPt_Pt_ITS5rebin_1->Draw("same");
			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS6rebin_1,0.8);
			Resolution_P_dPt_Pt_ITS6rebin_1->SetMarkerColor(kMagenta+2);
			Resolution_P_dPt_Pt_ITS6rebin_1->SetLineColor(kMagenta+2);
			Resolution_P_dPt_Pt_ITS6rebin_1->Draw("same");
			leg2->Draw("same");
		pad_c4->cd(1)->Update();	

		pad_c4->cd(2);
			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS0rebin_2,0.8);
			Resolution_P_dPt_Pt_ITS0rebin_2->SetMarkerColor(kRed+2);
			Resolution_P_dPt_Pt_ITS0rebin_2->SetLineColor(kRed+2);
			DrawResolutionGammaHisto( Resolution_P_dPt_Pt_ITS0rebin_2, 
								"", "p_{t} [GeV/c]", "#sigma dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, -1, 30.,
								kTRUE, 0., 18.);
/*			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS1rebin_2,1);
			Resolution_P_dPt_Pt_ITS1rebin_2->SetMarkerColor(kBlue+2);
			Resolution_P_dPt_Pt_ITS1rebin_2->SetLineColor(kBlue+2);
			Resolution_P_dPt_Pt_ITS1rebin_2->Draw("same");
*/			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS2rebin_2,0.8);
			Resolution_P_dPt_Pt_ITS2rebin_2->SetMarkerColor(kViolet+2);
			Resolution_P_dPt_Pt_ITS2rebin_2->SetLineColor(kViolet+2);
			Resolution_P_dPt_Pt_ITS2rebin_2->Draw("same");
/*			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS3rebin_2,1);
			Resolution_P_dPt_Pt_ITS3rebin_2->SetMarkerColor(kOrange+2);
			Resolution_P_dPt_Pt_ITS3rebin_2->SetLineColor(kOrange+2);
			Resolution_P_dPt_Pt_ITS3rebin_2->Draw("same");
*/			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS4rebin_2,0.8);
			Resolution_P_dPt_Pt_ITS4rebin_2->SetMarkerColor(kCyan+2);
			Resolution_P_dPt_Pt_ITS4rebin_2->SetLineColor(kCyan+2);
			Resolution_P_dPt_Pt_ITS4rebin_2->Draw("same");
			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS5rebin_2,1);
			Resolution_P_dPt_Pt_ITS5rebin_2->SetMarkerColor(kGreen+2);
			Resolution_P_dPt_Pt_ITS5rebin_2->SetLineColor(kGreen+2);
			Resolution_P_dPt_Pt_ITS5rebin_2->Draw("same");
			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS6rebin_2,0.8);
			Resolution_P_dPt_Pt_ITS6rebin_2->SetMarkerColor(kMagenta+2);
			Resolution_P_dPt_Pt_ITS6rebin_2->SetLineColor(kMagenta+2);
			Resolution_P_dPt_Pt_ITS6rebin_2->Draw("same");
			leg2->Draw("same");
		pad_c4->cd(2)->Update();	
		c4->Update();


		//************************************ Resolution dZAbs vs R **********************************
		
		ps_characteristics->NewPage();
		
		TCanvas * c6 = new TCanvas("c6","",700,1000);  // gives the page size
		pad_c6 = new TPad("pad_c6","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c6->SetFillColor(0);
		pad_c6->GetFrame()->SetFillColor(0);
		pad_c6->SetBorderMode(0);
		pad_c6->Divide(2,2);
		pad_c6->Draw();

		pad_c6->cd(1)->SetLogz(1);
		pad_c6->cd(1)->SetLogx(1);
		pad_c6->cd(2)->SetLogz(1);
		pad_c6->cd(2)->SetLogx(1);

		pad_c6->cd(1)->SetRightMargin(0.17); 		
		pad_c6->cd(1);
			DrawAutoGammaHisto2D(	Resolution_dZAbs_VS_R,
								"", "R [cm]", "dZ [cm]", "",
								kFALSE, 10., 140.,
								kTRUE, 2.5, 180.);
		pad_c6->cd(1)->Update();	

		pad_c6->cd(2)->SetRightMargin(0.17); 		
		pad_c6->cd(2);
			DrawAutoGammaHisto2D(	Resolution_dZAbs_VS_Rrebin,
								"", "R [cm]", "dZ [cm]", "",
								kFALSE, 10., 140.,
								kFALSE, 0., 180.);
		pad_c6->cd(2)->Update();	


		pad_c6->cd(3);
			StylingSliceHistos(Resolution_dZAbs_VS_Rrebin_1,0.8);
			Resolution_dZAbs_VS_Rrebin_1->SetMarkerColor(kRed+2);
			Resolution_dZAbs_VS_Rrebin_1->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_dZAbs_VS_Rrebin_1, 
								"", "R [cm]", "Peak pos. dZ [cm]", 
								kFALSE, 10., 140.,
								kFALSE, -0.2, 0.2,
								kTRUE, 0.8, 120.);
		pad_c6->cd(3)->Update();

		pad_c6->cd(4);
			StylingSliceHistos(Resolution_dZAbs_VS_Rrebin_2,0.8);
			Resolution_dZAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
			Resolution_dZAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_dZAbs_VS_Rrebin_2, 
							"", "R [cm]", "#sigma dZ [cm]", 
								kFALSE, 10., 140.,
								kTRUE, 0., 3.5,
								kTRUE, 0.8, 120.);
		pad_c6->cd(4)->Update();

		c6->Update();

		//************************************ Resolution dRAbs vs R **********************************
		
		ps_characteristics->NewPage();
		
		TCanvas * c7 = new TCanvas("c7","",700,1000);  // gives the page size
		pad_c7 = new TPad("pad_c7","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c7->SetFillColor(0);
		pad_c7->GetFrame()->SetFillColor(0);
		pad_c7->SetBorderMode(0);
		pad_c7->Divide(2,2);
		pad_c7->Draw();

		pad_c7->cd(1)->SetLogz(1);
		pad_c7->cd(1)->SetLogx(1);
		pad_c7->cd(2)->SetLogz(1);
		pad_c7->cd(2)->SetLogx(1);

		pad_c7->cd(1)->SetRightMargin(0.17); 		
		pad_c7->cd(1);
			DrawAutoGammaHisto2D(	Resolution_dRAbs_VS_R,
								"", "R [cm]", "dR [cm]", "",
								kFALSE, 10., 140.,
								kTRUE, 2.5, 180.);
		pad_c7->cd(1)->Update();	

		pad_c7->cd(2)->SetRightMargin(0.17); 		
		pad_c7->cd(2);
			DrawAutoGammaHisto2D(	Resolution_dRAbs_VS_Rrebin,
								"", "R [cm]", "dR [cm]", "",
								kFALSE, 10., 140.,
								kFALSE, 0., 180.);
		pad_c7->cd(2)->Update();	


		pad_c7->cd(3);
			StylingSliceHistos(Resolution_dRAbs_VS_Rrebin_1,0.8);
			Resolution_dRAbs_VS_Rrebin_1->SetMarkerColor(kRed+2);
			Resolution_dRAbs_VS_Rrebin_1->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_dRAbs_VS_Rrebin_1, 
								"", "R [cm]", "Peak pos. dR [cm]", 
								kFALSE, 10., 140.,
								kFALSE, -1, 1.5,
								kTRUE, 0.8, 120.);
		pad_c7->cd(3)->Update();

		pad_c7->cd(4);
			StylingSliceHistos(Resolution_dRAbs_VS_Rrebin_2,0.8);
			Resolution_dRAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
			Resolution_dRAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_dRAbs_VS_Rrebin_2, 
							"", "R [cm]", "#sigma dR [cm]", 
								kFALSE, 10., 140.,
								kTRUE, 0., 3.5,
								kTRUE, 0.8, 120.);
		pad_c7->cd(4)->Update();

		c7->Update();

		//************************************ Resolution dPhiAbs vs R **********************************
		
		ps_characteristics->NewPage();
		
		TCanvas * c8 = new TCanvas("c8","",700,1000);  // gives the page size
		pad_c8 = new TPad("pad_c8","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c8->SetFillColor(0);
		pad_c8->GetFrame()->SetFillColor(0);
		pad_c8->SetBorderMode(0);
		pad_c8->Divide(2,2);
		pad_c8->Draw();

		pad_c8->cd(1)->SetLogz(1);
		pad_c8->cd(1)->SetLogx(1);
		pad_c8->cd(2)->SetLogz(1);
		pad_c8->cd(2)->SetLogx(1);

		pad_c8->cd(1)->SetRightMargin(0.17); 		
		pad_c8->cd(1);
			DrawAutoGammaHisto2D(	Resolution_dPhiAbs_VS_R,
								"", "R [cm]", "d#phi[rad]", "",
								kFALSE, 10., 140.,
								kTRUE, 2.5, 180.);
		pad_c8->cd(1)->Update();	

		pad_c8->cd(2)->SetRightMargin(0.17); 		
		pad_c8->cd(2);
			DrawAutoGammaHisto2D(	Resolution_dPhiAbs_VS_Rrebin,
								"", "R [cm]", "d#phi[rad]", "",
								kFALSE, 10., 140.,
								kTRUE, 0., 180.);
		pad_c8->cd(2)->Update();	


		pad_c8->cd(3);
			StylingSliceHistos(Resolution_dPhiAbs_VS_Rrebin_1,0.8);
			Resolution_dPhiAbs_VS_Rrebin_1->SetMarkerColor(kRed+2);
			Resolution_dPhiAbs_VS_Rrebin_1->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Rrebin_1, 
								"", "R [cm]", "Peak pos. d#phi[rad]", 
								kFALSE, 10., 140.,
								kFALSE, -0.0005, 0.001,
								kTRUE, 0.8, 120.);
		pad_c8->cd(3)->Update();

		pad_c8->cd(4);
			StylingSliceHistos(Resolution_dPhiAbs_VS_Rrebin_2,0.8);
			Resolution_dPhiAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
			Resolution_dPhiAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Rrebin_2, 
							"", "R [cm]", "#sigma d#phi[rad]", 
								kFALSE, 10., 140.,
								kTRUE, 0.0, 0.01,
								kTRUE, 0.8, 120.);
		pad_c8->cd(4)->Update();

		c8->Update();

		//************************************ Resolution dpt vs R **********************************
		
		ps_characteristics->NewPage();
		
		TCanvas * c9 = new TCanvas("c9","",700,1000);  // gives the page size
		pad_c9 = new TPad("pad_c9","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c9->SetFillColor(0);
		pad_c9->GetFrame()->SetFillColor(0);
		pad_c9->SetBorderMode(0);
		pad_c9->Divide(2,2);
		pad_c9->Draw();

		pad_c9->cd(1)->SetLogz(1);
		pad_c9->cd(1)->SetLogx(1);
		pad_c9->cd(2)->SetLogz(1);
		pad_c9->cd(2)->SetLogx(1);

		pad_c9->cd(1)->SetRightMargin(0.17); 		
		pad_c9->cd(1);
			DrawAutoGammaHisto2D(	Resolution_R_dPt,
								"", "R [cm]", "dp_{t}/p_{t,MC} [%]", "",
								kFALSE, 10., 140.,
								kTRUE, 2.5, 180.);
		pad_c9->cd(1)->Update();	

		pad_c9->cd(2)->SetRightMargin(0.17); 		
		pad_c9->cd(2);
			DrawAutoGammaHisto2D(	Resolution_R_dPtrebin,
								"", "R [cm]", "dp_{t}/p_{t,MC} [%]", "",
								kFALSE, 10., 140.,
								kTRUE, 0., 180.);
		pad_c9->cd(2)->Update();	


		pad_c9->cd(3);
			StylingSliceHistos(Resolution_R_dPtrebin_1,0.8);
			Resolution_R_dPtrebin_1->SetMarkerColor(kRed+2);
			Resolution_R_dPtrebin_1->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_R_dPtrebin_1, 
								"", "R [cm]", "Peak pos. dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kFALSE, -0.001, 0.001,
								kTRUE, 0.8, 120.);
		pad_c9->cd(3)->Update();

		pad_c9->cd(4);
			StylingSliceHistos(Resolution_R_dPtrebin_2,0.8);
			Resolution_R_dPtrebin_2->SetMarkerColor(kRed+2);
			Resolution_R_dPtrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_R_dPtrebin_2, 
							"", "R [cm]", "#sigma dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, 0., 5,
								kTRUE, 0.8, 120.);
		pad_c9->cd(4)->Update();

		c9->Update();

//************************* Fitting for dPt in phi ********************************************************
		ps_characteristics->NewPage();
		
		TCanvas * c10 = new TCanvas("c10","",700,1000);  // gives the page size
		pad_c10 = new TPad("pad_c10","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c10->SetFillColor(0);
		pad_c10->GetFrame()->SetFillColor(0);
		pad_c10->SetBorderMode(0);
		pad_c10->Divide(2,3);
		pad_c10->Draw();


		pad_c10->cd(1);
			StylingSliceHistos(Resolution_Gamma_dPt_Phi_1,0.8);
			Resolution_Gamma_dPt_Phi_1->SetMarkerColor(kRed+2);
			Resolution_Gamma_dPt_Phi_1->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_Gamma_dPt_Phi_1, 
								"", "#phi [rad]", "Peak pos. dp_{t}/p_{t,MC} [%] Photon", 
								kFALSE, 10., 140.,
								kFALSE, 0., 1.8,
								kFALSE, 0.8, 180.);
		pad_c10->cd(1)->Update();	

		pad_c10->cd(2);
			StylingSliceHistos(Resolution_Gamma_dPt_Phi_2,0.8);
			Resolution_Gamma_dPt_Phi_2->SetMarkerColor(kRed+2);
			Resolution_Gamma_dPt_Phi_2->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_Gamma_dPt_Phi_2, 
								"", "#phi [rad]", "#sigma dp_{t}/p_{t,MC} [%] Photon", 
								kFALSE, 10., 140.,
								kTRUE, 0., 1.8,
								kFALSE, 0.8, 180.);		
		pad_c10->cd(2)->Update();	

		pad_c10->cd(3);
			StylingSliceHistos(Resolution_E_dPt_Phi_1,0.8);
			Resolution_E_dPt_Phi_1->SetMarkerColor(kRed+2);
			Resolution_E_dPt_Phi_1->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_E_dPt_Phi_1, 
								"", "#phi [rad]", "Peak pos. dp_{t}/p_{t,MC} [%] Electron", 
								kFALSE, 10., 140.,
								kFALSE, -0.001, 0.001,
								kFALSE, 0.8, 180.);
		pad_c10->cd(3)->Update();	

		pad_c10->cd(4);
			StylingSliceHistos(Resolution_E_dPt_Phi_2,0.8);
			Resolution_E_dPt_Phi_2->SetMarkerColor(kRed+2);
			Resolution_E_dPt_Phi_2->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_E_dPt_Phi_2, 
								"", "#phi [rad]", "#sigma dp_{t}/p_{t,MC} [%] Electron", 
								kFALSE, 10., 140.,
								kTRUE, 0., 1.8,
								kFALSE, 0.8, 180.);		
		pad_c10->cd(4)->Update();	

		pad_c10->cd(5);
			StylingSliceHistos(Resolution_P_dPt_Phi_1,0.8);
			Resolution_P_dPt_Phi_1->SetMarkerColor(kRed+2);
			Resolution_P_dPt_Phi_1->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_P_dPt_Phi_1, 
								"", "#phi [rad]", "Peak pos. dp_{t}/p_{t,MC} [%] Positron", 
								kFALSE, 10., 140.,
								kFALSE, -0.001, 0.001,
								kFALSE, 0.8, 180.);
		pad_c10->cd(5)->Update();	

		pad_c10->cd(6);
			StylingSliceHistos(Resolution_P_dPt_Phi_2,0.8);
			Resolution_P_dPt_Phi_2->SetMarkerColor(kRed+2);
			Resolution_P_dPt_Phi_2->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_P_dPt_Phi_2, 
								"", "#phi [rad]", "#sigma dp_{t}/p_{t,MC} [%] Positron", 
								kFALSE, 10., 140.,
								kTRUE, 0., 1.8,
								kFALSE, 0.8, 180.);		
		pad_c10->cd(6)->Update();	
		c10->Update();

		ps_characteristics->NewPage();

//******************************** TRD investigation ******************************************************		
		TCanvas * c11 = new TCanvas("c11","",700,1000);  // gives the page size
		pad_c11 = new TPad("pad_c11","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c11->SetFillColor(0);
		pad_c11->GetFrame()->SetFillColor(0);
		pad_c11->SetBorderMode(0);
		pad_c11->Divide(2,2);
		pad_c11->Draw();

		pad_c11->cd(1)->SetLogz(1);
		pad_c11->cd(1)->SetLogx(1);
		pad_c11->cd(2)->SetLogz(1);
		pad_c11->cd(2)->SetLogx(1);
		pad_c11->cd(3)->SetLogz(1);
		pad_c11->cd(3)->SetLogx(1);
		pad_c11->cd(4)->SetLogz(1);
		pad_c11->cd(4)->SetLogx(1);


		pad_c11->cd(1)->SetRightMargin(0.17); 		
		pad_c11->cd(1);
			DrawAutoGammaHisto2D(	Resolution_E_nTRDclusters_ESDPt,
								"", "p_{t} [GeV/c]", "N_{clusters} Electron", "",
								kFALSE, 10., 140.,
								kTRUE, 0.1, 180.);
		pad_c11->cd(1)->Update();	

		pad_c11->cd(2)->SetRightMargin(0.17); 		
		pad_c11->cd(2);
			DrawAutoGammaHisto2D(	Resolution_P_nTRDclusters_ESDPt,
								"", "p_{t} [GeV/c]", "N_{clusters} Positron", "",
								kFALSE, 10., 140.,
								kTRUE, 0.1, 180.);
		pad_c11->cd(2)->Update();	

		pad_c11->cd(3)->SetRightMargin(0.17); 		
		pad_c11->cd(3);
			DrawAutoGammaHisto2D(	Resolution_E_nTRDtracklets_ESDPt,
								"", "p_{t} [GeV/c]", "N_{tracklets} Electron", "",
								kFALSE, 10., 140.,
								kTRUE, 0.1, 180.);
		pad_c11->cd(3)->Update();	

		pad_c11->cd(4)->SetRightMargin(0.17); 		
		pad_c11->cd(4);
			DrawAutoGammaHisto2D(	Resolution_P_nTRDtracklets_ESDPt,
								"", "p_{t} [GeV/c]", "N_{tracklets} Positron", "",
								kFALSE, 10., 140.,
								kTRUE, 0.1, 180.);
		pad_c11->cd(4)->Update();	

	
		c11->Update();




		ps_characteristics->Close();
}

//****************************************************************************************************************
//************************************ single plotting for better resolution *************************************
//****************************************************************************************************************


if(SinglePlots){

		leg3 = new TLegend( 0.7,0.75,0.92,0.9);
		leg3->SetTextSize(0.03);			
		leg3->SetFillColor(0);			
		leg3->AddEntry(Resolution_Gamma_dPt_Phi_2,("Photons"));		
		leg3->AddEntry(Resolution_E_dPt_Phi_2,("Electrons"));
		leg3->AddEntry(Resolution_P_dPt_Phi_2,("Positrons"));


		
    // ----------------------- Resolution dPt vs Phi ------------------------------------------------------------
		TCanvas * c2_2 = new TCanvas("c2_2","",10,10,700,1000);  // gives the page size		
		c2_2->cd();
		pad_c11 = new TPad("pad_c11","",0,0,1,1,0);   // gives the size of the histo areas 
		pad_c11->SetFillColor(0);
		pad_c11->GetFrame()->SetFillColor(0);
		pad_c11->SetBorderMode(0);
		pad_c11->Divide(1,2);
		pad_c11->Draw();
		pad_c11->cd(1);
			pad_c11->cd(1)->SetBottomMargin(0.005);
			StylingSliceHistos(Resolution_E_dPt_Phi_1,0.8);
			Resolution_E_dPt_Phi_1->SetMarkerColor(kBlue+2);
			Resolution_E_dPt_Phi_1->SetLineColor(kBlue-9)	;	
			DrawResolutionGammaHisto( Resolution_E_dPt_Phi_1, 
								"", "#phi [rad]", "Peak pos. dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kFALSE, -0.001, 0.001,
								kFALSE, 0.8, 180.);		
			StylingSliceHistos(Resolution_P_dPt_Phi_1,0.8);
			Resolution_P_dPt_Phi_1->SetMarkerColor(kRed+2);
			Resolution_P_dPt_Phi_1->SetLineColor(kRed-9)	;	
			Resolution_P_dPt_Phi_1->Draw("same");
			StylingSliceHistos(Resolution_Gamma_dPt_Phi_1,0.8);
			Resolution_Gamma_dPt_Phi_1->SetMarkerColor(kCyan+3);
			Resolution_Gamma_dPt_Phi_1->SetLineColor(kCyan-9)	;	
			Resolution_Gamma_dPt_Phi_1->Draw("same");		
			leg3->Draw("same");
			DrawGammaLines(-TMath::Pi(), TMath::Pi(), 0,0,0.005);
		pad_c11->cd(2)->Update();
		pad_c11->cd(2);
			pad_c11->cd(2)->SetTopMargin(0+0.01);
			StylingSliceHistos(Resolution_E_dPt_Phi_2,0.8);
			Resolution_E_dPt_Phi_2->SetMarkerColor(kBlue+2);
			Resolution_E_dPt_Phi_2->SetLineColor(kBlue-9)	;	
			DrawResolutionGammaHisto( Resolution_E_dPt_Phi_2, 
								"", "#phi [rad]", "#sigma dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, 0, 6.,
								kFALSE, 0.8, 180.);		
			StylingSliceHistos(Resolution_P_dPt_Phi_2,0.8);
			Resolution_P_dPt_Phi_2->SetMarkerColor(kRed+2);
			Resolution_P_dPt_Phi_2->SetLineColor(kRed-9)	;	
			Resolution_P_dPt_Phi_2->Draw("same");
			StylingSliceHistos(Resolution_Gamma_dPt_Phi_2,0.8);
			Resolution_Gamma_dPt_Phi_2->SetMarkerColor(kCyan+3);
			Resolution_Gamma_dPt_Phi_2->SetLineColor(kCyan-9)	;	
			Resolution_Gamma_dPt_Phi_2->Draw("same");		
//			leg3->Draw("same");
		c2_2->Update();	
		c2_2->SaveAs(Form("%s%s/Resolution_dPt_Phi_1.%s",path,suffix,suffix));
		delete pad_c11;	
		delete c2_2;

	// ---------------------------- Resolution Electron, Gamma and Positron dPt s Pt -----------------------------------------	
		TCanvas * c2_2 = new TCanvas("c2_2","",10,10,700,1000);  // gives the page size		
		c2_2->cd();
		pad_c11 = new TPad("pad_c11","",0,0,1,1,0);   // gives the size of the histo areas 
		pad_c11->SetFillColor(0);
		pad_c11->GetFrame()->SetFillColor(0);
		pad_c11->SetBorderMode(0);
		pad_c11->Divide(1,2);
		pad_c11->Draw();

		pad_c11->cd(1);
			pad_c11->cd(1)->SetBottomMargin(0.005);
			StylingSliceHistos(Resolution_E_dPt_Ptrebin_1,0.8);
			Resolution_E_dPt_Ptrebin_1->SetMarkerColor(kBlue+2);
			Resolution_E_dPt_Ptrebin_1->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_E_dPt_Ptrebin_1, 
							"", "p_{t} [GeV/c]", " Peak position dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, -1.9, 3.,
								kTRUE, 0.8, 10.);
			StylingSliceHistos(Resolution_P_dPt_Ptrebin_1,0.8);
			Resolution_P_dPt_Ptrebin_1->SetMarkerColor(kRed+2);
			Resolution_P_dPt_Ptrebin_1->SetLineColor(kRed-9)	;	
			Resolution_P_dPt_Ptrebin_1->Draw("same");
			StylingSliceHistos(Resolution_Gamma_dPt_Ptrebin_1,0.8);
			Resolution_Gamma_dPt_Ptrebin_1->SetMarkerColor(kCyan+3);
			Resolution_Gamma_dPt_Ptrebin_1->SetLineColor(kCyan-9)	;	
			Resolution_Gamma_dPt_Ptrebin_1->Draw("same");
			DrawGammaLines(0.8, 50, 0,0,0.005);		
			leg3->Draw("same");
		pad_c11->cd(2)->Update();
		pad_c11->cd(2);
			pad_c11->cd(2)->SetTopMargin(0+0.01);
			StylingSliceHistos(Resolution_E_dPt_Ptrebin_2,0.8);
			Resolution_E_dPt_Ptrebin_2->SetMarkerColor(kBlue+2);
			Resolution_E_dPt_Ptrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_E_dPt_Ptrebin_2, 
							"", "p_{t} [GeV/c]", " #sigma dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, 0, 21.,
								kTRUE, 0.8, 10.);
			StylingSliceHistos(Resolution_P_dPt_Ptrebin_2,0.8);
			Resolution_P_dPt_Ptrebin_2->SetMarkerColor(kRed+2);
			Resolution_P_dPt_Ptrebin_2->SetLineColor(kRed-9)	;	
			Resolution_P_dPt_Ptrebin_2->Draw("same");
			StylingSliceHistos(Resolution_Gamma_dPt_Ptrebin_2,0.8);
			Resolution_Gamma_dPt_Ptrebin_2->SetMarkerColor(kCyan+3);
			Resolution_Gamma_dPt_Ptrebin_2->SetLineColor(kCyan-9)	;	
			Resolution_Gamma_dPt_Ptrebin_2->Draw("same");		
		pad_c11->cd(2)->Update();
		c2_2->Update();	
		c2_2->SaveAs(Form("%s%s/Resolution_dPt_Pt_all.%s",path,suffix,suffix));
		delete pad_c11;	
		delete c2_2;

	// ----------------------------- Resolution dPhiAbs vs R ---------------------------------------
		TCanvas * c2_2 = new TCanvas("c2_2","",10,10,700,1000);  // gives the page size		
		c2_2->cd();
		pad_c11 = new TPad("pad_c11","",0,0,1,1,0);   // gives the size of the histo areas 
		pad_c11->SetFillColor(0);
		pad_c11->GetFrame()->SetFillColor(0);
		pad_c11->SetBorderMode(0);
		pad_c11->Divide(1,2);
		pad_c11->Draw();

		pad_c11->cd(1);
			pad_c11->cd(1)->SetBottomMargin(0.005);
			StylingSliceHistos(Resolution_dPhiAbs_VS_Rrebin_1,0.8);
			Resolution_dPhiAbs_VS_Rrebin_1->SetMarkerColor(kRed+2);
			Resolution_dPhiAbs_VS_Rrebin_1->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Rrebin_1, 
								"", "R [cm]", "Peak position d#phi [rad]", 
								kFALSE, 10., 140.,
								kTRUE, -0.00038, 0.0004,
								kTRUE, 0.8, 120.);
			DrawGammaLines(0.8, 120, 0,0,0.005);
		pad_c11->cd(1)->Update();
		pad_c11->cd(2);
			pad_c11->cd(2)->SetTopMargin(0+0.04);
			StylingSliceHistos(Resolution_dPhiAbs_VS_Rrebin_2,0.8);
			Resolution_dPhiAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
			Resolution_dPhiAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Rrebin_2, 
							"", "R [cm]", "#sigma d#phi [rad]", 
								kFALSE, 10., 140.,
								kFALSE, 0.0, 0.00375,
								kTRUE, 0.8, 120.);
		pad_c11->cd(2)->Update();
		c2_2->Update();	
		c2_2->SaveAs(Form("%s%s/Resolution_dphiAbs_R.%s",path,suffix,suffix));
		delete pad_c11;	
		delete c2_2;

	// ----------------------------- Resolution dPt vs R -------------------------------
		TCanvas * c2_2 = new TCanvas("c2_2","",10,10,700,1000);  // gives the page size		
		c2_2->cd();
		pad_c11 = new TPad("pad_c11","",0,0,1,1,0);   // gives the size of the histo areas 
		pad_c11->SetFillColor(0);
		pad_c11->GetFrame()->SetFillColor(0);
		pad_c11->SetBorderMode(0);
		pad_c11->Divide(1,2);
		pad_c11->Draw();

		pad_c11->cd(1);
			pad_c11->cd(1)->SetBottomMargin(0.005);
			StylingSliceHistos(Resolution_R_dPtrebin_1,0.8);
			Resolution_R_dPtrebin_1->SetMarkerColor(kRed+2);
			Resolution_R_dPtrebin_1->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_R_dPtrebin_1, 
								"", "R [cm]", "Peak position dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, -1.7, 1.7,
								kTRUE, 0.8, 120.);
			DrawGammaLines(0.8, 120, 0,0,0.005);
		pad_c11->cd(2)->Update();
		pad_c11->cd(2);
			pad_c11->cd(2)->SetTopMargin(0+0.005);
			StylingSliceHistos(Resolution_R_dPtrebin_2,0.8);
			Resolution_R_dPtrebin_2->SetMarkerColor(kRed+2);
			Resolution_R_dPtrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_R_dPtrebin_2, 
							"", "R [cm]", "#sigma dp_{t}/p_{t,MC} [%]", 
								kFALSE, 10., 140.,
								kTRUE, 0,4.15,
								kTRUE, 0.8, 120.);
		pad_c11->cd(2)->Update();
		c2_2->Update();	
		c2_2->SaveAs(Form("%s%s/Resolution_dPt_R.%s",path,suffix,suffix));
		delete pad_c11;	
		delete c2_2;
     
	// ---------------------------------- Resolution dRAbs vs R -------------------------------
		TCanvas * c2_2 = new TCanvas("c2_2","",10,10,700,1000);  // gives the page size		
		c2_2->cd();
		pad_c11 = new TPad("pad_c11","",0,0,1,1,0);   // gives the size of the histo areas 
		pad_c11->SetFillColor(0);
		pad_c11->GetFrame()->SetFillColor(0);
		pad_c11->SetBorderMode(0);
		pad_c11->Divide(1,2);
		pad_c11->Draw();

		pad_c11->cd(1);
			pad_c11->cd(1)->SetBottomMargin(0.005);
			StylingSliceHistos(Resolution_dRAbs_VS_Rrebin_1,0.8);
			Resolution_dRAbs_VS_Rrebin_1->SetMarkerColor(kRed+2);
			Resolution_dRAbs_VS_Rrebin_1->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_dRAbs_VS_Rrebin_1, 
								"", "R [cm]", "Peak pos. dR [cm]", 
								kFALSE, 10., 140.,
								kFALSE, -1, 1.5,
								kTRUE, 0.8, 120.);
			DrawGammaLines(0.8, 120, 0,0,0.005);
		pad_c11->cd(2)->Update();
		pad_c11->cd(2);
			pad_c11->cd(2)->SetTopMargin(0+0.01);
			StylingSliceHistos(Resolution_dRAbs_VS_Rrebin_2,0.8);
			Resolution_dRAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
			Resolution_dRAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_dRAbs_VS_Rrebin_2, 
							"", "R [cm]", "#sigma dR [cm]", 
								kFALSE, 10., 140.,
								kTRUE, 0., 4.75,
								kTRUE, 0.8, 120.);
		pad_c11->cd(2)->Update();
		c2_2->Update();	
		c2_2->SaveAs(Form("%s%s/Resolution_dRAbs_R.%s",path,suffix,suffix));
		delete pad_c11;	
		delete c2_2;

	// -------------------------------- Resolution dZAbs vs R -----------------------------	
		TCanvas * c2_2 = new TCanvas("c2_2","",10,10,700,1000);  // gives the page size		
		c2_2->cd();
		pad_c11 = new TPad("pad_c11","",0,0,1,1,0);   // gives the size of the histo areas 
		pad_c11->SetFillColor(0);
		pad_c11->GetFrame()->SetFillColor(0);
		pad_c11->SetBorderMode(0);
		pad_c11->Divide(1,2);
		pad_c11->Draw();

		pad_c11->cd(1);
			pad_c11->cd(1)->SetBottomMargin(0.005);
			StylingSliceHistos(Resolution_dZAbs_VS_Rrebin_1,0.8);
			Resolution_dZAbs_VS_Rrebin_1->SetMarkerColor(kRed+2);
			Resolution_dZAbs_VS_Rrebin_1->SetLineColor(kBlue-8)	;	
			DrawResolutionGammaHisto( Resolution_dZAbs_VS_Rrebin_1, 
								"", "R [cm]", "Peak pos. dZ [cm]", 
								kFALSE, 10., 140.,
								kTRUE, -0.225, 0.2,
								kTRUE, 0.8, 120.);
			DrawGammaLines(0.8, 120, 0,0,0.005);
		pad_c11->cd(2)->Update();
		pad_c11->cd(2);
			pad_c11->cd(2)->SetTopMargin(0+0.005);
			StylingSliceHistos(Resolution_dZAbs_VS_Rrebin_2,0.8);
			Resolution_dZAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
			Resolution_dZAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_dZAbs_VS_Rrebin_2, 
							"", "R [cm]", "#sigma dZ [cm]", 
								kFALSE, 10., 140.,
								kTRUE, 0., 3.75,
								kTRUE, 0.8, 120.);
		pad_c11->cd(2)->Update();
		c2_2->Update();	
		c2_2->SaveAs(Form("%s%s/Resolution_dZAbs_R.%s",path,suffix,suffix));
		delete pad_c11;	
		delete c2_2;

	// -------------------------------- Resolution dPt vs Pt for different ITS clusters -----------
		TCanvas * c2_2 = new TCanvas("c2_2","",10,10,700,1000);  // gives the page size		
		c2_2->cd();
		pad_c11 = new TPad("pad_c11","",0,0,1,1,0);   // gives the size of the histo areas 
		pad_c11->SetFillColor(0);
		pad_c11->GetFrame()->SetFillColor(0);
		pad_c11->SetBorderMode(0);
		pad_c11->Divide(1,2);
		pad_c11->Draw();

		pad_c11->cd(1);
			pad_c11->cd(1)->SetBottomMargin(0.005);
			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS0rebin_2,0.8);
			Resolution_E_dPt_Pt_ITS0rebin_2->SetMarkerColor(kRed+2);
			Resolution_E_dPt_Pt_ITS0rebin_2->SetLineColor(kRed+2);
			DrawResolutionGammaHisto( Resolution_E_dPt_Pt_ITS0rebin_2, 
								"", "p_{t} [GeV/c]", "#sigma dp_{t}/p_{t,MC} [%] for Electrons", 
								kFALSE, 10., 140.,
								kTRUE, -1, 30.,
								kTRUE, 0, 10.);
/*			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS1rebin_2,1);
			Resolution_E_dPt_Pt_ITS1rebin_2->SetMarkerColor(kBlue+2);
			Resolution_E_dPt_Pt_ITS1rebin_2->SetLineColor(kBlue+2);
			Resolution_E_dPt_Pt_ITS1rebin_2->Draw("same");
*/			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS2rebin_2,0.8);
			Resolution_E_dPt_Pt_ITS2rebin_2->SetMarkerColor(kViolet+2);
			Resolution_E_dPt_Pt_ITS2rebin_2->SetLineColor(kViolet+2);
			Resolution_E_dPt_Pt_ITS2rebin_2->Draw("same");
/*			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS3rebin_2,1);
			Resolution_E_dPt_Pt_ITS3rebin_2->SetMarkerColor(kOrange+2);
			Resolution_E_dPt_Pt_ITS3rebin_2->SetLineColor(kOrange+2);
			Resolution_E_dPt_Pt_ITS3rebin_2->Draw("same");
*/			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS4rebin_2,0.8);
			Resolution_E_dPt_Pt_ITS4rebin_2->SetMarkerColor(kCyan+2);
			Resolution_E_dPt_Pt_ITS4rebin_2->SetLineColor(kCyan+2);
			Resolution_E_dPt_Pt_ITS4rebin_2->Draw("same");
			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS5rebin_2,0.8);
			Resolution_E_dPt_Pt_ITS5rebin_2->SetMarkerColor(kGreen+2);
			Resolution_E_dPt_Pt_ITS5rebin_2->SetLineColor(kGreen+2);
			Resolution_E_dPt_Pt_ITS5rebin_2->Draw("same");
			StylingSliceHistos(	Resolution_E_dPt_Pt_ITS6rebin_2,0.8);
			Resolution_E_dPt_Pt_ITS6rebin_2->SetMarkerColor(kMagenta+2);
			Resolution_E_dPt_Pt_ITS6rebin_2->SetLineColor(kMagenta+2);
			Resolution_E_dPt_Pt_ITS6rebin_2->Draw("same");
			leg2 = new TLegend( 0.7,0.72,0.92,0.9);
			leg2->SetTextSize(0.03);			
			leg2->SetFillColor(0);			
			leg2->AddEntry(Resolution_E_dPt_Pt_ITS0rebin_2,("only TPC"));
	//		leg2->AddEntry(Resolution_E_dPt_Pt_ITS1rebin_2,("1 ITS cluster"));
			leg2->AddEntry(Resolution_E_dPt_Pt_ITS2rebin_2,("2 ITS cluster"));
	//		leg2->AddEntry(Resolution_E_dPt_Pt_ITS3rebin_2,("3 ITS cluster"));
			leg2->AddEntry(Resolution_E_dPt_Pt_ITS4rebin_2,("4 ITS cluster"));
			leg2->AddEntry(Resolution_E_dPt_Pt_ITS5rebin_2,("5 ITS cluster"));
			leg2->AddEntry(Resolution_E_dPt_Pt_ITS6rebin_2,("6 ITS cluster"));		
			leg2->Draw("same");
		pad_c11->cd(1)->Update();
		pad_c11->cd(2);
			pad_c11->cd(2)->SetTopMargin(0+0.01);
			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS0rebin_2,0.8);
			Resolution_P_dPt_Pt_ITS0rebin_2->SetMarkerColor(kRed+2);
			Resolution_P_dPt_Pt_ITS0rebin_2->SetLineColor(kRed+2);
			DrawResolutionGammaHisto( Resolution_P_dPt_Pt_ITS0rebin_2, 
								"", "p_{t} [GeV/c]", "#sigma dp_{t}/p_{t,MC} [%] for Positrons", 
								kFALSE, 10., 140.,
								kTRUE, -1, 30.,
								kTRUE, 0., 10.);
/*			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS1rebin_2,1);
			Resolution_P_dPt_Pt_ITS1rebin_2->SetMarkerColor(kBlue+2);
			Resolution_P_dPt_Pt_ITS1rebin_2->SetLineColor(kBlue+2);
			Resolution_P_dPt_Pt_ITS1rebin_2->Draw("same");
*/			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS2rebin_2,0.8);
			Resolution_P_dPt_Pt_ITS2rebin_2->SetMarkerColor(kViolet+2);
			Resolution_P_dPt_Pt_ITS2rebin_2->SetLineColor(kViolet+2);
			Resolution_P_dPt_Pt_ITS2rebin_2->Draw("same");
/*			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS3rebin_2,1);
			Resolution_P_dPt_Pt_ITS3rebin_2->SetMarkerColor(kOrange+2);
			Resolution_P_dPt_Pt_ITS3rebin_2->SetLineColor(kOrange+2);
			Resolution_P_dPt_Pt_ITS3rebin_2->Draw("same");
*/			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS4rebin_2,0.8);
			Resolution_P_dPt_Pt_ITS4rebin_2->SetMarkerColor(kCyan+2);
			Resolution_P_dPt_Pt_ITS4rebin_2->SetLineColor(kCyan+2);
			Resolution_P_dPt_Pt_ITS4rebin_2->Draw("same");
			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS5rebin_2,1);
			Resolution_P_dPt_Pt_ITS5rebin_2->SetMarkerColor(kGreen+2);
			Resolution_P_dPt_Pt_ITS5rebin_2->SetLineColor(kGreen+2);
			Resolution_P_dPt_Pt_ITS5rebin_2->Draw("same");
			StylingSliceHistos(	Resolution_P_dPt_Pt_ITS6rebin_2,0.8);
			Resolution_P_dPt_Pt_ITS6rebin_2->SetMarkerColor(kMagenta+2);
			Resolution_P_dPt_Pt_ITS6rebin_2->SetLineColor(kMagenta+2);
			Resolution_P_dPt_Pt_ITS6rebin_2->Draw("same");
		pad_c11->cd(2)->Update();
		c2_2->Update();	
		c2_2->SaveAs(Form("%s%s/Resolution_dPt_PtITScluster.%s",path,suffix,suffix));
		delete pad_c11;	
		delete c2_2;

	// --------------------  TRD cluster and tracklet distributions for E and P of conversion --------------
		TCanvas * c2_2 = new TCanvas("c2_2","",10,10,2400,2000);  // gives the page size		
		c2_2->cd();
		pad_c11 = new TPad("pad_c11","",0,0,1,1,0);   // gives the size of the histo areas 
		pad_c11->SetFillColor(0);
		pad_c11->GetFrame()->SetFillColor(0);
		pad_c11->SetBorderMode(0);
		pad_c11->Divide(2,2);
		pad_c11->Draw();


		pad_c11->cd(1)->SetLogz(1);
		pad_c11->cd(1)->SetLogx(1);
		pad_c11->cd(2)->SetLogz(1);
		pad_c11->cd(2)->SetLogx(1);
		pad_c11->cd(3)->SetLogz(1);
		pad_c11->cd(3)->SetLogx(1);
		pad_c11->cd(4)->SetLogz(1);
		pad_c11->cd(4)->SetLogx(1);


		pad_c11->cd(1)->SetRightMargin(0.14); 		
		pad_c11->cd(1)->SetLeftMargin(0.14);
		pad_c11->cd(1)->SetTopMargin(0.01);
		pad_c11->cd(1)->SetBottomMargin(0.009); 		
		pad_c11->cd(1);
			DrawAutoGammaHisto2DRes(	Resolution_E_nTRDclusters_ESDPt,
								"", "p_{t} [GeV/c]", "N_{clusters} Electron", "",
								kFALSE, 10., 140.,
								kTRUE, 0.1, 180.);
		pad_c11->cd(1)->Update();	

		pad_c11->cd(2)->SetRightMargin(0.14);
		pad_c11->cd(2)->SetLeftMargin(0.14); 		
		pad_c11->cd(2)->SetTopMargin(0.01); 	
		pad_c11->cd(2)->SetBottomMargin(0.009);	
		pad_c11->cd(2);
			DrawAutoGammaHisto2DRes(	Resolution_P_nTRDclusters_ESDPt,
								"", "p_{t} [GeV/c]", "N_{clusters} Positron", "",
								kFALSE, 10., 140.,
								kTRUE, 0.1, 180.);

		pad_c11->cd(2)->Update();	

		pad_c11->cd(3)->SetRightMargin(0.14); 		
		pad_c11->cd(3)->SetLeftMargin(0.14);
		pad_c11->cd(3)->SetTopMargin(0.006); 		
		pad_c11->cd(3);
			DrawAutoGammaHisto2DRes(	Resolution_E_nTRDtracklets_ESDPt,
								"", "p_{t} [GeV/c]", "N_{tracklets} Electron", "",
								kFALSE, 10., 140.,
								kTRUE, 0.1, 180.);
		pad_c11->cd(3)->Update();	

		pad_c11->cd(4)->SetRightMargin(0.14); 		
		pad_c11->cd(4)->SetLeftMargin(0.14);
		pad_c11->cd(4)->SetTopMargin(0.006); 	
		pad_c11->cd(4);
			DrawAutoGammaHisto2DRes(	Resolution_P_nTRDtracklets_ESDPt,
								"", "p_{t} [GeV/c]", "N_{tracklets} Positron", "",
								kFALSE, 10., 140.,
								kTRUE, 0.1, 180.);
		pad_c11->cd(4)->Update();	
		c2_2->Update();	
		c2_2->SaveAs(Form("%s%s/TRDstudies.%s",path,suffix,suffix));
		delete pad_c11;	
		delete c2_2;

		leg4 = new TLegend( 0.7,0.75,0.92,0.9);
		leg4->SetTextSize(0.03);			
		leg4->SetFillColor(0);			
		leg4->AddEntry(ESD_ConvGamma_Pt,("Photons"));		
		leg4->AddEntry(ESD_E_Pt,("Electrons"));
		leg4->AddEntry(ESD_P_Pt,("Positrons"));

    // ------------------ Only \sigma for already created plots ------------------------------
		TCanvas * c2_2 = new TCanvas("c2_2","",10,10,700,500);  // gives the page size		
		c2_2->SetLogy(1);
		c2_2->cd();
			ESD_E_Pt->SetLineColor(kBlue+2)	;	
			DrawResolutionGammaHisto( ESD_E_Pt, 
								"", "p_{t} [GeV/c]", "N_{#gamma}", 
								kFALSE, 10., 140.,
								kFALSE, -0.001, 0.001,
								kFALSE, 0.8, 180.);		
			ESD_P_Pt->SetLineColor(kRed+2)	;	
			ESD_P_Pt->Draw("same");
			ESD_ConvGamma_Pt->SetLineColor(kCyan+1)	;	
			ESD_ConvGamma_Pt->Draw("same");		
			leg4->Draw("same");
		c2_2->Update();	
		c2_2->SaveAs(Form("%s%s/Resolution_Ptdistribution.%s",path,suffix,suffix));
		delete c2_2;

		TCanvas * c2_2 = new TCanvas("c2_2","",10,10,700,500);  // gives the page size		
		c2_2->cd();
			c2_2->SetTopMargin(0+0.005);
			StylingSliceHistos(Resolution_dZAbs_VS_Rrebin_2,0.8);
			Resolution_dZAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
			Resolution_dZAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_dZAbs_VS_Rrebin_2, 
							"", "R [cm]", "#sigma dZ [cm]", 
								kFALSE, 10., 140.,
								kTRUE, 0., 2.1,
								kTRUE, 0.8, 120.);
		c2_2->Update();	
		c2_2->SaveAs(Form("%s%s/Resolution_Zdistributionsingle.%s",path,suffix,suffix));
		delete c2_2;

		TCanvas * c2_2 = new TCanvas("c2_2","",10,10,700,500);  // gives the page size		
		c2_2->cd();
			c2_2->SetTopMargin(0+0.04);
			StylingSliceHistos(Resolution_dPhiAbs_VS_Rrebin_2,0.8);
			Resolution_dPhiAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
			Resolution_dPhiAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Rrebin_2, 
							"", "R [cm]", "#sigma d#phi [rad]", 
								kFALSE, 10., 140.,
								kTRUE, 0.0, 0.01,
								kTRUE, 0.8, 120.);
		c2_2->Update();	
		c2_2->SaveAs(Form("%s%s/Resolution_Phidistributionsingle.%s",path,suffix,suffix));
		delete c2_2;

		TCanvas * c2_2 = new TCanvas("c2_2","",10,10,700,500);  // gives the page size		
		c2_2->cd();
			c2_2->SetTopMargin(0+0.01);
			StylingSliceHistos(Resolution_dRAbs_VS_Rrebin_2,0.8);
			Resolution_dRAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
			Resolution_dRAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
			DrawResolutionGammaHisto( Resolution_dRAbs_VS_Rrebin_2, 
							"", "R [cm]", "#sigma dR [cm]", 
								kFALSE, 10., 140.,
								kTRUE, 0., 3.1,
								kTRUE, 0.8, 120.);
		c2_2->Update();	
		c2_2->SaveAs(Form("%s%s/Resolution_Rdistributionsingle.%s",path,suffix,suffix));
		delete c2_2;


	}
	
}


