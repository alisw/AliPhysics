/****************************************************************************************************************************
****** 	provided by Gamma Conversion Group, PWGGA, 														*****
******		Friederike Bock, friederike.bock@cern.ch													*****
*****************************************************************************************************************************
*** This macro can be used to observe the effects of the different cuts used in the GammaConversion analysis, it can 	*****
*** operate on root files created by the GammaConversionTask. It can take two root-files, where the second one should be ****
*** a MC file otherwise all histograms including MC need to be commented out.									*****
****************************************************************************************************************************/

#include <Riostream.h>
#include <fstream>
#include "PlottingGammaConversionHistos.h"
#include "PlottingGammaConversionAdditional.h"
using namespace std;


extern TRandom *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem *gSystem;


void  Cuts_Events_new(const char *data = "myOutput", const char *MCfile = "",const char *cutsel = "", const char *path = "", const char *output = "Reconstruction_Cuts", const char *Plots = "kTRUE", const char *suffix = "gif", const char *outlname = "Export_Cuts2.dat"){	
	
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
//	StyleSettingsThesis();	
	StyleSettings();	
	set_plot_style();	
	
	Bool_t SinglePlots = kFALSE;
	if(Plots == "kTRUE") SinglePlots = kTRUE;
	
	// textsize for legend
	Float_t ts=0.04;
	// textsize for label
	Float_t ls = 0.04;
	
	//Array defintion for printing Logo in right upper corner
	Float_t right_up[4]={0.7,0.63,0.15, 0.02};
	Float_t right_down[4]={0.7,0.23,0.15, 0.02};
	//Array defintion for printing Logo in left upper corner
	Float_t left_up[4]={0.13,0.73, 0.15, 0.02};
	//Array defintion for printing text in right upper corner
	Float_t right_up_text[4]={0.6,0.8,0.15, 0.02};

	/* If Full is set to kTRUE you will display all of the Cuts compared to the Data before any cut, 
	if it is set to kFAlSE you will display them compared to the Data before the Cut you made. */
	Bool_t Full = kFALSE;
		
	/* If MC is set to kTRUE the second input will be treated as Montecarlo and the true Photons will
	 be extracted as well and if not it stays as it was.*/
	Bool_t MC1 = kFALSE;
	Bool_t MC2 = kTRUE;
	
	char *StandardYAxis = "#gamma/ event scaled by multiplicity";

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

	char *Date = "25th June 2010";
// --------------------------------- end of self definitions --------------------------------
	
	
	TPostScript *ps_characteristics;
	if(!SinglePlots)ps_characteristics = new TPostScript(Form("%s%s.ps",path,output),111);	
	
	Char_t filename_data[200] = (Form("%s%s",path,data));	
	TFile f(filename_data);  

	 TDirectory *fPWGGAGammaConversion_data = new TDirectory(); // definition of first folder / list	
	 TList *fHistosGammaConversion_data = new TList(); // definition of first folder / list
	 TList *fESDContainer_data = new TList();  // definition of following folder / list
	 TList *fMappingContainer_data = new TList();	 

	 if(!(fPWGGAGammaConversion_data = (TDirectory*)f.Get(GammaDirectory))) cout <<"PWGGAGammConversion TList NOT loaded correctly"<<endl; 
	 if(!(fHistosGammaConversion_data = (TList*)fPWGGAGammaConversion_data->Get(GammaList))) cout<<"histogramsAliGammaConversion NOT loaded correctly!"<<endl; 
	 if(!(fESDContainer_data = (TList*)fHistosGammaConversion_data->FindObject("ESD histograms"))) cout<<"ESD histograms NOT loaded correctly!"<<endl; 
 	fMappingContainer_data = (TList*)fHistosGammaConversion_data->FindObject("Mapping histograms"); 
	if(MC1){
			fTRUEContainer_data = (TList*)fHistosGammaConversion_data->FindObject("True conversion histograms");  
			TH1F *TRUE_Conversion_R_data = fTRUEContainer_data->FindObject("ESD_TrueConversion_R");
			TH1F *ESD_Conversion_R_data= fESDContainer_data->FindObject("ESD_Conversion_R");        
		}
	
	TH1F * ESD_CutLikeSign_data = fESDContainer_data->FindObject("ESD_CutLikeSign_InvMass");
	TH1F * ESD_CutRefit_data = fESDContainer_data->FindObject("ESD_CutRefit_InvMass");
	TH1F * ESD_CutKink_data = fESDContainer_data->FindObject("ESD_CutKink_InvMass");
	TH1F * ESD_CutdEdxSigmaElectronLine_data = fESDContainer_data->FindObject("ESD_CutdEdxSigmaElectronLine_InvMass");
	TH1F * ESD_CutdEdxSigmaPionLine_data = fESDContainer_data->FindObject("ESD_CutdEdxSigmaPionLine_InvMass");
	TH1F * ESD_CutKaonRejectionLowP_data = fESDContainer_data->FindObject("ESD_CutKaonRejectionLowP_InvMass");
	TH1F * ESD_CutProtonRejectionLowP_data = fESDContainer_data->FindObject("ESD_CutProtonRejectionLowP_InvMass");
	TH1F * ESD_CutPionRejectionLowP_data = fESDContainer_data->FindObject("ESD_CutPionRejectionLowP_InvMass");
	TH1F * ESD_CutGetOnFly_data = fESDContainer_data->FindObject("ESD_CutGetOnFly_InvMass");	
	TH1F * ESD_CutNContributors_data = fESDContainer_data->FindObject("ESD_CutNContributors_InvMass");
	TH1F * ESD_CutPIDProb_data = fESDContainer_data->FindObject("ESD_CutPIDProb_InvMass");
	TH1F * ESD_CutR_data = fESDContainer_data->FindObject("ESD_CutR_InvMass");	
	TH1F * ESD_CutLine_data = fESDContainer_data->FindObject("ESD_CutLine_InvMass");
	TH1F * ESD_CutZ_data = fESDContainer_data->FindObject("ESD_CutZ_InvMass");
	TH1F * ESD_CutMinNClsTPC_data = fESDContainer_data->FindObject("ESD_CutMinNClsTPC_InvMass");
	TH1F * ESD_CutSinglePt_data = fESDContainer_data->FindObject("ESD_CutSinglePt_InvMass");
	TH1F * ESD_CutNDF_data = fESDContainer_data->FindObject("ESD_CutNDF_InvMass");	
	TH1F * ESD_CutChi2_data = fESDContainer_data->FindObject("ESD_CutChi2_InvMass");
	TH1F * ESD_CutEta_data = fESDContainer_data->FindObject("ESD_CutEta_InvMass");
	TH1F * ESD_CutPt_data = fESDContainer_data->FindObject("ESD_CutPt_InvMass");	
	TH1F * ESD_ConvGamma_Pt_data=fESDContainer_data->FindObject("ESD_ConvGamma_Pt");
	TH1F * ESD_NumberOfContributorsVtx_data=fESDContainer_data->FindObject("ESD_NumberOfContributorsVtx");
	TH1F * ESD_ConvGamma_Mass_data = fESDContainer_data->FindObject("ESD_ConvGamma_Mass");
	TH1F * ESD_NumberOfV0_data = fESDContainer_data->FindObject("ESD_NumberOfV0s");
	TH1F * ESD_AllV0_data = fESDContainer_data->FindObject("ESD_AllV0s_InvMass");
	TH1F * ESD_NumberOfGoodESDTracks_data=fESDContainer_data->FindObject("ESD_NumberOfGoodESDTracks");
	TH1F *ScalingDiagramm_data= fMappingContainer_data->FindObject("ESD_Conversion_Mapping_Phi_in_R_11");

	//get Entries
	ESD_NumberOfContributorsVtx_data->SetAxisRange(1.,100.);
	ESD_NumberOfGoodESDTracks_data->SetAxisRange(1.,100.);

	
	Float_t Scaling_data = ScalingDiagramm_data->Integral();
	Float_t nGoodEvents_data = ESD_NumberOfContributorsVtx_data->Integral();
	Float_t nGoodTrig_data = ESD_NumberOfContributorsVtx_data->GetEntries();
	Float_t nRecGamma_data = ESD_ConvGamma_Pt_data->Integral();
	cout<< data << "    Number of events::   " << nGoodEvents_data << "    Number of triggers::   " << nGoodTrig_data <<endl;

	Float_t normFac_data=1./nGoodEvents_data;
//	Float_t normFac_data=1./Scaling_data;
	Double_t mean_data = ESD_NumberOfGoodESDTracks_data->GetMean();

	//Scaling reconstr.
	ESD_CutLikeSign_data->Sumw2();
	ESD_CutRefit_data->Sumw2();
	ESD_CutKink_data->Sumw2();
	ESD_CutdEdxSigmaElectronLine_data->Sumw2();
	ESD_CutdEdxSigmaPionLine_data->Sumw2();
	ESD_CutKaonRejectionLowP_data->Sumw2();
	ESD_CutProtonRejectionLowP_data->Sumw2();
	ESD_CutPionRejectionLowP_data->Sumw2();
	ESD_CutGetOnFly_data->Sumw2();
	ESD_CutNContributors_data->Sumw2();	
	ESD_CutPIDProb_data->Sumw2();
	ESD_CutR_data->Sumw2();
	ESD_CutLine_data->Sumw2();
	ESD_CutZ_data->Sumw2();
	ESD_CutMinNClsTPC_data->Sumw2();
	ESD_CutSinglePt_data->Sumw2();
	ESD_CutNDF_data->Sumw2();
	ESD_CutChi2_data->Sumw2(); 
	ESD_CutEta_data->Sumw2();
	ESD_CutPt_data->Sumw2();
	ESD_ConvGamma_Mass_data->Sumw2();
	ESD_AllV0_data->Sumw2();
	// Rebinning, has to be done before merging otherwise it would effect everytime we add the file to another

	const int rebin = 20 ;
	
	ESD_CutLikeSign_data->Rebin(rebin);
	ESD_CutRefit_data->Rebin(rebin);
	ESD_CutKink_data->Rebin(rebin);
	ESD_CutdEdxSigmaElectronLine_data->Rebin(rebin);
	ESD_CutdEdxSigmaPionLine_data->Rebin(rebin);
	ESD_CutKaonRejectionLowP_data->Rebin(rebin);
	ESD_CutProtonRejectionLowP_data->Rebin(rebin);
	ESD_CutPionRejectionLowP_data->Rebin(rebin);
	ESD_CutGetOnFly_data->Rebin(rebin);
	ESD_CutNContributors_data->Rebin(rebin);
	ESD_CutPIDProb_data->Rebin(rebin);
	ESD_CutR_data->Rebin(rebin);
	ESD_CutLine_data->Rebin(rebin);
	ESD_CutZ_data->Rebin(rebin);
	ESD_CutMinNClsTPC_data->Rebin(rebin);
	ESD_CutSinglePt_data->Rebin(rebin);
	ESD_CutNDF_data->Rebin(rebin);
	ESD_CutChi2_data->Rebin(rebin);
	ESD_CutEta_data->Rebin(rebin);
	ESD_CutPt_data->Rebin(rebin);
	ESD_ConvGamma_Mass_data->Rebin(rebin);
	ESD_AllV0_data->Rebin(rebin);

	// Creation of histograms to compare cut with data before	

	TH1F* NoCutGetOnFly_data = (TH1F*) ESD_AllV0_data->Clone();

	TH1F* NoCutLikeSign_data = (TH1F*) NoCutGetOnFly_data->Clone();
	NoCutLikeSign_data->Add(ESD_CutGetOnFly_data,-1);

	TH1F* NoCutRefit_data = (TH1F*) NoCutLikeSign_data->Clone();		
	NoCutRefit_data->Add(ESD_CutLikeSign_data,-1);

	TH1F* NoCutKink_data = (TH1F*) NoCutRefit_data->Clone();		
	NoCutKink_data->Add(ESD_CutRefit_data,-1);

	TH1F* NoCutdEdxSigmaElectronLine_data = (TH1F*) NoCutKink_data->Clone();		
	NoCutdEdxSigmaElectronLine_data->Add(ESD_CutKink_data,-1);

	TH1F* NoCutdEdxSigmaPionLine_data = (TH1F*) NoCutdEdxSigmaElectronLine_data->Clone();		
	NoCutdEdxSigmaPionLine_data->Add(ESD_CutdEdxSigmaElectronLine_data,-1);

	TH1F* NoCutKaonRejectionLowP_data = (TH1F*) NoCutdEdxSigmaPionLine_data->Clone();		
	NoCutKaonRejectionLowP_data->Add(ESD_CutdEdxSigmaPionLine_data,-1);

	TH1F* NoCutProtonRejectionLowP_data = (TH1F*) NoCutKaonRejectionLowP_data->Clone();		
	NoCutProtonRejectionLowP_data->Add(ESD_CutKaonRejectionLowP_data,-1);

	TH1F* NoCutPionRejectionLowP_data = (TH1F*) NoCutProtonRejectionLowP_data->Clone();		
	NoCutPionRejectionLowP_data->Add(ESD_CutProtonRejectionLowP_data,-1);

	TH1F* NoCutNContributors_data = (TH1F*) NoCutPionRejectionLowP_data->Clone();		
	NoCutNContributors_data->Add(ESD_CutPionRejectionLowP_data,-1);

	TH1F* NoCutPIDProb_data = (TH1F*) NoCutNContributors_data->Clone();		
	NoCutPIDProb_data->Add(ESD_CutNContributors_data,-1);

	TH1F* NoCutR_data = (TH1F*) NoCutPIDProb_data->Clone();		
	NoCutR_data->Add(ESD_CutPIDProb_data,-1);

	TH1F* NoCutLine_data = (TH1F*) NoCutR_data->Clone();		
	NoCutLine_data->Add(ESD_CutR_data,-1);

	TH1F* NoCutZ_data = (TH1F*) NoCutLine_data->Clone();		
	NoCutZ_data->Add(ESD_CutLine_data,-1);

	TH1F* NoCutMinNClsTPC_data = (TH1F*) NoCutZ_data->Clone();		
	NoCutMinNClsTPC_data->Add(ESD_CutZ_data,-1);

	TH1F* NoCutSinglePt_data = (TH1F*) NoCutMinNClsTPC_data->Clone();		
	NoCutSinglePt_data->Add(ESD_CutMinNClsTPC_data,-1);

	TH1F* NoCutNDF_data = (TH1F*) NoCutSinglePt_data->Clone();		
	NoCutNDF_data->Add(ESD_CutSinglePt_data,-1);

	TH1F* NoCutChi2_data = (TH1F*) NoCutNDF_data->Clone();		
	NoCutChi2_data->Add(ESD_CutNDF_data,-1);

	TH1F* NoCutEta_data = (TH1F*) NoCutChi2_data->Clone();		
	NoCutEta_data->Add(ESD_CutChi2_data,-1);

	TH1F* NoCutPt_data = (TH1F*) NoCutEta_data->Clone();		
	NoCutPt_data->Add(ESD_CutEta_data,-1);

	TH1F* Final_data = (TH1F*) NoCutPt_data->Clone();			
	Final_data->Add(ESD_CutPt_data,-1);
	

	Double_t NumberConvGamma_data = ESD_ConvGamma_Mass_data->Integral();
	Double_t NumberV0_data= ESD_AllV0_data->Integral();
	//Scaling
	ESD_CutLikeSign_data->Scale(normFac_data);
	ESD_CutRefit_data->Scale(normFac_data);
	ESD_CutKink_data->Scale(normFac_data);
	ESD_CutdEdxSigmaElectronLine_data->Scale(normFac_data);
	ESD_CutdEdxSigmaPionLine_data->Scale(normFac_data);
	ESD_CutKaonRejectionLowP_data->Scale(normFac_data);
	ESD_CutProtonRejectionLowP_data->Scale(normFac_data);
	ESD_CutPionRejectionLowP_data->Scale(normFac_data);
	ESD_CutGetOnFly_data->Scale(normFac_data);
	ESD_CutNContributors_data->Scale(normFac_data);	
	ESD_CutPIDProb_data->Scale(normFac_data);
	ESD_CutR_data->Scale(normFac_data);
	ESD_CutLine_data->Scale(normFac_data);
	ESD_CutZ_data->Scale(normFac_data);
	ESD_CutMinNClsTPC_data->Scale(normFac_data);
	ESD_CutSinglePt_data->Scale(normFac_data);
	ESD_CutNDF_data->Scale(normFac_data);
	ESD_CutChi2_data->Scale(normFac_data); 
	ESD_CutEta_data->Scale(normFac_data);
	ESD_CutPt_data->Scale(normFac_data);
	ESD_ConvGamma_Mass_data->Scale(normFac_data);

	NoCutPt_data->Scale(normFac_data);
	NoCutEta_data->Scale(normFac_data);
	NoCutChi2_data->Scale(normFac_data);
	NoCutNDF_data->Scale(normFac_data);
	NoCutZ_data->Scale(normFac_data);	
	NoCutMinNClsTPC_data->Scale(normFac_data);
	NoCutSinglePt_data->Scale(normFac_data);
	NoCutLine_data->Scale(normFac_data);
	NoCutR_data->Scale(normFac_data);
	NoCutPIDProb_data->Scale(normFac_data);
	NoCutNContributors_data->Scale(normFac_data);
	NoCutGetOnFly_data->Scale(normFac_data);
	NoCutKaonRejectionLowP_data->Scale(normFac_data);
	NoCutProtonRejectionLowP_data->Scale(normFac_data);
	NoCutPionRejectionLowP_data->Scale(normFac_data);
	NoCutdEdxSigmaPionLine_data->Scale(normFac_data);
	NoCutdEdxSigmaElectronLine_data->Scale(normFac_data);
	NoCutKink_data->Scale(normFac_data);
	NoCutRefit_data->Scale(normFac_data);
	NoCutLikeSign_data->Scale(normFac_data);
	ESD_AllV0_data->Scale(normFac_data);
	Final_data->Scale(normFac_data);

	Double_t RealV0_data = 0. ;
	for(Int_t i=2;i<ESD_NumberOfV0_data->GetNbinsX()+1;i++){
	RealV0_data += ((ESD_NumberOfV0_data->GetBinCenter(i))* ESD_NumberOfV0_data->GetBinContent(i));
	}

	cout << "Number of V0's data:   " << NumberV0_data <<"   "  << RealV0_data << endl;

// Second input ------------- same as done for first input -------------------------------------------------------------------------------------
	TFile *montecarlo = 0x0;
	
	
	if(MCfile != ""){
		montecarlo = new TFile(Form("%s%s",path, MCfile));
		
		
		// for new versions
		 TDirectory *fPWGGAGammaConversion_montecarlo = new TDirectory(); // definition of first folder / list	
		 TList *fHistosGammaConversion_montecarlo = new TList(); // definition of first folder / list
		 TList *fESDContainer_montecarlo = new TList();  // definition of following folder / list
		 TList *fMappingContainer_montecarlo = new TList();  // definition of following folder / list

		 if(!(fPWGGAGammaConversion_montecarlo = (TDirectory*)montecarlo->Get(GammaDirectory))) cout <<"PWGGAGammConversion TList NOT loaded correctly"<<endl; 
		 if(!(fHistosGammaConversion_montecarlo = (TList*)fPWGGAGammaConversion_montecarlo->Get(GammaList))) cout<<"histogramsAliGammaConversion NOT loaded correctly!"<<endl; 
		 if(!(fESDContainer_montecarlo = (TList*)fHistosGammaConversion_montecarlo->FindObject("ESD histograms"))) cout<<"ESD histograms NOT loaded correctly!"<<endl; 
		fMappingContainer_montecarlo = (TList*)fHistosGammaConversion_montecarlo->FindObject("Mapping histograms");

		// If second Input is a Montecarlodata file
		if(MC2){
			fTRUEContainer_montecarlo = (TList*)fHistosGammaConversion_montecarlo->FindObject("True conversion histograms");  
			TH1F *TRUE_Conversion_R_montecarlo = fTRUEContainer_montecarlo->FindObject("ESD_TrueConversion_R");
			TH1F *ESD_Conversion_R_montecarlo= fESDContainer_montecarlo->FindObject("ESD_Conversion_R");        
		}
	
		TH1F * ESD_CutLikeSign_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutLikeSign_InvMass");
		TH1F * ESD_CutRefit_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutRefit_InvMass");
		TH1F * ESD_CutKink_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutKink_InvMass");
		TH1F * ESD_CutdEdxSigmaElectronLine_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutdEdxSigmaElectronLine_InvMass");
		TH1F * ESD_CutdEdxSigmaPionLine_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutdEdxSigmaPionLine_InvMass");
	TH1F * ESD_CutKaonRejectionLowP_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutKaonRejectionLowP_InvMass");
	TH1F * ESD_CutProtonRejectionLowP_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutProtonRejectionLowP_InvMass");
	TH1F * ESD_CutPionRejectionLowP_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutPionRejectionLowP_InvMass");

		TH1F * ESD_CutGetOnFly_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutGetOnFly_InvMass");	
		TH1F * ESD_CutNContributors_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutNContributors_InvMass");
		TH1F * ESD_CutPIDProb_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutPIDProb_InvMass");
		TH1F * ESD_CutR_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutR_InvMass");	
		TH1F * ESD_CutLine_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutLine_InvMass");
		TH1F * ESD_CutZ_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutZ_InvMass");
		TH1F * ESD_CutMinNClsTPC_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutMinNClsTPC_InvMass");
		TH1F * ESD_CutSinglePt_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutSinglePt_InvMass");
		TH1F * ESD_CutNDF_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutNDF_InvMass");	
		TH1F * ESD_CutChi2_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutChi2_InvMass");
		TH1F * ESD_CutEta_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutEta_InvMass");
		TH1F * ESD_CutPt_montecarlo = fESDContainer_montecarlo->FindObject("ESD_CutPt_InvMass");	
		TH1F * ESD_ConvGamma_Pt_montecarlo=fESDContainer_montecarlo->FindObject("ESD_ConvGamma_Pt");
		TH1F * ESD_NumberOfContributorsVtx_montecarlo=fESDContainer_montecarlo->FindObject("ESD_NumberOfContributorsVtx");
		TH1F * ESD_ConvGamma_Mass_montecarlo = fESDContainer_montecarlo->FindObject("ESD_ConvGamma_Mass");
		TH1F * ESD_NumberOfV0_montecarlo = fESDContainer_montecarlo->FindObject("ESD_NumberOfV0s");
		TH1F * ESD_AllV0_montecarlo = fESDContainer_montecarlo->FindObject("ESD_AllV0s_InvMass");
		TH1F * ESD_NumberOfGoodESDTracks_montecarlo=fESDContainer_montecarlo->FindObject("ESD_NumberOfGoodESDTracks");		
		TH1F *ScalingDiagramm_montecarlo= fMappingContainer_montecarlo->FindObject("ESD_Conversion_Mapping_Phi_in_R_11");
		//get Entries
		ESD_NumberOfContributorsVtx_montecarlo->SetAxisRange(1.,100.);
		ESD_NumberOfGoodESDTracks_montecarlo->SetAxisRange(1.,100.);	

		Float_t Scaling_montecarlo = ScalingDiagramm_montecarlo->Integral();
		Float_t nGoodEvents_montecarlo = ESD_NumberOfContributorsVtx_montecarlo->Integral();
		Float_t nGoodTrig_montecarlo = ESD_NumberOfContributorsVtx_montecarlo->GetEntries();
		Float_t nRecGamma_montecarlo = ESD_ConvGamma_Pt_montecarlo->Integral();
		cout<< MCfile << "    Number of events::   " << nGoodEvents_montecarlo << "    Number of triggers::   " << nGoodTrig_montecarlo << endl;
		Double_t mean_montecarlo = ESD_NumberOfGoodESDTracks_montecarlo->GetMean();
		Float_t normFac_montecarlo=1./nGoodEvents_montecarlo * mean_data/mean_montecarlo;
	//	Float_t normFac_montecarlo=1./Scaling_montecarlo;
		
		//Scaling reconstr.
		ESD_CutLikeSign_montecarlo->Sumw2();
		ESD_CutRefit_montecarlo->Sumw2();
		ESD_CutKink_montecarlo->Sumw2();
		ESD_CutdEdxSigmaElectronLine_montecarlo->Sumw2();
		ESD_CutdEdxSigmaPionLine_montecarlo->Sumw2();
		ESD_CutKaonRejectionLowP_montecarlo->Sumw2();
		ESD_CutProtonRejectionLowP_montecarlo->Sumw2();
		ESD_CutPionRejectionLowP_montecarlo->Sumw2();
		ESD_CutGetOnFly_montecarlo->Sumw2();
		ESD_CutNContributors_montecarlo->Sumw2();	
		ESD_CutPIDProb_montecarlo->Sumw2();
		ESD_CutR_montecarlo->Sumw2();
		ESD_CutLine_montecarlo->Sumw2();
		ESD_CutZ_montecarlo->Sumw2();
		ESD_CutMinNClsTPC_montecarlo->Sumw2();
		ESD_CutSinglePt_montecarlo->Sumw2();
		ESD_CutNDF_montecarlo->Sumw2();
		ESD_CutChi2_montecarlo->Sumw2(); 
		ESD_CutEta_montecarlo->Sumw2();
		ESD_CutPt_montecarlo->Sumw2();
		ESD_ConvGamma_Mass_montecarlo->Sumw2();
		ESD_AllV0_montecarlo->Sumw2();
		// Rebinning, has to be done before merging otherwise it would effect everytime we add the file to another
	
		const int rebin = 20 ;
		
		ESD_CutLikeSign_montecarlo->Rebin(rebin);
		ESD_CutRefit_montecarlo->Rebin(rebin);
		ESD_CutKink_montecarlo->Rebin(rebin);
		ESD_CutdEdxSigmaElectronLine_montecarlo->Rebin(rebin);
		ESD_CutdEdxSigmaPionLine_montecarlo->Rebin(rebin);
		ESD_CutKaonRejectionLowP_montecarlo->Rebin(rebin);
		ESD_CutProtonRejectionLowP_montecarlo->Rebin(rebin);
		ESD_CutPionRejectionLowP_montecarlo->Rebin(rebin);

		ESD_CutGetOnFly_montecarlo->Rebin(rebin);
		ESD_CutNContributors_montecarlo->Rebin(rebin);
		ESD_CutPIDProb_montecarlo->Rebin(rebin);
		ESD_CutR_montecarlo->Rebin(rebin);
		ESD_CutLine_montecarlo->Rebin(rebin);
		ESD_CutZ_montecarlo->Rebin(rebin);
		ESD_CutMinNClsTPC_montecarlo->Rebin(rebin);
		ESD_CutSinglePt_montecarlo->Rebin(rebin);
		ESD_CutNDF_montecarlo->Rebin(rebin);
		ESD_CutChi2_montecarlo->Rebin(rebin);
		ESD_CutEta_montecarlo->Rebin(rebin);
		ESD_CutPt_montecarlo->Rebin(rebin);
		ESD_ConvGamma_Mass_montecarlo->Rebin(rebin);
		ESD_AllV0_montecarlo->Rebin(rebin);
			
	// Creation of histograms to compare cut with data before	

		TH1F* NoCutGetOnFly_montecarlo = (TH1F*) ESD_AllV0_montecarlo->Clone();

		TH1F* NoCutLikeSign_montecarlo = (TH1F*) NoCutGetOnFly_montecarlo->Clone();
		NoCutLikeSign_montecarlo->Add(ESD_CutGetOnFly_montecarlo,-1);
	
		TH1F* NoCutRefit_montecarlo = (TH1F*) NoCutLikeSign_montecarlo->Clone();		
		NoCutRefit_montecarlo->Add(ESD_CutLikeSign_montecarlo,-1);
	
		TH1F* NoCutKink_montecarlo = (TH1F*) NoCutRefit_montecarlo->Clone();		
		NoCutKink_montecarlo->Add(ESD_CutRefit_montecarlo,-1);
	
		TH1F* NoCutdEdxSigmaElectronLine_montecarlo = (TH1F*) NoCutKink_montecarlo->Clone();		
		NoCutdEdxSigmaElectronLine_montecarlo->Add(ESD_CutKink_montecarlo,-1);
	
		TH1F* NoCutdEdxSigmaPionLine_montecarlo = (TH1F*) NoCutdEdxSigmaElectronLine_montecarlo->Clone();		
		NoCutdEdxSigmaPionLine_montecarlo->Add(ESD_CutdEdxSigmaElectronLine_montecarlo,-1);
	
		TH1F* NoCutNContributors_montecarlo = (TH1F*) NoCutdEdxSigmaPionLine_montecarlo->Clone();		
		NoCutNContributors_montecarlo->Add(ESD_CutdEdxSigmaPionLine_montecarlo,-1);
	
		TH1F* NoCutKaonRejectionLowP_montecarlo = (TH1F*) NoCutdEdxSigmaPionLine_montecarlo->Clone();		
		NoCutKaonRejectionLowP_montecarlo->Add(ESD_CutdEdxSigmaPionLine_montecarlo,-1);
	
		TH1F* NoCutProtonRejectionLowP_montecarlo = (TH1F*) NoCutKaonRejectionLowP_montecarlo->Clone();		
		NoCutProtonRejectionLowP_montecarlo->Add(ESD_CutKaonRejectionLowP_montecarlo,-1);

		TH1F* NoCutPionRejectionLowP_montecarlo = (TH1F*) NoCutProtonRejectionLowP_montecarlo->Clone();		
		NoCutPionRejectionLowP_montecarlo->Add(ESD_CutProtonRejectionLowP_montecarlo,-1);

		TH1F* NoCutNContributors_montecarlo = (TH1F*) NoCutPionRejectionLowP_montecarlo->Clone();		
		NoCutNContributors_montecarlo->Add(ESD_CutPionRejectionLowP_montecarlo,-1);

	TH1F* NoCutPIDProb_montecarlo = (TH1F*) NoCutNContributors_montecarlo->Clone();		
	NoCutPIDProb_montecarlo->Add(ESD_CutNContributors_montecarlo,-1);

	
		TH1F* NoCutR_montecarlo = (TH1F*) NoCutPIDProb_montecarlo->Clone();		
		NoCutR_montecarlo->Add(ESD_CutPIDProb_montecarlo,-1);
	
		TH1F* NoCutLine_montecarlo = (TH1F*) NoCutR_montecarlo->Clone();		
		NoCutLine_montecarlo->Add(ESD_CutR_montecarlo,-1);
	
		TH1F* NoCutZ_montecarlo = (TH1F*) NoCutLine_montecarlo->Clone();		
		NoCutZ_montecarlo->Add(ESD_CutLine_montecarlo,-1);
	
		TH1F* NoCutMinNClsTPC_montecarlo = (TH1F*) NoCutZ_montecarlo->Clone();		
		NoCutMinNClsTPC_montecarlo->Add(ESD_CutZ_montecarlo,-1);

		TH1F* NoCutSinglePt_montecarlo = (TH1F*) NoCutMinNClsTPC_montecarlo->Clone();		
		NoCutSinglePt_montecarlo->Add(ESD_CutMinNClsTPC_montecarlo,-1);
	
		TH1F* NoCutNDF_montecarlo = (TH1F*) NoCutSinglePt_montecarlo->Clone();		
		NoCutNDF_montecarlo->Add(ESD_CutSinglePt_montecarlo,-1);

		TH1F* NoCutChi2_montecarlo = (TH1F*) NoCutNDF_montecarlo->Clone();		
		NoCutChi2_montecarlo->Add(ESD_CutNDF_montecarlo,-1);
	
		TH1F* NoCutEta_montecarlo = (TH1F*) NoCutChi2_montecarlo->Clone();		
		NoCutEta_montecarlo->Add(ESD_CutChi2_montecarlo,-1);
	
		TH1F* NoCutPt_montecarlo = (TH1F*) NoCutEta_montecarlo->Clone();		
		NoCutPt_montecarlo->Add(ESD_CutEta_montecarlo,-1);
	
		TH1F* Final_montecarlo = (TH1F*) NoCutPt_montecarlo->Clone();			
		Final_montecarlo->Add(ESD_CutPt_montecarlo,-1);
		
		Double_t NumberConvGamma_montecarlo = ESD_ConvGamma_Mass_montecarlo->Integral();
		Double_t NumberV0_montecarlo= ESD_AllV0_montecarlo->Integral();
		//Scaling
		ESD_CutLikeSign_montecarlo->Scale(normFac_montecarlo);
		ESD_CutRefit_montecarlo->Scale(normFac_montecarlo);
		ESD_CutKink_montecarlo->Scale(normFac_montecarlo);
		ESD_CutdEdxSigmaElectronLine_montecarlo->Scale(normFac_montecarlo);
		ESD_CutdEdxSigmaPionLine_montecarlo->Scale(normFac_montecarlo);
		ESD_CutKaonRejectionLowP_montecarlo->Scale(normFac_montecarlo);
		ESD_CutProtonRejectionLowP_montecarlo->Scale(normFac_montecarlo);
		ESD_CutPionRejectionLowP_montecarlo->Scale(normFac_montecarlo);
		ESD_CutGetOnFly_montecarlo->Scale(normFac_montecarlo);
		ESD_CutNContributors_montecarlo->Scale(normFac_montecarlo);	
		ESD_CutPIDProb_montecarlo->Scale(normFac_montecarlo);
		ESD_CutR_montecarlo->Scale(normFac_montecarlo);
		ESD_CutLine_montecarlo->Scale(normFac_montecarlo);
		ESD_CutZ_montecarlo->Scale(normFac_montecarlo);
		ESD_CutMinNClsTPC_montecarlo->Scale(normFac_montecarlo);
		ESD_CutSinglePt_montecarlo->Scale(normFac_montecarlo);
		ESD_CutNDF_montecarlo->Scale(normFac_montecarlo);
		ESD_CutChi2_montecarlo->Scale(normFac_montecarlo); 
		ESD_CutEta_montecarlo->Scale(normFac_montecarlo);
		ESD_CutPt_montecarlo->Scale(normFac_montecarlo);
		ESD_ConvGamma_Mass_montecarlo->Scale(normFac_montecarlo);
		NoCutPt_montecarlo->Scale(normFac_montecarlo);
		NoCutEta_montecarlo->Scale(normFac_montecarlo);
		NoCutChi2_montecarlo->Scale(normFac_montecarlo);
		NoCutNDF_montecarlo->Scale(normFac_montecarlo);
		NoCutZ_montecarlo->Scale(normFac_montecarlo);	
		NoCutMinNClsTPC_montecarlo->Scale(normFac_montecarlo);
		NoCutSinglePt_montecarlo->Scale(normFac_montecarlo);
		NoCutLine_montecarlo->Scale(normFac_montecarlo);
		NoCutR_montecarlo->Scale(normFac_montecarlo);;
		NoCutPIDProb_montecarlo->Scale(normFac_montecarlo);
		NoCutNContributors_montecarlo->Scale(normFac_montecarlo);
		NoCutGetOnFly_montecarlo->Scale(normFac_montecarlo);
		NoCutKaonRejectionLowP_montecarlo->Scale(normFac_montecarlo);
		NoCutProtonRejectionLowP_montecarlo->Scale(normFac_montecarlo);
		NoCutPionRejectionLowP_montecarlo->Scale(normFac_montecarlo);
		NoCutdEdxSigmaPionLine_montecarlo->Scale(normFac_montecarlo);
		NoCutdEdxSigmaElectronLine_montecarlo->Scale(normFac_montecarlo);
		NoCutKink_montecarlo->Scale(normFac_montecarlo);
		NoCutRefit_montecarlo->Scale(normFac_montecarlo);
		NoCutLikeSign_montecarlo->Scale(normFac_montecarlo);
		ESD_AllV0_montecarlo->Scale(normFac_montecarlo);
		Final_montecarlo->Scale(normFac_montecarlo);

		Double_t RealV0_montecarlo = 0. ;
		for(Int_t i=2;i<ESD_NumberOfV0_montecarlo->GetNbinsX()+1;i++){
		RealV0_montecarlo += ((ESD_NumberOfV0_montecarlo->GetBinCenter(i))* ESD_NumberOfV0_montecarlo->GetBinContent(i));
		}
	
		cout << "Number of V0's montecarlo:   " << NumberV0_montecarlo <<"   "  << RealV0_montecarlo << endl;
	}

	//------------------------------Def. of legend for all ------------------------------
	leg1 = new TLegend(0.55,0.79,0.9,0.9);
	leg1->AddEntry(ESD_CutLikeSign_data,("Data"),"l");
	if(MCfile != ""){
	leg1->AddEntry(ESD_CutLikeSign_montecarlo,("MC"),"l");}
	if(!Full){
		leg1->AddEntry(NoCutLikeSign_data,("Data without last cut"),"l");
		if(MCfile!= "") {	
			leg1->AddEntry(NoCutLikeSign_montecarlo,("MC without last cut"),"l");}
	}	
	if(Full){
		leg1->AddEntry(NoCutLikeSign_data,("Data without cut"),"l");
		if(MCfile!= "") {	
			leg1->AddEntry(NoCutLikeSign_montecarlo,("MC without cut"),"l");}
	}
	leg1->SetTextSize(0.04);
	
	
	// ----------------------- Start  ps file -----------------------------------------------------

	// ----------------------- Titlepage - Page 1  -----------------------------------------
	if(!SinglePlots) { 	
//		ps_characteristics->NewPage();
		
//		TCanvas * c0 = new TCanvas("c0","",700,1000);  // gives the page size
		
		title0 = new TPaveLabel(0.05,0.92,0.95,0.96,(Form("data: %s",data)));	
		title0->SetFillColor(16);
		title0->SetTextSize(0.25);
//		title0->Draw();
		
		if(MCfile != ""){
			title1 = new TPaveLabel(0.05,0.87,0.95,0.91,(Form("montecarlo: %s",MCfile)));	
			title1->SetFillColor(16);
			title1->SetTextSize(0.25);
//			title1->Draw();
			
		}
		
//		c0->Update();
		
		//------------ Track Cuts - Page 2 ------------------------------------
		
		ps_characteristics->NewPage();
		
		TCanvas * c1 = new TCanvas("c1","",700,1000);  // gives the page size
		pad_c1 = new TPad("pad_c1","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c1->SetFillColor(0);
		pad_c1->GetFrame()->SetFillColor(0);
		pad_c1->SetBorderMode(0);
		pad_c1->Divide(2,2);
		pad_c1->Draw();
	
		title0->Draw();
		if(MCfile != ""){title1->Draw();}
		

		pad_c1->cd(1)->SetLogy(1);
		pad_c1->cd(2)->SetLogy(1);
		pad_c1->cd(3)->SetLogy(1);
		pad_c1->cd(4)->SetLogy(1);

		//*************** CutGetOnFly
		pad_c1->cd(1);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutGetOnFly_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutGetOnFly_data, NoCutGetOnFly_data,
									 ESD_CutGetOnFly_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutGetOnFly_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutGetOnFly_data, NoCutGetOnFly_data,
									 ESD_CutGetOnFly_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		//************** CutLikeSign 
		pad_c1->cd(2);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutLikeSign_data, NoCutLikeSign_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutLikeSign_data, NoCutLikeSign_data,
									 ESD_CutLikeSign_montecarlo, NoCutLikeSign_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutLikeSign_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutLikeSign_data, NoCutGetOnFly_data,
									 ESD_CutLikeSign_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}

		//******************* CutRefit 
		pad_c1->cd(3);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutRefit_data, NoCutRefit_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutRefit_data, NoCutRefit_data,
									 ESD_CutRefit_montecarlo, NoCutRefit_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutRefit_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutRefit_data, NoCutGetOnFly_data,
									 ESD_CutRefit_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}

		//*************  CutKink
		pad_c1->cd(4);
	if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutKink_data, NoCutKink_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutKink_data, NoCutKink_data,
									 ESD_CutKink_montecarlo, NoCutKink_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutKink_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutKink_data, NoCutGetOnFly_data,
									 ESD_CutKink_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		
		c1->Update();

		// ----------ElectronSigma Selection/ Pion rejection - Page 3-------------------------------------------------------------------
		
		ps_characteristics->NewPage();
		
		TCanvas * c3 = new TCanvas("c3","",700,1000);  // gives the page size
		
		pad_c3 = new TPad("pad_c3","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c3->SetFillColor(0);
		pad_c3->GetFrame()->SetFillColor(0);
		pad_c3->SetBorderMode(0);
		pad_c3->Divide(1,2);
		pad_c3->Draw();
		
		if(MCfile != ""){title1->Draw();}
		title0->Draw();
		
		pad_c3->cd(1)->SetLogy(1);
		pad_c3->cd(2)->SetLogy(1);

		//****************** CutdEdxSigmaElectronLine 
		pad_c3->cd(1);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutdEdxSigmaElectronLine_data, NoCutdEdxSigmaElectronLine_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutdEdxSigmaElectronLine_data, NoCutdEdxSigmaElectronLine_data,
									 ESD_CutdEdxSigmaElectronLine_montecarlo, NoCutdEdxSigmaElectronLine_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutdEdxSigmaElectronLine_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutdEdxSigmaElectronLine_data, NoCutGetOnFly_data,
									 ESD_CutdEdxSigmaElectronLine_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}

		//*************** CutdEdxSigmaPionLine
		pad_c3->cd(2);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutdEdxSigmaPionLine_data, NoCutdEdxSigmaPionLine_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutdEdxSigmaPionLine_data, NoCutdEdxSigmaPionLine_data,
									 ESD_CutdEdxSigmaPionLine_montecarlo, NoCutdEdxSigmaPionLine_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutdEdxSigmaPionLine_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutdEdxSigmaPionLine_data, NoCutGetOnFly_data,
									 ESD_CutdEdxSigmaPionLine_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		
		c3->Update();
		
		// -------Selection in converted gammas - Page 4--------------------------------------------------------
		
		ps_characteristics->NewPage();
		
		TCanvas * c4 = new TCanvas("c4","",700,1000);  // gives the page size
		pad_c4 = new TPad("pad_c4","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c4->SetFillColor(0);
		pad_c4->GetFrame()->SetFillColor(0);
		pad_c4->SetBorderMode(0);
		pad_c4->Divide(2,3);
		pad_c4->Draw();
		
		title0->Draw();
		if(MCfile != ""){title1->Draw();}
		
		pad_c4->cd(1)->SetLogy(1);
		pad_c4->cd(2)->SetLogy(1);
		pad_c4->cd(3)->SetLogy(1);
		pad_c4->cd(4)->SetLogy(1);
		pad_c4->cd(5)->SetLogy(1);
		pad_c4->cd(6)->SetLogy(1);
		
		//***************** CutNContributors
		pad_c4->cd(1);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutNContributors_data, NoCutNContributors_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutNContributors_data, NoCutNContributors_data,
									 ESD_CutNContributors_montecarlo, NoCutNContributors_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutNContributors_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutNContributors_data, NoCutGetOnFly_data,
									 ESD_CutNContributors_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}

		//***************** CutPIDProb 
		pad_c4->cd(2);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutPIDProb_data, NoCutPIDProb_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutPIDProb_data, NoCutPIDProb_data,
									 ESD_CutPIDProb_montecarlo, NoCutPIDProb_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutPIDProb_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutPIDProb_data, NoCutGetOnFly_data,
									 ESD_CutPIDProb_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}

		//*************** CutR
		pad_c4->cd(3);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutR_data, NoCutR_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutR_data, NoCutR_data,
									 ESD_CutR_montecarlo, NoCutR_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutR_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutR_data, NoCutGetOnFly_data,
									 ESD_CutR_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}

		//****************** CutLine
		pad_c4->cd(4);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutLine_data, NoCutLine_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutLine_data, NoCutLine_data,
									 ESD_CutLine_montecarlo, NoCutLine_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutLine_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutLine_data, NoCutGetOnFly_data,
									 ESD_CutLine_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		
		//*************** CutZ
		pad_c4->cd(5);
		if(!Full){	
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutZ_data, NoCutZ_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutZ_data, NoCutZ_data,
									 ESD_CutZ_montecarlo, NoCutZ_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutZ_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutZ_data, NoCutGetOnFly_data,
									 ESD_CutZ_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}


		//******************* CutMinNClsTPC
		pad_c4->cd(6);
		if(!Full){
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutMinNClsTPC_data, NoCutMinNClsTPC_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutMinNClsTPC_data, NoCutMinNClsTPC_data,
									 ESD_CutMinNClsTPC_montecarlo, NoCutMinNClsTPC_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutMinNClsTPC_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutMinNClsTPC_data, NoCutGetOnFly_data,
									 ESD_CutMinNClsTPC_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		
		c4->Update();
				
		// -----------------------------------------------------------------------------
		
		ps_characteristics->NewPage();
		
		TCanvas * c5 = new TCanvas("c5","",700,1000);  // gives the page size
		pad_c5 = new TPad("pad_c5","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c5->SetFillColor(0);
		pad_c5->GetFrame()->SetFillColor(0);
		pad_c5->SetBorderMode(0);
		pad_c5->Divide(2,3);
		pad_c5->Draw();
		
		title0->Draw();
		if(MCfile != ""){title1->Draw();}
		
		pad_c5->cd(1)->SetLogy(1);
		pad_c5->cd(2)->SetLogy(1);
		pad_c5->cd(3)->SetLogy(1);
		pad_c5->cd(4)->SetLogy(1);
		pad_c5->cd(5)->SetLogy(1);
		pad_c5->cd(6)->SetLogy(0);

		//******************* CutSinglePt
		pad_c5->cd(1);
		if(!Full){
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutSinglePt_data, NoCutSinglePt_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutSinglePt_data, NoCutSinglePt_data,
									 ESD_CutSinglePt_montecarlo, NoCutSinglePt_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutSinglePt_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutSinglePt_data, NoCutGetOnFly_data,
									 ESD_CutSinglePt_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}

		//******************* CutNDF
		pad_c5->cd(2);
		if(!Full){
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutNDF_data, NoCutNDF_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutNDF_data, NoCutNDF_data,
									 ESD_CutNDF_montecarlo, NoCutNDF_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutNDF_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutNDF_data, NoCutGetOnFly_data,
									 ESD_CutNDF_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}

		
				
		//***************** CutChi2 
		pad_c5->cd(3);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutChi2_data, NoCutChi2_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutChi2_data, NoCutChi2_data,
									 ESD_CutChi2_montecarlo, NoCutChi2_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutChi2_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutChi2_data, NoCutGetOnFly_data,
									 ESD_CutChi2_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		
		//**************** CutEta
		pad_c5->cd(4);
		if(!Full){
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutEta_data, NoCutEta_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutEta_data, NoCutEta_data,
									 ESD_CutEta_montecarlo, NoCutEta_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutEta_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutEta_data, NoCutGetOnFly_data,
									 ESD_CutNDF_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}

		//********************** CutPt		
		pad_c5->cd(5);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutPt_data, NoCutPt_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutPt_data, NoCutPt_data,
									 ESD_CutPt_montecarlo, NoCutPt_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutPt_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutPt_data, NoCutGetOnFly_data,
									 ESD_CutPt_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
				
		c5->Update();





	// Rescaling plots
		ESD_CutLikeSign_data->Scale(1./normFac_data);
		ESD_CutRefit_data->Scale(1./normFac_data);
		ESD_CutKink_data->Scale(1./normFac_data);
		ESD_CutdEdxSigmaElectronLine_data->Scale(1./normFac_data);
		ESD_CutdEdxSigmaPionLine_data->Scale(1./normFac_data);
		ESD_CutGetOnFly_data->Scale(1./normFac_data);
		ESD_CutNContributors_data->Scale(1./normFac_data);	
		ESD_CutPIDProb_data->Scale(1./normFac_data);
		ESD_CutR_data->Scale(1./normFac_data);
		ESD_CutLine_data->Scale(1./normFac_data);
		ESD_CutZ_data->Scale(1./normFac_data);
		ESD_CutSinglePt_data->Scale(1./normFac_data);
		ESD_CutMinNClsTPC_data->Scale(1./normFac_data);
		ESD_CutNDF_data->Scale(1./normFac_data);
		ESD_CutChi2_data->Scale(1./normFac_data); 
		ESD_CutEta_data->Scale(1./normFac_data);
		ESD_CutPt_data->Scale(1./normFac_data);
		ESD_ConvGamma_Mass_data->Scale(1./normFac_data);
		NoCutPt_data->Scale(1./normFac_data);
		NoCutEta_data->Scale(1./normFac_data);
		NoCutChi2_data->Scale(1./normFac_data);
		NoCutNDF_data->Scale(1./normFac_data);
		NoCutZ_data->Scale(1./normFac_data);	
		NoCutSinglePt_data->Scale(1./normFac_data);
		NoCutMinNClsTPC_data->Scale(1./normFac_data);
		NoCutLine_data->Scale(1./normFac_data);
		NoCutR_data->Scale(1./normFac_data);;
		NoCutPIDProb_data->Scale(1./normFac_data);
		NoCutNContributors_data->Scale(1./normFac_data);
		NoCutGetOnFly_data->Scale(1./normFac_data);
		NoCutdEdxSigmaPionLine_data->Scale(1./normFac_data);
		NoCutdEdxSigmaElectronLine_data->Scale(1./normFac_data);
		NoCutKink_data->Scale(1./normFac_data);
		NoCutRefit_data->Scale(1./normFac_data);
		NoCutLikeSign_data->Scale(1./normFac_data);
		ESD_AllV0_data->Scale(1./normFac_data);
		Final_data->Scale(1./normFac_data);
		if(MCfile != ""){
			ESD_CutLikeSign_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutRefit_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutKink_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutdEdxSigmaElectronLine_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutdEdxSigmaPionLine_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutGetOnFly_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutNContributors_montecarlo->Scale(1./normFac_montecarlo);	
			ESD_CutPIDProb_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutR_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutLine_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutZ_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutSinglePt_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutMinNClsTPC_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutNDF_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutChi2_montecarlo->Scale(1./normFac_montecarlo); 
			ESD_CutEta_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutPt_montecarlo->Scale(1./normFac_montecarlo);
			ESD_ConvGamma_Mass_montecarlo->Scale(1./normFac_montecarlo);
			NoCutPt_montecarlo->Scale(1./normFac_montecarlo);
			NoCutEta_montecarlo->Scale(1./normFac_montecarlo);
			NoCutChi2_montecarlo->Scale(1./normFac_montecarlo);
			NoCutNDF_montecarlo->Scale(1./normFac_montecarlo);
			NoCutZ_montecarlo->Scale(1./normFac_montecarlo);	
			NoCutSinglePt_montecarlo->Scale(1./normFac_montecarlo);
			NoCutMinNClsTPC_montecarlo->Scale(1./normFac_montecarlo);
			NoCutLine_montecarlo->Scale(1./normFac_montecarlo);
			NoCutR_montecarlo->Scale(1./normFac_montecarlo);;
			NoCutPIDProb_montecarlo->Scale(1./normFac_montecarlo);
			NoCutNContributors_montecarlo->Scale(1./normFac_montecarlo);
			NoCutGetOnFly_montecarlo->Scale(1./normFac_montecarlo);
			NoCutdEdxSigmaPionLine_montecarlo->Scale(1./normFac_montecarlo);
			NoCutdEdxSigmaElectronLine_montecarlo->Scale(1./normFac_montecarlo);
			NoCutKink_montecarlo->Scale(1./normFac_montecarlo);
			NoCutRefit_montecarlo->Scale(1./normFac_montecarlo);
			NoCutLikeSign_montecarlo->Scale(1./normFac_montecarlo);
			ESD_AllV0_montecarlo->Scale(1./normFac_montecarlo);
			Final_montecarlo->Scale(1./normFac_montecarlo);
		}

		ps_characteristics->NewPage()	;

		TCanvas * c6 = new TCanvas("c6","",700,1000);  // gives the page size
		pad_c6 = new TPad("pad_c6","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c6->SetFillColor(0);
		pad_c6->GetFrame()->SetFillColor(0);
		pad_c6->SetBorderMode(0);
		pad_c6->Divide(2,2);
		pad_c6->Draw();
		
		title0->Draw();
		if(MCfile != ""){title1->Draw();}

		pad_c6->cd(1)->SetLogy(1);
		//******************** Comparison of Ratio of Cuts
		pad_c6->cd(1)	;
		const Int_t nx = 17;
		char *cuts[nx] = {"V0Finder","LikeSign", "Refit", "Kink", "dEdxElectronLine", "dEdxPionLine" , "NContributors", "PIDProb", "R", "Line", "Z", "Min", "SinglePt","NDF", "Chi2", "Eta", "Pt"};
		leg2 = new TLegend(0.6,0.82,0.9,0.9);
			leg2->AddEntry(ESD_CutLikeSign_data,("Data"),"l");
			if(MCfile != ""){
			leg2->AddEntry(ESD_CutLikeSign_montecarlo,("MC"),"l");}
		leg2->SetFillColor(0);

		Double_t maxRange = ESD_CutdEdxSigmaElectronLine_data->Integral()/NumberV0_data *100;
		if(MCfile != ""){
			if(maxRange < ESD_CutdEdxSigmaElectronLine_montecarlo->Integral()/NumberV0_montecarlo *100){
				maxRange = ESD_CutdEdxSigmaElectronLine_montecarlo->Integral()/NumberV0_montecarlo *100;
			}
		}		

		Output_data = new TFile("Output_data.root","RECREATE");
		
		TH1F *ComCuts_data = new TH1F("ComCuts_data","Comparison of Cuts",3,0,3);
		ComCuts_data->SetStats(0);
		ComCuts_data->SetYTitle("ratio Cut/V0") ;
		ComCuts_data->SetBit(TH1::kCanRebin);
		ComCuts_data->GetYaxis()->SetRangeUser(0.001, 5 *maxRange);
		ComCuts_data->GetYaxis()->SetTitleOffset(1.35);
		ComCuts_data->Fill(cuts[0],ESD_CutGetOnFly_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[1],ESD_CutLikeSign_data->Integral()/NumberV0_data *100);
		ComCuts_data->Fill(cuts[2],ESD_CutRefit_data->Integral()/NumberV0_data *100);
		ComCuts_data->Fill(cuts[3],ESD_CutKink_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[4],ESD_CutdEdxSigmaElectronLine_data->Integral()/NumberV0_data *100);
		ComCuts_data->Fill(cuts[5],ESD_CutdEdxSigmaPionLine_data->Integral()/NumberV0_data *100);
		ComCuts_data->Fill(cuts[6],ESD_CutNContributors_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[7],ESD_CutPIDProb_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[8],ESD_CutR_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[9],ESD_CutLine_data->Integral()/NumberV0_data *100);		 		
		ComCuts_data->Fill(cuts[10],ESD_CutZ_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[11],ESD_CutMinNClsTPC_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[12],ESD_CutSinglePt_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[13],ESD_CutNDF_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[14],ESD_CutChi2_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[15],ESD_CutEta_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[16],ESD_CutPt_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->LabelsDeflate();
		ComCuts_data->Draw("e,hist");
                    ComCuts_data->Write();
                    Output_data->Write();

		if(MCfile != ""){
			Output_montecarlo = new TFile("Output_montecarlo.root","RECREATE");
			TH1F *ComCuts_montecarlo = new TH1F("ComCuts_montecarlo","Comparison of Cuts",3,0,3);
			ComCuts_montecarlo->SetStats(0);
			ComCuts_montecarlo->SetBit(TH1::kCanRebin);
			ComCuts_montecarlo->SetLineColor(2);
			ComCuts_montecarlo->Fill(cuts[0],ESD_CutGetOnFly_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[1],ESD_CutLikeSign_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts_montecarlo->Fill(cuts[2],ESD_CutRefit_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts_montecarlo->Fill(cuts[3],ESD_CutKink_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[4],ESD_CutdEdxSigmaElectronLine_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts_montecarlo->Fill(cuts[5],ESD_CutdEdxSigmaPionLine_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts_montecarlo->Fill(cuts[6],ESD_CutNContributors_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[7],ESD_CutPIDProb_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[8],ESD_CutR_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[9],ESD_CutLine_montecarlo->Integral()/NumberV0_montecarlo *100);		 		
			ComCuts_montecarlo->Fill(cuts[10],ESD_CutZ_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[11],ESD_CutMinNClsTPC_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[12],ESD_CutSinglePt_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[13],ESD_CutNDF_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[14],ESD_CutChi2_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[15],ESD_CutEta_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[16],ESD_CutPt_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->LabelsDeflate();
			ComCuts_montecarlo->Draw("e,hist,same");
			ComCuts_montecarlo->Write();
			Output_montecarlo->Write();
		}
		leg2->Draw();

		pad_c6->cd(2)->SetLogy(0);

		//************************* Remaining V0 candidates after Cuts
		pad_c6->cd(2)	;
		const Int_t nx = 18;
		char *cuts2[nx] = {"No Cut", "V0Finder", "LikeSign", "Refit", "Kink", "dEdxElectronLine", "dEdxPionLine",  "NContributors", "PIDProb", "R", "Line", "Z","MinClsTPC", "SinglePt","NDF", "Chi2", "Eta", "Pt"};
		Double_t maxRange = NoCutLikeSign_data->Integral()/NumberV0_data*100;
		if(MCfile != ""){
			if(maxRange < NoCutLikeSign_montecarlo->Integral()/NumberV0_data*100){
				maxRange = NoCutLikeSign_montecarlo->Integral()/NumberV0_data*100;
			}
		}		
		
		Output_data = new TFile("Output_data.root","UPDATE");
		
		TH1F *ComCuts2_data = new TH1F("ComCuts2_data","Remaining V0 candidates after Cuts ",3,0,3);
		ComCuts2_data->SetYTitle("ratio Cut/V0");
		ComCuts2_data->SetStats(0);
		ComCuts2_data->SetBit(TH1::kCanRebin);
		ComCuts2_data->GetYaxis()->SetRangeUser(0.1, 110);
		ComCuts2_data->GetYaxis()->SetTitleOffset(1.35);
		ComCuts2_data->Fill(cuts2[0],NoCutGetOnFly_data->Integral()/NumberV0_data *100);
		cout << "Number NoGetOnFly   " << NoCutGetOnFly_data->Integral()/NumberV0_data *100 << endl;		 
		ComCuts2_data->Fill(cuts2[1],NoCutLikeSign_data->Integral()/NumberV0_data *100);
		ComCuts2_data->Fill(cuts2[2],NoCutRefit_data->Integral()/NumberV0_data *100);
		ComCuts2_data->Fill(cuts2[3],NoCutKink_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[4],NoCutdEdxSigmaElectronLine_data->Integral()/NumberV0_data *100);
		ComCuts2_data->Fill(cuts2[5],NoCutdEdxSigmaPionLine_data->Integral()/NumberV0_data *100);
		ComCuts2_data->Fill(cuts2[6],NoCutNContributors_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[7],NoCutPIDProb_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[8],NoCutR_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[9],NoCutLine_data->Integral()/NumberV0_data *100);		 		
		ComCuts2_data->Fill(cuts2[10],NoCutZ_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[11],NoCutMinNClsTPC_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[12],NoCutSinglePt_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[13],NoCutNDF_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[14],NoCutChi2_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[15],NoCutEta_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[16],NoCutPt_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[17],Final_data->Integral()/NumberV0_data *100);			
		ComCuts2_data->Fill("RealLeftOver",ESD_ConvGamma_Mass_data->Integral()/NumberV0_data *100);	
		if(MC2){ ComCuts2_data->Fill("TruePhotons",0);}
		if(MC1){ ComCuts2_data->Fill("TruePhotons",TRUE_Conversion_R_data->Integral()/NumberV0_data *100);}
		ComCuts2_data->LabelsDeflate();
		ComCuts2_data->Draw("e,hist");

		ComCuts2_data->Write();
		Output_data->Write();

		if(MCfile != ""){
			Output_montecarlo = new TFile("Output_montecarlo.root","UPDATE");
			TH1F *ComCuts2_montecarlo = new TH1F("ComCuts2_montecarlo","Remaining V0 candidates after Cuts",3,0,3);
			ComCuts2_montecarlo->SetStats(0);
			ComCuts2_montecarlo->SetBit(TH1::kCanRebin);
			ComCuts2_montecarlo->SetLineColor(2);
			ComCuts2_montecarlo->Fill(cuts2[0],NoCutGetOnFly_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[1],NoCutLikeSign_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts2_montecarlo->Fill(cuts2[2],NoCutRefit_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts2_montecarlo->Fill(cuts2[3],NoCutKink_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[4],NoCutdEdxSigmaElectronLine_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts2_montecarlo->Fill(cuts2[5],NoCutdEdxSigmaPionLine_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts2_montecarlo->Fill(cuts2[6],NoCutNContributors_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[7],NoCutPIDProb_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[8],NoCutR_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[9],NoCutLine_montecarlo->Integral()/NumberV0_montecarlo *100);		 		
			ComCuts2_montecarlo->Fill(cuts2[10],NoCutZ_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[11],NoCutMinNClsTPC_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[12],NoCutSinglePt_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[13],NoCutNDF_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[14],NoCutChi2_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[15],NoCutEta_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[16],NoCutPt_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[17],Final_montecarlo->Integral()/NumberV0_montecarlo *100);			
			ComCuts2_montecarlo->Fill("RealLeftOver",ESD_ConvGamma_Mass_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			if(MC2){ ComCuts2_montecarlo->Fill("TruePhotons",TRUE_Conversion_R_montecarlo->Integral()/NumberV0_montecarlo *100);}
			ComCuts2_montecarlo->LabelsDeflate();
			ComCuts2_montecarlo->Draw("e,hist,same");
			ComCuts2_montecarlo->Write();
			Output_montecarlo->Write();
		}		
		leg2->Draw();
		
	if(MCfile != ""){
		if(MC2){	
		
		Double_t truePhotons = TRUE_Conversion_R_montecarlo->Integral();
		Double_t measPhotons = ESD_Conversion_R_montecarlo->Integral();
		cout << "True Photons     " << truePhotons << endl;
		cout << "Simulated Photons     " << measPhotons << endl;



		pad_c6->cd(3)->SetLogy(0);
	
		pad_c6->cd(3);
		DrawAutoGammaHistos( ESD_Conversion_R_montecarlo,
						 TRUE_Conversion_R_montecarlo,  
						 "Comparison True Gammas and MC Gammas","R [cm]",StandardYAxis,
						 kFALSE,2 ,0.0000002,
						 kFALSE,0. ,0.,
						 kTRUE, 0.,180.);		
		}
	}
		
		
		c6->Update();
		ps_characteristics->Close();		


// ------------------------------ Statistics -------------------------------------------------------------------------		
/* In this part the output-file for the Data is generated, but as it only makes sense for compared outputs, 
   this section is only executed if a second input is there. 
*/		
		if(MCfile != ""){
		fstream outl;
    		outl.open(outlname, ios::out);
		outl << "data :\t" << data  << endl;
		outl << "montecarlo :\t" << MCfile << endl;
		outl << "------------------------------------------------------------------------------------------" << endl;
		outl << endl;
		outl << "\t data \t montecarlo " << endl;
		outl << "Number of events" << "\t" << nGoodEvents_data << "\t" << nGoodEvents_montecarlo << endl;
		outl << "Number of triggers" << "\t" << nGoodTrig_data << "\t" << nGoodTrig_montecarlo << endl;
		outl << "Number reconstructed gammas"<< "\t" << nRecGamma_data << "\t" << nRecGamma_montecarlo <<endl;
		outl << "Number of V0's from AllV0s" << "\t" << NumberV0_data <<"\t" <<  NumberV0_montecarlo << endl;
		outl << "Number of V0's real" << "\t" << RealV0_data << "\t" << RealV0_montecarlo << endl;
		outl << endl;
		if(MC2){
		outl << "True Photons     " << truePhotons << endl;
		outl << "Simulated Photons     " << measPhotons << endl;
		}		
		outl << endl;
		outl << "------------------------------------------------------------------------------------------" << endl;
		outl << "Applied Cut \t" << "data entries \t" << "data norm \t" << "montecarlo entries \t" << "montecarlo norm \t" << endl;
		outl << "------------------------------------------------------------------------------------------" << endl;
		outl << "GetOnFly \t" << ESD_CutGetOnFly_data->Integral() << "\t" << ESD_CutGetOnFly_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutGetOnFly_montecarlo->Integral() << "\t" << ESD_CutGetOnFly_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "LikeSign \t" << ESD_CutLikeSign_data->Integral() << "\t" << ESD_CutLikeSign_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutLikeSign_montecarlo->Integral() << "\t" << ESD_CutLikeSign_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "Refit \t" << ESD_CutRefit_data->Integral() << "\t" << ESD_CutRefit_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutRefit_montecarlo->Integral() << "\t" << ESD_CutRefit_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "Kink \t" << ESD_CutKink_data->Integral() << "\t" << ESD_CutKink_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutKink_montecarlo->Integral() << "\t" << ESD_CutKink_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "dEdxSigmaElectronLine \t" << ESD_CutdEdxSigmaElectronLine_data->Integral() << "\t" << ESD_CutdEdxSigmaElectronLine_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutdEdxSigmaElectronLine_montecarlo->Integral() << "\t" << ESD_CutdEdxSigmaElectronLine_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "dEdxSigmaPionLine \t" << ESD_CutdEdxSigmaPionLine_data->Integral() << "\t" << ESD_CutdEdxSigmaPionLine_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutdEdxSigmaPionLine_montecarlo->Integral() << "\t" << ESD_CutdEdxSigmaPionLine_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NContributors \t" << ESD_CutNContributors_data->Integral() << "\t" << ESD_CutNContributors_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutNContributors_montecarlo->Integral() << "\t" << ESD_CutNContributors_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "PIDProb \t" << ESD_CutPIDProb_data->Integral() << "\t" << ESD_CutPIDProb_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutPIDProb_montecarlo->Integral() << "\t" << ESD_CutPIDProb_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "R \t" << ESD_CutR_data->Integral() << "\t" << ESD_CutR_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutR_montecarlo->Integral() << "\t" << ESD_CutR_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "Line \t" << ESD_CutLine_data->Integral() << "\t" << ESD_CutLine_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutLine_montecarlo->Integral() << "\t" << ESD_CutLine_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "Z \t" << ESD_CutZ_data->Integral() << "\t" << ESD_CutZ_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutZ_montecarlo->Integral() << "\t" << ESD_CutZ_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NDF \t" << ESD_CutNDF_data->Integral() << "\t" << ESD_CutNDF_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutNDF_montecarlo->Integral() << "\t" << ESD_CutNDF_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "Chi2 \t" << ESD_CutChi2_data->Integral() << "\t" << ESD_CutChi2_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutChi2_montecarlo->Integral() << "\t" << ESD_CutChi2_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "Eta \t" << ESD_CutEta_data->Integral() << "\t" << ESD_CutEta_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutEta_montecarlo->Integral() << "\t" << ESD_CutEta_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "Pt \t" << ESD_CutPt_data->Integral() << "\t" << ESD_CutPt_data->Integral()/NumberV0_data * 100 << "\t" << ESD_CutPt_montecarlo->Integral() << "\t" << ESD_CutPt_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "------------------------------------------------------------------------------------------" << endl;
		outl << endl;
		outl << endl;
		outl << "NoCutHistos" << endl;
		outl << "------------------------------------------------------------------------------------------" << endl;
		outl << "NoGetOnFly \t" << NoCutGetOnFly_data->Integral() << "\t" << NoCutGetOnFly_data->Integral()/NumberV0_data * 100 << "\t" << NoCutGetOnFly_montecarlo->Integral() << "\t" << NoCutGetOnFly_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NoLikeSign \t" << NoCutLikeSign_data->Integral() << "\t" << NoCutLikeSign_data->Integral()/NumberV0_data * 100 << "\t" << NoCutLikeSign_montecarlo->Integral() << "\t" << NoCutLikeSign_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NoRefit \t" << NoCutRefit_data->Integral() << "\t" << NoCutRefit_data->Integral()/NumberV0_data * 100 << "\t" << NoCutRefit_montecarlo->Integral() << "\t" << NoCutRefit_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NoKink \t" << NoCutKink_data->Integral() << "\t" << NoCutKink_data->Integral()/NumberV0_data * 100 << "\t" << NoCutKink_montecarlo->Integral() << "\t" << NoCutKink_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NodEdxSigmaElectronLine \t" << NoCutdEdxSigmaElectronLine_data->Integral() << "\t" << NoCutdEdxSigmaElectronLine_data->Integral()/NumberV0_data * 100 << "\t" << NoCutdEdxSigmaElectronLine_montecarlo->Integral() << "\t" << NoCutdEdxSigmaElectronLine_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NodEdxSigmaPionLine \t" << NoCutdEdxSigmaPionLine_data->Integral() << "\t" << NoCutdEdxSigmaPionLine_data->Integral()/NumberV0_data * 100 << "\t" << NoCutdEdxSigmaPionLine_montecarlo->Integral() << "\t" << NoCutdEdxSigmaPionLine_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NoNContributors \t" << NoCutNContributors_data->Integral() << "\t" << NoCutNContributors_data->Integral()/NumberV0_data * 100 << "\t" << NoCutNContributors_montecarlo->Integral() << "\t" << NoCutNContributors_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NoPIDProb \t" << NoCutPIDProb_data->Integral() << "\t" << NoCutPIDProb_data->Integral()/NumberV0_data * 100 << "\t" << NoCutPIDProb_montecarlo->Integral() << "\t" << NoCutPIDProb_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NoR \t" << NoCutR_data->Integral() << "\t" << NoCutR_data->Integral()/NumberV0_data * 100 << "\t" << NoCutR_montecarlo->Integral() << "\t" << NoCutR_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NoLine \t" << NoCutLine_data->Integral() << "\t" << NoCutLine_data->Integral()/NumberV0_data * 100 << "\t" << NoCutLine_montecarlo->Integral() << "\t" << NoCutLine_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NoZ \t" << NoCutZ_data->Integral() << "\t" << NoCutZ_data->Integral()/NumberV0_data * 100 << "\t" << NoCutZ_montecarlo->Integral() << "\t" << NoCutZ_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NoNDF \t" << NoCutNDF_data->Integral() << "\t" << NoCutNDF_data->Integral()/NumberV0_data * 100 << "\t" << NoCutNDF_montecarlo->Integral() << "\t" << NoCutNDF_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NoChi2 \t" << NoCutChi2_data->Integral() << "\t" << NoCutChi2_data->Integral()/NumberV0_data * 100 << "\t" << NoCutChi2_montecarlo->Integral() << "\t" << NoCutChi2_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NoEta \t" << NoCutEta_data->Integral() << "\t" << NoCutEta_data->Integral()/NumberV0_data * 100 << "\t" << NoCutEta_montecarlo->Integral() << "\t" << NoCutEta_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "NoPt \t" << NoCutPt_data->Integral() << "\t" << NoCutPt_data->Integral()/NumberV0_data * 100 << "\t" << NoCutPt_montecarlo->Integral() << "\t" << NoCutPt_montecarlo->Integral()/NumberV0_montecarlo *100 << endl;
		outl << "------------------------------------------------------------------------------------------" << endl;


		outl.close();
		}
	
	}
//------------------------------------ Single-Plots --------------------------------------------------------------------
	if(SinglePlots){
// -----------------------------------------CutLikeSign ---------------------
		TCanvas * c_SinglePlot_1 = new TCanvas("c_SinglePlot_1","",1000,1000);  // gives the page size
		c_SinglePlot_1->SetLogy(1);	
		c_SinglePlot_1->cd();
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutLikeSign_data, NoCutLikeSign_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutLikeSign_data, NoCutLikeSign_data,
									 ESD_CutLikeSign_montecarlo, NoCutLikeSign_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutLikeSign_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutLikeSign_data, NoCutGetOnFly_data,
									 ESD_CutLikeSign_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_1->Update();	
		c_SinglePlot_1->SaveAs(Form("%s%s/CutLikeSign2.%s",path,suffix,suffix));
		delete c_SinglePlot_1;
//-------------------------------------------- CutRefit -------------------------------------------------------

		TCanvas * c_SinglePlot_2 = new TCanvas("c_SinglePlot_2","",1000,1000);  // gives the page size
		c_SinglePlot_2->SetLogy(1);
		c_SinglePlot_2->cd();
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutRefit_data, NoCutRefit_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutRefit_data, NoCutRefit_data,
									 ESD_CutRefit_montecarlo, NoCutRefit_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutRefit_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutRefit_data, NoCutGetOnFly_data,
									 ESD_CutRefit_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_2->Update();	
		c_SinglePlot_2->SaveAs(Form("%s%s/CutRefit2.%s",path,suffix,suffix));
		delete c_SinglePlot_2;
//-------------------------------------------- CutKink --------------------------------------------------------
		TCanvas * c_SinglePlot_3 = new TCanvas("c_SinglePlot_3","",1000,1000);  // gives the page size
		c_SinglePlot_3->SetLogy(1);
		c_SinglePlot_3->cd();
	if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutKink_data, NoCutKink_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutKink_data, NoCutKink_data,
									 ESD_CutKink_montecarlo, NoCutKink_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutKink_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutKink_data, NoCutGetOnFly_data,
									 ESD_CutKink_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_3->Update();	
		c_SinglePlot_3->SaveAs(Form("%s%s/CutKink2.%s",path,suffix,suffix));
		delete c_SinglePlot_3;
//-------------------------------------------- CutdEdxElectronLine --------------------------------------------
		


		TCanvas * c_SinglePlot_4 = new TCanvas("c_SinglePlot_4","",1000,1000);  // gives the page size
		c_SinglePlot_4->SetLogy(1);
		c_SinglePlot_4->cd();
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutdEdxSigmaElectronLine_data, NoCutdEdxSigmaElectronLine_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutdEdxSigmaElectronLine_data, NoCutdEdxSigmaElectronLine_data,
									 ESD_CutdEdxSigmaElectronLine_montecarlo, NoCutdEdxSigmaElectronLine_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutdEdxSigmaElectronLine_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutdEdxSigmaElectronLine_data, NoCutGetOnFly_data,
									 ESD_CutdEdxSigmaElectronLine_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_4->Update();	
		c_SinglePlot_4->SaveAs(Form("%s%s/CutdEdxElectronLine2.%s",path,suffix,suffix));
		delete c_SinglePlot_4;
//-------------------------------------------- CutdEdxPionLine ------------------------------------------------
		TCanvas * c_SinglePlot_5 = new TCanvas("c_SinglePlot_5","",1000,1000);  // gives the page size
		c_SinglePlot_5->SetLogy(1);
		c_SinglePlot_5->cd();
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutdEdxSigmaPionLine_data, NoCutdEdxSigmaPionLine_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutdEdxSigmaPionLine_data, NoCutdEdxSigmaPionLine_data,
									 ESD_CutdEdxSigmaPionLine_montecarlo, NoCutdEdxSigmaPionLine_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutdEdxSigmaPionLine_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutdEdxSigmaPionLine_data, NoCutGetOnFly_data,
									 ESD_CutdEdxSigmaPionLine_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_5->Update();	
		c_SinglePlot_5->SaveAs(Form("%s%s/CutdEdxPionLine2.%s",path,suffix,suffix));
		delete c_SinglePlot_5;
//-------------------------------------------- CutGetOnFly ----------------------------------------------------
		TCanvas * c_SinglePlot_6 = new TCanvas("c_SinglePlot_6","",1000,1000);  // gives the page size
		c_SinglePlot_6->SetLogy(1);
		c_SinglePlot_6->cd();
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutGetOnFly_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutGetOnFly_data, NoCutGetOnFly_data,
									 ESD_CutGetOnFly_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutGetOnFly_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutGetOnFly_data, NoCutGetOnFly_data,
									 ESD_CutGetOnFly_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_6->Update();	
		c_SinglePlot_6->SaveAs(Form("%s%s/CutGetOnFly2.%s",path,suffix,suffix));
		delete c_SinglePlot_6;
//-------------------------------------------- CutNContributors ----------------------------------------------
		TCanvas * c_SinglePlot_7 = new TCanvas("c_SinglePlot_7","",1000,1000);  // gives the page size
		c_SinglePlot_7->cd();
		c_SinglePlot_7->SetLogy(1);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutNContributors_data, NoCutNContributors_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutNContributors_data, NoCutNContributors_data,
									 ESD_CutNContributors_montecarlo, NoCutNContributors_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutNContributors_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutNContributors_data, NoCutGetOnFly_data,
									 ESD_CutNContributors_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_7->Update();	
		c_SinglePlot_7->SaveAs(Form("%s%s/CutNContributors2.%s",path,suffix,suffix));
		delete c_SinglePlot_7;
//-------------------------------------------- CutPIDProb ----------------------------------------------------
		TCanvas * c_SinglePlot_8 = new TCanvas("c_SinglePlot_8","",1000,1000);  // gives the page size
		c_SinglePlot_8->SetLogy(1);
		c_SinglePlot_8->cd();
	if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutPIDProb_data, NoCutPIDProb_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutPIDProb_data, NoCutPIDProb_data,
									 ESD_CutPIDProb_montecarlo, NoCutPIDProb_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutPIDProb_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutPIDProb_data, NoCutGetOnFly_data,
									 ESD_CutPIDProb_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_8->Update();	
		c_SinglePlot_8->SaveAs(Form("%s%s/CutPIDProb2.%s",path,suffix,suffix));
		delete c_SinglePlot_8;
//-------------------------------------------- CutR ----------------------------------------------------------
		TCanvas * c_SinglePlot_9 = new TCanvas("c_SinglePlot_9","",1000,1000);  // gives the page size
		c_SinglePlot_9->SetLogy(1);
		c_SinglePlot_9->cd();
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutR_data, NoCutR_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutR_data, NoCutR_data,
									 ESD_CutR_montecarlo, NoCutR_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutR_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutR_data, NoCutGetOnFly_data,
									 ESD_CutR_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_9->Update();	
		c_SinglePlot_9->SaveAs(Form("%s%s/CutR2.%s",path,suffix,suffix));
		delete c_SinglePlot_9;
//-------------------------------------------- CutLine ------------------------------------------------------- 
		TCanvas * c_SinglePlot_10 = new TCanvas("c_SinglePlot_10","",1000,1000);  // gives the page size	
		c_SinglePlot_10->SetLogy(1);
		c_SinglePlot_10->cd();
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutLine_data, NoCutLine_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutLine_data, NoCutLine_data,
									 ESD_CutLine_montecarlo, NoCutLine_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutLine_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutLine_data, NoCutGetOnFly_data,
									 ESD_CutLine_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
				c_SinglePlot_10->Update();	
		c_SinglePlot_10->SaveAs(Form("%s%s/CutLine2.%s",path,suffix,suffix));
		delete c_SinglePlot_10;

//-------------------------------------------- CutZ ------------------------------------------------------------
		TCanvas * c_SinglePlot_11 = new TCanvas("c_SinglePlot_11","",1000,1000);  // gives the page size
		c_SinglePlot_11->SetLogy(1);
		c_SinglePlot_11->cd();
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutZ_data, NoCutZ_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutZ_data, NoCutZ_data,
									 ESD_CutZ_montecarlo, NoCutZ_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutZ_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutZ_data, NoCutGetOnFly_data,
									 ESD_CutZ_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_11->Update();	
		c_SinglePlot_11->SaveAs(Form("%s%s/CutZ.%s",path,suffix,suffix));
		delete c_SinglePlot_11;

//-------------------------------------------- CutMinNClsTPC ------------------------------------------------------------
		TCanvas * c_SinglePlot_11 = new TCanvas("c_SinglePlot_11","",1000,1000);  // gives the page siMinNClsTPCe
		c_SinglePlot_11->SetLogy(1);
		c_SinglePlot_11->cd();
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutMinNClsTPC_data, NoCutMinNClsTPC_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutMinNClsTPC_data, NoCutMinNClsTPC_data,
									 ESD_CutMinNClsTPC_montecarlo, NoCutMinNClsTPC_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutMinNClsTPC_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutMinNClsTPC_data, NoCutGetOnFly_data,
									 ESD_CutMinNClsTPC_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_11->Update();	
		c_SinglePlot_11->SaveAs(Form("%s%s/CutMinNClsTPC.%s",path,suffix,suffix));
		delete c_SinglePlot_11;

//-------------------------------------------- CutSinglePt ------------------------------------------------------------
		TCanvas * c_SinglePlot_11 = new TCanvas("c_SinglePlot_11","",1000,1000);  // gives the page siSinglePte
		c_SinglePlot_11->SetLogy(1);
		c_SinglePlot_11->cd();
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutSinglePt_data, NoCutSinglePt_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutSinglePt_data, NoCutSinglePt_data,
									 ESD_CutSinglePt_montecarlo, NoCutSinglePt_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutSinglePt_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutSinglePt_data, NoCutGetOnFly_data,
									 ESD_CutSinglePt_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_11->Update();	
		c_SinglePlot_11->SaveAs(Form("%s%s/CutSinglePt.%s",path,suffix,suffix));
		delete c_SinglePlot_11;

//-------------------------------------------- CutNDF ----------------------------------------------------------
		TCanvas * c_SinglePlot_12 = new TCanvas("c_SinglePlot_12","",1000,1000);  // gives the page size
		c_SinglePlot_12->cd();
		c_SinglePlot_12->SetLogy(1);
		if(!Full){				
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutNDF_data, NoCutNDF_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutNDF_data, NoCutNDF_data,
									 ESD_CutNDF_montecarlo, NoCutNDF_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutNDF_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutNDF_data, NoCutGetOnFly_data,
									 ESD_CutNDF_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_12->Update();	
		c_SinglePlot_12->SaveAs(Form("%s%s/CutNDF2.%s",path,suffix,suffix));
		delete c_SinglePlot_12;
//-------------------------------------------- CutChi2 ---------------------------------------------------------
		TCanvas * c_SinglePlot_13 = new TCanvas("c_SinglePlot_13","",1000,1000);  // gives the page size
		c_SinglePlot_13->SetLogy(1);
		c_SinglePlot_13->cd();
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutChi2_data, NoCutChi2_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutChi2_data, NoCutChi2_data,
									 ESD_CutChi2_montecarlo, NoCutChi2_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutChi2_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutChi2_data, NoCutGetOnFly_data,
									 ESD_CutChi2_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_13->Update();	
		c_SinglePlot_13->SaveAs(Form("%s%s/CutChi2_2.%s",path,suffix,suffix));
		delete c_SinglePlot_13;
//-------------------------------------------- CutEta ----------------------------------------------------------
		TCanvas * c_SinglePlot_14 = new TCanvas("c_SinglePlot_14","",1000,1000);  // gives the page size
		c_SinglePlot_14->cd();
		c_SinglePlot_14->SetLogy(1);
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutEta_data, NoCutEta_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutEta_data, NoCutEta_data,
									 ESD_CutEta_montecarlo, NoCutEta_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutEta_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutEta_data, NoCutGetOnFly_data,
									 ESD_CutNDF_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_14->Update();	
		c_SinglePlot_14->SaveAs(Form("%s%s/CutEta2.%s",path,suffix,suffix));
		delete c_SinglePlot_14;

//-------------------------------------------- CutPt ----------------------------------------------------------
		TCanvas * c_SinglePlot_15 = new TCanvas("c_SinglePlot_15","",1000,1000);  // gives the page size
		c_SinglePlot_15->SetLogy(1);
		c_SinglePlot_15->cd();
		if(!Full){		
			if(MCfile == ""){
					DrawCutGammaHisto( ESD_CutPt_data, NoCutPt_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutPt_data, NoCutPt_data,
									 ESD_CutPt_montecarlo, NoCutPt_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without last cut", "MC without last cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}
		}else{
			if(MCfile == ""){
					DrawCutGammaHisto( 	ESD_CutPt_data, NoCutGetOnFly_data,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}else{
					DrawCutGammaHistos( ESD_CutPt_data, NoCutGetOnFly_data,
									 ESD_CutPt_montecarlo, NoCutGetOnFly_montecarlo,
									 "","InvMass [GeV]",StandardYAxis, "Data without any cut", "MC without any cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kTRUE, 0.,0.5);
			}	
		}
		c_SinglePlot_15->Update();	
		c_SinglePlot_15->SaveAs(Form("%s%s/CutPt2.%s",path,suffix,suffix));
		delete c_SinglePlot_15;

//-------------------------------------------- InvMassdistribution of leftover ---------------------------------------------
		TCanvas * c_SinglePlot_15 = new TCanvas("c_SinglePlot_15","",1000,1000);  // gives the page size
		c_SinglePlot_15->SetLogy(1);
		c_SinglePlot_15->cd();
					DrawCutGammaHisto( 	NoCutGetOnFly_data, ESD_ConvGamma_Mass_data, 
									 "","InvMass [GeV]",StandardYAxis, "Data w/o cut",
									 kFALSE, 1.5,0.00001,
									 kTRUE,0.000001 ,1,
									 kFALSE, 0.,0);
					NoCutGetOnFly_data->SetLineWidth(2);
					NoCutGetOnFly_data->SetLineColor(1);
					NoCutGetOnFly_data->SetTitle();
					ESD_ConvGamma_Mass_data->SetTitle();
					ESD_ConvGamma_Mass_data->SetLineColor(2);
					ESD_ConvGamma_Mass_data->SetLineWidth(2);
					ESD_ConvGamma_Mass_data->Draw("same,hist");
				leg1 = new TLegend( 0.6,0.82,0.92,0.9);
				leg1->SetTextSize(0.03);			
				leg1->SetFillColor(0);
				leg1->AddEntry(NoCutGetOnFly_data,("before cuts"));
				leg1->AddEntry(ESD_ConvGamma_Mass_data,("after cuts"));
				leg1->Draw();
				DrawAliceLogoPerformance(right_up [0],right_up[1],right_up[2],right_up[3],0.03,Date);
		c_SinglePlot_15->Update();	
		c_SinglePlot_15->SaveAs(Form("%s%s/Leftover.%s",path,suffix,suffix));
		delete c_SinglePlot_15;







// ------------------------ True Photons for MonteCarlo Data ---------------------------------
	if(MCfile != ""){
	if(MC2){	
		
		Double_t truePhotons = TRUE_Conversion_R_montecarlo->Integral();
		Double_t measPhotons = ESD_Conversion_R_montecarlo->Integral();
		cout << "True Photons     " << truePhotons << endl;
		cout << "Simulated Photons     " << measPhotons << endl;

		TCanvas * c20 = new TCanvas("c3_1","",10,10,500,500);  // gives the page size
		c20->SetLogy(1);
	
		c20->cd();
		DrawAutoGammaHistos( ESD_Conversion_R_montecarlo,
						 TRUE_Conversion_R_montecarlo,  
						 "Comparison True Gammas and MC Gammas","R [cm]",StandardYAxis,
						 kFALSE,2 ,0.0000002,
						 kFALSE,0. ,0.,
						 kTRUE, 0.,180.);	
	
		c20->Update();
		c20->SaveAs("TRUE-SIM2.gif");
		delete c20;
	}
	}
	
//----------------------------- Comparison of Cuts ------------------------------------------------
	// Rescaling plots
		ESD_CutLikeSign_data->Scale(1./normFac_data);
		ESD_CutRefit_data->Scale(1./normFac_data);
		ESD_CutKink_data->Scale(1./normFac_data);
		ESD_CutdEdxSigmaElectronLine_data->Scale(1./normFac_data);
		ESD_CutdEdxSigmaPionLine_data->Scale(1./normFac_data);
		ESD_CutKaonRejectionLowP_data->Scale(1./normFac_data);
		ESD_CutProtonRejectionLowP_data->Scale(1./normFac_data);
		ESD_CutPionRejectionLowP_data->Scale(1./normFac_data);
		ESD_CutGetOnFly_data->Scale(1./normFac_data);
		ESD_CutNContributors_data->Scale(1./normFac_data);	
		ESD_CutPIDProb_data->Scale(1./normFac_data);
		ESD_CutR_data->Scale(1./normFac_data);
		ESD_CutLine_data->Scale(1./normFac_data);
		ESD_CutZ_data->Scale(1./normFac_data);
		ESD_CutSinglePt_data->Scale(1./normFac_data);
		ESD_CutMinNClsTPC_data->Scale(1./normFac_data);
		ESD_CutNDF_data->Scale(1./normFac_data);
		ESD_CutChi2_data->Scale(1./normFac_data); 
		ESD_CutEta_data->Scale(1./normFac_data);
		ESD_CutPt_data->Scale(1./normFac_data);
		ESD_ConvGamma_Mass_data->Scale(1./normFac_data);
		NoCutPt_data->Scale(1./normFac_data);
		NoCutEta_data->Scale(1./normFac_data);
		NoCutChi2_data->Scale(1./normFac_data);
		NoCutNDF_data->Scale(1./normFac_data);
		NoCutZ_data->Scale(1./normFac_data);	
		NoCutSinglePt_data->Scale(1./normFac_data);
		NoCutMinNClsTPC_data->Scale(1./normFac_data);
		NoCutLine_data->Scale(1./normFac_data);
		NoCutR_data->Scale(1./normFac_data);;
		NoCutPIDProb_data->Scale(1./normFac_data);
		NoCutNContributors_data->Scale(1./normFac_data);
		NoCutGetOnFly_data->Scale(1./normFac_data);
		NoCutKaonRejectionLowP_data->Scale(1./normFac_data);
		NoCutProtonRejectionLowP_data->Scale(1./normFac_data);
		NoCutPionRejectionLowP_data->Scale(1./normFac_data);
		NoCutdEdxSigmaPionLine_data->Scale(1./normFac_data);
		NoCutdEdxSigmaElectronLine_data->Scale(1./normFac_data);
		NoCutKink_data->Scale(1./normFac_data);
		NoCutRefit_data->Scale(1./normFac_data);
		NoCutLikeSign_data->Scale(1./normFac_data);
		ESD_AllV0_data->Scale(1./normFac_data);
		Final_data->Scale(1./normFac_data);
		if(MCfile != ""){
			ESD_CutLikeSign_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutRefit_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutKink_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutdEdxSigmaElectronLine_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutdEdxSigmaPionLine_montecarlo->Scale(1./normFac_montecarlo);
		ESD_CutKaonRejectionLowP_montecarlo->Scale(1./normFac_montecarlo);
		ESD_CutProtonRejectionLowP_montecarlo->Scale(1./normFac_montecarlo);
		ESD_CutPionRejectionLowP_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutGetOnFly_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutNContributors_montecarlo->Scale(1./normFac_montecarlo);	
			ESD_CutPIDProb_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutR_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutLine_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutZ_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutSinglePt_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutMinNClsTPC_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutNDF_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutChi2_montecarlo->Scale(1./normFac_montecarlo); 
			ESD_CutEta_montecarlo->Scale(1./normFac_montecarlo);
			ESD_CutPt_montecarlo->Scale(1./normFac_montecarlo);
			ESD_ConvGamma_Mass_montecarlo->Scale(1./normFac_montecarlo);
			NoCutPt_montecarlo->Scale(1./normFac_montecarlo);
			NoCutEta_montecarlo->Scale(1./normFac_montecarlo);
			NoCutChi2_montecarlo->Scale(1./normFac_montecarlo);
			NoCutNDF_montecarlo->Scale(1./normFac_montecarlo);
			NoCutZ_montecarlo->Scale(1./normFac_montecarlo);	
			NoCutSinglePt_montecarlo->Scale(1./normFac_montecarlo);
			NoCutMinNClsTPC_montecarlo->Scale(1./normFac_montecarlo);
			NoCutLine_montecarlo->Scale(1./normFac_montecarlo);
			NoCutR_montecarlo->Scale(1./normFac_montecarlo);;
			NoCutPIDProb_montecarlo->Scale(1./normFac_montecarlo);
			NoCutNContributors_montecarlo->Scale(1./normFac_montecarlo);
			NoCutGetOnFly_montecarlo->Scale(1./normFac_montecarlo);
		NoCutKaonRejectionLowP_montecarlo->Scale(1./normFac_montecarlo);
		NoCutProtonRejectionLowP_montecarlo->Scale(1./normFac_montecarlo);
		NoCutPionRejectionLowP_montecarlo->Scale(1./normFac_montecarlo);
			NoCutdEdxSigmaPionLine_montecarlo->Scale(1./normFac_montecarlo);
			NoCutdEdxSigmaElectronLine_montecarlo->Scale(1./normFac_montecarlo);
			NoCutKink_montecarlo->Scale(1./normFac_montecarlo);
			NoCutRefit_montecarlo->Scale(1./normFac_montecarlo);
			NoCutLikeSign_montecarlo->Scale(1./normFac_montecarlo);
			ESD_AllV0_montecarlo->Scale(1./normFac_montecarlo);
			Final_montecarlo->Scale(1./normFac_montecarlo);
		}
		TCanvas * c_SinglePlot_16 = new TCanvas("c_SinglePlot_16","",1000,1000);  // gives the page size
		c_SinglePlot_16->SetBottomMargin(1.5);
		c_SinglePlot_16->SetLogy(1);
		c_SinglePlot_16->cd();
		const Int_t nx = 20;
		char *cuts[nx] = {"V0Finder","LikeSign", "Refit", "Kink", "dEdxElectronLine", "dEdxPionLine", "KaonRejection", "Proton Rejection", "Pion Rejection" , "NContributors", "PIDProb", "R", "Line", "Z", "Min", "SinglePt","NDF", "Chi2", "Eta", "Pt"};
		leg2 = new TLegend(0.6,0.82,0.92,0.9);
			leg2->AddEntry(ESD_CutLikeSign_data,("Data"),"l");
			if(MCfile != ""){
			leg2->AddEntry(ESD_CutLikeSign_montecarlo,("MC"),"l");}
			leg2->SetFillColor(0);


		Double_t maxRange = ESD_CutdEdxSigmaElectronLine_data->Integral()/NumberV0_data * 100;
		if(MCfile != ""){
			if(maxRange < ESD_CutdEdxSigmaElectronLine_montecarlo->Integral()/NumberV0_montecarlo *100){
				maxRange = ESD_CutdEdxSigmaElectronLine_montecarlo->Integral()/NumberV0_montecarlo *100;
			}
		}		

		Output_data = new TFile("Output_data_2.root","RECREATE");
		
		TH1F *ComCuts_data = new TH1F("ComCuts_data","Comparison of Cuts",3,0,3);
		ComCuts_data->SetStats(0);
		ComCuts_data->SetYTitle("ratio Cut/V0");
		ComCuts_data->SetBit(TH1::kCanRebin);
		ComCuts_data->GetYaxis()->SetRangeUser(0.001, 5 *maxRange);
		ComCuts_data->GetYaxis()->SetTitleOffset(1.35);
		ComCuts_data->Fill(cuts[0],ESD_CutGetOnFly_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[1],ESD_CutLikeSign_data->Integral()/NumberV0_data *100);
		ComCuts_data->Fill(cuts[2],ESD_CutRefit_data->Integral()/NumberV0_data *100);
		ComCuts_data->Fill(cuts[3],ESD_CutKink_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[4],ESD_CutdEdxSigmaElectronLine_data->Integral()/NumberV0_data *100);
		ComCuts_data->Fill(cuts[5],ESD_CutdEdxSigmaPionLine_data->Integral()/NumberV0_data *100);
		ComCuts_data->Fill(cuts[6],ESD_CutKaonRejectionLowP_data->Integral()/NumberV0_data *100);
		ComCuts_data->Fill(cuts[7],ESD_CutProtonRejectionLowP_data->Integral()/NumberV0_data *100);
		ComCuts_data->Fill(cuts[8],ESD_CutPionRejectionLowP_data->Integral()/NumberV0_data *100);
		ComCuts_data->Fill(cuts[9],ESD_CutNContributors_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[10],ESD_CutPIDProb_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[11],ESD_CutR_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[12],ESD_CutLine_data->Integral()/NumberV0_data *100);		 		
		ComCuts_data->Fill(cuts[13],ESD_CutZ_data->Integral()/NumberV0_data *100);		 			 
		ComCuts_data->Fill(cuts[14],ESD_CutMinNClsTPC_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[15],ESD_CutSinglePt_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[16],ESD_CutNDF_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[17],ESD_CutChi2_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[18],ESD_CutEta_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->Fill(cuts[19],ESD_CutPt_data->Integral()/NumberV0_data *100);		 
		ComCuts_data->LabelsDeflate();
		ComCuts_data->Draw();
		ComCuts_data->Write();

		Output_data->Write();

		if(MCfile != ""){
			Output_montecarlo = new TFile("Output_montecarlo_2.root","RECREATE");
			TH1F *ComCuts_montecarlo = new TH1F("ComCuts_montecarlo","Comparison of Cuts",3,0,3);
			ComCuts_montecarlo->SetStats(0);
			ComCuts_montecarlo->SetBit(TH1::kCanRebin);
			ComCuts_montecarlo->SetLineColor(2);
			ComCuts_montecarlo->Fill(cuts[0],ESD_CutGetOnFly_montecarlo->Integral()/NumberV0_montecarlo *100);		
			ComCuts_montecarlo->Fill(cuts[1],ESD_CutLikeSign_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts_montecarlo->Fill(cuts[2],ESD_CutRefit_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts_montecarlo->Fill(cuts[3],ESD_CutKink_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[4],ESD_CutdEdxSigmaElectronLine_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts_montecarlo->Fill(cuts[5],ESD_CutdEdxSigmaPionLine_montecarlo->Integral()/NumberV0_montecarlo *100); 
			ComCuts_montecarlo->Fill(cuts[5],ESD_CutdEdxSigmaPionLine_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts_montecarlo->Fill(cuts[6],ESD_CutKaonRejectionLowP_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts_montecarlo->Fill(cuts[7],ESD_CutProtonRejectionLowP_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts_montecarlo->Fill(cuts[8],ESD_CutPionRejectionLowP_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts_montecarlo->Fill(cuts[9],ESD_CutNContributors_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[10],ESD_CutPIDProb_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[11],ESD_CutR_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[12],ESD_CutLine_montecarlo->Integral()/NumberV0_montecarlo *100);		 		
			ComCuts_montecarlo->Fill(cuts[13],ESD_CutZ_montecarlo->Integral()/NumberV0_montecarlo *100);		 			 
			ComCuts_montecarlo->Fill(cuts[14],ESD_CutMinNClsTPC_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[15],ESD_CutSinglePt_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[16],ESD_CutNDF_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[17],ESD_CutChi2_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[18],ESD_CutEta_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->Fill(cuts[19],ESD_CutPt_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts_montecarlo->LabelsDeflate();
			ComCuts_montecarlo->Draw("e,hist,same");
			ComCuts_montecarlo->Write();
			Output_montecarlo->Write();	
		}
		leg2->Draw();
		DrawAliceText(right_up_text[0],right_up_text[1], right_up_text[3]);				
		c_SinglePlot_16->Update();	
		c_SinglePlot_16->SaveAs(Form("%s%s/ComparisonCuts2.%s",path,suffix,suffix));


//------------------------------ Remaining V0 candidates after cut
		TCanvas * c_SinglePlot_17 = new TCanvas("c_SinglePlot_17","",1000,1000);  // gives the page size
		c_SinglePlot_17->SetLogy(1);		
		c_SinglePlot_17->cd();
		const Int_t nx = 21;
		char *cuts2[nx] = {"No Cut", "V0Finder", "LikeSign", "Refit", "Kink", "dEdxElectronLine", "dEdxPionLine",  "KaonRejection", "Proton Rejection", "Pion Rejection", "NContributors", "PIDProb", "R", "Line", "Z","MinClsTPC", "SinglePt","NDF", "Chi2", "Eta", "Pt"};
		Double_t maxRange = NoCutLikeSign_data->Integral()/NumberV0_data *100;
		if(MCfile != ""){
			if(maxRange < NoCutLikeSign_montecarlo->Integral()/NumberV0_montecarlo *100){
				maxRange = NoCutLikeSign_montecarlo->Integral()/NumberV0_montecarlo *100;
			}
		}		
		
		Output_data = new TFile("Output_data_2.root","UPDATE");
		
		TH1F *ComCuts2_data = new TH1F("ComCuts2_data","Remaining V0 candidates after Cuts ",3,0,3);	
		ComCuts2_data->SetYTitle("ratio Cut/V0");
		ComCuts2_data->GetYaxis()->SetTitleSize(0.03);
		ComCuts2_data->SetTitle();
		ComCuts2_data->SetStats(0);
		ComCuts2_data->SetBit(TH1::kCanRebin);
		ComCuts2_data->GetYaxis()->SetRangeUser(1, 300.);
		ComCuts2_data->GetYaxis()->SetTitleOffset(1.35);
		ComCuts2_data->Fill(cuts2[0],NoCutGetOnFly_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[1],NoCutLikeSign_data->Integral()/NumberV0_data *100);
		ComCuts2_data->Fill(cuts2[2],NoCutRefit_data->Integral()/NumberV0_data *100);
		ComCuts2_data->Fill(cuts2[3],NoCutKink_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[4],NoCutdEdxSigmaElectronLine_data->Integral()/NumberV0_data *100);
		ComCuts2_data->Fill(cuts2[5],NoCutdEdxSigmaPionLine_data->Integral()/NumberV0_data *100);
		ComCuts2_data->Fill(cuts2[6],NoCutKaonRejectionLowP_data->Integral()/NumberV0_data *100);
		ComCuts2_data->Fill(cuts2[7],NoCutProtonRejectionLowP_data->Integral()/NumberV0_data *100);
		ComCuts2_data->Fill(cuts2[8],NoCutPionRejectionLowP_data->Integral()/NumberV0_data *100);
		ComCuts2_data->Fill(cuts2[9],NoCutNContributors_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[10],NoCutPIDProb_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[11],NoCutR_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[12],NoCutLine_data->Integral()/NumberV0_data *100);		 		
		ComCuts2_data->Fill(cuts2[13],NoCutZ_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[14],NoCutMinNClsTPC_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[15],NoCutSinglePt_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[16],NoCutNDF_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[17],NoCutChi2_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[18],NoCutEta_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[19],NoCutPt_data->Integral()/NumberV0_data *100);		 
		ComCuts2_data->Fill(cuts2[20],Final_data->Integral()/NumberV0_data *100);			
		ComCuts2_data->Fill("RealLeftOver",ESD_ConvGamma_Mass_data->Integral()/NumberV0_data *100);	
		if(MC2){ ComCuts2_data->Fill("TruePhotons",0);}
		if(MC1){ ComCuts2_data->Fill("TruePhotons",TRUE_Conversion_R_data->Integral()/NumberV0_data *100);}
		ComCuts2_data->LabelsDeflate();
		ComCuts2_data->Draw();
		ComCuts2_data->Write();

		Output_data->Write();



		if(MCfile != ""){
			Output_montecarlo = new TFile("Output_montecarlo_2.root","UPDATE");
			TH1F *ComCuts2_montecarlo = new TH1F("ComCuts2_montecarlo","Remaining V0 candidates after Cuts",3,0,3);
			ComCuts2_montecarlo->SetTitle();
			ComCuts2_montecarlo->SetStats(0);
			ComCuts2_montecarlo->SetBit(TH1::kCanRebin);
			ComCuts2_montecarlo->SetLineColor(2);
			ComCuts2_montecarlo->Fill(cuts2[0],NoCutGetOnFly_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[1],NoCutLikeSign_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts2_montecarlo->Fill(cuts2[2],NoCutRefit_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts2_montecarlo->Fill(cuts2[3],NoCutKink_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[4],NoCutdEdxSigmaElectronLine_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts2_montecarlo->Fill(cuts2[5],NoCutdEdxSigmaPionLine_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts2_montecarlo->Fill(cuts2[6],NoCutKaonRejectionLowP_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts2_montecarlo->Fill(cuts2[7],NoCutProtonRejectionLowP_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts2_montecarlo->Fill(cuts2[8],NoCutPionRejectionLowP_montecarlo->Integral()/NumberV0_montecarlo *100);
			ComCuts2_montecarlo->Fill(cuts2[9],NoCutNContributors_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[10],NoCutPIDProb_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[11],NoCutR_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[12],NoCutLine_montecarlo->Integral()/NumberV0_montecarlo *100);		 		
			ComCuts2_montecarlo->Fill(cuts2[13],NoCutZ_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[14],NoCutMinNClsTPC_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[15],NoCutSinglePt_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[16],NoCutNDF_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[17],NoCutChi2_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[18],NoCutEta_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[19],NoCutPt_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			ComCuts2_montecarlo->Fill(cuts2[20],Final_montecarlo->Integral()/NumberV0_montecarlo *100);			
			ComCuts2_montecarlo->Fill("RealLeftOver",ESD_ConvGamma_Mass_montecarlo->Integral()/NumberV0_montecarlo *100);		 
			if(MC2){ ComCuts2_montecarlo->Fill("TruePhotons",TRUE_Conversion_R_montecarlo->Integral()/NumberV0_montecarlo *100);}
			ComCuts2_montecarlo->LabelsDeflate();
			ComCuts2_montecarlo->Draw("e,hist,same");
			ComCuts2_montecarlo->Write();
			Output_montecarlo->Write();
			
		}		
		leg2->Draw();
		DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);	
		
		c_SinglePlot_17->Update();	
		c_SinglePlot_17->SaveAs(Form("%s%s/RemainingV0s2.%s",path,suffix,suffix));
	
	}
	// end end end
	
	
}

