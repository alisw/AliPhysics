/****************************************************************************************************************************
****** 	provided by Gamma Conversion Group, PWGGA, 														*****
******		Kathrin Koch, kkoch@physi.uni-heidelberg.de 													*****
******		Friederike Bock, friederike.bock@cern.ch													*****
*****************************************************************************************************************************
*** This macro can be used to display the Photon Characteristics of the conversion method in ALICE, it can be operated  *****
*** on the output of the GammaConversionTask. It can take 2 input files, the second one should be MC, if this is not 	*****
*** the case all histograms including MC need to be commented out otherwise the running will crash.					*****
****************************************************************************************************************************/

#include <Riostream.h>
#include "PlottingGammaConversionHistos.h"
#include "PlottingGammaConversionAdditional.h"

extern TRandom *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem *gSystem;

void  Photon_Characteristics_Events(const char *data = "myOutput", const char *MCfile = "", const char *cutsel = "",const char *path = "", const char *output = "Photon-Characteristics", const char *Plots = "kTRUE", const char *suffix = "gif"){	
	
	gROOT->Reset();	
	gROOT->SetStyle("Plain");

	Bool_t SinglePlots = kFALSE;
	if(Plots == "kTRUE") SinglePlots = kTRUE;
	

	StyleSettingsThesis();	
	//StyleSettings();

	set_plot_style();
	
	// textsize for legend
	Float_t ts=0.04;
	// textsize for label
	Float_t ls = 0.04;
	
	//Array defintion for printing Logo in right upper corner
	Float_t right_up[4]={0.7,0.63,0.15, 0.02};
	Float_t right_up2D[4]={0.65,0.8,0.11, 0.02};
	Float_t right_down[4]={0.7,0.23,0.15, 0.02};
	//Array defintion for printing Logo in left upper corner
	Float_t left_up[4]={0.2,0.8, 0.15, 0.02};

	// which maximum pt in the plots
	Float_t maxPt=10;
	Float_t rebin = 4;	

	Char_t filename_data[100] = (Form("%s%s",path,data));
	
	char *StandardYAxis = "#gamma/ event scaled by multiplicity";
	char *Date = "25th June 2010";

	TLatex *Scaling = new TLatex(0.6,0.5,"MC scaled to mean of Data"); // Bo: this was modified
	            Scaling->SetNDC();
	            Scaling->SetTextColor(2);
	            Scaling->SetTextFont(42);
	            Scaling->SetTextSize(0.03);
	            Scaling->SetLineWidth(2);	

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

	Float_t linewidth = 1;

/*****************************************************************************************************************************
			Function defintions for plotting dEdx lines MC with Bethe Bloch parametrisation
				* parameters taken from TPC tender 
******************************************************************************************************************************
******************************************************************************************************************************/

		//Functions taken from TPC - for different line only change (x/0.000511) to (x/mass[TeV])
		TF1 foElecMC("foElecMC", "[0]*([1]*TMath::Power(TMath::Sqrt(1 + (x/0.000511)*(x/0.000511))/(x/0.000511) , [3]) - 1 - TMath::Power(TMath::Sqrt(1 + (x/0.000511)*(x/0.000511))/(x/0.000511) , [3])*TMath::Log([2] + 1/TMath::Power((x/0.000511), [4])))",0.005,100);
		TF1 foPionMC("foPionMC", "[0]*([1]*TMath::Power(TMath::Sqrt(1 + (x/0.13957)*(x/0.13957))/(x/0.13957) , [3]) - 1 - TMath::Power(TMath::Sqrt(1 + (x/0.13957)*(x/0.13957))/(x/0.13957) , [3])*TMath::Log([2] + 1/TMath::Power((x/0.13957), [4])))",0.005,100);
		TF1 foProtonMC("foProtonMC", "[0]*([1]*TMath::Power(TMath::Sqrt(1 + (x/0.93827)*(x/0.93827))/(x/0.93827) , [3]) - 1 - TMath::Power(TMath::Sqrt(1 + (x/0.93827)*(x/0.93827))/(x/0.93827) , [3])*TMath::Log([2] + 1/TMath::Power((x/0.93827), [4])))",0.005,100);
          TF1 foKaonMC("foKaronMC", "[0]*([1]*TMath::Power(TMath::Sqrt(1 + (x/0.493677)*(x/0.493677))/(x/0.493677) , [3]) - 1 - TMath::Power(TMath::Sqrt(1 + (x/0.493677)*(x/0.493677))/(x/0.493677) , [3])*TMath::Log([2] + 1/TMath::Power((x/0.493677), [4])))",0.005,100);

		// mip might change, see updates in TPC tender
		//Float_t mip = 50.3 ;
		foElecMC.SetParameters(2.15898e+00,1.75295e+01,3.40030e-09,1.96178e+00,3.91720e+00);
		foElecMC.SetLineColor(1);
		foElecMC.SetLineWidth(linewidth);	
		foPionMC.SetParameters(2.15898e+00,1.75295e+01,3.40030e-09,1.96178e+00,3.91720e+00);	
		foPionMC.SetLineColor(2);
		foPionMC.SetLineWidth(linewidth);
		foProtonMC.SetParameters(2.15898e+00,1.75295e+01,3.40030e-09,1.96178e+00,3.91720e+00);
		foProtonMC.SetLineColor(kGreen+3);
		foProtonMC.SetLineWidth(linewidth);
	     foKaonMC.SetParameters(2.15898e+00,1.75295e+01,3.40030e-09,1.96178e+00,3.91720e+00);
	     foKaonMC.SetLineColor(4);
		foKaonMC.SetLineWidth(linewidth);


/*****************************************************************************************************************************
			Function defintions for plotting dEdx lines DATA with Bethe Bloch Parametrisation 
				* parameters taken from TPC tender
******************************************************************************************************************************
******************************************************************************************************************************/

	     TF1 foElec("foElec", "(1.74648e+01*TMath::Power(TMath::Sqrt(1 +(x/0.000511)*(x/0.000511))/(x/0.000511) , 3.06481e+00) + 1.92177e+01 -1.90939e+00*TMath::Power(TMath::Sqrt(1 +(x/0.000511)*(x/0.000511))/(x/0.000511) ,3.06481e+00)*TMath::Log(2.78081e-07 + 1/TMath::Power((x/0.000511),3.26689e+00)))/4.77199500359375293e+01*50.3",0.005,100);
		foElec.SetLineColor(1);
		foElec.SetLineWidth(linewidth);

	     TF1 foPion("foPion", "(1.74648e+01*TMath::Power(TMath::Sqrt(1 +(x/0.13957)*(x/0.13957))/(x/0.13957) , 3.06481e+00) + 1.92177e+01 -1.90939e+00*TMath::Power(TMath::Sqrt(1 +(x/0.13957)*(x/0.13957))/(x/0.13957) ,3.06481e+00)*TMath::Log(2.78081e-07 + 1/TMath::Power((x/0.13957),3.26689e+00)))/4.77199500359375293e+01*50.3",0.005,100);
		foPion.SetLineColor(2);
		foPion.SetLineWidth(linewidth);

	    	TF1 foProton("foProton", "(1.74648e+01*TMath::Power(TMath::Sqrt(1 +(x/0.93827)*(x/0.93827))/(x/0.93827) , 3.06481e+00) + 1.92177e+01 -1.90939e+00*TMath::Power(TMath::Sqrt(1 +(x/0.93827)*(x/0.93827))/(x/0.93827) ,3.06481e+00)*TMath::Log(2.78081e-07 + 1/TMath::Power((x/0.93827),3.26689e+00)))/4.77199500359375293e+01*50.3",0.005,100);
		foProton.SetLineColor(kGreen+3);
		foProton.SetLineWidth(linewidth);

		TF1 foKaon("foKaron", "(1.74648e+01*TMath::Power(TMath::Sqrt(1 +(x/0.493677)*(x/0.493677))/(x/0.493677) , 3.06481e+00) + 1.92177e+01 -1.90939e+00*TMath::Power(TMath::Sqrt(1 +(x/0.493677)*(x/0.493677))/(x/0.493677) ,3.06481e+00)*TMath::Log(2.78081e-07 + 1/TMath::Power((x/0.493677),3.26689e+00)))/4.77199500359375293e+01*50.3",0.005,100);
		foKaon.SetLineColor(4);		
		foKaon.SetLineWidth(linewidth);

	
// ---------------------------- LOAD ALL FILES ---------------------------------------------------
	TFile f(filename_data);  
	
	//for new versions
	 TDirectory *fPWGGAGammaConversion_data = new TDirectory(); // definition of first folder / list	
	 TList *fHistosGammaConversion_data = new TList(); // definition of first folder / list
	 TList *fESDContainer_data = new TList();  // definition of following folder / list
	TList *fMappingContainer_data = new TList();
	 
	 if(!(fPWGGAGammaConversion_data = (TDirectory*)f.Get(GammaDirectory))) cout <<"PWGGAGammConversion TList NOT loaded correctly"<<endl; 
	 if(!(fHistosGammaConversion_data = (TList*)fPWGGAGammaConversion_data->Get(GammaList))) cout<<"histogramsAliGammaConversion NOT loaded correctly!"<<endl; 
	 if(!(fESDContainer_data = (TList*)fHistosGammaConversion_data->FindObject("ESD histograms"))) cout<<"ESD histograms NOT loaded correctly!"<<endl; 
	fMappingContainer_data = (TList*)fHistosGammaConversion_data->FindObject("Mapping histograms"); 

	TH2F * ESD_AllV0CurrentFinder_alfa_qt_data = fESDContainer_data->FindObject("ESD_AllV0sCurrentFinder_alfa_qt");	
	TH2F * ESD_ConvGamma_alfa_qt_data = fESDContainer_data->FindObject("ESD_ConvGamma_alfa_qt");	
	TH2F * ESD_ConvGamma_P_dEdxP_data = fESDContainer_data->FindObject("ESD_ConvGamma_P_dEdxP");
	TH2F * ESD_ConvGamma_E_dEdxP_data = fESDContainer_data->FindObject("ESD_ConvGamma_E_dEdxP");
	TH2F * ESD_ConvGamma_P_AsymmetryP_data = fESDContainer_data->FindObject("ESD_ConvGamma_P_AsymmetryP");
	TH1F * ESD_ConvGamma_DcaDaughters_data = fESDContainer_data->FindObject("ESD_ConvGamma_DcaDaughters");
	TH1F * ESD_ConvGamma_CosPointingAngle_data = fESDContainer_data->FindObject("ESD_ConvGamma_CosPointingAngle");
	TH1F * ESD_ConvGamma_NormDcaDistDaughters_data = fESDContainer_data->FindObject("ESD_ConvGamma_NormDcaDistDaughters");
	TH1F * ESD_ConvGamma_LikeLihood_data = fESDContainer_data->FindObject("ESD_ConvGamma_LikelihoodAP");
	TH1F * ESD_ConvGamma_Pt_data=fESDContainer_data->FindObject("ESD_ConvGamma_Pt");
	TH1F * ESD_ConvGamma_Eta_data=fESDContainer_data->FindObject("ESD_ConvGamma_Eta");
	TH1F * ESD_ConvGamma_Mass_data=fESDContainer_data->FindObject("ESD_ConvGamma_Mass");
	TH1F * ESD_ConvGamma_Chi2_data=fESDContainer_data->FindObject("ESD_ConvGamma_Chi2");
	TH1F * ESD_Mother_InvMass_data=fESDContainer_data->FindObject("ESD_Mother_InvMass");
	TH1F * ESD_NumberOfGoodESDTracks_data=fESDContainer_data->FindObject("ESD_NumberOfGoodESDTracksVtx");
	TH1F * ESD_NumberOfGoodESDTrackswith0_data=fESDContainer_data->FindObject("ESD_NumberOfGoodESDTracks");
	TH1F * ESD_E_Pt_data=fESDContainer_data->FindObject("ESD_E_Pt");
	TH1F * ESD_P_Pt_data=fESDContainer_data->FindObject("ESD_P_Pt");
	TH1F * ESD_E_nITSClusters_data=fESDContainer_data->FindObject("ESD_E_nITSClusters");
	TH1F * ESD_E_nTPCClusters_data=fESDContainer_data->FindObject("ESD_E_nTPCClusters");
	TH1F * ESD_P_nITSClusters_data=fESDContainer_data->FindObject("ESD_P_nITSClusters");
	TH1F * ESD_P_nTPCClusters_data=fESDContainer_data->FindObject("ESD_P_nTPCClusters");	
	TH1D *ESD_ConvGamma_P_dEdx_Y1_data=ESD_ConvGamma_P_dEdxP_data->ProjectionY("ESD_ConvGamma_dEdx_Y1_data",10,20);
	TH1D *ESD_ConvGamma_P_dEdx_Y2_data=ESD_ConvGamma_P_dEdxP_data->ProjectionY("ESD_ConvGamma_dEdx_Y2_data",20,30);
	TH1D *ESD_ConvGamma_P_dEdx_Y3_data=ESD_ConvGamma_P_dEdxP_data->ProjectionY("ESD_ConvGamma_dEdx_Y3_data",30,50);
	TH1D *ESD_ConvGamma_P_dEdx_Y4_data=ESD_ConvGamma_P_dEdxP_data->ProjectionY("ESD_ConvGamma_dEdx_Y4_data",50,100);
	TH1D *ESD_ConvGamma_P_AsymmetryP_Y1_data=ESD_ConvGamma_P_AsymmetryP_data->ProjectionY("ESD_ConvGamma_P_AsymmetryP_Y1_data",2,2);	
	TH1D *ESD_ConvGamma_P_AsymmetryP_Y2_data=ESD_ConvGamma_P_AsymmetryP_data->ProjectionY("ESD_ConvGamma_P_AsymmetryP_Y2_data",5,19);	
	TH1F * ESD_NumberOfContributorsVtx_data=fESDContainer_data->FindObject("ESD_NumberOfContributorsVtx");
	TH1F *ScalingDiagramm_data= fMappingContainer_data->FindObject("ESD_Conversion_Mapping_Phi_in_R_11");

	ESD_NumberOfContributorsVtx_data->SetAxisRange(1.,100.);
	ESD_NumberOfGoodESDTrackswith0_data->SetAxisRange(1.,100.);

	Float_t Scaling_data = ScalingDiagramm_data->Integral();
	Float_t nGoodEvents_data = ESD_NumberOfContributorsVtx_data->Integral();
	Float_t nGoodTrig_data = ESD_NumberOfContributorsVtx_data->GetEntries();
	cout<< data << "    Number of events::   " << nGoodEvents_data << "    Number of triggers::   " << nGoodTrig_data << endl;

	Float_t normFac_data=1./nGoodEvents_data;
//	Float_t normFac_data=1./Scaling_data;	//Float_t normFac_data=1.;
	Double_t mean_data = ESD_NumberOfGoodESDTracks_data->GetMean();

	//Scaling reconstr.

	GammaScalingHistogramm(ESD_ConvGamma_P_dEdxP_data,normFac_data);
	GammaScalingHistogramm(ESD_ConvGamma_E_dEdxP_data,normFac_data);
	GammaScalingHistogramm(ESD_ConvGamma_P_AsymmetryP_data,normFac_data);
	GammaScalingHistogramm(ESD_ConvGamma_Pt_data,normFac_data);
	GammaScalingHistogramm(ESD_ConvGamma_Eta_data,normFac_data);
	GammaScalingHistogramm(ESD_ConvGamma_Mass_data,normFac_data);
	GammaScalingHistogramm(ESD_ConvGamma_Chi2_data,normFac_data);
	GammaScalingHistogramm(ESD_E_Pt_data,normFac_data);
	GammaScalingHistogramm(ESD_P_Pt_data,normFac_data);
	GammaScalingHistogramm(ESD_E_nITSClusters_data,normFac_data);
	GammaScalingHistogramm(ESD_P_nITSClusters_data,normFac_data);
	GammaScalingHistogramm(ESD_E_nTPCClusters_data,normFac_data);
	GammaScalingHistogramm(ESD_P_nTPCClusters_data,normFac_data);
	GammaScalingHistogramm(ESD_Mother_InvMass_data,normFac_data);
	GammaScalingHistogramm(ESD_ConvGamma_DcaDaughters_data,normFac_data);
	GammaScalingHistogramm(ESD_ConvGamma_CosPointingAngle_data,normFac_data);
	GammaScalingHistogramm(ESD_ConvGamma_NormDcaDistDaughters_data,normFac_data);
	GammaScalingHistogramm(ESD_ConvGamma_LikeLihood_data,normFac_data); 
	GammaScalingHistogramm(ESD_ConvGamma_P_dEdx_Y1_data,normFac_data); 
	GammaScalingHistogramm(ESD_ConvGamma_P_dEdx_Y2_data,normFac_data);
	GammaScalingHistogramm(ESD_ConvGamma_P_dEdx_Y3_data,normFac_data); 
	GammaScalingHistogramm(ESD_ConvGamma_P_dEdx_Y4_data,normFac_data);
	GammaScalingHistogramm(ESD_ConvGamma_P_AsymmetryP_Y1_data,normFac_data); 
	GammaScalingHistogramm(ESD_ConvGamma_P_AsymmetryP_Y2_data,normFac_data); 
	GammaScalingHistogramm(ESD_NumberOfGoodESDTracks_data,normFac_data); 
	GammaScalingHistogramm(ESD_NumberOfGoodESDTrackswith0_data,normFac_data); 
	ESD_Mother_InvMass_data->Rebin(rebin);
	MultoutData = new TFile("Multout1.root","RECREATE");
	ESD_NumberOfGoodESDTracks_data->Write();
	MultoutData->Write();
	MultoutData->Close();
	//------------------------------- SECOND FILE ------------------------------
	TFile *montecarlo = 0x0;
	
	
	if(MCfile != ""){
		montecarlo = new TFile(Form("%s%s",path, MCfile));
		
		
		// for new versions
		 TDirectory *fPWGGAGammaConversion_montecarlo = new TDirectory(); // definition of first folder / list	
		 TList *fHistosGammaConversion_montecarlo = new TList(); // definition of first folder / list
		 TList *fESDContainer_montecarlo = new TList();  // definition of following folder / list
		 TList *fMappingContainer_montecarlo = new TList();  // definition of following folder / list
		 TList *fMCtruth_montecarlo = new TList();
		 if(!(fPWGGAGammaConversion_montecarlo = (TDirectory*)montecarlo->Get(GammaDirectory))) cout <<"PWGGAGammConversion TList NOT loaded correctly"<<endl; 
		 if(!(fHistosGammaConversion_montecarlo = (TList*)fPWGGAGammaConversion_montecarlo->Get(GammaList))) cout<<"histogramsAliGammaConversion NOT loaded correctly!"<<endl; 
		 if(!(fESDContainer_montecarlo = (TList*)fHistosGammaConversion_montecarlo->FindObject("ESD histograms"))) cout<<"ESD histograms NOT loaded correctly!"<<endl; 
		fMappingContainer_montecarlo = (TList*)fHistosGammaConversion_montecarlo->FindObject("Mapping histograms");		 
		fMCtruth_montecarlo = (TList*)fHistosGammaConversion_montecarlo->FindObject("MC histograms");		 

		TH2F * ESD_AllV0CurrentFinder_alfa_qt_montecarlo = fESDContainer_montecarlo->FindObject("ESD_AllV0sCurrentFinder_alfa_qt");	
		TH2F * ESD_ConvGamma_alfa_qt_montecarlo = fESDContainer_montecarlo->FindObject("ESD_ConvGamma_alfa_qt");			
		TH2F * ESD_ConvGamma_P_dEdxP_montecarlo = fESDContainer_montecarlo->FindObject("ESD_ConvGamma_P_dEdxP");
		TH2F * ESD_ConvGamma_E_dEdxP_montecarlo = fESDContainer_montecarlo->FindObject("ESD_ConvGamma_E_dEdxP");
		TH2F * ESD_ConvGamma_P_AsymmetryP_montecarlo = fESDContainer_montecarlo->FindObject("ESD_ConvGamma_P_AsymmetryP");
		TH2F * MC_ConvGamma_P_AsymmetryP_montecarlo = fMCtruth_montecarlo->FindObject("MC_ConvGamma_P_AsymmetryP");
		TH1F * ESD_ConvGamma_DcaDaughters_montecarlo = fESDContainer_montecarlo->FindObject("ESD_ConvGamma_DcaDaughters");
		TH1F * ESD_ConvGamma_CosPointingAngle_montecarlo = fESDContainer_montecarlo->FindObject("ESD_ConvGamma_CosPointingAngle");
		TH1F * ESD_ConvGamma_NormDcaDistDaughters_montecarlo = fESDContainer_montecarlo->FindObject("ESD_ConvGamma_NormDcaDistDaughters");
		TH1F * ESD_ConvGamma_LikeLihood_montecarlo = fESDContainer_montecarlo->FindObject("ESD_ConvGamma_LikelihoodAP");
		TH1F * ESD_ConvGamma_Pt_montecarlo=fESDContainer_montecarlo->FindObject("ESD_ConvGamma_Pt");
		TH1F * ESD_ConvGamma_Eta_montecarlo=fESDContainer_montecarlo->FindObject("ESD_ConvGamma_Eta");
		TH1F * ESD_ConvGamma_Mass_montecarlo=fESDContainer_montecarlo->FindObject("ESD_ConvGamma_Mass");
		TH1F * ESD_ConvGamma_Chi2_montecarlo=fESDContainer_montecarlo->FindObject("ESD_ConvGamma_Chi2");
		TH1F * ESD_E_Pt_montecarlo=fESDContainer_montecarlo->FindObject("ESD_E_Pt");
		TH1F * ESD_P_Pt_montecarlo=fESDContainer_montecarlo->FindObject("ESD_P_Pt");
		TH1F * ESD_E_nITSClusters_montecarlo=fESDContainer_montecarlo->FindObject("ESD_E_nITSClusters");
		TH1F * ESD_E_nTPCClusters_montecarlo=fESDContainer_montecarlo->FindObject("ESD_E_nTPCClusters");
		TH1F * ESD_P_nITSClusters_montecarlo=fESDContainer_montecarlo->FindObject("ESD_P_nITSClusters");
		TH1F * ESD_P_nTPCClusters_montecarlo=fESDContainer_montecarlo->FindObject("ESD_P_nTPCClusters");	
		TH1F * ESD_Mother_InvMass_montecarlo=fESDContainer_montecarlo->FindObject("ESD_Mother_InvMass");
		TH1F * ESD_NumberOfGoodESDTracks_montecarlo=fESDContainer_montecarlo->FindObject("ESD_NumberOfGoodESDTracksVtx");
		TH1F * ESD_NumberOfGoodESDTrackswith0_montecarlo=fESDContainer_montecarlo->FindObject("ESD_NumberOfGoodESDTracks");		
		TH1D *ESD_ConvGamma_P_dEdx_Y1_montecarlo=ESD_ConvGamma_P_dEdxP_montecarlo->ProjectionY("ESD_ConvGamma_dEdx_Y1_montecarlo",1,10);
		TH1D *ESD_ConvGamma_P_dEdx_Y2_montecarlo=ESD_ConvGamma_P_dEdxP_montecarlo->ProjectionY("ESD_ConvGamma_dEdx_Y2_montecarlo",10,20);	
		TH1D *ESD_ConvGamma_P_dEdx_Y3_montecarlo=ESD_ConvGamma_P_dEdxP_montecarlo->ProjectionY("ESD_ConvGamma_dEdx_Y3_montecarlo",20,50);
		TH1D *ESD_ConvGamma_P_dEdx_Y4_montecarlo=ESD_ConvGamma_P_dEdxP_montecarlo->ProjectionY("ESD_ConvGamma_dEdx_Y4_montecarlo",50,100);
		TH1D *ESD_ConvGamma_P_AsymmetryP_Y1_montecarlo=ESD_ConvGamma_P_AsymmetryP_montecarlo->ProjectionY("ESD_ConvGamma_P_AsymmetryP_Y1_montecarlo",2,2);
		TH1D *ESD_ConvGamma_P_AsymmetryP_Y2_montecarlo=ESD_ConvGamma_P_AsymmetryP_montecarlo->ProjectionY("ESD_ConvGamma_P_AsymmetryP_Y2_montecarlo",5,19);
		TH1D *MC_ConvGamma_P_AsymmetryP_Y1_montecarlo=MC_ConvGamma_P_AsymmetryP_montecarlo->ProjectionY("MC_ConvGamma_P_AsymmetryP_Y1_montecarlo",2,2);
		TH1D *MC_ConvGamma_P_AsymmetryP_Y2_montecarlo=MC_ConvGamma_P_AsymmetryP_montecarlo->ProjectionY("MC_ConvGamma_P_AsymmetryP_Y2_montecarlo",5,19);

		TH1F * ESD_NumberOfContributorsVtx_montecarlo=fESDContainer_montecarlo->FindObject("ESD_NumberOfContributorsVtx");
		TH1F *ScalingDiagramm_montecarlo= fMappingContainer_montecarlo->FindObject("ESD_Conversion_Mapping_Phi_in_R_11");

		ESD_NumberOfContributorsVtx_montecarlo->SetAxisRange(1.,100.);
		ESD_NumberOfGoodESDTrackswith0_montecarlo->SetAxisRange(1.,100.);

		Float_t Scaling_montecarlo = ScalingDiagramm_montecarlo->Integral();
		Float_t nGoodEvents_montecarlo = ESD_NumberOfContributorsVtx_montecarlo->Integral();
		Float_t nGoodTrig_montecarlo = ESD_NumberOfContributorsVtx_montecarlo->GetEntries();
		cout<< MCfile << "    Number of events::   " << nGoodEvents_montecarlo << "    Number of triggers::   " << nGoodTrig_montecarlo << endl;

		Double_t mean_montecarlo = ESD_NumberOfGoodESDTracks_montecarlo->GetMean();
		Float_t normFac_montecarlo=1./nGoodEvents_montecarlo * mean_data/mean_montecarlo;
	//	Float_t normFac_montecarlo=1./Scaling_montecarlo;
		Float_t nGoodEvents_MCtruth = 
		Float_t normFac_MCtruth = 1./
		
		//Scaling reconstr.
	
		GammaScalingHistogramm(ESD_ConvGamma_P_dEdxP_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_ConvGamma_E_dEdxP_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_ConvGamma_P_AsymmetryP_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_ConvGamma_Pt_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_ConvGamma_Eta_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_ConvGamma_Mass_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_ConvGamma_Chi2_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_E_Pt_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_P_Pt_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_E_nITSClusters_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_P_nITSClusters_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_E_nTPCClusters_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_P_nTPCClusters_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_Mother_InvMass_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_ConvGamma_DcaDaughters_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_ConvGamma_CosPointingAngle_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_ConvGamma_NormDcaDistDaughters_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_ConvGamma_LikeLihood_montecarlo,normFac_montecarlo); 
		GammaScalingHistogramm(ESD_ConvGamma_P_dEdx_Y1_montecarlo,normFac_montecarlo); 
		GammaScalingHistogramm(ESD_ConvGamma_P_dEdx_Y2_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_ConvGamma_P_dEdx_Y3_montecarlo,normFac_montecarlo); 
		GammaScalingHistogramm(ESD_ConvGamma_P_dEdx_Y4_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_ConvGamma_P_AsymmetryP_Y1_montecarlo,normFac_montecarlo); 
		GammaScalingHistogramm(MC_ConvGamma_P_AsymmetryP_Y1_montecarlo,normFac_montecarlo); 
		GammaScalingHistogramm(MC_ConvGamma_P_AsymmetryP_Y2_montecarlo,normFac_montecarlo); 
		GammaScalingHistogramm(ESD_ConvGamma_P_AsymmetryP_Y2_montecarlo,normFac_montecarlo); 
		GammaScalingHistogramm(ESD_NumberOfGoodESDTracks_montecarlo,normFac_montecarlo);
		GammaScalingHistogramm(ESD_NumberOfGoodESDTrackswith0_montecarlo,normFac_montecarlo);
		ESD_Mother_InvMass_montecarlo->Rebin(rebin);

		MultoutMC = new TFile("Multout2.root","RECREATE");
		ESD_NumberOfGoodESDTracks_montecarlo->Write();
		MultoutMC->Write();
		MultoutMC->Close();

	}
		//-----------------------------------------------------------------------------------------
		TF1* myGaus = new TF1("myGaus","gaus",50.,70.);
 		myGaus->SetParameters(0.010,65.,20);
		myGaus->SetParLimits(0,1e-3, 100e-3);
		myGaus->SetParLimits(1,60,90);
		ESD_ConvGamma_P_dEdx_Y1_data->Fit("myGaus","RME");
		Float_t mean_P_dEdx_Y1_data = myGaus->GetParameter(1);
		delete myGaus;

		TF1* myGaus = new TF1("myGaus","gaus",50.,70.);
 		myGaus->SetParameters(0.010,65.,20);
		myGaus->SetParLimits(0,1e-3, 100e-3);
		myGaus->SetParLimits(1,60,90);
		ESD_ConvGamma_P_dEdx_Y2_data->Fit("myGaus","RME");
		Float_t mean_P_dEdx_Y2_data = myGaus->GetParameter(1);
		delete myGaus;

		TF1* myGaus = new TF1("myGaus","gaus",60.,75.);
 		myGaus->SetParameters(0.010,65.,20);
		myGaus->SetParLimits(0,1e-3, 100e-3);
		myGaus->SetParLimits(1,60,90);
		ESD_ConvGamma_P_dEdx_Y3_data->Fit("myGaus","RME");
		Float_t mean_P_dEdx_Y3_data = myGaus->GetParameter(1);
		delete myGaus;

		TF1* myGaus = new TF1("myGaus","gaus",60.,79.);
 		myGaus->SetParameters(0.010,65.,20);
		myGaus->SetParLimits(0,1e-3, 100e-3);
		myGaus->SetParLimits(1,60,90);
		ESD_ConvGamma_P_dEdx_Y4_data->Fit("myGaus","RME");
		Float_t mean_P_dEdx_Y4_data = myGaus->GetParameter(1);
		delete myGaus;


	if(MCfile != ""){
		TF1* myGaus2 = new TF1("myGaus2","gaus",60.,79.);
 		myGaus2->SetParameters(0.010,80.,20);
		myGaus2->SetParLimits(0,1e-3, 100e-3);
		myGaus2->SetParLimits(1,60,100);
		ESD_ConvGamma_P_dEdx_Y1_montecarlo->Fit("myGaus2","RME");
		Float_t mean_P_dEdx_Y1_montecarlo = myGaus2->GetParameter(1);
		delete myGaus2;	
		TH1F* ESD_ConvGamma_P_dEdx_Y1_montecarlo_2 = (TH1F*) ESD_ConvGamma_P_dEdx_Y1_montecarlo->Clone();	
		Float_t nBins = ESD_ConvGamma_P_dEdx_Y1_montecarlo_2->GetNbinsX();
		for ( int i=1; i<nBins;i++) {	ESD_ConvGamma_P_dEdx_Y1_montecarlo_2->SetBinContent(i,0); }
		for ( int i=1; i<nBins;i++) {
			Float_t binwidth = ESD_ConvGamma_P_dEdx_Y1_montecarlo->GetBinWidth(i);
			Float_t yvalue=ESD_ConvGamma_P_dEdx_Y1_montecarlo->GetBinContent(i);
			Float_t xvalue=ESD_ConvGamma_P_dEdx_Y1_montecarlo->GetBinCenter(i);
			Float_t newX= xvalue*mean_P_dEdx_Y1_data/mean_P_dEdx_Y1_montecarlo; 
			Int_t j = (newX+0.5*binwidth)/binwidth;			
			ESD_ConvGamma_P_dEdx_Y1_montecarlo_2->AddBinContent(j,yvalue);
		}

		TF1* myGaus2 = new TF1("myGaus2","gaus",60.,79.);
 		myGaus2->SetParameters(0.010,80.,20);
		myGaus2->SetParLimits(0,1e-3, 100e-3);
		myGaus2->SetParLimits(1,60,100);
		ESD_ConvGamma_P_dEdx_Y2_montecarlo->Fit("myGaus2","RME");
		Float_t mean_P_dEdx_Y2_montecarlo = myGaus2->GetParameter(1);
		delete myGaus2;
		TH1F* ESD_ConvGamma_P_dEdx_Y2_montecarlo_2 = (TH1F*) ESD_ConvGamma_P_dEdx_Y2_montecarlo->Clone();	
		Float_t nBins = ESD_ConvGamma_P_dEdx_Y2_montecarlo_2->GetNbinsX();
		for ( int i=1; i<nBins;i++) {	ESD_ConvGamma_P_dEdx_Y2_montecarlo_2->SetBinContent(i,0);}
		for ( int i=1; i<nBins;i++) {
			Float_t binwidth = ESD_ConvGamma_P_dEdx_Y2_montecarlo->GetBinWidth(i);
			Float_t yvalue=ESD_ConvGamma_P_dEdx_Y2_montecarlo->GetBinContent(i);
			Float_t xvalue=ESD_ConvGamma_P_dEdx_Y2_montecarlo->GetBinCenter(i);
			Float_t newX= xvalue*mean_P_dEdx_Y2_data/mean_P_dEdx_Y2_montecarlo; 
			Int_t j = (newX+0.5*binwidth)/binwidth;			
			ESD_ConvGamma_P_dEdx_Y2_montecarlo_2->AddBinContent(newX,yvalue);
		}


		TF1* myGaus2 = new TF1("myGaus2","gaus",60.,79.);
 		myGaus2->SetParameters(0.010,81.,20);
		myGaus2->SetParLimits(0,1e-3, 100e-3);
		myGaus2->SetParLimits(1,60,100);
		ESD_ConvGamma_P_dEdx_Y3_montecarlo->Fit("myGaus2","RME");		
		Float_t mean_P_dEdx_Y3_montecarlo = myGaus2->GetParameter(1);
		delete myGaus2;
		TH1F* ESD_ConvGamma_P_dEdx_Y3_montecarlo_2 = (TH1F*) ESD_ConvGamma_P_dEdx_Y3_montecarlo->Clone();	
		Float_t nBins = ESD_ConvGamma_P_dEdx_Y3_montecarlo_2->GetNbinsX();
		for ( int i=1; i<nBins;i++) {	ESD_ConvGamma_P_dEdx_Y3_montecarlo_2->SetBinContent(i,0); }
		for ( int i=1; i<nBins;i++) {
			Float_t binwidth = ESD_ConvGamma_P_dEdx_Y3_montecarlo->GetBinWidth(i);
			Float_t yvalue=ESD_ConvGamma_P_dEdx_Y3_montecarlo->GetBinContent(i);
			Float_t xvalue=ESD_ConvGamma_P_dEdx_Y3_montecarlo->GetBinCenter(i);
			Float_t newX= xvalue*mean_P_dEdx_Y3_data/mean_P_dEdx_Y3_montecarlo; 
			Int_t j = (newX+0.5*binwidth)/binwidth;			
			ESD_ConvGamma_P_dEdx_Y3_montecarlo_2->AddBinContent(j,yvalue);
		}

		TF1* myGaus2 = new TF1("myGaus2","gaus",60.,79.);
		myGaus2->SetParameters(0.010,81.,20);
		myGaus2->SetParLimits(0,1e-4, 100e-3);
		myGaus2->SetParLimits(1,60,100);
		ESD_ConvGamma_P_dEdx_Y4_montecarlo->Fit("myGaus2","RME");
		Float_t mean_P_dEdx_Y4_montecarlo = myGaus2->GetParameter(1);
		delete myGaus2;
		TH1F* ESD_ConvGamma_P_dEdx_Y4_montecarlo_2 = (TH1F*) ESD_ConvGamma_P_dEdx_Y4_montecarlo->Clone();	
		Float_t nBins = ESD_ConvGamma_P_dEdx_Y4_montecarlo_2->GetNbinsX();
		for ( int i=1; i<nBins;i++) {	ESD_ConvGamma_P_dEdx_Y4_montecarlo_2->SetBinContent(i,0); }
		for ( int i=1; i<nBins;i++) {
			Float_t binwidth = ESD_ConvGamma_P_dEdx_Y4_montecarlo->GetBinWidth(i);
			Float_t yvalue=ESD_ConvGamma_P_dEdx_Y4_montecarlo->GetBinContent(i);
			Float_t xvalue=ESD_ConvGamma_P_dEdx_Y4_montecarlo->GetBinCenter(i);
			Float_t newX= xvalue*mean_P_dEdx_Y4_data/mean_P_dEdx_Y4_montecarlo; 
			Int_t j = (newX+0.5*binwidth)/binwidth;			
			ESD_ConvGamma_P_dEdx_Y4_montecarlo_2->AddBinContent(j,yvalue);
		}
	}

	
	// ----------------------------------------------------------------------------------------
	
	// ----------------------- Start  ps file -----------------------------------------------------
	
	TPostScript *ps_characteristics;
	if(!SinglePlots)ps_characteristics = new TPostScript(Form("%s%s.ps",path,output),111);	
		
	
	// ----------------------- dE/dx of pion/electron  -----------------------------------------

	// --------------------------- page 1 -------------------------------------------------------
	if(!SinglePlots) {
ps_characteristics->NewPage();
	
		TCanvas * c0 = new TCanvas("c0","",10,10,700,1000);  // gives the page size
		c0->Divide(2,2);
	
		title0 = new TPaveLabel(0.05,0.92,0.95,0.96,(Form("data: %s",data)));	
		title0->SetFillColor(16);
		title0->SetTextSize(0.25);
		title0->Draw();
	
		if(MCfile != ""){
			title1 = new TPaveLabel(0.05,0.87,0.95,0.91,(Form("montecarlo: %s",MCfile)));	
			title1->SetFillColor(16);
			title1->SetTextSize(0.25);
			title1->Draw();
		}
	

		c0->Update();
				
	//-------------------------- page 2 - dE/dx distribution-----------------------------------------------------------
		
		ps_characteristics->NewPage();

		TCanvas * c1 = new TCanvas("c1","",10,10,700,1000);  // gives the page size
		pad_c1 = new TPad("pad_c1","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c1->SetFillColor(0);
		pad_c1->GetFrame()->SetFillColor(0);
		pad_c1->SetBorderMode(0);
		pad_c1->Divide(2,3);
		pad_c1->Draw();
		
		title0->Draw();
		if(MCfile != ""){
			title1->Draw();
		}
				
		pad_c1->cd(1)->SetLogx(1);
		pad_c1->cd(1)->SetLogz(1);
		pad_c1->cd(2)->SetLogx(1);
		pad_c1->cd(2)->SetLogz(1);
		pad_c1->cd(3)->SetLogx(1);
		pad_c1->cd(4)->SetLogx(1);
		pad_c1->cd(3)->SetLogz(1);
		pad_c1->cd(4)->SetLogz(1);


		pad_c1->cd(1)->SetRightMargin(0.17); 
		pad_c1->cd(1);
			DrawAutoGammaHisto2D(	ESD_ConvGamma_P_dEdxP_data,
								"", "p [GeV/c]", "dE/dx Positron", "Data",
								kTRUE, 10., 140.,
								kTRUE, 0., 20.);
		foElec.Draw("e,hist,same");
		foPion.Draw("e,hist,same");
		foProton.Draw("e,hist,same");
		foKaon.Draw("e,hist,same");

		pad_c1->cd(2)->SetRightMargin(0.17); 
		if(MCfile != ""){
			pad_c1->cd(2);	
				DrawAutoGammaHisto2D(	ESD_ConvGamma_P_dEdxP_montecarlo,
									"", "p [GeV/c]", "dE/dx Positron", "MC",
									kTRUE, 20., 150.,
								kTRUE, 0., 20.);
			foElecMC.Draw("e,hist,same");
			foPionMC.Draw("e,hist,same");
			foProtonMC.Draw("e,hist,same");
			foKaonMC.Draw("e,hist,same");

		}

		pad_c1->cd(3)->SetRightMargin(0.17); 
		pad_c1->cd(3);
			DrawAutoGammaHisto2D(	ESD_ConvGamma_E_dEdxP_data,
								"", "p [GeV/c]", "dE/dx Electron", "Data",
								kTRUE, 10., 140.,
								kTRUE, 0., 20.);
		foElec.Draw("e,hist,same");
		foPion.Draw("e,hist,same");
		foProton.Draw("e,hist,same");
		foKaon.Draw("e,hist,same");

		pad_c1->cd(4)->SetRightMargin(0.17); 
		if(MCfile != ""){
			pad_c1->cd(4);	
				DrawAutoGammaHisto2D(	ESD_ConvGamma_E_dEdxP_montecarlo,
									"", "p [GeV/c]", "dE/dx Electron", "MC",
									kTRUE, 20., 150.,
								kTRUE, 0., 20.);
			foElecMC.Draw("e,hist,same");
			foPionMC.Draw("e,hist,same");
			foProtonMC.Draw("e,hist,same");
			foKaonMC.Draw("e,hist,same");

		}

		pad_c1->cd(5);
			if(MCfile == ""){
				ESD_ConvGamma_P_dEdx_Y3_data->Rebin(2);
				DrawAutoGammaHisto( ESD_ConvGamma_P_dEdx_Y3_data, 
								 "Projection of dE/dx 20th-50th pt-bins","dE/dx ",StandardYAxis,
								 kTRUE, 1.2,0,
								 kFALSE,0. ,0.,
								 kTRUE, 50.,100.);
			}else{
				ESD_ConvGamma_P_dEdx_Y3_data->Rebin(2);
				ESD_ConvGamma_P_dEdx_Y3_montecarlo_2->Rebin(2);
				DrawAutoGammaHistos( ESD_ConvGamma_P_dEdx_Y3_data, 
								 ESD_ConvGamma_P_dEdx_Y3_montecarlo_2, 
								 "Projection of dE/dx 20th-50th pt-bins","dE/dx ",StandardYAxis,
								 kTRUE, 1.2,0,
								 kFALSE,0. ,0.,
								 kTRUE, 50.,100.);
			}
			Scaling->Draw();		

		pad_c1->cd(6);
			if(MCfile == ""){
				ESD_ConvGamma_P_dEdx_Y4_data->Rebin(2);
				DrawAutoGammaHisto( ESD_ConvGamma_P_dEdx_Y4_data, 
								 "Projection of dE/dx in other pt-bins","dE/dx ",StandardYAxis,
								 kTRUE, 1.2,0,
								 kFALSE,0. ,0.,
								 kTRUE, 50.,100.);
			}else{
				ESD_ConvGamma_P_dEdx_Y4_data->Rebin(2);
				ESD_ConvGamma_P_dEdx_Y4_montecarlo_2->Rebin(2);
				DrawAutoGammaHistos( ESD_ConvGamma_P_dEdx_Y4_data, 
								 ESD_ConvGamma_P_dEdx_Y4_montecarlo_2, 
								 "Projection of dE/dx in other pt-bins","dE/dx ",StandardYAxis,
								 kTRUE, 1.2,0,
								 kFALSE,0. ,0.,
								 kTRUE, 50.,100.);
			}
			Scaling->Draw();		
		
		c1->Update();

 //----------------------------- page 3 - Asymmetry distribution ---------------------------------------------------------------
		ps_characteristics->NewPage();
		
		TCanvas * c2 = new TCanvas("c2","",10,10,700,1000);  // gives the page size
		pad_c2 = new TPad("pad_c2","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c2->SetFillColor(0);
		pad_c2->GetFrame()->SetFillColor(0);
		pad_c2->SetBorderMode(0);
		pad_c2->Divide(2,2);
		pad_c2->Draw();
		
		title0->Draw();
		if(MCfile != ""){
			title1->Draw();
		}

		pad_c2->cd(1)->SetLogz(1);
		pad_c2->cd(2)->SetLogz(1);

		pad_c2->cd(1);
			pad_c2->cd(1)->SetRightMargin(0.17);
			DrawAutoGammaHisto2D(	ESD_ConvGamma_P_AsymmetryP_data,
								"", "p [GeV/c]", "Asymmetry", "Data",
								kFALSE, 30., 100.,
								kTRUE, 0., 20.);

		if(MCfile != ""){
			pad_c2->cd(2);
				pad_c2->cd(2)->SetRightMargin(0.17);
				DrawAutoGammaHisto2D(	ESD_ConvGamma_P_AsymmetryP_montecarlo,
									"", "p [GeV/c]", "Asymmetry", "MC",
									kFALSE, 30., 100.,
									kTRUE, 0., 20.);
		}

		pad_c2->cd(3);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_P_AsymmetryP_Y1_data, 
								"Projection of Asymmetry in 1st -3rd bin pt-bin", "Asymmetry",StandardYAxis,
								kTRUE, 1.2,0,
								kFALSE,0. ,0., 
								kFALSE, 0.,0.9);

			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_P_AsymmetryP_Y1_data,
								 ESD_ConvGamma_P_AsymmetryP_Y1_montecarlo,
								 "Projection of Asymmetry in 1st -3rd bin pt-bin", "Asymmetry",StandardYAxis, 
								 kTRUE, 1.2,0,
								 kFALSE,0. ,0.,
								 kFALSE, 0.1,0.9);
			}

		pad_c2->cd(4);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_P_AsymmetryP_Y2_data,
								"Projection of Asymmetry in other pt-bins", "Asymmetry",StandardYAxis, 
								kTRUE, 1.2,0,
								kFALSE,0. ,0., 
								kFALSE, 0.,1.);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_P_AsymmetryP_Y2_data, 
								 ESD_ConvGamma_P_AsymmetryP_Y2_montecarlo, 
								 "Projection of Asymmetry in other pt-bins", "Asymmetry",StandardYAxis, 
								 kTRUE, 1.2,0,
								 kFALSE,0. ,0., 
								 kFALSE, 0.,1.);
			}
		c2->Update();
		

 //----------------------------- page 4 - Armenteros-plots ---------------------------------------------------------------
		ps_characteristics->NewPage();
		
		TCanvas * c20 = new TCanvas("c20","",10,10,700,1000);  // gives the page size
		pad_c20 = new TPad("pad_c20","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c20->SetFillColor(0);
		pad_c20->GetFrame()->SetFillColor(0);
		pad_c20->SetBorderMode(0);
		pad_c20->Divide(2,2);
		pad_c20->Draw();
		
		pad_c20->cd(1)->SetLogz(1);
		pad_c20->cd(2)->SetLogz(1);
		pad_c20->cd(3)->SetLogz(1);
		pad_c20->cd(4)->SetLogz(1);

		title0->Draw();
		if(MCfile != ""){
			title1->Draw();
		}

		pad_c20->cd(1);
			pad_c20->cd(1)->SetRightMargin(0.17);
			DrawAutoGammaHisto2D(	ESD_AllV0CurrentFinder_alfa_qt_data,
								"", "#alpha = (p^{+}_{L}-p^{-}_{L})/(p^{+}_{L}+p^{-}_{L})", "q_{T} [GeV/c]", "Data",
								kFALSE, 30., 100.,
								kTRUE, -1., 1.);

		if(MCfile != ""){
			pad_c20->cd(2);
				pad_c20->cd(2)->SetRightMargin(0.17);
				DrawAutoGammaHisto2D(	ESD_AllV0CurrentFinder_alfa_qt_montecarlo,
									"","#alpha = (p^{+}_{L}-p^{-}_{L})/(p^{+}_{L}+p^{-}_{L})", "q_{T} [GeV/c]", "MC",
									kFALSE, 30., 100.,
									kTRUE, -1., 1.);
		}

		pad_c20->cd(3);
			pad_c20->cd(3)->SetRightMargin(0.17);
			DrawAutoGammaHisto2D(	ESD_ConvGamma_alfa_qt_data,
								"", "#alpha = (p^{+}_{L}-p^{-}_{L})/(p^{+}_{L}+p^{-}_{L})", "q_{T} [GeV/c]", "Data",
								kFALSE, 30., 100.,
								kTRUE, -1., 1.);

		if(MCfile != ""){
			pad_c20->cd(4);
				pad_c20->cd(4)->SetRightMargin(0.17);
				DrawAutoGammaHisto2D(	ESD_ConvGamma_alfa_qt_montecarlo,
									"", "#alpha = (p^{+}_{L}-p^{-}_{L})/(p^{+}_{L}+p^{-}_{L})", "q_{T} [GeV/c]", "MC",
									kFALSE, 30., 100.,
									kTRUE, -1., 1.);
		}


		c20->Update();
		
		// -----------------------------Page 5 ------------------------------------------------
		
		ps_characteristics->NewPage();
		
		TCanvas * c3 = new TCanvas("c3","",10,10,700,1000);  // gives the page size
		pad_c3 = new TPad("pad_c3","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c3->SetFillColor(0);
		pad_c3->GetFrame()->SetFillColor(0);
		pad_c3->SetBorderMode(0);
		pad_c3->Divide(2,2);
		pad_c3->Draw();
		
		title0->Draw();
		if(MCfile != ""){
			title1->Draw();
		}
		
		pad_c3->cd(1)->SetLogy(1);
		pad_c3->cd(2)->SetLogy(1);
		pad_c3->cd(3)->SetLogy(1);
		pad_c3->cd(4)->SetLogy(1);
		
		pad_c3->cd(1);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_DcaDaughters_data,
								"", "DCA [cm]",StandardYAxis,
								kTRUE, 3,0.000001,								kFALSE,0 ,0, 
								kFALSE, 0,1);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_DcaDaughters_data,
								 ESD_ConvGamma_DcaDaughters_montecarlo, 
								 "", "DCA [cm]",StandardYAxis,
								 kTRUE, 3,0.000001,
								 kFALSE,0. ,0., 
								 kFALSE, 0.,1);
			}
		
		pad_c3->cd(2);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_NormDcaDistDaughters_data,
								"", "Normalized distance [cm]",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0 ,0, 
								kFALSE, 0,1);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_NormDcaDistDaughters_data,
								 ESD_ConvGamma_NormDcaDistDaughters_montecarlo, 
								"", "Normalized distance [cm]",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0. ,0., 
								kFALSE, 0.,1);
			}
		
		pad_c3->cd(3);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_CosPointingAngle_data,
								"", "Cos Pointing Angle",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0 ,0, 
								kTRUE, 0,1);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_CosPointingAngle_data,
								 ESD_ConvGamma_CosPointingAngle_montecarlo, 
								"", "Cos Pointing Angle",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0. ,0., 
								kTRUE, 0.,1);
			}
				
		pad_c3->cd(4);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_LikeLihood_data,
								"", "Likelihood",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0 ,0, 
								kTRUE, 0,1);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_LikeLihood_data,
								 ESD_ConvGamma_LikeLihood_montecarlo, 
								"", "Likelihood",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0. ,0., 
								kTRUE, 0.,1);
			}

		
		c3->Update();
	//-----------------------------Page 6 - Clusters---------------------------------
		ps_characteristics->NewPage();
		
		TCanvas * c6 = new TCanvas("c6","",10,10,700,1000);  // gives the page size
		
		pad_c6 = new TPad("pad_c6","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c6->SetFillColor(0);
		pad_c6->GetFrame()->SetFillColor(0);
		pad_c6->SetBorderMode(0);
		pad_c6->Divide(2,3);
		pad_c6->Draw();
		
		title0->Draw();
		if(MCfile != ""){
			title1->Draw();
		}
		
		pad_c6->cd(1)->SetLogy(1);
		pad_c6->cd(2)->SetLogy(1);
		TText *l = new TText();		
		
		pad_c6->cd(1);	
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_NumberOfGoodESDTracks_data, 
								 "","Number good ESDTracks ","",
								 kTRUE, 3.,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,100.);
				Double_t mean_data = ESD_NumberOfGoodESDTracks_data->GetMean();
				l->SetTextAlign(11);
				l->SetTextColor(1);			
				l->SetTextSize(0.03);
				l->DrawText(60,0.0002, "Mean Data:");
				l->DrawText(80,0.0002, Form("%lg",  mean_data));	

			}else{
				DrawAutoGammaHistos( ESD_NumberOfGoodESDTracks_data, 
								 ESD_NumberOfGoodESDTracks_montecarlo, 
								 "","Number good ESDTracks","",
								 kTRUE, 3.,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,100.);
				Double_t mean_data = ESD_NumberOfGoodESDTracks_data->GetMean();
				l->SetTextAlign(11);
				l->SetTextColor(1);			
				l->SetTextSize(0.03);
				l->DrawText(60,0.0002, "Mean Data:");
				l->DrawText(80,0.0002, Form("%lg",  mean_data));	
	
				Double_t mean_montecarlo = ESD_NumberOfGoodESDTracks_montecarlo->GetMean();
				l->SetTextAlign(11);
				l->SetTextColor(2);			
				l->SetTextSize(0.03);
				l->DrawText(60,0.0001, "Mean MC:");
				l->DrawText(80,0.0001, Form("%lg",  mean_montecarlo));		
	
			}

		pad_c6->cd(2);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_NumberOfGoodESDTrackswith0_data, 
								 "","Number good ESDTracks with 0 ","",
								 kTRUE, 3.,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,100.);
				Double_t mean_data = ESD_NumberOfGoodESDTrackswith0_data->GetMean();
				l->SetTextAlign(11);
				l->SetTextColor(1);			
				l->SetTextSize(0.03);
				l->DrawText(60,0.0002, "Mean Data:");
				l->DrawText(80,0.0002, Form("%lg",  mean_data));	

			}else{
				DrawAutoGammaHistos( ESD_NumberOfGoodESDTrackswith0_data, 
								 ESD_NumberOfGoodESDTrackswith0_montecarlo, 
								 "","Number good ESDTracks with 0","",
								 kTRUE, 3.,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,100.);
				Double_t mean_data = ESD_NumberOfGoodESDTrackswith0_data->GetMean();
				l->SetTextAlign(11);
				l->SetTextColor(1);			
				l->SetTextSize(0.03);
				l->DrawText(60,0.0002, "Mean Data:");
				l->DrawText(80,0.0002, Form("%lg",  mean_data));	
	
				Double_t mean_montecarlo = ESD_NumberOfGoodESDTrackswith0_montecarlo->GetMean();
				l->SetTextAlign(11);
				l->SetTextColor(2);			
				l->SetTextSize(0.03);
				l->DrawText(60,0.0001, "Mean MC:");
				l->DrawText(80,0.0001, Form("%lg",  mean_montecarlo));		
	
			}

		pad_c6->cd(3);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_E_nITSClusters_data,
								"", "Number of ITS Clusters for Electron ","Hits/event scaled by multiplicity ",
								kTRUE, 1.5,0,
								kFALSE,0 ,0, 
								kFALSE, 0,0);
			}else{
				DrawAutoGammaHistos( ESD_E_nITSClusters_data,
								 ESD_E_nITSClusters_montecarlo, 
								"","Number of ITS Clusters for Electron ","Hits/event scaled by multiplicity", 
								kTRUE, 1.5,0,
								kFALSE,0. ,0., 
								kFALSE, 0.,0.);
			}
	
		pad_c6->cd(4);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_E_nTPCClusters_data,
								"", "Number of TPC Clusters for Electron","Hits/event scaled by multiplicity", 
								kTRUE, 1.5,0,
								kFALSE,0. ,0., 
								kFALSE, 0.,0.);
			}else{
				DrawAutoGammaHistos( ESD_E_nTPCClusters_data, 
								 ESD_E_nTPCClusters_montecarlo,
								 "","Number of TPC Clusters for Electron ","Hits/event scaled by multiplicity",
								 kTRUE, 1.5,0,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,0.);
			}

		pad_c6->cd(5);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_P_nITSClusters_data,
								"", "Number of ITS Clusters for PositronS","Hits/event scaled by multiplicity", 
								kTRUE, 1.5, 0,
								kFALSE,0. ,0., 
								kFALSE, 0.,0.);
			}else{
				DrawAutoGammaHistos( ESD_P_nITSClusters_data, 
								 ESD_P_nITSClusters_montecarlo,
								 "","Number of ITS Clusters for Positron","Hits/event scaled by multiplicity", 
								 kTRUE, 1.5, 0,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,0.);
			}

		pad_c6->cd(6);
			if(MCfile != ""){	
				DrawAutoGammaHistos( ESD_P_nTPCClusters_data, 
								 ESD_P_nTPCClusters_montecarlo, 
								 "","Number of TPC Clusters for Positron","Hits/event scaled by multiplicity",
								 kTRUE, 1.5, 0,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,0.);
				
			}else{
				DrawAutoGammaHisto( ESD_P_nTPCClusters_data,
								"", "Number of TPC Clusters for PositronS","Hits/event scaled by multiplicity", 
								kTRUE, 1.5, 0,
								kFALSE,0. ,0., 
								kFALSE, 0.,0.);			
			}

		c6->Update();
	
		// --------------------------Page 7-------------------------------------
		
		ps_characteristics->NewPage();
		
		TCanvas * c4 = new TCanvas("c4","",10,10,700,1000);  // gives the page size
		pad_c4 = new TPad("pad_c4","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c4->SetFillColor(0);
		pad_c4->GetFrame()->SetFillColor(0);
		pad_c4->SetBorderMode(0);
		pad_c4->Divide(2,3);
		pad_c4->Draw();
		
		title0->Draw();
		if(MCfile != ""){
			title1->Draw();
		}

		pad_c4->cd(1)->SetLogy(1);
		pad_c4->cd(2)->SetLogy(1);
		pad_c4->cd(3)->SetLogy(1);
		pad_c4->cd(4)->SetLogy(1);
		pad_c4->cd(5)->SetLogy(1);
		pad_c4->cd(6)->SetLogy(1);
		
		
		pad_c4->cd(1);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_Pt_data,
								"", "p_{t} Conv #gamma [GeV]",StandardYAxis,
								kTRUE, 3,0.00000001,
								kFALSE,0 ,0, 
								kTRUE, 0,20.);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_Pt_data,
								 ESD_ConvGamma_Pt_montecarlo, 
								"", "p_{t} Conv #gamma [GeV]",StandardYAxis,
								kTRUE, 3,0.00000001,
								kFALSE,0. ,0., 
								kTRUE, 0.,20.);
			}

		pad_c4->cd(2);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_Eta_data,
								"", "#eta Conv  #gamma",StandardYAxis,
								kTRUE, 3,0.001,
								kFALSE,0 ,0, 
								kFALSE, 0,0.7);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_Eta_data,
								 ESD_ConvGamma_Eta_montecarlo, 
								"", "#eta Conv  #gamma",StandardYAxis,
								kTRUE, 3,0.001,
								kFALSE,0. ,0., 
								kFALSE, 0.,0.7);
			}
		
		pad_c4->cd(3);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_Mass_data,
								"", "InvMass e^{+}e^{-} [GeV]",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0 ,0, 
								kTRUE, 0,0.1);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_Mass_data,
								 ESD_ConvGamma_Mass_montecarlo, 
								"", "InvMass e^{+}e^{-} [GeV]",StandardYAxis,
								kFALSE, 3,0.000001,
								kFALSE,0. ,0., 
								kTRUE, 0.,0.1);
			}
		
		pad_c4->cd(4);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_Chi2_data,
								"", "#chi^{2} Conv  #gamma",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0 ,0, 
								kFALSE, 0,8.);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_Chi2_data,
								 ESD_ConvGamma_Chi2_montecarlo, 
								"", "#chi^{2} Conv  #gamma",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0. ,0., 
								kFALSE, 0.,8.);
			}
		
		pad_c4->cd(5);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_E_Pt_data,
								"", "p_{t} electron [GeV]",StandardYAxis,
								kTRUE, 3,0.00000001,
								kFALSE,0 ,0, 
								kTRUE, 0,20.);
			}else{
				DrawAutoGammaHistos( ESD_E_Pt_data,
								 ESD_E_Pt_montecarlo, 
								"", "p_{t} electron [GeV]",StandardYAxis,
								kTRUE, 3,0.00000001,
								kFALSE,0. ,0., 
								kTRUE, 0.,20.);
			}

		pad_c4->cd(6);
			if(MCfile == ""){
			DrawAutoGammaHisto( ESD_P_Pt_data,
							"", "p_{t} positron [GeV]",StandardYAxis,
							kTRUE, 3,0.00000001,
							kFALSE,0 ,0, 
							kTRUE, 0,20.);
			}else{
			DrawAutoGammaHistos( ESD_P_Pt_data,
							 ESD_P_Pt_montecarlo, 
							"", "p_{t} positron [GeV]",StandardYAxis,
							kTRUE, 3,0.00000001,
							kFALSE,0. ,0., 
							kTRUE, 0.,20.);
			}

		c4->Update();
		
		// --------------------------------Page 8 - Mother mass---------------------------------------------
		
		ps_characteristics->NewPage();
		
		TCanvas * c5 = new TCanvas("c5","",10,10,700,1000);  // gives the page size
		pad_c5 = new TPad("pad_c5","",0.05,0.05,0.95,0.85,0);   // gives the size of the histo areas 
		pad_c5->SetFillColor(0);
		pad_c5->GetFrame()->SetFillColor(0);
		pad_c5->SetBorderMode(0);
		pad_c5->Divide(2,2);
		pad_c5->Draw();
		
		title0->Draw();
		if(MCfile != ""){
			title1->Draw();
		}
	
		pad_c5->cd(1);
				if(MCfile == ""){
				DrawAutoGammaHisto( ESD_Mother_InvMass_data,
								"", "InvMass #gamma#gamma [GeV]","Counts/event scaled by multiplicity",
								kTRUE, 1.5,0,
								kFALSE,0 ,0, 
								kFALSE, 0,0);
			}else{
				DrawAutoGammaHistos( ESD_Mother_InvMass_data,
								 ESD_Mother_InvMass_montecarlo, 
								"","InvMass #gamma#gamma [GeV]","Counts/event scaled by multiplicity", 
								kTRUE, 1.5,0,
								kFALSE,0. ,0., 
								kFALSE, 0.,0.);
			}
				
		c5->Update();
		ps_characteristics->Close();
		
		
	}
	
	

	// ----------------------- Start Single Plots -----------------------------------------------------
	
	
	// ----------------------- dE/dx of pion/electron  -----------------------------------------
	if(SinglePlots){		
		//------------


		TCanvas * c1_1 = new TCanvas("c1_1","",10,10,500,500);  // gives the page size		
		c1_1->SetLogx(1);
		c1_1->SetLogz(1);

		c1_1->SetRightMargin(0.14); 		
		c1_1->cd();
			DrawAutoGammaHisto2D(	ESD_ConvGamma_P_dEdxP_data,
								"", "p [GeV/c]", "dE/dx", "",
								kTRUE, 10., 140.,
								kTRUE, 0., 20.);
			foElec.Draw("e,hist,same");
			foPion.Draw("e,hist,same");
			foProton.Draw("e,hist,same");
			foKaon.Draw("e,hist,same");
			DrawdEdxLabel();
		DrawAliceLogoPerformance(right_up2D[0],right_up2D[1],right_up2D[2],right_up2D[3],0.03, Date);	


		c1_1->Update();
		c1_1->SaveAs(Form("%s%s/dEdx_P_data.%s",path,suffix,suffix));
		delete c1_1;
		
		if(MCfile != ""){
			TCanvas * c1_2 = new TCanvas("c1_2","",10,10,500,500);  // gives the page size		
			
			c1_2->SetLogx(1);
			c1_2->SetLogz(1);
			c1_2->SetRightMargin(0.14); 			
			c1_2->cd();	
				DrawAutoGammaHisto2D(	ESD_ConvGamma_P_dEdxP_montecarlo,
									"", "p [GeV/c]", "dE/dx", "",
									kTRUE, 20., 150.,
									kTRUE, 0., 20.);
			foElecMC.Draw("e,hist,same");
			foPionMC.Draw("e,hist,same");
			foProtonMC.Draw("e,hist,same");
			foKaonMC.Draw("e,hist,same");
			DrawdEdxLabel();
			DrawAliceLogo(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3]);
			c1_2->Update();
			c1_2->SaveAs(Form("%s%s/dEdx_P_montecarlo.%s",path,suffix,suffix));
			delete c1_1;
		}

		TCanvas * c1_1 = new TCanvas("c1_1","",10,10,500,500);  // gives the page size		
		c1_1->SetLogx(1);
		c1_1->SetLogz(1);

		c1_1->SetRightMargin(0.14); 		
		c1_1->cd();
			DrawAutoGammaHisto2D(	ESD_ConvGamma_E_dEdxP_data,
								"", "p [GeV/c]", "dE/dx", "",
								kTRUE, 10., 140.,
								kTRUE, 0., 20.);
			foElec.Draw("e,hist,same");
			foPion.Draw("e,hist,same");
			foProton.Draw("e,hist,same");
			foKaon.Draw("e,hist,same");
			DrawdEdxLabel();
		DrawAliceLogoPerformance(right_up2D[0],right_up2D[1],right_up2D[2],right_up2D[3],0.03, Date);	


		c1_1->Update();
		c1_1->SaveAs(Form("%s%s/dEdx_E_data.%s",path,suffix,suffix));
		delete c1_1;
		
		if(MCfile != ""){
			TCanvas * c1_2 = new TCanvas("c1_2","",10,10,500,500);  // gives the page size		
			
			c1_2->SetLogx(1);
			c1_2->SetLogz(1);
			c1_2->SetRightMargin(0.14); 			
			c1_2->cd();	
				DrawAutoGammaHisto2D(	ESD_ConvGamma_E_dEdxP_montecarlo,
									"", "p [GeV/c]", "dE/dx", "",
									kTRUE, 20., 150.,
									kTRUE, 0., 20.);
			foElecMC.Draw("e,hist,same");
			foPionMC.Draw("e,hist,same");
			foProtonMC.Draw("e,hist,same");
			foKaonMC.Draw("e,hist,same");
			DrawdEdxLabel();
			DrawAliceLogo(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3]);
			c1_2->Update();
			c1_2->SaveAs(Form("%s%s/dEdx_E_montecarlo.%s",path,suffix,suffix));
			delete c1_1;
		}

		
		// ------------------------ asymetry ----------------------------------------------------
		
		
		
		TCanvas * c2_1 = new TCanvas("c2_1","",10,10,500,500);  // gives the page size
		c2_1->SetLogx(0);
		c2_1->SetLogz(1);
		c2_1->SetRightMargin(0.17); 			
		c2_1->cd();
			DrawAutoGammaHisto2D(	ESD_ConvGamma_P_AsymmetryP_data,
								"", "p [GeV/c]", "Asymmetry", "",
								kFALSE, 30., 100.,
								kTRUE, 0., 6.);
			DrawAliceLogo(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3]);
		c2_1->Update();
		c2_1->SaveAs(Form("%s%s/Asymetry_data.%s",path,suffix,suffix));
		delete c2_1;

		if(MCfile != ""){			
			TCanvas * c2_2 = new TCanvas("c2_2","",10,10,500,500);  // gives the page size			
			c2_2->SetLogx(0);
			c2_2->SetLogz(1);
			c2_2->SetRightMargin(0.17); 			
			c2_2->cd();
				DrawAutoGammaHisto2D(	ESD_ConvGamma_P_AsymmetryP_montecarlo,
									"", "p [GeV/c]", "Asymmetry", "",
									kFALSE, 30., 100.,
									kTRUE, 0., 6.);

			DrawAliceLogo(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3]);
			c2_2->Update();
			c2_2->SaveAs(Form("%s%s/Asymetry_montecarlo.%s",path,suffix,suffix));
			delete c2_2;	
		}
// ---------------------------- Armenteros single plots  ------------------------------------------------------------------
		
		TCanvas * c20_1 = new TCanvas("c20_1","",10,10,500,450);  // gives the page size
		c20_1->SetLogx(0);
		c20_1->SetLogz(1);
		c20_1->SetRightMargin(0.17); 			
		c20_1->cd();
			ESD_AllV0CurrentFinder_alfa_qt_data->SetMinimum(20);
			DrawAutoGammaHisto2D(	ESD_AllV0CurrentFinder_alfa_qt_data,
								"", "#alpha = (p^{+}_{L}-p^{-}_{L})/(p^{+}_{L}+p^{-}_{L})", "q_{T} [GeV/c]", "",
								kFALSE, 30., 100.,
								kTRUE, -1., 1.);
			DrawArmenteros();
			DrawAliceLogoPerformance(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3],0.03,Date);
		c20_1->Update();
		c20_1->SaveAs(Form("%s%s/Armenterosfull_data.%s",path,suffix,suffix));
		delete c20_1;

		if(MCfile != ""){			
			TCanvas * c20_2 = new TCanvas("c20_2","",10,10,500,450);  // gives the page size			
			c20_2->SetLogx(0);
			c20_2->SetLogz(1);
			c20_2->SetRightMargin(0.17); 			
			c20_2->cd();
				DrawAutoGammaHisto2D(	ESD_AllV0CurrentFinder_alfa_qt_montecarlo,
									"","#alpha = (p^{+}_{L}-p^{-}_{L})/(p^{+}_{L}+p^{-}_{L})", "q_{T} [GeV/c]", "",
									kFALSE, 30., 100.,
									kTRUE, -1., 1.);
				DrawArmenteros();
				DrawAliceLogoPerformance(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3],0.03, Date);
			c20_2->Update();
			c20_2->SaveAs(Form("%s%s/Armenterosfull_montecarlo.%s",path,suffix,suffix));
			delete c20_2;	
		}

		TCanvas * c20_1 = new TCanvas("c20_1","",10,10,500,450);  // gives the page size
		c20_1->SetLogx(0);
		c20_1->SetLogz(1);
		c20_1->SetRightMargin(0.17); 			
		c20_1->cd();
			DrawAutoGammaHisto2D(	ESD_ConvGamma_alfa_qt_data,
								"", "#alpha = (p^{+}_{L}-p^{-}_{L})/(p^{+}_{L}+p^{-}_{L})", "q_{T} [GeV/c]", "",
								kFALSE, 30., 100.,
								kTRUE, -1., 1.);
			DrawAliceLogoPerformance(right_up2D[0],right_up2D[1],right_up2D[2],right_up2D[3],0.03,Date);
		c20_1->Update();
		c20_1->SaveAs(Form("%s%s/ArmenterosAftCuts_data.%s",path,suffix,suffix));
		delete c20_1;

		if(MCfile != ""){			
			TCanvas * c20_2 = new TCanvas("c20_2","",10,10,500,450);  // gives the page size			
			c20_2->SetLogx(0);
			c20_2->SetLogz(1);
			c20_2->SetRightMargin(0.17); 			
			c20_2->cd();
				DrawAutoGammaHisto2D(	ESD_ConvGamma_alfa_qt_montecarlo,
									"", "#alpha = (p^{+}_{L}-p^{-}_{L})/(p^{+}_{L}+p^{-}_{L})", "q_{T} [GeV/c]", "",
									kFALSE, 30., 100.,
									kTRUE, -1., 1.);
				DrawAliceLogoPerformance(right_up2D [0],right_up2D[1],right_up2D[2],right_up2D[3],0.03,Date);
			c20_2->Update();
			c20_2->SaveAs(Form("%s%s/ArmenterosAftCuts_montecarlo.%s",path,suffix,suffix));
			delete c20_2;	
		}


		
		// -----------------------------------------------------------------------------
		
		
		TCanvas * c3_1 = new TCanvas("c3_1","",10,10,500,500);  // gives the page size
		TCanvas * c3_2 = new TCanvas("c3_2","",10,10,500,500);  // gives the page size
		TCanvas * c3_3 = new TCanvas("c3_3","",10,10,500,500);  // gives the page size
		TCanvas * c3_4 = new TCanvas("c3_4","",10,10,500,500);  // gives the page size
		
		
		c3_1->SetLogy(1);
		c3_2->SetLogy(1);
		c3_3->SetLogy(1);
		c3_4->SetLogy(1);
		
//--------------------------------------- DCA -------------------------------------------------------------		
		c3_1->cd();
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_DcaDaughters_data,
								"", "DCA [cm]",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0 ,0, 
								kFALSE, 0,1);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_DcaDaughters_data,
								 ESD_ConvGamma_DcaDaughters_montecarlo, 
								 "", "DCA [cm]",StandardYAxis,
								 kTRUE, 3,0.000001,
								 kFALSE,0. ,0., 
								 kFALSE, 0.,1);
			}
			DrawAliceLogoPerformance(right_up [0],right_up[1],right_up[2],right_up[3],0.03, Date);	

		c3_1->Update();
		c3_1->SaveAs(Form("%s%s/Photon_DcaDaughters.%s",path,suffix,suffix));
		delete c3_1;


//-------------------------------- Normalised Distance ----------------------------------------
		c3_2->cd();
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_NormDcaDistDaughters_data,
								"", "Normalized distance [cm]",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0 ,0, 
								kFALSE, 0,1);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_NormDcaDistDaughters_data,
								 ESD_ConvGamma_NormDcaDistDaughters_montecarlo, 
								"", "Normalized distance [cm]",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0. ,0., 
								kFALSE, 0.,1);
			}
			DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);

		c3_2->Update();
		c3_2->SaveAs(Form("%s%s/Photon_NormDcaDaughters.%s",path,suffix,suffix));
		delete c3_2;
//------------------------------------- Cos Pointing Angle -----------------------------------------		
		c3_3->cd();
		if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_CosPointingAngle_data,
								"", "Cos Pointing Angle",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0 ,0, 
								kTRUE, 0,1);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_CosPointingAngle_data,
								 ESD_ConvGamma_CosPointingAngle_montecarlo, 
								"", "Cos Pointing Angle",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0. ,0., 
								kTRUE, 0.,1);
			}
			DrawAliceLogoPerformance(left_up [0],left_up[1],left_up[2],left_up[3],0.03,Date);
		
		c3_3->Update();
		c3_3->SaveAs(Form("%s%s/Photon_CosPointingAngle.%s",path,suffix,suffix));
		delete c3_3;	

//-------------------------------------- Likelihood -----------------------------------------------------			
		c3_4->cd();
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_LikeLihood_data,
								"", "Likelihood",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0 ,0, 
								kTRUE, 0,1);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_LikeLihood_data,
								 ESD_ConvGamma_LikeLihood_montecarlo, 
								"", "Likelihood",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0. ,0., 
								kTRUE, 0.,1);
			}
			DrawAliceLogoPerformance(left_up [0],left_up[1],left_up[2],left_up[3],0.03,Date);

		c3_4->Update();
		c3_4->SaveAs(Form("%s%s/Photon_Likelihood.%s",path,suffix,suffix));
		delete c3_4;
		
// -------------------------------------------------------------------------------------------------
		TCanvas * c4_1 = new TCanvas("c4_1","",10,10,500,500);  // gives the page size
		TCanvas * c4_2 = new TCanvas("c4_2","",10,10,500,500);  // gives the page size
		TCanvas * c4_3 = new TCanvas("c4_3","",10,10,500,500);  // gives the page size
		TCanvas * c4_4 = new TCanvas("c4_4","",10,10,500,500);  // gives the page size
		TCanvas * c4_5 = new TCanvas("c4_5","",10,10,500,500);  // gives the page size
		TCanvas * c4_6 = new TCanvas("c4_6","",10,10,500,500);  // gives the page size
		
		
		c4_1->SetLogy(1);
		c4_2->SetLogy(1);
		c4_3->SetLogy(1);
		c4_4->SetLogy(1);
		c4_5->SetLogy(1);
		c4_6->SetLogy(1);
		
//------------------------------------ Photon Momentum -----------------------------------------		
		c4_1->cd();
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_Pt_data,
								"", "p_{t} Conv #gamma [GeV]",StandardYAxis,
								kTRUE, 3,0.00000001,
								kFALSE,0 ,0, 
								kTRUE, 0,40.);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_Pt_data,
								 ESD_ConvGamma_Pt_montecarlo, 
								"", "p_{t} Conv #gamma [GeV]",StandardYAxis,
								kTRUE, 3,0.0000000001,
								kFALSE,0. ,0., 
								kTRUE, 0.,40.);
			}
			DrawAliceLogoPerformance(right_up [0],right_up[1],right_up[2],right_up[3],0.03,Date);

		c4_1->Update();
		c4_1->SaveAs(Form("%s%s/Photon_Pt.%s",path,suffix,suffix));
		delete c4_1;
		
//----------------------------------- Photon eta ---------------------------------------------
		c4_2->cd();
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_Eta_data,
								"", "#eta Conv  #gamma",StandardYAxis,
								kTRUE, 3,0.001,
								kFALSE,0 ,0, 
								kFALSE, 0,0.7);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_Eta_data,
								 ESD_ConvGamma_Eta_montecarlo, 
								"", "#eta Conv  #gamma",StandardYAxis,
								kTRUE, 3,0.001,
								kFALSE,0. ,0., 
								kFALSE, 0.,0.7);
		}
		DrawAliceLogoPerformance(left_up [0],left_up[1],left_up[2],left_up[3],0.03,Date);           		
		c4_2->Update();
		c4_2->SaveAs(Form("%s%s/Photon_Eta.%s",path,suffix,suffix));
		delete c4_2;

//------------------------------------- Photon Mass ----------------------------------------
		c4_3->cd();
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_Mass_data,
								"", "InvMass e^{+}e^{-} [GeV]",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0 ,0, 
								kTRUE, 0,0.1);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_Mass_data,
								 ESD_ConvGamma_Mass_montecarlo, 
								"", "InvMass e^{+}e^{-} [GeV]",StandardYAxis,
								kTRUE, 3,0.00000001,
								kFALSE,0. ,0., 
								kTRUE, 0.,0.15);
			}
		DrawAliceLogoPerformance(right_up [0],right_up[1],right_up[2],right_up[3],0.03,Date);
		c4_3->Update();
		c4_3->SaveAs(Form("%s%s/Photon_Mass.%s",path,suffix,suffix));
		delete c4_3;		

//-------------------------------------- Photon Chi2 ----------------------------------------------
		c4_4->cd();
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_ConvGamma_Chi2_data,
								"", "#chi^{2} Conv  #gamma",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0 ,0, 
								kFALSE, 0,10.);
			}else{
				DrawAutoGammaHistos( ESD_ConvGamma_Chi2_data,
								 ESD_ConvGamma_Chi2_montecarlo, 
								"", "#chi^{2} Conv  #gamma",StandardYAxis,
								kTRUE, 3,0.000001,
								kFALSE,0. ,0., 
								kFALSE, 0.,10.);
			}

		DrawAliceLogoPerformance(right_up [0],right_up[1],right_up[2],right_up[3],0.03,Date);

		c4_4->Update();
		c4_4->SaveAs(Form("%s%s/Photon_Chi2.%s",path,suffix,suffix));
		delete c4_4;		
		

//---------------------------------- Electron Momentum -------------------------------------------		
		c4_5->cd();
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_E_Pt_data,
								"", "p_{t} electron [GeV]",StandardYAxis,
								kTRUE, 3,0.00000001,
								kFALSE,0 ,0, 
								kTRUE, 0,20.);
			}else{
				DrawAutoGammaHistos( ESD_E_Pt_data,
								 ESD_E_Pt_montecarlo, 
								"", "p_{t} electron [GeV]",StandardYAxis,
								kTRUE, 3,0.00000001,
								kFALSE,0. ,0., 
								kTRUE, 0.,20.);
			}
			DrawAliceLogoPerformance(right_up [0],right_up[1],right_up[2],right_up[3],0.03,Date);
		
		c4_5->Update();
		c4_5->SaveAs(Form("%s%s/Electron_Pt.%s",path,suffix,suffix));
		delete c4_5;
//---------------------------------- Positron Momentum -------------------------------------------
		c4_6->cd();
		if(MCfile == ""){
			DrawAutoGammaHisto( ESD_P_Pt_data,
							"", "p_{t} positron [GeV]",StandardYAxis,
							kTRUE, 3,0.00000001,
							kFALSE,0 ,0, 
							kTRUE, 0,20.);
		}else{
			DrawAutoGammaHistos( ESD_P_Pt_data,
							 ESD_P_Pt_montecarlo, 
							"", "p_{t} positron [GeV]",StandardYAxis,
							kTRUE, 3,0.00000001,
							kFALSE,0. ,0., 
							kTRUE, 0.,20.);
		}
		DrawAliceLogoPerformance(right_up [0],right_up[1],right_up[2],right_up[3],0.03,Date);

		c4_6->Update();
		c4_6->SaveAs(Form("%s%s/Positron_Pt.%s",path,suffix,suffix));
		delete c4_6;
		
		
		// -----------------------------------------------------------------------------
		
		
		TCanvas *c5_1 = new TCanvas("c5","",10,10,500,500);  // gives the page size
		
		
		c5_1->cd();
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_Mother_InvMass_data,
								"", "InvMass #gamma#gamma [GeV]","counts/event scaled by multiplicity",
								kTRUE, 1.5,0,
								kFALSE,0 ,0, 
								kFALSE, 0,0);
			}else{
				DrawAutoGammaHistos( ESD_Mother_InvMass_data,
								 ESD_Mother_InvMass_montecarlo, 
								"","InvMass #gamma#gamma [GeV]","counts/event scaled by multiplicity", 
								kTRUE, 1.5,0,
								kFALSE,0. ,0., 
								kFALSE, 0.,0.);
			}

			DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);		
		c5_1->Update();
		c5_1->SaveAs(Form("%s%s/Mother_InvMass.%s",path,suffix,suffix));
		delete c5_1;

// ---------------------------------------- Multiplicity Plot ---------------------------------------------------
		TCanvas * c6_1 = new TCanvas("c6_1","",10,10,500,500);  // gives the page size
		c6_1->SetLogy(1);		
		TText *l = new TText();

		c6_1->cd(1);
			if(MCfile == ""){
				DrawAutoGammaHisto( ESD_NumberOfGoodESDTracks_data, 
								 "","Number good ESDTracks ","",
								 kTRUE, 3.,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,100.);
				Double_t mean_data = ESD_NumberOfGoodESDTracks_data->GetMean();
				l->SetTextAlign(11);
				l->SetTextSize(0.03);
				l->DrawText(60,0.0002, "Mean Data:");
				l->DrawText(80,0.0002, Form("%lg",  mean_data));	

			}else{
				DrawAutoGammaHistos( ESD_NumberOfGoodESDTracks_data, 
								 ESD_NumberOfGoodESDTracks_montecarlo, 
								 "","Number good ESDTracks","",
								 kTRUE, 3.,0.000001,
								 kFALSE,0. ,0.,
								 kFALSE, 0.,100.);
				Double_t mean_data = ESD_NumberOfGoodESDTracks_data->GetMean();
				l->SetTextAlign(11);
				l->SetTextSize(0.03);
				l->DrawText(60,0.0002, "Mean Data:");
				l->DrawText(80,0.0002, Form("%lg",  mean_data));	
	
				Double_t mean_montecarlo = ESD_NumberOfGoodESDTracks_montecarlo->GetMean();
				l->SetTextAlign(11);
				l->SetTextColor(2);			
				l->SetTextSize(0.03);
				l->DrawText(60,0.0001, "Mean MC:");
				l->DrawText(80,0.0001, Form("%lg",  mean_montecarlo));		
	
			}
			DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);

		c6_1->Update();
		c6_1->SaveAs(Form("%s%s/NumberOfGoodESDTracks.%s",path,suffix,suffix));
		delete c6_1;

//------------------------------------ Projections dEdx -----------------------------------------------------------------		
		TCanvas * c9 = new TCanvas("c9","",10,10,1400,700);  // gives the page size
		pad_c9 = new TPad("pad_c9","",0.05,0.05,0.95,0.95,0);   // gives the size of the histo areas 
		pad_c9->SetFillColor(0);
		pad_c9->GetFrame()->SetFillColor(0);
		pad_c9->SetBorderMode(0);
		pad_c9->Divide(2,1);
		pad_c9->Draw();
		
		pad_c9->cd(1)->SetLogx(1);
		pad_c9->cd(1)->SetLogz(1);
		pad_c9->cd(2)->SetLogx(1);
		pad_c9->cd(2)->SetLogz(1);
			
		pad_c9->cd(1);
		if(MCfile == ""){
			ESD_ConvGamma_P_dEdx_Y1_data->Rebin(2);
			DrawAutoGammaHisto( ESD_ConvGamma_P_dEdx_Y1_data, 
							 "","dE/dx ",StandardYAxis,
							 kTRUE, 1.2,0,
							 kFALSE,0. ,0.,
							 kTRUE, 50.,100.);
		}else{
			ESD_ConvGamma_P_dEdx_Y1_data->Rebin(2);
			ESD_ConvGamma_P_dEdx_Y1_montecarlo_2->Rebin(2);
			DrawAutoGammaHistos( ESD_ConvGamma_P_dEdx_Y1_data, 
							 ESD_ConvGamma_P_dEdx_Y1_montecarlo_2, 
							 "","dE/dx ",StandardYAxis,
							 kTRUE, 1.2,0,
							 kFALSE,0. ,0.,
							 kTRUE, 50.,100.);

		}
		Scaling->Draw();		
		DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);					

 		pad_c9->cd(2);
		if(MCfile == ""){
			ESD_ConvGamma_P_dEdx_Y2_data->Rebin(2);
			DrawAutoGammaHisto( ESD_ConvGamma_P_dEdx_Y2_data, 
							 "","dE/dx ",StandardYAxis,
							 kTRUE, 1.2,0,
							 kFALSE,0. ,0.,
							 kTRUE, 50.,100.);
		}else{
			ESD_ConvGamma_P_dEdx_Y2_data->Rebin(2);
			ESD_ConvGamma_P_dEdx_Y2_montecarlo_2->Rebin(2);
			DrawAutoGammaHistos( ESD_ConvGamma_P_dEdx_Y2_data, 
							 ESD_ConvGamma_P_dEdx_Y2_montecarlo_2, 
							 "","dE/dx ",StandardYAxis,
							 kTRUE, 1.2,0,
							 kFALSE,0. ,0.,
							 kTRUE, 50.,100.);
		}
		Scaling->Draw();		
		DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);	

		
		c9->Update();
		c9->SaveAs(Form("%s%s/ProjectiondEdx.%s",path,suffix,suffix));
		delete pad_c9; 
		delete c9;


// ------------------------------- Asymmetry Projections -------------------------------------------------------------		
		TCanvas * c10 = new TCanvas("c10","",10,10,1400,700);  // gives the page size
		pad_c10 = new TPad("pad_c2","",0.0,0.0,1,1,0);   // gives the size of the histo areas 
		pad_c10->SetFillColor(0);
		pad_c10->GetFrame()->SetFillColor(0);
		pad_c10->SetBorderMode(0);
		pad_c10->Divide(2,1);
		pad_c10->Draw();

		leg1 = new TLegend( 0.65,0.87,0.98,0.97);
		leg1->SetTextSize(0.025);			
		leg1->SetFillColor(0);
		leg1->AddEntry(ESD_ConvGamma_P_AsymmetryP_Y1_data,("Data"));
		leg1->AddEntry(ESD_ConvGamma_P_AsymmetryP_Y1_montecarlo,("MC reconstructed"));
		leg1->AddEntry(MC_ConvGamma_P_AsymmetryP_Y1_montecarlo,("MC truth"));

		//Asymmetrie in 1st and 3rd bin		
		pad_c10->cd(1);
		pad_c10->cd(1)->SetTopMargin(0.03);
		pad_c10->cd(1)->SetBottomMargin(1);
		pad_c10->cd(1)->SetRightMargin(0.02);
		pad_c10->cd(1)->SetLeftMargin(0.16);
		if(MCfile == ""){
			DrawAutoGammaHisto( ESD_ConvGamma_P_AsymmetryP_Y1_data, 
							"", "x = E/k","(X_{0}N_{A}/A)d#sigma_{LPM}/dx",
							kTRUE, 2,0,
							kFALSE,0. ,0., 
							kFALSE, 0.,0.9);

		}else{
			DrawAutoGammaHistos( ESD_ConvGamma_P_AsymmetryP_Y1_data,
							 ESD_ConvGamma_P_AsymmetryP_Y1_montecarlo,
							 "", "x = E/k","(X_{0}N_{A}/A)d#sigma_{LPM}/dx", 
							 kTRUE, 2,0,
							 kFALSE,0. ,0.,
							 kFALSE, 0.1,0.9);
			MC_ConvGamma_P_AsymmetryP_Y1_montecarlo->SetLineColor(kBlue+1);
			MC_ConvGamma_P_AsymmetryP_Y1_montecarlo->Draw("same,hist");
//			leg1->Draw("same");
		}
		DrawAliceLogoPerformance(right_up[0]+0.1,right_up[1]+0.17,right_up[2],right_up[3],0.03, Date);

		//Asymmetry in other bins
		pad_c10->cd(2);
		pad_c10->cd(2)->SetTopMargin(0.03);
		pad_c10->cd(2)->SetBottomMargin(1);
		pad_c10->cd(2)->SetRightMargin(0.02);
		pad_c10->cd(2)->SetLeftMargin(0.16);
		if(MCfile == ""){
			DrawAutoGammaHisto( ESD_ConvGamma_P_AsymmetryP_Y2_data,
							"", "x = E/k","(X_{0}N_{A}/A)d#sigma_{LPM}/dx", 
							kTRUE, 2,0,
							kFALSE,0. ,0., 
							kFALSE, 0.,1.);
		}else{

			DrawAutoGammaHistos( ESD_ConvGamma_P_AsymmetryP_Y2_data, 
							 ESD_ConvGamma_P_AsymmetryP_Y2_montecarlo, 
							 "", "x = E/k","(X_{0}N_{A}/A)d#sigma_{LPM}/dx", 
							 kTRUE, 2,0,
							 kFALSE,0. ,0., 
							 kFALSE, 0.,1.);
			MC_ConvGamma_P_AsymmetryP_Y2_montecarlo->SetLineColor(kBlue+1);
			MC_ConvGamma_P_AsymmetryP_Y2_montecarlo->Draw("same,hist");
			leg1->Draw("same");
		}
//		DrawAliceLogo(right_up [0],right_up[1],right_up[2],right_up[3]);
		c10->Update();
		c10->SaveAs(Form("%s%s/ProjectionAsymmetry.%s",path,suffix,suffix));
		delete pad_c10;
		delete c10;		


// ---------------------------------------------- Clusters ---------------------------------------------------------
		TCanvas * c11 = new TCanvas("c11","",10,10,2000,1400);  // gives the page size
		pad_c11 = new TPad("pad_c11","",0.0,0.0,1,1,0);   // gives the size of the histo areas 
		pad_c11->SetFillColor(0);
		pad_c11->GetFrame()->SetFillColor(0);
		pad_c11->SetBorderMode(0);
		
		pad_c11->Divide(2,2);
		pad_c11->Draw();
		
		// ITS Clusters Electron
		pad_c11->cd(1);
		if(MCfile == ""){
			DrawAutoGammaHisto( ESD_E_nITSClusters_data,
							"", "Number of ITS Clusters for Electron ","Hits/event scaled by multiplicity",
							kTRUE, 1.5,0,
							kFALSE,0 ,0, 
							kFALSE, 0,0);
		}else{
			DrawAutoGammaHistos( ESD_E_nITSClusters_data,
							 ESD_E_nITSClusters_montecarlo, 
							"","Number of ITS Clusters for Electron ","Hits/event scaled by multiplicity", 
							kTRUE, 1.5,0,
							kFALSE,0. ,0., 
							kFALSE, 0.,0.);
		}
		DrawAliceLogoPerformance(right_up [0],right_up[1],right_up[2],right_up[3],0.03,Date);
		//TPC Clusters Electron
		pad_c11->cd(2);
		if(MCfile == ""){
			DrawAutoGammaHisto( ESD_E_nTPCClusters_data,
							"", "Number of TPC Clusters for Electron","Hits/event scaled by multiplicity", 
							kTRUE, 1.5,0,
							kFALSE,0. ,0., 
							kFALSE, 0.,0.);
		}else{
			DrawAutoGammaHistos( ESD_E_nTPCClusters_data, 
							 ESD_E_nTPCClusters_montecarlo,
							 "","Number of TPC Clusters for Electron ","Hits/event scaled by multiplicity",
							 kTRUE, 1.5,0,
							 kFALSE,0. ,0.,
							 kFALSE, 0.,0.);
		}
		//ITS Clusters Positron
		pad_c11->cd(3);
		if(MCfile == ""){
			DrawAutoGammaHisto( ESD_P_nITSClusters_data,
							"", "Number of ITS Clusters for PositronS","Hits/event scaled by multiplicity", 
							kTRUE, 1.5, 0,
							kFALSE,0. ,0., 
							kFALSE, 0.,0.);
		}else{
			DrawAutoGammaHistos( ESD_P_nITSClusters_data, 
							 ESD_P_nITSClusters_montecarlo,
							 "","Number of ITS Clusters for Positron","Hits/event scaled by multiplicity", 
							 kTRUE, 1.5, 0,
							 kFALSE,0. ,0.,
							 kFALSE, 0.,0.);
		}
		//TPC Clusters Positron
		pad_c11->cd(4);
		if(MCfile != ""){	
			DrawAutoGammaHistos( ESD_P_nTPCClusters_data, 
							 ESD_P_nTPCClusters_montecarlo, 
							 "","Number of TPC Clusters for Positron","Hits/event scaled by multiplicity",
							 kTRUE, 1.5, 0,
							 kFALSE,0. ,0.,
							 kFALSE, 0.,0.);
			
		}else{
			DrawAutoGammaHisto( ESD_P_nTPCClusters_data,
							"", "Number of TPC Clusters for PositronS","Hits/event scaled by multiplicity", 
							kTRUE, 1.5, 0,
							kFALSE,0. ,0., 
							kFALSE, 0.,0.);			
		}
			
		c11->Update();
		c11->SaveAs(Form("%s%s/Clusters.%s",path,suffix,suffix));
		delete pad_c11;
		delete c11;		
	//end drawing single Plots
	}	

}

