// provided by Gamma Conversion Group, PWG4, Kathrin Koch, kkoch@physi.uni-heidelberg.de

#include <fstream>
#include <Riostream.h>
//#include <PlottingGammaHistograms.h>

extern TRandom *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem *gSystem;

void Extract_IntegratedPi0Yield(const char* cutSelection="", const char *inputRootFile = "AnalysisResults",const char *path = "./",const char* outputDir="./Output/",const char* outputDatFileBase= "RB-data",const char* suffix=""){	

  Bool_t writePlotFiles=kFALSE; // if one wants to write f.ex .gif files from the plots
  if(suffix!=""){
    writePlotFiles=kTRUE;
  }

  gROOT->Reset();	
  gROOT->SetStyle("Plain");
	
  // rebinning
  Int_t RebinMass = 2;
	
  Double_t BGFit_range[2] = {0.17,0.3};   	
	
  Double_t Pi0Mass = 0.135;
  Double_t Pi0Mass_range[2] = {0.12,0.15};
  Double_t Pi0Width = 0.003;
  Double_t Pi0Width_range[2] = {0.001,0.007};
  Double_t Pi0Slope = 0.007;
  Double_t Pi0Slope_range[2] = {0.001,0.016};
	
  Double_t Pi0FitRange[2] = {0.1,0.16};
	
  Bool_t MC = kFALSE;

	
  // --------------------------------- end of self definitions --------------------------------

  TString filename = Form("%s%s.root",path,inputRootFile);	
  TFile f(filename.Data());  

  TCanvas* c_mass = new TCanvas("c_mass","",200,10,600,600);  // gives the page size
  c_mass->SetFillColor(0);

      
  TDirectory *pwg4dir =(TDirectory*)f.Get(Form("PWG4_GammaConversion_%s",cutSelection));
  TList *fHistosGammaConversion = (TList*)pwg4dir->Get(Form("histogramsAliGammaConversion_%s",cutSelection));
  TList *fESDContainer = (TList*)fHistosGammaConversion->FindObject("ESD histograms");
  TList *fMCContainer = (TList*)fHistosGammaConversion->FindObject("MC histograms");
  TList *fBackgroundContainer = (TList*)fHistosGammaConversion->FindObject("Background histograms");
  TH1F *ESD_Mother_Mass= (TH1F*)fESDContainer->FindObject("ESD_Mother_InvMass");
	    
  TH1F *ESD_Background_Mass=(TH1F*)fBackgroundContainer->FindObject("ESD_Background_InvMass");
	    
  TH1F * ESD_NumberOfContributorsVtx=(TH1F*)fESDContainer->FindObject("ESD_NumberOfContributorsVtx");
  ESD_NumberOfContributorsVtx->SetAxisRange(1.,100.);
	    
  Float_t nGoodEvents = ESD_NumberOfContributorsVtx->Integral();

  ESD_Mother_Mass->Sumw2();
  ESD_Background_Mass->Sumw2();

  // for Pi0 reco efficiency
  Double_t nPi0_MC = 0;
  Double_t nPi0_MC_error = 0;
  if(MC){
    TH1D *MC_Pi0_Pt=fMCContainer->FindObject("MC_Pi0_Pt");
    nPi0_MC = MC_Pi0_Pt->GetEntries();
    nPi0_MC_error = sqrt(nPi0_MC);
  }

  // get counts for reco and background within pi0 range
  TAxis *xaxis_reco = ESD_Mother_Mass->GetXaxis();							
  Int_t r_1 = xaxis_reco->FindBin(BGFit_range[0]);
  Int_t r_2 = xaxis_reco->FindBin(BGFit_range[1]);  
  Double_t r =  ESD_Mother_Mass->Integral(r_1,r_2); // Integral(75,125)
  TAxis *xaxis_back =ESD_Background_Mass->GetXaxis();							
  Int_t b_1 = xaxis_back->FindBin(BGFit_range[0]);
  Int_t b_2 = xaxis_back->FindBin(BGFit_range[1]);   
  Double_t b =  ESD_Background_Mass->Integral(b_1,b_2);
  Double_t norm = 1;
  if(b != 0) norm = r/b;
  ESD_Background_Mass->Sumw2();		
  ESD_Background_Mass->Scale(norm);
	    
  ESD_Mother_Mass->Rebin(RebinMass);
  ESD_Background_Mass->Rebin(RebinMass);

  ESD_Mother_Mass->SetTitle(Form("%s - cut: %s",inputRootFile,cutSelection));
  ESD_Mother_Mass->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  ESD_Mother_Mass->GetXaxis()->SetRangeUser(0.0,0.3);	
  ESD_Mother_Mass->SetLineColor(2);
  ESD_Mother_Mass->Draw("");
	    
  ESD_Background_Mass->GetXaxis()->SetRangeUser(0.0,0.3);	
  ESD_Background_Mass->SetLineColor(4);
  ESD_Background_Mass->Draw("same");
	    
  if(writePlotFiles){
    c_mass->Print(Form("%sSpectra-%s-%s.%s",outputDir,inputRootFile,cutSelection,suffix));
  }
  c_mass->Clear();

  TAxis *xaxis_Reco = ESD_Mother_Mass->GetXaxis();	
  Double_t ilowmassPi0_bin = xaxis_Reco->FindBin(0.1);  
  Double_t ihighmassPi0_bin = xaxis_Reco->FindBin(0.16);	
  Double_t Reco_error;
  Double_t Reco = ESD_Mother_Mass->IntegralAndError(ilowmassPi0_bin,ihighmassPi0_bin,Reco_error);
  //	    cout << "Reco:   " << Reco << "      error:  " << Reco_error << endl;
	    
  TAxis *xaxis_Back = ESD_Background_Mass->GetXaxis();	
  Double_t ilowmassPi0_bin = xaxis_Back->FindBin(0.1);  
  Double_t ihighmassPi0_bin = xaxis_Back->FindBin(0.16);	
  Double_t Back_error;
  Double_t Back = ESD_Background_Mass->IntegralAndError(ilowmassPi0_bin,ihighmassPi0_bin,Back_error);
  //	    cout << "Back:   " << Back << "      error:  " << Back_error << endl;
	    
		
  TH1D* Signal = (TH1D*)ESD_Mother_Mass->Clone();
  Signal->Sumw2();	
  Signal->Add(ESD_Background_Mass,-1.);
  Signal->Draw();
	    
  Double_t Pi0Amplitude = Signal->GetMaximum();
  Double_t Pi0Amplitude_Min = Pi0Amplitude*80/100;
  Double_t Pi0Amplitude_Max = Pi0Amplitude*100/100;

  TF1 *fReco = new TF1("Exp+Gauss","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",Pi0FitRange[0],Pi0FitRange[1]); 
	    
  // fit the peak around pi0 mass and draw		
  fReco->SetParameter(0,Pi0Amplitude);
  fReco->SetParameter(1,Pi0Mass);
  fReco->SetParameter(2,Pi0Width);
  fReco->SetParameter(3,Pi0Slope);
	    
  fReco->SetParLimits(0,Pi0Amplitude_Min,Pi0Amplitude_Max);
  fReco->SetParLimits(1,Pi0Mass_range[0],Pi0Mass_range[1]);
  fReco->SetParLimits(2,Pi0Width_range[0],Pi0Width_range[1]);
  fReco->SetParLimits(3,Pi0Slope_range[0],Pi0Slope_range[1]);
	    
  Signal->Fit(fReco,"RME");
  fReco->SetLineColor(3);
  fReco->SetLineWidth(0.7);
  fReco->SetLineStyle(4);		
  fReco->Draw("same");
	    
  Double_t parameter[4];
  fReco->GetParameters(parameter);
	    
  Double_t Mass = parameter[1];
  Double_t Mass_error = fReco->GetParError(1);
	    
  // calculation of FWHM + error
  Int_t Pi0Amplitude_fit = fReco->GetMaximum(Pi0FitRange[0],Pi0FitRange[1]);
  Int_t Pi0Amplitude_bin = fReco->GetXaxis()->FindBin(Pi0Amplitude_fit);
  Double_t MassFWHM1 = fReco->GetX(Pi0Amplitude_fit*0.5,0.1,0.135);
  Double_t MassFWHM2 = fReco->GetX(Pi0Amplitude_fit*0.5,0.135,0.16);		
  Double_t FWHM = MassFWHM2 - MassFWHM1;

  // + error
  TF1 *fReco_EPlus = new TF1("Exp+Gauss","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",Pi0FitRange[0],Pi0FitRange[1]); 
  fReco_EPlus->SetParameter(0,parameter[0] + fReco->GetParError(0));
  fReco_EPlus->SetParameter(1,parameter[1] + fReco->GetParError(1));
  fReco_EPlus->SetParameter(2,parameter[2] + fReco->GetParError(2));
  fReco_EPlus->SetParameter(3,parameter[3] + fReco->GetParError(3));

  Int_t Pi0Amplitude_EPlus_fit = fReco_EPlus->GetMaximum(Pi0FitRange[0],Pi0FitRange[1]);
  Int_t Pi0Amplitude_EPlus_bin = fReco_EPlus->GetXaxis()->FindBin(Pi0Amplitude_EPlus_fit);
  Double_t MassFWHM1_EPlus = fReco_EPlus->GetX(Pi0Amplitude_EPlus_fit*0.5,0.1,0.135);
  Double_t MassFWHM2_EPlus = fReco_EPlus->GetX(Pi0Amplitude_EPlus_fit*0.5,0.135,0.16);		
  Double_t FWHM_EPlus = MassFWHM2_EPlus - MassFWHM1_EPlus;
    
  // - error
  TF1 *fReco_EMinus = new TF1("Exp+Gauss","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",Pi0FitRange[0],Pi0FitRange[1]); 
  fReco_EMinus->SetParameter(0,parameter[0] - fReco->GetParError(0));
  fReco_EMinus->SetParameter(1,parameter[1] - fReco->GetParError(1));
  fReco_EMinus->SetParameter(2,parameter[2] - fReco->GetParError(2));
  fReco_EMinus->SetParameter(3,parameter[3] - fReco->GetParError(3));
	    
  Int_t Pi0Amplitude_EMinus_fit = fReco_EMinus->GetMaximum(Pi0FitRange[0],Pi0FitRange[1]);
  Int_t Pi0Amplitude_EMinus_bin = fReco_EMinus->GetXaxis()->FindBin(Pi0Amplitude_EMinus_fit);
  Double_t MassFWHM1_EMinus = fReco_EMinus->GetX(Pi0Amplitude_EMinus_fit*0.5,0.1,0.135);
  Double_t MassFWHM2_EMinus = fReco_EMinus->GetX(Pi0Amplitude_EMinus_fit*0.5,0.135,0.16);		
  Double_t FWHM_EMinus = MassFWHM2_EMinus - MassFWHM1_EMinus;		 
	    
  Double_t Error1 = TMath::Abs(FWHM-FWHM_EPlus);
  Double_t Error2 = TMath::Abs(FWHM-FWHM_EMinus);
	    
  Double_t FWHM_error;
  if(Error1>=Error2) FWHM_error = Error1;
  if(Error1 < Error2) FWHM_error = Error2;
	    
  fstream file;
  file.open(Form("%s%s-%s.dat",outputDir,outputDatFileBase,inputRootFile), ios::out|ios::app);
  file << inputRootFile << "	" << cutSelection << "    " << nGoodEvents<< "	" << Reco << "	" << Reco_error  << "	" << Back << "	" << Back_error << "	" << Mass << "	" << Mass_error << "	" << FWHM << "	" << FWHM_error << "    " << nPi0_MC<< "    " << nPi0_MC_error << endl;
  file.close();

  if(writePlotFiles){
    c_mass->Print(Form("%sInvMass-%s-%s.%s",outputDir,inputRootFile,cutSelection,suffix));
  }
  c_mass->Update();

  delete fReco;
  delete fReco_EPlus;
  delete fReco_EMinus;
  delete c_mass;
}
