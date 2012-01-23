
#include <fstream>
#include <Riostream.h>
/*
 *
 *
 *
 */

extern TRandom *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem *gSystem;

void Extract_Pi0_Characteristics(const char *cutSelection="", const char *inputRootFile = "AnalysisResults",const char *path = "./",const char*outputDir = "./Output/",Bool_t makeMappingPlots=kTRUE){	
	
  gROOT->Reset();	
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(0);
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
	
  // pt bins
  Int_t NBinsPt = 31;
  Double_t BinsPt[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4., 5., 6.,8.,10.,13.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.};  
	
  // input file
  TString filename = Form("%s%s.root",path,inputRootFile);

  TFile f(filename.Data());
  
  TFile *outputFile;
  if(makeMappingPlots==kTRUE){
    outputFile = new TFile(Form("%sPi0Characteristics.root",outputDir),"UPDATE");
  }
  TDirectory *pwg4dir =(TDirectory*)f.Get(Form("PWG4_GammaConversion_%s",cutSelection));
  TList *fHistosGammaConversion = (TList*)pwg4dir->Get(Form("histogramsAliGammaConversion_%s",cutSelection));
  TList *fESDContainer = (TList*)fHistosGammaConversion->FindObject("ESD histograms");
  TList *fMCContainer = (TList*)fHistosGammaConversion->FindObject("MC histograms");
  TList *fBackgroundContainer = (TList*)fHistosGammaConversion->FindObject("Background histograms");

  TH1F *ESD_Mother_InvMass =fESDContainer->FindObject("ESD_Mother_InvMass");
  ESD_Mother_InvMass->SetName(Form("ESD_Mother_InvMass_%s",cutSelection));
  ESD_Mother_InvMass->Write();

  TH1F *ESD_Background_InvMass =fBackgroundContainer->FindObject("ESD_Background_InvMass");
  ESD_Background_InvMass->SetName(Form("ESD_Background_InvMass_%s",cutSelection));
  ESD_Background_InvMass->Write();

  TH1F *ESD_dedx =fESDContainer->FindObject("ESD_ConvGamma_E_dEdxP");
  ESD_dedx->SetName(Form("ESD_ConvGamma_E_dEdxP_%s",cutSelection));
  ESD_dedx->Write();

  TH2F *ESD_alfaqt = fESDContainer->FindObject("ESD_ConvGamma_alfa_qt");
  ESD_alfaqt->SetName(Form("ESD_ConvGamma_alfa_qt_%s",cutSelection));
  ESD_alfaqt->Write();


  TH2F *ESD_Mother_InvMass_vs_Pt = fESDContainer->FindObject("ESD_Mother_InvMass_vs_Pt");
  ESD_Mother_InvMass_vs_Pt->Sumw2();	
  TH2F *ESD_Background_InvMass_vs_Pt = fBackgroundContainer->FindObject("ESD_Background_InvMass_vs_Pt");
  ESD_Background_InvMass_vs_Pt->Sumw2();
		
  //histos for comparison
  TH1F *histoNormYield_Pi0 = new TH1F(Form("Norm_Yield_Pi0_%s",cutSelection),"",NBinsPt,BinsPt);
  TH1F *histoRawYield_Pi0 = new TH1F(Form("Raw_Yield_Pi0_%s",cutSelection),"",NBinsPt,BinsPt);
  TH1F *histoMass_Pi0 = new TH1F(Form("Mass_Pi0_%s",cutSelection),"",NBinsPt,BinsPt);
  TH1F *histoFWHM_Pi0 = new TH1F(Form("FWHM_Pi0_%s",cutSelection),"",NBinsPt,BinsPt);
  TH1F *histoSignificance_Pi0 = new TH1F(Form("Significance_Pi0_%s",cutSelection),"",NBinsPt,BinsPt);
  TH1F *histoSB_Pi0 = new TH1F(Form("SB_Pi0_%s",cutSelection),"",NBinsPt,BinsPt);	    
  TH1F *histoSBRatio_Pi0 = new TH1F(Form("SBRatio_Pi0_%s",cutSelection),"",NBinsPt,BinsPt);	    

  // for yield correction in pt bins
  TH1F *deltaPt = new TH1F("deltaPt","",NBinsPt,BinsPt);
	    
  TH1F * ESD_NumberOfContributorsVtx = fESDContainer->FindObject("ESD_NumberOfContributorsVtx");
  ESD_NumberOfContributorsVtx->SetAxisRange(1.,100.);
  Float_t nGoodEvents = ESD_NumberOfContributorsVtx->Integral();
  cout <<"# good events: "<< nGoodEvents << endl;
		
		
  // binning / max pt from analysis
  Double_t maxPt = ESD_Mother_InvMass_vs_Pt->GetYaxis()->GetBinCenter(ESD_Mother_InvMass_vs_Pt->GetNbinsY()) + ESD_Mother_InvMass_vs_Pt->GetYaxis()->GetBinWidth(ESD_Mother_InvMass_vs_Pt->GetNbinsY())/2.;
  Double_t binningPt = ESD_Mother_InvMass_vs_Pt->GetNbinsY();
		
		
  // ---------------------------------------------------------------------------------------------
	    
  // extract yield, mass, significance and FWHM in pt bins

  TH1D* Mapping_Reco_InvMass_PtBin[32];    // aray of histos for pt slices
  TH1D* Mapping_Back_InvMass_PtBin[32];
  TH1D* Mapping_Signal_InvMass_PtBin[32];   
  TH1D* Mapping_SBRatio_InvMass_PtBin[32];   

  TCanvas *canvas = new TCanvas("canvas","",200,10,600,600);  // gives the page size
  canvas->SetFillColor(0);
		
  for(Int_t iPt = 2; iPt < NBinsPt+1; iPt++){
    TString histonameReco = Form("Mapping_Reco_InvMass_in_Pt_Bin%s%02d",cutSelection ,iPt);
    TString histonameBack = Form("Mapping_Back_InvMass_in_Pt_Bin%s%02d",cutSelection, iPt);

    Mapping_Reco_InvMass_PtBin[iPt]=new TH1D(histonameReco.Data(),histonameReco.Data(),ESD_Mother_InvMass_vs_Pt->GetNbinsX(),0.,1.);
    Mapping_Back_InvMass_PtBin[iPt]=new TH1D(histonameBack.Data(),histonameBack.Data(),ESD_Background_InvMass_vs_Pt->GetNbinsX(),0.,1.);
			
    Int_t startBin = ESD_Mother_InvMass_vs_Pt->GetYaxis()->FindBin(BinsPt[iPt-1])-1;
    Int_t endBin = ESD_Mother_InvMass_vs_Pt->GetYaxis()->FindBin(BinsPt[iPt])-1;
			
    Double_t startpt = startBin*maxPt/binningPt;	
    Double_t endpt = endBin*maxPt/binningPt;	
			
    ESD_Mother_InvMass_vs_Pt->ProjectionX(histonameReco.Data(),startBin,endBin);
    ESD_Background_InvMass_vs_Pt->ProjectionX(histonameBack.Data(),startBin,endBin);
			
    Mapping_Reco_InvMass_PtBin[iPt]=(TH1D*)gDirectory->Get(histonameReco.Data());
    Mapping_Reco_InvMass_PtBin[iPt]->Rebin(RebinMass);
    Mapping_Back_InvMass_PtBin[iPt]=(TH1D*)gDirectory->Get(histonameBack.Data());
    Mapping_Back_InvMass_PtBin[iPt]->Rebin(RebinMass);
			
			
    // normalisation of background
    TAxis *xaxis_reco = Mapping_Reco_InvMass_PtBin[iPt]->GetXaxis();							
    Int_t r_1 = xaxis_reco->FindBin(BGFit_range[0]);
    Int_t r_2 = xaxis_reco->FindBin(BGFit_range[1]);  
    Double_t r =  Mapping_Reco_InvMass_PtBin[iPt]->Integral(r_1,r_2); // Integral(75,125)
    TAxis *xaxis_back = Mapping_Back_InvMass_PtBin[iPt]->GetXaxis();							
    Int_t b_1 = xaxis_back->FindBin(BGFit_range[0]);
    Int_t b_2 = xaxis_back->FindBin(BGFit_range[1]);   
    Double_t b =  Mapping_Back_InvMass_PtBin[iPt]->Integral(b_1,b_2);
    Double_t norm = 1;
    if(b != 0) norm = r/b;
    Mapping_Back_InvMass_PtBin[iPt]->Sumw2();		
    Mapping_Back_InvMass_PtBin[iPt]->Scale(norm);
	      
    TString histonameSignal = Form("Mapping_Signal_InvMass_in_Pt_Bin%s%02d",cutSelection, iPt);
    Mapping_Signal_InvMass_PtBin[iPt] = (TH1D*)Mapping_Reco_InvMass_PtBin[iPt]->Clone();
    Mapping_Signal_InvMass_PtBin[iPt]->SetName(histonameSignal.Data());
    Mapping_Signal_InvMass_PtBin[iPt]->SetTitle(histonameSignal.Data());
    Mapping_Signal_InvMass_PtBin[iPt]->Sumw2();	
    Mapping_Signal_InvMass_PtBin[iPt]->Add(Mapping_Back_InvMass_PtBin[iPt],-1.);

    /*    TString histonameSBRatio= Form("Mapping_SBRatio_InvMass_in_Pt_Bin%02d", iPt);
    Mapping_SBRatio_InvMass_PtBin[iPt]=->ProjectionX(histonameSBRatio.Data(),startBin,endBin);;
    Mapping_SBRatio_InvMass_PtBin[iPt]->Add(Mapping_Reco_InvMass_PtBin[iPt],1);
    Mapping_SBRatio_InvMass_PtBin[iPt]->Divide(Mapping_Back_InvMass_PtBin[iPt]);
    Mapping_SBRatio_InvMass_PtBin[iPt]->Write();
    */
    Int_t bins =  Mapping_Signal_InvMass_PtBin[iPt]->GetNbinsX();
    for(Int_t ib = 1; ib < bins; ib++){
      Double_t error = Mapping_Signal_InvMass_PtBin[iPt]->GetBinError(ib);
      if(abs(error) < 0.1) error = Mapping_Signal_InvMass_PtBin[iPt]->GetBinError(ib) +1;
      Mapping_Signal_InvMass_PtBin[iPt]->SetBinError(ib,error);
    }
			
    Double_t Pi0Amplitude = Mapping_Signal_InvMass_PtBin[iPt]->GetMaximum();
    Double_t Pi0Amplitude_Min = Pi0Amplitude*80/100;
    Double_t Pi0Amplitude_Max = Pi0Amplitude*100/100;
	      
    TAxis *xaxis = Mapping_Signal_InvMass_PtBin[iPt]->GetXaxis();	
    Int_t ilowmassPi0_bin = xaxis->FindBin(Pi0FitRange[0]);  
    Int_t ihighmassPi0_bin =xaxis->FindBin(Pi0FitRange[1]);	
    if(makeMappingPlots){
      Mapping_Signal_InvMass_PtBin[iPt]->Write();
    }
    TF1 *fReco = new TF1(Form("Exp+Gauss PtBin[%02d]",iPt),"(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",Pi0FitRange[0],Pi0FitRange[1]); 
	      
    // fit the peak around pi0 mass and draw		
    fReco->SetParameter(0,Pi0Amplitude);
    fReco->SetParameter(1,Pi0Mass);
    fReco->SetParameter(2,Pi0Width);
    fReco->SetParameter(3,Pi0Slope);
	      
    fReco->SetParLimits(0,Pi0Amplitude_Min,Pi0Amplitude_Max);
    fReco->SetParLimits(1,Pi0Mass_range[0],Pi0Mass_range[1]);
    fReco->SetParLimits(2,Pi0Width_range[0],Pi0Width_range[1]);
    fReco->SetParLimits(3,Pi0Slope_range[0],Pi0Slope_range[1]);
	      
    Mapping_Signal_InvMass_PtBin[iPt]->Fit(fReco,"RME");
    //    Mapping_Signal_InvMass_PtBin[iPt]->Draw();
    fReco->SetLineColor(3);
    fReco->SetLineWidth(0.7);
    fReco->SetLineStyle(4);		
    fReco->Draw("same");
    if(makeMappingPlots){
      fReco->Write();
    }
    cout<<"Kenneth: Setting Bin content for bin: "<<iPt<<"  to "<<fReco->GetParameter(1)<<endl;
    histoMass_Pi0->SetBinContent(iPt, fReco->GetParameter(1));
    cout<<"Kenneth: Setting Bin content for bin: "<<iPt<<"  to "<<fReco->GetParError(1)<<endl;
    histoMass_Pi0->SetBinError(iPt, fReco->GetParError(1));

    TAxis *xaxis = Mapping_Signal_InvMass_PtBin[iPt]->GetXaxis();	
    Int_t ilowmassPi0_bin = xaxis->FindBin(Pi0FitRange[0]);  
    Int_t ihighmassPi0_bin =xaxis->FindBin(Pi0FitRange[1]);

    Double_t Reco_error;
    Double_t Background_error;
    
    Double_t IntReco = Mapping_Reco_InvMass_PtBin[iPt]->IntegralAndError(ilowmassPi0_bin,ihighmassPi0_bin,Reco_error);
    Double_t IntBackground = Mapping_Back_InvMass_PtBin[iPt]->IntegralAndError(ilowmassPi0_bin,ihighmassPi0_bin,Background_error);
    if(IntBackground>0){
      Double_t SB= IntReco/IntBackground;
      Double_t SB_error = sqrt( pow(Reco_error/IntBackground,2) + pow(IntReco/IntBackground/IntBackground*Background_error,2));		
      
      histoSB_Pi0->SetBinContent(iPt,SB);
      histoSB_Pi0->SetBinError(iPt, SB_error);
    }
    //    canvas->Update();
    //    canvas->Clear();
    
    Double_t parameter[4];
    fReco->GetParameters(parameter);
			
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
	      
    histoFWHM_Pi0->SetBinContent(iPt, FWHM);
    histoFWHM_Pi0->SetBinError(iPt,FWHM_error);
			
    Double_t Reco_error;
    Double_t Background_error;
	      
    Double_t Reco = Mapping_Reco_InvMass_PtBin[iPt]->IntegralAndError(ilowmassPi0_bin,ihighmassPi0_bin,Reco_error);
    Double_t Background = Mapping_Back_InvMass_PtBin[iPt]->IntegralAndError(ilowmassPi0_bin,ihighmassPi0_bin,Background_error);
	      
    if(Background>0){
      Double_t Significance = Reco/sqrt(Background);
      Double_t Significance_error = sqrt( pow(Reco_error/sqrt(Background),2) + pow( -0.5 * Reco/pow(Background,1.5) * Background_error,2) );
      
      histoSignificance_Pi0->SetBinContent(iPt,Significance);
      histoSignificance_Pi0->SetBinError(iPt, Significance_error);
    }
    histoNormYield_Pi0->SetBinContent(iPt,(Reco-Background)/nGoodEvents);
    histoNormYield_Pi0->SetBinError(iPt,sqrt(Reco_error*Reco_error + Background_error*Background_error)/nGoodEvents);

    histoRawYield_Pi0->SetBinContent(iPt,Reco-Background);
    histoRawYield_Pi0->SetBinError(iPt,sqrt(Reco_error*Reco_error + Background_error*Background_error));
	      
	      
    deltaPt->SetBinContent(iPt, endpt-startpt);
    deltaPt->SetBinError(iPt,0);
    if(iPt == 2){
      deltaPt->SetBinContent(1,startpt-0);
      deltaPt->SetBinError(1,0);
    } 
  }
		
  delete [] Mapping_Reco_InvMass_PtBin;
  delete [] Mapping_Back_InvMass_PtBin;
  delete [] Mapping_Signal_InvMass_PtBin;   
  delete [] Mapping_SBRatio_InvMass_PtBin;
  // end of extract yield, mass .....
		
  // write the created histos into the output root file	  

  if(makeMappingPlots == kFALSE){// open the root file here so the mapping histograms is not included
    outputFile = new TFile(Form("%sPi0Characteristics.root",outputDir),"UPDATE");
  }

  
  //    ESD_Mother_InvMass_vs_Pt->ProjectionX(histonameReco.Data(),startBin,endBin);
  //    ESD_Background_InvMass_vs_Pt->ProjectionX(histonameBack.Data(),startBin,endBin);
  /*
  histoSBRatio_Pi0->Add(ESD_Mother_InvMass_vs_Pt,1);
  histoSBRatio_Pi0->Divide(ESD_Background_InvMass_vs_Pt);
  histoSBRatio_Pi0->Write();
  */
  histoNormYield_Pi0->SetXTitle("p_{t} (GeV/c)");
  histoNormYield_Pi0->Divide(deltaPt);
  histoNormYield_Pi0->SetYTitle("#pi^{0} raw yield/per event");  
  histoNormYield_Pi0->DrawCopy("e1");  	
  histoNormYield_Pi0->Write();
  //  outputFile->Write();
  canvas->Update();
  canvas->Clear();

  histoRawYield_Pi0->SetXTitle("p_{t} (GeV/c)");
  histoRawYield_Pi0->Divide(deltaPt);
  histoRawYield_Pi0->SetYTitle("#pi^{0} raw yield");
  histoRawYield_Pi0->DrawCopy("e1");  	
  histoRawYield_Pi0->Write();
  //  outputFile->Write();
  canvas->Update();
  canvas->Clear();
	    
  histoMass_Pi0->SetXTitle("p_{t} (GeV/c)");
  //  histoMass_Pi0->Divide(deltaPt);
  histoMass_Pi0->SetYTitle("#pi^{0} mass");  
  histoMass_Pi0->DrawCopy("e1");  	
  histoMass_Pi0->Write();
  canvas->Update();
  canvas->Clear();
	   
  
  histoSB_Pi0->SetXTitle("p_{t} (GeV/c)");
  //  histoMass_Pi0->Divide(deltaPt);
  histoSB_Pi0->SetYTitle("S/B");  
  histoSB_Pi0->DrawCopy("e1");  	
  histoSB_Pi0->Write();
  canvas->Update();
  canvas->Clear(); 
  
 
  histoFWHM_Pi0->SetXTitle("p_{t} (GeV/c)");
  //  histoFWHM_Pi0->Divide(deltaPt);
  histoFWHM_Pi0->SetYTitle("#pi^{0} FWHM/2.36");  
  histoFWHM_Pi0->DrawCopy("e1");  	
  histoFWHM_Pi0->Write();
  canvas->Update();
  canvas->Clear();
	    
  histoSignificance_Pi0->SetXTitle("p_{t} (GeV/c)");
  //  histoSignificance_Pi0->Divide(deltaPt);
  histoSignificance_Pi0->SetYTitle("Significance");  
  histoSignificance_Pi0->DrawCopy("e1");  	
  histoSignificance_Pi0->Write();
  canvas->Update();
  canvas->Clear();
	    
  histoSBRatio_Pi0->Delete();
  histoNormYield_Pi0->Delete();
  histoRawYield_Pi0->Delete();
  histoMass_Pi0->Delete();
  histoFWHM_Pi0->Delete();
  histoSignificance_Pi0->Delete();
  histoSB_Pi0->Delete();
  deltaPt->Delete();
	    
  canvas->Delete();
  
  outputFile->Write();  
  outputFile->Close();  
}
