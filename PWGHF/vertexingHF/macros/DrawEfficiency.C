
// macro to calculate the pt RecAnCuts/MCAcc
// efficiency, from a CF container 
// Author: C. Zampolli 

// channel could be:
// D0: D0 -> Kpi, from old task
// D0_New: D0 -> Kpi, from new CF common implementation
// Dplus_New: D+ -> Kpipi, from new CF common implementation
//
// eff = sum of the efficiency indexes to compute
// MCAcc_over_MCLimAcc = 0x001
// RecPPR_over_MCAcc = 0x002
// RecPID_over_MCAcc = 0x004
// all = 0x007
//
// selection = D origin
// 0 --> from c only
// 1 --> from b only
// 2 --> from both c and b

#include <Riostream.h>

extern TRandom *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem *gSystem;

void DrawEfficiency(const char* channel, Int_t selection = 0, Int_t ieff = 7){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetOptTitle(0);
	
	gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include  -I$ROOTSYS/include");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");
	gSystem->Load("$ALICE_ROOT/CORRFW/libCORRFW") ;
	gSystem->Load("libPWGHFbase");
	gSystem->Load("libPWGHFvertexingHF");
		
	Int_t mcAcc_over_mcLimAcc = 0x001;
	Int_t recPPR_over_mcAcc = 0x002;
	Int_t recPID_over_mcAcc = 0x004;

	// pt index
	Int_t ipt =0;

	Int_t stepNum = -1;
	Int_t stepDen = -1;
	
	// Read the  container from file
	TFile* f = new TFile("AnalysisResults.root");	
	TString directoryName;
	TString containerName;
	TString cutName;
	TString outfileName;

	if (channel == "D0") {
		if (selection == 0){
			directoryName = "PWG3_D2H_CFtaskD0toKpi";
			containerName = "CFHFccontainer0";
			cutName = "Cuts";
			outfileName = "fileEff_D0_from_c.root";
		}
		else if (selection == 1){
			directoryName = "PWG3_D2H_CFtaskD0toKpiKeepD0fromBOnly";
			containerName = "CFHFccontainer0D0fromB";
			cutName = "Cuts";
			outfileName = "fileEff_D0_from_b.root";
		}
		else if (selection == 2){
			directoryName = "PWG3_D2H_CFtaskD0toKpiKeepD0fromB";
			containerName = "CFHFccontainer0allD0";
			cutName = "Cuts";
			outfileName = "fileEff_D0_from_c_and_b.root";
		}
		else {
			Printf("not a valid selection, return");
			return;
		}
	}
	else if (channel == "D0_New"){
		if (selection == 0){
			directoryName = "PWG3_D2H_CFtaskD0toKpi_NEW";
			containerName = "CFHFccontainer0_New";
			cutName = "Cuts_New";
			//directoryName = "PWG3_D2H_CFtaskD0toKpi_CommonFramework";
			//containerName = "CFHFccontainer0_CommonFramework";
			//cutName = "Cuts_CommonFramework";
			outfileName = "fileEff_D0_CommonFramework_from_c.root";
		}
		else if (selection == 1){
			directoryName = "PWG3_D2H_CFtaskD0toKpiKeepDfromBOnly";
			containerName = "CFHFccontainer0DfromB_New";
			cutName = "Cuts_New";
			//directoryName = "PWG3_D2H_CFtaskD0toKpiKeepDfromBOnly_CommonFramework";
			//containerName = "CFHFccontainer0DfromB_CommonFramework";	
			//cutName = "Cuts_CommonFramework";
			outfileName = "fileEff_D0_CommonFramework_from_b.root";
		}
		else if (selection == 2){
			directoryName = "PWG3_D2H_CFtaskD0toKpiKeepDfromB_NEW";
			containerName = "CFHFccontainer0allD_New";
			cutName = "Cuts_New";
			//directoryName = "PWG3_D2H_CFtaskD0toKpiKeepDfromB_CommonFramework";
			//containerName = "CFHFccontainer0allD_CommonFramework";
			//cutName = "Cuts_CommonFramework";
			outfileName = "fileEff_D0_CommonFramework_from_c_and_b.root";
		}
		else {
			Printf("not a valid selection, return");
			return;
		}
	}
	else if (channel == "Dplus_New"){
		if (selection == 0){
			directoryName = "PWG3_D2H_CFtaskDplustoKpipi_NEW";
			containerName = "CFHFccontainer0_New_3Prong";
			cutName = "Cuts_3Prong";
			//directoryName = "PWG3_D2H_CFtaskDplustoKpipi_CommonFramework";
			//containerName = "CFHFccontainer0_3Prong_CommonFramework";
			//cutName =  "Cuts_3Prong_CommonFramework";
			outfileName = "fileEff_Dplus_CommonFramework_from_c.root";
		}
		else if (selection == 1){
			directoryName = "PWG3_D2H_CFtaskDplustoKpipiKeepDfromBOnly";
			containerName = "CFHFccontainer0DfromB_New_3Prong";
			cutName = "Cuts_3Prong";
			//directoryName = "PWG3_D2H_CFtaskDplustoKpipiKeepDfromBOnly_CommonFramework";
			//containerName = "CFHFccontainer0DfromB_3Prong_CommonFramework";
			//cutName =  "Cuts_3Prong_CommonFramework";
			outfileName = "fileEff_Dplus_CommonFramework_from_b.root";
		}
		else if (selection == 2){
			directoryName = "PWG3_D2H_CFtaskDplustoKpipiKeepDfromB_NEW";
			containerName = "CFHFccontainer0allD_New_3Prong";
			cutName = "Cuts_3Prong";
			//directoryName = "PWG3_D2H_CFtaskDplustoKpipiKeepDfromB_CommonFramework";
			//containerName = "CFHFccontainer0allD_3Prong_CommonFramework";
			//cutName =  "Cuts_3Prong_CommonFramework";
			outfileName = "fileEff_Dplus_CommonFramework_from_c_and_b.root";
		}
		else {
			Printf("not a valid selection, return");
			return;
		}
	}
	else {
		Printf("not a valid channel, return");
		return;
	}

	Printf("Opening file Analysisresults.root");
	Printf("Reading Directory \"%s\"",directoryName.Data());
	Printf("Getting CF Container \"%s\"",containerName.Data());
	Printf("Getting Cut Object \"%s\"",cutName.Data());


	TDirectoryFile* d = (TDirectoryFile*)f->Get(directoryName.Data());
	if (!d){
		Printf("Directory does not exist! Check directory name (%s) in file AnalysisResults.root - returning...", directoryName.Data());
		return;
	}
	AliCFContainer *data = (AliCFContainer*) (d->Get(containerName.Data()));
	AliRDHFCuts *cutsRDHF = (AliRDHFCuts*)(d->Get(cutName.Data()));

	if (!data){
		Printf("Container does not exist! Check container name (%s) in directory %s - returning...", containerName.Data(), directoryName.Data());
		return;
	}
	
	TFile* fileEff = new TFile(outfileName.Data(), "RECREATE");
	TString plotDir(Form("EffPlots/%s",channel));
	gSystem->Exec(Form("mkdir EffPlots"));
	gSystem->Exec(Form("mkdir %s",plotDir.Data()));
	
	//construct the efficiency grid from the data container 
	AliCFEffGrid *eff = new AliCFEffGrid("eff"," The efficiency",*data);

	TCanvas *ceffpt = new TCanvas("ceffpt","Efficiency vs pt",50,50,550,550);
	ceffpt->cd();
	ceffpt->SetLeftMargin(0.15);
	ceffpt->SetRightMargin(0.05);
	TH1D *hpteffCF = 0x0; //the efficiency vs pt

	if (ieff & mcAcc_over_mcLimAcc){
		AliCFEffGrid *eff = new AliCFEffGrid("eff"," The efficiency",*data);
		stepDen = (Int_t)(AliCFTaskVertexingHF::kStepGeneratedLimAcc);	
		stepNum = (Int_t)(AliCFTaskVertexingHF::kStepAcceptance);	
		printf("Calculating efficiency for mcAcc_over_mcLimAcc: stepDen = %d, stepNum = %d\n",stepDen,stepNum);	
		eff->CalculateEfficiency(stepNum,stepDen); //eff= step1/step0
		
		//canvas
		ceffpt->cd();

		//The efficiency along the variables
		hpteffCF = eff->Project(ipt); 
		SetHistoEff(hpteffCF,8,20,"mcAcc_over_mcLimAcc");
		hpteffCF->Draw("hist");
		hpteffCF->Draw("err same");
		fileEff->cd();
		hpteffCF->Write("hpteffCF_mcAcc_over_mcLimAcc");
		
		// printing png files
		ceffpt->Print(Form("%s/effpt_mcAcc_over_mcLimAcc.png", plotDir.Data()));
		ceffpt->Print(Form("%s/effpt_mcAcc_over_mcLimAcc.eps", plotDir.Data()));
		ceffpt->Print(Form("%s/effpt_mcAcc_over_mcLimAcc.gif", plotDir.Data()));
		delete eff;
		eff = 0x0;
	}

	if (ieff & recPPR_over_mcAcc){
		AliCFEffGrid *eff = new AliCFEffGrid("eff"," The efficiency",*data);
		stepDen = (Int_t)(AliCFTaskVertexingHF::kStepAcceptance);	
		stepNum = (Int_t)(AliCFTaskVertexingHF::kStepRecoPPR);	
		printf("Calculating efficiency for RecPPR_over_mcAcc: stepDen = %d, stepNum = %d\n",stepDen,stepNum);	
		eff->CalculateEfficiency(stepNum,stepDen); //eff= step1/step0
		
		//canvas
		ceffpt->cd();
		
		//The efficiency along the variables
		hpteffCF = eff->Project(ipt); 
		SetHistoEff(hpteffCF,8,20, "recAnCuts_over_mcAcc");
		hpteffCF->Draw("hist");
		hpteffCF->Draw("err same");
		fileEff->cd();
		hpteffCF->Write("hpteffCF_RecAnCut_over_mcAcc");
		
		// printing png files
		ceffpt->Print(Form("%s/effpt_RecAnCut_over_mcAcc.png", plotDir.Data()));
		ceffpt->Print(Form("%s/effpt_RecAnCut_over_mcAcc.eps", plotDir.Data()));
		ceffpt->Print(Form("%s/effpt_RecAnCut_over_mcAcc.gif", plotDir.Data()));
		delete eff;
		eff = 0x0;
	}

	if (ieff & recPID_over_mcAcc){
		AliCFEffGrid *eff = new AliCFEffGrid("eff"," The efficiency",*data);
		stepDen = (Int_t)(AliCFTaskVertexingHF::kStepAcceptance);	
		stepNum = (Int_t)(AliCFTaskVertexingHF::kStepRecoPID);	
		printf("Calculating efficiency for RecPID_over_mcAcc: stepDen = %d, stepNum = %d\n",stepDen,stepNum);	
		eff->CalculateEfficiency(stepNum,stepDen); //eff= step1/step0
		
		//canvas
		ceffpt->cd();
		
		//The efficiency along the variables
		hpteffCF = eff->Project(ipt); 
		SetHistoEff(hpteffCF,8,20,"recPID_over_mcAcc");
		hpteffCF->Draw("hist");
		hpteffCF->Draw("err same");
		fileEff->cd();
		hpteffCF->Write("hpteffCF_RecPID_over_mcAcc");
		
		// printing png files
		ceffpt->Print(Form("%s/effpt_RecPID_over_mcAcc.png", plotDir.Data()));
		ceffpt->Print(Form("%s/effpt_RecPID_over_mcAcc.eps", plotDir.Data()));
		ceffpt->Print(Form("%s/effpt_RecPID_over_mcAcc.gif", plotDir.Data()));
		delete eff;
		eff = 0x0;
	}
	
	cutsRDHF->Write("Cuts");

	// writing single distributions
	TH1D *hMCAccpt = data->ShowProjection(ipt, AliCFTaskVertexingHF::kStepAcceptance);
	SetHistoDistribution(hMCAccpt, 1, 20);
	hMCAccpt->Draw();
	TH1D *hMCLimAccpt = data->ShowProjection(ipt, AliCFHeavyFlavourTaskMultiVarMultiStep::kStepGeneratedLimAcc);
	SetHistoDistribution(hMCLimAccpt, 4, 20);
	TH1D *hRecoAnCutspt = data->ShowProjection(ipt, AliCFTaskVertexingHF::kStepRecoPPR);
	SetHistoDistribution(hRecoAnCutspt, 8, 20);
	TH1D *hRecoPIDpt = data->ShowProjection(ipt, AliCFTaskVertexingHF::kStepRecoPID);
	SetHistoDistribution(hRecoPIDpt, 6, 20);
	hMCAccpt->Write("hMCAccpt");
	hMCLimAccpt->Write("hMCLimAccpt");
	hRecoAnCutspt->Write("hRecoAnCutspt");
	hRecoPIDpt->Write("hRecoPIDpt");

	//	fileEff->Close(); // commented out to see the canvas on the screen....

}

void SetHistoEff(TH1D* h, Int_t color, Int_t style, const char* effType){

	h->SetLineColor(color);
	h->SetLineWidth(3);
	h->SetMarkerStyle(style);
	h->SetMarkerColor(color);
	h->SetMarkerSize(1.2);
	h->GetYaxis()->SetTitleOffset(1.5);
	h->GetXaxis()->SetTitleOffset(1.2);
	h->GetXaxis()->SetTitle("p_{t} [GeV/c]");
	h->GetYaxis()->SetTitle(Form("%s, Efficiency",effType));
	return;
}
void SetHistoDistribution(TH1D* h, Int_t color, Int_t style){

	h->SetLineColor(color);
	h->SetLineWidth(3);
	h->SetMarkerStyle(style);
	h->SetMarkerColor(color);
	h->SetMarkerSize(1.2);
	h->GetXaxis()->SetTitleOffset(1.2);
	h->GetXaxis()->SetTitle("p_{t} [GeV/c]");
	return;
}
