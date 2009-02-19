#include <Riostream.h>

extern TRandom *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem *gSystem;

void ReadCFHeavyFlavourOutput(){

	// example macro for reading a container coming from a CF analysis
	// the generated and reconstructed distributions are got, 
	// the efficiencies are calculated,
	// the reco distributions are corrected,
	// the results plotted
	
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(1110);
	gStyle->SetPalette(1);
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetOptTitle(0);
	
	gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include  -I$ROOTSYS/include");
	gSystem->Load("libANALYSIS.so");
	gSystem->Load("libANALYSISalice.so");
	gSystem->Load("$ALICE_ROOT/CORRFW/libCORRFW.so") ;
	
	//Setting up the container grid... 
	
	const Int_t nstep=2; //number of selection steps  (just 2 in this ex)
	
	const Int_t nvar=6; //number of variables on the grid: pt, y, cosThetaStar, ptPi, ptK, cT
	
	// Flag the sel steps. In this example, we have two, may be any nstep
	Int_t stepGen=0;
	Int_t stepRec=1;
	
	// the sensitive variables, their indeces
	UInt_t ipt  =0;
	UInt_t iy   =1;
	UInt_t icTS = 2; // cTS stands for cosThetaStar
	UInt_t iptPi= 3;
	UInt_t iptK = 4;
	UInt_t icT  = 5;
	
	// Read the  container from file
	TFile *file = new TFile("output.root");
	AliCFContainer *data = (AliCFContainer*) (file->Get("container"));
	
	// Make some 1 & 2-D projections..
	// MC level
	TCanvas *cmc1 = new TCanvas("cmc1","The MC distributions",0,800,900,1200);
	TCanvas *cmc2 = new TCanvas("cmc2","The MC distributions",0,800,900,1200);
	TCanvas *cmcpt = new TCanvas("cmcpt","pt distribution from MC",50,50,550,550);
	TCanvas *cmcy = new TCanvas("cmcy","y distribution from MC",50,50,550,550);
	TCanvas *cmccTS = new TCanvas("cmcTS","cosThetaStar distribution from MC",50,50,550,550);  // cTS stands for cosThetaStar
	TCanvas *cmcptPi = new TCanvas("cmcptPi","pt_pi distribution from MC",50,50,550,550);
	TCanvas *cmcptK = new TCanvas("cmcptK","pt_K distribution from MC",50,50,550,550);
	TCanvas *cmccT = new TCanvas("cmccT","cT distribution from MC",50,50,550,550);

	// Reco-aod level
	TCanvas *cpt = new TCanvas("cpt","pt distribution from reco aod",50,50,550,550);
	TCanvas *cy = new TCanvas("cy","y distribution from reco aod",50,50,550,550);
	TCanvas *ccTS = new TCanvas("cTS","cosThetaStar distribution from reco aod",50,50,550,550); 
	TCanvas *cptPi = new TCanvas("cptPi","pt_pi distribution from reco aod",50,50,550,550);
	TCanvas *cptK = new TCanvas("cptK","pt_K distribution from reco aod",50,50,550,550);
	TCanvas *ccT = new TCanvas("ccT","cT distribution from reco aod",50,50,550,550);

	// some 2D distributions in pt and y
	TCanvas *cGen2d = new TCanvas("cGen2d","2D distribution from MC ",50,50,550,550);
	TCanvas *cRec2d = new TCanvas("cRec2d","2D distribution from reco aod",50,50,550,550);
	
	// 2D plots 
	TH2D *hMCpty2d = data->ShowProjection(ipt, iy, stepGen);
	hMCpty2d->SetLineColor(2);
	hMCpty2d->SetLineWidth(3);
	hMCpty2d->SetMarkerColor(2);
	hMCpty2d->SetMarkerStyle(20);
	hMCpty2d->GetXaxis()->SetTitleOffset(1.2);
	hMCpty2d->GetYaxis()->SetTitleOffset(1.5);
	hMCpty2d->GetXaxis()->SetTitle("p_{T} (GeV/c), MC data");
	hMCpty2d->GetYaxis()->SetTitle("y, MC data");
	hMCpty2d->Draw("text");
	cGen2d->cd();
	cGen2d->SetLeftMargin(0.15);
	cGen2d->SetRightMargin(0.05);
	cGen2d->Update();  

	TH2D *hRECpty2d = data->ShowProjection(ipt, iy, stepRec);
	hRECpty2d->Sumw2();
	hRECpty2d->SetLineColor(4);
	hRECpty2d->SetLineWidth(3);
	hRECpty2d->SetMarkerColor(4);
	hRECpty2d->SetMarkerStyle(20);
	hRECpty2d->GetXaxis()->SetTitleOffset(1.2);
	hRECpty2d->GetYaxis()->SetTitleOffset(1.5);
	hRECpty2d->GetXaxis()->SetTitle("p_{T} (GeV/c), AOD");
	hRECpty2d->GetYaxis()->SetTitle("y, AOD");
	hRECpty2d->Draw("text");
	cRec2d->cd();
	cRec2d->SetLeftMargin(0.15);
	cRec2d->SetRightMargin(0.05);
	cRec2d->Update();  

	// MC + REC 1D plots 
	// pt, y, cosThetaStar
	cmc1->Divide(2,3);   

	cmc1->cd(1);
	TH1D *hMCpt1D = data->ShowProjection(ipt, stepGen);
	Double_t maxpt = hMCpt1D->GetMaximum();
	hMCpt1D->GetYaxis()->SetRangeUser(0,maxpt*1.2);
	hMCpt1D->Sumw2();
	hMCpt1D->SetMinimum(0.01);
	hMCpt1D->SetLineColor(2);
	hMCpt1D->SetLineWidth(3);
	hMCpt1D->SetMarkerColor(2);
	hMCpt1D->SetFillColor(2);
	hMCpt1D->SetFillStyle(3005);
	hMCpt1D->SetMarkerStyle(20);
	hMCpt1D->GetXaxis()->SetTitleOffset(1.2);
	hMCpt1D->GetXaxis()->SetTitle("p_{T} (GeV/c), MC data");
	hMCpt1D->Draw("hist");
	cmcpt->cd();
	cmcpt->SetLeftMargin(0.15);
	cmcpt->SetRightMargin(0.05);
	hMCpt1D->Draw("hist");
	hMCpt1D->Draw("err same");
	cmcpt->Update();

	cmc1->cd(2);
	TH1D *hRECpt1D = data->ShowProjection(ipt, stepRec);
	hRECpt1D->GetYaxis()->SetRangeUser(0,maxpt*1.2);
	hRECpt1D->SetLineColor(4);
	hRECpt1D->SetLineWidth(3);
	hRECpt1D->SetMarkerColor(4);
	hRECpt1D->SetFillColor(4);
	hRECpt1D->SetFillStyle(3004);
	hRECpt1D->SetMarkerStyle(20);
	hRECpt1D->GetXaxis()->SetTitleOffset(1.2);
	hRECpt1D->Sumw2();
	hRECpt1D->SetMinimum(0.01);
	hRECpt1D->GetXaxis()->SetTitle("p_{T} (GeV/c), AOD");
	hRECpt1D->Draw("hist");
	cpt->cd();
	cpt->SetLeftMargin(0.15);
	cpt->SetRightMargin(0.05);
	hRECpt1D->Draw("hist");
	hRECpt1D->Draw("err same");
	cpt->Update();

	cmc1->cd(3);
	TH1D *hMCy1D = data->ShowProjection(iy, stepGen);
	Double_t maxy = hMCy1D->GetMaximum();
	hMCy1D->GetYaxis()->SetRangeUser(0,maxy*1.2);
	hMCy1D->SetLineColor(2);
	hMCy1D->SetLineWidth(3);
	hMCy1D->SetMarkerColor(2);
	hMCy1D->SetFillColor(2);
	hMCy1D->SetFillStyle(3005);
	hMCy1D->SetMarkerStyle(20);
	hMCy1D->GetXaxis()->SetTitleOffset(1.2);
	hMCy1D->Sumw2();
	hMCy1D->SetMinimum(0.01);
	hMCy1D->GetXaxis()->SetTitle("y, MC data");
	hMCy1D->Draw("hist");
	cmcy->cd();
	cmcy->SetLeftMargin(0.15);
	cmcy->SetRightMargin(0.05);
	hMCy1D->Draw("hist");
	hMCy1D->Draw("err same");
	cmcy->Update();

	cmc1->cd(4);
	TH1D *hRECy1D = data->ShowProjection(iy, stepRec);
	hRECy1D->GetYaxis()->SetRangeUser(0,maxy*1.2);
	hRECy1D->SetLineColor(4);
	hRECy1D->SetLineWidth(3);
	hRECy1D->SetMarkerColor(4);
	hRECy1D->SetFillColor(4);
	hRECy1D->SetFillStyle(3004);
	hRECy1D->SetMarkerStyle(20);
	hRECy1D->GetXaxis()->SetTitleOffset(1.2);
	hRECy1D->Sumw2();
	hRECy1D->SetMinimum(0.01);
	hRECy1D->Draw("hist");
	cy->cd();
	cy->SetLeftMargin(0.15);
	cy->SetRightMargin(0.05);
	hRECy1D->GetXaxis()->SetTitle("y, AOD");
	hRECy1D->Draw("hist");
	hRECy1D->Draw("err same");
	cy->Update();

	cmc1->cd(5);
	TH1D *hMCcTS1D = data->ShowProjection(icTS, stepGen);
	Double_t maxcTS = hMCcTS1D->GetMaximum();
	hMCcTS1D->GetYaxis()->SetRangeUser(0,maxcTS*1.2);
	hMCcTS1D->SetLineColor(2);
	hMCcTS1D->SetLineWidth(3);
	hMCcTS1D->SetMarkerColor(2);
	hMCcTS1D->SetFillColor(2);
	hMCcTS1D->SetFillStyle(3005);
	hMCcTS1D->SetMarkerStyle(20);
	hMCcTS1D->GetXaxis()->SetTitleOffset(1.2);
	hMCcTS1D->Sumw2();
	hMCcTS1D->SetMinimum(0.01);
	hMCcTS1D->GetXaxis()->SetTitle("cosThetaStar, MC data");
	hMCcTS1D->Draw("hist");
	cmccTS->cd();
	cmccTS->SetLeftMargin(0.15);
	cmccTS->SetRightMargin(0.05);
	hMCcTS1D->Draw("hist");
	hMCcTS1D->Draw("err same");
	cmccTS->Update();

	cmc1->cd(6);
	TH1D *hRECcTS1D = data->ShowProjection(icTS, stepRec);
	hRECcTS1D->GetYaxis()->SetRangeUser(0,maxcTS*1.2);
	hRECcTS1D->SetLineColor(4);
	hRECcTS1D->SetLineWidth(3);
	hRECcTS1D->SetMarkerColor(4);
	hRECcTS1D->SetFillColor(4);
	hRECcTS1D->SetFillStyle(3004);
	hRECcTS1D->SetMarkerStyle(20);
	hRECcTS1D->GetXaxis()->SetTitleOffset(1.2);
	hRECcTS1D->Sumw2();
	hRECcTS1D->SetMinimum(0.01);
	hRECcTS1D->Draw("hist");
	ccTS->cd();
	ccTS->SetLeftMargin(0.15);
	ccTS->SetRightMargin(0.05);
	hRECcTS1D->GetXaxis()->SetTitle("cosThetaStar, AOD");
	hRECcTS1D->Draw("hist");
	hRECcTS1D->Draw("err same");
	ccTS->Update();

	// ptPi, ptK, cT
	cmc2->Divide(2,3); 

	cmc2->cd(1);
	TH1D *hMCptPi1D = data->ShowProjection(iptPi, stepGen);
	Double_t maxptPi = hMCptPi1D->GetMaximum();
	hMCptPi1D->GetYaxis()->SetRangeUser(0,maxptPi*1.2);
	hMCptPi1D->Sumw2();
	hMCptPi1D->SetMinimum(0.01);
	hMCptPi1D->SetLineColor(2);
	hMCptPi1D->SetLineWidth(3);
	hMCptPi1D->SetMarkerColor(2);
	hMCptPi1D->SetFillColor(2);
	hMCptPi1D->SetFillStyle(3005);
	hMCptPi1D->SetMarkerStyle(20);
	hMCptPi1D->GetXaxis()->SetTitleOffset(1.2);
	hMCptPi1D->GetXaxis()->SetTitle("p_{T, #pi} (GeV/c), MC data");
	hMCptPi1D->Draw("hist");
	cmcptPi->cd();
	cmcptPi->SetLeftMargin(0.15);
	cmcptPi->SetRightMargin(0.05);
	hMCptPi1D->Draw("hist");
	hMCptPi1D->Draw("err same");
	cmcptPi->Update();

	cmc2->cd(2);
	TH1D *hRECptPi1D = data->ShowProjection(iptPi, stepRec);
	hRECptPi1D->GetYaxis()->SetRangeUser(0,maxptPi*1.2);
	hRECptPi1D->SetLineColor(4);
	hRECptPi1D->SetLineWidth(3);
	hRECptPi1D->SetMarkerColor(4);
	hRECptPi1D->SetFillColor(4);
	hRECptPi1D->SetFillStyle(3004);
	hRECptPi1D->SetMarkerStyle(20);
	hRECptPi1D->GetXaxis()->SetTitleOffset(1.2);
	hRECptPi1D->Sumw2();
	hRECptPi1D->SetMinimum(0.01);
	hRECptPi1D->GetXaxis()->SetTitle("p_{T, #pi} (GeV/c), AOD");
	hRECptPi1D->Draw("hist");
	cptPi->cd();
	cptPi->SetLeftMargin(0.15);
	cptPi->SetRightMargin(0.05);
	hRECptPi1D->Draw("hist");
	hRECptPi1D->Draw("err same");
	cptPi->Update();

	cmc2->cd(3);
	TH1D *hMCptK1D = data->ShowProjection(iptK, stepGen);
	Double_t maxptK = hMCptK1D->GetMaximum();
	hMCptK1D->GetYaxis()->SetRangeUser(0,maxptK*1.2);
	hMCptK1D->SetLineColor(2);
	hMCptK1D->SetLineWidth(3);
	hMCptK1D->SetMarkerColor(2);
	hMCptK1D->SetFillColor(2);
	hMCptK1D->SetFillStyle(3005);
	hMCptK1D->SetMarkerStyle(20);
	hMCptK1D->GetXaxis()->SetTitleOffset(1.2);
	hMCptK1D->Sumw2();
	hMCptK1D->SetMinimum(0.01);
	hMCptK1D->GetXaxis()->SetTitle("p_{T, K} (GeV/c), MC data");
	hMCptK1D->Draw("hist");
	cmcptK->cd();
	cmcptK->SetLeftMargin(0.15);
	cmcptK->SetRightMargin(0.05);
	hMCptK1D->Draw("hist");
	hMCptK1D->Draw("err same");
	cmcptK->Update();

	cmc2->cd(4);
	TH1D *hRECptK1D = data->ShowProjection(iptK, stepRec);
	hRECptK1D->GetYaxis()->SetRangeUser(0,maxptK*1.2);
	hRECptK1D->SetLineColor(4);
	hRECptK1D->SetLineWidth(3);
	hRECptK1D->SetMarkerColor(4);
	hRECptK1D->SetFillColor(4);
	hRECptK1D->SetFillStyle(3004);
	hRECptK1D->SetMarkerStyle(20);
	hRECptK1D->GetXaxis()->SetTitleOffset(1.2);
	hRECptK1D->Sumw2();
	hRECptK1D->SetMinimum(0.01);
	hRECptK1D->Draw("hist");
	cptK->cd();
	cptK->SetLeftMargin(0.15);
	cptK->SetRightMargin(0.05);
	hRECptK1D->GetXaxis()->SetTitle("p_{T, K} (GeV/c), AOD");
	hRECptK1D->Draw("hist");
	hRECptK1D->Draw("err same");
	cptK->Update();

	cmc2->cd(5);
	TH1D *hMCcT1D = data->ShowProjection(icT, stepGen);
	Double_t maxcT = hMCcT1D->GetMaximum();
	hMCcT1D->GetYaxis()->SetRangeUser(0,maxcT*1.2);
	hMCcT1D->SetLineColor(2);
	hMCcT1D->SetLineWidth(3);
	hMCcT1D->SetMarkerColor(2);
	hMCcT1D->SetFillColor(2);
	hMCcT1D->SetFillStyle(3005);
	hMCcT1D->SetMarkerStyle(20);
	hMCcT1D->GetXaxis()->SetTitleOffset(1.2);
	hMCcT1D->Sumw2();
	hMCcT1D->SetMinimum(0.01);
	hMCcT1D->GetXaxis()->SetTitle("ct (#mum), MC data");
	hMCcT1D->Draw("hist");
	cmccT->cd();
	cmccT->SetLeftMargin(0.15);
	cmccT->SetRightMargin(0.05);
	hMCcT1D->Draw("hist");
	hMCcT1D->Draw("err same");
	cmccT->Update();

	cmc2->cd(6);
	TH1D *hRECcT1D = data->ShowProjection(icT, stepRec);
	hRECcT1D->GetYaxis()->SetRangeUser(0,maxcT*1.2);
	hRECcT1D->SetLineColor(4);
	hRECcT1D->SetLineWidth(3);
	hRECcT1D->SetMarkerColor(4);
	hRECcT1D->SetFillColor(4);
	hRECcT1D->SetFillStyle(3004);
	hRECcT1D->SetMarkerStyle(20);
	hRECcT1D->GetXaxis()->SetTitleOffset(1.2);
	hRECcT1D->Sumw2();
	hRECcT1D->SetMinimum(0.01);
	hRECcT1D->Draw("hist");
	ccT->cd();
	ccT->SetLeftMargin(0.15);
	ccT->SetRightMargin(0.05);
	hRECcT1D->GetXaxis()->SetTitle("c#t (#mum), AOD");
	hRECcT1D->Draw("hist");
	hRECcT1D->Draw("err same");
	ccT->Update();

	/*
	// printing on eps files
	cmc1->Print("Plots/dataMC_pt_y_cTS.gif");
	cmc2->Print("Plots/dataMC_ptPi_ptK_cT.gif");
	cmcpt->Print("Plots/pt_Gen.eps");
	cmcy->Print("Plots/y_Gen.eps");
	cmccTS->Print("Plots/cTS_Gen.eps");
	cmcptPi->Print("Plots/ptPi_Gen.eps");
	cmcptK->Print("Plots/ptK_Gen.eps");
	cmccT->Print("Plots/cT_Gen.eps");
	cpt->Print("Plots/pt_Rec.eps");
	cy->Print("Plots/y_Rec.eps");
	ccTS->Print("Plots/cTS_Rec.eps");
	cptPi->Print("Plots/ptPi_Rec.eps");
	cptK->Print("Plots/ptK_Rec.eps");
	ccT->Print("Plots/cT_Rec.eps");
	cGen2d->Print("Plots/pt_y_Gen_2D.eps");  
	cRec2d->Print("Plots/pt_y_Rec_2D.eps");  

	// printing on gif files
	cmc1->Print("Plots/dataMC_pt_y_cTS.gif");
	cmc2->Print("Plots/dataMC_ptPi_ptK_cT.gif");
	cmcpt->Print("Plots/pt_Gen.eps");
	cmcy->Print("Plots/y_Gen.eps");
	cmccTS->Print("Plots/cTS_Gen.eps");
	cmcptPi->Print("Plots/ptPi_Gen.eps");
	cmcptK->Print("Plots/ptK_Gen.eps");
	cmccT->Print("Plots/cT_Gen.eps");
	cpt->Print("Plots/pt_Rec.eps");
	cy->Print("Plots/y_Rec.eps");
	ccTS->Print("Plots/cTS_Rec.eps");
	cptPi->Print("Plots/ptPi_Rec.eps");
	cptK->Print("Plots/ptK_Rec.eps");
	ccT->Print("Plots/cT_Rec.eps");
	cGen2d->Print("Plots/pt_y_Gen_2D.eps");  
	cRec2d->Print("Plots/pt_y_Rec_2D.eps");  
	*/

	//construct the efficiency grid from the data container 
	AliCFEffGrid *eff = new AliCFEffGrid("eff"," The efficiency",*data);
	eff->CalculateEfficiency(stepRec,stepGen); //eff= step1/step0

	//The efficiency along the variables, and some 2-D projections
	TCanvas *ceff =new TCanvas("ceff"," Efficiency",0,0,1600,1200);
	ceff->Divide(3,2);
	TCanvas *ceff2D =new TCanvas("ceff2D"," Efficiency for pt and y",50,50,550,550);
	TCanvas *ceffpt = new TCanvas("ceffpt","Efficiency vs pt",50,50,550,550);
	TCanvas *ceffy = new TCanvas("ceffy","Efficiency vs y",50,50,550,550);
	TCanvas *ceffcTS = new TCanvas("ceffcTS","Efficiency vs cosThetaStar",50,50,550,550);
	TCanvas *ceffptPi = new TCanvas("ceffptPi","Efficiency vs ptPi",50,50,550,550);
	TCanvas *ceffptK = new TCanvas("ceffptK","Efficiency vs ptK",50,50,550,550);
	TCanvas *ceffcT = new TCanvas("ceffcT","Efficiency vs cT",50,50,550,550);
	TCanvas *ceff2Dtext = new TCanvas("ceff2Dtext","Text plot for efficiency in pt and y",50,50,550,550);

	ceff->cd(1);
	TH1D *hpteffCF = eff->Project(ipt); //the efficiency vs pt
	hpteffCF->Sumw2();
	//hpteffCF->SetMinimum(0.01);
	hpteffCF->SetLineColor(8);
	hpteffCF->SetLineWidth(3);
	hpteffCF->SetMarkerColor(8);
	hpteffCF->SetMarkerStyle(20);
	hpteffCF->GetXaxis()->SetTitleOffset(1.2);
	hpteffCF->GetYaxis()->SetTitleOffset(1.5);
	hpteffCF->Draw("hist");
	ceffpt->cd();
	ceffpt->SetLeftMargin(0.15);
	ceffpt->SetRightMargin(0.05);
	hpteffCF->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hpteffCF->GetYaxis()->SetTitle("Efficiency");
	hpteffCF->Draw("hist");
	hpteffCF->Draw("err same");

	ceff->cd(2);
	TH1D *hyeffCF = eff->Project(iy); //the efficiency vs y
	hyeffCF->Sumw2();
	//hyeffCF->SetMinimum(0.01);
	hyeffCF->SetLineColor(8);
	hyeffCF->SetLineWidth(3);
	hyeffCF->SetMarkerColor(8);
	hyeffCF->SetMarkerStyle(20);
	hyeffCF->GetXaxis()->SetTitleOffset(1.2);
	hyeffCF->GetYaxis()->SetTitleOffset(1.5);
	hyeffCF->Draw("hist");
	ceffy->cd();
	ceffy->SetLeftMargin(0.15);
	ceffy->SetRightMargin(0.05);
	hyeffCF->GetXaxis()->SetTitle("y");
	hyeffCF->GetYaxis()->SetTitle("Efficiency");
	hyeffCF->Draw("hist");
	hyeffCF->Draw("err same");

	ceff->cd(3);
	TH1D *hcTSeffCF = eff->Project(icTS); //the efficiency vs cosThetaStar
	hcTSeffCF->Sumw2();
	//hcTSeffCF->SetMinimum(0.01);
	hcTSeffCF->SetLineColor(8);
	hcTSeffCF->SetLineWidth(3);
	hcTSeffCF->SetMarkerColor(8);
	hcTSeffCF->SetMarkerStyle(20);
	hcTSeffCF->GetXaxis()->SetTitleOffset(1.2);
	hcTSeffCF->GetYaxis()->SetTitleOffset(1.5);
	hcTSeffCF->Draw("hist");
	ceffcTS->cd();
	ceffcTS->SetLeftMargin(0.15);
	ceffcTS->SetRightMargin(0.05);
	hcTSeffCF->GetXaxis()->SetTitle("cosThetaStar");
	hcTSeffCF->GetYaxis()->SetTitle("Efficiency");
	hcTSeffCF->Draw("hist");
	hcTSeffCF->Draw("err same");

	ceff->cd(4);
	TH1D *hptPieffCF = eff->Project(iptPi); //the efficiency vs ptPi
	hptPieffCF->Sumw2();
	//hptPieffCF->SetMinimum(0.01);
	hptPieffCF->SetLineColor(8);
	hptPieffCF->SetLineWidth(3);
	hptPieffCF->SetMarkerColor(8);
	hptPieffCF->SetMarkerStyle(20);
	hptPieffCF->GetXaxis()->SetTitleOffset(1.2);
	hptPieffCF->GetYaxis()->SetTitleOffset(1.5);
	hptPieffCF->Draw("hist");
	ceffptPi->cd();
	ceffptPi->SetLeftMargin(0.15);
	ceffptPi->SetRightMargin(0.05);
	hptPieffCF->GetXaxis()->SetTitle("p_{T, #pi} (GeV/c)");
	hptPieffCF->GetYaxis()->SetTitle("Efficiency");
	hptPieffCF->Draw("hist");
	hptPieffCF->Draw("err same");

	ceff->cd(5);
	TH1D *hptKeffCF = eff->Project(iptK); //the efficiency vs ptK
	hptKeffCF->Sumw2();
	//hptKeffCF->SetMinimum(0.01);
	hptKeffCF->SetLineColor(8);
	hptKeffCF->SetLineWidth(3);
	hptKeffCF->SetMarkerColor(8);
	hptKeffCF->SetMarkerStyle(20);
	hptKeffCF->GetXaxis()->SetTitleOffset(1.2);
	hptKeffCF->GetYaxis()->SetTitleOffset(1.5);
	hptKeffCF->Draw("hist");
	ceffptK->cd();
	ceffptK->SetLeftMargin(0.15);
	ceffptK->SetRightMargin(0.05);
	hptKeffCF->GetXaxis()->SetTitle("p_{T, K} (GeV/c)");
	hptKeffCF->GetYaxis()->SetTitle("Efficiency");
	hptKeffCF->Draw("hist");
	hptKeffCF->Draw("err same");

	ceff->cd(6);
	TH1D *hcTeffCF = eff->Project(icT); //the efficiency vs cT
	hcTeffCF->Sumw2();
	//hcTeffCF->SetMinimum(0.01);
	hcTeffCF->SetLineColor(8);
	hcTeffCF->SetLineWidth(3);
	hcTeffCF->SetMarkerColor(8);
	hcTeffCF->SetMarkerStyle(20);
	hcTeffCF->GetXaxis()->SetTitleOffset(1.2);
	hcTeffCF->GetYaxis()->SetTitleOffset(1.5);
	hcTeffCF->Draw("hist");
	ceffcT->cd();
	ceffcT->SetLeftMargin(0.15);
	ceffcT->SetRightMargin(0.05);
	hcTeffCF->GetXaxis()->SetTitle("c#t (#mum)");
	hcTeffCF->GetYaxis()->SetTitle("Efficiency");
	hcTeffCF->Draw("hist");
	hcTeffCF->Draw("err same");

	ceff2D->cd();
	TH2D *hptyeffCF = eff->Project(ipt,iy); //look at the numerator
	//hptyeffCF->SetMinimum(0.01);
	hptyeffCF->SetMarkerColor(8);
	hptyeffCF->SetLineColor(8);
	hptyeffCF->SetMinimum(0.01);
	hptyeffCF->Draw("lego");
	ceff2Dtext->cd();
	hptyeffCF->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hptyeffCF->GetXaxis()->SetTitleOffset(1.2);
	hptyeffCF->GetYaxis()->SetTitle("y");
	hptyeffCF->GetYaxis()->SetTitleOffset(1.2);
	hptyeffCF->Draw("text");

	/*
	// printing eps files
	ceff->Print("Plots/efficiencies.eps");
	ceffpt->Print("Plots/effpt.eps");
	ceffy->Print("Plots/effy.eps");
	ceffcTS->Print("Plots/effcTS.eps");
	ceffptPi->Print("Plots/effptPi.eps");
	ceffptK->Print("Plots/effptK.eps");
	ceffcT->Print("Plots/effcT.eps");
	ceff2D->Print("Plots/eff2D_pt_y.eps");
	ceff2Dtext->Print("Plots/eff2d_text_pt_y.eps");

	// printing gif files
	ceff->Print("Plots/efficiencies.gif");
	ceffpt->Print("Plots/effpt.gif");
	ceffy->Print("Plots/effy.gif");
	ceffcTS->Print("Plots/effcTS.gif");
	ceffptPi->Print("Plots/effptPi.gif");
	ceffptK->Print("Plots/effptK.gif");
	ceffcT->Print("Plots/effcT.gif");
	ceff2D->Print("Plots/eff2D_pt_y.gif");
	ceff2Dtext->Print("Plots/eff2d_text_pt_y.gif");
	*/

	// applying efficiencies - using projections, not the ApplyEffCorrection method

	TCanvas *cmultpiplypt = new TCanvas("cmultiplypt","Reco From Eff in pt distribution",50,50,550,550);
	TCanvas *cmultpiplyy = new TCanvas("cmultiplyy","Reco From Eff in y distribution",50,50,550,550);
	TCanvas *cmultpiplycTS = new TCanvas("cmultiplycTS","Reco From Eff in cTS distribution",50,50,550,550);
	TCanvas *cmultpiplyptPi = new TCanvas("cmultiplyptPi","Reco From Eff in ptPi distribution",50,50,550,550);
	TCanvas *cmultpiplyptK = new TCanvas("cmultiplyptK","Reco From Eff in ptK distribution",50,50,550,550);
	TCanvas *cmultpiplycT = new TCanvas("cmultiplycT","Reco From Eff in cT distribution",50,50,550,550);

	TH1D *hmultiplypt = new TH1D("hmultiplypt","hmultiplypt",13,0,10);
	cout << " bin for histo MC = " << hMCpt1D->GetNbinsX() << " while for RECO histo = " << hRECpt1D->GetNbinsX() << " while for efficiency histo = " << hpteffCF->GetNbinsX() << endl;

	const Double_t ptmin_0_4 =  0.0 ;
	const Double_t ptmax_0_4 =  4.0 ;
	const Double_t ptmin_4_8 =  4.0 ;
	const Double_t ptmax_4_8 =  8.0 ;
	const Double_t ptmin_8_10 =  8.0 ;
	const Double_t ptmax_8_10 =  10.0 ;
	const Int_t nbin0_0_4  = 8 ; //bins in pt
	const Int_t nbin0_4_8  = 4 ; //bins in pt
	const Int_t nbin0_8_10  = 1 ; //bins in pt
	const Int_t nbins = nbin0_0_4 + nbin0_4_8 + nbin0_8_10;
	Double_t binLim0[nbins+1];
	for(Int_t i=0; i<=nbin0_0_4; i++) binLim0[i]=(Double_t)ptmin_0_4 + (ptmax_0_4-ptmin_0_4)/nbin0_0_4*(Double_t)i ; 
	if (binLim0[nbin0_0_4] != ptmin_4_8)  {
		printf("Calculated bin lim for pt - 1st range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin0_4_8; i++) binLim0[i+nbin0_0_4]=(Double_t)ptmin_4_8 + (ptmax_4_8-ptmin_4_8)/nbin0_4_8*(Double_t)i ; 
	if (binLim0[nbin0_0_4+nbin0_4_8] != ptmin_8_10)  {
		printf("Calculated bin lim for pt - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin0_8_10; i++) binLim0[i+nbin0_0_4+nbin0_4_8]=(Double_t)ptmin_8_10 + (ptmax_8_10-ptmin_8_10)/nbin0_8_10*(Double_t)i ; 

	hmultiplypt->SetBins(nbins,binLim0);
	hmultiplypt->Divide(hRECpt1D,hpteffCF,1,1,"b");
	hmultiplypt->SetMinimum(0);
	hmultiplypt->SetMaximum(hMCpt1D->GetMaximum());
	hmultiplypt->SetLineColor(38);
	hmultiplypt->SetLineWidth(3);
	hmultiplypt->SetMarkerColor(38);
	hmultiplypt->SetMarkerStyle(20);
	hmultiplypt->GetXaxis()->SetTitleOffset(1.2);
	hmultiplypt->GetYaxis()->SetTitleOffset(1.5);
	hmultiplypt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	cmultiplypt->SetLeftMargin(0.15);
	cmultiplypt->SetRightMargin(0.05);
	hmultiplypt->SetFillStyle(3004);
	hmultiplypt->SetFillColor(38);
	cmultiplypt->cd();
	hmultiplypt->Draw("hist err");
	hMCpt1D->Draw("hist same err");

	TH1D *hmultiplyy = new TH1D("hmultiplyy","hmultiplyy",8,-2,2);
	hmultiplyy->Divide(hRECy1D,hyeffCF,1,1,"b");
	hmultiplyy->SetMinimum(0);
	hmultiplyy->SetMaximum(hMCy1D->GetMaximum());
	hmultiplyy->SetLineColor(38);
	hmultiplyy->SetLineWidth(3);
	hmultiplyy->SetMarkerColor(38);
	hmultiplyy->SetMarkerStyle(20);
	hmultiplyy->GetXaxis()->SetTitleOffset(1.2);
	hmultiplyy->GetXaxis()->SetTitle("y");
	hmultiplyy->GetYaxis()->SetTitleOffset(1.5);
	cmultiplyy->SetLeftMargin(0.15);
	cmultiplyy->SetRightMargin(0.05);
	cmultiplyy->cd();
	hmultiplyy->SetFillStyle(3004);
	hmultiplyy->SetFillColor(38);
	hmultiplyy->Draw("hist err");
	hMCy1D->Draw("hist same err");

	TH1D *hmultiplycTS = new TH1D("hmultiplycTS","hmultiplycTS",8,-1,1);
	hmultiplycTS->Divide(hRECcTS1D,hcTSeffCF,1,1,"b");
	hmultiplycTS->SetMinimum(0);
	hmultiplycTS->SetMaximum(hMCcTS1D->GetMaximum());
	hmultiplycTS->SetLineColor(38);
	hmultiplycTS->SetLineWidth(3);
	hmultiplycTS->SetMarkerColor(38);
	hmultiplycTS->SetMarkerStyle(20);
	hmultiplycTS->GetXaxis()->SetTitleOffset(1.2);
	hmultiplycTS->GetXaxis()->SetTitle("cosThetaStar");
	hmultiplycTS->GetYaxis()->SetTitleOffset(1.5);
	cmultiplycTS->SetLeftMargin(0.15);
	cmultiplycTS->SetRightMargin(0.05);
	cmultiplycTS->cd();
	hmultiplycTS->SetFillStyle(3004);
	hmultiplycTS->SetFillColor(38);
	hmultiplycTS->Draw("hist err");
	hMCcTS1D->Draw("hist same err");

	TH1D *hmultiplyptPi = new TH1D("hmultiplyptPi","hmultiplyptPi",13,0,10);
	hmultiplyptPi->SetBins(nbins,binLim0);
	hmultiplyptPi->Divide(hRECptPi1D,hptPieffCF,1,1,"b");
	hmultiplyptPi->SetMinimum(0);
	hmultiplyptPi->SetMaximum(hMCptPi1D->GetMaximum());
	hmultiplyptPi->SetLineColor(38);
	hmultiplyptPi->SetLineWidth(3);
	hmultiplyptPi->SetMarkerColor(38);
	hmultiplyptPi->SetMarkerStyle(20);
	hmultiplyptPi->GetXaxis()->SetTitleOffset(1.2);
	hmultiplyptPi->GetXaxis()->SetTitle("p_{T, #pi} (GeV/c)");
	hmultiplyptPi->GetYaxis()->SetTitleOffset(1.5);
	cmultiplyptPi->SetLeftMargin(0.15);
	cmultiplyptPi->SetRightMargin(0.05);
	cmultiplyptPi->cd();
	hmultiplyptPi->SetFillStyle(3004);
	hmultiplyptPi->SetFillColor(38);
	hmultiplyptPi->Draw("hist err");
	hMCptPi1D->Draw("hist same err");

	TH1D *hmultiplyptK = new TH1D("hmultiplyptK","hmultiplyptK",13,0,10);
	hmultiplyptK->SetBins(nbins,binLim0);
	hmultiplyptK->Divide(hRECptK1D,hptKeffCF,1,1,"b");
	hmultiplyptK->SetMinimum(0);
	hmultiplyptK->SetMaximum(hMCptK1D->GetMaximum());
	hmultiplyptK->SetLineColor(38);
	hmultiplyptK->SetLineWidth(3);
	hmultiplyptK->SetMarkerColor(38);
	hmultiplyptK->SetMarkerStyle(20);
	hmultiplyptK->GetXaxis()->SetTitleOffset(1.2);
	hmultiplyptK->GetXaxis()->SetTitle("p_{T, K} (GeV/c)");
	hmultiplyptK->GetYaxis()->SetTitleOffset(1.5);
	cmultiplyptK->SetLeftMargin(0.15);
	cmultiplyptK->SetRightMargin(0.05);
	cmultiplyptK->cd();
	hmultiplyptK->SetFillStyle(3004);
	hmultiplyptK->SetFillColor(38);
	hmultiplyptK->Draw("hist err");
	hMCptK1D->Draw("hist same err");

	TH1D *hmultiplycT = new TH1D("hmultiplycT","hmultiplycT",24,0,500);
	hmultiplycT->Divide(hRECcT1D,hcTeffCF,1,1,"b");
	hmultiplycT->SetMinimum(0);
	hmultiplycT->SetMaximum(hMCcT1D->GetMaximum());
	hmultiplycT->SetLineColor(38);
	hmultiplycT->SetLineWidth(3);
	hmultiplycT->SetMarkerColor(38);
	hmultiplycT->SetMarkerStyle(20);
	hmultiplycT->GetXaxis()->SetTitleOffset(1.2);
	hmultiplycT->GetXaxis()->SetTitle("c#t (#mum)");
	hmultiplycT->GetYaxis()->SetTitleOffset(1.5);
	cmultiplycT->SetLeftMargin(0.15);
	cmultiplycT->SetRightMargin(0.05);
	cmultiplycT->cd();
	hmultiplycT->SetFillStyle(3004);
	hmultiplycT->SetFillColor(38);
	hmultiplycT->Draw("hist err");
	hMCcT1D->Draw("hist same err");

	/*
	TFile* file_histo = new TFile("fileHisto_180004.root","RECREATE");
	hMCpt1D->Write("hMCpt1D");
	hRECpt1D->Write("hRECpt1D");
	hMCy1D->Write("hMCy1D");
	hRECy1D->Write("hRECy1D");
	hMCcTS1D->Write("hMCcTS1D");
	hRECcTS1D->Write("hRECcTS1D");
	hMCptPi1D->Write("hMCptPi1D");
	hRECptPi1D->Write("hRECptPi1D");
	hMCptK1D->Write("hMCptK1D");
	hRECptK1D->Write("hRECptK1D");
	hMCcT1D->Write("hMCcT1D");
	hRECcT1D->Write("hRECcT1D");
	hpteffCF->Write("hpteffCF");
	hyeffCF->Write("hyeffCF");
	hcTSeffCF->Write("hcTSeffCF");
	hptPieffCF->Write("hptPieffCF");
	hptKeffCF->Write("hptKeffCF");
	hcTeffCF->Write("hcTeffCF");
	hmultiplypt->Write("");
	hmultiplyy->Write("");
	hmultiplycTS->Write("");
	hmultiplyptPi->Write("");
	hmultiplyptK->Write("");
	hmultiplycT->Write("");
	file_histo->Close();
	*/

	Int_t nsliceVars = 5;
	Int_t sliceVars[5];
	sliceVars[0] = 1; 
	sliceVars[1] = 2; 
	sliceVars[2] = 3; 
	sliceVars[3] = 4; 
	sliceVars[4] = 5; 
	Double_t varMin[6] = {0.5, -2, -1, 0., 0., 0. };
	Double_t varMax[6] = {1.5, 2, 1, 10., 10., 500.};
	AliCFContainer* data_sliced = data->MakeSlice(nsliceVars, sliceVars, varMin, varMax);
	cout << " the container now has " << data_sliced->GetNVar() << " dimensions " << endl;
	
}
