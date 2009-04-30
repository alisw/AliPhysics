// macro to applying slice method to calculate 
// efficiencies, on a CF container created with 
// AliCFHeavyFlavourTaskMultiVarMultiStep (commenting cut on
// cosPoitingAngle at REC PPR step)
// slicing on cosPointingAngle
// The efficiency are calculated by default using RecPPR / GenAcc steps


#include <Riostream.h>

extern TRandom *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem *gSystem;

void ReadCFHeavyFlavourOutputSlice(Int_t stepDen=1,
				   Int_t stepNum=5,
		   const char* plotsDir="./SlicePlots/GenAcc_RecPPR")
{

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetOptTitle(0);
	
	gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include  -I$ROOTSYS/include");
	gSystem->Load("libANALYSIS.so");
	gSystem->Load("libANALYSISalice.so");
	gSystem->Load("$ALICE_ROOT/CORRFW/libCORRFW.so") ;
	
	// Expected container: 
	// number of selection steps
	// 0 --> MC
        // 1 --> MC Acceptance
        // 2 --> Reco
        // 3 --> Reco in Acceptance
        // 4 --> Reco with n. of ITS clusters
	// 5 --> Reco PPR cuts		
	// variables on the grid: pt, y, cosThetaStar, ptPi, ptK, cTau, dca, d0pi, d0K, d0xd0, cosPointingAngle
		
	// Flag the sel steps
	Int_t stepGen=0;
	Int_t stepGenAcc=1;
	Int_t stepRec=2;
	Int_t stepRecAcc=3;
	Int_t stepRecITSClus=4;
	Int_t stepRecPPR=5;
	
	//the sensitive variables, their indeces
	Int_t ipt =0;
	Int_t iy=1;
	UInt_t icTS  = 2;   //cTS stands for cosThetaStar
	UInt_t iptPi  = 3;
	UInt_t iptK  = 4;
	UInt_t icT  = 5;
	UInt_t idca  = 6;
	UInt_t id0pi  = 7;
	UInt_t id0K  = 8;
	UInt_t id0xd0  = 9;
	UInt_t ipointing  = 10;
	
	// Read the  container from file
	TFile *file = new TFile("./output.root");
	AliCFContainer *data = (AliCFContainer*) (file->Get("container"));
	
	TH1D *hMCpt1D = data->ShowProjection(ipt, stepGen);
	TH1D *hMCAccpt1D = data->ShowProjection(ipt, stepGenAcc);
	TH1D *hRECpt1D = data->ShowProjection(ipt, stepRec);
	TH1D *hRECAccpt1D = data->ShowProjection(ipt, stepRecAcc);
	TH1D *hRECITSClpt1D = data->ShowProjection(ipt, stepRecITSClus);
	TH1D *hRECPPRpt1D = data->ShowProjection(ipt, stepRecPPR);

	TH1D *hMCy1D = data->ShowProjection(iy, stepGen);
	TH1D *hMCAccy1D = data->ShowProjection(iy, stepGenAcc);
	TH1D *hRECy1D = data->ShowProjection(iy, stepRec);
	TH1D *hRECAccy1D = data->ShowProjection(iy, stepRecAcc);
	TH1D *hRECITSCly1D = data->ShowProjection(iy, stepRecITSClus);
	TH1D *hRECPPRy1D = data->ShowProjection(iy, stepRecPPR);

	TH1D *hMCcTS1D = data->ShowProjection(icTS, stepGen);
	TH1D *hMCAcccTS1D = data->ShowProjection(icTS, stepGenAcc);
	TH1D *hRECcTS1D = data->ShowProjection(icTS, stepRec);
	TH1D *hRECAcccTS1D = data->ShowProjection(icTS, stepRecAcc);
	TH1D *hRECITSClcTS1D = data->ShowProjection(icTS, stepRecITSClus);
	TH1D *hRECPPRcTS1D = data->ShowProjection(icTS, stepRecPPR);

	TH1D *hMCptPi1D = data->ShowProjection(iptPi, stepGen);
	TH1D *hMCAccptPi1D = data->ShowProjection(iptPi, stepGen);
	TH1D *hRECptPi1D = data->ShowProjection(iptPi, stepRec);
	TH1D *hRECAccptPi1D = data->ShowProjection(iptPi, stepRecAcc);
	TH1D *hRECITSClptPi1D = data->ShowProjection(iptPi, stepRecITSClus);
	TH1D *hRECPPRptPi1D = data->ShowProjection(iptPi, stepRecPPR);

	TH1D *hMCptK1D = data->ShowProjection(iptK, stepGen);
	TH1D *hMCAccptK1D = data->ShowProjection(iptK, stepGenAcc);
	TH1D *hRECptK1D = data->ShowProjection(iptK, stepRec);
	TH1D *hRECAccptK1D = data->ShowProjection(iptK, stepRecAcc);
	TH1D *hRECITSClptK1D = data->ShowProjection(iptK, stepRecITSClus);
	TH1D *hRECPPRptK1D = data->ShowProjection(iptK, stepRecPPR);

	TH1D *hMCcT1D = data->ShowProjection(icT, stepGen);
	TH1D *hMCAcccT1D = data->ShowProjection(icT, stepGenAcc);
	TH1D *hRECcT1D = data->ShowProjection(icT, stepRec);
	TH1D *hRECAcccT1D = data->ShowProjection(icT, stepRecAcc);
	TH1D *hRECITSClcT1D = data->ShowProjection(icT, stepRecITSClus);
	TH1D *hRECPPRcT1D = data->ShowProjection(icT, stepRecPPR);

	TH1D *hMCdca1D = data->ShowProjection(idca, stepGen);
	TH1D *hMCAccdca1D = data->ShowProjection(idca, stepGenAcc);
	TH1D *hRECdca1D = data->ShowProjection(idca, stepRec);
	TH1D *hRECAccdca1D = data->ShowProjection(idca, stepRecAcc);
	TH1D *hRECITSCldca1D = data->ShowProjection(idca, stepRecITSClus);
	TH1D *hRECPPRdca1D = data->ShowProjection(idca, stepRecPPR);

	TH1D *hMCd0pi1D = data->ShowProjection(id0pi, stepGen);
	TH1D *hMCAccd0pi1D = data->ShowProjection(id0pi, stepGenAcc);
	TH1D *hRECd0pi1D = data->ShowProjection(id0pi, stepRec);
	TH1D *hRECAccd0pi1D = data->ShowProjection(id0pi, stepRecAcc);
	TH1D *hRECITSCld0pi1D = data->ShowProjection(id0pi, stepRecITSClus);
	TH1D *hRECPPRd0pi1D = data->ShowProjection(id0pi, stepRecPPR);

	TH1D *hMCd0K1D = data->ShowProjection(id0K, stepGen);
	TH1D *hMCAccd0K1D = data->ShowProjection(id0K, stepGenAcc);
	TH1D *hRECd0K1D = data->ShowProjection(id0K, stepRec);
	TH1D *hRECAccd0K1D = data->ShowProjection(id0K, stepRecAcc);
	TH1D *hRECITSCld0K1D = data->ShowProjection(id0K, stepRecITSClus);
	TH1D *hRECPPRd0K1D = data->ShowProjection(id0K, stepRecPPR);

	TH1D *hMCd0xd01D = data->ShowProjection(id0xd0, stepGen);
	TH1D *hMCAccd0xd01D = data->ShowProjection(id0xd0, stepGenAcc);
	TH1D *hRECd0xd01D = data->ShowProjection(id0xd0, stepRec);
	TH1D *hRECAccd0xd01D = data->ShowProjection(id0xd0, stepRecAcc);
	TH1D *hRECITSCld0xd01D = data->ShowProjection(id0xd0, stepRecITSClus);
	TH1D *hRECPPRd0xd01D = data->ShowProjection(id0xd0, stepRecPPR);

	TH1D *hMCpointing1D = data->ShowProjection(ipointing, stepGen);
	TH1D *hMCAccpointing1D = data->ShowProjection(ipointing, stepGenAcc);
	TH1D *hRECpointing1D = data->ShowProjection(ipointing, stepRec);
	TH1D *hRECAccpointing1D = data->ShowProjection(ipointing, stepRecAcc);
	TH1D *hRECITSClpointing1D = data->ShowProjection(ipointing, stepRecITSClus);
	TH1D *hRECPPRpointing1D = data->ShowProjection(ipointing, stepRecPPR);

	//construct the efficiency grid from the data container 
	AliCFEffGrid *eff = new AliCFEffGrid("eff"," The efficiency",*data);
	eff->CalculateEfficiency(stepNum,stepDen); //eff= step1/step0

	//The efficiency along the variables
	TCanvas *ceff =new TCanvas("ceff"," Efficiency",0,0,1600,1200);
	ceff->Divide(4,3);
	TCanvas *ceffpt = new TCanvas("ceffpt","Efficiency vs pt",50,50,550,550);
	TCanvas *ceffy = new TCanvas("ceffy","Efficiency vs y",50,50,550,550);
	TCanvas *ceffcTS = new TCanvas("ceffcTS","Efficiency vs cosThetaStar",50,50,550,550);
	TCanvas *ceffptPi = new TCanvas("ceffptPi","Efficiency vs ptPi",50,50,550,550);
	TCanvas *ceffptK = new TCanvas("ceffptK","Efficiency vs ptK",50,50,550,550);
	TCanvas *ceffcT = new TCanvas("ceffcT","Efficiency vs cT",50,50,550,550);
	TCanvas *ceffdca = new TCanvas("ceffdca","Efficiency vs dca",50,50,550,550);
	TCanvas *ceffd0pi = new TCanvas("ceffd0pi","Efficiency vs d0pi",50,50,550,550);
	TCanvas *ceffd0K = new TCanvas("ceffd0K","Efficiency vs d0K",50,50,550,550);
	TCanvas *ceffd0xd0 = new TCanvas("ceffd0xd0","Efficiency vs d0xd0",50,50,550,550);
	TCanvas *ceffpointing = new TCanvas("ceffpointing","Efficiency vs pointing",50,50,550,550);

	ceff->cd(1);
	TH1D *hpteffCF = eff->Project(ipt); //the efficiency vs pt
	//hpteffCF->Setw2();
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
	//hyeffCF->Setw2();
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
	//hcTSeffCF->Setw2();
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
	//hptPieffCF->Setw2();
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
	//hptKeffCF->Setw2();
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
	//hcTeffCF->Setw2();
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
	hcTeffCF->GetXaxis()->SetTitle("ct (#mum)");
	hcTeffCF->GetYaxis()->SetTitle("Efficiency");
	hcTeffCF->Draw("hist");
	hcTeffCF->Draw("err same");

	ceff->cd(7);
	TH1D *hdcaeffCF = eff->Project(idca); //the efficiency vs dca
	//hdcaeffCF->Setw2();
	//hdcaeffCF->SetMinimum(0.01);
	hdcaeffCF->SetLineColor(8);
	hdcaeffCF->SetLineWidth(3);
	hdcaeffCF->SetMarkerColor(8);
	hdcaeffCF->SetMarkerStyle(20);
	hdcaeffCF->GetXaxis()->SetTitleOffset(1.2);
	hdcaeffCF->GetYaxis()->SetTitleOffset(1.5);
	hdcaeffCF->Draw("hist");
	ceffdca->cd();
	ceffdca->SetLeftMargin(0.15);
	ceffdca->SetRightMargin(0.05);
	hdcaeffCF->GetXaxis()->SetTitle("dca (#mum)");
	hdcaeffCF->GetYaxis()->SetTitle("Efficiency");
	hdcaeffCF->Draw("hist");
	hdcaeffCF->Draw("err same");

	ceff->cd(8);
	TH1D *hd0pieffCF = eff->Project(id0pi); //the efficiency vs d0pi
	//hd0pieffCF->Setw2();
	//hd0pieffCF->SetMinimum(0.01);
	hd0pieffCF->SetLineColor(8);
	hd0pieffCF->SetLineWidth(3);
	hd0pieffCF->SetMarkerColor(8);
	hd0pieffCF->SetMarkerStyle(20);
	hd0pieffCF->GetXaxis()->SetTitleOffset(1.2);
	hd0pieffCF->GetYaxis()->SetTitleOffset(1.5);
	hd0pieffCF->Draw("hist");
	ceffd0pi->cd();
	ceffd0pi->SetLeftMargin(0.15);
	ceffd0pi->SetRightMargin(0.05);
	hd0pieffCF->GetXaxis()->SetTitle("d0_{#pi} (#mum)");
	hd0pieffCF->GetYaxis()->SetTitle("Efficiency");
	hd0pieffCF->Draw("hist");
	hd0pieffCF->Draw("err same");

	ceff->cd(9);
	TH1D *hd0KeffCF = eff->Project(id0K); //the efficiency vs d0K
	//hd0KeffCF->Setw2();
	//hd0KeffCF->SetMinimum(0.01);
	hd0KeffCF->SetLineColor(8);
	hd0KeffCF->SetLineWidth(3);
	hd0KeffCF->SetMarkerColor(8);
	hd0KeffCF->SetMarkerStyle(20);
	hd0KeffCF->GetXaxis()->SetTitleOffset(1.2);
	hd0KeffCF->GetYaxis()->SetTitleOffset(1.5);
	hd0KeffCF->Draw("hist");
	ceffd0K->cd();
	ceffd0K->SetLeftMargin(0.15);
	ceffd0K->SetRightMargin(0.05);
	hd0KeffCF->GetXaxis()->SetTitle("d0_{K} (#mum)");
	hd0KeffCF->GetYaxis()->SetTitle("Efficiency");
	hd0KeffCF->Draw("hist");
	hd0KeffCF->Draw("err same");

	ceff->cd(10);
	TH1D *hd0xd0effCF = eff->Project(id0xd0); //the efficiency vs d0xd0
	//hd0xd0effCF->Setw2();
	//hd0xd0effCF->SetMinimum(0.01);
	hd0xd0effCF->SetLineColor(8);
	hd0xd0effCF->SetLineWidth(3);
	hd0xd0effCF->SetMarkerColor(8);
	hd0xd0effCF->SetMarkerStyle(20);
	hd0xd0effCF->GetXaxis()->SetTitleOffset(1.2);
	hd0xd0effCF->GetYaxis()->SetTitleOffset(1.5);
	hd0xd0effCF->Draw("hist");
	ceffd0xd0->cd();
	ceffd0xd0->SetLeftMargin(0.15);
	ceffd0xd0->SetRightMargin(0.05);
	hd0xd0effCF->GetXaxis()->SetTitle("d0_{#pi}xd0_{K} (#mum^2)");
	hd0xd0effCF->GetYaxis()->SetTitle("Efficiency");
	hd0xd0effCF->Draw("hist");
	hd0xd0effCF->Draw("err same");

	ceff->cd(11);
	TH1D *hpointingeffCF = eff->Project(ipointing); //the efficiency vs pointing
	//hpointingeffCF->Setw2();
	//	hpointingeffCF->SetMinimum(0.01);
	hpointingeffCF->SetLineColor(8);
	hpointingeffCF->SetLineWidth(3);
	hpointingeffCF->SetMarkerColor(8);
	hpointingeffCF->SetMarkerStyle(20);
	hpointingeffCF->GetXaxis()->SetTitleOffset(1.2);
	hpointingeffCF->GetYaxis()->SetTitleOffset(1.5);
	hpointingeffCF->Draw("hist");
	ceffpointing->cd();
	ceffpointing->SetLeftMargin(0.15);
	ceffpointing->SetRightMargin(0.05);
	hpointingeffCF->GetXaxis()->SetTitle("cosPointingAngle");
	hpointingeffCF->GetYaxis()->SetTitle("Efficiency");
	hpointingeffCF->Draw("hist");
	hpointingeffCF->Draw("err same");

	// making slice

	// define vars on which to make the slice
	Int_t nvarSlice = 11;
	Int_t* ivarSlice = new Int_t[nvarSlice];
	ivarSlice[0] = ipt;
	ivarSlice[1] = iy;
	ivarSlice[2] = icTS;
	ivarSlice[3] = iptPi;
	ivarSlice[4] = iptK;
	ivarSlice[5] = icT;
	ivarSlice[6] = idca;
	ivarSlice[7] = id0pi;
	ivarSlice[8] = id0K;
	ivarSlice[9] = id0xd0;
	ivarSlice[10] = ipointing;

	// define limits for all vars
	Double_t mins[11];
	Double_t maxs[11];
	// pt
	mins[0] = 2;
	maxs[0] = 10;
	// y
	mins[1] = -2.1;
	maxs[1] = 2.1; 
	// cTS
	mins[2] = -1.05;
	maxs[2] = 1.05; 
	// ptPi
	mins[3] = 0;
	maxs[3] = 10; 
	// ptK
	mins[4] = 0;
	maxs[4] = 10; 
	// cT
	mins[5] = 0;
	maxs[5] = 500; 
	// dca
	mins[6] = 0;
	maxs[6] = 500; 
	// d0pi
	mins[7] = -1000;
	maxs[7] = 1000; 
	// d0K
	mins[8] = -1000;
	maxs[8] = 1000; 
	//d0xd0
	mins[9] = -100000;
	maxs[9] = 100000; 
	// pointing
	mins[10] = 0.5;
	maxs[10] = 1.05; 

	AliCFContainer* slicedCont = data->MakeSlice(nvarSlice, ivarSlice, mins, maxs);

	cout << "the new container has " << slicedCont->GetNStep() << " steps " << endl;

	TFile* fSlice = new TFile("fileSliceCont_pt_gt_2_pointing_gt_05.root","RECREATE");
	slicedCont->Write("container");
	fSlice->Close();

	//construct the efficiency grid from the sliced data container 
	AliCFEffGrid *effSlice = new AliCFEffGrid("effSlice"," The Sliced efficiency",*slicedCont);
	effSlice->CalculateEfficiency(stepNum,stepDen); //eff= step1/step0

	//The efficiency along the variables, and some 2-D projections
	TCanvas *ceffSlice =new TCanvas("ceffSlice"," Sliced Efficiency",0,0,1600,1200);
	ceffSlice->Divide(4,3);
	TCanvas *ceffSlicept = new TCanvas("ceffSlicept","Sliced Efficiency vs pt",50,50,550,550);
	TCanvas *ceffSlicey = new TCanvas("ceffSlicey","Sliced Efficiency vs y",50,50,550,550);
	TCanvas *ceffSlicecTS = new TCanvas("ceffSlicecTS","Sliced Efficiency vs cosThetaStar",50,50,550,550);
	TCanvas *ceffSliceptPi = new TCanvas("ceffSliceptPi","Sliced Efficiency vs ptPi",50,50,550,550);
	TCanvas *ceffSliceptK = new TCanvas("ceffSliceptK","Sliced Efficiency vs ptK",50,50,550,550);
	TCanvas *ceffSlicecT = new TCanvas("ceffSlicecT","Sliced Efficiency vs cT",50,50,550,550);
	TCanvas *ceffSlicedca = new TCanvas("ceffSlicedca","Sliced Efficiency vs dca",50,50,550,550);
	TCanvas *ceffSliced0pi = new TCanvas("ceffSliced0pi","Sliced Efficiency vs d0pi",50,50,550,550);
	TCanvas *ceffSliced0K = new TCanvas("ceffSliced0K","Sliced Efficiency vs d0K",50,50,550,550);
	TCanvas *ceffSliced0xd0 = new TCanvas("ceffSliced0xd0","Sliced Efficiency vs d0xd0",50,50,550,550);
	TCanvas *ceffSlicepointing = new TCanvas("ceffSlicepointing","Sliced Efficiency vs pointing",50,50,550,550);

	ceffSlice->cd(1);
	TH1D *hpteffSliceCF = effSlice->Project(ipt); //the efficiency vs pt
	//hpteffCF->Setw2();
	//hpteffSliceCF->SetMinimum(0.01);
	hpteffSliceCF->SetLineColor(8);
	hpteffSliceCF->SetLineWidth(3);
	hpteffSliceCF->SetMarkerColor(8);
	hpteffSliceCF->SetMarkerStyle(20);
	hpteffSliceCF->GetXaxis()->SetTitleOffset(1.2);
	hpteffSliceCF->GetYaxis()->SetTitleOffset(1.5);
	hpteffSliceCF->Draw("hist");
	ceffSlicept->cd();
	ceffSlicept->SetLeftMargin(0.15);
	ceffSlicept->SetRightMargin(0.05);
	hpteffSliceCF->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hpteffSliceCF->GetYaxis()->SetTitle("Sliced Efficiency");
	hpteffSliceCF->Draw("hist");
	hpteffSliceCF->Draw("err same");

	ceffSlice->cd(2);
	TH1D *hyeffSliceCF = effSlice->Project(iy); //the Sliced efficiency vs y
	//hyeffSliceCF->Setw2();
	//hyeffSliceCF->SetMinimum(0.01);
	hyeffSliceCF->SetLineColor(8);
	hyeffSliceCF->SetLineWidth(3);
	hyeffSliceCF->SetMarkerColor(8);
	hyeffSliceCF->SetMarkerStyle(20);
	hyeffSliceCF->GetXaxis()->SetTitleOffset(1.2);
	hyeffSliceCF->GetYaxis()->SetTitleOffset(1.5);
	hyeffSliceCF->Draw("hist");
	ceffSlicey->cd();
	ceffSlicey->SetLeftMargin(0.15);
	ceffSlicey->SetRightMargin(0.05);
	hyeffSliceCF->GetXaxis()->SetTitle("y");
	hyeffSliceCF->GetYaxis()->SetTitle("Sliced Efficiency");
	hyeffSliceCF->Draw("hist");
	hyeffSliceCF->Draw("err same");

	ceffSlice->cd(3);
	TH1D *hcTSeffSliceCF = effSlice->Project(icTS); //the Sliced efficiency vs cosThetaStar
	//hcTSeffSliceCF->Setw2();
	//hcTSeffSliceCF->SetMinimum(0.01);
	hcTSeffSliceCF->SetLineColor(8);
	hcTSeffSliceCF->SetLineWidth(3);
	hcTSeffSliceCF->SetMarkerColor(8);
	hcTSeffSliceCF->SetMarkerStyle(20);
	hcTSeffSliceCF->GetXaxis()->SetTitleOffset(1.2);
	hcTSeffSliceCF->GetYaxis()->SetTitleOffset(1.5);
	hcTSeffSliceCF->Draw("hist");
	ceffSlicecTS->cd();
	ceffSlicecTS->SetLeftMargin(0.15);
	ceffSlicecTS->SetRightMargin(0.05);
	hcTSeffSliceCF->GetXaxis()->SetTitle("cosThetaStar");
	hcTSeffSliceCF->GetYaxis()->SetTitle("EffSliceiciency");
	hcTSeffSliceCF->Draw("hist");
	hcTSeffSliceCF->Draw("err same");

	ceffSlice->cd(4);
	TH1D *hptPieffSliceCF = effSlice->Project(iptPi); //the Sliced efficiency vs ptPi
	//hptPieffSliceCF->Setw2();
	//hptPieffSliceCF->SetMinimum(0.01);
	hptPieffSliceCF->SetLineColor(8);
	hptPieffSliceCF->SetLineWidth(3);
	hptPieffSliceCF->SetMarkerColor(8);
	hptPieffSliceCF->SetMarkerStyle(20);
	hptPieffSliceCF->GetXaxis()->SetTitleOffset(1.2);
	hptPieffSliceCF->GetYaxis()->SetTitleOffset(1.5);
	hptPieffSliceCF->Draw("hist");
	ceffSliceptPi->cd();
	ceffSliceptPi->SetLeftMargin(0.15);
	ceffSliceptPi->SetRightMargin(0.05);
	hptPieffSliceCF->GetXaxis()->SetTitle("p_{T, #pi} (GeV/c)");
	hptPieffSliceCF->GetYaxis()->SetTitle("Sliced Efficiency");
	hptPieffSliceCF->Draw("hist");
	hptPieffSliceCF->Draw("err same");

	ceffSlice->cd(5);
	TH1D *hptKeffSliceCF = effSlice->Project(iptK); //the Sliced efficiency vs ptK
	//hptKeffSliceCF->Setw2();
	//hptKeffSliceCF->SetMinimum(0.01);
	hptKeffSliceCF->SetLineColor(8);
	hptKeffSliceCF->SetLineWidth(3);
	hptKeffSliceCF->SetMarkerColor(8);
	hptKeffSliceCF->SetMarkerStyle(20);
	hptKeffSliceCF->GetXaxis()->SetTitleOffset(1.2);
	hptKeffSliceCF->GetYaxis()->SetTitleOffset(1.5);
	hptKeffSliceCF->Draw("hist");
	ceffSliceptK->cd();
	ceffSliceptK->SetLeftMargin(0.15);
	ceffSliceptK->SetRightMargin(0.05);
	hptKeffSliceCF->GetXaxis()->SetTitle("p_{T, K} (GeV/c)");
	hptKeffSliceCF->GetYaxis()->SetTitle("Sliced Efficiency");
	hptKeffSliceCF->Draw("hist");
	hptKeffSliceCF->Draw("err same");

	ceffSlice->cd(6);
	TH1D *hcTeffSliceCF = effSlice->Project(icT); //the Sliced efficiency vs cT
	//hcTeffSliceCF->Setw2();
	//hcTeffSliceCF->SetMinimum(0.01);
	hcTeffSliceCF->SetLineColor(8);
	hcTeffSliceCF->SetLineWidth(3);
	hcTeffSliceCF->SetMarkerColor(8);
	hcTeffSliceCF->SetMarkerStyle(20);
	hcTeffSliceCF->GetXaxis()->SetTitleOffset(1.2);
	hcTeffSliceCF->GetYaxis()->SetTitleOffset(1.5);
	hcTeffSliceCF->Draw("hist");
	ceffSlicecT->cd();
	ceffSlicecT->SetLeftMargin(0.15);
	ceffSlicecT->SetRightMargin(0.05);
	hcTeffSliceCF->GetXaxis()->SetTitle("ct (#mum)");
	hcTeffSliceCF->GetYaxis()->SetTitle("Sliced Efficiency");
	hcTeffSliceCF->Draw("hist");
	hcTeffSliceCF->Draw("err same");

	ceffSlice->cd(7);
	TH1D *hdcaeffSliceCF = effSlice->Project(idca); //the Sliced efficiency vs dca
	//hdcaeffSliceCF->Setw2();
	//hdcaeffSliceCF->SetMinimum(0.01);
	hdcaeffSliceCF->SetLineColor(8);
	hdcaeffSliceCF->SetLineWidth(3);
	hdcaeffSliceCF->SetMarkerColor(8);
	hdcaeffSliceCF->SetMarkerStyle(20);
	hdcaeffSliceCF->GetXaxis()->SetTitleOffset(1.2);
	hdcaeffSliceCF->GetYaxis()->SetTitleOffset(1.5);
	hdcaeffSliceCF->Draw("hist");
	ceffSlicedca->cd();
	ceffSlicedca->SetLeftMargin(0.15);
	ceffSlicedca->SetRightMargin(0.05);
	hdcaeffSliceCF->GetXaxis()->SetTitle("dca (#mum)");
	hdcaeffSliceCF->GetYaxis()->SetTitle("Sliced Efficiency");
	hdcaeffSliceCF->Draw("hist");
	hdcaeffSliceCF->Draw("err same");

	ceffSlice->cd(8);
	TH1D *hd0pieffSliceCF = effSlice->Project(id0pi); //the Sliced efficiency vs d0pi
	//hd0pieffSliceCF->Setw2();
	//hd0pieffSliceCF->SetMinimum(0.01);
	hd0pieffSliceCF->SetLineColor(8);
	hd0pieffSliceCF->SetLineWidth(3);
	hd0pieffSliceCF->SetMarkerColor(8);
	hd0pieffSliceCF->SetMarkerStyle(20);
	hd0pieffSliceCF->GetXaxis()->SetTitleOffset(1.2);
	hd0pieffSliceCF->GetYaxis()->SetTitleOffset(1.5);
	hd0pieffSliceCF->Draw("hist");
	ceffSliced0pi->cd();
	ceffSliced0pi->SetLeftMargin(0.15);
	ceffSliced0pi->SetRightMargin(0.05);
	hd0pieffSliceCF->GetXaxis()->SetTitle("d0_{#pi} (#mum)");
	hd0pieffSliceCF->GetYaxis()->SetTitle("Efficiency");
	hd0pieffSliceCF->Draw("hist");
	hd0pieffSliceCF->Draw("err same");

	ceffSlice->cd(9);
	TH1D *hd0KeffSliceCF = effSlice->Project(id0K); //the Sliced efficiency vs d0K
	//hd0KeffSliceCF->Setw2();
	//hd0KeffSliceCF->SetMinimum(0.01);
	hd0KeffSliceCF->SetLineColor(8);
	hd0KeffSliceCF->SetLineWidth(3);
	hd0KeffSliceCF->SetMarkerColor(8);
	hd0KeffSliceCF->SetMarkerStyle(20);
	hd0KeffSliceCF->GetXaxis()->SetTitleOffset(1.2);
	hd0KeffSliceCF->GetYaxis()->SetTitleOffset(1.5);
	hd0KeffSliceCF->Draw("hist");
	ceffSliced0K->cd();
	ceffSliced0K->SetLeftMargin(0.15);
	ceffSliced0K->SetRightMargin(0.05);
	hd0KeffSliceCF->GetXaxis()->SetTitle("d0_{K} (#mum)");
	hd0KeffSliceCF->GetYaxis()->SetTitle("Sliced Efficiency");
	hd0KeffSliceCF->Draw("hist");
	hd0KeffSliceCF->Draw("err same");

	ceffSlice->cd(10);
	TH1D *hd0xd0effSliceCF = effSlice->Project(id0xd0); //the Sliced efficiency vs d0xd0
	//hd0xd0effSliceCF->Setw2();
	//hd0xd0effSliceCF->SetMinimum(0.01);
	hd0xd0effSliceCF->SetLineColor(8);
	hd0xd0effSliceCF->SetLineWidth(3);
	hd0xd0effSliceCF->SetMarkerColor(8);
	hd0xd0effSliceCF->SetMarkerStyle(20);
	hd0xd0effSliceCF->GetXaxis()->SetTitleOffset(1.2);
	hd0xd0effSliceCF->GetYaxis()->SetTitleOffset(1.5);
	hd0xd0effSliceCF->Draw("hist");
	ceffSliced0xd0->cd();
	ceffSliced0xd0->SetLeftMargin(0.15);
	ceffSliced0xd0->SetRightMargin(0.05);
	hd0xd0effSliceCF->GetXaxis()->SetTitle("d0_{#pi}xd0_{K} (#mum^2)");
	hd0xd0effSliceCF->GetYaxis()->SetTitle("Sliced Efficiency");
	hd0xd0effSliceCF->Draw("hist");
	hd0xd0effSliceCF->Draw("err same");

	ceffSlice->cd(11);
	TH1D *hpointingeffSliceCF = effSlice->Project(ipointing); //the efficiency vs pointing
	//hpointingeffSliceCF->Setw2();
	//	hpointingeffSliceCF->SetMinimum(0.01);
	hpointingeffSliceCF->SetLineColor(8);
	hpointingeffSliceCF->SetLineWidth(3);
	hpointingeffSliceCF->SetMarkerColor(8);
	hpointingeffSliceCF->SetMarkerStyle(20);
	hpointingeffSliceCF->GetXaxis()->SetTitleOffset(1.2);
	hpointingeffSliceCF->GetYaxis()->SetTitleOffset(1.5);
	hpointingeffSliceCF->Draw("hist");
	ceffSlicepointing->cd();
	ceffSlicepointing->SetLeftMargin(0.15);
	ceffSlicepointing->SetRightMargin(0.05);
	hpointingeffSliceCF->GetXaxis()->SetTitle("cosPointingAngle");
	hpointingeffSliceCF->GetYaxis()->SetTitle("Sliced Efficiency");
	hpointingeffSliceCF->Draw("hist");
	hpointingeffSliceCF->Draw("err same");


	TFile* fileEff = new TFile("fileEff_GenAcc_RecPPR_pt_gt_2_pointing_gt_05.root", "RECREATE");
	hpteffSliceCF->Write("hpteffCFGenAccRecPPR");
	hyeffSliceCF->Write("hyeffCFGenAccRecPPR");
	hcTSeffSliceCF->Write("hcTSeffCFGenAccRecPPR");
	hptPieffSliceCF->Write("hptPieffCFGenAccRecPPR");
	hptKeffSliceCF->Write("hptKeffCFGenAccRecPPR");
	hcTeffSliceCF->Write("hcTeffCFGenAccRecPPR");
	hdcaeffSliceCF->Write("hdcaeffCFGenAccRecPPR");
	hd0pieffSliceCF->Write("hd0pieffCFGenAccRecPPR");
	hd0KeffSliceCF->Write("hd0KeffCFGenAccRecPPR");
	hd0xd0effSliceCF->Write("hd0xd0effCFGenAccRecPPR");
	hpointingeffSliceCF->Write("hpointingeffCFGenAccRecPPR");
	fileEff->Close();
	

	TString dir(plotsDir);
        if (gSystem->Exec(Form("ls %s",dir.Data()))!=0){
                cout << " Creating directory for plots" << endl;

                gSystem->mkdir(Form("%s_pt_gt_2_pointing_gt_05",dir.Data()));
		gSystem->cd(Form("%s_pt_gt_2_pointing_gt_05",dir.Data()));
        }

        // printing eps files
        ceffSlice->Print(Form("efficiencies.eps", dir.Data()));
        ceffSlicept->Print(Form("effpt.eps", dir.Data()));
        ceffSlicey->Print(Form("effy.eps", dir.Data()));
        ceffSlicecTS->Print(Form("effcTS.eps", dir.Data()));
        ceffSliceptPi->Print(Form("effptPi.eps", dir.Data()));
        ceffSliceptK->Print(Form("effptK.eps", dir.Data()));
        ceffSlicecT->Print(Form("effcT.eps", dir.Data()));
        ceffSlicedca->Print(Form("effdca.eps", dir.Data()));
        ceffSliced0pi->Print(Form("effd0pi.eps", dir.Data()));
        ceffSliced0K->Print(Form("effd0K.eps", dir.Data()));
        ceffSliced0xd0->Print(Form("effd0xd0.eps", dir.Data()));
        ceffSlicepointing->Print(Form("effpointing.eps", dir.Data()));

	// printing gif files
        ceffSlice->Print(Form("efficiencies.gif", dir.Data()));
        ceffSlicept->Print(Form("effpt.gif", dir.Data()));
        ceffSlicey->Print(Form("effy.gif", dir.Data()));
        ceffSlicecTS->Print(Form("effcTS.gif", dir.Data()));
        ceffSliceptPi->Print(Form("effptPi.gif", dir.Data()));
        ceffSliceptK->Print(Form("effptK.gif", dir.Data()));
        ceffSlicecT->Print(Form("effcT.gif", dir.Data()));
        ceffSlicedca->Print(Form("effdca.gif", dir.Data()));
        ceffSliced0pi->Print(Form("effd0pi.gif", dir.Data()));
        ceffSliced0K->Print(Form("effd0K.gif", dir.Data()));
        ceffSliced0xd0->Print(Form("effd0xd0.gif", dir.Data()));
        ceffSlicepointing->Print(Form("effpointing.gif", dir.Data()));


}
