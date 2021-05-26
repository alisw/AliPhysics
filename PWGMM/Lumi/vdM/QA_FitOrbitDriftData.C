#include "CandRoot.h"

//-------------------------------------------------------
// Macro read the TCanvas we got from Martino to fit the orbit drift data
//-------------------------------------------------------

//-------------------------------------------------------
void ProcessFill6012(Int_t Fill)
{
	// time borders of individual scans
	// obtained using QA_print_timestamps.C
	Long_t time_start_x0 = 1501171527;
	Long_t time_end_x0   = 1501172628;  
	Long_t time_start_y0 = 1501172798;
	Long_t time_end_y0   = 1501173903;  
	Long_t time_start_x1 = 1501174067;
	Long_t time_end_x1   = 1501175170;  
	Long_t time_start_y1 = 1501175390;
	Long_t time_end_y1   = 1501176496;

	// open the file with the canvas and get the graphs
	TFile *Input = new TFile("../Fill-6012/OrbitDriftCorrection_6012.root");
	TCanvas *canvas = (TCanvas *) Input->Get("c1_n3");
	TList *ObjectsInCanvas = canvas->GetListOfPrimitives();
	// this prints the objects in the canvas
	// ObjectsInCanvas->ls();
	TGraph *gr_hor = (TGraph*) ObjectsInCanvas->At(2);
	TGraph *gr_ver = (TGraph*) ObjectsInCanvas->At(3);

	// loop over the graphs and and create new graph
	// with points outside the scan time
	// assuming same time coordinates in both graphs ...
	// --> variables to read the info from graph
	Double_t time;
	Double_t drift_hor;
	Double_t drift_ver;
	// --> new graphs to stored selected infor
	TGraph *gr_hor_new = new TGraph();
	TGraph *gr_ver_new = new TGraph();
	// --> variables to define the drawing frame
	Double_t min_time = 1e12;
	Double_t max_time = 0;  
	Double_t min_drift = 1e12;
	Double_t max_drift = -1e12;  
	for(Int_t i=0;i<gr_hor->GetN();i++)
	{
		// get the points
		gr_hor->GetPoint(i,time,drift_hor);
		gr_ver->GetPoint(i,time,drift_ver);
		Long_t timeL = (Long_t) time;

		//cout <<Form("Fill %i, index %5i, drift_hor: %8.4f, drift_ver: %8.4f\n", Fill, i, drift_hor, drift_ver);
		if (TMath::IsNaN(drift_hor)) continue;
		if (TMath::IsNaN(drift_ver)) continue;

		// update min/max
		if (time<min_time) min_time = time;
		if (time>max_time) max_time = time;    
		if (drift_hor<min_drift) min_drift = drift_hor;
		if (drift_hor>max_drift) max_drift = drift_hor;    
		if (drift_ver<min_drift) min_drift = drift_ver;
		if (drift_ver>max_drift) max_drift = drift_ver;
		// apply the selection
		if (timeL>time_start_x0 && timeL<time_end_x0) continue;
		if (timeL>time_start_y0 && timeL<time_end_y0) continue;     
		if (timeL>time_start_x1 && timeL<time_end_x1) continue;
		if (timeL>time_start_y1 && timeL<time_end_y1) continue;
		if (timeL>(time_end_y1+200)) continue;
		// fill new graph
		gr_hor_new->SetPoint(gr_hor_new->GetN(),timeL,drift_hor);
		gr_ver_new->SetPoint(gr_ver_new->GetN(),timeL,drift_ver);
	}

	TCanvas* c1 = new TCanvas("c1", "ODC fit", 1600, 600*3); c1->Divide(1, 3);

	// draw data
	min_drift*=1.5;
	max_drift*=1.5;
	//TCanvas *cgr = new TCanvas("cgr","cgr",0,0,900,600);
	c1->cd(1);
	TH1F* fr = gPad->DrawFrame(min_time-100, min_drift,max_time+100,max_drift);
	fr->SetTitle(Form("Fill %i;time (s);drift (#mum)", Fill));
	// original data
	gr_hor->SetMarkerStyle(20);gr_hor->SetMarkerColor(2);
	gr_hor->SetMarkerSize(0.2);
	gr_hor->Draw("p,same");
	gr_ver->SetMarkerStyle(20);gr_ver->SetMarkerColor(4);  
	gr_ver->SetMarkerSize(0.2);
	gr_ver->Draw("p,same");
	// selected data
	gr_hor_new->SetMarkerStyle(20);gr_hor_new->SetMarkerColor(2);
	gr_hor_new->Draw("p,same");
	gr_ver_new->SetMarkerStyle(20);gr_ver_new->SetMarkerColor(4);  
	gr_ver_new->Draw("p,same");

	// fit horizontal
	//  TF1 *f_hor = new TF1("f_hor","[1]+[2]*(x-[0])+[3]*pow(x-[0],2)+[4]*pow(x-[0],3)+[5]*pow(x-[0],4)",
	TF1 *f_hor = new TF1("f_hor","[1]+[2]*(x-[0])+[3]*pow(x-[0],3)+[4]*pow(x-[0],4)", min_time,time_end_y1+200);
	f_hor->SetLineColor(2);
	f_hor->SetParameter(0, time_start_x1);
	f_hor->SetParLimits(0,min_time,time_end_y1+200);
	f_hor->SetParameter(1,0);
	f_hor->SetParLimits(1,min_drift,max_drift);  
	f_hor->SetParameter(2,1);
	gr_hor_new->Fit("f_hor");

	// plot horizontal fit
	//TCanvas *ch = new TCanvas("ch","ch",0,0,900,600);
	c1->cd(2);
	TH1F* frh = gPad->DrawFrame(min_time-100,min_drift, time_end_y1+200,max_drift);
	frh->SetTitle("Horizontal;time (s);drift (#mum)");
	gr_hor->Draw("p,same");
	gr_hor_new->Draw("p,same");
	cout << " hor p0 = " << ((Long_t) f_hor->GetParameter(0)) << endl;

	// fit vertical
	TF1 *f_ver = new TF1("f_ver","[1]+[2]*(x-[0])",min_time,time_end_y1+200);
	f_ver->SetLineColor(4);
	f_ver->SetParameter(0, time_start_x1);
	f_ver->SetParLimits(0,min_time,time_end_y1+200);
	f_ver->SetParameter(1,0);
	f_ver->SetParLimits(1,min_drift,max_drift);  
	f_ver->SetParameter(2,1);
	gr_ver_new->Fit("f_ver");

	// plot vertical fit
	gStyle->SetOptFit(1);
	//TCanvas *cv = new TCanvas("cv","cv",0,0,900,600);
	c1->cd(3);
	TH1F* frv = gPad->DrawFrame(min_time-100,min_drift,time_end_y1+200,max_drift);
	frv->SetTitle("Vertical;time (s);drift (#mum)");
	gr_ver->Draw("p,same");
	gr_ver_new->Draw("p,same");
	cout << " ver p0 = " << ((Long_t) f_ver->GetParameter(0)) << endl;

	c1->Print(Form("c1d_ODC_Fill%i.png", Fill)); //kimc
	return;
}


//-------------------------------------------------------
void ProcessFill6864(Int_t Fill)
{
	// time borders of individual scans
	// obtained using QA_print_timestamps.C
	Long_t time_start_x0 = 1530314417;
	Long_t time_end_x0   = 1530315529;  
	Long_t time_start_y0 = 1530315682;
	Long_t time_end_y0   = 1530316794;  
	Long_t time_start_x1 = 1530316972;
	Long_t time_end_x1   = 1530318083;  
	Long_t time_start_y1 = 1530318233;
	Long_t time_end_y1   = 1530319344;

	// open the file with the canvas and get the graphs
	TFile *Input = new TFile("../Fill-6864/OrbitDriftCorrection_6864.root");
	TCanvas *canvas = (TCanvas *) Input->Get("c1_n3");
	TList *ObjectsInCanvas = canvas->GetListOfPrimitives();
	// this prints the objects in the canvas
	// ObjectsInCanvas->ls();
	TGraph *gr_hor = (TGraph*) ObjectsInCanvas->At(2);
	TGraph *gr_ver = (TGraph*) ObjectsInCanvas->At(3);

	// loop over the graphs and and create new graph
	// with points outside the scan time
	// assuming same time coordinates in both graphs ...
	// --> variables to read the info from graph
	Double_t time;
	Double_t drift_hor;
	Double_t drift_ver;
	// --> new graphs to stored selected infor
	TGraph *gr_hor_new = new TGraph();
	TGraph *gr_ver_new = new TGraph();
	// --> variables to define the drawing frame
	Double_t min_time = 1e12;
	Double_t max_time = 0;  
	Double_t min_drift = 1e12;
	Double_t max_drift = -1e12;  
	for(Int_t i=0;i<gr_hor->GetN();i++)
	{
		// get the points
		gr_hor->GetPoint(i,time,drift_hor);
		gr_ver->GetPoint(i,time,drift_ver);    
		Long_t timeL = (Long_t) time;
		// update min/max
		if (time<min_time) min_time = time;
		if (time>max_time) max_time = time;    
		if (drift_hor<min_drift) min_drift = drift_hor;
		if (drift_hor>max_drift) max_drift = drift_hor;    
		if (drift_ver<min_drift) min_drift = drift_ver;
		if (drift_ver>max_drift) max_drift = drift_ver;
		// apply the selection
		if (timeL>time_start_x0 && timeL<time_end_x0) continue;
		if (timeL>time_start_y0 && timeL<time_end_y0) continue;     
		if (timeL>time_start_x1 && timeL<time_end_x1) continue;
		if (timeL>time_start_y1 && timeL<time_end_y1) continue;
		if (timeL>(time_end_y1+200)) continue;
		// check for NaN
		if (TMath::IsNaN(drift_hor)) continue;
		if (TMath::IsNaN(drift_ver)) continue;    
		// fill new graph
		gr_hor_new->SetPoint(gr_hor_new->GetN(),timeL,drift_hor);
		gr_ver_new->SetPoint(gr_ver_new->GetN(),timeL,drift_ver);
		//  cout << timeL << " " << drift_hor << " " << drift_ver << endl;
	}

	TCanvas* c1 = new TCanvas("c1", "ODC fit", 1600, 600*3); c1->Divide(1, 3);

	// draw data
	min_drift*=1.5;
	max_drift*=1.5;
	//TCanvas *cgr = new TCanvas("cgr","cgr",0,0,900,600); //kimc
	c1->cd(1);
	TH1F* fr = gPad->DrawFrame(min_time-100,min_drift, max_time+100,max_drift);
	fr->SetTitle(Form("Fill %i;time (s);drift (#mum)", Fill));
	// original data
	gr_hor->SetMarkerStyle(20);gr_hor->SetMarkerColor(2);
	gr_hor->SetMarkerSize(0.2);
	gr_hor->Draw("p,same");
	gr_ver->SetMarkerStyle(20);gr_ver->SetMarkerColor(4);  
	gr_ver->SetMarkerSize(0.2);
	gr_ver->Draw("p,same");
	// selected data
	gr_hor_new->SetMarkerStyle(20);gr_hor_new->SetMarkerColor(2);
	gr_hor_new->Draw("p,same");
	gr_ver_new->SetMarkerStyle(20);gr_ver_new->SetMarkerColor(4);  
	gr_ver_new->Draw("p,same");

	// fit horizontal
	TF1 *f_hor = new TF1("f_hor","[1]+[2]*(x-[0])+[3]*pow(x-[0],2)+[4]*pow(x-[0],3)+[5]*pow(x-[0],4)",
			min_time,time_end_y1+200);
	f_hor->SetLineColor(2);
	f_hor->SetParameter(0, time_start_x1);
	f_hor->SetParLimits(0,min_time,time_end_y1+200);
	f_hor->SetParameter(1,0);
	f_hor->SetParLimits(1,min_drift,max_drift);  
	f_hor->SetParameter(2,1);
	gr_hor_new->Fit("f_hor");

	// plot horizontal fit
	//TCanvas *ch = new TCanvas("ch","ch",0,0,900,600);
	c1->cd(2);
	TH1F* frh = gPad->DrawFrame(min_time-100,min_drift, time_end_y1+200,max_drift);
	frh->SetTitle("Horizontal;time (s);drift (#mum)");
	gr_hor->Draw("p,same");
	gr_hor_new->Draw("p,same");
	cout << " hor p0 = " << ((Long_t) f_hor->GetParameter(0)) << endl;

	// fit vertical
	TF1 *f_ver = new TF1("f_ver","[1]+[2]*(x-[0])+[3]*pow(x-[0],2)", min_time,time_end_y1+200);
	f_ver->SetLineColor(4);
	f_ver->SetParameter(0, time_start_x1);
	f_ver->SetParLimits(0,min_time,time_end_y1+200);
	f_ver->SetParameter(1,2);
	// f_ver->SetParLimits(1,min_drift,max_drift);  
	f_ver->SetParameter(2,1);
	gr_ver_new->Fit("f_ver");

	// plot vertical fit
	gStyle->SetOptFit(1);
	//TCanvas *cv = new TCanvas("cv","cv",0,0,900,600); //kimc
	c1->cd(3);
	TH1F* frv = gPad->DrawFrame(min_time-100,min_drift,	time_end_y1+200,max_drift);
	frv->SetTitle("Vertical;time (s);drift (#mum)");
	gr_ver->Draw("p,same");
	gr_ver_new->Draw("p,same");
	cout << " ver p0 = " << ((Long_t) f_ver->GetParameter(0)) << endl;

	c1->Print(Form("c1d_ODC_Fill%i.png", Fill)); //kimc
	return;
}

//-------------------------------------------------------

void QA_FitOrbitDriftData(Int_t Fill)
{
	if (Fill == 6012) ProcessFill6012(Fill);
	if (Fill == 6864) ProcessFill6864(Fill);

	return;
}
