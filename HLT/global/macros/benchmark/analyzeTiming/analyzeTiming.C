#include <fstream>
#include <algorithm>
#include <string.h>

TCanvas * c;
TString folder;
TTree * t;
TString conversion[4] = {
	"esd",
	"flat",
	"esdToFlat",
	"flatVsEsd"
};
TString condition[3] = {
	"AliHLTGlobalEsdConverterComponent::DoEvent.Stop",
	"AliHLTGlobalFlatEsdConverterComponent::DoEvent.Stop",
	"AliHLTGlobalEsdToFlatConverterComponent::DoEvent.Stop"
};

// maximum values for plotting:
// inSize, outSize, time
Double_t maxs[4][3];



/**
 * Option mask:
 * 
 * 1: log
 * 2: fit
 * 4: bisector // Winkelhalbierende
 * 8: 2-dim
 * 
 * */


/*
 * 
 * values
 * 
 * id0 inSize
 * id1 outSize
 * id2 nTracks
 * id3 nV0s
 * 
 */




void createHistos( TString fo= "pp", bool individual=true, bool combined = true ){

	initialize(fo );
	

	
	
// create histograms for different types of conversion	
	
	
	if(individual){
		
		for(int i=0; i < 3 ; i++){

			// Draw input size distribution
			makeHist( i, "id0/1024>>h()", "inSize", "input size (kB)", "\# events", 1);
			// Draw output size distribution
			makeHist( i, "id1/1024>>h()", "outSize", "output size (kB)", "\# events", 1);
			// Draw track number distribution
			makeHist( i, "id2>>h()", "nTracks", "\# tracks", "\# events", 1);
			// Draw v0 number distribution
			makeHist(  i, "id3>>h()","nV0s", "\# V0s", "\# events", 1);
			// Draw clusters/track number distribution
		//	makeHist( i,  "id3/id2>>h()", "nClustersPerTrack", "\# clusters/track", "\# events", 1);
			// Draw output vs input size
			makeHist( i,  "id1/1024:id0/1024>>h()", "outVsIn", "input size (kB)", "output size (kB)", 2);

			createTimeHistos(i);

		}

		
		
	}

//  create combined histograms to compare conversions	

	if(combined){
		
		
		createTimeHistos(3);

		
		
		
	// loop over trees to get combined histos
		
		TH2D*flatVsNormalUserTime = new TH2D("userTime_flatVsEsd", "cpu time (user) : flat vs normal", 400,0,maxs[0][2], 150,0,maxs[1][2]);
		flatVsNormalUserTime->GetXaxis()->SetTitle("cpu time ESD converter (ms)");
		flatVsNormalUserTime->GetYaxis()->SetTitle("cpu time flatESD converter (ms)");
		
		TH2D*flatVsNormalInSize = new TH2D("inSize_flatVsEsd", "input size: flat vs normal", 200,0,maxs[0][0], 2000,0,maxs[1][0]);
		flatVsNormalInSize->GetXaxis()->SetTitle("input size ESD converter (kB)");
		flatVsNormalInSize->GetYaxis()->SetTitle("input size flatESD converter (kB)");
		
		TH2D*flatVsNormalOutSize = new TH2D("outSize_flatVsEsd", "output size: flat vs normal", 80,0,maxs[0][1], 800,0,maxs[1][1]);
		flatVsNormalOutSize->GetXaxis()->SetTitle("output size ESD converter (kB)");
		flatVsNormalOutSize->GetYaxis()->SetTitle("output size flatESD converter (kB)");
		
		Double_t cpu, cpuTmp, cpuOld, cpuOldTmp;
		Int_t inSize, inSizeTmp, outSize, outSizeTmp;
		Char_t sname[100];
		Char_t snameTmp[100];

		t->SetBranchAddress("pI.fCpuUser",&cpuTmp);
		t->SetBranchAddress("pIOld.fCpuUser",&cpuOldTmp);
		t->SetBranchAddress("sname",&snameTmp);
		t->SetBranchAddress("id0",&inSizeTmp);
		t->SetBranchAddress("id1",&outSizeTmp);

		t->SetBranchStatus("*",0); //disable all branches
		t->SetBranchStatus("sname",1);
		t->SetBranchStatus("pI.fCpuUser",1);
		t->SetBranchStatus("pIOld.fCpuUser",1);
		t->SetBranchStatus("id0",1);
		t->SetBranchStatus("id1",1);

		// loop over entries in tree
		for (Int_t i1=0, i2=0; i1 < t->GetEntries()  ; i1++) {
			t->GetEntry(i1);
		// when found correct stamp for ESD converter in tree	
			if( condition[0].EqualTo(snameTmp) ){
				cpu = cpuTmp;
				cpuOld = cpuOldTmp;
				inSize = inSizeTmp;
				outSize = outSizeTmp;
				strcpy(sname, snameTmp);
				// loop again over entries in tree until correct stamp for flat converter
				
				for( t->GetEntry(i2++); strcmp(snameTmp, condition[1].Data() ) != 0 && i2 < t->GetEntries(); t->GetEntry(i2++) );
				flatVsNormalUserTime->Fill( 1000 * (cpu - cpuOld), 1000 * (cpuTmp - cpuOldTmp) );
				flatVsNormalInSize->Fill( inSize / 1024, inSizeTmp / 1024);
				flatVsNormalOutSize->Fill( outSize / 1024, outSizeTmp / 1024);
			}
		}

		saveHist(flatVsNormalUserTime, 2, "colz");
		saveHist(flatVsNormalInSize, 0, "colz");
		saveHist(flatVsNormalOutSize, 2, "colz");
		delete flatVsNormalUserTime;
		delete flatVsNormalInSize;
		delete flatVsNormalOutSize;

	}
	
	return;
}

//________________________________________________
//
//           Drawing functions
//________________________________________________



void createTimeHistos(Int_t ic){

	/* 
	* cpu time histos
	* loop over
	*	- user time
	*	- sys time
	*	- user time + sys time
	*/


	TString  plotString[3] = {
		"1000*(pI.fCpuUser-pIOld.fCpuUser)",
		"1000*(pI.fCpuSys-pIOld.fCpuSys)",
		"1000*(pI.fCpuUser+pI.fCpuSys-pIOld.fCpuUser-pIOld.fCpuSys)"
	};
	

	TString  axisString[3] = {
		"cpu time (user)/event(ms)",
		"cpu time (sys)/event(ms)",
		"cpu time (sys + user)/event(ms)"
	};
	TString  names[3] = {
		"userTime",
		"sysTime",
		"totalTime"
	};
	
	TString legends[2] = ic==3?  {"esd", "flat"}: 0;

	for(int j=0; j<3; ++j){
		
		
		// cpu time distribution		
		makeHist( ic, Form("%s>>h(200,0,%f)", plotString[j].Data(),maxs[ic][2]), Form("%sDistribution", names[j].Data() ), axisString[j], "\# events", 1, legends);
		
		// cpu time vs input size
		makeHist( ic, Form("%s:id0/1024>>h(100,0,%f,100,0,%f)",plotString[j].Data(), maxs[ic][0], maxs[ic][2]), Form("%sVsInSize", names[j].Data() ), "input size (kB)",  axisString[j], 8, legends );
					
		// cpu time vs output size
		makeHist( ic, Form("%s:id1/1024>>h(100,0,%f,100,0,%f)",plotString[j].Data(), maxs[ic][1],maxs[ic][2]), Form("%sVsOutSize", names[j].Data() ), "output size (kB)",  axisString[j], 8, legends);
					
		// cpu time vs nTracks
		makeHist( ic, Form("%s:id2>>h()",plotString[j].Data()), Form("%sVsTracks", names[j].Data()), "\# tracks", axisString[j], 8, legends);	
		
		// cpu time vs nV0s
		makeHist( ic, Form("%s:id3>>h()",plotString[j].Data()), Form("%sVsV0s", names[j].Data() ), "\# V0s", axisString[j], 8, legends);		
		
		// cpu time vs timestamp
		makeHist( ic, Form("%s:stampSec>>h()",plotString[j].Data()), Form("%sVsT", names[j].Data() ), "timestamp (s)", axisString[j], 8 , legends );	
		
		// cpu time vs time in file
		makeHist( ic, Form("%s:T>>h()",plotString[j].Data()), Form("%sVsTFile", names[j].Data() ), "time in file (s)",  axisString[j], 8, legends );		
		
	}
}



void saveHist(TH1*h, UInt_t optionsMask = 0, TString  options = ""){
	TCanvas *c = new TCanvas();
	if(optionsMask & 1 ) c->SetLogy();
	h->SetStats(0);
	h->SetMarkerStyle(7);
	h->Draw(options);
	
	
	
	c->SaveAs( Form("%s/png/%s.png", folder.Data(), h->GetName() ) );
	c->SaveAs( Form("%s/root/%s.root", folder.Data(), h->GetName()) );
	
	
	TLegend l= TLegend(.68,0.78,.98,0.95);
		l.SetFillColor(0);
		l.SetBorderSize(0);
	if(optionsMask & 2 ){
		TF1 linear = TF1("linear","pol1(0)", 0,2000);
		linear.SetParameters(0.5,0.2);
		linear.SetLineColor(kRed);
		linear.SetLineWidth(2);
		TFitResultPtr r = h->Fit(&linear, "S");
		l.AddEntry(&linear, Form("%.2g + %.2g * x", r->Parameter(0), r->Parameter(1) ) , "l");
		h->Fit(&linear);
		l.Draw("same");
		c->SaveAs( Form("%s/png/%s_fit.png", folder.Data(), h->GetName()) );
		h->SaveAs( Form("%s/root/%s_fit.root", folder.Data(), h->GetName()) );
	}
	if( optionsMask & 4  ){
		TF1 linear = TF1("lin","x",0,8);
		linear.SetLineColor(kRed);
		linear.SetLineWidth(2);
		linear.Draw("same");
		c->SaveAs( Form("%s/png/%s_bis.png", folder.Data(), h->GetName()) );
		h->SaveAs( Form("%s/root/%s_bis.root", folder.Data(), h->GetName()) );
	}
}


void makeHist(Int_t ic, TString plotString, TString name, TString x, TString y, UInt_t optionsMask = 0, TString *legends =0x0, TString options = ""){
	c = new TCanvas();
	if(optionsMask & 1) c->SetLogy();
	
	TH1* h = 0 ;
	TH1* i = 0;
	
	
	
	if(ic == 3){
		// flat vs esd
		t->Draw(plotString, Form("sname==\"%s\"", condition[0].Data()) );
		TH1* tmp = (TH1*) gDirectory->Get("h");
		h= (TH1*) tmp->Clone();
		t->Draw(plotString, Form("sname==\"%s\"", condition[1].Data()) );
		tmp = (TH1*) gDirectory->Get("h");
		i = (TH1*)tmp->Clone();
	}
	else{
		t->Draw(plotString,Form("sname==\"%s\"", condition[ic].Data()));
		h = (TH1*) gDirectory->Get("h");
	}
		
	h->SetTitle(Form("%s_%s", name.Data(), conversion[ic].Data()) );
	h->GetXaxis()->SetTitle(x);
	h->GetYaxis()->SetTitle(y);
	h->SetStats(0);
	
	TLegend l= TLegend(.58,0.68,.88,0.85);
		l.SetFillColor(0);
		l.SetBorderSize(0);
	if(i){
		if(optionsMask & 8){
			h->SetMarkerColor(kRed);
			i->SetMarkerColor(kBlue);
			h->SetMarkerStyle(21);
			i->SetMarkerStyle(21);
			h->SetMarkerSize(1);
			i->SetMarkerSize(1);
			if(legends){
				l.AddEntry(h, legends[0],"p");
				l.AddEntry(i, legends[1],"p");
			}
		}
		else{
			h->SetFillColor(kRed);
			i->SetFillColor(kBlue);
			if(legends){
				l.AddEntry(h, legends[0], "f");
				l.AddEntry(i, legends[1],"f");
			}
		}
	}
	
	if(optionsMask & 2){
		TF1 linear = TF1("linear","pol1(0)", 0,2000);
		linear.SetParameters(0.5,0.2);
		linear.SetLineColor(kRed);
		linear.SetLineWidth(2);
		TFitResultPtr r = h->Fit(&linear, "S");
		l.AddEntry(&linear, Form("%.2g + %.2g * x", r->Parameter(0), r->Parameter(1) ) , "l");
	}
	h->Draw(options);
	if(i) i->Draw("same");
	if(legends || optionsMask & 2) l.Draw("same");
	if( optionsMask & 4  ){
		TF1 linear = TF1("lin","x",0,8);
		linear.SetLineColor(kRed);
		linear.SetLineWidth(2);
		linear.Draw("same");
	}

	c->SaveAs( Form("%s/png/%s_%s.png", folder.Data(), name.Data(), conversion[ic].Data()) );
	h->SaveAs( Form("%s/root/%s_%s.root", folder.Data(), name.Data(), conversion[ic].Data()) );
	gDirectory->DeleteAll();
}


void initialize( TString fo ){
	folder = fo;
	t = AliSysInfo::MakeTree( Form("%s/syswatch_merged.log" , fo.Data() ) );
	
	ifstream axes(fo + "/axisRanges.txt");
	if( axes.good() ){
		for ( int i=0 ; i<3 && axes >> maxs[i][0] >> maxs[i][1] >> maxs[i][2]; ++i);
		axes.close();
	}
	else{
		for ( int i=0 ; i<3 ; ++i){
			maxs[i][0] = 2000;
			maxs[i][1] = 1000;
			maxs[i][2] = 20;
		}
		
	}	
	
	maxs[3][0] = std::max(maxs[0][0], maxs[1][0] );
	maxs[3][1] = std::max(maxs[0][1], maxs[1][1] );
	maxs[3][2] = std::max(maxs[0][2], maxs[1][2] );
}

