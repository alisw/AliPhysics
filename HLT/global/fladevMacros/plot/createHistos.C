TCanvas * c;


void createHistos(bool single=true, bool combined = true, bool merged=true, int fromType=0, int toType=2){


	c = new TCanvas();
// 0 = raw->esd
// 1 = raw->flat
// 2= esd -> flat


	TTree * t[4];


	char*sname[4]={"DoEvent.Stop","DoEvent.Stop","DoEvent.Stop","DoEvent.Stop"};
	char*conversion[4]={"rawToEsd","rawToFlat","HltEsdToFlat","NormalEsdToFlat"};
	char *infolder="rawToHltEsd";
	char *subfolder="";
	Int_t outMax = 80;
	Int_t inMax = 150;
	Int_t timeMax = 20;


// create histograms for different types of conversion	

for(int type=fromType; type<=toType; ++type){

	if(type==1){
		timeMax=6;
		outMax = 700;
		inMax = 2000;
		infolder = "rawToFlat";
	}
	else if(type==2){
		timeMax = 1;
		inMax =200;
		infolder = "HltEsdToFlat";
	}
	else if(type==3){
		timeMax = 1;
		inMax =200;
		infolder = "NormalEsdToFlat";
	}


	if(merged){
		t[type] = AliSysInfo::MakeTree(Form("%s/syswatch_merged.log" , infolder) );
		subfolder = "/merged";
	}
	else{
		t[type] = AliSysInfo::MakeTree(Form("%s/syswatch.log" , infolder) );
		subfolder = "/single";

	}



if (single){

	cout<<"Creating plots for "<< infolder<<endl;

	// Draw cpu, sys, cpu+sys time distribution

	t[type]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser)>>h(100,0,15)", Form("sname==\"%s\"", sname[type]) );
	saveHist("cpuTimeUserDistribution",conversion[type] , infolder, subfolder, "cpu time (user)/event (ms)", "\#events","",kTRUE);

	t[type]->Draw("1000*(pI.fCpuSys-pIOld.fCpuSys)>>h()",Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeSysDistribution",conversion[type] , infolder, subfolder, "cpu time (sys)/event (ms)", "\#events","",kTRUE);

	t[type]->Draw("1000*(pI.fCpuUser+pI.fCpuSys-pIOld.fCpuUser-pIOld.fCpuSys)>>h()",Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeUserSysDistribution",conversion[type] , infolder, subfolder, "cpu time (user+sys)/event (ms)", "\#events","",kTRUE);



	// Draw cpu, sys, cpu+sys time vs input size

	//t[type]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser):id0/1000>>h(1000,0,120,1000,0,15)",Form("sname==\"%s\"", sname[type]));


	int scale = type==0? 1000:1;
	t[type]->Draw( Form("1000*(pI.fCpuUser-pIOld.fCpuUser):id0/1024/%d>>h(100,0,%d,100,0,%d)", scale,inMax,timeMax) ,Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeUserVsInputSize",conversion[type] , infolder, subfolder, "input size (kB)", "cpu time (user)/event (ms)", "");

	//t[type]->Draw("1000*(pI.fCpuSys-pIOld.fCpuSys):id0/1000>>h(1000,0,120,1000,0,15)",Form("sname==\"%s\"", sname[type]));
	t[type]->Draw( Form("1000*(pI.fCpuSys-pIOld.fCpuSys):id0/1024/%d>>h(100,0,%d,100,0,%d)",scale,inMax,timeMax) ,Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeSysVsInputSize",conversion[type] , infolder, subfolder, "input size (kB)", "cpu time (sys)/event (ms)","");

	//t[type]->Draw("1000*(pI.fCpuUser+pI.fCpuSys-pIOld.fCpuUser-pIOld.fCpuSys):id0/1000>>h(1000,0,120,1000,0,15)",Form("sname==\"%s\"", sname[type]));
	t[type]->Draw( Form("1000*(pI.fCpuUser+pI.fCpuSys-pIOld.fCpuUser-pIOld.fCpuSys):id0/1024/%d>>h(100,0,%d,100,0,%d)",scale, inMax,timeMax) ,Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeUserSysVsInputSize",conversion[type] , infolder, subfolder, "input size (kB)", "cpu time (user+sys)/event (ms)","");


	// Draw cpu, sys, cpu+sys time vs output size

	//t[type]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser):id0/1000>>h(1000,0,120,1000,0,15)",Form("sname==\"%s\"", sname[type]));
	t[type]->Draw(Form("1000*(pI.fCpuUser-pIOld.fCpuUser):id1/1024>>h(1000,0,%d,1000,0,%d)",outMax,timeMax) ,Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeUserVsOutSize",conversion[type] , infolder, subfolder, "output size (kB)", "cpu time (user)/event (ms)");
	//t[type]->Draw("1000*(pI.fCpuSys-pIOld.fCpuSys):id0/1000>>h(1000,0,120,1000,0,15)",Form("sname==\"%s\"", sname[type]));
	t[type]->Draw(Form("1000*(pI.fCpuSys-pIOld.fCpuSys):id1/1024>>h(1000,0,%d,1000,0,%d)",outMax,timeMax) ,Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeSysVsOutSize",conversion[type] , infolder, subfolder, "output size (kB)", "cpu time (sys)/event (ms)");

	//t[type]->Draw("1000*(pI.fCpuUser+pI.fCpuSys-pIOld.fCpuUser-pIOld.fCpuSys):id0/1000>>h(1000,0,120,1000,0,15)",Form("sname==\"%s\"", sname[type]));
	t[type]->Draw(Form("1000*(pI.fCpuUser+pI.fCpuSys-pIOld.fCpuUser-pIOld.fCpuSys):id1/1024>>h(1000,0,%d,1000,0,%d)",outMax,timeMax) ,Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeUserSysVsOutSize",conversion[type] , infolder, subfolder, "output size (kB)", "cpu time (user+sys)/event (ms)");


	// Draw cpu, sys, cpu+sys time vs timestamp

	//t[type]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser):stampSec>>h(1000,1402524e3,1402534e3,1000,0,15)",Form("sname==\"%s\"", sname[type]));
	t[type]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser):stampSec>>h()",Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeUserVsTimestamp",conversion[type] , infolder, subfolder, "timestamp (s)", "cpu time (user)/event (ms)");

	//t[type]->Draw("1000*(pI.fCpuSys-pIOld.fCpuSys):stampSec>>h(1000,1402524e3,1402534e3,1000,0,15)",Form("sname==\"%s\"", sname[type]));
	t[type]->Draw("1000*(pI.fCpuSys-pIOld.fCpuSys):stampSec>>h()",Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeSysVsTimestamp",conversion[type] , infolder, subfolder, "timestamp (s)", "cpu time (sys)/event (ms)");

	//t[type]->Draw("1000*(pI.fCpuUser+pI.fCpuSys-pIOld.fCpuUser-pIOld.fCpuSys):stampSec>>h(1000,1402524e3,1402534e3,1000,0,15)",Form("sname==\"%s\"", sname[type]));
	t[type]->Draw("1000*(pI.fCpuUser+pI.fCpuSys-pIOld.fCpuUser-pIOld.fCpuSys):stampSec>>h()",Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeUserSysVsTimestamp",conversion[type] , infolder, subfolder, "timestamp (s)", "cpu time (user+sys)/event (ms)");


	// Draw cpu, sys, cpu+sys time vs time in event

	t[type]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser):T>>h()",Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeUserVsT",conversion[type] , infolder, subfolder, "time in file (s)", "cpu time (user)/event (ms)");

	t[type]->Draw("1000*(pI.fCpuSys-pIOld.fCpuSys):T>>h()",Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeSysVsT",conversion[type] , infolder, subfolder, "time in file (s)", "cpu time (sys)/event (ms)");

	t[type]->Draw("1000*(pI.fCpuUser+pI.fCpuSys-pIOld.fCpuUser-pIOld.fCpuSys):T>>h()",Form("sname==\"%s\"", sname[type]));
	saveHist("cpuTimeUserSysVsT",conversion[type] , infolder, subfolder, "time in file (s)", "cpu time (user+sys)/event (ms)");

	
	// Draw input size distribution

	t[type]->Draw("id0/1024>>h()",Form("sname==\"%s\"", sname[type]));
	saveHist("inputSize",conversion[type] , infolder, subfolder, "input size (kB)", "\#events","",kTRUE);


	// Draw output size distribution

	t[type]->Draw("id1/1024>>h()",Form("sname==\"%s\"", sname[type]));
	saveHist("outputSize",conversion[type] , infolder, subfolder, "output size (kB)", "\#events","",kTRUE);


	// Draw output vs input size

	t[type]->Draw(Form("id1/1024/%d:id0/1024/%d>>h()",scale,scale),Form("sname==\"%s\"", sname[type]));
	saveHist("OutVsIn",conversion[type] , infolder, subfolder, "input size (kB)", "output size (kB)","", kFALSE, kTRUE);

}
}


// combined histo



// draw time distribution for flat and esd in one histo
if(combined){
	cout<<"Creating combined plots"<<endl;

	infolder = "combined";
	subfolder = "";
	char* legends[2] = {"ESD","flatESD"};
	char* options[2] = {"",""};



	t[0]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser)>>h1(100,0,15)", Form("sname==\"%s\"", sname[0]) );
	t[1]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser)>>h2(100,0,15)", Form("sname==\"%s\"", sname[1]) );
//	t[0]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser)>>h1()", Form("sname==\"%s\"", sname[0]) );
//	t[1]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser)>>h2()", Form("sname==\"%s\"", sname[1]) );
	save2Hists("cpuTimeUserDistribution","flatVsNormal" , infolder, subfolder, "cpu time (user)/event (ms)", "\#events",legends, options,kTRUE,1) ;

	t[0]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser):id0/1024/1000>>h1(1000,0,2000,1000,0,16)", Form("sname==\"%s\"", sname[0]) );
	t[1]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser):id0/1024>>h2(1000,0,2000,1000,0,16)", Form("sname==\"%s\"", sname[1]) );
	//t[0]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser):id0/1024/1000>>h1()", Form("sname==\"%s\"", sname[0]) );
	//t[1]->Draw("1000*(pI.fCpuUser-pIOld.fCpuUser):id0/1024>>h2()", Form("sname==\"%s\"", sname[1]) );
	save2Hists("cpuTimeUserVsInputSize","flatVsNormal" , infolder, subfolder, "input size (kB)", "cpu time (user)/event (ms)",legends, options,kFALSE,2);




// loop over trees to get combined histos



	TH2D*flatVsNormalCpuTime = new TH2D("flatVsNormalCpuTime", "cpu time: flat vs normal", 400,0,8,150,0,3);
	TH2D*flatVsNormalInSize = new TH2D("flatVsNormalInSize", "input size: flat vs normal", 200,0,200,2000,0,2000);
	TH2D*flatVsNormalOutSize = new TH2D("flatVsNormalOutSize", "output size: flat vs normal", 80,0,80,800,0,800);
	Double_t cpuEsd , cpuOldEsd, cpuFlat, cpuOldFlat;
	Int_t inputSizeEsd, inputSizeFlat,outputSizeEsd, outputSizeFlat;
	char snameEsd[100]=".";
	char snameFlat[100]=".";

	t[0]->SetBranchAddress("pI.fCpuUser",&cpuEsd);
	t[0]->SetBranchAddress("pIOld.fCpuUser",&cpuOldEsd);
	t[0]->SetBranchAddress("sname",&snameEsd);
	t[0]->SetBranchAddress("id0",&inputSizeEsd);
	t[0]->SetBranchAddress("id1",&outputSizeEsd);

	t[1]->SetBranchAddress("pI.fCpuUser",&cpuFlat);
	t[1]->SetBranchAddress("pIOld.fCpuUser",&cpuOldFlat);
	t[1]->SetBranchAddress("sname",&snameFlat);
	t[1]->SetBranchAddress("id0",&inputSizeFlat);
	t[1]->SetBranchAddress("id1",&outputSizeFlat);


     t[0]->SetBranchStatus("*",0); //disable all branches
     t[0]->SetBranchStatus("sname",1);
     t[0]->SetBranchStatus("pI.fCpuUser",1);
     t[0]->SetBranchStatus("pIOld.fCpuUser",1);
     t[0]->SetBranchStatus("id0",1);
     t[0]->SetBranchStatus("id1",1);

     t[1]->SetBranchStatus("*",0); //disable all branches
     t[1]->SetBranchStatus("sname",1);
     t[1]->SetBranchStatus("pI.fCpuUser",1);
     t[1]->SetBranchStatus("pIOld.fCpuUser",1);
     t[1]->SetBranchStatus("id0",1);
     t[1]->SetBranchStatus("id1",1);

Int_t i2=0;

// loop over timestamps in normalESD tree
for (Int_t i1 = 0; i1 < t[0]->GetEntries()  ; i1++) {
	
	t[0]->GetEntry(i1);
	if( strcmp(snameEsd,"DoEvent.Stop") == 0){
		t[1]->GetEntry(i2++);
		while(strcmp(snameFlat,"DoEvent.Stop") != 0 && i2 < t[1]->GetEntries()){
			t[1]->GetEntry(i2++);
		}
		if(i2 < t[1]->GetEntries()){
//cout<<outputSizeEsd<<" "<<outputSizeFlat<<endl;
		  flatVsNormalCpuTime->Fill( 1000*(cpuEsd-cpuOldEsd), 1000*(cpuFlat-cpuOldFlat) );
		  flatVsNormalInSize->Fill( inputSizeEsd/1000/1024, inputSizeFlat/1024);
		  flatVsNormalOutSize->Fill( outputSizeEsd/1000/1024, outputSizeFlat/1024);
		}
	}
	
}


	saveHist(flatVsNormalCpuTime, infolder, subfolder, "cpu time ESD converter (ms)", "cpu time flatESD converter (ms)", "colz", kFALSE, kTRUE,kTRUE);
	saveHist(flatVsNormalInSize, infolder, subfolder, "input size time ESD converter (kB)", "input size flatESD converter (kB)", "colz");
	saveHist(flatVsNormalOutSize, infolder, subfolder, "output size ESD converter (kB)", "output size flatESD converter (kB)", "colz",kFALSE,kTRUE);
delete flatVsNormalCpuTime;
delete flatVsNormalInSize;
delete flatVsNormalOutSize;

}
return;
}

void saveHist(TH1*h, char*infolder, char*subfolder, char*x="",char*y="",char*options="", Bool_t log=kFALSE, Bool_t fit = kFALSE, Bool_t xy=kFALSE){
	TCanvas *c = new TCanvas();
	if(log){
		c->SetLogy();
	}
	h->GetXaxis()->SetTitle(x);
	h->GetYaxis()->SetTitle(y);
	h->SetStats(0);
		h->SetMarkerStyle(6);
	h->Draw(options);
	//c->SaveAs( Form("%s%s/%s.png", infolder, subfolder, h->GetName()) );
	//c->SaveAs( Form("%s%s/%s.root", infolder, subfolder, h->GetName()) );
	if(fit){

        TF1 *linear = new TF1("linear","pol2(0)", 0,2000);

        linear->SetParameters(0.5,0.2);
        linear->SetLineColor(kRed);
        linear->SetLineWidth(2);
		h->Fit(linear);

	}
	if(xy){

        TF1 *linear = new TF1("lin","x",0,8);
        linear->SetLineColor(kRed);
        linear->SetLineWidth(2);
		linear->Draw("same");

	}
	c->SaveAs( Form("%s%s/%s_fit.png", infolder, subfolder, h->GetName()) );
	c->SaveAs( Form("%s%s/%s_fit.root", infolder, subfolder, h->GetName()) );
}





void saveHist(char* name, char* type, char*infolder, char*subfolder, char*x="",char*y="", char*options="", Bool_t log=kFALSE, Bool_t fit = kFALSE){
	TCanvas *c = new TCanvas();
	if(log){
		c->SetLogy();
	}
	TH1* h = (TH1*)gDirectory->Get("h");
	cout<<name<<" mean: "<< h->GetMean()<<endl;
	h->SetTitle(Form("%s_%s", name, type) );
	h->GetXaxis()->SetTitle(x);
	h->GetYaxis()->SetTitle(y);
	h->SetStats(0);
	//h->SetMarkerStyle(6);
	h->Draw(options);
	if(fit){

        TF1 *linear = new TF1("linear","pol2(0)", 0,2000);

        linear->SetParameters(0.5,0.2);
        linear->SetLineColor(kRed);
        linear->SetLineWidth(2);
		h->Fit(linear);

	}


	c->SaveAs( Form("%s%s/%s_%s.png", infolder, subfolder, name, type) );
	c->SaveAs( Form("%s%s/%s_%s.root", infolder, subfolder, name, type) );
	gDirectory->DeleteAll();
}


void save2Hists(char* name, char* infolder, char*subfolder, char*x="",char*y="", char**legends =0x0, char**options=0x0,  Bool_t log=kFALSE, Int_t dim=1){
	TCanvas *c = new TCanvas();
	if(log){
		c->SetLogy();
	}
	TH1* h1 = (TH1*)gDirectory->Get("h1");
	TH1* h2 = (TH1*)gDirectory->Get("h2");

	if(dim==1){
		h1->SetFillColor(kRed);
		h2->SetFillColor(kBlue);
	}
	else if(dim==2){
		h1->SetMarkerColor(kRed);
		h2->SetMarkerColor(kBlue);
		h1->SetMarkerStyle(6);
		h2->SetMarkerStyle(6);
	}

	h1->SetTitle(Form("%s_%s", name, type) );
	h1->GetXaxis()->SetTitle(x);
	h1->GetYaxis()->SetTitle(y);
	h1->SetStats(0);
	h1->Draw();
	h2->Draw("same");
    TLegend* l=new TLegend(.58,0.68,.88,0.85);
	l->SetFillColor(0);
	l->SetBorderSize(0);

	if(dim==1){
		l->AddEntry(h1, legends[0], "f");
		l->AddEntry(h2, legends[1],"f");
	}
	else if(dim==2){
		l->AddEntry(h1, legends[0],"p");
		l->AddEntry(h2, legends[1],"p");
	}

	l->Draw("same");
	c->SaveAs( Form("%s%s/%s_%s.png", infolder, subfolder, name, type) );
	c->SaveAs( Form("%s%s/%s_%s.root", infolder, subfolder, name,type) );
	gDirectory->DeleteAll();
	delete l;
}
