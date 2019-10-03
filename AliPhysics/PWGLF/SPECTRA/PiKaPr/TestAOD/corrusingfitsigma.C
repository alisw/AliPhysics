//This macro calculates the correction factor for the difference of the n sigma TOF distributions in data and mc  


class AliSpectraBothHistoManager;
AliSpectraBothHistoManager* managerdata=0x0;
AliSpectraBothHistoManager* managermc=0x0;

Double_t Gaus_plus_tail(Double_t* x, Double_t* par )
{

	Double_t mean=par[0];
	Double_t rms=par[1];
	Double_t c=par[2];
	Double_t slope=par[3]/par[1];//df/dx continues
	Double_t cut=par[3]; 



	Double_t one_over_sqrt_2pi=1.0/(TMath::Sqrt(2.0*TMath::Pi()));
	Double_t returnvalue=0.0;

	Double_t n=0.5*(1.0+TMath::Erf(cut/TMath::Sqrt(2.0)))+TMath::Exp(-cut*cut*0.5)*one_over_sqrt_2pi/(TMath::Abs(rms)*slope);
	if (x[0]<mean+cut*rms)
		returnvalue=TMath::Exp(-1.0*(x[0]-mean)*(x[0]-mean)/(2.0*rms*rms))*one_over_sqrt_2pi/TMath::Abs(rms);
	else
		returnvalue=TMath::Exp(slope*(mean+cut*rms-x[0]))*TMath::Exp(-cut*cut*0.5)*one_over_sqrt_2pi/TMath::Abs(rms);
                                                                                                         
	return c*returnvalue/n;	
}

//Double_t background(Double_t x, Double_t scale, Double_t slope, Double_t startpoint)
Double_t background(Double_t* x,Double_t* par)
{
	//at start ponit it is eqaul to scale
	Double_t scale=par[0];
	Double_t slope=par[1];
	Double_t startpoint=par[2];
	return scale*TMath::Exp((startpoint-x[0])*slope);
}


Double_t Gaus_on_background(Double_t* x,Double_t* par)
{
	//par[0] peak position 
	//par[1] sigma gaus
	//par[2] Normalization of signal
	//par[3] cut
	//par[4] scale of bg
	//par[5] slope of bg
	//par[6] start point of bg const 
	//cout<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]<<" "<<par[4]<<" "<<par[5]<<" "<<par[6]<<par[7]<<" "<<x[0]<< endl;
	Double_t returnvalue=0.0;
	 returnvalue+=Gaus_plus_tail(x,par);
	//returnvalue+=par[3]*TMath::Exp((par[5]-x[0])*par[4]);
	returnvalue+=background(x,&(par[4]));
	return returnvalue;
}

Int_t   OpenFile(TString nameFile,TString outputname,Int_t mg1)
{
	

	TFile *file = TFile::Open(nameFile.Data());
	if(!file)
	{
		cout<<"no file"<<endl;
		return -1;
	}	
	TDirectoryFile *dir=(TDirectoryFile*)file->Get(outputname.Data());
	if(!dir)
	{
		cout<<"no dir "<<outputname.Data()<<endl;
		return -1;
	}
	if(!dir->Get("SpectraHistos"))
		return -1;
	if(mg1==1)
	{	
		managerdata= (AliSpectraBothHistoManager*) dir->Get("SpectraHistos");
		return managerdata->GetGenMulvsRawMulHistogram("hHistGenMulvsRawMul")->GetEntries();

	}
	else 
	{
		managermc= (AliSpectraBothHistoManager*) dir->Get("SpectraHistos");
		return managermc->GetGenMulvsRawMulHistogram("hHistGenMulvsRawMul")->GetEntries();
	}
	return -1;
}


void corrusingfitsigma(TString nameFile1,TString outputname1,TString nameFile2,TString outputname2, Float_t nsigmacut=3.0, Float_t minpt=0.6,Float_t maxpt=1.4, Bool_t saveoption=kFALSE)
{	
	ofstream oufiletxt("TOFcorrPID.txt");
	Float_t nsigmacut=TMath::Abs(nsigmacut);
	TString pidmethods[3]={"TPC","TOF","TPCTOF"};
	TString Particle[]={"Pion","Kaon","Proton"};
	gStyle->SetOptStat(0);	
	TH1::AddDirectory(kFALSE);
	gSystem->Load("libCore");
	gSystem->Load("libPhysics");
	gSystem->Load("libTree");
	gSystem->Load("libMatrix");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libOADB");
	gSystem->Load("libANALYSISalice");
	gSystem->Load("libTender");
	gSystem->Load("libCORRFW");
	gSystem->Load("libPWGTools");
	gSystem->Load("libPWGLFspectra");
	Int_t ne1= OpenFile(nameFile1,outputname1,1);
	Int_t ne2= OpenFile(nameFile2,outputname2,2);
	Int_t imethod=1;
		oufiletxt<<"file 1 "<<nameFile1.Data()<<endl;
	oufiletxt<<"outputname 1 "<<outputname1.Data()<<endl;
	oufiletxt<<"file 2 "<<nameFile2.Data()<<endl;
	oufiletxt<<"outputname 2 "<<outputname2.Data()<<endl;
	oufiletxt<<"nsigma "<<nsigmacut<<endl;
	oufiletxt<<"minpt "<<minpt<<endl;
	oufiletxt<<"maxpt "<<maxpt<<endl;
	
	TList* lout=new TList();

	for(Int_t ipart=0;ipart<3;ipart++)
	{
		if(!managerdata->GetNSigHistogram(Form("hHistNSig%sPt%s",Particle[ipart].Data(),pidmethods[imethod].Data())))
			continue;		
		TH2F *nsig_data = (TH2F*)((TH2F*)managerdata->GetNSigHistogram(Form("hHistNSig%sPt%s",Particle[ipart].Data(),pidmethods[imethod].Data())))->Clone();
		if(!nsig_data)
			continue;

		if(!managermc->GetNSigHistogram(Form("hHistNSig%sPt%s",Particle[ipart].Data(),pidmethods[imethod].Data())))
			continue;
		TH2F *nsig_mc = (TH2F*)((TH2F*)managermc->GetNSigHistogram(Form("hHistNSig%sPt%s",Particle[ipart].Data(),pidmethods[imethod].Data())))->Clone();
		if(!nsig_mc)
			continue;
		Int_t firstbin=nsig_data->GetXaxis()->FindBin(minpt);
		Int_t lastbin=nsig_data->GetXaxis()->FindBin(maxpt);
		Float_t ptlow=nsig_data->GetXaxis()->GetBinLowEdge(firstbin);
		Float_t ptup=nsig_data->GetXaxis()->GetBinUpEdge(lastbin);

		TH1D* hist1data= new TH1D(Form("datahist3sigmafun4sigma%s",Particle[ipart].Data()),"",lastbin-firstbin+1,ptlow,ptup);
		TH1D* hist2data= new TH1D(Form("datahist3sigmahist4sigma%s",Particle[ipart].Data()),"",lastbin-firstbin+1,ptlow,ptup);
		TH1D* sigmadata= new TH1D(Form("sigmadata%s",Particle[ipart].Data()),"",lastbin-firstbin+1,ptlow,ptup);
	
		hist1data->SetLineColor(kBlack);
		hist2data->SetLineColor(kBlack);
		sigmadata->SetMarkerColor(kBlack);
		sigmadata->SetMarkerStyle(22);
		hist1data->SetLineStyle(1);
		hist2data->SetLineStyle(2);
		hist1data->SetLineWidth(2);
		hist2data->SetLineWidth(2);

	
		TH1D* hist1mc= new TH1D(Form("mchist3sigmafun4sigma%s",Particle[ipart].Data()),"",lastbin-firstbin+1,ptlow,ptup);
		TH1D* hist2mc= new TH1D(Form("mchist3sigmahist4sigma%s",Particle[ipart].Data()),"",lastbin-firstbin+1,ptlow,ptup);
		TH1D* sigmamc= new TH1D(Form("sigmamc%s",Particle[ipart].Data()),"",lastbin-firstbin+1,ptlow,ptup);

		hist1mc->SetLineColor(kRed);
		hist2mc->SetLineColor(kRed);
		sigmamc->SetMarkerColor(kRed);
		sigmamc->SetMarkerStyle(23);
		sigmamc->SetLineColor(kRed);
		hist1mc->SetLineStyle(1);
		hist2mc->SetLineStyle(2);
		hist1mc->SetLineWidth(2);
		hist2mc->SetLineWidth(2);



		for(int ibin=firstbin;ibin<lastbin;ibin++)
		{
				
			TH1F *nsig_data_Proj1=(TH1F*)nsig_data->ProjectionY(Form("%s%sdata[%.2f,%.2f]",Particle[ipart].Data(),pidmethods[imethod].Data(),nsig_data->GetXaxis()->GetBinLowEdge(ibin),nsig_data->GetXaxis()->GetBinUpEdge(ibin)),ibin,ibin);
			nsig_data_Proj1->Sumw2();
			nsig_data_Proj1->GetXaxis()->SetRangeUser(-6,6);
			nsig_data_Proj1->SetLineColor(kBlack);		
			TH1F *nsig_mc_Proj1=(TH1F*)nsig_mc->ProjectionY(Form("%s%smc[%.2f,%.2f]",Particle[ipart].Data(),pidmethods[imethod].Data(),nsig_mc->GetXaxis()->GetBinLowEdge(ibin),nsig_mc->GetXaxis()->GetBinUpEdge(ibin)),ibin,ibin);
                	nsig_mc_Proj1->Sumw2();
                	nsig_mc_Proj1->GetXaxis()->SetRangeUser(-6,6);
			nsig_mc_Proj1->SetLineColor(kRed);
			cout<<"data"<<endl;
			TF1* gfundata=new TF1("fun_for_fitdata",Gaus_on_background,-5,5,7);
			gfundata->SetParameter(0,0.0);
			gfundata->SetParameter(1,1.0);
			gfundata->SetParameter(2,nsig_data_Proj1->GetBinContent(nsig_data_Proj1->GetXaxis()->FindBin(0.0)));
			gfundata->FixParameter(3,1.0);
			gfundata->FixParameter(4,0.0);
			gfundata->FixParameter(5,1.0);
			gfundata->FixParameter(6,-5);
			gfundata->SetLineColor(kBlack);
			nsig_data_Proj1->Fit(gfundata,"SR+","E");
			Float_t sigmadatav=gfundata->GetParameter(1);
			Float_t meandatav=gfundata->GetParameter(0);
			sigmadata->SetBinContent(ibin-firstbin+1,gfundata->GetParameter(1));	
			sigmadata->SetBinError(ibin-firstbin+1,gfundata->GetParError(1));

			cout<<"MC"<<endl;
			TF1* gfunmc=new TF1("fun_for_fitmc",Gaus_on_background,-5,5,7);
			gfunmc->SetParameter(0,0.0);
			gfunmc->SetParameter(1,1.0);
			gfunmc->SetParameter(2,nsig_mc_Proj1->GetBinContent(nsig_mc_Proj1->GetXaxis()->FindBin(0.0)));
			gfunmc->FixParameter(3,1.0);
			gfunmc->FixParameter(4,0.0);
			gfunmc->FixParameter(5,1.0);
			gfunmc->FixParameter(6,-5);
			gfunmc->SetLineColor(kRed);
			nsig_mc_Proj1->Fit(gfunmc,"SR+","E");
			sigmamc->SetBinContent(ibin-firstbin+1,gfunmc->GetParameter(1));	
			sigmamc->SetBinError(ibin-firstbin+1,gfunmc->GetParError(1));
			Float_t sigmamcv=gfunmc->GetParameter(1);
			Float_t meanmcv=gfunmc->GetParameter(0);


			hist1data->SetBinContent(ibin-firstbin+1,nsig_data_Proj1->Integral(nsig_data_Proj1->GetXaxis()->FindBin(-1.0*nsigmacut),nsig_data_Proj1->GetXaxis()->FindBin(nsigmacut-0.0001)));			
			Float_t width=nsig_data_Proj1->GetXaxis()->GetBinWidth(-1.0*nsigmacut);
			Int_t lowbind=nsig_data_Proj1->GetXaxis()->FindBin(-1.0*nsigmacut*sigmadatav+meandatav);
			Int_t upbind=nsig_data_Proj1->GetXaxis()->FindBin(nsigmacut*sigmadatav+meandatav);
			Float_t lowbinpartd=TMath::Abs(nsig_data_Proj1->GetXaxis()->GetBinLowEdge(lowbind)+1.0*nsigmacut*sigmadatav-meandatav);
			Float_t upbinpartd=TMath::Abs(nsigmacut*sigmadatav+meandatav-nsig_data_Proj1->GetXaxis()->GetBinUpEdge(upbind));
			
//nsig_data_Proj1->Print("all");
//	nsig_mc_Proj1->Print("all");

			Float_t corrd=(lowbinpartd*nsig_data_Proj1->GetBinContent(lowbind)+upbinpartd*nsig_data_Proj1->GetBinContent(upbind))/width;

			hist2data->SetBinContent(ibin-firstbin+1,nsig_data_Proj1->Integral(nsig_data_Proj1->GetXaxis()->FindBin(-1.0*nsigmacut*sigmadatav+meandatav),nsig_data_Proj1->GetXaxis()->FindBin(nsigmacut*sigmadatav+meandatav))-corrd);
			cout<<corrd/nsig_data_Proj1->Integral(nsig_data_Proj1->GetXaxis()->FindBin(-1.0*nsigmacut*sigmadatav+meandatav),nsig_data_Proj1->GetXaxis()->FindBin(nsigmacut*sigmadatav+meandatav))<<endl;

			hist1mc->SetBinContent(ibin-firstbin+1,nsig_mc_Proj1->Integral(nsig_mc_Proj1->GetXaxis()->FindBin(-1.0*nsigmacut),nsig_mc_Proj1->GetXaxis()->FindBin(nsigmacut-0.0001)));
			width=nsig_mc_Proj1->GetXaxis()->GetBinWidth(-1.0*nsigmacut);
			Int_t lowbinmc=nsig_mc_Proj1->GetXaxis()->FindBin(-1.0*nsigmacut*sigmamcv+meanmcv);
			Int_t upbinmc=nsig_mc_Proj1->GetXaxis()->FindBin(nsigmacut*sigmamcv+meanmcv);
			Float_t lowbinpartmc=TMath::Abs(nsig_mc_Proj1->GetXaxis()->GetBinLowEdge(lowbinmc)+1.0*nsigmacut*sigmamcv-meanmcv);
			Float_t upbinpartmc=TMath::Abs(nsigmacut*sigmamcv+meanmcv-nsig_mc_Proj1->GetXaxis()->GetBinUpEdge(upbinmc));
			cout<<lowbinpartmc<<" "<<upbinpartmc<<endl;
			Float_t corrmc=(lowbinpartmc*nsig_mc_Proj1->GetBinContent(lowbinmc)+upbinpartmc*nsig_mc_Proj1->GetBinContent(upbinmc))/width;
			cout<<corrmc/nsig_mc_Proj1->Integral(nsig_mc_Proj1->GetXaxis()->FindBin(-1.0*nsigmacut*sigmamcv+meanmcv),nsig_mc_Proj1->GetXaxis()->FindBin(nsigmacut*sigmamcv+meanmcv))<<" aaa "<<endl;

			hist2mc->SetBinContent(ibin-firstbin+1,nsig_mc_Proj1->Integral(nsig_mc_Proj1->GetXaxis()->FindBin(-1.0*nsigmacut*sigmamcv+meanmcv),nsig_mc_Proj1->GetXaxis()->FindBin(nsigmacut*sigmamcv+meanmcv))-corrmc);
			TCanvas* ctmp=new TCanvas(Form("c%d%s",ibin,Particle[ipart].Data()),Form("c%d%s",ibin,Particle[ipart].Data()),800,800);
			ctmp->cd()->SetLogy();

			nsig_data_Proj1->GetXaxis()->SetTitle("#sigma_{TOF}");
			nsig_data_Proj1->GetYaxis()->SetTitle("N_{entries}");
			nsig_data_Proj1->SetTitle(Form("PionTOF[%.2f,%.2f]",nsig_data->GetXaxis()->GetBinLowEdge(ibin),nsig_data->GetXaxis()->GetBinUpEdge(ibin)));
			nsig_data_Proj1->Draw("E1");
			nsig_mc_Proj1->Draw("E1same");
			if(saveoption)
				ctmp->SaveAs(Form("c%d%s.png",ibin,Particle[ipart].Data()));
			delete ctmp;
			delete gfundata;
			delete gfunmc;
			delete nsig_data_Proj1;
			delete nsig_mc_Proj1;
			
		}
	
			hist1data->Sumw2();
			hist1data->Print("all");
			hist1data->Sumw2();
			hist1mc->Sumw2();
			hist2mc->Sumw2();

		hist2data->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		hist2data->GetYaxis()->SetTitle("N(-3fit,3fit)/N(-3,3)");
		/*hist1data->Print("all");
		hist2data->Print("all");
		hist1mc->Print("all");
		hist2mc->Print("all");
*/
		hist2mc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		hist2mc->GetYaxis()->SetTitle("N(-3fit,3fit)/N(-3,3)");
		hist2data->Divide(hist2data,hist1data,1,1,"B");
		hist2mc->Divide(hist2mc,hist1mc,1,1,"B");
		lout->Add(hist2data);	
		lout->Add(hist2mc);	

		sigmadata->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		sigmadata->GetYaxis()->SetTitle("#sigma_{TOF} fitted");
		sigmamc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		sigmamc->GetYaxis()->SetTitle("#sigma_{TOF} fitted");
		

		
		lout->Add(sigmadata);
		lout->Add(sigmamc);

		TH1D* corr=(TH1D*)hist2data->Clone(Form("corr%s",Particle[ipart].Data()));	
		corr->Divide(hist2mc);
		corr->SetLineColor(kBlue);
		corr->GetYaxis()->SetTitle("correction");
		corr->Fit("pol0","R0","",minpt,maxpt);
		
		oufiletxt<<Particle[ipart].Data()<<" "<<((TF1*)corr->GetListOfFunctions()->At(0))->GetParameter(0)<<endl;
		lout->Add(corr);
		TCanvas* corrcanvas=new TCanvas(Form("corrcanvas%s",Particle[ipart].Data()), Form("corrcanvas%s",Particle[ipart].Data()),600,600);
		corrcanvas->SetMargin(0.12,0.05,0.15,0.1);
		corr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		corr->GetYaxis()->SetTitle("TOF PID mismatch corr.");
		corr->GetYaxis()->SetTitleOffset(1.5*corr->GetYaxis()->GetTitleOffset());
		corr->GetYaxis()->SetRangeUser(0.98,1.04);
	
		corrcanvas->cd();
		corr->Draw("E1");
		corrcanvas->SaveAs(Form("Corr%s.eps",Particle[ipart].Data()));
		corrcanvas->SaveAs(Form("Corr%s.png",Particle[ipart].Data()));

	}
	TFile* fileout=TFile::Open("TOFcorrPID.root","Recreate");
	
	lout->Write("output",TObject::kSingleKey);
	fileout->Close();
}
