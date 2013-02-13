Float_t QAPlotsBoth( AliSpectraBothHistoManager* hman_data, AliSpectraBothHistoManager* hman_mc,
	      AliSpectraBothEventCuts* ecuts_data, AliSpectraBothEventCuts* ecuts_mc,
	      AliSpectraBothTrackCuts* tcuts_data, AliSpectraBothTrackCuts* tcuts_mc,
	      TList * flistqa,TList * flistcanvas)
{
TString pidmethods[3]={"TPC","TOF","TPCTOF"};	
	Double_t neventsdata =  ecutsdata->NumberOfPhysSelEvents();
	Double_t neventsmc =  ecutsmc->NumberOfPhysSelEvents();
	
	
	
	for(Int_t ipart=0;ipart<3;ipart++)
	{
			
			for(Int_t imethod=0;imethod<3;imethod++)
			{
				 TH2F *nsig_data = (TH2F*)((TH2F*)hman_data->GetNSigHistogram(Form("hHistNSig%sPt%s",Particle[ipart].Data(),pidmethods[imethod].Data())))->Clone();
				// nsig_data->RebinX(20);			 
				// nsig_data->RebinY(4);
				// nsig_data->Sumw2();

				 TH2F *nsig_mc = (TH2F*)((TH2F*)hman_mc->GetNSigHistogram(Form("hHistNSig%sPt%s",Particle[ipart].Data(),pidmethods[imethod].Data())))->Clone();
				 //nsig_mc->RebinX(20);			 
				// nsig_mc->RebinY(4);
				// nsig_mc->Sumw2();
				 
				Int_t ibin=1;
				Float_t binsize=nsig_mc->GetXaxis()->GetBinWidth(1);

				TH1F* maxposdata=(TH1F*)nsig_data->ProjectionX(Form("%s%sdatamaxpos",Particle[ipart].Data(),pidmethods[imethod].Data()),-1,-1));
				maxposdata->Reset();
				maxposdata->SetTitle(";p_{T} (GeV/c);max in (-2,2)");
				TH1F* maxposmc=(TH1F*)nsig_data->ProjectionX(Form("%s%smcmaxpos",Particle[ipart].Data(),pidmethods[imethod].Data()),-1,-1));
				maxposmc->Reset();
				maxposmc->SetTitle(";p_{T} (GeV/c);max in (-2,2)");

	
				 while (ibin*binsize<3.0)
				 {
					// TCanvas* c=new TCanvas(Form("canvas%s%s%d",Particle[ipart].Data(),pidmethods[imethod].Data(),ibin),Form("canvas%s%s%d",Particle[ipart].Data(),pidmethods[imethod].Data(),ibin),700,500);
					
					 TH1F *nsig_data_Proj1=(TH1F*)nsig_data->ProjectionY(Form("%s%sdata[%.2f,%.2f]",Particle[ipart].Data(),pidmethods[imethod].Data(),nsig_data->GetXaxis()->GetBinLowEdge(ibin),nsig_data->GetXaxis()->GetBinUpEdge(ibin)),ibin,ibin));
					 TH1F *nsig_mc_Proj1=(TH1F*)nsig_mc->ProjectionY(Form("%s%smc[%.2f,%.2f]",Particle[ipart].Data(),pidmethods[imethod].Data(),nsig_mc->GetXaxis()->GetBinLowEdge(ibin),nsig_mc->GetXaxis()->GetBinUpEdge(ibin)),ibin,ibin));
					 nsig_data_Proj1->GetXaxis()->SetRangeUser(-3,3);
					 nsig_data_Proj1->SetLineColor(kRed);
					 if(nsig_data_Proj1->Integral()<1&&nsig_mc_Proj1->Integral()<1)
					 {	
						ibin++;	
						continue;
					 } 
									//	nsig_data_Proj1->Sumw2();
				//	nsig_mc_Proj1->Sumw2();
					  //nsig_data_Proj1->GetXaxis()->SetRangeUser(-5,5);
					 
					 //c->cd()->SetLogy();
					 //nsig_data_Proj1->Draw();
					 //nsig_mc_Proj1->Draw("same");
					nsig_data_Proj1->GetXaxis()->SetRange(nsig_data_Proj1->GetXaxis()->FindBin(-2.0),nsig_data_Proj1->GetXaxis()->FindBin(2.0));
					if(nsig_data_Proj1->GetMaximumBin()<=nsig_data_Proj1->GetXaxis()->FindBin(2.0)&&nsig_data_Proj1->GetMaximumBin()>=nsig_data_Proj1->GetXaxis()->FindBin(-2.0))
					{
						maxposdata->SetBinContent(ibin,nsig_data_Proj1->GetXaxis()->GetBinCenter(nsig_data_Proj1->GetMaximumBin()));	
						maxposdata->SetBinError(ibin,nsig_data_Proj1->GetXaxis()->GetBinWidth(nsig_data_Proj1->GetMaximumBin())/2.0);	
	 				}
				cout<<Form("%s%sdatamaxpos",Particle[ipart].Data(),pidmethods[imethod].Data())<<" "<<nsig_data_Proj1->GetMaximumBin()<<" "<<nsig_data_Proj1->GetXaxis()->FindBin(2.0)<<" "<<nsig_data_Proj1->GetXaxis()->FindBin(-2.0)<<" "<<nsig_data_Proj1->GetXaxis()->GetBinCenter(nsig_data_Proj1->GetMaximumBin())<<" "<<ibin<<endl;

					nsig_data_Proj1->GetXaxis()->SetRange(0,nsig_data_Proj1->GetXaxis()->GetNbins());
									
					nsig_mc_Proj1->GetXaxis()->SetRange(nsig_mc_Proj1->GetXaxis()->FindBin(-2.0),nsig_mc_Proj1->GetXaxis()->FindBin(2.0));
					if(nsig_mc_Proj1->GetMaximumBin()<=nsig_mc_Proj1->GetXaxis()->FindBin(2.0)&&nsig_mc_Proj1->GetMaximumBin()>=nsig_mc_Proj1->GetXaxis()->FindBin(-2.0))
					{
						maxposmc->SetBinContent(ibin,nsig_mc_Proj1->GetXaxis()->GetBinCenter(nsig_mc_Proj1->GetMaximumBin()));	
						maxposmc->SetBinError(ibin,nsig_mc_Proj1->GetXaxis()->GetBinWidth(nsig_mc_Proj1->GetMaximumBin())/2.0);	
	 				}
						cout<<Form("%s%smcmaxpos",Particle[ipart].Data(),pidmethods[imethod].Data())<<" "<<nsig_mc_Proj1->GetMaximumBin()<<" "<<nsig_mc_Proj1->GetXaxis()->FindBin(2.0)<<" "<<nsig_mc_Proj1->GetXaxis()->FindBin(-2.0)<<" "<<ibin<<" "<<nsig_mc_Proj1->GetXaxis()->GetBinCenter(nsig_mc_Proj1->GetMaximumBin())<<endl;

					nsig_mc_Proj1->GetXaxis()->SetRange(0,nsig_mc_Proj1->GetXaxis()->GetNbins());


					nsig_data_Proj1->Scale(1.0/neventsdata);
					 nsig_mc_Proj1->Scale(1.0/neventsmc);

					
					flistqa->Add(nsig_data_Proj1);
					flistqa->Add(nsig_mc_Proj1);
					ibin++;
				 }
				flistqa->Add(maxposmc);
				flistqa->Add(maxposdata);
			}
	}
	TH1F* fHistoVtxAftSeldata=(TH1F*)ecuts_data->GetHistoVtxAftSel();
	TH1F* fHistoVtxAftSelmc=(TH1F*)ecuts_mc->GetHistoVtxAftSel();
	flistcanvas->Add(plot_on_canvas("vertex",fHistoVtxAftSeldata,fHistoVtxAftSelmc));
	TF1* fdata=new TF1("dataveretxfit","gausn");
	TF1* fmc=new TF1("mcveretxfit","gausn");
	fHistoVtxAftSeldata->Fit("dataveretxfit","0");
	fHistoVtxAftSelmc->Fit("mcveretxfit","0");
	Float_t datavertexratio=fHistoVtxAftSeldata->Integral(-1,-1,"width")/fdata->GetParameter(0);
	Float_t mcvertexratio=fHistoVtxAftSelmc->Integral(-1,-1,"width")/fmc->GetParameter(0);
	

	 TH1F* fHistoEtaAftSeldata=(TH1F*)ecuts_data->GetHistoEtaAftSel();
	 TH1F* fHistoEtaAftSelmc=(TH1F*)ecuts_mc->GetHistoEtaAftSel();
	flistcanvas->Add(plot_on_canvas("ETA",fHistoEtaAftSeldata,fHistoEtaAftSelmc));


	 TH1F* fITSclustershistdata=(TH1F*)tcuts_data->GetHistoNclustersITS();
 	  TH1F* fITSclustershistmc=(TH1F*)tcuts_mc->GetHistoNclustersITS();

	flistcanvas->Add(plot_on_canvas("NITS",fITSclustershistdata,fITSclustershistmc));
	cout<<" data "<<datavertexratio<<" mc "<<mcvertexratio<<endl;
	
	TH2F* hmul=(TH2F*)hman_mc->GetGenMulvsRawMulHistogram("hHistGenMulvsRawMul");	
	hmul->Sumw2();	
	TCanvas* cbc=new TCanvas("broken chunks","broken chunks",1200,600);
	cbc->Divide(2,1);
	cbc->cd(1);
	hmul->Draw();
	cbc->cd(2);
	TH1F* nonzero=(TH1F*)hmul->ProjectionX("nonzerotracks",2,-1);
	nonzero->SetMarkerColor(kRed);
	nonzero->SetMarkerStyle(21);
	TH1F* binzero=(TH1F*)hmul->ProjectionX("binzerotracks",1,1);
	binzero->SetMarkerColor(kBlack);
	binzero->SetMarkerStyle(22);
	binzero->Sumw2();
	nonzero->Sumw2();
	binzero->Divide(nonzero);
	TF1* badchunk=new TF1("badchunkfit","pol0",10,40);
	binzero->Fit("badchunkfit","R");
	Float_t badchunksfraction=badchunk->GetParameter(0);
	binzero->Draw("E1");
	flistcanvas->Add(cbc);
	
	return (1.0-badchunksfraction)*mcvertexratio/datavertexratio;

}
TCanvas* plot_on_canvas(TString name, TH1* h1,TH1* h2)
{
	TCanvas* cvrt=new TCanvas(name.Data(),name.Data(),600,600);
	cvrt->cd();
	h1->SetLineColor(kRed);
	h2->SetLineColor(kBlue);
	h1->Sumw2();
	h2->Sumw2();
	TLegend *lvtr=new TLegend(0.2,0.2,0.5,0.3,"","NDC");
	lvtr->SetLineColor(kWhite);
	lvtr->AddEntry(h1,"data","l");
	lvtr->AddEntry(h2,"MC","l");
	h1->Scale(1.0/h1->GetBinContent(h1->GetXaxis()->FindBin(0.0)));
	h2->Scale(1.0/h2->GetBinContent(h2->GetXaxis()->FindBin(0.0)));
	h1->DrawCopy("L");
	h2->DrawCopy("Lsame");
	lvtr->Draw();
	return cvrt;
}

