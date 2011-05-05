class AliOADBContainer;
class AliOADBPWG2Spectra;
Int_t colors[]={ kBlack,kRed,kBlue+2,kGreen+2, kSpring, kTeal, kAzure};
Int_t markers[]={21,22,23,24,25,26,27,28};
Float_t min[]={0.1,0.2,0.3,-1.0,-1.0,-1.0,0.5,0.5,0.8,-1.0,-1.0,-1.0};
Float_t max[]={0.5,0.5,0.5,-1.0,-1.0,-1.0,2.0,2.0,2.5,1.5,1.5,1.7};
Float_t min2[]={0.1,0.2,0.3,-1.0,-1.0,-1.0,0.5,0.5,0.8,-1.0,-1.0,-1.0};
Float_t max2[]={0.5,0.5,0.5,-1.0,-1.0,-1.0,2.0,2.0,2.5,0.1,0.1,0.1};
TString mcmodels[]={"Phojet","Pythia109","Pythia320"};
const char * fgkDetectorNames[] = {"ITS", "ITSTPC", "TPC", "TOF", "TOFTPC", "Dummy", "Dummy"};
const char * fgkPidTypeNames[]  = {"GaussFit", "NSigma", "Bayes", "Kinks"};
const char * partext[]={"#pi^{+}","K^{+}","p","#pi^{-}","K^{-}","#bar{p}"};
void ALICEWorkInProgress(TCanvas *c,TString today="11/05/2010", TString label = "ALICE performance");
TH1D* DrawHisto(AliOADBContainer* fOADBContainer, TCanvas* c1,Int_t ana,Int_t method,Int_t par,Int_t charge,Int_t color, Int_t marker,TLegend* leg, TString estimator ,Int_t bin);
Double_t myLevyPt(Double_t *pt, Double_t *par);
TH1D* CombineSpectra(TList* l,Int_t startpoint, Float_t min* , Float_t* max);
void DrawandFit(TH1D* hist,TCanvas* c1,Int_t par,TCanvas* c1ratio,Float_t* yileds,Float_t* yiledserr);
void GetMCratios(TString filename, TList* ktopi,TList* ptopi);
void plotRatios(TCanvas* c1, TH1D* data , TList* mc,TLegend* leg ,Int_t opt=0);
TH1D* Plotratiodatafit(TH1D* hist1, TF1* fun);
TCanvas* Ktopi(Float_t value,Float_t err);
void DrawMulRatios(TH1D** bin1,TH1D** bin2,TCanvas* c,TLegend* leg);


void PlotCombine7TeV()
{
	gSystem->Load("libCore.so");  
	gSystem->Load("libGeom.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libVMC");
	gSystem->Load("libTree");
	gSystem->Load("libProof");
	gSystem->Load("libMatrix");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libOADB");
	gSystem->Load("libANALYSISalice");
	gSystem->Load("libTENDER");
	gSystem->Load("libCORRFW");
	gSystem->Load("libPWG0base");
	gSystem->Load("libMinuit");
	gSystem->Load("libPWG2spectra");
	
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1,0);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	Int_t det[]={AliOADBPWG2Spectra::kITSsa,AliOADBPWG2Spectra::kITSTPC,AliOADBPWG2Spectra::kTOF,AliOADBPWG2Spectra::kTOFTPC} ;
	Int_t method[]={AliOADBPWG2Spectra::kGaussFit,AliOADBPWG2Spectra::kGaussFit,AliOADBPWG2Spectra::kGaussFit,AliOADBPWG2Spectra::kNSigma };
	
	
	 TString fileName = AliOADBPWG2Spectra::GetOADBPWG2SpectraFileName();
	TFile * f = new TFile (fileName);
	AliOADBContainer* fOADBContainer = (AliOADBContainer*) f->Get("Corrected");
	f->Close();
	if(!fOADBContainer)
		return;
	AliOADBPWG2Spectra* fOADBSpectra = (AliOADBPWG2Spectra*) fOADBContainer->GetObject(116562); 
	if(!fOADBSpectra)
		return;
	fOADBSpectra->Print();	
	TCanvas* c1= new TCanvas("pos","pos",1200,800);
	c1->cd()->SetLogy();
	TLegend* Leg1 = new TLegend(0.1,0.1,0.45,0.45,"","NDC");
	Leg1->SetNColumns(3);
	Leg1->SetTextSize(0.027);
	Leg1->SetFillStyle(kFALSE);
	Leg1->SetLineColor(kWhite);
	Leg1->SetBorderSize(0);
	
	
	TCanvas* c2= new TCanvas("neg","neg",1200,800);
	c2->cd()->SetLogy();
	TLegend* Leg2 = new TLegend(0.1,0.1,0.45,0.45,"","NDC");
	Leg2->SetNColumns(3);
	Leg2->SetTextSize(0.027);
	Leg2->SetFillStyle(kFALSE);
	Leg2->SetLineColor(kWhite);
	Leg2->SetBorderSize(0);
	
	
	TDatime date;
	TList* posparthist=new TList();
	TList* negparthist=new TList();
	
	TList* posparthistbin1=new TList();
	TList* negparthistbin1=new TList();
	
	TList* posparthistbin4=new TList();
	TList* negparthistbin4=new TList();
	
	for(int i=0;i<4;i++)// ITS , ITSTPC, TOF TPCTOF 
	{
		for(int j=0;j<3;j++) //pion kaon proton
		{	
			posparthist->Add(DrawHisto(fOADBSpectra,c1,det[i],method[i],j,0,colors[j],markers[i],Leg1));
			negparthist->Add(DrawHisto(fOADBSpectra,c2,det[i],method[i],j,1,colors[j],markers[i],Leg2));
			posparthistbin1->Add(DrawHisto(fOADBSpectra,0x0,det[i],method[i],j,0,colors[j],markers[i],0x0,"SPD2",1));
			negparthistbin1->Add(DrawHisto(fOADBSpectra,0x0,det[i],method[i],j,1,colors[j],markers[i],0x0,"SPD2",1));
			posparthistbin4->Add(DrawHisto(fOADBSpectra,0x0,det[i],method[i],j,0,colors[j],markers[i],0x0,"SPD2",4));
			negparthistbin4->Add(DrawHisto(fOADBSpectra,0x0,det[i],method[i],j,1,colors[j],markers[i],0x0,"SPD2",4));
		}	
	}	
	//negparthistbin4->ls();
	ALICEWorkInProgress(c1,Form("%d/%d/%d",date.GetDay(),date.GetMonth(),date.GetYear()),"ALICE preliminary");
	Leg1->Draw();
	ALICEWorkInProgress(c2,Form("%d/%d/%d",date.GetDay(),date.GetMonth(),date.GetYear()),"ALICE preliminary");
	Leg2->Draw();
	
	c1->SaveAs("Pos.eps");
	c2->SaveAs("Neg.eps");
	
	
	TH1D* combinepos[3];
	TH1D* combineneg[3];
	TH1D* combineposbin1[3];
	TH1D* combinenegbin1[3];
	TH1D* combineposbin4[3];
	TH1D* combinenegbin4[3];
	TCanvas* c1fit= new TCanvas("posfit","posfit",1200,800);
	c1fit->cd()->SetLogy();
	TCanvas* c2fit= new TCanvas("negfit","negfit",1200,800);
	c2fit->cd()->SetLogy();
	TCanvas* c1fitratio= new TCanvas("posfitratio","posfitrato",1200,800);
	TCanvas* c2fitratio= new TCanvas("negfitratio","negfitratio",1200,800);
	Float_t yileds[3]={0.0,0.0,0.0};
	Float_t yiledserr[3]={0.0,0.0,0.0};
	for(int j=0;j<3;j++) //pion kaon proton
	{
			combinepos[j]=CombineSpectra(posparthist,j, min,max);
			combineposbin1[j]=CombineSpectra(posparthistbin1,j,min2,max2);
			combineposbin4[j]=CombineSpectra(posparthistbin4,j,min2,max2);
			DrawandFit(combinepos[j],c1fit,j,c1fitratio,yileds,yiledserr);
			combinepos[j]->Sumw2();
			combineposbin1[j]->Sumw2();
			combineposbin4[j]->Sumw2();
			combineposbin1[j]->Divide(combinepos[j]);
			combineposbin4[j]->Divide(combinepos[j]);
			combineneg[j]=CombineSpectra(negparthist,j, min,max);
			combinenegbin1[j]=CombineSpectra(negparthistbin1,j,min2,max2);
			combinenegbin4[j]=CombineSpectra(negparthistbin4,j,min2,max2);
			DrawandFit(combineneg[j],c2fit,j,c2fitratio,yileds,yiledserr);
			combineneg[j]->Sumw2();
			combinenegbin1[j]->Sumw2();
			combinenegbin4[j]->Sumw2();
			combinenegbin1[j]->Divide(combinepos[j]);
			combinenegbin4[j]->Divide(combinepos[j]);
			combinepos[j]->Add(combineneg[j]);
				
	}
	ALICEWorkInProgress(c1fit,Form("%d/%d/%d",date.GetDay(),date.GetMonth(),date.GetYear()),"ALICE preliminary");
	ALICEWorkInProgress(c2fit,Form("%d/%d/%d",date.GetDay(),date.GetMonth(),date.GetYear()),"ALICE preliminary");
	c1fit->SaveAs("posfit.eps");
	c2fit->SaveAs("negfit.eps");
	c1fitratio->SaveAs("posfitratio.eps");
	c2fitratio->SaveAs("negfitratio.eps");
	combinepos[1]->Divide(combinepos[0]);
	combinepos[2]->Divide(combinepos[0]);
	
	TString mcfiles[3]={"/home/marek/Analysis/Spectra/7TeVLHC10b/MCmodels/AnalyseFastPhojet7000AccCut.root","/home/marek/Analysis/Spectra/7TeVLHC10b/MCmodels/AnalyseFastPythia7000AccCutTune109.root","/home/marek/Analysis/Spectra/7TeVLHC10b/MCmodels/AnalyseFastPythia7000AccCutTune320.root"};
	TList* ktopi=new TList();
	TList* ptopi=new TList();
	//ktopi->ls();
	for(int i=0;i<3;i++)
		GetMCratios(mcfiles[i],ktopi,ptopi);
	//ktopi->ls();	
	TCanvas* c1MC= new TCanvas("K/pi","K/pi",1200,800);
	c1MC->cd();
	
	TLegend* Leg3 = new TLegend(0.1,0.35,0.45,0.75,"","NDC");
	Leg3->SetTextSize(0.027);
	Leg3->SetFillStyle(kFALSE);
	Leg3->SetLineColor(kWhite);
	Leg3->SetBorderSize(0);
	plotRatios(c1MC,combinepos[1],ktopi,Leg3,1);
	Leg3->Draw();
	TCanvas* c2MC= new TCanvas("p/pi","p/pi",1200,800);
	TLegend* Leg4 = new TLegend(0.1,0.25,0.45,0.65,"","NDC");
	Leg4->SetTextSize(0.027);
	Leg4->SetFillStyle(kFALSE);
	Leg4->SetLineColor(kWhite);
	Leg4->SetBorderSize(0);
	plotRatios(c2MC,combinepos[2],ptopi,Leg4,0);
	Leg4->Draw();
	ALICEWorkInProgress(c1MC,Form("%d/%d/%d",date.GetDay(),date.GetMonth(),date.GetYear()),"ALICE preliminary");
	ALICEWorkInProgress(c2MC,Form("%d/%d/%d",date.GetDay(),date.GetMonth(),date.GetYear()),"ALICE preliminary");

	c1MC->SaveAs("spectraKtopiMC.eps");
	c2MC->SaveAs("spectraptopiMC.eps");
	
	TCanvas* cK2pi=0x0;
	if(yileds[0]>0.0)
	{
		cK2pi=Ktopi(yileds[1]/yileds[0],yileds[1]/yileds[0]*TMath::Sqrt((yiledserr[0]/yileds[0])*(yiledserr[0]/yileds[0])+(yiledserr[1]/yileds[1])*(yiledserr[1]/yileds[1])));
		ALICEWorkInProgress(cK2pi,Form("%d/%d/%d",date.GetDay(),date.GetMonth(),date.GetYear()),"ALICE preliminary");
	}
	cK2pi->SaveAs("K2pi.eps");
	TCanvas* c1mul= new TCanvas("posmul","posmul",1200,800);
	
	c1mul->cd();
	TLegend* Leg5 = new TLegend(0.2,0.2,0.45,0.45,"","NDC");
	Leg5->SetNColumns(2);
	Leg5->SetTextSize(0.027);
	Leg5->SetFillStyle(kFALSE);
	Leg5->SetLineColor(kWhite);
	Leg5->SetBorderSize(0);
	DrawMulRatios(combineposbin1,combineposbin4,c1mul,Leg5,0);
	Leg5->Draw();
	TCanvas* c2mul= new TCanvas("negmul","negmul",1200,800);
	TLegend* Leg6 = new TLegend(0.2,0.2,0.45,0.45,"","NDC");
	Leg6->SetNColumns(2);
	Leg6->SetTextSize(0.027);
	Leg6->SetFillStyle(kFALSE);
	Leg6->SetLineColor(kWhite);
	Leg6->SetBorderSize(0);
	DrawMulRatios(combinenegbin1,combinenegbin4,c2mul,Leg6,1);
	Leg6->Draw();
	c2mul->cd();
	ALICEWorkInProgress(c1mul,Form("%d/%d/%d",date.GetDay(),date.GetMonth(),date.GetYear()),"ALICE preliminary");
	ALICEWorkInProgress(c2mul,Form("%d/%d/%d",date.GetDay(),date.GetMonth(),date.GetYear()),"ALICE preliminary");
	c1mul->SaveAs("MulPos.eps");
	c2mul->SaveAs("MulNeg.eps");
}
void ALICEWorkInProgress(TCanvas *c,TString today, TString label)
{
		c->cd();
	  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.72,0.72,0.89,0.89);
	  //TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.72,0.62,0.85,0.75);
	  myPadLogo->SetFillColor(0); 
	  myPadLogo->SetBorderMode(0);
	  myPadLogo->SetBorderSize(2);
	  myPadLogo->SetFrameBorderMode(0);
	  myPadLogo->SetLeftMargin(0.0);
	  myPadLogo->SetTopMargin(0.0);
	  myPadLogo->SetBottomMargin(0.0);
	  myPadLogo->SetRightMargin(0.0);
	  myPadLogo->SetFillStyle(0);
	  myPadLogo->Draw();
	  myPadLogo->cd();
	  TASImage *myAliceLogo = new TASImage("alice_logo_trans.png");
	  myAliceLogo->Draw();
	  c->cd();  
	  TPaveText* t1=new TPaveText(0.418103, 0.837798, 0.656609, 0.888393,"NDC");
	  t1->SetFillStyle(0);
	  t1->SetBorderSize(0);
	  t1->AddText(0.,0.,label);
	  t1->SetTextColor(kRed);
	  t1->SetTextFont(42);
	  t1->SetTextSize(0.04);
	  t1->Draw();
	  TPaveText* t2=new TPaveText(0.418103, 0.80, 0.656609, 0.84,"NDC");
	  t2->SetFillStyle(0);
	  t2->SetBorderSize(0);
	  t2->AddText(0.,0.,"pp at #sqrt{s} = 7 TeV");
	  t2->SetTextColor(kRed);
	  t2->SetTextFont(42);
	  t2->SetTextSize(0.027);
	  t2->Draw();
	  TPaveText* t3=new TPaveText(0.418103, 0.76, 0.656609, 0.80,"NDC");
	  t3->SetFillStyle(0);
	  t3->SetBorderSize(0);
	 // t3->AddText(0.,0.,"Statistical and systematic errors");
	  t3->AddText(0.,0.,"Statistical  errors");
	  t3->SetTextColor(kRed);
	  t3->SetTextFont(42);
	  t3->SetTextSize(0.027);
	  t3->Draw(); 
	   TPaveText* t2=new TPaveText(0.65,0.65,0.89,0.7,"NDC");
	   t2->SetFillStyle(0);
	   t2->SetBorderSize(0);
	   t2->SetTextColor(kRed);
	   t2->SetTextFont(52);
	   t2->AddText(0.,0.,today.Data());
	   t2->Draw();
}

TH1D* DrawHisto(AliOADBPWG2Spectra* fOADBSpectra, TCanvas* c1, Int_t ana,Int_t method,Int_t par,Int_t charge,Int_t color, Int_t marker,TLegend* leg,TString estimator="MB" ,Int_t bin=-1)
{
 
   TH1D * h = 0;
    TString binname=estimator;
   if(estimator.CompareTo("MB")==0)
   {
		cout<<fOADBSpectra->GetHistoName(ana,method,par,charge,"MB")<<endl;
		h = fOADBSpectra->GetHisto(ana,method,par,charge,"MB");
   // Draw the selected histogram
	}
	else
	{
		cout<<fOADBSpectra->GetHistoName(ana,method,par,charge,estimator.Data(),bin)<<endl;
		h = fOADBSpectra->GetHisto(ana,method,par,charge,estimator.Data(),bin);
		binname+=bin;	
	}
   if(!h) 
   {
     cout << "Cannot get pointer to histo" << endl;
     return 0x0;
   }
   h->GetYaxis()->SetRangeUser(0.001,10.0);
    h->GetXaxis()->SetRangeUser(0.0,3.0);
   
    if(charge==0)
		h->SetTitle(Form("Positive_%s",binname.Data()));
    else
		h->SetTitle(Form("Negative_%s",binname.Data()));
    h->SetXTitle("p_{T} (GeV/c)");
     h->SetYTitle("1/N_{events} dN/dp_{T} |y|<0.5");
 
   h->SetMarkerColor(color);
   h->SetMarkerStyle(marker);
   
   
   
   
 if(c1&&leg)
 {
	TString opt = "E1";
	c1->cd();
    c1->GetListOfPrimitives()->Print();
    if (c1->GetListOfPrimitives()->GetEntries()>0) 
		opt += "same";
	c1->Update();
	c1->Modified();
	c1->Update();
   
	leg->AddEntry(h,Form("%s_{%s_%s}",partext[charge*3+par],fgkDetectorNames[ana],fgkPidTypeNames[method]),"p");
    h->Draw(opt.Data());
   // TCanvas::Update() draws the frame, after which it can be changed
}
   return h;
}
Double_t myLevyPt(Double_t *pt, Double_t *par)
{
  
  Double_t pdNdy  = par[0];
  Double_t pTemp = par[1];
  Double_t pPower = par[2];
  Double_t pMass  = par[3];
  Double_t pBigCoef = ((pPower-1)*(pPower-2)) / (pPower*pTemp*(pPower*pTemp+pMass*(pPower-2)));
 
  Double_t pInPower = 1.0 + (TMath::Sqrt(pt[0]*pt[0]+pMass*pMass)-pMass) / (pPower*pTemp);
  return pdNdy * pt[0] * pBigCoef * TMath::Power(pInPower,(-1.0)*pPower);
}
TH1D* CombineSpectra(TList* l,Int_t startpoint, Float_t* min, Float_t* max)
{
	//l->ls();
//	((TH1D*)l->At(0))->Print("all");
	TH1D* hcombine=new TH1D(*((TH1D*)l->At(0)));
	hcombine->SetName(Form("%s_combine",((TH1D*)l->At(startpoint))->GetName()));
	//cout<<hcombine->GetName()<<endl;
	for (int i=1; i<=hcombine->GetNbinsX();i++)
	{
		Float_t pt=hcombine->GetXaxis()->GetBinCenter(i);
		Int_t use[4]={1,1,1,1};
		Float_t values[4]={0.0,0.0,0.0,0.0};
		Float_t errors[4]={-1.0,-1.0,-1.0,-1.0};
	//	cout<<pt<<endl;
		for(int j=0;j<4;j++)
		{
			if(min[3*j+startpoint]>0.0&&min[3*j+startpoint]>pt)
				use[j]=0;
			if(max[3*j+startpoint]>0.0&&max[3*j+startpoint]<pt)
				use[j]=0;
			if(use[j])
			{
				values[j]=((TH1D*)l->At(3*j+startpoint))->GetBinContent(i);
				if(values[j]>0.0)
					errors[j]=((TH1D*)l->At(3*j+startpoint))->GetBinError(i);
				else
				{
					values[j]=0.0;
					use[j]=0;	
				}		
			}	
		}	
		if(use[0]+use[1]+use[2]+use[3]==0)
		{
				hcombine->SetBinContent(i,0.0);
				hcombine->SetBinError(i,0.0);
		}	
		else
		{
			hcombine->SetBinContent(i,(values[0]+values[1]+values[2]+values[3])/(use[0]+use[1]+use[2]+use[3]));
			hcombine->SetBinError(i,TMath::MaxElement(4,errors));	
		}
		//cout<<use[0]+use[1]+use[2]+use[3]<<endl;
	}	
	return hcombine;
} 
void DrawandFit(TH1D* hist,TCanvas* c1, Int_t par , TCanvas* cratio,Float_t* yileds,Float_t* yiledserr)
{
	c1->cd();
	TF1 *f1 = new TF1(Form("%s_fit",hist->GetName()),myLevyPt,0.1,3.0,4);
	f1->SetParNames("dN/dy","T","n","mass");
	f1->FixParameter(3,AliPID::ParticleMass(par+2));
	f1->SetParameter(0,hist->GetMaximum()/2.0);
						
	f1->SetParameter(1,0.1);
	f1->SetParLimits(1,0.01,10.0);
	f1->SetParameter(2,3.0);
	f1->SetParLimits(2,2.1,1000.0);	
	
	hist->SetMarkerColor(colors[par]);
    hist->SetMarkerStyle(markers[par]);
    
    
    hist->Fit(Form("%s_fit",hist->GetName()),"I0");
    
    TString opt = "E1";
  
     c1->GetListOfPrimitives()->Print();
     if (c1->GetListOfPrimitives()->GetEntries()>0) 
		opt += "same";	
		cout<<opt<<endl;
	hist->DrawCopy(opt.Data());
	f1->DrawCopy("Lsame");	
	cratio->cd();
	Plotratiodatafit(hist,f1)->Draw(opt.Data());
	yileds[par]+=f1->GetParameter(0);
	yiledserr[par]=TMath::Sqrt(yiledserr[par]*yiledserr[par]+f1->GetParError(0)*f1->GetParError(0));
}	
void GetMCratios(TString filename, TList* ktopi,TList* ptopi)
{
		TFile* f=TFile::Open(filename.Data());
		if(!f)
					return;
		TH1F* fpiplus=(TH1F*)f->Get("h1PrimariesPiPlusPtVar");		
		TH1F* fpiminus=(TH1F*)f->Get("h1PrimariesPiMinusPtVar");	
		TH1F* fKplus=(TH1F*)f->Get("h1PrimariesKPlusPtVar");		
		TH1F* fKminus=(TH1F*)f->Get("h1PrimariesKMinusPtVar");	
		TH1F* fpplus=(TH1F*)f->Get("h1PrimariesProtonPtVar");		
		TH1F* fpminus=(TH1F*)f->Get("h1PrimariesAntiProtonPtVar");
		
		fpiplus->Sumw2();
		fpiminus->Sumw2();
		fKplus->Sumw2();
		fKminus->Sumw2();
		fpplus->Sumw2();
		fpminus->Sumw2();
		
		fpiplus->Add(fpiminus);
		fKplus->Add(fKminus);
		fpplus->Add(fpminus);	

		fKplus->Divide(fpiplus);
		fpplus->Divide(fpiplus);
		
		ktopi->Add(fKplus);
		ptopi->Add(fpplus);
}
void plotRatios(TCanvas* c1, TH1D* data , TList* mc,TLegend* leg, Int_t opt)
{
	c1->cd();
	data->SetMarkerStyle(markers[0]);
	data->SetMarkerColor(colors[0]);
	data->Draw("E1");
	leg->AddEntry(data,"data","p");
	if(opt==1)
		data->SetTitle("K/#pi");
	else
		data->SetTitle("p/#pi");		
	data->GetYaxis()->SetTitle("");
	data->GetYaxis()->SetRangeUser(0.0,0.5);
	data->GetXaxis()->SetRangeUser(0.0,2.5);
	
	for (int i=0;i<3;i++)
	{
			TH1F* mchist=(TH1F*)mc->At(i);
			leg->AddEntry(mchist,mcmodels[i].Data(),"p");
			mchist->SetMarkerStyle(markers[i+1]);
			mchist->SetMarkerColor(colors[i+1]);
			mchist->Draw("E1same");
	}
}	
TH1D* Plotratiodatafit(TH1D* hist1, TF1* fun)
{
	TString name(hist1->GetName());
	TH1D* histratio=new TH1D(*hist1);
	
	histratio->SetName(Form("ratio%s",name.Data()));
	//histratio->SetTitle("ratio");
	histratio->GetYaxis()->SetTitle("ratio data/fit");
	for (int i=1;i<histratio->GetNbinsX();i++)
	{
		if(histratio->GetBinContent(i)>0.0)
		{
			histratio->SetBinContent(i,histratio->GetBinContent(i)/fun->Eval(histratio->GetXaxis()->GetBinCenter(i)));
			histratio->SetBinError(i,histratio->GetBinError(i)/fun->Eval(histratio->GetXaxis()->GetBinCenter(i)));	
		}
	}
	histratio->GetYaxis()->SetRangeUser(0.8,1.2);
	return histratio;
}
TCanvas* Ktopi(Float_t value,Float_t err)
{
//=========Macro generated from canvas: c1/c1
//=========  (Thu Oct 14 14:18:11 2010) by ROOT version5.26/00b
   TCanvas *c1 = new TCanvas("K2pi", "K2pi",0,44,1105,782);
   c1->Range(0.4639814,-0.024875,3.950796,0.223875);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogx();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH2F *__1 = new TH2F("__1","",4002,1.5,9000.5,199,0,0.199);
   __1->SetDirectory(0);
   __1->SetStats(0);
   __1->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
   __1->GetXaxis()->SetRangeUser(6,9002);
   __1->GetXaxis()->CenterTitle(true);
   __1->GetXaxis()->SetTitleSize(0.05);
   __1->GetXaxis()->SetTitleOffset(0.87);
   __1->GetYaxis()->SetTitle("K/#pi");
   __1->GetYaxis()->CenterTitle(true);
   __1->GetYaxis()->SetTitleSize(0.06);
   __1->GetYaxis()->SetTitleOffset(0.67);
   __1->Draw("");
   
   TGraphErrors *gre = new TGraphErrors(1);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerColor(3);
   gre->SetMarkerStyle(29);
   gre->SetMarkerSize(3.1);
   gre->SetPoint(0,199.9004,0.1030843);
   gre->SetPointError(0,0,0.008);
   
   TH1F *Graph1 = new TH1F("Graph1","Graph",100,199.9,201.1);
   Graph1->SetMinimum(0.0934);
   Graph1->SetMaximum(0.1126);
   Graph1->SetDirectory(0);
   Graph1->SetStats(0);
   gre->SetHistogram(Graph1);
   
   gre->Draw("p");
   
   gre = new TGraphErrors(1);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerColor(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(2.5);
   //gre->SetPoint(0,897.8542,0.1235148);
   //gre->SetPointError(0,0,0.008624);
   //gre->SetPoint(0,897.8542,0.122);  // new vale from MF: K/pi ratio = 0.122 +- 0.01  Oct. 14, 2010
   //    gre->SetPointError(0,0,0.01);

   gre->SetPoint(0,897.8542,0.12294);  // new vale from MF: K/pi ratio = 0.12294 +- 0.012  Feb.. 21, 2011
   gre->SetPointError(0,0,0.012);
   
   gre->SetPoint(1,7000,value);  // new vale from MF: K/pi ratio = 0.12294 +- 0.012  Feb.. 21, 2011
   gre->SetPointError(1,0,err);

   
   TH1F *Graph2 = new TH1F("Graph2","Graph",100,899.9,901.1);
   Graph2->SetMinimum(0.1128512);
   Graph2->SetMaximum(0.1335488);
   Graph2->SetDirectory(0);
   Graph2->SetStats(0);
   gre->SetHistogram(Graph2);
   
   gre->Draw("p");
   
   gre = new TGraphErrors(1);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerColor(4);
   gre->SetMarkerStyle(22);
   gre->SetMarkerSize(2.2);
   gre->SetPoint(0,17.24717,0.08050324);
   gre->SetPointError(0,0,0.0034);
   
   TH1F *Graph3 = new TH1F("Graph3","Graph",100,17.2,18.4);
   Graph3->SetMinimum(0.07634);
   Graph3->SetMaximum(0.0845);
   Graph3->SetDirectory(0);
   Graph3->SetStats(0);
   gre->SetHistogram(Graph3);
   
   gre->Draw("p");
   
   gre = new TGraphErrors(4);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(28);
   gre->SetMarkerSize(2.5);
   gre->SetPoint(0,300,0.105);
   gre->SetPointError(0,0,0.02);
   gre->SetPoint(1,534.998,0.112045);
   gre->SetPointError(1,0,0.01);
   gre->SetPoint(2,994.3588,0.1041596);
   gre->SetPointError(2,0,0.008);
   gre->SetPoint(3,1795.007,0.1131203);
   gre->SetPointError(3,0,0.005);
   
   TH1F *Graph4 = new TH1F("Graph4","Graph",100,150,1950);
   Graph4->SetMinimum(0.081);
   Graph4->SetMaximum(0.129);
   Graph4->SetDirectory(0);
   Graph4->SetStats(0);
   gre->SetHistogram(Graph4);
   
   gre->Draw("p");
   
   gre = new TGraphErrors(1);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(24);
   gre->SetMarkerSize(2.6);
   gre->SetPoint(0,546.8309,0.09519885);
   gre->SetPointError(0,0,0.0114);
   
   TH1F *Graph5 = new TH1F("Graph5","Graph",100,544.9,546.1);
   Graph5->SetMinimum(0.08132);
   Graph5->SetMaximum(0.10868);
   Graph5->SetDirectory(0);
   Graph5->SetStats(0);
   gre->SetHistogram(Graph5);
   
   gre->Draw("p");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
  // c1->ToggleToolBar();
  return c1;
}
void DrawMulRatios(TH1D** bin1,TH1D** bin2,TCanvas* c,TLegend* leg,Int_t charge=0)
{
	c->cd();
	for(int i=0;i<3;i++)
	{
		bin1[i]->GetYaxis()->SetTitle("ratio over MB spectra");
		bin1[i]->GetYaxis()->SetRangeUser(0.0,4.0);
		bin1[i]->SetMarkerStyle(markers[i]);
		bin1[i]->SetMarkerColor(colors[i]);
		leg->AddEntry(bin1[i],Form("%s_{bin1}",partext[i+charge*3],"p"));
		if(i==0)
			bin1[i]->Draw("E1");
		else
			bin1[i]->Draw("E1same");
		bin2[i]->GetYaxis()->SetTitle("ratio over MB spectra");
		bin2[i]->GetYaxis()->SetRangeUser(0.0,4.0);
		bin2[i]->SetMarkerStyle(markers[i+3]);
		bin2[i]->SetMarkerColor(colors[i]);	
		bin2[i]->Draw("E1same")	;	
		leg->AddEntry(bin2[i],Form("%s_{bin4}",partext[i+charge*3],"p"));
	}
}	
