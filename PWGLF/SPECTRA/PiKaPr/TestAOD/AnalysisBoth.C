class AliSpectraBothHistoManager;
class AliSpectraBothEventCuts; 
class AliSpectraBothTrackCuts;
TString Charge[]={"Pos","Neg"};
TString Sign[]={"Plus","Minus"};
TString Particle[]={"Pion","Kaon","Proton"};
AliSpectraBothHistoManager* managerdata=0x0;
AliSpectraBothEventCuts* ecutsdata=0x0; 
AliSpectraBothTrackCuts* tcutsdata=0x0;
	
AliSpectraBothHistoManager* managermc=0x0;
AliSpectraBothEventCuts* ecutsmc=0x0; 
AliSpectraBothTrackCuts* tcutsmc=0x0;

Float_t TOFMatchingScalling[2]={-1,-1};
Int_t Color[3]={1,2,4};
Int_t Marker[6]={20,21,22,24,25,26};
Double_t Range[3]={0.3,0.3,0.5}; // LowPt range for pi k p

enum ECharge_t {
  kPositive,
  kNegative,
  kNCharges
};
enum {
 kdodca=0x1,
 kgeantflukaKaon=0x2,
 kgeantflukaProton=0x4,
 knormalizationtoeventspassingPhySel=0x8,
 kveretxcorrectionandbadchunkscorr=0x10,
 kmcisusedasdata=0x20,
 kdonotusedcacuts=0x40		
};	

Bool_t OpenFile(TString dirname, TString outputname, Bool_t mcflag,Bool_t mcasdata=false);
void AnalysisBoth (UInt_t options=0xF,TString outdate, TString outnamedata, TString outnamemc="" )
{
	TH1::AddDirectory(kFALSE);
	gSystem->Load("libCore.so");  
	gSystem->Load("libPhysics.so");
	gSystem->Load("libTree");
	gSystem->Load("libMatrix");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libOADB");
	gSystem->Load("libANALYSISalice");
	gSystem->Load("libTENDER");
	gSystem->Load("libCORRFW");
	gSystem->Load("libPWGTools");
	gSystem->Load("libPWGLFspectra");
  	
  	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD/QAPlotsBoth.C");
	Double_t mass[3];
	mass[0]   = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
	mass[1]   = TDatabasePDG::Instance()->GetParticle("K+")->Mass();
	mass[2] = TDatabasePDG::Instance()->GetParticle("proton")->Mass();

	TFormula* dcacutxy=0x0;
	if(!(options&kdonotusedcacuts))
	{
	
		AliESDtrackCuts* esdtrackcuts= AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
		TString formulastring(esdtrackcuts->GetMaxDCAToVertexXYPtDep());
		formulastring.ReplaceAll("pt","x");
		dcacutxy=new TFormula("dcacutxy",formulastring.Data());
	}
	TList* lout=new TList();


	TString indirname=Form("/output/train%s",outdate.Data());
	//TString indirname("/output/train24072012");
	if(outnamemc.Length()==0)
	outnamemc=outnamedata;
	cout<<indirname.Data()<<" "<<outnamemc.Data()<<endl;
	// Do the job 


	OpenFile(indirname,outnamemc,true);
	OpenFile(indirname,outnamedata,false,((Bool_t)(options&kmcisusedasdata)));
	if(!managermc||!managerdata)
	{
		cout<<managermc<<" "<<managerdata<<endl;
		return;	
	}
	TH1F* rawspectradata[6];
	TH1F* rawspectramc[6];
	TH1F* MCTruth[6];
	TH1F* eff[6];
	TH1F* contallMC[6];
	TH1F* contPID[6];
	TH1F* contWD[6];
	TH1F* contMat[6];
	TH1F* confinal[6];
	
	TH1F* contfit[12];
	TH1F* contWDfit[12];
	TH1F* contMatfit[12];
	TH1F* primaryfit[12];

	
	
	TH1F* spectra[6];
	TH1F* spectraLeonardo[6];
	
	TH1F* corrLeonardo[6]; 
	//GetSpectra(managerdata,rawspectradata,true);
	//GetSpectra(managermc,rawspectramc,true,true);
	
	GetPtHistFromPtDCAhisto("hHistPtRecSigma","SpectraMC",managermc,rawspectramc,dcacutxy);
	GetPtHistFromPtDCAhisto("hHistPtRecSigma","SpectraDATA",managerdata,rawspectradata,dcacutxy);
	GetPtHistFromPtDCAhisto("hHistPtRecTruePrimary","eff",managermc,eff,dcacutxy);
	GetPtHistFromPtDCAhisto("hHistPtRecTrue","conPID",managermc,contPID,dcacutxy);
	GetPtHistFromPtDCAhisto("hHistPtRecSigmaSecondaryWeakDecay","conWD",managermc,contWD,dcacutxy);
	GetPtHistFromPtDCAhisto("hHistPtRecSigmaSecondaryMaterial","conMat",managermc,contMat,dcacutxy);
	
	
//	Double_t neventsdata =  ecutsdata->NumberOfEvents();
	Double_t neventsmcall = 1 ;  //if loop over MC is done after or befor events cuts this will be changed 
	Double_t neventsdata =  1;
	Double_t neventsmc =  1;

	if(options&knormalizationtoeventspassingPhySel)
	{
		neventsmcall= ecutsmc->NumberOfProcessedEvents();
		 neventsdata=ecutsdata->NumberOfPhysSelEvents();
		 neventsmc=ecutsmc->NumberOfPhysSelEvents();
	}
	else
	{
		neventsdata=ecutsdata->NumberOfEvents();
		 neventsmc=ecutsmc->NumberOfEvents();
		neventsmcall= ecutsmc->NumberOfEvents();


	}
	GetMCTruth(MCTruth);
	
	
	
	TH1F* allgen=((TH1F*)managermc->GetPtHistogram1D("hHistPtGen",1,1))->Clone();
	allgen->SetName("AllGen");
	TH1F* allrecMC=GetOneHistFromPtDCAhisto("hHistPtRec","rawallMC",managermc,dcacutxy);
	TH1F* alleff=GetOneHistFromPtDCAhisto("hHistPtRecPrimary","effall",managermc,dcacutxy);
	TH1F* allrecdata=GetOneHistFromPtDCAhisto("hHistPtRec","rawalldata",managerdata,dcacutxy);
	
  	
	
	
	TH1F* spectraall=(TH1F*)allrecdata->Clone("recNch");
	TH1F* contall=(TH1F*)allrecMC->Clone("contall");
	contall->Add(alleff,-1);
	alleff->Divide(alleff,allgen,1,1,"B");
	contall->Divide(contall,allrecMC,1,1,"B");
	
	GetCorrectedSpectra(spectraall,allrecdata,alleff,contall);
	Divideby2pipt(spectraall);

	allrecdata->Scale(1./neventsdata,"width");
	allgen->Scale(1./neventsmcall,"width");
	allrecMC->Scale(1./neventsmc,"width");
	spectraall->Scale(1./neventsdata,"width");


	lout->Add(allgen);
	lout->Add(allrecMC);
	lout->Add(alleff);
	lout->Add(allrecdata);
	lout->Add(spectraall);
	lout->Add(contall);
	
	for (int i=0;i<6;i++)
	{
	
		
		TString tmpname(rawspectramc[i]->GetTitle());
		tmpname.ReplaceAll("SpectraMC","%s");
		contallMC[i]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contallMC")); 
		contfit[i]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contfit"));
		contWDfit[i]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contWDfit"));
		contMatfit[i]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contMatfit"));
		primaryfit[i]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"primaryfit"));
		
		contfit[i+6]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contfitonMC"));
		contWDfit[i+6]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contWDfitonMC"));
		contMatfit[i+6]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contMatfitonMC"));
		primaryfit[i+6]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"primaryfitMC"));
		
		contfit[i]->Reset();
		contWDfit[i]->Reset();
		contMatfit[i]->Reset();
		primaryfit[i]->Reset();
		

		contfit[i+6]->Reset();
		contWDfit[i+6]->Reset();
		contMatfit[i+6]->Reset();
		primaryfit[i+6]->Reset();
		
		SetBintoOne(primaryfit[i]);
		SetBintoOne(primaryfit[i+6]);		
		spectra[i]=(TH1F*)rawspectradata[i]->Clone(Form(tmpname.Data(),"SpectraFinal"));
		
		spectraLeonardo[i]=(TH1F*)rawspectradata[i]->Clone(Form(tmpname.Data(),"SpectraFinalLeonardo"));
		corrLeonardo[i]=(TH1F*)MCTruth[i]->Clone(Form(tmpname.Data(),"CorrFactLeonardo"));
		
		corrLeonardo[i]->Divide(corrLeonardo[i],rawspectramc[i],1,1,"B");
		
		
		
		contallMC[i]->Add(eff[i],-1.0);
		RecomputeErrors(contallMC[i]);
		contallMC[i]->Sumw2(); 
		contallMC[i]->Divide(contallMC[i],rawspectramc[i],1,1,"B");
		
		eff[i]->Divide(eff[i],MCTruth[i],1,1,"B");
		
		
		contPID[i]->Sumw2();
		rawspectramc[i]->Sumw2();
		contPID[i]->Add(contPID[i],rawspectramc[i],-1,1);
		RecomputeErrors(contPID[i]);
		contPID[i]->ResetStats();
		contPID[i]->Sumw2();
		contPID[i]->Divide(contPID[i],rawspectramc[i],1,1,"B");
		
		confinal[i]=(TH1F*)contPID[i]->Clone(Form(tmpname.Data(),"confinal"));


		contWD[i]->Divide(contWD[i],rawspectramc[i],1,1,"B");
		contMat[i]->Divide(contMat[i],rawspectramc[i],1,1,"B");
	
	
	
		rawspectradata[i]->Scale(1./neventsdata,"width");
		rawspectramc[i]->Scale(1./neventsmc,"width");
		MCTruth[i]->Scale(1./neventsmcall,"width");
		spectraLeonardo[i]->Scale(1./neventsdata,"width");
	
	
	
		lout->Add(rawspectradata[i]);
		lout->Add(rawspectramc[i]);
		lout->Add(MCTruth[i]);
		lout->Add(eff[i]);
		lout->Add(contallMC[i]);
		lout->Add(contPID[i]);
		lout->Add(contWD[i]);
		lout->Add(contMat[i]);
		lout->Add(contfit[i]);
		lout->Add(contWDfit[i]);
		lout->Add(contMatfit[i]);
 		lout->Add(primaryfit[i]);	
		lout->Add(contfit[i+6]);
		lout->Add(contWDfit[i+6]);
		lout->Add(contMatfit[i+6]);
		lout->Add(primaryfit[i+6]);	
		lout->Add(spectra[i]);
		lout->Add(spectraLeonardo[i]);
		lout->Add(confinal[i]);
	}
	TFile* fout=new TFile(Form("./results/ResMY_%s_%s.root",outnamemc.Data(),outdate.Data()),"RECREATE");
	if (options&kdodca)
		DCACorrectionMarek(managerdata,managermc,dcacutxy,fout,contfit,contWDfit,contMatfit,primaryfit);
	for (int i=0;i<6;i++)
	{
			if(options&kdodca)
			{
				confinal[i]->Add(contfit[i]);
				GetCorrectedSpectra(spectra[i],rawspectradata[i],eff[i],confinal[i]);
			}
			else
			{
				GetCorrectedSpectra(spectra[i],rawspectradata[i],eff[i],contallMC[i]);	
			}
			GetCorrectedSpectraLeonardo(spectraLeonardo[i],corrLeonardo[i],primaryfit[i],primaryfit[i+6]);
			CleanHisto(spectra[i],-1,100,contPID[i]);
			CleanHisto(spectraLeonardo[i],-1,100,contPID[i]);				
	}
	
	GFCorrection(spectra,tcutsdata->GetPtTOFMatching(),options);
	GFCorrection(spectraLeonardo,tcutsdata->GetPtTOFMatching(),options);
	 MatchingTOFEff(spectra,lout);
	  MatchingTOFEff(spectraLeonardo);
	TH1F* allch=GetSumAllCh(spectra,mass);
	lout->Add(allch);	

//	lout->ls();
	fout->cd();	
	TList* listqa=new TList();
	TList* canvaslist=new TList();
	Float_t vertexcorrection=1.0;
	Float_t corrbadchunksvtx=QAPlotsBoth(managerdata,managermc,ecutsdata,ecutsmc,tcutsdata,tcutsmc,listqa,canvaslist);
	if (options&kveretxcorrectionandbadchunkscorr)
		vertexcorrection=corrbadchunksvtx;
	cout<<" VTX corr="<<vertexcorrection<<endl;
	Double_t ycut=tcutsdata->GetY();
	if(TMath::Abs(ycut)>0.0)
		vertexcorrection=vertexcorrection/(2.0*ycut);
	for (int i=0;i<6;i++)
	{
		spectra[i]->Scale(vertexcorrection);
		spectraLeonardo[i]->Scale(vertexcorrection);
		if(TMath::Abs(ycut)>0.0)
		{
			rawspectradata[i]->Scale(1.0/(2.0*ycut));
                	rawspectramc[i]->Scale(1.0/(2.0*ycut));
                	MCTruth[i]->Scale(1.0/(2.0*ycut));
		}
	}	
	allch->Scale(vertexcorrection);
	spectraall->Scale(vertexcorrection/1.6);

	//spectraall->Scale(1.0/1.6);
	lout->Write("output",TObject::kSingleKey);	
	listqa->Write("outputQA",TObject::kSingleKey);
	canvaslist->Write("outputcanvas",TObject::kSingleKey);

	fout->Close();

}

Bool_t   OpenFile(TString dirname,TString outputname, Bool_t mcflag, Bool_t mcasdata)
{
	

	TString nameFile = Form("./%s/AnalysisResults%s.root",dirname.Data(),(mcflag?"MC":"DATA"));
	TFile *file = TFile::Open(nameFile.Data());
	if(!file)
	{
		cout<<"no file"<<endl;
		return false;
	}	
	TString sname=Form("OutputBothSpectraTask_%s_%s",(mcflag?"MC":"Data"),outputname.Data());
	if(mcasdata)
	{
		cout<<"using MC as data "<<endl;
		sname=Form("OutputBothSpectraTask_%s_%s","MC",outputname.Data());
	}
	file->ls();
	TDirectoryFile *dir=(TDirectoryFile*)file->Get(sname.Data());
	if(!dir)
	{
	//	cout<<"no dir "<<sname.Data()<<endl;	
		if(mcasdata)
		{
			cout<<"using MC as data "<<endl;
			sname=Form("OutputAODSpectraTask_%s_%s","MC",outputname.Data());
		}
		else	
			sname=Form("OutputAODSpectraTask_%s_%s",(mcflag?"MC":"Data"),outputname.Data());
	//	cout<<"trying "<<sname.Data()<<endl;
		dir=(TDirectoryFile*)file->Get(sname.Data());
		if(!dir)
		{
			cout<<"no dir "<<sname.Data()<<endl;
			return false;
		}
	}
	cout << " -- Info about " <<(mcflag?"MC":"DATA") <<" -- "<< endl;
	if(mcflag)
	{
		managermc= (AliSpectraBothHistoManager*) dir->Get("SpectraHistos");	
		ecutsmc = (AliSpectraBothEventCuts*) dir->Get("Event Cuts");
		tcutsmc = (AliSpectraBothTrackCuts*) dir->Get("Track Cuts");
		ecutsmc->PrintCuts();
		tcutsmc->PrintCuts();
		if(!managermc||!ecutsmc||!tcutsmc)
			return false;
	}
	else
	{
		managerdata= (AliSpectraBothHistoManager*) dir->Get("SpectraHistos");	
		ecutsdata = (AliSpectraBothEventCuts*) dir->Get("Event Cuts");
		tcutsdata = (AliSpectraBothTrackCuts*) dir->Get("Track Cuts");
		ecutsdata->PrintCuts();
		tcutsdata->PrintCuts();
		if(!managerdata||!ecutsdata||!tcutsdata)
			return false;
	}
	return true;
}

 void GetMCTruth(TH1F** MCTruth)
 {
	for(Int_t icharge=0;icharge<2;icharge++)
	{
		for(Int_t ipart=0;ipart<3;ipart++)
		{
			Int_t index=ipart+3*icharge;
			TString hname=Form("hHistPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
			MCTruth[index]=(TH1F*)((TH1F*)managermc->GetPtHistogram1D(hname.Data(),1,1))->Clone();
			MCTruth[index]->SetName(Form("MCTruth_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
			MCTruth[index]->SetTitle(Form("MCTruth_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
			MCTruth[index]->Sumw2(); 
		}
	}
}

TH1F* GetOneHistFromPtDCAhisto(TString name,TString hnameout,AliSpectraBothHistoManager* hman,TFormula* dcacutxy)
{
			histo =(TH1F*)((TH1F*) hman->GetPtHistogram1D(name.Data(),-1,-1))->Clone();
			histo->SetName(hnameout.Data());
			histo->SetTitle(hnameout.Data());
		  
			if(dcacutxy)
			{
				for(int ibin=1;ibin<histo->GetNbinsX();ibin++)
				{
					Double_t lowedge=histo->GetBinLowEdge(ibin);
					Float_t cut=dcacutxy->Eval(lowedge);
					TH1F* dcahist=(TH1F*)hman->GetDCAHistogram1D(name.Data(),lowedge,lowedge));
					Float_t inyield=dcahist->Integral(dcahist->GetXaxis()->FindBin(-1.0*cut),dcahist->GetXaxis()->FindBin(cut));
					cout<<"corr data "<<histo->GetBinContent(ibin)<<" "<<inyield<<" "<<dcahist->Integral()<<" "<<hnameout.Data()<<endl;
					cout<<"test dca "<<lowedge<<" "<<dcacutxy->Eval(lowedge)<<" "<<dcacutxy->Eval(histo->GetXaxis()->GetBinUpEdge(ibin))<<" "<<dcahist->GetBinLowEdge(dcahist->GetXaxis()->FindBin(-1.0*cut))<<" "<<dcahist->GetXaxis()->GetBinUpEdge(dcahist->GetXaxis()->FindBin(-1.0*cut))<<endl;
					histo->SetBinContent(ibin,inyield);
					histo->SetBinError(ibin,TMath::Sqrt(inyield));
				}
			}
			histo->Sumw2();
			return histo;
}



	
void GetPtHistFromPtDCAhisto(TString hnamein, TString hnameout, AliSpectraBothHistoManager* hman,TH1F** histo,TFormula* dcacutxy)
{
	Float_t min[3]={0.3,0.3,0.4};
	Float_t max[3]={1.5,1.2,2.2};
	for(Int_t icharge=0;icharge<2;icharge++)
	{
		for(Int_t ipart=0;ipart<3;ipart++)
		{
			Int_t index=ipart+3*icharge;
			//TString hname=Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data());
			printf("Getting %s",hnamein.Data());
			TString nameinfinal=Form("%s%s%s",hnamein.Data(),Particle[ipart].Data(),Sign[icharge].Data());
			TString nameoutfinal=Form("%s%s%s",hnameout.Data(),Particle[ipart].Data(),Sign[icharge].Data());
			
			
			histo[index]=GetOneHistFromPtDCAhisto(nameinfinal,nameoutfinal,hman,dcacutxy);
			/*histo[index] =(TH1F*)((TH1F*) hman->GetPtHistogram1D(Form("%s%s%s",hnamein.Data(),Particle[ipart].Data(),Sign[icharge].Data()),-1,-1))->Clone();
			histo[index]->SetName(Form("%s%s%s",hnameout.Data(),Particle[ipart].Data(),Sign[icharge].Data()));
			histo[index]->SetTitle(Form("%s%s%s",hnameout.Data(),Particle[ipart].Data(),Sign[icharge].Data()));
		  
			if(dcacutxy)
			{
				for(int ibin=1;ibin<histo[index]->GetNbinsX();ibin++)
				{
					Double_t lowedge=histo[index]->GetBinLowEdge(ibin);
					Float_t cut=dcacutxy->Eval(lowedge);
					TH1F* dcahist=(TH1F*)hman->GetDCAHistogram1D(Form("%s%s%s",hnamein.Data(),Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge));
					Float_t inyield=dcahist->Integral(dcahist->GetXaxis()->FindBin(-1.0*cut),dcahist->GetXaxis()->FindBin(cut));
					cout<<"corr data "<<histo[index]->GetBinContent(ibin)<<" "<<inyield<<" "<<dcahist->Integral()<<endl;
					histo[index]->SetBinContent(ibin,inyield);
					histo[index]->SetBinError(ibin,TMath::Sqrt(inyield));
				}
			}*/
			CleanHisto(histo[index],min[ipart],max[ipart]);
			// histo[index]->Sumw2();
		}
	} 
}
void CleanHisto(TH1F* h, Float_t minV, Float_t maxV,TH1* contpid=0x0)
{
	for (int i=0;i<=h->GetNbinsX();i++)
	{	
		if(h->GetXaxis()->GetBinCenter(i)<minV||h->GetXaxis()->GetBinCenter(i)>maxV)
		{
			h->SetBinContent(i,0);
			h->SetBinError(i,0);
		}	
		if(contpid)
		{
			if(contpid->GetBinContent(i)>0.2)
			{
				h->SetBinContent(i,0);
				h->SetBinError(i,0);
			}
		}
	}
}


void DCACorrectionMarek(AliSpectraBothHistoManager* hman_data, AliSpectraBothHistoManager* hman_mc,TFormula* fun,TFile *fout,TH1F** hcon,TH1F** hconWD,TH1F** hconMat,TH1F** hprimary)
{
  printf("\n\n-> DCA Correction");  
  Double_t FitRange[2]={-2.0,2.0};
  Double_t CutRange[2]={-3.0,3.0};
  Double_t minptformaterial[6]={0.0,100.0,0.0,0.0,100.0,0.0};
  Double_t maxptformaterial[6]={0.0,-100.0,1.3,0.0,-100.0,0.0};
  Bool_t usefit[3]={true,false,true};
  Double_t maxbinforfit[6]={1.5,0,2.0,1.5,0,2.0};
  Printf("\DCACorr");
 // TH1F *hcorrection[2];
 /* TCanvas *ccorrection=new TCanvas("DCAcorrection","DCAcorrection",700,500);
  TCanvas *cRatiocorrection=new TCanvas("DCARatiocorrection","DCARatiocorrection",700,500);
  cRatiocorrection->Divide(2,1);
  ccorrection->Divide(2,1);*/
  TString sample[2]={"data","mc"};
  ofstream debug("debugDCA.txt");
  TList* listofdcafits=new TList();
  for(Int_t icharge=0;icharge<2;icharge++)
  {
		for(Int_t ipart=0;ipart<3;ipart++)
		{
			Int_t index=ipart+3*icharge;
			for(Int_t isample=0;isample<2;isample++)
			{

				/*hcorrection[isample]=(TH1F*)Spectra[index]->Clone();
				hcorrection[isample]->Reset("all");*/
				for(Int_t ibin_data=7;ibin_data<40;ibin_data++)
				{	
						
					Double_t lowedge=hcon[index]->GetBinLowEdge(ibin_data);
					Double_t binwidth=hcon[index]->GetBinWidth(ibin_data);
					debug<<"NEW "<<Particle[ipart].Data()<<" "<<Sign[icharge].Data()<<" "<<lowedge<<endl;
					if(fun)
					{						
						CutRange[1]=fun->Eval(lowedge);
						CutRange[0]=-1.0*CutRange[1];
					}	
					debug<<"cut  "<<CutRange[1]<<" "<<CutRange[0]<<endl;		
					Bool_t useMaterial=kFALSE;
					cout<<"try to fit "<<lowedge<<" "<<maxbinforfit[index]<<endl;
					if(lowedge>maxbinforfit[index])
						continue;
					if(lowedge>minptformaterial[index]&&lowedge<maxptformaterial[index])
						useMaterial=kTRUE;
	  
					TCanvas *cDCA=new TCanvas(Form("cDCA%d%s%s%sbin%d",index,sample[isample].Data(),Particle[ipart].Data(),Sign[icharge].Data(),ibin_data),Form("cDCA%d%s%s%sbin%d",index,sample[isample].Data(),Particle[ipart].Data(),Sign[icharge].Data(),ibin_data),1700,1500);
					if(isample==0)
						TH1F *hToFit =(TH1F*) ((TH1F*)hman_data->GetDCAHistogram1D(Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge))->Clone();
					if(isample==1)
						TH1F *hToFit =(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge))->Clone();
					debug<<Particle[ipart].Data()<<" "<<Sign[icharge].Data()<<" "<<lowedge<<endl;
					TH1F *hmc1=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecSigmaPrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge))->Clone();
					TH1F *hmc2=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecSigmaSecondaryWeakDecay%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge))->Clone();
					TH1F *hmc3=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecSigmaSecondaryMaterial%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge))->Clone();
					Double_t minentries=1;
					debug<<" Entries "<<isample<<" "<<hToFit->GetEntries()<<" "<<hmc1->GetEntries()<<" "<<hmc2->GetEntries()<<" "<<hmc3->GetEntries()<<endl;
					debug<<"2 Entries "<<isample<<" "<<hToFit->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))<<" "<<hmc1->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))<<" "<<hmc2->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))<<" "<<hmc3->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))<<endl;
					debug<< CutRange[0]<<" "<<CutRange[1]<<" "<<lowedge<<endl;
					if(hToFit->GetEntries()<=minentries || hmc1->GetEntries()<=minentries || hmc2->GetEntries()<=minentries || hmc3->GetEntries()<=minentries)
						continue;
					//Data and MC can have different stat
					hToFit->Sumw2();
					hmc1->Sumw2();
					hmc2->Sumw2();
					hmc3->Sumw2();
					
					Float_t corrforrebinning[4]={1.0,1.0,1.0,1.0};	
					
	
					if(hmc3->GetNbinsX()>300)
					{
					
						corrforrebinning[0]=hToFit->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]));
						corrforrebinning[1]=hmc1->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]));
						corrforrebinning[2]=hmc2->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]));
						corrforrebinning[3]=hmc3->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]));

						hToFit->Rebin(30);
						hmc1->Rebin(30);
						hmc2->Rebin(30);
						hmc3->Rebin(30);
						if(hToFit->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))>0.0)
							corrforrebinning[0]=corrforrebinning[0]/hToFit->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]));
						else	
							corrforrebinning[0]=1.0;
						if(hmc1->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))>0.0)
							corrforrebinning[1]=corrforrebinning[1]/hmc1->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]));
						else	
							corrforrebinning[1]=1.0;
						if(hmc2->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))>0.0)
							corrforrebinning[2]=corrforrebinning[2]/hmc2->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]));
						else	
							corrforrebinning[2]=1.0;
						if(hmc3->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))>0.0)
							corrforrebinning[3]=corrforrebinning[3]/hmc3->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]));
						else	
							corrforrebinning[3]=1.0;
							
						debug<<" cor bin "<<corrforrebinning[0]<<" "<<corrforrebinning[1]<<" "<<corrforrebinning[2]<<" "<<corrforrebinning[3]<<endl;


					}

					cDCA->cd();
					gPad->SetGridy();
					gPad->SetGridx();
					gPad->SetLogy();
	 
					TObjArray *mc=0x0;
					if(useMaterial)
						mc = new TObjArray(3);        // MC histograms are put in this array
					else
						mc = new TObjArray(2);
					mc->Add(hmc1);
					mc->Add(hmc2);
					if(useMaterial)
						mc->Add(hmc3);
					TFractionFitter* fit = new TFractionFitter(hToFit,mc); // initialise
					fit->Constrain(0,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
					fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
					if(useMaterial)
						fit->Constrain(2,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
					fit->SetRangeX(hToFit->GetXaxis()->FindBin(FitRange[0]),hToFit->GetXaxis()->FindBin(FitRange[1]));
					hToFit->GetXaxis()->SetRange(hToFit->GetXaxis()->FindBin(FitRange[0]),hToFit->GetXaxis()->FindBin(FitRange[1]));
					hToFit->SetTitle(Form("DCA distr - %s %s %s %lf",sample[isample].Data(),Particle[ipart].Data(),Sign[icharge].Data(),lowedge));
					Int_t status = fit->Fit();               // perform the fit
					cout << "fit status: " << status << endl;
					debug<<"fit status: " << status << endl;
		
					if (status == 0 && usefit[ipart])
					{ 	                   // check on fit status
						TH1F* result = (TH1F*) fit->GetPlot();
						TH1F* PrimMCPred=(TH1F*)fit->GetMCPrediction(0);
						TH1F* secStMCPred=(TH1F*)fit->GetMCPrediction(1);
					
						TH1F* secMCPred=0x0;
						if(useMaterial)
							secMCPred=(TH1F*)fit->GetMCPrediction(2);
	    
						Double_t v1=0,v2=0,v3=0;
						Double_t ev1=0,ev2=0,ev3=0;
						Double_t cov=0.0;
						//first method, use directly the fit result
						fit->GetResult(0,v1,ev1);
						fit->GetResult(1,v2,ev2);
						if(useMaterial)
						{
							fit->GetResult(2,v3,ev3);
							fit->GetFitter()->GetCovarianceMatrixElement(1,2);
						}
						debug<<v1<<" "<<ev1<<" "<<v2<<" "<<ev2<<" "<<v3<<" "<<ev3<<" "<<endl;
					
	    
	   
						Float_t normalizationdata=hToFit->Integral(hToFit->GetXaxis()->FindBin(CutRange[0]),hToFit->GetXaxis()->FindBin(CutRange[1]))/hToFit->Integral(hToFit->GetXaxis()->FindBin(FitRange[0]),hToFit->GetXaxis()->FindBin(FitRange[1]));
						
						Float_t normalizationmc1=hmc1->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))/hmc1->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1]))/normalizationdata;
						Float_t normalizationmc2=hmc2->Integral(hmc2->GetXaxis()->FindBin(CutRange[0]),hmc2->GetXaxis()->FindBin(CutRange[1]))/hmc2->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc2->GetXaxis()->FindBin(FitRange[1]))/normalizationdata;
						Float_t normalizationmc3=hmc3->Integral(hmc3->GetXaxis()->FindBin(CutRange[0]),hmc3->GetXaxis()->FindBin(CutRange[1]))/hmc3->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc3->GetXaxis()->FindBin(FitRange[1]))/normalizationdata;
						debug<<"After Nor"<<endl;
						debug<<v1*normalizationmc1<<" "<<ev1*normalizationmc1<<" "<<v2*normalizationmc2<<" "<<ev2*normalizationmc2<<" "<<v3*normalizationmc3<<" "<<ev3*normalizationmc3<<" "<<endl;
						debug<<1.0-v1*normalizationmc1<<" "<<ev1*normalizationmc1<<" "<<v2*normalizationmc2+v3*normalizationmc3<<" "<<TMath::Sqrt(ev2*ev2*normalizationmc2*normalizationmc2+ev3*ev3*normalizationmc3*normalizationmc3+cov*normalizationmc3*normalizationmc2)<<endl;
						debug<<"addtional info"<<endl;
						Float_t normalizationmc1b=(hmc1->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))/hmc1->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1])))/normalizationdata;

						debug<<normalizationmc1<<" "<<normalizationmc1b<<endl;
						debug<<hmc1->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))/hmc1->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1]))<<" "<<secStMCPred->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))/secStMCPred->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1]))<<endl;
						//debug<<hmc1->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1]))<<" "<<secStMCPred->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1]))<<endl;
						debug<<hmc2->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))/hmc2->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1]))<<" "<<PrimMCPred->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))/PrimMCPred->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1]))<<endl;
						//debug<<hmc2->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1]))<<" "<<PrimMCPred->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1]))<<endl;
						if(PrimMCPred->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))>0.0)
						{
						if(useMaterial)
						{
							debug<<" ITSsa "<<endl;
							debug<<PrimMCPred->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))/(PrimMCPred->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))+secStMCPred->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1]))+secMCPred->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1])))<<endl;
						}			
						else
						{
							debug<<" ITSsa "<<endl;
							debug<<PrimMCPred->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))/(PrimMCPred->Integral(hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1]))+secStMCPred->Integral(hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1])))<<endl;
		
						}
						}
						else 
						{
							debug<<" ITSsa "<<endl;
							debug<<" NO "<<endl;
						}



						Float_t normalizationdata1=result->Integral(result->GetXaxis()->FindBin(CutRange[0]),result->GetXaxis()->FindBin(CutRange[1]))/result->Integral(result->GetXaxis()->FindBin(FitRange[0]),result->GetXaxis()->FindBin(FitRange[1]));
						

						normalizationdata1*=corrforrebinning[0];


						Float_t normalizationmc11=PrimMCPred->Integral(PrimMCPred->GetXaxis()->FindBin(CutRange[0]),PrimMCPred->GetXaxis()->FindBin(CutRange[1]))/PrimMCPred->Integral(PrimMCPred->GetXaxis()->FindBin(FitRange[0]),PrimMCPred->GetXaxis()->FindBin(FitRange[1]))/normalizationdata1;
						Float_t normalizationmc21=secStMCPred->Integral(secStMCPred->GetXaxis()->FindBin(CutRange[0]),secStMCPred->GetXaxis()->FindBin(CutRange[1]))/secStMCPred->Integral(secStMCPred->GetXaxis()->FindBin(FitRange[0]),secStMCPred->GetXaxis()->FindBin(FitRange[1]))/normalizationdata1;
						Float_t normalizationmc31=0;
						if(useMaterial)
							normalizationmc31=secMCPred->Integral(secMCPred->GetXaxis()->FindBin(CutRange[0]),secMCPred->GetXaxis()->FindBin(CutRange[1]))/secMCPred->Integral(secMCPred->GetXaxis()->FindBin(FitRange[0]),secMCPred->GetXaxis()->FindBin(FitRange[1]))/normalizationdata1;
						
						normalizationmc11*=corrforrebinning[1];
						normalizationmc21*=corrforrebinning[2];
						normalizationmc31*=corrforrebinning[3];

				debug<<"After Nor 2"<<endl;
						debug<<v1*normalizationmc11<<" "<<ev1*normalizationmc11<<" "<<v2*normalizationmc21<<" "<<ev2*normalizationmc21<<" "<<v3*normalizationmc31<<" "<<ev3*normalizationmc31<<endl;
						
						debug<<1.0-v1*normalizationmc11<<" "<<ev1*normalizationmc11<<" "<<v2*normalizationmc21+v3*normalizationmc31<<" "<<TMath::Sqrt(ev2*ev2*normalizationmc21*normalizationmc21+ev3*ev3*normalizationmc31*normalizationmc31+cov*normalizationmc31*normalizationmc21)<<endl;
					
						debug<<CutRange[0]<<" "<<CutRange[1]<<endl;		
						debug<<" Entries "<<isample<<" "<<hToFit->GetEntries()<<" "<<hmc1->GetEntries()<<" "<<hmc2->GetEntries()<<" "<<hmc3->GetEntries()<<endl;
						debug<<"2 Entries "<<isample<<" "<<hToFit->Integral(result->GetXaxis()->FindBin(CutRange[0]),result->GetXaxis()->FindBin(CutRange[1]))<<" "<<hmc1->Integral(result->GetXaxis()->FindBin(CutRange[0]),result->GetXaxis()->FindBin(CutRange[1]))<<" "<<hmc2->Integral(result->GetXaxis()->FindBin(CutRange[0]),result->GetXaxis()->FindBin(CutRange[1]))<<" "<<hmc3->Integral(result->GetXaxis()->FindBin(CutRange[0]),result->GetXaxis()->FindBin(CutRange[1]))<<endl;

						
						hconWD[index+6*isample]->SetBinContent(ibin_data,v2*normalizationmc21);
						hconWD[index+6*isample]->SetBinError(ibin_data,ev2*normalizationmc21);
						hconMat[index+6*isample]->SetBinContent(ibin_data,v3*normalizationmc31);
						hconMat[index+6*isample]->SetBinError(ibin_data,ev3*normalizationmc31);
						hprimary[index+6*isample]->SetBinContent(ibin_data,v1*normalizationmc11);
						hprimary[index+6*isample]->SetBinError(ibin_data,ev1*normalizationmc11);
						if(useMaterial)
						{
							hcon[index+6*isample]->SetBinContent(ibin_data,v2*normalizationmc21+v3*normalizationmc31);
							hcon[index+6*isample]->SetBinError(ibin_data,TMath::Sqrt(ev2*ev2*normalizationmc21*normalizationmc21+ev3*ev3*normalizationmc31*normalizationmc31+cov*normalizationmc31*normalizationmc21));
						}
						else
						{
							hcon[index+6*isample]->SetBinContent(ibin_data,v2*normalizationmc21);
							hcon[index+6*isample]->SetBinError(ibin_data,ev2*normalizationmc21);
						}
						
						
						
						//Drawing section
						result->Scale(1.0/result->Integral(result->GetXaxis()->FindBin(FitRange[0]),result->GetXaxis()->FindBin(FitRange[1])));
						hToFit->Scale(1.0/hToFit->Integral(hToFit->GetXaxis()->FindBin(FitRange[0]),hToFit->GetXaxis()->FindBin(FitRange[1])));
						PrimMCPred->Scale(v1/PrimMCPred->Integral(PrimMCPred->GetXaxis()->FindBin(FitRange[0]),PrimMCPred->GetXaxis()->FindBin(FitRange[1])));
						secStMCPred->Scale(v2/secStMCPred->Integral(secStMCPred->GetXaxis()->FindBin(FitRange[0]),secStMCPred->GetXaxis()->FindBin(FitRange[1])));
						if(useMaterial)
							secMCPred->Scale(v3/secMCPred->Integral(secMCPred->GetXaxis()->FindBin(FitRange[0]),secMCPred->GetXaxis()->FindBin(FitRange[1])));	    
						   
						result->SetLineColor(kBlack);
						PrimMCPred->SetLineColor(kGreen+2);
						secStMCPred->SetLineColor(kRed);
						
						hToFit->SetMinimum(0.0001);
						hToFit->DrawClone("E1x0");
						result->SetTitle("Fit result");
						result->DrawClone("lhistsame");
						PrimMCPred->DrawClone("lhistsame");
						secStMCPred->DrawClone("lhistsame");
						if(useMaterial)
						{
							secMCPred->SetLineColor(kBlue);	
							secMCPred->DrawClone("lhistsame");
						}
					}
					else
					{
						hconWD[index+6*isample]->SetBinContent(ibin_data,0.0);
						hconWD[index+6*isample]->SetBinError(ibin_data,0.0);
						hconMat[index+6*isample]->SetBinContent(ibin_data,0.0);
						hconMat[index+6*isample]->SetBinError(ibin_data,0.0);					
						hcon[index+6*isample]->SetBinContent(ibin_data,0.0);
						hcon[index+6*isample]->SetBinError(ibin_data,0.0);
						hprimary[index+6*isample]->SetBinContent(ibin_data,1.0);
						hprimary[index+6*isample]->SetBinError(ibin_data,0.0);
					}
					listofdcafits->Add(cDCA);
					
					//cDCA->Write();
					delete hToFit;
				}
	

			}

		}
	}
	fout->cd();
	listofdcafits->Write("DCAfits",TObject::kSingleKey);	
}

void RecomputeErrors(TH1* h)
{
	for (int i=0; i<=h->GetXaxis()->GetNbins(); i++)
		h->SetBinError(i,TMath::Sqrt(h->GetBinContent(i)));
	h->Sumw2(); 	
}
void SetBintoOne(TH1* h)
{
	for (int i=0;i<=h->GetXaxis()->GetNbins(); i++) 
	{
		h->SetBinContent(i,1);
		h->SetBinError(i,0);
	}
}


void GetCorrectedSpectra(TH1F* corr,TH1F* raw,TH1F* eff, TH1F* con)
{
	for (int i=0;i<=corr->GetXaxis()->GetNbins(); i++) 
	{
		corr->SetBinContent(i,1);
		corr->SetBinError(i,0);
	}
	corr->Sumw2(); 
	corr->Add(con,-1);
	corr->Sumw2(); 	
	corr->Divide(eff);
	corr->Sumw2(); 
	corr->Multiply(raw);
	corr->Sumw2(); 
}
void GetCorrectedSpectraLeonardo(TH1F* spectra,TH1F* correction, TH1F* hprimaryData,TH1F* hprimaryMC)
{
	spectra->Sumw2(); 
	spectra->Multiply(correction);
	spectra->Sumw2(); 
	hprimaryData->Sumw2(); 
	spectra->Multiply(hprimaryData);
	hprimaryMC->Sumw2(); 
	spectra->Divide(hprimaryMC);
}

void GFCorrection(TH1F **Spectra,Float_t tofpt,UInt_t options)
{
	if (options&kgeantflukaKaon)
	{		
	 	 //Geant/Fluka Correction
	  	Printf("\nGF correction for Kaons");
	  	//Getting GF For Kaons in TPC
	  	TGraph *gGFCorrectionKaonPlus=new TGraph();
	  	gGFCorrectionKaonPlus->SetName("gGFCorrectionKaonPlus");
	  	gGFCorrectionKaonPlus->SetTitle("gGFCorrectionKaonPlus");
	  	TGraph *gGFCorrectionKaonMinus=new TGraph();
	  	gGFCorrectionKaonMinus->SetName("gGFCorrectionKaonMinus");
	  	gGFCorrectionKaonMinus->SetTitle("gGFCorrectionKaonMinus");
	 	 TString fnameGeanFlukaK="GFCorrection/correctionForCrossSection.321.root";
  	  	TFile *fGeanFlukaK= new TFile(fnameGeanFlukaK.Data());
	  	if (!fGeanFlukaK)		
			return;
	  	TH1F *hGeantFlukaKPos=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionParticles");
	  	TH1F *hGeantFlukaKNeg=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionAntiParticles");
	  	//getting GF func for Kaons with TOF
	  	TF1 *fGFKPosTracking;
	  	fGFKPosTracking = TrackingEff_geantflukaCorrection(3,kPositive);
	  	TF1 *fGFKNegTracking;
	 	 fGFKNegTracking = TrackingEff_geantflukaCorrection(3,kNegative);
	 	 TF1 *fGFKPosMatching;
	 	 fGFKPosMatching = TOFmatchMC_geantflukaCorrection(3,kPositive);
	 	 TF1 *fGFKNegMatching;
	 	 fGFKNegMatching = TOFmatchMC_geantflukaCorrection(3,kNegative);
	 	 for(Int_t binK=0;binK<=Spectra[1]->GetNbinsX();binK++)
	  	{
			if(Spectra[1]->GetBinCenter(binK)<tofpt)
			{//use TPC GeantFlukaCorrection
				Float_t FlukaCorrKPos=hGeantFlukaKPos->GetBinContent(hGeantFlukaKPos->FindBin(Spectra[1]->GetBinCenter(binK)));
				Float_t FlukaCorrKNeg=hGeantFlukaKNeg->GetBinContent(hGeantFlukaKNeg->FindBin(Spectra[4]->GetBinCenter(binK)));
			//	Printf("TPC Geant/Fluka: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrKPos,FlukaCorrKNeg);
				Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPos);
				Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNeg);
				Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPos);
				Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNeg);
				gGFCorrectionKaonPlus->SetPoint(binK,Spectra[1]->GetBinCenter(binK),FlukaCorrKPos);
				gGFCorrectionKaonMinus->SetPoint(binK,Spectra[4]->GetBinCenter(binK),FlukaCorrKNeg);
			}
			else
			{
				gGFCorrectionKaonPlus->SetPoint(binK,Spectra[1]->GetBinCenter(binK),0);
				gGFCorrectionKaonMinus->SetPoint(binK,Spectra[4]->GetBinCenter(binK),0);
				Float_t FlukaCorrKPosTracking=fGFKPosTracking->Eval(Spectra[1]->GetBinCenter(binK));
				Float_t FlukaCorrKNegTracking=fGFKNegTracking->Eval(Spectra[1]->GetBinCenter(binK));
			//	Printf("TPC/TOF Geant/Fluka Tracking: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrKPosTracking,FlukaCorrKNegTracking);
				Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPosTracking);
				Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNegTracking);
				Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPosTracking);
				Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNegTracking);
				Float_t FlukaCorrKPosMatching=fGFKPosMatching->Eval(Spectra[1]->GetBinCenter(binK));
				Float_t FlukaCorrKNegMatching=fGFKNegMatching->Eval(Spectra[1]->GetBinCenter(binK));
			//	Printf("TPC/TOF Geant/Fluka Matching: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrKPosMatching,FlukaCorrKNegMatching);
				Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPosMatching);
				Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNegMatching);
				Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPosMatching);
				Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNegMatching);
			}
	  	}
	  }
	  if(!(options&kgeantflukaProton))	
		return;
	  //Geant Fluka for P in TPC
	  Printf("\nGF correction for Protons");
	  const Int_t kNCharge=2;
	  Int_t kPos=0;
	  Int_t kNeg=1;
	  TFile* fGFProtons = new TFile ("GFCorrection/correctionForCrossSection.root");
	  TH2D * hCorrFluka[kNCharge];
	  TH2D * hCorrFluka[2];
	  hCorrFluka[kPos] = (TH2D*)fGFProtons->Get("gHistCorrectionForCrossSectionProtons");
	  hCorrFluka[kNeg] = (TH2D*)fGFProtons->Get("gHistCorrectionForCrossSectionAntiProtons");
	  //getting GF func for Kaons with TPCTOF
	  TF1 *fGFpPosTracking;
	  fGFpPosTracking = TrackingEff_geantflukaCorrection(4,kPositive);
	  TF1 *fGFpNegTracking;
	  fGFpNegTracking = TrackingEff_geantflukaCorrection(4,kNegative);
	  TF1 *fGFpPosMatching;
	  fGFpPosMatching = TOFmatchMC_geantflukaCorrection(4,kPositive);
	  TF1 *fGFpNegMatching;
	  fGFpNegMatching = TOFmatchMC_geantflukaCorrection(4,kNegative);
	  
	 
		Int_t nbins = Spectra[2]->GetNbinsX();
		
		for(Int_t ibin = 0; ibin < nbins; ibin++)
		{
			if(Spectra[2]->GetBinCenter(ibin)<tofpt)
			{//use TPC GeantFlukaCorrection
				for(Int_t icharge = 0; icharge < kNCharge; icharge++)
				{
					Int_t nbinsy=hCorrFluka[icharge]->GetNbinsY();
					Float_t pt = Spectra[2]->GetBinCenter(ibin);
					Float_t minPtCorrection = hCorrFluka[icharge]->GetYaxis()->GetBinLowEdge(1);
					Float_t maxPtCorrection = hCorrFluka[icharge]->GetYaxis()->GetBinLowEdge(nbinsy+1);
					if (pt < minPtCorrection) 
						pt = minPtCorrection+0.0001;
					if (pt > maxPtCorrection) 
						pt = maxPtCorrection;
					Float_t correction = hCorrFluka[icharge]->GetBinContent(1,hCorrFluka[icharge]->GetYaxis()->FindBin(pt));
				
					if (correction > 0.0) 
					{// If the bin is empty this is a  0
			//			cout<<icharge<<" "<<ibin<<" "<<correction<<endl;
						Spectra[icharge*3+2]->SetBinContent(ibin,Spectra[icharge*3+2]->GetBinContent(ibin)*correction);
						Spectra[icharge*3+2]->SetBinError(ibin,Spectra[icharge*3+2]->GetBinError(ibin)*correction);
					}
					else if (Spectra[icharge*3+2]->GetBinContent(ibin) > 0.0) 
					{ 
						// If we are skipping a non-empty bin, we notify the user
						cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for " <<"protons Pt:"<< pt<< endl;
						cout << " Bin content: " << Spectra[icharge*3+2]->GetBinContent(ibin)  << endl;
					}
				}
			}
			else
			{
				Float_t FlukaCorrpPosTracking=fGFpPosTracking->Eval(Spectra[2]->GetBinCenter(ibin));
				Float_t FlukaCorrpNegTracking=fGFpNegTracking->Eval(Spectra[2]->GetBinCenter(ibin));
			//	Printf("TPC/TOF Geant/Fluka Tracking: pt:%f  Pos:%f  Neg:%f",Spectra[2]->GetBinCenter(ibin),FlukaCorrpPosTracking,FlukaCorrpNegTracking);
				Spectra[2]->SetBinContent(ibin,Spectra[2]->GetBinContent(ibin)*FlukaCorrpPosTracking);
				Spectra[5]->SetBinContent(ibin,Spectra[5]->GetBinContent(ibin)*FlukaCorrpNegTracking);
				Spectra[2]->SetBinError(ibin,Spectra[2]->GetBinError(ibin)*FlukaCorrpPosTracking);
				Spectra[5]->SetBinError(ibin,Spectra[5]->GetBinError(ibin)*FlukaCorrpNegTracking);
				Float_t FlukaCorrpPosMatching=fGFpPosMatching->Eval(Spectra[2]->GetBinCenter(ibin));
				Float_t FlukaCorrpNegMatching=fGFpNegMatching->Eval(Spectra[2]->GetBinCenter(ibin));
			//	Printf("TPC/TOF Geant/Fluka Matching: pt:%f  Pos:%f  Neg:%f",Spectra[2]->GetBinCenter(ibin),FlukaCorrpPosMatching,FlukaCorrpNegMatching);
				Spectra[2]->SetBinContent(ibin,Spectra[2]->GetBinContent(ibin)*FlukaCorrpPosMatching);
				Spectra[5]->SetBinContent(ibin,Spectra[5]->GetBinContent(ibin)*FlukaCorrpNegMatching);
				Spectra[2]->SetBinError(ibin,Spectra[2]->GetBinError(ibin)*FlukaCorrpPosMatching);
				Spectra[5]->SetBinError(ibin,Spectra[5]->GetBinError(ibin)*FlukaCorrpNegMatching);
			}		
		}
	 fGeanFlukaK->Close();
	 delete fGeanFlukaK;
}


///////////
TF1 *
TrackingEff_geantflukaCorrection(Int_t ipart, Int_t icharge)
{

  if (ipart == 3 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "TrackingPtGeantFlukaCorrectionKaMinus(x)", 0., 5.);
    return f;
  }
  else if (ipart == 4 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "TrackingPtGeantFlukaCorrectionPrMinus(x)", 0., 5.);
  }
  else
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "TrackingPtGeantFlukaCorrectionNull(x)", 0., 5.);

  return f;
}

Double_t
TrackingPtGeantFlukaCorrectionNull(Double_t pTmc)
{
  return 1.;
}

Double_t
TrackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
  return (1 - 0.129758 *TMath::Exp(-pTmc*0.679612));
}

Double_t
TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
  return TMath::Min((0.972865 + 0.0117093*pTmc), 1.);
}
///////////////////////////////////////////
TF1 *
TOFmatchMC_geantflukaCorrection(Int_t ipart, Int_t icharge)
{

  if (ipart == 3 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "MatchingPtGeantFlukaCorrectionKaMinus(x)", 0., 5.);
    return f;
  }
  else if (ipart == 4 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "MatchingPtGeantFlukaCorrectionPrMinus(x)", 0., 5.);
  }
  else
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "MatchingPtGeantFlukaCorrectionNull(x)", 0., 5.);

  return f;
}


Double_t
MatchingPtGeantFlukaCorrectionNull(Double_t pTmc)
{
  return 1.;
}

Double_t 
MatchingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
  Float_t ptTPCoutP =pTmc*(1-6.81059e-01*TMath::Exp(-pTmc*4.20094));
  return (TMath::Power(1 - 0.129758*TMath::Exp(-ptTPCoutP*0.679612),0.07162/0.03471));
}

Double_t
MatchingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
  Float_t ptTPCoutK=pTmc*(1- 3.37297e-03/pTmc/pTmc - 3.26544e-03/pTmc);
  return TMath::Min((TMath::Power(0.972865 + 0.0117093*ptTPCoutK,0.07162/0.03471)), 1.);
}
void MatchingTOFEff(TH1F** Spectra, TList* list=0x0)
{
	  if(TOFMatchingScalling[0]<0.0&&TOFMatchingScalling[1]<0.0)
	  {
		  TH1F *hMatcEffPos_data=(TH1F*)tcutsdata->GetHistoNMatchedPos();
		  hMatcEffPos_data->Divide((TH1F*)tcutsdata->GetHistoNSelectedPos());
		  hMatcEffPos_data->SetTitle("Matching Eff Pos - data");
		  TH1F *hMatcEffNeg_data=(TH1F*)tcutsdata->GetHistoNMatchedNeg();
		  hMatcEffNeg_data->Divide((TH1F*)tcutsdata->GetHistoNSelectedNeg());
		  hMatcEffNeg_data->SetTitle("Matching Eff Neg - data");
		  TH1F *hMatcEffPos_mc=(TH1F*)tcutsmc->GetHistoNMatchedPos();
		  hMatcEffPos_mc->Divide((TH1F*)tcutsmc->GetHistoNSelectedPos());
		  hMatcEffPos_mc->SetTitle("Matching Eff Pos - mc");
		  TH1F *hMatcEffNeg_mc=(TH1F*)tcutsmc->GetHistoNMatchedNeg();
		  hMatcEffNeg_mc->Divide((TH1F*)tcutsmc->GetHistoNSelectedNeg());
		  hMatcEffNeg_mc->SetTitle("Matching Eff Neg - mc");


		  hMatcEffPos_data->Divide(hMatcEffPos_mc);
		  hMatcEffNeg_data->Divide(hMatcEffNeg_mc);
		  hMatcEffPos_data->SetName("MatchingTOFPos");
		  hMatcEffNeg_data->SetName("MatchingTOFNeg");
		  
		  
		  TF1 *pol0MatchPos_data=new TF1("pol0MatchPos_data","pol0",.6,5);
		  hMatcEffPos_data->Fit("pol0MatchPos_data","MNR");
		  TF1 *pol0MatchNeg_data=new TF1("pol0MatchNeg_data","pol0",.6,5);
		  hMatcEffNeg_data->Fit("pol0MatchNeg_data","MNR");
		
		list->Add(hMatcEffPos_data);
		  list->Add(hMatcEffNeg_data);
		
			
		  TOFMatchingScalling[0]=pol0MatchPos_data->GetParameter(0);
		  TOFMatchingScalling[1]=pol0MatchNeg_data->GetParameter(0);
	  }
	  //Correction spectra for matching efficiency
	  //For the moment I'm using the inclusive correction
	  for(Int_t ipart=0;ipart<3;ipart++)
	  {
		  
		for(Int_t ibin=1;ibin<Spectra[ipart]->GetNbinsX();ibin++)
		{
			Float_t ptspectra=Spectra[ipart]->GetBinCenter(ibin);
			if(ptspectra<tcutsdata->GetPtTOFMatching())
				continue;
		  //Spectra[ipart]->SetBinContent(ibin,( Spectra[ipart]->GetBinContent(ibin)/hMatcEffPos_data->GetBinContent(hMatcEffPos_data->FindBin(ptspectra))));
		  //Spectra[ipart+3]->SetBinContent(ibin,( Spectra[ipart+3]->GetBinContent(ibin)/hMatcEffNeg_data->GetBinContent(hMatcEffNeg_data->FindBin(ptspectra))));
			Spectra[ipart]->SetBinContent(ibin,( Spectra[ipart]->GetBinContent(ibin)/TOFMatchingScalling[0]));
			Spectra[ipart+3]->SetBinContent(ibin,( Spectra[ipart+3]->GetBinContent(ibin)/TOFMatchingScalling[1]));
		}
	  }
}
Double_t eta2y(Double_t pt, Double_t mass, Double_t eta)
{
  Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
  return TMath::ASinH(pt / mt * TMath::SinH(eta));
}

TH1* GetSumAllCh(TH1F** spectra, Double_t* mass  )
{
	TH1F* allch=(((TH1F*))spectra[0]->Clone("allCh"));
	allch->Reset();
	for (int i=0;i<6;i++)
	{
		Double_t masstmp=mass[i%3];
		for (int j=1;j<spectra[i]->GetXaxis()->GetNbins();j++)
		{
			Float_t value=spectra[i]->GetBinContent(j);
			Float_t error=spectra[i]->GetBinError(j);
			if(value>0.0)
			{
				Float_t pt=spectra[i]->GetXaxis()->GetBinCenter(j);
				Float_t mt2=masstmp*masstmp+pt*pt;		
				Float_t mt=TMath::Sqrt(mt2);
				Float_t maxy=eta2y(pt,masstmp,0.8);
				Float_t conver=maxy*(TMath::Sqrt(1-masstmp*masstmp/(mt2*TMath::CosH(maxy)*TMath::CosH(maxy)))+TMath::Sqrt(1-masstmp*masstmp/(mt2*TMath::CosH(0.0)*TMath::CosH(0.0))));
				conver=conver/1.6;
				cout<<maxy<<" "<<conver<<" "<<masstmp<<""<<spectra[i]->GetName()<<endl;
				Float_t bincontent=allch->GetBinContent(j);
				Float_t binerror=allch->GetBinError(j);
				bincontent+=conver*value;
				binerror=TMath::Sqrt(binerror*binerror+conver*conver*error*error);
				allch->SetBinContent(j,bincontent);
				allch->SetBinError(j,binerror);
			}
			
		}
	}
	Divideby2pipt(allch);
	allch->SetTitle("N_{ch};p_{T} (GeV/c);1/(2 #pi p_{T})dN/p_{T} |#eta|<0.8");		
	return allch;
}

void Divideby2pipt(TH1* hist)
{

	for (int i=1;i<hist->GetXaxis()->GetNbins();i++)
	{
		Float_t value=hist->GetBinContent(i);
		Float_t error=hist->GetBinError(i);
		Float_t pt=hist->GetXaxis()->GetBinCenter(i);
		hist->SetBinContent(i,value/(2.0*TMath::Pi()*pt));
		hist->SetBinError(i,error/(2.0*TMath::Pi()*pt));

	}
}
