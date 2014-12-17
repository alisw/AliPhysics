#include  "AliJEbeHistos.h"
#include  "AliJCard.h"

//______________________________________________________________________________
AliJEbeHistos::AliJEbeHistos(AliJCard* cardP):
	fcard(cardP),
	fmaxEtaRange(0),
	fmaxTriggEtaRange(0),
	ftriggFiducCut(0),
	fhistoList(NULL),
	fUseDirectory(true),
	fTopDirectory(NULL)
{   // constructor

	fcard=cardP;
	fmaxEtaRange = fcard->Get("EtaRange");
	ftriggFiducCut =  fcard->Get("TriggerFiducialEtaCut"); //FK// Fiduc cut 
	fmaxTriggEtaRange =  fmaxEtaRange - ftriggFiducCut; //FK// Trigger range

	fhistoList = new TList();
	fTopDirectory = gDirectory;

	//fhtyp[1] = "Real";
	//fhtyp[2] = "Mixed";
	//fhtyp[3] = "Rap. Gap";
}

//______________________________________________________________________________
AliJEbeHistos::AliJEbeHistos(const AliJEbeHistos& obj):
	fcard(obj.fcard),
	fmaxEtaRange(obj.fmaxEtaRange),
	fmaxTriggEtaRange(obj.fmaxTriggEtaRange),
	ftriggFiducCut(obj.ftriggFiducCut),
	fhistoList(obj.fhistoList),
	fUseDirectory(obj.fUseDirectory),
	fTopDirectory(obj.fTopDirectory)
{
	// copy constructor
}

//______________________________________________________________________________
AliJEbeHistos& AliJEbeHistos::operator=(const AliJEbeHistos& obj)
{

	JUNUSED(obj);
	// copy constructor
	return *this;
}


void AliJEbeHistos::CreateUnfoldingHistos(){
	// pi0mass histos

	int    hic;
	float  b1 = 0, b2 = 0.50;
	int    bins=300;
	double lbin = 0.0;
	double hbin = 0.5;
	const int NPTBins = 51;

	Double_t PTBINS[NPTBins+1] = {
		0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
		0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
		1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
		2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
		4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.5, 14.0,
		21.0, 32.0};

	if(fcard->GetNoOfBins(kCentrType) > kMaxNoCentrBin ){
		cout<<"ERROR: No of Centrality bins exceed max dim in AliJEbeHistos.cxx "<<endl;
		exit(0);
	}

	TDirectory * cwd = gDirectory;
	TDirectory * nwd = gDirectory->mkdir("UnfoldingHistos");
	if( nwd ) nwd->cd();

	for( hic=0;hic<fcard->GetNoOfBins(kCentrType);hic++){
		for(int ih = 1; ih < kNHarmonics; ih++){
			b1 = fcard->GetBinBorder(kCentrType, hic);
			b2 = fcard->GetBinBorder(kCentrType, hic+1);
			//inclusive fpt spectra================================================
			sprintf(fhname,"hVnObsVector%02d%02d", hic, ih);
			sprintf(fhtit, "Observed v_{%d}: %3.0f -%3.0f ", ih, b1, b2);
			fhVnObsVector[hic][ih] = new TH1D(fhname, fhtit, bins, lbin,hbin); fhVnObsVector[hic][ih]->Sumw2();
			fhistoList->Add(fhVnObsVector[hic][ih]);
			sprintf(fhname,"hResponseDist%02d%02d", hic,ih);
			sprintf(fhtit, "Response Distribution v_{%d}: %3.0f -%3.0f ", ih, b1, b2);
			fhResponseDist[hic][ih] = new TH2D(fhname, fhtit, bins, -0.5, 0.5, bins, -0.5,0.5); fhResponseDist[hic][ih]->Sumw2();
			fhistoList->Add(fhResponseDist[hic][ih]);
			sprintf(fhname,"hVnObsEP%02d%02d", hic,ih);
			sprintf(fhtit, "Observed EP v_{%d}: %3.0f -%3.0f ", ih, b1, b2);
			fhVnObsEP[hic][ih] = new TH1D(fhname, fhtit, bins, -0.2,0.5); fhVnObsEP[hic][ih]->Sumw2();
			fhistoList->Add(fhVnObsEP[hic][ih]);
			sprintf(fhname,"hCosndPhiPt%02d%02d", hic,ih);
			sprintf(fhtit, "Cos%d#Delta#phi Pt: %3.0f -%3.0f ", ih, b1, b2);
			fhCosndPhiPt[hic][ih] = new TH1D(fhname, fhtit, NPTBins, PTBINS); fhCosndPhiPt[hic][ih]->Sumw2();
			fhistoList->Add(fhCosndPhiPt[hic][ih]);
			sprintf(fhname,"heCosndPhiPt%02d%02d", hic,ih);
			sprintf(fhtit, "eCos%d#Delta#phi Pt: %3.0f -%3.0f ", ih, b1, b2);
			fheCosndPhiPt[hic][ih] = new TH1D(fhname, fhtit, NPTBins, PTBINS); fheCosndPhiPt[hic][ih]->Sumw2();
			fhistoList->Add(fheCosndPhiPt[hic][ih]);
			sprintf(fhname,"hCounter%02d%02d", hic,ih);
			sprintf(fhtit, "Counter Pt %02d: %3.0f -%3.0f ", ih, b1, b2);
			fhCounter[hic][ih] = new TH1D(fhname, fhtit, NPTBins, PTBINS); fhCounter[hic][ih]->Sumw2();
			fhistoList->Add(fhCounter[hic][ih]);

			sprintf(fhname,"hEventPlane%02d%02d", hic,ih);
			sprintf(fhtit, "Event Planes %02d: %3.0f -%3.0f ", ih, b1, b2);
			fhEventPlane[hic][ih] = new TH1D(fhname, fhtit, 100, -6.5,6.5); fhEventPlane[hic][ih]->Sumw2();
			fhistoList->Add(fhEventPlane[hic][ih]);

			sprintf(fhtit, "Q vector from V0 %02d: %3.0f -%3.0f ", ih, b1, b2);
			fhQvectorV0[hic][ih] = new TH1D(Form("hQvectorV0_%02d%02d",hic,ih),fhtit ,400, 0, 100); fhQvectorV0[hic][ih]->Sumw2();
			sprintf(fhtit, "Q vector from V0A %02d: %3.0f -%3.0f ", ih, b1, b2);
			fhQvectorV0A[hic][ih] = new TH1D(Form("hQvectorV0A_%02d%02d",hic,ih),fhtit ,400, 0, 100); fhQvectorV0A[hic][ih]->Sumw2();
			sprintf(fhtit, "Q vector from V0C %02d: %3.0f -%3.0f ", ih, b1, b2);
			fhQvectorV0C[hic][ih] = new TH1D(Form("hQvectorV0C_%02d%02d",hic,ih), fhtit ,400, 0, 100); fhQvectorV0C[hic][ih]->Sumw2();
			fhistoList->Add(fhQvectorV0[hic][ih]);
			fhistoList->Add(fhQvectorV0A[hic][ih]);
			fhistoList->Add(fhQvectorV0C[hic][ih]);

			sprintf(fhtit, "Q vector correlation V0A and V0C %02d: %3.0f -%3.0f ", ih, b1, b2);
			fhQvectorCorrelation[hic][ih] = new TH2D(Form("hQvectorCorrelation_%02d%02d",hic,ih),fhtit , 400, 0, 100, 400, 0 ,100); fhQvectorCorrelation[hic][ih]->Sumw2();
			fhistoList->Add(fhQvectorCorrelation[hic][ih]);

			sprintf(fhtit, "VnObs vs V0C Qvector correlation %02d: %3.0f -%3.0f ", ih, b1, b2);
			fhVnObsVsQvectorCorrelation[hic][ih] = new TH2D(Form("hVnObsVsQvectorCorrelation_%02d%02d",hic,ih),fhtit , 400, 0, 0.5, 400, 0 ,100); fhVnObsVsQvectorCorrelation[hic][ih]->Sumw2();
			fhistoList->Add(fhVnObsVsQvectorCorrelation[hic][ih]);
		}
		sprintf(fhname,"hMultiCount%02d", hic);
		sprintf(fhtit, "Multiplicity count: %3.0f -%3.0f ", b1, b2);
		fhMultiCount[hic] = new TH1D( fhname, fhtit, 1000, 0, 10000 );
		fhistoList->Add(fhMultiCount[hic]);
	}
	// tmp histogram for e-b-e
	for(int ih = 1; ih < kNHarmonics; ih++){
		fhEPCosndPhi[ih] = new TH1D(Form("hEPCosndPhi_%02d",ih), Form("cos%d#Delta#phi",ih), 400, -1.1, 1.1); fhEPCosndPhi[ih]->Sumw2();
		fhEPCosndPhi2[ih] = new TH1D(Form("hEPCosndPhi2_%02d",ih), Form("cos%d#Delta#phi",ih), 400, -1.1, 1.1); fhEPCosndPhi2[ih]->Sumw2();
		fhistoList->Add(fhEPCosndPhi[ih]);
		fhistoList->Add(fhEPCosndPhi2[ih]);


	}
	cwd->cd();
}

