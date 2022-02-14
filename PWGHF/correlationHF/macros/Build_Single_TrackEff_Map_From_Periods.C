void Build_Single_3DTrackEff_Map_From_Periods(Bool_t isHMmap = kTRUE /* choose if HM or MB */, TString inputfolder="./", TString filemapprefix="3D_TrackingEffMap") {

	// *******************************************************************************	
	// Shall be executed on the particle-reweighted maps, so after the macro
	// Reweight_3DTrackEff_By_Species.C
	// *******************************************************************************	
	// Tuned for 13 TeV, so HM sums all periods from 16k included, MB sums everything
	// *******************************************************************************

	TString periods[33] = {"16d","16e","16g","16h","16j","16k","16l","16o","16p","17c","17e","17f","17h","17i","17j","17k","17l","17m","17o","17r","18b","18d","18e","18f","18g","18h","18i","18k","18l","18m","18n","18o","18p"};
    Double_t events_HM[33] = { //the difference is that up to 16k we put artificially 0!! Skipping the bad periods
    	0, //16d - MANUAL 0
    	0, //16e - MANUAL 0
    	0, //16g - MANUAL 0
    	0, //16h - MANUAL 0
    	0, //16j - MANUAL 0
    	103377100, //16k
    	38980848, //16l
    	28561273, //16o
    	46008822, //16p
    	0, //17c
    	378907, //17e
    	170513, //17f
    	32733157, //17h
    	22228350, //17i
    	0, //17j
    	60490144, //17k
    	76619050, //17l
    	77119780, //17m
    	101078476, //17o
    	23369262, //17r
    	674821, //18b
    	26299842, //18d
    	30454814, //18e
    	34948496, //18f
    	0, //18g
    	2437652, //18h
    	0, //18i
    	4560894, //18k
    	39255375, //18l
    	108209700, //18m
    	46680508, //18n
    	21538523, //18o
    	0 //18p
    };
    Double_t events_MB[33] = { //all the periods are ok
    	15349249, //16d
    	49958203, //16e
    	26503561, //16g
    	69350174, //16h
    	44565728, //16j
    	116282427, //16k
    	28930354, //16l
    	32532599, //16o
    	20205783, //16p
    	9050465, //17c
    	9875989, //17e
    	9172750, //17f
    	117608605, //17h
    	43531035, //17i
    	39897714, //17j
    	90959876, //17k
    	68415835, //17l
    	96012749, //17m
    	95702356, //17o
    	25067722, //17r
    	171079652, //18b
    	35927318, //18d
    	46810077, //18e
    	51218125, //18f
    	8056857, //18g
    	3464784, //18h
    	50839588, //18i
    	8475072, //18k
    	57888862, //18l
    	176584688, //18m
    	3324104, //18n
    	26951519, //18o
    	60194743 //18p
    };
    Double_t sum_MB = 0;
    Double_t sum_HM = 0;
    Double_t weight_HM[33];
	Double_t weight_MB[33];

	//define weights for HM and MB
	for(int i=0; i<33; i++) {
		sum_HM += events_HM[i];
		sum_MB += events_MB[i];
	}
	for(int i=0; i<33; i++) {	
		weight_HM[i] = events_HM[i]/sum_HM;
		weight_MB[i] = events_MB[i]/sum_MB;
		if(!i) printf("Total number of events: HM = %d, MB = %d\n",(int)sum_HM,(int)sum_MB);
		printf("Weight of period %s: HM = %.3f (%.2e ev)\t MB = %.3f (%.2e ev)\n",periods[i].Data(),weight_HM[i],(int)events_HM[i],weight_MB[i],(int)events_MB[i]);
	}

	//define and load input maps
	TFile *fMap[33];
	TCanvas *cMap[33];
	TH3D *hMap[33];
	TH3D *hMapRebin[33];	

	for(int i=0; i<33; i++) {
		fMap[i] = TFile::Open(Form("%s%s_%s.root",inputfolder.Data(),filemapprefix.Data(),periods[i].Data()));
		cMap[i] = (TCanvas*)fMap[i]->Get("c");
		cMap[i]->SetName(Form("c_%d",i));
		hMap[i] = (TH3D*)cMap[i]->FindObject("heff");
		hMap[i]->SetName(Form("heff_%d",i));
		hMapRebin[i] = (TH3D*)cMap[i]->FindObject("heff_rebin");
		hMapRebin[i]->SetName(Form("heff_rebin_%d",i));
	}

	//create final map
	TH3D *hMapFinal = (TH3D*)hMap[32]->Clone("heff");
	TH3D *hMapFinalRebin = (TH3D*)hMapRebin[32]->Clone("heff_rebin");
	hMapFinal->Reset();
	hMapFinalRebin->Reset();

	//reweighting and building of final map
	for(int i=0; i<33; i++) {
		if(isHMmap) {
			hMapFinal->Add(hMap[i],weight_HM[i]); //the weight is the second argument of Add!!
			hMapFinalRebin->Add(hMapRebin[i],weight_HM[i]);
		} else {
			hMapFinal->Add(hMap[i],weight_MB[i]);
			hMapFinalRebin->Add(hMapRebin[i],weight_MB[i]);		
		}
	}	

	//create final canvas
	TCanvas *cOut = new TCanvas("c","pT, Eta, Zvtx Efficiency distribution",400,900);
	cOut->Divide(1,2);	
	cOut->cd(1);
	hMapFinal->DrawClone();
	cOut->cd(2);
	hMapFinalRebin->DrawClone();

	//save final map
	if(isHMmap) {
		cOut->SaveAs(Form("%s_All13TeV_ForHM.png",filemapprefix.Data()));
		cOut->SaveAs(Form("%s_All13TeV_ForHM.root",filemapprefix.Data()));		
	} else {
		cOut->SaveAs(Form("%s_All13TeV_ForMB.png",filemapprefix.Data()));
		cOut->SaveAs(Form("%s_All13TeV_ForMB.root",filemapprefix.Data()));				
	}

}

void Build_Single_1DTrackEff_Map_From_Periods(Bool_t isHMmap = kTRUE /* choose if HM or MB */, TString inputfolder="./", TString filemapprefix="1D_TrackingEffMap") {

	// *******************************************************************************	
	// Tuned for 13 TeV, so HM sums all periods from 16k included, MB sums everything
	// *******************************************************************************

	TString periods[33] = {"16d","16e","16g","16h","16j","16k","16l","16o","16p","17c","17e","17f","17h","17i","17j","17k","17l","17m","17o","17r","18b","18d","18e","18f","18g","18h","18i","18k","18l","18m","18n","18o","18p"};
    Double_t events_HM[33] = { //the difference is that up to 16k we put artificially 0!! Skipping the bad periods
    	0, //16d - MANUAL 0
    	0, //16e - MANUAL 0
    	0, //16g - MANUAL 0
    	0, //16h - MANUAL 0
    	0, //16j - MANUAL 0
    	103377100, //16k
    	38980848, //16l
    	28561273, //16o
    	46008822, //16p
    	0, //17c
    	378907, //17e
    	170513, //17f
    	32733157, //17h
    	22228350, //17i
    	0, //17j
    	60490144, //17k
    	76619050, //17l
    	77119780, //17m
    	101078476, //17o
    	23369262, //17r
    	674821, //18b
    	26299842, //18d
    	30454814, //18e
    	34948496, //18f
    	0, //18g
    	2437652, //18h
    	0, //18i
    	4560894, //18k
    	39255375, //18l
    	108209700, //18m
    	46680508, //18n
    	21538523, //18o
    	0 //18p
    };
    Double_t events_MB[33] = { //all the periods are ok
    	15349249, //16d
    	49958203, //16e
    	26503561, //16g
    	69350174, //16h
    	44565728, //16j
    	116282427, //16k
    	28930354, //16l
    	32532599, //16o
    	20205783, //16p
    	9050465, //17c
    	9875989, //17e
    	9172750, //17f
    	117608605, //17h
    	43531035, //17i
    	39897714, //17j
    	90959876, //17k
    	68415835, //17l
    	96012749, //17m
    	95702356, //17o
    	25067722, //17r
    	171079652, //18b
    	35927318, //18d
    	046810077, //18e
    	51218125, //18f
    	8056857, //18g
    	3464784, //18h
    	50839588, //18i
    	8475072, //18k
    	57888862, //18l
    	176584688, //18m
    	3324104, //18n
    	26951519, //18o
    	60194743 //18p
    };
    Double_t sum_MB = 0;
    Double_t sum_HM = 0;
    Double_t weight_HM[33];
	Double_t weight_MB[33];

	//define weights for HM and MB
	for(int i=0; i<33; i++) {
		sum_HM += events_HM[i];
		sum_MB += events_MB[i];
	}
	for(int i=0; i<33; i++) {	
		weight_HM[i] = events_HM[i]/sum_HM;
		weight_MB[i] = events_MB[i]/sum_MB;
		if(!i) printf("Total number of events: HM = %d, MB = %d\n",(int)sum_HM,(int)sum_MB);
		printf("Weight of period %s: HM = %.3f (%.2e ev)\t MB = %.3f (%.2e ev)\n",periods[i].Data(),weight_HM[i],(int)events_HM[i],weight_MB[i],(int)events_MB[i]);
	}

	//define and load input maps
	TFile *fMap[33];
	TCanvas *cMap[33];
	TH1D *hMap[33];
	TH1D *hMapRebin[33];	

	for(int i=0; i<33; i++) {
		fMap[i] = TFile::Open(Form("%s%s_%s.root",inputfolder.Data(),filemapprefix.Data(),periods[i].Data()));
		cMap[i] = (TCanvas*)fMap[i]->Get("c");
		cMap[i]->SetName(Form("c_%d",i));
		hMap[i] = (TH1D*)cMap[i]->FindObject("heff");
		hMap[i]->SetName(Form("heff_%d",i));
		hMapRebin[i] = (TH1D*)cMap[i]->FindObject("heff_rebin");
		hMapRebin[i]->SetName(Form("heff_rebin_%d",i));
	}

	//create final map
	TH1D *hMapFinal = (TH1D*)hMap[32]->Clone("heff");
	TH1D *hMapFinalRebin = (TH1D*)hMapRebin[32]->Clone("heff_rebin");
	hMapFinal->Reset();
	hMapFinalRebin->Reset();

	//reweighting and building of final map
	for(int i=0; i<33; i++) {
		if(isHMmap) {
			hMapFinal->Add(hMap[i],weight_HM[i]); //the weight is the second argument of Add!!
			hMapFinalRebin->Add(hMapRebin[i],weight_HM[i]);
		} else {
			hMapFinal->Add(hMap[i],weight_MB[i]);
			hMapFinalRebin->Add(hMapRebin[i],weight_MB[i]);		
		}
	}	

	//create final canvas
	TCanvas *cOut = new TCanvas("c","pT, Eta, Zvtx Efficiency distribution",400,900);
	cOut->Divide(1,2);	
	cOut->cd(1);
	hMapFinal->DrawClone();
	cOut->cd(2);
	hMapFinalRebin->DrawClone();

	//save final map
	if(isHMmap) {
		cOut->SaveAs(Form("%s_All13TeV_ForHM.png",filemapprefix.Data()));
		cOut->SaveAs(Form("%s_All13TeV_ForHM.root",filemapprefix.Data()));		
	} else {
		cOut->SaveAs(Form("%s_All13TeV_ForMB.png",filemapprefix.Data()));
		cOut->SaveAs(Form("%s_All13TeV_ForMB.root",filemapprefix.Data()));				
	}

}

void Ratio_All13TeV_Over_SinglePeriod_1D(TString inputfileFin = "1D_TrackingEffMap_Reweighted_All13TeV_ForMB", TString inputfileSin = "1D_TrackingEffMap_Reweighted", TString period="18m", Int_t color=1) {

	TFile *fFin = new TFile(Form("%s.root",inputfileFin.Data()),"read");
	TFile *fSin = new TFile(Form("%s_%s.root",inputfileSin.Data(),period.Data()),"read");

	TCanvas *cFin = (TCanvas*)fFin->Get("c");
	TCanvas *cSin = (TCanvas*)fSin->Get("c");
		
	TH1D* hFin = (TH1D*)cFin->FindObject("heff_rebin");
	TH1D* hSin = (TH1D*)cSin->FindObject("heff_rebin");

    TH1D* hRatio = (TH1D*)hFin->Clone("heff_ratio");		
    hRatio->Divide(hSin);

    TCanvas *cOut = new TCanvas("cOut","Ratio final/single period",900,900);
    hRatio->SetLineColor(color);
    hRatio->SetMarkerColor(color);
    hRatio->Draw();

	cOut->SaveAs(Form("RATIO_Final_SinglePeriod_1D_TrackingEffMap_%s.png",period.Data()));
	TFile *fOut = new TFile(Form("RATIO_Final_SinglePeriod_1D_TrackingEffMap_%s.root",period.Data()),"recreate");
	hRatio->Write();
	fOut->Close(); 
}

void Ratio_All13TeV_Over_SinglePeriod_3D(TString inputfileFin = "3D_TrackingEffMap_Reweighted_All13TeV_ForMB", TString inputfileSin = "3D_TrackingEffMap_Reweighted", TString period="18m", Int_t color=1) {

	TFile *fFin = new TFile(Form("%s.root",inputfileFin.Data()),"read");
	TFile *fSin = new TFile(Form("%s_%s.root",inputfileSin.Data(),period.Data()),"read");

	TCanvas *cFin = (TCanvas*)fFin->Get("c");
	TCanvas *cSin = (TCanvas*)fSin->Get("c");
		
	TH3D* hFin = (TH3D*)cFin->FindObject("heff_rebin");
	TH3D* hSin = (TH3D*)cSin->FindObject("heff_rebin");

    TH3D* hRatio = (TH3D*)hFin->Clone("heff_ratio");		
    hRatio->Divide(hSin);

    //DEBUG
    for(int i=1; i<=hRatio->GetNbinsX(); i++) { //loop on pT bins
		for(int j=1; j<=hRatio->GetNbinsY(); j++) { //loop on eta bins
			for(int k=1; k<=hRatio->GetNbinsZ(); k++) { //loop on zVtx bins
				printf("Ratio for (%d,%d,%d) -> pT %.2f, eta %2f, zVtx %2f = %.3f\n",i,j,k,hRatio->GetXaxis()->GetBinCenter(i),hRatio->GetYaxis()->GetBinCenter(j),hRatio->GetZaxis()->GetBinCenter(k),hRatio->GetBinContent(i,j,k));
			}
		}
	}

    TCanvas *cOut = new TCanvas("cOut","Ratio final/single period",900,900);
    hRatio->SetLineColor(color);
    hRatio->SetMarkerColor(color);
    hRatio->Draw();

	cOut->SaveAs(Form("RATIO_Final_SinglePeriod_3D_TrackingEffMap_%s.png",period.Data()));
	TFile *fOut = new TFile(Form("RATIO_Final_SinglePeriod_3D_TrackingEffMap_%s.root",period.Data()),"recreate");
	hRatio->Write();
	fOut->Close();
}

void Ratio_All13TeV_Over_SinglePeriod_1D_ForAll(TString inputfileFin = "1D_TrackingEffMap_Reweighted_All13TeV_ForMB", TString inputfileSin = "1D_TrackingEffMap_Reweighted") {

    TString file[33] = {"16d","16e","16g","16h","16j","16k","16l","16o","16p","17c","17e","17f","17h","17i","17j","17k","17l","17m","17o","17r","18b","18d","18e","18f","18g","18h","18i","18k","18l","18m","18n","18o","18p"};

	for(int i=0; i<33; i++) {
		Ratio_All13TeV_Over_SinglePeriod_1D(inputfileFin.Data(),inputfileSin.Data(),file[i].Data(),i+1);
	}

}

void Ratio_All13TeV_Over_SinglePeriod_3D_ForAll(TString inputfileFin = "3D_TrackingEffMap_Reweighted_All13TeV_ForMB", TString inputfileSin = "1D_TrackingEffMap_Reweighted") {

    TString file[33] = {"16d","16e","16g","16h","16j","16k","16l","16o","16p","17c","17e","17f","17h","17i","17j","17k","17l","17m","17o","17r","18b","18d","18e","18f","18g","18h","18i","18k","18l","18m","18n","18o","18p"};

	for(int i=0; i<33; i++) {
		Ratio_All13TeV_Over_SinglePeriod_3D(inputfileFin.Data(),inputfileSin.Data(),file[i].Data(),i+1);
	}

}
