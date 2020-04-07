void Reweight_3DTrackEff_By_Species(TString inputfolder="./", TString filemapprefix="3D_TrackingEffMap", TString filemapsuffix="18m", TString fileabund="pp_Monach13.root") {

	TString species[6] = {"pi","K","p","SigmaM","SigmaP","Rest"};

    //load input abundancies
	TFile *fAbund  = new TFile(fileabund.Data(),"read");
    TH1D *hAbund[6];
    hAbund[0] = (TH1D*)fAbund->Get("RelativeAbundancesData_kPion");    
    hAbund[1] = (TH1D*)fAbund->Get("RelativeAbundancesData_kKaon");
    hAbund[2] = (TH1D*)fAbund->Get("RelativeAbundancesData_kProton");
    hAbund[3] = (TH1D*)fAbund->Get("RelativeAbundancesData_kSigmaMinus");
    hAbund[4] = (TH1D*)fAbund->Get("RelativeAbundancesData_kSigmaPlus");
    hAbund[5] = (TH1D*)hAbund[0]->Clone("RelativeAbundancesData_kRest");
    //build Rest abundancies
    hAbund[5]->Reset();
    for(int i=1; i<=hAbund[5]->GetNbinsX(); i++) {
    	hAbund[5]->SetBinContent(i,1. - hAbund[0]->GetBinContent(i) - hAbund[1]->GetBinContent(i) - hAbund[2]->GetBinContent(i) - hAbund[3]->GetBinContent(i) - hAbund[4]->GetBinContent(i));
    	hAbund[5]->SetBinError(i,0.);
    	printf("DEBUG (LOAD STAGE): Bin %d), pT %.2f, abs (pi,K,p,Sm,Sp,Rest): %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",i,hAbund[0]->GetBinCenter(i),hAbund[0]->GetBinContent(i),hAbund[1]->GetBinContent(i),hAbund[2]->GetBinContent(i),hAbund[3]->GetBinContent(i),hAbund[4]->GetBinContent(i),hAbund[5]->GetBinContent(i));
    }

    //load input maps
	TFile *fMap[6];
	TCanvas *cMap[6];
	TH3D *hMap[6];
	TH3D *hMapRebin[6];	
	for(int i=0; i<6; i++) { //loop on species
		fMap[i] = new TFile(Form("%s%s_%s_%s.root",inputfolder.Data(),filemapprefix.Data(),species[i].Data(),filemapsuffix.Data()),"read");
		cMap[i] = (TCanvas*)fMap[i]->Get("c");
		cMap[i]->SetName(Form("c_%d",i));
		hMap[i] = (TH3D*)cMap[i]->FindObject("heff");
		hMap[i]->SetName(Form("heff_%d",i));
		hMapRebin[i] = (TH3D*)cMap[i]->FindObject("heff_rebin");
		hMapRebin[i]->SetName(Form("heff_rebin_%d",i));	
	}
	//create weightemap
	TH3D *hWeighMap = (TH3D*)hMap[0]->Clone("heff");
	TH3D *hWeighMapRebin = (TH3D*)hMapRebin[0]->Clone("heff_rebin");
	hWeighMap->Reset();
	hWeighMapRebin->Reset();

	//reweight the 3D maps (normal)
	for(int i=1; i<=hMap[0]->GetNbinsX(); i++) { //loop on pT bins
		//evaluate abundancies and for each bin
		Double_t abund[6];
		Double_t pT = hWeighMap->GetXaxis()->GetBinCenter(i);
		for(int s=0; s<6; s++) abund[s] = hAbund[s]->GetBinContent(hAbund[s]->FindBin(pT));
		printf("DEBUG (EVAL STAGE): Bin %d), pT %.2f, abs (pi,K,p,Sm,Sp,Rest): %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",i,pT,abund[0],abund[1],abund[2],abund[3],abund[4],abund[5]);

		for(int j=1; j<=hMap[0]->GetNbinsY(); j++) { //loop on eta bins
			for(int k=1; k<=hMap[0]->GetNbinsZ(); k++) { //loop on zVtx bins
				//evaluate abundancies and for each bin
				Double_t weightedEff = 0., weightedEffErr = 0.;
				for(int s=0; s<6; s++) weightedEff += hMap[s]->GetBinContent(i,j,k)*abund[s];
				for(int s=0; s<6; s++) weightedEffErr += hMap[s]->GetBinError(i,j,k)*hMap[s]->GetBinError(i,j,k)*abund[s]*abund[s];	
				weightedEffErr = TMath::Sqrt(weightedEffErr);
			    //set values in the map
				hWeighMap->SetBinContent(i,j,k,weightedEff);
				hWeighMap->SetBinError(i,j,k,weightedEffErr);
				//DEBUG LINE
				printf("--> WEIGHTING BIN %d,%d,%d): ALL %.3f - pi %.3f (w %.3f) - K %.3f (w %.3f) - p %.3f (w %.3f) - SM %.3f (w %.3f) - SP %.3f (w %.3f) - R %.3f (w %.3f)\n",
			    	i,j,k,weightedEff,hMap[0]->GetBinContent(i,j,k),abund[0],hMap[1]->GetBinContent(i,j,k),abund[1],hMap[2]->GetBinContent(i,j,k),abund[2],hMap[3]->GetBinContent(i,j,k),abund[3],hMap[4]->GetBinContent(i,j,k),abund[4],hMap[5]->GetBinContent(i,j,k),abund[5]);
			}
		}
	}

	//reweight the 3D maps (rebinned)
	for(int i=1; i<=hMapRebin[0]->GetNbinsX(); i++) { //loop on pT bins
		//evaluate abundancies and for each bin
		Double_t abund[6];
		Double_t pT = hWeighMapRebin->GetXaxis()->GetBinCenter(i);
		for(int s=0; s<6; s++) abund[s] = hAbund[s]->GetBinContent(hAbund[s]->FindBin(pT));
		printf("DEBUG (EVAL STAGE): Bin %d), pT %.2f, abs (pi,K,p,Sm,Sp,Rest): %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",i,pT,abund[0],abund[1],abund[2],abund[3],abund[4],abund[5]);

		for(int j=1; j<=hMapRebin[0]->GetNbinsY(); j++) { //loop on eta bins
			for(int k=1; k<=hMapRebin[0]->GetNbinsZ(); k++) { //loop on zVtx bins
				//evaluate abundancies and for each bin
				Double_t weightedEff = 0., weightedEffErr = 0.;
				for(int s=0; s<6; s++) weightedEff += hMapRebin[s]->GetBinContent(i,j,k)*abund[s];
				for(int s=0; s<6; s++) weightedEffErr += hMapRebin[s]->GetBinError(i,j,k)*hMapRebin[s]->GetBinError(i,j,k)*abund[s]*abund[s];	
				weightedEffErr = TMath::Sqrt(weightedEffErr);
			    //set values in the map
				hWeighMapRebin->SetBinContent(i,j,k,weightedEff);
				hWeighMapRebin->SetBinError(i,j,k,weightedEffErr);
			    //DEBUG LINE
			    printf("--> WEIGHT REB BIN %d,%d,%d): ALL %.3f - pi %.3f (w %.3f) - K %.3f (w %.3f) - p %.3f (w %.3f) - SM %.3f (w %.3f) - SP %.3f (w %.3f) - R %.3f (w %.3f)\n",
			    	i,j,k,weightedEff,hMapRebin[0]->GetBinContent(i,j,k),abund[0],hMapRebin[1]->GetBinContent(i,j,k),abund[1],hMapRebin[2]->GetBinContent(i,j,k),abund[2],hMapRebin[3]->GetBinContent(i,j,k),abund[3],hMapRebin[4]->GetBinContent(i,j,k),abund[4],hMapRebin[5]->GetBinContent(i,j,k),abund[5]);
			}
		}
	}

	TCanvas *cOut = new TCanvas("c","pT, Eta, Zvtx Efficiency distrubution",400,900);
	cOut->Divide(1,2);
	cOut->cd(1);
	hWeighMap->DrawClone();
	cOut->cd(2);
	hWeighMapRebin->DrawClone();

	cOut->SaveAs(Form("%s_Reweighted_%s.png",filemapprefix.Data(),filemapsuffix.Data()));
	cOut->SaveAs(Form("%s_Reweighted_%s.root",filemapprefix.Data(),filemapsuffix.Data()));
}


void Reweight_1DTrackEff_By_Species(TString inputfolder="./", TString filemapprefix="1D_TrackingEffMap", TString filemapsuffix="18m", TString fileabund="pp_Monach13.root") {

	TString species[6] = {"pi","K","p","SigmaM","SigmaP","Rest"};

    //load input abundancies
	TFile *fAbund  = new TFile(fileabund.Data(),"read");
    TH1D *hAbund[6];
    hAbund[0] = (TH1D*)fAbund->Get("RelativeAbundancesData_kPion");    
    hAbund[1] = (TH1D*)fAbund->Get("RelativeAbundancesData_kKaon");
    hAbund[2] = (TH1D*)fAbund->Get("RelativeAbundancesData_kProton");
    hAbund[3] = (TH1D*)fAbund->Get("RelativeAbundancesData_kSigmaMinus");
    hAbund[4] = (TH1D*)fAbund->Get("RelativeAbundancesData_kSigmaPlus");
    hAbund[5] = (TH1D*)hAbund[0]->Clone("RelativeAbundancesData_kRest");
    //build Rest abundancies
    hAbund[5]->Reset();
    for(int i=1; i<=hAbund[5]->GetNbinsX(); i++) {
    	hAbund[5]->SetBinContent(i,1. - hAbund[0]->GetBinContent(i) - hAbund[1]->GetBinContent(i) - hAbund[2]->GetBinContent(i) - hAbund[3]->GetBinContent(i) - hAbund[4]->GetBinContent(i));
    	hAbund[5]->SetBinError(i,0.);
    	printf("DEBUG (LOAD STAGE): Bin %d), pT %.2f, abs (pi,K,p,Sm,Sp,Rest): %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",i,hAbund[0]->GetBinCenter(i),hAbund[0]->GetBinContent(i),hAbund[1]->GetBinContent(i),hAbund[2]->GetBinContent(i),hAbund[3]->GetBinContent(i),hAbund[4]->GetBinContent(i),hAbund[5]->GetBinContent(i));
    }

    //load input maps
	TFile *fMap[6];
	TCanvas *cMap[6];
	TH1D *hMap[6];
	TH1D *hMapRebin[6];	
	for(int i=0; i<6; i++) { //loop on species
		fMap[i] = new TFile(Form("%s%s_%s_%s.root",inputfolder.Data(),filemapprefix.Data(),species[i].Data(),filemapsuffix.Data()),"read");
		cMap[i] = (TCanvas*)fMap[i]->Get("c");
		printf("cMap address %p\n",cMap[i]);
		cMap[i]->SetName(Form("c_%d",i));
		hMap[i] = (TH1D*)cMap[i]->FindObject("heff");
		hMap[i]->SetName(Form("heff_%d",i));
		hMapRebin[i] = (TH1D*)cMap[i]->FindObject("heff_rebin");
		hMapRebin[i]->SetName(Form("heff_rebin_%d",i));
	}
	//create weightemap
	TH1D *hWeighMap = (TH1D*)hMap[0]->Clone("heff");
	TH1D *hWeighMapRebin = (TH1D*)hMapRebin[0]->Clone("heff_rebin");
	hWeighMap->Reset();
	hWeighMapRebin->Reset();

	//reweight the 1D maps (normal)
	for(int i=1; i<=hMap[0]->GetNbinsX(); i++) { //loop on pT bins
		//evaluate abundancies and for each bin
		Double_t abund[6];
		Double_t pT = hWeighMap->GetXaxis()->GetBinCenter(i);
		for(int s=0; s<6; s++) abund[s] = hAbund[s]->GetBinContent(hAbund[s]->FindBin(pT));
		printf("DEBUG (EVAL STAGE): Bin %d), pT %.2f, abs (pi,K,p,Sm,Sp,Rest): %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",i,pT,abund[0],abund[1],abund[2],abund[3],abund[4],abund[5]);

		//evaluate abundancies and for each bin
		Double_t weightedEff = 0., weightedEffErr = 0.;
		for(int s=0; s<6; s++) weightedEff += hMap[s]->GetBinContent(i)*abund[s];
		for(int s=0; s<6; s++) weightedEffErr += hMap[s]->GetBinError(i)*hMap[s]->GetBinError(i)*abund[s]*abund[s];	
		weightedEffErr = TMath::Sqrt(weightedEffErr);
	    //set values in the map
		hWeighMap->SetBinContent(i,weightedEff);
		hWeighMap->SetBinError(i,weightedEffErr);
		//DEBUG LINE
		printf("--> WEIGHTING BIN %d): ALL %.3f - pi %.3f (w %.3f) - K %.3f (w %.3f) - p %.3f (w %.3f) - SM %.3f (w %.3f) - SP %.3f (w %.3f) - R %.3f (w %.3f)\n",
			    	i,weightedEff,hMap[0]->GetBinContent(i),abund[0],hMap[1]->GetBinContent(i),abund[1],hMap[2]->GetBinContent(i),abund[2],hMap[3]->GetBinContent(i),abund[3],hMap[4]->GetBinContent(i),abund[4],hMap[5]->GetBinContent(i),abund[5]);			
	}

	//reweight the 1D maps (rebinned)
	for(int i=1; i<=hMapRebin[0]->GetNbinsX(); i++) { //loop on pT bins
		//evaluate abundancies and for each bin
		Double_t abund[6];
		Double_t pT = hWeighMapRebin->GetXaxis()->GetBinCenter(i);
		for(int s=0; s<6; s++) abund[s] = hAbund[s]->GetBinContent(hAbund[s]->FindBin(pT));
		printf("DEBUG (EVAL STAGE): Bin %d), pT %.2f, abs (pi,K,p,Sm,Sp,Rest): %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",i,pT,abund[0],abund[1],abund[2],abund[3],abund[4],abund[5]);

		//evaluate abundancies and for each bin
		Double_t weightedEff = 0., weightedEffErr = 0.;
		for(int s=0; s<6; s++) weightedEff += hMapRebin[s]->GetBinContent(i)*abund[s];
		for(int s=0; s<6; s++) weightedEffErr += hMapRebin[s]->GetBinError(i)*hMapRebin[s]->GetBinError(i)*abund[s]*abund[s];	
		weightedEffErr = TMath::Sqrt(weightedEffErr);
	    //set values in the map
		hWeighMapRebin->SetBinContent(i,weightedEff);
		hWeighMapRebin->SetBinError(i,weightedEffErr);
		//DEBUG LINE
		printf("--> WEIGHTING BIN %d): ALL %.3f - pi %.3f (w %.3f) - K %.3f (w %.3f) - p %.3f (w %.3f) - SM %.3f (w %.3f) - SP %.3f (w %.3f) - R %.3f (w %.3f)\n",
			    	i,weightedEff,hMapRebin[0]->GetBinContent(i),abund[0],hMapRebin[1]->GetBinContent(i),abund[1],hMapRebin[2]->GetBinContent(i),abund[2],hMapRebin[3]->GetBinContent(i),abund[3],hMapRebin[4]->GetBinContent(i),abund[4],hMapRebin[5]->GetBinContent(i),abund[5]);				
	}

	TCanvas *cOut = new TCanvas("c","pT Efficiency distrubution",400,900);
	cOut->Divide(1,2);	
	cOut->cd(1);
	hWeighMap->DrawClone();
	cOut->cd(2);
	hWeighMapRebin->DrawClone();

	cOut->SaveAs(Form("%s_Reweighted_%s.png",filemapprefix.Data(),filemapsuffix.Data()));
	cOut->SaveAs(Form("%s_Reweighted_%s.root",filemapprefix.Data(),filemapsuffix.Data()));
}

void Ratio_Reweight_Over_All_1D(TString inputfileRew = "1D_TrackingEffMap_Reweighted", TString inputfileAll = "1D_TrackingEffMap_All", TString period="18m", Int_t color=1) {

	TFile *fRew = new TFile(Form("%s_%s.root",inputfileRew.Data(),period.Data()),"read");
	TFile *fAll = new TFile(Form("%s_%s.root",inputfileAll.Data(),period.Data()),"read");

	TCanvas *cRew = (TCanvas*)fRew->Get("c");
	TCanvas *cAll = (TCanvas*)fAll->Get("c");
		
	TH1D* hRew = (TH1D*)cRew->FindObject("heff_rebin");
	TH1D* hAll = (TH1D*)cAll->FindObject("heff_rebin");

    TH1D* hRatio = (TH1D*)hRew->Clone("heff_ratio");		
    hRatio->Divide(hAll);

    TCanvas *cOut = new TCanvas("cOut","Ratio reweighted/all species",900,900);
    hRatio->SetLineColor(color);
    hRatio->SetMarkerColor(color);
    hRatio->Draw();

	cOut->SaveAs(Form("RATIO_RewAll_1D_TrackingEffMap_%s.png",period.Data()));
	TFile *fOut = new TFile(Form("RATIO_RewAll_1D_TrackingEffMap_%s.root",period.Data()),"recreate");
	hRatio->Write();
	fOut->Close();	
}

void Ratio_Reweight_Over_All_3D(TString inputfileRew = "3D_TrackingEffMap_Reweighted", TString inputfileAll = "3D_TrackingEffMap_All", TString period="18m", Int_t color=1) {

	TFile *fRew = new TFile(Form("%s_%s.root",inputfileRew.Data(),period.Data()),"read");
	TFile *fAll = new TFile(Form("%s_%s.root",inputfileAll.Data(),period.Data()),"read");

	TCanvas *cRew = (TCanvas*)fRew->Get("c");
	TCanvas *cAll = (TCanvas*)fAll->Get("c");
		
	TH3D* hRew = (TH3D*)cRew->FindObject("heff_rebin");
	TH3D* hAll = (TH3D*)cAll->FindObject("heff_rebin");

    TH3D* hRatio = (TH3D*)hRew->Clone("heff_ratio");		
    hRatio->Divide(hAll);

    //DEBUG
    for(int i=1; i<=hRatio->GetNbinsX(); i++) { //loop on pT bins
		for(int j=1; j<=hRatio->GetNbinsY(); j++) { //loop on eta bins
			for(int k=1; k<=hRatio->GetNbinsZ(); k++) { //loop on zVtx bins
				printf("Ratio for (%d,%d,%d) -> pT %.2f, eta %2f, zVtx %2f = %.3f\n",i,j,k,hRatio->GetXaxis()->GetBinCenter(i),hRatio->GetYaxis()->GetBinCenter(j),hRatio->GetZaxis()->GetBinCenter(k),hRatio->GetBinContent(i,j,k));
			}
		}
	}

    TCanvas *cOut = new TCanvas("cOut","Ratio reweighted/all species",900,900);
    hRatio->SetLineColor(color);
    hRatio->SetMarkerColor(color);
    hRatio->Draw();

	cOut->SaveAs(Form("RATIO_RewAll_3D_TrackingEffMap_%s.png",period.Data()));
	TFile *fOut = new TFile(Form("RATIO_RewAll_3D_TrackingEffMap_%s.root",period.Data()),"recreate");
	hRatio->Write();
	fOut->Close();	
}

void Ratio_Reweight_Over_All_1D_Full13TeV(TString inputfileRew = "1D_TrackingEffMap_Reweighted", TString inputfileAll = "1D_TrackingEffMap_All") {

    TString file[33] = {"16d","16e","16g","16h","16j","16k","16l","16o","16p","17c","17e","17f","17h","17i","17j","17k","17l","17m","17o","17r","18b","18d","18e","18f","18g","18h","18i","18k","18l","18m","18n","18o","18p"};

	for(int i=0; i<33; i++) {
		Ratio_Reweight_Over_All_1D(inputfileRew.Data(),inputfileAll.Data(),file[i].Data(),i+1);
	}

}

void Ratio_Reweight_Over_All_3D_Full13TeV(TString inputfileRew = "3D_TrackingEffMap_Reweighted", TString inputfileAll = "3D_TrackingEffMap_All") {

    TString file[33] = {"16d","16e","16g","16h","16j","16k","16l","16o","16p","17c","17e","17f","17h","17i","17j","17k","17l","17m","17o","17r","18b","18d","18e","18f","18g","18h","18i","18k","18l","18m","18n","18o","18p"};

	for(int i=0; i<33; i++) {
		Ratio_Reweight_Over_All_3D(inputfileRew.Data(),inputfileAll.Data(),file[i].Data(),i+1);
	}

}