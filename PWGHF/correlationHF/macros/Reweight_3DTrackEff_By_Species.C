void Reweight_3DTrackEff_By_Species(TString inputfolder="./", TString filemapprefix="3D_TrackingEffMap", TString filemapsuffix="18m", TString fileabund_data="pp_Monach13.root", TString fileabund_MC="RelAbundances_5Species.root") {

	TString species[5] = {"pi","K","p","e","mu"};

    //load input abundancies - data (pi, K, p)
	TFile *fAbundData  = new TFile(fileabund_data.Data(),"read");
    TH1D *hAbund[5];
    hAbund[0] = (TH1D*)fAbundData->Get("RelativeAbundancesData_kPion");    
    hAbund[1] = (TH1D*)fAbundData->Get("RelativeAbundancesData_kKaon");
    hAbund[2] = (TH1D*)fAbundData->Get("RelativeAbundancesData_kProton");

    //load input abundancies - MC (e, mu)
	TFile *fAbundMC  = new TFile(fileabund_MC.Data(),"read");
	TCanvas *cMCin = fAbundMC->Get("c1_n2");
    hAbund[3] = (TH1D*)cMCin->FindObject("containerpp13TeV_e_SelStep1_proj_0");    
    hAbund[4] = (TH1D*)cMCin->FindObject("containerpp13TeV_mu_SelStep1_proj_0");

    //APPLY REWEIGHTING OF pi, K, p FROM DATA SO THAT THEIR ABUNCANCY SUM REMAINS THE SAME IN MC (I.E. abound(pi+K+p)_data = abound(pi+K+p)_MC) 
    //THIS BECAUSE WE USE MC ABUNDANCIES FOR e AND mu (NOT MEASURED), SO THE TOTAL SUM OF 5 SPECIES ABUNDANCIES SHALL BE 1
    //HERE OF COURSE abound(pi+K+p)_MC IS 1-abound(e+mu)_MC SINCE IN THE MC ABUNDANCY PLOT ONLY THOSE 5 SPECIES ARE CONSIDERED
	for(int i=1; i<=hAbund[0]->GetNbinsX(); i++) { //loop on pT bins of abundancy plots
		Double_t initialPi = hAbund[0]->GetBinContent(i);
		Double_t initialK = hAbund[1]->GetBinContent(i);
		Double_t initialP = hAbund[2]->GetBinContent(i);
		Double_t MCe = hAbund[3]->GetBinContent(hAbund[3]->FindBin(hAbund[0]->GetBinCenter(i))); //hAbund 3,4 and 0,1,2, have different binnings!!
		Double_t MCmu = hAbund[4]->GetBinContent(hAbund[4]->FindBin(hAbund[0]->GetBinCenter(i)));
		Double_t correctedPi = initialPi/(initialPi+initialK+initialP)*(1-MCe-MCmu); //equivalent to initialPi/(initialPi+initialK+initialP)*(MCpi+MCK+MCp)
		Double_t correctedK  = initialK/(initialPi+initialK+initialP)*(1-MCe-MCmu);
		Double_t correctedP  = initialP/(initialPi+initialK+initialP)*(1-MCe-MCmu);
		hAbund[0]->SetBinContent(i,correctedPi);
		hAbund[1]->SetBinContent(i,correctedK);
		hAbund[2]->SetBinContent(i,correctedP);
		printf("DEBUG (CHANGE WEIGHT pi): Bin %d), pT %.2f, initial p,K,pi,MCe,MCmu: %.3f, %.3f, %.3f, %.3f, %.3f, final pi %.3f \n",i,hAbund[0]->GetBinCenter(i),initialPi,initialK,initialP,MCe,MCmu,hAbund[0]->GetBinContent(i));
		printf("DEBUG (CHANGE WEIGHT K) : Bin %d), pT %.2f, initial p,K,pi,MCe,MCmu: %.3f, %.3f, %.3f, %.3f, %.3f, final K  %.3f \n",i,hAbund[0]->GetBinCenter(i),initialPi,initialK,initialP,MCe,MCmu,hAbund[1]->GetBinContent(i));
		printf("DEBUG (CHANGE WEIGHT p) : Bin %d), pT %.2f, initial p,K,pi,MCe,MCmu: %.3f, %.3f, %.3f, %.3f, %.3f, final p  %.3f \n",i,hAbund[0]->GetBinCenter(i),initialPi,initialK,initialP,MCe,MCmu,hAbund[2]->GetBinContent(i));
	}
	
    //load input maps
	TFile *fMap[5];
	TCanvas *cMap[5];
	TH3D *hMap[5];
	TH3D *hMapRebin[5];	
	for(int i=0; i<5; i++) { //loop on species
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
		Double_t abund[5];
		Double_t pT = hWeighMap->GetXaxis()->GetBinCenter(i);
		for(int s=0; s<5; s++) abund[s] = hAbund[s]->GetBinContent(hAbund[s]->FindBin(pT));
		printf("DEBUG (EVAL STAGE): Bin %d), pT %.2f, abs (pi,K,p,e,mu): %.3f, %.3f, %.3f, %.3f, %.3f\n",i,pT,abund[0],abund[1],abund[2],abund[3],abund[4]);

		for(int j=1; j<=hMap[0]->GetNbinsY(); j++) { //loop on eta bins
			for(int k=1; k<=hMap[0]->GetNbinsZ(); k++) { //loop on zVtx bins
				//evaluate abundancies and for each bin
				Double_t weightedEff = 0., weightedEffErr = 0.;
				for(int s=0; s<5; s++) weightedEff += hMap[s]->GetBinContent(i,j,k)*abund[s];
				for(int s=0; s<5; s++) weightedEffErr += hMap[s]->GetBinError(i,j,k)*hMap[s]->GetBinError(i,j,k)*abund[s]*abund[s];	
				weightedEffErr = TMath::Sqrt(weightedEffErr);
			    //set values in the map
				hWeighMap->SetBinContent(i,j,k,weightedEff);
				hWeighMap->SetBinError(i,j,k,weightedEffErr);
				//DEBUG LINE
				if(j<3 && k<3) printf("--> WEIGHTING BIN %d,%d,%d (pT %.3f)): ALL %.3f - pi %.3f (w %.3f) - K %.3f (w %.3f) - p %.3f (w %.3f) - e %.3f (w %.3f) - mu %.3f (w %.3f)\n",
			    	i,j,k,pT,weightedEff,hMap[0]->GetBinContent(i,j,k),abund[0],hMap[1]->GetBinContent(i,j,k),abund[1],hMap[2]->GetBinContent(i,j,k),abund[2],hMap[3]->GetBinContent(i,j,k),abund[3],hMap[4]->GetBinContent(i,j,k),abund[4]);
			}
		}
	}

	//reweight the 3D maps (rebinned)
	for(int i=1; i<=hMapRebin[0]->GetNbinsX(); i++) { //loop on pT bins
		//evaluate abundancies and for each bin
		Double_t abund[5];
		Double_t pT = hWeighMapRebin->GetXaxis()->GetBinCenter(i);
		for(int s=0; s<5; s++) abund[s] = hAbund[s]->GetBinContent(hAbund[s]->FindBin(pT));
		printf("DEBUG (EVAL STAGE): Bin %d), pT %.2f, abs (pi,K,p,e,mu): %.3f, %.3f, %.3f, %.3f, %.3f\n",i,pT,abund[0],abund[1],abund[2],abund[3],abund[4]);

		for(int j=1; j<=hMapRebin[0]->GetNbinsY(); j++) { //loop on eta bins
			for(int k=1; k<=hMapRebin[0]->GetNbinsZ(); k++) { //loop on zVtx bins
				//evaluate abundancies and for each bin
				Double_t weightedEff = 0., weightedEffErr = 0.;
				for(int s=0; s<5; s++) weightedEff += hMapRebin[s]->GetBinContent(i,j,k)*abund[s];
				for(int s=0; s<5; s++) weightedEffErr += hMapRebin[s]->GetBinError(i,j,k)*hMapRebin[s]->GetBinError(i,j,k)*abund[s]*abund[s];	
				weightedEffErr = TMath::Sqrt(weightedEffErr);
			    //set values in the map
				hWeighMapRebin->SetBinContent(i,j,k,weightedEff);
				hWeighMapRebin->SetBinError(i,j,k,weightedEffErr);
			    //DEBUG LINE
			    if(j<3 && k<3) printf("--> WEIGHTING BIN %d,%d,%d (pT %.3f)): ALL %.3f - pi %.3f (w %.3f) - K %.3f (w %.3f) - p %.3f (w %.3f) - e %.3f (w %.3f) - mu %.3f (w %.3f)\n",
			    	i,j,k,pT,weightedEff,hMapRebin[0]->GetBinContent(i,j,k),abund[0],hMapRebin[1]->GetBinContent(i,j,k),abund[1],hMapRebin[2]->GetBinContent(i,j,k),abund[2],hMapRebin[3]->GetBinContent(i,j,k),abund[3],hMapRebin[4]->GetBinContent(i,j,k),abund[4]);
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


void Reweight_1DTrackEff_By_Species(TString inputfolder="./", TString filemapprefix="1D_TrackingEffMap", TString filemapsuffix="18m", TString fileabund_data="pp_Monach13.root", TString fileabund_MC="RelAbundances_5Species.root") {

	TString species[5] = {"pi","K","p","e","mu"};

    //load input abundancies - data (pi, K, p)
	TFile *fAbundData  = new TFile(fileabund_data.Data(),"read");
    TH1D *hAbund[5];
    hAbund[0] = (TH1D*)fAbundData->Get("RelativeAbundancesData_kPion");    
    hAbund[1] = (TH1D*)fAbundData->Get("RelativeAbundancesData_kKaon");
    hAbund[2] = (TH1D*)fAbundData->Get("RelativeAbundancesData_kProton");

    //load input abundancies - MC (e, mu)
	TFile *fAbundMC  = new TFile(fileabund_MC.Data(),"read");
	TCanvas *cMCin = fAbundMC->Get("c1_n2");
    hAbund[3] = (TH1D*)cMCin->FindObject("containerpp13TeV_e_SelStep1_proj_0");    
    hAbund[4] = (TH1D*)cMCin->FindObject("containerpp13TeV_mu_SelStep1_proj_0");

    //APPLY REWEIGHTING OF pi, K, p FROM DATA SO THAT THEIR ABUNCANCY SUM REMAINS THE SAME IN MC (I.E. abound(pi+K+p)_data = abound(pi+K+p)_MC) 
    //THIS BECAUSE WE USE MC ABUNDANCIES FOR e AND mu (NOT MEASURED), SO THE TOTAL SUM OF 5 SPECIES ABUNDANCIES SHALL BE 1
    //HERE OF COURSE abound(pi+K+p)_MC IS 1-abound(e+mu)_MC SINCE IN THE MC ABUNDANCY PLOT ONLY THOSE 5 SPECIES ARE CONSIDERED
	for(int i=1; i<=hAbund[0]->GetNbinsX(); i++) { //loop on pT bins of abundancy plots
		Double_t initialPi = hAbund[0]->GetBinContent(i);
		Double_t initialK = hAbund[1]->GetBinContent(i);
		Double_t initialP = hAbund[2]->GetBinContent(i);
		Double_t MCe = hAbund[3]->GetBinContent(hAbund[3]->FindBin(hAbund[0]->GetBinCenter(i))); //hAbund 3,4 and 0,1,2, have different binnings!!
		Double_t MCmu = hAbund[4]->GetBinContent(hAbund[4]->FindBin(hAbund[0]->GetBinCenter(i)));
		Double_t correctedPi = initialPi/(initialPi+initialK+initialP)*(1-MCe-MCmu); //equivalent to initialPi/(initialPi+initialK+initialP)*(MCpi+MCK+MCp)
		Double_t correctedK  = initialK/(initialPi+initialK+initialP)*(1-MCe-MCmu);
		Double_t correctedP  = initialP/(initialPi+initialK+initialP)*(1-MCe-MCmu);
		hAbund[0]->SetBinContent(i,correctedPi);
		hAbund[1]->SetBinContent(i,correctedK);
		hAbund[2]->SetBinContent(i,correctedP);
		printf("DEBUG (CHANGE WEIGHT pi): Bin %d), pT %.2f, initial p,K,pi,MCe,MCmu: %.3f, %.3f, %.3f, %.3f, %.3f, final pi %.3f \n",i,hAbund[0]->GetBinCenter(i),initialPi,initialK,initialP,MCe,MCmu,hAbund[0]->GetBinContent(i));
		printf("DEBUG (CHANGE WEIGHT K) : Bin %d), pT %.2f, initial p,K,pi,MCe,MCmu: %.3f, %.3f, %.3f, %.3f, %.3f, final K  %.3f \n",i,hAbund[0]->GetBinCenter(i),initialPi,initialK,initialP,MCe,MCmu,hAbund[1]->GetBinContent(i));
		printf("DEBUG (CHANGE WEIGHT p) : Bin %d), pT %.2f, initial p,K,pi,MCe,MCmu: %.3f, %.3f, %.3f, %.3f, %.3f, final p  %.3f \n",i,hAbund[0]->GetBinCenter(i),initialPi,initialK,initialP,MCe,MCmu,hAbund[2]->GetBinContent(i));
	}

    //load input maps
	TFile *fMap[5];
	TCanvas *cMap[5];
	TH1D *hMap[5];
	TH1D *hMapRebin[5];	
	for(int i=0; i<5; i++) { //loop on species
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
		Double_t abund[5];
		Double_t pT = hWeighMap->GetXaxis()->GetBinCenter(i);
		for(int s=0; s<5; s++) abund[s] = hAbund[s]->GetBinContent(hAbund[s]->FindBin(pT));
		printf("DEBUG (EVAL STAGE): Bin %d), pT %.2f, abs (pi,K,p,e,mu): %.3f, %.3f, %.3f, %.3f, %.3f\n",i,pT,abund[0],abund[1],abund[2],abund[3],abund[4]);

		//evaluate abundancies and for each bin
		Double_t weightedEff = 0., weightedEffErr = 0.;
		for(int s=0; s<5; s++) weightedEff += hMap[s]->GetBinContent(i)*abund[s];
		for(int s=0; s<5; s++) weightedEffErr += hMap[s]->GetBinError(i)*hMap[s]->GetBinError(i)*abund[s]*abund[s];	
		weightedEffErr = TMath::Sqrt(weightedEffErr);
	    //set values in the map
		hWeighMap->SetBinContent(i,weightedEff);
		hWeighMap->SetBinError(i,weightedEffErr);
		//DEBUG LINE
		printf("--> WEIGHTING BIN %d (pT %.3f)): ALL %.3f - pi %.3f (w %.3f) - K %.3f (w %.3f) - p %.3f (w %.3f) - e %.3f (w %.3f) - mu %.3f (w %.3f)\n",
			    	i,pT,weightedEff,hMap[0]->GetBinContent(i),abund[0],hMap[1]->GetBinContent(i),abund[1],hMap[2]->GetBinContent(i),abund[2],hMap[3]->GetBinContent(i),abund[3],hMap[4]->GetBinContent(i),abund[4]);			
	}

	//reweight the 1D maps (rebinned)
	for(int i=1; i<=hMapRebin[0]->GetNbinsX(); i++) { //loop on pT bins
		//evaluate abundancies and for each bin
		Double_t abund[5];
		Double_t pT = hWeighMapRebin->GetXaxis()->GetBinCenter(i);
		for(int s=0; s<5; s++) abund[s] = hAbund[s]->GetBinContent(hAbund[s]->FindBin(pT));
		printf("DEBUG (EVAL STAGE): Bin %d), pT %.2f, abs (pi,K,p,e,mu): %.3f, %.3f, %.3f, %.3f, %.3f\n",i,pT,abund[0],abund[1],abund[2],abund[3],abund[4]);

		//evaluate abundancies and for each bin
		Double_t weightedEff = 0., weightedEffErr = 0.;
		for(int s=0; s<5; s++) weightedEff += hMapRebin[s]->GetBinContent(i)*abund[s];
		for(int s=0; s<5; s++) weightedEffErr += hMapRebin[s]->GetBinError(i)*hMapRebin[s]->GetBinError(i)*abund[s]*abund[s];	
		weightedEffErr = TMath::Sqrt(weightedEffErr);
	    //set values in the map
		hWeighMapRebin->SetBinContent(i,weightedEff);
		hWeighMapRebin->SetBinError(i,weightedEffErr);
		//DEBUG LINE
		printf("--> WEIGHTING BIN %d (pT %.3f)): ALL %.3f - pi %.3f (w %.3f) - K %.3f (w %.3f) - p %.3f (w %.3f) - e %.3f (w %.3f) - mu %.3f (w %.3f)\n",
			    	i,pT,weightedEff,hMap[0]->GetBinContent(i),abund[0],hMapRebin[1]->GetBinContent(i),abund[1],hMapRebin[2]->GetBinContent(i),abund[2],hMapRebin[3]->GetBinContent(i),abund[3],hMapRebin[4]->GetBinContent(i),abund[4]);			
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