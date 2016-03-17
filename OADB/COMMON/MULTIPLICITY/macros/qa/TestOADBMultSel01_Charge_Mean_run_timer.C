#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

#endif

void TestOADBMultSel01_Charge_Mean_run_timer(Int_t run = 0) {
     TStopwatch timer_all;
    TStopwatch timer_run;
    TStopwatch timer_event;
	timer_all.Start();
    TString gLibs[] =   {"STEER",
        "ANALYSIS", "ANALYSISalice", "ANALYSIScalib"};
    TString thislib = "lib";
    for(Int_t ilib = 0; ilib<4; ilib++){
        thislib="lib";
        thislib.Append(gLibs[ilib].Data());
        cout<<"Will load "<<thislib.Data()<<endl;
        gSystem->Load(thislib.Data());
    }
    gSystem->SetIncludePath("-I$ROOTSYS/include  -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    
    char name[10];
    char name_c[10];
    TObjArray Hlist_mean_value(0);
    TObjArray Hlist2(0);
	const Int_t nRuns = 3;
	
	Long64_t lNEv;
	Long_t lSelected = 0;
    Long_t lAll = 0;
    Long_t exeptVtxSelected = 0;
    
    Float_t Sel_dev_All = 0;
    Float_t Sel_dev_exeptVtx = 0;

    Float_t fAmplitude_V0A = 0;
    Float_t fAmplitude_V0C = 0;
    Float_t fAmplitude_V0Apartial = 0;
    Float_t fAmplitude_V0Cpartial = 0;
    Float_t fAmplitude_V0AEq = 0;
    Float_t fAmplitude_V0CEq = 0;
    Float_t fMultiplicity_AD=0;
    Float_t fMultiplicity_ADA=0;
    Float_t fMultiplicity_ADC=0;
    Int_t fEvSel_nSPDClusters=0; 
    Float_t fTriggerChargeA = 0;
    Float_t fTriggerChargeC = 0;
    Float_t fTriggerChargeV0M = 0;
    Float_t fZncEnergy =0;
    Float_t fZnaEnergy =0;
    Int_t fRunNumber;
    
    //Event Selection Variables
    Bool_t fEvSel_IsNotPileupInMultBins      = kFALSE ;
    TCut cut4="fEvSel_IsNotPileupInMultBins==1";
	
    Bool_t fEvSel_Triggered                  = kFALSE ;
    TCut cut1="fEvSel_Triggered==1";
	
    Bool_t fEvSel_INELgtZERO                 = kFALSE ;
    TCut cut5="fEvSel_INELgtZERO==1"; 
    
    Float_t fEvSel_VtxZ                      = 10.0 ;
    TCut cut3="fEvSel_VtxZCut=10.0";
	
    Int_t fRefMultEta8;

		
		//output added by T.Drozhzhova

	//run dependance 
	TH1F* hMeanValue_vs_Run_V0M = new TH1F("fMeanValue_vs_Run_V0M","fMeanValue_vs_Run_V0M",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_V0A = new TH1F("fMeanValue_vs_Run_V0A","fMeanValue_vs_Run_V0A",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_V0C = new TH1F("fMeanValue_vs_Run_V0C","fMeanValue_vs_Run_V0C",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_AD = new TH1F("fMeanValue_vs_Run_AD","fMeanValue_vs_Run_AD",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_ADA = new TH1F("fMeanValue_vs_Run_ADA","fMeanValue_vs_Run_ADA",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_ADC = new TH1F("fMeanValue_vs_Run_ADC","fMeanValue_vs_Run_ADC",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_V0MEq = new TH1F("fMeanValue_vs_Run_V0MEq","fMeanValue_vs_Run_V0MEq",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_V0AEq = new TH1F("fMeanValue_vs_Run_V0AEq","fMeanValue_vs_Run_V0AEq",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_V0CEq = new TH1F("fMeanValue_vs_Run_V0CEq","fMeanValue_vs_Run_V0CEq",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_V0Apartial = new TH1F("fMeanValue_vs_Run_V0Apartial","fMeanValue_vs_Run_V0Apartial",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_V0Cpartial = new TH1F("fMeanValue_vs_Run_V0Cpartial","fMeanValue_vs_Run_V0Cpartial",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_ChargeA = new TH1F("fMeanValue_vs_Run_ChargeA","fMeanValue_vs_Run_ChargeA",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_ChargeC = new TH1F("fMeanValue_vs_Run_ChargeC","fMeanValue_vs_Run_ChargeC",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_ChargeV0M = new TH1F("fMeanValue_vs_Run_ChargeV0M","fMeanValue_vs_Run_ChargeV0M",nRuns,0,nRuns);
	TH1F* hNevent_Vtx_vs_All= new TH1F("hNevent_Vtx_vs_All","hNevent_Vtx_vs_All",nRuns,0,nRuns);
	TH1F* hNevent_Vtx_exeptVtxSelected = new TH1F("hNevent_Vtx_exeptVtxSelected","hNevent_Vtx_vs_exeptVtxSelected",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_VtxZ= new TH1F("hMeanValue_vs_Run_VtxZ","hMeanValue_vs_Run_VtxZ",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_SPD= new TH1F("hMeanValue_vs_Run_nSPDClusters","hMeanValue_vs_Run_nSPDClusters",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_ZnaEnergy= new TH1F("hMeanValue_vs_Run_ZnaEnergy","hZnaEnergy",nRuns,0,nRuns);
	TH1F* hMeanValue_vs_Run_ZncEnergy= new TH1F("hMeanValue_vs_Run_ZncEnergy","hZncEnergy",nRuns,0,nRuns);
		
	Hlist_mean_value.Add(hNevent_Vtx_vs_All);
	Hlist_mean_value.Add(hNevent_Vtx_exeptVtxSelected);
	Hlist_mean_value.Add(hMeanValue_vs_Run_V0M); 
	Hlist_mean_value.Add(hMeanValue_vs_Run_V0A); 
	Hlist_mean_value.Add(hMeanValue_vs_Run_V0C); 
	Hlist_mean_value.Add(hMeanValue_vs_Run_AD);
	Hlist_mean_value.Add(hMeanValue_vs_Run_ADA);
	Hlist_mean_value.Add(hMeanValue_vs_Run_ADC);
	Hlist_mean_value.Add(hMeanValue_vs_Run_V0MEq);
	Hlist_mean_value.Add(hMeanValue_vs_Run_V0AEq);
	Hlist_mean_value.Add(hMeanValue_vs_Run_V0CEq);
	Hlist_mean_value.Add(hMeanValue_vs_Run_V0Apartial);
	Hlist_mean_value.Add(hMeanValue_vs_Run_V0Cpartial);
	Hlist_mean_value.Add(hMeanValue_vs_Run_ChargeA);
	Hlist_mean_value.Add(hMeanValue_vs_Run_ChargeC);
	Hlist_mean_value.Add(hMeanValue_vs_Run_ChargeV0M);
	Hlist_mean_value.Add(hMeanValue_vs_Run_SPD);
	Hlist_mean_value.Add(hMeanValue_vs_Run_VtxZ);
	Hlist_mean_value.Add(hMeanValue_vs_Run_ZncEnergy);
    Hlist_mean_value.Add(hMeanValue_vs_Run_ZnaEnergy);
      
   TString fFileName_temp;
   Int_t fLabel[nRuns]={225717,226062,225106};
   TString fFileNameCalibration = "Calibration_";
   TString fFileNameAnalysisResults  = "AnalysisResults_";
   TString fFileName[nRuns]={"225717_train39","226062_train39","225106_train_tatiana"};
   TFile * fCalibration;
   TFile * fAnalysisResults;
   TFile* output_file;
   
   //values
   
		Float_t valueV0M = 0;
		Float_t valueV0A = 0;
		Float_t valueV0C = 0;
		
		Float_t valueAD =0;
		Float_t valueADA =0;
		Float_t valueADC =0;
		
		Float_t valueV0Apartial = 0;
		Float_t valueV0Cpartial = 0;
		Float_t valueV0Mpartial = 0;
		
		Float_t valueV0AEq = 0;
		Float_t valueV0CEq = 0;
		Float_t valueV0MEq = 0;
		
		Int_t valueSPD = 0;
		Float_t valueVtxZ = 0;
		Float_t valueChargeA = 0;
		Float_t valueChargeC = 0;
		Float_t valueChargeV0M = 0;
		
		Float_t valueZncEnergy = 0;
		Float_t valueZnaEnergy = 0;
   
   //Runs loop
   for (Int_t i=0;i<nRuns;i++)
   {	timer_run.Start();
		cout << "fFileName[" << i << "] = " << fFileName[i] << endl;
		fFileName_temp = "/home/tatiana/home_GSI/CENTRALITY_GSI/output/" + fFileName[i] + "/" + fFileNameAnalysisResults + fFileName[i] + ".root";
		cout << "fFileName_temp = " << fFileName_temp << endl;
		fAnalysisResults = new TFile (fFileName_temp);
		
		sprintf(name,"%d",fLabel[i]);
		
		
		 // STEP 1: Basic I/O
		cout<<"(1) Opening File"<<endl;
		
		//Open File
		TFile *fInput = TFile::Open( fFileName_temp.Data(), "READ");
		if(!fInput){
			AliWarningF("File %s not found!", fFileName_temp.Data() );
			return kFALSE;
		}
		//Locate TTree object
		TTree* fTree = (TTree*)fInput->FindObjectAny("fTreeEvent");
		if(!fTree){
			AliWarning("fTreeEvent object not found!" );
			return kFALSE;
		}
				
		//SetBranchAddresses
		fTree->SetBranchAddress("fAmplitude_V0A",&fAmplitude_V0A);
		fTree->SetBranchAddress("fAmplitude_V0C",&fAmplitude_V0C);
		fTree->SetBranchAddress("fAmplitude_V0Apartial",&fAmplitude_V0Apartial);
		fTree->SetBranchAddress("fAmplitude_V0Cpartial",&fAmplitude_V0Cpartial);
		fTree->SetBranchAddress("fAmplitude_V0AEq",&fAmplitude_V0AEq);
		fTree->SetBranchAddress("fAmplitude_V0CEq",&fAmplitude_V0CEq);
		fTree->SetBranchAddress("fEvSel_IsNotPileupInMultBins",&fEvSel_IsNotPileupInMultBins);
		fTree->SetBranchAddress("fEvSel_Triggered",&fEvSel_Triggered);
		fTree->SetBranchAddress("fEvSel_INELgtZERO",&fEvSel_INELgtZERO);
		fTree->SetBranchAddress("fEvSel_VtxZ",&fEvSel_VtxZ);
		fTree->SetBranchAddress("fRunNumber",&fRunNumber);
		fTree->SetBranchAddress("fRefMultEta8",&fRefMultEta8);
		fTree->SetBranchAddress("fMultiplicity_AD",&fMultiplicity_AD);
		fTree->SetBranchAddress("fMultiplicity_ADA",&fMultiplicity_ADA);
		fTree->SetBranchAddress("fMultiplicity_ADC",&fMultiplicity_ADC);
		fTree->SetBranchAddress("fEvSel_nSPDClusters",&fEvSel_nSPDClusters);
		fTree->SetBranchAddress("fTriggerChargeA",&fTriggerChargeA);
		fTree->SetBranchAddress("fTriggerChargeC",&fTriggerChargeC);
		fTree->SetBranchAddress("fZncEnergy",&fZncEnergy);
		fTree->SetBranchAddress("fZnaEnergy",&fZnaEnergy);
		
		
				
		lNEv = fTree->GetEntries();
		cout<<"(1) File opened, event count is "<<lNEv<<endl;
		
		lSelected = 0;
		lAll = 0;
		exeptVtxSelected = 0;
			
		valueV0M = 0;
		valueV0A = 0;
		valueV0C = 0;
		
		valueAD =0;
		valueADA =0;
		valueADC =0;
		
		valueV0Apartial = 0;
		valueV0Cpartial = 0;
		valueV0Mpartial = 0;
		
		valueV0AEq = 0;
		valueV0CEq = 0;
		valueV0MEq = 0;
		
		valueSPD = 0;
		valueVtxZ = 0;
		
		valueChargeA = 0;
		valueChargeC = 0;
		valueChargeV0M = 0;
		
		valueZncEnergy =0;
		valueZnaEnergy =0;
		timer_event.Start();
		for(Long64_t iEv = 0; iEv<lNEv; iEv++)
		{
			fTree->GetEntry(iEv); //Look at next event
			lAll ++;
			
			//Check Selections as they are in the fMultSelectionCuts Object
			if( ! fEvSel_Triggered  ) continue;
			if( ! fEvSel_INELgtZERO ) continue;
			if( ! fEvSel_IsNotPileupInMultBins ) continue;
			if( fRefMultEta8 == -4  ) continue;
			exeptVtxSelected++;
			if( TMath::Abs(fEvSel_VtxZ) > 10 ) continue;
			lSelected ++;
			
			valueV0M = valueV0M + fAmplitude_V0A + fAmplitude_V0C;
			valueV0A = valueV0A + fAmplitude_V0A;
			valueV0C = valueV0C + fAmplitude_V0C;

			valueAD = valueAD + fMultiplicity_AD;
			valueADA = valueADA + fMultiplicity_ADA;
			valueADC = valueADC + fMultiplicity_ADC;

			valueV0Apartial = valueV0Apartial+fAmplitude_V0Apartial;
			valueV0Cpartial = valueV0Cpartial+fAmplitude_V0Cpartial;
			valueV0Mpartial = valueV0Mpartial+fAmplitude_V0Apartial+fAmplitude_V0Cpartial;

			valueV0AEq = valueV0AEq + fAmplitude_V0AEq;
			valueV0CEq = valueV0CEq + fAmplitude_V0CEq;
			valueV0MEq = valueV0MEq + fAmplitude_V0AEq+fAmplitude_V0CEq;

			valueSPD = valueSPD + fEvSel_nSPDClusters;
			valueVtxZ = valueVtxZ + fEvSel_VtxZ;

			valueChargeA = valueChargeA+fTriggerChargeA;
			valueChargeC = valueChargeC+fTriggerChargeC;
			valueChargeV0M = valueChargeV0M+fTriggerChargeA+fTriggerChargeC;
		
			valueZncEnergy = valueZncEnergy +fZncEnergy;
			valueZnaEnergy = valueZnaEnergy +fZnaEnergy;
			
		}
	timer_event.Stop();
	cout<<"Run:"<<fLabel[i]<<" timer run = "<<endl;
	timer_event.Print();
		
			cout<<"lSelected = "<<lSelected<< endl;
			cout<<valueV0M<<" = V0M "<<valueV0M/Float_t(lSelected)<<endl;
			cout<<valueAD<<" = AD "<<valueAD/Float_t(lSelected)<<endl;
			cout<<valueSPD<<" = SPD "<<valueSPD/Float_t(lSelected)<<endl;
			cout<<valueChargeV0M <<" = ChargeV0M "<< valueChargeV0M/Float_t(lSelected)<<endl;
			cout<<valueV0Cpartial<<"= valueV0Cpartial "<<valueV0Cpartial/Float_t(lSelected)<<endl;
			cout<< valueVtxZ<<" = valueVtxZ "<<valueVtxZ/Float_t(lSelected)<<endl;
			cout<<valueZncEnergy<<" = valueZncEnergy " <<valueZncEnergy/Float_t(lSelected)<<endl;
			cout<<valueZnaEnergy<<" = valueZnaEnergy " <<valueZnaEnergy/Float_t(lSelected)<<endl;
			
			//Mean 
			hMeanValue_vs_Run_SPD->Fill(i,valueSPD/Float_t(lSelected));
			hMeanValue_vs_Run_VtxZ->Fill(i,valueVtxZ/Float_t(lSelected));
			hMeanValue_vs_Run_V0M->Fill(i,valueV0M/Float_t(lSelected));
			hMeanValue_vs_Run_V0A->Fill(i,valueV0A/Float_t(lSelected));
			hMeanValue_vs_Run_V0C->Fill(i,valueV0C/Float_t(lSelected));
			hMeanValue_vs_Run_ZncEnergy->Fill(i,valueZncEnergy/Float_t(lSelected));
			hMeanValue_vs_Run_ZnaEnergy->Fill(i,valueZnaEnergy/Float_t(lSelected));
			
			if (fMultiplicity_AD)
			{
				hMeanValue_vs_Run_AD->Fill(i,valueAD/Float_t(lSelected));
				cout<<"AD"<<endl;
			}
			if (fMultiplicity_ADA)
			{
				hMeanValue_vs_Run_ADA->Fill(i,valueADA/Float_t(lSelected));
			}
			if (fMultiplicity_ADC)
			{
				hMeanValue_vs_Run_ADC->Fill(i,valueADC/Float_t(lSelected));
			}
			hMeanValue_vs_Run_V0MEq->Fill(i,valueV0MEq/Float_t(lSelected));
			cout<<"V0MEq" <<endl;
			hMeanValue_vs_Run_V0AEq->Fill(i,valueV0AEq/Float_t(lSelected));
			
			hMeanValue_vs_Run_V0CEq->Fill(i,valueV0CEq/Float_t(lSelected));
			
			hMeanValue_vs_Run_V0Apartial->Fill(i,valueV0Apartial/Float_t(lSelected));
			
			hMeanValue_vs_Run_V0Cpartial->Fill(i,valueV0Cpartial/Float_t(lSelected));
			cout<<"V0Cpartial" <<endl;
			if (fTriggerChargeA)
			{
				hMeanValue_vs_Run_ChargeA->Fill(i,valueChargeA/Float_t(lSelected));
			
			}
			if (fTriggerChargeC)
			{
				hMeanValue_vs_Run_ChargeC->Fill(i,valueChargeC/Float_t(lSelected));
						
			}
			if (fTriggerChargeA && fTriggerChargeC)
			{
				hMeanValue_vs_Run_ChargeV0M->Fill(i,valueChargeV0M/Float_t(lSelected));
			}
			cout<<"Fill all" <<endl;			
			fInput->Close();
			cout<<"Input close"<<endl;
			
			Sel_dev_All = Float_t(lSelected)/Float_t(lAll);
			cout<<"Sel_dev_All="<<Sel_dev_All<<endl;
			Sel_dev_exeptVtx = Float_t(lSelected)/Float_t(exeptVtxSelected);
			cout << "fFileName[" << i << "] : Sel_dev_All = " << Sel_dev_All<<"Sel_dev_exeptVtx ="<<Sel_dev_exeptVtx << endl;
				
			hNevent_Vtx_vs_All -> Fill(i,Sel_dev_All);
			hNevent_Vtx_exeptVtxSelected ->Fill(i,Sel_dev_exeptVtx);
		
		// Run label added to the histogram
				hMeanValue_vs_Run_ZncEnergy->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_ZnaEnergy->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_VtxZ->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_V0M->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_V0A->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_V0C->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_AD->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_ADA->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_ADC->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_V0MEq->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_V0AEq->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_V0CEq->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_V0Apartial->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_V0Cpartial->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_ChargeA->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_ChargeC->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_ChargeV0M->GetXaxis()->SetBinLabel(i+1,name);
				hMeanValue_vs_Run_SPD->GetXaxis()->SetBinLabel(i+1,name);
				hNevent_Vtx_vs_All->GetXaxis()->SetBinLabel(i+1,name);
				hNevent_Vtx_exeptVtxSelected->GetXaxis()->SetBinLabel(i+1,name);
	}
	timer_run.Stop();
		cout<<"timer event = "<<endl;
		timer_run.Print();
	TFile* output_file = new TFile("Percentile_results/Mean_run_distribution.root","RECREATE");
	Hlist_mean_value->Write(); 
	output_file->Close(); 
	
	TCanvas* myc0=new TCanvas("myc", "V0", 1350,1474);
    TCanvas* myc1=new TCanvas("myc1", "V0Eq", 1350,1474);
    TCanvas* myc2=new TCanvas("myc2", "V0partial", 1350,1474);
    TCanvas* myc3=new TCanvas("myc3", "AD", 1350,1474);
    TCanvas* myc4=new TCanvas("myc4", "SPD", 1350,1474);
    TCanvas* myc5=new TCanvas("myc5", "TriggerCharge", 1350,1474);
    TCanvas* myc6=new TCanvas("myc6", "VtxZ", 1350,1474);
    TCanvas* myc7=new TCanvas("myc7", "Vtx_vs_All", 1350,1474);
    TCanvas* myc8=new TCanvas("myc8", "Vtx_vs_exeptVtxSelected", 1350,1474);
    TCanvas* myc9=new TCanvas("myc9", "Znc", 1350,1474);
    
    myc2->Divide(1,2);
    myc9->Divide(1,2);
    
    myc6->cd();
				hMeanValue_vs_Run_VtxZ->Draw("e");
	myc0->cd();		
				hMeanValue_vs_Run_V0M->Draw("e");
	myc3->cd();			
				hMeanValue_vs_Run_AD->Draw("e");
	myc1->cd();			
				hMeanValue_vs_Run_V0MEq->Draw("e");
	myc2->cd(1);			
				hMeanValue_vs_Run_V0Apartial->Draw("e");
	myc2->cd(2);			
				hMeanValue_vs_Run_V0Cpartial->Draw("e");
	myc5->cd();			
				hMeanValue_vs_Run_ChargeV0M->Draw("e");
	myc4->cd();			
				hMeanValue_vs_Run_SPD->Draw("e");
	myc7->cd();			
				hNevent_Vtx_vs_All->Draw("e");
	myc8->cd();			
				hNevent_Vtx_exeptVtxSelected->Draw("e");
	myc9->cd(1);			
				hMeanValue_vs_Run_ZncEnergy->Draw("e");
	myc9->cd(2);			
				hMeanValue_vs_Run_ZnaEnergy->Draw("e");
	
	myc0->SaveAs("Percentile_results/Mean_value_hist_with_error/V0M.png");
	myc1->SaveAs("Percentile_results/Mean_value_hist_with_error/V0MEq.png");
	myc2->SaveAs("Percentile_results/Mean_value_hist_with_error/V0Apartial.png");
	myc3->SaveAs("Percentile_results/Mean_value_hist_with_error/AD.png");
	myc4->SaveAs("Percentile_results/Mean_value_hist_with_error/SPD.png");
	myc5->SaveAs("Percentile_results/Mean_value_hist_with_error/ChargeV0M.png");
	myc6->SaveAs("Percentile_results/Mean_value_hist_with_error/VtxZ.png");
	myc7->SaveAs("Percentile_results/Mean_value_hist_with_error/Nevent_Vtx_vs_All.png");
	myc8->SaveAs("Percentile_results/Mean_value_hist_with_error/Nevwnt_Vtx_vs_exeptVtxSelected.png");
	myc9->SaveAs("Percentile_results/Mean_value_hist_with_error/hMeanValue_vs_Run_ZncEnergy.png");
	timer_all.Stop();
	cout<<"Run:"<<fLabel[i]<<" timer run = "<<endl;
	timer_all.Print();
	
 return 0;   
}

    

    
    
