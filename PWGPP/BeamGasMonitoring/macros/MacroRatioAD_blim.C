// this is macro to get "efficiency", "purity" and "rejection efficiency".

Int_t bunchinputarray[7] = {201, 999, 999, 999, 999, 999, 999}; // this array contains our interesting output 2015.08.18. (blim)

void MacroRatioAD_blim(){
    for (Int_t irun=1; irun<2 ; irun++ ) {
        GetRatio(irun); // need to be modyfied. irun automatically read run number and produce output for each run
    }
}
GetRatio(Int_t split=-1){
    
    
    
    TString ss = Form("%03d",split);
    cout << ss << endl;
    
    
    TFile *filein = TFile::Open(ss+"/AnalysisResults.root","READ");
    if(filein==NULL) continue;
    cout << "open root file : " << ss << "/AnalysisResults.root"<< endl;

    TDirectoryFile *dir = (TDirectoryFile*)filein->Get("BeamGasMon");
    
    TList *list = (TList*)dir->Get("cOutputH_MB");
    TList *list2 = (TList*)dir->Get("cOutputH_HM"); //add for HM data
    /*
    TH2F *hTotalTrkVsClsSPID = (TH2F*)list->FindObject("hTotalTrkVsClsSPID");
    hTotalTrkVsClsSPID->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID->GetYaxis()->SetTitle("Cluster (fspdC1)");
    */
    TH2F *hTotalV0 = (TH2F*)list->FindObject("hTotalV0");
    hTotalV0->GetXaxis()->SetTitle("V0A-V0C");
    hTotalV0->GetYaxis()->SetTitle("V0A+V0C");
    
    TH2F *hTotalAD = (TH2F*)list->FindObject("hTotalAD");
    hTotalAD->GetXaxis()->SetTitle("ADA-ADC");
    hTotalAD->GetYaxis()->SetTitle("ADA+ADC");
  
    TH1F *hNumEffPurityBC[3][3][3];
    TH1F *hDenomEffBC[3][3][3];
    TH1F *hDenomPurityBC[3][3][3];
    TH1F *hDenomRejecEffBC[3][3][3];
    TH1F *hNumRejecEffBC[3][3][3];
    TH1F *hSPDNumBC[3][3][3];
    TH1F *hSPDDenomBC[3][3][3];
    
    TH2F *hNumTrkVsClsSPID[3][3][3];
    TH2F *hDenomTrkVsClsSPID[3][3][3];
    TH2F *hNumV0[3][3][3];
    TH2F *hDenomV0[3][3][3];
    
    
    TH1F *hADNumEffPurityBC[3][3][3];
    TH1F *hADDenomPurityBC[3][3][3];
    TH1F *hADNumRejecEffBC[3][3][3];
    TH2F *hADNumTrkVsClsSPD[3][3][3];
    TH2F *hNumAD[3][3][3];
    TH2F *hDenomAD[3][3][3];

    TH1F *hNumEffPurityBC_HM[3][3][3];
    TH1F *hDenomEffBC_HM[3][3][3];
    TH1F *hDenomPurityBC_HM[3][3][3];
    TH1F *hDenomRejecEffBC_HM[3][3][3];
    TH1F *hNumRejecEffBC_HM[3][3][3];
    TH1F *hSPDNumBC_HM[3][3][3];
    TH1F *hSPDDenomBC_HM[3][3][3];
    
    TH2F *hNumTrkVsClsSPID_HM[3][3][3];
    TH2F *hDenomTrkVsClsSPID_HM[3][3][3];
    TH2F *hNumV0_HM[3][3][3];
    TH2F *hDenomV0_HM[3][3][3];
    
    
    TH1F *hADNumEffPurityBC_HM[3][3][3];
    TH1F *hADDenomPurityBC_HM[3][3][3];
    TH1F *hADNumRejecEffBC_HM[3][3][3];
    TH2F *hADNumTrkVsClsSPD_HM[3][3][3];
    TH2F *hNumAD_HM[3][3][3];
    TH2F *hDenomAD_HM[3][3][3];
    
    
     TH2F *hTrkVsClsSPIDSlopeM = (TH2F*)list->FindObject("hTrkVsClsSPIDSlopeM");
     hTrkVsClsSPIDSlopeM->GetXaxis()->SetTitle("Tracklet");
     hTrkVsClsSPIDSlopeM->GetYaxis()->SetTitle("Cluster (fspdC1+fsdpC2)");
    TH2F *hTrkVsClsSPIDSlopeM2 = (TH2F*)list->FindObject("hTrkVsClsSPIDSlopeM2");
    hTrkVsClsSPIDSlopeM2->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM2->GetYaxis()->SetTitle("Cluster (fspdC1+fsdpC2)");
    
    TH2F *hTrkVsClsSPIDSlopeM_HM = (TH2F*)list2->FindObject("hTrkVsClsSPIDSlopeM_HM");
    hTrkVsClsSPIDSlopeM_HM->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM_HM->GetYaxis()->SetTitle("Cluster (fspdC1+fsdpC2)");
    TH2F *hTrkVsClsSPIDSlopeM_HM2 = (TH2F*)list2->FindObject("hTrkVsClsSPIDSlopeM_HM2");
    hTrkVsClsSPIDSlopeM_HM2->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM_HM2->GetYaxis()->SetTitle("Cluster (fspdC1+fsdpC2)");
//    hTrkVsClsSPIDSlopeM->SaveAs("./Fig/hTrkVsClsSPIDSlopeM.pdf");
    //addd
    
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
    for(int k=0; k<3; k++){
        //________ we want to make histograms which we use really. 2015.08.18. (blim)
        Int_t check000 = i*100+j*10+k;
        for ( int l=0; l<8;l++){
            if(check000==bunchinputarray[l]){
        //________
            
        hNumEffPurityBC[i][j][k] =(TH1F*)list->FindObject(Form("hNumEffPurityBC%d_V0%01d_Flag%d",i,j,k));
        hDenomEffBC[i][j][k]=(TH1F*)list->FindObject(Form("hDenomEffBC%d_V0%01d_Flag%d",i,j,k));
        hDenomPurityBC[i][j][k]=(TH1F*)list->FindObject(Form("hDenomPurityBC%d_V0%01d_Flag%d",i,j,k));
        hDenomRejecEffBC[i][j][k]=(TH1F*)list->FindObject(Form("hDenomRejecEffBC%d_V0%01d_Flag%d",i,j,k));
        hNumRejecEffBC[i][j][k]=(TH1F*)list->FindObject(Form("hNumRejecEffBC%d_V0%01d_Flag%d",i,j,k));
        hSPDNumBC[i][j][k]=(TH1F*)list->FindObject(Form("hSPDNumBC%d_V0%01d_Flag%d",i,j,k));
        hSPDDenomBC[i][j][k]=(TH1F*)list->FindObject(Form("hSPDDenomBC%d_V0%01d_Flag%d",i,j,k));
        
        hNumTrkVsClsSPID[i][j][k] = (TH2F*)list->FindObject(Form("hNumTrkVsClsSPID%d_V0%d_Flag%d",i,j,k));
        hDenomTrkVsClsSPID[i][j][k]=(TH2F*)list->FindObject(Form("hDenomTrkVsClsSPID%d_V0%d_Flag%d",i,j,k));
        hNumV0[i][j][k] = (TH2F*)list->FindObject(Form("hNumV0%d_V0%d_Flag%d",i,j,k));
        hDenomV0[i][j][k] = (TH2F*)list->FindObject(Form("hDenomV0%d_V0%d_Flag%d",i,j,k));
        
        hADNumEffPurityBC[i][j][k] =(TH1F*)list->FindObject(Form("hADNumEffPurityBC%d_AD%01d_Flag%d",i,j,k));
        hADDenomPurityBC[i][j][k] =(TH1F*)list->FindObject(Form("hADDenomPurityBC%d_AD%01d_Flag%d",i,j,k));
        hADNumRejecEffBC[i][j][k] =(TH1F*)list->FindObject(Form("hADNumRejecEffBC%d_AD%01d_Flag%d",i,j,k));
        hADNumTrkVsClsSPD[i][j][k] = (TH2F*)list->FindObject(Form("hADNumTrkVsClsSPD%d_AD%d_Flag%d",i,j,k));
        hNumAD[i][j][k] = (TH2F*)list->FindObject(Form("hNumAD%d_AD%d_Flag%d",i,j,k));
        hDenomAD[i][j][k] = (TH2F*)list->FindObject(Form("hDenomAD%d_AD%d_Flag%d",i,j,k));
            }
        }
        
    }
        }
    } // read histogram from output

    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                //________ we want to make histograms which we use really. 2015.08.18. (blim)
                Int_t check000 = i*100+j*10+k;
                for ( int l=0; l<8;l++){
                    if(check000==bunchinputarray[l]){
                        //________
                        
                        hNumEffPurityBC_HM[i][j][k] =(TH1F*)list2->FindObject(Form("hNumEffPurityBC_HM%d_V0%01d_Flag%d",i,j,k));
                        hDenomEffBC_HM[i][j][k]=(TH1F*)list2->FindObject(Form("hDenomEffBC_HM%d_V0%01d_Flag%d",i,j,k));
                        hDenomPurityBC_HM[i][j][k]=(TH1F*)list2->FindObject(Form("hDenomPurityBC_HM%d_V0%01d_Flag%d",i,j,k));
                        hDenomRejecEffBC_HM[i][j][k]=(TH1F*)list2->FindObject(Form("hDenomRejecEffBC_HM%d_V0%01d_Flag%d",i,j,k));
                        hNumRejecEffBC_HM[i][j][k]=(TH1F*)list2->FindObject(Form("hNumRejecEffBC_HM%d_V0%01d_Flag%d",i,j,k));
                        hSPDNumBC_HM[i][j][k]=(TH1F*)list2->FindObject(Form("hSPDNumBC_HM%d_V0%01d_Flag%d",i,j,k));
                        hSPDDenomBC_HM[i][j][k]=(TH1F*)list2->FindObject(Form("hSPDDenomBC_HM%d_V0%01d_Flag%d",i,j,k));
                        
                        hNumTrkVsClsSPID_HM[i][j][k] = (TH2F*)list2->FindObject(Form("hNumTrkVsClsSPID_HM%d_V0%d_Flag%d",i,j,k));
                        hDenomTrkVsClsSPID_HM[i][j][k]=(TH2F*)list2->FindObject(Form("hDenomTrkVsClsSPID_HM%d_V0%d_Flag%d",i,j,k));
                        hNumV0_HM[i][j][k] = (TH2F*)list2->FindObject(Form("hNumV0_HM%d_V0%d_Flag%d",i,j,k));
                        hDenomV0_HM[i][j][k] = (TH2F*)list2->FindObject(Form("hDenomV0_HM%d_V0%d_Flag%d",i,j,k));
                        
                        hADNumEffPurityBC_HM[i][j][k] =(TH1F*)list2->FindObject(Form("hADNumEffPurityBC_HM%d_AD%01d_Flag%d",i,j,k));
                        hADDenomPurityBC_HM[i][j][k] =(TH1F*)list2->FindObject(Form("hADDenomPurityBC_HM%d_AD%01d_Flag%d",i,j,k));
                        hADNumRejecEffBC_HM[i][j][k] =(TH1F*)list2->FindObject(Form("hADNumRejecEffBC_HM%d_AD%01d_Flag%d",i,j,k));
                        hADNumTrkVsClsSPD_HM[i][j][k] = (TH2F*)list2->FindObject(Form("hADNumTrkVsClsSPD_HM%d_AD%d_Flag%d",i,j,k));
                        hNumAD_HM[i][j][k] = (TH2F*)list2->FindObject(Form("hNumAD_HM%d_AD%d_Flag%d",i,j,k));
                        hDenomAD_HM[i][j][k] = (TH2F*)list2->FindObject(Form("hDenomAD_HM%d_AD%d_Flag%d",i,j,k));
                    }
                }
                
            }
        }
    } // read histogram from output
    
    
    TH1F *hPurityBC[3][3][3];
    TH1F *hEfficiencyBC[3][3][3];
    TH1F *hRejectEfficiencyBC[3][3][3];
    TH1F *hSPDEfficiencyBC[3][3][3];
    TH2F *h2DSPDEfficiency[3][3][3];
    TH2F *h2DV0Efficiency[3][3][3];


    TH1F *hADPurityBC[3][3][3];
    TH1F *hADEfficiencyBC[3][3][3];
    TH1F *hADRejectEfficiencyBC[3][3][3];
    TH2F *h2DADEfficiency[3][3][3];

    double denomeff = hDenomEffBC[2][0][1]->GetBinContent(10);
    double denomrejeceff = hDenomRejecEffBC[2][0][1]->GetBinContent(10);

    for(int i=0; i<3; i++){
       for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                Int_t check000 = i*100+j*10+k;
                for ( int l=0; l<8;l++){
                    if(check000==bunchinputarray[l]){
    hPurityBC[i][j][k] = new TH1F(Form("hPurityBC%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
                hPurityBC[i][j][k]->SetTitle("Purity");
    hEfficiencyBC[i][j][k] = new TH1F(Form("hEfficiencyBC%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
                hEfficiencyBC[i][j][k]->SetTitle("Efficiency");

    hRejectEfficiencyBC[i][j][k] = new TH1F(Form("hRejectEfficiencyBC%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
                hRejectEfficiencyBC[i][j][k]->SetTitle("Rejection-Efficiency");

                
    hSPDEfficiencyBC[i][j][k] = new TH1F(Form("hSPDEfficiencyBC%d_V0%d_Flag%d",i,j,k),"; Spd tracklet", 200, 0, 200);

    h2DSPDEfficiency[i][j][k] = new TH2F(Form("h2DSPDEfficiency%d_V0%d_Flag%d",i,j,k),"; Spd Cluster vs tracklet",140,0,140,500,0,500);
    h2DSPDEfficiency[i][j][k]->GetXaxis()->SetTitle("Tracklet");
    h2DSPDEfficiency[i][j][k]->GetYaxis()->SetTitle("Cluster (fspdC1)");
   
    h2DV0Efficiency[i][j][k] = new TH2F(Form("h2DV0Efficiency%d_V0%d_Flag%d",i,j,k),"; V0 timing information", 600,-300,300,2000,-1000,1000);
    h2DV0Efficiency[i][j][k]->GetXaxis()->SetTitle("V0A-V0C");
    h2DV0Efficiency[i][j][k]->GetYaxis()->SetTitle("V0A+V0C");
    
                
    
    hADPurityBC[i][j][k] = new TH1F(Form("hADPurityBC%d_V0%d_Flag%d",i,j,k),"; #ADflags in PF", 10, 0, 10);
    hADPurityBC[i][j][k]->SetTitle("AD : Purity");
    hADEfficiencyBC[i][j][k] = new TH1F(Form("hADEfficiencyBC%d_V0%d_Flag%d",i,j,k),"; #ADflags in PF", 10, 0, 10);
    hADEfficiencyBC[i][j][k]->SetTitle("AD : Efficiency");
    hADRejectEfficiencyBC[i][j][k] = new TH1F(Form("hADRejectEfficiencyBC%d_V0%d_Flag%d",i,j,k),"; #ADflags in PF", 10, 0, 10);
    hADRejectEfficiencyBC[i][j][k]->SetTitle("AD : Rejection-Efficiency");

                
          

                
    
    hPurityBC[i][j][k]->Divide(hNumEffPurityBC[i][j][k], hDenomPurityBC[i][j][k], 1, 1, "B");
                hPurityBC[i][j][k]->GetYaxis()->SetRangeUser(0,1.0);

    hEfficiencyBC[i][j][k]->Divide(hNumEffPurityBC[i][j][k], hDenomEffBC[i][j][k], 1, 1, "B");
                hEfficiencyBC[i][j][k]->GetYaxis()->SetRangeUser(0,1.0);

    hRejectEfficiencyBC[i][j][k]->Divide(hNumRejecEffBC[i][j][k], hDenomRejecEffBC[i][j][k], 1, 1, "B");
                hRejectEfficiencyBC[i][j][k]->GetYaxis()->SetRangeUser(0,1.0);

    hSPDEfficiencyBC[i][j][k]->Divide(hSPDNumBC[i][j][k], hSPDDenomBC[i][j][k], 1, 1, "B");
    hSPDEfficiencyBC[i][j][k]->GetYaxis()->SetRangeUser(0,1.1);
    
   h2DSPDEfficiency[i][j][k]->Divide(hNumTrkVsClsSPID[i][j][k], hDenomTrkVsClsSPID[i][j][k], 1, 1, "B");
   h2DV0Efficiency[i][j][k]->Divide(hNumV0[i][j][k], hDenomV0[i][j][k], 1, 1, "B");
        
                hADPurityBC[i][j][k]->Divide(hADNumEffPurityBC[i][j][k], hADDenomPurityBC[i][j][k], 1, 1, "B");
                
            for(int iad=0; iad<10; iad++){
                    hADEfficiencyBC[i][j][k]->SetBinContent(iad, (hADNumEffPurityBC[i][j][k]->GetBinContent(iad))/denomeff);
                    hADRejectEfficiencyBC[i][j][k]->SetBinContent(iad, (hADNumRejecEffBC[i][j][k]->GetBinContent(iad))/denomrejeceff);
                }
                    }
                }
            }
        }
    }

    TLatex *texbc1 = new TLatex(3. , 0.7,"BC after : 3-#color[4]{7}, before : 11-17");
    texbc1->SetTextSize(.06);
    texbc1->SetTextColor(1);
    TLatex *texbc2 = new TLatex(3. , 0.7,"BC after : 3-#color[8]{8}, before : 11-17");
    texbc2->SetTextSize(.06);
    texbc2->SetTextColor(1);
    TLatex *texbc3 = new TLatex(3. , 0.7,"BC after : 3-#color[2]{9}, before : 11-17");
    texbc3->SetTextSize(.06);
    texbc3->SetTextColor(1);
    
    
    TLatex *texv01 = new TLatex(3. , 0.8,"#color[8]{V0A + V0C}");
    texv01->SetTextSize(.06);
    texv01->SetTextColor(1);
    TLatex *texv02 = new TLatex(3. , 0.8,"#color[2]{V0A} only");
    texv02->SetTextSize(.06);
    texv02->SetTextColor(1);
    TLatex *texv03 = new TLatex(3. , 0.8,"#color[4]{V0C} only");
    texv03->SetTextSize(.06);
    texv03->SetTextColor(1);
    
    TLatex *texf1 = new TLatex(3. , 0.9,"#color[8]{BB + BG} flag");
    texf1->SetTextSize(.06);
    texf1->SetTextColor(1);
    TLatex *texf2 = new TLatex(3. , 0.9,"#color[2]{BB} flag only");
    texf2->SetTextSize(.06);
    texf2->SetTextColor(1);
    TLatex *texf3 = new TLatex(3. , 0.9,"#color[4]{BG} flag only");
    texf3->SetTextSize(.06);
    texf3->SetTextColor(1);
    
    //______Add some typical plot (blim)
    cOutCanvas_output= new TCanvas("cOutCanvas_output","cOutCanvas_output",1000,500);
    cOutCanvas_output->Divide(4,1);
    cOutCanvas_output->cd(1);
    hPurityBC[2][0][1]->Draw("same");
    //    hADPurityBC[2][0][1]->SetLineColor(2);
    //    hADPurityBC[2][0][1]->Draw("same");
    texbc3->Draw("same");
    texv01->Draw("same");
    texf2->Draw("same");
    
    cOutCanvas_output->cd(2);
    hEfficiencyBC[2][0][1]->Draw("same");
    //    hADEfficiencyBC[2][0][1]->SetLineColor(2);
    //    hADEfficiencyBC[2][0][1]->Draw("same");
    texbc3->Draw("same");
    texv01->Draw("same");
    texf2->Draw("same");
    //    texNum->Draw("same");
    
    cOutCanvas_output->cd(3);
    hRejectEfficiencyBC[2][0][1]->Draw("same");
    //    hADRejectEfficiencyBC[2][0][1]->SetLineColor(2);
    //    hADRejectEfficiencyBC[2][0][1]->Draw("same");
    texbc3->Draw("same");
    texv01->Draw("same");
    texf2->Draw("same");
    //    texNum->Draw("same");
    
    cOutCanvas_output->cd(4);
    hSPDEfficiencyBC[2][0][1]->Draw("same");
    texbc3->Draw("same");
    texv01->Draw("same");
    texf2->Draw("same");
    //    texNum->Draw("same");
    
    
    cOutCanvas_output->SaveAs("./Fig/cOutCanvas_output.pdf");
    
    //____________________________
    
    
    TH1F *hPurityBC_HM[3][3][3];
    TH1F *hEfficiencyBC_HM[3][3][3];
    TH1F *hRejectEfficiencyBC_HM[3][3][3];
    TH1F *hSPDEfficiencyBC_HM[3][3][3];
    TH2F *h2DSPDEfficiency_HM[3][3][3];
    TH2F *h2DV0Efficiency_HM[3][3][3];
    
    
    TH1F *hADPurityBC_HM[3][3][3];
    TH1F *hADEfficiencyBC_HM[3][3][3];
    TH1F *hADRejectEfficiencyBC_HM[3][3][3];
    TH2F *h2DADEfficiency_HM[3][3][3];
    
    double denomeff_HM = hDenomEffBC_HM[2][0][1]->GetBinContent(10);
    double denomrejeceff_HM = hDenomRejecEffBC_HM[2][0][1]->GetBinContent(10);
    
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                Int_t check000 = i*100+j*10+k;
                for ( int l=0; l<8;l++){
                    if(check000==bunchinputarray[l]){
                        hPurityBC_HM[i][j][k] = new TH1F(Form("hPurityBC_HM%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
                        hPurityBC_HM[i][j][k]->SetTitle("Purity");
                        hEfficiencyBC_HM[i][j][k] = new TH1F(Form("hEfficiencyBC_HM%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
                        hEfficiencyBC_HM[i][j][k]->SetTitle("Efficiency");
                        
                        hRejectEfficiencyBC_HM[i][j][k] = new TH1F(Form("hRejectEfficiencyBC_HM%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
                        hRejectEfficiencyBC_HM[i][j][k]->SetTitle("Rejection-Efficiency");
                        
                        
                        hSPDEfficiencyBC_HM[i][j][k] = new TH1F(Form("hSPDEfficiencyBC_HM%d_V0%d_Flag%d",i,j,k),"; Spd tracklet", 200, 0, 200);
                        
                        h2DSPDEfficiency_HM[i][j][k] = new TH2F(Form("h2DSPDEfficiency_HM%d_V0%d_Flag%d",i,j,k),"; Spd Cluster vs tracklet",140,0,140,500,0,500);
                        h2DSPDEfficiency_HM[i][j][k]->GetXaxis()->SetTitle("Tracklet");
                        h2DSPDEfficiency_HM[i][j][k]->GetYaxis()->SetTitle("Cluster (fspdC1)");
                        
                        h2DV0Efficiency_HM[i][j][k] = new TH2F(Form("h2DV0Efficiency_HM%d_V0%d_Flag%d",i,j,k),"; V0 timing information", 600,-300,300,2000,-1000,1000);
                        h2DV0Efficiency_HM[i][j][k]->GetXaxis()->SetTitle("V0A-V0C");
                        h2DV0Efficiency_HM[i][j][k]->GetYaxis()->SetTitle("V0A+V0C");
                        
                        
                        
                        hADPurityBC_HM[i][j][k] = new TH1F(Form("hADPurityBC_HM%d_V0%d_Flag%d",i,j,k),"; #ADflags in PF", 10, 0, 10);
                        hADPurityBC_HM[i][j][k]->SetTitle("AD : Purity");
                        hADEfficiencyBC_HM[i][j][k] = new TH1F(Form("hADEfficiencyBC_HM%d_V0%d_Flag%d",i,j,k),"; #ADflags in PF", 10, 0, 10);
                        hADEfficiencyBC_HM[i][j][k]->SetTitle("AD : Efficiency");
                        hADRejectEfficiencyBC_HM[i][j][k] = new TH1F(Form("hADRejectEfficiencyBC_HM%d_V0%d_Flag%d",i,j,k),"; #ADflags in PF", 10, 0, 10);
                        hADRejectEfficiencyBC_HM[i][j][k]->SetTitle("AD : Rejection-Efficiency");
                        
                        
                        
                        
                        
                        
                        hPurityBC_HM[i][j][k]->Divide(hNumEffPurityBC_HM[i][j][k], hDenomPurityBC_HM[i][j][k], 1, 1, "B");
                        hPurityBC_HM[i][j][k]->GetYaxis()->SetRangeUser(0,1.0);
                        
                        hEfficiencyBC_HM[i][j][k]->Divide(hNumEffPurityBC_HM[i][j][k], hDenomEffBC_HM[i][j][k], 1, 1, "B");
                        hEfficiencyBC_HM[i][j][k]->GetYaxis()->SetRangeUser(0,1.0);
                        
                        hRejectEfficiencyBC_HM[i][j][k]->Divide(hNumRejecEffBC_HM[i][j][k], hDenomRejecEffBC_HM[i][j][k], 1, 1, "B");
                        hRejectEfficiencyBC_HM[i][j][k]->GetYaxis()->SetRangeUser(0,1.0);
                        
                        hSPDEfficiencyBC_HM[i][j][k]->Divide(hSPDNumBC_HM[i][j][k], hSPDDenomBC_HM[i][j][k], 1, 1, "B");
                        hSPDEfficiencyBC_HM[i][j][k]->GetYaxis()->SetRangeUser(0,1.1);
                        
                        h2DSPDEfficiency_HM[i][j][k]->Divide(hNumTrkVsClsSPID_HM[i][j][k], hDenomTrkVsClsSPID_HM[i][j][k], 1, 1, "B");
                        h2DV0Efficiency_HM[i][j][k]->Divide(hNumV0_HM[i][j][k], hDenomV0_HM[i][j][k], 1, 1, "B");
                        
                        hADPurityBC_HM[i][j][k]->Divide(hADNumEffPurityBC_HM[i][j][k], hADDenomPurityBC_HM[i][j][k], 1, 1, "B");
                        
                        for(int iad=0; iad<10; iad++){
                            hADEfficiencyBC_HM[i][j][k]->SetBinContent(iad, (hADNumEffPurityBC_HM[i][j][k]->GetBinContent(iad))/denomeff_HM);
                            hADRejectEfficiencyBC_HM[i][j][k]->SetBinContent(iad, (hADNumRejecEffBC_HM[i][j][k]->GetBinContent(iad))/denomrejeceff_HM);
                        }
                    }
                }
            }
        }
    }
    
    
    //______Add some typical plot for HM (blim) 2015.11.09
    cOutCanvas_output_HM= new TCanvas("cOutCanvas_output_HM","cOutCanvas_output_HM",1000,500);
    cOutCanvas_output_HM->Divide(4,1);
    cOutCanvas_output_HM->cd(1);
    hPurityBC_HM[2][0][1]->Draw("same");
    //    hADPurityBC[2][0][1]->SetLineColor(2);
    //    hADPurityBC[2][0][1]->Draw("same");
    texbc3->Draw("same");
    texv01->Draw("same");
    texf2->Draw("same");
    
    cOutCanvas_output_HM->cd(2);
    hEfficiencyBC_HM[2][0][1]->Draw("same");
    //    hADEfficiencyBC[2][0][1]->SetLineColor(2);
    //    hADEfficiencyBC[2][0][1]->Draw("same");
    texbc3->Draw("same");
    texv01->Draw("same");
    texf2->Draw("same");
    //    texNum->Draw("same");
    
    cOutCanvas_output_HM->cd(3);
    hRejectEfficiencyBC_HM[2][0][1]->Draw("same");
    //    hADRejectEfficiencyBC[2][0][1]->SetLineColor(2);
    //    hADRejectEfficiencyBC[2][0][1]->Draw("same");
    texbc3->Draw("same");
    texv01->Draw("same");
    texf2->Draw("same");
    //    texNum->Draw("same");
    
    cOutCanvas_output_HM->cd(4);
    hSPDEfficiencyBC_HM[2][0][1]->Draw("same");
    texbc3->Draw("same");
    texv01->Draw("same");
    texf2->Draw("same");
    //    texNum->Draw("same");
    
    
    cOutCanvas_output_HM->SaveAs("./Fig/cOutCanvas_output_HM.pdf");
    
    //______

    /* I think that we don't need this part, now. 2015.08.18.(blim)
    TCanvas *cOutCanvas[3][3][3];
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
    cOutCanvas[i][j][k]= new TCanvas(Form("cOutCanvas%d_%d_%d",i,j,k),Form("cOutCanvas%d_%d_%d",i,j,k),1000,500);
    cOutCanvas[i][j][k]->Divide(4,1);
    cOutCanvas[i][j][k]->cd(1);
    hPurityBC[i][j][k]->Draw("same");
    hADPurityBC[i][j][k]->SetLineColor(2);
    hADPurityBC[i][j][k]->Draw("same");
            
    if(i==0)texbc1->Draw("same");
    if(i==1)texbc2->Draw("same");
    if(i==2)texbc3->Draw("same");
    if(j==0)texv01->Draw("same");
    if(j==1)texv02->Draw("same");
    if(j==2)texv03->Draw("same");
    if(k==0)texf1->Draw("same");
    if(k==1)texf2->Draw("same");
    if(k==2)texf3->Draw("same");
    cOutCanvas[i][j][k]->cd(2);
    hEfficiencyBC[i][j][k]->Draw("same");
                hADEfficiencyBC[i][j][k]->SetLineColor(2);
                hADEfficiencyBC[i][j][k]->Draw("same");
                if(i==0)texbc1->Draw("same");
                if(i==1)texbc2->Draw("same");
                if(i==2)texbc3->Draw("same");
                if(j==0)texv01->Draw("same");
                if(j==1)texv02->Draw("same");
                if(j==2)texv03->Draw("same");
                if(k==0)texf1->Draw("same");
                if(k==1)texf2->Draw("same");
                if(k==2)texf3->Draw("same");

    cOutCanvas[i][j][k]->cd(3);
    hRejectEfficiencyBC[i][j][k]->Draw("same");
                hADRejectEfficiencyBC[i][j][k]->SetLineColor(2);
                hADRejectEfficiencyBC[i][j][k]->Draw("same");
                if(i==0)texbc1->Draw("same");
                if(i==1)texbc2->Draw("same");
                if(i==2)texbc3->Draw("same");
                if(j==0)texv01->Draw("same");
                if(j==1)texv02->Draw("same");
                if(j==2)texv03->Draw("same");
                if(k==0)texf1->Draw("same");
                if(k==1)texf2->Draw("same");
                if(k==2)texf3->Draw("same");

    cOutCanvas[i][j][k]->cd(4);
    hSPDEfficiencyBC[i][j][k]->Draw("same");
                if(i==0)texbc1->Draw("same");
                if(i==1)texbc2->Draw("same");
                if(i==2)texbc3->Draw("same");
                if(j==0)texv01->Draw("same");
                if(j==1)texv02->Draw("same");
                if(j==2)texv03->Draw("same");
                if(k==0)texf1->Draw("same");
                if(k==1)texf2->Draw("same");
                if(k==2)texf3->Draw("same");

                
//     cOutCanvas[i][j][k]->SaveAs(Form("./Fig/cOutCanvas%d_%d_%d.pdf",i,j,k));
                    }
        }
    }*/
    
    TLatex *texNum = new TLatex(-15,25,"Numerator (!bgID & Good Event)");
    texNum->SetTextSize(.06);
    texNum->SetTextColor(1);
    TLatex *texDe = new TLatex(-15,25,"Denominator (!bgID)");
    texDe->SetTextSize(.06);
    texDe->SetTextColor(1);
    TLatex *texR = new TLatex(-15,25,"Ratio");
    texR->SetTextSize(.06);
    texR->SetTextColor(1);
    
    TLatex *texbc1 = new TLatex(-15,-15,"BC after : 3-#color[4]{7}, before : 11-17");
    texbc1->SetTextSize(.06);
    texbc1->SetTextColor(1);
    TLatex *texbc2 = new TLatex(-15,-15,"BC after : 3-#color[8]{8}, before : 11-17");
    texbc2->SetTextSize(.06);
    texbc2->SetTextColor(1);
    TLatex *texbc3 = new TLatex(-15,-15,"BC after : 3-#color[2]{9}, before : 11-17");
    texbc3->SetTextSize(.06);
    texbc3->SetTextColor(1);
    
    
    TLatex *texv01 = new TLatex(-15,-10,"#color[8]{V0A + V0C}");
    texv01->SetTextSize(.06);
    texv01->SetTextColor(1);
    TLatex *texv02 = new TLatex(-15,-10,"#color[2]{V0A} only");
    texv02->SetTextSize(.06);
    texv02->SetTextColor(1);
    TLatex *texv03 = new TLatex(-15,-10,"#color[4]{V0C} only");
    texv03->SetTextSize(.06);
    texv03->SetTextColor(1);
    
    TLatex *texf1 = new TLatex(-15,-5,"#color[8]{BB + BG} flag");
    texf1->SetTextSize(.06);
    texf1->SetTextColor(1);
    TLatex *texf2 = new TLatex(-15,-5,"#color[2]{BB} flag only");
    texf2->SetTextSize(.06);
    texf2->SetTextColor(1);
    TLatex *texf3 = new TLatex(-15,-5,"#color[4]{BG} flag only");
    texf3->SetTextSize(.06);
    texf3->SetTextColor(1);
    
    /* I think that we don't need this part, now. 2015.08.18.(blim)
    TCanvas *cOutCanvas1[3][3][3];
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                cOutCanvas1[i][j][k]= new TCanvas(Form("cOutCanvas1%d_%d_%d",i,j,k),Form("cOutCanvas1%d_%d_%d",i,j,k),1000,500);
                cOutCanvas1[i][j][k]->Divide(3,2);
                cOutCanvas1[i][j][k]->cd(1);
                hNumTrkVsClsSPID[i][j][k]->Draw("colz");
                cOutCanvas1[i][j][k]->cd(2);
                hDenomTrkVsClsSPID[i][j][k]->Draw("colz");
                cOutCanvas1[i][j][k]->cd(3);
                h2DSPDEfficiency[i][j][k]->Draw("colz");
                 cOutCanvas1[i][j][k]->cd(4);
                hNumV0[i][j][k]->Draw("colz");
                if(i==0)texbc1->Draw("same");
                if(i==1)texbc2->Draw("same");
                if(i==2)texbc3->Draw("same");
                if(j==0)texv01->Draw("same");
                if(j==1)texv02->Draw("same");
                if(j==2)texv03->Draw("same");
                if(k==0)texf1->Draw("same");
                if(k==1)texf2->Draw("same");
                if(k==2)texf3->Draw("same");
                texNum->Draw("same");
                cOutCanvas1[i][j][k]->cd(5);
                hDenomV0[i][j][k]->Draw("colz");
                if(i==0)texbc1->Draw("same");
                if(i==1)texbc2->Draw("same");
                if(i==2)texbc3->Draw("same");
                if(j==0)texv01->Draw("same");
                if(j==1)texv02->Draw("same");
                if(j==2)texv03->Draw("same");
                if(k==0)texf1->Draw("same");
                if(k==1)texf2->Draw("same");
                if(k==2)texf3->Draw("same");
                texDe->Draw("same");
                cOutCanvas1[i][j][k]->cd(6);
                h2DV0Efficiency[i][j][k]->Draw("colz");
                if(i==0)texbc1->Draw("same");
                if(i==1)texbc2->Draw("same");
                if(i==2)texbc3->Draw("same");
                if(j==0)texv01->Draw("same");
                if(j==1)texv02->Draw("same");
                if(j==2)texv03->Draw("same");
                if(k==0)texf1->Draw("same");
                if(k==1)texf2->Draw("same");
                if(k==2)texf3->Draw("same");
                texR->Draw("same");

//                cOutCanvas1[i][j][k]->SaveAs(Form("./Fig/cOutCanvas1%d_%d_%d.pdf",i,j,k));


            }
        }
    }*/
    


    
    TLatex *texbc11 = new TLatex(3,0.3,"BC after : 3-#color[4]{7}, before : 11-17");
    texbc11->SetTextSize(.06);
    texbc11->SetTextColor(1);
    TLatex *texbc21 = new TLatex(3,0.4,"BC after : 3-#color[2]{8}, before : 11-17");
    texbc21->SetTextSize(.06);
    texbc21->SetTextColor(1);
    TLatex *texbc31 = new TLatex(3,0.5,"BC after : 3-#color[1]{9}, before : 11-17");
    texbc31->SetTextSize(.06);
    texbc31->SetTextColor(1);
    TLatex *texf21 = new TLatex(3,0.1,"#color[1]{BB} flag only");
    texf21->SetTextSize(.06);
    texf21->SetTextColor(1);
    TLatex *texv011 = new TLatex(3,0.2,"#color[1]{V0A + V0C}");
    texv011->SetTextSize(.06);
    texv011->SetTextColor(1);
    
    TCanvas *cOutCanvasBC= new TCanvas("cOutCanvasBC","cOutCanvasBC",1000,500);
    cOutCanvasBC->Divide(4,1);
    cOutCanvasBC->cd(1);
    hPurityBC[2][0][1]->SetLineColor(1);
    hPurityBC[2][0][1]->Draw("same");
    //hPurityBC[1][0][1]->SetLineColor(2);
    //hPurityBC[1][0][1]->Draw("same");
    //hPurityBC[0][0][1]->SetLineColor(4);
    //hPurityBC[0][0][1]->Draw("same");
    texbc11->Draw("same");
    texbc21->Draw("same");
    texbc31->Draw("same");
    texf21->Draw("same");
    texv011->Draw("same");
    cOutCanvasBC->cd(2);
    hEfficiencyBC[2][0][1]->SetLineColor(1);
    hEfficiencyBC[2][0][1]->Draw("same");
    //hEfficiencyBC[1][0][1]->SetLineColor(2);
    //hEfficiencyBC[1][0][1]->Draw("same");
    //hEfficiencyBC[0][0][1]->SetLineColor(4);
    //hEfficiencyBC[0][0][1]->Draw("same");
    texbc11->Draw("same");
    texbc21->Draw("same");
    texbc31->Draw("same");
    texf21->Draw("same");
    texv011->Draw("same");
    cOutCanvasBC->cd(3);
    hRejectEfficiencyBC[2][0][1]->SetLineColor(1);
    hRejectEfficiencyBC[2][0][1]->Draw("same");
    //hRejectEfficiencyBC[1][0][1]->SetLineColor(2);
    //hRejectEfficiencyBC[1][0][1]->Draw("same");
    //hRejectEfficiencyBC[0][0][1]->SetLineColor(4);
    //hRejectEfficiencyBC[0][0][1]->Draw("same");
    texbc11->Draw("same");
    texbc21->Draw("same");
    texbc31->Draw("same");
    texf21->Draw("same");
    texv011->Draw("same");
    cOutCanvasBC->cd(4);
    hSPDEfficiencyBC[2][0][1]->SetLineColor(1);
    hSPDEfficiencyBC[2][0][1]->Draw("same");
    //hSPDEfficiencyBC[1][0][1]->SetLineColor(2);
    //hSPDEfficiencyBC[1][0][1]->Draw("same");
    //hSPDEfficiencyBC[0][0][1]->SetLineColor(4);
    //hSPDEfficiencyBC[0][0][1]->Draw("same");
    texbc11->Draw("same");
    texbc21->Draw("same");
    texbc31->Draw("same");
    texf21->Draw("same");
    texv011->Draw("same");

    TLatex *texbc32 = new TLatex(3,0.1,"BC after : 3-#color[1]{9}, before : 11-17");
    texbc32->SetTextSize(.06);
    texbc32->SetTextColor(1);
    TLatex *texf22 = new TLatex(3,0.5,"#color[1]{BB} flag only");
    texf22->SetTextSize(.06);
    texf22->SetTextColor(1);
    TLatex *texf12 = new TLatex(3,0.4,"#color[2]{BB + BG} flag");
    texf12->SetTextSize(.06);
    texf12->SetTextColor(1);
    TLatex *texf32 = new TLatex(3,0.3,"#color[4]{BG} flag only");
    texf32->SetTextSize(.06);
    texf32->SetTextColor(1);
    TLatex *texv012 = new TLatex(3,0.2,"#color[1]{V0A + V0C}");
    texv012->SetTextSize(.06);
    texv012->SetTextColor(1);
    
    TCanvas *cOutCanvasFlag= new TCanvas("cOutCanvasFlag","cOutCanvasFlag",1000,500);
    cOutCanvasFlag->Divide(4,1);
    cOutCanvasFlag->cd(1);
    hPurityBC[2][0][1]->SetLineColor(1);
    hPurityBC[2][0][1]->Draw("same");
    //hPurityBC[2][0][0]->SetLineColor(2);
    //hPurityBC[2][0][0]->Draw("same");
    //hPurityBC[2][0][2]->SetLineColor(4);
    //hPurityBC[2][0][2]->Draw("same");
    texbc32->Draw("same");
    texf22->Draw("same");
    texf12->Draw("same");
    texf32->Draw("same");
    texv012->Draw("same");
    cOutCanvasFlag->cd(2);
    hEfficiencyBC[2][0][1]->SetLineColor(1);
    hEfficiencyBC[2][0][1]->Draw("same");
    //hEfficiencyBC[2][0][0]->SetLineColor(2);
    //hEfficiencyBC[2][0][0]->Draw("same");
    //hEfficiencyBC[2][0][2]->SetLineColor(4);
    //hEfficiencyBC[2][0][2]->Draw("same");
    texbc32->Draw("same");
    texf22->Draw("same");
    texf12->Draw("same");
    texf32->Draw("same");
    texv012->Draw("same");
    cOutCanvasFlag->cd(3);
    hRejectEfficiencyBC[2][0][1]->SetLineColor(1);
    hRejectEfficiencyBC[2][0][1]->Draw("same");
    //hRejectEfficiencyBC[2][0][0]->SetLineColor(2);
    //hRejectEfficiencyBC[2][0][0]->Draw("same");
    //hRejectEfficiencyBC[2][0][2]->SetLineColor(4);
    //hRejectEfficiencyBC[2][0][2]->Draw("same");
    texbc32->Draw("same");
    texf22->Draw("same");
    texf12->Draw("same");
    texf32->Draw("same");
    texv012->Draw("same");
    cOutCanvasFlag->cd(4);
    hSPDEfficiencyBC[2][0][1]->SetLineColor(1);
    hSPDEfficiencyBC[2][0][1]->Draw("same");
    //hSPDEfficiencyBC[2][0][0]->SetLineColor(2);
    //hSPDEfficiencyBC[2][0][0]->Draw("same");
    //hSPDEfficiencyBC[2][0][2]->SetLineColor(4);
    //hSPDEfficiencyBC[2][0][2]->Draw("same");
    texbc32->Draw("same");
    texf22->Draw("same");
    texf12->Draw("same");
    texf32->Draw("same");
    texv012->Draw("same");
    

    TLatex *texbc33 = new TLatex(3,0.1,"BC after : 3-#color[1]{9}, before : 11-17");
    texbc33->SetTextSize(.06);
    texbc33->SetTextColor(1);
    TLatex *texf23 = new TLatex(3,0.2,"#color[1]{BB} flag only");
    texf23->SetTextSize(.06);
    texf23->SetTextColor(1);
    TLatex *texv013 = new TLatex(3,0.5,"#color[1]{V0A + V0C}");
    texv013->SetTextSize(.06);
    texv013->SetTextColor(1);
    TLatex *texv023 = new TLatex(3,0.4,"#color[2]{V0A} only");
    texv023->SetTextSize(.06);
    texv023->SetTextColor(1);
    TLatex *texv033 = new TLatex(3,0.3,"#color[4]{V0C} only");
    texv033->SetTextSize(.06);
    texv033->SetTextColor(1);
    
    TCanvas *cOutCanvasV0= new TCanvas("cOutCanvasV0","cOutCanvasV0",1000,500);
    cOutCanvasV0->Divide(4,1);
    cOutCanvasV0->cd(1);
    hPurityBC[2][0][1]->SetLineColor(1);
    hPurityBC[2][0][1]->Draw("same");
    //hPurityBC[2][1][1]->SetLineColor(2);
    //hPurityBC[2][1][1]->Draw("same");
    //hPurityBC[2][2][1]->SetLineColor(4);
    //hPurityBC[2][2][1]->Draw("same");
    texbc33->Draw("same");
    texf23->Draw("same");
    texv013->Draw("same");
    texv023->Draw("same");
    texv033->Draw("same");
    cOutCanvasV0->cd(2);
    hEfficiencyBC[2][0][1]->SetLineColor(1);
    hEfficiencyBC[2][0][1]->Draw("same");
    //hEfficiencyBC[2][1][1]->SetLineColor(2);
    //hEfficiencyBC[2][1][1]->Draw("same");
    //hEfficiencyBC[2][2][1]->SetLineColor(4);
    //hEfficiencyBC[2][2][1]->Draw("same");
    texbc33->Draw("same");
    texf23->Draw("same");
    texv013->Draw("same");
    texv023->Draw("same");
    texv033->Draw("same");    cOutCanvasV0->cd(3);
    hRejectEfficiencyBC[2][0][1]->SetLineColor(1);
    hRejectEfficiencyBC[2][0][1]->Draw("same");
    //hRejectEfficiencyBC[2][1][1]->SetLineColor(2);
    //hRejectEfficiencyBC[2][1][1]->Draw("same");
    //hRejectEfficiencyBC[2][2][1]->SetLineColor(4);
    //hRejectEfficiencyBC[2][2][1]->Draw("same");
    texbc33->Draw("same");
    texf23->Draw("same");
    texv013->Draw("same");
    texv023->Draw("same");
    texv033->Draw("same");
    cOutCanvasV0->cd(4);
    hSPDEfficiencyBC[2][0][1]->SetLineColor(1);
    hSPDEfficiencyBC[2][0][1]->Draw("same");
    //hSPDEfficiencyBC[2][1][1]->SetLineColor(2);
    //hSPDEfficiencyBC[2][1][1]->Draw("same");
    //hSPDEfficiencyBC[2][2][1]->SetLineColor(4);
    //hSPDEfficiencyBC[2][2][1]->Draw("same");
    texbc33->Draw("same");
    texf23->Draw("same");
    texv013->Draw("same");
    texv023->Draw("same");
    texv033->Draw("same");


    //// here star to AD information ////
 

    TLatex *ADtexv011 = new TLatex(3,0.2,"#color[1]{ADA + ADC}");
    ADtexv011->SetTextSize(.06);
    ADtexv011->SetTextColor(1);
    
    TCanvas *cADOutCanvasBC= new TCanvas("cADOutCanvasBC","cADOutCanvasBC",750,500);
    cADOutCanvasBC->Divide(3,1);
    cADOutCanvasBC->cd(1);
    hADPurityBC[2][0][1]->SetLineColor(1);
    hADPurityBC[2][0][1]->SetMaximum(1.);

    hADPurityBC[2][0][1]->Draw("same");
    //hADPurityBC[1][0][1]->SetLineColor(2);
    //hADPurityBC[1][0][1]->Draw("same");
    //hADPurityBC[0][0][1]->SetLineColor(4);
    //hADPurityBC[0][0][1]->Draw("same");
    texbc11->Draw("same");
    texbc21->Draw("same");
    texbc31->Draw("same");
    texf21->Draw("same");
    ADtexv011->Draw("same");
    cADOutCanvasBC->cd(2);
    hADEfficiencyBC[2][0][1]->SetLineColor(1);
    hADEfficiencyBC[2][0][1]->SetMaximum(1.);

    hADEfficiencyBC[2][0][1]->Draw("same");
    //hADEfficiencyBC[1][0][1]->SetLineColor(2);
    //hADEfficiencyBC[1][0][1]->Draw("same");
    //hADEfficiencyBC[0][0][1]->SetLineColor(4);
    //hADEfficiencyBC[0][0][1]->Draw("same");
    texbc11->Draw("same");
    texbc21->Draw("same");
    texbc31->Draw("same");
    texf21->Draw("same");
    ADtexv011->Draw("same");
    cADOutCanvasBC->cd(3);
    hADRejectEfficiencyBC[2][0][1]->SetLineColor(1);
    hADRejectEfficiencyBC[2][0][1]->SetMaximum(1.);

    hADRejectEfficiencyBC[2][0][1]->Draw("same");
    //hADRejectEfficiencyBC[1][0][1]->SetLineColor(2);
    //hADRejectEfficiencyBC[1][0][1]->Draw("same");
    //hADRejectEfficiencyBC[0][0][1]->SetLineColor(4);
    //hADRejectEfficiencyBC[0][0][1]->Draw("same");
    texbc11->Draw("same");
    texbc21->Draw("same");
    texbc31->Draw("same");
    texf21->Draw("same");
    ADtexv011->Draw("same");
   
    
   
    TLatex *ADtexv012 = new TLatex(3,0.2,"#color[1]{ADA + ADC}");
    ADtexv012->SetTextSize(.06);
    ADtexv012->SetTextColor(1);
    
    TCanvas *cADOutCanvasFlag= new TCanvas("cADOutCanvasFlag","cADOutCanvasFlag",750,500);
    cADOutCanvasFlag->Divide(3,1);
    cADOutCanvasFlag->cd(1);
    hADPurityBC[2][0][1]->SetLineColor(1);
    hADPurityBC[2][0][1]->SetMaximum(1.);
    hADPurityBC[2][0][1]->Draw("same");
    //hADPurityBC[2][0][0]->SetLineColor(2);
    //hADPurityBC[2][0][0]->Draw("same");
    //hADPurityBC[2][0][2]->SetLineColor(4);
    //hADPurityBC[2][0][2]->Draw("same");
    texbc32->Draw("same");
    texf22->Draw("same");
    texf12->Draw("same");
    texf32->Draw("same");
    ADtexv012->Draw("same");
    cADOutCanvasFlag->cd(2);
    hADEfficiencyBC[2][0][1]->SetLineColor(1);
    hADEfficiencyBC[2][0][1]->SetMaximum(1.);

    hADEfficiencyBC[2][0][1]->Draw("same");
    //hADEfficiencyBC[2][0][0]->SetLineColor(2);
    //hADEfficiencyBC[2][0][0]->Draw("same");
    //hADEfficiencyBC[2][0][2]->SetLineColor(4);
    //hADEfficiencyBC[2][0][2]->Draw("same");
    texbc32->Draw("same");
    texf22->Draw("same");
    texf12->Draw("same");
    texf32->Draw("same");
    ADtexv012->Draw("same");
    cADOutCanvasFlag->cd(3);
    hADRejectEfficiencyBC[2][0][1]->SetLineColor(1);
    hADRejectEfficiencyBC[2][0][1]->SetMaximum(1.);

    hADRejectEfficiencyBC[2][0][1]->Draw("same");
    //hADRejectEfficiencyBC[2][0][0]->SetLineColor(2);
    //hADRejectEfficiencyBC[2][0][0]->Draw("same");
    //hADRejectEfficiencyBC[2][0][2]->SetLineColor(4);
    //hADRejectEfficiencyBC[2][0][2]->Draw("same");
    texbc32->Draw("same");
    texf22->Draw("same");
    texf12->Draw("same");
    texf32->Draw("same");
    ADtexv012->Draw("same");
   
    
    
    TLatex *ADtexv013 = new TLatex(3,0.5,"#color[1]{ADA + ADC}");
    ADtexv013->SetTextSize(.06);
    ADtexv013->SetTextColor(1);
    TLatex *ADtexv023 = new TLatex(3,0.4,"#color[2]{ADA} only");
    ADtexv023->SetTextSize(.06);
    ADtexv023->SetTextColor(1);
    TLatex *ADtexv033 = new TLatex(3,0.3,"#color[4]{ADC} only");
    ADtexv033->SetTextSize(.06);
    ADtexv033->SetTextColor(1);
    
    TCanvas *cADOutCanvasV0= new TCanvas("cADOutCanvasV0","cADOutCanvasV0",750,500);
    cADOutCanvasV0->Divide(3,1);
    cADOutCanvasV0->cd(1);
    hADPurityBC[2][0][1]->SetLineColor(1);
    hADPurityBC[2][0][1]->SetMaximum(1.);

    hADPurityBC[2][0][1]->Draw("same");
    //hADPurityBC[2][1][1]->SetLineColor(2);
    //hADPurityBC[2][1][1]->Draw("same");
    //hADPurityBC[2][2][1]->SetLineColor(4);
    //hADPurityBC[2][2][1]->Draw("same");
    texbc33->Draw("same");
    texf23->Draw("same");
    ADtexv013->Draw("same");
    ADtexv023->Draw("same");
    ADtexv033->Draw("same");
    cADOutCanvasV0->cd(2);
    hADEfficiencyBC[2][0][1]->SetLineColor(1);
    hADEfficiencyBC[2][0][1]->SetMaximum(1.);

    hADEfficiencyBC[2][0][1]->Draw("same");
    //hADEfficiencyBC[2][1][1]->SetLineColor(2);
    //hADEfficiencyBC[2][1][1]->Draw("same");
    //hADEfficiencyBC[2][2][1]->SetLineColor(4);
    //hADEfficiencyBC[2][2][1]->Draw("same");
    texbc33->Draw("same");
    texf23->Draw("same");
    ADtexv013->Draw("same");
    ADtexv023->Draw("same");
    ADtexv033->Draw("same");
    cADOutCanvasV0->cd(3);
    hADRejectEfficiencyBC[2][0][1]->SetLineColor(1);
    hADRejectEfficiencyBC[2][0][1]->SetMaximum(1.);
    hADRejectEfficiencyBC[2][0][1]->Draw("same");
    //hADRejectEfficiencyBC[2][1][1]->SetLineColor(2);
    //hADRejectEfficiencyBC[2][1][1]->Draw("same");
    //hADRejectEfficiencyBC[2][2][1]->SetLineColor(4);
    //hADRejectEfficiencyBC[2][2][1]->Draw("same");
    texbc33->Draw("same");
    texf23->Draw("same");
    ADtexv013->Draw("same");
    ADtexv023->Draw("same");
    ADtexv033->Draw("same");
    cout << "end" << endl;

    hTrkVsClsSPID_Slope_MB01= new TCanvas("hTrkVsClsSPID_Slope_MB01","hTrkVsClsSPID_Slope_MB01",800,500);
    hTrkVsClsSPIDSlopeM->Draw("colz");
    hTrkVsClsSPID_Slope_MB01->SetLogz();
    hTrkVsClsSPID_Slope_MB01->SaveAs("./Fig/hTrkVsClsSPID_Slope_MB01.pdf");

    hTrkVsClsSPID_Slope_MB02= new TCanvas("hTrkVsClsSPID_Slope_MB02","hTrkVsClsSPID_Slope_MB02",800,500);
    hTrkVsClsSPIDSlopeM2->Draw("colz");
    hTrkVsClsSPID_Slope_MB02->SetLogz();
    hTrkVsClsSPID_Slope_MB02->SaveAs("./Fig/hTrkVsClsSPID_Slope_MB02.pdf");
    
    hTrkVsClsSPID_Slope_HM01= new TCanvas("hTrkVsClsSPID_Slope_HM01","hTrkVsClsSPID_Slope_HM01",800,500);
    hTrkVsClsSPIDSlopeM_HM->Draw("colz");
    hTrkVsClsSPID_Slope_HM01->SetLogz();
    hTrkVsClsSPID_Slope_HM01->SaveAs("./Fig/hTrkVsClsSPID_Slope_HM01.pdf");
    
    hTrkVsClsSPID_Slope_HM02= new TCanvas("hTrkVsClsSPID_Slope_HM02","hTrkVsClsSPID_Slope_HM02",800,500);
    hTrkVsClsSPIDSlopeM_HM2->Draw("colz");
    hTrkVsClsSPID_Slope_HM02->SetLogz();
    hTrkVsClsSPID_Slope_HM02->SaveAs("./Fig/hTrkVsClsSPID_Slope_HM02.pdf");
    
    TFile *fileout = TFile::Open(ss+"/BGMonitorHisto.root","RECREATE");
    cOutCanvasBC->Write("cOutCanvasBC");
    cOutCanvasFlag->Write("cOutCanvasFlag");
    cOutCanvasV0->Write("cOutCanvasV0");
    cADOutCanvasBC->Write("cADOutCanvasBC");
    cADOutCanvasFlag->Write("cADOutCanvasFlag");
    cADOutCanvasV0->Write("cADOutCanvasV0");
//    hTotalTrkVsClsSPID->Write("hTotalTrkVsClsSPID");
//    hTotalTrkVsClsSPID->SaveAs("./Fig/hTotalTrkVsClsSPID.pdf");
    //---------------------
    /*
    cOutCanvasBC->SaveAs("./Fig/cOutCanvasBC.pdf");
    cOutCanvasFlag->SaveAs("./Fig/cOutCanvasFlag.pdf");
    cOutCanvasV0->SaveAs("./Fig/cOutCanvasV0.pdf");
    cADOutCanvasBC->SaveAs("./Fig/cADOutCanvasBC.pdf");
    cADOutCanvasFlag->SaveAs("./Fig/cADOutCanvasFlag.pdf");
    cADOutCanvasV0->SaveAs("./Fig/cADOutCanvasV0.pdf");
    */
//    cOutCanvas_output->SaveAs("./Fig/cOutCanvas_output.pdf");  // it is on upper part
    //---------------------save output 2015.08.18.(blim)
    
    
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                Int_t check000 = i*100+j*10+k;
                for ( int l=0; l<8;l++){
                    if(check000==bunchinputarray[l]){
                hPurityBC[i][j][k]->Write();
                hEfficiencyBC[i][j][k]->Write();
                hRejectEfficiencyBC[i][j][k]->Write();
                hSPDEfficiencyBC[i][j][k]->Write();
                hNumTrkVsClsSPID[i][j][k]->Write();
                hDenomTrkVsClsSPID[i][j][k]->Write();
                h2DSPDEfficiency[i][j][k]->Write();
                hNumV0[i][j][k]->Write();
                hDenomV0[i][j][k]->Write();
                h2DV0Efficiency[i][j][k]->Write();
                    }
                }
            }
        }
    }
    fileout->Close();
 
}