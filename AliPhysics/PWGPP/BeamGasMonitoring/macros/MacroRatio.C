// this is macro to get "efficiency", "purity" and "rejection efficiency".
void MacroRatio(){
  for (Int_t irun=1; irun<2 ; irun++ ) {
    GetRatio(irun); // need to be modyfied. irun automatically read run number and produce output for each run
  }
}
GetRatio(Int_t split=-1)
{   
  TString ss = Form("%03d",split);
  cout << ss << endl;   
    
  TFile *filein = TFile::Open(ss+"/AnalysisResults.root","READ");
  if(filein==NULL) continue;
  cout << "open root file : " << ss << "/AnalysisResults.root"<< endl;

  TDirectoryFile *dir = (TDirectoryFile*)filein->Get("BeamGasMon");
    
  TList *list = (TList*)dir->Get("cOutputH");
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
                                                                 
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
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
	//  hNumV0[i][j][k]->GetXaxis()->SetTitle("V0A-V0C");
	//  hNumV0[i][j][k]->GetYaxis()->SetTitle("V0A+V0C");
        hDenomV0[i][j][k] = (TH2F*)list->FindObject(Form("hDenomV0%d_V0%d_Flag%d",i,j,k));
	//  hDenomV0[i][j][k]->GetXaxis()->SetTitle("V0A-V0C");
	//  hDenomV0[i][j][k]->GetYaxis()->SetTitle("V0A+V0C");
      }
    }
  } // read histogram from output

  TH1F *hPurityBC[3][3][3];
  TH1F *hEfficiencyBC[3][3][3];
  TH1F *hRejectEfficiencyBC[3][3][3];
  TH1F *hSPDEfficiencyBC[3][3][3];
  TH2F *h2DSPDEfficiency[3][3][3];
  TH2F *h2DV0Efficiency[3][3][3];

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
          
	hPurityBC[i][j][k] = new TH1F(Form("hPurityBC%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
	hEfficiencyBC[i][j][k] = new TH1F(Form("hEfficiencyBC%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
	hRejectEfficiencyBC[i][j][k] = new TH1F(Form("hRejectEfficiencyBC%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
	hSPDEfficiencyBC[i][j][k] = new TH1F(Form("hSPDEfficiencyBC%d_V0%d_Flag%d",i,j,k),"; Spd tracklet", 200, 0, 200);

	h2DSPDEfficiency[i][j][k] = new TH2F(Form("h2DSPDEfficiency%d_V0%d_Flag%d",i,j,k),"; Spd Cluster vs tracklet",140,0,140,1000,0,1000);
	h2DSPDEfficiency[i][j][k]->GetXaxis()->SetTitle("Tracklet");
	h2DSPDEfficiency[i][j][k]->GetYaxis()->SetTitle("Cluster (fspdC1)");
   
	h2DV0Efficiency[i][j][k] = new TH2F(Form("h2DV0Efficiency%d_V0%d_Flag%d",i,j,k),"; V0 timing information",1000,-20,30,1000,-20,30);
	h2DV0Efficiency[i][j][k]->GetXaxis()->SetTitle("V0A-V0C");
	h2DV0Efficiency[i][j][k]->GetYaxis()->SetTitle("V0A+V0C");
     
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
      }
    }
  }
  
  TLatex *texbc1 = new TLatex(0.4,0.1,"Bunch after : 3-#color[4]{7}, before : 11-17");
  texbc1->SetTextSize(.04);
  texbc1->SetTextColor(1);
  TLatex *texbc2 = new TLatex(0.4,0.1,"Bunch after : 3-#color[8]{8}, before : 11-17");
  texbc2->SetTextSize(.04);
  texbc2->SetTextColor(1);
  TLatex *texbc3 = new TLatex(0.4,0.1,"Bunch after : 3-#color[2]{9}, before : 11-17");
  texbc3->SetTextSize(.04);
  texbc3->SetTextColor(1);
     
  TLatex *texv01 = new TLatex(0.4,0.2,"#color[8]{V0A + V0C}");
  texv01->SetTextSize(.04);
  texv01->SetTextColor(1);
  TLatex *texv02 = new TLatex(0.4,0.2,"#color[2]{V0A} only");
  texv02->SetTextSize(.04);
  texv02->SetTextColor(1);
  TLatex *texv03 = new TLatex(0.4,0.2,"#color[4]{V0C} only");
  texv03->SetTextSize(.04);
  texv03->SetTextColor(1);
    
  TLatex *texf1 = new TLatex(0.4,0.3,"#color[8]{BB + BG} flag");
  texf1->SetTextSize(.04);
  texf1->SetTextColor(1);
  TLatex *texf2 = new TLatex(0.4,0.3,"#color[2]{BB} flag only");
  texf2->SetTextSize(.04);
  texf2->SetTextColor(1);
  TLatex *texf3 = new TLatex(0.4,0.3,"#color[4]{BG} flag only");
  texf3->SetTextSize(.04);
  texf3->SetTextColor(1);
    
  TCanvas *cOutCanvas[3][3][3];
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
	cOutCanvas[i][j][k]= new TCanvas(Form("cOutCanvas%d_%d_%d",i,j,k),Form("cOutCanvas%d_%d_%d",i,j,k),1000,500);
	cOutCanvas[i][j][k]->Divide(4,1);
	cOutCanvas[i][j][k]->cd(1);
	hPurityBC[i][j][k]->Draw("same");
	if(i==0)texbc1->Draw("same");
	if(i==1)texbc2->Draw("same");
	if(i==2)texbc3->Draw("same");
	if(j==0)texv01->Draw("same");
	if(j==1)texv02->Draw("same");
	if(j==2)texv03->Draw("same");
	if(k==0)texf1->Draw("same");
	if(k==1)texf1->Draw("same");
	if(k==2)texf1->Draw("same");
	cOutCanvas[i][j][k]->cd(2);
	hEfficiencyBC[i][j][k]->Draw("same");
	if(i==0)texbc1->Draw("same");
	if(i==1)texbc2->Draw("same");
	if(i==2)texbc3->Draw("same");
	if(j==0)texv01->Draw("same");
	if(j==1)texv02->Draw("same");
	if(j==2)texv03->Draw("same");
	if(k==0)texf1->Draw("same");
	if(k==1)texf1->Draw("same");
	if(k==2)texf1->Draw("same");

	cOutCanvas[i][j][k]->cd(3);
	hRejectEfficiencyBC[i][j][k]->Draw("same");
	if(i==0)texbc1->Draw("same");
	if(i==1)texbc2->Draw("same");
	if(i==2)texbc3->Draw("same");
	if(j==0)texv01->Draw("same");
	if(j==1)texv02->Draw("same");
	if(j==2)texv03->Draw("same");
	if(k==0)texf1->Draw("same");
	if(k==1)texf1->Draw("same");
	if(k==2)texf1->Draw("same");

	cOutCanvas[i][j][k]->cd(4);
	hSPDEfficiencyBC[i][j][k]->Draw("same");
	if(i==0)texbc1->Draw("same");
	if(i==1)texbc2->Draw("same");
	if(i==2)texbc3->Draw("same");
	if(j==0)texv01->Draw("same");
	if(j==1)texv02->Draw("same");
	if(j==2)texv03->Draw("same");
	if(k==0)texf1->Draw("same");
	if(k==1)texf1->Draw("same");
	if(k==2)texf1->Draw("same");

                
	cOutCanvas[i][j][k]->SaveAs(Form("./Fig/cOutCanvas%d_%d_%d.pdf",i,j,k));
      }
    }
  }
    
  TLatex *texNum = new TLatex(6,900,"Numerator (!bgID & Good Event)");
  texNum->SetTextSize(.04);
  texNum->SetTextColor(1);
  TLatex *texDe = new TLatex(6,900,"Denominator (!bgID)");
  texDe->SetTextSize(.04);
  texDe->SetTextColor(1);
  TLatex *texR = new TLatex(6,900,"Ratio");
  texR->SetTextSize(.04);
  texR->SetTextColor(1);
    
  TLatex *texbc1 = new TLatex(3,240,"Bunch after : 3-#color[4]{7}, before : 11-17");
  texbc1->SetTextSize(.04);
  texbc1->SetTextColor(1);
  TLatex *texbc2 = new TLatex(3,240,"Bunch after : 3-#color[8]{8}, before : 11-17");
  texbc2->SetTextSize(.04);
  texbc2->SetTextColor(1);
  TLatex *texbc3 = new TLatex(3,240,"Bunch after : 3-#color[2]{9}, before : 11-17");
  texbc3->SetTextSize(.04);
  texbc3->SetTextColor(1);
    
    
  TLatex *texv01 = new TLatex(3,80,"#color[8]{V0A + V0C}");
  texv01->SetTextSize(.04);
  texv01->SetTextColor(1);
  TLatex *texv02 = new TLatex(3,80,"#color[2]{V0A} only");
  texv02->SetTextSize(.04);
  texv02->SetTextColor(1);
  TLatex *texv03 = new TLatex(3,80,"#color[4]{V0C} only");
  texv03->SetTextSize(.04);
  texv03->SetTextColor(1);
    
  TLatex *texf1 = new TLatex(3,160,"#color[8]{BB + BG} flag");
  texf1->SetTextSize(.04);
  texf1->SetTextColor(1);
  TLatex *texf2 = new TLatex(3,160,"#color[2]{BB} flag only");
  texf2->SetTextSize(.04);
  texf2->SetTextColor(1);
  TLatex *texf3 = new TLatex(3,160,"#color[4]{BG} flag only");
  texf3->SetTextSize(.04);
  texf3->SetTextColor(1);
    
  TCanvas *cOutCanvas1[3][3][3];
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
	cOutCanvas1[i][j][k]= new TCanvas(Form("cOutCanvas1%d_%d_%d",i,j,k),Form("cOutCanvas1%d_%d_%d",i,j,k),1000,500);
	cOutCanvas1[i][j][k]->Divide(3,2);
	cOutCanvas1[i][j][k]->cd(1);
	hNumTrkVsClsSPID[i][j][k]->Draw("colz");
	if(i==0)texbc1->Draw("same");
	if(i==1)texbc2->Draw("same");
	if(i==2)texbc3->Draw("same");
	if(j==0)texv01->Draw("same");
	if(j==1)texv02->Draw("same");
	if(j==2)texv03->Draw("same");
	if(k==0)texf1->Draw("same");
	if(k==1)texf1->Draw("same");
	if(k==2)texf1->Draw("same");
	texNum->Draw("same");
	cOutCanvas1[i][j][k]->cd(2);
                
	hDenomTrkVsClsSPID[i][j][k]->Draw("colz");
	if(i==1)texbc2->Draw("same");
	if(i==2)texbc3->Draw("same");
	if(j==0)texv01->Draw("same");
	if(j==1)texv02->Draw("same");
	if(j==2)texv03->Draw("same");
	if(k==0)texf1->Draw("same");
	if(k==1)texf1->Draw("same");
	if(k==2)texf1->Draw("same");
	texDe->Draw("same");
	cOutCanvas1[i][j][k]->cd(3);

	h2DSPDEfficiency[i][j][k]->Draw("colz");
	if(i==1)texbc2->Draw("same");
	if(i==2)texbc3->Draw("same");
	if(j==0)texv01->Draw("same");
	if(j==1)texv02->Draw("same");
	if(j==2)texv03->Draw("same");
	if(k==0)texf1->Draw("same");
	if(k==1)texf1->Draw("same");
	if(k==2)texf1->Draw("same");
	texR->Draw("same");
	cOutCanvas1[i][j][k]->cd(4);

	hNumV0[i][j][k]->Draw("colz");
	cOutCanvas1[i][j][k]->cd(5);

	hDenomV0[i][j][k]->Draw("colz");
	cOutCanvas1[i][j][k]->cd(6);

	h2DV0Efficiency[i][j][k]->Draw("colz");
	cOutCanvas1[i][j][k]->SaveAs(Form("./Fig/cOutCanvas1%d_%d_%d.pdf",i,j,k));


      }
    }
  }
    
    
  TFile *fileout = TFile::Open(ss+"/BGMonitorHisto.root","RECREATE");
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
                
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
  fileout->Close();
 
}
