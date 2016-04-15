#include <stdio.h>

void ValueFromRoot(int i, const char *filename);
#define buf 128
Double_t ResultValue[500][7]; // [n][0]=RunNumber, [n][1]=Purity of MB [n][2]=Efficiency of MB, [n][3]=Rejection Eff. of MB, [n][4]=Purity of HM and so on..

void RunByRunAnal(){
    const char *line[500] = new char*;
    char str[buf];

    FILE *fp;
    char str[buf];
    char* line_p;

    fp = fopen("./inputlist.txt","r");

    cout << "check the inputlist file" << endl;
    cout << "________________________" << endl;

    int ii=0;
    while(fgets(str,buf,fp)){
        if((line_p = strchr(str, '\n')) != NULL)*line_p ='\0';
        //memcpy(line[ii], str, strlen(str));
        cout << ii+1 << " file is loaded at: " << str << endl;
        ValueFromRoot(ii,str);
        cout << "________________________" << endl;
        ii++;
    }
    fclose(fp);
 
    cout << "________________________" << endl;

    int j=0;
    
    //Test
    while(1){
        if(ResultValue[j][1]==0) break;
        if(j==3) break;
        cout << j << endl;
        cout << "RunNumber       :" << ResultValue[j][0] << endl;
        cout << "Purity(MB)      :" << ResultValue[j][1] << endl;
        cout << "Efficiency(MB)  :" << ResultValue[j][2] << endl;
        cout << "Reject Eff.(MB) :" << ResultValue[j][3] << endl;
        cout << "Purity(HM)      :" << ResultValue[j][4] << endl;
        cout << "Efficiency(HM)  :" << ResultValue[j][5] << endl;
        cout << "Reject Eff.(HM) :" << ResultValue[j][6] << endl;
        cout << endl;
        j++;
    }
    
    double noffiles = ii;
    
    TH1F *hRunbyRunPurity;
    TH1F *hRunbyRunPurity_HM;
    TH1F *hRunbyRunEff;
    TH1F *hRunbyRunEff_HM;
    TH1F *hRunbyRunRejectEff;
    TH1F *hRunbyRunRejectEff_HM;
    
    TH1F *hRunbyRunPurity = new TH1F("hRunbyRunPurity","Purity", noffiles, 0, noffiles);
    TH1F *hRunbyRunPurity_HM = new TH1F("hRunbyRunPurity_HM","Purity", noffiles, 0, noffiles);
    TH1F *hRunbyRunEff = new TH1F("hRunbyRunEff","Efficiency", noffiles, 0, noffiles);
    TH1F *hRunbyRunEff_HM = new TH1F("hRunbyRunEff_HM","Efficiency", noffiles, 0, noffiles);
    TH1F *hRunbyRunRejectEff = new TH1F("hRunbyRunRejectEff","Rejection Eff.", noffiles, 0, noffiles);
    TH1F *hRunbyRunRejectEff_HM = new TH1F("hRunbyRunRejectEff_HM","Rejection Eff.", noffiles, 0, noffiles);
    
    int k = 0;
    TString RunN;
    while(1){
        if(ResultValue[k][1]==0) break;
        hRunbyRunPurity->SetBinContent(k+1,ResultValue[k][1]);
        RunN = Form("%.0f",ResultValue[k][0]);
//        cout << RunN << endl;
        hRunbyRunPurity->GetXaxis()->SetBinLabel(k+1, RunN );
        
        hRunbyRunEff->SetBinContent(k+1,ResultValue[k][2]);
        hRunbyRunEff->GetXaxis()->SetBinLabel(k+1, RunN );
        
        hRunbyRunRejectEff->SetBinContent(k+1,ResultValue[k][3]);
        hRunbyRunRejectEff->GetXaxis()->SetBinLabel(k+1, RunN );
        
        hRunbyRunPurity_HM->SetBinContent(k+1,ResultValue[k][4]);
        hRunbyRunPurity_HM->GetXaxis()->SetBinLabel(k+1, RunN );
        
        hRunbyRunEff_HM->SetBinContent(k+1,ResultValue[k][5]);
        hRunbyRunEff_HM->GetXaxis()->SetBinLabel(k+1, RunN );
        
        hRunbyRunRejectEff_HM->SetBinContent(k+1,ResultValue[k][6]);
        hRunbyRunRejectEff_HM->GetXaxis()->SetBinLabel(k+1, RunN );
        
        k++;
    }
    
    hRunbyRunPurity->GetYaxis()->SetRangeUser(0,1.0);
    hRunbyRunEff->GetYaxis()->SetRangeUser(0,1.0);
    hRunbyRunRejectEff->GetYaxis()->SetRangeUser(0,1.0);
    hRunbyRunPurity->Draw();
    hRunbyRunEff->Draw();
    hRunbyRunRejectEff->Draw();
    hRunbyRunPurity_HM->GetYaxis()->SetRangeUser(0,1.0);
    hRunbyRunEff_HM->GetYaxis()->SetRangeUser(0,1.0);
    hRunbyRunRejectEff_HM->GetYaxis()->SetRangeUser(0,1.0);
    hRunbyRunPurity_HM->Draw();
    hRunbyRunEff_HM->Draw();
    hRunbyRunRejectEff_HM->Draw();
    
    TLatex *texbc1 = new TLatex(1 , 0.1,"#color[4]{Blue dot line} : Results from HM trigger");
    texbc1->SetTextSize(.04);
    texbc1->SetTextColor(1);
    TLatex *texbc2 = new TLatex(1 , 0.05,"#color[2]{Red line} : Results from MB trigger");
    texbc2->SetTextSize(.04);
    texbc2->SetTextColor(1);

    TLatex *texbc3 = new TLatex(1 , 0.95,"#color[4]{Blue dot line} : Results from HM trigger");
    texbc3->SetTextSize(.04);
    texbc3->SetTextColor(1);
    TLatex *texbc4 = new TLatex(1 , 0.9,"#color[2]{Red line} : Results from MB trigger");
    texbc4->SetTextSize(.04);
    texbc4->SetTextColor(1);
    
    
    RunByRunResult= new TCanvas("RunByRunResult","RunByRunResult",1000,500);
    gStyle->SetOptStat(0);
    RunByRunResult->Divide(3,1);
    
    RunByRunResult->cd(1);
    
    hRunbyRunPurity->SetLineColor(2); //red is from MB
    hRunbyRunPurity->SetLineWidth(2);
    hRunbyRunPurity_HM->SetLineWidth(2);
    hRunbyRunPurity_HM->SetLineColor(4); //blue is from HM
    hRunbyRunPurity_HM->SetLineStyle(2);
    hRunbyRunPurity->Draw("same");
    hRunbyRunPurity_HM->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    
    RunByRunResult->cd(2);
    
    hRunbyRunEff->SetLineColor(2);
    hRunbyRunEff->SetLineWidth(2);
    hRunbyRunEff_HM->SetLineWidth(2);
    hRunbyRunEff_HM->SetLineColor(4);
    hRunbyRunEff_HM->SetLineStyle(2);
    hRunbyRunEff->Draw("same");
    hRunbyRunEff_HM->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    
    RunByRunResult->cd(3);
    hRunbyRunRejectEff->SetLineColor(2);
    hRunbyRunRejectEff->SetLineWidth(2);
    hRunbyRunRejectEff_HM->SetLineWidth(2);
    hRunbyRunRejectEff_HM->SetLineColor(4);
    hRunbyRunRejectEff_HM->SetLineStyle(2);
    hRunbyRunRejectEff->Draw("same");
    hRunbyRunRejectEff_HM->Draw("same");
    texbc3->Draw("same");
    texbc4->Draw("same");
    
    RunByRunResult->SaveAs("./Fig/RunByRunResult.pdf");
    
}

void ValueFromRoot(int inputnum, const char *filename){

    int i = inputnum;
    TString fileroute = filename;
    TFile *filein = TFile::Open(fileroute+"/AnalysisResults.root","READ");
    if(filein==NULL) continue;
    cout << "open root file : " << "AnalysisResults.root"<< endl;
    
    TDirectoryFile *dir = (TDirectoryFile*)filein->Get("BeamGasMon");
    
    TList *list = (TList*)dir->Get("cOutputH_MB");
    TList *list2 = (TList*)dir->Get("cOutputH_HM");
    TTree *TTree = (TTree*)dir->Get("TreeTrack");

    Int_t Temp = 0;
    TTree->SetBranchAddress("runNumber",&Temp);
    TTree->GetEntry(0);
    
    
    //____________________MB________________________
    TH1F *hNumEffPurityBC[3][3][3];
    TH1F *hDenomEffBC[3][3][3];
    TH1F *hDenomPurityBC[3][3][3];
    TH1F *hDenomRejecEffBC[3][3][3];
    TH1F *hNumRejecEffBC[3][3][3];
    
    TH1F *hPurityBC[3][3][3];
    TH1F *hEfficiencyBC[3][3][3];
    TH1F *hRejectEfficiencyBC[3][3][3];
    
    TH1F *runNumber_hist;
    runNumber_hist = (TH1F*)list->FindObject("runNumber_hist");
    
    
    hNumEffPurityBC[2][0][1] =(TH1F*)list->FindObject("hNumEffPurityBC2_V00_Flag1");
    hDenomEffBC[2][0][1]=(TH1F*)list->FindObject("hDenomEffBC2_V00_Flag1");
    hDenomPurityBC[2][0][1]=(TH1F*)list->FindObject("hDenomPurityBC2_V00_Flag1");
    hDenomRejecEffBC[2][0][1]=(TH1F*)list->FindObject("hDenomRejecEffBC2_V00_Flag1");
    hNumRejecEffBC[2][0][1]=(TH1F*)list->FindObject("hNumRejecEffBC2_V00_Flag1");
    
    hPurityBC[2][0][1] = new TH1F("hPurityBC2_V00_Flag1","; #V0flags in PF", 35, 0, 35);
    hPurityBC[2][0][1]->SetTitle("Purity");
    hEfficiencyBC[2][0][1] = new TH1F("hEfficiencyBC2_V00_Flag1","; #V0flags in PF", 35, 0, 35);
    hEfficiencyBC[2][0][1]->SetTitle("Efficiency");
    hRejectEfficiencyBC[2][0][1] = new TH1F("hRejectEfficiencyBC2_V00_Flag1","; #V0flags in PF", 35, 0, 35);
    hRejectEfficiencyBC[2][0][1]->SetTitle("Rejection-Efficiency");
    
    //____________________HM________________________
    TH1F *hNumEffPurityBC_HM[3][3][3];
    TH1F *hDenomEffBC_HM[3][3][3];
    TH1F *hDenomPurityBC_HM[3][3][3];
    TH1F *hDenomRejecEffBC_HM[3][3][3];
    TH1F *hNumRejecEffBC_HM[3][3][3];
    
    TH1F *hPurityBC_HM[3][3][3];
    TH1F *hEfficiencyBC_HM[3][3][3];
    TH1F *hRejectEfficiencyBC_HM[3][3][3];
    
    hNumEffPurityBC_HM[2][0][1] =(TH1F*)list2->FindObject("hNumEffPurityBC_HM2_V00_Flag1");
    hDenomEffBC_HM[2][0][1]=(TH1F*)list2->FindObject("hDenomEffBC_HM2_V00_Flag1");
    hDenomPurityBC_HM[2][0][1]=(TH1F*)list2->FindObject("hDenomPurityBC_HM2_V00_Flag1");
    hDenomRejecEffBC_HM[2][0][1]=(TH1F*)list2->FindObject("hDenomRejecEffBC_HM2_V00_Flag1");
    hNumRejecEffBC_HM[2][0][1]=(TH1F*)list2->FindObject("hNumRejecEffBC_HM2_V00_Flag1");
    
    hPurityBC_HM[2][0][1] = new TH1F("hPurityBC_HM2_V00_Flag1","; #V0flags in PF", 35, 0, 35);
    hPurityBC_HM[2][0][1]->SetTitle("Purity");
    hEfficiencyBC_HM[2][0][1] = new TH1F("hEfficiencyBC_HM2_V00_Flag1","; #V0flags in PF", 35, 0, 35);
    hEfficiencyBC_HM[2][0][1]->SetTitle("Efficiency");
    hRejectEfficiencyBC_HM[2][0][1] = new TH1F("hRejectEfficiencyBC_HM2_V00_Flag1","; #V0flags in PF", 35, 0, 35);
    hRejectEfficiencyBC_HM[2][0][1]->SetTitle("Rejection-Efficiency");
    

    
    
    
    //____________________MB________________________
    double denomeff = hDenomEffBC[2][0][1]->GetBinContent(10);
    double denomrejeceff = hDenomRejecEffBC[2][0][1]->GetBinContent(10);
    
    hPurityBC[2][0][1]->Divide(hNumEffPurityBC[2][0][1], hDenomPurityBC[2][0][1], 1, 1, "B");
    hPurityBC[2][0][1]->GetYaxis()->SetRangeUser(0,1.0);
    
    hEfficiencyBC[2][0][1]->Divide(hNumEffPurityBC[2][0][1], hDenomEffBC[2][0][1], 1, 1, "B");
    hEfficiencyBC[2][0][1]->GetYaxis()->SetRangeUser(0,1.0);
    
    hRejectEfficiencyBC[2][0][1]->Divide(hNumRejecEffBC[2][0][1], hDenomRejecEffBC[2][0][1], 1, 1, "B");
    hRejectEfficiencyBC[2][0][1]->GetYaxis()->SetRangeUser(0,1.0);

    //____________________HM________________________
    double denomeff_HM = hDenomEffBC_HM[2][0][1]->GetBinContent(10);
    double denomrejeceff_HM = hDenomRejecEffBC_HM[2][0][1]->GetBinContent(10);
    
    hPurityBC_HM[2][0][1]->Divide(hNumEffPurityBC_HM[2][0][1], hDenomPurityBC_HM[2][0][1], 1, 1, "B");
    hPurityBC_HM[2][0][1]->GetYaxis()->SetRangeUser(0,1.0);
    
    hEfficiencyBC_HM[2][0][1]->Divide(hNumEffPurityBC_HM[2][0][1], hDenomEffBC_HM[2][0][1], 1, 1, "B");
    hEfficiencyBC_HM[2][0][1]->GetYaxis()->SetRangeUser(0,1.0);
    
    hRejectEfficiencyBC_HM[2][0][1]->Divide(hNumRejecEffBC_HM[2][0][1], hDenomRejecEffBC_HM[2][0][1], 1, 1, "B");
    hRejectEfficiencyBC_HM[2][0][1]->GetYaxis()->SetRangeUser(0,1.0);
    
    //____________________RunbyRun_Analysis_Result________________________
    ResultValue[i][0] = Temp;
    //    ResultValue[i][0] = runNumber_hist->GetBinContent(1);
//    ResultValue[i][0] = i; //temporary, will be modifyed soon.
    ResultValue[i][1] = hPurityBC[2][0][1]->GetBinContent(10);
    ResultValue[i][2] = hEfficiencyBC[2][0][1]->GetBinContent(10);
    ResultValue[i][3] = hRejectEfficiencyBC[2][0][1]->GetBinContent(10);
    ResultValue[i][4] = hPurityBC_HM[2][0][1]->GetBinContent(10);
    ResultValue[i][5] = hEfficiencyBC_HM[2][0][1]->GetBinContent(10);
    ResultValue[i][6] = hRejectEfficiencyBC_HM[2][0][1]->GetBinContent(10);

    filein->Close();
//    cout << "close root file" << endl;
}