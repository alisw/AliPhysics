#include <stdio.h>

void ValueFromRoot(int i, const char *filename);
Double_t GetFractionOfBG(TH2 * hClsVsTrkOrigin, TH2 * hClsVsTrkInput, Int_t typeOfSpdCut, Int_t typeOfResult);
#define buf 128
Double_t ResultValue[500][10][3]; // 

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
        cout << "RunNumber       :" << ResultValue[j][0][0] << endl;
        cout << "==CINT7==" << endl;
        cout << "Purity(None, PF2, PF10)        :" << ResultValue[j][1][0] << ", " << ResultValue[j][2][0] << ", " << ResultValue[j][3][0] << endl;
        cout << "Efficiency(None, PF2, PF10)  :" << ResultValue[j][4][0] << ", " << ResultValue[j][5][0] << ", " << ResultValue[j][6][0] << endl;
        cout << "Reject Eff.(None, PF2, PF10) :" << ResultValue[j][7][0] << ", " << ResultValue[j][8][0] << ", " << ResultValue[j][9][0] << endl;
        cout << "==V0M==" << endl;
        cout << "Purity(None, PF2, PF10)        :" << ResultValue[j][1][1] << ", " << ResultValue[j][2][1] << ", " << ResultValue[j][3][1] << endl;
        cout << "Efficiency(None, PF2, PF10)  :" << ResultValue[j][4][1] << ", " << ResultValue[j][5][1] << ", " << ResultValue[j][6][1] << endl;
        cout << "Reject Eff.(None, PF2, PF10) :" << ResultValue[j][7][1] << ", " << ResultValue[j][8][1] << ", " << ResultValue[j][9][1] << endl;
        cout << "==SH2==" << endl;
        cout << "Purity(None, PF2, PF10)        :" << ResultValue[j][1][2] << ", " << ResultValue[j][2][2] << ", " << ResultValue[j][3][2] << endl;
        cout << "Efficiency(None, PF2, PF10)  :" << ResultValue[j][4][2] << ", " << ResultValue[j][5][2] << ", " << ResultValue[j][6][2] << endl;
        cout << "Reject Eff.(None, PF2, PF10) :" << ResultValue[j][7][2] << ", " << ResultValue[j][8][2] << ", " << ResultValue[j][9][2] << endl;
        cout << endl;
        j++;
    }
    
    const char* fNames[] = {"CINT7","V0M","SH2"};
    const char* pfValues[] = {"","PF2","PF10"};
    const char* resultValues[] = {"Purity","Efficiency","RejectEffciency"};
    
    double noffiles = ii;
    
    TH1F *hRunbyRun[3][3][3];
    for(Int_t iValue = 0 ; iValue<3 ; iValue++){
        for(Int_t iPF = 0 ; iPF<3 ; iPF++){                
            for(Int_t iTrig = 0 ; iTrig<3 ; iTrig++){
                if(iPF ==0){
                    hRunbyRun[iValue][iPF][iTrig] =  new TH1F(Form("hRunbyRun%s_%s",resultValues[iValue],fNames[iTrig]),resultValues[iValue], noffiles, 0, noffiles);
                    hRunbyRun[iValue][iPF][iTrig]->GetYaxis()->SetRangeUser(0,1.0);
                }
                else{
                    hRunbyRun[iValue][iPF][iTrig] =  new TH1F(Form("hRunbyRun%s_%s_%s",resultValues[iValue],fNames[iTrig],pfValues[iPF]),resultValues[iValue], noffiles, 0, noffiles);
                    hRunbyRun[iValue][iPF][iTrig]->GetYaxis()->SetRangeUser(0,1.0);
                }
                hRunbyRun[iValue][iPF][iTrig]->SetLineWidth(2);
                hRunbyRun[iValue][iPF][iTrig]->SetLineColor(iTrig+1); // CINT7 = Black, V0M = Red, SH2 = Green
                hRunbyRun[iValue][iPF][iTrig]->SetLineStyle(iPF); // None = line, PF2 = dot line, PF10 = many dot line
            }
        }
    }
    int k = 0;
    TString RunN;
    while(1){
        RunN = Form("%.0f",ResultValue[k][0][0]);
        if(!ResultValue[k][1][0]) break;
        for(Int_t iTrig = 0 ; iTrig<3 ; iTrig++){
        hRunbyRun[0][0][iTrig]->SetBinContent(k+1,ResultValue[k][1][iTrig]);
        hRunbyRun[0][0][iTrig]->GetXaxis()->SetBinLabel(k+1, RunN);
        hRunbyRun[0][1][iTrig]->SetBinContent(k+1,ResultValue[k][2][iTrig]);
        hRunbyRun[0][1][iTrig]->GetXaxis()->SetBinLabel(k+1, RunN);
        hRunbyRun[0][2][iTrig]->SetBinContent(k+1,ResultValue[k][3][iTrig]);
        hRunbyRun[0][2][iTrig]->GetXaxis()->SetBinLabel(k+1, RunN);


        hRunbyRun[1][0][iTrig]->SetBinContent(k+1,ResultValue[k][4][iTrig]);
        hRunbyRun[1][0][iTrig]->GetXaxis()->SetBinLabel(k+1, RunN);
        hRunbyRun[1][1][iTrig]->SetBinContent(k+1,ResultValue[k][5][iTrig]);
        hRunbyRun[1][1][iTrig]->GetXaxis()->SetBinLabel(k+1, RunN);
        hRunbyRun[1][2][iTrig]->SetBinContent(k+1,ResultValue[k][6][iTrig]);
        hRunbyRun[1][2][iTrig]->GetXaxis()->SetBinLabel(k+1, RunN);


        hRunbyRun[2][0][iTrig]->SetBinContent(k+1,ResultValue[k][7][iTrig]);
        hRunbyRun[2][0][iTrig]->GetXaxis()->SetBinLabel(k+1, RunN);
        hRunbyRun[2][1][iTrig]->SetBinContent(k+1,ResultValue[k][8][iTrig]);
        hRunbyRun[2][1][iTrig]->GetXaxis()->SetBinLabel(k+1, RunN);
        hRunbyRun[2][2][iTrig]->SetBinContent(k+1,ResultValue[k][9][iTrig]);
        hRunbyRun[2][2][iTrig]->GetXaxis()->SetBinLabel(k+1, RunN );
        }

        k++;
    }
    
    
    TLatex *texbc1 = new TLatex(0.1 , 0.15,"#color[1]{Black line} : Results from CINT7 trigger");
    texbc1->SetTextSize(.04);
    texbc1->SetTextColor(1);
    TLatex *texbc2 = new TLatex(0.1 , 0.1,"#color[2]{Red line} : Results from V0M trigger");
    texbc2->SetTextSize(.04);
    texbc2->SetTextColor(1);
    TLatex *texbc3 = new TLatex(0.1 , 0.05,"#color[3]{Green line} : Results from SH2 trigger");
    texbc3->SetTextSize(.04);
    texbc3->SetTextColor(1);
    
    
    RunByRunResult01= new TCanvas("RunByRunResult01","RunByRunResult01",1000,500); // with Triggers, PF=2 condition
    gStyle->SetOptStat(0);
    RunByRunResult01->Divide(3,1);
    
    RunByRunResult01->cd(1);

    hRunbyRun[0][1][0]->Draw("same");
    hRunbyRun[0][1][1]->Draw("same");
    hRunbyRun[0][1][2]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult01->cd(2);
    
    hRunbyRun[1][1][0]->Draw("same");
    hRunbyRun[1][1][1]->Draw("same");
    hRunbyRun[1][1][2]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult01->cd(3);
    hRunbyRun[2][1][0]->Draw("same");
    hRunbyRun[2][1][1]->Draw("same");
    hRunbyRun[2][1][2]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult01->SaveAs("./Fig/CINT7_PF2_OverallRuns.pdf");
    
    RunByRunResult02= new TCanvas("RunByRunResult02","RunByRunResult02",1000,500); // with Triggers, PF=10 condition
    gStyle->SetOptStat(0);
    RunByRunResult02->Divide(3,1);
    
    RunByRunResult02->cd(1);

    hRunbyRun[0][2][0]->Draw("same");
    hRunbyRun[0][2][1]->Draw("same");
    hRunbyRun[0][2][2]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult02->cd(2);
    
    hRunbyRun[1][2][0]->Draw("same");
    hRunbyRun[1][2][1]->Draw("same");
    hRunbyRun[1][2][2]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult02->cd(3);
    hRunbyRun[2][2][0]->Draw("same");
    hRunbyRun[2][2][1]->Draw("same");
    hRunbyRun[2][2][2]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult02->SaveAs("./Fig/CINT7_PF10_OverallRuns.pdf");

    RunByRunResult03= new TCanvas("RunByRunResult03","RunByRunResult03",1000,500); // with Triggers, No PF condition
    gStyle->SetOptStat(0);
    RunByRunResult03->Divide(3,1);
    
    RunByRunResult03->cd(1);

    hRunbyRun[0][0][0]->Draw("same");
    hRunbyRun[0][0][1]->Draw("same");
    hRunbyRun[0][0][2]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult03->cd(2);
    
    hRunbyRun[1][0][0]->Draw("same");
    hRunbyRun[1][0][1]->Draw("same");
    hRunbyRun[1][0][2]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult03->cd(3);
    hRunbyRun[2][0][0]->Draw("same");
    hRunbyRun[2][0][1]->Draw("same");
    hRunbyRun[2][0][2]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult03->SaveAs("./Fig/CINT7_NOPF_OverallRuns.pdf");

    RunByRunResult04= new TCanvas("RunByRunResult04","RunByRunResult04",1000,500); // with NOPF, Trigger condition
    gStyle->SetOptStat(0);
    RunByRunResult04->Divide(3,1);
    
    RunByRunResult04->cd(1);

    hRunbyRun[0][0][0]->Draw("same");
    hRunbyRun[0][1][0]->Draw("same");
    hRunbyRun[0][2][0]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult04->cd(2);
    
    hRunbyRun[1][0][0]->Draw("same");
    hRunbyRun[1][1][0]->Draw("same");
    hRunbyRun[1][2][0]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult04->cd(3);
    hRunbyRun[2][0][0]->Draw("same");
    hRunbyRun[2][1][0]->Draw("same");
    hRunbyRun[2][2][0]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult04->SaveAs("./Fig/NOPF_Trigger_OverallRuns.pdf");

    RunByRunResult05= new TCanvas("RunByRunResult05","RunByRunResult05",1000,500); // with PF=2, Trigger condition
    gStyle->SetOptStat(0);
    RunByRunResult05->Divide(3,1);
    
    RunByRunResult05->cd(1);

    hRunbyRun[0][0][1]->Draw("same");
    hRunbyRun[0][1][1]->Draw("same");
    hRunbyRun[0][2][1]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult05->cd(2);
    
    hRunbyRun[1][0][1]->Draw("same");
    hRunbyRun[1][1][1]->Draw("same");
    hRunbyRun[1][2][1]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult05->cd(3);
    hRunbyRun[2][0][1]->Draw("same");
    hRunbyRun[2][1][1]->Draw("same");
    hRunbyRun[2][2][1]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult05->SaveAs("./Fig/PF2_Trigger_OverallRuns.pdf");


    RunByRunResult06= new TCanvas("RunByRunResult06","RunByRunResult06",1000,500); // with PF=10, Trigger condition
    gStyle->SetOptStat(0);
    RunByRunResult06->Divide(3,1);
    
    RunByRunResult06->cd(1);

    hRunbyRun[0][0][2]->Draw("same");
    hRunbyRun[0][1][2]->Draw("same");
    hRunbyRun[0][2][2]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult06->cd(2);
    
    hRunbyRun[1][0][2]->Draw("same");
    hRunbyRun[1][1][2]->Draw("same");
    hRunbyRun[1][2][2]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult06->cd(3);
    hRunbyRun[2][0][2]->Draw("same");
    hRunbyRun[2][1][2]->Draw("same");
    hRunbyRun[2][2][2]->Draw("same");
    texbc1->Draw("same");
    texbc2->Draw("same");
    texbc3->Draw("same");
    
    RunByRunResult06->SaveAs("./Fig/PF10_Trigger_OverallRuns.pdf");
}

void ValueFromRoot(int inputnum, const char *filename){

    int i = inputnum;
    TString fileroute = filename;
    TFile *filein = TFile::Open(fileroute+"/AnalysisResults.root","READ");
    if(filein==NULL) continue;
    cout << "open root file : " << "AnalysisResults.root"<< endl;
    
    TDirectoryFile *dir = (TDirectoryFile*)filein->Get("BeamGasMon");
    
    TList *list = (TList*)dir->Get("cOutputH");
    TTree *TTree = (TTree*)dir->Get("TreeTrack");

    Int_t Temp = 0;
    TTree->SetBranchAddress("runNumber",&Temp);
    TTree->GetEntry(0);

    const char* fNames[] = {"CINT7","V0M","SH2"};
    const char* pfValues[] = {"","PF2","PF10"};
    
    //____________________Open Histos________________________
    TH2F* hTotalTrkVsClsSPID[3][3];
    for(Int_t iTrig = 0 ; iTrig<3 ; iTrig++){
        for(Int_t iPF = 0 ; iPF<3 ; iPF++){
            if (iPF == 0){
                hTotalTrkVsClsSPID[iTrig][0]=(TH2F*)list->FindObject(Form("hTotalTrkVsClsSPID_%s",fNames[iTrig]));
            }
            else{
                hTotalTrkVsClsSPID[iTrig][iPF]=(TH2F*)list->FindObject(Form("hTotalTrkVsClsSPID_%s_%s",fNames[iTrig],pfValues[iPF]));    
            }
        }
    }
    cout << "Loading histograms complete" << endl;
    //____________________CINT7________________________
    TH2F *hTotalTrkVsClsSPID_CINT7_PF2_NP;
    TH2F *hTotalTrkVsClsSPID_CINT7_PF10_NP;
    hTotalTrkVsClsSPID_CINT7_PF2_NP = new TH2F("hTotalTrkVsClsSPID_CINT7_PF2_NP","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID_CINT7_PF10_NP = new TH2F("hTotalTrkVsClsSPID_CINT7_PF10_NP","; Spd : total",140,0,140,500,0,500);

    TH2F *hTotalTrkVsClsSPID_CINT7_PF10_NP;
    hTotalTrkVsClsSPID_CINT7_PF2_NP->Add(hTotalTrkVsClsSPID[0][0],hTotalTrkVsClsSPID[0][1],1,-1);
    hTotalTrkVsClsSPID_CINT7_PF10_NP->Add(hTotalTrkVsClsSPID[0][0],hTotalTrkVsClsSPID[0][2],1,-1);


    Double_t purityNone = GetFractionOfBG(hTotalTrkVsClsSPID[0][0],hTotalTrkVsClsSPID[0][0],0,0);
    Double_t purityPF2 = GetFractionOfBG(hTotalTrkVsClsSPID[0][0],hTotalTrkVsClsSPID[0][1],0,0);
    Double_t purityPF10 = GetFractionOfBG(hTotalTrkVsClsSPID[0][0],hTotalTrkVsClsSPID[0][2],0,0);

    Double_t efficiencyNone = GetFractionOfBG(hTotalTrkVsClsSPID[0][0],hTotalTrkVsClsSPID[0][0],0,1);
    Double_t efficiencyPF2 = GetFractionOfBG(hTotalTrkVsClsSPID[0][0],hTotalTrkVsClsSPID[0][1],0,1);
    Double_t efficiencyPF10 = GetFractionOfBG(hTotalTrkVsClsSPID[0][0],hTotalTrkVsClsSPID[0][2],0,1);
    Double_t rejectionEfficiencyNone = GetFractionOfBG(hTotalTrkVsClsSPID[0][0],hTotalTrkVsClsSPID[0][0],0,2);
    Double_t rejectionEfficiencyPF2 = GetFractionOfBG(hTotalTrkVsClsSPID[0][0],hTotalTrkVsClsSPID_CINT7_PF2_NP,0,2);
    Double_t rejectionEfficiencyPF10 = GetFractionOfBG(hTotalTrkVsClsSPID[0][0],hTotalTrkVsClsSPID_CINT7_PF10_NP,0,2);

    ResultValue[i][0][0] = Temp;
    //    ResultValue[i][0] = runNumber_hist->GetBinContent(1);
//    ResultValue[i][0] = i; //temporary, will be modifyed soon.
    ResultValue[i][1][0] = purityNone;
    ResultValue[i][2][0] = purityPF2;
    ResultValue[i][3][0] = purityPF10;
    ResultValue[i][4][0] = efficiencyNone;
    ResultValue[i][5][0] = efficiencyPF2;
    ResultValue[i][6][0] = efficiencyPF10;
    ResultValue[i][7][0] = rejectionEfficiencyNone;
    ResultValue[i][8][0] = rejectionEfficiencyPF2;
    ResultValue[i][9][0] = rejectionEfficiencyPF10;    


    //____________________V0M________________________
    TH2F *hTotalTrkVsClsSPID_V0M_PF2_NP;
    TH2F *hTotalTrkVsClsSPID_V0M_PF10_NP;
    hTotalTrkVsClsSPID_V0M_PF2_NP = new TH2F("hTotalTrkVsClsSPID_V0M_PF2_NP","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID_V0M_PF10_NP = new TH2F("hTotalTrkVsClsSPID_V0M_PF10_NP","; Spd : total",140,0,140,500,0,500);

    TH2F *hTotalTrkVsClsSPID_V0M_PF10_NP;
    hTotalTrkVsClsSPID_V0M_PF2_NP->Add(hTotalTrkVsClsSPID[1][0],hTotalTrkVsClsSPID[1][1],1,-1);
    hTotalTrkVsClsSPID_V0M_PF10_NP->Add(hTotalTrkVsClsSPID[1][0],hTotalTrkVsClsSPID[1][2],1,-1);


    Double_t purityNone = GetFractionOfBG(hTotalTrkVsClsSPID[1][0],hTotalTrkVsClsSPID[1][0],0,0);
    Double_t purityPF2 = GetFractionOfBG(hTotalTrkVsClsSPID[1][0],hTotalTrkVsClsSPID[1][1],0,0);
    Double_t purityPF10 = GetFractionOfBG(hTotalTrkVsClsSPID[1][0],hTotalTrkVsClsSPID[1][2],0,0);

    Double_t efficiencyNone = GetFractionOfBG(hTotalTrkVsClsSPID[1][0],hTotalTrkVsClsSPID[1][0],0,1);
    Double_t efficiencyPF2 = GetFractionOfBG(hTotalTrkVsClsSPID[1][0],hTotalTrkVsClsSPID[1][1],0,1);
    Double_t efficiencyPF10 = GetFractionOfBG(hTotalTrkVsClsSPID[1][0],hTotalTrkVsClsSPID[1][2],0,1);
    Double_t rejectionEfficiencyNone = GetFractionOfBG(hTotalTrkVsClsSPID[1][0],hTotalTrkVsClsSPID[1][0],0,2);
    Double_t rejectionEfficiencyPF2 = GetFractionOfBG(hTotalTrkVsClsSPID[1][0],hTotalTrkVsClsSPID_V0M_PF2_NP,0,2);
    Double_t rejectionEfficiencyPF10 = GetFractionOfBG(hTotalTrkVsClsSPID[1][0],hTotalTrkVsClsSPID_V0M_PF10_NP,0,2);

    ResultValue[i][0][1] = Temp;
    ResultValue[i][1][1] = purityNone;
    ResultValue[i][2][1] = purityPF2;
    ResultValue[i][3][1] = purityPF10;
    ResultValue[i][4][1] = efficiencyNone;
    ResultValue[i][5][1] = efficiencyPF2;
    ResultValue[i][6][1] = efficiencyPF10;
    ResultValue[i][7][1] = rejectionEfficiencyNone;
    ResultValue[i][8][1] = rejectionEfficiencyPF2;
    ResultValue[i][9][1] = rejectionEfficiencyPF10;    


    //____________________SH2________________________
    TH2F *hTotalTrkVsClsSPID_SH2_PF2_NP;
    TH2F *hTotalTrkVsClsSPID_SH2_PF10_NP;
    hTotalTrkVsClsSPID_SH2_PF2_NP = new TH2F("hTotalTrkVsClsSPID_SH2_PF2_NP","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID_SH2_PF10_NP = new TH2F("hTotalTrkVsClsSPID_SH2_PF10_NP","; Spd : total",140,0,140,500,0,500);

    TH2F *hTotalTrkVsClsSPID_SH2_PF10_NP;
    hTotalTrkVsClsSPID_SH2_PF2_NP->Add(hTotalTrkVsClsSPID[2][0],hTotalTrkVsClsSPID[2][1],1,-1);
    hTotalTrkVsClsSPID_SH2_PF10_NP->Add(hTotalTrkVsClsSPID[2][0],hTotalTrkVsClsSPID[2][2],1,-1);


    Double_t purityNone = GetFractionOfBG(hTotalTrkVsClsSPID[2][0],hTotalTrkVsClsSPID[2][0],0,0);
    Double_t purityPF2 = GetFractionOfBG(hTotalTrkVsClsSPID[2][0],hTotalTrkVsClsSPID[2][1],0,0);
    Double_t purityPF10 = GetFractionOfBG(hTotalTrkVsClsSPID[2][0],hTotalTrkVsClsSPID[2][2],0,0);

    Double_t efficiencyNone = GetFractionOfBG(hTotalTrkVsClsSPID[2][0],hTotalTrkVsClsSPID[2][0],0,1);
    Double_t efficiencyPF2 = GetFractionOfBG(hTotalTrkVsClsSPID[2][0],hTotalTrkVsClsSPID[2][1],0,1);
    Double_t efficiencyPF10 = GetFractionOfBG(hTotalTrkVsClsSPID[2][0],hTotalTrkVsClsSPID[2][2],0,1);
    Double_t rejectionEfficiencyNone = GetFractionOfBG(hTotalTrkVsClsSPID[2][0],hTotalTrkVsClsSPID[2][0],0,2);
    Double_t rejectionEfficiencyPF2 = GetFractionOfBG(hTotalTrkVsClsSPID[2][0],hTotalTrkVsClsSPID_SH2_PF2_NP,0,2);
    Double_t rejectionEfficiencyPF10 = GetFractionOfBG(hTotalTrkVsClsSPID[2][0],hTotalTrkVsClsSPID_SH2_PF10_NP,0,2);


    ResultValue[i][0][2] = Temp;
    ResultValue[i][1][2] = purityNone;
    ResultValue[i][2][2] = purityPF2;
    ResultValue[i][3][2] = purityPF10;
    ResultValue[i][4][2] = efficiencyNone;
    ResultValue[i][5][2] = efficiencyPF2;
    ResultValue[i][6][2] = efficiencyPF10;
    ResultValue[i][7][2] = rejectionEfficiencyNone;
    ResultValue[i][8][2] = rejectionEfficiencyPF2;
    ResultValue[i][9][2] = rejectionEfficiencyPF10;    

    filein->Close();
//    cout << "close root file" << endl;
}

Double_t GetFractionOfBG(TH2 * hClsVsTrkOrigin, TH2 * hClsVsTrkInput, Int_t typeOfSpdCut, Int_t typeOfResult) {
  // returns the fraction of events not bassing the SPD BG cut
  // typeOfSpdCut == 0: standard linear cut (as in physics selection)
  // typeOfSpdCut == 1: "Christoph's" cut
  // typeOfResult == 0: Purity
  // typeOfResult == 1: Efficiency
  // typeOfResult == 2: Rejection Efficiency

  //First, get the !BG and Total entry of Original (No PF) histogram
  Int_t nbinx = hClsVsTrkOrigin->GetNbinsX();
  Int_t nbiny = hClsVsTrkOrigin->GetNbinsY();
  Double_t nBGEvts = 0;

  for(Int_t ibinx = 0; ibinx < nbinx; ibinx++){
    for(Int_t ibiny = 0; ibiny < nbiny; ibiny++){
      Bool_t spdBg = false;
      Float_t trk = hClsVsTrkOrigin->GetXaxis()->GetBinCenter(ibinx);
      Float_t cls = hClsVsTrkOrigin->GetYaxis()->GetBinCenter(ibiny);
      if (typeOfSpdCut == 0 && (cls > 65.0+4.0*trk)) spdBg = true;
      if (typeOfSpdCut == 1) {
        std::cout << "Not yet implemented" << std::endl;
        exit(1);
      }

      if(spdBg) nBGEvts += hClsVsTrkOrigin->GetBinContent(ibinx, ibiny);
    }
  }

  std::cout << "Origin/// Total number of events: " << hClsVsTrkOrigin->GetEntries()
            << ", BG: " << nBGEvts
            << ", Not BG: " <<  hClsVsTrkOrigin->GetEntries()-nBGEvts << std::endl;

  //Second, Get the same values from Input value
  nbinx = hClsVsTrkInput->GetNbinsX();
  nbiny = hClsVsTrkInput->GetNbinsY();
  Double_t nBGEvtsInput = 0;

  for(Int_t ibinx = 0; ibinx < nbinx; ibinx++){
    for(Int_t ibiny = 0; ibiny < nbiny; ibiny++){
      Bool_t spdBg = false;
      Float_t trk = hClsVsTrkInput->GetXaxis()->GetBinCenter(ibinx);
      Float_t cls = hClsVsTrkInput->GetYaxis()->GetBinCenter(ibiny);
      if (typeOfSpdCut == 0 && (cls > 65.0+4.0*trk)) spdBg = true;
      if (typeOfSpdCut == 1) {
        std::cout << "Not yet implemented" << std::endl;
        exit(1);
      }

      if(spdBg) nBGEvtsInput += hClsVsTrkInput->GetBinContent(ibinx, ibiny);
    }
  }


  std::cout << "Input/// Total number of events: " << hClsVsTrkInput->GetEntries()
            << ", BG: " << nBGEvtsInput
            << ", Not BG: " <<  hClsVsTrkInput->GetEntries()-nBGEvtsInput  << std::endl;

  if(typeOfResult==0){
        if(hClsVsTrkInput->GetEntries() != 0){
            std::cout << "Purity: " << (hClsVsTrkInput->GetEntries()-nBGEvtsInput)/hClsVsTrkInput->GetEntries() << std::endl;
        }
        else{
                std::cout << "Purity: Null"<< std::endl;
        }   
    }     
  if(typeOfResult==1){
    if((hClsVsTrkOrigin->GetEntries()-nBGEvts) != 0){
        std::cout << "Efficiency: " << (hClsVsTrkInput->GetEntries()-nBGEvtsInput)/(hClsVsTrkOrigin->GetEntries()-nBGEvts) << std::endl;    
    }
    else{
        std::cout << "Efficiency: Null"<< std::endl;
    }
    }
  if(typeOfResult==2){
    if(nBGEvts != 0){
        std::cout << "Rejection Efficiency: " << nBGEvtsInput/nBGEvts << std::endl;          
    }
    else{
        std::cout << "Rejection Efficiency: Null"<< std::endl;
    }
    }



  if(typeOfResult==0){
    if(hClsVsTrkInput->GetEntries() != 0){
        return (hClsVsTrkInput->GetEntries()-nBGEvtsInput)/hClsVsTrkInput->GetEntries();
    }
    else{
        return 0.001;
    }
  }
  if(typeOfResult==1){
    if((hClsVsTrkOrigin->GetEntries()-nBGEvts) != 0){
            return (hClsVsTrkInput->GetEntries()-nBGEvtsInput)/(hClsVsTrkOrigin->GetEntries()-nBGEvts);
    }
    else{
        return 0.001;
    }
    }
  if(typeOfResult==2){
    if(nBGEvts != 0){
        return nBGEvtsInput/nBGEvts;
    }
    else{
       return 0.001;
    }
}
}
