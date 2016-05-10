// this is macro to get "efficiency", "purity" and "rejection efficiency".

void MacroHistosHLTTest(){
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

    TFile *fileout = TFile::Open(ss+"/BGMonitorHisto.root","RECREATE");

    TDirectoryFile *dir = (TDirectoryFile*)filein->Get("BeamGasMon");
    
    TList *list = (TList*)dir->Get("cOutputH_CINT7");
    TList *list2 = (TList*)dir->Get("cOutputH_V0MandSH2"); 
    TCanvas* canvas = new TCanvas("canvas");
    TH2F *hTotalTrkVsClsSPID_CINT7 = (TH2F*)list->FindObject("hTotalTrkVsClsSPID_CINT7");
    hTotalTrkVsClsSPID_CINT7->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_CINT7->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    hTotalTrkVsClsSPID_CINT7->Draw("colz");
    hTotalTrkVsClsSPID_CINT7->Write("hTotalTrkVsClsSPID_CINT7");
    canvas->SaveAs("./Fig/hTotalTrkVsClsSPID_CINT7.png");

    TH2F *hTotalTrkVsClsSPID_CINT7_PF2 = (TH2F*)list->FindObject("hTotalTrkVsClsSPID_CINT7_PF2");
    hTotalTrkVsClsSPID_CINT7_PF2->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_CINT7_PF2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    hTotalTrkVsClsSPID_CINT7_PF2->Draw("colz");
    hTotalTrkVsClsSPID_CINT7_PF2->Write("hTotalTrkVsClsSPID_CINT7");
    canvas->SaveAs("./Fig/hTotalTrkVsClsSPID_CINT7_PF2.png");

    TH2F *hTotalTrkVsClsSPID_CINT7_PF10 = (TH2F*)list->FindObject("hTotalTrkVsClsSPID_CINT7_PF10");
    hTotalTrkVsClsSPID_CINT7_PF10->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_CINT7_PF10->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    hTotalTrkVsClsSPID_CINT7_PF10->Draw("colz");
    hTotalTrkVsClsSPID_CINT7_PF10->Write("hTotalTrkVsClsSPID_CINT7");
    canvas->SaveAs("./Fig/hTotalTrkVsClsSPID_CINT7_PF10.png");

    TH2F *hTotalTrkVsClsSPID_V0M = (TH2F*)list2->FindObject("hTotalTrkVsClsSPID_V0M");
    hTotalTrkVsClsSPID_V0M->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_V0M->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    hTotalTrkVsClsSPID_V0M->Draw("colz");
    hTotalTrkVsClsSPID_V0M->Write("hTotalTrkVsClsSPID_V0M");
    canvas->SaveAs("./Fig/hTotalTrkVsClsSPID_V0M.png");

    TH2F *hTotalTrkVsClsSPID_V0M_PF2 = (TH2F*)list2->FindObject("hTotalTrkVsClsSPID_V0M_PF2");
    hTotalTrkVsClsSPID_V0M_PF2->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_V0M_PF2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    hTotalTrkVsClsSPID_V0M_PF2->Draw("colz");
    hTotalTrkVsClsSPID_V0M_PF2->Write("hTotalTrkVsClsSPID_V0M");
    canvas->SaveAs("./Fig/hTotalTrkVsClsSPID_V0M_PF2.png");

    TH2F *hTotalTrkVsClsSPID_V0M_PF10 = (TH2F*)list2->FindObject("hTotalTrkVsClsSPID_V0M_PF10");
    hTotalTrkVsClsSPID_V0M_PF10->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_V0M_PF10->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    hTotalTrkVsClsSPID_V0M_PF10->Draw("colz");
    hTotalTrkVsClsSPID_V0M_PF10->Write("hTotalTrkVsClsSPID_V0M");
    canvas->SaveAs("./Fig/hTotalTrkVsClsSPID_V0M_PF10.png");

    TH2F *hTotalTrkVsClsSPID_SH2 = (TH2F*)list2->FindObject("hTotalTrkVsClsSPID_SH2");
    hTotalTrkVsClsSPID_SH2->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_SH2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    hTotalTrkVsClsSPID_SH2->Draw("colz");
    hTotalTrkVsClsSPID_SH2->Write("hTotalTrkVsClsSPID_SH2");
    canvas->SaveAs("./Fig/hTotalTrkVsClsSPID_SH2.png");

    TH2F *hTotalTrkVsClsSPID_SH2_PF2 = (TH2F*)list2->FindObject("hTotalTrkVsClsSPID_SH2_PF2");
    hTotalTrkVsClsSPID_SH2_PF2->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_SH2_PF2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    hTotalTrkVsClsSPID_SH2_PF2->Draw("colz");
    hTotalTrkVsClsSPID_SH2_PF2->Write("hTotalTrkVsClsSPID_SH2");
    canvas->SaveAs("./Fig/hTotalTrkVsClsSPID_SH2_PF2.png");

    TH2F *hTotalTrkVsClsSPID_SH2_PF10 = (TH2F*)list2->FindObject("hTotalTrkVsClsSPID_SH2_PF10");
    hTotalTrkVsClsSPID_SH2_PF10->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_SH2_PF10->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    hTotalTrkVsClsSPID_SH2_PF10->Draw("colz");
    hTotalTrkVsClsSPID_SH2_PF10->Write("hTotalTrkVsClsSPID_SH2");
    canvas->SaveAs("./Fig/hTotalTrkVsClsSPID_SH2_PF10.png");

    fileout->Close();
 
}