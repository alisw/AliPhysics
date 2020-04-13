void MakeMassPlot(TString Name, TString PtRangeName,TH1F *Mass2sq,TF1 *Mass_Signal,TF1 *Mass_Background,Int_t PtBin, Double_t TheoMass){
    TLatex ptBinRange;
    TCanvas *MassPlot = new TCanvas("MassPlot","MassPlot",0,0,800,600);
    TString pdfMassName = Form("%s/MassFit/%s_MassFit_Bin_%d.pdf",Name.Data(),Name.Data(),PtBin);
    gPad->SetLogy();
    
    Mass2sq->Draw("ep");
    Mass2sq->GetXaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
    Mass2sq->GetXaxis()->SetTitleSize(0.045);
    Mass2sq->GetXaxis()->SetTitleOffset(0.75);
    Mass2sq->GetYaxis()->SetTitle("Counts");
    Mass2sq->GetXaxis()->SetRangeUser(TheoMass-2,TheoMass+1.5);
    ptBinRange.SetNDC(kTRUE);
    ptBinRange.SetTextSize(0.035);
    ptBinRange.DrawLatex(.15,.7,PtRangeName);
    Mass_Signal->DrawCopy("same");
    Mass_Background->DrawCopy("same");
    
    MassPlot->Print(pdfMassName);
    delete MassPlot;
}
void MassFitProton(TString fname, TH1F* histo,Int_t PtBin,TF1* Signal,TF1* BackGr, TH1F *hQA){
    
    TF1 *MassFit = new TF1("MassFit",MassFitFunctionProton,-1.5,5.5,8);
    //Parameters Deuteron
    
    if(PtBin<20){
        MassFit->SetParameters(30000,0.87,0.05,1,1000,-1000,10000000,-20);
    }else if(PtBin<33){
        MassFit->SetParameters(200000,0.87,0.06,0.1,1000,-1000,10000000,-10);
    }else{
        MassFit->SetParameters(200000,0.87,0.06,0.1,1000,-1000,10000000,-10);
        if(PtBin==39)MassFit->SetParameters(200000,0.87,0.06,0.1,1000,-1000,10000000,-5);
    }
    
    TFitter::SetPrecision(10);
    MassFit->SetLineColor(kBlack);
    TFitResultPtr Result;
    Result = histo->Fit(MassFit,"SB","",0.2,2.2);
    histo->GetYaxis()->SetTitle("Counts");
    Double_t par[8];
    
    MassFit->GetParameters(par);
    Signal->SetParameters(par);
    BackGr->SetParameters(&par[4]);
    Signal->SetLineColor(kBlue);
    
    Double_t MassFitInt = 0;
    Double_t SignalInt = 0;
    Double_t sigma;
    Double_t gaussmean;
    sigma = par[2];
    gaussmean = par[1];
    
    Double_t Intlow = gaussmean - (3*sigma);
    Double_t Intup = gaussmean + (3*sigma);
    
    //calculation of Signal and respective error
    Double_t PurityPtbin;
    SignalInt = Signal->Integral(Intlow,Intup,0.001)/(histo->GetXaxis()->GetBinWidth(1));
    MassFitInt = MassFit->Integral(Intlow,Intup,0.001)/(histo->GetXaxis()->GetBinWidth(1));
    if(!(MassFitInt==0))PurityPtbin = SignalInt/MassFitInt;
    
    Double_t errGlobal = MassFit->IntegralError(Intlow,Intup);
    Double_t errSignal1 = Signal->IntegralError(Intlow,Intup);
    Double_t RelError1;
    if(!(SignalInt==0)&&!(MassFitInt==0))RelError1 = TMath::Sqrt((errSignal1/SignalInt)*(errSignal1/SignalInt)+(errGlobal/MassFitInt)*(errGlobal/MassFitInt));
    hQA->SetBinContent(PtBin,PurityPtbin);
    hQA->SetBinError(PtBin,PurityPtbin*RelError1);
}



void MassFitDeuteron(TString fname, TH1F* histo,Int_t PtBin,TF1* Signal,TF1* BackGr, TH1F *hQA){
    
    TF1 *MassFit = new TF1("MassFit",MassFitFunction,-1.5,5.5,8);
    
    //Parameters Deuteron
    if(PtBin==9){
        MassFit->SetParameters(200,3.5,0.15,0.3,1,0,500,-10);
    }else if(PtBin<20){
        MassFit->SetParameters(60,3.6,0.15,0.3,1,-1,5000,-7);
    }else{
        MassFit->SetParameters(30,3.5,0.15,1.3,100,-15,40000,-3);
        if(strcmp(fname,"Deuteron")==0&&PtBin==25) MassFit->SetParameters(30,3.5,0.8,10,150,-30,400000,-4);
        if(strcmp(fname,"Antideuteron")==0&&PtBin==25) MassFit->SetParameters(30,3.5,0.4,10,150,-30,400000,-4);
        
        if(strcmp(fname,"Deuteron")==0&&PtBin==26) MassFit->SetParameters(30,3.5,0.31,10,100,-20,400000,-5);
        if(strcmp(fname,"Antideuteron")==0&&PtBin==26) MassFit->SetParameters(30,3.55,0.06,10,150,-40,400000,-3);
    }
    
    TFitter::SetPrecision(10);
    MassFit->SetLineColor(kBlack);
    TFitResultPtr Result;
    Result = histo->Fit(MassFit,"SB","",1.5,5.1);
    histo->GetYaxis()->SetTitle("Counts");
    Double_t par[8];
    
    MassFit->GetParameters(par);
    Signal->SetParameters(par);
    BackGr->SetParameters(&par[4]);
    Signal->SetLineColor(kBlue);
    
    Double_t MassFitInt = 0;
    Double_t SignalInt = 0;
    Double_t sigma;
    Double_t gaussmean;
    sigma = par[2];
    gaussmean = par[1];
    
    Double_t Intlow = gaussmean - (3*sigma);
    Double_t Intup = gaussmean + (3*sigma);
    
    //calculation of Signal and respective error
    Double_t PurityPtbin;
    SignalInt = Signal->Integral(Intlow,Intup,0.001)/(histo->GetXaxis()->GetBinWidth(1));
    MassFitInt = MassFit->Integral(Intlow,Intup,0.001)/(histo->GetXaxis()->GetBinWidth(1));
    if(!(MassFitInt==0))PurityPtbin = SignalInt/MassFitInt;
    
    Double_t errGlobal = MassFit->IntegralError(Intlow,Intup);
    Double_t errSignal1 = Signal->IntegralError(Intlow,Intup);
    Double_t RelError1;
    if(!(SignalInt==0)&&!(MassFitInt==0))RelError1 = TMath::Sqrt((errSignal1/SignalInt)*(errSignal1/SignalInt)+(errGlobal/MassFitInt)*(errGlobal/MassFitInt));
    hQA->SetBinContent(PtBin,PurityPtbin);
    hQA->SetBinError(PtBin,PurityPtbin*RelError1);
    return;
}
void MakeRawSpectra(TH1F *SpectrumHist,TF1* Signal,Double_t PtBin, Double_t binwidth_mass){
    Double_t SignalInt = Signal->Integral(-5,5,0.001)/(binwidth_mass);
    Double_t SignalError = TMath::Sqrt(SignalInt);
    SpectrumHist->SetBinContent(PtBin,SignalInt);
    SpectrumHist->SetBinError(PtBin,SignalError);
    return;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RawParticleSpectraProton(TString Name,TH1F *MomentumData, TH1F *hMomentum,TH2F *Mass2sqData, TH1F *QAMassFit){
    TH1F *Mass2sq[45];
    TString PtRangeName[45];
    TF1 *Mass_Signal[45];
    TF1 *Mass_Background[45];
    gSystem->Exec(Form("mkdir %s/MassFit",Name.Data()));
    
    for (Int_t PtBin=3; PtBin<8; PtBin++) {
        hMomentum->SetBinContent(PtBin,MomentumData->GetBinContent(PtBin));
    }
    
    for(Int_t PtBin=5; PtBin<45; PtBin++){
        
        std::cout << "*******************************" << std::endl;
        std::cout << "Mass fit of PtBin: " << PtBin << std::endl;
        std::cout << Name.Data() << std::endl;
        std::cout << "*******************************" << std::endl;
        
        TString MassName = Form("%s Bin %d",Name.Data(),PtBin);
        Double_t BinLowerEdge = Mass2sqData->GetXaxis()->GetBinLowEdge(PtBin);
        Double_t BinUpperEdge = Mass2sqData->GetXaxis()->GetBinLowEdge(PtBin+1);
        PtRangeName[PtBin] = Form("%1.2f< #it{p} <%1.2f",BinLowerEdge,BinUpperEdge);
        Mass2sq[PtBin] = (TH1F*)(Mass2sqData->ProjectionY(MassName.Data(),PtBin,PtBin));
        
        Mass_Signal[PtBin] = new TF1("Signal",gaussianSignal,-1.5,5.5,4);
        Mass_Background[PtBin] = new TF1("BackGr",background,-1.5,5.5,4);
        
        MassFitProton(Name.Data(),Mass2sq[PtBin],PtBin,Mass_Signal[PtBin],Mass_Background[PtBin],QAMassFit);
        if(PtBin>7)MakeRawSpectra(hMomentum,Mass_Signal[PtBin],PtBin,Mass2sq[PtBin]->GetXaxis()->GetBinWidth(10));
        
        MakeMassPlot(Name.Data(),PtRangeName[PtBin],Mass2sq[PtBin],Mass_Signal[PtBin],Mass_Background[PtBin],PtBin,0.88);
    }
    return;
}

void RawParticleSpectraDeuteron(TString Name,TH1F *MomentumData, TH1F *hMomentum,TH2F *Mass2sqData, TH1F *QAMassFit){
    TH1F *Mass2sq[45];
    TString PtRangeName[45];
    TF1 *Mass_Signal[45];
    TF1 *Mass_Background[45];
    gSystem->Exec(Form("mkdir %s",Name.Data()));
    gSystem->Exec(Form("mkdir %s/MassFit",Name.Data()));
    
    int bins_num[20]  = {8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 35, 40};
    for (Int_t PtBin=5; PtBin<15; PtBin++) {
        hMomentum->SetBinContent(PtBin,MomentumData->GetBinContent(PtBin));
    }
    
    for(Int_t PtBin=9; PtBin<28; PtBin++){
        
        std::cout << "*******************************" << std::endl;
        std::cout << "Mass fit of PtBin: " << PtBin << std::endl;
        std::cout << Name.Data() << std::endl;
        std::cout << "*******************************" << std::endl;
        
        TString MassName = Form("%s Bin %d",Name.Data(),PtBin);
        Double_t BinLowerEdge = Mass2sqData->GetXaxis()->GetBinLowEdge(bins_num[PtBin-9]+1);
        Double_t BinUpperEdge = Mass2sqData->GetXaxis()->GetBinLowEdge(bins_num[PtBin-8]+1);
        PtRangeName[PtBin] = Form("%1.2f< #it{p} <%1.2f",BinLowerEdge,BinUpperEdge);
        Mass2sq[PtBin] = (TH1F*)(Mass2sqData->ProjectionY(MassName.Data(),bins_num[PtBin-9]+1,bins_num[PtBin-8]));
        
        Mass_Signal[PtBin] = new TF1("Signal",gaussianSignal,-1.5,5.5,4);
        Mass_Background[PtBin] = new TF1("BackGr",background,-1.5,5.5,4);
        
        MassFitDeuteron(Name.Data(),Mass2sq[PtBin],PtBin,Mass_Signal[PtBin],Mass_Background[PtBin],QAMassFit);
        if(PtBin>14)MakeRawSpectra(hMomentum,Mass_Signal[PtBin],PtBin,Mass2sq[PtBin]->GetXaxis()->GetBinWidth(10));
        
        MakeMassPlot(Name.Data(),PtRangeName[PtBin],Mass2sq[PtBin],Mass_Signal[PtBin],Mass_Background[PtBin],PtBin,3.5);
    }
    hMomentum->Scale(0.1, "width");
    return;
}

void RawParticleSpectra(TH1F *hMomentum,TH2F *Mass2sqData){}
