void Summary(){}

void SetLabels(TH1F *Histo, TString xAxis, TString yAxis){
    
    Histo->GetXaxis()->SetTitle(xAxis.Data());
    Histo->GetYaxis()->SetTitle(yAxis.Data());
    Histo->GetXaxis()->SetTitleSize(0.045);
    Histo->GetXaxis()->SetTitleOffset(0.8);
    Histo->GetYaxis()->SetTitleOffset(1.4);
    
}

void MakeQAPlot(TString CName,TString xAxis,TString yAxis,TH1F *hProton,TH1F *hAntiproton,TH1F *hDeuteron,TH1F *hAntideuteron){
    TCanvas *QACanvas = new TCanvas(CName.Data(),CName.Data(),0,0,800,600);
    QACanvas->Divide(2,2);
    
    QACanvas->cd(1);
    SetLabels(hProton,xAxis.Data(),yAxis.Data());
    hProton->DrawCopy("EP");
    
    QACanvas->cd(2);
    SetLabels(hAntiproton,xAxis.Data(),yAxis.Data());
    hAntiproton->DrawCopy("EP");
    
    QACanvas->cd(3);
    SetLabels(hDeuteron,xAxis.Data(),yAxis.Data());
    hDeuteron->DrawCopy("EP");
    
    QACanvas->cd(4);
    SetLabels(hAntideuteron,xAxis.Data(),yAxis.Data());
    hAntideuteron->DrawCopy("EP");
    
}

TH1F *GetPrimarySpectra(TString Name, TH1F *RawSpectrum,TH1F *Fraction){
    
    RawSpectrum->GetXaxis()->SetRangeUser(0,4.0);
    Fraction->GetXaxis()->SetRangeUser(0,4.0);
    
    TH1F *PrimarySpectrum = (TH1F*)(RawSpectrum->Clone(Form("PrimarySpectrum %s",Name.Data())));
    PrimarySpectrum->SetTitle(Form("Primary %s",Name.Data()));
    if(Fraction->Integral()>0) PrimarySpectrum->Multiply(Fraction);
    
    return PrimarySpectrum;
}

TH1F *GetPrimaryRatios(TString Name, TH1F *Particle,TH1F *Antiparticle){
    
    TH1F *Ratio = (TH1F*)(Antiparticle->Clone(Name.Data()));
    Ratio->SetTitle(Form("Primary Ratio %s",Name.Data()));
    Ratio->Divide(Particle);
    Ratio->SetMarkerStyle(kFullSquare);
    return Ratio;
    
}
