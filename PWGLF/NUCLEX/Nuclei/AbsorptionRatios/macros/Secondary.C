void LabelDCAHisto(TH1F *Histo, TString Title, EColor color){
    
    Histo->SetMarkerColor(color);
    Histo->SetMarkerSize(0.8);
    Histo->SetMarkerStyle(20);
    Histo->SetTitle(Title.Data());
    Histo->GetXaxis()->SetRangeUser(-0.7, 0.7);
    Histo->GetXaxis()->SetTitle("DCAxy (cm)");
    Histo->GetYaxis()->SetTitle("Counts");
    Histo->GetXaxis()->SetTitleOffset(1.0);
    
}

void MakeDCAPlot(TString Name,TString PtRangeName, TH1F *ProjDCAData, TH1F *ProjPrimMC, TH1F *ProjSecMC, TH1F *ProjMatMC, TH1F *result,Int_t PtBin){
    
    
    LabelDCAHisto(ProjDCAData, Form("Data %s",Name.Data()) ,kBlack);
    LabelDCAHisto(ProjPrimMC, "Primary",kRed);
    LabelDCAHisto(ProjSecMC, "Secondary",kMagenta);
    LabelDCAHisto(ProjMatMC, "Material",kGreen);
    
    TCanvas *DCAPlot = new TCanvas("DCAPlot","DCAPlot",0,0,800,600);
    TString pdfDCAName = Form("%s/DCAFit/%s_DCAFit_Bin_%d.pdf",Name.Data(),Name.Data(),PtBin);
    TLatex PtBinRange;
    gPad->SetLogy();
    
    ProjDCAData->DrawCopy("Ep");
    ProjSecMC->DrawCopy("same");
    ProjPrimMC->DrawCopy("same");
    ProjMatMC->DrawCopy("same");
    if(result)result->DrawCopy("same");
    gPad->BuildLegend(0.35,0.9,0.15,0.7);
    PtBinRange.SetNDC(kTRUE);
    PtBinRange.SetTextSize(0.035);
    PtBinRange.DrawLatex(.15,.6,PtRangeName);
    
    DCAPlot->Print(pdfDCAName);
    delete DCAPlot;
    return;
}

void MakeDCAPlot(TString Name,TString PtRangeName, TH1F *ProjDCAData, TF1 *Signal, TF1 *Background, Int_t PtBin){
    
    LabelDCAHisto(ProjDCAData, Form("Data %s",Name.Data()) ,kBlack);
    
    TCanvas *DCAPlot = new TCanvas("DCAPlot","DCAPlot",0,0,800,600);
    TString pdfDCAName = Form("%s/DCAFit/%s_DCAFit_Bin_%d.pdf",Name.Data(),Name.Data(),PtBin);
    TLatex PtBinRange;
    gPad->SetLogy();
    gStyle->SetOptStat("nic");
    gStyle->SetOptFit(1100);
    gStyle->SetOptTitle(0);
    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.1);
    
    ProjDCAData->DrawCopy("Ep");
    Signal->DrawCopy("lsame");
    Background->DrawCopy("lsame");
    gPad->BuildLegend(0.35,0.9,0.15,0.7);
    PtBinRange.SetNDC(kTRUE);
    PtBinRange.SetTextSize(0.035);
    PtBinRange.DrawLatex(.15,.6,PtRangeName);
    
    DCAPlot->Print(pdfDCAName);
    delete DCAPlot;
    return;
}

void DCAFit(TH1F* histo,Int_t PtBin,Int_t ParticleType,TF1 *Signal,TF1 *BackGr,Double_t *ValueErrorArray)
{
    TF1 *DCAxyFit = new TF1("DCAxyFit",DCAFitFunction,-2.0,2.0,9);
    DCAxyFit->SetLineColor(kBlue+2);
    
    if(PtBin==5){
        DCAxyFit->SetParameters(-30,-1,20,200,0,0.08,700,0,0.02);
        DCAxyFit->SetParLimits(4,-0.1,0.1);
        DCAxyFit->SetParLimits(7,-0.1,0.1);
    }
    DCAxyFit->SetParameters(-10,-1,100,1500,0,0.1,500,0,0.05);
    DCAxyFit->SetParLimits(4,-0.1,0.1);
    DCAxyFit->SetParLimits(7,-0.1,0.1);
    if(PtBin>10){
        DCAxyFit->SetParameters(0,1,50,2400,0,0.07,200,0,0.04);
        DCAxyFit->SetParLimits(4,-0.1,0.1);
        DCAxyFit->SetParLimits(7,-0.1,0.1);
    }
    if(PtBin>17){
        DCAxyFit->SetParameters(0,0,3,1500,0,0.02,200,0,0.05);
        DCAxyFit->SetParLimits(4,-0.1,0.1);
        DCAxyFit->SetParLimits(7,-0.1,0.1);
    }
    
    auto Result = histo->Fit("DCAxyFit","SB","",-1.0,1.0);
    
    Signal->SetLineColor(kRed+2);
    BackGr->SetLineColor(kMagenta);
    
    Double_t par[9];
    DCAxyFit->GetParameters(par);
    BackGr->SetParameters(par);
    Signal->SetParameters(&par[3]);
    
    //calculation of fraction of primaries and respective error
    Double_t Binlow = -0.1;
    Double_t Binup = 0.1;
    
    Double_t DCAFitInt = DCAxyFit->Integral(Binlow,Binup,0.01);
    Double_t SignalInt = Signal->Integral(Binlow,Binup,0.01);
    Double_t errGlobal = DCAxyFit->IntegralError(Binlow,Binup);
    Double_t errSignal1 = Signal->IntegralError(Binlow,Binup);
    
    Double_t PrimaryFraction = 0;
    if(!(DCAFitInt==0))PrimaryFraction = SignalInt/DCAFitInt;
    Double_t RelError1;
    if(!(SignalInt==0)&&!(DCAFitInt==0))RelError1 = TMath::Sqrt((errSignal1/SignalInt)*(errSignal1/SignalInt)+(errGlobal/DCAFitInt)*(errGlobal/DCAFitInt));
    ValueErrorArray[0] = PrimaryFraction;
    ValueErrorArray[1] = PrimaryFraction*RelError1;
}


void CalculatePrimFrac(TH1F *primary_fraction,TH1F *ProjPrimMC,TH1F *ProjSecMC,TH1F *ProjMatMC,Int_t PtBin){
    
    //Find DCA integration range
    Double_t Binlow = ProjPrimMC->FindBin(-0.1);
    Double_t Binup = ProjPrimMC->FindBin(0.1);
    
    //Calculate Integral and respective Error
    Double_t ErrorPrim;
    Double_t ErrorSec;
    Double_t ErrorMat;
    Double_t Integral_primary = ProjPrimMC->IntegralAndError(Binlow,Binup,ErrorPrim);
    Double_t Integral_secondary = ProjSecMC->IntegralAndError(Binlow,Binup,ErrorSec);
    Double_t Integral_material = ProjMatMC->IntegralAndError(Binlow,Binup,ErrorMat);
    Double_t RelError;
    RelError = TMath::Sqrt( (ErrorPrim/Integral_primary)*(ErrorPrim/Integral_primary) + (ErrorPrim/(Integral_primary+Integral_secondary))*(ErrorPrim/(Integral_primary+Integral_secondary)) );
    
    //Calculate the primary fraction with error
    Double_t Val_PrimFraction = Integral_primary/(Integral_primary+Integral_secondary);
    if(PtBin<18)Val_PrimFraction = Integral_primary/(Integral_primary+Integral_secondary+Integral_material);
    primary_fraction->SetBinContent(PtBin,Val_PrimFraction);
    primary_fraction->SetBinError(PtBin,Val_PrimFraction*RelError);
    return;
}

void MakePrimFracProton(TString Name, TH1F *primary_fraction, TH2F *DCAData,TH2F *PrimaryMC,TH2F *SecondaryMC,TH2F *MaterialMC){
    
    TH1F *ProjDCAData[45];
    TH1F *ProjPrimMC[45];
    TH1F *ProjSecMC[45];
    TH1F *ProjMatMC[45];
    TH1F *result[45];
    TString PtRangeName[45];
    gSystem->Exec(Form("mkdir %s/DCAFit",Name.Data()));
    
    for (Int_t PtBin=4; PtBin<44; PtBin++) {
        std::cout << "*******************************" << std::endl;
        std::cout << "DCA fit of PtBin: " << PtBin << std::endl;
        std::cout << Name.Data() << std::endl;
        std::cout << "*******************************" << std::endl;
        
        Double_t BinLowerEdge = DCAData->GetXaxis()->GetBinLowEdge(PtBin);
        Double_t BinUpperEdge = DCAData->GetXaxis()->GetBinLowEdge(PtBin+1);
        PtRangeName[PtBin] = Form("%1.2f< #it{p} <%1.2f",BinLowerEdge,BinUpperEdge);
        
        //Data
        ProjDCAData[PtBin] = (TH1F*)(DCAData->ProjectionY(Form("DCABin %d %s",PtBin,Name.Data()),PtBin,PtBin,"e"));
        Double_t BinlowData = ProjDCAData[PtBin]->FindBin(-0.5);
        Double_t BinupData = ProjDCAData[PtBin]->FindBin(0.5);
        Double_t Scaling_Intgral = (ProjDCAData[PtBin]->Integral(BinlowData,BinupData));
        
        //Primary
        ProjPrimMC[PtBin] = (TH1F*)(PrimaryMC->ProjectionY(Form("Primbin %s %d",Name.Data(), PtBin),PtBin,PtBin,"e"));
        ProjPrimMC[PtBin]->Scale(Scaling_Intgral/(ProjPrimMC[PtBin]->Integral(BinlowData,BinupData)));
        
        //Secondary
        ProjSecMC[PtBin] = (TH1F*)(SecondaryMC->ProjectionY(Form("Secbin %s %d",Name.Data(), PtBin),PtBin,PtBin,"e"));
        ProjSecMC[PtBin]->Scale(Scaling_Intgral/(ProjSecMC[PtBin]->Integral(BinlowData,BinupData)));
        
        //Material
        ProjMatMC[PtBin] = (TH1F*)(MaterialMC->ProjectionY(Form("Matbin %s %d",Name.Data(), PtBin),PtBin,PtBin,"e"));
        ProjMatMC[PtBin]->Scale(Scaling_Intgral/(ProjMatMC[PtBin]->Integral(BinlowData,BinupData)));
        
        Double_t par0, errpar0, par1, errpar1, par2, errpar2;
        
        TObjArray *mcDCA_arrayLow = new TObjArray(3);  // MC histograms are put in this array
        mcDCA_arrayLow->Add(ProjPrimMC[PtBin]);
        mcDCA_arrayLow->Add(ProjSecMC[PtBin]);
        mcDCA_arrayLow->Add(ProjMatMC[PtBin]);
        mcDCA_arrayLow->ls();
        
        TObjArray *mcDCA_arrayHigh = new TObjArray(2);  // MC histograms are put in this array
        mcDCA_arrayHigh->Add(ProjPrimMC[PtBin]);
        mcDCA_arrayHigh->Add(ProjSecMC[PtBin]);
        mcDCA_arrayHigh->ls();
        
        TFractionFitter* fit;
        if(PtBin<18){
            fit = new TFractionFitter(ProjDCAData[PtBin],mcDCA_arrayLow);  // initialise
        }else{
            fit = new TFractionFitter(ProjDCAData[PtBin],mcDCA_arrayHigh);  // initialise
        }
        
        fit->Constrain(0,0.5,1.0);
        fit->Constrain(1,0.001,0.05);
        if(PtBin<18)fit->Constrain(2,0.0,0.05);
        
        if(PtBin==28||PtBin==22)fit->Constrain(0,0.5,1.0);
        if(PtBin==28||PtBin==22)fit->Constrain(1,0.001,0.07);
        if(strcmp(Name,"Antiproton")==0||PtBin==12)fit->Constrain(0,0.7,1.0);
        if(PtBin==25)fit->Constrain(0,0.6,1.0);
        if(strcmp(Name,"Antiproton")==0||PtBin==25)fit->Constrain(0,0.8,1.0);
        std::cout <<"I am here"<< std::endl;
        
        fit->SetRangeX(ProjDCAData[PtBin]->FindBin(-1.0),ProjDCAData[PtBin]->FindBin(1.0));
        Int_t status = -1;
        status = fit->Fit();// perform the fit
        
        result[PtBin] = 0;
        if (status == 0) {
            // obtain weights:
            fit->GetResult(0, par0, errpar0);
            fit->GetResult(1, par1, errpar1);
            fit->GetResult(2, par2, errpar2);
            
            ProjPrimMC[PtBin]->Scale(par0);
            ProjSecMC[PtBin]->Scale(par1);
            ProjMatMC[PtBin]->Scale(par2);
            result[PtBin] = (TH1F*)fit->GetPlot();
            result[PtBin]->SetTitle(Form("Results %s %d",Name.Data(),PtBin));
        }
        
        CalculatePrimFrac(primary_fraction,ProjPrimMC[PtBin],ProjSecMC[PtBin],ProjMatMC[PtBin],PtBin);
        MakeDCAPlot(Name.Data(),PtRangeName[PtBin],ProjDCAData[PtBin],ProjPrimMC[PtBin],ProjSecMC[PtBin],ProjMatMC[PtBin],result[PtBin],PtBin);
        
        delete mcDCA_arrayLow;
        delete mcDCA_arrayHigh;
        delete fit;
    }
    return;
}

void MakePrimFracDeuteron(TString Name, TH1F *primary_fraction, TH2F *DCAData){
    
    TH1F *ProjDCAData[45];
    TH1F *ProjPrimMC[45];
    TF1 *DCA_Signal[45];
    TF1 *DCA_Background[45];
    TString PurityDCAName[45];
    TString PtRangeName[45];
    gSystem->Exec(Form("mkdir %s/DCAFit",Name.Data()));
    int bins_num[23]  = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 35, 40};
    
    for (Int_t PtBin=6; PtBin<28; PtBin++) {
        
        std::cout << "*******************************" << std::endl;
        std::cout << "DCA fit of PtBin: " << PtBin << std::endl;
        std::cout << Name.Data() << std::endl;
        std::cout << "*******************************" << std::endl;
        
        Double_t BinLowerEdge = DCAData->GetXaxis()->GetBinLowEdge(PtBin);
        Double_t BinUpperEdge = DCAData->GetXaxis()->GetBinLowEdge(PtBin+1);
        PtRangeName[PtBin] = Form("%1.2f< #it{p} <%1.2f",BinLowerEdge,BinUpperEdge);
        
        //Data
        ProjDCAData[PtBin] = (TH1F*)(DCAData->ProjectionY(Form("DCABin %d %s",PtBin,Name.Data()),bins_num[PtBin-6]+1,bins_num[PtBin-5],"e"));
        Double_t BinlowData = ProjDCAData[PtBin]->FindBin(-0.5);
        Double_t BinupData = ProjDCAData[PtBin]->FindBin(0.5);
        Double_t Scaling_Intgral = (ProjDCAData[PtBin]->Integral(BinlowData,BinupData));
        
        Double_t PrimaryFractionErrorArray[2] = {0,0};
        DCA_Signal[PtBin] = new TF1("Signal",DCAGauss,-2.,2.,6);
        DCA_Background[PtBin] = new TF1("BackGr",DCAQuadraticFunc,-2.,2.,3);
        
        if(PtBin>14){
            primary_fraction->SetBinContent(PtBin,1);
        }else{
            DCAFit(ProjDCAData[PtBin],PtBin,1,DCA_Signal[PtBin],DCA_Background[PtBin],PrimaryFractionErrorArray);
            primary_fraction->SetBinContent(PtBin,PrimaryFractionErrorArray[0]);
            primary_fraction->SetBinError(PtBin,PrimaryFractionErrorArray[1]);
            PurityDCAName[PtBin] = Form("Purity: %1.3f",PrimaryFractionErrorArray[0]);
        }
        
        MakeDCAPlot(Name.Data(),PtRangeName[PtBin],ProjDCAData[PtBin], DCA_Signal[PtBin], DCA_Background[PtBin], PtBin);
        
    }
}

void Secondary(){}

