DrawMeanPt01(){
    //Example Macro to Compute Mean Pt evolution for a variety of particle species
    gStyle->SetOptStat(0);
    
    //Hopefully this will be expanded in the near future
    TString lModel = "Pythia6-Perugia0-New";
    
    TFile *fin = new TFile(Form( "AnalysisResults-%s.root", lModel.Data() ),"READ");
    fin -> cd("PWGLF_PPVsMultXCheck_MC");
    TList *flist = fin-> FindObjectAny( "cList" );
    
    //Variable binning in the raw estimator value
    //Will be converted into Mid-pseudorapidity density
    Double_t lBoundaries[] = { 0,3,5,7,9,11,14,18,22,26,30,35,40,45,50,58,66,75,85,95/*,105,115,125*/};

    Int_t lNBins = sizeof (lBoundaries) / sizeof(Double_t) - 1 ;
    
    //Particles for which one wants to draw <pT> evolution
    const int lNPart = 9;
    //Actually, I could have Xi and Omega too but Perugia 0 is HORRIBLE for them so stat is not there...
    //(histos are there though)
    
    //========================================================================
    //Good list
    //TString lPartName [lNPart] = {"Pion", "Kaon", "Proton","K0Short", "Lambda", "Phi"};
    //TString lPartNameLatex [lNPart] = {"#pi","K","p","K^{0}_{S}","#Lambda","#phi"};
    //Int_t lColors [lNPart] = {kRed, kBlue, kGreen, kCyan, kMagenta, kGray+1};
    //========================================================================
    
    //========================================================================
    //Full LF particle listing - Will result in a busy plot + requires lots of stat
    TString lPartName [lNPart] = {"Pion", "Kaon", "Proton","K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar"};
    TString lPartNameLatex [lNPart] = {"#pi","K","p","K^{0}_{S}","#Lambda","#Xi", "#Omega", "#phi", "K*"};
    Int_t lColors [lNPart] = {kRed, kBlue, kGreen, kCyan, kMagenta, kGray+1, kOrange, kTeal, kViolet};
    //========================================================================
    
    //Compute Event Counts
    TH2F* hV0MTrueVsNch = (TH2F*) flist -> FindObject("fHistV0MTrueVsMidRapidityTrue_MCSelection");
    
    //Project
    TH1D* hV0MTrue = hV0MTrueVsNch -> ProjectionX ( "hV0MTrue" );
    
    Double_t lEventCount[100];
    Double_t lEpsilon = 1e-6;
    Double_t lNchVal[100];
    TH1D *hNch[100];
    for( Int_t ib = 0; ib<lNBins; ib++){
        //Epsilon: ensures correct bins are always picked up
        lEventCount[ib] = hV0MTrue -> Integral ( hV0MTrue->FindBin( lBoundaries[ib] ), hV0MTrue->FindBin( lBoundaries[ib+1] - 1 ) );
        hNch[ib] = hV0MTrueVsNch -> ProjectionY ( Form( "hNch_%i", ib ), hV0MTrueVsNch->GetXaxis()->FindBin(lBoundaries[ib] ),
                                                 hV0MTrueVsNch->GetXaxis()->FindBin(lBoundaries[ib+1] - 1 ) );
        lNchVal[ib] = hNch[ib] -> GetMean(); //correct for bin center position, please
    }
    
    //2D Particle Histograms
    TH2F *hPartV0M_Data[lNPart];
    TH2F *hPartV0M_MC[lNPart];
    for( Int_t ip = 0; ip<lNPart; ip++){
        //Monte Carlo Selection of INEL>0
        hPartV0M_MC[ip]   = (TH2F*) flist -> FindObject(Form("fHistPtVsV0MTrue_MCSelection_%s"  ,lPartName[ip].Data()));
    }
    
    //Compute transverse momentum spectra
    //First Index: species
    //Second Index: Multiplicity Interval
    TH1D *hPtV0M_MC[lNPart][100];
    
    const Int_t NRGBs = 5;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.9*0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 0.9*1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red,
                                                green, blue, lNBins);
    
    TH1D *hPtV0M_Ratio[lNPart][100];
    Double_t lMeanPt[lNPart][100];
    Double_t lMeanPtError[lNPart][100];
    for( Int_t ib = 0; ib<lNBins; ib++){
        //Particle Species Loop
        for( Int_t ip = 0; ip<lNPart; ip++){
            hPtV0M_MC[ip][ib] = hPartV0M_MC[ip] -> ProjectionX( Form("hPtV0M_MC_%s_%i",lPartName[ip].Data(),ib ),
                    hPartV0M_MC[ip]->GetYaxis()->FindBin(lBoundaries[ib]   ),
                    hPartV0M_MC[ip]->GetYaxis()->FindBin(lBoundaries[ib+1] - 1)
                                                                   );
            hPtV0M_MC[ip][ib] -> Scale(1./lEventCount[ib]); // not strictly necessary for <pT>
            lMeanPt[ip][ib] = hPtV0M_MC[ip][ib] -> GetMean();
            lMeanPtError[ip][ib] = hPtV0M_MC[ip][ib] -> GetMeanError();
        }
    }
    
    //Draw TGraphErrors
    TCanvas *c1 = new TCanvas("c1", "",800,600);
    c1->SetTicks(1,1);
    c1->SetBottomMargin(0.135);
    c1->SetLeftMargin(0.135);
    c1->SetRightMargin(0.02);
    c1->SetTopMargin(0.02);
 
    TH1F *hdum = new TH1F("hdum", "",1000,0,50);
    hdum->GetYaxis()->SetRangeUser(0,2.1);
    hdum->GetXaxis()->SetRangeUser(0,36);
    hdum->Draw();
    hdum->GetYaxis()->SetTitle("#LTp_{T}#GT");
    hdum->GetXaxis()->SetTitle("dN_{ch}/d#eta (|#eta|<0.5)");
    hdum->GetYaxis()->SetTitleSize(0.055);
    hdum->GetXaxis()->SetTitleSize(0.055);
     
    TGraphErrors *gp[lNPart];
    for (Int_t ip =0; ip<lNPart; ip++){
        gp[ip] = new TGraphErrors ( lNBins, lNchVal, lMeanPt[ip],0, lMeanPtError[ip] );
        //gp[ip]->SetMarkerStyle(20);
        gp[ip]->SetMarkerColor(lColors[ip]);
        gp[ip]->SetLineColor(lColors[ip]);
        gp[ip]->SetLineWidth(2);
        //option
        gp[ip]->SetLineStyle(7);
        gp[ip]->Draw("L");
    }

    TLegend *leg = new TLegend(0.829, 0.202, 0.945, 0.575);
    leg->SetMargin(0.65);
    leg->SetBorderSize(0);
    for(Int_t ip=0; ip<lNPart; ip++) leg->AddEntry(gp[ip],Form("%s", lPartNameLatex[ip].Data()),"l");
    leg->Draw();
    
    TLatex Tl;
    Tl.SetTextSize(0.035);
    Tl.SetNDC();
    //Tl.DrawLatex(0.84,0.45, lModel.Data());
    
    c1->SaveAs( Form("MeanPt_%s.pdf", lModel.Data() ) );

    
}
