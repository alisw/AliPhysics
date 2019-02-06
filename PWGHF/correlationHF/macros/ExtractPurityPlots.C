void DrawFractionSecondariesAccepted(TList *l1);
void DrawFractionPrimariesRejected(TList *l1, TList *l2);
void DrawDeltaPhiDistrAndRatios_Secondaries(TList *l1);
void DrawDeltaPhiRatios_CharmAccepted(TList* l1, TList* l2);
void DrawDeltaPhiRatios_BeautyAccepted(TList* l1, TList* l2);
void ReflectPointsIn32Bin(TH1F* h);

void ExtractPurityPlots(TString filename="AnalysisResults.root", TString dirnameWithDCAcut="D0hCorrelOutput_pPb_Purity_0_20_DCA010_MultWeig", TString dirnameNoDCAcut="D0hCorrelOutput_pPb_Purity_0_20_NoDCA_MultWeig", TString listnameWithDCAcut="debugPlotsOutput_pPb_Purity_0_20_DCA010_MultWeig_020", TString listnameNoDCAcut="debugPlotsOutput_pPb_Purity_0_20_NoDCA_MultWeig_020", TString suffix="0_20_DCA010") {
 
  TFile f(filename.Data());
  TDirectoryFile *d1 = (TDirectoryFile*)f.Get(dirnameWithDCAcut);
  TDirectoryFile *d2 = (TDirectoryFile*)f.Get(dirnameNoDCAcut);
  TList *l1 = (TList*)d1->Get(listnameWithDCAcut);
  TList *l2 = (TList*)d2->Get(listnameNoDCAcut);

  gSystem->Exec(Form("rm -r Output_Plots_%s",suffix.Data()));
  gSystem->Exec(Form("rm -r Output_Root_%s",suffix.Data()));
  gSystem->Exec(Form("mkdir Output_Plots_%s",suffix.Data()));
  gSystem->Exec(Form("mkdir Output_Root_%s",suffix.Data()));

  DrawFractionSecondariesAccepted(l1,suffix);
  DrawFractionPrimariesRejected(l1,l2,suffix);
  DrawDeltaPhiDistrAndRatios_Secondaries(l1,suffix);
  DrawDeltaPhiDistrAndRatios_Primaries(l1,suffix);
  DrawDeltaPhiRatios_CharmAccepted(l1,l2,suffix);
  DrawDeltaPhiRatios_BeautyAccepted(l1,l2,suffix);

  return;
} 


void DrawFractionSecondariesAccepted(TList* l1, TString suffix) {  //secondary tracks which pass DCA cut over all tracks which pass DCA cut

  TString namebinD[5] = {"2to3","3to5","5to8","8to16","16to24"};
  TString namebinAss[5] = {"03to99","03to1","1to99","1to3","3to99"};

  for(int i=0; i<5; i++) {

    TH1F *hSecAcc = new TH1F("hSecAcc",Form("Fraction of secondary track accepted, %s",namebinD[i].Data()),5,-0.5,4.5);
    for(int j=0;j<5;j++) hSecAcc->GetXaxis()->SetBinLabel(j+1,namebinAss[j]);
    hSecAcc->Sumw2(); //in realtà no, andrebbe usato l'errore binomiale, ma vabbè tanto è useless...
	
    for(int j=0;j<5;j++) {						    
      TH1F *hInputSec = (TH1F*)l1->FindObject(Form("hPurityCount_SecAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));
      TH1F *hInputPrim = (TH1F*)l1->FindObject(Form("hPurityCount_PrimAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));
      Double_t num = hInputSec->GetBinContent(1);
      Double_t den = hInputPrim->GetBinContent(1)+hInputSec->GetBinContent(1);
      if(den) {
        hSecAcc->SetBinContent(j+1,num/den);
        hSecAcc->SetBinError(j+1,TMath::Sqrt(num/den*(1-num/den)*1/den));
      }
      else hSecAcc->SetBinContent(j+1,1);
    }

    TCanvas *cOut = new TCanvas("c1","cOut",100,100,1200,900);  
    hSecAcc->SetStats(kFALSE);
    hSecAcc->SetLineColor(kBlack);
    hSecAcc->SetMarkerStyle(20);
    hSecAcc->SetMarkerSize(1);
    hSecAcc->Draw();

    cOut->SaveAs(Form("Output_Plots_%s/FractOfSecOverTotal_%s.png",suffix.Data(),namebinD[i].Data()));
    cOut->SaveAs(Form("Output_Root_%s/FractOfSecOverTotal_%s.root",suffix.Data(),namebinD[i].Data()));
  }
}


void DrawFractionPrimariesRejected(TList* l1, TList* l2, TString suffix) {  //primary tracks which pass the DCA cut over primary tracks w/o any DCA cut

  TString namebinD[5] = {"2to3","3to5","5to8","8to16","16to24"};
  TString namebinAss[5] = {"03to99","03to1","1to99","1to3","3to99"};

  for(int i=0; i<5; i++) {

    TH1F *hPrimAcc = new TH1F("hPrimRej",Form("Fraction of primary tracks accepted by DCA cut, %s",namebinD[i].Data()),5,-0.5,4.5);
    TH1F *hcAcc = new TH1F("hcRej",Form("Fraction of charm tracks accepted by DCA cut, %s",namebinD[i].Data()),5,-0.5,4.5);
    TH1F *hbAcc = new TH1F("hbRej",Form("Fraction of beauty tracks accepted by DCA cut, %s",namebinD[i].Data()),5,-0.5,4.5);        
    for(int j=0;j<5;j++) {
    	hSecAcc->GetXaxis()->SetBinLabel(j+1,namebinAss[j]);
    	hcAcc->GetXaxis()->SetBinLabel(j+1,namebinAss[j]);
    	hbAcc->GetXaxis()->SetBinLabel(j+1,namebinAss[j]);
    }
    hSecAcc->Sumw2(); //in realtà no, andrebbe usato l'errore binomiale, ma vabbè tanto è useless...
    hcAcc->Sumw2(); //in realtà no, andrebbe usato l'errore binomiale, ma vabbè tanto è useless...
    hbAcc->Sumw2(); //in realtà no, andrebbe usato l'errore binomiale, ma vabbè tanto è useless...
	
    for(int j=0;j<5;j++) {
      TH1F *hInputPrStd = (TH1F*)l1->FindObject(Form("hPurityCount_PrimAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));
      hInputPrStd->SetName("h1");
      TH1F *hInputPrNo  = (TH1F*)l2->FindObject(Form("hPurityCount_PrimAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));
      hInputPrNo->SetName("h2");
      TH1F *hInputcStd = (TH1F*)l1->FindObject(Form("hPurityCount_CharmAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));
      hInputcStd->SetName("h1");
      TH1F *hInputcNo  = (TH1F*)l2->FindObject(Form("hPurityCount_CharmAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));
      hInputcStd->SetName("h2");
      TH1F *hInputbStd = (TH1F*)l1->FindObject(Form("hPurityCount_BeautyAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));
      hInputbStd->SetName("h1");
      TH1F *hInputbNo  = (TH1F*)l2->FindObject(Form("hPurityCount_BeautyAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));
      hInputbStd->SetName("h2");      

      Double_t numP = hInputPrStd->GetBinContent(1);
      Double_t denP = hInputPrNo->GetBinContent(1);
      Double_t numC = hInputcStd->GetBinContent(1);
      Double_t denC = hInputcNo->GetBinContent(1);
      Double_t numB = hInputbStd->GetBinContent(1);
      Double_t denB = hInputbNo->GetBinContent(1);

      if(denP) {
        hPrimAcc->SetBinContent(j+1,numP/denP); 
        hPrimAcc->SetBinError(j+1,TMath::Sqrt(numP/denP*(1-numP/denP)*1/denP)); 
      }
      else hPrimAcc->SetBinContent(j+1,1);
      if(denC) {
        hcAcc->SetBinContent(j+1,numC/denC); 
        hcAcc->SetBinError(j+1,TMath::Sqrt(numC/denC*(1-numC/denC)*1/denC)); 
      }
      else hcAcc->SetBinContent(j+1,1);
      if(denB) {
        hbAcc->SetBinContent(j+1,numB/denB);
        hbAcc->SetBinError(j+1,TMath::Sqrt(numB/denB*(1-numB/denB)*1/denB)); 
      }
      else hbAcc->SetBinContent(j+1,1);
    }

    TCanvas *cOut1 = new TCanvas("c1","cOutPr",100,100,1200,900);  
    hPrimAcc->SetStats(kFALSE);
    hPrimAcc->SetLineColor(kBlack);
    hPrimAcc->SetMarkerStyle(20);
    hPrimAcc->SetMarkerSize(1);
    hPrimAcc->Draw();
    TCanvas *cOut2 = new TCanvas("c2","cOutC",100,100,1200,900);  
    hcAcc->SetStats(kFALSE);
    hcAcc->SetLineColor(kBlack);
    hcAcc->SetMarkerStyle(20);
    hcAcc->SetMarkerSize(1);
    hcAcc->Draw();
    TCanvas *cOut3 = new TCanvas("c3","cOutB",100,100,1200,900);  
    hbAcc->SetStats(kFALSE);
    hbAcc->SetLineColor(kBlack);
    hbAcc->SetMarkerStyle(20);
    hbAcc->SetMarkerSize(1);
    hbAcc->Draw();    

    cOut1->SaveAs(Form("Output_Plots_%s/FractOfPrimAccepted_%s.png",suffix.Data(),namebinD[i].Data()));
    cOut2->SaveAs(Form("Output_Plots_%s/FractOfCharmAccepted_%s.png",suffix.Data(),namebinD[i].Data()));
    cOut3->SaveAs(Form("Output_Plots_%s/FractOfBeautyAccepted_%s.png",suffix.Data(),namebinD[i].Data()));
    cOut1->SaveAs(Form("Output_Root_%s/FractOfPrimAccepted_%s.root",suffix.Data(),namebinD[i].Data()));
    cOut2->SaveAs(Form("Output_Root_%s/FractOfCharmAccepted_%s.root",suffix.Data(),namebinD[i].Data()));
    cOut3->SaveAs(Form("Output_Root_%s/FractOfBeautyAccepted_%s.root",suffix.Data(),namebinD[i].Data()));
  }
}


void DrawDeltaPhiDistrAndRatios_Primaries(TList* l1, TString suffix) {  //secondary tracks which pass DCA cut over all tracks which pass DCA cut

  TString namebinD[5] = {"2to3","3to5","5to8","8to16","16to24"};
  TString namebinAss[5] = {"03to99","03to1","1to99","1to3","3to99"};

  for(int i=0; i<5; i++) {
    for(int j=0;j<5;j++) {

      TH1F *hInputSec = (TH1F*)l1->FindObject(Form("hPuritydPhi_SecAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));
      TH1F *hInputPrim = (TH1F*)l1->FindObject(Form("hPuritydPhi_PrimAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));

      TH1F *hSum = (TH1F*)hInputPrim->Clone("hSum");
      hSum->Add(hInputSec);
      TH1F *hRatio = (TH1F*)hInputPrim->Clone("hRatio");
      hRatio->Divide(hRatio,hSum,1,1,"B");
      TH1F *hRatioRefl = ReflectPointsIn32Bin(hRatio);

      TCanvas *cSec = new TCanvas("cSec","cSec",100,100,1200,900);  
      hInputSec->Draw();
      TCanvas *cPrim = new TCanvas("cPrim","cPrim",100,100,1200,900);  
      hInputPrim->Draw();
      TCanvas *cRatio = new TCanvas("cRatio","cRatio",100,100,1200,900);  
      hRatioRefl->Draw();      
 
      TF1 *fitf = new TF1("fitf","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)+[8]*pow(x,8)+[9]*pow(x,9)+[1]*(x-6.28)+[2]*pow(x-6.28,2)+[3]*pow(x-6.28,3)+[4]*pow(x-6.28,4)+[5]*pow(x-6.28,5)+[6]*pow(x-6.28,6)+[7]*pow(x-6.28,7)+[8]*pow(x-6.28,8)+[9]*pow(x-6.28,9)+[1]*(x+6.28)+[2]*pow(x+6.28,2)+[3]*pow(x+6.28,3)+[4]*pow(x+6.28,4)+[5]*pow(x+6.28,5)+[6]*pow(x+6.28,6)+[7]*pow(x+6.28,7)+[8]*pow(x+6.28,8)+[9]*pow(x+6.28,9)",-TMath::Pi()/2.,3.*TMath::Pi()/2.);
      //fitf->SetLineWidth(1);
      hRatioRefl->Fit(fitf);
      TPaveText *pt = new TPaveText(2.2,1.3,4.6,1.78);
      pt->SetFillColor(kWhite);

      TH1F *hMovMean=new TH1F("hRatioReflMovMean","hRatioReflMovMean",32,-TMath::Pi()/2.,3.*TMath::Pi()/2.);
   
      for(Int_t k=1;k<= (hRatioRefl->GetNbinsX());k++){
        if(k==1) {
            Double_t y = (hRatioRefl->GetBinContent(k)+hRatioRefl->GetBinContent(32)+hRatioRefl->GetBinContent(k+1))/3.;
            Double_t ey = TMath::Sqrt(hRatioRefl->GetBinError(k)*hRatioRefl->GetBinError(k)+hRatioRefl->GetBinError(32)*hRatioRefl->GetBinError(32)+hRatioRefl->GetBinError(k+1)*hRatioRefl->GetBinError(k+1))/3.;
            hMovMean->SetBinContent(k,y);
            hMovMean->SetBinError(k,ey);
        }
        if(k>=2 && k<=31) {
            Double_t y = (hRatioRefl->GetBinContent(k)+hRatioRefl->GetBinContent(k-1)+hRatioRefl->GetBinContent(k+1))/3.;
            Double_t ey = TMath::Sqrt(hRatioRefl->GetBinError(k)*hRatioRefl->GetBinError(k)+hRatioRefl->GetBinError(k-1)*hRatioRefl->GetBinError(k-1)+hRatioRefl->GetBinError(k+1)*hRatioRefl->GetBinError(k+1))/3.;
            hMovMean->SetBinContent(k,y);
            hMovMean->SetBinError(k,ey);
        }     
        if(k==32) {
            Double_t y = (hRatioRefl->GetBinContent(k)+hRatioRefl->GetBinContent(1)+hRatioRefl->GetBinContent(k-1))/3.;
            Double_t ey = TMath::Sqrt(hRatioRefl->GetBinError(k)*hRatioRefl->GetBinError(k)+hRatioRefl->GetBinError(k-1)*hRatioRefl->GetBinError(k-1)+hRatioRefl->GetBinError(1)*hRatioRefl->GetBinError(1))/3.;
            hMovMean->SetBinContent(k,y);
            hMovMean->SetBinError(k,ey);
        }
      }

      hMovMean->SetLineColor(kRed);
      hMovMean->SetMarkerColor(kRed);
      hMovMean->Draw("same");

      cSec->SaveAs(Form("Output_Plots_%s/DeltaPhi_%s_%s_TrackSec.png",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));
      cPrim->SaveAs(Form("Output_Plots_%s/DeltaPhi_%s_%s_TrackPrim.png",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));
      cRatio->SaveAs(Form("Output_Plots_%s/DeltaPhi_%s_%s_RatioPrimOverAll.png",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));
      cSec->SaveAs(Form("Output_Root_%s/DeltaPhi_%s_%s_TrackSec.root",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));
      cPrim->SaveAs(Form("Output_Root_%s/DeltaPhi_%s_%s_TrackPrim.root",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));
      cRatio->SaveAs(Form("Output_Root_%s/DeltaPhi_%s_%s_RatioPrimOverAll.root",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));

    }
  }
}

void DrawDeltaPhiDistrAndRatios_Secondaries(TList* l1, TString suffix) {  //secondary tracks which pass DCA cut over all tracks which pass DCA cut

  TString namebinD[5] = {"2to3","3to5","5to8","8to16","16to24"};
  TString namebinAss[5] = {"03to99","03to1","1to99","1to3","3to99"};

  for(int i=0; i<5; i++) {
    for(int j=0;j<5;j++) {

      TH1F *hInputSec = (TH1F*)l1->FindObject(Form("hPuritydPhi_SecAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));
      TH1F *hInputPrim = (TH1F*)l1->FindObject(Form("hPuritydPhi_PrimAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));

      TH1F *hSum = (TH1F*)hInputSec->Clone("hSum");
      hSum->Add(hInputPrim);
      TH1F *hRatio = (TH1F*)hInputSec->Clone("hRatio");
      hRatio->Divide(hRatio,hSum,1,1,"B");
      TH1F *hRatioRefl = ReflectPointsIn32Bin(hRatio);

      TCanvas *cSec = new TCanvas("cSec","cSec",100,100,1200,900);  
      hInputSec->Draw();
      TCanvas *cPrim = new TCanvas("cPrim","cPrim",100,100,1200,900);  
      hInputPrim->Draw();
      TCanvas *cRatio = new TCanvas("cRatio","cRatio",100,100,1200,900);  
      hRatioRefl->Draw();    

      TF1 *fitf = new TF1("fitf","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)+[8]*pow(x,8)+[9]*pow(x,9)+[1]*(x-6.28)+[2]*pow(x-6.28,2)+[3]*pow(x-6.28,3)+[4]*pow(x-6.28,4)+[5]*pow(x-6.28,5)+[6]*pow(x-6.28,6)+[7]*pow(x-6.28,7)+[8]*pow(x-6.28,8)+[9]*pow(x-6.28,9)+[1]*(x+6.28)+[2]*pow(x+6.28,2)+[3]*pow(x+6.28,3)+[4]*pow(x+6.28,4)+[5]*pow(x+6.28,5)+[6]*pow(x+6.28,6)+[7]*pow(x+6.28,7)+[8]*pow(x+6.28,8)+[9]*pow(x+6.28,9)",-TMath::Pi()/2.,3.*TMath::Pi()/2.);
      //fitf->SetLineWidth(1);
      hRatioRefl->Fit(fitf);
      TPaveText *pt = new TPaveText(2.2,1.3,4.6,1.78);
      pt->SetFillColor(kWhite);

      TH1F *hMovMean=new TH1F("hRatioReflMovMean","hRatioReflMovMean",32,-TMath::Pi()/2.,3.*TMath::Pi()/2.);
   
      for(Int_t k=1;k<= (hRatioRefl->GetNbinsX());k++){
        if(k==1) {
            Double_t y = (hRatioRefl->GetBinContent(k)+hRatioRefl->GetBinContent(32)+hRatioRefl->GetBinContent(k+1))/3.;
            Double_t ey = TMath::Sqrt(hRatioRefl->GetBinError(k)*hRatioRefl->GetBinError(k)+hRatioRefl->GetBinError(32)*hRatioRefl->GetBinError(32)+hRatioRefl->GetBinError(k+1)*hRatioRefl->GetBinError(k+1))/3.;
            hMovMean->SetBinContent(k,y);
            hMovMean->SetBinError(k,ey);
        }
        if(k>=2 && k<=31) {
            Double_t y = (hRatioRefl->GetBinContent(k)+hRatioRefl->GetBinContent(k-1)+hRatioRefl->GetBinContent(k+1))/3.;
            Double_t ey = TMath::Sqrt(hRatioRefl->GetBinError(k)*hRatioRefl->GetBinError(k)+hRatioRefl->GetBinError(k-1)*hRatioRefl->GetBinError(k-1)+hRatioRefl->GetBinError(k+1)*hRatioRefl->GetBinError(k+1))/3.;
            hMovMean->SetBinContent(k,y);
            hMovMean->SetBinError(k,ey);
        }     
        if(k==32) {
            Double_t y = (hRatioRefl->GetBinContent(k)+hRatioRefl->GetBinContent(1)+hRatioRefl->GetBinContent(k-1))/3.;
            Double_t ey = TMath::Sqrt(hRatioRefl->GetBinError(k)*hRatioRefl->GetBinError(k)+hRatioRefl->GetBinError(k-1)*hRatioRefl->GetBinError(k-1)+hRatioRefl->GetBinError(1)*hRatioRefl->GetBinError(1))/3.;
            hMovMean->SetBinContent(k,y);
            hMovMean->SetBinError(k,ey);
        }
      }

      hMovMean->SetLineColor(kRed);
      hMovMean->SetMarkerColor(kRed);
      hMovMean->Draw("same");

      cSec->SaveAs(Form("Output_Plots_%s/DeltaPhi_%s_%s_TrackSec.png",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));
      cPrim->SaveAs(Form("Output_Plots_%s/DeltaPhi_%s_%s_TrackPrim.png",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));
      cRatio->SaveAs(Form("Output_Plots_%s/DeltaPhi_%s_%s_RatioSecOverAll.png",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));
      cSec->SaveAs(Form("Output_Root_%s/DeltaPhi_%s_%s_TrackSec.root",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));
      cPrim->SaveAs(Form("Output_Root_%s/DeltaPhi_%s_%s_TrackPrim.root",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));
      cRatio->SaveAs(Form("Output_Root_%s/DeltaPhi_%s_%s_RatioSecOverAll.root",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));

    }
  }
}

void DrawDeltaPhiRatios_CharmAccepted(TList* l1, TList* l2, TString suffix) {  //secondary tracks which pass DCA cut over all tracks which pass DCA cut

  TString namebinD[5] = {"2to3","3to5","5to8","8to16","16to24"};
  TString namebinAss[5] = {"03to99","03to1","1to99","1to3","3to99"};

  for(int i=0; i<5; i++) {
    for(int j=0;j<5;j++) {

      TH1F *hInputCAcc = (TH1F*)l1->FindObject(Form("hPuritydPhi_CharmAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));
      TH1F *hInputCAll = (TH1F*)l2->FindObject(Form("hPuritydPhi_CharmAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));

      TH1F *hRatio = (TH1F*)hInputCAcc->Clone("hRatio");
      hRatio->Divide(hInputCAcc,hInputCAll,1,1,"B");
      TH1F *hRatioRefl = ReflectPointsIn32Bin(hRatio);

      TCanvas *cRatio = new TCanvas("cRatio","cRatio",100,100,1200,900);  
      hRatioRefl->Draw();

      TF1 *fitf = new TF1("fitf","[0]",-TMath::Pi()/2.,3*TMath::Pi()/2.);
      //fitf->SetLineWidth(1);
      hRatioRefl->Fit(fitf);
      TPaveText *pt = new TPaveText(2.2,1.3,4.6,1.78);
      pt->SetFillColor(kWhite);
      printf("Fit outcome %1.3f #pm %1.3f\n",fitf->GetParameter(0),fitf->GetParError(0));

      cRatio->SaveAs(Form("Output_Plots_%s/DeltaPhi_CharmAcceptedOverAll_%s_%s.png",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));
      cRatio->SaveAs(Form("Output_Root_%s/DeltaPhi_CharmAcceptedOverAll_%s_%s.root",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));

    }
  }
}

void DrawDeltaPhiRatios_BeautyAccepted(TList* l1, TList* l2, TString suffix) {  //secondary tracks which pass DCA cut over all tracks which pass DCA cut

  TString namebinD[5] = {"2to3","3to5","5to8","8to16","16to24"};
  TString namebinAss[5] = {"03to99","03to1","1to99","1to3","3to99"};

  for(int i=0; i<5; i++) {
    for(int j=0;j<5;j++) {

      TH1F *hInputBAcc = (TH1F*)l1->FindObject(Form("hPuritydPhi_BeautyAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));
      TH1F *hInputBAll = (TH1F*)l2->FindObject(Form("hPuritydPhi_BeautyAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data()));

      TH1F *hRatio = (TH1F*)hInputBAcc->Clone("hRatio");
      hRatio->Divide(hInputBAcc,hInputBAll,1,1,"B");
      TH1F *hRatioRefl = ReflectPointsIn32Bin(hRatio);

      TCanvas *cRatio = new TCanvas("cRatio","cRatio",100,100,1200,900);  
      hRatioRefl->Draw();

      TF1 *fitf = new TF1("fitf","[0]",-TMath::Pi()/2.,3*TMath::Pi()/2.);
      //fitf->SetLineWidth(1);
      hRatioRefl->Fit(fitf);
      TPaveText *pt = new TPaveText(2.2,1.3,4.6,1.78);
      pt->SetFillColor(kWhite);
      printf("Fit outcome %1.3f #pm %1.3f\n",fitf->GetParameter(0),fitf->GetParError(0));

      cRatio->SaveAs(Form("Output_Plots_%s/DeltaPhi_BeautyAcceptedOverAll_%s_%s.png",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));
      cRatio->SaveAs(Form("Output_Root_%s/DeltaPhi_BeautyAcceptedOverAll_%s_%s.root",suffix.Data(),namebinD[i].Data(),namebinAss[j].Data()));

    }
  }
}


//________________________________________________________________________________________________
TH1F* ReflectPointsIn32Bin(TH1F *h){
  
  TH1F *h2=new TH1F(Form("%sReflected",h->GetName()),Form("%sReflected",h->GetName()),h->GetNbinsX(),-TMath::Pi()/2.,3.*TMath::Pi()/2.);

  for(Int_t j=1;j<=h2->GetNbinsX();j++){
    Int_t k = -1;
    if(j<=16) k=16-j+1;
    else k=32-j+17;
    Double_t y = h->GetBinContent(j)+h->GetBinContent(k);
    Double_t ey = TMath::Sqrt(h->GetBinError(j)*h->GetBinError(j)+h->GetBinError(k)*h->GetBinError(k));
    h2->SetBinContent(j,y);
    h2->SetBinError(j,ey);
    printf("Combining bins %d-%d (%1.3f #pm %1.4f with %1.4f #pm %1.4f) to get %1.4f #pm %1.4f\n",j,k,h->GetBinContent(j),h->GetBinError(j),h->GetBinContent(k),h->GetBinError(k),y/2.,ey/2.);
  }
  h2->Scale(0.5);

  return h2;
}
