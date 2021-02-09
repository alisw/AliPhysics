void CompMassOne(Int_t cen=6){

 
  TString centext[8] = {"0-5%", "5-10%","10-20%","20-40%","40-60%","60-80%","0-20%", "0-10%"};
  Int_t col[8]={kRed,kOrange+1,kYellow+1,kGreen+2,kCyan+1,kBlue,kViolet,kRed+3} ;
  Int_t colBg[8]={kRed-8,kOrange-8,kYellow-6,kGreen-8,kCyan-8,kBlue-9,kViolet-8,kRed-10} ;  

  
  
  TFile * f1 = new TFile("mass_MB.root") ;
  
  TF1 * fit = new TF1("a","[0]",0.,20.) ;

  f1->ls() ; 
  
  const Int_t nPID=8 ;
  char cPID[14][12] ;
  sprintf(cPID[0],"All") ;
  sprintf(cPID[1],"Disp");
  sprintf(cPID[2],"CPV") ;
  sprintf(cPID[3],"Both"); 

  char key[255] ;
  TH1D * m1[8] ;
  for(Int_t iPID=0; iPID<4; iPID++){
    sprintf(key,"mass1_GS_Emin3_%s_cen%d",cPID[iPID],cen) ;
//printf("cen=%d, >%s< \n",cen,key) ;    
    m1[iPID]=(TH1D*)f1->Get(key) ;
    
    m1[iPID]->SetMarkerStyle(20) ;
    m1[iPID]->SetMarkerColor(col[iPID]) ;
    m1[iPID]->SetLineColor(col[iPID]) ;
  }


  TFile * fR = new TFile("../raw_MBmix.root") ;
  TH1D * mR[8] ; 
  for(Int_t iPID=0; iPID<4; iPID++){
    sprintf(key,"mass1_GS_Emin3_%s_cen%d",cPID[iPID],cen) ;
    mR[iPID]=(TH1D*)fR->Get(key) ;
    mR[iPID]->SetMarkerStyle(24) ;
    mR[iPID]->SetMarkerColor(col[iPID]) ;
    mR[iPID]->SetLineColor(col[iPID]) ;
  }

  TH1D * box = new TH1D("box","",100,1.,50.) ;
  box->SetMinimum(0.8) ;
  box->SetMaximum(1.2) ;
//   box->SetMinimum(0.125) ;
//   box->SetMaximum(0.145) ;
  box->SetXTitle("p_{t} (GeV/c)") ;
  box->GetYaxis()->SetTitleOffset(1) ;
  box->GetYaxis()->SetTitleSize(0.055) ;
//  box->SetYTitle("m (GeV/c^{2})") ;
  box->SetYTitle("m_{MC}/m_{Data}") ;
//  box->GetYaxis()->SetTitleOffset(1.2) ;
//  box->SetYTitle("m (GeV/c^{2})") ;
  box->SetStats(0) ;
  
  sprintf(key,"Mass_Cen_%s",cPID[iPID]) ; 
  TCanvas * cEff = new TCanvas(key,key,10,10,600,600) ;

  TAxis * xa = m1[0]->GetXaxis() ;
  for(Int_t iPID=0; iPID<1; iPID++){
    cEff->cd(iPID+1) ;
    gPad->SetLogx(1) ;
    gPad->SetGridx() ;
    gPad->SetGridy() ;
    box->SetTitle(cPID[iPID]) ;
    box->DrawClone() ;
    TH1D * tmp = (TH1D*)mR[iPID]->Clone(Form("Mass_%s_cen%d",cPID[iPID],cen)) ;
    tmp->SetMarkerStyle(20) ;
    
    m1[iPID]->Divide(mR[iPID]) ;

    m1[iPID]->Fit(fit,"","",1.,10.) ;   
    box->DrawClone() ;
    m1[iPID]->Draw("same") ;
//    mR[iPID]->Draw("same") ;

TLatex * text = new TLatex(6.,1.07,Form("%5.3f#pm%5.3f",fit->GetParameter(0),fit->GetParError(0))) ;
text->SetTextSize(0.075) ;
text->Draw() ;
  }

  
 
  cEff->cd(2) ;
  
  TLegend * l = new TLegend(0.2,0.7,0.5,0.9) ;
  l->AddEntry(mR[1],"Data","p") ;
  l->AddEntry(m1[1], "MC","p") ;

//   l->Draw() ;

}
