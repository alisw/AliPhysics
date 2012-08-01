void CompWidth(){

  Int_t nIter=18;
  TFile *f[20] ;
  TH3F * hReM1[20] ;
  TH3F * hMiM1[20] ;
  TH3F * hReM2[20] ;
  TH3F * hMiM2[20] ;
  TH3F * hReM3[20] ;
  TH3F * hMiM3[20] ;
  char key[155] ;

  for(Int_t i=0;i<nIter;i++){
    sprintf(key,"histos_pass%d.root",i) ;
    f[i] = TFile::Open(key) ;
    hReM1[i] = (TH3F*)f[i]->Get("hMggDispM1") ;
    sprintf(key,"hMggDispM1_iter%d",i) ;
    hReM1[i]->SetName(key) ;
    hMiM1[i] = (TH3F*)f[i]->Get("hMiMggDispM1") ;
    sprintf(key,"hMiMggDispM1_iter%d",i) ;
    hMiM1[i]->SetName(key) ;

    hReM2[i] = (TH3F*)f[i]->Get("hMggDispM2") ;
    sprintf(key,"hMggDispM2_iter%d",i) ;
    hReM2[i]->SetName(key) ;
    hMiM2[i] = (TH3F*)f[i]->Get("hMiMggDispM2") ;
    sprintf(key,"hMiMggDispM2_iter%d",i) ;
    hMiM2[i]->SetName(key) ;

    hReM3[i] = (TH3F*)f[i]->Get("hMggDispM3") ;
    sprintf(key,"hMggDispM3_iter%d",i) ;
    hReM3[i]->SetName(key) ;
    hMiM3[i] = (TH3F*)f[i]->Get("hMiMggDispM3") ;
    sprintf(key,"hMiMggDispM3_iter%d",i) ;
    hMiM3[i]->SetName(key) ;
  }

  TH1D * hW1 = new TH1D("WidthM1","Module 1",nIter,-0.5,nIter-0.5) ;
  TH1D * hW2 = new TH1D("WidthM2","Module 1",nIter,-0.5,nIter-0.5) ;
  TH1D * hW3 = new TH1D("WidthM3","Module 1",nIter,-0.5,nIter-0.5) ;

//  TF1 * gs = new TF1("gs","[0]*exp(-(x-[1])*(x-[1])/2./[2]/[2])+[3]+[4]*x",0.,1.) ;
  TF1 * gs = new TF1("gs","[0]*exp(-(x-[1])*(x-[1])/2./[2]/[2])+[3]",0.,1.) ;
  TCanvas * c = new TCanvas("c","c") ;

  TH1D * hD[100] ; 

  for(Int_t i=0;i<nIter;i++){
    TH1D * tmp = hReM1[i]->ProjectionZ("a") ;
    TH1D * tmpMi = hMiM1[i]->ProjectionZ("b") ;
    tmp->Sumw2() ;
    //Normalize
    Double_t nMi=tmpMi->Integral(70,120) ;
    Double_t nRe=tmp->Integral(70,120) ;
    if(nMi>0.)
      tmpMi->Scale(nRe/nMi) ;
    tmp->Add(tmpMi,-1.) ;
    gs->SetParameters(3500.,0.135,0.008,0.,0.) ;
    //    gs->SetParLimits(0,0.,1.e+6) ;
    gs->SetParLimits(1,0.09,0.2) ;
    gs->SetParLimits(2,0.003,0.015) ;
    
    c->cd() ;
//    tmp->Fit(gs,"q","",0.05,0.220) ;
    tmp->GetXaxis()->SetRangeUser(0.,0.3) ;
    tmp->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})") ;
    tmp->SetTitle("Module 2") ;
    tmp->Fit(gs,"","",0.125,0.145) ;
    hD[i]=(TH1D*)tmp->Clone(Form("Copy%d",i)) ;
//    c->Update() ; if(getchar()=='q')return  ;
//    if(i==12) return ;
    //    mass->SetBinContent(i,j,gs->GetParameter(1)) ;
    //	mass->SetBinError(i,j,gs->GetParError(1)) ;
    hW1->SetBinContent(i+1,gs->GetParameter(2)) ;
    hW1->SetBinError(i+1,gs->GetParError(2)) ;
    delete tmp ;
    delete tmpMi ;

    tmp = hReM2[i]->ProjectionZ("a") ;
    tmpMi = hMiM2[i]->ProjectionZ("b") ;
    tmp->Sumw2() ;
    //Normalize
    nMi=tmpMi->Integral(70,120) ;
    nRe=tmp->Integral(70,120) ;
    if(nMi>0.)
      tmpMi->Scale(nRe/nMi) ;
    tmp->Add(tmpMi,-1.) ;
    gs->SetParameters(10000.,0.150,0.008,0.,0.) ;
    //    gs->SetParLimits(0,0.,100.) ;
    gs->SetParLimits(1,0.09,0.2) ;
    gs->SetParLimits(2,0.003,0.015) ;
    
    tmp->Fit(gs,"q","",0.125,0.145) ;
    //    mass->SetBinContent(i,j,gs->GetParameter(1)) ;
    //	mass->SetBinError(i,j,gs->GetParError(1)) ;
    hW2->SetBinContent(i+1,gs->GetParameter(2)) ;
    hW2->SetBinError(i+1,gs->GetParError(2)) ;
    delete tmp ;
    delete tmpMi ;

    tmp = hReM3[i]->ProjectionZ("a") ;
    tmpMi = hMiM3[i]->ProjectionZ("b") ;
    tmp->Sumw2() ;
    //Normalize
    Double_t nMi=tmpMi->Integral(70,120) ;
    Double_t nRe=tmp->Integral(70,120) ;
    if(nMi>0.)
      tmpMi->Scale(nRe/nMi) ;
    tmp->Add(tmpMi,-1.) ;
    gs->SetParameters(10000.,0.150,0.008,0.,0.) ;
    //    gs->SetParLimits(0,0.,100.) ;
    gs->SetParLimits(1,0.09,0.2) ;
    gs->SetParLimits(2,0.003,0.015) ;
    
//    tmp->Fit(gs,"q","",0.05,0.220) ;
    tmp->Fit(gs,"q","",0.125,0.145) ;
    //    mass->SetBinContent(i,j,gs->GetParameter(1)) ;
    //	mass->SetBinError(i,j,gs->GetParError(1)) ;
    hW3->SetBinContent(i+1,gs->GetParameter(2)) ;
    hW3->SetBinError(i+1,gs->GetParError(2)) ;
    delete tmp ;
    delete tmpMi ;

  }

  hW1->SetMarkerStyle(20) ;
  hW1->SetMarkerColor(2) ;
  hW2->SetMarkerStyle(21) ;
  hW2->SetMarkerColor(4) ;
  hW3->SetMarkerStyle(24) ;
  hW3->SetMarkerColor(8) ;

  hW1->SetXTitle("Iteration") ;
  hW1->SetYTitle("#sigma (GeV/c^{2})") ;
  hW1->Draw() ;
  hW2->Draw("same") ;
  hW3->Draw("same") ;

  TLegend * l = new TLegend(0.7,0.7,0.85,0.85) ;
  l->AddEntry(hW3,"Module 2","p") ;
  l->AddEntry(hW2,"Module 3","p") ;
  l->AddEntry(hW1,"Module 4","p") ;
  l->Draw() ;

  TCanvas * cMinv = new TCanvas("cMinv") ;
  hD[0]->SetMarkerStyle(20) ;
  hD[0]->SetMarkerSize(0.8) ;
  hD[0]->Draw() ;
  for(Int_t i=1; i<nIter; i++){
   hD[i]->SetMarkerStyle(20+i) ;
   hD[i]->SetMarkerSize(0.8) ;
   hD[i]->Draw("same") ;
  }

}
