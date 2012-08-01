void MakeFits(Int_t mod=1){

  TFile *f = TFile::Open("histos_pass15.root") ;
  char key[55] ;
  sprintf(key,"hMggDispM%d",mod) ;
  TH3F * h3Re = (TH3F*)f->Get(key) ;
  sprintf(key,"hMggDispM%d",mod) ;
  TH3F * h3Re = (TH3F*)f->Get(key) ;
  sprintf(key,"hMiMggAllM%d",mod) ;
  TH3F * h3Mi = (TH3F*)f->Get(key) ;

  TH2D * mass = (TH2D*)h3Re->Project3D("yx") ;
  sprintf(key,"Mass_mod%d",mod) ;
  mass->SetName(key) ;
  mass->Reset();
  TH2D * width = (TH2D*)h3Re->Project3D("yx") ;
  sprintf(key,"Width_mod%d",mod) ;
  width->SetName(key) ;
  width->Reset();

  TH2D * yield = (TH2D*)h3Re->Project3D("yx") ;
  sprintf(key,"Yield_mod%d",mod) ;
  yield->SetName(key) ;
  yield->Reset();

  TSpectrum *s = new TSpectrum(1);
  TCanvas * c = new TCanvas("fit") ;
  TF1 * gs = new TF1("gs","[0]*exp(-(x-[1])*(x-[1])/2./[2]/[2])+[3]+0*[4]*x",0.,1.) ;
  for(Int_t i=1;i<=mass->GetNbinsX();i++){
    for(Int_t j=1;j<=mass->GetNbinsY();j++){
      printf("i=%d, j=%d \n",i,j) ;
      TH1D * tmp = h3Re->ProjectionZ("a",i,i,j,j) ;
      TH1D * tmpMi = h3Mi->ProjectionZ("b",i,i,j,j) ;
      tmp->Sumw2() ;
      //      tmp->Rebin(4) ;
      //      tmpMi->Rebin(4) ;
      //Normalize
      //      Double_t nMi=tmpMi->Integral(70/4,120/4) ;
      //      Double_t nRe=tmp->Integral(70/4,120/4) ;
      Double_t nMi=tmpMi->Integral(70,120) ;
      Double_t nRe=tmp->Integral(70,120) ;
      if(nMi>0.)
	tmpMi->Scale(nRe/nMi) ;
      tmp->Add(tmpMi,-1.) ;

      tmp->GetXaxis()->SetRangeUser(0.05,0.25) ;
      Int_t nfound = s->Search(tmp,1,"new");
      printf("Found %d candidate peaks to fitn",nfound);
      Float_t *xpeaks = s->GetPositionX();       

      gs->SetParameters(8.,xpeaks[0],0.008,0.,0.) ;
      gs->SetParLimits(0,0.,100.) ;
      gs->SetParLimits(1,0.05,0.2) ;
      gs->SetParLimits(2,0.003,0.015) ;
      if(tmp->Integral(40,80)){
	tmp->Fit(gs,"q","",0.05,0.220) ;
	tmp->Fit(gs,"qM","",0.05,0.220) ;
	mass->SetBinContent(i,j,gs->GetParameter(1)) ;
	mass->SetBinError(i,j,gs->GetParError(1)) ;
	width->SetBinContent(i,j,gs->GetParameter(2)) ;
	width->SetBinError(i,j,gs->GetParError(2)) ;
	yield->SetBinContent(i,j,TMath::Sqrt(6.28)*gs->GetParameter(0)*gs->GetParameter(2)/tmp->GetXaxis()->GetBinWidth(1)) ;
	
	if(gs->GetParameter(0)<3){//no visible peak
	  mass->SetBinContent(i,j,0.136) ;
	  mass->SetBinError(i,j,0.) ;
	  width->SetBinContent(i,j,0.) ;
	  width->SetBinError(i,j,0.) ;
	}	
	/*
	//	if(gs->GetParameter(0)<3){
	  tmp->Fit(gs,"","",0.05,0.220) ;
	  c->Update() ;
	  if(getchar()=='q')return ;
	  //	}
	  */

      }
      //      tmpMi->SetLineColor(2) ;
      //      tmpMi->Draw("same") ;
      //            c->Update() ;
      //            if(getchar()=='q')return ;
      delete tmp ;
      delete tmpMi ;
    }
  }

  TFile fout("Calibration_pass16.root","update") ;
  mass->Write(0,TObject::kOverwrite) ;
  width->Write(0,TObject::kOverwrite) ;
  yield->Write(0,TObject::kOverwrite) ;
  
  fout.Close() ;
 

}
