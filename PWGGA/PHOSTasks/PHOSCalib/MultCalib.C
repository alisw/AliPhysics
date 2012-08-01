void MultCalib(Int_t nIter=16){

  TH2D * c[6]={0} ;
  for(Int_t i=0; i<=nIter; i++){
    TFile * f = new TFile(Form("Calibration_pass%d.root",i)) ;
    for(Int_t m=1; m<4; m++){
       TH2D * tmp = (TH2D*)f->Get(Form("Mass_mod%d",m)) ;
       if(c[m]==0){
         c[m]=(TH2D*)tmp->Clone(Form("MassC_mod%d",m)) ;
	 for(Int_t ix=1; ix<=64;ix++)
	   for(Int_t iz=1; iz<=56; iz++){
	     Double_t clb= tmp->GetBinContent(ix,iz) ;
	     if(clb>0){
               if(clb<0.2)
	         c[m]->SetBinContent(ix,iz,0.136/clb) ;
               else
	         c[m]->SetBinContent(ix,iz,clb) ;
             }
	     else
	       c[m]->SetBinContent(ix,iz,1.) ;
	   }
       }	 
       else{
	 for(Int_t ix=1; ix<=64;ix++)
	   for(Int_t iz=1; iz<=56; iz++){
	     Double_t clb= tmp->GetBinContent(ix,iz) ;
	     if(clb>0)
               if(clb<0.2)
                 if(i<11)
                   clb=0.136/clb ;
                 else
                   clb=(0.136/clb)*(0.136/clb) ;
	       c[m]->SetBinContent(ix,iz,c[m]->GetBinContent(ix,iz)*clb) ;
	   }
       }
    }
  }

  TFile * fout = new TFile("calib.root","recreate") ;
  c[1]->Write("Mass_mod1") ;
  c[2]->Write("Mass_mod2") ;
  c[3]->Write("Mass_mod3") ;
}
