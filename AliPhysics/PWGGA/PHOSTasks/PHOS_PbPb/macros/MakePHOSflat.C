void MakePHOSflat(const Char_t *run = "167915")
{

  //Fills PHOS flattening parameters into OADB
  //Each run-dependent object contains list of 3 objects:
  //Flattening for TPC, V0A, V0C
  //Input histograms are read from the file Pi0Flow_*.root
  // which is the output of the analysis task AliAnalysisTaskPi0Flow
  // Produced root file should go to
  //"$ALICE_ROOT/OADB/PHOS/PHOSflat.root"

  /* $Id$ */

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWGGAPHOSTasks");


  TFile * f = TFile::Open(Form("Pi0Flow_%s.root",run)) ;
  TList * l = (TList*)f->Get("PHOSPi0Flow") ;
  TH2F *h1 = (TH2F*)l->FindObject("phiRP") ;
  TH2F *h2 = (TH2F*)l->FindObject("phiRPV0A") ;
  TH2F *h3 = (TH2F*)l->FindObject("phiRPV0C") ;
  
  TF1 * flat = new TF1("flat2","[0]*(1.+2.*[1]*cos(x)+2.*[2]*sin(x)+2.*[3]*cos(2.*x)+2.*[4]*sin(2.*x)+2.*[5]*cos(3.*x)+2.*[6]*sin(3.*x)+2.*[7]*cos(4.*x)+2.*[8]*sin(4.*x))",0.,5.) ;
  
  TH2D * hTPC = new TH2D("hTPC","TPC EP",20,0.,100.,8,0.,8.) ;
  TH2D * hV0A = new TH2D("hV0A","V0A EP",20,0.,100.,8,0.,8.) ;
  TH2D * hV0C = new TH2D("hV0C","V0C EP",20,0.,100.,8,0.,8.) ;
  
  Int_t nCent = h1->GetNbinsY();
  for(Int_t i=1; i<=nCent; i++){
    TH1D * h1D = h1->ProjectionX("a",i,i) ;
    if (h1D->Integral() > 0) {
      h1D->Sumw2() ;
      h1D->Fit(flat,"0") ;
      for(Int_t j=1; j<9; j++)
	hTPC->SetBinContent(i,j,flat->GetParameter(j)) ;
    }
    delete h1D ;
    
    h1D = h2->ProjectionX("a",i,i) ;
    if (h1D->Integral() > 0) {
      h1D->Sumw2() ;
      h1D->Fit(flat,"0") ;
      for(Int_t j=1; j<9; j++)
	hV0A->SetBinContent(i,j,flat->GetParameter(j)) ;   
    }
    delete h1D ;
    
    h1D = h3->ProjectionX("a",i,i) ;
    if (h1D->Integral() > 0) {
      h1D->Sumw2() ;
      h1D->Fit(flat,"0") ;
      for(Int_t j=1; j<9; j++)
	hV0C->SetBinContent(i,j,flat->GetParameter(j)) ;   
    }
    delete h1D ;
  }
   
  AliOADBContainer calibContainer("phosFlat");

  AliPHOSEPFlattener * tpc = new AliPHOSEPFlattener("TPC") ;
  tpc->SetParameterization(hTPC) ;
  AliPHOSEPFlattener * v0a = new AliPHOSEPFlattener("V0A") ;
  v0a->SetParameterization(hV0A) ;
  AliPHOSEPFlattener * v0c = new AliPHOSEPFlattener("V0C") ;
  v0c->SetParameterization(hV0C) ;
  

  // 
  TObjArray * lhc11h = new TObjArray(3); 
  lhc11h->SetName("PHOSflat_LHC11h");
  lhc11h->AddAt(tpc,0) ; //pass 3 reconstruction
  lhc11h->AddAt(v0a,1) ; //pass 3 reconstruction
  lhc11h->AddAt(v0c,2) ; //pass 3 reconstruction
  calibContainer.AppendObject(lhc11h,167693,170593) ;

  calibContainer.WriteToFile("PHOSCflat.root");


}
