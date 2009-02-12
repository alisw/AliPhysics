void TestMultiVector(){

  // Example of usage of AliMultiDimVector and AliSignificanceCalculator classes

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libAOD");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWG3");
  gSystem->Load("libPWG3vertexingHF");

  Int_t nptbins=4;
  const Int_t npars=3;
  Int_t nofcells[npars]={20,20,20};
  Float_t looses[npars]={0.,700.,0.8};
  Float_t tights[npars]={1500.,0.,1.};
  TString AxisTitle[npars];
  AxisTitle[0]="DistPS (micron)";
  AxisTitle[1]="TrackDispersion";
  AxisTitle[2]="CosPoint";

  AliMultiDimVector* mvsig=new AliMultiDimVector("Signal","Signal",npars,nptbins,nofcells,looses,tights,AxisTitle);
  AliMultiDimVector* mvsig2=new AliMultiDimVector("Signal","Signal",npars,nptbins,nofcells,looses,tights,AxisTitle);
  AliMultiDimVector* mvbkg=new AliMultiDimVector("Background","Background",npars,nptbins,nofcells,looses,tights,AxisTitle);


  TF1* dsig=new TF1("dsig","exp(-x/310.)",0.,5000.);
  TF1* dbak=new TF1("dbak","exp(-x/50)",0.,5000.);
  gStyle->SetOptStat(2);
  for(Int_t isignev=0;isignev<5000;isignev++){
    Int_t ptbin=gRandom->Integer(nptbins);
    Float_t dists=dsig->GetRandom();
    Float_t distb=dbak->GetRandom();
    Float_t cpas=1-TMath::Abs(gRandom->Gaus(0.,0.02));
    Float_t cpab=0.8+gRandom->Rndm()*0.2;
    Float_t sigverts=gRandom->Gaus(200.,25.);
    Float_t sigvertb=gRandom->Gaus(250.,25.);
    Float_t vals[npars]={dists,sigverts,cpas};
    Float_t valb[npars]={distb,sigvertb,cpab};
    mvsig->FillAndIntegrate(vals,ptbin);
    mvsig2->Fill(vals,ptbin);  // alternative way of filling
    mvbkg->FillAndIntegrate(valb,ptbin);
  }
  mvsig2->Integrate(); // mvsig2 is now equal to mvsig
  
  // Merge Pt bins
  AliMultiDimVector* mvsigallpt=mvsig->ShrinkPtBins(0,3);
  AliMultiDimVector* mvbkgallpt=mvbkg->ShrinkPtBins(0,3);
  
  //caluclate significance
  AliSignificanceCalculator* cal=new AliSignificanceCalculator(mvsigallpt,mvbkgallpt,1.,5.);
  AliMultiDimVector* mvsts=cal->GetSignificance();
  AliMultiDimVector* mvess=cal->GetSignificanceError();
  AliMultiDimVector* mvsob=cal->CalculateSOverB();
  AliMultiDimVector* mvesob=cal->CalculateSOverBError();
  AliMultiDimVector* mvpur=cal->CalculatePurity();
  AliMultiDimVector* mvepur=cal->CalculatePurityError();
  Int_t fixed[3]={0,0,0};
  gStyle->SetPalette(1);

  Int_t maxInd[3];
  Float_t sigMax=cal->GetMaxSignificance(maxInd,0);
  Float_t cut0=mvsigallpt->GetCutValue(0,maxInd[0]);
  Float_t cut1=mvsigallpt->GetCutValue(1,maxInd[1]);
  Float_t cut2=mvsigallpt->GetCutValue(2,maxInd[2]);

  printf("=========== Pt Integrated ==============\n");
  printf("Maximum of significance found in bin %d %d %d\n",maxInd[0],maxInd[1],maxInd[2]);
  printf("                                 cuts= %f %f %f\n",cut0,cut1,cut2);
  printf("Significance = %f +- %f\n",mvsts->GetElement(maxInd,0),mvess->GetElement(maxInd,0));
  printf("S/B          = %f +- %f\n",mvsob->GetElement(maxInd,0),mvesob->GetElement(maxInd,0));
  printf("Purity       = %f +- %f\n",mvpur->GetElement(maxInd,0),mvepur->GetElement(maxInd,0));


  // 2D projections
  TH2F* hsig1 = mvsigallpt->Project(0,2,fixed,0);
  TH2F* hbkg1 = mvbkgallpt->Project(0,2,fixed,0);
  TH2F* hsts1 = mvsts->Project(0,2,fixed,0);
  TH2F* hess1 = mvess->Project(0,2,fixed,0);
  TH2F* hpur1 = mvpur->Project(0,2,fixed,0);
  TH2F* hepur1 = mvepur->Project(0,2,fixed,0);
  TH2F* hsob1 = mvsob->Project(0,2,fixed,0);
  TH2F* hesob1 = mvesob->Project(0,2,fixed,0);

  TH2F* hsig2 = mvsigallpt->Project(1,2,fixed,0);
  TH2F* hbkg2 = mvbkgallpt->Project(1,2,fixed,0);
  TH2F* hsts2 = mvsts->Project(1,2,fixed,0);
  TH2F* hess2 = mvess->Project(1,2,fixed,0);
  TH2F* hpur2 = mvpur->Project(1,2,fixed,0);
  TH2F* hepur2 = mvepur->Project(1,2,fixed,0);
  TH2F* hsob2 = mvsob->Project(1,2,fixed,0);
  TH2F* hesob2 = mvesob->Project(1,2,fixed,0);

  TCanvas* c1=new TCanvas("c1","Var 0 vs. 2",1000,800);
  c1->Divide(4,2);
  c1->cd(1);
  hsig1->Draw("colz");
  c1->cd(2);
  hbkg1->Draw("colz");
  c1->cd(3);
  hsts1->Draw("colz");
  c1->cd(4);
  hess1->Draw("colz");
  c1->cd(5);
  hsob1->Draw("colz");
  c1->cd(6);
  hesob1->Draw("colz");
  c1->cd(7);
  hpur1->Draw("colz");
  c1->cd(8);
  hepur1->Draw("colz");

  TCanvas* c2=new TCanvas("c2","Var 1 vs. 2",1000,800);
  c2->Divide(4,2);
  c2->cd(1);
  hsig2->Draw("colz");
  c2->cd(2);
  hbkg2->Draw("colz");
  c2->cd(3);
  hsts2->Draw("colz");
  c2->cd(4);
  hess2->Draw("colz");
  c2->cd(5);
  hsob2->Draw("colz");
  c2->cd(6);
  hesob2->Draw("colz");
  c2->cd(7);
  hpur2->Draw("colz");
  c2->cd(8);
  hepur2->Draw("colz");

}
