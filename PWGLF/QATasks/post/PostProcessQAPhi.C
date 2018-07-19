void PostProcessQAPhi(char* name_fin="LHC17d1_B0p2.root",char* name_fout="QAphi_LHC17d1_B0p2",char* name_list="RsnQA_phi",char* name_hist_base="taskRsnQA_phi_",Bool_t* MC=kTRUE){

  //This macro was written by Anders Knospe (anders.knospe@cern.ch).
  //5 November 2013
 //Modified by Sourav Kundu (s.kundu@cern.ch).
  //4 october 2017

  //FOR DETAILED INSTRUCTIONS, SEE THE END OF THIS DOCUMENT.

  //Arguments:
  //name_fin: name of input file
  //name_fout: base name of output files (without suffix)
  //name_list: name of the list in fin that contains histograms, may be different depending on the cuts applied to the kaon daughters
  //name_hist_base: base name of the THnSparse histograms, may be different depending on the cuts applied to the kaon daughters
  //Bool_t* MC=kTRUE when we analyzed MC output and Bool_t* MC=kFALSE when we analyzed data
  
  //This Macro does the following:
  //1.) plots the standard histograms showing the number of events that pass selection of criteria and the number of events in each multiplicity or centrality bin
  
  //2.)Plots Accepted events vs. Multiplicity percentile

  //3.)Plots Resolution distribution, acceptance X efficiency vss. eta and pt, acceptance X efficiency vs, pT, acceptance X efficiency vs, eta, eta distribution and pT distribution for MC production.

  //4.) plots invariant-mass histograms (unlike-charge, like-charge, background-subtracted)

  //5.) fits the background-subtracted invariant-mass distribution.  In case there is a problem with the "default" fit, it does multiple fits using different residual backgrounds and fit regions.  All fits are plotted and all fit results are saved.

  //6.) prints out relevant information: pT range, phi yield per event, phi mass, and phi width (taken from the default fit)

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  //----- get input histograms -----

  TFile* fin=TFile::Open(name_fin);//open input file
  if(!fin) return;

  TList* l=(TList*) fin->Get(name_list);
  if(!l){cerr<<"Error in QAphi(): missing input list "<<name_list<<" in file "<<name_fin<<".  Stopping."<<endl; return;}

  TH1D *hevent,*hmult,*hu,*hm,*hp,*hl,*hs,*hs_plot,*hresolution,*efyeta,*efypt;
  TH2F *su,*sm,*sp,*resolution,*effypt,*effyeta,*gneta,*trueeta,*genpt,*truept,*ratio;
  float A,B,UA,UB;

  int lowptbin=3;
  int highptbin=4;



  hevent=(TH1D*) l->FindObject("hEventStat");//standard histogram in resonance package, shows number of events after various selections
  if(!hevent){cerr<<"Error in QAphi(): missing input histogram hEventStat in file "<<name_fin<<".  Stopping."<<endl; return;}

  hmult=(TH1D*) l->FindObject("hAEventsVsMulti");//standard histogram in resonance package, shows number of events in each multiplicity (for pp) or centrality (for Pb-Pb) bin
  if(!hmult){cerr<<"Error in QAphi(): missing input histogram hAEventsVsMulti in file "<<name_fin<<".  Stopping."<<endl; return;}
  double nevt=hmult->Integral();

  su=(TH2F*) l->FindObject(Form("%sUnlike_MPt",name_hist_base));//invariant-mass histogram, unlike-charge
  if(!su){cerr<<"Error in QAphi(): missing input histogram "<<name_hist_base<<"_Unlike in file "<<name_fin<<".  Stopping."<<endl; return;}

  sm=(TH2F*) l->FindObject(Form("%sLikeMM_MPt",name_hist_base));//invariant-mass histogram, like-charge (K-K-)
  if(!sm){cerr<<"Error in QAphi(): missing input histogram "<<name_hist_base<<"_LikeMM in file "<<name_fin<<".  Stopping."<<endl; return;}

  sp=(TH2F*) l->FindObject(Form("%sLikePP_MPt",name_hist_base));//invariant-mass histogram, like-charge (K+K+)
  if(!sp){cerr<<"Error in QAphi(): missing input histogram "<<name_hist_base<<"_LikePP in file "<<name_fin<<".  Stopping."<<endl; return;}

 if(MC)
   {
  
  resolution=(TH2F*) l->FindObject("taskRsnQA_phi_Res_ResPt");//resolution_pT
  genpt=(TH2F*) l->FindObject("taskRsnQA_phi_Gen_MPt");//resolution_pT
  recpt=(TH2F*) l->FindObject("taskRsnQA_phi_Trues_MPt");//resolution_pT
  
  genetapt=(TH2F*) l->FindObject("taskRsnQA_phi_Gen_EtaPt");//resolution_pT
  recetapt=(TH2F*) l->FindObject("taskRsnQA_phi_Trues_EtaPt");//resolution_pT
   }

  TFile* fout=new TFile(Form("%s.root",name_fout),"RECREATE","HistoFile");//open output file
 

  // TFile* fouteff=new TFile(Form("%seff.root",name_fout),"RECREATE","HistoFile");

  
  
  float dy=1.0;
  
  TCanvas* c=new TCanvas("c","",10,10,1500,500);
  c->SetFillColor(0);

  TPad* p1=new TPad("p1_0","",0.,0.,1./3,1.);
  p1->SetFillColor(0);
  
  TPad* p2=new TPad("p2_0","",1./3,0.,2./3,1.);
  p2->SetFillColor(0);

  TPad* p3=new TPad("p3_0","",2./3,0.,1.,1.);
  p3->SetFillColor(0);
  p3->SetLogy();

  p1->cd();
  hevent->SetLineColor(1);
  hevent->Draw();

  p2->cd();
  hmult->SetLineColor(1);
  int j,k,xmax=0;
  for(j=1;j<=hmult->GetNbinsX();j++) if(hmult->GetBinContent(j)>0.) xmax=j;
  xmax=(int) (1.1*xmax);
  if(xmax>hmult->GetNbinsX()) xmax=hmult->GetNbinsX();
  hmult->GetXaxis()->SetRange(1,xmax);
  hmult->Draw();
  // hmult->GetYaxis()->SetTitle("Events");
  hmult->GetXaxis()->SetTitle("multiplicity%");
  p3->cd();
  hmult->Draw();//same as p2, but log scale on y axis

  c->cd();
  p1->Draw();
  p2->Draw();
  p3->Draw();

  c->SaveAs(Form("%s.pdf(",name_fout));

  delete p1;
  delete p2;
  delete p3;
  delete c;
  
  fout->cd();
  hevent->Write();
  hmult->Write();




  if(MC)
    {
  TCanvas* c=new TCanvas("c","",10,10,1500,500);
  c->SetFillColor(0);
  
  TPad* p1=new TPad("p1_0","",0.,0.,1./3,1.);
  p1->SetFillColor(0);
  
  TPad* p2=new TPad("p2_0","",1./3,0.,2./3,1.);
  p2->SetFillColor(0);
  
  TPad* p3=new TPad("p3_0","",2./3,0.,1.,1.);
  p3->SetFillColor(0);

  p1->cd();
  TProfile *px = resolution->ProfileY("resolution",1,200,"e");
  px->SetMarkerColor(2);
  px->SetMarkerStyle(20);
  px->SetTitle("Resolution vs. #it{p}_{T} (GeV/#it{c})");
  px->SetXTitle("#it{p}_{T} (GeV/c)");
  px->Draw();
  
  p2->cd();


  TH2F *genetaptclone=(TH2F*)genetapt->Clone();
  TH2F *recetaptclone=(TH2F*)recetapt->Clone();
  recetaptclone->Divide(genetaptclone);
  recetaptclone->SetTitle("Acceptance #times Efficiency vs. #eta and #it{p}_{T}");
  recetaptclone->SetMarkerColor(2);
  recetaptclone->SetMarkerStyle(24);
  recetaptclone->GetXaxis()->SetTitle("#eta");
  recetaptclone->GetXaxis()->SetTitleSize(0.05);
  recetaptclone->GetXaxis()->SetTitleOffset(1.3);
  recetaptclone->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  recetaptclone->GetYaxis()->SetTitleSize(0.05);
  recetaptclone->GetYaxis()->SetTitleOffset(1.4);
  recetaptclone->Draw("surf");


  p3->cd();
  p3->SetLogy();
  TH1D* geneta2 = (TH1D*)genetapt->ProjectionX("geneta2",1,100,"E");
  TH1D* receta2 = (TH1D*)recetapt->ProjectionX("effyvseta",1,100,"E");
  TH1D* geneta2clone = (TH1D*) geneta2->Clone("geneta");
  TH1D* receta2clone = (TH1D*) receta2->Clone("receta");
  geneta2clone->SetMarkerColor(2);
  geneta2clone->SetMarkerStyle(20);
  receta2clone->SetMarkerColor(4);
  receta2clone->SetMarkerStyle(20);
  geneta2clone->SetTitle("#eta distribution");
  geneta2clone->SetXTitle("#eta");
  geneta2clone->GetYaxis()->SetRangeUser(2,1000000);
  geneta2clone->Draw("EP");
  receta2clone->Draw("same");
  TLegend* L=new TLegend(0.25,0.2,0.5,0.3);
  L->SetFillColor(0);
  L->AddEntry(geneta2clone,"genarated #phi","p");
  L->AddEntry(receta2clone,"reconstructed #phi","p");
  L->Draw();



  c->cd();
  p1->Draw();
  p2->Draw();
  p3->Draw();
  c->SaveAs(Form("%s.pdf(",name_fout));

  delete p1;
  delete p2;
  delete p3;
  delete c;
  
  
  TCanvas* c=new TCanvas("c","",10,10,1500,500);
  c->SetFillColor(0);
  
  TPad* p1=new TPad("p1_0","",0.,0.,1./4,1.);
  p1->SetFillColor(0);
  
  TPad* p2=new TPad("p2_0","",1./4,0.,2./4,1.);
  p2->SetFillColor(0);
  
  TPad* p3=new TPad("p3_0","",2./4,0.,3./4,1.);
  p3->SetFillColor(0);

  TPad* p4=new TPad("p3_0","",3./4,0.,1.,1.);
  p4->SetFillColor(0);


  p1->cd();
  receta2->Divide(geneta2);
  receta2->GetYaxis()->SetRangeUser(0.0,0.3);
  receta2->SetMarkerColor(2);
  receta2->SetMarkerStyle(20);
  receta2->SetTitle("Acceptance #times Efficiency vs. #eta");
  receta2->SetXTitle("#eta");
  receta2->Draw();
  
  p2->cd();
  p2->SetLogy();
  TH1D* genpt2 = (TH1D*)genetapt->ProjectionY("genpt2",1,100,"E");
  TH1D* recpt2 = (TH1D*)recetapt->ProjectionY("effyvspt",1,100,"E");
  TH1D* genpt2clone = (TH1D*) genpt2->Clone("genpt");
  TH1D* recpt2clone = (TH1D*) recpt2->Clone("recpt");
  
  genpt2clone->Scale(1/hevent->Integral(1,100));
  recpt2clone->Scale(1/hevent->Integral(1,100));
  genpt2clone->SetMarkerColor(2);
  genpt2clone->SetMarkerStyle(20);
  recpt2clone->SetMarkerColor(4);
  recpt2clone->SetMarkerStyle(20);
  genpt2clone->SetTitle("#it{p}_{T} distribution");
  genpt2clone->SetXTitle("#it{p}_{T} (GeV/c)");
  genpt2clone->Draw("EP");
  recpt2clone->Draw("Same");

  TLegend* L1=new TLegend(0.15,0.2,0.4,0.3);
  L1->SetFillColor(0);
  L1->AddEntry(genpt2clone,"genarated #phi","p");
  L1->AddEntry(recpt2clone,"reconstructed #phi","p");
  L1->Draw();

  p3->cd();
  recpt2->Rebin(2);
  genpt2->Rebin(2);
  recpt2->Divide(genpt2);
  recpt2->GetYaxis()->SetRangeUser(0.01,1.2);
  recpt2->SetMarkerColor(2);
  recpt2->SetMarkerStyle(20);
  recpt2->SetTitle("Acceptance #times Efficiency vs. #it{p}_{T}");
  recpt2->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  recpt2->Draw();


  p4->cd();
  TH1D* genpt22 = (TH1D*)genetapt->ProjectionY("genpt22",6,15,"E");
  TH1D* recpt22 = (TH1D*)recetapt->ProjectionY("effyvsptforetalessthanpoint5",6,15,"E");
  recpt22->Rebin(2);
  genpt22->Rebin(2);
  recpt22->Divide(genpt22);
  recpt22->GetYaxis()->SetRangeUser(0.01,1.2);
  recpt22->SetMarkerColor(2);
  recpt22->SetMarkerStyle(20);
  recpt22->SetTitle("Acceptance #times Efficiency vs. #it{p}_{T} (|#eta|#leq0.5)");
  recpt22->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  recpt22->Draw();


  c->cd();
  p1->Draw();
  p2->Draw();
  p3->Draw();
  p4->Draw();
  c->SaveAs(Form("%s.pdf(",name_fout));
  
  delete p1;
  delete p2;
  delete p3;
  delete c;

  
  fout->cd();
  genpt2clone->Write();
  recpt2clone->Write();
  receta2->Write();
  recpt2->Write();
  px->Write();
  
    }











  //----- get invariant-mass histograms -----

  //Project the TH2F histograms onto their x axes.
  hu=(TH1D*) su->ProjectionX("su",lowptbin,highptbin,"E");
  hu->SetName("mass_unlike");
  hm=(TH1D*) sm->ProjectionX("sm",lowptbin,highptbin,"E");
  hm->SetName("mass_likeMM");
  hp=(TH1D*) sp->ProjectionX("sp",lowptbin,highptbin,"E");
  hp->SetName("mass_likePP");

  float lowpt=(lowptbin-1)*0.2;
  float highpt=highptbin*0.2;

 hl=(TH1D*) hm->Clone("mass_like");//total like-charge histogram
  for(j=1;j<=hl->GetNbinsX();j++){
    A=hm->GetBinContent(j);
    B=hp->GetBinContent(j);
    hl->SetBinContent(j,2*sqrt(A*B));//note: 2sqrt(n{K-K-}*n{K+K+})
    hl->SetBinError(j,sqrt(A+B));
  }
  hl->SetTitle("Like-Charge Combinatorial Background");

  hs=(TH1D*) hu->Clone("mass_signal");

    for(j=1;j<=hs->GetNbinsX();j++){
    A=hu->GetBinContent(j); UA=hu->GetBinError(j);
    B=hl->GetBinContent(j); UB=hl->GetBinError(j);
    hs->SetBinContent(j,A-B);//subtract the combinatorial background
    hs->SetBinError(j,sqrt(UA*UA+UB*UB));
  }
  
  double Imin=1.01,Imax=1.03,dx=hs->GetXaxis()->GetBinWidth(1);

  hb=(TH1D*) hs->Clone("mass_background");//temporary histogram with the peak removed; used to fit the residual background
  for(j=1;j<=hb->GetNbinsX();j++){
    A=hb->GetXaxis()->GetBinCenter(j);
    if(A>Imin && A<Imax){hb->SetBinContent(j,0); hb->SetBinError(j,0);}
  }

  hs_plot=(TH1D*) hs->Clone("mass_signal_plot");//temporary histogram for plotting

  //----- plot the invariant-mass histograms -----

  c=new TCanvas("c","",10,10,1000,500);
  c->SetFillColor(0);

  p1=new TPad("p1_1","",0.,0.,0.5,1.);
  p1->SetFillColor(0);
  p1->SetRightMargin(0.05);

  p2=new TPad("p2_1","",0.5,0.,1.,1.);
  p2->SetFillColor(0);
  p2->SetRightMargin(0.05);

  p1->cd();
  hu->SetLineColor(1); hu->SetMarkerColor(1);
  hl->SetLineColor(2); hl->SetMarkerColor(2);
  hu->SetTitle("Invariant Mass (KK)     ");
  hu->SetXTitle("Invariant Mass (KK) [GeV/#it{c}^{2}]");
  hu->Draw();
  hl->Draw("same");

  TLatex* t1=new TLatex(0.7,0.96,"Unlike-Charge");
  t1->SetTextSize(0.04);
  t1->SetNDC();
  t1->Draw();

  TLatex* t2=new TLatex(0.7,0.91,"Like-Charge");
  t2->SetTextSize(0.04);
  t2->SetNDC();
  t2->SetTextColor(2);
  t2->Draw();

  p2->cd();
  hs->SetLineColor(1); hs->SetMarkerColor(1);
  hs->SetTitle("Combinatorial Background Subtracted");
  hs->SetXTitle("Invariant Mass (KK) [GeV/#it{c}^{2}]");
  hs->Draw();

  c->cd();
  p1->Draw();
  p2->Draw();

  c->SaveAs(Form("%s.pdf",name_fout));

  delete p1;
  delete p2;
  delete c;

  fout->cd();
  hu->Write();
  hl->Write();
  hs->Write();


  c=new TCanvas("c","",10,10,500,500);
  c->SetFillColor(0);

  //fit peak and extract yield, mass, and width
  //The peak is fit using a polynomial for the residual background, plus a Voigtian peak.  A Voigtian peak is the convolution of a (non-relativistic) Breit-Wigner peak and a Gaussian (to describe detector effects).  The resolution (Gaussian sigma) of the Voigtian peak will be fixed.
  //Ideally, we would only need to fit once.  However, the fits are not always stable so we try a variety of fits and record all results.  A similar procedure was followed in the analysis of phi mesons in Pb+Pb collisions (2010 data).

  int jBF;//index controlling the order of the residual background polynomial
  int nBF=3;//number of different orders used for the residual background polynomial; currently 3 (linear, quadratic, cubic)
  int jFR;//index controlling the fit region
  int nFR=4;//number of different fit regions
  double urb;

  //This TH2D is just an array of numbers to hold the different fit results and associated quantities.  The x-axis gives the name of the fit (order of residual background and number of the fit region), the y-axis gives the name of the quantity stored.
  TH2D* r=new TH2D("fit_results","",nBF*nFR,0,nBF*nFR, 16,0,16);
  r->GetYaxis()->SetBinLabel(1,"peak integral");
  r->GetYaxis()->SetBinLabel(2,"mass");
  r->GetYaxis()->SetBinLabel(3,"width");
  r->GetYaxis()->SetBinLabel(4,"resolution (fixed)");
  r->GetYaxis()->SetBinLabel(5,"p0");
  r->GetYaxis()->SetBinLabel(6,"p1");
  r->GetYaxis()->SetBinLabel(7,"p2");
  r->GetYaxis()->SetBinLabel(8,"p3");
  r->GetYaxis()->SetBinLabel(9,"chi2");
  r->GetYaxis()->SetBinLabel(10,"NDF");
  r->GetYaxis()->SetBinLabel(11,"yield (fit)");
  r->GetYaxis()->SetBinLabel(12,"yield/event (fit)");
  r->GetYaxis()->SetBinLabel(13,"peak correction factor");
  r->GetYaxis()->SetBinLabel(14,"yield (bin counting)");
  r->GetYaxis()->SetBinLabel(15,"yield/event (bin conunting)");
  r->GetYaxis()->SetBinLabel(16,"yield/event (bin conunting + tails)");

  TF1 *g[3][4],*gp[3][4],*gb[3][4];
  TFitResultPtr fr;
  int status;
  double pc;

  cerr<<"Please wait a moment while the peak fits are performed."<<endl;

  for(jBF=0;jBF<nBF;jBF++) for(jFR=0;jFR<nFR;jFR++){//loops: change the order of the residual background polynomial and the fit region
      if(!jFR){A=0.995; B=1.06;}//fr0
      else if(jFR==1){A=1.; B=1.06;}//fr1
      else if(jFR==2){A=0.995; B=1.07;}//fr2
      else if(jFR==3){A=1.; B=1.07;}//fr3

      if(!jBF){//pol1
	g[jBF][jFR]=new TF1(Form("fit_pol1_fr%i",jFR),"[0]*TMath::Voigt(x-[1],[3],[2])+[4]+[5]*x",A,B);
	gb[jBF][jFR]=new TF1(Form("fit_back_pol1_fr%i",jFR),"pol1",A,B);
	gp[jBF][jFR]=new TF1(Form("fit_peak_pol1_fr%i",jFR),"[0]*TMath::Voigt(x-[1],[3],[2])",A,B);
      }else if(jBF==1){//pol2
	g[jBF][jFR]=new TF1(Form("fit_pol2_fr%i",jFR),"[0]*TMath::Voigt(x-[1],[3],[2])+[4]+[5]*x+[6]*x*x",A,B);
	gb[jBF][jFR]=new TF1(Form("fit_back_pol2_fr%i",jFR),"pol2",A,B);
	gp[jBF][jFR]=new TF1(Form("fit_peak_pol2_fr%i",jFR),"[0]*TMath::Voigt(x-[1],[3],[2])",A,B);
      }else if(jBF==2){//pol3
	g[jBF][jFR]=new TF1(Form("fit_pol3_fr%i",jFR),"[0]*TMath::Voigt(x-[1],[3],[2])+[4]+[5]*x+[6]*x*x+[7]*x*x*x",A,B);
	gb[jBF][jFR]=new TF1(Form("fit_back_pol3_fr%i",jFR),"pol3",A,B);
	gp[jBF][jFR]=new TF1(Form("fit_peak_pol3_fr%i",jFR),"[0]*TMath::Voigt(x-[1],[3],[2])",A,B);
      }

      for(j=0;j<100;j++){//initial fit of residual background
	fr=hb->Fit(gb[jBF][jFR],"RSQ");
	status=fr->Status();
	if(!status) break;
      }
      urb=gb[jBF][jFR]->IntegralError(Imin,Imax)/dx;//estimated uncertainty of residual background integral

      j=hs->GetXaxis()->FindBin(Imin+1e-5); k=hs->GetXaxis()->FindBin(Imax-1e-5);
      g[jBF][jFR]->SetParameter(0,hs->Integral(j,k)*dx);//initial value of peak integral
      g[jBF][jFR]->SetParameter(1,1.019455); g[jBF][jFR]->FixParameter(1,1.019455);//fix mass to vacuum value
      g[jBF][jFR]->SetParameter(2,0.00426); g[jBF][jFR]->FixParameter(2,0.00426);//fix width to vacuum value
      g[jBF][jFR]->SetParameter(3,0.0011); g[jBF][jFR]->FixParameter(3,0.0011);//fix resolution to 1.1 MeV/c^2 (will not be changed)
      for(j=0;j<=jBF+1;j++){//fix residual background to the initial fit
	A=gb[jBF][jFR]->GetParameter(j);
	g[jBF][jFR]->SetParameter(j+4,A);
	g[jBF][jFR]->SetParError(j+4,gb[jBF][jFR]->GetParError(j));
	g[jBF][jFR]->FixParameter(j+4,A);
      }

      cerr<<"Fitting pol"<<jBF+1<<"_fr"<<jFR<<endl;

      for(j=0;j<100;j++){//fit with only the peak integral free
	fr=hs->Fit(g[jBF][jFR],"RSQ");
	status=fr->Status();
	if(!status) break;
      }

      g[jBF][jFR]->ReleaseParameter(2);//release width parameter
      for(j=0;j<100;j++){
	fr=hs->Fit(g[jBF][jFR],"RSQ");
	status=fr->Status();
	if(!status) break;
      }

      g[jBF][jFR]->ReleaseParameter(1);//release mass parameter
      for(j=0;j<100;j++){
	fr=hs->Fit(g[jBF][jFR],"RSQ");
	status=fr->Status();
	if(!status) break;
      }

      //If we want to use resolution as free parameter then we need to uncomment the below section//
      
      /*
 g[jBF][jFR]->ReleaseParameter(3);//release resolution parameter
      for(j=0;j<100;j++){
	fr=hs->Fit(g[jBF][jFR],"RSQ");
	status=fr->Status();
	if(!status) break;
      }

      */

      
      g[jBF][jFR]->ReleaseParameter(4);//release residual background constant parameter
      for(j=0;j<100;j++){
	fr=hs->Fit(g[jBF][jFR],"RSQ");
	status=fr->Status();
	if(!status) break;
      }

      for(j=1;j<=gb[jBF][jFR]->GetNpar();j++) g[jBF][jFR]->ReleaseParameter(4+j);//release other residual background parameters
      for(j=0;j<100;j++){
	fr=hs->Fit(g[jBF][jFR],"RSQ");
	status=fr->Status();
	if(!status) break;
      }

      for(j=0;j<100;j++){//redo fit using "I" option, typically does not affect the result very much, but best to do it anyway.  Comment this out if you are debugging and want to save time.
	fr=hs->Fit(g[jBF][jFR],"RSQI");
	status=fr->Status();
	if(!status) break;
      }

      fout->cd();
      g[jBF][jFR]->Write();//save the combined fit

      for(j=0;j<gb[jBF][jFR]->GetNpar();j++) gb[jBF][jFR]->SetParameter(j,g[jBF][jFR]->GetParameter(j+4));//put the background parameters into the function gb
      for(j=0;j<4;j++) gp[jBF][jFR]->SetParameter(j,g[jBF][jFR]->GetParameter(j));//put the peak parameters into the funciton gp

      j=nFR*jBF+jFR+1;//bin for the x-axis of results histogram
      r->GetXaxis()->SetBinLabel(j,Form("pol%i_fr%i",jBF+1,jFR));
      for(k=0;k<g[jBF][jFR]->GetNpar();k++){//store the values of the peak parameters
	r->SetBinContent(j,k+1,g[jBF][jFR]->GetParameter(k));
	r->SetBinError(j,k+1,g[jBF][jFR]->GetParError(k));
      }
      r->SetBinContent(j,9,g[jBF][jFR]->GetChisquare());//store chi^2
      r->SetBinContent(j,10,g[jBF][jFR]->GetNDF());//store number of degrees of freedom

      A=r->GetBinContent(j,1)-gp[jBF][jFR]->Integral(0.,2*0.493667);//subtract integral of peak below kinematic cutoff
      UA=A*r->GetBinError(j,1)/r->GetBinContent(j,1);
      r->SetBinContent(j,11,A/dx); r->SetBinError(j,11,UA/dx);//peak yield extracted from the fit
      r->SetBinContent(j,12,A/dx/nevt/dy); r->SetBinError(j,12,UA/dx/nevt/dy);//peak yield per event extracted from the fit, corrected for dy

      B=gp[jBF][jFR]->Integral(Imin,Imax);
      pc=B/A;//peak correction factor: the yield in the interval (Imin,Imax) divided by the total integral above the kinematic cutoff, used to correct the yield from bin counting to account for the yield in the tails
      r->SetBinContent(j,13,pc);//store peak correction factor

      A=UA=0;
      for(k=hs->GetXaxis()->FindBin(Imin+1e-5);k<=hs->GetXaxis()->FindBin(Imax-1e-5);k++){//yield from bin counting
	A+=hs->GetBinContent(k);
	UA+=pow(hs->GetBinError(k),2);
      }
      A-=gb[jBF][jFR]->Integral(Imin,Imax)/dx;//subtract residual background integral
      UA+=urb*urb;

      r->SetBinContent(j,14,A); r->SetBinError(j,14,sqrt(UA));//yield from bin counting
      r->SetBinContent(j,15,A/nevt/dy); r->SetBinError(j,15,sqrt(UA)/nevt/dy);//yield per event from bin counting, corrected for dy
      r->SetBinContent(j,16,A/nevt/pc/dy); r->SetBinError(j,16,sqrt(UA)/nevt/pc/dy);//yield per event from bin counting, corrected for dy, corrected to account for yield in tails
    }

  fout->cd();
  r->Write();//save the results histogram in case it is needed later

  delete c;

  //----- plot the fits -----

  c=new TCanvas("c","",10,10,1000,750);
  c->SetFillColor(0);

  p1=new TPad("p1_2","",0,0.5,0.5,1);
  p1->SetFillColor(0);
  p1->SetRightMargin(0.05);

  p2=new TPad("p2_2","",0.5,0.5,1,1);
  p2->SetFillColor(0);
  p2->SetRightMargin(0.05);

  p3=new TPad("p3_2","",0,0,0.5,0.5);
  p3->SetFillColor(0);
  p3->SetRightMargin(0.05);

  TPad* p4=new TPad("p4_2","",0.5,0,1,0.5);
  p4->SetFillColor(0);

  hs_plot->GetXaxis()->SetRangeUser(0.995,1.07-1e-5);
  hs_plot->SetLineColor(1); hs_plot->SetMarkerColor(1);
  hs_plot->SetTitle(hs->GetTitle());
  hs_plot->SetXTitle(hs->GetXaxis()->GetTitle());

  for(jBF=0;jBF<nBF;jBF++) for(jFR=0;jFR<nFR;jFR++){
      //assign styles and colors to the functions
      g[jBF][jFR]->SetNpx(300);
      if(!jFR) g[jBF][jFR]->SetLineColor(2);
      else if(jFR==1) g[jBF][jFR]->SetLineColor(TColor::GetColor("#ffaa00"));
      else if(jFR==2) g[jBF][jFR]->SetLineColor(TColor::GetColor("#238e23"));
      else if(jFR==3) g[jBF][jFR]->SetLineColor(4);
      gb[jBF][jFR]->SetLineColor(g[jBF][jFR]->GetLineColor());
      gb[jBF][jFR]->SetLineStyle(2);
    }

  p1->cd();
  hs_plot->Draw();
  for(jFR=0;jFR<nFR;jFR++) gb[0][jFR]->Draw("same");
  for(jFR=0;jFR<nFR;jFR++) g[0][jFR]->Draw("same");
  t1=new TLatex(0.94,0.87,"linear RB (pol1)"); t1->SetTextAlign(32);
  t1->SetNDC();
  t1->Draw();

  p2->cd();
  hs_plot->Draw();
  for(jFR=0;jFR<nFR;jFR++) gb[1][jFR]->Draw("same");
  for(jFR=0;jFR<nFR;jFR++) g[1][jFR]->Draw("same");
  t2=new TLatex(0.94,0.87,"quadratic RB (pol2)"); t2->SetTextAlign(32);
  t2->SetNDC();
  t2->Draw();

  p3->cd();
  hs_plot->Draw();
  for(jFR=0;jFR<nFR;jFR++) gb[2][jFR]->Draw("same");
  for(jFR=0;jFR<nFR;jFR++) g[2][jFR]->Draw("same");
  TLatex* t3=new TLatex(0.94,0.87,"cubic RB (pol3)"); t3->SetTextAlign(32);
  t3->SetNDC();
  t3->Draw();

//additional information in pad p4

  p4->cd();
  TLine* dummy1=new TLine(0,0,1,0);
  dummy1->SetLineColor(1); dummy1->SetLineWidth(2);
  TLine* dummy2=new TLine(0,0,1,0);
  dummy2->SetLineColor(1); dummy2->SetLineWidth(2); dummy2->SetLineStyle(gb[0][0]->GetLineStyle());
  TLegend* L=new TLegend(0.3,0.54,0.84,0.99);
  L->SetFillColor(0);
  L->AddEntry(g[0][0],"Fit Region 0 (0.995,1.06)","l");
  L->AddEntry(g[0][1],"Fit Region 1 (1,1.06)","l");
  L->AddEntry(g[0][2],"Fit Region 2 (0.995,1.07)","l");
  L->AddEntry(g[0][3],"Fit Region 3 (1,1.07)","l");
  L->AddEntry(dummy1,"Combined Fit","l");
  L->AddEntry(dummy2,"Residual Background Fit","l");
  L->AddEntry(dummy2,Form("%2.2f #leq #it{p}_{T} < %2.2f GeV/#it{c}",lowpt,highpt),"");
  L->Draw();




  double yield=r->GetBinContent(8,16);
  double mass=r->GetBinContent(8,2);
  double width=r->GetBinContent(8,3);
   double massresolution=r->GetBinContent(8,4);
   cout<<"resolution"<<massresolution<<"\n";

   
  double tx=0.005,ty=0.6,dty=0.06;
  TString s;

  TLatex* a1=new TLatex(tx,ty-2*dty,Form("#phi Yield/event = %1.5e #pm %1.5e (%s)",r->GetBinContent(8,16),r->GetBinError(8,16),r->GetXaxis()->GetBinLabel(8)));
  a1->Draw();

  TLatex* a2=new TLatex(tx,ty-3*dty,Form("#phi Mass = %1.5f #pm %1.5f GeV/#it{c}^{2}",r->GetBinContent(8,2),r->GetBinError(8,2)));
  a2->Draw();

  s.Form("PDG Mass = 1.019455");
  TLatex* a3=new TLatex(tx+0.1,ty-4*dty,s.Data());
  a3->Draw();

  TLatex* a4=new TLatex(tx,ty-5*dty,Form("#phi Width = %1.5f #pm %1.5f MeV/#it{c}^{2}",r->GetBinContent(8,3)*1000.,r->GetBinError(8,3)*1000.));
  a4->Draw();

  s.Form("PDG Width = 4.26");
  TLatex* a5=new TLatex(tx+0.1,ty-6*dty,s.Data());
  a5->Draw();

 s.Form("#phi mass resolution = 0.0011 (fixed during fitting)");
  TLatex* a6=new TLatex(tx+0.1,ty-7*dty,s.Data());
  a6->Draw();

 

  c->cd();
  p1->Draw();
  p2->Draw();
  p3->Draw();
  p4->Draw();

  c->SaveAs(Form("%s.pdf)",name_fout));



  //----- print out results -----

  j=su->GetXaxis()->GetFirst(); k=su->GetXaxis()->GetLast();
  printf("%2.2f < pT < %2.2f GeV/c, ",lowpt,highpt);
  j=su->GetYaxis()->GetFirst(); k=su->GetYaxis()->GetLast();
  printf("%2.2f < y < %2.2f\n",su->GetYaxis()->GetBinLowEdge(j),su->GetYaxis()->GetBinLowEdge(k+1));
  printf("phi Yield/event = %1.5e +/- %1.5e (%s)\n",r->GetBinContent(8,16),r->GetBinError(8,16),r->GetXaxis()->GetBinLabel(8));
  printf("phi Mass = %1.5f +/- %1.5f GeV/c^2 (PDG Value = 1.019455)\n",r->GetBinContent(8,2),r->GetBinError(8,2));
  printf("phi Width = %1.5f +/- %1.5f MeV/c^2 (PDG Value = 4.26)\n",r->GetBinContent(8,3)*1000.,r->GetBinError(8,3)*1000.);
 
  
  //----- finish -----

  cerr<<"\nMacro has finished successfully."<<endl;

  fin->Close();
  fout->Close();








  




  /*
  
 //----- check results -----

  double ex_yield=0;
  if(!strcmp(system,"pp276")) ex_yield=1.416e-3;
  else if(!strcmp(system,"pp7")) ex_yield=1.706e-3;
  else if(!strcmp(system,"PbPb276")) ex_yield=0.245;
  else{
    cerr<<"Warning in QAphi(): Unknown collision system "<<system<<".  Expected yield not known."<<endl;
    ex_yield=1e10;
  }
  double yield=r->GetBinContent(8,16);
  int status_yield;
  //if(yield/ex_yield<1.5 && yield/ex_yield>0.5) status_yield=0;//OK
  if(yield/ex_yield>0.0) status_yield=0;//OK
  //else if(yield/ex_yield<5 && yield/ex_yield>0.2) status_yield=1;//problem
  else if(yield/ex_yield<0.0) status_yield=1;//problem
  else status_yield=2;//big problem

  double ex_mass=1.019455;
  double mass=r->GetBinContent(8,2);
  int status_mass;
  if(fabs(mass-ex_mass)<1) status_mass=0;//OK
  else if(fabs(mass-ex_mass)<3) status_mass=1;//problem
  else status_mass=2;//big problem

  double ex_width=0.00426;
  double width=r->GetBinContent(8,3);
  int status_width;
  if(fabs(width-ex_width)<1) status_width=0;//OK
  else if(fabs(width-ex_width)<2) status_width=1;//problem
  else status_width=2;//big problem;

  double tx=0.005,ty=0.6,dty=0.06;

  TString s;
  bool ptOK=1;
  bool yOK=1;
  j=su->GetXaxis()->GetFirst(); k=su->GetXaxis()->GetLast();
  s.Form("%2.2f < #it{p}_{T} < %2.2f GeV/#it{c}, ",lowpt,highpt);
  j=su->GetYaxis()->GetFirst(); k=su->GetYaxis()->GetLast();
  s.Append(Form("%2.2f < #it{y} < %2.2f",su->GetYaxis()->GetBinLowEdge(j),su->GetYaxis()->GetBinLowEdge(k+1)));
  if(!ptOK || !yOK) s.Append(" [PROBLEM]");
  TLatex* a1=new TLatex(tx,ty,s.Data());
  SetText(a1,ptOK);
  // a1->Draw();
  
  TLatex* a2=new TLatex(tx,ty-2*dty,Form("#phi Yield/event = %1.5e #pm %1.5e (%s)",r->GetBinContent(8,16),r->GetBinError(8,16),r->GetXaxis()->GetBinLabel(8)));
  SetText(a2,status_yield);
  a2->Draw();

  s.Form("Expected Yield = %1.5e",ex_yield);
  if(status_yield) s.Append(" [PROBLEM]");
  TLatex* a3=new TLatex(tx+0.1,ty-2*dty,s.Data());
  SetText(a3,status_yield);
  // a3->Draw();

  TLatex* a4=new TLatex(tx,ty-3*dty,Form("#phi Mass = %1.5f #pm %1.5f GeV/#it{c}^{2}",r->GetBinContent(8,2),r->GetBinError(8,2)));
  SetText(a4,status_mass);
  a4->Draw();

  s.Form("PDG Mass = 1.019455");
  if(status_mass) s.Append(" [PROBLEM]");
  TLatex* a5=new TLatex(tx+0.1,ty-4*dty,s.Data());
  SetText(a5,status_mass);
  a5->Draw();

  TLatex* a6=new TLatex(tx,ty-5*dty,Form("#phi Width = %1.5f #pm %1.5f MeV/#it{c}^{2}",r->GetBinContent(8,3)*1000.,r->GetBinError(8,3)*1000.));
  SetText(a6,status_width);
  a6->Draw();

  s.Form("PDG Width = 4.26");
  if(status_width) s.Append(" [PROBLEM]");
  TLatex* a7=new TLatex(tx+0.1,ty-6*dty,s.Data());
  SetText(a7,status_width);
  a7->Draw();

  if(!ptOK || !yOK || status_yield || status_mass || status_width){
    status=2;
    s.Form("THERE ARE PROBLEMS!!!");
  }else{
    status=0;
    s.Form("Things look OK.");
  }

  TLatex* a8=new TLatex(tx,ty-7*dty,s.Data());
  SetText(a8,status); a8->SetTextSize(0.05);
  //a8->Draw();

  if(status) s.Form("PLEASE CHECK THE ERROR MESSAGES.");
  else s.Form("");

  TLatex* a9=new TLatex(tx,ty-8*dty,s.Data());
  SetText(a9,status); a9->SetTextSize(0.05);
  a9->Draw();

  c->cd();
  p1->Draw();
  p2->Draw();
  p3->Draw();
  p4->Draw();

  c->SaveAs(Form("%s.pdf)",name_fout));
  */
 
}


bool SetPtRange(TH2F* h){
  /*  TAxis* a=h->GetYaxis();
  double min=0.5,max=1.5;

  int b1=a->FindBin(1.00001*min);
  int b2=a->FindBin(0.99999*max);
  a->SetRange(b1,b2);
  b1=a->GetFirst();
  b2=a->GetLast();
  
  //min and/or max are not bin boundaries
  if(fabs(a->GetBinLowEdge(b1)-min)>1e-5*min || fabs(a->GetBinLowEdge(b2+1)-max)>1e-5*max){
    cerr<<"Error in QAphi(): pT range cannot be set to ("<<min<<","<<max<<").  The yield may therefore be different from the expected value."<<endl;
    return false;
  }
  */
  return true;
}

bool SetRapidityRange(TH2F* h,double& dy){
  /* TAxis* a=h->GetYaxis();
  double min=-0.5,max=0.5;

  int b1=a->FindBin(min+1e-5);
  int b2=a->FindBin(max-1e-5);
  a->SetRange(b1,b2);
  b1=a->GetFirst();
  b2=a->GetLast();
  dy=a->GetBinLowEdge(b2+1)-a->GetBinLowEdge(b1);

  //min and/or max are not bin boundaries
  if(fabs(a->GetBinLowEdge(b1)-min)>1e-5 || fabs(a->GetBinLowEdge(b2+1)-max)>1e-5){
    cerr<<"Error in QAphi(): rapidity range cannot be set to |y|<0.5.  Although this macro does correct for the rapidity range dy, using a different rapidity range may cause the yield to be different from the expected value."<<endl;
    return false;
  }
  */
  return true;
}


void SetText(TLatex* t,int flag){
  t->SetTextSize(0.04);
  if(!flag) t->SetTextColor(TColor::GetColor("#238e23"));
  else if(flag==1) t->SetTextColor(TColor::GetColor("#cc5500"));
  else t->SetTextColor(2);
  return;
}


void SetText(TLatex* t,bool flag){
  if(flag) SetText(t,0);
  else SetText(t,2);
  return;
}


/*-----------------------------------------------

INSTRUCTIONS:

*Introduction:

This macro is designed to read the output of the analysis task to check that the yield per event, mass, and width of the phi meson in the production are within acceptable ranges.  The analysis task constructs invariant-mass distributions of charged-kaon pairs in the vicinity of the mass of the phi meson (1.019455 GeV/c^2).  The branching ratio of the phi->K-K+ decay is 0.489.  In this macro, a combinatorial background is constructed using the like-charge distributions (K-K- and K+K+) and subtracted from the unlike-charge (K-K+) distribution.  This leaves a phi peak sitting on top of a residual background.  The background-subtracted distribution is fit with a Voigtian peak (a convolution of a non-relativistic Breit-Wigner peak with a Gaussian to account for detector effects) plus a polynomial to describe the residual background.  Because these fits sometimes fail, a total of 12 fits are performed.  There are 3 choices of residual background polynomial (linear, quadratic, or cubic) and 4 different fit regions (defined in the macro).  For more information on the fitting procedure, see refs. [2] and [3].  The yield, mass, and width of the phi meson extracted from the default fit (quadratic residual background, fit region 3) are compared to the expected values.

*Output:

This macro produces an output PDF file and an output ROOT file.

**Output PDF:

The output PDF consists of three canvases, organized into panels as follows:

Canvas 1: [[1a][1b][1c]]
Canvas 2: [[2a][2b]]
Canvas 3: [[3a][3b]]
          [[3c][3d]]

The contents of the panels is:
1a.) histogram hEventStat (see description below)
1b.) histogram hAEventsVsMulti (see description below)
1c.) histogram hAEventsVsMulti, logarithmic y-axis
2a.) unlike-charge invariant-mass distribution and like-charge combinatorial background
2b.) background-subtracted invariant-mass distribution
3a.) background-subtracted invariant-mass distribution with fits.  The fits shown here assume a linear residual background.  A fit is shown for each of the four different fit regions.  Fit region 3 (blue) is plotted on top.  The dashed lines are the residual background.
3b.) same as 3a, but with a quadratic residual background.  Note that the blue fit in this panel is the default fit, from which the yield, mass, and width are extracted.
3c.) same as 3a, but with a cubic residual background
3d.) a legend describing the fit functions and a printout of information.  The printed information should be green.  If it is not, there is a problem.

**Output ROOT file:

The output ROOT file contains the following:
1.) histogram hEventStat (see description below)
2.) histogram hAEventsVsMulti (see description below)
3.) histogram mass_unlike: the unlike-charge invariant-mass distribution
4.) histogram mass_like: the like-charge combinatorial background
5.) histogram mass_signal: background-subtracted invariant-mass distribution (mass_unlike - mass_like)
6.) functions fit_pol*_fr*: the fit functions for mass_signal.  pol1: linear residual background, pol2: quadratic residual background, pol3: cubic residual background.  fr* indicates the fit region.
7.) histogram fit_results:  contains the results of the different fits (see below)

***Histogram fit_results:

This is an array of numbers (with uncertainties) stored in a TH2D.  The x-axis is the name of the fit (e.g., "pol2_fr3") and the y-axis gives the name of the quantity stored.  The quantities stored (numbered along the y-axis) are:

1.) peak integral: integral of the phi peak in the fit function (parameter [0])
2.) peak mass: parameter [1]
3.) peak width: parameter [2]
4.) peak resolution: parameter [3]: currently fixed to 1.1 MeV/c^2, but stored anyway
5-8.) coefficients of the residual background polynomial (parameters [4] through [7])
9.) fit chi^2
10.) fit number of degrees of freedom
11.) yield extracted from the fit: parameter [0] minus the yield below the kinematic cutoff and corrected for the width of the invariant-mass bins
12.) yield/event (fit): the value above divided by the number of events and the width of the rapidity range (dy)
13.) peak correction factor: the fraction of the peak that lies inside the range (1.01,1.03) GeV/c^2
14.) yield (bin counting): the integral of the mass_signal histogram for (1.01,1.03) minus the integral of the residual background over the same interval
15.) yield/event (bin counting) the value above divided by the number of events and the width of the rapidity range (dy)
16.) yield/event (bin counting + tails) the value above further corrected by dividing by the peak correction factor (now accounting for the yield in the tails)

The default fit is "pol2_fr3".  Quantity 2 is reported as the mass, quantity 3 is reported as the width, and quantity 16 is reported as the yield.  So, to find the default mass (which is printed out) do fit_results->GetBinContent(8,2).  The default width is fit_results->GetBinContent(8,3) and the default yield is fit_results->GetBinContent(8,16).

*More Information:

In the analysis task, the kaons should be selected using track selection cuts and a 3-sigma cut (about the kaon mean) on TPC dE/dx.  See the analysis task for the exact cuts.

The output of the analysis task is expected to contain:
1.) the histogram hEventStat, which gives the number of events that pass different selection criteria
2.) the histogram hAEventsVsMulti, which gives the number of events in each multiplicity (for pp) or centrality (for Pb-Pb) bin
3.) three THnSparse histograms with three dimensions and the suffixes "Unlike", "LikePP", and "LikeMM"

For each THnSparse, axis 0 is the KK invariant mass, axis 1 is the transverse momentum, and axis 2 is the rapidity.  This macro attempts to select the range 0.5<pT<1.5 GeV/c; if this is not possible an error message will be generated, as the expected yields were calculated assuming that pT range.  This macro uses the full rapidity range on axis 2.  The expected yields were calculated assuming a rapidity range of |y|<0.5.  The macro corrects for the width of the rapidity range (dy), but using a different rapidity range may cause the yield to deviate from the expected value.

*Expected Values

**Expected Yields

If the yield is less than 50% or more than 150% of the expected value, it is flagged as an orange problem.  If the yield deviates from the expected value by more than a factor of 5, it is flagged as a red problem.

The expected yields differ depending on the collision system.

pp Collisions at 7 TeV (system="pp7"): The expected yield of 1.706e-3 is computed based on the published (ref. [1]) phi spectrum.  The yield per event for 0.5<pT<1.5 GeV/c is multiplied by the efficiency and the branching ratio.

pp Collisions at 2.76 TeV (system="pp276"): The expected yield of 1.416e-3 is computed by extrapolating from the pp 7 TeV yield, assuming that the phi yield scales as s^0.1 (energy dependence taken from ref. [2]).

Pb-Pb Collisions at 2.76 TeV (system="PbPb276"): The expected yield of 0.245 is computed based on the nearly published (ref. [2]) phi spectrum.

**Expected Mass

The expected phi mass is the PDG value: 1.019455 GeV/c^2.  A deviation from this value of more than 1 MeV/c^2 is flagged as an orange problem.  A deviation from this value of more than 3 MeV/c^2 is flagged as a red problem.

**Expected Width

The expected phi width is the PDG value: 4.26 MeV/c^2.  A deviation from this value of more than 1 MeV/c^2 is flagged as an orange problem.  A deviation from this value of more than 3 MeV/c^2 is flagged as a red problem.

*Troubleshooting

If a parameter does not have the expected value, it will be flagged.  In the output PDF file, a problematic parameter will be colored orange or red, with red indicating a larger deviation from the expected value.

If the yield, mass, or width deviates from the expected value, you should check to see if the default fit is bad.  Look at the blue fit in the upper right plot of the third canvas in the PDF file (panel 3b).  This is the default fit.  Does it describe the data?  Does the residual background (dashed line) appear to be reasonable.  If the default fit looks bad, look through the other fits and see if you can find a fit that looks OK.  You can read out the yield, mass, and width from the fit_results histogram.  Compare these values to the expected values.  It might also make sense to compare the results of all fits to the expected value.  Are the fit parameters always outside the acceptable range?

If the phi yield deviates from the expected value, the following should be noted.  The yield depends on
1.) collision system and energy
2.) the triggers used
3.) multiplicity or centrality range
4.) cuts used to select the decay daughters
5.) pT and rapidity range for phi mesons

The expected yields were calculated for
1.) pp and Pb-Pb collisions at 2.76 TeV, pp collisions at 7 TeV
2.) minimum-bias triggers
3.) all multiplicities for pp, centrality 0-90% for Pb-Pb
4.) the track selection cuts used in the original version of the analysis task
5.) 0.5<pT<1.5 GeV/c and |y|<0.5 for phi mesons

Changes to any of these criteria can cause the yield to change.  If the criteria you are using do not exactly match the criteria for the expected yield, you will need to correct for those differences.  The macro corrects for the rapidity range, so SMALL changes to the rapidity window will be accounted for to some extent.

If the width deviates from the expected value, it is possible that the resolution is not 1.1 MeV/c^2 (the value to which it is currently fixed).  Studying this issue is beyond the scope of these instructions, but you should be aware of the possibility.

*References

[1]: B. Abelev et al. (ALICE Collaboration), Eur. Phys. J. C 72, 2183 (2012)
[2]: Paper in preparation: "K*(892)^0 and phi(1020) resonances in Pb-Pb collisions at 2.76 TeV", by the ALICE Collaboration, intended for publication in Phys. Rev. C (2014)
[3]: Analysis Note: A. G. Knospe, "Yield of phi mesons at low pT in Pb-Pb collisions at 2.76 TeV (2010 data)", ALICE-ANA-2012-300, https://aliceinfo.cern.ch/Notes/node/42
 
  -----------------------------------------------*/


