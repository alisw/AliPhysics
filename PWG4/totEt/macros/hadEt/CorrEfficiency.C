//Christine Nattrass, University of Tennessee at Knoxville
//This macro is for calculating the single track efficiency for the hadronic transverse energy measurement
//Uses the output of AliAnalysisTaskHadEt
//This is not actually what gets used in the correction class AliAnalysisHadEtCorrections - that is done in the macro GetCorrections.C - but this is useful for making plots and playing around with different options

//this function calculates the efficiency propagating errors properly
TH1D* bayneseffdiv(TH1D* numerator, TH1D* denominator,Char_t* name) 
{
    if(!numerator){
      cerr<<"Error:  numerator does not exist!"<<endl;
      return NULL;
    }
    if(!denominator){
      cerr<<"Error:  denominator does not exist!"<<endl;
      return NULL;
    }
    TH1D* result = (TH1D*)numerator->Clone(name);
    Int_t nbins = numerator->GetNbinsX();
    for (Int_t ibin=0; ibin<= nbins+1; ++ibin) {
      Double_t numeratorVal = numerator->GetBinContent(ibin);
      Double_t denominatorVal = denominator->GetBinContent(ibin);
      // Check if the errors are right or the thing is scaled
      Double_t numeratorValErr = numerator->GetBinError(ibin);
      if (!(numeratorValErr==0. || numeratorVal ==0.) ) {
	Double_t rescale = numeratorValErr*numeratorValErr/numeratorVal;
	numeratorVal /= rescale;
      }
      Double_t denominatorValErr = denominator->GetBinError(ibin);
      if (!(denominatorValErr==0. || denominatorVal==0. )) {
	Double_t rescale = denominatorValErr*denominatorValErr/denominatorVal;
	denominatorVal /= rescale;
      }
      Double_t quotient = 0.;
      if (denominatorVal!=0.) {
	quotient = numeratorVal/denominatorVal;
      }
      Double_t quotientErr=0;
      quotientErr = TMath::Sqrt(
				(numeratorVal+1.0)/(denominatorVal+2.0)*
				((numeratorVal+2.0)/(denominatorVal+3.0)-(numeratorVal+1.0)/(denominatorVal+2.0)));
      result->SetBinContent(ibin,quotient);
      result->SetBinError(ibin,quotientErr);
      //cout<<"Setting bin "<<ibin<<" to "<<quotient<<" "<<numeratorVal<<"/"<<denominatorVal<<endl;
    }
    return result;
}


//This is a somewhat messy function that gets the efficiency for different particles
TH1D *GetHisto(float cut = 0.12, char *name, int mycase, bool eta, int color, int marker,bool TPC, bool ITS, bool PbPb, int cb, int cblast){
  //TFile *file = new TFile("Et.ESD.new.sim.merged.root");
  TFile *file = new TFile("rootFiles/LHC11a4_bis/Et.ESD.new.sim.LHC11a4_bis.root");
  //TFile *file = new TFile("rootFiles/LHC10d7/Et.ESD.new.sim.LHC10d7.root");
  TList *list = file->FindObject("out2");
  char *myname = "ITS";
  if(TPC&&!ITS) myname = "TPC";
  if(TPC&&ITS) myname = "TPCITS";
  cout<<"Using tracks from "<<myname<<" for efficiency"<<endl;
  TH2F *numeratorParent; 
  switch(mycase){
  case 0:
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(i==cb) numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiPlus",i)))->Clone("RecoHadron");
	else{numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiPlus",i)));}
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiMinus",i)));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KMinus",i)));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KPlus",i)));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"Proton",i)));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"AntiProton",i)));
      }
    }
    else{
      numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiPlus")))->Clone("RecoHadron");
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiMinus")));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KMinus")));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KPlus")));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"Proton")));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"AntiProton")));
    }
    break;
  case 1://pion
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(i==cb) numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiPlus",i)))->Clone("RecoPion");
	else{numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiPlus",i)));}
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiMinus",i)));
      }
    }
    else{
      numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiPlus")))->Clone("RecoPion");
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiMinus")));
    }
    break;
  case 2://kaon
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(i==cb) numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KPlus",i)))->Clone("RecoKaon");
	else{numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KPlus",i)));}
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KMinus",i)));
      }
    }
    else{
      cout<<"I am kaoning here !"<<endl;
      numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KPlus")))->Clone("RecoKaon");
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KMinus")));
      cout<<"Numerator "<<numeratorParent->GetEntries()<<endl;
    }
    break;
  case 3://proton
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(i==cb) numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"Proton",i)))->Clone("RecoProton");
	else{numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"Proton",i)));}
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"AntiProton",i)));
      }
    }
    else{
      cout<<"I am protoning here !"<<endl;
      numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"Proton")))->Clone("RecoProton");
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"AntiProton")));
      cout<<"Numerator "<<numeratorParent->GetEntries()<<endl;
    }
    break;
  case 4://electron
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s%s",myname,"EPlus",cbname)))->Clone("RecoElectron");
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s%s",myname,"EMinus",cbname)));
    break;
  }
  TH2F *denominatorParent; 
  switch(mycase){
  case 0:
    if(cblast != -1){//add more centrality bins
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNSimulatedChargedHadronCB%i",cb)))->Clone("RecoHadron");
      for(int i=cb+1;i<=cblast;i++){
	denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedChargedHadronCB%i",i)));
      }
    }
    else{
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedChargedHadron"))->Clone("RecoHadron");
    }
    break;
  case 1://pion
    if(cblast != -1){//add more centrality bins
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNSimulatedPiPlusCB%i",cb)))->Clone("RecoPion");
      denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedPiMinusCB%i",cb)));
      for(int i=cb+1;i<=cblast;i++){
	denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedPiPlusCB%i",i)));
	denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedPiMinusCB%i",i)));
      }
    }
    else{
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedPiPlus"))->Clone("RecoPion");
      denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedPiMinus"));
    }
    break;
  case 2://kaon
    if(cblast != -1){//add more centrality bins
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNSimulatedKPlusCB%i",cb)))->Clone("RecoKaon");
      denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedKMinusCB%i",cb)));
      for(int i=cb+1;i<=cblast;i++){
	denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedKPlusCB%i",i)));
	denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedKMinusCB%i",i)));
      }
    }
    else{
      cout<<"I am here kaoning"<<endl;
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedKPlus"))->Clone("RecoKaon");
      denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedKMinus"));
      cout<<"Denominator "<<denominatorParent->GetEntries()<<endl;
    }
    break;
  case 3://proton
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(cb==i)denominatorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNSimulatedProtonCB%i",i)))->Clone("RecoProton");
	else{denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedProtonCB%i",i)));}
	denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedAntiProtonCB%i",i)));
      }
    }
    else{
      cout<<"I am here protoning"<<endl;
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedProton"))->Clone("RecoProton");
      denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedAntiProton"));
      cout<<"Denominator "<<denominatorParent->GetEntries()<<endl;
    }
    break;
  case 4://electron
    denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedEPlus"))->Clone("RecoElectron");
    denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedEMinus"));
    break;
  }
      cout<<"Numerator "<<numeratorParent->GetEntries()<<endl;
      cout<<"Denominator "<<denominatorParent->GetEntries()<<endl;
  numeratorParent->Sumw2();
  denominatorParent->Sumw2();
  TH1D *denominator;
  TH1D *numerator;
  if(eta){
    int lowbin = numeratorParent->GetYaxis()->FindBin(-cut+.001);//make sure we don't accv0entally get the wrong bin
    int highbin = numeratorParent->GetYaxis()->FindBin(cut-.001);
    cout<<"Projecting from "<<numeratorParent->GetYaxis()->GetBinLowEdge(lowbin)<<" to "<<numeratorParent->GetYaxis()->GetBinLowEdge(highbin+1)<<endl;
    denominator = denominatorParent->ProjectionX(Form("garbage%s",name),lowbin,highbin);
    numerator = numeratorParent->ProjectionX(name,lowbin,highbin);
  }
  else{
    int lowbin = denominatorParent->GetXaxis()->FindBin(cut);//make sure we don't accidentally get the wrong bin
    int highbin = denominatorParent->GetXaxis()->GetNbins();
    cout<<"Here Projecting from "<<denominatorParent->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<denominatorParent->GetXaxis()->GetBinLowEdge(highbin+1)<<endl;
    numerator = numeratorParent->ProjectionY(name,lowbin,highbin);
    denominator = denominatorParent->ProjectionY(Form("denominator%s",name),lowbin,highbin);
  }
  delete numeratorParent;
  delete denominatorParent;
  //numerator->Divide(denominator);
  TH1D *result = bayneseffdiv((TH1D*) numerator,(TH1D*)denominator,name);
  //result->Rebin(2);
  //result->Scale(0.5);
  result->SetYTitle("Efficiency");
  result->GetYaxis()->SetTitleOffset(0.8);
  result->GetXaxis()->SetTitleOffset(0.8);
  result->GetYaxis()->SetLabelSize(0.05);
  result->GetXaxis()->SetLabelSize(0.05);
  result->GetYaxis()->SetTitleSize(0.05);
  result->GetXaxis()->SetTitleSize(0.05);
  result->SetMarkerColor(color);
  result->SetLineColor(color);
  result->SetMarkerStyle(marker);
  //result->Draw("e");
  return result;

}


//this is a method that makes pretty plots
void CorrEfficiency(char *shortprodname= "LHC11a4_bis", char *prodname = "LHC11a4_bis HIJING 2.76 TeV Pb+Pb", bool TPC = true,bool ITS = true, bool eta = true, bool PbPb = true, int cb = 0, int cblast = -1){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",600,400);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.02);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  //c->SetLogx();

  int colortotal = 1;
  int colorpi = 2;
  int colork = 3;
  int colorp = 4;
  int phosmarker = 20;
  int emcalmarker = 24;
  float ptcut1 = 0.05;
  float ptcut2 = 0.1;
  float phoscut = 0.12;
  float emcalcut = 0.7;
  if(!eta){
    phoscut = 0.1;
    emcalcut = 0.15;
  }
  TH1D *PHOStotal = GetHisto(phoscut,"PHOStotal",0,eta,colortotal,phosmarker,TPC,ITS,PbPb,cb,cblast);
  TH1D *PHOSpi = GetHisto(phoscut,"PHOSpi",1,eta,colorpi,phosmarker,TPC,ITS,PbPb,cb,cblast);
  TH1D *PHOSp = GetHisto(phoscut,"PHOSp",3,eta,colork,phosmarker,TPC,ITS,PbPb,cb,cblast);
  TH1D *PHOSk = GetHisto(phoscut,"PHOSk",2,eta,colorp,phosmarker,TPC,ITS,PbPb,cb,cblast);
  if(eta) PHOStotal->GetXaxis()->SetRange(PHOStotal->GetXaxis()->FindBin(0.05),PHOStotal->GetXaxis()->FindBin(1.0));
//if(ITS&&!TPC){PHOStotal->GetXaxis()->SetRange(PHOStotal->GetXaxis()->FindBin(0.05),PHOStotal->GetXaxis()->FindBin(1.0));}
//else{PHOStotal->GetXaxis()->SetRange(PHOStotal->GetXaxis()->FindBin(0.0),PHOStotal->GetXaxis()->FindBin(3.0));}
  PHOStotal->SetMinimum(0.0);
  PHOStotal->SetMaximum(1.0);
  //parameters[centbin][0]*exp(-pow(parameters[centbin][1]/pt,parameters[centbin][2]))
  TF1 *func = new TF1("func","[0]*exp(-pow([1]/x,[2]))");
  func->SetParameter(0,.9);
  func->SetParameter(1,.05);
  func->SetParLimits(1,1e-3,1);
  func->SetParameter(2,.1);
  //PHOStotal->Fit(func);
  PHOStotal->Draw();
  PHOSpi->Draw("same");
  PHOSp->Draw("same");
  PHOSk->Draw("same");
  cout<<"Hadrons"<<endl;
  TH1D *EMCALtotal = GetHisto(emcalcut,"EMCALtotal",0,eta,colortotal,emcalmarker,TPC,ITS,PbPb,cb,cblast);
  cout<<endl<<endl<<"=================================PIONS================================="<<endl;
  TH1D *EMCALpi = GetHisto(emcalcut,"EMCALpi",1,eta,colorpi,emcalmarker,TPC,ITS,PbPb,cb,cblast);
  cout<<endl<<endl<<"=================================PROTONS================================="<<endl;
  TH1D *EMCALp = GetHisto(emcalcut,"EMCALp",3,eta,colork,emcalmarker,TPC,ITS,PbPb,cb,cblast);
  cout<<endl<<endl<<"=================================KAONS================================="<<endl;
  TH1D *EMCALk = GetHisto(emcalcut,"EMCALk",2,eta,colorp,emcalmarker,TPC,ITS,PbPb,cb,cblast);
  EMCALtotal->Draw("same");
  EMCALpi->Draw("same");
  EMCALp->Draw("same");
  EMCALk->Draw("same");


  TLegend *leg = new  TLegend(0.22651,0.247312,0.370805,0.438172);
  leg->AddEntry(PHOStotal,"#pi,K,p");
  leg->AddEntry(PHOSpi,"#pi^{#pm}");
  leg->AddEntry(PHOSk,"K^{#pm}");
  leg->AddEntry(PHOSp,"p,#bar{p}");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);
 leg->Draw();

  TLine *line = new TLine(0.150,0.0,0.150,1.0);
  if(eta)line->Draw();
  line->SetLineWidth(3.0);
  //line->SetLineColor(TColor::kYellow);
  line->SetLineStyle(2);
  TLatex *tex = new TLatex(0.398954,0.0513196,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  TLatex *tex3 = new TLatex(1.16186,0.28348,"Closed symbols |#eta|<0.12 (PHOS)");
  tex3->SetTextSize(0.0537634);
  tex3->Draw();
  TLatex *tex4 = new TLatex(1.16186,0.213221,"Open symbols |#eta|<0.70 (EMCal)");
  tex4->SetTextSize(0.0537634);
  tex4->Draw();
  TLatex *tex2 = new TLatex(0.164016,0.860826,"TPC cut-off 150 MeV/c");
  tex2->SetTextSize(0.0537634);
  if(eta) tex2->Draw();


  TLine *line2 = new TLine(0.10,0.0,0.10,1.0);
  line2->SetLineWidth(3.0);
  TLatex *tex5 = new TLatex(0.10817,0.924976,"ITS cut-off 100 MeV/c");
  tex5->SetTextSize(0.0537634);
  line2->SetLineStyle(2);
  tex5->SetTextColor(4);
  line2->SetLineColor(4);
  if(!TPC && eta){
    line2->Draw();
    tex5->Draw();
  }
  if(!TPC){
    c->SaveAs("pics/CorrEfficiency.eps");
    c->SaveAs("pics/CorrEfficiency.png");
  }
  else{
    c->SaveAs("pics/CorrEfficiencyTPCITS.eps");
    c->SaveAs("pics/CorrEfficiencyTPCITS.png");
  }

  if(PbPb){//make one more plot
    //pions - no real centrality dependence
    //three centrality bins for efficiency 0-25%, 25-50%, 50-90%
    //use same for unidentified hadrons
    //kaons & protons - centrality dependence is more significant but I don't think I can do better on the binning
    int pid = 0;//h=0,pi=1,p=3,k=2
  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  c2->SetTopMargin(0.02);
  c2->SetRightMargin(0.02);
  c2->SetBorderSize(0);
  c2->SetFillColor(0);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetFrameFillColor(0);
  c2->SetFrameBorderMode(0);
  cout<<endl<<endl;

  TH1D *cb0 = GetHisto(phoscut,"cb0",pid,eta,1,20,TPC,ITS,PbPb,0,4);
  TH1D *cb4 = GetHisto(phoscut,"cb5",pid,eta,4,20,TPC,ITS,PbPb,5,9);
  //TH1D *cb9 = GetHisto(phoscut,"cb9",pid,eta,TColor::kGreen+4,20,TPC,ITS,PbPb,10,14);
  TH1D *cb14 = GetHisto(phoscut,"cb14",pid,eta,2,20,TPC,ITS,PbPb,10,18);
  cb0->GetXaxis()->SetRange(cb0->GetXaxis()->FindBin(0.05),cb0->GetXaxis()->FindBin(1.0));
  cb0->SetMaximum(1.0);
  cb0->Draw();
  cb4->Draw("same");
  // cb9->Draw("same");
  cb14->Draw("same");
  TLegend *leg = new TLegend(0.124161,0.747312,0.318792,0.959677);//(left,bottom,right,top)
   leg->SetTextSize(0.0537634);
   leg->SetBorderSize(0);
   leg->SetLineColor(0);
   leg->SetLineStyle(0);
   leg->SetLineWidth(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry(cb0,"0-25%");
   leg->AddEntry(cb4,"25-50%");
   //leg->AddEntry(cb9,"45-50%");
   leg->AddEntry(cb14,"50-90%");
   leg->Draw();
   c2->SaveAs(Form("pics/CorrEfficiency%i.png",pid));
  }
}
