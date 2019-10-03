void SetStyles(TH1 *histo,int marker, int color,char *xtitle = "p", char *ytitle = "E");

TH1* bayneseffdiv(TH1* numerator, TH1* denominator,Char_t* name);
void PlotMatchedTracksDepositsFromData(Bool_t isPhos = kFALSE,Bool_t allCells = kTRUE){
  TString datafilename;
  TString simfilename;
  TString det = "Phos";
  //Bool_t isPhos = kTRUE;
  if(!isPhos){
    det = "Emcal";
    isPhos = kFALSE;
    datafilename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.root";
    simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root";
  }
  else{
    simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS.LHC11a10a_bis.Run139465.root";
    datafilename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.PHOS.LHC10hPass2.Run139465.root";
  }
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TString outname = "/tmp/DataVsSimComparisons"+det;
  if(allCells){outname+="AllCells";}
  else{outname+="BkgdOnly";}
  TString textfilename = outname+".txt";
  TString textfilenameEffCorr = outname+"EffCorr.txt";
  TString outnameEffCorr = outname+"EffCorr.png";
  outname+=".png";
  TString outnameratios = "/tmp/DataVsSimComparisonsRatios"+det;
  if(allCells){outnameratios+="AllCells";}
  else{outnameratios+="BkgdOnly";}
  TString outnameratiosEffCorr = outnameratios+"EffCorr.png";
  outnameratios+=".png";
  float nbins = 18;
  TObjArray hMatchedTracksSimulation(nbins+1);
  TObjArray data(nbins+1);
  TObjArray ratios(nbins+1);
  TObjArray hMatchedTracksSimulationEffCorr(nbins+1);
  TObjArray dataEffCorr(nbins+1);
  TObjArray dataEffCorr500MeV(nbins+1);
  TObjArray ratiosEffCorr(nbins+1);
  TObjArray trackHit(nbins+1);
  TObjArray trackDep(nbins+1);
  TObjArray trackDep500MeV(nbins+1);
  TObjArray trackrecoDep(nbins+1);
  TObjArray trackEffVsPt(nbins+1);
  TObjArray trackrecoEffVsPt(nbins+1);
  TObjArray etForLowPtTracks(nbins+1);
  TObjArray etForLowPtTracksEffCorr(nbins+1);
  TObjArray etForLowPtTracks500MeV(nbins+1); 
  TObjArray etForLowPtTracksNoAntiProtons(nbins+1); 
  TObjArray etForLowPtTracksAntiProtons(nbins+1); 
  TObjArray nFound(nbins+1); 
  TObjArray nNotFound(nbins+1); 
  TObjArray nFound500MeV(nbins+1); 
  TObjArray nNotFound500MeV(nbins+1); 
//   TH1D *hMatchedTracksSimulation[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};
//   TH1D *data[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};
  int colors[] = {TColor::kRed, TColor::kOrange-3,TColor::kGreen+3, TColor::kBlue, TColor::kViolet,
		  TColor::kRed, TColor::kOrange-3,TColor::kGreen+3, TColor::kBlue, TColor::kViolet,
		  TColor::kRed, TColor::kOrange-3,TColor::kGreen+3, TColor::kBlue, TColor::kViolet,
		  TColor::kRed, TColor::kOrange-3,TColor::kGreen+3, TColor::kBlue, TColor::kViolet};
  int markersSim[] = {24,25,26,32,27,24,25,26,32,27,24,25,26,32,27,24,25,26,32,27};
  int markersData[] = {20,21,22,23,33,20,21,22,23,33,20,21,22,23,33,20,21,22,23,33};
  int colorsAlt[] = {TColor::kRed,TColor::kRed, TColor::kOrange-3,TColor::kOrange-3,TColor::kGreen+3,TColor::kGreen+3, TColor::kBlue, TColor::kBlue, TColor::kViolet, TColor::kViolet,TColor::kRed,TColor::kRed, TColor::kOrange-3,TColor::kOrange-3,TColor::kGreen+3,TColor::kGreen+3, TColor::kBlue, TColor::kBlue, TColor::kViolet, TColor::kViolet};
  int markersAlt[] = {20,24,21,25,22,26,23,32,33,27,20,24,21,25,22,26,23,32,33,27};
  float belowCut[] = {0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0,   0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0};
  float total[] = {0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0,   0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0};
  float belowCut500MeV[] = {0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0,   0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0};
  float total500MeV[] = {0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0,   0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0};
  float found[] = {0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0,   0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0};
  float notfound[] = {0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0,   0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0};
  float found500MeV[] = {0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0,   0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0};
  float notfound500MeV[] = {0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0,   0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,  0.0};

  Float_t npart[20] = {382.8,329.7,281.1,238.6,201.8,
                       169.5,141.6,116.2,94.09,75.45,
                       59.40,45.75,34.48,25.07,17.76,
                       12.49,8.810,6.150,4.351,3.060};



  float energycut = 0.0;

  TFile *fsim = TFile::Open(simfilename, "READ");
  TList *lsim = (TList*)fsim->Get("out1");

  TH3F  *fHistMatchedTracksEvspTBkgdvsCent = lsim->FindObject("fHistMatchedTracksEvspTBkgdvsCent");
  TH3F  *fHistMatchedTracksEvspTSignalvsCent = lsim->FindObject("fHistMatchedTracksEvspTSignalvsCent");
  TH3F  *fHistMatchedTracksEvspTBkgdvsCentEffCorr = lsim->FindObject("fHistMatchedTracksEvspTBkgdvsCentEffCorr");
  TH3F  *fHistMatchedTracksEvspTSignalvsCentEffCorr = lsim->FindObject("fHistMatchedTracksEvspTSignalvsCentEffCorr");
  TH2F *fHistHadronsAllCent = lsim->FindObject("fHistHadronsAllCent");//ALL Charged tracks
  TH2F *fHistHadronDepositsAllCent = lsim->FindObject("fHistHadronDepositsAllCent");//All charged tracks which leave deposits in the calorimeter
  TH2F *fHistHadronDepositsAllCent500MeV = lsim->FindObject("fHistHadronDepositsAllCent500MeV");//All charged tracks which leave deposits in the calorimeter
  TH2F *fHistHadronDepositsRecoCent = lsim->FindObject("fHistHadronDepositsRecoCent");//All charged tracks which leave deposits in the calorimeter
  TH2F *fHistChargedTracksAcceptedLowPtCent = lsim->FindObject("fHistChargedTracksAcceptedLowPtCent");//distribution of energies left by low pT tracks in simulation
  TH2F *fHistChargedTracksAcceptedLowPtCentEffCorr = lsim->FindObject("fHistChargedTracksAcceptedLowPtCentEffCorr");//distribution of energies left by low pT tracks in simulation
  TH2F *fHistChargedTracksAcceptedLowPtCent500MeV = lsim->FindObject("fHistChargedTracksAcceptedLowPtCent500MeV");//distribution of energies left by low pT tracks in simulation
  TH2F *fHistChargedTracksAcceptedLowPtCentNoAntiProtons = lsim->FindObject("fHistChargedTracksAcceptedLowPtCentNoAntiProtons");//distribution of energies left by low pT tracks in simulation
  TH2F *fHistChargedTracksAcceptedLowPtCentAntiProtons = lsim->FindObject("fHistChargedTracksAcceptedLowPtCentAntiProtons");//distribution of energies left by low pT tracks in simulation
  if(allCells){
    fHistMatchedTracksEvspTBkgdvsCent->Add(fHistMatchedTracksEvspTSignalvsCent);
    fHistMatchedTracksEvspTBkgdvsCentEffCorr->Add(fHistMatchedTracksEvspTSignalvsCentEffCorr);
  }
  int lowbin = fHistMatchedTracksEvspTBkgdvsCent->GetXaxis()->FindBin(energycut);
  int highbin = fHistMatchedTracksEvspTBkgdvsCent->GetXaxis()->GetNbins();
  fHistMatchedTracksEvspTBkgdvsCent->GetXaxis()->SetRange(lowbin,highbin);
  fHistMatchedTracksEvspTBkgdvsCentEffCorr->GetXaxis()->SetRange(lowbin,highbin);

  for(int bin = 1; bin<nbins;bin++){
    fHistMatchedTracksEvspTBkgdvsCent->GetZaxis()->SetRange(bin,bin);
    hMatchedTracksSimulation[bin] = fHistMatchedTracksEvspTBkgdvsCent->ProjectionY("y");
    ((TH1D*)hMatchedTracksSimulation[bin])->SetName(Form("Sim%i",bin));
    SetStyles((TH1D*)hMatchedTracksSimulation[bin],markersSim[bin-1],colors[bin-1]);
    fHistMatchedTracksEvspTBkgdvsCentEffCorr->GetZaxis()->SetRange(bin,bin);
    hMatchedTracksSimulationEffCorr[bin] = fHistMatchedTracksEvspTBkgdvsCentEffCorr->ProjectionY("y");
    ((TH1D*)hMatchedTracksSimulationEffCorr[bin])->SetName(Form("Sim%i",bin));
    SetStyles((TH1D*)hMatchedTracksSimulationEffCorr[bin],markersSim[bin-1],colors[bin-1]);
    //
    //   TObjArray trackHit(nbins+1);
    //   TObjArray trackDep(nbins+1);
    //   TObjArray trackEffVsPt(nbins+1);
    //fHistHadronDepositsRecoCent
    trackHit[bin] = fHistHadronsAllCent->ProjectionX(Form("trackAllCentVsPt%i",bin),bin,bin);
    trackDep[bin] = fHistHadronDepositsAllCent->ProjectionX(Form("trackDepVsPt%i",bin),bin,bin);
    trackDep500MeV[bin] = fHistHadronDepositsAllCent500MeV->ProjectionX(Form("trackDepVsPt500MeV%i",bin),bin,bin);
    trackrecoDep[bin] = fHistHadronDepositsRecoCent->ProjectionX(Form("trackrecoDepVsPt%i",bin),bin,bin);
    trackEffVsPt[bin] = bayneseffdiv((TH1*)trackDep[bin],(TH1*)trackHit[bin],Form("trackEffVsPt%i",bin));
    trackrecoEffVsPt[bin] = bayneseffdiv((TH1*)trackrecoDep[bin],(TH1*)trackDep[bin],Form("trackrecoEffVsPt%i",bin));
    SetStyles((TH1D*)trackEffVsPt[bin],markersAlt[bin-1],colorsAlt[bin-1]);
    SetStyles((TH1D*)trackrecoEffVsPt[bin],markersAlt[bin-1],colorsAlt[bin-1]);
    belowCut[bin] = ((TH1D*)trackDep[bin])->Integral(1,((TH1D*)trackDep[bin])->FindBin(0.49));
    total[bin] = ((TH1D*)trackDep[bin])->Integral();
    belowCut500MeV[bin] = ((TH1D*)trackDep500MeV[bin])->Integral(1,((TH1D*)trackDep500MeV[bin])->FindBin(0.49));
    total500MeV[bin] = ((TH1D*)trackDep500MeV[bin])->Integral();
    etForLowPtTracks[bin] = fHistChargedTracksAcceptedLowPtCent->ProjectionX(Form("etForLowPtTracks%i",bin),bin,bin);
    etForLowPtTracks500MeV[bin] = fHistChargedTracksAcceptedLowPtCent500MeV->ProjectionX(Form("etForLowPtTracks500MeV%i",bin),bin,bin);
    //etForLowPtTracks500MeV[bin] = fHistChargedTracksAcceptedLowPtCent500MeV->ProjectionX(Form("etForLowPtTracks500MeV%i",bin),bin,bin);
    etForLowPtTracksNoAntiProtons[bin] = fHistChargedTracksAcceptedLowPtCentNoAntiProtons->ProjectionX(Form("etForLowPtTracksNoAntiProtons%i",bin),bin,bin);
    etForLowPtTracksAntiProtons[bin] = fHistChargedTracksAcceptedLowPtCentAntiProtons->ProjectionX(Form("etForLowPtTracksAntiProtons%i",bin),bin,bin);
    cout<<"fraction bin "<<bin;
    //cout<<" "<<belowCut[bin]/total[bin]<<" below "<<belowCut[bin]<<" total "<<total[bin];
    //cout<<" fraction > 500 MeV bin "<<bin<<" "<<belowCut500MeV[bin]/total500MeV[bin];
    //cout<<" below cut "<<belowCut500MeV[bin]<<" total "<< total500MeV[bin];
    //cout<<" <et> low pT 500 MeV "<<((TH1D*)etForLowPtTracks500MeV[bin])->GetMean();
    //cout<<" <et> low pT "<<((TH1D*)etForLowPtTracks[bin])->GetMean();
    //cout<<" <et> high pT "<<((TH1D*)hMatchedTracksSimulation[bin])->GetMean();
    float low = ((TH1D*)etForLowPtTracks[bin])->GetMean();
    float high = ((TH1D*)hMatchedTracksSimulationEffCorr[bin])->GetMean();
    cout<<" low/total "<<low/(low+high)<<" low/high "<<low/high;
    //cout<<" low/high "<<((TH1D*)etForLowPtTracks[bin])->GetMean()/((TH1D*)hMatchedTracksSimulation[bin])->GetMean();
    //cout<<" <et> no antiprotons "<<((TH1D*)etForLowPtTracksNoAntiProtons[bin])->GetMean();
    //cout<<" mean low pT et antiprotons "<<((TH1D*)etForLowPtTracksAntiProtons[bin])->GetMean();
    cout<<endl;
    SetStyles((TH1D*)etForLowPtTracks[bin],markersAlt[bin-1],colorsAlt[bin-1]);
    SetStyles((TH1D*)etForLowPtTracksNoAntiProtons[bin],markersAlt[bin-1],colorsAlt[bin-1]);
    SetStyles((TH1D*)etForLowPtTracksAntiProtons[bin],markersAlt[bin-1],colorsAlt[bin-1]);
  }
  TFile *f = TFile::Open(datafilename, "READ");
  f->cd();  
  TList *l = (TList*)f->Get("out1");
  TH3F  *fHistMatchedTracksEvspTvsCent = l->FindObject("fHistMatchedTracksEvspTvsCent");
  int lowbin = fHistMatchedTracksEvspTvsCent->GetXaxis()->FindBin(energycut);
  int highbin = fHistMatchedTracksEvspTvsCent->GetXaxis()->GetNbins();
  fHistMatchedTracksEvspTvsCent->GetXaxis()->SetRange(lowbin,highbin);
  // TH3F  *fHistMatchedTracksEvspTvsCentEffCorr = l->FindObject("fHistMatchedTracksEvspTvsCentEffCorr");
  TH3F  *fHistMatchedTracksEvspTvsCentEffCorr = l->FindObject("fHistMatchedTracksEvspTvsCentEffTMCorr");
  TH3F  *fHistMatchedTracksEvspTvsCentEffCorr500MeV = l->FindObject("fHistMatchedTracksEvspTvsCentEffTMCorr500MeV");
  int lowbin = fHistMatchedTracksEvspTvsCentEffCorr->GetXaxis()->FindBin(energycut);
  int highbin = fHistMatchedTracksEvspTvsCentEffCorr->GetXaxis()->GetNbins();
  fHistMatchedTracksEvspTvsCentEffCorr->GetXaxis()->SetRange(lowbin,highbin);
  TH2F *fHistFoundHadronsvsCent = l->FindObject("fHistFoundHadronsvsCent");
  TH2F *fHistNotFoundHadronsvsCent = l->FindObject("fHistNotFoundHadronsvsCent");
  TH2F *fHistFoundHadronsvsCent500MeV = l->FindObject("fHistFoundHadronsvsCent500MeV");
  TH2F *fHistNotFoundHadronsvsCent500MeV = l->FindObject("fHistNotFoundHadronsvsCent500MeV");

  for(int bin = 1; bin<nbins;bin++){
    fHistMatchedTracksEvspTvsCent->GetZaxis()->SetRange(bin,bin);
    data[bin] = fHistMatchedTracksEvspTvsCent->Project3D("y");
    SetStyles((TH1D*)data[bin],markersData[bin-1],colors[bin-1]);
    ((TH1D*)data[bin])->SetName(Form("Sim%i",bin));
    fHistMatchedTracksEvspTvsCentEffCorr->GetZaxis()->SetRange(bin,bin);
    dataEffCorr[bin] = fHistMatchedTracksEvspTvsCentEffCorr->Project3D("y");
    fHistMatchedTracksEvspTvsCentEffCorr500MeV->GetZaxis()->SetRange(bin,bin);
    dataEffCorr500MeV[bin] = fHistMatchedTracksEvspTvsCentEffCorr500MeV->Project3D("y");
    SetStyles((TH1D*)dataEffCorr[bin],markersData[bin-1],colors[bin-1]);
    ((TH1D*)dataEffCorr[bin])->SetName(Form("DataEffCorr%i",bin));
    ((TH1D*)dataEffCorr500MeV[bin])->SetName(Form("DataEffCorr500MeV%i",bin));
    nFound[bin] = fHistFoundHadronsvsCent->ProjectionX(Form("Found%i",bin),bin,bin);
    nNotFound[bin] = fHistNotFoundHadronsvsCent->ProjectionX(Form("NotFound%i",bin),bin,bin);
    found[bin] = ((TH1D*)nFound[bin])->GetMean();//fHistFoundHadronsvsCent->GetBinContent(bin);
    notfound[bin] = ((TH1D*)nNotFound[bin])->GetMean();//fHistNotFoundHadronsvsCent->GetBinContent(bin);
    nFound500MeV[bin] = fHistFoundHadronsvsCent500MeV->ProjectionX(Form("Found%i",bin),bin,bin);
    nNotFound500MeV[bin] = fHistNotFoundHadronsvsCent500MeV->ProjectionX(Form("NotFound%i",bin),bin,bin);
    found500MeV[bin] = ((TH1D*)nFound500MeV[bin])->GetMean();//fHistFoundHadronsvsCent->GetBinContent(bin);
    notfound500MeV[bin] = ((TH1D*)nNotFound500MeV[bin])->GetMean();//fHistNotFoundHadronsvsCent->GetBinContent(bin);
    cout<<"bin "<<bin<<" low pT found "<<found[bin]<<" not found "<<notfound[bin]<<" high pt found "<<found500MeV[bin]<<" not found "<<notfound500MeV[bin]<<endl;
  }
  for(int bin = 1; bin<nbins;bin++){
    float scale = ((TH1D*)data[bin])->Integral();
    ((TH1D*)data[bin])->Scale(1.0/scale);
    scale = ((TH1D*)hMatchedTracksSimulation[bin])->Integral();
    ((TH1D*)hMatchedTracksSimulation[bin])->Scale(1.0/scale);
    ratios[bin] = data[bin]->Clone(Form("ratio%i",bin));
    ((TH1D*)ratios[bin])->Divide((TH1D*)hMatchedTracksSimulation[bin]);
    scale = ((TH1D*)dataEffCorr[bin])->Integral();
    ((TH1D*)dataEffCorr[bin])->Scale(1.0/scale);
    scale = ((TH1D*)dataEffCorr500MeV[bin])->Integral();
    ((TH1D*)dataEffCorr500MeV[bin])->Scale(1.0/scale);
    scale = ((TH1D*)hMatchedTracksSimulationEffCorr[bin])->Integral();
    ((TH1D*)hMatchedTracksSimulationEffCorr[bin])->Scale(1.0/scale);
    ratiosEffCorr[bin] = dataEffCorr[bin]->Clone(Form("ratio%i",bin));
    ((TH1D*)ratiosEffCorr[bin])->Divide((TH1D*)hMatchedTracksSimulationEffCorr[bin]);
  }

  TLegend *leg = new TLegend(0.631661,0.550143,0.731975,0.955587);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextSize(0.038682);
  TCanvas *c1 = new TCanvas("c1","Data",1200,600);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetLogy();
  ((TH1D *)data[1])->GetXaxis()->SetRange(((TH1D *)data[1])->FindBin(0.3),((TH1D *)data[1])->GetNbinsX());
  ((TH1D *)data[1])->Draw();
  ofstream myfile;
  myfile.open (textfilename.Data());
  myfile<<"data\tsimulation"<<endl;
  for(int bin = 1; bin<nbins;bin++){
    float scale = TMath::Power(2,bin);
    ((TH1D *)data[bin])->Scale(1.0/scale);
    ((TH1D *)hMatchedTracksSimulation[bin])->Scale(1.0/scale);
    leg->AddEntry(data[bin],Form("%i<N_{cluster}<%i times %2.1f",(bin-1)*10,bin*10,scale));
    data[bin]->Draw("same");
    hMatchedTracksSimulation[bin]->Draw("same");
    myfile<< ((TH1D *)data[bin])->GetMean() <<"\t"<< ((TH1D *)hMatchedTracksSimulation[bin])->GetMean()<<endl;

  }
  myfile.close();
  leg->Draw();
  cout<<"Saving "<<outname.Data()<<endl;
  c1->SaveAs(outname.Data());

  TCanvas *c2 = new TCanvas("c2","Data/simulation",1200,600);
  c2->SetTopMargin(0.02);
  c2->SetRightMargin(0.02);
  c2->SetBorderSize(0);
  c2->SetFillColor(0);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetFrameFillColor(0);
  c2->SetFrameBorderMode(0);
  ((TH1D *)ratios[1])->GetYaxis()->SetTitle("data/simulation");

  ((TH1D *)ratios[1])->GetXaxis()->SetRange(((TH1D *)ratios[1])->FindBin(0.3),((TH1D *)ratios[1])->GetNbinsX());
  ((TH1D *)ratios[1])->Draw();
  for(int bin = 1; bin<nbins;bin++){
    //float scale = 2-0.2*bin;
     //((TH1D *)ratios[bin])->Scale(1.0/scale);
    ratios[bin]->Draw("same");
  }
  ((TH1D *)ratios[1])->SetMaximum(1.5);
  ((TH1D *)ratios[1])->SetMinimum(0.0);
  //leg->Draw();
  c2->SaveAs(outnameratios.Data());
  TCanvas *c3 = new TCanvas("c3","Data eff corr",1200,600);
  c3->SetTopMargin(0.02);
  c3->SetRightMargin(0.02);
  c3->SetBorderSize(0);
  c3->SetFillColor(0);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetFrameFillColor(0);
  c3->SetFrameBorderMode(0);
  c3->SetLogy();
  ((TH1D *)dataEffCorr[1])->GetXaxis()->SetRange(((TH1D *)dataEffCorr[1])->FindBin(0.3),((TH1D *)dataEffCorr[1])->GetNbinsX());
  ((TH1D *)dataEffCorr[1])->Draw();
  ofstream myfile2;
  myfile2.open (textfilenameEffCorr.Data());
  myfile2<<"data\tsimulation"<<endl;
  for(int bin = 1; bin<nbins;bin++){
    float scale = TMath::Power(2,bin);
    ((TH1D *)dataEffCorr[bin])->Scale(1.0/scale);
    ((TH1D *)hMatchedTracksSimulationEffCorr[bin])->Scale(1.0/scale);
    //leg->AddEntry(dataEffCorr[bin],Form("%i<N_{cluster}<%i times %2.1f",(bin-1)*10,bin*10,scale));
    dataEffCorr[bin]->Draw("same");
    hMatchedTracksSimulationEffCorr[bin]->Draw("same");
    myfile2<< ((TH1D *)dataEffCorr[bin])->GetMean() <<"\t"<< ((TH1D *)hMatchedTracksSimulationEffCorr[bin])->GetMean()<<endl;
  }
  myfile2.close();
  leg->Draw();
  c3->SaveAs(outnameEffCorr.Data());

  TCanvas *c4 = new TCanvas("c4","data/simulation eff corr",1200,600);
  c4->SetTopMargin(0.02);
  c4->SetRightMargin(0.02);
  c4->SetBorderSize(0);
  c4->SetFillColor(0);
  c4->SetFillColor(0);
  c4->SetBorderMode(0);
  c4->SetFrameFillColor(0);
  c4->SetFrameBorderMode(0);
  ((TH1D *)ratiosEffCorr[1])->GetYaxis()->SetTitle("data/simulation");

  ((TH1D *)ratiosEffCorr[1])->GetXaxis()->SetRange(((TH1D *)ratiosEffCorr[1])->FindBin(0.3),((TH1D *)ratiosEffCorr[1])->GetNbinsX());
  ((TH1D *)ratiosEffCorr[1])->Draw();
  for(int bin = 1; bin<nbins;bin++){
    //float scale = 2-0.2*bin;
     //((TH1D *)ratios[bin])->Scale(1.0/scale);
    ratiosEffCorr[bin]->Draw("same");
  }
  ((TH1D *)ratiosEffCorr[1])->SetMaximum(1.5);
  ((TH1D *)ratiosEffCorr[1])->SetMinimum(0.0);
  //leg->Draw();
  c4->SaveAs(outnameratiosEffCorr.Data());

  TLegend *leg2 = new TLegend(0.151007,0.545699,0.251678,0.951613);
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->SetTextSize(0.038682);
  TCanvas *c5 = new TCanvas("c5","fraction depositing energy",600,400);
  c5->SetTopMargin(0.02);
  c5->SetRightMargin(0.02);
  c5->SetBorderSize(0);
  c5->SetFillColor(0);
  c5->SetFillColor(0);
  c5->SetBorderMode(0);
  c5->SetFrameFillColor(0);
  c5->SetFrameBorderMode(0);
  ((TH1D *)trackEffVsPt[1])->GetYaxis()->SetTitle("fraction depositing energy");

  ((TH1D *)trackEffVsPt[1])->GetXaxis()->SetRange(1,((TH1D *)trackEffVsPt[1])->FindBin(1.5));
  ((TH1D *)trackEffVsPt[1])->GetXaxis()->SetTitle("p_{T}");
  ((TH1D *)trackEffVsPt[1])->Draw();
  for(int bin = 1; bin<nbins;bin++){
    leg2->AddEntry(trackEffVsPt[bin],Form("%i<N_{cluster}<%i",(bin-1)*10,bin*10));
    //float scale = 2-0.2*bin;
     //((TH1D *)ratios[bin])->Scale(1.0/scale);
    trackEffVsPt[bin]->Draw("same");
  }
  leg2->Draw();
  ((TH1D *)trackEffVsPt[1])->SetMaximum(0.03);
  if(det.Contains("PHOS")) ((TH1D *)trackEffVsPt[1])->SetMaximum(0.005);
  ((TH1D *)trackEffVsPt[1])->SetMinimum(0.0);
  TString filename = "/tmp/FractionHadronsDeposted"+det+".eps";
  c5->SaveAs(filename.Data());

  TCanvas *c6 = new TCanvas("c6","fraction subtracted",600,400);
  c6->SetTopMargin(0.02);
  c6->SetRightMargin(0.02);
  c6->SetBorderSize(0);
  c6->SetFillColor(0);
  c6->SetFillColor(0);
  c6->SetBorderMode(0);
  c6->SetFrameFillColor(0);
  c6->SetFrameBorderMode(0);
  ((TH1D *)trackrecoEffVsPt[1])->GetYaxis()->SetTitle("fraction subtracted");

  ((TH1D *)trackrecoEffVsPt[1])->GetXaxis()->SetRange(1,((TH1D *)trackrecoEffVsPt[1])->FindBin(1.5));
  ((TH1D *)trackrecoEffVsPt[1])->Draw();
  for(int bin = 1; bin<nbins;bin++){
    //float scale = 2-0.2*bin;
     //((TH1D *)ratios[bin])->Scale(1.0/scale);
    trackrecoEffVsPt[bin]->Draw("same");
  }
  ((TH1D *)trackrecoEffVsPt[1])->SetMaximum(1.0);
  ((TH1D *)trackrecoEffVsPt[1])->SetMinimum(0.0);

  TCanvas *c7 = new TCanvas("c7","fraction subtracted",600,400);
  c7->SetTopMargin(0.02);
  c7->SetRightMargin(0.02);
  c7->SetBorderSize(0);
  c7->SetFillColor(0);
  c7->SetFillColor(0);
  c7->SetBorderMode(0);
  c7->SetFrameFillColor(0);
  c7->SetFrameBorderMode(0);
  ((TH1D *)etForLowPtTracks[1])->GetYaxis()->SetTitle("fraction subtracted");

  ((TH1D *)etForLowPtTracks[1])->GetXaxis()->SetRange(1,((TH1D *)etForLowPtTracks[1])->FindBin(1.5));
  ((TH1D *)etForLowPtTracks[1])->Draw();
  for(int bin = 1; bin<nbins;bin++){
    //float scale = 2-0.2*bin;
     //((TH1D *)ratios[bin])->Scale(1.0/scale);
    etForLowPtTracks[bin]->Draw("same");
  }
  //((TH1D *)etForLowPtTracks[1])->SetMaximum(1.0);
  //((TH1D *)etForLowPtTracks[1])->SetMinimum(0.0);

  ofstream myfile3;
  TString texfilename = "/tmp/energydepositstable"+det;
  if(allCells){texfilename+="AllCells";}
  else{texfilename+="BkgdOnly";}
  texfilename+=".tex";
  myfile3.open (texfilename.Data());
  float low = ((TH1D *)dataEffCorr[1])->GetMean();
  float low500MeV = ((TH1D *)dataEffCorr500MeV[1])->GetMean();
  for(int bin = 2; bin<nbins;bin++){
    if(((TH1D *)dataEffCorr[bin])->GetMean()<low){
      low = ((TH1D *)dataEffCorr[bin])->GetMean();
      low500MeV = ((TH1D *)dataEffCorr500MeV[bin])->GetMean();
    }
  }
  float  factor = 1-0.04;
  float corrfac = 1.275-1;
  float corrfacerr = 0.059 ;
  if(det.Contains("PHOS")){
    factor = 1-0.03;
    corrfac = 1.300-1;
    corrfacerr = 0.065;
  }
  low = factor*low;
  float lowEffCorr = ((TH1D *)hMatchedTracksSimulationEffCorr[1])->GetMean();
  myfile3<<"Simulation & "<<Form("%2.3f",((TH1D *)hMatchedTracksSimulation[1])->GetMean())<<" & "<<Form("%2.3f",((TH1D *)hMatchedTracksSimulationEffCorr[1])->GetMean())<<"\\\\ \\hline"<<endl;//simulation[1]<<" & "<<simulationEffCorr[1]<<"\\\\"<<endl;
  Float_t hadError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
  Float_t hadCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
  Float_t had500MeVError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
  Float_t had500MeVCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
  for(int bin = 1; bin<nbins;bin++){
    float high = ((TH1D *)dataEffCorr[bin])->GetMean();
    float high500MeV = ((TH1D *)dataEffCorr500MeV[bin])->GetMean();
    float avg = (high+low)/2.0;//average energy deposited by the mean hadron (E)
    float err = TMath::Abs((high-low)/2.0);
    float avg500MeV = (high500MeV+low500MeV)/2.0;//average energy deposited by the mean hadron (E)
    float err500MeV = TMath::Abs((high500MeV-low500MeV)/2.0);
    //hadronic correction:  y = f N1 E + N2 E
    //N1 = total number of hadrons
    //N2 = number of hadrons not found
    float percentEfficiencyError = 0.01;
    //fraction in low pT
    float eLow = corrfac * (found[bin]+notfound[bin])*avg;
    //float eLowErr = TMath::Sqrt(TMath::Power(corrfacerr*(found[bin]+notfound[bin])*avg,2)+TMath::Power(err*corrfac* (found[bin]+notfound[bin]),2)+TMath::Power(percentEfficiencyError*eLowErr,2));//error on the hadronic correction
    float eLowErr = TMath::Sqrt(TMath::Power(corrfacerr*(found[bin]+notfound[bin])*avg,2)+TMath::Power(err*corrfac* (found[bin]+notfound[bin]),2)+TMath::Power(eLow*percentEfficiencyError,2));//error on the hadronic correction
    //fraction in hadrons not found
    float eNotFound = notfound[bin]*avg;
    float eNotFoundErr = TMath::Sqrt(TMath::Power(err*notfound[bin],2)+TMath::Power(percentEfficiencyError*eNotFound,2));//error on the hadronic correction
    //cout<<"y = ("<<corrfac <<"*"<<(found[bin]+notfound[bin])<<"+"<<notfound[bin]<<")*"<<avg<<endl;
    float y = (corrfac * (found[bin]+notfound[bin]) +notfound[bin])*avg;
    float y500MeV =  notfound[bin]*avg500MeV;
    //float y = found[bin]+notfound[bin];
    //ey^2 = ef^2 N1^2 E^2 + eE^2 *(  f^2 N1^2 + N2^2)
    float finalerr = TMath::Sqrt(TMath::Power(corrfacerr*(found[bin]+notfound[bin])*avg,2)+err*err*(TMath::Power(corrfac* (found[bin]+notfound[bin]),2)+notfound[bin]*notfound[bin])+TMath::Power(percentEfficiencyError*y,2));//error on the hadronic correction
    float finalerr500MeV = TMath::Sqrt(TMath::Power(corrfacerr*(found[bin]+notfound[bin])*avg,2)+err500MeV*err500MeV*notfound[bin]*notfound[bin]+TMath::Power(percentEfficiencyError*y,2));//error on the hadronic correction
    //if(avg>0 && corrfac>0) finalerr = avg*corrfac*TMath::Sqrt(TMath::Power(corrfacerr/corrfac,2)+TMath::Power(err/avg,2));
    myfile3<<Form("%i-%i",(bin-1)*5,bin*5)<<"\\% & "<<Form("%2.3f",((TH1D *)data[bin])->GetMean())<<" & "<<Form("%2.3f",((TH1D *)dataEffCorr[bin])->GetMean())<<Form("& %2.3f $\\pm$ %2.3f",avg,err);
    myfile3<<" & "<< Form("%2.1f $\\pm$ %2.1f",eLow,eLowErr);
    myfile3<<" & "<< Form("%2.2f $\\pm$ %2.2f",eNotFound,eNotFoundErr);
    myfile3<<"& "<< Form("%2.2f $\\pm$ %2.2f",y,finalerr);
    myfile3<<" & "<< Form("%2.3f $\\pm$ %2.3f",y/npart[bin-1],finalerr/npart[bin-1]) <<"\\\\"<<endl;
    hadError[bin-1] = finalerr;
    hadCorr[bin-1] = y;
    had500MeVError[bin-1] = finalerr500MeV;
    had500MeVCorr[bin-1] = y500MeV;
    //<<Form("& %2.3f $\\pm$ %2.3f",avg*corrfac,finalerr)
  }
  myfile3.close();

  cout<<"Float_t hadError[20] = {";
  for(int i=0;i<19;i++){
    cout<<hadError[i];
    if(i!=19) cout<<",";
  }
  cout<<"};"<<endl;
  cout<<"Float_t hadCorr[20] = {";
  for(int i=0;i<19;i++){
    cout<<hadCorr[i];
    if(i!=19) cout<<",";
  }
  cout<<"};"<<endl;
  cout<<"Float_t hadError500MeV[20] = {";
  for(int i=0;i<19;i++){
    cout<<had500MeVError[i];
    if(i!=19) cout<<",";
  }
  cout<<"};"<<endl;
  cout<<"Float_t hadCorr500MeV[20] = {";
  for(int i=0;i<19;i++){
    cout<<had500MeVCorr[i];
    if(i!=19) cout<<",";
  }
  cout<<"};"<<endl;


  TCanvas *c8 = new TCanvas("c8","c8",600,400);
  c8->SetTopMargin(0.02);
  c8->SetRightMargin(0.02);
  c8->SetBorderSize(0);
  c8->SetFillColor(0);
  c8->SetFillColor(0);
  c8->SetBorderMode(0);
  c8->SetFrameFillColor(0);
  c8->SetFrameBorderMode(0);
  ((TH1D *)etForLowPtTracksNoAntiProtons[1])->GetYaxis()->SetTitle("fraction depositing energy");

  ((TH1D *)etForLowPtTracksNoAntiProtons[1])->GetXaxis()->SetRange(1,((TH1D *)etForLowPtTracksNoAntiProtons[1])->FindBin(1.5));
  ((TH1D *)etForLowPtTracksNoAntiProtons[1])->GetXaxis()->SetTitle("p_{T}");
  ((TH1D *)etForLowPtTracksNoAntiProtons[1])->Draw();
  for(int bin = 1; bin<nbins;bin++){
    //float scale = 2-0.2*bin;
     //((TH1D *)ratios[bin])->Scale(1.0/scale);
    etForLowPtTracksNoAntiProtons[bin]->Draw("same");
  }
  //((TH1D *)etForLowPtTracksNoAntiProtons[1])->SetMaximum(0.03);
  if(det.Contains("PHOS")) ((TH1D *)etForLowPtTracksNoAntiProtons[1])->SetMaximum(400);
  ((TH1D *)etForLowPtTracksNoAntiProtons[1])->SetMinimum(0.0);
  TString filename = "/tmp/FractionHadronsDeposted"+det+".eps";
  c8->SaveAs(filename.Data());

  TCanvas *c9 = new TCanvas("c9","c9",600,400);
  c9->SetTopMargin(0.02);
  c9->SetRightMargin(0.02);
  c9->SetBorderSize(0);
  c9->SetFillColor(0);
  c9->SetFillColor(0);
  c9->SetBorderMode(0);
  c9->SetFrameFillColor(0);
  c9->SetFrameBorderMode(0);
  ((TH1D *)etForLowPtTracksAntiProtons[1])->GetYaxis()->SetTitle("fraction depositing energy");

  ((TH1D *)etForLowPtTracksAntiProtons[1])->GetXaxis()->SetRange(1,((TH1D *)etForLowPtTracksAntiProtons[1])->FindBin(1.5));
  ((TH1D *)etForLowPtTracksAntiProtons[1])->GetXaxis()->SetTitle("p_{T}");
  ((TH1D *)etForLowPtTracksAntiProtons[1])->Draw();
  for(int bin = 1; bin<nbins;bin++){
    //float scale = 2-0.2*bin;
     //((TH1D *)ratios[bin])->Scale(1.0/scale);
    etForLowPtTracksAntiProtons[bin]->Draw("same");
  }
  //((TH1D *)etForLowPtTracksAntiProtons[1])->SetMaximum(0.03);
  if(det.Contains("PHOS")) ((TH1D *)etForLowPtTracksAntiProtons[1])->SetMaximum(14);
  ((TH1D *)etForLowPtTracksAntiProtons[1])->SetMinimum(0.0);
  TString filename = "/tmp/FractionHadronsDeposted"+det+".eps";
  c9->SaveAs(filename.Data());

  return;

}
TH1* bayneseffdiv(TH1* numerator, TH1* denominator,Char_t* name) 
{
    if(!numerator){
      cerr<<"Error:  numerator does not exist!"<<endl;
      return NULL;
    }
    if(!denominator){
      cerr<<"Error:  denominator does not exist!"<<endl;
      return NULL;
    }
    TH1F* result = (TH1F*)numerator->Clone(name);
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

void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle){
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetYaxis()->SetTitle(ytitle);
  histo->Sumw2();
}

