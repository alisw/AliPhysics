Int_t colors[] = {TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		  TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		  TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		  TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack};
Int_t markers[] = {20,21,22,23,33, 24,25,26,32,27, 20,21,22,23,33, 24,25,26,32,27};
void SetStyles(TGraph *graph, Int_t marker, Int_t color){
  graph->SetMarkerStyle(marker);
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetMarkerSize(1.5);
}

Float_t npartShort[10] =    {382.8,329.7,260.5,186.4,128.9, 85,52.8,30.0,15.8,7.48};
Float_t npartErrShort[10] = {    6,    6,  4.4,  3.9,  3.3,2.6,   2, 1.3, 0.6,0.29};
Float_t npart[20] = {382.7, 329.4, 281.2, 239, 202.1, 169.5, 141, 116, 94.11, 75.3, 59.24, 45.58, 34.33, 25.21, 17.96, 12.58, 8.812, 6.158, 4.376, 3.064};
Float_t npartErr[20] = {3, 4.3, 4.1, 3.5, 3.3, 3.3, 3.1, 2.8, 2.6, 2.3, 1.8, 1.4, 1.1, 0.87, 0.66, 0.45, 0.26, 0.19, 0.1, 0.059};

//========================Charged Pion Reference========================================
//Arrays for defining comparison plots
Float_t pionPlusEt[10] = {360.7,298.3,223.8,149.9,96.1, 58.1,32.4,16.4,7.3,2.7};
Float_t pionMinusEt[10] ={363.7,300.4,225.4,150.5,96.6, 58.4,32.5,16.5,7.4,2.8};
Float_t pionEtError[10] = {19.3, 15.3,11.3 ,7.5  , 4.8,  2.9, 1.6, 0.8,0.4,0.1};
Float_t pionEt[10] = {0,0,0,0,0, 0,0,0,0,0};
Float_t ypion[10]  = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t ypionerr[10]  = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t emEtFromPions[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t emEtFromPionsErr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t emEtFromPionsPerNpart[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t emEtFromPionsPerNpartErr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t npartAlt1[20],npartAlt2[20],npartAlt3[20];

TGraphErrors *GetPionEtGraph(){
  for(int i=0;i<10;i++){
    pionEt[i] = (pionPlusEt[i]+pionMinusEt[i])/2.0;
    emEtFromPions[i] = pionEt[i]*1.085;
    emEtFromPionsErr[i] = emEtFromPions[i]*TMath::Sqrt(TMath::Power(0.030/1.085,2)+TMath::Power(pionEtError[i]/pionEt[i],2));
    ypion[i] = pionEt[i]/(npartShort[i]/2);
    ypionerr[i] = pionEtError[i]/(npartShort[i]/2);
    emEtFromPionsPerNpart[i] = emEtFromPions[i]/(npartShort[i]/2);
    emEtFromPionsPerNpartErr[i] = emEtFromPionsErr[i]/(npartShort[i]/2);
    npartAlt1[i] = npartShort[i]+2;
    npartAlt2[i] = npartShort[i]-2;
    npartAlt3[i] = npartShort[i]+4;
  }
  TGraphErrors *gr2 = new TGraphErrors(10,npartShort,ypion,npartErrShort,ypionerr);
  gr2->GetYaxis()->SetTitle("dE_{T}/d#eta#frac{1}{0.5*N_{part}} [GeV]");
  gr2->GetXaxis()->SetTitle("N_{part}");
  gr2->SetTitle("");
  gr2->GetXaxis()->SetRangeUser(0, 400);
  SetStyles(gr2,30,TColor::kBlue);

  return gr2;
}

TGraphErrors *GetPionEmEtGraph(){
    TGraphErrors *gr3 = new TGraphErrors(10,npartAlt3,emEtFromPionsPerNpart,npartErrShort,emEtFromPionsPerNpartErr);
    SetStyles(gr3,29,TColor::kBlue);
    return gr3;

}
//========================Reading in corrections========================================
Float_t nonLinError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t nonLinErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t efficiencyError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t efficiencyErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronErrorShort[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0};
Float_t neutronCorrShort[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0};
Float_t neutronErrorPerNChShort[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0};
Float_t neutronCorrPerNChShort[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0};
Float_t neutronErrorPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t neutronCorrPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorrPerNChShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryErrorPerNChShort[10] =  {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorrPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryErrorPerNCh[20] =  {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonCorrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonErrorPerNChShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonCorrPerNChShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonErrorPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonCorrPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadCorrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadronErrorPerNChShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadronCorrPerNChShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadronErrorPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t hadronCorrPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t minEtErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t minEtCorrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t minEtError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t minEtCorr[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};

Float_t kaonYield[10] = {109,90.5,68,46,30,18.2,10.2,5.1,2.3,0.855};
Float_t kaonYieldStatErr[10] = {0.3,0.2,0.1,0.1,0.1, 0.06,0.04,0.03,0.02,0.01};
Float_t kaonYieldSysErr[10] = {9,7,5,4,2, 1.5,0.8,0.4,0.2,0.09};
Float_t kaonYieldTotErr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonYieldPerNCh[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonYieldPerNChErr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonEtPerNCh[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t kaonEtPerNChErr[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Double_t kaonPlusEt[2][10] = {{91.7712,75.8971,56.841,37.6962,23.9923,14.255,7.73469,3.74477,1.60505,0.578278},{6.61811,5.39337,3.9978,2.6337,1.67854,1.01849,0.557879,0.278199,0.125057,0.0592682}};
Double_t kaonMinusEt[2][10] = {{90.4723,74.9444,55.9463,37.286,23.6591,14.0413,7.63067,3.69337,1.59219,0.571019},{7.01588,5.76588,4.20933,2.80388,1.77983,1.06934,0.588003,0.292737,0.138191,0.0600075}};

void ReadMinEtCorrections(){
  cout<<"Reading in min et corrections..."<<endl;
  string inline;
  float value = 0;
  float error = 0;
  int i=0;
  TString minetInfileNameShort = "MinEt"+detector+"Short.dat";
  ifstream myminetfileShort (minetInfileNameShort.Data());
  if (myminetfileShort.is_open()){
    while ( myminetfileShort.good() )
      {
	getline (myminetfileShort,inline);
	istringstream tmp(inline);
	tmp >> value;
	tmp >> error;
	if(i<10){
	  minEtCorrShort[i] = value;
	  minEtErrorShort[i] = error;
	}
	//cout<<"min et corr cb "<<i<<" "<< minEtCorrShort[i]<<" +/- "<<minEtErrorShort[i]<<endl;
	i++;
      }
    myminetfileShort.close();
  }
  //doing the linear interpolation between the two data points we have
  TGraphErrors *graphMinEtCorrectionShort = GetMinEtCorrectionGraphShort();
  int shortbin = 0;
  for(int i=0;i<19;i++){
    //cout<<"Long bin "<<i<<" short bin "<<shortbin;
    if(i<2){//we have exact numbers so we don't need to interpolate
      minEtCorr[i] = minEtCorrShort[i];
      minEtError[i] = minEtErrorShort[i];
      shortbin++;
    }
    else{
      minEtCorr[i] = graphMinEtCorrectionShort->Eval(trackmultiplicity[i]);
      int altbin = shortbin-1;
      if(i%2==1){altbin = shortbin+1;}
      //cout<<" altbin "<<altbin;
      if(minEtErrorShort[shortbin]>minEtErrorShort[altbin]) minEtError[i] = minEtErrorShort[shortbin];
      else{minEtError[i] = minEtErrorShort[altbin];}
      if(i%2==1 && shortbin<10){shortbin++;}
    }
    // cout<<"min et corr cb "<<i<<" "<< minEtCorr[i]<<" +/- "<<minEtError[i]<<endl;
    //cout<<endl;
  }
  delete graphMinEtCorrectionShort;

}
void ReadInNeutronCorrections(){
  cout<<"Reading in neutron corrections..."<<endl;
  TString neutronInfileName = "Neutrons"+detector+".dat";
  ifstream myneutronfile (neutronInfileName.Data());
  string inline;
  float value = 0;
  float error = 0;
  int i=0;
  if (myneutronfile.is_open()){
    while ( myneutronfile.good() )
      {
	getline (myneutronfile,inline);
	istringstream tmp(inline);
	tmp >> value;
	tmp >> error;
	if(i<20){
	  neutronCorr[i] = value;
	  neutronError[i] = error;
	  if(trackmultiplicity[i]>0){
	    neutronCorrPerNCh[i] = value/(trackmultiplicity[i]);
	    neutronErrorPerNCh[i] = error/(trackmultiplicity[i]);
	  }
	}

	// cout<<"neutroncorr cb "<<i<<" "<<value<<" +/- "<<error<<endl;
	i++;
      }
    myneutronfile.close();

  }
  TString neutronInfileNameShort = "Neutrons"+detector+"Short.dat";
  ifstream myneutronfileShort (neutronInfileNameShort.Data());
  i=0;
  if (myneutronfileShort.is_open()){
    while ( myneutronfileShort.good() )
      {
	getline (myneutronfileShort,inline);
	istringstream tmp(inline);
	tmp >> value;
	tmp >> error;
	if(i<10){
	  neutronCorrShort[i] = value;
	  neutronErrorShort[i] = error;
	  neutronCorrPerNChShort[i] = value/(trackmultiplicityShort[i]);
	  neutronErrorPerNChShort[i] = error/(trackmultiplicityShort[i]);
	}
	//cout<<"neutroncorr cb "<<i<<" "<<neutronCorrShort[i]<<" +/- "<<neutronErrorShort[i]<<endl;
	i++;
      }
    myneutronfileShort.close();
  }

}
void ReadInSecondaryCorrections(){
  cout<<"Reading in secondary corrections..."<<endl;

    TString secondaryInfileName = "Secondaries"+detector+".dat";
    ifstream mysecondaryfile (secondaryInfileName.Data());
    string inline;
    float value = 0;
    float error = 0;
    int i=0;
    if (mysecondaryfile.is_open()){
      while ( mysecondaryfile.good() )
	{
	  getline (mysecondaryfile,inline);
	  istringstream tmp(inline);
	  tmp >> value;
	  tmp >> error;
	  if(i<20){
	    secondaryCorr[i] = value;
	    secondaryError[i] = error;
	    if(trackmultiplicity[i]>0){
	      secondaryCorrPerNCh[i] = value/(trackmultiplicity[i]);
	      secondaryErrorPerNCh[i] = error/(trackmultiplicity[i]);
	    }
	  }
	  //cout<<"secondarycorr cb "<<i<<" "<<value<<" +/- "<<error<<endl;
	  i++;
	}
        mysecondaryfile.close();
    }
    TString secondaryShortInfileName = "Secondaries"+detector+"Short.dat";
    ifstream mysecondaryShortfile (secondaryShortInfileName.Data());
    i=0;
    if (mysecondaryShortfile.is_open()){
      while ( mysecondaryShortfile.good() && i<10 )
	{
	  getline (mysecondaryShortfile,inline);
	  istringstream tmp(inline);
	  tmp >> value;
	  tmp >> error;
	  if(i<10){
	    secondaryCorrShort[i] = value;
	    secondaryErrorShort[i] = error;
	    secondaryCorrPerNChShort[i] = value/(trackmultiplicityShort[i]);
	    secondaryErrorPerNChShort[i] = error/(trackmultiplicityShort[i]);
	  }
	  //cout<<"secondarycorr cb "<<i<<" "<<value<<" +/- "<<error<<endl;
 	  i++;
	}
       mysecondaryShortfile.close();
    }

}
void ReadInKaonCorrections(){
  cout<<"Reading in kaon corrections..."<<endl;
  //junk.PHOS.CutNum6.txt
  TString kaonInfileName = "../spectrafits/KaonCut7EMCal.dat";
  if(detector.Contains("P")){
    kaonInfileName = "../spectrafits/KaonCut6PHOS.dat";
     }
    ifstream mykaonfile (kaonInfileName.Data());
    string inline;
    float value = 0;
    float error = 0;
    int i=0;
    if (mykaonfile.is_open()){
      while ( mykaonfile.good() )
	{
	  getline (mykaonfile,inline);
	  istringstream tmp(inline);
	  tmp >> value;
	  tmp >> error;
	  if(i<10){
	    kaonCorrShort[i] = value;
	    kaonErrorShort[i] = error;
	    kaonCorrPerNChShort[i] = value/(trackmultiplicityShort[i]);
	    kaonErrorPerNChShort[i] = error/(trackmultiplicityShort[i]);
	    //cout<<"kaoncorr cb "<<i<<" "<<value<<" +/- "<<error<<" "<<kaonCorrPerNPartShort[i]<<" +/- "<<kaonErrorPerNPartShort[i]<<endl;
	  }
	  i++;
	}
        mykaonfile.close();
    }


  TGraphErrors *graphKaonCorrectionShort = GetKaonCorrectionGraphShort();
  int shortbin = 0;
  for(int i=0;i<19;i++){
    //cout<<"Long bin "<<i<<" short bin "<<shortbin;
    if(i<2){//we have exact numbers so we don't need to interpolate
      kaonCorr[i] = kaonCorrShort[i];
      kaonError[i] = kaonErrorShort[i];
      shortbin++;
    }
    else{
      kaonCorr[i] = graphKaonCorrectionShort->Eval(trackmultiplicity[i]) * trackmultiplicity[i];
      int altbin = shortbin-1;
      if(i%2==1){altbin = shortbin+1;}
      //cout<<" altbin "<<altbin;
      if(kaonErrorPerNChShort[shortbin]>kaonErrorPerNChShort[altbin]) kaonError[i] = kaonErrorPerNChShort[shortbin] * trackmultiplicity[i];
      else{kaonError[i] =  kaonErrorPerNChShort[altbin] * trackmultiplicity[i];}
      if(i%2==1 && shortbin<10){shortbin++;}
    }
    //cout<<"kaoncorr cb "<<i<<" "<<kaonCorr[i]<<" +/- "<<kaonError[i]<<endl;
    kaonCorrPerNCh[i] = kaonCorr[i]/(trackmultiplicity[i]);
    kaonErrorPerNCh[i] = kaonError[i]/(trackmultiplicity[i]);
    //cout<<"min et corr cb "<<i<<" "<< kaonCorr[i]<<" +/- "<<kaonError[i]<<endl;
    //cout<<endl;
  }
  delete graphKaonCorrectionShort;


    for(int i=0;i<10;i++){
      kaonYieldTotErr[i] = TMath::Sqrt(kaonYieldStatErr[i]*kaonYieldStatErr[i]+kaonYieldSysErr[i]*kaonYieldSysErr[i]);
      kaonYieldPerNCh[i] = kaonYield[i]/(trackmultiplicityShort[i])/5;
      kaonYieldPerNChErr[i] = kaonYieldTotErr[i]/(trackmultiplicityShort[i])/5;
      float total = kaonPlusEt[0][i]+kaonMinusEt[0][i];
      float err = kaonPlusEt[1][i]+kaonMinusEt[1][i];
      kaonEtPerNCh[i] = total/(trackmultiplicityShort[i])/10;
      kaonEtPerNChErr[i] = err/(trackmultiplicityShort[i])/10;
    }


}
void CalculateHadronicCorrectionForOneBin(Int_t centbin1, Int_t centbin2, Bool_t isPhos, Bool_t isOver500MeV, Float_t &correction, Float_t &error){

    fHistMatchedTracksEvspTvsCentEffCorr->GetZaxis()->SetRange(centbin1,centbin2);
    fHistMatchedTracksEvspTvsCentEffCorr500MeV->GetZaxis()->SetRange(centbin1,centbin2);
    fHistMatchedTracksEvspTvsCent->GetZaxis()->SetRange(centbin1,centbin2);
    TH1D *dataEffCorrTmp = NULL;
    TH1D *dataEffCorrTmp2 = NULL;
    TH1D *dataTmp = NULL;
    TH1D *foundTmp = NULL;
    TH1D *notfoundTmp = NULL;
    dataTmp = (TH1D*)fHistMatchedTracksEvspTvsCent->Project3D("y");
    dataTmp->SetName(Form("dataTmp%i",centbin1));
    if(isOver500MeV){
      dataEffCorrTmp =(TH1D*) fHistMatchedTracksEvspTvsCentEffCorr500MeV->Project3D("y");
      dataEffCorrTmp2 =(TH1D*) fHistMatchedTracksEvspTvsCentEffCorr->Project3D("y");
      foundTmp = fHistFoundHadronsvsCent500MeV->ProjectionX(Form("Found%iTmp",centbin1),centbin1,centbin2);
      notfoundTmp = fHistNotFoundHadronsvsCent500MeV->ProjectionX(Form("NotFound%iTmp",centbin1),centbin1,centbin2);
    }
    else{
      dataEffCorrTmp = (TH1D*)fHistMatchedTracksEvspTvsCentEffCorr->Project3D("y");
      dataEffCorrTmp2 = dataEffCorrTmp;
      foundTmp = fHistFoundHadronsvsCent->ProjectionX(Form("Found%iTmp",centbin1),centbin1,centbin2);
      notfoundTmp = fHistNotFoundHadronsvsCent->ProjectionX(Form("NotFound%iTmp",centbin1),centbin1,centbin2);
    }
    float nfound = foundTmp->GetMean();//fHistFoundHadronsvsCent->GetBinContent(bin);
    float nnotfound = notfoundTmp->GetMean();//fHistNotFoundHadronsvsCent->GetBinContent(bin);
    float scaleLow = 0;
    float scaleHigh = 0;
    if(centbin1>=refBin){//for peripheral
      scaleHigh = 1.01;
      scaleLow = 0.97;
    }
    else{
      float scale1 = 1.0;
      float scale2 = 1.0;
      float refData = ((TH1D*)data[refBin])->GetMean();
      float myData = ((TH1D*)dataTmp)->GetMean();
      float refDataEffCorr = ((TH1D*)dataEffCorr[refBin])->GetMean();
      float myDataEffCorr = ((TH1D*)dataEffCorrTmp2)->GetMean();
      if(TMath::Abs(myData)>1e-5) scale1 = refData/myData;
      if(TMath::Abs(myDataEffCorr)>1e-5) scale2 = refDataEffCorr/myDataEffCorr;
      if(scale1<scale2){
	scaleLow = 0.97*scale1;
	scaleHigh = scale2;
      }
      else{
	scaleLow = 0.97*scale2;
	scaleHigh = scale1;
      }
    }
    float myavg = dataEffCorrTmp->GetMean();
    float avg = (scaleLow+scaleHigh)/2.0*myavg;
    float err = TMath::Abs((scaleLow-scaleHigh))/2.0*myavg;
    if(TMath::Abs(avg)<1e-3){
      avg = 1e-3;
      cerr<<"WARNING:  ERROR NOT CALCULATED CORRECTLY!!"<<endl;//prevents a crash
    }
    //factor is the fraction to reduce the track-matched ET by to get the true background ET
    //corrfac is the factor to multiply by in order to get the fraction of hadrons leaving deposits which come from low pT
    float percentEfficiencyError = 0.01;
    float  factor = 1-0.04;
    float corrfac = 1.275-1;
    float corrfacerr = 0.059 ;
    float eLowAverage = avg;
    float eLowAverageErr = err;
    if(isPhos){
      factor = 1-0.03;
      corrfac = 1.300-1;
      corrfacerr = 0.065;
    }
    if(isOver500MeV){
      eLowAverage = 1.0;
      eLowAverageErr = 0.05;
      if(isPhos){
	//fraction ranges from 5% - 26%
	corrfac = ((1.0/0.95 + 1.0/0.74)/2.0 - 1)*eLowAverage/avg;
	float corrfacerrtmp = (TMath::Abs((1.0/0.95 - 1.0/0.74)/2.0))*eLowAverage/avg;
	corrfacerr = corrfac * TMath::Sqrt(TMath::Power(corrfacerrtmp/corrfac,2)+TMath::Power(eLowAverageErr/eLowAverage,2));
      }
      else{
	//fraction ranges from 5% - 17%
	corrfac = ((1.0/0.95 + 1.0/0.83)/2.0 - 1)*eLowAverage/avg;
	float corrfacerrtmp = (TMath::Abs((1.0/0.95 - 1.0/0.83)/2.0))*eLowAverage/avg;
	corrfacerr = corrfac * TMath::Sqrt(TMath::Power(corrfacerrtmp/corrfac,2)+TMath::Power(eLowAverageErr/eLowAverage,2));
      }
    }
    //the energy from low pT is the fraction of tracks that come from low pT tracks times the average energy of these tracks
    float eLow = corrfac * (nfound+nnotfound)*eLowAverage;
    float eLowErr = TMath::Sqrt(TMath::Power(corrfacerr*(nfound+nnotfound)*eLowAverage,2)+TMath::Power(eLowAverageErr*corrfac* (nfound+nnotfound),2)+TMath::Power(eLow*percentEfficiencyError,2));//error on the hadronic correction

    float eNotFound = nnotfound*avg;
    float eNotFoundErr = TMath::Sqrt(TMath::Power(err*nnotfound,2)+TMath::Power(percentEfficiencyError*eNotFound,2));//error on the hadronic correction

    float y = (corrfac * (nfound+nnotfound) +nnotfound)*avg;
    float finalerr = TMath::Sqrt(TMath::Power(corrfacerr*(nfound+nnotfound)*avg,2)+err*err*(TMath::Power(corrfac* (nfound+nnotfound),2)+nnotfound*nnotfound)+TMath::Power(percentEfficiencyError*y,2));//error on the hadronic correction
    correction = y;
    error = finalerr;
    delete dataEffCorrTmp;
    delete foundTmp;
    delete notfoundTmp;
}
void CalculateHadronCorrections(Bool_t isPhos){
  float plotscale = 5.0;
  for(int i=0;i<19;i++){
    Float_t correction = 0;
    Float_t error = 0;
    CalculateHadronicCorrectionForOneBin(i+1,i+1,isPhos,kFALSE,correction,error);
    hadCorr[i] = correction;//hadCorrEmcal[i];
    hadError[i] = error;//hadErrorEmcal[i];
    hadronCorrPerNCh[i] = correction/(trackmultiplicity[i])/plotscale;//hadCorrEmcal[i];
    hadronErrorPerNCh[i] = error/(trackmultiplicity[i])/plotscale;//hadErrorEmcal[i];
    //cout<<"had cor "<<i<<" "<<correction<<" +/- "<<error<< "  "<<  correction/(trackmultiplicity[i])<< " +/- "<<  error/(trackmultiplicity[i]) <<endl;
  }

  int j=0;
  for(int i=0;i<10;i++){
    int centbinlow = i+1;
    int centbinhigh = i+1;
    if(i<2){//These bins are exactly what they should bin in the 20 bin binning
      j++;//i=0 j=0; i=1 j=1
      centbinlow = j;
      centbinhigh = j;
    }
    else{
      centbinlow = j+1;
      centbinhigh = j+2;
      j+=2;
    }
    Float_t correction = 0;
    Float_t error = 0;
    CalculateHadronicCorrectionForOneBin(centbinlow,centbinhigh,isPhos,kFALSE,correction,error);
    hadCorrShort[i] = correction;//hadCorrEmcal[i];
    hadErrorShort[i] = error;//hadErrorEmcal[i];
    hadronCorrPerNChShort[i] = correction/(trackmultiplicityShort[i])/plotscale;//hadCorrEmcal[i];
    hadronErrorPerNChShort[i] = error/(trackmultiplicityShort[i])/plotscale;//hadErrorEmcal[i];
    //cout<<"had cor "<<i<<" "<<correction<<" +/- "<<error<< "  "<<  correction/(trackmultiplicityShort[i])<< " +/- "<<  error/(trackmultiplicityShort[i]) <<endl;
  }
}

TGraphErrors *GetMinEtCorrectionGraphShort(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,secondaryCorrPerNPartShort,npartErrShort,secondaryErrorPerNPartShort);
  TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,minEtCorrShort,trackmultiplicityShortError,minEtErrorShort);
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,minEtCorrShort,npartErrShort,minEtErrorShort);
  SetStyles(gr3,29,TColor::kGreen+3);
    return gr3;

}
TGraphErrors *GetMinEtCorrectionGraph(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,secondaryCorrPerNPartShort,npartErrShort,secondaryErrorPerNPartShort);
  TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,minEtCorr,trackmultiplicityError,minEtError);
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,minEtCorrShort,npartErrShort,minEtErrorShort);
  SetStyles(gr3,30,TColor::kGreen+3);
    return gr3;

}


TGraphErrors *GetSecondaryCorrectionGraph(){
  TGraphErrors *gr3 = new TGraphErrors(18,trackmultiplicity,secondaryCorrPerNCh,trackmultiplicityError,secondaryErrorPerNCh);
  SetStyles(gr3,30,TColor::kGreen+3);
    return gr3;

}
TGraphErrors *GetSecondaryCorrectionGraphShort(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,secondaryCorrPerNPartShort,npartErrShort,secondaryErrorPerNPartShort);
  TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,secondaryCorrPerNChShort,trackmultiplicityShortError,secondaryErrorPerNChShort);
  SetStyles(gr3,29,TColor::kGreen+3);
    return gr3;

}
TGraphErrors *GetNeutronCorrectionGraph(){
    TGraphErrors *gr3 = new TGraphErrors(18,trackmultiplicity,neutronCorrPerNCh,trackmultiplicityError,neutronErrorPerNCh);
    SetStyles(gr3,24,TColor::kBlue);
    return gr3;

}
TGraphErrors *GetNeutronCorrectionGraphShort(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,neutronCorrPerNChShort,npartErrShort,neutronErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,neutronCorrPerNChShort,trackmultiplicityShortError,neutronErrorPerNChShort);
    SetStyles(gr3,20,TColor::kBlue);
    return gr3;

}
TGraphErrors *GetHadronCorrectionGraph(){
    TGraphErrors *gr3 = new TGraphErrors(18,trackmultiplicity,hadronCorrPerNCh,trackmultiplicityError,hadronErrorPerNCh);
    SetStyles(gr3,25,1);
    return gr3;

}
TGraphErrors *GetHadronCorrectionGraphShort(){
  //TGraphErrors *gr3 = new TGraphErrors(10,npart,neutronCorrPerNChShort,xpionerr,neutronErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,hadronCorrPerNChShort,trackmultiplicityShortError,hadronErrorPerNChShort);
    SetStyles(gr3,21,1);
    return gr3;

}
TGraphErrors *GetKaonCorrectionGraphShort(){
  //TGraphErrors *gr3 = new TGraphErrors(10,xpion,kaonCorrPerNChShort,xpionerr,kaonErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,kaonCorrPerNChShort,trackmultiplicityShortError,kaonErrorPerNChShort);
    SetStyles(gr3,33,TColor::kRed);
    return gr3;

}
TGraphErrors *GetKaonCorrectionGraph(){
  //TGraphErrors *gr3 = new TGraphErrors(10,xpion,kaonCorrPerNChShort,xpionerr,kaonErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,kaonCorrPerNCh,trackmultiplicityError,kaonErrorPerNCh);
    SetStyles(gr3,27,TColor::kRed);
    return gr3;

}
TGraphErrors *GetKaonGraph(){
  //TGraphErrors *gr3 = new TGraphErrors(10,xpion,kaonCorrPerNChShort,xpionerr,kaonErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,kaonYieldPerNCh,trackmultiplicityShortError,kaonYieldPerNChErr);
    SetStyles(gr3,33,TColor::kBlue);
    return gr3;
}
TGraphErrors *GetKaonEtGraph(){
  //TGraphErrors *gr3 = new TGraphErrors(10,xpion,kaonCorrPerNChShort,xpionerr,kaonErrorPerNChShort);
    TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,kaonEtPerNCh,trackmultiplicityShortError,kaonEtPerNChErr);
    SetStyles(gr3,27,TColor::kBlue);
    return gr3;
}


//=====================READ IN DATA===================================
Float_t arrayofzeros[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t trackmultiplicity[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t trackmultiplicityShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t clustermultiplicity[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t clustermultiplicityShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t trackmultiplicityError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t trackmultiplicityShortError[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t clustermultiplicityError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t clustermultiplicityShortError[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t matchedtrackmultiplicity[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t matchedtrackmultiplicityPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t matchedtrackmultiplicityPerNCl[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t notmatchedtrackmultiplicity[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t notmatchedtrackmultiplicityPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t notmatchedtrackmultiplicityPerNCl[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totaltrackmultiplicityPerNCh[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t totaltrackmultiplicityPerNCl[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNChValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNChValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNChError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNChErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNClValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNClValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNClError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNClErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNPartPairValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNPartPairValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNPartPairError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t partialCorrEtPerNPartPairErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};

Float_t rawEtNoEffCorrValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtNoEffCorrError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtNoEffCorrPerNChValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtNoEffCorrPerNChError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtNoEffCorrPerNClValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtNoEffCorrPerNClError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrPerNChValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrPerNChError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrPerNClValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllNoEffCorrPerNClError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllPerNChValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllPerNChError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllPerNClValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t rawEtAllPerNClError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};

Float_t corrEtValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtPerNPartPairValues[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtPerNPartPairValuesShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtPerNPartPairError[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t corrEtPerNPartPairErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
TH2F *fHistNominalRawEt;
TH2F *fHistNominalNonLinLowEt;
TH2F *fHistNominalNonLinHighEt;
TH2F *fHistNominalEffLowEt;
TH2F *fHistNominalEffHighEt;
TH2F *fHistTotRawEt;
TH2F *fHistTotAllRawEt;
TH2F *fHistTotRawEtNoEffCorr;
TH2F *fHistTotAllRawEtNoEffCorr;
TH2F *fHistTotRawEt500MeV;
TH3F  *fHistMatchedTracksEvspTvsCent;
TH3F  *fHistMatchedTracksEvspTvsCentEffCorr;
TH3F  *fHistMatchedTracksEvspTvsCentEffCorr500MeV;
TH2F *fHistFoundHadronsvsCent;
TH2F *fHistNotFoundHadronsvsCent;
TH2F *fHistFoundHadronsvsCent500MeV;
Int_t refBin = 15;//Reference bin for scaling
TObjArray data(21);
TObjArray dataEffCorr(21);
TObjArray dataEffCorr500MeV(21);
TObjArray rawEt(21);
TObjArray rawEtShort(11);
TObjArray rawEtNoEffCorr(21);
TObjArray rawEtAll(21);
TObjArray rawEtAllNoEffCorr(21);

void ReadInData(char *filename,TString det){
  cout<<"Reading in data..."<<endl;


    TFile *f = TFile::Open(filename, "READ");
    if (!f)
    {
        std::cerr << "Could not open file: " << filename << std::endl;
    }

    TList *l = dynamic_cast<TList*>(f->Get("out1"));
    if (!l)
    {
        std::cerr << "Could not get object list from: " << filename << std::endl;
    }
    TString prefix = "fHistNominal";
    fHistNominalRawEt =(TH2F*) l->FindObject((prefix+"RawEt"+det+"Rec").Data());
    fHistNominalNonLinLowEt = (TH2F*)l->FindObject((prefix+"NonLinLowEt"+det+"Rec").Data());
    fHistNominalNonLinHighEt = (TH2F*)l->FindObject((prefix+"NonLinHighEt"+det+"Rec").Data());
    fHistNominalEffLowEt = (TH2F*)l->FindObject((prefix+"EffLowEt"+det+"Rec").Data());
    fHistNominalEffHighEt = (TH2F*)l->FindObject((prefix+"EffHighEt"+det+"Rec").Data());
    fHistTotRawEt =(TH2F*) l->FindObject("fHistTotRawEtEffCorr");
    fHistTotAllRawEt = (TH2F*) l->FindObject("fHistTotAllRawEtEffCorr");
    fHistTotRawEtNoEffCorr = (TH2F*)l->FindObject("fHistTotRawEt");
    fHistTotAllRawEtNoEffCorr = (TH2F*)l->FindObject("fHistTotAllRawEt");
    fHistTotRawEt500MeV =(TH2F*) l->FindObject("fHistTotRawEtEffCorr500MeV");

    fHistCentVsNchVsNcl =(TH3F*) l->FindObject("fHistCentVsNchVsNclReco");


    fHistMatchedTracksEvspTvsCent =(TH3F*) l->FindObject("fHistMatchedTracksEvspTvsCent");
    fHistMatchedTracksEvspTvsCentEffCorr =(TH3F*) l->FindObject("fHistMatchedTracksEvspTvsCentEffTMCorr");
    fHistMatchedTracksEvspTvsCentEffCorr500MeV = (TH3F*) l->FindObject("fHistMatchedTracksEvspTvsCentEffTMCorr500MeV");
    fHistFoundHadronsvsCent = (TH2F*)l->FindObject("fHistFoundHadronsvsCent");
    fHistNotFoundHadronsvsCent = (TH2F*)l->FindObject("fHistNotFoundHadronsvsCent");
    fHistFoundHadronsvsCent500MeV = (TH2F*)l->FindObject("fHistFoundHadronsvsCent500MeV");
    fHistNotFoundHadronsvsCent500MeV = (TH2F*)l->FindObject("fHistNotFoundHadronsvsCent500MeV");
    int nbins = 20;
    for(int bin = 1; bin<=nbins;bin++){
      fHistMatchedTracksEvspTvsCent->GetZaxis()->SetRange(bin,bin);
      data[bin] = fHistMatchedTracksEvspTvsCent->Project3D("y");
      fHistMatchedTracksEvspTvsCentEffCorr->GetZaxis()->SetRange(bin,bin);
      dataEffCorr[bin] = fHistMatchedTracksEvspTvsCentEffCorr->Project3D("y");
      fHistMatchedTracksEvspTvsCentEffCorr500MeV->GetZaxis()->SetRange(bin,bin);
      dataEffCorr500MeV[bin] = fHistMatchedTracksEvspTvsCentEffCorr500MeV->Project3D("y");
      ((TH1D*)data[bin])->SetName(Form("DataEff%i",bin));
      ((TH1D*)dataEffCorr[bin])->SetName(Form("DataEffCorr%i",bin));
      ((TH1D*)dataEffCorr500MeV[bin])->SetName(Form("DataEffCorr500MeV%i",bin));

      rawEt[bin]= fHistTotRawEt->ProjectionX(Form("RawEt%i",bin),bin,bin);
      ((TH1D*)rawEt[bin])->SetName(Form("rawEt%i",bin));
      partialCorrEtValues[bin-1] = (Float_t)((TH1D*)rawEt[bin])->GetMean();
      partialCorrEtError[bin-1] = (Float_t) ((TH1D*)rawEt[bin])->GetMeanError();
      partialCorrEtPerNPartPairValues[bin-1] = partialCorrEtValues[bin-1]/(npart[bin-1])/2.0*10;
      partialCorrEtPerNPartPairError[bin-1]  =  partialCorrEtError[bin-1]/(npart[bin-1])/2.0*10;

      rawEtNoEffCorr[bin]= fHistTotRawEtNoEffCorr->ProjectionX(Form("RawEtNoEffCorr%i",bin),bin,bin);
      ((TH1D*)rawEtNoEffCorr[bin])->SetName(Form("rawEtNoEffCorr%i",bin));
      rawEtNoEffCorrValues[bin-1] = (Float_t)((TH1D*)rawEtNoEffCorr[bin])->GetMean();
      rawEtNoEffCorrError[bin-1] = (Float_t) ((TH1D*)rawEtNoEffCorr[bin])->GetMeanError();

      rawEtAllNoEffCorr[bin]= fHistTotAllRawEtNoEffCorr->ProjectionX(Form("RawEtAllNoEffCorr%i",bin),bin,bin);
      ((TH1D*)rawEtAllNoEffCorr[bin])->SetName(Form("rawEtAllNoEffCorr%i",bin));
      rawEtAllNoEffCorrValues[bin-1] = (Float_t)((TH1D*)rawEtAllNoEffCorr[bin])->GetMean();
      rawEtAllNoEffCorrError[bin-1] = (Float_t) ((TH1D*)rawEtAllNoEffCorr[bin])->GetMeanError();

      rawEtAll[bin]= fHistTotAllRawEt->ProjectionX(Form("RawEtAll%i",bin),bin,bin);
      ((TH1D*)rawEtAll[bin])->SetName(Form("rawEtAll%i",bin));
      rawEtAllValues[bin-1] = (Float_t)((TH1D*)rawEtAll[bin])->GetMean();
      rawEtAllError[bin-1] = (Float_t) ((TH1D*)rawEtAll[bin])->GetMeanError();
      //cout<<"bin "<<bin<<" "<<partialCorrEtValues[bin-1]<<" "<<rawEtNoEffCorrValues[bin-1]<<" "<<rawEtAllNoEffCorrValues[bin-1]<<" "<<rawEtAllValues[bin-1]<<endl;

      TH1D *temp = fHistNominalRawEt->ProjectionX("temp",bin,bin);
      float nominal = temp->GetMean();
      //cout<<" Mean "<<temp->GetMean()<<" nbins "<<temp->GetNbinsX()<<endl;
      delete temp;
      temp = fHistNominalNonLinLowEt->ProjectionX("temp",bin,bin);
      float nonlinlow = temp->GetMean();
      delete temp;
      temp = fHistNominalNonLinHighEt->ProjectionX("temp",bin,bin);
      float nonlinhigh = temp->GetMean();
      delete temp;
      temp = fHistNominalEffLowEt->ProjectionX("temp",bin,bin);
      float efflow = temp->GetMean();
      delete temp;
      temp = fHistNominalEffHighEt->ProjectionX("temp",bin,bin);
      float effhigh = temp->GetMean();
      delete temp;
      float nonlinfracerr = 0;
      if(nonlinhigh >0 || nonlinlow >0) nonlinfracerr = TMath::Abs(nonlinhigh-nonlinlow)/(nonlinhigh+nonlinlow);
      float efffracerr = 0;
      if(effhigh >0 || efflow>0)efffracerr = TMath::Abs(effhigh-efflow)/(effhigh+efflow);
      //cout<<"cb "<<bin-1<<" nonlinerr "<<nonlinfracerr<<" efficiencyError "<<efffracerr<<endl;
      nonLinError[bin-1] = nonlinfracerr;
      efficiencyError[bin-1] = efffracerr;

      //(not)matchedtrackmultiplicity, (not)matchedtrackmultiplicityPerNCh, (not)matchedtrackmultiplicityPerNCl
      temp = fHistFoundHadronsvsCent->ProjectionX(Form("Found%iTmp",bin),bin,bin);
      matchedtrackmultiplicity[bin-1] = temp->GetMean();
      delete temp;
      temp = fHistNotFoundHadronsvsCent->ProjectionX(Form("NotFound%iTmp",bin),bin,bin);
      notmatchedtrackmultiplicity[bin-1] = temp->GetMean();
      delete temp;



    }


    for(int cb=0;cb<20;cb++){
      fHistCentVsNchVsNcl->GetXaxis()->SetRange(cb+1,cb+1);
      TH1D *trackmultiplicityHist = fHistCentVsNchVsNcl->Project3D("y");
      TH1D *clustermultiplicityHist = fHistCentVsNchVsNcl->Project3D("z");
      trackmultiplicity[cb] = (Float_t) trackmultiplicityHist->GetMean();
      clustermultiplicity[cb] = (Float_t)clustermultiplicityHist->GetMean();
      trackmultiplicityError[cb] = (Float_t) trackmultiplicityHist->GetMeanError();
      clustermultiplicityError[cb] = (Float_t)clustermultiplicityHist->GetMeanError();
      delete trackmultiplicityHist;
      delete clustermultiplicityHist;
    }
    int cb1 = 0;
    for(int cb=0;cb<10;cb++){
      int cb2 = cb1+1;
      if(cb1<2) cb2 = cb1;
      //cout<<"From "<<cb1<<" to "<<cb2<<endl;
      fHistCentVsNchVsNcl->GetXaxis()->SetRange(cb1+1,cb2+1);
      TH1D *trackmultiplicityHistShort = fHistCentVsNchVsNcl->Project3D("y");
      TH1D *clustermultiplicityHistShort = fHistCentVsNchVsNcl->Project3D("z");
      trackmultiplicityShort[cb] = (Float_t) trackmultiplicityHistShort->GetMean();
      clustermultiplicityShort[cb] = (Float_t)clustermultiplicityHistShort->GetMean();
      trackmultiplicityShortError[cb] = (Float_t) trackmultiplicityHistShort->GetMeanError();
      clustermultiplicityShortError[cb] = (Float_t)clustermultiplicityHistShort->GetMeanError();
      delete trackmultiplicityHistShort;
      delete clustermultiplicityHistShort;
      if(cb1<2) cb1++;
      else{cb1+=2;}
      rawEtShort[cb]= fHistTotRawEt->ProjectionX(Form("RawEtShort%i",bin),cb1+1,cb2+1);
      ((TH1D*)rawEtShort[cb])->SetName(Form("rawEtShort%i",cb));
      partialCorrEtValuesShort[cb] = (Float_t)((TH1D*)rawEtShort[cb])->GetMean();
      partialCorrEtErrorShort[cb] = (Float_t) ((TH1D*)rawEtShort[cb])->GetMeanError();

      partialCorrEtPerNChValuesShort[cb] = partialCorrEtValuesShort[cb]/(trackmultiplicityShort[cb])/2.0;
      partialCorrEtPerNChErrorShort[cb]  =  partialCorrEtErrorShort[cb]/(trackmultiplicityShort[cb])/2.0;

      partialCorrEtPerNPartPairValuesShort[cb] = partialCorrEtValues[cb]/(npart[cb])/2.0*10;
      partialCorrEtPerNPartPairErrorShort[cb]  =  partialCorrEtError[cb]/(npart[cb])/2.0*10;

      TH1D *temp = fHistNominalRawEt->ProjectionX("temp",cb1+1,cb2+1);
      float nominal = temp->GetMean();
      //cout<<" Mean "<<temp->GetMean()<<" nbins "<<temp->GetNbinsX()<<endl;
      delete temp;
      temp = fHistNominalNonLinLowEt->ProjectionX("temp",cb1+1,cb2+1);
      float nonlinlow = temp->GetMean();
      delete temp;
      temp = fHistNominalNonLinHighEt->ProjectionX("temp",cb1+1,cb2+1);
      float nonlinhigh = temp->GetMean();
      delete temp;
      temp = fHistNominalEffLowEt->ProjectionX("temp",cb1+1,cb2+1);
      float efflow = temp->GetMean();
      delete temp;
      temp = fHistNominalEffHighEt->ProjectionX("temp",cb1+1,cb2+1);
      float effhigh = temp->GetMean();
      delete temp;
      float nonlinfracerr = 0;
      if(nonlinhigh >0 || nonlinlow >0) nonlinfracerr = TMath::Abs(nonlinhigh-nonlinlow)/(nonlinhigh+nonlinlow);
      float efffracerr = 0;
      if(effhigh >0 || efflow>0)efffracerr = TMath::Abs(effhigh-efflow)/(effhigh+efflow);
      nonLinErrorShort[cb] = nonlinfracerr;
      efficiencyErrorShort[cb] = efffracerr;



    }

}
void ApplyCorrections(Float_t scale){//scale takes into account the acceptance in eta and phi and the 1/etaacc
  //correlation of errors:
    //hadCorr - calculated using matched tracks, correcting with tracking efficiency, systematic error dominated by uncertainty in energy deposited in calorimeter
  //kaon correction - calculated using kaon spectra, systematic errors from spectra systematic errors
  //neutron correction - systematic errors a bit of a fudge.
  //secondary correction - systematic errors from Nch vs Ncl scaling
  //etmin - calculated from kinematics, simulation, and pion spectra
  //efficiency error - determined from different material budgets
  //nonlinearity error - difference between test beam data and simulation
  //arguably kaon correction and etmin correction are somewhat correlated

  for(int cb = 0;cb<19;cb++){
    partialCorrEtPerNChValues[cb] = partialCorrEtValues[cb]/(trackmultiplicity[cb]);
    partialCorrEtPerNChError[cb]  =  partialCorrEtError[cb]/(trackmultiplicity[cb]);
    partialCorrEtPerNClValues[cb] = partialCorrEtValues[cb]/(clustermultiplicity[cb]);
    partialCorrEtPerNClError[cb]  =  partialCorrEtError[cb]/(clustermultiplicity[cb]);


    rawEtNoEffCorrPerNChValues[cb] = rawEtNoEffCorrValues[cb]/(trackmultiplicity[cb]);
    rawEtNoEffCorrPerNChError[cb]  =  rawEtNoEffCorrError[cb]/(trackmultiplicity[cb]);
    rawEtNoEffCorrPerNClValues[cb] = rawEtNoEffCorrValues[cb]/(clustermultiplicity[cb]);
    rawEtNoEffCorrPerNClError[cb]  =  rawEtNoEffCorrError[cb]/(clustermultiplicity[cb]);


    rawEtAllNoEffCorrPerNChValues[cb] = rawEtAllNoEffCorrValues[cb]/(trackmultiplicity[cb]);
    rawEtAllNoEffCorrPerNChError[cb]  =  rawEtAllNoEffCorrError[cb]/(trackmultiplicity[cb]);
    rawEtAllNoEffCorrPerNClValues[cb] = rawEtAllNoEffCorrValues[cb]/(clustermultiplicity[cb]);
    rawEtAllNoEffCorrPerNClError[cb]  =  rawEtAllNoEffCorrError[cb]/(clustermultiplicity[cb]);


    rawEtAllPerNChValues[cb] = rawEtAllValues[cb]/(trackmultiplicity[cb]);
    rawEtAllPerNChError[cb]  =  rawEtAllError[cb]/(trackmultiplicity[cb]);
    rawEtAllPerNClValues[cb] = rawEtAllValues[cb]/(clustermultiplicity[cb]);
    rawEtAllPerNClError[cb]  =  rawEtAllError[cb]/(clustermultiplicity[cb]);

    matchedtrackmultiplicityPerNCh[cb] = matchedtrackmultiplicity[cb]/(trackmultiplicity[cb]);
    notmatchedtrackmultiplicityPerNCh[cb] = notmatchedtrackmultiplicity[cb]/(trackmultiplicity[cb]);
    totaltrackmultiplicityPerNCh[cb] = matchedtrackmultiplicityPerNCh[cb] + notmatchedtrackmultiplicityPerNCh[cb];
    matchedtrackmultiplicityPerNCl[cb] = matchedtrackmultiplicity[cb]/(clustermultiplicity[cb]);
    notmatchedtrackmultiplicityPerNCl[cb] = notmatchedtrackmultiplicity[cb]/(clustermultiplicity[cb]);
    totaltrackmultiplicityPerNCl[cb] = matchedtrackmultiplicityPerNCl[cb] + notmatchedtrackmultiplicityPerNCl[cb];

      cout<<"cb "<<cb<<" "<<partialCorrEtValues[cb]<<" "<<rawEtNoEffCorrValues[cb]<<" "<<rawEtAllNoEffCorrValues[cb]<<" "<<rawEtAllValues[cb]<<endl;
      cout<<"cb "<<cb<<" "<<partialCorrEtPerNChValues[cb]<<" "<<rawEtNoEffCorrPerNChValues[cb]<<" "<<rawEtAllNoEffCorrPerNChValues[cb]<<" "<<rawEtAllPerNChValues[cb]<<endl;
      //cout<<"cb "<<cb<<" "<<partialCorrEtPerNChValues[cb]<<" "<< partialCorrEtValues[cb]<<" "<< partialCorrEtError[cb]<<" "<<trackmultiplicity[cb] <<endl;
    corrEtValues[cb] = scale*(partialCorrEtValues[cb] - hadCorr[cb] - kaonCorr[cb] - neutronCorr[cb] - secondaryCorr[cb])/minEtCorr[cb];
    //this comes up enough in the error calculations we'll just make a variable for it
    float partialEt = scale*(partialCorrEtValues[cb])/minEtCorr[cb];
    //add up the error squared
    float err = 0;
    float partialerr = 0;
    bool writeerror = false;
    if(writeerror)cout<<"partialEt "<<partialEt<<" err^2 = ";
    //Et min correction
    //partialerr += TMath::Power(minEtError[cb]/minEtCorr[cb]*corrEtValues[cb],2);
    partialerr = TMath::Power(minEtError[cb]/minEtCorr[cb]*corrEtValues[cb],2.0);
    if(writeerror)cout<<partialerr<<"+";
    err+=partialerr;
    //nonlinearity correction - this is saved as a fractional error
    partialerr = TMath::Power(nonLinError[cb]*partialEt,2.0);
    if(writeerror)cout<<partialerr<<"+";
    err+=partialerr;
    //efficiency correction - this is also saved as a fractional error
    partialerr = TMath::Power(efficiencyError[cb]*partialEt,2);
    if(writeerror)cout<<partialerr<<"+";
    err+=partialerr;
    //hadron correction
    partialerr = TMath::Power(hadError[cb]*scale/minEtCorr[cb],2);
    if(writeerror)cout<<partialerr<<"+";
    err+=partialerr;
    //neutron correction
    partialerr = TMath::Power(neutronError[cb]*scale/minEtCorr[cb],2);
    if(writeerror)cout<<partialerr<<"+";
    err+=partialerr;
    //kaon correction
    partialerr = TMath::Power(kaonError[cb]*scale/minEtCorr[cb],2);
    if(writeerror)cout<<partialerr<<"+";
    err+=partialerr;
    //secondary correction
    partialerr = TMath::Power(secondaryError[cb]*scale/minEtCorr[cb],2);
    if(writeerror)cout<<partialerr;
    err+=partialerr;
    //And take the square root
    err = TMath::Sqrt(err);
    if(writeerror)cout<<"... = "<<err<<endl;
    corrEtError[cb] = err;
    corrEtPerNPartPairValues[cb] = corrEtValues[cb]/(npart[cb]/2.0);
    corrEtPerNPartPairError[cb]  =  corrEtError[cb]/(npart[cb]/2.0);
    //cout<<"cb "<<cb <<" et "<< corrEtPerNPartPairValues[cb] <<" +/- "<<corrEtPerNPartPairError[cb];

    //cout<<" = "<<scale<<"*"<<"("<<partialCorrEtValues[cb]<<"+/-"<<partialCorrEtError[cb]<<" - "<<hadCorr[cb]<<"+/-"<<hadError[cb]<<" - "<<kaonCorr[cb]<<"+/-"<<kaonError[cb]<<" - "<<neutronCorr[cb]<<"+/-"<<neutronError[cb]<<" - "<<secondaryCorr[cb]<<"+/-"<<secondaryError[cb]<<")/"<<minEtCorr[cb]<<"+/-"<<minEtError[cb];
    //cout<<endl;
  }
//   for(int cb = 0;cb<10;cb++){
//   }

}

TGraphErrors *GetEtGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,npart,corrEtPerNPartPairValues,npartErr,corrEtPerNPartPairError);
    SetStyles(gr3,25,1);
    return gr3;
}
TGraphErrors *GetPartialCorrEtPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,partialCorrEtPerNChValues,trackmultiplicityError,partialCorrEtPerNChError);
    //TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,kaonEtPerNCh,trackmultiplicityShortError,kaonEtPerNChErr);

    SetStyles(gr3,25,1);
    return gr3;
}
TGraphErrors *GetPartialCorrEtPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,partialCorrEtPerNClValues,clustermultiplicityError,partialCorrEtPerNClError);
    //TGraphErrors *gr3 = new TGraphErrors(10,trackmultiplicityShort,kaonEtPerNCh,trackmultiplicityShortError,kaonEtPerNChErr);

    SetStyles(gr3,25,1);
    return gr3;
}
TGraphErrors *GetRawEtNoEffCorrPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,rawEtNoEffCorrPerNChValues,trackmultiplicityError,rawEtNoEffCorrPerNChError);
    SetStyles(gr3,21,1);
    return gr3;
}
TGraphErrors *GetRawEtNoEffCorrPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,rawEtNoEffCorrPerNClValues,clustermultiplicityError,rawEtNoEffCorrPerNClError);
    for(int i=0;i<20;i++){
      cout<<"i "<<i<<" "<<clustermultiplicity[i]<<": "<<rawEtNoEffCorrPerNClValues[i]<<"+/-"<<rawEtNoEffCorrPerNClError[i]<<endl;
    }
    SetStyles(gr3,21,1);
    return gr3;
}
TGraphErrors *GetRawEtAllNoEffCorrPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,rawEtAllNoEffCorrPerNChValues,trackmultiplicityError,rawEtAllNoEffCorrPerNChError);
    SetStyles(gr3,26,1);
    return gr3;
}
TGraphErrors *GetRawEtAllNoEffCorrPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,rawEtAllNoEffCorrPerNClValues,clustermultiplicityError,rawEtAllNoEffCorrPerNClError);
    SetStyles(gr3,26,1);
    return gr3;
}
TGraphErrors *GetRawEtAllPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,rawEtAllPerNChValues,trackmultiplicityError,rawEtAllPerNChError);
    SetStyles(gr3,26,1);
    return gr3;
}
TGraphErrors *GetRawEtAllPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,rawEtAllPerNClValues,clustermultiplicityError,rawEtAllPerNClError);
    SetStyles(gr3,22,1);
    return gr3;
}
TGraphErrors *GetMatchedTracksPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,matchedtrackmultiplicityPerNCh,trackmultiplicityError,arrayofzeros);
    SetStyles(gr3,20,1);
    return gr3;
}
TGraphErrors *GetMatchedTracksPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,matchedtrackmultiplicityPerNCl,clustermultiplicityError,arrayofzeros);
    SetStyles(gr3,20,1);
    return gr3;
}
TGraphErrors *GetTotalTracksPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,totaltrackmultiplicityPerNCl,clustermultiplicityError,arrayofzeros);
    SetStyles(gr3,34,1);
    return gr3;
}
TGraphErrors *GetTotalTracksPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,totaltrackmultiplicityPerNCh,trackmultiplicityError,arrayofzeros);
    SetStyles(gr3,34,1);
    return gr3;
}
TGraphErrors *GetNotMatchedTracksPerNChGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,trackmultiplicity,notmatchedtrackmultiplicityPerNCh,trackmultiplicityError,arrayofzeros);
    SetStyles(gr3,24,1);
    return gr3;
}
TGraphErrors *GetNotMatchedTracksPerNClGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,clustermultiplicity,notmatchedtrackmultiplicityPerNCl,clustermultiplicityError,arrayofzeros);
    SetStyles(gr3,24,1);
    return gr3;
}
TGraphErrors *GetPartialCorrEtPerNPartPairGraph(){
    TGraphErrors *gr3 = new TGraphErrors(20,npart,partialCorrEtPerNPartPairValues,npartErr,partialCorrEtPerNPartPairError);
    SetStyles(gr3,22,1);
    return gr3;
}
//=================================Plotting code====================================
Bool_t sim = false;
//Bool_t isPhos = kFALSE;
TString detector = "";
//void PlotEmEtVer2(TString filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.root")
void PlotEmEtVer2(Bool_t isPhos = kFALSE, Int_t cutset = 0)
{
  TString filename, simfilename;
  if(cutset==0){
    if(isPhos){
      filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.PHOS.LHC10hPass2.Run139465.root";
      simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS.LHC11a10a_bis.Run139465.root";
    }
    else{
      //filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.LooseTrackMatchCuts.root";
      //simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.LooseTrackMatchCuts.root";
      filename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.root";
      simfilename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root";
    }
  }
  TString detector = "Emcal";
  Float_t scale = 360.0/40.0/1.2;//Azimuthal acceptance over eta range
  if(filename.Contains("PHOS")){
    detector = "Phos";
    //isPhos = kTRUE;
    scale = 360.0/60/0.24;//Azimuthal acceptance over eta range
  }
  TString det = detector;
  //gROOT->LoadMacro("macros/PlotSecondariesCorr.C");
  //PlotSecondariesCorr(simfilename,filename);
  ReadInData(filename,detector);
  ReadInNeutronCorrections();
  ReadInSecondaryCorrections();
  ReadInKaonCorrections();
  ReadMinEtCorrections();
  CalculateHadronCorrections(isPhos);
  ApplyCorrections(scale);

  TCanvas *c1 = new TCanvas("c1","Corrections",600,400);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetLogx();
  TGraphErrors *graphKaonCorrectionShort = GetKaonCorrectionGraphShort();
  graphKaonCorrectionShort->Draw("AP");
  //graphKaonCorrectionShort->SetMaximim(0.02);
  //graphKaonCorrectionShort->GetYaxis()->SetRange(1,graphKaonCorrectionShort->GetYaxis()->FindBin(0.02));
  float emcalmax = 0.01;
  graphKaonCorrectionShort->GetHistogram()->SetMaximum(emcalmax);
  //set scales the same within naive geometric scaling
  if(isPhos) graphKaonCorrectionShort->GetHistogram()->SetMaximum(emcalmax*0.24/1.2*60.0/40.0);
  graphKaonCorrectionShort->GetXaxis()->SetTitle("N_{Ch}");
  graphKaonCorrectionShort->GetYaxis()->SetTitle("Correction/N_{Ch}");
  TGraphErrors *graphKaonCorrection = GetKaonCorrectionGraph();
  graphKaonCorrection->Draw("P same");
  TGraphErrors *graphSecondaryCorrectionShort = GetSecondaryCorrectionGraphShort();
  graphSecondaryCorrectionShort->Draw("P same");
  TGraphErrors *graphNeutronCorrectionShort = GetNeutronCorrectionGraphShort();
  graphNeutronCorrectionShort->Draw("P same");
  TGraphErrors *graphSecondaryCorrection = GetSecondaryCorrectionGraph();
  graphSecondaryCorrection->Draw("P same");
  TGraphErrors *graphNeutronCorrection = GetNeutronCorrectionGraph();
  graphNeutronCorrection->Draw("P same");
  TGraphErrors *graphHadronCorrection = GetHadronCorrectionGraph();
  graphHadronCorrection->Draw("P same");
  TGraphErrors *graphHadronCorrectionShort = GetHadronCorrectionGraphShort();
  graphHadronCorrectionShort->Draw("P same");
  TLegend *leg = new TLegend(0.607383,0.704301,0.810403,0.927419);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextSize(0.038682);
  leg->AddEntry(graphKaonCorrectionShort,"Kaon correction/N_{ch}","P");
  leg->AddEntry(graphNeutronCorrectionShort,"Neutron correction/N_{ch}","P");
  leg->AddEntry(graphSecondaryCorrectionShort,"Secondary correction/N_{ch}","P");
  leg->AddEntry(graphHadronCorrectionShort,"Hadron correction/5.0 N_{ch}","P");
  leg->Draw();
  TString corrplotname1 = "/tmp/CorrectionsVsNch"+detector+".png";
  c1->SaveAs(corrplotname1.Data());
//   TH1F *frame = new TH1F("frame","frame",1,0,2);
//   frame->GetYaxis()->SetTitle("dE_{T}/d#eta");
//   frame->GetXaxis()->SetTitle("N_{part}");
//   //fPion->SetRange(0,2);
//   frame->SetMinimum(0);
//   frame->SetMaximum(10);
  TCanvas *c2 = new TCanvas("c2","Min Et Correction",600,400);
  c2->SetTopMargin(0.02);
  c2->SetRightMargin(0.02);
  c2->SetBorderSize(0);
  c2->SetFillColor(0);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetFrameFillColor(0);
  c2->SetFrameBorderMode(0);
  c2->SetLogx();
  TGraphErrors *graphMinEtCorrectionShort = GetMinEtCorrectionGraphShort();
  graphMinEtCorrectionShort->GetXaxis()->SetTitle("f_{ETmin}");
  graphMinEtCorrectionShort->GetYaxis()->SetTitle("Correction/N_{Ch}");
//   graphMinEtCorrectionShort->GetHistogram()->SetMaximum(1.0);
//   graphMinEtCorrectionShort->GetHistogram()->SetMinimum(0.0);
  graphMinEtCorrectionShort->Draw("AP");
  graphKaonCorrectionShort->GetXaxis()->SetTitle("N_{Ch}");
  graphKaonCorrectionShort->GetYaxis()->SetTitle("f_{EtMin}");
  TGraphErrors *graphMinEtCorrection = GetMinEtCorrectionGraph();
  graphMinEtCorrection->Draw("P same");
  TString corrplotname2 = "/tmp/EtMinVsNch"+detector+".png";
  c2->SaveAs(corrplotname2.Data());

  TCanvas *c3 = new TCanvas("c3","Kaon Correction",600,400);
  c3->SetTopMargin(0.02);
  c3->SetRightMargin(0.02);
  c3->SetBorderSize(0);
  c3->SetFillColor(0);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetFrameFillColor(0);
  c3->SetFrameBorderMode(0);
  c3->SetLogx();

  TGraphErrors *graphKaonGraph = GetKaonGraph();
  graphKaonGraph->Draw("AP");
  //graphKaonGraph->GetHistogram()->SetMaximum();
  graphKaonGraph->GetHistogram()->SetMinimum(0.0);
  graphKaonGraph->GetXaxis()->SetTitle("N_{Ch}");
  graphKaonGraph->GetYaxis()->SetTitle("val/N_{Ch}");
  graphKaonCorrectionShort->Draw("P same");
  TGraphErrors *graphKaonEtGraph = GetKaonEtGraph();
  graphKaonEtGraph->Draw("P same");
  TLegend *leg3 = new TLegend(0.602349,0.362903,0.805369,0.586022);
  leg3->SetFillStyle(0);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.03);
  leg3->SetTextSize(0.038682);
  leg3->AddEntry(graphKaonGraph,"N_{K^{#pm}}/N_{Ch}","P");
  leg3->AddEntry(graphKaonCorrectionShort,"E_{T}^{kaon}","P");
  leg3->AddEntry(graphKaonEtGraph,"E_{T}^{K^{#pm}}/N_{Ch}","P");
  leg3->Draw();
  TString corrplotname3 = "/tmp/KaonVsNch"+detector+".png";
  c3->SaveAs(corrplotname3.Data());


  TCanvas *c4 = new TCanvas("c4", "dE_{T}/d#eta#frac{1}{0.5*N_{part}} [GeV]",700, 600);
  c4->SetTopMargin(0.02);
  c4->SetRightMargin(0.02);
  c4->SetBorderSize(0);
  c4->SetFillColor(0);
  c4->SetFillColor(0);
  c4->SetBorderMode(0);
  c4->SetFrameFillColor(0);
  c4->SetFrameBorderMode(0);
  TGraphErrors *graphEt = GetEtGraph();
  graphEt->GetHistogram()->GetXaxis()->SetTitle("N_{part}");
  graphEt->GetHistogram()->GetYaxis()->SetTitle("(E_{T}/(N_{part}/2)");
  graphEt->GetHistogram()->SetMinimum(0.0);
  graphEt->GetHistogram()->SetMaximum(3.5);
  graphEt->Draw("AP");
  TGraphErrors *graphPartialCorr = GetPartialCorrEtPerNPartPairGraph();
  //graphPartialCorr->Draw("P same");
  TGraphErrors *graphPionEt = GetPionEtGraph();
  graphPionEt->Draw("P same");
  TGraphErrors *graphPionEmEt = GetPionEmEtGraph();
  graphPionEmEt->Draw("P same");
  TLegend *leg2 = new TLegend(0.607383,0.704301,0.810403,0.927419);
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->SetTextSize(0.038682);
  leg2->AddEntry(graphEt,"E_{T}^{EM}","P");
  leg2->AddEntry(graphPionEt,"E_{T}^{#pi^{#pm}}","P");
  leg2->AddEntry(graphPionEmEt,"E_{T}^{EM calc. from #pi}","P");
  leg2->Draw();
  TString corrplotname4 = "/tmp/EtVsNpart"+detector+".png";
  c4->SaveAs(corrplotname4.Data());

  TCanvas *c5 = new TCanvas("c5", "Raw Et Vs NCh",1200, 400);
  c5->SetTopMargin(0.02);
  c5->SetRightMargin(0.02);
  c5->SetBorderSize(0);
  c5->SetFillColor(0);
  c5->SetFillColor(0);
  c5->SetBorderMode(0);
  c5->SetFrameFillColor(0);
  c5->SetFrameBorderMode(0);
  c5->Divide(2);
  c5->cd(1);
  TGraphErrors *partialCorrEtPerNCh = GetPartialCorrEtPerNChGraph();
  partialCorrEtPerNCh->GetHistogram()->SetMaximum(0.1*7.5/scale);//so PHOS and EMCal scales are similar...  partialCorrEtPerNCh->GetHistogram()->SetMinimum(0.0);
  partialCorrEtPerNCh->GetHistogram()->GetXaxis()->SetTitle("N_{Ch}");
  partialCorrEtPerNCh->GetHistogram()->GetYaxis()->SetTitle("(E_{T}/#epsilon f_{nonlin})/N_{Ch}");
  partialCorrEtPerNCh->Draw("AP");
  TGraphErrors *graphRawNoEffCorrCh = GetRawEtNoEffCorrPerNChGraph();
  graphRawNoEffCorrCh->Draw("P same");
  TGraphErrors *graphRawAllNoEffCorrCh = GetRawEtAllNoEffCorrPerNChGraph();
  graphRawAllNoEffCorrCh->Draw("P same");
  TGraphErrors *graphRawAllCh = GetRawEtAllPerNChGraph();
  graphRawAllCh->Draw("P same");
  TLegend *leg5 = new TLegend(0.416667,0.173077,0.62069,0.396853);
  leg5->SetFillStyle(0);
  leg5->SetFillColor(0);
  leg5->SetBorderSize(0);
  leg5->SetTextSize(0.03);
  leg5->SetTextSize(0.038682);
  leg5->AddEntry(partialCorrEtPerNCh,"E_{T}#delta_{TM}/#epsilon f_{nonlin}","P");
  leg5->AddEntry(graphRawNoEffCorrCh,"E_{T}#delta_{TM}/ f_{nonlin}","P");
  leg5->AddEntry(graphRawAllCh,"E_{T}/#epsilon f_{nonlin}","P");
  leg5->AddEntry(graphRawAllNoEffCorrCh,"E_{T}/f_{nonlin}","P");
  leg5->Draw();
  TString corrplotname5 = "/tmp/RawEtVsNch"+detector+".png";
  c5->SaveAs(corrplotname5.Data());

//   TCanvas *c6 = new TCanvas("c6", "Raw Et Vs NCl",700, 600);
//   c6->SetTopMargin(0.02);
//   c6->SetRightMargin(0.02);
//   c6->SetBorderSize(0);
//   c6->SetFillColor(0);
//   c6->SetFillColor(0);
//   c6->SetBorderMode(0);
//   c6->SetFrameFillColor(0);
//   c6->SetFrameBorderMode(0);

  c5->cd(2);
  TGraphErrors *graphRawAllCl = GetRawEtAllPerNClGraph();
  graphRawAllCl->GetHistogram()->GetXaxis()->SetTitle("N_{Ch}");
  graphRawAllCl->GetHistogram()->GetYaxis()->SetTitle("(E_{T}/#epsilon f_{nonlin})/N_{Ch}");
  graphRawAllCl->Draw("AP");
  TGraphErrors *partialCorrEtPerNCl = GetPartialCorrEtPerNClGraph();
  partialCorrEtPerNCl->GetHistogram()->GetXaxis()->SetTitle("N_{Ch}");
  partialCorrEtPerNCl->GetHistogram()->GetYaxis()->SetTitle("(E_{T}/#epsilon f_{nonlin})/N_{Ch}");
  partialCorrEtPerNCl->Draw("P same");
  TGraphErrors *graphRawNoEffCorrCl = GetRawEtNoEffCorrPerNClGraph();
  graphRawNoEffCorrCl->Draw("P same");
  TGraphErrors *graphRawAllNoEffCorrCl = GetRawEtAllNoEffCorrPerNClGraph();
  graphRawAllNoEffCorrCl->Draw("P same");
  TLegend *leg6 = new TLegend(0.416667,0.173077,0.62069,0.396853);
  leg6->SetFillStyle(0);
  leg6->SetFillColor(0);
  leg6->SetBorderSize(0);
  leg6->SetTextSize(0.03);
  leg6->SetTextSize(0.038682);
  leg6->AddEntry(partialCorrEtPerNCl,"E_{T}#delta_{TM}/#epsilon f_{nonlin}","P");
  leg6->AddEntry(graphRawNoEffCorrCl,"E_{T}#delta_{TM}/ f_{nonlin}","P");
  leg6->AddEntry(graphRawAllCl,"E_{T}/#epsilon f_{nonlin}","P");
  leg6->AddEntry(graphRawAllNoEffCorrCl,"E_{T}/f_{nonlin}","P");
  leg6->Draw();
  TString corrplotname5 = "/tmp/RawEtVsNcl"+detector+".png";
  c5->SaveAs(corrplotname5.Data());

  TCanvas *c7 = new TCanvas("c7", "Matched tracks",1200, 400);
  c7->SetTopMargin(0.02);
  c7->SetRightMargin(0.02);
  c7->SetBorderSize(0);
  c7->SetFillColor(0);
  c7->SetFillColor(0);
  c7->SetBorderMode(0);
  c7->SetFrameFillColor(0);
  c7->SetFrameBorderMode(0);
  c7->Divide(2);
  c7->cd(1);
  TGraphErrors *graphTotalPerNCh = GetTotalTracksPerNChGraph();
  graphTotalPerNCh->GetHistogram()->GetXaxis()->SetTitle("N_{Ch}");
  graphTotalPerNCh->GetHistogram()->GetYaxis()->SetTitle("Number of tracks/N_{Ch}");
  graphTotalPerNCh->GetHistogram()->SetMaximum(0.06*7.5/scale);//so PHOS and EMCal scales are similar...
  graphTotalPerNCh->Draw("AP");
  TGraphErrors *graphMatchedPerNCh = GetMatchedTracksPerNChGraph();
  graphMatchedPerNCh->Draw("P same");
  TGraphErrors *graphNotMatchedPerNCh = GetNotMatchedTracksPerNChGraph();
  graphNotMatchedPerNCh->Draw("P same");
  TLegend *leg7 = new TLegend(0.176003,0.637152,0.379808,0.86208);
  leg7->SetFillStyle(0);
  leg7->SetFillColor(0);
  leg7->SetBorderSize(0);
  leg7->SetTextSize(0.03);
  leg7->SetTextSize(0.038682);
  leg7->AddEntry(graphMatchedPerNCh,"Number of tracks matched to clusters","P");
  leg7->AddEntry(graphNotMatchedPerNCh,"Number of tracks not matched to clusters","P");
  leg7->AddEntry(graphTotalPerNCh,"Total number of tracks","P");
  leg7->Draw();
  c7->cd(2);
  TGraphErrors *graphTotalPerNCl = GetTotalTracksPerNClGraph();
  graphTotalPerNCl->GetHistogram()->GetXaxis()->SetTitle("N_{Cl}");
  graphTotalPerNCl->GetHistogram()->GetYaxis()->SetTitle("Number of tracks/N_{Cl}");
  graphTotalPerNCl->GetHistogram()->SetMaximum(0.65);
  graphTotalPerNCl->Draw("AP");
  TGraphErrors *graphMatchedPerNCl = GetMatchedTracksPerNClGraph();
  graphMatchedPerNCl->Draw("P same");
  TGraphErrors *graphNotMatchedPerNCl = GetNotMatchedTracksPerNClGraph();
  graphNotMatchedPerNCl->Draw("P same");
  TLine *line = new TLine(graphTotalPerNCl->GetHistogram()->GetXaxis()->GetBinLowEdge(1),0.5,graphTotalPerNCl->GetHistogram()->GetXaxis()->GetBinLowEdge(graphTotalPerNCl->GetHistogram()->GetXaxis()->GetNbins()),0.5);
  line->Draw();
  TString corrplotname7 = "/tmp/TrackMatchCrossCheck"+detector+".png";
  c7->SaveAs(corrplotname7.Data());

}
