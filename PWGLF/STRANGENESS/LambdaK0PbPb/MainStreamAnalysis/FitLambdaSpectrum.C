static  int      myAliceBlue    = TColor::GetColor(0,0,192);
static  int      myAliceRed    = TColor::GetColor(192,0,0);
void FitLambdaSpectrum(const char* name, const char * listName = "clambdak0Histo_00", const char* HistName = "h2PtVsMassLambda", const char* SaveName = "SPECTRA_MK_Lambda_7TeV.root")  	
    									    		
////////////////////////zaciatok FitLambdaSpectrum//////////////////////////////////////////////////////////
{

  gROOT->SetStyle("Plain");
  gROOT->LoadMacro("run.C");
  InitAndLoadLibs();
 
  TFile *myFile = new TFile(name);
  myFile->ls();
  TList *myList = myFile->Get(listName); // FIXME: find suffix?

  if(!myList) {
    cout << "Cannot get list " << listName << " from file " << name << endl;
    return;
  }

  TCanvas *myCan = new TCanvas("spectra","lambda");
  myCan->Draw();
  myCan->cd();

  TLatex *myKind = new TLatex(0.15,0.92,"ALICE: LHC10b Pass1 7 TeV"); // FIXME
  myKind->SetNDC();
  myKind->SetTextSize(0.03);
  myKind->SetTextColor(myAliceBlue);
  myKind->Draw();


  TText *textTitle = new TText(0.6,0.86,"");
  textTitle->SetNDC();
  textTitle->SetTextSize(0.04);
  textTitle->Draw();

  TPad *myPad1 = new TPad("myPad1","myPad1",0,0,1,0.85);
  myPadSetUp(myPad1,0.12,0.02,0.02,0.15);
  myPad1->Draw();
  myPad1->cd();

  //    gStyle->SetOptTitle(0);
  //    gStyle->SetOptStat(0);

  TH1F * SaveParameters= new TH1F("SaveParameters","parameters",6,0,6);

  ///zaciatok binovania/////////////
  const Int_t range=33;

  Double_t b[range];

  b[0]=0;     
  b[1]=0.5;
  for (Int_t i = 2;i<17;i++)
    {
      b[i]=b[i-1]+0.1;
      cout<<b[i]<<endl;
    }

  for (Int_t i = 17;i<27;i++)
    {
      b[i]=b[i-1]+0.2;
      cout<<b[i]<<endl;
    }
  b[27]=4.5;
  b[28]=5.0;
  b[29]=5.5;
  b[30]=6.5;
  b[31]=8.0;
  b[32]=12.0;


  //koniec binovania
  TH1F *DrawSpectrumLambda = new TH1F("DrawSpectrumLambda","#Lambda Spectrum;p_{t} [GeV/c]; N",32,b);

  ///////////////////////////////////////////////run over all bins//////////////////////////////////////////////////////////////
  int iLoBin = 15, iHiBin = 16;
  //	int iLoBin = 71, iHiBin = 72, hMax = 250;

  for (Int_t rBin = 4;rBin<33;rBin++){
    //for (Int_t rBin = 36;rBin<51;rBin++){


    myBinCounting(myList, DrawSpectrumLambda, rBin,iLoBin ,iHiBin,HistName,SaveParameters);  
  
    if ((rBin>=1)&&(rBin<16)) {iLoBin = iLoBin+2; iHiBin = iHiBin+2;}
	
    if  (rBin==16) {iLoBin = iLoBin+2; iHiBin = iHiBin+4;}

    if ((rBin>=17)&&(rBin<26)){iLoBin = iLoBin+4; iHiBin = iHiBin+4;}

    if  (rBin==26) {iLoBin = iLoBin+4; iHiBin = iHiBin+10;}

    if  (rBin==27) {iLoBin = iLoBin+10; iHiBin = iHiBin+10;}

    if  (rBin==28) {iLoBin = iLoBin+10; iHiBin = iHiBin+10;}

    if  (rBin==29) {iLoBin = iLoBin+10; iHiBin = iHiBin+20;}

    if  (rBin==30) {iLoBin = iLoBin+20; iHiBin = iHiBin+30;}

    if  (rBin==31) {iLoBin = iLoBin+30; iHiBin = iHiBin+80;}



    char saveFileName[60]; 
    char saveFileName1[60];   
    sprintf(saveFileName,"hm_InvariantMassLambdafit_bin%d.gif",rBin);
    sprintf(saveFileName1,"hm_InvariantMassLambdafit_bin%d.root",rBin);
    myCan->SaveAs(saveFileName); 
    myCan->SaveAs(saveFileName1);
    myCan->Clear();
 	
  }
  //////////////////////////koniec cyklu cez jednotlive biny/////////////////////////////



  //////////////////////////vykreslenie a ulozenie histogramu///////////////////////////
  DrawSpectrumLambda->Draw();

  TFile *SPECTRA_MK_Lambda = new TFile(SaveName,"RECREATE");
  SPECTRA_MK_Lambda->WriteObject(DrawSpectrumLambda,"DrawSpectrumLambda");

  SPECTRA_MK_Lambda->Close();
  //////////////////////////koniec ulozenia a vykreslenia/////////////////////////////////


}
/////////end FitLambdaSpectrum///////////////////////////////////////////////////////////////////////////////  






//////////////////////funkcia pre fitovanie////////////////////////////////////////////////////////


void myBinCounting(TList *myList = 0, TH1F *DrawSpectrumLambda = 0, Int_t rBin = 0, Int_t iLoBin =0, Int_t iHiBin =0, const char* HistName =0,TH1F *SaveParameters = 0){


  TH2F *h2PtVsMassLambda = (TH2F*)myList->FindObject(HistName);
  TH1D *h1MassLambdaPtBin = h2PtVsMassLambda->ProjectionX("h1MassLambdaPtBin",iLoBin,iHiBin);
  h1MassLambdaPtBin->SetFillColor(0);
  //	h1MassLambdaPtBin->SetTitle("h1MassLambdaPtBin;invariant mass (GeV/c^{2});Counts");
  //h1MassLambdaPtBin->Draw("SAME");
	   

  ///////////////////////////////////////Gauss+2nd pol fit////////////////////////////////////////////

  TH1D *myHistGauss = h1MassLambdaPtBin;
  //Double_t hMax=myHistGauss->GetMaximum();
  //myHistGauss->SetMaximum(hMax);
  myHistGauss->Draw();
  //   double myLevel = 500, mySlope = 0, mySlope2 = 0 ;
  //   double myNorm  = 100, myMean = 0.497, myWidth = 0.003;
  double myRange[2]={0,0};  
 
  // Gaussain Fit

  double    myLevel = 0, mySlope = 0, mySlope2 = 0; 
  double     myNorm  = 0, myMean  = 1.116, myWidth = 0.003;
  //	myNorm=(myHistGauss->GetMaximum()*TMath::Sqrt(2*TMath::Pi()));
  myNorm=myHistGauss->GetMaximum();
  //     myRange[0]=1.09; myRange[1]=1.14;
  myRange[0]=1.095; myRange[1]=1.135;

  if(rBin==4)
    {
      //myLevel =-2.46E5; mySlope = 3.39E5; mySlope2 = -1.03E5; myNorm  = 5.12E3;
      myRange[0]=1.10; myRange[1]=1.14;
    }

  if(rBin==3)
    {
      myLevel =-2.91E6; mySlope = 5.21E6; mySlope2 = -2.33E6; myNorm  = 9.32E3;
      myRange[0]=1.104; myRange[1]=1.14;
    }

  char cRange[30],cInterval[30],cIntegration[30],cFitInfo[100], cFitInfo2[100], cFitInfo3[100];
     
  TF1 *fPolGauss = new TF1("fpolgauss",myPolGauss,myRange[0],myRange[1],6);
  fPolGauss->SetTitle("The gaussian function");
  fPolGauss->SetParNames("level","slope","slope2","norm","mean","width");
  fPolGauss->SetParameters(myLevel,mySlope,mySlope2,myNorm,myMean,myWidth);
  fPolGauss->SetLineStyle(2);
  fPolGauss->SetLineWidth(2);
  fPolGauss->SetLineColor(602);
  ///    fPolGauss->Draw("SAME");
  myHistGauss->Fit(fPolGauss,"IREM");
  fPolGauss->Draw("SAME");

  //myHistGauss->Fit(fPolGauss,"LIREM");
     
  myNorm  = fPolGauss->GetParameter(3);
  myMean  = fPolGauss->GetParameter(4);
  myWidth = fPolGauss->GetParameter(5);



  /////store parameters for next fit if own will fall////
  //if ((myMean-(4*myWidth))<(myMean+(4*myWidth)))
  //if (myWidth>0.0005)
  if ((myWidth>0.0005)&&(myWidth<0.004))
    {
      for (Int_t i=1;i<7;i++) SaveParameters->SetBinContent(i,0);
      SaveParameters->SetBinContent(1,fPolGauss->GetParameter(0));
      SaveParameters->SetBinContent(2,fPolGauss->GetParameter(1));
      SaveParameters->SetBinContent(3,fPolGauss->GetParameter(2));
      SaveParameters->SetBinContent(4,fPolGauss->GetParameter(3));
      SaveParameters->SetBinContent(5,fPolGauss->GetParameter(4));
      SaveParameters->SetBinContent(6,fPolGauss->GetParameter(5));


    }
  ///////////////////////////////////////////////////////

  ///////////if parameters are wrong use previous/////
  //if ((myMean-(4*myWidth))>=(myMean+(4*myWidth)))
  //if (myWidth<=0.0005)
  if ((myWidth<=0.0005)||(myWidth>=0.004))
    {
      cout<<"Previous parameters"<<endl;
      for (Int_t i=1;i<7;i++) cout<<SaveParameters->GetBinContent(i)<<endl;
      myLevel=SaveParameters->GetBinContent(1);
      mySlope=SaveParameters->GetBinContent(2);
      mySlope2=SaveParameters->GetBinContent(3);
      myNorm=SaveParameters->GetBinContent(4);
      myMean=SaveParameters->GetBinContent(5);
      myWidth=SaveParameters->GetBinContent(6);

      cout<<"nova sirka"<<myWidth<<endl;
    }
  ////////////////////////////////////////////////////

     
  TF1 *fGauss = new TF1("fgauss",myGauss,myRange[0],myRange[1],3);
  fGauss->SetTitle("The gaussian function");
  fGauss->SetParNames("norm","mean","width");
  fGauss->SetParameters(myNorm,myMean,myWidth);
  fGauss->SetLineColor(myAliceRed);
  fGauss->SetFillColor(myAliceRed);
  fGauss->SetFillStyle(3005);
  fGauss->SetLineWidth(1);
  fGauss->Draw("SAME");
     
  int    NumSigmaSignal      = 4;
  int    NumSigmaBackground   =6;
  double integralWidth = NumSigmaSignal*myWidth;
     
  double binWidth      = myHistGauss->GetBinWidth(1);
  double integralGauss = fGauss->Integral(myMean-integralWidth,myMean+integralWidth);
  sprintf(cRange,"%d#sigma range (GeV/c^{2})",NumSigmaSignal);
  sprintf(cInterval,"[%.3f;%.3f]",myMean-integralWidth,myMean+integralWidth);
  sprintf(cIntegration,"integral: %.0f",integralGauss/binWidth);
  sprintf(cFitInfo,"#Chi^{2}/ndf = %.1f/%d",fPolGauss->GetChisquare(),fPolGauss->GetNDF());
  //     sprintf(cFitInfo2,"mean = %.4f #pm %.4f",fPolGauss->GetParameter(4),fPolGauss->GetParError(4));
  //     sprintf(cFitInfo3,"width = %.5f #pm %.5f",fPolGauss->GetParameter(5),fPolGauss->GetParError(5));

  cout<<cRange<<endl;
  cout<<cInterval<<endl;
  cout<<cIntegration<<endl;
  cout<<cFitInfo<<endl;
  //cout<<cFitInfo2<<endl;
  //cout<<cFitInfo3<<endl;


  //break;

  ////////////////////////////////////////koniec Gauss+2nd pol fit////////////////////////////////////


  Double_t GaussMean =myMean; //fPolGauss->GetParameter(4); //stredna hodnota gausovskeho fitu 
  Double_t GaussSigma =myWidth; //fPolGauss->GetParameter(5); //sigma gausovskeho fitu 

  Double_t LeftSide = myHistGauss->FindBin(GaussMean - NumSigmaSignal*GaussSigma);  // hranice signalu
  Double_t RightSide = myHistGauss->FindBin(GaussMean + NumSigmaSignal*GaussSigma);

  Double_t NumberBinsSignal = (RightSide - LeftSide) + 1; //pocet binov ktore obsahuju signal

  Double_t LeftB_Side_Left = myHistGauss->FindBin(myRange[0]); //bacground on left side of peak
  Double_t LeftB_Side_Right = myHistGauss->FindBin(GaussMean - NumSigmaBackground*GaussSigma);

  Double_t RightB_Side_Left = myHistGauss->FindBin(GaussMean + NumSigmaBackground*GaussSigma); //bacground on Right side of peak
  Double_t RightB_Side_Right = myHistGauss->FindBin(myRange[1]);



  Double_t NumberBinsBackground = (LeftB_Side_Right - LeftB_Side_Left) + 1 + (RightB_Side_Right - RightB_Side_Left) + 1; // pocet binov co 


  //
  // Bin Counting Method 
  //

  // Compute the whole content for Signal in 4 sigma 
  Float_t SignalContentCount = h1MassLambdaPtBin->Integral(LeftSide,RightSide );


  // Create a histo with the Noise only.   

  Double_t  nCountNoise   = 0;
  TH1D *myHistNoise = new TH1D(*h1MassLambdaPtBin);
  myHistNoise->SetName("myHistNoise");
  myHistNoise->SetTitle("myHistNoise");
  myHistNoise->Scale(0);
  myHistNoise->SetFillColor(myAliceBlue);



  for (Int_t iBin = LeftB_Side_Left; iBin <= LeftB_Side_Right; iBin++) {
   
    nCountNoise      +=h1MassLambdaPtBin ->GetBinContent(iBin);
    myHistNoise->SetBinContent(iBin,h1MassLambdaPtBin->GetBinContent(iBin));
    myHistNoise->SetBinError(iBin,h1MassLambdaPtBin->GetBinError(iBin));
  }

  for (Int_t iBin = RightB_Side_Left; iBin <= RightB_Side_Right; iBin++) {
    
    nCountNoise      +=h1MassLambdaPtBin ->GetBinContent(iBin);
    myHistNoise->SetBinContent(iBin,h1MassLambdaPtBin->GetBinContent(iBin));
    myHistNoise->SetBinError(iBin,h1MassLambdaPtBin->GetBinError(iBin));
  }

  // Create a histo with the Signal in 4 sigma.
  int totInHistSignal = 0;
  TH1D *myHistSignal = new TH1D(*h1MassLambdaPtBin);
  myHistSignal->SetName("myHistSignal");
  myHistSignal->SetTitle("myHistSignal");
  myHistSignal->Scale(0);
  myHistSignal->SetFillStyle(3004);


  for (Int_t iBin = LeftSide; iBin <= RightSide; iBin++){
    myHistSignal->SetBinContent(iBin,h1MassLambdaPtBin->GetBinContent(iBin));
    myHistSignal->SetBinError(iBin,h1MassLambdaPtBin->GetBinError(iBin));    
    totInHistSignal += h1MassLambdaPtBin->GetBinContent(iBin);
   
  }


  //
  // Background Fit
  //


  cout<<"zaciatok background fit casti "<<endl;
  cout<<endl;
  cout<<endl;

  //   if  (rOrderPoly == 1 ) {
  //     TF1*fitNoise = new TF1("fitNoise",myPol1,LeftB_Side_Left,RightB_Side_Right,4);
  //
  //     fitNoise->FixParameter(2,(GaussMean - 6*GaussSigma));
  //     fitNoise->FixParameter(3,(GaussMean + 6*GaussSigma));   
  //   }


  TF1*fitNoise = new TF1("fitNoise",myPol2,myRange[0],myRange[1],5);
  fitNoise->SetParLimits(2,-1E8,0);
  fitNoise->SetParameter(0,fPolGauss->GetParameter(0));
  fitNoise->SetParameter(1,fPolGauss->GetParameter(1));
  fitNoise->SetParameter(2,fPolGauss->GetParameter(2));

  fitNoise->FixParameter(3,(GaussMean - NumSigmaBackground*GaussSigma));
  fitNoise->FixParameter(4,(GaussMean + NumSigmaBackground*GaussSigma));    



  Double_t inteGrapmyHistNoise = 0;
  Double_t inteGrapmyHistSignal = 0;

    
  myHistNoise->Fit("fitNoise","irem+");
  //myHistNoise->Fit("fitNoise","Lirem+");

  Double_t hMax=myHistGauss->GetMaximum();
  myHistNoise->SetMaximum(hMax);

  myHistNoise->GetYaxis()->SetTitle("counts");
  myHistNoise->GetYaxis()->SetTitleOffset(1.2);
  myHistNoise->SetMaximum(hMax);
  myHistNoise->Draw("HIST");
  myHistSignal->SetFillColor(myAliceRed);
  myHistSignal->Draw("A BAR E SAME");
  myHistGauss->SetLineWidth(1);
  myHistGauss->Draw("A H E SAME");

  //fitNoise->Draw("same");

  //****************************************
  // Bin Counting
  //****************************************

  // Noise Under Peak 

  //     TF1*fNoise = new TF1("fNoise","[0]+x*[1]",LeftB_Side_Left,RightB_Side_Right); 
  //    fNoise->FixParameter(0,fitNoise->GetParameter(0));
  //     fNoise->FixParameter(1,fitNoise->GetParameter(1));



  TF1*fNoise = new TF1("fNoise","[0]+x*[1]+x*x*[2]",myRange[0],myRange[1]); 
  fNoise->FixParameter(0,fitNoise->GetParameter(0));
  fNoise->FixParameter(1,fitNoise->GetParameter(1));
  fNoise->FixParameter(2,fitNoise->GetParameter(2));


  fNoise->Draw("same");

  //   Double_t nNoiseUnderPeak    = fNoise->Integral(GaussMean - NumSigmaSignal*GaussSigma,GaussMean + NumSigmaSignal*GaussSigma)/(myHistNoise->GetBinWidth(1));
  Double_t lBinWidth = myHistNoise->GetBinWidth(1);
  cout<<"sirka binu"<< lBinWidth<<endl;
  Double_t nNoiseUnderPeak    = fNoise->Integral(h1MassLambdaPtBin->GetBinCenter(LeftSide)-0.5*lBinWidth,h1MassLambdaPtBin->GetBinCenter(RightSide)+0.5*lBinWidth)/lBinWidth;


  Double_t nNoiseUnderPeakErr = TMath::Sqrt(nNoiseUnderPeak);


  // Signal
  Double_t   Signal[2]; // Signal[0]=signal and Signal[1]=error
  Signal[0] = Signal[1]=0.0;
  Signal[0] = totInHistSignal-nNoiseUnderPeak;
  Signal[1] = TMath::Sqrt(Signal[0]+nNoiseUnderPeak);

  //	DrawSpectrumLambda->SetBinContent(rBin,Signal[0]/DrawSpectrumLambda->GetBinWidth(rBin);
  //	DrawSpectrumLambda->SetBinError(rBin,Signal[1]/DrawSpectrumLambda->GetBinWidth(rBin));
  DrawSpectrumLambda->SetBinContent(rBin,Signal[0]);
  DrawSpectrumLambda->SetBinError(rBin,Signal[1]);
     
  printf("BoInfo: Noise intervals [%.3f;%.3f],[%.3f;%.3f]\n",myRange[0],(GaussMean - NumSigmaBackground*GaussSigma),(GaussMean +NumSigmaBackground *GaussSigma),myRange[1]);
  printf("BoInfo: nBinNoise(left+right)=%d  nCountNoise(left+right)=%f \n",NumberBinsBackground,nCountNoise);
  printf("BoInfo: nBinCentral=%d  totInHistSignal=%f \n",NumberBinsSignal,totInHistSignal);
  printf("BoInfo: Noise under peak %.2f +-%.2f \n",nNoiseUnderPeak, nNoiseUnderPeakErr);
  printf("BoInfo: Signal %f \n",Signal[0]);
   
   
  // printf only and keep only on plot the Signal, noise and S/N.
  Char_t     cIntervalSignal[50], cSignalOverNoise[50], cNoiseUnderPeak[50], cFitNoiseInfo[50];
  Double_t   SignalOverNoise    = 0.0;
  Double_t   SignalOverNoiseErr = 0.0;
   
  sprintf(cIntervalSignal,"Signal [%.3f;%.3f] = %.0f #pm %.0f",LeftSide,RightSide,Signal[0], Signal[1]);
  sprintf(cNoiseUnderPeak,"Noise under peak = %.0f #pm %.0f",nNoiseUnderPeak, nNoiseUnderPeakErr);
  SignalOverNoise    = Signal[0]/nNoiseUnderPeak;
  SignalOverNoiseErr = TMath::Sqrt(Signal[0]+2*nNoiseUnderPeak);
  sprintf(cSignalOverNoise,"S/N= %.2f #pm %.2f",SignalOverNoise,SignalOverNoiseErr);
  cout<<"cSignalOverNoise"<<cSignalOverNoise<<endl;
  sprintf(cFitNoiseInfo,"#Chi^{2}/ndf = %.1f/%d",fitNoise->GetChisquare(),fitNoise->GetNDF());
   
  cout<<cIntervalSignal<<endl;
  cout<<cNoiseUnderPeak<<endl;
  cout<<cFitNoiseInfo<<endl;

}

/////////////////////////////////////koniec funkcie////////////////////////////////////////



/////////////////////////////////////////funkcie///////////////////////////////////////////////////////


void myPadSetUp(TPad *currentPad=0, float rLeft=0, float rTop=0, float rRight=0, float rBottom=0){
  currentPad->SetLeftMargin(rLeft); 
  currentPad->SetTopMargin(rTop);  
  currentPad->SetRightMargin(rRight);
  currentPad->SetBottomMargin(rBottom);
  return;
}


double myGauss(double *x, double *par){
  double absc  = x[0];
  double norm  = par[0];
  double mean  = par[1];
  double width = par[2];
  //  double mygauss = (norm/TMath::Sqrt(2*TMath::Pi()))*TMath::Exp(-0.5*(absc-mean)*(absc-mean)/(width*width));
  double mygauss = (norm)*TMath::Exp(-0.5*(absc-mean)*(absc-mean)/(width*width));
  return mygauss;
}

double myPolGauss(double *x, double *pars){
  double absc   = x[0];
  double level  = pars[0];
  double slope  = pars[1];
  double slope2 = pars[2];
  double norm   = pars[3];
  double mean   = pars[4];
  double width  = pars[5];
  //  double ordo = level + slope * absc +slope2 * absc*absc + (norm/TMath::Sqrt(2*TMath::Pi()))*TMath::Exp(-0.5*(absc-mean)*(absc-mean)/(width*width));
  double ordo = level + slope * absc +slope2 * absc*absc + (norm)*TMath::Exp(-0.5*(absc-mean)*(absc-mean)/(width*width));
  return ordo;
}

double myPol1(double *x, double *pars) {

  double absc    = x[0];
  double level   = pars[0];
  double slope   = pars[1];
  double abscMin = pars[2];
  double abscMax = pars[3];

  if (absc >= abscMin && absc < abscMax) {
    TF1::RejectPoint();
    return 0;
  }
  
  double ordo = level + slope * absc;
  return ordo;

}

double myPol2(double *x, double *pars) {

  double absc    = x[0];
  double level   = pars[0];
  double slope   = pars[1];
  double slope2  = pars[2];
  double abscMin = pars[3];
  double abscMax = pars[4];

  if (absc >= abscMin && absc < abscMax) {
    TF1::RejectPoint();
    return 0;
  }
 
  double ordo = level + slope * absc +slope2 * absc*absc;
  return ordo;

}



void LoadLibs(){

  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");   
}
