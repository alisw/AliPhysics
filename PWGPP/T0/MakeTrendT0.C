#define NPMTs 24

double GetParameterGaus(TH1F *histo, int whichParameter);

int MakeTrendT0( char *infile, int run, char* ocdbStorage="raw://"){
  /*
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");
  */
  char *outfile = "trending.root";
 
  if(!infile) return -1;
  if(!outfile) return -1;

  TFile *f = TFile::Open(infile,"read");
  if (!f) {
    printf("File %s not available\n", infile);
    return -1;
  }
  TH1I* fTriggerCounter=0;
  //            LOAD HISTOGRAMS FROM QAresults.root   
  TObjArray   *fTzeroObject = (TObjArray*) f->Get("T0_Performance/QAT0chists");
  TObjArray   *fSPDObject   = (TObjArray*) f->Get("SPD_Performance/coutput1");
  
  TH1F* fTzeroORAplusORC =(TH1F*) ((TH1F*) fTzeroObject->FindObject("fTzeroORAplusORC"))->Clone("A"); 
  TH1F* fResolution      =(TH1F*) ((TH1F*) fTzeroObject->FindObject("fResolution"))->Clone("B");
  TH1F* fTzeroORA        =(TH1F*) ((TH1F*) fTzeroObject->FindObject("fTzeroORA"))->Clone("C");
  TH1F* fTzeroORC        =(TH1F*) ((TH1F*) fTzeroObject->FindObject("fTzeroORC"))->Clone("D");

  TH1F* fSPDVertexZ      =(TH1F*) ((TH1F*) fSPDObject->FindObject("hVertexZ"))->Clone("E");
  TH1F* f0TVX =0x0;
  if(fTzeroObject->FindObject("f0TVX")){
    f0TVX            =(TH1F*) ((TH1F*) fTzeroObject->FindObject("f0TVX"))->Clone("F");
  }
  if( fTzeroObject->FindObject("fTriggerCounter") )
    fTriggerCounter  =(TH1I*) ((TH1I*) fTzeroObject->FindObject("fTriggerCounter"))->Clone("F");
  TH2F *fTimeVSAmplitude[NPMTs];//counting PMTs from 0
  TH1D *fAmplitude[NPMTs];
  TH1D *fTime[NPMTs];
  //-----> add new histogram here  



  // sigma fit of resolution, mean T0A, T0C and T0C&A,  mean amplitude for each PMT
  double resolutionSigma = -9999; //dummy init
  double tzeroOrAPlusOrC = -9999; //dummy init
  double tzeroOrA        = -9999; //dummy init
  double tzeroOrC        = -9999; //dummy init
  double meanAmplitude[NPMTs]; //dummy init
  double meanTime[NPMTs]; //dummy init
  double timeDelayOCDB[NPMTs]; //dummy init
  double efficiency0TVX_SPD      = -9999.; //dummy init
  double efficiency0TVX_CINT7 = -9999.;//dummy init
  double efficiency0TVX_CADAND = -9999.;//dummy init
 
  //.................................................................................

  for(int ipmt=1; ipmt<=NPMTs; ipmt++){ //loop over all PMTs
    fTimeVSAmplitude[ipmt-1] =(TH2F*) ((TH2F*) fTzeroObject->FindObject(Form("fTimeVSAmplitude%d",ipmt)))->Clone(Form("E%d",ipmt));
    int nbinsX = fTimeVSAmplitude[ipmt-1]->GetNbinsX();
    int nbinsY = fTimeVSAmplitude[ipmt-1]->GetNbinsY();
    //Amplitude
    fAmplitude[ipmt-1] = (TH1D*) fTimeVSAmplitude[ipmt-1]->ProjectionX(Form("fAmplitude%d",ipmt), 1, nbinsY); 
    meanAmplitude[ipmt-1] = -9999; //dummy init 
    if(fAmplitude[ipmt-1]->GetEntries()>0){
      meanAmplitude[ipmt-1] = fAmplitude[ipmt-1]->GetMean();
    }

    //Time
    fTime[ipmt-1] = (TH1D*) fTimeVSAmplitude[ipmt-1]->ProjectionY(Form("fTime%d",ipmt), 1, nbinsX); 
    meanTime[ipmt-1] = -9999; //dummy init 
    if(fTime[ipmt-1]->GetEntries()>0){
      if(fTime[ipmt-1]->GetEntries()>20){
        meanTime[ipmt-1] = GetParameterGaus((TH1F*) fTime[ipmt-1], 1); // Mean Time
      }else if(fTime[ipmt-1]->GetEntries()>0){
        meanTime[ipmt-1] = fTime[ipmt-1]->GetMean();
      } 
    }
  }



  if(fResolution->GetEntries()>20){
    resolutionSigma = GetParameterGaus(fResolution, 2); //gaussian sigma 
  }else if(fResolution->GetEntries()>0){
    resolutionSigma = fResolution->GetRMS(); //gaussian sigma 
  }

  if(fTzeroORAplusORC->GetEntries()>20){
    tzeroOrAPlusOrC = GetParameterGaus(fTzeroORAplusORC, 1); //gaussian mean 
  }else if(fTzeroORAplusORC->GetEntries()>0){
    tzeroOrAPlusOrC = fTzeroORAplusORC->GetMean();
  }

  if(fTzeroORA->GetEntries()>20){ 
    tzeroOrA = GetParameterGaus(fTzeroORA, 1); //gaussian mean
  }else if(fTzeroORA->GetEntries()>0){
    tzeroOrA = fTzeroORA->GetMean();
  }

  if(fTzeroORC->GetEntries()>20){
    tzeroOrC = GetParameterGaus(fTzeroORC, 1); //gaussian mean 
  }else if(fTzeroORC->GetEntries()>0){
    tzeroOrC = fTzeroORC->GetMean(); //gaussian mean 
  }
  
  if(f0TVX && (fSPDVertexZ->GetEntries()-fSPDVertexZ->GetBinContent(0)-fSPDVertexZ->GetBinContent(fSPDVertexZ->GetNbinsX()+1))>0)
  efficiency0TVX_SPD = (f0TVX->GetEntries()-f0TVX->GetBinContent(0)-f0TVX->GetBinContent(f0TVX->GetNbinsX()+1)) /
  (fSPDVertexZ->GetEntries()-fSPDVertexZ->GetBinContent(0)-fSPDVertexZ->GetBinContent(fSPDVertexZ->GetNbinsX()+1));
  //efficiency of triggers
  double nspd   = double(fSPDVertexZ->GetEntries()-fSPDVertexZ->GetBinContent(0)-fSPDVertexZ->GetBinContent(fSPDVertexZ->GetNbinsX()+1));
  if (fTriggerCounter) {
    double otvx   = double(fTriggerCounter->GetBinContent(fTriggerCounter->GetXaxis()->FindBin("C0TVX-B")));cout << "0TVX: " << otvx << "\t";
    double cint7  = double(fTriggerCounter->GetBinContent(fTriggerCounter->GetXaxis()->FindBin("CINT7-B")));cout << "CINT7: " << cint7 << "\t";
    double cadand = double(fTriggerCounter->GetBinContent(fTriggerCounter->GetXaxis()->FindBin("CADAND-B")));cout << "CADAND: " << cadand << "\n";
    
    efficiency0TVX_SPD    = (nspd>0?otvx/nspd:-1.);    cout << "efficiency0TVX_SPD: " << efficiency0TVX_SPD << "\t";
    efficiency0TVX_CINT7  = (cint7>0?otvx/cint7:-1.);  cout << "efficiency0TVX_CINT7: " << efficiency0TVX_CINT7 << "\t";
    efficiency0TVX_CADAND = (cadand>0?otvx/cadand:-1); cout << "efficiency0TVX_CADAND: " << efficiency0TVX_CADAND << "\n";
  }

  //-----> analyze the new histogram here and set mean/sigma


  f->Close();

  //-------------------- READ OCDB TIME DELAYS ---------------------------
  // Arguments:
  /*  AliCDBManager* man = AliCDBManager::Instance();
  if (gSystem->Getenv("eocdbStorage")!=NULL){
    man->SetDefaultStorage(ocdbStorage);
  }else{
    man->SetDefaultStorage(gSystem->Getenv("ocdbStorage"));
  }
  man->SetRun(run);
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("T0/Calib/TimeDelay");
  AliT0CalibTimeEq *clb = (AliT0CalibTimeEq*)entry->GetObject(); */
  for (Int_t i=0; i<NPMTs; i++){
    timeDelayOCDB[i] = 0;
    //   if(clb)
    //    timeDelayOCDB[i] = 0;//clb->GetCFDvalue(i,0);
  }


 
  //--------------- write walues to the output ------------ 
  TTreeSRedirector* pcstream = NULL;
  pcstream = new TTreeSRedirector(outfile,"recreate");
  if (!pcstream) return -1;


  TFile *x =  pcstream->GetFile();
  x->cd();

  TObjString runType;
  Int_t startTimeGRP=0;
  Int_t stopTimeGRP=0;
  Int_t time=0;
  Int_t duration=0;

  time     = (startTimeGRP+stopTimeGRP)/2;
  duration = (stopTimeGRP-startTimeGRP);
  TVectorD vecAmp(24,meanAmplitude);
  TVectorD vecTime(24,meanTime);
  TVectorD vecDelay(24,timeDelayOCDB);

  (*pcstream)<<"trending"<<
    "run="<<run<<
    "vecAmp.="<<&vecAmp<<
    "vecTime.="<<&vecTime<<
    "vecDelay.="<<&vecDelay;

    (*pcstream)<<"trending"<<
      "resolution="<< resolutionSigma<<
      "tzeroOrAPlusOrC="<< tzeroOrAPlusOrC<<
      "tzeroOrA="<< tzeroOrA<<
      "tzeroOrC="<< tzeroOrC<<      
      "amplPMT1="<<meanAmplitude[0]<<
      "amplPMT2="<<meanAmplitude[1]<<
      "amplPMT3="<<meanAmplitude[2]<<
      "amplPMT4="<<meanAmplitude[3]<<
      "amplPMT5="<<meanAmplitude[4]<<
      "amplPMT6="<<meanAmplitude[5]<<
      "amplPMT7="<<meanAmplitude[6]<<
      "amplPMT8="<<meanAmplitude[7]<<
      "amplPMT9="<<meanAmplitude[8]<<
      "amplPMT10="<<meanAmplitude[9]<<
      "amplPMT11="<<meanAmplitude[10]<<
      "amplPMT12="<<meanAmplitude[11]<<
      "amplPMT13="<<meanAmplitude[12]<<
      "amplPMT14="<<meanAmplitude[13]<<
      "amplPMT15="<<meanAmplitude[14]<<
      "amplPMT16="<<meanAmplitude[15]<<
      "amplPMT17="<<meanAmplitude[16]<<
      "amplPMT18="<<meanAmplitude[17]<<
      "amplPMT19="<<meanAmplitude[18]<<
      "amplPMT20="<<meanAmplitude[19]<<
      "amplPMT21="<<meanAmplitude[20]<<
      "amplPMT22="<<meanAmplitude[21]<<
      "amplPMT23="<<meanAmplitude[22]<<
      "amplPMT24="<<meanAmplitude[23]<<
      "timePMT1="<<meanTime[0]<<
      "timePMT2="<<meanTime[1]<<
      "timePMT3="<<meanTime[2]<<
      "timePMT4="<<meanTime[3]<<
      "timePMT5="<<meanTime[4]<<
      "timePMT6="<<meanTime[5]<<
      "timePMT7="<<meanTime[6]<<
      "timePMT8="<<meanTime[7]<<
      "timePMT9="<<meanTime[8]<<
      "timePMT10="<<meanTime[9]<<
      "timePMT11="<<meanTime[10]<<
      "timePMT12="<<meanTime[11]<<
      "timePMT13="<<meanTime[12]<<
      "timePMT14="<<meanTime[13]<<
      "timePMT15="<<meanTime[14]<<
      "timePMT16="<<meanTime[15]<<
      "timePMT17="<<meanTime[16]<<
      "timePMT18="<<meanTime[17]<<
      "timePMT19="<<meanTime[18]<<
      "timePMT20="<<meanTime[19]<<
      "timePMT21="<<meanTime[20]<<
      "timePMT22="<<meanTime[21]<<
      "timePMT23="<<meanTime[22]<<
      "timePMT24="<<meanTime[23];
    (*pcstream)<<"trending"<<
      "timeDelayPMT1="<<timeDelayOCDB[0]<<
      "timeDelayPMT2="<<timeDelayOCDB[1]<<
      "timeDelayPMT3="<<timeDelayOCDB[2]<<
      "timeDelayPMT4="<<timeDelayOCDB[3]<<
      "timeDelayPMT5="<<timeDelayOCDB[4]<<
      "timeDelayPMT6="<<timeDelayOCDB[5]<<
      "timeDelayPMT7="<<timeDelayOCDB[6]<<
      "timeDelayPMT8="<<timeDelayOCDB[7]<<
      "timeDelayPMT9="<<timeDelayOCDB[8]<<
      "timeDelayPMT10="<<timeDelayOCDB[9]<<
      "timeDelayPMT11="<<timeDelayOCDB[10]<<
      "timeDelayPMT12="<<timeDelayOCDB[11]<<
      "timeDelayPMT13="<<timeDelayOCDB[12]<<
      "timeDelayPMT14="<<timeDelayOCDB[13]<<
      "timeDelayPMT15="<<timeDelayOCDB[14]<<
      "timeDelayPMT16="<<timeDelayOCDB[15]<<
      "timeDelayPMT17="<<timeDelayOCDB[16]<<
      "timeDelayPMT18="<<timeDelayOCDB[17]<<
      "timeDelayPMT19="<<timeDelayOCDB[18]<<
      "timeDelayPMT20="<<timeDelayOCDB[19]<<
      "timeDelayPMT21="<<timeDelayOCDB[20]<<
      "timeDelayPMT22="<<timeDelayOCDB[21]<<
      "timeDelayPMT23="<<timeDelayOCDB[22]<<
      "timeDelayPMT24="<<timeDelayOCDB[23];
    (*pcstream)<<"trending"<<
      "efficiency0TVX_SPD="<<efficiency0TVX_SPD<<
      "efficiency0TVX_CINT7="<<efficiency0TVX_CINT7<<
      "efficiency0TVX_CADAND="<<efficiency0TVX_CADAND;
    
    //-----> add the mean/sigma of the new histogram here      
 
 (*pcstream)<<"trending"<<"\n";
  
  pcstream->Close(); 
 
  delete pcstream; 

  return 0;

}

//_____________________________________________________________________________
double GetParameterGaus(TH1F *histo, int whichParameter){

   int    maxBin   =  histo->GetMaximumBin(); 
   double max      =  (histo->GetBinContent(maxBin-1) + histo->GetBinContent(maxBin) + histo->GetBinContent(maxBin+1))/3;
   double mean     =  histo->GetBinCenter(maxBin); //mean
   double lowfwhm  =  histo->GetBinCenter(histo->FindFirstBinAbove(max/2));
   double highfwhm =  histo->GetBinCenter(histo->FindLastBinAbove(max/2));
   double sigma    = (highfwhm - lowfwhm)/2.35482; //estimate fwhm  FWHM = 2.35482*sigma

   TF1 *gaussfit   = new TF1("gaussfit","gaus", mean - 4*sigma, mean + 4*sigma); // fit in +- 4 sigma window
   gaussfit->SetParameters(max, mean, sigma); 

   if(whichParameter==2)
     histo->Fit(gaussfit,"RQNI");
   else
     histo->Fit(gaussfit,"RQN");

   double parValue = gaussfit->GetParameter(whichParameter); 

   delete gaussfit; 
 
   return parValue;
}


double GetParameterNotGaus(TH1F *histo, int whichParameter){

    if(whichParameter==1) return double(histo->GetMean());
    if(whichParameter==2) return double(0.57735027 * histo->GetRMS());
    if(whichParameter!=1 && whichParameter!=2)return -99999.;
}
    
