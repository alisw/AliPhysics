void readDigits() 
{
 
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  char filename[100];
  sprintf(filename,"galice.root");
  AliRunLoader* rl = AliRunLoader::Open("galice.root",AliConfig::fgkDefaultEventFolderName,"read");
  if (rl == 0x0)
   {
     cerr<<"Can not open session for file galice.root\n";
     return;
   }

  rl->LoadgAlice();
  gAlice = rl->GetAliRun();
  
  AliSTART* START  = (AliSTART *)gAlice->GetDetector("START");
  
  rl->LoadHeader();
  rl->LoadKinematics("READ");
  Int_t retval;
  AliLoader* lstart = rl->GetLoader("STARTLoader");

  Int_t iNevents=rl->GetNumberOfEvents();
  cout<<"  nevents   "<<iNevents<<endl;

   char nameTD[8],nameTR[8];

  TH1F *hTimediff = new TH1F("hTimediff","Time difference",100,0,200);
  TH1F *hTimePs = new TH1F("hTimePs","Time in Ps",100,-200,200);
  TH1F *hMeanTime = new TH1F("hMeanTime","Time in Ps",100,2000,2500);
  TH1F *hallADC = new TH1F("hallright","ADC summary right",100,0,200);
  TH1F *hADCright = new TH1F("hADCright","ADC right",100,0,200);
  TH1F *hADCleft = new TH1F("hADCleft","ADC left",100,0,200);
  TH1F *hTimeright = new TH1F("hTimeright","Time right",100,2000.,3200.);
  TH1F *hTimeleft = new TH1F("hTimeleft","Time left",100,2000.,3200.);
  TH1F *hBestTimeright = new TH1F("hBestTimeright","First time right",
				  100,2000.,3200);
  TH1F *hBestTimeleft = new TH1F("hBestTimeleft","First time left",
				 100,11000.,13000.);
   
  //  digits = new AliSTARTdigit();
  AliSTARTdigit *digits   ; // digits
  digits = new AliSTARTdigit();
 
  timeRight = new TArrayI(12);
  timeLeft = new TArrayI(12);
  ADCRight = new TArrayI(12);
  ADCLeft = new TArrayI(12);

 // Event ------------------------- LOOP  
  for (Int_t j=0; j<iNevents; j++){
    //  gAlice->GetEvent(j);
    rl->GetEvent(j);
    //   lstart->Dump();
    lstart->LoadHits("READ");
    lstart->LoadDigits("READ");
    sprintf(nameTD,"START_D_%d",j);
    printf("%s\n",nameTD);
  TObject *td = (TObject*)gDirectory->Get(nameTD);
  // td->Dump();
    //   cout<<" td "<<td<<endl;
  td->Read(nameTD);
   digits->Read(nameTD);
   //   digits->Read(nameTD);
  // digits->Print();
  printf("time %d\n",digits->GetTimeDiff());
   
    if(digits->GetTimeDiff()!=999999){
      Int_t timediff = digits->GetTimeDiff();
      //     Double_t timePs=(timediff-128)*10.; // time in Ps channel_width =10ps
      Int_t timePs=(512-timediff)*2.5.; // time in Ps channel_width =10ps
      cout<<"timediff "<<timediff<<" timePs "<<timePs<<endl;
      hTimediff->Fill(timediff);
      hTimePs->Fill(timePs);
      Int_t mean=digits->GetMeanTime();
      cout<<" mean "<<mean<<endl;
      mean=mean*2.5;
      hMeanTime->Fill(mean);
      Int_t br=digits->GetBestTimeRight(); 
      Int_t bl=digits->GetBestTimeLeft();
      cout<<"BestTimeRight "<<br*2.5<<" BestTimeLeft "<<bl*2.5<<endl;
      hBestTimeright->Fill(br*2.5); 
      hBestTimeleft ->Fill(bl*2.5); 
      digits->GetTimeRight(*timeRight );
      digits->GetTimeLeft(*timeLeft );
      digits->GetADCRight(*ADCRight );
      digits->GetADCLeft(*ADCLeft );
	for (Int_t i=0; i<12; i++) 
	{
	  Int_t t=timeRight.At(i);
	  Int_t ADC=ADCRight.At(i);
	  hTimeright->Fill(t*2.5);
	  hADCright->Fill(ADC);
	}
	for (Int_t i=0; i<12; i++) 
	{
	  Int_t ADC=ADCRight.At(i);
	  Int_t t=timeLeft.At(i);
	  hTimeleft->Fill(t*2.5);
	  hADCleft->Fill(ADC);
	}
    }
  }
  
  Hfile = new TFile("Figdigits.root","RECREATE","Histograms for START digits");
  printf("Writting histograms to root file \n");
  Hfile->cd();
  //Create a canvas, set the view range, show histograms
  gStyle->SetOptStat(111111);
  //  TCanvas *c1 = new TCanvas("c1","Alice START Time ",400,10,600,600);
  hADCright->Write();
  hADCleft->Write();
  hTimeright->Write();
  hTimeleft->Write();
  hTimePs->SetXTitle("arriving time, ps");
  hTimePs->SetYTitle("number of events");
  hTimePs->Write();
  hMeanTime->Write();
  hBestTimeright->Write(); 
  hBestTimeleft ->Write(); 
  Hfile->Close();
    
 
} // end of macro




