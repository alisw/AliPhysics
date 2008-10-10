// aliroot macro to dump CC-USB raw data
// 

void dumpCCUSB(const int runno  = 117, const int evtnum = -1 )
{
	TH1F *hTDC[40], *hQDC[32], *hScalerCCUSB[2], *hScalerLecroy[12];
	
	for (Int_t i=0;i<40;i++) hTDC[i]          = new TH1F(Form("hTDC%d",i),Form("hTDC%d",i),4096,0.,4096.);
	for (Int_t i=0;i<32;i++) hQDC[i]          = new TH1F(Form("hQDC%d",i),Form("hQDC%d",i),4096,0.,4096.);
	for (Int_t i=0;i<2;i++)  hScalerCCUSB[i]  = new TH1F(Form("hScalerCCUSB%d",i),Form("hScalerCCUSB%d",i),4096,0.,4096.);
	for (Int_t i=0;i<12;i++) hScalerLecroy[i] = new TH1F(Form("hScalerLecroy%d",i),Form("hScalerLecroy%d",i),4096,0.,4096.);
	
	Char_t fname[256];
	sprintf(fname, "/local/data/Run_%09d.Seq_1A.Stream_0.root",runno);
	
	cout << "Reading file: " << fname << endl;
	AliRawReader *rawReader = new AliRawReaderRoot(fname);
	rawReader->SelectEvents( 7 ); // physics events
	
	int firstevent = evtnum;
	int lastevent  = evtnum;
	
	if (evtnum < 0) { // get a bunch of events
		firstevent = 0;
		lastevent  = - evtnum;
	}
	
	if (evtnum == 0) { // get all events
		firstevent = 0;
		lastevent  = 1000000;
	}
  
	Int_t iev = 1;
  
	while ( rawReader->NextEvent() && iev <= lastevent ) 
	{
		AliEMCALCCUSBRawStream *inCCUSB = new AliEMCALCCUSBRawStream(rawReader);
		
		while ( inCCUSB->Next() )
		{
			for (Int_t i=0;i<2;i++)  hScalerCCUSB[i]->Fill(inCCUSB->GetScalerCCUSB(i));
			for (Int_t i=0;i<12;i++) hScalerLecroy[i]->Fill(inCCUSB->GetScalerLecroy(i));
			for (Int_t i=0;i<40;i++) hTDC[i]->Fill(inCCUSB->GetTDC(i));
			for (Int_t i=0;i<32;i++) hQDC[i]->Fill(inCCUSB->GetQDC(i));
		}

		iev ++;
	}
	
	TCanvas *c1 = new TCanvas("c1","",600,600);
	c1->Divide(4,4);
	for (Int_t i=0;i<8;i++) {c1->cd(i+1); hTDC[i]->Draw();}
	c1->cd(9); hScalerCCUSB[0]->Draw(); 
}





