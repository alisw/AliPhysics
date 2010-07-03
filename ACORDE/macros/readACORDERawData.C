//////////////////////////////////////////////////////////////////
// 								//
//	readACORDERawData.C macro 				//
//								//
//	Reads the information of ACORDE from raw data file 	//
//	from the 4-32-words (SL0 & MCN mode)			//
//	and draws four histograms (2 by each mode)		//
// 								//
//	Author: Mario Rodriguez Cahuantzi			//
//	E-MaiL: mario.rocah@gmail.com, mrodrigu@mail.cern.ch	//
//								//
//	Created: July 1st. 2010 @ FCFM -BUAP, Puebla, MX	//
//	Last update: created from old AcoReco.C	macro		//
// 								//
//////////////////////////////////////////////////////////////////

void readACORDERawData(char* fileName)
{

	TStopwatch timer;
  	timer.Start();
	// Pointer to rawReader class

  	AliRawReader* rawReader = new AliRawReaderRoot(fileName); 
	
	// Pointer to rawData of Acorde

  	AliACORDERawStream* rawStream  = new AliACORDERawStream(rawReader);    
	
	// Create some histograms

  	TH1F *h1 = new TH1F("h1","ACORDE - Single Muon Hits (SL0)",60,-0.5,59.5);
  	TH1F *h2 = new TH1F("h2","ACORDE - Multiplicity of Acorde Modules (SL0)",61,-1,60);
  	TH1F *h3 = new TH1F("h3","ACORDE - Single Muon Hits (MCN)",60,-0.5,59.5);
  	TH1F *h4 = new TH1F("h4","ACORDE - Multiplicity of Acorde Modules (MCN)",61,-1,60);
	
	// Declare some counters
	
  	size_t contSingle=0;
  	size_t contMulti=0;
  	UInt_t acorde_word[4]; // array to store the 4-words
  	bool word_sl0[60],word_mcn[60]; // boolean array if some hit in module
  	UInt_t shiftword; // shift word

  	for(Int_t m=0;m<60;m++) {word_sl0[m]=0;word_mcn[m]=0;}
  
 
  	Int_t nEvents = rawStream->GetNEvents(fileName);
  
  	printf("File: %s, Number of events: %d \n",fileName,nEvents);


	// Loop over all the events

	for (Int_t i=1; i<=nEvents; i++) 
	{

       		if (!rawReader->NextEvent()) break;
		rawStream->Reset();
		if (!rawStream->Next()) continue;

		acorde_word[0] = rawStream->GetWord(0);
		acorde_word[1] = rawStream->GetWord(1);
		acorde_word[2] = rawStream->GetWord(2);
		acorde_word[3] = rawStream->GetWord(3);

		shiftword = acorde_word[0];
		for(Int_t iaco=0;iaco<30;iaco++)
		{
			word_sl0[iaco] = shiftword & 1;
			shiftword>>=1;
		}
	
		shiftword = acorde_word[1];
		for(Int_t iaco=30;iaco<60;iaco++)
		{
			word_sl0[iaco] = shiftword & 1;
			shiftword>>=1;
		}

		shiftword=acorde_word[2];
		for(Int_t iaco=0;iaco<30;iaco++)
		{
			word_mcn[iaco] = shiftword & 1;
			shiftword>>=1;
		}

		shiftword=acorde_word[3];
		for(Int_t iaco=30;iaco<60;iaco++)
		{
			word_mcn[iaco] = shiftword & 1;
			shiftword>>=1;
		}

		contSingle=0;	
		for(Int_t iaco=0;iaco<60;iaco++) 
		{
			if(word_sl0[iaco]==1) 
			{
				h1->Fill(iaco);
				contSingle++;
			}
		
		}h2->Fill(contSingle);
		contMulti=0;
		for(Int_t iaco=0;iaco<60;iaco++) 
		{
			if(word_mcn[iaco]==1) 
			{
				h3->Fill(iaco);
				contMulti++;
			}
		
		}h4->Fill(contMulti);

	}

	TCanvas *acorde = new TCanvas("ACORDE","ACORDE-Histograms from Raw-Data",1);
	acorde->Divide(2,2);
	acorde->cd(1);
	h1->GetXaxis()->SetTitle("No. of module");
	h1->GetYaxis()->SetTitle("No. of Hits");
	h1->SetFillColor(kRed);
	h1->Draw();

	acorde->cd(2);
	gPad->SetLogy();
	h2->GetXaxis()->SetTitle("No. of fired modules");
	h2->GetYaxis()->SetTitle("No. of events");
	h2->SetFillColor(kBlue);
	h2->Draw();

	acorde->cd(3);
	h3->GetXaxis()->SetTitle("No. of module");
	h3->GetYaxis()->SetTitle("No. of events");
	h3->SetFillColor(kRed);
	h3->Draw();

	acorde->cd(4);
	gPad->SetLogy();
	h4->GetXaxis()->SetTitle("No. of fired modules");
	h4->GetYaxis()->SetTitle("No. of events");
	h4->SetFillColor(kBlue);
	h4->Draw();
	

	delete rawReader;
	delete rawStream;

	timer.Stop();
	timer.Print();
}
