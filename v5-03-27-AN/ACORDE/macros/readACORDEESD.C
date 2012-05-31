/*************************************************************

	Macro used for reading the information
	recorded by ACORDE detector from ESD fles and test the
	correction of ACORDE ESDs
	
	Author: 
	
		Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch>

	Created: Sep. 24th 2009 @ CERN


**************************************************************/

void readACORDEESD()
{

	// Time's counter
	
	TStopwatch timer;
	
	timer.Start();

	TH1D *h1 = new TH1D("h1","ACORDE - Single Muon Trigger Hits",60,-0.5,59.5);
  	TH1D *h2 = new TH1D("h2","ACORDE - Single Muon Trigger Hit Multiplicity",60,-1,60);
  	TH1D *h3 = new TH1D("h3","ACORDE - Multi Muon Trigger Hits",60,-0.5,59.5);
  	TH1D *h4 = new TH1D("h4","ACORDE - Multi Muon Trigger Hit Multiplicity",60,-1,60);

	// Pointer to the ESD file
	
	TFile *ef = TFile::Open("AliESDs.root");
    
	// Check if the ESD file is OK
		    
	if (!ef || !ef->IsOpen()) 
	{
		cerr<<"Can't read AliESDs.root !\n"; 
		return 1;
	}
	
	// Pointer to the ESD event
	
	AliESDEvent* fESD = new AliESDEvent();
	
	// Pointer to the esdTree 
	
	TTree* tree = (TTree*) ef->Get("esdTree");
	
	// Check if the  esdTree is Ok
	
	if (!tree) 
	{
		cerr<<"no ESD tree found\n"; 
		return 1;
	}
	
	fESD->ReadFromTree(tree);

	Int_t n=1;

	// Loop over all events
	
	while (tree->GetEvent(n))
	{

		Int_t nMuons = 0;
	
		cout<<endl<<"Processing event number : "<<n++<<endl;
		
		// We select only events triggered by ACORDE
	
		TString ActiveTriggerDetector = fESD->GetFiredTriggerClasses();
		printf("Event:%d, Trigger:%s\n",fESD->GetEventNumberInFile(),ActiveTriggerDetector.Data());

		// if ACORDE trigger is on
		// else if ACORDE is working as readout just comment the condition
		if (ActiveTriggerDetector.Contains("AMU") || ActiveTriggerDetector.Contains("SL0") ) 		
		{

			cout<<endl<<"Processing event number : "<<n++<<endl;
			printf("Trigger Mask:%s\n",ActiveTriggerDetector.Data());

			AliESDACORDE *acordeESD = fESD->GetACORDEData();
			Int_t contMulti = 0;
			for(Int_t i=0;i<60;i++)
			{
				if (acordeESD->GetHitChannel(i)==kTRUE) 
				//if (acordeESD->GetACORDEBitPattern(i))
				{
					h1->Fill(i);
					contMulti++;
				}
			}		
			h2->Fill(contMulti);
		} //Only acorde trigger
	} // Loop over all events from AliESDs.root		

	TCanvas *c1 = new TCanvas();
	c1->Divide(2,1);
	c1->cd(1);
	h1->Draw();
	c1->cd(2);
	h2->Draw();

	timer.Stop();
	timer.Print();	



} // Main Procedure 

