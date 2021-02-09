AliGenerator* CreatePythia8Gen(
                Float_t e_cms,  
                Bool_t kCR, 
                Int_t kProcess
                );
AliGenerator* AddMCGenPythia8(
		Float_t e_cms= 13000.,
		Bool_t kCR= kTRUE,
		Int_t kProcess= 0 
		)
{
	// Add Pythia 6 generator: 
	// -- kProcess=0  MB generation
	//-- kProcess=1  Jet production, pthard generation
	//   -- Color reconnection = ON/OFF

	gSystem->Load("liblhapdf");

	AliGenerator *myGen  = NULL;
	myGen                = CreatePythia8Gen(e_cms, kCR, kProcess);

	return myGen;
}
AliGenerator* CreatePythia8Gen( 
		Float_t e_cms, 
		Bool_t kCR, 
		Int_t kProcess
		)
{


	gSystem->Load("libpythia6");
	gSystem->Load("libEGPythia6");
	gSystem->Load("libAliPythia6");
	gSystem->Load("libpythia8");
	gSystem->Load("libAliPythia8");
	gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
	gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
	gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets")); 

	AliPythia8 *pythia = AliPythia8::Instance();
	AliGenPythiaPlus* myGener = new AliGenPythiaPlus(pythia);

	if(kProcess==0)
		myGener->SetProcess(kPyMbDefault);

	if(kProcess==1){
		myGener->SetProcess(kPyJets);
	}

	//Center-of-mass energy 
	myGener->SetEnergyCMS(e_cms); // in GeV

	//Event list
	myGener->SetEventListRange(-1, -1);

	//Set tune
	// ... default (Monash 2013)// Perugia 2011

	//Set random seed (based on time)
	AliPythia8::Instance()->ReadString("Random:setSeed = on");
	AliPythia8::Instance()->ReadString("Random:seed = 0");

	//Set Color reconnection: ON / OFF
	if(kCR)
		AliPythia8::Instance()->ReadString("ColourReconnection:reconnect = on");
	else
		AliPythia8::Instance()->ReadString("ColourReconnection:reconnect = off");


	return myGener;
}


