// Test Macro, shows how to load Digits and Geometry, and how can we get 
// some of the parameters and variables.
// Author: Gustavo Conesa

void TestEMCALDigit()
{
  
  // Getting EMCAL Detector and Geometry.
  
  AliRunLoader *rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),"read");
  
  if (rl == 0x0)
    cout<<"Can not instatiate the Run Loader"<<endl;
  
  rl->LoadgAlice();//Needed to get geometry
  
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>
    (rl->GetDetectorLoader("EMCAL"));
  
  TGeoManager::Import("geometry.root");
  
  AliRun * alirun   = rl->GetAliRun(); // Needed to get Geometry
  AliEMCALGeometry * geom ;
  if(alirun){
    AliEMCAL * emcal  = (AliEMCAL*)alirun->GetDetector("EMCAL");
    geom = emcal->GetGeometry();
  }
  
  if (geom == 0) cout<<"Did not get geometry from EMCALLoader"<<endl;
  else   geom->PrintGeometry();
  
  //Load Digits
  rl->LoadDigits("EMCAL");
  
  //Get maximum number of events
  Int_t maxevent =  rl->GetNumberOfEvents();
  cout<<"Number of events "<<maxevent<<endl;
  //maxevent = 10 ;
  
  Int_t iEvent  = -1 ;
  Float_t amp   = -1 ;
  Float_t time  = -1 ;
  Int_t id      = -1 ;
  Int_t iSupMod =  0 ;
  Int_t iTower  =  0 ;
  Int_t iIphi   =  0 ;
  Int_t iIeta   =  0 ;
  Int_t iphi    =  0 ;
  Int_t ieta    =  0 ;
  
  AliEMCALDigit * dig;
  
  for ( iEvent=0; iEvent<maxevent; iEvent++)
    {
      cout <<  " ======> Event " << iEvent << endl ;
      //Load Event
      rl->GetEvent(iEvent);
      
      //Fill array of digits
      TClonesArray *digits = emcalLoader->Digits();	
      
      //Get digits from the list      
      for(Int_t idig = 0; idig< digits->GetEntries();idig++){
	//cout<<">> idig "<<idig<<endl;
	dig = static_cast<AliEMCALDigit *>(digits->At(idig)) ;
	
	if(dig != 0){
	  id   = dig->GetId() ; //cell (digit) label
	  amp  = dig->GetAmplitude(); //amplitude in cell (digit)
	  time = dig->GetTime();//time of creation of digit after collision
	  
	  cout<<"Cell ID "<<id<<" Amp "<<amp<<endl;//" time "<<time<<endl;
	  
	  //Geometry methods  
	  if(geom){
	    geom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta); 
	    //Gives SuperModule and Tower numbers
	    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
					      iIphi, iIeta,iphi,ieta);
	    //Gives label of cell in eta-phi position per each supermodule
	    cout<< "SModule "<<iSupMod<<"; Tower "<<iTower <<"; Eta "<<iIeta
		<<"; Phi "<<iIphi<<"; Cell Eta "<<ieta<<"; Cell Phi "<<iphi<<endl;
	  }
	}
	else
	  cout<<"Digit pointer 0x0"<<endl;
      }
      
    }


}

