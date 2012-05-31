// Test Macro, shows how to load Hits and Geometry, and how can we get 
// some of the parameters and variables.
// Author: Gustavo Conesa

void TestEMCALHit()
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
  //else   geom->PrintGeometry();
  
  //Load Hits
  rl->LoadHits("EMCAL");
  
  //Get maximum number of events
  Int_t maxevent =  rl->GetNumberOfEvents();
  cout<<"Number of events "<<maxevent<<endl;
  //maxevent = 8000 ;
  
  
  AliEMCALHit * hit;
  TClonesArray *hits = 0;
  
  for (Int_t iEvent=0; iEvent<maxevent; iEvent++)
    {
      //cout <<  " ======> Event " << iEvent <<endl ;  
      //Load Event
      rl->GetEvent(iEvent);
      Float_t elos=-1;
      Float_t time  = -1 ;
      Int_t id      = -1 ;
      Int_t iSupMod =  0 ;
      Int_t iTower  =  0 ;
      Int_t iIphi   =  0 ;
      Int_t iIeta   =  0 ;
      Int_t iphi    =  0 ;
      Int_t ieta    =  0 ;
      
      //Fill array of hits
      cout <<  " ======> Event " << iEvent << endl; 
      
      
      //Get hits from the list      
      
      //Hits are stored in different branches in the hits Tree, 
      //first get the branch and then access the hits in the branch
      TTree *treeH = emcalLoader->TreeH();	
      if (treeH) {
	// TreeH exists, get the branch
	Int_t nTrack = treeH->GetEntries();  // TreeH has array of hits for every primary
	TBranch * branchH = treeH->GetBranch("EMCAL");
	branchH->SetAddress(&hits);
	for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {
	  branchH->GetEntry(iTrack);
	  //Now get the hits in this branch
	  Int_t nHit = hits->GetEntriesFast();
	  for(Int_t ihit = 0; ihit< nHit;ihit++){
	    hit = static_cast<AliEMCALHit *>hits->At(ihit);//(hits->At(ihit)) ;
	    
	    if(hit != 0){
	      id   = hit->GetId() ; //cell (hit) label
	      elos = hit->GetEnergy(); //amplitude in cell (hit)
	      time = hit->GetTime();//time of creation of hit after collision
	      
	      cout<<"Hit ID "<<id<<" ELoss "<<elos;
	      
	      //Geometry methods  
	      if(geom){
		geom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta); 
		//Gives SuperModule and Tower numbers
		geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
						  iIphi, iIeta,iphi,ieta);
		//Gives label of cell in eta-phi position per each supermodule
		// cout<< "SModule "<<iSupMod<<"; Tower "<<iTower <<"; Eta "<<iIeta
		//<<"; Phi "<<iIphi<<"; Cell Eta "<<ieta<<"; Cell Phi "<<iphi<<endl;
		cout<< ";  SModule "<<iSupMod<<"; Cell Eta "<<ieta<<"; Cell Phi "<<iphi<<endl;
	      }//geom?
	    }//hit?
	    else
	      cout<<"Hit pointer 0x0"<<endl;
	  }//hit loop
	}// track loop
      }//treeH?
    }//event loop
}

