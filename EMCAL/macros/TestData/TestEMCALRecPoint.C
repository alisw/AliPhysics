// Test Macro, shows how to load RecPoints, and how can we get 
// some of the parameters and variables.
// Author: Gustavo Conesa

void TestEMCALRecPoint()
{
  
  // Getting EMCAL Detector and Geometry.
  
  AliRunLoader *rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),"read");
  
  if (rl == 0x0)
    cout<<"Can not instatiate the Run Loader"<<endl;
  
  rl->LoadgAlice();//Needed to get geometry
  
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>
    (rl->GetDetectorLoader("EMCAL"));
  
  //TGeoManager::Import("geometry.root");
  AliGeomManager::LoadGeometry("geometry.root");  
  AliRun * alirun   = rl->GetAliRun(); // Needed to get Geometry
  AliEMCALGeometry * geom ;
  if(alirun){
    AliEMCAL * emcal  = (AliEMCAL*)alirun->GetDetector("EMCAL");
    geom = emcal->GetGeometry();
  }
  else {
    cout<<"alirun not available, instantiate"<<endl;
    geom =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM") ;
  } 
  
  //Load RecPoints
  rl->LoadRecPoints("EMCAL");
  //Get maximum number of events
  Int_t maxevent =  rl->GetNumberOfEvents();
  cout<<"Number of events "<<maxevent<<endl;
  
  Int_t iEvent     = -1 ;
  Int_t iprim      = -1 ;
  Float_t energy   = -1 ;
  TVector3 gpos ;
  
  AliEMCALRecPoint * rp;
  for ( iEvent=0; iEvent<maxevent; iEvent++)
  {
    cout <<  " ======> Event " << iEvent << endl ;
    //Load Event
    rl->GetEvent(iEvent);
    //      AliStack *sta=rl->Stack();
    //Fill array of digits
    TObjArray *rpoints ;//= emcalLoader->RecPoints();
    
    TTree *treeR = emcalLoader->TreeR();
    TBranch * branchR = treeR->GetBranch("EMCALECARP");	
    branchR->SetAddress(&rpoints);
    branchR->GetEntry(0);
    
    if(!rpoints->GetEntries()) continue;
    cout <<  " ======> Event " << iEvent << endl ;
    
    cout<<">> Entries "<<rpoints->GetEntries()<<endl;
    
    //Get recpoints  from the list      
    for(Int_t irp = 0; irp< rpoints->GetEntries();irp++){
      rp = static_cast<AliEMCALRecPoint *>(rpoints->At(irp)) ;
      //rp = emcalLoader->RecPoint(irp);
      if(rp != 0){
        {
          energy  = rp->GetEnergy(); //cluster energy
          cout<<"Energy "<<energy<<" PointEnergy "<<rp->GetPointEnergy()<<endl;
          rp->GetGlobalPosition(gpos);//global ALICE xyz position
          Int_t  primMult  = 0;
          //Int_t *primInts =  rp->GetPrimaries(primMult);
          //for (Int_t ipr=0; ipr<primMult; ipr++) 
          //	cout<<"primlist "<<ipr<<" index "<< primInts[ipr]<<endl;
          //iprim=rp->GetPrimaryIndex() ;
          //    cout<<"Selected primary index "<<iprim<<endl;
          //  if(iprim>0){
          //TParticle *primary=sta->Particle(iprim);
          //cout<<"Primary phi "<<primary->Phi()*180/TMath::Pi()<<" Reconstructed phi "<<gpos.Phi()*180/TMath::Pi()<<"; Energy "<<primary->Energy()<<endl;
          //h->Fill((gpos.Phi()-primary->Phi())*180/TMath::Pi());
          cout<<"rec point "<<irp<<"; Energy "<<energy<<" Phi "<<gpos.Phi()*180/TMath::Pi()<<" Eta "<<gpos.Eta()<<" iprim "<<iprim<<endl; 
          //}
          Float_t lmb[2];
          rp->GetElipsAxis(lmb);
          cout<<"lmb0 "<<lmb[0]<<" lmb1 "<<lmb[1]<<endl;
        }
      }
      else
        cout<<"recpoint pointer 0x0"<<endl;
    } 
  }
}

