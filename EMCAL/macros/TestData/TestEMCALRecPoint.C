///
/// \file TestEMCALRecPoint.C
/// \brief RecPoints reading example
///
/// Test Macro, shows how to load EMCal RecPoints and Geometry, and how can we get 
/// some of the parameters and variables.
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

///
/// Main executing method
///
void TestEMCALRecPoint()
{
  //
  // Getting EMCAL Detector and Geometry.
  //
  AliRunLoader *rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),"read");
  
  if (rl == 0x0)
    cout<<"Can not instatiate the Run Loader"<<endl;
  
  rl->LoadgAlice(); // Needed to get geometry, kinematics
  
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>
    (rl->GetDetectorLoader("EMCAL"));
  
  //TGeoManager::Import("geometry.root");
  AliGeomManager::LoadGeometry("geometry.root");  
  
  AliRun * alirun   = rl->GetAliRun(); // Needed to get Geometry
  
  AliEMCALGeometry * geom ;
  if(alirun)
  {
    AliEMCAL * emcal  = (AliEMCAL*)alirun->GetDetector("EMCAL");
    geom = emcal->GetGeometry();
  }
  else 
  {
    cout<<"alirun not available, instantiate"<<endl;
    geom =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM") ;
  } 
  
  // Load RecPoints
  rl->LoadRecPoints("EMCAL");

  // Load Kinematics if MC
  //rl->LoadKinematics();
  
  // Get maximum number of events
  Int_t maxevent =  rl->GetNumberOfEvents();
  cout<<"Number of events "<<maxevent<<endl;
  
  Int_t iEvent     = -1 ;
  Int_t iprim      = -1 ;
  Float_t energy   = -1 ;
  TVector3 gpos ;
  Float_t lmb[2]   = {0,0};
  Int_t nlm        = -1;
  
  AliEMCALRecPoint * rp;
  for ( iEvent=0; iEvent<maxevent; iEvent++)
  {
    cout <<  "======> Event " << iEvent << endl ;
    
    // Load Event
    rl->GetEvent(iEvent);

    // In case of simulation get stack with kinematics
    //AliStack *sta=rl->Stack();
    
    //
    // Fill array of rec. points
    //
    TObjArray *rpoints ;//= emcalLoader->RecPoints();
    
    TTree *treeR = emcalLoader->TreeR();
    TBranch * branchR = treeR->GetBranch("EMCALECARP");	
    branchR->SetAddress(&rpoints);
    branchR->GetEntry(0);
    
    if(!rpoints->GetEntries()) continue;
    
    cout<<"Number of RecPoints "<<rpoints->GetEntries()<<endl;
    
    //
    // Get recpoints from the list      
    //
    for(Int_t irp = 0; irp< rpoints->GetEntries();irp++)
    {
      rp = static_cast<AliEMCALRecPoint *>(rpoints->At(irp)) ;
      //rp = emcalLoader->RecPoint(irp);
      
      if(!rp) 
      {
        cout<<"Null rec point!"<<endl;
        continue;
      }
      
      //
      // Basic parameters
      //
      energy  = rp->GetEnergy();      // Cluster energy
      rp->GetGlobalPosition(gpos);    // Global ALICE xyz position
      rp->GetElipsAxis(lmb);          // Shower shape 
      nlm     = rp->GetNExMax();      // Number of local maxima
      
      printf("\t RP %d, Energy %2.3f; Phi %2.1f, Eta %2.2f; "
             "Shower shape l0 %2.2f, l1 %2.2f;  NLM %d\n", 
             irp,energy,gpos.Phi()*TMath::RadToDeg(),gpos.Eta(),lmb[0],lmb[1],nlm);
      
      //
      // MC 
      //
      Int_t  primMult  = 0;
      Int_t *primaries =  rp->GetPrimaries(primMult);
      Int_t  parMult   = 0;
      Int_t *parents   =  rp->GetParents  (parMult );
      
      printf("\t MC: N primaries %d, N parents %d\n",primMult,parMult);
      
      //if ( !sta ) continue;
      
      for (Int_t ipr=0; ipr < primMult; ipr++) 
      {
        printf("\t \t primary %d, index %d\n",ipr,primaries[ipr]);

//        if ( primaries[ipr] < 0 ) continue ;
//        
//        TParticle *primary = sta->Particle(primaries[ipr]);
//        if(primary) printf("\t \t Primary: E %2.2f, Phi %2.2f, eta %2.2f\n",
//                           primary->Energy(),primary->Phi()*TMath::RadToDeg(),primary->Eta());
      } // primary
      
      for (Int_t ipa = 0; ipa < parMult; ipa++) 
      {
        printf("\t \t parent  %d, index %d\n",ipa,parents[ipa]);
        
//        if ( parents[ipa] < 0 ) continue ;
//        
//        TParticle *primary = sta->Particle(primaries[ipa]);
//        if(primary) printf("\t \t Primary: E %2.2f, Phi %2.2f, eta %2.2f\n",
//                           primary->Energy(),primary->Phi()*TMath::RadToDeg(),primary->Eta());
      } // parents
      
    } // rp loop
    
  } // event loop
}

