void testTrackSegment (Int_t evt = 0)  
{

//========= Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) 
  {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

//========== Opening galice.root file  
  TFile * file   = new TFile("galice.root"); 

//========== Get AliRun object from file 
  gAlice = (AliRun*) file->Get("gAlice");

//=========== Gets the PHOS object and associated geometry from the file 
  AliPHOSv0 * PHOS  = (AliPHOSv0 *)gAlice->GetDetector("PHOS");
  AliPHOSGeometry * Geom = AliPHOSGeometry::GetInstance(PHOS->GetGeometry()->GetName(),PHOS->GetGeometry()->GetTitle());
  
//========== Creates the Clusterizer
  AliPHOSClusterizerv1 clusterizer; 
  clusterizer.SetEmcEnergyThreshold(0.01) ; 
  clusterizer.SetEmcClusteringThreshold(0.1) ; 
  clusterizer.SetPpsdEnergyThreshold(0.0000001) ; 
  clusterizer.SetPpsdClusteringThreshold(0.0000002) ; 
  clusterizer.SetLocalMaxCut(0.03) ;
  clusterizer.SetCalibrationParameters(0.,0.0000001) ;


//========== Creates the track segment maker
  AliPHOSTrackSegmentMakerv1 tracksegmentmaker ;

//========== Creates the Reconstructioner  
  AliPHOSReconstructioner Reconstructioner(clusterizer,tracksegmentmaker);     

  //=========== Connects the various Tree's for evt
  gAlice->GetEvent(evt);
  //=========== Gets the Digit TTree
  gAlice->TreeD()->GetEvent(0) ;     
  //=========== Gets the number of entries in the Digits array
  Int_t nId = PHOS->Digits()->GetEntries();    
  //  printf("AnaPHOSv0.C> Number of entries in the Digit array is %d \n",nId);
  
  //=========== Do the reconstruction

  AliPHOSDigit * digit ;
  TIter next(PHOS->Digits()) ;
  Float_t Etot=0 ;
  while((digit = (AliPHOSDigit *)next())) Etot+=clusterizer.Calibrate(digit->GetAmp()) ;
  cout <<"Found   " << nId << "  digits with total energy " << Etot << endl ;



    PHOS->Reconstruction(Reconstructioner);  

  //================Make checks===========================
  AliPHOSDigit * digit ;
  TIter next(PHOS->Digits()) ;
  Float_t Etot=0 ;
  while((digit = (AliPHOSDigit *)next())) Etot+=clusterizer.Calibrate(digit->GetAmp() );
  cout <<"Found   " << nId << "  digits with total energy " << Etot << endl ;

  TClonesArray * EmcRP = PHOS->EmcClusters() ;
  Etot = 0.;
  TIter nextemc(EmcRP) ;
  AliPHOSEmcRecPoint * emc ;
  while((emc = (AliPHOSEmcRecPoint *)nextemc())) {
    Etot+=emc->GetTotalEnergy() ;
    TVector3 pos ;
    emc->GetLocalPosition(pos ) ;
  TMatrix Dummy ;  
    emc->GetGlobalPosition(pos,Dummy) ;
  }

  cout << "Found " << EmcRP->GetEntries() << " EMC Clusters with total energy  "<<Etot << endl ; 
  TClonesArray * PpsdRP = PHOS->PpsdClusters() ;
  cout << "Found " << PpsdRP->GetEntries() << " Ppsd Clusters " << endl ; 

  TObjArray * trsegl = PHOS->TrackSegments() ;
  AliPHOSTrackSegment trseg ;

  Int_t NTrackSegments = trsegl->GetEntries() ;
  Int_t index ;
  Etot = 0 ;
  for(index = 0; index < NTrackSegments ; index++){
    trseg = (AliPHOSTrackSegment * )trsegl->At(index) ;
    Etot+= trseg->GetEnergy() ;
    trseg->Print() ;
  } 
  cout << "Found " << trsegl->GetEntries() << " Track segments with total energy "<< Etot << endl ;
  
}


