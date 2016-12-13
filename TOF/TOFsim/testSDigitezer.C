testSDigitezer(){
  Int_t nevents = 100000;

 
  AliCDBManager *ocdb = AliCDBManager::Instance();
  ocdb->SetDefaultStorage("alien://folder=/alice/data/2010/OCDB");
  ocdb->SetRun(137161);
  AliRunLoader rl;
  rl.SetRunNumber(137161);
  AliTOFSDigitizer digitezer;
  digitezer.SetExternalRunLoader(&rl);
  digitezer.InitParameters();

  digitezer.PrintParameters();


  Float_t z0=0,x0=0; // hit position in the strip system (0,0 - center of the strip), cm
  
  Float_t geantTime = 0;


  // output
  Int_t nActivatedPads,nFiredPads;
  Int_t nActivatedPads2,nFiredPads2;
  Int_t isFired[6];
  Int_t isFired2[6];
  Bool_t isF[6];
  Int_t nPlace[6],ix[6],iz[6];
  Int_t nPlace2[6],ix2[6],iz2[6];
  Float_t qInduced[6];
  Float_t tofTime[6],tofTime2[6],averageTime;

  TFile *fout = new TFile("digitizerQA.root","RECREATE");
  TTree *tree = new TTree("tree","tree");
  tree->Branch("x",&x0,"x/F"); // position on the strip reference system (cm)
  tree->Branch("z",&z0,"z/F");
  tree->Branch("nFiredPads",&nFiredPads,"nFiredPads/I");
  tree->Branch("nFiredPads2",&nFiredPads2,"nFiredPads2/I");
//   tree->Branch("isFired",isFired,"isFired[nFiredPads]/I");
//   tree->Branch("isFired2",isFired2,"isFired2[nFiredPads2]/I");
  tree->Branch("tofTime",tofTime,"tofTime[nFiredPads]/F");
  tree->Branch("tofTime2",tofTime2,"tofTime2[nFiredPads2]/F");
  tree->Branch("ix",ix,"ix[nFiredPads]/I");
  tree->Branch("iz",iz,"iz[nFiredPads]/I");
  tree->Branch("ix2",ix2,"ix2[nFiredPads2]/I");
  tree->Branch("iz2",iz2,"iz2[nFiredPads2]/I");
  tree->Branch("qInduced",qInduced,"qInduced[nFiredPads2]/F");

  for(Int_t i=0;i < nevents;i++){

    isF[0]=isF[1]=isF[2]=isF[3]=isF[4]=isF[5]=0;

    digitezer.SetEdgeEffect(1);
    digitezer.SetTimeDelayFlag(0);

    z0 = (gRandom->Rndm())*3.5;
    x0 = (gRandom->Rndm())*2.5;

    digitezer.SimulateDetectorResponse(z0, x0, geantTime, nActivatedPads, nFiredPads, isF, nPlace, qInduced, tofTime, averageTime);
    isFired[0]=isF[0]>0;
    isFired[1]=isF[1]>0;
    isFired[2]=isF[2]>0;
    isFired[3]=isF[3]>0;
    isFired[4]=isF[4]>0;
    isFired[5]=isF[5]>0;

    convertNplace(nActivatedPads,nFiredPads,isFired,nPlace,ix,iz,tofTime,qInduced);
    // printf("nActivated = %i -- nFired = %i\n",nActivatedPads,nFiredPads);
       
    digitezer.SetEdgeEffect(2);
    digitezer.SetTimeDelayFlag(1);
    
    isF[0]=isF[1]=isF[2]=isF[3]=isF[4]=isF[5]=0;

    digitezer.SimulateDetectorResponse(z0, x0, geantTime, nActivatedPads2, nFiredPads2, isF, nPlace2, qInduced, tofTime2, averageTime);
    isFired2[0]=isF[0]>0;
    isFired2[1]=isF[1]>0;
    isFired2[2]=isF[2]>0;
    isFired2[3]=isF[3]>0;
    isFired2[4]=isF[4]>0;
    isFired2[5]=isF[5]>0;

    convertNplace(nActivatedPads2,nFiredPads2,isFired2,nPlace2,ix2,iz2,tofTime2,qInduced);
    //  printf("nActivated = %i -- nFired = %i\n",nActivatedPads,nFiredPads);

    tree->Fill();
  }

  fout->cd();
  tree->Write();
  fout->Close();
}

convertNplace(Int_t nActivatedPads,Int_t nFiredPads,Int_t isfired[6],Int_t nplace[6],Int_t ix[6],Int_t iz[6],Float_t tofTime[6],Float_t qInduced[6]){

  Int_t ipad=0;

  for(Int_t i=0;i<6;i++){
    if(isfired[i]==0)continue;
    switch(nplace[i]){
    case 73:
      ix[ipad]=0;
      iz[ipad]=0;
      break;
    case 25:
      ix[ipad]=0;
      iz[ipad]=-1;
      break;
    case 72:
      ix[ipad]=-1;
      iz[ipad]=0;
      break;
    case 74:
      ix[ipad]=1;
      iz[ipad]=0;
      break;
    case 24:
      ix[ipad]=-1;
      iz[ipad]=-1;
      break;
    case 26:
      ix[ipad]=1;
      iz[ipad]=-1;
      break;
    }
    tofTime[ipad]=tofTime[i];
    qInduced[ipad]=qInduced[i];
    ipad++;
  }


}
