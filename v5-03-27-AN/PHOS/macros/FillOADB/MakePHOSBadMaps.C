void MakePHOSBadMaps(){
  //Creates OADB entry with PHOS bad maps for different periods
  //You probably will need alien connection
  
  
  //Init Bad channels map
  AliOADBContainer badmapContainer(Form("phosBadMap"));

  TGrid::Connect("alien://") ;
  
  char key[55] ;
  //For the period LHC10b
  TObjArray * lhc10b = new TObjArray(5) ;
  lhc10b->SetName("LHC10b_pass1") ;
  TFile * fLHC10b = TFile::Open("alien:///alice/cern.ch/user/p/prsnko/BadMaps/BadMap_LHC10b.root") ;
  if(fLHC10b->IsOpen()){
    for(Int_t mod=1;mod<4; mod++){
      sprintf(key,"PHOS_BadMap_mod%d",mod) ;
      TH2I * h = (TH2I*)fLHC10b->Get(key) ;
      if(!h)
        lhc10b->AddAt(0x0,mod);
      else	
        lhc10b->AddAt(new TH2I(*h),mod) ;
    }
    badmapContainer.AppendObject(lhc10b,114737,117223) ;
  }

  //For the period LHC10h
  TObjArray * lhc10h1 = new TObjArray(5) ;
  lhc10h1->SetName("LHC10h_period1") ;
  TFile * fLHC10h1 = TFile::Open("alien:///alice/cern.ch/user/p/prsnko/BadMaps/BadMap_LHC10h_period1.root") ;
  if(fLHC10h1->IsOpen()){
    for(Int_t mod=1;mod<4; mod++){
      sprintf(key,"PHOS_BadMap_mod%d",mod) ;
      TH2I * h = (TH2I*)fLHC10h1->Get(key) ;
      if(!h)
        lhc10h1->AddAt(0x0,mod);
      else	
        lhc10h1->AddAt(new TH2I(*h),mod) ;
    }
    badmapContainer.AppendObject(lhc10h1,136851,137848) ;
  }
  TObjArray * lhc10h234 = new TObjArray(5) ;
  lhc10h234->SetName("LHC10h_period234") ;
  TFile * fLHC10h234 = TFile::Open("alien:///alice/cern.ch/user/p/prsnko/BadMaps/BadMap_LHC10h_period234.root") ;
  if(fLHC10h234->IsOpen()){
    for(Int_t mod=1;mod<4; mod++){
      sprintf(key,"PHOS_BadMap_mod%d",mod) ;
      TH2I * h = (TH2I*)fLHC10h234->Get(key) ;
      if(!h)
        lhc10h234->AddAt(0x0,mod);
      else	
        lhc10h234->AddAt(new TH2I(*h),mod) ;
    }
    badmapContainer.AppendObject(lhc10h234,138732,139517) ;
  }

  badmapContainer.WriteToFile("PHOSBadMaps.root");

}