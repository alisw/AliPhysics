void deadChannels(Int_t runnumber){
  
  TGrid::Connect("alien://",0,0,"t");
  AliCDBManager* cdb = AliCDBManager::Instance();
  //cdb->SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB/SHUTTLE/TestShuttle/TestCDB");
  cdb->SetDefaultStorage("local:///home/canute/ALICE/AliRoot/OCDB/SHUTTLE/TestShuttle/TestCDB");
  cdb->SetRun(runnumber);
  
  AliFMDParameters* pars = AliFMDParameters::Instance();
  pars->Init();
  Int_t nDead = 0;
  
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      
      for(UShort_t sec =0; sec < nsec;  sec++) {
	
	for(UShort_t strip = 0; strip < nstr; strip++) {
	    
	  if(pars->IsDead(det,ring,sec,strip)) {
	    std::cout<<Form("FMD%d%c[%d,%d] is dead with gain %f and pedestal %f and noise %f",det,ring,sec,strip,pars->GetPulseGain(det,ring,sec,strip),pars->GetPedestal(det,ring,sec,strip),pars->GetPedestalWidth(det,ring,sec,strip) )<<std::endl;
	    nDead++;
	        
	  }
	  //  if((sec%2 ==0 && strip%128 == 0) || ((strip+1)%128 == 0 && sec%2!=0))
	  //  std::cout<<pars->GetPulseGain(det,ring,sec,strip)<<std::endl;
	}
      }
    }
  }
  
  Float_t reldead = (Float_t)nDead / 51200.;
  
  std::cout<<Form("Found %d dead channels or %f percent",nDead,100*reldead)<<std::endl;

}
