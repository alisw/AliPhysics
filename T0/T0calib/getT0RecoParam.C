void getT0RecoParam(Int_t run)
{
  
  // Read calibration coefficients into the Calibration DB
  // Arguments:
  AliCDBManager* man = AliCDBManager::Instance();
   man->SetDefaultStorage("raw://");
  //  man->SetDefaultStorage("local:///home/alla/alice/Jul14/OCDB/");
  man->SetRun(run);
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("T0/Calib/RecoParam");
  AliT0RecoParam* recoParam = 0x0;
  cout<<" entry "<< entry<<endl;
  if(entry) {
    // load recoParam according OCDB content (single or array)
    //    if (!(recoParam = dynamic_cast<AliT0RecoParam*>(entry->GetObject()))) {
    
    TObjArray* recoParamArray = static_cast<TObjArray*>(entry->GetObject());     
    cout<<" TObjArray* recoParamArray "<<recoParamArray->GetEntriesFast()<<endl;
    for(Int_t ie = 0; ie < recoParamArray->GetEntriesFast()-1; ie++) {
      
      recoParam = static_cast<AliT0RecoParam*>(recoParamArray->UncheckedAt(ie));	
      cout<<ie<<endl;
      cout<<" eq "<<recoParam->GetEq()<<endl;
      //	 recoParam->Dump();
      cout<<" cfd range "<<recoParam->GetLow(300)<<" amplitude "<<recoParam->GetLow(200)<<" "<<recoParam->GetHigh(200)<<endl;
      
      for (int i=0; i<500; i++) 
	if( recoParam->GetLow(i) !=0) cout<<i<<" low "<<recoParam->GetLow(i)<<"   "<<endl;
      for (int i=0; i<500; i++) 
	if( recoParam->GetHigh(i) !=50000)  cout<<i<<" high "<<recoParam->GetHigh(i)<<endl;
      recoParam = 0x0;
    }
  }
  else 
    cout<<" no entry "<<endl;
  
  
}
