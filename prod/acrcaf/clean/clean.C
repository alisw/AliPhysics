void clean()
{
  // Login to CAF
  TProof::Open("aliprod@alicecaf");
  
  TMap* map = gProof->GetDataSets("/ALIREC/aliprod");
  
  Int_t max = 0;
  
  TIterator* iter = map->MakeIterator();
  TObjString* str = 0;
  while ((str = (TObjString*) iter->Next()))
  {
    //str->Print();
    
    TObjArray* tokens = str->String().Tokenize("/");
    if (tokens->GetEntries() < 3)
      continue;
      
    //tokens->Print();
    
    TString str2(((TObjString*) tokens->At(2))->String());
    str2.ReplaceAll("run", "");
      
    Int_t run = str2.Atoi();
    
    //Printf("%d", run);
    
    max = TMath::Max(max, run);
  }
  
  Printf("Maximum run number found is %d", max);
  
  Int_t removeUntil = max - 2000;
  
  Printf("Removing all datasets for runs with smaller run numbers than %d", removeUntil);
  
  iter->Reset();
  
  while ((str = (TObjString*) iter->Next()))
  {
    TObjArray* tokens = str->String().Tokenize("/");
    if (tokens->GetEntries() < 3)
      continue;
      
    //tokens->Print();
    
    TString str2(((TObjString*) tokens->At(2))->String());
    str2.ReplaceAll("run", "");
      
    Int_t run = str2.Atoi();
    
    if (run < removeUntil)
    {
      Printf("Removing run %d", run);
      
      gProof->RemoveDataSet(str->String());
    }
  }
}
