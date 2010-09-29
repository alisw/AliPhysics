alienmerge(const char* path, 
	   const char* pattern, 
	   const char* outfile=0, 
	   Int_t nFiles = 10000, 
	   const char* blacklist1 ="roskilde", 
	   const char* blacklist2 = "hvidovre"){
    
  TGrid::Connect("alien://",0,0,"t");
  TGridResult* result = gGrid->Query(path,pattern);
  result->Print();
  
  TFileMerger m;
  
   if (outfile) m.OutputFile(outfile);
  Int_t i=0;
  
  TObjArray* outArray = new TObjArray();
  //TList list;
  
  while (result->GetKey(i,"turl") && i<nFiles) {
    // TAlienFile* file = TAlienFile::Open(result->GetKey(i,"turl"),"READ");
    //file->ls();
    // TList* list = (TList*)file->Get("listOfhists");

    
    // std::cout << list;// << std::endl;
    
    //
    //    TObjArray* tmp = (TObjArray*)file.Get("FMD");
    // list.Add(tmp);
    
    //    if(i!=blacklist1 && i!=blacklist2) {
    TString test(result->GetKey(i,"turl"));
    test.ToLower();
    if(test.Contains(blacklist1) || test.Contains(blacklist2) ) {
      i++;
      continue;
    }
    m.AddFile(result->GetKey(i,"turl"));
    cout<<i<<"   "<<result->GetKey(i,"turl")<<endl;
      
      // }
      i++;
  }
  if (i)
    m.Merge();
  
  //outArray->Merge(list);
  //TFile fout("MergeTest.root","RECREATE");
  //outArray->Write();
  //fout.Close();
}
