Binary(char* in,int first, int last,char *path=""){
  char name[256];
  AliL3FileHandler *fFileHandler = new AliL3FileHandler(); 
  fFileHandler->SetAliInput(in);
  Int_t row[5][2] = {{ 0, 45},{46,77},{78,109},{110,141},{142,173}};
  for(int slice=first; slice<=last; slice++){
    for(int patch=4;patch>=0;patch--){
      cerr<<"reading slice: "<<slice<<" patch: "<<patch;
      fFileHandler->Init(slice,patch,row[patch]);      
      sprintf(name,"%sdigits_%d_%d.raw",path,slice,patch);
      fFileHandler->SetBinaryOutput(name);
      fFileHandler->AliDigits2CompBinary();
      fFileHandler->CloseBinaryOutput();      
      cerr<<" done"<<endl;
    }      
  }
  fFileHandler->CloseAliInput();
}












