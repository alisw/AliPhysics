Binary(char* in,int first, int last,char *path=""){
  char name[256];
  AliL3FileHandler *fFileHandler = new AliL3FileHandler(); 
  fFileHandler->SetAliInput(in);
//  Int_t row[5][2] = {{ 0, 45},{46,77},{78,109},{110,141},{142,173}};
  const Int_t npatch = 6;
  Int_t row[npatch][2] = {{0,31},{32,63},{64,91},{92,119},{120,143},{144,175}};
  for(int slice=first; slice<=last; slice++){
    for(int patch=0;patch<npatch;patch++){
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
