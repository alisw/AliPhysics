/* $Id$ */

/**
Macro for converting AliRoot digits into L3 RawData. Binary create for each patch its own file. singlepatch uses one file per slice. 
 */

Binary(char* in,int first, int last,char *path="."){
  char name[256];
  const Int_t npatch = 6;
  AliL3Logger l;
  //l.UnSet(AliL3Logger::kDebug);
  //l.UnSet(AliL3Logger::kAll);
  //l.Set(AliL3Logger::kInformational);
  l.UseStdout();
  //l.UseStream();

  AliL3FileHandler *fFileHandler = new AliL3FileHandler(); 
  fFileHandler->SetAliInput(in);
  AliL3Transform::Init(path);

  for(int slice=first; slice<=last; slice++){
    for(int patch=0;patch<npatch;patch++){
      cerr<<"reading slice: "<<slice<<" patch: "<<patch<<" and storing to: "<<path<<"digits_"<<slice<<"_"<<patch<<".raw"<<endl;
      fFileHandler->Free();
      fFileHandler->Init(slice,patch);      
      sprintf(name,"%s/digits_%d_%d.raw",path,slice,patch);
      fFileHandler->SetBinaryOutput(name);
      fFileHandler->AliDigits2CompBinary();
      fFileHandler->CloseBinaryOutput();      
      cerr<<" done"<<endl;
    }      
  }
  fFileHandler->CloseAliInput();
}
 
void singlepatch(char* in,int first, int last,char *path="",int event=0)
{
  AliL3Logger l;
  //l.UnSet(AliL3Logger::kDebug);
  //l.UnSet(AliL3Logger::kAll);
  //l.Set(AliL3Logger::kInformational);
  //l.UseStdout();
  l.UseStream();
  
  char fname[100];
  sprintf(fname,"%sevent_%d/",path,event);
  char name[256];
  AliL3FileHandler *fFileHandler = new AliL3FileHandler(); 
  fFileHandler->SetAliInput(in);
  Int_t srow[2] = {0,175};
  int patch=0;
  for(int slice=first; slice<=last; slice++)
    {
      cerr<<"reading slice: "<<slice;
      fFileHandler->Free();
      fFileHandler->Init(slice,patch,srow);
      sprintf(name,"%sdigits_%d_%d.raw",fname,slice,patch);
      fFileHandler->SetBinaryOutput(name);
      fFileHandler->AliDigits2CompBinary(event);
      fFileHandler->CloseBinaryOutput();      
      cerr<<" done"<<endl;
    }
  fFileHandler->CloseAliInput();
  
}
