// $Id$

/**
   Macro for converting AliRoot digits into L3 RawData. 
   Binary creates for each patch its own file. 
   Singlepatch uses one file per slice (sp=kTRUE). 
*/

binary(char* in,int first, int last,char *path=".",Bool_t sp=kFALSE){

  AliL3Transform::Init(path);
  
  if(sp) {
    singlepatch(in,first,last,path);
    return;
  }

  char name[256];
  const Int_t npatch = 6;
  
  AliL3FileHandler *fFileHandler = new AliL3FileHandler(); 
  fFileHandler->SetAliInput(in);

  for(int slice=first; slice<=last; slice++){
    for(int patch=0;patch<npatch;patch++){
      cerr<<"reading slice: "<<slice<<" patch: "<<patch<<" and storing to: "<<path<<"/digits_"<<slice<<"_"<<patch<<".raw"<<endl;
      fFileHandler->Init(slice,patch);      
      sprintf(name,"%s/digits_%d_%d.raw",path,slice,patch);
      fFileHandler->SetBinaryOutput(name);
      fFileHandler->AliDigits2CompBinary();
      fFileHandler->CloseBinaryOutput();      
      fFileHandler->Free();
      cerr<<" done"<<endl;
    }      
  }
  fFileHandler->CloseAliInput();
}

void write2rootfile(char *in,int first,int last,char *path)
{
  //Write new rootfile, using data from the binary files. 

  AliL3Transform::Init(path);
  char filename[1024];
  sprintf(filename,"%s/digitfile.root",path);
  file = TFile::Open(filename,"READ");
  if(file->IsOpen())
    {
      cout<<"Delete file "<<filename<<endl;
      return;
    }
  for(int slice=first; slice<=last; slice++)
    {
      for(int patch=0; patch<=5; patch++)
	{
	  c = new AliL3Compress(slice,patch,path);
	  c->WriteRootFile(filename,in);
	  delete c;
	}
    }
  
}
 
void singlepatch(char* in,int first=0, int last=0,char *path=".",int event=0)
{
  AliL3Logger l;
  //l.UnSet(AliL3Logger::kDebug);
  //l.UnSet(AliL3Logger::kAll);
  //l.Set(AliL3Logger::kInformational);
  l.UseStderr();
  //l.UseStream();
  
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
      sprintf(name,"%s/digits_%d_%d.raw",path,slice,patch);
      fFileHandler->SetBinaryOutput(name);
      fFileHandler->AliDigits2CompBinary(event);
      fFileHandler->CloseBinaryOutput();      
      cerr<<" done"<<endl;
    }
  fFileHandler->CloseAliInput();
  
}

void make_init_file(Char_t *f,Char_t *path="./"){
  AliL3Logger l;
  //l.UnSet(AliL3Logger::kDebug);
  //l.UnSet(AliL3Logger::kAll);
  //l.Set(AliL3Logger::kInformational);
  l.UseStderr();
  //l.UseStream();

  AliL3Transform::MakeInitFile(f,path);
}
