// $Id$

/**
   Macro for converting AliRoot digits into HLT RawData. 
   Binary creates for each patch its own file. 
   Singlepatch uses one file per slice (sp=kTRUE). 
   Use altro=kFALSE if you dont want to 
   filter out single timebins

   Run with ALIROOT (not ROOT)
*/

binary(Char_t* inpath,Char_t *outpath,Int_t first,Int_t last,Int_t event,Bool_t sp=kFALSE,Bool_t altro=kTRUE){

  AliHLTTransform::Init(inpath,kTRUE);

  if(sp) {
    singlepatch(inpath,outpath,first,last,event,altro);
    return;
  }

  Char_t name[256];
  const Int_t npatch = 6;

  sprintf(name,"%s/digitfile.root",inpath);
  AliHLTFileHandler *fFileHandler = new AliHLTFileHandler(); 
  fFileHandler->SetAliInput(name);

  for(Int_t slice=first; slice<=last; slice++){
    for(Int_t patch=0;patch<npatch;patch++){
      cerr<<"reading slice: "<<slice<<" patch: "<<patch<<" and storing to: "<<outpath<<"/digits_"<<slice<<"_"<<patch<<".raw"<<endl;
      fFileHandler->Init(slice,patch);      
      sprintf(name,"%s/digits_%d_%d_%d.raw",outpath,event,slice,patch);
      fFileHandler->SetBinaryOutput(name);
      fFileHandler->AliDigits2CompBinary(event,altro);
      fFileHandler->CloseBinaryOutput();      
      fFileHandler->Free();
      cerr<<" done"<<endl;
    }      
  }
  fFileHandler->CloseAliInput();
}

void singlepatch(Char_t* inpath,Char_t *outpath,Int_t first=0, Int_t last=0,Int_t event=0,Bool_t altro=kTRUE)
{
   
  Char_t name[256];
  AliHLTFileHandler *fFileHandler = new AliHLTFileHandler(); 
  sprintf(name,"%s/digitfile.root",inpath);
  fFileHandler->SetAliInput(name);
  
  Int_t patch=-1;
  for(Int_t slice=first; slice<=last; slice++)
    {
      cerr<<"reading slice: "<<slice;
      fFileHandler->Free();
      fFileHandler->Init(slice,patch);
      sprintf(name,"%s/digits_%d_%d_%d.raw",outpath,event,slice,patch);
      fFileHandler->SetBinaryOutput(name);
      fFileHandler->AliDigits2CompBinary(event,altro);
      fFileHandler->CloseBinaryOutput();      
      cerr<<" done"<<endl;
    }
  fFileHandler->CloseAliInput();
  
}

void write2rootfile(Char_t *in,Int_t first,Int_t last,Char_t *path)
{
  //Write new rootfile, using data from the binary files. 

  AliHLTTransform::Init(path);
  char filename[1024];
  sprintf(filename,"%s/digitfile.root",path);
  file = TFile::Open(filename,"READ");
  if(file->IsOpen())
    {
      cout<<"Delete file "<<filename<<endl;
      return;
    }
  for(Int_t slice=first; slice<=last; slice++)
    {
      for(Int_t patch=0; patch<=5; patch++)
	{
	  c = new AliHLTCompress(slice,patch,path);
	  c->WriteRootFile(filename,in);
	  delete c;
	}
    }
  
}

void make_init_file(Char_t *f,Char_t *path="./"){
  AliHLTTransform::MakeInitFile(f,path);
}
