// $Id$

/**
   Macro for converting DDL (digit) files into HLT RawData. 
   In "inpath" expected are the ddl files (called Ev0TPCslice 
   followed by number), the compressed HLT files will be
   in "outpath", and can be used as the input to the tracker.
   Singlepatch uses one file per slice (sp=kTRUE) (which 
   is probably what we want if we run it on the GDCs)

   Note: event is not used yet.
*/


ddl2binary(Char_t* inpath,Char_t *outpath,Int_t first=0,Int_t last=35,Bool_t sp=kTRUE,Int_t event=-1){

  AliHLTTransform::Init(inpath); //expect l3transform.config in "inpath"

  Int_t patchfrom = 0;
  Int_t patchend = 6;
  if(sp){
    patchfrom = -1;
    patchend = 0;
  }

  Char_t name[256];
  sprintf(name,"%s/Ev0TPCslice",inpath);

  //create the file handler
  AliHLTDDLDataFileHandler *fFileHandler = new AliHLTDDLDataFileHandler(); 
  fFileHandler->SetReaderInput(name);

  for(Int_t slice=first; slice<=last; slice++){
    for(Int_t patch=patchfrom;patch<patchend;patch++){
      cerr<<"reading slice: "<<slice<<" patch: "<<patch<<" and storing to: "<<outpath<<"/ddl_digits_"<<slice<<"_"<<patch<<".raw"<<endl;
      fFileHandler->Init(slice,patch);      
      sprintf(name,"%s/digits_%d_%d_%d.raw",outpath,event,slice,patch);
      fFileHandler->SetBinaryOutput(name);
      fFileHandler->DDLData2CompBinary(event);
      fFileHandler->CloseBinaryOutput();
      cerr<<" done"<<endl;
    }      
  }
  fFileHandler->CloseReaderInput();

  delete fFileHandler;
}
