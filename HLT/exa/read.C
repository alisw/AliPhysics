// $Id$

void read(Char_t *path="./",Int_t min=0,Int_t max=35)
{

  for(int slice=0; slice<35; slice++)
    {
      char fname[256];
      sprintf(fname,"%s/digits_%d_0.raw",path,slice);
      file = new AliL3FileHandler();
      if(!file->SetBinaryInput(fname))
	{
	  cerr<<"Error opening file "<<fname<<endl;
	  return;
	}

      int row[2]={0,175};
      file->Init(slice,0,row);

      UInt_t size;
      char name[256];
      AliL3DigitRowData *data = file->CompBinary2Memory(size);
      
      for(Int_t r=0; r<175; r++)
	{
	  UInt_t padrow=data->fRow;
	  AliL3DigitData *dPt = (AliL3DigitData*)data->fDigitData;
	  cout<<"padrow "<<padrow<<" ndigits "<<data->fNDigit<<endl;
	  
	  for(Int_t d=0; d<data->fNDigit; d++)
	    {
	      if(d>0 && dPt[d]->fPad == dPt[d-1]->fPad && dPt[d]->fTime == dPt[d-1]->fTime)
		cout<<"Slice "<<slice<<" padrow "<<padrow<<" pad "<<(int)dPt[d]->fPad<<" time "
		    <<(int)dPt[d]->fTime<<" charge "<<(int)dPt[d]->fCharge<<endl;
	    }
	  
	  file->UpdateRowPointer(data);
	}
      
      file->CloseBinaryInput();
      delete file;
    }
}

void read_ali(Char_t *fname, Int_t sl=0, Int_t sh=35)
{
  AliL3Logger l;
  //l.UnSet(AliL3Logger::kDebug);
  //l.UnSet(AliL3Logger::kAll);
  //l.Set(AliL3Logger::kInformational);
  l.UseStderr();
  //l.UseStream();

#if 0
  //need galice file or alirunfile.root link
  if(AliL3Transform::Init(fname,kTRUE))
  {
    cout << "created temp init file!" << endl;
  }
#endif

  AliL3FileHandler *fileHandler = new AliL3FileHandler();

  if(!fileHandler->SetAliInput(fname))
    {
      cerr<<"Error opening file "<<fname<<endl;
      return;
    }

  Int_t event=0;
  UInt_t nrow=0;

  for(Int_t slice=sl; slice<=sl; slice++){
    for(Int_t patch=0;patch<AliL3Transform::GetNPatches();patch++){

      cerr<<"reading slice: "<<slice<<" patch: "<<patch<<endl;

      fileHandler->Free();
      fileHandler->Init(slice,patch);      
      AliL3DigitRowData *data=fileHandler->AliDigits2Memory(nrow,event);

      cerr<<" found "<< nrow << " rows" <<endl;
    }      
  }
  fileHandler->CloseAliInput();
}
