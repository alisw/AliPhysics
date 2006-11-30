// $Id$

/**
   Macro for converting 10 bit rawdata into 8 bit rawdata.
*/

void convert(Int_t sl=0,Int_t sh=35,Char_t *path1="./", Char_t *path2="./")
{
  char fname[250];
  AliHLTLogger l;
  //l.UnSet(AliHLTLogger::kDebug);
  //l.UnSet(AliHLTLogger::kAll);
  //l.Set(AliHLTLogger::kInformational);
  l.UseStdout();
  //l.UseStream();

  for(int slice=sl; slice<=sh; slice++)
    {
      for(int patch=0; patch<6; patch++)
	{
	  file = new AliHLTDataHandler();
	  file->Init(slice,patch);
	  sprintf(fname,"%s/digits_%d_%d.raw",path1,slice,patch);
	  file->SetBinaryInput(fname);
	  sprintf(fname,"%s/digits_c8_%d_%d.raw",path2,slice,patch);
	  file->SetBinaryOutput(fname);
	  file->Convert10to8Bit();
	  file->CloseBinaryInput();
	  file->CloseBinaryOutput();
	  
	  delete file;
	}
    }
}

#if 0
void read()
{
  AliHLTMemHandler *file = new AliHLTDataHandler();
  file->SetBinaryInput("/prog/alice/data/new/fixed-slice0/rawdata_8bit/1part/digits_c8_0_5.raw");
  //file->SetBinaryInput("test.raw");
  UInt_t ndigits;
  AliHLTDigitRowData *data = (AliHLTDigitRowData*)file->CompBinary2Memory(ndigits);
  
  cout<<"Number of rows: "<<ndigits<<endl;
  
  AliHLTDigitData *digPt=0;
  for(int i=0; i<32; i++)
    {
      //cout<<(int)data->fRow<<endl;
      digPt = data->fDigitData;
      //cout<<"Number of digits: "<<(int)data->fNDigit<<endl;
      for(int j=0; j<data->fNDigit; j++)
	{
	  cout<<"Row "<<(int)data->fRow<<" pad "<<(int)digPt[j].fPad<<" time "<<(int)digPt[j].fTime<<" charge "<<(int)digPt[j].fCharge<<endl;
	}

      file->UpdateRowPointer(data);
      
    }
}
endif
